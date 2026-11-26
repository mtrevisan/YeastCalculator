/**
 * Copyright (c) 2026 Mauro Trevisan
 * <p>
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.KneadingInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.domain.YeastInput;
import io.github.mtrevisan.yeastcalculator.plugin.PrefermentPlugin;
import io.github.mtrevisan.yeastcalculator.simulation.DoughOdeSystem;
import io.github.mtrevisan.yeastcalculator.simulation.FoldEventHandler;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;

import java.util.Arrays;


/**
 * Conducts bounded search algorithms to match yeast velocity targets to structural volume plateaus.
 * Dynamically computes rehydration limits based on current cell moisture characteristics.
 */
public class YeastOptimizer{

	// Component densities used to build the dough volume from its composition [g/ml]
	// Search range of the fresh yeast fraction: from 0.01% (a pinch, as used in 24-72 h cold fermentations and in
	// pre-ferments) up to 3.5%
	private static final double MIN_YEAST_RATIO = 0.0001;
	private static final double MAX_YEAST_RATIO = 0.035;
	// Gluten development: ceiling reachable by kneading alone, and time constant of the kneading [min]
	private static final double MAX_KNEADING_DEVELOPMENT = 0.9;
	private static final double KNEADING_TIME_CONSTANT = 10.;
	// Time constant of the development reached by the autolyze [min]: 30 min close about 18% of the gap, 60 min about 33%
	// (assumed value, to be refined with measurements)
	private static final double AUTOLYZE_TIME_CONSTANT = 150.;
	// Fraction of the remaining development gap closed by a fully effective stretch-and-fold set
	private static final double FOLD_DEVELOPMENT_GAIN = 0.25;
	// Rest after which a fold is fully effective [min], and the effectiveness of a fold given right after the
	// previous one
	private static final double FOLD_REST_MINUTES = 30.;
	private static final double MIN_FOLD_EFFECTIVENESS = 0.25;
	// Rest needed after the last fold before the end of the timeline for it to be fully effective [min]
	private static final double FINAL_REST_MINUTES = 60.;
	// Dry flour solids (starch ~1.5, protein ~1.3, fibers ~1.4)
	private static final double DRY_SOLIDS_DENSITY = 1.45;
	// Free water (the flour moisture above the GAB monolayer, and the water added by the recipe)
	private static final double WATER_DENSITY = 1.;
	// Water strictly bound in the GAB monolayer: it is slightly contracted with respect to the free water (approximate
	// value, the effect on the total volume is below 1%)
	private static final double BOUND_WATER_DENSITY = 1.2;
	// Apparent density of salt once dissolved (partial specific volume of NaCl in water)
	private static final double SALT_DENSITY = 2.16;
	private static final double OIL_DENSITY = 0.92;
	// Flour strength W at which the gluten strength factor reaches half of its saturation value
	private static final double W_HALF_SATURATION = 150.;
	// Scale of the gas capacity per unit of dough volume
	private static final double EXPANSION_SCALE = 4.5;
	// Neutral gluten strength/development factor used to build the proofing target of a product (see
	// calculateReferenceGasPotential)
	private static final double REFERENCE_CAPACITY_FACTOR = 0.6;
	// Floor of the kneading contribution to the gas capacity: even a hand-kneaded dough builds a network
	// (folds and rest complete it), while a mixer saturates it at 1
	private static final double KNEADING_FACTOR_FLOOR = 0.6;

	/**
	 * Finds the optimal initial yeast ratio using a Brent Optimizer.
	 */
	public static double findOptimalYeast(final SimulationInputs in, final BakeryProduct targetProduct,
			final PrefermentPlugin activePlugin){
		// 1. If an upstream plugin exists, let it intercept and pre-condition the inputs (Recipe re-balancing)
		final SimulationInputs runtimeInputs;
		final PrefermentPlugin.StatePreConditioner preConditioner;
		if(activePlugin != null){
			final PrefermentPlugin.PrefermentDecorator decorator = activePlugin.interceptAndPreCondition(in);
			runtimeInputs = decorator.decoratedInputs();
			preConditioner = decorator.statePreConditioner();
		}
		else{
			runtimeInputs = in;
			preConditioner = null;
		}

		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(runtimeInputs);
		final double baseMaxGasPotential = calculateMaxGasPotential(inputsSummary(runtimeInputs), runtimeInputs,
			targetProduct);

		// The objective function defines the target for the optimizer to minimize
		// The search runs on the natural logarithm of the yeast ratio: the range spans more than two decades, and a
		// linear search would barely sample the lowest doses used in long cold fermentations
		final UnivariateFunction objective = logYeastAttempt -> {
			// Run the ODE simulation using the current yeast attempt value
			final double yeastAttempt = Math.exp(logYeastAttempt);
			final double[] y = runSimulation(runtimeInputs, gab, baseMaxGasPotential, yeastAttempt, preConditioner);
			// If the simulation crashed or hit a mathematical wall, return a maximum penalty score
			if(y == null)
				return Double.MAX_VALUE;

			return evaluateFitness(runtimeInputs, gab, baseMaxGasPotential, y, targetProduct);
		};

		final BrentOptimizer optimizer = new BrentOptimizer(1.e-4, 1.e-4);
		return Math.exp(optimizer.optimize(
			new MaxEval(200),
			new UnivariateObjectiveFunction(objective),
			GoalType.MINIMIZE,
			new SearchInterval(Math.log(MIN_YEAST_RATIO), Math.log(MAX_YEAST_RATIO))
		).getPoint());
	}

	/**
	 * Scores a final simulated state against the target product.
	 * <p>
	 * The end-of-timeline gas retention ceiling is rebuilt exactly as the integration saw it: the structural gain of
	 * every fold applied within the timeline and the surface skinning of the last stage.
	 * </p>
	 */
	static double evaluateFitness(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double baseMaxGasPotential, final double[] y, final BakeryProduct targetProduct){
		// Calculate total timeline execution [min]
		double totalMinutes = 0.;
		for(final StageInput stage : in.stages())
			totalMinutes += stage.duration();
		totalMinutes *= 60.;

		// Only the folds falling inside the timeline have been applied by the integrator.
		// Counting them manually here lets us bake the structural gain directly into the
		// retention ceiling, without having to re-run the ODE or rebuild a stateful handler.
		int appliedFolds = 0;
		for(final double fold : in.folds())
			if(fold < totalMinutes)
				appliedFolds ++;

		final FoldEventHandler foldHandler = new FoldEventHandler(new double[0]);

		// Pre-multiply the structural gain of every applied fold into the base potential.
		final double foldedMaxGasPotential = baseMaxGasPotential
			* Math.pow(FoldEventHandler.STRUCTURAL_GAIN, appliedFolds);

		final double[] yDotFinal = new double[7];
		// Final gas production speed. The structural gain is already inside foldedMaxGasPotential,
		// so an empty fold handler is used: it contributes a modifier of exactly 1.0 and does not
		// apply any further folds on top of the pre-multiplied ceiling.
		final StageInput lastStage = in.stages()[in.stages().length - 1];
		final double lastStageRH = lastStage.relativeHumidity();
		final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.temperature(), lastStageRH, gab,
			foldedMaxGasPotential, in, foldHandler);
		finalOde.computeDerivatives(totalMinutes, y, yDotFinal);

		// Apply skinning from the final ambient humidity
		final double dynamicMaxPotential = foldedMaxGasPotential
			* DoughOdeSystem.skinningModifier(lastStageRH);

		// Delegate the product-specific cost to the BakeryProduct enum. The degassing of the folds is already in the
		// simulated gas volume, so the proofing target is not lowered by the number of folds.
		final double fermentationFitness = targetProduct.computeFitness(y[2], dynamicMaxPotential,
			calculateReferenceGasPotential(in, targetProduct), y[1], y[0], yDotFinal[2]);

		// Hydration, salt, malt and oil cannot be judged by the fermentation alone: add their formulation cost
		final DoughRecipe recipe = in.recipe();
		// The number of folds is not rewarded for itself: it is judged through the gluten development it brings
		final double developmentIndex = calculateDevelopmentIndex(in, totalMinutes);
		return fermentationFitness
			+ targetProduct.computeFormulationCost(recipe.waterRatio(), recipe.saltRatio(), recipe.maltRatio(),
				recipe.oilRatio())
			+ targetProduct.computeStructuralCost(recipe.waterRatio(), developmentIndex, appliedFolds);
	}

	/**
	 * Gluten development reached at the end of the timeline, in [0, 1].
	 * <p>
	 * Kneading sets the starting value, depending on the mixer and on the kneading time. Every stretch-and-fold set
	 * applied within the timeline then closes a fraction of the remaining gap; a fold given too soon after the
	 * previous one (or after the mixing) finds a dough that has not relaxed and is less effective.
	 * </p>
	 */
	public static double calculateDevelopmentIndex(final SimulationInputs in, final double totalMinutes){
		final double initialDevelopment = getInitialDevelopment(in);

		final double[] folds = in.folds().clone();
		Arrays.sort(folds);
		double gap = 1. - initialDevelopment;
		double previousMinute = 0.;
		for(final double fold : folds){
			if(fold >= totalMinutes)
				break;

			final double relaxation = Math.clamp((fold - previousMinute) / FOLD_REST_MINUTES, MIN_FOLD_EFFECTIVENESS, 1.);
			// A fold given close to the end leaves no time to rest and recover the gas before the shaping: it counts less
			final double lateness = Math.clamp((totalMinutes - fold) / FINAL_REST_MINUTES, 0., 1.);
			gap *= 1. - FOLD_DEVELOPMENT_GAIN * relaxation * lateness;
			previousMinute = fold;
		}
		return 1. - gap;
	}

	private static double getInitialDevelopment(SimulationInputs in){
		final KneadingInput kneading = in.kneading();
		final double mixerDevelopment = Math.min(1., kneading.type().getDevelopmentEfficiency());
		final double kneadingDevelopment = Math.min(MAX_KNEADING_DEVELOPMENT, MAX_KNEADING_DEVELOPMENT * mixerDevelopment
			* (1. - Math.exp(-kneading.durationMinutes() / KNEADING_TIME_CONSTANT)));
		// Autolyze: the flour hydrates and the gluten develops with no work, closing part of the remaining gap
		return 1. - (1. - kneadingDevelopment)
			* Math.exp(-kneading.autolyseMinutes() / AUTOLYZE_TIME_CONSTANT);
	}

	/**
	 * High-level simulation orchestrator that encapsulates internal structural metrics
	 * to prevent Feature Envy in the calling client.
	 */
	public static double[] runSimulation(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double yeastRatio, final BakeryProduct targetProduct){
		return runSimulation(in, gab, yeastRatio, targetProduct, null);
	}

	public static double[] runSimulation(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double yeastRatio, final BakeryProduct targetProduct, final PrefermentPlugin activePlugin){
		final SimulationInputs runtimeInputs;
		final PrefermentPlugin.StatePreConditioner preConditioner;
		if(activePlugin != null){
			final PrefermentPlugin.PrefermentDecorator decorator = activePlugin.interceptAndPreCondition(in);
			runtimeInputs = decorator.decoratedInputs();
			preConditioner = decorator.statePreConditioner();
		}
		else{
			runtimeInputs = in;
			preConditioner = null;
		}

		// Internal parameters are resolved automatically inside the owner class
		final double[] summary = inputsSummary(runtimeInputs);
		final double baseMaxGasPotential = calculateMaxGasPotential(summary, runtimeInputs, targetProduct);

		// Delegates to the low-level simulation logic
		return runSimulation(runtimeInputs, gab, baseMaxGasPotential, yeastRatio, preConditioner);
	}

	/**
	 * The core numerical orchestrator. Executes sequential integration blocks via Dormand-Prince 8(5,3).
	 * Calculates the dynamic biological activation modifier as a direct function of cell moisture.
	 */
	private static double[] runSimulation(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double baseMaxGasPotential, final double yeastRatio,
			final PrefermentPlugin.StatePreConditioner preConditioner){
		final YeastInput yeastInput = in.yeastProperties();
		final double yeastMoisture = yeastInput.yeastMoisture();

		// Rehydration efficiency modifier: reduces the effective initial biomass if the cell has
		// not completed its osmotic swelling, or if it has been left in pure water for so long
		// that autolysis and starvation stress have started. The biologically acceptable window
		// is provided by YeastInput so the optimization layer stays free of domain knowledge.
		//
		// If a PrefermentPlugin preconditioner is present, the plugin has already handled the
		// rehydration (for example through a Poolish or Biga phase) and the initial state is
		// set directly by preConditioner.applyToInitialState. In that case the rehydration
		// penalty is skipped entirely to avoid double-counting the plugin's work.
		double rehydrationEfficiencyModifier = getRehydrationEfficiencyModifier(preConditioner, yeastInput, yeastMoisture);

		// Set initial conditions for state variables: [Biomass X, Sugars S, Gas Volume V, Ethanol EtOH, Damaged starch D, Physiological state h]
		final double initialX = yeastRatio * (1. - yeastMoisture) * rehydrationEfficiencyModifier;

		double totalFlourSugar = 0.;
		for(int i = 0; i < in.flourMatrix().length; i ++)
			totalFlourSugar += in.fractions()[i] * in.flourMatrix()[i].sugar();
		double initialS = totalFlourSugar + (in.recipe().maltRatio() * in.recipe().maltSugarContent());
		final double initialV = 0.;
		final double initialEtOH = 0.;
		// Build the direct dough baseline state array vector
		double initialDamagedStarch = damagedStarchOf(in);

		// Autolyze: the amylases work on the damaged starch while flour and water rest together (at the temperature of
		// the flour), before the yeast is added
		final double autolyseHydrolysis = initialDamagedStarch * DoughOdeSystem.starchHydrolysisFraction(
			in.kneading().autolyseMinutes(), in.flourTemperature(), in.recipe().waterRatio() + gab.flourActiveWater);
		initialDamagedStarch -= autolyseHydrolysis;
		initialS += DoughOdeSystem.STARCH_TO_MALTOSE_YIELD * autolyseHydrolysis;

		// Lag phase: the yeast of a pre-ferment is already adapted, the one added directly to the dough is not
		final double initialPhysiologicalState = (preConditioner != null?
			DoughOdeSystem.ADAPTED_PHYSIOLOGICAL_STATE: DoughOdeSystem.initialPhysiologicalState(yeastMoisture));
		final double[] y = new double[]{initialX, initialS, initialV, initialEtOH, initialDamagedStarch,
			initialPhysiologicalState, 0.};

		// --- INTERCEPTOR TRIGGER FIX ---
		// If an active plugin preconditioner hook is provided, let it execute now!
		// It intercepts y and modifies biomass, sugars, and ethanol based on Phase 1 outcomes.
		if(preConditioner != null){
			preConditioner.applyToInitialState(y, totalFlourSugar);
			// The growth of the yeast in the pre-ferment already used part of its doublings
			if(initialX > 0.)
				y[6] = Math.max(0., Math.log(y[0] / initialX) / Math.log(2.));
		}

		double currentTime = 0.;
		// Create the event handler to intercept stretch and fold timestamps
		final FoldEventHandler foldHandler = new FoldEventHandler(in.folds());
		for(final StageInput stage : in.stages()){
			final DoughOdeSystem ode = new DoughOdeSystem(stage.temperature(), stage.relativeHumidity(),
				gab, baseMaxGasPotential, in, foldHandler);
			// Variable-step size integrator setup
			final DormandPrince853Integrator integrator = new DormandPrince853Integrator(1.e-4, 2.,
				1.e-5, 1.e-5);
			integrator.addEventHandler(foldHandler, 0.01, 1.e-4, 100);

			final double stageEndTime = currentTime + stage.duration() * 60.;
			try{
				// Execute numerical step integration
				integrator.integrate(ode, currentTime, y, stageEndTime, y);
			}
			catch(final Exception e){
				// If a severe arithmetic error (like an un-guarded division by zero) occurs,
				// catch the exception and return null so the objective function applies a max penalty.
				return null;
			}
			currentTime = stageEndTime;
		}
		return y;
	}

	private static double getRehydrationEfficiencyModifier(final PrefermentPlugin.StatePreConditioner preConditioner,
		   final YeastInput yeastInput, final double yeastMoisture){
		double rehydrationEfficiencyModifier = 1.;
		if(preConditioner == null){
			final double rehydrationWindowHours = yeastInput.getRehydrationWindow() / 60.;
			final double[] optimalWindow = yeastInput.getOptimalRehydrationWindow();
			final double minHours = optimalWindow[0] / 60.;
			final double maxHours = optimalWindow[1] / 60.;

			if(yeastMoisture >= 0.65){
				// --- FRESH YEAST PATHWAY (e.g., 70% moisture panetto) ---
				// Cells are awake and active. If left floating in pure water without substrates
				// for too long, autolysis and early starvation acceleration triggers.
				if(rehydrationWindowHours > maxHours)
					rehydrationEfficiencyModifier = Math.max(0.4, 1. - (rehydrationWindowHours - maxHours) * 0.8);
			}
			else{
				// --- DEHYDRATED DRY YEAST PATHWAY (e.g., active dry yeast) ---
				// Structural mummified cells require a rigid timeline to safely re-swell membrane
				// proteins. Below the lower bound: shear trauma. Above the upper bound: starvation.
				if(rehydrationWindowHours < minHours)
					rehydrationEfficiencyModifier = 0.6 + (rehydrationWindowHours / minHours) * 0.4;
				else if(rehydrationWindowHours > maxHours)
					rehydrationEfficiencyModifier = Math.max(0.5, 1. - (rehydrationWindowHours - maxHours) * 0.4);
			}
		}
		return rehydrationEfficiencyModifier;
	}

	/**
	 * Damaged starch of the flour blend, weighted by the flour fractions [g/g_flour].
	 */
	private static double damagedStarchOf(final SimulationInputs in){
		double damagedStarch = 0.;
		for(int i = 0; i < in.flourMatrix().length; i ++)
			damagedStarch += in.fractions()[i] * in.flourMatrix()[i].estimatedDamagedStarch();
		return damagedStarch;
	}

	/**
	 * Volume of the dough before any fermentation [ml/g_flour].
	 */
	private static double doughVolume(final SimulationInputs in){
		final DoughRecipe recipe = in.recipe();
		// The moisture of the flour (which depends on the air humidity) is split as in the GAB model: the bound water
		// and the free water have different densities, and the rest of the gram of flour is dry solids
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in);
		final double drySolids = 1. - gab.flourStrictlyBoundWater - gab.flourActiveWater;
		// Malt is a flour-like powder; masses are per gram of flour, so the volume is the sum of the component volumes
		return (drySolids + recipe.maltRatio()) / DRY_SOLIDS_DENSITY
			+ gab.flourStrictlyBoundWater / BOUND_WATER_DENSITY
			+ (gab.flourActiveWater + recipe.waterRatio()) / WATER_DENSITY
			+ recipe.saltRatio() / SALT_DENSITY + recipe.oilRatio() / OIL_DENSITY;
	}

	/**
	 * Reference gas capacity [ml/g_flour] from which the proofing target of a product (underproof threshold) is derived.
	 * <p>
	 * It deliberately ignores flour strength, mixer, ash, lipids and ambient humidity: the volume a product must reach
	 * is a property of the product and of the dough volume, so a stronger flour or a better mixer must not make the
	 * target harder to reach.
	 * </p>
	 */
	public static double calculateReferenceGasPotential(final SimulationInputs in, final BakeryProduct targetProduct){
		return doughVolume(in) * targetProduct.getGlutenTearingLimit() * EXPANSION_SCALE * REFERENCE_CAPACITY_FACTOR;
	}

	/**
	 * Gas the gluten network can hold before it tears [ml/g_flour]: it drives the overproof (gas leakage) side of the
	 * fitness and the retention efficiency inside the ODE. It grows with the flour strength and the network
	 * development, so a stronger flour or a mixer widens the safe window without moving the proofing target.
	 */
	public static double calculateMaxGasPotential(final double[] summary, final SimulationInputs in,
			final BakeryProduct targetProduct){
		final double blendW = summary[0];
		final double blendPL = summary[1];
		final double blendFat = summary[2];
		final double blendAsh = summary[3];
		final double oilRatio = in.recipe().oilRatio();

		// 1. An unbalanced P/L (> 0.6) stiffens the dough, reducing its ability to extend without tearing.
		final double plModifier = (blendPL <= 0.5? 1.: Math.max(0.4, 1. - (blendPL - 0.5) * 0.8));

		// 2. Ashes (bran) mechanically cut the gluten network, lowering the structural strength
		final double ashModifier = Math.max(0.5, 1. - (blendAsh * 12.));

		// 3. The fats (recipe oil + flour lipids) within certain limits (up to 8%) lubricate the macro-bubbles expanding the mechanical potential
		final double totalLipids = oilRatio + blendFat;
		final double lipidModifier = 1. + (totalLipids <= 0.08? totalLipids * 1.5: 0.12 - (totalLipids - 0.08) * 2.);

		// 4. Gluten network strength saturates with W: a stronger flour holds more gas, with diminishing returns
		final double strengthFactor = blendW / (blendW + W_HALF_SATURATION);

		// 5. Network development reached by kneading alone (folds add to it during the simulation), saturating at 1
		final double kneadingEfficiency = in.kneading().type().getDevelopmentEfficiency();
		final double developmentFactor = KNEADING_FACTOR_FLOOR
			+ (1. - KNEADING_FACTOR_FLOOR) * Math.min(1., kneadingEfficiency);

		return doughVolume(in) * targetProduct.getGlutenTearingLimit() * EXPANSION_SCALE * strengthFactor
			* developmentFactor * plModifier * ashModifier * lipidModifier;
	}

	public static double[] inputsSummary(final SimulationInputs in){
		double blendW = 0;
		double blendFat = 0;
		double blendAsh = 0;
		double sumProductLnPL = 0.;

		final double[] fr = in.fractions();
		final FlourInput[] mx = in.flourMatrix();

		for(int i = 0; i < mx.length; i ++){
			blendW += fr[i] * mx[i].strength();
			blendFat += fr[i] * mx[i].fat();
			blendAsh += fr[i] * mx[i].ash();

			final double pl = mx[i].plRatio();
			if(pl > 0)
				sumProductLnPL += fr[i] * Math.log(pl);
		}
		final double blendPL = Math.exp(sumProductLnPL);

		return new double[]{blendW, blendPL, blendFat, blendAsh};
	}

}
