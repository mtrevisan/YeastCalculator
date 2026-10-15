package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.domain.YeastInput;
import io.github.mtrevisan.yeastcalculator.plugin.PrefermentPlugin;
import io.github.mtrevisan.yeastcalculator.simulation.DoughOdeSystem;
import io.github.mtrevisan.yeastcalculator.simulation.FoldEventHandler;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;


/**
 * Conducts bounded search algorithms to match yeast velocity targets to structural volume plateaus.
 * Dynamically computes rehydration limits based on current cell moisture characteristics.
 */
public class YeastOptimizer{

	/**
	 * Diagnostic utility to return the exact calculated optimal rehydration window bounds (in minutes).
	 * @return An array containing [MinOptimalMinutes, MaxOptimalMinutes]
	 */
	public static double[] calculateOptimalRehydrationWindow(final double yeastMoisture){
		if(yeastMoisture >= 0.65)
			// Fresh yeast block: cells are already active. Optimal use is immediate (0 mins).
			// Starvation stress triggers if floating in pure water without substrates for more than 15 mins.
			return new double[]{0.0, 15.0};
		else
			// Dry mummified yeast: requires strict structural re-swelling time to protect cell walls.
			// Optimal window is exactly between 9.6 minutes (0.16h) and 30 minutes (0.50h).
			return new double[]{0.16 * 60.0, 0.50 * 60.0};
	}

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
		final UnivariateFunction objective = yeastAttempt -> {
			// Run the ODE simulation using the current yeast attempt value
			final double[] y = runSimulation(runtimeInputs, gab, baseMaxGasPotential, yeastAttempt, preConditioner);
			// If the simulation crashed or hit a mathematical wall, return a maximum penalty score
			if(y == null)
				return Double.MAX_VALUE;

			return evaluateFitness(runtimeInputs, gab, baseMaxGasPotential, y, targetProduct);
		};

		final BrentOptimizer optimizer = new BrentOptimizer(1.e-6, 1.e-6);
		return optimizer.optimize(
			new MaxEval(200),
			new UnivariateObjectiveFunction(objective),
			GoalType.MINIMIZE,
			// Scans from 0.05% up to 3.5% fresh yeast fractions
			new SearchInterval(0.0005, 0.035)
		).getPoint();
	}

	/**
	 * Scores a final simulated state against the target product.
	 * <p>
	 * The end-of-timeline gas retention ceiling is rebuilt exactly as the integration saw it: the structural gain of
	 * every fold applied within the timeline and the surface skinning of the last stage.
	 * </p>
	 */
	private static double evaluateFitness(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double baseMaxGasPotential, final double[] y, final BakeryProduct targetProduct){
		// Calculate total timeline execution [min]
		double totalMinutes = 0.;
		for(final StageInput stage : in.stages())
			totalMinutes += stage.duration();
		totalMinutes *= 60.;

		// Only folds that actually occur within the timeline boost the final ceiling
		final FoldEventHandler foldHandler = new FoldEventHandler(in.folds());
		foldHandler.advanceFoldsUpTo(totalMinutes);
		final int appliedFolds = foldHandler.getAppliedFolds();

		final double[] yDotFinal = new double[4];
		// Final gas production speed: the fold gain is already in the potential, so use an exhausted fold handler
		final StageInput lastStage = in.stages()[in.stages().length - 1];
		final double lastStageRH = lastStage.relativeHumidity();
		final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.temperature(), lastStageRH, gab,
			baseMaxGasPotential, in, foldHandler);
		finalOde.computeDerivatives(totalMinutes, y, yDotFinal);

		final double skinningModifier = DoughOdeSystem.skinningModifier(lastStageRH);
		final double dynamicMaxPotential = baseMaxGasPotential
			* Math.pow(FoldEventHandler.STRUCTURAL_GAIN, appliedFolds)
			* skinningModifier;

		// Pass oilRatio and appliedFolds so the new reward terms can act on them
		return targetProduct.computeFitness(y[2], dynamicMaxPotential, y[1], y[0], yDotFinal[2],
			in.recipe().oilRatio(), appliedFolds, FoldEventHandler.DEGASSING_RETENTION);
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
		final double rehydrationWindow = (preConditioner == null? yeastInput.getRehydrationWindow() / 60.: 0.);
		final double yeastMoisture = yeastInput.yeastMoisture();
		double rehydrationEfficiencyModifier = 1.;
		if(yeastMoisture >= 0.65){
			// --- FRESH YEAST PATHWAY (e.g., 70% moisture panetto) ---
			// Cells are awake and active. If left floating in pure water without substrates for too long,
			// autolysis and early starvation acceleration triggers.
			// More than 15 minutes in pure water solvent matrix
			if(rehydrationWindow > 0.25)
				rehydrationEfficiencyModifier = Math.max(0.4, 1. - (rehydrationWindow - 0.25) * 0.8);
		}
		else{
			// --- DEHYDRATED DRY YEAST PATHWAY (e.g., active dry yeast) ---
			// Structural mummified cells require a rigid timeline to safely re-swell membrane proteins.
			// Approx 9.6 minutes threshold limit
			final double minRequiredRehydrationWindow = 0.16;
			if(rehydrationWindow < minRequiredRehydrationWindow)
				// Penalty for dry cell shear trauma
				rehydrationEfficiencyModifier = 0.6 + (rehydrationWindow / minRequiredRehydrationWindow) * 0.4;
			else if(rehydrationWindow > 0.50)
				// Starvation envelope after 30 mins
				rehydrationEfficiencyModifier = Math.max(0.5, 1. - (rehydrationWindow - 0.50) * 0.4);
		}

		// Set initial conditions for state variables: [Biomass X, Sugars S, Gas Volume V, Ethanol EtOH]
		final double initialX = yeastRatio * (1. - yeastMoisture) * rehydrationEfficiencyModifier;

		double totalFlourSugar = 0.;
		for(int i = 0; i < in.flourMatrix().length; i ++)
			totalFlourSugar += in.fractions()[i] * in.flourMatrix()[i].sugar();
		final double initialS = totalFlourSugar + (in.recipe().maltRatio() * in.recipe().maltSugarContent());
		final double initialV = 0.;
		final double initialEtOH = 0.;
		// Build the direct dough baseline state array vector
		final double[] y = new double[]{initialX, initialS, initialV, initialEtOH};

		// --- INTERCEPTOR TRIGGER FIX ---
		// If an active plugin preconditioner hook is provided, let it execute now!
		// It intercepts y and modifies biomass, sugars, and ethanol based on Phase 1 outcomes.
		if(preConditioner != null)
			preConditioner.applyToInitialState(y, totalFlourSugar);

		double currentTime = 0.;
		// Create the event handler to intercept stretch and fold timestamps
		final FoldEventHandler foldHandler = new FoldEventHandler(in.folds());
		for(final StageInput stage : in.stages()){
			final DoughOdeSystem ode = new DoughOdeSystem(stage.temperature(), stage.relativeHumidity(),
				gab, baseMaxGasPotential, in, foldHandler);
			// Variable-step size integrator setup
			final DormandPrince853Integrator integrator = new DormandPrince853Integrator(1.e-4, 0.1,
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

	public static double calculateMaxGasPotential(final double[] summary, final SimulationInputs in,
			final BakeryProduct targetProduct){
		// ACTIVATE FlourInput.getPlRatio, getFat, getAsh and DoughRecipe.getOilRatio
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

		final double kneadingEfficiency = in.kneading().type().getDevelopmentEfficiency();

		return blendW * targetProduct.getGlutenTearingLimit() * 0.05 * kneadingEfficiency * plModifier
			* ashModifier * lipidModifier;
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
