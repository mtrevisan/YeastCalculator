package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
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
	public static double findOptimalYeast(final SimulationInputs in, final BakeryProduct targetProduct){
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in);
		final double baseMaxGasPotential = calculateMaxGasPotential(inputsSummary(in), in, targetProduct);

		// The objective function defines the target for the optimizer to minimize
		final UnivariateFunction objective = yeastAttempt -> {
			// Run the ODE simulation using the current yeast attempt value
			final double[] y = runSimulation(in, gab, baseMaxGasPotential, yeastAttempt);
			// If the simulation crashed or hit a mathematical wall, return a maximum penalty score
			if(y == null)
				return Double.MAX_VALUE;

			// Calculate total timeline execution
			double currentTime = 0.;
			for(final StageInput stage : in.stages())
				currentTime += stage.duration();

			final FoldEventHandler foldHandler = new FoldEventHandler(in.folds());
			final double[] yDotFinal = new double[4];
			final StageInput lastStage = in.stages()[in.stages().length - 1];
			final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.temperature(), lastStage.relativeHumidity(),
				gab.flourActiveWater, baseMaxGasPotential, in, foldHandler);
			finalOde.computeDerivatives(currentTime * 60., y, yDotFinal);

			// Apply skinning and surface adjustments from the final ambient humidity
			final double skinningModifier = (lastStage.relativeHumidity() >= 0.70
				? 1.
				: Math.max(0.6, 1. - (0.70 - lastStage.relativeHumidity()) * 0.8));
			final double dynamicMaxPotential = baseMaxGasPotential * Math.pow(1.05, in.folds().length) * skinningModifier;

			// Delegates the specific cost target directly to your BakeryProduct enum structure
			return targetProduct.computeFitness(y[2], dynamicMaxPotential, y[1] ,y[0], yDotFinal[2]);
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
	 * High-level simulation orchestrator that encapsulates internal structural metrics
	 * to prevent Feature Envy in the calling client.
	 */
	public static double[] runSimulation(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double yeastRatio, final BakeryProduct targetProduct){
		// Internal parameters are resolved automatically inside the owner class
		final double[] summary = inputsSummary(in);
		final double baseMaxGasPotential = calculateMaxGasPotential(summary, in, targetProduct);

		// Delegates to the low-level simulation logic
		return runSimulation(in, gab, baseMaxGasPotential, yeastRatio);
	}

	/**
	 * The core numerical orchestrator. Executes sequential integration blocks via Dormand-Prince 8(5,3).
	 * Calculates the dynamic biological activation modifier as a direct function of cell moisture.
	 */
	private static double[] runSimulation(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double baseMaxGasPotential, final double yeastRatio){
		final double rehydrationDuration = in.yeastProperties().rehydrationDurationHours();
		final double yeastMoisture = in.yeastProperties().yeastMoisture();
		double rehydrationEfficiencyModifier = 1.;
		if(yeastMoisture >= 0.65){
			// --- FRESH YEAST PATHWAY (e.g., 70% moisture panetto) ---
			// Cells are awake and active. If left floating in pure water without substrates for too long,
			// autolysis and early starvation acceleration triggers.
			// More than 15 minutes in pure water solvent matrix
			if(rehydrationDuration > 0.25)
				rehydrationEfficiencyModifier = Math.max(0.4, 1. - (rehydrationDuration - 0.25) * 0.8);
		}
		else{
			// --- DEHYDRATED DRY YEAST PATHWAY (e.g., active dry yeast) ---
			// Structural mummified cells require a rigid timeline to safely re-swell membrane proteins.
			// Approx 9.6 minutes threshold limit
			final double minRequiredHours = 0.16;
			if(rehydrationDuration < minRequiredHours)
				// Penalty for dry cell shear trauma
				rehydrationEfficiencyModifier = 0.6 + (rehydrationDuration / minRequiredHours) * 0.4;
			else if(rehydrationDuration > 0.5)
				// Starvation envelope after 30 mins
				rehydrationEfficiencyModifier = Math.max(0.5, 1. - (rehydrationDuration - 0.50) * 0.4);
		}

		// Set initial conditions for state variables: [Biomass X, Sugars S, Gas Volume V, Ethanol EtOH]
		final double initialX = yeastRatio * (1. - yeastMoisture) * rehydrationEfficiencyModifier;

		double totalFlourSugar = 0.;
		for(int i = 0; i < in.flourMatrix().length; i ++)
			totalFlourSugar += in.fractions()[i] * in.flourMatrix()[i].sugar();
		final double initialS = totalFlourSugar + (in.recipe().maltRatio() * in.recipe().maltSugarContent());
		final double initialV = 0.;
		final double initialEtOH = 0.;
		final double[] y = new double[]{initialX, initialS, initialV, initialEtOH};
		double currentTime = 0.;

		// Create the event handler to intercept stretch and fold timestamps
		final FoldEventHandler foldHandler = new FoldEventHandler(in.folds());

		for(final StageInput stage : in.stages()){
			final DoughOdeSystem ode = new DoughOdeSystem(stage.temperature(), stage.relativeHumidity(),
				gab.flourActiveWater, baseMaxGasPotential, in, foldHandler);
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

		return blendW * targetProduct.getGlutenTearingLimit() * 0.05 * kneadingEfficiency * plModifier * ashModifier
			* lipidModifier;
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

	/**
	 * Calculates the optimal target rehydration duration range based on the moisture profile.
	 *
	 * @return A descriptive string indicating the unpenalized time window.
	 */
	public static String getOptimalRehydrationWindow(final double yeastMoisture){
		if(yeastMoisture >= 0.65)
			return "0 to 15 minutes (Instant mix recommended)";
		else
			return "10 to 30 minutes (Warm-water membrane rehydration mandatory)";
	}

}
