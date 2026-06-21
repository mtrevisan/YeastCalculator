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
 */
public class YeastOptimizer{

	public static double findOptimalYeast(final SimulationInputs in, final BakeryProduct targetProduct){
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in);
		final double baseMaxGasPotential = calculateMaxGasPotential(inputsSummary(in), in, targetProduct);

		final UnivariateFunction objective = yeastAttempt -> {
			final double[] y = runSimulation(in, gab, baseMaxGasPotential, yeastAttempt);
			if(y == null)
				return Double.MAX_VALUE;

			double currentTime = 0.;
			for(final StageInput stage : in.getStages())
				currentTime += stage.getDuration();

			final FoldEventHandler foldHandler = new FoldEventHandler(in.getFolds());
			final double[] yDotFinal = new double[4];
			final StageInput lastStage = in.getStages()[in.getStages().length - 1];
			final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.getTemperature(), lastStage.getRelativeHumidity(),
				gab.flourActiveWater, baseMaxGasPotential, in, foldHandler);
			finalOde.computeDerivatives(currentTime, y, yDotFinal);

			final double skinningModifier = (lastStage.getRelativeHumidity() >= 0.70
				? 1.
				: Math.max(0.6, 1. - (0.70 - lastStage.getRelativeHumidity()) * 0.8));
			final double dynamicMaxPotential = baseMaxGasPotential * Math.pow(1.05, in.getFolds().length) * skinningModifier;

			// CLEAN ENCAPSULATION: Delegates the specific target cost to the selected product enum
			return targetProduct.computeFitness(y[2], dynamicMaxPotential, y[1] ,y[0], yDotFinal[2]);
		};

		final BrentOptimizer optimizer = new BrentOptimizer(1.e-6, 1.e-6);
		return optimizer.optimize(
			new MaxEval(200),
			new UnivariateObjectiveFunction(objective),
			GoalType.MINIMIZE,
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

	private static double[] runSimulation(final SimulationInputs in, final GabMoistureModel.GabResult gab,
			final double baseMaxGasPotential, final double yeastRatio){
		// ACTIVATION YeastInput.getRehydrationDurationHours
		// Calculation of biological efficiency based on rehydration time
		// Optimal window assumed to be 10 minutes (0.166 hours). Below or above this window results in loss of cell viability.
		final double rehydrHours = in.getYeastProperties().getRehydrationDurationHours();
		double rehydrationEfficiencyModifier = 1.;
		if(rehydrHours < 0.16)
			// Penalty for failure to activate
			rehydrationEfficiencyModifier = 0.7 + (rehydrHours / 0.16) * 0.3;
		else if(rehydrHours > 0.5)
			// Autolysis/Early Starvation in Water
			rehydrationEfficiencyModifier = Math.max(0.5, 1. - (rehydrHours - 0.5) * 0.4);

		final double initialX = yeastRatio * (1. - in.getYeastProperties().getYeastMoisture()) * rehydrationEfficiencyModifier;

		double totalFlourSugar = 0.;
		for(int i = 0; i < in.getFlourMatrix().length; i ++)
			totalFlourSugar += in.getFractions()[i] * in.getFlourMatrix()[i].getSugar();
		final double initialS = totalFlourSugar + (in.getRecipe().getMaltRatio() * in.getRecipe().getMaltSugarContent());
		final double initialV = 0.;
		final double initialEtOH = 0.;
		final double[] y = new double[]{initialX, initialS, initialV, initialEtOH};
		double currentTime = 0.;

		final FoldEventHandler foldHandler = new FoldEventHandler(in.getFolds());

		for(final StageInput stage : in.getStages()){
			final DoughOdeSystem ode = new DoughOdeSystem(stage.getTemperature(), stage.getRelativeHumidity(),
				gab.flourActiveWater, baseMaxGasPotential, in, foldHandler);
			final DormandPrince853Integrator integrator = new DormandPrince853Integrator(1.e-4, 0.1,
				1.e-5, 1.e-5);
			integrator.addEventHandler(foldHandler, 0.01, 1.e-4, 100);

			final double stageEndTime = currentTime + stage.getDuration();
			try{
				integrator.integrate(ode, currentTime, y, stageEndTime, y);
			}
			catch(final Exception e){
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
		final double oilRatio = in.getRecipe().getOilRatio();

		// 1. An unbalanced P/L (> 0.6) stiffens the dough, reducing its ability to extend without tearing.
		final double plModifier = (blendPL <= 0.5? 1.: Math.max(0.4, 1. - (blendPL - 0.5) * 0.8));

		// 2. Ashes (bran) mechanically cut the gluten network, lowering the structural strength
		final double ashModifier = Math.max(0.5, 1. - (blendAsh * 12.));

		// 3. The fats (recipe oil + flour lipids) within certain limits (up to 8%) lubricate the macro-bubbles expanding the mechanical potential
		final double totalLipids = oilRatio + blendFat;
		final double lipidModifier = 1. + (totalLipids <= 0.08? totalLipids * 1.5: 0.12 - (totalLipids - 0.08) * 2.);

		final double kneadingEfficiency = in.getKneading().getType().getDevelopmentEfficiency();

		return blendW * targetProduct.getGlutenTearingLimit() * 0.05 * kneadingEfficiency * plModifier * ashModifier
			* lipidModifier;
	}

	private static double[] inputsSummary(final SimulationInputs in){
		double blendW = 0;
		double blendFat = 0;
		double blendAsh = 0;
		double sumProductLnPL = 0.;

		final double[] fr = in.getFractions();
		final FlourInput[] mx = in.getFlourMatrix();

		for(int i = 0; i < mx.length; i ++){
			blendW += fr[i] * mx[i].getStrength();
			blendFat += fr[i] * mx[i].getFat();
			blendAsh += fr[i] * mx[i].getAsh();

			final double pl = mx[i].getPlRatio();
			if(pl > 0)
				sumProductLnPL += fr[i] * Math.log(pl);
		}
		final double blendPL = Math.exp(sumProductLnPL);

		return new double[]{blendW, blendPL, blendFat, blendAsh};
	}

}
