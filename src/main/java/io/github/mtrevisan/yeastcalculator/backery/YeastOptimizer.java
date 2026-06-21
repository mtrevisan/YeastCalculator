package io.github.mtrevisan.yeastcalculator.backery;

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
		final double baseMaxGasPotential = calculateMaxGasPotential(in, targetProduct);

		UnivariateFunction objective = yeastAttempt -> {
			double[] y = runSimulation(in, gab, baseMaxGasPotential, yeastAttempt);
			if(y == null) return Double.MAX_VALUE;

			double currentTime = 0.0;
			for(StageInput stage : in.getStages()){
				currentTime += stage.getDurationHours();
			}

			FoldEventHandler foldHandler = new FoldEventHandler(in.getFolds());
			double[] yDotFinal = new double[3];
			StageInput lastStage = in.getStages()[in.getStages().length - 1];
			DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.getTemperature(), gab.flourActiveWater, baseMaxGasPotential, in, foldHandler);
			finalOde.computeDerivatives(currentTime, y, yDotFinal);

			double residualSugar = y[1];
			double sugarPenalty = residualSugar < targetProduct.getMinSafeSugarThreshold()? 500.0 * (targetProduct.getMinSafeSugarThreshold() - residualSugar): 0.0;

			return Math.abs(yDotFinal[2]) + sugarPenalty;
		};

		BrentOptimizer optimizer = new BrentOptimizer(1e-6, 1e-6);
		return optimizer.optimize(
			new MaxEval(150),
			new UnivariateObjectiveFunction(objective),
			GoalType.MINIMIZE,
			new SearchInterval(0.0005, 0.05)
		).getPoint();
	}

	public static double[] runSimulation(SimulationInputs in, GabMoistureModel.GabResult gab, double baseMaxGasPotential, double yeastRatio){
		double initialX = yeastRatio * (1.0 - in.getYeastProperties().getYeastMoisture());
		double totalFlourSugar = 0.0;
		for(int i = 0; i < in.getFlourMatrix().length; i++){
			totalFlourSugar += in.getFractions()[i] * in.getFlourMatrix()[i].getSugar();
		}
		double initialS = totalFlourSugar + (in.getRecipe().getMaltRatio() * in.getRecipe().getMaltSugarContent());
		double initialV = 0.0;

		double[] y = new double[]{initialX, initialS, initialV};
		double currentTime = 0.0;

		FoldEventHandler foldHandler = new FoldEventHandler(in.getFolds());

		for(StageInput stage : in.getStages()){
			DoughOdeSystem ode = new DoughOdeSystem(stage.getTemperature(), gab.flourActiveWater, baseMaxGasPotential, in, foldHandler);
			DormandPrince853Integrator integrator = new DormandPrince853Integrator(1.0e-4, 0.1, 1.0e-5, 1.0e-5);
			integrator.addEventHandler(foldHandler, 0.01, 1.0e-4, 100);

			double stageEndTime = currentTime + stage.getDurationHours();
			try{
				integrator.integrate(ode, currentTime, y, stageEndTime, y);
			}
			catch(Exception e){
				return null;
			}
			currentTime = stageEndTime;
		}
		return y;
	}

	public static double calculateMaxGasPotential(SimulationInputs in, BakeryProduct targetProduct){
		double totalW = 0.0;
		double[] fractions = in.getFractions();
		FlourInput[] matrix = in.getFlourMatrix();
		for(int i = 0; i < matrix.length; i++){
			totalW += fractions[i] * matrix[i].getStrengthW();
		}
		double kneadingEfficiency = in.getKneading().getType().getDevelopmentEfficiency();
		return totalW * targetProduct.getGlutenTearingLimit() * 0.05 * kneadingEfficiency;
	}

}
