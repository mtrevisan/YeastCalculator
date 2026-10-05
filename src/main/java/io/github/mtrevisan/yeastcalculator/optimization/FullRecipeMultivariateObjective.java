package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.simulation.DoughOdeSystem;
import io.github.mtrevisan.yeastcalculator.simulation.FoldEventHandler;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
import org.apache.commons.math3.analysis.MultivariateFunction;


/**
 * True multivariate objective function that maps Apache's continuous search vector
 * directly into the synchronized recipe space.
 */
public class FullRecipeMultivariateObjective implements MultivariateFunction{

	private final SimulationInputs baseInputs;
	private final BakeryProduct targetProduct;


	public FullRecipeMultivariateObjective(final SimulationInputs baseInputs, final BakeryProduct targetProduct){
		this.baseInputs = baseInputs;
		this.targetProduct = targetProduct;
	}

	@Override
	public double value(final double[] point){
		// Extract simultaneous multi-variable parameters from the optimizer vector
		final double hydration = point[0];
		final double salt = point[1];
		final double oil = point[2];
		final double malt = point[3];

		// Construct the candidate recipe payload
		final DoughRecipe candidateRecipe = new DoughRecipe(hydration, salt, oil, malt,
			baseInputs.getRecipe().getMaltSugarContent(), baseInputs.getRecipe().getMaltPollakUnit(),
			baseInputs.getRecipe().getMixerFrictionFactor());

		final SimulationInputs candidateInputs = new SimulationInputs(baseInputs.getFractions(),
			baseInputs.getFlourMatrix(), baseInputs.getFlourTemperature(), baseInputs.getAirRelativeHumidity(),
			baseInputs.getYeastProperties(), candidateRecipe, baseInputs.getKneading(), baseInputs.getStages(),
			baseInputs.getFolds());

		// 1. Solve the nested univariate problem for the optimal yeast target fraction
		final double optimalYeast = YeastOptimizer.findOptimalYeast(candidateInputs, targetProduct);

		// 2. Execute the full system ODE integration loop
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(candidateInputs);
		final double[] finalState = YeastOptimizer.runSimulation(candidateInputs, gab, optimalYeast, targetProduct);

		if(finalState == null)
			// Return maximum penalty if numerical system crashes
			return Double.MAX_VALUE;

		// 3. Reconstruct the dynamic Bloksma rheological constraints at the timeline endpoint
		final double[] summary = YeastOptimizer.inputsSummary(candidateInputs);
		final double baseMaxGasPotential = YeastOptimizer.calculateMaxGasPotential(summary, candidateInputs,
			targetProduct);
		final StageInput lastStage = baseInputs.getStages()[baseInputs.getStages().length - 1];

		final double skinningModifier = (lastStage.getRelativeHumidity() >= 0.70
			? 1.
			: Math.max(0.6, 1. - (0.70 - lastStage.getRelativeHumidity()) * 0.8));
		final double dynamicMaxPotential = baseMaxGasPotential * Math.pow(1.05, baseInputs.getFolds().length)
			* skinningModifier;

		double[] yDotFinal = new double[4];
		double currentTime = 0.;
		for(StageInput stage : baseInputs.getStages())
			currentTime += stage.getDuration();

		final FoldEventHandler foldHandler = new FoldEventHandler(baseInputs.getFolds());
		final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.getTemperature(), lastStage.getRelativeHumidity(),
			gab.flourActiveWater, baseMaxGasPotential, candidateInputs, foldHandler);
		finalOde.computeDerivatives(currentTime * 60., finalState, yDotFinal);

		// Evaluate and return global multi-variable fitness cost directly from the enum
		return targetProduct.computeFitness(finalState[2], dynamicMaxPotential, finalState[1], finalState[0], yDotFinal[2]);
	}

}
