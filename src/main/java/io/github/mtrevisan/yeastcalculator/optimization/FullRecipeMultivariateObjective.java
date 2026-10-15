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
 * directly into the synchronized recipe space with [0, 1] input normalization.
 */
public class FullRecipeMultivariateObjective implements MultivariateFunction{

	// The next fold interval is physically scaled longer to account for gluten relaxation
	public static final double RELAXATION_FACTOR = 1.4;


	private final SimulationInputs baseInputs;
	private final BakeryProduct targetProduct;
	private final double[] lowerBounds;
	private final double[] upperBounds;


	public FullRecipeMultivariateObjective(final SimulationInputs baseInputs, final BakeryProduct targetProduct,
			final double[] lowerBounds, final double[] upperBounds){
		this.baseInputs = baseInputs;
		this.targetProduct = targetProduct;
		this.lowerBounds = lowerBounds;
		this.upperBounds = upperBounds;
	}


	/**
	 * Maps a physical point onto the normalized [0, 1] search space.
	 */
	public double[] normalize(final double[] point){
		final int length = point.length;
		final double[] normalizedPoint = new double[length];
		for(int i = 0; i < length; i ++)
			normalizedPoint[i] = (point[i] - lowerBounds[i]) / (upperBounds[i] - lowerBounds[i]);
		return normalizedPoint;
	}

	/**
	 * Maps a normalized [0, 1] search vector onto its physical bounds.
	 */
	private double denormalize(final double normalizedValue, final int index){
		return lowerBounds[index] + normalizedValue * (upperBounds[index] - lowerBounds[index]);
	}

	@Override
	public double value(final double[] point){
		// Denormalize inputs from [0, 1] to physical ranges
		final double hydration = denormalize(point[0], 0);
		final double salt = denormalize(point[1], 1);
		final double oil = denormalize(point[2], 2);
		final double malt = denormalize(point[3], 3);
		final int numberOfFolds = (int)Math.round(denormalize(point[4], 4));
		final double initialFoldDelay = denormalize(point[5], 5);
		final double relaxationFactor = denormalize(point[6], 6);

		// Dynamically build the folds timestamp array in minutes
		final double[] dynamicFolds = new double[numberOfFolds];
		double localizedDelay = initialFoldDelay;
		double cumulativeTime = 0.;
		for(int i = 0; i < numberOfFolds; i ++){
			cumulativeTime += localizedDelay;
			dynamicFolds[i] = cumulativeTime;

			// The next interval is physically scaled by the optimized relaxation factor
			localizedDelay *= relaxationFactor;
		}

		// Construct the candidate recipe payload
		final DoughRecipe candidateRecipe = new DoughRecipe(hydration, salt, oil, malt,
			baseInputs.recipe().maltSugarContent(), baseInputs.recipe().maltPollakUnit(),
			baseInputs.recipe().mixerFrictionFactor());

		final SimulationInputs candidateInputs = new SimulationInputs(baseInputs.fractions(),
			baseInputs.flourMatrix(), baseInputs.flourTemperature(), baseInputs.airRelativeHumidity(),
			baseInputs.yeastProperties(), candidateRecipe, baseInputs.kneading(), baseInputs.stages(),
			dynamicFolds);

		// 1. Solve the nested univariate problem for the optimal yeast target fraction
		final double optimalYeast = YeastOptimizer.findOptimalYeast(candidateInputs, targetProduct, null);

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
		final StageInput lastStage = baseInputs.stages()[baseInputs.stages().length - 1];

		final double skinningModifier = (lastStage.relativeHumidity() >= 0.70
			? 1.
			: Math.max(0.6, 1. - (0.70 - lastStage.relativeHumidity()) * 0.8));
		final double dynamicMaxPotential = baseMaxGasPotential * Math.pow(1.05, baseInputs.folds().length)
			* skinningModifier;

		double[] yDotFinal = new double[4];
		double currentTime = 0.;
		for(StageInput stage : baseInputs.stages())
			currentTime += stage.duration();

		final FoldEventHandler foldHandler = new FoldEventHandler(baseInputs.folds());
		final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.temperature(), lastStage.relativeHumidity(),
			gab, baseMaxGasPotential, candidateInputs, foldHandler);
		finalOde.computeDerivatives(currentTime * 60., finalState, yDotFinal);

		// Evaluate and return global multi-variable fitness cost directly from the enum
		return targetProduct.computeFitness(finalState[2], dynamicMaxPotential, finalState[1], finalState[0], yDotFinal[2]);
	}

	public static double[] buildFoldSchedule(final int optimizedFoldCount, double currentDelay){
		final double[] optimizedFoldsArray = new double[optimizedFoldCount];
		// Serves as initial fold delay field variable mapping
		double cumulativeMinutes = 0.;
		for(int i = 0; i < optimizedFoldCount; i ++){
			cumulativeMinutes += currentDelay;
			optimizedFoldsArray[i] = cumulativeMinutes;
			currentDelay *= RELAXATION_FACTOR;
		}
		return optimizedFoldsArray;
	}

}
