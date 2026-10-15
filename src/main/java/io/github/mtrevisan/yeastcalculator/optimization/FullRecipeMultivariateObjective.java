package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
import org.apache.commons.math3.analysis.MultivariateFunction;


/**
 * Multivariate objective function for a FIXED number of stretch-and-fold sets.
 * <p>
 * The fold count is a discrete variable and cannot be optimized by BOBYQA, which assumes a
 * smooth objective. It is therefore injected at construction time, and the search vector
 * contains only the six continuous variables: [Hydration, Salt, Oil, Malt, InitialDelay,
 * RelaxationFactor]. The caller (GlobalMultivariateRecipeOptimizer) iterates over the fold
 * count externally and compares the resulting fitness values.
 * </p>
 */
public class FullRecipeMultivariateObjective implements MultivariateFunction{

	private final SimulationInputs baseInputs;
	private final BakeryProduct targetProduct;
	private final double[] lowerBounds;
	private final double[] upperBounds;
	private final int numberOfFolds;


	public FullRecipeMultivariateObjective(final SimulationInputs baseInputs, final BakeryProduct targetProduct,
			final double[] lowerBounds, final double[] upperBounds, final int numberOfFolds){
		this.baseInputs = baseInputs;
		this.targetProduct = targetProduct;
		this.lowerBounds = lowerBounds;
		this.upperBounds = upperBounds;
		this.numberOfFolds = numberOfFolds;
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
	public double[] denormalize(final double[] normalizedPoint){
		final int length = normalizedPoint.length;
		final double[] point = new double[length];
		for(int i = 0; i < length; i ++)
			point[i] = lowerBounds[i] + normalizedPoint[i] * (upperBounds[i] - lowerBounds[i]);
		return point;
	}

	/**
	 * Builds the candidate simulation inputs for a physical (not normalized) search point.
	 */
	public SimulationInputs buildCandidateInputs(final double[] point){
		final DoughRecipe candidateRecipe = new DoughRecipe(point[0], point[1], point[2], point[3],
			baseInputs.recipe().maltSugarContent(), baseInputs.recipe().maltPollakUnit(),
			baseInputs.recipe().mixerFrictionFactor());

		return new SimulationInputs(baseInputs.fractions(), baseInputs.flourMatrix(), baseInputs.flourTemperature(),
			baseInputs.airRelativeHumidity(), baseInputs.yeastProperties(), candidateRecipe, baseInputs.kneading(),
			baseInputs.stages(), buildFoldSchedule(numberOfFolds, point[4], point[5]));
	}

	@Override
	public double value(final double[] normalizedPoint){
		final SimulationInputs candidateInputs = buildCandidateInputs(denormalize(normalizedPoint));

		// 1. Solve the nested univariate problem for the optimal yeast target fraction
		final double optimalYeast = YeastOptimizer.findOptimalYeast(candidateInputs, targetProduct, null);

		// 2. Execute the full system ODE integration loop
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(candidateInputs);
		final double[] finalState = YeastOptimizer.runSimulation(candidateInputs, gab, optimalYeast, targetProduct);

		if(finalState == null)
			// Return maximum penalty if numerical system crashes
			return Double.MAX_VALUE;

		// 3. Evaluate the global fitness with the same end-of-timeline constraints used by the yeast optimizer.
		//    Delegating to YeastOptimizer.evaluateFitness guarantees that the outer BOBYQA search and the inner
		//    Brent search see exactly the same fitness landscape, with no risk of divergence.
		final double baseMaxGasPotential = YeastOptimizer.calculateMaxGasPotential(
			YeastOptimizer.inputsSummary(candidateInputs), candidateInputs, targetProduct);
		return YeastOptimizer.evaluateFitness(candidateInputs, gab, baseMaxGasPotential, finalState, targetProduct);
	}

	/**
	 * Builds the fold timestamps [min] as a geometric relaxation curve starting from the first delay.
	 * The relaxation factor is applied to the interval between consecutive folds (not the cumulative time),
	 * so a factor of 1.0 yields uniform intervals and a factor of 1.5 yields progressively longer intervals.
	 */
	public static double[] buildFoldSchedule(final int numberOfFolds, double initialFoldDelay,
			final double relaxationFactor){
		final double[] folds = new double[numberOfFolds];
		// Serves as initial fold delay field variable mapping
		double cumulativeMinutes = 0.;
		for(int i = 0; i < numberOfFolds; i ++){
			cumulativeMinutes += initialFoldDelay;
			folds[i] = cumulativeMinutes;
			initialFoldDelay *= relaxationFactor;
		}
		return folds;
	}

}
