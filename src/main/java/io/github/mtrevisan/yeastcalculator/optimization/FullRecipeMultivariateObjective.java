package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
import org.apache.commons.math3.analysis.MultivariateFunction;


/**
 * True multivariate objective function that maps Apache's continuous search vector
 * directly into the synchronized recipe space.
 * <p>
 * The search vector is normalized: every coordinate lives in [0, 1] and is mapped linearly onto its physical bounds.
 * BOBYQA uses a single trust region radius for all the variables, so without normalization the variables with a
 * wide range (minutes) would never move while those with a narrow range (salt) would dictate the step size.
 * </p>
 * <p>
 * The number of stretch & fold sets is an integer, so it is not part of the continuous search vector: it is fixed
 * per instance and enumerated by the caller.
 * </p>
 */
public class FullRecipeMultivariateObjective implements MultivariateFunction{

	// Search vector layout: [Hydration, Salt, Oil, Malt, InitialFoldDelayMinutes]
	public static final int DIMENSION = 5;

	// The next fold interval is physically scaled longer to account for gluten relaxation
	public static final double RELAXATION_FACTOR = 1.4;

	private final SimulationInputs baseInputs;
	private final BakeryProduct targetProduct;
	private final int numberOfFolds;
	private final double[] lowerBounds;
	private final double[] upperBounds;


	public FullRecipeMultivariateObjective(final SimulationInputs baseInputs, final BakeryProduct targetProduct,
			final int numberOfFolds, final double[] lowerBounds, final double[] upperBounds){
		this.baseInputs = baseInputs;
		this.targetProduct = targetProduct;
		this.numberOfFolds = numberOfFolds;
		this.lowerBounds = lowerBounds.clone();
		this.upperBounds = upperBounds.clone();
	}


	/**
	 * Maps a normalized [0, 1] search vector onto its physical bounds.
	 */
	public double[] denormalize(final double[] normalizedPoint){
		final double[] point = new double[normalizedPoint.length];
		for(int i = 0; i < normalizedPoint.length; i ++)
			point[i] = lowerBounds[i] + normalizedPoint[i] * (upperBounds[i] - lowerBounds[i]);
		return point;
	}

	/**
	 * Maps a physical point onto the normalized [0, 1] search space.
	 */
	public double[] normalize(final double[] point){
		final double[] normalizedPoint = new double[point.length];
		for(int i = 0; i < point.length; i ++)
			normalizedPoint[i] = (point[i] - lowerBounds[i]) / (upperBounds[i] - lowerBounds[i]);
		return normalizedPoint;
	}

	/**
	 * Builds the fold timestamps [min] as a geometric relaxation curve starting from the first delay.
	 */
	public static double[] buildFoldSchedule(final int numberOfFolds, final double initialFoldDelay){
		final double[] folds = new double[numberOfFolds];
		double localizedDelay = initialFoldDelay;
		double cumulativeTime = 0.;
		for(int i = 0; i < numberOfFolds; i ++){
			cumulativeTime += localizedDelay;
			folds[i] = cumulativeTime;
			localizedDelay *= RELAXATION_FACTOR;
		}
		return folds;
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
			baseInputs.stages(), buildFoldSchedule(numberOfFolds, point[4]));
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

		// 3. Evaluate the global fitness with the same end-of-timeline constraints used by the yeast optimizer
		final double baseMaxGasPotential = YeastOptimizer.calculateMaxGasPotential(
			YeastOptimizer.inputsSummary(candidateInputs), candidateInputs, targetProduct);
		return YeastOptimizer.evaluateFitness(candidateInputs, gab, baseMaxGasPotential, finalState, targetProduct);
	}

}
