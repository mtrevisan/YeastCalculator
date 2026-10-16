package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.PointValuePair;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;


public class GlobalMultivariateRecipeOptimizer{

	public record OptimizedRecipeResult(DoughRecipe recipe, double optimalYeastRatio, double finalFitness,
		int optimizedFoldCount, double optimizedFoldInterval){}


	// Maximum number of stretch & fold sets enumerated by the search
	private static final int MAX_FOLDS = 5;
	// Trust region radii, expressed in the normalized [0, 1] search space
	private static final double INITIAL_TRUST_REGION_RADIUS = 0.2;
	private static final double STOPPING_TRUST_REGION_RADIUS = 1.e-4;
	private static final int MAX_EVALUATIONS = 300;
	// Minimum width of the hydration search interval [g_water/g_flour]
	private static final double MIN_HYDRATION_SEARCH_WIDTH = 0.03;


	public static OptimizedRecipeResult optimizeFullRecipe(final SimulationInputs baseInputs,
			final BakeryProduct targetProduct){
		// --- AUTOMATIC HYDRATION CAP HANDLER (BIOCHEMICAL WAC MODEL) ---
		// Resolves the weighted average characteristics of the flour blend matrix
		double blendedProtein = 0.;
		double blendedFiber = 0.;
		double blendedStrengthW = 0.;
		final double[] fractions = baseInputs.fractions();
		final FlourInput[] flours = baseInputs.flourMatrix();

		for(int i = 0; i < flours.length; i ++){
			blendedProtein += fractions[i] * flours[i].protein();
			blendedFiber += fractions[i] * flours[i].fiber();
			blendedStrengthW += fractions[i] * flours[i].strength();
		}

		// Estimated damaged starch factor scaling derived from flour strength W
		final double estimatedDamagedStarch = 0.02 + blendedStrengthW / 15000.;
		// Predictive Water Absorption Capacity (WAC) baseline linear regression formula:
		// Proteins absorb ~1.3x their dry weight, fibers absorb ~2.1x, damaged starch absorbs ~2.5x.
		final double baselineWac = 0.46 + 1.35 * blendedProtein + 2.10 * blendedFiber + 2.50 * estimatedDamagedStarch;
		// Container-baked styles can carry extra water thanks to the mechanical support of the pan walls
		double dynamicHydrationCap = baselineWac + targetProduct.getPanSupportHydrationBonus();

		// Hard physical constraints boundaries clamp
		dynamicHydrationCap = Math.clamp(dynamicHydrationCap, 0.55, 0.88);

		// Diagnostic info track printout
		System.out.printf("Flour Blend WAC resolved. Dynamically clamping Hydration Cap to: %.1f%%\n", dynamicHydrationCap * 100.);

		// --- HYDRATION SEARCH RANGE: STYLE RANGE INTERSECTED WITH THE FLOUR CAPACITY ---
		double minHydration = targetProduct.getMinHydration();
		final double maxHydration = Math.min(targetProduct.getMaxHydration(), dynamicHydrationCap);
		if(maxHydration - minHydration < MIN_HYDRATION_SEARCH_WIDTH){
			// The flour cannot absorb the water the style needs: search just below the flour capacity instead
			minHydration = maxHydration - MIN_HYDRATION_SEARCH_WIDTH;
			System.out.printf("[WARNING] Flour blend too weak for %s: style range %.0f%%-%.0f%%, flour capacity %.1f%%.\n",
				targetProduct.name(), targetProduct.getMinHydration() * 100., targetProduct.getMaxHydration() * 100.,
				dynamicHydrationCap * 100.);
		}
		System.out.printf("Hydration search range resolved to: %.1f%% - %.1f%%\n", minHydration * 100., maxHydration * 100.);

		// Search Vector layout: [Hydration, Salt, Oil, Malt, InitialFoldDelayMinutes]
		// Start from the center of the hydration search range
		final double initialHydration = (minHydration + maxHydration) / 2.;
		final double[] initialGuess = new double[]{initialHydration, 0.022, 0.02, 0.002, 20.};
		// Strictly enforce lower and upper physical bounds to maintain culinary and rheological safety
		// Min: style/flour Hydr, 1.8% Salt, 0% Oil, 0% Malt, 15 min first fold delay
		final double[] lowerBounds = new double[]{minHydration, 0.018, 0., 0., 15.};
		// Max: style/flour Hydr, 3.2% Salt, 10% Oil, 1.5% Malt, 60 min first fold delay
		final double[] upperBounds = new double[]{maxHydration, 0.032, 0.1, 0.015, 60.};

		// The fold count is an integer (rounding it inside a continuous search makes the objective a step function that
		// BOBYQA cannot model): run an independent continuous search for each count and keep the best one
		return IntStream.rangeClosed(0, MAX_FOLDS)
			.parallel()
			.mapToObj(numberOfFolds -> optimizeForFoldCount(baseInputs, targetProduct, numberOfFolds, initialGuess,
				lowerBounds, upperBounds))
			.min(Comparator.comparingDouble(OptimizedRecipeResult::finalFitness))
			.orElseThrow();
	}

	private static OptimizedRecipeResult optimizeForFoldCount(final SimulationInputs baseInputs,
			final BakeryProduct targetProduct, final int numberOfFolds, final double[] initialGuess,
			final double[] lowerBounds, final double[] upperBounds){
		final FullRecipeMultivariateObjective objective = new FullRecipeMultivariateObjective(baseInputs, targetProduct,
			numberOfFolds, lowerBounds, upperBounds);

		// BOBYQA requires the number of interpolation points to be mathematically bounded:
		// Standard rule of thumb: 2 * numVariables + 1 points
		final int numVariables = FullRecipeMultivariateObjective.DIMENSION;
		final BOBYQAOptimizer optimizer = new BOBYQAOptimizer(2 * numVariables + 1, INITIAL_TRUST_REGION_RADIUS,
			STOPPING_TRUST_REGION_RADIUS);

		// The search runs in the normalized [0, 1] space
		final double[] normalizedLowerBounds = new double[numVariables];
		final double[] normalizedUpperBounds = new double[numVariables];
		Arrays.fill(normalizedUpperBounds, 1.);

		// Execute multivariate optimization execution pass
		final PointValuePair optimum = optimizer.optimize(new MaxEval(MAX_EVALUATIONS), new ObjectiveFunction(objective),
			GoalType.MINIMIZE, new InitialGuess(objective.normalize(initialGuess)),
			new SimpleBounds(normalizedLowerBounds, normalizedUpperBounds));

		final double[] optimizedPoint = objective.denormalize(optimum.getPoint());
		final SimulationInputs finalInputs = objective.buildCandidateInputs(optimizedPoint);

		// Find the absolute final synchronized optimal yeast dosage for the resolved recipe matrix
		final double bestYeast = YeastOptimizer.findOptimalYeast(finalInputs, targetProduct, null);

		return new OptimizedRecipeResult(finalInputs.recipe(), bestYeast, optimum.getValue(), numberOfFolds,
			optimizedPoint[4]);
	}

}
