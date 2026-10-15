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


public class GlobalMultivariateRecipeOptimizer{

	// Trust region radii, expressed in the normalized [0, 1] search space
	private static final double INITIAL_TRUST_REGION_RADIUS = 0.2;
	private static final double STOPPING_TRUST_REGION_RADIUS = 1.e-4;
	private static final int MAX_EVALUATIONS = 2000;

	// Number of continuous variables per inner optimization: [Hydration, Salt, Oil, Malt,
	// InitialFoldDelay, RelaxationFactor]. The fold count is handled by an outer loop.
	private static final int NUM_CONTINUOUS_VARIABLES = 6;

	// Exhaustive enumeration of the fold count. BOBYQA cannot optimize a discrete variable
	// because Math.round() makes the objective non-smooth; a local quadratic model built on
	// a flat plateau cannot generate a descent direction. The outer loop sidesteps the issue.
	private static final int MAX_FOLD_SETS = 5;


	public record OptimizedRecipeResult(DoughRecipe recipe, double optimalYeastRatio, double finalFitness,
		int optimizedFoldCount, double optimizedFoldInitialDelay, double optimizedRelaxationFactor){}


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
		// Adjust absorption cell ceiling slightly based on product style tolerances (e.g. Pan pizza handles more water than freestanding loaves)
		double dynamicHydrationCap = baselineWac;
		if(targetProduct == BakeryProduct.ROMAN_PAN_PIZZA)
			// High aeration pan styles can push boundaries via mechanical support
			dynamicHydrationCap += 0.08;
		else if(targetProduct == BakeryProduct.NEAPOLITAN_PIZZA)
			// Rigid upper bounds limit to safeguard hand stretching structural integrity
			dynamicHydrationCap = Math.min(dynamicHydrationCap, 0.68);

		// Hard physical constraints boundaries clamp
		dynamicHydrationCap = Math.clamp(dynamicHydrationCap, 0.55, 0.88);

		// Diagnostic info track printout
		System.out.printf("Flour Blend WAC resolved. Dynamically clamping Hydration Cap to: %.1f%%\n",
			dynamicHydrationCap * 100.);

		// Continuous variable bounds
		final double[] lowerBounds = new double[]{0.50, 0.015, 0., 0., 15., 1.};
		final double[] upperBounds = new double[]{dynamicHydrationCap, 0.030, 0.05, 0.015, 60., 1.5};

		// Baseline initial guess for the continuous variables
		final double safeInitialHydration = Math.min(0.62, dynamicHydrationCap - 0.02);
		final double[] initialGuess = new double[]{safeInitialHydration, 0.022, 0.02, 0.004, 20., 1.2};

		// Execute multivariate optimization execution pass.
		// The search runs in the normalized [0, 1] space.
		final double[] normalizedLowerBounds = new double[NUM_CONTINUOUS_VARIABLES];
		final double[] normalizedUpperBounds = new double[NUM_CONTINUOUS_VARIABLES];
		Arrays.fill(normalizedUpperBounds, 1.);

		OptimizedRecipeResult bestResult = null;
		double previousFitness = Double.MAX_VALUE;

// Outer loop over the discrete fold count. The fold reward is a Gaussian centered on a
// product-specific optimum, so the fitness as a function of the fold count is unimodal
// (it has a single minimum). Once the fitness starts increasing, every subsequent fold
// count is guaranteed to be worse, so the loop can stop early without losing the optimum.
		System.out.println("Scanning fold counts 0 to " + MAX_FOLD_SETS + "...");
		for(int folds = 0; folds <= MAX_FOLD_SETS; folds ++){
			final FullRecipeMultivariateObjective objective = new FullRecipeMultivariateObjective(baseInputs,
				targetProduct, lowerBounds, upperBounds, folds);

			final BOBYQAOptimizer optimizer = new BOBYQAOptimizer(2 * NUM_CONTINUOUS_VARIABLES + 1,
				INITIAL_TRUST_REGION_RADIUS, STOPPING_TRUST_REGION_RADIUS);

			final PointValuePair optimum = optimizer.optimize(new MaxEval(MAX_EVALUATIONS),
				new ObjectiveFunction(objective), GoalType.MINIMIZE,
				new InitialGuess(objective.normalize(initialGuess)),
				new SimpleBounds(normalizedLowerBounds, normalizedUpperBounds));

			final double currentFitness = optimum.getValue();
			final double[] p = objective.denormalize(optimum.getPoint());

			System.out.printf("  folds=%d -> fitness=%.4f (hydration=%.1f%%, salt=%.2f%%, oil=%.2f%%, malt=%.2f%%)\n",
				folds, currentFitness, p[0] * 100., p[1] * 100., p[2] * 100., p[3] * 100.);

			// Keep the best result seen so far
			if(bestResult == null || currentFitness < bestResult.finalFitness()){
				final DoughRecipe candidateRecipe = new DoughRecipe(p[0], p[1], p[2], p[3],
					baseInputs.recipe().maltSugarContent(), baseInputs.recipe().maltPollakUnit(),
					baseInputs.recipe().mixerFrictionFactor());

				final double candidateDelay = p[4];
				final double candidateRelaxation = p[5];
				final double[] candidateFolds = FullRecipeMultivariateObjective.buildFoldSchedule(folds,
					candidateDelay, candidateRelaxation);

				final SimulationInputs finalInputs = new SimulationInputs(baseInputs.fractions(),
					baseInputs.flourMatrix(), baseInputs.flourTemperature(), baseInputs.airRelativeHumidity(),
					baseInputs.yeastProperties(), candidateRecipe, baseInputs.kneading(), baseInputs.stages(),
					candidateFolds);

				final double bestYeast = YeastOptimizer.findOptimalYeast(finalInputs, targetProduct, null);

				bestResult = new OptimizedRecipeResult(candidateRecipe, bestYeast, currentFitness, folds,
					candidateDelay, candidateRelaxation);
			}

			// Early termination: once the fitness has increased relative to the previous fold count,
			// and we have already passed the optimum, every subsequent fold count will be worse.
			// The Gaussian reward makes the fold-fitness curve symmetric, so this is safe.
			if(currentFitness > previousFitness + 1.e-3)
				break;

			previousFitness = currentFitness;
		}

		return bestResult;
	}

}
