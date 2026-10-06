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


public class GlobalMultivariateRecipeOptimizer{

	public record OptimizedRecipeResult(DoughRecipe recipe, double optimalYeastRatio, double finalFitness,
		int optimizedFoldCount, double optimizedFoldInterval){}


	public static OptimizedRecipeResult optimizeFullRecipe(final SimulationInputs baseInputs,
			final BakeryProduct targetProduct){
		// Number of variables to optimize concurrently: [Hydration, Salt, Oil, Malt]
		final int numVariables = 6;

		// BOBYQA requires the number of interpolation points to be mathematically bounded:
		// Standard rule of thumb: 2 * numVariables + 1 = 9 points
		final BOBYQAOptimizer optimizer = new BOBYQAOptimizer(2 * numVariables + 1);

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
		System.out.printf("Flour Blend WAC resolved. Dynamically clamping Hydration Cap to: %.1f%%\n", dynamicHydrationCap * 100.);

		// Search Vector layout: [Hydration, Salt, Oil, Malt, NumberOfFolds, InitialFoldDelayMinutes]
		// Self-correct initial guess if it exceeds our newly calculated ceiling cap
		final double safeInitialHydration = Math.min(0.62, dynamicHydrationCap - 0.02);
		final double[] initialGuess = new double[]{safeInitialHydration, 0.022, 0.02, 0.002, 2., 20.};
		// Strictly enforce lower and upper physical bounds to maintain culinary and rheological safety
		// Min: 50% Hydr, 1.8% Salt, 0% Oil, 0% Malt, Min 0 folds, 15 min apart
		final double[] lowerBounds = new double[]{0.50, 0.018, 0., 0., 0., 15.};
		// Max: 85% Hydr, 3.2% Salt, 10% Oil, 1.5% Malt, Max 5 folds, 60 min apart
		final double[] upperBounds = new double[]{dynamicHydrationCap, 0.032, 0.1, 0.015, 5., 60.};

		final FullRecipeMultivariateObjective objective = new FullRecipeMultivariateObjective(baseInputs, targetProduct);

		// Execute multivariate optimization execution pass
		final PointValuePair optimum = optimizer.optimize(new MaxEval(300), new ObjectiveFunction(objective),
			GoalType.MINIMIZE, new InitialGuess(initialGuess), new SimpleBounds(lowerBounds, upperBounds));

		final double[] optimizedPoints = optimum.getPoint();
		final double finalFitness = optimum.getValue();

		final DoughRecipe finalRecipe = new DoughRecipe(optimizedPoints[0], optimizedPoints[1], optimizedPoints[2],
			optimizedPoints[3], baseInputs.recipe().maltSugarContent(), baseInputs.recipe().maltPollakUnit(),
			baseInputs.recipe().mixerFrictionFactor());

		final int finalFoldsCount = (int)Math.round(optimizedPoints[4]);
		final double finalInitialFoldDelay = optimizedPoints[5];
		final double relaxationFactor = 1.4;
		// Build the verification array using the exact same geometric progression as the objective function
		final double[] finalFolds = new double[finalFoldsCount];
		double localizedDelay = finalInitialFoldDelay;
		double cumulativeTime = 0.;
		for(int i = 0; i < finalFoldsCount; i ++){
			cumulativeTime += localizedDelay;
			finalFolds[i] = cumulativeTime;
			localizedDelay *= relaxationFactor;
		}

		final SimulationInputs finalInputs = new SimulationInputs(baseInputs.fractions(), baseInputs.flourMatrix(),
			baseInputs.flourTemperature(), baseInputs.airRelativeHumidity(), baseInputs.yeastProperties(),
			finalRecipe, baseInputs.kneading(), baseInputs.stages(), finalFolds);

		// Find the absolute final synchronized optimal yeast dosage for the resolved recipe matrix
		final double bestYeast = YeastOptimizer.findOptimalYeast(finalInputs, targetProduct, null);

		return new OptimizedRecipeResult(finalRecipe, bestYeast, finalFitness, finalFoldsCount, finalInitialFoldDelay);
	}

}
