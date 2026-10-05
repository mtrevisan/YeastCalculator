package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.PointValuePair;


public class GlobalMultivariateRecipeOptimizer{

	public record OptimizedRecipeResult(DoughRecipe recipe, double optimalYeastRatio, double finalFitness){}

	public static OptimizedRecipeResult optimizeFullRecipe(final SimulationInputs baseInputs,
			final BakeryProduct targetProduct){
		// Number of variables to optimize concurrently: [Hydration, Salt, Oil, Malt]
		final int numVariables = 4;

		// BOBYQA requires the number of interpolation points to be mathematically bounded:
		// Standard rule of thumb: 2 * numVariables + 1 = 9 points
		final BOBYQAOptimizer optimizer = new BOBYQAOptimizer(9);

		// Define initial guesses for the simultaneous search engine
		final double[] initialGuess = new double[]{0.62, 0.022, 0.02, 0.002};

		// Strictly enforce lower and upper physical bounds to maintain culinary and rheological safety
		// Min: 50% Hydr, 1.8% Salt, 0% Oil, 0% Malt
		final double[] lowerBounds = new double[]{0.50, 0.018, 0.0, 0.0};
		// Max: 85% Hydr, 3.2% Salt, 10% Oil, 1.5% Malt
		final double[] upperBounds = new double[]{0.85, 0.032, 0.10, 0.015};

		final FullRecipeMultivariateObjective objective = new FullRecipeMultivariateObjective(baseInputs, targetProduct);

		// Execute multivariate optimization execution pass
		final PointValuePair optimum = optimizer.optimize(new MaxEval(300), new ObjectiveFunction(objective),
			GoalType.MINIMIZE, new InitialGuess(initialGuess), new SimpleBounds(lowerBounds, upperBounds));

		final double[] optimizedPoints = optimum.getPoint();
		final double finalFitness = optimum.getValue();

		final DoughRecipe finalRecipe = new DoughRecipe(optimizedPoints[0], optimizedPoints[1], optimizedPoints[2],
			optimizedPoints[3], baseInputs.recipe().maltSugarContent(), baseInputs.recipe().maltPollakUnit(),
			baseInputs.recipe().mixerFrictionFactor());

		final SimulationInputs finalInputs = new SimulationInputs(baseInputs.fractions(), baseInputs.flourMatrix(),
			baseInputs.flourTemperature(), baseInputs.airRelativeHumidity(), baseInputs.yeastProperties(),
			finalRecipe, baseInputs.kneading(), baseInputs.stages(), baseInputs.folds());

		// Find the absolute final synchronized optimal yeast dosage for the resolved recipe matrix
		final double bestYeast = YeastOptimizer.findOptimalYeast(finalInputs, targetProduct);

		return new OptimizedRecipeResult(finalRecipe, bestYeast, finalFitness);
	}

}
