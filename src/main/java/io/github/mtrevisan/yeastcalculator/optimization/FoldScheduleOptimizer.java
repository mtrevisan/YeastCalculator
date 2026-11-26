/**
 * Copyright (c) 2026 Mauro Trevisan
 * <p>
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.stream.IntStream;


/**
 * Optimizes the stretch-and-fold schedule (count, first delay, relaxation factor) and the yeast amount for a recipe
 * whose ingredients are fixed. Counterpart of {@link GlobalMultivariateRecipeOptimizer}, where the recipe is free.
 */
public class FoldScheduleOptimizer{

	private static final double INITIAL_TRUST_REGION_RADIUS = 0.2;
	private static final double STOPPING_TRUST_REGION_RADIUS = 1.e-4;
	private static final int MAX_EVALUATIONS = 1000;

	// [InitialFoldDelay, RelaxationFactor]
	private static final int NUM_CONTINUOUS_VARIABLES = 2;
	private static final double[] LOWER_BOUNDS = {15., 1.};
	private static final double[] UPPER_BOUNDS = {60., 1.5};
	private static final double[] INITIAL_GUESS = {20., 1.2};

	private static final int MAX_FOLD_SETS = 8;


	private record FoldCandidate(int folds, double fitness, double[] point){}

	public record OptimizedFoldResult(double[] folds, double optimalYeastRatio, double finalFitness){}


	public static OptimizedFoldResult optimizeFolds(final SimulationInputs baseInputs, final BakeryProduct targetProduct){
		final List<FoldCandidate> candidates = IntStream.rangeClosed(0, MAX_FOLD_SETS)
			.parallel()
			.mapToObj(folds -> optimizeForFoldCount(baseInputs, targetProduct, folds))
			.filter(Objects::nonNull)
			.sorted(Comparator.comparingInt(FoldCandidate::folds))
			.toList();

		FoldCandidate best = null;
		for(final FoldCandidate candidate : candidates){
			System.out.printf("  folds=%d -> fitness=%.4f%n", candidate.folds(), candidate.fitness());

			if(best == null || candidate.fitness() < best.fitness())
				best = candidate;
		}
		if(best == null)
			throw new IllegalStateException("No fold count produced a valid schedule");

		// The folds scheduled after the end of the timeline are never given: only the applied ones are reported
		double totalMinutes = 0.;
		for(final StageInput stage : baseInputs.stages())
			totalMinutes += stage.duration() * 60.;
		final double[] schedule = buildSchedule(best, totalMinutes);

		final SimulationInputs finalInputs = withFolds(baseInputs, schedule);
		final double bestYeast = YeastOptimizer.findOptimalYeast(finalInputs, targetProduct, null);
		return new OptimizedFoldResult(schedule, bestYeast, best.fitness());
	}

	private static double[] buildSchedule(final FoldCandidate candidate, final double totalMinutes){
		if(candidate.folds() == 0)
			return new double[0];

		final double[] scheduled = FullRecipeMultivariateObjective.buildFoldSchedule(candidate.folds(),
			candidate.point()[0], candidate.point()[1]);
		int applied = 0;
		for(final double fold : scheduled)
			if(fold < totalMinutes)
				applied ++;
		return Arrays.copyOf(scheduled, applied);
	}

	private static FoldCandidate optimizeForFoldCount(final SimulationInputs baseInputs,
			final BakeryProduct targetProduct, final int folds){
		final MultivariateFunction objective = normalizedPoint -> {
			final double[] point = denormalize(normalizedPoint);
			return evaluate(withFolds(baseInputs, FullRecipeMultivariateObjective.buildFoldSchedule(folds, point[0],
				point[1])), targetProduct);
		};

		// No fold: no variable to search
		if(folds == 0)
			return new FoldCandidate(0, evaluate(withFolds(baseInputs, new double[0]), targetProduct),
				INITIAL_GUESS.clone());

		final double[] normalizedLowerBounds = new double[NUM_CONTINUOUS_VARIABLES];
		final double[] normalizedUpperBounds = new double[NUM_CONTINUOUS_VARIABLES];
		Arrays.fill(normalizedUpperBounds, 1.);
		final double[] normalizedGuess = new double[NUM_CONTINUOUS_VARIABLES];
		for(int i = 0; i < NUM_CONTINUOUS_VARIABLES; i ++)
			normalizedGuess[i] = (INITIAL_GUESS[i] - LOWER_BOUNDS[i]) / (UPPER_BOUNDS[i] - LOWER_BOUNDS[i]);

		final BOBYQAOptimizer optimizer = new BOBYQAOptimizer(2 * NUM_CONTINUOUS_VARIABLES + 1,
			INITIAL_TRUST_REGION_RADIUS, STOPPING_TRUST_REGION_RADIUS);
		try{
			final PointValuePair optimum = optimizer.optimize(new MaxEval(MAX_EVALUATIONS),
				new ObjectiveFunction(objective), GoalType.MINIMIZE, new InitialGuess(normalizedGuess),
				new SimpleBounds(normalizedLowerBounds, normalizedUpperBounds));
			return new FoldCandidate(folds, optimum.getValue(), denormalize(optimum.getPoint()));
		}
		catch(final MathIllegalStateException e){
			System.out.printf("  folds=%d -> skipped (%s)%n", folds, e.getMessage());
			return null;
		}
	}

	private static double evaluate(final SimulationInputs inputs, final BakeryProduct targetProduct){
		final double yeast = YeastOptimizer.findOptimalYeast(inputs, targetProduct, null);
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(inputs);
		final double[] finalState = YeastOptimizer.runSimulation(inputs, gab, yeast, targetProduct);
		if(finalState == null)
			return Double.MAX_VALUE;

		final double baseMaxGasPotential = YeastOptimizer.calculateMaxGasPotential(
			YeastOptimizer.inputsSummary(inputs), inputs, targetProduct);
		return YeastOptimizer.evaluateFitness(inputs, gab, baseMaxGasPotential, finalState, targetProduct);
	}

	private static SimulationInputs withFolds(final SimulationInputs inputs, final double[] folds){
		return new SimulationInputs(inputs.fractions(), inputs.flourMatrix(), inputs.flourTemperature(),
			inputs.airRelativeHumidity(), inputs.yeastProperties(), inputs.recipe(), inputs.kneading(),
			inputs.stages(), folds);
	}

	private static double[] denormalize(final double[] normalizedPoint){
		final double[] point = new double[normalizedPoint.length];
		for(int i = 0; i < point.length; i ++)
			point[i] = LOWER_BOUNDS[i] + normalizedPoint[i] * (UPPER_BOUNDS[i] - LOWER_BOUNDS[i]);
		return point;
	}

}
