package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.PointValuePair;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.stream.IntStream;


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
	private static final int MAX_FOLD_SETS = 8;

	// Margins by which the search ranges of the recipe exceed the envelope of the product
	private static final double ENVELOPE_MARGIN_HYDRATION = 0.03;
	private static final double ENVELOPE_MARGIN_SALT = 0.003;
	private static final double ENVELOPE_MARGIN_OIL = 0.01;
	private static final double ENVELOPE_MARGIN_MALT = 0.004;


	private record FoldCandidate(int folds, double fitness, double[] point){}

	public record OptimizedRecipeResult(DoughRecipe recipe, double optimalYeastRatio, double finalFitness,
		int optimizedFoldCount, double optimizedFoldInitialDelay, double optimizedRelaxationFactor){}


	public static OptimizedRecipeResult optimizeFullRecipe(final SimulationInputs baseInputs,
			final BakeryProduct targetProduct){
		// --- AUTOMATIC HYDRATION CAP HANDLER (BIOCHEMICAL WAC MODEL) ---
		// Resolves the weighted average characteristics of the flour blend matrix
		double blendedProtein = 0.;
		double blendedFiber = 0.;
		double estimatedDamagedStarch = 0.;
		final double[] fractions = baseInputs.fractions();
		final FlourInput[] flours = baseInputs.flourMatrix();

		for(int i = 0; i < flours.length; i ++){
			blendedProtein += fractions[i] * flours[i].protein();
			blendedFiber += fractions[i] * flours[i].fiber();
			estimatedDamagedStarch += fractions[i] * flours[i].estimatedDamagedStarch();
		}

		// The damaged starch is estimated from the flour strength W (see FlourInput.estimatedDamagedStarch)
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
		// Continuous variable bounds: the search stays around what the product is made with (a hydration of 50% is a
		// biscuit, not a bread), within the water the flour can absorb. The cost of formulation keeps the optimum inside
		// the envelope; the margin lets the search see its edges.
		final BakeryProduct.RecipeEnvelope envelope = targetProduct.getRecipeEnvelope();
		final double hydrationLowerBound = Math.min(Math.max(0.50, envelope.hydrationMin() - ENVELOPE_MARGIN_HYDRATION),
			dynamicHydrationCap - 0.02);
		final double hydrationUpperBound = Math.max(hydrationLowerBound + 0.02,
			Math.min(dynamicHydrationCap, envelope.hydrationMax() + ENVELOPE_MARGIN_HYDRATION));
		final double[] lowerBounds = new double[]{hydrationLowerBound, Math.max(0.01, envelope.saltMin() - ENVELOPE_MARGIN_SALT),
			Math.max(0., envelope.oilMin() - ENVELOPE_MARGIN_OIL), 0., 15., 1.};
		final double[] upperBounds = new double[]{hydrationUpperBound, envelope.saltMax() + ENVELOPE_MARGIN_SALT,
			envelope.oilMax() + ENVELOPE_MARGIN_OIL, envelope.maltMax() + ENVELOPE_MARGIN_MALT, 60., 1.5};

		// Initial guess: the middle of the search ranges for the recipe, a generic schedule for the folds
		final double[] initialGuess = new double[]{(lowerBounds[0] + upperBounds[0]) / 2.,
			(lowerBounds[1] + upperBounds[1]) / 2., (lowerBounds[2] + upperBounds[2]) / 2.,
			(lowerBounds[3] + upperBounds[3]) / 2., 20., 1.2};

		// Execute multivariate optimization execution pass.
		// The search runs in the normalized [0, 1] space.
		final double[] normalizedLowerBounds = new double[NUM_CONTINUOUS_VARIABLES];
		final double[] normalizedUpperBounds = new double[NUM_CONTINUOUS_VARIABLES];
		Arrays.fill(normalizedUpperBounds, 1.);

		// The discrete fold count cannot be handled by BOBYQA. The best count emerges from the gluten development
		// required by the product and the one reached by kneading; the fitness is not assumed to be unimodal, so all
		// counts are scanned. The scans are independent (every one has its own optimizer and the simulation is
		// stateless), so they run in parallel.
		System.out.println("Scanning fold counts 0 to " + MAX_FOLD_SETS + "...");
		final List<FoldCandidate> candidates = IntStream.rangeClosed(0, MAX_FOLD_SETS)
			.parallel()
			.mapToObj(folds -> {
				final FullRecipeMultivariateObjective objective = new FullRecipeMultivariateObjective(baseInputs,
					targetProduct, lowerBounds, upperBounds, folds);

				final BOBYQAOptimizer optimizer = new BOBYQAOptimizer(2 * NUM_CONTINUOUS_VARIABLES + 1,
					INITIAL_TRUST_REGION_RADIUS, STOPPING_TRUST_REGION_RADIUS);
				try{
					final PointValuePair optimum = optimizer.optimize(new MaxEval(MAX_EVALUATIONS),
						new ObjectiveFunction(objective), GoalType.MINIMIZE,
						new InitialGuess(objective.normalize(initialGuess)),
						new SimpleBounds(normalizedLowerBounds, normalizedUpperBounds));
					return new FoldCandidate(folds, optimum.getValue(), objective.denormalize(optimum.getPoint()));
				}
				catch(final MathIllegalStateException e){
					// A flat or discontinuous landscape can stop the trust region: this fold count is skipped
					System.out.printf("  folds=%d -> skipped (%s)%n", folds, e.getMessage());
					return null;
				}
			})
			.filter(Objects::nonNull)
			.sorted(Comparator.comparingInt(FoldCandidate::folds))
			.toList();

		FoldCandidate best = null;
		for(final FoldCandidate candidate : candidates){
			final double[] p = candidate.point();
			System.out.printf("  folds=%d -> fitness=%.4f (hydration=%.1f%%, salt=%.2f%%, oil=%.2f%%, malt=%.2f%%)\n",
				candidate.folds(), candidate.fitness(), p[0] * 100., p[1] * 100., p[2] * 100., p[3] * 100.);

			if(best == null || candidate.fitness() < best.fitness())
				best = candidate;
		}
		if(best == null)
			throw new IllegalStateException("No fold count produced a valid recipe");

		final double[] p = best.point();
		final DoughRecipe bestRecipe = new DoughRecipe(p[0], p[1], p[2], p[3],
			baseInputs.recipe().maltSugarContent(), baseInputs.recipe().maltPollakUnit(),
			baseInputs.recipe().mixerFrictionFactor());

		// The folds scheduled after the end of the timeline are never given: only the applied ones are reported
		double totalMinutes = 0.;
		for(final StageInput stage : baseInputs.stages())
			totalMinutes += stage.duration() * 60.;
		final double[] scheduledFolds = FullRecipeMultivariateObjective.buildFoldSchedule(best.folds(), p[4], p[5]);
		int appliedFolds = 0;
		for(final double fold : scheduledFolds)
			if(fold < totalMinutes)
				appliedFolds ++;
		final double[] bestFolds = Arrays.copyOf(scheduledFolds, appliedFolds);
		final SimulationInputs finalInputs = new SimulationInputs(baseInputs.fractions(),
			baseInputs.flourMatrix(), baseInputs.flourTemperature(), baseInputs.airRelativeHumidity(),
			baseInputs.yeastProperties(), bestRecipe, baseInputs.kneading(), baseInputs.stages(), bestFolds);

		final double bestYeast = YeastOptimizer.findOptimalYeast(finalInputs, targetProduct, null);

		return new OptimizedRecipeResult(bestRecipe, bestYeast, best.fitness(), appliedFolds, p[4], p[5]);
	}

}
