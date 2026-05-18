package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;

import io.github.mtrevisan.yeastcalculator.working.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.working.domain.GabMoistureModel;
import io.github.mtrevisan.yeastcalculator.working.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.working.domain.YeastFermentationModel;
import io.github.mtrevisan.yeastcalculator.working.optimization.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.working.simulation.DoughContext;
import io.github.mtrevisan.yeastcalculator.working.simulation.DoughOdeSystem;
import io.github.mtrevisan.yeastcalculator.working.simulation.FoldEventHandler;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariateOptimizer;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;


/**
 * Orchestrator class acting as the primary entry point for the simulation program.
 * <p>
 * This class coordinates the execution by mapping environmental fields, launching
 * the continuous 6D Dormand-Prince differential solver loop, and packaging the metrics
 * into a Brent scalar optimizer to locate the exact mathematical ceiling for yeast dosing.
 * </p>
 */
public final class Main2{

	private static final double YEAST_ABUSE_PENALTY_MULTIPLIER = 8.;
	private static final double AMYLASE_ACCESSIBLE_SUGAR_OFFSET = 0.01;
	private static final double MAX_SAFE_DRY_YEAST_LIMIT = 0.015;
	private static final double MIN_DRY_YEAST_LIMIT = 0.0001;

	/**
	 * Clean, boilerplate-free data transfer container modeling optimized output dimensions.
	 */
	public record SimulationResult(
		double optimalYeast,
		double peakVolume,
		double remainingSugar
	){}

	/**
	 * GC-Friendly static StepHandler wrapper.
	 * Eliminates high-frequency object allocations on the Heap inside optimization loops.
	 */
	private static final class DoughMetricsTracker implements StepHandler{
		private final double[] metricsRef;

		public DoughMetricsTracker(final double[] metricsRef){
			this.metricsRef = metricsRef;
		}

		@Override
		public void init(final double t0, final double[] y0, final double t){
			this.metricsRef[0] = y0[0];
			this.metricsRef[1] = y0[2];
		}

		@Override
		public void handleStep(final StepInterpolator interpolator, final boolean isLast){
			final double[] yInterp = interpolator.getInterpolatedState();
			// Retain peak maximum volume expansion encountered over any integration step
			if(yInterp[0] > this.metricsRef[0])
				this.metricsRef[0] = yInterp[0];
			// Keep tracking active terminal simple carbohydrates remaining
			this.metricsRef[1] = yInterp[2];
		}
	}


	/**
	 * CLI execution anchor method. Loads input models and logs structured baking results.
	 */
	public static void main(final String[] args){
		final SimulationInputs inputs = new SimulationInputs();
		final SimulationResult result = calculateOptimalYeast(inputs);

		System.out.printf("Optimal Yeast [%%]:   %.4f%n", result.optimalYeast());
		System.out.printf("Peak Volume Ratio:   %.1f%n", result.peakVolume());
		System.out.printf("Remaining Sugar [%%]: %.4f%n", result.remainingSugar());
	}

	/**
	 * Configures environmental setups, defines target optimization curves,
	 * and fires the numerical search algorithm.
	 *
	 * @param in	The unified SimulationInputs model instance.
	 * @return	Configured SimulationResult containing optimal coordinates.
	 */
	public static SimulationResult calculateOptimalYeast(final SimulationInputs in){
		// 1. Environmental & Material Setup via GAB Isotherms
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in);
		final double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;
		final double totalWaterContent = (in.getDoughWater() + flourMoisture) / (1. + in.getDoughWater());

		// 2. Timeline Aggregation
		final StageInput[] stages = Arrays.stream(in.getStages())
			.filter(s -> s != null && s.getDuration() > 0.)
			.toArray(StageInput[]::new);
		final double totalDuration = Arrays.stream(stages)
			.mapToDouble(StageInput::getDuration)
			.sum();
		final double[] folds = (in.getFolds() == null? new double[0]: in.getFolds());

		final double sugarInitial = in.getBlendedSugar() + AMYLASE_ACCESSIBLE_SUGAR_OFFSET;
		final BakeryProduct product = in.getTargetProduct();

		// 3. Context packaging using direct safe public record constructor access
		final DoughContext context = DoughContext.create(
			stages, folds, totalDuration,
			in.getHydratedStiffnessIndex(flourMoisture),
			YeastFermentationModel.calculateSaltInhibition(in.getDoughSalt()),
			YeastFermentationModel.calculateOilInhibition(in.getDoughOil()),
			totalWaterContent, sugarInitial
		);

		// 4. Optimization Engine Function Objective Blueprint
		final UnivariateFunction structuralObjectiveFunction = yDry -> {
			final double[] metrics = runDoughSimulation(yDry, context, product);
			final double peakVolume = metrics[0];
			final double remainingSugar = metrics[1];

			// Elastic structural failure threshold penalty
			final double structuralTearingPenalty = (peakVolume > product.getGlutenTearingLimit()
				? (peakVolume - product.getGlutenTearingLimit()) * product.getTearingPenaltyMultiplier()
				: 0.);

			// Chemical starch depletion safety barrier penalty
			final double sugarStarvationPenalty = (remainingSugar < product.getMinSafeSugarThreshold()
				? (product.getMinSafeSugarThreshold() - remainingSugar) * product.getStarvationPenaltyMultiplier()
				: 0.);

			final double yeastAbusePenalty = yDry * YEAST_ABUSE_PENALTY_MULTIPLIER;

			// Returns unified performance score targeted for maximize loop routes
			return peakVolume - structuralTearingPenalty - sugarStarvationPenalty - yeastAbusePenalty;
		};

		// 5. Brent Search Execution Loop
		final double yDryOptimal = executeBrentOptimization(structuralObjectiveFunction);
		final double[] finalMetrics = runDoughSimulation(yDryOptimal, context, product);

		final double peakVolume = finalMetrics[0];
		final double remainingSugar = finalMetrics[1];

		// 6. Rescale dry biomass scalar back to real commercial wet yeast numbers
		final double commercialYeast = yDryOptimal / (1. - in.getYeastMoisture());
		return new SimulationResult(commercialYeast, peakVolume, remainingSugar);
	}

	/**
	 * Configures and executes a single isolated continuous integration run of the 6D ODE network.
	 */
	private static double[] runDoughSimulation(final double yDry, final DoughContext ctx, final BakeryProduct product){
		final FirstOrderIntegrator integrator = new DormandPrince853Integrator(1e-6, ctx.getTotalDuration(),
			1e-6, 1e-6);
		if(ctx.getFolds().length > 0)
			integrator.addEventHandler(new FoldEventHandler(ctx.getFolds()), 0.01, 1e-4, 100);

		final double[] metrics = {1., ctx.getSugarInitial()};
		integrator.addStepHandler(new DoughMetricsTracker(metrics));

		// Boundary state vector initial conditions assignment array
		// y[0] = 1 (Initial Volume)
		// y[1] = 0 (Initial Lag Phase tracking coordinate)
		// y[2] = sugarInitial (Initial carbohydrate substrate pool)
		// y[3] = 0 (Initial dissolved carbon dioxide concentration)
		// y[4] = 0 (Initial macro-structural protein network degradation)
		// y[5] = 0.001 (Atmospheric/equilibrium basal micro-bubble gas pressure reference)
		final double[] y = {1., 0., ctx.getSugarInitial(), 0., 0., 0.001};

		final DoughOdeSystem ode = new DoughOdeSystem(yDry, ctx.getStages(), ctx.getStiffness(), ctx.getSaltK(),
			ctx.getOilK(), ctx.getWaterContent(), ctx.getSugarInitial(), product.getGlutenTearingLimit());
		try{
			integrator.integrate(ode, 0., y, ctx.getTotalDuration(), y);
		}
		catch(final Exception e){
			// Severe error tracking escape coordinate
			return new double[]{-100., 0.};
		}
		return metrics;
	}

	/**
	 * Evaluates and returns the highest efficiency coordinate using univariate Brent searching.
	 */
	private static double executeBrentOptimization(final UnivariateFunction objectiveFunction){
		final UnivariateOptimizer optimizer = new BrentOptimizer(1e-6, 1e-7);
		final UnivariatePointValuePair optimizationResult = optimizer.optimize(
			new MaxEval(150),
			new UnivariateObjectiveFunction(objectiveFunction),
			GoalType.MAXIMIZE,
			new SearchInterval(MIN_DRY_YEAST_LIMIT, MAX_SAFE_DRY_YEAST_LIMIT)
		);
		return optimizationResult.getPoint();
	}

}
