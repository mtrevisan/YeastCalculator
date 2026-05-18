package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;

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


public final class Main2{

	// Parameter Object to group and safely transport bio-mechanical dough properties
	private static final class DoughContext{
		private final StageInput[] stages;
		private final double[] folds;
		private final double totalDuration;
		private final double stiffness;
		private final double saltK;
		private final double oilK;
		private final double waterContent;
		private final double sugarInitial;

		public DoughContext(final StageInput[] stages, final double[] folds, final double totalDuration,
				final double stiffness, final double saltK, final double oilK, final double waterContent,
				final double sugarInitial){
			this.stages = stages;
			this.folds = folds;
			this.totalDuration = totalDuration;
			this.stiffness = stiffness;
			this.saltK = saltK;
			this.oilK = oilK;
			this.waterContent = waterContent;
			this.sugarInitial = sugarInitial;
		}

	}

	public static final class SimulationResult{
		private final double optimalYeast;
		private final double peakVolume;
		private final double remainingSugar;
		private final double structuralTearingPenalty;

		public SimulationResult(final double optimalYeast, final double peakVolume, final double remainingSugar,
				final double structuralTearingPenalty){
			this.optimalYeast = optimalYeast;
			this.peakVolume = peakVolume;
			this.remainingSugar = remainingSugar;
			this.structuralTearingPenalty = structuralTearingPenalty;
		}

		public double getOptimalYeast(){
			return optimalYeast;
		}

		public double getPeakVolume(){
			return peakVolume;
		}

		public double getRemainingSugar(){
			return remainingSugar;
		}

		public double getStructuralTearingPenalty(){
			return structuralTearingPenalty;
		}

	}

	public static void main(final String[] args){
		final SimulationInputs inputs = new SimulationInputs();
		final SimulationResult result = calculateOptimalYeast(inputs);

		System.out.printf("Optimal Commercial Yeast Quantity: %.6f%n", result.getOptimalYeast());
		System.out.printf("Peak Volume Ratio (V_max / V_0):   %.4f%n", result.getPeakVolume());
		System.out.printf("Remaining Sugar Substrate:          %.4f%n", result.getRemainingSugar());
		System.out.printf("Structural Tearing Penalty Applied: %.4f%n", result.getStructuralTearingPenalty());
	}

	public static SimulationResult calculateOptimalYeast(final SimulationInputs in){
		// 1. Environmental & Material Setup
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in);
		final double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;
		final double totalWaterContent = (in.getDoughWater() + flourMoisture) / (1. + in.getDoughWater());

		// 2. Extract Timeline Properties
		final StageInput[] stages = Arrays.stream(in.getStages())
			.filter(s -> s != null && s.getDuration() > 0.)
			.toArray(StageInput[]::new);
		final double totalDuration = Arrays.stream(stages)
			.mapToDouble(StageInput::getDuration)
			.sum();
		final double[] folds = (in.getFolds() == null? new double[0]: in.getFolds());

		// 3. Assemble the Parameter Object Context
		final DoughContext context = new DoughContext(
			stages, folds, totalDuration,
			calculateStiffnessIndex(in, flourMoisture),
			Math.exp(-15. * in.getDoughSalt()),
			(1. + 5. * in.getDoughOil()) * Math.exp(-8. * in.getDoughOil()),
			totalWaterContent,
			aggregateFlourSugar(in) + 0.01
		);

		// 4. Define Objective Function Strategy using the wrapped context
		final UnivariateFunction structuralObjectiveFunction = yDry -> {
			final double[] metrics = runDoughSimulation(yDry, context);
			final double peakVolume = metrics[0];
			final double remainingSugar = metrics[1];

			final double structuralTearingPenalty = (peakVolume > 2.2? (peakVolume - 2.2) * 15.: 0.);
			final double yeastAbusePenalty = yDry * 8.;

			return peakVolume + remainingSugar * 5. - structuralTearingPenalty - yeastAbusePenalty;
		};

		// 5. Execute Optimization Search Loop
		final double yDryOptimal = executeBrentOptimization(structuralObjectiveFunction);

		// 6. Diagnostic Run using the Discovered Optimal Yeast Coordinate
		final double[] finalMetrics = runDoughSimulation(yDryOptimal, context);
		final double peakVolume = finalMetrics[0];
		final double remainingSugar = finalMetrics[1];
		final double structuralTearingPenalty = (peakVolume > 2.2? (peakVolume - 2.2) * 15.: 0.);

		final double commercialYeast = yDryOptimal / (1. - in.getYeastMoisture());
		return new SimulationResult(commercialYeast, peakVolume, remainingSugar, structuralTearingPenalty);
	}


	private static double calculateStiffnessIndex(final SimulationInputs in, final double flourMoisture){
		double dotStrength = 0.;
		double dotPL = 0.;
		double dotProtein = 0.;
		double dotFiber = 0.;
		double dotAsh = 0.;
		double dotγ = 0.;
		final double[] fractions = in.getFractions();
		final FlourInput[] matrix = in.getFlourMatrix();
		for(int i = 0; i < in.getFlourCount(); i ++){
			final FlourInput flour = matrix[i];
			final double γ_i = flour.getBaseLookup() * Math.exp(-2. * flour.getFat()) * Math.exp(-0.5 * flour.getSugar());

			dotStrength += flour.getStrengthW() * fractions[i];
			dotPL += flour.getPlRatio() * fractions[i];
			dotProtein += flour.getProtein() * fractions[i];
			dotFiber += flour.getFiber() * fractions[i];
			dotAsh += flour.getAsh() * fractions[i];
			dotγ += γ_i * fractions[i];
		}

		final double waterEff = (in.getDoughWater() + flourMoisture) / (1. - flourMoisture);
		final double structuralHydrationModifier = 1. + Math.pow(waterEff - 0.6, 2);
		return (dotStrength * dotProtein / (dotPL * (1. + 2. * dotFiber + 5. * dotAsh)) * dotγ) / structuralHydrationModifier;
	}

	private static double aggregateFlourSugar(final SimulationInputs in){
		double dotSugar = 0.;
		final double[] fractions = in.getFractions();
		final FlourInput[] matrix = in.getFlourMatrix();
		for(int i = 0; i < in.getFlourCount(); i ++)
			dotSugar += matrix[i].getSugar() * fractions[i];
		return dotSugar;
	}

	// Signature is now incredibly clean: takes just the dry yeast input and the context object container
	private static double[] runDoughSimulation(double yDry, DoughContext ctx){
		final FirstOrderIntegrator integrator = new DormandPrince853Integrator(1e-6, ctx.totalDuration,
			1e-6, 1e-6);
		if(ctx.folds.length > 0)
			integrator.addEventHandler(new FoldEventHandler(ctx.folds), 0.01, 1e-4, 100);

		final double[] metrics = {1., ctx.sugarInitial};
		integrator.addStepHandler(new StepHandler(){
			@Override
			public void init(double t0, double[] y0, double t){
				metrics[0] = y0[0];
				metrics[1] = y0[2];
			}

			@Override
			public void handleStep(StepInterpolator interpolator, boolean isLast){
				double[] yInterp = interpolator.getInterpolatedState();
				if(yInterp[0] > metrics[0]){
					metrics[0] = yInterp[0];
				}
				metrics[1] = yInterp[2];
			}
		});

		final DoughOdeSystem ode = new DoughOdeSystem(yDry, ctx.stages, ctx.stiffness, ctx.saltK, ctx.oilK,
			ctx.waterContent);
		final double[] y = {1., 0., ctx.sugarInitial, 0.};
		try{
			integrator.integrate(ode, 0., y, ctx.totalDuration, y);
		}
		catch(Exception e){
			return new double[]{-100.0, 0.0};
		}
		return metrics;
	}

	private static double executeBrentOptimization(UnivariateFunction objectiveFunction){
		final double lowerBound = 0.0001;
		// Capped at 5% commercial fresh yeast equivalent
		final double upperBound = 0.015;
		final UnivariateOptimizer optimizer = new BrentOptimizer(1e-6, 1e-7);

		final UnivariatePointValuePair optimizationResult = optimizer.optimize(
			new MaxEval(150),
			new UnivariateObjectiveFunction(objectiveFunction),
			GoalType.MAXIMIZE,
			new SearchInterval(lowerBound, upperBound)
		);
		return optimizationResult.getPoint();
	}

}
