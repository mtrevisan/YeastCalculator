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

	public static final class SimulationResult{
		private final double optimalYeast;
		private final double peakVolume;
		private final double remainingSugar;
		private final double structuralTearingPenalty;

		public SimulationResult(double optimalYeast, double peakVolume, double remainingSugar, double structuralTearingPenalty){
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
		final int flours = in.getFlourCount();
		final double[] fractions = in.getFractions();

		// 1. GAB Moisture Modeling
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in);
		final double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;
		final double totalWaterContent = (in.getDoughWater() + flourMoisture) / (1. + in.getDoughWater());

		// 2. Weighted Flour Property Aggregation (Dot Product)
		double dotStrength = 0.;
		double dotPL = 0.;
		double dotSugar = 0.;
		double dotProtein = 0.;
		double dotFiber = 0.;
		double dotAsh = 0.;
		double dotγ = 0.;
		final FlourInput[] matrix = in.getFlourMatrix();
		for(int i = 0; i < flours; i++){
			final FlourInput flour = matrix[i];
			final double strength = flour.getStrengthW();
			final double pl = flour.getPlRatio();
			final double sugar = flour.getSugar();
			final double protein = flour.getProtein();
			final double fats = flour.getFat();
			final double fiber = flour.getFiber();
			final double ash = flour.getAsh();
			final double baseLookup = flour.getBaseLookup();
			final double γ_i = baseLookup * Math.exp(-2. * fats) * Math.exp(-0.5 * sugar);

			dotStrength += strength * fractions[i];
			dotPL += pl * fractions[i];
			dotSugar += sugar * fractions[i];
			dotProtein += protein * fractions[i];
			dotFiber += fiber * fractions[i];
			dotAsh += ash * fractions[i];
			dotγ += γ_i * fractions[i];
		}

		// 3. Bio-Mechanical Dough Environment Setup
		final double waterEff = (in.getDoughWater() + flourMoisture) / (1. - flourMoisture);
		final double structuralHydrationModifier = 1. + Math.pow(waterEff - 0.6, 2);
		final double stiffnessIndexBase = (dotStrength * dotProtein / (dotPL * (1. + 2. * dotFiber + 5. * dotAsh)) * dotγ) / structuralHydrationModifier;

		final double sugarInitial = dotSugar + 0.01;
		final double saltK = Math.exp(-15. * in.getDoughSalt());
		final double oilK = (1. + 5. * in.getDoughOil()) * Math.exp(-8. * in.getDoughOil());

		// 4. Extract Stages
		final StageInput[] stages = Arrays.stream(in.getStages())
			.filter(s -> s != null && s.getDuration() > 0.)
			.toArray(StageInput[]::new);

		final double totalDuration = Arrays.stream(stages)
			.mapToDouble(StageInput::getDuration)
			.sum();
		final double[] folds = (in.getFolds() == null) ? new double[0] : in.getFolds();

		// 5. Objective Function targeting maximum structural efficiency
		final UnivariateFunction structuralObjectiveFunction = yDry -> {
			final FirstOrderIntegrator integrator = new DormandPrince853Integrator(1e-6, totalDuration, 1e-6, 1e-6);

			if(folds.length > 0){
				integrator.addEventHandler(new FoldEventHandler(folds), 0.01, 1e-4, 100);
			}

			final double[] metrics = {1., sugarInitial};

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

			final DoughOdeSystem ode = new DoughOdeSystem(yDry, stages, stiffnessIndexBase, saltK, oilK, totalWaterContent);
			final double[] y = {1., 0., sugarInitial, 0.};

			try{
				integrator.integrate(ode, 0., y, totalDuration, y);
			}
			catch(Exception e){
				return -100.;
			}

			double peakVolume = metrics[0];
			double remainingSugar = metrics[1];

			double structuralTearingPenalty = 0.;
			if(peakVolume > 2.2){
				structuralTearingPenalty = (peakVolume - 2.2) * 15.;
			}

			double yeastAbusePenalty = yDry * 8.;

			return peakVolume + remainingSugar * 5. - structuralTearingPenalty - yeastAbusePenalty;
		};

		// 6. Brent Optimizer Configuration & Selection Execution
		// APPLIED PHYSICAL LIMIT: Lower bound is 0.01% and Upper Bound is strictly capped
		// at 1.5% dry yeast, which equates to exactly 5% standard commercial fresh yeast.
		final double lowerBound = 0.0001;
		final double upperBound = 0.015;
		final double relativeAccuracy = 1e-6;
		final double absoluteAccuracy = 1e-7;

		final UnivariateOptimizer optimizer = new BrentOptimizer(relativeAccuracy, absoluteAccuracy);

		final UnivariatePointValuePair optimizationResult = optimizer.optimize(
			new MaxEval(150),
			new UnivariateObjectiveFunction(structuralObjectiveFunction),
			GoalType.MAXIMIZE,
			new SearchInterval(lowerBound, upperBound)
		);

		final double yDryOptimal = optimizationResult.getPoint();

		// 7. Diagnostic Pass
		final FirstOrderIntegrator diagnosticIntegrator = new DormandPrince853Integrator(1e-6, totalDuration, 1e-6, 1e-6);
		if(folds.length > 0)
			diagnosticIntegrator.addEventHandler(new FoldEventHandler(folds), 0.01, 1e-4, 100);

		final double[] finalMetrics = {1., sugarInitial};
		diagnosticIntegrator.addStepHandler(new StepHandler(){
			@Override
			public void init(double t0, double[] y0, double t){
				finalMetrics[0] = y0[0];
				finalMetrics[1] = y0[2];
			}

			@Override
			public void handleStep(StepInterpolator interpolator, boolean isLast){
				double[] yInterp = interpolator.getInterpolatedState();
				if(yInterp[0] > finalMetrics[0]){
					finalMetrics[0] = yInterp[0];
				}
				finalMetrics[1] = yInterp[2];
			}
		});

		final DoughOdeSystem finalOde = new DoughOdeSystem(yDryOptimal, stages, stiffnessIndexBase, saltK, oilK, totalWaterContent);
		final double[] finalY = {1., 0., sugarInitial, 0.};

		try{
			diagnosticIntegrator.integrate(finalOde, 0., finalY, totalDuration, finalY);
		}
		catch(Exception ignored){
		}

		double peakVolume = finalMetrics[0];
		double remainingSugar = finalMetrics[1];
		double structuralTearingPenalty = 0.;
		if(peakVolume > 2.2){
			structuralTearingPenalty = (peakVolume - 2.2) * 15.;
		}

		final double commercialYeast = yDryOptimal / (1. - in.getYeastMoisture());

		return new SimulationResult(commercialYeast, peakVolume, remainingSugar, structuralTearingPenalty);
	}

}
