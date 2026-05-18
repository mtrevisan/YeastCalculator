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

	// Data carrier to return both optimal yeast and the highest volume ratio achieved
	public static final class SimulationResult{
		private final double optimalYeast;
		private final double maxVolumeRatio;

		public SimulationResult(double optimalYeast, double maxVolumeRatio){
			this.optimalYeast = optimalYeast;
			this.maxVolumeRatio = maxVolumeRatio;
		}

		public double getOptimalYeast(){
			return optimalYeast;
		}

		public double getMaxVolumeRatio(){
			return maxVolumeRatio;
		}

	}

	public static void main(final String[] args){
		final SimulationInputs inputs = new SimulationInputs();

		// Execute calculation and fetch the combined data result
		final SimulationResult result = calculateOptimalYeast(inputs);

		System.out.printf("Optimal Commercial Yeast Quantity for MAX Volume: %.6f%n", result.getOptimalYeast());
		System.out.printf("Max Volume / Initial Volume Ratio: %.2f%n", result.getMaxVolumeRatio());
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
		// RESTORED: Calculating the effective hydration variable representing unbound fluid mechanics
		final double waterEff = (in.getDoughWater() + flourMoisture) / (1. - flourMoisture);

		// Adjusting the structural stiffness base index using the restored water efficiency factor
		// Higher water deviations from the ideal 0.6 reference point soften the matrix structure
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
			.mapToDouble(s -> s.getDuration())
			.sum();
		final double[] folds = (in.getFolds() == null) ? new double[0] : in.getFolds();

		// 5. Objective Function targeting peak volume reached at any point during the timeline
		final UnivariateFunction volumeObjectiveFunction = new UnivariateFunction(){
			@Override
			public double value(double yDry){
				final FirstOrderIntegrator integrator = new DormandPrince853Integrator(1e-6, totalDuration, 1e-6, 1e-6);

				if(folds.length > 0){
					integrator.addEventHandler(new FoldEventHandler(folds), 0.01, 1e-4, 100);
				}

				// Container to track the maximum volume value encountered during simulation steps
				final double[] peakVolumeContainer = {1.0};

				integrator.addStepHandler(new StepHandler(){
					@Override
					public void init(double t0, double[] y0, double t){
						peakVolumeContainer[0] = y0[0];
					}

					@Override
					public void handleStep(StepInterpolator interpolator, boolean isLast){
						double[] yInterp = interpolator.getInterpolatedState();
						if(yInterp[0] > peakVolumeContainer[0]){
							peakVolumeContainer[0] = yInterp[0];
						}
					}
				});

				final DoughOdeSystem ode = new DoughOdeSystem(yDry, stages, stiffnessIndexBase, saltK, oilK, totalWaterContent);

				// Initial boundary conditions: Volume (y[0]) starts at 1.0
				final double[] y = {1., 0., sugarInitial, 0.};

				try{
					integrator.integrate(ode, 0., y, totalDuration, y);
				}
				catch(Exception e){
					return 1.0;
				}

				return peakVolumeContainer[0];
			}
		};

		// 6. Brent Optimizer Execution for MAXIMIZATION
		// NOTE: You can adjust the upperBound here (e.g., 0.015) to restrict search to realistic baking ranges
		final double lowerBound = 0.0001;
		final double upperBound = 0.015;
		final double relativeAccuracy = 1e-6;
		final double absoluteAccuracy = 1e-7;

		final UnivariateOptimizer optimizer = new BrentOptimizer(relativeAccuracy, absoluteAccuracy);

		final UnivariatePointValuePair optimizationResult = optimizer.optimize(
			new MaxEval(100),
			new UnivariateObjectiveFunction(volumeObjectiveFunction),
			GoalType.MAXIMIZE,
			new SearchInterval(lowerBound, upperBound)
		);

		final double yDryOptimal = optimizationResult.getPoint();
		final double maxVolumeReached = optimizationResult.getValue();
		final double commercialYeast = yDryOptimal / (1. - in.getYeastMoisture());

		return new SimulationResult(commercialYeast, maxVolumeReached);
	}

}
