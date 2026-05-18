package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;


public final class Main2{

	public static void main(final String[] args){
		final SimulationInputs inputs = new SimulationInputs();
		final double result = calculateOptimalYeast(inputs);
		System.out.printf("Optimal Commercial Yeast Quantity: %.6f%n", result);
	}

	public static double calculateOptimalYeast(final SimulationInputs in){
		final int flours = in.getFlourCount();
		final double[] fractions = in.getFractions();

		// 1. GAB Moisture Modeling
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in, fractions);
		final double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;

		// 2. Weighted Flour Property Aggregation (Dot Product)
		double dotStrength = 0.;
		double dotPL = 0.;
		double dotSugar = 0.;
		double dotProtein = 0.;
		double dotFiber = 0.;
		double dotAsh = 0.;
		double dotγ = 0.;
		final Object[][] matrix = in.flourMatrix;

		for(int i = 0; i < flours; i++){
			final double strength = ((Number)matrix[i][0]).doubleValue();
			final double pl = ((Number)matrix[i][1]).doubleValue();
			final double sugar = ((Number)matrix[i][2]).doubleValue();
			final double protein = ((Number)matrix[i][3]).doubleValue();
			final double fats = ((Number)matrix[i][4]).doubleValue();
			final double fiber = ((Number)matrix[i][5]).doubleValue();
			final double ash = ((Number)matrix[i][6]).doubleValue();
			final String type = (String)matrix[i][7];

			final FlourRegistry.FlourProperties props = FlourRegistry.resolveProperties(type);
			final double γ_i = props.baseLookup * Math.exp(-2. * fats) * Math.exp(-0.5 * sugar);

			dotStrength += strength * fractions[i];
			dotPL += pl * fractions[i];
			dotSugar += sugar * fractions[i];
			dotProtein += protein * fractions[i];
			dotFiber += fiber * fractions[i];
			dotAsh += ash * fractions[i];
			dotγ += γ_i * fractions[i];
		}

		// 3. Bio-Mechanical Dough Environment Setup
		final double waterEff = (in.doughWater + flourMoisture) / (1. - flourMoisture);
		final double stiffnessIndexBase = dotStrength * dotProtein / (dotPL * (1. + 2. * dotFiber + 5. * dotAsh)) * dotγ;
		final double vTarget = 1. + 0.005 * stiffnessIndexBase / (1. + Math.pow(waterEff - 0.6, 2));
		final double sugarInitial = dotSugar + 0.01;
		final double saltK = Math.exp(-15. * in.doughSalt);
		final double oilK = (1. + 5. * in.doughOil) * Math.exp(-8. * in.doughOil);

		// 4. Extract Profiles and Total Fermentation Lifespan (T max)
		final double[][] stages = Arrays.stream(in.stagesRaw)
			.filter(r -> r.length >= 3 && r[2] > 0)
			.toArray(double[][]::new);

		final double totalDurationHours = Arrays.stream(stages).mapToDouble(r -> r[2]).sum();
		final double[] folds = (in.foldsRaw == null) ? new double[0] : in.foldsRaw;

		// 5. Target Objective Function for the Bisection Root Finder
		final UnivariateFunction targetFunction = new UnivariateFunction(){
			@Override
			public double value(double yDry){
				// Adaptive step-size integrator setup (Dormand-Prince 8(5,3))
				final FirstOrderIntegrator integrator = new DormandPrince853Integrator(1e-6, totalDurationHours, 1e-6, 1e-6);

				// Add discrete event handler interceptor for folds
				if(folds.length > 0){
					integrator.addEventHandler(new FoldEventHandler(folds), 0.01, 1e-4, 100);
				}

				final DoughOdeSystem ode = new DoughOdeSystem(yDry, stages, folds, stiffnessIndexBase, saltK, oilK);

				// Initial boundary states at t = 0.
				final double[] y = {1., 0., sugarInitial}; // Volume=1., Lag=0., Initial Sugar

				// Solve the entire time profile analytically using dynamic adaptive steps
				integrator.integrate(ode, 0., y, totalDurationHours, y);

				// Return residual delta compared to target volume constraint
				return y[0] - vTarget;
			}
		};

		// 6. Solvers execution
		final double lowerBound = 0.0001;
		final double upperBound = 0.15;
		final double absoluteAccuracy = 1e-7;
		final int maxEvaluations = 100;

		final BisectionSolver solver = new BisectionSolver(absoluteAccuracy);
		final double yDryOptimal = solver.solve(maxEvaluations, targetFunction, lowerBound, upperBound);

		return yDryOptimal / (1. - in.yeastMoisture);
	}

}
