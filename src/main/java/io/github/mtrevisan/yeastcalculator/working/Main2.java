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
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in);
		final double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;
		final double totalWaterContent = (in.getDoughWater() + flourMoisture) / (1.0 + in.getDoughWater());

		// 2. Weighted Flour Property Aggregation (Dot Product)
		double dotStrength = 0.;
		double dotPL = 0.;
		double dotSugar = 0.;
		double dotProtein = 0.;
		double dotFiber = 0.;
		double dotAsh = 0.;
		double dotγ = 0.;
		final FlourInput[] matrix = in.getFlourMatrix();
		for(int i = 0; i < flours; i ++){
			final FlourInput flour = matrix[i];
			final double strength = flour.getStrengthW();
			final double pl = flour.getPlRatio();
			final double sugar = flour.getSugar();
			final double protein = flour.getProtein();
			final double fats = flour.getFat();
			final double fiber = flour.getFiber();
			final double ash = flour.getAsh();
			final String type = flour.getType();

			final FlourRegistry.FlourProperties props = FlourRegistry.resolveProperties(type);
			final double γ_i = props.baseLookup * Math.exp(-2.0 * fats) * Math.exp(-0.5 * sugar);

			dotStrength += strength * fractions[i];
			dotPL += pl * fractions[i];
			dotSugar += sugar * fractions[i];
			dotProtein += protein * fractions[i];
			dotFiber += fiber * fractions[i];
			dotAsh += ash * fractions[i];
			dotγ += γ_i * fractions[i];
		}

		// 3. Bio-Mechanical Dough Environment Setup
		final double waterEff = (in.getDoughWater() + flourMoisture) / (1.0 - flourMoisture);
		final double stiffnessIndexBase = dotStrength * dotProtein / (dotPL * (1.0 + 2.0 * dotFiber + 5.0 * dotAsh)) * dotγ;
		final double vTarget = 1.0 + 0.005 * stiffnessIndexBase / (1.0 + Math.pow(waterEff - 0.6, 2));
		final double sugarInitial = dotSugar + 0.01;
		final double saltK = Math.exp(-15.0 * in.getDoughSalt());
		final double oilK = (1.0 + 5.0 * in.getDoughOil()) * Math.exp(-8.0 * in.getDoughOil());

		// 4. Extract Stages
		final double[][] stages = Arrays.stream(in.getStages())
			.filter(r -> r.length >= 3 && r[2] > 0)
			.toArray(double[][]::new);

		final double totalDurationHours = Arrays.stream(stages).mapToDouble(r -> r[2]).sum();
		final double[] folds = (in.getFolds() == null) ? new double[0] : in.getFolds();

		// 5. Target Objective Function
		final UnivariateFunction targetFunction = new UnivariateFunction(){
			@Override
			public double value(double yDry){
				final FirstOrderIntegrator integrator = new DormandPrince853Integrator(1e-6, totalDurationHours, 1e-6, 1e-6);

				if(folds.length > 0){
					integrator.addEventHandler(new FoldEventHandler(folds), 0.01, 1e-4, 100);
				}

				final DoughOdeSystem ode = new DoughOdeSystem(yDry, stages, stiffnessIndexBase, saltK, oilK, totalWaterContent);

				// Initial boundary conditions:
				// Volume = 1.0, Lag = 0.0, Sugar = sugarInitial, Dissolved CO2 = 0.0
				final double[] y = {1.0, 0.0, sugarInitial, 0.0};

				integrator.integrate(ode, 0.0, y, totalDurationHours, y);
				return y[0] - vTarget;
			}
		};

		// 6. Solvers Execution
		final double lowerBound = 0.0001;
		final double upperBound = 0.15;
		final double absoluteAccuracy = 1e-7;
		final int maxEvaluations = 100;

		final BisectionSolver solver = new BisectionSolver(absoluteAccuracy);
		final double yDryOptimal = solver.solve(maxEvaluations, targetFunction, lowerBound, upperBound);

		return yDryOptimal / (1.0 - in.getYeastMoisture());
	}

}
