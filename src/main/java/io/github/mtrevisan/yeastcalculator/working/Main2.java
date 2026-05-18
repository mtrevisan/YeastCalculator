package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;
import java.util.stream.IntStream;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;


public final class Main2{

	public static void main(final String[] args){
		final SimulationInputs inputs = new SimulationInputs();
		final double dt = 2. / 60.;
		final double result = calculateOptimalYeast(inputs, dt);
		System.out.printf("Optimal Commercial Yeast Quantity: %.6f%n", result);
	}

	public static double calculateOptimalYeast(final SimulationInputs in, final double dt){
		final int flours = in.getFlourCount();
		final double[] fractions = in.getFractions();

		// 2. GAB Moisture Modeling
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(in, fractions);
		final double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;

		// 3. Weighted Flour Property Aggregation (Dot Product)
		double dotStrength = 0, dotPL = 0, dotSugar = 0, dotProtein = 0, dotFiber = 0, dotAsh = 0, dotγ = 0;

		for(int i = 0; i < flours; i++){
			final double strength = ((Number)in.flourMatrix[i][0]).doubleValue();
			final double pl = ((Number)in.flourMatrix[i][1]).doubleValue();
			final double sugar = ((Number)in.flourMatrix[i][2]).doubleValue();
			final double protein = ((Number)in.flourMatrix[i][3]).doubleValue();
			final double fats = ((Number)in.flourMatrix[i][4]).doubleValue();
			final double fiber = ((Number)in.flourMatrix[i][5]).doubleValue();
			final double ash = ((Number)in.flourMatrix[i][6]).doubleValue();
			final String type = (String)in.flourMatrix[i][7];

			// Isolated chemical lookup via domain registry
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

		// 4. Bio-Mechanical Dough Environment Setup
		final double waterEff = (in.doughWater + flourMoisture) / (1. - flourMoisture);
		final double stiffnessIndexBase = dotStrength * dotProtein / (dotPL * (1. + 2. * dotFiber + 5. * dotAsh)) * dotγ;
		final double vTarget = 1. + 0.005 * stiffnessIndexBase / (1. + Math.pow(waterEff - 0.6, 2));
		final double sugarInitial = dotSugar + 0.01;
		final double saltK = Math.exp(-15. * in.doughSalt);
		final double oilK = (1. + 5. * in.doughOil) * Math.exp(-8. * in.doughOil);

		// 5. Timeline Discretization & Stage Indexing
		final double[][] stages = Arrays.stream(in.stagesRaw)
			.filter(r -> r.length >= 3 && r[2] > 0)
			.toArray(double[][]::new);

		final int[] stepsPerStage = Arrays.stream(stages)
			.mapToInt(r -> (int)Math.ceil(r[2] / dt))
			.toArray();

		final int totalSteps = IntStream.of(stepsPerStage).sum();

		final int[] stagesIdx = new int[totalSteps];
		int currentStepGlobal = 0;
		for(int s = 0; s < stages.length; s++){
			for(int step = 0; step < stepsPerStage[s]; step++){
				stagesIdx[currentStepGlobal] = s;
				currentStepGlobal++;
			}
		}

		final double[] folds = (in.foldsRaw == null) ? new double[0] : in.foldsRaw;

		// 6. Root Finding via Apache Commons Bisection Solver
		final double lowerBound = 0.0001;
		final double upperBound = 0.15;
		final double absoluteAccuracy = 1e-7;
		final int maxEvaluations = 100;

		final UnivariateFunction targetFunction = new UnivariateFunction(){
			@Override
			public double value(double yDry){
				return simulateVolume(yDry, totalSteps, dt, stagesIdx, stages, folds,
					stiffnessIndexBase, sugarInitial, saltK, oilK) - vTarget;
			}
		};

		final BisectionSolver solver = new BisectionSolver(absoluteAccuracy);
		final double yDryOptimal = solver.solve(maxEvaluations, targetFunction, lowerBound, upperBound);

		return yDryOptimal / (1. - in.yeastMoisture);
	}

	// --- CONTINUOUS TIME ITERATIVE FERMENTATION SIMULATION ---
	private static double simulateVolume(final double yDry, final int totalSteps, final double dt,
		final int[] stagesIdx, final double[][] stages, final double[] folds,
		final double stiffnessIndexBase, final double sugarInitial, final double saltK, final double oilK){

		double vPrev = 1.;
		double stiffnessIndexPrev = stiffnessIndexBase;
		double lagPrev = 0.;
		double sugarPrev = sugarInitial;
		double tPrev = 0.;

		for(int step = 1; step <= totalSteps; step++){
			final double tCurrTime = step * dt;
			final int stageIdx = stagesIdx[step - 1];

			final double tCurr = stages[stageIdx][0];
			final double rhCurr = stages[stageIdx][1];

			// 1. Structural Dough Mechanics (Folds & Rheology)
			boolean isFoldStep = false;
			for(final double fold : folds){
				if(fold > tPrev && fold <= tCurrTime){
					isFoldStep = true;
					break;
				}
			}

			final double stiffnessTarget = stiffnessIndexBase / ((rhCurr < 0.60) ? (0.50 + 0.50 * (rhCurr / 0.60)) : 1.);
			final double proteolysisRate = 0.002 * Math.exp(0.07 * (tCurr - 20.));
			final double stiffnessRelaxed = stiffnessTarget * Math.exp(-proteolysisRate * tCurrTime);

			final double stiffnessIndexDynamic = (isFoldStep
				? stiffnessIndexPrev * 1.35
				: stiffnessRelaxed + (stiffnessIndexPrev - stiffnessRelaxed) * Math.exp(-1.8 * dt));

			// 2. Microorganism Biology (Delegated to YeastFermentationModel)
			final double alphaBio = YeastFermentationModel.calculateThermalEfficiency(tCurr);
			final double lagNew = lagPrev + dt * alphaBio;

			final double muBio = YeastFermentationModel.calculateBiomassGrowthRate(
				yDry, alphaBio, lagPrev, sugarPrev, saltK, oilK
			);
			final double sugarNew = YeastFermentationModel.consumeSugar(sugarPrev, muBio, dt);

			// 3. Bio-Mechanical Coupling (How yeast gas interacts with dough macrostructures)
			final double muEff = muBio * stiffnessIndexBase / stiffnessIndexPrev;

			final double pressureFactor = (vPrev - 1.) * stiffnessIndexDynamic / stiffnessIndexBase;
			final double leakingCoeff = 1. / (1. + Math.exp(-12. * (pressureFactor - 1.8)));
			final double leakingK = 1. - (leakingCoeff * 0.4);

			// 4. Volume State Update
			final double vGenerated = vPrev + muEff * vPrev * dt * leakingK;
			final double vNew = isFoldStep ? (1. + (vGenerated - 1.) * 0.75) : vGenerated;

			// State Migration for Next Iteration
			vPrev = vNew;
			stiffnessIndexPrev = stiffnessIndexDynamic;
			lagPrev = lagNew;
			sugarPrev = sugarNew;
			tPrev = tCurrTime;
		}

		return vPrev;
	}

}
