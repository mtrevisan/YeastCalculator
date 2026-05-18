package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;


public class Main2{

	// wet-basis moisture breakdown
	public static class GabResult{
		public double flourStrictlyBoundWater;
		public double flourActiveWater;

	}

	public static void main(final String[] args){
		final SimulationInputs inputs = new SimulationInputs();
		final double dt = 2. / 60.;
		final double result = calculateOptimalYeast(inputs, dt);

		System.out.printf("Optimal Commercial Yeast Quantity: %.6f%n", result);
	}

	public static double calculateOptimalYeast(final SimulationInputs in, final double dt){
		// 1. Flour Normalization
		final int flours = in.getFlourCount();
		final double[] fractions = in.getFractions();

		// 2. GAB Moisture Modeling
		final GabResult gab = calculateGabMoisture(in, fractions);
		final double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;

		// 3. Weighted Flour Property Aggregation (Dot Product)
		double dotStrength = 0, dotPL = 0, dotSugar = 0, dotProtein = 0, dotFiber = 0, dotAsh = 0, dotγ = 0;
		final List<String> flourTypesLookup = Arrays.asList(
			"wheat", "wheat semolina", "wheat semolina fine", "durum wheat", "durum wheat semolina",
			"durum wheat semolina fine", "rye", "einkorn", "emmer", "spelt", "buckwheat", "barley", "chestnut"
		);
		final double[] baseLookup = {1.0, 0.96, 0.98, 0.85, 0.82, 0.85, 0.35, 0.52, 0.68, 0.78, 0.05, 0.25, 0.02};

		for(int i = 0; i < flours; i++){
			final double strength = ((Number)in.flourMatrix[i][0]).doubleValue();
			final double pl = ((Number)in.flourMatrix[i][1]).doubleValue();
			final double sugar = ((Number)in.flourMatrix[i][2]).doubleValue();
			final double protein = ((Number)in.flourMatrix[i][3]).doubleValue();
			final double fats = ((Number)in.flourMatrix[i][4]).doubleValue();
			final double fiber = ((Number)in.flourMatrix[i][5]).doubleValue();
			final double ash = ((Number)in.flourMatrix[i][6]).doubleValue();
			final String type = (String)in.flourMatrix[i][7];

			final int fidx = flourTypesLookup.indexOf(type);
			final double base = (fidx != -1) ? baseLookup[fidx] : 1.0;
			final double γ_i = base * Math.exp(-2.0 * fats) * Math.exp(-0.5 * sugar);

			dotStrength += strength * fractions[i];
			dotPL += pl * fractions[i];
			dotSugar += sugar * fractions[i];
			dotProtein += protein * fractions[i];
			dotFiber += fiber * fractions[i];
			dotAsh += ash * fractions[i];
			dotγ += γ_i * fractions[i];
		}

		// 4. Bio-Mechanical Dough Environment Setup
		final double waterEff = (in.doughWater + flourMoisture) / (1.0 - flourMoisture);
		final double stiffnessIndexBase = dotStrength * dotProtein / (dotPL * (1.0 + 2.0 * dotFiber + 5.0 * dotAsh)) * dotγ;
		final double vTarget = 1.0 + 0.005 * stiffnessIndexBase / (1.0 + Math.pow(waterEff - 0.6, 2));
		final double sugarInitial = dotSugar + 0.01;
		final double saltK = Math.exp(-15.0 * in.doughSalt);
		final double oilK = (1.0 + 5.0 * in.doughOil) * Math.exp(-8.0 * in.doughOil);

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
		final int maxEvaluations = 100; // Increased safety margin ceiling

		// Wrap simulation into Apache's UnivariateFunction interface
		final UnivariateFunction targetFunction = new UnivariateFunction(){
			@Override
			public double value(double yDry){
				// Solver expects f(x) = 0. We want simulateVolume == vTarget, so we solve: simulateVolume - vTarget = 0
				return simulateVolume(yDry, totalSteps, dt, stagesIdx, stages, folds,
					stiffnessIndexBase, sugarInitial, saltK, oilK) - vTarget;
			}
		};

		// Configure the Bisection solver with the customized threshold
		final BisectionSolver solver = new BisectionSolver(absoluteAccuracy);
		final double yDryOptimal = solver.solve(maxEvaluations, targetFunction, lowerBound, upperBound);

		return yDryOptimal / (1.0 - in.yeastMoisture);
	}

	// --- GUGGENHEIM-ANDERSON-DE BOER (GAB) WATER ADSORPTION MODEL ---
	private static GabResult calculateGabMoisture(final SimulationInputs in, final double[] fractions){
		final List<String> flourTypesLookup = Arrays.asList(
			"wheat", "wheat semolina", "wheat semolina fine", "durum wheat", "durum wheat semolina",
			"durum wheat semolina fine", "rye", "einkorn", "emmer", "spelt", "buckwheat", "barley", "chestnut"
		);

		final double[] wmBaseLookup = {0.062, 0.06, 0.061, 0.063, 0.062, 0.062, 0.071, 0.064, 0.065, 0.063, 0.055, 0.068, 0.05};
		final double[] cGabLookup = {12.5, 11.8, 12.1, 13.2, 12.6, 12.8, 8.4, 13.0, 12.4, 11.9, 6.2, 9.5, 4.8};
		final double[] kGabLookup = {0.78, 0.79, 0.78, 0.76, 0.77, 0.77, 0.82, 0.75, 0.77, 0.79, 0.85, 0.81, 0.88};

		double mixTotalDB = 0.0;
		double mixBoundDB = 0.0;

		// Use Apache Commons Lang to clamp variables safely
		final double aw = StrictMath.max(StrictMath.min(in.airRelativeHumidity, 0.95), 0.1);
		final double temperatureCoeff = 1.0 - 0.0025 * (in.flourTemperature - 20.0);

		for(int i = 0; i < fractions.length; i++){
			final double protein = ((Number)in.flourMatrix[i][3]).doubleValue();
			final double fiber = ((Number)in.flourMatrix[i][6]).doubleValue();
			final String type = (String)in.flourMatrix[i][7];

			int fidx = flourTypesLookup.indexOf(type);
			if(fidx == -1){
				fidx = 0;
			}

			final double wmBase = wmBaseLookup[fidx];
			final double cGab = cGabLookup[fidx];
			final double kGab = kGabLookup[fidx];

			// Monolayer moisture content on a Dry Basis (DB)
			final double wmDb = wmBase + 0.085 * protein + 0.12 * fiber;

			// Thermodynamic equilibrium sorption formula
			final double uEquilibriumDB = (wmDb * cGab * kGab * aw) / ((1.0 - kGab * aw) * (1.0 - kGab * aw + cGab * kGab * aw));

			// Apply thermal drift scale and bound it between [0.08, 0.2] DB
			final double uTotalDB = StrictMath.max(StrictMath.min(uEquilibriumDB * temperatureCoeff, 0.2), 0.08);

			mixTotalDB += uTotalDB * fractions[i];
			mixBoundDB += wmDb * fractions[i];
		}

		// Dry Basis to Wet Basis conversions -> x / (1 + x)
		final double moistureTotal = mixTotalDB / (1.0 + mixTotalDB);
		final double flourStrictlyBoundWater = mixBoundDB / (1.0 + mixBoundDB);
		final double flourActiveWater = moistureTotal - flourStrictlyBoundWater;

		final GabResult res = new GabResult();
		res.flourStrictlyBoundWater = flourStrictlyBoundWater;
		res.flourActiveWater = flourActiveWater;
		return res;
	}

	// --- CONTINUOUS TIME ITERATIVE FERMENTATION SIMULATION ---
	private static double simulateVolume(final double yDry, final int totalSteps, final double dt,
		final int[] stagesIdx, final double[][] stages, final double[] folds,
		final double stiffnessIndexBase, final double sugarInitial, final double saltK, final double oilK){

		double vPrev = 1.0;
		double stiffnessIndexPrev = stiffnessIndexBase;
		double lagPrev = 0.0;
		double sugarPrev = sugarInitial;
		double tPrev = 0.0;

		for(int step = 1; step <= totalSteps; step++){
			final double tCurrTime = step * dt;
			final int stageIdx = stagesIdx[step - 1];

			final double tCurr = stages[stageIdx][0];
			final double rhCurr = stages[stageIdx][1];

			// Precise intersection check over continuous timeline slice window
			boolean isFoldStep = false;
			for(final double fold : folds){
				if(fold > tPrev && fold <= tCurrTime){
					isFoldStep = true;
					break;
				}
			}

			// Rheological structural relaxation models
			final double stiffnessTarget = stiffnessIndexBase / ((rhCurr < 0.60) ? (0.50 + 0.50 * (rhCurr / 0.60)) : 1.0);
			final double proteolysisRate = 0.002 * Math.exp(0.07 * (tCurr - 20.0));
			final double stiffnessRelaxed = stiffnessTarget * Math.exp(-proteolysisRate * tCurrTime);

			final double stiffnessIndexDynamic = isFoldStep
				? stiffnessIndexPrev * 1.35
				: stiffnessRelaxed + (stiffnessIndexPrev - stiffnessRelaxed) * Math.exp(-1.8 * dt);

			// Microbial growth activation kinetics
			final double alphaBio = calculateCtmi(tCurr, 4.0, 34.0, 42.0);
			final double lagNew = lagPrev + dt * alphaBio;

			final double sugarK = sugarPrev / (sugarPrev + 0.005);
			final double muBio = 205.0 * yDry * alphaBio * saltK * oilK * sugarK / (1.0 + Math.exp(-4.0 * (lagPrev - 0.5)));
			final double muEff = muBio * stiffnessIndexBase / stiffnessIndexPrev;

			final double sugarNew = Math.max(0.0, sugarPrev - (muBio * 0.015 * dt));

			// Viscoelastic pressure retention and degassing metrics
			final double pressureFactor = (vPrev - 1.0) * stiffnessIndexDynamic / stiffnessIndexBase;
			final double leakingCoeff = 1.0 / (1.0 + Math.exp(-12.0 * (pressureFactor - 1.8)));
			final double leakingK = 1.0 - (leakingCoeff * 0.4);

			final double vGenerated = vPrev + muEff * vPrev * dt * leakingK;
			final double vNew = isFoldStep ? (1.0 + (vGenerated - 1.0) * 0.75) : vGenerated;

			// Shift states forward
			vPrev = vNew;
			stiffnessIndexPrev = stiffnessIndexDynamic;
			lagPrev = lagNew;
			sugarPrev = sugarNew;
			tPrev = tCurrTime;
		}

		return vPrev;
	}

	// Cardinal Temperature Model with Inspection — Rosso et al. (1995)
	private static double calculateCtmi(final double temperature, final double yeastTemperatureMin,
			final double yeastTemperatureOpt, final double yeastTemperatureMax){
		if(temperature <= yeastTemperatureMin || temperature >= yeastTemperatureMax)
			return 0.;

		final double numerator = (temperature - yeastTemperatureMax) * Math.pow(temperature - yeastTemperatureMin, 2);
		final double denominator = (yeastTemperatureOpt - yeastTemperatureMin)
			* ((yeastTemperatureOpt - yeastTemperatureMin) * (temperature - yeastTemperatureOpt)
			- (yeastTemperatureOpt - yeastTemperatureMax) * (yeastTemperatureOpt + yeastTemperatureMin - 2.0 * temperature));
		return numerator / denominator;
	}

}