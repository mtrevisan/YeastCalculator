package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;


public class Main2{

	//wet-basis
	public static class GabResult{
		public double flourStrictlyBoundWater;
		public double flourActiveWater;
	}


	public static void main(final String[] args){
		final SimulationInputs inputs = new SimulationInputs();
		final double result = calculateOptimalYeast(inputs, 2. / 60.);
		System.out.printf("Quantità di Lievito Ottimale Commerciale: %.6f%n", result);
	}

	public static double calculateOptimalYeast(final SimulationInputs in, final double dt){
		// 1. Normalizzazione Farine
		final int flours = in.getFlourCount();
		final double[] fractions = in.getFractions();

		// 2. CALCOLO GAB
		GabResult gab = calculateGabMoisture(in, fractions);
		// Ricaviamo flourMoisture combinando l'acqua attiva e quella legata
		double flourMoisture = gab.flourActiveWater + gab.flourStrictlyBoundWater;

		// 3. Calcolo proprietà ponderate delle farine
		double dotStrength = 0, dotPL = 0, dotSugar = 0, dotProtein = 0, dotFiber = 0, dotAsh = 0, dotγ = 0;
		List<String> flourTypesLookup = Arrays.asList(
			"wheat", "wheat semolina", "wheat semolina fine", "durum wheat", "durum wheat semolina",
			"durum wheat semolina fine", "rye", "einkorn", "emmer", "spelt", "buckwheat", "barley", "chestnut"
		);
		double[] baseLookup = {1., 0.96, 0.98, 0.85, 0.82, 0.85, 0.35, 0.52, 0.68, 0.78, 0.05, 0.25, 0.02};

		for(int i = 0; i < flours; i++){
			double strength = ((Number)in.flourMatrix[i][0]).doubleValue();
			double pl = ((Number)in.flourMatrix[i][1]).doubleValue();
			double sugar = ((Number)in.flourMatrix[i][2]).doubleValue();
			double protein = ((Number)in.flourMatrix[i][3]).doubleValue();
			double fats = ((Number)in.flourMatrix[i][4]).doubleValue();
			double fiber = ((Number)in.flourMatrix[i][5]).doubleValue();
			double ash = ((Number)in.flourMatrix[i][6]).doubleValue();
			String type = (String)in.flourMatrix[i][7];

			int fidx = flourTypesLookup.indexOf(type);
			double base = (fidx != -1) ? baseLookup[fidx] : 1.;
			double γ_i = base * Math.exp(-2. * fats) * Math.exp(-0.5 * sugar);

			dotStrength += strength * fractions[i];
			dotPL += pl * fractions[i];
			dotSugar += sugar * fractions[i];
			dotProtein += protein * fractions[i];
			dotFiber += fiber * fractions[i];
			dotAsh += ash * fractions[i];
			dotγ += γ_i * fractions[i];
		}

		// 4. Parametri bio-meccanici (Usa il flourMoisture calcolato da GAB)
		double waterEff = (in.doughWater + flourMoisture) / (1. - flourMoisture);
		double stiffnessIndexBase = dotStrength * dotProtein / (dotPL * (1. + 2. * dotFiber + 5. * dotAsh)) * dotγ;
		double vTarget = 1. + 0.005 * stiffnessIndexBase / (1. + Math.pow(waterEff - 0.6, 2));
		double sugarInitial = dotSugar + 0.01;
		double saltK = Math.exp(-15. * in.doughSalt);
		double oilK = (1. + 5. * in.doughOil) * Math.exp(-8. * in.doughOil);

		// 5. Gestione Fasi discretizzate
		double[][] stages = Arrays.stream(in.stagesRaw)
			.filter(r -> r.length >= 3 && r[2] > 0)
			.toArray(double[][]::new);

		int[] stepsPerStage = Arrays.stream(stages)
			.mapToInt(r -> (int)Math.ceil(r[2] / dt))
			.toArray();

		int totalSteps = IntStream.of(stepsPerStage).sum();

		int[] stagesIdx = new int[totalSteps];
		int currentStepGlobal = 0;
		for(int s = 0; s < stages.length; s++){
			for(int step = 0; step < stepsPerStage[s]; step++){
				stagesIdx[currentStepGlobal] = s;
				currentStepGlobal++;
			}
		}

		double[] folds = (in.foldsRaw == null) ? new double[0] : in.foldsRaw;

		// 6. Bisezione ed esecuzione
		double a = 0.0001;
		double b = 0.15;
		int iterations = 22;

		final double finalVTarget = vTarget;
		double yDryOptimal = runBisection(a, b, iterations, finalVTarget, (yDry) ->
			simulateVolume(yDry, totalSteps, dt, stagesIdx, stages, folds,
				stiffnessIndexBase, sugarInitial, saltK, oilK)
		);

		return yDryOptimal / (1. - in.yeastMoisture);
	}

	// --- NUOVO METODO: IMPLEMENTAZIONE ISOTERMA DI ASSORBIMENTO GAB ---
	private static GabResult calculateGabMoisture(SimulationInputs in, double[] fractions){
		List<String> flourTypesLookup = Arrays.asList(
			"wheat", "wheat semolina", "wheat semolina fine", "durum wheat", "durum wheat semolina",
			"durum wheat semolina fine", "rye", "einkorn", "emmer", "spelt", "buckwheat", "barley", "chestnut"
		);

		// Costanti dei modelli GAB per tipo di farina
		double[] wmBaseLookup = {0.062, 0.06, 0.061, 0.063, 0.062, 0.062, 0.071, 0.064, 0.065, 0.063, 0.055, 0.068, 0.05};
		double[] cGabLookup = {12.5, 11.8, 12.1, 13.2, 12.6, 12.8, 8.4, 13.0, 12.4, 11.9, 6.2, 9.5, 4.8};
		double[] kGabLookup = {0.78, 0.79, 0.78, 0.76, 0.77, 0.77, 0.82, 0.75, 0.77, 0.79, 0.85, 0.81, 0.88};

		double mixTotalDB = 0.;
		double mixBoundDB = 0.;

		double aw = clamp(in.airRelativeHumidity, 0.1, 0.95);
		double temperatureCoeff = 1. - 0.0025 * (in.flourTemperature - 20.);

		for(int i = 0; i < fractions.length; i++){
			double protein = ((Number)in.flourMatrix[i][3]).doubleValue();
			double fiber = ((Number)in.flourMatrix[i][6]).doubleValue(); // Indice 6 per le fibre nella formula LET precedente
			String type = (String)in.flourMatrix[i][7];

			int fidx = flourTypesLookup.indexOf(type);
			if(fidx == -1){
				fidx = 0; // IFERROR defaults to 1 (0 in Java base-0)
			}

			double wmBase = wmBaseLookup[fidx];
			double cGab = cGabLookup[fidx];
			double kGab = kGabLookup[fidx];

			// WmDB: Contenuto idrico del monostrato su sostanza secca
			double wmDb = wmBase + 0.085 * protein + 0.12 * fiber;

			// Equazione GAB standard per u_EquilibriumDB
			double uEquilibriumDB = (wmDb * cGab * kGab * aw) / ((1. - kGab * aw) * (1. - kGab * aw + cGab * kGab * aw));

			// Correzione termica e clamp
			double uTotalDB = clamp(uEquilibriumDB * temperatureCoeff, 0.08, 0.2);

			// Somma ponderata (Dot Product)
			mixTotalDB += uTotalDB * fractions[i];
			mixBoundDB += wmDb * fractions[i];
		}

		// Conversione da Base Secca (Dry Basis) a Base Umida (Wet Basis) -> x / (1 + x)
		double moistureTotal = mixTotalDB / (1. + mixTotalDB);
		double flourStrictlyBoundWater = mixBoundDB / (1. + mixBoundDB);
		double flourActiveWater = moistureTotal - flourStrictlyBoundWater;

		GabResult res = new GabResult();
		res.flourStrictlyBoundWater = flourStrictlyBoundWater;
		res.flourActiveWater = flourActiveWater;
		return res;
	}

	// --- FUNZIONI DI UTILITÀ ---
	private static double clamp(double x, double lo, double hi){
		if(x < lo)
			return lo;
		if(x > hi)
			return hi;
		return x;
	}

	// --- SIMULAZIONE CON CONFRONTO DI INTERVALLO CONTINUO ---
	private static double simulateVolume(double yDry, int totalSteps, double dt,
		int[] stagesIdx, double[][] stages, double[] folds,
		double stiffnessIndexBase, double sugarInitial, double saltK, double oilK){
		double vPrev = 1.;
		double stiffnessIndexPrev = stiffnessIndexBase;
		double lagPrev = 0.;
		double sugarPrev = sugarInitial;

		// Teniamo traccia del tempo precedente e attuale (continuo, double)
		double tPrev = 0.;

		for(int step = 1; step <= totalSteps; step++){
			double tCurrTime = step * dt;
			int stageIdx = stagesIdx[step - 1];

			double tCurr = stages[stageIdx][0];
			double rhCurr = stages[stageIdx][1];

			// CONFRONTO ESATTO: Il fold è accaduto nell'intervallo di questo step?
			// Condizione: tPrev < fold <= tCurrTime
			boolean isFoldStep = false;
			for(double fold : folds){
				if(fold > tPrev && fold <= tCurrTime){
					isFoldStep = true;
					break; // Trovato, non serve controllare gli altri per questo step
				}
			}

			// Calcoli di rilassamento stiffness
			double stiffnessTarget = stiffnessIndexBase / ((rhCurr < 0.60) ? (0.50 + 0.50 * (rhCurr / 0.60)) : 1.);
			double proteolysisRate = 0.002 * Math.exp(0.07 * (tCurr - 20.));
			double stiffnessRelaxed = stiffnessTarget * Math.exp(-proteolysisRate * tCurrTime);

			double stiffnessIndexDynamic = (isFoldStep
				? stiffnessIndexPrev * 1.35
				: stiffnessRelaxed + (stiffnessIndexPrev - stiffnessRelaxed) * Math.exp(-1.8 * dt));

			// Crescita biologica e consumi
			double alphaBio = calculateCtmi(tCurr, 4., 34., 42.);
			double lagNew = lagPrev + dt * alphaBio;

			double sugarK = sugarPrev / (sugarPrev + 0.005);
			double muBio = 205. * yDry * alphaBio * saltK * oilK * sugarK / (1. + Math.exp(-4. * (lagPrev - 0.5)));
			double muEff = muBio * stiffnessIndexBase / stiffnessIndexPrev;

			double sugarNew = Math.max(0., sugarPrev - (muBio * 0.015 * dt));

			// Perdita di gas e calcolo volume
			double pressureFactor = (vPrev - 1.) * stiffnessIndexDynamic / stiffnessIndexBase;
			double leakingCoeff = 1. / (1. + Math.exp(-12. * (pressureFactor - 1.8)));
			double leakingK = 1. - (leakingCoeff * 0.4);

			double vGenerated = vPrev + muEff * vPrev * dt * leakingK;
			double vNew = isFoldStep ? (1. + (vGenerated - 1.) * 0.75) : vGenerated;

			// Aggiornamento dello stato per il prossimo ciclo
			vPrev = vNew;
			stiffnessIndexPrev = stiffnessIndexDynamic;
			lagPrev = lagNew;
			sugarPrev = sugarNew;

			// Il tempo attuale diventa il tempo precedente del prossimo step
			tPrev = tCurrTime;
		}

		return vPrev;
	}

	private static double runBisection(double a, double b, int iterations, double comp, java.util.function.Function<Double, Double> f){
		double lo = a;
		double hi = b;
		for(int i = 0; i < iterations; i++){
			double mid = (lo + hi) / 2.;
			if(f.apply(mid) > comp){
				hi = mid;
			}
			else{
				lo = mid;
			}
		}
		return (lo + hi) / 2.;
	}

	// Cardinal temperature model — Rosso et al. (1995)
	// Returns γ_T ∈ [0, 1]. γ_T = 1 at T_OPT.
	private static double calculateCtmi(final double temperature, final double yeastTemperatureMin,
			final double yeastTemperatureOpt, final double yeastTemperatureMax){
		if(temperature <= yeastTemperatureMin || temperature >= yeastTemperatureMax)
			return 0.;

		final double numerator = (temperature - yeastTemperatureMax) * Math.pow(temperature - yeastTemperatureMin, 2);
		final double denominator = (yeastTemperatureOpt - yeastTemperatureMin)
			* ((yeastTemperatureOpt - yeastTemperatureMin) * (temperature - yeastTemperatureOpt)
			- (yeastTemperatureOpt - yeastTemperatureMax) * (yeastTemperatureOpt + yeastTemperatureMin - 2. * temperature));
		return numerator / denominator;
	}

}
