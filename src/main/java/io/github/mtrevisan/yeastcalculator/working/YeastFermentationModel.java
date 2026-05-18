package io.github.mtrevisan.yeastcalculator.working;


final class YeastFermentationModel{

	// Saccharomyces cerevisiae thermal cardinal boundaries (Rosso et al., 1995)
	private static final double TEMP_MIN = 4.;
	private static final double TEMP_OPT = 34.;
	private static final double TEMP_MAX = 42.;

	// Kinetic constants
	private static final double POTENTIAL_MU_MAX = 205.;
	private static final double SUGAR_AFFINITY_K = 0.005;
	// Stoichiometric coefficient: kg of sugar consumed per kg of biomass generated
	private static final double YEAST_SUGAR_YIELD_Y = 0.015;
	// Maintenance coefficient: sugar consumed just to keep cells alive per hour, even without growth
	private static final double MAINTENANCE_COEFF_M = 0.0012;

	// Flour Alpha-Amylase maximum conversion velocity (starch -> maltose conversion)
	private static final double AMYLASE_VMAX_BASE = 0.008;


	private YeastFermentationModel(){}


	/**
	 * Cardinal Temperature Model with Inspection (CTMI) — Rosso et al. (1995).
	 * Evaluates the thermal scaling factor [0., 1.] for yeast metabolic activity.
	 */
	static double calculateThermalEfficiency(final double temperature){
		if(temperature <= TEMP_MIN || temperature >= TEMP_MAX){
			return 0.;
		}

		final double numerator = (temperature - TEMP_MAX) * Math.pow(temperature - TEMP_MIN, 2);
		final double denominator = (TEMP_OPT - TEMP_MIN)
			* ((TEMP_OPT - TEMP_MIN) * (temperature - TEMP_OPT)
			- (TEMP_OPT - TEMP_MAX) * (TEMP_OPT + TEMP_MIN - 2. * temperature));

		return numerator / denominator;
	}

	/**
	 * Calculates the specific biological growth/fermentation rate (muBio) based on
	 * cell concentration, enzyme activation state (lag), substrate availability, and biochemical inhibitors.
	 */
	public static double calculateBiomassGrowthRate(final double yDry, final double alphaBio,
		final double lagPrev, final double sugarPrev,
		final double saltK, final double oilK){
		// Substrate-limiting Monod kinetics parameter
		final double sugarK = sugarPrev / (sugarPrev + SUGAR_AFFINITY_K);

		// Sigmoidal lag phase transition factor representing enzyme adjustment
		final double lagTransition = 1. / (1. + Math.exp(-4. * (lagPrev - 0.5)));

		return POTENTIAL_MU_MAX * yDry * alphaBio * saltK * oilK * sugarK * lagTransition;
	}

	/**
	 * Calculates the net sugar rate (dSugar/dt).
	 * Accounts for concurrent enzymatic starch breakdown (generation) and yeast metabolism (consumption).
	 */
	static double calculateNetSugarRate(final double sugarCurr, final double muBio, final double yDry, final double tCurr){
		if(sugarCurr <= 0.0 && muBio <= 0.0){
			return 0.0;
		}

		// Arrhenius-like activation for flour amylase activity based on temperature
		final double amylaseThermalK = Math.exp(0.06 * (tCurr - 20.0)) * (1.0 - 0.005 * Math.pow(tCurr - 35.0, 2));
		final double sugarGeneration = AMYLASE_VMAX_BASE * Math.max(0.0, amylaseThermalK);

		// Total metabolic uptake (Growth requirements + cellular maintenance)
		final double sugarConsumption = (muBio * YEAST_SUGAR_YIELD_Y) + (yDry * MAINTENANCE_COEFF_M);

		return sugarGeneration - sugarConsumption;
	}

}
