package io.github.mtrevisan.yeastcalculator.working;

public final class YeastFermentationModel{

	// Microorganism biological constants (Saccharomyces cerevisiae)
	private static final double TEMP_MIN = 4.;
	private static final double TEMP_OPT = 34.;
	private static final double TEMP_MAX = 42.;

	private static final double POTENTIAL_MU_MAX = 205.;
	private static final double SUGAR_AFFINITY_K = 0.005; // Monod-like half-velocity constant
	private static final double SUGAR_CONSUMPTION_COEFF = 0.015;

	/**
	 * Cardinal Temperature Model with Inspection (CTMI) — Rosso et al. (1995).
	 * Evaluates the thermal scaling factor [0., 1.] for yeast metabolic activity.
	 */
	public static double calculateThermalEfficiency(final double temperature){
		if(temperature <= TEMP_MIN || temperature >= TEMP_MAX)
			return 0.;

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
	 * Calculates the sugar consumption during the time step using a linear metabolic yield relation.
	 */
	public static double consumeSugar(final double sugarPrev, final double muBio, final double dt){
		return Math.max(0., sugarPrev - (muBio * SUGAR_CONSUMPTION_COEFF * dt));
	}

}
