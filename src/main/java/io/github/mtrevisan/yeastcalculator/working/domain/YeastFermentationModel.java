package io.github.mtrevisan.yeastcalculator.working.domain;


public final class YeastFermentationModel{

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

	// Biochemical inhibition multipliers (Fixed Osmotic Code Smell constants)
	private static final double SALT_INHIBITION_MULTIPLIER = -15.;
	private static final double OIL_INHIBITION_LINEAR_COEFF = 5.;
	private static final double OIL_INHIBITION_EXP_COEFF = -8.;


	private YeastFermentationModel(){}


	/**
	 * Calculates the osmotic cell-stress inhibition coefficient driven by salt concentration.
	 */
	public static double calculateSaltInhibition(final double salt){
		return Math.exp(SALT_INHIBITION_MULTIPLIER * salt);
	}

	/**
	 * Calculates the structural lipid barrier and cell hydration inhibition coefficient driven by oil concentration.
	 */
	public static double calculateOilInhibition(final double oil){
		return (1. + OIL_INHIBITION_LINEAR_COEFF * oil) * Math.exp(OIL_INHIBITION_EXP_COEFF * oil);
	}

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
	public static double calculateBiomassGrowthRate(final double yeast, final double alphaBio,
			final double lag, final double sugar, final double saltK, final double oilK){
		final double sugarK = sugar / (sugar + SUGAR_AFFINITY_K);
		final double lagTransition = 1. / (1. + Math.exp(-4. * (lag - 0.5)));
		return POTENTIAL_MU_MAX * yeast * alphaBio * saltK * oilK * sugarK * lagTransition;
	}

	/**
	 * Calculates the net sugar rate (dSugar/dt).
	 * Accounts for concurrent enzymatic starch breakdown (generation) and yeast metabolism (consumption).
	 */
	public static double calculateNetSugarRate(final double sugar, final double muBio, final double yeas,
		final double temperature){
		if(sugar <= 0. && muBio <= 0.)
			return 0.;

		final double amylaseThermalK = Math.exp(0.06 * (temperature - 20.))
			* (1. - 0.005 * Math.pow(temperature - 35., 2));
		final double sugarGeneration = AMYLASE_VMAX_BASE * StrictMath.max(0., amylaseThermalK);
		final double sugarConsumption = muBio * YEAST_SUGAR_YIELD_Y + yeas * MAINTENANCE_COEFF_M;
		return sugarGeneration - sugarConsumption;
	}

}
