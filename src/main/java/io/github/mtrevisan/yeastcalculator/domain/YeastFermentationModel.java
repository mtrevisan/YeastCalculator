package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Core biochemical model representing Saccharomyces cerevisiae metabolic pathways.
 * <p>
 * This class hosts the functional equations for thermal efficiency, osmotic inhibitions,
 * amylase starch breakdown, and integrates the Baranyi & Roberts (1994) physiological
 * readiness state equations to govern the lag-to-exponential growth transition.
 * </p>
 */
public final class YeastFermentationModel {

	// Saccharomyces cerevisiae thermal cardinal boundaries (Rosso et al., 1995)
	private static final double TEMP_MIN = 4.;
	private static final double TEMP_OPT = 34.;
	private static final double TEMP_MAX = 42.;

	// Baranyi & Roberts (1994) / Rosso Biological Constants
	public static final double MU_MAX_REF = 0.52; // Maximum specific growth rate reference [h^-1]
	public static final double WB_ACTIVATION_THRESHOLD = 0.05; // Critical moisture below which cells are in anabiosis
	public static final double WB_FRESH_COMPRESSED = 0.70; // Reference moisture for fully active compressed yeast
	public static final double Q0_ANHYDROUS = 0.05; // Initial metabolic readiness pool for dry active cells
	public static final double Q0_FRESH_COMPRESSED = 0.45; // Initial metabolic readiness pool for fresh cells

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
	 * Computes the initial intracellular enzyme pool readiness (Q0) from yeast moisture content.
	 * Models cell anabiosis stabilization by trehalose in dry states (Gervais & Beney 2001).
	 *
	 * @param yeastMoistureFraction Mass fraction of water in the yeast [g_water / g_wet_yeast].
	 * @return Initial physiological state value Q0 [dimensionless].
	 */
	public static double calculateInitialQ0(final double yeastMoistureFraction) {
		if (yeastMoistureFraction <= WB_ACTIVATION_THRESHOLD) {
			return Q0_ANHYDROUS;
		}
		final double clipped = Math.min(yeastMoistureFraction, WB_FRESH_COMPRESSED);
		final double normalizedHydration = (clipped - WB_ACTIVATION_THRESHOLD) / (WB_FRESH_COMPRESSED - WB_ACTIVATION_THRESHOLD);

		// Power-law transition modeling membrane fluidity recovery during rehydration
		return Q0_ANHYDROUS + (Q0_FRESH_COMPRESSED - Q0_ANHYDROUS) * Math.pow(normalizedHydration, 1.5);
	}

	/**
	 * Simulates a warm water rehydration soak before mixing.
	 * Independent of external dough inhibitors, Q grows at its pure maximum reference velocity.
	 *
	 * @param q0Dry                 Initial physiological state computed from moisture.
	 * @param rehydrationHours      Duration of the water soak in hours.
	 * @return Blended effective Q0 coordinate ready for dough integration initiation.
	 */
	public static double activeQAfterRehydration(final double q0Dry, final double rehydrationHours) {
		if (rehydrationHours <= 0.0) return q0Dry;
		return q0Dry * Math.exp(MU_MAX_REF * rehydrationHours);
	}

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
	 * Evaluates the thermal scaling factor [0, 1] for yeast metabolic activity.
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
	 * Computes specific biomass fermentation velocity incorporating the rational Baranyi lag multiplier.
	 * Replaces old empirical sigmoid functions with the standard alpha_lag = Q / (Q + 1).
	 */
	public static double calculateBiomassGrowthRate(final double yeast, final double alphaBio,
		final double qCurr, final double sugar, final double saltK, final double oilK) {
		final double sugarK = sugar / (sugar + SUGAR_AFFINITY_K);

		// Baranyi rational adjustment factor alpha_lag. Highly robust, non-arbitrary.
		final double alphaLag = qCurr / (qCurr + 1.0);

		return POTENTIAL_MU_MAX * yeast * alphaBio * saltK * oilK * sugarK * alphaLag;
	}

	/**
	 * Calculates the net sugar rate (dSugar/dt).
	 * Accounts for concurrent enzymatic starch breakdown (generation) and yeast metabolism (consumption).
	 */
	public static double calculateNetSugarRate(final double sugar, final double muBio, final double yeast,
		final double temperature){
		if(sugar <= 0. && muBio <= 0.)
			return 0.;

		final double amylaseThermalK = Math.exp(0.06 * (temperature - 20.))
			* (1. - 0.005 * Math.pow(temperature - 35., 2));
		final double sugarGeneration = AMYLASE_VMAX_BASE * StrictMath.max(0., amylaseThermalK);
		final double sugarConsumption = muBio * YEAST_SUGAR_YIELD_Y + yeast * MAINTENANCE_COEFF_M;
		return sugarGeneration - sugarConsumption;
	}

}
