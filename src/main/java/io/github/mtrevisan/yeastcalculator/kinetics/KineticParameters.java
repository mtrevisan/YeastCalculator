package io.github.mtrevisan.yeastcalculator.kinetics;

/**
 * Global kinetic parameters for S. cerevisiae in dough.
 * Sources: Baranyi (1994), Hamad (2012), Rosso (1995), Ross (1993), Nagodawithana (1986)
 */
public final class KineticParameters{

	// ─ Growth ─
	/**
	 * Maximum specific growth rate [h⁻¹] at T_OPT and aw ≥ AW_OPT_DOUGH,
	 * with no lag and no substrate limitation.
	 * Broth reference from Baranyi (1994); the model applies environmental gamma-factors on top.
	 */
	public static final double MU_MAX_REF = 0.4;

	// ─ Temperature [°C] ─
	/**
	 * Cardinal temperatures [°C] for the dough system — Hamad et al. (2012).
	 * These differ from broth values (T_MIN = 0, T_OPT = 30, T_MAX = 45) because the
	 * dough matrix shifts the optimum and extends the viable range.
	 */
	public static final double TEMPERATURE_CARDINAL_MIN = 2.;
	public static final double TEMPERATURE_CARDINAL_OPT = 32.;
	public static final double TEMPERATURE_CARDINAL_MAX = 48.;

	// ─ Water Activity ─
	/**
	 * Minimum water activity for yeast growth (S. cerevisiae). Ross (1993).
	 * Below this threshold, growth is zero.
	 */
	public static final double AW_MIN = 0.94;
	/**
	 * Optimal water activity for S. cerevisiae in dough — Hamad et al. (2012).
	 * Note: the broth value (0.995) severely underestimates yeast activity in typical
	 * dough (aw ≈ 0.970–0.975). Using the broth value penalizes growth by ~36%.
	 */
	public static final double AW_OPT_DOUGH = 0.975;

	// ─ Ethanol ─
	/**
	 * Ethanol concentration at which growth is fully inhibited [g_EtOH / g_dough].
	 * S. cerevisiae: significant inhibition above ~8% v/v ≈ 0.063 g/g dough.
	 * Nagodawithana et al. (1986).
	 */
	public static final double ETHANOL_INHIBITION_MAX = 0.063;

	// ─ Sugar ─
	/**
	 * Minimum fermentable sugar fraction [g / g_dough] below which fermentation is
	 * considered depleted. Prevents division-by-zero when sugar is exhausted.
	 */
	public static final double SUGAR_DEPLETION_THRESHOLD = 0.001;

	// ─ Yields ─
	/**
	 * Biomass yield on glucose [g_biomass / g_glucose].
	 * Reflects aerobic and anaerobic metabolism mixed in dough.
	 */
	public static final double YIELD_BIOMASS_PER_SUGAR = 0.5;
	/**
	 * CO₂ yield on glucose [g_CO₂ / g_glucose] — stoichiometric.
	 * C₆H₁₂O₆ → 2 CO₂ + 2 C₂H₅OH: 2 × 44 / 180 = 0.4889.
	 */
	public static final double YIELD_CO2_PER_SUGAR = 0.4889;
	/**
	 * Ethanol yield on glucose [g_EtOH / g_glucose] — stoichiometric.
	 * 2 × 46 / 180 = 0.5111.
	 */
	public static final double YIELD_ETHANOL_PER_SUGAR = 0.5111;

	// ─────────────────────────────────────────────────────────────────────────
	// Ross aw model constants — logarithmic multi-solute form
	// ln(aw_i) = a1·m + a2·m²  where m = molality [mol/kg_water]
	// Source: Ross (1993), Chirife & Favetto (1992)
	// ─────────────────────────────────────────────────────────────────────────
	public static final double AW_NACL_COEFF_LINEAR = -0.0808;
	public static final double AW_NACL_COEFF_QUADRATIC = -0.000681;
	public static final double AW_SUCROSE_COEFF_LINEAR = -0.01867;
	public static final double AW_SUCROSE_COEFF_QUADRATIC = -0.000226;

	/** Molecular weight of NaCl [g/mol]. */
	public static final double MOLECULAR_WEIGHT_NACL = 58.44;
	/** Molecular weight of sucrose [g/mol]. */
	public static final double MOLECULAR_WEIGHT_SUCROSE = 342.3;

	// ─────────────────────────────────────────────────────────────────────────
	// Q0 model: physiological state as a function of yeast moisture content
	//
	// Biological rationale:
	//   Trehalose stabilizes dry yeast cells in a glassy state.
	//   Membrane fluidity and enzyme pools are restored progressively as the
	//   cell rehydrates. Below w_b ≈ 0.05 (critical aw ≈ 0.90) the membrane
	//   remains rigid and trehalose mobilization has not begun (Gervais &
	//   Beney (2001)). Above that threshold, recovery follows a power-law as a
	//   function of hydration.
	//
	// Calibration points:
	//   w_b = 0.05 (instant dry yeast, direct use):  Q0 → Q0_ANHYDROUS (cells
	//     hydrate rapidly in dough, but start below the activation threshold)
	//   w_b = 0.70 (fresh compressed yeast):         Q0 = Q0_FRESH_COMPRESSED
	// ─────────────────────────────────────────────────────────────────────────

	/** Q0 for fully anhydrous or near-anhydrous yeast (w_b ≤ WB_ACTIVATION_THRESHOLD). */
	public static final double Q0_ANHYDROUS = 0.02;
	/** Q0 for fresh compressed yeast at w_b = WB_FRESH_COMPRESSED. */
	public static final double Q0_FRESH_COMPRESSED = 2.5;
	/**
	 * Moisture content below which cells remain in near-anabiosis [g_water / g_wet_yeast].
	 * At this point membrane fluidity and trehalose mobilization have not yet begun.
	 */
	public static final double WB_ACTIVATION_THRESHOLD = 0.05;
	/** Moisture content of standard fresh compressed yeast [g_water / g_wet_yeast]. */
	public static final double WB_FRESH_COMPRESSED = 0.7;

	// ─ Amylase ─
	public static final double T_MIN_AMYLASE = 5.;
	public static final double T_OPT_AMYLASE = 52.5;
	public static final double T_MAX_AMYLASE = 75.;
	/** max 3% from starch */
	public static final double SUGAR_GENERATION_MAX = 0.03;
	/** [h⁻¹] — parameter to be calibrated */
	public static final double AMYLASE_ACTIVITY_RATE = 0.15;

	// ─ Solver tolerances ─
	public static final double ODE_RELATIVE_TOLERANCE = 1e-8;
	public static final double ODE_ABSOLUTE_TOLERANCE = 1e-10;
	public static final double BRENT_RELATIVE_TOLERANCE = 1e-9;
	public static final double BRENT_ABSOLUTE_TOLERANCE = 1e-12;
	/** Maximum Brent solver iterations (generous — convergence is fast). */
	public static final int BRENT_MAX_EVALUATIONS = 200;

	/** Target: fermentation reaches 99% of the substrate-limited maximum at t_final. */
	public static final double TARGET_FERMENTATION_PROGRESS = 0.99;

	// ─ State indices ─
	public static final int IDX_PHYSIOLOGICAL_STATE = 0;
	public static final int IDX_BIOMASS_DENSITY = 1;
	public static final int IDX_SUGAR_CONCENTRATION = 2;
	public static final int IDX_CO2_PRODUCED = 3;
	public static final int IDX_ETHANOL_PRODUCED = 4;
	public static final int IDX_SUGAR_ENZYMATIC = 5;
	public static final int STATE_DIMENSION = 6;


	private KineticParameters(){}

}