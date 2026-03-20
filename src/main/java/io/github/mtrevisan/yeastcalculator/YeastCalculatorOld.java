package io.github.mtrevisan.yeastcalculator;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;


/**
 * <h2>Overview</h2>
 * Predictive microbiology model for <i>Saccharomyces cerevisiae</i> in bread dough.
 * Computes the optimal initial yeast inoculum (g yeast / g dough) so that
 * fermentation reaches 99% of the substrate-limited maximum exactly at the end
 * of the last fermentation stage.
 *
 * <h2>What this model computes</h2>
 * Given a bread dough composition and a fermentation schedule, finds the minimum
 * quantity of <em>anhydrous</em> (dry) yeast [g_dry_yeast / g_dough] such that
 * biomass reaches 99% of the substrate-limited maximum exactly at the end of the
 * last fermentation stage.
 * <p>
 * The result is expressed as <strong>anhydrous</strong> yeast. To convert to the
 * actual wet-yeast quantity to add, divide by {@code (1 − yeastMoistureFraction)}:
 * <pre>
 *   wetYeastFraction = anhydrousYeastFraction / (1 − yeastMoistureFraction)
 * </pre>
 *
 * <h2>ODE system — Baranyi-Roberts (1994)</h2>
 * State vector: [physiologicalState, biomassDensity, sugarConcentration,
 * co2Produced, ethanolProduced]
 * <pre>
 *   dQ/dt = MU_MAX · Q                                       (Baranyi (1994), eq.2)
 *   dN/dt = MU_MAX · γ_T · γ_aw · γ_E · [Q/(Q+1)] · [1 - N/N_max] · N
 *   dS/dt = −(1/Y_XS) · dN/dt
 *   dP/dt = Y_PS · (−dS/dt)
 *   dE/dt = Y_ES · (−dS/dt)
 * </pre>
 *
 * <h2>Environmental factors</h2>
 * <ul>
 *   <li>γ_T  — Cardinal Model, Rosso et al. (1995); dough cardinal points from Hamad et al. (2012)</li>
 *   <li>γ_aw — linear ramp on [AW_MIN, AW_OPT_DOUGH = 0.975]; NOT the broth value 0.995</li>
 *   <li>γ_E  — linear ethanol inhibition, Nagodawithana et al. (1986)</li>
 * </ul>
 *
 * <h2>Q0 from yeast moisture content</h2>
 * Q0 (Baranyi physiological state) is not a property of a "yeast type" label —
 * it reflects how hydrated and metabolically active the cells are at the moment
 * they enter the dough. This model derives Q0 from the yeast moisture content
 * {@code yeastMoistureFraction} (g_water / g_wet_yeast):
 * <ul>
 *   <li>Anhydrous (w_b ≤ 0.05): Q0 = Q0_ANHYDROUS = 0.02 (cells in near-anabiosis)</li>
 *   <li>Fresh compressed (w_b ≈ 0.70): Q0 = Q0_FRESH_COMPRESSED = 2.50</li>
 *   <li>Intermediate: power-law interpolation (calibrated from rehydration kinetics)</li>
 * </ul>
 *
 * <h2>Literature</h2>
 * <ul>
 *   <li>Baranyi, J., & Roberts, T. A. (1994). A dynamic approach to predicting bacterial growth in food. International Journal of Food Microbiology, 23(3–4), 277–294.</li>
 *   <li>Rosso L. et al. (1995) — J Theor Biol 175:159-171</li>
 *   <li>Ross T. (1993) — J Appl Bacteriol 74:608-611</li>
 *   <li>Hamad S. et al. (2012) — Food Microbiol 32:228-236</li>
 *   <li>Nagodawithana T.W. et al. (1986) — Can J Microbiol 32:467-472</li>
 *   <li>Gervais P., Beney L. (2001) — J Appl Microbiol 90:579-585</li>
 * </ul>
 *
 * <h2>Units</h2>
 * All mass fractions are dimensionless (g / g_dough).
 * Time in hours.
 * Temperature in Celsius.
 */
public class YeastCalculatorOld{

	// ─────────────────────────────────────────────────────────────────────────
	// S. cerevisiae kinetic parameters — dough system
	// Sources: Hamad (2012), Rosso (1995), Baranyi (1994)
	// ─────────────────────────────────────────────────────────────────────────

	/**
	 * Maximum specific growth rate [h⁻¹] at T_OPT and aw ≥ AW_OPT_DOUGH,
	 * with no lag and no substrate limitation.
	 * Broth reference from Baranyi (1994); the model applies environmental gamma-factors on top.
	 */
	private static final double MU_MAX_REF = 0.4;

	/**
	 * Cardinal temperatures [°C] for the dough system — Hamad et al. (2012).
	 * These differ from broth values (T_MIN = 0, T_OPT = 30, T_MAX = 45) because the
	 * dough matrix shifts the optimum and extends the viable range.
	 */
	private static final double TEMPERATURE_CARDINAL_MIN = 2.;
	private static final double TEMPERATURE_CARDINAL_OPT = 32.;
	private static final double TEMPERATURE_CARDINAL_MAX = 48.;

	/**
	 * Minimum water activity for yeast growth (S. cerevisiae). Ross (1993).
	 * Below this threshold growth is zero.
	 */
	private static final double AW_MIN = 0.940;

	/**
	 * Optimal water activity for S. cerevisiae in dough — Hamad et al. (2012).
	 * Note: the broth value (0.995) severely underestimates yeast activity in typical
	 * dough (aw ≈ 0.970–0.975). Using the broth value penalizes growth by ~36%.
	 */
	private static final double AW_OPT_DOUGH = 0.975;

	/**
	 * Ethanol concentration at which growth is fully inhibited [g_EtOH / g_dough].
	 * S. cerevisiae: significant inhibition above ~8% v/v ≈ 0.063 g/g dough.
	 * Nagodawithana et al. (1986).
	 */
	private static final double ETHANOL_INHIBITION_MAX = 0.063;

	/**
	 * Minimum fermentable sugar fraction [g / g_dough] below which fermentation is
	 * considered depleted. Prevents division-by-zero when sugar is exhausted.
	 */
	private static final double SUGAR_DEPLETION_THRESHOLD = 0.001;

	// ─────────────────────────────────────────────────────────────────────────
	// Stoichiometric yield coefficients
	// ─────────────────────────────────────────────────────────────────────────

	/**
	 * Biomass yield on glucose [g_biomass / g_glucose].
	 * Reflects aerobic and anaerobic metabolism mix in dough.
	 */
	private static final double YIELD_BIOMASS_PER_SUGAR = 0.50;

	/**
	 * CO₂ yield on glucose [g_CO₂ / g_glucose] — stoichiometric.
	 * C₆H₁₂O₆ → 2 CO₂ + 2 C₂H₅OH: 2 × 44 / 180 = 0.4889.
	 */
	private static final double YIELD_CO2_PER_SUGAR = 0.4889;

	/**
	 * Ethanol yield on glucose [g_EtOH / g_glucose] — stoichiometric.
	 * 2 × 46 / 180 = 0.5111.
	 */
	private static final double YIELD_ETHANOL_PER_SUGAR = 0.5111;

	// ─────────────────────────────────────────────────────────────────────────
	// Ross aw model constants — logarithmic multi-solute form
	// ln(aw_i) = a1·m + a2·m²  where m = molality [mol/kg_water]
	// Source: Ross (1993), Chirife & Favetto (1992)
	// ─────────────────────────────────────────────────────────────────────────

	private static final double AW_NACL_COEFF_LINEAR = -0.0808;
	private static final double AW_NACL_COEFF_QUADRATIC = -0.000681;
	private static final double AW_SUCROSE_COEFF_LINEAR = -0.01867;
	private static final double AW_SUCROSE_COEFF_QUADRATIC = -0.000226;

	/** Molecular weight of NaCl [g/mol]. */
	private static final double MOLECULAR_WEIGHT_NACL = 58.44;
	/** Molecular weight of sucrose [g/mol]. */
	private static final double MOLECULAR_WEIGHT_SUCROSE = 342.3;

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
	private static final double Q0_ANHYDROUS = 0.02;

	/** Q0 for fresh compressed yeast at w_b = WB_FRESH_COMPRESSED. */
	private static final double Q0_FRESH_COMPRESSED = 2.50;

	/**
	 * Moisture content below which cells remain in near-anabiosis [g_water / g_wet_yeast].
	 * At this point membrane fluidity and trehalose mobilization have not yet begun.
	 */
	private static final double WB_ACTIVATION_THRESHOLD = 0.05;

	/** Moisture content of standard fresh compressed yeast [g_water / g_wet_yeast]. */
	private static final double WB_FRESH_COMPRESSED = 0.70;

	// ─────────────────────────────────────────────────────────────────────────
	// GAB (Guggenheim-Anderson-de Boer) flour moisture sorption isotherm
	//
	// W(aw) = Wm · C · Ka / [(1 − Ka) · (1 + (C−1)·Ka)]    Ka = min(K·aw, 0.999)
	//
	// Wm, C, K depend on flour composition (protein, fiber, ash) and temperature
	// via an Arrhenius correction.  Parameters below are for wheat flour; other
	// types are accepted via the constructor.
	//
	// Purpose: part of the recipe water is tightly bound to flour and is NOT
	// available as solvent for the Ross aw calculation.  Using free water only
	// lowers aw by ≈ 0.015–0.020 for typical dough, a physically significant
	// correction (much larger than the yeast-water effect, ≈ 0.001).
	//
	// The GAB reference aw is fixed at AW_GAB_REFERENCE (≈ 0.97) rather than
	// solved implicitly, because dW/d(aw) × 0.01 ≈ 0.003 g/g_dough — the
	// error introduced by the fixed-reference approximation is smaller than the
	// uncertainty in the flour composition parameters.
	//
	// Source: van den Berg & Bruin (1981); Chirife & Iglesias (1992);
	//   original Excel LET formula for Wm, C, K from flour composition.
	// ─────────────────────────────────────────────────────────────────────────

	/**
	 * Reference aw used as input to the GAB isotherm.
	 * Fixed at 0.97 (typical bread dough) to avoid an implicit equation between
	 * aw and bound water.  The error from this approximation is < 0.001 on aw.
	 */
	private static final double AW_GAB_REFERENCE = 0.97;

	// ─────────────────────────────────────────────────────────────────────────
	// Olive oil — kinetic modelling decision
	//
	// Olive oil is immiscible with water and at bread doses (0–10% baker's) has
	// no measurable effect on S. cerevisiae growth kinetics.  It is tracked only
	// as a mass component that dilutes all other dough fractions.
	//
	// Source: reviewed literature (Aleixandre-Tudo et al. 2015; Praphailong &
	//   Fleet 1997) — oleic acid inhibition requires aqueous concentrations
	//   > 0.5 g/L, far above what partitions into the water phase from 5% oil.
	// ─────────────────────────────────────────────────────────────────────────

	// ─────────────────────────────────────────────────────────────────────────
	// Diastatic malt powder (DMP) — fermentable sugar contribution
	//
	// DMP contributes fermentable sugar via two pathways:
	//   a) Endogenous simple sugars present in the malt (glucose, maltose, …):
	//        ≈ MALT_ENDOGENOUS_SUGAR_FRACTION × malt_mass
	//   b) Enzymatic conversion of damaged starch by α/β-amylase during
	//      fermentation:
	//        ≈ MALT_AMYLASE_SUGAR_FRACTION × malt_mass
	//
	// The sum (≈ 0.25 g sugar / g DMP) is added to sugarFraction at
	// construction time.  A time-resolved Michaelis-Menten model is not
	// warranted because starch substrate >> Km in dough (zero-order regime);
	// the total sugar generated over a 2–4 h fermentation at 28–32°C is well
	// approximated by the constant factor below.
	//
	// Source: typical barley DMP specifications; Goesaert et al. (2005)
	//   Trends in Food Science & Technology 16:12-30.
	// ─────────────────────────────────────────────────────────────────────────

	/** Endogenous simple sugars in DMP [g_sugar / g_DMP wet basis]. */
	private static final double MALT_ENDOGENOUS_SUGAR_FRACTION = 0.05;

	/**
	 * Amylase-generated sugars from flour starch per unit DMP
	 * [g_maltose_equiv / g_DMP wet basis].
	 * Calibrated for 2–4 h fermentation at 28–32°C with typical barley DMP
	 * (diastatic power 150–250 °Lintner).  Caller may override via the
	 * {@code maltDiastaticFactor} constructor parameter for non-standard malts.
	 */
	private static final double MALT_AMYLASE_SUGAR_FRACTION = 0.20;

	// ─────────────────────────────────────────────────────────────────────────
	// Solver tolerances
	// ─────────────────────────────────────────────────────────────────────────

	private static final double ODE_RELATIVE_TOLERANCE = 1e-8;
	private static final double ODE_ABSOLUTE_TOLERANCE = 1e-10;
	private static final double BRENT_RELATIVE_TOLERANCE = 1e-9;
	private static final double BRENT_ABSOLUTE_TOLERANCE = 1e-12;

	/** Maximum Brent solver iterations (generous — convergence is fast). */
	private static final int BRENT_MAX_EVALUATIONS = 200;

	/** Target: fermentation reaches 99% of the substrate-limited maximum at t_final. */
	private static final double TARGET_FERMENTATION_PROGRESS = 0.99;

	// ─────────────────────────────────────────────────────────────────────────
	// State vector indices
	// ─────────────────────────────────────────────────────────────────────────

	//Q — Baranyi physiological state
	private static final int IDX_PHYSIOLOGICAL_STATE = 0;
	//N — yeast biomass [g/g_dough]
	private static final int IDX_BIOMASS_DENSITY = 1;
	//S — fermentable sugar [g/g_dough]
	private static final int IDX_SUGAR_CONCENTRATION = 2;
	//P — CO₂ produced [g/g_dough]
	private static final int IDX_CO2_PRODUCED = 3;
	//E — ethanol [g/g_dough]
	private static final int IDX_ETHANOL_PRODUCED = 4;
	private static final int STATE_DIMENSION = 5;


	// ─────────────────────────────────────────────────────────────────────────
	// Dough composition — mass fractions on total dough weight
	// ─────────────────────────────────────────────────────────────────────────

	/** Water fraction [g_water / g_dough] */
	private final double waterFraction;
	/**
	 * Free water available as solvent [g / g_dough].
	 * = waterFraction − boundWaterFraction(GAB).
	 * Used in the Ross aw calculation instead of the total waterFraction,
	 * because bound water on flour is not osmotically active.
	 */
	private final double freeWaterFraction;
	/** Salt fraction [g_salt / g_dough] */
	private final double saltFraction;
	/**
	 * Effective fermentable sugar fraction [g / g_dough].
	 * = recipe sugar + DMP endogenous sugar + DMP amylase-generated sugar.
	 */
	private final double effectiveSugarFraction;

	// ─────────────────────────────────────────────────────────────────────────
	// Fermentation schedule — arrays of equal length, one entry per stage
	// ─────────────────────────────────────────────────────────────────────────

	/** Stage temperatures [°C]. */
	private final double[] stageTemperatures;
	/** Stage durations [h]. */
	private final double[] stageDurations;

	/**
	 * Moisture content of the yeast [g_water / g_wet_yeast].
	 * Used both to derive Q0 and to convert the anhydrous result to wet-yeast mass.
	 */
	private final double yeastMoistureFraction;

	/**
	 * Time the yeast spends in water before being added to the dough [h].
	 * During this phase Q grows at the maximum rate (dQ/dt = MU_MAX_REF · Q),
	 * so the physiological state at the start of dough fermentation is:
	 * <pre>
	 *   Q0_effective = Q0_dry × exp(MU_MAX_REF × rehydrationDurationHours)
	 * </pre>
	 * Use 0 for instant dry yeast added directly to flour, or for fresh compressed
	 * yeast that needs no pre-activation.
	 */
	private final double rehydrationDurationHours;

	/**
	 * Effective initial physiological state at the start of dough fermentation.
	 * Accounts for any pre-activation in water: Q0_eff = Q0_dry × exp(μ_max × t_rehydration).
	 * When rehydrationDurationHours = 0 this equals the raw Q0 from moisture content.
	 */
	private final double initialPhysiologicalState;

	/**
	 * Substrate-limited biomass ceiling [g_biomass / g_dough].
	 * = YIELD_BIOMASS_PER_SUGAR × sugarFraction.
	 * Using a fixed broth value (0.10 g/g) instead produces nonsensical results
	 * (~16% baker's) because the solver always finds N0 near that ceiling.
	 */
	private final double substrateLimitedBiomassCapacity;


	/**
	 * <p>All ingredient fractions are expressed as <strong>g / g_total_dough</strong>
	 * (including oil and malt in the denominator). The caller is responsible for
	 * dividing each ingredient by the total dough weight before passing it here.
	 * The {@link #main} method shows a concrete example.
	 *
	 * @param waterFraction	Water / total dough [g/g].
	 * @param saltFraction	Salt / total dough [g/g].
	 * @param sugarFraction	Recipe fermentable sugar / total dough [g/g]
	 * 	(endogenous flour sugars + added sugar, excluding malt).
	 * @param flourFraction	Flour / total dough [g/g].
	 * 	Used for the GAB bound-water calculation.
	 * @param flourProtein	Flour protein [g/g dry flour].
	 * 	Used for GAB Wm and C parameters.
	 * @param flourFiber	Flour dietary fibre [g/g dry flour].
	 * @param flourAsh	Flour ash [g/g dry flour].
	 * @param stageTemperatures	Stage temperatures [°C].
	 * @param stageDurations	Stage durations [h].
	 * @param yeastMoistureFraction	Moisture of the yeast [g_water / g_wet_yeast]:
	 * 	<ul>
	 * 		<li>≈ 0.04–0.06 — instant/rapid dry yeast (IDY/SDY)</li>
	 * 		<li>≈ 0.07–0.09 — active dry yeast (ADY)</li>
	 * 		<li>≈ 0.68–0.72 — fresh compressed yeast (CY)</li>
	 * 	</ul>
	 * @param rehydrationDurationHours	Time the yeast soaks in warm water before being added to the dough [h]. Use 0
	 * 	for yeast added directly to flour without pre-activation. During rehydration Q grows as
	 * 	Q(t) = Q0_dry × exp(MU_MAX_REF × t), so even a short soak significantly raises Q0 for very dry yeast.
	 * @param oliveoilFraction	Olive oil / total dough [g/g].
	 * 	Tracked for mass balance only; has no effect on yeast kinetics at bread doses (0–10% baker's).
	 * @param maltFraction	Diastatic malt powder / total dough [g/g].
	 * 	Contributes fermentable sugar:
	 * 	≈ (MALT_ENDOGENOUS + MALT_AMYLASE) × maltFraction.
	 * @param maltDiastaticFactor	Override for the amylase sugar contribution
	 * 	[g_sugar / g_DMP].  Pass{@link #MALT_AMYLASE_SUGAR_FRACTION} (0.20) for standard barley DMP (150–250 °Lintner).
	 */
	public YeastCalculatorOld(final double waterFraction, final double saltFraction, final double sugarFraction,
		final double flourFraction, final double flourProtein, final double flourFiber, final double flourAsh,
		final double[] stageTemperatures, final double[] stageDurations, final double yeastMoistureFraction,
		final double rehydrationDurationHours, final double oliveoilFraction, final double maltFraction,
		final double maltDiastaticFactor){
		this.waterFraction = waterFraction;
		this.saltFraction = saltFraction;
		this.stageTemperatures = stageTemperatures;
		this.stageDurations = stageDurations;
		this.yeastMoistureFraction = yeastMoistureFraction;
		this.rehydrationDurationHours = rehydrationDurationHours;

		// ── GAB bound water ────────────────────────────────────────────────
		// Wm, C, K for wheat flour at 20°C from van den Berg & Bruin (1981)
		// with protein/fiber/ash corrections from the LET formula.
		final double wm = 0.0682 + 0.112 * flourProtein + 0.184 * flourFiber;
		final double c = 7.21 + 41.2  * flourProtein + 92.1 * flourAsh - 112.4 * flourFiber / (1. + 8.2 * flourFiber);
		final double k  = 0.8124 - 0.724 * flourProtein - 0.801 * flourFiber - 1.423 * flourAsh;
		final double ka = Math.min(k * AW_GAB_REFERENCE, 0.999);
		final double wGab = wm * c * ka / ((1. - ka) * (1. + (c - 1.) * ka));
		// boundWater [g/g_dough] = W_GAB [g/g_dry_flour] × flourFraction
		final double boundWaterFraction = wGab * flourFraction;
		this.freeWaterFraction = Math.max(0., waterFraction - boundWaterFraction);

		// ── Malt sugar contribution ────────────────────────────────────────
		final double maltSugarContribution = maltFraction * (MALT_ENDOGENOUS_SUGAR_FRACTION + maltDiastaticFactor);
		this.effectiveSugarFraction = sugarFraction + maltSugarContribution;

		// ── Physiological state ────────────────────────────────────────────
		// Q0 of the dry yeast, then advanced by the rehydration soak.
		// dQ/dt = MU_MAX_REF · Q has the exact solution Q(t) = Q0 · exp(MU_MAX_REF · t).
		final double q0Dry = physiologicalStateFromMoisture(yeastMoistureFraction);
		this.initialPhysiologicalState  = q0Dry * Math.exp(MU_MAX_REF * rehydrationDurationHours);

		// ── Biomass substrate ceiling ──────────────────────────────────────
		//substrate-limited capacity: the maximum biomass the available sugar can support
		this.substrateLimitedBiomassCapacity = YIELD_BIOMASS_PER_SUGAR * effectiveSugarFraction;
	}


	// =========================================================================
	// Q0 from yeast moisture content — piecewise power model
	//
	// Below WB_ACTIVATION_THRESHOLD: cells remain in near-anabiosis → Q0_ANHYDROUS.
	// Above that: physiological recovery scales as t^1.5 toward Q0_FRESH_COMPRESSED.
	// =========================================================================
	static double physiologicalStateFromMoisture(final double yeastMoistureFraction){
		if(yeastMoistureFraction <= WB_ACTIVATION_THRESHOLD)
			return Q0_ANHYDROUS;

		final double clipped = Math.min(yeastMoistureFraction, WB_FRESH_COMPRESSED);
		final double normalizedHydration = (clipped - WB_ACTIVATION_THRESHOLD)
			/ (WB_FRESH_COMPRESSED - WB_ACTIVATION_THRESHOLD);
		return Q0_ANHYDROUS + (Q0_FRESH_COMPRESSED - Q0_ANHYDROUS) * Math.pow(normalizedHydration, 1.5);
	}


	// =========================================================================
	// Cardinal temperature model — Rosso et al. (1995)
	// Returns γ_T ∈ [0, 1]. γ_T = 1 at T_OPT.
	// =========================================================================
	private static double cardinalTemperatureFactor(final double temperatureCelsius){
		if(temperatureCelsius <= TEMPERATURE_CARDINAL_MIN || temperatureCelsius >= TEMPERATURE_CARDINAL_MAX)
			return 0.;

		final double numerator = (temperatureCelsius - TEMPERATURE_CARDINAL_MAX)
			* Math.pow(temperatureCelsius - TEMPERATURE_CARDINAL_MIN, 2.);
		final double denominator = (TEMPERATURE_CARDINAL_OPT - TEMPERATURE_CARDINAL_MIN)
			* ((TEMPERATURE_CARDINAL_OPT - TEMPERATURE_CARDINAL_MIN) * (temperatureCelsius - TEMPERATURE_CARDINAL_OPT)
			- (TEMPERATURE_CARDINAL_OPT - TEMPERATURE_CARDINAL_MAX)
			* (TEMPERATURE_CARDINAL_OPT + TEMPERATURE_CARDINAL_MIN - 2. * temperatureCelsius));
		return (denominator == 0. ? 0. : numerator / denominator);
	}

	// =========================================================================
	// Water activity — Ross multi-solute logarithmic model (Ross (1993))
	//
	// Uses freeWaterFraction (= total water − GAB bound water) as the solvent
	// mass, because water tightly bound to flour is not osmotically active.
	// At 65% hydration this lowers aw by ≈ 0.015–0.020 relative to using the
	// total water fraction — a correction larger than any other second-order
	// effect in this model.
	//
	// Solutes: NaCl and sucrose (dominant osmolytes in bread dough).
	// sugarConcentration is the current (dynamic) sugar fraction, because
	// fermentation consumes sugar and changes aw over time.
	// =========================================================================
	private double computeWaterActivity(final double sugarConcentration){
		if(freeWaterFraction <= 0.)
			return 0.9;

		//molality [mol/kg_freeWater] = (mass_solute/mass_freeWater) / (MW [kg/mol])
		final double naclMolality = (saltFraction / freeWaterFraction) / (MOLECULAR_WEIGHT_NACL / 1_000.);
		final double sucroseMolality = (sugarConcentration / freeWaterFraction) / (MOLECULAR_WEIGHT_SUCROSE / 1_000.);

		final double lnAw =
			AW_NACL_COEFF_LINEAR * naclMolality + AW_NACL_COEFF_QUADRATIC * naclMolality * naclMolality
				+ AW_SUCROSE_COEFF_LINEAR * sucroseMolality + AW_SUCROSE_COEFF_QUADRATIC * sucroseMolality
				* sucroseMolality;

		return Math.exp(lnAw);
	}

	// =========================================================================
	// Water activity inhibition factor — γ_aw ∈ [0, 1]
	//
	// Linear ramp between AW_MIN and AW_OPT_DOUGH (Hamad (2012) dough system).
	// Using the broth AW_OPT = 0.995 would penalize typical dough (aw ≈ 0.97)
	// by ~36%, producing unrealistically slow growth.
	// =========================================================================
	private static double waterActivityInhibitionFactor(final double waterActivity){
		if(waterActivity <= AW_MIN)
			return 0.;
		if(waterActivity >= AW_OPT_DOUGH)
			return 1.;
		return (waterActivity - AW_MIN) / (AW_OPT_DOUGH - AW_MIN);
	}


	// =========================================================================
	// Ethanol inhibition factor — γ_E ∈ [0, 1] — Nagodawithana et al. (1986)
	// =========================================================================
	private static double ethanolInhibitionFactor(final double ethanolConcentration){
		return Math.max(0., 1. - ethanolConcentration / ETHANOL_INHIBITION_MAX);
	}


	// =========================================================================
	// Baranyi-Roberts (1994) ODE system — implemented as FirstOrderDifferentialEquations
	//
	// State: [Q, N, S, P_CO2, E_ethanol]
	//
	// Growth equations:
	//   dQ/dt = MU_MAX · Q	← key: constant rate, not mu_eff! Q represents intracellular enzyme reserves that
	//   	replenish at the maximum rate, independently of external conditions - Baranyi (1994), eq. 2)
	//   dN/dt = MU_MAX · γ_T · γ_aw · γ_E · α_lag · α_cap · N
	//   dS/dt = −(1/Y_XS) · dN/dt
	//   dP/dt = Y_PS · (−dS/dt)
	//   dE/dt = Y_ES · (−dS/dt)
	//
	// Where:
	//   α_lag = Q / (Q + 1)                              — lag adjustment → 1 when Q >> 1
	//   α_cap = max(0, 1 − N / N_MAX)                    — capacity adjustment → 0 at saturation
	//   N_MAX = YIELD_BIOMASS_PER_SUGAR × sugarFraction  — substrate-limited ceiling
	// =========================================================================
	private class BaranyiRobertsODE implements FirstOrderDifferentialEquations{

		private final double stageTempCelsius;

		BaranyiRobertsODE(final double stageTempCelsius){
			this.stageTempCelsius = stageTempCelsius;
		}

		@Override
		public int getDimension(){
			return STATE_DIMENSION;
		}

		@Override
		public void computeDerivatives(final double time, final double[] state, final double[] derivatives){
			final double physiologicalState = Math.max(0., state[IDX_PHYSIOLOGICAL_STATE]);
			final double biomassDensity = Math.max(0., state[IDX_BIOMASS_DENSITY]);
			final double sugarConcentration = Math.max(0., state[IDX_SUGAR_CONCENTRATION]);
			final double ethanolConcentration = Math.max(0., state[IDX_ETHANOL_PRODUCED]);

			//── Environmental inhibition factors ──────────────────────────────
			final double waterActivity = computeWaterActivity(sugarConcentration);
			final double gammaTemperature = cardinalTemperatureFactor(stageTempCelsius);
			final double gammaWaterActivity = waterActivityInhibitionFactor(waterActivity);
			final double gammaEthanol = ethanolInhibitionFactor(ethanolConcentration);

			//── Baranyi structural factors ────────────────────────────────────
			//lagAdjustment → 0 during lag phase (Q small), → 1 when fully adapted (Q >> 1)
			final double lagAdjustment = physiologicalState / (physiologicalState + 1.);
			//capacityAdjustment → 0 when biomass saturates the substrate-limited ceiling
			final double capacityAdjustment = Math.max(0., 1. - biomassDensity / substrateLimitedBiomassCapacity);
			//binary sugar availability flag
			final double sugarAvailability = (sugarConcentration > SUGAR_DEPLETION_THRESHOLD ? 1. : 0.);

			//── Effective growth rate [h⁻¹] ──────────────────────────────────
			final double effectiveGrowthRate = MU_MAX_REF
				* gammaTemperature
				* gammaWaterActivity
				* gammaEthanol
				* lagAdjustment
				* capacityAdjustment
				* sugarAvailability;

			//── ODE derivatives ───────────────────────────────────────────────
			//dQ/dt uses MU_MAX_REF (not effectiveGrowthRate): Q evolves at maximum rate regardless of environment — it
			// represents intracellular enzyme reserves that replenish at a fixed rate (Baranyi (1994), eq. 2)
			derivatives[IDX_PHYSIOLOGICAL_STATE] = MU_MAX_REF * physiologicalState;

			final double biomassGrowthRate = effectiveGrowthRate * biomassDensity;
			derivatives[IDX_BIOMASS_DENSITY] = biomassGrowthRate;

			final double sugarConsumptionRate = -(1. / YIELD_BIOMASS_PER_SUGAR) * biomassGrowthRate;
			derivatives[IDX_SUGAR_CONCENTRATION] = sugarConsumptionRate;
			derivatives[IDX_CO2_PRODUCED] = YIELD_CO2_PER_SUGAR * (-sugarConsumptionRate);
			derivatives[IDX_ETHANOL_PRODUCED] = YIELD_ETHANOL_PER_SUGAR * (-sugarConsumptionRate);
		}

	}

	// =========================================================================
	// Clamp state variables to [0, ∞) to prevent the adaptive step-size
	// controller from exploring unphysical regions.
	// =========================================================================
	private static class NonNegativeClampHandler implements StepHandler{
		@Override
		public void init(final double t0, final double[] y0, final double t){}

		@Override
		public void handleStep(final StepInterpolator interpolator, final boolean isLast){
			final double[] state = interpolator.getInterpolatedState();
			for(int i = 0; i < state.length; i ++)
				state[i] = Math.max(0., state[i]);
		}

	}


	// =========================================================================
	// Integrate the ODE across all fermentation stages.
	//
	// @param initialAnhydrousYeastFraction [g_dry_yeast / g_dough]
	// @return	Final ODE state vector [Q, N, S, P_CO2, E_EtOH].
	// =========================================================================
	private double[] simulate(final double initialAnhydrousYeastFraction){
		final double[] state = new double[STATE_DIMENSION];
		state[IDX_PHYSIOLOGICAL_STATE] = initialPhysiologicalState;
		state[IDX_BIOMASS_DENSITY] = initialAnhydrousYeastFraction;
		state[IDX_SUGAR_CONCENTRATION] = effectiveSugarFraction;
		state[IDX_CO2_PRODUCED] = 0.;
		state[IDX_ETHANOL_PRODUCED] = 0.;

		double stageStartTime = 0.;
		for(int stage = 0; stage < stageTemperatures.length; stage ++){
			final double stageDuration = stageDurations[stage];
			if(stageDuration <= 0.)
				continue;

			final DormandPrince853Integrator integrator = new DormandPrince853Integrator(
				//minimum adaptive step size [h] — prevents stalling
				1e-6,
				//maximum adaptive step size [h] — ODE solver may use fewer steps
				stageDuration,
				ODE_ABSOLUTE_TOLERANCE,
				ODE_RELATIVE_TOLERANCE);

			integrator.addStepHandler(new NonNegativeClampHandler());
			integrator.integrate(
				new BaranyiRobertsODE(stageTemperatures[stage]),
				stageStartTime, state,
				stageStartTime + stageDuration, state);

			//hard clamp after each stage (belt-and-suspenders)
			for(int i = 0; i < STATE_DIMENSION; i ++)
				state[i] = Math.max(0., state[i]);

			stageStartTime += stageDuration;
		}
		return state;
	}


	// =========================================================================
	// Fermentation progress function
	//
	// progress(N0) = N(t_final) / N_MAX_substrate ∈ [0, 1]
	//
	// Monotone increasing in initialAnhydrousYeastFraction for values in (0, substrateLimitedBiomassCapacity):
	//   - More initial yeast → faster approach to substrate ceiling → higher progress.
	//   - Using N_MAX_substrate (not a fixed broth value) ensures the bracket is
	//     valid and results are physically realistic.
	// =========================================================================
	private double computeFermentationProgress(final double initialAnhydrousYeastFraction){
		return simulate(initialAnhydrousYeastFraction)[IDX_BIOMASS_DENSITY]
			/ substrateLimitedBiomassCapacity;
	}


	// =========================================================================
	// Estimated lag time in dough [h] — Baranyi & Roberts (1994), eq.(4)
	//
	// λ_dough = ln(1 + 1/Q0_effective) / MU_MAX_REF
	//
	// Q0_effective already incorporates any pre-dough rehydration soak
	// (see constructor). For the total wait from the start of rehydration:
	//   λ_total = rehydrationDurationHours + λ_dough
	// Rationale: dQ/dt = MU_MAX_REF · Q (constant, environment-independent).
	// Q exits the lag state — i.e. Q/(Q+1) → 1 — after time
	//   λ = ln(1 + 1/Q0) / MU_MAX_REF
	// regardless of temperature or water activity, because those factors only
	// affect the BIOMASS growth rate dN/dt, not the physiological recovery dQ/dt.
	//
	// The denominator is MU_MAX_REF, not μ_eff: dQ/dt = MU_MAX_REF · Q is
	// environment-independent (Baranyi 1994, eq. 2). Temperature and aw only
	// affect dN/dt, not dQ/dt, so the lag exit time is temperature-independent.
	// λ_dough is therefore a lower bound on visible leavening; the baker will
	// observe significant dough rise somewhat later, when N has grown enough.
	//
	// An infinite value is returned when stage-1 temperature is outside the
	// viable range [T_MIN, T_MAX].
	// =========================================================================
	public double estimatedLagTimeInDoughHours(){
		if(cardinalTemperatureFactor(stageTemperatures[0]) <= 0.)
			return Double.POSITIVE_INFINITY;
		return Math.log(1. + 1. / initialPhysiologicalState) / MU_MAX_REF;
	}


	// =========================================================================
	// Main entry point: find the optimal initial anhydrous yeast fraction
	//
	// Uses BrentSolver to find:
	//   N0* = argmin { |progress(N0) − TARGET_FERMENTATION_PROGRESS| }
	//       = root of f(N0) = progress(N0) − TARGET_FERMENTATION_PROGRESS
	//
	// Search domain: (0, N_MAX_substrate)
	//   Lower bound: effectively zero yeast — 1e-8 g/g
	//   Upper bound: substrate ceiling — N_MAX_substrate × 0.9999 (avoids capacity_adj = 0)
	//
	// @return	Optimal initial anhydrous yeast fraction [g_dry_yeast / g_dough].
	//          Use {@link #toWetYeastFraction} or {@link #toBakersWetPercent}
	//          to convert to the wet-yeast quantity to weigh on a scale.
	// =========================================================================
	public double findOptimalAnhydrousYeastFraction(){
		final String warning = scheduleWarning();
		if(warning != null)
			System.err.println("WARNING: " + warning);

		final double lowerBound = 1e-8;
		final double upperBound = substrateLimitedBiomassCapacity * 0.9999;

		final double progressAtLowerBound = computeFermentationProgress(lowerBound);
		final double progressAtUpperBound = computeFermentationProgress(upperBound);

		//check whether the bracket is valid for BrentSolver
		if(progressAtUpperBound < TARGET_FERMENTATION_PROGRESS){
			System.err.printf(Locale.US,
				"WARNING: even at the substrate ceiling N0 = %.6f g/g (%.4f%% baker's dry)%n"
					+ "  progress = %.4f < %.2f.%n"
					+ "  Schedule is too short for this yeast type (w_b = %.2f, Q0 = %.4f).%n"
					+ "  Returning the substrate ceiling as a conservative fallback.%n",
				upperBound, toBakersDryPercent(upperBound),
				progressAtUpperBound, TARGET_FERMENTATION_PROGRESS,
				yeastMoistureFraction, initialPhysiologicalState);
			return upperBound;
		}

		if(progressAtLowerBound >= TARGET_FERMENTATION_PROGRESS){
			System.err.printf("WARNING: near-zero yeast already reaches the target. Stages are too long.%n");
			return lowerBound;
		}

		//f(N0) = progress(N0) − TARGET is monotone increasing on [lo, hi]
		final UnivariateFunction objectiveFunction
			= N0 -> computeFermentationProgress(N0) - TARGET_FERMENTATION_PROGRESS;

		final BrentSolver solver = new BrentSolver(BRENT_RELATIVE_TOLERANCE, BRENT_ABSOLUTE_TOLERANCE);
		try{
			return solver.solve(BRENT_MAX_EVALUATIONS, objectiveFunction, lowerBound, upperBound);
		}
		catch(final NoBracketingException exception){
			System.err.printf(Locale.US,
				"BrentSolver failed: no sign change in [%.2e, %.2e]: %s%n",
				lowerBound, upperBound, exception.getMessage());
			return upperBound;
		}
	}

	// =========================================================================
	// Schedule feasibility check — Baranyi (1994): λ ≈ ln(1 + 1/Q0) / μ_max_eff
	//
	// Returns a warning string when the total estimated wait (rehydration +
	// lag in dough) exceeds 50% of the total fermentation schedule, or null.
	// Returning the message instead of printing it lets callers embed it inline
	// in tables — avoids interleaved stderr/stdout noise when called in loops.
	// =========================================================================
	private String scheduleWarning(){
		final double lagDoughHours = estimatedLagTimeInDoughHours();
		final double totalDoughHours = sumArray(stageDurations);

		if(Double.isInfinite(lagDoughHours))
			return String.format(
				"stage-1 temperature (%.1f°C) is outside the viable range [%.1f, %.1f°C].",
				stageTemperatures[0], TEMPERATURE_CARDINAL_MIN, TEMPERATURE_CARDINAL_MAX);

		if(lagDoughHours > 0.5 * totalDoughHours)
			return String.format(
				"total wait (lag %.2f h = %.2f h) exceeds 50%% of dough time (%.2f h) "
					+ "for w_b=%.2f — consider a longer schedule or more pre-hydrated yeast.",
				rehydrationDurationHours * 60., lagDoughHours, totalDoughHours,
				yeastMoistureFraction);

		return null;
	}


	// =========================================================================
	// Conversion helpers
	// =========================================================================

	/**
	 * Converts an anhydrous mass fraction to baker's percentage — dry-yeast basis.
	 * Result: g_dry_yeast / 100 g_flour.
	 */
	public double toBakersDryPercent(final double anhydrousMassFraction){
		final double flourFraction = 1. - waterFraction - saltFraction - effectiveSugarFraction;
		if(flourFraction <= 0.)
			throw new ArithmeticException("Flour fraction is zero or negative.");

		return anhydrousMassFraction / flourFraction * 100.;
	}

	/**
	 * Converts an anhydrous mass fraction [g/g_dough] to the equivalent
	 * wet-yeast fraction [g_wet_yeast / g_dough]:
	 * <pre>
	 *   wetFraction = anhydrousFraction / (1 − yeastMoistureFraction)
	 * </pre>
	 */
	public double toWetYeastFraction(final double anhydrousMassFraction){
		if(yeastMoistureFraction >= 1.)
			throw new ArithmeticException("yeastMoistureFraction must be < 1.");

		return anhydrousMassFraction / (1. - yeastMoistureFraction);
	}

	/**
	 * Converts an anhydrous mass fraction to baker's percentage — wet-yeast basis.
	 * Result: g_wet_yeast / 100 g_flour.
	 */
	public double toBakersWetPercent(final double anhydrousMassFraction){
		return toBakersDryPercent(toWetYeastFraction(anhydrousMassFraction));
	}


	// =========================================================================
	// Diagnostic trace: print time-series of key state variables
	// =========================================================================
	public void printSimulationTrace(final double initialAnhydrousYeastFraction){
		final int PRINT_POINTS_PER_STAGE = 20;

		final double[] state = new double[STATE_DIMENSION];
		state[IDX_PHYSIOLOGICAL_STATE] = initialPhysiologicalState;
		state[IDX_BIOMASS_DENSITY] = initialAnhydrousYeastFraction;
		state[IDX_SUGAR_CONCENTRATION] = effectiveSugarFraction;
		state[IDX_CO2_PRODUCED] = 0.;
		state[IDX_ETHANOL_PRODUCED] = 0.;

		System.out.printf(Locale.US,
			"%-8s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-8s%n",
			"time[h]", "stage", "Q", "N", "S", "CO2", "EtOH", "aw", "μ_eff");
		System.out.println("─".repeat(100));

		double stageStartTime = 0.;
		for(int stage = 0; stage < stageTemperatures.length; stage ++){
			final double stageDuration = stageDurations[stage];
			final double printInterval = stageDuration / PRINT_POINTS_PER_STAGE;

			for(int printPoint = 0; printPoint < PRINT_POINTS_PER_STAGE; printPoint ++){
				final double currentTime = stageStartTime + printPoint * printInterval;
				final double nextTime = Math.min(currentTime + printInterval, stageStartTime + stageDuration);

				//print state at currentTime before advancing
				final double physiologicalState = state[IDX_PHYSIOLOGICAL_STATE];
				final double biomassDensity = state[IDX_BIOMASS_DENSITY];
				final double sugarConcentration = state[IDX_SUGAR_CONCENTRATION];

				final double waterActivity = computeWaterActivity(sugarConcentration);
				final double sugarAvailability = (sugarConcentration > SUGAR_DEPLETION_THRESHOLD ? 1. : 0.);
				final double effectiveGrowthRate = MU_MAX_REF
					* cardinalTemperatureFactor(stageTemperatures[stage])
					* waterActivityInhibitionFactor(waterActivity)
					* ethanolInhibitionFactor(state[IDX_ETHANOL_PRODUCED])
					* (physiologicalState / (physiologicalState + 1.))
					* Math.max(0., 1. - biomassDensity / substrateLimitedBiomassCapacity)
					* sugarAvailability;

				System.out.printf(Locale.US,
					"%-8.2f %-6d %-6.2f %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-6.3f %-1.4f%n",
					currentTime, stage + 1,
					physiologicalState, biomassDensity * 100, sugarConcentration * 100,
					state[IDX_CO2_PRODUCED] * 100, state[IDX_ETHANOL_PRODUCED] * 100,
					waterActivity, effectiveGrowthRate);

				//advance integrator to next print point
				final DormandPrince853Integrator integrator = new DormandPrince853Integrator(
					1e-6, printInterval, ODE_ABSOLUTE_TOLERANCE, ODE_RELATIVE_TOLERANCE);
				integrator.integrate(
					new BaranyiRobertsODE(stageTemperatures[stage]),
					currentTime, state,
					nextTime, state);
				for(int i = 0; i < STATE_DIMENSION; i ++)
					state[i] = Math.max(0., state[i]);
			}
			stageStartTime += stageDuration;
		}

		System.out.printf(Locale.US,
			"%-8.2f %-6s %-6.2f %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-6.3f %-1s%n",
			stageStartTime, "FINAL",
			state[IDX_PHYSIOLOGICAL_STATE], state[IDX_BIOMASS_DENSITY] * 100,
			state[IDX_SUGAR_CONCENTRATION] * 100, state[IDX_CO2_PRODUCED] * 100,
			state[IDX_ETHANOL_PRODUCED] * 100, computeWaterActivity(state[IDX_SUGAR_CONCENTRATION]), "—");
	}


	private static double sumArray(final double[] values){
		double sum = 0.;
		for(final double value : values)
			sum += value;
		return sum;
	}


	// =========================================================================
	// Example usage
	// =========================================================================
	public static void main(final String[] args){

		//── Dough composition ──────────────────────────────────────────────
		//total = flour + water + salt + sugar
		final double flourUnit = 1.;
		final double waterUnit = 0.65;
		final double saltUnit = 0.012;
		//endogenous + added sugar (not malt)
		final double sugarUnit = 0.02;
		final double oilUnit = 0.03;
		final double maltUnit = 0.01;
		final double totalDoughWeight = flourUnit + waterUnit + saltUnit + sugarUnit + oilUnit + maltUnit;

		final double waterFrac = waterUnit / totalDoughWeight;
		final double saltFrac = saltUnit / totalDoughWeight;
		final double sugarFrac = sugarUnit / totalDoughWeight;
		final double flourFrac  = flourUnit / totalDoughWeight;
		final double oilFrac = oilUnit / totalDoughWeight;
		final double maltFrac = maltUnit / totalDoughWeight;

		// ── Flour composition (typical wheat, dry basis) ───────────────────
		final double flourProtein = 0.12;
		final double flourFiber = 0.025;
		final double flourAsh = 0.006;

		//── Fermentation schedule ──────────────────────────────────────────
		//two-stage: bulk fermentation + proofing
		final double[] temperatures = {28., 30.};
		final double[] durations = {2., 1.5};

		//── Yeast: fresh compressed ────────────────────────────────────────
		//the model computes Q0 from this value and returns an anhydrous mass
		final double yeastMoisture = 0.7;
		// rehydrationHours = time spent soaking in warm water before mixing
		final double rehydrationHours = 20. / 60.;

		//── Instantiate model ──────────────────────────────────────────────
		final YeastCalculatorOld model = new YeastCalculatorOld(
			waterFrac, saltFrac, sugarFrac,
			flourFrac, flourProtein, flourFiber, flourAsh,
			temperatures, durations,
			yeastMoisture, rehydrationHours,
			oilFrac, maltFrac, MALT_AMYLASE_SUGAR_FRACTION);

		//── Optimize ──────────────────────────────────────────────────────
		final double optimalAnhydrousYeast = model.findOptimalAnhydrousYeastFraction();

		System.out.printf("%n=== RESULT ===%n");
		System.out.printf(Locale.US,
			"Yeast moisture (w_b) : %.1f%%%n", yeastMoisture * 100);
		System.out.printf(Locale.US,
			"Free water fraction     : %.3f  (total %.3f, bound %.3f)%n",
			model.freeWaterFraction, waterFrac, waterFrac - model.freeWaterFraction);
		System.out.printf(Locale.US,
			"Rehydration soak     : %.2f hrs (Q0_dry=%.2f, Q0_eff=%.2f)%n",
			rehydrationHours,
			physiologicalStateFromMoisture(yeastMoisture),
			model.initialPhysiologicalState);
		System.out.printf(Locale.US,
			"Lag in dough         : %d min%n", (int)Math.ceil(model.estimatedLagTimeInDoughHours() * 60));
		System.out.printf(Locale.US,
			"Optimal yeast (WB)   : %.2f%%%n",
			model.toBakersWetPercent(optimalAnhydrousYeast));

		//── Diagnostic trace ──────────────────────────────────────────────
		System.out.printf("%n=== SIMULATION TRACE ===%n");
		model.printSimulationTrace(optimalAnhydrousYeast);

		//── Sensitivity analysis ──────────────────────────────────────────
		System.out.printf("%n=== SENSITIVITY ===%n");
		System.out.printf(Locale.US,
			"%-7s %-7s %-10s %-1s%n",
			"factor", "N0(DB)", "yeast(WB)", "progress");
		for(final double factor : new double[]{0.25, 0.50, 0.75, 1.00, 1.50, 2.00}){
			final double testFraction = Math.min(
				optimalAnhydrousYeast * factor,
				model.substrateLimitedBiomassCapacity * 0.9999);
			System.out.printf(Locale.US,
				"%-7.2f %-1.2f%%   %-1.2f%%      %-1.2f%n",
				factor, testFraction * 100,
				model.toBakersWetPercent(testFraction),
				model.computeFermentationProgress(testFraction));
		}

		// ── Effect of olive oil (varying %, all else fixed) ───────────────
		System.out.printf("%n=== EFFECT OF OLIVE OIL ===%n");
		System.out.printf("%-12s %-14s %-12s%n", "Oil%(baker)", "Yeast%(WB)", "ΔYeast%");
		final double yeastBase = model.toBakersWetPercent(optimalAnhydrousYeast);
		for(final double oilPct : new double[]{0., 2., 3., 5., 10.}){
			final double oilU2 = oilPct / 100.;
			final double total2 = flourUnit + waterUnit + saltUnit + sugarUnit + oilU2 + maltUnit;
			final YeastCalculatorOld mo = new YeastCalculatorOld(
				waterUnit / total2, saltUnit / total2, sugarUnit / total2,
				flourUnit / total2, flourProtein, flourFiber, flourAsh,
				temperatures, durations, yeastMoisture, rehydrationHours,
				oilU2 / total2, maltUnit / total2, MALT_AMYLASE_SUGAR_FRACTION);
			final boolean hw = mo.scheduleWarning() != null;
			if(hw)
				System.setErr(new PrintStream(OutputStream.nullOutputStream()));
			if(hw)
				System.setErr(System.err);
			final double bpo = mo.toBakersWetPercent(mo.findOptimalAnhydrousYeastFraction());
			System.out.printf(Locale.US,
				"%-12.1f %-14.4f %-12.4f%n",
				oilPct, bpo, bpo - yeastBase);
		}

		// ── Effect of diastatic malt powder ───────────────────────────────
		System.out.printf("%n=== EFFECT OF DIASTATIC MALT POWDER ===%n");
		System.out.printf("%-12s %-14s %-14s %-12s%n",
			"Malt%(baker)", "Yeast%(WB)", "+Sugar%", "ΔYeast%");
		for(final double maltPct : new double[]{0., 0.5, 1., 2., 5.}){
			final double maltU2 = maltPct / 100.;
			final double total2 = flourUnit + waterUnit + saltUnit + sugarUnit + oilUnit + maltU2;
			final double mf2 = maltU2 / total2;
			final YeastCalculatorOld mm = new YeastCalculatorOld(
				waterUnit / total2, saltUnit / total2, sugarUnit / total2,
				flourUnit / total2, flourProtein, flourFiber, flourAsh,
				temperatures, durations, yeastMoisture, rehydrationHours,
				oilUnit / total2, mf2, MALT_AMYLASE_SUGAR_FRACTION);
			final boolean hw = mm.scheduleWarning() != null;
			final java.io.PrintStream se = System.err;
			if(hw)
				System.setErr(new java.io.PrintStream(java.io.OutputStream.nullOutputStream()));
			final double n0m = mm.findOptimalAnhydrousYeastFraction();
			if(hw)
				System.setErr(se);
			final double bpm = mm.toBakersWetPercent(n0m);
			final double sugarAdded = mf2 * (MALT_ENDOGENOUS_SUGAR_FRACTION + MALT_AMYLASE_SUGAR_FRACTION)
				/ (flourUnit / total2) * 100.;
			System.out.printf(Locale.US,
				"%-12.1f %-14.4f %-14.4f %-12.4f%n",
				maltPct, bpm, sugarAdded, bpm - yeastBase);
		}
	}

}
