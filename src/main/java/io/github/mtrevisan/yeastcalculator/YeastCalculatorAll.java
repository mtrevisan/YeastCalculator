package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.kinetics.KineticParameters;
import io.github.mtrevisan.yeastcalculator.model.FlourMoistureModel;
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
 * Computes the optimal initial yeast inoculum (g yeast / g dough) so that fermentation reaches 99% of the
 * substrate-limited maximum exactly at the end of the last fermentation stage.
 *
 * <h2>What this model computes</h2>
 * Given a bread dough composition and a fermentation schedule, finds the minimum quantity of <em>anhydrous</em> (dry)
 * yeast [g_dry_yeast / g_dough] such that biomass reaches 99% of the substrate-limited maximum exactly at the end of
 * the last fermentation stage.
 * <p>
 * The result is expressed as <strong>anhydrous</strong> yeast. To convert to the actual wet-yeast quantity to add,
 * divide by {@code (1 − yeastMoistureFraction)}:
 * <pre>
 *   wetYeastFraction = anhydrousYeastFraction / (1 − yeastMoistureFraction)
 * </pre>
 *
 * <h2>ODE system — Baranyi-Roberts (1994)</h2>
 * State vector: [physiologicalState, biomassDensity, sugarConcentration, co2Produced, ethanolProduced, sugarEnzymatic]
 * <pre>
 *   dQ/dt = MU_MAX · Q   (Baranyi (1994), eq.2)
 *   dN/dt = MU_MAX · γ_T · γ_aw · γ_E · [Q/(Q+1)] · [1 - N/N_max] · N
 *   dS/dt = −(1/Y_XS) · dN/dt
 *   dP/dt = Y_PS · (−dS/dt)
 *   dE/dt = Y_ES · (−dS/dt)
 *   dS_enz/dt = k_amylase · (1 − S_enz/S_max) · γ_T_amylase
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
 * Q0 (Baranyi physiological state) is not a property of a "yeast type" label — it reflects how hydrated and
 * metabolically active the cells are at the moment they enter the dough. This model derives Q0 from the yeast moisture
 * content
 * {@code yeastMoistureFraction} (g_water / g_wet_yeast):
 * <ul>
 *   <li>Anhydrous (w_b ≤ 0.05): Q0 = Q0_ANHYDROUS = 0.02 (cells in near-anabiosis)</li>
 *   <li>Fresh compressed (w_b ≈ 0.70): Q0 = Q0_FRESH_COMPRESSED = 2.50</li>
 *   <li>Intermediate: power-law interpolation (calibrated from rehydration kinetics)</li>
 * </ul>
 *
 * <h2>Literature</h2>
 * <ul>
 *   <li>Baranyi, J., & Roberts, T. A. (1994). A dynamic approach to predicting bacterial growth in food. International Journal of Food Microbiology, 23 (3–4), 277–294.</li>
 *   <li>Rosso L. et al. (1995) — J Theor Biol 175:159-171</li>
 *   <li>Ross T. (1993) — J Appl Bacteriol 74:608-611</li>
 *   <li>Hamad S. et al. (2012) — Food Microbiol 32:228-236</li>
 *   <li>Nagodawithana T. W. et al. (1986) — Can J Microbiol 32:467-472</li>
 *   <li>Gervais P., Beney L. (2001) — J Appl Microbiol 90:579-585</li>
 * </ul>
 *
 * <h2>Units</h2>
 * All mass fractions are dimensionless (g / g_dough).
 * Time in hours.
 * Temperature in Celsius.
 */
public class YeastCalculatorAll{

	// ─────────────────────────────────────────────────────────────────────────
	// Dough composition — mass fractions on total dough weight
	// ─────────────────────────────────────────────────────────────────────────

	/** Water fraction [g_water / g_dough] */
	private final double water;
	/** Salt fraction [g_salt / g_dough] */
	private final double salt;
	/**
	 * Effective fermentable sugar fraction [g / g_dough].
	 * = recipe sugar + DMP endogenous sugar + DMP amylase-generated sugar.
	 */
	private final double sugar;

	// ─────────────────────────────────────────────────────────────────────────
	// Fermentation schedule — arrays of equal length, one entry per stage
	// ─────────────────────────────────────────────────────────────────────────

	/** Stage temperatures [°C]. */
	private final double[] stageTemperatures;
	/** Stage durations [h]. */
	private final double[] stageDurations;

	/**
	 * Free water content of the flour [g_water / g_wet_flour].
	 * = waterFraction − boundWaterFraction(GAB).
	 * Used in the Ross aw calculation instead of the total waterFraction,
	 * because bound water on flour is not osmotically active.
	 */
	private final double flourFreeWater;
	/**
	 * Moisture content of the yeast [g_water / g_wet_yeast].
	 * Used both to derive Q0 and to convert the anhydrous result to wet-yeast mass.
	 */
	private final double yeastMoisture;

	/**
	 * Time the yeast spends in water before being added to the dough [h].
	 * During this phase Q grows at the maximum rate (dQ/dt = MU_MAX_REF · Q),
	 * so the physiological state at the start of dough fermentation is:
	 * <pre>
	 *   Q0_effective = Q0_dry × exp(MU_MAX_REF × rehydrationDurationHours)
	 * </pre>
	 * Use 0 for instant dry yeast added directly to flour or for fresh compressed
	 * yeast that needs no pre-activation.
	 */
	private final double rehydrationDuration;

	/**
	 * Effective initial physiological state at the start of dough fermentation.
	 * Accounts for any pre-activation in water: Q0_eff = Q0_dry × exp(μ_max × t_rehydration).
	 * When rehydrationDurationHours is equal to zero, this equals the raw Q0 from moisture content.
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
	 * @param water	Water [baker's percentage].
	 * @param salt	Salt [baker's percentage].
	 * @param sugar	Recipe fermentable sugar [baker's percentage]
	 * 	(endogenous flour sugars plus added sugar, excluding malt).
	 * @param stageTemperatures	Stage temperatures [°C].
	 * @param stageDurations	Stage durations [h].
	 * @param yeastMoisture	Moisture of the yeast [g_water / g_wet_yeast]:
	 * 	<ul>
	 * 		<li>≈ 0.04–0.06 — instant/rapid dry yeast (IDY/SDY)</li>
	 * 		<li>≈ 0.07–0.09 — active dry yeast (ADY)</li>
	 * 		<li>≈ 0.68–0.72 — fresh compressed yeast (CY)</li>
	 * 	</ul>
	 * @param rehydrationDuration	Time the yeast soaks in warm water before being added to the dough [h]. Use 0
	 * 	for yeast added directly to flour without pre-activation. During rehydration Q grows as
	 * 	Q(t) = Q0_dry × exp(MU_MAX_REF × t), so even a short soak significantly raises Q0 for very dry yeast.
	 */
	public YeastCalculatorAll(final double water, final double salt, final double sugar,
			final double[] stageTemperatures, final double[] stageDurations,
			final double flourFreeWater, final double yeastMoisture,
			final double rehydrationDuration){
		final double totalDoughWeight = 1. + water + salt + sugar;
		this.water = water / totalDoughWeight;
		this.salt = salt / totalDoughWeight;
		this.sugar = sugar / totalDoughWeight;
		this.stageTemperatures = stageTemperatures;
		this.stageDurations = stageDurations;
		this.flourFreeWater = flourFreeWater / totalDoughWeight;
		this.yeastMoisture = yeastMoisture;
		this.rehydrationDuration = rehydrationDuration;

		// ── Physiological state ────────────────────────────────────────────
		// Q0 of the dry yeast, then advanced by the rehydration soak.
		// dQ/dt = MU_MAX_REF · Q has the exact solution Q(t) = Q0 · exp(MU_MAX_REF · t).
		final double q0Dry = physiologicalStateFromMoisture(yeastMoisture);
		this.initialPhysiologicalState = q0Dry * Math.exp(KineticParameters.MU_MAX_REF * rehydrationDuration);

		// ── Biomass substrate ceiling ──────────────────────────────────────
		//substrate-limited capacity: the maximum biomass the available sugar can support
		this.substrateLimitedBiomassCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * this.sugar;
	}


	// =========================================================================
	// Q0 from yeast moisture content — piecewise power model
	//
	// Below WB_ACTIVATION_THRESHOLD: cells remain in near-anabiosis → Q0_ANHYDROUS.
	// Above that: physiological recovery scales as t^1.5 toward Q0_FRESH_COMPRESSED.
	// =========================================================================
	static double physiologicalStateFromMoisture(final double yeastMoistureFraction){
		if(yeastMoistureFraction <= KineticParameters.WB_ACTIVATION_THRESHOLD)
			return KineticParameters.Q0_ANHYDROUS;

		final double clipped = Math.min(yeastMoistureFraction, KineticParameters.WB_FRESH_COMPRESSED);
		final double normalizedHydration = (clipped - KineticParameters.WB_ACTIVATION_THRESHOLD)
			/ (KineticParameters.WB_FRESH_COMPRESSED - KineticParameters.WB_ACTIVATION_THRESHOLD);
		return KineticParameters.Q0_ANHYDROUS
			+ (KineticParameters.Q0_FRESH_COMPRESSED - KineticParameters.Q0_ANHYDROUS) * Math.pow(normalizedHydration, 1.5);
	}


	// =========================================================================
	// Cardinal temperature model — Rosso et al. (1995)
	// Returns γ_T ∈ [0, 1]. γ_T = 1 at T_OPT.
	// =========================================================================
	private static double cardinalTemperatureFactor(final double temperature, final double Tmin, final double Topt,
			final double Tmax){
		if(temperature <= Tmin || temperature >= Tmax)
			return 0.;

		final double numerator = (temperature - Tmax) * Math.pow(temperature - Tmin, 2.);
		final double denominator = (Topt - Tmin) * ((Topt - Tmin) * (temperature - Topt)
			- (Topt - Tmax) * (Topt + Tmin - 2. * temperature));
		return (denominator == 0.? 0.: numerator / denominator);
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
		final double totalWater = water + flourFreeWater;
		if(totalWater <= 0.)
			return 0.9;

		//molality [mol/kg_freeWater] = (mass_solute/mass_freeWater) / (MW [kg/mol])
		final double naclMolality = (salt / totalWater)
			/ (KineticParameters.MOLECULAR_WEIGHT_NACL / 1_000.);
		final double sucroseMolality = (sugarConcentration / totalWater)
			/ (KineticParameters.MOLECULAR_WEIGHT_SUCROSE / 1_000.);

		final double lnAw = KineticParameters.AW_NACL_COEFF_LINEAR * naclMolality
			+ KineticParameters.AW_NACL_COEFF_QUADRATIC * naclMolality * naclMolality
			+ KineticParameters.AW_SUCROSE_COEFF_LINEAR * sucroseMolality
			+ KineticParameters.AW_SUCROSE_COEFF_QUADRATIC * sucroseMolality * sucroseMolality;

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
		if(waterActivity <= KineticParameters.AW_MIN)
			return 0.;
		if(waterActivity >= KineticParameters.AW_OPT_DOUGH)
			return 1.;
		return (waterActivity - KineticParameters.AW_MIN) / (KineticParameters.AW_OPT_DOUGH - KineticParameters.AW_MIN);
	}


	// =========================================================================
	// Ethanol inhibition factor — γ_E ∈ [0, 1] — Nagodawithana et al. (1986)
	// =========================================================================
	private static double ethanolInhibitionFactor(final double ethanolConcentration){
		return Math.max(0., 1. - ethanolConcentration / KineticParameters.ETHANOL_INHIBITION_MAX);
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

		private static final double T_MIN_AMYLASE = 5.;
		private static final double T_OPT_AMYLASE = (50. + 55.) / 2.;
		private static final double T_MAX_AMYLASE = 75.;


		private final double stageTemperature;

		BaranyiRobertsODE(final double stageTemperature){
			this.stageTemperature = stageTemperature;
		}

		@Override
		public int getDimension(){
			return KineticParameters.STATE_DIMENSION;
		}

		@Override
		public void computeDerivatives(final double time, final double[] state, final double[] derivatives){
			final double physiologicalState = Math.max(0., state[KineticParameters.IDX_PHYSIOLOGICAL_STATE]);
			final double biomassDensity = Math.max(0., state[KineticParameters.IDX_BIOMASS_DENSITY]);
			final double sugarConcentration = Math.max(0., state[KineticParameters.IDX_SUGAR_CONCENTRATION]);
			final double ethanolConcentration = Math.max(0., state[KineticParameters.IDX_ETHANOL_PRODUCED]);
			final double sugarEnzymatic = Math.max(0., state[KineticParameters.IDX_SUGAR_ENZYMATIC]);

			//── Environmental inhibition factors ──────────────────────────────
			final double totalSugarAvailable = sugarConcentration + sugarEnzymatic;
			final double waterActivity = computeWaterActivity(totalSugarAvailable);
			final double gammaTemperature = cardinalTemperatureFactor(stageTemperature,
				KineticParameters.TEMPERATURE_CARDINAL_MIN, KineticParameters.TEMPERATURE_CARDINAL_OPT,
				KineticParameters.TEMPERATURE_CARDINAL_MAX);
			final double gammaWaterActivity = waterActivityInhibitionFactor(waterActivity);
			final double gammaEthanol = ethanolInhibitionFactor(ethanolConcentration);

			//── Baranyi structural factors ────────────────────────────────────
			//lagAdjustment → 0 during lag phase (Q small), → 1 when fully adapted (Q >> 1)
			final double lagAdjustment = physiologicalState / (physiologicalState + 1.);
			//capacityAdjustment → 0 when biomass saturates the substrate-limited ceiling
			final double dynamicBiomassCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * totalSugarAvailable;
			final double capacityAdjustment = Math.max(0., 1. - biomassDensity / dynamicBiomassCapacity);
			//binary sugar availability flag
			final double sugarAvailability = (totalSugarAvailable > KineticParameters.SUGAR_DEPLETION_THRESHOLD? 1.: 0.);

			//── Effective growth rate [h⁻¹] ──────────────────────────────────
			final double effectiveGrowthRate = KineticParameters.MU_MAX_REF
				* gammaTemperature
				* gammaWaterActivity
				* gammaEthanol
				* lagAdjustment
				* capacityAdjustment
				* sugarAvailability;

			//── ODE derivatives ───────────────────────────────────────────────
			//dQ/dt uses MU_MAX_REF (not effectiveGrowthRate): Q evolves at the maximum rate regardless of environment — it
			// represents intracellular enzyme reserves that replenish at a fixed rate (Baranyi (1994), eq. 2)
			derivatives[KineticParameters.IDX_PHYSIOLOGICAL_STATE] = KineticParameters.MU_MAX_REF * physiologicalState;

			final double biomassGrowthRate = effectiveGrowthRate * biomassDensity;
			derivatives[KineticParameters.IDX_BIOMASS_DENSITY] = biomassGrowthRate;

			final double sugarConsumptionRate = -(1. / KineticParameters.YIELD_BIOMASS_PER_SUGAR) * biomassGrowthRate;
			derivatives[KineticParameters.IDX_SUGAR_CONCENTRATION] = sugarConsumptionRate;
			derivatives[KineticParameters.IDX_CO2_PRODUCED] = KineticParameters.YIELD_CO2_PER_SUGAR * (-sugarConsumptionRate);
			derivatives[KineticParameters.IDX_ETHANOL_PRODUCED] = KineticParameters.YIELD_ETHANOL_PER_SUGAR * (-sugarConsumptionRate);

			// ─ Amylase activity ─────────────────────────────────
			//max 3% from starch
			final double sugarGenerationMax = 0.03;
			//[h⁻¹] — parameter to be calibrated
			final double amylaseActivityRate = 0.15;
			// T_MIN, T_OPT, T_MAX for amylase
			final double gammaTemperatureAmylase = cardinalTemperatureFactor(stageTemperature,
				T_MIN_AMYLASE, T_OPT_AMYLASE, T_MAX_AMYLASE);

			// dS_enz/dt = k · (1 − S_enz/S_max) · γ_T
			final double sugarGenerationRate = amylaseActivityRate * (1. - sugarEnzymatic / sugarGenerationMax)
				* gammaTemperatureAmylase;
			derivatives[KineticParameters.IDX_SUGAR_ENZYMATIC] = sugarGenerationRate;
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
	// @param initialAnhydrousYeastFraction	[g_dry_yeast / g_dough]
	// @return	Final ODE state vector [Q, N, S, P_CO2, E_EtOH].
	// =========================================================================
	private double[] simulate(final double initialAnhydrousYeastFraction){
		final double[] state = new double[KineticParameters.STATE_DIMENSION];
		state[KineticParameters.IDX_PHYSIOLOGICAL_STATE] = initialPhysiologicalState;
		state[KineticParameters.IDX_BIOMASS_DENSITY] = initialAnhydrousYeastFraction;
		state[KineticParameters.IDX_SUGAR_CONCENTRATION] = sugar;
		state[KineticParameters.IDX_CO2_PRODUCED] = 0.;
		state[KineticParameters.IDX_ETHANOL_PRODUCED] = 0.;
		state[KineticParameters.IDX_SUGAR_ENZYMATIC] = 0.;

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
				KineticParameters.ODE_ABSOLUTE_TOLERANCE,
				KineticParameters.ODE_RELATIVE_TOLERANCE);

			integrator.addStepHandler(new NonNegativeClampHandler());
			integrator.integrate(
				new BaranyiRobertsODE(stageTemperatures[stage]),
				stageStartTime, state,
				stageStartTime + stageDuration, state);

			//hard clamp after each stage (belt-and-suspenders)
			for(int i = 0; i < KineticParameters.STATE_DIMENSION; i ++)
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
		final double[] finalState = simulate(initialAnhydrousYeastFraction);
		final double finalTotalSugar = finalState[KineticParameters.IDX_SUGAR_CONCENTRATION]
			+ finalState[KineticParameters.IDX_SUGAR_ENZYMATIC];
		final double dynamicCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * finalTotalSugar;
		return finalState[KineticParameters.IDX_BIOMASS_DENSITY] / dynamicCapacity;
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
	// Q exits the lag state — i.e., Q/(Q+1) → 1 — after time
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
		final double factor = cardinalTemperatureFactor(stageTemperatures[0],
			KineticParameters.TEMPERATURE_CARDINAL_MIN, KineticParameters.TEMPERATURE_CARDINAL_OPT,
			KineticParameters.TEMPERATURE_CARDINAL_MAX);
		return (factor <= 0.
			? Double.POSITIVE_INFINITY
			: Math.log(1. + 1. / initialPhysiologicalState) / KineticParameters.MU_MAX_REF);
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
	public double findOptimalAnhydrousYeastFraction(final boolean outputWarning){
		final String warning = scheduleWarning();
		if(warning != null)
			System.err.println("WARNING: " + warning);

		final double lowerBound = 1e-8;
		final double upperBound = substrateLimitedBiomassCapacity * 0.9999;

		final double progressAtLowerBound = computeFermentationProgress(lowerBound);
		final double progressAtUpperBound = computeFermentationProgress(upperBound);

		//check whether the bracket is valid for BrentSolver
		if(progressAtUpperBound < KineticParameters.TARGET_FERMENTATION_PROGRESS){
			if(outputWarning){
				final double requiredScheduleExtension
					= (KineticParameters.TARGET_FERMENTATION_PROGRESS - progressAtUpperBound)
					* sumArray(stageDurations) / progressAtUpperBound;
				System.err.printf(Locale.US,
					"%n" + "═".repeat(100) + "%n"
						+ "⚠️  SUBSTRATE-LIMITED CEILING REACHED%n"
						+ "═".repeat(100) + "%n%n"
						+ "SUMMARY:%n"
						+ "  Maximum yeast possible: %.2f%% WB%n"
						+ "  Fermentation achieved:  %.0f%%%n"
						+ "  Target fermentation:    %.0f%%%n"
						+ "ROOT CAUSE:%n"
						+ "  The limiting factor is SUBSTRATE (fermentable sugar), not yeast quantity.%n%n"
						+ "TO ACHIEVE TARGET:%n"
						+ "  - OPTION 1 (Recommended): Increase sugar%n"
						+ "  - OPTION 2: Extend fermentation schedule by ~%.1f hours%n"
						+ "RETURNING: Substrate ceiling (%.2f%% WB) as conservative estimate.%n"
						+ "═".repeat(100) + "%n%n",
					toBakersWetPercent(upperBound),
					progressAtUpperBound * 100.,
					KineticParameters.TARGET_FERMENTATION_PROGRESS * 100.,
					requiredScheduleExtension,
					toBakersWetPercent(upperBound));
			}
			return upperBound;
		}

		if(progressAtLowerBound >= KineticParameters.TARGET_FERMENTATION_PROGRESS){
			System.err.printf("WARNING: near-zero yeast already reaches the target. Stages are too long.%n");
			return lowerBound;
		}

		//f(N0) = progress(N0) − TARGET is monotone increasing on [lo, hi]
		final UnivariateFunction objectiveFunction
			= N0 -> computeFermentationProgress(N0) - KineticParameters.TARGET_FERMENTATION_PROGRESS;

		final BrentSolver solver = new BrentSolver(KineticParameters.BRENT_RELATIVE_TOLERANCE,
			KineticParameters.BRENT_ABSOLUTE_TOLERANCE);
		try{
			return solver.solve(KineticParameters.BRENT_MAX_EVALUATIONS, objectiveFunction, lowerBound, upperBound);
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
				stageTemperatures[0],
				KineticParameters.TEMPERATURE_CARDINAL_MIN,
				KineticParameters.TEMPERATURE_CARDINAL_MAX);

		if(lagDoughHours > 0.5 * totalDoughHours)
			return String.format(
				"total wait (lag %.2f h = %.2f h) exceeds 50%% of dough time (%.2f h) "
					+ "for w_b=%.2f — consider a longer schedule or more pre-hydrated yeast.",
				rehydrationDuration * 60., lagDoughHours, totalDoughHours,
				yeastMoisture);

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
		final double flourFraction = 1. - water - salt - sugar;
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
		if(yeastMoisture >= 1.)
			throw new ArithmeticException("yeastMoistureFraction must be < 1.");

		return anhydrousMassFraction / (1. - yeastMoisture);
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

		final double[] state = new double[KineticParameters.STATE_DIMENSION];
		state[KineticParameters.IDX_PHYSIOLOGICAL_STATE] = initialPhysiologicalState;
		state[KineticParameters.IDX_BIOMASS_DENSITY] = initialAnhydrousYeastFraction;
		state[KineticParameters.IDX_SUGAR_CONCENTRATION] = sugar;
		state[KineticParameters.IDX_CO2_PRODUCED] = 0.;
		state[KineticParameters.IDX_ETHANOL_PRODUCED] = 0.;
		state[KineticParameters.IDX_SUGAR_ENZYMATIC] = 0.;

		System.out.printf(Locale.US,
			"%-8s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-8s%n",
			"time[h]", "stage", "Q", "N", "S", "S_enz", "CO2", "EtOH", "aw", "μ_eff");

		double stageStartTime = 0.;
		for(int stage = 0; stage < stageTemperatures.length; stage ++){
			final double stageDuration = stageDurations[stage];
			final double printInterval = stageDuration / PRINT_POINTS_PER_STAGE;

			for(int printPoint = 0; printPoint < PRINT_POINTS_PER_STAGE; printPoint ++){
				final double currentTime = stageStartTime + printPoint * printInterval;
				final double nextTime = Math.min(currentTime + printInterval, stageStartTime + stageDuration);

				//print state at the currentTime before advancing
				final double physiologicalState = state[KineticParameters.IDX_PHYSIOLOGICAL_STATE];
				final double biomassDensity = state[KineticParameters.IDX_BIOMASS_DENSITY];
				final double sugarConcentration = state[KineticParameters.IDX_SUGAR_CONCENTRATION];
				final double sugarEnzymatic = state[KineticParameters.IDX_SUGAR_ENZYMATIC];
				final double totalSugarAvailable = sugarConcentration + sugarEnzymatic;
				final double dynamicBiomassCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * totalSugarAvailable;

				final double waterActivity = computeWaterActivity(totalSugarAvailable);
				final double sugarAvailability = (totalSugarAvailable > KineticParameters.SUGAR_DEPLETION_THRESHOLD? 1.: 0.);
				final double effectiveGrowthRate = KineticParameters.MU_MAX_REF
					* cardinalTemperatureFactor(stageTemperatures[stage],
						KineticParameters.TEMPERATURE_CARDINAL_MIN, KineticParameters.TEMPERATURE_CARDINAL_OPT,
						KineticParameters.TEMPERATURE_CARDINAL_MAX)
					* waterActivityInhibitionFactor(waterActivity)
					* ethanolInhibitionFactor(state[KineticParameters.IDX_ETHANOL_PRODUCED])
					* (physiologicalState / (physiologicalState + 1.))
					* Math.max(0., 1. - biomassDensity / dynamicBiomassCapacity)
					* sugarAvailability;

				System.out.printf(Locale.US,
					"%-8.2f %-6d %-6.2f %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-6.3f %-1.3f%n",
					currentTime,
					stage + 1,
					physiologicalState,
					biomassDensity * 100.,
					sugarConcentration * 100.,
					sugarEnzymatic * 100.,
					state[KineticParameters.IDX_CO2_PRODUCED] * 100.,
					state[KineticParameters.IDX_ETHANOL_PRODUCED] * 100.,
					waterActivity,
					effectiveGrowthRate);

				//advance integrator to the next print point
				final DormandPrince853Integrator integrator = new DormandPrince853Integrator(
					//minimum adaptive step size [h] — prevents stalling
					1e-6,
					//maximum adaptive step size [h] — ODE solver may use fewer steps
					stageDuration,
					KineticParameters.ODE_ABSOLUTE_TOLERANCE,
					KineticParameters.ODE_RELATIVE_TOLERANCE);
				integrator.integrate(
					new BaranyiRobertsODE(stageTemperatures[stage]),
					currentTime, state,
					nextTime, state);
				for(int i = 0; i < KineticParameters.STATE_DIMENSION; i ++)
					state[i] = Math.max(0., state[i]);
			}
			stageStartTime += stageDuration;
		}

		final double finalPhysiologicalState = state[KineticParameters.IDX_PHYSIOLOGICAL_STATE];
		final double finalBiomassDensity = state[KineticParameters.IDX_BIOMASS_DENSITY];
		final double finalSugarConcentration = state[KineticParameters.IDX_SUGAR_CONCENTRATION];
		final double finalSugarEnzymatic = state[KineticParameters.IDX_SUGAR_ENZYMATIC];
		final double finalTotalSugar = finalSugarConcentration + finalSugarEnzymatic;
		final double finalDynamicCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * finalTotalSugar;
		final double finalWaterActivity = computeWaterActivity(finalTotalSugar);
		final double finalSugarAvailability = (finalTotalSugar > KineticParameters.SUGAR_DEPLETION_THRESHOLD? 1.: 0.);
		final double finalEffectiveGrowthRate = KineticParameters.MU_MAX_REF
			* cardinalTemperatureFactor(stageTemperatures[stageTemperatures.length - 1],
				KineticParameters.TEMPERATURE_CARDINAL_MIN, KineticParameters.TEMPERATURE_CARDINAL_OPT,
				KineticParameters.TEMPERATURE_CARDINAL_MAX)
			* waterActivityInhibitionFactor(finalWaterActivity)
			* ethanolInhibitionFactor(state[KineticParameters.IDX_ETHANOL_PRODUCED])
			* (finalPhysiologicalState / (finalPhysiologicalState + 1.))
			* Math.max(0., 1. - finalBiomassDensity / finalDynamicCapacity)
			* finalSugarAvailability;

		System.out.printf(Locale.US,
			"%-8.2f %-6s %-6.2f %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-6.3f %-1.3f%n",
			stageStartTime,
			"FINAL",
			state[KineticParameters.IDX_PHYSIOLOGICAL_STATE],
			state[KineticParameters.IDX_BIOMASS_DENSITY] * 100.,
			state[KineticParameters.IDX_SUGAR_CONCENTRATION] * 100.,
			state[KineticParameters.IDX_SUGAR_ENZYMATIC] * 100.,
			state[KineticParameters.IDX_CO2_PRODUCED] * 100.,
			state[KineticParameters.IDX_ETHANOL_PRODUCED] * 100.,
			finalWaterActivity,
			finalEffectiveGrowthRate);
	}

	public void printSensitivityAnalysis(final double initialAnhydrousYeast){
		System.out.printf(Locale.US,
			"%-7s %-10s %-1s%n",
			"factor", "yeast(WB)", "progress");
		for(final double factor : new double[]{0.25, 0.5, 0.75, 1., 1.5, 2.}){
			final double testFraction = Math.min(
				initialAnhydrousYeast * factor,
				substrateLimitedBiomassCapacity * 0.9999);
			System.out.printf(Locale.US,
				"%-7.2f %-1.2f%%      %-1.2f%n",
				factor,
				toBakersWetPercent(testFraction),
				computeFermentationProgress(testFraction));
		}
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
		final double[] fractions = {1.};
		final double[][] flourMatrix = new double[][]{
			{170., 0.7, 0.72, 0.017, 0.11, 0.007, 0.022, 0.005, 0.}
		};
		final String[] flourTypes = {"wheat"};
		final double flourTemperature = 17.8;
		final double airRelativeHumidity = 0.54;
		final double[] moisture = FlourMoistureModel.compute(fractions, flourMatrix, flourTypes,
			flourTemperature, airRelativeHumidity);
		final double flourFreeWater = moisture[0] - moisture[1];

		//── Dough composition ──────────────────────────────────────────────
		//total = flour + water + salt + sugar
		final double water = 0.65;
		final double salt = 0.012;
		final double sugar = 0.02;

		//── Fermentation schedule ──────────────────────────────────────────
		//two-stage: bulk fermentation + proofing
		final double[] temperatures = {28., 30.};
		final double[] durations = {2., 2.};

		//── Yeast: fresh compressed ────────────────────────────────────────
		//the model computes Q0 from this value and returns an anhydrous mass
		final double yeastMoisture = 0.7;
		// rehydrationDuration = time spent soaking in warm water before mixing
		final double rehydrationDuration = 20. / 60.;

		//── Instantiate model ──────────────────────────────────────────────
		final YeastCalculatorAll model = new YeastCalculatorAll(
			water, salt, sugar,
			temperatures, durations,
			flourFreeWater, yeastMoisture,
			rehydrationDuration);

		//── Optimize ──────────────────────────────────────────────────────
		final double optimalAnhydrousYeast = model.findOptimalAnhydrousYeastFraction(true);

		System.out.printf("%n=== RESULT ===%n");
		System.out.printf(Locale.US,
			"Lag in dough       : %d min%n", (int)Math.ceil(model.estimatedLagTimeInDoughHours() * 60));
		System.out.printf(Locale.US,
			"Optimal yeast (WB) : %.2f%%%n",
			model.toBakersWetPercent(optimalAnhydrousYeast));

		//── Diagnostic trace ──────────────────────────────────────────────
		System.out.printf("%n=== SIMULATION TRACE ===%n");
		model.printSimulationTrace(optimalAnhydrousYeast);

		//── Sensitivity analysis ──────────────────────────────────────────
		System.out.printf("%n=== SENSITIVITY ===%n");
		model.printSensitivityAnalysis(optimalAnhydrousYeast);

		//── Effect of rehydration time on lag (same yeast, same schedule) ─────────
		System.out.printf("%n=== EFFECT OF REHYDRATION TIME ===%n");
		System.out.printf(Locale.US,
			"%-13s %-7s %-9s %-1s%n",
			"t_rehyd[min]", "Q0_eff", "lag[min]", "yeast(WB)");
		for(final double tMin : new double[]{0., 5., 10., 15., 20., 30., 45., 60.}){
			final YeastCalculatorAll mr = new YeastCalculatorAll(
				water, salt, sugar,
				temperatures, durations,
				flourFreeWater, yeastMoisture,
				tMin / 60.);
			final boolean hasWarning = (mr.scheduleWarning() != null);
			final PrintStream originalErr = System.err;
			if(hasWarning)
				System.setErr(new PrintStream(OutputStream.nullOutputStream()));
			final double n0r = mr.findOptimalAnhydrousYeastFraction(false);
			if(hasWarning)
				System.setErr(originalErr);
			System.out.printf(Locale.US,
				"%-13.0f %-7.2f %-9d %-1.2f%%%s%n",
				tMin, mr.initialPhysiologicalState,
				(int)Math.ceil(mr.estimatedLagTimeInDoughHours() * 60),
				mr.toBakersWetPercent(n0r), (hasWarning? " ⚠️": ""));
		}
	}

}
