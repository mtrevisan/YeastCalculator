package io.github.mtrevisan.yeastcalculator.physiology;

import io.github.mtrevisan.yeastcalculator.kinetics.KineticParameters;


/**
 * Baranyi physiological state (Q) — a measure of yeast metabolic readiness.
 * <p>
 * Q represents the intracellular enzyme reserve pool. During lag phase, Q grows
 * at a constant rate (dQ/dt = μ_max_ref · Q) independent of environment. Once
 * Q reaches ~1, lag phase ends and exponential growth begins.
 * <p>
 * The lag adjustment factor is: α_lag = Q / (Q + 1)
 * - Q → 0: α_lag → 0 (complete lag, no growth)
 * - Q → ∞: α_lag → 1 (fully adapted, exponential growth)
 * <p>
 * Source: Baranyi & Roberts (1994), eq. 2.
 */
public class PhysiologicalState{

	/**
	 * Compute Q0 (initial physiological state) from yeast moisture content.
	 * <p>
	 * Biological rationale:
	 * Trehalose stabilizes dry yeast cells in a glassy state (anabiosis).
	 * Membrane fluidity and enzyme pools are restored as the cell rehydrates.
	 * Below w_b ≈ 0.05 (critical aw ≈ 0.90), the membrane remains rigid
	 * and trehalose mobilization has not begun (Gervais & Beney 2001).
	 * Above that threshold, recovery follows a power-law.
	 * <p>
	 * Calibration points:
	 * w_b = 0.05 (instant dry): Q0 → Q0_ANHYDROUS (cells hydrate in dough)
	 * w_b = 0.70 (fresh compressed): Q0 = Q0_FRESH_COMPRESSED
	 *
	 * @param yeastMoistureFraction	[g_water / g_wet_yeast]
	 * @return	Q0 [dimensionless]
	 */
	public static double fromMoisture(final double yeastMoistureFraction){
		//below activation threshold: cells remain in near-anabiosis
		if(yeastMoistureFraction <= KineticParameters.WB_ACTIVATION_THRESHOLD)
			return KineticParameters.Q0_ANHYDROUS;

		//above threshold: interpolate via power-law (exponent 1.5)
		final double clipped = Math.min(yeastMoistureFraction, KineticParameters.WB_FRESH_COMPRESSED);
		final double normalizedHydration = (clipped - KineticParameters.WB_ACTIVATION_THRESHOLD)
			/ (KineticParameters.WB_FRESH_COMPRESSED - KineticParameters.WB_ACTIVATION_THRESHOLD);

		return KineticParameters.Q0_ANHYDROUS
			+ (KineticParameters.Q0_FRESH_COMPRESSED - KineticParameters.Q0_ANHYDROUS)
			* Math.pow(normalizedHydration, 1.5);
	}

	/**
	 * Compute effective Q after rehydration soak [h].
	 * <p>
	 * During rehydration in warm water, Q grows at maximum rate:
	 * dQ/dt = μ_max_ref · Q  (environment-independent)
	 * Solution: Q(t) = Q0_dry · exp(μ_max_ref · t)
	 * <p>
	 * Even a short soak significantly raises Q0 for very dry yeast.
	 *
	 * @param q0Dry	Initial physiological state [dimensionless]
	 * @param rehydrationDurationHours	Soak time [h]
	 * @return	Effective Q0 at start of dough fermentation [dimensionless]
	 */
	public static double afterRehydration(final double q0Dry, final double rehydrationDurationHours){
		return q0Dry * Math.exp(KineticParameters.MU_MAX_REF * rehydrationDurationHours);
	}

	/**
	 * Lag phase exit time (Baranyi & Roberts 1994, eq. 4).
	 * <p>
	 * λ = ln(1 + 1/Q0) / μ_max_ref
	 * <p>
	 * This is the time for Q to evolve from Q0 to the point where lag adjustment factor α_lag = Q / (Q+1) ≈ 0.95.
	 * It is INDEPENDENT of temperature and aw, because dQ/dt = μ_max_ref · Q is environment-independent.
	 *
	 * @param q0Effective	Physiological state at start of dough [dimensionless].
	 * @return	Lag phase duration [h].
	 */
	public static double lagPhaseHours(final double q0Effective){
		if(q0Effective <= 0.)
			throw new IllegalArgumentException("Q0 must be positive.");

		return Math.log(1. + 1. / q0Effective) / KineticParameters.MU_MAX_REF;
	}

	/**
	 * Lag adjustment factor: α_lag = Q / (Q + 1).
	 * <p>
	 * Used in the Baranyi ODE to modulate growth during lag phase.
	 *
	 * @param q	Physiological state [dimensionless].
	 * @return	Adjustment factor ∈ [0, 1].
	 */
	public static double lagAdjustmentFactor(final double q){
		return q / (q + 1.);
	}

}