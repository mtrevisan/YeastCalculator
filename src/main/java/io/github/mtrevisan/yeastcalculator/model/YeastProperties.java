package io.github.mtrevisan.yeastcalculator.model;

import io.github.mtrevisan.yeastcalculator.physiology.PhysiologicalState;


/**
 * Yeast moisture content, Q0 calculation, and rehydration effects.
 */
public class YeastProperties{

	/**
	 * Moisture content of the yeast [g_water / g_wet_yeast].
	 * Used both to derive Q0 and to convert the anhydrous result to wet-yeast mass.
	 */
	public final double moisture;
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
	public final double rehydrationDuration;
	/**
	 * Effective initial physiological state at the start of dough fermentation.
	 * Accounts for any pre-activation in water: Q0_eff = Q0_dry × exp(μ_max × t_rehydration).
	 * When rehydrationDurationHours is equal to 0, this equals the raw Q0 from moisture content.
	 */
	public final double initialPhysiologicalState;


	public YeastProperties(final double moisture, final double rehydrationDuration){
		this.moisture = moisture;
		this.rehydrationDuration = rehydrationDuration;

		// Q0 of the dry yeast, then advanced by the rehydration soak.
		// dQ/dt = MU_MAX_REF · Q has the exact solution Q(t) = Q0 · exp(MU_MAX_REF · t).
		final double q0Dry = PhysiologicalState.fromMoisture(moisture);
		this.initialPhysiologicalState = PhysiologicalState.afterRehydration(q0Dry, rehydrationDuration);
	}

	/**
	 * Converts an anhydrous mass fraction [g/g_dough] to the equivalent wet-yeast fraction [g_wet_yeast / g_dough]:
	 * <pre>
	 *   wetFraction = anhydrousFraction / (1 − yeastMoistureFraction)
	 * </pre>
	 *
	 * @param anhydrousMassFraction	[g_dry_yeast / g_dough]
	 * @return	[g_wet_yeast / g_dough]
	 */
	public double toWetYeastFraction(final double anhydrousMassFraction){
		if(moisture >= 1.)
			throw new ArithmeticException("yeastMoisture must be < 1.");

		return anhydrousMassFraction / (1. - moisture);
	}

}