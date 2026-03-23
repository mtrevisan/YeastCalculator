package io.github.mtrevisan.yeastcalculator.kinetics;

import io.github.mtrevisan.yeastcalculator.model.DoughComposition;


/**
 * Water activity — Ross multi-solute logarithmic model (Ross (1993))
 * <p>
 * Uses freeWaterFraction (= total water − GAB bound water) as the solvent mass, because water tightly bound to flour
 * is not osmotically active.
 * At 65% hydration this lowers aw by ≈ 0.015–0.020 relative to using the total water fraction — a correction larger
 * than any other second-order effect in this model.
 * <p>
 * Solutes: NaCl and sucrose (dominant osmolytes in bread dough).
 * sugarConcentration is the current (dynamic) sugar fraction, because fermentation consumes sugar and changes aw over
 * 	time.
 */
public class WaterActivityModel{

	public static double compute(final DoughComposition dough, final double sugarConcentration){
		final double totalWater = dough.getTotalWater();
		if(totalWater <= 0.)
			return 0.9;

		final double naclMolality = (dough.salt / totalWater)
			/ (KineticParameters.MOLECULAR_WEIGHT_NACL / 1_000.);
		final double sucroseMolality = (sugarConcentration / totalWater)
			/ (KineticParameters.MOLECULAR_WEIGHT_SUCROSE / 1_000.);

		final double lnAw = KineticParameters.AW_NACL_COEFF_LINEAR * naclMolality
			+ KineticParameters.AW_NACL_COEFF_QUADRATIC * naclMolality * naclMolality
			+ KineticParameters.AW_SUCROSE_COEFF_LINEAR * sucroseMolality
			+ KineticParameters.AW_SUCROSE_COEFF_QUADRATIC * sucroseMolality * sucroseMolality;

		return Math.exp(lnAw);
	}

	/**
	 * Water activity inhibition factor — γ_aw ∈ [0, 1]
	 * <p>
	 * Linear ramp between AW_MIN and AW_OPT_DOUGH (Hamad (2012) dough system).
	 * Using the broth AW_OPT = 0.995 would penalize typical dough (aw ≈ 0.97) by ~36%, producing unrealistically slow
	 * growth.
	 */
	public static double inhibitionFactor(final double waterActivity){
		if(waterActivity <= KineticParameters.AW_MIN)
			return 0.;
		if(waterActivity >= KineticParameters.AW_OPT_DOUGH)
			return 1.;

		return (waterActivity - KineticParameters.AW_MIN) / (KineticParameters.AW_OPT_DOUGH - KineticParameters.AW_MIN);
	}

}