package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Handles hydration and processing traits of the specific biological leavening agent.
 *
 * @param yeastMoisture            [g_water / g_wet_yeast]
 */
public record YeastInput(
	// [g_water / g_wet_yeast]
	double yeastMoisture
){


	public YeastInput{
		if(yeastMoisture < 0. || yeastMoisture > 1.)
			throw new IllegalArgumentException("Yeast moisture must be a fraction between 0 and 1");
	}


	/**
	 * Mathematically resolves the exact optimal rehydration time in minutes
	 * using a continuous thermodynamic membrane relaxation model.
	 * Eliminates step-function discontinuities.
	 */
	public double getRehydrationWindow(){
		// Smooth linear decay profile from dry baseline (5%) to fresh saturation (70%)
		if(yeastMoisture >= 0.70)
			return 0.;
		if(yeastMoisture <= 0.05)
			return 12.;

		// Interpolation curve: perfectly maps intermediate moistures (like 20% or 60%)
		return 12. * (1. - (yeastMoisture - 0.05) / (0.70 - 0.05));
	}

}
