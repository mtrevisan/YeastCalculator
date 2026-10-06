package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Handles hydration and processing traits of the specific biological leavening agent.
 *
 * @param yeastMoisture            [g_water / g_wet_yeast]
 * @param rehydrationDurationHours Rehydration window before mixing
 */
public record YeastInput(
	// [g_water / g_wet_yeast]
	double yeastMoisture,
	// Rehydration window before mixing
	double rehydrationDurationHours){


	public YeastInput{
		if(yeastMoisture < 0. || yeastMoisture > 1.)
			throw new IllegalArgumentException("Yeast moisture must be a fraction between 0 and 1");
	}


}
