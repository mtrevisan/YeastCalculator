package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Handles hydration and processing traits of the specific biological leavening agent.
 */
public class YeastInput{

	// [g_water / g_wet_yeast]
	private final double yeastMoisture;
	// Rehydration window before mixing
	private final double rehydrationDurationHours;


	public YeastInput(final double yeastMoisture, final double rehydrationDurationHours){
		if(yeastMoisture < 0. || yeastMoisture > 1.)
			throw new IllegalArgumentException("Yeast moisture must be a fraction between 0.0 and 1.0");

		this.yeastMoisture = yeastMoisture;
		this.rehydrationDurationHours = rehydrationDurationHours;
	}


	public double getYeastMoisture(){
		return yeastMoisture;
	}

	public double getRehydrationDurationHours(){
		return rehydrationDurationHours;
	}

}
