package io.github.mtrevisan.yeastcalculator.backery;


/**
 * Handles hydration and processing traits of the specific biological leavening agent.
 */
public class YeastInput{

	private final double yeastMoisture;             // [g_water / g_wet_yeast]
	private final double rehydrationDurationHours;  // Rehydration window before mixing


	public YeastInput(double yeastMoisture, double rehydrationDurationHours){
		if(yeastMoisture < 0.0 || yeastMoisture > 1.0){
			throw new IllegalArgumentException("Yeast moisture must be a fraction between 0.0 and 1.0");
		}
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
