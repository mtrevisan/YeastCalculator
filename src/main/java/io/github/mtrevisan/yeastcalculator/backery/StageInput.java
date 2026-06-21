package io.github.mtrevisan.yeastcalculator.backery;


/**
 * Identifies ambient boundaries for a single chronological resting or leavening stage.
 */
public class StageInput{

	private final double temperature;
	private final double relativeHumidity;
	private final double durationHours;


	public StageInput(double temperature, double relativeHumidity, double durationHours){
		this.temperature = temperature;
		this.relativeHumidity = relativeHumidity;
		this.durationHours = durationHours;
	}


	public double getTemperature(){
		return temperature;
	}

	public double getRelativeHumidity(){
		return relativeHumidity;
	}

	public double getDurationHours(){
		return durationHours;
	}

}
