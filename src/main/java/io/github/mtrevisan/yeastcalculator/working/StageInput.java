package io.github.mtrevisan.yeastcalculator.working;


public final class StageInput{

	private final double temperature;
	private final double relativeHumidity;
	private final double durationHours;


	public StageInput(final double temperature, final double relativeHumidity, final double durationHours){
		this.temperature = temperature;
		// Restrict RH between 0.0 and 1.0 just in case
		this.relativeHumidity = Math.clamp(relativeHumidity, 0., 1.);
		this.durationHours = Math.max(0., durationHours);
	}


	double getTemperature(){
		return temperature;
	}

	double getRelativeHumidity(){
		return relativeHumidity;
	}

	double getDuration(){
		return durationHours;
	}

}
