package io.github.mtrevisan.yeastcalculator.working.domain;


/**
 * Represents an immutable environmental timeline configuration block (a single proofing stage).
 * <p>
 * This class acts as a value object holding thermodynamic boundary conditions
 * (temperature, relative humidity, and duration) required to drive the chemical
 * and biological kinetic equations during a specific segment of the dough's life cycle.
 * </p>
 */
public final class StageInput{

	private final double temperature;
	private final double relativeHumidity;
	private final double durationHours;


	/**
	 * Constructs a new proofing stage with strictly validated environmental parameters.
	 *
	 * @param temperature	The ambient temperature in Celsius [°C].
	 * @param relativeHumidity	The ambient relative humidity as a decimal fraction [0.0, 1.0].
	 * @param durationHours	The total duration of this specific stage in hours (>= 0.0).
	 */
	public StageInput(final double temperature, final double relativeHumidity, final double durationHours){
		this.temperature = temperature;
		this.relativeHumidity = relativeHumidity;
		this.durationHours = Math.max(0., durationHours);
	}


	/**
	 * Gets the targeted stage temperature.
	 *
	 * @return	The temperature in Celsius [°C].
	 */
	public double getTemperature(){
		return temperature;
	}

	/**
	 * Gets the safely clamped relative humidity.
	 *
	 * @return	The relative humidity as a decimal fraction.
	 */
	public double getRelativeHumidity(){
		return relativeHumidity;
	}

	/**
	 * Gets the safely bounded duration of the stage.
	 *
	 * @return	The duration in hours.
	 */
	public double getDuration(){
		return durationHours;
	}

}
