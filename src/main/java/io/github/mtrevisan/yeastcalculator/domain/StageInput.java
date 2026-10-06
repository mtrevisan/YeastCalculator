package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Represents an immutable environmental timeline configuration block (a single proofing stage).
 * <p>
 * This class acts as a value object holding thermodynamic boundary conditions
 * (temperature, relative humidity, and duration) required to drive the chemical
 * and biological kinetic equations during a specific segment of the dough's life cycle.
 * </p>
 */
public record StageInput(
	// The ambient temperature in Celsius [°C].
	double temperature,
	// The ambient relative humidity as a decimal fraction [0, 1].
	double relativeHumidity,
	// The total duration of this specific stage in hours (>= 0).
	double duration){

}
