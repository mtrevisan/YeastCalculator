package io.github.mtrevisan.yeastcalculator.model;

import java.util.Arrays;


/** Fermentation schedule: temperatures and durations. */
public class FermentationSchedule{

	/** Stage temperatures [°C]. */
	public final double[] temperatures;
	/** Stage durations [h]. */
	public final double[] durations;


	public FermentationSchedule(final double[] temperatures, final double[] durations){
		if(temperatures.length != durations.length)
			throw new IllegalArgumentException("Temperature and duration arrays must have equal length.");

		this.temperatures = temperatures.clone();
		this.durations = durations.clone();
	}

	public int getStageCount(){
		return temperatures.length;
	}

	public double getTotalDuration(){
		return Arrays.stream(durations).sum();
	}

	public double getStageTemperature(final int stage){
		return temperatures[stage];
	}

	public double getStageDuration(final int stage){
		return durations[stage];
	}

}