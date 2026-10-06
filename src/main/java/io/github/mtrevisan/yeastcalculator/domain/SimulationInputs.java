package io.github.mtrevisan.yeastcalculator.domain;

import java.util.Arrays;


/**
 * Orchestrator payload containing the complete input matrix of the bake environment.
 */
public record SimulationInputs(
	// Flour Matrix inputs
	double[] fractions,
	FlourInput[] flourMatrix,
	// Environmental Baselines
	double flourTemperature,
	double airRelativeHumidity,
	// Yeast Target Specs
	YeastInput yeastProperties,
	// Process Timelines
	DoughRecipe recipe,
	KneadingInput kneading,
	StageInput[] stages,
	double[] folds){


	public SimulationInputs(final double[] fractions, final FlourInput[] flourMatrix,
			final double flourTemperature, final double airRelativeHumidity,
			final YeastInput yeastProperties, final DoughRecipe recipe,
			final KneadingInput kneading, final StageInput[] stages, final double[] folds){
		if(fractions == null || flourMatrix == null || fractions.length != flourMatrix.length)
			throw new IllegalArgumentException("Fractions and flourMatrix arrays must match in dimensions.");

		this.fractions = Arrays.copyOf(fractions, fractions.length);
		this.flourMatrix = Arrays.copyOf(flourMatrix, flourMatrix.length);
		this.flourTemperature = flourTemperature;
		this.airRelativeHumidity = airRelativeHumidity;
		this.yeastProperties = yeastProperties;
		this.recipe = recipe;
		this.kneading = kneading;
		this.stages = (stages != null? Arrays.copyOf(stages, stages.length): new StageInput[0]);
		this.folds = (folds != null? Arrays.copyOf(folds, folds.length): new double[0]);
	}

}
