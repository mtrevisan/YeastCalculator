package io.github.mtrevisan.yeastcalculator.backery;


import java.util.Arrays;


/**
 * Orchestrator payload containing the complete input matrix of the bake environment.
 */
public class SimulationInputs{

	private final double[] fractions;
	private final FlourInput[] flourMatrix;
	private final double flourTemperature;
	private final double airRelativeHumidity;
	private final YeastInput yeastProperties;
	private final DoughRecipe recipe;
	private final KneadingInput kneading;
	private final StageInput[] stages;
	private final double[] folds;


	public SimulationInputs(double[] fractions, FlourInput[] flourMatrix,
			double flourTemperature, double airRelativeHumidity,
			YeastInput yeastProperties, DoughRecipe recipe,
			KneadingInput kneading, StageInput[] stages, double[] folds){
		if(fractions == null || flourMatrix == null || fractions.length != flourMatrix.length){
			throw new IllegalArgumentException("Fractions and flourMatrix arrays must match in dimensions.");
		}

		this.fractions = Arrays.copyOf(fractions, fractions.length);
		this.flourMatrix = Arrays.copyOf(flourMatrix, flourMatrix.length);
		this.flourTemperature = flourTemperature;
		this.airRelativeHumidity = airRelativeHumidity;
		this.yeastProperties = yeastProperties;
		this.recipe = recipe;
		this.kneading = kneading;
		this.stages = stages != null? Arrays.copyOf(stages, stages.length): new StageInput[0];
		this.folds = folds != null? Arrays.copyOf(folds, folds.length): new double[0];
	}

	public double[] getFractions(){
		return Arrays.copyOf(fractions, fractions.length);
	}

	public FlourInput[] getFlourMatrix(){
		return Arrays.copyOf(flourMatrix, flourMatrix.length);
	}

	public double getFlourTemperature(){
		return flourTemperature;
	}

	public double getAirRelativeHumidity(){
		return airRelativeHumidity;
	}

	public YeastInput getYeastProperties(){
		return yeastProperties;
	}

	public DoughRecipe getRecipe(){
		return recipe;
	}

	public KneadingInput getKneading(){
		return kneading;
	}

	public StageInput[] getStages(){
		return Arrays.copyOf(stages, stages.length);
	}

	public double[] getFolds(){
		return Arrays.copyOf(folds, folds.length);
	}

}
