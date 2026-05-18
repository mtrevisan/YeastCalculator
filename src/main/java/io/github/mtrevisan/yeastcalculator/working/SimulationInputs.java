package io.github.mtrevisan.yeastcalculator.working;


public class SimulationInputs{

//	private double[] fractions = {0.70, 0.20, 0.10};
	private double[] fractions = {1.};

	private final FlourInput[] flourMatrix = {
//		new FlourInput(300., 0.55, 0.015, 0.13, 0.012, 0.02, 0.0055, FlourType.WHEAT),
//		new FlourInput(220., 0.60, 0.010, 0.12, 0.015, 0.025, 0.0065, FlourType.WHEAT_SEMOLINA_FINE),
//		new FlourInput(110., 0.40, 0.030, 0.10, 0.020, 0.06, 0.015, FlourType.RYE)
		new FlourInput(295., 0.55, 0.013, 0.13, 0.011, 0.019, 0.003, FlourType.WHEAT)
	};

	private StageInput[] stages = {
//		new StageInput(24., 0.75, 2.),
//		new StageInput(28., 0.80, 3.5},
//		new StageInput(4., 0.70, 12.)
		new StageInput(26., 0.55, 2.)
	};

//	private double[] folds = {1.2, 2.5};
	private double[] folds = {};

	private double doughWater = 0.615;
	private double doughSalt = 0.022;
	private double doughOil = 0.039;
	private double yeastMoisture = 0.70;

	private double flourTemperature = 19.;
	private double airRelativeHumidity = 0.54;


	int getFlourCount(){
		return fractions.length;
	}

	double[] getFractions(){
		final int flours = fractions.length;
		double sumFractions = 0.;
		for(final double f : fractions)
			sumFractions += f;
		final double[] fractions = new double[flours];
		for(int i = 0; i < flours; i++)
			fractions[i] = this.fractions[i] / sumFractions;
		return fractions;
	}

	FlourInput[] getFlourMatrix(){
		return flourMatrix;
	}

	double getDoughWater(){
		return doughWater;
	}

	double getDoughSalt(){
		return doughSalt;
	}

	double getDoughOil(){
		return doughOil;
	}

	StageInput[] getStages(){
		return stages;
	}

	double[] getFolds(){
		return folds;
	}

	double getYeastMoisture(){
		return yeastMoisture;
	}

	double getAirRelativeHumidity(){
		return airRelativeHumidity;
	}

	double getFlourTemperature(){
		return flourTemperature;
	}

}
