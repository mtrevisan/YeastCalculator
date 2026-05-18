package io.github.mtrevisan.yeastcalculator.working;


public class SimulationInputs{

	private double[] fractions = {0.70, 0.20, 0.10};

	private Object[][] flourMatrix = {
		// W,		P/L,	Sugar,	Prot,	Fat,		Fiber,	Ashes,	Type
		{300.,	0.55,	0.015,	0.13,	0.012,	0.02,		0.0055,	"wheat"},
		{220.,	0.60,	0.010,	0.12,	0.015,	0.025,	0.0065,	"wheat semolina fine"},
		{110.,	0.40,	0.030,	0.10,	0.020,	0.06,		0.015,	"rye"}
	};

	private double[][] stages = {
		// Temp,	RH,	Duration
		{24.,		0.75,	2.},
		{28.,		0.80,	3.5},
		{4.,		0.70,	12.}
	};

	private double[] folds = {1.2, 2.5};

	private double doughWater = 0.65;
	private double doughSalt = 0.02;
	private double doughOil = 0.03;
	private double yeastMoisture = 0.70;

	private double flourTemperature = 22.;
	private double airRelativeHumidity = 0.55;


	int getFlourCount(){
		return fractions.length;
	}

	double[] getFractions(){
		final int flours = fractions.length;
		double sumFractions = 0.;
		for(final double f : fractions)
			sumFractions += f;
		final double[] fractions = new double[flours];
		for(int i = 0; i < flours; i ++)
			fractions[i] = this.fractions[i] / sumFractions;
		return fractions;
	}

	double getDoughWater(){
		return doughWater;
	}

	Object[][] getFlourMatrix(){
		return flourMatrix;
	}

	double getDoughSalt(){
		return doughSalt;
	}

	double getDoughOil(){
		return doughOil;
	}

	double[][] getStages(){
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
