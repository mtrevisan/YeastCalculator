package io.github.mtrevisan.yeastcalculator.working;


public class SimulationInputs{
	private double[] fractionsRaw = {0.70, 0.20, 0.10};

	private Object[][] flourMatrix = {
		// W,		P/L,	Sugar,	Prot,	Fat,		Fiber,	Ashes,	Type
		{300.,	0.55,	0.015,	0.13,	0.012,	0.02,		0.0055,	"wheat"},
		{220.,	0.60,	0.010,	0.12,	0.015,	0.025,	0.0065,	"wheat semolina fine"},
		{110.,	0.40,	0.030,	0.10,	0.020,	0.06,		0.015,	"rye"}
	};

	private double[][] stagesRaw = {
		// Temp,	RH,	Duration
		{24.,		0.75,	2.},
		{28.,		0.80,	3.5},
		{4.,		0.70,	12.}
	};

	private double[] foldsRaw = {1.2, 2.5};

	private double doughWater = 0.65;
	private double doughSalt = 0.02;
	private double doughOil = 0.03;
	private double yeastMoisture = 0.70;

	private double flourTemperature = 22.;
	private double airRelativeHumidity = 0.55;


	int getFlourCount(){
		return fractionsRaw.length;
	}

	double[] getFractions(){
		final int flours = fractionsRaw.length;
		double sumFractions = 0.;
		for(final double f : fractionsRaw)
			sumFractions += f;
		final double[] fractions = new double[flours];
		for(int i = 0; i < flours; i ++)
			fractions[i] = fractionsRaw[i] / sumFractions;
		return fractions;
	}

	public double getDoughWater(){
		return doughWater;
	}

	public Object[][] getFlourMatrix(){
		return flourMatrix;
	}

	public double getDoughSalt(){
		return doughSalt;
	}

	public double getDoughOil(){
		return doughOil;
	}

	public double[][] getStagesRaw(){
		return stagesRaw;
	}

	public double[] getFoldsRaw(){
		return foldsRaw;
	}

	public double getYeastMoisture(){
		return yeastMoisture;
	}

	public double getAirRelativeHumidity(){
		return airRelativeHumidity;
	}

	public double getFlourTemperature(){
		return flourTemperature;
	}

}
