package io.github.mtrevisan.yeastcalculator.backery;


/**
 * Standard implementation of the Guggenheim-Anderson-de Boer (GAB) adsorption isotherm.
 */
public class GabMoistureModel{

	public static final class GabResult{
		public final double flourStrictlyBoundWater;
		public final double flourActiveWater;

		private GabResult(final double bound, final double active){
			this.flourStrictlyBoundWater = bound;
			this.flourActiveWater = active;
		}
	}


	public static GabResult calculateMoisture(final SimulationInputs in){
		double[] fractionsRaw = in.getFractions();
		FlourInput[] flourMatrix = in.getFlourMatrix();
		int numFlours = flourMatrix.length;

		double totalRawAmount = 0.0;
		for(double fraction : fractionsRaw){
			totalRawAmount += fraction;
		}

		double[] fractions = new double[numFlours];
		for(int i = 0; i < numFlours; i++){
			fractions[i] = totalRawAmount > 0? fractionsRaw[i] / totalRawAmount: 0.0;
		}

		double[] wmDB = new double[numFlours];
		double[] uTotalDB = new double[numFlours];
		double[] flourMoistureWB = new double[numFlours];
		double[] dryMassComponent = new double[numFlours];
		double totalDryMassComponent = 0.0;

		double aw = clamp(in.getAirRelativeHumidity(), 0.1, 0.95);
		double temperatureCoeff = 1.0 - 0.0025 * (in.getFlourTemperature() - 20.0);

		for(int i = 0; i < numFlours; i++){
			FlourInput flour = flourMatrix[i];
			FlourType type = flour.getType();

			wmDB[i] = type.getWmBase() + (0.085 * flour.getProtein()) + (0.12 * flour.getFiber());
			double tmp = type.getkGab() * aw;
			double uEquilibriumDB = (wmDB[i] * type.getcGab() * tmp) / ((1.0 - tmp) * (1.0 + (type.getcGab() - 1.0) * tmp));
			uTotalDB[i] = clamp(uEquilibriumDB * temperatureCoeff, 0.08, 0.2);
			flourMoistureWB[i] = uTotalDB[i] / (1.0 + uTotalDB[i]);
			dryMassComponent[i] = fractions[i] * (1.0 - flourMoistureWB[i]);
			totalDryMassComponent += dryMassComponent[i];
		}

		double mixTotalDB = 0.0;
		double mixBoundDB = 0.0;

		for(int i = 0; i < numFlours; i++){
			double dryFraction = totalDryMassComponent > 0? dryMassComponent[i] / totalDryMassComponent: 0.0;
			mixTotalDB += uTotalDB[i] * dryFraction;
			mixBoundDB += wmDB[i] * dryFraction;
		}

		double boundWater = mixBoundDB / (1.0 + mixTotalDB);
		double activeWater = (mixTotalDB - mixBoundDB) / (1.0 + mixTotalDB);

		return new GabResult(boundWater, activeWater);
	}

	private static double clamp(double x, double lo, double hi){
		return x < lo? lo: (x > hi? hi: x);
	}

}
