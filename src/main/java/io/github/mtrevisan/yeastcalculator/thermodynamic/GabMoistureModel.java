package io.github.mtrevisan.yeastcalculator.thermodynamic;

import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.FlourType;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;


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


	private GabMoistureModel(){}


	public static GabResult calculateMoisture(final SimulationInputs in){
		final double[] fractionsRaw = in.getFractions();
		final FlourInput[] flourMatrix = in.getFlourMatrix();
		final int numFlours = flourMatrix.length;

		double totalRawAmount = 0.;
		for(final double fraction : fractionsRaw)
			totalRawAmount += fraction;

		final double[] fractions = new double[numFlours];
		for(int i = 0; i < numFlours; i ++)
			fractions[i] = (totalRawAmount > 0? fractionsRaw[i] / totalRawAmount: 0.);

		final double[] wmDB = new double[numFlours];
		final double[] uTotalDB = new double[numFlours];
		final double[] flourMoistureWB = new double[numFlours];
		final double[] dryMassComponent = new double[numFlours];
		double totalDryMassComponent = 0.;

		final double aw = clamp(in.getAirRelativeHumidity(), 0.1, 0.95);
		final double temperatureCoeff = 1. - 0.0025 * (in.getFlourTemperature() - 20.);

		for(int i = 0; i < numFlours; i ++){
			final FlourInput flour = flourMatrix[i];
			final FlourType type = flour.getType();

			wmDB[i] = type.getWmBase() + (0.085 * flour.getProtein()) + (0.12 * flour.getFiber());
			final double tmp = type.getkGab() * aw;
			final double uEquilibriumDB = (wmDB[i] * type.getcGab() * tmp)
				/ ((1. - tmp) * (1. + (type.getcGab() - 1.) * tmp));
			uTotalDB[i] = clamp(uEquilibriumDB * temperatureCoeff, 0.08, 0.2);
			flourMoistureWB[i] = uTotalDB[i] / (1. + uTotalDB[i]);
			dryMassComponent[i] = fractions[i] * (1. - flourMoistureWB[i]);
			totalDryMassComponent += dryMassComponent[i];
		}

		double mixTotalDB = 0.;
		double mixBoundDB = 0.;
		for(int i = 0; i < numFlours; i ++){
			final double dryFraction = (totalDryMassComponent > 0.? dryMassComponent[i] / totalDryMassComponent: 0.);
			mixTotalDB += uTotalDB[i] * dryFraction;
			mixBoundDB += wmDB[i] * dryFraction;
		}

		final double boundWater = mixBoundDB / (1. + mixTotalDB);
		final double activeWater = (mixTotalDB - mixBoundDB) / (1. + mixTotalDB);
		return new GabResult(boundWater, activeWater);
	}

	private static double clamp(final double x, final double lo, final double hi){
		return (x < lo? lo: (x > hi? hi: x));
	}

}
