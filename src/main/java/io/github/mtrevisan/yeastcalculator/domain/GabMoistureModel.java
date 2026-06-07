package io.github.mtrevisan.yeastcalculator.domain;

import io.github.mtrevisan.yeastcalculator.optimization.SimulationInputs;


public final class GabMoistureModel{

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
		final double[] fractions = in.getFractions();
		double mixTotalDB = 0.;
		double mixBoundDB = 0.;

		final double aw = Math.clamp(in.getAirRelativeHumidity(), 0.1, 0.95);
		final double temperatureCoeff = 1. - 0.0025 * (in.getFlourTemperature() - 20.);
		final FlourInput[] matrix = in.getFlourMatrix();

		final int length = fractions.length;
		final double[] uTotalDBArray = new double[length];
		final double[] wmDbArray = new double[length];
		final double[] dryMassComponent = new double[length];
		double sumDryMass = 0.;
		for(int i = 0; i < length; i ++){
			final FlourInput flour = matrix[i];
			final double wmDb = flour.getWmBase() + 0.085 * flour.getProtein() + 0.12 * flour.getFiber();
			final double uEquilibriumDB = (wmDb * flour.getCGab() * flour.getKGab() * aw)
				/ ((1. - flour.getKGab() * aw) * (1. - flour.getKGab() * aw + flour.getCGab() * flour.getKGab() * aw));
			final double uTotalDB = Math.clamp(uEquilibriumDB * temperatureCoeff, 0.08, 0.2);

			wmDbArray[i] = wmDb;
			uTotalDBArray[i] = uTotalDB;

			// Temporary conversion to wet base moisture to extract the dry component (1 - moistureWB)
			final double flourMoistureWB = uTotalDB / (1. + uTotalDB);
			dryMassComponent[i] = fractions[i] * (1. - flourMoistureWB);
			sumDryMass += dryMassComponent[i];
		}

		// calculation of weighted averages based on DRY FRACTIONS
		for(int i = 0; i < length; i ++){
			final double dryFraction = dryMassComponent[i] / sumDryMass;

			mixTotalDB += uTotalDBArray[i] * dryFraction;
			mixBoundDB += wmDbArray[i] * dryFraction;
		}

		final double flourStrictlyBoundWater = mixBoundDB / (1. + mixTotalDB);
		final double flourActiveWater = (mixTotalDB - mixBoundDB) / (1. + mixTotalDB);
		return new GabResult(flourStrictlyBoundWater, flourActiveWater);
	}

}
