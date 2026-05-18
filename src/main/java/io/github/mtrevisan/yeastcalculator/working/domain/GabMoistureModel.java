package io.github.mtrevisan.yeastcalculator.working.domain;

import io.github.mtrevisan.yeastcalculator.working.optimization.SimulationInputs;


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

		for(int i = 0; i < fractions.length; i ++){
			final FlourInput flour = matrix[i];
			final double wmDb = flour.getWmBase() + 0.085 * flour.getProtein() + 0.12 * flour.getFiber();
			final double uEquilibriumDB = (wmDb * flour.getCGab() * flour.getKGab() * aw)
				/ ((1. - flour.getKGab() * aw) * (1. - flour.getKGab() * aw + flour.getCGab() * flour.getKGab() * aw));
			final double uTotalDB = Math.clamp(uEquilibriumDB * temperatureCoeff, 0.08, 0.2);

			mixTotalDB += uTotalDB * fractions[i];
			mixBoundDB += wmDb * fractions[i];
		}

		final double moistureTotal = mixTotalDB / (1. + mixTotalDB);
		final double flourStrictlyBoundWater = mixBoundDB / (1. + mixBoundDB);

		return new GabResult(flourStrictlyBoundWater, moistureTotal - flourStrictlyBoundWater);
	}

}
