package io.github.mtrevisan.yeastcalculator.working;


final class GabMoistureModel{

	static final class GabResult{
		final double flourStrictlyBoundWater;
		final double flourActiveWater;

		private GabResult(final double bound, final double active){
			this.flourStrictlyBoundWater = bound;
			this.flourActiveWater = active;
		}
	}


	private GabMoistureModel(){}


	static GabResult calculateMoisture(final SimulationInputs in){
		final double[] fractions = in.getFractions();

		double mixTotalDB = 0.;
		double mixBoundDB = 0.;

		final double aw = clamp(in.getAirRelativeHumidity(), 0.1, 0.95);
		final double temperatureCoeff = 1. - 0.0025 * (in.getFlourTemperature() - 20.);
		final FlourInput[] matrix = in.getFlourMatrix();

		for(int i = 0; i < fractions.length; i ++){
			final FlourInput flour = matrix[i];
			final double protein = flour.getProtein();
			final double fiber = flour.getFiber();
			final String type = flour.getType();

			final FlourRegistry.FlourProperties props = FlourRegistry.resolveProperties(type);

			final double wmDb = props.wmBase + 0.085 * protein + 0.12 * fiber;
			final double uEquilibriumDB = (wmDb * props.cGab * props.kGab * aw)
				/ ((1. - props.kGab * aw) * (1. - props.kGab * aw + props.cGab * props.kGab * aw));
			final double uTotalDB = clamp(uEquilibriumDB * temperatureCoeff, 0.08, 0.2);

			mixTotalDB += uTotalDB * fractions[i];
			mixBoundDB += wmDb * fractions[i];
		}

		final double moistureTotal = mixTotalDB / (1. + mixTotalDB);
		final double flourStrictlyBoundWater = mixBoundDB / (1. + mixBoundDB);

		return new GabResult(flourStrictlyBoundWater, moistureTotal - flourStrictlyBoundWater);
	}

	private static double clamp(final double x, final double min, final double max){
		return StrictMath.max(StrictMath.min(x, min), max);
	}

}
