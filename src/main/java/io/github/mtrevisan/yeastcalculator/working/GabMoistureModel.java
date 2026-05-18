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

	static GabResult calculateMoisture(final SimulationInputs in, final double[] fractions){
		double mixTotalDB = 0.;
		double mixBoundDB = 0.;

		final double aw = StrictMath.max(StrictMath.min(in.airRelativeHumidity, 0.95), 0.1);
		final double temperatureCoeff = 1. - 0.0025 * (in.flourTemperature - 20.);
		final Object[][] matrix = in.flourMatrix;

		for(int i = 0; i < fractions.length; i++){
			final double protein = ((Number)matrix[i][3]).doubleValue();
			final double fiber = ((Number)matrix[i][6]).doubleValue();
			final String type = (String)matrix[i][7];

			final FlourRegistry.FlourProperties props = FlourRegistry.resolveProperties(type);

			final double wmDb = props.wmBase + 0.085 * protein + 0.12 * fiber;
			final double uEquilibriumDB = (wmDb * props.cGab * props.kGab * aw)
				/ ((1. - props.kGab * aw) * (1. - props.kGab * aw + props.cGab * props.kGab * aw));

			final double uTotalDB = StrictMath.max(StrictMath.min(uEquilibriumDB * temperatureCoeff, 0.2), 0.08);

			mixTotalDB += uTotalDB * fractions[i];
			mixBoundDB += wmDb * fractions[i];
		}

		final double moistureTotal = mixTotalDB / (1. + mixTotalDB);
		final double flourStrictlyBoundWater = mixBoundDB / (1. + mixBoundDB);

		return new GabResult(flourStrictlyBoundWater, moistureTotal - flourStrictlyBoundWater);
	}

}
