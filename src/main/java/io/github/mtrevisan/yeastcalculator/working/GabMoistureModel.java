package io.github.mtrevisan.yeastcalculator.working;

public class GabMoistureModel{

	public static class GabResult{
		public double flourStrictlyBoundWater;
		public double flourActiveWater;

	}

	public static GabResult calculateMoisture(final SimulationInputs in, final double[] fractions){
		double mixTotalDB = 0.0;
		double mixBoundDB = 0.0;

		final double aw = StrictMath.max(StrictMath.min(in.airRelativeHumidity, 0.95), 0.1);
		final double temperatureCoeff = 1.0 - 0.0025 * (in.flourTemperature - 20.0);

		for(int i = 0; i < fractions.length; i++){
			final double protein = ((Number)in.flourMatrix[i][3]).doubleValue();
			final double fiber = ((Number)in.flourMatrix[i][6]).doubleValue();
			final String type = (String)in.flourMatrix[i][7];

			// Fetch properties directly from the registry
			final FlourRegistry.FlourProperties props = FlourRegistry.resolveProperties(type);

			final double wmDb = props.wmBase + 0.085 * protein + 0.12 * fiber;
			final double uEquilibriumDB = (wmDb * props.cGab * props.kGab * aw)
				/ ((1.0 - props.kGab * aw) * (1.0 - props.kGab * aw + props.cGab * props.kGab * aw));

			final double uTotalDB = StrictMath.max(StrictMath.min(uEquilibriumDB * temperatureCoeff, 0.2), 0.08);

			mixTotalDB += uTotalDB * fractions[i];
			mixBoundDB += wmDb * fractions[i];
		}

		final double moistureTotal = mixTotalDB / (1.0 + mixTotalDB);
		final double flourStrictlyBoundWater = mixBoundDB / (1.0 + mixBoundDB);

		final GabResult res = new GabResult();
		res.flourStrictlyBoundWater = flourStrictlyBoundWater;
		res.flourActiveWater = moistureTotal - flourStrictlyBoundWater;
		return res;
	}

}
