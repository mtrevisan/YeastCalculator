package io.github.mtrevisan.yeastcalculator.thermodynamic;

import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.FlourType;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import org.junit.jupiter.api.Test;

import java.util.EnumMap;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class GabMoistureModelTest{

	private static final double EPS = 1.e-12;

	/** Typical [protein, fiber] mass fractions of each flour (indicative values). */
	private static final Map<FlourType, double[]> COMPOSITION = new EnumMap<>(FlourType.class);
	static{
		COMPOSITION.put(FlourType.WHEAT, new double[]{0.12, 0.017});
		COMPOSITION.put(FlourType.WHEAT_SEMOLINA, new double[]{0.12, 0.03});
		COMPOSITION.put(FlourType.WHEAT_SEMOLINA_FINE, new double[]{0.12, 0.03});
		COMPOSITION.put(FlourType.DURUM_WHEAT, new double[]{0.13, 0.03});
		COMPOSITION.put(FlourType.DURUM_WHEAT_SEMOLINA, new double[]{0.13, 0.03});
		COMPOSITION.put(FlourType.DURUM_WHEAT_SEMOLINA_FINE, new double[]{0.13, 0.03});
		COMPOSITION.put(FlourType.KHORASAN, new double[]{0.14, 0.02});
		COMPOSITION.put(FlourType.SPELT, new double[]{0.14, 0.02});
		COMPOSITION.put(FlourType.EMMER, new double[]{0.13, 0.02});
		COMPOSITION.put(FlourType.EINKORN, new double[]{0.15, 0.02});
		COMPOSITION.put(FlourType.TRITICALE, new double[]{0.12, 0.02});
		COMPOSITION.put(FlourType.RYE, new double[]{0.08, 0.03});
		COMPOSITION.put(FlourType.WHOLE_WHEAT, new double[]{0.13, 0.10});
		COMPOSITION.put(FlourType.WHOLE_RYE, new double[]{0.10, 0.13});
		COMPOSITION.put(FlourType.WHOLE_SPELT, new double[]{0.15, 0.09});
		COMPOSITION.put(FlourType.WHOLE_EINKORN, new double[]{0.16, 0.09});
		COMPOSITION.put(FlourType.WHOLE_EMMER, new double[]{0.14, 0.09});
		COMPOSITION.put(FlourType.WHOLE_KHORASAN, new double[]{0.15, 0.09});
		COMPOSITION.put(FlourType.BARLEY, new double[]{0.10, 0.04});
		COMPOSITION.put(FlourType.OAT, new double[]{0.13, 0.10});
		COMPOSITION.put(FlourType.MAIZE, new double[]{0.08, 0.03});
		COMPOSITION.put(FlourType.RICE, new double[]{0.07, 0.01});
		COMPOSITION.put(FlourType.MILLET, new double[]{0.11, 0.04});
		COMPOSITION.put(FlourType.SORGHUM, new double[]{0.10, 0.06});
		COMPOSITION.put(FlourType.BUCKWHEAT, new double[]{0.13, 0.10});
		COMPOSITION.put(FlourType.QUINOA, new double[]{0.14, 0.07});
		COMPOSITION.put(FlourType.AMARANTH, new double[]{0.14, 0.07});
		COMPOSITION.put(FlourType.TEFF, new double[]{0.13, 0.08});
		COMPOSITION.put(FlourType.CHICKPEA, new double[]{0.22, 0.11});
		COMPOSITION.put(FlourType.LENTIL, new double[]{0.25, 0.11});
		COMPOSITION.put(FlourType.PEA, new double[]{0.23, 0.15});
		COMPOSITION.put(FlourType.SOY, new double[]{0.36, 0.09});
		COMPOSITION.put(FlourType.CHESTNUT, new double[]{0.06, 0.07});
		COMPOSITION.put(FlourType.ALMOND, new double[]{0.21, 0.12});
		COMPOSITION.put(FlourType.COCONUT, new double[]{0.18, 0.38});
		COMPOSITION.put(FlourType.HAZELNUT, new double[]{0.15, 0.10});
		COMPOSITION.put(FlourType.FLAXSEED, new double[]{0.20, 0.27});
		COMPOSITION.put(FlourType.HEMP, new double[]{0.30, 0.25});
		COMPOSITION.put(FlourType.SESAME, new double[]{0.25, 0.12});
	}


	private static FlourInput flour(final FlourType type){
		final double[] c = COMPOSITION.get(type);
		return new FlourInput(200., 0.5, 0.014, c[0], 0.01, c[1], 0.005, type);
	}

	private static SimulationInputs inputs(final double temperature, final double humidity, final double[] fractions,
			final FlourType... types){
		final FlourInput[] matrix = new FlourInput[types.length];
		for(int i = 0; i < types.length; i ++)
			matrix[i] = flour(types[i]);
		return new SimulationInputs(fractions, matrix, temperature, humidity, null, null, null, null, null);
	}

	/** @return total moisture, wet basis [g/g]. */
	private static double total(final GabMoistureModel.GabResult r){
		return r.flourStrictlyBoundWater + r.flourActiveWater;
	}

	private static GabMoistureModel.GabResult single(final FlourType type, final double temperature, final double humidity){
		return GabMoistureModel.calculateMoisture(inputs(temperature, humidity, new double[]{1.}, type));
	}


	@Test
	void compositionCoversEveryFlourType(){
		for(final FlourType type : FlourType.values())
			assertTrue(COMPOSITION.containsKey(type), "missing test composition for " + type);
	}

	@Test
	void everyFlourIsPhysicallyPlausible(){
		final double[][] points = {{20., 0.5}, {25., 0.75}, {5., 0.3}};
		for(final FlourType type : FlourType.values())
			for(final double[] p : points){
				final GabMoistureModel.GabResult r = single(type, p[0], p[1]);
				final String msg = type + " @ " + p[0] + "°C/" + p[1];

				assertTrue(r.flourStrictlyBoundWater > 0., "bound water must be positive: " + msg);
				assertTrue(r.flourActiveWater >= 0., "active water must not be negative: " + msg);
				assertTrue(total(r) >= 0.03 / 1.03 - EPS, "moisture below lower clamp: " + msg);
				assertTrue(total(r) < 0.20, "moisture above 20% w.b.: " + msg);
			}
	}

	@Test
	void flourMoistureAt50PercentHumidityIsInTheUsualRange(){
		for(final FlourType type : FlourType.values()){
			final double moisture = total(single(type, 20., 0.5));

			assertTrue(moisture > 0.04 && moisture < 0.14, type + " moisture " + moisture + " out of 4-14% at 20°C/50%");
		}
	}

	@Test
	void moistureDoesNotDecreaseWithHumidity(){
		for(final FlourType type : FlourType.values()){
			double previous = 0.;
			for(double rh = 0.2; rh <= 0.9; rh += 0.1){
				final double moisture = total(single(type, 20., rh));

				assertTrue(moisture >= previous - EPS, type + " moisture decreased at RH " + rh);
				previous = moisture;
			}
		}
	}

	@Test
	void hotterFlourHoldsLessMoistureWhenNotClamped(){
		final double cool = total(single(FlourType.WHEAT, 10., 0.5));
		final double warm = total(single(FlourType.WHEAT, 20., 0.5));
		final double hot = total(single(FlourType.WHEAT, 35., 0.5));

		assertTrue(cool > warm);
		assertTrue(warm > hot);
	}

	@Test
	void humidityIsClampedBetween10And95Percent(){
		assertEquals(total(single(FlourType.WHEAT, 20., 0.95)), total(single(FlourType.WHEAT, 20., 1.5)), EPS);
		assertEquals(total(single(FlourType.WHEAT, 20., 0.1)), total(single(FlourType.WHEAT, 20., 0.)), EPS);
	}

	@Test
	void wholeGrainFlourHoldsMoreWaterThanRefinedCounterpart(){
		final FlourType[][] pairs = {
			{FlourType.WHOLE_WHEAT, FlourType.WHEAT},
			{FlourType.WHOLE_RYE, FlourType.RYE},
			{FlourType.WHOLE_SPELT, FlourType.SPELT},
			{FlourType.WHOLE_EINKORN, FlourType.EINKORN},
			{FlourType.WHOLE_EMMER, FlourType.EMMER},
			{FlourType.WHOLE_KHORASAN, FlourType.KHORASAN}};
		for(final FlourType[] pair : pairs){
			final GabMoistureModel.GabResult whole = single(pair[0], 20., 0.5);
			final GabMoistureModel.GabResult refined = single(pair[1], 20., 0.5);

			assertTrue(total(whole) > total(refined), pair[0] + " should be wetter than " + pair[1]);
			assertTrue(whole.flourStrictlyBoundWater > refined.flourStrictlyBoundWater,
				pair[0] + " should bind more water than " + pair[1]);
		}
	}

	@Test
	void fatRichNutsAndSeedsAreDrierThanWheat(){
		final double wheat = total(single(FlourType.WHEAT, 20., 0.5));
		final FlourType[] dry = {FlourType.HAZELNUT, FlourType.ALMOND, FlourType.SESAME, FlourType.FLAXSEED, FlourType.HEMP};
		for(final FlourType type : dry)
			assertTrue(total(single(type, 20., 0.5)) < wheat, type + " should be drier than wheat");
	}

	@Test
	void flaxseedAndHempStayBelowNineAndAHalfPercentAtMidHumidity(){
		assertTrue(total(single(FlourType.FLAXSEED, 20., 0.5)) < 0.095);
		assertTrue(total(single(FlourType.HEMP, 20., 0.5)) < 0.095);
	}

	@Test
	void mixMoistureLiesBetweenComponents(){
		final double a = total(single(FlourType.WHEAT, 20., 0.5));
		final double b = total(single(FlourType.HAZELNUT, 20., 0.5));

		final double mix = total(GabMoistureModel.calculateMoisture(
			inputs(20., 0.5, new double[]{0.7, 0.3}, FlourType.WHEAT, FlourType.HAZELNUT)));

		assertTrue(mix > Math.min(a, b));
		assertTrue(mix < Math.max(a, b));
	}

	@Test
	void fractionsAreNormalized(){
		final GabMoistureModel.GabResult one = GabMoistureModel.calculateMoisture(
			inputs(20., 0.5, new double[]{0.3, 0.7}, FlourType.WHEAT, FlourType.WHOLE_RYE));
		final GabMoistureModel.GabResult scaled = GabMoistureModel.calculateMoisture(
			inputs(20., 0.5, new double[]{30., 70.}, FlourType.WHEAT, FlourType.WHOLE_RYE));

		assertEquals(one.flourStrictlyBoundWater, scaled.flourStrictlyBoundWater, EPS);
		assertEquals(one.flourActiveWater, scaled.flourActiveWater, EPS);
	}

	@Test
	void identicalFloursMixToTheSingleFlour(){
		final GabMoistureModel.GabResult single = single(FlourType.HEMP, 20., 0.5);
		final GabMoistureModel.GabResult mix = GabMoistureModel.calculateMoisture(
			inputs(20., 0.5, new double[]{0.4, 0.6}, FlourType.HEMP, FlourType.HEMP));

		assertEquals(single.flourStrictlyBoundWater, mix.flourStrictlyBoundWater, EPS);
		assertEquals(single.flourActiveWater, mix.flourActiveWater, EPS);
	}

}
