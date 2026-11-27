/**
 * Copyright (c) 2026 Mauro Trevisan
 * <p>
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
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
		COMPOSITION.put(FlourType.FOXTAIL_MILLET, new double[]{0.12, 0.07});
		COMPOSITION.put(FlourType.SORGHUM, new double[]{0.10, 0.06});
		COMPOSITION.put(FlourType.BUCKWHEAT, new double[]{0.13, 0.10});
		COMPOSITION.put(FlourType.QUINOA, new double[]{0.14, 0.07});
		COMPOSITION.put(FlourType.AMARANTH, new double[]{0.14, 0.07});
		COMPOSITION.put(FlourType.TEFF, new double[]{0.13, 0.08});
		COMPOSITION.put(FlourType.CHICKPEA, new double[]{0.22, 0.11});
		COMPOSITION.put(FlourType.LENTIL, new double[]{0.25, 0.11});
		COMPOSITION.put(FlourType.PEA, new double[]{0.23, 0.15});
		COMPOSITION.put(FlourType.FAVA_BEAN, new double[]{0.26, 0.13});
		COMPOSITION.put(FlourType.SOY, new double[]{0.36, 0.09});
		COMPOSITION.put(FlourType.CHESTNUT, new double[]{0.06, 0.07});
		COMPOSITION.put(FlourType.WATER_CHESTNUT, new double[]{0.04, 0.035});
		COMPOSITION.put(FlourType.ALMOND, new double[]{0.21, 0.12});
		COMPOSITION.put(FlourType.COCONUT, new double[]{0.18, 0.38});
		COMPOSITION.put(FlourType.HAZELNUT, new double[]{0.15, 0.10});
		COMPOSITION.put(FlourType.FLAXSEED, new double[]{0.20, 0.27});
		COMPOSITION.put(FlourType.HEMP, new double[]{0.30, 0.25});
		COMPOSITION.put(FlourType.SESAME, new double[]{0.25, 0.12});
		COMPOSITION.put(FlourType.WHOLE_BARLEY, new double[]{0.11, 0.12});
		COMPOSITION.put(FlourType.WHOLE_OAT, new double[]{0.13, 0.11});
		COMPOSITION.put(FlourType.MAIZE_FINE, new double[]{0.08, 0.02});
		COMPOSITION.put(FlourType.MAIZE_MEDIUM, new double[]{0.08, 0.03});
		COMPOSITION.put(FlourType.MAIZE_COARSE, new double[]{0.08, 0.04});
		COMPOSITION.put(FlourType.MAIZE_WHITE, new double[]{0.08, 0.03});
		COMPOSITION.put(FlourType.MAIZE_RED, new double[]{0.08, 0.03});
		COMPOSITION.put(FlourType.WHOLE_RICE, new double[]{0.08, 0.035});
		COMPOSITION.put(FlourType.WHOLE_MILLET, new double[]{0.11, 0.08});
		COMPOSITION.put(FlourType.WHOLE_BUCKWHEAT, new double[]{0.13, 0.11});
		COMPOSITION.put(FlourType.WHITE_BEAN, new double[]{0.22, 0.15});
		COMPOSITION.put(FlourType.BLACK_BEAN, new double[]{0.22, 0.16});
		COMPOSITION.put(FlourType.RED_LENTIL, new double[]{0.25, 0.11});
		COMPOSITION.put(FlourType.GREEN_PEA, new double[]{0.23, 0.15});
		COMPOSITION.put(FlourType.LUPIN, new double[]{0.38, 0.10});
		COMPOSITION.put(FlourType.ACORN, new double[]{0.05, 0.06});
		COMPOSITION.put(FlourType.WALNUT, new double[]{0.15, 0.07});
		COMPOSITION.put(FlourType.PISTACHIO, new double[]{0.21, 0.10});
		COMPOSITION.put(FlourType.PEANUT, new double[]{0.26, 0.09});
		COMPOSITION.put(FlourType.BEECH_MAST, new double[]{0.20, 0.10});
		COMPOSITION.put(FlourType.PUMPKIN_SEED, new double[]{0.30, 0.06});
		COMPOSITION.put(FlourType.SUNFLOWER_SEED, new double[]{0.21, 0.09});
		COMPOSITION.put(FlourType.CHIA, new double[]{0.17, 0.34});
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

	/**
	 * Expected equilibrium moisture [% w.b.] at 20°C and water activity 0.30, 0.50, 0.65, as [min, max] pairs.
	 * <p>
	 * The ranges are typical values of the sorption-isotherm literature on flours (type II isotherms; e.g. Moreira et al.,
	 * 2010, for chestnut and wheat flours): starchy flours hold about 7-9 % at a<sub>w</sub> 0.3, 10-12 % at 0.5 and 13-14.5 %
	 * at 0.65, whole-grain/fibrous flours about 1 point more, pulses about 1.5 points less at medium humidity, fat-rich nuts
	 * and seeds much less.
	 * </p>
	 */
	private enum Band{
		STARCHY(6.5, 9.5, 9.5, 12.5, 12., 15.),
		FIBROUS(7., 10., 10., 13., 12., 15.5),
		TEFF(5.5, 8.5, 7.5, 10.5, 9.5, 12.),
		/** Fat-rich kernels/seeds with a published monolayer well below 0.04 g/g d.b. (hazelnut Wm 0.0155, chia 0.024-0.029). */
		LOW_SORPTION(2.5, 4.5, 2.5, 5.5, 3.5, 8.),
		PULSE(6., 8.5, 8.5, 10.5, 11.5, 13.5),
		COCONUT(5.5, 8.5, 8., 10.5, 10.5, 13.5),
		NUT(2.5, 5., 4., 7., 6., 9.),
		SEED(4., 6.5, 6., 9., 8.5, 12.);

		private static final double[] WATER_ACTIVITY = {0.3, 0.5, 0.65};

		private final double[] min = new double[3];
		private final double[] max = new double[3];

		Band(final double min30, final double max30, final double min50, final double max50, final double min65,
				final double max65){
			min[0] = min30;
			max[0] = max30;
			min[1] = min50;
			max[1] = max50;
			min[2] = min65;
			max[2] = max65;
		}
	}

	private static final Map<FlourType, Band> BANDS = new EnumMap<>(FlourType.class);
	static{
		for(final FlourType type : FlourType.values())
			BANDS.put(type, Band.STARCHY);
		for(final FlourType type : new FlourType[]{FlourType.WHOLE_WHEAT, FlourType.WHOLE_RYE, FlourType.WHOLE_SPELT,
				FlourType.WHOLE_EINKORN, FlourType.WHOLE_EMMER, FlourType.WHOLE_KHORASAN, FlourType.OAT, FlourType.FOXTAIL_MILLET,
				FlourType.WHOLE_BARLEY, FlourType.WHOLE_OAT, FlourType.WHOLE_RICE, FlourType.WHOLE_MILLET, FlourType.WHOLE_BUCKWHEAT,
				FlourType.BUCKWHEAT, FlourType.QUINOA, FlourType.AMARANTH, FlourType.TEFF})
			BANDS.put(type, Band.FIBROUS);
		for(final FlourType type : new FlourType[]{FlourType.CHICKPEA, FlourType.LENTIL, FlourType.PEA, FlourType.FAVA_BEAN,
				FlourType.SOY, FlourType.WHITE_BEAN, FlourType.BLACK_BEAN, FlourType.RED_LENTIL, FlourType.GREEN_PEA,
				FlourType.LUPIN})
			BANDS.put(type, Band.PULSE);
		BANDS.put(FlourType.COCONUT, Band.COCONUT);
		BANDS.put(FlourType.TEFF, Band.TEFF);
		for(final FlourType type : new FlourType[]{FlourType.ALMOND, FlourType.HAZELNUT, FlourType.SESAME,
				FlourType.WALNUT, FlourType.PISTACHIO, FlourType.PEANUT, FlourType.BEECH_MAST})
			BANDS.put(type, Band.NUT);
		BANDS.put(FlourType.FLAXSEED, Band.SEED);
		BANDS.put(FlourType.HEMP, Band.SEED);
		BANDS.put(FlourType.PUMPKIN_SEED, Band.NUT);	// about 50 % fat: sorbs like a nut
		BANDS.put(FlourType.SUNFLOWER_SEED, Band.NUT);	// about 50 % fat: sorbs like a nut
		BANDS.put(FlourType.HAZELNUT, Band.LOW_SORPTION);
		BANDS.put(FlourType.CHIA, Band.LOW_SORPTION);

	}


	@Test
	void everyFlourHasALiteratureBand(){
		for(final FlourType type : FlourType.values())
			assertTrue(BANDS.containsKey(type), "missing literature band for " + type);
	}

	@Test
	void moistureMatchesLiteratureAtLowMediumAndHighHumidity(){
		for(final FlourType type : FlourType.values()){
			final Band band = BANDS.get(type);
			for(int i = 0; i < Band.WATER_ACTIVITY.length; i ++){
				final double moisture = total(single(type, 20., Band.WATER_ACTIVITY[i])) * 100.;

				assertTrue(moisture >= band.min[i] && moisture <= band.max[i],
					String.format("%s at a_w %.2f: %.2f %% w.b. outside literature range [%.1f, %.1f]", type,
						Band.WATER_ACTIVITY[i], moisture, band.min[i], band.max[i]));
			}
		}
	}

	@Test
	void chestnutMonolayerMatchesChenlo(){
		// monolayer 0.059 ± 0.007 g/g d.b. (GAB desorption, chestnut flour, Chenlo et al.; Moreira et al. 2010 print 0.062 for desorption and 0.049 for adsorption); the model monolayer is back-calculated from the wet-basis results
		final GabMoistureModel.GabResult r = single(FlourType.CHESTNUT, 20., 0.5);
		final double monolayerDB = r.flourStrictlyBoundWater / (1. - total(r));

		assertEquals(0.059, monolayerDB, 0.007);
	}

	@Test
	void legumesAreDrierThanCerealsAtMediumHumidity(){
		final double wheat = total(single(FlourType.WHEAT, 20., 0.65));
		for(final FlourType type : new FlourType[]{FlourType.CHICKPEA, FlourType.LENTIL, FlourType.PEA, FlourType.FAVA_BEAN,
				FlourType.SOY})
			assertTrue(total(single(type, 20., 0.65)) < wheat, type + " should be drier than wheat at a_w 0.65");
	}

	@Test
	void newFloursAreInTheRightNeighbourhood(){
		final double millet = total(single(FlourType.MILLET, 20., 0.5));
		final double foxtail = total(single(FlourType.FOXTAIL_MILLET, 20., 0.5));
		final double pea = total(single(FlourType.PEA, 20., 0.5));
		final double fava = total(single(FlourType.FAVA_BEAN, 20., 0.5));
		final double wheat = total(single(FlourType.WHEAT, 20., 0.5));
		final double waterChestnut = total(single(FlourType.WATER_CHESTNUT, 20., 0.5));

		assertEquals(millet, foxtail, 0.015);
		assertEquals(pea, fava, 0.01);
		assertEquals(wheat, waterChestnut, 0.015);
	}

	@Test
	void wholeGrainFloursHoldAtLeastAsMuchMoistureAsTheirRefinedCounterparts(){
		final FlourType[][] pairs = {{FlourType.BARLEY, FlourType.WHOLE_BARLEY}, {FlourType.OAT, FlourType.WHOLE_OAT},
			{FlourType.RICE, FlourType.WHOLE_RICE}, {FlourType.MILLET, FlourType.WHOLE_MILLET},
			{FlourType.BUCKWHEAT, FlourType.WHOLE_BUCKWHEAT}};
		for(final FlourType[] pair : pairs)
			assertTrue(total(single(pair[1], 20., 0.65)) >= total(single(pair[0], 20., 0.65)) - 0.005,
				pair[1] + " should not be drier than " + pair[0]);
	}

	@Test
	void maizeVariantsAreCloseToEachOther(){
		final double reference = total(single(FlourType.MAIZE_MEDIUM, 20., 0.5));
		for(final FlourType type : new FlourType[]{FlourType.MAIZE_FINE, FlourType.MAIZE_COARSE, FlourType.MAIZE_WHITE,
				FlourType.MAIZE_RED})
			assertEquals(reference, total(single(type, 20., 0.5)), 0.01, type + " too far from medium maize");
	}

	@Test
	void newPulsesAreDrierThanCerealsAndNutsAreDrierThanPulses(){
		final double wheat = total(single(FlourType.WHEAT, 20., 0.65));
		final double lentil = total(single(FlourType.LENTIL, 20., 0.65));
		for(final FlourType type : new FlourType[]{FlourType.WHITE_BEAN, FlourType.BLACK_BEAN, FlourType.GREEN_PEA,
				FlourType.LUPIN})
			assertTrue(total(single(type, 20., 0.65)) < wheat, type + " should be drier than wheat");
		assertEquals(lentil, total(single(FlourType.RED_LENTIL, 20., 0.65)), 0.005);
		for(final FlourType type : new FlourType[]{FlourType.WALNUT, FlourType.PISTACHIO, FlourType.PEANUT,
				FlourType.BEECH_MAST, FlourType.PUMPKIN_SEED, FlourType.SUNFLOWER_SEED, FlourType.CHIA})
			assertTrue(total(single(type, 20., 0.65)) < lentil, type + " should be drier than lentil");
	}

	@Test
	void acornIsStarchyLikeChestnut(){
		assertEquals(total(single(FlourType.CHESTNUT, 20., 0.5)), total(single(FlourType.ACORN, 20., 0.5)), 0.015);
	}

	@Test
	void flourMoistureAt50PercentHumidityIsInTheUsualRange(){
		for(final FlourType type : FlourType.values()){
			final double moisture = total(single(type, 20., 0.5));

			assertTrue(moisture > 0.025 && moisture < 0.14, type + " moisture " + moisture + " out of 2.5-14% at 20°C/50%");
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

			// 0.01 tolerance: published whole wheat (Martin-Santos 2012) and wheat flour (Moreira 2010) isotherms come from different studies and differ by less than 1 point
			assertTrue(total(whole) > total(refined) - 0.01, pair[0] + " should not be much drier than " + pair[1]);
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
