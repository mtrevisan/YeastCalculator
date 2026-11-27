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
		final double[] fractionsRaw = in.fractions();
		final FlourInput[] flourMatrix = in.flourMatrix();
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

		final double aw = Math.clamp(in.airRelativeHumidity(), 0.1, 0.95);
		final double temperatureCoeff = 1. - 0.0025 * (in.flourTemperature() - 20.);

		for(int i = 0; i < numFlours; i ++){
			final FlourInput flour = flourMatrix[i];
			final FlourType type = flour.type();

			// fat-rich materials (hazelnut, chia) bind far less water per unit of protein and fiber than cereals do
			wmDB[i] = type.getWmBase() + type.getCompositionFactor() * ((0.085 * flour.protein()) + (0.12 * flour.fiber()));
			final double tmp = type.getkGab() * aw;
			final double uEquilibriumDB = (wmDB[i] * type.getcGab() * tmp)
				/ ((1. - tmp) * (1. + (type.getcGab() - 1.) * tmp));
			// The total moisture can never be below the monolayer: protein- and fiber-rich flours (soy, pulses) have a
			// monolayer above the usual clamp bounds, which would otherwise give a negative active water
			uTotalDB[i] = Math.clamp(uEquilibriumDB * temperatureCoeff, Math.max(0.03, wmDB[i]), Math.max(0.2, wmDB[i]));
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


}
