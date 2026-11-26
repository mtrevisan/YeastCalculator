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
package io.github.mtrevisan.yeastcalculator.domain;

import java.util.Arrays;


/**
 * Orchestrator payload containing the complete input matrix of the bake environment.
 */
public record SimulationInputs(
	// Flour Matrix inputs
	double[] fractions,
	FlourInput[] flourMatrix,
	// Environmental Baselines
	double flourTemperature,
	double airRelativeHumidity,
	// Yeast Target Specs
	YeastInput yeastProperties,
	// Process Timelines
	DoughRecipe recipe,
	KneadingInput kneading,
	StageInput[] stages,
	double[] folds){


	public SimulationInputs(final double[] fractions, final FlourInput[] flourMatrix,
			final double flourTemperature, final double airRelativeHumidity,
			final YeastInput yeastProperties, final DoughRecipe recipe,
			final KneadingInput kneading, final StageInput[] stages, final double[] folds){
		if(fractions == null || flourMatrix == null || fractions.length != flourMatrix.length)
			throw new IllegalArgumentException("Fractions and flourMatrix arrays must match in dimensions.");

		this.fractions = Arrays.copyOf(fractions, fractions.length);
		this.flourMatrix = Arrays.copyOf(flourMatrix, flourMatrix.length);
		this.flourTemperature = flourTemperature;
		this.airRelativeHumidity = airRelativeHumidity;
		this.yeastProperties = yeastProperties;
		this.recipe = recipe;
		this.kneading = kneading;
		this.stages = (stages != null? Arrays.copyOf(stages, stages.length): new StageInput[0]);
		this.folds = (folds != null? Arrays.copyOf(folds, folds.length): new double[0]);
	}

}
