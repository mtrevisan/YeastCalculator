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
package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.FlourType;
import io.github.mtrevisan.yeastcalculator.domain.KneadingInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.domain.YeastInput;
import io.github.mtrevisan.yeastcalculator.optimization.FoldScheduleOptimizer;
import io.github.mtrevisan.yeastcalculator.optimization.YeastOptimizer;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;

import java.util.Locale;


/**
 * Inverse of {@link Main1}: water, salt, oil, malt and the rest of the process (flour, kneading, stages, folds) are
 * fixed by the user, and the amount of yeast and the stretch & fold schedule are computed.
 */
public class Main2{

	public static void main(final String[] args){
		// Enforce dot formatting for clean metric output readings
		Locale.setDefault(Locale.US);

		System.out.println("====================================================");
		System.out.println("            YEAST AMOUNT CALCULATOR                 ");
		System.out.println("====================================================\n");

		// 1. Flour blend
		final double[] fractions = {1.};
		final FlourInput[] flourMatrix = {
			new FlourInput(205., 0.523, 0.014, 0.12, 0.01, 0.017, 0.003, FlourType.WHEAT)
		};

		// 2. Environmental baselines
		final double flourTemperature = 19.6;
		final double airRelativeHumidity = 0.54;

		// 3. Yeast properties
		final double yeastMoisture = 0.70;
		final YeastInput yeastProps = new YeastInput(yeastMoisture);

		// 4. Mechanical Kneading profiles and autolyze (flour and water resting together before adding the yeast and
		// the salt, 0 for none: it develops the gluten with no work, so fewer folds are needed)
		final double autolyseMinutes = 30.;
		final KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.MANUAL, 15., autolyseMinutes);

		// 5. Fermentation schedule
		final StageInput[] stages = {
			new StageInput(28., 0.55, 5.29)
		};

		// 6. Stretch & Fold timestamps [min]: an output of the model, so none at the start
		final double[] noFolds = {};

		// 7. Recipe INPUT (ratios over the flour weight): water, salt, oil, malt
		final double waterRatio = 0.62;
		final double saltRatio = 0.022;
		final double oilRatio = 0.02;
		final double maltRatio = 0.002;
		final double maltSugarContent = 0.1;
		final DoughRecipe recipe = new DoughRecipe(waterRatio, saltRatio, oilRatio, maltRatio, maltSugarContent,
			(15_000. / 110.) * (1. - maltSugarContent), 0.05);

		final SimulationInputs baseInputs = new SimulationInputs(fractions, flourMatrix, flourTemperature,
			airRelativeHumidity, yeastProps, recipe, kneading, stages, noFolds);

		final BakeryProduct selectedProduct = BakeryProduct.GASTRONOMY_PAN_PIZZA;

		System.out.println("Selected Target Product: " + selectedProduct.name());
		System.out.printf("Water                : %.1f%%\n", waterRatio * 100.);
		System.out.printf("Salt                 : %.1f%%\n", saltRatio * 100.);
		System.out.printf("Oil                  : %.1f%%\n", oilRatio * 100.);
		System.out.printf("Malt                 : %.2f%%\n", maltRatio * 100.);
		System.out.printf("Autolyse             : %.0f min\n", autolyseMinutes);
		System.out.println("Computing folds and yeast...\n");

		// --- Fold schedule and yeast inversion ---
		final FoldScheduleOptimizer.OptimizedFoldResult result = FoldScheduleOptimizer.optimizeFolds(baseInputs,
			selectedProduct);
		final double optimalYeastRatio = result.optimalYeastRatio();
		final double[] folds = result.folds();
		final SimulationInputs inputs = new SimulationInputs(fractions, flourMatrix, flourTemperature,
			airRelativeHumidity, yeastProps, recipe, kneading, stages, folds);

		System.out.println("\n====================================================");
		System.out.printf("Optimal Yeast        : %.2f%% of flour weight\n", optimalYeastRatio * 100.);
		System.out.printf("Optimal S&F Count    : %d sets\n", folds.length);
		final StringBuilder schedule = new StringBuilder();
		for(final double fold : folds)
			schedule.append(schedule.isEmpty()? "": ", ").append(String.format("%.0f", fold));
		System.out.printf("S&F Timestamps [min] : %s\n", folds.length > 0? schedule: "none");
		System.out.println("====================================================");

		// Verification run at the computed yeast amount
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(inputs);
		final double[] finalState = YeastOptimizer.runSimulation(inputs, gab, optimalYeastRatio, selectedProduct);
		if(finalState != null){
			System.out.println("\n--- FINAL DOUGH STATE AT TIMELINE EXPIRATION ---");
			System.out.printf("Yeast (Dry Biomass)  : %.2f%%\n", finalState[0] * 100.);
			System.out.printf("Sugars Remaining     : %.1f%%\n", finalState[1] * 100.);
			System.out.printf("Retained Gas Volume  : %.1f ml/g_flour\n", finalState[2]);
			System.out.printf("Ethanol Accumulation : %.1f g/kg_dough\n", finalState[3] * 1000.);

			DoughDiagnostics.report(inputs, finalState, selectedProduct);
		}
		else
			System.err.println("Critical Error: verification simulation failed.");
	}

}
