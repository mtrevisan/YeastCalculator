package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.FlourType;
import io.github.mtrevisan.yeastcalculator.domain.KneadingInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.domain.YeastInput;
import io.github.mtrevisan.yeastcalculator.optimization.YeastOptimizer;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;

import java.util.Locale;


/**
 * Inverse of {@link Main1}: water, salt, oil, malt and the rest of the process (flour, kneading, stages, folds) are
 * fixed by the user, and only the amount of yeast is computed.
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

		// 4. Mechanical Kneading profiles and autolyse (flour and water resting together before adding the yeast and
		// the salt, 0 for none: it develops the gluten with no work, so fewer folds are needed)
		final double autolyseMinutes = 30.;
		final KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.MANUAL, 15., autolyseMinutes);

		// 5. Fermentation schedule
		final StageInput[] stages = {
			new StageInput(28., 0.55, 5.29)
		};

		// 6. Stretch & Fold timestamps [min], fixed by the user
		final double[] folds = {30., 60., 90.};

		// 7. Recipe INPUT (ratios over the flour weight): water, salt, oil, malt
		final double waterRatio = 0.62;
		final double saltRatio = 0.022;
		final double oilRatio = 0.02;
		final double maltRatio = 0.002;
		final double maltSugarContent = 0.1;
		final DoughRecipe recipe = new DoughRecipe(waterRatio, saltRatio, oilRatio, maltRatio, maltSugarContent,
			(15_000. / 110.) * (1. - maltSugarContent), 0.05);

		final SimulationInputs inputs = new SimulationInputs(fractions, flourMatrix, flourTemperature,
			airRelativeHumidity, yeastProps, recipe, kneading, stages, folds);

		final BakeryProduct selectedProduct = BakeryProduct.GASTRONOMY_PAN_PIZZA;

		System.out.println("Selected Target Product: " + selectedProduct.name());
		System.out.printf("Water                : %.1f%%\n", waterRatio * 100.);
		System.out.printf("Salt                 : %.1f%%\n", saltRatio * 100.);
		System.out.printf("Oil                  : %.1f%%\n", oilRatio * 100.);
		System.out.printf("Malt                 : %.2f%%\n", maltRatio * 100.);
		System.out.printf("Autolyse             : %.0f min\n", autolyseMinutes);
		System.out.printf("S&F Count            : %d sets\n", folds.length);
		System.out.println("Computing yeast...\n");

		// --- Yeast inversion ---
		final double optimalYeastRatio = YeastOptimizer.findOptimalYeast(inputs, selectedProduct, null);

		System.out.println("====================================================");
		System.out.printf("Optimal Yeast        : %.2f%% of flour weight\n", optimalYeastRatio * 100.);
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
