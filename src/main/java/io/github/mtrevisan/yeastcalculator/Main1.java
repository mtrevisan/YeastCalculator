package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.FlourType;
import io.github.mtrevisan.yeastcalculator.domain.KneadingInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.domain.YeastInput;
import io.github.mtrevisan.yeastcalculator.optimization.FullRecipeMultivariateObjective;
import io.github.mtrevisan.yeastcalculator.optimization.GlobalMultivariateRecipeOptimizer;
import io.github.mtrevisan.yeastcalculator.optimization.YeastOptimizer;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;

import java.util.Locale;


/* Ingredient optimization */
/*
 * 1. What is a Poolish? (The Liquid Pre-ferment)
 * A Poolish is a liquid pre-ferment mixed at a 1:1 ratio of flour to water (100% hydration) with a tiny fraction of yeast (usually 0.1% to 0.3% of the Poolish flour).
 * • Biochemical Reality: Because it is completely liquid, the yeast cells can swim freely. Mass transport is highly efficient, leading to a massive, rapid production of organic acids (primarily lactic acid) and ethanol. Protease enzymes break down proteins intensely, yielding a high concentration of free amino acids.
 * • Impact on the Final Dough: A Poolish drastically enhances the extensibility (L) of the dough, making it incredibly stretchy, soft, and aromatic. It is the gold standard for high-hydration products with large, irregular air pockets, like your ROMAN_PAN_PIZZA.
 *
 * 2. What is a Biga? (The Solid Pre-ferment)
 * A Biga is a dry, crumbly pre-ferment mixed at a 44% to 50% hydration with roughly 1% fresh yeast relative to the Biga flour. It is not kneaded; it is just turned over until no dry flour remains.
 * • Biochemical Reality: Because water is highly restricted (locked tightly by starch and fibers in the GAB isotherm), the yeast cells are structurally trapped and experience high osmotic stress. Under these conditions, the fermentation leans heavily toward acetic acid production rather than lactic acid.
 * • Impact on the Final Dough: Acetic acid strengthens the gluten network, drastically increasing the dough's tenacity (P) and structural elasticity. When mixed into the final dough, a Biga provides a massive "oven spring" (the explosive rise in the first minutes of baking) and creates a crispy, highly resilient crumb. It is ideal for structural loaves like BREAD.
 */
public class Main1{

	public static void main(final String[] args){
		// Enforce dot formatting for clean metric output readings
		Locale.setDefault(Locale.US);

		System.out.println("====================================================");
		System.out.println("   BOBYQA MULTIVARIATE RECIPE OPTIMIZATION ENGINE   ");
		System.out.println("====================================================\n");

		// 1. Definition of the Flour Blend
		final double[] fractions = {1.};
		final FlourInput[] flourMatrix = {
			new FlourInput(205., 0.523, 0.014, 0.12, 0.01, 0.017, 0.003, FlourType.WHEAT)
		};

		// 2. Structural & Storage environmental baselines
		final double flourTemperature = 19.6;
		final double airRelativeHumidity = 0.54;

		// 3. Yeast properties
		final double targetYeastMoisture = 0.70;
		final YeastInput yeastProps = new YeastInput(targetYeastMoisture);

		// 4. Mechanical Kneading profiles and autolyse (flour and water resting together before adding the yeast and
		// the salt, 0 for none: it develops the gluten with no work, so fewer folds are needed)
		final double autolyseMinutes = 30.;
		final KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.MANUAL, 15., autolyseMinutes);

		// 5. Multi-stage fermentation schedule (Hard constraints set by user)
		final StageInput[] stages = {
			new StageInput(28., 0.55, 5.29)
		};

		// 6. Physical structural interventions (Stretch & Fold timestamps in minutes)
		// The placeholder array inside initialInputs is sent as empty; BOBYQA will generate it dynamically
		final double[] placeholderFolds = new double[0];

		// 7. Establish the dynamic baseline context configuration
		final double maltSugarContent = 0.1;
		final DoughRecipe dummyBaseRecipe = new DoughRecipe(0.62, 0.022, 0.02, 0.002,
			maltSugarContent, (15_000. / 110.) * (1. - maltSugarContent), 0.05);

		final SimulationInputs initialInputs = new SimulationInputs(fractions, flourMatrix, flourTemperature,
			airRelativeHumidity, yeastProps, dummyBaseRecipe, kneading, stages, placeholderFolds);

		// Target profile parameter definitions from your custom enum
		final BakeryProduct selectedProduct = BakeryProduct.GASTRONOMY_PAN_PIZZA;
		System.out.println("Selected Target Product: " + selectedProduct.name());
		// PRINT CALCULATED TIMELINE INSIGHTS PRIOR TO EXECUTING OPTIMIZATION RAMP
		System.out.printf("Yeast Type Target Moisture : %.1f%%\n", targetYeastMoisture * 100.);
		System.out.printf("Autolyse                   : %.0f min\n", autolyseMinutes);
		final double[] window = yeastProps.getOptimalRehydrationWindow();
		System.out.printf("Rehydration Recommended    : %.1f min\n", yeastProps.getRehydrationWindow());
		System.out.printf("Rehydration Optimal Window : %.1f – %.1f min\n", window[0], window[1]);
		System.out.println("Scanning space...\n");

		// --- EXECUTE THE MULTIVARIATE SEARCH OVER CONTINUOUS BOUNDS ---
		final GlobalMultivariateRecipeOptimizer.OptimizedRecipeResult result =
			GlobalMultivariateRecipeOptimizer.optimizeFullRecipe(initialInputs, selectedProduct);

		// Reconstruct runtime fold vector array from optimized output points to run verification pass
		final int optimizedFoldCount = result.optimizedFoldCount();
		final double optimizedInitialDelay = result.optimizedFoldInitialDelay();
		final double optimizedRelaxationFactor = result.optimizedRelaxationFactor();
		final double[] optimizedFoldsArray = FullRecipeMultivariateObjective.buildFoldSchedule(optimizedFoldCount,
			optimizedInitialDelay, optimizedRelaxationFactor);

		// Reconstruct final optimized input payload for validation and diagnostic logging
		final SimulationInputs optimizedInputs = new SimulationInputs(fractions, flourMatrix, flourTemperature,
			airRelativeHumidity, yeastProps, result.recipe(), kneading, stages, optimizedFoldsArray);

		// Run final verification simulation at the optimum vector coordinate points
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(optimizedInputs);
		final double[] finalState = YeastOptimizer.runSimulation(optimizedInputs, gab, result.optimalYeastRatio(),
			selectedProduct);

		// --- GLOBAL OPTIMIZATION REPORT DISPLAY ---
		System.out.println("====================================================");
		System.out.println("           GLOBAL OPTIMIZED RECIPE                  ");
		System.out.println("====================================================");
		System.out.printf("Optimized Hydration  : %.1f%%\n", result.recipe().waterRatio() * 100.);
		System.out.printf("Optimized Salt Ratio : %.1f%%\n", result.recipe().saltRatio() * 100.);
		System.out.printf("Optimized Oil Ratio  : %.1f%%\n", result.recipe().oilRatio() * 100.);
		System.out.printf("Optimized Malt Ratio : %.2f%%\n", result.recipe().maltRatio() * 100.);
		System.out.printf("Optimal Starter Yeast: %.2f%%\n", result.optimalYeastRatio() * 100.);
		System.out.println("----------------------------------------------------");
		System.out.printf("Optimized S&F Count  : %d sets\n", optimizedFoldCount);
		if(optimizedFoldCount > 0){
			System.out.printf("Optimized First Delay: %.1f minutes after mixing\n", result.optimizedFoldInitialDelay());
			System.out.printf("Optimized Relaxation : %.2fx geometric scaling\n", optimizedRelaxationFactor);
			System.out.println("\n--- MECHANICAL STEP-BY-STEP WORKBOOK ---");
			for(int i = 0; i < optimizedFoldCount; i ++){
				BakeryProduct.HandlingInstruction stepInfo = selectedProduct.getHandlingInstructions((i + 1),
					optimizedFoldCount, optimizedFoldsArray[i]);
				System.out.printf("[TIMEMARK: %4.1f min] -> SET %d: %s (Intensity: %s)\n",
					optimizedFoldsArray[i], (i + 1), stepInfo.title(), stepInfo.intensity());
				System.out.printf("                    Action: %s\n\n", stepInfo.instructions());
			}
		}
		else
			System.out.println("Optimized S&F Timing : 0 sets (No mechanical folding required for this profile)");
		System.out.println("----------------------------------------------------");
		System.out.printf("Global Fitness Score : %.4f\n", result.finalFitness());

		if(finalState != null){
			System.out.println("\n--- FINAL DOUGH STATE AT TIMELINE EXPIRATION ---");
			System.out.printf("Yeast (Dry Biomass)  : %.2f%%\n", finalState[0] * 100.);
			System.out.printf("Sugars Remaining     : %.1f%%\n", finalState[1] * 100.);
			System.out.printf("Retained Gas Volume  : %.1f ml/g_flour\n", finalState[2]);
			System.out.printf("Ethanol Accumulation : %.1f g/kg_dough\n", finalState[3] * 1000.);

			DoughDiagnostics.report(optimizedInputs, finalState, selectedProduct);
		}
		else
			System.err.println("Critical Error: Core tracking verification integration run failed.");
		System.out.println("====================================================");
	}

}
