package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.FlourType;
import io.github.mtrevisan.yeastcalculator.domain.KneadingInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.domain.YeastInput;
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
public class Main2{

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
		final double flourTemperature = 19.;
		final double airRelativeHumidity = 0.54;

		// 3. Yeast properties
		final double targetYeastMoisture = 0.70;
		final YeastInput yeastProps = new YeastInput(targetYeastMoisture);

		// 4. Mechanical Kneading profiles
		final KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.MANUAL, 15.);

		// 5. Multi-stage fermentation schedule (Hard constraints set by user)
		final StageInput[] stages = {
			// 4.1 hours warm maturation profile at 28 °C
			new StageInput(28., 0.75, 4.1)
		};

		// 6. Physical structural interventions (Stretch & Fold timestamps in hours)
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
		System.out.printf("Rehydration Optimal Window : %s\n", yeastProps.getRehydrationWindow() * 60.);
		System.out.println("Scanning space...\n");

		// --- EXECUTE THE MULTIVARIATE SEARCH OVER CONTINUOUS BOUNDS ---
		final GlobalMultivariateRecipeOptimizer.OptimizedRecipeResult result =
			GlobalMultivariateRecipeOptimizer.optimizeFullRecipe(initialInputs, selectedProduct);

		// Reconstruct runtime fold vector array from optimized output points to run verification pass
		final double relaxationFactor = 1.4;
		final double[] optimizedFoldsArray = new double[result.optimizedFoldCount()];
		// Serves as initial fold delay field variable mapping
		double currentDelay = result.optimizedFoldInterval();
		double cumulativeMinutes = 0.;
		for(int i = 0; i < result.optimizedFoldCount(); i ++){
			cumulativeMinutes += currentDelay;
			optimizedFoldsArray[i] = cumulativeMinutes;
			currentDelay *= relaxationFactor;
		}

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
		System.out.printf("Optimized S&F Count  : %d sets\n", result.optimizedFoldCount());
		if(result.optimizedFoldCount() > 0){
			System.out.printf("Optimized First Delay: %.1f minutes after mixing\n", result.optimizedFoldInterval());
			System.out.println("Optimized S&F Style  : Geometric Relaxation Curve (1.4x factor scaling)");
			System.out.println("\n--- MECHANICAL STEP-BY-STEP WORKBOOK ---");
			for(int i = 0; i < result.optimizedFoldCount(); i ++){
				BakeryProduct.HandlingInstruction stepInfo = selectedProduct.getHandlingInstructions((i + 1),
					result.optimizedFoldCount(), optimizedFoldsArray[i]);
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
			System.out.printf("Yeast (Active Matter): %.2f%%\n", finalState[0] * 100.);
			System.out.printf("Sugars Remaining     : %.1f%%\n", finalState[1] * 100.);
			System.out.printf("Retained Gas Volume  : %.1f ml/g_flour\n", finalState[2]);
			System.out.printf("Ethanol Accumulation : %.1f g/kg_dough\n", finalState[3] * 1000.);

			// Calculate total timeline execution hours
			double totalDuration = 0.;
			for(final StageInput stage : stages)
				totalDuration += stage.duration();

			// Execute the dynamic rheological evaluation
			evaluateMaturationQuality(optimizedInputs, totalDuration);

			// Execute Biochemical & Sensory diagnostics
			System.out.println("\n--- BIOCHEMICAL & SENSORY DIAGNOSTICS ---");

			// 1. Maillard Reaction / Sugar Check
			if(finalState[1] < selectedProduct.getMinSafeSugarThreshold())
				System.out.println("[CRITICAL] Sugar levels dropped below safe baking charts thresholds! Crust will bake pale.");
			else
				System.out.println("[SUCCESS] Sugar levels safely satisfied target baking requirement benchmarks.");

			// 2. Yeast Residue / Off-Flavor Check
			final double finalYeastPercent = finalState[0] * 100.;
			final double maxAllowedYeastPercent = selectedProduct.getMaxAllowedFinalYeast() * 100.;

			if(finalState[0] > selectedProduct.getMaxAllowedFinalYeast()){
				System.out.printf("[CRITICAL] Sensory Alert: Residual active yeast (%.2f%%) exceeds the off-flavor profile threshold (Max: %.2f%%) for %s.\n",
					finalYeastPercent, maxAllowedYeastPercent, selectedProduct.name());
				System.out.println("           Result: The final product will have a strong, pungent chemical/yeasty smell and taste.");
			}
			else if(finalState[0] >= (selectedProduct.getMaxAllowedFinalYeast() * 0.85)){
				System.out.printf("[NOTICE] Pushing Sensory Limits: Residual active yeast (%.2f%%) is approaching the maximum allowed ceiling (%.2f%%).\n",
					finalYeastPercent, maxAllowedYeastPercent);
				System.out.println("         Result: Excellent oven-spring expected, but do not push the timeline any shorter to avoid taste degradation.");
			}
			else{
				System.out.printf("[SUCCESS] Clean Sensory Profile: Residual active yeast (%.2f%%) is well within safe bounds (Max: %.2f%%).\n",
					finalYeastPercent, maxAllowedYeastPercent);
				System.out.println("         Result: Clean, traditional fermentation aroma without heavy yeasty overtones.");
			}
		}
		else
			System.err.println("Critical Error: Core tracking verification integration run failed.");
		System.out.println("====================================================");
	}

	/**
	 * Evaluates the structural and enzymatic maturation quality of the dough
	 * based on flour strength (W), hydration ratio, and total process duration.
	 */
	private static void evaluateMaturationQuality(final SimulationInputs inputs, final double totalDuration){
		System.out.println("\n--- RHEOLOGICAL MATURATION DIAGNOSTICS ---");

		double blendW = 0.;
		final double[] fractions = inputs.fractions();
		final FlourInput[] flourMatrix = inputs.flourMatrix();
		for(int i = 0; i < fractions.length; i ++)
			blendW += fractions[i] * flourMatrix[i].strength();

		// [hours]
		double estimatedRequired = (blendW / 300.) * 6.;

		final double waterRatio = inputs.recipe().waterRatio();
		if(waterRatio >= 0.70)
			estimatedRequired *= 0.85;
		else if(waterRatio <= 0.55)
			estimatedRequired *= 1.15;

		if(totalDuration < estimatedRequired){
			System.out.printf("[WARNING] Insufficient maturation window for this flour strength (W: %.0f).\n", blendW);
			System.out.printf("          Required: ~%.1fh, Provided: %.1fh.\n", estimatedRequired, totalDuration);
			System.out.println("          Result: The gluten mesh will remain overly tense, causing high springback (elastic snap)");
			System.out.println("                  and potential gas retention instability during stretching.");
		}
		else if(totalDuration > (estimatedRequired * 2.5)){
			System.out.printf("[WARNING] Excessive room-temperature maturation window detected for W: %.0f.\n", blendW);
			System.out.printf("          Optimal window capped around ~%.1fh, Provided: %.1fh.\n", estimatedRequired * 2., totalDuration);
			System.out.println("          Result: Protease activity might over-degrade the gluten matrix structure,");
			System.out.println("                  leading to a sticky, fragile dough prone to tearing.");
		}
		else{
			System.out.printf("[SUCCESS] Maturation timeline (%.1fh) is well-proportioned for this blend (W: %.0f, Hydration: %.0f%%).\n",
				totalDuration, blendW, waterRatio * 100.);
			System.out.println("          Result: Optimal balance between gluten extensibility and gas holding matrix tenacity.");
		}
	}

}
