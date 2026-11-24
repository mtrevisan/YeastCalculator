package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;


/**
 * Console diagnostics shared by the entry points: they judge the final state of a simulated dough.
 */
final class DoughDiagnostics{

	private DoughDiagnostics(){}


	/**
	 * Prints the rheological and the biochemical/sensory diagnostics of the dough at the end of the timeline.
	 */
	static void report(final SimulationInputs inputs, final double[] finalState, final BakeryProduct product){
		double totalDuration = 0.;
		for(final StageInput stage : inputs.stages())
			totalDuration += stage.duration();

		evaluateMaturationQuality(inputs, totalDuration);
		evaluateBiochemistry(finalState, product);
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
			System.out.printf("          Required: ~%.2f h, Provided: %.2f h.\n", estimatedRequired, totalDuration);
			System.out.println("          Result: The gluten mesh will remain overly tense, causing high springback (elastic snap)");
			System.out.println("                  and potential gas retention instability during stretching.");
		}
		else if(totalDuration > estimatedRequired * 2.5){
			System.out.printf("[WARNING] Excessive room-temperature maturation window detected for W: %.0f.\n", blendW);
			System.out.printf("          Optimal window capped around ~%.1f h, Provided: %.1f h.\n", estimatedRequired * 2., totalDuration);
			System.out.println("          Result: Protease activity might over-degrade the gluten matrix structure,");
			System.out.println("                  leading to a sticky, fragile dough prone to tearing.");
		}
		else{
			System.out.printf("[SUCCESS] Maturation timeline (%.1f h) is well-proportioned for this blend (W: %.0f, Hydration: %.0f%%).\n",
				totalDuration, blendW, waterRatio * 100.);
			System.out.printf("          Estimated ~%.1f h\n", estimatedRequired);
			System.out.println("          Result: Optimal balance between gluten extensibility and gas holding matrix tenacity.");
		}
	}

	private static void evaluateBiochemistry(final double[] finalState, final BakeryProduct product){
		System.out.println("\n--- BIOCHEMICAL & SENSORY DIAGNOSTICS ---");

		// 1. Maillard Reaction / Sugar Check
		if(finalState[1] < product.getMinSafeSugarThreshold())
			System.out.println("[CRITICAL] Sugar levels dropped below safe baking charts thresholds! Crust will bake pale.");
		else
			System.out.println("[SUCCESS] Sugar levels safely satisfied target baking requirement benchmarks.");

		// 2. Yeast Residue / Off-Flavor Check
		final double finalYeastPercent = finalState[0] * 100.;
		final double maxAllowedYeastPercent = product.getMaxAllowedFinalYeast() * 100.;

		if(finalState[0] > product.getMaxAllowedFinalYeast()){
			System.out.printf("[CRITICAL] Sensory Alert: Residual active yeast (%.2f%%) exceeds the off-flavor profile threshold (Max: %.2f%%) for %s.\n",
				finalYeastPercent, maxAllowedYeastPercent, product.name());
			System.out.println("           Result: The final product will have a strong, pungent chemical/yeasty smell and taste.");
		}
		else if(finalState[0] >= (product.getMaxAllowedFinalYeast() * 0.85)){
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

}
