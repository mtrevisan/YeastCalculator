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


public class Main{

	public static void main(final String[] args){
		// Enforce dot formatting for clean metric output readings
		Locale.setDefault(Locale.US);

		System.out.println("====================================================");
		System.out.println("   BAKERY ODE ISOTHERM OPTIMIZATION ENGINE ENGINE   ");
		System.out.println("====================================================\n");

		// 1. Definition of the Flour Blend
		final double[] fractions = {1.};
		final FlourInput[] flourMatrix = {
//			new FlourInput(295., 0.55, 0.013, 0.13, 0.011, 0.019, 0.003, FlourType.WHEAT)
			new FlourInput(205., 0.523, 0.014, 0.12, 0.01, 0.017, 0.003, FlourType.WHEAT)
		};

		// 2. Structural & Storage baselines
		final double flourTemperature = 19.;
		final double airRelativeHumidity = 0.54;

		// 3. Yeast properties
		final YeastInput yeastProps = new YeastInput(0.70, 5. / 60.);

		// 4. Dough Recipe Configuration
		final double maltSugarContent = 0.1;
		final DoughRecipe recipe = new DoughRecipe(
			0.62,	// Water Ratio
			0.022,	// Salt Ratio
			0.07,		// Oil Ratio
			0.009,	// Malt Ratio
			maltSugarContent,
			(15_000. / 110.) * (1. - maltSugarContent),	// Diastatic power
			0.05		// Friction Factor
		);

		// 5. Mechanical Kneading profiles
		final KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.MANUAL, 15.);

		// 6. Multi-stage fermentation schedule
		final StageInput[] stages = {
			new StageInput(28., 0.75, 4.5)
		};

		// 7. Physical structural interventions (Stretch & Fold timestamps in hours)
		final double[] folds = {0.5, 1., 1.5}; // Folds execution at 45m and 90m

		// Build composite simulation payload object
		final SimulationInputs inputs = new SimulationInputs(
			fractions, flourMatrix, flourTemperature, airRelativeHumidity,
			yeastProps, recipe, kneading, stages, folds
		);

		// Target profile parameter definitions
		final BakeryProduct selectedProduct = BakeryProduct.GASTRONOMY_PAN_PIZZA;

		System.out.println("Selected Target Product: " + selectedProduct.name());

		// Execute Inversion Optimization Target Calculation
		final double optimalYeastRatio = YeastOptimizer.findOptimalYeast(inputs, selectedProduct);

		System.out.printf("Optimal Yeast: %.2f%%\n", optimalYeastRatio * 100.);

		// Run confirmation simulation at target values to extract chemical endpoints
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(inputs);
		final double[] finalState = YeastOptimizer.runSimulation(inputs, gab, optimalYeastRatio, selectedProduct);

		if(finalState != null){
			System.out.println("\n--- FINAL DOUGH STATE AT TIMELINE EXPIRATION ---");
			System.out.printf("Yeast                : %.2f%%\n", finalState[0] * 100.);
			System.out.printf("Sugars               : %.2f%%\n", finalState[1] * 100.);
			System.out.printf("Retained Gas Volume  : %.1f ml/g_flour\n", finalState[2]);
			System.out.printf("Ethanol Accumulation : %.1f g/kg_dough\n", finalState[3] * 1000.);

			// Calculate total timeline execution hours
			double totalDurationHours = 0.;
			for(final StageInput stage : stages)
				totalDurationHours += stage.getDuration();

			// Execute the dynamic rheological evaluation
			evaluateMaturationQuality(inputs, totalDurationHours);

			// Execute Maillard / Sugar verification
			System.out.println("\n--- BIOCHEMICAL & SENSORY DIAGNOSTICS ---");
			// 1. Maillard Reaction / Sugar Check
			if(finalState[1] < selectedProduct.getMinSafeSugarThreshold())
				System.out.println("[CRITICAL] Sugar levels dropped below safe baking charts thresholds! Crust will bake pale.");
//			else
//				System.out.println("[SUCCESS] Sugar levels safely satisfied target baking requirement benchmarks.");
			// 2. NEW: Yeast Residue / Off-Flavor Check
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
//			else{
//				System.out.printf("[SUCCESS] Clean Sensory Profile: Residual active yeast (%.2f%%) is well within safe bounds (Max: %.2f%%).\n",
//					finalYeastPercent, maxAllowedYeastPercent);
//				System.out.println("         Result: Clean, traditional fermentation aroma without heavy yeasty overtones.");
//			}
		}
		else
			System.err.println("Critical Error: Core tracking verification integration run failed.");
	}

	/**
	 * Evaluates the structural and enzymatic maturation quality of the dough
	 * based on flour strength (W), hydration ratio, and total process duration.
	 * * @param inputs             The composite simulation payload containing flour and recipe data.
	 *
	 * @param totalDurationHours The cumulative duration of all fermentation stages.
	 */
	private static void evaluateMaturationQuality(final SimulationInputs inputs, final double totalDurationHours){
		System.out.println("\n--- RHEOLOGICAL MATURATION DIAGNOSTICS ---");

		// 1. Compute the weighted average of Flour Strength (W) from the blend
		double blendW = 0.;
		final double[] fractions = inputs.getFractions();
		final FlourInput[] flourMatrix = inputs.getFlourMatrix();
		for(int i = 0; i < fractions.length; i ++)
			blendW += fractions[i] * flourMatrix[i].getStrength();

		// 2. Estimate base required maturation hours at room temperature as a function of W
		// Standard benchmark: A W300 flour requires roughly 6 hours at 24-28°C for full protease relaxation.
		double estimatedRequiredHours = (blendW / 300.) * 6.;

		// 3. Kinetic correction factor based on hydration (Water Ratio)
		// Higher hydration increases enzymatic mobility (amylases/proteases), accelerating maturation.
		final double waterRatio = inputs.getRecipe().getWaterRatio();
		if(waterRatio >= 0.70)
			// Accelerate by 15% due to high enzymatic diffusion
			estimatedRequiredHours *= 0.85;
		else if(waterRatio <= 0.55)
			// Deccelerate by 15% due to high osmotic/viscous restriction
			estimatedRequiredHours *= 1.15;

		// 4. Structural evaluation and logging output
		if(totalDurationHours < estimatedRequiredHours){
			System.out.printf("[WARNING] Insufficient maturation window for this flour strength (W: %.0f).\n", blendW);
			System.out.printf("          Required: ~%.1fh, Provided: %.1fh.\n", estimatedRequiredHours, totalDurationHours);
			System.out.println("          Result: The gluten mesh will remain overly tense, causing high springback (elastic snap)");
			System.out.println("                  and potential gas retention instability during stretching.");
		}
		else if(totalDurationHours > (estimatedRequiredHours * 2.5)){
			System.out.printf("[WARNING] Excessive room-temperature maturation window detected for W: %.0f.\n", blendW);
			System.out.printf("          Optimal window capped around ~%.1fh, Provided: %.1fh.\n", estimatedRequiredHours * 2., totalDurationHours);
			System.out.println("          Result: Protease activity might over-degrade the gluten matrix structure,");
			System.out.println("                  leading to a sticky, fragile dough prone to tearing.");
		}
//		else{
//			System.out.printf("[SUCCESS] Maturation timeline (%.1fh) is well-proportioned for this blend (W: %.0f, Hydration: %.0f%%).\n",
//				totalDurationHours, blendW, waterRatio * 100.);
//			System.out.println("          Result: Optimal balance between gluten extensibility and gas holding matrix tenacity.");
//		}
	}

}
