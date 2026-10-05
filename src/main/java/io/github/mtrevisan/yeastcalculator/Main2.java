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
import io.github.mtrevisan.yeastcalculator.simulation.DoughOdeSystem;
import io.github.mtrevisan.yeastcalculator.simulation.FoldEventHandler;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;

import java.util.Locale;


public class Main2{

	public static void main(final String[] args){
		// Enforce dot formatting for clean metric output readings
		Locale.setDefault(Locale.US);

		// 1. Definition of the Flour Blend
		final double[] fractions = {1.};
		final FlourInput[] flourMatrix = {
			new FlourInput(205., 0.523, 0.014, 0.12, 0.01, 0.017, 0.003, FlourType.WHEAT)
		};

		// 2. Structural & Storage baselines
		final double flourTemperature = 19.;
		final double airRelativeHumidity = 0.54;

		// 3. Yeast properties
		final YeastInput yeastProps = new YeastInput(0.70, 5. / 60.);

		// 4. Mechanical Kneading profiles
		final KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.MANUAL, 15.);

		// 5. Multi-stage fermentation schedule (Hard constraints set by user)
		final StageInput[] stages = {
			new StageInput(28., 0.75, 4.1)
		};

		// 6. Physical structural interventions (Stretch & Fold timestamps in hours)
		final double[] folds = {0.5, 1., 1.5};

		// Target profile parameter definitions
		final BakeryProduct selectedProduct = BakeryProduct.GASTRONOMY_PAN_PIZZA;
		System.out.println("Selected Target Product: " + selectedProduct.name());
		System.out.println("Scanning ingredient space for global biochemical & sensory optimum at 28 °C...");

		// --- GRID SEARCH EXPLORATION SPACES ---
		final double targetHydration = 0.62;
		final double[] saltOptions = {0.018, 0.020, 0.022, 0.024, 0.026, 0.028, 0.030};
		final double[] oilOptions = {0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090};
		final double[] maltOptions = {0.000, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010};

		double bestFitness = Double.MAX_VALUE;
		DoughRecipe bestRecipe = null;
		double bestYeastRatio = 0.;
		double[] bestFinalState = null;

		// Execute nested combinatorics scans
		for(final double salt : saltOptions){
			for(final double oil : oilOptions){
				for(final double malt : maltOptions){

					final DoughRecipe candidateRecipe = new DoughRecipe(targetHydration, salt, oil, malt, 0.5,
						15_000., 0.05);

					final SimulationInputs candidateInputs = new SimulationInputs(fractions, flourMatrix, flourTemperature,
						airRelativeHumidity, yeastProps, candidateRecipe, kneading, stages, folds);

					// Run the univariate optimization inside the grid loop
					final double candidateYeast = YeastOptimizer.findOptimalYeast(candidateInputs, selectedProduct);

					final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(candidateInputs);
					final double[] finalState = YeastOptimizer.runSimulation(candidateInputs, gab, candidateYeast,
						selectedProduct);

					if(finalState == null)
						continue;

					// Extract structural parameters for internal enum fitness calculation
					final double[] summary = YeastOptimizer.inputsSummary(candidateInputs);
					final double baseMaxGasPotential = YeastOptimizer.calculateMaxGasPotential(summary, candidateInputs,
						selectedProduct);

					final StageInput lastStage = stages[stages.length - 1];
					final double skinningModifier = (lastStage.getRelativeHumidity() >= 0.70
						? 1.
						: Math.max(0.6, 1. - (0.70 - lastStage.getRelativeHumidity()) * 0.8));
					final double dynamicMaxPotential = baseMaxGasPotential * Math.pow(1.05, folds.length) * skinningModifier;

					// Reconstruct derivatives endpoint (yDot)
					final double[] yDotFinal = new double[4];
					double currentTime = 0.;
					for(final StageInput stage : stages)
						currentTime += stage.getDuration();

					final FoldEventHandler foldHandler = new FoldEventHandler(folds);
					final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.getTemperature(),
						lastStage.getRelativeHumidity(), gab.flourActiveWater, baseMaxGasPotential, candidateInputs,
						foldHandler);
					finalOde.computeDerivatives(currentTime * 60., finalState, yDotFinal);

					final double currentFitness = selectedProduct.computeFitness(finalState[2], dynamicMaxPotential,
						finalState[1], finalState[0], yDotFinal[2]);

					// Retain the absolute global minimum of the loss/fitness function
					if(currentFitness < bestFitness){
						bestFitness = currentFitness;
						bestRecipe = candidateRecipe;
						bestYeastRatio = candidateYeast;
						bestFinalState = finalState;
					}
				}
			}
		}

		// Build final optimized input payload for diagnostic logging
		final SimulationInputs optimizedInputs = new SimulationInputs(fractions, flourMatrix, flourTemperature,
			airRelativeHumidity, yeastProps, bestRecipe, kneading, stages, folds);

		// --- GLOBAL OPTIMIZATION REPORT DISPLAY ---
		System.out.println("\n====================================================");
		System.out.println("              GLOBAL OPTIMIZED RECIPE               ");
		System.out.println("====================================================");
		System.out.printf("Optimized Salt Ratio : %.2f%%\n", bestRecipe.getSaltRatio() * 100.);
		System.out.printf("Optimized Oil Ratio  : %.2f%%\n", bestRecipe.getOilRatio() * 100.);
		System.out.printf("Optimized Malt Ratio : %.2f%%\n", bestRecipe.getMaltRatio() * 100.);
		System.out.printf("Optimal Starter Yeast: %.2f%%\n", bestYeastRatio * 100.);
		System.out.printf("Global Fitness Score : %.4f\n", bestFitness);

		System.out.println("\n--- FINAL DOUGH STATE AT TIMELINE EXPIRATION ---");
		System.out.printf("Yeast                : %.2f%%\n", bestFinalState[0] * 100.);
		System.out.printf("Sugars               : %.2f%%\n", bestFinalState[1] * 100.);
		System.out.printf("Retained Gas Volume  : %.1f ml/g_flour\n", bestFinalState[2]);
		System.out.printf("Ethanol Accumulation : %.1f g/kg_dough\n", bestFinalState[3] * 1000.);

		// Calculate total timeline execution hours
		double totalDuration = 0.;
		for(final StageInput stage : stages)
			totalDuration += stage.getDuration();

		// Execute the dynamic rheological evaluation
		evaluateMaturationQuality(optimizedInputs, totalDuration);

		// Execute Biochemical & Sensory diagnostics
		System.out.println("\n--- BIOCHEMICAL & SENSORY DIAGNOSTICS ---");

		// 1. Maillard Reaction / Sugar Check
		if(bestFinalState[1] < selectedProduct.getMinSafeSugarThreshold())
			System.out.println("[CRITICAL] Sugar levels dropped below safe baking charts thresholds! Crust will bake pale.");
		else
			System.out.println("[SUCCESS] Sugar levels safely satisfied target baking requirement benchmarks.");

		// 2. Yeast Residue / Off-Flavor Check
		final double finalYeastPercent = bestFinalState[0] * 100.;
		final double maxAllowedYeastPercent = selectedProduct.getMaxAllowedFinalYeast() * 100.;

		if(bestFinalState[0] > selectedProduct.getMaxAllowedFinalYeast()){
			System.out.printf("[CRITICAL] Sensory Alert: Residual active yeast (%.2f%%) exceeds the off-flavor profile threshold (Max: %.2f%%) for %s.\n",
				finalYeastPercent, maxAllowedYeastPercent, selectedProduct.name());
			System.out.println("           Result: The final product will have a strong, pungent chemical/yeasty smell and taste.");
		}
		else if(bestFinalState[0] >= (selectedProduct.getMaxAllowedFinalYeast() * 0.85)){
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

	/**
	 * Evaluates the structural and enzymatic maturation quality of the dough
	 * based on flour strength (W), hydration ratio, and total process duration.
	 */
	private static void evaluateMaturationQuality(final SimulationInputs inputs, final double totalDuration){
		System.out.println("\n--- RHEOLOGICAL MATURATION DIAGNOSTICS ---");

		double blendW = 0.;
		final double[] fractions = inputs.getFractions();
		final FlourInput[] flourMatrix = inputs.getFlourMatrix();
		for(int i = 0; i < fractions.length; i ++)
			blendW += fractions[i] * flourMatrix[i].getStrength();

		double estimatedRequired = (blendW / 300.) * 6.;

		final double waterRatio = inputs.getRecipe().getWaterRatio();
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
