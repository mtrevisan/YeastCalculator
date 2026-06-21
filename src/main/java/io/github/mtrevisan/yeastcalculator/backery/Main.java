package io.github.mtrevisan.yeastcalculator.backery;

import java.util.Locale;


public class Main{

	public static void main(String[] args){
		// Enforce dot formatting for clean metric output readings
		Locale.setDefault(Locale.US);

		System.out.println("====================================================");
		System.out.println("   BAKERY ODE ISOTHERM OPTIMIZATION ENGINE ENGINE   ");
		System.out.println("====================================================\n");

		// 1. Definition of the Flour Blend
		double[] fractions = {1.0};
		FlourInput[] flourMatrix = {
			new FlourInput(295.0, 0.55, 0.013, 0.13, 0.011, 0.019, 0.003, FlourType.WHEAT)
		};

		// 2. Structural & Storage baselines
		double flourTemperature = 19.0;
		double airRelativeHumidity = 0.54;

		// 3. Yeast properties
		YeastInput yeastProps = new YeastInput(0.70, 5.0 / 60.0);

		// 4. Dough Recipe Configuration
		DoughRecipe recipe = new DoughRecipe(
			0.60,  // Water Ratio
			0.022, // Salt Ratio
			0.05,  // Oil Ratio
			0.008, // Malt Ratio
			0.5,   // Malt Sugar content
			80.0,  // Diastatic power (Pollak Units)
			3.0    // Mixer Friction Factor
		);

		// 5. Mechanical Kneading profiles
		KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.SPIRAL_MIXER, 12.0);

		// 6. Multi-stage fermentation schedule
		StageInput[] stages = {
			new StageInput(28.0, 0.55, 3.0) // 3 hours at 28 degrees Celsius
		};

		// 7. Physical structural interventions (Stretch & Fold timestamps in hours)
		double[] folds = {0.75, 1.5}; // Folds execution at 45m and 90m

		// Build composite simulation payload object
		SimulationInputs inputs = new SimulationInputs(
			fractions, flourMatrix, flourTemperature, airRelativeHumidity,
			yeastProps, recipe, kneading, stages, folds
		);

		// Target target profile parameter definitions
		BakeryProduct selectedProduct = BakeryProduct.BREAD;

		System.out.println("Selected Target Product: " + selectedProduct.name());
		System.out.println("Processing optimization calculations via single-variable Brent inversion...");

		// Execute Inversion Optimization Target Calculation
		double optimalYeastRatio = YeastOptimizer.findOptimalYeast(inputs, selectedProduct);

		System.out.println("\nOptimization completed successfully.");
		System.out.printf("Optimal Yeast Ratio Target Required: %.4f%% (relative to total flour)\n", optimalYeastRatio * 100.0);

		// Run confirmation simulation at target values to extract chemical endpoints
		GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(inputs);
		double maxGasPotential = YeastOptimizer.calculateMaxGasPotential(inputs, selectedProduct);
		double[] finalState = YeastOptimizer.runSimulation(inputs, gab, maxGasPotential, optimalYeastRatio);

		if(finalState != null){
			System.out.println("\n--- FINAL DOUGH STATE AT TIMELINE EXPIRATION ---");
			System.out.printf("Active Yeast Biomass (Dry Basis) : %.6f g/g_flour\n", finalState[0]);
			System.out.printf("Residual Sugar Content           : %.4f%% (g/g_flour)\n", finalState[1] * 100.0);
			System.out.printf("Trapped Gas Volume               : %.2f mL/g_flour\n", finalState[2]);

			if(finalState[1] < selectedProduct.getMinSafeSugarThreshold()){
				System.out.println("WARNING: Sugar levels dropped below safe baking charts thresholds!");
			}
			else{
				System.out.println("Status: Sugar levels safely satisfied target baking requirement benchmarks.");
			}
		}
		else{
			System.err.println("Critical Error: Core tracking verification integration run failed.");
		}
	}

}
