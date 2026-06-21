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
			new FlourInput(295., 0.55, 0.013, 0.13, 0.011, 0.019, 0.003, FlourType.WHEAT)
		};

		// 2. Structural & Storage baselines
		final double flourTemperature = 19.;
		final double airRelativeHumidity = 0.54;

		// 3. Yeast properties
		final YeastInput yeastProps = new YeastInput(0.70, 5. / 60.);

		// 4. Dough Recipe Configuration
		final DoughRecipe recipe = new DoughRecipe(
			0.60,	// Water Ratio
			0.022,	// Salt Ratio
			0.05,		// Oil Ratio
			0.008,	// Malt Ratio
			0.5,		// Malt Sugar content
			80.,		// Diastatic power (Pollak Units)
			3.			// Mixer Friction Factor
		);

		// 5. Mechanical Kneading profiles
		final KneadingInput kneading = new KneadingInput(KneadingInput.KneadingType.MANUAL, 12.);

		// 6. Multi-stage fermentation schedule
		final StageInput[] stages = {
			new StageInput(28., 0.55, 3.) // 3 hours at 28 degrees Celsius
		};

		// 7. Physical structural interventions (Stretch & Fold timestamps in hours)
		final double[] folds = {0.75, 1.5}; // Folds execution at 45m and 90m

		// Build composite simulation payload object
		final SimulationInputs inputs = new SimulationInputs(
			fractions, flourMatrix, flourTemperature, airRelativeHumidity,
			yeastProps, recipe, kneading, stages, folds
		);

		// Target profile parameter definitions
		final BakeryProduct selectedProduct = BakeryProduct.GASTRONOMY_PAN_PIZZA;

		System.out.println("Selected Target Product: " + selectedProduct.name());
		System.out.println("Processing optimization calculations via single-variable Brent inversion...");

		// Execute Inversion Optimization Target Calculation
		final double optimalYeastRatio = YeastOptimizer.findOptimalYeast(inputs, selectedProduct);

		System.out.println("\nOptimization completed successfully.");
		System.out.printf("Optimal Yeast Ratio Target Required: %.3f%% (relative to total flour)\n", optimalYeastRatio * 100.);

		// Run confirmation simulation at target values to extract chemical endpoints
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(inputs);
		final double maxGasPotential = YeastOptimizer.calculateMaxGasPotential(new double[]{295., 0.55, 0.011, 0.003},
			inputs, selectedProduct);
		final double[] finalState = YeastOptimizer.runSimulation(inputs, gab, maxGasPotential, optimalYeastRatio);

		if(finalState != null){
			System.out.println("\n--- FINAL DOUGH STATE AT TIMELINE EXPIRATION ---");
			System.out.printf("Active Dry Yeast          : %.3f%%\n", finalState[0] * 100.);
			System.out.printf("Residual Sugars Remaining : %.2f%%\n", finalState[1] * 100.);
			System.out.printf("Final Retained Gas Volume : %.2f mL/g_flour\n", finalState[2]);

			if(finalState[1] < selectedProduct.getMinSafeSugarThreshold())
				System.out.println("WARNING: Sugar levels dropped below safe baking charts thresholds!");
			else
				System.out.println("Status: Sugar levels safely satisfied target baking requirement benchmarks.");
		}
		else
			System.err.println("Critical Error: Core tracking verification integration run failed.");
	}

}
