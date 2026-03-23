package io.github.mtrevisan.yeastcalculator.output;

import io.github.mtrevisan.yeastcalculator.optimization.YeastOptimizer;

import java.util.Locale;


/** Sensitivity analysis: fermentation progress vs. the initial yeast fraction. */
public class SensitivityAnalyzer{

	private final YeastOptimizer optimizer;


	public SensitivityAnalyzer(final YeastOptimizer optimizer){
		this.optimizer = optimizer;
	}

	/**
	 * Print sensitivity table: factor vs. progress.
	 */
	public void printSensitivity(final double optimalAnhydrousYeast){
		System.out.printf(Locale.US,
			"%-7s %-10s %-1s%n",
			"factor", "yeast(WB)", "progress");

		for(final double factor : new double[]{0.25, 0.5, 0.75, 1., 1.5, 2.}){
			final double testFraction = Math.min(
				optimalAnhydrousYeast * factor,
				optimizer.getSubstrateMaxBiomassCapacity());

			System.out.printf(Locale.US,
				"%-7.2f   %-1.2f%%      %-1.2f%n",
				factor,
				optimizer.toBakersWetPercent(testFraction),
				optimizer.computeFermentationProgress(testFraction));
		}
	}

}