package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.model.FlourMoistureModel;
import io.github.mtrevisan.yeastcalculator.output.SensitivityAnalyzer;
import io.github.mtrevisan.yeastcalculator.output.SimulationTracer;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;


public class Main{

	public static void main(final String[] args){
		final double[] fractions = {1.};
		//W	P/L	carbo.	sugar	prot.	fat	fiber	ashes	salt
		final double[][] flourMatrix = new double[][]{
			{170., 0.7, 0.72, 0.017, 0.11, 0.007, 0.022, 0.005, 0.}
		};
		final String[] flourTypes = {"wheat"};
		final double flourTemperature = 17.8;
		final double airRelativeHumidity = 0.54;
		final double[] moisture = FlourMoistureModel.compute(fractions, flourMatrix, flourTypes,
			flourTemperature, airRelativeHumidity);
		final double flourFreeWater = moisture[0] - moisture[1];

		//── Dough composition ──────────────────────────────────────────────
		//total = flour + water + salt + sugar
		final double water = 0.65;
		final double salt = 0.012;
		final double sugar = 0.02;
		//extra ingredients, such as oil, malt, etc
		final double extra = 0.;

		//── Fermentation schedule ──────────────────────────────────────────
		//two-stage: bulk fermentation + proofing
		final double[] temperatures = {28., 30.};
		final double[] durations = {2., 2.};

		//── Yeast: fresh compressed ────────────────────────────────────────
		//the model computes Q0 from this value and returns an anhydrous mass
		final double yeastMoisture = 0.7;
		// rehydrationDuration = time spent soaking in warm water before mixing
		final double rehydrationDuration = 20. / 60.;

		//── Instantiate model ──────────────────────────────────────────────
		final YeastCalculator model = new YeastCalculator(
			water, salt, sugar, extra,
			temperatures, durations,
			flourFreeWater, yeastMoisture,
			rehydrationDuration);

		//── Optimize ──────────────────────────────────────────────────────
		final double optimalAnhydrousYeast = model.findOptimal(true);

		System.out.printf("%n=== RESULT ===%n");
		System.out.printf(Locale.US,
			"Lag in dough       : %d min%n", model.estimatedLagTime());
		System.out.printf(Locale.US,
			"Optimal yeast (WB) : %.2f%%%n",
			model.yeastWB(optimalAnhydrousYeast));

		//── Diagnostic trace ──────────────────────────────────────────────
		System.out.printf("%n=== SIMULATION TRACE ===%n");
		final SimulationTracer tracer = model.getSimulationTracer();
		tracer.printTrace(optimalAnhydrousYeast);

		//── Sensitivity analysis ──────────────────────────────────────────
		System.out.printf("%n=== SENSITIVITY ===%n");
		final SensitivityAnalyzer sensitivityAnalyzer = model.getSensitivityAnalyzer();
		sensitivityAnalyzer.printSensitivity(optimalAnhydrousYeast);

		//── Effect of rehydration time on lag (same yeast, same schedule) ─────────
		System.out.printf("%n=== EFFECT OF REHYDRATION TIME ===%n");
		System.out.printf(Locale.US,
			"%-13s %-9s %-1s%n",
			"t_rehyd[min]", "lag[min]", "yeast(WB)");
		for(final double tMin : new double[]{0., 5., 10., 15., 20., 30., 45., 60.}){
			final YeastCalculator mr = new YeastCalculator(
				water, salt, sugar, extra,
				temperatures, durations,
				flourFreeWater, yeastMoisture,
				tMin / 60.);
			final boolean hasWarning = mr.hasWarning();
			final PrintStream originalErr = System.err;
			if(hasWarning)
				System.setErr(new PrintStream(OutputStream.nullOutputStream()));
			final double n0r = mr.findOptimal(false);
			if(hasWarning)
				System.setErr(originalErr);
			System.out.printf(Locale.US,
				"%-13.0f %-9d %-1.2f%%%s%n",
				tMin,
				mr.estimatedLagTime(),
				mr.yeastWB(n0r), (hasWarning? " ⚠️": ""));
		}
	}

}
