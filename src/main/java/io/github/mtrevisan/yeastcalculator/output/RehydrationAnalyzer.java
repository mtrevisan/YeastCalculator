package io.github.mtrevisan.yeastcalculator.output;

import io.github.mtrevisan.yeastcalculator.model.DoughComposition;
import io.github.mtrevisan.yeastcalculator.model.FermentationSchedule;
import io.github.mtrevisan.yeastcalculator.model.YeastProperties;
import io.github.mtrevisan.yeastcalculator.optimization.YeastOptimizer;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;


/** Effect of yeast rehydration time on lag phase and optimal yeast fraction. */
public class RehydrationAnalyzer{

	private final DoughComposition dough;
	private final FermentationSchedule schedule;
	private final double yeastMoisture;


	public RehydrationAnalyzer(final DoughComposition dough, final FermentationSchedule schedule,
			final double yeastMoisture){
		this.dough = dough;
		this.schedule = schedule;
		this.yeastMoisture = yeastMoisture;
	}

	/**
	 * Print table: rehydration time vs. Q0_eff, lag, and optimal yeast.
	 */
	public void printRehydrationEffect(){
		System.out.printf(Locale.US,
			"%-13s %-9s %-1s%n",
			"t_rehyd[min]", "lag[min]", "yeast(WB)");

		for(final double tMin : new double[]{0., 5., 10., 15., 20., 30., 45., 60.}){
			final YeastProperties yeastRehydrated = new YeastProperties(yeastMoisture, tMin / 60.);
			final YeastOptimizer optimizer = new YeastOptimizer(dough, schedule, yeastRehydrated);

			// Suppress warnings
			final boolean hasWarning = (optimizer.estimatedLagTime() * 60. > 0.5 * schedule.getTotalDuration());
			final PrintStream originalErr = System.err;
			if(hasWarning)
				System.setErr(new PrintStream(OutputStream.nullOutputStream()));

			final double optimalYeast = optimizer.findOptimal(false);

			if(hasWarning)
				System.setErr(originalErr);

			System.out.printf(Locale.US,
				"%-13.0f %-9d %-1.2f%%%s%n",
				tMin,
				(int)Math.ceil(optimizer.estimatedLagTime() * 60),
				optimizer.toBakersWetPercent(optimalYeast),
				(hasWarning? " ⚠️": ""));
		}
	}

}