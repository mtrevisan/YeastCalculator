package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.kinetics.CardinalTemperatureModel;
import io.github.mtrevisan.yeastcalculator.kinetics.KineticParameters;
import io.github.mtrevisan.yeastcalculator.model.DoughComposition;
import io.github.mtrevisan.yeastcalculator.model.FermentationSchedule;
import io.github.mtrevisan.yeastcalculator.model.YeastProperties;
import io.github.mtrevisan.yeastcalculator.ode.OdeSolver;
import io.github.mtrevisan.yeastcalculator.physiology.PhysiologicalState;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.NoBracketingException;

import java.util.Locale;


/**
 * Finds the optimal initial yeast inoculum using BrentSolver.
 * <p>
 * Target: fermentation reaches 99% of the substrate-limited maximum at t_final.
 */
public class YeastOptimizer{

	private final DoughComposition dough;
	private final FermentationSchedule schedule;
	private final YeastProperties yeast;

	private final double substrateLimitedBiomassCapacity;
	private final OdeSolver odeSolver;


	public YeastOptimizer(final DoughComposition dough, final FermentationSchedule schedule,
			final YeastProperties yeast){
		this.dough = dough;
		this.schedule = schedule;
		this.yeast = yeast;

		// N_MAX = YIELD_BIOMASS_PER_SUGAR × sugar [g_biomass / g_dough]
		this.substrateLimitedBiomassCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * dough.sugar;
		odeSolver = new OdeSolver(dough, schedule);
	}

	/**
	 * Find optimal anhydrous yeast inoculum.
	 * <p>
	 * Uses BrentSolver to find:
	 * 	N0* = argmin { |progress(N0) − TARGET_FERMENTATION_PROGRESS| }
	 * 		 = root of f(N0) = progress(N0) − TARGET_FERMENTATION_PROGRESS
	 * <p>
	 * Search domain: (0, N_MAX_substrate)
	 * 	Lower bound: effectively zero yeast — 1e-8 g/g
	 * 	Upper bound: substrate ceiling — N_MAX_substrate × 0.9999 (avoids capacity_adj = 0)
	 *
	 * @param outputWarning	If true, print diagnostic messages when the substrate ceiling is reached.
	 * @return	Optimal initial anhydrous yeast fraction [g_dry_yeast / g_dough].
	 * 	Use {@link #toBakersWetPercent} or {@link #toBakersWetPercent} to convert to the wet-yeast quantity to weigh
	 * 	on a scale.
	 */
	public double findOptimal(final boolean outputWarning){
		final String warning = scheduleWarning();
		if(warning != null)
			System.err.println("WARNING: " + warning);

		final double lowerBound = 1e-8;
		final double upperBound = getSubstrateMaxBiomassCapacity();

		final double progressAtLowerBound = computeFermentationProgress(lowerBound);
		final double progressAtUpperBound = computeFermentationProgress(upperBound);

		//check whether the bracket is valid for BrentSolver
		if(progressAtUpperBound < KineticParameters.TARGET_FERMENTATION_PROGRESS){
			if(outputWarning)
				printSubstrateCeilingDiagnostics(upperBound, progressAtUpperBound);
			return upperBound;
		}

		if(progressAtLowerBound >= KineticParameters.TARGET_FERMENTATION_PROGRESS){
			System.err.printf("WARNING: near-zero yeast already reaches the target. Stages are too long.%n");
			return lowerBound;
		}

		//f(N0) = progress(N0) − TARGET is monotone increasing on [lo, hi]
		final UnivariateFunction objectiveFunction
			= N0 -> computeFermentationProgress(N0) - KineticParameters.TARGET_FERMENTATION_PROGRESS;

		final BrentSolver brentSolver = new BrentSolver(
			KineticParameters.BRENT_RELATIVE_TOLERANCE,
			KineticParameters.BRENT_ABSOLUTE_TOLERANCE);
		try{
			return brentSolver.solve(KineticParameters.BRENT_MAX_EVALUATIONS, objectiveFunction, lowerBound, upperBound);
		}
		catch(final NoBracketingException exception){
			System.err.printf(Locale.US,
				"BrentSolver failed: no sign change in [%.2e, %.2e]: %s%n",
				lowerBound, upperBound, exception.getMessage());
			return upperBound;
		}
	}

	public double getSubstrateMaxBiomassCapacity(){
		return substrateLimitedBiomassCapacity * 0.9999;
	}

	/**
	 * Fermentation progress function
	 * <p>
	 * progress(N0) = N(t_final) / N_MAX_substrate ∈ [0, 1]
	 * <p>
	 * Monotone increasing in initialAnhydrousYeastFraction for values in (0, substrateLimitedBiomassCapacity):
	 * 	- More initial yeast → faster approach to substrate ceiling → higher progress.
	 * 	- Using N_MAX_substrate (not a fixed broth value) ensures the bracket is valid and results are physically
	 * 		realistic.
	 */
	public double computeFermentationProgress(final double initialAnhydrousYeastFraction){
		final double[] finalState = odeSolver.simulate(yeast.initialPhysiologicalState, initialAnhydrousYeastFraction);
		final double finalTotalSugar = finalState[KineticParameters.IDX_SUGAR_CONCENTRATION]
			+ finalState[KineticParameters.IDX_SUGAR_ENZYMATIC];
		final double dynamicCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * finalTotalSugar;
		return finalState[KineticParameters.IDX_BIOMASS_DENSITY] / dynamicCapacity;
	}

	/**
	 * Estimated lag time in dough [h] — Baranyi & Roberts (1994), eq.(4).
	 * <p>
	 *  λ_dough = ln(1 + 1/Q0_effective) / MU_MAX_REF
	 * <p>
	 *  Q0_effective already incorporates any pre-dough rehydration soak
	 *  (see constructor). For the total wait from the start of rehydration:
	 *    λ_total = rehydrationDurationHours + λ_dough
	 *  Rationale: dQ/dt = MU_MAX_REF · Q (constant, environment-independent).
	 *  Q exits the lag state — i.e., Q/(Q+1) → 1 — after time
	 *    λ = ln(1 + 1/Q0) / MU_MAX_REF
	 *  regardless of temperature or water activity, because those factors only
	 *  affect the BIOMASS growth rate dN/dt, not the physiological recovery dQ/dt.
	 * <p>
	 *  The denominator is MU_MAX_REF, not μ_eff: dQ/dt = MU_MAX_REF · Q is
	 *  environment-independent (Baranyi 1994, eq. 2). Temperature and aw only
	 *  affect dN/dt, not dQ/dt, so the lag exit time is temperature-independent.
	 *  λ_dough is therefore a lower bound on visible leavening; the baker will
	 *  observe significant dough rise somewhat later, when N has grown enough.
	 * <p>
	 *  An infinite value is returned when stage-1 temperature is outside the viable range [T_MIN, T_MAX].
	 */
	public double estimatedLagTime(){
		final double gammaTemp = CardinalTemperatureModel.factor(schedule.getStageTemperature(0),
			KineticParameters.TEMPERATURE_CARDINAL_MIN,
			KineticParameters.TEMPERATURE_CARDINAL_OPT,
			KineticParameters.TEMPERATURE_CARDINAL_MAX);
		return (gammaTemp <= 0.
			? Double.POSITIVE_INFINITY
			: PhysiologicalState.lagPhaseHours(yeast.initialPhysiologicalState));
	}

	/**
	 * Check if the schedule allows sufficient fermentation time.
	 * <p>
	 * Returns a warning string when the total estimated wait (rehydration plus lag in dough) exceeds 50% of the total
	 * fermentation schedule, or null.
	 * Returning the message instead of printing it lets callers embed it inline in tables — avoids interleaved
	 * stderr/stdout noise when called in loops.
	 */
	public String scheduleWarning(){
		final double lagDoughHours = estimatedLagTime();
		final double totalDoughHours = schedule.getTotalDuration();

		if(Double.isInfinite(lagDoughHours))
			return String.format(
				"stage-1 temperature (%.1f °C) is outside the viable range [%.1f, %.1f °C].",
				schedule.getStageTemperature(0),
				KineticParameters.TEMPERATURE_CARDINAL_MIN,
				KineticParameters.TEMPERATURE_CARDINAL_MAX);

		if(lagDoughHours > 0.5 * totalDoughHours)
			return String.format(
				"total wait (lag %.2f h = %.2f h) exceeds 50%% of dough time (%.2f h) "
					+ "for w_b=%.2f — consider a longer schedule or more pre-hydrated yeast.",
				yeast.rehydrationDuration * 60., lagDoughHours, totalDoughHours,
				yeast.moisture);

		return null;
	}

	/**
	 * Print detailed diagnostics when the substrate ceiling is reached.
	 */
	private void printSubstrateCeilingDiagnostics(final double upperBound, final double progressAtUpperBound){
		final double requiredScheduleExtension = (KineticParameters.TARGET_FERMENTATION_PROGRESS - progressAtUpperBound)
			* schedule.getTotalDuration() / progressAtUpperBound;

		System.err.printf(Locale.US,
			"%n" + "═".repeat(100) + "%n"
				+ "⚠️  SUBSTRATE-LIMITED CEILING REACHED%n"
				+ "═".repeat(100) + "%n%n"
				+ "SUMMARY:%n"
				+ "  Maximum yeast possible: %.2f%% WB%n"
				+ "  Fermentation achieved:  %.0f%%%n"
				+ "  Target fermentation:    %.0f%%%n"
				+ "ROOT CAUSE:%n"
				+ "  The limiting factor is SUBSTRATE (fermentable sugar), not yeast quantity.%n%n"
				+ "TO ACHIEVE TARGET:%n"
				+ "  - OPTION 1 (Recommended): Increase sugar%n"
				+ "  - OPTION 2: Extend fermentation schedule by ~%.1f hours%n"
				+ "RETURNING: Substrate ceiling (%.2f%% WB) as conservative estimate.%n"
				+ "═".repeat(100) + "%n%n",
			toBakersWetPercent(upperBound),
			progressAtUpperBound * 100.,
			KineticParameters.TARGET_FERMENTATION_PROGRESS * 100.,
			requiredScheduleExtension,
			toBakersWetPercent(upperBound));
	}


	// ─────────────────────────────────────────────────────────────────────
	// Conversion helpers
	// ─────────────────────────────────────────────────────────────────────

	/**
	 * Converts an anhydrous mass fraction to baker's percentage — dry-yeast basis.
	 * Result: g_dry_yeast / 100 g_flour.
	 */
	private double toBakersDryPercent(final double anhydrousMassFraction){
		final double flourFraction = dough.getFlourFraction();
		if(flourFraction <= 0.)
			throw new ArithmeticException("Flour fraction is zero or negative.");

		return anhydrousMassFraction / flourFraction * 100.;
	}

	/**
	 * Converts an anhydrous mass fraction to baker's percentage — wet-yeast basis.
	 * Result: g_wet_yeast / 100 g_flour.
	 */
	public double toBakersWetPercent(final double anhydrousMassFraction){
		return toBakersDryPercent(yeast.toWetYeastFraction(anhydrousMassFraction));
	}

}