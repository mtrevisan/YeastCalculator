package io.github.mtrevisan.yeastcalculator;

import io.github.mtrevisan.yeastcalculator.model.DoughComposition;
import io.github.mtrevisan.yeastcalculator.model.FermentationSchedule;
import io.github.mtrevisan.yeastcalculator.model.YeastProperties;
import io.github.mtrevisan.yeastcalculator.optimization.YeastOptimizer;
import io.github.mtrevisan.yeastcalculator.output.SensitivityAnalyzer;
import io.github.mtrevisan.yeastcalculator.output.SimulationTracer;


/**
 * <h2>Overview</h2>
 * Predictive microbiology model for <i>Saccharomyces cerevisiae</i> in bread dough.
 * Computes the optimal initial yeast inoculum (g yeast / g dough) so that fermentation reaches 99% of the
 * substrate-limited maximum exactly at the end of the last fermentation stage.
 *
 * <h2>What this model computes</h2>
 * Given a bread dough composition and a fermentation schedule, finds the minimum quantity of <em>anhydrous</em> (dry)
 * yeast [g_dry_yeast / g_dough] such that biomass reaches 99% of the substrate-limited maximum exactly at the end of
 * the last fermentation stage.
 * <p>
 * The result is expressed as <strong>anhydrous</strong> yeast. To convert to the actual wet-yeast quantity to add,
 * divide by {@code (1 − yeastMoistureFraction)}:
 * <pre>
 *   wetYeastFraction = anhydrousYeastFraction / (1 − yeastMoistureFraction)
 * </pre>
 *
 * <h2>ODE system — Baranyi-Roberts (1994)</h2>
 * State vector: [physiologicalState, biomassDensity, sugarConcentration, co2Produced, ethanolProduced, sugarEnzymatic]
 * <pre>
 *   dQ/dt = MU_MAX · Q   (Baranyi (1994), eq.2)
 *   dN/dt = MU_MAX · γ_T · γ_aw · γ_E · [Q/(Q+1)] · [1 - N/N_max] · N
 *   dS/dt = −(1/Y_XS) · dN/dt
 *   dP/dt = Y_PS · (−dS/dt)
 *   dE/dt = Y_ES · (−dS/dt)
 *   dS_enz/dt = k_amylase · (1 − S_enz/S_max) · γ_T_amylase
 * </pre>
 *
 * <h2>Environmental factors</h2>
 * <ul>
 *   <li>γ_T  — Cardinal Model, Rosso et al. (1995); dough cardinal points from Hamad et al. (2012)</li>
 *   <li>γ_aw — linear ramp on [AW_MIN, AW_OPT_DOUGH = 0.975]; NOT the broth value 0.995</li>
 *   <li>γ_E  — linear ethanol inhibition, Nagodawithana et al. (1986)</li>
 * </ul>
 *
 * <h2>Q0 from yeast moisture content</h2>
 * Q0 (Baranyi physiological state) is not a property of a "yeast type" label — it reflects how hydrated and
 * metabolically active the cells are at the moment they enter the dough. This model derives Q0 from the yeast moisture
 * content
 * {@code yeastMoistureFraction} (g_water / g_wet_yeast):
 * <ul>
 *   <li>Anhydrous (w_b ≤ 0.05): Q0 = Q0_ANHYDROUS = 0.02 (cells in near-anabiosis)</li>
 *   <li>Fresh compressed (w_b ≈ 0.70): Q0 = Q0_FRESH_COMPRESSED = 2.50</li>
 *   <li>Intermediate: power-law interpolation (calibrated from rehydration kinetics)</li>
 * </ul>
 *
 * <h2>Literature</h2>
 * <ul>
 *   <li>Baranyi, J., & Roberts, T. A. (1994). A dynamic approach to predicting bacterial growth in food. International Journal of Food Microbiology, 23(3–4), 277–294.</li>
 *   <li>Rosso L. et al. (1995) — J Theor Biol 175:159-171</li>
 *   <li>Ross T. (1993) — J Appl Bacteriol 74:608-611</li>
 *   <li>Hamad S. et al. (2012) — Food Microbiol 32:228-236</li>
 *   <li>Nagodawithana T.W. et al. (1986) — Can J Microbiol 32:467-472</li>
 *   <li>Gervais P., Beney L. (2001) — J Appl Microbiol 90:579-585</li>
 * </ul>
 *
 * <h2>Units</h2>
 * All mass fractions are dimensionless (g / g_dough).
 * Time in hours.
 * Temperature in Celsius.
 */
public class YeastCalculator{

	private final DoughComposition dough;
	private final FermentationSchedule schedule;
	private final YeastProperties yeast;
	private final YeastOptimizer yeastOptimizer;


	/**
	 * @param water	Water [baker's percentage].
	 * @param salt	Salt [baker's percentage].
	 * @param sugar	Recipe fermentable sugar [baker's percentage]
	 * 	(endogenous flour sugars + added sugar, excluding malt).
	 * @param stageTemperatures	Stage temperatures [°C].
	 * @param stageDurations	Stage durations [h].
	 * @param yeastMoisture	Moisture of the yeast [g_water / g_wet_yeast]:
	 * 	<ul>
	 * 		<li>≈ 0.04–0.06 — instant/rapid dry yeast (IDY/SDY)</li>
	 * 		<li>≈ 0.07–0.09 — active dry yeast (ADY)</li>
	 * 		<li>≈ 0.68–0.72 — fresh compressed yeast (CY)</li>
	 * 	</ul>
	 * @param rehydrationDuration	Time the yeast soaks in warm water before being added to the dough [h]. Use 0
	 * 	for yeast added directly to flour without pre-activation. During rehydration Q grows as
	 * 	Q(t) = Q0_dry × exp(MU_MAX_REF × t), so even a short soak significantly raises Q0 for very dry yeast.
	 */
	public YeastCalculator(final double water, final double salt, final double sugar, final double extra,
			final double[] stageTemperatures, final double[] stageDurations,
			final double flourFreeWater, final double yeastMoisture,
			final double rehydrationDuration){
		dough = new DoughComposition(water, salt, sugar, extra, flourFreeWater);
		schedule = new  FermentationSchedule(stageTemperatures, stageDurations);
		yeast = new  YeastProperties(yeastMoisture, rehydrationDuration);
		yeastOptimizer = new YeastOptimizer(dough, schedule, yeast);
	}

	double findOptimal(final boolean outputWarning){
		return yeastOptimizer.findOptimal(outputWarning);
	}

	int estimatedLagTime(){
		return (int)Math.ceil(yeastOptimizer.estimatedLagTime() * 60);
	}

	double yeastWB(final double yeastDB){
		return yeastOptimizer.toBakersWetPercent(yeastDB);
	}

	boolean hasWarning(){
		return (yeastOptimizer.scheduleWarning() != null);
	}

	SimulationTracer getSimulationTracer(){
		return new SimulationTracer(dough, schedule, yeast);
	}

	SensitivityAnalyzer getSensitivityAnalyzer(){
		return new SensitivityAnalyzer(yeastOptimizer);
	}

}
