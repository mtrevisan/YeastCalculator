package io.github.mtrevisan.yeastcalculator.output;

import io.github.mtrevisan.yeastcalculator.kinetics.CardinalTemperatureModel;
import io.github.mtrevisan.yeastcalculator.kinetics.KineticParameters;
import io.github.mtrevisan.yeastcalculator.kinetics.WaterActivityModel;
import io.github.mtrevisan.yeastcalculator.model.DoughComposition;
import io.github.mtrevisan.yeastcalculator.model.FermentationSchedule;
import io.github.mtrevisan.yeastcalculator.model.YeastProperties;
import io.github.mtrevisan.yeastcalculator.ode.BaranyiRobertsODE;
import io.github.mtrevisan.yeastcalculator.ode.NonNegativeClampHandler;
import io.github.mtrevisan.yeastcalculator.physiology.EnvironmentalFactors;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import java.util.Locale;


/**
 * Diagnostic trace: time-series of key state variables during fermentation.
 */
public class SimulationTracer{

	private final DoughComposition dough;
	private final FermentationSchedule schedule;
	private final YeastProperties yeast;
	private final EnvironmentalFactors factors;


	public SimulationTracer(final DoughComposition dough, final FermentationSchedule schedule,
			final YeastProperties yeast){
		this.dough = dough;
		this.schedule = schedule;
		this.yeast = yeast;
		this.factors = new EnvironmentalFactors(dough);
	}

	/** Print time-series trace with 20 points per stage. */
	public void printTrace(final double initialAnhydrousYeastFraction){
		final int PRINT_POINTS_PER_STAGE = 20;

		final double[] state = new double[KineticParameters.STATE_DIMENSION];
		state[KineticParameters.IDX_PHYSIOLOGICAL_STATE] = yeast.initialPhysiologicalState;
		state[KineticParameters.IDX_BIOMASS_DENSITY] = initialAnhydrousYeastFraction;
		state[KineticParameters.IDX_SUGAR_CONCENTRATION] = dough.sugar;
		state[KineticParameters.IDX_CO2_PRODUCED] = 0.;
		state[KineticParameters.IDX_ETHANOL_PRODUCED] = 0.;
		state[KineticParameters.IDX_SUGAR_ENZYMATIC] = 0.;

		printHeader();

		double stageStartTime = 0.;
		for(int stage = 0; stage < schedule.getStageCount(); stage ++){
			final double stageDuration = schedule.getStageDuration(stage);
			if(stageDuration <= 0.)
				continue;

			final double printInterval = stageDuration / PRINT_POINTS_PER_STAGE;

			for(int printPoint = 0; printPoint < PRINT_POINTS_PER_STAGE; printPoint ++){
				final double currentTime = stageStartTime + printPoint * printInterval;
				final double nextTime = Math.min(currentTime + printInterval, stageStartTime + stageDuration);

				printStateRow(currentTime, stage, state);

				// Advance integrator to next print point
				final DormandPrince853Integrator integrator = new DormandPrince853Integrator(
					//minimum adaptive step size [h] — prevents stalling
					1e-6,
					//maximum adaptive step size [h] — ODE solver may use fewer steps
					stageDuration,
					KineticParameters.ODE_ABSOLUTE_TOLERANCE,
					KineticParameters.ODE_RELATIVE_TOLERANCE);
				integrator.addStepHandler(new NonNegativeClampHandler());
				integrator.integrate(new BaranyiRobertsODE(factors, schedule.getStageTemperature(stage)),
					currentTime, state,
					nextTime, state);
				for(int i = 0; i < KineticParameters.STATE_DIMENSION; i ++)
					state[i] = Math.max(0., state[i]);
			}
			stageStartTime += stageDuration;
		}

		final EnvironmentalFactors environmentalFactors = new EnvironmentalFactors(dough);
		printFinalRow(environmentalFactors, stageStartTime, state);
	}

	private void printHeader(){
		System.out.printf(Locale.US,
			"%-8s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-8s%n",
			"time[h]", "stage", "Q", "N", "S", "S_enz", "CO2", "EtOH", "aw", "μ_eff");
	}

	private void printStateRow(final double currentTime, final int stage, final double[] state){
		final double physiologicalState = state[KineticParameters.IDX_PHYSIOLOGICAL_STATE];
		final double biomassDensity = state[KineticParameters.IDX_BIOMASS_DENSITY];
		final double sugarConcentration = state[KineticParameters.IDX_SUGAR_CONCENTRATION];
		final double sugarEnzymatic = state[KineticParameters.IDX_SUGAR_ENZYMATIC];
		final double totalSugarAvailable = sugarConcentration + sugarEnzymatic;

		final double waterActivity = WaterActivityModel.compute(dough, totalSugarAvailable);
		final double sugarAvailability = (totalSugarAvailable > KineticParameters.SUGAR_DEPLETION_THRESHOLD? 1.: 0.);

		final double dynamicBiomassCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * totalSugarAvailable;
		final double effectiveGrowthRate = KineticParameters.MU_MAX_REF
			* factors.gammaTemperature(schedule.getStageTemperature(stage))
			* factors.gammaWaterActivity(totalSugarAvailable)
			* factors.gammaEthanol(state[KineticParameters.IDX_ETHANOL_PRODUCED])
			* (physiologicalState / (physiologicalState + 1.))
			* Math.max(0., 1. - biomassDensity / dynamicBiomassCapacity)
			* sugarAvailability;

		System.out.printf(Locale.US,
			"%-8.2f %-6d %-6.2f %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-6.3f %-1.3f%n",
			currentTime,
			stage + 1,
			physiologicalState,
			biomassDensity * 100.,
			sugarConcentration * 100.,
			sugarEnzymatic * 100.,
			state[KineticParameters.IDX_CO2_PRODUCED] * 100.,
			state[KineticParameters.IDX_ETHANOL_PRODUCED] * 100.,
			waterActivity,
			effectiveGrowthRate);
	}

	private void printFinalRow(final EnvironmentalFactors environmentalFactors, final double finalTime,
			final double[] state){
		final double finalPhysiologicalState = state[KineticParameters.IDX_PHYSIOLOGICAL_STATE];
		final double finalBiomassDensity = state[KineticParameters.IDX_BIOMASS_DENSITY];
		final double finalSugarConcentration = state[KineticParameters.IDX_SUGAR_CONCENTRATION];
		final double finalSugarEnzymatic = state[KineticParameters.IDX_SUGAR_ENZYMATIC];
		final double finalTotalSugar = finalSugarConcentration + finalSugarEnzymatic;
		final double finalDynamicCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * finalTotalSugar;
		final double finalWaterActivity = WaterActivityModel.compute(dough, finalTotalSugar);
		final double finalSugarAvailability = (finalTotalSugar > KineticParameters.SUGAR_DEPLETION_THRESHOLD? 1.: 0.);
		final double finalEffectiveGrowthRate = KineticParameters.MU_MAX_REF
			* CardinalTemperatureModel.factor(schedule.getStageTemperature(schedule.getStageCount() - 1),
				KineticParameters.TEMPERATURE_CARDINAL_MIN, KineticParameters.TEMPERATURE_CARDINAL_OPT,
				KineticParameters.TEMPERATURE_CARDINAL_MAX)
			* WaterActivityModel.inhibitionFactor(finalWaterActivity)
			* environmentalFactors.gammaEthanol(state[KineticParameters.IDX_ETHANOL_PRODUCED])
			* (finalPhysiologicalState / (finalPhysiologicalState + 1.))
			* Math.max(0., 1. - finalBiomassDensity / finalDynamicCapacity)
			* finalSugarAvailability;

		System.out.printf(Locale.US,
			"%-8.2f %-6s %-6.2f %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-4.2f%%  %-6.3f %-1.3f%n",
			finalTime,
			"FINAL",
			state[KineticParameters.IDX_PHYSIOLOGICAL_STATE],
			state[KineticParameters.IDX_BIOMASS_DENSITY] * 100.,
			state[KineticParameters.IDX_SUGAR_CONCENTRATION] * 100.,
			state[KineticParameters.IDX_SUGAR_ENZYMATIC] * 100.,
			state[KineticParameters.IDX_CO2_PRODUCED] * 100.,
			state[KineticParameters.IDX_ETHANOL_PRODUCED] * 100.,
			finalWaterActivity,
			finalEffectiveGrowthRate);
	}

}