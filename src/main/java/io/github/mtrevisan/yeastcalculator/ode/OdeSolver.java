package io.github.mtrevisan.yeastcalculator.ode;

import io.github.mtrevisan.yeastcalculator.kinetics.KineticParameters;
import io.github.mtrevisan.yeastcalculator.model.DoughComposition;
import io.github.mtrevisan.yeastcalculator.model.FermentationSchedule;
import io.github.mtrevisan.yeastcalculator.physiology.EnvironmentalFactors;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;


/** ODE integration across fermentation stages. */
public class OdeSolver{

	private final DoughComposition dough;
	private final FermentationSchedule schedule;
	private final EnvironmentalFactors factors;


	public OdeSolver(final DoughComposition dough, final FermentationSchedule schedule){
		this.dough = dough;
		this.schedule = schedule;
		this.factors = new EnvironmentalFactors(dough);
	}

	/**
	 * Integrate the ODE across all fermentation stages.
	 *
	 * @param initialAnhydrousYeastFraction	[g_dry_yeast / g_dough]
	 * @return	Final ODE state vector [Q, N, S, P_CO2, E_EtOH].
	 */
	public double[] simulate(final double initialPhysiologicalState, final double initialAnhydrousYeastFraction){
		final double[] state = new double[KineticParameters.STATE_DIMENSION];
		state[KineticParameters.IDX_PHYSIOLOGICAL_STATE] = initialPhysiologicalState;
		state[KineticParameters.IDX_BIOMASS_DENSITY] = initialAnhydrousYeastFraction;
		state[KineticParameters.IDX_SUGAR_CONCENTRATION] = dough.sugar;
		state[KineticParameters.IDX_CO2_PRODUCED] = 0.;
		state[KineticParameters.IDX_ETHANOL_PRODUCED] = 0.;
		state[KineticParameters.IDX_SUGAR_ENZYMATIC] = 0.;

		double stageStartTime = 0.;
		for(int stage = 0; stage < schedule.getStageCount(); stage ++){
			final double stageDuration = schedule.getStageDuration(stage);
			if(stageDuration <= 0.)
				continue;

			final DormandPrince853Integrator integrator = new DormandPrince853Integrator(
				//minimum adaptive step size [h] — prevents stalling
				1e-6,
				//maximum adaptive step size [h] — ODE solver may use fewer steps
				stageDuration,
				KineticParameters.ODE_ABSOLUTE_TOLERANCE,
				KineticParameters.ODE_RELATIVE_TOLERANCE);
			integrator.addStepHandler(new NonNegativeClampHandler());
			integrator.integrate(new BaranyiRobertsODE(factors, schedule.getStageTemperature(stage)),
				stageStartTime, state,
				stageStartTime + stageDuration, state);

			//hard clamp after each stage (belt-and-suspenders)
			for(int i = 0; i < KineticParameters.STATE_DIMENSION; i ++)
				state[i] = Math.max(0., state[i]);

			stageStartTime += stageDuration;
		}
		return state;
	}

}