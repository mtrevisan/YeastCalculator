package io.github.mtrevisan.yeastcalculator.ode;

import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;


/**
 * Step handler that clamps state variables to [0, ∞) after each ODE step.
 * <p>
 * Purpose: Prevents the adaptive step-size controller from exploring unphysical regions
 * where biomass, sugar, or other concentrations would be negative.
 * <p>
 * This is a safety mechanism; properly calibrated tolerances should prevent negative
 * values naturally. However, roundoff errors near boundaries (e.g., sugar depleted,
 * biomass at ceiling) can produce small negatives that accumulate.
 */
public class NonNegativeClampHandler implements StepHandler{

	@Override
	public void init(final double t0, final double[] y0, final double t){}

	@Override
	public void handleStep(final StepInterpolator interpolator, final boolean isLast){
		final double[] state = interpolator.getInterpolatedState();
		for(int i = 0; i < state.length; i ++)
			if(state[i] < 0.)
				state[i] = 0.;
	}

}