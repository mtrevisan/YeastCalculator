package io.github.mtrevisan.yeastcalculator.working.simulation;

import org.apache.commons.math3.ode.events.EventHandler;


public final class FoldEventHandler implements EventHandler{

	private final double[] folds;


	public FoldEventHandler(final double[] folds){
		this.folds = folds;
	}


	@Override
	public void init(final double t0, final double[] y0, final double t){}

	@Override
	public double g(final double t, final double[] y){
		// Triggers event boundaries when time 't' crosses a targeted fold coordinate.
		double minDistance = Double.MAX_VALUE;
		for(final double fold : folds){
			final double dist = t - fold;
			if(Math.abs(dist) < Math.abs(minDistance))
				minDistance = dist;
		}
		return minDistance;
	}

	@Override
	public Action eventOccurred(final double t, final double[] y, final boolean increasing){
		// Structural discontinuity detected: state modification needed
		return Action.RESET_STATE;
	}

	@Override
	public void resetState(final double t, final double[] y){
		// Physical Degassing: Gas bubbles are mechanical crushed out of existence (V -> 1.)
		final double vGenerated = y[0];
		// 85% of pocket gas is expelled
		y[0] = 1. + (vGenerated - 1.) * 0.15;

		// Internal trapped micro-bubble pressure does not instantly vanish to zero.
		// Manipulation leaves an elastic residual stress footprint (~35% preserved energy)
		y[5] = y[5] * 0.35;

		// Sugar pools (y[2]), active dissolved CO2 (y[3]), and chemical degradation (y[4]) are untouched
	}

}
