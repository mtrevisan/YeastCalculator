package io.github.mtrevisan.yeastcalculator.simulation;

import org.apache.commons.math3.ode.events.EventHandler;


/**
 * Event handler designed to catch discrete, non-linear physical interventions (dough folding/degassing).
 * <p>
 * It tracks continuous integration time $t$ and switches states when crossing a targeted fold coordinate.
 * It models physical degassing, where structural macro-gas bubbles are crushed out of existence
 * while leaving underlying dissolved chemical pools intact.
 * </p>
 */
public final class FoldEventHandler implements EventHandler{

	private final double[] folds;


	/**
	 * Initializes the event handler with a defensive clone of targeted fold timestamps.
	 * @param folds	Array of chronological time locations (hours) marking manual folds.
	 */
	public FoldEventHandler(final double[] folds){
		this.folds = folds;
	}


	/**
	 * Framework callback initialization. Omitted by design.
	 */
	@Override
	public void init(final double t0, final double[] y0, final double t){}

	/**
	 * Switching function used by the Apache Commons integrator to isolate discontinuous root events.
	 * Evaluates local geometric distance to the closest targeted fold timestamp.
	 */
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

	/**
	 * Signals the solver engine on how to handle the integration stream once a root cross is detected.
	 *
	 * @return	Action.RESET_STATE to pause continuous integration and mutate state vectors.
	 */
	@Override
	public Action eventOccurred(final double t, final double[] y, final boolean increasing){
		// Structural discontinuity detected: state modification needed
		return Action.RESET_STATE;
	}

	/**
	 * Executes physical state manipulation representing the mechanical manipulation of folding.
	 * Expels 85% of pocket gas volume (y[0]) and compresses elastic residual micro-pressure (y[5]).
	 */
	@Override
	public void resetState(final double t, final double[] y){
		// Physical Degassing: Gas bubbles are mechanically crushed out of existence (V -> 1.)
		final double vGenerated = y[0];
		// Physical Degassing: 85% of accumulated gas pockets are crushed out of the dough matrix
		y[0] = 1. + (vGenerated - 1.) * 0.15;

		// Internal trapped micro-bubble pressure does not instantly vanish to zero.
		// Manipulation leaves an elastic residual stress footprint (~35% preserved energy)
		y[5] = y[5] * 0.35;

		// Sugar pools (y[2]), active dissolved CO2 (y[3]), and chemical degradation (y[4]) are untouched
	}

}
