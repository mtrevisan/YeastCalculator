package io.github.mtrevisan.yeastcalculator.simulation;

import org.apache.commons.math3.ode.events.EventHandler;


/**
 * Handles instant structural dough adaptations caused by manual stretch & fold techniques.
 */
public class FoldEventHandler implements EventHandler{

	private final double[] foldTimestamps;
	private double currentGasPotentialModifier = 1.;


	public FoldEventHandler(final double[] foldTimestamps){
		this.foldTimestamps = foldTimestamps;
	}


	public double getCurrentGasPotentialModifier(){
		return currentGasPotentialModifier;
	}

	@Override
	public void init(final double t0, final double[] y0, final double t){}

	@Override
	public double g(final double t, final double[] y){
		double minDistance = Double.MAX_VALUE;
		for(final double foldTime : foldTimestamps){
			final double distance = t - foldTime;
			if(Math.abs(distance) < Math.abs(minDistance))
				minDistance = distance;
		}
		return minDistance;
	}

	@Override
	public Action eventOccurred(final double t, final double[] y, final boolean increasing){
		return Action.RESET_DERIVATIVES;
	}

	@Override
	public void resetState(final double t, final double[] y){
		// Partial stable gas degassing
		y[2] = y[2] * 0.85;
		// Elastic structural gain
		this.currentGasPotentialModifier *= 1.05;
	}

}
