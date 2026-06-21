package io.github.mtrevisan.yeastcalculator.backery;


import org.apache.commons.math3.ode.events.EventHandler;


/**
 * Handles instant structural dough adaptations caused by manual stretch & fold techniques.
 */
public class FoldEventHandler implements EventHandler{

	private final double[] foldTimestamps;
	private double currentGasPotentialModifier = 1.0;


	public FoldEventHandler(double[] foldTimestamps){
		this.foldTimestamps = foldTimestamps;
	}


	public double getCurrentGasPotentialModifier(){
		return currentGasPotentialModifier;
	}

	@Override
	public void init(double t0, double[] y0, double t){}

	@Override
	public double g(double t, double[] y){
		double minDistance = Double.MAX_VALUE;
		for(double foldTime : foldTimestamps){
			double distance = t - foldTime;
			if(Math.abs(distance) < Math.abs(minDistance)){
				minDistance = distance;
			}
		}
		return minDistance;
	}

	@Override
	public Action eventOccurred(double t, double[] y, boolean increasing){
		return Action.RESET_DERIVATIVES;
	}

	@Override
	public void resetState(double t, double[] y){
		y[2] = y[2] * 0.85; // Partial stable gas degassing
		this.currentGasPotentialModifier *= 1.05; // Elastic structural gain
	}

}
