package io.github.mtrevisan.yeastcalculator.working;

import org.apache.commons.math3.ode.events.EventHandler;


final class FoldEventHandler implements EventHandler{

	private final double[] folds;


	FoldEventHandler(final double[] folds){
		this.folds = folds;
	}


	@Override
	public void init(double t0, double[] y0, double t){}

	@Override
	public double g(double t, double[] y){
		// Triggers event boundaries when time 't' crosses a targeted fold coordinate.
		double minDistance = Double.MAX_VALUE;
		for(final double fold : folds){
			final double dist = t - fold;
			if(Math.abs(dist) < Math.abs(minDistance)){
				minDistance = dist;
			}
		}
		return minDistance;
	}

	@Override
	public Action eventOccurred(double t, double[] y, boolean increasing){
		// Structural discontinuity detected: state modification needed
		return Action.RESET_STATE;
	}

	@Override
	public void resetState(double t, double[] y){
		// Physical Degassing: Gas bubbles are mechanical crushed out of existence (V -> 1.0)
		final double vGenerated = y[0];
		y[0] = 1.0 + (vGenerated - 1.0) * 0.15; // 85% of pocket gas is expelled

		// Dissolved chemical CO2 in the water is unaffected by physical folding mechanics
		// y[3] remains unchanged
	}

}
