package io.github.mtrevisan.yeastcalculator.simulation;

import org.apache.commons.math3.ode.events.EventHandler;

import java.util.Arrays;


/**
 * Handles instant structural dough adaptations caused by manual stretch & fold techniques.
 * <p>
 * Fold timestamps are expressed in minutes from mixing, the same time axis used by the ODE integration.
 * </p>
 */
public class FoldEventHandler implements EventHandler{

	// Fraction of the retained gas that survives a fold (partial degassing)
	public static final double DEGASSING_RETENTION = 0.85;
	// Multiplicative gain of the gluten gas-holding potential after each fold
	public static final double STRUCTURAL_GAIN = 1.05;

	private final double[] foldTimestamps;
	// Number of folds already applied; the handler is shared across stage integrators, so this survives stage changes
	private int appliedFolds;
	private double currentGasPotentialModifier = 1.;


	public FoldEventHandler(final double[] foldTimestamps){
		this.foldTimestamps = foldTimestamps;
		Arrays.sort(this.foldTimestamps);
	}


	public double getCurrentGasPotentialModifier(){
		return currentGasPotentialModifier;
	}

	@Override
	public void init(final double t0, final double[] y0, final double t){}

	/**
	 * Switching function that crosses zero only at the next pending fold.
	 * <p>
	 * The sign is flipped after every applied fold so that g stays continuous across the reset (it is positive right
	 * after an odd number of folds and negative right after an even one). A "nearest fold" distance would instead jump
	 * from +d to -d halfway between two folds, which the integrator detects as a spurious extra fold.
	 * </p>
	 */
	@Override
	public double g(final double t, final double[] y){
		final double sign = (appliedFolds % 2 == 0? 1.: -1.);
		if(appliedFolds >= foldTimestamps.length)
			// No more folds: keep the sign reached after the last one
			return -sign;
		return sign * (t - foldTimestamps[appliedFolds]);
	}

	@Override
	public Action eventOccurred(final double t, final double[] y, final boolean increasing){
		return Action.RESET_STATE;
	}

	@Override
	public void resetState(final double t, final double[] y){
		// Partial stable gas degassing
		y[2] = y[2] * DEGASSING_RETENTION;
		// Elastic structural gain
		currentGasPotentialModifier *= STRUCTURAL_GAIN;
		appliedFolds ++;
	}

}
