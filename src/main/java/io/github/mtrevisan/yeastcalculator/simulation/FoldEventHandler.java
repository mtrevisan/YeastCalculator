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

	// Fraction of the retained gas that survives a fold (partial degassing).
	// A manual stretch-and-fold loses 5-8% of the accumulated gas (Hamelman, Bread, 2nd ed.;
	// Reinhart, The Bread Baker's Apprentice). The previous value of 0.85 (15% loss) was the
	// upper end of the observed range and systematically penalized folds in the optimizer.
	public static final double DEGASSING_RETENTION = 0.92;
	// Multiplicative gain of the gluten gas-holding potential after each fold
	public static final double STRUCTURAL_GAIN = 1.05;

	private final double[] foldTimestamps;
	// Number of folds already applied; the handler is shared across stage integrators, so this survives stage changes
	private int appliedFolds;
	private double currentGasPotentialModifier = 1.;


	public FoldEventHandler(final double[] foldTimestamps){
		this.foldTimestamps = foldTimestamps.clone();
		Arrays.sort(this.foldTimestamps);
	}


	public double getCurrentGasPotentialModifier(){
		return currentGasPotentialModifier;
	}

	public int getAppliedFolds(){
		return appliedFolds;
	}

	/**
	 * Advances the internal potential modifier (and applied fold counter) as if all folds whose timestamp
	 * is less than or equal to {@code timeInMinutes} had already been triggered.
	 * <p>
	 * Does NOT modify the state vector y (no degassing is applied here), which is exactly what we want when
	 * recomputing the final gas retention ceiling for the fitness evaluation at the timeline endpoint.
	 * </p>
	 */
	public void advanceFoldsUpTo(final double timeInMinutes){
		while(appliedFolds < foldTimestamps.length && foldTimestamps[appliedFolds] <= timeInMinutes){
			currentGasPotentialModifier *= STRUCTURAL_GAIN;
			appliedFolds ++;
		}
	}

	@Override
	public void init(final double t0, final double[] y0, final double t){}

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
