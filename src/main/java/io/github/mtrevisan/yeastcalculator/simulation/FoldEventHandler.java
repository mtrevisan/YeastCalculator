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
	// Fraction of the residual gluten underdevelopment recovered by each fold set
	public static final double DEVELOPMENT_GAIN = 0.2;
	// Gas retention factor of a completely undeveloped network (a fully developed one has factor 1)
	private static final double MIN_GAS_RETENTION_FACTOR = 0.6;

	private final double[] foldTimestamps;
	// Gluten development index reached by kneading alone, before any fold
	private final double initialDevelopment;
	// Number of folds already applied; the handler is shared across stage integrators, so this survives stage changes
	private int appliedFolds;
	private double currentGasPotentialModifier = 1.;


	/**
	 * @param foldTimestamps	Fold times [min].
	 * @param initialDevelopment	Gluten development index reached by kneading alone, in [0, 1].
	 */
	public FoldEventHandler(final double[] foldTimestamps, final double initialDevelopment){
		this.foldTimestamps = foldTimestamps.clone();
		Arrays.sort(this.foldTimestamps);
		this.initialDevelopment = initialDevelopment;
	}


	/**
	 * Gluten development index after a number of fold sets: each set recovers a fixed fraction of the residual
	 * underdevelopment, so folds matter a lot after a gentle hand mix and very little after a full machine mix.
	 */
	public static double developmentAfterFolds(final double initialDevelopment, final int folds){
		return 1. - (1. - initialDevelopment) * Math.pow(1. - DEVELOPMENT_GAIN, folds);
	}

	/**
	 * Fraction of the full gas-holding capacity available to a gluten network with the given development index.
	 */
	public static double gasRetentionFactor(final double developmentIndex){
		return MIN_GAS_RETENTION_FACTOR + (1. - MIN_GAS_RETENTION_FACTOR) * developmentIndex;
	}

	/**
	 * Gain of the gas-holding potential produced by a number of fold sets, relative to the kneaded-only network.
	 */
	public static double foldGasPotentialGain(final double initialDevelopment, final int folds){
		return gasRetentionFactor(developmentAfterFolds(initialDevelopment, folds)) / gasRetentionFactor(initialDevelopment);
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
		// RESET_STATE is required: with RESET_DERIVATIVES the integrator never calls resetState
		return Action.RESET_STATE;
	}

	@Override
	public void resetState(final double t, final double[] y){
		// Partial stable gas degassing
		y[2] = y[2] * DEGASSING_RETENTION;
		// Elastic structural gain from the additional gluten development
		appliedFolds ++;
		currentGasPotentialModifier = foldGasPotentialGain(initialDevelopment, appliedFolds);
	}

}
