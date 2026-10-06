package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Contains ingredient addition ratios computed relative to the total wet-basis flour mass.
 */
public record DoughRecipe(
	double waterRatio,
	double saltRatio,
	double oilRatio,
	double maltRatio,
	double maltSugarContent,
	double maltPollakUnit,
	double mixerFrictionFactor){}
