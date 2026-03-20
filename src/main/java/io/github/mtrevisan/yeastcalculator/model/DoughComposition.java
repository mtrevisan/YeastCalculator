package io.github.mtrevisan.yeastcalculator.model;


/** Immutable dough composition (mass fractions on total dough). */
public class DoughComposition{

	/** Water fraction [g_water / g_dough] */
	public final double water;
	/** Salt fraction [g_salt / g_dough] */
	public final double salt;
	/**
	 * Effective fermentable sugar fraction [g / g_dough].
	 * = recipe sugar + DMP endogenous sugar + DMP amylase-generated sugar.
	 */
	public final double sugar;
	/** extra ingredients, such as oil, malt, etc. */
	public final double extra;

	/**
	 * Free water content of the flour [g_water / g_wet_flour].
	 * = waterFraction − boundWaterFraction(GAB).
	 * Used in the Ross aw calculation instead of the total waterFraction,
	 * because bound water on flour is not osmotically active.
	 */
	public final double flourFreeWater;


	public DoughComposition(final double water, final double salt, final double sugar, final double extra,
			final double flourFreeWater){
		final double totalDoughWeight = 1. + water + salt + sugar;
		this.water = water / totalDoughWeight;
		this.salt = salt / totalDoughWeight;
		this.sugar = sugar / totalDoughWeight;
		this.extra = extra;

		this.flourFreeWater = flourFreeWater / totalDoughWeight;
	}

	public double getFlourFraction(){
		return 1. - water - salt - sugar - extra;
	}

	public double getTotalWater(){
		return water + flourFreeWater;
	}

}