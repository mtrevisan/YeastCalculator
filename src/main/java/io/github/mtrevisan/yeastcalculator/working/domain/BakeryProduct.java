package io.github.mtrevisan.yeastcalculator.working.domain;


public enum BakeryProduct{
	/**
	 * NEAPOLITAN_PIZZA:
	 * Low volume in tray to preserve extensibility for manual stretching.
	 * Lower sugar threshold because extreme baking temperatures (450 °C+)
	 * will cause flash charring if residual sugars are too high.
	 */
	NEAPOLITAN_PIZZA(1.7, 0.009, 25.0, 30.),

	/**
	 * PAN_PIZZA_ROMAN:
	 * High hydration, large open alveoli. Maximum volume ceiling.
	 * Medium-high sugar required to sustain longer baking charts at 250 °C.
	 */
	ROMAN_PAN_PIZZA(2.2, 0.013, 15.0, 45.0),

	/**
	 * GASTRONOMY_TEGLIA:
	 * Spongy, soft, high thickness crumb with small, uniform, dense bubble structures.
	 * Controlled volume cap to prevent cell walls from thinning and merging into caves.
	 * High sugar residue target (1.5%) to feed the crumb during very long, gentle bake profiles (220 °C).
	 */
	GASTRONOMY_PAN_PIZZA(2.0, 0.015, 20., 50.),

	/**
	 * BREAD:
	 * Balanced freestanding three-dimensional expansion profile.
	 * Requires structural elasticity reserves for scoring cuts and steam oven spring.
	 */
	BREAD(1.9, 0.012, 20., 40.);


	private final double glutenTearingLimit;
	private final double minSafeSugarThreshold;
	private final double tearingPenaltyMultiplier;
	private final double starvationPenaltyMultiplier;


	BakeryProduct(final double vLimit, final double sLimit, final double pTear, final double pStarve){
		this.glutenTearingLimit = vLimit;
		this.minSafeSugarThreshold = sLimit;
		this.tearingPenaltyMultiplier = pTear;
		this.starvationPenaltyMultiplier = pStarve;
	}


	public double getGlutenTearingLimit(){
		return glutenTearingLimit;
	}

	public double getMinSafeSugarThreshold(){
		return minSafeSugarThreshold;
	}

	public double getTearingPenaltyMultiplier(){
		return tearingPenaltyMultiplier;
	}

	public double getStarvationPenaltyMultiplier(){
		return starvationPenaltyMultiplier;
	}

}
