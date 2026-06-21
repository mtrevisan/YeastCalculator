package io.github.mtrevisan.yeastcalculator.backery;


public enum BakeryProduct{
	NEAPOLITAN_PIZZA(1.7, 0.009, 25.0, 30.0),
	ROMAN_PAN_PIZZA(2.2, 0.013, 15.0, 45.0),
	GASTRONOMY_PAN_PIZZA(2.0, 0.015, 20.0, 50.0),
	BREAD(1.9, 0.012, 20.0, 40.0);


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
