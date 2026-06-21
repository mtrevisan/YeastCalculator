package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Contains ingredient addition ratios computed relative to the total wet-basis flour mass.
 */
public class DoughRecipe{

	private final double waterRatio;

	private final double saltRatio;

	private final double oilRatio;

	private final double maltRatio;
	private final double maltSugarContent;
	private final double maltPollakUnit;
	private final double mixerFrictionFactor;


	public DoughRecipe(final double waterRatio, final double saltRatio, final double oilRatio,
			final double maltRatio, final double maltSugarContent, final double maltPollakUnit,
			final double mixerFrictionFactor){
		this.waterRatio = waterRatio;
		this.saltRatio = saltRatio;
		this.oilRatio = oilRatio;
		this.maltRatio = maltRatio;
		this.maltSugarContent = maltSugarContent;
		this.maltPollakUnit = maltPollakUnit;
		this.mixerFrictionFactor = mixerFrictionFactor;
	}


	public double getWaterRatio(){
		return waterRatio;
	}

	public double getSaltRatio(){
		return saltRatio;
	}

	public double getOilRatio(){
		return oilRatio;
	}

	public double getMaltRatio(){
		return maltRatio;
	}

	public double getMaltSugarContent(){
		return maltSugarContent;
	}

	public double getMaltPollakUnit(){
		return maltPollakUnit;
	}

	public double getMixerFrictionFactor(){
		return mixerFrictionFactor;
	}

}
