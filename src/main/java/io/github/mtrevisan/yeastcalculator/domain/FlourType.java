package io.github.mtrevisan.yeastcalculator.domain;


public enum FlourType{
	WHEAT("wheat", 0.062, 12.5, 0.78, 1.),
	WHEAT_SEMOLINA("wheat semolina", 0.06, 11.8, 0.79, 0.96),
	WHEAT_SEMOLINA_FINE("wheat semolina fine", 0.061, 12.1, 0.78, 0.98),
	DURUM_WHEAT("durum wheat", 0.063, 13.2, 0.76, 0.85),
	DURUM_WHEAT_SEMOLINA("durum wheat semolina", 0.062, 12.6, 0.77, 0.82),
	DURUM_WHEAT_SEMOLINA_FINE("durum wheat semolina fine", 0.062, 12.8, 0.77, 0.85),
	RYE("rye", 0.071, 8.4, 0.82, 0.35),
	EINKORN("einkorn", 0.064, 13., 0.75, 0.52),
	EMMER("emmer", 0.065, 12.4, 0.77, 0.68),
	SPELT("spelt", 0.063, 11.9, 0.79, 0.78),
	BUCKWHEAT("buckwheat", 0.055, 6.2, 0.85, 0.05),
	BARLEY("barley", 0.068, 9.5, 0.81, 0.25),
	CHESTNUT("chestnut", 0.05, 4.8, 0.88, 0.02);


	private final String label;
	private final double wmBase;
	private final double cGab;
	private final double kGab;
	private final double baseLookup;


	FlourType(final String label, final double wmBase, final double cGab, final double kGab, final double baseLookup){
		this.label = label;
		this.wmBase = wmBase;
		this.cGab = cGab;
		this.kGab = kGab;
		this.baseLookup = baseLookup;
	}


	public String getLabel(){
		return label;
	}

	public double getWmBase(){
		return wmBase;
	}

	public double getCGab(){
		return cGab;
	}

	public double getKGab(){
		return kGab;
	}

	public double getBaseLookup(){
		return baseLookup;
	}

}
