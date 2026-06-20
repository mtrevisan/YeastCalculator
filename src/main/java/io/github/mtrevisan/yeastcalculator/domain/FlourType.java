package io.github.mtrevisan.yeastcalculator.domain;


public enum FlourType{
	WHEAT("wheat", 0.062, 11., 0.89, 1.),
	WHEAT_SEMOLINA("wheat semolina", 0.06, 10.5, 0.9, 2.),
	WHEAT_SEMOLINA_FINE("wheat semolina fine", 0.061, 10.8, 0.9, 3.),
	DURUM_WHEAT("durum wheat", 0.063, 12.5, 0.88, 4.),
	DURUM_WHEAT_SEMOLINA("durum wheat semolina", 0.062, 12., 0.89, 5.),
	DURUM_WHEAT_SEMOLINA_FINE("durum wheat semolina fine", 0.062, 12.2, 0.89, 6.),
	RYE("rye", 0.071, 13.5, 0.87, 7.),
	EINKORN("einkorn", 0.064, 11.5, 0.89, 8.),
	EMMER("emmer", 0.065, 11.2, 0.89, 9.),
	SPELT("spelt", 0.063, 11., 0.89, 10.),
	BUCKWHEAT("buckwheat", 0.055, 9., 0.91, 11.),
	BARLEY("barley", 0.068, 10.5, 0.88, 12.),
	CHESTNUT("chestnut", 0.05, 8., 0.92, 13.);


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
