package io.github.mtrevisan.yeastcalculator.backery;


/**
 * Enumeration representing distinct flour grains with their GAB isotherm baseline properties
 * and lookup indices mirroring the original Excel matrix model.
 */
public enum FlourType{
	WHEAT("wheat", 0.062, 11.0, 0.89),
	WHEAT_SEMOLINA("wheat semolina", 0.060, 10.5, 0.90),
	WHEAT_SEMOLINA_FINE("wheat semolina fine", 0.061, 10.8, 0.90),
	DURUM_WHEAT("durum wheat", 0.063, 12.5, 0.88),
	DURUM_WHEAT_SEMOLINA("durum wheat semolina", 0.062, 12.0, 0.89),
	DURUM_WHEAT_SEMOLINA_FINE("durum wheat semolina fine", 0.062, 12.2, 0.89),
	RYE("rye", 0.071, 13.5, 0.87),
	EINKORN("einkorn", 0.064, 11.5, 0.89),
	EMMER("emmer", 0.065, 11.2, 0.89),
	SPELT("spelt", 0.063, 11.0, 0.89),
	BUCKWHEAT("buckwheat", 0.055, 9.0, 0.91),
	BARLEY("barley", 0.068, 10.5, 0.88),
	CHESTNUT("chestnut", 0.050, 8.0, 0.92);

	private final String label;
	private final double wmBase;
	private final double cGab;
	private final double kGab;


	FlourType(String label, double wmBase, double cGab, double kGab){
		this.label = label;
		this.wmBase = wmBase;
		this.cGab = cGab;
		this.kGab = kGab;
	}


	public String getLabel(){
		return label;
	}

	public double getWmBase(){
		return wmBase;
	}

	public double getcGab(){
		return cGab;
	}

	public double getkGab(){
		return kGab;
	}

}
