package io.github.mtrevisan.yeastcalculator.working;


public enum FlourType{
	WHEAT("wheat"),
	WHEAT_SEMOLINA("wheat semolina"),
	WHEAT_SEMOLINA_FINE("wheat semolina fine"),
	DURUM_WHEAT("durum wheat"),
	DURUM_WHEAT_SEMOLINA("durum wheat semolina"),
	DURUM_WHEAT_SEMOLINA_FINE("durum wheat semolina fine"),
	RYE("rye"),
	EINKORN("einkorn"),
	EMMER("emmer"),
	SPELT("spelt"),
	BUCKWHEAT("buckwheat"),
	BARLEY("barley"),
	CHESTNUT("chestnut");


	private final String label;

	FlourType(final String label){
		this.label = label;
	}

	public String getLabel(){
		return label;
	}

}
