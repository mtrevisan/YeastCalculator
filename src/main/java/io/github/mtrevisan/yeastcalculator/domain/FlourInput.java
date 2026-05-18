package io.github.mtrevisan.yeastcalculator.domain;


public final class FlourInput{

	private final double strengthW;
	private final double plRatio;
	private final double sugar;
	private final double protein;
	private final double fat;
	private final double fiber;
	private final double ash;
	private final FlourType type;


	public FlourInput(final double strengthW, final double plRatio, final double sugar, final double protein,
			final double fat, final double fiber, final double ash, final FlourType type){
		this.strengthW = strengthW;
		this.plRatio = plRatio;
		this.sugar = sugar;
		this.protein = protein;
		this.fat = fat;
		this.fiber = fiber;
		this.ash = ash;
		this.type = type;
	}


	public double getStrengthW(){
		return strengthW;
	}

	public double getPlRatio(){
		return plRatio;
	}

	public double getSugar(){
		return sugar;
	}

	public double getProtein(){
		return protein;
	}

	public double getFat(){
		return fat;
	}

	public double getFiber(){
		return fiber;
	}

	public double getAsh(){
		return ash;
	}

	FlourType getType(){
		return type;
	}

	// Delegated parameters directly mapped to the safe Enum properties
	double getWmBase(){
		return type.getWmBase();
	}

	double getCGab(){
		return type.getCGab();
	}

	double getKGab(){
		return type.getKGab();
	}

	public double getBaseLookup(){
		return type.getBaseLookup();
	}

}