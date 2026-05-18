package io.github.mtrevisan.yeastcalculator.working;


public final class FlourInput{

	private final double strengthW;
	private final double plRatio;
	private final double sugar;
	private final double protein;
	private final double fat;
	private final double fiber;
	private final double ash;
	private final String type;


	public FlourInput(final double strengthW, final double plRatio, final double sugar, final double protein,
			final double fat, final double fiber, final double ash, final String type){
		this.strengthW = strengthW;
		this.plRatio = plRatio;
		this.sugar = sugar;
		this.protein = protein;
		this.fat = fat;
		this.fiber = fiber;
		this.ash = ash;
		this.type = type;
	}


	double getStrengthW(){
		return strengthW;
	}

	double getPlRatio(){
		return plRatio;
	}

	double getSugar(){
		return sugar;
	}

	double getProtein(){
		return protein;
	}

	double getFat(){
		return fat;
	}

	double getFiber(){
		return fiber;
	}

	double getAsh(){
		return ash;
	}

	String getType(){
		return type;
	}

}
