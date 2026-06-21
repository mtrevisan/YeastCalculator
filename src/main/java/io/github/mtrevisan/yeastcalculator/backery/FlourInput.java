package io.github.mtrevisan.yeastcalculator.backery;


/**
 * Encapsulates chemical and rheological specifications of a single flour instance.
 */
public class FlourInput{

	private final double strengthW;
	private final double plRatio;
	private final double sugar;
	private final double protein;
	private final double fat;
	private final double fiber;
	private final double ash;
	private final FlourType type;


	public FlourInput(double strengthW, double plRatio, double sugar, double protein,
			double fat, double fiber, double ash, FlourType type){
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

	public FlourType getType(){
		return type;
	}

}
