package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Encapsulates chemical and rheological specifications of a single flour instance.
 */
public class FlourInput{

	private final double strength;
	private final double plRatio;
	private final double sugar;
	private final double protein;
	private final double fat;
	private final double fiber;
	private final double ash;
	private final FlourType type;


	public FlourInput(final double strength, final double plRatio, final double sugar, final double protein,
			final double fat, final double fiber, final double ash, final FlourType type){
		this.strength = strength;
		this.plRatio = plRatio;
		this.sugar = sugar;
		this.protein = protein;
		this.fat = fat;
		this.fiber = fiber;
		this.ash = ash;
		this.type = type;
	}


	public double getStrength(){
		return strength;
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
