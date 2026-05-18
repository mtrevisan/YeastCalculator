package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;
import java.util.List;


public final class FlourInput{

	private static final double[] WM_BASE_LOOKUP = {0.062, 0.06, 0.061, 0.063, 0.062, 0.062, 0.071, 0.064, 0.065, 0.063, 0.055, 0.068, 0.05};
	private static final double[] C_GAB_LOOKUP = {12.5, 11.8, 12.1, 13.2, 12.6, 12.8, 8.4, 13., 12.4, 11.9, 6.2, 9.5, 4.8};
	private static final double[] K_GAB_LOOKUP = {0.78, 0.79, 0.78, 0.76, 0.77, 0.77, 0.82, 0.75, 0.77, 0.79, 0.85, 0.81, 0.88};
	private static final double[] BASE_LOOKUP = {1., 0.96, 0.98, 0.85, 0.82, 0.85, 0.35, 0.52, 0.68, 0.78, 0.05, 0.25, 0.02};

	// User-defined parameters
	private final double strengthW;
	private final double plRatio;
	private final double sugar;
	private final double protein;
	private final double fat;
	private final double fiber;
	private final double ash;
	private final FlourType type;

	// Auto-resolved botanical constants
	private final double wmBase;
	private final double cGab;
	private final double kGab;
	private final double baseLookup;


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

		// Automatically resolve static properties on instantiation
		final int idx = type.ordinal();
		this.wmBase = WM_BASE_LOOKUP[idx];
		this.cGab = C_GAB_LOOKUP[idx];
		this.kGab = K_GAB_LOOKUP[idx];
		this.baseLookup = BASE_LOOKUP[idx];
	}


	// Dynamic instance accessors combining both worlds natively
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

	FlourType getType(){
		return type;
	}

	double getWmBase(){
		return wmBase;
	}

	double getCGab(){
		return cGab;
	}

	double getKGab(){
		return kGab;
	}

	double getBaseLookup(){
		return baseLookup;
	}

}
