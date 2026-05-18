package io.github.mtrevisan.yeastcalculator.working.optimization;

import io.github.mtrevisan.yeastcalculator.working.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.working.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.working.domain.FlourType;
import io.github.mtrevisan.yeastcalculator.working.domain.StageInput;


public class SimulationInputs{

	private final BakeryProduct targetProduct = BakeryProduct.GASTRONOMY_PAN_PIZZA;

	private final double[] fractions = {1.};

	private final FlourInput[] flourMatrix = {
		new FlourInput(295., 0.55, 0.013, 0.13, 0.011, 0.019, 0.003, FlourType.WHEAT)
	};

	private final StageInput[] stages = {
		new StageInput(26., 0.55, 8.)
	};

	private final double[] folds = {};

	private final double doughWater = 0.615;
	private final double doughSalt = 0.022;
	private final double doughOil = 0.039;
	private final double yeastMoisture = 0.70;
	private final double flourTemperature = 19.;
	private final double airRelativeHumidity = 0.54;

	// Rheological empirical constant reference
	private static final double IDEAL_HYDRATION_REFERENCE = 0.6;


	public int getFlourCount(){
		return fractions.length;
	}

	public BakeryProduct getTargetProduct() { return targetProduct; }

	public double[] getFractions(){
		final int flours = fractions.length;
		double sumFractions = 0.;
		for(final double f : fractions)
			sumFractions += f;
		final double[] normalized = new double[flours];
		for(int i = 0; i < flours; i ++)
			normalized[i] = this.fractions[i] / sumFractions;
		return normalized;
	}

	/**
	 * Computes the complete, water-adjusted structural stiffness base index.
	 * Resolves previous architectural Data Envy leakages.
	 */
	public double getHydratedStiffnessIndex(final double flourMoisture){
		double dotStrength = 0.;
		double dotPL = 0.;
		double dotProtein = 0.;
		double dotFiber = 0.;
		double dotAsh = 0.;
		double dotγ = 0.;
		final double[] normalizedFractions = getFractions();
		for(int i = 0; i < fractions.length; i ++){
			final FlourInput flour = flourMatrix[i];
			final double γ_i = flour.getBaseLookup() * Math.exp(-2. * flour.getFat()) * Math.exp(-0.5 * flour.getSugar());
			final double f = normalizedFractions[i];

			dotStrength += flour.getStrengthW() * f;
			dotPL += flour.getPlRatio() * f;
			dotProtein += flour.getProtein() * f;
			dotFiber += flour.getFiber() * f;
			dotAsh += flour.getAsh() * f;
			dotγ += γ_i * f;
		}

		final double rawStiffnessBase = (dotStrength / dotPL) * (dotProtein / (1. + 2. * dotFiber + 5. * dotAsh)) * dotγ;

		// Viscoelastic hydration modifier logic integrated here
		final double waterEff = (doughWater + flourMoisture) / (1. - flourMoisture);
		final double structuralHydrationModifier = 1. + Math.pow(waterEff - IDEAL_HYDRATION_REFERENCE, 2);
		return rawStiffnessBase / structuralHydrationModifier;
	}

	public double getBlendedSugar(){
		double dotSugar = 0.;
		final double[] normalizedFractions = getFractions();
		for(int i = 0; i < fractions.length; i ++)
			dotSugar += flourMatrix[i].getSugar() * normalizedFractions[i];
		return dotSugar;
	}

	public FlourInput[] getFlourMatrix(){
		return flourMatrix;
	}

	public StageInput[] getStages(){
		return stages;
	}

	public double[] getFolds(){
		return folds;
	}

	public double getDoughWater(){
		return doughWater;
	}

	public double getDoughSalt(){
		return doughSalt;
	}

	public double getDoughOil(){
		return doughOil;
	}

	public double getYeastMoisture(){
		return yeastMoisture;
	}

	public double getAirRelativeHumidity(){
		return airRelativeHumidity;
	}

	public double getFlourTemperature(){
		return flourTemperature;
	}

}
