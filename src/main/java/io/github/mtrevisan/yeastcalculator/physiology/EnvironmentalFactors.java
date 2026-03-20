package io.github.mtrevisan.yeastcalculator.physiology;

import io.github.mtrevisan.yeastcalculator.kinetics.CardinalTemperatureModel;
import io.github.mtrevisan.yeastcalculator.kinetics.KineticParameters;
import io.github.mtrevisan.yeastcalculator.kinetics.WaterActivityModel;
import io.github.mtrevisan.yeastcalculator.model.DoughComposition;


/** Compute gamma factors (environmental inhibition) for Baranyi model. */
public class EnvironmentalFactors{

	private final DoughComposition dough;


	public EnvironmentalFactors(final DoughComposition dough){
		this.dough = dough;
	}

	public double gammaTemperature(final double temperature){
		return CardinalTemperatureModel.factor(temperature,
			KineticParameters.TEMPERATURE_CARDINAL_MIN,
			KineticParameters.TEMPERATURE_CARDINAL_OPT,
			KineticParameters.TEMPERATURE_CARDINAL_MAX);
	}

	public double gammaWaterActivity(final double sugarConcentration){
		final double aw = WaterActivityModel.compute(dough, sugarConcentration);
		return WaterActivityModel.inhibitionFactor(aw);
	}

	/** Ethanol inhibition factor — γ_E ∈ [0, 1] — Nagodawithana et al. (1986) */
	public double gammaEthanol(final double ethanolConcentration){
		return Math.max(0., 1. - ethanolConcentration / KineticParameters.ETHANOL_INHIBITION_MAX);
	}

	public double gammaAmylaseTemperature(final double temperature){
		return CardinalTemperatureModel.factor(temperature,
			KineticParameters.T_MIN_AMYLASE,
			KineticParameters.T_OPT_AMYLASE,
			KineticParameters.T_MAX_AMYLASE);
	}

}