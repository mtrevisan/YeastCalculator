package io.github.mtrevisan.yeastcalculator.working.simulation;

import io.github.mtrevisan.yeastcalculator.working.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.working.domain.YeastFermentationModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public final class DoughOdeSystem implements FirstOrderDifferentialEquations{

	private final double yDry;
	private final StageInput[] stages;
	private final double stiffnessIndexBase;
	private final double saltK;
	private final double oilK;
	private final double waterContent;

	// Domain Physical Constants
	private static final double CO2_SATURATION_LIMIT = 0.0015;
	private static final double GLUTEN_POROSITY_THRESHOLD = 2.2;
	private static final double CO2_DESORPTION_K = 12.5;
	private static final double GAS_CONSTANT_R = 0.08206;
	private static final double TEMPERATURE_KELVIN_OFFSET = 273.15;
	private static final double STOICHIOMETRIC_CO2_YIELD = 0.48;

	private static final double MIN_STIFFNESS_BOUND = 10.0;
	private static final double PERMEABILITY_DECAY_RATE = -3.5;

	private static final double INITIAL_GAS_POROSITY = 0.05;


	public DoughOdeSystem(final double yDry, final StageInput[] stages, final double stiffnessIndexBase,
		final double saltK, final double oilK, final double waterContent){
		this.yDry = yDry;
		this.stages = stages;
		this.stiffnessIndexBase = stiffnessIndexBase;
		this.saltK = saltK;
		this.oilK = oilK;
		this.waterContent = waterContent;
	}


	@Override
	public int getDimension(){
		return 4;
	}

	@Override
	public void computeDerivatives(final double t, final double[] y, final double[] yDot){
		final double vCurr = y[0];
		final double lagCurr = y[1];
		final double sugarCurr = y[2];
		final double co2Dissolved = y[3];

		// 1. Resolve environmental boundaries
		double stageStart = 0.;
		double tCurr = stages[0].getTemperature();
		double rhCurr = stages[0].getRelativeHumidity();
		for(final StageInput stage : stages){
			final double stageDuration = stage.getDuration();
			if(t >= stageStart && t <= (stageStart + stageDuration)){
				tCurr = stage.getTemperature();
				rhCurr = stage.getRelativeHumidity();
				break;
			}
			stageStart += stageDuration;
		}

		// Clamped temperature evaluation for safety
		final double tClamped = Math.clamp(tCurr, 4.0, 45.0);

		// 2. Henry's law temperature adjustments
		final double temperatureCorrection = Math.exp(2400. * (1. / (tClamped + TEMPERATURE_KELVIN_OFFSET) - 1. / 298.15));
		final double dynamicSaturationLimit = CO2_SATURATION_LIMIT * temperatureCorrection * (1. + 0.5 * saltK);

		// 3. Structural Stiffness relaxation
		final double humidityDeficitFactor = Math.exp(-1.2 * (rhCurr - 1.0));
		final double stiffnessTarget = stiffnessIndexBase * humidityDeficitFactor;
		final double proteolysisRate = 0.0015 * Math.exp(0.08 * (tClamped - 20.));
		final double dynamicStiffness = stiffnessTarget * Math.exp(-proteolysisRate * t);

		// 4. Biological Core Integration
		final double alphaBio = YeastFermentationModel.calculateThermalEfficiency(tClamped);
		final double muBio = YeastFermentationModel.calculateBiomassGrowthRate(yDry, alphaBio, lagCurr, sugarCurr,
			saltK, oilK);

		yDot[1] = alphaBio;
		yDot[2] = YeastFermentationModel.calculateNetSugarRate(sugarCurr, muBio, yDry, tClamped);

		// 5. Gas Partitioning Kinetics
		final double totalCo2ProductionRate = muBio * STOICHIOMETRIC_CO2_YIELD;
		double gasDesorptionRate;
		if(co2Dissolved > dynamicSaturationLimit){
			gasDesorptionRate = CO2_DESORPTION_K * (co2Dissolved - dynamicSaturationLimit) * waterContent;
			yDot[3] = totalCo2ProductionRate - gasDesorptionRate;
		}
		else{
			yDot[3] = totalCo2ProductionRate;
			gasDesorptionRate = 0.;
		}

		// 6. Volumetric Work vs Resistance
		if(gasDesorptionRate <= 0.)
			yDot[0] = 0.;
		else{
			final double netGasVolume = (vCurr - 1.0) + INITIAL_GAS_POROSITY;
			final double internalPressure = (gasDesorptionRate * GAS_CONSTANT_R * (tClamped + TEMPERATURE_KELVIN_OFFSET))
				/ netGasVolume;

			double matrixPermeabilityK = 1.;
			if(vCurr > GLUTEN_POROSITY_THRESHOLD){
				final double overstretch = vCurr - GLUTEN_POROSITY_THRESHOLD;
				matrixPermeabilityK = Math.exp(PERMEABILITY_DECAY_RATE * overstretch);
			}

			yDot[0] = (internalPressure / StrictMath.max(MIN_STIFFNESS_BOUND, dynamicStiffness)) * vCurr
				* matrixPermeabilityK;
		}
	}

}
