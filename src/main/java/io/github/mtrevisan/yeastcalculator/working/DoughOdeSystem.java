package io.github.mtrevisan.yeastcalculator.working;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


final class DoughOdeSystem implements FirstOrderDifferentialEquations{

	private final double yDry;
	private final double[][] stages;
	private final double stiffnessIndexBase;
	private final double saltK;
	private final double oilK;
	private final double waterContent;

	// Henry's Law Constant for CO2 in water at 25°C (mol / m^3 * Pa) converted to mass ratios
	private static final double HENRY_CONSTANT_CO2 = 0.034;
	// Critical saturation threshold before gas bubbles physically nucleate
	private static final double CO2_SATURATION_LIMIT = 0.0015;
	// Critical volume expansion ratio where gluten membranes thin out and fail (rupture point)
	private static final double GLUTEN_POROSITY_THRESHOLD = 2.2;


	DoughOdeSystem(final double yDry, final double[][] stages, final double stiffnessIndexBase,
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
		// y[0] = Dough Volume (V)
		// y[1] = Latent enzyme activation (Lag)
		// y[2] = Residual metabolic sugar substrate
		// y[3] = Dissolved CO2 mass fraction
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
		double tCurr = stages[0][0];
		double rhCurr = stages[0][1];
		for(final double[] stage : stages){
			final double stageDuration = stage[2];
			if(t >= stageStart && t <= (stageStart + stageDuration)){
				tCurr = stage[0];
				rhCurr = stage[1];
				break;
			}
			stageStart += stageDuration;
		}

		// 2. Thermodynamic Henry's law adjustment for temperature
		final double temperatureCorrection = Math.exp(2400. * (1. / (tCurr + 273.15) - 1. / 298.15));
		final double dynamicSaturationLimit = CO2_SATURATION_LIMIT * temperatureCorrection * (1. + 0.5 * saltK);

		// 3. Structural Viscoelastic Stiffness (Maxwell-like degradation)
		final double humidityDeficitFactor = (rhCurr < 0.60) ? (0.50 + 0.50 * (rhCurr / 0.60)) : 1.;
		final double stiffnessTarget = stiffnessIndexBase / humidityDeficitFactor;
		final double proteolysisRate = 0.0015 * Math.exp(0.08 * (tCurr - 20.)); // Enzymatic gluten break-down rate

		// True differential relaxation rate towards structural decay
		final double dynamicStiffness = stiffnessTarget * Math.exp(-proteolysisRate * t);

		// 4. Biological Core Execution
		final double alphaBio = YeastFermentationModel.calculateThermalEfficiency(tCurr);
		final double muBio = YeastFermentationModel.calculateBiomassGrowthRate(yDry, alphaBio, lagCurr, sugarCurr, saltK, oilK);

		yDot[1] = alphaBio; // dLag/dt
		yDot[2] = YeastFermentationModel.calculateNetSugarRate(sugarCurr, muBio, yDry, tCurr); // dSugar/dt

		// 5. Advanced Gas Partitioning Mechanics
		// Total gas mass produced by metabolic respiration (Stoichiometry: Glucose -> 2 CO2)
		final double totalCo2ProductionRate = muBio * 0.48;

		// Dissolution kinetics into liquid water phase vs desorption into gas bubble phase
		double gasDesorptionRate;
		if(co2Dissolved > dynamicSaturationLimit){
			// Mass transfer driving force rate proportional to oversaturation
			gasDesorptionRate = 12.5 * (co2Dissolved - dynamicSaturationLimit) * waterContent;
			yDot[3] = totalCo2ProductionRate - gasDesorptionRate; // dCO2_dissolved/dt
		}
		else{
			// All generated gas remains trapped/dissolved in the fluid matrix phase; zero expansion pressure
			yDot[3] = totalCo2ProductionRate;
			gasDesorptionRate = 0.;
		}

		// 6. Volumetric Work vs Viscous Rheological Resistance
		if(gasDesorptionRate <= 0.){
			yDot[0] = 0.; // No gas phase available to drive mechanical volume work
		}
		else{
			// Internal mechanical gas pressure balanced against matrix structural viscosity
			final double internalPressure = (gasDesorptionRate * 0.08206 * (tCurr + 273.15)) / StrictMath.max(0.1, vCurr - 1.);

			// Permeability factor modeling bubble cell-wall thinning and venting
			double matrixPermeabilityK = 1.;
			if(vCurr > GLUTEN_POROSITY_THRESHOLD){
				// Progressive microstructural tearing of cell walls
				final double overstretch = vCurr - GLUTEN_POROSITY_THRESHOLD;
				matrixPermeabilityK = Math.exp(-3.5 * overstretch);
			}

			// Viscoelastic Volumetric Expansion Velocity Equation (dV/dt)
			yDot[0] = (internalPressure / StrictMath.max(10., dynamicStiffness)) * vCurr * matrixPermeabilityK;
		}
	}

}
