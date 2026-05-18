package io.github.mtrevisan.yeastcalculator.working.simulation;

import io.github.mtrevisan.yeastcalculator.working.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.working.domain.YeastFermentationModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


/**
 * Core physical and bio-chemical integration engine modeling the dough as a
 * 6-Dimensional system of coupled non-linear First Order Differential Equations (ODEs).
 * <p>
 * State Vector Layout ($y$):
 * <ul>
 * <li>y[0] : Macro Volumetric Expansion Ratio ($V = \frac{V_t}{V_0}$)</li>
 * <li>y[1] : Lag Phase Enzymatic Adjustment State Coordinate ($\lambda$)</li>
 * <li>y[2] : Residual Simple Sugar Carbohydrate Substrate ($S$)</li>
 * <li>y[3] : Dissolved chemical $CO_2$ mass fraction pool ($C_{\text{aq}}$)</li>
 * <li>y[4] : Accumulated Macro-Structural Proteolytic Gluten Network Degradation Index ($D$)</li>
 * <li>y[5] : Micro-bubble Internal Volumetric Gas Pressure State ($P$)</li>
 * </ul>
 * </p>
 * <p>
 * This class accounts for non-linear physical interactions, including moisture sineresi,
 * ethanol-driven cell auto-intoxication, and progressive sigmoidal matrix permeability.
 * </p>
 */
public final class DoughOdeSystem implements FirstOrderDifferentialEquations{

	private final double yDry;
	private final StageInput[] stages;
	private final double stiffnessIndexBase;
	private final double saltK;
	private final double oilK;
	private final double waterContent;
	private final double sugarInitial;
	private final double glutenTearingLimit;

	// Domain Physical Constants
	private static final double CO2_SATURATION_LIMIT = 0.0015;
	private static final double CO2_DESORPTION_K = 12.5;
	private static final double GAS_CONSTANT_R = 0.08206;
	private static final double TEMPERATURE_KELVIN_OFFSET = 273.15;
	private static final double STOICHIOMETRIC_CO2_YIELD = 0.48;

	private static final double MIN_STIFFNESS_BOUND = 10.;

	private static final double INITIAL_GAS_POROSITY = 0.05;

	// Critical fermentation shutdown due to ethanol saturation (~6% mass fraction threshold)
	private static final double ETHANOL_CRITICAL_LIMIT = 0.06;


	/**
	 * Configures the multi-dimensional dynamic system variables.
	 */
	public DoughOdeSystem(final double yDry, final StageInput[] stages, final double stiffnessIndexBase,
			final double saltK, final double oilK, final double waterContent, final double sugarInitial,
			final double glutenTearingLimit){
		this.yDry = yDry;
		this.stages = stages;
		this.stiffnessIndexBase = stiffnessIndexBase;
		this.saltK = saltK;
		this.oilK = oilK;
		this.waterContent = waterContent;
		this.sugarInitial = sugarInitial;
		this.glutenTearingLimit = glutenTearingLimit;
	}


	@Override
	public int getDimension(){
		// y[0]=Volume, y[1]=Lag, y[2]=Sugar, y[3]=Dissolved CO2, y[4]=Gluten Degradation, y[5]=Internal Pressure
		return 6;
	}

	/**
	 * Primary derivative computation loop executed by the numerical integrator step handlers.
	 * Solves mass transfers, metabolic pathways, and viscoelastic expansions.
	 */
	@Override
	public void computeDerivatives(final double t, final double[] y, final double[] yDot){
		final double vCurr = y[0];
		final double lagCurr = y[1];
		final double sugarCurr = y[2];
		final double co2Dissolved = y[3];
		final double glutenDegradation = y[4];
		final double internalPressure = y[5];

		// 1. Environmental lookup mapping timeline coordinate 't' to its active stage
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

		final double tClamped = Math.clamp(tCurr, 4., 45.);

		// 2. Thermodynamic Henry's law adjustment for dynamic CO2 fluid saturation limits
		final double temperatureCorrection = Math.exp(2400. * (1. / (tClamped + TEMPERATURE_KELVIN_OFFSET) - 1. / 298.15));
		final double dynamicSaturationLimit = CO2_SATURATION_LIMIT * temperatureCorrection * (1. + 0.5 * saltK);

		// 3. IMPROVEMENT 1: Syneresis modeling. Gluten breakdown releases bound water into the dough matrix
		final double dynamicWaterContent = waterContent + 0.05 * glutenDegradation;

		// 3. REOLOGY MODELLING - Sineresi: Damaged proteins release bound water, fluidizing the matrix
		final double humidityDeficitFactor = Math.exp(-1.2 * (rhCurr - 1.));
		final double stiffnessTarget = stiffnessIndexBase * humidityDeficitFactor;

		final double proteolysisRate = 0.0015 * Math.exp(0.08 * (tClamped - 20.));
		// dDegradation/dt: tracks asymptotic protein structural decay pathway
		yDot[4] = proteolysisRate * (1. - glutenDegradation);

		final double dynamicStiffness = stiffnessTarget * (1. - glutenDegradation);

		// Dynamic structural pH trajectory proxy: organic acid accumulation tracks sugar depletion
		final double sugarConsumed = Math.max(0., sugarInitial - sugarCurr);
		// 4. MICROBIOLOGY - Ethanol growth inhibition feedback loop (Ghose & Tyagi approach)
		final double ethanolInhibitionFactor = Math.clamp(1. - (sugarConsumed / ETHANOL_CRITICAL_LIMIT), 0., 1.);

		// Biological Core Integration with added Ethanol feedback loop
		final double alphaBio = YeastFermentationModel.calculateThermalEfficiency(tClamped);
		final double muBio = YeastFermentationModel.calculateBiomassGrowthRate(yDry, alphaBio, lagCurr, sugarCurr,
			saltK, oilK) * ethanolInhibitionFactor;

		// dLag/dt
		yDot[1] = alphaBio;
		// dSugar/dt
		yDot[2] = YeastFermentationModel.calculateNetSugarRate(sugarCurr, muBio, yDry, tClamped);

		// 5. GAS KINETICS - Separation of dissolved aqueous CO2 vs gaseous pocket phase
		// Pasteur transition effect tracking
		final double anaerobicFactor = 1. - Math.exp(-3. * t);
		final double totalCo2ProductionRate = muBio * STOICHIOMETRIC_CO2_YIELD * anaerobicFactor;
		double gasDesorptionRate;
		if(co2Dissolved > dynamicSaturationLimit){
			gasDesorptionRate = CO2_DESORPTION_K * (co2Dissolved - dynamicSaturationLimit) * dynamicWaterContent;
			yDot[3] = totalCo2ProductionRate - gasDesorptionRate;
		}
		else{
			yDot[3] = totalCo2ProductionRate;
			gasDesorptionRate = 0.;
		}

		// 6. VISCOELASTIC VOLUMETRIC EXPANSION (Maxwell-like pneumatic work balance)
		final double netGasVolume = (vCurr - 1.) + INITIAL_GAS_POROSITY;
		// Ideal equilibrium reference pressure calculation
		final double equilibriumPressure = gasDesorptionRate * GAS_CONSTANT_R * (tClamped + TEMPERATURE_KELVIN_OFFSET)
			/ netGasVolume;

		// Dynamic mass relaxation pathway moving toward spatial pressure equilibrium
		yDot[5] = 25. * (equilibriumPressure - internalPressure);

		if (gasDesorptionRate <= 0. && internalPressure <= 0.)
			yDot[0] = 0.;
		else{
			// Sigmoid continuous microstructural venting modeling membrane micro-fessuration
			final double matrixPermeabilityK = 1. / (1. + Math.exp(5. * (vCurr - glutenTearingLimit)));
			// Viscous mechanical expansion velocity response equation (dV/dt)
			yDot[0] = (internalPressure / StrictMath.max(MIN_STIFFNESS_BOUND, dynamicStiffness)) * vCurr
				* matrixPermeabilityK;
		}
	}

}
