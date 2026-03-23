package io.github.mtrevisan.yeastcalculator.ode;

import io.github.mtrevisan.yeastcalculator.kinetics.KineticParameters;
import io.github.mtrevisan.yeastcalculator.physiology.EnvironmentalFactors;
import io.github.mtrevisan.yeastcalculator.physiology.PhysiologicalState;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


/**
 * Baranyi-Roberts (1994) ODE system for S. cerevisiae in bread dough.
 * <p>
 * State: [Q, N, S, P_CO2, E_ethanol]
 * Growth equations:
 * 	dQ/dt = MU_MAX · Q	← key: constant rate, not mu_eff! Q represents intracellular enzyme reserves that
 * 		replenish at the maximum rate, independently of external conditions - Baranyi (1994), eq. 2)
 * 	dN/dt = MU_MAX · γ_T · γ_aw · γ_E · α_lag · α_cap · N
 * 	dS/dt = −(1/Y_XS) · dN/dt
 * 	dP/dt = Y_PS · (−dS/dt)
 * 	dE/dt = Y_ES · (−dS/dt)
 * <p>
 * Where:
 * 	α_lag = Q / (Q + 1)                              — lag adjustment → 1 when Q >> 1
 * 	α_cap = max(0, 1 − N / N_MAX)                    — capacity adjustment → 0 at saturation
 * 	N_MAX = YIELD_BIOMASS_PER_SUGAR × sugarFraction  — substrate-limited ceiling
 */
public class BaranyiRobertsODE implements FirstOrderDifferentialEquations{

	private final EnvironmentalFactors factors;
	private final double temperature;


	public BaranyiRobertsODE(final EnvironmentalFactors factors, final double temperature){
		this.factors = factors;
		this.temperature = temperature;
	}

	@Override
	public int getDimension(){
		return KineticParameters.STATE_DIMENSION;
	}

	@Override
	public void computeDerivatives(final double time, final double[] state, final double[] derivatives){
		final double physiologicalState = Math.max(0., state[KineticParameters.IDX_PHYSIOLOGICAL_STATE]);
		final double biomassDensity = Math.max(0., state[KineticParameters.IDX_BIOMASS_DENSITY]);
		final double sugarConcentration = Math.max(0., state[KineticParameters.IDX_SUGAR_CONCENTRATION]);
		final double ethanolConcentration = Math.max(0., state[KineticParameters.IDX_ETHANOL_PRODUCED]);
		final double sugarEnzymatic = Math.max(0., state[KineticParameters.IDX_SUGAR_ENZYMATIC]);

		//── Environmental inhibition factors ──────────────────────────────
		final double totalSugarAvailable = sugarConcentration + sugarEnzymatic;
		final double gammaTemperature = factors.gammaTemperature(temperature);
		final double gammaWaterActivity = factors.gammaWaterActivity(totalSugarAvailable);
		final double gammaEthanol = factors.gammaEthanol(ethanolConcentration);

		//── Baranyi structural factors ────────────────────────────────────
		//lagAdjustment → 0 during lag phase (Q small), → 1 when fully adapted (Q >> 1)
		final double lagAdjustment = PhysiologicalState.lagAdjustmentFactor(physiologicalState);
		//capacityAdjustment → 0 when biomass saturates the substrate-limited ceiling
		final double dynamicBiomassCapacity = KineticParameters.YIELD_BIOMASS_PER_SUGAR * totalSugarAvailable;
		final double capacityAdjustment = Math.max(0., 1. - biomassDensity / dynamicBiomassCapacity);
		//binary sugar availability flag
		final double sugarAvailability = (totalSugarAvailable > KineticParameters.SUGAR_DEPLETION_THRESHOLD? 1.: 0.);

		//── Effective growth rate [h⁻¹] ──────────────────────────────────
		final double effectiveGrowthRate = KineticParameters.MU_MAX_REF
			* gammaTemperature
			* gammaWaterActivity
			* gammaEthanol
			* lagAdjustment
			* capacityAdjustment
			* sugarAvailability;

		//── ODE derivatives ───────────────────────────────────────────────
		//dQ/dt uses MU_MAX_REF (not effectiveGrowthRate): Q evolves at the maximum rate regardless of environment — it
		// represents intracellular enzyme reserves that replenish at a fixed rate (Baranyi (1994), eq. 2)
		derivatives[KineticParameters.IDX_PHYSIOLOGICAL_STATE] = KineticParameters.MU_MAX_REF * physiologicalState;

		final double biomassGrowthRate = effectiveGrowthRate * biomassDensity;
		derivatives[KineticParameters.IDX_BIOMASS_DENSITY] = biomassGrowthRate;

		final double sugarConsumptionRate = -(1. / KineticParameters.YIELD_BIOMASS_PER_SUGAR) * biomassGrowthRate;
		derivatives[KineticParameters.IDX_SUGAR_CONCENTRATION] = sugarConsumptionRate;
		derivatives[KineticParameters.IDX_CO2_PRODUCED] = KineticParameters.YIELD_CO2_PER_SUGAR * (-sugarConsumptionRate);
		derivatives[KineticParameters.IDX_ETHANOL_PRODUCED] = KineticParameters.YIELD_ETHANOL_PER_SUGAR
			* (-sugarConsumptionRate);

		// T_MIN, T_OPT, T_MAX for amylase
		final double gammaTemperatureAmylase = factors.gammaAmylaseTemperature(temperature);
		// dS_enz/dt = k · (1 − S_enz/S_max) · γ_T
		final double sugarGenerationRate = KineticParameters.AMYLASE_ACTIVITY_RATE
			* (1. - sugarEnzymatic / KineticParameters.SUGAR_GENERATION_MAX)
			* gammaTemperatureAmylase;
		derivatives[KineticParameters.IDX_SUGAR_ENZYMATIC] = sugarGenerationRate;
	}

}