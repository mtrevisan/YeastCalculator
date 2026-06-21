package io.github.mtrevisan.yeastcalculator.simulation;

import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


/**
 * System of equations applying the Rosso CTMI temperature function and enzymatic kinetics.
 */
public class DoughOdeSystem implements FirstOrderDifferentialEquations{

	// Cardinal temperature parameters (Rosso et al., 1995)
	private static final double T_MIN = 2.;
	private static final double T_OPT = 32.;
	private static final double T_MAX = 43.;

	// Biological constants for S. cerevisiae
	// Minimum aw limit for growth
	private static final double AW_MIN = 0.88;
	// Ethanol toxicity ceiling (~60 g/kg dough)
	private static final double ETHANOL_MAX = 0.06;
	// Ghose & Tyagi exponent
	private static final double ETHANOL_N = 0.5;
	// Gay-Lussac yield (~48%)
	private static final double Y_ETHANOL_S = 0.48;

	private static final double K_S = 0.005;
	private static final double Y_XS = 0.12;
	private static final double MAINTENANCE_M = 0.01;
	private static final double Y_VS = 350.;

	private final double currentTemperature;
	private final double currentStageRH;
	private final double activeWater;
	private final double maxGasPotential;
	private final SimulationInputs in;
	private final FoldEventHandler foldHandler;


	public DoughOdeSystem(final double currentTemperature, final double currentStageRH, final double activeWater,
			final double maxGasPotential, final SimulationInputs in, final FoldEventHandler foldHandler){
		this.currentTemperature = currentTemperature;
		this.currentStageRH = currentStageRH;
		this.activeWater = activeWater;
		this.maxGasPotential = maxGasPotential;
		this.in = in;
		this.foldHandler = foldHandler;
	}


	@Override
	public int getDimension(){
		return 4;
	}

	@Override
	public void computeDerivatives(final double t, final double[] y, final double[] yDot){
		final double x = Math.max(0., y[0]);
		final double s = Math.max(0., y[1]);
		final double vGas = Math.max(0., y[2]);
		final double EtOH = Math.max(0., y[3]);

		// --- 1. ROSS MODEL (1975) FOR WATER ACTIVITY (aw) ---
		final double waterRatio = in.getRecipe().getWaterRatio();
		final double saltRatio = in.getRecipe().getSaltRatio();

		// Molality of salt (moles of NaCl per kg of water in the mixture)
		final double molesSalt = saltRatio / 58.44;
		final double kgWater = waterRatio;
		final double molalitySalt = (kgWater > 0.? molesSalt / kgWater: 0.);

		// aw due to salt: exp(-0.018 * phi * nu * m) where phi=0.93, nu=2
		final double awSalt = Math.exp(-0.018 * 0.93 * 2. * molalitySalt);

		// The total aw is the product of the aw of the native flour and the salt
		// Typical value of an unsalted dough
		final double awFlourBase = 0.96;
		final double totalAw = Math.min(1., awFlourBase * awSalt);

		// Cardinal inhibition via aw (Rosso et al., 1995)
		double gammaAw = 0.;
		if(totalAw > AW_MIN)
			gammaAw = (totalAw - AW_MIN) / (1. - AW_MIN);

		// --- 2. ETHANOL TOXICITY INHIBITION (Ghose & Tyagi, 1979) ---
		double gammaEthanol = 1. - Math.pow(Math.min(1., EtOH / ETHANOL_MAX), ETHANOL_N);

		// Ash (nutrient mineral boosting) & Temperature corrections
		double blendAsh = 0.;
		// Simplified to blend 100% in main
		for(final FlourInput f : in.getFlourMatrix())
			blendAsh += f.getAsh();
		final double adjustedMuOpt = 0.45 * (1. + Math.min(0.15, blendAsh * 10.));
		final double gammaT = calculateRossoGammaT(currentTemperature);

		// Fat/Oil membrane screening penalty
		final double oilRatio = in.getRecipe().getOilRatio();
		final double gammaOilInhibition = Math.max(0.7, 1. - oilRatio * 1.2);

		// Combined growth kinetic rate
		final double muEff = adjustedMuOpt * (s / (K_S + s)) * gammaT * gammaAw * gammaEthanol * gammaOilInhibition;

		final double speedYeastConsumption = ((muEff / Y_XS) + MAINTENANCE_M) * x;

		// Cell mortality rate
		final double kd = 0.01 + (s <= 0.? 0.05: 0.) + (currentTemperature >= T_MAX? 0.5: 0.);

		// Diastatic malt sugar conversion rate
		final double enzymeActivity = Math.max(0., (currentTemperature - T_MIN) / (T_OPT - T_MIN));
		final double rMalt = 0.002 * in.getRecipe().getMaltRatio() * in.getRecipe().getMaltPollakUnit() * enzymeActivity;

		// Differential Equations
		// dX/dt
		yDot[0] = (muEff - kd) * x;
		// dS/dt
		yDot[1] = rMalt - speedYeastConsumption;
		// dEtOH/dt (Alcohol Accumulation)
		yDot[3] = Y_ETHANOL_S * speedYeastConsumption;

		// --- 3. GLUTEN TEARING: BLOKSMA / CONSIDÈRE SIGMOID RETENTION ---
		// Integrating retention efficiency with dynamic structural modifiers
		final double skinningModifier = (currentStageRH >= 0.70 ? 1.: Math.max(0.6, 1. - (0.70 - currentStageRH) * 0.8));
		final double dynamicMaxPotential = maxGasPotential * foldHandler.getCurrentGasPotentialModifier()
			* skinningModifier;

		// K_tear defines the slope of the bubble collapse (polymer rheology)
		final double kTear = 15.;
		final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (vGas - dynamicMaxPotential)));

		// dVgas/dt
		yDot[2] = Y_VS * speedYeastConsumption * retentionEfficiency;
	}

	private double calculateRossoGammaT(final double T){
		if(T <= T_MIN || T >= T_MAX)
			return 0.;
		final double num = (T - T_MAX) * Math.pow(T - T_MIN, 2.);

		final double den = (T_OPT - T_MIN) * ((T_OPT - T_MIN) * (T - T_OPT) - (T_OPT - T_MAX) * (T_OPT + T_MIN - 2. * T));
		return (den == 0.? 0.: num / den);
	}

}
