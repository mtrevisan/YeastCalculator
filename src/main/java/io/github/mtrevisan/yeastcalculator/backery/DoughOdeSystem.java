package io.github.mtrevisan.yeastcalculator.backery;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


/**
 * System of equations applying the Rosso CTMI temperature function and enzymatic kinetics.
 */
public class DoughOdeSystem implements FirstOrderDifferentialEquations{

	private static final double T_MIN = 2.0;
	private static final double T_OPT = 32.0;
	private static final double T_MAX = 43.0;

	private static final double MU_OPT = 0.45;
	private static final double K_S = 0.005;
	private static final double Y_XS = 0.12;
	private static final double MAINTENANCE_M = 0.01;
	private static final double Y_VS = 350.0;

	private final double currentTemperature;
	private final double activeWater;
	private final double maxGasPotential;
	private final SimulationInputs in;
	private final FoldEventHandler foldHandler;


	public DoughOdeSystem(double currentTemperature, double activeWater, double maxGasPotential,
		SimulationInputs in, FoldEventHandler foldHandler){
		this.currentTemperature = currentTemperature;
		this.activeWater = activeWater;
		this.maxGasPotential = maxGasPotential;
		this.in = in;
		this.foldHandler = foldHandler;
	}


	@Override
	public int getDimension(){
		return 3;
	}

	@Override
	public void computeDerivatives(double t, double[] y, double[] yDot){
		double X = Math.max(0.0, y[0]);
		double S = Math.max(0.0, y[1]);
		double Vgas = Math.max(0.0, y[2]);

		double gammaT = calculateRossoGammaT(currentTemperature);
		double saltConcentration = in.getRecipe().getWaterRatio() > 0? in.getRecipe().getSaltRatio() / in.getRecipe().getWaterRatio(): 0;
		double gammaSalt = Math.max(0.0, 1.0 - 4.0 * saltConcentration);
		double gammaWater = Math.min(1.0, activeWater / 0.14);

		double muEff = MU_OPT * (S / (K_S + S)) * gammaT * gammaSalt * gammaWater;
		double kd = 0.01 + (S <= 0.0? 0.05: 0.0) + (currentTemperature >= T_MAX? 0.5: 0.0);

		double enzymeActivity = Math.max(0.0, (currentTemperature - T_MIN) / (T_OPT - T_MIN));
		double rMalt = 0.002 * in.getRecipe().getMaltRatio() * in.getRecipe().getMaltPollakUnit() * enzymeActivity;

		yDot[0] = (muEff - kd) * X;
		double sugarConsumedByYeast = ((muEff / Y_XS) + MAINTENANCE_M) * X;
		yDot[1] = rMalt - sugarConsumedByYeast;

		double dynamicMaxPotential = maxGasPotential * foldHandler.getCurrentGasPotentialModifier();
		double retentionEfficiency = Math.max(0.0, 1.0 - Math.pow(Vgas / dynamicMaxPotential, 2));
		yDot[2] = Y_VS * sugarConsumedByYeast * retentionEfficiency;
	}

	private double calculateRossoGammaT(double T){
		if(T <= T_MIN || T >= T_MAX) return 0.0;
		double num = (T - T_MAX) * Math.pow(T - T_MIN, 2);
		double den = (T_OPT - T_MIN) * ((T_OPT - T_MIN) * (T - T_OPT) - (T_OPT - T_MAX) * (T_OPT + T_MIN - 2.0 * T));
		return den == 0? 0.0: num / den;
	}

}
