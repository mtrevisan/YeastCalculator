package io.github.mtrevisan.yeastcalculator.simulation;

import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


/**
 * System of equations applying the Rosso CTMI temperature function and enzymatic kinetics.
 */
public class DoughOdeSystem implements FirstOrderDifferentialEquations{

	private static final double T_MIN = 2.;
	private static final double T_OPT = 32.;
	private static final double T_MAX = 43.;

	private static final double K_S = 0.005;
	private static final double Y_XS = 0.12;
	private static final double MAINTENANCE_M = 0.01;
	private static final double Y_VS = 350.;

	private final double currentTemperature;
	private final double activeWater;
	private final double maxGasPotential;
	private final SimulationInputs in;
	private final FoldEventHandler foldHandler;


	public DoughOdeSystem(final double currentTemperature, final double activeWater, final double maxGasPotential,
			final SimulationInputs in, final FoldEventHandler foldHandler){
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
	public void computeDerivatives(final double t, final double[] y, final double[] yDot){
		final double x = Math.max(0., y[0]);
		final double s = Math.max(0., y[1]);
		final double vGas = Math.max(0., y[2]);

		// 1. Dynamic mu_opt calculation influenced by ash (nutrients) of the flour mixture
		double blendAsh = 0.;
		final double[] fractions = in.getFractions();
		final FlourInput[] matrix = in.getFlourMatrix();
		for(int i = 0; i < matrix.length; i ++)
			blendAsh += fractions[i] * matrix[i].getAsh();
		// Ash/bran acts as a mineral growth booster (up to +15% mu_opt)
		final double adjustedMuOpt = 0.45 * (1. + Math.min(0.15, blendAsh * 10.));

		// 2. Temperature Modifier (Rosso CTMI)
		final double gammaT = calculateRossoGammaT(currentTemperature);

		// 3. Salt Inhibition Modifier
		final double saltConcentration = (in.getRecipe().getWaterRatio() > 0.
			? in.getRecipe().getSaltRatio() / in.getRecipe().getWaterRatio()
			: 0.);
		final double gammaSalt = Math.max(0., 1. - 4. * saltConcentration);

		// 4. Water Activity Modifier
		final double gammaWater = Math.min(1., activeWater / 0.14);

		// 5. ACTIVATION Dough Recipe.getOil Ratio (Membrane inhibitor effect on yeast)
		// An excess of oil (e.g. > 5%) partially shields the cells by reducing the osmotic exchange
		final double oilRatio = in.getRecipe().getOilRatio();
		final double gammaOilInhibition = Math.max(0.7, 1. - (oilRatio * 1.2));

		// Combined effective growth kinetics
		final double muEff = adjustedMuOpt * (s / (K_S + s)) * gammaT * gammaSalt * gammaWater * gammaOilInhibition;

		// Cell mortality rate
		final double kd = 0.01 + (s <= 0.? 0.05: 0.) + (currentTemperature >= T_MAX? 0.5: 0.);

		// Enzymatic production of sugars from malt
		final double enzymeActivity = Math.max(0., (currentTemperature - T_MIN) / (T_OPT - T_MIN));
		final double rMalt = 0.002 * in.getRecipe().getMaltRatio() * in.getRecipe().getMaltPollakUnit() * enzymeActivity;

		// --- APPLICATION OF DIFFERENTIAL EQUATIONS ---
		yDot[0] = (muEff - kd) * x;

		final double sugarConsumedByYeast = ((muEff / Y_XS) + MAINTENANCE_M) * x;
		yDot[1] = rMalt - sugarConsumedByYeast;

		// Integrating retention efficiency with dynamic structural modifiers
		final double dynamicMaxPotential = maxGasPotential * foldHandler.getCurrentGasPotentialModifier();
		final double retentionEfficiency = Math.max(0., 1. - Math.pow(vGas / dynamicMaxPotential, 2.));

		yDot[2] = Y_VS * sugarConsumedByYeast * retentionEfficiency;
	}

	private double calculateRossoGammaT(final double T){
		if(T <= T_MIN || T >= T_MAX)
			return 0.;
		final double num = (T - T_MAX) * Math.pow(T - T_MIN, 2.);

		final double den = (T_OPT - T_MIN) * ((T_OPT - T_MIN) * (T - T_OPT) - (T_OPT - T_MAX) * (T_OPT + T_MIN - 2. * T));
		return (den == 0.? 0.: num / den);
	}

}
