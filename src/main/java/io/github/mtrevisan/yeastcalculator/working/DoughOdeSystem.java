package io.github.mtrevisan.yeastcalculator.working;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


final class DoughOdeSystem implements FirstOrderDifferentialEquations{

	private final double yDry;
	private final double[][] stages;
	private final double[] folds;
	private final double stiffnessIndexBase;
	private final double saltK;
	private final double oilK;

	DoughOdeSystem(final double yDry, final double[][] stages, final double[] folds,
		final double stiffnessIndexBase, final double saltK, final double oilK){
		this.yDry = yDry;
		this.stages = stages;
		this.folds = folds;
		this.stiffnessIndexBase = stiffnessIndexBase;
		this.saltK = saltK;
		this.oilK = oilK;
	}

	@Override
	public int getDimension(){
		// 3 continuous state variables:
		// y[0] = Dough Volume (V)
		// y[1] = Latent enzyme activation (Lag)
		// y[2] = Residual metabolic sugar substrate
		return 3;
	}

	@Override
	public void computeDerivatives(final double t, final double[] y, final double[] yDot){
		final double vCurr = y[0];
		final double lagCurr = y[1];
		final double sugarCurr = y[2];

		// 1. Identify current environmental stage based on continuous time 't' (hours)
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

		// 2. Structural Dough Rheology (Gluten viscoelastic relaxation)
		final double stiffnessTarget = stiffnessIndexBase / ((rhCurr < 0.60) ? (0.50 + 0.50 * (rhCurr / 0.60)) : 1.);
		final double proteolysisRate = 0.002 * Math.exp(0.07 * (tCurr - 20.));
		final double stiffnessRelaxed = stiffnessTarget * Math.exp(-proteolysisRate * t);
		final double stiffnessIndexDynamic = stiffnessRelaxed + (stiffnessIndexBase - stiffnessRelaxed) * Math.exp(-1.8 * t);

		// 3. Microorganism Biological Kinetics
		final double alphaBio = YeastFermentationModel.calculateThermalEfficiency(tCurr);
		final double muBio = YeastFermentationModel.calculateBiomassGrowthRate(yDry, alphaBio, lagCurr, sugarCurr, saltK, oilK);

		// Enzyme adjustment state rate (dLag/dt)
		yDot[1] = alphaBio;

		// Substrate exhaustion rate (dSugar/dt)
		if(sugarCurr <= 0.){
			yDot[2] = 0.;
		}
		else{
			yDot[2] = -(muBio * 0.015);
		}

		// 4. Bio-Mechanical Coupling & Gas Desorption (Macro-structural expansion)
		final double muEff = muBio * stiffnessIndexBase / stiffnessIndexDynamic;
		final double pressureFactor = (vCurr - 1.) * stiffnessIndexDynamic / stiffnessIndexBase;
		final double leakingCoeff = 1. / (1. + Math.exp(-12. * (pressureFactor - 1.8)));
		final double leakingK = 1. - (leakingCoeff * 0.4);

		// Volumetric structural expansion rate (dV/dt)
		yDot[0] = muEff * vCurr * leakingK;
	}

}
