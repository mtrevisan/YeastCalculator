package io.github.mtrevisan.yeastcalculator.simulation;

import io.github.mtrevisan.yeastcalculator.domain.FlourInput;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
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
	// Norrish-like constant for the aw depression caused by the hydrated flour matrix (starch, gluten, fibers),
	// calibrated so that an unsalted dough at 60% hydration sits at aw ~0.975
	private static final double K_MATRIX = 0.019;
	// Molar mass of NaCl [g/mol]
	private static final double NACL_MOLAR_MASS = 58.44;
	// Ethanol toxicity ceiling (~60 g/kg dough)
	private static final double ETHANOL_MAX = 0.06;
	// Ghose & Tyagi exponent
	private static final double ETHANOL_N = 0.5;
	// Gay-Lussac yield (~48%)
	private static final double Y_ETHANOL_S = 0.48;

	private static final double K_S = 0.005;
	private static final double Y_XS = 0.12;
	// Maintenance coefficient in [min^-1]. Converted from ~0.03 g_substrate/(g_X * h) (typical for
	// S. cerevisiae under anaerobic conditions, Verduyn et al., 1991, Antonie van Leeuwenhoek 60: 59-64)
	// to per-minute units: 0.03/60 = 0.0005.
	private static final double MAINTENANCE_M = 0.0005;
	// Calibrated gas yield factor per gram of metabolized sugar substrate [mL/g_sugar]
	private static final double Y_VS = 265.;
	// Base maximum specific growth rate adjusted from hourly limits to continuous minute scaling [min^-1]
	private static final double BASE_MU_OPT = 0.0085;

	// Diastatic malt sugar production coefficient.
	// Derivation: 1 Pollak unit (PK) produces approximately 1 mg of maltose per 10 g of flour
	// in 1 h at 30 °C (EBC standard method). With PK = 120 and 1.5% diastatic malt, the expected
	// sugar generation rate is therefore roughly 1.8e-5 g_sugar/g_flour/min. The previous
	// coefficient of 0.002 produced ~3.2e-3 g/g/min, i.e. roughly 180 times too much, which
	// saturated the dough almost instantly and pinned the residual sugar at the artificial
	// saturation cutoff. The corrected coefficient below yields a generation rate of ~3.2e-5
	// g/g/min at the reference temperature, which is the correct order of magnitude for a
	// 4-8 h direct fermentation.
	// References: Hoseney, Principles of Cereal Science and Technology, 3rd ed., chapter on
	// amylase activity; EBC Method 4.12 (Pollak units); Jayaram et al., 2014, J. Agric. Food
	// Chem. 62: 10312-10319 (maltose accumulation curves in fermented dough).
	private static final double MALT_RATE_COEFFICIENT = 2.e-5;
	// Product inhibition constant for alpha-amylase. Maltose and glucose exert weak product
	// inhibition on the enzyme, with reported Ki values in the range 5-20 mM (Greenwood &
	// Milne, 1968, Adv. Carbohydr. Chem. 23: 281). 5 mM of maltose corresponds to approximately
	// 0.005 g of sugar per gram of flour in a standard 60%-hydration dough. This constant is
	// used to build a Michaelis-like decay modifier in place of the previous hard cutoff, so
	// the enzyme keeps producing sugar (at a reduced rate) even when free sugars are abundant.
	private static final double MALT_PRODUCT_INHIBITION_KI = 0.005;

	private final double currentTemperature;
	private final double currentStageRH;
	// Flour moisture strictly bound in the GAB monolayer (non-solvent) [g/g_flour]
	private final double boundWater;
	// Flour moisture above the GAB monolayer (free to act as solvent) [g/g_flour]
	private final double activeWater;
	private final double maxGasPotential;
	private final SimulationInputs in;
	private final FoldEventHandler foldHandler;


	public DoughOdeSystem(final double currentTemperature, final double currentStageRH,
			final GabMoistureModel.GabResult gab, final double maxGasPotential, final SimulationInputs in,
			final FoldEventHandler foldHandler){
		this.currentTemperature = currentTemperature;
		this.currentStageRH = currentStageRH;
		this.boundWater = gab.flourStrictlyBoundWater;
		this.activeWater = gab.flourActiveWater;
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
		final double waterRatio = in.recipe().waterRatio();
		final double saltRatio = in.recipe().saltRatio();

		// Solvent water: added water plus the flour moisture above the GAB monolayer (monolayer water is strongly
		// bound and does not dissolve solutes) [g/g_flour]
		final double solventWater = waterRatio + activeWater;
		// Dry flour solids [g/g_flour]
		final double drySolids = 1. - boundWater - activeWater;

		// Molality of salt (moles of NaCl per kg of solvent water); all ratios are per gram of flour, hence the
		// grams-to-kilograms conversion of the water mass
		final double molesSalt = saltRatio / NACL_MOLAR_MASS;
		final double kgWater = solventWater / 1000.;
		final double molalitySalt = (kgWater > 0.? molesSalt / kgWater: 0.);

		// aw due to salt: exp(-0.018 * phi * nu * m) where phi=0.93, nu=2
		final double awSalt = Math.exp(-0.018 * 0.93 * 2. * molalitySalt);

		// aw due to the hydrated flour matrix: the more water per unit of dry solids, the closer to pure water
		final double awMatrix = (solventWater > 0.? Math.exp(-K_MATRIX * drySolids / solventWater): 0.);

		// The total aw is the product of the matrix and the salt contributions
		final double totalAw = Math.min(1., awMatrix * awSalt);

		// Cardinal inhibition via aw (Rosso et al., 1995)
		double gammaAw = 0.;
		if(totalAw > AW_MIN)
			gammaAw = (totalAw - AW_MIN) / (1. - AW_MIN);

		// --- 2. ETHANOL TOXICITY INHIBITION (Ghose & Tyagi, 1979) ---
		double gammaEthanol = 1. - Math.pow(Math.min(1., EtOH / ETHANOL_MAX), ETHANOL_N);

		// Mineral/Ash nutrient boosting scaling.
		// The ash contribution is weighted by the flour fraction in the blend.
		double blendAsh = 0.;
		final double[] fractions = in.fractions();
		final FlourInput[] flourMatrix = in.flourMatrix();
		for(int i = 0; i < flourMatrix.length; i ++)
			blendAsh += fractions[i] * flourMatrix[i].ash();
		final double adjustedMuOpt = BASE_MU_OPT * (1. + Math.min(0.25, blendAsh * 15.));
		final double gammaT = calculateRossoGammaT(currentTemperature);

		// Fat/Oil membrane screening penalty
		final double oilRatio = in.recipe().oilRatio();
		final double gammaOilInhibition = 1. / (1. + 2.5 * oilRatio);

		// Combined growth kinetic rate
		final double muEff = adjustedMuOpt * (s / (K_S + s)) * gammaT * gammaAw * gammaEthanol * gammaOilInhibition;

		final double speedYeastConsumption = ((muEff / Y_XS) + MAINTENANCE_M) * x;

		// Cell mortality rate
		final double kd = cellMortalityRate(s);

		// Diastatic malt sugar conversion rate.
		// The previous formulation used a hard cutoff at 0.001 g/g (0.1%) which forced the model
		// to always plateau at exactly 0.1% residual sugar, making it impossible for the optimizer
		// to reach any higher residual sugar target. The corrected formulation uses a Michaelis-
		// like product inhibition term (1 / (1 + s/Ki)) so the enzyme keeps producing sugar at a
		// reduced rate even at high sugar concentrations, which is what actually happens in a
		// fermenting dough (Jayaram et al., 2014).
		final double enzymeActivity = Math.max(0., (currentTemperature - T_MIN) / (T_OPT - T_MIN));
		final double sugarSaturationModifier = 1. / (1. + s / MALT_PRODUCT_INHIBITION_KI);
		final double rMalt = MALT_RATE_COEFFICIENT * in.recipe().maltRatio() * in.recipe().maltPollakUnit()
			* enzymeActivity * sugarSaturationModifier;

		// Calculate derivatives
		double dX_dt = (muEff - kd) * x;
		double dS_dt = rMalt - speedYeastConsumption;
		double dEtOH_dt = Y_ETHANOL_S * speedYeastConsumption;

		// Rheological volume evolution with smooth Bloksma leakage boundaries
		// Integrating retention efficiency with dynamic structural modifiers
		final double dynamicMaxPotential = maxGasPotential * foldHandler.getCurrentGasPotentialModifier()
			* skinningModifier(currentStageRH);

		// K_tear defines the slope of the bubble collapse (polymer rheology)
		final double kTear = 1.5;
		final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (vGas - dynamicMaxPotential)));
		double dVgas_dt = Y_VS * speedYeastConsumption * retentionEfficiency;

		if(x <= 0. && dX_dt < 0.)
			dX_dt = 0.;
		if(s <= 0. && dS_dt < 0.)
			dS_dt = 0.;

		yDot[0] = dX_dt;
		yDot[1] = dS_dt;
		yDot[2] = dVgas_dt;
		yDot[3] = dEtOH_dt;
	}

	/**
	 * Gas retention penalty caused by surface skinning (crust drying) under low ambient relative humidity.
	 */
	public static double skinningModifier(final double relativeHumidity){
		return (relativeHumidity >= 0.70? 1.: Math.max(0.6, 1. - (0.70 - relativeHumidity) * 0.8));
	}

	private double cellMortalityRate(double s){
		// --- CONTINUOUS CELLULAR MORTALITY RATE (MICROBIOLOGICAL LITERATURE) ---
		// 1. Basal physiological mortality rate under ideal conditions
		// Derived from standard S. cerevisiae baseline parameters (approx. 0.01% per minute) [min^-1]
		final double kd0 = 0.0001;

		// 2. Starvation-induced mortality modeling.
		// Values corrected against Boer, Amini & Botstein (2008, PNAS 105: 6936-6941) and
		// Gresham et al. (2011, Genetics 187: 209-220): death-phase onset delay of 60-80 h,
		// death rate 0.007-0.04 h^-1, half-lives 27.7 h (leucine) to 337 h (phosphate).
		// The previous 0.02 min^-1 (= 1.2 h^-1) was 30-170 times too fast.
		final double maxStarvationDeathRate = 0.0005;
		// Half-saturation below which starvation stress accelerates. Matches the Km of the
		// MALx1 maltose transporter (~5 mM, ~0.05% of flour mass) and the glucose derepression
		// threshold of the MAL regulon (0.1-0.4 C-mol/L).
		final double kStarveAffinity = 0.0005;
		final double kdStarvation = maxStarvationDeathRate * (1. - (s / (kStarveAffinity + s)));

		// 3. Thermal death/inactivation kinetics (Bigelow / Arrhenius extension)
		// Avoids discontinuous step changes at arbitrary thresholds like T_MAX.
		// Cell membrane denaturing and thermal protein degradation scale exponentially above 38-40 °C.
		double kdThermal = 0.;
		if(currentTemperature > 38.){
			// kD60 is the reference inactivation rate at 60 °C; zValue is the thermal sensitivity factor
			// Death velocity rate at base reference of 60 °C [min^-1]
			final double kD60 = 0.24;
			// Temperature change required to alter thermal death by one log factor [°C]
			final double zValue = 4.8;
			kdThermal = kD60 * Math.pow(10., (currentTemperature - 60.) / zValue);
		}

		// Combined continuous cellular death rate vector mapping
		return kd0 + kdStarvation + kdThermal;
	}

	private double calculateRossoGammaT(final double T){
		if(T <= T_MIN || T >= T_MAX)
			return 0.;
		final double num = (T - T_MAX) * Math.pow(T - T_MIN, 2.);

		final double den = (T_OPT - T_MIN) * ((T_OPT - T_MIN) * (T - T_OPT) - (T_OPT - T_MAX) * (T_OPT + T_MIN - 2. * T));
		return (den == 0.? 0.: num / den);
	}

}
