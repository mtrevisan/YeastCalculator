package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;
import java.util.List;


public class SimulationInputs{

	// Lookup structures for specific flour type GAB constants
	private static final List<String> FLOUR_TYPES_LOOKUP = Arrays.asList(
		"wheat", "wheat semolina", "wheat semolina fine", "durum wheat", "durum wheat semolina",
		"durum wheat semolina fine", "rye", "einkorn", "emmer", "spelt", "buckwheat", "barley", "chestnut"
	);


	public static class GabResult{
		public double flourStrictlyBoundWater;
		public double flourActiveWater;
	}

	// GAB Model constants associated to flour types
	private static final double[] WM_BASE_LOOKUP = {0.062, 0.06, 0.061, 0.063, 0.062, 0.062, 0.071, 0.064, 0.065, 0.063, 0.055, 0.068, 0.05};
	private static final double[] C_GAB_LOOKUP   = {12.5, 11.8, 12.1, 13.2, 12.6, 12.8, 8.4, 13.0, 12.4, 11.9, 6.2, 9.5, 4.8};
	private static final double[] K_GAB_LOOKUP   = {0.78, 0.79, 0.78, 0.76, 0.77, 0.77, 0.82, 0.75, 0.77, 0.79, 0.85, 0.81, 0.88};

	private static final double[] BASE_LOOKUP = {1.0, 0.96, 0.98, 0.85, 0.82, 0.85, 0.35, 0.52, 0.68, 0.78, 0.05, 0.25, 0.02};


	private double[] fractionsRaw = {
		0.70,
		0.20,
		0.10
	};
	public Object[][] flourMatrix = {
		//W		P/L	sugar		prot.	fat		fiber		ashes		type
		{300.,	0.55,	0.015,	0.13,	0.012,	0.02,		0.0055,	"wheat"},
		{220.,	0.60,	0.010,	0.12,	0.015,	0.025,	0.0065,	"wheat semolina fine"},
		{110.,	0.40,	0.030,	0.10,	0.020,	0.06,		0.015,	"rye"}
	};
	public double[][] stagesRaw = {
		//temp	RH		duration
		{24.,		0.75,	2.},
		{28.,		0.80,	3.5},
		{4.,		0.70,	12.}
	};
	public double[] foldsRaw = {1.2, 2.5};

	public double doughWater = 0.65;
	public double doughSalt = 0.02;
	public double doughOil = 0.03;
	public double yeastMoisture = 0.70;

	public double flourTemperature = 22.;
	public double airRelativeHumidity = 0.55;


	int getFlourCount(){
		return fractionsRaw.length;
	}

	double[] getFractions(){
		final int flours = fractionsRaw.length;
		double sumFractions = 0.;
		for(final double f : fractionsRaw)
			sumFractions += f;
		final double[] fractions = new double[flours];
		for(int i = 0; i < flours; i ++)
			fractions[i] = fractionsRaw[i] / sumFractions;
		return fractions;
	}

	/**
	 * Resolves the internal index for a given flour type name.
	 * Returns 0 (default fallback) if the type is unregistered.
	 */
	public int getFlourTypeIndex(final String type){
		final int index = FLOUR_TYPES_LOOKUP.indexOf(type);
		return (index != -1? index: 0);
	}

	/**
	 * Executes the Guggenheim-Anderson-de Boer (GAB) Water Adsorption Model over a multi-flour mixture.
	 */
	public static GabResult calculateMoisture(final SimulationInputs in, final double[] fractions){
		double mixTotalDB = 0.;
		double mixBoundDB = 0.;

		final double aw = StrictMath.max(StrictMath.min(in.airRelativeHumidity, 0.95), 0.1);
		final double temperatureCoeff = 1. - 0.0025 * (in.flourTemperature - 20.);

		for(int i = 0; i < fractions.length; i++){
			final double protein = ((Number)in.flourMatrix[i][3]).doubleValue();
			final double fiber = ((Number)in.flourMatrix[i][6]).doubleValue();
			final String type = (String)in.flourMatrix[i][7];

			// Ask SimulationInputs to resolve coefficients via internal lookup indexes
			final int fidx = in.getFlourTypeIndex(type);
			final double wmBase = in.getWmBase(fidx);
			final double cGab = in.getCGab(fidx);
			final double kGab = in.getKGab(fidx);

			// Monolayer moisture content on a Dry Basis (DB)
			final double wmDb = wmBase + 0.085 * protein + 0.12 * fiber;

			// Thermodynamic equilibrium sorption formula
			final double uEquilibriumDB = (wmDb * cGab * kGab * aw) / ((1. - kGab * aw) * (1. - kGab * aw + cGab * kGab * aw));

			// Apply thermal drift scale and bound it between [0.08, 0.2] DB
			final double uTotalDB = StrictMath.max(StrictMath.min(uEquilibriumDB * temperatureCoeff, 0.2), 0.08);

			mixTotalDB += uTotalDB * fractions[i];
			mixBoundDB += wmDb * fractions[i];
		}

		// Dry Basis to Wet Basis conversions -> x / (1 + x)
		final double moistureTotal = mixTotalDB / (1. + mixTotalDB);
		final double flourStrictlyBoundWater = mixBoundDB / (1. + mixBoundDB);
		final double flourActiveWater = moistureTotal - flourStrictlyBoundWater;

		final GabResult res = new GabResult();
		res.flourStrictlyBoundWater = flourStrictlyBoundWater;
		res.flourActiveWater = flourActiveWater;
		return res;
	}

	public double getWmBase(final int typeIndex){
		return WM_BASE_LOOKUP[typeIndex];
	}

	public double getCGab(final int typeIndex){
		return C_GAB_LOOKUP[typeIndex];
	}

	public double getKGab(final int typeIndex){
		return K_GAB_LOOKUP[typeIndex];
	}

	public double getBaseLookupValue(final int typeIndex) {
		return BASE_LOOKUP[typeIndex];
	}

}
