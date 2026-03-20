package io.github.mtrevisan.yeastcalculator.model;


public class FlourMoistureModel{

	private static final String[] FLOUR_TYPES = {
		"wheat",
		"wheat semolina",
		"wheat semolina fine",
		"durum wheat",
		"durum wheat semolina",
		"durum wheat semolina fine",
		"rye",
		"einkorn",
		"emmer",
		"spelt",
		"buckwheat",
		"barley",
		"chestnut"
	};

	private static final double[] WM_BASE = {
		0.0682, 0.0671, 0.0676, 0.0694, 0.0683, 0.0688, 0.0791, 0.0698, 0.0712, 0.0694, 0.0614, 0.0748, 0.0572
	};
	private static final double[] WM_PROTEIN_COEFF = {
		0.112, 0.108, 0.11, 0.104, 0.1, 0.102, 0.088, 0.118, 0.114, 0.116, 0.072, 0.081, 0.063
	};
	private static final double[] WM_FIBER_COEFF = {
		0.184, 0.162, 0.173, 0.156, 0.144, 0.15, 0.298, 0.201, 0.194, 0.189, 0.157, 0.241, 0.148
	};
	private static final double[] C_BASE = {
		7.21, 6.84, 7.02, 6.98, 6.61, 6.79, 4.87, 7.44, 7.28, 7.36, 5.12, 5.34, 4.63
	};
	private static final double[] D_HC = {
		38420, 37800, 38100, 39200, 38600, 38900, 35100, 39800, 39100, 39400, 34200, 35800, 32900
	};
	private static final double[] K_BASE = {
		0.8124, 0.8241, 0.8182, 0.8198, 0.8312, 0.8255, 0.8641, 0.8087, 0.8114, 0.8098, 0.7921, 0.8412, 0.7698
	};
	private static final double[] D_HK = {
		12840, 13120, 12980, 13210, 13450, 13330, 11240, 13640, 13420, 13510, 10980, 11760, 10210
	};

	private static final double R = 8.314;


	private static double[] normalize(final double[] v, final int n){
		double sum = 0.;
		for(int i = 0; i < n; i ++)
			sum += v[i];

		final double[] out = new double[n];
		for(int i = 0; i < n; i ++)
			out[i] = v[i] / sum;
		return out;
	}

	private static double arrhenius(final double A, final double T, final double Tref){
		return Math.exp(-(A / R) * (1 / (T + 273.15) - 1 / (Tref + 273.15)));
	}

	private static double arrhenius(final double A, final double T){
		return Math.exp(-A / (R * (T + 273.15)));
	}

	private static int matchIndex(final String type){
		for(int i = 0; i < FLOUR_TYPES.length; i ++)
			if(FLOUR_TYPES[i].equalsIgnoreCase(type))
				return i;

		throw new IllegalArgumentException("Unknown flour type: " + type);
	}

	public static double[] compute(final double[] fractionsRaw, final double[][] flourMatrix,
		final String[] flourTypes, final double flourTemperature, final double airRelativeHumidity){
		final int flours = fractionsRaw.length;
		final double[] fractions = normalize(fractionsRaw, flours);

		final double aw = airRelativeHumidity;

		double strictlyBoundWaterDB = 0.;
		double flourMoistureDB = 0.;
		for(int i = 0; i < flours; i ++){
			final int idx = matchIndex(flourTypes[i]);
			final double flourProtein = flourMatrix[i][4];
			final double flourFiber = flourMatrix[i][6];
			final double flourAsh = flourMatrix[i][7];

			final double Wm = WM_BASE[idx] + WM_PROTEIN_COEFF[idx] * flourProtein + WM_FIBER_COEFF[idx] * flourFiber;
			strictlyBoundWaterDB += Wm * fractions[i];

			final double C = (C_BASE[idx] + 41.2 * flourProtein + 92.1 * flourAsh
				- 112.4 * flourFiber / (1 + 8.2 * flourFiber)) * arrhenius(D_HC[idx], flourTemperature, 20.);
			final double K = (K_BASE[idx] - 0.724 * flourProtein - 0.801 * flourFiber - 1.423 * flourAsh)
				* arrhenius(D_HK[idx], flourTemperature, 20.);
			final double Ka = Math.min(K * aw, 0.999);
			flourMoistureDB += (Wm * C * Ka / ((1 - Ka) * (1 + (C - 1) * Ka))) * fractions[i];
		}

		return new double[]{
			flourMoistureDB,
			strictlyBoundWaterDB
		};
	}

}
