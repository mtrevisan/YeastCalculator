package io.github.mtrevisan.yeastcalculator.working;

import java.util.Arrays;
import java.util.List;


public final class FlourRegistry{

	private static final List<String> TYPES = Arrays.asList(
		"wheat", "wheat semolina", "wheat semolina fine", "durum wheat", "durum wheat semolina",
		"durum wheat semolina fine", "rye", "einkorn", "emmer", "spelt", "buckwheat", "barley", "chestnut"
	);

	private static final double[] WM_BASE = {0.062, 0.06, 0.061, 0.063, 0.062, 0.062, 0.071, 0.064, 0.065, 0.063, 0.055, 0.068, 0.05};
	private static final double[] C_GAB = {12.5, 11.8, 12.1, 13.2, 12.6, 12.8, 8.4, 13.0, 12.4, 11.9, 6.2, 9.5, 4.8};
	private static final double[] K_GAB = {0.78, 0.79, 0.78, 0.76, 0.77, 0.77, 0.82, 0.75, 0.77, 0.79, 0.85, 0.81, 0.88};
	private static final double[] BASE = {1.0, 0.96, 0.98, 0.85, 0.82, 0.85, 0.35, 0.52, 0.68, 0.78, 0.05, 0.25, 0.02};

	public static class FlourProperties{
		public final double wmBase;
		public final double cGab;
		public final double kGab;
		public final double baseLookup;

		public FlourProperties(final double wmBase, final double cGab, final double kGab, final double baseLookup){
			this.wmBase = wmBase;
			this.cGab = cGab;
			this.kGab = kGab;
			this.baseLookup = baseLookup;
		}
	}

	/**
	 * Resolves chemical coefficients for a flour type name.
	 * If type is unknown, returns the default parameters (index 0 - wheat).
	 */
	public static FlourProperties resolveProperties(final String type){
		int idx = TYPES.indexOf(type);
		if(idx == -1)
			idx = 0;
		return new FlourProperties(WM_BASE[idx], C_GAB[idx], K_GAB[idx], BASE[idx]);
	}

}
