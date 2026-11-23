package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Enumeration representing distinct flours with their GAB isotherm baseline properties
 * and lookup indices mirroring the original Excel matrix model.
 * <p>
 * The GAB parameters are:
 * </p>
 * <ul>
 * <li>{@code wmBase}: monolayer moisture of a flour with no protein and no fiber contribution [g_water/g_dry_flour];
 *    the model adds {@code 0.085 · protein + 0.12 · fiber} on top of it, so protein- and fiber-rich flours (pulses,
 *    nuts) have a lower base value;</li>
 * <li>{@code cGab}: Guggenheim constant, related to the heat of sorption of the monolayer;</li>
 * <li>{@code kGab}: correction factor of the heat of sorption of the multilayer (always below 1).</li>
 * </ul>
 * <p>
 * Flours are grouped by botanical family. The values of the non-wheat flours are approximations.
 * </p>
 */
public enum FlourType{

	// --- Wheat and its relatives (gluten-forming) ---
	WHEAT("wheat", 0.062, 11., 0.89),
	WHEAT_SEMOLINA("wheat semolina", 0.06, 10.5, 0.9),
	WHEAT_SEMOLINA_FINE("wheat semolina fine", 0.061, 10.8, 0.9),
	DURUM_WHEAT("durum wheat", 0.063, 12.5, 0.88),
	DURUM_WHEAT_SEMOLINA("durum wheat semolina", 0.062, 12., 0.89),
	DURUM_WHEAT_SEMOLINA_FINE("durum wheat semolina fine", 0.062, 12.2, 0.89),
	KHORASAN("khorasan", 0.063, 12., 0.89),
	SPELT("spelt", 0.063, 11., 0.89),
	EMMER("emmer", 0.065, 11.2, 0.89),
	EINKORN("einkorn", 0.064, 11.5, 0.89),
	TRITICALE("triticale", 0.063, 11., 0.89),
	RYE("rye", 0.071, 13.5, 0.87),
	BARLEY("barley", 0.068, 10.5, 0.88),

	// --- Other cereals (no gluten network) ---
	OAT("oat", 0.06, 10., 0.9),
	MAIZE("maize", 0.058, 9.5, 0.89),
	RICE("rice", 0.065, 12., 0.88),
	MILLET("millet", 0.058, 9., 0.9),
	SORGHUM("sorghum", 0.06, 10., 0.89),

	// --- Pseudocereals ---
	BUCKWHEAT("buckwheat", 0.055, 9., 0.91),
	QUINOA("quinoa", 0.062, 10., 0.89),
	AMARANTH("amaranth", 0.06, 9., 0.9),
	TEFF("teff", 0.062, 10., 0.89),

	// --- Legumes (protein- and fiber-rich: low base monolayer) ---
	CHICKPEA("chickpea", 0.04, 14., 0.88),
	LENTIL("lentil", 0.035, 13., 0.88),
	PEA("pea", 0.035, 13., 0.88),
	SOY("soy", 0.025, 12., 0.87),

	// --- Nuts and fruits (high fat: low moisture uptake) ---
	CHESTNUT("chestnut", 0.05, 8., 0.92),
	ALMOND("almond", 0.01, 7., 0.92),
	COCONUT("coconut", 0.01, 6., 0.92);


	private final String label;
	private final double wmBase;
	private final double cGab;
	private final double kGab;


	FlourType(final String label, final double wmBase, final double cGab, final double kGab){
		this.label = label;
		this.wmBase = wmBase;
		this.cGab = cGab;
		this.kGab = kGab;
	}


	public String getLabel(){
		return label;
	}

	public double getWmBase(){
		return wmBase;
	}

	public double getcGab(){
		return cGab;
	}

	public double getkGab(){
		return kGab;
	}

}
