/**
 * Copyright (c) 2026 Mauro Trevisan
 * <p>
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
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
 * Flours are grouped by botanical family. The values were tuned so that the equilibrium moisture at 20 °C and 30/50/65 % relative
 * humidity falls within the ranges reported for flours in the sorption-isotherm literature (see {@code GabMoistureModelTest});
 * the chestnut monolayer (about 0.06 g/g d.b.) follows Moreira et al., 2010. The values of the non-wheat flours remain approximations.
 * </p>
 */
public enum FlourType{

	// --- Wheat and its relatives (gluten-forming) ---
	WHEAT("wheat", 0.062, 11., 0.89),
	WHOLE_WHEAT("whole wheat", 0.057, 11.5, 0.88),	// NEW
	WHEAT_SEMOLINA("wheat semolina", 0.06, 10.5, 0.9),
	WHEAT_SEMOLINA_FINE("wheat semolina fine", 0.061, 10.8, 0.9),
	DURUM_WHEAT("durum wheat", 0.063, 12.5, 0.88),
	DURUM_WHEAT_SEMOLINA("durum wheat semolina", 0.062, 12., 0.89),
	DURUM_WHEAT_SEMOLINA_FINE("durum wheat semolina fine", 0.062, 12.2, 0.89),
	KHORASAN("khorasan", 0.063, 12., 0.89),
	WHOLE_KHORASAN("whole khorasan", 0.056, 12.5, 0.88),	// NEW
	SPELT("spelt", 0.063, 11., 0.89),
	WHOLE_SPELT("whole spelt", 0.055, 11.5, 0.89),	// NEW
	EMMER("emmer", 0.065, 11.2, 0.89),
	WHOLE_EMMER("whole emmer", 0.057, 11.7, 0.88),	// NEW
	EINKORN("einkorn", 0.06, 11.5, 0.89),	// NEW
	WHOLE_EINKORN("whole einkorn", 0.058, 12., 0.86),	// NEW
	TRITICALE("triticale", 0.063, 11., 0.89),
	RYE("rye", 0.067, 13.5, 0.87),	// NEW
	WHOLE_RYE("whole rye", 0.063, 14., 0.82),	// NEW
	BARLEY("barley", 0.064, 10.5, 0.88),	// NEW

	// --- Other cereals (no gluten network) ---
	OAT("oat", 0.055, 10., 0.9),	// NEW
	MAIZE("maize", 0.058, 9.5, 0.89),
	RICE("rice", 0.065, 12., 0.88),
	MILLET("millet", 0.058, 9., 0.9),
	FOXTAIL_MILLET("foxtail millet", 0.058, 9.5, 0.9),	// NEW
	SORGHUM("sorghum", 0.06, 10., 0.89),

	// --- Pseudocereals ---
	BUCKWHEAT("buckwheat", 0.055, 9., 0.91),
	QUINOA("quinoa", 0.062, 10., 0.89),
	AMARANTH("amaranth", 0.06, 9., 0.9),
	TEFF("teff", 0.062, 10., 0.89),

	// --- Legumes (protein- and fiber-rich: low base monolayer) ---
	CHICKPEA("chickpea", 0.036, 14., 0.88),	// NEW
	LENTIL("lentil", 0.035, 13., 0.88),
	PEA("pea", 0.032, 13., 0.87),	// NEW
	FAVA_BEAN("fava bean", 0.03, 13., 0.88),	// NEW
	SOY("soy", 0.025, 12., 0.87),

	// --- Nuts and fruits (high fat: low moisture uptake) ---
	CHESTNUT("chestnut", 0.049, 14., 0.92),	// NEW
	WATER_CHESTNUT("water chestnut", 0.065, 11., 0.9),	// NEW
	ALMOND("almond", 0.01, 7., 0.92),
	COCONUT("coconut", 0.01, 3.5, 0.92),	// NEW
	HAZELNUT("hazelnut", 0.01, 7., 0.92),

	// --- Seeds (high fat and protein: low moisture uptake) ---
	FLAXSEED("flaxseed", 0.005, 8., 0.91),
	HEMP("hemp", 0.005, 8., 0.91),
	SESAME("sesame", 0.007, 4.5, 0.92);	// NEW


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
