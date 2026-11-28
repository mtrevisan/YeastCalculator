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
 * <li>{@code kGab}: correction factor of the heat of sorption of the multilayer (always below 1);</li>
 * <li>{@code compositionFactor}: scale (default 1) of the protein and fiber contribution to the monolayer; fat-rich kernels and seeds
 *    (hazelnut, chia) bind far less water per unit of protein and fiber, so their published monolayers (Wm about 0.0155 and
 *    0.024-0.029 g/g d.b.) are reproduced with a factor of 0.4.</li>
 * </ul>
 * <p>
 * Flours are grouped by botanical family. The values were tuned so that the equilibrium moisture at 20 °C and 30/50/65 % relative
 * humidity falls within the ranges reported for flours in the sorption-isotherm literature (see {@code GabMoistureModelTest});
 * whole wheat and teff follow the GAB fits of Martin-Santos et al. 2012 and Abdissa et al. 2020, the chestnut monolayer (about 0.06 g/g d.b.) the desorption branch of Moreira et al. 2010 / Chenlo et al. The values of the non-wheat flours remain approximations.
 * </p>
 */
public enum FlourType{

	/*
	= LET(flourType; S12;
flourProtein; S13; flourFiber; S14; flourAsh; S15;
flourTemperature; S17; airRelativeHumidity; S18;
damageType; S16;

clamp; LAMBDA(x;lo;hi; IF(x < lo; lo; IF(x > hi; hi; x)));

fidx; IFERROR(XMATCH(flourType; {"wheat","whole wheat","wheat semolina","wheat semolina fine","durum wheat","durum wheat semolina","durum wheat semolina fine","khorasan","whole khorasan","spelt","whole spelt","emmer","whole emmer","einkorn","whole einkorn","triticale","rye","whole rye","barley","whole barley","oat","whole oat","maize","maize fine","maize medium","maize coarse","maize white","maize red","rice","whole rice","millet","whole millet","foxtail millet","sorghum","buckwheat","whole buckwheat","quinoa","amaranth","teff","chickpea","lentil","pea","fava bean","soy","white bean","black bean","red lentil","green pea","lupin","chestnut","acorn","water chestnut","almond","coconut","hazelnut","walnut","pistachio","peanut","beech mast","flaxseed","hemp","pumpkin seed","sunflower seed","chia","sesame"}); 1);
ff; IF(OR(flourType = "hazelnut"; flourType = "chia"); 0.4; 1);

wm_base; CHOOSE(fidx; 0.062; 0.057; 0.06; 0.061; 0.063; 0.062; 0.062; 0.063; 0.056; 0.063; 0.055; 0.065; 0.057; 0.06; 0.058; 0.063; 0.067; 0.063; 0.064; 0.058; 0.055; 0.054; 0.058; 0.059; 0.058; 0.057; 0.058; 0.058; 0.065; 0.06; 0.058; 0.055; 0.058; 0.06; 0.055; 0.052; 0.062; 0.06; 0.049; 0.036; 0.035; 0.032; 0.03; 0.025; 0.033; 0.032; 0.035; 0.032; 0.022; 0.049; 0.054; 0.065; 0.01; 0.01; 0.005; 0.01; 0.01; 0.01; 0.01; 0.005; 0.005; 0.005; 0.005; 0.005; 0.007) + ff * (0.085 * flourProtein + 0.12 * flourFiber);
C_GAB; CHOOSE(fidx; 11; 13.6; 10.5; 10.8; 12.5; 12; 12.2; 12; 12.5; 11; 11.5; 11.2; 11.7; 11.5; 12; 11; 13.5; 14; 10.5; 10.5; 10; 10; 9.5; 9.5; 9.5; 9.5; 9.5; 9.5; 12; 12; 9; 9; 9.5; 10; 9; 9.5; 10; 9; 10.4; 14; 13; 13; 13; 12; 13; 13; 13; 13; 12; 14; 14; 11; 7; 3.5; 7; 7; 7; 7; 7; 8; 8; 8; 8; 8; 4.5);
K_GAB; CHOOSE(fidx; 0.89; 0.733; 0.9; 0.9; 0.88; 0.89; 0.89; 0.89; 0.88; 0.89; 0.89; 0.89; 0.88; 0.89; 0.86; 0.89; 0.87; 0.82; 0.88; 0.88; 0.9; 0.9; 0.89; 0.89; 0.89; 0.89; 0.89; 0.89; 0.88; 0.88; 0.9; 0.9; 0.9; 0.89; 0.91; 0.9; 0.89; 0.9; 0.7255; 0.88; 0.88; 0.87; 0.88; 0.87; 0.88; 0.88; 0.88; 0.87; 0.87; 0.92; 0.92; 0.9; 0.92; 0.92; 0.98; 0.92; 0.92; 0.92; 0.92; 0.91; 0.91; 0.91; 0.91; 0.97; 0.92);

carbohydrates; MAX(0; 1 - flourProtein - flourFiber - flourAsh);
starchAdj; SWITCH(damageType; "high"; 0.006; "low"; -0.003; "legume"; 0.004; 0) * carbohydrates;

WmDB; wm_base + 0.085 * flourProtein + 0.12 * flourFiber + starchAdj;

aw; clamp(airRelativeHumidity; 0.1; 0.95);
tempFactor; EXP(-0.004 * (flourTemperature - 20));
Wm_t; WmDB * tempFactor;

tmp; K_GAB * aw;
X_dry_raw; Wm_t * C_GAB * tmp / ((1 - tmp) * (1 + (C_GAB - 1) * tmp));
X_dry; clamp(X_dry_raw; 0.05; 0.35);

M_wet; X_dry / (1 + X_dry);
dryMatter; 1 / (1 + X_dry);
boundWater; dryMatter * Wm_t;
activeWater; MAX(0; M_wet - boundWater);

activeWater
)

type	maize
protein	12.5%
fiber	3.0%
ashes	0.3%
damage type	low
flour temperature [°C]	18.0
air relative humidity	59.0%
flour humidity (wet basis)	6.7%

1. Select "low" for:
• Refined flours produced by industrial roller mills (e.g., "type 00" or "type 0" flours, or highly refined flours made from einkorn, spelt, or barley). The smooth steel rollers are designed to gently separate the bran, minimizing damage to the starch.
• Soft-textured grains (e.g., einkorn, spelt, oats). Because their endosperm structure is naturally friable, these grains break down easily without subjecting the starch to intense mechanical stress.
• Sifted flours intended for long fermentation processes, where excessive starch damage would be counterproductive (as it would lead to premature structural collapse due to the enzymatic activity of amylases).

2. Select "high" for:
• All stone-ground flours (whole-grain, "type 1," "type 2," or whole-grain flours made from rye, buckwheat, sorghum, or millet). Stone milling exerts significant pressure and friction, which generates heat and mechanically crushes—and ruptures—a very high percentage of starch granules.
• Very hard grains (e.g., corn, sorghum, emmer). These require greater milling force to be reduced into flour, resulting in a higher degree of structural starch damage.
• Whole-grain flours in general. The physical presence of bran fragments—which act like sharp edges during the milling process—promotes the mechanical rupture of the surrounding starch granules.

3. Select "legume" for:
• Fava bean flour (or any other flour derived from legumes, such as chickpea, lentil, or pea flour).
• Chestnut flour or water chestnuts.

	 */
	// Literature notes (details and caveats: docs/sorption-isotherms-research.md). "No source found" means that no published
	// sorption isotherm could be verified for that flour: the values are analogies with the closest verified flour.

	// --- Wheat and its relatives (gluten-forming) ---
	// Moreira et al. 2010, wheat flour, 20 °C : Wm 0.087, C 9.7, K 0.751 (adsorption). Not adopted: the current values already match its EMC within ~1 point, and its Wm is above that of whole wheat (Martin-Santos). https://doi.org/10.1016/j.indcrop.2010.04.021
	WHEAT("wheat", 0.062, 11., 0.89),
	// Martin-Santos et al. 2012, whole wheat flour, 20 °C: Wm 0.080, C 13.6, K 0.733. https://doi.org/10.1111/j.1365-2621.2012.02996.x
	WHOLE_WHEAT("whole wheat", 0.057, 13.6, 0.733),
	// Erbas et al. 2005, semolina, 20 °C (only m0 seen, basis to be checked). https://doi.org/10.1016/j.jfoodeng.2004.07.017
	WHEAT_SEMOLINA("wheat semolina", 0.06, 10.5, 0.9),
	// Erbas et al. 2005 (see WHEAT_SEMOLINA).
	WHEAT_SEMOLINA_FINE("wheat semolina fine", 0.061, 10.8, 0.9),
	// No source found for durum flour; closest: Erbas et al. 2005 (semolina). https://doi.org/10.1016/j.jfoodeng.2004.07.017
	DURUM_WHEAT("durum wheat", 0.063, 12.5, 0.88),
	// Erbas et al. 2005 (semolina). https://doi.org/10.1016/j.jfoodeng.2004.07.017
	DURUM_WHEAT_SEMOLINA("durum wheat semolina", 0.062, 12., 0.89),
	// Erbas et al. 2005 (semolina). https://doi.org/10.1016/j.jfoodeng.2004.07.017
	DURUM_WHEAT_SEMOLINA_FINE("durum wheat semolina fine", 0.062, 12.2, 0.89),
	// No source found (analogy with wheat).
	KHORASAN("khorasan", 0.063, 12., 0.89),
	// No source found (analogy with whole wheat).
	WHOLE_KHORASAN("whole khorasan", 0.056, 12.5, 0.88),
	// No source found (analogy with wheat).
	SPELT("spelt", 0.063, 11., 0.89),
	// No source found (analogy with whole wheat).
	WHOLE_SPELT("whole spelt", 0.055, 11.5, 0.89),
	// No source found (analogy with wheat).
	EMMER("emmer", 0.065, 11.2, 0.89),
	// No source found (analogy with whole wheat).
	WHOLE_EMMER("whole emmer", 0.057, 11.7, 0.88),
	// No GAB data found (einkorn flakes were fitted with other models only).
	EINKORN("einkorn", 0.06, 11.5, 0.89),
	// No GAB data found.
	WHOLE_EINKORN("whole einkorn", 0.058, 12., 0.86),
	// No source found (analogy with wheat).
	TRITICALE("triticale", 0.063, 11., 0.89),
	// No verified numbers; candidate: Ocieczek and Kostek, wheat and rye flour type 200, 20 °C (paper not accessible).
	RYE("rye", 0.067, 13.5, 0.87),
	// No verified numbers (see RYE).
	WHOLE_RYE("whole rye", 0.063, 14., 0.82),
	// Maia et al. 2016, barley grain (not flour), 15-50 °C, Wm 0.092-0.057 kg/kg d.b. (abstract only). https://doi.org/10.11301/jsfe.17.131
	BARLEY("barley", 0.064, 10.5, 0.88),
	// Maia et al. 2016 (see BARLEY).
	WHOLE_BARLEY("whole barley", 0.058, 10.5, 0.88),

	// --- Other cereals (no gluten network) ---
	// No verified numbers; candidate: Brett et al. 2009, oat flour, 5/23/45 °C. https://doi.org/10.1007/s11483-009-9112-0
	OAT("oat", 0.055, 10., 0.9),
	// Brett et al. 2009 (see OAT).
	WHOLE_OAT("whole oat", 0.054, 10., 0.9),
	// No verified numbers; candidate: Samapundo et al. 2007, yellow dent maize kernels (paper not accessible).
	MAIZE("maize", 0.058, 9.5, 0.89),
	// No source for the grind variants (see MAIZE).
	MAIZE_FINE("maize fine", 0.059, 9.5, 0.89),
	// No source for the grind variants (see MAIZE).
	MAIZE_MEDIUM("maize medium", 0.058, 9.5, 0.89),
	// No source for the grind variants (see MAIZE).
	MAIZE_COARSE("maize coarse", 0.057, 9.5, 0.89),
	// No source for the colour variants (see MAIZE).
	MAIZE_WHITE("maize white", 0.058, 9.5, 0.89),
	// No source for the colour variants (see MAIZE).
	MAIZE_RED("maize red", 0.058, 9.5, 0.89),
	// No verified numbers; candidate: Brett et al. 2009, rice flour, 5/23/45 °C. https://doi.org/10.1007/s11483-009-9112-0
	RICE("rice", 0.065, 12., 0.88),
	// Brett et al. 2009 (see RICE).
	WHOLE_RICE("whole rice", 0.06, 12., 0.88),
	// No GAB data found (pearl millet grain was fitted with modified Henderson, Goneli et al. 2010).
	MILLET("millet", 0.058, 9., 0.9),
	// No GAB data found (see MILLET).
	WHOLE_MILLET("whole millet", 0.055, 9., 0.9),
	// No usable data (the only hit, https://doi.org/10.1007/s11947-017-2007-z, reports implausible monolayers).
	FOXTAIL_MILLET("foxtail millet", 0.058, 9.5, 0.9),
	// Galdeano et al. 2018, whole sorghum flour, 25 °C: Wm 0.069, C 2.76, K 0.714 (outlier, not adopted). https://doi.org/10.1590/1981-6723.11817
	SORGHUM("sorghum", 0.06, 10., 0.89),

	// --- Pseudocereals ---
	// Arslan Tontul 2020, buckwheat grain (hulled), 15 °C: Wm 4.86 % (abstract only). https://doi.org/10.15237/gida.GD19110
	BUCKWHEAT("buckwheat", 0.055, 9., 0.91),
	// Arslan Tontul 2020 (see BUCKWHEAT).
	WHOLE_BUCKWHEAT("whole buckwheat", 0.052, 9.5, 0.9),
	// Contreras Lopez et al. 2026, quinoa, 20-40 °C, monolayer 6.97-11.55 g/100 g (abstract only). https://doi.org/10.15446/rfnam.v79.124831
	QUINOA("quinoa", 0.062, 10., 0.89),
	// No verified numbers; candidate: Pagano and Mascheroni 2003, amaranth seeds. https://doi.org/10.1081/JFP-120020116
	AMARANTH("amaranth", 0.06, 9., 0.9),
	// Abdissa et al. 2020, teff flour, 27 °C: Mo 6.97, C 10.4, K 0.7255; the parameters follow it. https://repository.rit.edu/japr/vol12/iss1/1
	// Solomon and Zewdu 2016, teff seed, 10-40 °C (parameters not seen). https://doi.org/10.1016/j.eaef.2016.05.001
	TEFF("teff", 0.049, 10.4, 0.7255),

	// --- Legumes (protein- and fiber-rich: low base monolayer) ---
	// No GAB values seen; candidates: Menkov 2000 (seeds) and Durakova and Menkov 2005 (flour), J. Food Eng. (abstracts only).
	CHICKPEA("chickpea", 0.036, 14., 0.88),
	// No source found.
	LENTIL("lentil", 0.035, 13., 0.88),
	// No usable numbers (Jacob et al. 2019, whole seed, parameter table unreadable).
	PEA("pea", 0.032, 13., 0.87),
	// Karatas and Arslan 2022, fava bean seed, 25 °C: Xm 7.15, C 5.67, K 0.81 (units/table layout to be checked). https://doi.org/10.35234/fumbd.1126606
	FAVA_BEAN("fava bean", 0.03, 13., 0.88),
	// No source found (Diaz-Contreras 2019 used Halsey/Oswin, not GAB).
	SOY("soy", 0.025, 12., 0.87),
	// No source found for the cultivar (a damaged-bean study gives a monolayer of 5.1-9.8 % d.b. over 20-50 °C).
	WHITE_BEAN("white bean", 0.033, 13., 0.88),
	// No source found.
	BLACK_BEAN("black bean", 0.032, 13., 0.88),
	// No source found (same as LENTIL).
	RED_LENTIL("red lentil", 0.035, 13., 0.88),
	// No usable numbers (see PEA).
	GREEN_PEA("green pea", 0.032, 13., 0.87),
	// Unverified (snippets only): Vazquez et al. 2003, J. Food Eng., monolayer about 0.0497; Briceno-Sevillano et al. 2025,
	// tarwi flour, 0.087 g/g at 20 °C. https://doi.org/10.31910/rudca.v28.n1.2025.2540
	LUPIN("lupin", 0.022, 12., 0.87),

	// --- Nuts and fruits (high fat: low moisture uptake) ---
	// Moreira et al. 2010, chestnut flour, 20 °C: adsorption Wm 0.049, C 60.5, K 0.902; desorption Wm 0.062, C 23.7, K 0.859
	// (the model follows the desorption branch). https://doi.org/10.1016/j.indcrop.2010.04.021
	// Chenlo et al., Food Bioprocess Technol 4:1497-1504, desorption Wm 0.059 ± 0.007. https://doi.org/10.1007/s11947-009-0239-2
	CHESTNUT("chestnut", 0.049, 14., 0.92),
	// No source found (analogy with chestnut).
	ACORN("acorn", 0.054, 14., 0.92),
	// No usable source (only a Trapa pericarp powder, not the kernel).
	WATER_CHESTNUT("water chestnut", 0.065, 11., 0.9),
	// Unverified: Taitano and Singh 2012, Int. J. Food Studies 1:61-75 (paper not accessible).
	ALMOND("almond", 0.01, 7., 0.92),
	// No usable source (green pulp at 30-70 °C; copra fitted with Henderson).
	COCONUT("coconut", 0.01, 3.5, 0.92),
	// Guzey et al. 2001, raw Tombul hazelnut, 25 °C: Wm 0.0155 kg/kg d.b. (abstract only; the model is too high). https://doi.org/10.1081/DRT-100103945
	// Full text not reachable (only the abstract: Wm 0.0155, C and K missing). ASABE Trans., hybrid hazelnut kernels, desorption, temperature-dependent modified GAB (A 3.55, B 0.717, C 120.5, % d.b.; not used, no single-temperature fit): https://www.midwesthazelnuts.org/uploads/3/8/3/5/38359971/asabe_trans-emc_of_hazelnut_husks_shells_and_kernels.pdf
	HAZELNUT("hazelnut", 0.005, 7., 0.98, 0.4),
	// Togrul and Arslan 2007, J. Stored Prod. Res., walnut kernels (BET monolayer, no GAB seen).
	WALNUT("walnut", 0.01, 7., 0.92),
	// Soleimanifard and Hamdami 2018, 25-70 °C (table not readable, values not trusted). https://doi.org/10.17221/460/2016-CJFS
	PISTACHIO("pistachio", 0.01, 7., 0.92),
	// Kaur et al. 2023, monolayer 3.1-5.6 % at 25-45 °C, consistent with the model (abstract only). J. Sci. Food Agric. 103(6):3077-3092
	PEANUT("peanut", 0.01, 7., 0.92),
	// No source found (all hits were beech wood).
	BEECH_MAST("beech mast", 0.01, 7., 0.92),

	// --- Seeds (high fat and protein: low moisture uptake) ---
	// Valente et al. 2020, whole seed, desorption at 40-80 °C only (no 15-25 °C data). https://doi.org/10.5935/1806-6690.20200021
	FLAXSEED("flaxseed", 0.005, 8., 0.91),
	// No numbers found.
	HEMP("hemp", 0.005, 8., 0.91),
	// No numbers found (Mayor et al. 2005, GAB values not seen).
	PUMPKIN_SEED("pumpkin seed", 0.005, 8., 0.91),
	// No numbers found (only Halsey/Smith fits).
	SUNFLOWER_SEED("sunflower seed", 0.005, 8., 0.91),
	// Arslan-Tontul 2020, whole seed, 15/25/35 °C: monolayer 2.39-2.91 g/100 g (abstract only; the model is about double).
	// https://doi.org/10.1016/j.lwt.2019.108859
	CHIA("chia", 0.005, 8., 0.97, 0.4),
	// Menkov and Durakova 2007, semi-defatted sesame flour (10.6 % fat), 25 °C, adsorption points only. https://hrcak.srce.hr/file/48177
	SESAME("sesame", 0.007, 4.5, 0.92);


	private final String label;
	private final double wmBase;
	private final double cGab;
	private final double kGab;
	private final double compositionFactor;


	FlourType(final String label, final double wmBase, final double cGab, final double kGab){
		this(label, wmBase, cGab, kGab, 1.);
	}

	FlourType(final String label, final double wmBase, final double cGab, final double kGab, final double compositionFactor){
		this.label = label;
		this.wmBase = wmBase;
		this.cGab = cGab;
		this.kGab = kGab;
		this.compositionFactor = compositionFactor;
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

	public double getCompositionFactor(){
		return compositionFactor;
	}

	public double getkGab(){
		return kGab;
	}

}
