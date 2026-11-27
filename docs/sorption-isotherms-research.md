# Moisture sorption isotherms / GAB parameters for the FlourType constants

Scope: compare the GAB parameters in `FlourType.java` (wmBase, C, K), the monolayer formula and the literature bands in `GabMoistureModelTest.java` against published sorption data. No code was modified.

## 0. Method, reliability and caveats

- Model: `Wm_total = wmBase + 0.085*protein + 0.12*fiber` [g/g d.b.], GAB with C and K, temperature factor `1 - 0.0025*(T-20)`. Total moisture is clamped to [max(0.03, Wm), max(0.2, Wm)] d.b.
- Protein and fiber come from the test's `COMPOSITION` map. "Model Wm_tot" below is the total monolayer, which is the only quantity comparable with published Wm. "Model EMC" is the model's equilibrium moisture [% w.b.] at 20 C and aw 0.3 / 0.5 / 0.65, computed offline with the model formulas, clamp included.
- Reliability tags:
  - **[T]** the number was seen in fetched source text or a table.
  - **[A]** the number was seen only in an abstract or search snippet. It is second-hand and not verified.
  - **[C]** the number was computed by me from printed GAB parameters (EMC = Wm*C*K*aw / ((1-K*aw)(1+(C-1)K*aw))). It is not printed in the paper.
- Many publisher, PMC and Springer pages returned 403, CAPTCHA or login walls, so coverage of primary full texts is poor. Nothing was invented. Where a number could not be verified, the entry says so.
- Published Wm is usually total monolayer in g/g d.b. (or % d.b. or g/100 g, as printed). Units were not always stated, and I flag this where relevant.
- No source was found that reports EMC at exactly aw 0.3, 0.5 or 0.65 for any flour. EMC values below are [C] (from printed GAB constants) or measured points as printed.
- Most data are at 25-27 C, not 20 C. The GAB Wm of flours typically falls slightly with temperature, so 25 C values are slightly conservative.

## 1. Sources with verified or partly verified numbers

### 1.1 Wheat flour and chestnut flour, 20 C [T]
Moreira R., Chenlo F., Torres M.D., Prieto D.M. (2010). *Water adsorption and desorption isotherms of chestnut and wheat flours.* Industrial Crops and Products 32:252-257. DOI 10.1016/j.indcrop.2010.04.021. Numbers come from the author manuscript at https://minerva.usc.gal/bitstreams/6efde658-0bb7-4fa9-b966-66f4e22b6fcf/download. Both are flours, aw 0.09-0.91. Xm is in kg/kg d.b.

| Material (20 C) | Branch | Wm | C | K | [C] EMC at aw 0.3 / 0.5 / 0.65, % w.b. |
|---|---|---|---|---|---|
| Wheat flour | adsorption | 0.087 | 9.70 | 0.751 | 7.7 / 10.6 / 13.3 |
| Wheat flour | desorption | 0.089 | 20.1 | 0.722 | 8.8 / 11.3 / 13.7 |
| Chestnut flour | adsorption | 0.049 | 60.5 | 0.902 | 6.0 / 8.1 / 10.5 |
| Chestnut flour | desorption | 0.062 | 23.7 | 0.859 | 7.0 / 9.3 / 12.0 |

For chestnut, the 35 / 50 / 65 C adsorption values were also printed: Xm 0.046 / 0.044 / 0.043, C 21.5 / 9.14 / 6.01, K 0.931 / 0.982 / 0.999.

Secondary source [A]: Chenlo F., Moreira R., Prieto D.M. et al., *Desorption isotherms and net isosteric heat of chestnut flour and starch.* Food Bioprocess Technol 4:1497-1504. DOI 10.1007/s11947-009-0239-2. A snippet gives Xm 0.059 +/- 0.007 for chestnut flour and 0.103 +/- 0.021 for starch (20-65 C desorption). The Springer page could not be read, and the DOI suffix (2009) conflicts with the 2011 volume, so check the year.

### 1.2 Whole wheat flour, 10/20/30 C [T]
Martin-Santos J., Vioque M., Gomez R. (2012). *Thermodynamic properties of moisture adsorption of whole wheat flour. Calculation of net isosteric heat.* Int J Food Sci Technol 47(7):1487-1495. DOI 10.1111/j.1365-2621.2012.02996.x (https://academic.oup.com/ijfst/article/47/7/1487/7866264). Flour, T. aestivum, 100 % extraction.
- GAB, Table 2. The unit is as printed, m0 in g H2O/100 g dry matter.
  - 20 C: m0 8.0052 (=0.0801 g/g), C 13.5918, K 0.7333, R2 0.9958.
  - 10 C: 8.1313 / 15.8784 / 0.7262.
  - 30 C: 7.8782 / 11.7330 / 0.7361.
- Measured adsorption EMC at 20 C, g/100 g d.m.: aw 0.331 -> 8.7; 0.432 -> 10.2; 0.544 -> 12.0; 0.591 -> 12.9; 0.725 -> 15.7. I did not interpolate these. As a rough read, 0.65 lies between 12.9 and 15.7 d.m., about 14 d.m., which is about 12.4 % w.b.
- [C] EMC from the printed constants: 7.5 / 10.1 / 12.4 % w.b. at aw 0.3 / 0.5 / 0.65.
- About 2,300 characters at the end of the page were not read.

### 1.3 Semolina and farina, 20 C [T, partial]
Erbas M., Ertugay M.F., Certel M. (2005). *Moisture adsorption behaviour of semolina and farina.* J Food Eng 69(2):191-198. DOI 10.1016/j.jfoodeng.2004.07.017.
- GAB m0 at 20 C: semolina (hard wheat) 11.8, farina (soft wheat) 9.16, both in g/100 g (basis presumably d.m.). C and K were not seen.
- Safe storage is below 75 % RH at 20 C (about 12.5 % moisture).
- m0 = 11.8 is far above every other wheat product (0.074-0.087). The unit or basis needs checking in the paper before it is used.

### 1.4 Sorghum flour, 25 C [T]
Galdeano M.C. et al. (2018). *Influence of milling and extrusion on the sorption properties of sorghum.* Braz J Food Technol 21:e2017118. DOI 10.1590/1981-6723.11817.
- Whole-grain milled flour, genotype SC319. Raw flour GAB: m0 0.069 g/g solids, C 2.762, K 0.714, R2 0.99.
- [C] EMC: 3.6 / 6.1 / 8.3 % w.b. at aw 0.3 / 0.5 / 0.65.
- These EMC are far below the usual 10-12 % w.b. at aw 0.5-0.6 for sorghum, and C = 2.76 is very low. Treat it as an outlier. It is not recommended for adoption without cross-checking.

### 1.5 Teff seed and flour, 27 C [T]
Abdissa B., Math R., Desham K., Korra S. (2020). *Moisture Sorption Behavior and Shelf Life Prediction of Teff Seed and Flour.* J Appl Packaging Res 12(1), Art. 1 (https://repository.rit.edu/japr/vol12/iss1/1). No DOI shown. Units of Mo are not stated (probably g/100 g d.b.).
- Flour: Mo 6.9656, C 10.4012, K 0.7255. [C] EMC 6.2 / 8.5 / 10.6 % w.b.
- Seed: Mo 7.5075, C 14.3886, K 0.7051. [C] EMC 7.0 / 9.3 / 11.4 % w.b.
- Teff seed, 10-40 C: Solomon W.K., Zewdu A.D. (2016). Eng Agric Environ Food 9(4):339-345. DOI 10.1016/j.eaef.2016.05.001. GAB R2 0.99, but the parameters were only in the PDF and were not seen.

### 1.6 Broad (fava) bean seed, 25/35/45 C [T, layout uncertain]
Karatas M., Arslan N. (2022). Firat Univ. Muh. Bil. Dergisi 34(2):719-731. DOI 10.35234/fumbd.1126606. Whole seed, Turkish text. The units are probably % d.b. The row labelling in the extracted table was shifted, so the GAB / BET assignment must be checked against the PDF.
- Adsorption at 25 C: Xm 7.154, C 5.674, K 0.810.
- Desorption at 25 C: Xm 9.198, C 6.304, K 0.757.
- [C] adsorption EMC at 25 C: 5.7 / 8.7 / 11.5 % w.b.
- Peleg and Chung-Pfost fitted best, not GAB. The isotherm is type III with hysteresis. The 45 C adsorption C value is a degenerate fit.

### 1.7 Sesame flour, 25 C [T]
Menkov N.D., Durakova A.G. (2007). *Moisture Sorption Isotherms of Sesame Flour at Several Temperatures.* Food Technol Biotechnol 45(1):96-100 (https://hrcak.srce.hr/file/48177). The material is semi-defatted commercial sesame flour, 10.6 % fat and 57.3 % protein, so it is **not** comparable with a full-fat sesame composition.
- Adsorption EMC at 25 C, % d.b., as printed:

| aw | 0.113 | 0.225 | 0.328 | 0.432 | 0.576 | 0.753 | 0.843 |
|---|---|---|---|---|---|---|---|
| EMC | 5.33 | 6.36 | 7.16 | 8.96 | 9.32 | 10.36 | 16.17 |

- Desorption at 25 C: 5.52 / 8.73 / 8.98 / 9.98 / 10.70 / 12.66 / 16.32 at the same aw.
- A temperature-dependent 5-parameter GAB form was fitted, not the classic Wm/C/K. The BET monolayer is about 5.5-5.9 % d.b., but its temperature assignment is garbled.
- Compared with the model's sesame (full-fat composition, 4.1 / 5.9 / 8.5 % w.b.), the measured adsorption at aw 0.328 is 7.16 % d.b., which is 6.7 % w.b.

### 1.8 Other items seen only as abstracts or snippets [A]
- **Hazelnut (Tombul, 25 C):** Guzey D. et al. (2001), Drying Technol 19:691-699, DOI 10.1081/DRT-100103945. GAB Wm 0.0155 kg/kg d.b. for raw nuts and 0.0100-0.0163 for roasted. C and K were not seen. Jung et al. (2018, J Food Sci; DOI not confirmed) gave a safe kernel moisture of about 4.2-5.0 % d.b. at aw <= 0.65 for Oregon kernels (cultivar values Barcelona 5.0, Yamhill 4.2, Jefferson 4.9).
- **Peanut:** Kaur et al. (2023), J Sci Food Agric 103(6):3077-3092. The GAB monolayer is 3.135-4.235 % (J87) and 3.906-5.640 % (TG37A) at 25-45 C. Kernel or pod was not stated, and C and K were not seen.
- **Chia (whole seed, 15/25/35 C):** Arslan-Tontul S. (2020), LWT 119:108859, DOI 10.1016/j.lwt.2019.108859. The monolayer is 2.39-2.91 g/100 g solid and EMC at the highest aw is about 18-20 g/100 g. The GAB table was not seen.
- **Flaxseed (whole, desorption, 40-80 C only):** Valente M.C.C. et al. (2020), Rev Cienc Agron 51(2), DOI 10.5935/1806-6690.20200021. Modified GAB Xm 3.78 / 2.64 / 2.33 at 40 / 60 / 80 C. The unit was not confirmed. [T] for the table, but outside the 15-25 C range.
- **Lupin:** Vazquez G., Chenlo F., Moreira R. (2003), *Sorption isotherms of lupine at different temperatures*, J Food Eng. A monolayer of about 0.0497 d.b. was quoted in a search summary only, at 25-45 C. The paper text was not seen. Briceno-Sevillano S.M. et al. (2025), Rev. U.D.C.A Act. & Div. Cient. 28(1), DOI 10.31910/rudca.v28.n1.2025.2540: tarwi (Lupinus mutabilis) flour, GAB at 20/30/40 C. The abstract says the monolayer is highest at 20 C, with 0.087 g/g d.m. for the defatted flour. C and K were not given.
- **Buckwheat grain, 15 C:** Arslan Tontul S. (2020), GIDA 45(2):253-261, DOI 10.15237/gida.GD19110. GAB Wm 4.86 % (basis not stated). C, K and EMC were not seen. The grain includes the hull.
- **Quinoa (Regalona), 20/30/40 C:** Contreras Lopez E. et al. (2026), Rev Fac Nac Agron Medellin 79:e124831, DOI 10.15446/rfnam.v79.124831. GAB was the best fit, with a monolayer of 6.97-11.55 g/100 g solids. The split by temperature, C and K were not seen. Seed or flour was not stated.
- **Barley seed (grain), 15-50 C:** Maia G.D. et al. (2016), Jpn J Food Eng 17(4):131-138, DOI 10.11301/jsfe.17.131. The GAB monolayer fell from 0.092 to 0.057 kg/kg d.b. as T rose from 15 to 50 C. The 20 C value, C and K were not seen.
- **Wheat flour (Hebrard et al. 2003, Powder Technol):** a snippet gave a GAB monolayer of about 7 % d.m., temperature unspecified.
- **Water caltrop (Trapa bispinosa) pericarp powder, not the kernel:** GAB monolayer 0.0561 / 0.0509 / 0.0409 g/g at 20 / 30 / 40 C. Authors, journal and DOI were not identified. This is a weak proxy at best.
- **Foxtail millet flour (Food Bioprocess Technol 2017, DOI 10.1007/s11947-017-2007-z):** a snippet gave GAB monolayers of 3.235-2.364 "g/g", which are implausibly high. Authors, title, C and K were not seen. Do not use.
- **Chickpea:** Menkov (2000, J Food Eng; seeds, 5-60 C) and Durakova and Menkov (2005, J Food Eng; flour, 10-40 C, BET monolayer) were seen as abstracts only. Neither fitted GAB parameters that I could see.
- **Pistachio:** Soleimanifard and Hamdami (2018), Czech J Food Sci 36(3):268-275, DOI 10.17221/460/2016-CJFS. Table 3 has GAB wm/C/K at 25-70 C, but the extraction was scrambled and the 25 C column could not be assigned. The values seen (kernels, 40-70 C, wm about 0.007-0.016) are not trusted.

### 1.9 Located but numbers not seen
Brett et al. (2009), Food Biophysics 4:151-157, DOI 10.1007/s11483-009-9112-0, oat and rice flour at 5/23/45 C. The Polish PTPS 2020 paper on corn, millet, oat and rice flours at 25 C. Ocieczek (2012), Acta Agrophysica 19(2):365, on wheat flour types 450/550/650 at 20 C. Ocieczek and Kostek on wheat and rye flour type 200 at 20 C. Murrieta-Pazos et al. on durum semolina fractions (HAL hal-01273456). Pagano and Mascheroni (2003), Int J Food Prop, DOI 10.1081/JFP-120020116, on amaranth seeds. Samapundo et al. (2007) on yellow dent maize kernels. Taitano and Singh (2012) on almonds (Int J Food Studies 1:61-75; J Food Process Eng 35:840-850). Toğrul and Arslan (2007), J Stored Prod Res, on walnut kernels (Peleg best, BET monolayer). These are the best targets if publisher access is available.

## 2. Coverage per FlourType constant

"Model" columns give the current parameters, the total monolayer and the model EMC at aw 0.3 / 0.5 / 0.65 (% w.b., 20 C). "Best source" is the best published data found. A dash means no source was found. Whole-grain and flour counterparts share a row where the same source (or the same absence) applies.

| Constant(s) | wmBase / C / K | Model Wm_tot | Model EMC 0.3 / 0.5 / 0.65 | Best source (material, T) | Published values | Verdict |
|---|---|---|---|---|---|---|
| WHEAT | 0.062 / 11 / 0.89 | 0.074 | 7.5 / 10.7 / 14.2 | Moreira 2010 [T], wheat flour, 20 C | Wm 0.087, C 9.7, K 0.751 (ads.); 0.089 / 20.1 / 0.722 (des.) | EMC agrees within about 1 point (0.65 is high by 0.5-0.9). Wm_tot is 0.013 lower and K much higher, and the two compensate. |
| WHOLE_WHEAT | 0.057 / 11.5 / 0.88 | 0.080 | 8.0 / 11.4 / 14.9 | Martin-Santos 2012 [T], whole flour, 20 C | Wm 0.0801, C 13.59, K 0.733 | Wm_tot matches. K and C do not, and the model overshoots EMC by about 1.3 (0.5) and 2.5 (0.65). |
| WHEAT_SEMOLINA, WHEAT_SEMOLINA_FINE | 0.060-0.061 / 10.5-10.8 / 0.90 | 0.074-0.075 | 7.4-7.6 / 10.7-10.9 / 14.3-14.5 | Erbas 2005 [T partial], semolina and farina, 20 C | m0 11.8 and 9.16 g/100 g | m0 is far above the model, and its basis is suspect. C and K are not seen. Unresolved. |
| DURUM_WHEAT, DURUM_WHEAT_SEMOLINA(_FINE) | 0.062-0.063 / 12-12.5 / 0.88-0.89 | 0.077-0.078 | 7.8-7.9 / 11.1-11.2 / 14.6 | Erbas 2005 (semolina) | as above | Not independently verified. Murrieta-Pazos and a Food Sci Technol Res 2012 paper (30-80 C) were not readable or are out of range. |
| KHORASAN, WHOLE_KHORASAN | 0.063 / 12 / 0.89; 0.056 / 12.5 / 0.88 | 0.077; 0.080 | 7.9 / 11.2 / 14.7; 8.1 / 11.4 / 14.9 | none | none | No source found. Values are analogies with wheat. |
| SPELT, WHOLE_SPELT | 0.063 / 11 / 0.89; 0.055 / 11.5 / 0.89 | 0.077; 0.079 | 7.8 / 11.1 / 14.7; 8.0 / 11.3 / 14.9 | none | none | No source found. |
| EMMER, WHOLE_EMMER | 0.065 / 11.2 / 0.89; 0.057 / 11.7 / 0.88 | 0.078; 0.080 | 7.9 / 11.3 / 14.9; 8.0 / 11.4 / 14.9 | none | none | No source found. |
| EINKORN, WHOLE_EINKORN | 0.060 / 11.5 / 0.89; 0.058 / 12 / 0.86 | 0.075; 0.082 | 7.6 / 10.9 / 14.4; 8.2 / 11.5 / 14.9 | einkorn flakes study (Bulgarian, 10-40 C), citation not confirmed | no GAB (Chung-Pfost, Oswin, Halsey, Henderson, BET) | No GAB data found. |
| TRITICALE | 0.063 / 11 / 0.89 | 0.076 | 7.6 / 10.9 / 14.4 | none | none | No source found. |
| RYE, WHOLE_RYE | 0.067 / 13.5 / 0.87; 0.063 / 14 / 0.82 | 0.077; 0.087 | 8.0 / 11.1 / 14.4; 8.7 / 11.8 / 14.9 | Ocieczek and Kostek (PDF blocked) | numbers not seen | No verified numbers. |
| BARLEY, WHOLE_BARLEY | 0.064 / 10.5 / 0.88; 0.058 / 10.5 / 0.88 | 0.077; 0.082 | 7.7 / 11.0 / 14.4; 8.1 / 11.5 / 15.1 | Maia 2016 [A], grain, 15-50 C | monolayer 0.092 (15 C) to 0.057 (50 C) kg/kg d.b. | Plausible, since the model sits inside that range. The 20 C value is not seen, and the paper is for grain. |
| OAT, WHOLE_OAT | 0.055 / 10 / 0.90; 0.054 / 10 / 0.90 | 0.078 | 7.8 / 11.2 / 14.9; 7.8 / 11.3 / 15.0 | Brett 2009 (not seen) | none | No verified numbers. |
| MAIZE, MAIZE_FINE/MEDIUM/COARSE/WHITE/RED | 0.057-0.059 / 9.5 / 0.89 | 0.068-0.069 | 6.7-6.8 / 9.8 / 13.1 | PTPS 2020 and Samapundo 2007 (not seen) | none | No verified numbers. There is no source for the variants (grind or colour). |
| RICE, WHOLE_RICE | 0.065 / 12 / 0.88; 0.060 / 12 / 0.88 | 0.072; 0.071 | 7.4 / 10.4 / 13.7; 7.3 / 10.3 / 13.5 | Brett 2009 (not seen); rough rice PMC4190263 (not opened) | none | No verified numbers. |
| MILLET, WHOLE_MILLET | 0.058 / 9 / 0.90; 0.055 / 9 / 0.90 | 0.072; 0.074 | 7.1 / 10.4 / 13.9; 7.2 / 10.6 / 14.2 | Goneli 2010 (pearl millet grain, modified Henderson, no GAB) | none | No GAB data found. |
| FOXTAIL_MILLET | 0.058 / 9.5 / 0.90 | 0.077 | 7.6 / 11.0 / 14.7 | FBT 2017 [A] | implausible "g/g" values | Unusable. No verified numbers. |
| SORGHUM | 0.060 / 10 / 0.89 | 0.076 | 7.5 / 10.8 / 14.3 | Galdeano 2018 [T], whole flour, 25 C | Wm 0.069, C 2.762, K 0.714 | Wm is close, but C and K are far off. The source is an outlier (see 1.4). |
| BUCKWHEAT, WHOLE_BUCKWHEAT | 0.055 / 9 / 0.91; 0.052 / 9.5 / 0.90 | 0.078; 0.076 | 7.7 / 11.2 / 15.1; 7.5 / 10.9 / 14.6 | Arslan Tontul 2020 [A], grain, 15 C | Wm 4.86 % (basis unknown) | The model is higher, but the source is hulled grain. C, K and EMC were not seen. |
| QUINOA | 0.062 / 10 / 0.89 | 0.082 | 8.1 / 11.6 / 15.4 | Contreras Lopez 2026 [A], 20-40 C | monolayer 6.97-11.55 g/100 g | The model is inside the range. Details unknown. |
| AMARANTH | 0.060 / 9 / 0.90 | 0.080 | 7.8 / 11.4 / 15.2 | Pagano and Mascheroni 2003 (not seen) | none | No verified numbers. |
| TEFF | 0.062 / 10 / 0.89 | 0.083 | 8.1 / 11.7 / 15.5 | Abdissa 2020 [T], flour, 27 C | Mo 6.97, C 10.40, K 0.7255 | The model overshoots EMC by 1.9 / 3.2 / 4.9 points. |
| CHICKPEA | 0.036 / 14 / 0.88 | 0.068 | 7.1 / 10.0 / 13.1 | Menkov 2000, Durakova and Menkov 2005 (abstracts) | no GAB values seen | No verified numbers. |
| LENTIL, RED_LENTIL | 0.035 / 13 / 0.88 | 0.069 | 7.2 / 10.1 / 13.3 | none (only lentil extrudates) | none | No source found. |
| PEA, GREEN_PEA | 0.032 / 13 / 0.87 | 0.070 | 7.2 / 10.1 / 13.1 | Jacob 2019 (whole seed; parameter table garbled) | unusable | No usable numbers. |
| FAVA_BEAN | 0.030 / 13 / 0.88 | 0.068 | 7.0 / 9.9 / 13.0 | Karatas and Arslan 2022 [T, layout uncertain], seed, 25 C | Xm 7.154, C 5.67, K 0.81 (ads.) | Wm agrees (0.0715 vs 0.068). C is lower, K is lower, and EMC is about 1.3 points lower than the model. |
| SOY | 0.025 / 12 / 0.87 | 0.066 | 6.8 / 9.6 / 12.6 | none (Diaz-Contreras 2019 used Halsey and Oswin) | none | No source found. |
| WHITE_BEAN, BLACK_BEAN | 0.033 / 13 / 0.88; 0.032 / 13 / 0.88 | 0.070 | 7.2-7.3 / 10.2 / 13.3-13.4 | Brazilian damaged-bean study [A] | monolayer 9.84 % to 5.10 % d.b. across 20-50 C (cultivar unknown) | Not applicable. Nothing found for black bean. |
| LUPIN | 0.022 / 12 / 0.87 | 0.066 | 6.8 / 9.6 / 12.5 | Vazquez 2003 [A]; Briceno-Sevillano 2025 [A] | 0.0497 (25-45 C, unverified); 0.087 defatted flour at 20 C | Inconclusive. The model lies between the two. |
| CHESTNUT | 0.049 / 14 / 0.92 | 0.063 | 6.8 / 9.6 / 12.9 | Moreira 2010 [T], flour, 20 C | see 1.1 | The model is close to the desorption branch. It is higher than adsorption by about 0.8 / 1.5 / 2.4 points. |
| ACORN | 0.054 / 14 / 0.92 | 0.065 | 7.1 / 10.1 / 13.4 | none | none | No source found. |
| WATER_CHESTNUT | 0.065 / 11 / 0.90 | 0.073 | 7.4 / 10.6 / 14.1 | Trapa pericarp [A] | 0.0561 at 20 C | Wrong tissue, and the authors are unknown. Not usable. |
| ALMOND | 0.010 / 7 / 0.92 | 0.042 | 4.1 / 6.3 / 8.8 | Taitano and Singh 2012 (not accessible); snippet 0.020-0.035 kg/kg [A, unverified] | none verified | Unverified. |
| COCONUT | 0.010 / 3.5 / 0.92 | 0.071 | 6.6 / 9.0 / 12.9 | green pulp, 30-70 C (out of range); copra (Henderson, no GAB) | none | No usable source. |
| HAZELNUT | 0.010 / 7 / 0.92 | 0.035 | 3.4 / 5.2 / 7.3 | Guzey 2001 [A], raw Tombul nuts, 25 C | Wm 0.0155 | The model total is more than double. |
| WALNUT | 0.010 / 7 / 0.92 | 0.031 | 3.0 / 4.7 / 6.6 | Togrul and Arslan 2007 [A] | BET monolayer, no GAB seen | No GAB numbers. |
| PISTACHIO | 0.010 / 7 / 0.92 | 0.040 | 3.8 / 5.9 / 8.3 | Soleimanifard and Hamdami 2018 (table garbled) | wm about 0.007-0.016 at 40-70 C, not trusted | Unverified. |
| PEANUT | 0.010 / 7 / 0.92 | 0.043 | 4.1 / 6.4 / 8.9 | Kaur 2023 [A], 25-45 C | monolayer 3.1-5.6 % | Consistent with the model (0.043). |
| BEECH_MAST | 0.010 / 7 / 0.92 | 0.039 | 3.8 / 5.8 / 8.1 | none | none | No source found (all hits were beech wood). |
| FLAXSEED | 0.005 / 8 / 0.91 | 0.054 | 5.3 / 8.0 / 10.9 | Valente 2020 [T], desorption at 40-80 C only | Xm 3.78 at 40 C | There is no 15-25 C data. |
| HEMP | 0.005 / 8 / 0.91 | 0.061 | 5.9 / 8.8 / 12.0 | none seen (PMC6145264 and a 2023 MDPI paper were not readable) | none | No numbers. |
| PUMPKIN_SEED | 0.005 / 8 / 0.91 | 0.038 | 3.7 / 5.7 / 7.8 | Mayor 2005 (GAB values not seen) | none | No numbers. |
| SUNFLOWER_SEED | 0.005 / 8 / 0.91 | 0.034 | 3.4 / 5.1 / 7.0 | none (Halsey and Smith fits only) | none | No numbers. |
| CHIA | 0.005 / 8 / 0.91 | 0.060 | 5.9 / 8.8 / 12.0 | Arslan-Tontul 2020 [A], whole seed | monolayer 2.39-2.91 g/100 g | The model monolayer is about double. |
| SESAME | 0.007 / 4.5 / 0.92 | 0.043 | 4.1 / 5.9 / 8.5 | Menkov and Durakova 2007 [T], semi-defatted flour, 25 C | EMC 7.16 / 8.96 / 9.32 % d.b. at aw 0.328 / 0.432 / 0.576 | The material is defatted, so it is not comparable with the model's full-fat composition. |

### Coverage summary
- 63 constants in total, in 40 independent families once whole/refined, maize and semolina variants are merged.
- With a verified, usable GAB source for the same material at roughly the right temperature: wheat flour, whole wheat, chestnut, teff, sorghum (outlier), fava bean (layout uncertain) = 6 constants.
- Verified for another aspect only (m0 only, EMC only, or other temperatures): semolina (2), sesame, flaxseed.
- Abstract/snippet only: hazelnut, peanut, chia, lupin, buckwheat, quinoa, barley, almond.
- No source found at all: khorasan (2), spelt (2), emmer (2), einkorn (2), triticale, rye (2), oat (2), maize (6), rice (2), millet (2), foxtail millet, amaranth, chickpea, lentil (2), pea (2), soy, white bean, black bean, acorn, water chestnut, coconut, walnut, pistachio, beech mast, hemp, pumpkin seed, sunflower seed.

## 3. Comparison with the test bands

The test's STARCHY band is 6.5-9.5 / 9.5-12.5 / 12-15 % w.b. at aw 0.3 / 0.5 / 0.65, and FIBROUS is 7-10 / 10-13 / 12.5-15.5.
- Wheat flour [C] 7.7 / 10.6 / 13.3 (ads.) and 8.8 / 11.3 / 13.7 (des.) are inside STARCHY.
- Whole wheat [C] 7.5 / 10.1 / 12.4 is inside the FIBROUS band at aw 0.3 and 0.5, but at aw 0.65 it is just below the lower bound (12.4 vs 12.5). The model's 14.9 at 0.65 is close to the upper edge of the band and above what the source supports.
- Teff [C] flour 6.2 / 8.5 / 10.6 is below the FIBROUS band at all three points (7 / 10 / 12.5). The seed 7.0 / 9.3 / 11.4 is also lower than FIBROUS. This suggests that teff belongs in a lower band than the other pseudocereals.
- Sorghum [C] 3.6 / 6.1 / 8.3 is far below STARCHY. This is probably a source artefact (see 1.4).
- Chestnut adsorption [C] 6.0 / 8.1 / 10.5 is below STARCHY at all three points. The desorption branch (7.0 / 9.3 / 12.0) is at the lower edge. The model (6.8 / 9.6 / 12.9) corresponds to desorption.
- Fava bean seed adsorption [C] 5.7 / 8.7 / 11.5 is below the PULSE band (6-8.5 / 8.5-10.5 / 11.5-13.5) at aw 0.3 and at its lower edge elsewhere. This is one source, with uncertain units.
- Hazelnut: Jung et al. (2018) [A] give about 4.2-5.0 % d.b. at aw 0.65, which is about 4.0-4.8 % w.b. The NUT band at 0.65 is 6-9 % w.b., so it looks too high. A single snippet source.
- Sesame: a measured 6.7 % w.b. at aw 0.328 sits above the NUT band (2.5-5). The material is defatted, so this is not a like-for-like comparison.
- The band comment in the test cites Moreira et al. 2010 for "starchy flours 7-9 % at 0.3, 10-12 % at 0.5 and 13-14.5 % at 0.65". That paper's printed 20 C wheat data support roughly 7.7-8.8 / 10.6-11.3 / 13.3-13.7, so the claim holds for wheat.
- The test `chestnutMonolayerMatchesMoreira2010` checks 0.059 +/- 0.007 g/g, which I could only see as a snippet from Chenlo et al. (desorption), not from Moreira et al. 2010. Moreira 2010 prints 0.049 (adsorption) and 0.062 (desorption) at 20 C. The test and the Javadoc attribute the value to Moreira et al., 2010, which should be rechecked.

## 4. Values that look like they should be corrected

Only changes with a source are proposed. All proposals use the model's convention, wmBase = Wm_total - 0.085*protein - 0.12*fiber, with the composition from the test.

1. **WHOLE_WHEAT: C 13.6, K 0.733** (Martin-Santos 2012 [T], 20 C). wmBase 0.057 is already consistent, since the source gives Wm 0.0801 and the model gives 0.080. With C 13.6 and K 0.733 the model gives 7.5 / 10.1 / 12.4 % w.b., matching the measured points. The FIBROUS band at 0.65 would need a lower bound of about 12 %. Confidence: high for this flour.
2. **WHEAT: K about 0.75, C about 9.7, wmBase about 0.075** (Moreira 2010 [T], adsorption: Wm 0.087 - 0.0102 - 0.0020). This changes the EMC little (about -0.5 to -0.9 point at aw 0.65). It is a refinement, and desorption would instead give C 20 / K 0.72 / wmBase 0.077. Lower K affects the high-aw end and the 0.2 d.b. clamp.
3. **TEFF: wmBase about 0.049, C 10.4, K 0.7255** (Abdissa 2020 [T] flour, 27 C, assuming Mo is % d.b.). The model is 2-5 points too wet at aw 0.3-0.65. The unit and the 27 C temperature are uncertain, and only one source supports this. Confidence: medium.
4. **CHESTNUT: re-document or adjust the branch.** The current model follows the desorption branch (Wm_tot 0.063 vs 0.062, EMC within 0.2-0.9 points). If the application is for adsorption (flour taking up moisture), the source supports Wm_tot 0.049 (wmBase about 0.036), C about 60, K about 0.90. Also correct the attribution of the 0.059 +/- 0.007 value, which is Chenlo et al., not Moreira et al. 2010 (see section 3).
5. **FAVA_BEAN and pulses: consider lower EMC.** A single seed source at 25 C (adsorption) gives about 1.3-1.5 points less than the model (5.7 / 8.7 / 11.5 vs 7.0 / 9.9 / 13.0). Do not change before the units and GAB assignment are checked in the Karatas and Arslan table.
6. **SORGHUM: do not adopt** Galdeano 2018 (C 2.76, K 0.714) without checking. It is an outlier that conflicts with all other cereals.
7. **Nuts and seeds, structural issue.** Published monolayers are far below the model's totals for hazelnut (0.0155 vs 0.035) and chia (0.024-0.029 vs 0.060), and flaxseed's published value is also lower. Because the additive terms 0.085*protein + 0.12*fiber alone give 0.0255 for hazelnut (protein 0.15, fiber 0.10) and about 0.054 for chia (protein 0.17, fiber 0.34), no non-negative wmBase can reproduce these. The additive protein/fiber formula looks unsuitable for fat-rich materials. Possible fixes are lower coefficients for them or a negative wmBase. Peanut (0.043 vs 3.1-5.6 %) is consistent. Confidence: medium, with snippet-level sources.
8. **SEMOLINA (and durum) wheat:** Erbas 2005 gives m0 11.8 g/100 g at 20 C, far above the model's 0.074. Verify the basis (and C and K) in the paper before changing anything.
9. **Unsupported analogies.** All khorasan, spelt, emmer, einkorn, triticale, rye, oat, maize, rice, millet, amaranth, chickpea, lentil, pea, soy, bean, lupin, acorn, coconut, walnut, pistachio, beech mast, hemp, pumpkin seed and sunflower seed values are not backed by any source seen in this study. They remain tuned approximations, as the `FlourType` Javadoc already states.

## 5. Best next steps for verification
1. Read the Ocieczek and Kostek papers (wheat and rye flour at 20 C) and Ocieczek 2012.
2. Read the Brett 2009 and PTPS 2020 papers (oat, rice, maize and millet flours).
3. Read the Erbas 2005 table (semolina C and K) and the Karatas and Arslan 2022 table (units).
4. Read the Arslan-Tontul 2020 chia table and the Guzey 2001 hazelnut table (C and K).
5. Confirm the Chenlo et al. citation and the 0.059 +/- 0.007 value.
