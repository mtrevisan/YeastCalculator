package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Encapsulates chemical and rheological specifications of a single flour instance.
 */
public record FlourInput(
	double strength,
	double plRatio,
	double sugar,
	double protein,
	double fat,
	double fiber,
	double ash,
	FlourType type){

	/**
	 * Estimated fraction of damaged starch granules [g/g_flour].
	 * <p>
	 * The damage is produced by milling and is higher in hard, strong flours: about 2% of the flour mass for a soft
	 * flour, growing with the strength W. It sets both the water absorption capacity and the starch the amylases can
	 * hydrolyze into fermentable sugars.
	 * </p>
	 */
	public double estimatedDamagedStarch(){
		return 0.02 + strength / 15000.;
	}

}
