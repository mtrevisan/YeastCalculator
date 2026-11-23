package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Describes mechanical parameters defining the matrix development profile.
 * <p>
 * The autolyse is the rest of the flour and the water alone (before adding the yeast and the salt): the flour hydrates
 * and the gluten starts to develop with no work, so fewer folds are needed afterwards. It is not part of the
 * fermentation timeline, which starts when the yeast is added.
 * </p>
 *
 * @param autolyseMinutes	Autolyse duration [min], 0 for none.
 */
public record KneadingInput(
	KneadingType type,
	double durationMinutes,
	double autolyseMinutes){

	public KneadingInput{
		if(autolyseMinutes < 0.)
			throw new IllegalArgumentException("Autolyse duration cannot be negative");
	}

	/**
	 * Kneading with no autolyse.
	 */
	public KneadingInput(final KneadingType type, final double durationMinutes){
		this(type, durationMinutes, 0.);
	}



	public enum KneadingType{
		MANUAL(0.1, 0.6),
		SPIRAL_MIXER(1., 1.),
		FORK_MIXER(0.6, 0.75),
		DIVING_ARMS(0.8, 1.1);

		private final double frictionModifier;
		private final double developmentEfficiency;

		KneadingType(final double frictionModifier, final double developmentEfficiency){
			this.frictionModifier = frictionModifier;
			this.developmentEfficiency = developmentEfficiency;
		}

		public double getFrictionModifier(){
			return frictionModifier;
		}

		public double getDevelopmentEfficiency(){
			return developmentEfficiency;
		}
	}

}
