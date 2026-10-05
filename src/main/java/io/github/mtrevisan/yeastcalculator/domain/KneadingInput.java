package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Describes mechanical parameters defining the matrix development profile.
 */
public class KneadingInput{

	public enum KneadingType{
		MANUAL(0.1, 0.4),
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


	private final KneadingType type;
	private final double durationMinutes;


	public KneadingInput(KneadingType type, double durationMinutes){
		this.type = type;
		this.durationMinutes = durationMinutes;
	}


	public KneadingType getType(){
		return type;
	}

	public double getDurationMinutes(){
		return durationMinutes;
	}

}
