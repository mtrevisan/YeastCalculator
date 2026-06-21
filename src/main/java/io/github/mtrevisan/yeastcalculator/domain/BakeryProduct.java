package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Defines the rheological targets, biochemical boundaries, and optimization penalty scales
 * specific to different core categories of baked goods.
 * <p>
 * Each product profile maps to a precise physical strategy, defining how far the gluten matrix
 * can safely stretch (volume limit) and how much residual carbohydrate substrate must be
 * preserved (sugar threshold) to feed the final baking phase and Maillard reaction.
 * </p>
 */
public enum BakeryProduct{
	/**
	 * NEAPOLITAN_PIZZA:
	 * Low volume in the tray to preserve extensibility for manual stretching.
	 * Lower sugar threshold because extreme baking temperatures (450 °C+)
	 * will cause flash charring if residual sugars are too high.
	 */
	NEAPOLITAN_PIZZA(1.7, 0.009, 25., 30.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeastDot){
			// Pizza targets a relaxed dough structure (~60% of max potential)
			final double targetVolume = maxPotential * 0.60;
			final double volumeError = Math.abs(targetVolume - finalVolume) * getTearingPenaltyMultiplier();

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);
			return volumeError + sugarPenalty;
		}
	},

	/**
	 * PAN_PIZZA_ROMAN:
	 * High hydration, large open alveoli. Maximum volume ceiling.
	 * Medium-high sugar required to sustain longer baking charts at 250 °C.
	 */
	ROMAN_PAN_PIZZA(2.2, 0.013, 15., 45.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeastDot){
			// Pan pizza targets higher volumetric development (~75% of max potential)
			final double targetVolume = maxPotential * 0.75;
			final double volumeError = Math.abs(targetVolume - finalVolume) * getTearingPenaltyMultiplier();

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);
			return volumeError + sugarPenalty;
		}
	},

	/**
	 * GASTRONOMY_TEGLIA:
	 * Spongy, soft, high-thickness crumb with small, uniform, dense bubble structures.
	 * Controlled volume cap to prevent cell walls from thinning and merging into caves.
	 * High-sugar residue target (1.5%) to feed the crumb during very long, gentle bake profiles (220 °C).
	 */
	GASTRONOMY_PAN_PIZZA(2., 0.015, 20., 50.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeastDot){
			// High-walled soft pan pizza targets high volume development (~80% of max potential)
			final double targetVolume = maxPotential * 0.80;
			final double volumeError = Math.abs(targetVolume - finalVolume) * getTearingPenaltyMultiplier();

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);
			return volumeError + sugarPenalty;
		}
	},

	/**
	 * BREAD:
	 * Balanced freestanding three-dimensional expansion profile.
	 * Requires structural elasticity reserves for scoring cuts and steam oven spring.
	 */
	BREAD(1.9, 0.012, 20., 40.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeastDot){
			// Bread structural targeting aims close to maximum expansion (~85% of max potential)
			final double targetVolume = maxPotential * 0.85;
			final double volumeError = Math.abs(targetVolume - finalVolume) * getTearingPenaltyMultiplier();

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);
			// Stability check for structural loaf layout
			final double stabilityPenalty = Math.abs(finalYeastDot) * 10.;
			return volumeError + sugarPenalty + stabilityPenalty;
		}
	};


	private final double glutenTearingLimit;
	private final double minSafeSugarThreshold;
	private final double tearingPenaltyMultiplier;
	private final double starvationPenaltyMultiplier;


	BakeryProduct(final double vLimit, final double sLimit, final double pTear, final double pStarve){
		this.glutenTearingLimit = vLimit;
		this.minSafeSugarThreshold = sLimit;
		this.tearingPenaltyMultiplier = pTear;
		this.starvationPenaltyMultiplier = pStarve;
	}


	public double getGlutenTearingLimit(){
		return glutenTearingLimit;
	}

	public double getMinSafeSugarThreshold(){
		return minSafeSugarThreshold;
	}

	public double getTearingPenaltyMultiplier(){
		return tearingPenaltyMultiplier;
	}

	public double getStarvationPenaltyMultiplier(){
		return starvationPenaltyMultiplier;
	}


	/**
	 * Computes the custom cost fitness value based on the final dough parameters.
	 * * @param finalVolume  Simulated gas volume at the timeline end (mL/g_flour).
	 * @param maxPotential Maximum theoretical expansion potential of the gluten network (mL/g_flour).
	 * @param residualSugar Leftover fermentable sugar percentage at execution end (g/g_flour).
	 * @param finalYeastDot Instantaneous final gas derivative (speed of volume growth).
	 * @return Fitness score. The closer to zero, the more optimal the yeast dosage is.
	 */
	public abstract double computeFitness(double finalVolume, double maxPotential, double residualSugar,
		double finalYeastDot);

}
