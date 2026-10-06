package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Defines the rheological targets, biochemical boundaries, and optimization penalty scales
 * specific to different core categories of baked goods.
 * <p>
 * Each product profile maps to a precise physical strategy, defining how far the gluten matrix
 * can safely stretch (volume limit), how much residual carbohydrate substrate must be
 * preserved (sugar threshold), and the maximum sensory threshold for residual biological biomass.
 * </p>
 */
public enum BakeryProduct{

	/**
	 * NEAPOLITAN_PIZZA:
	 * Targets a relaxed, highly extensible gluten structure to prevent elastic snapback during manual stretching.
	 * Can tolerate lower total gas volume (needs to stay within a flexible, workable envelope).
	 * Requires strict retention of residual sugars to prevent flash charring in high-temperature ovens (450 °C+).
	 * Strict yeast cell density caps to ensure a clean flavor profile in thin crusts.
	 */
	NEAPOLITAN_PIZZA(1.7, 0.009, 0.001, 25., 30.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			// Reconstruct the Bloksma polymer retention efficiency from DoughOdeSystem (kTear = 1.5)
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			// Neapolitan pizza targets an extensible, relaxed state (approx. 40% - 60% of max potential volume)
			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.92)
				// Heavily penalize if it enters the gas leakage phase (loses extensibility, tears during stretching)
				rheologicalPenalty += getTearingPenaltyMultiplier() * 120. * (0.92 - retentionEfficiency);
			else if(finalVolume < (maxPotential * 0.35))
				// Penalize if underproofed (dough will be too elastic and spring back)
				rheologicalPenalty += getTearingPenaltyMultiplier() * 15. * (maxPotential * 0.35 - finalVolume);

			// Biochemical check: Protect sugars for high-heat browning control
			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1200. * (getMinSafeSugarThreshold() - residualSugar);

			// Sensory check: Avoid yeasty off-flavors
			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 150.;

			// Kinetic check: Gas production rate should be stable and entering plateau
			final double kineticPenalty = Math.abs(finalYeastDot) * 8.;

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty;
		}
	},

	/**
	 * ROMAN_PAN_PIZZA:
	 * High hydration, large open alveoli structure. Pushes close to the structural threshold limit.
	 * Requires maximizing volume while maintaining a stable gluten mesh to support heavy aeration.
	 * Moderate yeast limit to tolerate intensive expansion without crumbling or structural collapse.
	 */
	ROMAN_PAN_PIZZA(2.2, 0.013, 0.002, 15., 45.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			// Roman pan pizza targets maximum structural development (approx. 65% - 80% of max potential volume)
			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.82)
				// Severe overproofing penalty: gas bubbles are merging, walls are thinning out, collapse imminent
				rheologicalPenalty += getTearingPenaltyMultiplier() * 150. * (0.82 - retentionEfficiency);
			else if(finalVolume < (maxPotential * 0.55))
				// Underproofing penalty: structure is dense, missing the classic open alveolar expansion
				rheologicalPenalty += getTearingPenaltyMultiplier() * 25. * (maxPotential * 0.55 - finalVolume);

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);

			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 100.;

			final double kineticPenalty = Math.abs(finalYeastDot) * 5.;

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty;
		}
	},

	/**
	 * GASTRONOMY_PAN_PIZZA:
	 * Spongy, highly resilient, fine crumb layout.
	 * Must optimize expansion to keep bubble structures uniform and prevent cells from merging into caverns.
	 * High sugar buffer target to sustain long baking charts at 220 °C.
	 * Higher yeast threshold tolerated due to the protective emulsifying effects of added oils/fats.
	 */
	GASTRONOMY_PAN_PIZZA(2., 0.015, 0.003, 20., 50.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			// Soft pan pizza targets highly uniform gas retention (approx. 60% - 75% of max potential volume)
			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.85)
				// Overproofed: Gluten walls are tearing, causing non-uniform large air pockets
				rheologicalPenalty += getTearingPenaltyMultiplier() * 130. * (0.85 - retentionEfficiency);
			else if(finalVolume < (maxPotential * 0.45))
				// Underproofed: The crumb will be overly dense, heavy, and lack soft sponginess
				rheologicalPenalty += getTearingPenaltyMultiplier() * 20. * (maxPotential * 0.45 - finalVolume);

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);

			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 100.;

			// Stability: Ensure the gas velocity curve has completely flattened out (plateau phase)
			final double kineticPenalty = Math.abs(finalYeastDot) * 6.;

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty;
		}
	},

	/**
	 * BREAD:
	 * Freestanding, three-dimensional balanced gas cell layout.
	 * Must maintain high elastic retention reserves to handle oven-spring expansions and scoring expansions.
	 * Strict kinetic constraints to prevent pocket channeling and uneven baking tunnels.
	 */
	BREAD(1.9, 0.012, 0.0015, 20., 40.){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			// Freestanding bread loaves target high retention reserves (approx. 70% - 82% of max potential volume)
			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.88)
				// Overproofed: Loaf will deflate or flat-line during scoring cuts or when hitting steam in the oven
				rheologicalPenalty += getTearingPenaltyMultiplier() * 140. * (0.88 - retentionEfficiency);
			else if(finalVolume < (maxPotential * 0.50))
				// Underproofed: Missing structural volume potential; dense core crumb expected
				rheologicalPenalty += getTearingPenaltyMultiplier() * 30. * (maxPotential * 0.50 - finalVolume);

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);

			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 100.;

			// Dynamic structural stabilization check for freestanding loaves
			final double stabilityPenalty = Math.abs(finalYeastDot) * 12.;

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + stabilityPenalty;
		}
	};


	private final double glutenTearingLimit;
	private final double minSafeSugarThreshold;
	private final double maxAllowedFinalYeast;
	private final double tearingPenaltyMultiplier;
	private final double starvationPenaltyMultiplier;


	BakeryProduct(final double glutenTearingLimit, final double minSafeSugarThreshold, final double maxAllowedFinalYeast,
			final double tearingPenaltyMultiplier, final double starvationPenaltyMultiplier){
		this.glutenTearingLimit = glutenTearingLimit;
		this.minSafeSugarThreshold = minSafeSugarThreshold;
		this.maxAllowedFinalYeast = maxAllowedFinalYeast;
		this.tearingPenaltyMultiplier = tearingPenaltyMultiplier;
		this.starvationPenaltyMultiplier = starvationPenaltyMultiplier;
	}


	public double getGlutenTearingLimit(){
		return glutenTearingLimit;
	}

	public double getMinSafeSugarThreshold(){
		return minSafeSugarThreshold;
	}

	public double getMaxAllowedFinalYeast(){
		return maxAllowedFinalYeast;
	}

	public double getTearingPenaltyMultiplier(){
		return tearingPenaltyMultiplier;
	}

	public double getStarvationPenaltyMultiplier(){
		return starvationPenaltyMultiplier;
	}


	/**
	 * Computes the custom cost fitness value based on the final dough parameters.
	 * @param finalVolume	Simulated gas volume at the timeline end (mL/g_flour).
	 * @param maxPotential	Maximum theoretical expansion potential of the gluten network (mL/g_flour).
	 * @param residualSugar	Leftover fermentable sugar percentage at execution end (g/g_flour).
	 * @param finalYeast	Instantaneous active yeast cell concentration at the execution end (g/g_flour).
	 * @param finalYeastDot	Instantaneous final gas derivative (speed of volume growth).
	 * @return	Fitness score. The closer to zero, the more optimal the yeast dosage is.
	 */
	public abstract double computeFitness(double finalVolume, double maxPotential, double residualSugar,
		double finalYeast, double finalYeastDot);

}
