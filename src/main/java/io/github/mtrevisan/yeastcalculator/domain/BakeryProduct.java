package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Defines the rheological targets, biochemical boundaries, and optimization penalty scales
 * specific to different core categories of baked goods.
 * <p>
 * Each product profile maps to a precise physical strategy, defining how far the gluten matrix
 * can safely stretch (volume limit), how much residual carbohydrate substrate must be
 * preserved (sugar threshold), and the maximum sensory threshold for residual biological biomass.
 * </p>
 * <p>
 * NEW: fold and oil rewards are now Gaussian functions centered on a product-specific optimum.
 * The previous linear reward was monotonically increasing, which drove the optimizer to always
 * select the maximum number of folds and the maximum amount of oil regardless of culinary
 * reality. A Gaussian reward peaks at the optimal value and decreases on both sides, which is
 * the correct shape for a "diminishing returns" contribution.
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
	NEAPOLITAN_PIZZA(1.7, 0.003, 0.001, 25., 30., 4.0, 1.5, 0.35, 15., 1, 0.01){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot, final double oilRatio, final int appliedFolds,
				final double foldDegassingRetention){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double foldDegassingScale = Math.pow(foldDegassingRetention, appliedFolds);
			final double effectiveUnderproofThreshold = maxPotential * getUnderproofThresholdFactor() * foldDegassingScale;

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.92)
				// Heavily penalize if it enters the gas leakage phase (loses extensibility, tears during stretching)
				rheologicalPenalty += getTearingPenaltyMultiplier() * 120. * (0.92 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Penalize if underproofed (dough will be too elastic and spring back)
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);

			// Biochemical check: Protect sugars for high-heat browning control
			// NOTE: threshold lowered from 0.009 to 0.003 (0.3% residual sugar) because a direct
			// dough fermented at 28-30 °C for 4-8 h retains only 0.1-0.3 g of fermentable sugars
			// per 100 g of flour (Jayaram et al., 2013, Food Chemistry 136: 301-308). The previous
			// threshold of 0.9% was physically unreachable in a standard recipe without added sucrose.
			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1200. * (getMinSafeSugarThreshold() - residualSugar);

			// Sensory check: Avoid yeasty off-flavors
			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 150.;

			// Kinetic check: Gas production rate should be stable and entering plateau
			final double kineticPenalty = Math.abs(finalYeastDot) * 8.;

			// Gaussian reward: peaks at the product-specific optimal fold count.
			final double deviation = appliedFolds - getOptimalFoldCount();
			final double sigma = (deviation < 0? 0.8 : 1.5);
			final double foldReward = -getFoldRewardMultiplier() * Math.exp(-deviation * deviation / (2. * sigma * sigma));

			// Gaussian reward: peaks at the product-specific optimal oil fraction.
			final double oilDeviation = (oilRatio - getOptimalOilRatio()) / OIL_REWARD_SIGMA;
			final double oilReward = -getOilRewardMultiplier() * Math.exp(-oilDeviation * oilDeviation / 2.);

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty + oilReward + foldReward;
		}

		@Override
		public HandlingInstruction getHandlingInstructions(int setNumber, int totalSets, double absoluteMinute){
			String title = "EXTENSIBILITY ALIGNMENT PHASE";
			String description = "Pull the dough edge upward with moderate tension and fold it over the center. "
				+ "Rotate 90° and repeat for all 4 sides. Work quickly to establish a baseline shape while "
				+ "keeping the gluten highly extensible for manual stretching later.";
			return new HandlingInstruction(title, description, "MODERATE");
		}
	},

	/**
	 * ROMAN_PAN_PIZZA:
	 * High hydration, large open alveoli structure. Pushes close to the structural threshold limit.
	 * Requires maximizing volume while maintaining a stable gluten mesh to support heavy aeration.
	 * Moderate yeast limit to tolerate intensive expansion without crumbling or structural collapse.
	 */
	ROMAN_PAN_PIZZA(2.2, 0.004, 0.002, 15., 45., 10.0, 2.5, 0.50, 15., 4, 0.025){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot, final double oilRatio, final int appliedFolds,
				final double foldDegassingRetention){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double foldDegassingScale = Math.pow(foldDegassingRetention, appliedFolds);
			final double effectiveUnderproofThreshold = maxPotential * getUnderproofThresholdFactor() * foldDegassingScale;

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.82)
				// Severe overproofing penalty: gas bubbles are merging, walls are thinning out, collapse imminent
				rheologicalPenalty += getTearingPenaltyMultiplier() * 150. * (0.82 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Underproofing penalty: structure is dense, missing the classic open alveolar expansion
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);

			// NOTE: threshold lowered from 0.013 to 0.004. Roman pan pizza is typically fermented
			// for 18-24 h at room temperature, so most fermentable sugars are consumed; a residual
			// of 0.3-0.4 g/100 g flour is the realistic ceiling (Jayaram et al., 2014, JAFC 62: 10312).
			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);

			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 100.;

			final double kineticPenalty = Math.abs(finalYeastDot) * 5.;

			final double foldDeviation = appliedFolds - getOptimalFoldCount();
			final double foldReward = -getFoldRewardMultiplier() * Math.exp(-foldDeviation * foldDeviation / 2.);

			final double oilDeviation = (oilRatio - getOptimalOilRatio()) / OIL_REWARD_SIGMA;
			final double oilReward = -getOilRewardMultiplier() * Math.exp(-oilDeviation * oilDeviation / 2.);

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty + oilReward + foldReward;
		}

		@Override
		public HandlingInstruction getHandlingInstructions(int setNumber, int totalSets, double absoluteMinute){
			if(setNumber <= 2)
				return new HandlingInstruction(
					"HIGH-TENSION COIL STRETCH",
					"Perform vigorous coil folds. Wet your hands, lift the entire dough mass from the center, "
						+ "allow it to stretch downward under its own weight, and tuck the ends underneath. "
						+ "Essential for driving vertical strength into highly hydrated dough networks.",
					"HIGH"
				);

			return new HandlingInstruction(
				"ALVEOLAR PRESERVATION COIL",
				"Gently slide hands under the dough and perform a light coil fold. Do not stretch to the breaking point. "
					+ "The cell walls are thinning; excessive force will pop large air pockets and cause gas weeping.",
				"GENTLE"
			);
		}
	},

	/**
	 * GASTRONOMY_PAN_PIZZA:
	 * Spongy, highly resilient, fine crumb layout.
	 * Must optimize expansion to keep bubble structures uniform and prevent cells from merging into caverns.
	 * High sugar buffer target to sustain long baking charts at 220 °C.
	 * Higher yeast threshold tolerated due to the protective emulsifying effects of added oils/fats.
	 */
	GASTRONOMY_PAN_PIZZA(2., 0.004, 0.003, 20., 50., 12.0, 8.0, 0.35, 8., 3, 0.030){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot, final double oilRatio, final int appliedFolds,
				final double foldDegassingRetention){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double foldDegassingScale = Math.pow(foldDegassingRetention, appliedFolds);
			final double effectiveUnderproofThreshold = maxPotential * getUnderproofThresholdFactor() * foldDegassingScale;

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.85)
				// Overproofed: Gluten walls are tearing, causing non-uniform large air pockets
				rheologicalPenalty += getTearingPenaltyMultiplier() * 130. * (0.85 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Underproofed: The crumb will be overly dense, heavy, and lack soft sponginess
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);

			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);

			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 100.;

			// Stability: Ensure the gas velocity curve has completely flattened out (plateau phase)
			final double kineticPenalty = Math.abs(finalYeastDot) * 6.;

			// Gaussian reward: peaks at 3 fold sets (Hamelman, Bread, 2nd ed.; Reinhart,
			// The Bread Baker's Apprentice). Both fewer and more folds reduce the reward.
			final double foldDeviation = appliedFolds - getOptimalFoldCount();
			final double foldReward = -getFoldRewardMultiplier() * Math.exp(-foldDeviation * foldDeviation / 2.);

			// Gaussian reward: peaks at 3% oil by flour mass, the typical upper target for a
			// gastronomy pan pizza where fat contributes tenderness and crust browning.
			final double oilDeviation = (oilRatio - getOptimalOilRatio()) / OIL_REWARD_SIGMA;
			final double oilReward = -getOilRewardMultiplier() * Math.exp(-oilDeviation * oilDeviation / 2.);

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty + oilReward + foldReward;
		}

		@Override
		public HandlingInstruction getHandlingInstructions(int setNumber, int totalSets, double absoluteMinute){
			if(setNumber == 1)
				return new HandlingInstruction(
					"STRUCTURAL ALIGNMENT PHASE",
					"Slide wet hands underneath the dough, pull firmly straight upward to maximize mechanical elongation, "
						+ "and fold it over the center. Repeat for all 4 quadrants. This forces the chaotic protein clumps "
						+ "to re-align into long, parallel sheets, building early structural tenacity (P).",
					"HIGH"
				);

			return new HandlingInstruction(
				"GAS POCKET SUBDIVISION PHASE",
				"Execute with extreme gentleness. Dip fingertips in water, slide them under the edge, lift slightly until "
					+ "you feel baseline resistance, and flip smoothly. Do not punch or stretch aggressively. This splits "
					+ "coarse air bubbles into thousands of micro-cavities, securing a sponge-like, soft crumb distribution.",
				"GENTLE"
			);
		}
	},

	/**
	 * BREAD:
	 * Freestanding, three-dimensional balanced gas cell layout.
	 * Must maintain high elastic retention reserves to handle oven-spring expansions and scoring expansions.
	 * Strict kinetic constraints to prevent pocket channeling and uneven baking tunnels.
	 */
	BREAD(1.9, 0.004, 0.0015, 20., 40., 8.0, 3.0, 0.40, 15., 3, 0.015){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double residualSugar,
				final double finalYeast, final double finalYeastDot, final double oilRatio, final int appliedFolds,
				final double foldDegassingRetention){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double foldDegassingScale = Math.pow(foldDegassingRetention, appliedFolds);
			final double effectiveUnderproofThreshold = maxPotential * getUnderproofThresholdFactor() * foldDegassingScale;

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.88)
				// Overproofed: Loaf will deflate or flat-line during scoring cuts or when hitting steam in the oven
				rheologicalPenalty += getTearingPenaltyMultiplier() * 140. * (0.88 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Underproofed: Missing structural volume potential; dense core crumb expected
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);

			// NOTE: threshold lowered from 0.012 to 0.004. For bread fermented for 4-6 h the residual
			// sugar is below 0.4 g/100 g flour unless the recipe contains added sucrose or honey.
			double sugarPenalty = 0.;
			if(residualSugar < getMinSafeSugarThreshold())
				sugarPenalty = getStarvationPenaltyMultiplier() * 1000. * (getMinSafeSugarThreshold() - residualSugar);

			double offFlavorPenalty = 0.;
			if(finalYeast > getMaxAllowedFinalYeast())
				offFlavorPenalty = Math.pow((finalYeast - getMaxAllowedFinalYeast()) * 1000., 2.) * 100.;

			// Dynamic structural stabilization check for freestanding loaves
			final double stabilityPenalty = Math.abs(finalYeastDot) * 12.;

			final double foldDeviation = appliedFolds - getOptimalFoldCount();
			final double foldReward = -getFoldRewardMultiplier() * Math.exp(-foldDeviation * foldDeviation / 2.);

			final double oilDeviation = (oilRatio - getOptimalOilRatio()) / OIL_REWARD_SIGMA;
			final double oilReward = -getOilRewardMultiplier() * Math.exp(-oilDeviation * oilDeviation / 2.);

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + stabilityPenalty + oilReward + foldReward;
		}

		@Override
		public HandlingInstruction getHandlingInstructions(int setNumber, int totalSets, double absoluteMinute){
			return new HandlingInstruction(
				"TENACITY OVERLAP FOLD",
				"Lift and fold the dough tightly onto itself. Aim to generate strong surface tension across the skin. "
					+ "This structural elasticity creates the pressure reservoir required to support clean scoring cuts "
					+ "and three-dimensional oven-spring without flattening out.",
				"HIGH"
			);
		}
	};


	// Width of the Gaussian oil reward, expressed in oil fraction units. A deviation of one
	// sigma from the optimum reduces the reward to 60.6% of its peak.
	private static final double OIL_REWARD_SIGMA = 0.02;


	private final double glutenTearingLimit;
	private final double minSafeSugarThreshold;
	private final double maxAllowedFinalYeast;
	private final double tearingPenaltyMultiplier;
	private final double starvationPenaltyMultiplier;
	private final double oilRewardMultiplier;
	private final double foldRewardMultiplier;
	private final double underproofThresholdFactor;
	private final double underproofPenaltyCoefficient;
	private final int optimalFoldCount;
	private final double optimalOilRatio;


	BakeryProduct(final double glutenTearingLimit, final double minSafeSugarThreshold, final double maxAllowedFinalYeast,
			final double tearingPenaltyMultiplier, final double starvationPenaltyMultiplier,
			final double oilRewardMultiplier, final double foldRewardMultiplier,
			final double underproofThresholdFactor, final double underproofPenaltyCoefficient,
			final int optimalFoldCount, final double optimalOilRatio){
		this.glutenTearingLimit = glutenTearingLimit;
		this.minSafeSugarThreshold = minSafeSugarThreshold;
		this.maxAllowedFinalYeast = maxAllowedFinalYeast;
		this.tearingPenaltyMultiplier = tearingPenaltyMultiplier;
		this.starvationPenaltyMultiplier = starvationPenaltyMultiplier;
		this.oilRewardMultiplier = oilRewardMultiplier;
		this.foldRewardMultiplier = foldRewardMultiplier;
		this.underproofThresholdFactor = underproofThresholdFactor;
		this.underproofPenaltyCoefficient = underproofPenaltyCoefficient;
		this.optimalFoldCount = optimalFoldCount;
		this.optimalOilRatio = optimalOilRatio;
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

	public double getOilRewardMultiplier(){
		return oilRewardMultiplier;
	}

	public double getFoldRewardMultiplier(){
		return foldRewardMultiplier;
	}

	public double getUnderproofThresholdFactor(){
		return underproofThresholdFactor;
	}

	public double getUnderproofPenaltyCoefficient(){
		return underproofPenaltyCoefficient;
	}

	public int getOptimalFoldCount(){
		return optimalFoldCount;
	}

	public double getOptimalOilRatio(){
		return optimalOilRatio;
	}


	/**
	 * Computes the custom cost fitness value based on the final dough parameters.
	 * @param finalVolume	Simulated gas volume at the timeline end (mL/g_flour).
	 * @param maxPotential	Maximum theoretical expansion potential of the gluten network (mL/g_flour).
	 * @param residualSugar	Leftover fermentable sugar percentage at execution end (g/g_flour).
	 * @param finalYeast	Instantaneous active yeast cell concentration at the execution end (g/g_flour).
	 * @param finalYeastDot	Instantaneous final gas derivative (speed of volume growth).
	 * @param oilRatio	Added oil ratio in the recipe (g/g_flour).
	 * @param appliedFolds	Number of stretch-and-fold sets actually applied within the timeline.
	 * @param foldDegassingRetention	Fraction of gas retained through each fold event.
	 * @return	Fitness score. The closer to negative infinity, the more optimal the recipe is.
	 */
	public abstract double computeFitness(double finalVolume, double maxPotential, double residualSugar,
		double finalYeast, double finalYeastDot, double oilRatio, int appliedFolds,
		double foldDegassingRetention);

	/**
	 * Generates specialized handling text advice mapped to the current process stage context.
	 */
	public abstract HandlingInstruction getHandlingInstructions(int setNumber, int totalSets, double absoluteMinute);

	/**
	 * Data transfer structure wrapping physical execution metrics.
	 */
	public record HandlingInstruction(
		String title,
		String instructions,
		// GENTLE, MODERATE, HIGH
		String intensity){}

}
