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
 * Defines the rheological targets, biochemical boundaries, and optimization penalty scales
 * specific to different core categories of baked goods.
 * <p>
 * Each product profile maps to a precise physical strategy, defining how far the gluten matrix
 * can safely stretch (volume limit), how much residual carbohydrate substrate must be
 * preserved (sugar threshold), and the maximum sensory threshold for residual biological biomass.
 * </p>
 * <p>
 * The number of stretch-and-fold
 * sets has no reward of its own: it emerges from the gluten development the product requires (see
 * {@link #computeStructuralCost(double, double, int)}).
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
	NEAPOLITAN_PIZZA(1.7, 0.003, 0.001, 25., 30., 0.35, 15., 0.60,
		new RecipeEnvelope(0.58, 0.65, 0.025, 0.030, 0., 0.005, 0., 0.02)){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double referencePotential,
				final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			double rheologicalPenalty = getRheologicalPenalty(finalVolume, maxPotential, referencePotential);

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

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty;
		}

		private double getRheologicalPenalty(final double finalVolume, final double maxPotential,
				final double referencePotential){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double effectiveUnderproofThreshold = referencePotential * getUnderproofThresholdFactor();

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.92)
				// Heavily penalize if it enters the gas leakage phase (loses extensibility, tears during stretching)
				rheologicalPenalty += getTearingPenaltyMultiplier() * 120. * (0.92 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Penalize if underproofed (dough will be too elastic and spring back)
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);
			return rheologicalPenalty;
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
	ROMAN_PAN_PIZZA(2.2, 0.004, 0.002, 15., 45., 0.50, 15., 0.85,
		new RecipeEnvelope(0.70, 0.85, 0.018, 0.024, 0., 0.010, 0.015, 0.04)){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double referencePotential,
				final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			double rheologicalPenalty = getRheologicalPenalty(finalVolume, maxPotential, referencePotential);

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

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + kineticPenalty;
		}

		private double getRheologicalPenalty(final double finalVolume, final double maxPotential,
				final double referencePotential){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double effectiveUnderproofThreshold = referencePotential * getUnderproofThresholdFactor();

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.82)
				// Severe overproofing penalty: gas bubbles are merging, walls are thinning out, collapse imminent
				rheologicalPenalty += getTearingPenaltyMultiplier() * 150. * (0.82 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Underproofing penalty: structure is dense, missing the classic open alveolar expansion
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);
			return rheologicalPenalty;
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
	GASTRONOMY_PAN_PIZZA(2., 0.004, 0.003, 20., 50., 0.35, 8., 0.80,
		new RecipeEnvelope(0.60, 0.75, 0.018, 0.024, 0., 0.010, 0.02, 0.05)){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double referencePotential,
				final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			double rheologicalPenalty = getRheologicalPenalty(finalVolume, maxPotential, referencePotential);

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

		private double getRheologicalPenalty(final double finalVolume, final double maxPotential,
				final double referencePotential){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double effectiveUnderproofThreshold = referencePotential * getUnderproofThresholdFactor();

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.85)
				// Overproofed: Gluten walls are tearing, causing non-uniform large air pockets
				rheologicalPenalty += getTearingPenaltyMultiplier() * 130. * (0.85 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Underproofed: The crumb will be overly dense, heavy, and lack soft sponginess
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);
			return rheologicalPenalty;
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
	BREAD(1.9, 0.004, 0.004, 20., 40., 0.40, 15., 0.90,
		new RecipeEnvelope(0.62, 0.78, 0.018, 0.022, 0., 0.010, 0., 0.03)){
		@Override
		public double computeFitness(final double finalVolume, final double maxPotential, final double referencePotential,
				final double residualSugar,
				final double finalYeast, final double finalYeastDot){
			double rheologicalPenalty = getRheologicalPenalty(finalVolume, maxPotential, referencePotential);

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

			return rheologicalPenalty + sugarPenalty + offFlavorPenalty + stabilityPenalty;
		}

		private double getRheologicalPenalty(final double finalVolume, final double maxPotential,
				final double referencePotential){
			final double kTear = 1.5;
			final double retentionEfficiency = 1. / (1. + Math.exp(kTear * (finalVolume - maxPotential)));

			final double effectiveUnderproofThreshold = referencePotential * getUnderproofThresholdFactor();

			double rheologicalPenalty = 0.;
			if(retentionEfficiency < 0.88)
				// Overproofed: Loaf will deflate or flat-line during scoring cuts or when hitting steam in the oven
				rheologicalPenalty += getTearingPenaltyMultiplier() * 140. * (0.88 - retentionEfficiency);
			else if(finalVolume < effectiveUnderproofThreshold)
				// Underproofed: Missing structural volume potential; dense core crumb expected
				rheologicalPenalty += getUnderproofPenaltyCoefficient()
					* (effectiveUnderproofThreshold - finalVolume);
			return rheologicalPenalty;
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



	// Weight of the formulation cost of hydration, salt and malt, and distance outside the envelope at which each
	// term reaches its weight
	private static final double HYDRATION_COST = 8.;
	private static final double HYDRATION_TOLERANCE = 0.04;
	private static final double SALT_COST = 6.;
	private static final double SALT_TOLERANCE = 0.003;
	private static final double MALT_COST = 4.;
	private static final double MALT_TOLERANCE = 0.004;
	// Oil has no point optimum: it is free inside the range the style is made with (fat tenderizes the crumb and
	// lengthens the shelf life, but slows the yeast and, above 8% with the flour lipids, tears the network)
	private static final double OIL_COST = 5.;
	private static final double OIL_TOLERANCE = 0.01;
	// Gluten development: gap from the requirement at which the cost reaches its weight, and the margin above the
	// requirement before the dough starts to be over-developed
	private static final double DEVELOPMENT_TOLERANCE = 0.1;
	private static final double DEVELOPMENT_MARGIN = 0.15;
	private static final double UNDERDEVELOPMENT_COST = 6.;
	private static final double OVERDEVELOPMENT_COST = 3.;
	// Cost of each stretch-and-fold set (work, the gas lost and the dough temperature)
	private static final double FOLD_LABOR_COST = 0.5;
	// Change of the required development per unit of hydration above the middle of the product's envelope
	private static final double HYDRATION_DEVELOPMENT_SLOPE = 1.5;
	// Weak pull toward the middle of the envelope, as a fraction of the weight (1 = at the envelope edge)
	private static final double CENTER_PULL = 0.05;
	// Stronger pull for the hydration: nothing in the fermentation rewards more water (the gluten development required
	// grows with it), so without it the hydration would always sit on the lowest edge of the envelope, while a baker
	// hydrates the dough toward the middle of the range the flour and the style allow
	private static final double HYDRATION_CENTER_PULL = 0.5;


	private final double glutenTearingLimit;
	private final double minSafeSugarThreshold;
	private final double maxAllowedFinalYeast;
	private final double tearingPenaltyMultiplier;
	private final double starvationPenaltyMultiplier;
	private final double underproofThresholdFactor;
	private final double underproofPenaltyCoefficient;
	private final double requiredDevelopment;
	private final RecipeEnvelope recipeEnvelope;


	BakeryProduct(final double glutenTearingLimit, final double minSafeSugarThreshold, final double maxAllowedFinalYeast,
			final double tearingPenaltyMultiplier, final double starvationPenaltyMultiplier,
			final double underproofThresholdFactor, final double underproofPenaltyCoefficient,
			final double requiredDevelopment, final RecipeEnvelope recipeEnvelope){
		this.glutenTearingLimit = glutenTearingLimit;
		this.minSafeSugarThreshold = minSafeSugarThreshold;
		this.maxAllowedFinalYeast = maxAllowedFinalYeast;
		this.tearingPenaltyMultiplier = tearingPenaltyMultiplier;
		this.starvationPenaltyMultiplier = starvationPenaltyMultiplier;
		this.underproofThresholdFactor = underproofThresholdFactor;
		this.underproofPenaltyCoefficient = underproofPenaltyCoefficient;
		this.requiredDevelopment = requiredDevelopment;
		this.recipeEnvelope = recipeEnvelope;
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



	public double getUnderproofThresholdFactor(){
		return underproofThresholdFactor;
	}

	public double getUnderproofPenaltyCoefficient(){
		return underproofPenaltyCoefficient;
	}


	public double getRequiredDevelopment(){
		return requiredDevelopment;
	}

	/**
	 * Gluten development this product requires for a given hydration, in [0, 1].
	 * <p>
	 * The product value refers to the middle of its hydration envelope; a wetter dough has less intrinsic structure
	 * and needs more development to hold its shape, a drier one less.
	 * </p>
	 */
	public double requiredDevelopment(final double hydration){
		final double hydrationMid = (recipeEnvelope.hydrationMin() + recipeEnvelope.hydrationMax()) / 2.;
		return Math.clamp(requiredDevelopment + HYDRATION_DEVELOPMENT_SLOPE * (hydration - hydrationMid), 0.3, 0.95);
	}

	/**
	 * Structural cost of the gluten development reached by the mixer and the stretch &amp; folds.
	 * <p>
	 * It replaces any fixed "best number of folds": that number emerges from the development the product requires
	 * and the development already reached by kneading. An under-developed dough is slack and holds the gas poorly
	 * (quadratic cost below the requirement), an over-developed one is tight and tears when stretched (milder cost
	 * above the requirement plus a margin), and every fold costs a little (work, and the gas it degasses).
	 * </p>
	 *
	 * @param hydration	Water ratio (g/g_flour).
	 * @param developmentIndex	Gluten development reached at the end of the timeline, in [0, 1].
	 * @param appliedFolds	Number of stretch-and-fold sets applied within the timeline.
	 * @return	Non-negative cost, to be added to the fitness (lower is better).
	 */
	public double computeStructuralCost(final double hydration, final double developmentIndex, final int appliedFolds){
		final double required = requiredDevelopment(hydration);
		final double shortfall = Math.max(0., required - developmentIndex) / DEVELOPMENT_TOLERANCE;
		final double excess = Math.max(0., developmentIndex - required - DEVELOPMENT_MARGIN) / DEVELOPMENT_TOLERANCE;
		return UNDERDEVELOPMENT_COST * shortfall * shortfall + OVERDEVELOPMENT_COST * excess * excess
			+ FOLD_LABOR_COST * appliedFolds;
	}

	public RecipeEnvelope getRecipeEnvelope(){
		return recipeEnvelope;
	}

	/**
	 * Formulation cost of the hydration, salt, malt and oil of a recipe against the envelope of this product.
	 * <p>
	 * Fermentation alone cannot judge these three: more water speeds up the yeast, less salt too, more malt feeds
	 * it, so without a counterweight the optimizer would push them to the search bounds. Each term is a dead-zone
	 * quadratic: free inside the envelope the style is made with, growing with the squared distance from it
	 * (normalized by a tolerance) outside, plus a weak pull toward the center of the envelope that keeps the optimum
	 * unique.
	 * </p>
	 *
	 * @param hydration	Water ratio (g/g_flour).
	 * @param salt	Salt ratio (g/g_flour).
	 * @param malt	Diastatic malt ratio (g/g_flour).
	 * @param oil	Added oil ratio (g/g_flour).
	 * @return	Non-negative cost, to be added to the fitness (lower is better).
	 */
	public double computeFormulationCost(final double hydration, final double salt, final double malt,
			final double oil){
		return HYDRATION_COST * envelopeCost(hydration, recipeEnvelope.hydrationMin(), recipeEnvelope.hydrationMax(),
			HYDRATION_TOLERANCE, HYDRATION_CENTER_PULL)
			+ SALT_COST * envelopeCost(salt, recipeEnvelope.saltMin(), recipeEnvelope.saltMax(), SALT_TOLERANCE, CENTER_PULL)
			+ MALT_COST * envelopeCost(malt, recipeEnvelope.maltMin(), recipeEnvelope.maltMax(), MALT_TOLERANCE, CENTER_PULL)
			+ OIL_COST * envelopeCost(oil, recipeEnvelope.oilMin(), recipeEnvelope.oilMax(), OIL_TOLERANCE, CENTER_PULL);
	}

	private static double envelopeCost(final double value, final double min, final double max, final double tolerance,
			final double centerPull){
		final double outside = (value < min? min - value: (value > max? value - max: 0.));
		final double center = (min + max) / 2.;
		final double halfWidth = Math.max((max - min) / 2., tolerance);
		final double fromCenter = (value - center) / halfWidth;
		final double outsideScaled = outside / tolerance;
		return outsideScaled * outsideScaled + centerPull * fromCenter * fromCenter;
	}

	/**
	 * Accepted range of hydration, salt, malt and oil (g/g_flour) for a product, from common baker's formulations.
	 */
	public record RecipeEnvelope(double hydrationMin, double hydrationMax, double saltMin, double saltMax,
		double maltMin, double maltMax, double oilMin, double oilMax){}



	/**
	 * Computes the custom cost fitness value based on the final dough parameters.
	 * @param finalVolume	Simulated gas volume at the timeline end (ml/g_flour).
	 * @param maxPotential	Gas the gluten network can hold before tearing (ml/g_flour): drives the overproof side.
	 * @param referencePotential	Strength-independent reference capacity (ml/g_flour) the proofing target (underproof
	 * threshold) is derived from.
	 * @param residualSugar	Leftover fermentable sugar percentage at execution end (g/g_flour).
	 * @param finalYeast	Instantaneous active yeast cell concentration at the execution end (g/g_flour).
	 * @param finalYeastDot	Instantaneous final gas derivative (speed of volume growth).
	 * @return	Fitness score. The closer to negative infinity, the more optimal the recipe is.
	 */
	public abstract double computeFitness(double finalVolume, double maxPotential, double referencePotential,
		double residualSugar,
		double finalYeast, double finalYeastDot);

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
