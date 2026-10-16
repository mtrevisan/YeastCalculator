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
	NEAPOLITAN_PIZZA(1.7, 0.009, 0.001, 25., 30., 0.58, 0.65, 0.){
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
	ROMAN_PAN_PIZZA(2.2, 0.013, 0.002, 15., 45., 0.75, 0.85, 0.08){
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
	GASTRONOMY_PAN_PIZZA(2., 0.015, 0.003, 20., 50., 0.65, 0.75, 0.04){
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
	BREAD(1.9, 0.012, 0.0015, 20., 40., 0.60, 0.72, 0.){
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


	// Cost of a hydration at the edge of the product range (quadratic pull toward the center of the range)
	private static final double HYDRATION_PREFERENCE_WEIGHT = 1.;
	// Scale of the penalty for a gluten network less developed than the hydration requires
	private static final double DEVELOPMENT_PENALTY_SCALE = 5.;
	// Required development index = BASE + SLOPE * (hydration - REFERENCE), clamped to [MIN, MAX]
	private static final double REQUIRED_DEVELOPMENT_BASE = 0.45;
	private static final double REQUIRED_DEVELOPMENT_SLOPE = 1.5;
	private static final double REQUIRED_DEVELOPMENT_REFERENCE_HYDRATION = 0.55;
	private static final double REQUIRED_DEVELOPMENT_MIN = 0.3;
	private static final double REQUIRED_DEVELOPMENT_MAX = 0.95;


	private final double glutenTearingLimit;
	private final double minSafeSugarThreshold;
	private final double maxAllowedFinalYeast;
	private final double tearingPenaltyMultiplier;
	private final double starvationPenaltyMultiplier;
	// Practical hydration range of the style [g_water/g_flour]
	private final double minHydration;
	private final double maxHydration;
	// Extra water a container-baked style can carry beyond the flour absorption capacity, thanks to the mechanical
	// support of the pan walls [g_water/g_flour]
	private final double panSupportHydrationBonus;


	BakeryProduct(final double glutenTearingLimit, final double minSafeSugarThreshold, final double maxAllowedFinalYeast,
			final double tearingPenaltyMultiplier, final double starvationPenaltyMultiplier, final double minHydration,
			final double maxHydration, final double panSupportHydrationBonus){
		this.glutenTearingLimit = glutenTearingLimit;
		this.minSafeSugarThreshold = minSafeSugarThreshold;
		this.maxAllowedFinalYeast = maxAllowedFinalYeast;
		this.tearingPenaltyMultiplier = tearingPenaltyMultiplier;
		this.starvationPenaltyMultiplier = starvationPenaltyMultiplier;
		this.minHydration = minHydration;
		this.maxHydration = maxHydration;
		this.panSupportHydrationBonus = panSupportHydrationBonus;
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

	public double getMinHydration(){
		return minHydration;
	}

	public double getMaxHydration(){
		return maxHydration;
	}

	public double getPanSupportHydrationBonus(){
		return panSupportHydrationBonus;
	}


	/**
	 * Minimum gluten development index a dough needs to hold its shape and gas at the given hydration.
	 * <p>
	 * Wetter doughs are slacker: the more free water dilutes the gluten network, the more mechanical development
	 * (kneading, stretch &amp; fold) is needed to obtain a cohesive, gas-tight structure.
	 * </p>
	 *
	 * @param hydration	Water ratio [g_water/g_flour].
	 * @return	Required development index in [0, 1].
	 */
	public static double requiredDevelopment(final double hydration){
		return Math.clamp(REQUIRED_DEVELOPMENT_BASE
				+ REQUIRED_DEVELOPMENT_SLOPE * (hydration - REQUIRED_DEVELOPMENT_REFERENCE_HYDRATION),
			REQUIRED_DEVELOPMENT_MIN, REQUIRED_DEVELOPMENT_MAX);
	}

	/**
	 * Computes the structural cost of the dough formulation, common to every product.
	 * <ul>
	 *    <li>Hydration preference: quadratic pull toward the center of the style hydration range (costing 1 at the
	 *    range edges), since crumb texture and handling are tied to the style hydration.</li>
	 *    <li>Gluten development: a network less developed than the hydration requires (slack, spreading dough) is
	 *    penalized proportionally to the product tearing sensitivity.</li>
	 * </ul>
	 *
	 * @param hydration	Water ratio [g_water/g_flour].
	 * @param developmentIndex	Gluten development index reached by kneading and folds, in [0, 1].
	 * @return	Structural penalty (zero at the center of the hydration range with a fully adequate development).
	 */
	public double computeStructuralPenalty(final double hydration, final double developmentIndex){
		final double targetHydration = (minHydration + maxHydration) / 2.;
		final double halfRange = (maxHydration - minHydration) / 2.;
		final double hydrationPenalty = HYDRATION_PREFERENCE_WEIGHT
			* Math.pow((hydration - targetHydration) / halfRange, 2.);

		final double developmentDeficit = Math.max(0., requiredDevelopment(hydration) - developmentIndex);
		final double developmentPenalty = tearingPenaltyMultiplier * DEVELOPMENT_PENALTY_SCALE * developmentDeficit;

		return hydrationPenalty + developmentPenalty;
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
