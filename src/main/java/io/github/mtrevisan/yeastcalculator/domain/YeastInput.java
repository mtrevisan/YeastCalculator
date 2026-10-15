package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Handles hydration and processing traits of the specific biological leavening agent.
 *
 * @param yeastMoisture            [g_water / g_wet_yeast]
 */
public record YeastInput(
	// [g_water / g_wet_yeast]
	double yeastMoisture
){

	public YeastInput{
		if(yeastMoisture < 0. || yeastMoisture > 1.)
			throw new IllegalArgumentException("Yeast moisture must be a fraction between 0 and 1");
	}


	/**
	 * Returns the recommended rehydration duration in minutes for this yeast type.
	 * <p>
	 * The value is a single scalar meant for display and for use as a default input.
	 * The biological tolerance around that value is expressed separately by
	 * {@link #getOptimalRehydrationWindow()}.
	 * </p>
	 * <p>
	 * Fresh yeast (moisture &ge; 0.65) is already hydrated and active, so it can be used
	 * immediately after crumbling. Dry yeast requires a short rehydration period to allow
	 * the cell wall to re-swell before mixing into the dough.
	 * </p>
	 */
	public double getRehydrationWindow(){
		// Fresh yeast block: cells are already active, immediate use is optimal.
		if(yeastMoisture >= 0.65)
			return 0.;

		// Dry mummified yeast: 12 minutes is the canonical rehydration time.
		return 12.;
	}

	/**
	 * Returns the biologically acceptable rehydration window as [min, max] in minutes.
	 * <p>
	 * The two boundaries represent distinct physiological failure modes:
	 * </p>
	 * <ul>
	 * <li>Below the lower bound the cell membrane has not completed its osmotic swelling.
	 *     Mixing at this stage exposes the partially rehydrated cell to shear damage and
	 *     leakage of intracellular solutes.</li>
	 * <li>Above the upper bound the cell begins to consume its internal trehalose reserves
	 *     in pure water without external substrates, triggering autolysis and early
	 *     starvation stress.</li>
	 * </ul>
	 * <p>
	 * Commercial yeasts are effectively bimodal: fresh blocks sit around 70% moisture and
	 * dry products (active dry, instant) sit between 4% and 8%. There is no meaningful
	 * commercial product in the 10-60% moisture range, so the window is modeled with two
	 * discrete regimes rather than a continuous interpolation.
	 * </p>
	 */
	public double[] getOptimalRehydrationWindow(){
		if(yeastMoisture >= 0.65)
			// Fresh yeast block: cells are already active, immediate use is optimal.
			// More than 15 minutes in pure water without substrates triggers starvation stress.
			return new double[]{0.0, 15.0};
		else
			// Dry mummified yeast: requires 9.6 to 30 minutes to safely re-swell membrane proteins.
			// Below 9.6 min: shear trauma. Above 30 min: starvation envelope.
			return new double[]{0.16 * 60.0, 0.50 * 60.0};
	}

}
