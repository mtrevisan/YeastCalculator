package io.github.mtrevisan.yeastcalculator.simulation;

import io.github.mtrevisan.yeastcalculator.domain.StageInput;

import java.util.Objects;


/**
 * Immutable Parameter Object container serving as a thread-safe transport bag
 * for blended biomechanical and biochemical dough properties.
 * <p>
 * This class isolates the raw analytical components calculated during the environmental setup phase,
 * grouping them into a read-only matrix context that safely feeds the underlying ODE solver
 * without risking external cross-thread mutations or parallel state contamination.
 * </p>
 */
public final class DoughContext{

	private final StageInput[] stages;
	private final double[] folds;
	private final double totalDuration;
	private final double stiffness;
	private final double saltK;
	private final double oilK;
	private final double waterContent;
	private final double sugarInitial;


	/**
	 * Constructs a new immutable dough simulation context.
	 * Automatically applies defensive deep cloning to any mutable internal array inputs.
	 *
	 * @param stages	Array of active environmental timelines.
	 * @param folds	Array of timeline coordinates (hours) where degassing events occur.
	 * @param totalDuration	The aggregated processing timeline duration in hours.
	 * @param stiffness	The unified initial viscoelastic matrix stiffness index.
	 * @param saltK	The static osmotic growth inhibition multiplier from salt.
	 * @param oilK	The structural lipid barrier inhibition multiplier from oil.
	 * @param waterContent	The global mass fraction of total active water.
	 * @param sugarInitial	The starting aggregate mass fraction of fermentable carbohydrates.
	 */
	public static DoughContext create(final StageInput[] stages, final double[] folds, final double totalDuration,
			final double stiffness, final double saltK, final double oilK, final double waterContent,
			final double sugarInitial) {
		return new DoughContext(stages, folds, totalDuration, stiffness, saltK, oilK, waterContent, sugarInitial);
	}


	DoughContext(final StageInput[] stages, final double[] folds, final double totalDuration,
			final double stiffness, final double saltK, final double oilK, final double waterContent,
			final double sugarInitial){
		this.stages = Objects.requireNonNull(stages);
		this.folds = Objects.requireNonNull(folds);
		this.totalDuration = totalDuration;
		this.stiffness = stiffness;
		this.saltK = saltK;
		this.oilK = oilK;
		this.waterContent = waterContent;
		this.sugarInitial = sugarInitial;
	}


	public StageInput[] getStages(){
		return stages;
	}

	public double[] getFolds(){
		return folds;
	}

	public double getTotalDuration(){
		return totalDuration;
	}

	public double getStiffness(){
		return stiffness;
	}

	public double getSaltK(){
		return saltK;
	}

	public double getOilK(){
		return oilK;
	}

	public double getWaterContent(){
		return waterContent;
	}

	public double getSugarInitial(){
		return sugarInitial;
	}

}
