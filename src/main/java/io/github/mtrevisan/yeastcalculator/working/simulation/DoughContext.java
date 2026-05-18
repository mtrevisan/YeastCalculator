package io.github.mtrevisan.yeastcalculator.working.simulation;

import io.github.mtrevisan.yeastcalculator.working.domain.StageInput;

import java.util.Objects;


public final class DoughContext{

	private final StageInput[] stages;
	private final double[] folds;
	private final double totalDuration;
	private final double stiffness;
	private final double saltK;
	private final double oilK;
	private final double waterContent;
	private final double sugarInitial;


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
