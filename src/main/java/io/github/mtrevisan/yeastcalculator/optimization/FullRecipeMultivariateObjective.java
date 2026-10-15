package io.github.mtrevisan.yeastcalculator.optimization;

import io.github.mtrevisan.yeastcalculator.domain.BakeryProduct;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;
import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.StageInput;
import io.github.mtrevisan.yeastcalculator.simulation.DoughOdeSystem;
import io.github.mtrevisan.yeastcalculator.simulation.FoldEventHandler;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;
import org.apache.commons.math3.analysis.MultivariateFunction;


/**
 * Multivariate objective function for a FIXED number of stretch-and-fold sets.
 * <p>
 * The fold count is a discrete variable and cannot be optimized by BOBYQA, which assumes a
 * smooth objective. It is therefore injected at construction time, and the search vector
 * contains only the six continuous variables: [Hydration, Salt, Oil, Malt, InitialDelay,
 * RelaxationFactor]. The caller (GlobalMultivariateRecipeOptimizer) iterates over the fold
 * count externally and compares the resulting fitness values.
 * </p>
 */
public class FullRecipeMultivariateObjective implements MultivariateFunction{

	private final SimulationInputs baseInputs;
	private final BakeryProduct targetProduct;
	private final double[] lowerBounds;
	private final double[] upperBounds;
	private final int numberOfFolds;


	public FullRecipeMultivariateObjective(final SimulationInputs baseInputs, final BakeryProduct targetProduct,
			final double[] lowerBounds, final double[] upperBounds, final int numberOfFolds){
		this.baseInputs = baseInputs;
		this.targetProduct = targetProduct;
		this.lowerBounds = lowerBounds;
		this.upperBounds = upperBounds;
		this.numberOfFolds = numberOfFolds;
	}


	/**
	 * Maps a physical point onto the normalized [0, 1] search space.
	 */
	public double[] normalize(final double[] point){
		final int length = point.length;
		final double[] normalizedPoint = new double[length];
		for(int i = 0; i < length; i ++)
			normalizedPoint[i] = (point[i] - lowerBounds[i]) / (upperBounds[i] - lowerBounds[i]);
		return normalizedPoint;
	}

	/**
	 * Maps a normalized [0, 1] search vector onto its physical bounds.
	 */
	public double[] denormalize(final double[] normalizedPoint){
		final int length = normalizedPoint.length;
		final double[] point = new double[length];
		for(int i = 0; i < length; i ++)
			point[i] = denormalize(normalizedPoint[i], i);
		return point;
	}

	private double denormalize(final double normalizedValue, final int index){
		return lowerBounds[index] + normalizedValue * (upperBounds[index] - lowerBounds[index]);
	}

	@Override
	public double value(final double[] point){
		// Denormalize inputs from [0, 1] to physical ranges
		final double[] p = denormalize(point);
		final double hydration = p[0];
		final double salt = p[1];
		final double oil = p[2];
		final double malt = p[3];
		final double initialFoldDelay = p[4];
		final double relaxationFactor = p[5];

		// The fold count is fixed externally
		final double[] dynamicFolds = buildFoldSchedule(numberOfFolds, initialFoldDelay, relaxationFactor);

		// Total timeline in minutes
		double totalMinutes = 0.;
		for(final StageInput stage : baseInputs.stages())
			totalMinutes += stage.duration() * 60.;

		// Construct the candidate recipe payload
		final DoughRecipe candidateRecipe = new DoughRecipe(hydration, salt, oil, malt,
			baseInputs.recipe().maltSugarContent(), baseInputs.recipe().maltPollakUnit(),
			baseInputs.recipe().mixerFrictionFactor());

		final SimulationInputs candidateInputs = new SimulationInputs(baseInputs.fractions(),
			baseInputs.flourMatrix(), baseInputs.flourTemperature(), baseInputs.airRelativeHumidity(),
			baseInputs.yeastProperties(), candidateRecipe, baseInputs.kneading(), baseInputs.stages(),
			dynamicFolds);

		// 1. Solve the nested univariate problem for the optimal yeast target fraction
		final double optimalYeast = YeastOptimizer.findOptimalYeast(candidateInputs, targetProduct, null);

		// 2. Execute the full system ODE integration loop
		final GabMoistureModel.GabResult gab = GabMoistureModel.calculateMoisture(candidateInputs);
		final double[] finalState = YeastOptimizer.runSimulation(candidateInputs, gab, optimalYeast, targetProduct);

		if(finalState == null)
			// Return maximum penalty if numerical system crashes
			return Double.MAX_VALUE;

		// 3. Reconstruct the dynamic Bloksma rheological constraints at the timeline endpoint
		final double[] summary = YeastOptimizer.inputsSummary(candidateInputs);
		final double baseMaxGasPotential = YeastOptimizer.calculateMaxGasPotential(summary, candidateInputs,
			targetProduct);
		final StageInput lastStage = baseInputs.stages()[baseInputs.stages().length - 1];
		final double lastStageRH = lastStage.relativeHumidity();
		final double skinningModifier = DoughOdeSystem.skinningModifier(lastStageRH);

		// Fold handler must reflect the dynamic folds; only those falling inside the timeline count
		final FoldEventHandler foldHandler = new FoldEventHandler(dynamicFolds);
		foldHandler.advanceFoldsUpTo(totalMinutes);
		final int appliedFolds = foldHandler.getAppliedFolds();

		final double dynamicMaxPotential = baseMaxGasPotential
			* Math.pow(FoldEventHandler.STRUCTURAL_GAIN, appliedFolds)
			* skinningModifier;

		// Compute the final derivative at the endpoint, using the fold-aware handler
		final double[] yDotFinal = new double[4];
		final DoughOdeSystem finalOde = new DoughOdeSystem(lastStage.temperature(), lastStageRH,
			gab, baseMaxGasPotential, candidateInputs, foldHandler);
		finalOde.computeDerivatives(totalMinutes, finalState, yDotFinal);

		// Evaluate and return global multi-variable fitness cost
		return targetProduct.computeFitness(finalState[2], dynamicMaxPotential, finalState[1], finalState[0],
			yDotFinal[2], candidateRecipe.oilRatio(), appliedFolds, FoldEventHandler.DEGASSING_RETENTION);
	}

	public static double[] buildFoldSchedule(final int optimizedFoldCount, double currentDelay,
			final double relaxationFactor){
		final double[] optimizedFoldsArray = new double[optimizedFoldCount];
		// Serves as initial fold delay field variable mapping
		double cumulativeMinutes = 0.;
		for(int i = 0; i < optimizedFoldCount; i ++){
			cumulativeMinutes += currentDelay;
			optimizedFoldsArray[i] = cumulativeMinutes;
			currentDelay *= relaxationFactor;
		}
		return optimizedFoldsArray;
	}

}
