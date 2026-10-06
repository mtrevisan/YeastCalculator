package io.github.mtrevisan.yeastcalculator.plugin;

import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;


public class PoolishPlugin implements PrefermentPlugin{

	// Fraction of total recipe flour used in Poolish (e.g., 0.30 for 30%)
	private final double flourFractionInPoolish;
	private final double durationHours;
	private final double temperature;


	public PoolishPlugin(final double flourFraction, final double durationHours, final double temperature){
		this.flourFractionInPoolish = flourFraction;
		this.durationHours = durationHours;
		this.temperature = temperature;
	}


	@Override
	public String getPluginId(){
		return "POOLISH_LIQUID_PREFERMENT_PLUGIN";
	}

	@Override
	public PrefermentDecorator interceptAndPreCondition(final SimulationInputs rawInputs){
		final DoughRecipe originalRecipe = rawInputs.recipe();

		// --- 1. RECIPE WATER RE-BALANCING MASS BALANCE ---
		// A traditional Poolish uses equal weights of flour and water (100% hydration relative to itself).
		// Water weight in Poolish = totalFlour * flourFractionInPoolish * 1.0.
		// Therefore, the water ratio subtracted from the main recipe pool is exactly equal to the flourFraction.
		final double waterUsedInPoolish = flourFractionInPoolish;
		double adjustedRefreshWaterRatio = originalRecipe.waterRatio() - waterUsedInPoolish;
		// Fallback protection to avoid unworkable negative hydration values if user configurations conflict
		if(adjustedRefreshWaterRatio < 0.35)
			adjustedRefreshWaterRatio = 0.35;

		final DoughRecipe decoratedRecipe = new DoughRecipe(
			adjustedRefreshWaterRatio,
			originalRecipe.saltRatio(),
			originalRecipe.oilRatio(),
			originalRecipe.maltRatio(),
			originalRecipe.maltSugarContent(),
			originalRecipe.maltPollakUnit(),
			// High water poolish lubricates initial mix, reducing friction heat
			originalRecipe.mixerFrictionFactor() - 0.5
		);

		final SimulationInputs decoratedInputs = new SimulationInputs(
			rawInputs.fractions(), rawInputs.flourMatrix(), rawInputs.flourTemperature(),
			rawInputs.airRelativeHumidity(), rawInputs.yeastProperties(), decoratedRecipe,
			rawInputs.kneading(), rawInputs.stages(), rawInputs.folds()
		);

		// --- 2. CLOSED STATE PRE-CONDITIONER INTERCEPTOR HOOK ---
		final StatePreConditioner preConditioner = (initialY, totalFlourSugarBudget) -> {
			// Read what YeastOptimizer set as the base baseline starting active biomass fraction (30% solids)
			final double baseActiveX = initialY[0];

			// Liquid poolish specific growth calculation kinetics (highly efficient mass transport)
			// Biomass multiplies exponentially over time under warm ambient temperatures
			final double growthRateConstant = 0.085 * (temperature / 28.0);
			final double theoreticalBiomassMultiplicationFactor = Math.exp(growthRateConstant * durationHours);

			// Theoretical sugar consumption required for this growth step
			final double theoreticalConsumedSugar = (baseActiveX * (theoreticalBiomassMultiplicationFactor - 1.)) / 0.12;

			// --- FIXED: ENFORCE THE LOCAL FLOUR SUGAR BUDGET SAFETY WALL ---
			// Calculate the absolute maximum simple sugar available inside the Poolish flour slice
			final double poolishFlourSugarBudget = totalFlourSugarBudget * flourFractionInPoolish;
			double actualConsumedSugar;
			double effectiveBiomassMultiplicationFactor;
			if(theoreticalConsumedSugar > poolishFlourSugarBudget){
				// Yeast starved! Cap consumption to the absolute maximum available in the Poolish partition
				actualConsumedSugar = poolishFlourSugarBudget;
				// Reverse-calculate the effective biomass cap supported by the empty budget
				effectiveBiomassMultiplicationFactor = 1. + ((actualConsumedSugar * 0.12) / baseActiveX);
			}
			else{
				// Growth completes uninhibited within safe substrate limits
				actualConsumedSugar = theoreticalConsumedSugar;
				effectiveBiomassMultiplicationFactor = theoreticalBiomassMultiplicationFactor;
			}

			// Cap the absolute maximum multiplication envelope to protect numerical stability bounds
			effectiveBiomassMultiplicationFactor = Math.min(3.5, effectiveBiomassMultiplicationFactor);

			// Re-calculate the initial active matter payload injected into the final dough matrix [Index 0 = Biomass]
			initialY[0] = baseActiveX * effectiveBiomassMultiplicationFactor;

			// Subtract the actually consumed sugars from the initial global pool [Index 1 = Sugars]
			initialY[1] = Math.max(0.002, initialY[1] - actualConsumedSugar);

			// Pre-generated ethanol injection accumulation (Gay-Lussac conversion yield ~48%) [Index 3 = Ethanol]
			final double preGeneratedEthanol = actualConsumedSugar * 0.48;
			// Diluted across total final dough mass
			initialY[3] = preGeneratedEthanol * flourFractionInPoolish;
		};

		return new PrefermentDecorator(decoratedInputs, preConditioner);
	}

}
