package io.github.mtrevisan.yeastcalculator.plugin;

import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.domain.DoughRecipe;


public class BigaPlugin implements PrefermentPlugin{

	private final double flourFractionInBiga; // Fraction of total recipe flour used in Biga (e.g., 0.50 for 50%)
	private final double durationHours;
	private final double temperature;


	public BigaPlugin(final double flourFraction, final double durationHours, final double temperature){
		this.flourFractionInBiga = flourFraction;
		this.durationHours = durationHours;
		this.temperature = temperature;
	}


	@Override
	public String getPluginId(){
		return "BIGA_SOLID_PREFERMENT_PLUGIN";
	}

	@Override
	public PrefermentDecorator interceptAndPreCondition(final SimulationInputs rawInputs){
		final DoughRecipe originalRecipe = rawInputs.recipe();

		// --- 1. RECIPE WATER RE-BALANCING MASS BALANCE ---
		// A standard Italian Biga is bound at a rigid 44% hydration relative to its specific flour weight.
		// Water weight in Biga = totalFlour * flourFractionInBiga * 0.44.
		final double waterUsedInBiga = flourFractionInBiga * 0.44;
		double adjustedRefreshWaterRatio = originalRecipe.waterRatio() - waterUsedInBiga;
		if(adjustedRefreshWaterRatio < 0.35)
			adjustedRefreshWaterRatio = 0.35;

		final DoughRecipe decoratedRecipe = new DoughRecipe(
			adjustedRefreshWaterRatio,
			originalRecipe.saltRatio(),
			originalRecipe.oilRatio(),
			originalRecipe.maltRatio(),
			originalRecipe.maltSugarContent(),
			originalRecipe.maltPollakUnit(),
			// Crumbly solid Biga pieces increase mixer friction resistance heat!
			originalRecipe.mixerFrictionFactor() + 1.5
		);

		final SimulationInputs decoratedInputs = new SimulationInputs(
			rawInputs.fractions(), rawInputs.flourMatrix(), rawInputs.flourTemperature(),
			rawInputs.airRelativeHumidity(), rawInputs.yeastProperties(), decoratedRecipe,
			rawInputs.kneading(), rawInputs.stages(), rawInputs.folds()
		);

		// --- 2. CLOSED STATE PRE-CONDITIONER INTERCEPTOR HOOK ---
		final StatePreConditioner preConditioner = (initialY, totalFlourSugarBudget) -> {
			final double baseActiveX = initialY[0];

			final double growthRateConstant = 0.042 * (temperature / 20.0);
			final double theoreticalBiomassMultiplicationFactor = Math.exp(growthRateConstant * durationHours);

			final double theoreticalConsumedSugar = (baseActiveX * (theoreticalBiomassMultiplicationFactor - 1.)) / 0.10;

			// --- FIXED: ENFORCE THE LOCAL BIGA SUGAR BUDGET SAFETY WALL ---
			final double bigaFlourSugarBudget = totalFlourSugarBudget * flourFractionInBiga;
			double actualConsumedSugar;
			double effectiveBiomassMultiplicationFactor;
			if(theoreticalConsumedSugar > bigaFlourSugarBudget){
				actualConsumedSugar = bigaFlourSugarBudget;
				effectiveBiomassMultiplicationFactor = 1. + (actualConsumedSugar * 0.1) / baseActiveX;
			}
			else{
				actualConsumedSugar = theoreticalConsumedSugar;
				effectiveBiomassMultiplicationFactor = theoreticalBiomassMultiplicationFactor;
			}

			effectiveBiomassMultiplicationFactor = Math.min(2.2, effectiveBiomassMultiplicationFactor);

			// Apply updates to the state vector coordinates
			// Biomass X
			initialY[0] = baseActiveX * effectiveBiomassMultiplicationFactor;
			// Sugars S
			initialY[1] = Math.max(0.002, initialY[1] - actualConsumedSugar);

			final double preGeneratedEthanol = actualConsumedSugar * 0.46;
			// Ethanol EtOH
			initialY[3] = preGeneratedEthanol * flourFractionInBiga;
		};

		return new PrefermentDecorator(decoratedInputs, preConditioner);
	}

}
