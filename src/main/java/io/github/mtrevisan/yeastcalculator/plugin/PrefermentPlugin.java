package io.github.mtrevisan.yeastcalculator.plugin;

import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;
import io.github.mtrevisan.yeastcalculator.thermodynamic.GabMoistureModel;


/**
 * Plugin blueprint for indirect fermentation strategies (Biga, Poolish, Sourdough).
 * Acts as an isolated upstream interceptor that pre-conditions the state vector
 * and rheological limits before entering the core direct ODE system.
 */
public interface PrefermentPlugin{

	/**
	 * Unique identifier for the plugin registry mapping.
	 */
	String getPluginId();

	/**
	 * Intercepts and transforms the inputs, returning the modified recipe
	 * and a hook to pre-condition the initial state vector variables.
	 */
	PrefermentDecorator interceptAndPreCondition(SimulationInputs rawInputs);


	interface StatePreConditioner{
		/**
		 * Modifies the initial state array vector y = [X, S, V, EtOH]
		 * immediately prior to Phase 2 integration.
		 */
		void applyToInitialState(double[] initialY, double totalFlourSugarBudget);
	}

	record PrefermentDecorator(SimulationInputs decoratedInputs, StatePreConditioner statePreConditioner){}

}
