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
package io.github.mtrevisan.yeastcalculator.plugin;

import io.github.mtrevisan.yeastcalculator.domain.SimulationInputs;


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
		 * Modifies the initial state array vector y = [X, S, V, EtOH, D, h, d]; D (damaged starch), h (physiological state) and d (doublings, derived from the biomass change) are set by the optimizer
		 * immediately prior to Phase 2 integration.
		 */
		void applyToInitialState(double[] initialY, double totalFlourSugarBudget);
	}

	record PrefermentDecorator(SimulationInputs decoratedInputs, StatePreConditioner statePreConditioner){}

}
