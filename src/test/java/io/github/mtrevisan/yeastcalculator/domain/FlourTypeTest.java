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

import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;


class FlourTypeTest{

	@Test
	void labelsAreUniqueAndNotEmpty(){
		final Set<String> labels = new HashSet<>();
		for(final FlourType type : FlourType.values()){
			assertNotNull(type.getLabel());
			assertTrue(!type.getLabel().isBlank(), type + " has an empty label");
			assertTrue(labels.add(type.getLabel()), "duplicate label " + type.getLabel());
		}
	}

	@Test
	void gabConstantsArePhysicallyPlausible(){
		for(final FlourType type : FlourType.values()){
			// base monolayer (without the protein and fiber contributions): 0.5-10 g/100 g d.b.
			assertTrue(type.getWmBase() > 0.004 && type.getWmBase() < 0.10, type + " wmBase " + type.getWmBase());
			// type II isotherm (sigmoid) requires C > 2; flours are well below 30
			assertTrue(type.getcGab() > 2. && type.getcGab() < 30., type + " C " + type.getcGab());
			// K < 1 (multilayer water less bound than bulk water, and the isotherm must stay finite for a_w up to 0.95);
			// food GAB fits are typically in 0.7-1.0
			assertTrue(type.getkGab() > 0.7 && type.getkGab() < 1., type + " K " + type.getkGab());
			assertTrue(type.getkGab() * 0.95 < 1., type + " K·a_w must be below 1");
		}
	}

	@Test
	void gabConstantsKeepTheModelMonolayerInTheCerealToNutRange(){
		// monolayer = wmBase + 0.085·protein + 0.12·fiber, with extreme protein/fiber of 40 % and 40 %: below 15 g/100 g d.b.
		for(final FlourType type : FlourType.values())
			assertTrue(type.getWmBase() + 0.085 * 0.4 + 0.12 * 0.4 < 0.15, type + " monolayer too large");
	}

}
