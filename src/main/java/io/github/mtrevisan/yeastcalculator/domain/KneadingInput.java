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


/**
 * Describes mechanical parameters defining the matrix development profile.
 * <p>
 * The autolyze is the rest of the flour and the water alone (before adding the yeast and the salt): the flour hydrates
 * and the gluten starts to develop with no work, so fewer folds are needed afterward. It is not part of the
 * fermentation timeline, which starts when the yeast is added.
 * </p>
 *
 * @param autolyseMinutes	Autolyze duration [min], 0 for none.
 */
public record KneadingInput(
	KneadingType type,
	double durationMinutes,
	double autolyseMinutes){

	public KneadingInput{
		if(autolyseMinutes < 0.)
			throw new IllegalArgumentException("Autolyse duration cannot be negative");
	}

	/**
	 * Kneading with no autolyze.
	 */
	public KneadingInput(final KneadingType type, final double durationMinutes){
		this(type, durationMinutes, 0.);
	}



	public enum KneadingType{
		MANUAL(0.1, 0.6),
		SPIRAL_MIXER(1., 1.),
		FORK_MIXER(0.6, 0.75),
		DIVING_ARMS(0.8, 1.1);

		private final double frictionModifier;
		private final double developmentEfficiency;

		KneadingType(final double frictionModifier, final double developmentEfficiency){
			this.frictionModifier = frictionModifier;
			this.developmentEfficiency = developmentEfficiency;
		}

		public double getFrictionModifier(){
			return frictionModifier;
		}

		public double getDevelopmentEfficiency(){
			return developmentEfficiency;
		}
	}

}
