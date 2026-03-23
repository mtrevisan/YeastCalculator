package io.github.mtrevisan.yeastcalculator.kinetics;


/**
 * Cardinal temperature model — Rosso et al. (1995).
 * <p>
 * Returns γ_T ∈ [0, 1]. γ_T = 1 at T_OPT.
 */
public class CardinalTemperatureModel{

	public static double factor(final double temperature, final double Tmin, final double Topt,
			final double Tmax){
		if(temperature <= Tmin || temperature >= Tmax)
			return 0.;

		final double numerator = (temperature - Tmax) * Math.pow(temperature - Tmin, 2.);
		final double denominator = (Topt - Tmin) * ((Topt - Tmin) * (temperature - Topt)
			- (Topt - Tmax) * (Topt + Tmin - 2. * temperature));
		return (denominator == 0.? 0.: numerator / denominator);
	}

}