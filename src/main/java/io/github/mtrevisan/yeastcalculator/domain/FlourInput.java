package io.github.mtrevisan.yeastcalculator.domain;


/**
 * Encapsulates chemical and rheological specifications of a single flour instance.
 */
public record FlourInput(
	double strength,
	double plRatio,
	double sugar,
	double protein,
	double fat,
	double fiber,
	double ash,
	FlourType type){}
