/* Copyright 2016. Martin Uecker.
 * All rights reserved. Use of this source code is governed by
 * a BSD-style license which can be found in the LICENSE file.
 * 
 * Authors:
 *	Martin Uecker <martin.uecker@med.uni-goettingen.de> 
 *	Nick Scholand
 */

#include <math.h>

#include "misc/misc.h"

#include "num/chebfun.h"

#include "specfun.h"


/* FIXME: improve precision
 * (but should be good enough for our purposes...)
 */

static const float coeff_0to8[] = {

	0.143432, 0.144372, 0.147260, 0.152300, 0.159883,
	0.170661, 0.185731, 0.207002, 0.238081, 0.286336,
	0.366540, 0.501252, 0.699580, 0.906853, 1.000000,
};

static const float coeff_8toinf[] = {

	0.405687, 0.405664, 0.405601, 0.405494, 0.405349,
	0.405164, 0.404945, 0.404692, 0.404413, 0.404107,
	0.403782, 0.403439, 0.403086, 0.402724, 0.402359,
	0.401995, 0.401637, 0.401287, 0.400951, 0.400631,
	0.400332, 0.400055, 0.399805, 0.399582, 0.399391,
	0.399231, 0.399106, 0.399012, 0.398998, 0.399001
};



/*
 * modified bessel function
 */
double bessel_i0(double x)
{
	if (x < 0.)
		return bessel_i0(-x);

	if (x < 8.)
		return exp(x) * chebeval(x  / 4. - 1., ARRAY_SIZE(coeff_0to8), coeff_0to8);

	return exp(x) * chebeval(16. / x - 1., ARRAY_SIZE(coeff_8toinf), coeff_8toinf) / sqrt(x);
}



static long factorial(long k)
{
	return (0 == k) ? 1 : (k * factorial(k - 1));
}



// approximate sine integral with power series (only for small x)
double Si_power(double x)
{
	int k_max = 10;
	double sum = 0;

	for (int k = 1; k < k_max; k++)
		sum += pow(-1. , (k - 1)) * pow(x, (2 * k - 1)) / ((2 * k - 1) * factorial(2 * k - 1));

	return sum;
}


// Efficient and accurate calculation of Sine Integral using Padé approximants of the convergent Taylor series
// For details see:
// 	Rowe, B., et al.
//	"GALSIM: The modular galaxy image simulation toolkit".
//	Astronomy and Computing. 10: 121. 2015.
//	arXiv:1407.7676
//  		-> Appendix B: Efficient evaluation of the Sine and Cosine integrals

// helper function to calculate Si accurate for large arguments (> 4)
static double Si_help_f(double x)
{
	double num =	1
			+ 7.44437068161936700618 * 10E2		* pow(x,-2)
			+ 1.96396372895146869801 * 10E5		* pow(x,-4)
			+ 2.37750310125431834034 * 10E7		* pow(x,-6)
			+ 1.43073403821274636888 * 10E9		* pow(x,-8)
			+ 4.33736238870432522765 * 10E10	* pow(x,-10)
			+ 6.40533830574022022911 * 10E11	* pow(x,-12)
			+ 4.20968180571076940208 * 10E12	* pow(x,-14)
			+ 1.00795182980368574617 * 10E13	* pow(x,-16)
			+ 4.94816688199951963482 * 10E12	* pow(x,-18)
			- 4.94701168645415959931 * 10E11	* pow(x,-20);

	double denum =	1
			+ 7.46437068161927678031 * 10E2		* pow(x,-2)
			+ 1.97865247031583951450 * 10E5		* pow(x,-4)
			+ 2.41535670165126845144 * 10E7		* pow(x,-6)
			+ 1.47478952192985464958 * 10E9		* pow(x,-8)
			+ 4.58595115847765779830 * 10E10	* pow(x,-10)
			+ 7.08501308149515401563 * 10E11	* pow(x,-12)
			+ 5.06084464593475076774 * 10E12	* pow(x,-14)
			+ 1.43468549171581016479 * 10E13	* pow(x,-16)
			+ 1.11535493509914254097 * 10E13	* pow(x,-18);

	return num / denum / x;
}

// helper function to calculate Si accurate for large arguments (> 4)
static double Si_help_g(double x)
{
	double num =	1
			+ 8.1359520115168615 * 10E2	* pow(x,-2)
			+ 2.35239181626478200 * 10E5	* pow(x,-4)
			+ 3.12557570795778731 * 10E7	* pow(x,-6)
			+ 2.06297595146763354 * 10E9	* pow(x,-8)
			+ 6.83052205423625007 * 10E10	* pow(x,-10)
			+ 1.09049528450362786 * 10E12	* pow(x,-12)
			+ 7.57664583257834349 * 10E12	* pow(x,-14)
			+ 1.81004487464664575 * 10E13	* pow(x,-16)
			+ 6.43291613143049485 * 10E12	* pow(x,-18)
			- 1.36517137670871689 * 10E12	* pow(x,-20);

	double denum =	1
			+ 8.19595201151451564 * 10E2	* pow(x,-2)
			+ 2.40036752835578777 * 10E5	* pow(x,-4)
			+ 3.26026661647090822 * 10E7	* pow(x,-6)
			+ 2.23355543278099360 * 10E9	* pow(x,-8)
			+ 7.87465017341829930 * 10E10	* pow(x,-10)
			+ 1.39866710696414565 * 10E12	* pow(x,-12)
			+ 1.17164723371736605 * 10E13	* pow(x,-14)
			+ 4.01839087307656620 * 10E13	* pow(x,-16)
			+ 3.99653257887490811 * 10E13	* pow(x,-18);

	return ((x < 0) ? -1 : 1) * num / denum / x / x;
}

static double Si_large_x(double x)
{
	return M_PI / 2 - Si_help_f(x) * cos(x) - Si_help_g(x) * sin(x);
}

double Si(double x)
{
	// Definition of Si_large_x just for x > 0,
	// therefore use Si(-z) = -Si(z)
	// For more information compare:
	// 	Abramowitz, M., & Stegun, I. 1964, Handbook of Mathematical Functions, 5th edn. (New York: Dover)
	// 		-> Chapter 5.2

	return (fabs(x) <= 4) ? Si_power(x) : ((x < 0) ? -1 : 1) * Si_large_x(fabs(x));
}
