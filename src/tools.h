/*
 * tools.h
 *
 *  Created on: 6 oct. 2022
 *      Author: Kévin Carrier
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpfr.h" // in "/opt/homebrew/Cellar/mpfr/4.1.0/include/"
#include "gmp.h" // in "/opt/homebrew/Cellar/gmp/6.2.1_1/include/"
#include "flint.h" // in "/opt/homebrew/Cellar/flint/2.9.0/include/flint/"
#include "fmpz_poly.h" // in "/opt/homebrew/Cellar/flint/2.9.0/include/flint/"

#define PRECISION 0.000000001

// Verbosity flags
#define VERB 1
#define NOVERB 0

// Asymptotic flags
#define ASYMPT 1
#define NOASYMPT 0


// Model for computing the optimal sieving dimension for Solving SVP for lattices.
#define ASYMPTOTIC_MODEL 1
#define G6K_MODEL 2


#define C_PROG 5.45759151825911 // 1.0 / (1.0 - pow(2.0,-0.292))  Matzov p.37
#define C_MUL 1024.0 // Matzov p.37 : size_word * size_word
#define C_ADD 160.0 //guessing based on C_mul : 5.0 * size_word


/*
 *  @brief The time complexity of the FFT
 *  @param q Modulus >= 2
 *  @param n Dimension
 *  @return The time complexity of the FFT in dimension n with modulus p.
 */
double Tfft_f(const double q, const double d){
	double log2_q = log2(q);
	if (fabs(round(log2_q) - log2_q) > PRECISION){
		return C_MUL * d * pow(q , d + 1);  // Matzov Theorem 7.6, p.38
	} else {
		log2_q = round(log2_q);
		return C_MUL * log2_q * d * pow(q , d);  // [Knu97, PP.306-311]
	}
}


/*
 *  @brief  Time complexity of updating the table in each iteration.
 *  @param N Number of nonzero entries
 *  @return  Time complexity of updating the table of size N in each iteration.
 */
double Ttable_f(const double N){
	return 4.0 * C_ADD * N; // Matzov Theorem 7.6, p.39
}

/*
 *  @brief  Time complexity for decoding all the samples.
 *  @param N Number of nonzero entries
 *  @return  Time complexity for decoding all the samples.
 */
double Tdecode_f(const double N, const double p, const double n){
	return 2.0 * C_MUL * N * p * log2(p) * n * log2(n);
}

/*
 * @brief Code de Martin Albretch. Compute root-Hermite factor.
 * @param beta Block size
 * @return root-Hermite factor delta from block size beta.
 */
double delta_f(const double beta){
	//double beta_ = beta; // in Matzov
	double beta_ = round(beta);

	double small[8][2] = {
			{2, 1.02190},  // noqa
			{5, 1.01862},  // noqa
			{10, 1.01616},
			{15, 1.01485},
			{20, 1.01420},
			{25, 1.01342},
			{28, 1.01331},
			{40, 1.01295},
	};

	if (beta_ <= 2){
		return 1.0219;
	} else if (beta_ < 40){
		for (int i = 1 ; i < 8 ; ++i){
			if (beta_ < small[i][0]) {
				return small[i - 1][1];
			}
		}
	} else if (beta_ == 40) {
		return 1.01295;
	}
	return pow( (beta_ / (2 * M_PI * M_E)) * pow(M_PI * beta_ , 1.0 / beta_) , 1.0 / (2.0 * (beta_ - 1.0)) );
}

double coth(const double x){
	return cosh(x)/sinh(x);
}

/*
 * @brief Code de Martin Albretch. Entropy.
 * @param sigma standard deviation of the random variable
 * @return Entropy of the Gaussian of standard deviation sigma.
 */
double H_f(const double sigma){
	return (
			0.5
			+ log(sqrt(2 * M_PI) * sigma)
			+ log(coth(M_PI * M_PI * sigma * sigma))
	) / log(2.0);
}

/*
 *  @brief Optimal sieving dimension for Solving SVP for lattices.
 *  @param beta dimension
 *  @param dimension_for_free_model The Dimension-For-Free model : ASYMPTOTIC_MODEL from [Duc18] or G6K_MODEL from [ADH+19]
 *  @return The optimal sieving dimension for Solving SVP for lattices in dimension beta.
 */
double beta_eff_f(const double beta, const int dimension_for_free_model){
	if (beta < 25){
		fprintf(stderr, "beta_eff_f error: beta must be greater than 25 for the ASYMPTOTIC_MODEL");
		exit(EXIT_FAILURE);
	}
	if (dimension_for_free_model == ASYMPTOTIC_MODEL){
		//return beta - (beta * log(4.0 / 3.0) / log(beta / (2 * M_PI * M_E))); // in Matzov
		return beta - floor( fmax(beta * log(4.0 / 3.0) / log(beta / (2 * M_PI * M_E)), 0.0) );
	} else if (dimension_for_free_model == G6K_MODEL) {
		//return beta - (11.46 + (0.0757*beta)); // in Matzov
		return beta - floor(11.46 + (0.0757*beta));
	} else {
		fprintf(stderr, "beta_eff_f error: the given dimension_for_free_model is not available");
		exit(EXIT_FAILURE);
	}
}

double CDF(const double x){
	return 0.5 + 0.5*erf(x/sqrt(2));
}

/*
 * @brief from https://gist.github.com/kmpm/1211922/6b7fcd0155b23c3dc71e6f4969f2c48785371292
 */
double invCDF(const double y)
{
    if (y < -PRECISION || y > 1.0 + PRECISION)
    {
    	fprintf(stderr, "invCDF error: p = %.12f must be between 0 and 1", y);
    	exit(EXIT_FAILURE);
    }
    double p = y;
    if (y < PRECISION){
    	p = PRECISION;
    }
    if (y > 1.0 - PRECISION){
    	p = 1.0 - PRECISION;
    }

    double r, val;

    const double q = p - 0.5;

    if (fabs(q) <= .425) {
        r = .180625 - q * q;
        val =
            q * (((((((r * 2509.0809287301226727 +
                33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r +
                1971.5909503065514427) * r + 133.14166789178437745) * r +
                3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r +
                687.1870074920579083) * r + 42.313330701600911252) * r + 1);
    }
    else {
        if (q > 0) {
            r = 1 - p;
        }
        else {
            r = p;
        }

        r = sqrt(-log(r));

        if (r <= 5)
        {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                .0227238449892691845833) * r + .24178072517745061177) *
                r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) *
                r + 4.6303378461565452959) * r +
                1.42343711074968357734)
                / (((((((r *
                    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                    .14810397642748007459) * r + .68976733498510000455) *
                    r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1);
        }
        else { /* very close to  0 or 1 */
            r += -5;
            val = (((((((r * 2.01033439929228813265e-7 +
                2.71155556874348757815e-5) * r +
                .0012426609473880784386) * r + .026532189526576123093) *
                r + .29656057182850489123) * r +
                1.7848265399172913358) * r + 5.4637849111641143699) *
                r + 6.6579046435011037772)
                / (((((((r *
                    2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    r + 1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r + .0148753612908506148525)
                    * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1);
        }

        if (q < 0.0) {
            val = -val;
        }
    }

    return val;
}






/*
 * @brief Example for fmpz_poly
 */
void example_fmpz_poly(){
	unsigned int t = 2, n = 5;
	fmpz_poly_t R, P, Q, MOD;
	fmpz_t N, UN;
	fmpz_init_set_ui(UN, 1);
	fmpz_init(N);
	fmpz_poly_init(R);
	fmpz_poly_init(Q);
	fmpz_poly_init(P);
	fmpz_poly_init(MOD);

	fmpz_poly_set_coeff_ui(P, 0, 1);
	fmpz_poly_set_coeff_ui(P, 1, 2);
	fmpz_poly_set_coeff_ui(P, 4, 2);


	fmpz_poly_pow(Q, P, n);

	fmpz_poly_set_coeff_ui(MOD, t+1, 1);
	fmpz_poly_rem(R, Q, MOD);

	fmpz_poly_evaluate_fmpz(N, R, UN);

	char *N_str = fmpz_get_str(NULL, 10, N);
	printf("%s\n", N_str);

	char *s = "X";
	fmpz_poly_print_pretty(P, s);

	fmpz_poly_clear(R);
	fmpz_poly_clear(Q);
	fmpz_poly_clear(P);
	fmpz_poly_clear(MOD);
	fmpz_clear(UN);
	fmpz_clear(N);
	free(N_str);
}


#endif /* TOOLS_H_ */
