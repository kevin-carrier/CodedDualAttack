/*
 ============================================================================
 Name        : CodedDualAttack.c
 Author      : Kévin Carrier
 Version     :
 Copyright   : open source
 Description : Coded Dual Attack in C, Ansi-style
 ============================================================================
 */

#include "tools.h"
#include "centered_binomial.h"
#include "nns.h"
#include "lwe.h"


typedef struct DualAttackParams{
	double m; // the number of used LWE samples
	double n_enum; // size of the enumerating part
	double n_fft; // length of the code (size of the zone for the modulo switching, decoding and FFT)
	double n_lat; // size of the lattice part: n_lat = n - n_enum - n_fft
	double n_code; // code dimension
	double p; // new modulus after modulo switching
	double beta_bkz; // BKZ block size
	double beta_sieve; // sieving dimension
	double mu; // failure probability of inner loop
	double nu; // failure probability of outer loop
} DualAttackParams_t;

void print_DualAttackParams(DualAttackParams_t *p){
	printf("Parameters of the algorithm:\n");
	printf("m = %.0f\n", p->m);
	printf("n_enum = %.0f\n", p->n_enum);
	printf("n_fft = %.0f\n", p->n_fft);
	printf("n_lat = %.0f\n", p->n_lat);
	printf("n_code = %.0f\n", p->n_code);
	printf("p = %.0f\n", p->p);
	printf("beta_bkz = %.0f\n", p->beta_bkz);
	printf("beta_sieve = %.0f\n", p->beta_sieve);
	printf("mu = %.40f\n", p->mu);
	printf("nu = %.40f\n\n", p->nu);
}

void copy_parameters(DualAttackParams_t *dest, DualAttackParams_t *src){
	dest->m = src->m;
	dest->n_enum = src->n_enum;
	dest->n_fft = src->n_fft;
	dest->n_lat = src->n_lat;
	dest->n_code = src->n_code;
	dest->p = src->p;
	dest->beta_bkz = src->beta_bkz;
	dest->beta_sieve = src->beta_sieve;
	dest->mu = src->mu;
}

/*
 *  @brief Required number of samples to distinguish with advantage.
 *  @param LWE The LWE parameters
 *  @param ALGO The parameters of the algorithm
 *  @param Nenum Number of enumerated vector s_enum
 *  @param flag_asympt Bollean for asymptotic computation or not
 *  @param flag_verb Boolean for verbosity mode
 *  @return Required number of samples to distinguish with advantage.
 */
double N_f(LWEparams_t *LWE, DualAttackParams_t *ALGO, double Nenum, int flag_asympt, int flag_verb){
	double lsigma_s = (
			pow(LWE->sigma_e , ALGO->m / (ALGO->m + ALGO->n_lat))
			* pow(LWE->sigma_s * LWE->q , ALGO->n_lat / (ALGO->m + ALGO->n_lat))
			* sqrt(4.0 / 3.0)
			* sqrt(ALGO->beta_sieve / (2.0 * M_PI * M_E))
			* pow(delta_f(ALGO->beta_bkz) , ALGO->m + ALGO->n_lat - ALGO->beta_sieve)
	);
	double Deq = exp(4 * pow(lsigma_s * M_PI / LWE->q , 2));

	double var_code = ALGO->n_fft * LWE->sigma_s * LWE->sigma_s * pow(ALGO->p , 2.0 - 2.0 * ALGO->n_code/ALGO->n_fft) / (2.0 * M_PI * M_E);
	double Dcode = exp(4 * var_code * pow(M_PI / ALGO->p , 2));
	if (ALGO->n_code > ALGO->n_fft - PRECISION){
		ALGO->n_code = ALGO->n_fft;
		var_code = 0.0;
		Dcode = 1.0;
	}

	double var_round = ALGO->n_fft * LWE->sigma_s * LWE->sigma_s / 12.0;

	double Dround, Darg, Dfpfn;
	if (flag_asympt) {
		Dround = exp(4 * var_round * pow(M_PI / ALGO->p , 2));
		Darg = 0.5;
		Dfpfn = (log2(Nenum) + ALGO->n_code * log(ALGO->p) + log(1.0 / ALGO->mu));
	} else {
		Dround = 1.0;
		for (int i = 0 ; i < LWE->b_s ; ++i){
			double s = CenteredBinomial[LWE->b_s][i][0];
			double p = CenteredBinomial[LWE->b_s][i][1];
			Dround *= pow( sin(M_PI * s / ALGO->p) / (M_PI * s / ALGO->p) , -4.0 * ALGO->n_fft * p);
		}

		/*
		// in Matzov but it seems to be wrong
		double Darg_tmp1 = CenteredBinomial[LWE.b_e][LWE.b_e][1] + exp(-8.0 * M_PI * M_PI * lsigma_s * lsigma_s / (LWE.sigma_e * LWE.sigma_e * LWE.q * LWE.q * (ALGO.m + ALGO.n_lat)));
		double Darg_tmp2 = CenteredBinomial[LWE.b_s][LWE.b_s][1] + exp(-8.0 * M_PI * M_PI * lsigma_s * lsigma_s / (LWE.sigma_s * LWE.sigma_s * LWE.q * LWE.q * (ALGO.m + ALGO.n_lat)));
		Darg = 0.5 * exp(2.0 * pow(Darg_tmp1, ALGO.m) * pow(Darg_tmp2, ALGO.n_lat));
		 */
		Darg = 0.5 + exp(-8.0 * M_PI * M_PI * ( pow(lsigma_s / LWE->q ,2.0) + ((var_round + var_code) / pow(ALGO->p , 2.0)) ) );

		Dfpfn = pow(invCDF( 1.0 - ( ALGO->mu / (2.0 * Nenum * pow(ALGO->p , ALGO->n_code)) ) ) + invCDF(1.0 - (ALGO->mu / 2.0)) , 2.0);
	}

	if (flag_verb){
		printf("lsigma_s^2 = %.12f\n", lsigma_s * lsigma_s);
		printf("var_code = %.12f\n", var_code);
		printf("var_round = %.12f\n", var_round);
		printf("Deq = 2^%f\nDround = 2^%f\nDcode = 2^%f\nDarg = 2^%f\nDfpfn = 2^%f\n\n", log2(Deq), log2(Dround), log2(Dcode), log2(Darg), log2(Dfpfn));
	}
	return Deq * Dround * Dcode * Darg * Dfpfn;
}


/*
 *  @brief Required number of samples to distinguish with advantage with Prange approach.
 *  @param LWE The LWE parameters
 *  @param ALGO The parameters of the algorithm
 *  @param Nenum Number of enumerated vector s_enum
 *  @param distrib Distribution of s_fft
 *  @param flag_asympt Bollean for asymptotic computation or not
 *  @param flag_verb Boolean for verbosity mode
 *  @return Required number of samples to distinguish with Prange approach.
 */
double N_with_bet_f(LWEparams_t *LWE, DualAttackParams_t *ALGO, double Nenum, double **distrib, int flag_asympt, int flag_verb){
	double lsigma_s = (
			pow(LWE->sigma_e , ALGO->m / (ALGO->m + ALGO->n_lat))
			* pow(LWE->sigma_s * LWE->q , ALGO->n_lat / (ALGO->m + ALGO->n_lat))
			* sqrt(4.0 / 3.0)
			* sqrt(ALGO->beta_sieve / (2.0 * M_PI * M_E))
			* pow(delta_f(ALGO->beta_bkz) , ALGO->m + ALGO->n_lat - ALGO->beta_sieve)
	);
	double Deq = exp(4 * pow(lsigma_s * M_PI / LWE->q , 2));

	double squared_norm_s_fft = 0;
	for (int i = 0 ; i < LWE->b_s ; ++i){
		squared_norm_s_fft += (2.0 * distrib[i][1] * distrib[i][0] * distrib[i][0]);
	}
	squared_norm_s_fft *= ALGO->n_fft;

	double var_code = squared_norm_s_fft * pow(ALGO->p , 2.0 - 2.0 * ALGO->n_code/ALGO->n_fft) / (2.0 * M_PI * M_E);
	double Dcode = exp(4 * var_code * pow(M_PI / ALGO->p , 2));
	if (ALGO->n_code > ALGO->n_fft - PRECISION){
		ALGO->n_code = ALGO->n_fft;
		var_code = 0.0;
		Dcode = 1.0;
	}

	double var_round = squared_norm_s_fft / 12.0;

	double Dround, Darg, Dfpfn;
	if (flag_asympt) {
		Dround = exp(4 * var_round * pow(M_PI / ALGO->p , 2));
		Darg = 0.5;
		Dfpfn = (log2(Nenum) + ALGO->n_code * log(ALGO->p) + log(1.0 / ALGO->mu));
	} else {
		Dround = 1.0;
		for (int i = 0 ; i < LWE->b_s ; ++i){
			double s = distrib[i][0];
			double p = distrib[i][1];
			Dround *= pow( sin(M_PI * s / ALGO->p) / (M_PI * s / ALGO->p) , -4.0 * ALGO->n_fft * p);
		}

		Darg = 0.5 + exp(-8.0 * M_PI * M_PI * ( pow(lsigma_s / LWE->q ,2.0) + ((var_round + var_code) / pow(ALGO->p , 2.0)) ) );

		Dfpfn = pow(invCDF( 1.0 - ( ALGO->mu / (2.0 * Nenum * pow(ALGO->p , ALGO->n_code)) ) ) + invCDF(1.0 - (ALGO->mu / 2.0)) , 2.0);
	}

	if (flag_verb){
		printf("lsigma_s^2 = %.12f\n", lsigma_s * lsigma_s);
		printf("var_code = %.12f\n", var_code);
		printf("var_round = %.12f\n", var_round);
		printf("Deq = 2^%f\nDround = 2^%f\nDcode = 2^%f\nDarg = 2^%f\nDfpfn = 2^%f\n\n", log2(Deq), log2(Dround), log2(Dcode), log2(Darg), log2(Dfpfn));
	}
	return Deq * Dround * Dcode * Darg * Dfpfn;
}


/*
 *  @brief Cost for the short vectors sampling.
 *  @param LWE The LWE parameters
 *  @param ALGO The parameters of the algorithm
 *  @param N Number of short vectors required
 *  @param red_cost_model Model for NNS
 *  @param dimension_for_free_model Model for dimension sieving for solving SVP
 *  @return The cost for short vectors sampling.
 */
double Tsample_f(LWEparams_t *LWE, DualAttackParams_t *ALGO, double N, int red_cost_model, int dimension_for_free_model){
	double beta_eff = beta_eff_f(ALGO->beta_bkz, dimension_for_free_model);
	double N_sieve = floor( pow(1.1547005383792515 , ALGO->beta_sieve) ); // 1.1547005383792515 = sqrt(4/3) --> it is an asymptotic estimation

	double res =  ceil(N / N_sieve) * (
			(C_PROG * Tnns_f(ALGO->beta_sieve , red_cost_model)) +
			(C_PROG * C_PROG * (ALGO->m + ALGO->n_lat - ALGO->beta_bkz + 1) * Tnns_f(beta_eff , red_cost_model))
	);
	//res = ceil(N / N_sieve) * pow(2.0 , 0.2920 * ALGO->beta_sieve); // (C0) Model classical asymptotic sieving
	if (res < 0){
		return INFINITY;
	}
	return res;
}

/*
 *  @brief Running time for the original Matzov Algorithm.
 *  		fixe ALGO->n_code == ALGO->n_fft for considering no code
 *  @param LWE The LWE parameters
 *  @param ALGO The parameters of the algorithm
 *  @param red_cost_model Model for NNS
 *  @param dimension_for_free_model Model for dimension sieving for solving SVP
 *  @param flag_asympt Bollean for asymptotic computation or not
 *  @param flag_verb Boolean for verbosity mode
 *  @return Running time for the original Matzov Algorithm.
 */
double Cost_codedMatzov(LWEparams_t *LWE, DualAttackParams_t *ALGO, int red_cost_model, int dimension_for_free_model, int flag_asympt, int flag_verb){
	ALGO->n_lat = LWE->n - ALGO->n_enum - ALGO->n_fft;
	if ((ALGO->m < 1) || (ALGO->m > LWE->m)) {
		ALGO->m = LWE->m;
	}

	// Correction of Martin Albretch for the cost of the enumeration
	double coeff = 1.0 / (1.0 - exp(-1.0 / (2.0 * LWE->sigma_s * LWE->sigma_s)));
	double tmp_alpha = M_PI * M_PI * LWE->sigma_s * LWE->sigma_s;
	double tmp_a = exp(8.0 * tmp_alpha * exp(-2.0 * tmp_alpha) * tanh(tmp_alpha));
	double Nenum = coeff * pow(2 * tmp_a / sqrt(M_E) , ALGO->n_enum) * pow(2 , ALGO->n_enum * H_f(LWE->sigma_s));
	// Nenum = pow(2 , ALGO->n_enum * H_f(LWE->sigma_s)); // in Matzov

	double N = N_f(LWE, ALGO, Nenum, flag_asympt, flag_verb);

	double T_sample = Tsample_f(LWE, ALGO, N, red_cost_model, dimension_for_free_model);
	double T_fft = Tfft_f(ALGO->p, ALGO->n_code) + Ttable_f(N);
	double T_decode = 1.0;
	if (ALGO->n_fft > ALGO->n_code) {
		T_decode = Tdecode_f(N, ALGO->p, ALGO->n_fft);
	} else {
		ALGO->n_code = ALGO->n_fft;
	}
	double T_guess = Nenum * (T_fft + T_decode);
	// T_guess = pow(2 , ALGO.n_enum * H) * (T_fftf(ALGO.n_fft, ALGO.p) + T_tablef(N)); //wrong formula from [Matzov22]

	if (flag_verb){
		printf("Nenum = 2^%f\n", log2(Nenum));
		printf("N = 2^%f\n", log2(N));
		printf("Nsieve = 2^%f\n", log2(floor( pow(1.1547005383792515 , ALGO->beta_sieve) )));
		printf("Tfft = 2^%f\n", log2(T_fft));
		printf("Tdecode = 2^%f\n", log2(T_decode));
		printf("Tsample = 2^%f\n", log2(T_sample));
		printf("Tguess = 2^%f\n", log2(T_guess));
	}
	return log2(T_sample + T_guess);
}



/*
 *  @brief Running time for the original Matzov Algorithm.
 *  		fixe ALGO->n_code == ALGO->n_fft for considering no code
 *  @param LWE The LWE parameters
 *  @param ALGO The parameters of the algorithm
 *  @param red_cost_model Model for NNS
 *  @param dimension_for_free_model Model for dimension sieving for solving SVP
 *  @param flag_asympt Bollean for asymptotic computation or not
 *  @param flag_verb Boolean for verbosity mode
 *  @return Running time for the original Matzov Algorithm.
 */
double Cost_codedMatzov_with_bet(LWEparams_t *LWE, DualAttackParams_t *ALGO, int red_cost_model, int dimension_for_free_model, int flag_asympt, int flag_verb){
	ALGO->n_lat = LWE->n - ALGO->n_enum - ALGO->n_fft;
	if ((ALGO->m < 1) || (ALGO->m > LWE->m)) {
		ALGO->m = LWE->m;
	}

	double Nenum = 1.0;

	double number_zeros = f[LWE->b_s][lround(ALGO->n_enum + ALGO->n_fft)];

	//double number_zeros = floor((ALGO->n_enum + ALGO->n_fft) * CenteredBinomial[LWE->b_s][LWE->b_s][1]);
	if (number_zeros < ALGO->n_enum){
		return INFINITY;
	}

	double **distrib = (double **)malloc( (2*LWE->b_s + 1) * sizeof(double *) );
	//double rate = (ALGO->n_enum + ALGO->n_fft) / ALGO->n_fft;
	for (int i = 0 ; i < LWE->b_s ; ++i){
		distrib[i] = (double *)malloc( 2 * sizeof(double) );
		distrib[2*LWE->b_s - i] = (double *)malloc( 2 * sizeof(double) );
		distrib[i][0] = CenteredBinomial[LWE->b_s][i][0];
		distrib[2*LWE->b_s - i][0] = CenteredBinomial[LWE->b_s][2*LWE->b_s - i][0];
		//distrib[i][1] = rate * CenteredBinomial[LWE->b_s][i][1];
		//distrib[2*LWE->b_s - i][1] = rate * CenteredBinomial[LWE->b_s][2*LWE->b_s - i][1];
		distrib[i][1] = CenteredBinomial[LWE->b_s][i][1];
		distrib[2*LWE->b_s - i][1] = CenteredBinomial[LWE->b_s][2*LWE->b_s - i][1];
	}

	double Psucc = 1.0;
	for (double i = 1 ; i <= ALGO->n_enum + PRECISION ; i += 1.0){
		Psucc *= ((number_zeros + 1 - i) / (ALGO->n_fft + i));
	}

	ALGO->mu = (ALGO->nu * Psucc)/((ALGO->nu * Psucc) + 2.0 - ALGO->nu);
	double R = log(ALGO->nu / 2.0) / (log(1.0 - Psucc) + log(1.0 - ALGO->mu));


	double N = N_with_bet_f(LWE, ALGO, Nenum, distrib, flag_asympt, flag_verb);

	double T_sample = Tsample_f(LWE, ALGO, N, red_cost_model, dimension_for_free_model);
	double T_fft = Tfft_f(ALGO->p, ALGO->n_code) + Ttable_f(N);
	double T_decode = 1.0;
	if (ALGO->n_fft > ALGO->n_code) {
		T_decode = Tdecode_f(N, ALGO->p, ALGO->n_fft);
	} else {
		ALGO->n_code = ALGO->n_fft;
	}

	double T_guess = R * Nenum * (T_fft + T_decode);


	if (flag_verb){
		printf("Psucc = 2^%f\n", log2(Psucc));
		printf("number of zeros = %f\n", number_zeros);
		printf("mu = %f\n", ALGO->mu);
		printf("R = 2^%f\n", log2(R));
		printf("Nenum = 2^%f\n", log2(Nenum));
		printf("N = 2^%f\n", log2(N));
		printf("Nsieve = 2^%f\n", log2(floor( pow(1.1547005383792515 , ALGO->beta_sieve) )));
		printf("Tfft = 2^%f\n", log2(T_fft));
		printf("Tdecode = 2^%f\n", log2(T_decode));
		printf("Tsample = 2^%f\n", log2(T_sample));
		printf("Tguess = 2^%f\n", log2(T_guess));
	}

	for (int i = 0 ; i < LWE->b_s ; ++i){
		free(distrib[i]);
	}
	free(distrib);

	return log2(T_sample + T_guess);
}



double optim_Matzov(LWEparams_t *LWE, DualAttackParams_t *ALGO, int red_cost_model, int dimension_for_free_model, int flag_asympt, double radius, double step){

	double Ccur, Cmin = INFINITY;

	DualAttackParams_t ALGOmin;
	copy_parameters(&ALGOmin, ALGO);


	double m_inf=fmax(2, ALGO->m - radius);
	double m_sup=fmin(LWE->m, ALGO->m + radius);

	double n_enum_inf=fmax(1, ALGO->n_enum - radius);
	double n_enum_sup=fmin(400, ALGO->n_enum + radius);

	double n_fft_inf=fmax(1, ALGO->n_fft - radius);
	double n_fft_sup=ALGO->n_fft + radius;

	double p_inf=fmax(2, ALGO->p - radius);;
	double p_sup=fmin(LWE->q, ALGO->p + radius);

	double beta_bkz_inf=fmax(40, ALGO->beta_bkz - radius);
	double beta_bkz_sup=fmin(LWE->n, ALGO->beta_bkz + radius);

	double beta_sieve_inf=fmax(40, ALGO->beta_sieve - radius);
	double beta_sieve_sup=fmin(LWE->n, ALGO->beta_sieve + radius);

	for (ALGO->m = m_inf ; ALGO->m <= m_sup; ALGO->m += step){
		for (ALGO->n_enum = n_enum_inf ; ALGO->n_enum <= n_enum_sup; ALGO->n_enum += step){
			for (ALGO->n_fft = n_fft_inf ; ALGO->n_fft <= fmin(LWE->n - ALGO->n_enum, n_fft_sup); ALGO->n_fft += step){
				ALGO->n_lat = LWE->n - ALGO->n_enum - ALGO->n_fft;
				for (ALGO->p = p_inf ; ALGO->p <= p_sup; ALGO->p += 1){
					ALGO->n_code = ALGO->n_fft;
					for (ALGO->beta_bkz = beta_bkz_inf ; ALGO->beta_bkz <= beta_bkz_sup; ALGO->beta_bkz += step){
						//beta_sieve_inf = beta_eff_f(ALGO->beta_bkz, dimension_for_free_model);
						//beta_sieve_sup = beta_sieve_inf;
						for (ALGO->beta_sieve = beta_sieve_inf ; ALGO->beta_sieve <= fmin(ALGO->m + ALGO->n_lat, beta_sieve_sup); ALGO->beta_sieve += step){
							Ccur = Cost_codedMatzov(LWE, ALGO, red_cost_model, dimension_for_free_model, flag_asympt, NOVERB);
							if (Ccur < Cmin){
								Cmin = Ccur;
								copy_parameters(&ALGOmin, ALGO);
							}
						}
					}
				}
			}
		}
	}

	copy_parameters(ALGO, &ALGOmin);
	return Cmin;
}

double optim_codedMatzov(LWEparams_t *LWE, DualAttackParams_t *ALGO, int red_cost_model, int dimension_for_free_model, int flag_asympt, double radius, double step){
	double Ccur, Cmin = INFINITY;

	DualAttackParams_t ALGOmin;
	copy_parameters(&ALGOmin, ALGO);


	double m_inf=fmax(2, ALGO->m - radius);
	double m_sup=fmin(LWE->m, ALGO->m + radius);

	double n_enum_inf=fmax(1, ALGO->n_enum - radius);
	double n_enum_sup=fmin(400, ALGO->n_enum + radius);



	double n_fft_inf=fmax(1, ALGO->n_fft - radius);
	double n_fft_sup=ALGO->n_fft + radius;

	double log2_p_inf=1;
	double log2_p_sup=floor(log2(LWE->q));

	double n_code_inf=fmax(1, ALGO->n_code - radius);
	double n_code_sup=ALGO->n_code + radius;

	double beta_bkz_inf=fmax(40, ALGO->beta_bkz - radius);
	double beta_bkz_sup=fmin(LWE->n, ALGO->beta_bkz + radius);

	double beta_sieve_inf=fmax(40, ALGO->beta_sieve - radius);
	double beta_sieve_sup=fmin(LWE->n, ALGO->beta_sieve + radius);

	for (ALGO->m = m_inf ; ALGO->m <= m_sup; ALGO->m += step){
		for (ALGO->n_enum = n_enum_inf ; ALGO->n_enum <= n_enum_sup; ALGO->n_enum += step){
			for (ALGO->n_fft = n_fft_inf ; ALGO->n_fft <= fmin(LWE->n - ALGO->n_enum, n_fft_sup); ALGO->n_fft += step){
				ALGO->n_lat = LWE->n - ALGO->n_enum - ALGO->n_fft;
				for (double log2_p = log2_p_inf ; log2_p <= log2_p_sup; log2_p += 1){
					ALGO->p = pow(2 , log2_p);
					for (ALGO->n_code = n_code_inf ; ALGO->n_code <= fmin(ALGO->n_fft, n_code_sup); ALGO->n_code += step){
						for (ALGO->beta_bkz = beta_bkz_inf ; ALGO->beta_bkz <= fmin(ALGO->m + ALGO->n_lat, beta_bkz_sup); ALGO->beta_bkz += step){
							//beta_sieve_inf = beta_eff_f(ALGO->beta_bkz, dimension_for_free_model);
							//beta_sieve_sup = beta_sieve_inf;
							for (ALGO->beta_sieve = beta_sieve_inf ; ALGO->beta_sieve <= beta_sieve_sup; ALGO->beta_sieve += step){
								Ccur = Cost_codedMatzov(LWE, ALGO, red_cost_model, dimension_for_free_model, flag_asympt, NOVERB);
								if (Ccur < Cmin){
									Cmin = Ccur;
									copy_parameters(&ALGOmin, ALGO);
								}
							}
						}
					}
				}
			}
		}
	}

	copy_parameters(ALGO, &ALGOmin);
	return Cmin;
}

double optim_codedMatzov_with_bet(LWEparams_t *LWE, DualAttackParams_t *ALGO, int red_cost_model, int dimension_for_free_model, int flag_asympt, double radius, double step){
	double Ccur, Cmin = INFINITY;

	DualAttackParams_t ALGOmin;
	copy_parameters(&ALGOmin, ALGO);

	double m_inf=fmax(2, ALGO->m - radius);
	double m_sup=fmin(LWE->m, ALGO->m + radius);

	double n_enum_inf=fmax(1, ALGO->n_enum - radius);
	double n_enum_sup=fmin(400, ALGO->n_enum + radius);

	double n_fft_inf=fmax(1, ALGO->n_fft - radius);
	double n_fft_sup=ALGO->n_fft + radius;

	double log2_p_inf=1;
	double log2_p_sup=floor(log2(LWE->q));

	double n_code_inf=fmax(1, ALGO->n_code - radius);
	double n_code_sup=ALGO->n_code + radius;

	double beta_bkz_inf=fmax(40, ALGO->beta_bkz - radius);
	double beta_bkz_sup=fmin(LWE->n, ALGO->beta_bkz + radius);

	double beta_sieve_inf=fmax(40, ALGO->beta_sieve - radius);
	double beta_sieve_sup=fmin(LWE->n, ALGO->beta_sieve + radius);

	ALGO->nu = 0.5;
	for (ALGO->m = m_inf ; ALGO->m <= m_sup; ALGO->m += step){
		for (ALGO->n_enum = n_enum_inf ; ALGO->n_enum <= n_enum_sup; ALGO->n_enum += fmax(1, step/4)){
			for (ALGO->n_fft = n_fft_inf ; ALGO->n_fft <= fmin(LWE->n - ALGO->n_enum, n_fft_sup); ALGO->n_fft += step){
				ALGO->n_lat = LWE->n - ALGO->n_enum - ALGO->n_fft;
				for (double log2_p = log2_p_inf ; log2_p <= log2_p_sup; log2_p += 1){
					ALGO->p = pow(2 , log2_p);
					for (ALGO->n_code = n_code_inf ; ALGO->n_code <= fmin(ALGO->n_fft, n_code_sup); ALGO->n_code += step){
						for (ALGO->beta_bkz = beta_bkz_inf ; ALGO->beta_bkz <= fmin(ALGO->m + ALGO->n_lat, beta_bkz_sup); ALGO->beta_bkz += step){
							//beta_sieve_inf = beta_eff_f(ALGO->beta_bkz, dimension_for_free_model);
							//beta_sieve_sup = beta_sieve_inf;
							for (ALGO->beta_sieve = beta_sieve_inf ; ALGO->beta_sieve <= beta_sieve_sup; ALGO->beta_sieve += step){
								Ccur = Cost_codedMatzov_with_bet(LWE, ALGO, red_cost_model, dimension_for_free_model, flag_asympt, NOVERB);
								if (Ccur < Cmin){
									Cmin = Ccur;
									copy_parameters(&ALGOmin, ALGO);
								}
							}
						}
					}
				}
			}
		}
	}

	copy_parameters(ALGO, &ALGOmin);
	return Cmin;
}

int main(void) {
	int rc_model = list_decoding_classical;  //list_decoding_naive_classical (CN) , list_decoding_classical (CC)
	int dff_model = G6K_MODEL; // ASYMPTOTIC_MODEL or G6K_MODEL
	int asympt = NOASYMPT;

	DualAttackParams_t ALGOparams = {
			652, //m
			27, //n_enum
			51, //n_fft
			690, //n_lat
			51, //n_code
			5, // p
			580, // beta_bkz
			574, // beta_sieve
			0.5, // mu
			0.5 // nu
	};


	print_LWEparams(&Saber);

	double radius = 128, step = 32;
	double Cmin_pred = INFINITY;
	double Cmin_cur;
	while (step >= 1 - PRECISION){
		Cmin_cur = optim_codedMatzov(&Saber, &ALGOparams, rc_model, dff_model, asympt, radius, step);
		printf("radius = %f ,  step = %f   --> %.12f\n", radius, step, Cmin_cur);
		if (Cmin_cur >= Cmin_pred){
			radius /= 2.0;
			step /= 2.0;
		}
		Cmin_pred = Cmin_cur;
	}
	printf("Complexity = 2^%.12f\n", Cost_codedMatzov(&Saber, &ALGOparams, rc_model, dff_model, asympt, VERB));
	print_DualAttackParams(&ALGOparams);


	/*
	DualAttackParams_t MATZOV_Kyber1024 = {
			905, //m
			32, //n_enum
			72, //n_fft
			690, //n_lat
			72, //n_code
			5, // p
			846, // beta_bkz
			814, // beta_sieve
			0.5, // mu
			0.5 // nu
	};

	DualAttackParams_t MATZOV_Kyber512 = {
			485, //m
			17, //n_enum
			33, //n_fft
			462, //n_lat
			33, //n_code
			5, // p
			379, // beta_bkz
			383, // beta_sieve
			0.5, // mu
			0.5  // nu
	};
	print_LWEparams(&Kyber512);
	printf("with the given parameters : %.12f\n", Cost_codedMatzov(&Kyber512, &MATZOV_Kyber512, rc_model, dff_model, NOASYMPT, NOVERB));
	print_DualAttackParams(&MATZOV_Kyber512);


	DualAttackParams_t MATZOV_Kyber768 = {
			652, //m
			27, //n_enum
			51, //n_fft
			690, //n_lat
			51, //n_code
			5, // p
			580, // beta_bkz
			574, // beta_sieve
			0.5, // mu
			0.5 // nu
	};
	print_LWEparams(&Kyber768);
	printf("with the given parameters : %.12f\n", Cost_codedMatzov(&Kyber768, &MATZOV_Kyber768, rc_model, dff_model, NOASYMPT, NOVERB));
	print_DualAttackParams(&MATZOV_Kyber768);
	*/
	return EXIT_SUCCESS;
}
