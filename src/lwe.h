/*
 * lwe.h
 *
 *  Created on: 6 oct. 2022
 *      Author: Kévin Carrier
 */

#ifndef LWE_H_
#define LWE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct LWEparams{
	double q; // over ZZ/qZZ
	double n; // LWE dimension
	double m; // number of sample max
	double sigma_s; // variance of the secret vector
	double sigma_e; // variance of the error vector
	int b_s; // centered binomial parameter of the secret vector
	int b_e; // centered binomial parameter of the error vector
	char * name;
} LWEparams_t;

LWEparams_t Kyber512 = {3329, 512, 512, 1.224744871391589, 1.224744871391589, 3, 3, "Kyber512"};
LWEparams_t Kyber768 = {3329, 768, 768, 1.00, 1.00, 2, 2, "Kyber768"};
LWEparams_t Kyber1024 = {3329, 1024, 1024, 1.00, 1.00, 2, 2, "Kyber1024"};

LWEparams_t LightSaber = {8192, 512, 512, 1.5811388300841898, 2.00, 5, 8, "LightSaber"};
LWEparams_t Saber = {8192, 768, 768, 1.4142135623730951, 2.00, 4, 8, "Saber"};
LWEparams_t FireSaber = {8192, 1024, 1024, 1.224744871391589, 2.00, 3, 8, "FireSaber"};

void print_LWEparams(LWEparams_t *LWE){
	printf("LWE parameters for %s:\n", LWE->name);
	printf("q = %.0f\n", LWE->q);
	printf("n = %.0f\n", LWE->n);
	printf("m = %.0f\n", LWE->m);
	printf("sigma_s = %.2f\n", LWE->sigma_s);
	printf("sigma_e = %.2f\n", LWE->sigma_e);
	printf("b_s = %d\n", LWE->b_s);
	printf("b_e = %d\n\n", LWE->b_e);
}

#endif /* LWE_H_ */
