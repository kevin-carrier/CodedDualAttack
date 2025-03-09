#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "polar_code.h"

polar_t* _random_polar_code(int q, int n, int k, unsigned long a, unsigned long b, int list_size, int nb_test_decodage);
void _free_polar_code(polar_t *code);
//double _mean_error(polar_t* code,unsigned long a,unsigned long b);
stat_t _mean_error(polar_t* code,unsigned long a,unsigned long b);
void _gen_random_codeword(polar_t* code, long* codeword);
void _decode(polar_t* code, long* received_punctured, long* decoded_punctured,unsigned long a);
// Matrix is vectorized line by line [line1 line2 ...]
void _generator_matrix(polar_t* code, long* matrix);
