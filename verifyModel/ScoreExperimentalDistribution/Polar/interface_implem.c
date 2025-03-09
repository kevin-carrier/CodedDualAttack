#include "interface.h"

void _init_random_seed(unsigned int s){
    srand(s);
}
void _init_random(){
    srand(time(NULL));
}

/*
polar_t* _random_polar_code(int q, int n, int k){
    //srand(time(NULL));
	int nb_samples = 100;
	polar_t *code = gen_polar(q, n, k, nb_samples, NULL);
	return code;
}*/
void _free_polar_code(polar_t *code){
	free_polar(code);
}
polar_t *_copy_polar(polar_t* c){
    return copy_polar(c);
}

void _print_polar(polar_t *code){
    printf("q=%lu\n", code->ring->q);
    printf("m=%lu\n", code->m);
    printf("k=%lu\n", code->k);
    printf("nb_punc=%lu\n\n", code->nb_punc);

    printf( "coeffs=\n");
    for (unsigned long i = 0 ; i < code->m - 1 ; ++i){
        for (unsigned long j = 0 ; j < (unsigned long)(1 << i) - 1 ; ++j){
            printf("%lu\t", code->coeffs[i][j]);
        }
        printf( "%lu\n", code->coeffs[i][(unsigned long)(1 << i) - 1]);
    }
    for (unsigned long j = 0 ; j < (code->n / 2) - 1 ; ++j){
        printf( "%lu\t", code->coeffs[code->m - 1][j]);
    }
    printf("%lu\n\n", code->coeffs[code->m - 1][(code->n / 2) - 1]);

    printf( "frozen_probs=");
    for (unsigned long i = 0 ; i < code->n - 1 ; ++i){
        printf( "%.12f\t", code->frozen_probs[i]);
    }
    printf( "%.12f\n", code->frozen_probs[code->n - 1]);

    printf("frozen_bits=");
    for (unsigned long i = 0 ; i < code->n - 1 ; ++i){
        printf( "%d\t", code->frozen_bits[i]);
    }
    printf("%d\n", code->frozen_bits[code->n - 1]);
}

stat_t _mean_error(polar_t *code, unsigned long list_size,unsigned long nb_tests){
    //unsigned long list_size = 1;
    //unsigned long nb_tests = 100;
    //stat_t s = mean_error(code,list_size,nb_tests);
    //return s.mean;
    stat_t s =mean_error(code,list_size,nb_tests);
    return s;
}

void _decode(polar_t* code, long* received_punctured, long* decoded_punctured, unsigned long list_size){

    //int list_size = 1;

    long *received = (long *)malloc(sizeof(long) * code->n); // negative integer for the punctured positions
    double **probs = (double **)malloc(sizeof(double *) * code->n);


    
  	for (unsigned long i = 0 ; i < code->nb_punc ; ++i){
        received[i] = -1;
    }
    for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
        received[i] = received_punctured[i - code->nb_punc];
    }
    /*
     * Probs received
     */
    for (unsigned long i = 0 ; i < code->n ; ++i){
        probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
        if (received[i] < 0){
            for (unsigned long u = 0 ; u < code->ring->q ; ++u){
                probs[i][u] = 1.0/(double)(code->ring->q);
            }
        } else {
            for (unsigned long u = 0 ; u < code->ring->q ; ++u){
                probs[i][u] = code->channel->distrib[(code->ring->q + received[i] - u ) % code->ring->q];
            }
        }
    }


    /*
     * List decoding using Probabilistic Success Cancelation
     */
    long *res_decoded = psc_list_decoder(code, probs, received, list_size);

	for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
		decoded_punctured[i - code->nb_punc] = res_decoded[i];
	}


    for (unsigned long i = 0 ; i < code->n ; ++i){
        free(probs[i]);
    }
    free(probs);
    free(received);
    free(res_decoded);
}

/*
void _decode(polar_t* code, long* received_punctured, long* decoded_punctured){
	int list_size = 1;

	double **probs = (double **)malloc(sizeof(double *) * code->n);

	long *received = (long *)malloc(sizeof(long) * code->n);
  	for (unsigned long i = 0 ; i < code->nb_punc ; ++i){
        received[i] = -1;
    }
    for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
        received[i] = received_punctured[i - code->nb_punc];
    }

    for (unsigned long i = 0 ; i < code->n ; ++i){
        probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
        if (received[i] < 0){
            for (unsigned long u = 0 ; u < code->ring->q ; ++u){
                probs[i][u] = 1.0/(double)(code->ring->q);
            }
        } else {
            for (unsigned long u = 0 ; u < code->ring->q ; ++u){
                probs[i][u] = code->channel->distrib[(code->ring->q + received[i] - u ) % code->ring->q];
            }
        }
    }

	long *res_decoded= psc_list_decoder(code, probs, received, list_size);

	for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
		decoded_punctured[i - code->nb_punc] = res_decoded[i];
	}


	free(res_decoded);
	for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
		free(probs[i]);
	}
	free(probs);
	free(received);

}
*/
void _gen_random_codeword(polar_t* code, long* codeword){
	int length_unpuctured = code->n;
    /*
     * generate a random word then applying the U|U+V tree
     */
	long *codeword_unpuctured = (long *)malloc(sizeof(long) * length_unpuctured);

    for (unsigned long i = 0 ; i < code->n ; ++i){
        if (code->frozen_bits[i] == -1){
            codeword_unpuctured[i] = rand_lu(code->ring->q);
        } else {
            codeword_unpuctured[i] = code->frozen_bits[i];
        }
    }
    for (long i = code->m - 1 ; i >= 0 ; --i){
        for (unsigned long j = 0 ; j < (unsigned long)(1 << i) ; ++j){
            unsigned long n_tmp = (unsigned long)(1 << (code->m - i - 1));
            for (unsigned long jj = j * 2 * n_tmp ; jj < (j * 2 * n_tmp + n_tmp) ; ++jj){
                codeword_unpuctured[jj] = ( ( codeword_unpuctured[jj] + codeword_unpuctured[jj + n_tmp] ) % code->ring->q );
                codeword_unpuctured[jj + n_tmp] = ( ( code->coeffs[i][j] * codeword_unpuctured[jj + n_tmp] ) % code->ring->q );
            }
        }
    }
	
	for(int i = code->nb_punc; i < length_unpuctured; ++i){
		codeword[i - code->nb_punc] = codeword_unpuctured[i];
	}

	free(codeword_unpuctured);

}
// matrix line by line
void _generator_matrix(polar_t* code, long* matrix){
	int length_punctured = code->n - code->nb_punc;
	for(int offset = 0; offset < (code->k)*(length_punctured); offset+= length_punctured){
		_gen_random_codeword(code,&matrix[offset]);
	}
}

polar_t* _random_polar_code(int q, int n, int k,unsigned long test_code, unsigned long nb_samples, int list_size, int nb_test_decodage){
    //srand(time(NULL));
    //int test_code = 1;
    //int nb_samples = 100;

    //int list_size = 1;
    //int nb_test_mean = 100;
    int nb_test_mean = nb_test_decodage;
    polar_t **L_code = (polar_t **)malloc(sizeof(polar_t *) * test_code);
    for (int i = 0; i < test_code; ++i){
        L_code[i] = gen_polar(q, n, k, nb_samples, NULL);
    }

    //double min = _mean_error(L_code[0],list_size,nb_test_mean);
    stat_t _min = _mean_error(L_code[0],list_size,nb_test_mean);
    double min = _min.mean;
    int i_min = 0;
    for (int i = 1; i < test_code; ++i){
        //double d = _mean_error(L_code[i],list_size,nb_test_mean);
        stat_t d_ = _mean_error(L_code[i],list_size,nb_test_mean);
        double d = d_.mean;
        if(d <= min){
            min = d;
            i_min = i;
        }
    }
	
	polar_t *code_min = L_code[i_min];

    for(int i = 0; i < test_code; ++i){
        if(i != i_min){
            free(L_code[i]);
        }
    }
    free(L_code);
	return code_min;
}
/*
int main(){
    polar_t* c = _random_polar_code(47,8,4);
    printf("%f\n",_mean_error(c));
    _free_polar_code(c);
}*/