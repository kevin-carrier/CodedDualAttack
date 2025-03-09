/*
 ============================================================================
 Name        : PolarCode_over_Zq.c
 Author      : Kévin Carrier
 Version     :
 Copyright   : open source
 Description : Distortion of a Polar Code over ZZ/qZZ.
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "polar_code.h"


/*
*    @brief Test of Probabilistic Success Cancellation decoder
*    @param q size of the ring ZZ_q
*    @param n length of the polar code
*    @param k dimension of the polar code
*    @param list_size size of the decoding list
*/
void psc_list_decoder_test(unsigned long q, unsigned long n, unsigned long k, unsigned long list_size){
    printf("Generating the polar code...\n");
    polar_t *code = gen_polar(q, n, k, 100, NULL);
    print_polar(stdout, code);
    printf("\n\n");

    printf("number of punctured positions: %lu\n\n", code->nb_punc);

    long *word = (long *)malloc(sizeof(long) * code->n);
    long *codeword = (long *)malloc(sizeof(long) * code->n);
    long *received = (long *)malloc(sizeof(long) * code->n); // negative integer for the punctured positions
    double **probs = (double **)malloc(sizeof(double *) * code->n);

    /*
     * generate a random word then applying the U|U+V tree
     */
    for (unsigned long i = 0 ; i < code->n ; ++i){
        if (code->frozen_bits[i] == -1){
            word[i] = rand_lu(code->ring->q);
        } else {
            word[i] = code->frozen_bits[i];
        }
    }
    printf("    word =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", word[i]);
    }
    printf("\n");

    memcpy((void *)codeword, (void *)word, sizeof(long) * code->n);
    for (long i = code->m - 1 ; i >= 0 ; --i){
        for (unsigned long j = 0 ; j < (unsigned long)(1 << i) ; ++j){
            unsigned long n_tmp = (unsigned long)(1 << (code->m - i - 1));
            for (unsigned long jj = j * 2 * n_tmp ; jj < (j * 2 * n_tmp + n_tmp) ; ++jj){
                codeword[jj] = ( ( codeword[jj] + codeword[jj + n_tmp] ) % code->ring->q );
                codeword[jj + n_tmp] = ( ( code->coeffs[i][j] * codeword[jj + n_tmp] ) % code->ring->q );
            }
        }
    }
    printf("codeword =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", codeword[i]);
    }
    printf("\n");

    /*
     * Gaussian noise + puncturing
     */
    for (unsigned long i = 0 ; i < code->nb_punc ; ++i){
        received[i] = -1;
    }
    for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
        received[i] = ( (codeword[i] + gen_error(code->channel)) % code->ring->q );
    }


    printf("received =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", received[i]);
    }
    printf("\n");


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
    long *decoded = psc_list_decoder(code, probs, received, list_size);

    printf(" decoded =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", decoded[i]);
    }
    printf("\n\n");

    printf("dist(codeword, received) = %.12f\n", dist(codeword, received, code->ring->q, code->n));
    printf("dist(decoded,  received) = %.12f\n", dist(decoded, received, code->ring->q, code->n));



    for (unsigned long i = 0 ; i < code->n ; ++i){
        free(probs[i]);
    }
    free(probs);
    free(word);
    free(codeword);
    free(received);
    free(decoded);
}




/*
*    @brief Test of Probabilistic Success Cancellation decoder
*    @param q size of the ring ZZ_q
*    @param n length of the polar code
*    @param k dimension of the polar code
*    @param list_size size of the decoding list
*/
void psc_list_decoder_soft_test(unsigned long q, unsigned long n, unsigned long k, unsigned long list_size){
    printf("Generating the polar code...\n");
    polar_t *code = gen_polar(q, n, k, 100, NULL);
    print_polar(stdout, code);
    printf("\n\n");

    printf("number of punctured positions: %lu\n\n", code->nb_punc);

    long *word = (long *)malloc(sizeof(long) * code->n);
    long *codeword = (long *)malloc(sizeof(long) * code->n);
    double *randword = (double *)malloc(sizeof(double) * code->n);
    double *received = (double *)malloc(sizeof(double) * code->n); // negative integer for the punctured positions
    double **probs = (double **)malloc(sizeof(double *) * code->n);

    /*
     * generate a random word then applying the U|U+V tree
     */
    for (unsigned long i = 0 ; i < code->n ; ++i){
        if (code->frozen_bits[i] == -1){
            word[i] = rand_lu(code->ring->q);
        } else {
            word[i] = code->frozen_bits[i];
        }
    }
    printf("    word =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", word[i]);
    }
    printf("\n");

    memcpy((void *)codeword, (void *)word, sizeof(long) * code->n);
    for (long i = code->m - 1 ; i >= 0 ; --i){
        for (unsigned long j = 0 ; j < (unsigned long)(1 << i) ; ++j){
            unsigned long n_tmp = (unsigned long)(1 << (code->m - i - 1));
            for (unsigned long jj = j * 2 * n_tmp ; jj < (j * 2 * n_tmp + n_tmp) ; ++jj){
                codeword[jj] = ( ( codeword[jj] + codeword[jj + n_tmp] ) % code->ring->q );
                codeword[jj + n_tmp] = ( ( code->coeffs[i][j] * codeword[jj + n_tmp] ) % code->ring->q );
            }
        }
    }
    printf("codeword =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", codeword[i]);
    }
    printf("\n");

    /*
     * Gaussian noise + puncturing
     */
    for (unsigned long i = 0 ; i < code->nb_punc ; ++i){
        received[i] = -1;
    }
    for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
        received[i] = fmod( (double)(codeword[i]) + gen_error_soft(code->channel) , code->ring->q );
    }


    printf("received =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%f\t", received[i]);
    }
    printf("\n");


    /*
     * Probs received
     */
    double p_tmp;
    for (unsigned long i = 0 ; i < code->n ; ++i){
        probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
        if (received[i] < 0){
            for (unsigned long u = 0 ; u < code->ring->q ; ++u){
                probs[i][u] = 1.0/(double)(code->ring->q);
            }
        } else {
            for (unsigned long u = 0 ; u < code->ring->q ; ++u){
                p_tmp = fmod(received[i], 1);
                probs[i][u] = (1.0 - p_tmp) * code->channel->distrib[(code->ring->q + (long)floor(received[i]) - u) % code->ring->q]
                              + p_tmp * code->channel->distrib[(code->ring->q + (long)floor(received[i]) + 1 - u) % code->ring->q];
            }
        }

    }


    /*
     * List decoding using Probabilistic Success Cancelation
     */
    long *decoded = psc_list_decoder_soft(code, probs, received, list_size);

    printf(" decoded =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", decoded[i]);
    }
    printf("\n\n");

    printf("dist(codeword, received) = %.12f\n", dist_double(codeword, received, code->ring->q, code->n));
    printf("dist(decoded,  received) = %.12f\n", dist_double(decoded, received, code->ring->q, code->n));


    // generate a punctured random word
    for (unsigned long i = 0 ; i < code->nb_punc ; ++i){
        randword[i] = -1;
    }
    for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
        randword[i] = rand_double(0, code->ring->q);
    }

    printf("\nrandword =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%f\t", randword[i]);
    }
    printf("\n");

    // Probs
    for (unsigned long i = 0 ; i < code->nb_punc ; ++i){
        probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
        for (unsigned long u = 0 ; u < code->ring->q ; ++u){
            probs[i][u] = 1.0 / (double)(code->ring->q);
        }
    }

    for (unsigned long i = code->nb_punc ; i < code->n ; ++i){
        probs[i] = (double *)malloc(sizeof(double) * code->ring->q);
        for (unsigned long u = 0 ; u < code->ring->q ; ++u){
            p_tmp = fmod(randword[i], 1);
            probs[i][u] = (1.0 - p_tmp) * code->channel->distrib[(code->ring->q + (long)floor(randword[i]) - u) % code->ring->q]
                                    + (p_tmp) * code->channel->distrib[(code->ring->q + (long)floor(randword[i]) + 1 - u) % code->ring->q];
        }
    }


    // List decoding using Probabilistic Success Cancelation
    long *decoded_rand = psc_list_decoder_soft(code, probs, randword, list_size);

    printf("decoded =\t");
    for (unsigned long i = 0 ; i < code->n ; ++i){
        printf("%ld\t", decoded_rand[i]);
    }
    printf("\n");

    printf("dist(decoded,  randword) = %.12f\n", dist_double(decoded_rand, randword, code->ring->q, code->n));

    for (unsigned long i = 0 ; i < code->n ; ++i){
        free(probs[i]);
    }
    free(probs);
    free(word);
    free(codeword);
    free(randword);
    free(received);
    free(decoded);
    free(decoded_rand);
}

unsigned long q = 0, n = 0, k = 0 ;

unsigned long nb_samples = 100;
unsigned long nb_codes = 100;
unsigned long nb_randwords = 100;
int flag_real = 0, flag_mse = 0;

int parse_cmdline(int argc, char *argv[]){
    int i;

    for (i = 1 ; i < argc ; ++i){
        if (argv[i][0] == '-') {
            switch (argv[i][1]){
            case 'q':
                if (!argv[i][2]){
                    ++i;
                    q = atoi(argv[i]);
                }else{
                    q = atoi(argv[i]+2);
                }
                break;
            case 'n':
                if (!argv[i][2]){
                    ++i;
                    n = atoi(argv[i]);
                }else{
                    n = atoi(argv[i]+2);
                }
                break;
            case 'k':
                if (!argv[i][2]){
                    ++i;
                    k = atoi(argv[i]);
                }else{
                    k = atoi(argv[i]+2);
                }
                break;
            case 'S':
                if (!argv[i][2]){
                    ++i;
                    nb_samples = atoi(argv[i]);
                }else{
                    nb_samples = atoi(argv[i]+2);
                }
                break;
            case 'N':
                if (!argv[i][2]){
                    ++i;
                    nb_codes = atoi(argv[i]);
                }else{
                    nb_codes = atoi(argv[i]+2);
                }
                break;
            case 'D':
                if (!argv[i][2]){
                    ++i;
                    nb_randwords = atoi(argv[i]);
                }else{
                    nb_randwords = atoi(argv[i]+2);
                }
                break;
            case 'R':
                flag_real = 1;
                break;
            case 'E':
                    if (!argv[i][2]){
                        ++i;
                        if (strcmp(argv[i],"MSE") == 0){
                            flag_mse = 1;
                        }
                    }else{
                        if (strcmp(argv[i]+2,"MSE") == 0){
                            flag_mse = 1;
                        }
                    }
                    break;
            case 'h' :
                printf("Usage : %s <arguments>\n", argv[0]);
                printf("\t\"-h\" : help\n");
                printf("\t\"-q\" : q (unsigned long REQUIRED)\n");
                printf("\t\"-n\" : length of the polar code (unsigned long REQUIRED)\n");
                printf("\t\"-k\" : dimension of the polar code (unsigned long REQUIRED)\n");
                printf("\t\"-S\" : number of samples for estimating the frozen positions (unsigned long OPTIONAL : DEFAULT=100)\n");
                printf("\t\"-N\" : number of codes to test (unsigned long OPTIONAL : DEFAULT=100)\n");
                printf("\t\"-D\" : number of random words to decode for each distortion computation (unsigned long OPTIONAL : DEFAULT=100)\n");
                printf("\t\"-R\" : put this option if you want the vectors to decode are real (OPTIONAL)\n");
                printf("\t\"-E\" : MSE for mean-square-error and ME mean-error (string OPTIONAL : DEFAULT='ME')\n");
                return 2;
                break;
            default :
                fprintf(stderr, "Arguments parse error : argument \"-%c\" unknown !\n", argv[i][1]);
                fprintf(stderr, " -h for help !\n\n");
                return 1;
            }
        }
    }
    if (!q || !n || !k){
        fprintf(stderr, "Arguments parse error : \n");
        if (!q){
            fprintf(stderr, "\t argument -q is required !\n");

        }
        if (!n){
            fprintf(stderr, "\t argument -n is required !\n");

        }
        if (-!k){
            fprintf(stderr, "\t argument -k is required !\n");

        }
        fprintf(stderr, " -h for help !\n\n");
        return 1;
    }

    return 0;
}


int main(int argc, char * argv[]){
    srand((unsigned int)time(NULL));
    //polar_t *code = gen_polar(47, 8, 2, 1, NULL);
    //free_polar(code);

    polar_t * c = gen_polar(101,8,2,100,NULL);
    printf("%f\n",mean_error(c,1,100).mean);
    free_polar(c);
    /*
    // parsing arguments
    switch(parse_cmdline(argc, argv)){
    case 0 :
        break;
    case 1 :
        exit(EXIT_FAILURE);
    case 2 :
        exit(EXIT_SUCCESS);
    }

    if (flag_real){
        printf("We decode real vectors\n");
    } else {
        printf("We decode integer vectors\n");
    }
    if (flag_mse){
        printf("We consider the normalized Mean-Square-Error\n");
    } else {
        printf("We consider the Mean-Error\n");
    }

    //psc_list_decoder_soft_test(q, n, k, 4);

    printf("\n");
    stat_t d;
    double d_min = n * q * q;
    polar_t *best_code = NULL;
    polar_t *code = NULL;
    for ( int i = 0 ; i < nb_codes ; ++i){
        if (flag_real && flag_mse){
            code = gen_polar_soft(q, n, k, nb_samples, NULL);
            d = mean_square_error_soft(code, 1, nb_randwords/10);
            printf("mean-square-error and standard deviation of square-error of the code %d: (%.5f, %.5f)\n",i, d.mean, d.confidence);
        } else if (!flag_real && flag_mse){
            code = gen_polar(q, n, k, nb_samples, NULL);
            d = mean_square_error(code, 1, nb_randwords/10);
            printf("mean-square-error and standard deviation of square-error of the code %d: (%.5f, %.5f)\n",i, d.mean, d.confidence);
        } else if (flag_real && !flag_mse){
            code = gen_polar_soft(q, n, k, nb_samples, NULL);
            d = mean_error_soft(code, 1, nb_randwords/10);
            printf("mean-square-error and standard deviation of square-error of the code %d: (%.5f, %.5f)\n",i, d.mean, d.confidence);
        } else {
            code = gen_polar(q, n, k, nb_samples, NULL);
            d = mean_error(code, 1, nb_randwords/10);
            printf("mean-square-error and standard deviation of square-error of the code %d: (%.5f, %.5f)\n",i, d.mean, d.confidence);
        }

        if (d.mean < d_min){
            d_min = d.mean;
            if (best_code != NULL){
                free_polar(best_code);
            }
            best_code = copy_polar(code);
        }
        free_polar(code);
    }

    printf("\n\n ************************ BEGIN: The polar code parameters ************************\n");
    print_polar(stdout, best_code);
    printf("\n ************************  END: The polar code parameters  ************************\n\n");

    double dgv = sqrt((double)n) * pow(q, 1.0 - (double)k/((double)(n))) / sqrt(2.0 * M_PI * M_E);
    double msegv = 1.0 / (2.0 * M_PI * M_E);

    if (flag_real && flag_mse){
        printf("[q, n, k, msegv, list size, mean of square-errors, std-dev of square-errors]\nBEGIN RESULTS\n");
        for (unsigned long list_size = 1 ; list_size < 9 ; list_size *= 2){
            d = mean_square_error_soft(best_code, list_size, nb_randwords);
            printf("[%lu, %lu, %lu, %.12f, %lu, %.12f, %.12f, %.12f]\n", q, n, k, msegv, list_size, d.mean, d.std_dev, d.confidence);
        }
    } else if (!flag_real && flag_mse){
        printf("[q, n, k, msegv, list size, mean of square-errors, std-dev of square-errors]\nBEGIN RESULTS\n");
        for (unsigned long list_size = 1 ; list_size < 9 ; list_size *= 2){
            d = mean_square_error(best_code, list_size, nb_randwords);
            printf("[%lu, %lu, %lu, %.12f, %lu, %.12f, %.12f, %.12f]\n", q, n, k, msegv, list_size, d.mean, d.std_dev, d.confidence);
        }
    } else if (flag_real && !flag_mse){
        printf("[q, n, k, dgv, list size, mean of errors, std-dev of errors]\nBEGIN RESULTS\n");
        for (unsigned long list_size = 1 ; list_size < 9 ; list_size *= 2){
            d = mean_error_soft(best_code, list_size, nb_randwords);
            printf("[%lu, %lu, %lu, %.12f, %lu, %.12f, %.12f, %.12f]\n", q, n, k, dgv, list_size, d.mean, d.std_dev, d.confidence);
        }
    } else {
        printf("[q, n, k, dgv, list size, mean of errors, std-dev of errors]\nBEGIN RESULTS\n");
        for (unsigned long list_size = 1 ; list_size < 9 ; list_size *= 2){
            d = mean_error(best_code, list_size, nb_randwords);
            printf("[%lu, %lu, %lu, %.12f, %lu, %.12f, %.12f, %.12f]\n", q, n, k, dgv, list_size, d.mean, d.std_dev, d.confidence);
        }
    }
    printf("END RESULTS\n");

    free_polar(best_code);
    
    */
    return EXIT_SUCCESS;
}
