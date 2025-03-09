/*
 * gaussian_channel.h
 *
 *  Created on: 12 may 2024
 *      Author: Kévin Carrier
 *
 * Description: Simulation of a Gaussian channel
 */

#ifndef GAUSSIAN_CHANNEL_H_
#define GAUSSIAN_CHANNEL_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"

#define RAND_PRECISION 1000000
#define LAYOVER 1000

typedef struct Channel {
    double std_dev;
    unsigned long q;
    double *distrib;
    unsigned long *rand_set;
} channel_t;

void free_channel(channel_t *channel){
    free(channel->distrib);
    free(channel->rand_set);
    free(channel);
}

channel_t *gen_channel(unsigned long q, double std_dev){
    channel_t *channel = (channel_t *)malloc(sizeof(channel_t));
    channel->q = q;
    channel->std_dev = std_dev;
    channel->distrib = (double *)malloc(sizeof(double) * q);
    double prob;
    long layover = (long)ceil(LAYOVER/(double)q);

    double variance = std_dev * std_dev;

    double norm_factor = 0;
    for ( unsigned long u = 0 ; u < q ; ++u){
        channel->distrib[u] = 0;
        for ( long x = -layover ; x <= layover ; ++x){
            for (double step = -0.5 ; step < 0.5 ; step += 0.01){
                channel->distrib[u] += exp(-((double)u + step + (double)x*q) * ((double)u + step + (double)x*q)/(2 * variance));
            }
        }
        channel->distrib[u] /= 100.0;
        norm_factor += channel->distrib[u];
    }
    for ( unsigned long u = 0 ; u < q ; ++u){
        channel->distrib[u] /= norm_factor;
    }


    channel->rand_set = (unsigned long *)malloc(sizeof(unsigned long)*RAND_PRECISION);

    unsigned long value = 0;
    double bound = RAND_PRECISION * channel->distrib[value];

    for (long i = 0 ; i < RAND_PRECISION ; ){
        if (i > bound){
            value += 1;
            if (value >= q){
                break;
            }
            bound += RAND_PRECISION * channel->distrib[value];
        } else {
            channel->rand_set[i++] = value;
        }
    }
    return channel;
}

void print_channel(FILE *fd, channel_t *channel){
    double sum_probs = 0;

    fprintf(fd, "Distribution:\n");
    for (long i = 0 ; i < channel->q ; ++i){
        sum_probs += channel->distrib[i];
        printf("[%ld , %.12f],\n", i, channel->distrib[i]);
    }

    fprintf(fd, "\n\nSum of probabilities: %.12f\n\n", sum_probs);
}

unsigned long gen_error(channel_t *channel){
    return channel->rand_set[rand()%RAND_PRECISION];
}


double gen_error_soft(channel_t *channel){
    double e = fmod(rand_normal(0, channel->std_dev), channel->q);
    return (e < 0 ? channel->q + e : e);
}

#endif /* GAUSSIAN_CHANNEL_H_ */
