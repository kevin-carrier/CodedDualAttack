/*
 * ring.h
 *
 *  Created on: 24 nov. 2022
 *      Author: Kévin Carrier
 *
 * Description: Ring ZZ/qZZ
 */

#ifndef RING_H_
#define RING_H_

#include <stdio.h>
#include <stdlib.h>
#include "tools.h"

typedef struct Ring {
    unsigned long q;
    long *inverts;
} ring_t;

void free_ring(ring_t *ring){
    free(ring->inverts);
    free(ring);
}

ring_t *copy_ring(ring_t *r){
    ring_t *ring = (ring_t *)malloc(sizeof(ring_t));
    ring->q = r->q;
    ring->inverts = (long *)malloc(sizeof(long) * ring->q);
    for (unsigned long x = 0; x < ring->q; ++x){
        ring->inverts[x] = r->inverts[x];
    }
    return ring;
}
ring_t *gen_ring(unsigned long q){
    ring_t *ring = (ring_t *)malloc(sizeof(ring_t));
    ring->q = q;
    ring->inverts = (long *)malloc(sizeof(long) * ring->q);
    for (unsigned long x = 0; x < ring->q; ++x){
        ring->inverts[x] = inv_mod(x, q);
    }
    return ring;
}

int is_invertible(ring_t *ring, unsigned long x){
    if (ring->inverts[x] < 0){
        return 0;
    } else {
        return 1;
    }
}

unsigned long invert(ring_t *ring, unsigned long x){
    if (ring->inverts[x] < 0){
        fprintf(stderr, "invert error: %lu cannot be inverted!\n", x);
        exit(EXIT_FAILURE);
    }
    return ring->inverts[x];
}

#endif /* RING_H_ */
