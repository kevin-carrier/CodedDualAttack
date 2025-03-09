/*
 * tools.h
 *
 *  Created on: 12 may 2024
 *      Author: Kévin Carrier
 */

#ifndef TOOLS_H_
#define TOOLS_H_

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <complex.h>

typedef struct Stat {
    double mean;
    double std_dev;
    double confidence;
} stat_t;


unsigned long rand_lu(unsigned long val_max){
    return (((unsigned long)rand() << 32) + (unsigned long)rand()) % val_max;
}

/*
 * @brief generate a random floating point number from min to max
 */
double rand_double(double min, double max)
{
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

double rand_normal(double mu, double sigma) {
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1) {
        call = !call;
        return (mu + sigma * (double) X2);
    }

    do {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    } while (W >= 1 || W == 0);

    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (double) X1);
}

/*
 * @brief Compute the Euclidean norm by ignoring the punctured positions
 */
double norm(long *v, unsigned long q, unsigned long length){
    double res = 0.0, tmp;
    for (unsigned long i = 0 ; i < length ; ++i){
        if (v[i] >= 0){
            tmp = (v[i] > q/2.0 ? (double)(v[i]) - (double)q : (double)(v[i]));
            res += ( tmp * tmp );
        }
    }
    return sqrt(res);
}


/*
 * @brief Compute the Euclidean distance by ignoring the punctured positions
 */
double dist(long *u, long *v, unsigned long q, unsigned long length){
    double res = 0.0, tmp;
    for (unsigned long i = 0 ; i < length ; ++i){
        if ((v[i] >= 0) && (u[i] >= 0)){
            tmp = fabs((double)(v[i]) - (double)(u[i]));
            if (tmp > q/2.0){
                tmp -= (double)q;
            }
            res += ( tmp * tmp );
        }
    }
    return sqrt(res);
}


/*
 * @brief Compute the Euclidean distance by ignoring the punctured positions. Contrary to dist function, here, one vector can be real.
 */
double dist_double(long *u, double *v, unsigned long q, unsigned long length){
    double res = 0.0, tmp;
    for (unsigned long i = 0 ; i < length ; ++i){
        if ((v[i] >= 0) && (u[i] >= 0)){
            tmp = fabs(v[i] - (double)(u[i]));
            if (tmp > q/2.0){
                tmp -= (double)q;
            }
            res += ( tmp * tmp );
        }
    }
    return sqrt(res);
}


long partition(double *arr, long l, long u) {
    long i, j;
    double v, tmp;
    v = arr[l];
    i = l;
    j = u + 1;
    do {
        do {
            ++i;
        } while ( (arr[i] < v) && (i <= u) );
        do {
            --j;
        } while (v < arr[j]);
        if (i < j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    } while (i < j);
    arr[l] = arr[j];
    arr[j] = v;
    return j;
}

void quick_sort(double *arr, long l, long u){
    long j;
    if (l < u) {
        j = partition(arr, l, u);
        quick_sort(arr, l, j-1);
        quick_sort(arr, j+1, u);
    }
}
/*

// function to swap elements
void swap(double *a, double *b) {
  double t = *a;
  *a = *b;
  *b = t;
}

// function to find the partition position
long partition(double array[], long low, long high) {
  
  // select the rightmost element as pivot
  long pivot = array[high];
  
  // pointer for greater element
  long i = (low - 1);

  // traverse each element of the array
  // compare them with the pivot
  for (long j = low; j < high; j++) {
    if (array[j] <= pivot) {
        
      // if element smaller than pivot is found
      // swap it with the greater element pointed by i
      i++;
      
      // swap element at i with element at j
      swap(&array[i], &array[j]);
    }
  }

  // swap the pivot element with the greater element at i
  swap(&array[i + 1], &array[high]);
  
  // return the partition point
  return (i + 1);
}

void quick_sort(double array[], long low, long high) {
  if (low < high) {
    
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on right of pivot
    long pi = partition(array, low, high);
    
    // recursive call on the left of pivot
    quick_sort(array, low, pi - 1);
    
    // recursive call on the right of pivot
    quick_sort(array, pi + 1, high);
  }
}
*/
long inv_mod(unsigned long x, unsigned long q){
    unsigned long mod = q;
    // extended euclidean algorithm to find Bezout identity q*a + x*b = gcd(x,q)
    long tmp;
    long a_prev = 1;
    long b_prev = 0;
    long a_cur = 0;
    long b_cur = 1;

    long quotient, remainder;
    while (x > 0){
        quotient = q/x;
        remainder = q%x;
        tmp = a_prev - quotient*a_cur;
        a_prev = a_cur;
        a_cur = tmp;

        tmp = b_prev - quotient*b_cur;
        b_prev = b_cur;
        b_cur = tmp;

        q = x;
        x = remainder;
    }
    if (q!=1){
        return -1;
    } else {
        while (b_prev < 0){
            b_prev += mod;
        }
        return b_prev % mod;
    }
}


// Cooley-Tukey method to compute Discrete Fourier Transform
void fft(complex long double * f, unsigned long n){
    if(n != 1){
        long n_over_2 = n/2;
        complex long double * tmp1 = (complex long double*)malloc(sizeof(complex long double)*n_over_2);
        complex long double * tmp2 = (complex long double*)malloc(sizeof(complex long double)*n_over_2);
        for (unsigned long k = 0 ; k < n_over_2 ; ++k){
            tmp1[k] = f[2*k];
            tmp2[k] = f[2*k + 1];
        }
        fft(tmp1, n_over_2);
        fft(tmp2, n_over_2);

        complex long double W=cexp(-2.0*I*M_PI/n), W_tmp=1;
        for(unsigned long k = 0; k<n_over_2 ; ++k){
            f[k+n_over_2] = tmp1[k] - W_tmp * tmp2[k];
            f[k] = tmp1[k] + W_tmp * tmp2[k];
            W_tmp *= W;
        }
        free(tmp1);
        free(tmp2);
    }
}

// Cooley-Tukey method to compute the invert of Discrete Fourier Transform
void ifft_tmp(complex long double * f, unsigned long n){
    if(n != 1){
        unsigned long n_over_2 = n/2;
        complex long double * tmp1 = (complex long double*)malloc(sizeof(complex long double)*n_over_2);
        complex long double * tmp2 = (complex long double*)malloc(sizeof(complex long double)*n_over_2);
        for (unsigned long k = 0 ; k < n_over_2 ; ++k){
            tmp1[k] = f[2*k];
            tmp2[k] = f[2*k + 1];
        }
        ifft_tmp(tmp1, n_over_2);
        ifft_tmp(tmp2, n_over_2);

        complex long double W=cexp(2.0*I*M_PI/n), W_tmp=1;
        for(unsigned long k = 0; k<n_over_2 ; ++k){
            f[k+n_over_2] = tmp1[k] - W_tmp * tmp2[k];
            f[k] = tmp1[k] + W_tmp * tmp2[k];
            W_tmp *= W;
        }
        free(tmp1);
        free(tmp2);
    }
}

void ifft(complex long double * f, unsigned long n){
    ifft_tmp(f, n);
    for (unsigned long k = 0 ; k < n ; ++k){
        f[k] /= (long double)n;
    }
}


/*
 * @brief At the end, out_ is the circular convolution of f and g
 */
void convolution(double * f, double * g, double * out_, unsigned long n){
    if (n & (n - 1)){ // n is not a power of 2
        unsigned long N = (1 << ((unsigned long)ceil(log2(n)) + 1));
        unsigned long i = 0;

        complex long double * f_dft = (complex long double*)calloc(N, sizeof(complex long double));
        complex long double * g_dft = (complex long double*)calloc(N, sizeof(complex long double));
        for ( ; i < n ; ++i){
            f_dft[i] = (complex long double)(f[i]);
            g_dft[i] = (complex long double)(g[i]);
        }

        fft(f_dft, N);
        fft(g_dft, N);

        for ( i = 0 ; i < N ; ++i){
            f_dft[i] *= g_dft[i];
        }

        ifft(f_dft, N);

        for ( i = 0 ; i < n ; ++i){
            out_[i] = creal(f_dft[i]) + creal(f_dft[i + n]);
        }
        free(f_dft);
        free(g_dft);
    } else { // n is a power of 2
        unsigned long i = 0;

        complex long double * f_dft = (complex long double*)malloc(n*sizeof(complex long double));
        complex long double * g_dft = (complex long double*)malloc(n*sizeof(complex long double));
        for ( ; i < n ; ++i){
            f_dft[i] = (complex long double)(f[i]);
            g_dft[i] = (complex long double)(g[i]);
        }

        fft(f_dft, n);
        fft(g_dft, n);

        for ( i = 0 ; i < n ; ++i){
            f_dft[i] *= g_dft[i];
        }

        ifft(f_dft, n);

        for ( i = 0 ; i < n ; ++i){
            out_[i] = creal(f_dft[i]);
        }
        free(f_dft);
        free(g_dft);
    }


}


void convolution_(double * f, double * g, double * out_, unsigned long n){
    unsigned int i, j, u;
    for (u = 0 ; u < n ; ++u){
        out_[u] = 0;
    }
    for (i = 0 ; i < n ; ++i){
        for (j = 0 ; j < n ; ++j){
            u = (i + j)%n;
            out_[u] += f[i]*g[j];
        }
    }
}

#endif /* TOOLS_H_ */
