/* ==================================================================== *
 *								        *
 *  block_dgemm.c -- Implemant a block matrix multiplication routine    *
 *                                                                      *
 * ==================================================================== */

#include "square_dgemm.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* block parameter ... */
#ifndef BLOCK_SIZE
#  define BLOCK_SIZE ((unsigned) 1)
#endif





void square_naive (const double *A, const double *B,
	      double *C,
	      const unsigned dim, int s)
{
    unsigned  i, j, k;

    for (i = 0; i < s; i++)
    {
    	const double *A_ix = A + i*dim;
    	for (j = 0; j < s; j++)
	    {
            const double *B_xj = B + j;

            double cij = C[i*dim + j];
	        unsigned k_ind = 0;

            for (k = 0; k < s; k++)
	        {
    	       cij += A_ix[k] * B_xj[k_ind];
	           k_ind += dim;
            }

            C[i*dim + j] = cij;
    	}
    }
}


/**
 *  square_dgemm -- multiply two block matrices A and B adding result to C, result is C = C + A*B
 */
void square_dgemm (const double  *A, const double  *B,  double  *C, const unsigned  M)
{
	int s = BLOCK_SIZE;
	for (int i =0; i < M/s; i++) {
		for (int j =0; j < M/s; j++) {
			double *C_ij = C + i*s*M + j*s;	
			
			for (int k =0; k < M/s; k++) {
				const double *A_ik = A + i*s*M + k*s;
				const double *B_kj = B + k*s*M + j*s;
				square_naive(A_ik, B_kj,C_ij,M,s);
			}
		}
	}
}



