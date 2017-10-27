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
#  define BLOCK_SIZE ((unsigned) 20)
#endif





void square_naive (const double *A, const double *B,
	      double *C,
	      const unsigned dim, unsigned s)
{
    unsigned  i, j, k, k_ind;
    double cij;

    for (i = 0; i < s; i++)
    {
    	const double *A_ix = A + i*dim;
    	for (j = 0; j < s; j++)
	    {
            const double *B_xj = B + j;

            cij = C[i*dim + j];
	        k_ind = 0;

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
	unsigned i,j,k,s;
	
	if (BLOCK_SIZE < M) s = (BLOCK_SIZE);
	else s = M;
	unsigned rem = (M%s);

	for (i =0; i < M/s; i++) {
		for (j =0; j < M/s; j++) {
			double *C_ij = C + i*s*M + j*s;
			
			for (k =0; k < M/s; k++) {
				const double *A_ik = A + i*s*M + k*s;
				const double *B_kj = B + k*s*M + j*s;
				square_naive(A_ik, B_kj,C_ij,M,s);
			}
			if (rem ==0) continue;
			//for (int ii = i*s; ii < (i+1)*s; ii++ ) {
			for (int ii = 0; ii < s; ii++){
				const double *A_ix = A + i*s*M + ii*M;
				//for (int jj = j*s; jj < (j+1)*s; jj++) {
				for (int jj = 0; jj < s; jj++){
					const double *B_xj = B + j*s + jj;
					double cij = C_ij[ii*M + jj];
					unsigned k_ind = (M-rem)*M;
					for (int kk = M-rem; kk<M; kk++ ){
						cij += A_ix[kk]*B_xj[k_ind];
						k_ind += M;
					}
					C_ij[ii*M + jj] = cij;
				}
			}
		}
	}

	// Remainder
	s = (M%s);
	if (s==0) {
		return;
	}

	for (i=M-s; i < M;i++){
		const double *A_ix = A + i*M;
		for(j=0; j<M; j++){
			const double *B_xj = B+j;
			double cij = C[i*M +j];
			unsigned k_ind =0;
			for(k=0; k<M; k++) {
				cij += A_ix[k] * B_xj[k_ind];
				k_ind += M;
			}

			C[i*M + j] = cij;
		}
	}
	for (i=0; i < M-s; i++) {
		const double *A_ix = A + i*M;
		for(j=M-s; j< M; j++) {
			const double *B_xj = B + j;
			double cij = C[i*M + j];
			unsigned k_ind = 0;
			for (k=0; k<M; k++) {
				cij += A_ix[k] * B_xj[k_ind];
				k_ind += M;
			}
			C[i*M + j] = cij;
		}
	}

}



