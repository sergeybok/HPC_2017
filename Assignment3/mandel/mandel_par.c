#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>
#include <time.h>
#include <sys/time.h>

#include "pngwriter.h"
#include "consts.h"

#include <math.h>


unsigned long get_time ()
{
    struct timeval tp;
    gettimeofday (&tp, NULL);
    return tp.tv_sec * 1000000 + tp.tv_usec;
}

int main (int argc, char** argv)
{
	png_data* pPng = png_create (IMAGE_WIDTH, IMAGE_HEIGHT);
	
	double z_re, z_im, z_re2, z_im2, cx, cy;
	cy = MIN_Y;
	
	double fDeltaX = (MAX_X - MIN_X) / (double) IMAGE_WIDTH;
	double fDeltaY = (MAX_Y - MIN_Y) / (double) IMAGE_HEIGHT;
	
	long nTotalIterationsCount = 0;
	unsigned long nTimeStart = get_time ();
	
	long i, j, n;
	long nTotalIterationsCount_local = 0;
        
    n = 0;
    //printf(MAX_ITERS);
	// do the calculation
	
	#pragma omp parallel firstprivate(cx,cy,i,j,n,z_re,z_im,z_re2,z_im2,nTotalIterationsCount_local) shared(pPng,nTotalIterationsCount,fDeltaX,fDeltaY) 
	{
	#pragma omp for
	for (j = 0; j < IMAGE_HEIGHT; j++)
	{
		cy = MIN_Y + fDeltaY*j;
		for (i = 0; i < IMAGE_WIDTH; i++)
		{
			cx = MIN_X + fDeltaX*i;

			z_re = cx;
			z_im = cy;
			z_re2 = z_re * z_re;
			z_im2 = z_im * z_im;
			for (n=0; n < MAX_ITERS && z_re2 + z_im2 < 4; n++) {
				z_im = z_re*z_im*2 + cy;
				z_re = z_re2 - z_im2 + cx;
				z_re2 = z_re * z_re;
				z_im2 = z_im * z_im;
			}

			// plot the number of iterations at point (i, j)
			int c = ((long) n * 255) / MAX_ITERS;
			
			#pragma omp critical
			png_plot (pPng, i, j, c, c, c);

			nTotalIterationsCount_local += n;
			
		}

	}
	#pragma omp critical
	nTotalIterationsCount += nTotalIterationsCount_local;
	}
	
	unsigned long nTimeEnd = get_time ();
	
	// print benchmark data
	printf ("Total time:                 %g ms\n", (nTimeEnd - nTimeStart) / 1000.0);
	printf ("Image size:                 %ld x %ld = %ld Pixels\n",
		(long) IMAGE_WIDTH, (long) IMAGE_HEIGHT, (long) (IMAGE_WIDTH * IMAGE_HEIGHT));
	printf ("Total number of iterations: %ld\n", nTotalIterationsCount);
	printf ("Avg. time per pixel:        %g µs\n", (nTimeEnd - nTimeStart) / (double) (IMAGE_WIDTH * IMAGE_HEIGHT));
	printf ("Avg. time per iteration:    %g µs\n", (nTimeEnd - nTimeStart) / (double) nTotalIterationsCount);
	printf ("Iterations/second:          %g\n", nTotalIterationsCount / (double) (nTimeEnd - nTimeStart) * 1e6);
	// assume there are 8 floating point operations per iteration
	printf ("MFlop/s:                    %g\n", nTotalIterationsCount * 8.0 / (double) (nTimeEnd - nTimeStart));
	
	png_write (pPng, "mandel.png");
	return 0;
}
