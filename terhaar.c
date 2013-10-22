/*
 * terhaar.c
 *
 * driver for a one-dimensional terhaar wavelet
 *
 * basic idea follows the conte-deboor fft (uugh)
 * where a buffer is used so that the bit-swapped output
 * happens automagically
 *
 * UNLIKE THE FFT n really does have to be a power of two
 *
 * rwh 2001
 *
 * 
 */

#include <stdio.h>
#include <stdlib.h>

int  terhaar( float x[], float scratch[], int n)
{

	int i,j,delta,stride;
	float xsum,xdiff;
	
	for( i=2; i <= n ; i *= 2)
	{
		if( i == n) goto N_IS_A_POWER_OF_TWO;
	}	

	fprintf(stderr,"terhaar requires a power of 2 %d is not\n", n);
	return 1==0;

N_IS_A_POWER_OF_TWO: ;

	for( i=1; i<n; i*=2)
	{
		delta = 2;
		stride = n/2;
		for( j=0; j< n/2; j++)
		{
			xsum = x[j+j] + x[j+j+1];
			xdiff = x[j+j] - x[j+j+1];

			scratch[j] = xsum;
			scratch[stride+j] = xdiff;

		}/*inc */
		for( j=0; j< n; j++)
		{
			x[j] = scratch[j];
		}
	} /* i */




}
