/*
*/
#include <math.h>

int main()
{

	float x[128],scratch[128];
	float y[128];


	int i,n;
	
	n = 16;

	for( i=0; i< n; i+=1)
		{x[i] = cos(2.*3.14159265*i/16 ); 
			x[i] += cos(2.*3.14159265*i/8);}



	for( i=0; i< n; i++)
		printf("%f ",x[i]);
	printf("\n");

	terhaar(x,scratch,n);
	terhaar(y,scratch,n);
	for( i=0; i< n; i++)
	{
		printf("%f ",x[i]);
/*		x[i] *= y[i];
*/
	}
	x[1] = 0.;
	printf("\n");

	terhaar(x,scratch,n);
	for( i=0; i< n; i++)
		printf("%f ",x[i]);
	printf("\n");

}
