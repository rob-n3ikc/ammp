/*
* semi.c
*
*  semi-classical virtual potential
*
*  should be the last potential in the use list.
*
*  as it is a function of the velocities don't use it in
* minimization.
*
*
*  the virtualness is  controlled by ksemi
*
* basically off-load some of the kinetic energy into the potential
* and make the classical path more alike the quantum one.
*
*
*  (c) 2001 RWH
*  GPL yada yada yada.
*
*/

#include "ammp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int v_semi(float *V, float lamda)
{

	float get_f_variable(char *);
	float k,accum,to_kcal;
	int i,numatm, a_number(void);
	ATOM *ap, *a_next(int);

	k = get_f_variable("ksemi"); 
	if( k < 0.) k = 0.;
	numatm = a_number();

	if(numatm < 1) return 0;
	
	to_kcal = 0.5/(4.184*1000.*1000.);

	accum = 0.;
	for( i=0; i< numatm; i++)
	{
		ap = a_next(i);  /* can use a_index in the c++ */	
		accum -= to_kcal*k*ap->mass*(ap->vx*ap->vx);
		accum -= to_kcal*k*ap->mass*(ap->vy*ap->vy);
		accum -= to_kcal*k*ap->mass*(ap->vz*ap->vz);
		
	}
	
	*V += accum;
	
	return 0;

}

/*
*  f = -dV/dx
*
*  V = -k * sum( dx/dt *dx/dt*mass)  (m v**2 )
*
*
*  dV/dx = -2*k* dx/dt * d(dx/dt)/dx
*        = -2*k* dx/dt * d2x/dt2 /(dx/dt)
*        = -2*k*f_x/mass
*
*/

int f_semi( float lambda)
{
	float get_f_variable(char *);
	float k,to_kcal;
	int i,numatm, a_number(void);
	ATOM *ap, *a_next(int);

	k = get_f_variable("ksemi"); 
	if( k < 0.) k = 0.;
	numatm = a_number();

	if(numatm < 1) return 0;

	to_kcal = 0.5/(4.184*1000.*1000.);
	/* use to_kcal*2 */
	to_kcal += to_kcal;

	for( i=0; i< numatm; i++)
	{
		ap = a_next(i);  /* can use a_index in the c++ */	
		ap->fx += to_kcal*k*ap->fx;
		ap->fy += to_kcal*k*ap->fy;
		ap->fz += to_kcal*k*ap->fz;
		
	}


	return 0;
}

