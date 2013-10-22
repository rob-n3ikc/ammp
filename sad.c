
/* dgeom()
*
* standard distance geometry with the full eigen expansion
* of the whole big matrix.
* 
*  UUGH.
*  maybe if we can get this to work we can go to the power method
*
*/
/* experiments with deflation */
/*
*  copyright 1993,1994,1995,2001 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/

#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"

/*#define RMIN 4. */
#define RMIN 7.
/* remember the RMIN is the SQUARE of the DISTANCE */

/* buggered finite in windows */
#ifdef WINDOWS
#define finite _finite
#endif


/* RMIN is the minimum radius squared */

#define  SQUARE
/* SQUARE uses the square of the xxt matrix and makes the problem
* parabolic (rather than hyperbolic ) */

/* this is based on an older iterative version that never
* quite worked.
* but to save time i've kept the call structure
* -this means it can be added to the other versions of AMMP without
* -a big headache
*/
void sa_dgeom(echo,op, vfs,nfs,niter,origin,eigshift )
FILE *op;
int (*vfs[])(),nfs;
int niter,echo;
int origin;
float eigshift;
{
ATOM *ap,*bp,*a_next(),*aorig,*a_m_serial();
int numatm , a_number();
int  iter,i,j,k,l;
float x,y,z,xo,yo,zo,x1,y1,z1,r;
int v_bond(),v_angle(),u_v_nonbon(),v_noel();
int v_mmbond(),v_mmangle(),v_c_angle(),v_nonbon();
int v_hard(),v_abc();
int v_torsion(), v_step(), v_ttarget();
int v_hybrid(),gsdg_hybrid();

int jacobi(); /* in normal.c */
int a_v_zero();

int use_bond, use_angle, use_noel;
int use_torsion, use_ttarget, use_step;
int use_hybrid;

int local_niter;

int box,itemper;
float temperature,wait,w;
float randf(),randg();
float xs,ys,zs, rold,rnew;

float clone_average_distance(); /* int, int, int, float,float */

float (*matrix)[];
float (*eig)[];
float (*raw)[];
float (*floyd)[]; /* keep the current floyded matrix */
float (*scratch)[]; /* keep the current floyded matrix */
float (*bat)[];/* deflation matrix (die fleidermaus) */
ATOM *(*allatom)[];

numatm = a_number();
if( numatm < 2 ) return; /* the solution for 1 atom is trivial (0,0,0) */

matrix = (float (*)[])malloc( numatm*numatm*sizeof(float));
raw = (float (*)[])malloc( numatm*numatm*sizeof(float));
eig = (float (*)[])malloc( numatm*numatm*sizeof(float));
floyd = (float (*)[])malloc( numatm*numatm*sizeof(float));
scratch = (float (*)[])malloc( numatm*numatm*sizeof(float));
bat = (float (*)[])malloc( numatm*numatm*sizeof(float));

allatom = (ATOM *(*)[])malloc( numatm *sizeof(ATOM *));


for( i=0; i< numatm; i++)
{ (*allatom)[i] = a_next(i);}


use_bond = (1==0);
use_angle = (1==0);
use_noel = (1==0);
use_step = (1==0);
use_torsion = (1==0);
use_ttarget = (1==0);
use_hybrid = (1==0);
for( i=0; i< nfs; i++)
{
	if( vfs[i] == v_bond) use_bond = (1==1);
	if( vfs[i] == v_angle) use_angle = (1==1);
	if( vfs[i] == v_noel) use_noel = (1==1);
	if( vfs[i] == v_step) use_step = (1==1);
	if( vfs[i] == v_ttarget) use_ttarget = (1==1);
	if( vfs[i] == v_torsion) use_torsion = (1==1);
	if( vfs[i] == v_hybrid) use_hybrid = (1==1);
}

	for( i=0; i< numatm*numatm; i++)
		(*bat)[i] = 0.;

	iter = 0;
	for( i=0; i< numatm; i++)
	{
	a_v_zero();
	if( use_bond) gsdg_bond((*allatom)[i]);
	if( use_angle) gsdg_angle((*allatom)[i]);
//	if( use_torsion) gsdg_torsion((*allatom)[i]);
//	if( use_ttarget) gsdg_ttarget((*allatom)[i]);
	if( use_noel) gsdg_noel((*allatom)[i]);
	if( use_step) gsdg_step((*allatom)[i]);

	(*eig)[i*numatm+i]  = 0.;
	(*raw)[i*numatm+i] = 0.;
	ap = (*allatom)[i];
	for( j=i+1; j< numatm; j++)
	{
		bp= (*allatom)[j];
/* if i and j are both inactive calculate the distance and use it */
/* otherwise use the vx value for the j'th atom (which was set by gsdg_) */
	
	if( !ap->active && !bp->active)
	{
		xo = ap->x -bp->x;
		yo = ap->y -bp->y;
		zo = ap->z -bp->z;
		r = xo*xo + yo*yo +zo*zo;
		if( bp->vx > 0.001)
			fprintf(op,"%d %d %f here %f theoretical\n",ap->serial, bp->serial, r, bp->vx);
		iter +=1;
	}else{
		r = bp->vx;

		if( r < 0.001 ) r = 10.e10;
		else iter += 1;
	}
	(*raw)[i*numatm + j] = r;
	(*eig)[i*numatm + j] = r;
	(*raw)[j*numatm + i] = r;
	(*eig)[j*numatm + i] = r;
	}/* j */
	}/* i */
	fprintf(op, "%d distances in %d unique %f fraction\n", iter,
		numatm*(numatm-1)/2,(float)iter/(numatm*(numatm-1)/2));
/*
 * before floyd algorithm
 * we should fill the distances in from the current
 * structure.   
 * floyd will catch inconsistencies between experimental and
 * theoretical distances.
 */
/* setup some defaults
*/
	local_niter = niter;
	if( local_niter < 0) local_niter = -local_niter;
	if( local_niter < 1) local_niter = 100;

	if( niter < 0) 
	{
	niter = -niter;

	for( i=0; i< numatm; i++)
	for( j=i; j< numatm; j++)
	{
		/* skip the theoretical terms */
		if( (*raw)[i*numatm+j] < 10.e9) continue;
		ap = (*allatom)[i];
		bp = (*allatom)[j];
		xo = ap->x -bp->x;
		yo = ap->y -bp->y;
		zo = ap->z -bp->z;
		r = xo*xo + yo*yo +zo*zo;
		(*eig)[i*numatm +j] = r;	
		(*eig)[j*numatm +i] = r;	
			
	}


	for( i=0; i< numatm*numatm; i++)
		(*floyd)[i] = (*eig)[i];

	}else{
#define FLOYD
#ifdef FLOYD
/* now use Floyd's algorithm to fill in the matrix */
	l = 0;
REFLOYD:;
		for( k=0; k< numatm; k++)
		{
	for( i=0; i<numatm; i++)
	for( j=0; j<numatm; j++)
	{

/* if we don't know the distance */	
	if( (*raw)[i*numatm + j] > 10.e6)
	{
		r = (*eig)[i*numatm + j];
		x = (*eig)[i*numatm +k] + (*eig)[k*numatm + j];
		if( x < r ) r = x;
		(*eig)[i*numatm + j] = r;
		(*eig)[j*numatm + i] = r;

	}


	}
		}
	if( l< 1)
	for( i=0; i<numatm*numatm; i++)
	{
		r = (*eig)[i];
		(*floyd)[i] = r;
		if( r > 10.e6)
		{
		l ++; /* put a limit on it */
		goto REFLOYD;
		} 
	}
	for( i=0; i< numatm; i++)
	for( j=i; j< numatm; j++)
	{
		if( (*eig)[numatm*i+j] < 10.e8) continue; /* skip real distances */
		 r = clone_average_distance( i,j,3,7.*7.,.2);
		if( r < 0.0001) continue; /* skip unobserved */
		if( (*floyd)[numatm*i+j] > r) (*floyd)[numatm*i+j] = r;
	}
	
#else
	for( i=0; i< numatm; i++)
	for( j=0; j< numatm; j++)
		if( (*eig)[i*numatm +j] > 10.e6 ) (*eig)[i*numatm + j] = 0.;
	 for( i=0; i< numatm; i++)
		 for( j=0; j< numatm; j++)
		  (*floyd)[i*numatm+j] = (*eig)[i*numatm+j];

#endif
	}
/*
// now we have the distances in raw,  and eig  raw < 0 is unobserved
//
// floyd has appx distances based on raw,eig
//
// scratch and matrix are fossils from dgeom
*/

	for( itemper = 0; itemper < local_niter+1; itemper ++)
	{
	temperature = (float)(local_niter-itemper)*0.01;
//	temperature *= 10;
	temperature *= 2;

	wait = (float)(local_niter - itemper)*0.01 + 1.;
	for( iter=0; iter< (100*numatm < 10000? 100*numatm: 10000); iter ++)
	{
	
//		box = (float)(101-itemper)*0.1*randf()+1;
//		if( box > 5) box = 5;
		box = 1;
//		box = (int)(randf()*10.) + 1;
		if( niter >= 0)
		{
		if( itemper < 4)
			box = (int)(randf()*numatm*0.5) + 1;
		else if( itemper < 6)
			box = (int)(randf()*numatm*0.25) + 1;
		else if( itemper < 8)
			box = (int)(randf()*numatm*0.125) + 1;
		else if( itemper < 10)
			box = (int)(randf()*numatm*0.0625) + 1;
		else if( itemper < 12)
			box = (int)(randf()*numatm*0.03125) + 1;
		else
			box = 1;
		} else {
			box = 1;
		}/* niter */



		i = numatm;
		while(  i == numatm)
			i = numatm*randf();
		
		if( box + i  > numatm) box = numatm -i-1;

		if( use_hybrid)
			for( j=i; j< box+i; j++)
			{
			gsdg_hybrid((*allatom)[j]);		
			}

//		xs = (1.*randf()-.5);
//		ys = (1.*randf()-.5);
//		zs = (1.*randf()-.5);
		xs = randg()*2;
		ys = randg()*2;
		zs = randg()*2;
		rold = 0.;
		rnew = 0.;
		for( j=i;j < i+box; j++)
		{
		ap = (*allatom)[j];
		for( k=0; k< numatm; k++)
		{
		w = 1.;
/* use floyd for the early steps */
//		if( (*raw)[j*numatm + k ] > 10.e7) continue;
//		if( (*floyd)[j*numatm + k ] > 10.e7) continue;
//		if(itemper > 4 && (*raw)[j*numatm + k ] > 10.e7) continue;
//		if( (*raw)[j*numatm + k ] > 10.e7) w = wait;
		bp = (*allatom)[k];
		xo = ap->x -bp->x;
		yo = ap->y -bp->y;
		zo = ap->z -bp->z;
		r  = (xo*xo + yo*yo +zo*zo);
		r -= ((*floyd)[j*numatm + k]);
		if( (*raw)[j*numatm+k] > 10.e7 && r > 0.) w = 0.;
		rold += fabs(r)*w;
		if( k >=i && k < i+box)
		{
			rnew +=fabs(r)*w;
		}else{
		w = 1.;
		xo = ap->x + xs -bp->x;
		yo = ap->y + ys -bp->y;
		zo = ap->z + zs -bp->z;
		r  = (xo*xo + yo*yo +zo*zo);
		r -= ((*floyd)[j*numatm + k]);
		if( (*raw)[j*numatm+k] > 10.e7 && r > 0.) w = 0.;
			rnew += fabs(r)*w ;
		}
		}/* k */
		}/* j */
		if( rnew < 1.e-7) continue;
		if( rnew < rold  || randf() < exp( (rold-rnew)/temperature/box)) 
//		if( rnew < rold  ) 
		{
	//printf("%d %d %f %f\n",itemper,iter,rold,rnew);
		for( j=i; j < i + box; j++)
		{
		ap = (*allatom)[j];
		ap->x += xs;
		ap->y += ys;
		ap->z += zs;
		}
		}
		

	}/* iter */
	}/* itemper */

/* clean up at the end of the routine */
free(allatom);
free(scratch);
free(bat);
free(floyd);
free(raw);
free(eig);
free(matrix);
}/* end of dgeom */

