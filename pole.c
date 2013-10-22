/*  AMMP pole correction */
/*   add a factored pole term to the system
//  form of
//  pole atom_id  constant  n
//
//   for   (constant *constant)/r^n
//
// only applied when the n values match.

copyright (C) 2006  Robert W Harrison
  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice

*/

#include "ammp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
ATOM *atom;
float k;
int n;
void *next;
}  POLE;

#define poleLONG sizeof(POLE)
static POLE *pole_first = NULL;
static POLE *pole_last = NULL;

int pole( p1, konst, n)
int p1,n;
float konst;
{
	ATOM *ap1, *a_m_serial(int);
	POLE *new;
	ap1 = a_m_serial( p1);
	if( ap1 == (ATOM*)NULL) return 1;
	
	new = pole_first;
	if( new != (POLE*)NULL)
	{
	while( 1)
	{
		if( new->next == NULL) break;
		if( new->atom == ap1 && n == new->n)
		{
			new->k = konst;
			return 1;
		}
		if( new->next == new) break;
		new = new->next;
	}
	}/* new is not defined */
	if( (new = malloc( poleLONG) ) == NULL)
	{
		aaerror("cannot allocate memory in pole\n");
		return 0;
	}
	if( pole_first == NULL) pole_first = new;
	if( pole_last == NULL) pole_last = new;
	new->atom = ap1;
	new->k = konst;
	new->n = n;
	new->next = new;
	pole_last -> next = new;
	pole_last = new;
	return 1;
}/* end of pole */

/*  this gets used a bunch so we'll wrap it
*   remember static means that it only has scope in this file */
/* probably better to use pow( double, double) for this */
static float my_pow( r,n)
float r;
int n;
{
	float accum = 1.;
	float mult ;
	int  i;
	i = n;
	mult = r;
	while( i > 0)
	{
	
		if( (i & 1) != 0)  /* low order bit */
		accum *= mult;
		mult *= mult;
		i /= 2;
	}
	return accum;
}
int v_pole( V, lambda)
float *V;
float lambda; 
{
	POLE *apole;
	POLE *bpole;
	ATOM *ap,*bp;
	

	if( pole_first == NULL ) return 0 ;
	apole = pole_first;
/* check for trivial calculations here */
	bpole = apole->next;
	if( apole == bpole) return 0;
	if( bpole == NULL ) return 0;
/* loop over every pair of poles */
	while( 1==1)
	{
	int i;
	float dx,dy,dz,dr;
	ap = apole->atom;
	bpole = apole->next;
	while( 1==1)
	{	
	bp = bpole->atom;
	if( bpole->n != apole->n) goto SKIP_THIS_ATOM;
	for(  i=0; i< ap->dontuse; i++)
		{
		if( bp == ap->excluded[i]) goto SKIP_THIS_ATOM;
		}

		dx = ap->x - bp->x + lambda*(ap->dx - bp->dx);
		dy = ap->y - bp->y + lambda*(ap->dy - bp->dy);
		dz = ap->z - bp->z + lambda*(ap->dz - bp->dz);
		dr = dx*dx + dy*dy + dz*dz;
		if( dr > 1.)
		{
		dr = sqrt(dr);
		dr = my_pow( dr, apole->n);	
		} else { dr = 1.;}
		*V += apole->k*bpole->k/dr;

SKIP_THIS_ATOM: ;
	if( bpole == bpole->next) break;
	if( bpole->next == NULL) break;
	bpole = bpole->next;
	}/* inner loop */
	if( apole == apole->next) return 0;
	if( apole->next == NULL ) return 0;
	apole = apole->next;
	}

}/* end of v_pole */
int f_pole( lambda)
float lambda; 
{
	POLE *apole;
	POLE *bpole;
	ATOM *ap,*bp;
	float accum;
	accum = 0.;

	if( pole_first == NULL ) return 0.;
	apole = pole_first;
/* check for trivial calculations here */
	bpole = apole->next;
	if( apole == bpole) return 0;
	if( bpole == NULL ) return 0;
/* loop over every pair of poles */
	while( 1==1)
	{
	int i;
	float konst;
	float dx,dy,dz,dr;
	ap = apole->atom;
	bpole = apole->next;
	while( 1==1)
	{	
	bp = bpole->atom;
	if( bpole->n != apole->n) goto SKIP_THIS_ATOM;
	for(  i=0; i< ap->dontuse; i++)
		{
		if( bp == ap->excluded[i]) goto SKIP_THIS_ATOM;
		}

		dx = ap->x - bp->x + lambda*(ap->dx - bp->dx);
		dy = ap->y - bp->y + lambda*(ap->dy - bp->dy);
		dz = ap->z - bp->z + lambda*(ap->dz - bp->dz);
		dr = dx*dx + dy*dy + dz*dz;
		if( dr > 1.)
		{
		dr = sqrt(dr);
		dr = my_pow( dr, apole->n+2);	
		} else { dr = 1.;}
		konst =  apole->n*apole->k*bpole->k/dr;
		if( ap->active)
		{
		ap->fx += dx*konst;
		ap->fy += dy*konst;
		ap->fz += dz*konst;
		}
		if( bp->active)
		{
		bp->fx -= dx*konst;
		bp->fy -= dy*konst;
		bp->fz -= dz*konst;
		}

SKIP_THIS_ATOM: ;
	if( bpole == bpole->next) break;
	if( bpole->next == NULL ) break;
	bpole = bpole->next;
	}/* inner loop */
	if( apole == apole->next) return 0;
	if( apole->next == NULL ) return 0;
	apole = apole->next;
	}

}/* end of f_pole */

int dump_pole( op)
FILE *op;
{
	POLE *apole;
	ATOM *ap;
	if( pole_first == NULL ) return 0 ;
	apole = pole_first;
	while( 1==1)
	{
		ap = apole->atom;
		fprintf(op,"pole %d %f %d;\n", ap->serial, apole->k, apole->n);
		if( apole->next == NULL) return 0;
		if( apole == apole->next) return 0;
		apole = apole->next;
	}
}
