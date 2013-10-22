/*
 *  polarization - use a  -qa/r^3 term.
 */
/* polar.c
*
* collection of routines to service polared non-locaL nonbonded potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are symmetric - so we don't have to fuck around with
* s matrices and the like.
*
* note that the non-bonded information is in the ATOM structures 
*
*
*/
/*
*  copyright 1992, 1993, 1994 Robert W. Harrison
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
/* ATOM structure contains a serial number for indexing into
* arrays and the like (a Hessian)
* but otherwise is self-contained. Note the hooks for Non-nonboned potentials
*/

#ifdef GRACELESS
#define exp(x) ((x) > -100.? exp(x): 0.)
#endif

/* f_polar()
*
* f_polar increments the forces in the atom structures by the force
* due to the polared-nonbonded components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_polar(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
	float ux,uy,uz;
	float k,r,r4,r0,xt,yt,zt;
	float lcutoff,cutoff,get_f_variable();
	int inbond,inangle,i,test;
	ATOM *a1,*a2,*bonded[10],*angled[10];
	ATOM *a_next( ); /* returns first ATOM when called with -1 */
	int a_number(),inbuffer,imax;
	float (*buffer)[];
	int invector,atomsused,ii,jj;
	float (*vector)[];
	ATOM *(*atms)[],*(*atomall)[];
        float dielectric;
	float fx,fy,fz;
	float xt2,xt3,xt4;
	float yt2,yt3,yt4;
	float zt2,zt3,zt4;

/* first update the lists 
*  this routine checks if any atom has
*   broken the mxdq barrier and updates the
* forces, potentials and expansions thereof */
	fv_update_polar( lambda);

        dielectric = get_f_variable("dielec");
        if( dielectric < 1.) dielectric = 1.;
        dielectric = 332.17752/dielectric;

/*  get the number of atoms and allocate the memory for the array space */
        i = a_number();
        atomall = malloc( i*sizeof(ATOM *) );
        if( atomall == NULL)
        {aaerror("cannot allocate memory in f_nonbon"); return 0;}

        imax = a_number();
        for( i=0; i< imax; i++)
        {
                (*atomall)[i] = a_next(i);
        }
        for( i= 0; i< imax; i++)
        {
	fx = 0.; fy = 0.; fz = 0.;
                a1 = (*atomall)[i];
        xt = a1->dx*lambda +a1->x - a1->px;
        yt = a1->dy*lambda +a1->y - a1->py;
        zt = a1->dz*lambda +a1->z - a1->pz;


        fx = (a1->qxx*xt + a1->qxy*yt
                        + a1->qxz*zt) ;
        fy = (a1->qxy*xt + a1->qyy*yt
                        + a1->qyz*zt) ;
        fz = (a1->qxz*xt + a1->qyz*yt
                        + a1->qzz*zt) ;
#ifdef CUBIC
	xt2 = xt*xt; yt2 = yt*yt; zt2 = zt*zt;
	fx += a1->qxxx*xt2/2. + a1->qxxy*xt*yt + a1->qxxz*xt*zt
		+ a1->qxyy*yt/2. + a1->qxyz*yt*zt + a1->qxzz*zt2/2.;
	fy += a1->qxxy*xt2/2. + a1->qxyy*xt*yt + a1->qxyz*xt*zt
		+ a1->qyyy*yt/2. + a1->qyyz*yt*zt + a1->qyzz*zt2/2.;
	fz += a1->qxxz*xt2/2. + a1->qxyz*xt*yt + a1->qxzz*xt*zt
		+ a1->qyyz*yt/2. + a1->qyzz*yt*zt + a1->qzzz*zt2/2.;
#endif
#ifdef QUARTIC
	xt3 = xt*xt2; yt3 = yt*yt2; zt3 = zt*zt2;
	fx +=  a1->qxxxx*xt3/6. + a1->qxxxy*xt2*yt/2. + a1->qxxxz*xt2*zt/2.
		+ a1->qxxyy*xt*yt/2. + a1->qxxyz*xt*yt*zt + a1->qxxzz*xt*zt2/2.
		+ a1->qxyyy*yt3/6. + a1->qxyyz*yt2*zt/2. + a1->qxyzz*yt*zt2/2.
		+ a1->qxzzz*zt3/6.;
	fy +=  a1->qxxxy*xt3/6. + a1->qxxyy*xt2*yt/2. + a1->qxxyz*xt2*zt/2.
		+ a1->qxyyy*xt*yt/2. + a1->qxyyz*xt*yt*zt + a1->qxyzz*xt*zt2/2.
		+ a1->qyyyy*yt3/6. + a1->qyyyz*yt2*zt/2. + a1->qyyzz*yt*zt2/2.
		+ a1->qyzzz*zt3/6.;
	fz +=  a1->qxxxz*xt3/6. + a1->qxxyz*xt2*yt/2. + a1->qxxzz*xt2*zt/2.
		+ a1->qxyyz*xt*yt/2. + a1->qxyzz*xt*yt*zt + a1->qxzzz*xt*zt2/2.
		+ a1->qyyyz*yt3/6. + a1->qyyzz*yt2*zt/2. + a1->qyzzz*yt*zt2/2.
		+ a1->qzzzz*zt3/6.;
#endif
#ifdef QUINTIC
	xt4 = xt*xt3; yt4 = yt*yt3; zt4 = zt*zt3;
	fx += a1->qxxxxx*xt4/24. + a1->qxxxxy*xt3*yt/6. + a1->qxxxxz*xt3*zt/6.
		+ a1->qxxxyy*xt2*yt2/4. + a1->qxxxyz*xt2*yt*zt/2. + a1->qxxxzz*xt2*zt2/4.
		+ a1->qxxyyy*xt*yt3/6. + a1->qxxyyz*xt*yt2*zt/2. + a1->qxxyzz*xt*yt*zt2/2.
		+ a1->qxxzzz*xt*zt3/6. + a1->qxyyyy*yt4/24. + a1->qxyyyz*yt3*zt/6. 
		+ a1->qxyyzz*yt2*zt2/4. + a1->qxyzzz*yt*zt3/6. + a1->qxzzzz*zt4/24.;
	fy += a1->qxxxxy*xt4/24. + a1->qxxxyy*xt3*yt/6. + a1->qxxxyz*xt3*zt/6.
		+ a1->qxxyyy*xt2*yt2/4. + a1->qxxyyz*xt2*yt*zt/2. + a1->qxxyzz*xt2*zt2/4.
		+ a1->qxyyyy*xt*yt3/6. + a1->qxyyyz*xt*yt2*zt/2. + a1->qxyyzz*xt*yt*zt2/2.
		+ a1->qxyzzz*xt*zt3/6. + a1->qyyyyy*yt4/24. + a1->qyyyyz*yt3*zt/6. 
		+ a1->qyyyzz*yt2*zt2/4. + a1->qyyzzz*yt*zt3/6. + a1->qyzzzz*zt4/24.;
	fz += a1->qxxxxz*xt4/24. + a1->qxxxyz*xt3*yt/6. + a1->qxxxzz*xt3*zt/6.
		+ a1->qxxyyz*xt2*yt2/4. + a1->qxxyzz*xt2*yt*zt/2. + a1->qxxzzz*xt2*zt2/4.
		+ a1->qxyyyz*xt*yt3/6. + a1->qxyyzz*xt*yt2*zt/2. + a1->qxyzzz*xt*yt*zt2/2.
		+ a1->qxzzzz*xt*zt3/6. + a1->qyyyyz*yt4/24. + a1->qyyyzz*yt3*zt/6. 
		+ a1->qyyzzz*yt2*zt2/4. + a1->qyzzzz*yt*zt3/6. + a1->qzzzzz*zt4/24.;
#endif
	a1->fx += fx  + a1->dpx;
	a1->fy += fy  + a1->dpy;
	a1->fz += fz  + a1->dpz;
/* do the close atoms */
	for( jj=0; jj< NCLOSE; jj++)
	{ if( a1->close[jj] == NULL) break; }
        for( ii=0; ii< jj;ii++)
        {
        a2 = a1->close[ii]; 
/* note ux is backwards from below */
        ux = (a2->dx -a1->dx)*lambda +(a2->x -a1->x);
        uy = (a2->dy -a1->dy)*lambda +(a2->y -a1->y);
        uz = (a2->dz -a1->dz)*lambda +(a2->z -a1->z);
        r = ux*ux + uy*uy + uz*uz; r0 = sqrt(r);
        ux = ux/r0; uy = uy/r0; uz = uz/r0;
/* standard point atoms 
        k = -dielectric*a1->q*a2->q/r;
*/
        k = -dielectric*a1->q*a2->q/r;
	r4 = r*r;
	k += 0.5*(a1->a*a2->q*a2->q + a2->a*a1->q*a1->q)/r4/r0*4;
	r = r*r*r;
        k += a1->a*a2->a/r/r0*6;
        k -= a1->b*a2->b/r/r/r0*12;
        a1->fx += ux*k;
        a1->fy += uy*k;
        a1->fz += uz*k;
        a2->fx -= ux*k;
        a2->fy -= uy*k;
        a2->fz -= uz*k;
        }
        } 

	a_inactive_f_zero();
        free( atomall); return 1;

}
/* v_polar()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_polar( V, lambda )
        float *V,lambda;
{
        float r,r4,r0,xt,yt,zt;
        float get_f_variable();
        float k;
        int inbond,inangle,i;
        ATOM *a1,*a2,*bonded[10],*angled[10];
        ATOM *a_next( ); /* returns first ATOM when called with -1 */
        int a_number(),inbuffer;
        int invector,atomsused,ii,jj,imax;
        float (*vector)[];
        float vx;
        float k2;
        ATOM *(*atomall)[];
        float dielectric;
	float xt2,xt3,xt4,xt5;
	float yt2,yt3,yt4,yt5;
	float zt2,zt3,zt4,zt5;


	fv_update_polar( lambda);

        dielectric = get_f_variable("dielec");
        if( dielectric < 1.) dielectric = 1.;
        dielectric = 332.17752/dielectric;

/*  get the number of atoms and allocate the memory for the array space */
        i = a_number();
        atomall = malloc( i*sizeof(ATOM *) );
        if( atomall == NULL)
        {aaerror("cannot allocate memory in v_nonbon"); return 0;}

        imax = a_number();
        for( i=0; i< imax; i++)
        {
                (*atomall)[i] = a_next(i);
        }
        for( i= 0; i< imax; i++)
        {
                a1 = (*atomall)[i];
		vx = a1->VP;
                xt = a1->dx*lambda +a1->x - a1->px;
                yt = a1->dy*lambda +a1->y - a1->py;
                zt = a1->dz*lambda +a1->z - a1->pz;
                vx -= (a1->dpx*xt + a1->dpy*yt
                        + a1->dpz*zt) ;
                vx -= ( (xt*(.5*a1->qxx*xt + a1->qxy*yt + a1->qxz*zt)
                + yt*(.5*a1->qyy*yt + a1->qyz*zt) + .5*zt*a1->qzz*zt));
#ifdef CUBIC
		xt2 = xt*xt; yt2 = yt*yt;  zt2 = zt*zt;
		xt3 = xt2*xt; yt3 = yt2*yt; zt3 = zt2*zt;

	vx -= a1->qxxx*xt3/6. + a1->qxxy*xt2*yt/2 + a1->qxxz*xt2*zt/2
		+ a1->qxyy*xt*yt2/2 + a1->qxyz*xt*yt*zt + a1->qxzz*xt*zt2/2
		+ a1->qyyy*yt3/6 + a1->qyyz*yt2*zt/2 + a1->qyzz*yt*zt2/2 
		+ a1->qzzz*zt3/6.;
#endif
#ifdef QUARTIC
		xt4 = xt3*xt; yt4 = yt3*yt; zt4 = zt3*zt;
	vx -= a1->qxxxx*xt4/24. + a1->qxxxy*xt3*yt/6. + a1->qxxxz*xt3*yt/6. + a1->qxxyy*xt2*yt2/4.
		+ a1->qxxyz*xt2*yt*zt/2. + a1->qxxzz*xt2*zt2/4. + a1->qxyyy*xt*yt3/6.
		+ a1->qxyyz*xt*yt2*zt/2. + a1->qxyzz*xt*yt*zt2/2. + a1->qxzzz*xt*zt3/6.
		+ a1->qyyyy*yt4/24. + a1->qyyyz*yt3*zt/6. + a1->qyyzz*yt2*zt2/4. + a1->qyzzz*yt*zt3/6.
		+ a1->qzzzz*zt4/24.;
#endif
#ifdef QUINTIC
		xt5 = xt4*xt; yt5 = yt4*yt; zt5 = zt4*zt;
	vx -= a1->qxxxxx*xt5/120. + a1->qxxxxy*xt4*yt/24. + a1->qxxxxz*xt4*zt/24.
		+ a1->qxxxyy*xt3*yt2/12. + a1->qxxxyz*xt3*yt*zt/6. + a1->qxxxzz*xt3*zt2/12.
		+ a1->qxxyyy*xt2*yt3/12. + a1->qxxyyz*xt2*yt2*zt/4. + a1->qxxyzz*xt2*yt*zt2/4.
		+ a1->qxxzzz*xt2*zt3/12. + a1->qxyyyy*xt*yt4/24.  + a1->qxyyyz*xt*yt3*zt/6.
		+ a1->qxyyzz*xt*yt2*zt2/4. + a1->qxyzzz*xt*yt*zt3/6. + a1->qxzzzz*xt*zt4/24.
		+ a1->qyyyyy*yt5/120. + a1->qyyyyz*yt4*zt/24 + a1->qyyyzz*yt3*zt2/12.
		+ a1->qyyzzz*yt2*zt3/12. + a1->qyzzzz*yt*zt4/24. + a1->qzzzzz*zt5/120.;

#endif
/* do the close atoms */
	for( jj=0; jj< NCLOSE; jj++)
	{ if( a1->close[jj] == NULL) break; }
        for( ii=0; ii< jj;ii++)
        {
        a2 = a1->close[ii]; 
        xt = (a2->dx -a1->dx)*lambda +(a2->x -a1->x);
        yt = (a2->dy -a1->dy)*lambda +(a2->y -a1->y);
        zt = (a2->dz -a1->dz)*lambda +(a2->z -a1->z);
        r = xt*xt + yt*yt + zt*zt; r0 = sqrt(r);
/*      xt = xt/r0; yt = yt/r0; zt = zt/r0;
*/
/* standard point atom
        k = dielectric*a1->q*a2->q/r0;
*/
/* 1s 1s  polared atoms */
	k = dielectric*a1->q*a2->q/r0;
	
	r4 = r*r;
	k -= 0.5*(a1->a*a2->q*a2->q + a2->a*a1->q*a1->q)/r4;
        r = r*r*r;
        k -= a1->a*a2->a/r;
        k += a1->b*a2->b/r/r;
	vx += k;
        }
        *V += vx;
	}
	a_inactive_f_zero();
        free( atomall); return 1;
}

/*zone_polar()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int zone_polar( V, lambda ,alist, inalist )
	float *V,lambda;
	ATOM *( *alist)[] ;
	int inalist;
{
	float r,r0,xt,yt,zt;
	float lcutoff,cutoff,get_f_variable();
	int inbond,inangle,i,ii;
	ATOM *a1,*a2;
	ATOM *a_next( ); /* returns first ATOM when called with -1 */
	float dielectric,ve,va,vh;
	ATOM *a_m_serial();

/* nonbonded potentials 
* do a double loop starting from the first atom to the 
* last 
* then from the second to the last 
* etc
*
* also check to avoid bonded and 1-3 bonded atoms
*/
	if( inalist <= 0 ) return 1;
	dielectric = get_f_variable("dielec");
	if( dielectric < 1.) dielectric = 1.;
	dielectric = 332.17752/dielectric;
	cutoff = get_f_variable("cutoff");
	if( cutoff < 1.) cutoff = 1.e10;
	lcutoff = -cutoff;
	for( ii=0; ii< inalist; ii++)
	{
	a1 = (*alist)[ii];
	if( a1 == NULL ) goto NOTANATOM;
	ve = 0.; va = 0.; vh = 0.;
	a2 = a_next(-1);
/*	
*	for(i = 0; i< a1->dontuse; i++)
*	printf("%d ",a1->excluded[i]->serial);
*	printf("\n");
*/
/*
	while(  (a2->next != a2) && (a2->next != NULL))
	*/
	while(  (a2 != NULL) && (a2->next != NULL) && a2->next != a2)
	{
/* goto SKIP is used because this is one case where it makes sense */
/*	if( a2 == a1) break;  */
	if( a2 == a1) goto SKIP;
	for(i = 0; i< a1->dontuse; i++)
		if( a2 == a1->excluded[i]) goto SKIP;
/* non - bonded are only used when the atoms arent bonded */

	xt = (a1->x -a2->x) + lambda*(a1->dx - a2->dx);
	if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
	yt = (a1->y -a2->y) + lambda*(a1->dy - a2->dy);
	if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
	zt = (a1->z -a2->z) + lambda*(a1->dz - a2->dz);
	if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;

	r = xt*xt+yt*yt+zt*zt;
	if( r < 1.) r = 1.;  

	r0 = sqrt(r); r = r*r*r ;
/* debugging
	printf(" %d %d %f %f %f \n", a1->serial,a2->serial,a1->q,a2->q,
	332.17752*a1->q*a2->q/r0);
*/
	ve += dielectric*a1->q*a2->q/r0; 
	va -= a1->a*a2->a/r;
	vh += a1->b*a2->b/r/r; 

SKIP:
/*	if( a2->next == a1) break; */
	if( a2->next == a2) break;  
	a2 = a2->next;
	}
	*V += ve + va + vh;
NOTANATOM:
	i = i;
	}
	return 1;

}


/* a_polar()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_polar( V, lambda,ilow,ihigh,op )
	float *V,lambda;
	int ilow,ihigh;
	FILE *op;
{
	float r,r0,r4,xt,yt,zt;
	float lcutoff,cutoff,get_f_variable();
	int inbond,inangle,i,ii;
	ATOM *a1,*a2;
	ATOM *a_next( ); /* returns first ATOM when called with -1 */
	float dielectric,ve,va,vh,vp;
	float vel,val,vhl,vpl;
	float vtint,vtout, vtt;
	ATOM *a_m_serial();
	float k;


/* nonbonded potentials 
* do a double loop starting from the first atom to the 
* last 
* then from the second to the last 
* etc
*
* also check to avoid bonded and 1-3 bonded atoms
*/
	dielectric = get_f_variable("dielec");
	if( dielectric < 1.) dielectric = 1.;
	dielectric = 332.17752/dielectric;
	cutoff = get_f_variable("cutoff");
	if( cutoff < 1.) cutoff = 1.e10;
	lcutoff = -cutoff;
	vtint = 0.; vtout = 0.; vtt = 0.;
	for( ii=ilow; ii<=ihigh; ii++)
	{
	a1 = a_m_serial(ii);
	if( a1 == NULL ) goto NOTANATOM;
	ve = 0.; va = 0.; vh = 0.; vp = 0.;
	vel = 0.; val = 0.; vhl = 0.;vpl = 0.;
	a2 = a_next(-1);
/*	
*	for(i = 0; i< a1->dontuse; i++)
*	printf("%d ",a1->excluded[i]->serial);
*	printf("\n");
*/
/*
	while(  (a2->next != a2) && (a2->next != NULL))
	*/
	while(   (a2->next != NULL))
	{
/* goto SKIP is used because this is one case where it makes sense */
/*	if( a2 == a1) break;  */
/*	if( a2 == a1) goto SKIP;  */
	for(i = 0; i< a1->dontuse; i++)
		if( a2 == a1->excluded[i]) goto SKIP;
/* non - bonded are only used when the atoms arent bonded */

	if( lambda == 0.)
	{
	xt = (a1->x - a2->x);
	if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
	yt =  (a1->y - a2->y);
	if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
	zt =  (a1->z - a2->z);
	if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;
	} else
	{
	xt = (a1->x - a2->x +lambda*(a1->dx - a2->dx));
	if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
	yt = (a1->y - a2->y +lambda*(a1->dy - a2->dy));
	if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
	zt = (a1->z - a2->z +lambda*(a1->dz - a2->dz));
	if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;
	}
	r = xt*xt+yt*yt+zt*zt;
	if( r < .01) r = .01;  

	r0 = sqrt(r); r4 = r*r; r = r*r*r ;
/* debugging
*	printf(" %d %d %f %f %f \n", a1->serial,a2->serial,a1->q,a2->q,
*	332.17752*a1->q*a2->q/r0);
*/
/*
*	ve += dielectric*a1->q*a2->q/r0; 
*/
/* 1s 1s  polared atoms */
	k = dielectric*a1->q*a2->q/r0;
/*
*/
	ve += k;
	va -= a1->a*a2->a/r;
	vh += a1->b*a2->b/r/r; 
	vp -= 0.5*( a1->q*a1->q*a2->a + a1->a*a2->q*a2->q)/r4;
	if( a2->serial < ilow || a2->serial > ihigh)
	{
/*
	vel += dielectric*a1->q*a2->q/r0; 
*/
	vel += k;
	val -= a1->a*a2->a/r;
	vhl += a1->b*a2->b/r/r; 
	vpl -= 0.5*( a1->q*a1->q*a2->a + a1->a*a2->q*a2->q)/r4;
	}

SKIP:
/*	if( a2->next == a1) break; */
	if( a2->next == a2) break;  
	a2 = a2->next;
	}
	fprintf(op,"Vpolar internal %d Eq %f Ep %f E6 %f E12 %f\n",
	ii,ve-vel,vp-vpl,va-val,vh-vhl);
	fprintf(op,"Vpolar external %d Eq %f Ep  %f E6 %f E12 %f\n",ii,vel,vpl,val,vhl);
	fprintf(op,"Vpolar total    %d Eq %f Ep %f E6 %f E12 %f\n",ii,ve,vp,va,vh);
	*V += ve + va + vh+ vp;
	vtint += ve -vel+ va -val + vh -vhl + vp - vpl;
	vtout += vel + val + vhl + vpl;
	vtt  += ve + va + vh +vp;
NOTANATOM:
	i = i;
	}
	fprintf(op," Vpolar total internal %f \n",vtint);
	fprintf(op," Vpolar total external %f \n",vtout);
	fprintf(op," Vpolar total          %f \n",vtt);
	return 1;

}


int fv_update_polar(  lambda )
	float lambda;
{
	float r,r4,r0,xt,yt,zt,pt;
	float k,k1,k2,k3,k4,k5;
	float ka3,ka2;
	float kb3,kb2;
	float kp3,kp2;
	float get_f_variable();
	int inbond,inangle,i;
	ATOM *a1,*a2,*bonded[10],*angled[10];
	ATOM *a_next( ); /* returns first ATOM when called with -1 */
	int (*indexes)[],inindex,in;
	int a_number();
	int ii,j,jj,imax,inclose;
	float (*vector)[];
/*	float (*vecold)[];
*/
	ATOM *close[NCLOSE],*(*atomall)[];
	float mxdq,dielectric,mxcut; 
	static float dielecold = -1.;


	mxdq = get_f_variable("mxdq");
/*	if( mxdq < 0.) mxdq = 0.;
*/
	mxcut = get_f_variable("mxcut");
	if( mxcut < 0.) mxcut= 5.;
	if( mxdq > 0.) mxdq = mxdq*mxdq;

	dielectric = get_f_variable("dielec");
	if( dielectric < 1.) dielectric = 1.;
	if( dielecold != dielectric)
	{
	 dielecold = dielectric;
	mxdq = -1.;
	}
	dielectric = 332.17752/dielectric;
	
/*  get the number of atoms and allocate the memory for the array space */
	i = a_number();
	vector = malloc( 4*i*sizeof(float) );
	if( vector == NULL) 
	{ aaerror("cannot allocate memory in v_nonbon\n"); return 0;}
/*
	vecold = malloc( 4*i*sizeof(float) );
	if( vecold == NULL) 
	{ aaerror("cannot allocate memory in v_nonbon\n"); return 0;}
*/
	atomall = malloc( i*sizeof(ATOM *) );
	if( atomall == NULL)
	{aaerror("cannot allocate memory in v_nonbon\n"); return 0;}

	imax = a_number();
	for( i=0; i< imax; i++)
	{
		(*atomall)[i] = a_next(i);
	}
/* first check if anyone's moved and update the lists */
/* note that this must be a look-ahead rather than
*  look back search because
* we cannot update ->px until we've used that atom !!! */
	for( ii=0; ii< imax; ii++)
	{
        a1 = (*atomall)[ii];
        xt = a1->dx*lambda +a1->x - a1->px;
        yt = a1->dy*lambda +a1->y - a1->py;
        zt = a1->dz*lambda +a1->z - a1->pz;
        r = xt*xt + yt*yt + zt*zt;
	if( r > mxdq ) goto DOIT;
	}
	
	free( vector);
/*	free( vecold);
*/
	free (atomall);
	return 1;
DOIT:
	xt = get_f_variable("mmbox");
	if( xt > 0.)
	{ free(vector); free(atomall); mm_fv_update_polar(lambda); return 1;}
	indexes = malloc( imax* sizeof(int) );
	if( indexes == NULL ){
	aaerror(" cannot allocate memory in fv_update\n");
	return 0;}
	for( ii=0; ii< imax; ii++)
	{
	a1 = (*atomall)[ii];
	a1 -> VP = 0.;
	a1 -> dpx = 0.;
	a1 -> dpy = 0.;
	a1 -> dpz = 0.;
	a1 -> qxx = 0.;
	a1 -> qxy = 0.;
	a1 -> qxz = 0.;
	a1 -> qyy = 0.;
	a1 -> qyz = 0.;
	a1 -> qzz = 0.;
#ifdef CUBIC
	a1 -> qxxx = 0.;
	a1 -> qxxy = 0.;
	a1 -> qxxz = 0.;
	a1 -> qxyy = 0.;
	a1 -> qxyz = 0.;
	a1 -> qxzz = 0.;
	a1 -> qyyy = 0.;
	a1 -> qyyz = 0.;
	a1 -> qyzz = 0.;
	a1 -> qzzz = 0.;
#endif
	for( j=0; j< NCLOSE; j++) 
		a1->close[j] = NULL;

	}
	for( ii=0; ii<  imax; ii++)
	{ /* if this is met we update the expansion for this atom */
	a1 = (*atomall)[ii];
	inclose = 0;
	if( lambda != 0.)
	{
#pragma _CNX no_recurrence
	for( i=ii+1; i< imax; i++)
	{
	a2 = (*atomall)[i];
	j = i*4;
	(*vector)[j  ] = a2->x - a1->x + lambda*(a2->dx -a1->dx);
	(*vector)[j+1] = a2->y - a1->y + lambda*(a2->dy -a1->dy);
	(*vector)[j+2] = a2->z - a1->z + lambda*(a2->dz -a1->dz);
	}
	}else {
#pragma _CNX no_recurrence
	for( i=ii+1; i< imax; i++)
	{
	a2 = (*atomall)[i];
	j = i*4;
	(*vector)[j  ] = a2->x - a1->x ;
	(*vector)[j+1] = a2->y - a1->y ;
	(*vector)[j+2] = a2->z - a1->z ;
	}
	} /* end of difference position into vector loops */
#pragma _CNX no_recurrence
	for( i=ii+1; i< imax; i++)
	{
		j = i*4;
		(*vector)[j+3] = sqrt((*vector)[j]*(*vector)[j] +
				 (*vector)[j+1]*(*vector)[j+1] +
				 (*vector)[j+2]*(*vector)[j+2]);
	}
/* add the new components */
/* first extract indexes */
	inindex = 0;
	for( i=ii+1; i< imax; i++)
	{
	a2 = (*atomall)[i];
	for( j=0; j< a1->dontuse; j++)
	{ if( a2 == a1->excluded[j]) goto SKIPNEW;}
	j = i*4;
	if( (*vector)[j+3] > mxcut || inclose > NCLOSE)
	{
	(*indexes)[inindex++] = i;
	}else {
	a1->close[inclose++] = (*atomall)[i];
	}
	if( inclose == NCLOSE)
	{
	aaerror(
	" fv_update_nonbon> too many atoms increase NCLOSE or decrease mxcut");
	}
SKIPNEW:   i = i;
	}
	/* and then use them */
#pragma _CNX no_recurrence
	for( in=0; in< inindex; in++)
	{
	i = (*indexes)[in];
	a2 = (*atomall)[i];
	j = i*4;
	r0 = one/(*vector)[j+3];  
	r = r0*r0;
	r4 = r*r;
	r = r*r*r; /* r0^-6 */
	xt = a1->q*a2->q*dielectric*r0;
	pt = 0.5*(a1->a*a2->q*a2->q + a1->q*a1->q*a2->q)*r4;
	yt = a1->a*a2->a*r;
	zt = a1->b*a2->b*r*r;
	k = xt - yt + zt -pt;
	xt = xt*r0; yt = yt*r0; zt = zt*r0; pt = pt*r0;
	k1 = xt - yt*six + zt*twelve -pt*4;
	xt = xt*r0; yt = yt*r0; zt = zt*r0; pt = pt*r0;
	/*
	k2 = xt*three; ka2 = - yt*6*8; kb2 =  zt*12*14;  kp2 = -pt*4*6;
	*/
	k2 = xt*three; ka2 = - yt*48.; kb2 =  zt*168.; kp2 = -pt*24;
#ifdef CUBIC
	xt = xt*r0; yt = yt*r0; zt = zt*r0;
	k3 = -xt*5*3; ka3 =   yt*6*8*10 ; kb3 =  -zt*12*14*16;  kp3 = pt*4*6*8;
#endif
	k1 = -k1;
	xt = (*vector)[j]*r0 ;
	yt = (*vector)[j+1]*r0 ;
	zt = (*vector)[j+2] *r0;
	/*
	xt = (*vector)[j] ;
	yt = (*vector)[j+1] ;
	zt = (*vector)[j+2] ;
	*/
	a1->VP += k;
	a2->dpx -= k1*xt;
	a1->dpx += k1*xt;
	a2->dpy -= k1*yt;
	a1->dpy += k1*yt;
	a2->dpz -= k1*zt;
	a1->dpz += k1*zt;
/*  note that xt has the 1/r in it so k2*xt*xt is 1/r^5 */
	a2->qxx -= k2*(xt*xt - third) + ka2*(xt*xt - eightth) + kb2*(xt*xt-fourteenth) +kp2*(xt*xt - 0.166666666666666666666);
	a1->qxx -= k2*(xt*xt - third) + ka2*(xt*xt - eightth) + kb2*(xt*xt-fourteenth) +kp2*(xt*xt -0.1666666666666666666666);
	a2->qxy -= (k2+ka2+kb2 +kp2)*yt*xt;
	a1->qxy -= (k2+ka2+kb2 +kp2)*yt*xt;
	a2->qxz -= (k2+ka2+kb2 +kp2)*zt*xt;
	a1->qxz -= (k2+ka2+kb2 +kp2)*zt*xt;
	a2->qyy -= k2*(yt*yt - third) + ka2*(yt*yt - eightth) + kb2*(yt*yt-fourteenth) +kp2*(yt*yt-0.166666666666666666666);
	a1->qyy -= k2*(yt*yt - third) + ka2*(yt*yt - eightth) + kb2*(yt*yt-fourteenth) +kp2*(yt*yt-0.16666666666666666666666);
	a2->qyz -= (k2+ka2+kb2 +kp2)*yt*zt;
	a1->qyz -= (k2+ka2+kb2 +kp2)*yt*zt;
	a2->qzz -= k2*(zt*zt - third) + ka2*(zt*zt - eightth) + kb2*(zt*zt-fourteenth) +kp2*(zt*zt-0.1666666666666666666666);
	a1->qzz -= k2*(zt*zt - third) + ka2*(zt*zt - eightth) + kb2*(zt*zt-fourteenth) +kp2*(zt*zt-0.16666666666666666666666);
#ifdef CUBIC
	a2->qxxx -= k3*(xt*xt*xt - xt*( 9./15 )) ;
	a2->qxxx -= ka3*(xt*xt*xt - xt*( 24./80 )) ;
	a2->qxxx -= kb3*(xt*xt*xt - xt*( 42./(14*16))); 
	a1->qxxx += k3*(xt*xt*xt - xt*( 9./15 )) ;
	a1->qxxx += ka3*(xt*xt*xt - xt*( 24./80 )) ;
	a1->qxxx += kb3*(xt*xt*xt - xt*( 42./(14*16))); 
	a2->qxxy -= k3*(yt*xt*xt - yt*( 6./ 15));
	a2->qxxy -= ka3*(yt*xt*xt - yt*( 11./ 80));
	a2->qxxy -= kb3*(yt*xt*xt - yt*( 17./ (14*16)));
	a1->qxxy += k3*(yt*xt*xt - yt*( 6./ 15));
	a1->qxxy += ka3*(yt*xt*xt - yt*( 11./ 80));
	a1->qxxy += kb3*(yt*xt*xt - yt*( 17./ (14*16)));
	a2->qxxz -= k3*(zt*xt*xt - zt*( 6./ 15));
	a2->qxxz -= ka3*(zt*xt*xt - zt*( 11./ 80));
	a2->qxxz -= kb3*(zt*xt*xt - zt*( 17./ (14*16)));
	a1->qxxz += k3*(zt*xt*xt - zt*( 6./ 15));
	a1->qxxz += ka3*(zt*xt*xt - zt*( 11./ 80));
	a1->qxxz += kb3*(zt*xt*xt - zt*( 17./ (14*16)));
	a2->qxyy -= k3*(yt*yt*xt - xt*( 6./ 15));
	a2->qxyy -= ka3*(yt*yt*xt - xt*( 11./ 80));
	a2->qxyy -= kb3*(yt*yt*xt - xt*( 17./ (14*16)));
	a1->qxyy += k3*(yt*yt*xt - xt*( 6./ 15));
	a1->qxyy += ka3*(yt*yt*xt - xt*( 11./ 80));
	a1->qxyy += kb3*(yt*yt*xt - xt*( 17./ (14*16)));
	a2->qxyz -= (k3+ka3+kb3)*yt*zt*xt;
	a1->qxyz += (k3+ka3+kb3)*yt*zt*xt;
	a2->qxzz -= k3*(zt*zt*xt - xt*( 6./ 15));
	a2->qxzz -= ka3*(zt*zt*xt - xt*( 11./ 80));
	a2->qxzz -= kb3*(zt*zt*xt - xt*( 17./ (14*16)));
	a1->qxzz += k3*(zt*zt*xt - xt*( 6./ 15));
	a1->qxzz += ka3*(zt*zt*xt - xt*( 11./ 80));
	a1->qxzz += kb3*(zt*zt*xt - xt*( 17./ (14*16)));
	a2->qyyy -= k3*(yt*yt*yt - yt*( 9./15 )) ;
	a2->qyyy -= ka3*(yt*yt*yt - yt*( 24./80 )) ;
	a2->qyyy -= kb3*(yt*yt*yt - yt*( 42./(14*16))); 
	a1->qyyy += k3*(yt*yt*yt - yt*( 9./15 )) ;
	a1->qyyy += ka3*(yt*yt*yt - yt*( 24./80 )) ;
	a1->qyyy += kb3*(yt*yt*yt - yt*( 42./(14*16))); 
	a2->qyyz -= k3*(yt*yt*zt - zt*( 6./ 15));
	a2->qyyz -= ka3*(yt*yt*zt - zt*( 11./ 80));
	a2->qyyz -= kb3*(yt*yt*zt - zt*( 17./ (14*16)));
	a1->qyyz += k3*(yt*yt*zt - zt*( 6./ 15));
	a1->qyyz += ka3*(yt*yt*zt - zt*( 11./ 80));
	a1->qyyz += kb3*(yt*yt*zt - zt*( 17./ (14*16)));
	a2->qyzz -= k3*(zt*zt*yt - yt*( 6./ 15));
	a2->qyzz -= ka3*(zt*zt*yt - yt*( 11./ 80));
	a2->qyzz -= kb3*(zt*zt*yt - yt*( 17./ (14*16)));
	a1->qyzz += k3*(zt*zt*yt - yt*( 6./ 15));
	a1->qyzz += ka3*(zt*zt*yt - yt*( 11./ 80));
	a1->qyzz += kb3*(zt*zt*yt - yt*( 17./ (14*16)));
	a2->qzzz -= k3*(zt*zt*zt - zt*( 9./15 )) ;
	a2->qzzz -= ka3*(zt*zt*zt - zt*( 24./80 )) ;
	a2->qzzz -= kb3*(zt*zt*zt - zt*( 42./(14*16))); 
	a1->qzzz += k3*(zt*zt*zt - zt*( 9./15 )) ;
	a1->qzzz += ka3*(zt*zt*zt - zt*( 24./80 )) ;
	a1->qzzz += kb3*(zt*zt*zt - zt*( 42./(14*16))); 
#endif

/* debugging
	j = i *4;
	fprintf(stderr," mxcut %f %f inclose %d who %d \n",mxcut,(*vector)[j+3],inclose,(*atomall)[i]->serial);
	fprintf(stderr," vector %f %f %f \n", (*vector)[j],(*vector)[j+1],(*vector)[j+2]);
*/
	}/* end of loop i */
/* merge the non-bond mxcut lists */
/*	a1->close[inclose] = NULL; */
/* set the position */
	a1->px = a1->dx*lambda + a1->x;
	a1->py = a1->dy*lambda + a1->y;
	a1->pz = a1->dz*lambda + a1->z;

	}  /* end of ii loop */
	
	a_inactive_f_zero();

	free( indexes);
	free( vector);
/*	free( vecold);
*/
	free (atomall);
	return 1;

}



/*
 * *
 * polarmm.c   multipole method with polarization
 * *
 * v 0.5 add the polarization to the neighbor terms (closest boxes
 *   since it is 1/r^3 this should converge.  
 *   if this doesn't work then we'll add it to the other boxes.
 *    
*/
/* rectmm.c
*
* integrated multipole method with the amortized
* standard nonbonded program.
*  rectangular multipole expansion to the r^-6 order (5th order expansion)
*
*  from Eyges "The Classical Electromagnetic Field"
*  note that we use the oposite convention for sign of
*  expansion so use + for all the cumulants, while he uses
*   -1^n.   This is solely due to choice of origin and for
*   other applications (rxn feild ) the -1^n is correct.
*
*
* collection of routines to service nonbonded potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are symmetric - so we don't have to fuck around with
* s matrices and the like.
*
* note that the non-bonded information is in the ATOM structures 
*
*
* attempts at vectorization
*/
/*
*  copyright 1992, 1993, 1994, 1995 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/

/*
#define FOURTH
#define  FIFTH
#ifdef FIFTH
#define FOURTH
#endif
*/
typedef struct {
  float xc,yc,zc;
  float sqp;
  /* dipole appx for r^-6 */
  float sa,xa,ya,za;
  float q100,q010,q001;
  float q200,q020,q002,q110,q101,q011;
  float q300,q030,q003,q210,q201,q120,q021,q102,q012,q111;
#ifdef FOURTH
  float q400,q040,q004,q310,q301,q130,q031,q103,q013,q220,q202,q022,q211,q121,q112;
#endif
#ifdef FIFTH
  float q500,q050,q005,q410,q401,q140,q041,q104,q014,q320,q230,q302,q203,q032,q023,q311,q131,q113,q221,q212,q122;
#endif
  int first,last,innode; } MMNODE;
typedef struct {
 ATOM *who;
 int next,which; } MMATOM;
int mm_fv_update_polar(  lambda )
	float lambda;
{
	float r,r0,xt,yt,zt;
/* polarization terms */
	float r4, pt,kp2,kp3;
	float xt2,xt3,xt4,yt2,yt3,yt4,zt2,zt3,zt4;
	float k,k1,k2,k3,k4,k5;
	float ka2,ka3;
	float kb2,kb3;
	float c1,c2,c3,c4,c5; /* constants for the mm expansion */
	float get_f_variable();
	int inbond,inangle,i;
	ATOM *ap,*a1,*a2,*bonded[10],*angled[10];
	ATOM *a_next( ); /* returns first ATOM when called with -1 */
	int a_number();
	int ii,j,jj,imax,inclose;
	float (*vector)[];
/* */
	ATOM *close[NCLOSE],*(*atomall)[];
	float mxdq,dielectric,mxcut; 
	float mmbox;
	float xmax,xmin,ymax,ymin,zmax,zmin;
	int nx,ny,nz;
	int ix,iy,iz,inode;
	int naybor[27];
	char line[80];
	MMNODE (*nodelist)[];
	MMATOM (*atomlist)[];


	mmbox = get_f_variable("mmbox");
	mxcut = get_f_variable("mxcut");
	if( mxcut < 0.) mxcut= 5.;

	dielectric = get_f_variable("dielec");
	if( dielectric <= 0.) dielectric = 1.;
	dielectric = 332.17752/dielectric;
	
/*  get the number of atoms and allocate the memory for the array space */
	i = a_number();
	vector = malloc( 4*i*sizeof(float) );
	if( vector == NULL) 
	{ aaerror("cannot allocate memory in mm_fv_update\n"); return 0;}
	atomall = malloc( i*sizeof(ATOM *) );
	if( atomall == NULL)
	{aaerror("cannot allocate memory in mm_fv_update\n"); return 0;}
	atomlist = malloc( i * sizeof( MMATOM ));
	if( atomlist == NULL)
	{ aaerror("cannot allocate memory in mm_fv_update\n"); return 0;}

	imax = a_number();
	jj = imax;

	for( i=0; i< imax; i++)
	{
		(*atomall)[i] = a_next(i);
		(*atomlist)[i].who = (*atomall)[i];
	}
/* first check if anyone's moved and update the lists */
/* note that this must be a look-ahead rather than
*  look back search because
* we cannot update ->px until we've used that atom !!! */
#pragma _CNX no_recurrence
	for( ii=0; ii< imax; ii++)
	{
        a1 = (*atomall)[ii];
	j = ii*4;
        (*vector)[j] = a1->dx*lambda +a1->x ;
        (*vector)[j+1] = a1->dy*lambda +a1->y;
        (*vector)[j+2] = a1->dz*lambda +a1->z;
	}
/* determine the bounds of box which surrounds all of the atoms */
	xmax = -10e10;
	ymax = -10e10;
	zmax = -10e10;
	xmin =  10e10;
	ymin =  10e10;
	zmin =  10e10;
	for( ii= 0; ii< imax; ii++)
	{
		j = ii*4;
	  if( xmax < (*vector)[j] ) xmax = (*vector)[j];
	  if( ymax < (*vector)[j+1] ) ymax = (*vector)[j+1];
	  if( zmax < (*vector)[j+2] ) zmax = (*vector)[j+2];
	  if( xmin > (*vector)[j] ) xmin = (*vector)[j];
	  if( ymin > (*vector)[j+1] ) ymin = (*vector)[j+1];
	  if( zmin > (*vector)[j+2] ) zmin = (*vector)[j+2];
	  }
	nx = (xmax - xmin)/mmbox + 1 ;
	ny = (ymax - ymin)/mmbox + 1 ;
	nz = (zmax - zmin)/mmbox + 1 ;

/* DEBUG	
	sprintf(line,"before allocation mmbox %f nx %d ny %d nz %d \n",mmbox,nx,ny,nz);
	aaerror( line);
	sprintf(line," xmin xmax %f %f ymin ymax %f %f zmin zmax %f %f\n",
	xmin,xmax,ymin,ymax,zmin,zmax);
	aaerror( line);
end of DEBUG */

/* now try to malloc the mmnodes */
	nodelist = malloc( nx*ny*nz * sizeof( MMNODE ));
	if( nodelist == NULL)
	{ aaerror("cannot allocate node memory in mm_fv_update (doubling grid )\n"); 
	sprintf(line,"mmbox %f nx %d ny %d nz %d ",mmbox,nx,ny,nz);
	aaerror( line);
	sprintf(line," xmin xmax %f %f ymin ymax %f %f zmin zmax %f %f",
	xmin,xmax,ymin,ymax,zmin,zmax);
	aaerror( line);
	mmbox = mmbox *2;
	set_f_variable( "mmbox",mmbox);
	nx = (xmax - xmin)/mmbox + 1;
	ny = (ymax - ymin)/mmbox + 1;
	nz = (zmax - zmin)/mmbox + 1;
	nodelist = malloc( nx*ny*nz * sizeof( MMNODE ));
	if( nodelist == NULL)
	{ aaerror("cannot allocate node memory in mm_fv_update (cannot do it)\n"); return 0; }
	}
	for( ix=0; ix< nx; ix++)
	for( iy=0; iy< ny; iy++)
	for( iz=0; iz< nz; iz++)
	{
	inode = ((iz*ny)+iy)*nx + ix;
	(*nodelist)[inode].xc = ix*mmbox + .5*mmbox + xmin;
	(*nodelist)[inode].yc = iy*mmbox + .5*mmbox + ymin;
	(*nodelist)[inode].zc = iz*mmbox + .5*mmbox + zmin;
	}
#pragma _CNX no_recurrence
	for( ii=0; ii < nx*ny*nz; ii++)
	{
	(*nodelist)[ii].sqp = 0.;
	(*nodelist)[ii].sa = 0.;
	(*nodelist)[ii].xa = 0.;
	(*nodelist)[ii].ya = 0.;
	(*nodelist)[ii].za = 0.;
	(*nodelist)[ii].q100 = 0.;
	(*nodelist)[ii].q010 = 0.;
	(*nodelist)[ii].q001 = 0.;
	(*nodelist)[ii].q200 = 0.;
	(*nodelist)[ii].q020 = 0.;
	(*nodelist)[ii].q002 = 0.;
	(*nodelist)[ii].q101 = 0.;
	(*nodelist)[ii].q110 = 0.;
	(*nodelist)[ii].q011 = 0.;
	(*nodelist)[ii].q300 = 0.;
	(*nodelist)[ii].q030 = 0.;
	(*nodelist)[ii].q003 = 0.;
	(*nodelist)[ii].q210 = 0.;
	(*nodelist)[ii].q120 = 0.;
	(*nodelist)[ii].q201 = 0.;
	(*nodelist)[ii].q102 = 0.;
	(*nodelist)[ii].q021 = 0.;
	(*nodelist)[ii].q012 = 0.;
	(*nodelist)[ii].q111 = 0.;
#ifdef FOURTH
	(*nodelist)[ii].q400 = 0.;
	(*nodelist)[ii].q040 = 0.;
	(*nodelist)[ii].q004 = 0.;
	(*nodelist)[ii].q310 = 0.;
	(*nodelist)[ii].q130 = 0.;
	(*nodelist)[ii].q301 = 0.;
	(*nodelist)[ii].q103 = 0.;
	(*nodelist)[ii].q031 = 0.;
	(*nodelist)[ii].q013 = 0.;
	(*nodelist)[ii].q220 = 0.;
	(*nodelist)[ii].q202 = 0.;
	(*nodelist)[ii].q022 = 0.;
	(*nodelist)[ii].q211 = 0.;
	(*nodelist)[ii].q121 = 0.;
	(*nodelist)[ii].q112 = 0.;
#endif
#ifdef FIFTH
	(*nodelist)[ii].q500 = 0.;
	(*nodelist)[ii].q050 = 0.;
	(*nodelist)[ii].q005 = 0.;
	(*nodelist)[ii].q410 = 0.;
	(*nodelist)[ii].q140 = 0.;
	(*nodelist)[ii].q401 = 0.;
	(*nodelist)[ii].q104 = 0.;
	(*nodelist)[ii].q041 = 0.;
	(*nodelist)[ii].q014 = 0.;
	(*nodelist)[ii].q320 = 0.;
	(*nodelist)[ii].q230 = 0.;
	(*nodelist)[ii].q302 = 0.;
	(*nodelist)[ii].q203 = 0.;
	(*nodelist)[ii].q032 = 0.;
	(*nodelist)[ii].q023 = 0.;
	(*nodelist)[ii].q221 = 0.;
	(*nodelist)[ii].q212 = 0.;
	(*nodelist)[ii].q122 = 0.;
	(*nodelist)[ii].q311 = 0.;
	(*nodelist)[ii].q131 = 0.;
	(*nodelist)[ii].q113 = 0.;
#endif
	(*nodelist)[ii].first = -1;
	(*nodelist)[ii].last = -1;
	(*nodelist)[ii].innode = 0;
	}
/* now decide for each atom who he belongs to */	
#pragma _CNX no_recurrence
	for( ii=0; ii< imax; ii++)
	{
	j = ii*4;
	ix = ((*vector)[j] - xmin )/mmbox;
	iy = ((*vector)[j+1] - ymin )/mmbox;
	iz = ((*vector)[j+2] - zmin )/mmbox;
	inode = ((iz*ny)+iy)*nx + ix;
	(*atomlist)[ii].which = inode;
/* DEBUG
	printf(" error %f %f %f %d %d %d %d\n",
		(*vector)[j],(*vector)[j+1],(*vector)[j+2],ix,iy,iz,inode);

ENDDEBUG */
	}
/* and generate the links */ 
	for( inode = 0; inode < nx*ny*nz; inode++)
	{
	/* first find the first atom which belongs to me */
	for( ii = 0; ii< imax; ii++)
	{
		if( (*atomlist)[ii].which == inode)
		{
		(*nodelist)[inode].first = ii;
		(*nodelist)[inode].last = ii;
		(*nodelist)[inode].innode += 1;
		ap = (*atomlist)[ii].who;
		break;
		}
	}
	/* only if i'm not null */
	if( ii != imax )
	{
	for( ii= (*nodelist)[inode].first; ii < imax; ii++)
	{
		if( (*atomlist)[ii].which == inode)
		{
		(*atomlist)[(*nodelist)[inode].last].next  = ii;
		(*nodelist)[inode].last = ii;
		(*nodelist)[inode].innode += 1;
		ap = (*atomlist)[ii].who;
		xt = ap->x + lambda*ap->dx - (*nodelist)[inode].xc;
		yt = ap->y + lambda*ap->dy - (*nodelist)[inode].yc;
		zt = ap->z + lambda*ap->dz - (*nodelist)[inode].zc;
		(*nodelist)[inode].sqp +=  ap->q;
		(*nodelist)[inode].sa +=  ap->a;
		(*nodelist)[inode].xa +=  ap->a*xt;
		(*nodelist)[inode].ya +=  ap->a*yt;
		(*nodelist)[inode].za +=  ap->a*zt;
		xt2 = xt*xt;
		xt3 = xt2*xt;
		xt4 = xt3*xt;
		yt2 = yt*yt;
		yt3 = yt2*yt;
		yt4 = yt3*yt;
		zt2 = zt*zt;
		zt3 = zt2*zt;
		zt4 = zt3*zt;
	(*nodelist)[inode].q100 += ap->q*xt;
	(*nodelist)[inode].q010 += ap->q*yt;
	(*nodelist)[inode].q001 += ap->q*zt;
	(*nodelist)[inode].q200 += ap->q*xt2;
	(*nodelist)[inode].q020 += ap->q*yt2;
	(*nodelist)[inode].q002 += ap->q*zt2;
	(*nodelist)[inode].q101 += ap->q*xt*zt;
	(*nodelist)[inode].q110 += ap->q*xt*yt;
	(*nodelist)[inode].q011 += ap->q*yt*zt;
	(*nodelist)[inode].q300 += ap->q*xt3;
	(*nodelist)[inode].q030 += ap->q*yt3;
	(*nodelist)[inode].q003 += ap->q*zt3;
	(*nodelist)[inode].q210 += ap->q*xt2*yt;
	(*nodelist)[inode].q120 += ap->q*xt*yt2;
	(*nodelist)[inode].q201 += ap->q*xt2*zt;
	(*nodelist)[inode].q102 += ap->q*xt*zt2;
	(*nodelist)[inode].q021 += ap->q*yt2*zt;
	(*nodelist)[inode].q012 += ap->q*yt*zt2;
	(*nodelist)[inode].q111 += ap->q*xt*yt*zt;
#ifdef FOURTH
	(*nodelist)[inode].q400 += ap->q*xt4;
	(*nodelist)[inode].q040 += ap->q*yt4;
	(*nodelist)[inode].q004 += ap->q*zt4;
	(*nodelist)[inode].q310 += ap->q*xt3*yt;
	(*nodelist)[inode].q130 += ap->q*xt*yt3;
	(*nodelist)[inode].q301 += ap->q*xt3*zt;
	(*nodelist)[inode].q103 += ap->q*xt*zt3;
	(*nodelist)[inode].q031 += ap->q*yt3*zt;
	(*nodelist)[inode].q013 += ap->q*yt*zt3;
	(*nodelist)[inode].q220 += ap->q*xt2*yt2;
	(*nodelist)[inode].q202 += ap->q*xt2*zt2;
	(*nodelist)[inode].q022 += ap->q*yt2*zt2;
	(*nodelist)[inode].q211 += ap->q*xt2*yt*zt;
	(*nodelist)[inode].q121 += ap->q*xt*yt2*zt;
	(*nodelist)[inode].q112 += ap->q*xt*yt*zt2;
#endif
#ifdef FIFTH
	(*nodelist)[inode].q500 += ap->q*xt4*xt;
	(*nodelist)[inode].q050 += ap->q*yt4*yt;
	(*nodelist)[inode].q005 += ap->q*zt4*zt;
	(*nodelist)[inode].q410 += ap->q*xt4*yt;
	(*nodelist)[inode].q140 += ap->q*yt4*xt;
	(*nodelist)[inode].q401 += ap->q*xt4*zt;
	(*nodelist)[inode].q104 += ap->q*zt4*xt;
	(*nodelist)[inode].q041 += ap->q*yt4*zt;
	(*nodelist)[inode].q014 += ap->q*zt4*yt;
	(*nodelist)[inode].q320 += ap->q*xt3*yt2;
	(*nodelist)[inode].q230 += ap->q*yt3*xt2;
	(*nodelist)[inode].q302 += ap->q*xt3*zt2;
	(*nodelist)[inode].q203 += ap->q*zt3*xt2;
	(*nodelist)[inode].q032 += ap->q*yt3*zt2;
	(*nodelist)[inode].q023 += ap->q*zt3*yt2;
	(*nodelist)[inode].q221 += ap->q*xt2*yt2*zt;
	(*nodelist)[inode].q212 += ap->q*xt2*yt*zt2;
	(*nodelist)[inode].q122 += ap->q*xt*yt2*zt2;
	(*nodelist)[inode].q311 += ap->q*xt3*yt*zt;
	(*nodelist)[inode].q131 += ap->q*xt*yt3*zt;
	(*nodelist)[inode].q113 += ap->q*xt*yt*zt3;
#endif
		}
	}/* ii */
	}/* checking if ii != imax */
	}/* inode */
/* and now (almost done with the MM setup)
* normalize the accumulated nodal data */
#pragma _CNX no_recurrence
	/* multiplied by .5 to correct for double counting */
	k = dielectric *.5;
	xt = .5/3.;
	yt = xt/4.;
	zt = yt/5.;
	for( ii = 0; ii < nx*ny*nz; ii ++)
	{
	(*nodelist)[ii].sqp *= k;
	(*nodelist)[ii].q100 *= k;
	(*nodelist)[ii].q010 *= k;
	(*nodelist)[ii].q001 *= k;
	(*nodelist)[ii].q200 *= .5*k;
	(*nodelist)[ii].q020 *= .5*k;
	(*nodelist)[ii].q002 *= .5*k;
	(*nodelist)[ii].q101 *= k;
	(*nodelist)[ii].q110 *= k;
	(*nodelist)[ii].q011 *= k;
	(*nodelist)[ii].q300 *= xt*k;
	(*nodelist)[ii].q030 *= xt*k;
	(*nodelist)[ii].q003 *= xt*k;
	(*nodelist)[ii].q210 *= 0.5*k;
	(*nodelist)[ii].q120 *= 0.5*k;
	(*nodelist)[ii].q201 *= 0.5*k;
	(*nodelist)[ii].q102 *= 0.5*k;
	(*nodelist)[ii].q021 *= 0.5*k;
	(*nodelist)[ii].q012 *= 0.5*k;
	(*nodelist)[ii].q111 *= k;
#ifdef FOURTH
	(*nodelist)[ii].q400 *= yt*k;
	(*nodelist)[ii].q040 *= yt*k;
	(*nodelist)[ii].q004 *= yt*k;
	(*nodelist)[ii].q310 *= xt*k;
	(*nodelist)[ii].q130 *= xt*k;
	(*nodelist)[ii].q301 *= xt*k;
	(*nodelist)[ii].q103 *= xt*k;
	(*nodelist)[ii].q031 *= xt*k;
	(*nodelist)[ii].q013 *= xt*k;
	(*nodelist)[ii].q220 *= .25*k;
	(*nodelist)[ii].q202 *= .25*k;
	(*nodelist)[ii].q022 *= .25*k;
	(*nodelist)[ii].q211 *= .5*k;
	(*nodelist)[ii].q121 *= .5*k;
	(*nodelist)[ii].q112 *= .5*k;
#endif
#ifdef FIFTH
	(*nodelist)[ii].q500 *= zt*k;
	(*nodelist)[ii].q050 *= zt*k;
	(*nodelist)[ii].q005 *= zt*k;
	(*nodelist)[ii].q410 *= yt*k;
	(*nodelist)[ii].q140 *= yt*k;
	(*nodelist)[ii].q401 *= yt*k;
	(*nodelist)[ii].q104 *= yt*k;
	(*nodelist)[ii].q041 *= yt*k;
	(*nodelist)[ii].q014 *= yt*k;
	(*nodelist)[ii].q320 *= .5*xt*k;
	(*nodelist)[ii].q230 *= .5*xt*k;
	(*nodelist)[ii].q302 *= .5*xt*k;
	(*nodelist)[ii].q203 *= .5*xt*k;
	(*nodelist)[ii].q032 *= .5*xt*k;
	(*nodelist)[ii].q023 *= .5*xt*k;
	(*nodelist)[ii].q221 *= .25*k;
	(*nodelist)[ii].q212 *= .25*k;
	(*nodelist)[ii].q122 *= .25*k;
	(*nodelist)[ii].q311 *= xt*k;
	(*nodelist)[ii].q131 *= xt*k;
	(*nodelist)[ii].q113 *= xt*k;
#endif
	/*debug
	printf("%d %f %f\n",ii,(*nodelist)[ii].sqp,(*nodelist)[ii].q100);
	*/
	if( (*nodelist)[ii].sa != 0.)
	{
	(*nodelist)[ii].xa = (*nodelist)[ii].xa/(*nodelist)[ii].sa;
	(*nodelist)[ii].ya = (*nodelist)[ii].ya/(*nodelist)[ii].sa;
	(*nodelist)[ii].za = (*nodelist)[ii].za/(*nodelist)[ii].sa;
	}
	(*nodelist)[ii].xa += (*nodelist)[ii].xc;
	(*nodelist)[ii].ya += (*nodelist)[ii].yc;
	(*nodelist)[ii].za += (*nodelist)[ii].zc;
/* correct for double counting */
	(*nodelist)[ii].sa  *= .5;
	}

/* initiallization of the mmnodes is done !!! */

/*  initialize the data for every atom */
	for( ii=0; ii< jj; ii++)
	{
	a1 = (*atomall)[ii];
	a1-> px = a1->x + lambda*a1->dx;
	a1-> py = a1->y + lambda*a1->dy;
	a1-> pz = a1->z + lambda*a1->dz;
	a1 -> VP = 0.;
	a1 -> dpx = 0.;
	a1 -> dpy = 0.;
	a1 -> dpz = 0.;
	a1 -> qxx = 0.;
	a1 -> qxy = 0.;
	a1 -> qxz = 0.;
	a1 -> qyy = 0.;
	a1 -> qyz = 0.;
	a1 -> qzz = 0.;
#ifdef CUBIC
	a1 -> qxxx = 0.;
	a1 -> qxxy = 0.;
	a1 -> qxxz = 0.;
	a1 -> qxyy = 0.;
	a1 -> qxyz = 0.;
	a1 -> qxzz = 0.;
	a1 -> qyyy = 0.;
	a1 -> qyyz = 0.;
	a1 -> qyzz = 0.;
	a1 -> qzzz = 0.;
#endif
	for( j=0; j< NCLOSE; j++) 
		a1->close[j] = NULL;

	}/* end of initializations */


	for( ii=0; ii<  jj; ii++)
	{ /* if this is met we update the expansion for this atom */
/*	a1 = (*atomall)[ii];
	atomall will be reused in this loop so we refer to atomlist
	*/
	a1 = (*atomlist)[ii].who;
	inclose = 0;
/* loop over the nodes 
   if the node is mine or a neighbor then use an
   explicit summation
   otherwise use the MM node */
	ix = (a1->px  - xmin )/mmbox ;
	iy = (a1->py  - ymin )/mmbox ;
	iz = (a1->pz  - zmin )/mmbox ;
	naybor[0] = ((iz*ny)+iy)*nx + ix;
	naybor[1] = ((iz*ny)+iy)*nx + ix+1;
	naybor[2] = ((iz*ny)+iy)*nx + ix-1;
	naybor[3] = ((iz*ny)+iy)*nx+nx + ix;
	naybor[4] = ((iz*ny)+iy)*nx-nx + ix;
	naybor[5] = ((iz*ny)+iy)*nx+nx + ix+1;
	naybor[6] = ((iz*ny)+iy)*nx+nx + ix-1;
	naybor[7] = ((iz*ny)+iy)*nx-nx + ix+1;
	naybor[8] = ((iz*ny)+iy)*nx-nx + ix-1;
	naybor[9] = ((iz*ny)+ny+iy)*nx + ix;
	naybor[10] = ((iz*ny)+ny+iy)*nx + ix+1;
	naybor[11] = ((iz*ny)+ny+iy)*nx + ix-1;
	naybor[12] = ((iz*ny)+ny+iy)*nx+nx + ix;
	naybor[13] = ((iz*ny)+ny+iy)*nx-nx + ix;
	naybor[14] = ((iz*ny)+ny+iy)*nx+nx + ix+1;
	naybor[15] = ((iz*ny)+ny+iy)*nx+nx + ix-1;
	naybor[16] = ((iz*ny)+ny+iy)*nx-nx + ix+1;
	naybor[17] = ((iz*ny)+ny+iy)*nx-nx + ix-1;
	naybor[18] = ((iz*ny)-ny+iy)*nx + ix;
	naybor[19] = ((iz*ny)-ny+iy)*nx + ix+1;
	naybor[20] = ((iz*ny)-ny+iy)*nx + ix-1;
	naybor[21] = ((iz*ny)-ny+iy)*nx+nx + ix;
	naybor[22] = ((iz*ny)-ny+iy)*nx-nx + ix;
	naybor[23] = ((iz*ny)-ny+iy)*nx+nx + ix+1;
	naybor[24] = ((iz*ny)-ny+iy)*nx+nx + ix-1;
	naybor[25] = ((iz*ny)-ny+iy)*nx-nx + ix+1;
	naybor[26] = ((iz*ny)-ny+iy)*nx-nx + ix-1;

	for( inode = 0; inode < nx*ny*nz; inode ++)
	{/* loop over all mm nodes */
	/* check the origin */
	for(j=0; j< 27; j++)
	{  
	if( inode == naybor[j]) break; }
	if( j == 27  )
	{ /* then use mm */
	if( (*nodelist)[inode].innode > 0 )
	{
/* $%$%$%$%$  the expansion for f(r) goes here */
/* first the a terms */
	/*
	xt = (*nodelist)[inode].xa - a1->px;
	yt = (*nodelist)[inode].ya - a1->py;
	zt = (*nodelist)[inode].za - a1->pz;
	r = one/(xt*xt + yt*yt + zt*zt);
	r0 = sqrt(r);
	r = r*r*r;
	k = -(*nodelist)[inode].sa *a1->a*r;
	a1->VP  += k/r;
	k *= six*r0; 
	xt *= r0;
	yt *= r0;
	zt *= r0;
	a1->dpx += k*xt;
	a1->dpy += k*yt;
	a1->dpz += k*zt;
	k *= eight*r0;
	a1->qxx -= k *(xt*xt -eightth);
	a1->qxy -= k*xt*yt;
	a1->qxz -= k*xt*zt;
	a1->qyy -= k *(yt*yt -eightth);
	a1->qyz -= k*yt*zt;
	a1->qzz -= k *(zt*zt -eightth);
	*/
/* now do the multipole expansion for the electrostatic terms */
/* note that dielectric is included in the multipole expansion */
	xt = (*nodelist)[inode].xc - a1->px;
	yt = (*nodelist)[inode].yc - a1->py;
	zt = (*nodelist)[inode].zc - a1->pz;
	r = one/(xt*xt + yt*yt + zt*zt);
	r0 = sqrt(r);
	c1 =  -r*r0;
	c2 = -three*c1*r;
	c3 = -five*c2*r;
	c4 = -seven*c3*r;
	c5 = -nine*c4*r;
	xt2 = xt*xt;
	xt3 = xt2*xt;
	xt4 = xt3*xt;
	yt2 = yt*yt;
	yt3 = yt2*yt;
	yt4 = yt3*yt;
	zt2 = zt*zt;
	zt3 = zt2*zt;
	zt4 = zt3*zt;
	a1->VP += (*nodelist)[inode].sqp*a1->q*r0;
	k = c1*a1->q*xt;
	a1->VP += k*(*nodelist)[inode].q100;
	a1->dpx += k*(*nodelist)[inode].sqp;
	k = c1*a1->q*yt;
	a1->VP += k*(*nodelist)[inode].q010;
	a1->dpy += k*(*nodelist)[inode].sqp;
	k = c1*a1->q*zt;
	a1->VP += k*(*nodelist)[inode].q001;
	a1->dpz += k*(*nodelist)[inode].sqp;
/* n=2 */
	k = (c2*xt2 +c1)*a1->q;
	a1->VP += k*(*nodelist)[inode].q200;
	a1->dpx += k*(*nodelist)[inode].q100;
	a1->qxx += k*(*nodelist)[inode].sqp;
	k = (c2*yt2 +c1)*a1->q;
	a1->VP += k*(*nodelist)[inode].q020;
	a1->dpy += k*(*nodelist)[inode].q010;
	a1->qyy += k*(*nodelist)[inode].sqp;
	k = (c2*zt2 +c1)*a1->q;
	a1->VP += k*(*nodelist)[inode].q002;
	a1->dpz += k*(*nodelist)[inode].q001;
	a1->qzz += k*(*nodelist)[inode].sqp;
	k = c2*xt*yt*a1->q;
	a1->VP += k*(*nodelist)[inode].q110;
	a1->dpx += k*(*nodelist)[inode].q010;
	a1->dpy += k*(*nodelist)[inode].q100;
	a1->qxy += k*(*nodelist)[inode].sqp;
	k = c2*xt*zt*a1->q;
	a1->VP += k*(*nodelist)[inode].q101;
	a1->dpx += k*(*nodelist)[inode].q001;
	a1->dpz += k*(*nodelist)[inode].q100;
	a1->qxz += k*(*nodelist)[inode].sqp;
	k = c2*yt*zt*a1->q;
	a1->VP += k*(*nodelist)[inode].q011;
	a1->dpy += k*(*nodelist)[inode].q001;
	a1->dpz += k*(*nodelist)[inode].q010;
	a1->qyz += k*(*nodelist)[inode].sqp;
/* n=3 */
	k = (c3*xt3 +3*c2*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q300;
	a1->dpx += k*(*nodelist)[inode].q200;
	a1->qxx += k*(*nodelist)[inode].q100;
	k = (c3*yt3 +3*c2*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q030;
	a1->dpy += k*(*nodelist)[inode].q020;
	a1->qyy += k*(*nodelist)[inode].q010;
	k = (c3*zt3 +3*c2*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q003;
	a1->dpz += k*(*nodelist)[inode].q002;
	a1->qzz += k*(*nodelist)[inode].q001;
	k = (c3*xt2*yt+c2*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q210;
	a1->dpx += k*(*nodelist)[inode].q110;
	a1->dpy += k*(*nodelist)[inode].q200;
	a1->qxx += k*(*nodelist)[inode].q010;
	a1->qxy += k*(*nodelist)[inode].q100;
	k = (c3*yt2*xt+c2*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q120;
	a1->dpx += k*(*nodelist)[inode].q020;
	a1->dpy += k*(*nodelist)[inode].q110;
	a1->qyy += k*(*nodelist)[inode].q100;
	a1->qxy += k*(*nodelist)[inode].q010;
	k = (c3*xt2*zt+c2*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q201;
	a1->dpx += k*(*nodelist)[inode].q101;
	a1->dpz += k*(*nodelist)[inode].q200;
	a1->qxx += k*(*nodelist)[inode].q001;
	a1->qxz += k*(*nodelist)[inode].q100;
	k = (c3*zt2*xt+c2*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q102;
	a1->dpx += k*(*nodelist)[inode].q002;
	a1->dpz += k*(*nodelist)[inode].q101;
	a1->qzz += k*(*nodelist)[inode].q100;
	a1->qxz += k*(*nodelist)[inode].q001;
	k = (c3*yt2*zt+c2*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q021;
	a1->dpy += k*(*nodelist)[inode].q011;
	a1->dpz += k*(*nodelist)[inode].q020;
	a1->qyy += k*(*nodelist)[inode].q001;
	a1->qyz += k*(*nodelist)[inode].q010;
	k = (c3*zt2*yt+c2*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q012;
	a1->dpy += k*(*nodelist)[inode].q002;
	a1->dpz += k*(*nodelist)[inode].q011;
	a1->qzz += k*(*nodelist)[inode].q010;
	a1->qyz += k*(*nodelist)[inode].q001;
	k = (c3*zt*yt*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q111;
	a1->dpx += k*(*nodelist)[inode].q011;
	a1->dpy += k*(*nodelist)[inode].q101;
	a1->dpz += k*(*nodelist)[inode].q110;
/* n=4 */
#ifdef FOURTH
	k = (c4*xt4 +six*c3*(xt2) +three*c2)*a1->q;
	a1->VP += k*(*nodelist)[inode].q400;
	a1->dpx += k*(*nodelist)[inode].q300;
	a1->qxx += k*(*nodelist)[inode].q200;
	k = (c4*yt4 +six*c3*(yt2) +three*c2)*a1->q;
	a1->VP += k*(*nodelist)[inode].q040;
	a1->dpy += k*(*nodelist)[inode].q030;
	a1->qyy += k*(*nodelist)[inode].q020;
	k = (c4*zt4 +six*c3*(zt2) +three*c2)*a1->q;
	a1->VP += k*(*nodelist)[inode].q004;
	a1->dpz += k*(*nodelist)[inode].q003;
	a1->qzz += k*(*nodelist)[inode].q002;
	k = (c4*xt3*yt + three*c3*xt*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q310;
	a1->dpx += k*(*nodelist)[inode].q210;
	a1->dpy += k*(*nodelist)[inode].q300;
	a1->qxx += k*(*nodelist)[inode].q110;
	a1->qxy += k*(*nodelist)[inode].q200;
	k = (c4*yt3*xt + three*c3*xt*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q130;
	a1->dpx += k*(*nodelist)[inode].q030;
	a1->dpy += k*(*nodelist)[inode].q120;
	a1->qyy += k*(*nodelist)[inode].q110;
	a1->qxy += k*(*nodelist)[inode].q020;
	k = (c4*xt3*zt + three*c3*xt*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q301;
	a1->dpx += k*(*nodelist)[inode].q201;
	a1->dpz += k*(*nodelist)[inode].q300;
	a1->qxx += k*(*nodelist)[inode].q101;
	a1->qxz += k*(*nodelist)[inode].q200;
	k = (c4*zt3*yt + three*c3*xt*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q103;
	a1->dpx += k*(*nodelist)[inode].q003;
	a1->dpz += k*(*nodelist)[inode].q102;
	a1->qzz += k*(*nodelist)[inode].q101;
	a1->qxz += k*(*nodelist)[inode].q002;
	k = (c4*yt3*zt + three*c3*zt*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q031;
	a1->dpz += k*(*nodelist)[inode].q030;
	a1->dpy += k*(*nodelist)[inode].q021;
	a1->qyy += k*(*nodelist)[inode].q011;
	a1->qyz += k*(*nodelist)[inode].q020;
	k = (c4*zt3*yt + three*c3*zt*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q013;
	a1->dpz += k*(*nodelist)[inode].q012;
	a1->dpy += k*(*nodelist)[inode].q003;
	a1->qzz += k*(*nodelist)[inode].q011;
	a1->qyz += k*(*nodelist)[inode].q002;
	k = (c4*xt2*yt2 + c3*(xt2+yt2) +c2)*a1->q;
	a1->VP += k*(*nodelist)[inode].q220;
	a1->dpx += k*(*nodelist)[inode].q120;
	a1->dpy += k*(*nodelist)[inode].q210;
	a1->qxx += k*(*nodelist)[inode].q020;
	a1->qyy += k*(*nodelist)[inode].q200;
	a1->qxy += k*(*nodelist)[inode].q110;
	k = (c4*xt2*zt2 + c3*(xt2+zt2) +c2)*a1->q;
	a1->VP += k*(*nodelist)[inode].q202;
	a1->dpx += k*(*nodelist)[inode].q102;
	a1->dpz += k*(*nodelist)[inode].q201;
	a1->qxx += k*(*nodelist)[inode].q002;
	a1->qzz += k*(*nodelist)[inode].q200;
	a1->qxz += k*(*nodelist)[inode].q101;
	k = (c4*zt2*yt2 + c3*(zt2+yt2) +c2)*a1->q;
	a1->VP += k*(*nodelist)[inode].q022;
	a1->dpz += k*(*nodelist)[inode].q021;
	a1->dpy += k*(*nodelist)[inode].q012;
	a1->qzz += k*(*nodelist)[inode].q020;
	a1->qyy += k*(*nodelist)[inode].q002;
	a1->qyz += k*(*nodelist)[inode].q011;
	k = (c4*xt2*yt*zt +c3*yt*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q211;
	a1->dpz += k*(*nodelist)[inode].q210;
	a1->dpy += k*(*nodelist)[inode].q201;
	a1->dpx += k*(*nodelist)[inode].q111;
	a1->qxx += k*(*nodelist)[inode].q011;
	a1->qxy += k*(*nodelist)[inode].q101;
	a1->qyz += k*(*nodelist)[inode].q200;
	a1->qxz += k*(*nodelist)[inode].q110;
	k = (c4*xt*yt2*zt +c3*xt*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q121;
	a1->dpz += k*(*nodelist)[inode].q120;
	a1->dpy += k*(*nodelist)[inode].q111;
	a1->dpx += k*(*nodelist)[inode].q021;
	a1->qyy += k*(*nodelist)[inode].q101;
	a1->qxy += k*(*nodelist)[inode].q011;
	a1->qyz += k*(*nodelist)[inode].q110;
	a1->qxz += k*(*nodelist)[inode].q020;
	k = (c4*xt*yt*zt2 +c3*yt*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q112;
	a1->dpz += k*(*nodelist)[inode].q111;
	a1->dpy += k*(*nodelist)[inode].q102;
	a1->dpx += k*(*nodelist)[inode].q012;
	a1->qzz += k*(*nodelist)[inode].q110;
	a1->qxy += k*(*nodelist)[inode].q002;
	a1->qxz += k*(*nodelist)[inode].q011;
	a1->qyz += k*(*nodelist)[inode].q101;
#endif
/* n=5 */
#ifdef FIFTH
	k = ((c5*xt+9*c4)*xt4  +15*c3*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q500;
	a1->dpx += k*(*nodelist)[inode].q400;
	a1->qxx += k*(*nodelist)[inode].q300;
	k = ((c5*yt+9*c4)*yt4  +15*c3*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q050;
	a1->dpy += k*(*nodelist)[inode].q040;
	a1->qyy += k*(*nodelist)[inode].q030;
	k = ((c5*zt+9*c4)*zt4  +15*c3*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q005;
	a1->dpz += k*(*nodelist)[inode].q004;
	a1->qzz += k*(*nodelist)[inode].q003;
	k = (c5*xt4+six*c4*xt2 +three*c3)*yt*a1->q;
	a1->VP += k*(*nodelist)[inode].q410;
	a1->dpx += k*(*nodelist)[inode].q310;
	a1->dpy += k*(*nodelist)[inode].q400;
	a1->qxx += k*(*nodelist)[inode].q210;
	a1->qxy += k*(*nodelist)[inode].q300;
	k = (c5*yt4+six*c4*yt2 +three*c3)*xt*a1->q;
	a1->VP += k*(*nodelist)[inode].q140;
	a1->dpx += k*(*nodelist)[inode].q040;
	a1->dpy += k*(*nodelist)[inode].q130;
	a1->qyy += k*(*nodelist)[inode].q120;
	a1->qxy += k*(*nodelist)[inode].q030;
	k = (c5*xt4+six*c4*xt2 +three*c3)*zt*a1->q;
	a1->VP += k*(*nodelist)[inode].q401;
	a1->dpx += k*(*nodelist)[inode].q301;
	a1->dpz += k*(*nodelist)[inode].q400;
	a1->qxx += k*(*nodelist)[inode].q201;
	a1->qxz += k*(*nodelist)[inode].q300;
	k = (c5*zt4+six*c4*zt2 +three*c3)*xt*a1->q;
	a1->VP += k*(*nodelist)[inode].q104;
	a1->dpx += k*(*nodelist)[inode].q004;
	a1->dpz += k*(*nodelist)[inode].q103;
	a1->qzz += k*(*nodelist)[inode].q102;
	a1->qxz += k*(*nodelist)[inode].q003;
	k = (c5*yt4+six*c4*yt2 +three*c3)*zt*a1->q;
	a1->VP += k*(*nodelist)[inode].q041;
	a1->dpy += k*(*nodelist)[inode].q031;
	a1->dpz += k*(*nodelist)[inode].q040;
	a1->qyy += k*(*nodelist)[inode].q021;
	a1->qyz += k*(*nodelist)[inode].q030;
	k = (c5*zt4+six*c4*zt2 +three*c3)*yt*a1->q;
	a1->VP += k*(*nodelist)[inode].q014;
	a1->dpy += k*(*nodelist)[inode].q004;
	a1->dpz += k*(*nodelist)[inode].q013;
	a1->qzz += k*(*nodelist)[inode].q012;
	a1->qyz += k*(*nodelist)[inode].q003;
	k = (c5*xt3*yt2 +c4*(three*xt*yt2-xt3) +three*c3*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q320;
	a1->dpx += k*(*nodelist)[inode].q220;
	a1->dpy += k*(*nodelist)[inode].q310;
	a1->qxx += k*(*nodelist)[inode].q120;
	a1->qxy += k*(*nodelist)[inode].q210;
	a1->qyy += k*(*nodelist)[inode].q300;
	k = (c5*yt3*xt2 +c4*(three*yt*xt2-yt3) +three*c3*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q230;
	a1->dpx += k*(*nodelist)[inode].q130;
	a1->dpy += k*(*nodelist)[inode].q220;
	a1->qxx += k*(*nodelist)[inode].q030;
	a1->qxy += k*(*nodelist)[inode].q120;
	a1->qyy += k*(*nodelist)[inode].q210;
	k = (c5*xt3*zt2 +c4*(three*xt*zt2-xt3) +three*c3*xt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q302;
	a1->dpx += k*(*nodelist)[inode].q202;
	a1->dpz += k*(*nodelist)[inode].q301;
	a1->qxx += k*(*nodelist)[inode].q102;
	a1->qxz += k*(*nodelist)[inode].q201;
	a1->qzz += k*(*nodelist)[inode].q300;
	k = (c5*zt3*xt2 +c4*(three*zt*xt2-zt3) +three*c3*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q203;
	a1->dpx += k*(*nodelist)[inode].q103;
	a1->dpz += k*(*nodelist)[inode].q202;
	a1->qxx += k*(*nodelist)[inode].q003;
	a1->qxz += k*(*nodelist)[inode].q102;
	a1->qzz += k*(*nodelist)[inode].q201;
	k = (c5*yt3*zt2 +c4*(three*yt*zt2-yt3) +three*c3*yt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q032;
	a1->dpy += k*(*nodelist)[inode].q022;
	a1->dpz += k*(*nodelist)[inode].q031;
	a1->qyy += k*(*nodelist)[inode].q012;
	a1->qyz += k*(*nodelist)[inode].q021;
	a1->qzz += k*(*nodelist)[inode].q030;
	k = (c5*zt3*yt2 +c4*(three*zt*yt2-zt3) +three*c3*zt)*a1->q;
	a1->VP += k*(*nodelist)[inode].q023;
	a1->dpy += k*(*nodelist)[inode].q013;
	a1->dpz += k*(*nodelist)[inode].q022;
	a1->qyy += k*(*nodelist)[inode].q003;
	a1->qyz += k*(*nodelist)[inode].q012;
	a1->qzz += k*(*nodelist)[inode].q021;
	k = (c5*xt2*yt2 +c4*(xt2+yt2) +c3)*zt*a1->q;
	a1->VP += k*(*nodelist)[inode].q221;
	a1->dpx += k*(*nodelist)[inode].q121;
	a1->dpy += k*(*nodelist)[inode].q211;
	a1->dpz += k*(*nodelist)[inode].q220;
	a1->qxx += k*(*nodelist)[inode].q021;
	a1->qyy += k*(*nodelist)[inode].q201;
	a1->qxz += k*(*nodelist)[inode].q120;
	a1->qyz += k*(*nodelist)[inode].q210;
	k = (c5*xt2*zt2 +c4*(xt2+zt2) +c3)*yt*a1->q;
	a1->VP += k*(*nodelist)[inode].q212;
	a1->dpx += k*(*nodelist)[inode].q112;
	a1->dpy += k*(*nodelist)[inode].q202;
	a1->dpz += k*(*nodelist)[inode].q211;
	a1->qxx += k*(*nodelist)[inode].q012;
	a1->qzz += k*(*nodelist)[inode].q210;
	a1->qxz += k*(*nodelist)[inode].q111;
	a1->qyz += k*(*nodelist)[inode].q201;
	k = (c5*zt2*yt2 +c4*(zt2+yt2) +c3)*xt*a1->q;
	a1->VP += k*(*nodelist)[inode].q122;
	a1->dpx += k*(*nodelist)[inode].q022;
	a1->dpy += k*(*nodelist)[inode].q112;
	a1->dpz += k*(*nodelist)[inode].q121;
	a1->qzz += k*(*nodelist)[inode].q120;
	a1->qyy += k*(*nodelist)[inode].q102;
	a1->qxy += k*(*nodelist)[inode].q022;
	a1->qxz += k*(*nodelist)[inode].q022;
	k = (c5*xt3+three*c4*xt)*yt*zt*a1->q;
	a1->VP += k*(*nodelist)[inode].q311;
	a1->dpx += k*(*nodelist)[inode].q211;
	a1->dpy += k*(*nodelist)[inode].q301;
	a1->dpz += k*(*nodelist)[inode].q310;
	a1->qxx += k*(*nodelist)[inode].q211;
	a1->qxy += k*(*nodelist)[inode].q201;
	a1->qxz += k*(*nodelist)[inode].q210;
	k = (c5*yt3+three*c4*yt)*xt*zt*a1->q;
	a1->VP += k*(*nodelist)[inode].q131;
	a1->dpx += k*(*nodelist)[inode].q031;
	a1->dpy += k*(*nodelist)[inode].q121;
	a1->dpz += k*(*nodelist)[inode].q130;
	a1->qyy += k*(*nodelist)[inode].q111;
	a1->qxy += k*(*nodelist)[inode].q021;
	a1->qyz += k*(*nodelist)[inode].q120;
	k = (c5*zt3+three*c4*zt)*yt*xt*a1->q;
	a1->VP += k*(*nodelist)[inode].q113;
	a1->dpx += k*(*nodelist)[inode].q013;
	a1->dpy += k*(*nodelist)[inode].q103;
	a1->dpz += k*(*nodelist)[inode].q112;
	a1->qzz += k*(*nodelist)[inode].q111;
	a1->qyz += k*(*nodelist)[inode].q102;
	a1->qxz += k*(*nodelist)[inode].q012;
#endif

	} /* if innode > 0  end if */
	} else if( (*nodelist)[inode].innode > 0) 
		{ /* if not mm use explicit */
/* first load the atoms onto atomall */
	imax = 0;
	i = (*nodelist)[inode].first;
	 if( (*nodelist)[inode].innode > 0  && 
	    ((*atomlist)[i].who)->serial > a1->serial)
	 {  (*atomall)[imax++] = (*atomlist)[i].who;}
	for( j=1; j< (*nodelist)[inode].innode -1 ; j++)
	{
	 i = (*atomlist)[i].next;
	 if( ((*atomlist)[i].who)->serial > a1->serial)
	 {  (*atomall)[imax++] = (*atomlist)[i].who;}
	}
	/*
	for( j=0; j< jj; j++)
	{
		if( (*atomlist)[j].which == inode && 
		    (*atomlist)[j].who->serial > a1->serial)
		    {(*atomall)[imax] = (*atomlist)[j].who;
		    imax+= 1;}
	}
	*/
#pragma _CNX no_recurrence
	for( i=0; i< imax; i++)
	{
	a2 = (*atomall)[i];
	j = i*4;
	(*vector)[j  ] = a2->px - a1->px ;
	(*vector)[j+1] = a2->py - a1->py ;
	(*vector)[j+2] = a2->pz - a1->pz ;
	}
#pragma _CNX no_recurrence
	for( i=0; i< imax; i++)
	{
		j = i*4;
		(*vector)[j+3] = sqrt((*vector)[j]*(*vector)[j] +
				 (*vector)[j+1]*(*vector)[j+1] +
				 (*vector)[j+2]*(*vector)[j+2]);
	}
/* add the new components */
	for( i=0; i< imax; i++)
	{
	a2 = (*atomall)[i];
	for( j=0; j< a1->dontuse; j++)
	{ if( a2 == a1->excluded[j]) goto SKIPNEW;}
	j = i*4;
	if( (*vector)[j+3] > mxcut || inclose > NCLOSE )
	{
	r0 = one/(*vector)[j+3];  
	r = r0*r0;
	r4 = r*r;
	r = r*r*r; /* r0^-6 */
	xt = a1->q*a2->q*dielectric*r0;
	pt = 0.5*(a1->a*a2->q*a2->q +a1->q*a1->q*a2->a)*r4;
	yt = a1->a*a2->a*r;
	zt = a1->b*a2->b*r*r;
	k = xt - yt + zt -pt;
	xt = xt*r0; yt = yt*r0; zt = zt*r0; pt = pt*r0;
	k1 = xt - yt*six + zt*twelve -pt*3;
	xt = xt*r0; yt = yt*r0; zt = zt*r0; pt = pt*r0;
	k2 = xt*three; ka2 = - yt*six*eight; kb2 =  zt*twelve*14; kp2 = -pt*24;
#ifdef CUBIC
	xt = xt*r0; yt = yt*r0; zt = zt*r0;
	k3 = -xt*5*3; ka3 =   yt*6*8*10 ; kb3 =  -zt*12*14*16; kp3 = pt*4*6*8;
#endif
	k1 = -k1;
	xt = (*vector)[j]*r0 ;
	yt = (*vector)[j+1]*r0 ;
	zt = (*vector)[j+2] *r0;
	/*
	xt = (*vector)[j] ;
	yt = (*vector)[j+1] ;
	zt = (*vector)[j+2] ;
	*/
	a1->VP += k;
	a2->dpx -= k1*xt;
	a1->dpx += k1*xt;
	a2->dpy -= k1*yt;
	a1->dpy += k1*yt;
	a2->dpz -= k1*zt;
	a1->dpz += k1*zt;
	xt2 = xt*xt; yt2 = yt*yt; zt2 = zt*zt;
	a2->qxx -= k2*(xt2 - third) + ka2*(xt2 - eightth)+kb2*(xt2-fourteenth) +kp2*(xt*xt - 0.16666666666);
	a1->qxx -= k2*(xt2 - third) + ka2*(xt2 - eightth)+kb2*(xt2-fourteenth) +kp2*(xt*xt - 0.16666666666);
	a2->qxy -= (k2+ka2+kb2+kp2)*yt*xt;
	a1->qxy -= (k2+ka2+kb2+kp2)*yt*xt;
	a2->qxz -= (k2+ka2+kb2+kp2)*zt*xt;
	a1->qxz -= (k2+ka2+kb2+kp2)*zt*xt;
	a2->qyy -= k2*(yt2 - third) + ka2*(yt2 - eightth)+kb2*(yt2-fourteenth) +kp2*(yt2-0.1666666666666666666);
	a1->qyy -= k2*(yt2 - third) + ka2*(yt2 - eightth)+kb2*(yt2-fourteenth) +kp2*(yt2-0.16666666666666666);
	a2->qyz -= (k2+ka2+kb2+kp2)*yt*zt;
	a1->qyz -= (k2+ka2+kb2+kp2)*yt*zt;
	a2->qzz -= k2*(zt2 - third) + ka2*(zt2 - eightth)+kb2*(zt2-fourteenth) +kp2*(zt2-0.16666666666666666666);
	a1->qzz -= k2*(zt2 - third) + ka2*(zt2 - eightth)+kb2*(zt2-fourteenth) +kp2*(zt2-0.166666666666666666666);
#ifdef CUBIC
	a2->qxxx -= k3*(xt*xt*xt - xt*( 9./15 )) ;
	a2->qxxx -= ka3*(xt*xt*xt - xt*( 24./80 )) ;
	a2->qxxx -= kb3*(xt*xt*xt - xt*( 42./(14*18))); 
	a1->qxxx += k3*(xt*xt*xt - xt*( 9./15 )) ;
	a1->qxxx += ka3*(xt*xt*xt - xt*( 24./80 )) ;
	a1->qxxx += kb3*(xt*xt*xt - xt*( 42./(14*18))); 
	a2->qxxy -= k3*(yt*xt*xt - yt*( 6./ 15));
	a2->qxxy -= ka3*(yt*xt*xt - yt*( 11./ 80));
	a2->qxxy -= kb3*(yt*xt*xt - yt*( 17./ (14*18)));
	a1->qxxy += k3*(yt*xt*xt - yt*( 6./ 15));
	a1->qxxy += ka3*(yt*xt*xt - yt*( 11./ 80));
	a1->qxxy += kb3*(yt*xt*xt - yt*( 17./ (14*18)));
	a2->qxxz -= k3*(zt*xt*xt - zt*( 6./ 15));
	a2->qxxz -= ka3*(zt*xt*xt - zt*( 11./ 80));
	a2->qxxz -= kb3*(zt*xt*xt - zt*( 17./ (14*18)));
	a1->qxxz += k3*(zt*xt*xt - zt*( 6./ 15));
	a1->qxxz += ka3*(zt*xt*xt - zt*( 11./ 80));
	a1->qxxz += kb3*(zt*xt*xt - zt*( 17./ (14*18)));
	a2->qxyy -= k3*(yt*yt*xt - xt*( 6./ 15));
	a2->qxyy -= ka3*(yt*yt*xt - xt*( 11./ 80));
	a2->qxyy -= kb3*(yt*yt*xt - xt*( 17./ (14*18)));
	a1->qxyy += k3*(yt*yt*xt - xt*( 6./ 15));
	a1->qxyy += ka3*(yt*yt*xt - xt*( 11./ 80));
	a1->qxyy += kb3*(yt*yt*xt - xt*( 17./ (14*18)));
	a2->qxyz -= (k3+ka3+kb3)*yt*zt*xt;
	a1->qxyz += (k3+ka3+kb3)*yt*zt*xt;
	a2->qxzz -= k3*(zt*zt*xt - xt*( 6./ 15));
	a2->qxzz -= ka3*(zt*zt*xt - xt*( 11./ 80));
	a2->qxzz -= kb3*(zt*zt*xt - xt*( 17./ (14*18)));
	a1->qxzz += k3*(zt*zt*xt - xt*( 6./ 15));
	a1->qxzz += ka3*(zt*zt*xt - xt*( 11./ 80));
	a1->qxzz += kb3*(zt*zt*xt - xt*( 17./ (14*18)));
	a2->qyyy -= k3*(yt*yt*yt - yt*( 9./15 )) ;
	a2->qyyy -= ka3*(yt*yt*yt - yt*( 24./80 )) ;
	a2->qyyy -= kb3*(yt*yt*yt - yt*( 42./(14*18))); 
	a1->qyyy += k3*(yt*yt*yt - yt*( 9./15 )) ;
	a1->qyyy += ka3*(yt*yt*yt - yt*( 24./80 )) ;
	a1->qyyy += kb3*(yt*yt*yt - yt*( 42./(14*18))); 
	a2->qyyz -= k3*(yt*yt*zt - zt*( 6./ 15));
	a2->qyyz -= ka3*(yt*yt*zt - zt*( 11./ 80));
	a2->qyyz -= kb3*(yt*yt*zt - zt*( 17./ (14*18)));
	a1->qyyz += k3*(yt*yt*zt - zt*( 6./ 15));
	a1->qyyz += ka3*(yt*yt*zt - zt*(11./ 80));
	a1->qyyz += kb3*(yt*yt*zt - zt*( 17./ (14*18)));
	a2->qyzz -= k3*(zt*zt*yt - yt*( 6./ 15));
	a2->qyzz -= ka3*(zt*zt*yt - yt*( 11./ 80));
	a2->qyzz -= kb3*(zt*zt*yt - yt*( 17./ (14*18)));
	a1->qyzz += k3*(zt*zt*yt - yt*( 6./ 15));
	a1->qyzz += ka3*(zt*zt*yt - yt*( 11./ 80));
	a1->qyzz += kb3*(zt*zt*yt - yt*( 17./ (14*18)));
	a2->qzzz -= k3*(zt*zt*zt - zt*( 9./15 )) ;
	a2->qzzz -= ka3*(zt*zt*zt - zt*( 24./80 )) ;
	a2->qzzz -= kb3*(zt*zt*zt - zt*( 42./(14*18))); 
	a1->qzzz += k3*(zt*zt*zt - zt*( 9./15 )) ;
	a1->qzzz += ka3*(zt*zt*zt - zt*( 24./80 )) ;
	a1->qzzz += kb3*(zt*zt*zt - zt*( 42./(14*18))); 
#endif
	}else {
	a1->close[inclose++] = (*atomall)[i];
/* debugging
	j = i *4;
	fprintf(stderr," mxcut %f %f inclose %d who %d \n",mxcut,(*vector)[j+3],inclose,(*atomall)[i]->serial);
	fprintf(stderr," vector %f %f %f \n", (*vector)[j],(*vector)[j+1],(*vector)[j+2]);
*/
	if( inclose == NCLOSE)
	{
	aaerror(
	" fv_update_nonbon> too many atoms increase NCLOSE or decrease mxcut");
/*		exit(0); 
*/
	}
	}
SKIPNEW:  j =  j;
	}/* end of loop i */
	}/* end of if in close MM node */
	}/* end of loop inode */
/* merge the non-bond mxcut lists */
/*	a1->close[inclose] = NULL; */
/* set the position */
	a1->px = a1->dx*lambda + a1->x;
	a1->py = a1->dy*lambda + a1->y;
	a1->pz = a1->dz*lambda + a1->z;

	}  /* end of ii loop */
	
	a_inactive_f_zero();

	free( atomlist);
	free( nodelist);
	free( vector);
	free (atomall);
	return 1;

}


