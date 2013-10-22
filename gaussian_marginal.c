/* gaussian_marginal.c
*
* collection of routines to service gaussian_marginal length potentials
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of gaussian_marginal, with the associated length and force constant
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are gaussian_marginal wise symmetric - so we don't have to fuck around with
* s matrices and the like.
*/
/*
*  copyright 1992,1993,1994,1995 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-gaussian_marginaled potentials
*/
/* this struct defines upto 3 gaussians in the sum */
typedef struct{
    ATOM *atom1,*atom2;
    float radius[3],sigma[3],k[3]; // target distance, width, konstant
    int nterms;
    void *next;
}  GAUSSIAN_MARGINAL;
#define BLONG sizeof(GAUSSIAN_MARGINAL)

GAUSSIAN_MARGINAL *gaussian_marginal_first = NULL;
GAUSSIAN_MARGINAL *gaussian_marginal_last = NULL;
/* function gaussian_marginal adds a gaussian_marginal to the gaussian_marginal list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int gaussian_marginal( p1,p2,n,k1,r1,s1,k2,r2,s2,k3,r3,s3)
int p1,p2;
int n;
float k1,k2,k3;
float r1,r2,r3;
float s1,s2,s3;
{
    ATOM *ap1,*ap2,*a_m_serial();
    GAUSSIAN_MARGINAL *new;
    int i;
    char line[80];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        sprintf( line,"undefined atom in gaussian_marginal %d %d \0",p1,p2);
        aaerror( line );
        return 0;
    }

    if( ( new = malloc( BLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( gaussian_marginal_first == NULL) gaussian_marginal_first = new;
    if( gaussian_marginal_last == NULL) gaussian_marginal_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> nterms = n;
// not worth trying to be clever here
    new -> k[0] = k1;
    new -> k[1] = k2;
    new -> k[2] = k3;
    new -> radius[0] = r1;
    new -> radius[1] = r2;
    new -> radius[2] = r3;
    new -> sigma[0] = s1;
    new -> sigma[1] = s2;
    new -> sigma[2] = s3;
    new -> next = new;

    gaussian_marginal_last -> next = new;
    gaussian_marginal_last = new;
    return 1;
}


/* v_gaussian_marginal()
* this function sums up the potentials
* for the atoms defined in the GAUSSIAN_MARGINAL data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_gaussian_marginal( V, lambda )
float *V,lambda;
{
    GAUSSIAN_MARGINAL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;
    int i;


    bp = gaussian_marginal_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
        {
            if( lambda == 0.)
            {
                r = (a1->x - a2->x)*(a1->x - a2->x);
                r = r + (a1->y - a2->y)*(a1->y - a2->y);
                r = r + (a1->z - a2->z)*(a1->z - a2->z);
            } else
            {
                xt = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
                yt = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
                zt = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
                r = xt*xt+yt*yt+zt*zt;
            }
            r = sqrt(r); 
            for( i=0; i< bp->nterms ; i++)
                 *V += bp->k[i]*exp( -( r-bp->radius[i])*(r-bp->radius[i])/bp->sigma[i]);	
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_gaussian_marginal()
*
* f_gaussian_marginal increments the forces in the atom structures by the force
* due to the gaussian_marginal components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_gaussian_marginal(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    GAUSSIAN_MARGINAL *bp;
    float r,k,ux,uy,uz;
    float dux,duy,duz;
    ATOM *a1,*a2;
    int i;


    bp = gaussian_marginal_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
        {
            if( lambda == 0.)
            {
                ux = (a2->x - a1->x);
                uy = (a2->y - a1->y);
                uz = (a2->z - a1->z);
            }else{
                ux = (a2->x -a1->x +lambda*(a2->dx-a1->dx));
                uy = (a2->y -a1->y +lambda*(a2->dy-a1->dy));
                uz = (a2->z -a1->z +lambda*(a2->dz-a1->dz));
            }
            r = ux*ux + uy*uy + uz*uz;
            /* watch for FP errors*/
            if( r <= 1.e-15)
            { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            dux = 0.;
            duy = 0.;
            duz = 0.;
            for( i=0; i< bp->nterms ; i++)
		{
                 float dr;
                 dr = bp->k[i]*exp( -( r-bp->radius[i])*(r-bp->radius[i])/bp->sigma[i]);	
                 dr *= -(r-bp->radius[i])/bp->sigma[i]*2.;
                 dux += dr*ux;
                 duy += dr*uy;
                 duz += dr*uz;
		}
            if( a1->active ){
                a1->fx += dux;
                a1->fy += duy;
                a1->fz += duz;
            }
            if( a2->active ){
                a2->fx -= dux;
                a2->fy -= duy;
                a2->fz -= duz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}

/* routine dump_gaussian_marginal
* this function outputs the gaussian_marginal parameters
* and does it in a simple form
* gaussian_marginal ser1,ser2,k,req
* the rest is just free format
*/
void dump_gaussian_marginal( where )
FILE *where;
{
    GAUSSIAN_MARGINAL *b;
    ATOM *a1,*a2;
    int i;
    b = gaussian_marginal_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"marginal %d %d %d ",a1->serial,a2->serial,b->nterms);
        for( i=0; i< b->nterms; i++)
           fprintf(where,"%f %f %f ",b->k[i],b->radius[i],b->sigma[i]);
        fprintf(where,";\n"); 
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"marginal %d %d %d ",a1->serial,a2->serial,b->nterms);
        for( i=0; i< b->nterms; i++)
           fprintf(where,"%f %f %f ",b->k[i],b->radius[i],b->sigma[i]);
        fprintf(where,";\n"); 
}


/* a_gaussian_marginal()
* this function sums up the potentials
* for the atoms defined in the GAUSSIAN_MARGINAL data structure.
* only does gaussian_marginal in the given range
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_gaussian_marginal( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    GAUSSIAN_MARGINAL *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;
    int i;


    bp = gaussian_marginal_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active )
            if(( a1->serial >= ilow && a1->serial <=ihigh)
                    ||( a2->serial >= ilow && a2->serial <=ihigh))
            {
                if( lambda == 0.)
                {
                    r = (a1->x - a2->x)*(a1->x - a2->x);
                    r = r + (a1->y - a2->y)*(a1->y - a2->y);
                    r = r + (a1->z - a2->z)*(a1->z - a2->z);
                } else
                {
                    xt = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
                    yt = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
                    zt = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
                    r = xt*xt+yt*yt+zt*zt;
                }
                r = sqrt(r);  zt = 0.;
            for( i=0; i< bp->nterms ; i++)
                 zt += bp->k[i]*exp( -( r-bp->radius[i])*(r-bp->radius[i])/bp->sigma[i]);	
                *V += zt;
                fprintf(op,"Bond %s %d %s %d E %f value \n"
                        ,a1->name,a1->serial,a2->name,a2->serial,zt,r);
            }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}


/* these are inherited from bonds.c but may be useful so I've kept them */

GAUSSIAN_MARGINAL *get_gaussian_marginal_pointer( a1,a2)
ATOM *a1, *a2;
{
    GAUSSIAN_MARGINAL *bp;

    bp = gaussian_marginal_first;
    while(1)
    {
        if( bp == NULL ) return;
        if( a1 == bp->atom1 && a2 == bp->atom2) return bp;
        if( a1 == bp->atom2 && a2 == bp->atom1) return bp;
        if( bp == bp->next) return NULL;
        bp = bp->next;
    }
    return NULL;
}

int gaussian_marginal_next( i,ap1,ap2)
int i;
ATOM **ap1,**ap2;
{
    static GAUSSIAN_MARGINAL *bp;
    if( i <= 0){ bp = gaussian_marginal_first;
        if( bp == NULL ){*ap1 = NULL; *ap2 == NULL; return 0; }}
    *ap1 = bp->atom1;
    *ap2 = bp->atom2;
    if( bp->next == bp || bp == NULL) return 0;
    if( bp->next != NULL)bp = bp->next;

    return 1;
}
