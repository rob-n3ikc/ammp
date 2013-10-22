/*app.c
*
*  implement an "arbitrary Polar Potential"
*    V = k/r^n where k,N are user specified.
*
* app terms are not re-writable by default.
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of app bonds, with the associated length and force constants
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are bond wise symmetric - so we don't have to fuck around with
* s matrices and the like.
*/
/*
*  copyright 1992,2003 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-apped potentials
*/
typedef struct{
    ATOM *atom1,*atom2;
    int power;
    float konstant;
    void *next;
}  APP;
#define MLONG sizeof(APP)

static APP *app_first = NULL;
static APP *app_last = NULL;
/* function app adds a app to the app list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial numbers, length and constant
* allocates the papp memory, initializes it and
* returns
*/
#ifdef newcall
int app_rewrite( int p1,int p2,float fk, int order)
#else
int app_rewrite(  p1, p2, fk,  order)
	int p1,p2,order;
	float fk;
#endif
{
    ATOM *ap1,*ap2,*a_m_serial(int);
    APP *papp;
//    int app( int, int, float, int);
    ap1 = a_m_serial(p1);
    ap2 = a_m_serial(p2);
    papp = app_first;

    while(  papp != (APP*)NULL)
    {
	    if( ap1 == papp->atom1 && ap2 == papp->atom2)
	    {
		    papp->konstant = fk;
	    	    papp->power = order;
		    return 0;
	    }
	    if( papp == papp->next) break;
	    papp = papp->next;
    }
    /* if we didn't need to re-write create the papp one. */
    return app( p1,p2,fk,order);
}
#ifdef newcall
int app( int p1,int p2,float fk, int order)
#else
int app(  p1, p2, fk,  order)
	int p1,p2,order;
	float fk;
#endif
{
    ATOM *ap1,*ap2,*a_m_serial(int);
    APP *papp;
    char line[80];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    if( (ap1 == NULL) || (ap2 == NULL) )
    {
        sprintf( line,"undefined atom in app %d %d \0",p1,p2);
        aaerror( line );
        return 0;
    }

    if( ( papp = malloc( MLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( app_first == NULL) app_first = papp;
    if( app_last == NULL) app_last = papp;
    papp -> atom1 = ap1;
    papp -> atom2 = ap2;
    papp -> konstant = fk;
    papp ->  power = order;
    papp -> next = papp;
    app_last -> next = papp;
    app_last = papp;
    return 1;
}


/* v_app()
* this function sums up the potentials
* for the atoms defined in the APP data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
#ifdef newcall
int v_app(float  *V, float lambda )
#else
int v_app(V,  lambda )
	float *V, lambda;
#endif
{
    APP *bp;
    float r,xt,yt,zt;
    ATOM *a1,*a2;
    int i;


    bp = app_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active){
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
	    xt = 1.;
	    r = 1./r;
	    for( i=0; i< bp->power; i++)
		    xt = xt*r;
            /* *V += bp->k*( r - bp->length)*(r - bp->length);
            */
            *V += bp->konstant*xt;
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_app()
*
* f_app increments the forces in the atom structures by the force
* due to the app components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
#ifdef newcall
int f_app(float lambda)
#else
int f_app(lambda)
	float lambda;
#endif
/*  returns 0 if error, 1 if OK */
{
    APP *bp;
    float r,k,ux,uy,uz;
    float xt;
    int i;
    ATOM *a1,*a2;


    bp = app_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active){
            if( lambda == 0.)
            {
                ux = (a2->x - a1->x);
                uy = (a2->y - a1->y);
                uz = (a2->z - a1->z);
            }else{
                ux = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
                uy = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
                uz = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            }
            r = ux*ux + uy*uy + uz*uz;
            /* watch for FP errors*/
            if( r <= 1.e-5)
            { r = 0; ux = 1.; uy = 0.; uz = 0.; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
	    	r = 1./r;
            }
	    xt = 1.;
	    for( i=0; i< bp->power; i++)
		    xt = xt*r;
            ux =   -bp->power*bp->konstant*r*ux*xt;
            uy =   -bp->power*bp->konstant*r*uy*xt;
            uz =   -bp->power*bp->konstant*r*uz*xt;
            if( a1->active){
                a1->fx += ux;
                a1->fy += uy;
                a1->fz += uz;
            }
            if( a2->active){
                a2->fx -= ux;
                a2->fy -= uy;
                a2->fz -= uz;
            }
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* routine dump_apps
* this function outputs the app parameters
* and does it in a simple form
* app ser1,ser2,k,req
* the rest is just free format
*/
#ifdef newcall
void dump_app(FILE *where )
#else
void dump_app(where )
	FILE *where;
#endif
{
    APP *b;
    ATOM *a1,*a2;
    b = app_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;
        fprintf( where,"app %d %d %f %d ;\n",a1->serial,a2->serial,
                b->konstant,b->power);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;
    fprintf( where,"app %d %d %f %d ;\n",a1->serial,a2->serial,
             b->konstant,b->power);
}
