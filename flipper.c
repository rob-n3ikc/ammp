/* flipper.c is an experimental (right now hard coded for protein)
*  simulated annealing solver for peptide (and inter monomer in general)
*  configuration searching.
*
*  It will flip along the (CA-CA) axis, when possible (NOT PROLINE) the peptide
* groups followed by minimization (as in mcmin etc).
*
*/

/* include everything you might need? */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include "ammp.h"

// external definitions we will use
int math_match_atom( char*, ATOM* );
void get_bond(ATOM*, ATOM *[], int, int*);
/* clone routines from clone.c */
int clone( int);
int restore_from_clone(int);


int rotate_around_axis( a,b,what,how_many,how_much)
ATOM *a,*b, *what[];
float how_much;
int how_many;
{


    float nx,ny,nz;
    float phi,cphi,sphi,xphi;
    float rx,ry,rz, nnrx,nnry,nnrz, rnx,rny,rnz;
    ATOM *ap;
    int i;

    nx = b->x - a->x;
    ny = b->y - a->y;
    nz = b->z - a->z;
    rx = sqrt(nx*nx + ny*ny + nz*nz);
    if( rx < 1.e-6)
    {aaerror(" bad torsion radius in rotate_around_axis \n"); return 0 ;}
    nx = nx/rx;
    ny = ny/rx;
    nz = nz/rx;
    /*
            phi = 2.*3.141592653589793 /(float)tgp->ntry * (float)num;
    */
    cphi = cos(how_much); sphi = sin(how_much);
    phi = (1.-cphi);
    for( i=0; i< how_many; i++)
    {
        ap = what[i];
            rx = ap->x - a->x;
            ry = ap->y - a->y;
            rz = ap->z - a->z;
            xphi = nx*rx + ny*ry + nz*rz;
            nnrx = xphi*nx;
            nnry = xphi*ny;
            nnrz = xphi*nz;
            rnx = ny*rz - nz*ry;
            rny = -nx*rz + nz*rx;
            rnz = nx*ry - ny*rx;
            rx = cphi*rx + phi*nnrx + sphi*rnx;
            ry = cphi*ry + phi*nnry + sphi*rny;
            rz = cphi*rz + phi*nnrz + sphi*rnz;
            ap->x = rx + a->x;
            ap->y = ry + a->y;
            ap->z = rz + a->z;
	}
	return 0;
}





int flipper( op,vfs,ffs, nfs, at1,at2, temper, nstep)
	FILE *op;
	int (*vfs[])(), (*ffs[])(), nfs, nstep;
	float temper;
	char *at1,*at2;
{
	ATOM *(*aps)[], *a_next();
	ATOM *a1,*a2,*a3,*a4;
	ATOM *bonded[10];
	int inbonded;
	int a_number(), numatom;
	int i,istep;
	float ddr,randf();
	float V,Vbest, dV,dt;	
	int ica,jca,ic,io,in,ihn,istart;
	int ires, jres;

	numatom = a_number();
	if( numatom < 4) return 0;
	aps = malloc( numatom * sizeof( ATOM *));
	if( aps == NULL) { aaerror("FLIPPER error exhausted\n"); exit(0);}
	for(i = 0; i< numatom; i++)
		(*aps)[i] = a_next(i);
// dt will be pi.  Period.  (we can try others, but this should work)
	dt = 3.1415926535;

	printf("entering flipper\n"); fflush(stdout);
	clone( 1);  /* arbitrary integer label */
	Vbest = 10.e10;

	for( istep=0; istep< nstep; istep++)
	{
		int in_bond;
		int j;
		float dx;
	V = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&V,0.);
	if( Vbest > V) Vbest = V;
	clone( 2);  /* arbitrary integer label */
//	printf("flipper past clone\n"); fflush(stdout);
//
//		tset_bond_build( a1,a2,a3,a4);
//		set_torsion( a1,a2,a3,a4, dV );
//	fprintf(op,"%f %f \n",dV,get_torsion_value( a1,a2,a3,a4 ));
//	}
//
// pick a random atom, if it isn't ca go up until it is
// if we get to the end then do it again
// if it is the last one (there isn't one after) do it again.
//
// won't make the assumption that the atoms are consequitive.
// ARRGH!
RECHOOSE: ;;
	istart = numatom *randf();
	if( istart < 0) goto RECHOOSE;
	if( istart >= numatom) goto RECHOOSE;
//	printf("flipper choose %d %d\n", istart,(*aps)[istart]->serial); fflush(stdout);
	ires = (*aps)[istart]->serial/100;
	jres = ires +1; // but will assume that monomers are
// now look for the ca, c, o in residue i
	ica = -1;
	ic = -1;
	io = -1;
	in = -1;
	ihn = -1;
	jca = -1;
	for( i=0;i< numatom; i++)
	{
	int k;
	k = (*aps)[i]->serial/100;
	if( k == ires)
	{
	if( math_match_atom( "ca", (*aps)[i]) != 0){ ica = i;}
	if( math_match_atom( "c", (*aps)[i]) != 0){ ic = i;}
	if( math_match_atom( "o", (*aps)[i]) != 0){ io = i;}
	}
	if( k == jres)
	{
	if( math_match_atom( "ca", (*aps)[i]) != 0){ jca = i;}
	if( math_match_atom( "n", (*aps)[i]) != 0){ in = i;}
	if( math_match_atom( "hn", (*aps)[i]) != 0){ ihn = i;}
	}

	}
	if(  ica < 0) goto RECHOOSE;
	if(  ic < 0) goto RECHOOSE;
	if(  io < 0) goto RECHOOSE;
	if(  in < 0) goto RECHOOSE;
	if(  ihn < 0) goto RECHOOSE;
	if(  jca < 0) goto RECHOOSE;

// so now we can flip it.
	bonded[0] = (*aps)[ic];
	bonded[1] = (*aps)[io];
	bonded[2] = (*aps)[in];
	bonded[3] = (*aps)[ihn];
	inbonded = 4;

//	printf("flipper rotating\n"); fflush(stdout);
//	ddr = randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr += randf() -0.5;
//	ddr *= 2./10.;
//	rotate_around_axis( (*aps)[ica], (*aps)[jca], bonded,inbonded,dt*ddr);
	rotate_around_axis( (*aps)[ica], (*aps)[jca], bonded,inbonded,dt);
//	printf("flipper rotated\n"); fflush(stdout);
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer( vfs,ffs,nfs, 100,30, (1==0));
//		block_optimizer( vfs,ffs,nfs, 10,30, (1==0));
//standard		block_optimizer( vfs,ffs,nfs, 10,a_number()/20, (1==0));
		pc_cngdel( vfs,ffs,nfs, 100,100,0., (1==0));
//		pc_cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
//was 10,10		cngdel( vfs,ffs,nfs, 100,100,0., (1==0));
	dV = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&dV,0.);
	if( isnan(dV) )
	{
	fprintf(op,"error condition - restoring from previous best\n");

		restore_from_clone(1);
		V = 0.;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
		dV = V;
	} else{
	dx = exp(  (V-dV)/temper*1.987);
	fprintf(op, "V %f, dV, %f, Delta %f %f\n", V,dV, -V+dV,dx);
	if( dV < Vbest){ clone(1); Vbest = dV;}
	if( dV < V){ clone(2);}
	else {
		if( dx < randf(-1)){fprintf(op,"rejected\n"); restore_from_clone(2);}
	}
	}
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
	}/* istep */

	restore_from_clone(1);
	free( aps);
}
