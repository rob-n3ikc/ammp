/* implementation of Monte-Carlo /minimization
 * i.e. score for simulated annealling after minimization
 *  depends on the clone routines
 *  
* also depends on tset  actually using tset_bond_build and set_torsion
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ammp.h"

#ifdef graphic
void force_screen_update();
#endif

/* clone routines from clone.c */

int clone( int);
int restore_from_clone(int);
int tset_bond_build( ATOM*, ATOM* , ATOM*, ATOM*);
//int set_torsion( ATOM* , ATOM*, ATOM*, ATOM*, float);
float get_torsion_value( ATOM* , ATOM*, ATOM*, ATOM*);
int math_match_atom( char*, ATOM*);
void get_bond( ATOM*, ATOM *[], int, int*);
float get_f_variable(char *);

int mcmin( op,vfs,ffs, nfs, at1,at2, temper, nstep)
	FILE *op;
	int (*vfs[])(), (*ffs[])(), nfs, nstep;
	float temper;
	char *at1,*at2;
{
	ATOM *(*aps)[], *a_next();
	ATOM *a1,*a2,*a3,*a4;
	ATOM *bonded[10];
	int a_number(), numatom;
	int i,istep;
	float randf();
	float V,Vbest, dV,dt;	
	float mxdq;

	numatom = a_number();
	if( numatom < 4) return 0;
	aps = malloc( numatom * sizeof( ATOM *));
	if( aps == NULL) { aaerror("MCMIN error exhausted\n"); exit(0);}
	for(i = 0; i< numatom; i++)
		(*aps)[i] = a_next(i);
	dt = 3.1415926535*0.5;
	// special version for 2pi/3 torsion shifts;
//	dt = 3.1415926535*2./3.;
	// special version for backbone etc
	dt = 3.1415926535*0.1;

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
	if( Vbest > V)	Vbest = V;
	for( j = 0; j< 1; j++)
	{
	clone( 2);  /* arbitrary integer label */
RE_CHOOSE: ;
		i = numatom *randf();
		if( i < 0 ) i = 0; 
		if( i >= numatom) i = numatom-1;
		a2 = NULL;
		for( ; i< numatom; i++)
			if( math_match_atom( at1,(*aps)[i]) != 0)
			{a2 = (*aps)[i]; break;}
		if( a2 == NULL) goto RE_CHOOSE;
/* RARE PROTEIN DEPENDENCY HERE */
/* skip disulfides */
                if( strcmp( a2->name, "css.ca") == 0) goto RE_CHOOSE;
		get_bond( a2,bonded, 10, &in_bond);
		a3 = NULL;
		for( i =0; i< in_bond; i++)
			if( math_match_atom( at2,bonded[i]) != 0)
			{a3 = bonded[i]; break;}
		if( a3 == NULL) goto RE_CHOOSE;
		a1 = NULL;
		for( i=0; i< in_bond; i++)
			if( math_match_atom( at2,bonded[i]) == 0)
			{a1 = bonded[i]; break;}
		if( a1 == NULL) goto RE_CHOOSE;
		get_bond( a3,bonded, 10, &in_bond);
		a4 = NULL;
		for( i=0; i< in_bond; i++)
			if( math_match_atom( at1,bonded[i]) == 0 )
			{a4 = bonded[i]; break;}
			if( a4 == a1) goto RE_CHOOSE; /* its a loop */
			if( a4 == NULL) goto RE_CHOOSE;
		get_bond( a4,bonded, 10, &in_bond);
		for( i=0; i< in_bond; i++)
			if( bonded[i] == a1) goto RE_CHOOSE;
		fprintf(op, "%s %s %s %s\n", a1->name, a2->name, a3->name, a4->name);
		fflush(op);
/* so now we have a1,a2,a3,a4 which define a valid torsion */
	fprintf(op,"%f ",get_torsion_value( a1,a2,a3,a4 ));
	dV = dt*(randf()-0.5);
//	dV = dt; if (randf() < 0.5) dV = -dV;
		tset_bond_build( a1,a2,a3,a4);
		set_torsion( a1,a2,a3,a4, dV );
	fprintf(op,"%f %f \n",dV,get_torsion_value( a1,a2,a3,a4 ));
	}
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer( vfs,ffs,nfs, 100,30, (1==0));
//		block_optimizer( vfs,ffs,nfs, 10,30, (1==0));
//standard		block_optimizer( vfs,ffs,nfs, 10,a_number()/20, (1==0));
		pc_cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
//WORKED		cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
	dV = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&dV,0.);
	if( isnan(dV) )
	{
	fprintf(op,"error condition - restoring from previous best\n");
		restore_from_clone(1);
	mxdq = get_f_variable("mxdq");
	set_f_variable("mxdq", -1.);
		V = 0.;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
		dV = V;
	set_f_variable("mxdq", mxdq);
	} else{
	dx = exp(  (V-dV)/temper*1.987);
	fprintf(op, "V %f, dV, %f, Delta %f %f\n", V,dV, -V+dV,dx);
	if( dV < Vbest){ clone(1); Vbest = dV;}
	if( dV < V){ clone(2);}
	else {
		if( dx < randf(-1)){
			fprintf(op,"rejected\n"); 
			restore_from_clone(2);
			mxdq = get_f_variable("mxdq");
			set_f_variable("mxdq", -1.);
				V = 0.;
				for( i=0; i< nfs; i++)
					(*vfs[i])(&V,0.);
			set_f_variable("mxdq", mxdq);
			}
	}
	}
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
	}/* istep */

	restore_from_clone(1);
	free( aps);
}
