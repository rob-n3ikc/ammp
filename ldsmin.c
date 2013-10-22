// two changes over time
// 1) dV is now a sum of 10 small rng
// 2) block_opt() taken from 10*depth,30 to 10,100
//
/*
*  possible improvment on Monte-Carlo minimization
*
*  calculate local dominating sets on 'cb' - 'cb' distance
*  then vary only torsions in the lds
*
*  an lds is a connected dominating set 
* (for the distance graph and given cutoff)
*  however, the depth of the tree is restricted to (one, two, three) a 
*  small integer so that completely connected sets are not tested 
*  (simply because if everything is connected then we're back to the mcmin
*   algorithm with a sub-optimal step choice)
*
*/
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

//#define CUTOFF 5.
#define CUTOFF 5.
/* clone routines from clone.c */

int clone( int);
int restore_from_clone(int);
int tset_bond_build( ATOM*, ATOM* , ATOM*, ATOM*);
//int set_torsion( ATOM* , ATOM*, ATOM*, ATOM*, float);
float get_torsion_value( ATOM* , ATOM*, ATOM*, ATOM*);
int math_match_atom( char*, ATOM*);
void get_bond( ATOM*, ATOM *[], int, int*);
float get_f_variable(char *);
int ldsmin( op,vfs,ffs, nfs, at1,at2, temper, nstep)
	FILE *op;
	int (*vfs[])(), (*ffs[])(), nfs, nstep;
	float temper;
	char *at1,*at2;
{
	ATOM *(*aps)[], *a_next();
	ATOM *a1,*a2,*a3,*a4;
	ATOM *bonded[10];
#define MAXSET 100
	ATOM *ld_list[MAXSET];
	int in_the_list;
	int a_number(), numatom;
	int i,istep,ilayer;
	float randf();
	float V,Vbest, dV,dt;	
	int in_d_matrix;
	float *d_matrix;
        float mxdq;
	ATOM **atoms_in_matrix;
	
	int math_match_atom(); /* math_match_atom( char*, ATOM*)  defined in math.c*/

	numatom = a_number();
	if( numatom < 4) return 0;
	aps = malloc( numatom * sizeof( ATOM *));
	if( aps == NULL) { aaerror("LDSMIN error MEMORY exhausted\n"); exit(0);}
	in_d_matrix = 0;
	for(i = 0; i< numatom; i++)
	{
		(*aps)[i] = a_next(i);
		if( math_match_atom( at2, (*aps)[i])  != 0)
		{
			in_d_matrix++;
		}
	}
	d_matrix = (float *)malloc( in_d_matrix*in_d_matrix*sizeof(float));
	atoms_in_matrix = (ATOM **)malloc( in_d_matrix*sizeof(ATOM *));
	in_d_matrix = 0;
	for( i=0; i< numatom; i++)
	{
		if( math_match_atom( at2, (*aps)[i])  != 0)
		{
			atoms_in_matrix[in_d_matrix] = (*aps)[i];
			in_d_matrix++;
		}
	}
	dt = 3.1415926535*0.5;
	dt = 3.1415926535;
	dt = dt*0.1;
	// special version for 2pi/3 torsion shifts;
//	dt = 3.1415926535*2./3.;

	clone( 1);  /* arbitrary integer label */
	clone( 2);  /* arbitrary integer label */
//
	V = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&V,0.);
	Vbest = V;

	for( istep=0; istep< nstep; istep++)
	{
		int in_bond;
		int j,k;
		int listed[MAXSET]; /* use this to keep track of what is put in ld_list */
		float dx;
	V = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&V,0.);
/* make the distance matrix */
	for( i=0; i< in_d_matrix*in_d_matrix; i++) d_matrix[i] = 0.;
	for( i=0; i< in_d_matrix; i++)
	for( j= 0; j < in_d_matrix; j++)
	{
		float xx,yy,zz,r;
		xx = atoms_in_matrix[i]->x - atoms_in_matrix[j]->x;
		yy = atoms_in_matrix[i]->y - atoms_in_matrix[j]->y;
		zz = atoms_in_matrix[i]->z - atoms_in_matrix[j]->z;
		r = xx*xx +yy*yy + zz*zz;
		d_matrix[ i*in_d_matrix + j] = sqrt(r);
	}
/* now pick an atom from the list and build the lds list */
	k = in_d_matrix*randf();
	if( k < 0) k = 0; 
	if( k >= in_d_matrix) k = in_d_matrix -1;

/* make sure you have the center atom */
	in_the_list = 1;
	ld_list[0] = atoms_in_matrix[k];
	listed[0] = k;
	 i = k*in_d_matrix;
	for( j=0; j< in_d_matrix; j++)
	{
		if( j == k ) continue;
		if( d_matrix[i+j] < CUTOFF){ ld_list[in_the_list] = atoms_in_matrix[j]; listed[in_the_list++] = j;}
		if( in_the_list == MAXSET) break;
	}
/* second, 3,4 layer of lds */
	for( ilayer=2; ilayer < 7; ilayer++)
	if( in_the_list < MAXSET){

	in_bond = in_the_list;
	for( k=0; k< in_bond; k++)
	{
	i = listed[k] *in_d_matrix;
	for( j=0; j< in_d_matrix; j++)
	{
		int ll;
//		if( j == listed[k] ) continue;
		for(ll = 0; ll < in_bond; ll++)
			if( j == listed[ll]) goto SKIP; 
		if( d_matrix[i+j] < CUTOFF){ ld_list[in_the_list] = atoms_in_matrix[j]; listed[in_the_list++] = j;}
		if( in_the_list == MAXSET) break;
SKIP: ;
	}
	}
	}/* second layer */
	for( j = 0; j< 1; j++)
	{
	int list_index;
	clone( 2);  /* arbitrary integer label */
		for( list_index = 0; list_index < in_the_list; list_index++)
		{
		a3 = ld_list[list_index];
/* find those bonded to me who are of type at1 */
		get_bond( a3,bonded, 10, &in_bond);
		a2 = NULL;
		for( i =0; i< in_bond; i++)
			if( math_match_atom( at1,bonded[i]) != 0)
			{a2 = bonded[i]; break;}
		if( a2 == NULL) goto RE_CHOOSE;
/* now find those bonded to my mate who are not of type at2 */
		get_bond( a2,bonded, 10, &in_bond);
		a1 = NULL;
		for( i=0; i< in_bond; i++)
			if( math_match_atom( at2,bonded[i]) == 0)
			{a1 = bonded[i]; break;}
		if( a1 == NULL) goto RE_CHOOSE;
/* find those bound to me mate who are not of type at1 (my type) */
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
/*
	dV = dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
*/
	dV = dt*(2.*randf()-1.);
//	for( i=0; i< 10;i++)
	for( i=0; i< 20;i++)
	dV += dt*(2.*randf()-1.);
/* not quite right
//	dV = 3.1415926535*2./3; if (randf() < 0.5) dV = -dV;
*/
/*
this is the way to do discrete steps
	if( in_the_list > 1)
	{
	i = (int)(3. *randf());
	dV = 3.1415926535*2./3.*i;
	}
	else
	{ dV = 3.1415926535*2./3; if (randf() < 0.5) dV = -dV; }
*/
/* skip CSS (RARE PROTEIN SPECIFIC CODE- THESE ARE DISULFIDES) */
/* skip PRO  because torsion never works on ca-cb */
	if( strcmp(a2->name, "css.ca") != 0)
	if( strcmp(a2->name, "pro.ca") != 0)
	{
	fflush(stdout);
		tset_bond_build( a1,a2,a3,a4);
/* set_torsion uses (or should use a) delta_value */
		set_torsion( a1,a2,a3,a4, dV );
	}
	fprintf(op,"%f %f \n",dV,get_torsion_value( a1,a2,a3,a4 ));
RE_CHOOSE: ;
		}// in the list
	}
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer(vfs,ffs,nfs,10*in_the_list,30, (1==0));
		block_optimizer(vfs,ffs,nfs,10,100, (1==0));
// CASP7 USED THIS		block_optimizer(vfs,ffs,nfs,10,100, (1==0));
//		pc_cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
//		block_optimizer(vfs,ffs,nfs,100,30, (1==0));
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer(vfs,ffs,nfs,10,a_number()/20, (1==0));
//orig		cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
	dV = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&dV,0.);
	dx = exp(  (V-dV)/temper*1.987);
	fprintf(op, "V %f, dV, %f, Delta %f %f\n", V,dV, -V+dV,dx);
	if( isnan(dV) )
	{
		fprintf(op,"ERROR CONDITION -restore from best\n");
		restore_from_clone(1);
		mxdq = get_f_variable("mxdq");
		set_f_variable("mxdq",-1.);
		V = 0;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
		dV = V;
		set_f_variable("mxdq",mxdq);
	}
	if( dV < Vbest){ clone(1); Vbest = dV;}
	if( dV < V){ 
			clone(2);
		mxdq = get_f_variable("mxdq");
		set_f_variable("mxdq",-1.);
		V = 0;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
		set_f_variable("mxdq",mxdq);
			}
	else {
		if( dx < randf(-1)) {
		fprintf(op,"REJECTED\n");
		restore_from_clone(2);
		mxdq = get_f_variable("mxdq");
		set_f_variable("mxdq",-1.);
		V = 0;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
		set_f_variable("mxdq",mxdq);
		}
	}
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
	}/* istep */

	restore_from_clone(1);
	free( atoms_in_matrix);
	free( d_matrix);
	free( aps);
}
/* restricted version uses inactive/active to restrict search */
int restricted_ldsmin( op,vfs,ffs, nfs, at1,at2, temper, nstep)
	FILE *op;
	int (*vfs[])(), (*ffs[])(), nfs, nstep;
	float temper;
	char *at1,*at2;
{
	ATOM *(*aps)[], *a_next();
	ATOM *a1,*a2,*a3,*a4;
	ATOM *bonded[10];
#define MAXSET 100
	ATOM *ld_list[MAXSET];
	int in_the_list;
	int a_number(), numatom;
	int i,istep,ilayer;
	float randf();
	float V,Vbest, dV,dt;	
	int in_d_matrix;
	float *d_matrix;
	ATOM **atoms_in_matrix;
	
	int math_match_atom(); /* math_match_atom( char*, ATOM*)  defined in math.c*/

	numatom = a_number();
	if( numatom < 4) return 0;
	aps = malloc( numatom * sizeof( ATOM *));
	if( aps == NULL) { aaerror("LDSMIN error MEMORY exhausted\n"); exit(0);}
	in_d_matrix = 0;
	for(i = 0; i< numatom; i++)
	{
		(*aps)[i] = a_next(i);
		if( math_match_atom( at2, (*aps)[i])  != 0)
		{
			in_d_matrix++;
		}
	}
	d_matrix = (float *)malloc( in_d_matrix*in_d_matrix*sizeof(float));
	atoms_in_matrix = (ATOM **)malloc( in_d_matrix*sizeof(ATOM *));
	in_d_matrix = 0;
	for( i=0; i< numatom; i++)
	{
		if( math_match_atom( at2, (*aps)[i])  != 0)
		{
			atoms_in_matrix[in_d_matrix] = (*aps)[i];
			in_d_matrix++;
		}
	}
	dt = 3.1415926535*0.5;
	dt = 3.1415926535;
	dt = dt*0.1;
	// special version for 2pi/3 torsion shifts;
//	dt = 3.1415926535*2./3.;

	clone( 1);  /* arbitrary integer label */
	Vbest = 10.e10;
//

	for( istep=0; istep< nstep; istep++)
	{
		int in_bond;
		int j,k;
		int listed[MAXSET]; /* use this to keep track of what is put in ld_list */
		float dx;
	activate( 0, 10000000);
	V = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&V,0.);
	if( Vbest > V)Vbest = V;
	inactivate( 0, 10000000);
/* make the distance matrix */
	for( i=0; i< in_d_matrix*in_d_matrix; i++) d_matrix[i] = 0.;
	for( i=0; i< in_d_matrix; i++)
	for( j= 0; j < in_d_matrix; j++)
	{
		float xx,yy,zz,r;
		xx = atoms_in_matrix[i]->x - atoms_in_matrix[j]->x;
		yy = atoms_in_matrix[i]->y - atoms_in_matrix[j]->y;
		zz = atoms_in_matrix[i]->z - atoms_in_matrix[j]->z;
		r = xx*xx +yy*yy + zz*zz;
		d_matrix[ i*in_d_matrix + j] = sqrt(r);
	}
/* now pick an atom from the list and build the lds list */
	k = in_d_matrix*randf();
	if( k < 0) k = 0; 
	if( k >= in_d_matrix) k = in_d_matrix -1;

/* make sure you have the center atom */
	in_the_list = 1;
	ld_list[0] = atoms_in_matrix[k];
	listed[0] = k;
	 i = k*in_d_matrix;
	for( j=0; j< in_d_matrix; j++)
	{
		if( j == k ) continue;
		if( d_matrix[i+j] < CUTOFF){ ld_list[in_the_list] = atoms_in_matrix[j]; listed[in_the_list++] = j;}
		if( in_the_list == MAXSET) break;
	}
/* second, 3,4 layer of lds */
	for( ilayer=2; ilayer < 7; ilayer++)
	if( in_the_list < MAXSET){

	in_bond = in_the_list;
	for( k=0; k< in_bond; k++)
	{
	i = listed[k] *in_d_matrix;
	for( j=0; j< in_d_matrix; j++)
	{
		int ll;
//		if( j == listed[k] ) continue;
		for(ll = 0; ll < in_bond; ll++)
			if( j == listed[ll]) goto SKIP; 
		if( d_matrix[i+j] < CUTOFF){ ld_list[in_the_list] = atoms_in_matrix[j]; listed[in_the_list++] = j;}
		if( in_the_list == MAXSET) break;
SKIP: ;
	}
	}
	}/* second layer */
	for( j = 0; j< 1; j++)
	{
	int list_index;
	clone( 2);  /* arbitrary integer label */
		for( list_index = 0; list_index < in_the_list; list_index++)
		{
		a3 = ld_list[list_index];
		int ilow;
		ilow = a3->serial;
		ilow = ilow - ilow%100;
		activate( ilow, ilow+99);
/* find those bonded to me who are of type at1 */
		get_bond( a3,bonded, 10, &in_bond);
		a2 = NULL;
		for( i =0; i< in_bond; i++)
			if( math_match_atom( at1,bonded[i]) != 0)
			{a2 = bonded[i]; break;}
		if( a2 == NULL) goto RE_CHOOSE;
/* now find those bonded to my mate who are not of type at2 */
		get_bond( a2,bonded, 10, &in_bond);
		a1 = NULL;
		for( i=0; i< in_bond; i++)
			if( math_match_atom( at2,bonded[i]) == 0)
			{a1 = bonded[i]; break;}
		if( a1 == NULL) goto RE_CHOOSE;
/* find those bound to me mate who are not of type at1 (my type) */
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
/*
	dV = dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
*/
	dV = dt*(2.*randf()-1.);
//	for( i=0; i< 10;i++)
	for( i=0; i< 20;i++)
	dV += dt*(2.*randf()-1.);
/* not quite right
//	dV = 3.1415926535*2./3; if (randf() < 0.5) dV = -dV;
*/
/*
this is the way to do discrete steps
	if( in_the_list > 1)
	{
	i = (int)(3. *randf());
	dV = 3.1415926535*2./3.*i;
	}
	else
	{ dV = 3.1415926535*2./3; if (randf() < 0.5) dV = -dV; }
*/
/* skip CSS (RARE PROTEIN SPECIFIC CODE- THESE ARE DISULFIDES) */
/* skip PRO  because torsion never works on ca-cb */
	if( strcmp(a2->name, "css.ca") != 0)
	if( strcmp(a2->name, "pro.ca") != 0)
	{
	fflush(stdout);
		tset_bond_build( a1,a2,a3,a4);
/* set_torsion uses (or should use a) delta_value */
		set_torsion( a1,a2,a3,a4, dV );
	}
	fprintf(op,"%f %f \n",dV,get_torsion_value( a1,a2,a3,a4 ));
RE_CHOOSE: ;
		}// in the list
	}
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer(vfs,ffs,nfs,10*in_the_list,30, (1==0));
		block_optimizer(vfs,ffs,nfs,10,100, (1==0));
//		block_optimizer(vfs,ffs,nfs,100,30, (1==0));
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer(vfs,ffs,nfs,10,a_number()/20, (1==0));
//orig		cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
	dV = 0;
	activate( 0, 10000000);
	for( i=0; i< nfs; i++)
		(*vfs[i])(&dV,0.);
	dx = exp(  (V-dV)/temper*1.987);
	fprintf(op, "V %f, dV, %f, Delta %f %f\n", V,dV, -V+dV,dx);
	if( isnan(dV) )
	{
		fprintf(op,"ERROR CONDITION -restore from best\n");
		restore_from_clone(1);
		V = 0;
		activate( 0, 10000000);
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
		dV = V;
	}
	if( dV < Vbest){ clone(1); Vbest = dV;}
	if( dV < V){ clone(2);}
	else {
		if( dx < randf(-1)) {fprintf(op,"REJECTED\n");restore_from_clone(2);}
	}
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
	}/* istep */

	restore_from_clone(1);
	free( atoms_in_matrix);
	free( d_matrix);
	free( aps);
}


//
// use multiple search (i.e. tgroup/tsearch) on the l.d.set
//

// two changes over time
// 1) dV is now a sum of 10 small rng
// 2) block_opt() taken from 10*depth,30 to 10,100
//
/*
*  possible improvment on Monte-Carlo minimization
*
*  calculate local dominating sets on 'cb' - 'cb' distance
*  then vary only torsions in the lds
*
*  an lds is a connected dominating set 
* (for the distance graph and given cutoff)
*  however, the depth of the tree is restricted to (one, two, three) a 
*  small integer so that completely connected sets are not tested 
*  (simply because if everything is connected then we're back to the mcmin
*   algorithm with a sub-optimal step choice)
*
*/
/* implementation of Monte-Carlo /minimization
 * i.e. score for simulated annealling after minimization
 *  depends on the clone routines
 *  
* also depends on tset  actually using tset_bond_build and set_torsion
 */
/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ammp.h"

#ifdef graphic
void force_screen_update();
#endif

*/
//#define CUTOFF 5.
#define CUTOFF 5.
/* clone routines from clone.c */

int clone( int);
int restore_from_clone(int);
int tset_bond_build( ATOM*, ATOM* , ATOM*, ATOM*);
//int set_torsion( ATOM* , ATOM*, ATOM*, ATOM*, float);
float get_torsion_value( ATOM* , ATOM*, ATOM*, ATOM*);
int math_match_atom( char*, ATOM*);
void get_bond( ATOM*, ATOM *[], int, int*);
int ldssearch( op,vfs,ffs, nfs, at1,at2, temper, nstep)
	FILE *op;
	int (*vfs[])(), (*ffs[])(), nfs, nstep;
	float temper;
	char *at1,*at2;
{
	ATOM *(*aps)[], *a_next();
	ATOM *a1,*a2,*a3,*a4;
	ATOM *bonded[10];
#define MAXSET 100
	ATOM *ld_list[MAXSET];
	int in_the_list;
	int a_number(), numatom;
	int i,istep,ilayer;
	float randf();
	float V,Vbest, dV,dt;	
	int in_d_matrix;
	float *d_matrix;
	ATOM **atoms_in_matrix;
	
	int math_match_atom(); /* math_match_atom( char*, ATOM*)  defined in math.c*/

	numatom = a_number();
	if( numatom < 4) return 0;
	aps = malloc( numatom * sizeof( ATOM *));
	if( aps == NULL) { aaerror("LDSMIN error MEMORY exhausted\n"); exit(0);}
	in_d_matrix = 0;
	for(i = 0; i< numatom; i++)
	{
		(*aps)[i] = a_next(i);
		if( math_match_atom( at2, (*aps)[i])  != 0)
		{
			in_d_matrix++;
		}
	}
	d_matrix = (float *)malloc( in_d_matrix*in_d_matrix*sizeof(float));
	atoms_in_matrix = (ATOM **)malloc( in_d_matrix*sizeof(ATOM *));
	in_d_matrix = 0;
	for( i=0; i< numatom; i++)
	{
		if( math_match_atom( at2, (*aps)[i])  != 0)
		{
			atoms_in_matrix[in_d_matrix] = (*aps)[i];
			in_d_matrix++;
		}
	}
	dt = 3.1415926535*0.5;
	dt = 3.1415926535;
	dt = dt*0.1;
	// special version for 2pi/3 torsion shifts;
//	dt = 3.1415926535*2./3.;

	clone( 1);  /* arbitrary integer label */
//

	for( istep=0; istep< nstep; istep++)
	{
		int in_bond;
		int j,k;
		int listed[MAXSET]; /* use this to keep track of what is put in ld_list */
		float dx;
	V = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&V,0.);
	Vbest = V;
/* make the distance matrix */
	for( i=0; i< in_d_matrix*in_d_matrix; i++) d_matrix[i] = 0.;
	for( i=0; i< in_d_matrix; i++)
	for( j= 0; j < in_d_matrix; j++)
	{
		float xx,yy,zz,r;
		xx = atoms_in_matrix[i]->x - atoms_in_matrix[j]->x;
		yy = atoms_in_matrix[i]->y - atoms_in_matrix[j]->y;
		zz = atoms_in_matrix[i]->z - atoms_in_matrix[j]->z;
		r = xx*xx +yy*yy + zz*zz;
		d_matrix[ i*in_d_matrix + j] = sqrt(r);
	}
/* now pick an atom from the list and build the lds list */
	k = in_d_matrix*randf();
	if( k < 0) k = 0; 
	if( k >= in_d_matrix) k = in_d_matrix -1;

/* make sure you have the center atom */
	in_the_list = 1;
	ld_list[0] = atoms_in_matrix[k];
	listed[0] = k;
	 i = k*in_d_matrix;
	for( j=0; j< in_d_matrix; j++)
	{
		if( j == k ) continue;
		if( d_matrix[i+j] < CUTOFF){ ld_list[in_the_list] = atoms_in_matrix[j]; listed[in_the_list++] = j;}
		if( in_the_list == MAXSET) break;
	}
/* second, 3,4 layer of lds */
	for( ilayer=2; ilayer < 7; ilayer++)
	if( in_the_list < MAXSET){

	in_bond = in_the_list;
	for( k=0; k< in_bond; k++)
	{
	i = listed[k] *in_d_matrix;
	for( j=0; j< in_d_matrix; j++)
	{
		int ll;
//		if( j == listed[k] ) continue;
		for(ll = 0; ll < in_bond; ll++)
			if( j == listed[ll]) goto SKIP; 
		if( d_matrix[i+j] < CUTOFF){ ld_list[in_the_list] = atoms_in_matrix[j]; listed[in_the_list++] = j;}
		if( in_the_list == MAXSET) break;
SKIP: ;
	}
	}
	}/* second layer */
	for( j = 0; j< 1; j++)
	{
	int list_index;
	clone( 2);  /* arbitrary integer label */
		for( list_index = 0; list_index < in_the_list; list_index++)
		{
		a3 = ld_list[list_index];
/* find those bonded to me who are of type at1 */
		get_bond( a3,bonded, 10, &in_bond);
		a2 = NULL;
		for( i =0; i< in_bond; i++)
			if( math_match_atom( at1,bonded[i]) != 0)
			{a2 = bonded[i]; break;}
		if( a2 == NULL) goto RE_CHOOSE;
/* now find those bonded to my mate who are not of type at2 */
		get_bond( a2,bonded, 10, &in_bond);
		a1 = NULL;
		for( i=0; i< in_bond; i++)
			if( math_match_atom( at2,bonded[i]) == 0)
			{a1 = bonded[i]; break;}
		if( a1 == NULL) goto RE_CHOOSE;
/* find those bound to me mate who are not of type at1 (my type) */
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
/*
	dV = dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
	dV += dt*(2.*randf()-1.);
*/
	dV = dt*(2.*randf()-1.);
//	for( i=0; i< 10;i++)
	for( i=0; i< 20;i++)
	dV += dt*(2.*randf()-1.);
/* not quite right
//	dV = 3.1415926535*2./3; if (randf() < 0.5) dV = -dV;
*/
/*
this is the way to do discrete steps
	if( in_the_list > 1)
	{
	i = (int)(3. *randf());
	dV = 3.1415926535*2./3.*i;
	}
	else
	{ dV = 3.1415926535*2./3; if (randf() < 0.5) dV = -dV; }
*/
/* skip CSS (RARE PROTEIN SPECIFIC CODE- THESE ARE DISULFIDES) */
/* skip PRO  because torsion never works on ca-cb */
	if( strcmp(a2->name, "css.ca") != 0)
	if( strcmp(a2->name, "pro.ca") != 0)
	{
	fflush(stdout);
/* set_torsion uses (or should use a) delta_value */
/*		tset_bond_build( a1,a2,a3,a4);
		set_torsion( a1,a2,a3,a4, dV );
*/
		tmin( op,1==1,a1->serial,a2->serial,a3->serial,a4->serial, 6, vfs,nfs);
	}
	fprintf(op,"%f %f \n",dV,get_torsion_value( a1,a2,a3,a4 ));
RE_CHOOSE: ;
		}// in the list
	}
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer(vfs,ffs,nfs,10*in_the_list,30, (1==0));
		block_optimizer(vfs,ffs,nfs,10,100, (1==0));
// CASP7 USED THIS		block_optimizer(vfs,ffs,nfs,10,100, (1==0));
//		pc_cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
//		block_optimizer(vfs,ffs,nfs,100,30, (1==0));
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
//		block_optimizer(vfs,ffs,nfs,10,a_number()/20, (1==0));
//orig		cngdel( vfs,ffs,nfs, 10,10,0., (1==0));
	dV = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&dV,0.);
	dx = exp(  (V-dV)/temper*1.987);
	fprintf(op, "V %f, dV, %f, Delta %f %f\n", V,dV, -V+dV,dx);
	if( isnan(dV) )
	{
		fprintf(op,"ERROR CONDITION -restore from best\n");
		restore_from_clone(1);
		V = 0;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
		dV = V;
	}
	if( dV < Vbest){ clone(1); Vbest = dV;}
	if( dV < V){ clone(2);}
	else {
		if( dx < randf(-1)) {fprintf(op,"REJECTED\n");restore_from_clone(2);}
	}
#ifdef graphic
 force_screen_update();
#endif
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
	}/* istep */

	restore_from_clone(1);
	free( atoms_in_matrix);
	free( d_matrix);
	free( aps);
}
