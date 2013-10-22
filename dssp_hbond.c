/* this program implements the dssp hbond potential as used in tasser
*
* it needs the tuple  constant (q1,q2) names of four atoms;
*
* it then calculates
*  332.17752 * (q1q2) *( 1/rON + 1/rCH - 1/rOH -1/rCN)
*  the potential and derivative are 0 if rON > 5.2
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "ammp.h"

int math_match_atom( char*, ATOM* );
void get_bond(ATOM*, ATOM *[], int, int*);
void a_inactive_f_zero(void);

typedef struct{
	float k;
	char O[10];
	char N[10];
	char C[10];
	char H[10];
	void *next;
}  DSSP_HBOND;

DSSP_HBOND *dssp_first = NULL;

//int dssp_hbond( float k, char *o, char *n, char *c, char *h)
int dssp_hbond(  k, o,n,c,h)
float k;
char *o, *n, *c,*h;
{
	DSSP_HBOND *d,*dlast;
	d = dssp_first;
	dlast = dssp_first;
	while( d != NULL)
	{dlast = d;  d = d->next;}
	
	d = malloc( sizeof( DSSP_HBOND));
	if( d == NULL) return 1;

	d->k = k;
	strncpy( d->O,o,9);
	strncpy( d->N,n,9);
	strncpy( d->C,c,9);
	strncpy( d->H,h,9);
	d->next = NULL;

	if( dlast == NULL)
	{ dssp_first = d; }
	else
	{  dlast->next = d;}
}

//int v_dssp_hbond( float *V, float lambda)
int v_dssp_hbond( V,lambda)
float *V;
float lambda;
{
	DSSP_HBOND *dp;
	float ron,rcn,roh,rch;
	ATOM *O,*N,*C,*H;
	ATOM *bonded[10];
	ATOM *a1,*a2;
	int in_bond;
	ATOM *a_next(int);
	float r,x,y,z,x1,y1,z1,dx,dy,dz;
	int numatom,i,j,k, a_number(void);

	numatom = a_number();
	if( numatom < 4) return 0;
	if( dssp_first == NULL) return 0;
	dp = dssp_first;
DOIT_AGAIN:;	
	a1 = a_next(-1);
	for(  i=0; i<numatom; i++)
	{
	if( math_match_atom(dp->N,a1) == 0 ){a1 = a1->next; continue;}
	x = a1->x + lambda *a1->dx;
	y = a1->y + lambda *a1->dy;
	z = a1->z + lambda *a1->dz;
	N = a1;
	H = NULL;
	get_bond( N,bonded,10,&in_bond);
	for( j = 0; j< in_bond; j++)
		if( math_match_atom( dp->H,bonded[j]) !=0)
		{  H = bonded[j]; break;}
	if( H == NULL){a1= a1->next; continue;}	
	a2 = a_next(-1);
	for(  j = 0 ;j < numatom; j++)
	{
	if( a1 == a2) { a2 = a2->next; continue;}
	if( a1->serial/100 == a2->serial/100) { a2 = a2->next; continue;}
	if( math_match_atom(dp->O,a2) == 0 ){a2 = a2->next; continue;}
	x1 = a2->x + lambda *a2->dx;
	y1 = a2->y + lambda *a2->dy;
	z1 = a2->z + lambda *a2->dz;
	dx = x-x1;
	dy = y-y1;
	dz = z-z1;
	r = dx*dx+dy*dy+dz*dz;
// when commented out it is continuous, but SLOW
	if( r > 5.2*5.2) { a2= a2->next; continue;}
// ok so now we have a candidate 

	O = a2;
	C = NULL;
	get_bond( O,bonded,10,&in_bond);
	for( k = 0; k< in_bond; k++)
		if( math_match_atom( dp->C,bonded[k]) !=0)
		{  C = bonded[k]; break;}
	if( C == NULL)  {a2 = a2->next; continue;}
	ron = sqrt(r);

	dx = O->x -H->x + lambda*( O->dx - H->dx);
	dy = O->y -H->y + lambda*( O->dy - H->dy);
	dz = O->z -H->z + lambda*( O->dz - H->dz);
	roh =sqrt( dx*dx+dy*dy+dz*dz);

	dx = C->x -H->x + lambda*( C->dx - H->dx);
	dy = C->y -H->y + lambda*( C->dy - H->dy);
	dz = C->z -H->z + lambda*( C->dz - H->dz);
	rch =sqrt( dx*dx+dy*dy+dz*dz);

	dx = C->x -N->x + lambda*( C->dx - N->dx);
	dy = C->y -N->y + lambda*( C->dy - N->dy);
	dz = C->z -N->z + lambda*( C->dz - N->dz);
	rcn =sqrt( dx*dx+dy*dy+dz*dz);
	
//printf("%s %s %s %s %f %f %f %f\n", O->name, N->name, C->name, H->name, ron,rcn,rch,roh);

// chosen to be repelling if unphysical
	if( ron < 0.5) ron = 0.5;
	if( rch < 0.5) rch = 0.5;
	if( rcn < 1.) rcn  = 1.;
	if( roh < 1.) roh  = 1.;
//*  332.17752 * (q1q2) *( 1/rON + 1/rCH - 1/rOH -1/rCN)
	*V+= 332.17752*dp->k*(1./ron + 1./rch -1./roh - 1./rcn);


	a2 = a2->next;
	}//j
	a1 = a1->next;
	}//i 
	
	if( dp->next == NULL) return 0;
	goto DOIT_AGAIN;
}


int f_dssp_hbond(  lambda)
float lambda;
{
	DSSP_HBOND *dp;
	float ron,rcn,roh,rch;
	float V;
	ATOM *O,*N,*C,*H;
	ATOM *bonded[10];
	ATOM *a1,*a2;
	int in_bond;
	ATOM *a_next(int);
	float r,x,y,z,x1,y1,z1,dx,dy,dz;
	float d_on_x,d_on_y,d_on_z;
	float d_cn_x,d_cn_y,d_cn_z;
	float d_ch_x,d_ch_y,d_ch_z;
	float d_oh_x,d_oh_y,d_oh_z;
	int numatom,i,j,k, a_number(void);

	numatom = a_number();
	if( numatom < 4) return 0;
	if( dssp_first == NULL) return 0;
	dp = dssp_first;
DOIT_AGAIN:;	
	a1 = a_next(-1);
	for(  i=0; i<numatom; i++)
	{
	if( math_match_atom(dp->N,a1) == 0 ){a1 = a1->next; continue;}
	x = a1->x + lambda *a1->dx;
	y = a1->y + lambda *a1->dy;
	z = a1->z + lambda *a1->dz;
	N = a1;
	H = NULL;
	get_bond( N,bonded,10,&in_bond);
	for( j = 0; j< in_bond; j++)
		if( math_match_atom( dp->H,bonded[j]) !=0)
		{  H = bonded[j]; break;}
	if( H == NULL){a1= a1->next; continue;}	
	a2 = a_next(-1);
	for(  j = 0 ;j < numatom; j++)
	{
	if( a1 == a2) { a2 = a2->next; continue;}
	if( a1->serial/100 == a2->serial/100) { a2 = a2->next; continue;}
	if( math_match_atom(dp->O,a2) == 0 ){a2 = a2->next; continue;}
	x1 = a2->x + lambda *a2->dx;
	y1 = a2->y + lambda *a2->dy;
	z1 = a2->z + lambda *a2->dz;
	dx = x-x1;
	dy = y-y1;
	dz = z-z1;
	
	r = dx*dx+dy*dy+dz*dz;
// when commented out it is continous, but SLOW
	if( r > 5.2*5.2) { a2= a2->next; continue;}
// ok so now we have a candidate 

	O = a2;
	C = NULL;
	get_bond( O,bonded,10,&in_bond);
	for( k = 0; k< in_bond; k++)
		if( math_match_atom( dp->C,bonded[k]) !=0)
		{  C = bonded[k]; break;}
	if( C == NULL)  {a2 = a2->next; continue;}
	ron = sqrt(r);
	d_on_x = -dx/(r*ron);
	d_on_y = -dy/(r*ron);
	d_on_z = -dz/(r*ron);

	dx = O->x -H->x + lambda*( O->dx - H->dx);
	dy = O->y -H->y + lambda*( O->dy - H->dy);
	dz = O->z -H->z + lambda*( O->dz - H->dz);
	r = dx*dx+dy*dy+dz*dz;
	roh =sqrt( r);
	d_oh_x = dx/(r*roh);
	d_oh_y = dy/(r*roh);
	d_oh_z = dz/(r*roh);

	dx = C->x -H->x + lambda*( C->dx - H->dx);
	dy = C->y -H->y + lambda*( C->dy - H->dy);
	dz = C->z -H->z + lambda*( C->dz - H->dz);
	r = dx*dx+dy*dy+dz*dz;
	rch =sqrt( r);
	d_ch_x = dx/(r*rch);
	d_ch_y = dy/(r*rch);
	d_ch_z = dz/(r*rch);

	dx = C->x -N->x + lambda*( C->dx - N->dx);
	dy = C->y -N->y + lambda*( C->dy - N->dy);
	dz = C->z -N->z + lambda*( C->dz - N->dz);
	r = dx*dx+dy*dy+dz*dz;
	rcn =sqrt( r);
	d_cn_x = dx/(r*rcn);
	d_cn_y = dy/(r*rcn);
	d_cn_z = dz/(r*rcn);
	
// chosen to be repelling if unphysical
//	if( ron < 0.5) ron = 0.5;
//	if( rch < 0.5) rch = 0.5;
//	if( rcn < 1.) rcn  = 1.;
//	if( roh < 1.) roh  = 1.;
//*  332.17752 * (q1q2) *( 1/rON + 1/rCH - 1/rOH -1/rCN)
//	V = 332.17752*dp->k*(1./ron + 1./rch -1./roh - 1./rcn);
	V =    332.11752*dp->k;;
	O->fx += V*(d_on_x - d_oh_x);
	O->fy += V*(d_on_y - d_oh_y);
	O->fz += V*(d_on_z - d_oh_z);
	N->fx += V*(-d_on_x + d_cn_x);
	N->fy += V*(-d_on_y + d_cn_y);
	N->fz += V*(-d_on_z + d_cn_z);
	C->fx += V*(d_ch_x  - d_cn_x);
	C->fy += V*(d_ch_y  - d_cn_y);
	C->fz += V*(d_ch_z  - d_cn_z);
	H->fx += V*(-d_ch_x + d_oh_x);
	H->fy += V*(-d_ch_y + d_oh_y);
	H->fz += V*(-d_ch_z + d_oh_z);
	

	a2 = a2->next;
	}//j
	a1 = a1->next;
	}//i 
	
	if( dp->next == NULL) { a_inactive_f_zero();return 0;}
	goto DOIT_AGAIN;
}


int dump_dssp_hbond( where)
FILE *where;
{
	DSSP_HBOND *d;
	d = dssp_first;
	while( d != NULL)
	{
	fprintf( where,"dssp_hbond %f %s %s %s %s ;\n",
		d->k, d->O, d->N, d->C, d->H);
	  d = d->next;}

	return 0;
}

int gsdg_dssp_hbond( who)
ATOM *who;
{
	DSSP_HBOND *dp;
	float ron,rcn,roh,rch;
	ATOM *O,*N,*C,*H;
	ATOM *bonded[10];
	ATOM *a1,*a2;
	int in_bond;
	ATOM *a_next(int);
	float r,x,y,z,x1,y1,z1,dx,dy,dz;
	int numatom,i,j,k, a_number(void);

	numatom = a_number();
	if( numatom < 4) return 0;

	dp = dssp_first;
	while( dp != NULL)
	{
		if( math_match_atom( dp->N,who) != 0)
		{
		N = who;
		get_bond( N,bonded,10,&in_bond);
		H = NULL;
		for( j = 0; j< in_bond; j++)
			if( math_match_atom( dp->H,bonded[j]) !=0)
			{  H = bonded[j]; break;}
		if( H == NULL){return 0;}
// so now find the O and check the distance;
		a2 = a_next(-1);
		for( j=0; j< numatom; j++)
		{
	if( who == a2) { a2 = a2->next; continue;}
	if( who->serial/100 == a2->serial/100) { a2 = a2->next; continue;}
	if( math_match_atom(dp->O,a2) == 0 ){a2 = a2->next; continue;}
		O = a2;
		dx = N->x - O->x;
		dy = N->y - O->y;
		dz = N->z - O->z;
		r = dx*dx + dy*dy + dz*dz;
	if( r > 5.2*5.2) { a2= a2->next; continue;}
		ron = sqrt(r);
		C = NULL;
		get_bond( O,bonded,10,&in_bond);
		for( k = 0; k< in_bond; k++)
			if( math_match_atom( dp->C,bonded[k]) !=0)
				{  C = bonded[k]; break;}
		if( C == NULL)  {a2 = a2->next; continue;}
		dx = N->x - C->x;
		dy = N->y - C->y;
		dz = N->z - C->z;
		rcn = sqrt(dx*dx+dy*dy+dz*dz);
		dx = H->x - C->x;
		dy = H->y - C->y;
		dz = H->z - C->z;
		rch = sqrt(dx*dx+dy*dy+dz*dz);
		dx = H->x - O->x;
		dy = H->y - O->y;
		dz = H->z - O->z;
		roh = sqrt(dx*dx+dy*dy+dz*dz);
	r = 332.17752*dp->k*(1./ron + 1./rch -1./roh - 1./rcn);
		if( r < -0.5){ // dssp cutoff criterion
			O->vx = 1.6*1.6; O->vy = 10.;
			H->vx = 1.6*1.6; H->vy = 10.;
			return 0;
		} 

	
		a2 = a2->next;
		} // for j
		}// if an N

		dp = dp->next;
	}
	return 0;
}
