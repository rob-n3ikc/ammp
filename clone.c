/* clone.c
*  keep copies of the current coordinates
*  return one or an average to the working set
*  sort of does what tether does, but somewhat
*  more elegantly
*  clone
*  statclone
*  subclone
*
*  nnclone
*
*  repclone  - find the most representative (lowest total rms with set) clone
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ammp.h"
/*#define graphic
*/

#ifdef graphic
void force_screen_update();
#endif

typedef struct {
int inme, serial;
float (*x)[],(*y)[],(*z)[];
char (*active_flag)[];
void *next;
}  CLONE;

/* nobody else should use firstCLONE, lastCLONE */
static CLONE *firstCLONE = NULL;
static CLONE *lastCLONE = NULL;
#define CLONG sizeof(CLONE)

static int N_clones = 0;

int clone(  who)
int who;
{
ATOM *ap,*a_next();
int i,numatom,a_number();
CLONE *cp,*cp2;
float (*x)[],(*y)[],(*z)[];
char (*active_flag)[];

	numatom = a_number();
	if( numatom < 1) return 0;
	cp = NULL;
	if( firstCLONE != NULL)
	{
		cp2 = firstCLONE;
		while( 1==1)
		{
			if( cp2->serial == who) {cp = cp2; break;}
			if( cp2->next == cp2) break;
			if( cp2->next == NULL) break;
			cp2 = cp2->next;
		}
		if( cp != NULL && cp->inme != numatom)
		{cp->inme = 0;free(cp->z); free(cp->y); free(cp->x);}
	}
	if( cp == NULL) {cp = malloc( CLONG); 
		if( cp == NULL) 
		{aaerror("cannot allocate memory in clone"); return 0;}
			N_clones += 1;
			cp->inme = 0; cp->next = cp;
			cp->serial = who;
			if( firstCLONE == NULL) firstCLONE = cp;
			if( lastCLONE == NULL) lastCLONE = cp;
			lastCLONE->next = cp;
			lastCLONE = cp; 
				}


	if( cp->inme == 0){
			cp->inme = numatom;
			cp->x = malloc( numatom *sizeof(float));
			cp->y = malloc( numatom *sizeof(float));
			cp->z = malloc( numatom *sizeof(float));
			cp->active_flag = malloc( numatom *sizeof(char));
			}

	x = cp->x; y = cp->y; z = cp->z;
	active_flag = cp->active_flag;
 	for( i=0; i< numatom ; i++)
	{
		ap = a_next(i);
		(*x)[i] = ap->x;
		(*y)[i] = ap->y;
		(*z)[i] = ap->z;
		(*active_flag)[i] = 0==1;
		if( ap->active) (*active_flag)[i] = (1==1);
	}

	return 0;
}/* end of clone */

/*  destructor 
*/
int declone( who)
int who;
{
if( firstCLONE == NULL) return 0;
}/* end of declone */

int statclone(op)
FILE *op;
{
int i,j,numclone, numatom , a_number();
ATOM *ap,*a_next();
float (*x)[],(*y)[],(*z)[];
float dx,dy,dz,dr;
int indx;
CLONE *cp;

numatom = a_number();
if( numatom < 1) return 0;
if( firstCLONE == NULL) return 0;


cp = firstCLONE;
i = 0;
while( 1==1)
{	
	i+= 1;
	if( cp ->next == cp) break;
	if( cp ->next == NULL) break;	
	cp = cp->next;
}
	fprintf(op," There are %d cloned coordinate sets\n",i);
numclone = i;

	cp = firstCLONE ;
	for( j=0; j< numclone; j++)
	{
	
fprintf(op," Clone %d label %d number of atoms %d",j+1,cp->serial,cp->inme);
		indx = numatom;
		if( indx < cp->inme) indx = cp->inme;
		x = cp->x; y = cp->y; z = cp->z;
		dr = 0.;
		for( i=0; i< indx; i++)
	 	{
		ap = a_next(i);
		dx = ap->x-(*x)[i];
		dy = ap->y-(*y)[i];
		dz = ap->z-(*z)[i];
		dr += dx*dx + dy*dy + dz*dz;
		}/* i */
		fprintf(op," rms %f\n", sqrt(dr/indx) );
	cp = cp->next;
	}/* j */


}/* end of statclone */

int restore_from_clone( who)
	int who;
{
	CLONE *cp; ATOM *ap, *a_next(int);
	int i;
	int a_number();
	cp = firstCLONE;
	if( cp == NULL) return 0;
	for( i=0; i< N_clones; i++)
	{
		if( cp->serial == who)
		{
			float (*x)[],(*y)[],(*z)[];
			int j,a;
			a = a_number();
			x = cp->x;
			y = cp->y;
			z = cp->z;
			for( j=0; j< a; j++)
			{
				ap = a_next(j);
				ap->x = (*x)[j];
				ap->y = (*y)[j];
				ap->z = (*z)[j];
			}
			return 1;
		}else{
			cp = cp->next;
		}
	}
	return 0;
}
int repclone( op,vfs,nfs, low, high)
int nfs,(*vfs[])();
int low,high;
FILE *op;
{
int i,j,k,numatom,a_number();
ATOM *ap,*a_next();
float (*x)[],(*y)[],(*z)[];
float (*xa)[],(*ya)[],(*za)[];
float bstrot(); /* defined in tether.c */
CLONE *cp,*cp1, *(*cps)[];
float matrix[3][3],delta[3];
float (*rmstable)[];
int nclone;
int ib;
float rb;
/* some basic error checks */
numatom = a_number();
if( numatom < 1) return 0;
if( firstCLONE == NULL) return 0;
// malloc memory that I know we'll need
cps = malloc( N_clones *sizeof( CLONE *));
rmstable = malloc( N_clones *sizeof( float ));
// if a range is given , fill the cps table with things in the range
// otherwise just fill it will all the CLONE* that are defined.
if(low == 0 && high == 0)
{
nclone = N_clones;
cp = firstCLONE;
	for( i =0; i< nclone; i++)
	{
		(*cps)[i] = cp;
		cp = cp->next;
	}
}else{
	nclone = 0;
	cp = firstCLONE;
	for( i =0; i< N_clones; i++)
	{
	if( cp->serial >=low && cp->serial <= high)
		(*cps)[nclone++] = cp;
	cp = cp->next;
	}
}
if( nclone == 0) return 0;
for( i=0; i< nclone; i++)
	(*rmstable)[i] = 0.;
/*
for( i=0; i< nclone; i++)
{
	cp = (*cps)[i];
	x = cp->x; y = cp->y; z = cp->z;
	if( cp->inme != numatom) { (*rmstable)[i] = 10.e10; continue;}
	for( j=0; j<nclone; j++)
	{
		if( j == i) continue;
		cp1 = (*cps)[j];
		if( cp1->inme != numatom) continue;
		xa = cp1->x; ya = cp1->y; za = cp1->z;
	(*rmstable)[i] += bstrot(*xa,*ya,*za, *x,*y,*z, numatom, matrix,delta);
	}// j
}//i 
*/
for( i=0; i<nclone; i++)
{
	cp = (*cps)[i];
	x = cp->x; y = cp->y; z = cp->z;
	if( cp->inme != numatom) { (*rmstable)[i] = 10.e10; continue;}
	for( j=0; j< numatom; j++)
	{
		ap = a_next(j);
		ap->x  = (*x)[j];
		ap->y  = (*y)[j];
		ap->z  = (*z)[j];
	}
	for( j=0; j< nfs; j++)
		 (*vfs[j])(&(*rmstable)[i],0.);
}

rb = (*rmstable)[0];
ib = 0;
for( i=1; i< nclone; i++)
   if( (*rmstable)[i] < rb){ ib = i;rb = (*rmstable)[i];}

cp = (*cps)[ib];
x = cp->x; y = cp->y; z = cp->z;
for( j=0; j< numatom; j++)
{
	ap = a_next(j);
	ap->x  = (*x)[j];
	ap->y  = (*y)[j];
	ap->z  = (*z)[j];
}



free(rmstable);
free(cps);
return 0;
}


int subclone(op,  who, how_many)
int how_many,who[];
FILE *op;
{
int i,j,k, numatom,a_number();
ATOM *ap,*a_next();
float (*x)[],(*y)[],(*z)[];
float (*xa)[],(*ya)[],(*za)[];
char (*act)[];
float bstrot(); /* defined in tether.c */
float matrix[3][3],delta[3];
CLONE *cp;

/*
printf("%d %d\n",how_many,who[0]);
*/
numatom = a_number();
if( numatom < 1) return 0;
if( firstCLONE == NULL) return 0;

/* first check that the clone is valid */
	if( how_many > 0)
	{
		j = -1;
		cp = firstCLONE;
		while(1==1)
		{
			for(i=0; i< how_many; i++)
			if( cp->serial == who[i]){j=1; break;}
			if( j == 1) break;
			if(cp->next == cp) break;
			if( cp->next == NULL) break;
			cp = cp->next;
		}
		if( j < 0){
		aaerror("no valid Clones\n"); return 0;
		}
	} 
	xa  = malloc( numatom*sizeof(float));
	ya  = malloc( numatom*sizeof(float));
	za  = malloc( numatom*sizeof(float));
	if( xa == NULL || ya == NULL || za == NULL){
		aaerror( "memory allocation error in subclone");
		exit(0);}
/*  zero the atoms */
	for( i=0; i< numatom ; i++)
	{
		ap = a_next(i);
		(*xa)[i] = ap->x;
		(*ya)[i] = ap->y;
		(*za)[i] = ap->z;
		ap->x = 0.;
		ap->y = 0.;
		ap->z = 0.;
		ap->dx = 0.;
	}
/* if how_many == 0 average them all */
	if( how_many == 0)
	{
		cp = firstCLONE;
		while( 1==1)
		{
		x = cp->x; y = cp->y; z = cp->z;
		j = numatom ; if( j < cp->inme) j = cp->inme;
	fprintf(op,"Clone %d %d atoms RMS %f\n", cp->serial,j,
		bstrot(*xa,*ya,*za, 
			*x,*y,*z, j, matrix,delta));
		for(i=0; i< j; i++)
		{
			ap = a_next(i);
ap->x += matrix[0][0]*(*x)[i]+matrix[1][0]*(*y)[i]+matrix[2][0]*(*z)[i];
ap->y += matrix[0][1]*(*x)[i]+matrix[1][1]*(*y)[i]+matrix[2][1]*(*z)[i];
ap->z += matrix[0][2]*(*x)[i]+matrix[1][2]*(*y)[i]+matrix[2][2]*(*z)[i];
		ap->x += delta[0];
		ap->y += delta[1];
		ap->z += delta[2];
			ap->dx += 1.;
		}/* for i */ 
		if(cp->next == cp) break;
		if( cp->next == NULL) break;
		cp = cp->next;
		}/* end of the while */

	} else{


		cp = firstCLONE;
		while( 1==1)
		{
		x = cp->x; y = cp->y; z = cp->z;
		for( k= 0; k< how_many; k++)
		{
		if( cp->serial == who[k]){
			j = numatom ; if( j < cp->inme) j = cp->inme;
			for(i=0; i< j; i++)
			{
				ap = a_next(i);
				ap->x += (*x)[i];
				ap->y += (*y)[i];
				ap->z += (*z)[i];
				ap->dx += 1.;
			}/* for i */ 
		}
		}/* for k */	
		if(cp->next == cp) break;
		if( cp->next == NULL) break;
		cp = cp->next;
		}/* end of the while */
	}/* end of how_many == 0 */
/* now average out the data */
	for( i=0; i< numatom; i++)
	{
		ap = a_next(i);
		if( ap->dx > 0.) ap->dx = 1./(ap->dx);
		ap->x *= ap->dx;
		ap->y *= ap->dx;
		ap->z *= ap->dx;
	}
/* update the active status */
	act = cp->active_flag;
	for( i=0; i< numatom; i++)
	{
		ap = a_next(i);
		ap->active = 1==0;
		if( (*act)[i] ) ap->active = 1==1;
	}
#ifdef graphic
	force_screen_update();
#endif
	free(za); free(ya); free(xa); return 0;

}/* end of subclone */

int nnclone(op,  who, how_many, rmin, rmax)
int how_many,who[];
float rmin,rmax;
FILE *op;
{
int i,j,k, numatom,a_number();
int iatom,jatom;
int incpa;
ATOM *ap,*ap1,*a_next();
CLONE *cp;
CLONE *(*cpa)[];
float (*xp1)[];
float (*yp1)[];
float (*zp1)[];
float x,y,z, r;
float rup,rdown, rbar;
float kup,kdown;
float (*rtable)[];

int noel();
float get_f_variable();

printf("%d %d\n",how_many,who[0]);
numatom = a_number();
if( numatom < 1) return 0;
if( firstCLONE == NULL) return 0;

if( rmin == 0.) rmin = 4.; /* short distances are covalent */
if( rmax == 0.) rmax = rmin+rmin;
if( rmin > rmax) 
{   x = rmin; rmin = rmax; rmax = x;}
kup = get_f_variable("kup");
kdown = get_f_variable("kdown");
if( kup < 0.001) kup = 10.;
if( kdown < 0.001) kdown = 10.;

cpa = malloc( how_many* sizeof( CLONE *) );
incpa = 0;
/* first check that the clone is valid */
for( i=0; i< how_many; i++)
	(*cpa)[i] = NULL;

	if( how_many > 0)
	{
		j = -1;
		cp = firstCLONE;
		while(1==1)
		{
			for(i=0; i< how_many; i++)
			if( cp->serial == who[i])
				{j=1; break;}
			if( j == 1) break;
			if(cp->next == cp) break;
			if( cp->next == NULL) break;
			cp = cp->next;
		}
		if( j < 0){ free(cpa);
		aaerror("no valid Clones\n"); return 0;
		}
	} 
	cp = firstCLONE;
	while(1==1)
	{
			for( i=0; i< how_many; i++)
					if( cp->serial == who[i])
						(*cpa)[incpa++] = cp;
			if( cp->next == cp) break;
			if( incpa == how_many) break;
			if( cp->next == NULL) break;
			cp = cp->next;
	}
	if( incpa < 3)
	{ free(cpa);  aaerror(" too few clones for nn_clone \n");
	   return 0;
	   }
/* so now cpa stores the pointers to the clones 
* 
*  loop over all atom pairs and if the pairs of atoms
*  and build a table
*  of distances 
*/
/* the following approach is not memory access efficient
* and will thrash like hell
* fortunately you only need to do i once in a while
*/
rtable = malloc( incpa * sizeof( float ));
	ap = a_next(-1);
for( iatom = 0; iatom < numatom-1; iatom++)
{

	ap1 = ap->next;
	for( jatom = iatom +1; jatom < numatom; jatom++)
	{
	
	for( i=0; i< incpa; i++)
	{
		xp1 = (*cpa)[i]->x;
		yp1 = (*cpa)[i]->y;
		zp1 = (*cpa)[i]->z;

		x = (*xp1)[iatom]-(*xp1)[jatom];
		y = (*yp1)[iatom]-(*yp1)[jatom];
		z = (*zp1)[iatom]-(*zp1)[jatom];

		r = sqrt( x*x +y*y+z*z);
		(*rtable)[i] = r;
	}/* i */
/* now we need to analyze the table of radii
* and by majority vote extract rmin,r,rmax
*/
	k = incpa/4 ;
	for( i=0; i< incpa-1; i++)
		for( j=i; j< incpa; j++)
			if( (*rtable)[i] > (*rtable)[j])
			{   x = (*rtable)[i]; (*rtable)[i] = (*rtable)[j];
			(*rtable)[j] = x; }
/*
			for( i=0; i< incpa; i++)
				printf("%f ",(*rtable)[i]);
			printf("\n");
*/
			rup = (*rtable)[incpa -k];
			rdown = (*rtable)[k];
			rbar = 0.5*(rup+rdown);
			if( rbar >= rmin && rbar <= rmax)
			noel( ap->serial, ap1->serial,
			rbar, rbar-rdown,rup -rbar, kdown,kup);
	ap1 = ap1->next;
	}/* jatom */
	ap = ap->next;
}/* iatom */


#ifdef graphic
	force_screen_update();
#endif
	free( rtable);
	free(cpa);
return 0;
}/* end of nnclone */

/* clone_average_distance(iat,jat min, radius, error)
*  iat, jat specify the atom pair
*  min is the minimum number of clones to use
*  radius is the maxium distance
*  error is the maximum deviation as a fraction o
*
*  return -1 if not allowed, or too few clones
*
*  return the average distance if ok
*/
float clone_average_distance(iat,jat, min, radius, error)
int min,iat,jat;
float radius,error;
{
CLONE *cp;
float (*xp)[];
float (*yp)[];
float (*zp)[];
float x,y,z, r;
float  rbar,r2bar;
int i;
	if( N_clones < min) return -1.;
	
	
	cp = firstCLONE;
	rbar = 0.;
	r2bar = 0.;
	
	for( i=0;i< N_clones; i++)
	{
		xp = cp->x;
		yp = cp->y;
		zp = cp->z;
		
		x = (*xp)[iat]-(*xp)[jat];
		y = (*yp)[iat]-(*yp)[jat];
		z = (*zp)[iat]-(*zp)[jat];
		r = x*x + y*y + z*z;
		

		rbar += r; 
/*
		r2bar += r*r;
*/

		cp = cp->next;
	}

	rbar = rbar/N_clones;
	if( rbar > radius) return -1.;
/*	
	r2bar = r2bar/N_clones - rbar*rbar;
	printf("%f %f\n", rbar,r2bar);
	if( r2bar > error*error*rbar*rbar) return -1.;
*/

	return rbar;

		
}/* end of clone_average_distance */
