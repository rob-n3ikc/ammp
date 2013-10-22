/* more experiments
 *  can we use the 4th eigenvector to correct for the distances that 
 *  we don't observe?
 *   
 *   in particular can we avoid the Floyd step since we need to anyway.
 *   if we work in the Fourier Domain then we can go directly from distances to
 *   the eigenvector squared. 
 */
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
void dgeom(echo,op, vfs,nfs,niter,origin,eigshift )
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
int v_restrain(), f_restrain();
int v_mmbond(),v_mmangle(),v_c_angle(),v_nonbon();
int v_hard(),v_abc();
int v_torsion(), v_step(), v_ttarget();
int f_bond(),f_angle(),f_noel();
int f_mmbond(),f_mmangle(),f_c_angle(),f_nonbon();
int f_hard(),f_abc();
int f_torsion(), f_step(), f_ttarget();
int v_hybrid(),f_hybrid(),use_hybrid;
int v_dssp_hbond(), f_dssp_hbond(), use_dssp;
float right,left;
int v_dill(),f_dill();

int jacobi(); /* in normal.c */
int a_v_zero();

int use_bond, use_angle, use_noel;
int use_torsion, use_ttarget, use_step;
int use_restrain;

int (*vuse[20])(), (*fuse[20])(),nuse;

float (*matrix)[];
float (*eig)[];
float (*raw)[];
float (*floyd)[]; /* keep the current floyded matrix */
float (*rawfloyd)[]; /* keep the current floyded matrix */
float (*scratch)[]; /* keep the current floyded matrix */
float (*bat)[];/* deflation matrix (die fleidermaus) */
ATOM *(*allatom)[];

numatm = a_number();
if( numatm < 2 ) return; /* the solution for 1 atom is trivial (0,0,0) */

matrix = (float (*)[])malloc( numatm*numatm*sizeof(float));
raw = (float (*)[])malloc( numatm*numatm*sizeof(float));
eig = (float (*)[])malloc( numatm*numatm*sizeof(float));
floyd = (float (*)[])malloc( numatm*numatm*sizeof(float));
rawfloyd = (float (*)[])malloc( numatm*numatm*sizeof(float));
scratch = (float (*)[])malloc( numatm*numatm*sizeof(float));
/* bat will now hold the basis */
bat = (float (*)[])malloc( numatm*numatm*sizeof(float));

allatom = (ATOM *(*)[])malloc( numatm *sizeof(ATOM *));


for( i=0; i< numatm; i++)
{ (*allatom)[i] = a_next(i);}


use_bond = (1==0);
use_angle = (1==0);
use_noel = (1==0);
use_step = (1==0);
use_torsion = (1==0);
use_hybrid = (1==0);
use_restrain = (1==0);
use_ttarget = (1==0);
use_dssp = (1==0);
nuse = 0;
for( i=0; i< nfs; i++)
{
	if( vfs[i] == v_bond) {use_bond = (1==1); vuse[nuse] = v_bond; fuse[nuse] = f_bond; nuse++;}
	if( vfs[i] == v_restrain) {use_restrain = (1==1); vuse[nuse] = v_restrain; fuse[nuse] = f_restrain; nuse++;}
	if( vfs[i] == v_angle) {use_angle = (1==1); vuse[nuse] = v_angle; fuse[nuse] = f_angle; nuse++;}
	if( vfs[i] == v_noel) {use_noel = (1==1); vuse[nuse] = v_noel; fuse[nuse] = f_noel; nuse++;}
	if( vfs[i] == v_step) use_step = (1==1);
	if( vfs[i] == v_ttarget) use_ttarget = (1==1);
	if( vfs[i] == v_torsion) {use_torsion = (1==1); vuse[nuse] = v_torsion; fuse[nuse] = f_torsion; nuse++;}
	if( vfs[i] == v_hybrid) {use_hybrid = (1==1); vuse[nuse] = v_hybrid; fuse[nuse] = f_hybrid; nuse++;}
	if( vfs[i] == v_dill) { vuse[nuse] = v_dill; fuse[nuse] = f_dill; nuse++;}
	if( vfs[i] == v_dssp_hbond) { vuse[nuse] = v_dssp_hbond; fuse[nuse] = f_dssp_hbond; nuse++;}
}
vuse[nuse] = v_nonbon; fuse[nuse++] = f_nonbon;

	for( i=0; i< numatm*numatm; i++)
		(*bat)[i] = 0.;

	iter = 0;
	for( i=0; i< numatm; i++)
	{
	a_v_zero();
	if( use_bond) gsdg_bond((*allatom)[i]);
	if( use_angle) gsdg_angle((*allatom)[i]);
	if( use_restrain) gsdg_restrain((*allatom)[i]);
/*
//	if( use_torsion) gsdg_torsion((*allatom)[i]);
//	if( use_ttarget) gsdg_ttarget((*allatom)[i]);
*/
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
//#define KEEP_CURRENT
*/
/*
 * before floyd algorithm
 * we should fill the distances in from the current
 * structure.   
 * floyd will catch inconsistencies between experimental and
 * theoretical distances.
 */
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
		(*raw)[i*numatm +j] = r;	
		(*raw)[j*numatm +i] = r;	
			
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
	if( l< numatm)
	for( i=0; i<numatm*numatm; i++)
	{
		r = (*eig)[i];
		(*floyd)[i] = r;
		if( r > 10.e6)
		{
		goto REFLOYD;
		} 
	}
	
/*
//	dgeom_triangle( floyd,raw,numatm);
*/
	for( i=0; i< numatm*numatm; i++)
		(*rawfloyd)[i] = (*eig)[i]; /* a backup copy for later on */
#else
	for( i=0; i< numatm; i++)
	for( j=0; j< numatm; j++)
		if( (*eig)[i*numatm +j] > 10.e6 ) {
			(*eig)[i*numatm + j] = (float)(i+j);
			(*rawfloyd)[i*numatm + j] = 0.;
		}
/*		if( (*eig)[i*numatm +j] > 10.e6 ) (*eig)[i*numatm + j] = 0.;
 *		*/
	 for( i=0; i< numatm; i++)
		 for( j=0; j< numatm; j++)
		  (*floyd)[i*numatm+j] = (*eig)[i*numatm+j];

/*	for( i=0; i< numatm*numatm; i++)
		(*rawfloyd)[i] = (*eig)[i]; /* a backup copy for later on */
#endif
	}
/* generate the basis set in (*bat)  it will be sin/cos (n pi /numatm) n = 1,2,...*/
	for( i=0; i< numatm; i++)
		(*bat)[i] = 0.;
	for( i=1; i< numatm; i++)
	{
		float konst;
		if( i%2 == 0)
		konst = 2*3.14159265*(i)/2;
		else
		konst = 2*3.14159265*(i+1)/2;
	for( j = 0; j< numatm; j++)
	{
		float kkk;
			kkk = (konst *j )/numatm;
			if( i%2 == 0)
		(*bat)[i*numatm+j] = sin(kkk);
			else
		(*bat)[i*numatm+j] = cos(kkk);
	}
	}
/* check orthogonal (as a test )*/
	for( i=0; i< numatm; i++)
	{
		float xxx;
		for( j=0; j< numatm; j++)
		{
		xxx = 0.;
		for( k=0; k< numatm; k++)
			xxx += (*bat)[i*numatm+k]*(*bat)[j*numatm+k];
		}
		j = i;
		xxx = 0.;
		for( k=0; k< numatm; k++)
			xxx += (*bat)[i*numatm+k]*(*bat)[j*numatm+k];
		if( xxx <= 0.00001) continue;
		xxx = sqrt(1./xxx);
		for( k=0; k< numatm; k++)
			(*bat)[i*numatm+k] *= xxx;



	}
/* so now we have distances in floyd, eig
 * and a basis set in bat
 * raw possesses a copy of the input distances and the constrained ones from 
 * the inactive atoms.
 */
	for( iter=0; iter< niter; iter++)
	{
#ifdef TEST_DISTANCES
		/* test if the distance is ok */
		/* doesn't seem to make much difference */
        for( i=0; i< numatm; i++)
	        for( j=0; j< numatm; j++)
	        {
			if( (*raw)[i*numatm+j] < 10.e10) continue;
			// was 25
			if( (*eig)[i*numatm+j] < 1.8*1.8) (*eig)[i*numatm+j] = 100.;
		}
#endif
#define DISTANCE_CLASSIC
#ifdef DISTANCE_CLASSIC
        r = 0.;
        x = 0.;
        for( i=0; i< numatm; i++)
	        for( j=0; j< numatm; j++)
	        {
                   x1 = (*eig)[i*numatm +j];
                    {r += x1; x += 1.;}
            }
        r /= (x);

	        for( i=0; i< numatm; i++)
	        for( j=0; j< numatm; j++)
	        {
	                xo = 0.;
	                yo = 0.;
	                x = 0.;
	                y = 0.;
	                for( k=0; k< numatm; k++)
	                {
                          x1 = (*eig)[i*numatm +k];
                          y1 = (*eig)[k*numatm +j];
                          { xo += x1; x += 1.;}
                          { yo += y1; y += 1.;}
                         }
                        if( x > 0.) xo /= x;
                        if( y > 0.) yo /= y;
              (*matrix)[i*numatm + j] = -0.5*((*eig)[i*numatm+j] -xo - yo+r);
               }
	for( i=0; i< numatm*numatm; i++)
		(*floyd)[i] = (*matrix)[i];
#endif
/* transform floyd into matrix */

	dgeom_apply_basis( matrix, scratch, floyd,bat, numatm);	

	simultaneous_iteration(   matrix,eig, scratch, NULL,(1==1), numatm,4,100,0.001);
	printf("%f %f %f %f\n",(*scratch)[0],(*scratch)[1], (*scratch)[2], (*scratch)[3]);
/*	dgeom_deflate( matrix,eig,scratch,numatm,3);
//`	dgeom_unapply_basis( floyd,scratch,matrix,bat, numatm);
//`	dgeom_restore_constraint( floyd,raw,numatm);
*/
	dgeom_update( allatom, eig,scratch,bat,numatm);

	if( nuse == 0) continue;
#define hybrid_test
#ifdef hybrid_test
	if( use_hybrid)
	{/* check the hand of the solution.  right hand solution should have minimal hybrids */
		v_hybrid(&right,0.);
		for( i=0;i< numatm; i++)
		{
			(*allatom)[i]->z *= -1.;
		}
		v_hybrid(&left,0.);
		if( left > right)
		for( i=0;i< numatm; i++)
		{
			(*allatom)[i]->z *= -1.;
		}
	}
#endif
	// was 2,0
	a_build(nuse, vuse,fuse, 10,0,(*allatom)[numatm-1]->serial+1,stdout,(1==0));
	cngdel( vuse,fuse, nuse, 100, 100,0.,(1==0));
	for( i=0; i< numatm; i++)
		for( j=0; j<numatm; j++)
		{
			ATOM *aap,*bbp;
			float xxx,yyy,zzz,rrr;
			if( i == j) (*eig)[i*numatm+j] = 0.;
			else if( (*raw)[i*numatm+j] < 10.e6) 
			{/* distance is known */
				(*eig)[i*numatm+j] = (*raw)[i*numatm+j];
			}else {
				aap = (*allatom)[i];
				bbp = (*allatom)[j];
				xxx = aap->x -bbp->x;
				yyy = aap->y -bbp->y;
				zzz = aap->z -bbp->z;
				rrr = xxx*xxx + yyy*yyy + zzz*zzz;
				if( rrr < (*rawfloyd)[i*numatm+j])
				       	rrr = (*rawfloyd)[i*numatm+j];

/*				(*rawfloyd)[i*numatm+j] = rrr;
 */
				(*eig)[i*numatm+j] = rrr;
			
			}
		}

	} /* iter */
/* clean up at the end of the routine */
free(allatom);
free(scratch);
free(bat);
free(rawfloyd);
free(floyd);
free(raw);
free(eig);
free(matrix);
}/* end of dgeom */
int	dgeom_update( allatom,eig,scratch,bat,numatm)
	int numatm;
	ATOM* (*allatom)[];
	float  (*eig)[], (*scratch)[], (*bat)[];
	/* bat is the basis */
	/* the eigen values are the first n elements of scratch   */
	/* the eigenvectors are in eig  */
	/* the eigen vectors are the rows of eig */
{
	int i,j,k;
	float e;
	ATOM *ap;
	
	for( i=0; i< numatm; i++)
	{ ap = (*allatom)[i]; ap->x = 0.; ap->y = 0.; ap->z = 0.;}
/* x coordinate first */
	e = -0.5*sqrt( (*scratch)[0])/numatm;
	e = sqrt((*scratch)[0]);
	for( i=0; i< numatm; i++)
	{
		ap = (*allatom)[i];
		for( j=0; j< numatm; j++)
			ap->x += e*(*eig)[j]*(*bat)[i*numatm + j];
	}/* i */

/* y coordinate  */
	e = -0.5*sqrt( (*scratch)[1])/numatm;
	e = sqrt((*scratch)[1]);
	for( i=0; i< numatm; i++)
	{
		ap = (*allatom)[i];
		for( j=0; j< numatm; j++)
			ap->y += e*(*eig)[numatm+j]*(*bat)[i*numatm + j];
	}/* i */

/* z coordinate  */
	e = -0.5*sqrt( (*scratch)[2])/numatm;
	e = sqrt((*scratch)[2]);
	for( i=0; i< numatm; i++)
	{
		ap = (*allatom)[i];
		for( j=0; j< numatm; j++)
			ap->z += e*(*eig)[numatm+numatm+j]*(*bat)[i*numatm + j];
	}/* i */

	

	return 0;
}/* end of dgeom_update() */

int dgeom_apply_basis( out, scratch, in, basis, n)
	int n;
	float (*out)[],(*scratch)[], (*in)[], (*basis)[];
{
	int i,j,k;
	for( i=0; i< n*n; i++)
	{ (*out)[i]  = 0.; (*scratch)[i] = 0.;}

	for( i=0; i< n; i++)
		for( j = 0; j < n; j++)
		for( k = 0; k < n; k++)
		{ (*scratch)[i*n+j] += (*in)[i*n+k]*(*basis)[k*n+j];}

	/* remember to use transpose!!! */
	for( i=0; i< n; i++)
		for( j = 0; j < n; j++)
		for( k = 0; k < n; k++)
		{ (*out)[i*n+j] += (*basis)[k*n+i]*(*scratch)[k*n+j];}
	return 1;
}/* end of dgeom_apply_basis */
int	dgeom_deflate( matrix,eig,scratch,n,which)
	float (*matrix)[], (*eig)[],  (*scratch)[] ;
	int n,which;
{
	/* the eigen values are the first n elements of scratch   */
	/* the eigen vectors are the rows of eig */
	float e ;
	int i,j,offset;
	e = (*scratch)[which];
	/* e may need to be normalized  */
	e = e/n;
	printf("%f %f\n", e, (*scratch)[which]);
	offset = n*which;
	for( i=0; i< n; i++)
	       for( j = 0; j< n; j++)
	       (*matrix)[i*n+j] -= e*(*eig)[i*n+which]*(*eig)[offset+j];
	return 0;	
}/* end of dgeom deflate */
int	dgeom_unapply_basis( floyd,scratch,matrix,bat, n)
	float (*floyd)[], (*scratch)[], (*matrix)[], (*bat)[];
	int n;
{
	int i,j,k;
	for( i=0; i< n*n; i++)
		{(*scratch)[i] = 0.; (*floyd)[i] = 0.;}
	for( i=0; i< n; i++)
		for( j=0; j< n; j++)
			for( k = 0; k< n; k++)
				(*scratch)[i*n+j] += (*matrix)[i*n+k]*(*bat)[k*n+j];
	/* remember to use transpose!!! */
	for( i=0; i< n; i++)
		for( j=0; j< n; j++)
			for( k = 0; k< n; k++)
				(*floyd)[i*n+j] += (*bat)[k*n+i]*(*scratch)[k*n+j];

	for( i=0; i< n*n; i++)
		(*floyd)[i] /= n;
	return 0;
}/* end of dgeom_unapply_basis */
int	dgeom_restore_constraint( floyd,raw,numatm)
	float (*floyd)[],  (*raw)[];
	int numatm;
{
	int i;
	for( i=0; i< numatm*numatm; i++)
	{
		if( (*floyd)[i] < 2.) (*floyd)[i] = 2.; /* miniminum vdw */
		if( (*raw)[i] < 1.e10) (*floyd)[i] = (*raw)[i];  
	}
	return 0;
}/* end of dgeom_restore_constriant */

int simultaneous_iteration( matrix, eig, values,scratch,initialize, n,nfit, niter, toler)
	float (*matrix)[], (*eig)[], (*values)[], (*scratch)[], toler;
	int n,nfit,niter;
{
	int iter,i,j,k,free_scratch;
	float deig,x,y;
	float randf();

		free_scratch = (1==0);
		if( scratch  == NULL)
		{
			free_scratch = (1==1);
			scratch = ( float (*)[])malloc( n* sizeof(float));
		}		
		if( initialize)
		{
			for( i=0; i< nfit*n; i++)
				(*eig)[i] = randf();
		}
			for( i=0; i< nfit; i++)
				(*values)[i] = 0.; /* this ensures that it iterate at least once */
		
		/* normalize */
		for( i=0; i< nfit; i++)
		{
			x = 0.;
			for( j=0; j< n; j++)
				x += (*eig)[i*n+j] * (*eig)[i*n+j];
			if( x < 1.e-7) aaerror("suspect zero eigenvector in simultaneous_iteration input");
			x = 1./sqrt(x);
			for( j=0; j< n; j++)
				(*eig)[i*n+j] *= x;

		}
		for( iter=0; iter < niter; iter++)
		{

			 /* BE_NICE */ ;
			/*first find the products */
			for( k=0; k< nfit; k++)
			{
			for( i=0; i< n; i++)
				{
				(*scratch)[i] = 0.;
				for( j=0; j< n; j++)
				 (*scratch)[i] += (*matrix)[i*n+j]*(*eig)[k*n+j];
				} /*  i */

			for( i=0; i< n; i++)
				(*eig)[k*n+i] = (*scratch)[i];
			}/* k */
			/* then ortho-normalize them with Graham-schmidt orthogonalization */
			deig = 0.;
			for( i=0; i< nfit; i++)
			{
				x = 0.;
				for( j=0; j< n; j++)
					x += (*eig)[i*n+j] * (*eig)[i*n+j];
				x = sqrt(x);
				deig += fabs( x - (*values)[i]);	
				(*values)[i] = x;
				x = 1./x;
				for( j=0; j< n; j++)
					(*eig)[i*n+j] *= x;
			}/* i */
			/* now normal, but not orthogonal */
			for( k=0; k< nfit-1; k++)
			for( i=k+1; i < nfit; i++)
			{
				/* dot product */
				x = 0.;
				for( j=0; j< n; j++)
					x += (*eig)[k*n+j] * (*eig)[i*n+j];
				/* subtract */
				for( j=0; j< n; j++)
					(*eig)[i*n+j] -= x* (*eig)[k*n+j];
				/* re-normalize */
				x = 0.;
				for( j=0; j< n; j++)
					x += (*eig)[i*n+j] * (*eig)[i*n+j];
				x = 1./sqrt(x);
				for( j=0; j< n; j++)
					(*eig)[i*n+j] *= x;
			}/* k i */
			if( deig < toler ) break;

		}/* iter */	


		if( free_scratch) free(scratch);
}/* end of simultaneous_iteration */


/*
 * terhaar.c
 *
 * driver for a one-dimensional terhaar wavelet
 *
 * basic idea follows the conte-deboor fft (uugh)
 * where a buffer is used so that the bit-swapped output
 * happens automagically
 *
 * UNLIKE THE FFT n really does have to be a power of two
 *
 * rwh 2001
 *
 * 
 */


#ifdef SINGLE_LINE
int  terhaar( float x[], float scratch[], int n)
#else
int  terhaar( x,scratch,n)
float x[];
float scratch[];
int n;
#endif
{

	int i,j,delta,stride;
	float xsum,xdiff;
	
	for( i=2; i <= n ; i *= 2)
	{
		if( i == n) goto N_IS_A_POWER_OF_TWO;
	}	

	fprintf(stderr,"terhaar requires a power of 2 %d is not\n", n);
	return 1==0;

N_IS_A_POWER_OF_TWO: ;

	for( i=1; i<n; i*=2)
	{
		delta = 2;
		stride = n/2;
		for( j=0; j< n/2; j++)
		{
			xsum = x[j+j] + x[j+j+1];
			xdiff = x[j+j] - x[j+j+1];

			scratch[j] = xsum;
			scratch[stride+j] = xdiff;

		}/*inc */
		for( j=0; j< n; j++)
		{
			x[j] = scratch[j];
		}
	} /* i */
} /* end of terharr */

