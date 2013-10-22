/*
* optimist.c
*
* optimiser routines for AMMP
*
*  steepest descents
*  conjugate gradients
*  line search    routines
*
*/
/*
*  copyright 1992 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/

#define ANSI 1
#include <stdio.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif

#include "ammp.h"
/*
exclude outliers in a second pass of block_optimizer

*/

int block_optimizer( vfs,ffs,nfs,nstep, nblock, echo)
int nfs,(*vfs[])(),(*ffs[])();
int nstep,nblock,echo;
{
	float *xp,*yp,*zp, *xa,*ya,*za;
	ATOM *ap,*a_next(int), **alist;
	int natom,a_number(void);
	int i,j,k;  
	int aaerror(char*);
	float m[3][3],v[3];
	float bstrot();

	natom =  a_number();
	if( natom < nblock *3){ aaerror(" blockopt  nblock too big\n"); return 1==0;}
	if( natom < 1){ aaerror("need some atoms to minimize\n"); return 1==0;}

	alist = malloc( natom*sizeof(ATOM *));
	xp = malloc( natom*sizeof(float));
	yp = malloc( natom*sizeof(float));
	zp = malloc( natom*sizeof(float));
	xa = malloc( natom*sizeof(float));
	ya = malloc( natom*sizeof(float));
	za = malloc( natom*sizeof(float));

	for( i=0; i< natom; i++)
	{
	ap = a_next(i);
	alist[i] = ap;
	xp[i] = ap->x;
	yp[i] = ap->y;
	zp[i] = ap->z;
	}
	cngdel( vfs,ffs,nfs, nstep, nstep, 0., echo);
/*  not quite as good on (too simple) a test
	pc_cngdel( vfs,ffs,nfs, nstep, nstep, 0., echo);
*/

	for( i=0; i< natom; i++)
	{
	ap = alist[i];
	xa[i] = ap->x;
	ya[i] = ap->y;
	za[i] = ap->z;
	}

/* first one */
	i = 0;
	bstrot_out(&xp[i],&yp[i],&zp[i], &xa[i],&ya[i],&za[i],3*nblock,m,v);
	for( j=0; j< nblock; j++)
	{
	ap = alist[j];
	if( !ap->active) continue;
	ap->x = v[0]+m[0][0]*xa[j]+m[1][0]*ya[j]+m[2][0]*za[j];
	ap->y = v[1]+m[0][1]*xa[j]+m[1][1]*ya[j]+m[2][1]*za[j];
	ap->z = v[2]+m[0][2]*xa[j]+m[1][2]*ya[j]+m[2][2]*za[j];
	}
/* middle ones */
	for( i=0; i< natom; i+= nblock)
	{
	k = nblock*3;
	if( i+k > natom) k = natom-i;
	bstrot_out(&xp[i],&yp[i],&zp[i], &xa[i],&ya[i],&za[i],k,m,v);
	for( j=nblock+i; j< nblock+nblock+i; j++)
	{
		ap = alist[j];
		if( !ap->active) continue;
		ap->x = v[0]+m[0][0]*xa[j]+m[1][0]*ya[j]+m[2][0]*za[j];
		ap->y = v[1]+m[0][1]*xa[j]+m[1][1]*ya[j]+m[2][1]*za[j];
		ap->z = v[2]+m[0][2]*xa[j]+m[1][2]*ya[j]+m[2][2]*za[j];
	}/* j */
	if( natom == i+k) break;
	}/* i */
/* last ones */
	bstrot_out(&xp[i],&yp[i],&zp[i], &xa[i],&ya[i],&za[i],k,m,v);
	for( j=i+nblock+nblock; j< natom; j++)
	{
	ap = alist[j];
	if( !ap->active) continue;
	ap->x = v[0]+m[0][0]*xa[j]+m[1][0]*ya[j]+m[2][0]*za[j];
	ap->y = v[1]+m[0][1]*xa[j]+m[1][1]*ya[j]+m[2][1]*za[j];
	ap->z = v[2]+m[0][2]*xa[j]+m[1][2]*ya[j]+m[2][2]*za[j];
	}

	free(za);
	free(ya);
	free(xa);
	free(zp);
	free(yp);
	free(xp);
	free(alist);
	return 1==1;
}

int bstrot_out(x,y,z, xx,yy,zz, k, m,v)
float x[],y[],z[];
float xx[],yy[],zz[];
int k;
float m[3][3],v[3];
{
	float *x1,*y1,*z1,*xx1,*yy1,*zz1;
	int i,j;
	float me2;
	float dx;
	float bstrot();
	x1 = malloc( k*sizeof(float));
	y1 = malloc( k*sizeof(float));
	z1 = malloc( k*sizeof(float));
	xx1 = malloc( k*sizeof(float));
	yy1 = malloc( k*sizeof(float));
	zz1 = malloc( k*sizeof(float));
	bstrot(x,y,z,xx,yy,zz,k,m,v);
	
	me2 = 0.;
	for( j=0; j< k; j++)
	{
	x1[j] = v[0]+m[0][0]*xx[j]+m[1][0]*yy[j]+m[2][0]*zz[j];
	y1[j] = v[1]+m[0][1]*xx[j]+m[1][1]*yy[j]+m[2][1]*zz[j];
	z1[j] = v[2]+m[0][2]*xx[j]+m[1][2]*yy[j]+m[2][2]*zz[j];
	me2 += (x1[j]-x[j])*(x1[j]-x[j]);
	me2 += (y1[j]-y[j])*(y1[j]-y[j]);
	me2 += (z1[j]-z[j])*(z1[j]-z[j]);
	}
	me2 /=k;
/* 3 sigma (squared)
//	me2 *=3.;
//	me2 *=3.;
*/
	me2 *= 4.;
	i = 0;
	for( j=0; j< k; j++)
	{
	dx = 0.;
	dx += (x1[j]-x[j])*(x1[j]-x[j]);
	dx += (y1[j]-y[j])*(y1[j]-y[j]);
	dx += (z1[j]-z[j])*(z1[j]-z[j]);
	if( dx > me2) continue;
	x1[i] = x[j];
	y1[i] = y[j];
	z1[i] = z[j];
	xx1[i] = xx[j];
	yy1[i] = yy[j];
	zz1[i] = zz[j];
	i++;
	}
	
	if( i > 6) bstrot(x1,y1,z1,xx1,yy1,zz1,i,m,v);
	free(zz1);
	free(yy1);
	free(xx1);
	free(z1);
	free(y1);
	free(x1);
	return 0;
}
/* function steep()
* steepest descents optimiser
*
* passed a set of functions vfs for the potential
* passed a set of functions ffs for the force
*  how many functions  nfs
*  how many steps to try
*  when to stop
*/
int steep( vfs,ffs,nfs, nstep,toler )
int nfs,(*vfs[])(),(*ffs[])();
int nstep;
float toler;
{
    int i,j,ifs;
    float vb,vt,vto;
    float a_max_f(),lam;
    /*  do at most nstep steps */
    for( i=0; i< nstep ; i++)
    {
        vb = 0.;
        for( ifs = 0; ifs < nfs; ifs++)
            (*vfs[ifs])( &vb,  0.);
        a_f_zero();
        for( ifs = 0; ifs < nfs; ifs++)
            (*ffs[ifs])( 0.);
        lam = a_max_f();
        if( lam <= toler) return 1;
        if( lam <= 1. )
        { a_ftodx(1.,0.); } else { a_ftodx( 1./sqrt(lam),0. ); }
        vto = vb;
        lam = 0.;
        /* this is a fibonacci line search
        * rather odd, the steps keep getting
        * bigger until
        * either the end
        * or we go up
        */
        for( j=0; j< 200; j++)
        {
            vt = 0.;
            lam = lam + .01*j;
            for( ifs = 0; ifs < nfs; ifs++)
                (*vfs[ifs])(&vt, lam);
            if( vt > vto && j == 0 )
            { a_ftodx(0., .1);   }
            else if( vt > vto )
            { lam = lam -.01*j; break; }
            vto = vt;
        }
        a_inc_d( lam );
    }
    return 0;
}
/* function cngdel()
* conjugate gradients optimiser
*
* passed a set of functions vfs for the potential
* passed a set of functions ffs for the force
*  how many functions  nfs
*  how many steps to try
*  when to stop
*/
int cngdel( vfs,ffs,nfs, nstep,nreset,toler,echo )
int nfs,(*vfs[])(),(*ffs[])();
int nstep,nreset;
float toler;
int echo;
{
    int i,j,ifs;
    float vb,vt,vto;
    float a_max_f(),lam;
    float linmin(),a_max_d();
    float a_pr_beta();
    float beta,betad,a_l2_f(),a_l2_g();
    int insane();/* check status */
    /*  do at most nstep steps */
    if( nreset < 1) nreset = nstep;
    a_g_zero();
    a_d_zero();
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
    for( i=0; i< nstep ; i++)
    {
        if( insane() ) return 0;
        vb = 0.;
        for( ifs = 0; ifs < nfs; ifs++)
            (*vfs[ifs])( &vb,  0.);
        a_f_zero();
        for( ifs = 0; ifs < nfs; ifs++)
            (*ffs[ifs])( 0.);
        lam = a_max_f();
        /* make up the beta *
        /* use just the simple one */
        /*	beta = a_l2_f();
        	betad = a_l2_g();
        	if( betad < 1.e-4) {betad = 1.;beta = 0.;}
        	beta = -beta/betad;
        */
        beta = a_pr_beta();

        if( (i%nreset) == 0) beta = 0.;
        if( echo) printf(" vb, maxf %f %f %f \n",vb,lam,beta);
        if( lam <= toler) return 1;
        /*	printf(" max f max d  %f %f \n",a_max_f(),a_max_d() ); */
        /* and make up the descent direction */
        lam = a_max_f();
        a_ftodx(1.,beta);
        a_ftogx( 1.,0.);
        lam = linmin( vfs,nfs, sqrt(a_max_d()) );
        if( lam < 1.e-6)
        {
            a_f_zero();
            for( ifs = 0; ifs < nfs; ifs++)
                (*ffs[ifs])(0. );
            lam = sqrt(a_max_f()); if( lam < 1.) lam = 1.;
            a_ftodx( 1./(lam), 0.);
            lam = linmin(vfs,nfs,1. );
            if( lam < 1.e-6) return 0;
        }
        a_inc_d( lam );
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
    }
    return 0;
}



/*line minimization routine linmin( x, search, n, f )
*
* returns step size for approximate line minimizer
*  recursive line minimizer
*/
float linmin(ffs,nfs ,damp )
int (*ffs[])();
int nfs;
float damp;

{
    int i,iter,jter,imin;
    float alpha[401],fval[401],dstep,step,stpmin,fvt,fmin,fold;
    int nostep,get_i_variable();
    float mxdq, get_f_variable();

    nostep = get_i_variable("nostep");
    if( nostep < 1) nostep = 8;
    mxdq = get_f_variable("mxdq");
    set_f_variable("mxdq" , 100.);

    dstep = 1.;
    step = 0;
    stpmin = 0.;
    imin = 0;
    alpha[0] = 0;
    fval[0] = 0;
    /*	 if( damp < 1.) {dstep = .25;}else{ dstep = .01/sqrt(damp);} */
if( damp < 1.) {dstep = .25;}else{ dstep = 1./damp;}

    for( i=0; i< nfs ; i++)
        (*ffs[i])(&fval[0],0.);
    fmin = fval[0]; fold = fmin;
    /* recheck is to find a descent stp length */
    imin++;
recheck:
    alpha[imin] = dstep;
    fval[imin] = 0.;
    for( i=0; i< nfs ; i++)
        (*ffs[i])(&fval[imin],dstep);
    if( fval[imin] > fval[0])
    {
        dstep = dstep*.25;
        if( dstep > 1.e-8) goto recheck;
        goto cleanup;
    }
    /* if here we have dstep small enough to use */
    /*	for( iter=0; iter< 8; iter++)
    */
    for( iter=0; iter< nostep; iter++)
    {
        for(jter=1; jter<100; jter++)
        {
            step += dstep;
            /* get the function value */
            for( i=0; i< imin; i++)
            {
                if( alpha[i] == step)
                {
                    fvt = fval[i]; goto FOUND ;
                }
            }
            fvt = 0.;
            for( i=0; i< nfs ; i++)
                (*ffs[i])(&fvt,step);
            alpha[imin] = step;  fval[imin++] = fvt;
            if( imin > 400) imin = 400;
FOUND:
            /*		printf("\n %f %f %f %f", step,fvt,fmin,fold);  */
            if( fvt< fmin)
            {
                fmin = fvt; stpmin = step;
            }
            if( fvt > fold)
            {
                dstep = -dstep/2; break;
            }
            fold = fvt;
        }

    }
cleanup:
    /* do not allow 0 steps */
    /*	if( stpmin == 0.) stpmin = dstep/2; */
    set_f_variable("mxdq" , mxdq);
    return(stpmin);
}
/* function pc_cngdel()
* conjugate gradients optimiser
*
* pc_cngdel uses a preconditioner (basically a floating tether)
* to avoid certain nasties.
*
* passed a set of functions vfs for the potential
* passed a set of functions ffs for the force
*  how many functions  nfs
*  how many steps to try
*  when to stop
*/
int pc_cngdel( vfs,ffs,nfs, nstep,nreset,toler,echo )
int nfs,(*vfs[])(),(*ffs[])();
int nstep,nreset;
float toler;
int echo;
{
    int i,j,ifs;
    float vb,vt,vto;
    float a_max_f(),lam;
    float pc_linmin(),a_max_d();
/* pc_linmin() actually impliments the preconditioner */
    float a_pr_beta();
    float beta,betad,a_l2_f(),a_l2_g();
    int insane();/* check status */
    /*  do at most nstep steps */
    if( nreset < 1) nreset = nstep;
    a_g_zero();
    a_d_zero();
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
    for( i=0; i< nstep ; i++)
    {
        if( insane() ) return 0;
        vb = 0.;
        for( ifs = 0; ifs < nfs; ifs++)
            (*vfs[ifs])( &vb,  0.);
        a_f_zero();
        for( ifs = 0; ifs < nfs; ifs++)
            (*ffs[ifs])( 0.);
        lam = a_max_f();
        /* make up the beta *
        /* use just the simple one */
        /*	beta = a_l2_f();
        	betad = a_l2_g();
        	if( betad < 1.e-4) {betad = 1.;beta = 0.;}
        	beta = -beta/betad;
        */
        beta = a_pr_beta();

        if( (i%nreset) == 0) beta = 0.;
        if( echo) printf(" vb, maxf %f %f %f \n",vb,lam,beta);
        if( lam <= toler) return 1;
        /*	printf(" max f max d  %f %f \n",a_max_f(),a_max_d() ); */
        /* and make up the descent direction */
        lam = a_max_f();
        a_ftodx(1.,beta);
        a_ftogx( 1.,0.);
        lam = pc_linmin( vfs,nfs, sqrt(a_max_d()) );
        if( lam < 1.e-6)
        {
            a_f_zero();
            for( ifs = 0; ifs < nfs; ifs++)
                (*ffs[ifs])(0. );
            lam = sqrt(a_max_f()); if( lam < 1.) lam = 1.;
            a_ftodx( 1./(lam), 0.);
            lam = pc_linmin(vfs,nfs,1. );
            if( lam < 1.e-6) return 0;
        }
        a_inc_d( lam );
#ifdef GTK_VERSION
AMMP_whiz_bang();
#endif
    }
    return 0;
}



/*line minimization routine linmin( x, search, n, f )
*
* returns step size for approximate line minimizer
*  recursive line minimizer
*/
float pc_linmin(ffs,nfs ,damp )
int (*ffs[])();
int nfs;
float damp;

{
    int i,iter,jter,imin;
    float alpha[401],fval[401],dstep,step,stpmin,fvt,fmin,fold;
    int nostep,get_i_variable();
    float mxdq, get_f_variable();
	float a_l2_d(),sumdx;
	int natom,a_number();

    nostep = get_i_variable("nostep");
    if( nostep < 1) nostep = 8;
    mxdq = get_f_variable("mxdq");
    set_f_variable("mxdq" , 100.);

    dstep = 1.;
    step = 0;
    stpmin = 0.;
    imin = 0;
    alpha[0] = 0;
    fval[0] = 0;
    /*	 if( damp < 1.) {dstep = .25;}else{ dstep = .01/sqrt(damp);} */
if( damp < 1.) {dstep = .25;}else{ dstep = 1./damp;}

	natom = a_number();
	sumdx = a_l2_d();
	sumdx *=100.;
	

    for( i=0; i< nfs ; i++)
        (*ffs[i])(&fval[0],0.);
    fmin = fval[0]; fold = fmin;
    /* recheck is to find a descent stp length */
    imin++;
recheck:
    alpha[imin] = dstep;
    fval[imin] = 0.;
    for( i=0; i< nfs ; i++)
        (*ffs[i])(&fval[imin],dstep);
	fval[imin] += sumdx*dstep*dstep;
    if( fval[imin] > fval[0])
    {
        dstep = dstep*.25;
        if( dstep > 1.e-8) goto recheck;
        goto cleanup;
    }
    /* if here we have dstep small enough to use */
    /*	for( iter=0; iter< 8; iter++)
    */
    for( iter=0; iter< nostep; iter++)
    {
        for(jter=1; jter<100; jter++)
        {
            step += dstep;
            /* get the function value */
            for( i=0; i< imin; i++)
            {
                if( alpha[i] == step)
                {
                    fvt = fval[i]; goto FOUND ;
                }
            }
            fvt = 0.;
            for( i=0; i< nfs ; i++)
                (*ffs[i])(&fvt,step);
	fvt += sumdx*step*step;
            alpha[imin] = step;  fval[imin++] = fvt;
            if( imin > 400) imin = 400;
FOUND:
            /*		printf("\n %f %f %f %f", step,fvt,fmin,fold);  */
            if( fvt< fmin)
            {
                fmin = fvt; stpmin = step;
            }
            if( fvt > fold)
            {
                dstep = -dstep/2; break;
            }
            fold = fvt;
        }

    }
cleanup:
    /* do not allow 0 steps */
    /*	if( stpmin == 0.) stpmin = dstep/2; */
    set_f_variable("mxdq" , mxdq);
    return(stpmin);
}

