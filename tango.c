/* tango  (Torsion ANgle Global Optimizer)
*
*  basically reads 4 atom numbers then 10 floats that describe the distribution of
*  torsion angles for the 4 atoms.
*  tango terms are not directly a potential (well may convert them?) but are useful
*
*  the individual tango operations then use this to search the feasible space of the 
*  coordinates.
*/

#include "ammp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct{
	ATOM *a1,*a2,*a3,*a4;
	float raw[10];
	float cumulative[10];
	float current;  /* where to put the current torsion value if needed*/
	float sigma; /* an information status */
	void *next;
	}  TANGO;
#define TLONG sizeof(TANGO)

TANGO *tango_first = NULL;
TANGO *tango_last = NULL;
static int number_of_tangos = 0;
/* function tango sets up the tango database
* returns 1 if ok
* returns 1 if not
*
* undefined atoms are not lethal, but the term is not
* initialized if all the atoms are not defined
*
*
*  is passed the atom serial numbers and data
* allocates the new memory, initializes it and
* returns
*/
int tango( p1,p2,p3,p4,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9)
	int p1,p2,p3,p4;
	float k0,k1,k2,k3,k4,k5,k6,k7,k8,k9;
	{
	ATOM *ap1,*ap2,*ap3,*ap4,*a_m_serial();
	TANGO *new;
	int i;
	float s,s2;
	char line[80];
/* get the atom pointers for the two serial numbers */
	ap1 = a_m_serial( p1 );
	ap2 = a_m_serial( p2 );
	ap3 = a_m_serial( p3 );
	ap4 = a_m_serial( p4 );
	if( (ap1 == NULL) || (ap2 == NULL) || (ap3==NULL) || (ap4 == NULL) ) 
	{
/*	sprintf( line,"undefined atom in tango %d %d %d %d \0",p1,p2,p3,p4);
	aaerror( line ); */
	return 1;
	}

	if( ( new = malloc( TLONG ) ) == NULL)
	{ /* not being able to allocate memory is lethal */
	return 0;
	}
	/* initialize the pointers */
	if( tango_first == NULL) tango_first = new;
	if( tango_last == NULL) tango_last = new;
    number_of_tangos++;
	new -> a1 = ap1;
	new -> a2 = ap2;
	new -> a3 = ap3;
	new -> a4 = ap4;
	new -> raw[0] = k0;
	new -> raw[1] = k1;
	new -> raw[2] = k2;
	new -> raw[3] = k3;
	new -> raw[4] = k4;
	new -> raw[5] = k5;
	new -> raw[6] = k6;
	new -> raw[7] = k7;
	new -> raw[8] = k8;
	new -> raw[9] = k9;
	new -> next = new;

	new -> cumulative[0] = k0;
	/* while they are expected to be normalized input they might not be */
	for( i=1;i< 10; i++)
	      new -> cumulative[i] = new->cumulative[i-1] + new->raw[i];

	if( new->cumulative[9] < 0.00001)
	{
		for( i=0; i< 10; i++)
			new->cumulative[i] = (float)(i+1)*0.1;
	} else {
		for( i=0; i< 10; i++)
			new->cumulative[i] /= new->cumulative[9];
	}
	new->sigma = 0.;
	s2 = 0.; s = 0.;
	for( i = 0; i< 10; i++)
		  { s += new->raw[i]; s2 += new->raw[i]*new->raw[i];}
	new->sigma = s2/10. - s*s/100.;
           
	tango_last -> next = new;
	tango_last = new;
	return 1;
	}

	TANGO *t_iterator(i)
		int i;
	{
		static TANGO *t;
		if( i <= 0)
			t = tango_first;
		else
		{
//			if( t == tango_last) return NULL;
			if( t == t->next) return NULL;
			t = t->next;
		}
		return t;
	}



int clone( int);
int restore_from_clone(int);
float randf();
float randg();

float get_torsion_value( ATOM*, ATOM*, ATOM*,ATOM*);


int v_tango( V, dl)
float *V;
float dl;
{
	int i;
	TANGO *t;
	i = 0;
	while ( 1 ==1 )
		{
			int ir;
			float tv,lv; 
			t = t_iterator(i++);
			if( t == NULL) break;
			tv = get_torsion_value( t->a1,t->a2,t->a3,t->a4);
			while( tv < 0.) tv += 2.*3.14159265;
			while( tv > 2.*3.14159265) tv -= 2.*3.14159265;
			ir = (int)( 10.*tv /(2.*3.14159265));
			if( ir > 9) ir = 9;
			if( ir < 0) ir = 0;
			lv = t->raw[ir];
			lv += 1.e-3;  // bias off zero 
			lv = 10.*lv;   // correction  for mean
			*V += 300.*0.2*log( lv);

		}
	
return 1;
}
// dummy to keep eval lists of potential/forces happy
// since this is  a table there is no force
int f_tango( dl)
float dl;
{
	return 1;
}


int tango_search( op,vfs, ffs, nfs, how_many)
FILE *op;
int nfs, how_many;
int (*vfs[])();
int (*ffs[])();
{
	float V,Vbest,Vold;
	int i,iterations;
	float rv;
	float dstep, T, H;
	TANGO *t;

	
	/*  useful code snippets

 	clone( 1);  
	clone( 2);  

	restore_from_clone(1);

	V = 0;
	for( i=0; i< nfs; i++)
		(*vfs[i])(&V,0.);
		*/


	clone(1);
    clone(2);
	Vbest = 0.;
	for( i =0 ; i< nfs; i++)
	     (*vfs[i])(&Vbest, 0.);
	Vold = Vbest;
	V = Vbest;

     dstep = 0.7;
	 T = Vbest * 0.1; 
	 H = Vbest;
	 dstep = 0.1;
	for( iterations = 0; iterations < how_many; iterations++)
	{
		int local_iteration;
		float vbar;
		int best_found;

		vbar = 0.;
		best_found = 0;
//	dstep -= 0.01;
//	if( dstep < 0.1) dstep = 0.1;
		for( local_iteration =0; local_iteration < 1; local_iteration++)
		{
			int who_to_do;

			who_to_do = randf()* number_of_tangos;
			if(who_to_do < 0) who_to_do = 0;
			if( who_to_do >= number_of_tangos) who_to_do = number_of_tangos-1;
		i = 0;
//		if( V < H+T && V > H-T)
	//		clone(2);
		while ( 1 ==1 )
		{
			int ir;
			float tv; 
			float original,delta;
			float local_best;
			int ilocal_best;
			t = t_iterator(i++);
			if( t == NULL) break;
		//	rv = randf();  
	      //  if( rv > 0.01 ) continue;

     		if( i-1 != who_to_do) continue;
//#define pure
#define cannonical

			tset_bond_build( t->a1, t->a2, t->a3, t->a4);
		//	set_torsion( t->a1 , t-> a2, t->a3, t->a4, 0. );

	
		    for( ilocal_best = 0; ilocal_best < 1000; ilocal_best++)
	//		for( ilocal_best = 0; ilocal_best < 1; ilocal_best++)
			{
				float ltv;
			     int index;
			tv = get_torsion_value( t->a1, t->a2, t->a3, t->a4);
			//	if( t->sigma < 0.01)
			//	 delta = dstep*randg();
			//	else
			//	 delta = randf()*2.* 3.14159265;
			delta = 0.1*randg();	
				 rv = randf();
					ltv = tv + delta;
					while( ltv > 2.*3.14159265) ltv -= 2.*3.14159265;
					while( ltv < 0.) ltv += 2.*3.14159265;
					index = (int)(9.*ltv/(2.*3.14159265));

					if( index < 0) index = 0;
					if( index > 9) index = 9;
					if( rv > t->raw[index]) continue;

			set_torsion( t->a1,t->a2,t->a3,t->a4,delta );
#ifdef cannonical
			V = 0;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);

				if(V < H+T && V >H-T) break;
	//		break;
#endif	
#ifdef pure
				break;
#endif		
			}

//			break;

         
//		force_screen_update();
		}
#ifdef pure
//				cngdel( vfs,ffs,nfs, 100,100,0., (1==0));
#endif
//		force_screen_update();
#ifdef windoze
		BE_NICE();
#endif
		V = 0;
		for( i=0; i< nfs; i++)
			(*vfs[i])(&V,0.);
#ifdef windoze
		BE_NICE();
#endif
		fprintf(op,"%f %f %f \n", V, H,Vbest );
		vbar += V;
#ifdef pure
		if( V < Vbest) {best_found ++; clone(1);  clone(2); Vbest = V; Vold = V;}
		else if( V < Vold) { Vold = V; clone(2);}
		else {
			if( randf() > exp( -(V-Vold)/(T)))
			{ /*Vold = Vbest; */ restore_from_clone(2); /*force_screen_update();*/}
			else { Vold = V;}
		}
#endif	
#ifdef cannonical
//1		if( V < H+T && V > H-T)
//1		{
			if( V < Vbest) {clone(1); Vbest = V; Vold = V;}
//			if( V < Vold ) { clone(2); Vold = V;}
//1		}
		 if( V > H+T  ||V < H-T){ 
			 restore_from_clone(2); /*force_screen_update()*/;} else{ clone(2);}
#endif		
		} /* local_iteration */



		H = Vbest - T*0.95;
		T = Vbest * 0.1;
		if( T < 100.) T = 100.;


	}/* iterations*/
	restore_from_clone(1);
	return 1;
}
