/* routines to read and contour and display density
*
*  we'll use XPLOR format
*  
*  v0.0 will expect 90 90 90 grids.  This is not complete 
*       but ALL ammp xplor maps are 90 90 90
*
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gtk/gtk.h>
#include <gdk/gdk.h>
#include "ammp.h"

typedef  struct
{
float xo,yo,zo; /* the coordinates of the point*/
int ixo,iyo;  /* transformed copies in pixels */
float txo,tyo,tzo; /* tranformed copy */
}  CONTOUR_POINT ;

#define MAX_BUNCH 10000
typedef struct
{
 CONTOUR_POINT c[MAX_BUNCH];
GdkColor *pen;
int inme;
void *next;
}  BUNCH ;

static BUNCH *first_bunch  = NULL;
static BUNCH *last_bunch  = NULL;


void DENSITY_free_density(  GdkColor *pen)
{
	BUNCH *bp,*bp_old,*bp_new;
	bp = first_bunch;
	if( bp == NULL) return;
	if( bp->next == bp)/* only one level */
	{
		if( bp->pen == pen)
		{ free(bp); first_bunch = NULL; last_bunch = NULL; return;}
	}
	if( bp->pen == pen)/* free the first */
	{
		first_bunch  = bp->next;
		free( bp);
	}
	bp_old = bp; bp = bp->next; bp_new = bp->next;
	while( 1==1)
	{
	if( bp->next == bp) break;
	if( bp == NULL) return; /* error condition */
	if( bp->pen == pen)
	{
		bp_old->next = bp_new;
		free(bp);
		bp = bp_new;
	}
	bp_old = bp; bp = bp->next; bp_new = bp->next;
	}
	if( bp->pen == pen) /* free the last */
	{
		last_bunch=bp_old;
		bp_old->next = bp_old;
		free(bp);
	}
}
void DENSITY_scatter_scan_density( FILE *ip,float level, GdkColor *pen, int np)
{
	int ix,imin,imax,jx,jmin,jmax,kx,kmin,kmax;
	float dx,dy,dz,alpha,beta,gamma;
	float xorig,yorig,zorig, xstep,ystep,zstep;
	float x,y,z;
	char type[4];
	int i,j,k;
	int isec,irow;
	BUNCH *bunch,*bpo;
	float (*sec)[];
	float randf();

	char line[80];
	
	if( np > MAX_BUNCH) np = MAX_BUNCH;


/* malloc an array and keep it */
	bunch = malloc( sizeof(BUNCH));
	if( bunch == NULL) 
		{aaerror("cannot allocate memory for contour");return;}
	if( first_bunch  == NULL){ first_bunch = bunch; last_bunch=bunch;}
	bpo = last_bunch;
	last_bunch->next = bunch;
	bunch->next = bunch;
/* skip the header */
	fgets( line,80,ip); printf("%s\n",line);
	fgets( line,80,ip); printf("%s\n",line);
	fgets( line,80,ip); printf("%s\n",line);
	fgets( line,80,ip); printf("%s\n",line);
/* get the grid info */
	fgets(line,80,ip); printf("%s\n",line);
	sscanf(line,"%8i%8i%8i%8i%8i%8i%8i%8i%8i",
	&ix,&imin,&imax,&jx,&jmin,&jmax,&kx,&kmin,&kmax);
	fgets(line,80,ip); printf("%s\n",line);
	sscanf(line,"%12.5e%12.5e%12.5e%12.5e%12.5e", &dx,&dy,&dz,&alpha,&beta,&gamma);
/* the type */
	fgets(line,80,ip); printf("%s\n",line);

	xstep = dx/ix;
	ystep = dy/jx;
	zstep = dz/kx;
	xorig = imin *xstep;
	yorig = jmin *ystep;
	zorig = kmin *zstep;

/* assuming ZYX maps because we're in a hurry */
	sec = malloc( ix*jx*kx*sizeof(float));
	if( sec == NULL  )
	{
		aaerror("Cannot allocate memory for contors");
		return;
	}

/* read the whole thing */
	for( k=0; k< kx; k++)
	{
	fgets(line,80,ip);
	sscanf(line,"%d",&isec);
	for(j=0; j< jx; j++)
	for(i=0; i< ix; i++)
	{
/* this occurs in rows of 6 plus remainder, I'm hopping fscanf is clever
enough */
		fscanf(ip,"%12.5e",&(*sec)[(k*jx+j)*ix + i]);
	}
	}
/* if we get here we have the whole damn thing in the core */
	bunch->inme = 0;
	bunch->pen = pen;
	for( k=0; k< kx-1; k++)
	for( j=0; j< jx-1; j++)
	for( i=0; i< ix-1; i++)
	{
		x = (*sec)[ (k*jx+j)*ix+i];
		y = (*sec)[ (k*jx+j)*ix+i +1];
		if( (x > level && y  < level ) || (x< level && y > level))
		{
	bunch->c[bunch->inme].xo = i*xstep + xorig+ (level-x)/(y-x)*xstep;
		bunch->c[bunch->inme].yo = j*ystep + yorig;
		bunch->c[bunch->inme].zo = k*zstep + zorig;
		bunch->inme++;
		if( bunch->inme == MAX_BUNCH) break;
		}
		x = (*sec)[ (k*jx+j)*ix+i];
		y = (*sec)[ (k*jx+j)*ix+i+ix];
		if( (x > level && y  < level ) || (x< level && y > level))
		{
		bunch->c[bunch->inme].xo = i*xstep + xorig;
	bunch->c[bunch->inme].yo = j*ystep + yorig+ (level-x)/(y-x)*ystep;
		bunch->c[bunch->inme].zo = k*zstep + zorig;
		bunch->inme++;
		if( bunch->inme == MAX_BUNCH) break;
		}
		x = (*sec)[ (k*jx+j)*ix+i];
		y = (*sec)[ (k*jx+j+jx)*ix+i];
		if( (x > level && y  < level ) || (x< level && y > level))
		{
		bunch->c[bunch->inme].xo = i*xstep + xorig;
		bunch->c[bunch->inme].yo = j*ystep + yorig;
	bunch->c[bunch->inme].zo = k*zstep + zorig+ (level-x)/(y-x)*zstep;
		bunch->inme++;
		if( bunch->inme == MAX_BUNCH) break;
		}
	}
	free(sec);
}/* end of DENSITY_scan_density */
