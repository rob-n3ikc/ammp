/* this draws the molecule
*
*/
#include <stdio.h>
#include <math.h>
#include <gtk/gtk.h>
#include <gdk/gdk.h>

#include "ammp.h"

/* should be passed at the signal connect time (but isn't always)*/
extern GtkWidget *AMMP_draw_area;
extern int AMMP_READY_go_for_it;
extern int AMMP_NOOP_im_busy;
static int AMMP_DRAW_im_busy = 1==0;

/* the pixmap which is drawn */
static GdkPixmap *pix = NULL;
static GdkGC *Gcontext = NULL;
static GdkColor redpen ;
static GdkColor yellowpen;
static GdkColor greenpen;
static GdkColor bluepen;
static GdkColor blackpen;
static GdkColor greypen;
static GdkColor whitepen;
static GdkColor purplepen;
static GdkColor anypen;
static GdkColor redgradpen[64];
static GdkFont *AMMPfont = NULL;

typedef struct { int inme;  GtkWidget *bunch[10];} BUNCH_OF_WIDGETS;

static ATOM *(picked_queue)[10] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

gint set_anypen( GtkWidget *the_widget, float r,float g,float b)
{
	if( r < 0.) r = 0.; if( r > 1. ) r = 1.;
	if( g < 0.) g = 0.; if( g > 1. ) g = 1.;
	if( b < 0.) b = 0.; if( b > 1. ) b = 1.;
	anypen.red = 65535*r;
	anypen.green = 65535*g;
	anypen.blue = 65535*b;
/*
 1.2 only?
	gdk_color_change( gtk_widget_get_colormap(the_widget),&anypen);
*/
	anypen.pixel = (gulong)(anypen.red*65535*65535 + anypen.green*65535 + anypen.blue);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&anypen);

	return TRUE;
}

gint configure_pens(GtkWidget *the_widget)
{
	int i;
	redpen.red = 65535; redpen.blue = 0; redpen.green = 0; redpen.pixel = (gulong)(redpen.red*65536*65536 + redpen.green*65536 + redpen.blue);
	yellowpen.red = 65535; yellowpen.blue = 0; yellowpen.green = 65535; yellowpen.pixel = (gulong)(yellowpen.red*65536*65536 + yellowpen.green*65536 + yellowpen.blue);
	greenpen.red = 0; greenpen.blue = 0; greenpen.green = 65535; greenpen.pixel = (gulong)(greenpen.red*65536*65536 + greenpen.green*65536 + greenpen.blue);
	bluepen.red = 0; bluepen.blue = 65535; bluepen.green = 0; bluepen.pixel = (gulong)(bluepen.red*65536*65536 + bluepen.green*65536 + bluepen.blue);
	blackpen.red = 0; blackpen.blue = 0; blackpen.green = 0; blackpen.pixel = (gulong)(blackpen.red*65536*65536 + blackpen.green*65536 + blackpen.blue);
	whitepen.red = 65535; whitepen.blue = 65535; whitepen.green = 65535; whitepen.pixel = (gulong)(whitepen.red*65536*65536 + whitepen.green*65536 + whitepen.blue);
	greypen.red = 65535*0.65; greypen.blue = 65535*0.65; greypen.green = 65535*0.65; greypen.pixel = (gulong)(greypen.red*65536*65536 + greypen.green*65536 + greypen.blue);
	anypen.red = 65535; anypen.green = 65535; anypen.blue = 65535; anypen.pixel = (gulong) 1;
	purplepen.red = 65535/2; purplepen.blue = 65535; purplepen.green = 0;
	purplepen.pixel = (gulong)(purplepen.red*65535*65535 + purplepen.green*65535 + purplepen.blue);

	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&redpen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&yellowpen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&greenpen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&bluepen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&purplepen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&blackpen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&whitepen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&greypen);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&anypen);
	for( i=0; i< 64; i++)
	{
	redgradpen[i].red = 65535;	
	redgradpen[i].green = 65535*i/63;	
	redgradpen[i].blue = 65535*i/63;	
	redgradpen[i].pixel = (gulong)( redgradpen[i].red*65535*65535 +
		redgradpen[i].green*65535 + redgradpen[i].blue);
	gdk_color_alloc( gtk_widget_get_colormap(the_widget),&redgradpen[i]);
	}
	return TRUE;
}
#include "ammp.xpm" /* initial screen */
gint configure_event_handler( GtkWidget *wind, GdkEventConfigure *event)
{
	gint draw_the_molecule();
	GtkStyle *style;
	if( pix != NULL)
		gdk_pixmap_unref(pix);
	pix = gdk_pixmap_new( wind->window, 
		wind->allocation.width,wind->allocation.height,-1);
	if( Gcontext  == NULL){ 
	Gcontext = gdk_gc_new( wind->window); 
	style =gtk_widget_get_default_style();
	AMMPfont = style->font;
	if( AMMPfont == NULL ) printf(" font not loaded\n");
	configure_pens( wind);
/* fill it up with a test pattern */
/*
	gdk_gc_set_background( Gcontext,&greypen);
	gdk_gc_set_foreground( Gcontext,&greypen);
	gdk_draw_rectangle( pix,Gcontext,TRUE,0,0, wind->allocation.width,wind->allocation.height);
	gdk_gc_set_foreground( Gcontext,&redpen);
	gdk_draw_rectangle( pix, Gcontext,TRUE,10,10,30,30);				
	gdk_gc_set_foreground( Gcontext,&greenpen);
	gdk_draw_rectangle( pix, Gcontext,TRUE,40,10,30,30);				
	gdk_gc_set_foreground( Gcontext,&bluepen);
	gdk_draw_rectangle( pix, Gcontext,TRUE,70,10,30,30);				
	gdk_gc_set_foreground( Gcontext,&blackpen);
	gdk_draw_rectangle( pix, Gcontext,TRUE,100,10,30,30);				
	gdk_gc_set_foreground( Gcontext,&whitepen);
	gdk_draw_rectangle( pix, Gcontext,TRUE,130,10,30,30);				
	gdk_gc_set_foreground( Gcontext,&yellowpen);
	gdk_draw_rectangle( pix, Gcontext,TRUE,160,10,30,30);				
*/
/* wind->window is the GDK window see gtk includes */
	pix = gdk_pixmap_create_from_xpm_d( wind->window, NULL, NULL, ammp_xpm);
	} else {
	gdk_gc_set_background( Gcontext,&greypen);
	gdk_gc_set_foreground( Gcontext,&greypen);
	gdk_draw_rectangle( pix,Gcontext,TRUE,0,0, wind->allocation.width,wind->allocation.height);
	}
	return TRUE;
}

gint expose_event_handler( GtkWidget *wind, GdkEventExpose *event)
{
	gint draw_the_molecule();

	gdk_draw_pixmap( wind->window,
	Gcontext,
	pix,
	event->area.x, event->area.y,
	event->area.x, event->area.y,
	event->area.width,event->area.height);
	return 1==1;
}


static int  do_scale=TRUE;
static float AMMPcx,AMMPcy,AMMPcz;
static float AMMPslab=100.;
static float rotate[3][3] = {1.,0.,0., 0., 1.,0., 0., 0., 1.};
static float scale = 1.;
static int orient_flag = 0; 
static int center_flag = 0;
static int noel_flag = (1==0);
static int noel_violation_flag = (1==0);
static int tether_flag= (1==0);
static int draw_color_type = 1;
#ifndef PI
#define PI 3.14159265
#endif
gint set_orientation( GtkWidget *wind, gpointer *what)
{
	char *cp;
	gint draw_the_molecule();
	cp = (char *)what;
	if( *cp == 'X') {
	rotate[0][0] = 0.; rotate[0][1] = 1.; rotate[0][2] = 0.;
	rotate[1][0] = 0.; rotate[1][1] = 0.; rotate[1][2] = 1.;
	rotate[2][0] = 1.; rotate[2][1] = 0.; rotate[2][2] = 0.;
	}
	if( *cp == 'Y') {
	rotate[0][0] = 0.; rotate[0][1] = 0.; rotate[0][2] = 1.;
	rotate[1][0] = 1.; rotate[1][1] = 0.; rotate[1][2] = 0.;
	rotate[2][0] = 0.; rotate[2][1] = 1.; rotate[2][2] = 0.;
	}
	if( *cp == 'Z') {
	rotate[0][0] = 1.; rotate[0][1] = 0.; rotate[0][2] = 0.;
	rotate[1][0] = 0.; rotate[1][1] = 1.; rotate[1][2] = 0.;
	rotate[2][0] = 0.; rotate[2][1] = 0.; rotate[2][2] = 1.;
	}
	draw_the_molecule( AMMP_draw_area,NULL);
	return 1==1;
}

gint clear_picked_list( GtkWidget *wind, gpointer *dummy)
{
int i;
ATOM *ap,*a_first();
	for( i=0; i< 10; i++)
	picked_queue[i] = NULL;

	ap = a_first();
	while(1 == 1)
	{
		ap->is_label = (1==0);
		if( ap == ap->next) break;
		ap = ap->next;
		if( ap == NULL) break;
	}

}
gint AMMP_analyze_handler( GtkWidget *wind, int (*func)())
{
	int (*f)();
	int i;
	float x;
	ATOM *ap;

	f = func;
	ap = picked_queue[0];
	if( ap == NULL ) return 1==1;
	i = ap->serial;
	(*f)( &x,0., i,i,stdout);	
	return 1==1;

}
gint set_draw_color_type( GtkWidget *wind, gpointer *flag)
{
	draw_color_type = (int)*(char *)flag;
}
gint toggle_tether_flag( GtkWidget *wind, gpointer *dummy)
{   if( !tether_flag)  tether_flag = (1==1); 
	else tether_flag = (1==0);  return 1==1;}

gint toggle_noel_violation_flag( GtkWidget *wind, gpointer *dummy)
{   if( !noel_violation_flag)  noel_violation_flag = (1==1); 
	else noel_violation_flag = (1==0);  return 1==1;}
gint toggle_noel_flag( GtkWidget *wind, gpointer *dummy)
{   if( !noel_flag) noel_flag = (1==1); 
	else noel_flag = (1==0); return 1==1;}

gint set_orient_flag_xy( GtkWidget *wind, gpointer *dummy)
{   orient_flag = 0; return 1==1;}
gint set_orient_flag_yz( GtkWidget *wind, gpointer *dummy)
{   orient_flag = 1; return 1==1;}
gint set_orient_flag_scale( GtkWidget *wind, gpointer *dummy)
{   orient_flag = 2; return 1==1;}
gint set_center_flag_mean( GtkWidget *wind, gpointer *dummy)
{   center_flag = 2; return 1==1;}
gint set_center_flag_atom( GtkWidget *wind, gpointer *dummy)
{   center_flag = 1; return 1==1;}
gint set_center_flag_track( GtkWidget *wind, gpointer *dummy)
{   center_flag = 0; return 1==1;}
gint set_scale_flag( GtkWidget *wind, gpointer *dummy)
{ if( do_scale) do_scale = (1==0); else do_scale  = (1==1); return 1==1;}

void draw_matrix_update( GtkWidget *wind, int x,int y ,int flag)
{
	gint draw_the_molecule();
	static int oldx = -100;
	static int oldy ;
	float m[3][3],m2[3][3];
	float wtor,htor;
	float dx,dy;
	float cx,cy,sx,sy;

	if( orient_flag < 0 || orient_flag > 2) return ;

	if(oldx < 0) { /* avoid an initial jump */
	oldx = x; oldy = y;
	return ;
	}
	wtor = wind->allocation.width;
	htor = wind->allocation.height;
	if( orient_flag == 2)
	{
	dx = (x - oldx); dy = (y -oldy);
	oldx = x; oldy = y;
	if( !do_scale) scale = scale * exp( dx/wtor);
	AMMPslab = (y/htor)*(y/htor)*100.;
	if( AMMPslab < .1) AMMPslab = .1;
	draw_the_molecule( wind,NULL);
	return ;
	}
	wtor = 2.*PI/wtor;
	htor = 2.*PI/htor;
	dy = wtor*( x - oldx); dx = htor*(y -oldy);
	oldx = x; oldy = y;
	if( !flag) return;

	cx = cos(dx); sx = sin(dx);
	cy = cos(dy); sy = sin(dy);
/* build the rotation matrix */
/* put the cases for x,y y,z here */
/* pre-multiplied y,x */
	if( orient_flag  == 0){
	m[0][0] = cy;     m[0][1] = 0.; m[0][2] = -sy;
	m[1][0] = -sx*sy; m[1][1] = cx; m[1][2] = -sx*cy;
	m[2][0] = cx*sy;  m[2][1] = sx; m[2][2] = cx*cy;
		}
	if( orient_flag  == 1){
	m[0][0] = cx*cy;    m[0][1] = -sx; m[0][2] = -cx*sy;
	m[1][0] = sx*cy;    m[1][1] = cx;  m[1][2] = -sx*sy;
	m[2][0] = sy;       m[2][1] = 0.;  m[2][2] = cy;
		}
/* multiply by the rotation matrix (again write out the mults for speed) */
	m2[0][0] = m[0][0]*rotate[0][0] +m[0][1]*rotate[1][0] + m[0][2]*rotate[2][0];
	m2[0][1] = m[1][0]*rotate[0][0] +m[1][1]*rotate[1][0] + m[1][2]*rotate[2][0];
	m2[0][2] = m[2][0]*rotate[0][0] +m[2][1]*rotate[1][0] + m[2][2]*rotate[2][0];
	m2[1][0] = m[0][0]*rotate[0][1] +m[0][1]*rotate[1][1] + m[0][2]*rotate[2][1];
	m2[1][1] = m[1][0]*rotate[0][1] +m[1][1]*rotate[1][1] + m[1][2]*rotate[2][1];
	m2[1][2] = m[2][0]*rotate[0][1] +m[2][1]*rotate[1][1] + m[2][2]*rotate[2][1];
	m2[2][0] = m[0][0]*rotate[0][2] +m[0][1]*rotate[1][2] + m[0][2]*rotate[2][2];
	m2[2][1] = m[1][0]*rotate[0][2] +m[1][1]*rotate[1][2] + m[1][2]*rotate[2][2];
	m2[2][2] = m[2][0]*rotate[0][2] +m[2][1]*rotate[1][2] + m[2][2]*rotate[2][2];
/* and graham schmit it */
	m[0][0] = m2[0][0]*m2[0][0] + m2[0][1]*m2[0][1] + m2[0][2]*m2[0][2] ;
	m[0][1] = m2[1][0]*m2[1][0] + m2[1][1]*m2[1][1] + m2[1][2]*m2[1][2] ;
	m[0][2] = m2[2][0]*m2[2][0] + m2[2][1]*m2[2][1] + m2[2][2]*m2[2][2] ;
	dx = 1./sqrt(m[0][0] );
	sx = 1./sqrt( m[0][0]*m[0][1]);
	sy = 1./sqrt( m[0][0]*m[0][2]);
	m[0][1] = m2[0][0]*m2[1][0] + m2[0][1]*m2[1][1] + m2[0][2]*m2[1][2] ;
	m[0][2] = m2[0][0]*m2[2][0] + m2[0][1]*m2[2][1] + m2[0][2]*m2[2][2] ;
	
	rotate[0][0] = m2[0][0]*dx;
	rotate[1][0] = m2[0][1]*dx;
	rotate[2][0] = m2[0][2]*dx;
	m2[1][0] = m2[1][0] - m[0][1]*sx*m2[0][0];
	m2[1][1] = m2[1][1] - m[0][1]*sx*m2[0][1];
	m2[1][2] = m2[1][2] - m[0][1]*sx*m2[0][2];
	m2[2][0] = m2[2][0] - m[0][2]*sy*m2[0][0];
	m2[2][1] = m2[2][1] - m[0][2]*sy*m2[0][1];
	m2[2][2] = m2[2][2] - m[0][2]*sy*m2[0][2];
	m[0][1] = m2[1][0]*m2[1][0] + m2[1][1]*m2[1][1] + m2[1][2]*m2[1][2] ;
	m[0][2] = m2[2][0]*m2[2][0] + m2[2][1]*m2[2][1] + m2[2][2]*m2[2][2] ;
	dx = 1./sqrt(m[0][1] );
	sx = 1./sqrt( m[0][1]*m[0][2]);
	rotate[0][1] = m2[1][0]*dx;
	rotate[1][1] = m2[1][1]*dx;
	rotate[2][1] = m2[1][2]*dx;
	m[0][2] = m2[1][0]*m2[2][0] + m2[1][1]*m2[2][1] + m2[1][2]*m2[2][2] ;
	m2[2][0] = m2[2][0] - m[0][2]*sx*m2[1][0];
	m2[2][1] = m2[2][1] - m[0][2]*sx*m2[1][1];
	m2[2][2] = m2[2][2] - m[0][2]*sx*m2[1][2];
	m[0][2] = m2[2][0]*m2[2][0] + m2[2][1]*m2[2][1] + m2[2][2]*m2[2][2] ;
	dx = 1./sqrt(m[0][2] );
	rotate[0][2] = m2[2][0]*dx;
	rotate[1][2] = m2[2][1]*dx;
	rotate[2][2] = m2[2][2]*dx;
	draw_the_molecule( wind,NULL);

}
gint draw_the_molecule( GtkWidget *wind ,gpointer *dummy)
{
	int natoms,a_number();
	ATOM *ap,*a_first();
	GdkRectangle rect;
	

	void center_the_molecule();
	void scale_the_molecule();
	void atom_walker();
	void draw_bond_list();
	void draw_angle_list();
	void draw_torsion_list();
	void draw_noel();
	void draw_tether();
	void bond_walker();
	void bond_atom_color_walker();
	void bond_force_color_walker();
	void block_unblock();

	natoms = a_number();
	if( natoms < 1) return;
/*	if( !AMMP_READY_go_for_it ) return ;  don't go if the atom structure is mucky */
	if( AMMP_DRAW_im_busy) return;
	AMMP_DRAW_im_busy = 1==1;
	block_unblock(TRUE);

	if( pix != NULL)
		gdk_pixmap_unref(pix);
	pix = gdk_pixmap_new( wind->window, 
		wind->allocation.width,wind->allocation.height,-1);
	gdk_gc_set_background( Gcontext,&greypen);
	gdk_gc_set_foreground( Gcontext,&greypen);
	gdk_draw_rectangle( pix,Gcontext,TRUE,0,0, wind->allocation.width,wind->allocation.height);

/* scale the beast first */
	center_the_molecule(natoms);

	scale_the_molecule( wind, natoms);

/* walk the atoms */
	atom_walker(pix);
/* walk the distance list */
	draw_bond_list(pix);
	draw_angle_list(pix);
	draw_torsion_list(pix);
	if( noel_flag) draw_noel();
	if( tether_flag) draw_tether(wind);
/* walk the bonds */
/* monochrome or polycolor gets black background, all others white */
	if( draw_color_type == 0 || draw_color_type > 2)	
	bond_walker(pix,&blackpen);

/* CPK colors */
	if( draw_color_type == 1)
	bond_atom_color_walker( wind);
/* force colors */
	if( draw_color_type == 2)
	bond_force_color_walker( wind);

/* and draw the map */
	rect.x = 0;
	rect.y = 0;
	rect.width = wind->allocation.width;
	rect.height = wind->allocation.height;
	gtk_widget_draw(wind,&rect);
	block_unblock(FALSE);
	AMMP_DRAW_im_busy = 0==1;

}
void center_the_molecule(int na)
{
	ATOM *ap,*a_first();
	float inv_na;

	if( center_flag == 0){
	AMMPcx= 0.; AMMPcy = 0.; AMMPcz = 0.;
	ap = a_first();
	while( 1==1)
	{
	AMMPcx += ap->x;	
	AMMPcy += ap->y;	
	AMMPcz += ap->z;	
	if( ap == ap->next) break;
	ap = ap->next;
	if( ap == NULL ) break;
	}
	inv_na = 1./(float)na;
	AMMPcx  *= inv_na;
	AMMPcy  *= inv_na;
	AMMPcz  *= inv_na;
	}
	if( center_flag  == 1)
	{
	if( picked_queue[0] == NULL ){center_flag = 0; return;}
	AMMPcx = picked_queue[0]->x;
	AMMPcy = picked_queue[0]->y;
	AMMPcz = picked_queue[0]->z;
	center_flag   = 3;
	}
	if( center_flag   == 2)
	{
	AMMPcx= 0.; AMMPcy = 0.; AMMPcz = 0.;
	ap = a_first();
	while( 1==1)
	{
	AMMPcx += ap->x;	
	AMMPcy += ap->y;	
	AMMPcz += ap->z;	
	if( ap == ap->next) break;
	ap = ap->next;
	if( ap == NULL ) break;
	}
	inv_na = 1./(float)na;
	AMMPcx  *= inv_na;
	AMMPcy  *= inv_na;
	AMMPcz  *= inv_na;
	center_flag = 3;
	}
	
}

void scale_the_molecule( GtkWidget *window, int na)
{
	ATOM *ap,*a_first();
	float width,height;
	float xmax,xmin,ymax,ymin;	
	float x,y,z;
	
/* transform the molecule */
	xmax = -10e10;
	ymax = -10e10;
	xmin =  10e10;
	ymin =  10e10;
	ap = a_first();
	while(1==1)
	{
	x = ap->x - AMMPcx;
	y = ap->y - AMMPcy;
	z = ap->z - AMMPcz;
	ap->tx = rotate[0][0] *x + rotate[0][1]*y + rotate[0][2]*z;
	ap->ty = rotate[1][0] *x + rotate[1][1]*y + rotate[1][2]*z;
	ap->tz = rotate[2][0] *x + rotate[2][1]*y + rotate[2][2]*z;
	if( xmax < ap->tx) xmax = ap->tx;
	if( ymax < ap->ty) ymax = ap->ty;
	if( xmin > ap->tx) xmin = ap->tx;
	if( ymin > ap->ty) ymin = ap->ty;
	if( ap == ap->next) break;
	ap = ap ->next;
	if( ap == NULL) break;
	}

	width = (float)(window->allocation.width);
	height = (float)(window->allocation.height);
	if( do_scale)
	{
	x =  (xmax-xmin)/(0.9*width);
	y =  (ymax-ymin)/(0.9*height);
	if( y > x) scale = 1./y;
	if( x > y) scale = 1./x;
	if( x <= 1.e-7 && y <= 1.e-7) scale = 1.;
	}
	width *= 0.5;
	height *= 0.5;
	ap = a_first();
	while(1==1)
	{
	ap->ix = scale*ap->tx + width;
	ap->iy = scale*ap->ty + height;
	if( ap == ap->next) break;
	ap = ap ->next;
	if( ap == NULL) break;
	}
	
}
 

/* atom walker will walk the list and display labels */
void atom_walker( GtkWidget *where)
{
	ATOM *ap,*a_first();
	float upper,lower;
	char line[25];
	ap = a_first();
	upper =  AMMPslab;
	lower =  - AMMPslab;
	while( 1==1)
	{
/* check the slab now */
	ap->is_visible = (1==1);
	if( ap->tz < lower) ap->is_visible = (1==0);
	if( ap->tz > upper) ap->is_visible = (1==0);
/* draw a point */
	if( !ap->is_bonded){
	gdk_gc_set_foreground( Gcontext,&greenpen);
	gdk_draw_line(pix, Gcontext, ap->ix+2,ap->iy+2,ap->ix-2,ap->iy-2);
	gdk_draw_line(pix, Gcontext, ap->ix+2,ap->iy-2,ap->ix-2,ap->iy+2);
	}
/* check the label now */
	if( ap->is_visible && ap->is_label)
	{
		gdk_gc_set_foreground( Gcontext,&redpen);
		sprintf(line,"%d %s",ap->serial,ap->name);
		gdk_draw_string( pix, AMMPfont, Gcontext, ap->ix,ap->iy,line);
	}
	if( ap == ap->next) break;
	ap = ap->next;
	if( ap == NULL ) break;
	}
	return;
}
void bond_walker( GtkWidget *where, GdkColor *pen)
{
	ATOM *ap1,*ap2;
	int i;
	int bond_next(); /* in bonds.c */


	gdk_gc_set_foreground( Gcontext,pen);
	i = -1;
	while( bond_next(i,&ap1,&ap2) == 1)
	{
	if( ap1 == NULL) return ;
	if( ap2 == NULL) return ;
	i = 1;
/* now draw it */
	if( ap1->is_visible && ap2->is_visible)
	gdk_draw_line( pix,Gcontext, ap1->ix,ap1->iy, ap2->ix,ap2->iy);
	}
	if( ap1->is_visible && ap2->is_visible)
	gdk_draw_line( pix,Gcontext, ap1->ix,ap1->iy, ap2->ix,ap2->iy);

}/* end of bond_walker */

void bond_atom_color_walker( GtkWidget *where)
{
	ATOM *ap1,*ap2;
	int i,j,k;
	int bond_next(); /* in bonds.c */
	char *cp;
	GdkColor *pen;


/*	gdk_gc_set_foreground( Gcontext,&blackpen); */
	i = -1;
	while( 1 == 1)
	{
	k = bond_next(i,&ap1,&ap2);
	if( ap1 == NULL) return ;
	if( ap2 == NULL) return ;
	i = 1;
/* now draw it */
	if( ap1->is_visible && ap2->is_visible)
	{
	cp = ap1->name;
	for( j=0; j< 9; j++)
	{if( *cp == '.') {cp++; goto DOIT1;} if(*cp == '\0') goto SKIP1; cp++;}	
	goto SKIP1;
DOIT1: ;
	pen = &greenpen;
	if( *cp == 'h' || *cp < '9')  pen = &whitepen;
	if( *cp == 'c' && *(cp+1) != 'l')  pen=&blackpen;
	if( *cp == 'o' )  pen = &redpen;
	if( *cp == 'n' )  pen = &bluepen;
	if( *cp == 's' )  pen = &yellowpen;
	gdk_gc_set_foreground( Gcontext,pen);
	gdk_draw_line( pix,Gcontext, ap1->ix,ap1->iy, (ap2->ix+ap1->ix)/2,(ap2->iy+ap1->iy)/2);
SKIP1: ;
	cp = ap2->name;
	for( j=0; j< 9; j++)
	{if( *cp == '.') {cp++; goto DOIT2;} if(*cp == '\0') goto SKIP2; cp++;}	
	goto SKIP2;
DOIT2: ;
	pen = &greenpen;
	if( *cp == 'h' || *cp < '9')  pen = &whitepen;
	if( *cp == 'c' && *(cp+1) != 'l')  pen=&blackpen;
	if( *cp == 'o' )  pen = &redpen;
	if( *cp == 'n' )  pen = &bluepen;
	if( *cp == 's' )  pen = &yellowpen;
 	gdk_gc_set_foreground( Gcontext,pen);
	gdk_draw_line( pix,Gcontext, ap2->ix,ap2->iy, (ap2->ix+ap1->ix)/2,(ap2->iy+ap1->iy)/2);
SKIP2: ;
	}
	if( k != 1) break;
	}

}/* end of bond_atom_color_walker */


void bond_force_color_walker( GtkWidget *where)
{
	ATOM *ap1,*ap2;
	int i,j,k;
	float f;
	int bond_next(); /* in bonds.c */
	char *cp;
	GdkColor *pen;


/*	gdk_gc_set_foreground( Gcontext,&blackpen); */
	i = -1;
	while( 1 == 1)
	{
	k = bond_next(i,&ap1,&ap2);
	if( ap1 == NULL) return ;
	if( ap2 == NULL) return ;
	i = 1;
/* now draw it */
	if( ap1->is_visible && ap2->is_visible)
	{
	f = ap1->fx*ap1->fx;
	f += ap1->fy*ap1->fy;
	f += ap1->fz*ap1->fz;
	f = sqrt(f);
	j = 63*f*0.01;
	if( j > 63) j = 63;
	if( j < 0) j = 0;
	j = 63-j;
	gdk_gc_set_foreground( Gcontext,&redgradpen[j]); 
	gdk_draw_line( pix,Gcontext, ap1->ix,ap1->iy, (ap2->ix+ap1->ix)/2,(ap2->iy+ap1->iy)/2);
	f = ap2->fx*ap2->fx;
	f += ap2->fy*ap2->fy;
	f += ap2->fz*ap2->fz;
	f = sqrt(f);
	j = 63*f*0.01;
	if( j > 63) j = 63;
	if( j < 0) j = 0;
	j = 63-j;
	gdk_gc_set_foreground( Gcontext,&redgradpen[j]); 
	gdk_draw_line( pix,Gcontext, ap2->ix,ap2->iy, (ap2->ix+ap1->ix)/2,(ap2->iy+ap1->iy)/2);
	}
	if( k != 1) break;
	}

}/* end of bond_force_color_walker */


void  add_atom_to_picked_queue( ATOM *ap)
{
	int i; 
	if( ap == NULL) return;
	for( i=9; i> 0; i--)
		picked_queue[i] = picked_queue[i-1];
	picked_queue[0] = ap;
	return;
}
ATOM *get_atom_from_picked_queue( int i)
{
	if(i >=0 && i < 10) return picked_queue[i];
	return NULL;
}

/* these routines support distance lists */
typedef struct
{ATOM *a1,*a2;} ATOM_PAIR;
typedef struct
{ATOM *a1,*a2,*a3;} ATOM_TRIPLE;
typedef struct
{ATOM *a1,*a2,*a3,*a4;} ATOM_QUARTET;

static GSList *bondlist = NULL;
static GSList *anglelist = NULL;
static GSList *torsionlist = NULL;

gint add_pair_to_bond_list( GtkWidget *wind, gpointer *dummy)
{
	ATOM_PAIR *app;
	if( picked_queue[0] == NULL || picked_queue[1] == NULL) return 1==1;

	app = malloc( sizeof(ATOM_PAIR) );
	app->a1 = picked_queue[0];
	app->a2 = picked_queue[1];
	bondlist = g_slist_append( bondlist,app);
	draw_the_molecule( AMMP_draw_area,NULL);
	return 1==1;
} 
gint add_pair_to_angle_list( GtkWidget *wind, gpointer *dummy)
{
	ATOM_TRIPLE *app;
	if( picked_queue[0] == NULL || picked_queue[1] == NULL 
		|| picked_queue[2] == NULL) return 1==1;

	app = malloc( sizeof(ATOM_TRIPLE) );
	app->a1 = picked_queue[0];
	app->a2 = picked_queue[1];
	app->a3 = picked_queue[2];
	anglelist = g_slist_append( anglelist,app);
	draw_the_molecule( AMMP_draw_area,NULL);
	return 1==1;
} 
gint add_pair_to_torsion_list( GtkWidget *wind, gpointer *dummy)
{
	ATOM_QUARTET *app;
	if( picked_queue[0] == NULL || picked_queue[1] == NULL
	|| picked_queue[2] == NULL || picked_queue[3] == NULL) return 1==1;

	app = malloc( sizeof(ATOM_QUARTET) );
	app->a1 = picked_queue[0];
	app->a2 = picked_queue[1];
	app->a3 = picked_queue[2];
	app->a4 = picked_queue[3];
	torsionlist = g_slist_append( torsionlist,app);
	draw_the_molecule( AMMP_draw_area,NULL);
	return 1==1;
} 
gint clear_bond_list( GtkWidget *wind, gpointer *dummy)
{
	GSList *node;
	ATOM_PAIR *app;
	/* walk the list and free each object */
	node = bondlist;
	if( node == NULL) return 1==1;
	while( node != NULL)
	{
		app = (ATOM_PAIR *) node->data;
		free( app);
		node = node->next;
	}
	/* then free the list */
	g_slist_free( bondlist);
	bondlist = NULL;
	return 1==1;
} 
gint clear_angle_list( GtkWidget *wind, gpointer *dummy)
{
	GSList *node;
	ATOM_TRIPLE *app;
	/* walk the list and free each object */
	node = anglelist;
	if( node == NULL) return 1==1;
	while( node != NULL)
	{
		app = (ATOM_TRIPLE *) node->data;
		free( app);
		node = node->next;
	}
	/* then free the list */
	g_slist_free( anglelist);
	anglelist = NULL;
	return 1==1;
} 
gint clear_torsion_list( GtkWidget *wind, gpointer *dummy)
{
	GSList *node;
	ATOM_QUARTET *app;
	/* walk the list and free each object */
	node = torsionlist;
	if( node == NULL) return 1==1;
	while( node != NULL)
	{
		app = (ATOM_QUARTET *) node->data;
		free( app);
		node = node->next;
	}
	/* then free the list */
	g_slist_free( torsionlist);
	torsionlist = NULL;
	return 1==1;
} 
void draw_bond_list(  GtkWidget *where)
{
	GSList *node;
	float r,dx,dy,dz;
	char line[25];
	ATOM_PAIR *app;
	ATOM *a1,*a2;
	node = bondlist;
	if( node == NULL) return ;
	gdk_gc_set_foreground( Gcontext,&redpen);
	while( node != NULL)
	{
		
		app = (ATOM_PAIR *)node->data;
		a1 =  app->a1;
		a2 =  app->a2;
		if( a1->is_visible && a2->is_visible){
		dx = a1->x -a2->x;
		dy = a1->y -a2->y;
		dz = a1->z -a2->z;
		r = sqrt( dx*dx + dy*dy + dz*dz);
		sprintf(line,"%f",r);
	gdk_draw_line( pix,Gcontext, a1->ix,a1->iy, a2->ix,a2->iy);
	gdk_draw_string( pix, AMMPfont, Gcontext, (a1->ix+a2->ix)/2,(a1->iy+a2->iy)/2,line);
				}
		node = node->next;
	}
}
void draw_angle_list(  GtkWidget *where)
{
	GSList *node;
	float r,dx,dy,dz;
	float cth,ddx,ddy,ddz;
	char line[25];
	ATOM_TRIPLE *app;
	ATOM *a1,*a2,*a3;
	node = anglelist;
	if( node == NULL) return ;
	gdk_gc_set_foreground( Gcontext,&greenpen);
	while( node != NULL)
	{
		
		app = (ATOM_TRIPLE *)node->data;
		a1 =  app->a1;
		a2 =  app->a2;
		a3 =  app->a3;
		if( a1->is_visible && a2->is_visible && a3->is_visible){
		dx = a1->x -a2->x;
		dy = a1->y -a2->y;
		dz = a1->z -a2->z;
		ddx = a3->x -a2->x;
		ddy = a3->y -a2->y;
		ddz = a3->z -a2->z;
		r = sqrt( (dx*dx + dy*dy + dz*dz)*(ddx*ddx+ddy*ddy+ddz*ddz));
		if( r  >1.e-7)
		{
		cth = dx*ddx+dy*ddy + dz*ddz;
		cth /= r;
		cth   = acos( cth) * 180./PI;
		} else { cth = 0.;}
		sprintf(line,"%f",cth);
	gdk_draw_line( pix,Gcontext, a1->ix,a1->iy, a3->ix,a3->iy);
	gdk_draw_string( pix, AMMPfont, Gcontext, (a1->ix+a3->ix)/2,(a1->iy+a3->iy)/2,line);
				}
		node = node->next;
	}
}
void draw_torsion_list(  GtkWidget *where)
{
	GSList *node;
	float r;
	char line[25];
	ATOM_QUARTET *app;
	ATOM *a1,*a2,*a3,*a4;
	float get_torsion_value(); /* in tset.c */
	node = torsionlist;
	if( node == NULL) return ;
	gdk_gc_set_foreground( Gcontext,&bluepen);
	while( node != NULL)
	{
		
		app = (ATOM_QUARTET *)node->data;
		a1 =  app->a1;
		a2 =  app->a2;
		a3 =  app->a3;
		a4 =  app->a4;
		if( a1->is_visible && a2->is_visible && 
			a3->is_visible && a4->is_visible){
		r = 180.*get_torsion_value( a1,a2,a3,a4)/PI;
		sprintf(line,"%f",r);
	gdk_draw_line( pix,Gcontext, a2->ix,a2->iy, a3->ix,a3->iy);
	gdk_draw_string( pix, AMMPfont, Gcontext, (a2->ix+a3->ix)/2,(a2->iy+a3->iy)/2,line);
				}
		node = node->next;
	}
}
void draw_noel()
{
	int noel_next(); /* steps through the noel list */
	ATOM *a1,*a2;
	int i,ismore,ismet;

	i = 0;
	while(1==1)
	{
	ismore = noel_next( i++, &a1,&a2,&ismet);
	if( a1 == NULL || a2 == NULL) break;
	if( a1->is_visible && a2->is_visible){
	if( ismet )
	{
	if( !noel_violation_flag){
	gdk_gc_set_foreground( Gcontext,&greenpen);
	gdk_draw_line( pix,Gcontext, a1->ix,a1->iy, a2->ix,a2->iy);
		}
	}else{	
	gdk_gc_set_foreground( Gcontext,&redpen);
	gdk_draw_line( pix,Gcontext, a1->ix,a1->iy, a2->ix,a2->iy);
	}
	}
	if( !ismore ) break;	
	}
	return;
}
void draw_tether(GtkWidget *window)
{
	int tether_next(); /* steps through the noel list */
	ATOM *a1;
	float x,y,z;
	float tx,ty, width,height;
	int ix,iy;
	int i,ismore,ismet;

	ismore = tether_next( 0, &a1,&x,&y,&z);
	if( !ismore ) return ;
	gdk_gc_set_foreground( Gcontext,&purplepen);
	width = (float)(window->allocation.width)*0.5;
	height = (float)(window->allocation.height)*0.5;
	i = 0;
	while(1==1)
	{
	ismore = tether_next( i++, &a1,&x,&y,&z);
	if( a1 == NULL ) break;
	if( a1->is_visible ){
	x -= AMMPcx;
	y -= AMMPcy;
	z -= AMMPcz;

	tx = rotate[0][0]*x + rotate[0][1]*y + rotate[0][2]*z;
	ty = rotate[1][0]*x + rotate[1][1]*y + rotate[1][2]*z;

	ix = scale*tx + width;
	iy = scale*ty + height;
	gdk_draw_line( pix,Gcontext, a1->ix,a1->iy, ix,iy);

	}
	if( !ismore ) break;	
	}
	return;
}

gint torsion_search_handler( GtkWidget *root, gpointer *data)
{
	GSList *node;
	ATOM_QUARTET *app;
	ATOM *a1,*a2,*a3,*a4;
	char *cp,line[256];
	gint AMMP_indirect_handler();

	node = torsionlist;
	if( node == NULL ) return 1==1;

	while( node != NULL)
	{
		app = (ATOM_QUARTET *)node->data;
		a1 =  app->a1;
		a2 =  app->a2;
		a3 =  app->a3;
		a4 =  app->a4;
		node = node->next;
	}
/* at this stage a1,...,a4 point to the last torsion defined */
	if( a1 == NULL ) return 1==1;
	if( a2 == NULL ) return 1==1;
	if( a3 == NULL ) return 1==1;
	if( a4 == NULL ) return 1==1;
	cp = (char *)data;

	if( strcmp( cp,"tmin" ) == 0)
	{
	sprintf(line,"tmin %d %d %d %d 12 ;",a1->serial,a2->serial,a3->serial,a4->serial);
	}else{
	sprintf(line,"tgroup 1 %d %d %d %d 0 12 ; tsearch 1 ;",a1->serial,a2->serial,a3->serial,a4->serial);
	}
	AMMP_indirect_handler( NULL, line);
	return 1==1;
}
gint use_set_torsion_dialog( GtkWidget *root, gpointer *data)
{
	GSList *node;
	ATOM_QUARTET *app;
	ATOM *a1,*a2,*a3,*a4;
	BUNCH_OF_WIDGETS *wp;
	char *cp;
	float get_torsion_value();
	float target,original,delta;


	node = torsionlist;
	if( node == NULL ) return 1==1;

	wp = (BUNCH_OF_WIDGETS *)data;	
	while( node != NULL)
	{
		app = (ATOM_QUARTET *)node->data;
		a1 =  app->a1;
		a2 =  app->a2;
		a3 =  app->a3;
		a4 =  app->a4;
		node = node->next;
	}
/* at this stage a1,...,a4 point to the last torsion defined */
	if( a1 == NULL ) return 1==1;
	if( a2 == NULL ) return 1==1;
	if( a3 == NULL ) return 1==1;
	if( a4 == NULL ) return 1==1;
/* now get the values from the  dialog box */
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[2]));
	sscanf( cp,"%f", &target);
	
	original = get_torsion_value(a1,a2,a3,a4);
	tset_bond_build(a1,a2,a3,a4);
	delta = target*3.14159265/180.-original ;
	set_torsion(a1,a2,a3,a4,delta);
	AMMP_whiz_bang();	
	
	return 1==1;
}
gint use_define_noel_dialog( GtkWidget *caller, gpointer *data)
{
	GSList *node;
	ATOM_PAIR *app;
	ATOM *a1,*a2;
	BUNCH_OF_WIDGETS *wp;
	float target,dm,du,km,ku;
	char *cp;
	int noel();


	node = bondlist;
	if( node == NULL ) return 1==1;
	wp = (BUNCH_OF_WIDGETS *)data;	

	while( node != NULL)
	{
		app = (ATOM_PAIR *)node->data;
		a1 =  app->a1;
		a2 =  app->a2;
		node = node->next;
	}
	if( a1 == NULL ) return 1==1;
	if( a2 == NULL ) return 1==1;
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[2]));
	sscanf( cp,"%f", &target);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[4]));
	sscanf( cp,"%f", &dm);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[6]));
	sscanf( cp,"%f", &du);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[8]));
	sscanf( cp,"%f", &km);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[10]));
	sscanf( cp,"%f", &ku);
	noel(a1->serial, a2->serial, target,dm,du,km,ku);
	
	return 1==1;
}
gint use_define_restrain_dialog( GtkWidget *caller, gpointer *data)
{
	GSList *node;
	ATOM_PAIR *app;
	ATOM *a1,*a2;
	BUNCH_OF_WIDGETS *wp;
	float target,k;
	char *cp;
	int restrain();


	node = bondlist;
	if( node == NULL ) return 1==1;
	wp = (BUNCH_OF_WIDGETS *)data;	

	while( node != NULL)
	{
		app = (ATOM_PAIR *)node->data;
		a1 =  app->a1;
		a2 =  app->a2;
		node = node->next;
	}
	if( a1 == NULL ) return 1==1;
	if( a2 == NULL ) return 1==1;
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[2]));
	sscanf( cp,"%f", &target);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[4]));
	sscanf( cp,"%f", &k);
	restrain(a1->serial, a2->serial, target,k);
	return 1==1;
}
