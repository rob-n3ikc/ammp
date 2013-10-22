/* these routines build the menus
*  and hook them into dummy routines
*  which will be changed into real routines as needed
*/
#include <stdio.h>
#include <stdlib.h>
#include <gtk/gtk.h>
#include <signal.h>

/* external file which is where to write commands for 
*  buffered AMMP level access */
extern FILE *AMMP_QUEUED_COMMAND_FILE;
int AMMP_QUEUE_flag = (1==0);
int AMMP_QUEUE_ready = (1==0);
int AMMP_in_fileselect = (1==0);

/* external routines in eval */
void EVAL_assert_block( void );
void EVAL_assert_free( void );

/* structure used to pass a flag and a widget to a routine */
typedef struct
{ GtkWidget *which;
  void *flag;}  FANCY_WIDGET_DATA;

typedef struct
{
	int inme;
	GtkWidget *bunch[10];
} BUNCH_OF_WIDGETS;

/* flag definitions used as 0,1,2,3,4   */
static char f0,f1,f2,f3,f4;

/* block_unblock( flag)
* block if flag
* free if not
*
* manipulate the sigusr1 signal
*/
void block_unblock(int flag)
{

	struct sigaction wakeup;
	void no_operation(); /* in ammp.c */
/*block signals so fgets works */
	sigemptyset(&wakeup.sa_mask);
	sigaddset(&wakeup.sa_mask, SIGUSR1);

	sigprocmask(SIG_SETMASK, &wakeup.sa_mask,NULL);
	if( flag)
	sigprocmask(SIG_BLOCK, &wakeup.sa_mask,NULL);
	else
	sigprocmask(SIG_UNBLOCK, &wakeup.sa_mask,NULL);

	wakeup.sa_handler =  no_operation;
	wakeup.sa_flags =   SA_NOMASK;
	sigaction( SIGUSR1, &wakeup,NULL);
	return;
}

gint dummy_handler( GtkWidget *what, gpointer *data)
{
/* living dangerously aren't we? */
if( data != NULL) g_print("%s\n",data);
return 0;
}
extern int AMMP_EVAL_BLOCK;
gint AMMP_blocked_direct_handler( GtkWidget *what, gpointer *data)
{
int eval();
FILE *op, *get_current_op();
if( data == NULL ) return 0;
if( AMMP_EVAL_BLOCK){
g_print("the command\n %s\n blocked for program stablity\n",data); return 0;}
EVAL_assert_block();
	block_unblock(TRUE);
op = get_current_op();
if( op == NULL ) op = stdout;
eval( stdin, op, (char *)data);
	block_unblock(FALSE);
EVAL_assert_free();
return 1==1;
}

gint AMMP_indirect_handler( GtkWidget *what, gpointer *data )
{
	char *cp;
	
	cp = (char *) data;
	if( AMMP_QUEUE_flag) return 1==1;
AMMP_QUEUE_flag = (1==1);
	fprintf(AMMP_QUEUED_COMMAND_FILE,"%s;",cp);
	fflush(AMMP_QUEUED_COMMAND_FILE);
AMMP_QUEUE_flag = (1==0);
AMMP_QUEUE_ready = (1==1);
	return 1==1;
}
/* parse strings ; strings; ... into 
* individual ammp calls 
*/
gint AMMP_script_handler( GtkWidget *what, gpointer *data)
{
	int i;
	char line[256];
	char *cp;
	
	cp = (char *) data;
AGAIN: ;
	i = 0;
	while( *cp != ';' && *cp != '\0' && i < 256)
	{
		line[i++] = *(cp++);
	}
	line[i] = '\0';
	if( i > 0) AMMP_indirect_handler( what,(gpointer *)&line[0]);
	if( *cp !=  '\0') cp ++;
	if( *cp != '\0') goto AGAIN;
	return TRUE;
}

gint AMMP_direct_handler( GtkWidget *what, gpointer *data)
{
int eval();
FILE *op, *get_current_op();
if( data == NULL ) return 0;
op = get_current_op();
if( op == NULL ) op = stdout;
return eval( stdin, op, (char *)data);
}

/* File_read_write is linked to the OK button on a
* fileselectwidget
*  the data is what is passed at the gtk_signal_connect
*  and that's what i'll use to decide what to do
*  (i.e. what to dump or to read)
*
*  the data structure MUST contain a pointer to the file selector
*/
gint File_read_write( GtkWidget *what, gpointer *data)
{
	FANCY_WIDGET_DATA *ap;
	GtkWidget *wp;
	char *cp,*datap,line[256],touse[256];
	FILE *ip,*op,*fopen(),*get_current_op();
	void read_eval_do();
	int eval();

	AMMP_in_fileselect = (1==1);
	ap = (FANCY_WIDGET_DATA *)data;
	wp = (ap->which);
	datap = ap->flag;

	cp = gtk_file_selection_get_filename( GTK_FILE_SELECTION(wp) );

	if( *datap == (char)0)
	{
/*
	ip = stdin;
	op  = get_current_op();
	sprintf(line,"read %s",cp);
	eval(ip,op,"echo off");
	eval(ip,op,line);
	eval(ip,op,"echo ");
*/
/*
    the following should work but requires flagging on the draw side
*/
	sprintf(line,"echo off; read %s ; echo",cp);
	AMMP_indirect_handler( what, (gpointer *)&line[0]);

	}else{

/* block */
	block_unblock(TRUE);
	op  = fopen( cp,"w");
	if( op == NULL) {g_print("cannot open %s\n",cp); return (1==0);}
	ip = stdin;
	if( *datap == (char)1 || *datap == (char)3 )
	{
	if( *datap == (char)1)dump_atoms(op);
	dump_bonds(op);
	dump_angles(op);
	dump_hybrids(op);
	dump_torsions(op);
	dump_tethers(op);
	dump_noels(op);
	dump_orbit(op);
	dump_steps(op);
	dump_restrains(op);
	dump_table(op);
/*	dump_centrals(op); */
	dump_morse(op);
	dump_tbond(op);
	dump_tgroup(op);
	dump_variable(op);
	dump_swarms(op);
	dump_abcs(op);	
	dump_av5s(op);
	dump_ttargets(op);
/*	eval(ip,op,"dump atom bond angle hybrid torsion tether  ");
	eval(ip,op,"dump noel orbit step restrain table central ");
	eval(ip,op,"dump  morse tbond tgroup variable swarm ");
*/
	}
	if( *datap == (char)2)
	dump_atoms(op);
	if( *datap == (char)4)
	dump_pdb(op,100);

	fclose(op);
/* unblock */
	block_unblock(FALSE);
	}/* *datap selection */	
	AMMP_in_fileselect = (0==1);

/*	gtk_grab_remove(wp);
*/
	gtk_widget_destroy(wp);
return (1==1);
} /* File_read_write() */
void deselect_and_destroy( GtkWidget *dp, gpointer *dd)
{
	gtk_grab_remove(dp);
	gtk_widget_destroy( dp);
}

void destroy_widget(GtkWidget *wp, GtkWidget *dd)
{
	gtk_widget_destroy( dd);
}

/* File_selector is called by the menu and constructs
*  the file selector widget
* and destroyes it too 
*/
gint File_selector(GtkWidget *caller, gpointer *data)
{
static FANCY_WIDGET_DATA mine;
static GtkWidget *select;
char *lp;

	lp = (char *)data;

	if( *lp == (char)0)
	select = gtk_file_selection_new("Choose an Input File");
	else
	select = gtk_file_selection_new("Choose an Output File");

	mine.which = select;
	mine.flag = data;
	block_unblock(TRUE);

	gtk_signal_connect( GTK_OBJECT(GTK_FILE_SELECTION(select)->ok_button),
	"clicked", (GtkSignalFunc)File_read_write,&mine);
	gtk_signal_connect_object( GTK_OBJECT(GTK_FILE_SELECTION(select)->cancel_button),
	"clicked", (GtkSignalFunc)gtk_widget_destroy,GTK_OBJECT(select));


	gtk_signal_connect( GTK_OBJECT(select), "destroy", 
		(GtkSignalFunc)destroy_widget, select);

	gtk_widget_show(select);
	gtk_grab_add(select);
	block_unblock(FALSE);
}/* endof File_selector */


void build_file_menu( GtkWidget *main_menu)
{
static GtkWidget *read_ammp,*write_ammp,*write_atom,*write_geom,*write_pdb;
static GtkWidget *suicide;
GtkWidget *hline;

static GtkWidget *sub_menu;
sub_menu = gtk_menu_new();

read_ammp = gtk_menu_item_new_with_label("Read AMMP file");
write_ammp = gtk_menu_item_new_with_label("Write AMMP file");
write_atom = gtk_menu_item_new_with_label("Write AMMP atom file");
write_geom = gtk_menu_item_new_with_label("Write geometry file");
write_pdb = gtk_menu_item_new_with_label("Write pdb file");
suicide = gtk_menu_item_new_with_label("Exit");
hline = gtk_menu_item_new();


gtk_menu_append( GTK_MENU(sub_menu), read_ammp);
gtk_menu_append( GTK_MENU(sub_menu), write_ammp);
gtk_menu_append( GTK_MENU(sub_menu), write_atom);
gtk_menu_append( GTK_MENU(sub_menu), write_geom);
gtk_menu_append( GTK_MENU(sub_menu), write_pdb);
gtk_menu_append( GTK_MENU(sub_menu), hline);
gtk_menu_append( GTK_MENU(sub_menu), suicide);

f0 = 0;
f1 = 1;
f2 = 2;
f3 = 3;
f4 = 4;
gtk_signal_connect( GTK_OBJECT(read_ammp), "activate",
			  GTK_SIGNAL_FUNC(File_selector),&f0);
gtk_signal_connect( GTK_OBJECT(write_ammp), "activate",
			  GTK_SIGNAL_FUNC(File_selector),&f1);
gtk_signal_connect( GTK_OBJECT(write_atom), "activate",
			  GTK_SIGNAL_FUNC(File_selector),&f2);
gtk_signal_connect( GTK_OBJECT(write_geom), "activate",
			  GTK_SIGNAL_FUNC(File_selector),&f3);
gtk_signal_connect( GTK_OBJECT(write_pdb), "activate",
			  GTK_SIGNAL_FUNC(File_selector),&f4);
gtk_signal_connect( GTK_OBJECT(suicide), "activate",
			  GTK_SIGNAL_FUNC(AMMP_direct_handler),"exit");

gtk_widget_show( read_ammp);
gtk_widget_show( write_ammp);
gtk_widget_show( write_atom);
gtk_widget_show( write_geom);
gtk_widget_show( write_pdb);
gtk_widget_show( suicide);
gtk_widget_show(hline);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(main_menu), sub_menu);

}

void build_show_menu( GtkWidget *main_menu)
{
static GtkWidget *distance,*angle,*torsion,*noel,*tether,*clone;
static GtkWidget *add_distance,*clear_distance;
static GtkWidget *add_torsion,*clear_torsion;
static GtkWidget *add_angle,*clear_angle;
static GtkWidget *show_noel,*show_noel_violation;
static GtkWidget *show_tether;
static GtkWidget *clear_labels;

static GtkWidget *sub_menu;
static GtkWidget *sub_distance;
static GtkWidget *sub_angle;
static GtkWidget *sub_torsion;
static GtkWidget *sub_noel;
static GtkWidget *sub_tether;
gint add_pair_to_bond_list(), clear_bond_list();
gint add_pair_to_angle_list(), clear_angle_list();
gint add_pair_to_torsion_list(), clear_torsion_list();
gint toggle_noel_flag(), toggle_noel_violation_flag();
gint toggle_tether_flag();
gint clear_picked_list();

sub_menu = gtk_menu_new();
sub_distance = gtk_menu_new();
sub_angle = gtk_menu_new();
sub_torsion = gtk_menu_new();
sub_noel = gtk_menu_new();
sub_tether = gtk_menu_new();
f0 = 0;
f1 = 1;
f2 = 2;
f3 = 3;
f4 = 4;

show_noel  = gtk_check_menu_item_new_with_label("show Noels");
gtk_widget_show(show_noel);
gtk_check_menu_item_set_state( GTK_CHECK_MENU_ITEM( show_noel),FALSE);
show_noel_violation  = gtk_check_menu_item_new_with_label("show violations");
gtk_widget_show(show_noel_violation);
gtk_check_menu_item_set_state( GTK_CHECK_MENU_ITEM( show_noel_violation),FALSE);

show_tether  = gtk_check_menu_item_new_with_label("show tethers");
gtk_widget_show(show_tether);
gtk_check_menu_item_set_state( GTK_CHECK_MENU_ITEM( show_tether),FALSE);

distance = gtk_menu_item_new_with_label("Distance");
angle = gtk_menu_item_new_with_label("Angle");
torsion = gtk_menu_item_new_with_label("torsion");
noel = gtk_menu_item_new_with_label("Noel");
tether = gtk_menu_item_new_with_label("Tether");
clone = gtk_menu_item_new_with_label("Clones");

add_distance = gtk_menu_item_new_with_label("Show it");
clear_distance = gtk_menu_item_new_with_label("clear");
add_angle = gtk_menu_item_new_with_label("Show it");
clear_angle = gtk_menu_item_new_with_label("clear");
add_torsion = gtk_menu_item_new_with_label("Show it");
clear_torsion = gtk_menu_item_new_with_label("clear");
clear_labels = gtk_menu_item_new_with_label("clear labels");


gtk_menu_append( GTK_MENU(sub_menu), distance);
gtk_menu_append( GTK_MENU(sub_menu), angle);
gtk_menu_append( GTK_MENU(sub_menu), torsion);
gtk_menu_append( GTK_MENU(sub_menu), noel);
gtk_menu_append( GTK_MENU(sub_menu), tether);
gtk_menu_append( GTK_MENU(sub_menu), clone);
gtk_menu_append( GTK_MENU(sub_menu), clear_labels);

gtk_menu_append( GTK_MENU(sub_distance),add_distance);
gtk_menu_append( GTK_MENU(sub_distance),clear_distance);
gtk_menu_append( GTK_MENU(sub_angle),add_angle);
gtk_menu_append( GTK_MENU(sub_angle),clear_angle);
gtk_menu_append( GTK_MENU(sub_torsion),add_torsion);
gtk_menu_append( GTK_MENU(sub_torsion),clear_torsion);

gtk_menu_append( GTK_MENU(sub_noel),show_noel);
gtk_menu_append( GTK_MENU(sub_noel),show_noel_violation);
gtk_menu_append( GTK_MENU(sub_tether),show_tether);
/*
gtk_signal_connect( GTK_OBJECT(distance), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"Distance");
gtk_signal_connect( GTK_OBJECT(angle), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"angle");
gtk_signal_connect( GTK_OBJECT(noel), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"noel");
*/
gtk_signal_connect( GTK_OBJECT(show_noel),"activate",
			GTK_SIGNAL_FUNC( toggle_noel_flag),NULL);
gtk_signal_connect( GTK_OBJECT(show_noel_violation),"activate",
			GTK_SIGNAL_FUNC( toggle_noel_violation_flag),NULL);
gtk_signal_connect( GTK_OBJECT(show_tether), "activate",
			  GTK_SIGNAL_FUNC(toggle_tether_flag),"NULL");
gtk_signal_connect( GTK_OBJECT(clone), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"clone");
gtk_signal_connect( GTK_OBJECT(add_distance), "activate",
			  GTK_SIGNAL_FUNC(add_pair_to_bond_list),NULL);
gtk_signal_connect( GTK_OBJECT(clear_distance), "activate",
			  GTK_SIGNAL_FUNC(clear_bond_list),NULL);
gtk_signal_connect( GTK_OBJECT(add_angle), "activate",
			  GTK_SIGNAL_FUNC(add_pair_to_angle_list),NULL);
gtk_signal_connect( GTK_OBJECT(clear_angle), "activate",
			  GTK_SIGNAL_FUNC(clear_angle_list),NULL);
gtk_signal_connect( GTK_OBJECT(add_torsion), "activate",
			  GTK_SIGNAL_FUNC(add_pair_to_torsion_list),NULL);
gtk_signal_connect( GTK_OBJECT(clear_torsion), "activate",
			  GTK_SIGNAL_FUNC(clear_torsion_list),NULL);
gtk_signal_connect( GTK_OBJECT(clear_labels), "activate",
			  GTK_SIGNAL_FUNC(clear_picked_list),NULL);

gtk_widget_show( distance );
gtk_widget_show( angle );
gtk_widget_show( torsion );
gtk_widget_show( noel );
gtk_widget_show( tether );
gtk_widget_show( clone );
gtk_widget_show( add_distance );
gtk_widget_show( clear_distance );
gtk_widget_show( add_angle );
gtk_widget_show( clear_angle );
gtk_widget_show( clear_labels );
gtk_widget_show( add_torsion );
gtk_widget_show( clear_torsion );

gtk_menu_item_set_submenu( GTK_MENU_ITEM(main_menu), sub_menu);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(distance), sub_distance);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(angle), sub_angle);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(torsion), sub_torsion);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(noel), sub_noel);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(tether), sub_tether);

}/* build_show_menu */

void build_view_menu( GtkWidget *main_menu)
{
static GtkWidget *orient,*center,*scale,*color;
static GtkWidget *sub_menu,*sub_center,*sub_orient,*sub_scale;
static GtkWidget *center_on_mean,*center_on_atom,*track_center;
static GtkWidget *rotxy,*rotyz,*scaleslab,*track_scale;
static GtkWidget *sub_down,*down,*down_x,*down_y,*down_z;
static GtkWidget *sub_color,*color_mono, *color_atom,*color_force;
static  GSList *g_orient=NULL;
static  GSList *g_center=NULL;
static  GSList *g_color=NULL;
gint set_orient_flag_xy();
gint set_orient_flag_yz();
gint set_orient_flag_scale();
gint set_center_flag_mean();
gint set_center_flag_atom();
gint set_center_flag_track();
gint set_scale_flag();
gint set_orientation();
gint set_draw_color_type();
f0 = 0;
f1 = 1;
f2 = 2;
f3 = 3;
f4 = 4;

sub_menu = gtk_menu_new();
sub_center = gtk_menu_new();
sub_orient = gtk_menu_new();
sub_scale = gtk_menu_new();
sub_down = gtk_menu_new();
sub_color = gtk_menu_new();


orient = gtk_menu_item_new_with_label("Orientation Control");
center = gtk_menu_item_new_with_label("Center");
scale = gtk_menu_item_new_with_label("Scale");
color = gtk_menu_item_new_with_label("Color");
down = gtk_menu_item_new_with_label("Look Down");
down_x = gtk_menu_item_new_with_label("Down X axis");
down_y = gtk_menu_item_new_with_label("Down Y axis");
down_z = gtk_menu_item_new_with_label("Down Z axis");

rotxy = gtk_radio_menu_item_new_with_label(g_orient,"X-Y");
g_orient =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(rotxy));
rotyz = gtk_radio_menu_item_new_with_label(g_orient,"Y-Z");
g_orient =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(rotyz));
scaleslab= gtk_radio_menu_item_new_with_label(g_orient,"Scale/Slab");
g_orient =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(scaleslab));

gtk_menu_append( GTK_MENU( sub_orient), rotxy);
gtk_menu_append( GTK_MENU( sub_orient), rotyz);
gtk_menu_append( GTK_MENU( sub_orient), scaleslab);

gtk_menu_append( GTK_MENU(sub_down),down_x);
gtk_menu_append( GTK_MENU(sub_down),down_y);
gtk_menu_append( GTK_MENU(sub_down),down_z);

track_center = gtk_radio_menu_item_new_with_label(g_center,"track the center");
g_center =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(track_center));
center_on_mean = gtk_radio_menu_item_new_with_label(g_center,"Center on Mean");
g_center =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(center_on_mean));
gtk_widget_show(center_on_mean);
center_on_atom = gtk_radio_menu_item_new_with_label(g_center,"Center on last atom");
g_center =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(center_on_atom));
gtk_widget_show(center_on_atom);
gtk_widget_show(track_center);
gtk_menu_append( GTK_MENU( sub_center), track_center);
gtk_menu_append( GTK_MENU( sub_center), center_on_atom);
gtk_menu_append( GTK_MENU( sub_center), center_on_mean);

track_scale   = gtk_check_menu_item_new_with_label("track the scale");
gtk_widget_show(track_scale);
gtk_check_menu_item_set_state( GTK_CHECK_MENU_ITEM( track_scale),TRUE);
gtk_menu_append( GTK_MENU( sub_scale), track_scale);


color_atom = gtk_radio_menu_item_new_with_label(g_color,"CPK");
g_color =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(color_atom));
color_mono = gtk_radio_menu_item_new_with_label(g_color,"Mono");
g_color =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(color_mono));
color_force = gtk_radio_menu_item_new_with_label(g_color,"Force");
g_color =  gtk_radio_menu_item_group( GTK_RADIO_MENU_ITEM(color_force));
gtk_widget_show( color_mono);
gtk_widget_show( color_atom);
 gtk_widget_show( color_force); /* not yet */
gtk_menu_append( GTK_MENU(sub_color), color_mono);
gtk_menu_append( GTK_MENU(sub_color), color_atom);
gtk_menu_append( GTK_MENU(sub_color), color_force);

gtk_menu_append( GTK_MENU(sub_menu), orient);
gtk_menu_append( GTK_MENU(sub_menu), center);
gtk_menu_append( GTK_MENU(sub_menu), scale);
gtk_menu_append( GTK_MENU(sub_menu), color);
gtk_menu_append( GTK_MENU(sub_menu), down);

/*
gtk_signal_connect( GTK_OBJECT(scale), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"scale");
gtk_signal_connect( GTK_OBJECT(color), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"color");
*/
gtk_signal_connect( GTK_OBJECT(color_mono), "activate",
			  GTK_SIGNAL_FUNC(set_draw_color_type),&f0);
gtk_signal_connect( GTK_OBJECT(color_atom), "activate",
			  GTK_SIGNAL_FUNC(set_draw_color_type),&f1);
gtk_signal_connect( GTK_OBJECT(color_force), "activate",
			  GTK_SIGNAL_FUNC(set_draw_color_type),&f2);
gtk_signal_connect( GTK_OBJECT(rotxy),"activate",
			GTK_SIGNAL_FUNC(set_orient_flag_xy),NULL);
gtk_signal_connect( GTK_OBJECT(rotyz),"activate",
			GTK_SIGNAL_FUNC(set_orient_flag_yz),NULL);
gtk_signal_connect( GTK_OBJECT(scaleslab),"activate",
			GTK_SIGNAL_FUNC(set_orient_flag_scale),NULL);
gtk_signal_connect( GTK_OBJECT(center_on_mean),"activate",
			GTK_SIGNAL_FUNC(set_center_flag_mean),NULL);
gtk_signal_connect( GTK_OBJECT(center_on_atom),"activate",
			GTK_SIGNAL_FUNC(set_center_flag_atom),NULL);
gtk_signal_connect( GTK_OBJECT(track_center),"activate",
			GTK_SIGNAL_FUNC(set_center_flag_track),NULL);
gtk_signal_connect( GTK_OBJECT(track_scale),"activate",
			GTK_SIGNAL_FUNC(set_scale_flag),NULL);
gtk_signal_connect( GTK_OBJECT(down_x),"activate",
			GTK_SIGNAL_FUNC(set_orientation),"X");
gtk_signal_connect( GTK_OBJECT(down_y),"activate",
			GTK_SIGNAL_FUNC(set_orientation),"Y");
gtk_signal_connect( GTK_OBJECT(down_z),"activate",
			GTK_SIGNAL_FUNC(set_orientation),"Z");

gtk_widget_show( rotxy);
gtk_widget_show( rotyz);
gtk_widget_show( scaleslab);
gtk_widget_show( orient);
gtk_widget_show( center);
gtk_widget_show( scale);
gtk_widget_show( color);
gtk_widget_show( down);
gtk_widget_show( down_x);
gtk_widget_show( down_y);
gtk_widget_show( down_z);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(main_menu), sub_menu);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(orient), sub_orient);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(center), sub_center);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(scale), sub_scale);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(down), sub_down);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(color), sub_color);
}/* build_view */

void build_analyze_menu( GtkWidget *main_menu)
{
static GtkWidget *sub_menu;
static GtkWidget *bond,*angle,*nonbon,*hybrid,*noel,*tether,*torsion;
/* the functions */
int a_noel();
int a_bond(),a_angle();
int a_nonbon(),a_torsion(),a_hybrid();
int a_tether();
gint AMMP_analyze_handler();

sub_menu = gtk_menu_new();

bond = gtk_menu_item_new_with_label("Bond");
angle = gtk_menu_item_new_with_label("Angle");
nonbon = gtk_menu_item_new_with_label("Non-Bonded");
hybrid = gtk_menu_item_new_with_label("hybrid");
noel = gtk_menu_item_new_with_label("Noel");
tether = gtk_menu_item_new_with_label("tether");
torsion = gtk_menu_item_new_with_label("torsion");

gtk_menu_append( GTK_MENU(sub_menu),bond);
gtk_menu_append( GTK_MENU(sub_menu),angle);
gtk_menu_append( GTK_MENU(sub_menu),hybrid);
gtk_menu_append( GTK_MENU(sub_menu),torsion);
gtk_menu_append( GTK_MENU(sub_menu),nonbon);
gtk_menu_append( GTK_MENU(sub_menu),noel);
gtk_menu_append( GTK_MENU(sub_menu),tether);

gtk_signal_connect( GTK_OBJECT(bond),"activate",
			GTK_SIGNAL_FUNC(AMMP_analyze_handler),a_bond);
gtk_signal_connect( GTK_OBJECT(angle),"activate",
			GTK_SIGNAL_FUNC(AMMP_analyze_handler),a_angle);
gtk_signal_connect( GTK_OBJECT(nonbon),"activate",
			GTK_SIGNAL_FUNC(AMMP_analyze_handler),a_nonbon);
gtk_signal_connect( GTK_OBJECT(tether),"activate",
			GTK_SIGNAL_FUNC(AMMP_analyze_handler),a_tether);
gtk_signal_connect( GTK_OBJECT(noel),"activate",
			GTK_SIGNAL_FUNC(AMMP_analyze_handler),a_noel);
gtk_signal_connect( GTK_OBJECT(torsion),"activate",
			GTK_SIGNAL_FUNC(AMMP_analyze_handler),a_torsion);
gtk_signal_connect( GTK_OBJECT(hybrid),"activate",
			GTK_SIGNAL_FUNC(AMMP_analyze_handler),a_hybrid);

gtk_widget_show( bond);
gtk_widget_show( angle);
gtk_widget_show( nonbon);
gtk_widget_show( hybrid);
gtk_widget_show( noel);
gtk_widget_show( tether);
gtk_widget_show( torsion);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(main_menu), sub_menu);

}

void build_task_menu( GtkWidget *main_menu)
{
static GtkWidget *potentials,*dynamics,*minimize,*build,*set_charge;
static GtkWidget *fixem,*clone,*monitor;
static GtkWidget *freeze_all,*freeze_nz,*tether_all,*thaw_all;
static GtkWidget *build_new,*build_part;
static GtkWidget *sub_menu_search,*torsion_ops, *tgroup, *tmin,*tset;
static GtkWidget *hline, *torsion_group_label;

static GtkWidget *sub_menu,*sub_menu_freeze,*sub_menu_build;
gint build_dynamic_dialog();
gint build_minimize_dialog();
gint build_set_charge_dialog();
gint build_set_torsion_dialog();
gint torsion_search_handler();

f0 = 0;
f1 = 1;
f2 = 2;
f3 = 3;
f4 = 4;

sub_menu = gtk_menu_new();
sub_menu_freeze = gtk_menu_new();
sub_menu_build = gtk_menu_new();
sub_menu_search = gtk_menu_new();

torsion_group_label = gtk_menu_item_new_with_label("Pick and display the torsion first");
hline = gtk_menu_item_new();

potentials = gtk_menu_item_new_with_label("Choose Potentials");
dynamics = gtk_menu_item_new_with_label("Run Molecular dynamics");
minimize = gtk_menu_item_new_with_label("Minimize the energy");
build = gtk_menu_item_new_with_label("Build");
set_charge = gtk_menu_item_new_with_label("set charges");
fixem = gtk_menu_item_new_with_label("Freeze/Thaw atoms");
clone = gtk_menu_item_new_with_label("Clone operations");
torsion_ops = gtk_menu_item_new_with_label("Manipulate Torsions");
monitor = gtk_menu_item_new_with_label("Monitor");

freeze_all = gtk_menu_item_new_with_label("freeze All atoms");
freeze_nz = gtk_menu_item_new_with_label("freeze built atoms");
tether_all = gtk_menu_item_new_with_label("tether built atoms");
thaw_all = gtk_menu_item_new_with_label("Thaw All atoms");

build_new  = gtk_menu_item_new_with_label("Build from scratch");
build_part  = gtk_menu_item_new_with_label("Build from Partial structure");

tgroup = gtk_menu_item_new_with_label("side chain <200 atoms move (fast)");
tmin = gtk_menu_item_new_with_label("Any torsion (slow but general)");
tset = gtk_menu_item_new_with_label("Set the torsion");

gtk_menu_append( GTK_MENU(sub_menu), potentials);
gtk_menu_append( GTK_MENU(sub_menu), build);
gtk_menu_append( GTK_MENU(sub_menu), torsion_ops);
gtk_menu_append( GTK_MENU(sub_menu), minimize);
gtk_menu_append( GTK_MENU(sub_menu), dynamics);
gtk_menu_append( GTK_MENU(sub_menu), set_charge);
gtk_menu_append( GTK_MENU(sub_menu), fixem);
gtk_menu_append( GTK_MENU(sub_menu), clone);
gtk_menu_append( GTK_MENU(sub_menu), monitor);

gtk_menu_append( GTK_MENU(sub_menu_freeze),freeze_all);
gtk_menu_append( GTK_MENU(sub_menu_freeze),freeze_nz);
gtk_menu_append( GTK_MENU(sub_menu_freeze),tether_all);
gtk_menu_append( GTK_MENU(sub_menu_freeze),thaw_all);

gtk_menu_append( GTK_MENU(sub_menu_build), build_new);
gtk_menu_append( GTK_MENU(sub_menu_build), build_part);

gtk_menu_append( GTK_MENU(sub_menu_search), torsion_group_label);
gtk_menu_append( GTK_MENU(sub_menu_search), hline);
gtk_menu_append( GTK_MENU(sub_menu_search), tset);
gtk_menu_append( GTK_MENU(sub_menu_search), tgroup);
gtk_menu_append( GTK_MENU(sub_menu_search), tmin);

gtk_signal_connect( GTK_OBJECT(potentials) ,"activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"potentials");
gtk_signal_connect( GTK_OBJECT(minimize), "activate",
			  GTK_SIGNAL_FUNC(build_minimize_dialog),"minimize");
gtk_signal_connect( GTK_OBJECT(dynamics), "activate",
			  GTK_SIGNAL_FUNC(build_dynamic_dialog),"dynamics");

gtk_signal_connect( GTK_OBJECT(set_charge), "activate",
			  GTK_SIGNAL_FUNC(build_set_charge_dialog),"set charge");
gtk_signal_connect( GTK_OBJECT(clone), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"clone");
gtk_signal_connect( GTK_OBJECT(monitor), "activate",
			  GTK_SIGNAL_FUNC(AMMP_indirect_handler),"monitor");
gtk_signal_connect( GTK_OBJECT(freeze_all), "activate",
	  GTK_SIGNAL_FUNC(AMMP_blocked_direct_handler),"inactive  0 100000000");
gtk_signal_connect( GTK_OBJECT(freeze_nz), "activate",
	  GTK_SIGNAL_FUNC(AMMP_blocked_direct_handler),"nzinactive  0 100000000");
gtk_signal_connect( GTK_OBJECT(tether_all), "activate",
	  GTK_SIGNAL_FUNC(AMMP_blocked_direct_handler),"tether all 100");
gtk_signal_connect( GTK_OBJECT(thaw_all), "activate",
	  GTK_SIGNAL_FUNC(AMMP_blocked_direct_handler),"active  0 10000000");
gtk_signal_connect( GTK_OBJECT(build_new),"activate",
	GTK_SIGNAL_FUNC(AMMP_script_handler),"gsdg 20");
gtk_signal_connect( GTK_OBJECT(build_part),"activate",
	GTK_SIGNAL_FUNC(AMMP_script_handler),"nzinactive 0 10000000; abuild  10");

gtk_signal_connect( GTK_OBJECT(tset),"activate",
	GTK_SIGNAL_FUNC( build_set_torsion_dialog),NULL);
gtk_signal_connect( GTK_OBJECT(tgroup),"activate",
		GTK_SIGNAL_FUNC(torsion_search_handler),"tsearch");
gtk_signal_connect( GTK_OBJECT(tmin),"activate",
		GTK_SIGNAL_FUNC(torsion_search_handler),"tmin");

gtk_widget_show( potentials);
gtk_widget_show( build);
gtk_widget_show( minimize);
gtk_widget_show( dynamics);
gtk_widget_show( clone);
gtk_widget_show( set_charge);
gtk_widget_show( fixem);
gtk_widget_show( monitor);
gtk_widget_show(freeze_all);
gtk_widget_show(freeze_nz);
gtk_widget_show(tether_all);
gtk_widget_show(thaw_all);
gtk_widget_show( build_new);
gtk_widget_show( build_part);
gtk_widget_show( torsion_ops);
gtk_widget_show( hline);
gtk_widget_show( tset);
gtk_widget_show( tmin);
gtk_widget_show( tgroup);
gtk_widget_show( torsion_group_label);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(main_menu), sub_menu);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(fixem), sub_menu_freeze);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(build), sub_menu_build);
gtk_menu_item_set_submenu( GTK_MENU_ITEM(torsion_ops), sub_menu_search);

}
void build_define_menu( GtkWidget *main_menu)
{
	static GtkWidget *instruction,*hline,*define_noel,*define_restrain;
	static GtkWidget *sub_menu;
	gint build_define_noel();
	gint build_define_restrain();

	sub_menu = gtk_menu_new();
instruction  = gtk_menu_item_new_with_label("Pick and Display the two atoms first");
	hline = gtk_menu_item_new();
define_noel = gtk_menu_item_new_with_label("Noel term ");
define_restrain = gtk_menu_item_new_with_label("restrain term ");
	gtk_menu_append( GTK_MENU(sub_menu), instruction);
	gtk_menu_append( GTK_MENU(sub_menu), hline);
	gtk_menu_append( GTK_MENU(sub_menu), define_noel);
	gtk_menu_append( GTK_MENU(sub_menu), define_restrain);
	gtk_widget_show( instruction);
	gtk_widget_show( hline);
	gtk_widget_show( define_noel);
	gtk_widget_show( define_restrain);
	gtk_menu_item_set_submenu( GTK_MENU_ITEM(main_menu), sub_menu);
gtk_signal_connect( GTK_OBJECT(define_noel),"activate",
	GTK_SIGNAL_FUNC( build_define_noel),NULL);
gtk_signal_connect( GTK_OBJECT(define_restrain),"activate",
	GTK_SIGNAL_FUNC( build_define_restrain),NULL);
}/* end of build_define_menu */

void build_help_menu( GtkWidget *main_menu)
{

gtk_signal_connect( GTK_OBJECT(main_menu), "activate",
			  GTK_SIGNAL_FUNC(dummy_handler),"Help");
}

