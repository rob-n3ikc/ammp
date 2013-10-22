/* these routines build the dialogboxes
*
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <gtk/gtk.h>

/* defined in 1.2 and above but not in 1.0 */
gboolean gtk_toggle_button_get_active( GtkToggleButton *what)
{
	if(  what->active) return TRUE;
	return FALSE;
}

typedef struct
{ int inme;
 GtkWidget *bunch[12];
}  BUNCH_OF_WIDGETS ;

gint use_set_charge_dialog( GtkWidget *root, gpointer *data)
{
	BUNCH_OF_WIDGETS *wp;
	char line[256],*cp;
	int i,j,k;
	float tq;


	wp = (BUNCH_OF_WIDGETS *) data;
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[2]) );
	sscanf(cp,"%d", &i);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[4]) );
	sscanf(cp,"%d", &j);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[6]) );
	sscanf(cp,"%f", &tq);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[9]) );
	sscanf(cp,"%d", &k);
	 sprintf(line,"momadd %d %d ; mom %f %d",i,j,tq,k);
	 AMMP_indirect_handler(root,(gpointer *)&line[0]);

	return (1==1);
}
gint build_set_torsion_dialog( GtkWidget *root, gpointer *data)
{
	static GtkWidget *the_box,*label;
	static GtkWidget *ok, *cancel;
	static BUNCH_OF_WIDGETS mine;
	int i;
	gint use_set_torsion_dialog();  /* this will be in draw
                                        * because it needs access to the
                                        * clicked list */
	the_box = gtk_dialog_new();
	label = gtk_label_new("Set a Torsion Angle ");
	gtk_box_pack_start( GTK_BOX( GTK_DIALOG(the_box)->vbox ),label,
		TRUE,TRUE,0);
	gtk_widget_show(label);
	ok = gtk_button_new_with_label("Run");
	cancel = gtk_button_new_with_label("Close");
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),ok,
	TRUE,TRUE,0);
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),cancel,
	TRUE,TRUE,0);
/* now build the rest of the box */
	mine.bunch[0] = gtk_hseparator_new();
	mine.bunch[1] = gtk_label_new("Value (in degrees)");
	mine.bunch[2] = gtk_entry_new();
	mine.inme = 3;
	for( i=0; i< mine.inme; i++)
	{
		gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->vbox),
			mine.bunch[i],TRUE,TRUE,0);
		gtk_widget_show( mine.bunch[i]);
	}
	gtk_widget_show(ok);
	gtk_widget_show(cancel);
/* connect to signals here */
	gtk_signal_connect_object( GTK_OBJECT(cancel),
	"clicked", (GtkSignalFunc)gtk_widget_destroy,GTK_OBJECT(the_box));
	gtk_signal_connect( GTK_OBJECT(ok),
	"clicked", (GtkSignalFunc)use_set_torsion_dialog,(gpointer *)&mine);

	gtk_widget_show(the_box);	
	
	return 1==1;
}
gint build_set_charge_dialog( GtkWidget *root , gpointer *data)
{
	static GtkWidget *the_box,*label;
	static GtkWidget *ok,*cancel;
	static BUNCH_OF_WIDGETS mine;
	int i;

	the_box = gtk_dialog_new();
	label = gtk_label_new("Set up Charge Generation");
	gtk_box_pack_start( GTK_BOX( GTK_DIALOG(the_box)->vbox ),label,
		TRUE,TRUE,0);
	gtk_widget_show(label);
	ok = gtk_button_new_with_label("Run");
	cancel = gtk_button_new_with_label("Close");
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),ok,
	TRUE,TRUE,0);
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),cancel,
	TRUE,TRUE,0);
/* now set up the text/input parts of the box */
	mine.bunch[0] = gtk_hseparator_new();
	mine.bunch[1] = gtk_label_new("First atom");
	mine.bunch[2] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[2]),"0");
	mine.bunch[3] = gtk_label_new("Last atom");
	mine.bunch[4] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[4]),"100");
	mine.bunch[5] = gtk_label_new("Total Charge");
	mine.bunch[6] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[6]),"0.");
	mine.bunch[7] = gtk_hseparator_new();
	mine.bunch[8] = gtk_label_new("Number of iterations");
	mine.bunch[9] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[9]),"100");
	mine.inme = 10;
	for( i=0; i< mine.inme; i++)
	{
		gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->vbox),
			mine.bunch[i],TRUE,TRUE,0);
		gtk_widget_show( mine.bunch[i]);
	}
	gtk_widget_show(ok);
	gtk_widget_show(cancel);
/* connect to signals here */
	gtk_signal_connect_object( GTK_OBJECT(cancel),
	"clicked", (GtkSignalFunc)gtk_widget_destroy,GTK_OBJECT(the_box));
	gtk_signal_connect( GTK_OBJECT(ok),
	"clicked", (GtkSignalFunc)use_set_charge_dialog,(gpointer *)&mine);

	gtk_widget_show(the_box);	
	
	return 1==1;
} /* build_minimize_dialog */

gint use_minimize_dialog( GtkWidget *root, gpointer *data)
{
	BUNCH_OF_WIDGETS *wp;
	char line[256],*cp;
	int numatom,a_number(); /* in the main ammp section */
	gboolean min_cng,min_bfgs,min_steep,min_trunc;

	numatom = a_number();
	if( numatom < 1) return TRUE;
	if( numatom > 100) AMMP_indirect_handler(root,(gpointer *)"setf mxdq 0.75");
	if( numatom > 1000) AMMP_indirect_handler(root,(gpointer *)"setf mmbox 10.");

	wp = (BUNCH_OF_WIDGETS *) data;
	min_cng = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(wp->bunch[6]));
	min_bfgs = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(wp->bunch[7]));
	min_steep = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(wp->bunch[8]));
	min_trunc = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(wp->bunch[9]));
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[2]) );
	sscanf(cp,"%d", &numatom);
	if( min_cng) sprintf(line,"cngdel %d",numatom);
	if( min_bfgs) sprintf(line,"bfgs %d",numatom);
	if( min_steep) sprintf(line,"steep %d",numatom);
	if( min_trunc) sprintf(line,"trunc %d 0.01 0.01",numatom);
	 AMMP_indirect_handler(root,(gpointer *)&line[0]);

	return (1==1);
}
gint build_minimize_dialog( GtkWidget *root , gpointer *data)
{
	static GtkWidget *the_box,*label;
	static GtkWidget *ok,*cancel;
	static BUNCH_OF_WIDGETS mine;
	static GSList *group ;
	int i;

	group = NULL; /* since this dialog may be called multiple times it must initialize the group */
	the_box = gtk_dialog_new();
	label = gtk_label_new("Set up Energy Minimization");
	gtk_box_pack_start( GTK_BOX( GTK_DIALOG(the_box)->vbox ),label,
		TRUE,TRUE,0);
	gtk_widget_show(label);
	ok = gtk_button_new_with_label("Run");
	cancel = gtk_button_new_with_label("Close");
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),ok,
	TRUE,TRUE,0);
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),cancel,
	TRUE,TRUE,0);
/* now set up the text/input parts of the box */
	mine.bunch[0] = gtk_hseparator_new();
	mine.bunch[1] = gtk_label_new("Number of steps");
	mine.bunch[2] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[2]),"100");
	mine.bunch[3] = gtk_hseparator_new();
	mine.bunch[4] = gtk_label_new("Algorithm");
	mine.bunch[5] = gtk_hseparator_new();
	mine.bunch[6] = gtk_radio_button_new_with_label(group,"conjugate gradients");
	group = gtk_radio_button_group( GTK_RADIO_BUTTON(mine.bunch[6]) );
	mine.bunch[7] = gtk_radio_button_new_with_label(group,"Quasi-Newton");
	group = gtk_radio_button_group( GTK_RADIO_BUTTON(mine.bunch[7]) );
	mine.bunch[8] = gtk_radio_button_new_with_label(group,"Steepest Descents");
	group = gtk_radio_button_group( GTK_RADIO_BUTTON(mine.bunch[8]) );
	mine.bunch[9] = gtk_radio_button_new_with_label(group,"Truncated Newton");
	group = gtk_radio_button_group( GTK_RADIO_BUTTON(mine.bunch[9]) );
	mine.inme = 10;
	for( i=0; i< mine.inme; i++)
	{
		gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->vbox),
			mine.bunch[i],TRUE,TRUE,0);
		gtk_widget_show( mine.bunch[i]);
	}
	gtk_widget_show(ok);
	gtk_widget_show(cancel);
/* connect to signals here */
	gtk_signal_connect_object( GTK_OBJECT(cancel),
	"clicked", (GtkSignalFunc)gtk_widget_destroy,GTK_OBJECT(the_box));
	gtk_signal_connect( GTK_OBJECT(ok),
	"clicked", (GtkSignalFunc)use_minimize_dialog,(gpointer *)&mine);

	gtk_widget_show(the_box);	
	
	return 1==1;
} /* build_minimize_dialog */


gint use_dynamic_dialog( GtkWidget *root, gpointer *data)
{
	BUNCH_OF_WIDGETS *wp;
	char line[256],*cp;
	int numatom,a_number(); /* in the main ammp section */
	gboolean dyn_init,dyn_temp;
	float dt,T;

	numatom = a_number();
	if( numatom < 1) return TRUE;
	if( numatom > 100) AMMP_indirect_handler(root,(gpointer *)"setf mxdq 0.75");
	if( numatom > 1000) AMMP_indirect_handler(root,(gpointer *)"setf mmbox 10.");

	wp = (BUNCH_OF_WIDGETS *) data;
	
	dyn_init = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(wp->bunch[7]));
	dyn_temp = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON(wp->bunch[8]));

	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[2]) );
	sscanf(cp,"%d", &numatom);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[4]) );
	sscanf(cp,"%f",&dt);
	cp = gtk_entry_get_text( GTK_ENTRY(wp->bunch[6]) );
	sscanf(cp,"%f",&T);
	if( dyn_init){
	sprintf(line,"v_maxwell %f",T);
	AMMP_indirect_handler( root, (gpointer*)line);
	}
	if( dt > .00002 || dt < 0.0000001) dt = 0.00001;
	if( dyn_temp)
	{
	sprintf(line,"tpac %d %f %f",numatom,dt,T);
	}else{
	sprintf(line,"pac %d %f",numatom,dt);
	}
	AMMP_indirect_handler( root, (gpointer*)line);
	return (1==1);
}
gint build_dynamic_dialog( GtkWidget *root , gpointer *data)
{
	static GtkWidget *the_box,*label;
	static GtkWidget *ok,*cancel;
	static BUNCH_OF_WIDGETS mine;
	int i;

	the_box = gtk_dialog_new();
	label = gtk_label_new("Set up Molecular Dynamics");
	gtk_box_pack_start( GTK_BOX( GTK_DIALOG(the_box)->vbox ),label,
		TRUE,TRUE,0);
	gtk_widget_show(label);
	ok = gtk_button_new_with_label("Run");
	cancel = gtk_button_new_with_label("Close");
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),ok,
	TRUE,TRUE,0);
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),cancel,
	TRUE,TRUE,0);
/* now set up the text/input parts of the box */
	mine.bunch[0] = gtk_hseparator_new();
	mine.bunch[1] = gtk_label_new("Number of steps");
	mine.bunch[2] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[2]),"100");
	mine.bunch[3] = gtk_label_new("Time step");
	mine.bunch[4] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[4]),"0.00001");
	mine.bunch[5] = gtk_label_new("Temperature");
	mine.bunch[6] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[6]),"300");
	mine.bunch[7] = gtk_check_button_new_with_label("Initialize");
	mine.bunch[8] = gtk_check_button_new_with_label("Constrain Temperature");
	mine.inme = 9;
	for( i=0; i< mine.inme; i++)
	{
		gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->vbox),
			mine.bunch[i],TRUE,TRUE,0);
		gtk_widget_show( mine.bunch[i]);
	}
	gtk_widget_show(ok);
	gtk_widget_show(cancel);
/* connect to signals here */
	gtk_signal_connect_object( GTK_OBJECT(cancel),
	"clicked", (GtkSignalFunc)gtk_widget_destroy,GTK_OBJECT(the_box));
	gtk_signal_connect( GTK_OBJECT(ok),
	"clicked", (GtkSignalFunc)use_dynamic_dialog,(gpointer *)&mine);

	gtk_widget_show(the_box);	
	return 1==1;
} /* build_dynamic_dialog */

gint build_define_noel( GtkWidget *root, gpointer *data)
{
	static GtkWidget *the_box,*label;
	static GtkWidget *ok,*cancel;
	static BUNCH_OF_WIDGETS mine;
	gint use_define_noel_dialog();
	int i;

	the_box = gtk_dialog_new();
	label = gtk_label_new("Set up a NOEl term");
	gtk_box_pack_start( GTK_BOX( GTK_DIALOG(the_box)->vbox ),label,
		TRUE,TRUE,0);
	gtk_widget_show(label);
	ok = gtk_button_new_with_label("Define");
	cancel = gtk_button_new_with_label("Close");
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),ok,
	TRUE,TRUE,0);
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),cancel,
	TRUE,TRUE,0);
/* now set up the text/input parts of the box */
	mine.bunch[0] = gtk_hseparator_new();
	mine.bunch[1] = gtk_label_new("mean distance");
	mine.bunch[2] = gtk_entry_new();
	mine.bunch[3] = gtk_label_new("lower Deviation");
	mine.bunch[4] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[4]),"0");
	mine.bunch[5] = gtk_label_new("upper Deviation");
	mine.bunch[6] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[6]),"0");
	mine.bunch[7] = gtk_label_new("lower force constant");
	mine.bunch[8] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[8]),"100");
	mine.bunch[9] = gtk_label_new("upper force constant");
	mine.bunch[10] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[10]),"100");


	mine.inme = 11;
	for( i=0; i< mine.inme; i++)
	{
		gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->vbox),
			mine.bunch[i],TRUE,TRUE,0);
		gtk_widget_show( mine.bunch[i]);
	}

	gtk_widget_show(ok);
	gtk_widget_show(cancel);
/* connect to signals here */
	gtk_signal_connect_object( GTK_OBJECT(cancel),
	"clicked", (GtkSignalFunc)gtk_widget_destroy,GTK_OBJECT(the_box));
	gtk_signal_connect( GTK_OBJECT(ok),
	"clicked", (GtkSignalFunc)use_define_noel_dialog,(gpointer *)&mine);

	gtk_widget_show(the_box);	
	return 1==1;
}
gint build_define_restrain( GtkWidget *root, gpointer *data)
{
	static GtkWidget *the_box,*label;
	static GtkWidget *ok,*cancel;
	static BUNCH_OF_WIDGETS mine;
	gint use_define_restrain_dialog();
	int i;

	the_box = gtk_dialog_new();
	label = gtk_label_new("Set up a Restrain term");
	gtk_box_pack_start( GTK_BOX( GTK_DIALOG(the_box)->vbox ),label,
		TRUE,TRUE,0);
	gtk_widget_show(label);
	ok = gtk_button_new_with_label("Define");
	cancel = gtk_button_new_with_label("Close");
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),ok,
	TRUE,TRUE,0);
	gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->action_area),cancel,
	TRUE,TRUE,0);
/* now set up the text/input parts of the box */
	mine.bunch[0] = gtk_hseparator_new();
	mine.bunch[1] = gtk_label_new(" distance");
	mine.bunch[2] = gtk_entry_new();
	mine.bunch[3] = gtk_label_new("force constant");
	mine.bunch[4] = gtk_entry_new();
	gtk_entry_set_text(GTK_ENTRY(mine.bunch[4]),"100");


	mine.inme = 5;
	for( i=0; i< mine.inme; i++)
	{
		gtk_box_pack_start( GTK_BOX(GTK_DIALOG(the_box)->vbox),
			mine.bunch[i],TRUE,TRUE,0);
		gtk_widget_show( mine.bunch[i]);
	}

	gtk_widget_show(ok);
	gtk_widget_show(cancel);
/* connect to signals here */
	gtk_signal_connect_object( GTK_OBJECT(cancel),
	"clicked", (GtkSignalFunc)gtk_widget_destroy,GTK_OBJECT(the_box));
	gtk_signal_connect( GTK_OBJECT(ok),
	"clicked", (GtkSignalFunc)use_define_restrain_dialog,(gpointer *)&mine);

	gtk_widget_show(the_box);	
	return 1==1;
}
