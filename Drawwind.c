
/*

 routines
* starting with the AMMPdrawwndproc()
* and including the atomwalker and bondwalker routines
*/

#include <windows.h>
#include <shellapi.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <setjmp.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"
#define MAXPICK 40
#ifndef WIN32
#define WIN32
#endif
/* external definitions of important stuff */
extern HWND  AMMPmainhWnd;
extern HWND  AMMPdrawhWnd ;
extern int AMMPexecineval;
extern struct {int useit, atom_id; ATOM *ap;} choose_atom_returns;
extern struct {char atom[10]; float r; int useit;} select_draw_atom_returns;

extern int in_AMMP_mouse_proc;

BOOL CALLBACK choose_atom_dialog(HWND,UINT,UINT,LONG);
BOOL CALLBACK select_draw_atom_dialog(HWND,UINT,UINT,LONG);

long FAR PASCAL AMMPdrawwndproc( hWnd, wMessage,wParam,lParam)
HWND hWnd;
UINT wMessage,wParam;
LONG lParam;
{
HDC hdc,hdcm,hdcr;
static int pageH,pageV,slab=1000;
static float xr=0,yr=0,zr=0,scale=1,xc=0,yc=0,zc=0;
static float orient[3][3]={1.,0.,0.,0.,1.,0.,0.,0.,1.};
static void *picked[MAXPICK];
static char labels[MAXPICK][20],inlabels[MAXPICK];
static int inpicked = -1;
static void *distlist[MAXPICK][2], *anglist[MAXPICK][3];
static int indistlist = -1, inanglist = -1;
static int autoscale = 1; /* scale automatically */
static int autocenter = 1; /* center automatically */
static int hisy = 1, hisscale = 0; /* horizontal scroller meaning */
static int visx = 1, visz = 0, visslab = 0, vistor = 0 ; /* vertical scroller meaning */
static int nsV=0,nsH=0; /* scroller positions */
static int shownoel = 0; /* show noe's */
static int just_noel_violate = 0; /* if true only violations shown */
static int color_scheme = 1;
static int draw_tether = (0==1);
int i,ix,iy,icl,draw_itemp,numatm,a_number();
int it1,it2,it3,it4;
char line[80];
float x;
ATOM *ap,*apc,*a_next(),*first_AMMP_ATOM(); /* ammp ATOMs */
PAINTSTRUCT ps;
HPEN oldpen,redpen,greenpen;
HBITMAP newbitmap,oldbitmap;
BITMAP bdata;
RECT rect;
GLOBALHANDLE hGMem;
LPMETAFILEPICT lpMFP;
HMETAFILE hmf;
SIZE draw_size;
void draw_noel();
/*BOOL CALLBACK choose_atom_dialog();*/
void AMMP_mouse_proc();
void do_the_drawing();

switch( wMessage)
{
case WM_CREATE:
	if( inpicked < 0 )
	 {for( inpicked = 0; inpicked < MAXPICK; inpicked ++)
		picked[inpicked] = NULL;

					 inpicked = 0;
	 }
	 if( indistlist < 0 )
	 {for( indistlist = 0; indistlist < 20; indistlist++)
		distlist[indistlist][0] = NULL;
		indistlist = 0;
	 }
	 if( inanglist < 0)
	 { for( inanglist = 0; inanglist < 20; inanglist ++)
		anglist[inanglist][0] = NULL;
		inanglist = 0;
			}
	hdc = GetDC( hWnd );
	SetScrollRange( hWnd,SB_VERT,0,1000,FALSE);
		  SetScrollPos( hWnd,SB_VERT,0,TRUE);
	SetScrollRange( hWnd,SB_HORZ,0,1000,FALSE);
	SetScrollPos( hWnd,SB_HORZ,0,TRUE);
	CheckMenuItem(GetMenu(hWnd),201, MF_BYCOMMAND|MF_CHECKED);
			CheckMenuItem(GetMenu(hWnd),207, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),208, MF_BYCOMMAND|MF_CHECKED);
	EnableMenuItem(GetMenu(hWnd),209,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),301,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),302,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),306,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),307,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),310,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),311,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),312,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),313,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),314,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	EnableMenuItem(GetMenu(hWnd),315,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	DrawMenuBar(hWnd);
	ReleaseDC( hWnd,hdc);
	 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0);

	return 0;
case WM_SIZE:
	pageH = LOWORD( lParam);
	pageV = HIWORD( lParam);
	 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0);

	return 0;
case WM_LBUTTONDOWN:
case WM_MOUSEMOVE:
case WM_LBUTTONUP:
	AMMP_mouse_proc( hWnd,wMessage,wParam,lParam,nsH,nsV);
	return 0;
case WM_PAINT:
		 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0);

	hdc = BeginPaint(hWnd,&ps) ;
	hdcr = CreateCompatibleDC( hdc);
	if( hdcr != NULL ){ hdcm = hdcr;

	newbitmap = CreateCompatibleBitmap( hdc,pageH,pageV);
		  oldbitmap = SelectObject( hdcm,newbitmap);
	} else {hdcm = hdc;}

/* insert call to draw */
	do_the_drawing( 
		hdcm,pageV,pageH,slab,shownoel,color_scheme,just_noel_violate);
	if( draw_tether) 
		tetherwalkerpainter( orient,xr,yr,zr,&xc,&yc,&zc,&scale,
						(1==0),autocenter,pageH,pageV,0,hdcm,slab);

/* still do the picking and bond-angle stuff here */
	for(i=0; i< 20; i++)
	{ if( picked[i] == NULL) break;
	ap = picked[i];
	if( ap->iz >= -slab && ap->iz <= slab)
	  TextOut( hdcm,ap->ix,ap->iy,
	  &labels[i][0],inlabels[i]);
	  }

		  if( distlist[0][0] != NULL)
	{
	redpen = CreatePen( PS_DOT,1,RGB(255,0,0));
		  oldpen = SelectObject( hdcm,redpen);
		  SetTextColor( hdcm,RGB(255,0,0));
	distdraw( hdcm, distlist,slab );
	SetTextColor( hdcm,RGB(0,0,0));
	SelectObject( hdcm,oldpen);
	DeleteObject( redpen);
		  }
		  if( anglist[0][0] != NULL)
	{
	redpen = CreatePen( PS_DOT,1,RGB(124,0,255));
	oldpen = SelectObject( hdcm,redpen);
		  SetTextColor( hdcm,RGB(124,0,255));
	angledraw( hdcm,anglist,slab );
	SelectObject( hdcm,oldpen);
	DeleteObject( redpen);
		  }

	SetBkMode(hdcm,OPAQUE);
	if( hdcr != NULL)
	{/* bitblt it if memory draw worked */
	BitBlt( hdc,0,0,pageH,pageV,hdcm,0,0,SRCCOPY);
	SelectObject( hdcm,oldbitmap);
	DeleteObject( newbitmap );
	DeleteDC( hdcm);

		  }
	EndPaint( hWnd,&ps);
		 return 0;
case WM_LBUTTONDBLCLK:
	if( wParam == MK_LBUTTON){
	ix = LOWORD(lParam);iy = HIWORD(lParam);
	numatm = a_number();
		  apc = first_AMMP_ATOM(); icl = abs( apc->ix -ix)+abs(apc->iy-iy);
		  ap = first_AMMP_ATOM();
		  ap = ap->next;
	for( i=1;i<numatm; i++)
	{
					 if( ap->iz <= slab && ap->iz >=-slab){
		draw_itemp = abs(ap->ix-ix)+abs(ap->iy -iy);
		if( draw_itemp < icl)
					  {apc = ap; icl = draw_itemp;}
					  }
			ap = ap->next;
	}   /*
	printf(" the closest atom is %d %s\n",apc->serial,&apc->name[0]);
				 */
		 if( icl < 10 ){
		inlabels[inpicked] = sprintf( &labels[inpicked][0],
				 "%s %d",&apc->name[0],apc->serial);
							  picked[inpicked] = apc;
				 inpicked +=1; if(inpicked >= MAXPICK) inpicked = 0;
 if( inpicked == 1 )
		  {
	EnableMenuItem(GetMenu(hWnd),310,MF_BYCOMMAND|MF_ENABLED);
	EnableMenuItem(GetMenu(hWnd),311,MF_BYCOMMAND|MF_ENABLED);
	EnableMenuItem(GetMenu(hWnd),312,MF_BYCOMMAND|MF_ENABLED);
	EnableMenuItem(GetMenu(hWnd),313,MF_BYCOMMAND|MF_ENABLED);
	EnableMenuItem(GetMenu(hWnd),314,MF_BYCOMMAND|MF_ENABLED);
	EnableMenuItem(GetMenu(hWnd),315,MF_BYCOMMAND|MF_ENABLED);
	 }
 if( inpicked == 2 )
	EnableMenuItem(GetMenu(hWnd),301,MF_BYCOMMAND|MF_ENABLED);
 if( inpicked == 3 )
	EnableMenuItem(GetMenu(hWnd),302,MF_BYCOMMAND|MF_ENABLED);
if( inpicked == 4)
 {
	EnableMenuItem(GetMenu(hWnd),209,MF_BYCOMMAND|MF_ENABLED);
//	EnableMenuItem(GetMenu(hWnd),306,MF_BYCOMMAND|MF_ENABLED);
  }
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd)  ;

					}
						}
	return 0;
case WM_COMMAND:
	switch( wParam)
	{
	case 99: /* center on an atom */
		DialogBox( NULL,MAKEINTRESOURCE(6),hWnd,choose_atom_dialog);
		if( !choose_atom_returns.useit) break;
		ap = choose_atom_returns.ap;
		if( ap != NULL) {
		CheckMenuItem(GetMenu(hWnd),208, MF_BYCOMMAND|MF_UNCHECKED);
		autocenter = 0 ;
			xc = ap->x; yc = ap->y; zc = ap->z;
		atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0);  }
		InvalidateRect( hWnd,NULL,FALSE);
		UpdateWindow( hWnd);

	break;
	case 100: /* send to clipboard */
	hdcm = CreateMetaFile(NULL);
	if( hdcm == NULL ) break;
	draw_size.cx = pageH;
	draw_size.cy = pageV;
	SetWindowExtEx( hdcm,pageH, pageV,&draw_size);
/* code for drawing from WM_PAINT */

	/* insert call to draw */
	do_the_drawing( 
		hdcm,pageV,pageH,slab,shownoel,color_scheme,just_noel_violate);
	if( draw_tether) 
		tetherwalkerpainter( orient,xr,yr,zr,&xc,&yc,&zc,&scale,
						(1==0),autocenter,pageH,pageV,0,hdcm,slab);

/* still do the picking and bond-angle stuff here */


	SelectObject(hdcm,GetStockObject(SYSTEM_FIXED_FONT));
		  SetBkMode( hdcm,TRANSPARENT);
	for(i=0; i< 20; i++)
	{ if( picked[i] == NULL) break;
	ap = picked[i];
	if( ap->iz >= -slab && ap->iz <= slab)
	  TextOut( hdcm,ap->ix,ap->iy,
	  &labels[i][0],inlabels[i]);
	  }
		  if( distlist[0][0] != NULL)
	{
	redpen = CreatePen( PS_DOT,1,RGB(255,0,0));
		  oldpen = SelectObject( hdcm,redpen);
		  SetTextColor( hdcm,RGB(255,0,0));
	distdraw( hdcm, distlist,slab );
	SetTextColor( hdcm,RGB(0,0,0));
	SelectObject( hdcm,oldpen);
	DeleteObject( redpen);
		  }
		  if( anglist[0][0] != NULL)
	{
	redpen = CreatePen( PS_DOT,1,RGB(124,0,255));
	oldpen = SelectObject( hdcm,redpen);
		  SetTextColor( hdcm,RGB(124,0,255));
	angledraw( hdcm,anglist,slab );
	SelectObject( hdcm,oldpen);
	DeleteObject( redpen);
		  }

	hGMem = GlobalAlloc( GHND, (DWORD) sizeof(METAFILEPICT));
	lpMFP = (LPMETAFILEPICT) GlobalLock(hGMem);
	lpMFP->mm = MM_ISOTROPIC|MM_HIMETRIC;
	lpMFP->xExt = pageH*35;
	lpMFP->yExt = pageV*35;
	lpMFP->hMF = CloseMetaFile(hdcm);
	hmf = lpMFP->hMF;
	GlobalUnlock( hGMem);
	OpenClipboard( hWnd);
	EmptyClipboard();
	SetClipboardData(CF_METAFILEPICT, hGMem);
	CloseClipboard();


	break;
	case 101:
	orient[0][0] = 0;
	orient[0][1] = 0;
	orient[0][2] = -1;
	orient[1][0] = 0;
	orient[1][1] = 1;
	orient[1][2] = 0;
	orient[2][0] = 1;
	orient[2][1] = 0;
		  orient[2][2] = 0;
	xr = 0; yr = 3.14159265/2; zr = 0;
	/* call atomwalker()*/
 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
			 autoscale,autocenter,pageH,pageV,1);
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd)  ;
	break;
	case 102:
	orient[0][0] = 1;
	orient[0][1] = 0;
	orient[0][2] = 0;
	orient[1][0] = 0;
	orient[1][1] = 0;
	orient[1][2] = -1;
	orient[2][0] = 0;
	orient[2][1] = 1;
		  orient[2][2] = 0;
	xr = 3.14159265/2; yr = 0; zr = 0;
		/* call atomwalker()*/
 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,1)  ;
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd) ;

		  break;
	case 103:
	orient[0][0] = 1;
	orient[0][1] = 0;
	orient[0][2] = 0;
	orient[1][0] = 0;
	orient[1][1] = 1;
	orient[1][2] = 0;
	orient[2][0] = 0;
	orient[2][1] = 0;
		  orient[2][2] = 1;
	xr = 0; yr = 0; zr = 0;
		/* call atomwalker()*/
 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,1) ;
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd) ;

	break;
	case 104:
	color_scheme = 0;
	CheckMenuItem(GetMenu(hWnd),104, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),105, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),106, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),107, MF_BYCOMMAND|MF_UNCHECKED);
	InvalidateRect( hWnd,NULL,FALSE);
	UpdateWindow( hWnd) ;

	break;
	case 105:
		color_scheme = 1;
	CheckMenuItem(GetMenu(hWnd),105, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),104, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),106, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),107, MF_BYCOMMAND|MF_UNCHECKED);
	InvalidateRect( hWnd,NULL,FALSE);
	UpdateWindow( hWnd) ;

	break;
	case 106:
	color_scheme = 2;
	CheckMenuItem(GetMenu(hWnd),106, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),104, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),105, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),107, MF_BYCOMMAND|MF_UNCHECKED);
	InvalidateRect( hWnd,NULL,FALSE);
	UpdateWindow( hWnd) ;

	break;
	case 107:
	color_scheme = 3;
	CheckMenuItem(GetMenu(hWnd),107, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),104, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),106, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),105, MF_BYCOMMAND|MF_UNCHECKED);
	InvalidateRect( hWnd,NULL,FALSE);
	UpdateWindow( hWnd) ;

	break;
	
	case 108:
		DialogBox( NULL,MAKEINTRESOURCE(7),hWnd,select_draw_atom_dialog);
/*		if( !select_draw_atom_returns.useit) break;*/
/*		printf("%f %s\n",select_draw_atom_returns.r,select_draw_atom_returns.atom);*/
		InvalidateRect( hWnd,NULL,FALSE);
		UpdateWindow( hWnd);
	break;

	case 201: /* y-x*/
	hisy = 1; hisscale = 0;
	visx = 1; visz = 0; visslab = 0;  vistor = 0;
	nsH = 1000.*yr/6.2831853;
	nsV = 1000.*xr/6.2831853;
	CheckMenuItem(GetMenu(hWnd),201, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),202, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),203, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),209, MF_BYCOMMAND|MF_UNCHECKED);

		  SetScrollPos( hWnd,SB_HORZ,nsH,TRUE);
		  SetScrollPos( hWnd,SB_VERT,nsV,TRUE);
	break;
	case 202: /*y-z*/
	hisy = 1; hisscale = 0;
	visx = 0; visz = 1; visslab = 0; vistor = 0;
	nsH = 1000.*yr/6.2831853;
	nsV = 1000.*zr/6.2831853;
	CheckMenuItem(GetMenu(hWnd),201, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),202, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),203, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),209, MF_BYCOMMAND|MF_UNCHECKED);
		  SetScrollPos( hWnd,SB_HORZ,nsH,TRUE);
		  SetScrollPos( hWnd,SB_VERT,nsV,TRUE);
	break;
	case 209: /* y-Torsion */
	CheckMenuItem(GetMenu(hWnd),201, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),202, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),203, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),209, MF_BYCOMMAND|MF_CHECKED);
	visx = 0; visz = 0; visslab = 0;
	vistor = 1;
	hisy = 1; hisscale = 0;
	nsH = 1000.*yr/6.2831853;
	nsV = 1000.*zr/6.2831853;


	break;
	case 203: /*scale-slab*/
	hisy = 0; hisscale = 1;
	visx = 0; visz = 0; visslab = 1; vistor = 0;
	nsH = 10.*scale;
	nsV = slab;
	CheckMenuItem(GetMenu(hWnd),201, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),202, MF_BYCOMMAND|MF_UNCHECKED);
	CheckMenuItem(GetMenu(hWnd),203, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),209, MF_BYCOMMAND|MF_UNCHECKED);
		  SetScrollPos( hWnd,SB_HORZ,nsH,TRUE);
		  SetScrollPos( hWnd,SB_VERT,nsV,TRUE);
	break;
	case 204: /* center on last picked atom */
	 CheckMenuItem(GetMenu(hWnd),208, MF_BYCOMMAND|MF_UNCHECKED);
	 autocenter = 0;
	ix = inpicked;
	ix -= 1; if( ix < 0 ) ix = MAXPICK -1;
	if( picked[ix] != NULL){
	ap = picked[ix];
	xc = ap->x; yc = ap->y; zc = ap->z;
		atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd);
	  }
	break;
	case 205: /* center on mean position */
	if( !autocenter){
	autocenter = 1;
	atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd);
	autocenter = 0;
		  }
	break;
	case 206: /* clear picked list */
	 if( inpicked > 0 )
	 {for( inpicked = 0; inpicked < 20; inpicked ++)
		picked[inpicked] = NULL;

		inpicked = 0;
	 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd);
 if( vistor) {
			vistor = 0; visx = 1;
				nsH = 1000.*yr/6.2831853;
	nsV = 1000.*xr/6.2831853;
	CheckMenuItem(GetMenu(hWnd),201, MF_BYCOMMAND|MF_CHECKED);
	CheckMenuItem(GetMenu(hWnd),209, MF_BYCOMMAND|MF_UNCHECKED);
				}
 EnableMenuItem(GetMenu(hWnd),209,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),301,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),302,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),306,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),307,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),310,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),311,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),312,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),313,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),314,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
 EnableMenuItem(GetMenu(hWnd),315,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
	 }
	break;
	case 207:
	if( autoscale){
	CheckMenuItem(GetMenu(hWnd),207, MF_BYCOMMAND|MF_UNCHECKED);
	 autoscale = 0;}
	else {
	CheckMenuItem(GetMenu(hWnd),207, MF_BYCOMMAND|MF_CHECKED);
	autoscale = 1 ;  }
		  break;
	case 208:
	if( autocenter){
	 CheckMenuItem(GetMenu(hWnd),208, MF_BYCOMMAND|MF_UNCHECKED);
	 autocenter = 0;}
	else{
	 CheckMenuItem(GetMenu(hWnd),208, MF_BYCOMMAND|MF_CHECKED);
	 autocenter = 1 ;}
		  break;
	case 301: /* add atoms to distance list */
	EnableMenuItem(GetMenu(hWnd),306,MF_BYCOMMAND|MF_ENABLED);
	i = inpicked -1;
	if( i < 0 ) i = MAXPICK-1;
	if( picked[i] == NULL ) break;
	distlist[indistlist][1] = picked[i];
	i -= 1;
	if( i < 0 ) i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	distlist[indistlist][0] = picked[i];
	indistlist += 1; if( indistlist == 20) indistlist = 0;
	atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
	 InvalidateRect( hWnd,NULL,FALSE);
	 UpdateWindow( hWnd);
	break;
	case 302: /* add atoms to angle list */
	EnableMenuItem(GetMenu(hWnd),307,MF_BYCOMMAND|MF_ENABLED);
	i = inpicked -1;
	if( i< 0 ) i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	anglist[inanglist][2] = picked[i];
	i -= 1;
	if( i < 0 ) i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	anglist[inanglist][1] = picked[i];
	i -= 1;
	if( i < 0 ) i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	anglist[inanglist][0] = picked[i];
	inanglist += 1; if( inanglist > 19) inanglist = 0;
	atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
	 InvalidateRect( hWnd,NULL,FALSE);
	 UpdateWindow( hWnd);

		  break;
	case 303: /* toggle noe display */
	if( shownoel)
	{
	CheckMenuItem(GetMenu(hWnd),303, MF_BYCOMMAND|MF_UNCHECKED);
	 shownoel = 0;
	 }else {
	CheckMenuItem(GetMenu(hWnd),303, MF_BYCOMMAND|MF_CHECKED);
	shownoel = 1 ;
	}
	 InvalidateRect( hWnd,NULL,FALSE);
	 UpdateWindow( hWnd);

	break;
	case 304: /* toggle noe violations only */
		if( just_noel_violate)
		{
		CheckMenuItem(GetMenu(hWnd),304,MF_BYCOMMAND|MF_UNCHECKED);
		just_noel_violate = 0;
		}else{
		CheckMenuItem(GetMenu(hWnd),304,MF_BYCOMMAND|MF_CHECKED);
		just_noel_violate = 1;
		}
	 InvalidateRect( hWnd,NULL,FALSE);
	 UpdateWindow( hWnd);

	break;
	case 305:
		if( draw_tether) {
			CheckMenuItem(GetMenu(hWnd),305, MF_BYCOMMAND|MF_UNCHECKED);
			draw_tether = (1==0); 
		}else{ 
			CheckMenuItem(GetMenu(hWnd),305,MF_BYCOMMAND|MF_CHECKED);
			draw_tether = (1==1);}
	 InvalidateRect( hWnd,NULL,FALSE);
	 UpdateWindow( hWnd);
	break;
	case 306: /* clear distance list */
	for( indistlist = 0; indistlist < 20; indistlist++)
		distlist[indistlist][0] = NULL;

		indistlist = 0;
EnableMenuItem(GetMenu(hWnd),306,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);

		atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
		 InvalidateRect( hWnd,NULL,FALSE);
		 UpdateWindow( hWnd);

		  break;
	case 307: /* clear angle list */
	 for( inanglist = 0; inanglist < 20; inanglist ++)
		anglist[inanglist][0] = NULL;

		inanglist = 0;
EnableMenuItem(GetMenu(hWnd),307,MF_BYCOMMAND|MF_GRAYED|MF_DISABLED);
		atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
		 InvalidateRect( hWnd,NULL,FALSE);
		 UpdateWindow( hWnd);
	break;
	case 321: /* setup torsion search */
	case 322: /* setup torsion search */
	case 323: /* set to 0 */
	case 325: /* set to 180 */
	if( AMMPexecineval)
	{MessageBox(NULL,"Please wait for the main process to finish",
			  "AMMP caution", MB_OK|MB_ICONEXCLAMATION);
			  break;
	}else{
	i = inpicked -1; if( i< 0 ) i = MAXPICK -1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it4 = ap->serial;
	i -= 1; if( i < 0 ) i = MAXPICK-1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it3 = ap->serial;
	i -= 1; if( i < 0 ) i = MAXPICK-1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it2 = ap->serial;
	i -= 1; if( i < 0 ) i = MAXPICK-1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it1 = ap->serial;
	if( it1 == it2 || it1 == it3 || it1 == it4 ||
	it2 == it3 || it2 == it4 || it3 == it4)
	{
	sprintf( line,"Please pick 4 unique atoms %d %d %d %d",it1,it2,it3
		  ,it4);
	MessageBox( NULL,line,"AMMP caution",MB_OK|MB_ICONEXCLAMATION);
		break;
	}
	if( wParam == 321 )
	{
	sprintf(line,"tgroup 1 %d %d %d %d 0. 24 ",it1,it2,it3,it4);
	eval( stdin,stdout,line);
	sprintf(line,"tsearch 1 ");
	eval( stdin,stdout,line);
	} else if ( wParam == 322 ){
	sprintf( line,"tmin %d %d %d %d 12 ",it1,it2,it3,it4);
	eval( stdin,stdout,line);
	} else if (wParam == 323 ){
	sprintf( line,"tset %d %d %d %d 0 ",it1,it2,it3,it4);
	eval( stdin,stdout,line);
	} else if (wParam == 325 ){
	sprintf( line,"tset %d %d %d %d 180 ",it1,it2,it3,it4);
	eval( stdin,stdout,line);
	}
	atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
	 InvalidateRect( hWnd,NULL,FALSE);
	 UpdateWindow( hWnd);

		  break;
	}

	case 310: /*analyze the bonds on the current picked atom */
	i = inpicked -1; if( i < 0 ) i = MAXPICK -1;
	if( picked[i] == NULL ) break;
		  ap = picked[i];
	a_bond( &x,0.,ap->serial,ap->serial,stdout);
	break;
	case 311: /* analyze the angles on the current picked atom */
	i = inpicked-1; if( i<0 )i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	ap = picked[i];
	a_angle( &x,0., ap->serial,ap->serial,stdout);
		  break;
	case 312: /* analyze the torsion on the current picked atom */
	i = inpicked-1; if( i<0 )i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	ap = picked[i];
	a_torsion( &x,0., ap->serial,ap->serial,stdout);
		  break;
	case 313: /* analyze the hybrid on the current picked atom */
	i = inpicked-1; if( i<0 )i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	ap = picked[i];
	a_hybrid( &x,0., ap->serial,ap->serial,stdout);
		  break;
	case 314: /* analyze the nonbonded on the current picked atom */
	i = inpicked-1; if( i<0 )i = MAXPICK -1;
	if( picked[i] == NULL ) break;
	ap = picked[i];
	a_nonbon( &x,0., ap->serial,ap->serial,stdout);
		  break;
	case 315: /* analyze noel terms */
	i = inpicked-1; if( i<0) i=MAXPICK -i;
	if( picked[i] == NULL ) break;
	ap = picked[i];
	a_noel( &x,0., ap->serial,ap->serial,stdout);
	default:
	break;
		  }

		  return 0;
case WM_SYSCOMMAND:
	if( wParam == SC_CLOSE)
	{
	ShowWindow( AMMPdrawhWnd, SW_HIDE);
		  return 0;
	}
	return DefWindowProc( hWnd,wMessage,wParam,lParam);

case WM_QUIT:
	 /*	PostQuitMessage(hWnd);  */
		  return 0;

case WM_KEYDOWN:
	switch( wParam)
	{
	case VK_HOME:
	if( nsV != 0)
		  SendMessage( hWnd,WM_VSCROLL,SB_TOP,0L);
	if( nsH != 0)
	SendMessage( hWnd,WM_HSCROLL,SB_TOP,0L);
	break;
	case VK_PRIOR:
		  SendMessage( hWnd,WM_VSCROLL,SB_PAGEDOWN,0L);
		  break;
	case VK_NEXT:
	SendMessage( hWnd,WM_VSCROLL,SB_PAGEUP,0L);
		  break;
	case VK_UP:
	SendMessage( hWnd,WM_VSCROLL,SB_LINEDOWN,0L);
	break;
	case VK_DOWN:
	SendMessage( hWnd,WM_VSCROLL,SB_LINEUP,0L);
	break;
	case VK_RIGHT:
	SendMessage( hWnd,WM_HSCROLL,SB_PAGEDOWN,0L);
	break;
	case VK_LEFT:
	SendMessage( hWnd,WM_HSCROLL,SB_PAGEUP,0L);
		  break;
		  }
	  return 0;


case WM_HSCROLL:
#ifdef WIN32
	switch( LOWORD(wParam))
#else
	switch(wparam)
#endif
	{
	case SB_LINEUP:
			  i=-1;
	break;
	case SB_LINEDOWN:
		i=1;
	break;
	case SB_PAGEUP:
			  i = -10;
	break;
	case SB_PAGEDOWN:
			  i = 10;
	break;
	case SB_THUMBPOSITION:
#ifdef WIN32
			  i = HIWORD(wParam) - nsH;
#else
				i = LOWORD(lParam) - nsH;
#endif
	break;
	case SB_TOP:
	i = -nsH;
		  break;
		  default:
		  i = 0;
	}
	nsH += i;
	if( i != 0)
	{
		if(hisy)
	{     	if( nsH < 0 ){  nsH += 1000;}
		if( nsH >= 1000 ){   nsH -= 1000;}
		 yr = (float) nsH/1000. *6.283153;

		 }

	else if (hisscale) {
		scale += (float)i/10; if(scale < 0 ){
				nsH = 0; scale = 0;}
					 }

	SetScrollPos( hWnd,SB_HORZ,(short int)nsH,TRUE);
	/* call atomwalker()*/
	if( !in_AMMP_mouse_proc) {
 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
						  autoscale,autocenter,pageH,pageV,0)   ;
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd);
									 }


	}

	return 0 ;

case WM_VSCROLL:
#ifdef WIN32
	switch( LOWORD(wParam))
#else
	switch( wParam)
#endif
	{
	case SB_LINEUP:
		  i = -1;
	break;
	case SB_LINEDOWN:
	i = 1;
	break;
	case SB_PAGEUP:
		  i =  -10;
	break;
	case SB_PAGEDOWN:
		  i = +10;
	break;
	case SB_THUMBPOSITION:
#ifdef WIN32
			  i = HIWORD(wParam) - nsV;
#else
				i = LOWORD(lParam) - nsV;
#endif

		  break;
	case SB_TOP:
	i = 1;
		  nsV = -1;
		  break;
	default:
		  i = 0;

	}
	nsV += i;
	if( i != 0)
	{

	if( visx)
	{
		if( nsV < 0 ){  nsV += 1000;}
		if( nsV > 1000 ) {nsV -= 1000;}
		 xr = (float) nsV/1000. *6.283153;
		 }
	else if( visz)
		{
		if( nsV < 0 ){ nsV += 1000;}
		if( nsV > 1000 ) { nsV -= 1000;}
		zr = (float) nsV/1000. *6.283153;
		  }
	else if( vistor)
		{

	if( AMMPexecineval)
	{MessageBox(NULL,"Please wait for the main process to finish",
			  "AMMP caution", MB_OK|MB_ICONEXCLAMATION);
			  break;
	}else{
	i = inpicked -1; if( i< 0 ) i = MAXPICK -1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it4 = ap->serial;
	i -= 1; if( i < 0 ) i = MAXPICK-1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it3 = ap->serial;
	i -= 1; if( i < 0 ) i = MAXPICK-1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it2 = ap->serial;
	i -= 1; if( i < 0 ) i = MAXPICK-1;
	if( picked[i] == NULL) break;
	ap = picked[i];  it1 = ap->serial;
	if( it1 == it2 || it1 == it3 || it1 == it4 ||
	it2 == it3 || it2 == it4 || it3 == it4)
	{
	sprintf( line,"Please pick 4 unique atoms %d %d %d %d",it1,it2,it3
		  ,it4);
	MessageBox( NULL,line,"AMMP caution",MB_OK|MB_ICONEXCLAMATION);
		break;
	}
/* just call tset directly */

	tset( stdout,0,it1,it2,it3,it4, (float)nsV/1000.*6.283153);
	}
		 }
	else if (visslab) {
		slab += i; if( slab < 0 ) slab = 0;
		}
		SetScrollPos( hWnd,SB_VERT,(short int)nsV,TRUE);
		/* call atomwalker()*/
 atomwalker(orient, xr,yr,zr,&xc,&yc,&zc,&scale,
			 autoscale,autocenter,pageH,pageV,0) ;
 InvalidateRect( hWnd,NULL,FALSE);
 UpdateWindow( hWnd) ;


		}
	return 0;





}
return DefWindowProc( hWnd,wMessage,wParam,lParam);


}


int atomwalker(orient, xr,yr,zr,xc,yc,zc,scale,
                    autoscale,autocenter,pageH,pageV,force)
float orient[3][3];
float xr,yr,zr,*xc,*yc,*zc,*scale;
int autoscale, autocenter,force;
int pageH,pageV;
{
ATOM *ap,*a_next(),*first_AMMP_ATOM();
int numatm,a_number(),i;
static float xrk=0,yrk=0,zrk=0;
float r1[3][3];/* incremental rotation matrix */
float x,y,z,xx,yy,zz;
int phc,pvc;  /* center of page in device units*/
int xdmi,xdma,ydmi,ydma; /* for scaling */
	numatm = a_number();
	if( numatm <= 0) return;
		  phc = pageH/2; pvc = pageV/2;
	if( autocenter )
	{ /*this is done here because we need the center*/
	*xc = 0; *yc = 0; *zc = 0;
	ap = first_AMMP_ATOM();
	for( i=0; i< numatm; i++)
	{
	  *xc += ap->x; *yc += ap->y; *zc += ap->z;
	 ap = ap->next;
	  }
	  *xc /= numatm; *yc /= numatm; *zc /= numatm;
	  }
			 if( autoscale){
			 *scale = 100; /* too big, but will scale down */
	 }
	 if( force){ xrk = xr; yrk = yr; zrk = zr;}
/* if xr has changed update the matrixes */
	if( xr != xrk ){
		  xrk -= xr;
	r1[0][0] = 1;
	r1[0][1] = 0;
	r1[0][2] = 0;
	r1[1][0] = 0;
	r1[1][1] = cos(xrk);
	r1[1][2] = sin(xrk);
	r1[2][0] = 0;
	r1[2][1] = -r1[1][2];
	r1[2][2] = r1[1][1];
		  concatview( r1,orient);
	xrk = xr;
	}/* endif for xr update */
/* if yr has changed update the matrixes */
	if( yr != yrk ){
		  yrk -= yr;
	r1[0][0] = cos(yrk);
	r1[0][1] = 0;
	r1[0][2] = sin(yrk);
	r1[1][0] = 0;
	r1[1][1] = 1;
	r1[1][2] = 0;
	r1[2][0] = -r1[0][2];
	r1[2][1] = 0;
	r1[2][2] = r1[0][0];
	concatview( r1,orient);
	yrk = yr;
	}/* endif for yr update */
/* if zr has changed update the matrixes */
	if( zr != zrk ){
		  zrk -= zr;
	r1[0][0] = cos(zrk);
	r1[0][1] = sin(zrk);
	r1[0][2] = 0;
	r1[1][0] = -r1[0][1];
	r1[1][1] = r1[0][0];
	r1[1][2] = 0;
	r1[2][0] = 0;
	r1[2][1] = 0;
		  r1[2][2] = 1;
	concatview( r1,orient);
	zrk = zr;
	}/* endif for zr update */

/* apply the view transformation */
	ap = first_AMMP_ATOM();
	for( i=0; i< numatm; i++)
	{
	  x = ap->x - *xc;
	  y = ap->y - *yc;
	  z = ap->z - *zc;
	  xx = orient[0][0]*x + orient[0][1]*y + orient[0][2]*z;
	  yy = orient[1][0]*x + orient[1][1]*y + orient[1][2]*z;
	  zz = orient[2][0]*x + orient[2][1]*y + orient[2][2]*z;
	  ap->ix = -((int)(*scale*xx)-phc);
	  ap->iy = -((int)(*scale*yy)-pvc);
	  ap->iz = (int)(*scale*zz);
	  ap = ap->next;
	 }

		if(!autoscale)  return;
		xdma = -32767;
		xdmi = 32767;
		ydma = -32767;
		ydmi = 32767;
		ap = first_AMMP_ATOM();
		for( i=0; i< numatm; i++)
		{
	if( xdma < ap->ix) xdma = ap->ix;
	if( xdmi > ap->ix) xdmi = ap->ix;
	if( ydma < ap->iy) ydma = ap->iy;
		  if( ydmi > ap->iy) ydmi = ap->iy;
		  ap = ap->next;
		}
		xdma -= xdmi; ydma -= ydmi;
		if( xdma  < 2 && ydma < 2) return;
		if( xdma > ydma)
		 { *scale *= (float)phc/xdma ;  }else{ *scale *= (float)pvc/ydma;}
		 xx = *scale/100; /* set above in if (autoscale) */
		ap = first_AMMP_ATOM();
		 for( i=0; i< numatm; i++)
		 {
		 ap->ix = (ap->ix - phc)*xx + phc;
		 ap->iy = (ap->iy - pvc)*xx + pvc;
		 ap->iz *= xx;
		 ap = ap->next;
		 }
		 return;
}/* end of routine atomwalker */
int tetherwalkerpainter(orient, xr,yr,zr,xc,yc,zc,scale,
                    autoscale,autocenter,pageH,pageV,force,Hdc,slab)
HDC Hdc;
float orient[3][3];
float xr,yr,zr,*xc,*yc,*zc,*scale;
int autoscale, autocenter,force,slab;
int pageH,pageV;
{
ATOM *ap,*a_next(),*first_AMMP_ATOM();
int numatm,a_number(),i;
static float xrk=0,yrk=0,zrk=0;
float r1[3][3];/* incremental rotation matrix */
float x,y,z,xx,yy,zz,fk;
int phc,pvc;  /* center of page in device units*/
int xdmi,xdma,ydmi,ydma; /* for scaling */
int itx,ity,itz; /* for the tether values */
int get_tether_values();
int are_there_any_tethers();
HPEN oldpen,the_newpen;
	numatm = a_number();
	if( numatm <= 0) return 0;
	if( !are_there_any_tethers) return 0;
    the_newpen = CreatePen(PS_SOLID,1,RGB(255,0,255) );
	oldpen = SelectObject( Hdc,the_newpen);
		  phc = pageH/2; pvc = pageV/2;
	if( autocenter )
	{ /*this is done here because we need the center*/
	*xc = 0; *yc = 0; *zc = 0;
	ap = first_AMMP_ATOM();
	for( i=0; i< numatm; i++)
	{
	  *xc += ap->x; *yc += ap->y; *zc += ap->z;
	 ap = ap->next;
	  }
	  *xc /= numatm; *yc /= numatm; *zc /= numatm;
	  }
			 if( autoscale){
			 *scale = 100; /* too big, but will scale down */
	 }
	 if( force){ xrk = xr; yrk = yr; zrk = zr;}
/* if xr has changed update the matrixes */
	if( xr != xrk ){
		  xrk -= xr;
	r1[0][0] = 1;
	r1[0][1] = 0;
	r1[0][2] = 0;
	r1[1][0] = 0;
	r1[1][1] = cos(xrk);
	r1[1][2] = sin(xrk);
	r1[2][0] = 0;
	r1[2][1] = -r1[1][2];
	r1[2][2] = r1[1][1];
		  concatview( r1,orient);
	xrk = xr;
	}/* endif for xr update */
/* if yr has changed update the matrixes */
	if( yr != yrk ){
		  yrk -= yr;
	r1[0][0] = cos(yrk);
	r1[0][1] = 0;
	r1[0][2] = sin(yrk);
	r1[1][0] = 0;
	r1[1][1] = 1;
	r1[1][2] = 0;
	r1[2][0] = -r1[0][2];
	r1[2][1] = 0;
	r1[2][2] = r1[0][0];
	concatview( r1,orient);
	yrk = yr;
	}/* endif for yr update */
/* if zr has changed update the matrixes */
	if( zr != zrk ){
		  zrk -= zr;
	r1[0][0] = cos(zrk);
	r1[0][1] = sin(zrk);
	r1[0][2] = 0;
	r1[1][0] = -r1[0][1];
	r1[1][1] = r1[0][0];
	r1[1][2] = 0;
	r1[2][0] = 0;
	r1[2][1] = 0;
		  r1[2][2] = 1;
	concatview( r1,orient);
	zrk = zr;
	}/* endif for zr update */

/* apply the view transformation */
	ap = first_AMMP_ATOM();
	for( i=0; i< numatm; i++)
	{
	  x = ap->x - *xc;
	  y = ap->y - *yc;
	  z = ap->z - *zc;
	  xx = orient[0][0]*x + orient[0][1]*y + orient[0][2]*z;
	  yy = orient[1][0]*x + orient[1][1]*y + orient[1][2]*z;
	  zz = orient[2][0]*x + orient[2][1]*y + orient[2][2]*z;
	  ap->ix = -((int)(*scale*xx)-phc);
	  ap->iy = -((int)(*scale*yy)-pvc);
	  ap->iz = (int)(*scale*zz);
	  if( !autoscale ){
	  if( get_tether_values(ap,&x,&y,&z,&fk) ){
		x -=   *xc;
		y -=  *yc;
		z -=  *zc;
	  xx = orient[0][0]*x + orient[0][1]*y + orient[0][2]*z;
	  yy = orient[1][0]*x + orient[1][1]*y + orient[1][2]*z;
//	  zz = orient[2][0]*x + orient[2][1]*y + orient[2][2]*z;
	  itx = -((int)(*scale*xx)-phc);
	  ity = -((int)(*scale*yy)-pvc);
//	  itz = (int)(*scale*zz);
	  if( ap->iz <= slab && ap->iz >= -slab )
	  {
		  MoveToEx(Hdc,ap->ix,ap->iy,NULL);
		  LineTo( Hdc,itx,ity);
	  }
		 } /* if get_tether_values() */
	  } /* if !autoscale */
	  ap = ap->next;
	 }

	if(autoscale) {
		xdma = -32767;
		xdmi = 32767;
		ydma = -32767;
		ydmi = 32767;
		ap = first_AMMP_ATOM();
		for( i=0; i< numatm; i++)
		{
	if( xdma < ap->ix) xdma = ap->ix;
	if( xdmi > ap->ix) xdmi = ap->ix;
	if( ydma < ap->iy) ydma = ap->iy;
		  if( ydmi > ap->iy) ydmi = ap->iy;
		  ap = ap->next;
		}
		xdma -= xdmi; ydma -= ydmi;
		if( xdma  < 2 && ydma < 2) return;
		if( xdma > ydma)
		 { *scale *= (float)phc/xdma ;  }else{ *scale *= (float)pvc/ydma;}
		 xx = *scale/100; /* set above in if (autoscale) */
		ap = first_AMMP_ATOM();
		 for( i=0; i< numatm; i++)
		 {
		 ap->ix = (ap->ix - phc)*xx + phc;
		 ap->iy = (ap->iy - pvc)*xx + pvc;
		 ap->iz *= xx;
		 if( get_tether_values(ap,&x,&y,&z,&fk) ){
	  x -=   *xc;
	  y -=  *yc;
	  z -=  *zc;
	  fk = orient[0][0]*x + orient[0][1]*y + orient[0][2]*z;
	  yy = orient[1][0]*x + orient[1][1]*y + orient[1][2]*z;
//	  zz = orient[2][0]*x + orient[2][1]*y + orient[2][2]*z;
	  itx = -((int)(*scale*fk)-phc);
	  ity = -((int)(*scale*yy)-pvc);
//	  itz = (int)(*scale*zz);
	  itx = (itx-phc)*xx + phc;
	  ity = (ity -pvc)*xx + pvc;
//	  itz *= xx;
	  if( ap->iz <= slab && ap->iz >= -slab )
	  {
		  MoveToEx(Hdc,ap->ix,ap->iy,NULL);
		  LineTo( Hdc,itx,ity);
	  }
		 } /* if get_tether_values() */
		 ap = ap->next;
		 }
	} /* if autoscale */
	SelectObject(Hdc,oldpen);
	DeleteObject( the_newpen);
		 return;
}/* end of routine tetherwalker */


/* concatview makes the r1.rout product
*  for  xview = rout x ( rout is last left concat)
* and graham-scmidts it to be
* an orthogonal matrix
*/
int concatview( r1,rout)
float r1[3][3],rout[3][3];
{
float rs[3][3],a,b,c;
int i ,j,k;

for( i=0; i< 3; i++)
for( j = 0; j<3; j++)
{
rs[i][j] = 0.;
for( k=0; k<3; k++)
rs[i][j] += r1[i][k]*rout[k][j];
}

/* now orthogonalize */
/* normalize the first row */
 a = rs[0][0]*rs[0][0] + rs[1][0]*rs[1][0] + rs[2][0]*rs[2][0];
 c = 1./sqrt(a);
 rs[0][0]*=c; rs[1][0] *= c; rs[2][0] *= c;
for( i=0; i<2; i++)
{
a = rs[0][i]*rs[0][i] + rs[1][i]*rs[1][i] + rs[2][i]*rs[2][i];
for( j=i+1; j< 3; j++)
{
b = rs[0][j]*rs[0][j] + rs[1][j]*rs[1][j] + rs[2][j]*rs[2][j];
c = rs[0][j]*rs[0][i] + rs[1][j]*rs[1][i] + rs[2][j]*rs[2][i];
c = c/ sqrt(a*b);
rs[0][j] -= c*rs[0][i]; rs[1][j] -= c*rs[1][i]; rs[2][j] -= c*rs[2][i];
b = rs[0][j]*rs[0][j] + rs[1][j]*rs[1][j] + rs[2][j]*rs[2][j];
c = 1./sqrt(b);
rs[0][j] *= c; rs[1][j] *= c; rs[0][j] *= c;
}/* j */
}/* i */
rout[0][0] = rs[0][0];
rout[0][1] = rs[0][1];
rout[0][2] = rs[0][2];
rout[1][0] = rs[1][0];
rout[1][1] = rs[1][1];
rout[1][2] = rs[1][2];
rout[2][0] = rs[2][0];
rout[2][1] = rs[2][1];
rout[2][2] = rs[2][2];

}/* end of concatview */

/* bondwalker( HDC hdc)
*  walk the bonds and draw the atoms
*/
int bondwalker( hdc,slab)
HDC hdc;
int slab;
{
int i;
ATOM *ap1,*ap2,*first_AMMP_ATOM();
int bond_next();
int ix,iy,jx,jy,na,a_number();
float fz;
na = a_number();
if( na <= 0 ) return;
ap1 = first_AMMP_ATOM();
if( ap1 == NULL ) return;
for( i=0; i< na; i++)
{  if( ap1->dontuse <= 1){/* dontuse == 1 is just ap1 */
	if( ap1->iz > -slab && ap1->iz < slab){
	MoveToEx(hdc,ap1->ix-1,ap1->iy-1,NULL);
	LineTo( hdc,ap1->ix+1,ap1->iy+1);
	MoveToEx(hdc,ap1->ix+1,ap1->iy-1,NULL);
	LineTo( hdc,ap1->ix-1,ap1->iy+1);
	}  /* if in slab*/
	} /* if i dont have any bonds */
	ap1 = ap1->next;
}   /* i */

i = 0;
while(bond_next(i,&ap1,&ap2) )
{
/* check both end points */
 if( ap1->iz < -slab && ap2->iz <= -slab) goto SKIP;
 if( ap1->iz > slab && ap2->iz >= slab) goto SKIP;
 iy = ap1->iy;ix = ap1->ix;
 jy = ap2->iy;jx = ap2->ix;
 if( ap1->iz < -slab || ap1->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap2->iz)/( ap1->iz-ap2->iz);
 else fz = (float)(slab - ap2->iz)/( ap1->iz-ap2->iz);
 ix = fz*( ix -jx) +jx;
 iy = fz*( iy - jy) +jy  ;
}
 if( ap2->iz < -slab && ap2->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap1->iz)/( ap2->iz-ap1->iz);
 else fz = (float)(slab - ap1->iz)/( ap2->iz-ap1->iz);
 jx = fz*(jx - ap1->ix )+ap1->ix;
 jy = fz*(jy - ap1->iy )+ap1->iy ;

 }
MoveToEx( hdc, ix,iy, NULL);
LineTo( hdc,jx,jy);

SKIP:
i = 1;
}
if( ap1 != NULL && ap2 != NULL)
{

/* check both end points */
 if( ap1->iz < -slab && ap2->iz <= -slab) return;
 if( ap1->iz > slab && ap2->iz >= slab) return;
 iy = ap1->iy;ix = ap1->ix;
 jy = ap2->iy;jx = ap2->ix;
 if( ap1->iz < -slab || ap1->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap2->iz)/( ap1->iz-ap2->iz);
 else fz = (float)(slab - ap2->iz)/( ap1->iz-ap2->iz);
 ix = fz*( ix -jx) +jx;
 iy = fz*( iy - jy) +jy  ;
}
 if( ap2->iz < -slab && ap2->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap1->iz)/( ap2->iz-ap1->iz);
 else fz = (float)(slab - ap1->iz)/( ap2->iz-ap1->iz);
 jx = fz*(jx - ap1->ix )+ap1->ix;
 jy = fz*(jy - ap1->iy )+ap1->iy ;

 }


MoveToEx( hdc, ix,iy, NULL);
LineTo( hdc,jx,jy);

}

}/* end of routine */
distdraw( hdc, distlist,slab )
 HDC hdc;
 ATOM *distlist[20][2];
int slab;
 {
 int i;
 float x,rad;
 int inmess;
 char mess[20];
 for( i=0; i< 20; i++)
 { if( distlist[i][0] == NULL ) return;
 if( distlist[i][0]->iz >= -slab &&
	distlist[i][0]->iz <=slab) {
 if( distlist[i][1]->iz >= -slab &&
	distlist[i][1]->iz <=slab) {
 x = distlist[i][0]->x - distlist[i][1]->x;
 rad = x*x;
x = distlist[i][0]->y - distlist[i][1]->y;
 rad += x*x;
x = distlist[i][0]->z - distlist[i][1]->z;
 rad += x*x;
 rad = sqrt(rad);
 inmess = sprintf( mess,"%f",rad);
 if( inmess > 5) inmess = 5;
 MoveToEx(hdc, distlist[i][0]->ix, distlist[i][0]->iy,NULL);
 LineTo(hdc, distlist[i][1]->ix, distlist[i][1]->iy);
 TextOut( hdc,.5*(distlist[i][0]->ix+distlist[i][1]->ix),
	.5*(distlist[i][0]->iy+distlist[i][1]->iy),
	  mess,inmess);

	}   }
 }
 }/* end of routine */
angledraw( hdc,anglist,slab )
HDC hdc;
ATOM *anglist[20][3];
int slab;
{
int i ;
float x1,y1,z1,x2,y2,z2,r1,r2,r12;
int inmess;
char  mess[20];
for( i=0; i< 20; i++)
{ if( anglist[i][0] == NULL ) return;
if( anglist[i][0]->iz >= -slab && anglist[i][0]->iz <= slab)
if( anglist[i][1]->iz >= -slab && anglist[i][1]->iz <= slab)
if( anglist[i][2]->iz >= -slab && anglist[i][2]->iz <= slab)
 {
	x1 = anglist[i][0]->x - anglist[i][1]->x;
	y1 = anglist[i][0]->y - anglist[i][1]->y;
	z1 = anglist[i][0]->z - anglist[i][1]->z;
	x2 = anglist[i][2]->x - anglist[i][1]->x;
	y2 = anglist[i][2]->y - anglist[i][1]->y;
	z2 = anglist[i][2]->z - anglist[i][1]->z;
	r1 = x1*x1 + y1*y1 + z1*z1;
	r2 = x2*x2 + y2*y2 + z2*z2;
	r12 = x1*x2 + y1*y2 + z1*z2;
	r1 = sqrt( r1*r2);
	if( r1 > 1.e-3)
	{
	  r12 = r12/r1;
	  r12 = acos(r12)* 180./3.14159265;
	  inmess = sprintf( mess,"%f",r12);
	  if( inmess > 7) inmess = 7;
	  MoveToEx( hdc,anglist[i][0]->ix, anglist[i][0]->iy,NULL);
	  LineTo( hdc,anglist[i][2]->ix,anglist[i][2]->iy);
 TextOut( hdc,.5*(anglist[i][0]->ix+anglist[i][2]->ix),
	.5*(anglist[i][0]->iy+anglist[i][2]->iy),
	  mess,inmess);

	 }/* end of checking for r1 r2 */
  }  /* end of slab check*/
}     /* for( i */
}/* end of routine */
/* draw noel */
int noel_next( int , void *,void *,int*);
void draw_noel( hdc,slab,redpen,greenpen,just_violate )
HDC hdc;
int slab;
HPEN redpen,greenpen;
int just_violate;
{
	  void *np;
	int i,ix,iy,jx,jy;
	int is_met;
	float fz;
	  ATOM *ap1,*ap2;
		 i = 0;
		while( noel_next(i,&ap1,&ap2, &is_met) ){
	/* check both end points */
	if( is_met && just_violate) goto SKIP;
	if( is_met) SelectObject(hdc,greenpen);
	else SelectObject(hdc,redpen);
	
	if( ap1 == NULL )return;
	if( ap2 == NULL ) return;
 if( ap1->iz < -slab && ap2->iz <= -slab) goto SKIP;
 if( ap1->iz > slab && ap2->iz >= slab) goto SKIP;
 iy = ap1->iy;ix = ap1->ix;
 jy = ap2->iy;jx = ap2->ix;
 if( ap1->iz < -slab || ap1->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap2->iz)/( ap1->iz-ap2->iz);
 else fz = (float)(slab - ap2->iz)/( ap1->iz-ap2->iz);
 ix = fz*( ix -jx) +jx;
 iy = fz*( iy - jy) +jy  ;
}
 if( ap2->iz < -slab && ap2->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap1->iz)/( ap2->iz-ap1->iz);
 else fz = (float)(slab - ap1->iz)/( ap2->iz-ap1->iz);
 jx = fz*(jx - ap1->ix )+ap1->ix;
 jy = fz*(jy - ap1->iy )+ap1->iy ;

 }
MoveToEx( hdc, ix,iy, NULL);
LineTo( hdc,jx,jy);

SKIP:
i = 1;
}
if( ap1 == NULL || ap2 == NULL ) return ; /* can exit the above loop
either if there are no NOELS or you're on the last one */
if( is_met && just_violate) return; /* dont draw if only violations */
/* check both end points */
 if( ap1->iz < -slab && ap2->iz <= -slab) return;
 if( ap1->iz > slab && ap2->iz >= slab) return;
 iy = ap1->iy;ix = ap1->ix;
 jy = ap2->iy;jx = ap2->ix;
 if( ap1->iz < -slab || ap1->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap2->iz)/( ap1->iz-ap2->iz);
 else fz = (float)(slab - ap2->iz)/( ap1->iz-ap2->iz);
 ix = fz*( ix -jx) +jx;
 iy = fz*( iy - jy) +jy  ;
														}
 if( ap2->iz < -slab && ap2->iz > slab){
 if( ap1->iz < 0) fz = (float)(-slab - ap1->iz)/( ap2->iz-ap1->iz);
 else fz = (float)(slab - ap1->iz)/( ap2->iz-ap1->iz);
 jx = fz*(jx - ap1->ix )+ap1->ix;
 jy = fz*(jy - ap1->iy )+ap1->iy ;        }




MoveToEx( hdc, ix,iy, NULL);
LineTo( hdc,jx,jy);




}/* end of draw_noel() */
