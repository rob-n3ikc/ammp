
    static int (*serials)[] = NULL;
    static ATOM * (*atoms)[] = NULL;
    static int mynumber = 0;
    static int myalloc = 0;
static void atom_add( ATOM *ap)
{
	int (*stemp)[];
	ATOM * (*atemp)[];
	int i;
	if( serials == NULL)
	{
		serials = (int (*)[])malloc( 1000*sizeof(int));
		atoms = (ATOM * (*)[])malloc( 1000*sizeof(ATOM*));
		myalloc = 1000;
	}
	if( mynumber+1 >= myalloc)
	{
		stemp = (int  (*)[])malloc( (1000+myalloc)*sizeof(int));
		atemp = (ATOM*  (*)[])malloc( (1000+myalloc)*sizeof(ATOM*));
		for( i=0; i< myalloc; i++)
		{
			(*stemp)[i] = (*serials)[i];
			(*atemp)[i] = (*atoms)[i];
		}
		for( i=myalloc; i< myalloc+1000; i++)
			(*stemp)[i] = -1;/* -1 is not allowed as a serial */
		free(atoms);
		free(serials);
		serials = stemp;
		atoms = atemp;
		myalloc += 1000;
	
	}	
		(*serials)[mynumber] = ap->serial;
		(*atoms)[mynumber] = ap;
		mynumber += 1;


}
/* function a_m_serial( serial )
* returns NULL on error or returns the address of the ATOM
* which matches serial
* cute?
*/

ATOM *a_m_serial( serial )
int serial;
{

	int i;
	if( serials == NULL)
		return NULL;

	for( i=0; i< mynumber; i++)
	if( (*serials)[i] == serial) return (*atoms)[i];		

	return NULL;
}
