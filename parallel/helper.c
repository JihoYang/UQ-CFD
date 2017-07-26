#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include "helper.h"

/* ----------------------------------------------------------------------- */
/*                             auxiliary functions                         */
/* ----------------------------------------------------------------------- */
int min( int a, int b)           
{
    if( a < b ) return a;
    return b;
}

int max( int a, int b)
{
    if( a > b ) return a;
    return b;
}

double fmin( double a, double b)
{
    if( a < b ) return a;
    return b;
}

double fmax( double a, double b)
{
    if( a > b ) return a;
    return b;
}


/* ----------------------------------------------------------------------- */
/*                         local auxiliary functions                       */
/* ----------------------------------------------------------------------- */

clock_t last_timer_reset;

int min_int( const int n1, const int n2 )
{
    if( n1 < n2 ) return n1;
    return n2;
}



/* ----------------------------------------------------------------------- */
/*                             read datafile                               */
/* ----------------------------------------------------------------------- */

void errhandler( int nLine, const char *szFile, const char *szString )
{
    int err = errno;

    fprintf( ERROUT, "%s:%d Error : %s", szFile, nLine, szString );
    fprintf( ERROUT, "\n" );
    
    /* if an error within the c-library occured, an error code can be   */
    /* found in the global variable err                                 */
    if( err != 0 )
    {
	fprintf( ERROUT, "C-Lib   errno    = %d\n", err);
	fprintf( ERROUT, "C-Lib   strerror = %s\n", strerror( err ) );
    }
    exit(1);
}


/*  for comfort */
#define READ_ERROR(szMessage, szVarName, szFileName, nLine) \
  { char szTmp[80]; \
    if( nLine ) \
	sprintf( szTmp, " %s  File: %s   Variable: %s  Line: %d", szMessage, szFileName, szVarName, nLine ); \
    else \
	sprintf( szTmp, " %s  File: %s   Variable: %s ", szMessage, szFileName, szVarName); \
    ERROR( szTmp ); \
  }
    

/* --------------------------------------------------------------------------*/
/* The function searches the datafile fh for the line defining the variable  */
/* szVarName and returns the respctive string including the value of the     */
/* variable. If there's no appropriate line within the datafile, the program */
/* stops with an error messsage.                                             */
/* ATTENTION: The pointer returned refers to a static variable within the    */
/* function. To maintain the string over several program calls, it has to be */
/* copied!!!                                                                 */
/*                                                                           */
char* find_string( const char* szFileName, const char *szVarName )
{ 
    int nLine = 0;
    int i;
    FILE *fh = NULL;
    
    static char szBuffer[MAX_LINE_LENGTH];	/* containes the line read  */
                                               /* from the datafile        */

    char* szLine = szBuffer;
    char* szValue = NULL;
    char* szName = NULL;

    /* open file */
    fh = fopen( szFileName, "rt" );
    if( fh == 0 ) 
	READ_ERROR("Could not open file", szVarName, szFileName, 0);

    /* searching */
    while( ! feof(fh) )
    {
	fgets( szLine, MAX_LINE_LENGTH, fh );
	++nLine;

	/* remove comments */
	for( i = 0; i < strlen(szLine); i++)
	    if( szLine[i] == '#' )
	    {
		szLine[i] = '\0'; /* Stringende setzen */
		break;
	    }

	/* remove empty lines */
	while( isspace( (int)*szLine ) && *szLine) ++szLine;
	if( strlen( szLine ) == 0) continue; 

	/* now, the name can be extracted */
	szName = szLine;
	szValue = szLine;
	while( (isalnum( (int)*szValue ) || *szValue == '_') && *szValue) ++szValue;
	
	/* is the value for the respective name missing? */
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	*szValue = 0;		/* complete szName! at the right place */
	++szValue;
        
	/* read next line if the correct name wasn't found */
	if( strcmp( szVarName, szName)) continue;

	/* remove all leading blnkets and tabs from the value string  */
	while( isspace( (int)*szValue) ) ++szValue;
	if( *szValue == '\n' || strlen( szValue) == 0)  
	    READ_ERROR("wrong format", szName, szFileName, nLine);
	
	fclose(fh);
	return szValue;
    }  
   
    READ_ERROR("variable not found", szVarName, szFileName, nLine);
    
    return NULL;		/* dummy to satisfy the compiler  */
} 

void read_string( const char* szFileName, const char* szVarName, char*   pVariable)
{
    char* szValue = NULL;	/* string containg the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as variable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%s", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName,0);

    printf( "File: %s\t\t%s%s= %s\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              pVariable );
}

void read_int( const char* szFileName, const char* szVarName, int* pVariable)
{
    char* szValue = NULL;	/* string containing the read variable value */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%d", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName, 0);

    printf( "File: %s\t\t%s%s= %d\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              *pVariable );
}

void read_double( const char* szFileName, const char* szVarName, double* pVariable)
{
    char* szValue = NULL;	/* String mit dem eingelesenen Variablenwert */

    if( szVarName  == 0 )  ERROR("null pointer given as varable name" );
    if( szFileName == 0 )  ERROR("null pointer given as filename" );
    if( pVariable  == 0 )  ERROR("null pointer given as variable" );

    if( szVarName[0] == '*' )
	szValue = find_string( szFileName, szVarName +1 );
    else
	szValue = find_string( szFileName, szVarName );
    
    if( sscanf( szValue, "%lf", pVariable) == 0)
	READ_ERROR("wrong format", szVarName, szFileName, 0);

    printf( "File: %s\t\t%s%s= %f\n", szFileName, 
	                              szVarName,
	                              &("               "[min_int( strlen(szVarName), 15)]), 
	                              *pVariable );
}


/* ----------------------------------------------------------------------- */
/*                   write matrices to a file                              */
/* ----------------------------------------------------------------------- */

void write_matrix( const char* szFileName,       /* filename */
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch,		       /* last row */
		 double xlength,	       /* size of the geometry in */
                                               /* x-direction */
		 double ylength,	       /* size of the geometry in */
                                               /* y-direction  */
		   int fFirst ) 	       /* 0 == append, else overwrite*/
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
   int k = 0;

   if( fFirst )				/* first call of the function ? */
   {
       fh = fopen( szFileName, "w");	/* overwrite file/write new file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Outputfile %s cannot be created", szFileName );
	   ERROR( szBuff );
       }
       
/*       fprintf( fh,"%f\n%f\n%d\n%d\n%d\n%d\n", xlength, ylength, nrl, nrh, ncl, nch ); */
   }
   else
   {
       fh = fopen( szFileName ,"a");	/* append to the file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Outputfile %s cannot be opened", szFileName );
	   ERROR( szBuff );
       }
   } 

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
	   tmp[k++] = (float)m[i][j];

   fwrite( tmp, sizeof(float), nSize, fh);

   if( fclose(fh) )
   {
       char szBuff[80];
       sprintf( szBuff, "Outputfile %s cannot be closed", szFileName );
       ERROR( szBuff );
   };

   free( tmp );
}


void read_matrix( const char* szFileName,       /* filename */
		   double **m,		       /* matrix */
		   int nrl,		       /* first column */
		   int nrh,		       /* last column */
		   int ncl,		       /* first row */
		   int nch		       /* last row */
                  ) 	  
{
   int i, j;
   FILE * fh = 0;
   int nSize = (nrh-nrl+1) * (nch-ncl+1);
   float *tmp = (float *)malloc( (size_t)(nSize * sizeof(float)));
   int k = 0;

       fh = fopen( szFileName, "r");	/* overwrite file/write new file */
       if( fh == NULL )			/* opening failed ? */
       {
	   char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", szFileName );
	   ERROR( szBuff );
       }


   fread( tmp, sizeof(float), nSize, fh);

   for( j = ncl; j <= nch; j++)
       for( i = nrl; i <= nrh; i++)
	   m[i][j]=tmp[k++];

   if( fclose(fh) )
   {
       char szBuff[80];
       /*orig bug:
       sscanf( szBuff, "Inputfile %s cannot be closed", szFileName );*/
       sprintf( szBuff, "Inputfile %s cannot be closed", szFileName );
       ERROR( szBuff );
   };

   free( tmp );
}


/* ----------------------------------------------------------------------- */
/*                      general matrix functions                           */
/* ----------------------------------------------------------------------- */

/*  allocates storage for a matrix                                         */
double **matrix( int nrl, int nrh, int ncl, int nch )
{
   int i;
   int nrow = nrh - nrl + 1;	/* compute number of lines */
   int ncol = nch - ncl + 1;	/* compute number of columns */
   
   double **pArray  = (double **) malloc((size_t)( nrow * sizeof(double*)) );
   double  *pMatrix = (double *)  malloc((size_t)( nrow * ncol * sizeof( double )));

   if( pArray  == 0)  ERROR("Storage cannot be allocated");
   if( pMatrix == 0)  ERROR("Storage cannot be allocated");

   /* first entry of the array points to the value corrected by the 
      beginning of the column */
   pArray[0] = pMatrix - ncl; 

   /* compute the remaining array entries */
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   /* return the value corrected by the beginning of a line */
   return pArray - nrl;
}

//Allocates storage for a 3D matrix
double ***matrix_3d(int nrl, int nrh, int ncl, int nch, int npl, int nph){

    int i, j;
    int nrow = nrh - nrl + 1;    /* compute number of lines */
    int ncol = nch - ncl + 1;    /* compute number of columns */
    int np   = nph - npl + 1;    /* compute number of polynomials */
         
    double ***m3d = (double ***) malloc((size_t)( np * sizeof(double**)));
         
    if( m3d == 0)  ERROR("Storage cannot be allocated");
                 
    for (i = 0; i <= np; i++){
         
         m3d[i] = (double **) malloc((size_t)( nrow * sizeof(double*)));
      
         for (j = 0; j <= nrow; j++){
    
                 m3d[i][j] = (double*) malloc((size_t)( ncol * sizeof(double)));
    
         }   
    
    } 
    
    return m3d;

}

//Allocates storage for a 4D matrix
double ****matrix_4d(int nrl, int nrh, int ncl, int nch, int nzl, int nzh, int npl, int nph){

    int i, j, k;
    int nrow = nrh - nrl + 1;   //number of lines
    int ncol = nch - ncl + 1;   //number of columns
    int nz   = nzh - nzl + 1;   //number of z direction lines
    int np   = nph - npl + 1;   //number of polynomials

    double ****m4d = (double ****) malloc((size_t)( np * sizeof(double***)));
         
    if( m4d == 0)  ERROR("Storage cannot be allocated");
                 
    for (i = 0; i <= np; i++){
         
         m4d[i] = (double ***) malloc((size_t)( nrow * sizeof(double**)));

        for (j = 0; j <= nrow; j++){

            m4d[i][j] = (double**) malloc((size_t)( ncol * sizeof(double*)));

            for (k = 0; k <= ncol; k++){

                m4d[i][j][k] = (double*) malloc((size_t)( nz * sizeof(double)));

            }
    
        }   
    
    } 
    
    return m4d;

}


/* deallocates the storage of a matrix  */
void free_matrix( double **m, int nrl, int nrh, int ncl, int nch )
{
   double **pArray  = m + nrl;
   double  *pMatrix = m[nrl]+ncl;

   free( pMatrix );
   free( pArray );
}

//Deallocates the storage of a 3D matrix
void free_matrix_3d(double ***m, int nrl, int nrh, int ncl, int nch, int npl, int nph){

    double ***m3d = m;
    
    free( m3d );
 
}

//Deallocates the storage of a 4D matrix
void free_matrix_4d(double ****m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh, int npl, int nph){
    
    double ****m4d = m;
    
    free( m4d);
    
}

void init_matrix( double **m, int nrl, int nrh, int ncl, int nch, double a){

    int i, j;

    for( i = 0; i <= nrh; i++){

       for( j = ncl; j <= nch; j++){

	        m[i][j] = a;

        }

    }

}

void init_matrix_3d( double ***m, int nrl, int nrh, int ncl, int nch, int npl, int nph, double a){
    
    int i, j, k;
    
    for (i = npl; i <= nph; i++){
        
        for (j = max(0, nrl); j <= nrh; j++){
            
            for (k = max(0, ncl); k <= nch; k++){
                
                m[i][j][k] = a;
                
            }
            
        }
    }
    
}

void init_matrix_4d( double ****m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh, int npl, int nph, double a){
    
    int i, j, k, p;
    
    for ( p = npl; p <= nph; p++){
        
        for (i = nrl; i <= nrh; i++){
            
            for (j = ncl; j <= nch; j++){
                
                for (k = nzl; k <= nzh; k++){
                    
                    m[p][i][j][k] = a;
                    
                }
                
            }
            
        }
        
    }

}


/* allocates storage for a matrix */
int **imatrix( int nrl, int nrh, int ncl, int nch )
{
   int i;

   int nrow = nrh - nrl + 1;	/* compute number of rows */
   int ncol = nch - ncl + 1;	/* compute number of columns */
   
   int **pArray  = (int **) malloc((size_t)( nrow * sizeof( int* )) );
   int  *pMatrix = (int *)  malloc((size_t)( nrow * ncol * sizeof( int )));


   if( pArray  == 0)  ERROR("Storage cannot be allocated");
   if( pMatrix == 0)  ERROR("Storage cannot be allocated");

   /* first entry of the array points to the value corrected by the 
      beginning of the column */
   pArray[0] = pMatrix - ncl; 

   /* compute the remaining array entries */
   for( i = 1; i < nrow; i++ )
   {
       pArray[i] = pArray[i-1] + ncol;
   }

   /* return the value corrected by the beginning of a line */
   return pArray - nrl;
}

//Allocates storage for a integer based 3D matrix
int ***imatrix_3d( int nrl, int nrh, int ncl, int nch, int nzl, int nzh )
{
    
    int i, j;
    int nrow = nrh - nrl + 1;    /* compute number of lines */
    int ncol = nch - ncl + 1;    /* compute number of columns */
    int nz   = nzh - nzl + 1;    /* compute number of z direction elements */
    
    int ***m3d = (int ***) malloc((size_t)( nrow * sizeof(int**)));
    
    if( m3d == 0)  ERROR("Storage cannot be allocated");
    
    for (i = 0; i < nrow; i++){
        
        m3d[i] = (int **) malloc((size_t)( ncol * sizeof(int*)));
        
        for (j = 0; j < ncol; j++){
            
            m3d[i][j] = (int*) malloc((size_t)( nz * sizeof(int)));
            
        }
        
    }
    
    return m3d;
    
}

//Allocates storage for a integer based 4D matrix
int ****imatrix_4d( int nrl, int nrh, int ncl, int nch, int nzl, int nzh, int npl, int nph )
{
    
    int i, j, k;
    int nrow = nrh - nrl + 1;   //number of lines
    int ncol = nch - ncl + 1;   //number of columns
    int nz   = nzh - nzl + 1;   //number of z direction lines
    int np   = nph - npl + 1;   //number of polynomials
    
    int ****m4d = (int ****) malloc((size_t)( np * sizeof(int***)));
    
    if( m4d == 0)  ERROR("Storage cannot be allocated");
    
    for (i = 0; i < np; i++){
        
        m4d[i] = (int ***) malloc((size_t)( nrow * sizeof(int**)));
        
        for (j = 0; j < nrow; j++){
            
            m4d[i][j] = (int**) malloc((size_t)( ncol * sizeof(int*)));
            
            for (k = 0; k < ncol; k++){
                
                m4d[i][j][k] = (int*) malloc((size_t)( nz * sizeof(int)));
                
            }
            
        }   
        
    } 
    
    return m4d;
    
}


/* deallocates the storage of a matrix  */
void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch )
{
   int **pArray  = m + nrl;
   int  *pMatrix = m[nrl]+ncl;

   free( pMatrix );
   free( pArray );
}

/* deallocates the storage of a 3D matrix  */
void free_imatrix_3d(int ***m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh){
    
    int ***m3d = m;
    
    free( m3d );
    
}

//Deallocates the storage of a 4D matrix
void free_imatrix_4d(int ****m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh, int npl, int nph){
    
    int ****m4d = m;
    
    free( m4d);
    
}

void init_imatrix( int **m, int nrl, int nrh, int ncl, int nch, int a)
{
   int i,j;
   for( i = nrl; i <= nrh; i++)
       for( j = ncl; j <= nch; j++)
	   m[i][j] = a;
}

void init_imatrix_3d( int ***m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh, int a){
    
    int i, j, k;
    
    for ( i = nrl; i <= nrh; i++){
        
        for (j = ncl; j <= nch; j++){
            
            for (k = nzl; k <= nzh; k++){
                
                m[i][j][k] = a;
                
            }
            
        }
    }
    
}

void init_imatrix_4d( int ****m, int nrl, int nrh, int ncl, int nch, int nzl, int nzh, int npl, int nph, int a){
    
    int i, j, k, p;
    
    for ( p = npl; p <= nph; p++){
        
        for (i = nrl; i <= nrh; i++){
            
            for (j = ncl; j <= nch; j++){
                
                for (k = nzl; k <= nzh; k++){
                
                    m[p][i][j][k] = a;
                
                }
            
            }
            
        }
        
    }
    
}

int **read_pgm(const char *filename)
{
    FILE *input = NULL;
    char line[1024];
    int levels;
    int xsize, ysize;
    int i1, j1;
    int **pic = NULL;
    

    if ((input=fopen(filename,"rb"))==0)
    {
       char szBuff[80];
	   sprintf( szBuff, "Can not read file %s !!!", filename );
	   ERROR( szBuff );
    }

    /* check for the right "magic number" */
    if ( fread(line,1,3,input)!=3 )
    {
	    fclose(input);
	    ERROR("Error Wrong Magic field!");
    }

    /* skip the comments */
    do
    fgets(line,sizeof line,input);
    while(*line=='#');

    /* read the width and height */
    sscanf(line,"%d %d\n",&xsize,&ysize);

    printf("Image size: %d x %d\n", xsize,ysize);

    /* read # of gray levels */
    fgets(line,sizeof line,input);
    sscanf(line,"%d\n",&levels);

    /* allocate memory for image */
    pic = imatrix(0,xsize+1,0,ysize+1);
    printf("Image initialised...\n");

    /* read pixel row by row */
    for(j1=1; j1 < ysize+1; j1++)
    {
	    for (i1=1; i1 < xsize+1; i1++)
	    {
	        int byte;
            fscanf(input, "%d", &byte);

	        if (byte==EOF)
	        {
		        fclose(input);
		        ERROR("read failed");
	        }
	        else
	        {
		        pic[i1][ysize+1-j1] = byte;
		        /*printf("%d,%d: %d\n", i1,ysize+1-j1,byte);*/
	        }
	     }
    }
    for (i1 = 0; i1 < xsize+2; i1++)
    {
        pic[i1][0] = 0;
    }
    for (i1 = 0; i1 < xsize+2; i1++)
    {
        pic[i1][ysize+1] = 0;
    }
    for (j1 = 0; j1 < ysize+2; j1++)
    {
        pic[0][j1] = 0;
        pic[xsize+1][j1] = 0;
    }

    /* close file */
    fclose(input);
    
    return pic;
}
