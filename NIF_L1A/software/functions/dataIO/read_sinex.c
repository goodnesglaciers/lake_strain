/* read_sinex.c, Version 3.1
  
   Purpose: Reads data from sinex files into Matlab arrays
  
   Record of revisions:

   Date          Programmer            Description of Change
   ====          ==========            =====================
   10/03/2001    Peter Cervelli        Added support for a priori (constraint) covariance
   09/17/2001    Peter Cervelli        Fixed a bug that caused segmentation fault when trying to open non-
                                       existent files.
   04/20/2001    Peter Cervelli        Fixed a bug: invalid files had been left open.
   04/17/2001    Peter Cervelli        Major revision: (1) Generalized code to read sinex files
                                       containing arbitrary entries in the TYPE__ field. (2) 'sites'
                                       output now contains N rows, where N is the total number of
                                       data in the sinex file. (3) 'epochs' output is now Nx3, with
                                       the columns containing year:doy:second. (4) 'd' output is now
                                       Nx1.
   10/31/2000    Peter Cervelli        Added a 'types' output.
   04/23/2000    Peter Cervelli        Added check of SINEX header to verify valid SINEX file.
   04/13/2000    Peter Cervelli        Fixed a bug in determining when the 'SOLUTION/ESTIMATE' block begins.
   04/11/2000    Peter Cervelli        Original Code

*/

#include "mex.h"

FILE *fIn;

double **MatlabArray(double *M, int i, int j) {

    /* Permits indexing into a matlab array with matlab style
       row and column indices (i.e., 1 is first index, not 0).
       Note: because of the way Mathworks defines the matlab
       array, the column index comes first.
    */

     double **A;
     int n;
    
     A = (double **)mxCalloc(j+1,sizeof(double *));
    
     for(n=0;n<j;n++)
        A[n+1]=&M[n*i]-1;
    
     return A;

}

int ReadHeader(FILE *fIn) {

    /* ReadHeader queries the file header and returns the number of observations it contains. */
    
    char CheckString[BUFSIZ];
    int count;
    
    fscanf(fIn,"%s %*s %*s %*s %*s %*s %*s %*s %d",&CheckString,&count);

    if (strcmp(CheckString,"%=SNX") == 0)
        return count/1.0;
    else
        return 0;       
          
}

void ReadSinexData(double **d, double **dcov, double **epoch, char **sites, char **types, double **apcov, int n, FILE *fIn) {

     int i,j,k,c=0;
     double cov[3];

     
     char Line[BUFSIZ];
     
     while (fgets(Line,sizeof(Line),fIn)) {

        if (Line[0] == 43 && strncmp(Line, "+SOLUTION/EST",13)==0) {
            
            while (c<n) {
                fgets(Line, sizeof(Line), fIn);
                if (Line[0] == 42)
                    continue;

                sites[c]=sites[0]+(c)*5*sizeof(char);
                types[c]=types[0]+(c)*5*sizeof(char);
                                
                c++;
                
                sscanf(Line, "%*s %s %s %*s %*s %2lf %*c %3lf %*c %5lf %*s %*s %lf %*s\n",
                             types[c-1],sites[c-1],&epoch[1][c],&epoch[2][c],&epoch[3][c],&d[1][c]);
   
            }
        }
        
        if (Line[0] == 43 && strncmp(Line, "+SOLUTION/MATRIX_ESTIMATE",25)==0) {
                                       
            while (1) {                
                fgets(Line, sizeof(Line), fIn);  
                if (Line[0] == 42)
                    continue;
                if (Line[0] == 45)
                    break;
               
                c=sscanf(Line,"%d %d %lf %lf %lf",&i,&j,&cov[0],&cov[1],&cov[2]);

                for(k=0;k<c-2;k++) {
                    dcov[i][j+k]=cov[k];
                    dcov[j+k][i]=cov[k];
                }
                 
            }
        }
        
        if (Line[0] == 43 && strncmp(Line, "+SOLUTION/MATRIX_APRIORI L COVA",31)==0) {
                                       
            while (1) {                
                fgets(Line, sizeof(Line), fIn);  
                if (Line[0] == 42)
                    continue;
                if (Line[0] == 45)
                    break;
               
                c=sscanf(Line,"%d %d %lf %lf %lf",&i,&j,&cov[0],&cov[1],&cov[2]);

                for(k=0;k<c-2;k++) {
                    apcov[i][j+k]=cov[k];
                    apcov[j+k][i]=cov[k];
                }
                 
            }
        }
    }     
} 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double **d, **dcov, **apcov, **epoch;
    char *filename, **sites, **types;
    int n;

    if (nrhs != 1 || nlhs>6) {
        mexPrintf("\nread_sinex, v3.1. Reads sinex files.\n\n");
        mexPrintf("Usage: [d,dcov,epoch,sites,types,apcov]=read_sinex('filename')\n\n");
        return;
    }
    
    n= (mxGetM(prhs[0])*mxGetN(prhs[0]))+1;
    filename =mxCalloc(n,sizeof(char));
    mxGetString(prhs[0],filename,n);

    fIn=fopen(filename, "r");
    
    if (fIn != NULL)
        n=ReadHeader(fIn);
    else
        n=0;
    
    if (n == 0) {
        mexPrintf("Problem reading file: %s.\n",filename);
        for (n=0;n<nlhs;n++)
            plhs[n]=mxCreateDoubleMatrix(0,0,mxREAL);
        if (fIn != NULL)
            fclose(fIn);   
        return;
    }
    
    plhs[0]=mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(n,n,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n,3,mxREAL);
    plhs[5]=mxCreateDoubleMatrix(n,n,mxREAL);
    sites=(char **)mxCalloc(n,sizeof(char *));
    sites[0]=(char *)mxCalloc(5*n,sizeof(char));
    types=(char **)mxCalloc(n,sizeof(char *));
    types[0]=(char *)mxCalloc(5*n,sizeof(char));
    
    d=MatlabArray(mxGetPr(plhs[0]),n,1);
    dcov=MatlabArray(mxGetPr(plhs[1]),n,n);
    apcov=MatlabArray(mxGetPr(plhs[5]),n,n);
    epoch=MatlabArray(mxGetPr(plhs[2]),n,3);

    ReadSinexData(d, dcov, epoch, sites, types, apcov, n, fIn);

    fclose(fIn);    
   
    plhs[3]=mxCreateCharMatrixFromStrings(n,(const char**)sites);
    plhs[4]=mxCreateCharMatrixFromStrings(n,(const char**)types);

}
