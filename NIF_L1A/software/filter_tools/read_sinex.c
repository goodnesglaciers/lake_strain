/* read_sinex.c, Version 2.00
  
   Purpose: Reads data from sinex files into Matlab arrays

   
   Record of revisions:

   Date          Programmer            Description of Change
   ====          ==========            =====================
   04/23/2000    Peter Cervelli		   Added check of SINEX header to verify valid SINEX file.
   04/13/2000    Peter Cervelli        Fixed a bug in determining when the 'SOLUTION/ESTIMATE' block begins.
   04/11/2000    Peter Cervelli        Original Code


*/

#include "mex.h"

FILE *fIn;

/* FindNumber queries the header of the sinex file to determine how
   many observations it contains. */

int FindNumber(char filename[]) {

	char header[BUFSIZ];
	char CheckString[BUFSIZ];
	int count = 0;

	if ((fIn = fopen(filename, "r")) == NULL)
		return -1;

	fgets(header, sizeof(header), fIn);
	sscanf(header,"%s %*s %*s %*s %*s %*s %*s %*s %d",&CheckString,&count);
	fclose(fIn);

	if (strcmp(CheckString,"%=SNX") != 0)
		return 0;

	return count/3;
}

/* SinexTime2DecimalYear converts from year:doy:second format
   to decimal years.  It handles Y2K with a pivot year, 1985.
   Leap years are handled correctly. */

double SinexTime2DecimalYear(int year, int doy, int second) {

	int pivot = 85;
	double divisor=365;
	if (year > 85)
		year += 1900;
	else
		year += 2000;

	if ((year % 4 == 0 & year % 100 !=0 ) || year % 400 == 0)
		divisor=366;
	else
		divisor=365;

	return (double)year + (doy - 1 + second/86400.0) / divisor;
}

void ReadSinexData(double *pEstimate, double *pCov, double *pEpoch, char **pSites, int NumSites, char FileName[]) {

	char Line[BUFSIZ];
	int i, r, c, n, year, doy, second;
	double cov1,cov2,cov3;

	n=NumSites*3;
	
	fIn = fopen(FileName, "r");

	while (fgets(Line,sizeof(Line),fIn)) {

		if (Line[0] == 43 && strcmp(Line, "+SOLUTION/ESTIMATE\n")==0) { /* '+' is ASCII 43 */
			
			i = 0;
			fgets(Line, sizeof(Line), fIn);

			while(Line[0] != 45) {  /* '-' is ASCII 45 */

				if (Line[0] != 42) { /* '*' (comment) is ASCII 42 */

					pSites[i]=(char *)mxCalloc(5,sizeof(char));
					sscanf(Line, "%*s %*s %s %*s %*s %2d %*c %3d %*c %5d %*s %*s %lf",pSites[i],&year,&doy,&second,&pEstimate[i*3]);
					pEpoch[i]=SinexTime2DecimalYear(year,doy,second);
					fgets(Line, sizeof(Line), fIn);
					sscanf(Line, "%*s %*s %*s %*s %*s %*s %*s %*s %lf",&pEstimate[i*3+1]);
					fgets(Line, sizeof(Line), fIn);
					sscanf(Line, "%*s %*s %*s %*s %*s %*s %*s %*s %lf",&pEstimate[i*3+2]);
					i++;

				}
				fgets(Line, sizeof(Line), fIn);


			}
		}

		if (Line[0] == 43 && strcmp(Line, "+SOLUTION/MATRIX_ESTIMATE L COVA\n")==0) { /* '+' is ASCII 43 */
			
			fgets(Line, sizeof(Line), fIn);

			while(Line[0] != 45) { /* '-' is ASCII 45 */
				
				if (Line[0] != 42) { /* '*' (comment) is ASCII 42 */

					sscanf(Line, "%d %d %lf %lf %lf", &r, &c, &cov1, &cov2, &cov3);

					pCov[(c-1)*n+(r-1)]=cov1;
					pCov[(r-1)*n+(c-1)]=cov1;
					pCov[(c)*n+(r-1)]=cov2;
					pCov[(r-1)*n+(c)]=cov2;
					pCov[(c+1)*n+(r-1)]=cov3;
					pCov[(r-1)*n+(c+1)]=cov3;
					
				}

			fgets(Line, sizeof(Line), fIn);	
			}
		}
	}

	fclose(fIn);

}

				


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *pEpoch;
    double *pEstimate;
	double *pCov;

	char **pSites;
	char *pFileName;

	int FileNameLength;
	int NumSites;
	int dims[2];

	if (nrhs != 1) {
		printf("read_sinex, v2.0.\n\n");
		printf("Usage: [d,dcov,epoch,sites]=read_sinex('filename')\n");
		return;
	}
	FileNameLength = (mxGetM(prhs[0])*mxGetN(prhs[0]))+1;
	pFileName =mxCalloc(FileNameLength,sizeof(char));
	mxGetString(prhs[0],pFileName,FileNameLength);

	NumSites = FindNumber(pFileName);
	
	if (NumSites == -1)
		mexErrMsgTxt("Could not open file.");

	if (NumSites == 0)
		mexErrMsgTxt("Could not read data (invalid SINEX header).");

	
	plhs[0]=mxCreateDoubleMatrix(3,NumSites,mxREAL);
	plhs[1]=mxCreateDoubleMatrix(3*NumSites,3*NumSites,mxREAL);
	plhs[2]=mxCreateDoubleMatrix(NumSites,1,mxREAL);

	pSites=(char **)mxCalloc(NumSites,sizeof(char *));

	pEstimate=mxGetPr(plhs[0]);
	pCov=mxGetPr(plhs[1]);
	pEpoch=mxGetPr(plhs[2]);

	ReadSinexData(pEstimate,pCov,pEpoch,pSites,NumSites,pFileName);

	plhs[3]=mxCreateCharMatrixFromStrings(NumSites,(const char**)pSites);
}