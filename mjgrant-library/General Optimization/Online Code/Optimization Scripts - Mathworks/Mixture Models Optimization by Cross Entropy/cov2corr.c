/* corr2cov.c 

  usage: S = cov2corr(R); 
  
		
  Author : Sébastien PARIS © (sebastien.paris@lsis.org) 
		  
Example

 S   = cat(3 , [2 0; 0 1] , [2 -.2; -.2 2] , [1 .9; .9 1]  ); %(d x d x M)
 R   = corr2cov(S);
  
  
  Dans matlab il faut compiler de la manière suivante  : 
  
  mex -f mexopts_intelamd.bat cov2corr.c

  mex -f mexopts_intelamd.bat cov2corr.c
				  
*/


#include <math.h>
#include "mex.h"



/*--------------------------------------------------------------- */


void corr2cov(double * , double * , int  , int , 
			  double *);



/*--------------------------------------------------------------- */



void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[] )
{
	
	
	double *S ;
	
	double *R;
	
	double *sqrtS;
	
	
	
	int  *dimsS;
	
	
	int  numdimsS;
		
	int  i , d , N = 1;
	
	
	/* Check input */
	
	if((nrhs < 1) && (nrhs > 1))
		
	{     
		mexErrMsgTxt("1 inputs argument are required for cov2corr ");	
	}
	
	
	/* Input 1 */
	
	
	S            = mxGetPr(prhs[0]);
	
	numdimsS     = mxGetNumberOfDimensions(prhs[0]);
	
	dimsS        = mxGetDimensions(prhs[0]);
	
	if ( (dimsS[1] != dimsS[0]) && ( mxIsEmpty(prhs[0]) != 1) )
		
	{
		mexErrMsgTxt("S must be (d x d x n1 , .... x nl)");	
	}
	
	d              = dimsS[0];
	
	if (numdimsS > 2)
		
	{
		
		for (i = 2 ; i < numdimsS ; i++)
			
		{
			
			N *= dimsS[i];
			
		}
		
	}
	
	
	plhs[0]        = mxCreateNumericArray(numdimsS, dimsS, mxDOUBLE_CLASS, mxREAL);
	
	R              = mxGetPr(plhs[0]);
	
	/* vecteur temporaire */
	
	
	sqrtS          = (double *)mxMalloc(d*sizeof(double));
	
    /* Main call */
	
	
	corr2cov(S , R , d , N , sqrtS);
	
	
	/* Free ressources */
	
	
	mxFree(sqrtS);
}


/*----------------------------------------------------------------------------------------------*/


void corr2cov(double *S , double *R , int d , int N , 
			  double *sqrtS)
{
	
	int i , j , t , d2 = d*d , td2 , id , itd2 , jid;
	
	
	
	/* Main Loop */ 
	
	
	for(t = 0 ; t < N ; t++)
		
	{
		td2     = t*d2;
		
		
		for (i = 0 ; i < d ; i++)
			
		{
			
			sqrtS[i] = sqrt(S[i*(d + 1) + td2]);
			
		}
		
		for (i = 0 ; i < d ; i++)
			
		{
			id        = i*d + td2;
			
			itd2      = i + td2;
			
			R[id + i] = S[id + i];
			
			for (j = i + 1 ; j < d ; j++)
				
			{
				jid           = j + id;
				
				R[jid]        = S[jid]/(sqrtS[i]*sqrtS[j]);
				
				R[j*d + itd2] = R[jid];
				
			}
		}
	}
}


