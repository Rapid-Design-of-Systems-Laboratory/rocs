/* corr2cov.c 

  usage: S = corr2cov(R);
  
		
  Author : Sébastien PARIS © (sebastien.paris@lsis.org) 
		  
Example

 R   = cat(3 , [2 0; 0 1] , [2 -.2; -.2 2] , [1 .9; .9 1]  ); %(d x d x M)
 S   = corr2cov(R);
  
  
  Dans matlab il faut compiler de la manière suivante  : 
  
  mex (-f mexopts_intelamd.bat) corr2cov.c

  mex (-f mexopts_intelamd.bat) corr2cov.c
				  
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
	
	
	double *R ;
	
	double *S;
	
	double *sqrtR;
	
	
	
	int  *dimsR;
	
	
	int  numdimsR;
		
	int  i , d , N = 1;
	
	
	/* Check input */
	
	if((nrhs < 1) && (nrhs > 1))
		
	{     
		mexErrMsgTxt("1 inputs argument are required for corr2cov");	
	}
	
	
	/* Input 1 */
	
	
	R            = mxGetPr(prhs[0]);
	
	numdimsR     = mxGetNumberOfDimensions(prhs[0]);
	
	dimsR        = mxGetDimensions(prhs[0]);
	
	if ( (dimsR[1] != dimsR[0]) && ( mxIsEmpty(prhs[0]) != 1) )
		
	{
		mexErrMsgTxt("R must be (d x d x n1 , .... x nl)");	
	}
	
	d              = dimsR[0];
	
	if (numdimsR > 2)
		
	{
		
		for (i = 2 ; i < numdimsR ; i++)
			
		{
			
			N *= dimsR[i];
			
		}
		
	}
	
	
	plhs[0]        = mxCreateNumericArray(numdimsR, dimsR, mxDOUBLE_CLASS, mxREAL);
	
	S              = mxGetPr(plhs[0]);
	
	/* vecteur temporaire */
	
	
	sqrtR              = (double *)mxMalloc(d*sizeof(double));
	
    /* Main call */
	
	
	corr2cov(R , S , d , N , sqrtR);
	
	
	/* Free ressources */
	
	
	mxFree(sqrtR);
}


/*----------------------------------------------------------------------------------------------*/


void corr2cov(double *R , double *S , int d , int N , 
			  double *sqrtR)
{
	
	int i , j , t , d2 = d*d , td2 , id , itd2 , jid;
	
	
	
	/* Main Loop */ 
	
	
	for(t = 0 ; t < N ; t++)
		
	{
		td2     = t*d2;
		
		
		for (i = 0 ; i < d ; i++)
			
		{
			
			sqrtR[i] = sqrt(R[i*(d + 1) + td2]);
			
		}
		
		for (i = 0 ; i < d ; i++)
			
		{
			id        = i*d + td2;
			
			itd2      = i + td2;
			
			S[id + i] = R[id + i];
			
			for (j = i + 1 ; j < d ; j++)
				
			{
				
                jid           = j + id;

				S[jid]        = R[jid]*sqrtR[i]*sqrtR[j];
				
				S[j*d + itd2] = S[jid];
				
			}
		}
	}
}


