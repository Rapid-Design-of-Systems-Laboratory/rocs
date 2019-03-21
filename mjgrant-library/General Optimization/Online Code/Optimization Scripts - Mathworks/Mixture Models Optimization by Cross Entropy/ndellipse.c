/*

   Return coordinates of a ellipsoïds.
	  
  
	Author : Sébastien PARIS © (sebastien.paris@lsis.org) 
		

  Exemple



	nbe     = 3000;	
    M       = randn(3 , 1 , nbe);
    S       = repmat(([3.2 , 0.2 , 0.3; 0.2 , 1.1 , 0.25 ; 0.3 , 0.25 , 2.2]) , [1 , 1 , nbe]);
	H       = [1 0 0 ; 0 1 0];
	R       = sqrt(chi2inv(0.95 , 2));
	N       = 30;

    [x , y] = ndellipse(M , S , H , R , N);
    
	plot(x , y);

  To compile : 

  mex ndellipse.c

  Myself, I use Intel CPP compiler as : 
  
  mex -f mexopts_intelamd.bat ndellipse.c


  Ver 1.0 (03/04/05)

			
			  
*/

#include <math.h>
#include "mex.h"


#define max(a , b) ((a) >= (b) ? (a) : (b))


#define PI 3.14159265358979323846

#define TINY 1.0E-14



/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/


void ndellipse(double * , double * , double * , double  , int  , 
			   double * , double *,
			   int , int ,
			   double * , double * , double*);

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/


void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )

{
	
	
	double *M , *S , *H;
	
	double R;
	
	int  N;
	
	double *x , *y;

	double *cos_angle , *sin_angle , *Stemp;
	
	int  *dimsM , *dimsS , *dimsH;
	
	int  *dimsx;
	
	int i , j , d , slice = 1 , numdimsM , numdimsS , numdimsH , numdimsx;
	
	
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	/* -------------------------- Parse INPUT  -------------------------------------- */
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	
	if(nrhs < 2)
		
	{     
		mexErrMsgTxt("At least 2 inputs argument are required for ndellipse");	
	}
	
	
	/* ----- Input 1 ----- */
	
	M           = mxGetPr(prhs[0]);
	
	numdimsM    = mxGetNumberOfDimensions(prhs[0]);
	
	dimsM       = mxGetDimensions(prhs[0]);
	
	if ( (dimsM[0] <= 1) || (dimsM[1] != 1) )
		
	{
		mexErrMsgTxt("M must be at least(d x  1 x s1 x .... x sp), d >= 2");	
	}
	
	d              = dimsM[0];
	
	for (i = 2; i <numdimsM ; i++)
		
	{
		
		slice *= dimsM[i];
		
		
	}
		
	
    /* ----- Input 2 ----- */
	
	
	S           = mxGetPr(prhs[1]);
	
	numdimsS    = mxGetNumberOfDimensions(prhs[1]);
	
	dimsS       = mxGetDimensions(prhs[1]);
	
	if ( (dimsS[0] != dimsS[1]) || (numdimsS != (numdimsM)) || (dimsS[0] != d) )
		
	{
		mexErrMsgTxt("S must be (d x d x s1 x .... x sp)");	
	}
	
	
	/* ----- Input 3 ----- */
	
	if(nrhs < 3)
		
	{

			
		H        = (double *)mxMalloc((2*d)*sizeof(double));

		for(i = 0 ; i < 2*d ; i++)

		{

			H[i] = 0.0;

		}

		H[0]     = 1.0;

		H[3]     = 1.0;

						
	}
	
	else
	{
		
		H           = mxGetPr(prhs[2]);
		
		numdimsH    = mxGetNumberOfDimensions(prhs[2]);
		
		dimsH       = mxGetDimensions(prhs[2]);
		
		if ( (numdimsH != 2) || (dimsH[0] != 2) ||  (dimsH[1] != d) )
			
		{
			mexErrMsgTxt("H must be (2 x d )");	
		}
		
		
	}
	/* ----- Input 4 ----- */
	
	if(nrhs < 4)
	{
	
		R = 2.44774683; // sqrt(chi2inv(0.95 , 2))
	}
	
	else
	
	{
		
		R         = mxGetScalar(prhs[3]);
	}
	/* ----- Input 5 ----- */
	
	if(nrhs < 5)
		
	{
		
		N = 50;
		
	}
	
	else
	
	{
	
		N         = (int) mxGetScalar(prhs[4]);
		
	}
	
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	/* -------------------------- Parse OUTPUT  ------------------------------------- */
	/*--------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------*/
	
	
	
	numdimsx  = numdimsM - 1;
	
	
	dimsx     = (int *)mxMalloc(numdimsx*sizeof(int));
	
	
	dimsx[0]  = N;
	
	
	for (i = 2; i < numdimsM ; i++)
		
	{
		
		dimsx[i - 1] = dimsM[i];
		
	}
	
	/* ----- output 1 ----- */
	
	
	plhs[0]    = mxCreateNumericArray(numdimsx, dimsx, mxDOUBLE_CLASS, mxREAL);
	
	x          = mxGetPr(plhs[0]);
	
	
	/* ----- output 2 ----- */
	
	
	plhs[1]    = mxCreateNumericArray(numdimsx, dimsx, mxDOUBLE_CLASS, mxREAL);
	
	y          = mxGetPr(plhs[1]);
	
	
	/* ----- Usefull Matrices & vectors ----- */
	
	
	cos_angle        = (double *)mxMalloc(N*sizeof(double));
	
	sin_angle        = (double *)mxMalloc(N*sizeof(double));
	
	Stemp            = (double *)mxMalloc((2*d)*sizeof(double));
	
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/
	/* ----------------------- MAIN CALL  -------------------------------------------- */
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/
	/*---------------------------------------------------------------------------------*/
	
	
	ndellipse(M , S , H , R , N , 
		x , y,
		slice, d ,
		cos_angle, sin_angle , Stemp);
	
	
	
	/*-----------------------------------------------*/
	/*-----------------------------------------------*/
	/* ------------ END of Mex File ---------------- */
	/*-----------------------------------------------*/
	/*-----------------------------------------------*/
	
	
	mxFree(dimsx);
	
	mxFree(cos_angle);

	mxFree(sin_angle);

	mxFree(Stemp);

	if(nrhs < 3)
		
	{

		mxFree(H);

	}

	
}


/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/




void ndellipse(double *M , double *S , double *H , double R , int N , 
			   double *x , double *y,
			   int slice, int d, 
			   double *cos_angle , double *sin_angle , double *Stemp)
			   
			   
{
	
	register double twoPI = 2.0*PI , pas_PI , val ;
	
	double mx , my , mtemp , stemp , a , b , c , Hx , Hy;
	
	
	double cte1 , cte2 , theta , lambda1 , lambda2 , costheta , sintheta ;
	
	double l1costheta , l2costheta , l1sintheta , l2sintheta;
	
	int i , j , v ;
	
	register int  vd , vdd , ivdd , i2 , j2 , jd , d2 = d*2 , vN , ivN;
	



	pas_PI = twoPI/(N - 1);
	
	for(i = 0 ; i < N ; i++)
		
	{
		
		cos_angle[i] = cos(i*pas_PI);
		
		sin_angle[i] = sin(i*pas_PI);
			
	}
	
	
	for (v = 0 ; v < slice ; v++)
		
	{
		
		vd  = v*d;
		
		vdd = vd*d; 
		
		vN  = v*N;
		
		
		mx  = 0.0;
		
		my  = 0.0;
		
		a   = 0.0;
		
		b   = 0.0;
		
		c   = 0.0;
		
		//Mean = H*M & H*S*H'
		
		for(i = 0 ; i < d2; i ++)
			
		{
			
			Stemp[i] = 0.0;
			
		}
		
		
		for (i = 0 ; i < d ; i++)
			
		{
			
			i2    = i*2;
			
			ivdd  = i + vdd;
			
			mtemp = M[i + vd];
			
			Hx    = H[0 + i2];
			
			Hy    = H[1 + i2];
			
						
			for (j = 0 ; j < d ; j++)
				
			{
				
				j2             = j*2;

				val            = S[j*d + ivdd];
								
				Stemp[0 + j2] += Hx*val;
				
				Stemp[1 + j2] += Hy*val;			
				
			}

			mx   += Hx*mtemp;
			
			my   += Hy*mtemp;


			a    += Stemp[0 + i2]*Hx;
			
			b    += Stemp[0 + i2]*Hy;
			
			c    += Stemp[1 + i2]*Hy;

		}

		
		cte1     = 0.5*(a + c);
		
		cte2     = 0.5*sqrt(a*a - 2*a*c + c*c + 4*b*b);
		
		lambda1  = R*sqrt(cte1 + cte2);
		
		lambda2  = R*sqrt(cte1 - cte2);


       if (fabs(b) < TINY)

	   {

		   b = TINY;

	   }

		theta       = atan((b/(cte2 - 0.5*(c - a))));
		
		costheta    = cos(theta);
		
		sintheta    = sin(theta);
		
		l1costheta  = lambda1*costheta;
		
		l1sintheta  = lambda1*sintheta;
		
		l2costheta  = lambda2*costheta;
		
		l2sintheta  = lambda2*sintheta;
		
		
		for (i = 0 ; i < N ; i++)
			
		{
					
			ivN       = i + vN;
			
			
			x[ivN]    = mx +  l1costheta*cos_angle[i] - l2sintheta*sin_angle[i]; 
			
			
			y[ivN]    = my +  l1sintheta*cos_angle[i] + l2costheta*sin_angle[i]; 
			
		}
		
}


}



