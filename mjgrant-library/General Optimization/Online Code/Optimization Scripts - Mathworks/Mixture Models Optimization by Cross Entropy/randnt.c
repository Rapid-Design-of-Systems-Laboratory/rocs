/*  
   Truncated Gaussian number generator

   
    y = randnt(m , sigma² , [left] , [right] , [n1] , ... , [np]);


  Author : Sébastien PARIS © (sebastien.paris@lsis.org) 


    Examples

   
	y = randnt(0 , 1 , 0 , 3 , 2 , 3);

	y = randnt(0 , 1 , [] , 3 , 2 , 3 , 2);

	y = randnt(0 , 1 , [] , [] , 2 , 3 , 2); % randn(2 , 3 , 2)


	 
	y = randnt(ones(2,2) , ones(2));  N(1 , 1)*I_{2,2}


  	To compile :
	
	
	mex -DranSHR3  randnt.c or mex -DranKISS  randnt.c  (little bit slower)

    Myself, I use Intel CPP compiler as : 

    mex -DranSHR3 -f mexopts_intelAMD.bat randnt.c

    Ver 1.0 (03/04/05)

		 
*/

#include <math.h>
#include <time.h>

#include "mex.h"



/*---------------- Basic generators definition ------------------- */

#define mix(a , b , c) \
{ \
	a -= b; a -= c; a ^= (c>>13); \
	b -= c; b -= a; b ^= (a<<8); \
	c -= a; c -= b; c ^= (b>>13); \
	a -= b; a -= c; a ^= (c>>12);  \
	b -= c; b -= a; b ^= (a<<16); \
	c -= a; c -= b; c ^= (b>>5); \
	a -= b; a -= c; a ^= (c>>3);  \
	b -= c; b -= a; b ^= (a<<10); \
	c -= a; c -= b; c ^= (b>>15); \
}



#define znew   (z = 36969*(z&65535) + (z>>16) )

#define wnew   (w = 18000*(w&65535) + (w>>16) )

#define MWC    ((znew<<16) + wnew )

#define SHR3   ( jsr ^= (jsr<<17), jsr ^= (jsr>>13), jsr ^= (jsr<<5) )

#define CONG   (jcong = 69069*jcong + 1234567)

#define KISS   ((MWC^CONG) + SHR3)




#ifdef ranKISS

#define randint KISS

#define rand() (randint*2.328306e-10)

#endif 



#ifdef ranSHR3

#define randint SHR3

#define rand() (0.5 + (signed)randint*2.328306e-10)

#endif 

/*--------------------------------------------------------------- */


typedef unsigned long UL;


/*--------------------------------------------------------------- */


static UL jsrseed = 31340134 , jsr;

#ifdef ranKISS

static UL z=362436069, w=521288629, jcong=380116160;

#endif




/*-------------------------------------------------------------------------------------------------*/


void randini(void);  


double erf(double );


double icdfn(double );


double randnt(double , double , double , double );




/*-------------------------------------------------------------------------------------------------*/


void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs,  const mxArray *prhs[])
{
	
	
	double *a , *b , *left , *right;
	
	
	double *y;
	
	
	int  *dimsa , *dimsb , *dimsleft , *dimsright ;
	
	int  *dimsy;
	
	
	int proda = 1 , prodb = 1 , prodleft = 1 , prodright = 1;
	
	int numdimsa , numdimsb , numdimsleft , numdimsright;
	
	int numdimsy;
	
	int i , N = 1;
	
	
	
	
	if((nrhs < 2))    
		
	{
		mexErrMsgTxt("Usage: y = randnt(m , sigma² , [left] , [right] , [n1] , ... , [np]);");
		
	}
	
	/* --- Input 1 ---*/
	
	a           = mxGetPr(prhs[0]);
	
	numdimsa    = mxGetNumberOfDimensions(prhs[0]);
	
	dimsa       = mxGetDimensions(prhs[0]);
	
	
	/* --- Input 2 ---*/
	
	
	b           = mxGetPr(prhs[1]);
	
	numdimsb    = mxGetNumberOfDimensions(prhs[1]);
	
	dimsb       = mxGetDimensions(prhs[1]);
	
	
	if ( numdimsb != numdimsa)
		
	{
		
		mexErrMsgTxt("m and sigma got different size");
		
	}
	
	for (i = 0 ; i < numdimsa ; i++)
		
	{
		
		proda *= dimsa[i];
		
		prodb *= dimsb[i];
		
	}
	
	if ( proda != prodb)
		
	{
		
		mexErrMsgTxt("inner dimensions of m and sigma are different");
		
	}
	

	randini();


	
	/* --- Input 3 ---*/
	
	if((nrhs < 5))   
		
	{
		
		
		if((nrhs > 2))    
			
		{
			
			if(mxIsEmpty(prhs[2]))
				
			{
				
				left        = (double *)mxMalloc(proda*sizeof(double));
				
				for (i = 0 ; i < proda ; i++)
					
				{
					
					left[i] = a[i] - 5*b[i];
					
					
				}
				
			}
			
			else
			{
				
				left        = mxGetPr(prhs[2]);
				
				numdimsleft = mxGetNumberOfDimensions(prhs[2]);
				
				dimsleft    = mxGetDimensions(prhs[2]);
				
				if ( numdimsleft != numdimsa)
					
				{
					
					mexErrMsgTxt("left don't have the same size of m and sigma");
					
				}
				
				
				for (i = 0 ; i < numdimsleft ; i++)
					
				{
					
					prodleft *= dimsleft[i];
					
					
				}
				
				if ( proda != prodleft)
					
				{
					
					mexErrMsgTxt("inner dimensions of left are not matching with m and sigma");
					
				}
				
				
			}
			
		}
		
		
		else
		{
			
			left        = (double *)mxMalloc(proda*sizeof(double));
			
			for (i = 0 ; i < proda ; i++)
				
			{
				
				left[i] = a[i] - 5*b[i];
				
				
			}
			
		}
		
		if((nrhs > 3))    
			
		{
			
			if(mxIsEmpty(prhs[3]))
				
			{
				
				right        = (double *)mxMalloc(proda*sizeof(double));
				
				for (i = 0 ; i < proda ; i++)
					
				{
					
					right[i] = a[i] + 5*b[i];
					
					
				}
				
			}
			
			else
			{
				
				right        = mxGetPr(prhs[3]);
				
				numdimsright = mxGetNumberOfDimensions(prhs[3]);
				
				dimsright    = mxGetDimensions(prhs[3]);
				
				if ( numdimsright != numdimsa)
					
				{
					
					mexErrMsgTxt("right don't have the same size of m and sigma");
					
				}
				
				
				for (i = 0 ; i < numdimsright ; i++)
					
				{
					
					prodright *= dimsright[i];
					
					
				}
				
				if ( proda != prodright)
					
				{
					
					mexErrMsgTxt("inner dimensions of right are not matching with m and sigma");
					
				}
				
				
			}
			
		}
		
		else
		{
			
			right        = (double *)mxMalloc(proda*sizeof(double));
			
			for (i = 0 ; i < proda ; i++)
				
			{
				
				right[i] = a[i] + 5*b[i];
				
				
			}
			
		}
		
		
		plhs[0]        = mxCreateNumericArray(numdimsa, dimsa, mxDOUBLE_CLASS, mxREAL);
		
		y              = mxGetPr(plhs[0]);
		
		N              = proda;
		
		
		for (i = 0 ; i < N ; i++)
			
		{
			
						y[i] = randnt(a[i] , b[i] , left[i] , right[i]);
			
			
		}
		
		if((nrhs > 2))    
			
		{	
			if(mxIsEmpty(prhs[2]))
				
			{
				
				
				mxFree(left);
				
			}
		}
		
		else
			
		{
			
			mxFree(left);
			
		}
		
		
		
		if((nrhs > 3))    
			
		{	
			
			if(mxIsEmpty(prhs[3]))
				
			{
				
				
				mxFree(right);
				
			}
			
		}
		else
			
		{
			
			mxFree(right);
			
			
		}
		
		
	}
	
	else
		
	{
		
		if ((dimsa[0] != 1) && (dimsa[1] != 1) && (dimsb[0] != 1) && (dimsb[1] != 1))
			
		{
			
			mexErrMsgTxt("m, sigma² must be scalar if randnt is called with more than 4 inputs argument");
			
		}
		
		
		if((nrhs > 2))    
			
		{
			
			if(mxIsEmpty(prhs[2]))
				
			{
				
				left        = (double *)mxMalloc(sizeof(double));
				
				
				left[0]     = a[0] - 5*b[0];
				
				
			}
			
			else
			{
				
				left        = mxGetPr(prhs[2]);
				
				numdimsleft = mxGetNumberOfDimensions(prhs[2]);
				
				dimsleft    = mxGetDimensions(prhs[2]);
				
				if ( numdimsleft != numdimsa)
					
				{
					
					mexErrMsgTxt("left don't have the same size of m and sigma");
					
				}
				
				
				for (i = 0 ; i < numdimsleft ; i++)
					
				{
					
					prodleft *= dimsleft[i];
					
					
				}
				
				if ( proda != prodleft)
					
				{
					
					mexErrMsgTxt("inner dimensions of left are not matching with m and sigma");
					
				}
				
				
			}
			
		}
		
		if((nrhs > 3))    
			
		{
			
			if(mxIsEmpty(prhs[3]))
				
			{
				
				right        = (double *)mxMalloc(sizeof(double));
				
				right[0]     = a[0] + 5*b[0];
				
			}
			
			else
			{
				
				right        = mxGetPr(prhs[3]);
				
				numdimsright = mxGetNumberOfDimensions(prhs[3]);
				
				dimsright    = mxGetDimensions(prhs[3]);
				
				if ( numdimsright != numdimsa)
					
				{
					
					mexErrMsgTxt("right don't have the same size of m and sigma");
					
				}
				
				
				for (i = 0 ; i < numdimsright ; i++)
					
				{
					
					prodright *= dimsright[i];
					
					
				}
				
				if ( proda != prodright)
					
				{
					
					mexErrMsgTxt("inner dimensions of right are not matching with m and sigma");
					
				}
				
				
			}
			
		}
		
		
		
		numdimsy = (nrhs - 4);
		
		dimsy    = (int *)mxMalloc(numdimsy*sizeof(int)); 
		
		for (i = 4 ; i < nrhs  ; i++)
			
		{
			dimsy[i - 4] = (int) mxGetScalar(prhs[i]);
			
			N           *=dimsy[i - 4];
			
		}
		
		
		
		
		plhs[0]        = mxCreateNumericArray(numdimsy, dimsy, mxDOUBLE_CLASS, mxREAL);
		
		y              = mxGetPr(plhs[0]);
		
		
		for (i = 0 ; i < N ; i++)
			
		{
			
						y[i] = randnt(a[0] , b[0] , left[0] , right[0]);
			
			
		}
		
		mxFree(dimsy);
		
		if((nrhs > 3))    
			
		{	
			
			if(mxIsEmpty(prhs[3]))
				
			{
				
				
				mxFree(right);
				
			}
			
		}
		else
			
		{
			
			mxFree(right);
			
			
		}
		
		
	}
	
	
	
	
}

/*-------------------------------------------------------------------------------------------------*/

double randnt(double a , double b , double left , double right)

{
	
	double std = sqrt(b) , lowerProb , upperProb , u;


	lowerProb = 0.5*(1 + erf( ((left - a)/(1.414213562373095*std)) ) );

	upperProb = 0.5*(1 + erf( ((right - a)/(1.414213562373095*std)) ) );

	u         = lowerProb + (upperProb-lowerProb)*rand();

	
	return  (a + icdfn(u)*std);

		
	
}




/* ----------------------------------------------------------------------- */



double icdfn(double x)

{
	/* Coefficients in rational approximations. */
	
	const double a[4] = { 0.886226899, -1.645349621,  0.914624893, -0.140543331};
	const double b[4] = {-2.118377725,  1.442710462, -0.329097515,  0.012229801};
	const double c[4] = {-1.970840454, -1.624906493,  3.429567803,  1.641345311};
	const double d[2] = { 3.543889200,  1.637067800};
	
	int  i;
	double z, x0;
    
	/* Central range. */
    
	x  = 2*x - 1;
	
	x0 = .7;
	
	/* Near end points of range. */
	
	if (x0 < x)
	{
		z = sqrt(-log((1-x)*0.5));
		
		z = (((c[3]*z+c[2])*z+c[1])*z+c[0]) / ((d[1]*z+d[0])*z+1);
	}
	else if (x<-x0)
		
	{
		
		z = sqrt(-log((1+x)*0.5));
		
		z = -(((c[3]*z+c[2])*z+c[1])*z+c[0]) / ((d[1]*z+d[0])*z+1);
	}
	
	else
		
	{
		
		z = x*x;
		
		z = x*(((a[3]*z+a[2])*z+a[1])*z+a[0]) / ((((b[3]*z+b[2])*z+b[1])*z+b[0])*z+1);
	}
	
	/* Newton-Raphson correction to full accuracy.
	Without these steps, erfinv() would be about 3 times
	faster to compute, but accurate to only about 6 digits.  */
	
	for (i = 0 ; i < 2 ; i++)
	{
		
		z = z - (erf(z) - x) / (1.128379167095513 * exp(-z*z));
	}
	
	x = 1.414213562373095*z;
	
    return(x);
}


/* ----------------------------------------------------------------------- */


double erf(double x)

{
	const double p1[5] = 
	{
		3.20937758913846947e03,
			3.77485237685302021e02,
			1.13864154151050156e02,
			3.16112374387056560e00,
			1.85777706184603153e-1
	};
	const double q1[4] = 
	{
		2.84423683343917062e03,
			1.28261652607737228e03,
			2.44024637934444173e02,
			2.36012909523441209e01
	};
	const double p2[9] = 
	{ 
		1.23033935479799725e03,
			2.05107837782607147e03,
			1.71204761263407058e03,
			8.81952221241769090e02,
			2.98635138197400131e02,
			6.61191906371416295e01,
			8.88314979438837594e00,
			5.64188496988670089e-1,
			2.15311535474403846e-8
	};
	const double q2[8] = 
	{ 
		1.23033935480374942e03,
			3.43936767414372164e03,
			4.36261909014324716e03,
			3.29079923573345963e03,
			1.62138957456669019e03,
			5.37181101862009858e02,
			1.17693950891312499e02,
			1.57449261107098347e01
	};
	const double p3[6] = 
	{
		6.58749161529837803e-4,
			1.60837851487422766e-2,
			1.25781726111229246e-1,
			3.60344899949804439e-1,
			3.05326634961232344e-1,
			1.63153871373020978e-2
	};
	const double q3[5] = 
	{ 
		2.33520497626869185e-3,
			6.05183413124413191e-2,
			5.27905102951428412e-1,
			1.87295284992346047e00,
			2.56852019228982242e00
	};
	
	int i;
	
	double xval, xx, p, q;
	
	bool NegativeValue;
	
    xval=x;
	
	if (xval<0) 
		
	{
		xval          = -xval; 
		
		NegativeValue = true;
	}
    
	else 
	{
		
		NegativeValue=false;
	}
    
	if (xval<=0.46875)
    {
		xx = xval*xval;
		
		p  = p1[4];
		
		q  = 1.0;
		
		for (i = 3 ; i>=0 ; i--) 
		{
			p = p*xx + p1[i]; 
			
			q = q*xx + q1[i];
		}
		
		xx=p/q;
		
		return(x*xx);
    }
    
	else if (xval<=4)
		
	{
		xx = xval;
		
		p  = p2[8];
		
		q  = 1.0;
		
		for (i=7 ; i>=0 ; i--) 
			
		{
			
			p = p*xx + p2[i]; 
			
			q = q*xx + q2[i];
			
		}
		
		xx = p/q;
		
		xx = exp(-xval*xval)*xx;
		
		if (NegativeValue) 
		{ 
			
			return((xx - 0.5)-0.5);
		} 
		
		else 
		{
			
			return((0.5 - xx)+0.5);
		}
    }  
    else if (xval<10)
		
    {
		
		xx = 1.0/(xval*xval);
		
		p  = p3[5];
		
		q  = 1.0;
		
		for (i = 4 ; i >= 0 ; i--) 
			
		{
			
			p = p*xx+p3[i]; 
			
			q = q*xx+q3[i];
		}
		
		xx = p/q;
		
		xx = exp(-xval*xval)*(0.7071067811865475 - xx)/(xval);
		
		
		if (mxIsNaN(xx)) 
			
		{
			
			xx = 0;
			
		}
		
		if (NegativeValue) 
		{
			
			return((xx - 0.5) - 0.5); 
			
		}
		
		else 
			
		{
			
			return((0.5 - xx) + 0.5);
			
		}
    }
    else
		
	{
		
		if (NegativeValue) 
		{
			
			return(-1); 
		}
		
		else 
			
		{
			return(1);
			
		}
	}
}



/* ----------------------------------------------------------------------- */



void randini(void) 

{
	
	/* SHR3 Seed initialization */
	
	jsrseed  = (UL) time( NULL );
	
	jsr     ^= jsrseed;
	
	
	/* KISS Seed initialization */
	
#ifdef ranKISS
	
	z        = (UL) time( NULL );
	
	w        = (UL) time( NULL ); 
	
	jcong    = (UL) time( NULL );
	
	mix(z , w , jcong);
	
#endif 
	
	
}
