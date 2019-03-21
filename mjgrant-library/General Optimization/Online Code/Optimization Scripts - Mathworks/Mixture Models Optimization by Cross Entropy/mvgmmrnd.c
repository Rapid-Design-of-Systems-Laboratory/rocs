/* mvgmmrnd.c 


  Draw samples from a mixture of multivariate Gaussian PDF

  usage: [Z , index] = mvgmmrnd(mu , sigma , p , N)
  
		
   Author : Sébastien PARIS (sebastien.paris@lsis.org) 
		  
Example
			
 mu      = cat(3 , [-5 ; -5] , [0 ; 0] ,[ 5 ; 5]);                %(d x 1 x M)
 sigma   = cat(3 , [2 0; 0 1] , [2 -.2; -.2 2] , [1 .9; .9 1]  ); %(d x d x M)
 p       = cat(3 , [0.3] , [0.1]  , [0.6]);                       %(1 x 1 x M)
 N       = 20000;
 Z       = mvgmmrnd(mu , sigma , p , N);
 [x , y] = ndellipse(mu , sigma);
 plot(Z(1 , :) , Z(2 , :) , 'k+', x , y , 'g' , 'markersize' , 2 , 'linewidth' , 2);
 hold on
 plot(reshape(mu(1 , : , :) , 1 , 3) , reshape(mu(2 , : , :) , 1 , 3) , 'r+' , 'markersize' , 6);
 hold off
				
  To compile
  
  mex -DranSHR3 mvgmmrnd.c or mex -DranKISS mvgmmrnd.c

  Myself, I use Intel CPP compiler as : 

  
  mex -DranKISS -f mexopts_intelamd.bat mvgmmrnd.c

  or

  mex -DranSHR3 -f mexopts_intelamd.bat mvgmmrnd.c

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

#define zigstep 128 // Number of Ziggurat'Steps


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


static UL jz , iz , kn[zigstep];
		
static long hz;
		
static float wn[zigstep] , fn[zigstep];



 /*--------------------------------------------------------------- */
 
 
 void randini(void);  

 void randnini(void);


 float nfix(void);

 double randn(void); 


 void matvect(double * , double * , double *, int , int , int); 


 void chol(double * , double * , int  , int); 


  
 void mvgmmrnd(double * , double * , double * , int , int , int ,
               double * , double * ,
			   double * , double * , double *);


 
 /*--------------------------------------------------------------- */
 
 
 
 void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[] )
 {
	 
	 
	 double *mu , *sigma , *p;
	 
	 double *Z , *index;

	 double *choles ,*v , *b;

 	 int  *dimsmu , *dimssigma , *dimsp;

   	 int  numdimsmu , numdimssigma , numdimsp;
	 
	 int  d , N , M;


	 /* Check input */
	 
	if(nrhs < 4)
		
	{     
		mexErrMsgTxt("At least 4 inputs argument are required for mvgmmrnd");	
	}
	 

	 /* Input 1 */
	
	
	mu           = mxGetPr(prhs[0]);
	
	numdimsmu    = mxGetNumberOfDimensions(prhs[0]);
	
	dimsmu       = mxGetDimensions(prhs[0]);
	
	if ( (numdimsmu >3) && (dimsmu[1] != 1))
		
	{
		mexErrMsgTxt("mu must be (d x 1 x M)");	
	}
	
	d              = dimsmu[0];

	M              = dimsmu[2];

	 
	 /* Input 2 */
	
	
	sigma           = mxGetPr(prhs[1]);
	
	numdimssigma    = mxGetNumberOfDimensions(prhs[1]);
	
	dimssigma       = mxGetDimensions(prhs[1]);
	
	if ( (numdimssigma > 3) && (dimssigma[0] !=d) && (dimssigma[1] != d) && (dimssigma[2] != M))
		
	{
		mexErrMsgTxt("sigma must be (d x d x M)");	
	}
		 


	 /* Input 3 */
	
	
	p               = mxGetPr(prhs[2]);
	
	numdimsp        = mxGetNumberOfDimensions(prhs[2]);
	
	dimsp           = mxGetDimensions(prhs[2]);
	
	if ( (numdimsp > 3) && (dimsp[0] !=1) && (dimsp[1] !=1) && (dimsp[2] != M))
		
	{
		mexErrMsgTxt("p must be (1 x 1 xM)");	
	}


	 /* Input 4 */


	N               = (int) mxGetScalar(prhs[3]);

	 
	 /* Output 1 */
	 
	 
	 
	 plhs[0] = mxCreateDoubleMatrix(d , N ,  mxREAL);
	 
	 Z       = mxGetPr(plhs[0]);
	 


	 /* Output 2 */
	 
	 
	 
	 plhs[1] = mxCreateDoubleMatrix(1 , N ,  mxREAL);
	 
	 index   = mxGetPr(plhs[1]);

	 
	 /* vecteur temporaire */
	 
	 
	 
	 choles     = (double *)mxMalloc((d*d*M)*sizeof(double)); 


	 v          = (double *)mxMalloc((d)*sizeof(double)); 

	 b          = (double *)mxMalloc((d)*sizeof(double)); 

	 	 
	 
	 /* Rand ~U[0,1] Seed initialization */
	 
	 
	 randini();	

    /* Initialize Ziggurat Table with zigstep steps for Normal(0,1) */

     randnini(); 

    /* Main call */
	 

     mvgmmrnd(mu , sigma , p , d , M , N ,
               Z , index ,
			   choles , v , b);

	 
	 
	 /* Free ressources */
	 
	 
	 mxFree(choles);
	 
	 mxFree(v);

	 mxFree(b);
	 
	 
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
 
 
 /* ----------------------------------------------------------------------- */
 
 void mvgmmrnd(double *mu , double *sigma , double *p , int d , int M , int N ,
	 double *Z , double *index ,
	 double *choles , double *v , double *b)
	 
	 
 {
	 
	 int i , j , jd , d2 = d*d, val;
	 
	 double temp , cP;
	 
	 
	 // Compute choles=chol(sigma)'; //
	 
	 
     chol(sigma , choles , d , M); 
	 
	 for (j = 0 ; j < N ; j++)
		 
	 {
		 
		 jd   = j*d;
		 
		 temp = rand();
		 
		 val  = 1;
		 
		 cP   = p[0];
		 
		 while( (temp > cP) && (val < M))
		 {
			 
			 cP +=p[val];
			 
			 val++;
			 
		 }
		 
		 index[j] = val;
		 
		 
		 for (i = 0; i < d ; i++)
			 
		 {
			 
			 v[i] = randn();
			 
		 }
		 
         matvect(choles , v , b , d , d , val - 1); 
		 
		 
		 for (i = 0; i < d ; i++)
			 
		 {
			 
			 Z[i + jd] = b[i] + mu[i + (val - 1)*d];
			 
		 }
			 
	 }
	  
 }



/*----------------------------------------------------------*/

void matvect(double *A , double *v , double *w, int d , int n , int val) 

/*

  w = Av, A(d x n), v(n x 1)


*/


{
	
	
	int t , i , off = val*d*d;
	
	register double temp;
	
	
	for (t = 0 ; t < d ; t++)
	{
		
		temp   = 0.0;
		
		for(i = 0 ; i < n ; i++)
		{				
			
			temp += A[t + i*d + off]*v[i];
		}
		
		w[t] = temp;
		
	}
	
	
}




/*----------------------------------------------------------*/

void chol(double *Q , double *D , int d , int M) 

{
	
	int i , j , r , d2=d*d; 
	
	int id , d1 = d - 1 , i1 , l , knnn , jd , v;
	
	double sum , p , inv_p;

	
	for (r = 0 ; r  < M ; r++)
		
	{
		
        v = r*d2;
		
		for (i = 0 ; i < d2 ; i++)
			
		{
			
			D[i + v]    = Q[i + v];
			
		}
		
		
		p           = sqrt(D[0 + v]);
		
		inv_p       = 1.0/p;
		
		D[0 + v]    = p;
		
		
		for(i = 1 ; i < d; i++)
			
		{
			
			D[d*i + v]  *= inv_p;
			
		}
		
		
		for(i = 1 ; i < d; i++)
			
		{
			id   = i*d;
			
			i1   = i - 1;
			
			sum  = D[i + id + v];    //sum = B[i][i]
			
			for(l = 0; l < i; ++l)
				
			{			
				knnn = id + l;
				
				sum -= D[knnn + v]*D[knnn + v];
			}
			
			p     = sqrt(sum);
			
			inv_p = 1.0/p;
			
			for(j = d1; j > i ; --j)
			{
				jd   = j*d;
				
				sum  = D[jd + i + v];
				
				for(l = 0; l < i ; ++l)
					
				{				
					sum   -= D[jd + l + v]*D[id + l + v];
				}
				
				
				D[jd + i + v] = sum*inv_p;
				
			}
			
			D[i + id + v] = p;
			
			for(l = d1  ; l>i1 ; l--)
			{
				
				D[l + i1*d + v] = 0.0;	 
			}
		}
		
		// D = D';
		
		for (j = 0 ; j < d ; j++)
			
		{			
			jd = j*d;
			
			for(i = j + 1 ; i < d ; i++)
				
			{
				
				D[i + jd + v]  = D[j + i*d + v];

				D[j + i*d + v] = 0.0;
				
			}
		}
		
	}
	
}



/* --------------------------------------------------------------------------- */

void randnini(void) 
{	  
	register const double m1 = 2147483648.0, m2 = 4294967296.0 ;

	register double  invm1;

	register double dn = 3.442619855899 , tn = dn , vn = 9.91256303526217e-3 , q; 
	
	int i;

	
	/* Ziggurat tables for randn */	 
	
	invm1             = 1.0/m1;
	
	q                 = vn/exp(-0.5*dn*dn);  
	
	kn[0]             = (dn/q)*m1;	  
	
	kn[1]             = 0;
		  
	wn[0]             = q*invm1;	  
	
	wn[zigstep - 1 ]  = dn*invm1;
	
	fn[0]             = 1.0;	  
	
	fn[zigstep - 1]   = exp(-0.5*dn*dn);		
	
	for(i = (zigstep - 2) ; i >= 1 ; i--)      
	{   
		dn              = sqrt(-2.*log(vn/dn + exp(-0.5*dn*dn)));          
	
		kn[i+1]         = (dn/tn)*m1;		  
		
		tn              = dn;          
		
		fn[i]           = exp(-0.5*dn*dn);          
		
		wn[i]           = dn*invm1;      
	}

}


/* --------------------------------------------------------------------------- */


float nfix(void) 
{	
	const float r = 3.442620f; 	/* The starting of the right tail */	
	
	static float x, y;
	
	for(;;)

	{
		
		x = hz*wn[iz];
					
		if(iz == 0)
		
		{	/* iz==0, handle the base strip */
			
			do
			{	
				x = -log(rand())*0.2904764;  /* .2904764 is 1/r */  
								
				y = -log(rand());			
			} 
			
			while( (y + y) < (x*x));
			
			return (hz > 0) ? (r + x) : (-r - x);	
		}
		
		if( (fn[iz] + rand()*(fn[iz-1] - fn[iz])) < ( exp(-0.5*x*x) ) ) 
		
		{
		
			return x;

		}

		
		hz = randint;		

		iz = (hz & (zigstep - 1));		
		
		if(abs(hz) < kn[iz]) 

		{
			return (hz*wn[iz]);	

		}


	}

}


/* --------------------------------------------------------------------------- */


double randn(void) 
	
{ 

		hz = randint;
		
		iz = (hz & (zigstep - 1));

		return (abs(hz) < kn[iz]) ? (hz*wn[iz]) : ( nfix() );
	
};



/* --------------------------------------------------------------------------- */
