#include "cusp/complex.h"
#define EPS 1e-50
#define complex_t cusp::complex<double>
__device__ double calc_F(const int i, const int j, const int k, double *d_data, 
						double *d_params, double *d_const, const int N)
{
%JAC1_FUNCTIONS%
}



/****************************************************************************

Optimized version of STM propagator using multiple sub-arcs and multiple arcs

*****************************************************************************/
#include <stdio.h>
extern "C" __global__ 
	void stt1_Fstar_kernel(double *d_x0, double *d_F, double *d_params, double *d_const, 
		double* d_phi0, double *d_k, double *d_k2, double *d_phi1, const double h,
		const double k_mul1, const double k_mul2, const int N, const int P, double *arcSequence, int x0Size, int part)
{
	double x0_[%NUMSTATES%];
	int j;
	int realcol, realrow;

	if(part == 1)
	{
		// Copy states for every sub-arc into shared memory for current sub-step
		// N values per thread, P threads
		if(blockIdx.x<%NUMSTATES% && blockIdx.y<%NUMSTATES% && threadIdx.x < P)
		{
			
			for(j=0;j<%NUMSTATES%;j++)
			{
				realcol = blockIdx.z*x0Size + j*P+ threadIdx.x;
				x0_[j] = d_x0[realcol];
			}
	
			realcol = P*blockIdx.y + threadIdx.x;
			d_F[blockIdx.z*N*N*P + %NUMSTATES%*P*blockIdx.x + realcol] = calc_F(blockIdx.x,blockIdx.y,arcSequence(blockIdx.z),x0_, d_params, d_const, %NUMSTATES%);
		}
	}
	else if(part == 2)
	{
		// Part 2
		double phidot;
		if(blockIdx.x<%NUMSTATES% && blockIdx.y<%NUMSTATES% && threadIdx.x < P)
		{
			// Dot product of F and phi
			phidot = 0;
			for(j=0;j<%NUMSTATES%;j++)
			{
				// realrow = j*P + threadIdx.x;
				realcol = j*P + threadIdx.x;
				phidot += d_F[blockIdx.z*N*N*P + %NUMSTATES%*P*blockIdx.x + realcol]*
							(d_phi0[blockIdx.z*N*N*P + blockIdx.y*%NUMSTATES%*P + realcol] + d_k[blockIdx.z*N*N*P + blockIdx.y*%NUMSTATES%*P + realcol]*k_mul2);
			}
			phidot = h*phidot;
			
			// Store into "k" matrix in col-major form  k_i = h*phidot
			realrow = blockIdx.x*P + threadIdx.x;
			
			d_k2[blockIdx.z*N*N*P + blockIdx.y*%NUMSTATES%*P + realrow] = phidot;			
			d_phi1[blockIdx.z*N*N*P + blockIdx.y*%NUMSTATES%*P + realrow] += k_mul1*phidot;
		}
	}
}
