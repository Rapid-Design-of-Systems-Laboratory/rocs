#include <stdio.h>
#define _num_t %DATATYPE%
#define EPS 1e-8

__device__ __forceinline__ int compare(double a, double b, double tol)
{
	if(a > b && abs(a-b) > tol)
		return 1;
	else if(a < b && abs(a-b) > tol)
		return -1;
	else
		return 0;
}
__device__ __forceinline__ bool isNaN(double a)
{
	return isnan(a);
}

#include "computeControl.h"
__device__ double calc_xdot(const int _i, const int _arcId, const int _arcType, double *d_data, double *d_const, double *d_constraints, const int N)
{
%STATE_FUNCTIONS%
}

/****************************************************************************

Optimized version of STM propagator using multiple sub-arcs and multiple arcs

*****************************************************************************/

extern "C" __global__ 
	void stt1_Fstar_kernel(double *d_x0, double *d_F, double *d_const, double *d_constraints, 
		double* d_phi0, double *d_k, double *d_k2, double *d_phi1, const double h,
		const double k_mul1, const double k_mul2, const int N, const int P, double *arcSequence, int x0Size, int part)
{
	double x0_[%NUMSTATES%+%MAXARCS%];
	// complex_t x0_[N];
	int j;
	int realcol, realrow;

	if(part == 1)
	{
		// Copy states for every sub-arc into shared memory for current sub-step
		// N values per thread, P threads
		if(blockIdx.x<N && blockIdx.y<N && threadIdx.x < P)
		{
			
			for(j=0;j<N;j++)
			{
				realcol = blockIdx.z*x0Size + j*P+ threadIdx.x;
				x0_[j] = d_x0[realcol];
			}			
			realcol = P*blockIdx.y + threadIdx.x;
			
			double fx = calc_xdot(blockIdx.x,blockIdx.z, arcSequence[blockIdx.z], x0_, d_const, d_constraints, N);
			x0_[blockIdx.y] += EPS;
			fx = calc_xdot(blockIdx.x,blockIdx.z, arcSequence[blockIdx.z], x0_, d_const, d_constraints, N)/EPS - fx/EPS;
			d_F[blockIdx.z*N*N*P + N*P*blockIdx.x + realcol] = fx;
			// if(blockIdx.z == 1 && blockIdx.y == 14 && blockIdx.x == 6)
			// {
			// 	double bank, alfa,hamiltonian;
			// 	computeControlConstraint1(d_x0,d_const,d_constraints,&bank,&alfa,&hamiltonian);
			// 	printf("%%d, %%d, %%lf, alfa = %%lf, bank = %%lf\\n",blockIdx.x, blockIdx.y, fx/EPS, alfa, bank);
			// }
			
		}
	}
	else if(part == 2)
	{
		// Part 2
		double phidot;
		if(blockIdx.x<N && blockIdx.y<N && threadIdx.x < P)
		{
			// Dot product of F and phi
			phidot = 0;
			for(j=0;j<N;j++)
			{
				// realrow = j*P + threadIdx.x;
				realcol = j*P + threadIdx.x;
				phidot += d_F[blockIdx.z*N*N*P + N*P*blockIdx.x + realcol]*(d_phi0[blockIdx.z*N*N*P + blockIdx.y*N*P + realcol] + d_k[blockIdx.z*N*N*P + blockIdx.y*N*P + realcol]*k_mul2);
			}
			phidot = h*phidot;
			// Store into "k" matrix in col-major form  k_i = h*phidot
			realrow = blockIdx.x*P + threadIdx.x;
			
			d_k2  [blockIdx.z*N*N*P + blockIdx.y*N*P + realrow] = phidot;
			d_phi1[blockIdx.z*N*N*P + blockIdx.y*N*P + realrow] += k_mul1*phidot;		
		}
	}
}
