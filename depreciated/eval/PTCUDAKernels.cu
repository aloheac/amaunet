/*
 * PTCUDAKernels.cu
 *
 *  Created on: Jan 26, 2016
 *      Author: loheac
 */

#include <cuda_runtime.h>
#include <math.h>
#include <cuComplex.h>
#include "PTSymbolicObjects.h"

__global__ void evaluatePerturbationTheoryTerm( int* dorder, int* fourierSumIndices, double* result ) {

}

__device__ cuDoubleComplex cexpf( cuDoubleComplex z ) {
	cuDoubleComplex result;
	double real = expf( z.x );  // Field x of cuDoubleComplex refers to Re{z}.
								// Field y of cuDoubleComplex refers to Im{z}.
	sincos( z.y, &result.y, &result.x );

	result.x *= real;
	result.y *= real;

	return result;
}

__device__ void evaluateDUdaggerU_1D( int* kBasisSpace, int* kBasisTime,
								      int* xBasisSpaceA, int* xBasisTimeA,
								      int* xBasisSpaceB, int* xBasisTimeB,
								      double* tau, double* mu,
								      int* NX, int* NTAU, cuDoubleComplex* result ) {

	using namespace std;

	double omega = ( 2.0 * (double)(*kBasisTime) + 1.0 ) / (double)(*NTAU);
	double k = ( 2.0 * (double)(*kBasisSpace) ) / (double)(*NX);

	cuDoubleComplex exparg;
	cuDoubleComplex D;
	cuDoubleComplex UdaggerU;
	cuDoubleComplex expresult;

	exparg.x = -(*tau) * ( ( k * k / 2.0 ) - (*mu) );
	exparg.y = -omega;
	expresult = cexpf( exparg );

	D.x = ( expresult.x * ( 1.0 + expresult.x ) + expresult.y * expresult.y ) / ( ( 1.0 + expresult.x ) * ( 1.0 + expresult.x ) + expresult.y * expresult.y );
	D.y = -expresult.x * expresult.y / ( ( 1.0 + expresult.x ) * ( 1.0 + expresult.x ) + expresult.y * expresult.y );

	exparg.x = -k * ( (*xBasisSpaceA) + (*xBasisSpaceB) );
	exparg.y = omega * ( (*xBasisTimeA) + (*xBasisTimeB) );
	UdaggerU = cexpf( exparg );

	(*result).x = D.x * UdaggerU.x - D.y * UdaggerU.y;
	(*result).y = D.x * UdaggerU.y + D.y * UdaggerU.x;
}
