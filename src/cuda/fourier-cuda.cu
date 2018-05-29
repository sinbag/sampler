#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <utils_cuda.h>
#include <point2d.h>
#include <vector>
#include <thrust/complex.h>
//#include <VarianceFormulationAnalyzer.h>


// __global__ void atomic_continuous_fourier_transform(float *d_outReal, float *d_outImag,
//                                             std::vector<Point2d> &points, double twopi){
//     int myId = threadIdx.x + blockDim.x * blockIdx.x;
//     int tid = threadIdx.x;
//
//     if(m >= points.size()){
//       return;
//     }
//     else{
//       double exp = -twopi * (wx * points[myId].x + wy * points[myId].y);
//       atomicAdd(d_outReal, cos(exp));
//       atomicAdd(d_outImag, sin(exp));
//     }
// }

__global__
void fourier_compute_kernel(double * d_outReal, double * d_outImag,
                              double* d_inPoints, double* d_inFuncVals,
                              int npts, int numCols,
                              int numRows, double twopi, float frequencyStep){
  int row = threadIdx.x + blockDim.x * blockIdx.x;
  int col = threadIdx.y + blockDim.y * blockIdx.y;
  //int s = threadIdx.z + blockDim.z * blcokIdx.z;
  //printf("thread indices %d , %d \n", m , l);

  if ( col >= numCols ||  row >= numRows)// || s >= points.size())
    {
        return;
    }
 else{
      double realCoeff = 0.0, imagCoeff = 0.0;
      __syncthreads();
      int half_xRes = numCols * 0.5;
      int half_yRes = numRows * 0.5;
      int wy = (row - half_yRes) * frequencyStep;
      int wx = (col - half_xRes) * frequencyStep;
      for(int i = 0; i < npts; i++){
          double exp = -twopi * (wx * d_inPoints[2*i] + wy * d_inPoints[2*i+1]);
          realCoeff += d_inFuncVals[i] * cosf(exp);
          imagCoeff += d_inFuncVals[i] * sinf(exp);
      }
      ///Do not uncomment!!!
      ///Division by N for real and Imaginaru coeffs is done in the
      /// C++ Host machine FourierAnalyzer code
      //realCoeffs /= N; imagCoeffs /= N;

    __syncthreads();
    int index = row * numCols + col;
    //printf("I am thread %d , %d in block %d , %d \n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
    d_outReal[index] = realCoeff;
    d_outImag[index] = imagCoeff;
    __syncthreads();
    // if(row == 384 && col == 384)
    //     printf("\nreal %f and imag %f : resolution: %d %d", realCoeff, imagCoeff, numRows, numCols);
    // __syncthreads();
  }
}

void continuous_fourier_transform_cuda(double * d_outReal, double * d_outImag,
                              double* d_inPoints,  double* d_inFuncVals,
                              int npts, int numCols,
                              int numRows, double twopi, float frequencyStep){
    cudaDeviceSynchronize();
    //printf("Image Size: %d %d \n", numRows, numCols);

    const dim3 blockSize(16, 16);

    int bx = (numCols + blockSize.x - 1 ) / blockSize.x;
    int by = (numRows + blockSize.y - 1 ) / blockSize.y;
    dim3 gridSize( bx, by, 1);

    fourier_compute_kernel<<<gridSize, blockSize>>>(d_outReal, d_outImag, d_inPoints,
      d_inFuncVals, npts, numCols, numRows, twopi, frequencyStep);

    cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
}
