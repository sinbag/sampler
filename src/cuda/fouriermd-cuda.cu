#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <utils_cuda.h>
#include <point2d.h>
#include <vector>
#include <thrust/complex.h>

__global__
void fouriermd_compute_kernel(double * d_outReal, double * d_outImag,
                              double* d_inPoints, double* d_inFuncVals,
                              int npts, int resolution, int ndims,
                              double twopi, float frequencyStep){
    int row = threadIdx.x + blockDim.x * blockIdx.x;
    int col = threadIdx.y + blockDim.y * blockIdx.y;

    //printf("thread indices %d , %d \n", m , l);

    if ( col >= resolution ||  row >= resolution)// || s >= points.size())
    {
        return;
    }
    else{
        double realCoeff = 0.0, imagCoeff = 0.0;
        __syncthreads();
        int halfRes = resolution * 0.5;

        if(ndims == 2){
            int wy = (row - halfRes) * frequencyStep;
            int wx = (col - halfRes) * frequencyStep;

            for(int i = 0; i < npts; i++){
                double exp = -twopi * (wx * d_inPoints[2*i] + wy * d_inPoints[2*i+1]);
                realCoeff += d_inFuncVals[i] * cosf(exp);
                imagCoeff += d_inFuncVals[i] * sinf(exp);
            }
            __syncthreads();
            int index = row * resolution + col;
            //printf("I am thread %d , %d in block %d , %d \n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
            d_outReal[index] = realCoeff;
            d_outImag[index] = imagCoeff;
        }
        else if(ndims == 3){
            int s = threadIdx.z + blockDim.z * blockIdx.z;
            int wy = (row - halfRes) * frequencyStep;
            int wx = (col - halfRes) * frequencyStep;
            int wz = (s - halfRes) * frequencyStep;

            for(int i = 0; i < npts; i++){
                double exp = -twopi * (wx * d_inPoints[3*i] + wy * d_inPoints[3*i+1] + wz * d_inPoints[3*i+2]);
                realCoeff += d_inFuncVals[i] * cosf(exp);
                imagCoeff += d_inFuncVals[i] * sinf(exp);
            }
            __syncthreads();
            int index = s*resolution*resolution + row * resolution + col;
            //printf("I am thread %d , %d in block %d , %d \n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
            d_outReal[index] = realCoeff;
            d_outImag[index] = imagCoeff;
        }
        else if(ndims == 4){
            int s = threadIdx.z + blockDim.z * blockIdx.z;
            int wy = (row - halfRes) * frequencyStep;
            int wx = (col - halfRes) * frequencyStep;
            int wz = (s - halfRes) * frequencyStep;

            __syncthreads();
            for(int u = 0; u < resolution; u++)
            //int u = halfRes;
            {
                int wu = (u - halfRes) * frequencyStep;
                //printf("\n %d, %d, %d, %d %d: ", wx, wy, wz, wu, s);
                for(int i = 0; i < npts; i++){
                    double exp = -twopi * (wx * d_inPoints[4*i] +
                                           wy * d_inPoints[4*i+1] +
                                           wz * d_inPoints[4*i+2] +
                                           wu * d_inPoints[4*i+3]);

                    realCoeff += d_inFuncVals[i] * cosf(exp);
                    imagCoeff += d_inFuncVals[i] * sinf(exp);
                }
                __syncthreads();
                int index = u*resolution*resolution*resolution +
                            s*resolution*resolution +
                            row * resolution + col;
                //printf("I am thread %d , %d in block %d , %d \n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
                d_outReal[index] = realCoeff;
                d_outImag[index] = imagCoeff;
                __syncthreads();
            }
        }
        else if(ndims == 5){
            int s = threadIdx.z + blockDim.z * blockIdx.z;
            int wy = (row - halfRes) * frequencyStep;
            int wx = (col - halfRes) * frequencyStep;
            int wz = (s - halfRes) * frequencyStep;

            for(int u = 0; u< resolution; u++){
                for(int v = 0; v < resolution; v++){
                    int wu = (u - halfRes) * frequencyStep;
                    int wv = (v - halfRes) * frequencyStep;
                    for(int i = 0; i < npts; i++){
                        double exp = -twopi * (wx * d_inPoints[5*i] +
                                               wy * d_inPoints[5*i+1] +
                                               wz * d_inPoints[5*i+2] +
                                               wu * d_inPoints[5*i+3] +
                                               wv * d_inPoints[5*i+4]);

                        realCoeff += d_inFuncVals[i] * cosf(exp);
                        imagCoeff += d_inFuncVals[i] * sinf(exp);
                    }
                    __syncthreads();
                    int index = v*resolution*resolution*resolution*resolution +
                                u*resolution*resolution*resolution +
                                s*resolution*resolution +
                                row * resolution + col;
                    //printf("I am thread %d , %d in block %d , %d \n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
                    d_outReal[index] = realCoeff;
                    d_outImag[index] = imagCoeff;
                }
            }
        }
      ///Do not uncomment!!!
      ///Division by N for real and Imaginaru coeffs is done in the
      /// C++ Host machine FourierAnalyzer code
      //realCoeffs /= N; imagCoeffs /= N;

    __syncthreads();
    // if(row == 384 && col == 384)
    //     printf("\nreal %f and imag %f : resolution: %d %d", realCoeff, imagCoeff, numRows, numCols);
    // __syncthreads();
  }
}

void cftmd_cuda(double * d_outReal, double * d_outImag,
              double* d_inPoints,  double* d_inFuncVals,
              int npts, int resolution, int ndims,
              double twopi, float frequencyStep){
    cudaDeviceSynchronize();
    //printf("Image Size: %d %d \n", numRows, numCols);

    if(ndims == 2)
    {
    const dim3 blockThreadSize(16, 16);

    int bx = (resolution + blockThreadSize.x - 1 ) / blockThreadSize.x;
    int by = (resolution + blockThreadSize.y - 1 ) / blockThreadSize.y;
    dim3 gridBlockSize( bx, by, 1);

    fouriermd_compute_kernel<<<gridBlockSize, blockThreadSize>>>(d_outReal, d_outImag, d_inPoints,
      d_inFuncVals, npts, resolution, ndims, twopi, frequencyStep);
    }
    else{
      const dim3 blockThreadSize(8,8,8);

      int bx = (resolution + blockThreadSize.x - 1 ) / blockThreadSize.x;
      int by = (resolution + blockThreadSize.y - 1 ) / blockThreadSize.y;
      int bz = (resolution + blockThreadSize.z - 1 ) / blockThreadSize.z;
      dim3 gridBlockSize( bx, by, bz);

      fouriermd_compute_kernel<<<gridBlockSize, blockThreadSize>>>(d_outReal, d_outImag, d_inPoints,
        d_inFuncVals, npts, resolution, ndims, twopi, frequencyStep);
  }
    cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
}
