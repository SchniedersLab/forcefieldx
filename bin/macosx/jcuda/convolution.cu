#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cufft.h>
#include <cutil_inline.h>
#include <jni.h>
#include "ffx_numerics_fft_Complex3DCuda.h"

__global__ void recipSummation(float* data, float* recip, int len)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < len) { 
       const int j = 2 * i;
       data[j]     *= recip[i];
       data[j + 1] *= recip[i];
    }
    __syncthreads();
}

/*
 * Class:     ffx_numerics_fft_Complex3DCuda
 * Method:    init
 * Signature: (III[F[F[J)I
 */
JNIEXPORT jint JNICALL Java_ffx_numerics_fft_Complex3DCuda_init
  (JNIEnv *env, jobject obj, jint nx, jint ny, jint nz, 
   jfloatArray dataArray, jfloatArray recipArray, jlongArray pointerArray) {

   // Init the CUDA device.
   cudaSetDevice( cutGetMaxGflopsDeviceId() );

   // Compute the needed device memory.
   int len = nx * ny * nz;
   int dataSize = len * 2 * sizeof(float);
   int recipSize = len * sizeof(float);

   // Create the FFT plan and allocate device memory.
   float *d_data, *d_recip;
   cufftHandle plan;
   cutilSafeCall(cudaMalloc((void**)&d_data, dataSize));
   cutilSafeCall(cudaMalloc((void**)&d_recip, recipSize));
   cufftSafeCall(cufftPlan3d(&plan, nx, ny, nz, CUFFT_C2C));

   // Save the FFT plan ID and GPU memory addresses.
   jlong *pointers = (jlong*) env->GetPrimitiveArrayCritical(pointerArray, 0);
   if (pointers == NULL) {
       return -1;
   }
   pointers[0] = plan;
   pointers[1] = (jlong) d_data;
   pointers[2] = (jlong) d_recip; 
   env->ReleasePrimitiveArrayCritical(pointerArray, pointers, 0);

   // Copy the data and reciprocal vectors to the GPU.
   jfloat *data = (jfloat*) env->GetPrimitiveArrayCritical(dataArray, 0);
   jfloat *recip = (jfloat*) env->GetPrimitiveArrayCritical(recipArray, 0);
   if (recip == NULL || data == NULL) {
       return -1;
   }
   cutilSafeCall(cudaMemcpy(d_data, data, dataSize, cudaMemcpyHostToDevice));
   cutilSafeCall(cudaMemcpy(d_recip, recip, recipSize, cudaMemcpyHostToDevice));
   env->ReleasePrimitiveArrayCritical(dataArray, data, 0);
   env->ReleasePrimitiveArrayCritical(recipArray, recip, 0);
   return 1;
}

/*
 * Class:     ffx_numerics_fft_Complex3DCuda
 * Method:    convolution
 * Signature: ([F[J)I
 */
JNIEXPORT jint JNICALL Java_ffx_numerics_fft_Complex3DCuda_convolution
  (JNIEnv *env, jobject obj, jfloatArray dataArray, jlongArray pointerArray) {

   // Get a reference to the data array
   jint len = env->GetArrayLength(dataArray) / 2; 
   jfloat *data = (jfloat*) env->GetPrimitiveArrayCritical(dataArray, 0);
   if (data == NULL) {
       return -1;
   }

   // Get the FFT plan ID and GPU memory addresses.
   jlong *pointers = (jlong*) env->GetPrimitiveArrayCritical(pointerArray, 0);
   if (pointers == NULL) {
       return -1;
   }
   cufftHandle plan = (cufftHandle) pointers[0];
   float *d_data = (float*) pointers[1];
   float *d_recip = (float*) pointers[2];

   // Compute the needed GPU memory and numbers of thread blocks.
   int threads = 512;
   int blocks = (len + threads - 1) / threads;
   int dataSize = 2 * len * sizeof(float);

   // Copy the data to the GPU and do the convolution.
   cutilSafeCall(cudaMemcpy(d_data, data, dataSize, cudaMemcpyHostToDevice));
   cufftSafeCall(cufftExecC2C(plan, (cufftComplex *) d_data, (cufftComplex *) d_data, CUFFT_FORWARD));
   recipSummation<<<blocks,threads>>>(d_data, d_recip, len);
   cufftSafeCall(cufftExecC2C(plan, (cufftComplex *) d_data, (cufftComplex *) d_data, CUFFT_INVERSE));
   cutilSafeCall(cudaMemcpy(data, d_data, dataSize, cudaMemcpyDeviceToHost));

   env->ReleasePrimitiveArrayCritical(pointerArray, pointers, 0);
   env->ReleasePrimitiveArrayCritical(dataArray, data, 0);
   return 1;
}

/*
 * Class:     ffx_numerics_fft_Complex3DCuda
 * Method:    free
 * Signature: ([J)I
 */
JNIEXPORT jint JNICALL Java_ffx_numerics_fft_Complex3DCuda_free
  (JNIEnv *env, jobject obj, jlongArray pointerArray) {

   jlong *pointers = (jlong*) env->GetPrimitiveArrayCritical(pointerArray, 0);
   if (pointers == NULL) {
       return -1;
   }
   cufftHandle plan = (cufftHandle) pointers[0];
   float *d_data = (float*) pointers[1];
   float *d_recip = (float*) pointers[2];
   env->ReleasePrimitiveArrayCritical(pointerArray, pointers, 0);

   cufftSafeCall(cufftDestroy(plan));
   cutilSafeCall(cudaFree(d_data));
   cutilSafeCall(cudaFree(d_recip));
   cudaThreadExit();
   return 1;
}

