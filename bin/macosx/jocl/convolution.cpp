//
// File:       convolution.cpp
//

#include <jni.h>
#include <string.h>
#include <stdlib.h>
#include <OpenCL/opencl.h>
#include "clFFT.h"
#include "fft_base_kernels.h"
#include "ffx_numerics_fft_Complex3DOpenCL.h"

typedef struct
{
	double *real;
	double *imag;
} clFFT_SplitComplexDouble;

cl_device_type globalDeviceType()
{
	char *force_cpu = getenv( "CL_DEVICE_TYPE" );
	if( force_cpu != NULL )
	{
		if( strcmp( force_cpu, "gpu" ) == 0 || strcmp( force_cpu, "CL_DEVICE_TYPE_GPU" ) == 0 )
			return CL_DEVICE_TYPE_GPU;
		else if( strcmp( force_cpu, "cpu" ) == 0 || strcmp( force_cpu, "CL_DEVICE_TYPE_CPU" ) == 0 )
			return CL_DEVICE_TYPE_CPU;
		else if( strcmp( force_cpu, "accelerator" ) == 0 || strcmp( force_cpu, "CL_DEVICE_TYPE_ACCELERATOR" ) == 0 )
			return CL_DEVICE_TYPE_ACCELERATOR;
		else if( strcmp( force_cpu, "CL_DEVICE_TYPE_DEFAULT" ) == 0 )
			return CL_DEVICE_TYPE_DEFAULT;
	}
	// default
	return CL_DEVICE_TYPE_GPU;
}

int
checkMem(clFFT_Dim3 n, cl_ulong gMemSize)
{
	cl_ulong memReq = 2;
	memReq *= n.x*n.y*n.z*sizeof(clFFT_Complex);
	memReq = memReq/1024/1024;
	if(memReq >= gMemSize)
		return -1;
	return 0;
}

void 
notify(const char *errinfo, const void *private_info, size_t cb, void *user_data)
{
    printf( "%s\n", errinfo );
}

cl_device_id       device_id;
cl_context         context;
cl_command_queue   queue;
clFFT_Dim3         n;
clFFT_SplitComplex data_split;
float*             data_recip; 
clFFT_Plan         plan;
cl_mem             data_in_real;
cl_mem             data_in_imag;
cl_mem             data_in_recip;
cl_program         recipProgram;
cl_kernel          kernel;
cl_int             n3;

/*
 * Class:     ffx_numerics_fft_Complex3DOpenCL
 * Method:    init
 * Signature: (IIILjava/nio/FloatBuffer;Ljava/nio/FloatBuffer;Ljava/nio/FloatBuffer;)I
 */
JNIEXPORT jint JNICALL Java_ffx_numerics_fft_Complex3DOpenCL_init
  (JNIEnv *env, jobject obj, jint x, jint y, jint z, 
                jobject rBuffer, jobject iBuffer, jobject recipBuffer) {

        cl_ulong gMemSize;
        cl_int err;
        unsigned int num_devices;
        n.x = x;
        n.y = y;
        n.z = z;
        n3 = x * y * z;

        cl_device_type device_type = globalDeviceType();
        if(device_type != CL_DEVICE_TYPE_GPU)
        {
                printf("Only supported on DEVICE_TYPE_GPU\n");
                return -1;
        }

        err = clGetDeviceIDs(NULL, device_type, 1, &device_id, &num_devices);
        if(err)
        {
                printf("clGetDeviceIDs failed\n");
                return -1;
        }

        context = clCreateContext(0, 1, &device_id, notify, NULL, &err);
        if(!context || err)
        {
                printf("clCreateContext failed\n");
                return -1;
        }

    	queue = clCreateCommandQueue(context, device_id, 0, &err);
    	if(!queue || err)
        {
        	printf("clCreateQueue() failed.\n");
                clReleaseContext(context);
       		return -1;
        } 

        err = clGetDeviceInfo(device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &gMemSize, NULL);
        if(err)
        {
                printf("Failed to get global mem size\n");
                clReleaseContext(context);
                clReleaseCommandQueue(queue);
                return -2;
        }

        gMemSize /= (1024*1024);

        if(checkMem(n, gMemSize)) {
                printf("Memory requirements canot be met by the available device\n");
                return -1;
        }

        err = CL_SUCCESS;
        int length = n.x * n.y * n.z;

        data_split = (clFFT_SplitComplex) { NULL, NULL };
        plan = NULL;
        data_in_real = NULL;
        data_in_imag = NULL;
        data_in_recip = NULL;
        data_split.real = (float *) env->GetDirectBufferAddress(rBuffer);
        data_split.imag = (float *) env->GetDirectBufferAddress(iBuffer);
        data_recip = (float *) env->GetDirectBufferAddress(recipBuffer);

        plan = clFFT_CreatePlan( context, n, clFFT_3D, clFFT_SplitComplexFormat, &err );
        if(!plan || err)
        {
                printf("clFFT_CreatePlan failed\n");
                return -1;
        }

        //clFFT_DumpPlan(plan, stdout);

        data_in_real = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, length*sizeof(float), data_split.real, &err);
        data_in_imag = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, length*sizeof(float), data_split.imag, &err);
        data_in_recip = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, length*sizeof(float), data_recip, &err);
        if(!data_in_real || !data_in_imag || !data_in_recip || err)
        {
                printf("clCreateBuffer failed\n");
                return -1;
        }

        recipProgram = clCreateProgramWithSource(context, 1, (const char**) &recipKernel, NULL, &err);
        if(!recipProgram || err)
        {
                printf("clCreateProgram failed for reciprocal summation kernel\n");
                return -1;
        }

        err = clBuildProgram(recipProgram, 1, &device_id, NULL, NULL, NULL);
        if(err)
        {
                printf("clBuildProgram failed for reciprocal summation kernel\n");
                return -1;
        }
        
        kernel = clCreateKernel(recipProgram, "multiply", &err);
        if(!kernel || err)
        {
                printf("clCreateKernel failed for reciprocal summation kernel\n");
                return -1;
        }

        err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &data_in_real);
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &data_in_imag);
        err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &data_in_recip);
        err |= clSetKernelArg(kernel, 3, sizeof(cl_int), &n3);

        if(err)
        {
                printf("clSetKernelArg failed for reciprocal summation kernel\n");
                return -1;
        }

        return 1;
}


/*
 * Class:     ffx_numerics_fft_Complex3DOpenCL
 * Method:    convolution
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_ffx_numerics_fft_Complex3DOpenCL_convolution
  (JNIEnv *env, jobject obj)
{	
	cl_int err = CL_SUCCESS;
        cl_mem data_out_real = data_in_real;
        cl_mem data_out_imag = data_in_imag;
        size_t local[1] = { 256 };
        size_t r = n3 % 256;
        if (r == 0) { 
            r = n3;
        } else {
            r = n3 + 256 - r;
        }
        size_t global[1] =  { r };

        int length = n.x * n.y * n.z;
	
        /* Forward + Inverse */
	err |= clEnqueueWriteBuffer(queue, data_in_real, CL_TRUE, 0, length*sizeof(float), data_split.real, 0, NULL, NULL);
	err |= clEnqueueWriteBuffer(queue, data_in_imag, CL_TRUE, 0, length*sizeof(float), data_split.imag, 0, NULL, NULL);
	err |= clEnqueueWriteBuffer(queue, data_in_recip, CL_TRUE, 0, length*sizeof(float), data_recip, 0, NULL, NULL);
	err |= clFFT_ExecutePlannar(queue, plan, 1, clFFT_Forward, data_in_real, data_in_imag, data_out_real, data_out_imag, 0, NULL, NULL);
        err |= clEnqueueNDRangeKernel(queue, kernel, 1, NULL, global, local, 0, NULL, NULL);
	err |= clFFT_ExecutePlannar(queue, plan, 1, clFFT_Inverse, data_in_real, data_in_imag, data_out_real, data_out_imag, 0, NULL, NULL);
	err |= clEnqueueReadBuffer(queue, data_out_real, CL_TRUE, 0, length*sizeof(float), data_split.real, 0, NULL, NULL);
	err |= clEnqueueReadBuffer(queue, data_out_imag, CL_TRUE, 0, length*sizeof(float), data_split.imag, 0, NULL, NULL);
	err |= clFinish(queue);
	if(err) 
	{
		printf("clFFT_Execute\n");
                return -1;
	}
        return 1;
}
	
/*
 * Class:     ffx_numerics_fft_Complex3DOpenCL
 * Method:    free
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_ffx_numerics_fft_Complex3DOpenCL_free
  (JNIEnv *env, jobject obj) {
	
	clFFT_DestroyPlan(plan);

	if(data_in_real)
		clReleaseMemObject(data_in_real);
	if(data_in_imag)
		clReleaseMemObject(data_in_imag);
	if(data_in_recip)
		clReleaseMemObject(data_in_recip);
	
	clReleaseContext(context);
	clReleaseCommandQueue(queue);
	
	return 1;
}


