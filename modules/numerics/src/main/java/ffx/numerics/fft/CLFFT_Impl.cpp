
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <jni.h>
#include "clFFT.h"
#include "CL/cl.h"
#include "ffx_numerics_fft_CLFFT.h"

JNIEXPORT jlong JNICALL Java_ffx_numerics_fft_CLFFT_clfftSetupNative
(JNIEnv *env, jclass object) {
    clfftStatus_ err;
    clfftSetupData fftSetup;
    clfftSetupData* fftSetupPtr;
    err = clfftInitSetupData(&fftSetup);
    //printf(" clfftInitSetupData: %d\n", err);
    err = clfftSetup(&fftSetup);
    //printf(" clfftSetup: %d\n", err);
    fftSetupPtr = &fftSetup;
    return ((jlong) fftSetupPtr);
}

JNIEXPORT jlong JNICALL Java_ffx_numerics_fft_CLFFT_clfftCreateDefaultPlanNative
(JNIEnv *env, jclass object, jlong jContext, jint dimension, jint dimX, jint dimY, jint dimZ) {
    clfftStatus_ err;
    clfftDim dim;
    clfftPlanHandle planHandle;
    size_t clLengths[(int) dimension];
    cl_context context = (cl_context) jContext;

    switch ((int) dimension) {
        case 3:
            dim = CLFFT_3D;
            clLengths[0] = (size_t) dimX;
            clLengths[1] = (size_t) dimY;
            clLengths[2] = (size_t) dimZ;
            break;
        case 2:
            dim = CLFFT_2D;
            clLengths[0] = (size_t) dimX;
            clLengths[1] = (size_t) dimY;
            break;
        case 1:
        default:
            dim = CLFFT_1D;
            clLengths[0] = (size_t) dimX;
            break;
    }

    err = clfftCreateDefaultPlan(&planHandle, context, dim, clLengths);
    //printf(" Create %zd %d\n", planHandle, err);
    return ((jlong) planHandle);
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftSetPlanPrecisionNative
(JNIEnv *env, jclass object, jlong jPlanHandle, jint precisionType) {
    clfftStatus_ err;
    clfftPlanHandle planHandle = (clfftPlanHandle) jPlanHandle;
    clfftPrecision precision;

    switch ((int) precisionType) {
        case 0:
            precision = CLFFT_SINGLE;
            break;
        default:
        case 1:
            precision = CLFFT_DOUBLE;
            break;
    }

    err = clfftSetPlanPrecision(planHandle, precision);
    //printf(" Precision %zd %d\n", planHandle, err);
    return (int) err;
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftSetLayoutNative
(JNIEnv *env, jclass object, jlong jPlanHandle, jint inLayoutType, jint outLayoutType) {
    clfftStatus_ err;
    clfftPlanHandle planHandle = (clfftPlanHandle) jPlanHandle;
    clfftLayout inLayout;
    clfftLayout outLayout;

    switch ((int) inLayoutType) {
        default:
        case 0:
            inLayout = CLFFT_COMPLEX_INTERLEAVED;
            break;
        case 1:
            inLayout = CLFFT_COMPLEX_PLANAR;
            break;
        case 2:
            inLayout = CLFFT_REAL;
            break;
    }
    switch ((int) outLayoutType) {
        default:
        case 0:
            outLayout = CLFFT_COMPLEX_INTERLEAVED;
            break;
        case 1:
            outLayout = CLFFT_COMPLEX_PLANAR;
            break;
        case 2:
            outLayout = CLFFT_REAL;
            break;
    }
    err = clfftSetLayout(planHandle, inLayout, outLayout);
    //printf(" Layout %zd %d\n", planHandle, err);
    return (int) err;
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftExecuteTransformNative
(JNIEnv *env, jclass object, jlong jPlanHandle, int direction, jlong jQueue, jlong jrBuffer, jlong jcBuffer) {
    clfftStatus_ err;
    int status;
    clfftPlanHandle planHandle = (clfftPlanHandle) jPlanHandle;
    cl_command_queue queue = (cl_command_queue) jQueue;
    cl_mem rBuffer = (cl_mem) jrBuffer;
    cl_mem cBuffer = (cl_mem) jcBuffer;
    cl_mem buffers[] = {rBuffer, cBuffer};
    clfftDirection dir;

    switch ((int) direction) {
        default:
        case 1:
            dir = CLFFT_FORWARD;
            break;
        case -1:
            dir = CLFFT_BACKWARD;
            break;
    }
    err = clfftBakePlan(planHandle, 1, &queue, NULL, NULL);
    //printf(" Bake %zd %d\n", planHandle, err);
    err = clfftEnqueueTransform(planHandle, dir, 1, &queue, 0, NULL, NULL, buffers, NULL, NULL);
    //printf(" Enqueue %zd %d\n", planHandle, err);
    status = clFinish(queue);
    //printf(" Finish %zd %d\n", planHandle, status);
    return (int) err;
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftDestroyPlanNative
(JNIEnv *env, jclass object, jlong jPlanHandle) {
    clfftStatus_ err;
    clfftPlanHandle planHandle = (clfftPlanHandle) jPlanHandle;
    err = clfftDestroyPlan(&planHandle);
    //printf(" Destroy %zd %d\n", planHandle, err);
    return (int) err;
}

JNIEXPORT void JNICALL Java_ffx_numerics_fft_CLFFT_clfftTeardownNative
(JNIEnv *env, jclass object) {
    clfftTeardown();
}

