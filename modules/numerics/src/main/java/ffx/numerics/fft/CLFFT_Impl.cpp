
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
    clfftSetupData fftSetup;
    clfftSetupData* fftSetupPtr;
    clfftInitSetupData(&fftSetup);
    clfftSetup(&fftSetup);
    fftSetupPtr = &fftSetup;
    return ((jlong) fftSetupPtr);
}

JNIEXPORT jlong JNICALL Java_ffx_numerics_fft_CLFFT_clfftCreateDefaultPlanNative
(JNIEnv *env, jclass object, jlong jContext, jint dimension, jint dimX, jint dimY, jint dimZ) {
    int err = 0;
    clfftDim dim;
    clfftPlanHandle planHandle;
    clfftPlanHandle* planHandlePtr;
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
    planHandlePtr = &planHandle;
    return ((jlong) planHandlePtr);
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftSetPlanPrecisionNative
(JNIEnv *env, jclass object, jlong jPlanHandle, jint precisionType) {
    clfftPlanHandle* planHandlePtr = (clfftPlanHandle*) jPlanHandle;
    clfftPlanHanlde planHandle = *planHandlePtr;
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
    return clfftSetPlanPrecision(planHandle, precision);
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftSetLayoutNative
(JNIEnv *env, jclass object, jlong jPlanHandle, jint inLayoutType, jint outLayoutType) {
    clfftPlanHandle* planHandlePtr = (clfftPlanHandle*) jPlanHandle;
    clfftPlanHanlde planHandle = *planHandlePtr;
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
    return clfftSetLayout(planHandle, inLayout, outLayout);
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftExecuteTransformNative
(JNIEnv *env, jclass object, jlong jPlanHandle, int direction, jlong jQueue, jlong jrBuffer, jlong jcBuffer) {
    clfftPlanHandle* planHandlePtr = (clfftPlanHandle*) jPlanHandle;
    clfftPlanHanlde planHandle = *planHandlePtr;
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
    clfftBakePlan(planHandle, 1, &queue, NULL, NULL);
    int ret = clfftEnqueueTransform(planHandle, dir, 1, &queue, 0, NULL, NULL, buffers, NULL, NULL);
    clFinish(queue);
    return ret;
}

JNIEXPORT jint JNICALL Java_ffx_numerics_fft_CLFFT_clfftDestroyPlanNative
(JNIEnv *env, jclass object, jlong jPlanHandle) {
    clfftPlanHandle* planHandlePtr = (clfftPlanHandle*) jPlanHandle;
    clfftPlanHanlde planHandle = *planHandlePtr;
    return clfftDestroyPlan(&planHandle);
}

JNIEXPORT void JNICALL Java_ffx_numerics_fft_CLFFT_clfftTeardownNative
(JNIEnv *env, jclass object) {
    clfftTeardown();
}
