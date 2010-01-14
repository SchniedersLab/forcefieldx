#!/bin/sh

g++ -O3 -arch x86_64 -dynamiclib -I/System/Library/Frameworks/JavaVM.framework/Headers convolution.cpp fft_execute.cpp fft_kernelstring.cpp fft_setup.cpp  -o libconvolution.jnilib -framework JavaVM -framework OpenCL
