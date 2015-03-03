#! /bin/bash

export JAVA_HOME="/Library/Java/JavaVirtualMachines/jdk1.7.0_75.jdk/Contents/Home"
export CLFFT_HOME="/Data/opencl/clFFT/src"
export OPENCL_HOME="/Developer/GPU\ Computing/OpenCL/common"

gcc CLFFT_Impl.cpp -shared -fPIC -o libJclFFT.so -framework OpenCL -I$JAVA_HOME/include -I$JAVA_HOME/include/darwin -I$CLFFT_HOME/include -I/Developer/GPU\ Computing/OpenCL/common/inc -L$CLFFT_HOME/library -lclFFT
