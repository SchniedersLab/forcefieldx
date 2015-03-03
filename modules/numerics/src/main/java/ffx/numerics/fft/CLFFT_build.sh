#! /bin/bash
export JAVA_HOME="/Dedicated/schnieders/programs/jdk1.8.0_31"
export CLFFT_HOME="/Dedicated/schnieders/programs/clMath/clFFT-current/src"
export OPENCL_HOME="/opt/intel/opencl-1.2-3.2.1.16712"

icc CLFFT_Impl.cpp -shared -fPIC -o libJclFFT.so -I$JAVA_HOME/include -I$JAVA_HOME/include/linux -I$CLFFT_HOME/include -I$OPENCL_HOME/include
