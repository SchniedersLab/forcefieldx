#!/bin/sh
rm libcudafft.so

export CUDA_HOME="/usr/local/cuda"
export JAVA_HOME="/usr/lib/jvm/java-openjdk"
export CUDA_SDK="/usr/local/NVIDIA_GPU_Computing_SDK"
export PATH="$JAVA_HOME/bin:$CUDA_HOME/bin:$PATH"

nvcc --keep --compiler-options -fno-inline -Xcompiler -fPIC -c ../../macosx/jcuda/convolution.cu -I../../macosx/jcuda -I$CUDA_HOME/include -I$JAVA_HOME/include -I$JAVA_HOME/include/linux -I$CUDA_SDK/C/common/inc 

sed -i "s|__builtin_stdarg_start|__builtin_va_start|g" convolution.cu.cpp

g++ -fPIC --shared -o libcudafft.so -L/usr/lib64/nvidia -L$CUDA_HOME/lib64 -lcudart -lcuda -lcufft -I$CUDA_HOME/include convolution.cu.cpp 

rm convolution.*
