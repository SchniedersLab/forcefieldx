#!/bin/sh
rm libcudafft.so
nvcc -Xcompiler -fPIC --shared ../../macosx/jcuda/convolution.cu -o libcudafft.so -I../../macosx/jcuda -I/home/schnied/Software/jdk1.6.0_17/include -I/home/schnied/Software/jdk1.6.0_17/include/linux -I/home/schnied/NVIDIA_GPU_Computing_SDK/C/common/inc -L/usr/local/cuda/lib -lcudart -lcuda -lcufft
