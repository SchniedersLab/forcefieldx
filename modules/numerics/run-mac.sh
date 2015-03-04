#!/bin/bash

export MVN="/Users/mjschnie/.m2/repository"
export CP="$MVN/edu/rit/pj/pj/1.8/pj-1.8.jar:$MVN/commons-math/commons-math3/3.4.1/commons-math3-3.4.1.jar:target/numerics-1.0.0-beta.jar:$MVN/org/jogamp/jocl/jocl/2.2.0/jocl-2.2.0.jar:$MVN/org/jogamp/gluegen/gluegen-rt/2.2.0/gluegen-rt-2.2.0.jar"

export SOFT="/Data/opencl/clFFT/src/library"
export LP="/Data/ffx-project/ffx/modules/numerics/src/main/java/ffx/numerics/fft:$SOFT"

java -Xmx6G -cp $CP -Djava.library.path=$LP ffx.numerics.fft.Complex3DOpenCL 256 5

