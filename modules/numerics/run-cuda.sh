#!/bin/bash

export MVN="/Users/mjschnie/.m2/repository"
export SOFT="/Dedicated/schnieders/programs"
export CP="$MVN/commons-io/commons-io/2.4/commons-io-2.4.jar:$MVN/edu/rit/pj/pj/1.8/pj-1.8.jar:$MVN/commons-math/commons-math3/3.4.1/commons-math3-3.4.1.jar:target/numerics-1.0.0-beta.jar:$MVN/jcuda/jcuda-all/6.5.0/jcuda-all-6.5.0.jar:$MVN/org/jogamp/gluegen/gluegen-rt/2.2.0/gluegen-rt-2.2.0.jar"

#export LP="$SOFT/ffx/modules/numerics/src/main/java/ffx/numerics/fft:$SOFT/clMath/clFFT-current/src/library"
#export LD_PRELOAD="/Dedicated/schnieders/programs/clMath/clFFT-current/src/library/libclFFT.so"

export LP="$SOFT/ffx/modules/numerics/src/main/java/ffx/numerics/fft:$SOFT/ffx/bin/64-bit"
export LD_PRELOAD="/Dedicated/schnieders/programs/clMath/clFFT-nVidia/src/library/libclFFT.so"

java -cp $CP -Djava.library.path=$LP ffx.numerics.fft.Complex3DCuda 384 10

export LD_PRELOAD=""
