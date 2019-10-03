#!/bin/bash

# These paths are specific to the Schnieders Lab Argon Cluster" 
export SOFT="/Dedicated/schnieders/programs"
export MVN="$SOFT/maven-repo"

# The CP is not specific to any platform, as long as the MVN variable is set above.

/Dedicated/schnieders/programs/maven-repo/org/apache/commons/commons-math3/3.6.1/commons-math3-3.6.1.jar

export CP="$MVN/org/apache/commons/commons-io/2.6/commons-io-2.6.jar:$MVN/edu/uiowa/eng/ffx/pj/1.0.1/pj-1.0.1.jar:$MVN/org/apache/commons/commons-math3/3.6.1/commons-math3-3.6.1.jar:target/numerics-1.0.0-beta.jar:$MVN/org/jcuda/jcuda/0.9.2/jcuda-0.9.2.jar:$MVN/org/jcuda/jcufft/0.9.2/jcufft-0.9.2.jar:$MVN/org/jogamp/gluegen/gluegen-rt/2.4.0/gluegen-rt-2.4.0.jar"


export LP="$SOFT/ffx/modules/numerics/src/main/java/ffx/numerics/fft:$SOFT/ffx/bin/64-bit"
# export LD_PRELOAD="/Dedicated/schnieders/programs/clMath/clFFT-nVidia/src/library/libclFFT.so"

java -cp $CP -Djava.library.path=$LP ffx.numerics.fft.Complex3DCuda 128 10
java -cp $CP -Djava.library.path=$LP ffx.numerics.fft.Complex3DCuda 256 10

# export LD_PRELOAD=""
