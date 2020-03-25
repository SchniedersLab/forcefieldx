#!/bin/bash

# These paths are specific to the Schnieders Lab Argon Cluster"
export SOFT="/Dedicated/schnieders/programs"
export MVN="$SOFT/maven-repo"

# These paths are do not need to be changed if the above are correct.
export JCUDA="$MVN/org/jcuda/jcuda/10.1.0/jcuda-10.1.0.jar:$MVN/org/jcuda/jcuda-natives/10.1.0/jcuda-natives-10.1.0-linux-x86_64.jar"
export JCUFFT="$MVN/org/jcuda/jcufft/10.1.0/jcufft-10.1.0.jar:$MVN/org/jcuda/jcufft-natives/10.1.0/jcufft-natives-10.1.0-linux-x86_64.jar"
export CP="$MVN/org/apache/commons/commons-io/2.6/commons-io-2.6.jar:$MVN/edu/uiowa/eng/ffx/pj/1.0.1/pj-1.0.1.jar:$MVN/org/apache/commons/commons-math3/3.6.1/commons-math3-3.6.1.jar:target/numerics-1.0.0-beta.jar:$JCUDA:$JCUFFT:$MVN/org/jogamp/gluegen/gluegen-rt/2.4.0/gluegen-rt-2.4.0.jar"

export LP="$SOFT/ffx/modules/numerics/src/main/java/ffx/numerics/fft"

java -cp $CP -Djava.library.path=$LP ffx.numerics.fft.Complex3DCuda 128 10
java -cp $CP -Djava.library.path=$LP ffx.numerics.fft.Complex3DCuda 256 10
