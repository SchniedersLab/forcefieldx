#! /bin/bash

# This example only works with GraalPy when its part of a GraalVM.
# Otherwise, the following error will result:
# ERROR: '--jvm' is only supported when this launcher is part of a GraalVM.

# Uncomment to set the path to your FFX_HOME directory.
# export FFX_HOME=/Users/mjschnie/Data/ffx-project/forcefieldx

export FFX_JVM_FLAGS="--vm.Dj3d.rend=noop --vm.Djava.awt.headless=true"
export CLASSPATH="$FFX_HOME/lib/*"

graalpy --jvm --vm.Djava.class.path=$CLASSPATH $FFX_JVM_FLAGS ffx-script.py 

