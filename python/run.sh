#! /bin/bash

# Instructions to install the GraalVM Python module: 
# https://www.graalvm.org/22.2/reference-manual/python

# Uncomment to set the path to your FFX_HOME directory.
# export FFX_HOME=/Users/mjschnie/Data/ffx-project/forcefieldx

export FFX_JVM_FLAGS="--vm.Dj3d.rend=noop --vm.Djava.awt.headless=true"

graalpython --jvm --vm.classpath=$FFX_HOME/bin/ffx-all-1.0.0-beta.jar $FFX_JVM_FLAGS ffx-script.py 

