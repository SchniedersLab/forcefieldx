
Installing GraalVM, Python and Torch
====================================

## Download GraalVM
Use the following command:
  
    bash <(curl -sL https://get.graalvm.org/jdk) -c python graalvm-ce-java17-22.3.1

Or see the following [website](https://www.graalvm.org/downloads) for other versions.

## Set the JAVA_HOME and PATH environment variables 
See the logging from the curl command or the following example:
 
    export JAVA_HOME="/iahome/m/mj/mjschnie/software/graalvm-ce-java17-22.3.1"
    export PATH="$JAVA_HOME/bin:$PATH"

## Install Python for the GraalVM.
This step can be skipped if the curl command above was used with the "-c pythom" flag.
  
    gu install python

## Create a virtual environment and activate it.

    graalpy -m venv ffx_venv
    source ffx_venv/bin/activate

## Install Torch.
  
    graalpy -m ginstall install torch

## Download the ANI-2x Torch script.

    wget https://ffx.biochem.uiowa.edu/ANI2x.pt.gz
    gunzip ANI2x.pt.gz

## Evalute the ANI-2x energy & gradient.
  
    ffxc ANI.groovy ../examples/water-dimer.xyz  

---

