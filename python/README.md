
Installation of Torch for use of ANI-2x
=======================================

## Download and install GraalPy from Github

    Installation instructions are available [here](https://www.graalvm.org/latest/reference-manual/python).
    GraalPy can be downloaded from [GitHub](https://github.com/oracle/graalpython/releases).
    
    For example, the download for Linux on Intel CPUs is called:
    graalpy-24.0.1-linux-amd64.tar.gz

## Update your PATH environment variables 
 
    export GRAALPY_HOME="/iahome/m/mj/mjschnie/software/graalpy-24.0.1-linux-amd64"
    export PATH="$GRAALPY_HOME/bin:$PATH"

## Create a virtual environment within the Force Field X directory and activate it

    cd /iahome/m/mj/mjschnie/forcefieldx
    graalpy -m venv ffx_venv
    source ffx_venv/bin/activate

## Install Torch Using Pip
  
    pip install torch

## Download the ANI-2x Torch script

    wget https://ffx.biochem.uiowa.edu/ANI2x.pt.gz
    gunzip ANI2x.pt.gz

## Evaluate the ANI-2x energy & gradient
  
    ffxc ANI.groovy ../examples/water-dimer.xyz  

---

