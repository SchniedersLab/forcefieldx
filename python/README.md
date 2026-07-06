
Installation of Torch for use of ANI-2x
=======================================

## Download and install GraalPy from Github

    Installation instructions are available [here](https://www.graalvm.org/latest/reference-manual/python).
    GraalPy can be downloaded from [GitHub](https://github.com/oracle/graalpython/releases).
    
    For example, the download for Linux (Intel) is:
    graalpy3.12-25.1.3-linux-amd64.tar.gz

## Update your PATH environment variables 
 
    export GRAALPY_HOME="/iahome/m/mj/mjschnie/software/graalpy3.12-25.1.3-linux-amd64"
    export PATH="$GRAALPY_HOME/bin:$PATH"

## Create a virtual environment on Linux within the Force Field X directory and activate it
    
    cd /iahome/m/mj/mjschnie/forcefieldx
    mvn -Ppython -Dgraalpy.vfs.venvLauncher=/iahome/m/mj/mjschnie/software/graalpy3.12-25.1.3-linux-amd64/bin/graalpy
    source python-resources/venv/bin/activate

## Download the ANI-2x Torch script

    wget https://ffx.biochem.uiowa.edu/ANI2x.pt.gz
    gunzip ANI2x.pt.gz

## Evaluate the ANI-2x energy & gradient
  
    ffxc ANI ../examples/water-dimer.xyz  

---

