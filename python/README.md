
Installing GraalV
==================

* Download [GraalVM] (https://www.graalvm.org/downloads)

* Set environment variables.
  export JAVA_HOME="/iahome/m/mj/mjschnie/software/graalvm-community-openjdk-17.0.7+4.1"
  export PATH="$JAVA_HOME/bin:$PATH"

* Install Python for the GraalVM.
  gu install python

* Create a virtual environment and activate it.
  graalpy -m venv ffx_venv
  source ffx_venv/bin/activate

* Install Torch.
  graalpy -m ginstall install torch

* Download the ANI-2x Torch script.
  wget https://ffx.biochem.uiowa.edu/ANI2x.pt.gz
  gunzip ANI2x.pt.gz

* Evalute the ANI-2x energy & gradient.
  ffxc ANI.groovy ../examples/water-dimer.xyz  

---

