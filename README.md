
Force Field X 
=============

Force Field X is an atomic resolution molecular modeling application that targets open research questions in the areas of:
* predicting the structure, thermodynamic stability and solubility of organic polymer crystals,
* predicting the effect of missense mutations on protein structure, thermodynamics and molecular phenotype
* computational design of biomolecules in both soluble and crystalline environments

Please see the Force Field X [website](http://ffx.biochem.uiowa.edu) for more details.

[![Build Status](https://travis-ci.org/mjschnie/ffx.svg?branch=master)](https://travis-ci.org/mjschnie/ffx/)
[![codecov](https://codecov.io/gh/mjschnie/ffx/branch/master/graph/badge.svg)](https://codecov.io/gh/mjschnie/ffx)

---

## Clone Force Field X

To clone the Force Field X source using GIT:

    git clone git@github.com:mjschnie/ffx.git ffx

In the future, your clone of Force Field X can be updated to the latest version using the command:

    git pull origin master

---

## Build Using Maven

A Maven project file (pom.xml) is provided to build Force Field X on any platform. After cloing the Force Field X git repository, change directoies into the base project directory. Then execute:

    mvn

This requires Maven v. 3.2 or later to be installed with its bin directory included in your $PATH environment variable. The first time this command is executed, Maven will download build dependencies and Force Field X runtime dependecies. Future executions are quicker. Force Field X will self-test its modules and report failures. Only code that passes all testing should be pushed to the GitHub repository, so if any test fails it may be due to a local configuration issue. To execute the tests:

    mvn -DskipTests=false

Additional tests, ordinarily skipped due to length of running them (~15 minutes on a single core of a 2013 CPU) can be accessed via the ffx.ci property, as such:

   mvn -DskipTests=false -Dffx.ci=true

Currently, only JDK 1.8 is fully supported. Versions less than 1.8 do not support 1.8+ syntax such as lambda expressions. Versions above 1.8 are problematic due to changes in JDK policy, such as forbidding illegal reflective access operations. Purely experimental support for this has been added to the ffx/ffxc Bash scripts, though not the Windows batch scripts.

For most systems, simply install JDK 1.8, point the JAVA\_HOME environment variable to the JDK directory, and then add its bin directory to your path. For Mac OS X, it is recommended to add this to your .bash\_profile after installing JDK 1.8:

   export JAVA\_HOME=$(/usr/libexec/java\_home -v 1.8)

---

## Execute Force Field X

Once the Maven build succeeds, Force Field X can be executed using platform dependent start-up scripts located in the bin. On Mac OS X or Linux:

    bin/ffxc Energy examples/alamet.xyz

On Windows:

    bin/ffxc.bat Energy examples/alamet.xyz

The ffx/bin directory should be appended to your $PATH environment variable. The "Energy" command refers to an internal version of the Energy.groovy script that can be found in the ffx/scripts directory. To embed your own script within FFX, place it into the scripts directory and rebuild FFX.

---

