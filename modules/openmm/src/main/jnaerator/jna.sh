#!/bin/bash

java -jar jnaerator-0.12-shaded.jar config.jnaerator

tar -xvf openmm-7.2.0.jar

cp simtk/openmm/*java ../java/simtk/openmm/.

rm openmm-7.2.0.jar

rm -rf simtk
