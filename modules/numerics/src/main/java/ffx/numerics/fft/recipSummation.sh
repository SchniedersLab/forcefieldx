#!/bin/sh

nvcc -m64 -arch sm_13 -cubin recipSummation.cu -o recipSummation-64.cubin
