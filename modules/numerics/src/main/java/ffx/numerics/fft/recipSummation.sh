#!/bin/sh

nvcc -m32 -arch sm_11 -cubin recipSummation.cu -o recipSummation-32.cubin
nvcc -m64 -arch sm_11 -cubin recipSummation.cu -o recipSummation-64.cubin
