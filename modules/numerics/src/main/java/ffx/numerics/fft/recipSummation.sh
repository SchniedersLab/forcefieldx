#!/bin/sh

nvcc -m64 -gencode arch=compute_30,code=sm_30 -cubin recipSummation.cu -o recipSummation-64.cubin
