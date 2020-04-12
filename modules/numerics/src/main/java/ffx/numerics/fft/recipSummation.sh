#!/bin/sh

nvcc -ptx -Xptxas="-v" recipSummation.cu -o recipSummation.ptx
