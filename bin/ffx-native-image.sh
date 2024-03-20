#! /bin/bash

# export NATIVE_FLAGS="-H:+PrintClassInitialization -H:+ReportExceptionStackTraces -H:TraceClassInitialization=true --report-unsupported-elements-at-runtime"
export NATIVE_FLAGS="--report-unsupported-elements-at-runtime"

# native-image $NATIVE_FLAGS --report-unsupported-elements-at-runtime --no-fallback -cp ".:lib/*" ffx.numerics.fft.Complex3DParallel Complex3D

native-image $NATIVE_FLAGS --no-fallback -cp ".:lib/*" ffx.potential.commands.Energy Energy 

