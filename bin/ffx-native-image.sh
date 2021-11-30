#! /bin/bash

native-image -H:+PrintClassInitialization -H:+ReportExceptionStackTraces -H:TraceClassInitialization=true --initialize-at-run-time=ffx.potential.bonded.RendererCache,jcuda.jcufft.JCufft,org.jogamp.java3d.SoundRetained,org.jogamp.java3d.Screen3D,org.jogamp.java3d.VirtualUniverse,org.jogamp.java3d.Canvas3D,edu.uiowa.jopenmm.OpenMMAmoebaLibrary,edu.uiowa.jopenmm.OpenMMLibrary,jcuda.driver.JCudaDriver --allow-incomplete-classpath --report-unsupported-elements-at-runtime --no-fallback --no-server -cp "bin/ffx-all-1.0.0-beta.jar" ffx.Main

