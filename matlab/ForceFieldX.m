%
% This Force Field X example should be run after completing the
% setup process described in the SetupFFX sript.
%
% Matlab Java syntax is available here:
% https://www.mathworks.com/help/matlab/using-java-libraries-in-matlab.html
%
% Both JDK 10 and FFX classes can now be imported.

import java.util.ArrayList
import groovy.lang.Binding
import ffx.potential.groovy.Energy

% To avoid absolute paths, its possible (but optional) to set an FFX_HOME 
% environment variable to the root FFX directoy.

FFX_HOME = getenv("FFX_HOME");

% Script input Arguments
A = ArrayList;
A.add("-m");
A.add(FFX_HOME + "/examples/watersmall.xyz");

% Create a Groovy Script Binding
binding = Binding;
binding.setVariable("args",A);

% Run the Energy Groovy Script
script = Energy;
script.setBinding(binding);
script.run();

% Get back the ForceFieldEnergy instance
forceFieldEnergy = script.forceFieldEnergy;
forceFieldEnergy.toString()
