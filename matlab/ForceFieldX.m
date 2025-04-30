%
% This Force Field X example should be run after completing the
% setup process described in the SetupFFX sript.
%
% Matlab Java syntax is available here:
% https://www.mathworks.com/help/matlab/using-java-libraries-in-matlab.html

import java.util.ArrayList
import groovy.lang.Binding
import ffx.potential.groovy.Energy

% Set the FFX_HOME variable to the root FFX directoy.

FFX_HOME = "/Users/mjschnie/Data/ffx-project/forcefieldx";

% Script input Arguments
A = ArrayList;
A.add("-m");
A.add(FFX_HOME + "/examples/peptide.xyz");

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
