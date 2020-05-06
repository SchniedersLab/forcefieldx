% 
% JDK VERSION
%
% As of Matlab 2020a, Force Field X support depends on using JDK/JRE version 10.
% Please install JDK 10 from Oracle or from Adopt Open JDK:
% https://www.oracle.com/java/technologies/java-archive-javase10-downloads.html
% - or -
% https://adoptopenjdk.net/releases.html?variant=openjdk10
%
% Please then set the MATLAB_JAVA environment variable
% For csh / tcsh
% setenv MATLAB_JAVA <JRE path>
% For bash
% export MATLAB_JAVA=<JRE path>
% 
% For example, using bash on MacOS with for Oracle JDK 10:
% export MATLAB_JAVA=/Library/Java/JavaVirtualMachines/jdk-10.0.2.jdk/Contents/Home
%
% The restart Matlab from the command line (the GUI shortcut will not pick up
% the environment variable). 
% For example, on MacOS
% /Applications/MATLAB_R2020a.app/bin/matlab
%
% To check the current Matlab JRE, use the command

version -java

% ans =
% 'Java 10.0.2+13 with "Oracle Corporation" Java HotSpot(TM) 64-Bit Server VM mixed mode'
%
% ADDING FORCE FIELD X TO THE CLASSPATH
%
% Once Matlab is running on JRE/JDK 10, the ffx-all.jar needs to be added
% to the "static" Matlab classpath.
% https://www.mathworks.com/help/matlab/matlab_external/static-path.html
%
% To avoid absolute paths, its possible (but optional) to set an FFX_HOME 
% environment variable to the root FFX directoy.

FFX_HOME = getenv("FFX_HOME");

% Append /bin/ffx-all.jar to FFX_HOME

FFX_JAR = FFX_HOME + "/bin/ffx-all-1.0.0-beta.jar"

% CD into the preferences directory

cd(prefdir)

% Add ffx-all.jar to the file "javaclasspath.txt". The value of the 
% FFX_JAR variable was printed to the command window. 

edit javaclasspath.txt

% Matlab must be restarted for the new static classpath to take effect.