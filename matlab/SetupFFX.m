% 
% JDK VERSION
%
% As of Matlab 2024b, both Matlab and Force Field X support JDK/JRE version 21.
% Please install a JDK 21 LTS version from Adopt Open JDK:
% https://adoptium.net/temurin/releases
%
% Matlab can be configured to use this JDK use the 'jenv' command:
jenv("/Library/Java/JavaVirtualMachines/jdk-21.0.5+11/Contents/Home")

% Restart Matlab and then check the current JDK version:
version -java
% ans =
%   'Java 21.0.5+11-LTS with Eclipse Adoptium OpenJDK 64-Bit Server VM mixed mode, emulated-client, sharing'

%
% ADDING FORCE FIELD X TO THE CLASSPATH
%
% Once Matlab is running on JRE/JDK 21, the FFX jar files need to be added
% to the "static" Matlab classpath.

FFX_HOME = "/Users/mjschnie/Data/ffx-project/forcefieldx";
ffxlib = FFX_HOME + "/lib/"
jarfiles = dir(ffxlib)

% Add JAR files "javaclasspath.txt" in the Matlab pref directory. 
cd(prefdir)
fid = fopen('javaclasspath.txt','wt');
for n = 1:length(jarfiles)
    file = jarfiles(n);
    if (length(file.name) > 4)
        jar = ffxlib + file.name;
        fprintf(fid, '%s\n', jar);
%        javaaddpath(jar)
    end
end
fclose(fid);

% Uncomment the line below to check javaclasspath.txt.
% edit javaclasspath.txt

% Matlab must be restarted for the new static classpath to take effect.
