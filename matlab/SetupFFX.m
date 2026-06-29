% 
% JDK VERSION
%
% As of Matlab 2026b, both Matlab and Force Field X support JDK/JRE version 25.
% Please install a JDK 25 LTS version from Adopt Open JDK:
% https://adoptium.net/temurin/releases
%
% Matlab can be configured to use this JDK use the 'jenv' command:
jenv("/Library/Java/JavaVirtualMachines/jdk-25.0.3+9/Contents/Home")

% Force Java to be loaded
java.lang.String('Java has been loaded.')

% Check the current JDK version:
version -java
% ans =
%   'Java 25.0.3+9-LTS with Eclipse Adoptium OpenJDK 64-Bit Server VM mixed mode, emulated-client, sharing'


% ADDING FORCE FIELD X TO THE CLASSPATH
%
% Once Matlab is running on JRE/JDK 25, the FFX jar files need to be added
% to the "static" Matlab classpath.
FFX_HOME = "/Users/mjschnie/Data/ffx-project/forcefieldx";
ffxlib = FFX_HOME + "/lib/";
jarfiles = dir(ffxlib);

% Add JAR files to "javaclasspath.txt" in the Matlab pref directory. 
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

% Configure Java Start-Up Options to enable Vector Calculations
startupDir = matlabroot + "/bin/" + computer('arch');
cd(startupDir);
fid = fopen('java.opts','wt');
fprintf(fid, '%s\n', '--add-modules=jdk.incubator.vector');
fclose(fid);

% Uncomment the line below to check java.opts.
% edit java.opts

% Matlab must be restarted for the new static classpath and start-up options to take effect.