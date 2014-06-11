/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */

// Run a Parallel Java Job Scheduler on a cluster with multiple nodes.

// Groovy Imports
import groovy.util.CliBuilder;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Determine the number of CPUs per node
int CPUs = Runtime.getRuntime().availableProcessors();
String cpuString = CPUs.toString();
String memory = "1G";

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc scheduler dummy ');
cli.e(longOpt:'hostfile', args:1, argName:'PE_HOSTFILE', 'Environment variable that points to the host file.');
cli.p(longOpt:'cpus', args:1, argName:cpuString, 'The number of processor cores (threads) per process.');
cli.h(longOpt:'help', 'Print this help message.');
cli.m(longOpt:'memory', args:1, argName:'1G', 'String value of -Xmx to pass to worker nodes.');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h) {
    return cli.usage();
}

// The default is 1 process per node.
int processes = 1;
// More than 1 process per node can be selected if no cores are wasted.
if (options.p) {
    int n = Integer.parseInt(options.p);
    if (n < CPUs && CPUs % n == 0) {
        processes = CPUs / n;
        CPUs = n;
    }
}

if (options.m) {
    memory = options.m;
}

// Allow specification of an alternative to PE_HOSTFILE.
String hostfileName = "PE_HOSTFILE";
if (options.e) {
    hostfileName = options.e;
}
// Read in the Parallel Environment Host File environment variable.
String hostsFile = System.getenv(hostfileName);
if (hostsFile == null) {
    logger.severe("The " + hostfileName + " environment variable is empty.");
}

// Check that the supplied file exists and is readable.
File host = new File(hostsFile);
if (!host.exists() || !host.canRead()) {
    logger.severe("The file path specified by the " + hostfileName
        + " environment variable does not exist or cannot be read.");
}

// Read in the hosts.
List nodes = host.readLines();
String[] hostnames = new String[nodes.size()];
int i=0;
for (line in nodes) {
    hostnames[i] = line.split(" +")[0];
    i++;
}

// Create the Parallel Java cluster configuration file.
String frontend = hostnames[0];
StringBuffer sb = new StringBuffer();
sb.append("# Force Field X Cluster Configuration File\n");
sb.append("cluster Force Field X TACC Cluster\n");
sb.append("logfile ffx-scheduler.log\n");
sb.append("webhost 127.0.0.1\n");
sb.append("webport 8080\n");
sb.append("schedulerhost localhost\n");
sb.append("schedulerport 20617\n");
sb.append("frontendhost " + frontend + "\n");

// Locate the JRE being used.
String javaHome = System.getProperty("java.home");
// Locate the version of FFX being used.
String ffxHome = System.getProperty("basedir");

java = javaHome + "/bin/java";
ffx = ffxHome + "/bin/ffx-all.jar";
args = "-Xmx" + memory + " -Djava.system.class.loader='ffx.FFXClassLoader'";

// Create an entry for each process
i = 0;
for (p = 0; p < processes; p++) {
    for (node in hostnames) {
        sb.append("backend node" + i + " "
            + CPUs + " "
            + node + " "
            + java + " "
            + ffx + " "
            + args + "\n");
        i++;
    }
}

// Write the Parallel Java config file.
String pjConfig = "cluster.txt";
File config = new File(pjConfig);
config.write(sb.toString());

// Run the Parallel Java scheduler.
String command = ffxHome + "/bin/scheduler " + pjConfig;
def process = command.execute();
process.waitFor();

