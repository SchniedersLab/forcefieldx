
package ffx.algorithms

import groovy.cli.Option
import groovy.cli.Unparsed

/**
 * Run a Parallel Java Job Scheduler on a cluster with multiple nodes.
 *
 * <br>
 * Usage:
 * <br>
 * ffxc Scheduler [options] &lt;filename&gt;

 * @author Michael Schnieders
 */
class Scheduler extends Script {

    /**
     * Options for the Scheduler Script.
     * <br>
     * Usage:
     * <br>
     * ffxc Scheduler [options]
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName='h', longName='help', defaultValue='false',
                description='Print this help message.') boolean help;
        /**
         * -p or --CPUs to define the processor cores (threads) per process (default is all available cores).
         */
        @Option(shortName='p', longName='CPUs', 
                description='The number of processor cores (threads) per process.') int p;

        /**
         * -e or --hostfile to define the environment variable that points to the host file (default is PE_HOSTFILE).
         */
        @Option(shortName='e', longName='hostfile', defaultValue='PE_HOSTFILE',
                description='Environment variable that points to the host file.') String e;

        /**
         * -m or --memory to define the string value of -Xmx to pass to worker nodes (default is '2G').
         */
        @Option(shortName='m', longName='memory', defaultValue='2G',
                description='String value of -Xmx to pass to worker nodes.') String m;

        /**
         * There should be no final argument(s) to the Scheduler.
         */
        @Unparsed List<String> filenames;
    }

    def run() {
        def cli = new CliBuilder(usage: ' ffxc Scheduler [options]', header: ' Options:');

        AlgorithmFunctions afuncts;
        try {
            afuncts = getAlgorithmUtils();
        } catch (MissingMethodException ex) {
            afuncts = new AlgorithmUtils();
        }

        def options = new Options();
        cli.parseFromInstance(options, afuncts.getArguments());

        if (options.help == true) {
            return cli.usage();
        }

        // Determine the number of CPUs per node
        int CPUs = Runtime.getRuntime().availableProcessors();
        String memory = "2G";

        // The default is 1 process per node.
        int processes = 1;
        // More than 1 process per node can be selected if no cores are wasted.
        if (options.p) {
            int n = options.p;
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

        // Default to using only the current node (i.e. "localhost")
        String[] hostnames = new String[1];
        hostnames[0] = "localhost";

        if (hostsFile == null) {
            logger.info(" The " + hostfileName + " environment variable is empty.");
            logger.info(" Only localhost will be used.\n");
        } else {
            // Check that the supplied file exists and is readable.
            File host = new File(hostsFile);
            if (
            !host.exists() || !host.canRead()) {
                logger.info(" The file path specified by the " + hostfileName
                        + " environment variable does not exist or cannot be read.");
                logger.info(" Only localhost will be used.\n");
            } else {
                // Read in the hosts.
                List nodes = host.readLines();
                hostnames = new String[nodes.size()];
                int i = 0;
                for (line in nodes) {
                    hostnames[i] = line.split(" +")[0];
                    i++;
                }
            }
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

        // Run the Parallel Java Scheduler.
        String command = ffxHome + "/bin/scheduler " + pjConfig;
        def process = command.execute();
        process.waitFor();
    }

}
