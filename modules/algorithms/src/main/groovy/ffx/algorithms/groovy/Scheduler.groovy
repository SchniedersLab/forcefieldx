package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import picocli.CommandLine.Command
import picocli.CommandLine.Option

/**
 * Run a Parallel Java Job Scheduler on a cluster with multiple nodes.
 *
 * <br>
 * Usage:
 * <br>
 * ffxc Scheduler [options] &lt;filename&gt;

 * @author Michael Schnieders
 */
@Command(description = " The Scheduler runs Parallel jobs over nodes.", name = "ffxc Scheduler")
class Scheduler extends AlgorithmsScript {


    /**
     * -v or --verbose to turn on verbose backend Parallel Java logging.
     */
    @Option(names = ['-v', '--verbose'],
            description = 'Turn on verbose backend Parallel Java logging.')
    boolean v = false

    /**
     * -p or --CPUs to define the processor cores (threads) per process (default is all available cores).
     */
    @Option(names = ['-p', '--CPUs'],
            description = 'The number of processor cores (threads) per process.')
    int p = -1

    /**
     * -P or --port to define the port the Server will use.
     */
    @Option(names = ['-P', '--port'], paramLabel = "20617",
            description = 'Set the port the Front End server will listen on.')
    int port = 20617

    /**
     * -e or --hostfile to define the environment variable that points to the host file (default is PE_HOSTFILE).
     */
    @Option(names = ['-e', '--hostfile'], paramLabel = 'PE_HOSTFILE',
            description = 'Environment variable that points to the host file.')
    String hostfileName = "PE_HOSTFILE"

    /**
     * -m or --memory to define the string value of -Xmx to pass to worker nodes (default is '2G').
     */
    @Option(names = ['-m', '--memory'], paramLabel = '2G',
            description = 'String value of -Xmx to pass to worker nodes.')
    String memory = "2G"

    @Override
    Scheduler run() {

        if (!init()) {
            return
        }

        // Determine the number of CPUs per node
        int CPUs = Runtime.getRuntime().availableProcessors()

        // The default is 1 process per node.
        int processes = 1
        // More than 1 process per node can be selected if no cores are wasted.
        if (p > 0) {
            int n = p
            if (n < CPUs && CPUs % n == 0) {
                processes = CPUs / n
                CPUs = n
            }
        }

        // Read in the Parallel Environment Host File environment variable.
        String hostsFile = System.getenv(hostfileName)

        // Default to using only the current node (i.e. "localhost")
        String[] hostnames = new String[1]
        hostnames[0] = "localhost"

        if (hostsFile == null) {
            logger.info(" The " + hostfileName + " environment variable is empty.")
            logger.info(" Only localhost will be used.\n")
        } else {
            // Check that the supplied file exists and is readable.
            File host = new File(hostsFile)
            if (
            !host.exists() || !host.canRead()) {
                logger.info(" The file path specified by the " + hostfileName
                        + " environment variable does not exist or cannot be read.")
                logger.info(" Only localhost will be used.\n")
            } else {
                // Read in the hosts.
                List nodes = host.readLines()
                hostnames = new String[nodes.size()]
                int i = 0
                for (line in nodes) {
                    hostnames[i] = line.split(" +")[0]
                    i++
                }
            }
        }

        // Create the Parallel Java cluster configuration file.
        String frontend = hostnames[0]
        StringBuffer sb = new StringBuffer()
        sb.append("# Force Field X Cluster Configuration File\n")
        sb.append("cluster Force Field X TACC Cluster\n")
        sb.append("logfile ffx-scheduler.log\n")
        sb.append("webhost 127.0.0.1\n")
        sb.append("webport 8080\n")
        sb.append("schedulerhost localhost\n")
        sb.append("schedulerport " + port + "\n")
        sb.append("frontendhost " + frontend + "\n")

        // Locate the JRE being used.
        String javaHome = System.getProperty("java.home")
        // Locate the version of FFX being used.
        String ffxHome = System.getProperty("basedir")

        String java = javaHome + "/bin/java"
        String ffx = ffxHome + "/bin/ffx-all.jar"
        args = "-Xmx" + memory + " -Djava.system.class.loader='ffx.FFXClassLoader'"

        if (v) {
            args = args + " -Dpj.verbose='true'"
        }

        // Create an entry for each process
        int i = 0
        for (p = 0; p < processes; p++) {
            for (node in hostnames) {
                sb.append("backend node" + i + " "
                        + CPUs + " "
                        + node + " "
                        + java + " "
                        + ffx + " "
                        + args + "\n")
                i++
            }
        }

        // Write the Parallel Java config file.
        String pjConfig = "cluster.txt"
        File config = new File(pjConfig)
        config.write(sb.toString())

        // Run the Parallel Java Scheduler.
        String command = ffxHome + "/bin/scheduler " + pjConfig
        def process = command.execute()
        process.waitFor()

        return this
    }

    @Override
    public List<Potential> getPotentials() {
        return new ArrayList<>();
    }
}
