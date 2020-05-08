//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.utilities.PortUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option

import static java.lang.String.format

/**
 * Run a Parallel Java Job Scheduler on a cluster with multiple nodes.
 *
 * <br>
 * Usage:
 * <br>
 * ffxc Scheduler [options] &lt;filename&gt;

 * @author Michael Schnieders
 */
@Command(description = " The Scheduler runs parallel jobs over nodes.", name = "ffxc Scheduler")
class Scheduler extends AlgorithmsScript {

  /**
   * -v or --verbose to turn on verbose backend Parallel Java logging.
   */
  @Option(names = ['-v', '--verbose'],
      description = 'Turn on verbose back end logging.')
  boolean v = false

  /**
   * -p or --pj.nt to define the processor cores (threads) per process (default is all available cores).
   */
  @Option(names = ['-p', '--pj.nt'], paramLabel = 'all',
      description = 'The number of processor cores (threads) per process.')
  int p = -1

  /**
   * -P or --port to define the port the Server will use.
   */
  @Option(names = ['-P', '--port'], paramLabel = '20617',
      description = 'Set the port the front end server will listen on.')
  int port = 20617

  /**
   * -W or --webPort to define the port the Server will serve a webpage to (generally not used).
   */
  @Option(names = ['-W', '--webPort'], paramLabel = '8080',
      description = 'Set the port the server will serve a webpage to.')
  int webPort = 8080

  /**
   * --ib or --infiniband Replace the "hpc" domain with the "ipoib" domain to use the Argon high-speed network.
   */
  @Option(names = ['--ib', '--ipoib'],
      description = 'Replace the "hpc" domain with the "ipoib" domain to use the Argon high-speed network.')
  boolean ipoib = false

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
      return null
    }

    // Determine the number of CPUs per node
    int CPUs = Runtime.getRuntime().availableProcessors()

    // The default is 1 process per node.
    int processes = 1
    // More than 1 process per node can be selected if no cores are wasted.
    if (p > 0) {
      int n = p
      if (n < CPUs && CPUs % n == 0) {
        processes = (int) (CPUs / n)
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
          if (ipoib) {
            hostnames[i] = hostnames[i].replace("hpc", "ipoib")
          }
          i++
        }
      }
    }

    // Check for an invalid requested port.
    if (port <= 1) {
      logger.info(" The scheduler port must be greater than 1; the default of 20617 will be used.")
      port = 20617
    }

    // Check the availability of the desired Scheduler port.
    while (!PortUtils.isTcpPortAvailable(port)) {
      logger.info(format(" Scheduler port %d is not available.", port))
      if (++port > PortUtils.MAX_TCP_PORT) {
        logger.severe(" Reached port 65535 without finding an open scheduler port!")
      }
    }
    logger.info(format(" Scheduler port: %d.", port))
    String logFile = format("scheduler.%d.log", port)

    // Check the availability of the desired web port.
    while (!PortUtils.isTcpPortAvailable(webPort) || webPort == port) {
      logger.info(format(" Web port %d is not available.", webPort))
      if (++webPort > PortUtils.MAX_TCP_PORT) {
        logger.severe(" Reached port 65535 without finding an open web port!")
      }
    }
    logger.info(format(" Web port: %d.", webPort))

    // Create the Parallel Java cluster configuration file.
    String frontend = hostnames[0]
    StringBuffer sb = new StringBuffer()
    sb.append("# Force Field X Cluster Configuration File\n")
    sb.append("cluster Force Field X Cluster\n")
    sb.append(format("logfile %s\n", logFile))
    sb.append("webhost 127.0.0.1\n")
    sb.append(format("webport %d\n", webPort))
    sb.append("schedulerhost localhost\n")
    sb.append(format("schedulerport %d\n", port))
    sb.append(format("frontendhost %s\n", frontend))

    // Locate the JRE being used.
    String javaHome = System.getProperty("java.home")
    // Locate the version of FFX being used.
    String ffxHome = System.getProperty("basedir")

    String java = javaHome + "/bin/java"
    String ffx = ffxHome + "/bin/ffx-all-1.0.0-beta.jar"
    args = "-Xmx" + memory

    if (v) {
      args = args + " -Dpj.log='true'"
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

}
