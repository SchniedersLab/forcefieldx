//******************************************************************************
//
// File:    Configuration.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.Configuration
//
// This Java source file is copyright (C) 2012 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj.cluster;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Scanner;

/**
 * Class Configuration provides configuration information about a parallel
 * computer running Parallel Java. The configuration information is read from a
 * plain text file. Each configuration file entry is on a single line. Lines
 * beginning with <TT>#</TT> and blank lines are ignored. The order of the
 * entries in the file does not matter (unless stated otherwise below). The
 * items in each entry are separated by whitespace; there cannot be any
 * whitespace within an item (unless stated otherwise below). The configuration
 * file entries are:
 * <UL>
 * <P>
 * <LI>
 * <TT>cluster &lt;name&gt;</TT>
 * <BR>The name of the cluster is <TT>&lt;name&gt;</TT>. The name may contain
 * whitespace. This entry must be specified; there is no default.
 * <P>
 * <LI>
 * <TT>logfile &lt;file&gt;</TT>
 * <BR>The Job Scheduler will append log entries to the log file named
 * <TT>&lt;file&gt;</TT>. This entry must be specified; there is no default.
 * <P>
 * <LI>
 * <TT>webhost &lt;host&gt;</TT>
 * <BR>The host name for the Job Scheduler's web interface is
 * <TT>&lt;host&gt;</TT>. This entry must be specified; there is no default.
 * <P>
 * <LI>
 * <TT>webport &lt;port&gt;</TT>
 * <BR>The port number for the Job Scheduler's web interface is
 * <TT>&lt;port&gt;</TT>. If not specified, the default port number is 8080.
 * <P>
 * <LI>
 * <TT>schedulerhost &lt;host&gt;</TT>
 * <BR>The host name to which the Job Scheduler listens for connections from job
 * frontend processes is <TT>&lt;host&gt;</TT>. If not specified, the default is
 * <TT>"localhost"</TT>.
 * <P>
 * <LI>
 * <TT>schedulerport &lt;port&gt;</TT>
 * <BR>The port number to which the Job Scheduler listens for connections from
 * job frontend processes is <TT>&lt;port&gt;</TT>. If not specified, the
 * default port number is 20617.
 * <P>
 * <LI>
 * <TT>frontendhost &lt;host&gt;</TT>
 * <BR>The host name to which job frontend processes listen for connections from
 * job backend processes is <TT>&lt;host&gt;</TT>. This entry must be specified;
 * there is no default.
 * <P>
 * <LI>
 * <TT>backend &lt;name&gt; &lt;cpus&gt; &lt;host&gt; &lt;jvm&gt;
 * &lt;classpath&gt; [&lt;jvmflag&gt; ...]</TT>
 * <BR>The parallel computer includes a backend node named
 * <TT>&lt;name&gt;</TT> with <TT>&lt;cpus&gt;</TT> CPUs. The host name for SSH
 * remote logins to the backend node is <TT>&lt;host&gt;</TT>. The full pathname
 * for executing the Java Virtual Machine (JVM) on the backend node is
 * <TT>&lt;jvm&gt;</TT>. The Java class path for the Parallel Java Library on
 * the backend node is <TT>&lt;classpath&gt;</TT>. Each
 * <TT>&lt;jvmflag&gt;</TT> (zero or more) gives a flag passed to the JVM on the
 * command line. At least one of this entry must be specified.
 * <P>
 * <LI>
 * <TT>backendshell &lt;name&gt; &lt;shell command&gt;</TT>
 * <BR>On the backend node named <TT>&lt;name&gt;</TT>, use the given shell
 * command string when starting a job backend process. This entry, if present,
 * must appear after the corresponding <TT>backend &lt;name&gt;</TT> entry. If
 * this entry is omitted, the default shell command string is
 * <TT>"bash&nbsp;-l&nbsp;-c"</TT>.
 * <P>
 * <LI>
 * <TT>jobtime &lt;time&gt;</TT>
 * <BR>The maximum time in seconds any Parallel Java job is allowed to run. The
 * Job Scheduler will abort a job if it runs for this many seconds. If not
 * specified, the default is not to impose a maximum time on jobs. <I>Note:</I>
 * If the Job Scheduler is configured with a maximum job time and a particular
 * job is given a maximum time with the <TT>-Dpj.jobtime</TT> property, the
 * smaller of the Job Scheduler's maximum job time and the job's maximum time
 * will be used for that job.
 * </UL>
 * <P>
 * Here is an example of a configuration file:
 * <P>
 * <TABLE BORDER=1 CELLPADDING=10 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * <FONT SIZE="-1">
 * <PRE> # Parallel Java Job Scheduler configuration file
 * # Frontend node: tardis.cs.rit.edu
 * # Backend nodes: dr00-dr09
 *
 * cluster RIT CS Tardis Hybrid SMP Cluster
 * logfile /var/tmp/parajava/scheduler.log
 * webhost tardis.cs.rit.edu
 * webport 8080
 * schedulerhost localhost
 * schedulerport 20617
 * frontendhost 10.10.221.1
 * backend dr00 4 10.10.221.10 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr01 4 10.10.221.11 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr02 4 10.10.221.12 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr03 4 10.10.221.13 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr04 4 10.10.221.14 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr05 4 10.10.221.15 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr06 4 10.10.221.16 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr07 4 10.10.221.17 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr08 4 10.10.221.18 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar
 * backend dr09 4 10.10.221.19 /usr/local/versions/jdk-1.5.0_15/bin/java /var/tmp/parajava/pj.jar</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 *
 * @author Alan Kaminsky
 * @version 20-Jun-2012
 */
public class Configuration {

// Hidden data members.
    // Cluster name.
    private String myClusterName;

    // Log file.
    private String myLogFile;

    // Web interface host and port.
    private String myWebHost = Constants.ALL_NETWORK_INTERFACES;
    private int myWebPort = Constants.WEB_PORT;

    // Job Scheduler host and port.
    private String mySchedulerHost = "localhost";
    private int mySchedulerPort = Constants.PJ_PORT;

    // Frontend host.
    private String myFrontendHost;

    // List of backend information objects.
    private ArrayList<BackendInfo> myBackendInfo
            = new ArrayList<BackendInfo>();

    // Default shell comand string.
    private static final String DEFAULT_SHELL_COMMAND = "bash -l -c";

    // Maximum job time. 0 means no maximum.
    private int myJobTime;

// Exported constructors.
    /**
     * Construct a new configuration. The configuration information is read from
     * the given file.
     *
     * @param configfile Configuration file name.
     * @exception IOException Thrown if an I/O error occurred while reading the
     * configuration file. Thrown if there was an error in the configuration
     * file.
     * @throws java.io.IOException if any.
     */
    public Configuration(String configfile)
            throws IOException {
        parseConfigFile(configfile);
    }

// Exported operations.
    /**
     * Returns the cluster name.
     *
     * @return Cluster name.
     */
    public String getClusterName() {
        return myClusterName;
    }

    /**
     * Returns the Job Scheduler's log file name.
     *
     * @return Log file name.
     */
    public String getLogFile() {
        return myLogFile;
    }

    /**
     * Returns the Job Scheduler's web interface host name.
     *
     * @return Host name.
     */
    public String getWebHost() {
        return myWebHost;
    }

    /**
     * Returns the Job Scheduler's web interface port number.
     *
     * @return Port number.
     */
    public int getWebPort() {
        return myWebPort;
    }

    /**
     * Returns the Job Scheduler's channel group host name. To send messages to
     * the Job Scheduler, a job frontend connects a channel to this host.
     *
     * @return Host name.
     */
    public String getSchedulerHost() {
        return mySchedulerHost;
    }

    /**
     * Returns the Job Scheduler's channel group port number. To send messages
     * to the Job Scheduler, a job frontend connects a channel to this port.
     *
     * @return Port number.
     */
    public int getSchedulerPort() {
        return mySchedulerPort;
    }

    /**
     * Returns the host name of the cluster's frontend processor.
     *
     * @return Host name.
     */
    public String getFrontendHost() {
        return myFrontendHost;
    }

    /**
     * Returns the number of backend processors.
     *
     * @return Count.
     */
    public int getBackendCount() {
        return myBackendInfo.size();
    }

    /**
     * Returns information about the given backend processor.
     *
     * @param i Index in the range 0 .. <TT>getBackendCount()-1</TT>.
     * @return Backend information object.
     */
    public BackendInfo getBackendInfo(int i) {
        return myBackendInfo.get(i);
    }

    /**
     * Returns information about all backend processors.
     *
     * @return List of backend information objects.
     */
    public List<BackendInfo> getBackendInfoList() {
        return myBackendInfo;
    }

    /**
     * Returns the maximum job time.
     *
     * @return Maximum job time (seconds), or 0 if no maximum.
     */
    public int getJobTime() {
        return myJobTime;
    }

// Hidden operations.
    /**
     * Parse the configuration file.
     *
     * @param configfile Configuration file name.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private void parseConfigFile(String configfile)
            throws IOException {
        Scanner scanner = null;
        String line = null;
        long now = System.currentTimeMillis();
        try {
            scanner = new Scanner(new File(configfile));
            lineloop:
            while (scanner.hasNextLine()) {
                line = scanner.nextLine();
                Scanner linescanner = new Scanner(line);
                if (!linescanner.hasNext()) {
                    continue lineloop;
                }
                String command = linescanner.next();
                if (command.charAt(0) == '#') {
                } else if (command.equals("cluster")) {
                    myClusterName = linescanner.nextLine().trim();
                } else if (command.equals("logfile")) {
                    myLogFile = linescanner.next();
                } else if (command.equals("webhost")) {
                    myWebHost = linescanner.next();
                } else if (command.equals("webport")) {
                    myWebPort = Integer.parseInt(linescanner.next());
                } else if (command.equals("schedulerhost")) {
                    mySchedulerHost = linescanner.next();
                } else if (command.equals("schedulerport")) {
                    mySchedulerPort = Integer.parseInt(linescanner.next());
                } else if (command.equals("frontendhost")) {
                    myFrontendHost = linescanner.next();
                } else if (command.equals("backend")) {
                    String name = linescanner.next();
                    int cpus = linescanner.nextInt();
                    if (cpus < 1) {
                        throw new IOException("Invalid backend command, <cpus> must be >= 1: "
                                + line);
                    }
                    String host = linescanner.next();
                    String jvm = linescanner.next();
                    String classpath = linescanner.next();
                    ArrayList<String> jvmflags = new ArrayList<String>();
                    while (linescanner.hasNext()) {
                        jvmflags.add(linescanner.next());
                    }
                    BackendInfo backendinfo
                            = new BackendInfo(name,
                                    cpus,
                                    BackendInfo.State.IDLE,
                                    now,
                                    host,
                                    jvm,
                                    classpath,
                                    jvmflags.toArray(new String[jvmflags.size()]),
                                    DEFAULT_SHELL_COMMAND);
                    myBackendInfo.add(backendinfo);
                } else if (command.equals("backendshell")) {
                    String name = linescanner.next();
                    String shellCommand = linescanner.nextLine().trim();
                    BackendInfo backendinfo = backendInfoForName(name);
                    if (backendinfo == null) {
                        throw new IOException("Invalid backendshell command, no backend named \""
                                + name + "\"");
                    }
                    backendinfo.shellCommand = shellCommand;
                } else if (command.equals("jobtime")) {
                    int time = linescanner.nextInt();
                    if (time < 1) {
                        throw new IOException("Invalid configuration command: " + line);
                    }
                    myJobTime = time;
                } else {
                    throw new IOException("Invalid configuration command: " + line);
                }
            }
            if (myClusterName == null) {
                throw new IOException("Missing configuration command: cluster <name>");
            }
            if (myLogFile == null) {
                throw new IOException("Missing configuration command: logfile <file>");
            }
            if (myWebHost == null) {
                throw new IOException("Missing configuration command: webhost <host>");
            }
            if (myFrontendHost == null) {
                throw new IOException("Missing configuration command: frontendhost <host>");
            }
            if (myBackendInfo.isEmpty()) {
                throw new IOException("Missing configuration command: backend <name> <host> <port>");
            }
        } catch (NoSuchElementException exc) {
            throw new IOException("Invalid configuration command: " + line);
        } catch (NumberFormatException exc) {
            throw new IOException("Invalid configuration command: " + line);
        } finally {
            if (scanner != null) {
                scanner.close();
            }
        }
    }

// Hidden operations.
    /**
     * Returns the backend info object for the given backend name.
     *
     * @param name Backend name.
     *
     * @return Backend info, or null if <TT>name</TT> does not exist.
     */
    private BackendInfo backendInfoForName(String name) {
        for (BackendInfo backendinfo : myBackendInfo) {
            if (backendinfo.name.equals(name)) {
                return backendinfo;
            }
        }
        return null;
    }

}
