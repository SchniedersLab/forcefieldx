//******************************************************************************
//
// File:    PJProperties.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.PJProperties
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj;

/**
 * Class PJProperties provides static methods for reading Java system properties
 * that control the behavior of Parallel Java. The PJ properties are:
 * <UL>
 * <LI>
 * <B>pj.nn</B> -- The number of backend nodes in a parallel program.
 * ({@link #getPjNn()})
 * <LI>
 * <B>pj.np</B> -- The number of processes in a parallel program; also the size
 * of the world communicator. ({@link #getPjNp()})
 * <LI>
 * <B>pj.nt</B> -- The number of CPUs per process in a parallel program; also
 * the default number of threads in a parallel team. ({@link #getPjNt()})
 * <p>
 * When running a Parallel Java program via the job queue on a cluster parallel
 * computer (see package {@linkplain edu.rit.pj.cluster edu.rit.pj.cluster} for
 * further information), the Job Scheduler Daemon uses the <I>nn,</I>
 * <I>np,</I> and <I>nt</I> settings to assign resources to the job.
 * <p>
 * If neither <I>nn</I> nor <I>np</I> is specified, the Job Scheduler will run
 * the job with one process on one node.
 * <p>
 * If <I>nn</I> or <I>np</I> is specified but not both, the Job Scheduler will
 * run the job on <I>nn</I> (or <I>np</I>) nodes with one process on each node.
 * <p>
 * If <I>nn</I> and <I>np</I> are both specified and <I>nn</I> &gt;= <I>np</I>,
 * the Job Scheduler will run the job on <I>np</I> nodes with one process on
 * each node.
 * <p>
 * If <I>nn</I> and <I>np</I> are both specified and <I>nn</I> &lt; <I>np</I>,
 * the Job Scheduler will run the job on <I>nn</I> nodes with with more than one
 * process on some or all of the nodes, apportioning the <I>np</I> processes as
 * equally as possible among the <I>nn</I> nodes. Note that in this case,
 * different nodes may be assigned different numbers of processes.
 * <p>
 * On each node, the Job Scheduler will assign <I>nt</I> CPUs to each process.
 * If <I>nt</I> is not specified, the default is to use all the CPUs in the
 * node, apportioning the CPUs as equally as possible among the processes on the
 * node. Note that in this case, different processes may be assigned different
 * numbers of CPUs.
 * <LI>
 * <B>pj.schedule</B> -- The schedule for a parallel loop in a parallel program.
 * ({@link #getPjSchedule()})
 * <LI>
 * <B>pj.host</B> -- The host name of the Job Scheduler Daemon to use when
 * running a cluster parallel program. ({@link #getPjHost()})
 * <LI>
 * <B>pj.port</B> -- The port number of the Job Scheduler Daemon to use when
 * running a cluster parallel program. ({@link #getPjPort()})
 * <LI>
 * <B>pj.jobtime</B> -- The maximum amount of time (seconds) the job is allowed
 * to run when running a cluster parallel program. ({@link #getPjJobTime()})
 * <LI>
 * <B>pj.jvmflags</B> -- JVM flags to include on the Java command line when
 * running a backend process in a cluster parallel program. ({@link
 * #getPjJvmFlags()})
 * <LI>
 * <B>pj.prng</B> -- The fully-qualified class name of the default pseudorandom
 * number generator (PRNG) class. ({@link #getPjPrng()})
 * </UL>
 * <p>
 * You can specify a PJ property on the Java command line like this:
 * <p>
 * <code>&nbsp;&nbsp;&nbsp;&nbsp;java -Dpj.nt=4 . . .</code>
 *
 * @author Alan Kaminsky
 * @version 21-May-2008
 */
public class PJProperties {

    // Prevent construction.
    private PJProperties() {
    }

// Exported operations.

    /**
     * Determine the number of backend nodes in a parallel program.
     * <p>
     * If the <code>"pj.nn"</code> Java system property is specified, it must be an
     * integer greater than or equal to 1.
     *
     * @return Number of backend nodes for a parallel program.
     * @throws IllegalArgumentException (unchecked exception) Thrown if the
     *                                  <code>"pj.nn"</code> property value is not an integer greater than or equal
     *                                  to 1.
     */
    public static int getPjNn() {
        int k = 1;
        String pj_nn = System.getProperty("pj.nn");
        if (pj_nn != null) {
            try {
                k = Integer.parseInt(pj_nn);
            } catch (NumberFormatException exc) {
                throw new IllegalArgumentException("pj.nn system property is not an integer >= 1");
            }
            if (k < 1) {
                throw new IllegalArgumentException("pj.nn system property is not an integer >= 1");
            }
        } else // (pj_nn == null)
        {
            String pj_np = System.getProperty("pj.np");
            if (pj_np != null) {
                try {
                    k = Integer.parseInt(pj_np);
                } catch (NumberFormatException exc) {
                    throw new IllegalArgumentException("pj.np system property is not an integer >= 1");
                }
                if (k < 1) {
                    throw new IllegalArgumentException("pj.np system property is not an integer >= 1");
                }
            }
        }
        return k;
    }

    /**
     * Determine the number of processes in a parallel program. This is the
     * number of backend processes set up when the <code>Comm.init()</code> method
     * is executed (see class {@linkplain Comm}). This is also the size of the
     * world communicator.
     * <p>
     * If the <code>"pj.np"</code> Java system property is specified, it must be an
     * integer greater than or equal to 1.
     *
     * @return Number of processes for a parallel program.
     * @throws IllegalArgumentException (unchecked exception) Thrown if the
     *                                  <code>"pj.np"</code> property value is not an integer greater than or equal
     *                                  to 1.
     */
    public static int getPjNp() {
        int k = 1;
        String pj_np = System.getProperty("pj.np");
        if (pj_np != null) {
            try {
                k = Integer.parseInt(pj_np);
            } catch (NumberFormatException exc) {
                throw new IllegalArgumentException("pj.np system property is not an integer >= 1");
            }
            if (k < 1) {
                throw new IllegalArgumentException("pj.np system property is not an integer >= 1");
            }
        } else // (pj_np == null)
        {
            String pj_nn = System.getProperty("pj.nn");
            if (pj_nn != null) {
                try {
                    k = Integer.parseInt(pj_nn);
                } catch (NumberFormatException exc) {
                    throw new IllegalArgumentException("pj.nn system property is not an integer >= 1");
                }
                if (k < 1) {
                    throw new IllegalArgumentException("pj.nn system property is not an integer >= 1");
                }
            }
        }
        return k;
    }

    /**
     * Determine the number of CPUs per process in a parallel program. This is
     * the number of threads a {@linkplain ParallelTeam} will have if the number
     * of threads is not specified as a constructor argument.
     * <p>
     * If the <code>"pj.nt"</code> Java system property is specified, it must be an
     * integer greater than or equal to 1.
     * <p>
     * If the <code>"pj.nt"</code> Java system property is not specified, this
     * method returns 0 to signify that all available CPUs should be used.
     *
     * @return Number of CPUs per process for a parallel program, or 0 if not
     * specified.
     * @throws IllegalArgumentException (unchecked exception) Thrown if the
     *                                  <code>"pj.nt"</code> property value is not an integer greater than or equal
     *                                  to 1.
     */
    public static int getPjNt() {
        int k = 0;
        String pj_nt = System.getProperty("pj.nt");
        if (pj_nt != null) {
            try {
                k = Integer.parseInt(pj_nt);
            } catch (NumberFormatException exc) {
                throw new IllegalArgumentException("pj.nt system property is not an integer >= 1");
            }
            if (k < 1) {
                throw new IllegalArgumentException("pj.nt system property is not an integer >= 1");
            } else {
                System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", pj_nt);
            }
        }
        return k;
    }

    /**
     * Determine the schedule for a parallel loop in an SMP parallel program.
     * This is the schedule that will be used if the parallel for loop has a
     * runtime schedule. For further information, see class {@linkplain
     * IntegerSchedule} and class {@linkplain LongSchedule}.
     * <p>
     * If the <code>"pj.schedule"</code> Java property is specified, that value
     * gives the type of schedule, which must be one of the following:
     * <UL>
     * <LI><code>"fixed"</code> -- Fixed schedule.
     * <LI><code>"dynamic"</code> -- Dynamic schedule with a chunk size of 1.
     * <LI><code>"dynamic(&lt;n&gt;)"</code> -- Dynamic schedule with a chunk size
     * of <code>&lt;n&gt;</code>, an integer &gt;= 1.
     * <LI><code>"guided"</code> -- Self-guided schedule with a minimum chunk size
     * of 1.
     * <LI><code>"guided(&lt;n&gt;)"</code> -- Self-guided schedule with a minimum
     * chunk size of <code>&lt;n&gt;</code>, an integer &gt;= 1.
     * <LI><code>"<I>classname</I>"</code> -- Schedule that is an instance of the
     * given class. <I>classname</I> is the fully-qualified class name of the
     * schedule class. The instance is constructed using the subclass's
     * no-argument constructor.
     * <LI><code>"<I>classname</I>(<I>arg</I>,<I>arg</I>,...)"</code> -- Schedule
     * that is an instance of the given class. <I>classname</I> is the
     * fully-qualified class name of the schedule class. The arguments between
     * the parentheses are split into separate strings separated by commas.
     * There cannot be parentheses or commas within the arguments themselves.
     * The instance is constructed using the subclass's constructor whose
     * argument is an array of Strings, namely the individual arguments between
     * the parentheses.
     * </UL>
     * <p>
     * If the <code>"pj.schedule"</code> Java property is not specified, the default
     * schedule for a runtime schedule is used. Normally this is a fixed
     * schedule, but a program can specify a different default.
     *
     * @return Schedule for a parallel for loop (one of the above strings), or
     * null if the <code>"pj.schedule"</code> Java property is not specified.
     */
    public static String getPjSchedule() {
        return System.getProperty("pj.schedule");
    }

    /**
     * Determine the host name of the Job Scheduler Daemon to use when running a
     * cluster parallel program. The program contacts the Job Scheduler Daemon
     * when the <code>Comm.init()</code> method is executed (see class {@linkplain
     * Comm}). For further information, see package {@linkplain
     * edu.rit.pj.cluster edu.rit.pj.cluster} and class {@linkplain
     * edu.rit.pj.cluster.JobScheduler}.
     * <p>
     * If the <code>"pj.host"</code> Java system property is specified, it gives the
     * Job Scheduler Daemon's host name (or IP address).
     * <p>
     * If the <code>"pj.host"</code> Java system property is not specified, a host
     * name of <code>"localhost"</code> is returned.
     *
     * @return Job Scheduler Daemon host name.
     */
    public static String getPjHost() {
        return System.getProperty("pj.host", "localhost");
    }

    /**
     * Determine the port number of the Job Scheduler Daemon to use when running
     * a cluster parallel program. The program contacts the Job Scheduler Daemon
     * when the <code>Comm.init()</code> method is executed (see class {@linkplain
     * Comm}). For further information, see package {@linkplain
     * edu.rit.pj.cluster edu.rit.pj.cluster} and class {@linkplain
     * edu.rit.pj.cluster.JobScheduler}.
     * <p>
     * If the <code>"pj.port"</code> Java system property is specified, it gives the
     * Job Scheduler Daemon's port number.
     * <p>
     * If the <code>"pj.port"</code> Java system property is not specified, the
     * well-known Parallel Java port number (20617) is returned.
     *
     * @return Job Scheduler Daemon port number.
     * @throws IllegalArgumentException (unchecked exception) Thrown if the
     *                                  <code>"pj.port"</code> property value is not an integer.
     */
    public static int getPjPort() {
        try {
            String pj_port = System.getProperty("pj.port");
            int port
                    = pj_port == null
                    ? // edu.rit.pj.cluster.Constants.PJ_PORT :
                    20617
                    : Integer.parseInt(pj_port);
            return port;
        } catch (NumberFormatException exc) {
            throw new IllegalArgumentException("pj.port system property is not an integer");
        }
    }

    /**
     * Determine the maximum amount of time (seconds) the job is allowed to run
     * when running a cluster parallel program. If a program running in the
     * Parallel Java job queue has not finished in the given amount of time, the
     * job frontend process automatically termintes the job when the given
     * number of seconds have elapsed since the job started executing. If a
     * program is not running in the Parallel Java job queue (if there is no Job
     * Scheduler Daemon present), then the job time setting is ignored and the
     * program will not time out.
     * <p>
     * If the <code>"pj.jobtime"</code> Java system property is specified, it must
     * be an integer greater than or equal to 1, and that gives the the maximum
     * job timeout in seconds.
     * <p>
     * If the <code>"pj.jobtime"</code> Java system property is not specified, a
     * value of 0 is returned to signify that there is no job timeout.
     *
     * @return Job timeout (seconds), or 0 if no job timeout.
     * @throws IllegalArgumentException (unchecked exception) Thrown if the
     *                                  <code>"pj.jobtime"</code> property value is not an integer greater than or
     *                                  equal to 1.
     */
    public static int getPjJobTime() {
        try {
            String pj_jobtime = System.getProperty("pj.jobtime");
            if (pj_jobtime == null) {
                return 0;
            } else {
                int jobtime = Integer.parseInt(pj_jobtime);
                if (jobtime < 1) {
                    throw new IllegalArgumentException("pj.jobtime system property is not an integer >= 1");
                }
                return jobtime;
            }
        } catch (NumberFormatException exc) {
            throw new IllegalArgumentException("pj.jobtime system property is not an integer >= 1");
        }
    }

    /**
     * Determine the JVM flags to include on the Java command line when running
     * a backend process in a cluster parallel program. When a job backend
     * process is started, the JVM flags are included on the command line
     * immediately after the <code>"java"</code> command. These flags then control
     * the job backend process's JVM. For further information, see package
     * {@linkplain edu.rit.pj.cluster edu.rit.pj.cluster}.
     * <p>
     * If the <code>"pj.jvmflags"</code> Java system property is specified, it gives
     * the JVM flags exactly as they are to appear on the Java command line. If
     * there are multiple flags separated by whitespace, the
     * <code>"pj.jvmflags"</code> Java system property must be enclosed in quotation
     * marks.
     * <p>
     * If the <code>"pj.jvmflags"</code> Java system property is not specified,
     * there are no JVM flags, and an empty string is returned.
     * <p>
     * <B>Example.</B> To cause the job backend processes' JVMs to use an
     * initial heap size of 4 MB and a maximum heap size of 128 MB, specify the
     * <code>"pj.jvmflags"</code> property as follows when running the program:
     *
     * <code>&nbsp;&nbsp;&nbsp;&nbsp;java -Dpj.jvmflags="-Xms4m -Xmx128m" . .
     * .</code>
     * <p>
     * Note that quotation marks are needed around the property value because of
     * the embedded whitespace. This property value causes the Job Launcher
     * Daemon to launch each job backend process's JVM with this command:
     *
     * <code>&nbsp;&nbsp;&nbsp;&nbsp;java -Xms4m -Xmx128m . . .</code>
     * <p>
     * which in turn tells each job backend process's JVM to use initial and
     * maximum heap sizes of 4 MB and 128 MB.
     *
     * @return JVM flags.
     */
    public static String getPjJvmFlags() {
        return System.getProperty("pj.jvmflags", "");
    }

    /**
     * Determine the fully-qualified class name of the default pseudorandom
     * number generator (PRNG) class.
     * <p>
     * If the <code>"pj.prng"</code> Java system property is specified, it gives the
     * fully-qualified class name of the PRNG class that the static
     * <code>getInstance(long)</code> method in class {@linkplain
     * edu.rit.util.Random} will construct. Specifying the <code>"pj.prng"</code>
     * property will substitute a different PRNG algorithm into a program
     * without needing to recompile. See class {@linkplain edu.rit.util.Random}
     * for further information.
     * <p>
     * If the <code>"pj.prng"</code> Java system property is not specified, a PRNG
     * class name of <code>"edu.rit.util.DefaultRandom"</code> is returned. See
     * class {@linkplain edu.rit.util.DefaultRandom} for further information.
     *
     * @return Default PRNG class name.
     */
    public static String getPjPrng() {
        return System.getProperty("pj.prng", "edu.rit.util.DefaultRandom");
    }

}
