//******************************************************************************
//
// File:    JobFrontend.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobFrontend
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
package edu.rit.pj.cluster;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.InetSocketAddress;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;

import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;
import edu.rit.mp.ChannelGroupClosedException;
import edu.rit.mp.ObjectBuf;
import edu.rit.mp.Status;
import edu.rit.mp.buf.ObjectItemBuf;
import edu.rit.pj.PJProperties;
import edu.rit.util.ByteSequence;
import edu.rit.util.Timer;
import edu.rit.util.TimerTask;
import edu.rit.util.TimerThread;

/**
 * Class JobFrontend provides the message handler for the PJ job frontend
 * process.
 *
 * @author Alan Kaminsky
 * @version 20-Jun-2012
 */
public class JobFrontend
        implements Runnable, JobFrontendRef {

// Hidden data members.
    // User name.
    private String username;

    // Job number.
    private int jobnum;

    // Job resources.
    private int Nn;
    private int Np;
    private int Nt;

    // Whether the frontend communicator exists, true or false.
    private boolean hasFrontendComm;

    // Main class name.
    private String myMainClassName;

    // Command line arguments.
    private String[] myArgs;

    // Rank of next backend process to be assigned.
    private int myNextRank;

    // Timer thread for lease renewals and expirations.
    private TimerThread myLeaseTimerThread;

    // Timers for the lease with the Job Scheduler.
    private Timer mySchedulerRenewTimer;
    private Timer mySchedulerExpireTimer;

    // Timer for the job timeout if any.
    private Timer myJobTimer;

    // Array of job backend process info records, indexed by rank.
    private ProcessInfo[] myProcessInfo;

    // Mapping from job backend reference to job backend process info record.
    private Map<JobBackendRef, ProcessInfo> myProcessMap
            = new HashMap<JobBackendRef, ProcessInfo>();

    // Number of running job backend processes.
    private int myRunningCount;

    // Number of finished job backend processes.
    private int myFinishedCount;

    // Middleware channel group and address array.
    private ChannelGroup myMiddlewareChannelGroup;
    private InetSocketAddress[] myMiddlewareAddress;

    // Proxy for Job Scheduler Daemon.
    private JobSchedulerRef myJobScheduler;

    // World communicator channel group address array.
    private InetSocketAddress[] myWorldAddress;

    // Frontend communicator channel group and address array.
    private ChannelGroup myFrontendChannelGroup;
    private InetSocketAddress[] myFrontendAddress;

    // JVM flags.
    private String userJvmFlags = PJProperties.getPjJvmFlags();

    // Resource contents that have been reported to job backend processes.
    private ResourceCache myResourceCache = new ResourceCache();

    // Flag for shutting down the run() method.
    private boolean continueRun = true;

    // State of this job frontend.
    private State myState = State.RUNNING;

    private static enum State {

        RUNNING,
        TERMINATE_CANCEL_JOB,
        TERMINATING
    };

    // Error message if job canceled, or null if job finished normally.
    private String myCancelMessage = "User canceled job";

    // For writing and reading files on the job frontend's node.
    private FrontendFileWriter myFrontendFileWriter;
    private FrontendFileReader myFrontendFileReader;

// Exported constructors.
    /**
     * Construct a new job frontend object. The job frontend object will contact
     * the Job Scheduler Daemon specified by the <TT>"pj.host"</TT> and
     * <TT>"pj.port"</TT> Java system properties. See class {@linkplain
     * edu.rit.pj.PJProperties} for further information.
     *
     * @param username User name.
     * @param Nn Number of backend nodes (&gt;= 1).
     * @param Np Number of processes (&gt;= 1).
     * @param Nt Number of CPUs per process (&gt;= 0). 0 means "all CPUs."
     * @param hasFrontendComm True if the job has the frontend communicator,
     * false if it doesn't.
     * @param mainClassName Main class name.
     * @param args Command line arguments.
     * @exception JobSchedulerException (subclass of IOException) Thrown if the
     * job frontend object could not contact the Job Scheduler Daemon.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public JobFrontend(String username,
            int Nn,
            int Np,
            int Nt,
            boolean hasFrontendComm,
            String mainClassName,
            String[] args)
            throws IOException {
        // Record arguments.
        this.username = username;
        this.Nn = Nn;
        this.Np = Np;
        this.Nt = Nt;
        this.hasFrontendComm = hasFrontendComm;
        this.myMainClassName = mainClassName;
        this.myArgs = args;

        // Set up shutdown hook.
        Runtime.getRuntime().addShutdownHook(new Thread() {
            public void run() {
                shutdown();
            }
        });

        // Set up lease timer thread.
        myLeaseTimerThread = new TimerThread();
        myLeaseTimerThread.setDaemon(true);
        myLeaseTimerThread.start();

        // Set up Job Scheduler lease timers.
        mySchedulerRenewTimer
                = myLeaseTimerThread.createTimer(new TimerTask() {
                    public void action(Timer timer) {
                        try {
                            schedulerRenewTimeout();
                        } catch (Throwable exc) {
                        }
                    }
                });
        mySchedulerExpireTimer
                = myLeaseTimerThread.createTimer(new TimerTask() {
                    public void action(Timer timer) {
                        try {
                            schedulerExpireTimeout();
                        } catch (Throwable exc) {
                        }
                    }
                });

        // Set up job timer.
        myJobTimer
                = myLeaseTimerThread.createTimer(new TimerTask() {
                    public void action(Timer timer) {
                        try {
                            jobTimeout();
                        } catch (Throwable exc) {
                        }
                    }
                });

        // Set up array of job backend process info records.
        myProcessInfo = new ProcessInfo[Np];
        for (int i = 0; i < Np; ++i) {
            final int rank = i;
            ProcessInfo processinfo
                    = new ProcessInfo(/*state            */ProcessInfo.State.NOT_STARTED,
                            /*name             */ null,
                            /*rank             */ rank,
                            /*backend          */ null,
                            /*middlewareAddress*/ null,
                            /*worldAddress     */ null,
                            /*frontendAddress  */ null,
                            /*renewTimer       */
                            myLeaseTimerThread.createTimer(new TimerTask() {
                                public void action(Timer timer) {
                                    try {
                                        backendRenewTimeout(rank);
                                    } catch (Throwable exc) {
                                    }
                                }
                            }),
                            /*expireTimer      */
                            myLeaseTimerThread.createTimer(new TimerTask() {
                                public void action(Timer timer) {
                                    try {
                                        backendExpireTimeout(rank);
                                    } catch (Throwable exc) {
                                    }
                                }
                            }),
                            /*Nt               */ 0);
            myProcessInfo[rank] = processinfo;
        }

        // Set up middleware channel group and address array.
        myMiddlewareChannelGroup = new ChannelGroup();
        myMiddlewareAddress = new InetSocketAddress[Np + 1];

        // Set up world communicator address array.
        myWorldAddress = new InetSocketAddress[Np];

        // Set up frontend communicator channel group and address array.
        if (hasFrontendComm) {
            myFrontendChannelGroup = new ChannelGroup();
            myFrontendAddress = new InetSocketAddress[Np + 1];
        }

        // Set up frontend file writer and reader.
        myFrontendFileWriter = new FrontendFileWriter(this);
        myFrontendFileReader = new FrontendFileReader(this);

        // Set up Job Scheduler proxy.
        InetSocketAddress js_address = null;
        Channel js_channel = null;
        try {
            js_address
                    = new InetSocketAddress(PJProperties.getPjHost(),
                            PJProperties.getPjPort());
            js_channel = myMiddlewareChannelGroup.connect(js_address);
        } catch (IOException exc) {
            throw new JobSchedulerException("JobFrontend(): Cannot contact Job Scheduler Daemon at "
                    + js_address,
                    exc);
        }
        myJobScheduler
                = new JobSchedulerProxy(myMiddlewareChannelGroup, js_channel);

        // Start Job Scheduler lease timers.
        mySchedulerRenewTimer.start(Constants.LEASE_RENEW_INTERVAL,
                Constants.LEASE_RENEW_INTERVAL);
        mySchedulerExpireTimer.start(Constants.LEASE_EXPIRE_INTERVAL);

        // Kick off the job!
        myJobScheduler.requestJob(this, username, Nn, Np, Nt);
    }

// Exported operations.
    /**
     * Run this Job Frontend.
     */
    public void run() {
        ObjectItemBuf<JobFrontendMessage> buf
                = ObjectBuf.buffer((JobFrontendMessage) null);
        Status status = null;
        JobFrontendMessage message = null;
        JobBackendRef backend = null;

        try {
            while (continueRun) {
                // Receive a message from any channel.
                status = myMiddlewareChannelGroup.receive(null, null, buf);
                message = buf.item;

                // Process a message from the Job Scheduler.
                if (status.tag == Message.FROM_JOB_SCHEDULER) {
                    message.invoke(this, myJobScheduler);
                } // Process a message from a job backend.
                else if (status.tag == Message.FROM_JOB_BACKEND) {
                    // Get job backend associated with channel. If none, set up
                    // a new job backend proxy.
                    backend = (JobBackendRef) status.channel.info();
                    if (backend == null) {
                        backend
                                = new JobBackendProxy(myMiddlewareChannelGroup, status.channel);
                        status.channel.info(backend);
                    }

                    // Process message.
                    message.invoke(this, backend);
                }

                // Enable garbage collection of no-longer-needed objects while
                // waiting to receive next message.
                buf.item = null;
                status = null;
                message = null;
                backend = null;
            }
        } catch (ChannelGroupClosedException exc) {
        } catch (Throwable exc) {
            terminateCancelJob(exc);
        }

        // Exit process if necessary.
        switch (myState) {
            case TERMINATE_CANCEL_JOB:
                System.exit(1);
                break;
            case RUNNING:
            case TERMINATING:
                break;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Assign a backend process to the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void assignBackend(JobSchedulerRef theJobScheduler,
            String name,
            String host,
            String jvm,
            String classpath,
            String[] jvmflags,
            String shellCommand,
            int Nt)
            throws IOException {
        // Record backend name and number of CPUs.
        int rank = myNextRank++;
        ProcessInfo processinfo = myProcessInfo[rank];
        processinfo.name = name;
        processinfo.Nt = Nt;

        // Display backend.
        System.err.print(", ");
        System.err.print(name);
        System.err.flush();
        if (myNextRank == Np) {
            System.err.println();
        }

        try {
            // Build a command to run on the backend node.
            StringBuilder command = new StringBuilder();
            command.append(shellCommand);
            command.append(" \"");
            String cwd = System.getProperty("user.dir");
            if (cwd != null) {
                command.append("cd '");
                command.append(cwd);
                command.append("'; ");
            }
            command.append("nohup ");
            command.append(jvm);
            command.append(" -classpath '");
            command.append(classpath);
            command.append("'");
            for (String flag : jvmflags) {
                command.append(" ");
                command.append(flag);
            }
            command.append(" ");
            command.append(userJvmFlags);
            command.append(" edu.rit.pj.cluster.JobBackend '");
            command.append(username);
            command.append("' ");
            command.append(jobnum);
            command.append(" ");
            command.append(Np);
            command.append(" ");
            command.append(rank);
            command.append(" ");
            command.append(hasFrontendComm);
            command.append(" '");
            command.append(myMiddlewareChannelGroup.listenAddress().getHostName());
            command.append("' ");
            command.append(myMiddlewareChannelGroup.listenAddress().getPort());
            command.append(" '");
            command.append(host);

            String pjLogging = System.getProperty("pj.log");
            if (pjLogging != null && Boolean.parseBoolean(pjLogging)) {
                command.append("' &\"");
            } else {
                command.append("' >/dev/null 2>/dev/null &\"");
            }

            // So an SSH remote login and execute the above command.
            Process ssh
                    = Runtime.getRuntime().exec(new String[]{"ssh", host, command.toString()});

            // Start lease timers for the backend node.
            processinfo.renewTimer.start(Constants.LEASE_RENEW_INTERVAL,
                    Constants.LEASE_RENEW_INTERVAL);
            processinfo.expireTimer.start(Constants.LEASE_EXPIRE_INTERVAL);
        } // If an I/O error occurs, treat it as a backend node failure.
        catch (IOException exc) {
            if (myNextRank != Np) {
                System.err.println();
            }
            System.err.println(" Exception executing SSH command:\n" + exc.toString());
            terminateCancelJob(backendFailed(processinfo));
        }
    }

    /**
     * {@inheritDoc}
     *
     * Assign a job number to the job. The host name for the job frontend's
     * middleware channel group is also specified.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void assignJobNumber(JobSchedulerRef theJobScheduler,
            int jobnum,
            String pjhost)
            throws IOException {
        // Record job number.
        this.jobnum = jobnum;

        // Start listening for connections to the middleware channel group.
        myMiddlewareChannelGroup.listen(new InetSocketAddress(pjhost, 0));
        myMiddlewareChannelGroup.startListening();
        myMiddlewareAddress[Np] = myMiddlewareChannelGroup.listenAddress();

        // Start listening for connections to the frontend communicator channel
        // group.
        if (hasFrontendComm) {
            myFrontendChannelGroup.listen(new InetSocketAddress(pjhost, 0));
            myFrontendChannelGroup.startListening();
            myFrontendAddress[Np] = myFrontendChannelGroup.listenAddress();
        }

        // Report job number.
        System.err.print("Job " + jobnum);
        System.err.flush();
    }

    /**
     * {@inheritDoc}
     *
     * Cancel the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void cancelJob(JobSchedulerRef theJobScheduler,
            String errmsg)
            throws IOException {
        terminateCancelJob(errmsg);
    }

    /**
     * Renew the lease on the job.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public synchronized void renewLease(JobSchedulerRef theJobScheduler)
            throws IOException {
        mySchedulerExpireTimer.start(Constants.LEASE_EXPIRE_INTERVAL);
    }

    /**
     * {@inheritDoc}
     *
     * Report that a backend process has finished executing the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void backendFinished(JobBackendRef theJobBackend)
            throws IOException {
        ProcessInfo processinfo = myProcessMap.get(theJobBackend);
        if (processinfo == null) {
            return;
        }

        // Verify that this backend has not finished already.
        if (processinfo.state != ProcessInfo.State.RUNNING) {
            terminateCancelJob("Unexpected \"backend finished\" message, rank="
                    + processinfo.rank);
        }

        // Update job backend process state.
        processinfo.state = ProcessInfo.State.FINISHED;

        // Increase count of finished processes.
        ++myFinishedCount;

        // If all job backend processes have finished, terminate the run()
        // method. This will cause the job frontend process to exit when all
        // other non-daemon threads have also terminated.
        if (myFinishedCount == Np) {
            continueRun = false;
            myCancelMessage = null;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Report that a backend process is ready to commence executing the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void backendReady(JobBackendRef theJobBackend,
            int rank,
            InetSocketAddress middlewareAddress,
            InetSocketAddress worldAddress,
            InetSocketAddress frontendAddress)
            throws IOException {
        // Verify that rank is in range.
        if (0 > rank || rank >= Np) {
            terminateCancelJob("Illegal \"backend ready\" message, rank=" + rank);
        }

        // Verify that this backend has not started already.
        ProcessInfo processinfo = myProcessInfo[rank];
        if (processinfo.state != ProcessInfo.State.NOT_STARTED) {
            terminateCancelJob("Unexpected \"backend ready\" message, rank=" + rank);
        }

        // Record information in job backend process info record.
        processinfo.state = ProcessInfo.State.RUNNING;
        processinfo.backend = theJobBackend;
        processinfo.middlewareAddress = middlewareAddress;
        processinfo.worldAddress = worldAddress;
        processinfo.frontendAddress = frontendAddress;
        myProcessMap.put(theJobBackend, processinfo);

        // Record channel group addresses.
        myMiddlewareAddress[rank] = middlewareAddress;
        myWorldAddress[rank] = worldAddress;
        if (hasFrontendComm) {
            myFrontendAddress[rank] = frontendAddress;
        }

        // Increase count of running processes.
        ++myRunningCount;

        // If all job backend processes have reported ready, commence job.
        if (myRunningCount == Np) {
            // Start job timer if necessary.
            int jobtime = PJProperties.getPjJobTime();
            if (jobtime > 0) {
                myJobTimer.start(jobtime * 1000L);
            }

            // Get the system properties.
            Properties props = System.getProperties();

            // Send "commence job" message to each job backend, with system
            // property "pj.nt" set to the proper number of CPUs.
            for (ProcessInfo info : myProcessMap.values()) {
                props.setProperty("pj.nt", "" + info.Nt);
                info.backend.commenceJob(/*theJobFrontend   */this,
                        /*middlewareAddress*/ myMiddlewareAddress,
                        /*worldAddress     */ myWorldAddress,
                        /*frontendAddress  */ myFrontendAddress,
                        /*properties       */ props,
                        /*mainClassName    */ myMainClassName,
                        /*args             */ myArgs);
            }
        }
    }

    /**
     * {@inheritDoc}
     *
     * Cancel the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void cancelJob(JobBackendRef theJobBackend,
            String errmsg)
            throws IOException {
        terminateCancelJob(errmsg);
    }

    /**
     * {@inheritDoc}
     *
     * Renew the lease on the job.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void renewLease(JobBackendRef theJobBackend)
            throws IOException {
        ProcessInfo processinfo = myProcessMap.get(theJobBackend);
        if (processinfo != null) {
            processinfo.expireTimer.start(Constants.LEASE_EXPIRE_INTERVAL);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Request the given resource from this job frontend's class loader.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void requestResource(JobBackendRef theJobBackend,
            String resourceName)
            throws IOException {
        // To hold resource content.
        byte[] content = null;

        // Get resource content. If resource not found, content is null.
        if (myResourceCache.contains(resourceName)) {
            // Get resource content from cache.
            content = myResourceCache.getNoWait(resourceName);
        } else {
            // Get resource content from class loader, save it in cache.
            InputStream stream
                    = getClass().getClassLoader().getResourceAsStream(resourceName);
            if (stream != null) {
                content = new ByteSequence(stream).toByteArray();
            }
            myResourceCache.put(resourceName, content);
        }

        // Send resource to job backend.
        theJobBackend.reportResource(this, resourceName, content);
    }

    /**
     * {@inheritDoc}
     *
     * Open the given output file for writing or appending.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void outputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file,
            boolean append)
            throws IOException {
        myFrontendFileWriter.outputFileOpen(theJobBackend, bfd, file, append);
    }

    /**
     * {@inheritDoc}
     *
     * Write the given bytes to the given output file. <TT>ffd</TT> = 1 refers
     * to the job's standard output stream; <TT>ffd</TT> = 2 refers to the job's
     * standard error stream; other values refer to a previously opened file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void outputFileWrite(JobBackendRef theJobBackend,
            int ffd,
            byte[] buf,
            int off,
            int len)
            throws IOException {
        myFrontendFileWriter.outputFileWrite(theJobBackend, ffd, len);
    }

    /**
     * {@inheritDoc}
     *
     * Flush accumulated bytes to the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void outputFileFlush(JobBackendRef theJobBackend,
            int ffd)
            throws IOException {
        myFrontendFileWriter.outputFileFlush(theJobBackend, ffd);
    }

    /**
     * {@inheritDoc}
     *
     * Close the given output file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void outputFileClose(JobBackendRef theJobBackend,
            int ffd)
            throws IOException {
        myFrontendFileWriter.outputFileClose(theJobBackend, ffd);
    }

    /**
     * {@inheritDoc}
     *
     * Open the given input file for reading.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void inputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file)
            throws IOException {
        myFrontendFileReader.inputFileOpen(theJobBackend, bfd, file);
    }

    /**
     * {@inheritDoc}
     *
     * Read bytes from the given input file. <TT>ffd</TT> = 1 refers to the
     * job's standard input stream; other values refer to a previously opened
     * file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void inputFileRead(JobBackendRef theJobBackend,
            int ffd,
            int len)
            throws IOException {
        myFrontendFileReader.inputFileRead(theJobBackend, ffd, len);
    }

    /**
     * {@inheritDoc}
     *
     * Skip bytes from the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void inputFileSkip(JobBackendRef theJobBackend,
            int ffd,
            long len)
            throws IOException {
        myFrontendFileReader.inputFileSkip(theJobBackend, ffd, len);
    }

    /**
     * {@inheritDoc}
     *
     * Close the given input file.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void inputFileClose(JobBackendRef theJobBackend,
            int ffd)
            throws IOException {
        myFrontendFileReader.inputFileClose(theJobBackend, ffd);
    }

    /**
     * {@inheritDoc}
     *
     * Report a comment for a process.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public synchronized void reportComment(JobBackendRef theJobBackend,
            int rank,
            String comment)
            throws IOException {
        myJobScheduler.reportComment(this, rank, comment);
    }

    /**
     * Close communication with this Job Frontend.
     */
    public void close() {
    }

// Hidden operations.
    /**
     * Take action when the Job Scheduler's lease renewal timer times out.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private synchronized void schedulerRenewTimeout()
            throws IOException {
        if (mySchedulerRenewTimer.isTriggered()) {
            myJobScheduler.renewLease(this);
        }
    }

    /**
     * Take action when the Job Scheduler's lease expiration timer times out.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private void schedulerExpireTimeout()
            throws IOException {
        boolean doExit = false;
        synchronized (this) {
            if (mySchedulerExpireTimer.isTriggered()) {
                continueRun = false;
                if (myState == State.RUNNING) {
                    myState = State.TERMINATE_CANCEL_JOB;
                    myCancelMessage = "Job Scheduler failed";
                    System.err.println(myCancelMessage);
                    doExit = true;
                }
            }
        }

        // Cannot hold the synchronization lock while calling System.exit(),
        // otherwise a deadlock can occur between this thread (the timer thread)
        // and the shutdown hook thread.
        if (doExit) {
            System.exit(1);
        }
    }

    /**
     * Take action when the job timer times out.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private void jobTimeout()
            throws IOException {
        boolean doExit = false;
        synchronized (this) {
            if (myJobTimer.isTriggered()) {
                continueRun = false;
                if (myState == State.RUNNING) {
                    myState = State.TERMINATE_CANCEL_JOB;
                    myCancelMessage = "Job exceeded maximum running time";
                    System.err.println(myCancelMessage);
                    doExit = true;
                }
            }
        }

        // Cannot hold the synchronization lock while calling System.exit(),
        // otherwise a deadlock can occur between this thread (the timer thread)
        // and the shutdown hook thread.
        if (doExit) {
            System.exit(1);
        }
    }

    /**
     * Take action when a job backend process's lease renewal timer times out.
     *
     * @param rank Job backend process's rank.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private synchronized void backendRenewTimeout(int rank)
            throws IOException {
        ProcessInfo processinfo = myProcessInfo[rank];
        if (processinfo.renewTimer.isTriggered()) {
            processinfo.backend.renewLease(this);
        }
    }

    /**
     * Take action when a job backend process's lease expiration timer times
     * out.
     *
     * @param rank Job backend process's rank.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private void backendExpireTimeout(int rank)
            throws IOException {
        boolean doExit = false;
        synchronized (this) {
            ProcessInfo processinfo = myProcessInfo[rank];
            if (processinfo.expireTimer.isTriggered()) {
                // Terminate the Job Frontend.
                String msg = backendFailed(processinfo);
                continueRun = false;
                if (myState == State.RUNNING) {
                    myState = State.TERMINATE_CANCEL_JOB;
                    myCancelMessage = msg;
                    System.err.println(myCancelMessage);
                    doExit = true;
                }
            }
        }

        // Cannot hold the synchronization lock while calling System.exit(),
        // otherwise a deadlock can occur between this thread (the timer thread)
        // and the shutdown hook thread.
        if (doExit) {
            System.exit(1);
        }
    }

    /**
     * Take action when a backend process fails.
     *
     * @param processinfo Process info.
     *
     * @return Error message.
     */
    private String backendFailed(ProcessInfo processinfo) {
        // Mark the backend process as failed.
        processinfo.state = ProcessInfo.State.FAILED;

        // Tell the Job Scheduler that the backend process failed.
        try {
            myJobScheduler.backendFailed(this, processinfo.name);
        } catch (IOException exc) {
        }

        // Set up error message.
        return "Job backend process failed, node "
                + processinfo.name + ", rank " + processinfo.rank;
    }

    /**
     * Terminate this Job Frontend immediately, sending a "cancel job" message
     * to the Job Scheduler and all Job Backends. The error message is
     * <TT>msg</TT>. This method must only be called by the thread calling
     * <TT>run()</TT>.
     *
     * @param msg Error message.
     */
    private void terminateCancelJob(String msg) {
        continueRun = false;
        if (myState == State.RUNNING) {
            myState = State.TERMINATE_CANCEL_JOB;
            myCancelMessage = msg;
            System.err.println(myCancelMessage);
        }
    }

    /**
     * Terminate this Job Frontend immediately, sending a "cancel job" message
     * to the Job Scheduler and all Job Backends. The error message comes from
     * the given exception. This method must only be called by the thread
     * calling <TT>run()</TT>.
     *
     * @param exc Exception.
     */
    private void terminateCancelJob(Throwable exc) {
        continueRun = false;
        if (myState == State.RUNNING) {
            myCancelMessage = exc.getClass().getName();
            String msg = exc.getMessage();
            if (msg != null) {
                myCancelMessage = myCancelMessage + ": " + msg;
            }
            System.err.println(myCancelMessage);
            exc.printStackTrace(System.err);
        }
    }

    /**
     * Terminate this Job Frontend immediately, sending a "cancel job" message
     * to the Job Scheduler and all Job Backends. The error message comes from
     * the given exception. This method must only be called by a thread other
     * than the thread calling <TT>run()</TT>.
     *
     * @param exc Exception.
     */
    void terminateCancelJobOther(Throwable exc) {
        boolean doExit = false;
        synchronized (this) {
            continueRun = false;
            if (myState == State.RUNNING) {
                myCancelMessage = exc.getClass().getName();
                String msg = exc.getMessage();
                if (msg != null) {
                    myCancelMessage = myCancelMessage + ": " + msg;
                }
                System.err.println(myCancelMessage);
                exc.printStackTrace(System.err);
                doExit = true;
            }
        }

        // Cannot hold the synchronization lock while calling System.exit(),
        // otherwise a deadlock can occur between this thread and the shutdown
        // hook thread.
        if (doExit) {
            System.exit(1);
        }
    }

    /**
     * Shut down this Job Frontend.
     */
    private void shutdown() {
        synchronized (this) {
            // Stop all lease timers.
            mySchedulerRenewTimer.stop();
            mySchedulerExpireTimer.stop();
            for (ProcessInfo processinfo : myProcessInfo) {
                processinfo.renewTimer.stop();
                processinfo.expireTimer.stop();
            }

            // If state is RUNNING but myCancelMessage is not null, it means the
            // user canceled the job (e.g., by hitting CTRL-C).
            if (myState == State.RUNNING && myCancelMessage != null) {
                myState = State.TERMINATE_CANCEL_JOB;
            }

            // Inform Job Scheduler and Job Backends.
            switch (myState) {
                case RUNNING:
                    // Send "job finished" messages.
                    for (ProcessInfo processinfo : myProcessInfo) {
                        if (processinfo.backend != null) {
                            try {
                                processinfo.backend.jobFinished(this);
                            } catch (IOException exc) {
                            }
                        }
                    }
                    if (myJobScheduler != null) {
                        try {
                            myJobScheduler.jobFinished(this);
                        } catch (IOException exc) {
                        }
                    }
                    break;
                case TERMINATE_CANCEL_JOB:
                    // Send "cancel job" messages.
                    for (ProcessInfo processinfo : myProcessInfo) {
                        if (processinfo.backend != null
                                && processinfo.state != ProcessInfo.State.FAILED) {
                            try {
                                processinfo.backend.cancelJob(this, myCancelMessage);
                            } catch (IOException exc) {
                            }
                        }
                    }
                    if (myJobScheduler != null) {
                        try {
                            myJobScheduler.cancelJob(this, myCancelMessage);
                        } catch (IOException exc) {
                        }
                    }
                    break;
                case TERMINATING:
                    // Send nothing.
                    break;
            }

            // Record that we are terminating.
            myState = State.TERMINATING;
        }

        // All proxies, channels, and channel groups will close when the process
        // exits.
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.pj.cluster.JobFrontend <I>username</I> <I>K</I>
//	 * <I>hasFrontendComm</I> <I>mainClassName</I> [ <I>arg</I> . . . ]
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		if (args.length < 4) usage();
//		String username = args[0];
//		int K = Integer.parseInt (args[1]);
//		boolean hasFrontendComm = Boolean.parseBoolean (args[2]);
//		String mainClassName = args[3];
//		int n = args.length - 4;
//		String[] cmdargs = new String [n];
//		System.arraycopy (args, 4, cmdargs, 0, n);
//
//		new JobFrontend (username, K, hasFrontendComm, mainClassName, cmdargs)
//					.run();
//		}
//
//	/**
//	 * Print a usage message and exit.
//	 */
//	private static void usage()
//		{
//		System.err.println ("Usage: java edu.rit.pj.cluster.JobFrontend <username> <K> <hasFrontendComm> <mainClassName> [<arg>...]");
//		System.exit (1);
//		}
}
