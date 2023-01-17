//******************************************************************************
//
// File:    JobBackend.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobBackend
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

import static java.lang.String.format;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.CountDownLatch;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import edu.rit.mp.ChannelGroup;
import edu.rit.mp.ChannelGroupClosedException;
import edu.rit.mp.ObjectBuf;
import edu.rit.mp.Status;
import edu.rit.mp.buf.ObjectItemBuf;
import edu.rit.util.ByteSequence;
import edu.rit.util.Timer;
import edu.rit.util.TimerThread;

/**
 * Class JobBackend is the main program for a job backend process in the PJ
 * cluster middleware. The job backend process is launched by an SSH remote
 * login from the job frontend process (class {@linkplain JobFrontend}).
 * <p>
 * The command line for the job backend main program is:
 * <p>
 * java edu.rit.pj.cluster.JobBackend <I>username</I> <I>jobnum</I> <I>K</I>
 * <I>rank</I> <I>hasFrontendComm</I> <I>frontendHost</I> <I>frontendPort</I>
 * <I>backendHost</I>
 * <BR><I>username</I> = User name
 * <BR><I>jobnum</I> = Job number
 * <BR><I>K</I> = Number of backend processes (&gt;= 1)
 * <BR><I>rank</I> = Rank of this backend process (0 .. <I>K</I>-1)
 * <BR><I>hasFrontendComm</I> = Whether the frontend communicator exists
 * (<code>true</code> or <code>false</code>)
 * <BR><I>frontendHost</I> = Job frontend's middleware channel group host name
 * <BR><I>frontendPort</I> = Job frontend's middleware channel group port number
 * <BR><I>backendHost</I> = Job backend's middleware channel group host name
 *
 * @author Alan Kaminsky
 * @version 24-Jan-2012
 */
public class JobBackend
        implements Runnable, JobBackendRef {

    // Logger
    private static final Logger logger = Logger.getLogger(JobBackend.class.getName());

    // Hidden class-wide data members.
    private static JobBackend theJobBackend;

    // Hidden data members.
    // Command line arguments.
    private final String username;
    private final int jobnum;
    private final int K;
    private final int rank;
    private final boolean hasFrontendComm;
    private final String frontendHost;
    private final int frontendPort;
    private final String backendHost;

    // Logging
    FileHandler fileHandler = null;

    // Timer thread for lease renewals and expirations.
    private final TimerThread myLeaseTimerThread;

    // Timers for the lease with the job frontend.
    private final Timer myFrontendRenewTimer;
    private final Timer myFrontendExpireTimer;

    // Middleware channel group and address array.
    private final ChannelGroup myMiddlewareChannelGroup;
    private InetSocketAddress[] myMiddlewareAddress;

    // Job frontend proxy.
    private final JobFrontendRef myJobFrontend;

    // For loading classes from the job frontend process.
    private final ResourceCache myResourceCache;
    private final BackendClassLoader myClassLoader;

    // World communicator channel group and address array.
    private final ChannelGroup myWorldChannelGroup;
    private InetSocketAddress[] myWorldAddress;

    // Frontend communicator channel group and address array.
    private ChannelGroup myFrontendChannelGroup;
    private InetSocketAddress[] myFrontendAddress;

    // Java system properties.
    private Properties myProperties;

    // Main class name.
    private String myMainClassName;

    // Command line arguments.
    private String[] myArgs;

    // Flag set true to commence job.
    private boolean commence;

    // Buffer for receiving job backend messages.
    private final ObjectItemBuf<JobBackendMessage> myBuffer
            = ObjectBuf.buffer((JobBackendMessage) null);

    // Flags for shutting down the run() method.
    private boolean continueRun = true;
    private final CountDownLatch runFinished = new CountDownLatch(1);

    // State of this job backend.
    private State myState = State.RUNNING;

    private static enum State {
        RUNNING,
        TERMINATE_CANCEL_JOB,
        TERMINATE_NO_REPORT,
        TERMINATING
    }

    // Error message if job canceled, or null if job finished normally.
    private String myCancelMessage;

    // Original standard error stream; goes to the Job Launcher's log file.
    private final PrintStream myJobLauncherLog;

    // For writing and reading files in the job frontend.
    private final BackendFileWriter myFileWriter;
    private final BackendFileReader myFileReader;

// Hidden constructors.

    /**
     * Construct a new Job Backend.
     *
     * @param username        User name.
     * @param jobnum          Job number.
     * @param K               Number of backend processes.
     * @param rank            Rank of this backend process.
     * @param hasFrontendComm Whether the frontend communicator exists.
     * @param frontendHost    Host name of job frontend's middleware channel group.
     * @param frontendPort    Port number of job frontend's middleware channel group.
     * @param backendHost     Host name of job backend's middleware channel group.
     * @throws IOException Thrown if an I/O error occurred.
     */
    private JobBackend(String username,
                       int jobnum,
                       int K,
                       int rank,
                       boolean hasFrontendComm,
                       String frontendHost,
                       int frontendPort,
                       String backendHost)
            throws IOException {

        // Record command line arguments.
        this.username = username;
        this.jobnum = jobnum;
        this.K = K;
        this.rank = rank;
        this.hasFrontendComm = hasFrontendComm;
        this.frontendHost = frontendHost;
        this.frontendPort = frontendPort;
        this.backendHost = backendHost;

        // Turn on verbose start-up logging.
        StringBuilder sb = new StringBuilder();
        boolean verbose = Boolean.parseBoolean(System.getProperty("pj.verbose", "false"));
        if (verbose) {
            try {
                // Remove all log handlers from the default logger.
                Logger defaultLogger = LogManager.getLogManager().getLogger("");
                Handler[] defaultHandlers = defaultLogger.getHandlers();
                for (Handler h : defaultHandlers) {
                    defaultLogger.removeHandler(h);
                }

                // Create a FileHandler that logs messages with a SimpleFormatter.
                File dir = new File(Integer.toString(this.rank));
                if (!dir.exists()) {
                    boolean success = dir.mkdir();
                    sb.append(format("\n Creation of logging directory %s: %b", dir, success));
                }

                String logFile = dir.getAbsolutePath() + File.separator + "backend.log";
                fileHandler = new FileHandler(logFile);
                sb.append(format("\n Log file: %s", logFile));
                fileHandler.setFormatter(new SimpleFormatter());
                logger.addHandler(fileHandler);
                logger.setLevel(Level.INFO);
            } catch (Exception e) {
                logger.setLevel(Level.OFF);
            }
        } else {
            logger.setLevel(Level.OFF);
        }

        // Note that logging is OFF if the "pj.verbose" property was not true.
        sb.append(format("\n Username:          %s", this.username));
        sb.append(format("\n Job number:        %d", this.jobnum));
        sb.append(format("\n Nodes:             %d", this.K));
        sb.append(format("\n Rank:              %d", this.rank));
        sb.append(format("\n Has Frontend Comm: %b", this.hasFrontendComm));
        sb.append(format("\n Frontend Host:     %s", this.frontendHost));
        sb.append(format("\n Frontend Port:     %d", this.frontendPort));
        sb.append(format("\n Backend Host:      %s", this.backendHost));
        logger.log(Level.INFO, sb.toString());

        // Set up shutdown hook.
        Runtime.getRuntime().addShutdownHook(new Thread(() -> shutdown()));

        // Set up lease timer thread.
        logger.log(Level.INFO, " Set up lease timer thread.");
        myLeaseTimerThread = new TimerThread();
        myLeaseTimerThread.setDaemon(true);
        myLeaseTimerThread.start();

        // Set up job frontend lease timers.
        logger.log(Level.INFO, format(" Create frontend renew timer (%d sec)." , Constants.LEASE_RENEW_INTERVAL / 1000));
        myFrontendRenewTimer = myLeaseTimerThread.createTimer(timer -> {
                    try {
                        frontendRenewTimeout();
                    } catch (Throwable exc) {
                    }
                });
        logger.log(Level.INFO, format(" Create frontend expire timer (%d sec)." , Constants.LEASE_EXPIRE_INTERVAL / 1000));
        myFrontendExpireTimer = myLeaseTimerThread.createTimer(timer -> {
                    try {
                        frontendExpireTimeout();
                    } catch (Throwable exc) {
                    }
                });

        // Start job frontend lease expiration timer regardless of whether the
        // job frontend proxy gets set up.
        myFrontendExpireTimer.start(Constants.LEASE_EXPIRE_INTERVAL);

        // Set up middleware channel group.
        logger.log(Level.INFO, " Set up middleware channel group.");
        myMiddlewareChannelGroup = new ChannelGroup(new InetSocketAddress(backendHost, 0));
        myMiddlewareChannelGroup.startListening();

        // Set up job frontend proxy.
        logger.log(Level.INFO, " Set up job frontend proxy.");
        myJobFrontend = new JobFrontendProxy(myMiddlewareChannelGroup,
                myMiddlewareChannelGroup.connect(new InetSocketAddress(frontendHost, frontendPort)));

        // If we get here, the job frontend proxy has been set up.
        logger.log(Level.INFO, " The job frontend proxy has been set up.");

        // Start job frontend lease renewal timer.
        logger.log(Level.INFO, " Start frontend lease renewal timer.");
        myFrontendRenewTimer.start(Constants.LEASE_RENEW_INTERVAL,
                Constants.LEASE_RENEW_INTERVAL);

        // Set up backend class loader.
        myResourceCache = new ResourceCache();
        myClassLoader = new BackendClassLoader(
                /*parent        */ getClass().getClassLoader(),
                /*theJobBackend */ this,
                /*theJobFrontend*/ myJobFrontend,
                /*theCache      */ myResourceCache);

        // Set up world communicator channel group.
        logger.log(Level.INFO, " Set up world communicator channel group.");
        myWorldChannelGroup= new ChannelGroup(new InetSocketAddress(backendHost, 0));
        myWorldChannelGroup.setAlternateClassLoader(myClassLoader);

        // Set up frontend communicator channel group.
        if (hasFrontendComm) {
            logger.log(Level.INFO, " Set up frontend communicator channel group.");
            myFrontendChannelGroup = new ChannelGroup(new InetSocketAddress(backendHost, 0));
            myFrontendChannelGroup.setAlternateClassLoader(myClassLoader);
        }

        // Set up backend file writer and reader.
        logger.log(Level.INFO, " Set up backend file writer and reader.");
        myFileWriter = new BackendFileWriter(myJobFrontend, this);
        myFileReader = new BackendFileReader(myJobFrontend, this);

        // Redirect standard input, standard output, and standard error to job frontend.
        logger.log(Level.INFO, " Redirect standard input, standard output, and standard error to job frontend.");
        System.in.close();
        System.out.close();
        myJobLauncherLog = System.err;
        System.setIn(myFileReader.in);
        System.setOut(myFileWriter.out);
        System.setErr(myFileWriter.err);

        // Tell job frontend we're ready!
        logger.log(Level.INFO, " Tell job frontend we're ready.");
        myJobFrontend.backendReady(
                /*theJobBackend    */ this,
                /*rank             */ rank,
                /*middlewareAddress*/ myMiddlewareChannelGroup.listenAddress(),
                /*worldAddress     */ myWorldChannelGroup.listenAddress(),
                /*frontendAddress  */ hasFrontendComm ? myFrontendChannelGroup.listenAddress() : null);
    }
// Exported operations.

    /**
     * Run this Job Backend.
     */
    public void run() {
        Status status = null;
        JobBackendMessage message = null;

        try {
            while (continueRun) {
                // Receive a message from any channel.
                status
                        = myMiddlewareChannelGroup.receive(null, null, myBuffer);
                message = myBuffer.item;

                // Process message.
                message.invoke(this, myJobFrontend);

                // Enable garbage collection of no-longer-needed objects while
                // waiting to receive next message.
                myBuffer.item = null;
                status = null;
                message = null;
            }

            // Allow shutdown hook to proceed.
            reportRunFinished();
        } catch (ChannelGroupClosedException exc) {
            // Allow shutdown hook to proceed.
            reportRunFinished();
        } catch (Throwable exc) {
            // Allow shutdown hook to proceed.
            reportRunFinished();
            terminateCancelJob(exc);
        }

        // Exit process if necessary.
        switch (myState) {
            case TERMINATE_CANCEL_JOB:
            case TERMINATE_NO_REPORT:
                System.exit(1);
                break;
            case RUNNING:
            case TERMINATING:
                break;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Cancel the job.
     */
    public synchronized void cancelJob(JobFrontendRef theJobFrontend,
                                       String errmsg)
            throws IOException {
        terminateNoReport();
    }

    /**
     * Commence the job.
     *
     * @param theJobFrontend    Job Frontend that is calling this method.
     * @param middlewareAddress Array of hosts/ports for middleware messages.
     *                          The first <I>K</I>
     *                          elements are for the job backend processes in rank order, the
     *                          <I>K</I>+1st element is for the job frontend process. If the
     * @param worldAddress      Array of hosts/ports for the world communicator. The
     *                          <I>K</I>
     *                          elements are for the job backend processes in rank order.
     * @param frontendAddress   Array of hosts/ports for the frontend
     *                          communicator. The first
     *                          <I>K</I> elements are for the job backend processes in rank order, the
     *                          <I>K</I>+1st element is for the job frontend process. If the frontend
     *                          communicator does not exist, <code>frontendAddress</code> is null.
     * @param properties        Java system properties.
     * @param mainClassName     Fully qualified class name of the Java main program
     *                          class to execute.
     * @param args              Array of 0 or more Java command line arguments.
     * @throws java.io.IOException         Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public synchronized void commenceJob(JobFrontendRef theJobFrontend,
                                         InetSocketAddress[] middlewareAddress,
                                         InetSocketAddress[] worldAddress,
                                         InetSocketAddress[] frontendAddress,
                                         Properties properties,
                                         String mainClassName,
                                         String[] args)
            throws IOException {
        // Record information.
        myMiddlewareAddress = middlewareAddress;
        myWorldAddress = worldAddress;
        myFrontendAddress = frontendAddress;
        myProperties = properties;
        myMainClassName = mainClassName;
        myArgs = args;

        // Notify main program to commence job.
        commence = true;
        notifyAll();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report that the job finished.
     */
    public synchronized void jobFinished(JobFrontendRef theJobFrontend)
            throws IOException {
        continueRun = false;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Renew the lease on the job.
     */
    public synchronized void renewLease(JobFrontendRef theJobFrontend)
            throws IOException {
        myFrontendExpireTimer.start(Constants.LEASE_EXPIRE_INTERVAL);
    }

    /**
     * Report the content for a previously-requested resource.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param resourceName   Resource name.
     * @param content        Resource content, or null if resource not found.
     * @throws java.io.IOException         Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public synchronized void reportResource(JobFrontendRef theJobFrontend,
                                            String resourceName,
                                            byte[] content)
            throws IOException {
        myResourceCache.put(resourceName, content);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the content for a previously-requested resource.
     *
     * @throws java.io.IOException Thrown if an I/O error occurred.
     * @param theJobFrontend a {@link edu.rit.pj.cluster.JobFrontendRef} object
     * @param resourceName a {@link java.lang.String} object
     * @param content a {@link edu.rit.util.ByteSequence} object
     */
    public synchronized void reportResource(JobFrontendRef theJobFrontend,
                                            String resourceName,
                                            ByteSequence content)
            throws IOException {
        myResourceCache.put(resourceName,
                content == null ? null : content.toByteArray());
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of opening the given output file.
     */
    public synchronized void outputFileOpenResult(JobFrontendRef theJobFrontend,
                                                  int bfd,
                                                  int ffd,
                                                  IOException exc)
            throws IOException {
        myFileWriter.outputFileOpenResult(theJobFrontend, bfd, ffd, exc);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of writing the given output file.
     */
    public synchronized void outputFileWriteResult(JobFrontendRef theJobFrontend,
                                                   int ffd,
                                                   IOException exc)
            throws IOException {
        myFileWriter.outputFileWriteResult(theJobFrontend, ffd, exc);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of flushing the given output file.
     */
    public synchronized void outputFileFlushResult(JobFrontendRef theJobFrontend,
                                                   int ffd,
                                                   IOException exc)
            throws IOException {
        myFileWriter.outputFileFlushResult(theJobFrontend, ffd, exc);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of closing the given output file.
     */
    public synchronized void outputFileCloseResult(JobFrontendRef theJobFrontend,
                                                   int ffd,
                                                   IOException exc)
            throws IOException {
        myFileWriter.outputFileCloseResult(theJobFrontend, ffd, exc);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of opening the given input file.
     */
    public synchronized void inputFileOpenResult(JobFrontendRef theJobFrontend,
                                                 int bfd,
                                                 int ffd,
                                                 IOException exc)
            throws IOException {
        myFileReader.inputFileOpenResult(theJobFrontend, bfd, ffd, exc);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of reading the given input file.
     */
    public synchronized void inputFileReadResult(JobFrontendRef theJobFrontend,
                                                 int ffd,
                                                 byte[] buf,
                                                 int len,
                                                 IOException exc)
            throws IOException {
        myFileReader.inputFileReadResult(theJobFrontend, ffd, len, exc);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of skipping the given input file.
     */
    public synchronized void inputFileSkipResult(JobFrontendRef theJobFrontend,
                                                 int ffd,
                                                 long len,
                                                 IOException exc)
            throws IOException {
        myFileReader.inputFileSkipResult(theJobFrontend, ffd, len, exc);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Report the result of closing the given input file.
     */
    public synchronized void inputFileCloseResult(JobFrontendRef theJobFrontend,
                                                  int ffd,
                                                  IOException exc)
            throws IOException {
        myFileReader.inputFileCloseResult(theJobFrontend, ffd, exc);
    }

    /**
     * Close communication with this Job Backend.
     */
    public synchronized void close() {
    }

    /**
     * Obtain this job's user name.
     *
     * @return User name.
     */
    public String getUserName() {
        return username;
    }

    /**
     * Obtain this job's job number.
     *
     * @return Job number.
     */
    public int getJobNumber() {
        return jobnum;
    }

    /**
     * Obtain the number of backend processes in this job.
     *
     * @return <I>K</I>, the number of backend processes.
     */
    public int getK() {
        return K;
    }

    /**
     * Obtain the rank of this backend process in this job.
     *
     * @return Rank.
     */
    public int getRank() {
        return rank;
    }

    /**
     * Obtain the backend host name on which this job is running.
     *
     * @return Host name.
     */
    public String getBackendHost() {
        return backendHost;
    }

    /**
     * Determine whether the frontend communicator exists in this job.
     *
     * @return True if the frontend communicator exists, false if it doesn't.
     */
    public boolean hasFrontendCommunicator() {
        return hasFrontendComm;
    }

    /**
     * Obtain this job's backend class loader.
     *
     * @return Class loader.
     */
    public ClassLoader getClassLoader() {
        return myClassLoader;
    }

    /**
     * Obtain this job's backend file writer.
     *
     * @return Backend file writer.
     */
    public BackendFileWriter getFileWriter() {
        return myFileWriter;
    }

    /**
     * Obtain this job's backend file reader.
     *
     * @return Backend file reader.
     */
    public BackendFileReader getFileReader() {
        return myFileReader;
    }

    /**
     * Wait until this job commences.
     */
    public synchronized void waitForCommence() {
        while (!commence) {
            try {
                wait();
            } catch (InterruptedException exc) {
            }
        }
    }

    /**
     * Obtain this job's world communicator channel group. If this job has not
     * commenced yet, null is returned.
     *
     * @return Channel group.
     */
    public ChannelGroup getWorldChannelGroup() {
        return myWorldChannelGroup;
    }

    /**
     * Obtain this job's array of hosts/ports for the world communicator. The
     * <I>K</I> elements are for the job backend processes in rank order. If
     * this job has not commenced yet, null is returned.
     *
     * @return Array of world communicator addresses.
     */
    public InetSocketAddress[] getWorldAddress() {
        return myWorldAddress;
    }

    /**
     * Obtain this job's frontend communicator channel group. If the frontend
     * communicator does not exist, or if this job has not commenced yet, null
     * is returned.
     *
     * @return Channel group.
     */
    public ChannelGroup getFrontendChannelGroup() {
        return myFrontendChannelGroup;
    }

    /**
     * Obtain this job's array of hosts/ports for the frontend communicator. The
     * first <I>K</I> elements are for the job backend processes in rank order,
     * the <I>K</I>+1st element is for the job frontend process. If the frontend
     * communicator does not exist, or if this job has not commenced yet, null
     * is returned.
     *
     * @return Array of frontend communicator addresses.
     */
    public InetSocketAddress[] getFrontendAddress() {
        return myFrontendAddress;
    }

    /**
     * Obtain this job's Java system properties. If this job has not commenced
     * yet, null is returned.
     *
     * @return Properties.
     */
    public Properties getProperties() {
        return myProperties;
    }

    /**
     * Obtain this job's main class name. If this job has not commenced yet,
     * null is returned.
     *
     * @return Fully qualified class name of the Java main program class to
     * execute.
     */
    public String getMainClassName() {
        return myMainClassName;
    }

    /**
     * Obtain this job's command line arguments. If this job has not commenced
     * yet, null is returned.
     *
     * @return Array of 0 or more Java command line arguments.
     */
    public String[] getArgs() {
        return myArgs;
    }

    /**
     * Set the comment string for this job backend process. The comment string
     * appears in the detailed job status display in the Job Scheduler's web
     * interface. Each job backend process (rank) has its own comment string. If
     * <code>setComment()</code> is never called, the comment string is empty. The
     * comment string is typically used to display this job backend process's
     * progress. The comment string is rendered by a web browser and can contain
     * HTML tags.
     * <p>
     * Calling <code>setComment()</code> causes a message to be sent to the job
     * frontend process, which in turn causes a message to be sent to the Job
     * Scheduler. (Any I/O errors during message sending are ignored.)
     * Consequently, don't call <code>setComment()</code> too frequently, or the
     * program's performance will suffer.
     *
     * @param comment Comment string.
     */
    public void setComment(String comment) {
        try {
            myJobFrontend.reportComment(this, rank, comment);
        } catch (IOException exc) {
        }
    }

    /**
     * Obtain the Job Backend object. If the Job Backend main program is
     * running, the job backend object for the job is returned. If some other
     * main program is running, null is returned.
     *
     * @return Job backend object, or null.
     */
    public static JobBackend getJobBackend() {
        return theJobBackend;
    }

// More hidden operations.

    /**
     * Take action when the job frontend's lease renewal timer times out.
     *
     * @throws IOException Thrown if an I/O error occurred.
     */
    private synchronized void frontendRenewTimeout()
            throws IOException {
        if (myFrontendRenewTimer.isTriggered()) {
            myJobFrontend.renewLease(this);
        }
    }

    /**
     * Take action when the job frontend's lease expiration timer times out.
     *
     * @throws IOException Thrown if an I/O error occurred.
     */
    private void frontendExpireTimeout()
            throws IOException {
        boolean doExit = false;
        synchronized (this) {
            if (myFrontendExpireTimer.isTriggered()) {
                reportRunFinished();
                if (myState == State.RUNNING) {
                    myState = State.TERMINATE_NO_REPORT;
                    doExit = true;
                }
            }
        }

        // Cannot hold the synchronization lock while calling System.exit(),
        // otherwise a deadlock can occur between this thread (the timer thread)
        // and the shutdown hook thread.
        myJobLauncherLog.println("Job frontend lease expired");
        if (doExit) {
            System.exit(1);
        }
    }

    /**
     * Terminate this Job Backend immediately, sending a "cancel job" message to
     * the Job Frontend. The error message comes from the given exception.
     *
     * @param exc Exception.
     */
    private void terminateCancelJob(Throwable exc) {
        continueRun = false;
        if (myState == State.RUNNING) {
            myState = State.TERMINATE_CANCEL_JOB;
            myCancelMessage = exc.getClass().getName();
            String msg = exc.getMessage();
            if (msg != null) {
                myCancelMessage = myCancelMessage + ": " + msg;
            }
            //System.err.println (myCancelMessage);
            //exc.printStackTrace (System.err);
        }
    }

    /**
     * Terminate this Job Backend immediately, with no report to the Job
     * Frontend.
     */
    private void terminateNoReport() {
        continueRun = false;
        if (myState == State.RUNNING) {
            myState = State.TERMINATE_NO_REPORT;
        }
    }

    /**
     * Shut down this Job Backend.
     */
    private void shutdown() {
        synchronized (this) {
            // Tell job frontend that we are terminating.
            if (myJobFrontend != null) {
                try {
                    switch (myState) {
                        case RUNNING:
                            // Tell job frontend we finished normally.
                            myJobFrontend.backendFinished(this);
                            break;
                        case TERMINATE_CANCEL_JOB:
                            // Tell job frontend we're canceling.
                            myJobFrontend.cancelJob(this, myCancelMessage);
                            break;
                        case TERMINATE_NO_REPORT:
                        case TERMINATING:
                            // Tell job frontend nothing.
                            break;
                    }
                } catch (IOException exc) {
                }
            }

            // Record that we are terminating.
            myState = State.TERMINATING;
        }

        // Wait until the run() method thread terminates.
        waitForRunFinished();

        // Shut down job frontend lease timers.
        synchronized (this) {
            myFrontendRenewTimer.stop();
            myFrontendExpireTimer.stop();
        }

        // All proxies, channels, and channel groups will close when the process
        // exits.
    }

    /**
     * Wait for the run() method to finish.
     */
    private void waitForRunFinished() {
        for (; ; ) {
            try {
                runFinished.await();
                break;
            } catch (InterruptedException exc) {
            }
        }
    }

    /**
     * Report that the run() method finished.
     */
    private void reportRunFinished() {
        runFinished.countDown();
    }

    /**
     * Dump this job backend to the standard output, for debugging.
     */
    private synchronized void dump() {
        System.out.println("********************************");
        System.out.println("username = " + username);
        System.out.println("jobnum = " + jobnum);
        System.out.println("K = " + K);
        System.out.println("rank = " + rank);
        System.out.println("hasFrontendComm = " + hasFrontendComm);
        for (int i = 0; i <= K; ++i) {
            System.out.println("myMiddlewareAddress[" + i + "] = " + myMiddlewareAddress[i]);
        }
        for (int i = 0; i < K; ++i) {
            System.out.println("myWorldAddress[" + i + "] = " + myWorldAddress[i]);
        }
        if (hasFrontendComm) {
            for (int i = 0; i <= K; ++i) {
                System.out.println("myFrontendAddress[" + i + "] = " + myFrontendAddress[i]);
            }
        }
        myProperties.list(System.out);
        System.out.println("myMainClassName = " + myMainClassName);
        for (int i = 0; i < myArgs.length; ++i) {
            System.out.println("myArgs[" + i + "] = \"" + myArgs[i] + "\"");
        }
    }

// Main program.

    /**
     * Job Backend main program.
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.lang.Exception if any.
     */
    public static void main(String[] args)
            throws Exception {
        try {
            // Parse command line arguments.
            if (args.length != 8) {
                usage();
            }
            String username = args[0];
            int jobnum = Integer.parseInt(args[1]);
            int K = Integer.parseInt(args[2]);
            int rank = Integer.parseInt(args[3]);
            boolean hasFrontendComm = Boolean.parseBoolean(args[4]);
            String frontendHost = args[5];
            int frontendPort = Integer.parseInt(args[6]);
            String backendHost = args[7];

            // Set up job backend object.
            theJobBackend
                    = new JobBackend(username, jobnum, K, rank, hasFrontendComm,
                    frontendHost, frontendPort, backendHost);
        } catch (Throwable exc) {
            exc.printStackTrace(System.err);
            System.exit(1);
        }

        // Set the main thread's context class loader to be the job backend's
        // class loader.
        Thread.currentThread().setContextClassLoader(theJobBackend.getClassLoader());

        // Run job backend object in a separate thread.
        logger.log(Level.INFO, " Starting backend Daemon thread.");
        Thread thr = new Thread(theJobBackend);
        thr.setDaemon(true);
        thr.start();

        // Wait until job commences.
        logger.log(Level.INFO, " Waiting for commence.");
        theJobBackend.waitForCommence();
        logger.log(Level.INFO, " Commencing.");

        // Add any Java system properties from the job frontend process that do
        // not exist in this job backend process.
        Properties backendProperties = System.getProperties();
        Properties frontendProperties = theJobBackend.getProperties();
        for (Map.Entry<Object, Object> entry : frontendProperties.entrySet()) {
            String name = (String) entry.getKey();
            String value = (String) entry.getValue();
            if (backendProperties.getProperty(name) == null) {
                backendProperties.setProperty(name, value);
            }
        }

        // Turn on headless mode. This allows graphics drawing operations (that
        // do not require a screen, keyboard, or mouse) to work.
        System.setProperty("java.awt.headless", "true");

        // Call the job's main() method, passing in the job's command line
        // arguments.
        Class<?> mainclass
                = Class.forName(theJobBackend.getMainClassName(),
                true,
                // Force Field X modification to use our SystemClassLoader.
                ClassLoader.getSystemClassLoader());
        // Original:
        // theJobBackend.getClassLoader());
        Method mainmethod = mainclass.getMethod("main", String[].class);

        StringBuilder stringBuilder = new StringBuilder(format(" Preparing to invoke main method:\n %s\n With args: ", mainmethod.toString()));

        for (String arg : theJobBackend.getArgs()) {
            stringBuilder.append(format(" %s", arg));
        }
        stringBuilder.append("\n Backend start-up completed.");
        logger.log(Level.INFO, stringBuilder.toString());

        // Close down the FileHandler if it exists.
        if (theJobBackend.fileHandler != null) {
            theJobBackend.fileHandler.flush();
            logger.setLevel(Level.OFF);
            theJobBackend.fileHandler.close();
        }

        mainmethod.invoke(null, (Object) theJobBackend.getArgs());

        // After the main() method returns and all non-daemon threads have
        // terminated, the process will exit, and the shutdown hook will call
        // the shutdown() method.
    }

    /**
     * Print a usage message and exit.
     */
    private static void usage() {
        System.err.println("Usage: java edu.rit.pj.cluster.JobBackend <username> <jobnum> <K> <rank> <hasFrontendComm> <frontendHost> <frontendPort> <backendHost>");
        System.exit(1);
    }

}
