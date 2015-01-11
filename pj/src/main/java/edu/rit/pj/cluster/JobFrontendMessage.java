//******************************************************************************
//
// File:    JobFrontendMessage.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobFrontendMessage
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
import java.io.Externalizable;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.net.InetSocketAddress;

/**
 * Class JobFrontendMessage provides a message sent to a Job Frontend process
 * (interface {@linkplain JobFrontendRef}) in the PJ cluster middleware.
 *
 * @author Alan Kaminsky
 * @version 20-Jun-2012
 */
public abstract class JobFrontendMessage
        extends Message
        implements Externalizable {

// Hidden data members.
    private static final long serialVersionUID = -6601793901631997673L;

// Exported constructors.
    /**
     * Construct a new job frontend message.
     */
    public JobFrontendMessage() {
    }

    /**
     * Construct a new job frontend message with the given message tag.
     *
     * @param theTag Message tag to use when sending this message.
     */
    public JobFrontendMessage(int theTag) {
        super(theTag);
    }

// Exported operations.
    /**
     * Construct a new "assign backend" message.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @param name Backend node name.
     * @param host Host name for SSH remote login.
     * @param jvm Full pathname of Java Virtual Machine.
     * @param classpath Java class path for PJ Library.
     * @param jvmflags Array of JVM command line flags.
     * @param shellCommand Shell command string.
     * @param Nt Number of CPUs assigned to the process.
     *
     * @return "Assign backend" message.
     */
    public static JobFrontendMessage assignBackend(JobSchedulerRef theJobScheduler,
            String name,
            String host,
            String jvm,
            String classpath,
            String[] jvmflags,
            String shellCommand,
            int Nt) {
        return new AssignBackendMessage(theJobScheduler, name, host, jvm, classpath, jvmflags,
                shellCommand, Nt);
    }

    /**
     * Construct a new "assign job number" message.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @param jobnum Job number.
     * @param pjhost Host name for middleware channel group.
     *
     * @return "Assign job number" message.
     */
    public static JobFrontendMessage assignJobNumber(JobSchedulerRef theJobScheduler,
            int jobnum,
            String pjhost) {
        return new AssignJobNumberMessage(theJobScheduler, jobnum, pjhost);
    }

    /**
     * Construct a new "cancel job" message.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     * @param errmsg Error message string.
     *
     * @return "Cancel job" message.
     */
    public static JobFrontendMessage cancelJob(JobSchedulerRef theJobScheduler,
            String errmsg) {
        return new CancelJobMessage(theJobScheduler, errmsg);
    }

    /**
     * Construct a new "renew lease" message.
     *
     * @param theJobScheduler Job Scheduler that is calling this method.
     *
     * @return "Renew lease" message.
     */
    public static JobFrontendMessage renewLease(JobSchedulerRef theJobScheduler) {
        return new RenewLeaseMessage(theJobScheduler);
    }

    /**
     * Construct a new "backend finished" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     *
     * @return "Backend finished" message.
     */
    public static JobFrontendMessage backendFinished(JobBackendRef theJobBackend) {
        return new BackendFinishedMessage(theJobBackend);
    }

    /**
     * Construct a new "backend ready" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param rank Rank of the job backend process.
     * @param middlewareAddress Host/port to which the job backend process is
     * listening for middleware messages.
     * @param worldAddress Host/port to which the job backend process is
     * listening for the world communicator.
     * @param frontendAddress Host/port to which the job backend process is
     * listening for the frontend communicator, or null if the frontend
     * communicator does not exist.
     *
     * @return "Backend ready" message.
     */
    public static JobFrontendMessage backendReady(JobBackendRef theJobBackend,
            int rank,
            InetSocketAddress middlewareAddress,
            InetSocketAddress worldAddress,
            InetSocketAddress frontendAddress) {
        return new BackendReadyMessage(theJobBackend, rank, middlewareAddress,
                worldAddress, frontendAddress);
    }

    /**
     * Construct a new "cancel job" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param errmsg Error message string.
     *
     * @return "Cancel job" message.
     */
    public static JobFrontendMessage cancelJob(JobBackendRef theJobBackend,
            String errmsg) {
        return new CancelJobMessage(theJobBackend, errmsg);
    }

    /**
     * Construct a new "renew lease" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     *
     * @return "Renew lease" message.
     */
    public static JobFrontendMessage renewLease(JobBackendRef theJobBackend) {
        return new RenewLeaseMessage(theJobBackend);
    }

    /**
     * Construct a new "request resource" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param resourceName Resource name.
     *
     * @return "Request resource" message.
     */
    public static JobFrontendMessage requestResource(JobBackendRef theJobBackend,
            String resourceName) {
        return new RequestResourceMessage(theJobBackend, resourceName);
    }

    /**
     * Construct a new "output file open" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param bfd Backend file descriptor.
     * @param file File.
     * @param append True to append, false to overwrite.
     *
     * @return "Output file open" message.
     */
    public static JobFrontendMessage outputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file,
            boolean append) {
        return new OutputFileOpenMessage(theJobBackend, bfd, file, append);
    }

    /**
     * Construct a new "output file write" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes to write.
     *
     * @return "Output file write" message.
     */
    public static JobFrontendMessage outputFileWrite(JobBackendRef theJobBackend,
            int ffd,
            int len) {
        return new OutputFileWriteMessage(theJobBackend, ffd, len);
    }

    /**
     * Construct a new "output file flush" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     *
     * @return "Flush file" message.
     */
    public static JobFrontendMessage outputFileFlush(JobBackendRef theJobBackend,
            int ffd) {
        return new OutputFileFlushMessage(theJobBackend, ffd);
    }

    /**
     * Construct a new "output file close" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     *
     * @return "Output file close" message.
     */
    public static JobFrontendMessage outputFileClose(JobBackendRef theJobBackend,
            int ffd) {
        return new OutputFileCloseMessage(theJobBackend, ffd);
    }

    /**
     * Construct a new "input file open" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param bfd Backend file descriptor.
     * @param file File.
     *
     * @return "Input file open" message.
     */
    public static JobFrontendMessage inputFileOpen(JobBackendRef theJobBackend,
            int bfd,
            File file) {
        return new InputFileOpenMessage(theJobBackend, bfd, file);
    }

    /**
     * Construct a new "input file read" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes to read.
     *
     * @return "Input file read" message.
     */
    public static JobFrontendMessage inputFileRead(JobBackendRef theJobBackend,
            int ffd,
            int len) {
        return new InputFileReadMessage(theJobBackend, ffd, len);
    }

    /**
     * Construct a new "input file skip" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes to skip.
     *
     * @return "Input file skip" message.
     */
    public static JobFrontendMessage inputFileSkip(JobBackendRef theJobBackend,
            int ffd,
            long len) {
        return new InputFileSkipMessage(theJobBackend, ffd, len);
    }

    /**
     * Construct a new "input file close" message.
     *
     * @param theJobBackend Job Backend that is calling this method.
     * @param ffd Frontend file descriptor.
     *
     * @return "Input file close" message.
     */
    public static JobFrontendMessage inputFileClose(JobBackendRef theJobBackend,
            int ffd) {
        return new InputFileCloseMessage(theJobBackend, ffd);
    }

    /**
     * Construct a new "report comment" message.
     *
     * @param theJobBackend Job backend that is calling this method.
     * @param rank Process rank.
     * @param comment Comment string.
     *
     * @return "Report comment" message.
     */
    public static JobFrontendMessage reportComment(JobBackendRef theJobBackend,
            int rank,
            String comment) {
        return new ReportCommentMessage(theJobBackend, rank, comment);
    }

    /**
     * Invoke the method corresponding to this job frontend message on the given
     * Job Frontend object. The method arguments come from the fields of this
     * job frontend message object.
     *
     * @param theJobFrontend Job Frontend on which to invoke the method.
     * @param theJobScheduler Job Scheduler that is calling the method.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void invoke(JobFrontendRef theJobFrontend,
            JobSchedulerRef theJobScheduler)
            throws IOException {
        throw new UnsupportedOperationException();
    }

    /**
     * Invoke the method corresponding to this job frontend message on the given
     * Job Frontend object. The method arguments come from the fields of this
     * job frontend message object.
     *
     * @param theJobFrontend Job Frontend on which to invoke the method.
     * @param theJobBackend Job Backend that is calling the method.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void invoke(JobFrontendRef theJobFrontend,
            JobBackendRef theJobBackend)
            throws IOException {
        throw new UnsupportedOperationException();
    }

    /**
     * Write this job frontend message to the given object output stream.
     *
     * @param out Object output stream.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void writeExternal(ObjectOutput out)
            throws IOException {
    }

    /**
     * Read this job frontend message from the given object input stream.
     *
     * @param in Object input stream.
     *
     * @exception IOException Thrown if an I/O error occurred.
     * @exception ClassNotFoundException Thrown if a class needed to read this
     * job backend message could not be found.
     */
    public void readExternal(ObjectInput in)
            throws IOException, ClassNotFoundException {
    }

// Hidden subclasses.
    /**
     * Class AssignBackendMessage provides the Job Frontend "assign backend"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 20-Jun-2008
     */
    private static class AssignBackendMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 4426930548078115598L;

        private String name;
        private String host;
        private String jvm;
        private String classpath;
        private String[] jvmflags;
        private String shellCommand;
        private int Nt;

        public AssignBackendMessage() {
        }

        public AssignBackendMessage(JobSchedulerRef theJobScheduler,
                String name,
                String host,
                String jvm,
                String classpath,
                String[] jvmflags,
                String shellCommand,
                int Nt) {
            super(Message.FROM_JOB_SCHEDULER);
            this.name = name;
            this.host = host;
            this.jvm = jvm;
            this.classpath = classpath;
            this.jvmflags = jvmflags;
            this.shellCommand = shellCommand;
            this.Nt = Nt;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobSchedulerRef theJobScheduler)
                throws IOException {
            theJobFrontend.assignBackend(theJobScheduler, name, host, jvm, classpath, jvmflags,
                    shellCommand, Nt);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeUTF(name);
            out.writeUTF(host);
            out.writeUTF(jvm);
            out.writeUTF(classpath);
            int n = jvmflags.length;
            out.writeInt(n);
            for (int i = 0; i < n; ++i) {
                out.writeUTF(jvmflags[i]);
            }
            out.writeUTF(shellCommand);
            out.writeInt(Nt);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            name = in.readUTF();
            host = in.readUTF();
            jvm = in.readUTF();
            classpath = in.readUTF();
            int n = in.readInt();
            jvmflags = new String[n];
            for (int i = 0; i < n; ++i) {
                jvmflags[i] = in.readUTF();
            }
            shellCommand = in.readUTF();
            Nt = in.readInt();
        }
    }

    /**
     * Class AssignJobNumberMessage provides the Job Frontend "assign job
     * number" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 12-Oct-2006
     */
    private static class AssignJobNumberMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -197099388467750972L;

        private int jobnum;
        private String pjhost;

        public AssignJobNumberMessage() {
        }

        public AssignJobNumberMessage(JobSchedulerRef theJobScheduler,
                int jobnum,
                String pjhost) {
            super(Message.FROM_JOB_SCHEDULER);
            this.jobnum = jobnum;
            this.pjhost = pjhost;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobSchedulerRef theJobScheduler)
                throws IOException {
            theJobFrontend.assignJobNumber(theJobScheduler, jobnum, pjhost);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(jobnum);
            out.writeUTF(pjhost);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            jobnum = in.readInt();
            pjhost = in.readUTF();
        }
    }

    /**
     * Class BackendFinishedMessage provides the Job Frontend "backend finished"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 12-Oct-2006
     */
    private static class BackendFinishedMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 759872032212034107L;

        public BackendFinishedMessage() {
        }

        public BackendFinishedMessage(JobBackendRef theJobBackend) {
            super(Message.FROM_JOB_BACKEND);
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.backendFinished(theJobBackend);
        }
    }

    /**
     * Class BackendReadyMessage provides the Job Frontend "backend ready"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 20-Oct-2006
     */
    private static class BackendReadyMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -8872540352133660209L;

        private int rank;
        private InetSocketAddress middlewareAddress;
        private InetSocketAddress worldAddress;
        private InetSocketAddress frontendAddress;

        public BackendReadyMessage() {
        }

        public BackendReadyMessage(JobBackendRef theJobBackend,
                int rank,
                InetSocketAddress middlewareAddress,
                InetSocketAddress worldAddress,
                InetSocketAddress frontendAddress) {
            super(Message.FROM_JOB_BACKEND);
            this.rank = rank;
            this.middlewareAddress = middlewareAddress;
            this.worldAddress = worldAddress;
            this.frontendAddress = frontendAddress;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.backendReady(theJobBackend, rank, middlewareAddress,
                    worldAddress, frontendAddress);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(rank);
            out.writeObject(middlewareAddress);
            out.writeObject(worldAddress);
            out.writeObject(frontendAddress);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            rank = in.readInt();
            middlewareAddress = (InetSocketAddress) in.readObject();
            worldAddress = (InetSocketAddress) in.readObject();
            frontendAddress = (InetSocketAddress) in.readObject();
        }
    }

    /**
     * Class CancelJobMessage provides the Job Frontend "cancel job" message in
     * the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 12-Oct-2006
     */
    private static class CancelJobMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -2009079595329812496L;

        private String errmsg;

        public CancelJobMessage() {
        }

        public CancelJobMessage(JobSchedulerRef theJobScheduler,
                String errmsg) {
            super(Message.FROM_JOB_SCHEDULER);
            this.errmsg = errmsg;
        }

        public CancelJobMessage(JobBackendRef theJobBackend,
                String errmsg) {
            super(Message.FROM_JOB_BACKEND);
            this.errmsg = errmsg;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobSchedulerRef theJobScheduler)
                throws IOException {
            theJobFrontend.cancelJob(theJobScheduler, errmsg);
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.cancelJob(theJobBackend, errmsg);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeUTF(errmsg);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            errmsg = in.readUTF();
        }
    }

    /**
     * Class RenewLeaseMessage provides the Job Frontend "renew lease" message
     * in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 12-Oct-2006
     */
    private static class RenewLeaseMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -9030468939160436622L;

        public RenewLeaseMessage() {
        }

        public RenewLeaseMessage(JobSchedulerRef theJobScheduler) {
            super(Message.FROM_JOB_SCHEDULER);
        }

        public RenewLeaseMessage(JobBackendRef theJobBackend) {
            super(Message.FROM_JOB_BACKEND);
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobSchedulerRef theJobScheduler)
                throws IOException {
            theJobFrontend.renewLease(theJobScheduler);
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.renewLease(theJobBackend);
        }
    }

    /**
     * Class RequestResourceMessage provides the Job Frontend "request resource"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 26-Oct-2006
     */
    private static class RequestResourceMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 9184806604713339263L;

        private String resourceName;

        public RequestResourceMessage() {
        }

        public RequestResourceMessage(JobBackendRef theJobBackend,
                String resourceName) {
            super(Message.FROM_JOB_BACKEND);
            this.resourceName = resourceName;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.requestResource(theJobBackend, resourceName);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeUTF(resourceName);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            resourceName = in.readUTF();
        }
    }

    /**
     * Class OutputFileOpenMessage provides the Job Frontend "output file open"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class OutputFileOpenMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -4987754930269039852L;

        private int bfd;
        private File file;
        private boolean append;

        public OutputFileOpenMessage() {
        }

        public OutputFileOpenMessage(JobBackendRef theJobBackend,
                int bfd,
                File file,
                boolean append) {
            super(Message.FROM_JOB_BACKEND);
            this.bfd = bfd;
            this.file = file;
            this.append = append;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.outputFileOpen(theJobBackend, bfd, file, append);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(bfd);
            out.writeBoolean(append);
            out.writeObject(file);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            bfd = in.readInt();
            append = in.readBoolean();
            file = (File) in.readObject();
        }
    }

    /**
     * Class OutputFileWriteMessage provides the Job Frontend "output file
     * write" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 20-Nov-2006
     */
    private static class OutputFileWriteMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -4460426636308841602L;

        private int ffd;
        private int len;

        public OutputFileWriteMessage() {
        }

        public OutputFileWriteMessage(JobBackendRef theJobBackend,
                int ffd,
                int len) {
            super(Message.FROM_JOB_BACKEND);
            this.ffd = ffd;
            this.len = len;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.outputFileWrite(theJobBackend, ffd, null, 0, len);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeInt(len);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            ffd = in.readInt();
            len = in.readInt();
        }
    }

    /**
     * Class OutputFileFlushMessage provides the Job Frontend "output file
     * flush" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class OutputFileFlushMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 7074849708663078210L;

        private int ffd;

        public OutputFileFlushMessage() {
        }

        public OutputFileFlushMessage(JobBackendRef theJobBackend,
                int ffd) {
            super(Message.FROM_JOB_BACKEND);
            this.ffd = ffd;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.outputFileFlush(theJobBackend, ffd);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            ffd = in.readInt();
        }
    }

    /**
     * Class OutputFileCloseMessage provides the Job Frontend "output file
     * close" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class OutputFileCloseMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -5637017577338427153L;

        private int ffd;

        public OutputFileCloseMessage() {
        }

        public OutputFileCloseMessage(JobBackendRef theJobBackend,
                int ffd) {
            super(Message.FROM_JOB_BACKEND);
            this.ffd = ffd;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.outputFileClose(theJobBackend, ffd);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            ffd = in.readInt();
        }
    }

    /**
     * Class InputFileOpenMessage provides the Job Frontend "input file open"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class InputFileOpenMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = -791306998166025239L;

        private int bfd;
        private File file;

        public InputFileOpenMessage() {
        }

        public InputFileOpenMessage(JobBackendRef theJobBackend,
                int bfd,
                File file) {
            super(Message.FROM_JOB_BACKEND);
            this.bfd = bfd;
            this.file = file;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.inputFileOpen(theJobBackend, bfd, file);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(bfd);
            out.writeObject(file);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            bfd = in.readInt();
            file = (File) in.readObject();
        }
    }

    /**
     * Class InputFileReadMessage provides the Job Frontend "input file read"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class InputFileReadMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 7727558874589005187L;

        private int ffd;
        private int len;

        public InputFileReadMessage() {
        }

        public InputFileReadMessage(JobBackendRef theJobBackend,
                int ffd,
                int len) {
            super(Message.FROM_JOB_BACKEND);
            this.ffd = ffd;
            this.len = len;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.inputFileRead(theJobBackend, ffd, len);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeInt(len);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            ffd = in.readInt();
            len = in.readInt();
        }
    }

    /**
     * Class InputFileSkipMessage provides the Job Frontend "input file skip"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class InputFileSkipMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 7867427744015166954L;

        private int ffd;
        private long len;

        public InputFileSkipMessage() {
        }

        public InputFileSkipMessage(JobBackendRef theJobBackend,
                int ffd,
                long len) {
            super(Message.FROM_JOB_BACKEND);
            this.ffd = ffd;
            this.len = len;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.inputFileSkip(theJobBackend, ffd, len);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeLong(len);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            ffd = in.readInt();
            len = in.readLong();
        }
    }

    /**
     * Class InputFileCloseMessage provides the Job Frontend "input file close"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class InputFileCloseMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 1223549294646718409L;

        private int ffd;

        public InputFileCloseMessage() {
        }

        public InputFileCloseMessage(JobBackendRef theJobBackend,
                int ffd) {
            super(Message.FROM_JOB_BACKEND);
            this.ffd = ffd;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.inputFileClose(theJobBackend, ffd);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            ffd = in.readInt();
        }
    }

    /**
     * Class ReportCommentMessage provides the Job Frontend "report comment"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 24-Jan-2012
     */
    private static class ReportCommentMessage
            extends JobFrontendMessage {

        private static final long serialVersionUID = 1092254806177598252L;

        private int rank;
        private String comment;

        public ReportCommentMessage() {
        }

        public ReportCommentMessage(JobBackendRef theJobBackend,
                int rank,
                String comment) {
            super(Message.FROM_JOB_BACKEND);
            this.rank = rank;
            this.comment = comment == null ? "" : comment;
        }

        public void invoke(JobFrontendRef theJobFrontend,
                JobBackendRef theJobBackend)
                throws IOException {
            theJobFrontend.reportComment(theJobBackend, rank, comment);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(rank);
            out.writeUTF(comment);
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            rank = in.readInt();
            comment = in.readUTF();
        }
    }

}
