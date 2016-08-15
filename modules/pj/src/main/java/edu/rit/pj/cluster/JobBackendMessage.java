//******************************************************************************
//
// File:    JobBackendMessage.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.JobBackendMessage
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

import edu.rit.util.ByteSequence;

import java.io.IOException;
import java.io.Externalizable;
import java.io.ObjectInput;
import java.io.ObjectOutput;

import java.net.InetSocketAddress;

import java.util.Map;
import java.util.Properties;

/**
 * Class JobBackendMessage provides a message sent to a Job Backend process
 * (interface {@linkplain JobBackendRef}) in the PJ cluster middleware.
 *
 * @author Alan Kaminsky
 * @version 02-Jul-2007
 */
public abstract class JobBackendMessage
        extends Message
        implements Externalizable {

// Hidden data members.
    private static final long serialVersionUID = 3747140854797048519L;

// Exported constructors.
    /**
     * Construct a new job backend message.
     */
    public JobBackendMessage() {
    }

    /**
     * Construct a new job backend message with the given message tag.
     *
     * @param theTag Message tag to use when sending this message.
     */
    public JobBackendMessage(int theTag) {
        super(theTag);
    }

// Exported operations.
    /**
     * Construct a new "cancel job" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param errmsg Error message string.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage cancelJob(JobFrontendRef theJobFrontend,
            String errmsg) {
        return new CancelJobMessage(theJobFrontend, errmsg);
    }

    /**
     * Construct a new "commence job" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param middlewareAddress Array of hosts/ports for middleware messages.
     * The first <I>K</I>
     * elements are for the job backend processes in rank order, the
     * <I>K</I>+1st element is for the job frontend process. If the
     * @param worldAddress Array of hosts/ports for the world communicator. The
     * <I>K</I>
     * elements are for the job backend processes in rank order.
     * @param frontendAddress Array of hosts/ports for the frontend
     * communicator. The first
     * <I>K</I> elements are for the job backend processes in rank order, the
     * <I>K</I>+1st element is for the job frontend process. If the frontend
     * communicator does not exist, <TT>frontendAddress</TT> is null.
     * @param properties Java system properties.
     * @param mainClassName Fully qualified class name of the Java main program
     * class to execute.
     * @param args Array of 0 or more Java command line arguments.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage commenceJob(JobFrontendRef theJobFrontend,
            InetSocketAddress[] middlewareAddress,
            InetSocketAddress[] worldAddress,
            InetSocketAddress[] frontendAddress,
            Properties properties,
            String mainClassName,
            String[] args) {
        return new CommenceJobMessage(theJobFrontend, middlewareAddress, worldAddress, frontendAddress,
                properties, mainClassName, args);
    }

    /**
     * Construct a new "job finished" message.
     *
     * @param theJobFrontend Job frontend that is calling this method.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage jobFinished(JobFrontendRef theJobFrontend) {
        return new JobFinishedMessage(theJobFrontend);
    }

    /**
     * Construct a new "renew lease" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage renewLease(JobFrontendRef theJobFrontend) {
        return new RenewLeaseMessage(theJobFrontend);
    }

    /**
     * Construct a new "report resource" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param resourceName Resource name.
     * @param content Resource content, or null if resource not found.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage reportResource(JobFrontendRef theJobFrontend,
            String resourceName,
            byte[] content) {
        return new ReportResourceMessage(theJobFrontend, resourceName, content);
    }

    /**
     * Construct a new "report resource" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param resourceName Resource name.
     * @param content Resource content, or null if resource not found.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage reportResource(JobFrontendRef theJobFrontend,
            String resourceName,
            ByteSequence content) {
        return new ReportResourceMessage(theJobFrontend, resourceName, content);
    }

    /**
     * Construct a new "output file open result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param bfd Backend file descriptor.
     * @param ffd Frontend file descriptor if success.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage outputFileOpenResult(JobFrontendRef theJobFrontend,
            int bfd,
            int ffd,
            IOException exc) {
        return new OutputFileOpenResultMessage(theJobFrontend, bfd, ffd, exc);
    }

    /**
     * Construct a new "output file write result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage outputFileWriteResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc) {
        return new OutputFileWriteResultMessage(theJobFrontend, ffd, exc);
    }

    /**
     * Construct a new "output file flush result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage outputFileFlushResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc) {
        return new OutputFileFlushResultMessage(theJobFrontend, ffd, exc);
    }

    /**
     * Construct a new "output file close result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage outputFileCloseResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc) {
        return new OutputFileCloseResultMessage(theJobFrontend, ffd, exc);
    }

    /**
     * Construct a new "input file open result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param bfd Backend file descriptor.
     * @param ffd Frontend file descriptor if success.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage inputFileOpenResult(JobFrontendRef theJobFrontend,
            int bfd,
            int ffd,
            IOException exc) {
        return new InputFileOpenResultMessage(theJobFrontend, bfd, ffd, exc);
    }

    /**
     * Construct a new "input file read result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes read, or -1 if EOF.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage inputFileReadResult(JobFrontendRef theJobFrontend,
            int ffd,
            int len,
            IOException exc) {
        return new InputFileReadResultMessage(theJobFrontend, ffd, len, exc);
    }

    /**
     * Construct a new "input file skip result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param len Number of bytes skipped.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage inputFileSkipResult(JobFrontendRef theJobFrontend,
            int ffd,
            long len,
            IOException exc) {
        return new InputFileSkipResultMessage(theJobFrontend, ffd, len, exc);
    }

    /**
     * Construct a new "input file close result" message.
     *
     * @param theJobFrontend Job Frontend that is calling this method.
     * @param ffd Frontend file descriptor.
     * @param exc Null if success, exception if failure.
     * @return a {@link edu.rit.pj.cluster.JobBackendMessage} object.
     */
    public static JobBackendMessage inputFileCloseResult(JobFrontendRef theJobFrontend,
            int ffd,
            IOException exc) {
        return new InputFileCloseResultMessage(theJobFrontend, ffd, exc);
    }

    /**
     * Invoke the method corresponding to this job backend message on the given
     * Job Backend object. The method arguments come from the fields of this job
     * backend message object.
     *
     * @param theJobBackend Job Backend on which to invoke the method.
     * @param theJobFrontend Job Frontend that is calling the method.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void invoke(JobBackendRef theJobBackend,
            JobFrontendRef theJobFrontend)
            throws IOException {
        throw new UnsupportedOperationException();
    }

    /**
     * {@inheritDoc}
     *
     * Write this job backend message to the given object output stream.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void writeExternal(ObjectOutput out)
            throws IOException {
    }

    /**
     * {@inheritDoc}
     *
     * Read this job backend message from the given object input stream.
     * @exception IOException Thrown if an I/O error occurred.
     * @exception ClassNotFoundException Thrown if a class needed to read this
     * job backend message could not be found.
     */
    public void readExternal(ObjectInput in)
            throws IOException, ClassNotFoundException {
    }

// Hidden subclasses.
    /**
     * Class CancelJobMessage provides the Job Backend "cancel job" message in
     * the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 25-Oct-2006
     */
    private static class CancelJobMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = -1706674774429384654L;

        private String errmsg;

        public CancelJobMessage() {
        }

        public CancelJobMessage(JobFrontendRef theJobFrontend,
                String errmsg) {
            super(Message.FROM_JOB_FRONTEND);
            this.errmsg = errmsg;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.cancelJob(theJobFrontend, errmsg);
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
     * Class CommenceJobMessage provides the Job Backend "commence job" message
     * in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 25-Oct-2006
     */
    private static class CommenceJobMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = -8262872991140404870L;

        private InetSocketAddress[] middlewareAddress;
        private InetSocketAddress[] worldAddress;
        private InetSocketAddress[] frontendAddress;
        private Properties properties;
        private String mainClassName;
        private String[] args;

        public CommenceJobMessage() {
        }

        public CommenceJobMessage(JobFrontendRef theJobFrontend,
                InetSocketAddress[] middlewareAddress,
                InetSocketAddress[] worldAddress,
                InetSocketAddress[] frontendAddress,
                Properties properties,
                String mainClassName,
                String[] args) {
            super(Message.FROM_JOB_FRONTEND);
            this.middlewareAddress = middlewareAddress;
            this.worldAddress = worldAddress;
            this.frontendAddress = frontendAddress;
            this.properties = properties;
            this.mainClassName = mainClassName;
            this.args = args;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.commenceJob(theJobFrontend, middlewareAddress, worldAddress,
                    frontendAddress, properties, mainClassName, args);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(properties.size());
            for (Map.Entry<Object, Object> entry : properties.entrySet()) {
                out.writeUTF((String) entry.getKey());
                out.writeUTF((String) entry.getValue());
            }
            out.writeUTF(mainClassName);
            out.writeInt(args.length);
            for (String arg : args) {
                out.writeUTF(arg);
            }
            int n1 = middlewareAddress.length;
            int n2 = worldAddress.length;
            int n3 = frontendAddress == null ? 0 : frontendAddress.length;
            out.writeInt(n1);
            out.writeInt(n2);
            out.writeInt(n3);
            for (int i = 0; i < n1; ++i) {
                out.writeObject(middlewareAddress[i]);
            }
            for (int i = 0; i < n2; ++i) {
                out.writeObject(worldAddress[i]);
            }
            for (int i = 0; i < n3; ++i) {
                out.writeObject(frontendAddress[i]);
            }
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            int n = in.readInt();
            properties = new Properties();
            for (int i = 0; i < n; ++i) {
                properties.setProperty(in.readUTF(), in.readUTF());
            }
            mainClassName = in.readUTF();
            n = in.readInt();
            args = new String[n];
            for (int i = 0; i < n; ++i) {
                args[i] = in.readUTF();
            }
            int n1 = in.readInt();
            int n2 = in.readInt();
            int n3 = in.readInt();
            middlewareAddress = new InetSocketAddress[n1];
            for (int i = 0; i < n1; ++i) {
                middlewareAddress[i] = (InetSocketAddress) in.readObject();
            }
            worldAddress = new InetSocketAddress[n2];
            for (int i = 0; i < n2; ++i) {
                worldAddress[i] = (InetSocketAddress) in.readObject();
            }
            if (n3 > 0) {
                frontendAddress = new InetSocketAddress[n3];
                for (int i = 0; i < n3; ++i) {
                    frontendAddress[i] = (InetSocketAddress) in.readObject();
                }
            }
        }
    }

    /**
     * Class JobFinishedMessage provides the Job Backend "job finished" message
     * in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 26-Oct-2006
     */
    private static class JobFinishedMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = 1363549433797859519L;

        public JobFinishedMessage() {
        }

        public JobFinishedMessage(JobFrontendRef theJobFrontend) {
            super(Message.FROM_JOB_FRONTEND);
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.jobFinished(theJobFrontend);
        }
    }

    /**
     * Class RenewLeaseMessage provides the Job Backend "renew lease" message in
     * the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 25-Oct-2006
     */
    private static class RenewLeaseMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = -5146916532326180730L;

        public RenewLeaseMessage() {
        }

        public RenewLeaseMessage(JobFrontendRef theJobFrontend) {
            super(Message.FROM_JOB_FRONTEND);
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.renewLease(theJobFrontend);
        }
    }

    /**
     * Class ReportResourceMessage provides the Job Backend "report resource"
     * message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 26-Oct-2006
     */
    private static class ReportResourceMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = 8709827200639757997L;

        private String resourceName;
        private ByteSequence contentSeq;
        private byte[] content;

        public ReportResourceMessage() {
        }

        public ReportResourceMessage(JobFrontendRef theJobFrontend,
                String resourceName,
                byte[] content) {
            super(Message.FROM_JOB_FRONTEND);
            this.resourceName = resourceName;
            this.content = content;
        }

        public ReportResourceMessage(JobFrontendRef theJobFrontend,
                String resourceName,
                ByteSequence content) {
            super(Message.FROM_JOB_FRONTEND);
            this.resourceName = resourceName;
            this.contentSeq = content;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            if (content != null) {
                theJobBackend.reportResource(theJobFrontend, resourceName, content);
            } else {
                theJobBackend.reportResource(theJobFrontend, resourceName, contentSeq);
            }
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeUTF(resourceName);
            if (content != null) {
                out.writeInt(content.length);
                out.write(content);
            } else if (contentSeq != null) {
                out.writeInt(contentSeq.length());
                contentSeq.write(out);
            } else {
                out.writeInt(-1);
            }
        }

        public void readExternal(ObjectInput in)
                throws IOException {
            resourceName = in.readUTF();
            int n = in.readInt();
            if (n < 0) {
                content = null;
            } else {
                content = new byte[n];
                in.readFully(content);
            }
        }
    }

    /**
     * Class OutputFileOpenResultMessage provides the Job Backend "output file
     * open result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class OutputFileOpenResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = 1460222094425830409L;

        private int bfd;
        private int ffd;
        private IOException exc;

        public OutputFileOpenResultMessage() {
        }

        public OutputFileOpenResultMessage(JobFrontendRef theJobFrontend,
                int bfd,
                int ffd,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.bfd = bfd;
            this.ffd = ffd;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.outputFileOpenResult(theJobFrontend, bfd, ffd, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(bfd);
            out.writeInt(ffd);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            bfd = in.readInt();
            ffd = in.readInt();
            exc = (IOException) in.readObject();
        }
    }

    /**
     * Class OutputFileWriteResultMessage provides the Job Backend "output file
     * write result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class OutputFileWriteResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = -4734876024127370851L;

        private int ffd;
        private IOException exc;

        public OutputFileWriteResultMessage() {
        }

        public OutputFileWriteResultMessage(JobFrontendRef theJobFrontend,
                int ffd,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.ffd = ffd;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.outputFileWriteResult(theJobFrontend, ffd, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            ffd = in.readInt();
            exc = (IOException) in.readObject();
        }
    }

    /**
     * Class OutputFileFlushResultMessage provides the Job Backend "output file
     * flush result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class OutputFileFlushResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = 8921581871627030306L;

        private int ffd;
        private IOException exc;

        public OutputFileFlushResultMessage() {
        }

        public OutputFileFlushResultMessage(JobFrontendRef theJobFrontend,
                int ffd,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.ffd = ffd;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.outputFileFlushResult(theJobFrontend, ffd, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            ffd = in.readInt();
            exc = (IOException) in.readObject();
        }
    }

    /**
     * Class OutputFileCloseResultMessage provides the Job Backend "output file
     * close result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class OutputFileCloseResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = 5976163600237430235L;

        private int ffd;
        private IOException exc;

        public OutputFileCloseResultMessage() {
        }

        public OutputFileCloseResultMessage(JobFrontendRef theJobFrontend,
                int ffd,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.ffd = ffd;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.outputFileCloseResult(theJobFrontend, ffd, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            ffd = in.readInt();
            exc = (IOException) in.readObject();
        }
    }

    /**
     * Class InputFileOpenResultMessage provides the Job Backend "input file
     * open result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class InputFileOpenResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = -1082499193559062581L;

        private int bfd;
        private int ffd;
        private IOException exc;

        public InputFileOpenResultMessage() {
        }

        public InputFileOpenResultMessage(JobFrontendRef theJobFrontend,
                int bfd,
                int ffd,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.bfd = bfd;
            this.ffd = ffd;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.inputFileOpenResult(theJobFrontend, bfd, ffd, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(bfd);
            out.writeInt(ffd);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            bfd = in.readInt();
            ffd = in.readInt();
            exc = (IOException) in.readObject();
        }
    }

    /**
     * Class InputFileReadResultMessage provides the Job Backend "input file
     * read result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 20-Nov-2006
     */
    private static class InputFileReadResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = 4542261695089387333L;

        private int ffd;
        private int len;
        private IOException exc;

        public InputFileReadResultMessage() {
        }

        public InputFileReadResultMessage(JobFrontendRef theJobFrontend,
                int ffd,
                int len,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.ffd = ffd;
            this.len = len;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.inputFileReadResult(theJobFrontend, ffd, null, len, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeInt(len);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            ffd = in.readInt();
            len = in.readInt();
            exc = (IOException) in.readObject();
        }
    }

    /**
     * Class InputFileSkipResultMessage provides the Job Backend "input file
     * skip result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class InputFileSkipResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = 5050948958179612039L;

        private int ffd;
        private long len;
        private IOException exc;

        public InputFileSkipResultMessage() {
        }

        public InputFileSkipResultMessage(JobFrontendRef theJobFrontend,
                int ffd,
                long len,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.ffd = ffd;
            this.len = len;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.inputFileSkipResult(theJobFrontend, ffd, len, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeLong(len);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            ffd = in.readInt();
            len = in.readLong();
            exc = (IOException) in.readObject();
        }
    }

    /**
     * Class InputFileCloseResultMessage provides the Job Backend "input file
     * close result" message in the PJ cluster middleware.
     *
     * @author Alan Kaminsky
     * @version 05-Nov-2006
     */
    private static class InputFileCloseResultMessage
            extends JobBackendMessage {

        private static final long serialVersionUID = -5645040374361899565L;

        private int ffd;
        private IOException exc;

        public InputFileCloseResultMessage() {
        }

        public InputFileCloseResultMessage(JobFrontendRef theJobFrontend,
                int ffd,
                IOException exc) {
            super(Message.FROM_JOB_FRONTEND);
            this.ffd = ffd;
            this.exc = exc;
        }

        public void invoke(JobBackendRef theJobBackend,
                JobFrontendRef theJobFrontend)
                throws IOException {
            theJobBackend.inputFileCloseResult(theJobFrontend, ffd, exc);
        }

        public void writeExternal(ObjectOutput out)
                throws IOException {
            out.writeInt(ffd);
            out.writeObject(exc);
        }

        public void readExternal(ObjectInput in)
                throws IOException, ClassNotFoundException {
            ffd = in.readInt();
            exc = (IOException) in.readObject();
        }
    }

}
