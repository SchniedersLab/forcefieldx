//******************************************************************************
//
// File:    IORequest.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.IORequest
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
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.mp;

import java.io.IOException;
import java.io.InterruptedIOException;

import java.lang.reflect.Constructor;

/**
 * Class IORequest encapsulates the state of a message being sent or received in
 * the Message Protocol (MP).
 * <P>
 * Class {@linkplain ChannelGroup}'s non-blocking <TT>sendNoWait()</TT> method
 * includes an IORequest argument. This allows the caller to initiate the send
 * operation and continue processing while the channel group sends the message
 * in a separate thread. To wait for the message to be sent, the caller must
 * call the IORequest object's <TT>waitForFinish()</TT> method.
 * <P>
 * Class {@linkplain ChannelGroup}'s non-blocking <TT>receiveNoWait()</TT>
 * method includes an IORequest argument. This allows the caller to initiate the
 * receive operation and continue processing while the channel group receives
 * the message in a separate thread. To wait for the message to be received, the
 * caller must call the IORequest object's <TT>waitForFinish()</TT> method,
 * which returns a {@linkplain Status} object giving the results of the receive
 * operation.
 *
 * @author Alan Kaminsky
 * @version 12-Oct-2008
 */
public class IORequest {

// Hidden data members.
    /**
     * Channel on which to send or receive message. Null means "any channel"
     * (only valid for receiving).
     */
    protected Channel myChannel;

    /**
     * Message tag range lower bound. For sending, myTagLb gives the outgoing
     * message tag, and myTagUb is not used. For receiving, an incoming message
     * tag in the range myTagLb..myTabUb inclusive will match this I/O request.
     */
    protected int myTagLb;

    /**
     * Message tag range upper bound.
     */
    protected int myTagUb;

    /**
     * Source or destination of message items.
     */
    protected Buf myBuf;

    /**
     * Status of a received message.
     */
    protected Status myStatus;

    // Exception that occurred, or null if none.
    IOException myIOException;
    RuntimeException myRuntimeException;
    Error myError;

    // State: PENDING = request still in progress; SUCCEEDED = request finished
    // successfully; FAILED = request finished unsuccessfully and myIOException,
    // myRuntimeException, or myError contains the exception object to throw.
    int myState = PENDING;
    static final int PENDING = 0;
    static final int SUCCEEDED = 1;
    static final int FAILED = 2;

// Exported constructors.
    /**
     * Construct a new I/O request object.
     */
    public IORequest() {
    }

    /**
     * Initialize this I/O request object.
     *
     * @param theChannel Channel on which to send or receive message. Null
     * denotes "any channel."
     * @param theTagLb Message tag range lower bound.
     * @param theTagUb Message tag range upper bound.
     * @param theBuf Source or destination of message items.
     */
    void initialize(Channel theChannel,
            int theTagLb,
            int theTagUb,
            Buf theBuf) {
        myChannel = theChannel;
        myTagLb = theTagLb;
        myTagUb = theTagUb;
        myBuf = theBuf;
        myStatus = null;
        myIOException = null;
        myRuntimeException = null;
        myError = null;
        myState = PENDING;
    }

// Exported operations.
    /**
     * Determine if this I/O request has finished.
     *
     * @return False if this I/O request has not finished, true if this I/O
     * request has finished successfully.
     * @exception IOException Thrown if this I/O request has finished and an I/O
     * error occurred.
     * @throws java.io.IOException if any.
     */
    public synchronized boolean isFinished()
            throws IOException {
        if (myState == PENDING) {
            return false;
        }
        rethrow("IORequest: Exception during send or receive");
        return true;
    }

    /**
     * Wait until the send or receive operation corresponding to this I/O
     * request has finished. For a receive operation, a {@linkplain Status}
     * object containing the results of the receive operation is returned; for a
     * send operation, null is returned.
     *
     * @return Receive status for a receive operation, or null for a send
     * operation.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public synchronized Status waitForFinish()
            throws IOException {
        try {
            while (myState == PENDING) {
                wait();
            }
            rethrow("IORequest: Exception during send or receive");
            return myStatus;
        } catch (InterruptedException exc) {
            IOException exc2
                    = new InterruptedIOException("IORequest: waitForFinish() interrupted");
            exc2.initCause(exc);
            throw exc2;
        }
    }

    /**
     * Returns a string version of this I/O request.
     *
     * @return String version.
     */
    public String toString() {
        StringBuilder b = new StringBuilder();
        b.append("IORequest(myChannel=");
        b.append(myChannel);
        b.append(",myTagLb=");
        b.append(myTagLb);
        b.append(",myTagUb=");
        b.append(myTagUb);
        b.append(",myBuf=");
        b.append(myBuf);
        b.append(",myState=");
        b.append(myState);
        b.append(")");
        return b.toString();
    }

// Hidden operations.
    /**
     * Determine if this I/O request matches the given IORequest.
     *
     * @param that I/O request to match.
     *
     * @return True if this I/O request matches <TT>theIORequest</TT>, false
     * otherwise.
     */
    boolean match(IORequest that) {
        return (this.myChannel == null
                || that.myChannel == null
                || this.myChannel == that.myChannel)
                && (this.myTagLb <= that.myTagLb)
                && (that.myTagLb <= this.myTagUb)
                && (this.myBuf.myMessageType == that.myBuf.myMessageType);
    }

    /**
     * Determine if this I/O request matches the given channel, message tag, and
     * message type.
     *
     * @param channel Channel.
     * @param tag Message tag.
     * @param type Message type.
     *
     * @return True if this I/O request matches the given information, false
     * otherwise.
     */
    boolean match(Channel channel,
            int tag,
            byte type) {
        return (this.myChannel == null
                || channel == null
                || this.myChannel == channel)
                && (this.myTagLb <= tag)
                && (tag <= this.myTagUb)
                && (this.myBuf.myMessageType == type);
    }

    /**
     * Report that this I/O request succeeded.
     */
    protected synchronized void reportSuccess() {
        myState = SUCCEEDED;
        notifyAll();
    }

    /**
     * Report that this I/O request failed with an I/O exception.
     *
     * @param theIOException I/O exception.
     */
    protected synchronized void reportFailure(IOException theIOException) {
        myIOException = theIOException;
        myState = FAILED;
        notifyAll();
    }

    /**
     * Report that this I/O request failed with a runtime exception.
     *
     * @param theRuntimeException Runtime exception.
     */
    synchronized void reportFailure(RuntimeException theRuntimeException) {
        myRuntimeException = theRuntimeException;
        myState = FAILED;
        notifyAll();
    }

    /**
     * Report that this I/O request failed with an error.
     *
     * @param theError Error.
     */
    synchronized void reportFailure(Error theError) {
        myError = theError;
        myState = FAILED;
        notifyAll();
    }

    /**
     * Rethrow the exception reported to this I/O request if any.
     *
     * @param msg Detail message for rethrown exception.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private void rethrow(String msg)
            throws IOException {
        if (myIOException != null) {
            rethrowIOException(msg);
        } else if (myRuntimeException != null) {
            rethrowRuntimeException(msg);
        } else if (myError != null) {
            rethrowError(msg);
        }
    }

    /**
     * Rethrow the I/O exception reported to this I/O request.
     *
     * @param msg Detail message for rethrown exception.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private void rethrowIOException(String msg)
            throws IOException {
        IOException exc2 = null;
        try {
            Class<? extends IOException> excClass
                    = myIOException.getClass();
            Constructor<? extends IOException> excConstructor
                    = excClass.getConstructor(String.class);
            exc2 = excConstructor.newInstance(msg);
        } catch (Throwable exc) {
            exc2 = new IOException(msg);
        }
        exc2.initCause(myIOException);
        throw exc2;
    }

    /**
     * Rethrow the runtime exception reported to this I/O request.
     *
     * @param msg Detail message for rethrown exception.
     */
    private void rethrowRuntimeException(String msg) {
        RuntimeException exc2 = null;
        try {
            Class<? extends RuntimeException> excClass
                    = myRuntimeException.getClass();
            Constructor<? extends RuntimeException> excConstructor
                    = excClass.getConstructor(String.class);
            exc2 = excConstructor.newInstance(msg);
        } catch (Throwable exc) {
            exc2 = new RuntimeException(msg);
        }
        exc2.initCause(myRuntimeException);
        throw exc2;
    }

    /**
     * Rethrow the error reported to this I/O request.
     *
     * @param msg Detail message for rethrown exception.
     */
    private void rethrowError(String msg) {
        Error exc2 = null;
        try {
            Class<? extends Error> excClass
                    = myError.getClass();
            Constructor<? extends Error> excConstructor
                    = excClass.getConstructor(String.class);
            exc2 = excConstructor.newInstance(msg);
        } catch (Throwable exc) {
            exc2 = new Error(msg);
        }
        exc2.initCause(myError);
        throw exc2;
    }

}
