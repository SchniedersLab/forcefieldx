//******************************************************************************
//
// File:    Channel.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.Channel
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
package edu.rit.mp;

import java.io.IOException;
import java.io.PrintStream;
import java.net.InetSocketAddress;

/**
 * Class Channel provides a channel for sending and receiving messages in the
 * Message Protocol (MP).
 * <P>
 * A channel object is not constructed directly. Rather, a channel group object
 * is used to connect a channel to another computer, accept a channel connection
 * from another computer, or obtain a "loopback" channel within the same
 * computer. See class {@linkplain ChannelGroup} for further information.
 * <P>
 * An information object may be associated with a channel. An application using
 * MP can use the information object to hold additional data about the channel.
 *
 * @author Alan Kaminsky
 * @version 13-May-2008
 */
public abstract class Channel {

// Hidden data members.
    ChannelGroup myChannelGroup;
    IORequestList myIORequestList;
    Object myInfo;

    int myWriteState = WRITE_OPEN;
    static final int WRITE_OPEN = 0;
    static final int WRITE_CLOSED = 1;

    int myReadState = READ_OPEN;
    static final int READ_OPEN = 0;
    static final int READ_CLOSED = 1;

// Hidden constructors.
    /**
     * Construct a new channel.
     *
     * @param theChannelGroup Enclosing channel group.
     */
    Channel(ChannelGroup theChannelGroup) {
        myChannelGroup = theChannelGroup;
        myIORequestList = theChannelGroup.myIORequestList;
    }

// Exported operations.
    /**
     * Obtain the channel group that created this channel.
     *
     * @return Channel group.
     */
    public ChannelGroup getChannelGroup() {
        return myChannelGroup;
    }

    /**
     * Obtain the channel group ID of this channel's near end channel group.
     *
     * @return Near end channel group ID.
     */
    public abstract int nearEndChannelGroupId();

    /**
     * Obtain the channel group ID of this channel's far end channel group.
     *
     * @return Far end channel group ID.
     */
    public abstract int farEndChannelGroupId();

    /**
     * Obtain this channel's near end address. This is the host and port of the
     * near end of this channel's connection.
     *
     * @return Near end address.
     */
    public abstract InetSocketAddress nearEndAddress();

    /**
     * Obtain this channel's far end address. This is the host and port of the
     * far end of this channel's connection.
     *
     * @return Far end address.
     */
    public abstract InetSocketAddress farEndAddress();

    /**
     * Obtain this channel's information object.
     *
     * @return Information object, or null if none.
     */
    public Object info() {
        return myInfo;
    }

    /**
     * Set this channel's information object.
     *
     * @param theInfo Information object, or null if none.
     */
    public void info(Object theInfo) {
        myInfo = theInfo;
    }

    /**
     * Close this channel. Any pending send requests will fail with a
     * {@linkplain ChannelClosedException}.
     */
    public void close() {
        // To avoid deadlock, synchronize first on the channel group, then on
        // this channel.
        synchronized (myChannelGroup) {
            synchronized (this) {
                myChannelGroup.removeChannel(this);
                myWriteState = WRITE_CLOSED;
                myReadState = READ_CLOSED;
                subclassClose();
            }
        }
    }

    /**
     * Returns a string version of this channel.
     *
     * @return a {@link java.lang.String} object.
     */
    public String toString() {
        return getClass().getName()
                + "(near=" + nearEndAddress()
                + ",far=" + farEndAddress() + ")";
    }

    /**
     * Dump the state of this channel on the given print stream. For debugging.
     *
     * @param out Print stream.
     * @param prefix String to print at the beginning of each line.
     */
    public void dump(PrintStream out,
            String prefix) {
        out.println(prefix + getClass().getName() + "@" + Integer.toHexString(System.identityHashCode(this)));
        out.println(prefix + "myChannelGroup = " + myChannelGroup);
        out.println(prefix + "myIORequestList = " + myIORequestList);
        out.println(prefix + "myInfo = " + myInfo);
        out.println(prefix + "myWriteState = " + (myWriteState == WRITE_OPEN ? "WRITE_OPEN" : "WRITE_CLOSED"));
        out.println(prefix + "myReadState = " + (myReadState == READ_OPEN ? "READ_OPEN" : "READ_CLOSED"));
        out.println(prefix + "nearEndChannelGroupId() = " + nearEndChannelGroupId());
        out.println(prefix + "farEndChannelGroupId() = " + farEndChannelGroupId());
        out.println(prefix + "nearEndAddress() = " + nearEndAddress());
        out.println(prefix + "farEndAddress() = " + farEndAddress());
    }

// Hidden operations.
    /**
     * Start sending and receiving messages via this channel.
     */
    void start() {
    }

    /**
     * Send a message via this channel. The I/O request object must be newly
     * constructed with the message tag and source buffer fields filled in. This
     * method is allowed to return immediately and let the message be sent in a
     * separate thread. The calling thread should use the I/O request object to
     * wait for the message send to complete.
     *
     * @param theIORequest I/O request object.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    abstract void send(IORequest theIORequest)
            throws IOException;

    /**
     * Perform additional close actions in a subclass.
     */
    void subclassClose() {
    }

}
