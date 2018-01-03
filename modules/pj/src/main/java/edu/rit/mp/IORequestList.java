//******************************************************************************
//
// File:    IORequestList.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.IORequestList
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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
import java.io.PrintStream;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Class IORequestList provides a list of pending I/O requests in the Message
 * Protocol (MP).
 *
 * @author Alan Kaminsky
 * @version 26-Oct-2006
 */
class IORequestList {

// Hidden data members.
    LinkedList<IORequest> myList = new LinkedList<IORequest>();

// Exported constructors.
    /**
     * Construct a new I/O request list.
     */
    public IORequestList() {
    }

// Exported operations.
    /**
     * Add the given I/O request to this list. This list maintains the I/O
     * requests in FIFO order.
     *
     * @param theIORequest I/O request.
     */
    public synchronized void add(IORequest theIORequest) {
        myList.add(theIORequest);
        notifyAll();
    }

    /**
     * Remove an I/O request from this list that matches the given I/O request.
     * The first matching I/O request, if any, is returned.
     *
     * @param theIORequest I/O request.
     * @return Matching I/O request removed from this list, or null if there was
     * no match.
     */
    public synchronized IORequest removeMatch(IORequest theIORequest) {
        Iterator<IORequest> iter = myList.iterator();
        while (iter.hasNext()) {
            IORequest iorequest = iter.next();
            if (iorequest.match(theIORequest)) {
                iter.remove();
                return iorequest;
            }
        }
        return null;
    }

    /**
     * Remove an I/O request from this list that matches the given channel,
     * message tag, and message type. The first matching I/O request, if any, is
     * returned.
     *
     * @param channel Channel.
     * @param tag Message tag.
     * @param type Message type.
     * @return Matching I/O request removed from this list, or null if there was
     * no match.
     */
    public synchronized IORequest removeMatch(Channel channel,
            Integer tag,
            byte type) {
        Iterator<IORequest> iter = myList.iterator();
        while (iter.hasNext()) {
            IORequest iorequest = iter.next();
            if (iorequest.match(channel, tag, type)) {
                iter.remove();
                return iorequest;
            }
        }
        return null;
    }

    /**
     * Wait for an I/O request that matches the given I/O request to be added to
     * this list, then remove it. The first matching I/O request is returned.
     *
     * @param theIORequest I/O request.
     * @return Matching I/O request removed from this list.
     * @exception InterruptedException Thrown if the calling thread was
     * interrupted while blocked in this method.
     * @throws java.lang.InterruptedException if any.
     */
    public synchronized IORequest waitForMatch(IORequest theIORequest)
            throws InterruptedException {
        IORequest iorequest;
        while ((iorequest = removeMatch(theIORequest)) == null) {
//System.out.println ("IORequest.waitForMatch("+theIORequest+") waiting, myList="+this);
            wait();
        }
        return iorequest;
    }

    /**
     * Wait for an I/O request that matches the given channel, message tag, and
     * message type to be added to this list, then remove it. The first matching
     * I/O request is returned.
     *
     * @return Matching I/O request removed from this list.
     * @exception InterruptedException Thrown if the calling thread was
     * interrupted while blocked in this method.
     * @param channel a {@link edu.rit.mp.Channel} object.
     * @param tag a {@link java.lang.Integer} object.
     * @param type a byte.
     * @throws java.lang.InterruptedException if any.
     */
    public synchronized IORequest waitForMatch(Channel channel,
            Integer tag,
            byte type)
            throws InterruptedException {
        IORequest iorequest;
        while ((iorequest = removeMatch(channel, tag, type)) == null) {
//System.out.println ("IORequest.waitForMatch("+channel+","+tag+","+type+") waiting, myList="+this);
            wait();
        }
        return iorequest;
    }

    /**
     * Report that every I/O request in this list failed with an I/O exception.
     *
     * @param theIOException I/O exception.
     */
    public synchronized void reportFailure(IOException theIOException) {
        IORequest iorequest = null;
        while (!myList.isEmpty()) {
            iorequest = myList.remove(0);
            iorequest.reportFailure(theIOException);
        }
    }

    /**
     * Dump the state of this I/O request list on the given print stream. For
     * debugging.
     *
     * @param out Print stream.
     * @param prefix String to print at the beginning of each line.
     */
    public synchronized void dump(PrintStream out,
            String prefix) {
        out.println(prefix + getClass().getName() + "@" + Integer.toHexString(System.identityHashCode(this)));
        out.println(prefix + "\t" + myList.size() + " entries");
        for (IORequest r : myList) {
            out.println(prefix + "\t" + r);
        }
    }

}
