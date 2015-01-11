//******************************************************************************
//
// File:    CommRequest.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.CommRequest
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
package edu.rit.pj;

import edu.rit.mp.IORequest;

import java.io.IOException;

/**
 * Class CommRequest provides an object for doing a non-blocking message passing
 * operation in a PJ cluster parallel program. A non-blocking message passing
 * operation in class {@linkplain Comm} immediately returns a CommRequest object
 * while performing the message passing operation in a separate thread. This
 * allows the calling thread to do other work while the message passing
 * operation is in progress. When the calling thread needs to wait for the
 * message passing operation to finish, the calling thread calls the CommRequest
 * object's <TT>waitForFinish()</TT> method.
 *
 * @author Alan Kaminsky
 * @version 18-Sep-2007
 */
public class CommRequest {

// Hidden data members.
    IORequest mySendRequest;
    IORequest myRecvRequest;

// Exported constructors.
    /**
     * Construct a new CommRequest object.
     */
    public CommRequest() {
    }

// Exported operations.
    /**
     * Determine if the message passing operation associated with this
     * CommRequest object has finished.
     * <P>
     * <I>Note:</I> If the <TT>isFinished()</TT> method is called on a
     * newly-created CommRequest object that has not been passed to or returned
     * from a communicator's non-blocking send or receive method, the
     * <TT>isFinished()</TT> method immediately returns true.
     *
     * @return False if the message passing operation has not finished, true if
     * the message passing operation has finished successfully.
     *
     * @exception IOException Thrown if the message passing operation has
     * finished and an I/O error occurred.
     */
    public boolean isFinished()
            throws IOException {
        return (mySendRequest == null || mySendRequest.isFinished())
                && (myRecvRequest == null || myRecvRequest.isFinished());
    }

    /**
     * Wait for the message passing operation associated with this CommRequest
     * object to finish. If the message passing operation involved a receive, a
     * {@linkplain CommStatus} object is returned giving the results of the
     * receive, otherwise null is returned.
     * <P>
     * For a receive operation, the returned status object gives the actual rank
     * of the process that sent the message, the actual message tag that was
     * received, and the actual number of data items in the message. If the
     * actual number of data items in the message is less than the length of the
     * receive buffer, nothing is stored into the extra data items at the end of
     * the buffer. If the actual number of data items in the message is greater
     * than the length of the receive buffer, the extra data items at the end of
     * the message are discarded.
     * <P>
     * <I>Note:</I> If the <TT>waitForFinish()</TT> method is called on a
     * newly-created CommRequest object that has not been passed to or returned
     * from a communicator's non-blocking send or receive method, the
     * <TT>waitForFinish()</TT> method immediately returns null.
     *
     * @return Status object for a receive operation, otherwise null.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus waitForFinish()
            throws IOException {
        if (mySendRequest != null) {
            mySendRequest.waitForFinish();
        }
        if (myRecvRequest != null) {
            edu.rit.mp.Status status = myRecvRequest.waitForFinish();
            return new CommStatus(Comm.getFarRank(status.channel),
                    status.tag,
                    status.length);
        } else {
            return null;
        }
    }

}
