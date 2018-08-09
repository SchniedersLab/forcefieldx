//******************************************************************************
//
// File:    Proxy.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.Proxy
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

import java.io.IOException;

import edu.rit.mp.Buf;
import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;
import edu.rit.mp.ObjectBuf;
import edu.rit.mp.Status;

/**
 * Class Proxy is the abstract base class for a proxy object for sending
 * messages to a PJ process.
 *
 * @author Alan Kaminsky
 * @version 20-Nov-2006
 */
public abstract class Proxy {

// Hidden data members.
    private ChannelGroup myChannelGroup;
    private Channel myChannel;

// Exported constructors.
    /**
     * Construct a new proxy. The proxy will use the given channel in the given
     * channel group to send messages to the far end process.
     *
     * @param theChannelGroup Channel group.
     * @param theChannelGroup Channel group.
     * @param theChannel Channel.
     */
    public Proxy(ChannelGroup theChannelGroup,
            Channel theChannel) {
        myChannelGroup = theChannelGroup;
        myChannel = theChannel;
    }

// Exported operations.
    /**
     * Send the given message to this proxy's far end process.
     *
     * @param theMessage Message.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void send(Message theMessage)
            throws IOException {
        myChannelGroup.send(myChannel,
                theMessage.getTag(),
                ObjectBuf.buffer(theMessage));
    }

    /**
     * Send a message with the given tag and items to this proxy's far end
     * process.
     *
     * @param theTag Message tag.
     * @param theSrc Item source buffer.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void send(int theTag,
            Buf theSrc)
            throws IOException {
        myChannelGroup.send(myChannel, theTag, theSrc);
    }

    /**
     * Receive a message from this proxy's far end process. A message will be
     * received from this proxy's channel in this proxy's channel group. The
     * message must have a tag of 0. The message items are stored in the given
     * item destination buffer.
     * <P>
     * The <TT>receive()</TT> method does not return until the message has been
     * fully received.
     *
     * @param theDst Item destination buffer.
     * @return Status object giving the outcome of the message reception.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public Status receive(Buf theDst)
            throws IOException {
        return myChannelGroup.receive(myChannel, theDst);
    }

    /**
     * Receive a message with the given tag from this proxy's far end process. A
     * message will be received from this proxy's channel in this proxy's
     * channel group. If <TT>theTag</TT> is null, a message will be received
     * with any tag. The message items are stored in the given item destination
     * buffer.
     * <P>
     * The <TT>receive()</TT> method does not return until the message has been
     * fully received.
     *
     * @param theTag Message tag, or null to receive any tag.
     * @param theDst Item destination buffer.
     * @return Status object giving the outcome of the message reception.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public Status receive(Integer theTag,
            Buf theDst)
            throws IOException {
        return myChannelGroup.receive(myChannel, theTag, theDst);
    }

    /**
     * Close communication with this proxy's far end process.
     */
    public void close() {
        myChannel.close();
    }

}
