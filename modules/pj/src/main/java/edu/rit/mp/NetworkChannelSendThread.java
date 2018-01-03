//******************************************************************************
//
// File:    NetworkChannelSendThread.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.NetworkChannelSendThread
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
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;
import java.util.concurrent.LinkedBlockingQueue;

/**
 * Class NetworkChannelSendThread provides a thread for sending outgoing
 * messages for a {@linkplain NetworkChannel}.
 *
 * @author Alan Kaminsky
 * @version 12-Oct-2008
 */
class NetworkChannelSendThread
        extends Thread {

// Hidden data members.
    // Enclosing network channel.
    private NetworkChannel myNetworkChannel;

    // Underlying socket channel.
    private SocketChannel mySocketChannel;

    // Queue of outgoing I/O requests.
    private LinkedBlockingQueue<IORequest> myOutgoingQueue;

    // Byte buffer.
    private ByteBuffer myByteBuffer;

// Hidden constructors.
    /**
     * Construct a new network channel send thread.
     *
     * @param theNetworkChannel Enclosing network channel.
     * @param theSocketChannel Underlying socket channel.
     * @param theOutgoingQueue Queue of outgoing I/O requests.
     */
    NetworkChannelSendThread(NetworkChannel theNetworkChannel,
            SocketChannel theSocketChannel,
            LinkedBlockingQueue<IORequest> theOutgoingQueue) {
        myNetworkChannel = theNetworkChannel;
        mySocketChannel = theSocketChannel;
        myOutgoingQueue = theOutgoingQueue;
        myByteBuffer = ByteBuffer.allocateDirect(Constants.BUFFER_SIZE);
        setDaemon(true);
        start();
    }

// Exported operations.
    /**
     * Run this network channel send thread.
     */
    public void run() {
        IORequest iorequest = null;
        Buf buf = null;

        // Processing loop.
        sendloop:
        for (;;) {
            try {
                // Wait for an I/O request to show up in the outgoing queue.
                iorequest = myOutgoingQueue.take();
                buf = iorequest.myBuf;

                // Message preprocessing.
                buf.preSend();
                myByteBuffer.clear();
                int i = 0;
                int msglength = buf.myMessageLength;

                // Write message header.
                myByteBuffer.putInt(Constants.MAGIC_NUMBER);
                myByteBuffer.putInt(iorequest.myTagLb);
                myByteBuffer.put(buf.myMessageType);
                myByteBuffer.putInt(msglength);

                // Repeatedly transfer items from source buffer to byte buffer,
                // then from byte buffer to socket channel.
                while (i < msglength) {
                    i += buf.sendItems(i, myByteBuffer);
                    myByteBuffer.flip();
                    mySocketChannel.write(myByteBuffer);
                    myByteBuffer.compact();
                }
                myByteBuffer.flip();
                while (myByteBuffer.hasRemaining()) {
                    mySocketChannel.write(myByteBuffer);
                    myByteBuffer.compact();
                    myByteBuffer.flip();
                }

                // Message postprocessing.
                buf.postSend();

                // Report success of current I/O request.
                iorequest.reportSuccess();

                iorequest = null;
                buf = null;
            } catch (IOException exc) {
                // Report failure of current I/O request.
                if (iorequest != null) {
                    iorequest.reportFailure(exc);
                }
            } catch (InterruptedException exc) {
                ChannelClosedException exc2
                        = new ChannelClosedException("Channel closed");

                // Report failure of current I/O request.
                if (iorequest != null) {
                    iorequest.reportFailure(exc2);
                }

                // Report failure of any pending I/O requests.
                while ((iorequest = myOutgoingQueue.poll()) != null) {
                    iorequest.reportFailure(exc2);
                }

                // Terminate this thread.
                break sendloop;
            } catch (RuntimeException exc) {
                // Report failure of current I/O request.
                if (iorequest != null) {
                    iorequest.reportFailure(exc);
                }
            } catch (Error exc) {
                // Report failure of current I/O request.
                if (iorequest != null) {
                    iorequest.reportFailure(exc);
                }
            }
        }

        // This thread is terminating. Enable garbage collection of data
        // members.
        myNetworkChannel = null;
        mySocketChannel = null;
        myOutgoingQueue = null;
        myByteBuffer = null;
    }

}
