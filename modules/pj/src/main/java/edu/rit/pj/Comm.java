//******************************************************************************
//
// File:    Comm.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.Comm
//
// This Java source file is copyright (C) 2009 by Alan Kaminsky. All rights
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

import java.io.IOException;
import java.io.InterruptedIOException;
import java.io.PrintStream;
import java.net.InetSocketAddress;
import java.util.LinkedList;

import edu.rit.mp.Buf;
import edu.rit.mp.Channel;
import edu.rit.mp.ChannelGroup;
import edu.rit.mp.ConnectListener;
import edu.rit.mp.IORequest;
import edu.rit.mp.IntegerBuf;
import edu.rit.mp.ObjectBuf;
import edu.rit.mp.Status;
import edu.rit.pj.cluster.CommPattern;
import edu.rit.pj.cluster.JobBackend;
import edu.rit.pj.cluster.JobFrontend;
import edu.rit.pj.cluster.JobSchedulerException;
import edu.rit.pj.reduction.IntegerOp;
import edu.rit.pj.reduction.Op;
import edu.rit.util.Range;

/**
 * Class Comm provides a communicator for a PJ cluster parallel program. Class
 * Comm provides a method to initialize the PJ message passing middleware and
 * run the parallel program on multiple processors of a cluster parallel
 * computer. Class Comm also provides methods for passing messages between the
 * processes of the parallel program.
 * <P>
 * <HR>
 * <P>
 * <B>BASIC CONCEPTS</B>
 * <P>
 * A <B>cluster parallel computer</B> typically consists of a <B>frontend
 * processor</B> and a number of <B>backend processors</B> connected via a
 * dedicated high-speed network. A user logs into the frontend processor and
 * runs a PJ program there. The PJ message passing middleware causes copies of
 * the PJ program to run in a number of separate processes, each process on a
 * different backend processor. The backend processes run the PJ program, using
 * the PJ middleware to send messages amongst themselves. The PJ middleware
 * redirects the backend processes' standard output and standard error streams
 * to the frontend process. The frontend process does not actually execute the
 * PJ program, but merely displays all the backend processes' standard output
 * and standard error streams on the frontend process's own standard output and
 * standard error.
 * <P>
 * For the PJ message passing middleware to work, certain server processes must
 * be running. See package {@linkplain edu.rit.pj.cluster} for further
 * information.
 * <P>
 * To initialize the PJ message passing middleware, the program must first call
 * the static <TT>Comm.init()</TT> method, passing in the command line
 * arguments.
 * <P>
 * A <B>communicator</B> is associated with a group of backend processes. The
 * communicator's <B>size</B> is the number of processes in the group. Each
 * process in the communicator has a different <B>rank</B> in the range 0 ..
 * <I>size</I>-1. A process may obtain the size and rank by calling the
 * communicator's <TT>size()</TT> and <TT>rank()</TT> methods.
 * <P>
 * There is one predefined communicator, the <B>world communicator,</B>
 * consisting of all the backend processes in the parallel program. A process
 * may obtain a reference to the world communicator by calling the static
 * <TT>Comm.world()</TT> method. Typically, the first few lines in a PJ cluster
 * parallel program look like this:
 * <PRE>
 *     public class AParallelProgram
 *         {
 *         public static void main
 *             (String[] args)
 *             throws Exception
 *             {
 *             Comm.init (args);
 *             Comm world = Comm.world();
 *             int size = world.size();
 *             int rank = world.rank();
 *             . . .</PRE>
 * <P>
 * The number of processes in the parallel program -- that is, the size of the
 * world communicator -- is specified by the <TT>"pj.np"</TT> property, which
 * must be an integer greater than or equal to 1. You can specify the number of
 * processes on the Java command line like this:
 * <PRE>
 *     java -Dpj.np=4 . . .</PRE>
 * <P>
 * If the <TT>"pj.np"</TT> property is not specified, the default is 1.
 * <P>
 * The PJ program will run on the specified number of backend processors as
 * described above. To determine which backend processors to use, the PJ program
 * interacts with a <B>Job Scheduler</B> server process running on the frontend
 * processor. When the PJ program starts and calls the <TT>Comm.init()</TT>
 * method, the middleware first prints the job number on the standard error. The
 * middleware then waits until the required number of backend processors are
 * ready to run a job. As each backend processor becomes ready, the middleware
 * prints on the standard error the name of each backend processor assigned to
 * the job. Once all are ready, the PJ program starts running on those backend
 * processors, and all further output comes from the PJ program. Since each PJ
 * program interacts with the Job Scheduler, the Job Scheduler can ensure that
 * each backend processor is running a backend process for only one job at a
 * time.
 * <P>
 * Depending on the system load, your PJ program may have to wait in the Job
 * Scheduler's queue for a while until enough backend processors become ready.
 * If you get tired of waiting, you can kill your PJ program (e.g., by typing
 * CTRL-C), which will remove your PJ program from the Job Scheduler's queue.
 * <P>
 * The Job Scheduler has a web interface that lets you examine the cluster
 * status. Just point your web browser at this URL:
 * <P>
 * <TT>&nbsp;&nbsp;&nbsp;&nbsp;http://&lt;hostname&gt;:8080/</TT>
 * <P>
 * where <TT>&lt;hostname&gt;</TT> is replaced by the host name of the frontend
 * processor. The default port for the cluster status web interface is port
 * 8080. The Job Scheduler can be configured to use a different port. For
 * further information, see package {@linkplain edu.rit.pj.cluster}.
 * <P>
 * If the PJ program is executed on a host where no Job Scheduler is running,
 * the PJ program will run in <I>one</I> process on that host (i.e., the machine
 * you're logged into), rather than on the backend processors. The message
 * passing methods in class Comm will still work, though. This option can be
 * useful for debugging a PJ program's logic on a non-parallel machine before
 * running the PJ program on a cluster.
 * <P>
 * <HR>
 * <P>
 * <B>MESSAGE PASSING</B>
 * <P>
 * PJ provides two categories of communication, <B>point-to-point
 * communication</B> and <B>collective communication.</B> The following methods
 * of class Comm are used for point-to-point communication:
 * <UL>
 * <LI><TT>send()</TT>
 * <LI><TT>receive()</TT>
 * <LI><TT>sendReceive()</TT>
 * <LI><TT>floodSend()</TT>
 * <LI><TT>floodReceive()</TT>
 * </UL>
 * The following methods are used for collective communication:
 * <UL>
 * <LI><TT>broadcast()</TT>
 * <LI><TT>scatter()</TT>
 * <LI><TT>gather()</TT>
 * <LI><TT>allGather()</TT>
 * <LI><TT>reduce()</TT>
 * <LI><TT>allReduce()</TT>
 * <LI><TT>allToAll()</TT>
 * <LI><TT>scan()</TT>
 * <LI><TT>exclusiveScan()</TT>
 * <LI><TT>barrier()</TT>
 * </UL>
 * These methods are described further in the sections below.
 * <P>
 * In addition, you can create a new communicator consisting of all, or a subset
 * of, the processes in an existing communicator. Message passing in the new
 * communicator is completely independent of message passing in the original
 * communicator. The following method creates a new communicator:
 * <UL>
 * <LI><TT>createComm()</TT>
 * </UL>
 * <P>
 * <HR>
 * <P>
 * <B>POINT-TO-POINT COMMUNICATION</B>
 * <P>
 * One process in a PJ cluster parallel program, the <B>source process</B>, may
 * use a communicator to send a message to another process in the program, the
 * <B>destination process</B>. This is called a <B>point-to-point
 * communication</B> because just the two processes are involved (as opposed to
 * a collective communication, which involves all the processes). Five
 * point-to-point communication methods are available in this release: send,
 * receive, send-receive, flood-send, and flood-receive.
 * <P>
 * <B>Send and Receive</B>
 * <P>
 * To do a point-to-point communication, the source process calls the
 * <TT>send()</TT> method on a certain communicator, such as the world
 * communicator. The source process specifies the destination process's rank,
 * the <B>tag</B> for the message, and a <B>buffer</B> containing the data items
 * to be sent (type {@linkplain edu.rit.mp.Buf}). Likewise, the destination
 * process calls the <TT>receive()</TT> method on the same communicator as the
 * source process. The destination process specifies the source process's rank,
 * the message tag which must be the same as in the source process, and the
 * buffer for the data items to be received.
 * <P>
 * A <TT>send()</TT> method call and a <TT>receive()</TT> method call are said
 * to <B>match</B> if (a) the rank passed to the <TT>send()</TT> method equals
 * the rank of the process calling <TT>receive()</TT>, (b) the rank passed to
 * the <TT>receive()</TT> method equals the rank of the process calling
 * <TT>send()</TT>, (c) the item data type in the source buffer is the same as
 * the item data type in the destination buffer, and (d) the send message tag
 * equals the receive message tag. A <TT>receive()</TT> method call will block
 * until a matching <TT>send()</TT> method call occurs. If more than one
 * <TT>send()</TT> method call matches a <TT>receive()</TT> method call, one of
 * the matching <TT>send()</TT> method calls is picked in an unspecified manner.
 * A <TT>send()</TT> method call <I>may</I> block until a matching
 * <TT>receive()</TT> method call occurs due to flow control in the underlying
 * network communication.
 * <P>
 * The message tag can be used to distinguish different kinds of messages. A
 * <TT>receive()</TT> method call will only match a <TT>send()</TT> method call
 * with the same tag. If there is no need to distinguish different kinds of
 * messages, omit the tag (it will default to 0).
 * <P>
 * Once a <TT>send()</TT> method call and a <TT>receive()</TT> method call have
 * been matched together, the actual message data transfer takes place. Each
 * item in the source buffer, starting at index 0 and continuing for the entire
 * length of the source buffer, is written to the message. At the other end,
 * each item in the destination buffer, starting at index 0, is read from the
 * message.
 * <P>
 * The <TT>receive()</TT> method returns a {@linkplain CommStatus} object. The
 * status object gives the actual rank of the process that sent the message, the
 * actual message tag that was received, and the actual number of data items in
 * the message. If the actual number of data items in the message is less than
 * the length of the destination buffer, nothing is stored into the extra data
 * items at the end of the destination buffer. If the actual number of data
 * items in the message is greater than the length of the destination buffer,
 * the extra data items at the end of the message are discarded.
 * <P>
 * The <TT>send()</TT> method does not return until all the message elements
 * have been written from the source buffer. Likewise, the <TT>receive()</TT>
 * method does not return until all the message elements have been read into the
 * destination buffer. However, you cannot assume that because the
 * <TT>send()</TT> method has returned, the matching <TT>receive()</TT> method
 * has also returned. Because of buffering in the underlying network
 * communication, not all the destination items might have been received even
 * though all the source items have been sent.
 * <P>
 * The destination process, instead of specifying a particular source process,
 * can declare that it will receive a message from any source process by
 * specifying null for the source process rank in the <TT>receive()</TT> method
 * call. This is called a <B>wildcard source</B>. In this case the
 * <TT>receive()</TT> method call's returned status object will indicate the
 * actual source process that sent the message.
 * <P>
 * The destination process, instead of specifying a particular message tag, can
 * declare that it will receive a message with any tag by specifying null for
 * the tag in the <TT>receive()</TT> method call. This is called a <B>wildcard
 * tag</B>. Alternatively, the destination process can specify a range of
 * message tags, and it will receive a message with any tag in the given range.
 * In these cases the <TT>receive()</TT> method call's returned status object
 * will indicate the actual message tag that was sent.
 * <P>
 * A process can send a message to itself. In this case one thread must call
 * <TT>send()</TT> (specifying the process's own rank as the destination) and a
 * different thread must call <TT>receive()</TT> (specifying the process's own
 * rank as the source), otherwise a deadlock will ensue.
 * <P>
 * <B>Send-Receive</B>
 * <P>
 * By calling the <TT>sendReceive()</TT> method, a process can send a buffer of
 * outgoing message items to a destination process while simultaneously
 * receiving a buffer of incoming message items from a source process. The
 * destination process may be the same as the source process, or different from
 * the source process. The outgoing message items must come from a different
 * place than where the incoming message items will be stored, otherwise the
 * incoming message items may overwrite the outgoing message items before they
 * can be sent. When the <TT>sendReceive()</TT> method returns, the outgoing
 * message items have been fully sent, but they may not yet have been fully
 * received; and the incoming message items have been fully received.
 * <P>
 * With the <TT>sendReceive()</TT> method, a process cannot receive a message
 * from a wildcard source, and a process cannot receive a message with a
 * wildcard tag or a range of tags. The process calling <TT>sendReceive()</TT>
 * must know the rank of the source process and the message tag (if not 0). The
 * <TT>sendReceive()</TT> method does return a status object giving the outcome
 * of the receive half of the send-receive operation, just as the
 * <TT>receive()</TT> method does.
 * <P>
 * A process can send-receive messages with itself. In this case one thread must
 * call <TT>sendReceive()</TT> (specifying the process's own rank as the source
 * and destination) and a different thread must also call <TT>sendReceive()</TT>
 * (specifying the process's own rank as the source and destination), otherwise
 * a deadlock will ensue.
 * <P>
 * <B>Non-Blocking Communication</B>
 * <P>
 * The <TT>send()</TT>, <TT>receive()</TT>, and <TT>sendReceive()</TT> methods
 * each have a non-blocking version. A non-blocking communication method
 * initiates the communication operation and immediately returns, storing the
 * state of the communication operation in a {@linkplain CommRequest} object.
 * The communicator then performs the communication operation in a separate
 * thread. This allows the calling thread to do other work while the
 * communication operation is in progress. To wait for the send and receive
 * operations to finish, call the CommRequest object's <TT>waitForFinish()</TT>
 * method.
 * <P>
 * <B>Flood-Send and Flood-Receive</B>
 * <P>
 * Any process can send a message to all processes in the communicator. This is
 * called "flooding" the message. First, all processes must start a
 * flood-receive operation, either by calling the non-blocking
 * <TT>floodReceive()</TT> method, or by having a separate thread call the
 * blocking <TT>floodReceive()</TT> method. Then, one process (any process) must
 * call the <TT>floodSend()</TT> method. The data items in the flood-send
 * operation's outgoing buffer are copied into the flood-receive operation's
 * incoming buffer in all processes.
 * <P>
 * Message flooding is similar to the "broadcast" collective communication
 * operation (see below). The differences are these: Broadcasting combines
 * sending and receiving in a single operation; flooding uses separate send and
 * receive operations. For broadcasting, every process must know which process
 * is sending the outgoing data items; for flooding, the receiving processes do
 * not need to know which process is sending (any process can send).
 * <P>
 * <HR>
 * <P>
 * <B>COLLECTIVE COMMUNICATION</B>
 * <P>
 * A PJ cluster parallel program may use a communicator to send a message among
 * all the processes in the program at the same time. This is called a
 * <B>collective communication</B> because all the processes in the communicator
 * are involved (as opposed to a point-to-point communication). Ten collective
 * communication methods are available in this release: broadcast, scatter,
 * gather, all-gather, reduce, all-reduce, all-to-all, scan, exclusive-scan, and
 * barrier. Further collective communication methods will be added to class Comm
 * in a later release.
 * <P>
 * <B>Broadcast</B>
 * <P>
 * One process in the communicator, the <B>root</B> process, has a source buffer
 * (type {@linkplain edu.rit.mp.Buf Buf}) filled with data. The other processes
 * in the communicator each have a destination buffer with the same length and
 * the same item data type as the source buffer. Each process calls the
 * communicator's <TT>broadcast()</TT> method. Afterwards, all the destination
 * buffers contain the same data as the source buffer.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |    |    |    |    |    |
 * |  2 |    |    |    |    |    |    |
 * |  3 |    |    |    |    |    |    |
 * |  4 |    |    |    |    |    |    |
 * |  5 |    |    |    |    |    |    |
 * |  6 |    |    |    |    |    |    |
 * |  7 |    |    |    |    |    |    |
 * |  8 |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  1 |    |  1 |    |  1 |
 * |  2 |    |  2 |    |  2 |    |  2 |
 * |  3 |    |  3 |    |  3 |    |  3 |
 * |  4 |    |  4 |    |  4 |    |  4 |
 * |  5 |    |  5 |    |  5 |    |  5 |
 * |  6 |    |  6 |    |  6 |    |  6 |
 * |  7 |    |  7 |    |  7 |    |  7 |
 * |  8 |    |  8 |    |  8 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * <I>Note:</I> Any process can be the root of the broadcast. The above is only
 * one example with process 0 as the root.
 * <P>
 * <B>Scatter</B>
 * <P>
 * One process in the communicator, the root process, has <I>K</I> source
 * buffers (type {@linkplain edu.rit.mp.Buf Buf}) filled with data, where
 * <I>K</I> is the size of the communicator. For example, the source buffers
 * could be different portions of an array. Each process in the communicator
 * (including the root process) has a destination buffer with the same length
 * and the same item data type as the corresponding source buffer. Each process
 * calls the communicator's <TT>scatter()</TT> method. Afterwards, each
 * process's destination buffer contains the same data as the corresponding
 * source buffer in the root process.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+
 * |  1 |
 * |  2 |
 * +----+
 * |  3 |
 * |  4 |
 * +----+
 * |  5 |
 * |  6 |
 * +----+
 * |  7 |
 * |  8 |
 * +----+
 *
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+
 * |  1 |
 * |  2 |
 * +----+
 * |  3 |
 * |  4 |
 * +----+
 * |  5 |
 * |  6 |
 * +----+
 * |  7 |
 * |  8 |
 * +----+
 *
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * In the root process, the destination buffer can be the same as the source
 * buffer:
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+
 * |  1 |
 * |  2 |
 * +----+    +----+
 * |  3 |    |    |
 * |  4 |    |    |
 * +----+    +----+    +----+
 * |  5 |              |    |
 * |  6 |              |    |
 * +----+              +----+    +----+
 * |  7 |                        |    |
 * |  8 |                        |    |
 * +----+                        +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+
 * |  1 |
 * |  2 |
 * +----+    +----+
 * |  3 |    |  3 |
 * |  4 |    |  4 |
 * +----+    +----+    +----+
 * |  5 |              |  5 |
 * |  6 |              |  6 |
 * +----+              +----+    +----+
 * |  7 |                        |  7 |
 * |  8 |                        |  8 |
 * +----+                        +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * <I>Note:</I> Any process can be the root of the scatter. The above is only
 * one example with process 0 as the root.
 * <P>
 * <B>Gather</B>
 * <P>
 * Gather is the opposite of scatter. One process in the communicator, the root
 * process, has <I>K</I> destination buffers (type {@linkplain edu.rit.mp.Buf
 * Buf}), where <I>K</I> is the size of the communicator. For example, the
 * destination buffers could be different portions of an array. Each process in
 * the communicator (including the root process) has a source buffer with the
 * same length and the same item data type as the corresponding destination
 * buffer, filled with data. Each process calls the communicator's
 * <TT>gather()</TT> method. Afterwards, each destination buffer in the root
 * process contains the same data as the corresponding source buffer.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+
 *
 * +----+
 * |    |
 * |    |
 * +----+
 * |    |
 * |    |
 * +----+
 * |    |
 * |    |
 * +----+
 * |    |
 * |    |
 * +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+
 *
 * +----+
 * |  1 |
 * |  2 |
 * +----+
 * |  3 |
 * |  4 |
 * +----+
 * |  5 |
 * |  6 |
 * +----+
 * |  7 |
 * |  8 |
 * +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * In the root process, the destination buffer can be the same as the source
 * buffer:
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+
 * |  1 |
 * |  2 |
 * +----+    +----+
 * |    |    |  3 |
 * |    |    |  4 |
 * +----+    +----+    +----+
 * |    |              |  5 |
 * |    |              |  6 |
 * +----+              +----+    +----+
 * |    |                        |  7 |
 * |    |                        |  8 |
 * +----+                        +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+
 * |  1 |
 * |  2 |
 * +----+    +----+
 * |  3 |    |  3 |
 * |  4 |    |  4 |
 * +----+    +----+    +----+
 * |  5 |              |  5 |
 * |  6 |              |  6 |
 * +----+              +----+    +----+
 * |  7 |                        |  7 |
 * |  8 |                        |  8 |
 * +----+                        +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * <I>Note:</I> Any process can be the root of the gather. The above is only one
 * example with process 0 as the root.
 * <P>
 * <B>All-Gather</B>
 * <P>
 * All-gather is the same as gather, except that every process has an array of
 * destination buffers, and every process receives the results of the gather.
 * Each process in the communicator has a source buffer (type {@linkplain
 * edu.rit.mp.Buf Buf}) filled with data. Each process in the communicator has
 * <I>K</I> destination buffers, where <I>K</I> is the size of the communicator.
 * For example, the destination buffers could be different portions of an array.
 * Each destination buffer has the same length and the same item data type as
 * the corresponding source buffer. Each process calls the communicator's
 * <TT>allGather()</TT> method. Afterwards, each destination buffer in every
 * process contains the same data as the corresponding source buffer.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+
 *
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+
 *
 * +----+    +----+    +----+    +----+
 * |  1 |    |  1 |    |  1 |    |  1 |
 * |  2 |    |  2 |    |  2 |    |  2 |
 * +----+    +----+    +----+    +----+
 * |  3 |    |  3 |    |  3 |    |  3 |
 * |  4 |    |  4 |    |  4 |    |  4 |
 * +----+    +----+    +----+    +----+
 * |  5 |    |  5 |    |  5 |    |  5 |
 * |  6 |    |  6 |    |  6 |    |  6 |
 * +----+    +----+    +----+    +----+
 * |  7 |    |  7 |    |  7 |    |  7 |
 * |  8 |    |  8 |    |  8 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * The destination buffer can be the same as the source buffer in each process:
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |    |    |    |    |    |
 * |  2 |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |  3 |    |    |    |    |
 * |    |    |  4 |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |  5 |    |    |
 * |    |    |    |    |  6 |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |  7 |
 * |    |    |    |    |    |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  1 |    |  1 |    |  1 |
 * |  2 |    |  2 |    |  2 |    |  2 |
 * +----+    +----+    +----+    +----+
 * |  3 |    |  3 |    |  3 |    |  3 |
 * |  4 |    |  4 |    |  4 |    |  4 |
 * +----+    +----+    +----+    +----+
 * |  5 |    |  5 |    |  5 |    |  5 |
 * |  6 |    |  6 |    |  6 |    |  6 |
 * +----+    +----+    +----+    +----+
 * |  7 |    |  7 |    |  7 |    |  7 |
 * |  8 |    |  8 |    |  8 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * <P>
 * <B>Reduce</B>
 * <P>
 * Reduce is like gather, except the buffers' contents are combined together
 * instead of just copied. Each process in the communicator has a buffer (type
 * {@linkplain edu.rit.mp.Buf Buf}) filled with data. Each process calls the
 * communicator's <TT>reduce()</TT> method, specifying some binary operation
 * (type {@linkplain edu.rit.pj.reduction.Op Op}) for combining the data.
 * Afterwards, each element of the buffer in the root process contains the
 * result of combining all the corresponding elements in all the buffers using
 * the specified binary operation. For example, if the operation is addition,
 * each buffer element in the root process ends up being the sum of the
 * corresponding buffer elements in all the processes. In the non-root
 * processes, the buffers' contents may be changed from their original contents.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 * 0 (root)     1         2         3
 * +----+    +----+    +----+    +----+
 * | 16 |    | ?? |    | ?? |    | ?? |
 * | 20 |    | ?? |    | ?? |    | ?? |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * <I>Note:</I> Any process can be the root of the reduction. The above is only
 * one example with process 0 as the root.
 * <P>
 * <B>All-Reduce</B>
 * <P>
 * All-reduce is the same as reduce, except that every process receives the
 * results of the reduction. Each process in the communicator has a buffer (type
 * {@linkplain edu.rit.mp.Buf Buf}) filled with data. Each process calls the
 * communicator's <TT>allReduce()</TT> method, specifying some binary operation
 * (type {@linkplain edu.rit.pj.reduction.Op Op}) for combining the data.
 * Afterwards, each element of the buffer in each process contains the result of
 * combining all the corresponding elements in all the buffers using the
 * specified binary operation. For example, if the operation is addition, each
 * buffer element ends up being the sum of the corresponding buffer elements.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * | 16 |    | 16 |    | 16 |    | 16 |
 * | 20 |    | 20 |    | 20 |    | 20 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * <P>
 * <B>All-to-All</B>
 * <P>
 * Every process in the communicator has <I>K</I> source buffers (type
 * {@linkplain edu.rit.mp.Buf Buf}) filled with data, where <I>K</I> is the size
 * of the communicator. Every process in the communicator also has <I>K</I>
 * destination buffers (type {@linkplain edu.rit.mp.Buf Buf}). The source
 * buffers and the destination buffers must refer to different storage. For
 * example, the source buffers could be portions of an array, and the
 * destination buffers could be portions of a different array. Each process
 * calls the communicator's <TT>allToAll()</TT> method. Afterwards, for each
 * process rank <I>k</I>, 0 &lt;= <I>k</I> &lt;= <I>K</I>-1, and each buffer
 * index <I>i</I>, 0 &lt;= <I>i</I> &lt;= <I>K</I>-1, destination buffer
 * <I>i</I> in process <I>k</I> contains the same data as source buffer <I>k</I>
 * in process <I>i</I>.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  9 |    | 17 |    | 25 |
 * |  2 |    | 10 |    | 18 |    | 26 |
 * +----+    +----+    +----+    +----+
 * |  3 |    | 11 |    | 19 |    | 27 |
 * |  4 |    | 12 |    | 20 |    | 28 |
 * +----+    +----+    +----+    +----+
 * |  5 |    | 13 |    | 21 |    | 29 |
 * |  6 |    | 14 |    | 22 |    | 30 |
 * +----+    +----+    +----+    +----+
 * |  7 |    | 15 |    | 23 |    | 31 |
 * |  8 |    | 16 |    | 24 |    | 32 |
 * +----+    +----+    +----+    +----+
 *
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+
 * |    |    |    |    |    |    |    |
 * |    |    |    |    |    |    |    |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  9 |    | 17 |    | 25 |
 * |  2 |    | 10 |    | 18 |    | 26 |
 * +----+    +----+    +----+    +----+
 * |  3 |    | 11 |    | 19 |    | 27 |
 * |  4 |    | 12 |    | 20 |    | 28 |
 * +----+    +----+    +----+    +----+
 * |  5 |    | 13 |    | 21 |    | 29 |
 * |  6 |    | 14 |    | 22 |    | 30 |
 * +----+    +----+    +----+    +----+
 * |  7 |    | 15 |    | 23 |    | 31 |
 * |  8 |    | 16 |    | 24 |    | 32 |
 * +----+    +----+    +----+    +----+
 *
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+
 * |  9 |    | 11 |    | 13 |    | 15 |
 * | 10 |    | 12 |    | 14 |    | 16 |
 * +----+    +----+    +----+    +----+
 * | 17 |    | 19 |    | 21 |    | 23 |
 * | 18 |    | 20 |    | 22 |    | 24 |
 * +----+    +----+    +----+    +----+
 * | 25 |    | 27 |    | 29 |    | 31 |
 * | 26 |    | 28 |    | 30 |    | 32 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * <P>
 * <B>Scan</B>
 * <P>
 * Each process in the communicator has a buffer (type {@linkplain
 * edu.rit.mp.Buf Buf}) filled with data. Each process calls the communicator's
 * <TT>scan()</TT> method, specifying some binary operation (type {@linkplain
 * edu.rit.pj.reduction.Op Op}) for combining the data. Afterwards, each element
 * of the buffer in a particular process contains the result of combining all
 * the corresponding elements in its own and all lower-ranked processes' buffers
 * using the specified binary operation. For example, if the operation is
 * addition, each buffer element ends up being the sum of its own and all
 * lower-ranked processes' buffer elements.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  4 |    |  9 |    | 16 |
 * |  2 |    |  6 |    | 12 |    | 20 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * The scan operation is also known as "prefix scan" or "inclusive prefix scan"
 * -- "inclusive" because the process's own element is included in the result.
 * <P>
 * <B>Exclusive-Scan</B>
 * <P>
 * The exclusive-scan operation is a variation of the scan operation. Each
 * process in the communicator has a buffer (type {@linkplain edu.rit.mp.Buf
 * Buf}) filled with data. Each process calls the communicator's
 * <TT>exclusiveScan()</TT> method, specifying some binary operation (type
 * {@linkplain edu.rit.pj.reduction.Op Op}) for combining the data, and
 * specifying an initial data value. Afterwards, each element of the buffer in a
 * particular process contains the result of combining all the corresponding
 * elements in all lower-ranked processes' buffers using the specified binary
 * operation, except in process 0 each element of the buffer contains the
 * initial data value. For example, if the operation is addition and the initial
 * data value is 0, each buffer element ends up being the sum of all
 * lower-ranked processes' buffer elements.
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR>
 * <TD ALIGN="left" VALIGN="top">
 * Before:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  1 |    |  3 |    |  5 |    |  7 |
 * |  2 |    |  4 |    |  6 |    |  8 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * <TD WIDTH=50> </TD>
 * <TD ALIGN="left" VALIGN="top">
 * After:
 * <FONT SIZE="-1">
 * <PRE>
 * Process   Process   Process   Process
 *    0         1         2         3
 * +----+    +----+    +----+    +----+
 * |  0 |    |  1 |    |  4 |    |  9 |
 * |  0 |    |  2 |    |  6 |    | 12 |
 * +----+    +----+    +----+    +----+</PRE>
 * </FONT>
 * </TD>
 * </TR>
 * </TABLE>
 * This version of the scan operation is also known as "exclusive prefix scan"
 * -- "exclusive" because the process's own element is excluded from the result.
 * <P>
 * <B>Barrier</B>
 * <P>
 * The barrier operation causes all the processes to synchronize with each
 * other. Each process calls the communicator's <TT>barrier()</TT> method. The
 * calling thread blocks until all processes in the communicator have called the
 * <TT>barrier()</TT> method. Then the calling thread unblocks and returns from
 * the <TT>barrier()</TT> method call.
 *
 * @author Alan Kaminsky
 * @version 21-Jan-2009
 */
public class Comm {

// Hidden data members.
    // Predefined communicators.
    private static Comm theWorldCommunicator;
    private static Comm theFrontendCommunicator;

    // This communicator's size, rank, and host.
    private int mySize;
    private int myRank;
    private String myHost;

    // The largest power of 2 less than or equal to this communicator's size.
    private int mySizePowerOf2;

    // Channel group for message passing in this communicator.
    private ChannelGroup myChannelGroup;

    // Map from rank (array index) to channel group address (array element).
    private InetSocketAddress[] myAddressForRank;

    // Map from rank (array index) to channel for communicating with the process
    // at that rank (array element).
    private Channel[] myChannelForRank;

    // Broadcast trees for flood-send, flood-receive, broadcast, and reduce
    // operations, indexed by root.
    private int[][] myBroadcastTree;

// Hidden constructors.
    /**
     * Construct a new communicator.
     *
     * @param size Communicator's size.
     * @param rank Current process's rank in the communicator.
     * @param host Host name.
     * @param channelgroup Channel group for message passing in this
     * communicator.
     * @param address Map from rank (array index) to channel group address
     * (array element).
     */
    private Comm(int size,
            int rank,
            String host,
            ChannelGroup channelgroup,
            InetSocketAddress[] address) {
        // Record size, rank, channel group.
        mySize = size;
        myRank = rank;
        myHost = host;
        myChannelGroup = channelgroup;

        // Determine the largest power of 2 less than or equal to this
        // communicator's size.
        int p2 = 1;
        while (p2 <= size) {
            p2 <<= 1;
        }
        mySizePowerOf2 = p2 >>> 1;

        // Set channel group ID equal to rank.
        myChannelGroup.setChannelGroupId(rank);

        // Set up connect listener.
        myChannelGroup.setConnectListener(new ConnectListener() {
            public void nearEndConnected(ChannelGroup theChannelGroup,
                    Channel theChannel)
                    throws IOException {
            }

            public void farEndConnected(ChannelGroup theChannelGroup,
                    Channel theChannel)
                    throws IOException {
                doFarEndConnected(theChannel);
            }
        });

        // Record socket address for each process rank.
        myAddressForRank = address;

        // Set up channel for each process rank.
        myChannelForRank = new Channel[size];

        // Populate channel at my own rank with the loopback channel.
        myChannelForRank[myRank] = channelgroup.loopbackChannel();

        // If there's more than one process, start listening for incoming
        // connections.
        if (mySize > 1) {
            myChannelGroup.startListening();
        }
    }

// Exported operations.
    /**
     * Initialize the PJ message passing middleware. Certain Java system
     * properties specify the middleware's behavior; these properties are
     * typically given on the Java command line with the <TT>"-D"</TT> flag. For
     * further information, see class {@linkplain PJProperties}.
     *
     * @param args Command line arguments.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>args</TT> is null.
     * @exception IllegalArgumentException (unchecked exception) Thrown if the
     * value of one of the Java system properties is illegal.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public static void init(String[] args)
            throws IOException {
        // Verify preconditions.
        if (args == null) {
            throw new NullPointerException("Comm.init(): args is null");
        }

        // Get the job backend object.
        JobBackend backend = JobBackend.getJobBackend();

        if (backend == null) {
			// We're running on the frontend processor.

            // Prepare constructor arguments for the Job Frontend object.
            String username = System.getProperty("user.name");
            int Nn = PJProperties.getPjNn();
            int Np = PJProperties.getPjNp();
            int Nt = PJProperties.getPjNt();
            boolean hasFrontendComm = false;

            // Examine the call stack to find the main program class name.
            StackTraceElement[] stack
                    = Thread.currentThread().getStackTrace();
            StackTraceElement bottom = stack[stack.length - 1];
            if (!bottom.getMethodName().equals("main")) {
                throw new IllegalStateException("Comm.init(): Not called from main program");
            }
            String mainClassName = bottom.getClassName();

            // Set up the Job Frontend object.
            JobFrontend frontend = null;
            try {
                frontend
                        = new JobFrontend(username, Nn, Np, Nt, hasFrontendComm, mainClassName,
                                args);

                // We were able to contact the Job Scheduler.
                // Run the job frontend in this process, then exit.
                frontend.run();
                System.exit(0);
            } catch (JobSchedulerException exc) {
                // We were not able to contact the Job Scheduler.
                System.err.println(" No Job Scheduler at "
                        + PJProperties.getPjHost() + ":"
                        + PJProperties.getPjPort()
                        + ", running in this (one) process");

                // Set up world communicator.
                theWorldCommunicator
                        = new Comm(/*size        */1,
                                /*rank        */ 0,
                                /*host        */ "<unknown>",
                                /*channelgroup*/ new ChannelGroup(),
                                /*address     */
                                new InetSocketAddress[]{new InetSocketAddress(0)});
            }
        } else {
			// We're running on a backend processor.

            // Set up world communicator.
            theWorldCommunicator
                    = new Comm(/*size        */backend.getK(),
                            /*rank        */ backend.getRank(),
                            /*host        */ backend.getBackendHost(),
                            /*channelgroup*/ backend.getWorldChannelGroup(),
                            /*address     */ backend.getWorldAddress());
        }
    }

    /**
     * Obtain a reference to the world communicator.
     *
     * @return World communicator.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if
     * <TT>Comm.init()</TT> has not been called. Thrown if <TT>world()</TT> is
     * called in the job frontend process; the world communicator does not exist
     * in the job frontend process.
     */
    public static Comm world() {
        if (theWorldCommunicator != null) {
            return theWorldCommunicator;
        } else if (JobBackend.getJobBackend() != null) {
            throw new IllegalStateException("Comm.world(): Didn't call Comm.init()");
        } else {
            throw new IllegalStateException("Comm.world(): World communicator doesn't exist in job frontend process");
        }
    }

    /**
     * Obtain the number of processes in this communicator.
     *
     * @return Size.
     */
    public int size() {
        return mySize;
    }

    /**
     * Obtain the current process's rank in this communicator.
     *
     * @return Rank.
     */
    public int rank() {
        return myRank;
    }

    /**
     * Obtain the host name of this communicator's backend processor. If this
     * communicator is not running on a cluster backend processor, the host name
     * is <TT>"&lt;unknown&gt;"</TT>.
     *
     * @return Host name.
     */
    public String host() {
        return myHost;
    }

    /**
     * Create a new communicator. <I>Every</I> process in this communicator must
     * call the <TT>createComm()</TT> method. Each process passes true or false
     * for the <TT>participate</TT> argument to state whether the process will
     * participate in the new communicator. At least one process must
     * participate in the new communicator. Messages to set up the new
     * communicator are sent to all processes in this communicator, using a
     * message tag of 0.
     * <P>
     * In processes participating in the new communicator, the new communicator
     * is returned. The participating processes appear in the same order by rank
     * in the new communicator as in this communicator. The process can call the
     * new communicator's <TT>rank()</TT> method to determine the process's rank
     * in the new communicator.
     * <P>
     * In processes not participating in the new communicator, null is returned.
     *
     * @param participate True if this process will participate in the new
     * communicator; false otherwise.
     *
     * @return New communicator if this process will participate in the new
     * communicator; null otherwise.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public Comm createComm(boolean participate)
            throws IOException {
        return createComm(participate, 0);
    }

    /**
     * Create a new communicator. <I>Every</I> process in this communicator must
     * call the <TT>createComm()</TT> method. Each process passes true or false
     * for the <TT>participate</TT> argument to state whether the process will
     * participate in the new communicator. At least one process must
     * participate in the new communicator. Messages to set up the new
     * communicator are sent to all processes in this communicator, using the
     * given message tag.
     * <P>
     * In processes participating in the new communicator, the new communicator
     * is returned. The participating processes appear in the same order by rank
     * in the new communicator as in this communicator. The process can call the
     * new communicator's <TT>rank()</TT> method to determine the process's rank
     * in the new communicator.
     * <P>
     * In processes not participating in the new communicator, null is returned.
     *
     * @param participate True if this process will participate in the new
     * communicator; false otherwise.
     * @param tag Message tag.
     *
     * @return New communicator if this process will participate in the new
     * communicator; null otherwise.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public Comm createComm(boolean participate,
            int tag)
            throws IOException {
        // Set up array of socket addresses for all processes.
        InetSocketAddress[] address = new InetSocketAddress[mySize];
        ObjectBuf<InetSocketAddress>[] addressbuf
                = ObjectBuf.sliceBuffers(address,
                        new Range(0, mySize - 1).subranges(mySize));

        // Create channel group for new communicator, if participating.
        ChannelGroup channelgroup = null;
        InetSocketAddress myaddress = null;
        if (participate) {
            channelgroup
                    = new ChannelGroup(new InetSocketAddress(myChannelGroup.listenAddress().getAddress(), 0));
            myaddress = channelgroup.listenAddress();
            address[myRank] = myaddress;
        }

        // All-gather channel group socket addresses into every process.
        allGather(tag, addressbuf[myRank], addressbuf);

        // Close up gaps in the socket address array if any.
        int off = 0;
        int newsize = 0;
        int newrank = -1;
        for (int i = 0; i < mySize; ++i) {
            if (address[i] == null) {
                ++off;
            } else {
                if (i == myRank) {
                    newrank = i - off;
                }
                address[i - off] = address[i];
                ++newsize;
            }
        }

        // Verify size of new communicator.
        if (newsize == 0) {
            throw new IOException("Comm.createComm(): No processes in communicator");
        }

        // Return new communicator if participating.
        if (participate) {
            return new Comm(newsize, newrank, myHost, channelgroup, address);
        } // Return null if not participating.
        else {
            return null;
        }
    }

    /**
     * Send a message to the process at the given rank in this communicator. The
     * message uses a tag of 0. The message items come from the given buffer. To
     * receive the message, the destination process must call the
     * <TT>receive()</TT> method. When the <TT>send()</TT> method returns, the
     * message has been fully sent, but it may not yet have been fully received.
     * <P>
     * A process can send a message to itself; in this case a different thread
     * must call the <TT>receive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param buffer Buffer of data items to be sent.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void send(int toRank,
            Buf buffer)
            throws IOException {
        send(toRank, 0, buffer);
    }

    /**
     * Send a message to the process at the given rank in this communicator with
     * the given message tag. The message items come from the given buffer. To
     * receive the message, the destination process must call the
     * <TT>receive()</TT> method. When the <TT>send()</TT> method returns, the
     * message has been fully sent, but it may not yet have been fully received.
     * <P>
     * A process can send a message to itself; in this case a different thread
     * must call the <TT>receive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param tag Message tag.
     * @param buffer Buffer of data items to be sent.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void send(int toRank,
            int tag,
            Buf buffer)
            throws IOException {
        myChannelGroup.send(getChannel(toRank), tag, buffer);
    }

    /**
     * Send a message to the process at the given rank in this communicator
     * (non-blocking). A message tag of 0 is used. The message items come from
     * the given buffer. To receive the message, the destination process must
     * call the <TT>receive()</TT> method.
     * <P>
     * The <TT>send()</TT> method initiates the send operation and immediately
     * returns a {@linkplain CommRequest} object. The send operation is
     * performed by a separate thread. To wait for the send operation to finish,
     * call the returned {@linkplain CommRequest} object's
     * <TT>waitForFinish()</TT> method. When that method returns, the message
     * has been fully sent, but it may not yet have been fully received.
     * <P>
     * A process can send a message to itself; in this case a different thread
     * must call the <TT>receive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param buffer Buffer of data items to be sent.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest send(int toRank,
            Buf buffer,
            CommRequest request)
            throws IOException {
        return send(toRank, 0, buffer, request);
    }

    /**
     * Send a message to the process at the given rank in this communicator with
     * the given message tag (non-blocking). The message items come from the
     * given buffer. To receive the message, the destination process must call
     * the <TT>receive()</TT> method.
     * <P>
     * The <TT>send()</TT> method initiates the send operation and immediately
     * returns a {@linkplain CommRequest} object. The send operation is
     * performed by a separate thread. To wait for the send operation to finish,
     * call the returned {@linkplain CommRequest} object's
     * <TT>waitForFinish()</TT> method. When that method returns, the message
     * has been fully sent, but it may not yet have been fully received.
     * <P>
     * A process can send a message to itself; in this case a different thread
     * must call the <TT>receive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param tag Message tag.
     * @param buffer Buffer of data items to be sent.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest send(int toRank,
            int tag,
            Buf buffer,
            CommRequest request)
            throws IOException {
        // Set up CommRequest object.
        CommRequest req = request == null ? new CommRequest() : request;

        // Send message (non-blocking).
        req.mySendRequest = new IORequest();
        req.myRecvRequest = null;
        myChannelGroup.sendNoWait(getChannel(toRank), tag, buffer, req.mySendRequest);

        // Return CommRequest object.
        return req;
    }

    /**
     * Receive a message from the process at the given rank in this
     * communicator. If <TT>rank</TT> is null, a message will be received from
     * any process in this communicator. The message must have a tag of 0. The
     * received message items are stored in the given buffer. To send the
     * message, the source process must call the <TT>send()</TT> method. When
     * the <TT>receive()</TT> method returns, the message has been fully
     * received.
     * <P>
     * A {@linkplain CommStatus} object is returned. The status object gives the
     * actual rank of the process that sent the message, the actual message tag
     * that was received, and the actual number of data items in the message. If
     * the actual number of data items in the message is less than the length of
     * the buffer, nothing is stored into the extra data items at the end of the
     * buffer. If the actual number of data items in the message is greater than
     * the length of the buffer, the extra data items at the end of the message
     * are discarded.
     * <P>
     * A process can receive a message from itself; in this case a different
     * thread must call the <TT>send()</TT> method on this communicator.
     *
     * @param fromRank Source process's rank in this communicator, or null to
     * receive from any process.
     * @param buffer Buffer of data items to be received.
     *
     * @return Status object giving the outcome of the message reception.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>fromRank</TT> is not null and is not in the range 0 ..
     * <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus receive(Integer fromRank,
            Buf buffer)
            throws IOException {
        return receive(fromRank, 0, buffer);
    }

    /**
     * Receive a message from the process at the given rank in this communicator
     * with the given message tag. If <TT>rank</TT> is null, a message will be
     * received from any process in this communicator. The received message
     * items are stored in the given buffer. To send the message, the source
     * process must call the <TT>send()</TT> method. When the <TT>receive()</TT>
     * method returns, the message has been fully received.
     * <P>
     * A {@linkplain CommStatus} object is returned. The status object gives the
     * actual rank of the process that sent the message, the actual message tag
     * that was received, and the actual number of data items in the message. If
     * the actual number of data items in the message is less than the length of
     * the buffer, nothing is stored into the extra data items at the end of the
     * buffer. If the actual number of data items in the message is greater than
     * the length of the buffer, the extra data items at the end of the message
     * are discarded.
     * <P>
     * A process can receive a message from itself; in this case a different
     * thread must call the <TT>send()</TT> method on this communicator.
     *
     * @param fromRank Source process's rank in this communicator, or null to
     * receive from any process.
     * @param tag Message tag.
     * @param buffer Buffer of data items to be received.
     *
     * @return Status object giving the outcome of the message reception.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>fromRank</TT> is not null and is not in the range 0 ..
     * <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus receive(Integer fromRank,
            int tag,
            Buf buffer)
            throws IOException {
        Status status;

        // If source is a wildcard, ensure a channel to every process, then
        // receive from any process.
        if (fromRank == null) {
            for (int src = 0; src < mySize; ++src) {
                ensureChannel(src);
            }
            status = myChannelGroup.receive(null, tag, buffer);
        } // If source is not a wildcard, receive from that process.
        else {
            status
                    = myChannelGroup.receive(getChannel(fromRank), tag, buffer);
        }

        // Return CommStatus object.
        return new CommStatus(getFarRank(status.channel),
                status.tag,
                status.length);
    }

    /**
     * Receive a message from the process at the given rank in this communicator
     * with the given message tag range. If <TT>rank</TT> is null, a message
     * will be received from any process in this communicator. If
     * <TT>tagRange</TT> is null, a message will be received with any tag. If
     * <TT>tagRange</TT> is not null, a message will be received with any tag in
     * the given range. The received message items are stored in the given
     * buffer. To send the message, the source process must call the
     * <TT>send()</TT> method. When the <TT>receive()</TT> method returns, the
     * message has been fully received.
     * <P>
     * A {@linkplain CommStatus} object is returned. The status object gives the
     * actual rank of the process that sent the message, the actual message tag
     * that was received, and the actual number of data items in the message. If
     * the actual number of data items in the message is less than the length of
     * the buffer, nothing is stored into the extra data items at the end of the
     * buffer. If the actual number of data items in the message is greater than
     * the length of the buffer, the extra data items at the end of the message
     * are discarded.
     * <P>
     * A process can receive a message from itself; in this case a different
     * thread must call the <TT>send()</TT> method on this communicator.
     *
     * @param fromRank Source process's rank in this communicator, or null to
     * receive from any process.
     * @param tagRange Message tag range, or null to receive any tag.
     * @param buffer Buffer of data items to be received.
     *
     * @return Status object giving the outcome of the message reception.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>fromRank</TT> is not null and is not in the range 0 ..
     * <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus receive(Integer fromRank,
            Range tagRange,
            Buf buffer)
            throws IOException {
        Status status;

        // If source is a wildcard, ensure a channel to every process, then
        // receive from any process.
        if (fromRank == null) {
            for (int src = 0; src < mySize; ++src) {
                ensureChannel(src);
            }
            status = myChannelGroup.receive(null, tagRange, buffer);
        } // If source is not a wildcard, receive from that process.
        else {
            status = myChannelGroup.receive(getChannel(fromRank), tagRange, buffer);
        }

        // Return CommStatus object.
        return new CommStatus(getFarRank(status.channel),
                status.tag,
                status.length);
    }

    /**
     * Receive a message from the process at the given rank in this communicator
     * (non-blocking). If <TT>rank</TT> is null, a message will be received from
     * any process in this communicator. The message must have a tag of 0. The
     * received message items are stored in the given buffer. To send the
     * message, the source process must call the <TT>send()</TT> method.
     * <P>
     * The <TT>receive()</TT> method initiates the receive operation and
     * immediately returns a {@linkplain CommRequest} object. The receive
     * operation is performed by a separate thread. To wait for the receive
     * operation to finish, call the returned {@linkplain CommRequest} object's
     * <TT>waitForFinish()</TT> method. When that method returns, the incoming
     * message items have been fully received.
     * <P>
     * A process can receive a message from itself; in this case a different
     * thread must call the <TT>send()</TT> method on this communicator.
     *
     * @param fromRank Source process's rank in this communicator, or null to
     * receive from any process.
     * @param buffer Buffer of data items to be received.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>fromRank</TT> is not null and is not in the range 0 ..
     * <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest receive(Integer fromRank,
            Buf buffer,
            CommRequest request)
            throws IOException {
        return receive(fromRank, 0, buffer, request);
    }

    /**
     * Receive a message from the process at the given rank in this communicator
     * with the given message tag (non-blocking). If <TT>rank</TT> is null, a
     * message will be received from any process in this communicator. If
     * <TT>tag</TT> is null, a message will be received with any tag. The
     * received message items are stored in the given buffer. To send the
     * message, the source process must call the <TT>send()</TT> method.
     * <P>
     * The <TT>receive()</TT> method initiates the receive operation and
     * immediately returns a {@linkplain CommRequest} object. The receive
     * operation is performed by a separate thread. To wait for the receive
     * operation to finish, call the returned {@linkplain CommRequest} object's
     * <TT>waitForFinish()</TT> method. When that method returns, the incoming
     * message items have been fully received.
     * <P>
     * A process can receive a message from itself; in this case a different
     * thread must call the <TT>send()</TT> method on this communicator.
     *
     * @param fromRank Source process's rank in this communicator, or null to
     * receive from any process.
     * @param tag Message tag.
     * @param buffer Buffer of data items to be received.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>fromRank</TT> is not null and is not in the range 0 ..
     * <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest receive(Integer fromRank,
            int tag,
            Buf buffer,
            CommRequest request)
            throws IOException {
        // Set up CommRequest object.
        CommRequest req = request == null ? new CommRequest() : request;
        req.mySendRequest = null;
        req.myRecvRequest = new IORequest();

        // If source is a wildcard, ensure a channel to every process, then
        // receive (non-blocking) from any process.
        if (fromRank == null) {
            for (int src = 0; src < mySize; ++src) {
                ensureChannel(src);
            }
            myChannelGroup.receiveNoWait(null, tag, buffer, req.myRecvRequest);
        } // If source is not a wildcard, receive (non-blocking) from that
        // process.
        else {
            myChannelGroup.receiveNoWait(getChannel(fromRank), tag, buffer, req.myRecvRequest);
        }

        // Return CommRequest object.
        return req;
    }

    /**
     * Receive a message from the process at the given rank in this communicator
     * with the given message tag range (non-blocking). If <TT>rank</TT> is
     * null, a message will be received from any process in this communicator.
     * If <TT>tagRange</TT> is null, a message will be received with any tag. If
     * <TT>tagRange</TT> is not null, a message will be received with any tag in
     * the given range. The received message items are stored in the given
     * buffer. To send the message, the source process must call the
     * <TT>send()</TT> method.
     * <P>
     * The <TT>receive()</TT> method initiates the receive operation and
     * immediately returns a {@linkplain CommRequest} object. The receive
     * operation is performed by a separate thread. To wait for the receive
     * operation to finish, call the returned {@linkplain CommRequest} object's
     * <TT>waitForFinish()</TT> method. When that method returns, the incoming
     * message items have been fully received.
     * <P>
     * A process can receive a message from itself; in this case a different
     * thread must call the <TT>send()</TT> method on this communicator.
     *
     * @param fromRank Source process's rank in this communicator, or null to
     * receive from any process.
     * @param tagRange Message tag range, or null to receive any tag.
     * @param buffer Buffer of data items to be received.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>fromRank</TT> is not null and is not in the range 0 ..
     * <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest receive(Integer fromRank,
            Range tagRange,
            Buf buffer,
            CommRequest request)
            throws IOException {
        // Set up CommRequest object.
        CommRequest req = request == null ? new CommRequest() : request;
        req.mySendRequest = null;
        req.myRecvRequest = new IORequest();

        // If source is a wildcard, ensure a channel to every process, then
        // receive (non-blocking) from any process.
        if (fromRank == null) {
            for (int src = 0; src < mySize; ++src) {
                ensureChannel(src);
            }
            myChannelGroup.receiveNoWait(null, tagRange, buffer, req.myRecvRequest);
        } // If source is not a wildcard, receive (non-blocking) from that
        // process.
        else {
            myChannelGroup.receiveNoWait(getChannel(fromRank), tagRange, buffer, req.myRecvRequest);
        }

        // Return CommRequest object.
        return req;
    }

    /**
     * Send a message to the process at the given rank in this communicator, and
     * receive a message from the process at the given rank in this
     * communicator. A message tag of 0 is used. The outgoing message items come
     * from the buffer <TT>sendbuf</TT>. The incoming message items go into the
     * buffer <TT>recvbuf</TT>. The outgoing message items must come from a
     * different place than where the incoming message items will be stored. The
     * destination process (process <TT>toRank</TT>) must call a method to
     * receive this process's outgoing message items. The source process
     * (process <TT>fromRank</TT>) must call a method to send this process's
     * incoming message items. When the <TT>sendReceive()</TT> method returns,
     * the outgoing message items have been fully sent, but they may not yet
     * have been fully received; and the incoming message items have been fully
     * received.
     * <P>
     * A {@linkplain CommStatus} object is returned giving the results of the
     * receive half of the operation. The status object gives the rank of the
     * process that sent the incoming message, the message tag that was
     * received, and the actual number of data items in the message. If the
     * actual number of data items in the message is less than the length of the
     * receive buffer, nothing is stored into the extra data items at the end of
     * the receive buffer. If the actual number of data items in the message is
     * greater than the length of the receive buffer, the extra data items at
     * the end of the message are discarded.
     * <P>
     * A process can send-receive messages with itself; in this case a different
     * thread must call the <TT>sendReceive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param sendBuf Buffer of data items to be sent.
     * @param fromRank Source process's rank in this communicator.
     * @param recvBuf Buffer of data items to be received.
     *
     * @return Status object giving the outcome of the message reception.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> or <TT>fromRank</TT>
     * is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>sendBuf</TT> or <TT>recvBuf</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus sendReceive(int toRank,
            Buf sendBuf,
            int fromRank,
            Buf recvBuf)
            throws IOException {
        return sendReceive(toRank, 0, sendBuf, fromRank, 0, recvBuf);
    }

    /**
     * Send a message to the process at the given rank in this communicator with
     * the given message tag, and receive a message from the process at the
     * given rank in this communicator with the given message tag. The outgoing
     * message items come from the buffer <TT>sendbuf</TT>. The incoming message
     * items go into the buffer <TT>recvbuf</TT>. The outgoing message items
     * must come from a different place than where the incoming message items
     * will be stored. The destination process (process <TT>toRank</TT>) must
     * call a method to receive this process's outgoing message items. The
     * source process (process <TT>fromRank</TT>) must call a method to send
     * this process's incoming message items. When the <TT>sendReceive()</TT>
     * method returns, the outgoing message items have been fully sent, but they
     * may not yet have been fully received; and the incoming message items have
     * been fully received.
     * <P>
     * A {@linkplain CommStatus} object is returned giving the results of the
     * receive half of the operation. The status object gives the rank of the
     * process that sent the incoming message, the message tag that was
     * received, and the actual number of data items in the message. If the
     * actual number of data items in the message is less than the length of the
     * receive buffer, nothing is stored into the extra data items at the end of
     * the receive buffer. If the actual number of data items in the message is
     * greater than the length of the receive buffer, the extra data items at
     * the end of the message are discarded.
     * <P>
     * A process can send-receive messages with itself; in this case a different
     * thread must call the <TT>sendReceive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param sendTag Message tag for outgoing message.
     * @param sendBuf Buffer of data items to be sent.
     * @param fromRank Source process's rank in this communicator.
     * @param recvTag Message tag for incoming message.
     * @param recvBuf Buffer of data items to be received.
     *
     * @return Status object giving the outcome of the message reception.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> or <TT>fromRank</TT>
     * is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>sendBuf</TT> or <TT>recvBuf</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus sendReceive(int toRank,
            int sendTag,
            Buf sendBuf,
            int fromRank,
            int recvTag,
            Buf recvBuf)
            throws IOException {
        // Send the outgoing message (non-blocking).
        IORequest sendRequest = new IORequest();
        myChannelGroup.sendNoWait(getChannel(toRank), sendTag, sendBuf, sendRequest);

        // Receive the outgoing message (non-blocking).
        IORequest recvRequest = new IORequest();
        myChannelGroup.receiveNoWait(getChannel(fromRank), recvTag, recvBuf, recvRequest);

        // Wait for both messages to finish.
        sendRequest.waitForFinish();
        Status status = recvRequest.waitForFinish();

        return new CommStatus(getFarRank(status.channel),
                status.tag,
                status.length);
    }

    /**
     * Send a message to the process at the given rank in this communicator, and
     * receive a message from the process at the given rank in this communicator
     * (non-blocking). A message tag of 0 is used. The outgoing message items
     * come from the buffer <TT>sendbuf</TT>. The incoming message items go into
     * the buffer <TT>recvbuf</TT>. The outgoing message items must come from a
     * different place than where the incoming message items will be stored. The
     * destination process (process <TT>toRank</TT>) must call a method to
     * receive this process's outgoing message items. The source process
     * (process <TT>fromRank</TT>) must call a method to send this process's
     * incoming message items.
     * <P>
     * The <TT>sendReceive()</TT> method initiates the send and receive
     * operations and immediately returns a {@linkplain CommRequest} object. The
     * send and receive operations are performed by a separate thread. To wait
     * for the send and receive operations to finish, call the returned
     * {@linkplain CommRequest} object's <TT>waitForFinish()</TT> method. When
     * that method returns, the outgoing message items have been fully sent, but
     * they may not yet have been fully received; and the incoming message items
     * have been fully received.
     * <P>
     * A process can send-receive messages with itself; in this case a different
     * thread must call the <TT>sendReceive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param sendBuf Buffer of data items to be sent.
     * @param fromRank Source process's rank in this communicator.
     * @param recvBuf Buffer of data items to be received.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> or <TT>fromRank</TT>
     * is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>sendBuf</TT> or <TT>recvBuf</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest sendReceive(int toRank,
            Buf sendBuf,
            int fromRank,
            Buf recvBuf,
            CommRequest request)
            throws IOException {
        return sendReceive(toRank, 0, sendBuf, fromRank, 0, recvBuf, request);
    }

    /**
     * Send a message to the process at the given rank in this communicator with
     * the given message tag, and receive a message from the process at the
     * given rank in this communicator with the given message tag
     * (non-blocking). The outgoing message items come from the buffer
     * <TT>sendbuf</TT>. The incoming message items go into the buffer
     * <TT>recvbuf</TT>. The outgoing message items must come from a different
     * place than where the incoming message items will be stored. The
     * destination process (process <TT>toRank</TT>) must call a method to
     * receive this process's outgoing message items. The source process
     * (process <TT>fromRank</TT>) must call a method to send this process's
     * incoming message items.
     * <P>
     * The <TT>sendReceive()</TT> method initiates the send and receive
     * operations and immediately returns a {@linkplain CommRequest} object. The
     * send and receive operations are performed by a separate thread. To wait
     * for the send and receive operations to finish, call the returned
     * {@linkplain CommRequest} object's <TT>waitForFinish()</TT> method. When
     * that method returns, the outgoing message items have been fully sent, but
     * they may not yet have been fully received; and the incoming message items
     * have been fully received.
     * <P>
     * A process can send-receive messages with itself; in this case a different
     * thread must call the <TT>sendReceive()</TT> method on this communicator.
     *
     * @param toRank Destination process's rank in this communicator.
     * @param sendTag Message tag for outgoing message.
     * @param sendBuf Buffer of data items to be sent.
     * @param fromRank Source process's rank in this communicator.
     * @param recvTag Message tag for incoming message.
     * @param recvBuf Buffer of data items to be received.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>toRank</TT> or <TT>fromRank</TT>
     * is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>sendBuf</TT> or <TT>recvBuf</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest sendReceive(int toRank,
            int sendTag,
            Buf sendBuf,
            int fromRank,
            int recvTag,
            Buf recvBuf,
            CommRequest request)
            throws IOException {
        // Set up CommRequest object.
        CommRequest req = request == null ? new CommRequest() : request;

        // Send the outgoing message (non-blocking).
        req.mySendRequest = new IORequest();
        myChannelGroup.sendNoWait(getChannel(toRank), sendTag, sendBuf, req.mySendRequest);

        // Receive the outgoing message (non-blocking).
        req.myRecvRequest = new IORequest();
        myChannelGroup.receiveNoWait(getChannel(fromRank), recvTag, recvBuf, req.myRecvRequest);

        // Return CommRequest object.
        return req;
    }

    /**
     * Flood-send a message to all processes in this communicator. The message
     * uses a tag of 0. The message items come from the given buffer. To receive
     * the message, every process (including the sending process) must call the
     * <TT>floodReceive()</TT> method. When the <TT>floodSend()</TT> method
     * returns, the message has been fully sent, but it may not yet have been
     * fully received in all processes.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param buffer Buffer of data items to be sent.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void floodSend(Buf buffer)
            throws IOException {
        floodSend(0, buffer, null).waitForFinish();
    }

    /**
     * Flood-send a message to all processes in this communicator with the given
     * message tag. The message items come from the given buffer. To receive the
     * message, every process (including the sending process) must call the
     * <TT>floodReceive()</TT> method. When the <TT>floodSend()</TT> method
     * returns, the message has been fully sent, but it may not yet have been
     * fully received in all processes.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param tag Message tag.
     * @param buffer Buffer of data items to be sent.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void floodSend(int tag,
            Buf buffer)
            throws IOException {
        floodSend(tag, buffer, null).waitForFinish();
    }

    /**
     * Flood-send a message to all processes in this communicator
     * (non-blocking). A message tag of 0 is used. The message items come from
     * the given buffer. To receive the message, every process (including the
     * sending process) must call the <TT>floodReceive()</TT> method.
     * <P>
     * The <TT>floodSend()</TT> method initiates the flood-send operation and
     * immediately returns a {@linkplain CommRequest} object. The flood-send
     * operation is performed by a separate thread. To wait for the flood-send
     * operation to finish, call the returned {@linkplain CommRequest} object's
     * <TT>waitForFinish()</TT> method. When that method returns, the message
     * has been fully sent, but it may not yet have been fully received in all
     * processes.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param buffer Buffer of data items to be sent.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest floodSend(Buf buffer,
            CommRequest request)
            throws IOException {
        return floodSend(0, buffer, request);
    }

    /**
     * Flood-send a message to all processes in this communicator with the given
     * message tag (non-blocking). The message items come from the given buffer.
     * To receive the message, every process (including the sending process)
     * must call the <TT>floodReceive()</TT> method.
     * <P>
     * The <TT>floodSend()</TT> method initiates the flood-send operation and
     * immediately returns a {@linkplain CommRequest} object. The flood-send
     * operation is performed by a separate thread. To wait for the flood-send
     * operation to finish, call the returned {@linkplain CommRequest} object's
     * <TT>waitForFinish()</TT> method. When that method returns, the message
     * has been fully sent, but it may not yet have been fully received in all
     * processes.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param tag Message tag.
     * @param buffer Buffer of data items to be sent.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest floodSend(int tag,
            Buf buffer,
            CommRequest request)
            throws IOException {
        // Set up CommRequest object.
        CommRequest req = request == null ? new CommRequest() : request;
        req.mySendRequest = new IORequest();
        req.myRecvRequest = null;

        // Send data to process 0. Process 0's flood-receive I/O request object
        // will forward the data to all the processes.
        myChannelGroup.sendNoWait(getChannel(0), tag, buffer, req.mySendRequest);

        // Return CommRequest object.
        return req;
    }

    /**
     * Flood-receive a message from any process in this communicator. The
     * message must have a tag of 0. The received message items are stored in
     * the given buffer. To send the message, the source process must call the
     * <TT>floodSend()</TT> method. When the <TT>floodReceive()</TT> method
     * returns, the message has been fully received.
     * <P>
     * A {@linkplain CommStatus} object is returned. The status object gives the
     * actual rank of the process that sent the message, the actual message tag
     * that was received, and the actual number of data items in the message.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param buffer Buffer of data items to be received.
     *
     * @return Status object giving the outcome of the message reception.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus floodReceive(Buf buffer)
            throws IOException {
        return floodReceive(0, buffer, null).waitForFinish();
    }

    /**
     * Flood-receive a message from any process in this communicator with the
     * given message tag. If <TT>tag</TT> is null, a message will be received
     * with any tag. The received message items are stored in the given buffer.
     * To send the message, the source process must call the
     * <TT>floodSend()</TT> method. When the <TT>floodReceive()</TT> method
     * returns, the message has been fully received.
     * <P>
     * A {@linkplain CommStatus} object is returned. The status object gives the
     * actual rank of the process that sent the message, the actual message tag
     * that was received, and the actual number of data items in the message.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param tag Message tag, or null to receive any tag.
     * @param buffer Buffer of data items to be received.
     *
     * @return Status object giving the outcome of the message reception.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommStatus floodReceive(Integer tag,
            Buf buffer)
            throws IOException {
        return floodReceive(tag, buffer, null).waitForFinish();
    }

    /**
     * Flood-receive a message from any process in this communicator
     * (non-blocking). A message tag of 0 is used. If <TT>tag</TT> is null, a
     * message will be received with any tag. The received message items are
     * stored in the given buffer. To send the message, the source process must
     * call the <TT>floodSend()</TT> method.
     * <P>
     * The <TT>floodReceive()</TT> method initiates the flood-receive operation
     * and immediately returns a {@linkplain CommRequest} object. The
     * flood-receive operation is performed by a separate thread. To wait for
     * the flood-receive operation to finish, call the returned {@linkplain
     * CommRequest} object's <TT>waitForFinish()</TT> method. When that method
     * returns, the incoming message items have been fully received.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param buffer Buffer of data items to be received.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest floodReceive(Buf buffer,
            CommRequest request)
            throws IOException {
        return floodReceive(0, buffer, request);
    }

    /**
     * Flood-receive a message from any process in this communicator with the
     * given message tag (non-blocking). If <TT>tag</TT> is null, a message will
     * be received with any tag. The received message items are stored in the
     * given buffer. To send the message, the source process must call the
     * <TT>floodSend()</TT> method.
     * <P>
     * The <TT>floodReceive()</TT> method initiates the flood-receive operation
     * and immediately returns a {@linkplain CommRequest} object. The
     * flood-receive operation is performed by a separate thread. To wait for
     * the flood-receive operation to finish, call the returned {@linkplain
     * CommRequest} object's <TT>waitForFinish()</TT> method. When that method
     * returns, the incoming message items have been fully received.
     * <P>
     * <I>Note:</I> The length of the incoming buffer in the
     * <TT>floodReceive()</TT> method call must be the same as the length of the
     * outgoing buffer in the <TT>floodSend()</TT> method call.
     *
     * @param tag Message tag, or null to receive any tag.
     * @param buffer Buffer of data items to be received.
     * @param request CommRequest object to use to wait for the operation to
     * finish; in this case <TT>request</TT> is returned. If
     * <TT>request</TT> is null, a new CommRequest object is created and
     * returned.
     *
     * @return CommRequest object to use to wait for the operation to finish.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public CommRequest floodReceive(Integer tag,
            Buf buffer,
            CommRequest request)
            throws IOException {
        // Get broadcast tree for root=0.
        int[] tree = getBroadcastTree(0);

        // Set up CommRequest object with a special I/O request object that
        // forwards the message down the broadcast tree.
        CommRequest req = request == null ? new CommRequest() : request;
        req.mySendRequest = null;
        req.myRecvRequest = new FloodReceiveIORequest(tree);

        // In process 0, ensure a channel to every process, then receive
        // (non-blocking) a message from any process.
        if (myRank == 0) {
            for (int src = 0; src < mySize; ++src) {
                ensureChannel(src);
            }
            myChannelGroup.receiveNoWait(null, tag, buffer, req.myRecvRequest);
        } // In other processes, ensure a channel to the child processes in the
        // broadcast tree, then receive (non-blocking) a message from the parent
        // process in the broadcast tree.
        else {
            for (int i = 1; i < tree.length; ++i) {
                ensureChannel(tree[i]);
            }
            myChannelGroup.receiveNoWait(getChannel(tree[0]), tag, buffer, req.myRecvRequest);
        }

        // Return CommRequest object.
        return req;
    }

    /**
     * Class FloodReceiveIORequest overrides the methods of class IORequest with
     * additional processing to forward the message when a message is received.
     */
    private class FloodReceiveIORequest
            extends IORequest {

        // Broadcast tree.
        private int[] tree;

        // List of zero or more additional I/O requests to forward copies of the
        // received message.
        private LinkedList<IORequest> myForwardedIORequests
                = new LinkedList<IORequest>();

        /**
         * Construct a new I/O request object.
         *
         * @param tree Broadcast tree.
         */
        public FloodReceiveIORequest(int[] tree) {
            super();
            this.tree = tree;
        }

        /**
         * Determine if this I/O request has finished.
         *
         * @return False if this I/O request has not finished, true if this I/O
         * request has finished successfully.
         *
         * @exception IOException Thrown if this I/O request has finished and an
         * I/O error occurred.
         */
        public synchronized boolean isFinished()
                throws IOException {
            if (!super.isFinished()) {
                return false;
            }
            for (IORequest req : myForwardedIORequests) {
                if (!req.isFinished()) {
                    return false;
                }
            }
            return true;
        }

        /**
         * Wait until the send or receive operation corresponding to this I/O
         * request has finished. For a receive operation, a {@linkplain Status}
         * object containing the results of the receive operation is returned;
         * for a send operation, null is returned.
         *
         * @return Receive status for a receive operation, or null for a send
         * operation.
         *
         * @exception IOException Thrown if an I/O error occurred.
         */
        public synchronized Status waitForFinish()
                throws IOException {
            Status status = super.waitForFinish();
            for (IORequest req : myForwardedIORequests) {
                req.waitForFinish();
            }
            return status;
        }

        /**
         * Report that this I/O request succeeded.
         */
        protected synchronized void reportSuccess() {
            try {
                super.reportSuccess();

                // Get message tag.
                int msgtag = myStatus.tag;

                // Flood the message to every child process in the broadcast
                // tree.
                int n = tree.length;
                for (int i = 1; i < n; ++i) {
                    IORequest req = new IORequest();
                    myForwardedIORequests.add(req);
                    myChannelGroup.sendNoWait(getChannel(tree[i]), msgtag, myBuf, req);
                }
            } catch (IOException exc) {
                reportFailure(exc);
            }
        }
    }

    /**
     * Broadcast a message to all processes in this communicator. The broadcast
     * uses a message tag of 0. All processes must call <TT>broadcast()</TT>
     * with the same value for <TT>root</TT> and with a buffer of the same
     * length and the same item data type.
     * <P>
     * The root process (the process whose rank in this communicator is
     * <TT>root</TT>) sends the message items. The message items come from the
     * given buffer. When the <TT>broadcast()</TT> method returns, the message
     * has been fully sent, but it may not yet have been fully received by all
     * processes.
     * <P>
     * Each non-root process receives the message items. The message items are
     * stored in the given buffer. When the <TT>broadcast()</TT> method returns,
     * the message has been fully received.
     *
     * @param root Root process's rank in this communicator.
     * @param buffer Buffer of data items to be sent (root process) or received
     * (non-root processes).
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void broadcast(int root,
            Buf buffer)
            throws IOException {
        broadcast(root, 0, buffer);
    }

    /**
     * Broadcast a message to all processes in this communicator using the given
     * message tag. All processes must call <TT>broadcast()</TT> with the same
     * values for <TT>root</TT> and <TT>tag</TT> and with a buffer of the same
     * length and the same item data type.
     * <P>
     * The root process (the process whose rank in this communicator is
     * <TT>root</TT>) sends the message items. The message items come from the
     * given buffer. When the <TT>broadcast()</TT> method returns, the message
     * has been fully sent, but it may not yet have been fully received by all
     * processes.
     * <P>
     * Each non-root process receives the message items. The message items are
     * stored in the given buffer. When the <TT>broadcast()</TT> method returns,
     * the message has been fully received.
     *
     * @param root Root process's rank in this communicator.
     * @param tag Message tag.
     * @param buffer Buffer of data items to be sent (root process) or received
     * (non-root processes).
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buffer</TT> is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void broadcast(int root,
            int tag,
            Buf buffer)
            throws IOException {
        // Verify preconditions.
        if (0 > root || root >= mySize) {
            throw new IndexOutOfBoundsException("Comm.broadcast(): root = " + root + " out of bounds");
        }

        // Early return if only one process.
        if (mySize == 1) {
            return;
        }

        // A broadcast is done as a series of point-to-point messages. The
        // messages are organized into rounds. The number of rounds is
        // ceil(log_2(mySize)). In each round, processes send messages to other
        // processes in parallel. Here is the message pattern for a communicator
        // with 8 processes doing a broadcast from root process 0:
        //
        //     Process
        //     0     1     2     3     4     5     6     7
        //     |     |     |     |     |     |     |     |
        //     |---->|     |     |     |     |     |     |  Round 1
        //     |     |     |     |     |     |     |     |
        //     |---------->|     |     |     |     |     |  Round 2
        //     |     |---------->|     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |---------------------->|     |     |     |  Round 3
        //     |     |---------------------->|     |     |
        //     |     |     |---------------------->|     |
        //     |     |     |     |---------------------->|
        //     |     |     |     |     |     |     |     |
        //
        // If a process other than process 0 is the root, the message pattern is
        // the same, except the process ranks are circularly rotated.
        // Get array of process ranks in the broadcast tree.
        int[] broadcasttree = getBroadcastTree(root);
        int n = broadcasttree.length;

        // Receive data from parent if any (blocking).
        int parent = broadcasttree[0];
        if (parent != -1) {
            myChannelGroup.receive(getChannel(parent), tag, buffer);
        }

        // Send data to children if any (non-blocking).
        IORequest[] iorequest = new IORequest[n];
        for (int i = 1; i < n; ++i) {
            int child = broadcasttree[i];
            iorequest[i] = new IORequest();
            myChannelGroup.sendNoWait(getChannel(child), tag, buffer, iorequest[i]);
        }

        // Wait for sends to finish if any.
        for (int i = 1; i < n; ++i) {
            iorequest[i].waitForFinish();
        }
    }

    /**
     * Scatter messages to all processes in this communicator. The scatter uses
     * a message tag of 0. All processes must call <TT>scatter()</TT>
     * with the same value for <TT>root</TT>.
     * <P>
     * The root process (the process whose rank in this communicator is
     * <TT>root</TT>) sends the message items. The message items sent to process
     * <I>i</I> come from the source buffer at index <I>i</I> in the given array
     * of source buffers. When the <TT>scatter()</TT> method returns, the
     * messages have been fully sent, but they may not yet have been fully
     * received by all processes.
     * <P>
     * Each process, including the root process, receives the message items. The
     * message items are stored in the given destination buffer. This must have
     * the same length and the same item data type as the corresponding source
     * buffer. When the <TT>scatter()</TT> method returns, the message has been
     * fully received.
     * <P>
     * In the non-root processes, the source buffer array is ignored and may be
     * null.
     *
     * @param root Root process's rank in this communicator.
     * @param srcarray Array of source buffers to be sent by the root process.
     * Ignored in the non-root processes.
     * @param dst Destination buffer to be received.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1. Thrown if
     * <TT>srcarray</TT>'s length does not equal the size of this communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>srcarray</TT> or any element thereof is null. Thrown if <TT>dst</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void scatter(int root,
            Buf[] srcarray,
            Buf dst)
            throws IOException {
        scatter(root, 0, srcarray, dst);
    }

    /**
     * Scatter messages to all processes in this communicator using the given
     * message tag. All processes must call <TT>scatter()</TT> with the same
     * values for <TT>root</TT> and <TT>tag</TT>.
     * <P>
     * The root process (the process whose rank in this communicator is
     * <TT>root</TT>) sends the message items. The message items sent to process
     * <I>i</I> come from the source buffer at index <I>i</I> in the given array
     * of source buffers. When the <TT>scatter()</TT> method returns, the
     * messages have been fully sent, but they may not yet have been fully
     * received by all processes.
     * <P>
     * Each process, including the root process, receives the message items. The
     * message items are stored in the given destination buffer. This must have
     * the same length and the same item data type as the corresponding source
     * buffer. When the <TT>scatter()</TT> method returns, the message has been
     * fully received.
     * <P>
     * In the non-root processes, the source buffer array is ignored and may be
     * null.
     *
     * @param root Root process's rank in this communicator.
     * @param tag Message tag.
     * @param srcarray Array of source buffers to be sent by the root process.
     * Ignored in the non-root processes.
     * @param dst Destination buffer to be received.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1. Thrown if
     * <TT>srcarray</TT>'s length does not equal the size of this communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>srcarray</TT> or any element thereof is null. Thrown if <TT>dst</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void scatter(int root,
            int tag,
            Buf[] srcarray,
            Buf dst)
            throws IOException {
        // Verify preconditions.
        if (0 > root || root >= mySize) {
            throw new IndexOutOfBoundsException("Comm.scatter(): root = " + root + " out of bounds");
        }

        // A scatter is done as a series of point-to-point messages. The root
        // process sends a separate message to every other process. Here is the
        // message pattern for a communicator with 8 processes scattering from
        // root process 0:
        //
        //     Process
        //     0     1     2     3     4     5     6     7
        //     |     |     |     |     |     |     |     |
        //     |---->|     |     |     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |---------->|     |     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |---------------->|     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |---------------------->|     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |---------------------------->|     |     |
        //     |     |     |     |     |     |     |     |
        //     |---------------------------------->|     |
        //     |     |     |     |     |     |     |     |
        //     |---------------------------------------->|
        //     |     |     |     |     |     |     |     |
        // Root process sends all messages.
        if (myRank == root) {
            // Array of IORequest objects for non-blocking sends.
            IORequest[] iorequest = new IORequest[mySize];

            // Initiate sends to lower-ranked processes.
            for (int rank = 0; rank < myRank; ++rank) {
                iorequest[rank] = new IORequest();
                myChannelGroup.sendNoWait(getChannel(rank), tag, srcarray[rank], iorequest[rank]);
            }

            // Initiate sends to higher-ranked processes.
            for (int rank = myRank + 1; rank < mySize; ++rank) {
                iorequest[rank] = new IORequest();
                myChannelGroup.sendNoWait(getChannel(rank), tag, srcarray[rank], iorequest[rank]);
            }

            // Copy to itself.
            dst.copy(srcarray[myRank]);

            // Wait for completion of sends to lower-ranked processes.
            for (int rank = 0; rank < myRank; ++rank) {
                iorequest[rank].waitForFinish();
            }

            // Wait for completion of sends to higher-ranked processes.
            for (int rank = myRank + 1; rank < mySize; ++rank) {
                iorequest[rank].waitForFinish();
            }
        } // Non-root process receives one message.
        else {
            myChannelGroup.receive(getChannel(root), tag, dst);
        }
    }

    /**
     * Gather messages from all processes in this communicator. The gather uses
     * a message tag of 0. All processes must call <TT>gather()</TT> with the
     * same value for <TT>root</TT>.
     * <P>
     * The root process (the process whose rank in this communicator is
     * <TT>root</TT>) receives the message items. The message items received
     * from process <I>i</I> are stored in the destination buffer at index
     * <I>i</I> in the given array of destination buffers. When the
     * <TT>gather()</TT> method returns, all the messages have been fully
     * received.
     * <P>
     * Each process, including the root process, sends the message items. The
     * message items come from the given source buffer. This must have the same
     * length and the same item data type as the corresponding destination
     * buffer. When the <TT>gather()</TT> method returns, the message has been
     * fully sent, but it may not yet have been fully received by the root
     * process.
     * <P>
     * In the non-root processes, the destination buffer array is ignored and
     * may be null.
     *
     * @param root Root process's rank in this communicator.
     * @param src Source buffer to be sent.
     * @param dstarray Array of destination buffers to be received by the root
     * process. Ignored in the non-root processes.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1. Thrown if
     * <TT>dstarray</TT>'s length does not equal the size of this communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>dstarray</TT> or any element thereof is null. Thrown if <TT>src</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void gather(int root,
            Buf src,
            Buf[] dstarray)
            throws IOException {
        gather(root, 0, src, dstarray);
    }

    /**
     * Gather messages from all processes in this communicator using the given
     * message tag. All processes must call <TT>gather()</TT> with the same
     * values for <TT>root</TT> and <TT>tag</TT>.
     * <P>
     * The root process (the process whose rank in this communicator is
     * <TT>root</TT>) receives the message items. The message items received
     * from process <I>i</I> are stored in the destination buffer at index
     * <I>i</I> in the given array of destination buffers. When the
     * <TT>gather()</TT> method returns, all the messages have been fully
     * received.
     * <P>
     * Each process, including the root process, sends the message items. The
     * message items come from the given source buffer. This must have the same
     * length and the same item data type as the corresponding destination
     * buffer. When the <TT>gather()</TT> method returns, the message has been
     * fully sent, but it may not yet have been fully received by the root
     * process.
     * <P>
     * In the non-root processes, the destination buffer array is ignored and
     * may be null.
     *
     * @param root Root process's rank in this communicator.
     * @param tag Message tag.
     * @param src Source buffer to be sent.
     * @param dstarray Array of destination buffers to be received by the root
     * process. Ignored in the non-root processes.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1. Thrown if
     * <TT>dstarray</TT>'s length does not equal the size of this communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>dstarray</TT> or any element thereof is null. Thrown if <TT>src</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void gather(int root,
            int tag,
            Buf src,
            Buf[] dstarray)
            throws IOException {
        // Verify preconditions.
        if (0 > root || root >= mySize) {
            throw new IndexOutOfBoundsException("Comm.gather(): root = " + root + " out of bounds");
        }

        // A gather is done as a series of point-to-point messages. The root
        // process receives a separate message from every other process. Here is
        // the message pattern for a communicator with 8 processes gathering
        // into root process 0:
        //
        //     Process
        //     0     1     2     3     4     5     6     7
        //     |     |     |     |     |     |     |     |
        //     |<----|     |     |     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |<----------|     |     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |<----------------|     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |<----------------------|     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |<----------------------------|     |     |
        //     |     |     |     |     |     |     |     |
        //     |<----------------------------------|     |
        //     |     |     |     |     |     |     |     |
        //     |<----------------------------------------|
        //     |     |     |     |     |     |     |     |
        // Root process receives all messages.
        if (myRank == root) {
            // Array of IORequest objects for non-blocking receives.
            IORequest[] iorequest = new IORequest[mySize];

            // Initiate receives from lower-ranked processes.
            for (int rank = 0; rank < myRank; ++rank) {
                iorequest[rank] = new IORequest();
                myChannelGroup.receiveNoWait(getChannel(rank), tag, dstarray[rank], iorequest[rank]);
            }

            // Initiate receives from higher-ranked processes.
            for (int rank = myRank + 1; rank < mySize; ++rank) {
                iorequest[rank] = new IORequest();
                myChannelGroup.receiveNoWait(getChannel(rank), tag, dstarray[rank], iorequest[rank]);
            }

            // Copy to itself.
            dstarray[myRank].copy(src);

            // Wait for completion of receives from lower-ranked processes.
            for (int rank = 0; rank < myRank; ++rank) {
                iorequest[rank].waitForFinish();
            }

            // Wait for completion of receives from higher-ranked processes.
            for (int rank = myRank + 1; rank < mySize; ++rank) {
                iorequest[rank].waitForFinish();
            }
        } // Non-root process sends one message.
        else {
            myChannelGroup.send(getChannel(root), tag, src);
        }
    }

    /**
     * All-gather messages from each process to all processes in this
     * communicator. A message tag of 0 is used. All processes must call
     * <TT>allGather()</TT>.
     * <P>
     * Each process sends the message items in the given source buffer. When the
     * <TT>allGather()</TT> method returns, the source buffer has been fully
     * sent.
     * <P>
     * Each process receives message items from the other processes. The message
     * items received from process <I>i</I> are stored in the destination buffer
     * at index <I>i</I> in the given array of destination buffers. This
     * destination buffer must have the same length and the same item data type
     * as the source buffer in process <I>i</I>. When the <TT>allGather()</TT>
     * method returns, all the destination buffers have been fully received.
     * <P>
     * All-gather is the same as gather, except that every process has an array
     * of destination buffers, and every process receives the results of the
     * gather.
     *
     * @param src Source buffer to be sent.
     * @param dstarray Array of destination buffers to be received.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>dstarray</TT>'s length does not equal the size of this communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>dstarray</TT> or any element thereof is null. Thrown if <TT>src</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void allGather(Buf src,
            Buf[] dstarray)
            throws IOException {
        allGather(0, src, dstarray);
    }

    /**
     * All-gather messages from each process to all processes in this
     * communicator using the given message tag. All processes must call
     * <TT>allGather()</TT> with the same value for <TT>tag</TT>.
     * <P>
     * Each process sends the message items in the given source buffer. When the
     * <TT>allGather()</TT> method returns, the source buffer has been fully
     * sent.
     * <P>
     * Each process receives message items from the other processes. The message
     * items received from process <I>i</I> are stored in the destination buffer
     * at index <I>i</I> in the given array of destination buffers. This
     * destination buffer must have the same length and the same item data type
     * as the source buffer in process <I>i</I>. When the <TT>allGather()</TT>
     * method returns, all the destination buffers have been fully received.
     * <P>
     * All-gather is the same as gather, except that every process has an array
     * of destination buffers, and every process receives the results of the
     * gather.
     *
     * @param tag Message tag.
     * @param src Source buffer to be sent.
     * @param dstarray Array of destination buffers to be received.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>dstarray</TT>'s length does not equal the size of this communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>dstarray</TT> or any element thereof is null. Thrown if <TT>src</TT>
     * is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void allGather(int tag,
            Buf src,
            Buf[] dstarray)
            throws IOException {
        // Get ranks of predecessor and successor processes.
        int pred = (myRank - 1 + mySize) % mySize;
        int succ = (myRank + 1) % mySize;

        // Copy source buffer into destination buffer at my own rank.
        dstarray[myRank].copy(src);

        // Do (mySize-1) message rounds. Messages are sent in a pipelined
        // fashion from each process to its predecessor until each process's
        // source data has arrived in every process. Each outgoing message is
        // overlapped with an incoming message.
        for (int i = 1; i < mySize; ++i) {
            sendReceive(/*toRank  */pred,
                    /*sendTag */ tag,
                    /*sendBuf */ dstarray[(myRank + i - 1) % mySize],
                    /*fromRank*/ succ,
                    /*recvTag */ tag,
                    /*recvBuf */ dstarray[(myRank + i) % mySize]);
        }
    }

    /**
     * Perform a reduction on all processes in this communicator. The reduction
     * uses a message tag of 0. All processes must call <TT>reduce()</TT> with
     * the same value for <TT>root</TT>, with a buffer of the same length and
     * the same item data type, and with the same binary operation (class
     * {@linkplain edu.rit.pj.reduction.Op Op}).
     * <P>
     * Before calling <TT>reduce()</TT>, each process has a buffer filled with
     * data items. After <TT>reduce()</TT> returns, each data item in the root
     * process's buffer has been set to the <B>reduction</B> of the
     * corresponding data items in all the processes' buffers. The reduction is
     * calculated by this formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process rank 0, 1, 2, and so on. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     * <P>
     * In the root process, the reduce operation always changes the buffer's
     * contents as described above. In the non-root processes, the reduce
     * operation may or may not change the buffer's contents; the final contents
     * of the buffer in the non-root processes is not specified.
     * <P>
     * When the <TT>reduce()</TT> method returns in the root process, the
     * reduction has been fully performed as described above. When the
     * <TT>reduce()</TT> method returns in a non-root process, the non-root
     * process has sent all its data items into the reduction, but the reduction
     * may not be fully complete in the root process yet.
     *
     * @param root Root process's rank in this communicator.
     * @param buffer Buffer of data items to be reduced.
     * @param op Binary operation.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void reduce(int root,
            Buf buffer,
            Op op)
            throws IOException {
        reduce(root, 0, buffer, op);
    }

    /**
     * Perform a reduction on all processes in this communicator using the given
     * message tag. All processes must call <TT>reduce()</TT> with the same
     * value for <TT>root</TT>, with a buffer of the same length and the same
     * item data type, and with the same binary operation (class {@linkplain
     * edu.rit.pj.reduction.Op Op}).
     * <P>
     * Before calling <TT>reduce()</TT>, each process has a buffer filled with
     * data items. After <TT>reduce()</TT> returns, each data item in the root
     * process's buffer has been set to the <B>reduction</B> of the
     * corresponding data items in all the processes' buffers. The reduction is
     * calculated by this formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process rank 0, 1, 2, and so on. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     * <P>
     * In the root process, the reduce operation always changes the buffer's
     * contents as described above. In the non-root processes, the reduce
     * operation may or may not change the buffer's contents; the final contents
     * of the buffer in the non-root processes is not specified.
     * <P>
     * When the <TT>reduce()</TT> method returns in the root process, the
     * reduction has been fully performed as described above. When the
     * <TT>reduce()</TT> method returns in a non-root process, the non-root
     * process has sent all its data items into the reduction, but the reduction
     * may not be fully complete in the root process yet.
     *
     * @param root Root process's rank in this communicator.
     * @param tag Message tag.
     * @param buffer Buffer of data items to be reduced.
     * @param op Binary operation.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>root</TT> is not in the range 0 .. <TT>size()</TT>-1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void reduce(int root,
            int tag,
            Buf buffer,
            Op op)
            throws IOException {
        // Verify preconditions.
        if (0 > root || root >= mySize) {
            throw new IndexOutOfBoundsException("Comm.reduce(): root = " + root + " out of bounds");
        }

        // Early return if only one process.
        if (mySize == 1) {
            return;
        }

        // A reduction is done as a series of point-to-point messages. The
        // messages are organized into rounds. The number of rounds is
        // ceil(log_2(mySize)). The message pattern is the reverse of the
        // broadcast message pattern. In each round, processes receive messages
        // from other processes and reduce the data items into their accumulator
        // buffers in parallel. When a process has received all messages, it
        // sends the reduced results on. Here is the message pattern for a
        // communicator with 8 processes doing a reduction into root process 0:
        //
        //     Process
        //     0     1     2     3     4     5     6     7
        //     |     |     |     |     |     |     |     |
        //     |<----------------------|     |     |     |  Round 1
        //     |     |<----------------------|     |     |
        //     |     |     |<----------------------|     |
        //     |     |     |     |<----------------------|
        //     |     |     |     |     |     |     |     |
        //     |<----------|     |     |     |     |     |  Round 2
        //     |     |<----------|     |     |     |     |
        //     |     |     |     |     |     |     |     |
        //     |<----|     |     |     |     |     |     |  Round 3
        //     |     |     |     |     |     |     |     |
        //
        // If a process other than process 0 is the root, the message pattern is
        // the same, except the process ranks are circularly rotated.
        // Get array of process ranks in the broadcast tree.
        int[] broadcasttree = getBroadcastTree(root);
        int n = broadcasttree.length;

        // Set up reduction buffer on top of source buffer.
        Buf reductionbuf = buffer.getReductionBuf(op);

        // Receive data from children if any, one at a time in reverse order.
        for (int i = n - 1; i >= 1; --i) {
            int child = broadcasttree[i];
            myChannelGroup.receive(getChannel(child), tag, reductionbuf);
        }

        // Send data to parent if any.
        int parent = broadcasttree[0];
        if (parent != -1) {
            myChannelGroup.send(getChannel(parent), tag, buffer);
        }
    }

    /**
     * Perform an all-reduce on all processes in this communicator. The
     * all-reduce uses a message tag of 0. All processes must call
     * <TT>allReduce()</TT> with a buffer of the same length and the same item
     * data type, and with the same binary operation (class {@linkplain
     * edu.rit.pj.reduction.Op Op}).
     * <P>
     * Before calling <TT>allReduce()</TT>, each process has a buffer filled
     * with data items. After <TT>allReduce()</TT> returns, each data item in
     * the calling process's buffer has been set to the <B>reduction</B> of the
     * corresponding data items in all the processes' buffers. The reduction is
     * calculated by this formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process rank 0, 1, 2, and so on. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     * <P>
     * The <TT>allReduce()</TT> method is similar to the <TT>reduce()</TT>
     * method, except the results are stored in all the processes' buffers, not
     * just the one root process's buffer.
     *
     * @param buffer Buffer of data items to be reduced.
     * @param op Binary operation.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void allReduce(Buf buffer,
            Op op)
            throws IOException {
        allReduce(0, buffer, op);
    }

    /**
     * Perform an all-reduce on all processes in this communicator using the
     * given message tag. All processes must call <TT>allReduce()</TT> with a
     * buffer of the same length and the same item data type, and with the same
     * binary operation (class {@linkplain edu.rit.pj.reduction.Op Op}).
     * <P>
     * Before calling <TT>allReduce()</TT>, each process has a buffer filled
     * with data items. After <TT>allReduce()</TT> returns, each data item in
     * the calling process's buffer has been set to the <B>reduction</B> of the
     * corresponding data items in all the processes' buffers. The reduction is
     * calculated by this formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process rank 0, 1, 2, and so on. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     * <P>
     * The <TT>allReduce()</TT> method is similar to the <TT>reduce()</TT>
     * method, except the results are stored in all the processes' buffers, not
     * just the one root process's buffer.
     *
     * @param tag Message tag.
     * @param buffer Buffer of data items to be reduced.
     * @param op Binary operation.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void allReduce(int tag,
            Buf buffer,
            Op op)
            throws IOException {
		// An all-reduce is done using a "butterfly" message passing pattern.
        // Consider the case of K=8 processes. In the first round, processes one
        // rank apart exchange data, then each processes accumulates the data
        // from the other process using the reduction operator. In the second
        // round, processes two ranks apart exchange and accumulate data. In the
        // third round, processes four ranks apart exchange and accumulate data.
        //
        //     Process
        //     0     1     2     3     4     5     6     7
        //     |     |     |     |     |     |     |     |
        //     |<--->|     |<--->|     |<--->|     |<--->|  Round 1
        //     |     |     |     |     |     |     |     |
        //     |<--------->|     |     |<--------->|     |  Round 2
        //     |     |<--------->|     |     |<--------->|
        //     |     |     |     |     |     |     |     |
        //     |<--------------------->|     |     |     |  Round 3
        //     |     |<--------------------->|     |     |
        //     |     |     |<--------------------->|     |
        //     |     |     |     |<--------------------->|
        //     |     |     |     |     |     |     |     |
        //
        // The butterfly pattern works only if the number of processes is a
        // power of two. If this is not the case, there are two extra message
        // rounds. Each process outside the butterfly pattern sends its data to
        // its counterpart inside the butterfly pattern, which accumulates the
        // data using the reduction operator. Then the butterfly pattern takes
        // place. Afterwards, each process outside the butterfly pattern
        // receives the final result from its counterpart inside the butterfly
        // pattern. For the case of K=10 processes:
        //
        //     Process
        //     0     1     2     3     4     5     6     7     8     9
        //     |     |     |     |     |     |     |     |     |     |
        //     |<----------------------------------------------|     | Pre
        //     |     |<----------------------------------------------|
        //     |     |     |     |     |     |     |     |     |     |
        //     |<--->|     |<--->|     |<--->|     |<--->|     |     | Round 1
        //     |     |     |     |     |     |     |     |     |     |
        //     |<--------->|     |     |<--------->|     |     |     | Round 2
        //     |     |<--------->|     |     |<--------->|     |     |
        //     |     |     |     |     |     |     |     |     |     |
        //     |<--------------------->|     |     |     |     |     | Round 3
        //     |     |<--------------------->|     |     |     |     |
        //     |     |     |<--------------------->|     |     |     |
        //     |     |     |     |<--------------------->|     |     |
        //     |     |     |     |     |     |     |     |     |     |
        //     |---------------------------------------------->|     | Post
        //     |     |---------------------------------------------->|
        //     |     |     |     |     |     |     |     |     |     |
        //
        // If K is a power of two, the all-reduce takes (log_2 K) rounds. If K
        // is not a power of two, the all-reduce takes floor(log_2 K)+2 rounds.

        // Early exit if only one process.
        if (mySize == 1) {
            return;
        }

        // Determine the highest power of 2 less than or equal to this
        // communicator's size. Processes at this rank and above will be outside
        // the butterfly message passing pattern.
        int outside = mySizePowerOf2;

        // For processes outside the butterfly:
        if (myRank >= outside) {
            int insideRank = myRank - outside;

            // Send initial data to counterpart inside.
            send(insideRank, tag, buffer);

            // Receive reduced result from counterpart inside.
            receive(insideRank, tag, buffer);
        } // For processes inside the butterfly:
        else {
            // Set up temporary receive buffer.
            Buf receiveBuf = buffer.getTemporaryBuf();

            // Set up reduction buffer on top of data buffer.
            Buf reductionBuf = buffer.getReductionBuf(op);

            // If there is a counterpart process outside, receive and accumulate
            // its initial data.
            int outsideRank = myRank + outside;
            if (outsideRank < mySize) {
                receive(outsideRank, tag, reductionBuf);
            }

            // Perform butterfly message passing rounds with counterpart
            // processes inside.
            int round = 1;
            while (round < outside) {
                int otherRank = myRank ^ round;
                sendReceive(otherRank, tag, buffer, otherRank, tag, receiveBuf);
                reductionBuf.copy(receiveBuf);
                round <<= 1;
            }

            // If there is a counterpart process outside, send the reduced
            // result.
            if (outsideRank < mySize) {
                send(outsideRank, tag, buffer);
            }
        }
    }

    /**
     * Do an all-to-all among all processes in this communicator. A message tag
     * of 0 is used.
     * <P>
     * <TT>srcarray</TT> must be an array of <I>K</I> buffers, where <I>K</I> is
     * the size of this communicator. <TT>dstarray</TT> must be an array of
     * <I>K</I> buffers referring to different storage from the source buffers.
     * For each process rank <I>k</I>, 0 &lt;= <I>k</I> &lt;= <I>K</I>, and each
     * buffer index <I>i</I>, 0 &lt;= <I>i</I> &lt;= <I>K</I>, the contents of
     * <TT>srcarray[k]</TT> in process <I>i</I> are sent to <TT>dstarray[i]</TT>
     * in process <I>k</I>.
     *
     * @param srcarray Array of source buffers to be sent by this process.
     * @param dstarray Array of destination buffers to be received by this
     * process.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>srcarray</TT>'s length does not equal the size of this communicator.
     * Thrown if <TT>dstarray</TT>'s length does not equal the size of this
     * communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>srcarray</TT> or any element thereof is null. Thrown if
     * <TT>dstarray</TT> or any element thereof is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void allToAll(Buf[] srcarray,
            Buf[] dstarray)
            throws IOException {
        allToAll(0, srcarray, dstarray);
    }

    /**
     * Do an all-to-all among all processes in this communicator using the given
     * message tag. All processes must call <TT>allToAll()</TT> with the same
     * value for <TT>tag</TT>.
     * <P>
     * <TT>srcarray</TT> must be an array of <I>K</I> buffers, where <I>K</I> is
     * the size of this communicator. <TT>dstarray</TT> must be an array of
     * <I>K</I> buffers referring to different storage from the source buffers.
     * For each process rank <I>k</I>, 0 &lt;= <I>k</I> &lt;= <I>K</I>, and each
     * buffer index <I>i</I>, 0 &lt;= <I>i</I> &lt;= <I>K</I>, the contents of
     * <TT>srcarray[k]</TT> in process <I>i</I> are sent to <TT>dstarray[i]</TT>
     * in process <I>k</I>.
     *
     * @param tag Message tag.
     * @param srcarray Array of source buffers to be sent by this process.
     * @param dstarray Array of destination buffers to be received by this
     * process.
     *
     * @exception IndexOutOfBoundsException (unchecked exception) Thrown if
     * <TT>srcarray</TT>'s length does not equal the size of this communicator.
     * Thrown if <TT>dstarray</TT>'s length does not equal the size of this
     * communicator.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>srcarray</TT> or any element thereof is null. Thrown if
     * <TT>dstarray</TT> or any element thereof is null.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void allToAll(int tag,
            Buf[] srcarray,
            Buf[] dstarray)
            throws IOException {
		// An all-to-all is done as a series of send-receives. Each process
        // sends the appropriate buffer to the process one ahead and receives
        // the appropriate buffer from the process one behind. Then each process
        // sends the appropriate buffer to the process two ahead and receives
        // the appropriate buffer from the process two behind. And so on. Here
        // is the message pattern for a communicator with 4 processes doing an
        // all-to-all:
        //
        //          Process
        //          0     1     2     3
        //          |     |     |     |
        //          |---->|     |     | Round 1
        //          |     |---->|     |
        //          |     |     |---->|
        //     - -->|     |     |     |--- -
        //          |     |     |     |
        //          |     |     |     |
        //          |---------->|     | Round 2
        //          |     |---------->|
        //     - -->|     |     |--------- -
        //     - -------->|     |     |--- -
        //          |     |     |     |
        //          |     |     |     |
        //          |---------------->| Round 3
        //     - -->|     |--------------- -
        //     - -------->|     |--------- -
        //     - -------------->|     |--- -
        //          |     |     |     |

        // Copy source to destination at this process's own rank.
        dstarray[myRank].copy(srcarray[myRank]);

        // Initiate K-1 non-blocking send-receives.
        CommRequest[] commrequest = new CommRequest[mySize];
        for (int i = 1; i < mySize; ++i) {
            int toRank = (myRank + i) % mySize;
            int fromRank = (myRank - i + mySize) % mySize;
            commrequest[i]
                    = sendReceive(toRank, tag, srcarray[toRank],
                            fromRank, tag, dstarray[fromRank],
                            (CommRequest) null);
        }

        // Wait for completion of all send-receives.
        for (int i = 1; i < mySize; ++i) {
            commrequest[i].waitForFinish();
        }
    }

    /**
     * Perform a scan on all processes in this communicator. A message tag of 0
     * is used. All processes must call <TT>scan()</TT> with a buffer of the
     * same length and the same item data type, and with the same binary
     * operation (class {@linkplain edu.rit.pj.reduction.Op Op}).
     * <P>
     * Before calling <TT>scan()</TT>, each process has a buffer filled with
     * data items. After <TT>scan()</TT> returns, each data item in the buffer
     * of process rank <I>i</I> has been set to the <B>reduction</B> of the
     * corresponding data items in the buffers of process ranks 0 through
     * <I>i</I>. The reduction is calculated by this formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process ranks 0 through <I>i</I>. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     *
     * @param buf Buffer of data items to be scanned.
     * @param op Binary operation.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void scan(Buf buf,
            Op op)
            throws IOException {
        scan(0, buf, op);
    }

    /**
     * Perform a scan on all processes in this communicator using the given
     * message tag. All processes must call <TT>scan()</TT> with the same value
     * for <TT>tag</TT>, with a buffer of the same length and the same item data
     * type, and with the same binary operation (class {@linkplain
     * edu.rit.pj.reduction.Op Op}).
     * <P>
     * Before calling <TT>scan()</TT>, each process has a buffer filled with
     * data items. After <TT>scan()</TT> returns, each data item in the buffer
     * of process rank <I>i</I> has been set to the <B>reduction</B> of the
     * corresponding data items in the buffers of process ranks 0 through
     * <I>i</I>. The reduction is calculated by this formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process ranks 0 through <I>i</I>. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     *
     * @param tag Message tag.
     * @param buf Buffer of data items to be scanned.
     * @param op Binary operation.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void scan(int tag,
            Buf buf,
            Op op)
            throws IOException {
        // Early return if only one process.
        if (mySize == 1) {
            return;
        }

        // A scan is done as a series of point-to-point messages. The messages
        // are organized into rounds. The number of rounds is
        // ceil(log_2(mySize)). In the first round, each process sends its data
        // to the process one rank ahead, and the incoming data is combined with
        // the process's data using the reduction operator. In the second round,
        // each process sends its data to the process two ranks ahead. In the
        // third round, each process sends its data to process four ranks ahead.
        // And so on. Here is the message pattern for a communicator with 8
        // processes:
        //
        //     Process
        //     0     1     2     3     4     5     6     7
        //     |     |     |     |     |     |     |     |
        //     |---->|---->|---->|---->|---->|---->|---->|  Round 1
        //     |     |     |     |     |     |     |     |
        //     |---------->|     |     |     |     |     |  Round 2
        //     |     |---------->|     |     |     |     |
        //     |     |     |---------->|     |     |     |
        //     |     |     |     |---------->|     |     |
        //     |     |     |     |     |---------->|     |
        //     |     |     |     |     |     |---------->|
        //     |     |     |     |     |     |     |     |
        //     |---------------------->|     |     |     |  Round 3
        //     |     |---------------------->|     |     |
        //     |     |     |---------------------->|     |
        //     |     |     |     |---------------------->|
        //
        // Get temporary buffer for holding incoming data items.
        Buf tempbuf = buf.getTemporaryBuf();

        // Get reduction buffer for combining data items.
        Buf reductionbuf = buf.getReductionBuf(op);

        // Do rounds of message passing and reduction.
        int skip = 1;
        for (;;) {
            int toRank = myRank + skip;
            int fromRank = myRank - skip;
            boolean toExists = 0 <= toRank && toRank < mySize;
            boolean fromExists = 0 <= fromRank && fromRank < mySize;
            if (toExists && fromExists) {
                sendReceive(toRank, tag, buf, fromRank, tag, tempbuf);
                reductionbuf.copy(tempbuf);
            } else if (fromExists) {
                receive(fromRank, tag, reductionbuf);
            } else if (toExists) {
                send(toRank, tag, buf);
            } else {
                break;
            }
            skip <<= 1;
        }
    }

    /**
     * Perform an exclusive scan on all processes in this communicator. A
     * message tag of 0 is used. All processes must call
     * <TT>exclusiveScan()</TT> with a buffer of the same length and the same
     * item data type, with the same binary operation (class {@linkplain
     * edu.rit.pj.reduction.Op Op}), and with the same initial data value.
     * <P>
     * Before calling <TT>exclusiveScan()</TT>, each process has a buffer filled
     * with data items. After <TT>exclusiveScan()</TT> returns, each data item
     * in the buffer of process rank <I>i</I> &gt; 0 has been set to the
     * <B>reduction</B> of the corresponding data items in the buffers of
     * process ranks 0 through <I>i</I>-1. The reduction is calculated by this
     * formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process ranks 0 through <I>i</I>-1. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     * <P>
     * In process 0, each data item in the buffer has been set to the initial
     * data value using the buffer's <TT>fill()</TT> method.
     * <P>
     * If the buffer's item data type is a primitive type, the <TT>item</TT>
     * must be an instance of the corresponding primitive wrapper class -- class
     * Integer for type <TT>int</TT>, class Double for type <TT>double</TT>, and
     * so on. If the <TT>item</TT> is null, the item data type's default initial
     * value is assigned to each element in the buffer.
     * <P>
     * If the buffer's item data type is a nonprimitive type, the <TT>item</TT>
     * must be an instance of the item class or a subclass thereof. The
     * <TT>item</TT> may be null. Note that since <TT>item</TT> is
     * <I>assigned</I> to every buffer element, every buffer element ends up
     * referring to the same <TT>item</TT>.
     *
     * @param buf Buffer of data items to be scanned.
     * @param op Binary operation.
     * @param item Initial data value.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void exclusiveScan(Buf buf,
            Op op,
            Object item)
            throws IOException {
        exclusiveScan(0, buf, op, item);
    }

    /**
     * Perform an exclusive scan on all processes in this communicator using the
     * given message tag. All processes must call <TT>exclusiveScan()</TT> with
     * the same value for <TT>tag</TT>, with a buffer of the same length and the
     * same item data type, with the same binary operation (class {@linkplain
     * edu.rit.pj.reduction.Op Op}), and with the same initial data value.
     * <P>
     * Before calling <TT>exclusiveScan()</TT>, each process has a buffer filled
     * with data items. After <TT>exclusiveScan()</TT> returns, each data item
     * in the buffer of process rank <I>i</I> &gt; 0 has been set to the
     * <B>reduction</B> of the corresponding data items in the buffers of
     * process ranks 0 through <I>i</I>-1. The reduction is calculated by this
     * formula:
     * <P>
     * &nbsp;&nbsp;&nbsp;&nbsp;<I>item</I><SUB>0</SUB> <I>op</I>
     * <I>item</I><SUB>1</SUB> <I>op</I> <I>item</I><SUB>2</SUB> <I>op</I> . . .
     * <P>
     * where <I>op</I> is the binary operation passed in as an argument and
     * <I>item</I><SUB>0</SUB>, <I>item</I><SUB>1</SUB>,
     * <I>item</I><SUB>2</SUB>, and so on are the data items in the buffers of
     * process ranks 0 through <I>i</I>-1. However, the order in which the data
     * items actually are combined is not specified. Therefore, the binary
     * operation must be such that the answer will be the same regardless of the
     * order in which the data items are combined; that is, the binary operation
     * must be commutative and associative.
     * <P>
     * In process 0, each data item in the buffer has been set to the initial
     * data value using the buffer's <TT>fill()</TT> method.
     * <P>
     * If the buffer's item data type is a primitive type, the <TT>item</TT>
     * must be an instance of the corresponding primitive wrapper class -- class
     * Integer for type <TT>int</TT>, class Double for type <TT>double</TT>, and
     * so on. If the <TT>item</TT> is null, the item data type's default initial
     * value is assigned to each element in the buffer.
     * <P>
     * If the buffer's item data type is a nonprimitive type, the <TT>item</TT>
     * must be an instance of the item class or a subclass thereof. The
     * <TT>item</TT> may be null. Note that since <TT>item</TT> is
     * <I>assigned</I> to every buffer element, every buffer element ends up
     * referring to the same <TT>item</TT>.
     *
     * @param tag Message tag.
     * @param buf Buffer of data items to be scanned.
     * @param op Binary operation.
     * @param item Initial data value.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>buf</TT> is null or <TT>op</TT>
     * is null.
     * @exception ClassCastException (unchecked exception) Thrown if
     * <TT>buf</TT> and <TT>op</TT> do not use the same item data type.
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void exclusiveScan(int tag,
            Buf buf,
            Op op,
            Object item)
            throws IOException {
        // An exclusive scan begins with each process sending its buffer to the
        // next higher process. Then process 0 fills its buffer with the initial
        // data value, while processes 1 and higher do an inclusive scan.

        int toRank;
        int fromRank;
        boolean toExists;
        boolean fromExists;

        // Process 0 does this.
        if (myRank == 0) {
            // Send buffer to next higher process.
            toRank = 1;
            toExists = toRank < mySize;
            if (toExists) {
                send(toRank, tag, buf);
            }

            // Fill buffer with initial data value.
            buf.fill(item);
        } // Processes 1 and higher do this.
        else {
            // Get temporary buffer for holding incoming data items.
            Buf tempbuf = buf.getTemporaryBuf();

            // Get reduction buffer for combining data items.
            Buf reductionbuf = buf.getReductionBuf(op);

            // Send buffer to next higher process.
            toRank = myRank + 1;
            fromRank = myRank - 1;
            toExists = 0 <= toRank && toRank < mySize;
            if (toExists) {
                sendReceive(toRank, tag, buf, fromRank, tag, tempbuf);
                buf.copy(tempbuf);
            } else {
                receive(fromRank, tag, buf);
            }

            // Do rounds of message passing and reduction.
            int skip = 1;
            for (;;) {
                toRank = myRank + skip;
                fromRank = myRank - skip;
                toExists = 1 <= toRank && toRank < mySize;
                fromExists = 1 <= fromRank && fromRank < mySize;
                if (toExists && fromExists) {
                    sendReceive(toRank, tag, buf, fromRank, tag, tempbuf);
                    reductionbuf.copy(tempbuf);
                } else if (fromExists) {
                    receive(fromRank, tag, reductionbuf);
                } else if (toExists) {
                    send(toRank, tag, buf);
                } else {
                    break;
                }
                skip <<= 1;
            }
        }
    }

    /**
     * Cause all processes in this communicator to wait at a barrier. The
     * barrier uses a message tag of 0. All processes must call
     * <TT>barrier()</TT>. The calling thread blocks until every process has
     * called <TT>barrier()</TT>, then the calling thread unblocks and returns
     * from the <TT>barrier()</TT> call.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void barrier()
            throws IOException {
        barrier(0);
    }

    /**
     * Cause all processes in this communicator to wait at a barrier, using the
     * given message tag. All processes must call <TT>barrier()</TT> with the
     * same tag. The calling thread blocks until every process has called
     * <TT>barrier()</TT>, then the calling thread unblocks and returns from the
     * <TT>barrier()</TT> call.
     *
     * @param tag Message tag.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    public void barrier(int tag)
            throws IOException {
        // A barrier is done as an all-reduce of an empty buffer.
        allReduce(tag, IntegerBuf.emptyBuffer(), IntegerOp.SUM);
    }

    /**
     * Returns a string version of this communicator. The string includes the
     * communicator's size, the current process's rank, and the host and port of
     * each backend process.
     *
     * @return String version.
     */
    public String toString() {
        StringBuilder buf = new StringBuilder();
        buf.append("Comm(size=");
        buf.append(mySize);
        buf.append(",rank=");
        buf.append(myRank);
        buf.append(",backend");
        for (int i = 0; i < mySize; ++i) {
            if (i > 0) {
                buf.append(',');
            }
            buf.append('[');
            buf.append(i);
            buf.append("]=");
            buf.append(myAddressForRank[i]);
        }
        buf.append(')');
        return buf.toString();
    }

    /**
     * Dump the state of this communicator on the given print stream. For
     * debugging.
     *
     * @param out Print stream.
     * @param prefix String to print at the beginning of each line.
     */
    public void dump(PrintStream out,
            String prefix) {
        out.println();
        out.println(prefix + getClass().getName() + "@" + Integer.toHexString(System.identityHashCode(this)));
        out.println(prefix + "mySize = " + mySize);
        out.println(prefix + "myRank = " + myRank);
        out.println(prefix + "myHost = " + myHost);
        out.println(prefix + "mySizePowerOf2 = " + mySizePowerOf2);
        out.println(prefix + "myChannelGroup = " + myChannelGroup);
        out.println(prefix + "myAddressForRank:");
        for (int i = 0; i < myAddressForRank.length; ++i) {
            out.println(prefix + "\t[" + i + "] " + myAddressForRank[i]);
        }
        out.println(prefix + "myChannelForRank:");
        for (int i = 0; i < myChannelForRank.length; ++i) {
            out.println(prefix + "\t[" + i + "] " + myChannelForRank[i]);
        }
        out.println(prefix + "myBroadcastTree:");
        for (int i = 0; i < myBroadcastTree.length; ++i) {
            out.print(prefix + "\t[" + i + "]");
            int[] tree = myBroadcastTree[i];
            if (tree == null) {
                out.print(" null");
            } else {
                for (int j = 0; j < tree.length; ++j) {
                    out.print(" " + tree[j]);
                }
            }
            out.println();
        }
        out.println();
        myChannelGroup.dump(out, prefix);
    }

// Hidden operations.
    /**
     * Notify that another process connected a channel to this process.
     *
     * @param theChannel Channel.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private synchronized void doFarEndConnected(Channel theChannel)
            throws IOException {
        // Record channel and rank.
        myChannelForRank[getFarRank(theChannel)] = theChannel;

        // Notify any threads waiting in getChannel().
        notifyAll();
    }

    /**
     * Ensure that a channel for communicating with the process at the given
     * rank is or will be set up.
     *
     * @param farrank Rank of far end process.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private synchronized void ensureChannel(int farrank)
            throws IOException {
        // Get channel from channel array.
        Channel channel = myChannelForRank[farrank];

        // If the channel does not exist:
        if (channel == null) {
            // If this is the lower-ranked process, set up the connection.
            if (myRank < farrank) {
                myChannelForRank[farrank]
                        = myChannelGroup.connect(myAddressForRank[farrank]);
            }

            // If this is the higher-ranked process, the lower-ranked process
            // will set up the connection.
        }
    }

    /**
     * Get the channel for communicating with the process at the given rank.
     *
     * @param farrank Rank of far end process.
     *
     * @return Channel.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    private synchronized Channel getChannel(int farrank)
            throws IOException {
        // Get channel from channel array.
        Channel channel = myChannelForRank[farrank];

        // If the channel does not exist:
        if (channel == null) {
            // If this is the lower-ranked process, set up the connection.
            if (myRank < farrank) {
                channel = myChannelGroup.connect(myAddressForRank[farrank]);
                myChannelForRank[farrank] = channel;
            } // If this is the higher-ranked process, wait for the lower-ranked
            // process to set up the connection.
            else {
                try {
                    while (channel == null) {
                        wait();
                        channel = myChannelForRank[farrank];
                    }
                } catch (InterruptedException exc) {
                    IOException exc2 = new InterruptedIOException();
                    exc2.initCause(exc);
                    throw exc2;
                }
            }
        }

        return channel;
    }

    /**
     * Get the rank of the process at the far end of the given channel.
     *
     * @param channel Channel.
     *
     * @return Far end process rank.
     */
    static int getFarRank(Channel channel) {
        return channel.farEndChannelGroupId();
    }

    /**
     * Get an array of process ranks in the broadcast tree for the given root.
     * The broadcast tree is cached in the field myBroadcastTree for later use.
     *
     * @param root Root process's rank.
     *
     * @return Broadcast tree.
     */
    private synchronized int[] getBroadcastTree(int root) {
        if (myBroadcastTree == null) {
            myBroadcastTree = new int[mySize][];
        }
        int[] broadcasttree = myBroadcastTree[root];
        if (broadcasttree == null) {
            broadcasttree = CommPattern.broadcastPattern(mySize, myRank, root);
            myBroadcastTree[root] = broadcasttree;
        }
        return broadcasttree;
    }

}
