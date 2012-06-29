//******************************************************************************
//
// File:    WorkerIteration.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.WorkerIteration
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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

import edu.rit.mp.ObjectBuf;

import edu.rit.mp.buf.ObjectItemBuf;

import edu.rit.util.Range;

import java.io.IOException;

/**
 * Class WorkerIteration is the abstract base class for a worker iteration that
 * is executed inside a {@linkplain WorkerRegion}. The worker iteration lets you
 * iterate over a group of items, with a separate worker team thread processing
 * each item. The generic type parameter T specifies the items' data type. The
 * items can be the elements of an array, the items obtained from an {@linkplain
 * java.util.Iterator Iterator}, or the items contained in an {@linkplain
 * java.lang.Iterable Iterable} collection.
 * <P>
 * To execute a worker iteration, create a {@linkplain WorkerRegion} object;
 * create an instance of a concrete subclass of class WorkerIteration; and pass
 * this instance to the worker region's <TT>execute()</TT> method. Either every
 * worker team thread must call the worker region's <TT>execute()</TT> method
 * with identical arguments, or every thread must not call the
 * <TT>execute()</TT> method. You can do all this using an anonymous inner
 * class; for example:
 * <PRE>
 *     new WorkerRegion()
 *         {
 *         ArrayList&lt;String&gt; list = new ArrayList&lt;String&gt;();
 *         . . .
 *         public void run()
 *             {
 *             . . .
 *             execute (list, new WorkerIteration&lt;String&gt;()
 *                 {
 *                 // Thread local variable declarations
 *                 . . .
 *                 public void start()
 *                     {
 *                     // Per-thread pre-loop initialization code
 *                     . . .
 *                     }
 *                 public void run (String item)
 *                     {
 *                     // Loop code
 *                     . . .
 *                     }
 *                 public void finish()
 *                     {
 *                     // Per-thread post-loop finalization code
 *                     . . .
 *                     }
 *                 });
 *             }
 *         . . .
 *         }
 * </PRE>
 * <P>
 * In each process of a cluster parallel program, the worker team has one or
 * more worker threads. Every worker thread in every process has a unique worker
 * index, going from index 0 for the first worker thread in the first process to
 * index <I>K</I>&minus;1 for the last worker thread in the last process, where
 * <I>K</I> is the total number of worker threads in all the processes. In
 * addition, in one process there is a master thread. The worker and master
 * threads all call the worker region's <TT>execute()</TT> method to execute the
 * worker for loop. However, the worker and master threads differ in their
 * actions.
 * <P>
 * The master thread does the following. The master sets up the source of items
 * to be iterated over -- either an array's elements, an iterator's items, or an
 * iterable collection's contents. The master repeatedly sends "tasks" to the
 * workers and receives "responses" from the workers. To send a task to a
 * particular worker, the master (1) sends a message containing the next item to
 * the worker's process; and (2) calls the worker iteration's
 * <TT>sendTaskInput()</TT> method. This method's default implementation does
 * nothing, but it can be overridden to send additional task input data to the
 * worker. To receive a response from a particular worker, the master (1)
 * receives a message from the worker's process containing the item that was
 * processed (whose state may have changed); and (2) calls the worker
 * iteration's <TT>receiveTaskOutput()</TT> method. This method's default
 * implementation does nothing, but it can be overridden to receive additional
 * task output data from the worker. Once all tasks have been sent to the
 * workers and all responses have been received from the workers, the master
 * returns from the worker region's <TT>execute()</TT> method.
 * <P>
 * Each worker thread does the following. The worker calls the worker
 * iteration's <TT>start()</TT> method once before processing any items. The
 * worker repeatedly receives tasks from the master and sends responses to the
 * master. To receive a task from the master, the worker (1) receives a message
 * containing the next item from the master's process; and (2) calls the worker
 * iteration's <TT>receiveTaskInput()</TT> method. This method's default
 * implementation does nothing, but it can be overridden to receive additional
 * task input data from the master. The worker now calls the worker iteration's
 * <TT>run()</TT> method, passing in the item to be processed. When the
 * <TT>run()</TT> method returns, the worker sends the response to the master.
 * To send the response, the worker (1) sends a message to the master's process
 * containing the item that was processed (whose state may have changed); and
 * (2) calls the worker iteration's <TT>sendTaskOutput()</TT> method. This
 * method's default implementation does nothing, but it can be overridden to
 * send additional task output data to the master. Once all tasks have been
 * received from the master and all responses have been sent to the master, the
 * worker calls the worker iteration's <TT>finish()</TT> method. (Unlike a
 * {@linkplain ParallelIteration}'s threads, the workers do <I>not</I>
 * synchronize with each other at a barrier at this point.) The worker then
 * returns from the worker region's <TT>execute()</TT> method.
 * <P>
 * Each message described above is sent with a message tag equal to
 * <I>W</I>+<I>T</I>, where <I>W</I> is the worker index and <I>T</I> is the
 * "tag offset." The tag offset is <TT>Integer.MIN_VALUE</TT> by default, but
 * this can be changed by overriding the <TT>tagOffset()</TT> method. Thus, the
 * message tags fall in the range <I>T</I> .. <I>K</I>&minus;1+<I>T</I>, where
 * <I>K</I> is the total number of workers in all the processes. The program
 * should not use message tags in this range except to send and receive the
 * messages described above.
 * <P>
 * Note that each worker team thread actually creates its own instance of the
 * worker iteration class and passes that instance to the worker region's
 * <TT>execute()</TT> method. Thus, any fields declared in the worker iteration
 * class will <I>not</I> be shared by all the workers, but instead will be
 * private to each worker.
 * <P>
 * The <TT>start()</TT> method is intended for performing per-thread
 * initialization before starting the loop iterations. If no such initialization
 * is needed, omit the <TT>start()</TT> method.
 * <P>
 * The <TT>run()</TT> method contains the code for the loop body. It does
 * whatever processing is needed on the one item passed in as an argument. Note
 * that, unlike a worker for loop (class {@linkplain WorkerForLoop}), a worker
 * iteration is not "chunked;" each worker team thread always processes just one
 * item at a time.
 * <P>
 * The <TT>finish()</TT> method is intended for performing per-thread
 * finalization after finishing the loop iterations. If no such finalization is
 * needed, omit the <TT>finish()</TT> method.
 * <P>
 * If the worker iteration's <TT>start()</TT>, <TT>run()</TT>, or
 * <TT>finish()</TT> method throws an exception in one of the worker threads,
 * then that worker thread executes no further code in the loop, and the worker
 * region's <TT>execute()</TT> method throws that same exception in that thread.
 * However, the other worker threads in the worker team continue to execute.
 *
 * @param  <T>  Data type of the items iterated over.
 *
 * @author  Alan Kaminsky
 * @version 07-Oct-2010
 */
public abstract class WorkerIteration<T>
	extends WorkerConstruct
	{

// Exported constructors.

	/**
	 * Construct a new worker iteration.
	 */
	public WorkerIteration()
		{
		super();
		}

// Exported operations.

	/**
	 * Perform per-thread initialization actions before starting the loop
	 * iterations. Called by a worker thread.
	 * <P>
	 * The <TT>start()</TT> method may be overridden in a subclass. If not
	 * overridden, the <TT>start()</TT> method does nothing.
	 *
	 * @exception  Exception
	 *     The <TT>start()</TT> method may throw any exception.
	 */
	public void start()
		throws Exception
		{
		}

	/**
	 * Send additional input data associated with a task. Called by the master
	 * thread. The task is denoted by the given item to be processed. The input
	 * data must be sent using the given communicator, to the given worker
	 * process rank, with the given message tag.
	 * <P>
	 * The <TT>sendTaskInput()</TT> method may be overridden in a subclass. If
	 * not overridden, the <TT>sendTaskInput()</TT> method does nothing.
	 *
	 * @param  item   Item to be processed.
	 * @param  comm   Communicator.
	 * @param  wRank  Worker process rank.
	 * @param  tag    Message tag.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void sendTaskInput
		(T item,
		 Comm comm,
		 int wRank,
		 int tag)
		throws IOException
		{
		}

	/**
	 * Receive input data associated with a task. Called by a worker thread. The
	 * task is denoted by the given item to be processed. The input data must be
	 * received using the given communicator, from the given master process
	 * rank, with the given message tag.
	 * <P>
	 * The <TT>receiveTaskInput()</TT> method may be overridden in a subclass.
	 * If not overridden, the <TT>receiveTaskInput()</TT> method does nothing.
	 *
	 * @param  item   Item to be processed.
	 * @param  comm   Communicator.
	 * @param  mRank  Master process rank.
	 * @param  tag    Message tag.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void receiveTaskInput
		(T item,
		 Comm comm,
		 int mRank,
		 int tag)
		throws IOException
		{
		}

	/**
	 * Process one item in this worker iteration. The <TT>run()</TT> method must
	 * perform the loop body for the given item.
	 * <P>
	 * The <TT>run()</TT> method must be overridden in a subclass.
	 *
	 * @param  item  Item.
	 *
	 * @exception  Exception
	 *     The <TT>run()</TT> method may throw any exception.
	 */
	public abstract void run
		(T item)
		throws Exception;

	/**
	 * Send additional output data associated with a task. Called by a worker
	 * thread. The task is denoted by the given item that was processed (whose
	 * state may have changed during processing). The output data must be sent
	 * using the given communicator, to the given master process rank, with the
	 * given message tag.
	 * <P>
	 * The <TT>sendTaskOutput()</TT> method may be overridden in a subclass. If
	 * not overridden, the <TT>sendTaskOutput()</TT> method does nothing.
	 *
	 * @param  item   Item that was processed.
	 * @param  comm   Communicator.
	 * @param  mRank  Master process rank.
	 * @param  tag    Message tag.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void sendTaskOutput
		(T item,
		 Comm comm,
		 int mRank,
		 int tag)
		throws IOException
		{
		}

	/**
	 * Receive additional output data associated with a task. Called by the
	 * master thread. The task is denoted by the given item that was processed
	 * (whose state may have changed during processing). The output data must be
	 * received using the given communicator, from the given worker process
	 * rank, with the given message tag.
	 * <P>
	 * The <TT>receiveTaskOutput()</TT> method may be overridden in a subclass.
	 * If not overridden, the <TT>receiveTaskOutput()</TT> method does nothing.
	 *
	 * @param  item   Item that was processed.
	 * @param  comm   Communicator.
	 * @param  wRank  Worker process rank.
	 * @param  tag    Message tag.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void receiveTaskOutput
		(T item,
		 Comm comm,
		 int wRank,
		 int tag)
		throws IOException
		{
		}

	/**
	 * Perform per-thread finalization actions after finishing the loop
	 * iterations. Called by a worker thread.
	 * <P>
	 * The <TT>finish()</TT> method may be overridden in a subclass. If not
	 * overridden, the <TT>finish()</TT> method does nothing.
	 *
	 * @exception  Exception
	 *     The <TT>finish()</TT> method may throw any exception.
	 */
	public void finish()
		throws Exception
		{
		}

	/**
	 * Returns the tag offset for this worker for loop. Each message between the
	 * master and worker threads is sent with a message tag equal to
	 * <I>W</I>+<I>T</I>, where <I>W</I> is the worker index and <I>T</I> is the
	 * tag offset.
	 * <P>
	 * The <TT>tagOffset()</TT> method may be overridden in a subclass. If not
	 * overridden, the <TT>tagOffset()</TT> returns a default tag offset of
	 * <TT>Integer.MIN_VALUE</TT>.
	 *
	 * @return  Tag offset.
	 */
	public int tagOffset()
		{
		return Integer.MIN_VALUE;
		}

// Hidden operations.

	/**
	 * Execute this worker iteration in the master thread.
	 *
	 * @param  generator  Item generator.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	void masterExecute
		(ItemGenerator<T> generator)
		throws IOException
		{
		int count = myTeam.count;
		int remaining = count;
		ObjectItemBuf<ItemHolder<T>> buf = ObjectBuf.buffer();
		Range tagRange = new Range (tagFor (0), tagFor (count - 1));
		Comm comm = myTeam.comm;

		// Send initial task to each worker.
		for (int w = 0; w < count; ++ w)
			{
			ItemHolder<T> holder = generator.nextItem();
			buf.item = holder;
			buf.reset();
			int r = myTeam.workerRank (w);
			int tag = tagFor (w);
			comm.send (r, tag, buf);
			if (holder == null)
				{
				-- remaining;
				}
			else
				{
				sendTaskInput (holder.myItem, comm, r, tag);
				}
			}

		// Repeatedly receive a response from a worker and send next task to
		// that worker.
		while (remaining > 0)
			{
			CommStatus status = comm.receive (null, tagRange, buf);
			ItemHolder<T> holder = buf.item;
			int r = status.fromRank;
			int tag = status.tag;
			int w = workerFor (tag);
			receiveTaskOutput (holder.myItem, comm, r, tag);
			holder = generator.nextItem();
			buf.item = holder;
			buf.reset();
			comm.send (r, tag, buf);
			if (holder == null)
				{
				-- remaining;
				}
			else
				{
				sendTaskInput (holder.myItem, comm, r, tag);
				}
			}
		}

	/**
	 * Execute this worker for loop in a worker thread.
	 *
	 * @param  w  Worker index.
	 *
	 * @exception  Exception
	 *     This method may throw any exception.
	 */
	void workerExecute
		(int w)
		throws Exception
		{
		Comm comm = myTeam.comm;
		int r = myTeam.masterRank();
		int tag = tagFor (w);
		start();
		ObjectItemBuf<ItemHolder<T>> buf = ObjectBuf.buffer();
		for (;;)
			{
			comm.receive (r, tag, buf);
			ItemHolder<T> holder = buf.item;
			if (holder == null) break;
			receiveTaskInput (holder.myItem, comm, r, tag);
			run (holder.myItem);
			buf.reset();

			// The next two statements constitute a critical section; other
			// workers in this team must not send messages in between these two
			// messages, or the master can deadlock.
			synchronized (myTeam)
				{
				comm.send (r, tag, buf);
				sendTaskOutput (holder.myItem, comm, r, tag);
				}
			}
		finish();
		}

	/**
	 * Returns the message tag for the given worker index.
	 *
	 * @param  w  Worker index.
	 *
	 * @return  Message tag.
	 */
	private int tagFor
		(int w)
		{
		return w + tagOffset();
		}

	/**
	 * Returns the worker index for the given message tag.
	 *
	 * @param  tag  Message tag.
	 *
	 * @return  Worker index.
	 */
	private int workerFor
		(int tag)
		{
		return tag - tagOffset();
		}

	}
