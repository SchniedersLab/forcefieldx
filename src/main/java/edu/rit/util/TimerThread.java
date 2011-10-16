//******************************************************************************
//
// File:    TimerThread.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.TimerThread
//
// This Java source file is copyright (C) 2002-2004 by Alan Kaminsky. All rights
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

package edu.rit.util;

import java.util.Iterator;
import java.util.Vector;

/**
 * Class TimerThread encapsulates a thread that does the timing for {@linkplain
 * Timer}s and performs {@linkplain TimerTask}s' actions when timeouts occur.
 * <P>
 * A timer is created by calling a timer thread's <TT>createTimer()</TT> method,
 * giving a timer task to associate with the timer. Multiple timers may be
 * created under the control of the same timer thread. The timer thread will
 * perform the actions of its timers' timer tasks one at a time at the proper
 * instants (by causing each timer to call its timer task's <TT>action()</TT>
 * method). Note that the timer tasks of a single timer thread are performed
 * sequentially, which may cause some timer tasks' actions to be delayed
 * depending on how long other timer tasks' actions take to execute. If
 * necessary, consider running separate timers in separate timer threads so the
 * timer tasks' actions run concurrently.
 * <P>
 * Class TimerThread does <I>not</I> offer real-time guarantees. It merely makes
 * a best-effort attempt to perform each timer task's actions as soon as
 * possible after the timeouts occur.
 * <P>
 * A TimerThread is just a {@linkplain java.lang.Thread Thread}. After
 * constructing a timer thread, you can mark it as a daemon thread if you want.
 * You must also call the timer thread's <TT>start()</TT> method, or no timeouts
 * will occur. To gracefully stop a timer thread, call the <TT>shutdown()</TT>
 * method.
 * <P>
 * To simplify writing programs with multiple objects that all use the same
 * timer thread, the static <TT>TimerThread.getDefault()</TT> method returns a
 * single shared instance of class TimerThread. The default timer thread is
 * marked as a daemon thread and is started automatically. The default timer
 * thread is not created until the first call to
 * <TT>TimerThread.getDefault()</TT>.
 * <P>
 * Classes {@linkplain Timer}, {@linkplain TimerTask}, and TimerThread provide
 * capabilities similar to classes java.util.Timer and java.util.TimerTask.
 * Unlike the latter, they also provide the ability to stop and restart a timer
 * and the ability to deal with race conditions in multithreaded programs.
 *
 * @author  Alan Kaminsky
 * @version 18-Jun-2003
 */
public class TimerThread
	extends Thread
	{

// Hidden helper classes.

	/**
	 * Class TimerThread.TimeoutInfo is a record used to keep track of a
	 * timeout.
	 *
	 * @author  Alan Kaminsky
	 * @version 18-Sep-2002
	 */
	private static class TimeoutInfo
		{
		// Instant at which the timeout will occur (milliseconds since midnight
		// 01-Jan-1970 UTC).
		public long myTimeout;

		// Timer to trigger.
		public Timer myTimer;

		public TimeoutInfo
			(long theTimeout,
			 Timer theTimer)
			{
			myTimeout = theTimeout;
			myTimer = theTimer;
			}
		}

// Hidden data members.

	/**
	 * Queue of timeout info records, organized as a heap in order of timeout.
	 * Index 0 is unused. Indexes 1 .. mySize contain the heap. The queue's
	 * length is increased by INCR when necessary to add a new record.
	 */
	private static final int INCR = 4;
	private TimeoutInfo[] myQueue = new TimeoutInfo [INCR+1];

	/**
	 * Number of records in the queue.
	 */
	private int mySize = 0;

	/**
	 * True if this timer thread is running, false if it's shut down.
	 */
	private boolean iamRunning = true;

// Hidden static data members.

	/**
	 * The default timer thread.
	 */
	private static TimerThread theDefaultTimerThread = null;

// Exported constructors.

	/**
	 * Construct a new timer thread. After constructing it, you must call the
	 * timer thread's <TT>start()</TT> method, or no timeouts will occur.
	 */
	public TimerThread()
		{
		super();
		}

// Exported operations.

	/**
	 * Get the default timer thread, a single shared instance of class
	 * TimerThread. The default timer thread is marked as a daemon thread and is
	 * started automatically. The default timer thread is not created until the
	 * first call to <TT>TimerThread.getDefault()</TT>.
	 *
	 * @return  Default timer thread.
	 */
	public static synchronized TimerThread getDefault()
		{
		if (theDefaultTimerThread == null)
			{
			theDefaultTimerThread = new TimerThread();
			theDefaultTimerThread.setDaemon (true);
			theDefaultTimerThread.start();
			}
		return theDefaultTimerThread;
		}

	/**
	 * Create a new timer associated with the given timer task and under the
	 * control of this timer thread. When the timer is triggered, this timer
	 * thread will cause the timer to call the given timer task's
	 * <TT>action()</TT> method.
	 *
	 * @param  theTimerTask  Timer task.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theTimerTask</TT> is null.
	 */
	public Timer createTimer
		(TimerTask theTimerTask)
		{
		return new Timer (this, theTimerTask);
		}

	/**
	 * Shut down this timer thread.
	 */
	public void shutdown()
		{
		synchronized (this)
			{
			iamRunning = false;
			notifyAll();
			}
		}

	/**
	 * Perform this timer thread's processing. (Never call the <TT>run()</TT>
	 * method yourself!)
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if some thread other than this timer
	 *     thread called the <TT>run()</TT> method.
	 */
	public void run()
		{
		// Only this timer thread itself can call the run() method.
		if (Thread.currentThread() != this)
			{
			throw new IllegalStateException
				("Wrong thread called the run() method");
			}

		try
			{
			while (iamRunning)
				{
				long now = System.currentTimeMillis();
				Vector theTriggeredTimeouts = null;

				synchronized (this)
					{
					// If timeout queue is empty, wait until notified.
					if (mySize == 0)
						{
						wait();
						now = System.currentTimeMillis();
						}

					// If timeout queue is not empty and first timeout is in the
					// future, wait until timeout or until notified.
					else
						{
						long waitTime = myQueue[1].myTimeout - now;
						if (waitTime > 0L)
							{
							wait (waitTime);
							now = System.currentTimeMillis();
							}
						}

					// Pull all timeouts that have occurred out of the timeout
					// queue into a separate list.
					theTriggeredTimeouts = new Vector();
					while (mySize > 0 && myQueue[1].myTimeout <= now)
						{
						theTriggeredTimeouts.add (myQueue[1]);
						myQueue[1] = myQueue[mySize];
						myQueue[mySize] = null;
						-- mySize;
						siftDown (mySize);
						}
					}

				// Perform the action of each triggered timeout. Do this outside
				// the synchronized block, or a deadlock may happen if a timer
				// is restarted.
				Iterator iter = theTriggeredTimeouts.iterator();
				while (iter.hasNext())
					{
					((TimeoutInfo) iter.next()).myTimer.trigger (now);
					}
				}
			}

		catch (InterruptedException exc)
			{
			System.err.println ("TimerThread interrupted");
			exc.printStackTrace (System.err);
			}
		}

// Hidden operations.

	/**
	 * Schedule a timer.
	 *
	 * @param  theTimeout  Timeout instant.
	 * @param  theTimer    Timer.
	 */
	synchronized void schedule
		(long theTimeout,
		 Timer theTimer)
		{
		// Increase queue allocation if necessary.
		if (mySize == myQueue.length - 1)
			{
			TimeoutInfo[] newQueue = new TimeoutInfo [myQueue.length + INCR];
			System.arraycopy (myQueue, 1, newQueue, 1, mySize);
			myQueue = newQueue;
			}

		// Add a new timeout info record to the queue.
		++ mySize;
		myQueue[mySize] = new TimeoutInfo (theTimeout, theTimer);
		siftUp (mySize);

		// Wake up the timer thread.
		notifyAll();
		}

	/**
	 * Sift up the last element in the heap. Precondition: HEAP(1,n-1) and n
	 * &gt; 0. Postcondition: HEAP(1,n).
	 */
	private void siftUp
		(int n)
		{
		int i = n;
		int p;
		TimeoutInfo temp;
		for (;;)
			{
			// Invariant: HEAP(1,n) except perhaps between i and its parent.
			if (i == 1) break;
			p = i / 2;
			if (myQueue[p].myTimeout <= myQueue[i].myTimeout) break;
			temp = myQueue[p];
			myQueue[p] = myQueue[i];
			myQueue[i] = temp;
			i = p;
			}
		}

	/**
	 * Sift down the first element in the heap. Precondition: HEAP(2,n) and n
	 * &gt;= 0. Postcondition: HEAP(1,n).
	 */
	private void siftDown
		(int n)
		{
		int i = 1;
		int c;
		TimeoutInfo temp;
		for (;;)
			{
			// Invariant: HEAP(1,n) except perhaps between i and its 0, 1, or 2
			// children.
			c = 2 * i;
			if (c > n) break;
			// c is the left child of i.
			if (c+1 <= n)
				{
				// c+1 is the right child of i.
				if (myQueue[c+1].myTimeout < myQueue[c].myTimeout)
					{
					c = c + 1;
					}
				}
			// c is the least child of i.
			if (myQueue[i].myTimeout <= myQueue[c].myTimeout) break;
			temp = myQueue[c];
			myQueue[c] = myQueue[i];
			myQueue[i] = temp;
			i = c;
			}
		}

// Unit test main program.

//	/**
//	 * Helper class for unit test main program.
//	 */
//	private static class Message
//		implements TimerTask
//		{
//		private String myMessage;
//
//		public Message
//			(String theMessage)
//			{
//			myMessage = theMessage;
//			}
//
//		public void action
//			(Timer theTimer)
//			{
//			System.out.println (myMessage);
//			}
//		}
//
//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		{
//		try
//			{
//			TimerThread theTimerThread = new TimerThread();
//			theTimerThread.start();
//
//			Timer timer1 =
//				theTimerThread.createTimer
//					(new Message ("Timer 1 triggered"));
//
//			Timer timer2 =
//				theTimerThread.createTimer
//					(new Message ("Timer 2 triggered"));
//
//			Timer timer3 =
//				theTimerThread.createTimer
//					(new Message ("Timer 3 triggered"));
//
//			timer1.start (1000L, 1000L);
//			timer2.start (5000L, 5000L);
//			timer3.start (10000L);
//			}
//		catch (Throwable exc)
//			{
//			System.err.println ("edu.rit.util.TimerThread: Uncaught exception");
//			exc.printStackTrace (System.err);
//			}
//		}

	}
