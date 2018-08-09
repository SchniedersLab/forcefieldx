//******************************************************************************
//
// File:    Timer.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Timer
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
package edu.rit.util;

import java.util.Date;

/**
 * Class Timer controls the execution of a {@linkplain TimerTask}'s timed
 * actions.
 * <P>
 * A timer is in one of three states as shown by this state diagram:
 * <P>
 * <FONT SIZE="-1">
 * <PRE>
 *       +---+                   +---+
 *       |   | stop              |   | start
 *       |   v                   |   v
 *    +---------+    start    +---------+   timeout   +-----------+
 *    |         |------------&gt;|         |------------&gt;|           |
 * --&gt;| STOPPED |             | STARTED |             | TRIGGERED |
 *    |         |&lt;------------|         |&lt;------------|           |
 *    +---------+     stop    +---------+    start    +-----------+
 *       ^   ^                     ^                     |     |
 *       |   |                     |            action() |     | stop
 *       |   |                     |           performed |     |
 *       |   |                     |                     |     |
 *       |   |                     | periodic timeout    |     |
 *       |   |                     +---------------------+     |
 *       |   |                       one-shot timeout    |     |
 *       |   +-------------------------------------------+     |
 *       |                                                     |
 *       +-----------------------------------------------------+
 * </PRE>
 * </FONT>
 * <P>
 * When a timer is created, it is associated with a {@linkplain TimerTask}. The
 * timer is initially <B>stopped.</B> The timer can then be <B>started</B> in
 * various ways. The timer will then <B>time out</B> at some instant, depending
 * on how it was started. When the timer times out, the timer becomes
 * <B>triggered,</B> and the timer's timer task's <TT>action()</TT> method is
 * called. The timer task can stop the timer, or it can start the timer again
 * for a different timeout, or it can leave the timer alone. If the timer task
 * does not start or stop the timer, the timer goes into a state determined by
 * the way the timer was originally started. If the timer was started with a
 * <B>periodic timeout,</B> the timer automatically starts itself to time out
 * again after a certain interval. If the timer was started with a <B>one-shot
 * timeout,</B> the timer automatically stops itself.
 * <P>
 * By calling the appropriate method, a timer may be started to do timeouts in
 * any of the following ways:
 * <UL>
 * <LI>
 * <B>One-shot timeout.</B> The timer task's action is performed just once at a
 * certain point in time, specified either as an absolute time or as an interval
 * from now.
 * <BR>&nbsp;
 * <LI>
 * <B>Periodic fixed-rate timeout.</B> The timer task's action is performed
 * repeatedly at a given interval. The first action happens at a certain point
 * in time, specified either as an absolute time or as an interval from now.
 * Thereafter, each subsequent action starts at the given repetition interval
 * after the <I>scheduled start</I> of the previous action. Thus, the interval
 * between the scheduled starts of successive actions is always the same,
 * regardless of how long each action takes to complete, and the interval
 * between successive actions may vary, depending on how long each action takes
 * to complete.
 * <BR>&nbsp;
 * <LI>
 * <B>Periodic fixed-interval timeout.</B> The timer task's action is performed
 * repeatedly at a given interval. The first action happens at a certain point
 * in time, specified either as an absolute time or as an interval from now.
 * Thereafter, each subsequent action starts at the given repetition interval
 * after the <I>end</I> of the previous action. Thus, the interval between the
 * starts of successive actions may vary, depending on how long each action
 * takes to complete, and the interval between successive actions is always the
 * same, regardless of how long each action takes to complete.
 * </UL>
 * <P>
 * A Timer object is not constructed directly. Rather, a timer object is created
 * by calling a {@linkplain TimerThread}'s <TT>createTimer()</TT> method, giving
 * a timer task to associate with the timer. A single timer thread can control
 * any number of timers. The timer thread is responsible for doing all the
 * timers' timeouts and causing the timers to call the timer tasks'
 * <TT>action()</TT> methods at the proper times. Since a timer thread does not
 * provide real-time guarantees, the time at which a timer task's action
 * actually starts may be later than when the timeout occurred.
 * <P>
 * Classes Timer, {@linkplain TimerTask}, and {@linkplain TimerThread} provide
 * capabilities similar to classes java.util.Timer and java.util.TimerTask.
 * Unlike the latter, they also provide the ability to stop and restart a timer
 * and the ability to deal with race conditions in multithreaded programs.
 *
 * @author Alan Kaminsky
 * @version 15-Dec-2002
 */
public class Timer {

// Hidden data members.
    // Timer thread which created this timer.
    private TimerThread myTimerThread;

    // Associated timer task.
    private TimerTask myTimerTask;

    // State of this timer.
    private int myState = STOPPED;
    private static final int STOPPED = 0;
    private static final int STARTED = 1;
    private static final int TRIGGERED = 2;

    // Kind of timeout.
    private int myKind;
    private static final int ONE_SHOT_TIMEOUT = 0;
    private static final int FIXED_RATE_TIMEOUT = 1;
    private static final int FIXED_INTERVAL_TIMEOUT = 2;

    // Time at which this timer is next scheduled to time out (milliseconds
    // since midnight 01-Jan-1970 UTC).
    private long myTimeout;

    // Interval for periodic timeouts (milliseconds).
    private long myInterval;

// Hidden constructors.
    /**
     * Construct a new timer.
     *
     * @param theTimerThread Timer thread.
     * @param theTimerTask Timer task.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theTimerTask</TT> is null.
     */
    Timer(TimerThread theTimerThread,
            TimerTask theTimerTask) {
        if (theTimerTask == null) {
            throw new NullPointerException();
        }
        myTimerThread = theTimerThread;
        myTimerTask = theTimerTask;
    }

// Exported operations.
    /**
     * Start this timer with a one-shot timeout at the given absolute time. If
     * the time denotes a point in the past, the action is performed
     * immediately.
     *
     * @param theTime Absolute time at which to perform the action.
     */
    public synchronized void start(Date theTime) {
        myState = STARTED;
        myKind = ONE_SHOT_TIMEOUT;
        myTimeout = theTime.getTime();
        myTimerThread.schedule(myTimeout, this);
    }

    /**
     * Start this timer with a one-shot timeout at the given interval from now.
     * If the interval is less than or equal to zero, the action is performed
     * immediately.
     *
     * @param theInterval Timeout interval before performing the action
     * (milliseconds).
     */
    public synchronized void start(long theInterval) {
        myState = STARTED;
        myKind = ONE_SHOT_TIMEOUT;
        myTimeout = System.currentTimeMillis() + theInterval;
        myTimerThread.schedule(myTimeout, this);
    }

    /**
     * Start this timer with a periodic fixed-rate timeout starting at the given
     * absolute time. If the start time denotes a point in the past, the first
     * action is performed immediately. Each subsequent action is started at the
     * given repetition interval after the scheduled start of the previous
     * action.
     *
     * @param theFirstTime Absolute time at which to perform the first action.
     * @param theRepetitionInterval Interval between successive actions
     * (milliseconds).
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>theRepetitionInterval</TT> &lt;= 0.
     */
    public synchronized void start(Date theFirstTime,
            long theRepetitionInterval) {
        if (theRepetitionInterval <= 0) {
            throw new IllegalArgumentException();
        }
        myState = STARTED;
        myKind = FIXED_RATE_TIMEOUT;
        myTimeout = theFirstTime.getTime();
        myInterval = theRepetitionInterval;
        myTimerThread.schedule(myTimeout, this);
    }

    /**
     * Start this timer with a periodic fixed-rate timeout starting at the given
     * interval from now. If the interval is less than or equal to zero, the
     * first action is performed immediately. Each subsequent action is started
     * at the given repetition interval after the scheduled start of the
     * previous action.
     *
     * @param theFirstInterval Timeout interval before performing the first
     * action (milliseconds).
     * @param theRepetitionInterval Interval between successive actions
     * (milliseconds).
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>theRepetitionInterval</TT> &lt;= 0.
     */
    public synchronized void start(long theFirstInterval,
            long theRepetitionInterval) {
        if (theRepetitionInterval <= 0) {
            throw new IllegalArgumentException();
        }
        myState = STARTED;
        myKind = FIXED_RATE_TIMEOUT;
        myTimeout = System.currentTimeMillis() + theFirstInterval;
        myInterval = theRepetitionInterval;
        myTimerThread.schedule(myTimeout, this);
    }

    /**
     * Start this timer with a periodic fixed-interval timeout starting at the
     * given absolute time. If the start time denotes a point in the past, the
     * first action is performed immediately. Each subsequent action is started
     * at the given repetition interval after the end of the previous action.
     *
     * @param theFirstTime Absolute time at which to perform the first action.
     * @param theRepetitionInterval Interval between successive actions
     * (milliseconds).
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>theRepetitionInterval</TT> &lt;= 0.
     */
    public synchronized void startFixedIntervalTimeout(Date theFirstTime,
            long theRepetitionInterval) {
        if (theRepetitionInterval <= 0) {
            throw new IllegalArgumentException();
        }
        myState = STARTED;
        myKind = FIXED_INTERVAL_TIMEOUT;
        myTimeout = theFirstTime.getTime();
        myInterval = theRepetitionInterval;
        myTimerThread.schedule(myTimeout, this);
    }

    /**
     * Start this timer with a periodic fixed-interval timeout starting at the
     * given interval from now. If the interval is less than or equal to zero,
     * the first action is performed immediately. Each subsequent action is
     * started at the given repetition interval after the end of the previous
     * action.
     *
     * @param theFirstInterval Timeout interval before performing the first
     * action (milliseconds).
     * @param theRepetitionInterval Interval between successive actions
     * (milliseconds).
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>theRepetitionInterval</TT> &lt;= 0.
     */
    public synchronized void startFixedIntervalTimeout(long theFirstInterval,
            long theRepetitionInterval) {
        if (theRepetitionInterval <= 0) {
            throw new IllegalArgumentException();
        }
        myState = STARTED;
        myKind = FIXED_INTERVAL_TIMEOUT;
        myTimeout = System.currentTimeMillis() + theFirstInterval;
        myInterval = theRepetitionInterval;
        myTimerThread.schedule(myTimeout, this);
    }

    /**
     * Stop this timer.
     */
    public synchronized void stop() {
        myState = STOPPED;
    }

    /**
     * Determine whether this timer is stopped.
     *
     * @return True if this timer is in the stopped state, false otherwise.
     */
    public synchronized boolean isStopped() {
        return myState == STOPPED;
    }

    /**
     * Determine whether this timer is started.
     *
     * @return True if this timer is in the started state, false otherwise.
     */
    public synchronized boolean isStarted() {
        return myState == STARTED;
    }

    /**
     * Determine whether this timer is triggered.
     *
     * @return True if this timer is in the triggered state, false otherwise.
     */
    public synchronized boolean isTriggered() {
        return myState == TRIGGERED;
    }

    /**
     * Determine the time when this timer is or was scheduled to time out.
     * <P>
     * If the <TT>getTimeout()</TT> method is called when this timer is in the
     * stopped state, then <TT>Long.MAX_VALUE</TT> is returned.
     * <P>
     * If the <TT>getTimeout()</TT> method is called when this timer is in the
     * started state, then the <TT>getTimeout()</TT> method returns the time
     * when the next timeout will occur.
     * <P>
     * If the <TT>getTimeout()</TT> method is called when this timer is in the
     * triggered state, such as by this timer's timer task's <TT>action()</TT>
     * method, then the <TT>getTimeout()</TT> method returns the time when the
     * present timeout was <I>scheduled</I> to occur. Since a timer thread
     * provides no real-time guarantees, the time when the timeout
     * <I>actually</I> occurred may be some time later. If desired, the caller
     * can compare the scheduled timeout to the actual time and decide whether
     * it's too late to perform the action.
     *
     * @return Time of scheduled timeout (milliseconds since midnight
     * 01-Jan-1970 UTC).
     */
    public synchronized long getTimeout() {
        return myState == STOPPED ? Long.MAX_VALUE : myTimeout;
    }

    /**
     * Returns this timer's timer task.
     *
     * @return a {@link edu.rit.util.TimerTask} object.
     */
    public TimerTask getTimerTask() {
        return myTimerTask;
    }

// Hidden operations.
    /**
     * Trigger this timer.
     *
     * @param theTriggerTime Time at which the trigger occurred (milliseconds
     * since midnight 01-Jan-1970 UTC).
     */
    void trigger(long theTriggerTime) {
        synchronized (this) {
            // Make sure we're started and it really is time to trigger.
            if (myState != STARTED || myTimeout > theTriggerTime) {
                return;
            }

            // Switch to the triggered state.
            myState = TRIGGERED;
        }

        // Call the timer task's action() method. Do it outside the synchronized
        // block, otherwise a deadlock may happen if another thread tries to
        // start or stop the timer.
        myTimerTask.action(this);

        synchronized (this) {
            // Decide whether to do an automatic restart.
            if (myState != TRIGGERED) {
                // Someone already stopped or restarted us. Do nothing.
            } else if (myKind == FIXED_RATE_TIMEOUT) {
                myState = STARTED;
                myTimeout += myInterval;
                myTimerThread.schedule(myTimeout, this);
            } else if (myKind == FIXED_INTERVAL_TIMEOUT) {
                myState = STARTED;
                myTimeout = System.currentTimeMillis() + myInterval;
                myTimerThread.schedule(myTimeout, this);
            } else // ONE_SHOT_TIMEOUT
            {
                myState = STOPPED;
            }
        }
    }

}
