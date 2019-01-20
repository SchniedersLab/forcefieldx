//******************************************************************************
//
// File:    ParallelRegion.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelRegion
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
package edu.rit.pj;

import java.util.Iterator;

import edu.rit.util.LongRange;
import edu.rit.util.Range;

/**
 * Class ParallelRegion is the abstract base class for a parallel region that is
 * executed by a {@linkplain ParallelTeam} of threads.
 * <P>
 * To execute a parallel region, create a {@linkplain ParallelTeam} object;
 * create an instance of a concrete subclass of class ParallelRegion; and pass
 * this instance to the parallel team's <code>execute()</code> method. You can do
 * all this using an anonymous inner class; for example:
 * <PRE>
 *     new ParallelTeam().execute (new ParallelRegion()
 *         {
 *         // Shared variable declarations
 *         . . .
 *         public void start()
 *             {
 *             // Initialization code
 *             . . .
 *             }
 *         public void run()
 *             {
 *             // Thread private variable declarations
 *             // Parallel code
 *             . . .
 *             }
 *         public void finish()
 *             {
 *             // Finalization code
 *             . . .
 *             }
 *         });
 * </PRE>
 * <P>
 * The parallel team's <code>execute()</code> method does the following. The
 * parallel team has a certain number of threads <I>K</I>, where <I>K</I> was
 * specified when the parallel team was constructed. The main thread is the
 * thread calling the parallel team's <code>execute()</code> method. The main thread
 * calls the parallel region's <code>start()</code> method. When the
 * <code>start()</code> method returns, all the team threads call the parallel
 * region's <code>run()</code> method concurrently. When all the team threads have
 * returned from the <code>run()</code> method, the main thread calls the parallel
 * region's <code>finish()</code> method. When the <code>finish()</code> method returns,
 * the main thread returns from the parallel team's <code>execute()</code> method.
 * <P>
 * Variables to be shared by all threads in the team may be declared as fields
 * of the ParallelRegion subclass. The <code>start()</code> method is intended for
 * performing initialization in a single thread before parallel execution
 * begins. If no such initialization is needed, omit the <code>start()</code>
 * method. The <code>run()</code> method contains code to be executed in parallel by
 * all threads in the team. Variables that are private to each thread may be
 * declared inside the <code>run()</code> method. The <code>finish()</code> method is
 * intended for performing finalization in a single thread after parallel
 * execution ends. If no such finalization is needed, omit the <code>finish()</code>
 * method.
 * <P>
 * If the parallel region's <code>start()</code> method throws an exception, the
 * parallel team's <code>execute()</code> method throws that same exception, and the
 * <code>run()</code> method is not called.
 * <P>
 * If the parallel region's <code>run()</code> method throws an exception in one of
 * the team threads, the exception's stack trace is printed on the standard
 * error, the parallel team waits until all the other team threads have returned
 * from the <code>run()</code> method, then the parallel team's <code>execute()</code>
 * method throws that same exception, and the parallel region's
 * <code>finish()</code> method is not called. If the parallel region's
 * <code>run()</code> method throws an exception in more than one of the team
 * threads, each exception's stack trace is printed on the standard error, the
 * parallel team waits until all the other team threads have returned from the
 * <code>run()</code> method, then the parallel team's <code>execute()</code> method
 * throws a {@linkplain MultipleParallelException} wrapping all the thrown
 * exceptions, and the parallel region's <code>finish()</code> method is not called.
 * <P>
 * If the parallel region's <code>finish()</code> method throws an exception, the
 * parallel team's <code>execute()</code> method throws that same exception.
 *
 * @author Alan Kaminsky
 * @version 11-Nov-2007
 */
public abstract class ParallelRegion
        extends ParallelConstruct {

// Hidden data members.
    // Default lock for critical() and criticalNonexclusive() methods.
    private Lock myLock = new Lock();

// Exported constructors.
    /**
     * Construct a new parallel region.
     */
    public ParallelRegion() {
        super();
    }

// Exported operations.
    /**
     * Perform initialization actions before parallel execution begins. Only one
     * thread calls the <code>start()</code> method.
     * <P>
     * The <code>start()</code> method may be overridden in a subclass. If not
     * overridden, the <code>start()</code> method does nothing.
     *
     * @exception Exception The <code>start()</code> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public void start()
            throws Exception {
    }

    /**
     * Execute parallel code. All threads of the parallel team call the
     * <code>run()</code> method concurrently.
     * <P>
     * The <code>run()</code> method must be implemented in a subclass.
     *
     * @exception Exception The <code>run()</code> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public abstract void run()
            throws Exception;

    /**
     * Perform finalization actions after parallel execution ends. Only one
     * thread calls the <code>finish()</code> method.
     * <P>
     * The <code>finish()</code> method may be overridden in a subclass. If not
     * overridden, the <code>finish()</code> method does nothing.
     *
     * @exception Exception The <code>finish()</code> method may throw any
     * exception.
     * @throws java.lang.Exception if any.
     */
    public void finish()
            throws Exception {
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain IntegerForLoop}. The loop index goes
     * from <code>first</code> (inclusive) to <code>last</code> (inclusive) in steps of
     * +1. If <code>first</code> is greater than <code>last</code>, then no loop
     * iterations are performed. At the end of the parallel for loop, the
     * parallel team threads wait for each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param theLoop Parallel for loop.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(int first,
            int last,
            IntegerForLoop theLoop)
            throws Exception {
        execute(first, last, theLoop, BarrierAction.WAIT);
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain IntegerForLoop}. The loop index goes
     * from <code>first</code> (inclusive) to <code>last</code> (inclusive) in steps of
     * +1. If <code>first</code> is greater than <code>last</code>, then no loop
     * iterations are performed. At the end of the parallel for loop, the
     * parallel team threads encounter a barrier, and their behavior depends on
     * the given {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param theLoop Parallel for loop.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null. Thrown if
     * <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(int first,
            int last,
            IntegerForLoop theLoop,
            BarrierAction action)
            throws Exception {
        // Verify preconditions.
        if (theLoop == null) {
            throw new NullPointerException("ParallelRegion.execute(): Parallel for loop is null");
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            ParallelTeamThread currentThread = getCurrentThread();
            int currentIndex = currentThread.myIndex;

            // Do top-of-parallel-construct processing.
            IntegerSchedule schedule = null;
            if (currentThread.arriveAtParallelConstruct()) {
                // First thread to arrive sets up the shared parallel for loop
                // schedule object and stores it (or an exception if any) in
                // each team thread.
                try {
                    schedule = theLoop.schedule();
                    schedule.commonStart(myTeam.K, new Range(first, last));
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setIntegerSchedule(schedule);
                    }
                } catch (Throwable exc) {
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setConstructException(exc);
                    }
                }
            }

            // Get the shared parallel for loop schedule object.
            schedule = currentThread.getIntegerSchedule();
            theLoop.mySchedule = schedule;

            // Prepare to catch exceptions thrown by the parallel for loop body.
            Throwable runException = null;
            try {
                // Perform per-thread initialization.
                theLoop.start();

                // Repeatedly get and process a chunk of loop iterations.
                Range chunk;
                while ((chunk = schedule.commonNext(currentIndex)) != null) {
                    theLoop.commonRun(chunk.lb(), chunk.ub());
                }

                // Perform per-thread finalization.
                theLoop.finish();
            } catch (Throwable exc) {
                runException = exc;
                schedule.myBreak = true;
            }

            // Barrier synchronization.
            action.doBarrier(currentThread);

            // Propagate any exception thrown by the run() method.
            ParallelTeam.rethrow(runException);
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
            theLoop.mySchedule = null;
        }
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain IntegerStrideForLoop}. The loop index
     * goes from <code>first</code> (inclusive) to <code>last</code> (inclusive) in
     * steps of <code>stride</code>. The stride must be positive. If <code>first</code>
     * is greater than <code>last</code>, then no loop iterations are performed. At
     * the end of the parallel for loop, the parallel team threads wait for each
     * other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param stride Loop index stride, &gt;= 1.
     * @param theLoop Parallel for loop.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>stride</code> &lt; 1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(int first,
            int last,
            int stride,
            IntegerStrideForLoop theLoop)
            throws Exception {
        execute(first, last, stride, theLoop, BarrierAction.WAIT);
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain IntegerStrideForLoop}. The loop index
     * goes from <code>first</code> (inclusive) to <code>last</code> (inclusive) in
     * steps of <code>stride</code>. The stride must be positive. If <code>first</code>
     * is greater than <code>last</code>, then no loop iterations are performed. At
     * the end of the parallel for loop, the parallel team threads encounter a
     * barrier, and their behavior depends on the given {@linkplain
     * BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param stride Loop index stride, &gt;= 1.
     * @param theLoop Parallel for loop.
     * @param action Barrier action.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>stride</code> &lt; 1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null. Thrown if
     * <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(int first,
            int last,
            int stride,
            IntegerStrideForLoop theLoop,
            BarrierAction action)
            throws Exception {
        // Verify preconditions.
        if (stride <= 0) {
            throw new IllegalArgumentException("ParallelRegion.execute(): Stride = " + stride + " illegal");
        }
        if (theLoop == null) {
            throw new NullPointerException("ParallelRegion.execute(): Parallel for loop is null");
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            ParallelTeamThread currentThread = getCurrentThread();
            int currentIndex = currentThread.myIndex;

            // Do top-of-parallel-construct processing.
            IntegerSchedule schedule = null;
            if (currentThread.arriveAtParallelConstruct()) {
                // First thread to arrive sets up the shared parallel for loop
                // schedule object and stores it (or an exception if any) in
                // each team thread.
                try {
                    schedule = theLoop.schedule();
                    schedule.commonStart(myTeam.K, new Range(first, last, stride));
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setIntegerSchedule(schedule);
                    }
                } catch (Throwable exc) {
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setConstructException(exc);
                    }
                }
            }

            // Get the shared parallel for loop schedule object.
            schedule = currentThread.getIntegerSchedule();
            theLoop.mySchedule = schedule;

            // Prepare to catch exceptions thrown by the parallel for loop body.
            Throwable runException = null;
            try {
                // Perform per-thread initialization.
                theLoop.start();

                // Repeatedly get and process a chunk of loop iterations.
                Range chunk;
                while ((chunk = schedule.commonNext(currentIndex)) != null) {
                    theLoop.commonRun(chunk.lb(), chunk.ub(), chunk.stride());
                }

                // Perform per-thread finalization.
                theLoop.finish();
            } catch (Throwable exc) {
                runException = exc;
                schedule.myBreak = true;
            }

            // Barrier synchronization.
            action.doBarrier(currentThread);

            // Propagate any exception thrown by the run() method.
            ParallelTeam.rethrow(runException);
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
            theLoop.mySchedule = null;
        }
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain LongForLoop}. The loop index goes from
     * <code>first</code> (inclusive) to <code>last</code> (inclusive) in steps of +1.
     * If <code>first</code> is greater than <code>last</code>, then no loop iterations
     * are performed. At the end of the parallel for loop, the parallel team
     * threads wait for each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param theLoop Parallel for loop.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(long first,
            long last,
            LongForLoop theLoop)
            throws Exception {
        execute(first, last, theLoop, BarrierAction.WAIT);
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain LongForLoop}. The loop index goes from
     * <code>first</code> (inclusive) to <code>last</code> (inclusive) in steps of +1.
     * If <code>first</code> is greater than <code>last</code>, then no loop iterations
     * are performed. At the end of the parallel for loop, the parallel team
     * threads encounter a barrier, and their behavior depends on the given
     * {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param theLoop Parallel for loop.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null. Thrown if
     * <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(long first,
            long last,
            LongForLoop theLoop,
            BarrierAction action)
            throws Exception {
        // Verify preconditions.
        if (theLoop == null) {
            throw new NullPointerException("ParallelRegion.execute(): Parallel for loop is null");
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            ParallelTeamThread currentThread = getCurrentThread();
            int currentIndex = currentThread.myIndex;

            // Do top-of-parallel-construct processing.
            LongSchedule schedule = null;
            if (currentThread.arriveAtParallelConstruct()) {
                // First thread to arrive sets up the shared parallel for loop
                // schedule object and stores it (or an exception if any) in
                // each team thread.
                try {
                    schedule = theLoop.schedule();
                    schedule.commonStart(myTeam.K, new LongRange(first, last));
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setLongSchedule(schedule);
                    }
                } catch (Throwable exc) {
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setConstructException(exc);
                    }
                }
            }

            // Get the shared parallel for loop schedule object.
            schedule = currentThread.getLongSchedule();
            theLoop.mySchedule = schedule;

            // Prepare to catch exceptions thrown by the parallel for loop body.
            Throwable runException = null;
            try {
                // Perform per-thread initialization.
                theLoop.start();

                // Repeatedly get and process a chunk of loop iterations.
                LongRange chunk;
                while ((chunk = schedule.commonNext(currentIndex)) != null) {
                    theLoop.commonRun(chunk.lb(), chunk.ub());
                }

                // Perform per-thread finalization.
                theLoop.finish();
            } catch (Throwable exc) {
                runException = exc;
                schedule.myBreak = true;
            }

            // Barrier synchronization.
            action.doBarrier(currentThread);

            // Propagate any exception thrown by the run() method.
            ParallelTeam.rethrow(runException);
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
            theLoop.mySchedule = null;
        }
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain LongStrideForLoop}. The loop index
     * goes from <code>first</code> (inclusive) to <code>last</code> (inclusive) in
     * steps of <code>stride</code>. The stride must be positive. If <code>first</code>
     * is greater than <code>last</code>, then no loop iterations are performed. At
     * the end of the parallel for loop, the parallel team threads wait for each
     * other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param stride Loop index stride, &gt;= 1.
     * @param theLoop Parallel for loop.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>stride</code> &lt; 1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(long first,
            long last,
            long stride,
            LongStrideForLoop theLoop)
            throws Exception {
        execute(first, last, stride, theLoop, BarrierAction.WAIT);
    }

    /**
     * Execute a parallel for loop within this parallel region. For further
     * information, see class {@linkplain LongStrideForLoop}. The loop index
     * goes from <code>first</code> (inclusive) to <code>last</code> (inclusive) in
     * steps of <code>stride</code>. The stride must be positive. If <code>first</code>
     * is greater than <code>last</code>, then no loop iterations are performed. At
     * the end of the parallel for loop, the parallel team threads encounter a
     * barrier, and their behavior depends on the given {@linkplain
     * BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param stride Loop index stride, &gt;= 1.
     * @param theLoop Parallel for loop.
     * @param action Barrier action.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <code>stride</code> &lt; 1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLoop</code> is null. Thrown if
     * <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theLoop</code>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(long first,
            long last,
            long stride,
            LongStrideForLoop theLoop,
            BarrierAction action)
            throws Exception {
        // Verify preconditions.
        if (stride <= 0) {
            throw new IllegalArgumentException("ParallelRegion.execute(): Stride = " + stride + " illegal");
        }
        if (theLoop == null) {
            throw new NullPointerException("ParallelRegion.execute(): Parallel for loop is null");
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            ParallelTeamThread currentThread = getCurrentThread();
            int currentIndex = currentThread.myIndex;

            // Do top-of-parallel-construct processing.
            LongSchedule schedule = null;
            if (currentThread.arriveAtParallelConstruct()) {
                // First thread to arrive sets up the shared parallel for loop
                // schedule object and stores it (or an exception if any) in
                // each team thread.
                try {
                    schedule = theLoop.schedule();
                    schedule.commonStart(myTeam.K, new LongRange(first, last, stride));
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setLongSchedule(schedule);
                    }
                } catch (Throwable exc) {
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setConstructException(exc);
                    }
                }
            }

            // Get the shared parallel for loop schedule object.
            schedule = currentThread.getLongSchedule();
            theLoop.mySchedule = schedule;

            // Prepare to catch exceptions thrown by the parallel for loop body.
            Throwable runException = null;
            try {
                // Perform per-thread initialization.
                theLoop.start();

                // Repeatedly get and process a chunk of loop iterations.
                LongRange chunk;
                while ((chunk = schedule.commonNext(currentIndex)) != null) {
                    theLoop.commonRun(chunk.lb(), chunk.ub(), chunk.stride());
                }

                // Perform per-thread finalization.
                theLoop.finish();
            } catch (Throwable exc) {
                runException = exc;
                schedule.myBreak = true;
            }

            // Barrier synchronization.
            action.doBarrier(currentThread);

            // Propagate any exception thrown by the run() method.
            ParallelTeam.rethrow(runException);
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
            theLoop.mySchedule = null;
        }
    }

    /**
     * Execute a parallel iteration within this parallel region. For further
     * information, see class {@linkplain ParallelIteration}. The items
     * processed by the iteration are the elements of the given array. The
     * iteration order is from index 0 upwards. At the end of the parallel
     * iteration, the parallel team threads wait for each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theArray Array containing the items.
     * @param theIteration Parallel iteration.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null or
     * <code>theIteration</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theIteration</code>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(T[] theArray,
            ParallelIteration<T> theIteration)
            throws Exception {
        execute(theArray, theIteration, BarrierAction.WAIT);
    }

    /**
     * Suppress warnings for casts of ItemGenerator.
     * @param obj The ItemGenerator instance.
     * @param <T> Data type of the items iterated over.
     * @return
     */
    @SuppressWarnings("unchecked")
    private static <T> ItemGenerator<T> castItemGenerator(Object obj) {
        return (ItemGenerator<T>) obj;
    }

    /**
     * Execute a parallel iteration within this parallel region. For further
     * information, see class {@linkplain ParallelIteration}. The items
     * processed by the iteration are the elements of the given array. The
     * iteration order is from index 0 upwards. At the end of the parallel
     * iteration, the parallel team threads encounter a barrier, and their
     * behavior depends on the given {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theArray Array containing the items.
     * @param theIteration Parallel iteration.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theArray</code> is null. Thrown if
     * <code>theIteration</code> is null. Thrown if <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theIteration</code>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(T[] theArray,
            ParallelIteration<T> theIteration,
            BarrierAction action)
            throws Exception {
        // Verify preconditions.
        if (theArray == null) {
            throw new NullPointerException("ParallelRegion.execute(): Array is null");
        }
        if (theIteration == null) {
            throw new NullPointerException("ParallelRegion.execute(): Parallel iteration is null");
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theIteration.myTeam = this.myTeam;

            // Get current parallel team thread.
            ParallelTeamThread currentThread = getCurrentThread();

            // Do top-of-parallel-construct processing.
            ItemGenerator<T> generator = null;
            if (currentThread.arriveAtParallelConstruct()) {
                // First thread to arrive sets up the shared item generator
                // object and stores it (or an exception if any) in each team
                // thread.
                try {
                    generator = new ArrayItemGenerator<T>(theArray);
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setItemGenerator(generator);
                    }
                } catch (Throwable exc) {
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setConstructException(exc);
                    }
                }
            }

            // Get the shared item generator object.
            generator = castItemGenerator(currentThread.getItemGenerator());
            theIteration.myItemGenerator = generator;

            // Prepare to catch exceptions thrown by the parallel iteration
            // body.
            Throwable runException = null;
            try {
                // Perform per-thread initialization.
                theIteration.start();

                // Repeatedly get and process an item.
                ItemHolder<T> itemholder;
                while ((itemholder = generator.nextItem()) != null) {
                    theIteration.commonRun(itemholder.mySequenceNumber, itemholder.myItem);
                }

                // Perform per-thread finalization.
                theIteration.finish();
            } catch (Throwable exc) {
                runException = exc;
                generator.myBreak = true;
            }

            // Barrier synchronization.
            action.doBarrier(currentThread);

            // Propagate any exception thrown by the run() method.
            ParallelTeam.rethrow(runException);
        } finally {
            // Forget parallel team.
            theIteration.myTeam = null;
            theIteration.myItemGenerator = null;
        }
    }

    /**
     * Execute a parallel iteration within this parallel region. For further
     * information, see class {@linkplain ParallelIteration}. The items
     * processed by the iteration are the items returned by the given iterator.
     * The iteration order is that of the given iterator. At the end of the
     * parallel iteration, the parallel team threads wait for each other at a
     * barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theIterator Iterator over the items.
     * @param theIteration Parallel iteration.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theIterator</code> is null or
     * <code>theIteration</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theIteration</code>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(Iterator<T> theIterator,
            ParallelIteration<T> theIteration)
            throws Exception {
        execute(theIterator, theIteration, BarrierAction.WAIT);
    }

    /**
     * Execute a parallel iteration within this parallel region. For further
     * information, see class {@linkplain ParallelIteration}. The items
     * processed by the iteration are the items returned by the given iterator.
     * The iteration order is that of the given iterator. At the end of the
     * parallel iteration, the parallel team threads encounter a barrier, and
     * their behavior depends on the given {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theIterator Iterator over the items.
     * @param theIteration Parallel iteration.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theIterator</code> is null. Thrown if <code>theIteration</code> is null.
     * Thrown if <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theIteration</code>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(Iterator<T> theIterator,
            ParallelIteration<T> theIteration,
            BarrierAction action)
            throws Exception {
        // Verify preconditions.
        if (theIterator == null) {
            throw new NullPointerException("ParallelRegion.execute(): Iterator is null");
        }
        if (theIteration == null) {
            throw new NullPointerException("ParallelRegion.execute(): Parallel iteration is null");
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theIteration.myTeam = this.myTeam;

            // Get current parallel team thread.
            ParallelTeamThread currentThread = getCurrentThread();

            // Do top-of-parallel-construct processing.
            ItemGenerator<T> generator = null;
            if (currentThread.arriveAtParallelConstruct()) {
                // First thread to arrive sets up the shared item generator
                // object and stores it (or an exception if any) in each team
                // thread.
                try {
                    generator = new IteratorItemGenerator<T>(theIterator);
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setItemGenerator(generator);
                    }
                } catch (Throwable exc) {
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setConstructException(exc);
                    }
                }
            }

            // Get the shared item generator object.
            generator = castItemGenerator(currentThread.getItemGenerator());
            theIteration.myItemGenerator = generator;

            // Prepare to catch exceptions thrown by the parallel iteration
            // body.
            Throwable runException = null;
            try {
                // Perform per-thread initialization.
                theIteration.start();

                // Repeatedly get and process an item.
                ItemHolder<T> itemholder;
                while ((itemholder = generator.nextItem()) != null) {
                    theIteration.commonRun(itemholder.mySequenceNumber, itemholder.myItem);
                }

                // Perform per-thread finalization.
                theIteration.finish();
            } catch (Throwable exc) {
                runException = exc;
                generator.myBreak = true;
            }

            // Barrier synchronization.
            action.doBarrier(currentThread);

            // Propagate any exception thrown by the run() method.
            ParallelTeam.rethrow(runException);
        } finally {
            // Forget parallel team.
            theIteration.myTeam = null;
            theIteration.myItemGenerator = null;
        }
    }

    /**
     * Execute a parallel iteration within this parallel region. For further
     * information, see class {@linkplain ParallelIteration}. The items
     * processed by the iteration are the items contained in the given iterable
     * collection. The iteration order is that of the given iterable
     * collection's iterator. At the end of the parallel iteration, the parallel
     * team threads wait for each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theIterable Iterable collection containing the items.
     * @param theIteration Parallel iteration.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theIterable</code> is null or
     * <code>theIteration</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theIteration</code>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(Iterable<T> theIterable,
            ParallelIteration<T> theIteration)
            throws Exception {
        execute(theIterable, theIteration, BarrierAction.WAIT);
    }

    /**
     * Execute a parallel iteration within this parallel region. For further
     * information, see class {@linkplain ParallelIteration}. The items
     * processed by the iteration are the items contained in the given iterable
     * collection. The iteration order is that of the given iterable
     * collection's iterator. At the end of the parallel iteration, the parallel
     * team threads encounter a barrier, and their behavior depends on the given
     * {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theIterable Iterable collection containing the items.
     * @param theIteration Parallel iteration.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theIterable</code> is null. Thrown if <code>theIteration</code> is null.
     * Thrown if <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of <code>theIteration</code>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(Iterable<T> theIterable,
            ParallelIteration<T> theIteration,
            BarrierAction action)
            throws Exception {
        // Verify preconditions.
        if (theIterable == null) {
            throw new NullPointerException("ParallelRegion.execute(): Iterable collection is null");
        }
        if (theIteration == null) {
            throw new NullPointerException("ParallelRegion.execute(): Parallel iteration is null");
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theIteration.myTeam = this.myTeam;

            // Get current parallel team thread.
            ParallelTeamThread currentThread = getCurrentThread();

            // Do top-of-parallel-construct processing.
            ItemGenerator<T> generator = null;
            if (currentThread.arriveAtParallelConstruct()) {
                // First thread to arrive sets up the shared item generator
                // object and stores it (or an exception if any) in each team
                // thread.
                try {
                    generator
                            = new IteratorItemGenerator<T>(theIterable.iterator());
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setItemGenerator(generator);
                    }
                } catch (Throwable exc) {
                    for (ParallelTeamThread thread : myTeam.myThread) {
                        thread.setConstructException(exc);
                    }
                }
            }

            // Get the shared item generator object.
            generator = castItemGenerator(currentThread.getItemGenerator());
            theIteration.myItemGenerator = generator;

            // Prepare to catch exceptions thrown by the parallel iteration
            // body.
            Throwable runException = null;
            try {
                // Perform per-thread initialization.
                theIteration.start();

                // Repeatedly get and process an item.
                ItemHolder<T> itemholder;
                while ((itemholder = generator.nextItem()) != null) {
                    theIteration.commonRun(itemholder.mySequenceNumber, itemholder.myItem);
                }

                // Perform per-thread finalization.
                theIteration.finish();
            } catch (Throwable exc) {
                runException = exc;
                generator.myBreak = true;
            }

            // Barrier synchronization.
            action.doBarrier(currentThread);

            // Propagate any exception thrown by the run() method.
            ParallelTeam.rethrow(runException);
        } finally {
            // Forget parallel team.
            theIteration.myTeam = null;
            theIteration.myItemGenerator = null;
        }
    }

    /**
     * Execute a parallel section within this parallel region. The parallel
     * section's <code>run()</code> method is called by one of the parallel team
     * threads. For further information, see class {@linkplain ParallelSection}.
     * At the end of the parallel section, the parallel team threads wait for
     * each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param section Parallel section.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>section</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if the parallel section's <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection section)
            throws Exception {
        execute(new ParallelSection[]{section}, BarrierAction.WAIT);
    }

    /**
     * Execute a parallel section within this parallel region. The parallel
     * section's <code>run()</code> method is called by one of the parallel team
     * threads. For further information, see class {@linkplain ParallelSection}.
     * At the end of the parallel section, the parallel team threads encounter a
     * barrier, and their behavior depends on the given {@linkplain
     * BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param section Parallel section.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>section</code> is null. Thrown if
     * <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if the parallel section's <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection section,
            BarrierAction action)
            throws Exception {
        execute(new ParallelSection[]{section}, action);
    }

    /**
     * Execute a group of two parallel sections concurrently within this
     * parallel region. Each parallel section's <code>run()</code> method is called
     * by a different parallel team thread. For further information, see class
     * {@linkplain ParallelSection}. At the end of the parallel section group,
     * the parallel team threads wait for each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param section1 First parallel section.
     * @param section2 Second parallel section.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>section1</code> is null. Thrown if
     * <code>section2</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of the parallel sections'
     * <code>run()</code> methods throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection section1,
            ParallelSection section2)
            throws Exception {
        execute(new ParallelSection[]{section1, section2},
                BarrierAction.WAIT);
    }

    /**
     * Execute a group of two parallel sections concurrently within this
     * parallel region. Each parallel section's <code>run()</code> method is called
     * by a different parallel team thread. For further information, see class
     * {@linkplain ParallelSection}. At the end of the parallel section group,
     * the parallel team threads encounter a barrier, and their behavior depends
     * on the given {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param section1 First parallel section.
     * @param section2 Second parallel section.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>section1</code> is null. Thrown if
     * <code>section2</code> is null. Thrown if <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of the parallel sections'
     * <code>run()</code> methods throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection section1,
            ParallelSection section2,
            BarrierAction action)
            throws Exception {
        execute(new ParallelSection[]{section1, section2},
                action);
    }

    /**
     * Execute a group of three parallel sections concurrently within this
     * parallel region. Each parallel section's <code>run()</code> method is called
     * by a different parallel team thread. For further information, see class
     * {@linkplain ParallelSection}. At the end of the parallel section group,
     * the parallel team threads wait for each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param section1 First parallel section.
     * @param section2 Second parallel section.
     * @param section3 Third parallel section.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>section1</code> is null. Thrown if
     * <code>section2</code> is null. Thrown if <code>section3</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of the parallel sections'
     * <code>run()</code> methods throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection section1,
            ParallelSection section2,
            ParallelSection section3)
            throws Exception {
        execute(new ParallelSection[]{section1, section2, section3},
                BarrierAction.WAIT);
    }

    /**
     * Execute a group of three parallel sections concurrently within this
     * parallel region. Each parallel section's <code>run()</code> method is called
     * by a different parallel team thread. For further information, see class
     * {@linkplain ParallelSection}. At the end of the parallel section group,
     * the parallel team threads encounter a barrier, and their behavior depends
     * on the given {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param section1 First parallel section.
     * @param section2 Second parallel section.
     * @param section3 Third parallel section.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>section1</code> is null. Thrown if
     * <code>section2</code> is null. Thrown if <code>section3</code> is null. Thrown if
     * <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of the parallel sections'
     * <code>run()</code> methods throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection section1,
            ParallelSection section2,
            ParallelSection section3,
            BarrierAction action)
            throws Exception {
        execute(new ParallelSection[]{section1, section2, section3},
                action);
    }

    /**
     * Execute a group of parallel sections concurrently within this parallel
     * region. Each parallel section's <code>run()</code> method is called by a
     * different parallel team thread. For further information, see class
     * {@linkplain ParallelSection}. At the end of the parallel section group,
     * the parallel team threads wait for each other at a barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param sections Parallel sections.
     * @exception NullPointerException (unchecked exception) Thrown if any of
     * the <code>sections</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of the parallel sections'
     * <code>run()</code> methods throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection[] sections)
            throws Exception {
        execute(sections, BarrierAction.WAIT);
    }

    /**
     * Execute a group of parallel sections concurrently within this parallel
     * region. Each parallel section's <code>run()</code> method is called by a
     * different parallel team thread. For further information, see class
     * {@linkplain ParallelSection}. At the end of the parallel section group,
     * the parallel team threads encounter a barrier, and their behavior depends
     * on the given {@linkplain BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>execute()</code> method with identical arguments, or none of the
     * threads must call the <code>execute()</code> method.
     *
     * @param sections Parallel sections.
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if any of
     * the <code>sections</code> is null. Thrown if <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if one of the parallel sections'
     * <code>run()</code> methods throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelSection[] sections,
            BarrierAction action)
            throws Exception {
        if (sections == null) {
            throw new NullPointerException("ParallelRegion.execute(): sections is null");
        }
        for (ParallelSection section : sections) {
            if (section == null) {
                throw new NullPointerException("ParallelRegion.execute(): A parallel section is null");
            }
        }
        if (action == null) {
            throw new NullPointerException("ParallelRegion.execute(): Barrier action is null");
        }

        execute(sections, new ParallelIteration<ParallelSection>() {
            public void run(ParallelSection section) throws Exception {
                try {
                    section.myTeam = this.myTeam;
                    section.run();
                } finally {
                    section.myTeam = null;
                }
            }
        },
                action);
    }

    /**
     * Perform a section of code in a critical region with exclusive locking.
     * The locking is performed using the default lock, a hidden {@linkplain
     * Lock} variable shared by all the parallel team threads. The thread
     * calling the <code>critical()</code> method waits until no other thread is
     * executing a critical region with exclusive locking using the default lock
     * and no other thread is executing a critical region with nonexclusive
     * locking using the default lock. The thread then calls
     * <code>theSection</code>'s <code>run()</code> method with exclusive locking using
     * the default lock. When the <code>run()</code> method returns, the thread
     * unlocks the lock and returns from the <code>critical()</code> method.
     * <P>
     * If the parallel section's <code>run()</code> method throws an exception, the
     * <code>critical()</code> method throws that same exception in the thread that
     * called the <code>run()</code> method (after unlocking the lock).
     *
     * @param theSection Parallel section to execute in the critical region.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theSection</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if <code>theSection</code>'s <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void critical(ParallelSection theSection)
            throws Exception {
        critical(myLock, theSection);
    }

    /**
     * Perform a section of code in a critical region with nonexclusive locking.
     * The locking is performed using the default lock, a hidden {@linkplain
     * Lock} variable shared by all the parallel team threads. The thread
     * calling the <code>critical()</code> method waits until no other thread is
     * executing a critical region with exclusive locking using the default
     * lock. However, any number of other threads may be executing a critical
     * region with nonexclusive locking using the default lock. The thread then
     * calls <code>theSection</code>'s <code>run()</code> method with nonexclusive
     * locking using the default lock. When the <code>run()</code> method returns,
     * the thread unlocks the lock and returns from the <code>critical()</code>
     * method.
     * <P>
     * If the parallel section's <code>run()</code> method throws an exception, the
     * <code>critical()</code> method throws that same exception in the thread that
     * called the <code>run()</code> method (after unlocking the lock).
     *
     * @param theSection Parallel section to execute in the critical region.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theSection</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if <code>theSection</code>'s <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void criticalNonexclusive(ParallelSection theSection)
            throws Exception {
        criticalNonexclusive(myLock, theSection);
    }

    /**
     * Perform a section of code in a critical region with exclusive locking
     * using the given lock. The thread calling the <code>critical()</code> method
     * waits until no other thread is executing a critical region with exclusive
     * locking using the given lock and no other thread is executing a critical
     * region with nonexclusive locking using the given lock. The thread then
     * calls <code>theSection</code>'s <code>run()</code> method with exclusive locking
     * using the given lock. When the <code>run()</code> method returns, the thread
     * unlocks the lock and returns from the <code>critical()</code> method.
     * <P>
     * If the parallel section's <code>run()</code> method throws an exception, the
     * <code>critical()</code> method throws that same exception in the thread that
     * called the <code>run()</code> method (after unlocking the lock).
     *
     * @param theLock Lock.
     * @param theSection Parallel section to execute in the critical region.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLock</code> is null or
     * <code>theSection</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if <code>theSection</code>'s <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void critical(Lock theLock,
            ParallelSection theSection)
            throws Exception {
        // Verify preconditions.
        if (theLock == null) {
            throw new NullPointerException("ParallelRegion.critical(): Lock is null");
        }
        if (theSection == null) {
            throw new NullPointerException("ParallelRegion.critical(): Parallel section is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.critical(): No parallel team executing");
        }

        // Lock the lock.
        theLock.lockExclusive();

        // Process the parallel section.
        try {
            theSection.myTeam = this.myTeam;
            theSection.run();
        } // Unlock the lock.
        finally {
            theSection.myTeam = null;
            theLock.unlockExclusive();
        }
    }

    /**
     * Perform a section of code in a critical region with nonexclusive locking
     * using the given lock. The thread calling the <code>critical()</code> method
     * waits until no other thread is executing a critical region with exclusive
     * locking using the given lock. However, any number of other threads may be
     * executing a critical region with nonexclusive locking using the given
     * lock. The thread then calls <code>theSection</code>'s <code>run()</code> method
     * with nonexclusive locking using the given lock. When the <code>run()</code>
     * method returns, the thread unlocks the lock and returns from the
     * <code>critical()</code> method.
     * <P>
     * If the parallel section's <code>run()</code> method throws an exception, the
     * <code>critical()</code> method throws that same exception in the thread that
     * called the <code>run()</code> method (after unlocking the lock).
     *
     * @param theLock Lock.
     * @param theSection Parallel section to execute in the critical region.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theLock</code> is null or
     * <code>theSection</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if <code>theSection</code>'s <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void criticalNonexclusive(Lock theLock,
            ParallelSection theSection)
            throws Exception {
        // Verify preconditions.
        if (theLock == null) {
            throw new NullPointerException("ParallelRegion.criticalNonexclusive(): Lock is null");
        }
        if (theSection == null) {
            throw new NullPointerException("ParallelRegion.criticalNonexclusive(): Parallel section is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("ParallelRegion.criticalNonexclusive(): No parallel team executing");
        }

        // Lock the lock.
        theLock.lockNonexclusive();

        // Process the parallel section.
        try {
            theSection.myTeam = this.myTeam;
            theSection.run();
        } // Unlock the lock.
        finally {
            theSection.myTeam = null;
            theLock.unlockNonexclusive();
        }
    }

    /**
     * Perform a barrier. The parallel team threads wait for each other at a
     * barrier.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>barrier()</code> method, or none of the threads must call the
     * <code>barrier()</code> method.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     */
    public final void barrier() {
        getCurrentThread().barrier();
    }

    /**
     * Perform a barrier, with a barrier action. The parallel team threads
     * encounter a barrier, and their behavior depends on the given {@linkplain
     * BarrierAction}.
     * <P>
     * <I>Note:</I> Either all threads in the parallel team must call the
     * <code>barrier()</code> method, or none of the threads must call the
     * <code>barrier()</code> method.
     *
     * @param action Barrier action.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>action</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel region.
     * @exception Exception Thrown if <code>theSection</code>'s <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final void barrier(BarrierAction action)
            throws Exception {
        if (action == null) {
            throw new NullPointerException("ParallelRegion.barrier(): Barrier action is null");
        }
        action.doBarrier(getCurrentThread());
    }

}
