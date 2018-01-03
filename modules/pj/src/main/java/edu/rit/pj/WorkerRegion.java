//******************************************************************************
//
// File:    WorkerRegion.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.WorkerRegion
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

import java.util.Iterator;

import edu.rit.util.LongRange;
import edu.rit.util.Range;

/**
 * Class WorkerRegion is the abstract base class for a worker region that is
 * executed by a {@linkplain WorkerTeam} of threads distributed across the
 * processes of a cluster parallel program.
 * <P>
 * To execute a worker region, create a {@linkplain WorkerTeam} object; create
 * an instance of a concrete subclass of class WorkerRegion; and pass this
 * instance to the worker team's <TT>execute()</TT> method. You can do all this
 * using an anonymous inner class; for example:
 * <PRE>
 *     new WorkerTeam().execute (new WorkerRegion()
 *         {
 *         public void start()
 *             {
 *             // Initialization code
 *             . . .
 *             }
 *         public void run()
 *             {
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
 * The worker team's <TT>execute()</TT> method does the following. In each
 * process, the worker team has a certain number of <B>worker threads</B>
 * <I>K</I>, where <I>K</I> was specified when the worker team was constructed.
 * In the highest-ranked process of the communicator, there is a <B>master
 * thread</B> in addition to the worker threads. The <B>main thread</B> is the
 * thread calling the worker team's <TT>execute()</TT> method. The main thread
 * calls the worker region's <TT>start()</TT> method. When the <TT>start()</TT>
 * method returns, all the worker threads, plus the master thread if any, call
 * the worker region's <TT>run()</TT> method concurrently. When all the team
 * threads have returned from the <TT>run()</TT> method, the main thread calls
 * the worker region's <TT>finish()</TT> method. When the <TT>finish()</TT>
 * method returns, the main thread returns from the worker team's
 * <TT>execute()</TT> method.
 * <P>
 * The chief purpose of a worker team is to execute a work-sharing parallel loop
 * in a cluster parallel program, partitioning the loop iterations among the
 * worker threads in all the processes. The worker team uses the
 * <B>master-worker pattern</B> to partition the iterations. The master thread
 * partitions the loop iterations and sends tasks to the worker threads; the
 * worker threads send results back to the master thread. The worker team uses a
 * certain <B>communicator</B> to do this message passing. The communicator was
 * specified when the worker team was constructed. For further information, see
 * class {@linkplain WorkerIntegerForLoop} and {@linkplain WorkerLongForLoop}.
 * <P>
 * Within each process, variables to be shared by all threads in the team may be
 * declared as fields of the WorkerRegion subclass. (Variables cannot be shared
 * between processes.) The <TT>start()</TT> method is intended for performing
 * initialization in a single thread before parallel execution begins. If no
 * such initialization is needed, omit the <TT>start()</TT> method. The
 * <TT>run()</TT> method contains code to be executed in parallel by all threads
 * in the team. Variables that are private to each thread may be declared inside
 * the <TT>run()</TT> method. The <TT>finish()</TT> method is intended for
 * performing finalization in a single thread after parallel execution ends. If
 * no such finalization is needed, omit the <TT>finish()</TT> method.
 * <P>
 * If the worker region's <TT>start()</TT> method throws an exception, the
 * worker team's <TT>execute()</TT> method throws that same exception, and the
 * <TT>run()</TT> method is not called.
 * <P>
 * If the worker region's <TT>run()</TT> method throws an exception in one of
 * the team threads, the exception's stack trace is printed on the standard
 * error, the worker team waits until all the other team threads have returned
 * from the <TT>run()</TT> method, then the worker team's <TT>execute()</TT>
 * method throws that same exception, and the worker region's <TT>finish()</TT>
 * method is not called. If the worker region's <TT>run()</TT> method throws an
 * exception in more than one of the team threads, each exception's stack trace
 * is printed on the standard error, the worker team waits until all the other
 * team threads have returned from the <TT>run()</TT> method, then the worker
 * team's <TT>execute()</TT> method throws a {@linkplain
 * MultipleParallelException} wrapping all the thrown exceptions, and the worker
 * region's <TT>finish()</TT> method is not called.
 * <P>
 * If the worker region's <TT>finish()</TT> method throws an exception, the
 * worker team's <TT>execute()</TT> method throws that same exception.
 *
 * @author Alan Kaminsky
 * @version 07-Oct-2010
 */
public abstract class WorkerRegion
        extends WorkerConstruct {

// Exported constructors.
    /**
     * Construct a new worker region.
     */
    public WorkerRegion() {
        super();
    }

// Exported operations.
    /**
     * Perform initialization actions before parallel execution begins. Only one
     * thread in each process calls the <TT>start()</TT> method.
     * <P>
     * The <TT>start()</TT> method may be overridden in a subclass. If not
     * overridden, the <TT>start()</TT> method does nothing.
     *
     * @exception Exception The <TT>start()</TT> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public void start()
            throws Exception {
    }

    /**
     * Execute parallel code. All threads of the worker team in each process
     * call the <TT>run()</TT> method concurrently.
     * <P>
     * The <TT>run()</TT> method must be implemented in a subclass.
     *
     * @exception Exception The <TT>run()</TT> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public abstract void run()
            throws Exception;

    /**
     * Perform finalization actions after parallel execution ends. Only one
     * thread in each process calls the <TT>finish()</TT> method.
     * <P>
     * The <TT>finish()</TT> method may be overridden in a subclass. If not
     * overridden, the <TT>finish()</TT> method does nothing.
     *
     * @exception Exception The <TT>finish()</TT> method may throw any
     * exception.
     * @throws java.lang.Exception if any.
     */
    public void finish()
            throws Exception {
    }

    /**
     * Execute a worker for loop within this worker region. For further
     * information, see class {@linkplain WorkerIntegerForLoop}. The loop index
     * goes from <TT>first</TT> (inclusive) to <TT>last</TT> (inclusive) in
     * steps of +1. If <TT>first</TT> is greater than <TT>last</TT>, then no
     * loop iterations are performed.
     * <P>
     * <I>Note:</I> Either all threads in the worker team must call the
     * <TT>execute()</TT> method with identical arguments, or none of the
     * threads must call the <TT>execute()</TT> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param theLoop Worker for loop.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theLoop</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker region.
     * @exception Exception Thrown if one of <TT>theLoop</TT>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(int first,
            int last,
            WorkerIntegerForLoop theLoop)
            throws Exception {
        // Verify preconditions.
        if (theLoop == null) {
            throw new NullPointerException("WorkerRegion.execute(): Worker for loop is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("WorkerRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            WorkerTeamThread currentThread = getCurrentThread();
            int w = currentThread.myIndex;

            // Do master or worker thread processing.
            Range range = new Range(first, last);
            if (w == -1) {
                theLoop.masterExecute(range);
            } else {
                theLoop.workerExecute(w, range);
            }
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
        }
    }

    /**
     * Execute a worker for loop within this worker region. For further
     * information, see class {@linkplain WorkerIntegerStrideForLoop}. The loop
     * index goes from <TT>first</TT> (inclusive) to <TT>last</TT> (inclusive)
     * in steps of <TT>stride</TT>. The stride must be positive. If
     * <TT>first</TT> is greater than <TT>last</TT>, then no loop iterations are
     * performed.
     * <P>
     * <I>Note:</I> Either all threads in the worker team must call the
     * <TT>execute()</TT> method with identical arguments, or none of the
     * threads must call the <TT>execute()</TT> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param stride Loop index stride, &gt;= 1.
     * @param theLoop Worker for loop.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>stride</TT> &lt; 1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theLoop</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker region.
     * @exception Exception Thrown if one of <TT>theLoop</TT>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(int first,
            int last,
            int stride,
            WorkerIntegerStrideForLoop theLoop)
            throws Exception {
        // Verify preconditions.
        if (stride <= 0) {
            throw new IllegalArgumentException("WorkerRegion.execute(): Stride = " + stride + " illegal");
        }
        if (theLoop == null) {
            throw new NullPointerException("WorkerRegion.execute(): Worker for loop is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("WorkerRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            WorkerTeamThread currentThread = getCurrentThread();
            int w = currentThread.myIndex;

            // Do master or worker thread processing.
            Range range = new Range(first, last, stride);
            if (w == -1) {
                theLoop.masterExecute(range);
            } else {
                theLoop.workerExecute(w, range);
            }
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
        }
    }

    /**
     * Execute a worker for loop within this worker region. For further
     * information, see class {@linkplain WorkerLongForLoop}. The loop index
     * goes from <TT>first</TT> (inclusive) to <TT>last</TT> (inclusive) in
     * steps of +1. If <TT>first</TT> is greater than <TT>last</TT>, then no
     * loop iterations are performed.
     * <P>
     * <I>Note:</I> Either all threads in the worker team must call the
     * <TT>execute()</TT> method with identical arguments, or none of the
     * threads must call the <TT>execute()</TT> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param theLoop Worker for loop.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theLoop</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker region.
     * @exception Exception Thrown if one of <TT>theLoop</TT>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(long first,
            long last,
            WorkerLongForLoop theLoop)
            throws Exception {
        // Verify preconditions.
        if (theLoop == null) {
            throw new NullPointerException("WorkerRegion.execute(): Worker for loop is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("WorkerRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            WorkerTeamThread currentThread = getCurrentThread();
            int w = currentThread.myIndex;

            // Do master or worker thread processing.
            LongRange range = new LongRange(first, last);
            if (w == -1) {
                theLoop.masterExecute(range);
            } else {
                theLoop.workerExecute(w, range);
            }
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
        }
    }

    /**
     * Execute a worker for loop within this worker region. For further
     * information, see class {@linkplain WorkerLongStrideForLoop}. The loop
     * index goes from <TT>first</TT> (inclusive) to <TT>last</TT> (inclusive)
     * in steps of <TT>stride</TT>. The stride must be positive. If
     * <TT>first</TT> is greater than <TT>last</TT>, then no loop iterations are
     * performed.
     * <P>
     * <I>Note:</I> Either all threads in the worker team must call the
     * <TT>execute()</TT> method with identical arguments, or none of the
     * threads must call the <TT>execute()</TT> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param stride Loop index stride, &gt;= 1.
     * @param theLoop Worker for loop.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>stride</TT> &lt; 1.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theLoop</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker region.
     * @exception Exception Thrown if one of <TT>theLoop</TT>'s methods throws
     * an exception.
     * @throws java.lang.Exception if any.
     */
    public final void execute(long first,
            long last,
            long stride,
            WorkerLongStrideForLoop theLoop)
            throws Exception {
        // Verify preconditions.
        if (stride <= 0) {
            throw new IllegalArgumentException("WorkerRegion.execute(): Stride = " + stride + " illegal");
        }
        if (theLoop == null) {
            throw new NullPointerException("WorkerRegion.execute(): Worker for loop is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("WorkerRegion.execute(): No parallel team executing");
        }

        try {
            // Record parallel team.
            theLoop.myTeam = this.myTeam;

            // Get current parallel team thread.
            WorkerTeamThread currentThread = getCurrentThread();
            int w = currentThread.myIndex;

            // Do master or worker thread processing.
            LongRange range = new LongRange(first, last, stride);
            if (w == -1) {
                theLoop.masterExecute(range);
            } else {
                theLoop.workerExecute(w, range);
            }
        } finally {
            // Forget parallel team.
            theLoop.myTeam = null;
        }
    }

    /**
     * Execute a worker iteration within this worker region. For further
     * information, see class {@linkplain WorkerIteration}. The items processed
     * by the iteration are the elements of the given array. The iteration order
     * is from index 0 upwards.
     * <P>
     * <I>Note:</I> Either all threads in the worker team must call the
     * <TT>execute()</TT> method with identical arguments, or none of the
     * threads must call the <TT>execute()</TT> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theArray Array containing the items.
     * @param theIteration Worker iteration.
     * @exception NullPointerException (unchecked exception) Thrown if this is
     * the master process and
     * <TT>theArray</TT> is null. Thrown if <TT>theIteration</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker region.
     * @exception Exception Thrown if one of <TT>theIteration</TT>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(T[] theArray,
            WorkerIteration<T> theIteration)
            throws Exception {
        // Verify preconditions.
        if (myTeam == null) {
            throw new IllegalStateException("WorkerRegion.execute(): No parallel team executing");
        }
        if (myTeam.rank == myTeam.masterRank() && theArray == null) {
            throw new NullPointerException("WorkerRegion.execute(): Array is null");
        }
        if (theIteration == null) {
            throw new NullPointerException("WorkerRegion.execute(): Worker iteration is null");
        }

        try {
            // Record parallel team.
            theIteration.myTeam = this.myTeam;

            // Get current parallel team thread.
            WorkerTeamThread currentThread = getCurrentThread();
            int w = currentThread.myIndex;

            // Do master or worker thread processing.
            if (w == -1) {
                theIteration.masterExecute(new ArrayItemGenerator<T>(theArray));
            } else {
                theIteration.workerExecute(w);
            }
        } finally {
            // Forget parallel team.
            theIteration.myTeam = null;
        }
    }

    /**
     * Execute a worker iteration within this worker region. For further
     * information, see class {@linkplain WorkerIteration}. The items processed
     * by the iteration are the items returned by the given iterator. The
     * iteration order is that of the given iterator.
     * <P>
     * <I>Note:</I> Either all threads in the worker team must call the
     * <TT>execute()</TT> method with identical arguments, or none of the
     * threads must call the <TT>execute()</TT> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theIterator Iterator over the items.
     * @param theIteration Worker iteration.
     * @exception NullPointerException (unchecked exception) Thrown if this is
     * the master process and
     * <TT>theIterator</TT> is null. Thrown if <TT>theIteration</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker region.
     * @exception Exception Thrown if one of <TT>theIteration</TT>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(Iterator<T> theIterator,
            WorkerIteration<T> theIteration)
            throws Exception {
        // Verify preconditions.
        if (myTeam == null) {
            throw new IllegalStateException("WorkerRegion.execute(): No parallel team executing");
        }
        if (myTeam.rank == myTeam.masterRank() && theIterator == null) {
            throw new NullPointerException("WorkerRegion.execute(): Iterator is null");
        }
        if (theIteration == null) {
            throw new NullPointerException("WorkerRegion.execute(): Worker iteration is null");
        }

        try {
            // Record parallel team.
            theIteration.myTeam = this.myTeam;

            // Get current parallel team thread.
            WorkerTeamThread currentThread = getCurrentThread();
            int w = currentThread.myIndex;

            // Do master or worker thread processing.
            if (w == -1) {
                theIteration.masterExecute(new IteratorItemGenerator<T>(theIterator));
            } else {
                theIteration.workerExecute(w);
            }
        } finally {
            // Forget parallel team.
            theIteration.myTeam = null;
        }
    }

    /**
     * Execute a worker iteration within this worker region. For further
     * information, see class {@linkplain WorkerIteration}. The items processed
     * by the iteration are the items contained in the given iterable
     * collection. The iteration order is that of the given iterable
     * collection's iterator.
     * <P>
     * <I>Note:</I> Either all threads in the worker team must call the
     * <TT>execute()</TT> method with identical arguments, or none of the
     * threads must call the <TT>execute()</TT> method.
     *
     * @param <T> Data type of the items iterated over.
     * @param theIterable Iterable collection containing the items.
     * @param theIteration Worker iteration.
     * @exception NullPointerException (unchecked exception) Thrown if this is
     * the master process and
     * <TT>theIterable</TT> is null. Thrown if <TT>theIteration</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker region.
     * @exception Exception Thrown if one of <TT>theIteration</TT>'s methods
     * throws an exception.
     * @throws java.lang.Exception if any.
     */
    public final <T> void execute(Iterable<T> theIterable,
            WorkerIteration<T> theIteration)
            throws Exception {
        // Verify preconditions.
        if (myTeam == null) {
            throw new IllegalStateException("WorkerRegion.execute(): No parallel team executing");
        }
        if (myTeam.rank == myTeam.masterRank() && theIterable == null) {
            throw new NullPointerException("WorkerRegion.execute(): Iterable collection is null");
        }
        if (theIteration == null) {
            throw new NullPointerException("WorkerRegion.execute(): Worker iteration is null");
        }

        try {
            // Record parallel team.
            theIteration.myTeam = this.myTeam;

            // Get current parallel team thread.
            WorkerTeamThread currentThread = getCurrentThread();
            int w = currentThread.myIndex;

            // Do master or worker thread processing.
            if (w == -1) {
                theIteration.masterExecute(new IteratorItemGenerator<T>(theIterable.iterator()));
            } else {
                theIteration.workerExecute(w);
            }
        } finally {
            // Forget parallel team.
            theIteration.myTeam = null;
        }
    }

}
