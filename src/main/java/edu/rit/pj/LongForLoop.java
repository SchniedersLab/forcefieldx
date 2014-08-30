//******************************************************************************
//
// File:    LongForLoop.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.LongForLoop
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

/**
 * Class LongForLoop is the abstract base class for one variation of a parallel
 * for loop that is executed inside a {@linkplain ParallelRegion}. The loop
 * index data type is <TT>long</TT>. The loop stride is implicit (+1).
 * <P>
 * To execute a parallel for loop, create a {@linkplain ParallelRegion} object;
 * create an instance of a concrete subclass of class LongForLoop; and pass this
 * instance to the parallel region's <TT>execute()</TT> method. Either every
 * parallel team thread must call the parallel region's <TT>execute()</TT>
 * method with identical arguments, or every thread must not call the
 * <TT>execute()</TT> method. You can do all this using an anonymous inner
 * class; for example:
 * <PRE>
 *     new ParallelRegion()
 *         {
 *         . . .
 *         public void run()
 *             {
 *             . . .
 *             execute (0L, 99L, new LongForLoop()
 *                 {
 *                 // Thread local variable declarations
 *                 . . .
 *                 public void start()
 *                     {
 *                     // Per-thread pre-loop initialization code
 *                     . . .
 *                     }
 *                 public void run (long first, long last)
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
 * The parallel region's <TT>execute()</TT> method does the following. Each
 * parallel team thread calls the parallel for loop's <TT>start()</TT> method
 * once before beginning any loop iterations. The range of loop indexes is
 * divided into "chunks" and the chunks are apportioned among the threads, in a
 * manner determined by the parallel for loop's schedule as returned by the
 * <TT>schedule()</TT> method. Each thread repeatedly calls the parallel for
 * loop's <TT>run()</TT> method, passing in a different chunk on each call,
 * until all the chunks assigned to that thread have been performed. When a
 * thread has finished calling <TT>run()</TT>, the thread calls the parallel for
 * loop's <TT>finish()</TT> method. Then the thread waits at an implicit
 * barrier. When all the threads have reached the barrier, the
 * <TT>execute()</TT> method returns.
 * <P>
 * Note that each parallel team thread actually creates its own instance of the
 * parallel for loop class and passes that instance to the parallel region's
 * <TT>execute()</TT> method. Thus, any fields declared in the parallel for loop
 * class will <I>not</I> be shared by all the threads, but instead will be
 * private to each thread.
 * <P>
 * The <TT>start()</TT> method is intended for performing per-thread
 * initialization before starting the loop iterations. If no such initialization
 * is needed, omit the <TT>start()</TT> method.
 * <P>
 * The <TT>run()</TT> method contains the code for the loop. The first and last
 * indexes for a chunk of loop iterations are passed in as arguments. The loop
 * stride is implicit (+1). The parallel for loop's <TT>run()</TT> method must
 * be coded this way:
 * <PRE>
 *     public void run (long first, long last)
 *         {
 *         for (long i = first; i &lt;= last; ++ i)
 *             {
 *             // Loop body code
 *             . . .
 *             }
 *         }
 * </PRE> with the loop indexes running from <TT>first</TT> to <TT>last</TT>
 * inclusive and increasing by +1 on each iteration.
 * <P>
 * The <TT>finish()</TT> method is intended for performing per-thread
 * finalization after finishing the loop iterations. If no such finalization is
 * needed, omit the <TT>finish()</TT> method.
 * <P>
 * Sometimes a portion of a parallel for loop has to be executed sequentially in
 * the order of the loop indexes, while the rest of the parallel for loop can be
 * executed concurrently. For example, the loop body is performing some
 * computation that can be executed in parallel for different loop indexes, but
 * the results of each computation must be written to a file sequentially in the
 * order of the loop indexes. The <TT>ordered()</TT> method is provided for this
 * purpose. A call to the <TT>ordered()</TT> method may appear once in the
 * parallel for loop's <TT>run()</TT> method, like so:
 * <PRE>
 *     public void run (long first, long last)
 *         {
 *         for (long i = first; i &lt;= last; ++ i)
 *             {
 *             // This portion executed concurrently
 *             . . .
 *             ordered (new ParallelSection()
 *                 {
 *                 public void run()
 *                     {
 *                     // This portion executed sequentially
 *                     // in the order of the loop indexes
 *                     . . .
 *                     }
 *                 });
 *             // This portion executed concurrently again
 *             . . .
 *             }
 *         }
 * </PRE> When called, the <TT>ordered()</TT> method waits until the
 * <TT>ordered()</TT>
 * method has been called and has returned in all loop iterations prior to the
 * current loop iteration. Then the <TT>ordered()</TT> method calls the given
 * parallel section's <TT>run()</TT> method. When the parallel section's
 * <TT>run()</TT> method returns, the <TT>ordered()</TT> method returns. If the
 * parallel section's <TT>run()</TT> method throws an exception, the
 * <TT>ordered()</TT> method throws that same exception.
 * <P>
 * It is possible to stop a parallel for loop using the <TT>stopLoop()</TT>
 * method, like this:
 * <PRE>
 *     public void run (long first, long last)
 *         {
 *         for (long i = first; i &lt;= last; ++ i)
 *             {
 *             // Loop body
 *             . . .
 *             if (/&#42;time to stop the loop&#42;/)
 *                 {
 *                 stopLoop();
 *                 break;
 *                 }
 *             // More loop body
 *             . . .
 *             }
 *         }
 * </PRE> Once <TT>stopLoop()</TT> is called, after each parallel team thread
 * finishes executing its current chunk of iterations, each thread will execute
 * no further chunks and will proceed to finish the parallel for loop. Note well
 * that stopping a parallel for loop is not the same as executing a
 * <TT>break</TT> statement in a regular for loop. The parallel for loop does
 * not stop until each thread, <I>including the thread that called
 * <TT>stopLoop()</TT></I>, has finished its current <I>chunk</I> of iterations.
 * Thus, depending on the parallel for loop's schedule, additional iterations
 * may be executed after <TT>stopLoop()</TT> is called. (The <TT>break</TT>
 * statement in the above example causes the thread that called
 * <TT>stopLoop()</TT> to finish its chunk of iterations early.)
 * <P>
 * Normally, at the end of the parallel for loop, the parallel team threads wait
 * for each other at a barrier. To eliminate this barrier wait, include
 * {@link BarrierAction#NO_WAIT BarrierAction.NO_WAIT} in the <TT>execute()</TT>
 * method call:
 * <PRE>
 *     new ParallelRegion()
 *         {
 *         . . .
 *         public void run()
 *             {
 *             . . .
 *             execute (0L, 99L, new LongForLoop()
 *                 {
 *                 . . .
 *                 },
 *             BarrierAction.NO_WAIT);
 *             . . .
 *             }
 *         }
 * </PRE> To execute a section of code in a single thread as part of the barrier
 * synchronization, include an instance of class {@linkplain BarrierAction} in
 * the <TT>execute()</TT> method call. The barrier action object's
 * <TT>run()</TT> method contains the code to be executed in a single thread
 * while the other threads wait:
 * <PRE>
 *     new ParallelRegion()
 *         {
 *         . . .
 *         public void run()
 *             {
 *             . . .
 *             execute (0L, 99L, new LongForLoop()
 *                 {
 *                 . . .
 *                 },
 *             new BarrierAction()
 *                 {
 *                 public void run()
 *                     {
 *                     // Single-threaded code goes here
 *                     . . .
 *                     }
 *                 });
 *             . . .
 *             }
 *         }
 * </PRE> For further information, see class {@linkplain BarrierAction}.
 * <P>
 * If the parallel for loop's <TT>start()</TT>, <TT>run()</TT>, or
 * <TT>finish()</TT> method throws an exception in one of the threads, then that
 * thread executes no further code in the loop, and the parallel region's
 * <TT>execute()</TT> method throws that same exception in that thread.
 * Furthermore, the other threads in the parallel team also execute no further
 * code in the loop after finishing their current chunks. Thus, if one thread
 * throws an exception, the whole parallel for loop exits with some (perhaps
 * none) of the iterations unperformed.
 *
 * @author Alan Kaminsky
 * @version 11-Nov-2007
 */
public abstract class LongForLoop
        extends ParallelForLoop {

// Hidden data members.
    // Parallel for loop schedule.
    LongSchedule mySchedule;

    // Loop index for ordered() construct.
    long myOrderedIndex;

// Exported constructors.
    /**
     * Construct a new parallel for loop.
     */
    public LongForLoop() {
        super();
    }

// Exported operations.
    /**
     * Determine this parallel for loop's schedule. The schedule determines how
     * the loop iterations are apportioned among the parallel team threads. For
     * further information, see class {@linkplain LongSchedule}.
     * <P>
     * The <TT>schedule()</TT> method may be overridden in a subclass to return
     * the desired schedule. If not overridden, the default is a runtime
     * schedule (see {@link LongSchedule#runtime()}).
     *
     * @return Schedule for this parallel for loop.
     */
    public LongSchedule schedule() {
        return LongSchedule.runtime();
    }

    /**
     * Perform per-thread initialization actions before starting the loop
     * iterations.
     * <P>
     * The <TT>start()</TT> method may be overridden in a subclass. If not
     * overridden, the <TT>start()</TT> method does nothing.
     *
     * @exception Exception The <TT>start()</TT> method may throw any exception.
     */
    public void start()
            throws Exception {
    }

    /**
     * Execute one chunk of iterations of this parallel for loop. The
     * <TT>run()</TT> method must perform the loop body for indexes
     * <TT>first</TT> through <TT>last</TT> inclusive, increasing the loop index
     * by +1 after each iteration.
     * <P>
     * The <TT>run()</TT> method must be overridden in a subclass.
     *
     * @param first First loop index.
     * @param last Last loop index.
     *
     * @exception Exception The <TT>run()</TT> method may throw any exception.
     */
    public abstract void run(long first,
            long last)
            throws Exception;

    /**
     * Perform per-thread finalization actions after finishing the loop
     * iterations.
     * <P>
     * The <TT>finish()</TT> method may be overridden in a subclass. If not
     * overridden, the <TT>finish()</TT> method does nothing.
     *
     * @exception Exception The <TT>finish()</TT> method may throw any
     * exception.
     */
    public void finish()
            throws Exception {
    }

    /**
     * Execute the given section of code in order of the loop indexes. A call to
     * the <TT>ordered()</TT> method may appear in this parallel for loop's
     * <TT>run()</TT> method. When called, the <TT>ordered()</TT> method waits
     * until the <TT>ordered()</TT> method has been called and has returned in
     * all loop iterations prior to the current loop iteration. Then the
     * <TT>ordered()</TT> method calls the <TT>run()</TT> method of
     * <TT>theParallelSection</TT>. When the parallel section's <TT>run()</TT>
     * method returns, the <TT>ordered()</TT> method returns. If the parallel
     * section's <TT>run()</TT> method throws an exception, the
     * <TT>ordered()</TT> method throws that same exception.
     * <P>
     * The <TT>ordered()</TT> method is used when a portion of a parallel for
     * loop has to be executed sequentially in the order of the loop indexes,
     * while the rest of the parallel for loop can be executed concurrently.
     * <P>
     * <I>Note:</I> Either the <TT>ordered()</TT> method must be called exactly
     * once during each call of the parallel for loop's <TT>run()</TT> method,
     * or the <TT>ordered()</TT> method must not be called at all.
     *
     * @param theSection Parallel section to execute in order.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theSection</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel for loop.
     * @exception Exception Thrown if <TT>theSection</TT>'s <TT>run()</TT>
     * method throws an exception.
     */
    public final void ordered(ParallelSection theSection)
            throws Exception {
        // Verify preconditions.
        if (theSection == null) {
            throw new IllegalStateException("LongForLoop.ordered(): Parallel section is null");
        }
        if (myTeam == null) {
            throw new IllegalStateException("LongForLoop.ordered(): No parallel team executing");
        }

        // Wait until the ordered() construct has finished for all previous
        // iterations.
        if (mySchedule.myOrderedIndex.get() != myOrderedIndex) {
            Spinner spinner = new Spinner();
            while (mySchedule.myOrderedIndex.get() != myOrderedIndex) {
                spinner.spin();
            }
        }

        // Execute parallel section. Propagate any exception.
        theSection.myTeam = this.myTeam;
        try {
            theSection.run();
        } finally {
            theSection.myTeam = null;

            // Notify that the ordered construct has finished for this
            // iteration.
            ++this.myOrderedIndex;
            mySchedule.myOrderedIndex.set(this.myOrderedIndex);
        }
    }

    /**
     * Stop this parallel for loop. Once <TT>stopLoop()</TT> is called, after
     * each parallel team thread finishes executing its current chunk of
     * iterations, each thread will execute no further chunks and will proceed
     * to finish this parallel for loop.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel for loop.
     */
    public final void stopLoop() {
        if (myTeam == null) {
            throw new IllegalStateException("ParallelForLoop.stopLoop(): No parallel team executing");
        }
        mySchedule.myBreak = true;
    }

// Hidden operations.
    /**
     * Execute one chunk of iterations of this parallel for loop. This method
     * performs common processing, then calls the <TT>run()</TT> method.
     *
     * @param first First loop index.
     * @param last Last loop index.
     *
     * @exception Exception This method may throw any exception.
     */
    void commonRun(long first,
            long last)
            throws Exception {
        myOrderedIndex = first;
        run(first, last);
    }

}
