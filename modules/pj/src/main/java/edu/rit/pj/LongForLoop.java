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

/**
 * Class LongForLoop is the abstract base class for one variation of a parallel
 * for loop that is executed inside a {@linkplain ParallelRegion}. The loop
 * index data type is <code>long</code>. The loop stride is implicit (+1).
 * <P>
 * To execute a parallel for loop, create a {@linkplain ParallelRegion} object;
 * create an instance of a concrete subclass of class LongForLoop; and pass this
 * instance to the parallel region's <code>execute()</code> method. Either every
 * parallel team thread must call the parallel region's <code>execute()</code>
 * method with identical arguments, or every thread must not call the
 * <code>execute()</code> method. You can do all this using an anonymous inner
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
 * The parallel region's <code>execute()</code> method does the following. Each
 * parallel team thread calls the parallel for loop's <code>start()</code> method
 * once before beginning any loop iterations. The range of loop indexes is
 * divided into "chunks" and the chunks are apportioned among the threads, in a
 * manner determined by the parallel for loop's schedule as returned by the
 * <code>schedule()</code> method. Each thread repeatedly calls the parallel for
 * loop's <code>run()</code> method, passing in a different chunk on each call,
 * until all the chunks assigned to that thread have been performed. When a
 * thread has finished calling <code>run()</code>, the thread calls the parallel for
 * loop's <code>finish()</code> method. Then the thread waits at an implicit
 * barrier. When all the threads have reached the barrier, the
 * <code>execute()</code> method returns.
 * <P>
 * Note that each parallel team thread actually creates its own instance of the
 * parallel for loop class and passes that instance to the parallel region's
 * <code>execute()</code> method. Thus, any fields declared in the parallel for loop
 * class will <I>not</I> be shared by all the threads, but instead will be
 * private to each thread.
 * <P>
 * The <code>start()</code> method is intended for performing per-thread
 * initialization before starting the loop iterations. If no such initialization
 * is needed, omit the <code>start()</code> method.
 * <P>
 * The <code>run()</code> method contains the code for the loop. The first and last
 * indexes for a chunk of loop iterations are passed in as arguments. The loop
 * stride is implicit (+1). The parallel for loop's <code>run()</code> method must
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
 * </PRE> with the loop indexes running from <code>first</code> to <code>last</code>
 * inclusive and increasing by +1 on each iteration.
 * <P>
 * The <code>finish()</code> method is intended for performing per-thread
 * finalization after finishing the loop iterations. If no such finalization is
 * needed, omit the <code>finish()</code> method.
 * <P>
 * Sometimes a portion of a parallel for loop has to be executed sequentially in
 * the order of the loop indexes, while the rest of the parallel for loop can be
 * executed concurrently. For example, the loop body is performing some
 * computation that can be executed in parallel for different loop indexes, but
 * the results of each computation must be written to a file sequentially in the
 * order of the loop indexes. The <code>ordered()</code> method is provided for this
 * purpose. A call to the <code>ordered()</code> method may appear once in the
 * parallel for loop's <code>run()</code> method, like so:
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
 * </PRE> When called, the <code>ordered()</code> method waits until the
 * <code>ordered()</code>
 * method has been called and has returned in all loop iterations prior to the
 * current loop iteration. Then the <code>ordered()</code> method calls the given
 * parallel section's <code>run()</code> method. When the parallel section's
 * <code>run()</code> method returns, the <code>ordered()</code> method returns. If the
 * parallel section's <code>run()</code> method throws an exception, the
 * <code>ordered()</code> method throws that same exception.
 * <P>
 * It is possible to stop a parallel for loop using the <code>stopLoop()</code>
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
 * </PRE> Once <code>stopLoop()</code> is called, after each parallel team thread
 * finishes executing its current chunk of iterations, each thread will execute
 * no further chunks and will proceed to finish the parallel for loop. Note well
 * that stopping a parallel for loop is not the same as executing a
 * <code>break</code> statement in a regular for loop. The parallel for loop does
 * not stop until each thread, <I>including the thread that called
 * <code>stopLoop()</code></I>, has finished its current <I>chunk</I> of iterations.
 * Thus, depending on the parallel for loop's schedule, additional iterations
 * may be executed after <code>stopLoop()</code> is called. (The <code>break</code>
 * statement in the above example causes the thread that called
 * <code>stopLoop()</code> to finish its chunk of iterations early.)
 * <P>
 * Normally, at the end of the parallel for loop, the parallel team threads wait
 * for each other at a barrier. To eliminate this barrier wait, include
 * {@link edu.rit.pj.BarrierAction#NO_WAIT BarrierAction.NO_WAIT} in the <code>execute()</code>
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
 * the <code>execute()</code> method call. The barrier action object's
 * <code>run()</code> method contains the code to be executed in a single thread
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
 * If the parallel for loop's <code>start()</code>, <code>run()</code>, or
 * <code>finish()</code> method throws an exception in one of the threads, then that
 * thread executes no further code in the loop, and the parallel region's
 * <code>execute()</code> method throws that same exception in that thread.
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
     * The <code>schedule()</code> method may be overridden in a subclass to return
     * the desired schedule. If not overridden, the default is a runtime
     * schedule (see {@link edu.rit.pj.LongSchedule#runtime()}).
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
     * Execute one chunk of iterations of this parallel for loop. The
     * <code>run()</code> method must perform the loop body for indexes
     * <code>first</code> through <code>last</code> inclusive, increasing the loop index
     * by +1 after each iteration.
     * <P>
     * The <code>run()</code> method must be overridden in a subclass.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @exception Exception The <code>run()</code> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public abstract void run(long first,
            long last)
            throws Exception;

    /**
     * Perform per-thread finalization actions after finishing the loop
     * iterations.
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
     * Execute the given section of code in order of the loop indexes. A call to
     * the <code>ordered()</code> method may appear in this parallel for loop's
     * <code>run()</code> method. When called, the <code>ordered()</code> method waits
     * until the <code>ordered()</code> method has been called and has returned in
     * all loop iterations prior to the current loop iteration. Then the
     * <code>ordered()</code> method calls the <code>run()</code> method of
     * <code>theParallelSection</code>. When the parallel section's <code>run()</code>
     * method returns, the <code>ordered()</code> method returns. If the parallel
     * section's <code>run()</code> method throws an exception, the
     * <code>ordered()</code> method throws that same exception.
     * <P>
     * The <code>ordered()</code> method is used when a portion of a parallel for
     * loop has to be executed sequentially in the order of the loop indexes,
     * while the rest of the parallel for loop can be executed concurrently.
     * <P>
     * <I>Note:</I> Either the <code>ordered()</code> method must be called exactly
     * once during each call of the parallel for loop's <code>run()</code> method,
     * or the <code>ordered()</code> method must not be called at all.
     *
     * @param theSection Parallel section to execute in order.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theSection</code> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * parallel team is executing this parallel for loop.
     * @exception Exception Thrown if <code>theSection</code>'s <code>run()</code>
     * method throws an exception.
     * @throws java.lang.Exception if any.
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
     * Stop this parallel for loop. Once <code>stopLoop()</code> is called, after
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
     * performs common processing, then calls the <code>run()</code> method.
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

    // Kludge to avert false sharing in multithreaded programs.
    // Padding fields.
    volatile long p0 = 1000L;
    volatile long p1 = 1001L;
    volatile long p2 = 1002L;
    volatile long p3 = 1003L;
    volatile long p4 = 1004L;
    volatile long p5 = 1005L;
    volatile long p6 = 1006L;
    volatile long p7 = 1007L;
    volatile long p8 = 1008L;
    volatile long p9 = 1009L;
    volatile long pa = 1010L;
    volatile long pb = 1011L;
    volatile long pc = 1012L;
    volatile long pd = 1013L;
    volatile long pe = 1014L;
    volatile long pf = 1015L;

    // Method to prevent the JDK from optimizing away the padding fields.
    long preventOptimization() {
        return p0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 +
            p8 + p9 + pa + pb + pc + pd + pe + pf;
    }
    
}
