//******************************************************************************
//
// File:    ParallelIteration.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelIteration
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
 * Class ParallelIteration is the abstract base class for a parallel iteration
 * that is executed inside a {@linkplain ParallelRegion}. The parallel iteration
 * lets you iterate over a group of items, with a separate parallel team thread
 * processing each item. The generic type parameter T specifies the items' data
 * type. The items can be the elements of an array, the items obtained from an
 * {@linkplain java.util.Iterator Iterator}, or the items contained in an
 * {@linkplain java.lang.Iterable Iterable} collection.
 * <P>
 * To execute a parallel iteration, create a {@linkplain ParallelRegion} object;
 * create an instance of a concrete subclass of class ParallelIteration; and
 * pass this instance to the parallel region's <TT>execute()</TT> method. Either
 * every parallel team thread must call the parallel region's <TT>execute()</TT>
 * method with identical arguments, or every thread must not call the
 * <TT>execute()</TT> method. You can do all this using an anonymous inner
 * class; for example:
 * <PRE>
 *     new ParallelRegion()
 *         {
 *         ArrayList&lt;String&gt; list = new ArrayList&lt;String&gt;();
 *         . . .
 *         public void run()
 *             {
 *             . . .
 *             execute (list, new ParallelIteration&lt;String&gt;()
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
 *                     // Loop body code
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
 * The parallel region's <TT>execute()</TT> method does the following. One of
 * the parallel team threads sets up the source of the items to be iterated over
 * -- either an array's elements, an iterator's items, or an iterable
 * collection's contents. (Note that only <I>one</I> thread does this setup; but
 * because all threads must call the parallel region's <TT>execute()</TT> method
 * with identical arguments, it doesn't matter which thread does the setup.)
 * Each parallel team thread calls the parallel iteration's <TT>start()</TT>
 * method once before beginning any loop iterations. Each thread repeatedly
 * calls the parallel iteration's <TT>run()</TT> method, passing in a different
 * item on each call, until all the items have been processed. When a thread has
 * finished calling <TT>run()</TT>, the thread calls the parallel iteration's
 * <TT>finish()</TT> method. Then the thread waits at a barrier. When all the
 * threads have reached the barrier, the <TT>execute()</TT> method returns.
 * <P>
 * Note that each parallel team thread actually creates its own instance of the
 * parallel iteration class and passes that instance to the parallel region's
 * <TT>execute()</TT> method. Thus, any fields declared in the parallel
 * iteration class will <I>not</I> be shared by all the threads, but instead
 * will be private to each thread.
 * <P>
 * The <TT>start()</TT> method is intended for performing per-thread
 * initialization before starting the loop iterations. If no such initialization
 * is needed, omit the <TT>start()</TT> method.
 * <P>
 * The <TT>run()</TT> method contains the code for the loop body. It does
 * whatever processing is needed on the one item passed in as an argument. Note
 * that, unlike a parallel for loop (class {@linkplain ParallelForLoop}), a
 * parallel iteration is not "chunked;" each parallel team thread always
 * processes just one item at a time.
 * <P>
 * The <TT>finish()</TT> method is intended for performing per-thread
 * finalization after finishing the loop iterations. If no such finalization is
 * needed, omit the <TT>finish()</TT> method.
 * <P>
 * Sometimes a portion of a parallel iteration has to be executed sequentially
 * in the same order as the items' iteration order, while the rest of the
 * parallel iteration can be executed concurrently. For example, the loop body
 * is performing some computation that can be executed in parallel for different
 * items, but the results of each computation must be written to a file
 * sequentially in the items' iteration order. The <TT>ordered()</TT> method is
 * provided for this purpose. A call to the <TT>ordered()</TT> method may appear
 * once in the parallel iteration's <TT>run()</TT> method, like so:
 * <PRE>
 *     public void run (String item)
 *         {
 *         // This portion executed concurrently
 *         . . .
 *         ordered (new ParallelSection()
 *             {
 *             public void run()
 *                 {
 *                 // This portion executed sequentially
 *                 // in the items' iteration order
 *                 . . .
 *                 }
 *             });
 *         // This portion executed concurrently again
 *         . . .
 *         }
 * </PRE>
 * When called, the <TT>ordered()</TT> method waits until the <TT>ordered()</TT>
 * method has been called and has returned for all items prior to the current
 * item. Then the <TT>ordered()</TT> method calls the given parallel section's
 * <TT>run()</TT> method. When the parallel section's <TT>run()</TT> method
 * returns, the <TT>ordered()</TT> method returns. If the parallel section's
 * <TT>run()</TT> method throws an exception, the <TT>ordered()</TT> method
 * throws that same exception.
 * <P>
 * It is possible to stop a parallel iteration using the <TT>stopLoop()</TT>
 * method, like this:
 * <PRE>
 *     public void run (String item)
 *         {
 *         // Loop body
 *         . . .
 *         if (/&#42;time to stop the loop&#42;/)
 *             {
 *             stopLoop();
 *             return;
 *             }
 *         // More loop body
 *         . . .
 *         }
 * </PRE>
 * Once <TT>stopLoop()</TT> is called, after each parallel team thread finishes
 * processing its current item, each thread will process no further items and
 * will proceed to finish the parallel iteration. Note well that stopping a
 * parallel iteration is not the same as executing a <TT>break</TT> statement in
 * a regular loop. The parallel iteration does not stop until each thread,
 * <I>including the thread that called <TT>stopLoop()</TT></I>, has finished
 * processing its current item. Thus, processing may continue for a while after
 * <TT>stopLoop()</TT> is called. (The <TT>return</TT> statement in the above
 * example causes the thread that called <TT>stopLoop()</TT> to stop its
 * processing early.)
 * <P>
 * Normally, at the end of the parallel iteration, the parallel team threads
 * wait for each other at a barrier. To eliminate this barrier wait, include
 * {@link BarrierAction#NO_WAIT BarrierAction.NO_WAIT} in the <TT>execute()</TT>
 * method call:
 * <PRE>
 *     new ParallelRegion()
 *         {
 *         . . .
 *         public void run()
 *             {
 *             . . .
 *             execute (list, new ParallelIteration&lt;String&gt;()
 *                 {
 *                 . . .
 *                 },
 *             BarrierAction.NO_WAIT);
 *             . . .
 *             }
 *         }
 * </PRE>
 * To execute a section of code in a single thread as part of the barrier
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
 *             execute (list, new ParallelIteration&lt;String&gt;()
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
 * </PRE>
 * For further information, see class {@linkplain BarrierAction}.
 * <P>
 * If the parallel iteration's <TT>start()</TT>, <TT>run()</TT>, or
 * <TT>finish()</TT> method throws an exception in one of the threads, then that
 * thread executes no further code in the loop, and the parallel region's
 * <TT>execute()</TT> method throws that same exception in that thread.
 * Furthermore, the other threads in the parallel team also process no further
 * items after finishing their current items. Thus, if one thread throws an
 * exception, the whole parallel iteration exits with some (perhaps none) of the
 * iterations unperformed.
 *
 * @param  <T>  Data type of the items iterated over.
 *
 * @author  Alan Kaminsky
 * @version 11-Nov-2007
 */
public abstract class ParallelIteration<T>
	extends ParallelConstruct
	{

// Hidden data members.

	// Item generator, used to obtain the items, to break this parallel
	// iteration, and to implement the ordered() construct.
	ItemGenerator<T> myItemGenerator;

	// Iteration index for the ordered() construct.
	int myOrderedIndex;

// Exported constructors.

	/**
	 * Construct a new parallel iteration.
	 */
	public ParallelIteration()
		{
		super();
		}

// Exported operations.

	/**
	 * Perform per-thread initialization actions before starting the loop
	 * iterations.
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
	 * Process one item in this parallel iteration. The <TT>run()</TT> method
	 * must perform the loop body for the given item.
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
	 * Perform per-thread finalization actions after finishing the loop
	 * iterations.
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
	 * Execute the given section of code in the items' iteration order. A call
	 * to the <TT>ordered()</TT> method may appear in this parallel iteration's
	 * <TT>run()</TT> method. When called, the <TT>ordered()</TT> method waits
	 * until the <TT>ordered()</TT> method has been called and has returned for
	 * all items prior to the current item. Then the <TT>ordered()</TT> method
	 * calls the <TT>run()</TT> method of <TT>theParallelSection</TT>. When the
	 * parallel section's <TT>run()</TT> method returns, the <TT>ordered()</TT>
	 * method returns. If the parallel section's <TT>run()</TT> method throws an
	 * exception, the <TT>ordered()</TT> method throws that same exception.
	 * <P>
	 * The <TT>ordered()</TT> method is used when a portion of a parallel
	 * iteration has to be executed sequentially in the items' iteration order,
	 * while the rest of the parallel iteration can be executed concurrently.
	 * <P>
	 * <I>Note:</I> Either the <TT>ordered()</TT> method must be called exactly
	 * once during each call of the parallel iteration's <TT>run()</TT> method,
	 * or the <TT>ordered()</TT> method must not be called at all.
	 *
	 * @param  theSection  Parallel section to execute in order.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theSection</TT> is null.
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if no parallel team is executing this
	 *     parallel iteration.
	 * @exception  Exception
	 *     Thrown if <TT>theSection</TT>'s <TT>run()</TT> method throws an
	 *     exception.
	 */
	public final void ordered
		(ParallelSection theSection)
		throws Exception
		{
		// Verify preconditions.
		if (theSection == null)
			{
			throw new IllegalStateException
				("ParallelIteration.ordered(): Parallel section is null");
			}
		if (myTeam == null)
			{
			throw new IllegalStateException
				("ParallelIteration.ordered(): No parallel team executing");
			}

		// Wait until the ordered() construct has finished for all previous
		// iterations.
		if (myItemGenerator.myOrderedIndex != this.myOrderedIndex)
			{
			Spinner spinner = new Spinner();
			while (myItemGenerator.myOrderedIndex != this.myOrderedIndex)
				{
				spinner.spin();
				}
			}

		// Execute parallel section. Propagate any exception.
		theSection.myTeam = this.myTeam;
		try
			{
			theSection.run();
			}
		finally
			{
			theSection.myTeam = null;

			// Notify that the ordered construct has finished for this
			// iteration.
			++ this.myOrderedIndex;
			myItemGenerator.myOrderedIndex = this.myOrderedIndex;
			}
		}

	/**
	 * Stop this parallel iteration. Once <TT>stopLoop()</TT> is called, after
	 * each parallel team thread finishes processing its current item, each
	 * thread will process no further items and will proceed to finish this
	 * parallel iteration.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if no parallel team is executing this
	 *     parallel iteration.
	 */
	public final void stopLoop()
		{
		if (myTeam == null)
			{
			throw new IllegalStateException
				("ParallelIteration.stopLoop(): No parallel team executing");
			}
		myItemGenerator.myBreak = true;
		}

// Hidden operations.

	/**
	 * Execute one chunk of iterations of this parallel for loop. This method
	 * performs common processing, then calls the <TT>run()</TT> method.
	 *
	 * @param  index  Iteration index.
	 * @param  item   Item.
	 *
	 * @exception  Exception
	 *     This method may throw any exception.
	 */
	void commonRun
		(int index,
		 T item)
		throws Exception
		{
		myOrderedIndex = index;
		run (item);
		}

	}
