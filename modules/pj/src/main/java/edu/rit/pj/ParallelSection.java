//******************************************************************************
//
// File:    ParallelSection.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelSection
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
 * Class ParallelSection is the abstract base class for a section of code to be
 * executed in parallel.
 * <P>
 * A group of parallel sections may be executed concurrently by calling a
 * {@linkplain ParallelRegion}'s <code>execute()</code> method. Each section in the
 * group is executed by a different thread in the parallel thread team. Here is
 * one way to code a parallel section group with two parallel sections:
 * <PRE>
 *     new ParallelTeam(2).execute (new ParallelRegion()
 *         {
 *         public void run() throws Exception
 *             {
 *             execute
 *                 (new ParallelSection()
 *                     {
 *                     public void run()
 *                         {
 *                         // Code for first concurrent parallel section
 *                         . . .
 *                         }
 *                     },
 *                  new ParallelSection()
 *                     {
 *                     public void run()
 *                         {
 *                         // Code for second concurrent parallel section
 *                         . . .
 *                         }
 *                     });
 *             }
 *         });
 * </PRE> Here is another way to code a parallel section group with two parallel
 * sections:
 * <PRE>
 *     ParallelSection section_1 = new ParallelSection()
 *         {
 *         public void run()
 *             {
 *             // Code for first concurrent parallel section
 *             . . .
 *             }
 *         };
 *     ParallelSection section_2 = new ParallelSection()
 *         {
 *         public void run()
 *             {
 *             // Code for second concurrent parallel section
 *             . . .
 *             }
 *         };
 *     new ParallelTeam(2).execute (new ParallelRegion()
 *         {
 *         public void run() throws Exception
 *             {
 *             execute (section_1, section_2);
 *             }
 *         });
 * </PRE> A parallel section group may contain any number of parallel sections.
 * There are overloaded <code>execute()</code> methods that take one, two, or three
 * parallel section arguments. If there are four or more parallel sections, put
 * them in an array (type <code>ParallelSection[]</code>) and pass the array to the
 * <code>execute()</code> method.
 * <P>
 * Normally, at the end of the parallel section group, the parallel team threads
 * wait for each other at a barrier. To eliminate this barrier wait, include
 * {@link edu.rit.pj.BarrierAction#NO_WAIT BarrierAction.NO_WAIT} in the <code>execute()</code>
 * method call:
 * <PRE>
 *     new ParallelTeam(2).execute (new ParallelRegion()
 *         {
 *         public void run() throws Exception
 *             {
 *             execute (section_1, section_2, BarrierAction.NO_WAIT);
 *             }
 *         });
 * </PRE> To execute a section of code in a single thread as part of the barrier
 * synchronization, include an instance of class {@linkplain BarrierAction} in
 * the <code>execute()</code> method call. The barrier action object's
 * <code>run()</code> method contains the code to be executed in a single thread
 * while the other threads wait:
 * <PRE>
 *     new ParallelTeam(2).execute (new ParallelRegion()
 *         {
 *         public void run() throws Exception
 *             {
 *             execute (section_1, section_2, new BarrierAction()
 *                 {
 *                 public void run()
 *                     {
 *                     // Single-threaded code goes here
 *                     . . .
 *                     }
 *                 });
 *             }
 *         });
 * </PRE> For further information, see class {@linkplain BarrierAction}.
 * <P>
 * A parallel section may be executed by one thread in the parallel thread team
 * by executing a parallel section group consisting of the one parallel section,
 * as shown above:
 * <PRE>
 *     new ParallelTeam().execute (new ParallelRegion()
 *         {
 *         public void run() throws Exception
 *             {
 *             . . .
 *             execute (new ParallelSection()
 *                 {
 *                 public void run()
 *                     {
 *                     // Code to be executed by one parallel team thread
 *                     . . .
 *                     }
 *                 });
 *             . . .
 *             }
 *         });
 * </PRE>
 * <P>
 * A parallel section may be executed in a mutually exclusive fashion by calling
 * a {@linkplain ParallelRegion}'s <code>critical()</code> or
 * <code>criticalNonexclusive()</code> method. For example:
 * <PRE>
 *     new ParallelTeam().execute (new ParallelRegion()
 *         {
 *         public void run() throws Exception
 *             {
 *             . . .
 *             critical (new ParallelSection()
 *                 {
 *                 public void run()
 *                     {
 *                     // Mutually exclusive code
 *                     . . .
 *                     }
 *                 });
 *             . . .
 *             }
 *         });
 * </PRE>
 * <P>
 * By calling a {@linkplain ParallelForLoop}'s <code>ordered()</code> method, a
 * parallel section may be executed in the order of the loop indexes within a
 * parallel for loop body, while the rest of the parallel for loop body executes
 * concurrently. See classes {@linkplain IntegerForLoop}, {@linkplain
 * IntegerStrideForLoop}, {@linkplain LongForLoop}, and {@linkplain
 * LongStrideForLoop} for further information.
 * <P>
 * By calling a {@linkplain ParallelIteration}'s <code>ordered()</code> method, a
 * parallel section may be executed in the order of the items within a parallel
 * iteration body, while the rest of the parallel iteration body executes
 * concurrently. See class {@linkplain ParallelIteration} for further
 * information.
 *
 * @author Alan Kaminsky
 * @version 11-Nov-2007
 */
public abstract class ParallelSection
        extends ParallelConstruct {

// Exported constructors.
    /**
     * Construct a new parallel section.
     */
    public ParallelSection() {
        super();
    }

// Exported operations.
    /**
     * Execute this parallel section.
     * <P>
     * The <code>run()</code> method must be implemented in a subclass.
     *
     * @exception Exception The <code>run()</code> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public abstract void run()
            throws Exception;

}
