//******************************************************************************
//
// File:    BarrierAction.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.BarrierAction
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
 * Class BarrierAction is the abstract base class for an object containing code
 * that is executed as part of a barrier wait. A barrier wait occurs in these
 * situations within a {@linkplain ParallelRegion}:
 * <UL>
 * <LI>
 * When the {@linkplain ParallelTeam} threads finish executing all the
 * iterations of a parallel for loop.
 * </LI>
 * <LI>
 * When the {@linkplain ParallelTeam} threads finish executing all the
 * {@linkplain ParallelSection}s in a group of parallel sections.
 * </LI>
 * <LI>
 * When the {@linkplain ParallelTeam} threads call the {@linkplain
 * ParallelRegion}'s <code>barrier()</code> method explicitly.
 * </LI>
 * </UL>
 * As each thread finishes executing one of the above constructs, each thread
 * encounters a barrier. What happens next depends on the barrier action
 * specified for that construct. There are three possibilities:
 * <UL>
 * <LI>
 * If the barrier action is omitted, or if the barrier action is {@link #WAIT
 * BarrierAction.WAIT}, each thread stops and waits at the barrier. When all
 * threads have arrived at the barrier, each thread resumes and proceeds to
 * execute whatever comes after the construct.
 * </LI>
 * <LI>
 * If the barrier action is {@link #NO_WAIT BarrierAction.NO_WAIT}, nothing
 * happens. The threads do not wait for each other. Each thread immediately
 * proceeds to execute whatever comes after the construct.
 * </LI>
 * <LI>
 * If the barrier action is an instance of class BarrierAction with the
 * <code>run()</code> method overridden, each thread stops and waits at the barrier.
 * When all threads have arrived at the barrier, <I>one</I> thread calls the
 * BarrierAction object's <code>run()</code> method. The particular thread that
 * calls the <code>run()</code> method is not specified. During this time the other
 * threads remain stopped. When the <code>run()</code> method returns, each thread
 * resumes and proceeds to execute whatever comes after the construct.
 * </LI>
 * </UL>
 * <p>
 * Thus, the barrier serves to synchronize all the threads at the end of a
 * parallel construct and possibly to execute a section of code in a single
 * thread.
 *
 * @author Alan Kaminsky
 * @version 11-Nov-2007
 */
public abstract class BarrierAction
        extends ParallelConstruct {

// Exported constructors.

    /**
     * Construct a new barrier action.
     */
    public BarrierAction() {
    }

// Exported operations.

    /**
     * Execute this barrier action. The <code>run()</code> method is called by a
     * single thread after all threads have arrived at the barrier.
     * <p>
     * The <code>run()</code> method must be implemented in a subclass.
     *
     * @throws Exception           The <code>run()</code> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public abstract void run()
            throws Exception;

// Hidden operations.

    /**
     * Execute a barrier.
     *
     * @param currentThread Parallel team thread calling <code>doBarrier()</code>.
     * @throws Exception The <code>run()</code> method may throw any exception.
     */
    void doBarrier(ParallelTeamThread currentThread)
            throws Exception {
        // Default is to do a barrier wait with this as the barrier action.
        currentThread.barrier(this);
    }

// Exported constants.
    /**
     * Do a barrier wait, without executing any code in a single thread.
     */
    public static final BarrierAction WAIT = new BarrierAction() {
        public void run() {
        }

        void doBarrier(ParallelTeamThread currentThread)
                throws Exception {
            currentThread.barrier();
        }
    };

    /**
     * Do not do a barrier wait.
     */
    public static final BarrierAction NO_WAIT = new BarrierAction() {
        public void run() {
        }

        void doBarrier(ParallelTeamThread currentThread)
                throws Exception {
        }
    };

}
