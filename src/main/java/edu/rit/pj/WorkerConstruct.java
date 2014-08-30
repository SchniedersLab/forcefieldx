//******************************************************************************
//
// File:    WorkerConstruct.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.WorkerConstruct
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

import edu.rit.pj.reduction.*;

/**
 * Class WorkerConstruct is the common base class for all worker constructs that
 * are executed by a {@linkplain WorkerTeam}.
 *
 * @author Alan Kaminsky
 * @version 19-Jan-2010
 */
public abstract class WorkerConstruct {

// Hidden data members.
    // 128 bytes of extra padding to avert cache interference.
    private long p0, p1, p2, p3, p4, p5, p6, p7;
    private long p8, p9, pa, pb, pc, pd, pe, pf;

    // Worker team that is executing this worker construct, or null if none.
    WorkerTeam myTeam;

// Exported constructors.
    /**
     * Construct a new worker construct.
     */
    public WorkerConstruct() {
    }

// Exported operations.
    /**
     * Determine if a worker team is executing this worker construct.
     *
     * @return True if a worker team is executing this worker construct, false
     * otherwise.
     */
    public final boolean isExecutingInParallel() {
        return myTeam != null;
    }

    /**
     * Returns the worker team that is executing this worker construct.
     *
     * @return Worker team.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker construct.
     */
    public final WorkerTeam team() {
        if (myTeam == null) {
            throw new IllegalStateException("WorkerConstruct.team(): No worker team executing");
        }
        return myTeam;
    }

    /**
     * Returns the worker region of code within which a worker team is executing
     * this worker construct.
     *
     * @return Worker region.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker construct.
     */
    public final WorkerRegion region() {
        if (myTeam == null) {
            throw new IllegalStateException("WorkerConstruct.region(): No worker team executing");
        }
        return myTeam.myRegion;
    }

    /**
     * Determine the number of worker threads in the current process in the
     * worker team executing this worker construct. This does not include the
     * master thread if any.
     *
     * @return Number of worker threads in the current process.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker construct.
     */
    public final int getThreadCount() {
        if (myTeam == null) {
            throw new IllegalStateException("WorkerConstruct.getThreadCount(): No worker team executing");
        }
        return myTeam.K;
    }

    /**
     * Determine the total number of worker threads in all processes in the
     * worker team executing this worker construct. This does not include the
     * master thread.
     *
     * @return Number of worker threads in all processes.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker construct.
     */
    public final int getTotalThreadCount() {
        if (myTeam == null) {
            throw new IllegalStateException("WorkerConstruct.getTotalThreadCount(): No worker team executing");
        }
        return myTeam.count;
    }

    /**
     * Determine the index of the calling thread in the worker team executing
     * this worker construct. Every worker thread in every process has a unique
     * index, going from index 0 for the first thread in the first process to
     * index <I>K</I>&minus;1 for the last thread in the last process, where
     * <I>K</I> is the total number of worker threads in all the processes. The
     * master thread's index is &minus;1.
     *
     * @return Index of the calling thread.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker construct. Thrown if the thread
     * calling
     * <TT>getThreadIndex()</TT> is not part of the worker team executing this
     * worker construct.
     */
    public final int getThreadIndex() {
        return getCurrentThread().myIndex;
    }

    /**
     * Determine if the calling thread is the master thread in the worker team
     * executing this worker construct.
     *
     * @return True if the calling thread is the master thread, false if the
     * calling thread is a worker thread.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if no
     * worker team is executing this worker construct. Thrown if the thread
     * calling
     * <TT>getThreadIndex()</TT> is not part of the worker team executing this
     * worker construct.
     */
    public final boolean isMasterThread() {
        return getThreadIndex() == -1;
    }

// Hidden operations.
    /**
     * Get the worker team thread that is calling this method.
     *
     * @return Worker team thread.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if the
     * calling thread is not one of the worker team threads executing this
     * worker construct.
     */
    WorkerTeamThread getCurrentThread() {
        if (myTeam == null) {
            throw new IllegalStateException("WorkerConstruct.getCurrentThread(): No worker team executing");
        }
        try {
            WorkerTeamThread current = (WorkerTeamThread) Thread.currentThread();
            if (current.myTeam != this.myTeam) {
                throw new IllegalStateException("WorkerConstruct.getCurrentThread(): Current thread is not executing this worker construct");
            }
            return current;
        } catch (ClassCastException exc) {
            throw new IllegalStateException("WorkerConstruct.getCurrentThread(): Current thread is not a worker team thread",
                    exc);
        }
    }

}
