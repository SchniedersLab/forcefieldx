//******************************************************************************
//
// File:    WorkerTeam.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.WorkerTeam
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

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Semaphore;

/**
 * Class WorkerTeam provides a team of threads, distributed across the processes
 * of a cluster parallel program, for executing a {@linkplain WorkerRegion} in
 * parallel.
 * <P>
 * A worker team uses a communicator for message passing. The communicator is
 * specified as a constructor argument; if not specified, the world communicator
 * is used. Every process that is part of the communicator must create the
 * worker team. In class WorkerTeam, there is one worker thread per process. (To
 * get more than one worker thread per process, use class {@linkplain
 * HybridTeam}.) Every worker thread in every process has a unique index, going
 * from index 0 for the first thread in the first process to index
 * <I>K</I>&minus;1 for the last thread in the last process, where <I>K</I> is
 * the total number of worker threads in all the processes. In process rank 0,
 * there is an additional master thread.
 * <P>
 * To execute a worker region, create a WorkerTeam object; create an instance of
 * a concrete subclass of class {@linkplain WorkerRegion}; and pass this
 * instance to the worker team's <TT>execute()</TT> method. For further
 * information, see class {@linkplain WorkerRegion}.
 *
 * @author Alan Kaminsky
 * @version 19-Jan-2010
 */
public class WorkerTeam {

// Hidden data members.
    // Number of worker threads in this process.
    int K;

    // Communicator for message passing, its size, this process's rank.
    Comm comm;
    int size;
    int rank;

    // Number of worker threads in all processes.
    int count;

    // Array of worker and master team threads in this process. There are K
    // worker threads. There is an additional master thread in the last process
    // of the communicator.
    WorkerTeamThread[] myThread;

    // Worker region being executed, or null if none is being executed.
    WorkerRegion myRegion;

    // Semaphore for synchronizing threads at the end of a worker region.
    Semaphore myRegionEndSemaphore = new Semaphore(0);

    // Exception map for worker region, or null if none is being executed.
    ConcurrentHashMap<Integer, Throwable> myExceptionMap;

// Hidden constructors.
    /**
     * Construct a new, uninitialized worker team.
     *
     * @param flag To distinguish this constructor from the others.
     */
    WorkerTeam(boolean flag) {
    }

// Exported constructors.
    /**
     * Construct a new worker team with one thread per process and using the
     * world communicator for message passing.
     */
    public WorkerTeam() {
        this(Comm.world());
    }

    /**
     * Construct a new worker team with one thread per process and using the
     * given communicator for message passing.
     *
     * @param comm Communicator to use for message passing.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>comm</TT> is null.
     */
    public WorkerTeam(Comm comm) {
        if (comm == null) {
            throw new NullPointerException("WorkerTeam(): comm is null");
        }
        initialize(/*K    */1,
                /*comm */ comm,
                /*size */ comm.size(),
                /*rank */ comm.rank(),
                /*count*/ comm.size(),
                /*wlb  */ comm.rank());
    }

// Hidden initializers.
    /**
     * Initialize a new worker team.
     *
     * @param K Number of worker threads in this process.
     * @param comm Communicator to use for message passing.
     * @param size Communicator's size.
     * @param rank This process's rank in the communicator.
     * @param count Number of worker threads in all processes.
     * @param wlb First worker index in this process.
     */
    void initialize(int K,
            Comm comm,
            int size,
            int rank,
            int count,
            int wlb) {
        // Record parameters.
        this.K = K;
        this.comm = comm;
        this.size = size;
        this.rank = rank;
        this.count = count;

        // Set up worker team threads. Additional master thread in process 0.
        int WM = K + (rank == 0 ? 1 : 0);
        myThread = new WorkerTeamThread[WM];
        for (int i = 0; i < K; ++i) {
            myThread[i] = new WorkerTeamThread(this, wlb + i);
        }
        if (WM > K) {
            myThread[K] = new WorkerTeamThread(this, -1);
        }
    }

// Exported operations.
    /**
     * Execute the given worker region.
     *
     * @param theRegion Worker region.
     *
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theRegion</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * worker team is already executing a worker region. Thrown if
     * <TT>theRegion</TT> is already being executed by a worker team.
     * @exception Exception Exception thrown by the worker region's
     * <TT>start()</TT>,
     * <TT>run()</TT>, or <TT>finish()</TT> methods.
     */
    public final void execute(WorkerRegion theRegion)
            throws Exception {
        // Verify preconditions.
        if (theRegion == null) {
            throw new NullPointerException("WorkerTeam.execute(): theRegion is null");
        }
        if (myRegion != null) {
            throw new IllegalStateException("WorkerTeam.execute(): Already executing a worker region");
        }
        if (theRegion.myTeam != null) {
            throw new IllegalStateException("WorkerTeam.execute(): theRegion already being executed by a worker team");
        }

        // Record worker region.
        myRegion = theRegion;
        myExceptionMap = new ConcurrentHashMap<Integer, Throwable>(K, 0.75f, K);
        theRegion.myTeam = this;

        try {
            // Perform the worker region's start() method. Any exception aborts
            // the execute() method.
            myRegion.start();

            // Release the team threads to perform the worker region's run()
            // method.
            for (WorkerTeamThread thread : myThread) {
                thread.myRegionBeginSemaphore.release();
            }

            // Wait until all team threads have returned from the worker
            // region's run() method.
            myRegionEndSemaphore.acquireUninterruptibly(myThread.length);

            // Propagate any exceptions thrown by the run() method.
            if (myExceptionMap.isEmpty()) {
            } else if (myExceptionMap.size() == 1) {
                rethrow(myExceptionMap.values().iterator().next());
            } else {
                throw new MultipleParallelException("WorkerTeam.execute(): Multiple threads threw exceptions",
                        myExceptionMap);
            }

            // Perform the worker region's finish() method. Any exception aborts
            // the execute() method.
            myRegion.finish();
        } finally {
            // Clean up.
            myRegion.myTeam = null;
            myExceptionMap = null;
            myRegion = null;
        }
    }

    /**
     * Determine if this worker team is executing a worker region.
     *
     * @return True if this worker team is executing a worker region, false
     * otherwise.
     */
    public final boolean isExecutingInParallel() {
        return myRegion != null;
    }

    /**
     * Returns the worker region of code that this worker team is executing.
     *
     * @return Worker region.
     *
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * worker team is not executing a worker region.
     */
    public final WorkerRegion region() {
        if (myRegion == null) {
            throw new IllegalStateException("WorkerTeam.region(): Not executing a worker region");
        }
        return myRegion;
    }

    /**
     * Determine the number of worker threads in this worker team in this
     * process. This does not include the master thread if any.
     *
     * @return Number of worker threads in this process.
     */
    public final int getThreadCount() {
        return K;
    }

    /**
     * Determine the total number of worker threads in this worker team in all
     * processes. This does not include the master thread.
     *
     * @return Number of worker threads in all processes.
     */
    public final int getTotalThreadCount() {
        return count;
    }

    /**
     * Determine the rank of the process that contains the master thread. At
     * present, this is always rank 0.
     *
     * @return Master process rank.
     */
    public int masterRank() {
        return 0;
    }

    /**
     * Determine the rank of the process that contains the worker thread with
     * the given index.
     *
     * @param w Worker index.
     *
     * @return Worker process rank.
     *
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>w</TT> is not in the range 0 ..
     * <TT>getTotalThreadCount()</TT>&minus;1.
     */
    public int workerRank(int w) {
        if (0 > w || w >= count) {
            throw new IllegalArgumentException("WorkerTeam.workerRank(): w (= " + w + ") illegal");
        }
        return w;
    }

// Hidden operations.
    /**
     * Re-throw the given object as a checked or unchecked exception. If the
     * given object is null or is not throwable, do nothing.
     */
    static void rethrow(Object exc)
            throws Exception {
        if (exc instanceof RuntimeException) {
            throw (RuntimeException) exc;
        } else if (exc instanceof Exception) {
            throw (Exception) exc;
        } else if (exc instanceof Error) {
            throw (Error) exc;
        }
    }

}
