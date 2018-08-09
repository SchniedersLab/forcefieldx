//******************************************************************************
//
// File:    ParallelTeam.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelTeam
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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

//******************************************************************************
// File modified 10/8/2014 by Jacob Litman to enable garbage collection of
// ParallelTeamThreads.
//******************************************************************************
package edu.rit.pj;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Class ParallelTeam provides a team of threads for executing a {@linkplain
 * ParallelRegion} in parallel.
 * <P>
 * To execute a parallel region, create a ParallelTeam object; create an
 * instance of a concrete subclass of class {@linkplain ParallelRegion}; and
 * pass this instance to the parallel team's <TT>execute()</TT> method. For
 * further information, see class {@linkplain ParallelRegion}.
 *
 * @author Alan Kaminsky
 * @version 19-May-2008
 */
public class ParallelTeam {

// Hidden data members.
    // Number of threads.
    int K;

    // Array of threads in the team.
    ParallelTeamThread[] myThread;

    // Parallel region being executed, or null if none is being executed.
    ParallelRegion myRegion;

    // Semaphore for synchronizing threads at the end of a parallel region.
    Semaphore myRegionEndSemaphore = new Semaphore(0);

    // Exception map for parallel region, or null if none is being executed.
    ConcurrentHashMap<Integer, Throwable> myExceptionMap;

    // Team barrier flag. Used by the ParallelRegion.barrier() method.
    volatile int myBarrierFlag;

    // Parallel construct counter. Counts how many parallel constructs have been
    // encountered.
    AtomicInteger myConstructCount = new AtomicInteger(0);
    
    // Set false if the ParallelTeam is shut down.
    boolean isActive = true;

// Exported constructors.
    /**
     * Construct a new parallel team with the default number of threads. If the
     * <TT>"pj.nt"</TT> Java property is specified, that property gives the
     * default number of threads, which must be an integer greater than or equal
     * to 1. If the <TT>"pj.nt"</TT> Java property is not specified, the default
     * number of threads is the value returned by the
     * <TT>Runtime.availableProcessors()</TT> method. You can specify the
     * default number of threads on the Java command line like this:
     * <PRE>
     *     java -Dpj.nt=4 . . .
     * </PRE>
     *
     * @exception IllegalArgumentException (unchecked exception) Thrown if the
     * <TT>"pj.nt"</TT> property value is not an integer greater than or equal
     * to 1.
     */
    public ParallelTeam() {
        this(getDefaultThreadCount());
    }

    /**
     * Construct a new parallel team with the given number of threads.
     *
     * @param K Number of threads.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <I>K</I> is less than 1.
     */
    public ParallelTeam(int K) {
        if (K < 1) {
            throw new IllegalArgumentException("ParallelTeam(): K must be >= 1");
        }
        this.K = K;

        myThread = new ParallelTeamThread[K];
        myThread[0] = new ParallelTeamThread_0(this, 0);
        for (int i = 1; i < K; ++i) {
            myThread[i] = new ParallelTeamThread(this, i);
        }
    }

// Exported operations.
    /**
     * Execute the given parallel region.
     *
     * @param theRegion Parallel region.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <TT>theRegion</TT> is null.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * parallel team is already executing a parallel region. Thrown if
     * <TT>theRegion</TT> is already being executed by a parallel team.
     * @exception Exception Exception thrown by the parallel region's
     * <TT>start()</TT>,
     * <TT>run()</TT>, or <TT>finish()</TT> methods.
     * @throws java.lang.Exception if any.
     */
    public final void execute(ParallelRegion theRegion)
            throws Exception {
        // Verify preconditions.
        if (theRegion == null) {
            throw new NullPointerException("ParallelTeam.execute(): theRegion is null");
        }
        if (myRegion != null) {
            throw new IllegalStateException("ParallelTeam.execute(): Already executing a parallel region");
        }
        if (theRegion.myTeam != null) {
            throw new IllegalStateException("ParallelTeam.execute(): theRegion already being executed by a parallel team");
        }
        if (!isActive) {
            throw new IllegalStateException("ParallelTeam.execute(): The team has been shut down.");
        }

        // Record parallel region.
        myRegion = theRegion;
        myExceptionMap = new ConcurrentHashMap<Integer, Throwable>(K, 0.75f, K);
        theRegion.myTeam = this;

        try {
            // Perform the parallel region's start() method. Any exception
            // aborts the execute() method.
            myRegion.start();

            // Release the team threads to perform the parallel region's run()
            // method.
            for (ParallelTeamThread thread : myThread) {
                thread.myRegionBeginSemaphore.release();
            }

            // Wait until all team threads have returned from the parallel
            // region's run() method.
            myRegionEndSemaphore.acquireUninterruptibly(K);

            // Propagate any exceptions thrown by the run() method.
            if (myExceptionMap.isEmpty()) {
            } else if (myExceptionMap.size() == 1) {
                rethrow(myExceptionMap.values().iterator().next());
            } else {
                throw new MultipleParallelException("ParallelTeam.execute(): Multiple threads threw exceptions",
                        myExceptionMap);
            }

            // Perform the parallel region's finish() method. Any exception
            // aborts the execute() method.
            myRegion.finish();
        } finally {
            // Clean up.
            myRegion.myTeam = null;
            myExceptionMap = null;
            myRegion = null;
        }
    }

    /**
     * Determine if this parallel team is executing a parallel region.
     *
     * @return True if this parallel team is executing a parallel region, false
     * otherwise.
     */
    public final boolean isExecutingInParallel() {
        return myRegion != null;
    }

    /**
     * Returns the parallel region of code that this parallel team is executing.
     *
     * @return Parallel region.
     * @exception IllegalStateException (unchecked exception) Thrown if this
     * parallel team is not executing a parallel region.
     */
    public final ParallelRegion region() {
        if (myRegion == null) {
            throw new IllegalStateException("ParallelTeam.region(): Not executing a parallel region");
        }
        return myRegion;
    }

    /**
     * Determine the number of threads in this parallel team.
     *
     * @return Number of threads in the team.
     */
    public final int getThreadCount() {
        return K;
    }

    /**
     * Determine the default number of threads for a parallel team. If the
     * <TT>"pj.nt"</TT> Java property is specified, that property gives the
     * default number of threads, which must be an integer greater than or equal
     * to 1. If the <TT>"pj.nt"</TT> Java property is not specified, the default
     * number of threads is the value returned by the
     * <TT>Runtime.availableProcessors()</TT> method. You can specify the
     * default number of threads on the Java command line like this:
     * <PRE>
     *     java -Dpj.nt=4 . . .
     * </PRE>
     *
     * @return Default number of threads for a parallel team.
     * @exception IllegalArgumentException (unchecked exception) Thrown if the
     * <TT>"pj.nt"</TT> property value is not an integer greater than or equal
     * to 1.
     */
    public static int getDefaultThreadCount() {
        int k = PJProperties.getPjNt();
        if (k == 0) {
            k = Runtime.getRuntime().availableProcessors();
        }
        return k;
    }

// Hidden operations.
    /**
     * Do the thread-0 portion of a barrier with no barrier action. This method
     * is called by thread 0 of the parallel team.
     */
    void barrier() {
        // Get the new team barrier flag.
        int newBarrierFlag = myBarrierFlag ^ 1;

        // Wait until each team thread 1 .. K-1 has switched to the new
        // barrier flag.
        for (int i = 1; i < K; ++i) {
            ParallelTeamThread thread_i = myThread[i];
            if (thread_i.myBarrierFlag != newBarrierFlag) {
                Spinner spinner = new Spinner();
                while (thread_i.myBarrierFlag != newBarrierFlag) {
                    spinner.spin();
                }
            }
        }

        // Switch to the new team barrier flag.
        myBarrierFlag = newBarrierFlag;
    }

    /**
     * Do the thread-0 portion of a barrier with a barrier action. This method
     * is called by thread 0 of the parallel team.
     *
     * @param action Barrier action.
     *
     * @exception Exception Thrown if the <TT>action</TT>'s <TT>run()</TT>
     * method throws an exception.
     */
    void barrier(BarrierAction action)
            throws Exception {
        // Get the new team barrier flag.
        int newBarrierFlag = myBarrierFlag ^ 1;

        // Wait until each team thread 1 .. K-1 has switched to the new
        // barrier flag.
        for (int i = 1; i < K; ++i) {
            ParallelTeamThread thread_i = myThread[i];
            if (thread_i.myBarrierFlag != newBarrierFlag) {
                Spinner spinner = new Spinner();
                while (thread_i.myBarrierFlag != newBarrierFlag) {
                    spinner.spin();
                }
            }
        }

        try {
            // Do the barrier action.
            action.myTeam = this;
            action.run();
        } finally {
            action.myTeam = null;

            // Switch to the new team barrier flag.
            myBarrierFlag = newBarrierFlag;
        }
    }

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
    
    /**
     * Kills the team's threads run() methods so that they are no longer GC roots.
     * Useful if you are repetitively constructing ParallelTeam objects, although
     * it is slightly more elegant just to keep the same ParallelTeam objects through
     * the entire execution of a program.
     *
     * @throws java.lang.Exception if any.
     */
    public void shutdown() throws Exception {
        if (isActive) {
            this.execute(new KillRegion());
            isActive = false;
        }
    }
}
