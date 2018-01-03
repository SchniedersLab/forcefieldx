//******************************************************************************
//
// File:    WorkerTeamThread.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.WorkerTeamThread
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

import java.util.concurrent.Semaphore;

/**
 * Class WorkerTeamThread provides one thread in a {@linkplain WorkerTeam} of
 * threads for executing a {@linkplain WorkerRegion} in parallel.
 *
 * @author Alan Kaminsky
 * @version 19-Jan-2010
 */
class WorkerTeamThread
        extends Thread {

// Hidden data members.
    // Reference to the worker team.
    WorkerTeam myTeam;

    // Index of this thread within the worker team.
    int myIndex;

    // Semaphore for synchronizing threads at the beginning of a worker region.
    Semaphore myRegionBeginSemaphore = new Semaphore(0);

    // 128 bytes of extra padding to avert cache interference.
    private long p0, p1, p2, p3, p4, p5, p6, p7;
    private long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.
    /**
     * Construct a new worker team thread.
     *
     * @param theTeam Worker team to which this thread belongs.
     * @param theIndex Index of this thread within the team.
     */
    public WorkerTeamThread(WorkerTeam theTeam,
            int theIndex) {
        myTeam = theTeam;
        myIndex = theIndex;
        setDaemon(true);
        start();
    }

// Exported operations.
    /**
     * Run this worker team thread.
     */
    public void run() {
        for (;;) {
            // Wait until released by the main thread.
            myRegionBeginSemaphore.acquireUninterruptibly();

            // Call the worker region's run() method. Save any
            // exception for later.
            try {
                myTeam.myRegion.run();
            } catch (Throwable exc) {
                synchronized (System.err) {
                    System.err.println("Worker team thread " + myIndex
                            + ": WorkerRegion.run() threw an exception");
                    exc.printStackTrace(System.err);
                }
                myTeam.myExceptionMap.put(myIndex, exc);
            }

            // Tell the main thread we're done.
            myTeam.myRegionEndSemaphore.release();
        }
    }

}
