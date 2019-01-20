//******************************************************************************
//
// File:    ParallelTeamThread.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelTeamThread
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

//******************************************************************************
// File modified 10/8/2014 by Jacob Litman to enable garbage collection of
// ParallelTeamThreads.
//******************************************************************************
package edu.rit.pj;

import java.util.concurrent.Semaphore;

/**
 * Class ParallelTeamThread provides one thread in a {@linkplain ParallelTeam}
 * of threads for executing a {@linkplain ParallelRegion} in parallel.
 *
 * @author Alan Kaminsky
 * @version 20-Dec-2007
 */
class ParallelTeamThread
        extends Thread {

// Hidden data members.
    // Reference to the parallel team.
    ParallelTeam myTeam;

    // Index of this thread within the parallel team.
    int myIndex;

    // Semaphore for synchronizing threads at the beginning of a parallel
    // region.
    Semaphore myRegionBeginSemaphore = new Semaphore(0);

    // Thread barrier flag. Used by the ParallelRegion.barrier() method.
    volatile int myBarrierFlag;

    // Parallel construct counter. Counts how many times this thread has arrived
    // at the top of a parallel construct.
    volatile int myConstructCount;

    // IntegerSchedule object for a parallel for loop.
    volatile IntegerSchedule myIntegerSchedule;

    // LongSchedule object for a parallel for loop.
    volatile LongSchedule myLongSchedule;

    // Item generator object for a parallel iteration.
    volatile ItemGenerator<?> myItemGenerator;

    // Exception thrown while setting up a parallel construct, or null if none.
    volatile Throwable myConstructException;

    // 128 bytes of extra padding to avert cache interference.
    private long p0, p1, p2, p3, p4, p5, p6, p7;
    private long p8, p9, pa, pb, pc, pd, pe, pf;

// Exported constructors.
    /**
     * Construct a new parallel team thread.
     *
     * @param theTeam Parallel team to which this thread belongs.
     * @param theIndex Index of this thread within the team.
     */
    public ParallelTeamThread(ParallelTeam theTeam,
            int theIndex) {
        myTeam = theTeam;
        myIndex = theIndex;
        setDaemon(true);
        start();
    }

// Exported operations.
    /**
     * Run this parallel team thread.
     */
    public void run() {
        for (;;) {
            // Wait until released by the main thread.
            myRegionBeginSemaphore.acquireUninterruptibly();

            if (myTeam.myRegion instanceof KillRegion) {
                try {
                    myTeam.myRegion.run();
                } catch (Exception ex) {
                    System.err.println(String.format("Error exiting parallel "
                            + "team thread: %s", ex.toString()));
                }
                break;
            }
            
            // Call the parallel region's run() method. Save any
            // exception for later.
            try {
                myTeam.myRegion.run();
            } catch (Throwable exc) {
                synchronized (System.err) {
                    System.err.println("Parallel team thread " + myIndex
                            + ": ParallelRegion.run() threw an exception");
                    exc.printStackTrace(System.err);
                    if (exc.getCause() != null) {
                        exc.getCause().printStackTrace(System.err);
                    }
                }
                myTeam.myExceptionMap.put(myIndex, exc);
            }

            // Tell the main thread we're done.
            myTeam.myRegionEndSemaphore.release();
        }
        myTeam.myRegionEndSemaphore.release();
    }

// Hidden operations.
    /**
     * Do this thread's portion of a barrier with no barrier action. This method
     * is called by thread 1 through thread K-1 of the parallel team.
     */
    void barrier() {
        // Get the new team barrier flag.
        int newBarrierFlag = myTeam.myBarrierFlag ^ 1;

        // Switch this thread to the new barrier flag.
        myBarrierFlag = newBarrierFlag;

        // Wait until thread 0 has switched to the new team barrier flag.
        if (myTeam.myBarrierFlag != newBarrierFlag) {
            Spinner spinner = new Spinner();
            while (myTeam.myBarrierFlag != newBarrierFlag) {
                spinner.spin();
            }
        }
    }

    /**
     * Do this thread's portion of a barrier with a barrier action. This method
     * is called by thread 1 through thread K-1 of the parallel team.
     *
     * @param action Barrier action.
     *
     * @exception Exception Thrown if the <code>action</code>'s <code>run()</code>
     * method throws an exception.
     */
    void barrier(BarrierAction action)
            throws Exception {
        barrier();
    }

    /**
     * Do processing when this thread arrives at the top of a parallel
     * construct.
     *
     * @return True if this thread is the first to arrive, false otherwise.
     */
    boolean arriveAtParallelConstruct() {
        // Wait until every thread's construct count is greater than or equal to
        // this thread's construct count.
        int K = myTeam.K;
        for (int i = 0; i < K; ++i) {
            ParallelTeamThread thread_i = myTeam.myThread[i];
            if (thread_i.myConstructCount < this.myConstructCount) {
                Spinner spinner = new Spinner();
                while (thread_i.myConstructCount < this.myConstructCount) {
                    spinner.spin();
                }
            }
        }

        // Determine if this thread is the first to arrive.
        for (;;) {
            int oldConstructCount = myTeam.myConstructCount.get();
            if (oldConstructCount != myConstructCount) {
                // This thread is not the first to arrive.
                myConstructCount = oldConstructCount;
                return false;
            }
            // This thread may be the first to arrive.
            if (myTeam.myConstructCount.compareAndSet(oldConstructCount, oldConstructCount + 1)) {
                // This thread is the first to arrive.
                myConstructCount = oldConstructCount + 1;
                return true;
            }
        }
    }

    /**
     * Set this thread's IntegerSchedule.
     *
     * @param theSchedule Integer schedule.
     */
    void setIntegerSchedule(IntegerSchedule theSchedule) {
        // Wait until myIntegerSchedule is null.
        if (myIntegerSchedule != null) {
            Spinner spinner = new Spinner();
            while (myIntegerSchedule != null) {
                spinner.spin();
            }
        }

        myIntegerSchedule = theSchedule;
    }

    /**
     * Get this thread's IntegerSchedule.
     *
     * @return IntegerSchedule.
     *
     * @exception Exception This method can throw any exception.
     */
    IntegerSchedule getIntegerSchedule()
            throws Exception {
        // Wait until myIntegerSchedule or myConstructException is set to a
        // non-null value.
        if (myIntegerSchedule == null && myConstructException == null) {
            Spinner spinner = new Spinner();
            while (myIntegerSchedule == null && myConstructException == null) {
                spinner.spin();
            }
        }

        // Make temporary copies and null out originals.
        IntegerSchedule schedule = myIntegerSchedule;
        myIntegerSchedule = null;
        Throwable exc = myConstructException;
        myConstructException = null;

        // Re-throw exception if necessary.
        ParallelTeam.rethrow(exc);

        return schedule;
    }

    /**
     * Set this thread's LongSchedule.
     *
     * @param theSchedule Long schedule.
     */
    void setLongSchedule(LongSchedule theSchedule) {
        // Wait until myLongSchedule is null.
        if (myLongSchedule != null) {
            Spinner spinner = new Spinner();
            while (myLongSchedule != null) {
                spinner.spin();
            }
        }

        myLongSchedule = theSchedule;
    }

    /**
     * Get this thread's LongSchedule.
     *
     * @return LongSchedule.
     *
     * @exception Exception This method can throw any exception.
     */
    LongSchedule getLongSchedule()
            throws Exception {
        // Wait until myLongSchedule or myConstructException is set to a
        // non-null value.
        if (myLongSchedule == null && myConstructException == null) {
            Spinner spinner = new Spinner();
            while (myLongSchedule == null && myConstructException == null) {
                spinner.spin();
            }
        }

        // Make temporary copies and null out originals.
        LongSchedule schedule = myLongSchedule;
        myLongSchedule = null;
        Throwable exc = myConstructException;
        myConstructException = null;

        // Re-throw exception if necessary.
        ParallelTeam.rethrow(exc);

        return schedule;
    }

    /**
     * Set this thread's ItemGenerator.
     *
     * @param theItemGenerator Item generator.
     */
    void setItemGenerator(ItemGenerator<?> theItemGenerator) {
        // Wait until myItemGenerator is null.
        if (myItemGenerator != null) {
            Spinner spinner = new Spinner();
            while (myItemGenerator != null) {
                spinner.spin();
            }
        }

        myItemGenerator = theItemGenerator;
    }

    /**
     * Get this thread's ItemGenerator.
     *
     * @return Item generator.
     *
     * @exception Exception This method can throw any exception.
     */
    ItemGenerator<?> getItemGenerator()
            throws Exception {
        // Wait until myItemGenerator or myConstructException is set to a
        // non-null value.
        if (myItemGenerator == null && myConstructException == null) {
            Spinner spinner = new Spinner();
            while (myItemGenerator == null && myConstructException == null) {
                spinner.spin();
            }
        }

        // Make temporary copies and null out originals.
        ItemGenerator<?> generator = myItemGenerator;
        myItemGenerator = null;
        Throwable exc = myConstructException;
        myConstructException = null;

        // Re-throw exception if necessary.
        ParallelTeam.rethrow(exc);

        return generator;
    }

    /**
     * Set this thread's construct exception.
     *
     * @param theException Construct exception.
     */
    void setConstructException(Throwable theException) {
        // Wait until myConstructException is null.
        if (myConstructException != null) {
            Spinner spinner = new Spinner();
            while (myConstructException != null) {
                spinner.spin();
            }
        }

        myConstructException = theException;
    }

}
