//******************************************************************************
//
// File:    KillRegion.java
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
// Additional file added 10/8/2014 by Jacob Litman to address the inability of 
// the Java garbage collector to collect ParallelTeamThreads, whose run methods 
// are infinite loops: since running threads are always GC roots, programs which
// set up many ParallelTeam objects will accumulate inactive threads and crash.
//
// Other modifications in ParallelTeam and ParallelTeamThread. 
//******************************************************************************
package edu.rit.pj;

/**
 * Provides a mechanism to shut down a ParallelTeam's threads, enabling garbage
 * collection.
 *
 * @author Jacob Litman
 */
public class KillRegion extends ParallelRegion {

    /** {@inheritDoc} */
    @Override
    public void run() throws Exception {
        // Empty run method.
        ParallelTeamThread currentThread = getCurrentThread();
        Thread thread = currentThread.getThread();
        long totalThreads = ParallelTeamThread.totalThreads.getAndDecrement();
        // System.out.printf(" Killing team %s thread %s of %d\n", currentThread.myTeam, thread, totalThreads);
    }
}
