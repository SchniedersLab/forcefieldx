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
 * @author JacobLitman
 */
public class KillRegion extends ParallelRegion {
    @Override
    public void run() throws Exception {
        // Does precisely nothing save exist.
        // Post-modern poetry to be implemented.
    }
}
