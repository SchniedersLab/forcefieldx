//******************************************************************************
//
// File:    ParallelTeamThread_0.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelTeamThread_0
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
 * Class ParallelTeamThread_0 provides one thread in a {@linkplain ParallelTeam}
 * of threads for executing a {@linkplain ParallelRegion} in parallel. This
 * class is used for thread index 0.
 *
 * @author  Alan Kaminsky
 * @version 03-Oct-2007
 */
class ParallelTeamThread_0
	extends ParallelTeamThread
	{

// Exported constructors.

	/**
	 * Construct a new parallel team thread.
	 *
	 * @param  theTeam   Parallel team to which this thread belongs.
	 * @param  theIndex  Index of this thread within the team.
	 */
	public ParallelTeamThread_0
		(ParallelTeam theTeam,
		 int theIndex)
		{
		super (theTeam, theIndex);
		}

// Hidden operations.

	/**
	 * Do this thread's portion of a barrier with no barrier action. This method
	 * is called by thread 0 of the parallel team.
	 */
	void barrier()
		{
		myTeam.barrier();
		}

	/**
	 * Do this thread's portion of a barrier with a barrier action. This method
	 * is called by thread 0 of the parallel team.
	 *
	 * @param  action  Barrier action.
	 *
	 * @exception  Exception
	 *     Thrown if the <TT>action</TT>'s <TT>run()</TT> method throws an
	 *     exception.
	 */
	void barrier
		(BarrierAction action)
		throws Exception
		{
		myTeam.barrier (action);
		}

	}
