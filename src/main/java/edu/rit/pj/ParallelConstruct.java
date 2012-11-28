//******************************************************************************
//
// File:    ParallelConstruct.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ParallelConstruct
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

import edu.rit.pj.reduction.*;

/**
 * Class ParallelConstruct is the common base class for all parallel constructs
 * that are executed by a {@linkplain ParallelTeam}.
 *
 * @author  Alan Kaminsky
 * @version 20-Dec-2007
 */
public abstract class ParallelConstruct
	{

// Hidden data members.

	// 128 bytes of extra padding to avert cache interference.
	private long p0, p1, p2, p3, p4, p5, p6, p7;
	private long p8, p9, pa, pb, pc, pd, pe, pf;

	// Parallel team that is executing this parallel construct, or null if none.
	ParallelTeam myTeam;

// Exported constructors.

	/**
	 * Construct a new parallel construct.
	 */
	public ParallelConstruct()
		{
		}

// Exported operations.

	/**
	 * Determine if a parallel team is executing this parallel construct.
	 *
	 * @return  True if a parallel team is executing this parallel construct,
	 *          false otherwise.
	 */
	public final boolean isExecutingInParallel()
		{
		return myTeam != null;
		}

	/**
	 * Returns the parallel team that is executing this parallel construct.
	 *
	 * @return  Parallel team.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if no parallel team is executing this
	 *     parallel construct.
	 */
	public final ParallelTeam team()
		{
		if (myTeam == null)
			{
			throw new IllegalStateException
				("ParallelConstruct.team(): No parallel team executing");
			}
		return myTeam;
		}

	/**
	 * Returns the parallel region of code within which a parallel team is
	 * executing this parallel construct.
	 *
	 * @return  Parallel region.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if no parallel team is executing this
	 *     parallel construct.
	 */
	public final ParallelRegion region()
		{
		if (myTeam == null)
			{
			throw new IllegalStateException
				("ParallelConstruct.region(): No parallel team executing");
			}
		return myTeam.myRegion;
		}

	/**
	 * Determine the number of threads in the parallel team executing this
	 * parallel construct.
	 *
	 * @return  Number of threads in the thread team.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if no parallel team is executing this
	 *     parallel construct.
	 */
	public final int getThreadCount()
		{
		if (myTeam == null)
			{
			throw new IllegalStateException
				("ParallelConstruct.getThreadCount(): No parallel team executing");
			}
		return myTeam.K;
		}

	/**
	 * Determine the index of the calling thread in the parallel team executing
	 * this parallel construct.
	 *
	 * @return  Index of the calling thread in the range 0 ..
	 *          <TT>getThreadCount()-1</TT>.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if no parallel team is executing this
	 *     parallel construct. Thrown if the thread calling
	 *     <TT>getThreadIndex()</TT> is not part of the parallel team executing
	 *     this parallel construct.
	 */
	public final int getThreadIndex()
		{
		return getCurrentThread().myIndex;
		}

// Hidden operations.

	/**
	 * Get the parallel team thread that is calling this method.
	 *
	 * @return  Parallel team thread.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if the calling thread is not one of the
	 *     parallel team threads executing this parallel construct.
	 */
	ParallelTeamThread getCurrentThread()
		{
		if (myTeam == null)
			{
			throw new IllegalStateException
				("ParallelConstruct.getCurrentThread(): No parallel team executing");
			}
		try
			{
			ParallelTeamThread current = (ParallelTeamThread)
				Thread.currentThread();
			if (current.myTeam != this.myTeam)
				{
				throw new IllegalStateException
					("ParallelConstruct.getCurrentThread(): Current thread is not executing this parallel construct");
				}
			return current;
			}
		catch (ClassCastException exc)
			{
			throw new IllegalStateException
				("ParallelConstruct.getCurrentThread(): Current thread is not a parallel team thread",
				 exc);
			}
		}

	}
