//******************************************************************************
//
// File:    HybridTeam.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.HybridTeam
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

import java.io.IOException;

import edu.rit.mp.IntegerBuf;
import edu.rit.util.Range;

/**
 * Class HybridTeam provides a team of threads, distributed across the processes
 * of a cluster parallel program, for executing a {@linkplain WorkerRegion} in
 * parallel.
 * <P>
 * A hybrid team uses a communicator for message passing. The communicator is
 * specified as a constructor argument; if not specified, the world communicator
 * is used. Every process that is part of the communicator must create the
 * hybrid team. In class HybridTeam, there are one or more worker threads per
 * process; typically, the number of threads in each process equals the number
 * of CPUs on the node executing the process. (To get just one worker thread per
 * process, use class {@linkplain WorkerTeam}.) Every worker thread in every
 * process has a unique index, going from index 0 for the first thread in the
 * first process to index <I>K</I>&minus;1 for the last thread in the last
 * process, where <I>K</I> is the total number of worker threads in all the
 * processes. In process rank 0, there is an additional master thread.
 * <P>
 * When a hybrid team is constructed, the processes in the communicator send
 * messages amongst themselves to make every process aware of the number of
 * worker threads in each process. These messages use a tag of
 * <TT>Integer.MIN_VALUE</TT>. If an I/O error occurs during this message
 * passing, the constructor throws an IOException.
 * <P>
 * To execute a worker region, create a HybridTeam object; create an instance of
 * a concrete subclass of class {@linkplain WorkerRegion}; and pass this
 * instance to the hybrid team's <TT>execute()</TT> method. For further
 * information, see class {@linkplain WorkerRegion}.
 *
 * @author  Alan Kaminsky
 * @version 18-Nov-2009
 */
public class HybridTeam
	extends WorkerTeam
	{

// Hidden data members.

	// Array of lowest worker indexes, indexed by process rank.
	int[] lowestIndex;

// Exported constructors.

	/**
	 * Construct a new hybrid team with the default number of threads per
	 * process and using the world communicator for message passing. If the
	 * <TT>"pj.nt"</TT> Java property is specified, that property gives the
	 * default number of threads per process, which must be an integer greater
	 * than or equal to 1. If the <TT>"pj.nt"</TT> Java property is not
	 * specified, the default number of threads in each process is the value
	 * returned by the <TT>Runtime.availableProcessors()</TT> method. You can
	 * specify the default number of threads per process on the Java command
	 * line like this:
	 * <PRE>
	 *     java -Dpj.nt=4 . . .
	 * </PRE>
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the <TT>"pj.nt"</TT> property value
	 *     is not an integer greater than or equal to 1.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public HybridTeam()
		throws IOException
		{
		this (getDefaultThreadCount(), Comm.world());
		}

	/**
	 * Construct a new hybrid team with the given number of threads per
	 * process and using the world communicator for message passing.
	 *
	 * @param  K  Number of threads per process.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>K</I> is less than 1.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public HybridTeam
		(int K)
		throws IOException
		{
		this (K, Comm.world());
		}

	/**
	 * Construct a new hybrid team with the default number of threads per
	 * process and using the given communicator for message passing. If the
	 * <TT>"pj.nt"</TT> Java property is specified, that property gives the
	 * default number of threads per process, which must be an integer greater
	 * than or equal to 1. If the <TT>"pj.nt"</TT> Java property is not
	 * specified, the default number of threads in each process is the value
	 * returned by the <TT>Runtime.availableProcessors()</TT> method. You can
	 * specify the default number of threads per process on the Java command
	 * line like this:
	 * <PRE>
	 *     java -Dpj.nt=4 . . .
	 * </PRE>
	 *
	 * @param  comm  Communicator to use for message passing.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the <TT>"pj.nt"</TT> property value
	 *     is not an integer greater than or equal to 1.
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>comm</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public HybridTeam
		(Comm comm)
		throws IOException
		{
		this (getDefaultThreadCount(), comm);
		}

	/**
	 * Construct a new hybrid team with the given number of threads per
	 * process and using the given communicator for message passing.
	 *
	 * @param  K     Number of threads per process.
	 * @param  comm  Communicator to use for message passing.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <I>K</I> is less than 1.
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>comm</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public HybridTeam
		(int K,
		 Comm comm)
		throws IOException
		{
		// Construct superclass but don't initialize it yet.
		super (false);

		// Verify preconditions.
		if (K <= 0)
			{
			throw new IllegalArgumentException
				("HybridTeam(): K = "+K+" illegal");
			}
		if (comm == null)
			{
			throw new NullPointerException
				("HybridTeam(): comm is null");
			}

		// Communicate with the other processes to determine the number of
		// worker threads in each process; store these in lowestIndex.
		int size = comm.size();
		int rank = comm.rank();
		lowestIndex = new int [size + 1];
		lowestIndex[rank] = K;
		IntegerBuf[] bufs =
			IntegerBuf.sliceBuffers
				(lowestIndex,
				 new Range (0, size - 1) .subranges (size));
		IntegerBuf myBuf = bufs[rank];
		comm.allGather (Integer.MIN_VALUE, myBuf, bufs);

		// Scan lowestIndex to determine the lowest worker index in each process
		// and the total number of worker threads.
		int count = 0;
		for (int i = 0; i < size; ++ i)
			{
			int tmp = lowestIndex[i];
			lowestIndex[i] = count;
			count += tmp;
			}
		lowestIndex[size] = count;

		// Initialize superclass.
		initialize
			(/*K    */ K,
			 /*comm */ comm,
			 /*size */ size,
			 /*rank */ rank,
			 /*count*/ count,
			 /*wlb  */ lowestIndex[rank]);
		}

// Exported operations.

	/**
	 * Determine the default number of threads per process for a hybrid team. If
	 * the <TT>"pj.nt"</TT> Java property is specified, that property gives the
	 * default number of threads per process, which must be an integer greater
	 * than or equal to 1. If the <TT>"pj.nt"</TT> Java property is not
	 * specified, the default number of threads in each process is the value
	 * returned by the <TT>Runtime.availableProcessors()</TT> method. You can
	 * specify the default number of threads per process on the Java command
	 * line like this:
	 * <PRE>
	 *     java -Dpj.nt=4 . . .
	 * </PRE>
	 *
	 * @return  Default number of threads per process for a hybrid team.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if the <TT>"pj.nt"</TT> property value
	 *     is not an integer greater than or equal to 1.
	 */
	public static int getDefaultThreadCount()
		{
		int k = PJProperties.getPjNt();
		if (k == 0) k = Runtime.getRuntime().availableProcessors();
		return k;
		}

	/**
	 * Determine the rank of the process that contains the worker thread with
	 * the given index.
	 *
	 * @param  w  Worker index.
	 *
	 * @return  Worker process rank.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>w</TT> is not in the range 0 ..
	 *     <TT>getTotalThreadCount()</TT>&minus;1.
	 */
	public int workerRank
		(int w)
		{
		// Verify preconditions.
		if (0 > w || w >= count)
			{
			throw new IllegalArgumentException
				("HybridTeam.workerRank(): w (= "+w+") illegal");
			}

		// Do a binary search to find w in lowestIndex.
		int L = 0;
		int U = size;
		while (U - L > 1)
			{
			// Invariant: w >= lowestIndex[L] and w < lowestIndex[U].
			int M = (L + U)/2;
			if (w >= lowestIndex[M])
				{
				L = M;
				}
			else
				{
				U = M;
				}
			}

		// Assert: L = U-1 and w >= lowestIndex[L] and w < lowestIndex[U].
		return L;
		}

	}
