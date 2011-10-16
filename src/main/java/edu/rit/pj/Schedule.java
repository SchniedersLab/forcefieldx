//******************************************************************************
//
// File:    Schedule.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.Schedule
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
 * Class Schedule provides an object that determines how to schedule the
 * iterations of a {@linkplain ParallelForLoop} among the threads in a
 * {@linkplain ParallelTeam}. To create a schedule object, call one of the
 * following static methods in class {@linkplain IntegerSchedule} or class
 * {@linkplain LongSchedule}:
 * <UL>
 * <LI><TT>IntegerSchedule.fixed()</TT>
 * <LI><TT>IntegerSchedule.dynamic()</TT>
 * <LI><TT>IntegerSchedule.guided()</TT>
 * <LI><TT>IntegerSchedule.runtime()</TT>
 * <LI><TT>IntegerSchedule.parse()</TT>
 * <LI><TT>LongSchedule.fixed()</TT>
 * <LI><TT>LongSchedule.dynamic()</TT>
 * <LI><TT>LongSchedule.guided()</TT>
 * <LI><TT>LongSchedule.runtime()</TT>
 * <LI><TT>LongSchedule.parse()</TT>
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 20-Dec-2007
 */
public abstract class Schedule
	{

// Hidden data members.

	// 128 bytes of extra padding to avert cache interference.
	private long p0, p1, p2, p3, p4, p5, p6, p7;
	private long p8, p9, pa, pb, pc, pd, pe, pf;

	// True to break out of the parallel for loop.
	boolean myBreak;

// Hidden constructors.

	/**
	 * Construct a new schedule object.
	 */
	Schedule()
		{
		}

	}
