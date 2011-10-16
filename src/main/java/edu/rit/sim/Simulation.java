//******************************************************************************
//
// File:    Simulation.java
// Package: edu.rit.sim
// Unit:    Class edu.rit.sim.Simulation
//
// This Java source file is copyright (C) 2011 by Alan Kaminsky. All rights
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

package edu.rit.sim;

/**
 * Class Simulation provides a discrete event simulation. To write a discrete
 * event simulation program:
 * <OL TYPE=1>
 * <P><LI>
 * Create a Simulation object.
 * <P><LI>
 * Create one or more {@linkplain Event}s and add them to the simulation (by
 * calling the <TT>doAt()</TT> or <TT>doAfter()</TT> methods).
 * <P><LI>
 * Run the simulation (by calling the <TT>run()</TT> method). The simulation
 * performs events, by calling each event's <TT>perform()</TT> method, in order
 * according to the events' simulation times, as returned by each event's
 * <TT>time()</TT> method. Performing an event may cause further events to be
 * created and added to the simulation.
 * <P><LI>
 * When there are no more events, the simulation is finished. At this point the
 * simulation's <TT>run()</TT> method returns.
 * </OL>
 *
 * @author  Alan Kaminsky
 * @version 29-Jul-2011
 */
public class Simulation
	{

// Hidden data members.

	// Minimum-priority queue of events. Uses a heap data structure. The entry
	// at index 0 is a sentinel with time = 0.0.
	private Event[] heap = new Event [1024];

	// Number of entries in the heap (including the sentinel).
	private int N = 1;

	// Simulation time.
	private double T = 0.0;

// Exported constructors.

	/**
	 * Construct a new simulation.
	 */
	public Simulation()
		{
		heap[0] = new Event() { public void perform() { } };
		heap[0].sim = this;
		heap[0].time = 0.0;
		}

// Exported operations.

	/**
	 * Returns the current simulation time.
	 *
	 * @return  Simulation time.
	 */
	public double time()
		{
		return T;
		}

	/**
	 * Schedule the given event to be performed at the given time in this
	 * simulation.
	 *
	 * @param  t      Simulation time for <TT>event</TT>.
	 * @param  event  Event to be performed.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>t</TT> is less than the current
	 *     simulation time.
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>event</TT> is null.
	 */
	public void doAt
		(double t,
		 Event event)
		{
		// Verify preconditions.
		if (t < T)
			{
			throw new IllegalArgumentException
				("Simulation.doAt(): t = "+t+" less than simulation time ="+T+
				 ", illegal");
			}
		if (event == null)
			{
			throw new NullPointerException
				("Simulation.doAt(): event = null");
			}

		// Set event fields.
		event.sim = this;
		event.time = t;

		// Grow heap if necessary.
		if (N == heap.length)
			{
			Event[] newheap = new Event [N + 1024];
			System.arraycopy (heap, 0, newheap, 0, N);
			heap = newheap;
			}

		// Insert event into heap in min-priority order.
		heap[N] = event;
		siftUp (N);
		++ N;
		}

	/**
	 * Schedule the given event to be performed at a time <TT>dt</TT> in the
	 * future (at current simulation time + <TT>dt</TT>) in this simulation.
	 *
	 * @param  dt     Simulation time delta for <TT>event</TT>.
	 * @param  event  Event to be performed.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>dt</TT> is less than zero.
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>event</TT> is null.
	 */
	public void doAfter
		(double dt,
		 Event event)
		{
		doAt (T + dt, event);
		}

	/**
	 * Run the simulation. At the start of the simulation, the simulation time
	 * is 0. The <TT>run()</TT> method returns when there are no more events.
	 */
	public void run()
		{
		while (N > 1)
			{
			// Extract minimum event from heap.
			Event event = heap[1];
			-- N;
			heap[1] = heap[N];
			heap[N] = null;
			if (N > 1) siftDown (1);

			// Advance simulation time and perform event.
			T = event.time;
			event.perform();
			}
		}

// Hidden operations.

	/**
	 * Sift up the heap entry at the given index.
	 *
	 * @param  c  Index.
	 */
	private void siftUp
		(int c)
		{
		double c_time = heap[c].time;
		int p = c >> 1;
		double p_time = heap[p].time;
		while (c_time < p_time)
			{
			Event temp = heap[c];
			heap[c] = heap[p];
			heap[p] = temp;
			c = p;
			p = c >> 1;
			p_time = heap[p].time;
			}
		}

	/**
	 * Sift down the heap entry at the given index.
	 *
	 * @param  p  Index.
	 */
	private void siftDown
		(int p)
		{
		double p_time = heap[p].time;
		int lc = (p << 1);
		double lc_time = lc < N ? heap[lc].time : Double.POSITIVE_INFINITY;
		int rc = (p << 1) + 1;
		double rc_time = rc < N ? heap[rc].time : Double.POSITIVE_INFINITY;
		int c;
		double c_time;
		if (lc_time < rc_time)
			{
			c = lc;
			c_time = lc_time;
			}
		else
			{
			c = rc;
			c_time = rc_time;
			}
		while (c_time < p_time)
			{
			Event temp = heap[c];
			heap[c] = heap[p];
			heap[p] = temp;
			p = c;
			lc = (p << 1);
			lc_time = lc < N ? heap[lc].time : Double.POSITIVE_INFINITY;
			rc = (p << 1) + 1;
			rc_time = rc < N ? heap[rc].time : Double.POSITIVE_INFINITY;
			if (lc_time < rc_time)
				{
				c = lc;
				c_time = lc_time;
				}
			else
				{
				c = rc;
				c_time = rc_time;
				}
			}
		}

	}
