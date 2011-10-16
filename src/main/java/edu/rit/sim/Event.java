//******************************************************************************
//
// File:    Event.java
// Package: edu.rit.sim
// Unit:    Class edu.rit.sim.Event
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
 * Class Event is the abstract base class for an event in a discrete event
 * simulation.
 *
 * @author  Alan Kaminsky
 * @version 29-Jul-2011
 */
public abstract class Event
	{

// Hidden data members.

	// Simulation in which this event occurs.
	Simulation sim;

	// Simulation time of this event.
	double time;

// Exported constructors.

	/**
	 * Construct a new event.
	 */
	public Event()
		{
		}

// Exported operations.

	/**
	 * Returns the simulation in which this event occurs.
	 *
	 * @return  Simulation.
	 */
	public final Simulation simulation()
		{
		return sim;
		}

	/**
	 * Returns this event's simulation time, the time when this event is
	 * scheduled to take place.
	 *
	 * @return  Simulation time.
	 */
	public final double time()
		{
		return time;
		}

	/**
	 * Schedule the given event to be performed at the given time in this
	 * event's simulation.
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
	public final void doAt
		(double t,
		 Event event)
		{
		sim.doAt (t, event);
		}

	/**
	 * Schedule the given event to be performed at a time <TT>dt</TT> in the
	 * future (at current simulation time + <TT>dt</TT>) in this event's
	 * simulation.
	 *
	 * @param  dt     Simulation time delta for <TT>event</TT>.
	 * @param  event  Event to be performed.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>dt</TT> is less than zero.
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>event</TT> is null.
	 */
	public final void doAfter
		(double dt,
		 Event event)
		{
		sim.doAfter (dt, event);
		}

	/**
	 * Perform this event. Called by the {@linkplain Simulation} when the
	 * simulation time equals the time when this event is scheduled to take
	 * place.
	 */
	public abstract void perform();

	}
