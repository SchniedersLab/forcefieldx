//******************************************************************************
//
// File:    Status.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.Status
//
// This Java source file is copyright (C) 2005 by Alan Kaminsky. All rights
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

package edu.rit.mp;

/**
 * Class Status provides the result of receiving a message in the Message
 * Protocol (MP).
 *
 * @author  Alan Kaminsky
 * @version 29-Dec-2005
 */
public class Status
	{

// Exported data members.

	/**
	 * The channel from which the message was received.
	 */
	public Channel channel;

	/**
	 * The tag from the message that was received.
	 */
	public int tag;

	/**
	 * The actual number of items in the message that was received.
	 */
	public int length;

// Hidden constructors.

	/**
	 * Construct a new status object.
	 *
	 * @param  channel  Channel.
	 * @param  tag      Tag.
	 * @param  length   Length.
	 */
	Status
		(Channel channel,
		 int tag,
		 int length)
		{
		this.channel = channel;
		this.tag = tag;
		this.length = length;
		}

	}
