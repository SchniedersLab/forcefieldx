//******************************************************************************
//
// File:    Message.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.Message
//
// This Java source file is copyright (C) 2006 by Alan Kaminsky. All rights
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

package edu.rit.pj.cluster;

import java.io.IOException;
import java.io.Externalizable;
import java.io.ObjectInput;
import java.io.ObjectOutput;

/**
 * Class Message is the abstract base class for a message sent to a process in
 * the PJ cluster middleware.
 *
 * @author  Alan Kaminsky
 * @version 20-Nov-2006
 */
public abstract class Message
	implements Externalizable
	{

// Exported constants.

	/**
	 * The message tag for a message from a job backend process.
	 */
	public static final int FROM_JOB_BACKEND = 1;

	/**
	 * The message tag for a message from a job frontend process.
	 */
	public static final int FROM_JOB_FRONTEND = 2;

	/**
	 * The message tag for a message from a job launcher process.
	 */
	public static final int FROM_JOB_LAUNCHER = 3;

	/**
	 * The message tag for a message from a job scheduler process.
	 */
	public static final int FROM_JOB_SCHEDULER = 4;

	/**
	 * The message tag for a message containing data to write to a file.
	 */
	public static final int FILE_WRITE_DATA = 5;

	/**
	 * The message tag for a message containing data read from a file.
	 */
	public static final int FILE_READ_DATA = 6;

// Hidden data members.

	private static final long serialVersionUID = -3891573184096499571L;

	private int myTag;

// Exported constructors.

	/**
	 * Construct a new message.
	 */
	public Message()
		{
		}

	/**
	 * Construct a new message with the given message tag.
	 *
	 * @param  theTag  Message tag to use when sending this message.
	 */
	public Message
		(int theTag)
		{
		myTag = theTag;
		}

// Exported operations.

	/**
	 * Get the message tag to use when sending this message.
	 *
	 * @return  Message tag.
	 */
	public int getTag()
		{
		return myTag;
		}

	}
