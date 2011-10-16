//******************************************************************************
//
// File:    ConnectListener.java
// Package: edu.rit.mp
// Unit:    Interface edu.rit.mp.ConnectListener
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

package edu.rit.mp;

import java.io.IOException;

/**
 * Interface ConnectListener specifies the interface for an object that is
 * notified whenever a {@linkplain Channel} is connected in a {@linkplain
 * ChannelGroup}.
 *
 * @author  Alan Kaminsky
 * @version 20-Dec-2007
 */
public interface ConnectListener
	{

// Exported operations.

	/**
	 * Report that a channel was connected in the given channel group, initiated
	 * by the near end. The channel group calls the <TT>nearEndConnected()</TT>
	 * method of a registered connect listener when the channel group's
	 * <TT>connect()</TT> method is called by code in the same process.
	 *
	 * @param  theChannelGroup  Channel group that is calling this method.
	 * @param  theChannel       Newly created channel.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void nearEndConnected
		(ChannelGroup theChannelGroup,
		 Channel theChannel)
		throws IOException;

	/**
	 * Report that a channel was connected in the given channel group, initiated
	 * by the far end. The channel group calls the <TT>farEndConnected()</TT>
	 * method of a registered connect listener when an incoming connection
	 * request is received from another process.
	 *
	 * @param  theChannelGroup  Channel group that is calling this method.
	 * @param  theChannel       Newly created channel.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void farEndConnected
		(ChannelGroup theChannelGroup,
		 Channel theChannel)
		throws IOException;

	}
