//******************************************************************************
//
// File:    ChannelClosedException.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.ChannelClosedException
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
package edu.rit.mp;

import java.io.IOException;

/**
 * Class ChannelClosedException is thrown to indicate that an I/O operation
 * failed because the channel was closed.
 *
 * @author Alan Kaminsky
 * @version 26-Oct-2006
 */
public class ChannelClosedException
        extends IOException {

// Exported constructors.
    /**
     * Create a new channel closed exception with no detail message and no
     * cause.
     */
    public ChannelClosedException() {
        super();
    }

    /**
     * Create a new channel closed exception with the given detail message and
     * no cause.
     *
     * @param msg Detail message.
     */
    public ChannelClosedException(String msg) {
        super(msg);
    }

    /**
     * Create a new channel closed exception with no detail message and the
     * given cause.
     *
     * @param cause Cause.
     */
    public ChannelClosedException(Throwable cause) {
        super();
        initCause(cause);
    }

    /**
     * Create a new channel closed exception with the given detail message and
     * the given cause.
     *
     * @param msg Detail message.
     * @param cause Cause.
     */
    public ChannelClosedException(String msg,
            Throwable cause) {
        super(msg);
        initCause(cause);
    }

}
