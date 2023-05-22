//******************************************************************************
//
// File:    InvalidMatrixFileException.java
// Package: edu.rit.io
// Unit:    Class edu.rit.io.InvalidMatrixFileException
//
// This Java source file is copyright (C) 2008 by Alan Kaminsky. All rights
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
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.io;

import java.io.IOException;
import java.io.Serial;

/**
 * Class InvalidMatrixFileException provides an exception thrown when the
 * contents of a matrix file are invalid. The detail message and/or chained
 * exception give further information about the problem.
 *
 * @author Alan Kaminsky
 * @version 07-Jan-2008
 */
public class InvalidMatrixFileException
        extends IOException {

    @Serial
    private static final long serialVersionUID = 1L;

// Exported constructors.
    /**
     * Construct a new invalid matrix file exception with no detail message and
     * no cause.
     */
    public InvalidMatrixFileException() {
        super();
    }

    /**
     * Construct a new invalid matrix file exception with the given detail
     * message and no cause.
     *
     * @param message Detail message.
     */
    public InvalidMatrixFileException(String message) {
        super(message);
    }

    /**
     * Construct a new invalid matrix file exception with the given cause and
     * the default detail message.
     *
     * @param cause Cause.
     */
    public InvalidMatrixFileException(Throwable cause) {
        super();
        initCause(cause);
    }

    /**
     * Construct a new invalid matrix file exception with the given detail
     * message and the given cause.
     *
     * @param message Detail message.
     * @param cause Cause.
     */
    public InvalidMatrixFileException(String message,
            Throwable cause) {
        super(message);
        initCause(cause);
    }

}
