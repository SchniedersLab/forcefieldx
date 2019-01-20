//******************************************************************************
//
// File:    PrintStreamLogger.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.PrintStreamLogger
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
package edu.rit.util;

import java.io.PrintStream;
import java.util.Date;

/**
 * Class PrintStreamLogger provides an object that logs messages to a print
 * stream.
 *
 * @author Alan Kaminsky
 * @version 16-Apr-2008
 */
public class PrintStreamLogger
        implements Logger {

// Hidden operations.
    private PrintStream out;

// Exported constructors.
    /**
     * Construct a new print stream logger that logs to <code>System.err</code>.
     */
    public PrintStreamLogger() {
        this.out = System.err;
    }

    /**
     * Construct a new print stream logger that logs to the given print stream.
     *
     * @param out Print stream.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>out</code> is null.
     */
    public PrintStreamLogger(PrintStream out) {
        if (out == null) {
            throw new NullPointerException("PrintStreamLogger(): Print stream is null");
        }
        this.out = out;
    }

// Exported operations.
    /**
     * Log the given message.
     *
     * @param msg Message.
     */
    public void log(String msg) {
        log(System.currentTimeMillis(), msg, null);
    }

    /**
     * {@inheritDoc}
     *
     * Log the given exception.
     */
    public void log(Throwable exc) {
        log(System.currentTimeMillis(), null, exc);
    }

    /**
     * {@inheritDoc}
     *
     * Log the given message and exception.
     */
    public void log(String msg,
            Throwable exc) {
        log(System.currentTimeMillis(), msg, exc);
    }

    /**
     * Log the given date and message.
     *
     * @param date Date and time in milliseconds since midnight 01-Jan-1970 UTC.
     * @param msg Message.
     */
    public void log(long date,
            String msg) {
        log(date, msg, null);
    }

    /**
     * {@inheritDoc}
     *
     * Log the given date and exception.
     */
    public void log(long date,
            Throwable exc) {
        log(date, null, exc);
    }

    /**
     * {@inheritDoc}
     *
     * Log the given date, message, and exception.
     */
    public void log(long date,
            String msg,
            Throwable exc) {
        synchronized (out) {
            out.print(new Date(date));
            if (msg != null) {
                out.print(' ');
                out.print(msg);
            }
            out.println();
            if (exc != null) {
                exc.printStackTrace(out);
            }
        }
    }

}
