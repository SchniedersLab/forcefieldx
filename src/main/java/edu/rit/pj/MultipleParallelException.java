//******************************************************************************
//
// File:    MultipleParallelException.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.MultipleParallelException
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

import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Map;

/**
 * Class MultipleParallelException is thrown to indicate that multiple threads
 * in a parallel team threw exceptions while executing a parallel construct. The
 * individual exceptions thrown by the threads are contained in a {@linkplain
 * java.util.Map Map} mapping the thread index (type Integer) to the exception
 * thrown by that thread (type Throwable).
 *
 * @author  Alan Kaminsky
 * @version 16-May-2007
 */
public class MultipleParallelException
	extends Exception
	{

// Hidden data members.

	private Map<Integer,Throwable> myMap;

// Exported constructors.

	/**
	 * Create a new multiple parallel exception with no detail message and no
	 * exception map.
	 */
	public MultipleParallelException()
		{
		super();
		}

	/**
	 * Create a new multiple parallel exception with the given detail message
	 * and no exception map.
	 *
	 * @param  theMessage  Detail message.
	 */
	public MultipleParallelException
		(String theMessage)
		{
		super (theMessage);
		}

	/**
	 * Create a new multiple parallel exception with no detail message and the
	 * given exception map.
	 *
	 * @param  theMap  Exception map.
	 */
	public MultipleParallelException
		(Map<Integer,Throwable> theMap)
		{
		super();
		myMap = theMap;
		}

	/**
	 * Create a new multiple parallel exception with the given detail message
	 * and the given exception map.
	 *
	 * @param  theMessage  Detail message.
	 * @param  theMap      Exception map.
	 */
	public MultipleParallelException
		(String theMessage,
		 Map<Integer,Throwable> theMap)
		{
		super (theMessage);
		myMap = theMap;
		}

// Exported operations.

	/**
	 * Obtain this multiple parallel exception's exception map.
	 *
	 * @return  Exception map, or null if none.
	 */
	public Map<Integer,Throwable> getExceptionMap()
		{
		return myMap;
		}

	/**
	 * Print this throwable and its backtrace to the specified print stream. The
	 * stack traces of this multiple parallel exception itself and of all the
	 * wrapped exceptions in the exception map are printed.
	 *
	 * @param  s  Print stream to use for output.
	 */
	public void printStackTrace
		(PrintStream s)
		{
		synchronized (s)
			{
			super.printStackTrace (s);
			if (myMap != null)
				{
				for
					(Map.Entry<Integer,Throwable> entry : myMap.entrySet())
					{
					s.println ("Parallel team thread " + entry.getKey() + ":");
					entry.getValue().printStackTrace (s);
					}
				}
			}
		}

	/**
	 * Print this throwable and its backtrace to the specified print writer. The
	 * stack traces of this multiple parallel exception itself and of all the
	 * wrapped exceptions in the exception map are printed.
	 *
	 * @param  s  Print writer to use for output.
	 */
	public void printStackTrace
		(PrintWriter s)
		{
		synchronized (s)
			{
			super.printStackTrace (s);
			if (myMap != null)
				{
				for
					(Map.Entry<Integer,Throwable> entry : myMap.entrySet())
					{
					s.println ("Parallel team thread " + entry.getKey() + ":");
					entry.getValue().printStackTrace (s);
					}
				}
			}
		}

	}
