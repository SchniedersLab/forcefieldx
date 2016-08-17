//******************************************************************************
//
// File:    Instance.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Instance
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
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
package edu.rit.util;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

/**
 * Class Instance provides static methods for creating instances of classes.
 *
 * @author Alan Kaminsky
 * @version 09-Oct-2010
 */
public class Instance {

// Prevent construction.
    private Instance() {
    }

// Exported operations.
    /**
     * Create a new instance of a class as specified by the given string. The
     * string must consist of a fully-qualified class name, a left parenthesis,
     * zero or more comma-separated arguments, and a right parenthesis. No
     * whitespace is allowed. This method attempts to find a constructor for the
     * given class as follows, where <I>N</I> is the number of arguments:
     * <UL>
     * <P>
     * <LI>
     * If <I>N</I> = 0, use a no-argument constructor.
     * <P>
     * <LI>
     * Else if all arguments are integers, use a constructor with <I>N</I>
     * arguments of type <TT>int</TT>.
     * <P>
     * <LI>
     * Else if all arguments are integers and there is no such constructor, use
     * a constructor with one argument of type <TT>int[]</TT>.
     * <P>
     * <LI>
     * Else if not all arguments are integers, use a constructor with <I>N</I>
     * arguments of type <TT>String</TT>.
     * <P>
     * <LI>
     * Else if not all arguments are integers and there is no such constructor,
     * use a constructor with one argument of type <TT>String[]</TT>.
     * <P>
     * <LI>
     * Else throw a NoSuchMethodException.
     * </UL>
     * <P>
     * This method invokes the chosen constructor, passing in the given argument
     * values, and returns a reference to the newly-created instance.
     * <P>
     * <I>Note:</I> To find the given class, the calling thread's context class
     * loader is used.
     *
     * @param s Constructor expression string.
     * @return New instance.
     * @exception IllegalArgumentException (unchecked exception) Thrown if
     * <TT>s</TT> does not obey the required syntax.
     * @exception ClassNotFoundException Thrown if the given class cannot be
     * found.
     * @exception NoSuchMethodException Thrown if a suitable constructor cannot
     * be found in the given class.
     * @exception InstantiationException Thrown if an instance cannot be created
     * because the given class is an interface or an abstract class.
     * @exception IllegalAccessException Thrown if an instance cannot be created
     * because the calling method does not have access to the given constructor.
     * @exception InvocationTargetException Thrown if the given constructor
     * throws an exception.
     * @throws java.lang.ClassNotFoundException if any.
     * @throws java.lang.NoSuchMethodException if any.
     * @throws java.lang.InstantiationException if any.
     * @throws java.lang.IllegalAccessException if any.
     * @throws java.lang.reflect.InvocationTargetException if any.
     */
    public static Object newInstance(String s)
            throws
            ClassNotFoundException,
            NoSuchMethodException,
            InstantiationException,
            IllegalAccessException,
            InvocationTargetException {
        int i, j;

        // Parse class name.
        i = 0;
        j = s.indexOf('(', i);
        if (j == -1) {
            throw new IllegalArgumentException("Instance.newInstance(\"" + s + "\"): Missing '('");
        }
        if (i == j) {
            throw new IllegalArgumentException("Instance.newInstance(\"" + s + "\"): Missing class name");
        }
        String classname = s.substring(i, j);

        // Parse arguments.
        i = j + 1;
        j = s.indexOf(')', i);
        if (j == -1) {
            throw new IllegalArgumentException("Instance.newInstance(\"" + s + "\"): Missing ')'");
        }
        if (j != s.length() - 1) {
            throw new IllegalArgumentException("Instance.newInstance(\"" + s
                    + "\"): Extraneous characters after ')'");
        }
        String arguments = s.substring(i, j);

        // Parse individual arguments.
        String[] args;
        if (i == j) {
            args = new String[0];
        } else {
            args = arguments.split(",", -1);
        }
        Integer[] intargs = new Integer[args.length];
        boolean allAreInts = true;
        for (i = 0; i < args.length; ++i) {
            try {
                intargs[i] = new Integer(args[i]);
            } catch (NumberFormatException exc) {
                allAreInts = false;
            }
        }

        // Get class.
        Class<?> theClass = Class.forName(classname,
                true,
                Thread.currentThread().getContextClassLoader());

        // Get constructor and create instance.
        Constructor<?> ctor = null;
        Class<?>[] argtypes = null;

        // No-argument constructor.
        if (args.length == 0) {
            try {
                ctor = theClass.getConstructor();
                return ctor.newInstance();
            } catch (NoSuchMethodException exc) {
            }
        }

        // Constructor(int,int,...,int).
        if (allAreInts) {
            try {
                argtypes = new Class<?>[args.length];
                for (i = 0; i < args.length; ++i) {
                    argtypes[i] = Integer.TYPE;
                }
                ctor = theClass.getConstructor(argtypes);
                return ctor.newInstance((Object[]) intargs);
            } catch (NoSuchMethodException exc) {
            }
        }

        // Constructor(int[]).
        if (allAreInts) {
            try {
                ctor = theClass.getConstructor(int[].class);
                return ctor.newInstance((Object) intargs);
            } catch (NoSuchMethodException exc) {
            }
        }

        // Constructor(String,String,...,String).
        try {
            argtypes = new Class<?>[args.length];
            for (i = 0; i < args.length; ++i) {
                argtypes[i] = String.class;
            }
            ctor = theClass.getConstructor(argtypes);
            return ctor.newInstance((Object[]) args);
        } catch (NoSuchMethodException exc) {
        }

        // Constructor(String[]).
        try {
            ctor = theClass.getConstructor(String[].class);
            return ctor.newInstance((Object) args);
        } catch (NoSuchMethodException exc) {
        }

        // Could not find suitable constructor.
        throw new NoSuchMethodException("Instance.newInstance(\"" + s
                + "\"): Cannot find suitable constructor");
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 */
//	public static void main
//		(String[] args)
//		throws Exception
//		{
//		System.out.println (Instance.newInstance (args[0]));
//		}
}
