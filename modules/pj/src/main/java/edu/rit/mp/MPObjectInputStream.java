//******************************************************************************
//
// File:    MPObjectInputStream.java
// Package: edu.rit.mp
// Unit:    Class edu.rit.mp.MPObjectInputStream
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
package edu.rit.mp;

import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectStreamClass;

/**
 * Class MPObjectInputStream provides an object input stream that can load
 * classes from an alternate class loader.
 *
 * @author Alan Kaminsky
 * @version 31-Jan-2006
 */
class MPObjectInputStream
        extends ObjectInputStream {

// Hidden data members.
    private ClassLoader myClassLoader;

// Exported constructors.
    /**
     * Create a new MP object input stream. An alternate class loader will not
     * be used.
     *
     * @param in Underlying input stream.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public MPObjectInputStream(InputStream in)
            throws IOException {
        super(in);
    }

    /**
     * Create a new MP object input stream. The given class loader will be used
     * to load classes; if null, an alternate class loader will not be used.
     *
     * @param in Underlying input stream.
     * @param cl Alternate class loader, or null.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public MPObjectInputStream(InputStream in,
            ClassLoader cl)
            throws IOException {
        super(in);
        myClassLoader = cl;
    }

// Hidden operations.
    /**
     * {@inheritDoc}
     *
     * Load the local class equivalent of the specified stream class
     * description.
     * @exception IOException Thrown if an I/O error occurred.
     * @exception ClassNotFoundException Thrown if the local class could not be
     * found.
     */
    protected Class<?> resolveClass(ObjectStreamClass desc)
            throws IOException, ClassNotFoundException {
        try {
            return super.resolveClass(desc);
        } catch (ClassNotFoundException exc) {
            if (myClassLoader != null) {
                return Class.forName(desc.getName(), false, myClassLoader);
            } else {
                throw exc;
            }
        }
    }

}
