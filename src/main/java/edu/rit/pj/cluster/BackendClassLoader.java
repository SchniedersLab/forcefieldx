//******************************************************************************
//
// File:    BackendClassLoader.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.BackendClassLoader
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

/**
 * Class BackendClassLoader provides a class loader for a job backend process in
 * the PJ cluster middleware. If a backend class loader is requested to load a
 * class and the parent class loaders cannot do so, the backend class loader
 * sends a request for the class file to the job frontend process, waits for the
 * job frontend process to send the class file, and loads the class.
 * <P>
 * Class BackendClassLoader supports loading of resources that are class files.
 * Support for loading of resources that are not class files is TBD.
 *
 * @author  Alan Kaminsky
 * @version 26-Oct-2006
 */
public class BackendClassLoader
	extends ClassLoader
	{

// Hidden data members.

	private JobBackendRef myJobBackend;
	private JobFrontendRef myJobFrontend;
	private ResourceCache myCache;

// Exported constructors.

	/**
	 * Construct a new backend class loader. The parent class loader is the one
	 * returned by <TT>ClassLoader.getSystemClassLoader()</TT>. Class files will
	 * be requested from <TT>theJobFrontend</TT>. Class files will be stored in
	 * <TT>theCache</TT>.
	 *
	 * @param  theJobBackend   Reference to job backend.
	 * @param  theJobFrontend  Reference to job frontend.
	 * @param  theCache        Resource cache.
	 */
	public BackendClassLoader
		(JobBackendRef theJobBackend,
		 JobFrontendRef theJobFrontend,
		 ResourceCache theCache)
		{
		super();
		myJobBackend = theJobBackend;
		myJobFrontend = theJobFrontend;
		myCache = theCache;
		}

	/**
	 * Construct a new backend class loader. The parent class loader is
	 * <TT>parent</TT>. Class files will be requested from
	 * <TT>theJobFrontend</TT>. Class files will be stored in <TT>theCache</TT>.
	 *
	 * @param  parent          Parent class loader.
	 * @param  theJobBackend   Reference to job backend.
	 * @param  theJobFrontend  Reference to job frontend.
	 * @param  theCache        Resource cache.
	 */
	public BackendClassLoader
		(ClassLoader parent,
		 JobBackendRef theJobBackend,
		 JobFrontendRef theJobFrontend,
		 ResourceCache theCache)
		{
		super (parent);
		myJobBackend = theJobBackend;
		myJobFrontend = theJobFrontend;
		myCache = theCache;
		}

// Hidden operations.

	/**
	 * Find the class with the given name.
	 *
	 * @param  className  Fully-qualified class name.
	 *
	 * @return  Class object.
	 *
	 * @exception  ClassNotFoundException
	 *     Thrown if the class could not be found.
	 */
	protected Class<?> findClass
		(String className)
		throws ClassNotFoundException
		{
		try
			{
			// Convert class name to resource name.
			String resourceName = className.replace ('.', '/') + ".class";

			// If the resource is not in the cache, ask the Job Frontend for it.
			if (! myCache.contains (resourceName))
				{
				myJobFrontend.requestResource (myJobBackend, resourceName);
				}

			// Wait until the resource shows up in the cache.
			byte[] content = myCache.get (resourceName);
			if (content == null)
				{
				throw new ClassNotFoundException
					("Class " + className + " not found");
				}

			// Load the class.
			return defineClass (className, content, 0, content.length);
			}

		catch (IOException exc)
			{
			throw new ClassNotFoundException
				("Class " + className + " not found due to I/O error",
				 exc);
			}

		catch (InterruptedException exc)
			{
			throw new ClassNotFoundException
				("Class " + className + " not found because thread interrupted",
				 exc);
			}
		}

	}
