//******************************************************************************
//
// File:    BackendClassLoader.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.BackendClassLoader
//
// This Java source file is copyright (C) 2012 by Alan Kaminsky. All rights
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
package edu.rit.pj.cluster;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.net.URL;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Class BackendClassLoader provides a class loader for a job backend process in
 * the PJ cluster middleware. If a backend class loader is requested to load a
 * class and the parent class loaders cannot do so, the backend class loader
 * sends a request for the class file to the job frontend process, waits for the
 * job frontend process to send the class file, and loads the class.
 *
 * @author Alan Kaminsky
 * @version 15-Jun-2012
 */
public class BackendClassLoader
        extends ClassLoader {

// Hidden data members.
    private JobBackendRef myJobBackend;
    private JobFrontendRef myJobFrontend;
    private ResourceCache myCache;

    // Map from resource name to resource URL for non-class-file resources.
    private Map<String, URL> myResourceURLMap
            = Collections.synchronizedMap(new HashMap<String, URL>());

// Exported constructors.
    /**
     * Construct a new backend class loader. The parent class loader is the one
     * returned by <TT>ClassLoader.getSystemClassLoader()</TT>. Class files will
     * be requested from <TT>theJobFrontend</TT>. Class files will be stored in
     * <TT>theCache</TT>.
     *
     * @param theJobBackend Reference to job backend.
     * @param theJobFrontend Reference to job frontend.
     * @param theCache Resource cache.
     */
    public BackendClassLoader(JobBackendRef theJobBackend,
            JobFrontendRef theJobFrontend,
            ResourceCache theCache) {
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
     * @param parent Parent class loader.
     * @param theJobBackend Reference to job backend.
     * @param theJobFrontend Reference to job frontend.
     * @param theCache Resource cache.
     */
    public BackendClassLoader(ClassLoader parent,
            JobBackendRef theJobBackend,
            JobFrontendRef theJobFrontend,
            ResourceCache theCache) {
        super(parent);
        myJobBackend = theJobBackend;
        myJobFrontend = theJobFrontend;
        myCache = theCache;
    }

// Hidden operations.
    /**
     * {@inheritDoc}
     *
     * Find the class with the given name.
     * @exception ClassNotFoundException Thrown if the class could not be found.
     */
    protected Class<?> findClass(String className)
            throws ClassNotFoundException {
        try {
            // Convert class name to resource name.
            String resourceName = className.replace('.', '/') + ".class";

            // If the resource is not in the cache, ask the Job Frontend for it.
            if (!myCache.contains(resourceName)) {
                myJobFrontend.requestResource(myJobBackend, resourceName);
            }

            // Wait until the resource shows up in the cache.
            byte[] content = myCache.get(resourceName);
            if (content == null) {
                throw new ClassNotFoundException("Class " + className + " not found");
            }

            // Load the class.
            return defineClass(className, content, 0, content.length);
        } catch (IOException exc) {
            throw new ClassNotFoundException("Class " + className + " not found due to I/O error",
                    exc);
        } catch (InterruptedException exc) {
            throw new ClassNotFoundException("Class " + className + " not found because thread interrupted",
                    exc);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Find the resource with the given name.
     */
    protected URL findResource(String name) {
        try {
            URL url = myResourceURLMap.get(name);

            if (url == null) {
                // If the resource is not in the cache, ask the Job Frontend for
                // it.
                if (!myCache.contains(name)) {
                    myJobFrontend.requestResource(myJobBackend, name);
                }

                // Wait until the resource shows up in the cache.
                byte[] content = myCache.get(name);
                if (content == null) {
                    return null;
                }

                // Store the resource contents in a temporary file, which will
                // be deleted when the JVM exits.
                String fname = new File(name).getName();
                int i = fname.lastIndexOf('.');
                String fprefix = i == -1 ? fname : fname.substring(0, i);
                String fsuffix = i == -1 ? null : fname.substring(i);
                File file = File.createTempFile(fprefix + "_tmp", fsuffix);
                OutputStream out
                        = new BufferedOutputStream(new FileOutputStream(file));
                out.write(content);
                out.close();
                file.deleteOnExit();

                // Map resource name to temporary file URL.
                url = file.toURI().toURL();
                myResourceURLMap.put(name, url);
            }

            return url;
        } catch (IOException exc) {
            return null;
        } catch (InterruptedException exc) {
            return null;
        }
    }

}
