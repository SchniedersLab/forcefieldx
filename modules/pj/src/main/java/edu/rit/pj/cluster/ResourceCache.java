//******************************************************************************
//
// File:    ResourceCache.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.ResourceCache
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
package edu.rit.pj.cluster;

import java.util.HashMap;
import java.util.Map;

/**
 * Class ResourceCache provides a cache of resources, indexed by resource name.
 * A resource is a piece of content (sequence of bytes) obtained from a class
 * loader.
 * <P>
 * <I>Note:</I> Class ResourceCache is multiple thread safe.
 *
 * @author Alan Kaminsky
 * @version 26-Oct-2006
 */
public class ResourceCache {

// Hidden data members.
    private Map<String, ResourceInfo> myMap
            = new HashMap<String, ResourceInfo>();

// Hidden helper classes.
    private static class ResourceInfo {

        public byte[] content;

        public ResourceInfo(byte[] content) {
            this.content = content;
        }
    }

// Exported constructors.
    /**
     * Construct a new resource cache.
     */
    public ResourceCache() {
    }

// Exported operations.
    /**
     * Determine if this resource cache contains resource information for the
     * given resource name. If the answer is yes, the resource content may or
     * may not have been found.
     *
     * @param name Resource name.
     * @return True if this resource cache contains resource information for
     * <TT>name</TT>, false otherwise.
     */
    public synchronized boolean contains(String name) {
        ResourceInfo info = myMap.get(name);
        return info != null;
    }

    /**
     * Determine if this resource cache contains the resource content for the
     * given resource name. If the answer is yes, the resource content was
     * found, otherwise the resource content was not found or no resource
     * information is available.
     *
     * @param name Resource name.
     * @return True if this resource cache contains the resource content for
     * <TT>name</TT>, false otherwise.
     */
    public synchronized boolean containsContent(String name) {
        ResourceInfo info = myMap.get(name);
        return info == null ? false : info.content != null;
    }

    /**
     * Obtain the resource content for the given resource name from this
     * resource cache (blocking). This method will block if necessary until this
     * resource cache contains the content for <TT>name</TT> or until this
     * resource cache knows the content was not found.
     *
     * @param name Resource name.
     * @return Resource content, or null if not found.
     * @exception InterruptedException Thrown if the calling thread is
     * interrupted while blocked in this method.
     * @throws java.lang.InterruptedException if any.
     */
    public synchronized byte[] get(String name)
            throws InterruptedException {
        ResourceInfo info = myMap.get(name);
        while (info == null) {
            wait();
            info = myMap.get(name);
        }
        return info.content;
    }

    /**
     * Obtain the resource content for the given resource name from this
     * resource cache (non-blocking). This method will return null if the
     * resource content was not found or no resource information is available.
     *
     * @param name Resource name.
     * @return Resource content, or null if not found or no information is
     * available.
     */
    public synchronized byte[] getNoWait(String name) {
        ResourceInfo info = myMap.get(name);
        return info == null ? null : info.content;
    }

    /**
     * Store the given resource content under the given resource name in this
     * resource cache. Any existing content for <TT>name</TT> is overwritten.
     * <P>
     * <I>Note:</I> The resource cache assumes that the contents of
     * <TT>content</TT> are not changed after <TT>put()</TT> is called.
     *
     * @param name Resource name.
     * @param content Resource content, or null if not found.
     */
    public synchronized void put(String name,
            byte[] content) {
        myMap.put(name, new ResourceInfo(content));
        notifyAll();
    }

    /**
     * Remove the resource content for the given resource name from this
     * resource cache. If there is no content for <TT>name</TT>, the
     * <TT>remove()</TT> method does nothing.
     *
     * @param name Resource name.
     */
    public synchronized void remove(String name) {
        myMap.remove(name);
    }

}
