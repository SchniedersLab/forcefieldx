//******************************************************************************
//
// File:    FixedIntegerSchedule.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.FixedIntegerSchedule
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
package edu.rit.pj;

import edu.rit.util.Range;

/**
 * Class FixedIntegerSchedule provides a fixed schedule object. The loop index
 * is type <code>int</code>. The loop iterations are apportioned among the parallel
 * team threads once at the beginning of the parallel for loop, with each thread
 * getting a fixed number of iterations, the same number of iterations for each
 * thread (plus or minus one).
 *
 * @author Alan Kaminsky
 * @version 27-Jan-2010
 */
class FixedIntegerSchedule
        extends IntegerSchedule {

// Hidden data members.
    // Chunk for each thread.
    private Range[] myChunk;

// Exported constructors.
    /**
     * Construct a new fixed schedule object.
     */
    public FixedIntegerSchedule() {
        super();
    }

    /**
     * Construct a new fixed schedule object. This constructor is for use by the
     * <code>IntegerSchedule.parse()</code> method.
     *
     * @param args Array of argument strings.
     */
    public FixedIntegerSchedule(String[] args) {
        super();
        throw new IllegalArgumentException("FixedIntegerSchedule(): Usage: -Dpj.schedule=fixed");
    }

// Exported operations.
    /**
     * Determine if this schedule is a fixed schedule. For a parallel team with
     * <I>K</I> threads, a fixed schedule partitions the loop index range into
     * exactly <I>K</I> chunks, one chunk for each thread, each chunk with
     * predetermined upper and lower bounds.
     *
     * @return True if this is a fixed schedule, false otherwise.
     */
    public boolean isFixedSchedule() {
        return true;
    }

// Hidden operations.
    /**
     * {@inheritDoc}
     *
     * Start generating chunks of iterations for a parallel for loop using this
     * schedule.
     * <P>
     * The <code>start()</code> method is only called by a single thread in the
     * Parallel Java middleware.
     */
    public void start(int K,
            Range theLoopRange) {
        myChunk = theLoopRange.subranges(K);
        for (int i = 0; i < K; ++i) {
            if (myChunk[i].length() == 0) {
                myChunk[i] = null;
            }
        }
    }

    /**
     * {@inheritDoc}
     *
     * Obtain the next chunk of iterations for the given thread index. If there
     * are more iterations, a range object is returned whose lower bound, upper
     * bound, and stride specify the chunk of iterations to perform. The
     * returned range object's stride is the same as that given to the
     * <code>start()</code> method. The returned range object's lower bound and
     * upper bound are contained within the range given to the <code>start()</code>
     * method. If there are no more iterations, null is returned.
     * <P>
     * The <code>next()</code> method is called by multiple parallel team threads in
     * the Parallel Java middleware. The <code>next()</code> method must be multiple
     * thread safe.
     */
    public Range next(int theThreadIndex) {
        Range chunk = myChunk[theThreadIndex];
        myChunk[theThreadIndex] = null;
        return chunk;
    }

}
