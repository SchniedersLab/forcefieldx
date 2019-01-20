//******************************************************************************
//
// File:    JobGenerator.java
// Package: edu.rit.pj.job
// Unit:    Class edu.rit.pj.job.JobGenerator
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
package edu.rit.pj.job;

import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * Class JobGenerator is the abstract base class for an object that generates a
 * group of {@linkplain Job}s.
 * <P>
 * Jobs are numbered from 0 to <I>N</I>&minus;1, where <I>N</I> is the number of
 * jobs in the group. A subclass must override the <code>jobCount()</code> method to
 * return <I>N</I>. A subclass must override the <code>createJob()</code> method to
 * create and return the job corresponding to a given job number. The job
 * generator need not create all the jobs in the group, and it need not create
 * them in any particular order.
 * <P>
 * Class JobGenerator provides the <code>omit()</code> method to omit generating
 * certain job numbers. This is used for checkpointing. For further information,
 * see class {@linkplain Runner}.
 *
 * @author Alan Kaminsky
 * @version 08-Oct-2010
 */
public abstract class JobGenerator
        implements Iterable<Job> {

// Hidden data members.
    private Set<Integer> myOmittedJobNumbers;

// Exported constructors.
    /**
     * Construct a new job generator.
     */
    public JobGenerator() {
        myOmittedJobNumbers = new HashSet<Integer>();
    }

// Exported operations.
    /**
     * Omit the job numbers in the given set when generating jobs. To be
     * effective, <code>omit()</code> must be called before calling
     * <code>iterator()</code>. A snapshot of the given set is taken; changing the
     * set's contents thereafter will not affect the job numbers to be omitted.
     *
     * @param theOmittedJobNumbers Set of job numbers to be omitted.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>theOmittedJobNumbers</code> is null.
     */
    public void omit(Set<Integer> theOmittedJobNumbers) {
        myOmittedJobNumbers.clear();
        myOmittedJobNumbers.addAll(theOmittedJobNumbers);
    }

    /**
     * Get an iterator for generating the jobs in the job group.
     *
     * @return Iterator.
     */
    public Iterator<Job> iterator() {
        return new Iterator<Job>() {
            private int N = jobCount();
            private int myJobNumber = -1;
            private boolean generated = true;

            public boolean hasNext() {
                advance();
                return myJobNumber < N;
            }

            public Job next() {
                advance();
                if (myJobNumber >= N) {
                    throw new NoSuchElementException();
                }
                generated = true;
                return createJob(myJobNumber);
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }

            private void advance() {
                if (generated) {
                    generated = false;
                    do {
                        ++myJobNumber;
                    } while (myJobNumber < N
                            && myOmittedJobNumbers.contains(myJobNumber));
                }
            }
        };
    }

// Hidden operations.
    /**
     * Returns the number of jobs in the job group, <I>N</I>.
     *
     * @return Number of jobs.
     */
    protected abstract int jobCount();

    /**
     * Create the job with the given job number. This method must create and
     * return an instance of class {@linkplain Job} whose job number is
     * <code>theJobNumber</code>.
     *
     * @param theJobNumber Job number (0 .. <I>N</I>&minus;1).
     * @return a {@link edu.rit.pj.job.Job} object.
     */
    protected abstract Job createJob(int theJobNumber);

}
