//******************************************************************************
//
// File:    ProcessInfo.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.ProcessInfo
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
package edu.rit.pj.cluster;

import java.net.InetSocketAddress;

import edu.rit.util.Timer;

/**
 * Class ProcessInfo provides a record of information about one job backend
 * process in the PJ cluster middleware.
 *
 * @author Alan Kaminsky
 * @version 21-May-2008
 */
public class ProcessInfo {

// Exported enumerations.
    /**
     * The state of a job backend process.
     */
    public static enum State {

        /**
         * The job backend process has not started yet.
         */
        NOT_STARTED("Not started"),
        /**
         * The job backend process is running.
         */
        RUNNING("Running"),
        /**
         * The job backend process has finished.
         */
        FINISHED("Finished"),
        /**
         * The job backend process has failed.
         */
        FAILED("Failed");

        private final String stringForm;

        /**
         * Construct a new State value.
         *
         * @param stringForm String form.
         */
        State(String stringForm) {
            this.stringForm = stringForm;
        }

        /**
         * Returns a string version of this State value.
         *
         * @return String version.
         */
        public String toString() {
            return stringForm;
        }
    }

// Exported data members.
    /**
     * The job backend process's state.
     */
    public State state;

    /**
     * The job backend node's name.
     */
    public String name;

    /**
     * The job backend process's rank.
     */
    public int rank;

    /**
     * Reference to the job backend process.
     */
    public JobBackendRef backend;

    /**
     * Host/port to which the job backend process is listening for middleware
     * messages.
     */
    public InetSocketAddress middlewareAddress;

    /**
     * Host/port to which the job backend process is listening for the world
     * communicator.
     */
    public InetSocketAddress worldAddress;

    /**
     * Host/port to which the job backend process is listening for the frontend
     * communicator, or null if the frontend communicator does not exist.
     */
    public InetSocketAddress frontendAddress;

    /**
     * Lease renewal timer.
     */
    public Timer renewTimer;

    /**
     * Lease expiration timer.
     */
    public Timer expireTimer;

    /**
     * Number of CPUs assigned to the job backend process.
     */
    public int Nt;

// Exported constructors.
    /**
     * Construct a new job information record.
     *
     * @param state The job backend process's state.
     * @param name The job backend processor's name.
     * @param rank The job backend process's rank.
     * @param backend Reference to the job backend process.
     * @param middlewareAddress Host/port to which the job backend process is
     * listening for middleware messages.
     * @param worldAddress Host/port to which the job backend process is
     * listening for the world communicator.
     * @param frontendAddress Host/port to which the job backend process is
     * listening for the frontend communicator, or null if the frontend
     * communicator does not exist.
     * @param renewTimer Lease renewal timer.
     * @param expireTimer Lease expiration timer.
     * @param Nt Number of CPUs assigned to the job backend process.
     */
    public ProcessInfo(State state,
            String name,
            int rank,
            JobBackendRef backend,
            InetSocketAddress middlewareAddress,
            InetSocketAddress worldAddress,
            InetSocketAddress frontendAddress,
            Timer renewTimer,
            Timer expireTimer,
            int Nt) {
        this.state = state;
        this.name = name;
        this.rank = rank;
        this.backend = backend;
        this.middlewareAddress = middlewareAddress;
        this.worldAddress = worldAddress;
        this.frontendAddress = frontendAddress;
        this.renewTimer = renewTimer;
        this.expireTimer = expireTimer;
        this.Nt = Nt;
    }

}
