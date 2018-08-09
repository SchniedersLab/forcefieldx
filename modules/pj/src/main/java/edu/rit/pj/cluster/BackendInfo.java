//******************************************************************************
//
// File:    BackendInfo.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.BackendInfo
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

/**
 * Class BackendInfo provides a record of information about one backend node in
 * the PJ cluster middleware.
 *
 * @author Alan Kaminsky
 * @version 20-Jun-2012
 */
public class BackendInfo {

// Exported enumerations.
    /**
     * The state of a backend node.
     */
    public static enum State {

        /**
         * The backend is available for jobs.
         */
        IDLE("Idle"),
        /**
         * The backend is reserved for a job that has not yet started running.
         */
        RESERVED("Reserved"),
        /**
         * The backend is running a job.
         */
        RUNNING("Running"),
        /**
         * The backend has failed.
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
     * The backend's name.
     */
    public String name;

    /**
     * The total number of CPUs in the backend.
     */
    public int totalCpus;

    /**
     * The backend's state.
     */
    public State state;

    /**
     * The time when the backend entered its current state (milliseconds since
     * midnight 01-Jan-1970 GMT).
     */
    public long stateTime;

    /**
     * The host name for SSH remote logins to the backend.
     */
    public String host;

    /**
     * The full pathname for executing the Java Virtual Machine (JVM) on the
     * backend.
     */
    public String jvm;

    /**
     * The Java class path for the Parallel Java Library on the backend.
     */
    public String classpath;

    /**
     * Array of command line flags for the JVM (zero or more).
     */
    public String[] jvmflags;

    /**
     * Shell command string on the backend.
     */
    public String shellCommand;

    /**
     * The job that has reserved or is running on the backend.
     */
    public JobInfo job;

// Exported constructors.
    /**
     * Construct a new backend information record.
     *
     * @param name The backend's name.
     * @param totalCpus The total number of CPUs in the backend.
     * @param state The backend's state.
     * @param stateTime The time when the backend entered its current state.
     * @param stateTime The time when the backend entered its current state.
     * @param host The host name for SSH remote logins to the backend.
     * @param jvm The full pathname for executing the Java Virtual Machine (JVM)
     * on the backend.
     * @param jvmflags Array of command line flags for the JVM (zero or more).
     * @param classpath The Java class path for the Parallel Java Library on the
     * backend.
     * @param jvmflags Array of command line flags for the JVM (zero or more).
     * @param shellCommand Shell command string.
     */
    public BackendInfo(String name,
            int totalCpus,
            State state,
            long stateTime,
            String host,
            String jvm,
            String classpath,
            String[] jvmflags,
            String shellCommand) {
        this.name = name;
        this.totalCpus = totalCpus;
        this.state = state;
        this.stateTime = stateTime;
        this.host = host;
        this.jvm = jvm;
        this.classpath = classpath;
        this.jvmflags = jvmflags;
        this.shellCommand = shellCommand;
    }

}
