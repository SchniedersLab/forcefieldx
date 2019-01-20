//******************************************************************************
//
// File:    Runner.java
// Package: edu.rit.pj.job
// Unit:    Class edu.rit.pj.job.Runner
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

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Date;
import java.util.HashSet;
import java.util.Scanner;

import edu.rit.mp.BooleanBuf;
import edu.rit.mp.buf.BooleanItemBuf;
import edu.rit.pj.Comm;
import edu.rit.pj.WorkerIteration;
import edu.rit.pj.WorkerRegion;
import edu.rit.pj.WorkerTeam;
import edu.rit.util.Instance;

/**
 * Class Runner is a parallel program that runs, in parallel, a group of
 * {@linkplain Job}s created by a {@linkplain JobGenerator}. The job generator
 * is specified on the command line as a constructor expression. An instance of
 * the class specified in the constructor expression is constructed, with the
 * constructor arguments specified in the constructor expression. For further
 * information, see class {@linkplain edu.rit.util.Instance Instance}.
 * <P>
 * The Runner program is targeted at three use cases:
 * <UL>
 *
 * <LI>
 * <B>Sequential jobs on a cluster parallel computer.</B> Each job is a
 * sequential (single-threaded) program. The Runner program is running on a
 * cluster parallel computer with <I>N</I> nodes and one CPU per node. Run the
 * Runner program as follows:
 * <PRE>
 *     java -Dpj.nn=<I>N</I> edu.rit.pj.job.Runner . . .
 * </PRE> The Runner program runs with one process per node and one thread per
 * process.
 *
 * <LI>
 * <B>Sequential jobs on a hybrid parallel computer.</B> Each job is a
 * sequential (single-threaded) program. The Runner program is running on a
 * hybrid SMP cluster parallel computer with <I>N</I> nodes and <I>C</I> total
 * CPUs. (For example, on a hybrid parallel computer with 10 nodes and 4 CPUs
 * per node, <I>C</I> = 40.) Run the Runner program as follows:
 * <PRE>
 *     java -Dpj.nn=<I>N</I> -Dpj.np=<I>C</I> edu.rit.pj.job.Runner . . .
 * </PRE> The Runner program runs with multiple processes per node and one
 * thread per process.
 *
 * <LI>
 * <B>SMP parallel jobs on a hybrid parallel computer.</B> Each job is an SMP
 * parallel (multi-threaded) program. The Runner program is running on a hybrid
 * SMP cluster parallel computer with <I>N</I> nodes and multiple CPUs per node.
 * Run the Runner program as follows:
 * <PRE>
 *     java -Dpj.nn=<I>N</I> edu.rit.pj.job.Runner . . .
 * </PRE> The Runner program runs with one process per node and multiple threads
 * per process, typically as many threads as there are CPUs on the node.
 * </UL>
 * <P>
 * All these processes form a <I>worker team.</I> The Runner program uses the
 * job generator specified on the command line to create jobs and sends each job
 * to a worker team process to be executed.
 * <P>
 * When the Runner program starts, it prints the job generator constructor
 * expression on the standard output. Whenever a job starts or finishes, the
 * Runner program prints a log message on the standard output consisting of the
 * job's number and description.
 * <P>
 * <B>Checkpointing.</B> It is recommended to redirect the Runner program's
 * standard output into a <I>checkpoint file.</I> If a failure occurs before the
 * Runner program finishes running all the jobs, the checkpoint file contains a
 * record of the job generator that was used as well as which jobs did and did
 * not finish. To resume the Runner program where it left off, specify the
 * checkpoint file name on the command line instead of a job generator
 * constructor expression. The Runner program reads the checkpoint file to
 * determine the job generator and the jobs that finished. The Runner program
 * then generates and runs the jobs that did not finish.
 * <P>
 * Usage: java edu.rit.pj.job.Runner { <I>generator</I> | <I>file</I> }
 * <BR><I>generator</I> = Job generator constructor expression
 * <BR><I>file</I> = Checkpoint file name
 *
 * @author Alan Kaminsky
 * @version 22-Oct-2010
 */
public class Runner {

// Prevent construction.
    private Runner() {
    }

// Global variables.
    private static PrintStream stdout = System.out;
    private static PrintStream stderr = System.err;

    private static Comm world;
    private static int rank;

    private static WorkerTeam team;

    private static String generatorExpression;
    private static HashSet<Integer> omitted;
    private static JobGenerator generator;

// Main program.
    /**
     * Main program.
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.lang.Exception if any.
     */
    public static void main(String[] args)
            throws Exception {
        if (args.length != 1) {
            usage();
        }

        // Initialize world communicator.
        Comm.init(args);
        world = Comm.world();
        rank = world.rank();

        // Set up worker team.
        team = new WorkerTeam();

        // Master process sets up job generator.
        if (rank == team.masterRank()) {
            omitted = new HashSet<Integer>();

            // Assume argument is a checkpoint file name and try to read it.
            Scanner scanner = null;
            try {
                scanner = new Scanner(new File(args[0]));
            } catch (IOException exc) {
            }

            // Read checkpoint file.
            if (scanner != null) {
                while (scanner.hasNextLine()) {
                    Scanner linescanner = new Scanner(scanner.nextLine());
                    String word;
                    int jobnum;
                    if (!linescanner.hasNext()) {
                        continue;
                    }
                    word = linescanner.next();
                    if (!word.equals("***")) {
                        continue;
                    }
                    if (!linescanner.hasNext()) {
                        continue;
                    }
                    word = linescanner.next();
                    if (word.equals("Generator")) {
                        if (!linescanner.hasNext()) {
                            continue;
                        }
                        if (generatorExpression == null) {
                            generatorExpression = linescanner.next();
                        }
                    } else if (word.equals("Job")) {
                        if (!linescanner.hasNextInt()) {
                            continue;
                        }
                        jobnum = linescanner.nextInt();
                        if (!linescanner.hasNext()) {
                            continue;
                        }
                        word = linescanner.next();
                        if (word.equals("finished")) {
                            omitted.add(jobnum);
                        }
                    }
                }
                scanner.close();
            } // Assume argument is a job generator constructor expression.
            else {
                generatorExpression = args[0];
            }

            // Create job generator.
            if (generatorExpression == null) {
                stderr.printf("Runner: No job generator in checkpoint file %s%n",
                        args[0]);
            } else {
                try {
                    generator = (JobGenerator) Instance.newInstance(generatorExpression);
                    stdout.printf("*** Generator %s%n", generatorExpression);
                    generator.omit(omitted);
                } catch (Throwable exc) {
                    stderr.printf("Runner: Could not create job generator %s%n",
                            generatorExpression);
                    exc.printStackTrace(stderr);
                }
            }
        }

        // Abort every process if job generator was not created.
        BooleanItemBuf buf = BooleanBuf.buffer(generator != null);
        world.broadcast(team.masterRank(), buf);
        if (!buf.item) {
            System.exit(1);
        }

        // Generate and run jobs.
        team.execute(new WorkerRegion() {
            public void run() throws Exception {
                execute(generator, new WorkerIteration<Job>() {
                    public void sendTaskInput(Job job, Comm comm, int wRank, int tag) {
                        stdout.printf("*** Job %d started %s %s%n",
                                job.getJobNumber(),
                                new Date(),
                                job.getDescription());
                    }

                    public void run(Job job) {
                        job.run();
                    }

                    public void receiveTaskOutput(Job job, Comm comm, int wRank, int tag) {
                        stdout.printf("*** Job %d finished %s%n",
                                job.getJobNumber(),
                                new Date());
                    }
                });
            }
        });

        if (rank == team.masterRank()) {
            stdout.printf("*** All jobs finished%n");
        }
    }

// Hidden operations.
    /**
     * Print a usage message and exit.
     */
    private static void usage() {
        stderr.println("Usage: java edu.rit.pj.job.Runner {<generator>|<file>}");
        stderr.println("<generator> = Job generator constructor expression");
        stderr.println("<file> = Checkpoint file name");
        System.exit(1);
    }

}
