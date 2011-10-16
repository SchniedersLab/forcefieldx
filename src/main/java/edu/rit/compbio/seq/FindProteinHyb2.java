//******************************************************************************
//
// File:    FindProteinHyb2.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.FindProteinHyb2
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
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************

package edu.rit.compbio.seq;

import edu.rit.mp.ObjectBuf;

import edu.rit.mp.buf.ObjectItemBuf;

import edu.rit.pj.Comm;
import edu.rit.pj.CommStatus;
import edu.rit.pj.LongForLoop;
import edu.rit.pj.LongSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;

import edu.rit.util.LongRange;

import java.io.File;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Class FindProteinHyb2 is a hybrid parallel program that finds matches for a
 * protein in a protein sequence database. The query sequence is stored in a
 * file in FASTA format; for further information, see class {@linkplain
 * ProteinSequence}. The protein sequence database is stored in two files, a
 * database file in FASTA format and an index file; for further information, see
 * class {@linkplain ProteinDatabase}. The program uses the Smith-Waterman
 * algorithm to compute a local alignment between the query sequence and each
 * subject sequence in the database. The program uses the BLOSUM-62 protein
 * substitution matrix. The program uses affine gap penalties with a gap
 * existence penalty of &minus;11 and a gap extension penalty of &minus;1. The
 * program prints on the standard output the resulting alignments from highest
 * to lowest score. The program only prints alignments with an <I>E</I>-value
 * below the given threshold; if not specified, the default is 10.
 * <P>
 * The program uses the master-worker pattern for load balancing. The program
 * partitions the alignments among the processes using the schedule specified by
 * the <TT>-Dpj.schedule</TT> flag. The master sends a range of database indexes
 * to a worker; the worker sends back to the master a list of alignments for
 * those database indexes; and the process repeats. Within each worker, the
 * alignments are partitioned among multiple threads using the parallel loop
 * schedule specified by the fourth command line argument. If this argument is
 * missing, the default is to divide the alignments evenly among the threads.
 * <P>
 * Usage: java -Dpj.np=<I>Kp</I> -Dpj.nt=<I>Kt</I> [
 * -Dpj.schedule=<I>procschedule</I> ] edu.rit.compbio.seq.FindProteinHyb2
 * <I>queryfile</I> <I>databasefile</I> <I>indexfile</I> [ <I>thrschedule</I> [
 * <I>expect</I> ] ]
 * <BR><I>Kp</I> = Number of parallel processes
 * <BR><I>Kt</I> = Number of parallel threads per process
 * <BR><I>procschedule</I> = Load balancing schedule for processes (default:
 * fixed schedule)
 * <BR><I>queryfile</I> = Query sequence file
 * <BR><I>databasefile</I> = Protein sequence database file
 * <BR><I>indexfile</I> = Protein sequence index file
 * <BR><I>thrschedule</I> = Load balancing schedule for threads (default: fixed
 * schedule)
 * <BR><I>expect</I> = <I>E</I>-value threshold (default: 10)
 *
 * @author  Alan Kaminsky
 * @version 07-Jul-2008
 */
public class FindProteinHyb2
	{

// Prevent construction.

	private FindProteinHyb2()
		{
		}

// Global variables.

	// World communicator.
	static Comm world;
	static int size;
	static int rank;

	// Command line arguments.
	static File queryfile;
	static File databasefile;
	static File indexfile;
	static LongSchedule thrschedule;
	static double expect;

	// Query sequence.
	static ProteinSequence query;

	// Protein sequence database.
	static ProteinDatabase database;

	// Object to compute alignment statistics.
	static AlignmentStats stats;

	// List of alignments found.
	static List<Alignment> alignmentsFound;

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Start timing.
		long t1 = System.currentTimeMillis();

		// Initialize world communicator.
		Comm.init (args);
		world = Comm.world();
		size = world.size();
		rank = world.rank();

		// Parse command line arguments.
		if (3 > args.length || args.length > 5) usage();
		queryfile = new File (args[0]);
		databasefile = new File (args[1]);
		indexfile = new File (args[2]);
		thrschedule = LongSchedule.fixed();
		if (args.length >= 4) thrschedule = LongSchedule.parse (args[3]);
		expect = 10.0;
		if (args.length >= 5) expect = Double.parseDouble (args[4]);

		// Set up query sequence.
		query = new ProteinSequence (queryfile);

		// Set up protein sequence database.
		database = new ProteinDatabase (databasefile, indexfile);

		// Set up object to compute alignment statistics.
		stats = new DefaultAlignmentStats (database.getDatabaseLength());

		// Set up list of alignments found.
		alignmentsFound = new ArrayList<Alignment>();

		long t2 = System.currentTimeMillis();

		// In process 0, run the master and the worker in separate threads.
		if (rank == 0)
			{
			new ParallelTeam(2).execute (new ParallelRegion()
				{
				public void run() throws Exception
					{
					execute (new ParallelSection()
						{
						public void run() throws Exception
							{
							masterSection();
							}
						},
					new ParallelSection()
						{
						public void run() throws Exception
							{
							workerSection();
							}
						});
					}
				});
			}

		// In processes 1 and up, run just the worker section.
		else
			{
			workerSection();
			}

		long t3 = System.currentTimeMillis();

		// Process 0 does the postprocessing.
		if (rank == 0)
			{
			// Sort alignments into descending order of score.
			Collections.sort (alignmentsFound);

			// Set up alignment printer.
			AlignmentPrinter printer = new AlignmentPrinter (System.out, stats);

			// Print query sequence.
			System.out.println ("Query Description:");
			System.out.println (query.description());
			System.out.println ("Length = "+query.length());
			System.out.println();

			// Print summary of each alignment.
			System.out.println ("                                                                Bit  E-");
			System.out.println ("Subject Description                                           Score  Value");
			for (Alignment a : alignmentsFound)
				{
				printer.printSummary
					(a, database.getProteinSequence (a.getSubjectId()));
				}
			System.out.println();

			// Print details of each alignment.
			for (Alignment a : alignmentsFound)
				{
				printer.printDetails
					(a, query, database.getProteinSequence (a.getSubjectId()));
				}

			// Print various information about the alignment procedure.
			System.out.println
				("Query file: "+queryfile);
			System.out.println
				("Database file: "+databasefile);
			System.out.println
				("Database index file: "+indexfile);
			System.out.println
				("Number of sequences: "+database.getProteinCount());
			System.out.println
				("Number of matches: "+alignmentsFound.size());
			System.out.println
				("Query length: "+query.length());
			System.out.println
				("Database length: "+database.getDatabaseLength());
			stats.print (System.out);
			System.out.println();

			long t4 = System.currentTimeMillis();
			System.out.println ((t2-t1)+" msec pre");
			System.out.println ((t3-t2)+" msec calc");
			System.out.println ((t4-t3)+" msec post");
			System.out.println ((t4-t1)+" msec total");
			}

		// Done using the protein sequence database.
		database.close();
		}

// Hidden operations.

	/**
	 * Perform the master section.
	 *
	 * @exception  Exception
	 *     Thrown if an I/O error occurred.
	 */
	private static void masterSection()
		throws IOException
		{
		int worker;
		LongRange range;

		long t2 = System.currentTimeMillis();

		// Set up a schedule object.
		LongSchedule schedule = LongSchedule.runtime();
		schedule.start (size, new LongRange (0, database.getProteinCount()-1));

		// Send initial database index range to each worker. If range is null,
		// no more work for that worker. Keep count of active workers.
		int activeWorkers = size;
		for (worker = 0; worker < size; ++ worker)
			{
			range = schedule.next (worker);
			world.send (worker, ObjectBuf.buffer (range));
			if (range == null) -- activeWorkers;
			}

		// Repeat until all workers have finished.
		while (activeWorkers > 0)
			{
			// Receive a message containing a list of zero or more alignments
			// from any worker.
			ObjectItemBuf<List<Alignment>> buf = ObjectBuf.buffer();
			CommStatus status = world.receive (null, buf);
			worker = status.fromRank;

			// Send next database index range to that specific worker. If null,
			// no more work.
			range = schedule.next (worker);
			world.send (worker, ObjectBuf.buffer (range));
			if (range == null) -- activeWorkers;

			// Add alignments to list.
			alignmentsFound.addAll (buf.item);
			}
		}

	/**
	 * Perform the worker section.
	 *
	 * @exception  Exception
	 *     Thrown if an I/O error occurred.
	 */
	private static void workerSection()
		throws Exception
		{
		// Set up parallel team.
		ParallelTeam team = new ParallelTeam();
		int K = team.getThreadCount();

		// Set up per-thread objects to perform alignments.
		final ProteinLocalAlignment[] aligner = new ProteinLocalAlignment [K];
		for (int i = 0; i < K; ++ i)
			{
			aligner[i] = new ProteinLocalAlignmentSeq();
			aligner[i].setQuerySequence (query, 0);
			}

		// Process chunks from master.
		for (;;)
			{
			// Receive database index range from master. If null, no more work.
			ObjectItemBuf<LongRange> rangeBuf = ObjectBuf.buffer();
			world.receive (0, rangeBuf);
			LongRange range = rangeBuf.item;
			if (range == null) break;
			final long lb = range.lb();
			final long ub = range.ub();

			// Set up list to hold alignments.
			final List<Alignment> alignments = new ArrayList<Alignment>();

			// Align query sequence against every subject sequence.
			team.execute (new ParallelRegion()
				{
				public void run() throws Exception
					{
					execute (lb, ub, new LongForLoop()
						{
						// Per-thread variables plus padding.
						ProteinLocalAlignment thrAligner;
						List<Alignment> thrAlignments;
						long p0, p1, p2, p3, p4, p5, p6, p7;
						long p8, p9, pa, pb, pc, pd, pe, pf;

						// Initialize per-thread variables.
						public void start()
							{
							thrAligner = aligner[getThreadIndex()];
							thrAlignments = new ArrayList<Alignment>();
							}

						// Use thread-level loop schedule.
						public LongSchedule schedule()
							{
							return thrschedule;
							}

						// Do alignments.
						public void run (long first, long last) throws Exception
							{
							for (long id = first; id <= last; ++ id)
								{
								ProteinSequence subject =
									database.getProteinSequence (id);
								thrAligner.setSubjectSequence (subject, id);
								Alignment a = thrAligner.align();
								if (stats.eValue (a) <= expect)
									{
									thrAlignments.add (a);
									}
								}
							}

						// Reduce per-thread alignments into global alignments.
						public void finish() throws Exception
							{
							region().critical (new ParallelSection()
								{
								public void run()
									{
									alignments.addAll (thrAlignments);
									}
								});
							}
						});
					}
				});

			// Send alignments back to master.
			world.send (0, ObjectBuf.buffer (alignments));
			}
		};

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java -Dpj.np=<Kp> -Dpj.nt=<Kt> [-Dpj.schedule=<procschedule>] edu.rit.compbio.seq.FindProteinHyb2 <queryfile> <databasefile> <indexfile> [<thrschedule> [<expect>]]");
		System.err.println ("<Kp> = Number of parallel processes");
		System.err.println ("<Kt> = Number of parallel threads per process");
		System.err.println ("<procschedule> = Load balancing schedule for processes (default: fixed schedule)");
		System.err.println ("<queryfile> = Query sequence file");
		System.err.println ("<databasefile> = Protein sequence database file");
		System.err.println ("<indexfile> = Protein sequence index file");
		System.err.println ("<thrschedule> = Load balancing schedule for threads (default: fixed schedule)");
		System.err.println ("<expect> = <E>-value threshold (default: 10)");
		System.exit (1);
		}

	}
