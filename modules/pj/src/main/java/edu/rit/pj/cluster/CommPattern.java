//******************************************************************************
//
// File:    CommPattern.java
// Package: edu.rit.pj.cluster
// Unit:    Class edu.rit.pj.cluster.CommPattern
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

import java.util.ArrayList;

/**
 * Class CommPattern provides static methods for calculating communication
 * patterns for collective communication operations.
 *
 * @author Alan Kaminsky
 * @version 15-Mar-2008
 */
public class CommPattern {

// Prevent construction.
    private CommPattern() {
    }

// Exported operations.
    /**
     * Calculate the communication pattern for a parallel broadcast tree. This
     * is also used in reverse for a parallel reduction tree.
     *
     * @param size Size of the communicator. Must be &gt;= 1.
     * @param rank Rank of this process in the communicator. Must be in the
     * range 0 .. <TT>size</TT>-1.
     * @param root Rank of the root process for the broadcast. Must be in the
     * range 0 .. <TT>size</TT>-1.
     * @return Array of process ranks for the parallel broadcast pattern. The
     * element at index 0 is the parent process rank, or -1 if there is no
     * parent process. The elements at indexes 1 and above, if any, are the
     * child process ranks.
     * @exception IllegalArgumentException (unchecked exception) Thrown if any
     * argument is illegal.
     */
    public static int[] broadcastPattern(int size,
            int rank,
            int root) {
        // Verify preconditions.
        if (size < 1) {
            throw new IllegalArgumentException("broadcastPattern(): size must be >= 1");
        }
        if (0 > rank || rank >= size) {
            throw new IllegalArgumentException("broadcastPattern(): rank must be in the range 0 .. "
                    + (size - 1));
        }
        if (0 > root || root >= size) {
            throw new IllegalArgumentException("broadcastPattern(): root must be in the range 0 .. "
                    + (size - 1));
        }

        // Imagine for the moment that the processes are numbered 1 through K,
        // where K = size. The broadcast communication pattern takes place in a
        // number of rounds. The rounds are numbered 1, 2, 4, 8, and so on. In
        // round 1, process 1 sends to process 2. In round 2, processes 1 and 2
        // send to processes 3 and 4 in parallel. In round 4, processes 1, 2, 3,
        // and 4 send to processes 5, 6, 7, and 8 in parallel:
        //
        // Process
        // 1    2    3    4    5    6    7    8
        // |    |    |    |    |    |    |    |
        // |    |    |    |    |    |    |    |
        // |--->|    |    |    |    |    |    |  Round 1
        // |    |    |    |    |    |    |    |
        // |    |    |    |    |    |    |    |
        // |-------->|    |    |    |    |    |  Round 2
        // |    |-------->|    |    |    |    |
        // |    |    |    |    |    |    |    |
        // |    |    |    |    |    |    |    |
        // |------------------>|    |    |    |  Round 4
        // |    |------------------>|    |    |
        // |    |    |------------------>|    |
        // |    |    |    |------------------>|
        // |    |    |    |    |    |    |    |
        // |    |    |    |    |    |    |    |
        //
        // In general, in round i, processes 1 through i send to processes 1+i
        // through i+i in parallel. This continues until there are no more
        // processes to send to.
        //
        // After calculating the above pattern, the process numbers are shifted
        // such that process 1 becomes process root, process 2 becomes process
        // root+1, and so on. In general, process i becomes process (i+root-1)
        // (mod size).
        // This process's rank relative to the root, in the range 1 .. size.
        int thisrank
                = rank >= root
                ? rank - root + 1
                : rank + size - root + 1;

        // Parent process.
        int parent = -1;

        // List of child processes.
        ArrayList<Integer> childlist = new ArrayList<Integer>();

        // Do all rounds.
        int round = 1;
        while (round < size) {
            // Do all messages within this round.
            for (int src = 1; src <= round; ++src) {
                int dst = src + round;

                // If this process is the destination, record source as parent.
                if (thisrank == dst) {
                    parent = src;
                } // If this process is the source, record destination as child.
                else if (thisrank == src && dst <= size) {
                    childlist.add(dst);
                }
            }

            // Next round.
            round <<= 1;
        }

        // Make an array to hold parent and child processes.
        int n = childlist.size();
        int[] result = new int[n + 1];

        // Record parent, offsetting rank.
        result[0] = parent == -1 ? -1 : (parent + root - 1) % size;

        // Record children, offsetting ranks.
        for (int i = 0; i < n; ++i) {
            result[i + 1] = (childlist.get(i) + root - 1) % size;
        }

        // All done!
        return result;
    }

// Unit test main program.
//	/**
//	 * Unit test main program.
//	 * <P>
//	 * Usage: java edu.rit.pj.cluster.CommPattern <I>size</I> <I>rank</I>
//	 * <I>root</I>
//	 */
//	public static void main
//		(String[] args)
//		{
//		if (args.length != 3) usage();
//		int size = Integer.parseInt (args[0]);
//		int rank = Integer.parseInt (args[1]);
//		int root = Integer.parseInt (args[2]);
//		int[] pattern = CommPattern.broadcastPattern (size, rank, root);
//		System.out.print ("broadcastPattern (size=");
//		System.out.print (size);
//		System.out.print (", rank=");
//		System.out.print (rank);
//		System.out.print (", root=");
//		System.out.print (root);
//		System.out.print (") returns {");
//		System.out.print (pattern[0]);
//		for (int i = 1; i < pattern.length; ++ i)
//			{
//			System.out.print (", ");
//			System.out.print (pattern[i]);
//			}
//		System.out.println ("}");
//		}
//
//	/**
//	 * Print a usage message and exit.
//	 */
//	private static void usage()
//		{
//		System.err.println ("Usage: java edu.rit.pj.cluster.CommPattern <size> <rank> <root>");
//		System.exit (1);
//		}
}
