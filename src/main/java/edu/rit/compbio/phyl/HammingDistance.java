//******************************************************************************
//
// File:    HammingDistance.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.HammingDistance
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

package edu.rit.compbio.phyl;

/**
 * Class HammingDistance provides an object that computes the Hamming distance
 * between two {@linkplain DnaSequence}s. The Hamming distance is the number of
 * positions at which the sequences differ.
 *
 * @author  Alan Kaminsky
 * @version 23-Jul-2008
 */
public class HammingDistance
	implements Distance
	{

// Exported constructors.

	/**
	 * Construct a new Hamming distance object.
	 */
	public HammingDistance()
		{
		}

// Exported operations.

	/**
	 * Compute the distance between the two given DNA sequences. It is assumed
	 * that the DNA sequences are the same length.
	 *
	 * @param  seq1  First DNA sequence.
	 * @param  seq2  Second DNA sequence.
	 *
	 * @return  Distance.
	 */
	public double distance
		(DnaSequence seq1,
		 DnaSequence seq2)
		{
		return seq1.distance (seq2);
		}

	}
