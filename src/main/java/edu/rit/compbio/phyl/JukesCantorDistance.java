//******************************************************************************
//
// File:    JukesCantorDistance.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.JukesCantorDistance
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
 * Class JukesCantorDistance provides an object that computes the Jukes-Cantor
 * distance between two {@linkplain DnaSequence}s. This is the corrected
 * distance under the Jukes-Cantor model of DNA sequence evolution. The formula
 * is
 * <CENTER>
 * <I>D</I><SUB><I>JC</I></SUB> = &minus;3/4 <I>N</I> ln (1 &minus; 4/3 <I>D</I><SUB><I>H</I></SUB>/<I>N</I>)
 * </CENTER>
 * where <I>D</I><SUB><I>JC</I></SUB> is the Jukes-Cantor distance,
 * <I>D</I><SUB><I>H</I></SUB> is the Hamming distance (number of differing
 * sites), and <I>N</I> is the number of sites. For further information, see:
 * <UL>
 * <LI>
 * T. Jukes and C. Cantor. Evolution of protein molecules. In M. Munro,
 * editor. <I>Mammalian Protein Metabolism, Volume III.</I> Academic Press,
 * 1969, pages 21-132.
 * <LI>
 * J. Felsenstein. <I>Inferring Phylogenies.</I> Sinauer Associates, 2004,
 * pages 156-158.
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 23-Jul-2008
 */
public class JukesCantorDistance
	implements Distance
	{

// Exported constructors.

	/**
	 * Construct a new Jukes-Cantor distance object.
	 */
	public JukesCantorDistance()
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
		double D = seq1.distance (seq2);
		double N = seq1.length();
		double x = 1.0 - D/N/0.75;
		return
			x <= 0.0 ?
				Double.POSITIVE_INFINITY :
				Math.abs (-0.75*N*Math.log(x));
		}

	}
