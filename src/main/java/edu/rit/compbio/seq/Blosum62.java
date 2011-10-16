//******************************************************************************
//
// File:    Blosum62.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.Blosum62
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

/**
 * Class Blosum62 provides the BLOSUM-62 substitution matrix for protein
 * sequence alignment. The expression <TT>Blosum62.matrix[x][y]</TT> is the
 * score when amino acid <TT>x</TT> is aligned with amino acid <TT>y</TT>, where
 * <TT>x</TT> and <TT>y</TT> are amino acids in the range 0..27 from a
 * {@linkplain ProteinSequence}.
 * <P>
 * The BLOSUM-62 substitution matrix is taken from
 * <A HREF="http://www.ncbi.nlm.nih.gov/Class/BLAST/BLOSUM62.txt" TARGET="_top">http://www.ncbi.nlm.nih.gov/Class/BLAST/BLOSUM62.txt</A>.
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class Blosum62
	{

// Prevent construction.

	private Blosum62()
		{
		}

// Exported constants.

	/**
	 * The BLOSUM-62 protein substitution matrix.
	 */
	public static final int[][] matrix = new int[][]
		{
      /*A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z  *  -*/
/*A*/ { 4,-2, 0,-2,-1,-2, 0,-2,-1,-4,-1,-1,-1,-2,-4,-1,-1,-1, 1, 0, 0, 0,-3, 0,-2,-1,-4,-4},
/*B*/ {-2, 4,-3, 4, 1,-3,-1, 0,-3,-4, 0,-4,-3, 3,-4,-2, 0,-1, 0,-1,-3,-3,-4,-1,-3, 1,-4,-4},
/*C*/ { 0,-3, 9,-3,-4,-2,-3,-3,-1,-4,-3,-1,-1,-3,-4,-3,-3,-3,-1,-1, 9,-1,-2,-2,-2,-3,-4,-4},
/*D*/ {-2, 4,-3, 6, 2,-3,-1,-1,-3,-4,-1,-4,-3, 1,-4,-1, 0,-2, 0,-1,-3,-3,-4,-1,-3, 1,-4,-4},
/*E*/ {-1, 1,-4, 2, 5,-3,-2, 0,-3,-4, 1,-3,-2, 0,-4,-1, 2, 0, 0,-1,-4,-2,-3,-1,-2, 4,-4,-4},
/*F*/ {-2,-3,-2,-3,-3, 6,-3,-1, 0,-4,-3, 0, 0,-3,-4,-4,-3,-3,-2,-2,-2,-1, 1,-1, 3,-3,-4,-4},
/*G*/ { 0,-1,-3,-1,-2,-3, 6,-2,-4,-4,-2,-4,-3, 0,-4,-2,-2,-2, 0,-2,-3,-3,-2,-1,-3,-2,-4,-4},
/*H*/ {-2, 0,-3,-1, 0,-1,-2, 8,-3,-4,-1,-3,-2, 1,-4,-2, 0, 0,-1,-2,-3,-3,-2,-1, 2, 0,-4,-4},
/*I*/ {-1,-3,-1,-3,-3, 0,-4,-3, 4,-4,-3, 2, 1,-3,-4,-3,-3,-3,-2,-1,-1, 3,-3,-1,-1,-3,-4,-4},
/*J*/ {-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4},
/*K*/ {-1, 0,-3,-1, 1,-3,-2,-1,-3,-4, 5,-2,-1, 0,-4,-1, 1, 2, 0,-1,-3,-2,-3,-1,-2, 1,-4,-4},
/*L*/ {-1,-4,-1,-4,-3, 0,-4,-3, 2,-4,-2, 4, 2,-3,-4,-3,-2,-2,-2,-1,-1, 1,-2,-1,-1,-3,-4,-4},
/*M*/ {-1,-3,-1,-3,-2, 0,-3,-2, 1,-4,-1, 2, 5,-2,-4,-2, 0,-1,-1,-1,-1, 1,-1,-1,-1,-1,-4,-4},
/*N*/ {-2, 3,-3, 1, 0,-3, 0, 1,-3,-4, 0,-3,-2, 6,-4,-2, 0, 0, 1, 0,-3,-3,-4,-1,-2, 0,-4,-4},
/*O*/ {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4},
/*P*/ {-1,-2,-3,-1,-1,-4,-2,-2,-3,-4,-1,-3,-2,-2,-4, 7,-1,-2,-1,-1,-3,-2,-4,-2,-3,-1,-4,-4},
/*Q*/ {-1, 0,-3, 0, 2,-3,-2, 0,-3,-4, 1,-2, 0, 0,-4,-1, 5, 1, 0,-1,-3,-2,-2,-1,-1, 3,-4,-4},
/*R*/ {-1,-1,-3,-2, 0,-3,-2, 0,-3,-4, 2,-2,-1, 0,-4,-2, 1, 5,-1,-1,-3,-3,-3,-1,-2, 0,-4,-4},
/*S*/ { 1, 0,-1, 0, 0,-2, 0,-1,-2,-4, 0,-2,-1, 1,-4,-1, 0,-1, 4, 1,-1,-2,-3, 0,-2, 0,-4,-4},
/*T*/ { 0,-1,-1,-1,-1,-2,-2,-2,-1,-4,-1,-1,-1, 0,-4,-1,-1,-1, 1, 5,-1, 0,-2, 0,-2,-1,-4,-4},
/*U*/ { 0,-3, 9,-3,-4,-2,-3,-3,-1,-4,-3,-1,-1,-3,-4,-3,-3,-3,-1,-1, 9,-1,-2,-2,-2,-3,-4,-4},
/*V*/ { 0,-3,-1,-3,-2,-1,-3,-3, 3,-4,-2, 1, 1,-3,-4,-2,-2,-3,-2, 0,-1, 4,-3,-1,-1,-2,-4,-4},
/*W*/ {-3,-4,-2,-4,-3, 1,-2,-2,-3,-4,-3,-2,-1,-4,-4,-4,-2,-3,-3,-2,-2,-3,11,-2, 2,-3,-4,-4},
/*X*/ { 0,-1,-2,-1,-1,-1,-1,-1,-1,-4,-1,-1,-1,-1,-4,-2,-1,-1, 0, 0,-2,-1,-2,-1,-1,-1,-4,-4},
/*Y*/ {-2,-3,-2,-3,-2, 3,-3, 2,-1,-4,-2,-1,-1,-2,-4,-3,-1,-2,-2,-2,-2,-1, 2,-1, 7,-2,-4,-4},
/*Z*/ {-1, 1,-3, 1, 4,-3,-2, 0,-3,-4, 1,-3,-1, 0,-4,-1, 3, 0, 0,-1,-3,-2,-3,-1,-2, 4,-4,-4},
/***/ {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4},
/*-*/ {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1},
		};

	}
