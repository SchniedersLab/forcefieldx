//*****************************************************************************
//
// File:    DnaSequenceList.java
// Package: edu.rit.compbio.phyl
// Unit:    Class edu.rit.compbio.phyl.DnaSequenceList
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

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Scanner;

/**
 * Class DnaSequenceList provides a list of {@linkplain DnaSequence}s. Methods
 * for reading and writing textual files of DNA sequences are provided.
 * <P>
 * Each DNA sequence consists of a sequence of <B>sites</B>. Each site has a
 * <B>state,</B> which is a set of <B>bases</B>. The four bases are adenine,
 * cytosine, guanine, and thymine. For textual I/O, each state is represented by
 * a single character as follows:
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR><TD><I>Char.</I></TD><TD WIDTH=20> </TD>
 * <TD><I>Meaning</I></TD><TD WIDTH=20> </TD>
 * <TD><I>Set</I></TD></TR>
 * <TR><TD>A</TD><TD WIDTH=20> </TD>
 * <TD>Adenine</TD><TD WIDTH=20> </TD>
 * <TD>(A)</TD></TR>
 * <TR><TD>C</TD><TD WIDTH=20> </TD>
 * <TD>Cytosine</TD><TD WIDTH=20> </TD>
 * <TD>(C)</TD></TR>
 * <TR><TD>G</TD><TD WIDTH=20> </TD>
 * <TD>Guanine</TD><TD WIDTH=20> </TD>
 * <TD>(G)</TD></TR>
 * <TR><TD>T</TD><TD WIDTH=20> </TD>
 * <TD>Thymine</TD><TD WIDTH=20> </TD>
 * <TD>(T)</TD></TR>
 * <TR><TD>Y</TD><TD WIDTH=20> </TD>
 * <TD>pYrimidine</TD><TD WIDTH=20> </TD>
 * <TD>(C or T)</TD></TR>
 * <TR><TD>R</TD><TD WIDTH=20> </TD>
 * <TD>puRine</TD><TD WIDTH=20> </TD>
 * <TD>(A or G)</TD></TR>
 * <TR><TD>W</TD><TD WIDTH=20> </TD>
 * <TD>"Weak"</TD><TD WIDTH=20> </TD>
 * <TD>(A or T)</TD></TR>
 * <TR><TD>S</TD><TD WIDTH=20> </TD>
 * <TD>"Strong"</TD><TD WIDTH=20> </TD>
 * <TD>(C or G)</TD></TR>
 * <TR><TD>K</TD><TD WIDTH=20> </TD>
 * <TD>"Keto"</TD><TD WIDTH=20> </TD>
 * <TD>(G or T)</TD></TR>
 * <TR><TD>M</TD><TD WIDTH=20> </TD>
 * <TD>"aMino"</TD><TD WIDTH=20> </TD>
 * <TD>(A or C)</TD></TR>
 * <TR><TD>B</TD><TD WIDTH=20> </TD>
 * <TD>not A</TD><TD WIDTH=20> </TD>
 * <TD>(C or G or T)</TD></TR>
 * <TR><TD>D</TD><TD WIDTH=20> </TD>
 * <TD>not C</TD><TD WIDTH=20> </TD>
 * <TD>(A or G or T)</TD></TR>
 * <TR><TD>H</TD><TD WIDTH=20> </TD>
 * <TD>not G</TD><TD WIDTH=20> </TD>
 * <TD>(A or C or T)</TD></TR>
 * <TR><TD>V</TD><TD WIDTH=20> </TD>
 * <TD>not T</TD><TD WIDTH=20> </TD>
 * <TD>(A or C or G)</TD></TR>
 * <TR><TD>X</TD><TD WIDTH=20> </TD>
 * <TD>unknown</TD><TD WIDTH=20> </TD>
 * <TD>(A or C or G or T)</TD></TR>
 * <TR><TD>-</TD><TD WIDTH=20> </TD>
 * <TD>deletion</TD><TD WIDTH=20> </TD>
 * <TD>()</TD></TR>
 * </TABLE>
 * <P>
 * The DNA sequence file format is that used by Joseph Felsenstein's Phylogeny
 * Inference Package (PHYLIP). While the file is a plain text file, it often has
 * the extension <TT>".phy"</TT> to indicate that it is in PHYLIP format. For
 * further information, see:
 * <UL>
 * <LI>
 * PHYLIP -- <A HREF="http://evolution.genetics.washington.edu/phylip/phylip.html" TARGET="_top">http://evolution.genetics.washington.edu/phylip/phylip.html</A>
 * <LI>
 * Input file format -- <A HREF="http://evolution.genetics.washington.edu/phylip/doc/sequence.html" TARGET="_top">http://evolution.genetics.washington.edu/phylip/doc/sequence.html</A>
 * </UL>
 * <P>
 * Here is an example of an input file:
 * <P>
 * <TABLE BORDER=1 CELLPADDING=4 CELLSPACING=0>
 * <TR>
 * <TD>
 * <PRE>  5    42
 * Turkey     AAGCTNGGGC ATTTCAGGGT 
 * Salmo gair AAGCCTTGGC AGTGCAGGGT 
 * H. Sapiens ACCGGTTGGC CGTTCAGGGT 
 * Chimp      AAACCCTTGC CGTTACGCTT 
 * Gorilla    AAACCCTTGC CGGTACGCTT 
 * 
 * GAGCCCGGGC AATACAGGGT AT
 * GAGCCGTGGC CGGGCACGGT AT
 * ACAGGTTGGC CGTTCAGGGT AA
 * AAACCGAGGC CGGGACACTC AT
 * AAACCATTGC CGGTACGCTT AA</PRE>
 * </TD>
 * </TR>
 * </TABLE>
 * <P>
 * The first line contains the number of species <I>S</I> and the number of
 * sites <I>N</I> in each sequence. <I>S</I> must be &gt;= 2. <I>N</I> must be
 * &gt;= 1.
 * <P>
 * The next <I>S</I> lines contain the initial data for each species. The first
 * ten characters contain the sequence name. This must be exactly ten
 * characters, padded with blanks if necessary. Then comes one character for
 * each site in the sequence. Uppercase and lowercase are considered the same.
 * Characters other than those for the states listed above are ignored. Often, a
 * blank is inserted every ten characters for readability, but this is not
 * necessary. After these <I>S</I> lines come zero or more blank lines for
 * readability, which are ignored. If there is more sequence data, the next
 * <I>S</I> lines give the states for the next sites in the sequences. This
 * continues for the rest of the file.
 * <P>
 * This is known as the "interleaved" file format. There is also a "sequential"
 * file format, but the sequential file format is not supported.
 * <P>
 * Thus, the complete sequence for each species in the example is:
 * <P>
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR><TD><I>Species</I></TD><TD WIDTH=20> </TD>
 * <TD><I>Sequence</I></TD></TR>
 * <TR><TD>Turkey</TD><TD WIDTH=20> </TD>
 * <TD><TT>AAGCTNGGGCATTTCAGGGTGAGCCCGGGCAATACAGGGTAT</TT></TD></TR>
 * <TR><TD>Salmo gair</TD><TD WIDTH=20> </TD>
 * <TD><TT>AAGCCTTGGCAGTGCAGGGTGAGCCGTGGCCGGGCACGGTAT</TT></TD></TR>
 * <TR><TD>H. Sapiens</TD><TD WIDTH=20> </TD>
 * <TD><TT>ACCGGTTGGCCGTTCAGGGTACAGGTTGGCCGTTCAGGGTAA</TT></TD></TR>
 * <TR><TD>Chimp</TD><TD WIDTH=20> </TD>
 * <TD><TT>AAACCCTTGCCGTTACGCTTAAACCGAGGCCGGGACACTCAT</TT></TD></TR>
 * <TR><TD>Gorilla</TD><TD WIDTH=20> </TD>
 * <TD><TT>AAACCCTTGCCGGTACGCTTAAACCATTGCCGGTACGCTTAA</TT></TD></TR>
 * </TABLE>
 * <P>
 * In the input file, the following alternate characters can be used: X, N, and
 * ? all mean "unknown." O (capital letter O) and - (hyphen) both mean
 * "deletion." The character . (period) means "the same as the corresponding
 * site in the first species." Here is another input file with the same
 * sequences as the one above:
 * <P>
 * <TABLE BORDER=1 CELLPADDING=4 CELLSPACING=0>
 * <TR>
 * <TD>
 * <PRE>  5    42
 * Turkey     AAGCTNGGGC ATTTCAGGGT 
 * Salmo gair ..G.CTT... AG.G...... 
 * H. Sapiens .CCGGTT... .G........ 
 * Chimp      ..A.CCTT.. .G..AC.CT. 
 * Gorilla    ..A.CCTT.. .GG.AC.CT. 
 * 
 * GAGCCCGGGC AATACAGGGT AT
 * .....GT... CGGG..C... ..
 * ACAGGTT... CG.T...... .A
 * A.A..GA... CGGGACACTC ..
 * A.A..ATT.. CGGTAC.CT. .A</PRE>
 * </TD>
 * </TR>
 * </TABLE>
 * <P>
 * Here are some more example DNA sequence files:
 * <UL>
 * <LI><A HREF="doc-files/example.phy">example.phy</A>
 * <LI><A HREF="doc-files/iguana16.phy">iguana16.phy</A>
 * <LI><A HREF="doc-files/iguana18.phy">iguana18.phy</A>
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 20-Jul-2008
 */
public class DnaSequenceList
	implements Iterable<DnaSequence>
	{

// Hidden data members.

	// DNA sequences.
	DnaSequence[] mySequence;

	// Mapping from site (index) to whether site is informative (true/false). If
	// null, must be recomputed.
	private boolean[] isInformative;

	// Number of informative sites.
	private int nInformative;

	// Number of state changes in uninformative sites.
	private int nChanges;

// Hidden constructors.

	/**
	 * Construct a new DNA sequence list.
	 */
	DnaSequenceList()
		{
		}

	/**
	 * Construct a new DNA sequence list that is a copy of the given DNA
	 * sequence list.
	 * <P>
	 * <I>Note:</I> The DNA sequences in the new list are copies of (not
	 * references to) the DNA sequences in the given list.
	 *
	 * @param  list  DNA sequence list to copy.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>list</TT> is null.
	 */
	public DnaSequenceList
		(DnaSequenceList list)
		{
		int N = list.mySequence.length;
		this.mySequence = new DnaSequence [N];
		for (int i = 0; i < N; ++ i)
			{
			this.mySequence[i] = new DnaSequence (list.mySequence[i]);
			}
		if (list.isInformative != null)
			{
			this.isInformative = (boolean[]) list.isInformative.clone();
			}
		this.nInformative = list.nInformative;
		this.nChanges = list.nChanges;
		}

// Exported operations.

	/**
	 * Obtain this DNA sequence list's length.
	 *
	 * @return  Length <I>N</I> (number of DNA sequences).
	 */
	public int length()
		{
		return mySequence.length;
		}

	/**
	 * Get the DNA sequence at the given index in this DNA sequence list.
	 *
	 * @param  i  Index, 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 *
	 * @return  DNA sequence.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public DnaSequence seq
		(int i)
		{
		return mySequence[i];
		}

	/**
	 * Read a DNA sequence list from the given input file. The input file must
	 * be in interleaved PHYLIP format.
	 * <P>
	 * The DNA sequences' sites and names are read from the input file. The DNA
	 * sequences' scores are set to 0.
	 *
	 * @param  file  File.
	 *
	 * @return  DNA sequence list.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred. Thrown if the input file's contents
	 *     were invalid.
	 */
	public static DnaSequenceList read
		(File file)
		throws IOException
		{
		Scanner filescanner = new Scanner (file);
		Scanner linescanner;
		int S, N;
		DnaSequenceList list;
		int[] sitecount;
		String line;

		try
			{
			// Read number of species and number of sites from first line.
			if (! filescanner.hasNextLine())
				{
				throw new IOException
					("DnaSequenceList.read(\"" + file + "\"): " +
					 "Empty file");
				}
			linescanner = new Scanner (filescanner.nextLine());
			if (! linescanner.hasNextInt())
				{
				throw new IOException
					("DnaSequenceList.read(\"" + file + "\"): " +
					 "Number of species invalid or missing");
				}
			S = linescanner.nextInt();
			if (S < 2)
				{
				throw new IOException
					("DnaSequenceList.read(\"" + file + "\"): " +
					 "Number of species must be >= 2");
				}
			if (! linescanner.hasNextInt())
				{
				throw new IOException
					("DnaSequenceList.read(\"" + file + "\"): " +
					 "Number of sites invalid or missing");
				}
			N = linescanner.nextInt();
			if (N < 1)
				{
				throw new IOException
					("DnaSequenceList.read(\"" + file + "\"): " +
					 "Number of sites must be >= 1");
				}

			// Set up DNA sequence list and site count array.
			list = new DnaSequenceList();
			list.mySequence = new DnaSequence [S];
			sitecount = new int [S];

			// Read sequence data from groups of S lines until EOF.
			fileloop: for (;;)
				{
				speciesloop: for (int s = 0; s < S; ++ s)
					{
					// Get a line of sequence data for species s.
					if (filescanner.hasNextLine())
						{
						}
					else if (s != 0 || sitecount[s] == 0)
						{
						throw new IOException
							("DnaSequenceList.read(\"" + file + "\"): " +
							 "Missing a line of sequence data for species " +
							 (s+1));
						}
					else
						{
						break fileloop;
						}
					line = filescanner.nextLine();

					// Ignore blank lines.
					if (line.trim().equals (""))
						{
						-- s;
						continue;
						}

					// The first time, extract sequence name and create
					// DnaSequence object.
					if (sitecount[s] == 0)
						{
						if (line.length() < 10)
							{
							throw new IOException
								("DnaSequenceList.read(\"" + file + "\"): " +
								 "Name must be 10 characters for species " +
								 (s+1));
							}
						list.mySequence[s] =
							new DnaSequence
								(N, 0, line.substring (0, 10) .trim());
						line = line.substring (10);
						}

					// Parse characters in sequence data.
					int len = line.length();
					byte[] seq = list.mySequence[s].mySites;
					byte[] seq0 = list.mySequence[0].mySites;
					int count = sitecount[s];
					for (int i = 0; i < len; ++ i)
						{
						switch (line.charAt(i))
							{
							case 'O': case 'o': case '-':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  0; // ----
								++ count;
								break;
							case 'A': case 'a':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  1; // ---A
								++ count;
								break;
							case 'C': case 'c':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  2; // --C-
								++ count;
								break;
							case 'M': case 'm':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  3; // --CA
								++ count;
								break;
							case 'G': case 'g':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  4; // -G--
								++ count;
								break;
							case 'R': case 'r':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  5; // -G-A
								++ count;
								break;
							case 'S': case 's':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  6; // -GC-
								++ count;
								break;
							case 'V': case 'v':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  7; // -GCA
								++ count;
								break;
							case 'T': case 't':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  8; // T---
								++ count;
								break;
							case 'W': case 'w':
								verifyCount (count, N, file, s);
								seq[count] = (byte)  9; // T--A
								++ count;
								break;
							case 'Y': case 'y':
								verifyCount (count, N, file, s);
								seq[count] = (byte) 10; // T-C-
								++ count;
								break;
							case 'H': case 'h':
								verifyCount (count, N, file, s);
								seq[count] = (byte) 11; // T-CA
								++ count;
								break;
							case 'K': case 'k':
								verifyCount (count, N, file, s);
								seq[count] = (byte) 12; // TG--
								++ count;
								break;
							case 'D': case 'd':
								verifyCount (count, N, file, s);
								seq[count] = (byte) 13; // TG-A
								++ count;
								break;
							case 'B': case 'b':
								verifyCount (count, N, file, s);
								seq[count] = (byte) 14; // TGC-
								++ count;
								break;
							case 'X': case 'x': case 'N': case 'n': case '?':
								verifyCount (count, N, file, s);
								seq[count] = (byte) 15; // TGCA
								++ count;
								break;
							case '.':
								verifyCount (count, N, file, s);
								if (s == 0)
									{
									throw new IOException
										("DnaSequenceList.read(\"" + file +
										 "\"): " +
										 "'.' not allowed in species 1");
									}
								if (count >= sitecount[0])
									{
									throw new IOException
										("DnaSequenceList.read(\"" + file +
										 "\"): " +
										 "'.' in species " + (s+1) +
										 " has no corresponding site in species 1");
									}
								seq[count] = seq0[count];
								++ count;
								break;
							}
						}
					sitecount[s] = count;
					}
				}

			// Verify correct site count for all species.
			for (int s = 0; s < S; ++ s)
				{
				if (sitecount[s] < N)
					{
					throw new IOException
						("DnaSequenceList.read(\"" + file + "\"): " +
						 "Too few sites for species " + (s+1));
					}
				else if (sitecount[s] > N)
					{
					throw new IOException
						("DnaSequenceList.read(\"" + file + "\"): " +
						 "Too many sites for species " + (s+1));
					}
				}

			// Return DNA sequence list.
			return list;
			}

		finally
			{
			filescanner.close();
			}
		}

	private static void verifyCount
		(int count,
		 int N,
		 File file,
		 int s)
		throws IOException
		{
		if (count >= N)
			{
			throw new IOException
				("DnaSequenceList.read(\"" + file + "\"): " +
				 "Too many sites for species " + (s+1));
			}
		}

	/**
	 * Write this DNA sequence list to the given output file. The output file is
	 * in interleaved PHYLIP format. There are 70 sites on each output line.
	 * Periods are not used. Informative sites are not marked in bold.
	 *
	 * @param  file  File.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void write
		(File file)
		throws IOException
		{
		write (file, 70, false, false);
		}

	/**
	 * Write this DNA sequence list to the given output file. The output file is
	 * in interleaved PHYLIP format.
	 *
	 * @param  file     File.
	 * @param  sites    Number of sites per output line.
	 * @param  periods  True to use periods, false not to use periods.
	 * @param  bold     True to mark informative sites in bold, false not to.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>file</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>sites</TT> &lt;= 10.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void write
		(File file,
		 int sites,
		 boolean periods,
		 boolean bold)
		throws IOException
		{
		PrintStream ps =
			new PrintStream
				(new BufferedOutputStream
					(new FileOutputStream (file)));
		try
			{
			write (ps, sites, periods, bold);
			}
		finally
			{
			ps.close();
			}
		}

	/**
	 * Write this DNA sequence list to the given print stream in interleaved
	 * PHYLIP format.
	 *
	 * @param  ps       Print stream.
	 * @param  sites    Number of sites per output line.
	 * @param  periods  True to use periods, false not to.
	 * @param  bold     True to mark informative sites in bold, false not to.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>ps</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>sites</TT> &lt;= 10.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void write
		(PrintStream ps,
		 int sites,
		 boolean periods,
		 boolean bold)
		throws IOException
		{
		if (sites <= 10)
			{
			throw new IllegalArgumentException
				("DnaSequenceList.write(): sites = " + sites + " illegal");
			}

		// Determine informative sites if necessary.
		if (bold) computeInformativeSites();

		// Print number of species and number of sites.
		int S = mySequence.length;
		int N = mySequence[0].myLength;
		ps.print (S);
		ps.print (' ');
		ps.print (N);
		ps.println();

		// Print groups of sites for each species. On the first line, print
		// sequence name, padded or truncated to 10 characters.
		int lb = 0;
		int ub = Math.min (sites-10, N);
		byte[] seq0 = mySequence[0].mySites;
		while (lb < N)
			{
			for (int s = 0; s < S; ++ s)
				{
				byte[] seq = mySequence[s].mySites;
				if (lb == 0) ps.print (padName (mySequence[s].myName));
				for (int i = lb; i < ub; ++ i)
					{
					if ((lb == 0 || i > lb) && i % 10 == 0)
						{
						ps.print (' ');
						}
					if (periods && s > 0 && seq[i] == seq0[i])
						{
						printSite (ps, i, '.', bold);
						}
					else
						{
						printSite
							(ps, i, DnaSequence.state2char[seq[i]], bold);
						}
					}
				ps.println();
				}
			ps.println();
			lb = ub;
			ub = Math.min (ub+sites, N);
			}

		// Check for I/O errors.
		if (ps.checkError())
			{
			throw new IOException ("DnaSequenceList.write(): I/O error");
			}
		}

	private static String padName
		(String name)
		{
		if (name == null) return "<unnamed> ";
		int len = name.length();
		if (len == 10)
			{
			return name;
			}
		else if (len > 10)
			{
			return name.substring (0, 10);
			}
		else
			{
			return name + padding[len];
			}
		}

	private static String[] padding = new String[]
		{/*0*/ "          ",
		 /*1*/ "         ",
		 /*2*/ "        ",
		 /*3*/ "       ",
		 /*4*/ "      ",
		 /*5*/ "     ",
		 /*6*/ "    ",
		 /*7*/ "   ",
		 /*8*/ "  ",
		 /*9*/ " "};

	private void printSite
		(PrintStream ps,
		 int i,
		 char c,
		 boolean bold)
		{
		if (bold && isInformative[i])
			{
			ps.print ("<B>");
			ps.print (c);
			ps.print ("</B>");
			}
		else
			{
			ps.print (c);
			}
		}

	/**
	 * Truncate this DNA sequence list to the given length. If this list is
	 * already shorter than <TT>len</TT>, the <TT>truncate()</TT> method does
	 * nothing.
	 *
	 * @param  len  Length.
	 *
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>len</TT> &lt; 0.
	 */
	public void truncate
		(int len)
		{
		if (len < mySequence.length)
			{
			DnaSequence[] newSequence = new DnaSequence [len];
			System.arraycopy (mySequence, 0, newSequence, 0, len);
			mySequence = newSequence;
			}
		}

	/**
	 * Excise uninformative sites from the DNA sequences in this DNA sequence
	 * list.
	 * <P>
	 * Each site in the DNA sequences is either "uninformative" or
	 * "informative," defined as follows:
	 * <UL>
	 * <LI>
	 * If the site has the same state (A, C, G, or T) in all sequences, the
	 * site is uninformative. This site will contribute no state changes to the
	 * parsimony score in every possible phylogenetic tree.
	 * <P><LI>
	 * If the site has the same state in all sequences, except for one or more
	 * sequences that have a unique state at that site (i.e., a state that
	 * appears in no other sequences at that site), the site is uninformative.
	 * The site will contribute the same number of state changes to the
	 * parsimony score in every possible phylogenetic tree, namely the number of
	 * different states that appear at that site, minus 1.
	 * <P><LI>
	 * Otherwise, the site is informative. There are at least two different
	 * states at that site, and each state appears in at least two different
	 * sequences. The site will contribute a different number of state changes
	 * to the parsimony score, depending on where the sequences appear in the
	 * phylogenetic tree.
	 * </UL>
	 * <P>
	 * Since the uninformative sites do not affect the outcome of a maximum
	 * parsimony phylogenetic tree search, the uninformative sites can be
	 * omitted from the tree scoring process to save time. The informative sites
	 * do affect the outcome and must be included in the tree scoring process.
	 * <P>
	 * The <TT>exciseUninformativeSites()</TT> removes the uninformative sites
	 * from the DNA sequences in this list. The DNA sequences' scores and names
	 * are unchanged.
	 *
	 * @return  Number of state changes the (excised) uninformative sites
	 *          contribute to the parsimony score.
	 */
	public int exciseUninformativeSites()
		{
		int S = mySequence.length;
		int N = mySequence[0].length();

		// Determine which sites are informative.
		computeInformativeSites();

		// Excise uninformative sites from sequences.
		for (int s = 0; s < S; ++ s)
			{
			byte[] oldSites = mySequence[s].mySites;
			mySequence[s] =
				new DnaSequence
					(nInformative,
					 mySequence[s].myScore,
					 mySequence[s].myName);
			byte[] excSites = mySequence[s].mySites;
			int j = 0;
			for (int i = 0; i < N; ++ i)
				{
				if (isInformative[i])
					{
					excSites[j++] = oldSites[i];
					}
				}
			}

		// Mark all sites as informative.
		isInformative = new boolean [nInformative];
		Arrays.fill (isInformative, true);

		// Return number of state changes.
		return nChanges;
		}

	/**
	 * Returns the number of informative sites in this DNA sequence list.
	 *
	 * @return  Number of informative sites.
	 */
	public int informativeSiteCount()
		{
		computeInformativeSites();
		return nInformative;
		}

	/**
	 * Determine the number of absent states after adding each sequence in this
	 * DNA sequence list to a tree. The return value <I>A</I> is an
	 * <I>N</I>-element array, where <I>N</I> is the length of this DNA sequence
	 * list. As sequences from this list are added to a tree in order from
	 * <I>i</I> = 0 to <I>N</I>&minus;1, <I>A</I>[<I>i</I>] is the number of
	 * character states that do not yet appear in the tree. Thus, the number of
	 * state changes in the tree must increase by at least <I>A</I>[<I>i</I>]
	 * when the sequences after sequence <I>i</I> are added to the tree. This
	 * can be used to prune a branch-and-bound search.
	 *
	 * @return  Array <I>A</I>.
	 */
	public int[] countAbsentStates()
		{
		int N = mySequence.length;
		int L = mySequence[0].length();
		int[] A = new int [N];

		// Compute the union of all the DNA sequences.
		byte[] sites = new byte [L];
		for (int i = 0; i < N; ++ i)
			{
			byte[] mysites_i = mySequence[i].mySites;
			for (int j = 0; j < L; ++ j)
				{
				sites[j] |= mysites_i[j];
				}
			}

		// Subtract each sequence from the union, count and record states.
		for (int i = 0; i < N; ++ i)
			{
			byte[] mysites_i = mySequence[i].mySites;
			int count = 0;
			for (int j = 0; j < L; ++ j)
				{
				sites[j] &= ~ mysites_i[j];
				count += DnaSequence.state2bitCount [sites[j]];
				}
			A[i] = count;
			}

		return A;
		}

	/**
	 * Create a DNA sequence tree from this DNA sequence list and the given tree
	 * signature. The tree signature is an array of indexes of length <I>N</I>,
	 * where <I>N</I> is the length of this list. To construct the tree, for all
	 * <I>i</I> from 0 to <I>N</I>&minus;1, the DNA sequence at index <I>i</I>
	 * in this list is added to the tree at index <TT>signature[i]</TT> using
	 * the <TT>DnaSequenceTree.add()</TT> method. For all <I>i</I>,
	 * <TT>signature[i]</TT> must be in the range 0 ..
	 * 2(<I>i</I>&nbsp;&minus;&nbsp;1), except <TT>signature[0]</TT> is 0.
	 * <P>
	 * <I>Note:</I> The returned tree has references to (not copies of) the DNA
	 * sequences in this list.
	 *
	 * @param  signature  Tree signature (array of tree indexes).
	 *
	 * @return  Tree.
	 */
	public DnaSequenceTree toTree
		(int[] signature)
		{
		int N = mySequence.length;
		DnaSequenceTree tree = new DnaSequenceTree (2*N - 1);
		for (int i = 0; i < N; ++ i)
			{
			tree.add (signature[i], mySequence[i]);
			}
		return tree;
		}

	/**
	 * Returns an iterator for the DNA sequences in this list.
	 *
	 * @return  Iterator.
	 */
	public Iterator<DnaSequence> iterator()
		{
		return new Iterator<DnaSequence>()
			{
			int i = 0;

			public boolean hasNext()
				{
				return i < mySequence.length;
				}

			public DnaSequence next()
				{
				return mySequence[i++];
				}

			public void remove()
				{
				throw new UnsupportedOperationException();
				}
			};
		}

// Hidden operations.

	/**
	 * Compute information about informative sites.
	 */
	private void computeInformativeSites()
		{
		if (isInformative != null) return;

		int S = mySequence.length;
		int N = mySequence[0].length();

		// Allocate storage to remember each site's category: true =
		// informative, false = uninformative. Also count number of informative
		// sites and number of state changes in uninformative sites.
		isInformative = new boolean [N];
		nInformative = 0;
		nChanges = 0;

		// Allocate storage to count states at each site.
		int[] stateCount = new int [16];

		// Examine all sites.
		for (int i = 0; i < N; ++ i)
			{
			Arrays.fill (stateCount, 0);

			// Examine current site in all sequences.
			for (int s = 0; s < S; ++ s)
				{
				++ stateCount[mySequence[s].mySites[i]];
				}

			// Count how many values in stateCount are 2 or greater.
			int x = 0;
			for (int j = 0; j < 16; ++ j)
				{
				if (stateCount[j] >= 2) ++ x;
				}

			// Categorize current site.
			if (x >= 2)
				{
				// Informative site.
				isInformative[i] = true;
				++ nInformative;
				}
			else
				{
				// Uninformative site. Increase number of state changes by
				// (number of different states - 1).
				isInformative[i] = false;
				for (int j = 0; j < 16; ++ j)
					{
					if (stateCount[j] > 0) ++ nChanges;
					}
				-- nChanges;
				}
			}
		}

	}
