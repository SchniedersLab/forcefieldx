//******************************************************************************
//
// File:    ProteinSequence.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.ProteinSequence
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

import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Class ProteinSequence encapsulates a protein sequence. A protein sequence
 * object may be constructed from a string or read from a file. Protein sequence
 * objects may also be read from a protein sequence database using class
 * {@linkplain ProteinDatabase}.
 * <P>
 * In a file, a protein sequence is stored in FASTA format. A FASTA format
 * protein sequence consists of one description line and one or more sequence
 * lines. The description line consists of an initial <TT>'&gt;'</TT> character
 * followed by zero or more characters (the protein's description). A sequence
 * line consists of one or more characters <TT>'A'</TT> through <TT>'Z'</TT>,
 * <TT>'a'</TT> through <TT>'z'</TT>, <TT>'*'</TT>, or <TT>'-'</TT>.
 * <P>
 * In a program, a protein sequence is represented as a byte array (type
 * <TT>byte[]</TT>). For a protein sequence of length <I>L</I>, the byte array
 * contains <I>L</I>+1 bytes. The byte at index 0 is unused and contains a value
 * of -1. The bytes at indexes 1 through <I>L</I> contain the amino acids. Amino
 * acids <TT>'A'</TT> through <TT>'Z'</TT> (case insensitive) are represented by
 * the values 0 through 25; <TT>'*'</TT> is represented as 26; <TT>'-'</TT> is
 * represented as 27.
 * <P>
 * The amino acid letters and values are:
 * <TABLE BORDER=0 CELLPADDING=0 CELLSPACING=0>
 * <TR><TD><I>Letter</I>&nbsp;&nbsp;</I></TD><TD><I>Value</I>&nbsp;&nbsp;</TD><TD><I>Amino Acid</I></TD></TR>
 * <TR><TD>A</TD><TD>0</TD><TD>Alanine</TD></TR>
 * <TR><TD>B</TD><TD>1</TD><TD>Aspartate or asparagine</TD></TR>
 * <TR><TD>C</TD><TD>2</TD><TD>Cysteine</TD></TR>
 * <TR><TD>D</TD><TD>3</TD><TD>Aspartate</TD></TR>
 * <TR><TD>E</TD><TD>4</TD><TD>Glutamate</TD></TR>
 * <TR><TD>F</TD><TD>5</TD><TD>Phenylalanine</TD></TR>
 * <TR><TD>G</TD><TD>6</TD><TD>Glycine</TD></TR>
 * <TR><TD>H</TD><TD>7</TD><TD>Histidine</TD></TR>
 * <TR><TD>I</TD><TD>8</TD><TD>Isoleucine</TD></TR>
 * <TR><TD>J</TD><TD>9</TD><TD><I>unused</I></TD></TR>
 * <TR><TD>K</TD><TD>10</TD><TD>Lysine</TD></TR>
 * <TR><TD>L</TD><TD>11</TD><TD>Leucine</TD></TR>
 * <TR><TD>M</TD><TD>12</TD><TD>Methionine</TD></TR>
 * <TR><TD>N</TD><TD>13</TD><TD>Asparagine</TD></TR>
 * <TR><TD>O</TD><TD>14</TD><TD><I>unused</I></TD></TR>
 * <TR><TD>P</TD><TD>15</TD><TD>Proline</TD></TR>
 * <TR><TD>Q</TD><TD>16</TD><TD>Glutamine</TD></TR>
 * <TR><TD>R</TD><TD>17</TD><TD>Arginine</TD></TR>
 * <TR><TD>S</TD><TD>18</TD><TD>Serine</TD></TR>
 * <TR><TD>T</TD><TD>19</TD><TD>Threonine</TD></TR>
 * <TR><TD>U</TD><TD>20</TD><TD>Selenocysteine</TD></TR>
 * <TR><TD>V</TD><TD>21</TD><TD>Valine</TD></TR>
 * <TR><TD>W</TD><TD>22</TD><TD>Tryptophan</TD></TR>
 * <TR><TD>X</TD><TD>23</TD><TD>Any, unknown</TD></TR>
 * <TR><TD>Y</TD><TD>24</TD><TD>Tyrosine</TD></TR>
 * <TR><TD>Z</TD><TD>25</TD><TD>Glutamate or glutamine</TD></TR>
 * <TR><TD>*</TD><TD>26</TD><TD>Translation stop</TD></TR>
 * <TR><TD>-</TD><TD>27</TD><TD>Gap of indeterminate length</TD></TR>
 * </TABLE>
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class ProteinSequence
	extends Sequence
	{

// Exported constructors.

	/**
	 * Construct a new protein sequence from the given string.
	 *
	 * @param  description
	 *     Description string. Must start with a <TT>'&gt;'</TT> character.
	 * @param  sequence
	 *     Sequence string. Must consist of the characters <TT>'A'</TT> through
	 *     <TT>'Z'</TT>, <TT>'a'</TT> through <TT>'z'</TT>, <TT>'*'</TT>, and
	 *     <TT>'-'</TT>.
	 */
	public ProteinSequence
		(String description,
		 String sequence)
		{
		// Make sure description starts with '>'.
		if (description.charAt(0) != '>')
			{
			throw new IllegalArgumentException
				("ProteinSequence(): Invalid description");
			}
		myDescription = description;

		// Read characters of the protein sequence.
		myLength = sequence.length();
		mySequence = new byte [myLength + 1];
		mySequence[0] = (byte)(-1);
		for (int i = 0; i < myLength; ++ i)
			{
			char b = sequence.charAt(i);
			if ('A' <= b && b <= 'Z')
				{
				mySequence[i+1] = (byte)(b - 'A');
				}
			else if ('a' <= b && b <= 'z')
				{
				mySequence[i+1] = (byte)(b - 'a');
				}
			else if (b == '*')
				{
				mySequence[i+1] = (byte)(26);
				}
			else if (b == '-')
				{
				mySequence[i+1] = (byte)(27);
				}
			else
				{
				throw new IllegalArgumentException
					("ProteinSequence(): Invalid amino acid '"+b+"'");
				}
			}
		}

	/**
	 * Construct a new protein sequence read from the given file.
	 *
	 * @param  file  File.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public ProteinSequence
		(File file)
		throws IOException
		{
		InputStream in = null;
		try
			{
			in = new BufferedInputStream (new FileInputStream (file));
			read (in);
			}
		finally
			{
			if (in != null)
				{
				try { in.close(); } catch (IOException exc) {}
				}
			}
		}

	/**
	 * Construct a new protein sequence read from the given input stream.
	 *
	 * @param  in  Input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	ProteinSequence
		(InputStream in)
		throws IOException
		{
		read (in);
		}

	/**
	 * Read this protein sequence from the given input stream.
	 *
	 * @param  in  Input stream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void read
		(InputStream in)
		throws IOException
		{
		// Make sure description starts with '>'.
		int b = in.read();
		if (b != '>')
			{
			throw new IOException
				("ProteinSequence(): Invalid description line");
			}

		// Read bytes of the description up until end of line or end of file.
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		baos.write (b);
		while ((b = in.read()) != -1 && b != '\r' && b != '\n') baos.write (b);
		myDescription = new String (baos.toByteArray());

		// Read bytes of the protein sequence.
		baos.reset();
		baos.write (-1);
		for (;;)
			{
			b = in.read();
			if (b == -1 || b == '>')
				{
				break;
				}
			else if ('A' <= b && b <= 'Z')
				{
				baos.write (b - 'A');
				}
			else if ('a' <= b && b <= 'z')
				{
				baos.write (b - 'a');
				}
			else if (b == '*')
				{
				baos.write (26);
				}
			else if (b == '-')
				{
				baos.write (27);
				}
			else if (b == '\r' || b == '\n')
				{
				}
			else
				{
				throw new IOException
					("ProteinSequence(): Invalid amino acid '"+((char) b)+"'");
				}
			}
		mySequence = baos.toByteArray();
		myLength = mySequence.length - 1;
		}

// Exported operations.

	/**
	 * Returns a character version of this protein sequence's element at the
	 * given index.
	 *
	 * @param  i  Index in the range 1 .. <I>L</I>.
	 *
	 * @return  Character corresponding to element <TT>i</TT>.
	 */
	public char charAt
		(int i)
		{
		if (1 > i || i > myLength)
			{
			throw new IndexOutOfBoundsException
				("ProteinSequence.charAt(): Index "+i+" out of bounds");
			}
		int aa = mySequence[i];
		if (0 <= aa && aa <= 25) return (char)('A' + aa);
		else if (aa == 26) return '*';
		else return '-';
		}

	/**
	 * Returns a string version of this protein sequence. The string is
	 * <TT>"ProteinSequence(<I>description</I>)"</TT>.
	 *
	 * @return  String version.
	 */
	public String toString()
		{
		return "ProteinSequence(" + myDescription + ")";
		}

	}
