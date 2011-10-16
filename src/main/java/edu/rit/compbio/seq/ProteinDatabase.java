//******************************************************************************
//
// File:    ProteinDatabase.java
// Package: edu.rit.compbio.seq
// Unit:    Class edu.rit.compbio.seq.ProteinDatabase
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
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;

/**
 * Class ProteinDatabase encapsulates a protein sequence database stored in a
 * file. Class ProteinDatabase is used to read {@linkplain ProteinSequence}
 * objects from the database.
 * <P>
 * A protein sequence database file consists of one or more protein sequences
 * stored in FASTA format. A FASTA format protein sequence consists of one
 * description line and one or more sequence lines. The description line
 * consists of an initial <TT>'&gt;'</TT> character followed by zero or more
 * characters (the protein's description). A sequence line consists of one or
 * more characters <TT>'A'</TT> through <TT>'Z'</TT>, <TT>'a'</TT> through
 * <TT>'z'</TT>, <TT>'*'</TT>, or <TT>'-'</TT>. For further information, see
 * class {@linkplain ProteinSequence}.
 * <P>
 * Along with the protein sequence database file, there is a protein sequence
 * index file. The index file tells where each protein is located in the
 * database file. Class ProteinDatabase's constructor requires both the database
 * file and the index file to be supplied as arguments.
 * <P>
 * Class ProteinDatabase includes a main program that reads a protein sequence
 * database file and creates the requisite protein sequence index file. If the
 * <I>n</I> argument is specified, the program creates an index file for just
 * the first <I>n</I> sequences in the database. If the <I>n</I> argument is
 * omitted, the program creates an index file for all the sequences in the
 * database.
 * <P>
 * Usage: java edu.rit.compbio.seq.ProteinDatabase <I>databasefile</I>
 * <I>indexfile</I> [ <I>n</I> ]
 * <BR><I>databasefile</I> = Input protein sequence database file
 * <BR><I>indexfile</I> = Output protein sequence index file
 * <BR><I>n</I> = Number of sequences in the index (default: all)
 *
 * @author  Alan Kaminsky
 * @version 01-Jul-2008
 */
public class ProteinDatabase
	{

// Hidden data members.

	private File myDatabaseFile;
	private File myIndexFile;
	private RandomAccessFile myIndex;
	private long myProteinCount;
	private long myDatabaseLength;
	private long myFileLength;

// Exported constructors.

	/**
	 * Construct a new protein sequence database.
	 *
	 * @param  theDatabaseFile  Protein sequence database file.
	 * @param  theIndexFile     Protein sequence index file.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDatabaseFile</TT> is null.
	 *     Thrown if <TT>theIndexFile</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public ProteinDatabase
		(File theDatabaseFile,
		 File theIndexFile)
		throws IOException
		{
		// Verify preconditions.
		if (theDatabaseFile == null)
			{
			throw new NullPointerException
				("ProteinDatabase(): theDatabaseFile is null");
			}
		if (theIndexFile == null)
			{
			throw new NullPointerException
				("ProteinDatabase(): theIndexFile is null");
			}

		// Record files.
		myDatabaseFile = theDatabaseFile;
		myIndexFile = theIndexFile;

		// Open index file.
		myIndex = new RandomAccessFile (theIndexFile, "r");

		// Determine number of entries in index file.
		long n = myIndex.length();
		if ((n & 0x7L) != 0L)
			{
			throw new IOException
				("ProteinDatabase(): Index file \""+theIndexFile+
				 "\" has invalid length (= "+n+")");
			}
		myProteinCount = (n >> 3) - 1;

		// Determine total length of all sequences in database.
		myDatabaseLength = myIndex.readLong();

		// Determine length of database file.
		myFileLength = theDatabaseFile.length();
		}

// Exported operations.

	/**
	 * Get the sum of the lengths of the protein sequences in this protein
	 * sequence database.
	 *
	 * @return  Total protein sequence length.
	 */
	public long getDatabaseLength()
		{
		return myDatabaseLength;
		}

	/**
	 * Get the number of protein sequences in this protein sequence database.
	 *
	 * @return  Number of protein sequences, <I>N</I>.
	 */
	public long getProteinCount()
		{
		return myProteinCount;
		}

	/**
	 * Get the protein sequence at the given index in this protein sequence
	 * database.
	 *
	 * @param  i  Index in the range 0 &le; <TT>i</TT> &le; <I>N</I>&minus;1.
	 *
	 * @return  Protein sequence.
	 *
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public ProteinSequence getProteinSequence
		(long i)
		throws IOException
		{
		// Verify preconditions.
		if (0 > i || i >= myProteinCount)
			{
			throw new IndexOutOfBoundsException
				("ProteinDatabase.getProteinSequence(): i (= "+i+
				 ") out of bounds");
			}

		// Get offset to protein sequence from index file, in a critical section
		// for multiple thread safety.
		long offset;
		synchronized (this)
			{
			myIndex.seek ((i + 1) << 3);
			offset = myIndex.readLong();
			}
		if (0 > offset || offset >= myFileLength)
			{
			throw new IOException
				("ProteinDatabase.getProteinSequence("+i+
				 "): Invalid offset (= "+offset+")");
			}

		// Get protein sequence from database file.
		InputStream fis = null;
		try
			{
			// Open database file.
			fis =
				(new BufferedInputStream
					(new FileInputStream (myDatabaseFile)));

			// Seek to offset.
			long remaining = offset;
			long skipped = 0;
			while (remaining > 0)
				{
				skipped = fis.skip (remaining);
				if (skipped == 0)
					{
					throw new IOException
						("ProteinDatabase.getProteinSequence("+i+
						 "): Unexpected end of file");
					}
				remaining -= skipped;
				}

			// Read protein sequence.
			return new ProteinSequence (fis);
			}
		finally
			{
			// Close database file.
			if (fis != null)
				{
				try { fis.close(); } catch (IOException exc) {}
				}
			}
		}

	/**
	 * Close this protein sequence database.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void close()
		throws IOException
		{
		myIndex.close();
		}

	/**
	 * Finalize this protein sequence database.
	 */
	protected void finalize()
		{
		try { close(); } catch (IOException exc) {}
		}

	/**
	 * Main program that reads a protein sequence database file and creates the
	 * requisite protein sequence index file. If the <I>n</I> argument is
	 * specified, the program creates an index file for just the first <I>n</I>
	 * sequences in the database. If the <I>n</I> argument is omitted, the
	 * program creates an index file for all the sequences in the database.
	 * <P>
	 * Usage: java edu.rit.compbio.seq.ProteinDatabase <I>databasefile</I>
	 * <I>indexfile</I> [ <I>n</I> ]
	 * <BR><I>databasefile</I> = Input protein sequence database file
	 * <BR><I>indexfile</I> = Output protein sequence index file
	 * <BR><I>n</I> = Number of sequences in the index (default: all)
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Parse command line arguments.
		if (args.length < 2 || args.length > 3) usage();
		File databasefile = new File (args[0]);
		File indexfile = new File (args[1]);
		long n = Long.MAX_VALUE;
		if (args.length == 3) n = Long.parseLong (args[2]);

		// Open database file for reading.
		InputStream databasein =
			new BufferedInputStream
				(new FileInputStream (databasefile));

		// Open index file for writing. Set aside space for database length.
		RandomAccessFile indexout = new RandomAccessFile (indexfile, "rw");
		indexout.writeLong (0L);

		// Scan the database file. Whenever a '>' is encountered at the
		// beginning of a line, store its offset in the index file.
		long offset = 0L;
		long dblength = 0L;
		int state = 0;
		int b;
		readloop: while ((b = databasein.read()) != -1)
			{
			switch (state)
				{
				case 0: // Beginning of line
					if (b == '>')
						{
						if (n == 0) break readloop;
						indexout.writeLong (offset);
						-- n;
						state = 1;
						}
					else if (b == '\r' || b == '\n')
						{
						state = 0;
						}
					else
						{
						++ dblength;
						state = 2;
						}
					break;
				case 1: // In a description line
					if (b == '\r' || b == '\n')
						{
						state = 0;
						}
					else
						{
						state = 1;
						}
					break;
				case 2: // In a sequence line
					if (b == '\r' || b == '\n')
						{
						state = 0;
						}
					else
						{
						++ dblength;
						state = 2;
						}
					break;
				}
			++ offset;
			}

		// Write database length.
		indexout.seek (0L);
		indexout.writeLong (dblength);

		// All done.
		databasein.close();
		indexout.close();
		}

// Hidden operations.

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.compbio.seq.ProteinDatabase <databasefile> <indexfile> [<length>]");
		System.err.println ("<databasefile> = Input protein sequence database file");
		System.err.println ("<indexfile> = Output protein sequence index file");
		System.err.println ("<length> = Number of sequences in the index (default: all)");
		System.exit (1);
		}

	}
