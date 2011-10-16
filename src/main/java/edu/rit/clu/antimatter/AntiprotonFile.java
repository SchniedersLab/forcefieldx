//******************************************************************************
//
// File:    AntiprotonFile.java
// Package: edu.rit.clu.antimatter
// Unit:    Class edu.rit.clu.antimatter.AntiprotonFile
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

package edu.rit.clu.antimatter;

import edu.rit.vector.Vector2D;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Class AntiprotonFile provides an object for reading or writing a series of
 * snapshots of antiproton positions from or to a file.
 * <P>
 * To write the snapshots to a file, call the <TT>prepareToWrite()</TT> method,
 * specifying the output stream to write. The <TT>prepareToWrite()</TT> method
 * returns an instance of class {@linkplain AntiprotonFile.Writer}. Call the
 * methods of the writer object to write the antiproton position snapshots to
 * the output stream. When finished, close the writer.
 * <P>
 * To read the snapshots from a file, call the <TT>prepareToRead()</TT> method,
 * specifying the input stream to read. The <TT>prepareToRead()</TT> method
 * returns an instance of class {@linkplain AntiprotonFile.Reader}. Call the
 * methods of the reader object to read the antiproton position snapshots from
 * the input stream. When finished, close the reader.
 * <P>
 * Class AntiprotonFile includes a main program to combine a group of antiproton
 * files into one antiproton file. You might use this main program when the
 * processes of a parallel program have computed slices of the antiprotons and
 * stored the slices in separate files, and you want to combine the slices
 * together into one file.
 * <P>
 * <B>Antiproton File Format</B>
 * <P>
 * An antiproton file is a binary file containing the following items. Each
 * primitive item is written as though by java.io.DataOutput (<TT>int</TT>, four
 * bytes; <TT>long</TT>, eight bytes; <TT>double</TT>, eight bytes; most
 * significant byte first).
 * <UL>
 * <LI>
 * Random seed used to generate the initial antiproton positions, <I>seed</I>
 * (<TT>long</TT>).
 * <LI>
 * Side of the square used to generate the initial antiproton positions,
 * <I>R</I> (<TT>double</TT>).
 * <LI>
 * Time step size, <I>dt</I> (<TT>double</TT>).
 * <LI>
 * Number of time steps between snapshots, <I>steps</I> (<TT>int</TT>).
 * <LI>
 * Number of snapshots, <I>snaps</I> (<TT>int</TT>).
 * <LI>
 * Total number of antiprotons, <I>N</I> (<TT>int</TT>).
 * <LI>
 * Lower antiproton index for each snapshot, <I>L</I> (<TT>int</TT>). <I>L</I>
 * &gt;= 0.
 * <LI>
 * Number of antiprotons in each snapshot, <I>M</I> (<TT>int</TT>). <I>M</I>
 * &gt;= 0. <I>L</I>+<I>M</I> &lt;= <I>N</I>.
 * <LI>
 * The snapshots.
 * </UL>
 * <P>
 * Each snapshot contains the following items. A snapshot contains information
 * about a slice of the antiprotons.
 * <UL>
 * <LI>
 * The antiproton positions for antiproton indexes <I>L</I> ..
 * <I>L</I>+<I>M</I>-1. Each position is stored as the X coordinate
 * (<TT>double</TT>) followed by the Y coordinate (<TT>double</TT>).
 * <LI>
 * The total momentum for the antiprotons at antiproton indexes <I>L</I> ..
 * <I>L</I>+<I>M</I>-1. Stored as the X coordinate (<TT>double</TT>) followed by
 * the Y coordinate (<TT>double</TT>).
 * </UL>
 *
 * @author  Alan Kaminsky
 * @version 09-Feb-2008
 */
public class AntiprotonFile
	{

// Hidden data members.

	// Attributes.
	private long seed;
	private double R;
	private double dt;
	private int steps;
	private int snaps;
	private int N;
	private int L;
	private int M;

	// True if attributes are initialized, false otherwise.
	private boolean iamInitialized;

// Exported constructors.

	/**
	 * Construct a new antiproton file object. The antiproton file's attributes
	 * are uninitialized. Before using the antiproton file, specify the
	 * attributes by calling the <TT>setAttributes()</TT> method or by reading
	 * the antiproton file from an input stream.
	 */
	public AntiprotonFile()
		{
		}

	/**
	 * Construct a new antiproton file object with the given attributes.
	 *
	 * @param  seed   Random seed used to generate the initial antiproton
	 *                positions.
	 * @param  R      Side of the square used to generate the initial antiproton
	 *                positions.
	 * @param  dt     Time step size.
	 * @param  steps  Number of time steps between snapshots.
	 * @param  snaps  Number of snapshots.
	 * @param  N      Total number of antiprotons.
	 * @param  L      Lower antiproton index for each snapshot.
	 * @param  M      Number of antiprotons in each snapshot.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>R</TT> &lt;= 0. Thrown if
	 *     <TT>dt</TT> &lt;= 0. Thrown if <TT>steps</TT> &lt;= 0. Thrown if
	 *     <TT>snaps</TT> &lt; 0. Thrown if <TT>N</TT> &lt; 0.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>L</TT> &lt; 0. Thrown if
	 *     <TT>M</TT> &lt; 0. Thrown if <TT>L+M</TT> &gt; <TT>N</TT>.
	 */
	public AntiprotonFile
		(long seed,
		 double R,
		 double dt,
		 int steps,
		 int snaps,
		 int N,
		 int L,
		 int M)
		{
		setAttributes (seed, R, dt, steps, snaps, N, L, M);
		}

// Exported operations.

	/**
	 * Get the random seed used to generate the initial antiproton positions.
	 *
	 * @return  Random seed, <I>seed</I>.
	 */
	public long getSeed()
		{
		verifyInitialized();
		return this.seed;
		}

	/**
	 * Get the side of the square used to generate the initial antiproton
	 * positions.
	 *
	 * @return  Side, <I>R</I>.
	 */
	public double getR()
		{
		verifyInitialized();
		return this.R;
		}

	/**
	 * Get the time step size.
	 *
	 * @return  Time step size, <I>dt</I>.
	 */
	public double getDt()
		{
		verifyInitialized();
		return this.dt;
		}

	/**
	 * Get the number of time steps between snapshots.
	 *
	 * @return  Number of time steps, <I>steps</I>.
	 */
	public int getSteps()
		{
		verifyInitialized();
		return this.steps;
		}

	/**
	 * Get the number of snapshots.
	 *
	 * @return  Number of snapshots, <I>snaps</I>.
	 */
	public int getSnaps()
		{
		verifyInitialized();
		return this.snaps;
		}

	/**
	 * Get the total number of antiprotons.
	 *
	 * @return  Total number of antiprotons, <I>N</I>.
	 */
	public int getN()
		{
		verifyInitialized();
		return this.N;
		}

	/**
	 * Get the lower antiproton index for each snapshot.
	 *
	 * @return  Lower antiproton index for each snapshot, <I>L</I>.
	 */
	public int getL()
		{
		verifyInitialized();
		return this.L;
		}

	/**
	 * Get the number of antiprotons in each snapshot.
	 *
	 * @return  Number of antiprotons in each snapshot, <I>M</I>.
	 */
	public int getM()
		{
		verifyInitialized();
		return this.M;
		}

	/**
	 * Set this antiproton file's attributes.
	 *
	 * @param  seed   Random seed used to generate the initial antiproton
	 *                positions.
	 * @param  R      Side of the square used to generate the initial antiproton
	 *                positions.
	 * @param  dt     Time step size.
	 * @param  steps  Number of time steps between snapshots.
	 * @param  snaps  Number of snapshots.
	 * @param  N      Total number of antiprotons.
	 * @param  L      Lower antiproton index for each snapshot.
	 * @param  M      Number of antiprotons in each snapshot.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>R</TT> &lt;= 0. Thrown if
	 *     <TT>dt</TT> &lt;= 0. Thrown if <TT>steps</TT> &lt;= 0. Thrown if
	 *     <TT>snaps</TT> &lt; 0. Thrown if <TT>N</TT> &lt; 0.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>L</TT> &lt; 0. Thrown if
	 *     <TT>M</TT> &lt; 0. Thrown if <TT>L+M</TT> &gt; <TT>N</TT>.
	 */
	public void setAttributes
		(long seed,
		 double R,
		 double dt,
		 int steps,
		 int snaps,
		 int N,
		 int L,
		 int M)
		{
		if (R <= 0.0)
			{
			throw new IllegalArgumentException
				("AntiprotonFile.setAttributes(): R <= 0");
			}
		if (dt <= 0.0)
			{
			throw new IllegalArgumentException
				("AntiprotonFile.setAttributes(): dt <= 0");
			}
		if (steps < 0)
			{
			throw new IllegalArgumentException
				("AntiprotonFile.setAttributes(): steps < 0");
			}
		if (snaps < 0)
			{
			throw new IllegalArgumentException
				("AntiprotonFile.setAttributes(): snaps < 0");
			}
		if (snaps < 0)
			{
			throw new IllegalArgumentException
				("AntiprotonFile.setAttributes(): N < 0");
			}
		if (L < 0)
			{
			throw new IndexOutOfBoundsException
				("AntiprotonFile.setAttributes(): L < 0");
			}
		if (M < 0)
			{
			throw new IndexOutOfBoundsException
				("AntiprotonFile.setAttributes(): M < 0");
			}
		if (L+M > N)
			{
			throw new IndexOutOfBoundsException
				("AntiprotonFile.setAttributes(): L+M > N");
			}

		this.seed = seed;
		this.R = R;
		this.dt = dt;
		this.steps = steps;
		this.snaps = snaps;
		this.N = N;
		this.L = L;
		this.M = M;

		iamInitialized = true;
		}

	/**
	 * Prepare to write this antiproton file to the given output stream. This
	 * matrix file's attributes are written to the output stream at this time.
	 * To write this antiproton file's snapshots, call methods on the returned
	 * writer, then close the writer.
	 * <P>
	 * For improved performance, specify an output stream with buffering, such
	 * as an instance of class java.io.BufferedOutputStream.
	 *
	 * @param  theStream  Output stream.
	 *
	 * @return  Writer object with which to write this antiproton file.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if this antiproton file object is
	 *     uninitialized.
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Writer prepareToWrite
		(OutputStream theStream)
		throws IOException
		{
		verifyInitialized();
		return new Writer (theStream);
		}

	/**
	 * Prepare to read this antiproton file from the given input stream. This
	 * antiproton file's attributes are read from the input stream at this time.
	 * To read this antiproton file's snapshots, call methods on the returned
	 * reader, then close the reader.
	 * <P>
	 * For improved performance, specify an input stream with buffering, such
	 * as an instance of class java.io.BufferedInputStream.
	 *
	 * @param  theStream  Input stream.
	 *
	 * @return  Reader object with which to read this antiproton file.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public Reader prepareToRead
		(InputStream theStream)
		throws IOException
		{
		return new Reader (theStream);
		}

	/**
	 * Main program to combine a group of antiproton files into one antiproton
	 * file. Each input file must be for an antiproton array with the same
	 * number of elements. Each input file must contain the same number of
	 * snapshots.
	 * <P>
	 * You might use this main program when the processes of a parallel program
	 * have computed slices of antiprotons and stored the slices in separate
	 * files, and you want to combine the slices together into one file.
	 * <P>
	 * Usage: java edu.rit.io.AntiprotonFile <I>outfile</I> <I>infile1</I>
	 * [ <I>infile2</I> . . . ]
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		// Validate command line arguments.
		if (args.length < 2)
			{
			System.err.println ("Usage: java edu.rit.io.AntiprotonFile <outfile> <infile1> [<infile2> ...]");
			System.exit (1);
			}

		// Set up input antiproton file objects.
		int size = args.length - 1;
		AntiprotonFile[] infile = new AntiprotonFile [size];
		AntiprotonFile.Reader[] reader = new AntiprotonFile.Reader [size];
		for (int i = 0; i < size; ++ i)
			{
			infile[i] = new AntiprotonFile();
			reader[i] = infile[i].prepareToRead
				(new BufferedInputStream
					(new FileInputStream (args[i+1])));
			}

		// Verify that all input files have the same attributes.
		long seed = infile[0].getSeed();
		double R = infile[0].getR();
		double dt = infile[0].getDt();
		int steps = infile[0].getSteps();
		int snaps = infile[0].getSnaps();
		int N = infile[0].getN();
		for (int i = 1; i < size; ++ i)
			{
			if (seed != infile[i].getSeed())
				{
				System.err.println
					("AntiprotonFile: " + args[1] + " seed (" + seed +
					 ") != " + args[i+1] + " seed (" + infile[i].getSeed());
				System.exit (1);
				}
			if (R != infile[i].getR())
				{
				System.err.println
					("AntiprotonFile: " + args[1] + " R (" + R +
					 ") != " + args[i+1] + " R (" + infile[i].getR());
				System.exit (1);
				}
			if (dt != infile[i].getDt())
				{
				System.err.println
					("AntiprotonFile: " + args[1] + " dt (" + dt +
					 ") != " + args[i+1] + " dt (" + infile[i].getDt());
				System.exit (1);
				}
			if (steps != infile[i].getSteps())
				{
				System.err.println
					("AntiprotonFile: " + args[1] + " steps (" + steps +
					 ") != " + args[i+1] + " steps (" + infile[i].getSteps());
				System.exit (1);
				}
			if (snaps != infile[i].getSnaps())
				{
				System.err.println
					("AntiprotonFile: " + args[1] + " snaps (" + snaps +
					 ") != " + args[i+1] + " snaps (" + infile[i].getSnaps());
				System.exit (1);
				}
			if (N != infile[i].getN())
				{
				System.err.println
					("AntiprotonFile: " + args[1] + " N (" + N +
					 ") != " + args[i+1] + " N (" + infile[i].getN());
				System.exit (1);
				}
			}

		// Set up position array and total momentum.
		Vector2D[] p = new Vector2D [N];
		for (int i = 0; i < N; ++ i)
			{
			p[i] = new Vector2D();
			}
		Vector2D totalmv = new Vector2D();
		Vector2D mv = new Vector2D();

		// Set up output file.
		AntiprotonFile outfile =
			new AntiprotonFile (seed, R, dt, steps, snaps, N, 0, N);
		AntiprotonFile.Writer writer =
			outfile.prepareToWrite
				(new BufferedOutputStream
					(new FileOutputStream (args[0])));

		// Do all snapshots.
		for (int s = 0; s < snaps; ++ s)
			{
			totalmv.clear();

			// Read snapshot from all input files.
			for (int i = 0; i < size; ++ i)
				{
				reader[i].readSnapshot (p, infile[i].getL(), mv);
				totalmv.add (mv);
				}

			// Write snapshot to output file.
			writer.writeSnapshot (p, 0, totalmv);
			}

		// Close all files.
		for (int i = 0; i < size; ++ i)
			{
			reader[i].close();
			}
		writer.close();
		}

// Exported helper classes.

	/**
	 * Class AntiprotonFile.Writer provides an object with which to write an
	 * {@linkplain AntiprotonFile} to an output stream.
	 * <P>
	 * When a writer is created, the antiproton file's attributes are written to
	 * the output stream. Each time the <TT>writeSnapshot()</TT> method is
	 * called, one snapshot is written to the output stream. When finished, call
	 * the <TT>close()</TT> method.
	 * <P>
	 * <I>Note:</I> Class AntiprotonFile.Writer is not multiple thread safe.
	 *
	 * @author  Alan Kaminsky
	 * @version 05-Feb-2008
	 */
	public class Writer
		{

	// Hidden data members.

		private OutputStream myOs;
		private DataOutputStream myDos;
		private int mySnaps;

	// Hidden constructors.

		/**
		 * Construct a new antiproton file writer.
		 *
		 * @param  theStream  Output stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private Writer
			(OutputStream theStream)
			throws IOException
			{
			if (theStream == null)
				{
				throw new NullPointerException
					("AntiprotonFile.Writer(): theStream is null");
				}
			myOs = theStream;
			myDos = new DataOutputStream (theStream);
			myDos.writeLong (seed);
			myDos.writeDouble (R);
			myDos.writeDouble (dt);
			myDos.writeInt (steps);
			myDos.writeInt (snaps);
			myDos.writeInt (N);
			myDos.writeInt (L);
			myDos.writeInt (M);
			}

	// Exported operations.

		/**
		 * Write a snapshot to the output stream.
		 * <P>
		 * The <I>M</I> antiproton positions starting at index <TT>off</TT> in
		 * the array <TT>p</TT> are written. Note that, depending on how the
		 * antiproton position array <TT>p</TT> has been allocated, the argument
		 * <TT>off</TT> need not be the same as the attribute <I>L</I>. However,
		 * the element <TT>p[off]</TT> must contain the position of antiproton
		 * <I>L</I>, the element <TT>p[off+1]</TT> must contain the position of
		 * antiproton <I>L</I>+1, and so on.
		 * <P>
		 * The given total momentum is also written; this must be the total
		 * momentum of the antiprotons <I>L</I> .. <I>L</I>+<I>M</I>-1.
		 *
		 * @param  p    Array of antiproton positions.
		 * @param  off  Index of first array element to write.
		 * @param  mv   Total momentum.
		 *
		 * @exception  IllegalStateException
		 *     (unchecked exception) Thrown if <I>snaps</I> snapshots have
		 *     already been written.
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>p</TT> is null. Thrown if
		 *     <TT>mv</TT> is null.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if <TT>off</TT> &lt; 0. Thrown if
		 *     <TT>off</TT>+<I>M</I> &gt; <TT>p.length</TT>.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void writeSnapshot
			(Vector2D[] p,
			 int off,
			 Vector2D mv)
			throws IOException
			{
			if (mySnaps >= snaps)
				{
				throw new IllegalStateException
					("AntiprotonFile.Writer.writeSnapshot(): Too many snapshots");
				}
			if (off < 0 || off+M > p.length)
				{
				throw new IndexOutOfBoundsException();
				}
			for (int i = 0; i < M; ++ i)
				{
				myDos.writeDouble (p[i+off].x);
				myDos.writeDouble (p[i+off].y);
				}
			myDos.writeDouble (mv.x);
			myDos.writeDouble (mv.y);
			++ mySnaps;
			}

		/**
		 * Close the output stream.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void close()
			throws IOException
			{
			myDos.close();
			}
		}

	/**
	 * Class AntiprotonFile.Reader provides an object with which to read an
	 * {@linkplain AntiprotonFile} from an input stream.
	 * <P>
	 * When a reader is created, the antiproton file's attributes are read from
	 * the input stream. Each time the <TT>readSnapshot()</TT> method is
	 * called, one snapshot is read from the input stream. When finished, call
	 * the <TT>close()</TT> method.
	 * <P>
	 * <I>Note:</I> Class AntiprotonFile.Reader is not multiple thread safe.
	 *
	 * @author  Alan Kaminsky
	 * @version 05-Feb-2008
	 */
	public class Reader
		{

	// Hidden data members.

		private InputStream myIs;
		private DataInputStream myDis;
		private int mySnaps;

	// Hidden constructors.

		/**
		 * Construct a new antiproton file reader.
		 *
		 * @param  theStream  Input stream.
		 *
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>theStream</TT> is null.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		private Reader
			(InputStream theStream)
			throws IOException
			{
			if (theStream == null)
				{
				throw new NullPointerException
					("AntiprotonFile.Reader(): theStream is null");
				}
			myIs = theStream;
			myDis = new DataInputStream (theStream);
			setAttributes
				(/*seed */ myDis.readLong(),
				 /*R    */ myDis.readDouble(),
				 /*dt   */ myDis.readDouble(),
				 /*steps*/ myDis.readInt(),
				 /*snaps*/ myDis.readInt(),
				 /*N    */ myDis.readInt(),
				 /*L    */ myDis.readInt(),
				 /*M    */ myDis.readInt());
			}

	// Exported operations.

		/**
		 * Read a snapshot from the input stream.
		 * <P>
		 * The <I>M</I> antiproton positions are stored starting at index
		 * <TT>off</TT> in the array <TT>p</TT>. Note that, depending on how the
		 * antiproton position array <TT>p</TT> has been allocated, the argument
		 * <TT>off</TT> need not be the same as the attribute <I>L</I>. However,
		 * the element <TT>p[off]</TT> must contain the position of antiproton
		 * <I>L</I>, the element <TT>p[off+1]</TT> must contain the position of
		 * antiproton <I>L</I>+1, and so on.
		 * <P>
		 * The total momentum of the antiprotons <I>L</I> .. <I>L</I>+<I>M</I>-1
		 * is also read and is stored in <TT>mv</TT>.
		 *
		 * @param  p    Array of antiproton positions.
		 * @param  off  Index of first array element to store.
		 * @param  mv   Total momentum (output).
		 *
		 * @return  The number of the just-read snapshot. Snapshot numbers start
		 *          at 0.
		 *
		 * @exception  IllegalStateException
		 *     (unchecked exception) Thrown if <I>snaps</I> snapshots have
		 *     already been read.
		 * @exception  NullPointerException
		 *     (unchecked exception) Thrown if <TT>p</TT> is null. Thrown if
		 *     <TT>mv</TT> is null.
		 * @exception  IndexOutOfBoundsException
		 *     (unchecked exception) Thrown if <TT>off</TT> &lt; 0. Thrown if
		 *     <TT>off</TT>+<I>M</I> &gt; <TT>p.length</TT>.
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public int readSnapshot
			(Vector2D[] p,
			 int off,
			 Vector2D mv)
			throws IOException
			{
			if (mySnaps >= snaps)
				{
				throw new IllegalStateException
					("AntiprotonFile.Reader.readSnapshot(): Too many snapshots");
				}
			if (off < 0 || off+M > p.length)
				{
				throw new IndexOutOfBoundsException();
				}
			for (int i = 0; i < M; ++ i)
				{
				p[i+off].x = myDis.readDouble();
				p[i+off].y = myDis.readDouble();
				}
			mv.x = myDis.readDouble();
			mv.y = myDis.readDouble();
			return mySnaps ++;
			}

		/**
		 * Close the input stream.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void close()
			throws IOException
			{
			myDis.close();
			}
		}

// Hidden operations.

	/**
	 * Verify whether this antiproton file is initialized.
	 *
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if this antiproton file is not
	 *     initialized.
	 */
	private void verifyInitialized()
		{
		if (! iamInitialized)
			{
			throw new IllegalStateException
				("AntiprotonFile: Not initialized");
			}
		}
		
	}
