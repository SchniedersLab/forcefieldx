//******************************************************************************
//
// File:    CreateSignalDataSet.java
// Package: edu.rit.mri
// Unit:    Class edu.rit.mri.CreateSignalDataSet
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

package edu.rit.mri;

import edu.rit.numeric.ListSeries;

import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import java.nio.channels.FileChannel;

import java.util.ArrayList;
import java.util.Scanner;

/**
 * Class CreateSignalDataSet is a preprocessing program that reads a magnetic
 * resonance image stored in a group of files in a certain format and creates a
 * spin signal data set file. This spin signal data set file is then used by the
 * other programs in package edu.rit.mri.
 * <P>
 * Usage: java edu.rit.mri.CreateSignalDataSet <I>timesfile</I> <I>maskfile</I>
 * <I>signalfile</I>
 * <BR><I>timesfile</I> = Input times file
 * <BR><I>maskfile</I> = Input mask file
 * <BR><I>signalfile</I> = Output spin signal data set file
 * <P>
 * The program's error handling is rudimentary. If an error is encountered, the
 * program throws an exception and terminates.
 * <P>
 * The format of the output spin signal data set file is described in class
 * {@linkplain SignalDataSetWriter}. The formats of the input files are
 * described below.
 * <P>
 * <B>Times file.</B> The times file is a plain text file giving the image
 * dimensions in pixels and the list of discrete time values at which the images
 * were taken. Each line of the times file consists of fields separated by white
 * space. The first line of the times file has two fields, the image height
 * <I>H</I> (integer) and the image width <I>W</I> (integer). The rest of the
 * times file consists of <I>M</I> lines. The first field on each line gives
 * <I>t</I><SUB><I>i</I></SUB>, the time (sec, real number) at which measurement
 * <I>i</I> was taken, 0 &le; <I>i</I> &le; <I>M</I>&minus;1. The second field
 * on each line gives the name of the data file for measurement <I>i</I>.
 * <P>
 * <B>Data files.</B> Each data file is a binary file containing the pixel data
 * values for a certain measurement <I>i</I> taken at time
 * <I>t</I><SUB><I>i</I></SUB>. The file contains <I>H</I>&times;<I>W</I> spin
 * signal values <I>S</I><SUB><I>rc</I></SUB>(<I>t</I><SUB><I>i</I></SUB>),
 * where <I>r</I> is the pixel's row in the range 0 &le; <I>r</I> &le;
 * <I>H</I>&minus;1 and <I>c</I> is the pixel's column in the range 0 &le;
 * <I>c</I> &le; <I>W</I>&minus;1. The pixel data values are stored in row major
 * order -- first the pixels in the topmost row from left to right, then the
 * pixels in the second row from left to right, and so on.
 * <P>
 * Each pixel data value is stored in binary as a two-byte integer in the range
 * &minus;32768 &le; <I>S</I><SUB><I>rc</I></SUB>(<I>t</I><SUB><I>i</I></SUB>)
 * &le; 32767. The first two pixels (<I>S</I><SUB>00</SUB> and
 * <I>S</I><SUB>01</SUB>) do not contain actual data, but instead are set to the
 * values <I>S</I><SUB>00</SUB> = 0 and <I>S</I><SUB>01</SUB> = 1. These are
 * used to determine the byte order or <I>endianness</I> with which the pixel
 * data values are stored. If the first four bytes of the data file are 0, 0, 0,
 * 1, then the pixel data values are stored in <I>big-endian order</I> -- most
 * significant byte first. If the first four bytes of the data file are 0, 0, 1,
 * 0, then the pixel values are stored in <I>little-endian order</I> -- least
 * significant byte first. It is assumed that all data files use the same byte
 * order.
 * <P>
 * <B>Mask file.</B> The mask file is a binary file containing
 * <I>H</I>&times;<I>W</I> mask values <I>M</I><SUB><I>rc</I></SUB>.
 * <I>M</I><SUB><I>rc</I></SUB> is 1 if pixel (<I>r,c</I>) has data to be
 * analyzed. <I>M</I><SUB><I>rc</I></SUB> is 0 if pixel (<I>r,c</I>) is not to
 * be analyzed.
 * <P>
 * Each mask value is stored in binary as a two-byte integer. The first pixel
 * (<I>M</I><SUB>00</SUB>) does not contain actual mask data, but instead is set
 * to the value <I>M</I><SUB>00</SUB> = 1. This is used to determine the mask
 * values' byte order. If the first two bytes of the mask file are 0, 1, then
 * the mask values are stored in big-endian order. If the first two bytes of the
 * mask file are 1, 0, then the mask values are stored in little-endian order.
 *
 * @author  Alan Kaminsky
 * @version 21-Jun-2008
 */
public class CreateSignalDataSet
	{

// Prevent construction.

	private CreateSignalDataSet()
		{
		}

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Exception
		{
		byte b1, b2, b3, b4;

		// Parse command line arguments.
		if (args.length != 3) usage();
		File timesfile = new File (args[0]);
		File maskfile = new File (args[1]);
		File signalfile = new File (args[2]);

		// Read times file: height, width, time series, and list of data files.
		Scanner timesScanner = new Scanner (timesfile);
		Scanner lineScanner = new Scanner (timesScanner.nextLine());
		int H = lineScanner.nextInt();
		int W = lineScanner.nextInt();
		int P = H*W;
		ListSeries t_series = new ListSeries();
		ArrayList<File> datafile_list = new ArrayList<File>();
		while (timesScanner.hasNextLine())
			{
			lineScanner = new Scanner (timesScanner.nextLine());
			t_series.add (lineScanner.nextDouble());
			datafile_list.add (new File (lineScanner.next()));
			}
		timesScanner.close();
		int M = t_series.length();

		// Read mask file.
		boolean[] mask = new boolean [P];
		ByteBuffer buf = ByteBuffer.allocateDirect (2*P);
		fillBuffer (buf, maskfile);
		b1 = buf.get();
		b2 = buf.get();
		if (b1 == 0 && b2 == 1)
			{
			buf.order (ByteOrder.BIG_ENDIAN);
			}
		else if (b1 == 1 && b2 == 0)
			{
			buf.order (ByteOrder.LITTLE_ENDIAN);
			}
		else
			{
			throw new IOException
				("Illegal start of mask file \""+maskfile+"\": "+b1+", "+b2);
			}
		for (int i = 1; i < P; ++ i)
			{
			short val = buf.getShort();
			if (val == 0)
				{
				mask[i] = false;
				}
			else if (val == 1)
				{
				mask[i] = true;
				}
			else
				{
				throw new IOException
					("Illegal value "+val+" for pixel "+i+" in mask file \""+
					 maskfile+"\"");
				}
			}

		// Set up arrays to hold pixel signal data.
		short[][] S_array = new short [P] [];
		for (int i = 0; i < P; ++ i)
			{
			if (mask[i])
				{
				S_array[i] = new short [M];
				}
			}

		// Read data files.
		for (int j = 0; j < M; ++ j)
			{
			File datafile = datafile_list.get(j);
			fillBuffer (buf, datafile);
			b1 = buf.get();
			b2 = buf.get();
			b3 = buf.get();
			b4 = buf.get();
			if (b1 == 0 && b2 == 0 && b3 == 0 && b4 == 1)
				{
				buf.order (ByteOrder.BIG_ENDIAN);
				}
			if (b1 == 0 && b2 == 0 && b3 == 1 && b4 == 0)
				{
				buf.order (ByteOrder.LITTLE_ENDIAN);
				}
			else
				{
				throw new IOException
					("Illegal start of data file \""+datafile+"\": "+b1+", "+
					 b2+", "+b3+", "+b4);
				}
			for (int i = 2; i < P; ++ i)
				{
				short val = buf.getShort();
				if (mask[i])
					{
					S_array[i][j] = val;
					}
				}
			}

		// Write height, width, and time series to spin signal data set.
		SignalDataSetWriter writer =
			new SignalDataSetWriter (signalfile, H, W, t_series);

		// Write pixel signal data to spin signal data set.
		for (int i = 0; i < P; ++ i)
			{
			if (mask[i])
				{
				writer.addPixelSignal (new PixelSignal (0, i, S_array[i]));
				}
			}

		// All done!
		writer.close();
		}

// Hidden operations.

	/**
	 * Fill the given buffer with the contents of the given file.
	 */
	private static void fillBuffer
		(ByteBuffer buf,
		 File file)
		throws IOException
		{
		FileInputStream fis = new FileInputStream (file);
		FileChannel fchan = fis.getChannel();
		buf.clear();
		while (buf.hasRemaining())
			{
			if (fchan.read (buf) == -1)
				{
				throw new EOFException
					("Unexpected EOF in file \""+file+"\"");
				}
			}
		fis.close();
		buf.flip();
		}

	/**
	 * Print a usage message and exit.
	 */
	private static void usage()
		{
		System.err.println ("Usage: java edu.rit.mri.CreateSignalDataSet <timesfile> <maskfile> <signalfile>");
		System.err.println ("<timesfile> = Input times file");
		System.err.println ("<maskfile> = Input mask file");
		System.err.println ("<signalfile> = Output spin signal data set file");
		System.exit (1);
		}

	}
