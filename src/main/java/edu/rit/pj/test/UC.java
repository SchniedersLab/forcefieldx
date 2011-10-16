//******************************************************************************
//
// File:    UC.java
// Package: edu.rit.pj.test
// Unit:    Class edu.rit.pj.test.UC
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
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

package edu.rit.pj.test;

import edu.rit.pj.ParallelIteration;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import java.io.IOException;
import java.io.InputStream;

import java.net.MalformedURLException;
import java.net.URL;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Class UC is a main program that counts lines, words, and bytes in the URLs
 * given on the command line. It is a multi-threaded program using Parallel
 * Java.
 * <P>
 * Usage: java edu.rit.pj.test.UC <I>url</I> . . .
 *
 * @author  Alan Kaminsky
 * @version 04-Jun-2007
 */
public class UC
	{

// Prevent construction.

	private UC()
		{
		}

// Global variables.

	static String[] urls;

	static AtomicLong totalLines = new AtomicLong (0L);
	static AtomicLong totalWords = new AtomicLong (0L);
	static AtomicLong totalBytes = new AtomicLong (0L);

// Main program.

	/**
	 * Main program.
	 */
	public static void main
		(String[] args)
		throws Throwable
		{
		long stopwatch = -System.currentTimeMillis();

		urls = args;

		new ParallelTeam (urls.length) .execute (new ParallelRegion()
			{
			public void run() throws Exception
				{
				// Process all URLs.
				execute (urls, new ParallelIteration<String>()
					{
					// Process one URL.
					public void run (String url)
						{
						try
							{
							// Open a data stream from the URL.
							InputStream stream = new URL (url) .openStream();

							// Initialize counters.
							long lines = 0L;
							long words = 0L;
							long bytes = 0L;
							boolean inaword = false;

							// Read data and update counters.
							int b;
							while ((b = stream.read()) != -1)
								{
								++ bytes;
								if (b == '\n') ++ lines;
								if (0x21 <= b && b <= 0x7e) // Non-whitespace
									{
									if (! inaword) ++ words;
									inaword = true;
									}
								else // Whitespace
									{
									inaword = false;
									}
								}

							// Close data stream.
							stream.close();

							// Accumulate totals.
							totalLines.addAndGet (lines);
							totalWords.addAndGet (words);
							totalBytes.addAndGet (bytes);

							// Print results.
							System.out.println
								(lines + "\t" +
								 words + "\t" +
								 bytes + "\t" +
								 url);
							}

						catch (MalformedURLException exc)
							{
							System.err.println
								("\"" + url + "\": Malformed URL");
							}
						catch (IOException exc)
							{
							System.err.println
								(url + ": I/O error");
							}
						}
					});
				}
			});

		// Print final totals.
		System.out.println
			(totalLines + "\t" +
			 totalWords + "\t" +
			 totalBytes + "\tTotal");

		stopwatch += System.currentTimeMillis();
		System.out.println (stopwatch + " msec");
		}

	}
