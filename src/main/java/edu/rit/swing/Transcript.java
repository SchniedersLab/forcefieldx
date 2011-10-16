//******************************************************************************
//
// File:    Transcript.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.Transcript
//
// This Java source file is copyright (C) 2005 by the Rochester Institute of
// Technology. All rights reserved. For further information, contact the author,
// Alan Kaminsky, at ark@it.rit.edu.
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

package edu.rit.swing;

import java.awt.Font;
import java.awt.Insets;

import java.io.IOException;
import java.io.OutputStream;

import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;

/**
 * Class Transcript is a Java Swing component that displays plain text written
 * to an output stream. It is a {@link javax.swing.JScrollPane
 * </CODE>JScrollPane<CODE>} containing a {@link javax.swing.JTextArea
 * </CODE>JTextArea<CODE>}.
 *
 * @author  Alan Kaminsky
 * @version 02-Oct-2005
 */
public class Transcript
	extends JScrollPane
	{

// Exported constants.

	/**
	 * Default font: sanserif, plain, 9 point.
	 */
	public static final Font DEFAULT_FONT =
		new Font ("monospaced", Font.PLAIN, 9);

	/**
	 * Default maximum number of lines: 24.
	 */
	public static final int DEFAULT_MAX_LINES = 24;

	/**
	 * Default visible number of lines: 24.
	 */
	public static final int DEFAULT_VISIBLE_LINES = 24;

	/**
	 * Default visible number of columns: 80.
	 */
	public static final int DEFAULT_VISIBLE_COLUMNS = 80;

// Hidden constants.

	private static final int MARGIN = 3;
	private static final int INCR = 1024;

// Hidden data members.

	private JTextArea myTextArea;
	private JScrollBar myVerticalScrollBar;

	private int myMaxLines;
	private int myLineCount;
	private int myMaxBytes;
	private int myByteCount;
	private byte[] myBuffer;

	private TranscriptOutputStream myOutputStream;

// Hidden helper classes.

	/**
	 * Class Transcript.TranscriptOutputStream is an {@link java.io.OutputStream
	 * </CODE>OutputStream<CODE>} used to write into a {@link Transcript
	 * </CODE>Transcript<CODE>}.
	 *
	 * @author  Alan Kaminsky
	 * @version 13-Apr-2003
	 */
	private class TranscriptOutputStream
		extends OutputStream
		{
		private boolean flushed = true;

		/**
		 * Write the given byte to this output stream.
		 */
		public synchronized void write
			(int b)
			throws IOException
			{
			// Expand myBuffer if necessary.
			if (myByteCount == myMaxBytes)
				{
				int newMaxBytes = myMaxBytes + INCR;
				byte[] newBuffer = new byte [newMaxBytes];
				System.arraycopy (myBuffer, 0, newBuffer, 0, myByteCount);
				myMaxBytes = newMaxBytes;
				myBuffer = newBuffer;
				}

			// Append byte to buffer.
			myBuffer[myByteCount++] = (byte)(b & 0xFF);
			flushed = false;

			// Handle a newline.
			if (b == '\n')
				{
				if (myLineCount < myMaxLines)
					{
					// Maximum lines not reached yet.
					++ myLineCount;
					}
				else
					{
					// Get rid of the oldest line.

					// Find the index of the byte after the first newline in the
					// buffer.
					int i = 0;
					while (myBuffer[i++] != '\n');

					// Shift buffer bytes backwards.
					myByteCount -= i;
					System.arraycopy (myBuffer, i, myBuffer, 0, myByteCount);
					}

				// Flush the buffer into the display.
				flush();
				}
			}

		/**
		 * Flush this output stream.
		 */
		public synchronized void flush()
			throws IOException
			{
			// If no bytes were written since the last flush, do nothing.
			if (! flushed)
				{
				flushed = true;
				displayBuffer();
				}
			}

		}

// Exported constructors.

	/**
	 * Construct a new transcript. The transcript will use the default font (9
	 * point monospaced). The transcript will display the default maximum number
	 * of lines (24). The transcript will have the default number of lines and
	 * columns visible (24x80).
	 */
	public Transcript()
		{
		this
			(DEFAULT_FONT,
			 DEFAULT_MAX_LINES,
			 DEFAULT_VISIBLE_LINES,
			 DEFAULT_VISIBLE_COLUMNS);
		}

	/**
	 * Construct a new transcript that will use the given font. The transcript
	 * will display the default maximum number of lines (24). The transcript
	 * will have the default number of lines and columns visible (24x80).
	 *
	 * @param  theFont  Font to use for the transcript display.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFont</TT> is null.
	 */
	public Transcript
		(Font theFont)
		{
		this
			(theFont,
			 DEFAULT_MAX_LINES,
			 DEFAULT_VISIBLE_LINES,
			 DEFAULT_VISIBLE_COLUMNS);
		}

	/**
	 * Construct a new transcript that will use the given font and that will
	 * display the given maximum number of lines. The transcript will have the
	 * default number of lines and columns visible (24x80).
	 *
	 * @param  theFont      Font to use for the transcript display.
	 * @param  theMaxLines  Maximum number of lines displayed.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFont</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theMaxLines</TT> is less than or
	 *     equal to 0.
	 */
	public Transcript
		(Font theFont,
		 int theMaxLines)
		{
		this
			(theFont,
			 theMaxLines,
			 DEFAULT_VISIBLE_LINES,
			 DEFAULT_VISIBLE_COLUMNS);
		}

	/**
	 * Construct a new transcript that will use the given font, that will
	 * display the given maximum number of lines, and that will have the given
	 * number of lines and columns visible.
	 *
	 * @param  theFont            Font to use for the transcript display.
	 * @param  theMaxLines        Maximum number of lines displayed.
	 * @param  theVisibleLines    Number of lines visible.
	 * @param  theVisibleColumns  Number of columns visible.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFont</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theMaxLines</TT> is less than or
	 *     equal to 0. Thrown if <TT>theVisibleLines</TT> is less than or equal
	 *     to 0. Thrown if <TT>theVisibleColumns</TT> is less than or equal to
	 *     0.
	 */
	public Transcript
		(Font theFont,
		 int theMaxLines,
		 int theVisibleLines,
		 int theVisibleColumns)
		{
		super
			(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
			 JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

		if (theFont == null)
			{
			throw new NullPointerException();
			}
		if (theMaxLines <= 0 || theVisibleLines <= 0 || theVisibleColumns <= 0)
			{
			throw new IllegalArgumentException();
			}

		myTextArea = new JTextArea (theVisibleLines, theVisibleColumns);
		myTextArea.setEditable (false);
		myTextArea.setFont (theFont);
		myTextArea.setLineWrap (true);
		myTextArea.setWrapStyleWord (false);
		myTextArea.setMargin (new Insets (MARGIN, MARGIN, MARGIN, MARGIN));
		setViewportView (myTextArea);

		myVerticalScrollBar = getVerticalScrollBar();

		myMaxLines = theMaxLines;
		myLineCount = 0;
		myMaxBytes = INCR;
		myByteCount = 0;
		myBuffer = new byte [INCR];

		myOutputStream = new TranscriptOutputStream();
		}

// Exported operations.

	/**
	 * Returns the {@link javax.swing.JTextArea </CODE>JTextArea<CODE>}
	 * component used to display the text of this transcript.
	 */
	public JTextArea getTextArea()
		{
		return myTextArea;
		}

	/**
	 * Returns the {@link java.io.OutputStream </CODE>OutputStream<CODE>} used
	 * to write into this transcript.
	 */
	public OutputStream getOutputStream()
		{
		return myOutputStream;
		}

	/**
	 * Clear this transcript.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	public void clear()
		throws IOException
		{
		synchronized (myOutputStream)
			{
			myLineCount = 0;
			myByteCount = 0;
			displayBuffer();
			}
		}

// Hidden operations.

	/**
	 * Display this transcript's buffer. Assumes the calling thread is
	 * synchronized on myOutputStream.
	 *
	 * @exception  IOException
	 *     Thrown if an I/O error occurred.
	 */
	private void displayBuffer()
		throws IOException
		{
		// Display buffered text.
		myTextArea.setText (new String (myBuffer, 0, myByteCount));

		// Scroll to bottom of text.
		try
			{
			SwingUtilities.invokeLater
				(new Runnable()
					{
					public void run()
						{
						myVerticalScrollBar.setValue
							(myVerticalScrollBar.getMaximum());
						}
					});
			}
		catch (Exception exc)
			{
			IOException exc2 = new IOException();
			exc2.initCause (exc);
			throw exc2;
			}
		}

	}
