//******************************************************************************
//
// File:    TranscriptFrame.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.TranscriptFrame
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

import java.io.IOException;
import java.io.OutputStream;

import javax.swing.JFrame;
import javax.swing.JTextArea;

/**
 * Class TranscriptFrame is a Java Swing frame (window) that displays plain text
 * written to an output stream. It is a {@link javax.swing.JFrame
 * </CODE>JFrame<CODE>} containing a {@link Transcript </CODE>Transcript<CODE>}.
 *
 * @author  Alan Kaminsky
 * @version 20-Sep-2005
 */
public class TranscriptFrame
	extends JFrame
	{

// Exported constants.

	/**
	 * Default font: sanserif, plain, 9 point.
	 */
	public static final Font DEFAULT_FONT =
		Transcript.DEFAULT_FONT;

	/**
	 * Default maximum number of lines: 24.
	 */
	public static final int DEFAULT_MAX_LINES =
		Transcript.DEFAULT_MAX_LINES;

	/**
	 * Default visible number of lines: 24.
	 */
	public static final int DEFAULT_VISIBLE_LINES =
		Transcript.DEFAULT_VISIBLE_LINES;

	/**
	 * Default visible number of columns: 80.
	 */
	public static final int DEFAULT_VISIBLE_COLUMNS =
		Transcript.DEFAULT_VISIBLE_COLUMNS;

// Hidden data members.

	private Transcript myTranscript;

// Exported constructors.

	/**
	 * Construct a new transcript frame. The transcript will use the default
	 * font (9 point monospaced). The transcript will display the default
	 * maximum number of lines (24). The transcript will have the default number
	 * of lines and columns visible (24x80).
	 *
	 * @param  theTitle  Frame title.
	 */
	public TranscriptFrame
		(String theTitle)
		{
		this
			(theTitle,
			 DEFAULT_FONT,
			 DEFAULT_MAX_LINES,
			 DEFAULT_VISIBLE_LINES,
			 DEFAULT_VISIBLE_COLUMNS);
		}

	/**
	 * Construct a new transcript frame that will use the given font. The
	 * transcript will display the default maximum number of lines (24). The
	 * transcript will have the default number of lines and columns visible
	 * (24x80).
	 *
	 * @param  theTitle  Frame title.
	 * @param  theFont   Font to use for the transcript display.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFont</TT> is null.
	 */
	public TranscriptFrame
		(String theTitle,
		 Font theFont)
		{
		this
			(theTitle,
			 theFont,
			 DEFAULT_MAX_LINES,
			 DEFAULT_VISIBLE_LINES,
			 DEFAULT_VISIBLE_COLUMNS);
		}

	/**
	 * Construct a new transcript frame that will use the given font and that
	 * will display the given maximum number of lines. The transcript will have
	 * the default number of lines and columns visible (24x80).
	 *
	 * @param  theTitle     Frame title.
	 * @param  theFont      Font to use for the transcript display.
	 * @param  theMaxLines  Maximum number of lines displayed.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theFont</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theMaxLines</TT> is less than or
	 *     equal to 0.
	 */
	public TranscriptFrame
		(String theTitle,
		 Font theFont,
		 int theMaxLines)
		{
		this
			(theTitle,
			 theFont,
			 theMaxLines,
			 DEFAULT_VISIBLE_LINES,
			 DEFAULT_VISIBLE_COLUMNS);
		}

	/**
	 * Construct a new transcript frame that will use the given font, that will
	 * display the given maximum number of lines, and that will have the given
	 * number of lines and columns visible.
	 *
	 * @param  theTitle           Frame title.
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
	public TranscriptFrame
		(String theTitle,
		 Font theFont,
		 int theMaxLines,
		 int theVisibleLines,
		 int theVisibleColumns)
		{
		super (theTitle);
		myTranscript =
			new Transcript
				(theFont,
				 theMaxLines,
				 theVisibleLines,
				 theVisibleColumns);
		getContentPane().add (myTranscript);
		pack();
		}

// Exported operations.

	/**
	 * Returns the {@link javax.swing.JTextArea </CODE>JTextArea<CODE>}
	 * component used to display the text of this transcript.
	 */
	public JTextArea getTextArea()
		{
		return myTranscript.getTextArea();
		}

	/**
	 * Returns the {@link java.io.OutputStream </CODE>OutputStream<CODE>} used
	 * to write into this transcript.
	 */
	public OutputStream getOutputStream()
		{
		return myTranscript.getOutputStream();
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
		myTranscript.clear();
		}

	}
