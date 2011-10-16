//******************************************************************************
//
// File:    DrawingFrame.java
// Package: edu.rit.draw
// Unit:    Class edu.rit.draw.DrawingFrame
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

package edu.rit.draw;

import edu.rit.swing.Displayable;
import edu.rit.swing.DisplayableFrame;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;

import java.io.File;
import java.io.IOException;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;

/**
 * Class DrawingFrame is a {@linkplain edu.rit.swing.DisplayableFrame
 * DisplayableFrame} for viewing a {@linkplain Drawing} object on the screen.
 * The drawing object is specified when the viewer is constructed.
 * <P>
 * To display a drawing in a window on the screen, construct an instance of
 * class DrawingFrame; call <TT>setSize()</TT> to set the viewer window's size
 * if desired; and call <TT>setVisible(true)</TT> to display the viewer window.
 * <P>
 * The viewer window includes a "File" menu with menu items for saving the plot
 * in a PNG file, saving the plot in a PostScript file, and closing the window.
 * <P>
 * The viewer window includes a "View" menu with menu items for zooming in and
 * out. When the drawing is saved in a PNG file, the drawing is saved at its
 * zoomed size. When the drawing is saved in a PostScript file, the drawing is
 * saved at its normal (unzoomed) size.
 * <P>
 * When the viewer window is closed, either by using the "File" menu or by
 * clicking the window's close box, one of three possible operations is
 * performed: hide the window; hide and dispose the window; or exit the
 * application. The close operation is specified when the viewer window is
 * constructed. The default close operation is to hide and dispose the window.
 *
 * @author  Alan Kaminsky
 * @version 20-Dec-2007
 */
class DrawingFrame
	extends DisplayableFrame
	{

// Hidden data members.

	private Drawing myDrawing;

// Exported constructors.

	/**
	 * Construct a new drawing viewer. The frame has no title. The frame is
	 * initially invisible. When the frame is closed, the frame is hidden and
	 * disposed.
	 *
	 * @param  theDrawing  Drawing object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public DrawingFrame
		(Drawing theDrawing)
		{
		super (theDrawing);
		initialize (theDrawing);
		}

	/**
	 * Construct a new drawing viewer with the given title. The frame is
	 * initially invisible. When the frame is closed, the frame is hidden and
	 * disposed.
	 *
	 * @param  theTitle    Title.
	 * @param  theDrawing  Drawing object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 */
	public DrawingFrame
		(String theTitle,
		 Drawing theDrawing)
		{
		super (theTitle, theDrawing);
		initialize (theDrawing);
		}

	/**
	 * Construct a new drawing viewer. The frame has no title. The frame is
	 * initially invisible. When the frame is closed, <TT>theCloseOperation</TT>
	 * specifies what to do.
	 *
	 * @param  theDrawing         Drawing object.
	 * @param  theCloseOperation  What to do when the frame is closed, one of
	 *                            the following: <TT>JFrame.HIDE_ON_CLOSE</TT>,
	 *                            <TT>JFrame.DISPOSE_ON_CLOSE</TT>, or
	 *                            <TT>JFrame.EXIT_ON_CLOSE</TT>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theCloseOperation</TT> is
	 *     invalid.
	 */
	public DrawingFrame
		(Drawing theDrawing,
		 int theCloseOperation)
		{
		super (theDrawing, theCloseOperation);
		initialize (theDrawing);
		}

	/**
	 * Construct a new drawing viewer with the given title. The frame is
	 * initially invisible. When the frame is closed, <TT>theCloseOperation</TT>
	 * specifies what to do.
	 *
	 * @param  theTitle           Title.
	 * @param  theDrawing         Drawing object.
	 * @param  theCloseOperation  What to do when the frame is closed, one of
	 *                            the following: <TT>JFrame.HIDE_ON_CLOSE</TT>,
	 *                            <TT>JFrame.DISPOSE_ON_CLOSE</TT>, or
	 *                            <TT>JFrame.EXIT_ON_CLOSE</TT>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDrawing</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theCloseOperation</TT> is
	 *     invalid.
	 */
	public DrawingFrame
		(String theTitle,
		 Drawing theDrawing,
		 int theCloseOperation)
		{
		super (theTitle, theDrawing, theCloseOperation);
		initialize (theDrawing);
		}

	/**
	 * Initialize a newly constructed drawing viewer.
	 */
	private void initialize
		(Drawing theDrawing)
		{
		// Record drawing.
		if (theDrawing == null)
			{
			throw new NullPointerException
				("DrawingFrame(): theDrawing is null");
			}
		myDrawing = theDrawing;

		// Add item to File menu.
		JMenuBar menubar = getJMenuBar();
		JMenu fileMenu = menubar.getMenu (0);
		JMenuItem saveAsDrawingItem =
			new JMenuItem ("Save as Drawing...", KeyEvent.VK_D);
		saveAsDrawingItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_D, InputEvent.CTRL_MASK));
		saveAsDrawingItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doSaveAsDrawing();
					}
				});
		fileMenu.insert (saveAsDrawingItem, 2);

		// Add item to View menu.
		JMenu viewMenu = menubar.getMenu (1);
		JMenuItem windowTitleItem =
			new JMenuItem ("Window Title...", KeyEvent.VK_T);
		windowTitleItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_T, InputEvent.CTRL_MASK));
		windowTitleItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doWindowTitle();
					}
				});
		viewMenu.add (windowTitleItem);
		}

// Exported operations.

	/**
	 * Tell this displayable object viewer to display the given displayable
	 * object.
	 *
	 * @param  theDisplayable  Displayable object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 * @exception  ClassCastException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is not an
	 *     instance of class {@linkplain Drawing}.
	 */
	public void display
		(Displayable theDisplayable)
		{
		myDrawing = (Drawing) theDisplayable;
		super.display (theDisplayable);
		}

// Hidden operations.

	/**
	 * "Save as Drawing..." menu item.
	 */
	private void doSaveAsDrawing()
		{
		saveFile
			(".dwg",
			 "Serialized drawing object files",
			 new DisplayableFrame.FileSaver()
				{
				public void saveFile (File file) throws IOException
					{
					Drawing.write (myDrawing, file);
					}
				});
		}

	/**
	 * "Window Title..." menu item.
	 */
	private void doWindowTitle()
		{
		String newTitle =
			JOptionPane.showInputDialog
				(this,
				 "Window title:",
				 myDrawing.getTitle());
		if (newTitle != null)
			{
			myDrawing.setTitle (newTitle);
			setTitle (newTitle);
			}
		}

	}
