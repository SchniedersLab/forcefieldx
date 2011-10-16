//******************************************************************************
//
// File:    DisplayableFrame.java
// Package: edu.rit.swing
// Unit:    Class edu.rit.swing.DisplayableFrame
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

package edu.rit.swing;

import java.awt.Cursor;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import java.awt.geom.Rectangle2D;

import java.io.File;
import java.io.IOException;

import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.KeyStroke;

import javax.swing.filechooser.FileFilter;

/**
 * Class DisplayableFrame is a Swing {@link javax.swing.JFrame
 * </CODE>JFrame<CODE>} for viewing a {@link Displayable
 * </CODE>Displayable<CODE>} object on the screen. The displayable object is
 * specified when the viewer is constructed.
 * <P>
 * To display a displayable object in a window on the screen, construct an
 * instance of class DisplayableFrame; call <TT>setSize()</TT> to set the viewer
 * window's size if desired; and call <TT>setVisible(true)</TT> to display the
 * viewer window.
 * <P>
 * The viewer window includes a "File" menu with menu items for saving the
 * display in a PNG file, saving the display in a PostScript file, and closing
 * the window.
 * <P>
 * The viewer window includes a "View" menu with menu items for zooming in and
 * out. When the display is saved in a PNG file, the display is saved at its
 * zoomed size. When the display is saved in a PostScript file, the display is
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
public class DisplayableFrame
	extends JFrame
	{

// Hidden data members.

	private Displayable myDisplayable;
	private DisplayablePanel myDisplayablePanel;
	private JFileChooser myChooser;
	private FileTypeFilter myFilter;
	private JMenuItem myCloseMenuItem;
	private int myCloseOperation;

	private int myZoomCount = 0;
	private double myZoomFactor = 1.0;
	private JCheckBoxMenuItem myZoom13Item;
	private JCheckBoxMenuItem myZoom25Item;
	private JCheckBoxMenuItem myZoom50Item;
	private JCheckBoxMenuItem myZoom100Item;
	private JCheckBoxMenuItem myZoom200Item;
	private JCheckBoxMenuItem myZoom400Item;
	private JCheckBoxMenuItem myZoom800Item;

// Exported constructors.

	/**
	 * Construct a new displayable object viewer. The frame has no title. The
	 * frame is initially invisible. When the frame is closed, the frame is
	 * hidden and disposed.
	 *
	 * @param  theDisplayable  Displayable object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 */
	public DisplayableFrame
		(Displayable theDisplayable)
		{
		super();
		initialize (theDisplayable, JFrame.DISPOSE_ON_CLOSE);
		}

	/**
	 * Construct a new displayable object viewer with the given title. The frame
	 * is initially invisible. When the frame is closed, the frame is hidden and
	 * disposed.
	 *
	 * @param  theTitle        Title.
	 * @param  theDisplayable  Displayable object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 */
	public DisplayableFrame
		(String theTitle,
		 Displayable theDisplayable)
		{
		super (theTitle);
		initialize (theDisplayable, JFrame.DISPOSE_ON_CLOSE);
		}

	/**
	 * Construct a new displayable object viewer. The frame has no title. The
	 * frame is initially invisible. When the frame is closed,
	 * <TT>theCloseOperation</TT> specifies what to do.
	 *
	 * @param  theDisplayable     Displayable object.
	 * @param  theCloseOperation  What to do when the frame is closed, one of
	 *                            the following: <TT>JFrame.HIDE_ON_CLOSE</TT>,
	 *                            <TT>JFrame.DISPOSE_ON_CLOSE</TT>, or
	 *                            <TT>JFrame.EXIT_ON_CLOSE</TT>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theCloseOperation</TT> is
	 *     invalid.
	 */
	public DisplayableFrame
		(Displayable theDisplayable,
		 int theCloseOperation)
		{
		super();
		initialize (theDisplayable, theCloseOperation);
		}

	/**
	 * Construct a new displayable object viewer with the given title. The frame
	 * is initially invisible. When the frame is closed,
	 * <TT>theCloseOperation</TT> specifies what to do.
	 *
	 * @param  theTitle           Title.
	 * @param  theDisplayable     Displayable object.
	 * @param  theCloseOperation  What to do when the frame is closed, one of
	 *                            the following: <TT>JFrame.HIDE_ON_CLOSE</TT>,
	 *                            <TT>JFrame.DISPOSE_ON_CLOSE</TT>, or
	 *                            <TT>JFrame.EXIT_ON_CLOSE</TT>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theCloseOperation</TT> is
	 *     invalid.
	 */
	public DisplayableFrame
		(String theTitle,
		 Displayable theDisplayable,
		 int theCloseOperation)
		{
		super (theTitle);
		initialize (theDisplayable, theCloseOperation);
		}

	/**
	 * Initialize a newly constructed displayable object viewer.
	 *
	 * @param  theDisplayable     Displayable object.
	 * @param  theCloseOperation  What to do when the frame is closed, one of
	 *                            the following: <TT>JFrame.HIDE_ON_CLOSE</TT>,
	 *                            <TT>JFrame.DISPOSE_ON_CLOSE</TT>, or
	 *                            <TT>JFrame.EXIT_ON_CLOSE</TT>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theCloseOperation</TT> is
	 *     invalid.
	 */
	private void initialize
		(Displayable theDisplayable,
		 int theCloseOperation)
		{
		// Verify preconditions.
		if (theDisplayable == null)
			{
			throw new NullPointerException();
			}
		switch (theCloseOperation)
			{
			case JFrame.HIDE_ON_CLOSE:
			case JFrame.DISPOSE_ON_CLOSE:
			case JFrame.EXIT_ON_CLOSE:
				break;
			default:
				throw new IllegalArgumentException();
			}

		// Save arguments.
		myDisplayable = theDisplayable;
		myCloseOperation = theCloseOperation;

		// Set up file chooser.
		myChooser = new JFileChooser (System.getProperty ("user.dir"));
		myFilter = new FileTypeFilter();
		myChooser.setFileFilter (myFilter);

		// Create File menu.
		JMenu fileMenu = new JMenu ("File");
		fileMenu.setMnemonic (KeyEvent.VK_F);

		JMenuItem saveAsPngItem =
			new JMenuItem ("Save as PNG...", KeyEvent.VK_S);
		saveAsPngItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_S, InputEvent.CTRL_MASK));
		saveAsPngItem.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					doSaveAsPng();
					}
				});
		fileMenu.add (saveAsPngItem);

		JMenuItem saveAsPostScriptItem =
			new JMenuItem ("Save as PostScript...", KeyEvent.VK_P);
		saveAsPostScriptItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_P, InputEvent.CTRL_MASK));
		saveAsPostScriptItem.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					doSaveAsPostScript();
					}
				});
		fileMenu.add (saveAsPostScriptItem);

		fileMenu.addSeparator();

		myCloseMenuItem = new JMenuItem ("Quit", KeyEvent.VK_Q);
		setDefaultCloseOperation (theCloseOperation);
		myCloseMenuItem.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					doClose();
					}
				});
		fileMenu.add (myCloseMenuItem);

		// Create View menu.
		JMenu viewMenu = new JMenu ("View");
		viewMenu.setMnemonic (KeyEvent.VK_V);

		JMenuItem zoomOutItem = new JMenuItem ("Zoom Out", KeyEvent.VK_Z);
		zoomOutItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_MINUS, InputEvent.CTRL_MASK));
		zoomOutItem.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (myZoomCount - 1);
					}
				});
		viewMenu.add (zoomOutItem);

		JMenuItem zoomInItem = new JMenuItem ("Zoom In", KeyEvent.VK_O);
		zoomInItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_EQUALS, InputEvent.CTRL_MASK));
		zoomInItem.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (myZoomCount + 1);
					}
				});
		viewMenu.add (zoomInItem);

		JMenuItem fitWindowItem = new JMenuItem ("Fit Window", KeyEvent.VK_F);
		fitWindowItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_E, InputEvent.CTRL_MASK));
		fitWindowItem.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					doFitWindow();
					}
				});
		viewMenu.add (fitWindowItem);

		viewMenu.addSeparator();

		myZoom13Item = new JCheckBoxMenuItem ("13%");
		myZoom13Item.setAccelerator
			(KeyStroke.getKeyStroke
				(KeyEvent.VK_8,
				 InputEvent.CTRL_MASK+InputEvent.SHIFT_MASK));
		myZoom13Item.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (-18);
					}
				});
		viewMenu.add (myZoom13Item);

		myZoom25Item = new JCheckBoxMenuItem ("25%");
		myZoom25Item.setAccelerator
			(KeyStroke.getKeyStroke
				(KeyEvent.VK_4,
				 InputEvent.CTRL_MASK+InputEvent.SHIFT_MASK));
		myZoom25Item.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (-12);
					}
				});
		viewMenu.add (myZoom25Item);

		myZoom50Item = new JCheckBoxMenuItem ("50%");
		myZoom50Item.setAccelerator
			(KeyStroke.getKeyStroke
				(KeyEvent.VK_2,
				 InputEvent.CTRL_MASK+InputEvent.SHIFT_MASK));
		myZoom50Item.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (-6);
					}
				});
		viewMenu.add (myZoom50Item);

		myZoom100Item = new JCheckBoxMenuItem ("100%");
		myZoom100Item.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_1, InputEvent.CTRL_MASK));
		myZoom100Item.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (0);
					}
				});
		viewMenu.add (myZoom100Item);

		myZoom200Item = new JCheckBoxMenuItem ("200%");
		myZoom200Item.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_2, InputEvent.CTRL_MASK));
		myZoom200Item.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (6);
					}
				});
		viewMenu.add (myZoom200Item);

		myZoom400Item = new JCheckBoxMenuItem ("400%");
		myZoom400Item.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_4, InputEvent.CTRL_MASK));
		myZoom400Item.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (12);
					}
				});
		viewMenu.add (myZoom400Item);

		myZoom800Item = new JCheckBoxMenuItem ("800%");
		myZoom800Item.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_8, InputEvent.CTRL_MASK));
		myZoom800Item.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					setZoom (18);
					}
				});
		viewMenu.add (myZoom800Item);

		myZoom100Item.setSelected (true);

		viewMenu.addSeparator();

		JMenuItem imageInfoItem =
			new JMenuItem ("Image Info...", KeyEvent.VK_I);
		imageInfoItem.setAccelerator
			(KeyStroke.getKeyStroke (KeyEvent.VK_I, InputEvent.CTRL_MASK));
		imageInfoItem.addActionListener
			(new ActionListener()
				{
				public void actionPerformed
					(ActionEvent e)
					{
					doImageInfo();
					}
				});
		viewMenu.add (imageInfoItem);

		// Create menu bar.
		JMenuBar menubar = new JMenuBar();
		menubar.add (fileMenu);
		menubar.add (viewMenu);
		setJMenuBar (menubar);

		// Create panel to show the displayable object.
		myDisplayablePanel = new DisplayablePanel (myDisplayable);
		JScrollPane theScrollPane =
			new JScrollPane
				(myDisplayablePanel,
				 JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
				 JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		getContentPane().add (theScrollPane);

		// Arrange to handle the window closing.
		addWindowListener
			(new WindowAdapter()
				{
				public void windowClosing
					(WindowEvent e)
					{
					doClose();
					}
				});

		// Lay out widgets.
		pack();
		}

// Exported operations.

	/**
	 * Set the operation to be performed when the frame is closed.
	 *
	 * @param  operation
	 *     One of the following constants defined in class javax.swing.JFrame:
	 *     <TT>HIDE_ON_CLOSE</TT>, <TT>DISPOSE_ON_CLOSE</TT>,
	 *     <TT>EXIT_ON_CLOSE</TT>.
	 */
	public void setDefaultCloseOperation
		(int operation)
		{
		switch (operation)
			{
			case JFrame.HIDE_ON_CLOSE:
			case JFrame.DISPOSE_ON_CLOSE:
			case JFrame.EXIT_ON_CLOSE:
				break;
			default:
				throw new IllegalArgumentException();
			}

		super.setDefaultCloseOperation (operation);

		myCloseOperation = operation;
		if (operation == JFrame.EXIT_ON_CLOSE)
			{
			myCloseMenuItem.setText ("Quit");
			myCloseMenuItem.setMnemonic (KeyEvent.VK_Q);
			myCloseMenuItem.setAccelerator
				(KeyStroke.getKeyStroke (KeyEvent.VK_Q, InputEvent.CTRL_MASK));
			}
		else
			{
			myCloseMenuItem.setText ("Close");
			myCloseMenuItem.setMnemonic (KeyEvent.VK_C);
			myCloseMenuItem.setAccelerator
				(KeyStroke.getKeyStroke (KeyEvent.VK_W, InputEvent.CTRL_MASK));
			}
		}

	/**
	 * Tell this displayable object viewer to display the given displayable
	 * object.
	 *
	 * @param  theDisplayable  Displayable object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theDisplayable</TT> is null.
	 */
	public void display
		(Displayable theDisplayable)
		{
		if (theDisplayable == null)
			{
			throw new NullPointerException();
			}
		myDisplayable = theDisplayable;
		myDisplayablePanel.display (theDisplayable);
		}

	/**
	 * Perform a file save operation. This method displays a file chooser dialog
	 * listing files with the given extension; lets the user pick one; if the
	 * chosen file exists, displays a dialog asking the user if it's okay to
	 * overwrite the file; if all is well, calls the given file saver object's
	 * <TT>saveFile()</TT> method, passing in the chosen file; if an IOException
	 * occurs, displays a dialog with an error message.
	 *
	 * @param  theExtension    File extension for the file type to save.
	 * @param  theDescription  Description of the file type to save.
	 * @param  theFileSaver    File saver object.
	 */
	public void saveFile
		(String theExtension,
		 String theDescription,
		 FileSaver theFileSaver)
		{
		File file = null;

		try
			{
			myFilter.setType (theExtension, theDescription);
			if (myChooser.showSaveDialog (this) == JFileChooser.APPROVE_OPTION)
				{
				file = myChooser.getSelectedFile();
				if (okayToWrite (file))
					{
					setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
					theFileSaver.saveFile (file);
					setCursor (Cursor.getDefaultCursor());
					}
				}
			}

		catch (Throwable exc)
			{
			setCursor (Cursor.getDefaultCursor());
			displayErrorMessage (file, exc);
			}
		}

// Exported helper classes.

	/**
	 * Interface DisplayableFrame.FileSaver specifies the interface for an
	 * object that saves information in a file.
	 *
	 * @author  Alan Kaminsky
	 * @version 04-Nov-2007
	 */
	public static interface FileSaver
		{
		/**
		 * Save information in the given file.
		 *
		 * @param  file  File.
		 *
		 * @exception  IOException
		 *     Thrown if an I/O error occurred.
		 */
		public void saveFile
			(File file)
			throws IOException;
		}

// Hidden operations.

	/**
	 * "Save as PNG..." menu item.
	 */
	private void doSaveAsPng()
		{
		saveFile
			(".png",
			 "PNG image files",
			 new FileSaver()
				{
				public void saveFile (File file) throws IOException
					{
					DisplayableIO.writeColorPNGFile
						(myDisplayable, file, myZoomFactor);
					}
				});
		}

	/**
	 * "Save as PostScript..." menu item.
	 */
	private void doSaveAsPostScript()
		{
		saveFile
			(".ps",
			 "PostScript files",
			 new FileSaver()
				{
				public void saveFile (File file) throws IOException
					{
					DisplayableIO.writePostScriptFile (myDisplayable, file);
					}
				});
		}

	/**
	 * Check if it's okay to write a file.
	 */
	private boolean okayToWrite
		(File file)
		{
		if (file.exists())
			{
			// File exists. Ask to overwrite.
			try
				{
				return
					JOptionPane.showConfirmDialog
						(this,
						 "Overwrite \"" + file.getName() + "\"?",
						 "File Exists",
						 JOptionPane.YES_NO_OPTION,
						 JOptionPane.QUESTION_MESSAGE)
						== JOptionPane.YES_OPTION;
				}
			catch (Throwable exc)
				{
				return false;
				}
			}
		else
			{
			// File does not exist.
			return true;
			}
		}

	/**
	 * Display an error message.
	 */
	private void displayErrorMessage
		(File file,
		 Throwable exc)
		{
		try
			{
			StringBuffer buf = new StringBuffer();
			buf.append ("Cannot write file");
			if (file != null)
				{
				buf.append (" \"");
				buf.append (file.getName());
				buf.append ("\"");
				}
			buf.append (" -- ");
			buf.append (exc.getClass().getName());
			if (exc.getMessage() != null)
				{
				buf.append (" -- ");
				buf.append (exc.getMessage());
				}
			JOptionPane.showMessageDialog
				(this,
				 buf.toString(),
				 "Error",
				 JOptionPane.ERROR_MESSAGE);
			}
		catch (Throwable exc2)
			{
			}
		}

	/**
	 * "Close" or "Quit" menu item.
	 */
	private void doClose()
		{
		switch (myCloseOperation)
			{
			case JFrame.HIDE_ON_CLOSE:
				setVisible (false);
				break;
			case JFrame.DISPOSE_ON_CLOSE:
				setVisible (false);
				dispose();
				break;
			case JFrame.EXIT_ON_CLOSE:
				System.exit (0);
				break;
			}
		}

	/**
	 * Set the zoom factor.
	 */
	private void setZoom
		(int zoomCount)
		{
		myZoomCount = Math.max (-18, Math.min (zoomCount, +18));
		myZoomFactor = Math.pow (2.0, myZoomCount / 6.0);
		myDisplayablePanel.zoom (myZoomFactor);
		myZoom13Item .setSelected (myZoomCount == -18);
		myZoom25Item .setSelected (myZoomCount == -12);
		myZoom50Item .setSelected (myZoomCount ==  -6);
		myZoom100Item.setSelected (myZoomCount ==   0);
		myZoom200Item.setSelected (myZoomCount ==   6);
		myZoom400Item.setSelected (myZoomCount ==  12);
		myZoom800Item.setSelected (myZoomCount ==  18);
		}

	/**
	 * "Fit Window" menu item.
	 */
	private void doFitWindow()
		{
		myDisplayablePanel.display (myDisplayable);
		pack();
		}

	/**
	 * "Image Info" menu item.
	 */
	private void doImageInfo()
		{
		Rectangle2D bb = myDisplayable.getBoundingBox();
		double w = bb.getWidth();
		double h = bb.getHeight();
		JOptionPane.showMessageDialog
			(this,
			 new String[]
				{"Width = " + w + " pt = " + (w/72.0) + " in",
				 "Height = " + h + " pt = " + (h/72.0) + " in"},
			 "Image Info",
			 JOptionPane.INFORMATION_MESSAGE);
		}

	}
