//******************************************************************************
//
// File:    PlotFrame.java
// Package: edu.rit.numeric.plot
// Unit:    Class edu.rit.numeric.plot.PlotFrame
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

package edu.rit.numeric.plot;

import edu.rit.swing.Displayable;
import edu.rit.swing.DisplayableFrame;

import java.awt.Cursor;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import java.io.File;
import java.io.IOException;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;

/**
 * Class PlotFrame is a {@linkplain edu.rit.swing.DisplayableFrame
 * DisplayableFrame} for viewing a {@linkplain Plot} object on the screen. The
 * plot object is specified when the viewer is constructed.
 * <P>
 * To display a plot in a window on the screen, construct an instance of class
 * PlotFrame; call <TT>setSize()</TT> to set the viewer window's size if
 * desired; and call <TT>setVisible(true)</TT> to display the viewer window.
 * <P>
 * The viewer window includes a "File" menu with menu items for saving the plot
 * in a PNG file, saving the plot in a PostScript file, and closing the window.
 * <P>
 * The viewer window includes a "View" menu with menu items for zooming in and
 * out. When the plot is saved in a PNG file, the plot is saved at its zoomed
 * size. When the plot is saved in a PostScript file, the plot is saved at its
 * normal (unzoomed) size.
 * <P>
 * The viewer window includes a "Format" menu with menu items for setting the
 * various plot attributes interactively.
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
class PlotFrame
	extends DisplayableFrame
	{

// Hidden data members.

	private Plot myPlot;
	private TitleDialog myPlotTitleDialog;
	private PlotMarginsDialog myPlotMarginsDialog;
	private PlotAreaDialog myPlotAreaDialog;
	private AxisDialog myXAxisDialog;
	private TitleDialog myXAxisTitleDialog;
	private AxisDialog myYAxisDialog;
	private TitleDialog myYAxisTitleDialog;

// Exported constructors.

	/**
	 * Construct a new plot viewer. The frame has no title. The frame is
	 * initially invisible. When the frame is closed, the frame is hidden and
	 * disposed.
	 *
	 * @param  thePlot  Plot object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePlot</TT> is null.
	 */
	public PlotFrame
		(Plot thePlot)
		{
		super (thePlot);
		initialize (thePlot);
		}

	/**
	 * Construct a new plot viewer with the given title. The frame is initially
	 * invisible. When the frame is closed, the frame is hidden and disposed.
	 *
	 * @param  theTitle  Title.
	 * @param  thePlot   Plot object.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePlot</TT> is null.
	 */
	public PlotFrame
		(String theTitle,
		 Plot thePlot)
		{
		super (theTitle, thePlot);
		initialize (thePlot);
		}

	/**
	 * Construct a new plot viewer. The frame has no title. The frame is
	 * initially invisible. When the frame is closed, <TT>theCloseOperation</TT>
	 * specifies what to do.
	 *
	 * @param  thePlot            Plot object.
	 * @param  theCloseOperation  What to do when the frame is closed, one of
	 *                            the following: <TT>JFrame.HIDE_ON_CLOSE</TT>,
	 *                            <TT>JFrame.DISPOSE_ON_CLOSE</TT>, or
	 *                            <TT>JFrame.EXIT_ON_CLOSE</TT>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePlot</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theCloseOperation</TT> is
	 *     invalid.
	 */
	public PlotFrame
		(Plot thePlot,
		 int theCloseOperation)
		{
		super (thePlot, theCloseOperation);
		initialize (thePlot);
		}

	/**
	 * Construct a new plot viewer with the given title. The frame is initially
	 * invisible. When the frame is closed, <TT>theCloseOperation</TT> specifies
	 * what to do.
	 *
	 * @param  theTitle           Title.
	 * @param  thePlot            Plot object.
	 * @param  theCloseOperation  What to do when the frame is closed, one of
	 *                            the following: <TT>JFrame.HIDE_ON_CLOSE</TT>,
	 *                            <TT>JFrame.DISPOSE_ON_CLOSE</TT>, or
	 *                            <TT>JFrame.EXIT_ON_CLOSE</TT>.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePlot</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theCloseOperation</TT> is
	 *     invalid.
	 */
	public PlotFrame
		(String theTitle,
		 Plot thePlot,
		 int theCloseOperation)
		{
		super (theTitle, thePlot, theCloseOperation);
		initialize (thePlot);
		}

	/**
	 * Initialize a newly constructed plot viewer.
	 *
	 * @param  thePlot  Plot.
	 */
	private void initialize
		(Plot thePlot)
		{
		// Record plot.
		if (thePlot == null)
			{
			throw new NullPointerException
				("PlotFrame(): thePlot is null");
			}
		myPlot = thePlot;

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

		// Create Format menu.
		JMenu formatMenu = new JMenu ("Format");
		formatMenu.setMnemonic (KeyEvent.VK_O);

		JMenuItem titleItem =
			new JMenuItem ("Plot Title...", KeyEvent.VK_P);
		titleItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doPlotTitle();
					}
				});
		formatMenu.add (titleItem);

		JMenuItem marginsItem =
			new JMenuItem ("Margins...", KeyEvent.VK_M);
		marginsItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doMargins();
					}
				});
		formatMenu.add (marginsItem);

		JMenuItem areaItem =
			new JMenuItem ("Plot Area...", KeyEvent.VK_A);
		areaItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doPlotArea();
					}
				});
		formatMenu.add (areaItem);

		JMenuItem xAxisItem =
			new JMenuItem ("X Axis...", KeyEvent.VK_X);
		xAxisItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doXAxis();
					}
				});
		formatMenu.add (xAxisItem);

		JMenuItem xAxisTitleItem =
			new JMenuItem ("X Axis Title...", KeyEvent.VK_I);
		xAxisTitleItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doXAxisTitle();
					}
				});
		formatMenu.add (xAxisTitleItem);

		JMenuItem yAxisItem =
			new JMenuItem ("Y Axis...", KeyEvent.VK_Y);
		yAxisItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doYAxis();
					}
				});
		formatMenu.add (yAxisItem);

		JMenuItem yAxisTitleItem =
			new JMenuItem ("Y Axis Title...", KeyEvent.VK_T);
		yAxisTitleItem.addActionListener (new ActionListener()
				{
				public void actionPerformed (ActionEvent e)
					{
					doYAxisTitle();
					}
				});
		formatMenu.add (yAxisTitleItem);

		// Add Format menu to menu bar.
		menubar.add (formatMenu);

		// Arrange to handle the window closing.
		addWindowListener (new WindowAdapter()
			{
			public void windowClosing (WindowEvent e)
				{
				doClose();
				}
			});
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
	 *     instance of class {@linkplain Plot}.
	 */
	public void display
		(Displayable theDisplayable)
		{
		myPlot = (Plot) theDisplayable;
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
					Plot.write (myPlot, file);
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
				 myPlot.getTitle());
		if (newTitle != null)
			{
			myPlot.frameTitle (newTitle);
			setTitle (newTitle);
			}
		}

	/**
	 * "Plot Title..." menu item.
	 */
	private void doPlotTitle()
		{
		// Create plot title dialog if necessary.
		if (myPlotTitleDialog == null)
			{
			setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
			myPlotTitleDialog = new TitleDialog (this, "Plot");
			setCursor (Cursor.getDefaultCursor());
			}
		if (myPlot.formatPlotTitle (myPlotTitleDialog))
			{
			setTitle (myPlot.getTitle());
			repaint();
			}
		}

	/**
	 * "Margins..." menu item.
	 */
	private void doMargins()
		{
		// Create plot margins dialog if necessary.
		if (myPlotMarginsDialog == null)
			{
			setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
			myPlotMarginsDialog = new PlotMarginsDialog (this);
			setCursor (Cursor.getDefaultCursor());
			}
		if (myPlot.formatPlotMargins (myPlotMarginsDialog))
			{
			repaint();
			}
		}

	/**
	 * "Plot Area..." menu item.
	 */
	private void doPlotArea()
		{
		// Create plot area dialog if necessary.
		if (myPlotAreaDialog == null)
			{
			setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
			myPlotAreaDialog = new PlotAreaDialog (this);
			setCursor (Cursor.getDefaultCursor());
			}
		if (myPlot.formatPlotArea (myPlotAreaDialog))
			{
			repaint();
			}
		}

	/**
	 * "X Axis..." menu item.
	 */
	private void doXAxis()
		{
		// Create X axis dialog if necessary.
		if (myXAxisDialog == null)
			{
			setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
			myXAxisDialog = new AxisDialog (this, "X");
			setCursor (Cursor.getDefaultCursor());
			}
		if (myPlot.formatXAxis (myXAxisDialog))
			{
			repaint();
			}
		}

	/**
	 * "X Axis Title..." menu item.
	 */
	private void doXAxisTitle()
		{
		// Create X axis title dialog if necessary.
		if (myXAxisTitleDialog == null)
			{
			setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
			myXAxisTitleDialog = new TitleDialog (this, "X Axis");
			setCursor (Cursor.getDefaultCursor());
			}
		if (myPlot.formatXAxisTitle (myXAxisTitleDialog))
			{
			repaint();
			}
		}

	/**
	 * "Y Axis..." menu item.
	 */
	private void doYAxis()
		{
		// Create Y axis dialog if necessary.
		if (myYAxisDialog == null)
			{
			setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
			myYAxisDialog = new AxisDialog (this, "Y");
			setCursor (Cursor.getDefaultCursor());
			}
		if (myPlot.formatYAxis (myYAxisDialog))
			{
			repaint();
			}
		}

	/**
	 * "Y Axis Title..." menu item.
	 */
	private void doYAxisTitle()
		{
		// Create Y axis title dialog if necessary.
		if (myYAxisTitleDialog == null)
			{
			setCursor (Cursor.getPredefinedCursor (Cursor.WAIT_CURSOR));
			myYAxisTitleDialog = new TitleDialog (this, "Y Axis");
			setCursor (Cursor.getDefaultCursor());
			}
		if (myPlot.formatYAxisTitle (myYAxisTitleDialog))
			{
			repaint();
			}
		}

	/**
	 * "Close" or "Quit" menu item.
	 */
	private void doClose()
		{
		if (myPlotTitleDialog != null)
			{
			myPlotTitleDialog.dispose();
			myPlotTitleDialog = null;
			}
		if (myPlotMarginsDialog != null)
			{
			myPlotMarginsDialog.dispose();
			myPlotMarginsDialog = null;
			}
		if (myPlotAreaDialog != null)
			{
			myPlotAreaDialog.dispose();
			myPlotAreaDialog = null;
			}
		if (myXAxisTitleDialog != null)
			{
			myXAxisTitleDialog.dispose();
			myXAxisTitleDialog = null;
			}
		if (myYAxisTitleDialog != null)
			{
			myYAxisTitleDialog.dispose();
			myYAxisTitleDialog = null;
			}
		}

	}
