/**
 * <p>Title: Force Field X</p>
 * <p>Description: Force Field X is a Molecular Biophysics Environment</p>
 * <p>Copyright: Copyright (c) Michael J. Schnieders 2002-2009</p>
 *
 * @author Michael J. Schnieders
 * @version 0.1
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.ui;

import java.awt.Component;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import javax.swing.JSplitPane;

/**
 * The GraphicsSplitPane is an early attempt at working around issues caused by
 * the heavyweight Canvas3D inside a lightweight Swing SplitPane. Specifically,
 * users cannot drag the slider toward the heavyweight Canvas3D.
 *
 * @author Michael J. Schnieders
 * @version $Id: $
 */
public class GraphicsSplitPane extends JSplitPane implements MouseListener,
		MouseMotionListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	boolean mouseClicked = false;
	int currentPos = 0;

	/**
	 * <p>Constructor for GraphicsSplitPane.</p>
	 */
	public GraphicsSplitPane() {
		super();
		addMouseListener(this);
		addMouseMotionListener(this);
	}

	/**
	 * <p>Constructor for GraphicsSplitPane.</p>
	 *
	 * @param orient a int.
	 * @param b a boolean.
	 * @param left a {@link java.awt.Component} object.
	 * @param right a {@link java.awt.Component} object.
	 */
	public GraphicsSplitPane(int orient, boolean b, Component left,
			Component right) {
		super(orient, b, left, right);
		addMouseListener(this);
		addMouseMotionListener(this);
	}

	/** {@inheritDoc} */
	public void mouseClicked(MouseEvent e) {
		mouseClicked = true;
		currentPos = e.getX();
		// System.out.println("Split Clicked");
	}

	/** {@inheritDoc} */
	public void mouseDragged(MouseEvent e) {
		// System.out.println("Split Dragged");
		if (mouseClicked) {
			int change = e.getX() - currentPos;
			setDividerLocation(getDividerLocation() + change);
		}
	}

	/** {@inheritDoc} */
	public void mouseEntered(MouseEvent e) {
	}

	/** {@inheritDoc} */
	public void mouseExited(MouseEvent e) {
	}

	/** {@inheritDoc} */
	public void mouseMoved(MouseEvent e) {
	}

	/** {@inheritDoc} */
	public void mousePressed(MouseEvent e) {
	}

	/** {@inheritDoc} */
	public void mouseReleased(MouseEvent e) {
		mouseClicked = false;
	}
}
