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
 */
public class GraphicsSplitPane extends JSplitPane implements MouseListener,
		MouseMotionListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	boolean mouseClicked = false;
	int currentPos = 0;

	public GraphicsSplitPane() {
		super();
		addMouseListener(this);
		addMouseMotionListener(this);
	}

	public GraphicsSplitPane(int orient, boolean b, Component left,
			Component right) {
		super(orient, b, left, right);
		addMouseListener(this);
		addMouseMotionListener(this);
	}

	public void mouseClicked(MouseEvent e) {
		mouseClicked = true;
		currentPos = e.getX();
		// System.out.println("Split Clicked");
	}

	public void mouseDragged(MouseEvent e) {
		// System.out.println("Split Dragged");
		if (mouseClicked) {
			int change = e.getX() - currentPos;
			setDividerLocation(getDividerLocation() + change);
		}
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	public void mouseMoved(MouseEvent e) {
	}

	public void mousePressed(MouseEvent e) {
	}

	public void mouseReleased(MouseEvent e) {
		mouseClicked = false;
	}
}
