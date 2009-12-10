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

import java.awt.Frame;
import java.awt.Window;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

/**
 * The FullScreenWindow class controls full screen graphics.
 */
public class GraphicsFullScreen extends Window implements KeyListener {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private boolean fullScreen = false;

	public GraphicsFullScreen(Frame f, GraphicsCanvas graphics) {
		super(f);
		/*
		 * setLayout(new BorderLayout()); screenSize =
		 * Toolkit.getDefaultToolkit().getScreenSize(); setSize(screenSize);
		 * fullScreenCanvas = new Canvas3D(graphics.getGraphicsConfiguration());
		 * fullScreenCanvas.stopRenderer();
		 * graphics.getView().addCanvas3D(fullScreenCanvas);
		 * addKeyListener(this); fullScreenCanvas.addKeyListener(this);
		 * setFocusable(true); fullScreenCanvas.setFocusable(true);
		 * add(fullScreenCanvas, BorderLayout.CENTER);
		 */
	}

	public void enterFullScreen() {
		/*
		 * if (fullScreen) { return; } fullScreenCanvas.startRenderer();
		 * setVisible(true); fullScreenCanvas.requestFocus(); fullScreen = true;
		 */
	}

	public void exitFullScreen() {
		/*
		 * if (!fullScreen) { return; } setVisible(false);
		 * fullScreenCanvas.stopRenderer(); fullScreen = false;
		 */
	}

	public void keyPressed(KeyEvent evt) {
		if (evt.getKeyCode() == KeyEvent.VK_ESCAPE) {
			exitFullScreen();
		} else if (evt.getKeyChar() == 'e') {
			exitFullScreen();
		} else if (evt.getKeyChar() == 'x') {
			exitFullScreen();
		}
	}

	public void keyReleased(KeyEvent evt) {
		keyPressed(evt);
	}

	public void keyTyped(KeyEvent evt) {
		keyPressed(evt);
	}

	public void toggleFullScreen() {
		if (fullScreen) {
			exitFullScreen();
		} else {
			enterFullScreen();
		}
	}
}
