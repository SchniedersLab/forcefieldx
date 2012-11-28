/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
package ffx.ui;

import java.awt.Frame;
import java.awt.Window;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

/**
 * The FullScreenWindow class controls full screen graphics.
 *
 * @author Michael J. Schnieders
 *
 */
public class GraphicsFullScreen extends Window implements KeyListener {

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private boolean fullScreen = false;

    /**
     * <p>Constructor for GraphicsFullScreen.</p>
     *
     * @param f a {@link java.awt.Frame} object.
     * @param graphics a {@link ffx.ui.GraphicsCanvas} object.
     */
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

    /**
     * <p>enterFullScreen</p>
     */
    public void enterFullScreen() {
        /*
         * if (fullScreen) { return; } fullScreenCanvas.startRenderer();
         * setVisible(true); fullScreenCanvas.requestFocus(); fullScreen = true;
         */
    }

    /**
     * <p>exitFullScreen</p>
     */
    public void exitFullScreen() {
        /*
         * if (!fullScreen) { return; } setVisible(false);
         * fullScreenCanvas.stopRenderer(); fullScreen = false;
         */
    }

    /**
     * {@inheritDoc}
     */
    public void keyPressed(KeyEvent evt) {
        if (evt.getKeyCode() == KeyEvent.VK_ESCAPE) {
            exitFullScreen();
        } else if (evt.getKeyChar() == 'e') {
            exitFullScreen();
        } else if (evt.getKeyChar() == 'x') {
            exitFullScreen();
        }
    }

    /**
     * {@inheritDoc}
     */
    public void keyReleased(KeyEvent evt) {
        keyPressed(evt);
    }

    /**
     * {@inheritDoc}
     */
    public void keyTyped(KeyEvent evt) {
        keyPressed(evt);
    }

    /**
     * <p>toggleFullScreen</p>
     */
    public void toggleFullScreen() {
        if (fullScreen) {
            exitFullScreen();
        } else {
            enterFullScreen();
        }
    }
}
