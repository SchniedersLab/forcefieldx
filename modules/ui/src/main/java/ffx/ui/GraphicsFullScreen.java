//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.ui;

import java.awt.Frame;
import java.awt.Window;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

/**
 * The FullScreenWindow class controls full screen graphics.
 *
 * @author Michael J. Schnieders
 */
public class GraphicsFullScreen extends Window implements KeyListener {

    private boolean fullScreen = false;

    /**
     * <p>
     * Constructor for GraphicsFullScreen.</p>
     *
     * @param f        a {@link java.awt.Frame} object.
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
     * <p>
     * toggleFullScreen</p>
     */
    public void toggleFullScreen() {
        if (fullScreen) {
            exitFullScreen();
        } else {
            enterFullScreen();
        }
    }

    /**
     * <p>
     * enterFullScreen</p>
     */
    void enterFullScreen() {
        /*
         * if (fullScreen) { return; } fullScreenCanvas.startRenderer();
         * setVisible(true); fullScreenCanvas.requestFocus(); fullScreen = true;
         */
    }

    /**
     * <p>
     * exitFullScreen</p>
     */
    private void exitFullScreen() {
        /*
         * if (!fullScreen) { return; } setVisible(false);
         * fullScreenCanvas.stopRenderer(); fullScreen = false;
         */
    }
}
