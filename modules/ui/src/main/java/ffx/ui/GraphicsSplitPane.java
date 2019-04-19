//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

import javax.swing.JSplitPane;

import java.awt.Component;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * The GraphicsSplitPane is an early attempt at working around issues caused by
 * the heavyweight Canvas3D inside a lightweight Swing SplitPane. Specifically,
 * users cannot drag the slider toward the heavyweight Canvas3D.
 *
 * @author Michael J. Schnieders
 */
public class GraphicsSplitPane extends JSplitPane implements MouseListener,
        MouseMotionListener {


    private boolean mouseClicked = false;
    private int currentPos = 0;

    /**
     * <p>
     * Constructor for GraphicsSplitPane.</p>
     */
    public GraphicsSplitPane() {
        super();
        addMouseListener(this);
        addMouseMotionListener(this);
    }

    /**
     * <p>
     * Constructor for GraphicsSplitPane.</p>
     *
     * @param orient a int.
     * @param b      a boolean.
     * @param left   a {@link java.awt.Component} object.
     * @param right  a {@link java.awt.Component} object.
     */
    public GraphicsSplitPane(int orient, boolean b, Component left,
                             Component right) {
        super(orient, b, left, right);
        addMouseListener(this);
        addMouseMotionListener(this);
    }

    /**
     * {@inheritDoc}
     */
    public void mouseClicked(MouseEvent e) {
        mouseClicked = true;
        currentPos = e.getX();
    }

    /**
     * {@inheritDoc}
     */
    public void mouseDragged(MouseEvent e) {
        if (mouseClicked) {
            int change = e.getX() - currentPos;
            setDividerLocation(getDividerLocation() + change);
        }
    }

    /**
     * {@inheritDoc}
     */
    public void mouseEntered(MouseEvent e) {
    }

    /**
     * {@inheritDoc}
     */
    public void mouseExited(MouseEvent e) {
    }

    /**
     * {@inheritDoc}
     */
    public void mouseMoved(MouseEvent e) {
    }

    /**
     * {@inheritDoc}
     */
    public void mousePressed(MouseEvent e) {
    }

    /**
     * {@inheritDoc}
     */
    public void mouseReleased(MouseEvent e) {
        mouseClicked = false;
    }
}
