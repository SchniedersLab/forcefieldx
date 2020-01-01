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

import javax.swing.JPanel;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;

/**
 * The GraphicsPanel class contains the 3D Canvas and its status box.
 *
 * @author Michael J. Schnieders
 */
public class GraphicsPanel extends JPanel {

    private GraphicsCanvas graphicsCanvas;

    /**
     * <p>
     * Constructor for GraphicsPanel.</p>
     *
     * @param graphicsCanvas a {@link ffx.ui.GraphicsCanvas} object.
     * @param statusPanel    a {@link javax.swing.JPanel} object.
     */
    public GraphicsPanel(GraphicsCanvas graphicsCanvas, JPanel statusPanel) {
        super();
        setLayout(new BorderLayout());
        if (graphicsCanvas != null) {
            JPanel canvasPanel = new JPanel(new GridLayout(1, 1));
            canvasPanel.add(graphicsCanvas);
            add(canvasPanel, BorderLayout.CENTER);
        } else {
            setBackground(Color.BLACK);
        }
        add(statusPanel, BorderLayout.SOUTH);
        this.graphicsCanvas = graphicsCanvas;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setVisible(boolean v) {
        super.setVisible(v);
        if (graphicsCanvas != null) {
            graphicsCanvas.setVisible(v);
        }
    }
}
