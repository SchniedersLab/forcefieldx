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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;

import javax.swing.JPanel;

/**
 * The GraphicsPanel class contains the 3D Canvas and its status box.
 */
public class GraphicsPanel extends JPanel {
	private GraphicsCanvas graphics;
	private JPanel canvasPanel;
	private JPanel statusPanel;

	public GraphicsPanel(GraphicsCanvas g, JPanel s) {
		super();
		setLayout(new BorderLayout());
		statusPanel = s;
		if (g != null) {
			canvasPanel = new JPanel(new GridLayout(1, 1));
			canvasPanel.add(g);
			add(canvasPanel, BorderLayout.CENTER);
		} else {
			setBackground(Color.BLACK);
		}
		add(statusPanel, BorderLayout.SOUTH);
	}

	public void setVisible(boolean v) {
		super.setVisible(v);
		if (graphics != null) {
			graphics.setVisible(v);
		}
	}
}
