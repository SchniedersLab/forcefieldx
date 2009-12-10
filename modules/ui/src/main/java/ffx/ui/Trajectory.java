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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Timer;

import ffx.potential.bonded.MolecularAssembly;

/**
 * The Trajectory class controls playback of a TINKER trajectory.
 */
public class Trajectory implements ActionListener {
	private MolecularAssembly molecularSystem;
	private MainPanel mainPanel;
	private Timer timer;
	private int delay = 50;
	private int desiredspeed = 20;
	private int cycle = 1;
	private int sign = 1;
	private int skip = 1;
	private boolean oscillate = false;

	public Trajectory(MolecularAssembly mol, MainPanel f) {
		molecularSystem = mol;
		mainPanel = f;
		timer = new Timer(delay, this);
		timer.setCoalesce(true);
	}

	public void actionPerformed(ActionEvent e) {
		if (mainPanel.getGraphics3D().isCacheFull()) {
			return;
		}
		cycle = advance(skip * sign);
		setFrame(cycle);
	}

	private int advance(int adv) {
		if (molecularSystem != null) {
			cycle = molecularSystem.getCurrentCycle();
			int frame = cycle + adv;
			if ((frame) <= 0) {
				sign = 1;
				if (oscillate) {
					frame = -adv - cycle;
				} else {
					frame = molecularSystem.getCycles() + (adv + cycle);
				}
			} else if ((frame) > molecularSystem.getCycles()) {
				if (oscillate) {
					frame = molecularSystem.getCycles() + (-adv + cycle);
					sign = -1;
				} else {
					sign = 1;
					frame = cycle - molecularSystem.getCycles() + adv;
				}
			}
			return frame;
		}
		return 0;
	}

	public void back() {
		setFrame(getFrame() - 1);
	}

	public void forward() {
		setFrame(getFrame() + 1);
	}

	public int getFrame() {
		return molecularSystem.getCurrentCycle();
	}

	public MolecularAssembly getFSystem() {
		return molecularSystem;
	}

	public boolean getOscillate() {
		return oscillate;
	}

	public int getRate() {
		return desiredspeed;
	}

	public int getSkip() {
		return skip;
	}

	public void rewind() {
		setFrame(1);
	}

	public void setFrame(int f) {
		if (molecularSystem != null) {
			molecularSystem.setCurrentCycle(f);
			mainPanel.getGraphics3D().updateScene(molecularSystem, true, false,
					null, false, null);
			mainPanel.getHierarchy().updateStatus();
		}
	}

	public void setOscillate(boolean o) {
		oscillate = o;
	}

	public void setRate(int s) {
		if (s > 0 && s <= 100) {
			desiredspeed = s;
			delay = 1000 / s;
			timer.setDelay(delay);
		}
	}

	public void setSkip(int s) {
		if (s < 1) {
			return;
		}
		skip = s % molecularSystem.getAtomList().size();
	}

	public void start() {
		timer.start();
	}

	public void stop() {
		timer.stop();
	}
}
