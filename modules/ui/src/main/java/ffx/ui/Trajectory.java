/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.ui;

import javax.swing.Timer;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import ffx.potential.MolecularAssembly;

/**
 * The Trajectory class controls playback of a TINKER trajectory.
 *
 * @author Michael J. Schnieders
 *
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

    /**
     * <p>
     * Constructor for Trajectory.</p>
     *
     * @param mol a {@link ffx.potential.MolecularAssembly} object.
     * @param f a {@link ffx.ui.MainPanel} object.
     */
    public Trajectory(MolecularAssembly mol, MainPanel f) {
        molecularSystem = mol;
        mainPanel = f;
        timer = new Timer(delay, this);
        timer.setCoalesce(true);
    }

    /**
     * {@inheritDoc}
     */
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

    /**
     * <p>
     * back</p>
     */
    public void back() {
        setFrame(getFrame() - 1);
    }

    /**
     * <p>
     * forward</p>
     */
    public void forward() {
        setFrame(getFrame() + 1);
    }

    /**
     * <p>
     * getFrame</p>
     *
     * @return a int.
     */
    public int getFrame() {
        return molecularSystem.getCurrentCycle();
    }

    /**
     * <p>
     * getFSystem</p>
     *
     * @return a {@link ffx.potential.MolecularAssembly} object.
     */
    public MolecularAssembly getFSystem() {
        return molecularSystem;
    }

    /**
     * <p>
     * Getter for the field <code>oscillate</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getOscillate() {
        return oscillate;
    }

    /**
     * <p>
     * getRate</p>
     *
     * @return a int.
     */
    public int getRate() {
        return desiredspeed;
    }

    /**
     * <p>
     * Getter for the field <code>skip</code>.</p>
     *
     * @return a int.
     */
    public int getSkip() {
        return skip;
    }

    /**
     * <p>
     * rewind</p>
     */
    public void rewind() {
        setFrame(1);
    }

    /**
     * <p>
     * setFrame</p>
     *
     * @param f a int.
     */
    public void setFrame(int f) {
        if (molecularSystem != null) {
            molecularSystem.setCurrentCycle(f);
            mainPanel.getGraphics3D().updateScene(molecularSystem, true, false,
                    null, false, null);
            mainPanel.getHierarchy().updateStatus();
        }
    }

    /**
     * <p>
     * Setter for the field <code>oscillate</code>.</p>
     *
     * @param o a boolean.
     */
    public void setOscillate(boolean o) {
        oscillate = o;
    }

    /**
     * <p>
     * setRate</p>
     *
     * @param s a int.
     */
    public void setRate(int s) {
        if (s > 0 && s <= 100) {
            desiredspeed = s;
            delay = 1000 / s;
            timer.setDelay(delay);
        }
    }

    /**
     * <p>
     * Setter for the field <code>skip</code>.</p>
     *
     * @param s a int.
     */
    public void setSkip(int s) {
        if (s < 1) {
            return;
        }
        skip = s % molecularSystem.getAtomList().size();
    }

    /**
     * <p>
     * start</p>
     */
    public void start() {
        timer.start();
    }

    /**
     * <p>
     * stop</p>
     */
    public void stop() {
        timer.stop();
    }
}
