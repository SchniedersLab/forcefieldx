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
package ffx.potential.nonbonded.pme;

import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;

import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import static ffx.potential.nonbonded.ParticleMeshEwald.DEFAULT_ELECTRIC;

/**
 * Parallel computation of the polarization energy as sum over atomic contributions (-1/2 u.E).
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PolarizationEnergyRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(PolarizationEnergyRegion.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    private double[] polarizability;
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    private double[][][] inducedDipole;
    /**
     * Direct induced dipoles.
     */
    private double[][] directDipoleCR;
    private double polarizationScale;

    private final double electric;
    private final PolarizationEnergyLoop[] polarizationLoop;
    private final SharedDouble polarizationEnergy = new SharedDouble();

    public PolarizationEnergyRegion(int nt, ForceField forceField) {
        electric = forceField.getDouble("ELECTRIC", DEFAULT_ELECTRIC);
        polarizationLoop = new PolarizationEnergyLoop[nt];
    }

    public void init(Atom[] atoms, double[] polarizability,
                     double[][][] inducedDipole, double[][] directDipoleCR,
                     double polarizationScale) {
        this.atoms = atoms;
        this.polarizability = polarizability;
        this.inducedDipole = inducedDipole;
        this.directDipoleCR = directDipoleCR;
        this.polarizationScale = polarizationScale;
    }

    @Override
    public void start() {
        polarizationEnergy.set(0.0);
    }

    /**
     * Return the final polarization energy.
     *
     * @return The polarization energy.
     */
    public double getPolarizationEnergy() {
        return polarizationEnergy.get();
    }

    /**
     * Set the current polarization energy.
     *
     * @param energy Value to set the polarization energy to.
     */
    public void setPolarizationEnergy(double energy) {
        polarizationEnergy.set(energy);
    }

    @Override
    public void run() throws Exception {
        try {
            int ti = getThreadIndex();
            if (polarizationLoop[ti] == null) {
                polarizationLoop[ti] = new PolarizationEnergyLoop();
            }
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, polarizationLoop[ti]);
        } catch (RuntimeException ex) {
            logger.warning("Runtime exception computing polarization energy in thread " + getThreadIndex());
            throw ex;
        } catch (Exception e) {
            String message = "Fatal exception computing polarization energy in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    @Override
    public void finish() {
        double energy = polarizationEnergy.get();
        polarizationEnergy.set(energy * polarizationScale * electric);
    }

    private class PolarizationEnergyLoop extends IntegerForLoop {

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            double energy = 0.0;
            for (int i = lb; i <= ub; i++) {
                if (polarizability[i] == 0.0) {
                    continue;
                }
                double uix = inducedDipole[0][i][0];
                double uiy = inducedDipole[0][i][1];
                double uiz = inducedDipole[0][i][2];
                double pix = directDipoleCR[i][0];
                double piy = directDipoleCR[i][1];
                double piz = directDipoleCR[i][2];
                energy += (uix * pix + uiy * piy + uiz * piz) / polarizability[i];
            }
            polarizationEnergy.addAndGet(-0.5 * energy);
        }
    }
}