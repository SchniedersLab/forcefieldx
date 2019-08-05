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

import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.PCGVariables;

public class PCGInitRegion1 extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(PCGInitRegion1.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    private double[] polarizability;
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public double[][][] inducedDipole;
    public double[][][] inducedDipoleCR;
    /**
     * Direct induced dipoles.
     */
    public double[][] directDipole;
    public double[][] directDipoleCR;
    /**
     * Field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] field;
    /**
     * Chain rule field array for each thread. [threadID][X/Y/Z][atomID]
     */
    private double[][][] fieldCR;
    /**
     * PCG Variables.
     */
    private double[][] rsd;
    private double[][] rsdCR;
    private double[][] vec;
    private double[][] vecCR;

    private final PCGInitLoop[] pcgLoop;

    public PCGInitRegion1(int nt) {
        pcgLoop = new PCGInitLoop[nt];
    }

    public void init(Atom[] atoms, double[] polarizability,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     double[][] directDipole, double[][] directDipoleCR,
                     double[][][] field, double[][][] fieldCR,
                     PCGVariables pcgVariables) {
        this.atoms = atoms;
        this.polarizability = polarizability;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.directDipole = directDipole;
        this.directDipoleCR = directDipoleCR;
        this.field = field;
        this.fieldCR = fieldCR;
        this.rsd = pcgVariables.rsd;
        this.rsdCR = pcgVariables.rsdCR;
        this.vec = pcgVariables.vec;
        this.vecCR = pcgVariables.vecCR;
    }

    @Override
    public void run() throws Exception {
        try {
            int ti = getThreadIndex();
            if (pcgLoop[ti] == null) {
                pcgLoop[ti] = new PCGInitLoop();
            }
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, pcgLoop[ti]);
        } catch (Exception e) {
            String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }

    }

    private class PCGInitLoop extends IntegerForLoop {

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void run(int lb, int ub) throws Exception {

            for (int i = lb; i <= ub; i++) {
                // Set initial conjugate gradient residual (a field).
                double ipolar;
                if (polarizability[i] > 0) {
                    ipolar = 1.0 / polarizability[i];
                    rsd[0][i] = (directDipole[i][0] - inducedDipole[0][i][0]) * ipolar + field[0][0][i];
                    rsd[1][i] = (directDipole[i][1] - inducedDipole[0][i][1]) * ipolar + field[0][1][i];
                    rsd[2][i] = (directDipole[i][2] - inducedDipole[0][i][2]) * ipolar + field[0][2][i];
                    rsdCR[0][i] = (directDipoleCR[i][0] - inducedDipoleCR[0][i][0]) * ipolar + fieldCR[0][0][i];
                    rsdCR[1][i] = (directDipoleCR[i][1] - inducedDipoleCR[0][i][1]) * ipolar + fieldCR[0][1][i];
                    rsdCR[2][i] = (directDipoleCR[i][2] - inducedDipoleCR[0][i][2]) * ipolar + fieldCR[0][2][i];
                } else {
                    rsd[0][i] = 0.0;
                    rsd[1][i] = 0.0;
                    rsd[2][i] = 0.0;
                    rsdCR[0][i] = 0.0;
                    rsdCR[1][i] = 0.0;
                    rsdCR[2][i] = 0.0;
                }
                // Store the current induced dipoles and load the residual induced dipole
                double polar = polarizability[i];
                vec[0][i] = inducedDipole[0][i][0];
                vec[1][i] = inducedDipole[0][i][1];
                vec[2][i] = inducedDipole[0][i][2];
                vecCR[0][i] = inducedDipoleCR[0][i][0];
                vecCR[1][i] = inducedDipoleCR[0][i][1];
                vecCR[2][i] = inducedDipoleCR[0][i][2];
                inducedDipole[0][i][0] = polar * rsd[0][i];
                inducedDipole[0][i][1] = polar * rsd[1][i];
                inducedDipole[0][i][2] = polar * rsd[2][i];
                inducedDipoleCR[0][i][0] = polar * rsdCR[0][i];
                inducedDipoleCR[0][i][1] = polar * rsdCR[1][i];
                inducedDipoleCR[0][i][2] = polar * rsdCR[2][i];
            }
        }
    }
}