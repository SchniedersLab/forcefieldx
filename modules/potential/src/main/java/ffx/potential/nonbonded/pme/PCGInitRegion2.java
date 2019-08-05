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
import ffx.potential.nonbonded.ParticleMeshEwaldCart;

public class PCGInitRegion2 extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(PCGInitRegion2.class.getName());

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
    private double[][] rsdPre;
    private double[][] rsdPreCR;
    private double[][] conj;
    private double[][] conjCR;
    private double[][] vec;
    private double[][] vecCR;

    private final PCGInitLoop[] pcgLoop;

    public PCGInitRegion2(int nt) {
        pcgLoop = new PCGInitLoop[nt];
    }

    public void init(Atom[] atoms, double[] polarizability,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     double[][][] field, double[][][] fieldCR, ParticleMeshEwaldCart.PCGVariables pcgVectors) {
        this.atoms = atoms;
        this.polarizability = polarizability;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.field = field;
        this.fieldCR = fieldCR;
        this.rsd = pcgVectors.rsd;
        this.rsdCR = pcgVectors.rsdCR;
        this.rsdPre = pcgVectors.rsdPre;
        this.rsdPreCR = pcgVectors.rsdPreCR;
        this.conj = pcgVectors.conj;
        this.conjCR = pcgVectors.conjCR;
        this.vec = pcgVectors.vec;
        this.vecCR = pcgVectors.vecCR;
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

                // Revert to the stored induce dipoles.
                inducedDipole[0][i][0] = vec[0][i];
                inducedDipole[0][i][1] = vec[1][i];
                inducedDipole[0][i][2] = vec[2][i];
                inducedDipoleCR[0][i][0] = vecCR[0][i];
                inducedDipoleCR[0][i][1] = vecCR[1][i];
                inducedDipoleCR[0][i][2] = vecCR[2][i];

                // Set initial conjugate vector (induced dipoles).
                double udiag = 2.0;
                double polar = polarizability[i];
                rsdPre[0][i] = polar * (field[0][0][i] + udiag * rsd[0][i]);
                rsdPre[1][i] = polar * (field[0][1][i] + udiag * rsd[1][i]);
                rsdPre[2][i] = polar * (field[0][2][i] + udiag * rsd[2][i]);
                rsdPreCR[0][i] = polar * (fieldCR[0][0][i] + udiag * rsdCR[0][i]);
                rsdPreCR[1][i] = polar * (fieldCR[0][1][i] + udiag * rsdCR[1][i]);
                rsdPreCR[2][i] = polar * (fieldCR[0][2][i] + udiag * rsdCR[2][i]);
                conj[0][i] = rsdPre[0][i];
                conj[1][i] = rsdPre[1][i];
                conj[2][i] = rsdPre[2][i];
                conjCR[0][i] = rsdPreCR[0][i];
                conjCR[1][i] = rsdPreCR[1][i];
                conjCR[2][i] = rsdPreCR[2][i];
            }
        }
    }
}
