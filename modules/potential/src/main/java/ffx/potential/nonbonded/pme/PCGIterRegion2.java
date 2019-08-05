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
import ffx.potential.nonbonded.ParticleMeshEwaldCart.PCGVariables;

public class PCGIterRegion2 extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(PCGIterRegion2.class.getName());

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

    private final PCGIterLoop1[] iterLoop1;
    private final PCGIterLoop2[] iterLoop2;
    private final SharedDouble dotShared;
    private final SharedDouble dotCRShared;
    private final SharedDouble epsShared;
    private final SharedDouble epsCRShared;
    public double sum;
    public double sumCR;

    public PCGIterRegion2(int nt) {
        iterLoop1 = new PCGIterLoop1[nt];
        iterLoop2 = new PCGIterLoop2[nt];
        dotShared = new SharedDouble();
        dotCRShared = new SharedDouble();
        epsShared = new SharedDouble();
        epsCRShared = new SharedDouble();
    }

    public void init(Atom[] atoms, double[] polarizability,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     double[][][] field, double[][][] fieldCR, PCGVariables pcgVariables) {
        this.atoms = atoms;
        this.polarizability = polarizability;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.field = field;
        this.fieldCR = fieldCR;
        this.rsd = pcgVariables.rsd;
        this.rsdCR = pcgVariables.rsdCR;
        this.rsdPre = pcgVariables.rsdPre;
        this.rsdPreCR = pcgVariables.rsdPreCR;
        this.conj = pcgVariables.conj;
        this.conjCR = pcgVariables.conjCR;
        this.vec = pcgVariables.vec;
        this.vecCR = pcgVariables.vecCR;
    }

    public double getEps() {
        return epsShared.get();
    }

    public double getEpsCR() {
        return epsCRShared.get();
    }

    @Override
    public void start() {
        dotShared.set(0.0);
        dotCRShared.set(0.0);
        epsShared.set(0.0);
        epsCRShared.set(0.0);
        if (sum == 0.0) {
            sum = 1.0;
        }
        if (sumCR == 0.0) {
            sumCR = 1.0;
        }
    }

    @Override
    public void run() throws Exception {
        try {
            int ti = getThreadIndex();
            if (iterLoop1[ti] == null) {
                iterLoop1[ti] = new PCGIterLoop1();
                iterLoop2[ti] = new PCGIterLoop2();
            }
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, iterLoop1[ti]);
            execute(0, nAtoms - 1, iterLoop2[ti]);
        } catch (Exception e) {
            String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }

    }

    private class PCGIterLoop1 extends IntegerForLoop {

        public double dot;
        public double dotCR;

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void start() {
            dot = 0.0;
            dotCR = 0.0;
        }

        @Override
        public void finish() {
            dotShared.addAndGet(dot / sum);
            dotCRShared.addAndGet(dotCR / sumCR);
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            double udiag = 2.0;
            for (int i = lb; i <= ub; i++) {

                // Revert the induced dipoles to the saved values.
                inducedDipole[0][i][0] = vec[0][i];
                inducedDipole[0][i][1] = vec[1][i];
                inducedDipole[0][i][2] = vec[2][i];
                inducedDipoleCR[0][i][0] = vecCR[0][i];
                inducedDipoleCR[0][i][1] = vecCR[1][i];
                inducedDipoleCR[0][i][2] = vecCR[2][i];

                // Compute the dot product of the residual and preconditioner.
                double polar = polarizability[i];
                rsdPre[0][i] = polar * (field[0][0][i] + udiag * rsd[0][i]);
                rsdPre[1][i] = polar * (field[0][1][i] + udiag * rsd[1][i]);
                rsdPre[2][i] = polar * (field[0][2][i] + udiag * rsd[2][i]);
                rsdPreCR[0][i] = polar * (fieldCR[0][0][i] + udiag * rsdCR[0][i]);
                rsdPreCR[1][i] = polar * (fieldCR[0][1][i] + udiag * rsdCR[1][i]);
                rsdPreCR[2][i] = polar * (fieldCR[0][2][i] + udiag * rsdCR[2][i]);
                dot += rsd[0][i] * rsdPre[0][i]
                        + rsd[1][i] * rsdPre[1][i]
                        + rsd[2][i] * rsdPre[2][i];
                dotCR += rsdCR[0][i] * rsdPreCR[0][i]
                        + rsdCR[1][i] * rsdPreCR[1][i]
                        + rsdCR[2][i] * rsdPreCR[2][i];
            }
        }
    }

    private class PCGIterLoop2 extends IntegerForLoop {

        public double eps;
        public double epsCR;

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void start() {
            eps = 0.0;
            epsCR = 0.0;
        }

        @Override
        public void finish() {
            epsShared.addAndGet(eps);
            epsCRShared.addAndGet(epsCR);
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            double dot = dotShared.get();
            double dotCR = dotCRShared.get();
            for (int i = lb; i <= ub; i++) {
                // Update the conjugate vector and sum the square of the residual field.
                conj[0][i] = rsdPre[0][i] + dot * conj[0][i];
                conj[1][i] = rsdPre[1][i] + dot * conj[1][i];
                conj[2][i] = rsdPre[2][i] + dot * conj[2][i];
                conjCR[0][i] = rsdPreCR[0][i] + dotCR * conjCR[0][i];
                conjCR[1][i] = rsdPreCR[1][i] + dotCR * conjCR[1][i];
                conjCR[2][i] = rsdPreCR[2][i] + dotCR * conjCR[2][i];
                eps += rsd[0][i] * rsd[0][i]
                        + rsd[1][i] * rsd[1][i]
                        + rsd[2][i] * rsd[2][i];
                epsCR += rsdCR[0][i] * rsdCR[0][i]
                        + rsdCR[1][i] * rsdCR[1][i]
                        + rsdCR[2][i] * rsdCR[2][i];
            }
        }
    }
}


