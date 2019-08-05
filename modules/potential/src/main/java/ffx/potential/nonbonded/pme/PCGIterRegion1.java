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
import ffx.potential.nonbonded.ParticleMeshEwaldCart;

public class PCGIterRegion1 extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(PCGIterRegion1.class.getName());

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
    private final SharedDouble sumShared;
    private final SharedDouble sumCRShared;

    public PCGIterRegion1(int nt) {
        iterLoop1 = new PCGIterLoop1[nt];
        iterLoop2 = new PCGIterLoop2[nt];
        dotShared = new SharedDouble();
        dotCRShared = new SharedDouble();
        sumShared = new SharedDouble();
        sumCRShared = new SharedDouble();
    }

    public double getSum() {
        return sumShared.get();
    }

    public double getSumCR() {
        return sumCRShared.get();
    }

    public void init(Atom[] atoms, double[] polarizability,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     double[][][] field, double[][][] fieldCR,
                     ParticleMeshEwaldCart.PCGVariables pcgVectors) {
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
    public void start() {
        dotShared.set(0.0);
        dotCRShared.set(0.0);
        sumShared.set(0.0);
        sumCRShared.set(0.0);
    }

    @Override
    public void run() throws Exception {
        try {
            int ti = getThreadIndex();
            int nAtoms = atoms.length;
            if (iterLoop1[ti] == null) {
                iterLoop1[ti] = new PCGIterLoop1();
                iterLoop2[ti] = new PCGIterLoop2();
            }
            execute(0, nAtoms - 1, iterLoop1[ti]);
            if (ti == 0) {
                if (dotShared.get() != 0.0) {
                    dotShared.set(sumShared.get() / dotShared.get());
                }
                if (dotCRShared.get() != 0.0) {
                    dotCRShared.set(sumCRShared.get() / dotCRShared.get());
                }
            }
            barrier();
            execute(0, nAtoms - 1, iterLoop2[ti]);
        } catch (Exception e) {
            String message = "Fatal exception computing the mutual induced dipoles in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }

    }

    private class PCGIterLoop1 extends IntegerForLoop {

        public double dot;
        public double dotCR;
        public double sum;
        public double sumCR;

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void start() {
            dot = 0.0;
            dotCR = 0.0;
            sum = 0.0;
            sumCR = 0.0;
        }

        @Override
        public void finish() {
            dotShared.addAndGet(dot);
            dotCRShared.addAndGet(dotCR);
            sumShared.addAndGet(sum);
            sumCRShared.addAndGet(sumCR);
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            for (int i = lb; i <= ub; i++) {
                if (polarizability[i] > 0) {
                    double ipolar = 1.0 / polarizability[i];
                    inducedDipole[0][i][0] = vec[0][i];
                    inducedDipole[0][i][1] = vec[1][i];
                    inducedDipole[0][i][2] = vec[2][i];
                    vec[0][i] = conj[0][i] * ipolar - field[0][0][i];
                    vec[1][i] = conj[1][i] * ipolar - field[0][1][i];
                    vec[2][i] = conj[2][i] * ipolar - field[0][2][i];
                    inducedDipoleCR[0][i][0] = vecCR[0][i];
                    inducedDipoleCR[0][i][1] = vecCR[1][i];
                    inducedDipoleCR[0][i][2] = vecCR[2][i];
                    vecCR[0][i] = conjCR[0][i] * ipolar - fieldCR[0][0][i];
                    vecCR[1][i] = conjCR[1][i] * ipolar - fieldCR[0][1][i];
                    vecCR[2][i] = conjCR[2][i] * ipolar - fieldCR[0][2][i];
                } else {
                    inducedDipole[0][i][0] = 0.0;
                    inducedDipole[0][i][1] = 0.0;
                    inducedDipole[0][i][2] = 0.0;
                    vec[0][i] = 0.0;
                    vec[1][i] = 0.0;
                    vec[2][i] = 0.0;
                    inducedDipoleCR[0][i][0] = 0.0;
                    inducedDipoleCR[0][i][1] = 0.0;
                    inducedDipoleCR[0][i][2] = 0.0;
                    vecCR[0][i] = 0.0;
                    vecCR[1][i] = 0.0;
                    vecCR[2][i] = 0.0;
                }

                // Compute dot product of the conjugate vector and new residual.
                dot += conj[0][i] * vec[0][i]
                        + conj[1][i] * vec[1][i]
                        + conj[2][i] * vec[2][i];
                dotCR += conjCR[0][i] * vecCR[0][i]
                        + conjCR[1][i] * vecCR[1][i]
                        + conjCR[2][i] * vecCR[2][i];
                // Compute dot product of the previous residual and preconditioner.
                sum += rsd[0][i] * rsdPre[0][i]
                        + rsd[1][i] * rsdPre[1][i]
                        + rsd[2][i] * rsdPre[2][i];
                sumCR += rsdCR[0][i] * rsdPreCR[0][i]
                        + rsdCR[1][i] * rsdPreCR[1][i]
                        + rsdCR[2][i] * rsdPreCR[2][i];
            }

        }
    }

    private class PCGIterLoop2 extends IntegerForLoop {

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            double dot = dotShared.get();
            double dotCR = dotCRShared.get();
            for (int i = lb; i <= ub; i++) {
                    /*
                      Reduce the residual field, add to the induced dipoles
                      based on the scaled conjugate vector and finally set the
                      induced dipoles to the polarizability times the residual
                      field.
                     */
                rsd[0][i] -= dot * vec[0][i];
                rsd[1][i] -= dot * vec[1][i];
                rsd[2][i] -= dot * vec[2][i];
                rsdCR[0][i] -= dotCR * vecCR[0][i];
                rsdCR[1][i] -= dotCR * vecCR[1][i];
                rsdCR[2][i] -= dotCR * vecCR[2][i];
                vec[0][i] = inducedDipole[0][i][0] + dot * conj[0][i];
                vec[1][i] = inducedDipole[0][i][1] + dot * conj[1][i];
                vec[2][i] = inducedDipole[0][i][2] + dot * conj[2][i];
                vecCR[0][i] = inducedDipoleCR[0][i][0] + dotCR * conjCR[0][i];
                vecCR[1][i] = inducedDipoleCR[0][i][1] + dotCR * conjCR[1][i];
                vecCR[2][i] = inducedDipoleCR[0][i][2] + dotCR * conjCR[2][i];
                double polar = polarizability[i];
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
