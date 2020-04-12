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
package ffx.potential.nonbonded.pme;

import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.EwaldParameters;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t100;

/**
 * Parallel computation of the OPT induced dipoles.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class OPTRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(OPTRegion.class.getName());
    public final double[] optCoefficients;
    /**
     * Induced dipoles for extrapolated perturbation theory.
     */
    public final int optOrder = 2;
    private final OPTLoop[] optLoop;
    private final double[] optCoefficientsSum;
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public double[][][] inducedDipole;
    public double[][][] inducedDipoleCR;
    public double[][][] optDipole;
    public double[][][] optDipoleCR;
    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    private double[] polarizability;
    private double[][] cartesianDipolePhi;
    private double[][] cartesianDipolePhiCR;
    /**
     * Field array.
     */
    private AtomicDoubleArray3D field;
    /**
     * Chain rule field array.
     */
    private AtomicDoubleArray3D fieldCR;
    /**
     * Flag to indicate use of generalized Kirkwood.
     */
    private boolean generalizedKirkwoodTerm;
    private GeneralizedKirkwood generalizedKirkwood;
    private double aewald;
    private double aewald3;
    private int currentOptOrder;

    public OPTRegion(int nt) {
        optLoop = new OPTLoop[nt];
        optCoefficients = new double[optOrder + 1];
        optCoefficientsSum = new double[optOrder + 1];
        switch (optOrder) {
            case 1:
                optCoefficients[0] = 0.530;
                optCoefficients[1] = 0.604;
                break;
            case 2:
                optCoefficients[0] = 0.042;
                optCoefficients[1] = 0.635;
                optCoefficients[2] = 0.414;
                break;
            case 3:
                optCoefficients[0] = -0.132;
                optCoefficients[1] = 0.218;
                optCoefficients[2] = 0.637;
                optCoefficients[3] = 0.293;
                break;
            case 4:
                optCoefficients[0] = -0.071;
                optCoefficients[1] = -0.096;
                optCoefficients[2] = 0.358;
                optCoefficients[3] = 0.587;
                optCoefficients[4] = 0.216;
                break;
            case 5:
                optCoefficients[0] = -0.005;
                optCoefficients[1] = -0.129;
                optCoefficients[2] = -0.026;
                optCoefficients[3] = 0.465;
                optCoefficients[4] = 0.528;
                optCoefficients[5] = 0.161;
                break;
            case 6:
                optCoefficients[0] = 0.014;
                optCoefficients[1] = -0.041;
                optCoefficients[2] = -0.172;
                optCoefficients[3] = 0.073;
                optCoefficients[4] = 0.535;
                optCoefficients[5] = 0.467;
                optCoefficients[6] = 0.122;
                break;
            default:
                logger.severe(" Unsupported OPT order.");
        }

        for (int i = 0; i <= optOrder; i++) {
            for (int j = optOrder; j >= i; j--) {
                optCoefficientsSum[i] += optCoefficients[j];
            }
        }
    }

    public void init(int currentOptOrder, Atom[] atoms, double[] polarizability,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     double[][] cartesianDipolePhi, double[][] cartesianDipolePhiCR,
                     AtomicDoubleArray3D field, AtomicDoubleArray3D fieldCR,
                     boolean generalizedKirkwoodTerm, GeneralizedKirkwood generalizedKirkwood,
                     EwaldParameters ewaldParameters) {
        this.currentOptOrder = currentOptOrder;
        this.atoms = atoms;
        this.polarizability = polarizability;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.cartesianDipolePhi = cartesianDipolePhi;
        this.cartesianDipolePhiCR = cartesianDipolePhiCR;
        this.field = field;
        this.fieldCR = fieldCR;
        this.generalizedKirkwoodTerm = generalizedKirkwoodTerm;
        this.generalizedKirkwood = generalizedKirkwood;
        this.aewald = ewaldParameters.aewald;
        this.aewald3 = ewaldParameters.aewald3;
    }

    @Override
    public void run() throws Exception {
        try {
            int ti = getThreadIndex();
            if (optLoop[ti] == null) {
                optLoop[ti] = new OPTRegion.OPTLoop();
            }
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, optLoop[ti]);
        } catch (RuntimeException ex) {
            logger.warning("Fatal exception computing the opt induced dipoles in thread " + getThreadIndex());
            throw ex;
        } catch (Exception e) {
            String message = "Fatal exception computing the opt induced dipoles in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private class OPTLoop extends IntegerForLoop {

        @Override
        public void run(int lb, int ub) throws Exception {
            int threadID = getThreadIndex();
            final double[][] induced0 = inducedDipole[0];
            final double[][] inducedCR0 = inducedDipoleCR[0];

            if (aewald > 0.0) {
                // Add the self and reciprocal space fields to the real space field.
                for (int i = lb; i <= ub; i++) {
                    double[] dipolei = induced0[i];
                    double[] dipoleCRi = inducedCR0[i];
                    final double[] phii = cartesianDipolePhi[i];
                    final double[] phiCRi = cartesianDipolePhiCR[i];
                    double fx = aewald3 * dipolei[0] - phii[t100];
                    double fy = aewald3 * dipolei[1] - phii[t010];
                    double fz = aewald3 * dipolei[2] - phii[t001];
                    double fxCR = aewald3 * dipoleCRi[0] - phiCRi[t100];
                    double fyCR = aewald3 * dipoleCRi[1] - phiCRi[t010];
                    double fzCR = aewald3 * dipoleCRi[2] - phiCRi[t001];
                    field.add(threadID, i, fx, fy, fz);
                    fieldCR.add(threadID, i, fxCR, fyCR, fzCR);
                }
            }
            if (generalizedKirkwoodTerm) {
                AtomicDoubleArray3D fieldGK = generalizedKirkwood.getFieldGK();
                AtomicDoubleArray3D fieldGKCR = generalizedKirkwood.getFieldGKCR();
                // Add the GK reaction field to the intramolecular field.
                for (int i = lb; i <= ub; i++) {
                    field.add(threadID, i, fieldGK.getX(i), fieldGK.getY(i), fieldGK.getZ(i));
                    fieldCR.add(threadID, i, fieldGKCR.getX(i), fieldGKCR.getY(i), fieldGKCR.getZ(i));
                }
            }

            // Reduce the real space field.
            field.reduce(lb, ub);
            fieldCR.reduce(lb, ub);

            // Collect the current Opt Order induced dipole.
            for (int i = lb; i <= ub; i++) {
                final double[] ind = induced0[i];
                final double[] indCR = inducedCR0[i];
                final double polar = polarizability[i];
                for (int j = 0; j < 3; j++) {
                    optDipole[currentOptOrder][i][j] = polar * field.get(j, i);
                    optDipoleCR[currentOptOrder][i][j] = polar * fieldCR.get(j, i);
                    ind[j] = optDipole[currentOptOrder][i][j];
                    indCR[j] = optDipoleCR[currentOptOrder][i][j];
                }
            }
        }

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }
    }
}
