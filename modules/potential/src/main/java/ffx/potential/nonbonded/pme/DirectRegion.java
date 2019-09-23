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
import edu.rit.pj.ParallelTeam;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.EwaldParameters;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t100;

/**
 * Parallel computation of induced dipoles due to the direct field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class DirectRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(DirectRegion.class.getName());

    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    private double[] polarizability;
    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    private double[][][] globalMultipole;
    private double[][] cartMultipolePhi;
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
    private final DirectLoop[] directLoop;

    public DirectRegion(int nt) {
        directLoop = new DirectLoop[nt];
    }

    public void init(Atom[] atoms, double[] polarizability, double[][][] globalMultipole, double[][] cartMultipolePhi,
                     AtomicDoubleArray3D field, AtomicDoubleArray3D fieldCR,
                     boolean generalizedKirkwoodTerm, GeneralizedKirkwood generalizedKirkwood,
                     EwaldParameters ewaldParameters,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     double[][] directDipole, double[][] directDipoleCR) {
        // Input
        this.atoms = atoms;
        this.polarizability = polarizability;
        this.globalMultipole = globalMultipole;
        this.cartMultipolePhi = cartMultipolePhi;
        this.field = field;
        this.fieldCR = fieldCR;
        this.generalizedKirkwoodTerm = generalizedKirkwoodTerm;
        this.generalizedKirkwood = generalizedKirkwood;
        this.aewald = ewaldParameters.aewald;
        this.aewald3 = ewaldParameters.aewald3;
        // Output
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.directDipole = directDipole;
        this.directDipoleCR = directDipoleCR;
    }

    /**
     * Execute the DirectRegion with the passed ParallelTeam.
     *
     * @param parallelTeam The ParallelTeam instance to execute with.
     */
    public void executeWith(ParallelTeam parallelTeam) {
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = " Exception computing direct induced dipoles.\n";
            logger.log(Level.WARNING, message, e);
        }
    }

    @Override
    public void run() throws Exception {
        int ti = getThreadIndex();
        if (directLoop[ti] == null) {
            directLoop[ti] = new DirectLoop();
        }
        try {
            int nAtoms = atoms.length;
            execute(0, nAtoms - 1, directLoop[ti]);
        } catch (Exception e) {
            String message = "Fatal exception computing the direct induced dipoles in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private class DirectLoop extends IntegerForLoop {

        @Override
        public IntegerSchedule schedule() {
            return IntegerSchedule.fixed();
        }

        @Override
        public void run(int lb, int ub) throws Exception {
            int threadID = getThreadIndex();
            if (aewald > 0.0) {
                // Add the self and reciprocal space contributions.
                for (int i = lb; i <= ub; i++) {
                    double[] mpolei = globalMultipole[0][i];
                    double[] phii = cartMultipolePhi[i];
                    double fx = aewald3 * mpolei[t100] - phii[t100];
                    double fy = aewald3 * mpolei[t010] - phii[t010];
                    double fz = aewald3 * mpolei[t001] - phii[t001];
                    field.add(threadID, i, fx, fy, fz);
                    fieldCR.add(threadID, i, fx, fy, fz);
                }
            }
            if (generalizedKirkwoodTerm) {
                // Set the electric field to the direct field plus the permanent GK reaction field.
                AtomicDoubleArray3D fieldGK = generalizedKirkwood.getFieldGK();
                for (int i = lb; i <= ub; i++) {
                    double fx = fieldGK.getX(i);
                    double fy = fieldGK.getY(i);
                    double fz = fieldGK.getZ(i);
                    field.add(threadID, i, fx, fy, fz);
                    fieldCR.add(threadID, i, fx, fy, fz);
                }
            }
            // Reduce the direct field.
            field.reduce(lb, ub);
            fieldCR.reduce(lb, ub);

            // Set the direct induced dipoles to the polarizability multiplied by the direct field.
            final double[][] induced0 = inducedDipole[0];
            final double[][] inducedCR0 = inducedDipoleCR[0];
            for (int i = lb; i <= ub; i++) {
                final double polar = polarizability[i];
                final double[] ind = induced0[i];
                final double[] directi = directDipole[i];
                ind[0] = polar * field.getX(i);
                ind[1] = polar * field.getY(i);
                ind[2] = polar * field.getZ(i);
                directi[0] = ind[0];
                directi[1] = ind[1];
                directi[2] = ind[2];
                final double[] indCR = inducedCR0[i];
                final double[] directCRi = directDipoleCR[i];
                indCR[0] = polar * fieldCR.getX(i);
                indCR[1] = polar * fieldCR.getY(i);
                indCR[2] = polar * fieldCR.getZ(i);
                directCRi[0] = indCR[0];
                directCRi[1] = indCR[1];
                directCRi[2] = indCR[2];
            }
        }
    }
}