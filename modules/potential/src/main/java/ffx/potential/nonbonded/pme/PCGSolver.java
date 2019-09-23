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

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwaldCart;
import ffx.potential.nonbonded.ParticleMeshEwaldCart.EwaldParameters;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.utils.EnergyException;
import static ffx.numerics.special.Erf.erfc;

/**
 * Parallel pre-conditioned conjugate gradient solver for the self-consistent field.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PCGSolver {

    private static final Logger logger = Logger.getLogger(PCGSolver.class.getName());

    private final InducedDipolePreconditionerRegion inducedDipolePreconditionerRegion;
    private final PCGInitRegion1 pcgInitRegion1;
    private final PCGInitRegion2 pcgInitRegion2;
    private final PCGIterRegion1 pcgIterRegion1;
    private final PCGIterRegion2 pcgIterRegion2;
    private final double poleps;

    /**
     * Residual vector.
     */
    private double[][] rsd;
    /**
     * Residual vector for the chain-rule dipoles.
     */
    private double[][] rsdCR;
    /**
     * Preconditioner residual.
     */
    private double[][] rsdPre;
    /**
     * Preconditioner residual for the chain-rule dipoles.
     */
    private double[][] rsdPreCR;
    /**
     * Conjugate search direction.
     */
    private double[][] conj;
    /**
     * Conjugate search direction for the chain-rule dipoles.
     */
    private double[][] conjCR;
    /**
     * Work vector.
     */
    private double[][] vec;
    /**
     * Work vector for the chain-rule dipoles.
     */
    private double[][] vecCR;

    /**
     * Neighbor lists, without atoms beyond the preconditioner cutoff.
     * [nSymm][nAtoms][nIncludedNeighbors]
     */
    int[][][] preconditionerLists;
    /**
     * Number of neighboring atoms within the preconditioner cutoff.
     * [nSymm][nAtoms]
     */
    int[][] preconditionerCounts;
    public double preconditionerCutoff;
    public double preconditionerEwald = 0.0;
    private final int preconditionerListSize = 50;

    /**
     * An ordered array of atoms in the system.
     */
    private Atom[] atoms;
    /**
     * Dimensions of [nsymm][xyz][nAtoms].
     */
    private double[][][] coordinates;
    private double[] polarizability;
    private double[] ipdamp;
    private double[] thole;
    /**
     * When computing the polarization energy at Lambda there are 3 pieces.
     * <p>
     * 1.) Upol(1) = The polarization energy computed normally (ie. system with
     * ligand).
     * <p>
     * 2.) Uenv = The polarization energy of the system without the ligand.
     * <p>
     * 3.) Uligand = The polarization energy of the ligand by itself.
     * <p>
     * Upol(L) = L*Upol(1) + (1-L)*(Uenv + Uligand)
     * <p>
     * Set the "use" array to true for all atoms for part 1. Set the "use" array
     * to true for all atoms except the ligand for part 2. Set the "use" array
     * to true only for the ligand atoms for part 3.
     * <p>
     * The "use" array can also be employed to turn off atoms for computing the
     * electrostatic energy of sub-structures.
     */
    private boolean[] use;
    /**
     * Unit cell and spacegroup information.
     */
    private Crystal crystal;
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    private double[][][] inducedDipole;
    private double[][][] inducedDipoleCR;
    /**
     * Direct induced dipoles.
     */
    private double[][] directDipole;
    private double[][] directDipoleCR;

    /**
     * Field array.
     */
    private AtomicDoubleArray3D field;
    /**
     * Chain rule field array.
     */
    private AtomicDoubleArray3D fieldCR;

    private EwaldParameters ewaldParameters;
    /**
     * The default ParallelTeam encapsulates the maximum number of threads used
     * to parallelize the electrostatics calculation.
     */
    private ParallelTeam parallelTeam;
    /**
     * Pairwise schedule for load balancing.
     */
    private IntegerSchedule realSpaceSchedule;
    private long[] realSpaceSCFTime;

    /**
     * Constructor the PCG solver.
     *
     * @param maxThreads Number of threads.
     * @param poleps     Convergence criteria (RMS Debye).
     * @param forceField Force field in use.
     * @param nAtoms     Initial number of atoms.
     */
    public PCGSolver(int maxThreads, double poleps, ForceField forceField, int nAtoms) {
        this.poleps = poleps;
        inducedDipolePreconditionerRegion = new InducedDipolePreconditionerRegion(maxThreads);
        pcgInitRegion1 = new PCGInitRegion1(maxThreads);
        pcgInitRegion2 = new PCGInitRegion2(maxThreads);
        pcgIterRegion1 = new PCGIterRegion1(maxThreads);
        pcgIterRegion2 = new PCGIterRegion2(maxThreads);

        // The size of the preconditioner neighbor list depends on the size of the preconditioner cutoff.
        boolean preconditioner = forceField.getBoolean(ForceField.ForceFieldBoolean.USE_SCF_PRECONDITIONER, true);
        if (preconditioner) {
            preconditionerCutoff = forceField.getDouble(ForceField.ForceFieldDouble.CG_PRECONDITIONER_CUTOFF, 4.5);
            preconditionerEwald = forceField.getDouble(ForceField.ForceFieldDouble.CG_PRECONDITIONER_EWALD, 0.0);
        } else {
            preconditionerCutoff = 0.0;
        }

        allocateVectors(nAtoms);
    }

    /**
     * Allocate PCG vectors.
     *
     * @param nAtoms The number of atoms.
     */
    public void allocateVectors(int nAtoms) {
        if (rsd == null || rsd[0].length != nAtoms) {
            rsd = new double[3][nAtoms];
            rsdCR = new double[3][nAtoms];
            rsdPre = new double[3][nAtoms];
            rsdPreCR = new double[3][nAtoms];
            conj = new double[3][nAtoms];
            conjCR = new double[3][nAtoms];
            vec = new double[3][nAtoms];
            vecCR = new double[3][nAtoms];
        }
    }

    /**
     * Allocate storage for pre-conditioner neighbor list.
     *
     * @param nSymm  Number of symmetry operators.
     * @param nAtoms Number of atoms.
     */
    public void allocateLists(int nSymm, int nAtoms) {
        preconditionerLists = new int[nSymm][nAtoms][preconditionerListSize];
        preconditionerCounts = new int[nSymm][nAtoms];
    }

    public void init(Atom[] atoms, double[][][] coordinates, double[] polarizability,
                     double[] ipdamp, double[] thole, boolean[] use, Crystal crystal,
                     double[][][] inducedDipole, double[][][] inducedDipoleCR,
                     double[][] directDipole, double[][] directDipoleCR,
                     AtomicDoubleArray3D field, AtomicDoubleArray3D fieldCR,
                     EwaldParameters ewaldParameters, ParallelTeam parallelTeam,
                     IntegerSchedule realSpaceSchedule, long[] realSpaceSCFTime) {
        this.atoms = atoms;
        this.coordinates = coordinates;
        this.polarizability = polarizability;
        this.ipdamp = ipdamp;
        this.thole = thole;
        this.use = use;
        this.crystal = crystal;
        this.inducedDipole = inducedDipole;
        this.inducedDipoleCR = inducedDipoleCR;
        this.directDipole = directDipole;
        this.directDipoleCR = directDipoleCR;
        this.field = field;
        this.fieldCR = fieldCR;
        this.ewaldParameters = ewaldParameters;
        this.parallelTeam = parallelTeam;
        this.realSpaceSchedule = realSpaceSchedule;
        this.realSpaceSCFTime = realSpaceSCFTime;
    }

    public int scfByPCG(boolean print, long startTime, ParticleMeshEwaldCart pme) {
        long directTime = System.nanoTime() - startTime;
        // A request of 0 SCF cycles simplifies mutual polarization to direct polarization.
        StringBuilder sb = null;
        if (print) {
            sb = new StringBuilder("\n Self-Consistent Field\n Iter  RMS Change (Debye)  Time\n");
        }

        // Find the induced dipole field due to direct dipoles (or predicted induced dipoles from previous steps).
        pme.computeInduceDipoleField();

        try {
            // Set initial conjugate gradient residual (a field).
            // Store the current induced dipoles and load the residual induced dipole.
            parallelTeam.execute(pcgInitRegion1);

            // Compute preconditioner.
            pme.expandInducedDipoles();
            // Use a special Ewald coefficient for the pre-conditioner.
            double aewaldTemp = ewaldParameters.aewald;
            ewaldParameters.setEwaldParameters(ewaldParameters.off, preconditionerEwald);
            int nAtoms = atoms.length;
            field.reset(parallelTeam, 0, nAtoms - 1);
            fieldCR.reset(parallelTeam, 0, nAtoms - 1);
            parallelTeam.execute(inducedDipolePreconditionerRegion);
            field.reduce(parallelTeam, 0, nAtoms - 1);
            fieldCR.reduce(parallelTeam, 0, nAtoms - 1);
            ewaldParameters.setEwaldParameters(ewaldParameters.off, aewaldTemp);
            // Revert to the stored induce dipoles.
            // Set initial conjugate vector (induced dipoles).
            parallelTeam.execute(pcgInitRegion2);
        } catch (Exception e) {
            String message = "Exception initializing preconditioned CG.";
            logger.log(Level.SEVERE, message, e);
        }

        // Conjugate gradient iteration of the mutual induced dipoles.
        int completedSCFCycles = 0;
        int maxSCFCycles = 1000;
        double eps = 100.0;
        double previousEps;
        boolean done = false;
        while (!done) {
            long cycleTime = -System.nanoTime();

            // Store a copy of the current induced dipoles, then set the induced dipoles to the conjugate vector.
            int nAtoms = atoms.length;
            for (int i = 0; i < nAtoms; i++) {
                vec[0][i] = inducedDipole[0][i][0];
                vec[1][i] = inducedDipole[0][i][1];
                vec[2][i] = inducedDipole[0][i][2];
                inducedDipole[0][i][0] = conj[0][i];
                inducedDipole[0][i][1] = conj[1][i];
                inducedDipole[0][i][2] = conj[2][i];
                vecCR[0][i] = inducedDipoleCR[0][i][0];
                vecCR[1][i] = inducedDipoleCR[0][i][1];
                vecCR[2][i] = inducedDipoleCR[0][i][2];
                inducedDipoleCR[0][i][0] = conjCR[0][i];
                inducedDipoleCR[0][i][1] = conjCR[1][i];
                inducedDipoleCR[0][i][2] = conjCR[2][i];
            }

            // Find the induced dipole field.
            pme.computeInduceDipoleField();

            try {
                /*
                 * Revert the induced dipoles to the saved values, then save the new residual field.
                 * Compute dot product of the conjugate vector and new residual.
                 * Reduce the residual field, add to the induced dipoles based
                 * on the scaled conjugate vector and finally set the induced
                 * dipoles to the polarizability times the residual field.
                 */
                parallelTeam.execute(pcgIterRegion1);

                // Compute preconditioner.
                pme.expandInducedDipoles();
                // Use a special Ewald coefficient for the pre-conditioner.
                double aewaldTemp = ewaldParameters.aewald;
                ewaldParameters.setEwaldParameters(ewaldParameters.off, preconditionerEwald);
                field.reset(parallelTeam, 0, nAtoms - 1);
                fieldCR.reset(parallelTeam, 0, nAtoms - 1);
                parallelTeam.execute(inducedDipolePreconditionerRegion);
                field.reduce(parallelTeam, 0, nAtoms - 1);
                fieldCR.reduce(parallelTeam, 0, nAtoms - 1);
                ewaldParameters.setEwaldParameters(ewaldParameters.off, aewaldTemp);

                /*
                 * Revert the induced dipoles to the saved values.
                 * Compute the dot product of the residual and preconditioner.
                 * Update the conjugate vector and sum the square of the residual field.
                 */
                pcgIterRegion2.sum = pcgIterRegion1.getSum();
                pcgIterRegion2.sumCR = pcgIterRegion1.getSumCR();
                parallelTeam.execute(pcgIterRegion2);
            } catch (Exception e) {
                String message = "Exception in first CG iteration region.";
                logger.log(Level.SEVERE, message, e);
            }

            previousEps = eps;
            eps = max(pcgIterRegion2.getEps(), pcgIterRegion2.getEpsCR());
            completedSCFCycles++;
            eps = MultipoleType.DEBYE * sqrt(eps / (double) nAtoms);
            cycleTime += System.nanoTime();
            if (print) {
                sb.append(format(
                        " %4d     %15.10f %7.4f\n", completedSCFCycles, eps, cycleTime * TO_SECONDS));
            }

            // If the RMS Debye change increases, fail the SCF process.
            if (eps > previousEps) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Fatal SCF convergence failure: (%10.5f > %10.5f)\n", eps, previousEps);
                throw new EnergyException(message, false);
            }

            // The SCF should converge well before the max iteration check. Otherwise, fail the SCF process.
            if (completedSCFCycles >= maxSCFCycles) {
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = format("Maximum SCF iterations reached: (%d)\n", completedSCFCycles);
                throw new EnergyException(message, false);
            }

            // Check if the convergence criteria has been achieved.
            if (eps < poleps) {
                done = true;
            }
        }
        if (print) {
            sb.append(format(" Direct:                  %7.4f\n", TO_SECONDS * directTime));
            startTime = System.nanoTime() - startTime;
            sb.append(format(" Total:                   %7.4f", startTime * TO_SECONDS));
            logger.info(sb.toString());
        }

        // Find the final induced dipole field.
        pme.computeInduceDipoleField();

        return completedSCFCycles;
    }

    private class PCGInitRegion1 extends ParallelRegion {

        private final PCGInitLoop[] pcgLoop;

        public PCGInitRegion1(int nt) {
            pcgLoop = new PCGInitLoop[nt];
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
                        rsd[0][i] = (directDipole[i][0] - inducedDipole[0][i][0]) * ipolar + field.getX(i);
                        rsd[1][i] = (directDipole[i][1] - inducedDipole[0][i][1]) * ipolar + field.getY(i);
                        rsd[2][i] = (directDipole[i][2] - inducedDipole[0][i][2]) * ipolar + field.getZ(i);
                        rsdCR[0][i] = (directDipoleCR[i][0] - inducedDipoleCR[0][i][0]) * ipolar + fieldCR.getX(i);
                        rsdCR[1][i] = (directDipoleCR[i][1] - inducedDipoleCR[0][i][1]) * ipolar + fieldCR.getY(i);
                        rsdCR[2][i] = (directDipoleCR[i][2] - inducedDipoleCR[0][i][2]) * ipolar + fieldCR.getZ(i);
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

    private class PCGInitRegion2 extends ParallelRegion {

        private final PCGInitLoop[] pcgLoop;

        public PCGInitRegion2(int nt) {
            pcgLoop = new PCGInitLoop[nt];
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
                    rsdPre[0][i] = polar * (field.getX(i) + udiag * rsd[0][i]);
                    rsdPre[1][i] = polar * (field.getY(i) + udiag * rsd[1][i]);
                    rsdPre[2][i] = polar * (field.getZ(i) + udiag * rsd[2][i]);
                    rsdPreCR[0][i] = polar * (fieldCR.getX(i) + udiag * rsdCR[0][i]);
                    rsdPreCR[1][i] = polar * (fieldCR.getY(i) + udiag * rsdCR[1][i]);
                    rsdPreCR[2][i] = polar * (fieldCR.getZ(i) + udiag * rsdCR[2][i]);
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

    private class PCGIterRegion1 extends ParallelRegion {

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
                        vec[0][i] = conj[0][i] * ipolar - field.getX(i);
                        vec[1][i] = conj[1][i] * ipolar - field.getY(i);
                        vec[2][i] = conj[2][i] * ipolar - field.getZ(i);
                        inducedDipoleCR[0][i][0] = vecCR[0][i];
                        inducedDipoleCR[0][i][1] = vecCR[1][i];
                        inducedDipoleCR[0][i][2] = vecCR[2][i];
                        vecCR[0][i] = conjCR[0][i] * ipolar - fieldCR.getX(i);
                        vecCR[1][i] = conjCR[1][i] * ipolar - fieldCR.getY(i);
                        vecCR[2][i] = conjCR[2][i] * ipolar - fieldCR.getZ(i);
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

    private class PCGIterRegion2 extends ParallelRegion {

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
                    rsdPre[0][i] = polar * (field.getX(i) + udiag * rsd[0][i]);
                    rsdPre[1][i] = polar * (field.getY(i) + udiag * rsd[1][i]);
                    rsdPre[2][i] = polar * (field.getZ(i) + udiag * rsd[2][i]);
                    rsdPreCR[0][i] = polar * (fieldCR.getX(i) + udiag * rsdCR[0][i]);
                    rsdPreCR[1][i] = polar * (fieldCR.getY(i) + udiag * rsdCR[1][i]);
                    rsdPreCR[2][i] = polar * (fieldCR.getZ(i) + udiag * rsdCR[2][i]);
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

    /**
     * Evaluate the real space field due to induced dipoles using a short cutoff
     * (~3-4 A).
     */
    private class InducedDipolePreconditionerRegion extends ParallelRegion {
        private final InducedPreconditionerFieldLoop[] inducedPreconditionerFieldLoop;

        InducedDipolePreconditionerRegion(int threadCount) {
            inducedPreconditionerFieldLoop = new InducedPreconditionerFieldLoop[threadCount];
        }

        @Override
        public void run() {
            int threadIndex = getThreadIndex();
            if (inducedPreconditionerFieldLoop[threadIndex] == null) {
                inducedPreconditionerFieldLoop[threadIndex] = new InducedPreconditionerFieldLoop();
            }
            try {
                int nAtoms = atoms.length;
                execute(0, nAtoms - 1, inducedPreconditionerFieldLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing the induced real space field in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class InducedPreconditionerFieldLoop extends IntegerForLoop {

            private int threadID;
            private double[] x, y, z;
            private double[][] ind, indCR;

            InducedPreconditionerFieldLoop() {
            }

            @Override
            public IntegerSchedule schedule() {
                return realSpaceSchedule;
            }

            @Override
            public void start() {
                threadID = getThreadIndex();
                realSpaceSCFTime[threadID] -= System.nanoTime();
                x = coordinates[0][0];
                y = coordinates[0][1];
                z = coordinates[0][2];
                ind = inducedDipole[0];
                indCR = inducedDipoleCR[0];
            }

            @Override
            public void finish() {
                realSpaceSCFTime[threadID] += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                final double[] dx = new double[3];
                final double[][] transOp = new double[3][3];

                // Loop over a chunk of atoms.
                int[][] lists = preconditionerLists[0];
                int[] counts = preconditionerCounts[0];
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    double fx = 0.0;
                    double fy = 0.0;
                    double fz = 0.0;
                    double px = 0.0;
                    double py = 0.0;
                    double pz = 0.0;
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double[] dipolei = ind[i];
                    final double uix = dipolei[0];
                    final double uiy = dipolei[1];
                    final double uiz = dipolei[2];
                    final double[] dipoleCRi = indCR[i];
                    final double pix = dipoleCRi[0];
                    final double piy = dipoleCRi[1];
                    final double piz = dipoleCRi[2];
                    final double pdi = ipdamp[i];
                    final double pti = thole[i];

                    // Loop over the neighbor list.
                    final int[] list = lists[i];
                    final int npair = counts[i];
                    for (int j = 0; j < npair; j++) {
                        final int k = list[j];
                        if (!use[k]) {
                            continue;
                        }
                        final double pdk = ipdamp[k];
                        final double ptk = thole[k];
                        dx[0] = x[k] - xi;
                        dx[1] = y[k] - yi;
                        dx[2] = z[k] - zi;
                        final double r2 = crystal.image(dx);

                        // Calculate the error function damping terms.
                        final double r = sqrt(r2);
                        final double rr1 = 1.0 / r;
                        final double rr2 = rr1 * rr1;
                        final double ralpha = ewaldParameters.aewald * r;
                        final double exp2a = exp(-ralpha * ralpha);
                        final double bn0 = erfc(ralpha) * rr1;
                        // final double exp2a = 1.0;
                        // final double bn0 = rr1;
                        final double bn1 = (bn0 + ewaldParameters.an0 * exp2a) * rr2;
                        final double bn2 = (3.0 * bn1 + ewaldParameters.an1 * exp2a) * rr2;
                        double scale3 = 1.0;
                        double scale5 = 1.0;
                        double damp = pdi * pdk;
                        final double pgamma = min(pti, ptk);
                        final double rdamp = r * damp;
                        damp = -pgamma * rdamp * rdamp * rdamp;
                        if (damp > -50.0) {
                            final double expdamp = exp(damp);
                            scale3 = 1.0 - expdamp;
                            scale5 = 1.0 - expdamp * (1.0 - damp);
                        }
                        double rr3 = rr1 * rr2;
                        double rr5 = 3.0 * rr3 * rr2;
                        rr3 *= (1.0 - scale3);
                        rr5 *= (1.0 - scale5);
                        final double xr = dx[0];
                        final double yr = dx[1];
                        final double zr = dx[2];
                        final double[] dipolek = ind[k];
                        final double ukx = dipolek[0];
                        final double uky = dipolek[1];
                        final double ukz = dipolek[2];
                        final double ukr = ukx * xr + uky * yr + ukz * zr;
                        final double bn2ukr = bn2 * ukr;
                        final double fimx = -bn1 * ukx + bn2ukr * xr;
                        final double fimy = -bn1 * uky + bn2ukr * yr;
                        final double fimz = -bn1 * ukz + bn2ukr * zr;
                        final double rr5ukr = rr5 * ukr;
                        final double fidx = -rr3 * ukx + rr5ukr * xr;
                        final double fidy = -rr3 * uky + rr5ukr * yr;
                        final double fidz = -rr3 * ukz + rr5ukr * zr;
                        fx += (fimx - fidx);
                        fy += (fimy - fidy);
                        fz += (fimz - fidz);
                        final double[] dipolepk = indCR[k];
                        final double pkx = dipolepk[0];
                        final double pky = dipolepk[1];
                        final double pkz = dipolepk[2];
                        final double pkr = pkx * xr + pky * yr + pkz * zr;
                        final double bn2pkr = bn2 * pkr;
                        final double pimx = -bn1 * pkx + bn2pkr * xr;
                        final double pimy = -bn1 * pky + bn2pkr * yr;
                        final double pimz = -bn1 * pkz + bn2pkr * zr;
                        final double rr5pkr = rr5 * pkr;
                        final double pidx = -rr3 * pkx + rr5pkr * xr;
                        final double pidy = -rr3 * pky + rr5pkr * yr;
                        final double pidz = -rr3 * pkz + rr5pkr * zr;
                        px += (pimx - pidx);
                        py += (pimy - pidy);
                        pz += (pimz - pidz);
                        final double uir = uix * xr + uiy * yr + uiz * zr;
                        final double bn2uir = bn2 * uir;
                        final double fkmx = -bn1 * uix + bn2uir * xr;
                        final double fkmy = -bn1 * uiy + bn2uir * yr;
                        final double fkmz = -bn1 * uiz + bn2uir * zr;
                        final double rr5uir = rr5 * uir;
                        final double fkdx = -rr3 * uix + rr5uir * xr;
                        final double fkdy = -rr3 * uiy + rr5uir * yr;
                        final double fkdz = -rr3 * uiz + rr5uir * zr;
                        field.add(threadID, k, fkmx - fkdx, fkmy - fkdy, fkmz - fkdz);
                        final double pir = pix * xr + piy * yr + piz * zr;
                        final double bn2pir = bn2 * pir;
                        final double pkmx = -bn1 * pix + bn2pir * xr;
                        final double pkmy = -bn1 * piy + bn2pir * yr;
                        final double pkmz = -bn1 * piz + bn2pir * zr;
                        final double rr5pir = rr5 * pir;
                        final double pkdx = -rr3 * pix + rr5pir * xr;
                        final double pkdy = -rr3 * piy + rr5pir * yr;
                        final double pkdz = -rr3 * piz + rr5pir * zr;
                        fieldCR.add(threadID, k, pkmx - pkdx, pkmy - pkdy, pkmz - pkdz);
                    }
                    field.add(threadID, i, fx, fy, fz);
                    fieldCR.add(threadID, i, px, py, pz);
                }

                // Loop over symmetry mates.
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                int nSymm = symOps.size();
                for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                    SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                    crystal.getTransformationOperator(symOp, transOp);
                    lists = preconditionerLists[iSymm];
                    counts = preconditionerCounts[iSymm];
                    final double[] xs = coordinates[iSymm][0];
                    final double[] ys = coordinates[iSymm][1];
                    final double[] zs = coordinates[iSymm][2];
                    final double[][] inds = inducedDipole[iSymm];
                    final double[][] indCRs = inducedDipoleCR[iSymm];

                    // Loop over a chunk of atoms.
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        double fx = 0.0;
                        double fy = 0.0;
                        double fz = 0.0;
                        double px = 0.0;
                        double py = 0.0;
                        double pz = 0.0;
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final double[] dipolei = ind[i];
                        final double uix = dipolei[0];
                        final double uiy = dipolei[1];
                        final double uiz = dipolei[2];
                        final double[] dipoleCRi = indCR[i];
                        final double pix = dipoleCRi[0];
                        final double piy = dipoleCRi[1];
                        final double piz = dipoleCRi[2];
                        final double pdi = ipdamp[i];
                        final double pti = thole[i];

                        // Loop over the neighbor list.
                        final int[] list = lists[i];
                        final int npair = counts[i];
                        for (int j = 0; j < npair; j++) {
                            final int k = list[j];
                            if (!use[k]) {
                                continue;
                            }
                            double selfScale = 1.0;
                            if (i == k) {
                                selfScale = 0.5;
                            }
                            final double pdk = ipdamp[k];
                            final double ptk = thole[k];
                            dx[0] = xs[k] - xi;
                            dx[1] = ys[k] - yi;
                            dx[2] = zs[k] - zi;
                            final double r2 = crystal.image(dx);

                            // Calculate the error function damping terms.
                            final double r = sqrt(r2);
                            final double rr1 = 1.0 / r;
                            final double rr2 = rr1 * rr1;
                            final double ralpha = ewaldParameters.aewald * r;
                            final double exp2a = exp(-ralpha * ralpha);
                            final double bn0 = erfc(ralpha) * rr1;
                            //final double exp2a = 1.0;
                            //final double bn0 = rr1;
                            final double bn1 = (bn0 + ewaldParameters.an0 * exp2a) * rr2;
                            final double bn2 = (3.0 * bn1 + ewaldParameters.an1 * exp2a) * rr2;
                            double scale3 = 1.0;
                            double scale5 = 1.0;
                            double damp = pdi * pdk;
                            final double pgamma = min(pti, ptk);
                            final double rdamp = r * damp;
                            damp = -pgamma * rdamp * rdamp * rdamp;
                            if (damp > -50.0) {
                                final double expdamp = exp(damp);
                                scale3 = 1.0 - expdamp;
                                scale5 = 1.0 - expdamp * (1.0 - damp);
                            }
                            double rr3 = rr1 * rr2;
                            double rr5 = 3.0 * rr3 * rr2;
                            rr3 *= (1.0 - scale3);
                            rr5 *= (1.0 - scale5);
                            final double xr = dx[0];
                            final double yr = dx[1];
                            final double zr = dx[2];
                            final double[] dipolek = inds[k];
                            final double ukx = dipolek[0];
                            final double uky = dipolek[1];
                            final double ukz = dipolek[2];
                            final double[] dipolepk = indCRs[k];
                            final double pkx = dipolepk[0];
                            final double pky = dipolepk[1];
                            final double pkz = dipolepk[2];
                            final double ukr = ukx * xr + uky * yr + ukz * zr;
                            final double bn2ukr = bn2 * ukr;
                            final double fimx = -bn1 * ukx + bn2ukr * xr;
                            final double fimy = -bn1 * uky + bn2ukr * yr;
                            final double fimz = -bn1 * ukz + bn2ukr * zr;
                            final double rr5ukr = rr5 * ukr;
                            final double fidx = -rr3 * ukx + rr5ukr * xr;
                            final double fidy = -rr3 * uky + rr5ukr * yr;
                            final double fidz = -rr3 * ukz + rr5ukr * zr;
                            fx += selfScale * (fimx - fidx);
                            fy += selfScale * (fimy - fidy);
                            fz += selfScale * (fimz - fidz);
                            final double pkr = pkx * xr + pky * yr + pkz * zr;
                            final double bn2pkr = bn2 * pkr;
                            final double pimx = -bn1 * pkx + bn2pkr * xr;
                            final double pimy = -bn1 * pky + bn2pkr * yr;
                            final double pimz = -bn1 * pkz + bn2pkr * zr;
                            final double rr5pkr = rr5 * pkr;
                            final double pidx = -rr3 * pkx + rr5pkr * xr;
                            final double pidy = -rr3 * pky + rr5pkr * yr;
                            final double pidz = -rr3 * pkz + rr5pkr * zr;
                            px += selfScale * (pimx - pidx);
                            py += selfScale * (pimy - pidy);
                            pz += selfScale * (pimz - pidz);
                            final double uir = uix * xr + uiy * yr + uiz * zr;
                            final double bn2uir = bn2 * uir;
                            final double fkmx = -bn1 * uix + bn2uir * xr;
                            final double fkmy = -bn1 * uiy + bn2uir * yr;
                            final double fkmz = -bn1 * uiz + bn2uir * zr;
                            final double rr5uir = rr5 * uir;
                            final double fkdx = -rr3 * uix + rr5uir * xr;
                            final double fkdy = -rr3 * uiy + rr5uir * yr;
                            final double fkdz = -rr3 * uiz + rr5uir * zr;
                            double xc = selfScale * (fkmx - fkdx);
                            double yc = selfScale * (fkmy - fkdy);
                            double zc = selfScale * (fkmz - fkdz);
                            double fkx = (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                            double fky = (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                            double fkz = (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                            field.add(threadID, k, fkx, fky, fkz);
                            final double pir = pix * xr + piy * yr + piz * zr;
                            final double bn2pir = bn2 * pir;
                            final double pkmx = -bn1 * pix + bn2pir * xr;
                            final double pkmy = -bn1 * piy + bn2pir * yr;
                            final double pkmz = -bn1 * piz + bn2pir * zr;
                            final double rr5pir = rr5 * pir;
                            final double pkdx = -rr3 * pix + rr5pir * xr;
                            final double pkdy = -rr3 * piy + rr5pir * yr;
                            final double pkdz = -rr3 * piz + rr5pir * zr;
                            xc = selfScale * (pkmx - pkdx);
                            yc = selfScale * (pkmy - pkdy);
                            zc = selfScale * (pkmz - pkdz);
                            fkx = (xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0]);
                            fky = (xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1]);
                            fkz = (xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2]);
                            fieldCR.add(threadID, k, fkx, fky, fkz);
                        }
                        field.add(threadID, i, fx, fy, fz);
                        fieldCR.add(threadID, i, px, py, pz);
                    }
                }
            }
        }
    }

    private static final double TO_SECONDS = 1.0e-9;
}