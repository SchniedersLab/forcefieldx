/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
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
package ffx.potential.nonbonded;

import static java.lang.Math.min;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldInteger;
import ffx.potential.parameters.ForceField.ForceFieldString;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The PMEWisdom class searches through Ewald coefficients and
 * {@link ReciprocalSpace} grid spacings to find the most efficient combination
 * that maintains an RMS gradient error below a specified tolarance
 * (typically 10-4 RMS Kcal/mole/Angstrom or better).
 * This is especially useful for {@link Crystal} space groups with a large number
 * of symmetry operators. The gold standard gradients are computed using:
 * <ol>
 * <li>An Ewald coefficient of 0.42.</li>
 * <li>A real space cutoff of 10.0 Angstroms.</li>
 * <li>A grid spacing of 0.5 Angstroms.</li>
 * <li>A b-Spline order of 10.</li>
 * </ol>
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PMEWisdom {

    private static final Logger logger = Logger.getLogger(PMEWisdom.class.getName());
    
    private final MolecularAssembly molecularAssembly;
    private final ForceField forceField;
    private final Crystal crystal;
    private final int nAtoms;
    private final double highAccuracyGradients[][];
    private final double hx[], hy[], hz[];
    private final double gradients[][];
    private final double gx[], gy[], gz[];
    private final Atom atoms[];
    private final int neighborLists[][][];
    private final double coordinates[][];
    private final double buffer;
    private final ParallelTeam parallelTeam;

    /**
     * The PMEWisdom constructor requires a MolecularAssembly that is a
     * periodic.
     *
     * @param molecularAssembly
     * @since 1.0
     */
    public PMEWisdom(MolecularAssembly molecularAssembly) {
        this.molecularAssembly = molecularAssembly;
        this.forceField = molecularAssembly.getForceField();
        this.parallelTeam = new ParallelTeam();
        nAtoms = molecularAssembly.getAtomList().size();
        highAccuracyGradients = new double[3][nAtoms];
        hx = highAccuracyGradients[0];
        hy = highAccuracyGradients[1];
        hz = highAccuracyGradients[2];
        gradients = new double[3][nAtoms];
        gx = gradients[0];
        gy = gradients[1];
        gz = gradients[2];
        final double a = forceField.getDouble(ForceFieldDouble.A_AXIS, 10.0);
        final double b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
        final double c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
        final double alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
        final double beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
        final double gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
        final String spacegroup = forceField.getString(
                ForceFieldString.SPACEGROUP, "P1");
        crystal = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
        logger.info(crystal.toString());

        ArrayList<Atom> atomList = molecularAssembly.getAtomList();
        atoms = atomList.toArray(new Atom[nAtoms]);
        Arrays.sort(atoms);

        // Create the neighbor list.
        int nSymm = crystal.spaceGroup.symOps.size();
        neighborLists = new int[nSymm][][];
        coordinates = new double[nSymm][nAtoms * 3];
        buffer = 2.0;

        /**
         * Update the local coordinate array.
         */
        for (int i = 0; i < nAtoms; i++) {
            final double xyz[] = atoms[i].getXYZ();
            int i3 = i * 3;
            coordinates[0][i3++] = xyz[0];
            coordinates[0][i3++] = xyz[1];
            coordinates[0][i3] = xyz[2];
        }
        Vector<SymOp> symOps = crystal.spaceGroup.symOps;
        double in[] = new double[3];
        double out[] = new double[3];
        for (int iSymOp = 1; iSymOp < nSymm; iSymOp++) {
            SymOp symOp = symOps.get(iSymOp);
            double xyz[] = coordinates[iSymOp];
            for (int i = 0; i < nAtoms; i++) {
                int i3 = i * 3;
                int iX = i3;
                int iY = i3 + 1;
                int iZ = i3 + 2;
                in[0] = coordinates[0][iX];
                in[1] = coordinates[0][iY];
                in[2] = coordinates[0][iZ];
                crystal.applySymOp(in, out, symOp);
                xyz[iX] = out[0];
                xyz[iY] = out[1];
                xyz[iZ] = out[2];
            }
        }
    }
    private static double toSeconds = 0.000000001;

    public void run() {
        final double maxCutoff = min(min(crystal.a, crystal.b), crystal.c) /
                2.0 - buffer;

        /**
         * Following Essmann et al. (1995), we choose a real space Ewald cutoff
         * of 9.0 Angstroms and require that erfc(Beta * r) / r < 10^(-8) at
         * the cutoff. A beta value (aka the Ewald coefficient) of 0.42
         * satisfies this criteria.
         * <p>
         * In the limit of infinite b-Spline order a Gaussian is
         * achieved. Here we use order 10. Finally, a reciprocal space
         * grid spacing of 0.66 A is chosen.
         */
        double beta = 0.42;
        double cutoff = ParticleMeshEwald.ewaldCutoff(beta, maxCutoff, 1e-10);

        forceField.addForceFieldDouble(ForceFieldDouble.EWALD_CUTOFF, cutoff);
        forceField.addForceFieldDouble(ForceFieldDouble.EWALD_ALPHA, beta);
        forceField.addForceFieldDouble(ForceFieldDouble.PME_MESH_DENSITY, 0.5);
        forceField.addForceFieldInteger(ForceFieldInteger.PME_ORDER, 10);

        NeighborList neighborList = new NeighborList(null, crystal, atoms, cutoff, buffer,
                parallelTeam);
        neighborList.buildList(coordinates, neighborLists, null, true, true);
        ParticleMeshEwald particleMeshEwald = new ParticleMeshEwald(forceField,
                atoms, crystal, parallelTeam, neighborLists, neighborList.getPairwiseSchedule());

        /**
         * Time the high accuracy energy gradients.
         */
        long bestTime = System.nanoTime();
        double highAccuracyEnergy = particleMeshEwald.energy(true, false);
        bestTime = System.nanoTime() - bestTime;
        logger.info(String.format("High Accuracy Time:           %5.3f\n", bestTime *
                toSeconds));
        particleMeshEwald.getGradients(highAccuracyGradients);

        double realSpaceTolerance = 1.0e-6;
        double rmsGradientTolerance = 1.0e-5;
        beta = 0.42;
        cutoff = ParticleMeshEwald.ewaldCutoff(beta, maxCutoff, realSpaceTolerance);

        while (true) {
            long newTime = findPMEGridSpacing(cutoff, beta, 6, 0.6,
                    rmsGradientTolerance);
            if (newTime < bestTime) {
                bestTime = newTime;
                beta -= 0.02;
                cutoff = ParticleMeshEwald.ewaldCutoff(beta, maxCutoff,
                        realSpaceTolerance);
                if (cutoff > maxCutoff) {
                    logger.info("Breaking due to large real space cutoff.");
                    break;
                }
            } else {
                logger.info("Breaking due to slow time.\n" +
                        String.format("Best Time:                   %5.3f\n", bestTime *
                        toSeconds) +
                        String.format("New Time:                    %5.3f\n", newTime *
                        toSeconds));
                break;
            }
        }
        beta += 0.02;
        cutoff = ParticleMeshEwald.ewaldCutoff(beta, maxCutoff,
                realSpaceTolerance);

    }

    private long findPMEGridSpacing(double cutoff, double alpha, int order,
            double initial, double gradientTolerance) {
        double spacing = initial - 0.1;
        double previousRMS = 0.0;
        double rms = 0.0;
        NeighborList neighborList = new NeighborList(null, crystal, atoms, cutoff, buffer,
                    parallelTeam);
        neighborList.buildList(coordinates, neighborLists, null, true, false);

        logger.setLevel(Level.INFO);
        logger.info(String.format("RMS Gradient Tolerance: %5.3f\n",
                gradientTolerance));

        while (rms < gradientTolerance) {
            previousRMS = rms;
            // Increase the grid spacing until the rms force error is too large.
            spacing += 0.1;
            logger.info(
                    String.format("Evaluating Grid Spacing: %5.3f\n", spacing));
            forceField.addForceFieldDouble(ForceFieldDouble.EWALD_CUTOFF, cutoff);
            forceField.addForceFieldDouble(ForceFieldDouble.EWALD_ALPHA, alpha);
            forceField.addForceFieldDouble(ForceFieldDouble.PME_MESH_DENSITY, spacing);
            forceField.addForceFieldInteger(ForceFieldInteger.PME_ORDER, order);
            ParticleMeshEwald particleMeshEwald = new ParticleMeshEwald(forceField, atoms, crystal,
                    parallelTeam, neighborLists, neighborList.getPairwiseSchedule());
            System.gc();
            particleMeshEwald.energy(true, false);
            particleMeshEwald.getGradients(gradients);
            /* Calculate the RMS gradient error. */
            double denom = 0.0;
            for (int i = 0; i < nAtoms; i++) {
                double hxi = hx[i];
                double hyi = hy[i];
                double hzi = hz[i];
                double dx = hxi - gx[i];
                double dy = hyi - gy[i];
                double dz = hzi - gz[i];
                rms += dx * dx;
                rms += dy * dy;
                rms += dz * dz;
                denom += hxi * hxi + hyi * hyi + hzi * hzi;
            }
            rms = Math.sqrt(rms / denom);
            logger.info(String.format(
                    "RMS Gradient Error:     %6.3e (%6.3e percent of maximum)\n",
                    rms, rms / gradientTolerance * 100));
        }
        spacing -= 0.1;

        // Find the best timing for these PME parameters.
        forceField.addForceFieldDouble(ForceFieldDouble.EWALD_CUTOFF, cutoff);
        forceField.addForceFieldDouble(ForceFieldDouble.EWALD_ALPHA, alpha);
        forceField.addForceFieldDouble(ForceFieldDouble.PME_MESH_DENSITY, spacing);
        forceField.addForceFieldInteger(ForceField.ForceFieldInteger.PME_ORDER,
                order);
        ParticleMeshEwald particleMeshEwald = new ParticleMeshEwald(forceField,
                atoms, crystal, parallelTeam, neighborLists, neighborList.getPairwiseSchedule());
        System.gc();
        long bestTime = System.nanoTime();
        particleMeshEwald.energy(true, false);
        bestTime = System.nanoTime() - bestTime;
        long newTime = bestTime;
        for (int i = 0; i < 5; i++) {
            newTime = System.nanoTime();
            if (i == 4) {
                logger.setLevel(Level.FINE);
            }
            particleMeshEwald.energy(true, false);
            if (i == 4) {
                logger.setLevel(Level.INFO);
            }
            newTime = System.nanoTime() - newTime;
            if (newTime < bestTime) {
                bestTime = newTime;
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append(String.format("\nRMS Gradient Tolerance: %6.3e\n",
                gradientTolerance));
        sb.append(String.format("RMS Gradient Error:     %6.3e\n", previousRMS));
        sb.append(String.format("Ewald Coefficient:      %5.3f\n", alpha));
        sb.append(String.format("Real Space Cutoff:      %5.3f\n", cutoff));
        sb.append(String.format("Grid Spacing:           %5.3f\n", spacing));
        sb.append(String.format("Time:                   %5.3f\n", bestTime *
                toSeconds));
        logger.info(sb.toString());
        return bestTime;
    }
}
