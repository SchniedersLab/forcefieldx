/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
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
package ffx.potential.utils;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.Random;
import java.util.logging.Logger;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import static org.junit.Assert.assertEquals;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;

/**
 * Test the ForceFieldEnergy class for small crystals.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class OPLSCrystalEnergyTest {

    private static final Logger logger = Logger.getLogger(OPLSCrystalEnergyTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                false,
                "OPLS Acetanilide Benchmark",
                "ffx/potential/structures/acetanilide-oplsaa.xyz",
                0.45009773, 19,
                1.99659439, 30,
                3.12713717, 38,
                0.03415131, 12,
                -10.22874283, 7293,
                -24.11087720, 2229},
            {
                false,
                "OPLS Ethylparaben Benchmark",
                "ffx/potential/structures/ethylparaben-oplsaa.xyz",
                1.23483793, 44,
                1.17950432, 70,
                0.58946840, 88,
                0.01347012, 22,
                -23.22505628, 16759,
                -28.45972542, 5029},
            {
                false,
                "OPLS Methylparaben Benchmark",
                "ffx/potential/structures/methylparaben-oplsaa.xyz",
                0.48582171, 19,
                0.47331501, 29,
                0.53855046, 35,
                0.02598434, 11,
                -8.71533364, 7647,
                -16.384524645738537, 2302},
            {
                false,
                "OPLS Paracetamol Benchmark",
                "ffx/potential/structures/paracetamol-oplsaa.xyz",
                0.63722563, 20,
                2.63545869, 31,
                3.01547224, 40,
                0.06569712, 10,
                -9.45895598, 7831,
                -41.34298872, 3334},
            {
                false,
                "OPLS Phenacetin Benchmark",
                "ffx/potential/structures/phenacetin-oplsaa.xyz",
                0.53495810, 26,
                3.68523677, 43,
                4.12682233, 52,
                0.03932507, 10,
                -14.31979877, 10338,
                -31.561832525438113, 3116}
        });
    }
    private final String info;
    private final File structure;
    private final MolecularAssembly molecularAssembly;
    private final int nBonds;
    private final int nAngles;
    private final int nTorsions;
    private final int nImproperTorsions;
    private final int nVanDerWaals;
    private final int nPermanent;
    private final double bondEnergy;
    private final double angleEnergy;
    private final double torsionEnergy;
    private final double improperTorsionEnergy;
    private final double vanDerWaalsEnergy;
    private final double permanentEnergy;
    private final double tolerance = 1.0e-2;
    private final double gradientTolerance = 1.0e-2;
    private final boolean ci;
    private final boolean ciOnly;
    private boolean mpoleTerm;
    private final ForceFieldEnergy forceFieldEnergy;

    public OPLSCrystalEnergyTest(
            boolean ciOnly,
            String info, String filename,
            double bondEnergy, int nBonds,
            double angleEnergy, int nAngles,
            double torsionEnergy, int nTorsions,
            double improperTorsionEnergy, int nImproperTorsions,
            double vanDerWaalsEnergy, int nVanDerWaals,
            double permanentEnergy, int nPermanent) {
        this.ciOnly = ciOnly;
        this.info = info;
        this.bondEnergy = bondEnergy;
        this.nBonds = nBonds;
        this.angleEnergy = angleEnergy;
        this.nAngles = nAngles;
        this.torsionEnergy = torsionEnergy;
        this.nTorsions = nTorsions;
        this.improperTorsionEnergy = improperTorsionEnergy;
        this.nImproperTorsions = nImproperTorsions;
        this.vanDerWaalsEnergy = vanDerWaalsEnergy;
        this.nVanDerWaals = nVanDerWaals;
        this.permanentEnergy = permanentEnergy;
        this.nPermanent = nPermanent;

        ci = System.getProperty("ffx.ci","false").equalsIgnoreCase("true");
        if (!ci && ciOnly) {
            structure = null;
            molecularAssembly = null;
            forceFieldEnergy = null;
            return;
        }

        /**
         * Load the test system.
         */
        ClassLoader cl = this.getClass().getClassLoader();
        structure = new File(cl.getResource(filename).getPath());
        PotentialsUtils potentialUtils = new PotentialsUtils();
        molecularAssembly = potentialUtils.open(structure.getAbsolutePath())[0];

        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        mpoleTerm = molecularAssembly.getForceField().getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, true);

        if (ci) {
            testGradient();
            testSoftCore();
        }
    }

    /**
     * Test of energy method, of class PotentialEnergy.
     */
    @Test
    public void testEnergy() {

        /**
         * Skip expensive tests for normal builds.
         */
        if (!ci && ciOnly) {
            return;
        }

        boolean gradient = false;
        boolean print = true;
        double total = forceFieldEnergy.energy(gradient, print);
        // Bond Energy
        assertEquals(info + " Bond Energy", bondEnergy, forceFieldEnergy.getBondEnergy(), tolerance);
        assertEquals(info + " Bond Count", nBonds, forceFieldEnergy.getNumberofBonds());
        // Angle Energy
        assertEquals(info + " Angle Energy", angleEnergy, forceFieldEnergy.getAngleEnergy(), tolerance);
        assertEquals(info + " Angle Count", nAngles, forceFieldEnergy.getNumberofAngles());
        // Torsional Angle
        assertEquals(info + " Torsion Energy", torsionEnergy, forceFieldEnergy.getTorsionEnergy(), tolerance);
        assertEquals(info + " Torsion Count", nTorsions, forceFieldEnergy.getNumberofTorsions());
        // Improper Torsional Angle
        assertEquals(info + " Improper Torsion Energy", improperTorsionEnergy, forceFieldEnergy.getImproperTorsionEnergy(), tolerance);
        assertEquals(info + " Improper Torsion Count", nImproperTorsions, forceFieldEnergy.getNumberofImproperTorsions());
        // van Der Waals
        assertEquals(info + " van Der Waals Energy", vanDerWaalsEnergy, forceFieldEnergy.getVanDerWaalsEnergy(), tolerance);
        assertEquals(info + " van Der Waals Count", nVanDerWaals, forceFieldEnergy.getVanDerWaalsInteractions());
        // Permanent Multipoles
        if (mpoleTerm) {
            assertEquals(info + " Permanent Multipole Energy", permanentEnergy,
                    forceFieldEnergy.getPermanentMultipoleEnergy(), tolerance);
            assertEquals(info + " Permanent Multipole Count", nPermanent, forceFieldEnergy.getPermanentInteractions());
        }
    }

    /**
     * Test of energy gradient, of class PotentialEnergy.
     */
    public void testGradient() {
        boolean gradient = true;
        boolean print = true;
        forceFieldEnergy.energy(gradient, print);
        gradient = false;
        print = false;
        double step = 0.0001;
        double analytic[] = new double[3];
        double numeric[] = new double[3];
        double xyz[] = new double[3];

        Atom[] atoms = molecularAssembly.getAtomArray();
        int n = atoms.length;
        Random random = new Random();
        int i = random.nextInt(n);
        Atom a0 = atoms[i];

        a0.getXYZGradient(analytic);
        a0.getXYZ(xyz);
        // Find numeric dX
        xyz[0] += step;
        a0.moveTo(xyz);
        double e = forceFieldEnergy.energy(gradient, print);
        xyz[0] -= 2.0 * step;
        a0.moveTo(xyz);
        e -= forceFieldEnergy.energy(gradient, print);
        numeric[0] = e / (2.0 * step);
        xyz[0] += step;
        // Find numeric dY
        xyz[1] += step;
        a0.moveTo(xyz);
        e = forceFieldEnergy.energy(gradient, print);
        xyz[1] -= 2.0 * step;
        a0.moveTo(xyz);
        e -= forceFieldEnergy.energy(gradient, print);
        numeric[1] = e / (2.0 * step);
        xyz[1] += step;
        // Find numeric dZ
        xyz[2] += step;
        a0.moveTo(xyz);
        e = forceFieldEnergy.energy(gradient, print);
        xyz[2] -= 2.0 * step;
        a0.moveTo(xyz);
        e -= forceFieldEnergy.energy(gradient, print);
        numeric[2] = e / (2.0 * step);
        xyz[2] += step;
        a0.moveTo(xyz);
        double dx = analytic[0] - numeric[0];
        double dy = analytic[1] - numeric[1];
        double dz = analytic[2] - numeric[2];
        double len = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (len > gradientTolerance) {
            logger.severe("\n" + a0.toString() + String.format(" failed: %10.6f.", len) + String.format(
                    "\nAnalytic: (%12.4f, %12.4f, %12.4f)\n",
                    analytic[0], analytic[1], analytic[2]) + String.format("Numeric:  (%12.4f, %12.4f, %12.4f)\n",
                            numeric[0], numeric[1], numeric[2]));
        } else {
            logger.info("\n" + a0.toString() + String.format(" passed: %10.6f.", len) + String.format(
                    "\nAnalytic: (%12.4f, %12.4f, %12.4f)\n",
                    analytic[0], analytic[1], analytic[2]) + String.format("Numeric:  (%12.4f, %12.4f, %12.4f)\n",
                            numeric[0], numeric[1], numeric[2]));
        }
        assertEquals(a0.toString(), 0.0, len, gradientTolerance);
    }

    public void testSoftCore() {
        boolean gradient = false;
        boolean print = true;
        double e = forceFieldEnergy.energy(gradient, print);
        Atom atoms[] = molecularAssembly.getAtomArray();
        int n = atoms.length;
        // Make the 3 atoms soft
        for (int i = n; i > n - 3; i--) {
            Atom ai = atoms[i - 1];
            ai.setApplyLambda(true);
        }
        // Compute the energy with Lambda = 1.0;
        double lambda = 1.0;
        forceFieldEnergy.setLambda(lambda);
        double e2 = forceFieldEnergy.energy(gradient, print);
        assertEquals(e, e2, tolerance);
    }
}
