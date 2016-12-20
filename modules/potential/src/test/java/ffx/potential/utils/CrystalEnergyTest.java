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
public class CrystalEnergyTest {

    private static final Logger logger = Logger.getLogger(CrystalEnergyTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                false,
                "Acetanilide Benchmark",
                "ffx/potential/structures/acetanilide.xyz",
                0.67814832, 19,
                1.06817134, 30,
                -0.06290501, 26,
                0.0, 0,
                0.13033970, 24,
                -8.01262794, 38,
                0.0, 0,
                0.0, 0,
                0.0, 0,
                -3.62572228, 7295,
                -26.0720131242538, 2231,
                -1.9856761645305, 2231},
            {
                false,
                "Ethylparaben Benchmark",
                "ffx/potential/structures/ethylparaben.xyz",
                1.37991170, 44,
                2.43600964, 70,
                0.20006782, 58,
                0.0, 0,
                0.04434686, 42,
                -7.14298308, 88,
                0.0, 0,
                0.0, 0,
                0.0, 0,
                -4.52561611, 16755,
                -45.55012417320137, 5012,
                -3.97079858673400, 5012},
            {
                false,
                "Methylparaben Benchmark",
                "ffx/potential/structures/methylparaben.xyz",
                0.42501694, 19,
                1.23611203, 29,
                0.00351228, 24,
                0.0, 0,
                0.01058196, 21,
                -3.72750990, 35,
                0.0, 0,
                0.0, 0,
                0.0, 0,
                0.20573979, 7639,
                -26.558577471708134, 2314,
                -2.35904592, 2314},
            {
                false,
                "Paracetamol Benchmark",
                "ffx/potential/structures/paracetamol.xyz",
                0.55316710, 20,
                1.51848706, 31,
                0.05096337, 27,
                0.0, 0,
                0.46651316, 24,
                -7.40125527, 40,
                0.0, 0,
                0.0, 0,
                0.0, 0,
                0.76621272, 7832,
                -32.02011297249507, 2357,
                 -4.96649992387851, 2357},
            {
                false,
                "Phenacetin Benchmark",
                "ffx/potential/structures/phenacetin.xyz",
                0.56492623, 26,
                2.68273715, 43,
                0.00409248, 35,
                0.0, 0,
                0.33385711, 24,
                -6.49787814, 52,
                0.0, 0,
                0.0, 0,
                0.0, 0,
                -5.62144406, 10340,
                -25.42173449818894, 3138,
                 -2.09385416786188, 3138}
        });
    }
    private final String info;
    private final File structure;
    private final MolecularAssembly molecularAssembly;
    private final int nBonds;
    private final int nAngles;
    private final int nStretchBends;
    private final int nUreyBradleys;
    private final int nOutOfPlaneBends;
    private final int nTorsions;
    private final int nImproperTorsions;
    private final int nPiOrbitalTorsions;
    private final int nTorsionTorsions;
    private final int nVanDerWaals;
    private final int nPermanent;
    private final int nPolar;
    private final double bondEnergy;
    private final double angleEnergy;
    private final double stretchBendEnergy;
    private final double ureyBradleyEnergy;
    private final double outOfPlaneBendEnergy;
    private final double torsionEnergy;
    private final double improperTorsionEnergy;
    private final double piOrbitalTorsionEnergy;
    private final double torsionTorsionEnergy;
    private final double vanDerWaalsEnergy;
    private final double permanentEnergy;
    private final double polarizationEnergy;
    private final double tolerance = 1.0e-3;
    private final double gradientTolerance = 1.0e-3;
    private final boolean ci;
    private final boolean ciOnly;
    private boolean mpoleTerm;
    private final ForceFieldEnergy forceFieldEnergy;

    public CrystalEnergyTest(
            boolean ciOnly,
            String info, String filename,
            double bondEnergy, int nBonds,
            double angleEnergy, int nAngles,
            double stretchBendEnergy, int nStretchBends,
            double ureyBradleyEnergy, int nUreyBradleys,
            double outOfPlaneBendEnergy, int nOutOfPlaneBends,
            double torsionEnergy, int nTorsions,
            double improperTorsionEnergy, int nImproperTorsions,
            double piOrbitalTorsionEnergy, int nPiOrbitalTorsions,
            double torsionTorsionEnergy, int nTorsionTorsions,
            double vanDerWaalsEnergy, int nVanDerWaals,
            double permanentEnergy, int nPermanent,
            double polarizationEnergy, int nPolar) {
        this.ciOnly = ciOnly;
        this.info = info;
        this.bondEnergy = bondEnergy;
        this.nBonds = nBonds;
        this.angleEnergy = angleEnergy;
        this.nAngles = nAngles;
        this.stretchBendEnergy = stretchBendEnergy;
        this.nStretchBends = nStretchBends;
        this.ureyBradleyEnergy = ureyBradleyEnergy;
        this.nUreyBradleys = nUreyBradleys;
        this.outOfPlaneBendEnergy = outOfPlaneBendEnergy;
        this.nOutOfPlaneBends = nOutOfPlaneBends;
        this.torsionEnergy = torsionEnergy;
        this.nTorsions = nTorsions;
        this.improperTorsionEnergy = improperTorsionEnergy;
        this.nImproperTorsions = nImproperTorsions;
        this.piOrbitalTorsionEnergy = piOrbitalTorsionEnergy;
        this.nPiOrbitalTorsions = nPiOrbitalTorsions;
        this.torsionTorsionEnergy = torsionTorsionEnergy;
        this.nTorsionTorsions = nTorsionTorsions;
        this.vanDerWaalsEnergy = vanDerWaalsEnergy;
        this.nVanDerWaals = nVanDerWaals;
        this.permanentEnergy = permanentEnergy;
        this.nPermanent = nPermanent;
        this.polarizationEnergy = polarizationEnergy;
        this.nPolar = nPolar;

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
        // Stretch-Bend Energy
        assertEquals(info + " Stretch-Bend Energy", stretchBendEnergy, forceFieldEnergy.getStrenchBendEnergy(), tolerance);
        assertEquals(info + " Stretch-Bend Count", nStretchBends, forceFieldEnergy.getNumberofStretchBends());
        // Urey-Bradley Energy
        assertEquals(info + " Urey-Bradley Energy", ureyBradleyEnergy, forceFieldEnergy.getUreyBradleyEnergy(), tolerance);
        assertEquals(info + " Urey-Bradley Count", nUreyBradleys, forceFieldEnergy.getNumberofUreyBradleys());
        // Out-of-Plane Bend
        assertEquals(info + " Out-of-Plane Bend Energy", outOfPlaneBendEnergy, forceFieldEnergy.getOutOfPlaneBendEnergy(), tolerance);
        assertEquals(info + " Out-of-Plane Bend Count", nOutOfPlaneBends, forceFieldEnergy.getNumberofOutOfPlaneBends());
        // Torsional Angle
        assertEquals(info + " Torsion Energy", torsionEnergy, forceFieldEnergy.getTorsionEnergy(), tolerance);
        assertEquals(info + " Torsion Count", nTorsions, forceFieldEnergy.getNumberofTorsions());
        // Improper Torsional Angle
        assertEquals(info + " Improper Torsion Energy", improperTorsionEnergy, forceFieldEnergy.getImproperTorsionEnergy(), tolerance);
        assertEquals(info + " Improper Torsion Count", nImproperTorsions, forceFieldEnergy.getNumberofImproperTorsions());
        // Pi-Orbital Torsion
        assertEquals(info + " Pi-OrbitalTorsion Energy", piOrbitalTorsionEnergy,
                forceFieldEnergy.getPiOrbitalTorsionEnergy(), tolerance);
        assertEquals(info + " Pi-OrbitalTorsion Count", nPiOrbitalTorsions, forceFieldEnergy.getNumberofPiOrbitalTorsions());
        // Torsion-Torsion
        assertEquals(info + " Torsion-Torsion Energy", torsionTorsionEnergy, forceFieldEnergy.getTorsionTorsionEnergy(), tolerance);
        assertEquals(info + " Torsion-Torsion Count", nTorsionTorsions, forceFieldEnergy.getNumberofTorsionTorsions());
        // van Der Waals
        assertEquals(info + " van Der Waals Energy", vanDerWaalsEnergy, forceFieldEnergy.getVanDerWaalsEnergy(), tolerance);
        assertEquals(info + " van Der Waals Count", nVanDerWaals, forceFieldEnergy.getVanDerWaalsInteractions());
        // Permanent Multipoles
        if (mpoleTerm) {
            assertEquals(info + " Permanent Multipole Energy", permanentEnergy,
                    forceFieldEnergy.getPermanentMultipoleEnergy(), tolerance);
            assertEquals(info + " Permanent Multipole Count", nPermanent, forceFieldEnergy.getPermanentInteractions());
        }
        // Polarization
        assertEquals(info + " Polarization Energy", polarizationEnergy, forceFieldEnergy.getPolarizationEnergy(), tolerance);
        assertEquals(info + " Polarization Count", nPolar, forceFieldEnergy.getPermanentInteractions());

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
        double step = 0.00001;
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
