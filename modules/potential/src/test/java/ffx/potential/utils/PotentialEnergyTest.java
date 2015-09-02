/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
 * Test the PotentialEnergy class.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@RunWith(Parameterized.class)
public class PotentialEnergyTest {

    private static final Logger logger = Logger.getLogger(PotentialEnergyTest.class.getName());

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
            {
                false,
                "Ubiquitin Benchmark",
                "ffx/potential/structures/ubiquitin.xyz",
                2673.37683484, 6908,
                1637.34919041, 5094,
                -11.04350364, 1958,
                279.64162198, 2835,
                67.64798284, 651,
                215.14214012, 3297,
                24.69060350, 106,
                -29.43681349, 71,
                13183.92864934, 1483768,
                -33012.66179952, 623490,
                -13041.30955459, 623490},
            {true,
                "DHFR Benchmark",
                "ffx/potential/structures/dhfr.xyz",
                6423.84579926, 16569,
                3746.31506290, 11584,
                -21.85553039, 4031,
                687.46861123, 7023,
                198.72886589, 1566,
                426.23738971, 6701,
                48.26628393, 292,
                -41.71473465, 147,
                32630.94057333, 3480445,
                -79396.71166429, 1463353,
                -32141.39930772, 1463353},
            {true,
                "SNARE P1",
                "ffx/potential/structures/1n7s.P1.xyz",
                1405.28569930, 20160,
                2976.77005458, 33020,
                25.68062976, 27696,
                10.00655326, 1288,
                540.99677465, 9948,
                1671.56977674, 45796,
                159.42575736, 1480,
                -2243.98305878, 1072,
                16013.08734188, 2966572,
                -49215.72628076, 1328456,
                -11245.82734685, 1328456},
            {true,
                "SNARE P212121",
                "ffx/potential/structures/1n7s.P212121.xyz",
                351.32142483, 5040,
                744.19251364, 8255,
                6.42015744, 6924,
                2.50163831, 322,
                135.24919366, 2487,
                417.89244418, 11449,
                39.85643934, 370,
                -560.99576469, 268,
                4003.27183547, 741643,
                -12303.93157019, 332114,
                -2811.45683671, 332114}});
    }
    private final String info;
    private final File structure;
    private final MolecularAssembly molecularAssembly;
    private final int nAtoms;
    private final int nBonds;
    private final int nAngles;
    private final int nStretchBends;
    private final int nUreyBradleys;
    private final int nOutOfPlaneBends;
    private final int nTorsions;
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
    private final double piOrbitalTorsionEnergy;
    private final double torsionTorsionEnergy;
    private final double vanDerWaalsEnergy;
    private final double permanentEnergy;
    private final double polarizationEnergy;
    private final double tolerance = 1.0e-3;
    private final double gradientTolerance = 1.0e-4;
    private final boolean ci;
    private final boolean ciOnly;
    private boolean mpoleTerm;
    private final ForceFieldEnergy forceFieldEnergy;

    public PotentialEnergyTest(
            boolean ciOnly,
            String info, String filename,
            double bondEnergy, int nBonds,
            double angleEnergy, int nAngles,
            double stretchBendEnergy, int nStretchBends,
            double ureyBradleyEnergy, int nUreyBradleys,
            double outOfPlaneBendEnergy, int nOutOfPlaneBends,
            double torsionEnergy, int nTorsions,
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

        String ffxCi = System.getProperty("ffx.ci", "false");
        ci = ffxCi.equalsIgnoreCase("true");

        if (!ci && ciOnly) {
            structure = null;
            molecularAssembly = null;
            nAtoms = 0;
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

        nAtoms = molecularAssembly.getAtomArray().length;
        forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        mpoleTerm = molecularAssembly.getForceField().getBoolean(ForceField.ForceFieldBoolean.MPOLETERM, true);

        /*
         if (ci) {
         testGradient();
         testSoftCore();
         } */
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
