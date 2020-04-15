/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2020.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.potential.groovy;

import java.util.Arrays;
import java.util.Collection;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.random;
import static org.junit.Assert.assertEquals;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.groovy.test.Gradient;
import ffx.potential.groovy.test.LambdaGradient;
import ffx.utilities.BaseFFXTest;

import groovy.lang.Binding;

/**
 * Test the Energy script.
 */
@RunWith(Parameterized.class)
public class EnergyTest extends BaseFFXTest {

    private final String info;
    private final String filename;
    private final int nAtoms;
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
    private final int nGK;
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
    private final double gkEnergy;
    private final double totalEnergy;
    private final boolean testOpenMM;
    private final double tolerance = 1.0e-2;
    private final double openMMTolerance = 0.5;
    private Binding binding;
    private Energy energy;
    private Gradient gradient;
    private LambdaGradient lambdaGradient;

    public EnergyTest(String info, String filename, int nAtoms,
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
                      double polarizationEnergy, int nPolar,
                      double gkEnergy, int nGK, boolean testOpenMM) {
        this.filename = filename;
        this.info = info;
        this.nAtoms = nAtoms;
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
        this.gkEnergy = gkEnergy;
        this.nGK = nGK;
        this.testOpenMM = testOpenMM;

        totalEnergy = bondEnergy + angleEnergy + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy
                + torsionEnergy + improperTorsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy
                + vanDerWaalsEnergy + permanentEnergy + polarizationEnergy + gkEnergy;
    }

    @Before
    public void before() {
        binding = new Binding();
        System.clearProperty("platform");
    }

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {
                        "SNARE P212121",
                        "ffx/potential/structures/1n7s.P212121.xyz",
                        5357,
                        351.32142483, 5040,
                        744.19251364, 8255,
                        6.42015744, 6924,
                        2.50163831, 322,
                        135.24919366, 2487,
                        417.89244418, 11449,
                        0.0, 0,
                        39.85643934, 370,
                        -560.99576469, 268,
                        4003.27183547, 741643,
                        -12303.93157019, 332114,
                        -2811.45683671, 332114,
                        0.0, 0, false},
                {
                        // OpenMM does not handle this correctly yet.
                        "SNARE P1",
                        "ffx/potential/structures/1n7s.P1.xyz",
                        21428,
                        1405.28569930, 20160,
                        2976.77005458, 33020,
                        25.68062976, 27696,
                        10.00655326, 1288,
                        540.99677465, 9948,
                        1671.56977674, 45796,
                        0.0, 0,
                        159.42575736, 1480,
                        -2243.98305878, 1072,
                        16013.08734188, 2966572,
                        -49215.72628076, 1328456,
                        -11245.82734685, 1328456,
                        0.0, 0, false},
                {
                        "DHFR Benchmark",
                        "ffx/potential/structures/dhfr.xyz",
                        23558,
                        6623.70878543, 16569,
                        4748.78155255, 11584,
                        -19.78973386, 4031,
                        -132.74781407, 7023,
                        110.84355342, 1566,
                        433.69593807, 6701,
                        0.0, 0,
                        58.75073806, 292,
                        -41.99981847, 147,
                        31718.84803063, 3480619,
                        -78751.36662993, 1463589,
                        -31790.15161303, 1463589,
                        0.0, 0, true},
                {
                        "AMBER99SB GB (no dispersion) Capped DMHD",
                        "ffx/potential/structures/dmhd-amber99sb.xyz",
                        71,
                        1.56331971, 71,
                        25.07578077, 124,
                        0.0, 0,
                        0.0, 0,
                        0.0, 0,
                        27.92653816, 169,
                        0.15162117, 20,
                        0.0, 0,
                        0.0, 0,
                        -4.31922323, 2290,
                        -71.00737570, 2290,
                        0.0, 2290,
                        -161.32855355, 2556, true},
                {
                        "AMOEBA Bio 2018 Crambin GK/SA",
                        "ffx/potential/structures/crambin.xyz",
                        642,
                        134.83946761, 652,
                        133.94053912, 1183,
                        -4.81704334, 1044,
                        0.0, 0,
                        29.92634330, 375,
                        18.82386314, 1747,
                        0.0, 0,
                        12.01761948, 66,
                        -23.33097186, 4,
                        312.82576484, 203926,
                        -1208.54202510, 203926,
                        -127.41736039, 203926,
                        -237.20593090, 206403, false},
                {
                        "AMOEBA Protein 2013 GK Capped DMHD",
                        "ffx/potential/structures/dmhd-amoebapro13.xyz",
                        71,
                        4.00030221, 71,
                        15.27574588, 124,
                        -0.23989418, 110,
                        0.0, 0,
                        0.39337245, 51,
                        3.40745668, 169,
                        0.0, 0,
                        0.09663007, 10,
                        0.0, 0,
                        20.45106694, 2290,
                        -169.03958434, 2290,
                        -25.85878386, 2290,
                        -197.35407295, 2556, true},
                {
                        "AMOEBA 2009 Tetra-Alanine GK",
                        "ffx/potential/structures/alatet.int",
                        42,
                        1.66598227, 41,
                        3.41765862, 72,
                        -0.07860298, 57,
                        0.0, 0,
                        0.29369677, 24,
                        7.99325266, 91,
                        0.0, 0,
                        0.01292169, 4,
                        1.95288364, 3,
                        20.39805900, 748,
                        -78.53945229, 748,
                        -8.62910979, 748,
                        -58.7112682143682, 903, false},
                {
                        "Amber99sb Peptide",
                        "ffx/potential/structures/peptide-amber99sb.xyz",
                        328,
                        41.19699756, 333,
                        54.88751601, 596,
                        0.0, 0,
                        0.0, 0,
                        0.0, 0,
                        189.34195768, 875,
                        1.59062194, 101,
                        0.0, 0,
                        0.0, 0,
                        111362.79687915, 52696,
                        -413.52897447, 52699,
                        0.0, 52699,
                        0.0, 0, true},
                {
                        "OPLS-AA Peptide",
                        "ffx/potential/structures/peptide-oplsaa.xyz",
                        328,
                        72.08575480, 333,
                        32.38121260, 596,
                        0.0, 0,
                        0.0, 0,
                        0.0, 0,
                        65.56283245, 875,
                        0.85331537, 186,
                        0.0, 0,
                        0.0, 0,
                        91224.14951970, 40511,
                        -665.41688158, 52699,
                        0.0, 52699,
                        0.0, 0, false},
                {
                        "OPLS-AA/L Peptide",
                        "ffx/potential/structures/peptide-oplsaal.xyz",
                        328,
                        39.69175722, 333,
                        41.54908436, 596,
                        0.0, 0,
                        0.0, 0,
                        0.0, 0,
                        41.53603210, 875,
                        0.92410847, 97,
                        0.0, 0,
                        0.0, 0,
                        112122.04255274, 40511,
                        -671.66812023, 52699,
                        0.0, 52699,
                        0.0, 0, false},
                {
                        "Ubiquitin Benchmark",
                        "ffx/potential/structures/ubiquitin.xyz",
                        9737,
                        2856.59378034, 6908,
                        2018.16571190, 5094,
                        -6.98794687, 1958,
                        -55.63461732, 2835,
                        65.00809012, 651,
                        224.37608631, 3297,
                        0.0, 0,
                        20.43063520, 106,
                        -31.48011891, 71,
                        12862.34983956, 1482946,
                        -32736.92207773, 623627,
                        -12934.93512829, 623627,
                        0.0, 0, true},
                {
                        "Ubiquitin Amber 1999 Benchmark",
                        "ffx/potential/structures/ubiquitin-amber99.xyz",
                        9737,
                        3496.26834781, 6908,
                        2056.15089674, 5094,
                        0.0, 0,
                        0.0, 0,
                        0.0, 0,
                        705.79420086, 3297,
                        27.59366005, 287,
                        0.0, 0,
                        0.0, 0,
                        4860.41072437, 2279196,
                        -40036.04804903, 3562034,
                        0.0, 3562034,
                        0.0, 0, true},
                {
                        "Acetanilide Benchmark",
                        "ffx/potential/structures/acetanilide.xyz",
                        19,
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
                        -1.9856761645305, 2231,
                        0.0, 0, false},
                {
                        "Ethylparaben Benchmark",
                        "ffx/potential/structures/ethylparaben.xyz",
                        44,
                        1.37991170, 44,
                        2.43600964, 70,
                        0.20006782, 58,
                        0.0, 0,
                        0.04434686, 42,
                        -7.14298308, 88,
                        0.0, 0,
                        0.0, 0,
                        0.0, 0,
                        -4.816366774263829, 32644,
                        -45.55012417320137, 5012,
                        -3.97079858673400, 5012,
                        0.0, 0, false},
                {
                        "Methylparaben Benchmark",
                        "ffx/potential/structures/methylparaben.xyz",
                        19,
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
                        -2.35904592, 2314,
                        0.0, 0, false},
                {
                        "Paracetamol Benchmark",
                        "ffx/potential/structures/paracetamol.xyz",
                        20,
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
                        -4.96649992387851, 2357,
                        0.0, 0, false},
                {
                        "Phenacetin Benchmark",
                        "ffx/potential/structures/phenacetin.xyz",
                        26,
                        0.56492623, 26,
                        2.68273715, 43,
                        0.00409248, 35,
                        0.0, 0,
                        0.33385711, 24,
                        -6.49787814, 52,
                        0.0, 0,
                        0.0, 0,
                        0.0, 0,
                        -5.792072251098926, 20384,
                        -25.42173449818894, 3138,
                        -2.09385416786188, 3138,
                        0.0, 0, false}
        });
    }

    @Test
    public void testEnergy() {
        if (nAtoms > 10000 && !ffxCI) {
            return;
        }
        logger.info(" Testing energy for " + info);

        energy = new Energy();
        energy.setBinding(binding);

        // Set-up the input arguments for the Energy script.
        String[] args = {"src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate the script.
        energy.run();

        ForceFieldEnergy forceFieldEnergy = energy.forceFieldEnergy;

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
        assertEquals(info + " Pi-OrbitalTorsion Energy", piOrbitalTorsionEnergy, forceFieldEnergy.getPiOrbitalTorsionEnergy(), tolerance);
        assertEquals(info + " Pi-OrbitalTorsion Count", nPiOrbitalTorsions, forceFieldEnergy.getNumberofPiOrbitalTorsions());
        // Torsion-Torsion
        assertEquals(info + " Torsion-Torsion Energy", torsionTorsionEnergy, forceFieldEnergy.getTorsionTorsionEnergy(), tolerance);
        assertEquals(info + " Torsion-Torsion Count", nTorsionTorsions, forceFieldEnergy.getNumberofTorsionTorsions());
        // van Der Waals
        assertEquals(info + " van Der Waals Energy", vanDerWaalsEnergy, forceFieldEnergy.getVanDerWaalsEnergy(), tolerance);
        assertEquals(info + " van Der Waals Count", nVanDerWaals, forceFieldEnergy.getVanDerWaalsInteractions());
        // Permanent Multipoles
        assertEquals(info + " Permanent Multipole Energy", permanentEnergy, forceFieldEnergy.getPermanentMultipoleEnergy(), tolerance);
        assertEquals(info + " Permanent Multipole Count", nPermanent, forceFieldEnergy.getPermanentInteractions());
        // Polarization Energy
        assertEquals(info + " Polarization Energy", polarizationEnergy, forceFieldEnergy.getPolarizationEnergy(), tolerance);
        assertEquals(info + " Polarization Count", nPolar, forceFieldEnergy.getPermanentInteractions());
        // GK Energy
        assertEquals(info + " Solvation", gkEnergy, forceFieldEnergy.getSolvationEnergy(), tolerance);
        assertEquals(info + " Solvation Count", nGK, forceFieldEnergy.getSolvationInteractions());
        energy.destroyPotentials();
        System.gc();
    }

    @Test
    public void testGradient() {
        if (nAtoms > 5000 && !ffxCI) {
            return;
        }
        logger.info(" Testing Cartesian gradient(s) for " + info);

        gradient = new Gradient();
        gradient.setBinding(binding);

        // Choose a random atom to test.
        int atomID = (int) floor(random() * nAtoms) + 1;
        double stepSize = 1.0e-5;
        // Set-up the input arguments for the Gradient script.
        String[] args = {"--ga", Integer.toString(atomID),
                "--dx", Double.toString(stepSize),
                "--tol", Double.toString(tolerance),
                "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate the script.
        gradient.run();

        assertEquals(info + " gradient failures: ", 0, gradient.nFailures);
        gradient.destroyPotentials();
        System.gc();
    }

    @Test
    public void testLambdaGradient() {
        if (nAtoms > 5000 && !ffxCI) {
            return;
        }
        logger.info(" Testing lambda gradient(s) for " + info);

        lambdaGradient = new LambdaGradient();
        lambdaGradient.setBinding(binding);

        // Choose a random atom to test dEdX gradient.
        int atomID = (int) floor(random() * nAtoms) + 1;

        double stepSize = 1.0e-5;
        // Set-up the input arguments for the Lambda Gradient script.
        String[] args = {"--ga", Integer.toString(atomID),
                "--dx", Double.toString(stepSize),
                "--tol", Double.toString(tolerance),
                "--ac", "1" + "-" + nAtoms,
                "-l", "0.5",
                "src/main/java/" + filename};
        binding.setVariable("args", args);

        // Evaluate the script.
        lambdaGradient.run();

        System.clearProperty("lambdaterm");

        assertEquals(info + " dEdL failures: ", 0, lambdaGradient.ndEdLFailures);
        assertEquals(info + " d2EdL2 failures: ", 0, lambdaGradient.nd2EdL2Failures);
        assertEquals(info + " dEdXdL failures: ", 0, lambdaGradient.ndEdXdLFailures);
        assertEquals(info + " dEdX failures: ", 0, lambdaGradient.ndEdXFailures);
        lambdaGradient.destroyPotentials();
        System.gc();
    }

    @Test
    public void testOpenMMEnergy() {
        if (!testOpenMM || !ffxOpenMM) {
            return;
        }
        logger.info(" Testing OpenMM energy for " + info);

        energy = new Energy();
        energy.setBinding(binding);

        // Set-up the input arguments for the Biotype script.
        String[] args = {"src/main/java/" + filename};
        binding.setVariable("args", args);

        System.setProperty("platform", "OMM");

        // Evaluate the script.
        energy.run();

        System.clearProperty("platform");

        assertEquals(info + " OpenMM Energy", totalEnergy, energy.energy, openMMTolerance);
        energy.destroyPotentials();
        System.gc();
    }
}
