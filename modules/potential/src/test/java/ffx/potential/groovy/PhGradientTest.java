// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.groovy;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.groovy.test.Gradient;
import ffx.potential.groovy.test.PhGradient;
import ffx.potential.utils.PotentialTest;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.random;
import static org.junit.Assert.assertEquals;

/** JUnit Tests for the PhGradient Script */

@RunWith(Parameterized.class)
public class PhGradientTest extends PotentialTest{
    private final String info;
    private final String filename;
    private final String key;
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
    private final double discretizerEnergy;
    private final double acidostatEnergy;
    private final double extendedSystemBias;
    private final double totalEnergy;
    private final double tolerance = 1.0e-2;

    public PhGradientTest(
            String info,
            String filename,
            String key,
            double bondEnergy,
            int nBonds,
            double angleEnergy,
            int nAngles,
            double stretchBendEnergy,
            int nStretchBends,
            double ureyBradleyEnergy,
            int nUreyBradleys,
            double outOfPlaneBendEnergy,
            int nOutOfPlaneBends,
            double torsionEnergy,
            int nTorsions,
            double improperTorsionEnergy,
            int nImproperTorsions,
            double piOrbitalTorsionEnergy,
            int nPiOrbitalTorsions,
            double torsionTorsionEnergy,
            int nTorsionTorsions,
            double vanDerWaalsEnergy,
            int nVanDerWaals,
            double permanentEnergy,
            int nPermanent,
            double polarizationEnergy,
            int nPolar,
            double discretizerEnergy,
            double acidostatEnergy){
        this.filename = filename;
        this.info = info;
        this.key = key;
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
        this.discretizerEnergy = discretizerEnergy;
        this.acidostatEnergy = acidostatEnergy;

        extendedSystemBias = discretizerEnergy + acidostatEnergy;

        totalEnergy =
                bondEnergy
                        + angleEnergy
                        + stretchBendEnergy
                        + ureyBradleyEnergy
                        + outOfPlaneBendEnergy
                        + torsionEnergy
                        + improperTorsionEnergy
                        + piOrbitalTorsionEnergy
                        + torsionTorsionEnergy
                        + vanDerWaalsEnergy
                        + permanentEnergy
                        + polarizationEnergy
                        + discretizerEnergy
                        + acidostatEnergy;
    }

    @Parameterized.Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(
                new Object[][] {
                        {
                                "Lys-Lys Dipeptide",
                                "src/main/java/ffx/potential/structures/lys-lys.pdb",
                                "0.0000, 0.0000",
                                3.54438920,
                                105,
                                4.92151429,
                                190,
                                -0.34239425,
                                161,
                                0.0,
                                0,
                                0.00008300,
                                28,
                                7.05487554,
                                271,
                                0.0,
                                0,
                                0.00000080,
                                4,
                                0.0,
                                0,
                                38.25724479,
                                1960,
                                -56.53581530,
                                2145,
                                -2.58637980,
                                2145,
                                0.00000000,
                                0.00000000
                        },
                        {
                                "Lys-Lys Dipeptide",
                                "src/main/java/ffx/potential/structures/lys-lys.pdb",
                                "1.0000, 0.0000",
                                3.55340652,
                                105,
                                5.09658486,
                                190,
                                -0.34819224,
                                161,
                                0.0,
                                0,
                                0.00008300,
                                28,
                                6.99474604,
                                271,
                                0.0,
                                0,
                                0.00000080,
                                4,
                                0.0,
                                0,
                                38.51986213,
                                1960,
                                -38.93435694,
                                2145,
                                -9.00775857,
                                2145,
                                0.00000000,
                                -46.16990685
                        },
                        {
                                "Lys-Lys Dipeptide",
                                "src/main/java/ffx/potential/structures/lys-lys.pdb",
                                "0.0000, 1.0000",
                                3.55408365,
                                105,
                                5.09836859,
                                190,
                                -0.34831497,
                                161,
                                0.0,
                                0,
                                0.00008300,
                                28,
                                6.99477443,
                                271,
                                0.0,
                                0,
                                0.00000080,
                                4,
                                0.0,
                                0,
                                38.51992476,
                                1960,
                                -31.63779852,
                                2145,
                                -9.95024270,
                                2145,
                                0.00000000,
                                -46.16990685
                        },
                        {
                                "Lys-Lys Dipeptide",
                                "src/main/java/ffx/potential/structures/lys-lys.pdb",
                                "1.0000, 1.0000",
                                3.56310097,
                                105,
                                5.27343916,
                                190,
                                -0.35411297,
                                161,
                                0.0,
                                0,
                                0.00008300,
                                28,
                                6.93464493,
                                271,
                                0.0,
                                0,
                                0.00000080,
                                4,
                                0.0,
                                0,
                                38.78251486,
                                1960,
                                33.73717085,
                                2145,
                                -18.42522652,
                                2145,
                                0.00000000,
                                -92.33981370
                        }

                });
    }


    @Test
    public void testEndstateEnergy() {
        // Configure input arguments for the PhGradient script.
        logger.info(" Testing endstate energy for " + info + "at lambdas: " + key);
        String[] args = {"-v", filename};
        binding.setVariable("args", args);
        // Construct and evaluate the Volume script.
        PhGradient pHGradient = new PhGradient(binding).run();
        potentialScript = pHGradient;
        
        double[] energyAndInteractionList = pHGradient.endstateEnergyMap.get(key);
        // Bond Energy
        assertEquals(info + " Bond Energy", bondEnergy, energyAndInteractionList[0], tolerance);
        assertEquals(info + " Bond Count", nBonds, (int) energyAndInteractionList[1]);
        // Angle Energy
        assertEquals(info + " Angle Energy", angleEnergy, energyAndInteractionList[2], tolerance);
        assertEquals(info + " Angle Count", nAngles, (int) energyAndInteractionList[3]);
        // Stretch-Bend Energy
        assertEquals(
                info + " Stretch-Bend Energy",
                stretchBendEnergy,
                energyAndInteractionList[4],
                tolerance);
        assertEquals(
                info + " Stretch-Bend Count", nStretchBends, (int) energyAndInteractionList[5]);
        // Urey-Bradley Energy
        assertEquals(
                info + " Urey-Bradley Energy",
                ureyBradleyEnergy,
                energyAndInteractionList[6],
                tolerance);
        assertEquals(
                info + " Urey-Bradley Count", nUreyBradleys, (int) energyAndInteractionList[7]);
        // Out-of-Plane Bend
        assertEquals(
                info + " Out-of-Plane Bend Energy",
                outOfPlaneBendEnergy,
                energyAndInteractionList[8],
                tolerance);
        assertEquals(
                info + " Out-of-Plane Bend Count",
                nOutOfPlaneBends,
                (int) energyAndInteractionList[9]);
        // Torsional Angle
        assertEquals(
                info + " Torsion Energy", torsionEnergy, energyAndInteractionList[10], tolerance);
        assertEquals(info + " Torsion Count", nTorsions, (int) energyAndInteractionList[11]);
        // Improper Torsional Angle
        assertEquals(
                info + " Improper Torsion Energy",
                improperTorsionEnergy,
                energyAndInteractionList[12],
                tolerance);
        assertEquals(
                info + " Improper Torsion Count",
                nImproperTorsions,
                (int) energyAndInteractionList[13]);
        // Pi-Orbital Torsion
        assertEquals(
                info + " Pi-OrbitalTorsion Energy",
                piOrbitalTorsionEnergy,
                energyAndInteractionList[14],
                tolerance);
        assertEquals(
                info + " Pi-OrbitalTorsion Count",
                nPiOrbitalTorsions,
                (int) energyAndInteractionList[15]);
        // Torsion-Torsion
        assertEquals(
                info + " Torsion-Torsion Energy",
                torsionTorsionEnergy,
                energyAndInteractionList[16],
                tolerance);
        assertEquals(
                info + " Torsion-Torsion Count",
                nTorsionTorsions,
                (int) energyAndInteractionList[17]);
        // van Der Waals
        assertEquals(
                info + " van Der Waals Energy",
                vanDerWaalsEnergy,
                energyAndInteractionList[18],
                tolerance);
        assertEquals(
                info + " van Der Waals Count", nVanDerWaals, (int) energyAndInteractionList[19]);
        // Permanent Multipoles
        assertEquals(
                info + " Permanent Multipole Energy",
                permanentEnergy,
                energyAndInteractionList[20],
                tolerance);
        assertEquals(
                info + " Permanent Multipole Count",
                nPermanent,
                (int) energyAndInteractionList[21]);
        // Polarization Energy
        assertEquals(
                info + " Polarization Energy",
                polarizationEnergy,
                energyAndInteractionList[22],
                tolerance);
        assertEquals(info + " Polarization Count", nPolar, (int) energyAndInteractionList[23]);
        // Extended System Bias
        assertEquals(info + " ExtendedSystemBias", extendedSystemBias, energyAndInteractionList[24],tolerance);
        // Total Energy
        assertEquals(info + " Total Energy", totalEnergy, energyAndInteractionList[25],tolerance);
    }

    @Test
    public void testGradients() {
        logger.info(" Testing gradients for " + info + "at lambdas: " + key);
        String[] args = {"-v", filename};
        binding.setVariable("args", args);
        // Construct and evaluate the PhGradient script.
        PhGradient pHGradient = new PhGradient(binding).run();
        potentialScript = pHGradient;
        assertEquals("Lys-Lys gradient failures: ", 0, pHGradient.nFailures);
        assertEquals("Lys-Lys ESV gradient failures: ", 0, pHGradient.nESVFailures);
    }
}

