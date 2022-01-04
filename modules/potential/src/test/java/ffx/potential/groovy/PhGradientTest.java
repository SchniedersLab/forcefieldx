// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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

import static org.junit.Assert.assertEquals;

import ffx.potential.utils.PotentialTest;
import ffx.potential.groovy.test.PhGradient;
import java.util.Arrays;
import java.util.Collection;
import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

/** JUnit Tests for the PhGradient Script */

@RunWith(Parameterized.class)
public class PhGradientTest extends PotentialTest {

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
      double acidostatEnergy) {
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
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "\n Titration Lambdas: 0.0000, 0.0000, 0.0000, 0.0000, " +
                        "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.36542345,
                6016,
                -110.21708124,
                5478,
                -28.24011977,
                5478,
                0.00000000,
                -20.56768333
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 0.0000, 0.0000, 0.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.96855658,
                6016,
                -176.50397580,
                5691,
                -15.18071794,
                5691,
                0.00000000,
                -80.79416948
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 0.0000, 1.0000, 0.0000, 0.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.05692022,
                6016,
                -170.33703658,
                5691,
                -17.82995367,
                5691,
                0.00000000,
                -86.86494229
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 1.0000, 0.0000, 0.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.65993599,
                6016,
                -185.41105516,
                5908,
                -8.40685849,
                5908,
                0.00000000,
                -147.09142844
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 0.0000, 0.0000, 1.0000, 0.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.31240247,
                6016,
                -165.84274841,
                5478,
                -32.92926880,
                5478,
                0.00000000,
                13.28255800
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 0.0000, 1.0000, 0.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.91553489,
                6016,
                -208.18172444,
                5691,
                -19.10501950,
                5691,
                0.00000000,
                -46.94392815
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 0.0000, 1.0000, 1.0000, 0.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.00386493,
                6016,
                -184.91785146,
                5691,
                -24.96658110,
                5691,
                0.00000000,
                -53.01470096
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 1.0000, 1.0000, 0.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.60687999,
                6016,
                -176.04395152,
                5908,
                -14.81292194,
                5908,
                0.00000000,
                -113.24118712
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 0.0000, 0.0000, 0.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.62854247,
                6016,
                -124.60312166,
                5582,
                -36.41919614,
                5582,
                0.00000000,
                20.64654243
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 0.0000, 0.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.23167558,
                6016,
                -175.24747721,
                5797,
                -22.40684819,
                5797,
                0.00000000,
                -39.57994372
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 0.0000, 1.0000, 0.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.32003898,
                6016,
                -162.41187049,
                5797,
                -25.39207665,
                5797,
                0.00000000,
                -45.65071653
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 1.0000, 0.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.92305472,
                6016,
                -161.84335007,
                6016,
                -15.02343905,
                6016,
                0.00000000,
                -105.87720268
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 0.0000, 0.0000, 1.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                    27.29020028,
                    113,
                    17.85143041,
                    199,
                    -1.35218784,
                    168,
                    0.0,
                    0,
                    0.00017447,
                    69,
                    26.86603133,
                    282,
                    0.0,
                    0,
                    0.00000716,
                    13,
                0.0,
                0,
                85.57543191,
                6016,
                -134.50074439,
                5582,
                -42.26053086,
                5582,
                0.00000000,
                54.49678376
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 0.0000, 1.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.17856431,
                6016,
                -161.19718142,
                5797,
                -27.45638385,
                5797,
                0.00000000,
                -5.72970239
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 0.0000, 1.0000, 1.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.26689412,
                6016,
                -131.26464093,
                5797,
                -33.62322895,
                5797,
                0.00000000,
                -11.80047520
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                    "\n Titration Lambdas: 1.0000, 1.0000, 1.0000, 1.0000, " +
                            "\n Tautomer Lambdas: 0.5000, 0.5000, 0.5000",
                27.29020028,
                113,
                17.85143041,
                199,
                -1.35218784,
                168,
                0.0,
                0,
                0.00017447,
                69,
                26.86603133,
                282,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.86990915,
                6016,
                -106.74820198,
                6016,
                -22.49713513,
                6016,
                0.00000000,
                -72.02696135
            }

        });
  }

  @Test
  public void testEndStateEnergyAndGradient() {
    // Configure input arguments for the PhGradient script.
    logger.info(" Testing endstate energy for " + info + "at lambdas: " + key);

    // Choose a random Atom (1..111) to test XYZ gradient.
    Random random = new Random();
    String randomAtom = Integer.toString(random.nextInt(110) + 1);

    String[] args = {"-v", "--esvLambda", "0.5", "-d", "0.000001", "--ga", randomAtom, filename};
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
    assertEquals(info + " Stretch-Bend Energy", stretchBendEnergy, energyAndInteractionList[4],
        tolerance);
    assertEquals(info + " Stretch-Bend Count", nStretchBends, (int) energyAndInteractionList[5]);
    // Urey-Bradley Energy
    assertEquals(info + " Urey-Bradley Energy", ureyBradleyEnergy, energyAndInteractionList[6],
        tolerance);
    assertEquals(info + " Urey-Bradley Count", nUreyBradleys, (int) energyAndInteractionList[7]);
    // Out-of-Plane Bend
    assertEquals(info + " Out-of-Plane Bend Energy", outOfPlaneBendEnergy,
        energyAndInteractionList[8], tolerance);
    assertEquals(info + " Out-of-Plane Bend Count", nOutOfPlaneBends,
        (int) energyAndInteractionList[9]);
    // Torsional Angle
    assertEquals(info + " Torsion Energy", torsionEnergy, energyAndInteractionList[10], tolerance);
    assertEquals(info + " Torsion Count", nTorsions, (int) energyAndInteractionList[11]);
    // Improper Torsional Angle
    assertEquals(info + " Improper Torsion Energy", improperTorsionEnergy,
        energyAndInteractionList[12], tolerance);
    assertEquals(info + " Improper Torsion Count", nImproperTorsions,
        (int) energyAndInteractionList[13]);
    // Pi-Orbital Torsion
    assertEquals(info + " Pi-OrbitalTorsion Energy", piOrbitalTorsionEnergy,
        energyAndInteractionList[14], tolerance);
    assertEquals(info + " Pi-OrbitalTorsion Count", nPiOrbitalTorsions,
        (int) energyAndInteractionList[15]);
    // Torsion-Torsion
    assertEquals(info + " Torsion-Torsion Energy", torsionTorsionEnergy,
        energyAndInteractionList[16], tolerance);
    assertEquals(info + " Torsion-Torsion Count", nTorsionTorsions,
        (int) energyAndInteractionList[17]);
    // van Der Waals
    assertEquals(info + " van Der Waals Energy", vanDerWaalsEnergy, energyAndInteractionList[18],
        tolerance);
    assertEquals(info + " van Der Waals Count", nVanDerWaals, (int) energyAndInteractionList[19]);
    // Permanent Multipoles
    assertEquals(info + " Permanent Multipole Energy", permanentEnergy, energyAndInteractionList[20],
        tolerance);
    assertEquals(info + " Permanent Multipole Count", nPermanent,
        (int) energyAndInteractionList[21]);
    // Polarization Energy
    assertEquals(info + " Polarization Energy", polarizationEnergy, energyAndInteractionList[22],
        tolerance);
    assertEquals(info + " Polarization Count", nPolar, (int) energyAndInteractionList[23]);
    // Extended System Bias
    assertEquals(info + " ExtendedSystemBias", extendedSystemBias, energyAndInteractionList[24],
        tolerance);
    // Total Energy
    assertEquals(info + " Total Energy", totalEnergy, energyAndInteractionList[25], tolerance);

    // Check for a gradient failure.
    assertEquals("DEHK gradient failures: ", 0, pHGradient.nFailures);
    assertEquals("DEHK ESV gradient failures: ", 0, pHGradient.nESVFailures);
  }
}

