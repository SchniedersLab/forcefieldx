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

import static org.junit.Assert.assertEquals;

import ffx.potential.groovy.test.PhGradient;
import ffx.potential.utils.PotentialTest;
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
                "0.0000, 0.0000, 0.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.15767805,
                5797,
                -117.49018161,
                6105,
                -28.28571436,
                6105,
                0.00000000,
                -6.42212098
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 0.0000, 0.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.45199032,
                5797,
                -184.64811367,
                6105,
                -15.15183374,
                6105,
                0.00000000,
                -54.83525644
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "0.0000, 1.0000, 0.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.58761295,
                5797,
                -168.76678648,
                6105,
                -18.69501522,
                6105,
                0.00000000,
                -61.51474289
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 1.0000, 0.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.88178733,
                5797,
                -184.27924195,
                6105,
                -8.39983635,
                6105,
                0.00000000,
                -109.92787835
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "0.0000, 0.0000, 1.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.09124761,
                5797,
                -169.62270694,
                6105,
                -32.83203586,
                6105,
                0.00000000,
                38.41082483
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 0.0000, 1.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.38555964,
                5797,
                -213.78774728,
                6105,
                -18.75451109,
                6105,
                0.00000000,
                -10.00231062
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "0.0000, 1.0000, 1.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.52117997,
                5797,
                -183.09872963,
                6105,
                -25.32650991,
                6105,
                0.00000000,
                -16.68179707
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 1.0000, 1.0000, 0.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.81535412,
                5797,
                -175.61829338,
                6105,
                -14.14018274,
                6105,
                0.00000000,
                -65.09493253
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "0.0000, 0.0000, 0.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.42088059,
                5797,
                -136.92341075,
                6105,
                -36.13604092,
                6105,
                0.00000000,
                34.82842021
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 0.0000, 0.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.71519284,
                5797,
                -188.40693546,
                6105,
                -21.93473964,
                6105,
                0.00000000,
                -13.58471525
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "0.0000, 1.0000, 0.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.85081533,
                5797,
                -166.31281253,
                6105,
                -25.85245332,
                6105,
                0.00000000,
                -20.26420170
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 1.0000, 0.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.14498969,
                5797,
                -166.15086065,
                6105,
                -14.50012773,
                6105,
                0.00000000,
                -68.67733716
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "0.0000, 0.0000, 1.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.35427749,
                5797,
                -138.30069801,
                6105,
                -42.20880342,
                6105,
                0.00000000,
                79.66136602
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 0.0000, 1.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.64858951,
                5797,
                -166.79133100,
                6105,
                -27.00713890,
                6105,
                0.00000000,
                31.24823057
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "0.0000, 1.0000, 1.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                85.78420970,
                5797,
                -129.88951761,
                6105,
                -33.89597126,
                6105,
                0.00000000,
                24.56874411
            },
            {
                "DEHK peptide",
                "src/main/java/ffx/potential/structures/DEHK.pdb",
                "1.0000, 1.0000, 1.0000, 1.0000",
                17.05831258,
                111,
                16.60445139,
                197,
                -1.23617487,
                170,
                0.0,
                0,
                0.00017436,
                69,
                24.76559922,
                278,
                0.0,
                0,
                0.00000716,
                13,
                0.0,
                0,
                86.07838383,
                5797,
                -106.73467401,
                6105,
                -21.59591450,
                6105,
                0.00000000,
                -23.84439134
            }

        });
  }

  @Test
  public void testEndStateEnergyAndGradient() {
    // Configure input arguments for the PhGradient script.
    logger.info(" Testing endstate energy for " + info + "at lambdas: " + key);

    // Choose a random Atom to test XYZ gradient.
    Random random = new Random();
    String randomAtom = Integer.toString(random.nextInt(111));
    String[] args = {"-v", "-d", "0.000001", "--ga", randomAtom, filename};
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
    assertEquals(info + " ExtendedSystemBias", extendedSystemBias, energyAndInteractionList[24],
        tolerance);
    // Total Energy
    assertEquals(info + " Total Energy", totalEnergy, energyAndInteractionList[25], tolerance);

    // Check for a gradient failure.
    assertEquals("DEHK gradient failures: ", 0, pHGradient.nFailures);
    assertEquals("DEHK ESV gradient failures: ", 0, pHGradient.nESVFailures);
  }
}

