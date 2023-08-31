// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
import ffx.potential.groovy.test.LambdaGradient;
import ffx.potential.utils.PotentialTest;
import org.junit.Test;

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.random;
import static org.junit.Assert.assertEquals;

/**
 * Test the Energy script.
 */
public class ParentEnergyTest extends PotentialTest {

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
  private final int nStretchTorsions;
  private final int nAngleTorsions;
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
  private final double stretchTorsionEnergy;
  private final double angleTorsionEnergy;
  private final double vanDerWaalsEnergy;
  private final double permanentEnergy;
  private final double polarizationEnergy;
  private final double gkEnergy;
  private final double totalEnergy;
  private final boolean testOpenMM;
  private final double tolerance = 1.0e-2;

  public ParentEnergyTest(String info, String filename, int nAtoms, double bondEnergy, int nBonds,
                          double angleEnergy, int nAngles, double stretchBendEnergy, int nStretchBends,
                          double ureyBradleyEnergy, int nUreyBradleys, double outOfPlaneBendEnergy,
                          int nOutOfPlaneBends, double torsionEnergy, int nTorsions, double improperTorsionEnergy,
                          int nImproperTorsions, double piOrbitalTorsionEnergy, int nPiOrbitalTorsions,
                          double torsionTorsionEnergy, int nTorsionTorsions, double stretchTorsionEnergy,
                          int nStretchTorsions, double angleTorsionEnergy, int nAngleTorsions, double vanDerWaalsEnergy,
                          int nVanDerWaals, double permanentEnergy, int nPermanent, double polarizationEnergy,
                          int nPolar, double gkEnergy, int nGK, boolean testOpenMM) {
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
    this.stretchTorsionEnergy = stretchTorsionEnergy;
    this.nStretchTorsions = nStretchTorsions;
    this.angleTorsionEnergy = angleTorsionEnergy;
    this.nAngleTorsions = nAngleTorsions;
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
        + torsionEnergy + improperTorsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy + stretchTorsionEnergy
        + angleTorsionEnergy + vanDerWaalsEnergy + permanentEnergy + polarizationEnergy + gkEnergy;
  }

  @Test
  public void testEnergy() {
    if (nAtoms > 10000 && !ffxCI) {
      return;
    }
    logger.info(" Testing energy for " + info);

    String[] args = {"src/main/java/" + filename};
    binding.setVariable("args", args);

    // Create and evaluate the script.
    Energy energy = new Energy(binding).run();
    potentialScript = energy;
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
    // Stretch-Torsion
    assertEquals(info + " Stretch-Torsion Energy", stretchTorsionEnergy, forceFieldEnergy.getStretchTorsionEnergy(), tolerance);
    assertEquals(info + " Stretch-Torsion Count", nStretchTorsions, forceFieldEnergy.getNumberofStretchTorsions());
    // Angle-Torsion
    assertEquals(info + " Angle-Torsion Energy", angleTorsionEnergy, forceFieldEnergy.getAngleTorsionEnergy(), tolerance);
    assertEquals(info + " Angle-Torsion Count", nAngleTorsions, forceFieldEnergy.getNumberofAngleTorsions());
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
  }

  @Test
  public void testGradient() {
    if (!ffxCI) {
      if (nAtoms > 5000) {
        return;
      } else if (nGK > 0 && nAtoms > 500) {
        return;
      }
    }


    logger.info(" Testing Cartesian gradient(s) for " + info);

    // Set up the input arguments for the Gradient script.
    // Choose a random atom to test.
    int atomID = (int) floor(random() * nAtoms) + 1;
    double stepSize = 1.0e-5;
    String[] args = {"--ga", Integer.toString(atomID),
        "--dx", Double.toString(stepSize),
        "--tol", Double.toString(tolerance),
        "src/main/java/" + filename};
    binding.setVariable("args", args);

    // Create and evaluate the script.
    Gradient gradient = new Gradient(binding).run();
    potentialScript = gradient;
    assertEquals(info + " gradient failures: ", 0, gradient.nFailures);
  }

  @Test
  public void testLambdaGradient() {
    if (!ffxCI) {
      if (nAtoms > 5000) {
        return;
      } else if (nGK > 0 && nAtoms > 500) {
        return;
      }
    }

    logger.info(" Testing lambda gradient(s) for " + info);

    // Set up the input arguments for the Lambda Gradient script.
    // Choose a random atom to test dEdX gradient.
    int atomID = (int) floor(random() * nAtoms) + 1;
    double stepSize = 1.0e-5;
    String[] args = {"--ga", Integer.toString(atomID),
        "--dx", Double.toString(stepSize),
        "--tol", Double.toString(tolerance),
        "--ac", "ALL",
        "-l", "0.9", "src/main/java/" + filename};
    binding.setVariable("args", args);

    // Create and evaluate the script.
    LambdaGradient lambdaGradient = new LambdaGradient(binding).run();
    potentialScript = lambdaGradient;
    assertEquals(info + " dEdL failures: ", 0, lambdaGradient.ndEdLFailures);
    assertEquals(info + " d2EdL2 failures: ", 0, lambdaGradient.nd2EdL2Failures);
    assertEquals(info + " dEdXdL failures: ", 0, lambdaGradient.ndEdXdLFailures);
    assertEquals(info + " dEdX failures: ", 0, lambdaGradient.ndEdXFailures);
  }

  @Test
  public void testOpenMMEnergy() {
    if (!testOpenMM || !ffxOpenMM) {
      return;
    }
    logger.info(" Testing OpenMM energy for " + info);

    // Set up the input arguments for the Energy script.
    String[] args = {"src/main/java/" + filename};
    binding.setVariable("args", args);
    System.setProperty("platform", "OMM");

    // Construct and evaluate the script.
    Energy energy = new Energy(binding).run();
    potentialScript = energy;
    double openMMTolerance = 0.5;
    assertEquals(info + " OpenMM Energy", totalEnergy, energy.energy, openMMTolerance);
  }
}
