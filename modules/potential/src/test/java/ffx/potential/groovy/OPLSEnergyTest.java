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

import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.random;
import static org.junit.Assert.assertEquals;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.groovy.Energy;
import ffx.potential.groovy.test.Gradient;
import ffx.potential.groovy.test.LambdaGradient;
import ffx.potential.utils.PotentialTest;
import java.util.Arrays;
import java.util.Collection;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

/** Test OPLS-AA energy and gradient. */
@RunWith(Parameterized.class)
public class OPLSEnergyTest extends PotentialTest {

  private final String info;
  private final String filename;
  private final int nAtoms;
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

  public OPLSEnergyTest(
      String info,
      String filename,
      int nAtoms,
      double bondEnergy,
      int nBonds,
      double angleEnergy,
      int nAngles,
      double torsionEnergy,
      int nTorsions,
      double improperTorsionEnergy,
      int nImproperTorsions,
      double vanDerWaalsEnergy,
      int nVanDerWaals,
      double permanentEnergy,
      int nPermanent) {

    this.info = info;
    this.filename = filename;
    this.nAtoms = nAtoms;
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
  }

  @Parameters
  public static Collection<Object[]> data() {
    return Arrays.asList(
        new Object[][] {
            {
                "OPLS Acetanilide Benchmark",
                "ffx/potential/structures/acetanilide-oplsaa.xyz",
                19,
                0.45009773,
                19,
                1.99659439,
                30,
                3.12713717,
                38,
                0.03415131,
                12,
                -10.22874283,
                7293,
                -24.11087720,
                2229
            },
            {
                "OPLS Ethylparaben Benchmark",
                "ffx/potential/structures/ethylparaben-oplsaa.xyz",
                44,
                1.23483793,
                44,
                1.17950432,
                70,
                0.58946840,
                88,
                0.01347012,
                22,
                -23.22505628,
                16759,
                -28.45972542,
                5029
            },
            {
                "OPLS Methylparaben Benchmark",
                "ffx/potential/structures/methylparaben-oplsaa.xyz",
                19,
                0.48582171,
                19,
                0.47331501,
                29,
                0.53855046,
                35,
                0.02598434,
                11,
                -8.71533364,
                7647,
                -16.384524645738537,
                2302
            },
            {
                "OPLS Paracetamol Benchmark",
                "ffx/potential/structures/paracetamol-oplsaa.xyz",
                20,
                0.63722563,
                20,
                2.63545869,
                31,
                3.01547224,
                40,
                0.06569712,
                10,
                -9.45895598,
                7831,
                -41.34298872,
                3334
            },
            {
                "OPLS Phenacetin Benchmark",
                "ffx/potential/structures/phenacetin-oplsaa.xyz",
                26,
                0.53495810,
                26,
                3.68523677,
                43,
                4.12682233,
                52,
                0.03932507,
                10,
                -14.31979877,
                10338,
                -31.561832525438113,
                3116
            }
        });
  }

  @Test
  public void testEnergy() {
    logger.info(" Testing energy for " + info);

    // Set-up the input arguments for the Energy script.
    String[] args = {"src/main/java/" + filename};
    binding.setVariable("args", args);

    // Evaluate the script.
    Energy energy = new Energy(binding).run();
    potentialScript = energy;

    ForceFieldEnergy forceFieldEnergy = energy.forceFieldEnergy;

    // Bond Energy
    assertEquals(info + " Bond Energy", bondEnergy, forceFieldEnergy.getBondEnergy(), tolerance);
    assertEquals(info + " Bond Count", nBonds, forceFieldEnergy.getNumberofBonds());
    // Angle Energy
    assertEquals(info + " Angle Energy", angleEnergy, forceFieldEnergy.getAngleEnergy(), tolerance);
    assertEquals(info + " Angle Count", nAngles, forceFieldEnergy.getNumberofAngles());
    // Torsional Angle
    assertEquals(
        info + " Torsion Energy", torsionEnergy, forceFieldEnergy.getTorsionEnergy(), tolerance);
    assertEquals(info + " Torsion Count", nTorsions, forceFieldEnergy.getNumberofTorsions());
    // Improper Torsional Angle
    assertEquals(
        info + " Improper Torsion Energy",
        improperTorsionEnergy,
        forceFieldEnergy.getImproperTorsionEnergy(),
        tolerance);
    assertEquals(
        info + " Improper Torsion Count",
        nImproperTorsions,
        forceFieldEnergy.getNumberofImproperTorsions());
    // van Der Waals
    assertEquals(
        info + " van Der Waals Energy",
        vanDerWaalsEnergy,
        forceFieldEnergy.getVanDerWaalsEnergy(),
        tolerance);
    assertEquals(
        info + " van Der Waals Count", nVanDerWaals, forceFieldEnergy.getVanDerWaalsInteractions());
    // Permanent Multipoles
    assertEquals(
        info + " Permanent Multipole Energy",
        permanentEnergy,
        forceFieldEnergy.getPermanentMultipoleEnergy(),
        tolerance);
    assertEquals(
        info + " Permanent Multipole Count",
        nPermanent,
        forceFieldEnergy.getPermanentInteractions());
  }

  @Test
  public void testGradient() {
    logger.info(" Testing Cartesian gradient(s) for " + info);

    // Set-up the input arguments for the Gradient script.
    // Choose a random atom to test.
    int atomID = (int) floor(random() * nAtoms) + 1;
    double stepSize = 1.0e-5;
    String[] args = {"--ga", Integer.toString(atomID),
        "--dx", Double.toString(stepSize),
        "src/main/java/" + filename
    };
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    Gradient gradient = new Gradient(binding).run();
    potentialScript = gradient;

    assertEquals(info + " gradient failures: ", 0, gradient.nFailures);
  }

  @Test
  public void testLambdaGradient() {
    logger.info(" Testing lambda gradient(s) for " + info);

    // Set-up the input arguments for the LambdaGradient script.
    // Choose a random atom to test dEdX gradient.
    int atomID = (int) floor(random() * nAtoms) + 1;
    double stepSize = 1.0e-5;
    String[] args = {
        "--ga", Integer.toString(atomID),
        "--dx", Double.toString(stepSize),
        "--tol", Double.toString(tolerance),
        "--ac", "1" + "-" + nAtoms,
        "-l", "0.5",
        "src/main/java/" + filename
    };
    binding.setVariable("args", args);

    // Construct and evaluate the script.
    LambdaGradient lambdaGradient = new LambdaGradient(binding).run();
    potentialScript = lambdaGradient;

    assertEquals(info + " dEdL failures: ", 0, lambdaGradient.ndEdLFailures);
    assertEquals(info + " d2EdL2 failures: ", 0, lambdaGradient.nd2EdL2Failures);
    assertEquals(info + " dEdXdL failures: ", 0, lambdaGradient.ndEdXdLFailures);
    assertEquals(info + " dEdX failures: ", 0, lambdaGradient.ndEdXFailures);
  }
}
