//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
//******************************************************************************
package ffx.algorithms.commands.test;

import ffx.algorithms.cli.AlgorithmsCommand;
import ffx.utilities.FFXBinding;
import org.apache.commons.math3.special.Erf;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;
import picocli.CommandLine.Parameters;

import java.util.List;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

/**
 * The Freefix script calculates analytical free energy, entropy, and enthalpy corrections
 * for flat-bottom harmonic restraints with inner and outer radii.
 * <br>
 * Usage:
 * <br>
 * ffxc test.Freefix [options] &lt;filename&gt;
 */
@Command(description = " Calculate restraint free energy corrections.", name = "test.Freefix")
public class Freefix extends AlgorithmsCommand {

  /**
   * --ri Constant Inner Radius value.
   */
  @Option(names = {"--ri"}, paramLabel = "0.00", defaultValue = "0.00",
      description = "Constant Inner Radius value.")
  private double innerRadius;

  /**
   * --fi Constant Inner Force value.
   */
  @Option(names = {"--fi"}, paramLabel = "0.0", defaultValue = "0.0",
      description = "Constant Inner Force value.")
  private double innerForce;

  /**
   * --ro Constant Outer Radius value.
   */
  @Option(names = {"--ro"}, paramLabel = "0.00", defaultValue = "0.00",
      description = "Constant Outer Radius value.")
  private double outerRadius;

  /**
   * --fo Constant Outer Force value.
   */
  @Option(names = {"--fo"}, paramLabel = "15.0", defaultValue = "15.0",
      description = "Constant Outer Force value.")
  private double outerForce;

  /**
   * --temp Constant Temperature value.
   */
  @Option(names = {"--temp"}, paramLabel = "298.0", defaultValue = "298.0",
      description = "Constant Temperature value.")
  private double temperature;

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = "Atomic coordinate files in PDB or XYZ format.")
  private List<String> filenames;

  public static final double AVOGADRO = 6.02214076e23;
  public static final double STD_CONVERSION = 1.0e27 / AVOGADRO;

  /**
   * Freefix Constructor.
   */
  public Freefix() {
    super();
  }

  /**
   * Freefix Constructor.
   * @param binding The Binding to use.
   */
  public Freefix(FFXBinding binding) {
    super(binding);
  }

  /**
   * Freefix constructor that sets the command line arguments.
   * @param args Command line arguments.
   */
  public Freefix(String[] args) {
    super(args);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Freefix run() {

    // Init the context and bind variables.
    if (!init()) {
      return this;
    }

    if (innerForce == 0.0) {
      innerForce = 1.0;
      innerRadius = 0.0;
    }
    if (outerForce == 0.0) {
      outerForce = 1.0;
    }

    System.out.printf("%-35s %.4f Ang%n", "Inner Flat-Bottom Radius:", innerRadius);
    System.out.printf("%-35s %.4f Kcal/mole/Ang^2%n", "Inner Force Constant:", innerForce);
    System.out.printf("%-35s %.4f Ang%n", "Outer Flat-Bottom Radius:", outerRadius);
    System.out.printf("%-35s %.4f Kcal/mole/Ang^2%n", "Outer Force Constant:", outerForce);
    System.out.printf("%-35s %.4f Kelvin%n%n", "System Temperature Value:", temperature);

    double kB = 0.0019872041;
    double kt = kB * temperature;

    double[] volumeResults = computeVolumeIntegrals(innerRadius, innerForce, outerRadius, outerForce, temperature);
    double vol = volumeResults[0];
    double dvol = volumeResults[1];

    logger.info(
        format("%-35s %.4f Ang^3%n", "Analytical Volume Integral:", vol));
    logger.info(
        format("%-35s %.4f Ang^3/K%n", "Analytical dVol/dT Value:", dvol));

    double dg = -kt * Math.log(vol / STD_CONVERSION);
    double ds = -dg / temperature + kt * dvol / vol;
    double dh = dg + temperature * ds;

    logger.info(
        format("%-35s %.4f Kcal/mole%n", "Restraint Free Energy:", dg));
    logger.info(
        format("%-35s %.4f Kcal/mole/K%n", "Restraint Entropy Value:", ds));
    logger.info(
        format("%-35s %.4f Kcal/mole%n", "Restraint Enthalpy Value:", dh));
    logger.info(
        format("%-35s %.4f Kcal/mole%n", "Restraint -T deltaS Value:", -temperature * ds));

    return this;
  }

  /**
   * Compute volume integrals and their temperature derivatives for restraint thermodynamics.
   *
   * @param innerRadius Inner flat-bottom radius
   * @param innerForce Inner force constant
   * @param outerRadius Outer flat-bottom radius
   * @param outerForce Outer force constant
   * @param temperature System temperature
   * @return Array containing [volume, dvolume/dT]
   */
  public static double[] computeVolumeIntegrals(double innerRadius, double innerForce,
      double outerRadius, double outerForce, double temperature) {
    double kB = 0.0019872041;
    double kt = kB * temperature;
    double volume1 = 0.0, volume2, volume3 = 0.0;

    if (innerRadius != 0.0) {
      double term1_v1 = 2.0 * PI * innerRadius * (-2.0 + exp(-pow(innerRadius, 2) * innerForce / kt)) * kt / innerForce;
      double term2_v1 = sqrt(kt * pow(PI / innerForce, 3)) *
          (2.0 * innerForce * pow(innerRadius, 2) + kt) * Erf.erf(innerRadius * sqrt(innerForce / kt));
      volume1 = term1_v1 + term2_v1;
    }

    volume2 = (4.0 / 3.0) * PI * (pow(outerRadius, 3) - pow(innerRadius, 3));

    if (outerRadius == 0.0) {
      double term1_v3 = sqrt(kt * pow(PI / outerForce, 3)) * kt;
      volume3 = term1_v3;
    } else {
      double term1_v3 = sqrt(kt * pow(PI / outerForce, 3)) *
          (2.0 * outerForce * pow(outerRadius, 2) + kt + 4.0 * outerRadius * sqrt(kt * outerForce / PI));
      volume3 = term1_v3;
    }

    double volume = volume1 + volume2 + volume3;

    double dv1 = (innerRadius != 0.0) ?
        computeDv1(innerRadius, innerForce, kt, temperature) : 0.0;
    double dv3 = (outerRadius != 0.0) ?
        computeDv3(outerRadius, outerForce, kt, temperature) :
        1.5 * pow(kt * PI / outerForce, 1.5) / temperature;

    double dvolume = dv1 + dv3;
    return new double[]{volume, dvolume};
  }

  /**
   * Compute temperature derivative of inner volume integral.
   *
   * @param innerRadius Inner flat-bottom radius
   * @param innerForce Inner force constant
   * @param kt Boltzmann constant times temperature
   * @param temperature System temperature
   * @return Temperature derivative of inner volume
   */
  private static double computeDv1(double innerRadius, double innerForce, double kt, double temperature) {
    double term1_dv1 = 2.0 * PI * pow(innerRadius, 3) * exp(-pow(innerRadius, 2) * innerForce / kt) / temperature;
    double term2_dv1 = 2.0 * PI * innerRadius * (-2.0 + exp(-pow(innerRadius, 2) * innerForce / kt)) * kt / (innerForce * temperature);
    double term3_dv1 = 0.5 * sqrt(pow(PI / innerForce, 3)) * sqrt(kt) *
        (2.0 * innerForce * pow(innerRadius, 2) + kt) * Erf.erf(innerRadius * sqrt(innerForce / kt)) / temperature;
    double term4_dv1 = -PI * innerRadius * exp(-pow(innerRadius, 2) * innerForce / kt) *
        (2.0 * innerForce * pow(innerRadius, 2) + kt) / (innerForce * temperature);
    double term5_dv1 = sqrt(pow(kt * PI / innerForce, 3)) * Erf.erf(innerRadius * sqrt(innerForce / kt)) / temperature;
    return term1_dv1 + term2_dv1 + term3_dv1 + term4_dv1 + term5_dv1;
  }

  /**
   * Compute temperature derivative of outer volume integral.
   *
   * @param outerRadius Outer flat-bottom radius
   * @param outerForce Outer force constant
   * @param kt Boltzmann constant times temperature
   * @param temperature System temperature
   * @return Temperature derivative of outer volume
   */
  private static double computeDv3(double outerRadius, double outerForce, double kt, double temperature) {
    double term1_dv3 = sqrt(kt * pow(PI / outerForce, 3)) * outerForce * pow(outerRadius, 2) / temperature;
    double term2_dv3 = 4.0 * kt * (PI / outerForce) * outerRadius / temperature;
    double term3_dv3 = 1.5 * sqrt(pow(kt * PI / outerForce, 3)) / temperature;
    return term1_dv3 + term2_dv3 + term3_dv3;
  }
}
