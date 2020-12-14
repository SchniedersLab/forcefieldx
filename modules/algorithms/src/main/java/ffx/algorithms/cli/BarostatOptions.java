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
package ffx.algorithms.cli;

import static java.lang.String.format;

import ffx.algorithms.dynamics.Barostat;
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;
import java.util.logging.Logger;
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that use a barostat/NPT.
 *
 * @author Michael J. Schnieders
 * @author Hernan V. Bernabe
 * @since 1.0
 */
public class BarostatOptions {

  /**
   * Default maximum density constraint on the barostat that prevents reduction in unit cell
   * (particularly at or near vapor states).
   */
  public static final String DEFAULT_MAX_DENSITY = "1.6";
  /**
   * Default "tin box" constraint on the barostat that prevents expansion of the unit cell
   * (particularly at or near vapor states).
   */
  public static final String DEFAULT_MIN_DENSITY = "0.75";
  /** Default size of proposed unit cell volume moves (uniformly distributed) in Angstroms^3. */
  public static final String DEFAULT_MAX_VOLUME_MOVE = "1.0";
  /** Default mean number of MD steps (Poisson distribution) between barostat move proposals. */
  public static final String DEFAULT_BAROSTAT_INTERVAL = "10";
  /**
   * Default number of Barostat moves between print statements.
   */
  public static final String DEFAULT_BAROSTAT_PRINT_INTERVAL = "1000";

  private static final Logger logger = Logger.getLogger(BarostatOptions.class.getName());

  /**
   * The ArgGroup keeps the BarostatOptions together when printing help.
   */
  @ArgGroup(heading = "%n Monte Carlo Pressure Options%n", validate = false)
  public BarostatOptionGroup group = new BarostatOptionGroup();

  /**
   * If pressure has been set &gt; 0, creates a Barostat around a CrystalPotential, else returns the
   * original, unmodified CrystalPotential.
   *
   * @param molecularAssembly Primary molecularAssembly of the CrystalPotential.
   * @param crystalPotential A CrystalPotential.
   * @return Either a Barostat (NPT enabled) or cpot.
   */
  public CrystalPotential checkNPT(
      MolecularAssembly molecularAssembly, CrystalPotential crystalPotential) {
    if (group.pressure > 0) {
      return createBarostat(molecularAssembly, crystalPotential);
    } else {
      return crystalPotential;
    }
  }

  /**
   * Creates a Barostat around a CrystalPotential.
   *
   * @param assembly Primary assembly of the CrystalPotential.
   * @param crystalPotential A CrystalPotential.
   * @return An NPT potential.
   * @throws IllegalArgumentException If this BarostatOptions has pressure &lt;= 0.
   */
  public Barostat createBarostat(MolecularAssembly assembly, CrystalPotential crystalPotential)
      throws IllegalArgumentException {
    if (group.pressure > 0) {
      Barostat barostat = new Barostat(assembly, crystalPotential);
      barostat.setPressure(group.pressure);
      barostat.setIsotropic(group.isotropic);
      barostat.setMaxDensity(group.maxD);
      barostat.setMinDensity(group.minD);
      barostat.setBarostatPrintFrequency(group.barPrint);
      double dens = barostat.density();
      if (dens < group.minD) {
        logger.info(
            format(
                " Barostat: initial density %9.4g < minimum density %9.4g, resetting to minimum density",
                dens, group.minD));
        barostat.setDensity(group.minD);
      } else if (dens > group.maxD) {
        logger.info(
            format(
                " Barostat: initial density %9.4g > maximum density %9.4g, resetting to maximum density",
                dens, group.maxD));
        barostat.setDensity(group.maxD);
      }
      barostat.setMaxVolumeMove(group.maxV);
      barostat.setMeanBarostatInterval(group.barInt);
      return barostat;
    } else {
      throw new IllegalArgumentException(" Pressure is <= 0; cannot create a Barostat!");
    }
  }

  /**
   * -p or --npt Specify use of a MC Barostat at the given pressure (default 0 = constant volume).
   *
   * @return Returns the pressure (atm).
   */
  public double getPressure() {
    return group.pressure;
  }

  public void setPressure(double pressure) {
    group.pressure = pressure;
  }

  /**
   * Restrict the MC Barostat to isotropic moves. The lattice angles are * held fixed, and lattice
   * lengths are scaled equally.
   *
   * @return Returns true if the Barostat is isotropic.
   */
  public boolean isIsotropic() {
    return group.isotropic;
  }

  public void setIsotropic(boolean isotropic) {
    group.isotropic = isotropic;
  }

  /**
   * The maximum density accepted by the MC Barostat (g/cc).
   *
   * @return Returns the maximum density.
   */
  public double getMaxD() {
    return group.maxD;
  }

  public void setMaxD(double maxD) {
    group.maxD = maxD;
  }

  /**
   * The minimum density accepted by the MC Barostat (g/cc).
   *
   * @return Returns the minimum density.
   */
  public double getMinD() {
    return group.minD;
  }

  public void setMinD(double minD) {
    group.minD = minD;
  }

  /**
   * The volume of a proposed unit cell volume moves (uniformly distributed) in Angstroms^3.
   *
   * @return Returns the maximum volume move.
   */
  public double getMaxV() {
    return group.maxV;
  }

  public void setMaxV(double maxV) {
    group.maxV = maxV;
  }

  /**
   * The mean number of MD steps (Poisson distribution) between barostat move proposals.
   *
   * @return Returns the mean number of MD steps between barostat MC trials.
   */
  public int getBarInt() {
    return group.barInt;
  }

  public void setBarInt(int barInt) {
    group.barInt = barInt;
  }

  /**
   * --bpi or --barostatPrintInterval Sets the number of Barostat MC cycles between print
   * statements.
   */
  public int getPrintInt() {
    return group.barPrint;
  }

  public void setPrintInt(int printInterval) {
    group.barPrint = printInterval;
  }

  /**
   * Collection of Barostat Options.
   */
  private static class BarostatOptionGroup {

    /**
     * -p or --npt Specify use of a MC Barostat at the given pressure (default 0 = constant volume).
     */
    @Option(
        names = {"-p", "--npt"},
        paramLabel = "0",
        defaultValue = "0",
        description =
            "Specify use of a MC Barostat at the given pressure; the default 0 disables NPT (atm).")
    private double pressure;
    /**
     * -iso or --isotropic Restrict the MC Barostat to isotropic moves. The lattice angles are held
     * fixed, and lattice lengths are scaled equally.
     */
    @Option(
        names = {"--iso", "--isotropic"},
        defaultValue = "false",
        description = "Restrict the MC Barostat to isotropic moves.")
    private boolean isotropic;
    /** --maxD or --maxDensity Specify the maximum density accepted by the MC Barostat (g/cc). */
    @Option(
        names = {"--maxD", " --maxDensity"},
        paramLabel = DEFAULT_MAX_DENSITY,
        defaultValue = DEFAULT_MAX_DENSITY,
        description = "Specify the maximum density accepted by the MC Barostat (g/cc).")
    private double maxD;
    /** --minD or --minDensity Specify the minimum density accepted by the MC Barostat (g/cc). */
    @Option(
        names = {"--minD", " --minDensity"},
        paramLabel = DEFAULT_MIN_DENSITY,
        defaultValue = DEFAULT_MIN_DENSITY,
        description = "Specify the minimum density accepted by the MC Barostat (g/cc).")
    private double minD;
    /**
     * --maxSM or --maxSideMove Sets the width of proposed unit cell side length moves (uniformly
     * distributed) in Angstroms.
     */
    @Option(
        names = {"--maxV", "--maxVolumeMove"},
        paramLabel = DEFAULT_MAX_VOLUME_MOVE,
        defaultValue = DEFAULT_MAX_VOLUME_MOVE,
        description =
            "Default width of proposed unit cell side length moves (uniformly distributed) in Angstroms.")
    private double maxV;
    /**
     * --barInt or --meanBarostatInterval Sets the mean number of MD steps (Poisson distribution)
     * between barostat move proposals.
     */
    @Option(
        names = {"--barInt", "--meanBarostatInterval"},
        paramLabel = DEFAULT_BAROSTAT_INTERVAL,
        defaultValue = DEFAULT_BAROSTAT_INTERVAL,
        description = "Sets the mean number of MD steps between barostat move proposals.")
    private int barInt;
    /**
     * --bpi or --barostatPrintInterval Sets the number of Barostat MC cycles between print
     * statements.
     */
    @Option(
        names = {"--bpi", "--barostatPrintInterval"},
        paramLabel = DEFAULT_BAROSTAT_PRINT_INTERVAL,
        defaultValue = DEFAULT_BAROSTAT_PRINT_INTERVAL,
        description = "Sets the number of Barostat MC cycles between print statements.")
    private int barPrint;
  }
}
