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
package ffx.potential.nonbonded.implicit;

import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.switching.MultiplicativeSwitch;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.parameters.ForceField;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The ChandlerCavitation class smoothly switches between a volume based dependence for small
 * solutes to a surface area dependence for large solutes.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ChandlerCavitation {

  private static final Logger logger = Logger.getLogger(ChandlerCavitation.class.getName());
  /** Volume / Surface Area are switched on / off over 7.0 Angstroms; half is 3.5 A. */
  private final double HALF_SWITCH_RANGE = 3.5;
  /**
   * A value of 0.2 A smoothly moves from Volume dependence to Surface Area dependence.
   *
   * <p>However, a smaller offset of 0.0 to ~0.1 A overshoots the limiting surface tension.
   */
  private final double SA_SWITCH_OFFSET = 0.2;

  private final boolean doVolume;
  private final boolean doSurfaceArea;

  private final GaussVol gaussVol;
  private final ConnollyRegion connollyRegion;
  /** Surface area (Ang^2). */
  private double surfaceArea;
  /** Surface area energy (kcal/mol). */
  private double surfaceAreaEnergy;
  /** Volume (Ang^3). */
  private double volume;
  /** Volume energy (kcal/mol). */
  private double volumeEnergy;
  /** Cavitation energy, which is a function of volume and/or surface area. */
  private double cavitationEnergy;
  /**
   * Effective radius probe.
   *
   * <p>In cavitation volume scaling regime, approximate solvent excluded volume and effective
   * radius are computed as follow.
   *
   * <p>1) GaussVol vdW volume is computed from defined radii. 2) Effective radius is computed as
   * Reff = cbrt(3.0 * volume / (4.0 * PI)) + effectiveRadiusProbe. 3) Solvent Excluded Volume = 4/3
   * * Pi * Reff^3
   */
  private double effectiveRadius;
  /**
   * Solvent pressure in kcal/mol/Ang^3. Original value from Schnieders thesis work: 0.0327 Value
   * based on testing with Schnieders thesis test set, Sept 2019: 0.11337
   */
  private double solventPressure = GeneralizedKirkwood.DEFAULT_SOLVENT_PRESSURE;
  /** Surface tension in kcal/mol/Ang^2. */
  private double surfaceTension = GeneralizedKirkwood.DEFAULT_CAVDISP_SURFACE_TENSION;
  /**
   * Radius where volume dependence crosses over to surface area dependence (approximately at 1 nm).
   */
  private double crossOver = GeneralizedKirkwood.DEFAULT_CROSSOVER;
  /** Begin turning off the Volume term. */
  private double beginVolumeOff = crossOver - HALF_SWITCH_RANGE;
  /** Volume term is zero at the cut-off. */
  private double endVolumeOff = beginVolumeOff + 2.0 * HALF_SWITCH_RANGE;
  /** Begin turning off the SA term. */
  private double beginSurfaceAreaOff = crossOver + SA_SWITCH_OFFSET + HALF_SWITCH_RANGE;
  /** SA term is zero at the cut-off. */
  private double endSurfaceAreaOff = beginSurfaceAreaOff - 2.0 * HALF_SWITCH_RANGE;
  /** Volume multiplicative switch. */
  private MultiplicativeSwitch volumeSwitch =
      new MultiplicativeSwitch(beginVolumeOff, endVolumeOff);
  /** Surface area multiplicative switch. */
  private MultiplicativeSwitch surfaceAreaSwitch =
      new MultiplicativeSwitch(beginSurfaceAreaOff, endSurfaceAreaOff);

  private final int nAtoms;
  private final Atom[] atoms;

  public ChandlerCavitation(Atom[] atoms, GaussVol gaussVol, ForceField forceField) {
    this.atoms = atoms;
    this.nAtoms = atoms.length;
    this.gaussVol = gaussVol;
    this.connollyRegion = null;

    doVolume = forceField.getBoolean("VOLUMETERM", true);
    doSurfaceArea = forceField.getBoolean("AREATERM", true);
  }

  public ChandlerCavitation(Atom[] atoms, ConnollyRegion connollyRegion, ForceField forceField) {
    this.atoms = atoms;
    this.nAtoms = atoms.length;
    this.connollyRegion = connollyRegion;
    this.gaussVol = null;

    doVolume = forceField.getBoolean("VOLUMETERM", true);
    doSurfaceArea = forceField.getBoolean("AREATERM", true);
  }

  /**
   * Compute molecular volume and surface area.
   *
   * @param positions Atomic positions to use.
   * @param gradient Atomic coordinate gradient.
   * @return The cavitation energy.
   */
  public double energyAndGradient(double[][] positions, AtomicDoubleArray3D gradient) {
    if (gaussVol != null) {
      return energyAndGradientGausVol(positions, gradient);
    } else {
      return energyAndGradientConnolly(gradient);
    }
  }

  /**
   * Compute the cavitation energy.
   *
   * @param gradient Add the gradient to this AtomicDoubleArray3D.
   * @return Returns the cavitation energy.
   */
  public double energyAndGradientConnolly(AtomicDoubleArray3D gradient) {

    connollyRegion.init(atoms, true);
    connollyRegion.runVolume();

    // Calculate a purely surface area based cavitation energy.
    surfaceArea = connollyRegion.getSurfaceArea();
    surfaceAreaEnergy = surfaceArea * surfaceTension;

    // Calculate a purely volume based cavitation energy.
    volume = connollyRegion.getVolume();
    volumeEnergy = volume * solventPressure;

    // effectiveRadius = 0.5 * sqrt(surfaceArea / PI);
    // double reff = effectiveRadius;
    // double reff2 = reff * reff;
    // double reff3 = reff2 * reff;
    // double reff4 = reff3 * reff;
    // double reff5 = reff4 * reff;
    // double dReffdvdW = reff / (2.0 * surfaceArea);

    // Use Volume to find an effective cavity radius.
    effectiveRadius = cbrt(3.0 * volume / (4.0 * PI));
    double reff = effectiveRadius;
    double reff2 = reff * reff;
    double reff3 = reff2 * reff;
    double reff4 = reff3 * reff;
    double reff5 = reff4 * reff;
    double vdWVolPI23 = pow(volume / PI, 2.0 / 3.0);
    double dReffdvdW = 1.0 / (pow(6.0, 2.0 / 3.0) * PI * vdWVolPI23);

    // Find the cavitation energy using a combination of volume and surface area dependence.
    if (doVolume) {
      if (!doSurfaceArea || reff < beginVolumeOff) {
        // Find cavity energy from only the molecular volume.
        cavitationEnergy = volumeEnergy;
        double[][] volumeGradient = connollyRegion.getVolumeGradient();
        for (int i = 0; i < nAtoms; i++) {
          double dx = solventPressure * volumeGradient[0][i];
          double dy = solventPressure * volumeGradient[1][i];
          double dz = solventPressure * volumeGradient[2][i];
          gradient.add(0, i, dx, dy, dz);
        }
      } else if (reff <= endVolumeOff) {
        // Include a tapered molecular volume.
        double taper = volumeSwitch.taper(reff, reff2, reff3, reff4, reff5);
        cavitationEnergy = taper * volumeEnergy;
        double dtaper = volumeSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
        double factor = dtaper * volumeEnergy + taper * solventPressure;
        double[][] volumeGradient = connollyRegion.getVolumeGradient();
        for (int i = 0; i < nAtoms; i++) {
          double dx = factor * volumeGradient[0][i];
          double dy = factor * volumeGradient[1][i];
          double dz = factor * volumeGradient[2][i];
          gradient.add(0, i, dx, dy, dz);
        }
      }
    }

    // TODO: We need to include the surface area contribution to the gradient.
    if (doSurfaceArea) {
      if (!doVolume || reff > beginSurfaceAreaOff) {
        // Find cavity energy from only SA.
        cavitationEnergy = surfaceAreaEnergy;
        // addSurfaceAreaGradient(0.0, surfaceTension, gradient);
      } else if (reff >= endSurfaceAreaOff) {
        // Include a tapered surface area term.
        double taperSA = surfaceAreaSwitch.taper(reff, reff2, reff3, reff4, reff5);
        cavitationEnergy += taperSA * surfaceAreaEnergy;
        // double dtaperSA = surfaceAreaSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
        // addSurfaceAreaGradient(dtaperSA * surfaceAreaEnergy, taperSA * surfaceTension, gradient);
      }
    }

    if (logger.isLoggable(Level.FINE)) {
      logger.fine("\n Connolly");
      logger.fine(format(" Volume:              %8.3f (Ang^3)", volume));
      logger.fine(format(" Volume Energy:       %8.3f (kcal/mol)", volumeEnergy));
      logger.fine(format(" Surface Area:        %8.3f (Ang^2)", surfaceArea));
      logger.fine(format(" Surface Area Energy: %8.3f (kcal/mol)", surfaceAreaEnergy));
      logger.fine(format(" Volume + SA Energy:  %8.3f (kcal/mol)", cavitationEnergy));
      logger.fine(format(" Effective Radius:    %8.3f (Ang)", reff));
    }

    return cavitationEnergy;
  }

  /**
   * Compute molecular volume and surface area.
   *
   * @param positions Atomic positions to use.
   * @param gradient Atomic coordinate gradient.
   * @return The cavitation energy.
   */
  public double energyAndGradientGausVol(double[][] positions, AtomicDoubleArray3D gradient) {

    gaussVol.computeVolumeAndSA(positions);

    // Save the total molecular volume.
    volume = gaussVol.getVolume();
    double[] volumeGradient = gaussVol.getVolumeGradient();

    // Calculate a purely volume based cavitation energy.
    volumeEnergy = volume * solventPressure;

    // Calculate the surface area.
    surfaceArea = gaussVol.getSurfaceArea();
    double[] surfaceAreaGradient = gaussVol.getSurfaceAreaGradient();

    // Calculate a purely surface area based cavitation energy.
    surfaceAreaEnergy = surfaceArea * surfaceTension;

    // Use SA to find an effective cavity radius.
    effectiveRadius = 0.5 * sqrt(surfaceArea / PI);
    double reff = effectiveRadius;
    double reff2 = reff * reff;
    double reff3 = reff2 * reff;
    double reff4 = reff3 * reff;
    double reff5 = reff4 * reff;
    double dRdSA = reff / (2.0 * surfaceArea);

    // Use Volume to find an effective cavity radius.
    //        effectiveRadius = cbrt(3.0 * volume / (4.0 * PI));
    //        double reff = effectiveRadius;
    //        double reff2 = reff * reff;
    //        double reff3 = reff2 * reff;
    //        double reff4 = reff3 * reff;
    //        double reff5 = reff4 * reff;
    //        double volPI23 = pow(volume / PI, 2.0 / 3.0);
    //        double dReffdvdW = 1.0 / (pow(6.0, 2.0 / 3.0) * PI * volPI23);

    // Find the cavitation energy using a combination of volume and surface area dependence.
    if (doVolume) {
      if (!doSurfaceArea || reff < beginVolumeOff) {
        // Find cavity energy from only the molecular volume.
        cavitationEnergy = volumeEnergy;
        addVolumeGradient(0.0, solventPressure, surfaceAreaGradient, volumeGradient, gradient);
      } else if (reff <= endVolumeOff) {
        // Include a tapered molecular volume.
        double taper = volumeSwitch.taper(reff, reff2, reff3, reff4, reff5);
        double dtaper = volumeSwitch.dtaper(reff, reff2, reff3, reff4) * dRdSA;
        cavitationEnergy = taper * volumeEnergy;
        double factorSA = dtaper * volumeEnergy;
        double factorVol = taper * solventPressure;
        addVolumeGradient(factorSA, factorVol, surfaceAreaGradient, volumeGradient, gradient);
      }
    }

    if (doSurfaceArea) {
      if (!doVolume || reff > beginSurfaceAreaOff) {
        // Find cavity energy from only SA.
        cavitationEnergy = surfaceAreaEnergy;
        addSurfaceAreaGradient(surfaceTension, surfaceAreaGradient, gradient);
      } else if (reff >= endSurfaceAreaOff) {
        // Include a tapered surface area term.
        double taperSA = surfaceAreaSwitch.taper(reff, reff2, reff3, reff4, reff5);
        double dtaperSA = surfaceAreaSwitch.dtaper(reff, reff2, reff3, reff4) * dRdSA;
        cavitationEnergy += taperSA * surfaceAreaEnergy;
        double factor = dtaperSA * surfaceAreaEnergy + taperSA * surfaceTension;
        addSurfaceAreaGradient(factor, surfaceAreaGradient, gradient);
      }
    }

    if (logger.isLoggable(Level.FINE)) {
      logger.fine("\n GaussVol");
      logger.fine(format(" Volume:              %8.3f (Ang^3)", volume));
      logger.fine(format(" Volume Energy:       %8.3f (kcal/mol)", volumeEnergy));
      logger.fine(format(" Surface Area:        %8.3f (Ang^2)", surfaceArea));
      logger.fine(format(" Surface Area Energy: %8.3f (kcal/mol)", surfaceAreaEnergy));
      logger.fine(format(" Volume + SA Energy:  %8.3f (kcal/mol)", cavitationEnergy));
      logger.fine(format(" Effective Radius:    %8.3f (Ang)", reff));
    }

    return cavitationEnergy;
  }

  public ConnollyRegion getConnollyRegion() {
    return connollyRegion;
  }

  public double getCrossOver() {
    return crossOver;
  }

  public void setCrossOver(double crossOver) {
    if (crossOver < 3.5) {
      logger.severe(
          format(" The cross-over point (%8.6f A) must be greater than 3.5 A", crossOver));
      return;
    }
    this.crossOver = crossOver;
    beginVolumeOff = crossOver - HALF_SWITCH_RANGE;
    endVolumeOff = beginVolumeOff + 2.0 * HALF_SWITCH_RANGE;
    beginSurfaceAreaOff = crossOver + SA_SWITCH_OFFSET + HALF_SWITCH_RANGE;
    endSurfaceAreaOff = beginSurfaceAreaOff - 2.0 * HALF_SWITCH_RANGE;
    volumeSwitch = new MultiplicativeSwitch(beginVolumeOff, endVolumeOff);
    surfaceAreaSwitch = new MultiplicativeSwitch(beginSurfaceAreaOff, endSurfaceAreaOff);
  }

  public double getEffectiveRadius() {
    return effectiveRadius;
  }

  public double getEnergy() {
    return cavitationEnergy;
  }

  public GaussVol getGaussVol() {
    return gaussVol;
  }

  public double getSolventPressure() {
    return solventPressure;
  }

  public void setSolventPressure(double solventPressure) {
    double newCrossOver = 3.0 * surfaceTension / solventPressure;
    if (newCrossOver < 3.5) {
      logger.severe(
          format(
              " The solvent pressure (%8.6f kcal/mol/A^3)"
                  + " and surface tension (%8.6f kcal/mol/A^2) combination is not supported.",
              solventPressure, surfaceTension));
      return;
    }
    this.solventPressure = solventPressure;
    this.setCrossOver(newCrossOver);
  }

  /**
   * Return Surface Area (A^2).
   *
   * @return Surface Area (A^2).
   */
  public double getSurfaceArea() {
    return surfaceArea;
  }

  /**
   * Return Surface Area based cavitation energy.
   *
   * @return Surface Area based cavitation energy.
   */
  public double getSurfaceAreaEnergy() {
    return surfaceAreaEnergy;
  }

  public double getSurfaceTension() {
    return surfaceTension;
  }

  public void setSurfaceTension(double surfaceTension) {
    double newCrossOver = 3.0 * surfaceTension / solventPressure;
    if (newCrossOver < 3.5) {
      logger.severe(
          format(
              " The solvent pressure (%8.6f kcal/mol/A^3)"
                  + " and surface tension (%8.6f kcal/mol/A^2) combination is not supported.",
              solventPressure, surfaceTension));
      return;
    }
    this.surfaceTension = surfaceTension;
    this.setCrossOver(newCrossOver);
  }

  /**
   * Return Volume (A^3).
   *
   * @return Volume (A^3).
   */
  public double getVolume() {
    return volume;
  }

  /**
   * Return Volume based cavitation energy.
   *
   * @return Volume based cavitation energy.
   */
  public double getVolumeEnergy() {
    return volumeEnergy;
  }

  /**
   * Collect the volume base cavitation energy and gradient contributions.
   *
   * @param factorSA Factor to multiply surface area gradient by.
   * @param factorVol Factor to multiply surface area gradient by.
   * @param surfaceAreaGradient Surface area gradient.
   * @param volumeGradient Volume gradient.
   * @param gradient Array to accumulate derivatives.
   */
  private void addVolumeGradient(
      double factorSA,
      double factorVol,
      double[] surfaceAreaGradient,
      double[] volumeGradient,
      AtomicDoubleArray3D gradient) {
    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      double gx = factorSA * surfaceAreaGradient[index] + factorVol * volumeGradient[index++];
      double gy = factorSA * surfaceAreaGradient[index] + factorVol * volumeGradient[index++];
      double gz = factorSA * surfaceAreaGradient[index] + factorVol * volumeGradient[index];
      gradient.add(0, i, gx, gy, gz);
    }
  }

  /**
   * Collect the surface area based cavitation energy and gradient contributions.
   *
   * @param factor Factor to multiply surface area gradient by.
   * @param surfaceAreaGradient Surface area gradient.
   * @param gradient Array to accumulate derivatives.
   */
  private void addSurfaceAreaGradient(
      double factor, double[] surfaceAreaGradient, AtomicDoubleArray3D gradient) {
    for (int i = 0; i < nAtoms; i++) {
      int index = i * 3;
      double gx = factor * surfaceAreaGradient[index++];
      double gy = factor * surfaceAreaGradient[index++];
      double gz = factor * surfaceAreaGradient[index];
      gradient.add(0, i, gx, gy, gz);
    }
  }
}
