// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.algorithms.dynamics;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.crystal.SpaceGroup;
import ffx.numerics.Potential;
import ffx.numerics.math.RunningStatistics;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parsers.XYZFilter;
import org.apache.commons.io.FilenameUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.crystal.LatticeSystem.check;
import static ffx.numerics.math.ScalarMath.mirrorDegrees;
import static ffx.utilities.Constants.*;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.*;

/**
 * The Barostat class maintains constant pressure using random trial moves in lattice parameters,
 * which are consistent with the space group.
 *
 * <p>Frenkel and Smit, "Understanding Molecular Simulation, 2nd Edition", Academic Press, San
 * Diego, CA, 2002; Section 5.4
 *
 * @author Michael J. Schnieders
 */
public class Barostat implements CrystalPotential {

  private static final Logger logger = Logger.getLogger(Barostat.class.getName());

  /**
   * MolecularAssembly being simulated.
   */
  private final MolecularAssembly molecularAssembly;
  /**
   * Boundary conditions and symmetry operators (could be a ReplicatedCrystal).
   */
  private Crystal crystal;
  /**
   * The unit cell.
   */
  private Crystal unitCell;
  /**
   * Mass of the system.
   */
  private final double mass;
  /**
   * ForceFieldEnergy that describes the system.
   */
  private final CrystalPotential potential;
  /**
   * Atomic coordinates.
   */
  private final double[] x;
  /**
   * The number of space group symmetry operators.
   */
  private final int nSymm;
  /**
   * Number of independent molecules in the simulation cell.
   */
  private final SpaceGroup spaceGroup;
  /**
   * Number of molecules in the systems.
   */
  private final int nMolecules;
  /**
   * Ideal gas constant * temperature (kcal/mol).
   */
  private double kT;
  /**
   * Sampling pressure (atm)
   */
  private double pressure = 1.0;
  /**
   * Only isotropic MC moves.
   */
  private boolean isotropic = false;
  /**
   * Flag to turn the Barostat on or off. If false, MC moves will not be tried.
   */
  private boolean active = true;
  /**
   * Default angle move (degrees).
   */
  private double maxAngleMove = 0.5;
  /**
   * Maximum volume move (Angstroms^3).
   */
  private double maxVolumeMove = 1.0;
  /**
   * Minimum density constraint.
   */
  private double minDensity = 0.75;
  /**
   * Maximum density constraint.
   */
  private double maxDensity = 1.60;
  /**
   * Number of energy evaluations between application of MC moves.
   */
  private int meanBarostatInterval = 1;
  /**
   * A counter for the number of barostat calls.
   */
  private int barostatCount = 0;
  /**
   * Number of unit cell Monte Carlo moves attempted.
   */
  private final RunningStatistics unitCellMoves = new RunningStatistics();
  /**
   * Number of lattice axis Monte Carlo moves attempted.
   */
  private final RunningStatistics sideMoves = new RunningStatistics();
  /**
   * Number of Monte Carlo moves attempted.
   */
  private final RunningStatistics angleMoves = new RunningStatistics();
  /**
   * Energy STATE.
   */
  private STATE state = STATE.BOTH;
  /**
   * Barostat move type.
   */
  private MoveType moveType = MoveType.SIDE;
  /**
   * True when a Monte Carlo move is accepted.
   */
  private boolean moveAccepted = false;
  /**
   * Current density value.
   */
  private double currentDensity = 0;
  /**
   * Current A-axis value.
   */
  private double a = 0;
  /**
   * Current B-axis value.
   */
  private double b = 0;
  /**
   * Current C-axis value.
   */
  private double c = 0;
  /**
   * Current alpha value.
   */
  private double alpha = 0;
  /**
   * Current beta value.
   */
  private double beta = 0;
  /**
   * Current gamma value.
   */
  private double gamma = 0;
  /**
   * A-axis statistics.
   */
  private final RunningStatistics aStats = new RunningStatistics();
  /**
   * B-axis statistics.
   */
  private final RunningStatistics bStats = new RunningStatistics();
  /**
   * C-axis statistics.
   */
  private final RunningStatistics cStats = new RunningStatistics();
  /**
   * Alpha statistics.
   */
  private final RunningStatistics alphaStats = new RunningStatistics();
  /**
   * Beta statistics.
   */
  private final RunningStatistics betaStats = new RunningStatistics();
  /**
   * Gamma statistics.
   */
  private final RunningStatistics gammaStats = new RunningStatistics();
  /**
   * Density statistics.
   */
  private final RunningStatistics densityStats = new RunningStatistics();
  /**
   * Barostat print frequency.
   */
  private int printFrequency = 1000;

  /**
   * Initialize the Barostat.
   *
   * @param molecularAssembly The molecular assembly to apply the MC barostat to.
   * @param potential         a {@link ffx.crystal.CrystalPotential} object.
   */
  public Barostat(MolecularAssembly molecularAssembly, CrystalPotential potential) {
    this(molecularAssembly, potential, 298.15);
  }

  /**
   * Initialize the Barostat.
   *
   * @param molecularAssembly The molecular assembly to apply the MC barostat to.
   * @param potential         a {@link ffx.crystal.CrystalPotential} object.
   * @param temperature       The Metropolis Monte Carlo temperature (Kelvin).
   */
  public Barostat(MolecularAssembly molecularAssembly, CrystalPotential potential,
                  double temperature) {

    this.molecularAssembly = molecularAssembly;
    this.potential = potential;
    this.kT = temperature * R;

    crystal = potential.getCrystal();
    unitCell = crystal.getUnitCell();
    spaceGroup = unitCell.spaceGroup;
    // Atoms in the system.
    Atom[] atoms = molecularAssembly.getAtomArray();
    // Number of atoms.
    int nAtoms = atoms.length;
    nSymm = spaceGroup.getNumberOfSymOps();
    mass = molecularAssembly.getMass();
    x = new double[3 * nAtoms];
    nMolecules = molecularAssembly.fractionalCount();
  }

  /**
   * density.
   *
   * @return a double.
   */
  public double density() {
    return (mass * nSymm / AVOGADRO) * (1.0e24 / unitCell.volume);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public boolean destroy() {
    // Nothing at this level to destroy.
    return potential.destroy();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {
    // Do not apply the Barostat for energy only evaluations.
    return potential.energy(x);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energyAndGradient(double[] x, double[] g) {

    // Calculate the energy and gradient as usual.
    double energy = potential.energyAndGradient(x, g);

    // Apply the barostat during computation of slowly varying forces.
    if (active && state != STATE.FAST) {
      if (random() < (1.0 / meanBarostatInterval)) {

        // Attempt to change the unit cell parameters.
        moveAccepted = false;

        applyBarostat(energy);

        // Collect Statistics.
        collectStats();

        // If a move was accepted, then re-calculate the gradient so
        // that it's consistent with the current unit cell parameters.
        if (moveAccepted) {
          energy = potential.energyAndGradient(x, g);
        }
      }
    }
    return energy;
  }

  /**
   * Restrict the MC Barostat to isotropic moves. The lattice angles are held fixed, and lattice
   * lengths are scaled equally.
   *
   * @return Returns true if only isotropic moves are allowed.
   */
  public boolean isIsotropic() {
    return isotropic;
  }

  /**
   * Restrict the MC Barostat to isotropic moves. The lattice angles are held fixed, and lattice
   * lengths are scaled equally.
   *
   * @param isotropic If true, if only isotropic moves are allowed.
   */
  public void setIsotropic(boolean isotropic) {
    this.isotropic = isotropic;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getAcceleration(double[] acceleration) {
    return potential.getAcceleration(acceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getCoordinates(double[] parameters) {
    return potential.getCoordinates(parameters);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public Crystal getCrystal() {
    return potential.getCrystal();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setCrystal(Crystal crystal) {
    potential.setCrystal(crystal);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public STATE getEnergyTermState() {
    return state;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setEnergyTermState(STATE state) {
    this.state = state;
    potential.setEnergyTermState(state);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getMass() {
    return potential.getMass();
  }

  /**
   * Returns the mean number of steps between barostat applications.
   *
   * @return Mean steps between barostat applications.
   */
  public int getMeanBarostatInterval() {
    return meanBarostatInterval;
  }

  /**
   * Setter for the field <code>meanBarostatInterval</code>.
   *
   * @param meanBarostatInterval The mean number of steps between barostat applications.
   */
  public void setMeanBarostatInterval(int meanBarostatInterval) {
    this.meanBarostatInterval = meanBarostatInterval;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public int getNumberOfVariables() {
    return potential.getNumberOfVariables();
  }

  /**
   * Gets the pressure of this Barostat in atm.
   *
   * @return Pressure in atm.
   */
  public double getPressure() {
    return pressure;
  }

  /**
   * Setter for the field <code>pressure</code>.
   *
   * @param pressure a double.
   */
  public void setPressure(double pressure) {
    this.pressure = pressure;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getPreviousAcceleration(double[] previousAcceleration) {
    return potential.getPreviousAcceleration(previousAcceleration);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getScaling() {
    return potential.getScaling();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setScaling(double[] scaling) {
    potential.setScaling(scaling);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return potential.getTotalEnergy();
  }

  @Override
  public List<Potential> getUnderlyingPotentials() {
    List<Potential> underlying = new ArrayList<>();
    underlying.add(potential);
    underlying.addAll(potential.getUnderlyingPotentials());
    return underlying;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public VARIABLE_TYPE[] getVariableTypes() {
    return potential.getVariableTypes();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double[] getVelocity(double[] velocity) {
    return potential.getVelocity(velocity);
  }

  public boolean isActive() {
    return active;
  }

  /**
   * Setter for the field <code>active</code>.
   *
   * @param active a boolean.
   */
  public void setActive(boolean active) {
    this.active = active;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setAcceleration(double[] acceleration) {
    potential.setAcceleration(acceleration);
  }

  /**
   * setDensity.
   *
   * @param density a double.
   */
  public void setDensity(double density) {
    molecularAssembly.computeFractionalCoordinates();
    crystal.setDensity(density, mass);
    potential.setCrystal(crystal);
    molecularAssembly.moveToFractionalCoordinates();
  }

  /**
   * Setter for the field <code>maxAngleMove</code>.
   *
   * @param maxAngleMove a double.
   */
  public void setMaxAngleMove(double maxAngleMove) {
    this.maxAngleMove = maxAngleMove;
  }

  /**
   * Setter for the field <code>maxVolumeMove</code>.
   *
   * @param maxVolumeMove a double.
   */
  public void setMaxVolumeMove(double maxVolumeMove) {
    this.maxVolumeMove = maxVolumeMove;
  }

  /**
   * Setter for the field <code>maxDensity</code>.
   *
   * @param maxDensity a double.
   */
  public void setMaxDensity(double maxDensity) {
    this.maxDensity = maxDensity;
  }

  /**
   * Setter for the field <code>minDensity</code>.
   *
   * @param minDensity a double.
   */
  public void setMinDensity(double minDensity) {
    this.minDensity = minDensity;
  }

  /**
   * Set the Barostat print frequency.
   *
   * @param frequency The number of Barostat moves between print statements.
   */
  public void setBarostatPrintFrequency(int frequency) {
    this.printFrequency = frequency;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setPreviousAcceleration(double[] previousAcceleration) {
    potential.setPreviousAcceleration(previousAcceleration);
  }

  /**
   * Set the Metropolis Monte Carlo temperature.
   *
   * @param temperature Temperature (Kelvin).
   */
  public void setTemperature(double temperature) {
    this.kT = temperature * R;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setVelocity(double[] velocity) {
    potential.setVelocity(velocity);
  }

  private double mcStep(double currentE, double currentV) {

    // Enforce minimum & maximum density constraints.
    double den = density();
    if(Double.isNaN(den) || Double.isNaN(a) || Double.isNaN(b) || Double.isNaN(c) || Double.isNaN(alpha) || Double.isNaN(beta) || Double.isNaN(gamma)){
      logger.warning(format(" Found value of NaN: density: %7.4f a: %7.4f b: %7.4f c: %7.4f alpha: %7.4f beta: %7.4f gamma: %7.4f",
              den, a, b, c, alpha, beta, gamma));
      File errorFile = new File(FilenameUtils.removeExtension(molecularAssembly.getFile().getName()) + "_err.xyz");
      XYZFilter.version(errorFile);
      XYZFilter writeFilter = new XYZFilter(errorFile, molecularAssembly, molecularAssembly.getForceField(), molecularAssembly.getProperties());
      writeFilter.writeFile(errorFile, true, null);
    }
    if (den < minDensity || den > maxDensity) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(
            format(" MC Density %10.6f is outside the range %10.6f - %10.6f.", den, minDensity,
                maxDensity));
      }
      // Reject moves outside the specified density range.
      crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
      return currentE;
    }

    // A carbon atom cannot fit into a unit cell without interfacial radii greater than ~1.2
    // Angstroms.
    double minInterfacialRadius = 1.2;
    double currentMinInterfacialRadius = min(min(unitCell.interfacialRadiusA, unitCell.interfacialRadiusB),
        unitCell.interfacialRadiusC);
    // Enforce minimum interfacial radii of 1.2 Angstroms.
    if (currentMinInterfacialRadius < minInterfacialRadius) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" MC An interfacial radius (%10.6f,%10.6f,%10.6f) is below the minimum %10.6f",
            unitCell.interfacialRadiusA, unitCell.interfacialRadiusB,
            unitCell.interfacialRadiusC, minInterfacialRadius));
      }
      // Reject due to small interfacial radius.
      crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
      return currentE;
    }

    // Apply the boundary condition for the proposed move. If the move is
    // rejected, then the previous boundary conditions must be restored.
    potential.setCrystal(crystal);

    // Update atomic coordinates to maintain molecular fractional centers of mass.
    molecularAssembly.moveToFractionalCoordinates();

    // Save the new volume
    double newV = unitCell.volume / nSymm;

    potential.getCoordinates(x);

    // Compute the new energy
    double newE = potential.energy(x);

    // Compute the change in potential energy
    double dPE = newE - currentE;

    // Compute the pressure-volume work for the asymmetric unit.
    double dEV = pressure * (newV - currentV) / PRESCON;

    // Compute the volume entropy
    double dES = -nMolecules * kT * log(newV / currentV);

    // Add up the contributions
    double dE = dPE + dEV + dES;

    // Energy decreased so the move is accepted.
    if (dE < 0.0) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" MC Accept: %12.6f (dPE) + %12.6f (dEV) + %12.6f (dES) = %12.6f with (V=%12.6f, E=%12.6f)",
            dPE, dEV, dES, dE, newV, newE));
      }
      moveAccepted = true;
      switch (moveType) {
        case SIDE -> sideMoves.addValue(1.0);
        case ANGLE -> angleMoves.addValue(1.0);
        case UNIT -> unitCellMoves.addValue(1.0);
      }
      return newE;
    }

    // Apply the Metropolis criteria.
    double acceptanceProbability = exp(-dE / kT);

    // Draw a pseudorandom double greater than or equal to 0.0 and less than 1.0
    // from an (approximately) uniform distribution from that range.
    double metropolis = random();

    // Energy increase without the Metropolis criteria satisfied.
    if (metropolis > acceptanceProbability) {
      rejectMove();
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" MC Reject: %12.6f (dPE) + %12.6f (dEV) + %12.6f (dES) = %12.6f with (V=%12.6f, E=%12.6f)",
            dPE, dEV, dES, dE, currentV, currentE));
      }
      switch (moveType) {
        case SIDE -> sideMoves.addValue(0.0);
        case ANGLE -> angleMoves.addValue(0.0);
        case UNIT -> unitCellMoves.addValue(0.0);
      }
      return currentE;
    }

    // Energy increase with Metropolis criteria satisfied.
    if (logger.isLoggable(Level.FINE)) {
      logger.fine(format(" MC Accept: %12.6f (dPE) + %12.6f (dEV) + %12.6f (dES) = %12.6f with (V=%12.6f, E=%12.6f)",
          dPE, dEV, dES, dE, newV, newE));
    }

    moveAccepted = true;
    switch (moveType) {
      case SIDE -> sideMoves.addValue(1.0);
      case ANGLE -> angleMoves.addValue(1.0);
      case UNIT -> unitCellMoves.addValue(1.0);
    }

    return newE;
  }

  /**
   * Reset the to state prior to trial move.
   */
  private void rejectMove() {
    // Reset the unit cell parameters
    crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);

    // Reset the potential PBC.
    potential.setCrystal(crystal);

    // Reset the atomic coordinates to maintain molecular fractional centers of mass.
    molecularAssembly.moveToFractionalCoordinates();
  }

  /**
   * Propose to change the A-axis length.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcA(double currentE) {
    moveType = MoveType.SIDE;
    double currentAUV = unitCell.volume / nSymm;
    double dAUVolume = maxVolumeMove * (2.0 * random() - 1.0);
    double dVdA = unitCell.dVdA / nSymm;
    double dA = dAUVolume / dVdA;
    boolean succeed = crystal.changeUnitCellParameters(a + dA, b, c, alpha, beta, gamma);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change of the a-axis to (%6.3f) with volume change %6.3f (A^3)", a, dAUVolume));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Propose to change the B-axis length.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcB(double currentE) {
    moveType = MoveType.SIDE;
    double currentAUV = unitCell.volume / nSymm;
    double dAUVolume = maxVolumeMove * (2.0 * random() - 1.0);
    double dVdB = unitCell.dVdB / nSymm;
    double dB = dAUVolume / dVdB;
    boolean succeed = crystal.changeUnitCellParameters(a, b + dB, c, alpha, beta, gamma);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change of the b-axis to (%6.3f) with volume change %6.3f (A^3)", b, dAUVolume));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Propose to change the C-axis length.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcC(double currentE) {
    moveType = MoveType.SIDE;
    double currentAUV = unitCell.volume / nSymm;
    double dAUVolume = maxVolumeMove * (2.0 * random() - 1.0);
    double dVdC = unitCell.dVdC / nSymm;
    double dC = dAUVolume / dVdC;
    boolean succeed = crystal.changeUnitCellParameters(a, b, c + dC, alpha, beta, gamma);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change to the c-axis to (%6.3f) with volume change %6.3f (A^3)", c, dAUVolume));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Propose the same change to the A-axis and B-axis lengths.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcAB(double currentE) {
    moveType = MoveType.SIDE;
    double currentAUV = unitCell.volume / nSymm;
    double dAUVolume = maxVolumeMove * (2.0 * random() - 1.0);
    double dVdAB = (unitCell.dVdA + unitCell.dVdB) / nSymm;
    double dAB = dAUVolume / dVdAB;
    boolean succeed = crystal.changeUnitCellParametersAndVolume(a + dAB, b + dAB, c, alpha, beta,
        gamma, currentAUV + dAUVolume);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change to the a,b-axes to (%6.3f) with volume change %6.3f (A^3)", a, dAUVolume));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  private double mcABC(double currentE) {
    moveType = MoveType.SIDE;
    double currentAUV = unitCell.volume / nSymm;
    double dAUVolume = maxVolumeMove * (2.0 * random() - 1.0);
    double dVdABC = (unitCell.dVdA + unitCell.dVdB + unitCell.dVdC) / nSymm;
    double dABC = dAUVolume / dVdABC;
    boolean succeed = crystal.changeUnitCellParametersAndVolume(a + dABC, b + dABC, c + dABC, alpha,
        beta, gamma, currentAUV + dAUVolume);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change to the a,b,c-axes to (%6.3f) with volume change %6.3f (A^3)", a, dAUVolume));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Propose to change the Alpha angle.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcAlpha(double currentE) {
    moveType = MoveType.ANGLE;
    double move = maxAngleMove * (2.0 * random() - 1.0);
    double currentAUV = unitCell.volume / nSymm;
    double newAlpha = mirrorDegrees(alpha + move);
    boolean succeed = crystal.changeUnitCellParametersAndVolume(a, b, c, newAlpha, beta, gamma, currentAUV);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change of alpha to %6.3f.", newAlpha));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Propose to change the Beta angle.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcBeta(double currentE) {
    moveType = MoveType.ANGLE;
    double move = maxAngleMove * (2.0 * random() - 1.0);
    double currentAUV = unitCell.volume / nSymm;
    double newBeta = mirrorDegrees(beta + move);
    boolean succeed = crystal.changeUnitCellParametersAndVolume(a, b, c, alpha, newBeta, gamma, currentAUV);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change of beta to %6.3f.", newBeta));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Propose to change the Gamma angle.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcGamma(double currentE) {
    moveType = MoveType.ANGLE;
    double move = maxAngleMove * (2.0 * random() - 1.0);
    double currentAUV = unitCell.volume / nSymm;
    double newGamma = mirrorDegrees(gamma + move);
    boolean succeed = crystal.changeUnitCellParametersAndVolume(a, b, c, alpha, beta, newGamma, currentAUV);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change of gamma to %6.3f.", newGamma));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Propose the same change to all 3 angles.
   *
   * @param currentE Current energy.
   * @return Energy after the MC trial.
   */
  private double mcABG(double currentE) {
    moveType = MoveType.ANGLE;
    double move = maxAngleMove * (2.0 * random() - 1.0);
    double currentAUV = unitCell.volume / nSymm;
    double newAlpha = mirrorDegrees(alpha + move);
    double newBeta = mirrorDegrees(beta + move);
    double newGamma = mirrorDegrees(gamma + move);
    boolean succeed = crystal.changeUnitCellParametersAndVolume(a, b, c, newAlpha, newBeta, newGamma, currentAUV);
    if (succeed) {
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(" Proposing MC change to all angles to %6.3f.", newAlpha));
      }
      return mcStep(currentE, currentAUV);
    }
    return currentE;
  }

  /**
   * Attempt an MC move on a lattice parameter.
   *
   * @param currentE The current potential energy.
   * @return The potential energy after attempting the MC move.
   */
  private double applyBarostat(double currentE) {
    // Determine the current molecular centers of mass in fractional coordinates.
    molecularAssembly.computeFractionalCoordinates();

    // Collect the current unit cell parameters.
    crystal = potential.getCrystal();
    unitCell = crystal.getUnitCell();
    a = unitCell.a;
    b = unitCell.b;
    c = unitCell.c;
    alpha = unitCell.alpha;
    beta = unitCell.beta;
    gamma = unitCell.gamma;

    if (isotropic) {
      currentE = mcABC(currentE);
    } else {
      switch (spaceGroup.latticeSystem) {
        case MONOCLINIC_LATTICE -> {
          // alpha = gamma = 90
          int move = (int) floor(random() * 4.0);
          switch (move) {
            case 0 -> currentE = mcA(currentE);
            case 1 -> currentE = mcB(currentE);
            case 2 -> currentE = mcC(currentE);
            case 3 -> currentE = mcBeta(currentE);
            default -> logger.severe(" Barostat programming error.");
          }
        }
        case ORTHORHOMBIC_LATTICE -> {
          // alpha = beta = gamma = 90
          int move = (int) floor(random() * 3.0);
          switch (move) {
            case 0 -> currentE = mcA(currentE);
            case 1 -> currentE = mcB(currentE);
            case 2 -> currentE = mcC(currentE);
            default -> logger.severe(" Barostat programming error.");
          }
        }
        case TETRAGONAL_LATTICE -> {
          // a = b, alpha = beta = gamma = 90
          int move = (int) floor(random() * 2.0);
          switch (move) {
            case 0 -> currentE = mcAB(currentE);
            case 1 -> currentE = mcC(currentE);
            default -> logger.severe(" Barostat programming error.");
          }
        }
        case RHOMBOHEDRAL_LATTICE -> {
          // a = b = c, alpha = beta = gamma
          int move = (int) floor(random() * 2.0);
          switch (move) {
            case 0 -> currentE = mcABC(currentE);
            case 1 -> currentE = mcABG(currentE);
            default -> logger.severe(" Barostat programming error.");
          }
        }
        case HEXAGONAL_LATTICE -> {
          // a = b, alpha = beta = 90, gamma = 120
          int move = (int) floor(random() * 2.0);
          switch (move) {
            case 0 -> currentE = mcAB(currentE);
            case 1 -> currentE = mcC(currentE);
            default -> logger.severe(" Barostat programming error.");
          }
        }
        case CUBIC_LATTICE ->
            // a = b = c, alpha = beta = gamma = 90
            currentE = mcABC(currentE);
        case TRICLINIC_LATTICE -> {
          if (check(a, b) && check(b, c) && check(alpha, 90.0) && check(beta, 90.0) && check(gamma, 90.0)) {
            currentE = mcABC(currentE);
          } else {
            int move = (int) floor(random() * 6.0);
            switch (move) {
              case 0 -> currentE = mcA(currentE);
              case 1 -> currentE = mcB(currentE);
              case 2 -> currentE = mcC(currentE);
              case 3 -> currentE = mcAlpha(currentE);
              case 4 -> currentE = mcBeta(currentE);
              case 5 -> currentE = mcGamma(currentE);
              default -> logger.severe(" Barostat programming error.");
            }
          }
        }
      }
    }

    currentDensity = density();
    if (moveAccepted) {
      if (logger.isLoggable(Level.FINE)) {
        StringBuilder sb = new StringBuilder(" MC Barostat Acceptance:");
        if (sideMoves.getCount() > 0) {
          sb.append(format(" Side %5.1f%%", sideMoves.getMean() * 100.0));
        }
        if (angleMoves.getCount() > 0) {
          sb.append(format(" Angle %5.1f%%", angleMoves.getMean() * 100.0));
        }
        if (unitCellMoves.getCount() > 0) {
          sb.append(format(" UC %5.1f%%", unitCellMoves.getMean() * 100.0));
        }
        sb.append(format("\n Density: %5.3f  UC: %s", currentDensity, unitCell.toShortString()));
        logger.fine(sb.toString());
      }
    } else {
      // Check that the unit cell parameters have not changed.
      if (unitCell.a != a || unitCell.b != b || unitCell.c != c || unitCell.alpha != alpha
          || unitCell.beta != beta || unitCell.gamma != gamma) {
        logger.severe(" Reversion of unit cell parameters did not succeed after failed Barostat MC move.");
      }
    }

    return currentE;
  }

  /**
   * Collect statistics on lattice parameters and density.
   */
  private void collectStats() {
    // Collect statistics.
    barostatCount++;

    // Sanity check for values.
    if(Double.isNaN(currentDensity)||Double.isNaN(a)||Double.isNaN(b)||Double.isNaN(c)||Double.isNaN(alpha)||Double.isNaN(beta)||Double.isNaN(gamma)){
      logger.warning(format(" Statistic Value was NaN: Density: %5.3f A: %5.3f B: %5.3f C: %5.3f Alpha: %5.3f Beta: %5.3f Gamma: %5.3f",
              currentDensity, a, b, c, alpha, beta, gamma));
    }
    densityStats.addValue(currentDensity);
    aStats.addValue(a);
    bStats.addValue(b);
    cStats.addValue(c);
    alphaStats.addValue(alpha);
    betaStats.addValue(beta);
    gammaStats.addValue(gamma);
    if (barostatCount % printFrequency == 0) {
      logger.info(format(" Barostat statistics for the last %d moves:", printFrequency));
      // Density
      logger.info(format("  Density: %5.3f±%5.3f with range %5.3f .. %5.3f (g/cc)", densityStats.getMean(),
          densityStats.getStandardDeviation(), densityStats.getMin(), densityStats.getMax()));
      densityStats.reset();
      // Lattice Parameters
      logger.info(format("  Lattice a-axis: %5.2f±%3.2f b-axis: %5.2f±%3.2f c-axis: %5.2f±%3.2f",
          aStats.getMean(), aStats.getStandardDeviation(), bStats.getMean(), bStats.getStandardDeviation(), cStats.getMean(),
          cStats.getStandardDeviation()));
      logger.info(format("          alpha:  %5.2f±%3.2f beta:   %5.2f±%3.2f gamma:  %5.2f±%3.2f",
          alphaStats.getMean(), alphaStats.getStandardDeviation(), betaStats.getMean(),
          betaStats.getStandardDeviation(), gammaStats.getMean(), gammaStats.getStandardDeviation()));
      aStats.reset();
      bStats.reset();
      cStats.reset();
      alphaStats.reset();
      betaStats.reset();
      gammaStats.reset();
      // MC Move Statistics
      if (sideMoves.getCount() > 0) {
        logger.info(format("  Axis length moves: %4d (%5.2f%%)", sideMoves.getCount(), sideMoves.getMean() * 100.0));
        sideMoves.reset();
      }
      if (angleMoves.getCount() > 0) {
        logger.info(format("  Angle moves:       %4d (%5.2f%%)", angleMoves.getCount(), angleMoves.getMean() * 100.0));
        angleMoves.reset();
      }
      if (unitCellMoves.getCount() > 0) {
        logger.info(format("  Unit cell moves:   %4d (%5.2f%%)", unitCellMoves.getCount(), unitCellMoves.getMean() * 100.0));
        unitCellMoves.reset();
      }
    }
  }

  /**
   * The type of Barostat move.
   */
  private enum MoveType {
    SIDE, ANGLE, UNIT
  }
}
