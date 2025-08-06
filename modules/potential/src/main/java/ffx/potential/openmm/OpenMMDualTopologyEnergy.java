// ******************************************************************************
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
// ******************************************************************************
package ffx.potential.openmm;

import ffx.crystal.Crystal;
import ffx.numerics.switching.UnivariateSwitchingFunction;
import ffx.potential.DualTopologyEnergy;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Platform;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.utils.EnergyException;

import javax.annotation.Nullable;
import java.util.logging.Level;
import java.util.logging.Logger;

import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static java.lang.Double.isFinite;
import static java.lang.String.format;

public class OpenMMDualTopologyEnergy extends DualTopologyEnergy implements OpenMMPotential {

  private static final Logger logger = Logger.getLogger(OpenMMDualTopologyEnergy.class.getName());

  /**
   * OpenMM Context.
   */
  private final OpenMMContext openMMContext;
  /**
   * OpenMMDualTopologySystem.
   */
  private final OpenMMDualTopologySystem openMMDualTopologySystem;
  /**
   * The atoms this OpenMMEnergy operates on.
   */
  private final Atom[] atoms;
  /**
   * MolecularAssembly for topology 1.
   */
  private final MolecularAssembly molecularAssembly1;
  /**
   * MolecularAssembly for topology 2.
   */
  private final MolecularAssembly molecularAssembly2;

  /**
   * Constructor for DualTopologyEnergy.
   *
   * @param topology1         a {@link MolecularAssembly} object.
   * @param topology2         a {@link MolecularAssembly} object.
   * @param switchFunction    a {@link UnivariateSwitchingFunction} object.
   * @param requestedPlatform a {@link Platform} object.
   */
  public OpenMMDualTopologyEnergy(MolecularAssembly topology1, MolecularAssembly topology2,
                                  UnivariateSwitchingFunction switchFunction, Platform requestedPlatform) {
    super(topology1, topology2, switchFunction);

    logger.info("\n Initializing an OpenMM Dual Topology System.");

    molecularAssembly1 = topology1;
    molecularAssembly2 = topology2;

    // Check that the two topologies do not use symmetry operators.
    Crystal crystal1 = molecularAssembly1.getCrystal();
    int symOps1 = crystal1.spaceGroup.getNumberOfSymOps();
    Crystal crystal2 = molecularAssembly2.getCrystal();
    int symOps2 = crystal2.spaceGroup.getNumberOfSymOps();
    if (symOps1 > 1 || symOps2 > 1) {
      logger.severe(" OpenMM does not support symmetry operators.");
    }

    // Load the OpenMM plugins
    ForceField forceField = topology1.getForceField();
    ffx.openmm.Platform openMMPlatform = OpenMMContext.loadPlatform(requestedPlatform, forceField);

    // Create the OpenMM System.
    openMMDualTopologySystem = new OpenMMDualTopologySystem(this);
    openMMDualTopologySystem.addForces();

    atoms = getDualTopologyAtoms(0);

    // Create the Context.
    openMMContext = new OpenMMContext(openMMPlatform, openMMDualTopologySystem, atoms);
  }


  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x) {
    return energy(x, false);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double energy(double[] x, boolean verbose) {
    // Make sure the context has been created.
    openMMContext.update();

    updateParameters(atoms);

    // Unscale and set the coordinates.
    unscaleCoordinates(x);
    setCoordinates(x);

    OpenMMState openMMState = openMMContext.getOpenMMState(OpenMM_State_Energy);
    double e = openMMState.potentialEnergy;
    openMMState.destroy();

    if (!isFinite(e)) {
      String message = format(" OpenMMDualTopologyEnergy was a non-finite %8g", e);
      logger.warning(message);
      throw new EnergyException(message);
    }

    if (verbose) {
      logger.log(Level.INFO, format("\n OpenMM Energy: %14.10g", e));
    }

    return e;
  }

  /**
   * Compute the energy using the pure Java code path.
   *
   * @param x Atomic coordinates.
   * @return The energy (kcal/mol)
   */
  public double energyFFX(double[] x) {
    return super.energy(x, false);
  }

  /**
   * Compute the energy using the pure Java code path.
   *
   * @param x       Input atomic coordinates
   * @param verbose Use verbose logging.
   * @return The energy (kcal/mol)
   */
  public double energyFFX(double[] x, boolean verbose) {
    return super.energy(x, verbose);
  }

  /**
   * Set FFX and OpenMM coordinates for active atoms.
   *
   * @param x Atomic coordinates.
   */
  public void setCoordinates(double[] x) {
    // Set both OpenMM and FFX coordinates to x.
    openMMContext.setPositions(x);
  }

  /**
   * Get the MolecularAssembly.
   *
   * @param topology The topology index (0 for topology 1, 1 for topology 2).
   * @return a {@link MolecularAssembly} object for topology 1.
   */
  public MolecularAssembly getMolecularAssembly(int topology) {
    if (topology == 0) {
      return molecularAssembly1;
    } else if (topology == 1) {
      return molecularAssembly2;
    } else {
      throw new IllegalArgumentException(" Invalid topology index: " + topology);
    }
  }

  /**
   * Get the OpenMMEnergy for the specified topology.
   *
   * @param topology The topology index (0 for topology 1, 1 for topology 2).
   * @return The OpenMMEnergy for the specified topology.
   */
  public ForceFieldEnergy getForceFieldEnergy(int topology) {
    if (topology == 0) {
      return getForceFieldEnergy1();
    } else if (topology == 1) {
      return getForceFieldEnergy2();
    } else {
      throw new IllegalArgumentException(" Invalid topology index: " + topology);
    }
  }

  /**
   * Update parameters if the Use flags and/or Lambda value has changed.
   *
   * @param atoms Atoms in this list are considered.
   */
  @Override
  public void updateParameters(@Nullable Atom[] atoms) {
    if (atoms == null) {
      atoms = this.atoms;
    }
    if (openMMDualTopologySystem != null) {
      openMMDualTopologySystem.updateParameters(atoms);
    }
  }

  /**
   * Returns the Context instance.
   *
   * @return context
   */
  @Override
  public OpenMMContext getContext() {
    return openMMContext;
  }

  /**
   * Create an OpenMM Context.
   *
   * <p>Context.free() must be called to free OpenMM memory.
   *
   * @param integratorName Integrator to use.
   * @param timeStep       Time step.
   * @param temperature    Temperature (K).
   * @param forceCreation  Force a new Context to be created, even if the existing one matches the
   *                       request.
   */
  @Override
  public void updateContext(String integratorName, double timeStep, double temperature, boolean forceCreation) {
    openMMContext.update(integratorName, timeStep, temperature, forceCreation);
  }

  /**
   * Create an immutable OpenMM State.
   *
   * <p>State.free() must be called to free OpenMM memory.
   *
   * @param mask The State mask.
   * @return Returns the State.
   */
  @Override
  public OpenMMState getOpenMMState(int mask) {
    return openMMContext.getOpenMMState(mask);
  }

  /**
   * Get a reference to the System instance.
   *
   * @return a reference to the OpenMMSystem.
   */
  @Override
  public OpenMMSystem getSystem() {
    return openMMDualTopologySystem;
  }

  /**
   * Update active atoms.
   */
  @Override
  public void setActiveAtoms() {
    openMMDualTopologySystem.updateAtomMass();
    // Tests show reinitialization of the OpenMM Context is not necessary to pick up mass changes.
    // context.reinitContext();
  }

}
