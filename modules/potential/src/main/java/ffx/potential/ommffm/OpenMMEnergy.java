// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
package ffx.potential.ommffm;

import ffx.crystal.Crystal;
import ffx.openmm.ffm.bindings.OpenMMNative;
import ffx.potential.FiniteDifferenceUtils;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Platform;
import ffx.potential.bonded.Atom;
import ffx.potential.parameters.ForceField;
import ffx.potential.utils.EnergyException;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Double.isFinite;
import static java.lang.String.format;

/**
 * Compute the potential energy and derivatives using OpenMM via the FFM backend.
 */
public class OpenMMEnergy extends ForceFieldEnergy implements OpenMMPotential {

  private static final Logger logger = Logger.getLogger(OpenMMEnergy.class.getName());

  private final Platform platform;
  private OpenMMContext openMMContext;
  private OpenMMSystem openMMSystem;
  private final Atom[] atoms;
  private final boolean computeDEDL;

  public OpenMMEnergy(MolecularAssembly molecularAssembly, Platform requestedPlatform, int nThreads) {
    super(molecularAssembly, nThreads);

    Crystal crystal = getCrystal();
    int symOps = crystal.spaceGroup.getNumberOfSymOps();
    if (symOps > 1) {
      logger.severe(" OpenMM does not support symmetry operators.");
    }

    logger.info("\n Initializing OpenMM (FFM)");

    ForceField forceField = molecularAssembly.getForceField();
    atoms = molecularAssembly.getAtomArray();

    this.platform = requestedPlatform;
    ffx.openmm.ffm.Platform openMMPlatform = OpenMMContext.loadPlatform(platform, forceField);

    openMMSystem = new OpenMMSystem(this);
    openMMSystem.addForces();

    openMMContext = new OpenMMContext(openMMPlatform, openMMSystem, atoms);
    computeDEDL = forceField.getBoolean("OMM_DUDL", false);
  }

  @Override
  public double energy(double[] x) {
    return energy(x, false);
  }

  @Override
  public double energy(double[] x, boolean verbose) {
    if (lambdaBondedTerms) {
      return 0.0;
    }

    openMMContext.update();
    updateParameters(atoms);

    unscaleCoordinates(x);
    setCoordinates(x);

    OpenMMState openMMState = openMMContext.getOpenMMState(OpenMMNative.OpenMM_State_Energy());
    double e = openMMState.potentialEnergy;
    openMMState.destroy();

    if (!isFinite(e)) {
      String message = String.format(" Energy from OpenMM was a non-finite %8g", e);
      logger.warning(message);
      throw new EnergyException(message);
    }

    if (verbose) {
      logger.log(Level.INFO, String.format("\n OpenMM Energy: %14.10g", e));
    }

    scaleCoordinates(x);
    return e;
  }

  public double energyFFX(double[] x) {
    return super.energy(x, false);
  }

  public double energyFFX(double[] x, boolean verbose) {
    return super.energy(x, verbose);
  }

  public double energyAndGradientFFX(double[] x, double[] g) {
    return super.energyAndGradient(x, g, false);
  }

  public double energyAndGradientFFX(double[] x, double[] g, boolean verbose) {
    return super.energyAndGradient(x, g, verbose);
  }

  @Override
  public double energyAndGradient(double[] x, double[] g) {
    return energyAndGradient(x, g, false);
  }

  @Override
  public double energyAndGradient(double[] x, double[] g, boolean verbose) {
    if (lambdaBondedTerms) {
      return 0.0;
    }

    unscaleCoordinates(x);
    openMMContext.update();
    updateParameters(atoms);
    setCoordinates(x);

    OpenMMState openMMState = openMMContext.getOpenMMState(
        OpenMMNative.OpenMM_State_Energy() | OpenMMNative.OpenMM_State_Forces());
    double e = openMMState.potentialEnergy;
    g = openMMState.getGradient(g);
    openMMState.destroy();

    if (!isFinite(e)) {
      String message = format(" Energy from OpenMM was a non-finite %8g", e);
      logger.warning(message);
      throw new EnergyException(message);
    }

    scaleCoordinates(x);
    return e;
  }

  @Override
  public OpenMMContext getContext() {
    return openMMContext;
  }

  @Override
  public void updateContext(String integratorName, double timeStep, double temperature, boolean forceCreation) {
    openMMContext.update(integratorName, timeStep, temperature, forceCreation);
  }

  @Override
  public OpenMMState getOpenMMState(int mask) {
    return openMMContext.getOpenMMState(mask);
  }

  @Override
  public OpenMMSystem getSystem() {
    return openMMSystem;
  }

  @Override
  public double getd2EdL2() {
    return 0.0;
  }

  @Override
  public double getdEdL() {
    if (!lambdaTerm || !computeDEDL) {
      return 0.0;
    }
    return FiniteDifferenceUtils.computedEdL(this, this, molecularAssembly.getForceField());
  }

  @Override
  public void getdEdXdL(double[] gradients) {
  }

  @Override
  public boolean setActiveAtoms() {
    return openMMSystem.updateAtomMass();
  }

  @Override
  public void setCoordinates(double[] x) {
    super.setCoordinates(x);

    int n = atoms.length * 3;
    double[] xall = new double[n];
    int i = 0;
    for (Atom atom : atoms) {
      xall[i] = atom.getX();
      xall[i + 1] = atom.getY();
      xall[i + 2] = atom.getZ();
      i += 3;
    }
    openMMContext.setPositions(xall);
  }

  @Override
  public void setVelocity(double[] v) {
    super.setVelocity(v);

    int n = atoms.length * 3;
    double[] vall = new double[n];
    double[] v3 = new double[3];
    int i = 0;
    for (Atom atom : atoms) {
      atom.getVelocity(v3);
      vall[i] = v3[0];
      vall[i + 1] = v3[1];
      vall[i + 2] = v3[2];
      i += 3;
    }
    openMMContext.setVelocities(vall);
  }

  @Override
  public void setCrystal(Crystal crystal) {
    super.setCrystal(crystal);
    openMMContext.setPeriodicBoxVectors(crystal);
  }

  @Override
  public void setLambda(double lambda) {
    if (!lambdaTerm) {
      logger.fine(" Attempting to set lambda for an OpenMMEnergy with lambdaterm false.");
      return;
    }

    super.setLambda(lambda);

    if (atoms != null) {
      List<Atom> atomList = new ArrayList<>();
      for (Atom atom : atoms) {
        if (atom.applyLambda()) {
          atomList.add(atom);
        }
      }
      updateParameters(atomList.toArray(new Atom[0]));
    } else {
      updateParameters(null);
    }
  }

  @Override
  public void updateParameters(@Nullable Atom[] atoms) {
    if (atoms == null) {
      atoms = this.atoms;
    }
    if (openMMSystem != null) {
      openMMSystem.updateParameters(atoms);
    }
  }

  @Override
  public boolean destroy() {
    free();
    return super.destroy();
  }

  public void free() {
    if (openMMContext != null) {
      openMMContext.free();
      openMMContext = null;
    }
    if (openMMSystem != null) {
      openMMSystem.free();
      openMMSystem = null;
    }
  }
}
