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
package ffx.potential.terms;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import ffx.numerics.Potential;
import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.parameters.ForceField;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

import static ffx.potential.parameters.ForceField.toEnumForm;
import static java.lang.String.format;

public class EnergyTermRegion extends ParallelRegion {

  private static final Logger logger = Logger.getLogger(EnergyTermRegion.class.getName());

  /**
   * The atoms in the MolecularAssembly.
   */
  private final Atom[] atoms;
  /**
   * The number of atoms in the MolecularAssembly.
   */
  private final int nAtoms;
  /**
   * If true, the gradient will be computed.
   */
  private boolean gradient = false;
  /**
   * If true, the gradient for each atom instance will be initialized to zero.
   */
  private boolean initAtomGradients = true;

  /**
   * The state of the potential energy term, which can be either
   * BOTH, FAST or SLOW. If SLOW, no bonded terms are evaluated.
   */
  private Potential.STATE state = Potential.STATE.BOTH;

  /*
   * The logic here deals with how DualTopologyEnergy (DTE) works.
   * If this ForceFieldEnergy (FFE) is not part of a DTE, lambdaBondedTerms is always false, and the term is always evaluated.
   *
   * Most bonded terms should be at full-strength, regardless of lambda, "complemented" to 1.0*strength.
   * Outside FFE, DTE scales the initial, !lambdaBondedTerms evaluation by f(lambda).
   *
   * In the case of a bonded term lacking any softcore atoms, this will be externally complemented by the other FFE topology.
   * If there is internal lambda-scaling, that will separately apply at both ends.
   *
   * In the case of a bonded term with a softcore atom, it's a bit trickier.
   * If it's unscaled by lambda, it needs to be internally complemented; we re-evaluate it with lambdaBondedTerms true.
   * This second evaluation is scaled by DTE by a factor of f(1-lambda), and becomes "restraintEnergy".
   * If it is scaled internally by lambda, we assume that the energy term is not meant to be internally complemented.
   * In that case, we skip evaluation into restraintEnergy.
   */

  /**
   * Alchemical atoms will not be checked for restraints.
   */
  protected boolean checkAlchemicalAtoms = true;
  /**
   * Indicates only bonded energy terms effected by Lambda should be evaluated.
   */
  protected boolean lambdaBondedTerms = false;
  /**
   * Indicates all bonded energy terms should be evaluated if lambdaBondedTerms is true.
   */
  protected boolean lambdaAllBondedTerms = false;
  /**
   * List of energy terms that are based on local bonded interactions.
   */
  private final List<EnergyTerm> energyTerms = new ArrayList<>();
  /**
   * The bonded energy terms will be executed in parallel using these loops.
   */
  private final BondedTermLoop[] bondedTermLoops;
  /**
   * Loops to initialize the gradient to zero.
   */
  private final GradInitLoop[] gradInitLoops;
  /**
   * Loops to reduce the gradient after computation.
   */
  private final GradReduceLoop[] gradReduceLoops;
  /**
   * The gradient for each atom, indexed by atom index and thread index.
   */
  private final AtomicDoubleArray3D grad;
  /**
   * The lambda gradient for each atom, indexed by atom index and thread index.
   * This is only used if lambdaTerm is true.
   */
  private final AtomicDoubleArray3D lambdaGrad;

  /**
   * Constructor for BondedRegion.
   *
   * @param molecularAssembly The MolecularAssembly that this region will operate on.
   * @param lambdaTerm        If true, the lambda gradient will be computed.
   * @param parallelTeam      The ParallelTeam that this region will run in.
   */
  public EnergyTermRegion(ParallelTeam parallelTeam,
                          MolecularAssembly molecularAssembly,
                          boolean lambdaTerm) {
    this.atoms = molecularAssembly.getAtomArray();
    this.nAtoms = atoms.length;

    // Define how the gradient will be accumulated.
    AtomicDoubleArray.AtomicDoubleArrayImpl atomicDoubleArrayImpl = AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI;
    ForceField forceField = molecularAssembly.getForceField();

    int nThreads = parallelTeam.getThreadCount();
    bondedTermLoops = new BondedTermLoop[nThreads];
    gradInitLoops = new GradInitLoop[nThreads];
    gradReduceLoops = new GradReduceLoop[nThreads];
    for (int i = 0; i < nThreads; i++) {
      bondedTermLoops[i] = new BondedTermLoop();
      gradInitLoops[i] = new GradInitLoop();
      gradReduceLoops[i] = new GradReduceLoop();
    }

    String value = forceField.getString("ARRAY_REDUCTION", "MULTI");
    try {
      atomicDoubleArrayImpl = AtomicDoubleArray.AtomicDoubleArrayImpl.valueOf(toEnumForm(value));
    } catch (Exception e) {
      logger.info(format(" Unrecognized ARRAY-REDUCTION %s; defaulting to %s", value,
          atomicDoubleArrayImpl));
    }
    logger.fine(format("  Bonded using %s arrays.", atomicDoubleArrayImpl));

    int nAtoms = molecularAssembly.getAtomArray().length;
    grad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, nThreads);
    if (lambdaTerm) {
      lambdaGrad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, nThreads);
    } else {
      lambdaGrad = null;
    }
  }

  public void setCheckAlchemicalAtoms(boolean checkAlchemicalAtoms) {
    this.checkAlchemicalAtoms = checkAlchemicalAtoms;
  }

  public boolean getCheckAlchemicalAtoms() {
    return checkAlchemicalAtoms;
  }

  public void setLambdaBondedTerms(boolean lambdaBondedTerms) {
    this.lambdaBondedTerms = lambdaBondedTerms;
  }

  public boolean getLambdaBondedTerms() {
    return lambdaBondedTerms;
  }

  public void setLambdaAllBondedTerms(boolean lambdaAllBondedTerms) {
    this.lambdaAllBondedTerms = lambdaAllBondedTerms;
  }

  public boolean getLambdaAllBondedTerms() {
    return lambdaAllBondedTerms;
  }

  public void setInitAtomGradients(boolean initAtomGradients) {
    this.initAtomGradients = initAtomGradients;
  }

  public boolean getInitAtomGradients() {
    return initAtomGradients;
  }

  /**
   * Get the total energy of all EnergyTerms.
   *
   * @return The total energy computed by summing the energies of all EnergyTerms.
   */
  public double getEnergy() {
    double energy = 0.0;
    for (EnergyTerm term : energyTerms) {
      energy += term.getEnergy();
    }
    return energy;
  }

  /**
   * Log the details of the energy terms in this bonded region.
   */
  public void log() {
    for (EnergyTerm term : energyTerms) {
      term.log();
    }
  }

  /**
   * String representation of this EnergyTermRegion.
   *
   * @return A string containing the details of all energy terms.
   */
  public String toString() {
    StringBuilder sb = new StringBuilder();
    for (EnergyTerm term : energyTerms) {
      sb.append(term.toString());
    }
    return sb.toString();
  }

  /**
   * String representation for PDB Headers
   *
   * @return A string containing the details of all energy terms.
   */
  public String toPDBString() {
    StringBuilder sb = new StringBuilder();
    for (EnergyTerm term : energyTerms) {
      sb.append(term.toPDBString());
    }
    return sb.toString();
  }

  /**
   * Set the state of the potential energy term.
   *
   * @param state The new state, which can be either BOTH, FAST or SLOW.
   */
  public void setState(Potential.STATE state) {
    if (state == null) {
      throw new IllegalArgumentException("Potential state cannot be null.");
    }
    this.state = state;
  }

  /**
   * Set whether the gradient will be computed.
   *
   * @param gradient If true, compute the gradient.
   */
  public void setGradient(boolean gradient) {
    this.gradient = gradient;
  }

  /**
   * Get whether the gradient will be computed.
   *
   * @return true if the gradient will be computed.
   */
  public boolean getGradient() {
    return gradient;
  }

  /**
   * Add a EnergyTerm to this bonded region.
   *
   * @param term EnergyTerm to add (ignored if null).
   * @return true if the term was added.
   */
  public boolean addEnergyTerm(EnergyTerm term) {
    if (term == null) {
      return false;
    }
    return energyTerms.add(term);
  }

  /**
   * Remove an EnergyTerm from this bonded region.
   *
   * @param term EnergyTerm to remove (ignored if null).
   * @return true if the term was present and removed.
   */
  public boolean removeEnergyTerm(EnergyTerm term) {
    if (term == null) {
      return false;
    }
    return energyTerms.remove(term);
  }

  /**
   * Get the BondedEnergyTerm at a specific index.
   *
   * @param index Index into the bonded energy terms list.
   * @return EnergyTerm at the specified index.
   * @throws IndexOutOfBoundsException if index is invalid.
   */
  public EnergyTerm getEnergyTerm(int index) {
    return energyTerms.get(index);
  }

  /**
   * Get an unmodifiable view of the energy terms.
   *
   * @return Unmodifiable list of EnergyTerm.
   */
  public List<EnergyTerm> getEnergyTerms() {
    return Collections.unmodifiableList(energyTerms);
  }

  @Override
  public void start() {
    // Initialize the BondedEnergyTerms.
    for (EnergyTerm term : energyTerms) {
      term.setEnergy(0.0);
      term.setRMSD(0.0);
    }
  }

  @Override
  public void run() throws Exception {
    int threadIndex = getThreadIndex();

    // Initialize the Gradient to zero.
    if (gradient) {
      execute(0, nAtoms - 1, gradInitLoops[threadIndex]);
    }

    // Determine if the bonded energy terms should be evaluated.
    if (state == Potential.STATE.BOTH || state == Potential.STATE.FAST) {
      // Loop over all BondedEnergyTerms and execute them in parallel.
      BondedTermLoop bondedTermLoop = bondedTermLoops[threadIndex];
      for (EnergyTerm term : energyTerms) {
        if (threadIndex == 0) {
          term.startTime();
        }
        bondedTermLoop.setEnergyTerm(term);
        execute(0, term.getNumberOfTerms() - 1, bondedTermLoop);
        if (threadIndex == 0) {
          term.stopTime();
        }
      }
    }

    // Reduce the gradients if they were computed.
    if (gradient) {
      execute(0, nAtoms - 1, gradReduceLoops[threadIndex]);
    }

  }

  private class BondedTermLoop extends IntegerForLoop {

    private EnergyTerm energyTerm;

    public void setEnergyTerm(EnergyTerm terms) {
      this.energyTerm = terms;
    }

    @Override
    public void run(int first, int last) {
      int threadIndex = getThreadIndex();
      double localEnergy = 0.0;
      double localRMSD = 0.0;

      BondedTerm[] bondedTerms = energyTerm.getBondedTermsArray();
      for (int i = first; i <= last; i++) {
        BondedTerm term = bondedTerms[i];
        boolean used = true;
        /*
         * The logic here deals with how DualTopologyEnergy (DTE) works.
         * If this ForceFieldEnergy (FFE) is not part of a DTE, lambdaBondedTerms is always false, and the term is always evaluated.
         *
         * Most bonded terms should be at full-strength, regardless of lambda, "complemented" to 1.0*strength.
         * Outside FFE, DTE scales the initial, !lambdaBondedTerms evaluation by f(lambda).
         *
         * In the case of a bonded term lacking any softcore atoms, this will be externally complemented by the other FFE topology.
         * If there is internal lambda-scaling, that will separately apply at both ends.
         *
         * In the case of a bonded term with a softcore atom, it's a bit trickier.
         * If it's unscaled by lambda, it needs to be internally complemented; we re-evaluate it with lambdaBondedTerms true.
         * This second evaluation is scaled by DTE by a factor of f(1-lambda), and becomes "restraintEnergy".
         * If it is scaled internally by lambda, we assume that the energy term is not meant to be internally complemented.
         * In that case, we skip evaluation into restraintEnergy.
         */
        if (checkAlchemicalAtoms) {
          used = !lambdaBondedTerms || lambdaAllBondedTerms || (term.applyLambda() && !term.isLambdaScaled());
        }
        if (used) {
          localEnergy += term.energy(gradient, threadIndex, grad, lambdaGrad);
          double value = term.getValue();
          localRMSD += value * value;
        }
      }

      energyTerm.addAndGetEnergy(localEnergy);
      energyTerm.addAndGetRMSD(localRMSD);
    }

  }

  private class GradInitLoop extends IntegerForLoop {

    @Override
    public void run(int first, int last) {
      int threadID = getThreadIndex();
      if (gradient) {
        grad.reset(threadID, first, last);
        if (initAtomGradients) {
          for (int i = first; i <= last; i++) {
            atoms[i].setXYZGradient(0.0, 0.0, 0.0);
          }
        }
      }
      if (lambdaGrad != null) {
        lambdaGrad.reset(threadID, first, last);
        if (initAtomGradients) {
          for (int i = first; i <= last; i++) {
            atoms[i].setLambdaXYZGradient(0.0, 0.0, 0.0);
          }
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }
  }

  private class GradReduceLoop extends IntegerForLoop {

    @Override
    public void run(int first, int last) {
      if (gradient) {
        grad.reduce(first, last);
        for (int i = first; i <= last; i++) {
          Atom a = atoms[i];
          a.addToXYZGradient(grad.getX(i), grad.getY(i), grad.getZ(i));
        }
      }
      if (lambdaGrad != null) {
        lambdaGrad.reduce(first, last);
        for (int i = first; i <= last; i++) {
          Atom a = atoms[i];
          a.addToLambdaXYZGradient(lambdaGrad.getX(i), lambdaGrad.getY(i), lambdaGrad.getZ(i));
        }
      }
    }

    @Override
    public IntegerSchedule schedule() {
      return IntegerSchedule.fixed();
    }
  }
}
