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
package ffx.algorithms.ph;

import static ffx.potential.extended.TitrationUtils.findResiduePolymer;
import static ffx.potential.extended.TitrationUtils.inactivateResidue;
import static java.lang.String.format;

import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.MultiTerminus;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.extended.ExtendedVariable;
import ffx.potential.extended.TitrationESV;
import ffx.potential.extended.TitrationUtils;
import ffx.potential.extended.TitrationUtils.TitrationConfig;
import ffx.potential.parameters.ForceField;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * PhMD class.
 *
 * @author Stephen D. LuCore
 */
public class PhMD  {
  private static final Logger logger = Logger.getLogger(PhMD.class.getName());

  /** Disallow code paths other than that which is being validated. DEBUG: reset to false. */
  private static final boolean CAUTIOUS = true;

  private static final double NS_TO_SEC = 0.000000001;
  /** Boltzmann's constant is kcal/mol/Kelvin. */
  private static final double BOLTZMANN = 0.0019872041;
  /** The MolecularAssembly. */
  private final MolecularAssembly molecularAssembly;
  /** Simulation pH. */
  private final double pH;
  /** The forcefield being used. Needed by MultiResidue constructor. */
  private final ForceField forceField;
  /**
   * The ForceFieldEnergy object being used by MD. Needed by MultiResidue constructor and for
   * reinitializing after a chemical change.
   */
  private final ForceFieldEnergy forceFieldEnergy;

  private final TitrationConfig titrationConfig;
  /** The MolecularDynamics object controlling the simulation. */
  private MolecularDynamics molecularDynamics;
  /** The current MD step. */
  private int stepCount = 0;
  /** Residues selected by user. */
  private final List<Residue> chosenResidues;

  private final List<MultiTerminus> titratingTermini = new ArrayList<>();
  private final List<ExtendedVariable> titratingESVs = new ArrayList<>();

  private ExtendedSystem esvSystem;

  /**
   * Construct a Monte-Carlo protonation state switching mechanism.
   *
   * @param molecularAssembly the molecular assembly
   * @param pH the simulation pH
   * @param molecularDynamics a {@link MolecularDynamics} object.
   * @param titrating a {@link java.util.List} object.
   */
  public PhMD(
      MolecularAssembly molecularAssembly,
      MolecularDynamics molecularDynamics,
      List<Residue> titrating,
      double pH) {
    this.titrationConfig = new TitrationConfig();
    this.molecularAssembly = molecularAssembly;
    this.molecularDynamics = molecularDynamics;
    this.forceField = molecularAssembly.getForceField();
    this.forceFieldEnergy = molecularAssembly.getPotentialEnergy();
    this.pH = pH;
    this.chosenResidues = titrating;

    reInitialize(true, false);
    readyup();
  }

  private void readyup() {
    // Create MultiTerminus objects to wrap termini.
    if (titrationConfig.titrateTermini) {
      for (Residue res : molecularAssembly.getResidueList()) {
        if (res.getPreviousResidue() == null || res.getNextResidue() == null) {
          MultiTerminus multiTerminus =
              new MultiTerminus(res, forceField, forceFieldEnergy, molecularAssembly);
          Polymer polymer = findResiduePolymer(res, molecularAssembly);
          polymer.addMultiTerminus(res, multiTerminus);
          reInitialize(true, false);
          titratingTermini.add(multiTerminus);
          logger.info(String.format(" Titrating: %s", multiTerminus));
        }
      }
    }

    /* Create containers for titratables: MultiResidues for discrete, ExtendedVariables for continuous. */
    esvSystem = new ExtendedSystem(molecularAssembly);
    esvSystem.setConstantPh(pH);
    for (Residue res : chosenResidues) {
        MultiResidue multi = TitrationUtils.titratingMultiresidueFactory(molecularAssembly, res);
        TitrationESV esv = new TitrationESV(esvSystem, multi);
        titratingESVs.add(esv);
        for (Residue background : multi.getInactive()) {
          inactivateResidue(background);
        }
        esvSystem.addVariable(esv);
      }
    forceFieldEnergy.attachExtendedSystem(esvSystem);
    logger.info(format(" Continuous pHMD readied with %d residues.", titratingESVs.size()));
    molecularDynamics.attachExtendedSystem(esvSystem, 1000);
  }

  /** Wraps reinitialization calls so as to provide un-fucked atom numbering. */
  private void reInitialize(boolean initFFE, boolean initMolDyn) {
    //        renumberAtoms(mola);	// TODO Determine if+why this is necessary.
    if (initFFE) {
      forceFieldEnergy.reInit();
    }
    if (initMolDyn) {
      molecularDynamics.reInit();
    }
  }



}
