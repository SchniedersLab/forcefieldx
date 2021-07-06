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
package ffx.potential.extended;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.*;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.PDBFilter.Mutation;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.Constants;
import org.apache.commons.io.FilenameUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalDouble;
import java.util.logging.Logger;

import static ffx.potential.extended.ExtUtils.prop;
import static java.lang.String.format;

/**
 * Helper methods to define titration-specific phenomena.
 *
 * @author Stephen LuCore
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("serial")
public class TitrationUtils {

  /** Constant <code>heavyStrandedDynamics=prop("phmd-heavyStrandedDynamics", false)</code> */
  public static final boolean heavyStrandedDynamics = prop("phmd-heavyStrandedDynamics", false);

  /** Constant <code>threeStateHistidines=prop("phmd-threeState", false)</code> */
  public static final boolean threeStateHistidines =
      prop("phmd-threeState", false); // not yet implemented

  private static final Logger logger = Logger.getLogger(TitrationUtils.class.getName());

  /** Utility class */
  private TitrationUtils() {}

  /**
   * activateResidue.
   *
   * @param addDoF a {@link ffx.potential.bonded.Residue} object.
   */
  public static void activateResidue(Residue addDoF) {
    List<Atom> atomList = addDoF.getAtomList();
    for (Atom atom : atomList) {
      atom.setActive(true);
    }
  }

  /**
   * Identify titratable residues and choose them all.
   *
   * @param searchMe a {@link ffx.potential.MolecularAssembly} object.
   * @return a {@link java.util.List} object.
   */
  public static List<Residue> chooseTitratables(MolecularAssembly searchMe) {
    List<Residue> chosen = new ArrayList<>();
    Polymer polymers[] = searchMe.getChains();
    for (int i = 0; i < polymers.length; i++) {
      List<Residue> residues = polymers[i].getResidues();
      for (int j = 0; j < residues.size(); j++) {
        Residue res = residues.get(j);
        Titration[] avail = Titration.multiLookup(res);
        if (avail != null) {
          chosen.add(residues.get(j));
        }
      }
    }
    return chosen;
  }

  /**
   * Choose titratables with intrinsic pKa inside (pH-window,pH+window).
   *
   * @param pH a double.
   * @param window a double.
   * @param searchMe a {@link ffx.potential.MolecularAssembly} object.
   * @return a {@link java.util.List} object.
   */
  public static List<Residue> chooseTitratables(
      double pH, double window, MolecularAssembly searchMe) {
    List<Residue> chosen = new ArrayList<>();
    Polymer polymers[] = searchMe.getChains();
    for (int i = 0; i < polymers.length; i++) {
      List<Residue> residues = polymers[i].getResidues();
      for (int j = 0; j < residues.size(); j++) {
        Residue res = residues.get(j);
        Titration[] avail = Titration.multiLookup(res);
        for (Titration titration : avail) {
          double pKa = titration.pKa;
          if (pKa >= pH - window && pKa <= pH + window) {
            chosen.add(residues.get(j));
          }
        }
      }
    }
    return chosen;
  }

  /**
   * chooseTitratables.
   *
   * @param residueIDs a {@link java.lang.String} object.
   * @param searchMe a {@link ffx.potential.MolecularAssembly} object.
   * @return a {@link java.util.List} object.
   */
  public static List<Residue> chooseTitratables(String residueIDs, MolecularAssembly searchMe) {
    String[] tokens =
        (residueIDs.split(".").length > 1)
            ? residueIDs.split(".")
            : (residueIDs.split(",").length > 1)
                ? residueIDs.split(",")
                : new String[] {residueIDs};
    return chooseTitratables(Arrays.asList(tokens), searchMe);
  }

  /**
   * Select titrating residues by amino acid.
   *
   * @param aa a {@link ffx.potential.bonded.ResidueEnumerations.AminoAcid3} object.
   * @param searchMe a {@link ffx.potential.MolecularAssembly} object.
   * @return a {@link java.util.List} object.
   */
  public static List<Residue> chooseTitratables(AminoAcid3 aa, MolecularAssembly searchMe) {
    List<Residue> chosen = new ArrayList<>();
    Polymer polymers[] = searchMe.getChains();
    for (Polymer polymer : polymers) {
      List<Residue> residues = polymer.getResidues();
      for (Residue res : residues) {
        if (res.getAminoAcid3() == aa) {
          Titration[] avail = Titration.multiLookup(res);
          if (avail != null) {
            chosen.add(res);
          }
        }
      }
    }
    return chosen;
  }

  /**
   * chooseTitratables.
   *
   * @param crIDs a {@link java.util.List} object.
   * @param searchMe a {@link ffx.potential.MolecularAssembly} object.
   * @return a {@link java.util.List} object.
   */
  public static List<Residue> chooseTitratables(List<String> crIDs, MolecularAssembly searchMe) {
    List<Residue> chosen = new ArrayList<>();
    for (String crID : crIDs) {
      char chain = crID.charAt(0);
      int num = Integer.parseInt(crID.substring(1));
      boolean found = false;
      List<Residue> allRes = searchMe.getResidueList();
      for (Residue res : allRes) {
        if (res.getChainID() == chain && res.getResidueNumber() == num) {
          chosen.add(res);
          found = true;
          break;
        }
      }
      if (!found) {
        logger.severe(format("Couldn't find residue for crID %c,%d.", chain, num));
      }
    }
    return chosen;
  }

  /**
   * chooseTitratables.
   *
   * @param chain a char.
   * @param resID a int.
   * @param searchMe a {@link ffx.potential.MolecularAssembly} object.
   * @return a {@link java.util.List} object.
   */
  public static List<Residue> chooseTitratables(char chain, int resID, MolecularAssembly searchMe) {
    List<Residue> chosen = new ArrayList<>();
    Polymer polymers[] = searchMe.getChains();
    for (Polymer polymer : polymers) {
      if (polymer.getChainID() == chain) {
        List<Residue> residues = polymer.getResidues();
        for (Residue residue : residues) {
          if (residue.getResidueNumber() == resID) {
            chosen.add(residue);
            logger.info(String.format(" Chosen: %s", residue));
          }
        }
      }
    }
    return chosen;
  }

  /**
   * Locate to which Polymer in a MolecularAssembly the given Residue belongs.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param mola a {@link ffx.potential.MolecularAssembly} object.
   * @return a {@link ffx.potential.bonded.Polymer} object.
   */
  public static Polymer findResiduePolymer(Residue residue, MolecularAssembly mola) {
    if (residue.getChainID() == null) {
      logger.severe("No chain ID for residue " + residue);
    }
    Polymer polymers[] = mola.getChains();
    Polymer location = null;
    for (Polymer p : polymers) {
      if (p.getChainID().equals(residue.getChainID())) {
        location = p;
      }
    }
    if (location == null) {
      logger.severe("Couldn't find polymer for residue " + residue);
    }
    return location;
  }

  /**
   * inactivateResidue.
   *
   * @param killDoF a {@link ffx.potential.bonded.Residue} object.
   */
  public static void inactivateResidue(Residue killDoF) {
    List<Atom> atomList = killDoF.getAtomList();
    for (Atom atom : atomList) {
      atom.setActive(false);
    }
  }

  /**
   * initDiscountPreloadProperties.
   *
   * @param cutoffs a {@link java.lang.Double} object.
   */
  public static void initDiscountPreloadProperties(Double cutoffs) {
    initEsvPreloadProperties(cutoffs);
  }

  /** initDiscountPreloadProperties. */
  public static void initDiscountPreloadProperties() {
    initDiscountPreloadProperties(null);
  }

  /**
   * Note that this must (generally) be called before loading the input file or instantiating
   * titration classes.
   *
   * @param cutoffs a {@link java.lang.Double} object.
   */
  public static void initEsvPreloadProperties(Double cutoffs) {
    // Active Potential
    // System.setProperty("forcefield", "AMOEBA_PROTEIN_2013");
    System.setProperty("esvterm", "true");
    //System.setProperty("lambdaterm", "true");
    // System.setProperty("bondterm", "true");
    // System.setProperty("angleterm", "true");
    // System.setProperty("strbndterm", "true");
    // System.setProperty("ureyterm", "true");
    // System.setProperty("opbendterm", "true");
    // System.setProperty("torsionterm", "true");
    // System.setProperty("pitorsterm", "true");
    // System.setProperty("tortorterm", "true");
    // System.setProperty("improperterm", "true");

    // Optional Potential
    // System.setProperty("vdwterm", "true"); // van der Waals
    System.setProperty("esv.vdw", "true");
    // System.setProperty("mpoleterm", "true"); // permanent real space
    System.setProperty("pme-qi", "true");
    System.setProperty("esv.pme", "true");
    // System.setProperty("recipterm", "true"); // permanent reciprocal space

    // Inactive Potential
    // System.setProperty("polarizeterm", "false"); // polarization
    // System.setProperty("polarization", "NONE");
    // System.setProperty("gkterm", "false");
    // System.setProperty("restrainterm", "false");
    // System.setProperty("comrestrainterm", "false");
    // System.setProperty("lambda_torsions", "false");

    // Potential Settings
    System.setProperty("permanent-lambda-alpha", "2.0");
    System.setProperty("permanent-lambda-exponent", "3.0");

    // Polarize on the whole range [0,1]
    System.setProperty("polarization-lambda-start", "0.0");

    // Polarization is not soft-cored, only a factor of lambda is applied.
    System.setProperty("polarization-lambda-exponent", "0.0");

    // No special SCF between proton atoms that are being turned off.
    System.setProperty("ligand-vapor-elec", "false");

    // No special SCF for the protein without protons.
    System.setProperty("no-ligand-condensed-scf", "false");

    // Protons on different amino acids do not feel each other when turned off.
    System.setProperty("intramolecular-softcore", "true");

    // Protons on different proteins do not feel each other when turned off.
    System.setProperty("intermolecular-softcore", "true");

    // ESV Settings
    //        System.setProperty("esv.biasTerm", "true");             // include discretization and
    // pH biases
    //        System.setProperty("esv.scaleBonded", "true");          // include effects on bonded
    // terms
    //        System.setProperty("esv.backgroundBonded", "true");     // hook up BG bonded terms to
    // FG node
    //        System.setProperty("esv.scaleUnshared", "true");        // use multipole scaling in
    // all cases (eliminates softcoring)
  }

  /**
   * isTitratableHydrogen.
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @return a boolean.
   */
  public static boolean isTitratableHydrogen(Atom atom) {
    String name = atom.getName();
    switch (atom.getResidueName()) {
      case "LYS":
        if (name.equals("HZ3")) {
          return true;
        }
        break;
      case "TYR":
        if (name.equals("HH")) {
          return true;
        }
        break;
      case "CYS":
        if (name.equals("HG")) {
          return true;
        }
        break;
      case "HIS":
        if (name.equals("HD1") || name.equals("HE2")) {
          return true;
        }
        break;
      case "HID":
        if (name.equals("HE2")) {  //HD1?
          return true;
        }
        break;
      case "HIE":
        if (name.equals("HD1")) {  //HE2?
          return true;
        }
        break;
      case "ASH":
        if (name.equals("HD2")) {
          return true;
        }
        break;
      case "GLH":
        if (name.equals("HE2")) {
          return true;
        }
        break;
    }
    return false;
  }

  /**
   * openFullyProtonated.
   *
   * @param structure a {@link java.io.File} object.
   * @return a {@link ffx.potential.MolecularAssembly} object.
   */
  public static MolecularAssembly openFullyProtonated(File structure) {
    String name = format("%s-prot", FilenameUtils.removeExtension(structure.getName()));
    MolecularAssembly mola = new MolecularAssembly(name);
    mola.setFile(structure);

    List<Mutation> mutations = new ArrayList<>();
    List<Residue> residues = mola.getResidueList();
    for (Residue res : residues) {
      char chain = res.getChainID();
      int resID = res.getResidueNumber();
      Titration titration = Titration.lookup(res);
      if (res.getAminoAcid3() != titration.protForm) {
        String protName = titration.protForm.name();
        mutations.add(new PDBFilter.Mutation(chain, resID, protName));
      }
    }

    PotentialsUtils utils = new PotentialsUtils();
    return utils.openWithMutations(structure, mutations);
  }

  /**
   * openFullyProtonated.
   *
   * @param filename a {@link java.lang.String} object.
   * @return a {@link ffx.potential.MolecularAssembly} object.
   */
  public static MolecularAssembly openFullyProtonated(String filename) {
    return openFullyProtonated(new File(filename));
  }

  /**
   * propagateInactiveResidues.
   *
   * @param multiResidues a {@link java.util.List} object.
   * @param propagateDynamics a boolean.
   */
  public static void propagateInactiveResidues(
      List<MultiResidue> multiResidues, boolean propagateDynamics) {
    for (MultiResidue multiRes : multiResidues) {
      propagateInactiveResidues(multiRes, propagateDynamics);
    }
  }

  /**
   * propagateInactiveResidues.
   *
   * @param multiResidues a {@link java.util.List} object.
   */
  public static void propagateInactiveResidues(List<MultiResidue> multiResidues) {
    for (MultiResidue multiRes : multiResidues) {
      propagateInactiveResidues(multiRes, true);
    }
  }

  /**
   * propagateInactiveResidues.
   *
   * @param multiResidue a {@link ffx.potential.bonded.MultiResidue} object.
   */
  public static void propagateInactiveResidues(MultiResidue multiResidue) {
    propagateInactiveResidues(multiResidue, true);
  }

  /**
   * Copies atomic coordinates from each active residue to its inactive counterparts. Inactive
   * hydrogen coordinates are updated by geometry with the propagated heavies.
   *
   * @param multiRes a {@link ffx.potential.bonded.MultiResidue} object.
   * @param propagateDynamics a boolean.
   */
  public static void propagateInactiveResidues(MultiResidue multiRes, boolean propagateDynamics) {
    // Propagate all atom coordinates from active residues to their inactive counterparts.
    Residue active = multiRes.getActive();
    List<Residue> inactives = multiRes.getInactive();
    for (Atom activeAtom : active.getAtomList()) {
      String activeName = activeAtom.getName();
      for (Residue inactive : inactives) {
        Atom inactiveAtom = (Atom) inactive.getAtomNode(activeName);
        if (inactiveAtom != null) {
          // Propagate position and gradient.
          double[] activeXYZ = activeAtom.getXYZ(null);
          inactiveAtom.setXYZ(activeXYZ);
          double[] grad = new double[3];
          activeAtom.getXYZGradient(grad);
          inactiveAtom.setXYZGradient(grad[0], grad[1], grad[2]);
          if (propagateDynamics) {
            // Propagate velocity, acceleration, and previous acceleration.
            double[] activeVelocity = new double[3];
            activeAtom.getVelocity(activeVelocity);
            inactiveAtom.setVelocity(activeVelocity);
            double[] activeAccel = new double[3];
            activeAtom.getAcceleration(activeAccel);
            inactiveAtom.setAcceleration(activeAccel);
            double[] activePrevAcc = new double[3];
            activeAtom.getPreviousAcceleration(activePrevAcc);
            inactiveAtom.setPreviousAcceleration(activePrevAcc);
          }
        } else {
          if (activeName.equals("C")
              || activeName.equals("O")
              || activeName.equals("N")
              || activeName.equals("CA")
              || activeName.equals("H")
              || activeName.equals("HA")) {
            // Backbone atoms aren't supposed to exist in inactive multiResidue components; so no
            // problem.
          } else if (isTitratableHydrogen(activeAtom)) {
            /**
             * i.e. ((activeResName.equals("LYS") && activeName.equals("HZ3")) ||
             * (activeResName.equals("TYR") && activeName.equals("HH")) ||
             * (activeResName.equals("CYS") && activeName.equals("HG")) ||
             * (activeResName.equals("HIS") && (activeName.equals("HD1") ||
             * activeName.equals("HE2"))) || (activeResName.equals("HID") &&
             * activeName.equals("HD1")) || (activeResName.equals("HIE") &&
             * activeName.equals("HE2")) || (activeResName.equals("ASH") &&
             * activeName.equals("HD2")) || (activeResName.equals("GLH") &&
             * activeName.equals("HE2")))
             */
            // These titratable protons are handled below; so no problem.
          } else {
            // Now we have a problem.
            logger.warning(
                format(
                    "Couldn't propagate inactive MultiResidue atom: %s: %s, %s",
                    multiRes, activeName, activeAtom));
          }
        }
      }
    }
    rebuildStrandedProtons(multiRes);
  }

  /**
   * renumberAtoms.
   *
   * @param mola a {@link ffx.potential.MolecularAssembly} object.
   */
  public static void renumberAtoms(MolecularAssembly mola) {
    Atom[] atoms = mola.getAtomArray();
    int setter = 0;
    for (Atom atom : atoms) {
      atom.setXyzIndex(setter++);
    }
  }

  /**
   * Create a MultiResidue from the given Residue by adding its alternated protonation state(s) as
   * alternate possibilities.
   *
   * @param mola a {@link ffx.potential.MolecularAssembly} object.
   * @param res a {@link ffx.potential.bonded.Residue} object.
   * @return a {@link ffx.potential.bonded.MultiResidue} object.
   */
  public static MultiResidue titratingMultiresidueFactory(MolecularAssembly mola, Residue res) {
    ForceField ff = mola.getForceField();
    Potential potential = mola.getPotentialEnergy();
    if (!(potential instanceof ForceFieldEnergy)) {
      logger.warning(
          String.format("TitrationFactory only supported by ForceFieldEnergy potentials."));
      throw new IllegalStateException();
    }
    ForceFieldEnergy ffe = (ForceFieldEnergy) potential;

    /* Create new titration state. */
    Titration titration = Titration.lookup(res);
    String targetName =
        (titration.protForm != res.getAminoAcid3())
            ? titration.protForm.toString()
            : titration.deprotForm.toString();
    int resNumber = res.getResidueNumber();
    Residue.ResidueType resType = res.getResidueType();
    Residue newRes = new Residue(targetName, resNumber, resType);

    /* Wrap both states in a MultiResidue. */
    MultiResidue multiRes = new MultiResidue(res, ff, ffe);
    Polymer polymer = findResiduePolymer(res, mola);
    polymer.addMultiResidue(multiRes);
    multiRes.addResidue(newRes);

    /* Begin in protonated state by default. */
    multiRes.setActiveResidue(titration.protForm);
    propagateInactiveResidues(multiRes, false);
    ffe.reInit();
    return multiRes;
  }

  /**
   * Rebuild stranded titratable protons from ideal geometry. "Stranded protons" are titrating H+
   * atoms on inactive MultiRes members; when propagating new coordinates to inactive residues, no
   * coords/velocity exist for them.
   */
  private static void rebuildStrandedProtons(MultiResidue multiRes) {
    // If inactive residue is a protonated form, move the stranded hydrogen to new coords (based on
    // propagated heavies).
    // Also give the stranded hydrogen a maxwell velocity and remove its accelerations.
    List<Residue> inactives = multiRes.getInactive();
    for (Residue inactive : inactives) {
      List<Atom> resetMe = new ArrayList<>();
      switch (inactive.getName()) {
        case "LYS":
          {
            Atom HZ3 = (Atom) inactive.getAtomNode("HZ3");
            Atom NZ = (Atom) inactive.getAtomNode("NZ");
            Atom CE = (Atom) inactive.getAtomNode("CE");
            Atom HZ1 = (Atom) inactive.getAtomNode("HZ1");
            BondedUtils.intxyz(HZ3, NZ, 1.02, CE, 109.5, HZ1, 109.5, -1);
            resetMe.add(HZ3);
            break;
          }
        case "ASH":
          {
            Atom HD2 = (Atom) inactive.getAtomNode("HD2");
            Atom OD2 = (Atom) inactive.getAtomNode("OD2");
            Atom CG = (Atom) inactive.getAtomNode("CG");
            Atom OD1 = (Atom) inactive.getAtomNode("OD1");
            BondedUtils.intxyz(HD2, OD2, 0.98, CG, 108.7, OD1, 0.0, 0);
            resetMe.add(HD2);
            break;
          }
        case "GLH":
          {
            Atom HE2 = (Atom) inactive.getAtomNode("HE2");
            Atom OE2 = (Atom) inactive.getAtomNode("OE2");
            Atom CD = (Atom) inactive.getAtomNode("CD");
            Atom OE1 = (Atom) inactive.getAtomNode("OE1");
            BondedUtils.intxyz(HE2, OE2, 0.98, CD, 108.7, OE1, 0.0, 0);
            resetMe.add(HE2);
            break;
          }
        case "HIS":
          {
            Atom HE2 = (Atom) inactive.getAtomNode("HE2");
            Atom NE2 = (Atom) inactive.getAtomNode("NE2");
            Atom CD2 = (Atom) inactive.getAtomNode("CD2");
            Atom CE1 = (Atom) inactive.getAtomNode("CE1");
            Atom HD1 = (Atom) inactive.getAtomNode("HD1");
            Atom ND1 = (Atom) inactive.getAtomNode("ND1");
            Atom CG = (Atom) inactive.getAtomNode("CG");
            Atom CB = (Atom) inactive.getAtomNode("CB");
            BondedUtils.intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
            BondedUtils.intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
            resetMe.add(HE2);
            resetMe.add(HD1);
            break;
          }
        case "HID":
          {
            Atom HD1 = (Atom) inactive.getAtomNode("HD1");
            Atom ND1 = (Atom) inactive.getAtomNode("ND1");
            Atom CG = (Atom) inactive.getAtomNode("CG");
            Atom CB = (Atom) inactive.getAtomNode("CB");
            BondedUtils.intxyz(HD1, ND1, 1.02, CG, 126.0, CB, 0.0, 0);
            resetMe.add(HD1);
            break;
          }
        case "HIE":
          {
            Atom HE2 = (Atom) inactive.getAtomNode("HE2");
            Atom NE2 = (Atom) inactive.getAtomNode("NE2");
            Atom CD2 = (Atom) inactive.getAtomNode("CD2");
            Atom CE1 = (Atom) inactive.getAtomNode("CE1");
            BondedUtils.intxyz(HE2, NE2, 1.02, CD2, 126.0, CE1, 126.0, 1);
            resetMe.add(HE2);
            break;
          }
        case "CYS":
          {
            Atom HG = (Atom) inactive.getAtomNode("HG");
            Atom SG = (Atom) inactive.getAtomNode("SG");
            Atom CB = (Atom) inactive.getAtomNode("CB");
            Atom CA = (Atom) inactive.getAtomNode("CA");
            BondedUtils.intxyz(HG, SG, 1.34, CB, 96.0, CA, 180.0, 0);
            resetMe.add(HG);
            break;
          }
        case "TYR":
          {
            Atom HH = (Atom) inactive.getAtomNode("HH");
            Atom OH = (Atom) inactive.getAtomNode("OH");
            Atom CZ = (Atom) inactive.getAtomNode("CZ");
            Atom CE2 = (Atom) inactive.getAtomNode("CE2");
            BondedUtils.intxyz(HH, OH, 0.97, CZ, 108.0, CE2, 0.0, 0);
            resetMe.add(HH);
            break;
          }
        default:
      }
      for (Atom a : resetMe) {
        if (heavyStrandedDynamics) {
          // Use of heavy atom dynamics properties is in testing.
          a.setXYZGradient(0, 0, 0);
          double[] heavyVelocity = new double[3];
          double[] heavyAccel = new double[3];
          double[] heavyPrevAccel = new double[3];
          Atom heavy = a.getBonds().get(0).get1_2(a);
          heavy.getVelocity(heavyVelocity);
          heavy.getAcceleration(heavyAccel);
          heavy.getPreviousAcceleration(heavyPrevAccel);
          a.setVelocity(heavyVelocity);
          a.setAcceleration(heavyAccel);
          a.setPreviousAcceleration(heavyPrevAccel);
        } else {
          // PREVIOUSLY: draw vel from maxwell and set accel to zero
          a.setXYZGradient(0, 0, 0);
          a.setVelocity(ExtUtils.maxwellVelocity(a.getMass(), Constants.ROOM_TEMPERATURE));
          a.setAcceleration(new double[] {0, 0, 0});
          a.setPreviousAcceleration(new double[] {0, 0, 0});
        }
      }
    }
  }

  /** How DISCOUNT initializes lambda values at outset of continuous dynamics. */
  public enum ContinuousSeedDistribution {
    FLAT,
    BETA,
    BOLTZMANN,
    DIRAC_CURRENT,
    DIRAC_POINTFIVE;
  }

  /** Global override of MC acceptance criteria. */
  public enum MCOverride {
    NONE,
    ACCEPT,
    REJECT,
    ONCE;
  }

  /**
   * Writes .s-[num] and .f-[num] files representing before/after MC move structures. Note: The
   * 'after' snapshot represents the change that was PROPOSED, regardless of accept/reject.
   */
  public enum Snapshots {
    INTERLEAVED,
    SEPARATE,
    NONE;
  }

  public enum HistidineMode {
    HIE_ONLY,
    HID_ONLY,
    SINGLE,
    DOUBLE;
  }

  /**
   * Amino acid protonation reactions. Constructors below specify intrinsic pKa and reference free
   * energy of protonation, obtained via (OST) metadynamics on capped monomers.
   */
  public enum Titration {
    ctoC(8.18, 60.168, 0.0, AminoAcid3.CYD, AminoAcid3.CYS),
    Dtod(3.90, 53.188, 0.0, AminoAcid3.ASP, AminoAcid3.ASH),
    Etoe(4.25, 59.390, 0.0, AminoAcid3.GLU, AminoAcid3.GLH),
    ktoK(10.53, -52.807,0.077955, AminoAcid3.LYD, AminoAcid3.LYS),
    ytoY(10.07, 34.961, 0.0, AminoAcid3.TYD, AminoAcid3.TYR),
    UtoH(6.00, -42.923, 0.0, AminoAcid3.HID, AminoAcid3.HIS),
    ZtoH(6.00, 00.000, 0.0, AminoAcid3.HIE, AminoAcid3.HIS),
    TerminalNH3toNH2(8.23, 0.0, 0.0, AminoAcid3.UNK, AminoAcid3.UNK),
    TerminalCOOHtoCOO(3.55, 0.0, 0.0, AminoAcid3.UNK, AminoAcid3.UNK);

    public final double pKa;
    public final double refEnergy;
    public final double lambdaIntercept;
    public final AminoAcid3 protForm;
    public final AminoAcid3 deprotForm;

    /** Invoked by Enum; use the factory method to obtain instances. */
    private Titration(double pKa, double refEnergy, double lambdaIntercept, AminoAcid3 deprotForm, AminoAcid3 protForm) {
      this.pKa = pKa;
      this.refEnergy = refEnergy;
      this.lambdaIntercept = lambdaIntercept;
      this.deprotForm = deprotForm;
      this.protForm = protForm;
    }

    public static Titration lookup(Residue res) {
      Titration[] titrations = multiLookup(res);
      if (titrations.length > 1) {
        logger.warning(
            "Titration::lookup returned more results than expected. Did you mean to invoke multi-state?");
      }
      return (titrations != null) ? titrations[0] : null;
    }

    /**
     * Return a Titration object for the given Residue. TODO: Add support for multi-state titrations
     * (HIS,ASP,GLU).
     *
     * @param res a Residue.
     * @return a Titration instance.
     */
    public static Titration[] multiLookup(Residue res) {
      AminoAcid3 current = AminoAcid3.valueOf(res.getName());

      if (threeStateHistidines) {
        if (current == AminoAcid3.HIS || current == AminoAcid3.HID || current == AminoAcid3.HIE) {
          return new Titration[] {ZtoH, UtoH};
        }
      }

      for (Titration titration : Titration.values()) {
        if (current == titration.protForm || current == titration.deprotForm) {
          return new Titration[] {titration};
        }
      }

      logger.warning(format("No titration lookup found for residue %s", res));
      return null;
    }
  }

  public enum TitrationType {
    PROT,
    DEPROT;
  }

  /** Advanced options to both DiscreteMCMD and DiscountPh. */
  public static class TitrationConfig {

    public final ContinuousSeedDistribution seedDistribution =
        prop("phmd-seedMode", ContinuousSeedDistribution.class, ContinuousSeedDistribution.FLAT);
    public final Snapshots snapshots = prop("phmd-snapshots", Snapshots.class, Snapshots.NONE);
    public final HistidineMode histidineMode =
        prop("phmd-histidineMode", HistidineMode.class, HistidineMode.HIE_ONLY);
    public final OptionalDouble referenceOverride =
        prop("phmd-referenceOverride", OptionalDouble.empty());
    public final double meltdownTemperature = prop("phmd-meltdownTemp", 6000.0);
    public final double warningTemperature = prop("phmd-warningTemp", 1000.0);
    public final boolean logTimings = prop("phmd-logTimings", false);
    public final boolean titrateTermini = prop("phmd-termini", false);
    public final boolean zeroReferences = prop("phmd-zeroReferences", true);
    public final int debugLogLevel = prop("phmd-debugLog", 0);
    public final boolean useConformationalBias = prop("phmd-cbmcRotamerMoves", false);
    public final boolean inactivateBackground = prop("phmd-inactivateBackground", false);
    public final boolean zeroReferenceEnergies =
        prop("phmd-zeroReferences", false, "Zeroing all reference energies!");
    public final OptionalDouble refOverride =
        prop(
            "phmd-refOverride",
            OptionalDouble.empty(),
            "Reference protonation energies overridden!");
    public MCOverride mcOverride = prop("phmd-override", MCOverride.class, MCOverride.NONE);

    public void print() {
      ExtUtils.printConfigSet("Titration Config:", System.getProperties(), "phmd");
    }
  }
}
