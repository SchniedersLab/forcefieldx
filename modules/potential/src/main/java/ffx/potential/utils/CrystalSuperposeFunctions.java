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
package ffx.potential.utils;

import static ffx.numerics.math.DoubleMath.dist;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.DoubleIndexPair;
import ffx.utilities.FFXScript;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * PACCOMFunctions contains methods utilized in the CrystalSuperpose.groovy file.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 */
public class CrystalSuperposeFunctions {

  /**
   * The logger for this class.
   */
  protected static final Logger logger = Logger.getLogger(FFXScript.class.getName());

  /**
   * Generate a molecular assembly of a sphere containing a specified number of molecules.
   *
   * @param molecularAssembly Base molecular assembly to expand.
   * @param nMolecules Number of molecules (asymmetric unit copies) in the expansion.
   * @param radius Approximate radius of the expanded system.
   * @param original If following the original version of PACCOM (by Okimasa Okada)
   * @param symOps List of SymOps used to generate the expanded system.
   * @return MolecularAssembly of the desired sphere size.
   */
  public static MolecularAssembly generateBaseSphere(MolecularAssembly molecularAssembly,
      int nMolecules, boolean original, double radius, List<SymOp> symOps) {

    // Move molecules into the unit cell.
    molecularAssembly.moveAllIntoUnitCell();

    // When the system was read in, a replicates crystal may have been created to satisfy the cutoff.
    // Retrieve a reference to the unit cell (not the replicates crystal).
    Crystal unitCell = molecularAssembly.getCrystal().getUnitCell();
    Crystal replicatesCrystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, radius * 2.0);

    logger.fine(" Base replicates crystal " + replicatesCrystal.toString());

    // Update the Molecular Assembly with the ReplicatesCrystal (used outside this method).
    molecularAssembly.setCrystal(replicatesCrystal);

    Atom[] atoms = molecularAssembly.getAtomArray();
    int nAtoms = atoms.length;
    // Collect asymmetric unit atomic coordinates.
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      x[i] = atom.getX();
      y[i] = atom.getY();
      z[i] = atom.getZ();
    }

    // Symmetry coordinates for each molecule in replicates crystal
    int nSymm = replicatesCrystal.getNumSymOps();
    double[][] xS = new double[nSymm][nAtoms];
    double[][] yS = new double[nSymm][nAtoms];
    double[][] zS = new double[nSymm][nAtoms];
    // Cartesian center of each molecule
    double[][] centerMolsCart = new double[nSymm][3];

    // Loop over replicate crystal SymOps
    List<SymOp> inflatedSymOps = replicatesCrystal.spaceGroup.symOps;
    for (int iSym = 0; iSym < nSymm; iSym++) {
      SymOp symOp = inflatedSymOps.get(iSym);
      // Apply SymOp to the asymmetric unit atoms Cartesian Coordinates.
      replicatesCrystal.applySymOp(nAtoms, x, y, z, xS[iSym], yS[iSym], zS[iSym], symOp);
      // Compute center-of-mass (CoM) for Cartesian coordinates
      double[] centerOfMass = new double[3];
      int index = 0;
      if (original) {
        // PACCOM conformity -- original algorithm did not use masses.
        for (int i = 0; i < nAtoms; i++) {
          centerOfMass[0] += xS[iSym][index];
          centerOfMass[1] += yS[iSym][index];
          centerOfMass[2] += zS[iSym][index++];
        }
        centerOfMass[0] /= nAtoms;
        centerOfMass[1] /= nAtoms;
        centerOfMass[2] /= nAtoms;
      } else {
        double totalMass = 0.0;
        for (Atom atom : atoms) {
          double m = atom.getMass();
          centerOfMass[0] += xS[iSym][index] * m;
          centerOfMass[1] += yS[iSym][index] * m;
          centerOfMass[2] += zS[iSym][index++] * m;
          totalMass += m;
        }
        centerOfMass[0] /= totalMass;
        centerOfMass[1] /= totalMass;
        centerOfMass[2] /= totalMass;
      }

      // Save CoM cartesian coordinates
      centerMolsCart[iSym] = centerOfMass;
    }

    // Save (mark) a molecule as being closest to the center of the replicates crystal (0.5, 0.5, 0.5)
    // Convert (0.5, 0.5, 0.5) to Cartesian Coordinates
    double[] fracCenter = {0.5, 0.5, 0.5};
    double[] cartCenter = new double[3];
    replicatesCrystal.toCartesianCoordinates(fracCenter, cartCenter);

    DoubleIndexPair[] moleculesDistances = new DoubleIndexPair[nSymm];
    for (int i = 0; i < nSymm; i++) {
      // Then compute Euclidean distance from Cartesian center of the replicates cell
      moleculesDistances[i] = new DoubleIndexPair(i, dist(cartCenter, centerMolsCart[i]));
    }

    // Sort the molecules by their distance from the center.
    // Note that the smallest distances are first in the array after the sort.
    Arrays.sort(moleculesDistances);

    // Create expanded assembly
    MolecularAssembly expandedAssembly = new MolecularAssembly(molecularAssembly.getName());
    List<Bond> bondList = molecularAssembly.getBondList();
    ArrayList<Atom> newAtomList = new ArrayList<>();
    int atomIndex = 0;
    logger.info("\n Copy  SymOp        Distance");
    for (int n = 0; n < nMolecules; n++) {
      // Current molecule
      int iSym = moleculesDistances[n].getIndex();
      double distance = moleculesDistances[n].getDoubleValue();
      logger.info(format(" %4d  %5d  %16.8f", n, iSym, distance));

      symOps.add(inflatedSymOps.get(iSym));

      ArrayList<Atom> atomList = new ArrayList<>();
      // Create a new set of Atoms for each SymOp molecule
      for (int i = 0; i < nAtoms; i++) {
        Atom a = atoms[i];
        double[] xyz = new double[3];
        xyz[0] = xS[iSym][i];
        xyz[1] = yS[iSym][i];
        xyz[2] = zS[iSym][i];
        Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
        atomList.add(atom);
      }
      // Create a new set of Bonds for each SymOp molecule
      for (Bond bond : bondList) {
        Atom a1 = bond.getAtom(0);
        Atom a2 = bond.getAtom(1);
        Atom newA1 = atomList.get(a1.getIndex() - 1);
        Atom newA2 = atomList.get(a2.getIndex() - 1);
        Bond b = new Bond(newA1, newA2);
        b.setBondType(bond.getBondType());
      }
      newAtomList.addAll(atomList);
    }

    // Construct the force field for the expanded set of molecules
    ForceField forceField = molecularAssembly.getForceField();

    // Clear all periodic boundary keywords to achieve an aperiodic system.
    forceField.clearProperty("a-axis");
    forceField.clearProperty("b-axis");
    forceField.clearProperty("c-axis");
    forceField.clearProperty("alpha");
    forceField.clearProperty("beta");
    forceField.clearProperty("gamma");
    forceField.clearProperty("spacegroup");

    expandedAssembly.setForceField(forceField);

    // The biochemistry method is designed to load chemical entities into the
    // Polymer, Molecule, Water and Ion data structure.
    Utilities.biochemistry(expandedAssembly, newAtomList);
    expandedAssembly.finalize(true, forceField);
    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(expandedAssembly);
    expandedAssembly.setPotential(energy);

    expandedAssembly.setFile(molecularAssembly.getFile());

    return expandedAssembly;
  }

  /**
   * Generate new coordinates based on SymOps from system 1
   *
   * @param baseCrystal ReplicatesCrystal for System 1.
   * @param molecularAssembly Base molecular assembly to expand.
   * @param radius Approximate radius of the inflated system.
   * @return Coordinates of expanded system 2.
   */
  public static double[] mimicSphere(Crystal baseCrystal,
      MolecularAssembly molecularAssembly, List<SymOp> symOps, double radius) {

    // Move molecules into the unit cell.
    molecularAssembly.moveAllIntoUnitCell();

    // Collect asymmetric unit atomic coordinates.
    Atom[] atoms = molecularAssembly.getAtomArray();
    int nAtoms = atoms.length;
    double[] x = new double[nAtoms];
    double[] y = new double[nAtoms];
    double[] z = new double[nAtoms];
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      x[i] = atom.getX();
      y[i] = atom.getY();
      z[i] = atom.getZ();
    }

    // Create a ReplicatesCrystal with the unit cell and space group of System 2,
    // but the replicated size [l,m,n] of System 1.
    int l = 1;
    int m = 1;
    int n = 1;
    if (baseCrystal instanceof ReplicatesCrystal) {
      ReplicatesCrystal baseReplicatesCrystal = (ReplicatesCrystal) baseCrystal;
      l = baseReplicatesCrystal.getL();
      m = baseReplicatesCrystal.getM();
      n = baseReplicatesCrystal.getN();
    }
    Crystal unitCell = molecularAssembly.getCrystal().getUnitCell();
    Crystal replicatesCrystal = new ReplicatesCrystal(unitCell, l, m, n, radius * 2.0);
    logger.fine(" Mimic replicates crystal " + replicatesCrystal.toString());

    // Symmetry coordinates
    int nSymm = symOps.size();
    double[] xS = new double[nAtoms];
    double[] yS = new double[nAtoms];
    double[] zS = new double[nAtoms];
    // Output coordinates
    double[] xM = new double[3 * nSymm * nAtoms];

    // Loop over replicate crystal SymOps
    for (int iSym = 0; iSym < nSymm; iSym++) {
      SymOp symOp = symOps.get(iSym);
      // Apply SymOp to the asymmetric unit atoms Cartesian Coordinates.
      replicatesCrystal.applySymOp(nAtoms, x, y, z, xS, yS, zS, symOp);

      // Load the 1D return array.
      int offset = iSym * 3 * nAtoms;
      for (int i = 0; i < nAtoms; i++) {
        int index = offset + i * 3;
        xM[index] = xS[i];
        xM[index + 1] = yS[i];
        xM[index + 2] = zS[i];
      }
    }

    return xM;
  }

//  /**
//   * Cut a molecular sub-sphere from a larger sphere (inflatedAssembly) that matches a smaller
//   * molecular sphere (baseAssembly).
//   *
//   * @param asymmetricUnitAssembly Assembly of the original molecule.
//   * @param inflatedAssembly Assembly of the current sphere.
//   * @param similarMolecules Indices for desired molecules.
//   * @return MolecularAssembly cut from molecularAssembly that is similar to baseAssembly.
//   */
//  public static MolecularAssembly cutSubSphere(MolecularAssembly asymmetricUnitAssembly,
//      MolecularAssembly inflatedAssembly, DoubleIndexPair[] similarMolecules, int nMolecules) {
//
//    int nAtoms = asymmetricUnitAssembly.getAtomArray().length;
//    int nInflatedMolecules = inflatedAssembly.getMolecules().size();
//
//    Atom[] inflatedAtoms = inflatedAssembly.getAtomArray();
//    double[][] xInflated = new double[nInflatedMolecules][nAtoms];
//    double[][] yInflated = new double[nInflatedMolecules][nAtoms];
//    double[][] zInflated = new double[nInflatedMolecules][nAtoms];
//    // This assumes the atoms of each molecule are contiguous. This is not enforced for XYZ files.
//    for (int i = 0; i < nInflatedMolecules; i++) {
//      for (int j = 0; j < nAtoms; j++) {
//        Atom atom = inflatedAtoms[nAtoms * i + j];
//        xInflated[i][j] = atom.getX();
//        yInflated[i][j] = atom.getY();
//        zInflated[i][j] = atom.getZ();
//      }
//    }
//
//    // Create assembly for second system based on optimal nMolecule matches
//    MolecularAssembly expandedAssembly = new MolecularAssembly(asymmetricUnitAssembly.getName());
//    List<Bond> bondList = asymmetricUnitAssembly.getBondList();
//    ArrayList<Atom> newAtomList = new ArrayList<>();
//    int atomIndex = 0;
//    for (int n = 0; n < nMolecules; n++) {
//      DoubleIndexPair molecule = similarMolecules[n];
//      ArrayList<Atom> atomList = new ArrayList<>();
//      // Create a new set of Atoms for each SymOp molecule
//      int index = molecule.getIndex();
//      int atomOffset = index * nAtoms;
//      for (int i = 0; i < nAtoms; i++) {
//        Atom a = inflatedAtoms[atomOffset + i];
//        double[] xyz = new double[3];
//        xyz[0] = xInflated[index][i];
//        xyz[1] = yInflated[index][i];
//        xyz[2] = zInflated[index][i];
//        Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz);
//        atomList.add(atom);
//      }
//      // Create a new set of Bonds for each SymOp molecule
//      for (Bond bond : bondList) {
//        Atom a1 = bond.getAtom(0);
//        Atom a2 = bond.getAtom(1);
//        Atom newA1 = atomList.get(a1.getIndex() - 1);
//        Atom newA2 = atomList.get(a2.getIndex() - 1);
//        Bond b = new Bond(newA1, newA2);
//        b.setBondType(bond.getBondType());
//      }
//      newAtomList.addAll(atomList);
//    }
//
//    // Construct the force field for the expanded set of molecules
//    ForceField forceField = inflatedAssembly.getForceField();
//
//    // Clear all periodic boundary keywords to achieve an aperiodic system.
//    forceField.clearProperty("a-axis");
//    forceField.clearProperty("b-axis");
//    forceField.clearProperty("c-axis");
//    forceField.clearProperty("alpha");
//    forceField.clearProperty("beta");
//    forceField.clearProperty("gamma");
//    forceField.clearProperty("spacegroup");
//
//    expandedAssembly.setForceField(forceField);
//    Utilities.biochemistry(expandedAssembly, newAtomList);
//    expandedAssembly.finalize(true, forceField);
//    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(expandedAssembly);
//    expandedAssembly.setPotential(energy);
//
//    expandedAssembly.setFile(inflatedAssembly.getFile());
//    return expandedAssembly;
//  }

  /**
   * Reorients the given atoms so n2 and n3 are on the X axis and X-Y plane. n1 should be at origin.
   *
   * @param n1 Atom located at the origin
   * @param n2 Atom to go to X axis
   * @param n3 Atom to go to X-Y plane
   * @return double[] New coordinates for the atoms
   * @author Okimasa OKADA
   * @author Aaron J. Nessler
   */
  static double[] standardOrientation(Atom n1, Atom n2, Atom n3) {
    double[] orientedCoords = new double[9];
    double p1n2 = n2.getX();
    double q1n2 = n2.getY();
    double r1n2 = n2.getZ();
    double p2n2;
    double q2n2;

    double p1n3 = n3.getX();
    double q1n3 = n3.getY();
    double r1n3 = n3.getZ();
    double p2n3;
    double q2n3;

    logger.finer(format("START: N1: %16.8f %16.8f %16.8f", n1.getX(), n1.getY(), n1.getZ()));
    logger.finer(format("START: N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
    logger.finer(format("START: N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));

    // Calculation of sigma, phai, and cita angles needed to get specified atoms to desired loci
    double cita0 = acos(p1n2 / sqrt(pow(p1n2, 2) + pow(q1n2, 2)));
    double phai0 = acos(
        sqrt(pow(p1n2, 2) + pow(q1n2, 2)) / sqrt(pow(p1n2, 2) + pow(q1n2, 2) + pow(r1n2, 2)));
    if (q1n2 < 0.0) {
      cita0 = -cita0;
    }

    logger.finer(format("cita: %16.8f", cita0));
    logger.finer(format("cos(cita): %16.8f", cos(cita0)));
    logger.finer(format("sin(cita): %16.8f", sin(cita0)));

    p2n2 = p1n2 * cos(cita0) + q1n2 * sin(cita0);
    q1n2 = -p1n2 * sin(cita0) + q1n2 * cos(cita0);
    p2n3 = p1n3 * cos(cita0) + q1n3 * sin(cita0);
    q1n3 = -p1n3 * sin(cita0) + q1n3 * cos(cita0);

    p1n2 = p2n2;
    p1n3 = p2n3;

    logger.finer(format("Step 1 N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
    logger.finer(format("Step 1 N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));

    if (r1n2 > 0.0) {
      phai0 = -phai0;
    }

    logger.finer(format("phai: %16.8f", phai0));

    p2n2 = p1n2 * cos(phai0) - r1n2 * sin(phai0);
    r1n2 = p1n2 * sin(phai0) + r1n2 * cos(phai0);
    p2n3 = p1n3 * cos(phai0) - r1n3 * sin(phai0);
    r1n3 = p1n3 * sin(phai0) + r1n3 * cos(phai0);

    p1n2 = p2n2;
    p1n3 = p2n3;

    logger.finer(format("Step 2 N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
    logger.finer(format("Step 2 N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));

    double sigma0 = acos(q1n3 / sqrt(q1n3 * q1n3 + r1n3 * r1n3));
    if (r1n3 < 0.0) {
      sigma0 = -sigma0;
    }

    logger.finer(format("sigma: %16.8f", sigma0));

    q2n2 = q1n2 * cos(sigma0) + r1n2 * sin(sigma0);
    r1n2 = -q1n2 * sin(sigma0) + r1n2 * cos(sigma0);
    q2n3 = q1n3 * cos(sigma0) + r1n3 * sin(sigma0);
    r1n3 = -q1n3 * sin(sigma0) + r1n3 * cos(sigma0);

    q1n2 = q2n2;
    q1n3 = q2n3;

    logger.finer(format("DONE N1: %16.8f %16.8f %16.8f", n1.getX(), n1.getY(), n1.getZ()));
    logger.finer(format("DONE N2: %16.8f %16.8f %16.8f", p1n2, q1n2, r1n2));
    logger.finer(format("DONE N3: %16.8f %16.8f %16.8f", p1n3, q1n3, r1n3));

    orientedCoords[0] = n1.getX();
    orientedCoords[1] = n1.getY();
    orientedCoords[2] = n1.getZ();
    orientedCoords[3] = p1n2;
    orientedCoords[4] = q1n2;
    orientedCoords[5] = r1n2;
    orientedCoords[6] = p1n3;
    orientedCoords[7] = q1n3;
    orientedCoords[8] = r1n3;

    return orientedCoords;
  }

  /**
   * Save the current assembly as a PDB file.
   *
   * @param currentAssembly Current assembly of which to save a PDB.
   * @param description Prefix for the file being saved.
   */
  public static void saveAssemblyPDB(String fileName, MolecularAssembly currentAssembly,
      String description) {
    File saveLocationPDB = new File(fileName + description + ".pdb");
    // Save aperiodic system of n_mol closest atoms for visualization
    logger.info(" Saving " + saveLocationPDB.getName());
    PDBFilter pdbfilter = new PDBFilter(saveLocationPDB, currentAssembly, null, null);
    pdbfilter.writeFile(saveLocationPDB, false);
  }
}