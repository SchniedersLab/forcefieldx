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

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.cbrt;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SymOp;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSNode;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.FFXScript;
import org.apache.commons.math3.ml.distance.EuclideanDistance;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.logging.Logger;

/**
 * PACCOMFunctions contains methods utilized in the PACCOM.groovy file.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 */
public class PACCOMFunctions {

  /** The logger for this class. */
  protected static final Logger logger = Logger.getLogger(FFXScript.class.getName());

  /** Distance calculator from Apache Commons. */
  protected static final EuclideanDistance distance = new EuclideanDistance();
  /**
   * Generate a molecular assembly of a sphere containing a specified number of molecules.
   *
   * @param molecularAssembly Base molecular assembly to expand.
   * @param nMolecules Number of molecules to include in sphere.
   * @param original If following the original version of PACCOM (by Okimasa Okada)
   * @return MolecularAssembly of the desired sphere size.
   */
  public static MolecularAssembly generateBaseSphere(MolecularAssembly molecularAssembly,
                                                     int nMolecules,
                                                     boolean original) {

    // Move molecules into the unit cell.
    molecularAssembly.moveAllIntoUnitCell();

    Crystal crystal = molecularAssembly.getCrystal();
    double asymmetricUnitVolume = crystal.volume / crystal.getNumSymOps();

    // Add wiggle room for boundary cutoffs
    double inflationFactor = 2.0;

    // Estimate a radius that will include "nMolecules".
    double radius = cbrt((3.0 / (4.0 * PI) * nMolecules * asymmetricUnitVolume)) + inflationFactor;

    // Atoms in the asymmetric unit.
    Atom[] atoms = molecularAssembly.getAtomArray();
    int nAtoms = atoms.length;

    Crystal replicatesCrystal =
        ReplicatesCrystal.replicatesCrystalFactory(molecularAssembly.getCrystal(), radius * 2.0);
    molecularAssembly.setCrystal(replicatesCrystal);

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

    int nSymm = replicatesCrystal.spaceGroup.symOps.size();
    // Symmetry coordinates for each molecule in replicates crystal
    double[][] xS = new double[nSymm][nAtoms];
    double[][] yS = new double[nSymm][nAtoms];
    double[][] zS = new double[nSymm][nAtoms];
    double[][] centerMolsFrac = new double[nSymm][3]; //Fractional center of each molecule
    double[][] centerMolsCart = new double[nSymm][3]; //Cartesian center of each molecule

    // Loop over replicate crystal SymOps
    for (int iSym = 0; iSym < nSymm; iSym++) {
      SymOp symOp = replicatesCrystal.spaceGroup.getSymOp(iSym);
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

      // Convert CoM to fractional coordinates
      replicatesCrystal.toFractionalCoordinates(centerOfMass, centerOfMass);

      // Save CoM fractional coordinates
      centerMolsFrac[iSym] = centerOfMass;
    }

    int[] desiredMols = new int[nMolecules];

    // Save (mark) a molecule as being closest to the center of the replicates crystal (0.5, 0.5, 0.5)
    // Convert (0.5, 0.5, 0.5) to Cartesian Coordinates
    double[] fracCenter = {0.5, 0.5, 0.5};
    double[] cartCenter = new double[3];
    replicatesCrystal.toCartesianCoordinates(fracCenter, cartCenter);

    // Molecule index closest to center of crystal.
    int minIndex = -1;
    double minDist = Double.MAX_VALUE;
    for (int i = 0; i < nSymm; i++) {
      // Convert each CoM to cartesian coordinates
      replicatesCrystal.toCartesianCoordinates(centerMolsFrac[i], centerMolsCart[i]);

      // Then compute euclidean distance from cartesian center of the replicates cell
      double value = distance.compute(cartCenter, centerMolsCart[i]);

      // Save index of molecule closest to replicates cell center
      if (value < minDist) {
        minDist = value;
        minIndex = i;
      }
    }

    // Hashmap with index and distance
    HashMap<Integer, Double> closestMols = new HashMap<>();

    // Loop over saved CoM coordinates to find the nMolecules closest to the marked molecule
    for (int i = 0; i < nSymm; i++) {
      double value = distance.compute(centerMolsCart[minIndex], centerMolsCart[i]);

      if (closestMols.size() < nMolecules) {
        closestMols.put(i, value);
      } else {
        double maxValue = Collections.max(closestMols.values());
        if (value < maxValue) {
          // Remove most distant molecule
          Iterator<Entry<Integer, Double>> it = closestMols.entrySet().iterator();
          while (it.hasNext()) {
            Entry<Integer, Double> item = it.next();
            // Hopefully precision issues do not cause any issues with a direct comparison of double values.
            if (item.getValue() == maxValue) {
              it.remove();
              break;
            }
          }
          closestMols.put(i, value);
        }
      }
    }

    // Final collection of nMolcules closest to the marked molecule.
    int index = 0;
    for (Integer key : closestMols.keySet()) {
      desiredMols[index++] = key;
    }

    // Create expanded assemblies containing nMolecules
    MolecularAssembly expandedAssembly = new MolecularAssembly(molecularAssembly.getName());
    List<Bond> bondList = molecularAssembly.getBondList();
    ArrayList<Atom> newAtomList = new ArrayList<>();
    int atomIndex = 0;
    for (int iSym : desiredMols) {
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
    Utilities.biochemistry(expandedAssembly, newAtomList);
    expandedAssembly.finalize(true, forceField);
    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(expandedAssembly);
    expandedAssembly.setPotential(energy);

    // This energy evaluation could be removed in the future. Is it useful as a consistency check?
    energy.energy(true, true);

    expandedAssembly.setFile(molecularAssembly.getFile());

    return expandedAssembly;
  }

  /**
   * Cut a molecular sub-sphere from a larger sphere (molecularAssembly) that matches a smaller molecular sphere
   * (baseAssembly).
   *
   * @param asymmetricUnitAssembly Assembly of the original molecule.
   * @param targetAssembly Assembly containing the desired sphere size/orientation.
   * @param inflatedAssembly Assembly of the current sphere.
   * @param original If following the original version of PACCOM (by Okimasa Okada)
   * @return MolecularAssembly cut from molecularAssembly that is similar to baseAssembly.
   */
  public static MolecularAssembly cutSubSphere(MolecularAssembly asymmetricUnitAssembly,
                                               MolecularAssembly targetAssembly,
                                               MolecularAssembly inflatedAssembly, boolean original) {

    int nAtoms = asymmetricUnitAssembly.getAtomArray().length;
    int nMolecules = targetAssembly.getMolecules().size();
    int nInflatedMolecules = inflatedAssembly.getMolecules().size();

    Atom[] inflatedAtoms = inflatedAssembly.getAtomArray();
    double[][] xInflated = new double[nInflatedMolecules][nAtoms];
    double[][] yInflated = new double[nInflatedMolecules][nAtoms];
    double[][] zInflated = new double[nInflatedMolecules][nAtoms];
    // This assumes the atoms of each molecule are contiguous. This is not enforced for XYZ files.
    for (int i = 0; i < nInflatedMolecules; i++) {
      for (int j = 0; j < nAtoms; j++) {
        Atom atom = inflatedAtoms[nAtoms * i + j];
        xInflated[i][j] = atom.getX();
        yInflated[i][j] = atom.getY();
        zInflated[i][j] = atom.getZ();
      }
    }

    // Compute center-of-mass (CoM) for Cartesian coordinates
    double[][] centerOfMols = new double[nInflatedMolecules][3];

    if (original) {
      // PACCOM conformity. Original PACCOM did not use masses.
      for (int molIndex = 0; molIndex < nInflatedMolecules; molIndex++) {
        double[] centerOfMass = {0.0, 0.0, 0.0};
        for (int atomIndex = 0; atomIndex < nAtoms; atomIndex++) {
          centerOfMass[0] += xInflated[molIndex][atomIndex];
          centerOfMass[1] += yInflated[molIndex][atomIndex];
          centerOfMass[2] += zInflated[molIndex][atomIndex];
        }
        centerOfMass[0] /= inflatedAtoms.length;
        centerOfMass[1] /= inflatedAtoms.length;
        centerOfMass[2] /= inflatedAtoms.length;
        centerOfMols[molIndex] = centerOfMass;
      }
    } else {
      for (int molIndex = 0; molIndex < nInflatedMolecules; molIndex++) {
        double[] centerOfMass = {0.0, 0.0, 0.0};
        double totalMass = 0.0;
        for (int atomIndex = 0; atomIndex < nAtoms; atomIndex++) {
          double m = inflatedAtoms[molIndex * nAtoms + atomIndex].getMass();
          centerOfMass[0] += xInflated[molIndex][atomIndex] * m;
          centerOfMass[1] += yInflated[molIndex][atomIndex] * m;
          centerOfMass[2] += zInflated[molIndex][atomIndex] * m;
          totalMass += m;
        }
        centerOfMass[0] /= totalMass;
        centerOfMass[1] /= totalMass;
        centerOfMass[2] /= totalMass;
        centerOfMols[molIndex] = centerOfMass;
      }
    }

    int[] desiredMols = new int[nMolecules];
    Arrays.fill(desiredMols, -1);
    int index = 0;
    for (MSNode baseMolecule : targetAssembly.getMolecules()) {
      double[] baseCenter = baseMolecule.getCenter(true);
      int key = -1;
      double minValue = Double.MAX_VALUE;
      for (int i = 0; i < nInflatedMolecules; i++) {
        boolean duplicate = false;
        double value = distance.compute(baseCenter, centerOfMols[i]);
        if (value < minValue) {
          // Avoid using a molecule from the inflated system twice.
          for (int prevKey : desiredMols) {
            if (prevKey == i) {
              // Perhaps log a warning.
              duplicate = true;
              break;
            }
          }
          if (!duplicate) {
            minValue = value;
            key = i;
          }
        }
      }
      desiredMols[index++] = key;
    }

    // Create assembly for second system based on optimal nMolecule matches
    MolecularAssembly expandedAssembly = new MolecularAssembly(asymmetricUnitAssembly.getName());
    List<Bond> bondList = asymmetricUnitAssembly.getBondList();
    ArrayList<Atom> newAtomList = new ArrayList<>();
    int atomIndex = 0;
    for (int molecule : desiredMols) {
      ArrayList<Atom> atomList = new ArrayList<>();
      // Create a new set of Atoms for each SymOp molecule
      for (int i = 0; i < nAtoms; i++) {
        Atom a = inflatedAtoms[molecule * nAtoms + i];
        double[] xyz = new double[3];
        xyz[0] = xInflated[molecule][i];
        xyz[1] = yInflated[molecule][i];
        xyz[2] = zInflated[molecule][i];
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
    ForceField forceField = inflatedAssembly.getForceField();

    // Clear all periodic boundary keywords to achieve an aperiodic system.
    forceField.clearProperty("a-axis");
    forceField.clearProperty("b-axis");
    forceField.clearProperty("c-axis");
    forceField.clearProperty("alpha");
    forceField.clearProperty("beta");
    forceField.clearProperty("gamma");
    forceField.clearProperty("spacegroup");

    expandedAssembly.setForceField(forceField);
    Utilities.biochemistry(expandedAssembly, newAtomList);
    expandedAssembly.finalize(true, forceField);
    ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(expandedAssembly);
    expandedAssembly.setPotential(energy);

    // This energy evaluation could be removed in the future. Is it useful as a consistency check?
    energy.energy(true, true);

    expandedAssembly.setFile(inflatedAssembly.getFile());
    return expandedAssembly;
  }

  /**
   * Determine the index for the central molecule based on center of mass for an assembly.
   *
   * @param assembly Assembly of which the center is desired.
   * @return int index of the molecule closest to the center of mass.
   */
  public static int centerMoleculeIndex(MolecularAssembly assembly) {
    int minIndex = -1; // Molecule index closest to center of crystal.
    double minDist = Double.MAX_VALUE;
    int nMolecules = assembly.getMolecules().size();
    double[] centerOfMass = new double[3];
    double totalMass = 0.0;

    for (Atom atom : assembly.getAtomArray()) {
      double m = atom.getMass();
      centerOfMass[0] += atom.getX() * m;
      centerOfMass[1] += atom.getY() * m;
      centerOfMass[2] += atom.getZ() * m;
      totalMass += m;
    }
    centerOfMass[0] /= totalMass;
    centerOfMass[1] /= totalMass;
    centerOfMass[2] /= totalMass;

    assembly.getMolecules().get(0).getCenter(true);
    for (int i = 0; i < nMolecules; i++) {
      double[] centerMol = assembly.getMolecules().get(i).getCenter(true);
      double value = distance.compute(centerOfMass, centerMol);
      if (value < minDist) {
        minDist = value;
        minIndex = i;
      }
    }
    return minIndex;
  }

  /**
   * Sort the given HashMap (Integer, Double) by values (Double).
   *
   * @param hm Hashmap to be sorted.
   * @return HashMap(Integer, Double) that has been sorted.
   */
  public static HashMap<Integer, Double> sortHashMapByValue(HashMap<Integer, Double> hm) {
    // Create a list from elements of HashMap
    List<Entry<Integer, Double>> list = new LinkedList<>(hm.entrySet());

    // Sort the list
    list.sort(Entry.comparingByValue());

    // put data from sorted list to hashMap
    HashMap<Integer, Double> temp = new LinkedHashMap<>();
    for (Entry<Integer, Double> aa : list) {
      temp.put(aa.getKey(), aa.getValue());
    }
    return temp;
  }

  /**
   * Save the current assembly as a PDB file.
   *
   * @param currentAssembly Current assembly of which to save a PDB.
   * @param saveLocationPDB A pre-made file to save the PDB.
   */
  public static void saveAssemblyPDB(MolecularAssembly currentAssembly,
      File saveLocationPDB) {
    //Save aperiodic system of n_mol closest atoms for visualization
    logger.info(" Saving " + saveLocationPDB.getName());
    PDBFilter pdbfilter = new PDBFilter(saveLocationPDB, currentAssembly, null, null);
    pdbfilter.writeFile(saveLocationPDB, false);
  }
}
