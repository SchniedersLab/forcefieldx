// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.potential.bonded;

import static java.lang.String.format;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The ResidueState class encodes the current chemical and coordinate state of a Residue,
 * particularly a MultiResidue, for ease of reverting coordinates. Should not be too difficult to get
 * it to also store velocities, etc.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class ResidueState {

  private static final Logger logger = Logger.getLogger(ResidueState.class.getName());

  /**
   * For a MultiResidue, parent is the MultiResidue, while res is the chemical state to which this
   * ResidueState belongs. For a regular Residue, both are just the Residue. In general, use res to
   * store a chemical or other state, and parent to store the Residue to which this belongs.
   */
  private final Residue parent;

  private final Residue residue;
  private final HashMap<Atom, double[]> atomMap;
  private final Atom[] atoms;

  /**
   * Constructor for ResidueState.
   *
   * @param parent a {@link ffx.potential.bonded.Residue} object.
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   */
  public ResidueState(Residue parent, Residue residue) {
    this.parent = parent;
    this.residue = residue;

    List<Atom> atomList = residue.getAtomList();
    int nAtoms = atomList.size();
    atomMap = new HashMap<>(nAtoms);
    atoms = new Atom[nAtoms];

    atomList.toArray(atoms);
    for (Atom atom : atoms) {
      double[] atXYZ = new double[3];
      atom.getXYZ(atXYZ);
      atomMap.put(atom, atXYZ);
    }
  }

  /**
   * Constructor for ResidueState.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   */
  public ResidueState(Residue residue) {
    this(residue, residue);
    if (residue instanceof MultiResidue) {
      throw new IllegalArgumentException(format(
          " Residue %s is a MultiResidue: this ResidueState has been incorrectly constructed!",
          residue));
    }
  }

  /**
   * revertAllCoordinates.
   *
   * @param residueList a {@link java.util.List} object.
   * @param states an array of {@link ffx.potential.bonded.ResidueState} objects.
   */
  public static void revertAllCoordinates(List<Residue> residueList, ResidueState[] states) {
    revertAllCoordinates(residueList.toArray(new Residue[0]), states);
  }

  /**
   * Uses a double[nAtoms][3] to revert the coordinates of an array of atoms. Does not check to see
   * if the arrays are properly ordered; make sure you are using the correct arrays.
   *
   * @param atoms To revert coordinates of
   * @param coords Coordinates
   */
  public static void revertAtomicCoordinates(Atom[] atoms, double[][] coords) {
    int nAtoms = atoms.length;
    if (coords.length != nAtoms) {
      throw new IllegalArgumentException(
          format(" Length %d of atoms " + "array does not match length %d of coordinates array",
              nAtoms, coords.length));
    }
    for (int i = 0; i < nAtoms; i++) {
      atoms[i].setXYZ(coords[i]);
    }
  }

  /**
   * storeAllCoordinates.
   *
   * @param residueList a {@link java.util.List} object.
   * @return an array of {@link ffx.potential.bonded.ResidueState} objects.
   */
  public static ResidueState[] storeAllCoordinates(List<Residue> residueList) {
    return storeAllCoordinates(residueList.toArray(new Residue[0]));
  }

  /**
   * storeAllCoordinates.
   *
   * @param residues an array of {@link ffx.potential.bonded.Residue} objects.
   * @return an array of {@link ffx.potential.bonded.ResidueState} objects.
   */
  public static ResidueState[] storeAllCoordinates(Residue[] residues) {
    int nResidues = residues.length;
    ResidueState[] states = new ResidueState[nResidues];
    for (int i = 0; i < nResidues; i++) {
      states[i] = residues[i].storeState();
    }
    return states;
  }

  /**
   * Returns a new double[nAtoms][3] with the coordinates of an array of atoms.
   *
   * @param atoms To store coordinates of
   * @return Coordinates
   */
  public static double[][] storeAtomicCoordinates(Atom[] atoms) {
    int nAtoms = atoms.length;
    double[][] coords = new double[nAtoms][3];
    for (int i = 0; i < nAtoms; i++) {
      atoms[i].getXYZ(coords[i]);
    }
    return coords;
  }

  /**
   * revertAllCoordinates.
   *
   * @param residueArray an array of {@link ffx.potential.bonded.Residue} objects.
   * @param states an array of {@link ffx.potential.bonded.ResidueState} objects.
   */
  private static void revertAllCoordinates(Residue[] residueArray, ResidueState[] states) {
    int nResidues = residueArray.length;
    if (nResidues != states.length) {
      throw new IllegalArgumentException(
          format("Length of residue " + "array %d and residue state array %d do not match.",
              nResidues, states.length));
    }

    for (int i = 0; i < nResidues; i++) {
      Residue resi = residueArray[i];
      // If indexing did not get thrown off:
      if (resi.equals(states[i].getParent())) {
        resi.revertState(states[i]);
      } else {
        // Else must search states[]
        boolean matchFound = false;
        for (int j = 0; j < nResidues; j++) {
          if (resi.equals(states[j].getParent())) {
            matchFound = true;
            resi.revertState(states[j]);
            break;
          }
        }

        if (!matchFound) {
          throw new IllegalArgumentException(
              format("Could not " + "find match for residue %s among residue states array.", resi));
        }
      }
    }
  }

  public double compareTo(ResidueState residueState) {
    logger.info("Comparing rotamers using the compareTo method in ResidueState");
    double[][] tempx1 = new double[this.atoms.length][3];
    double[][] tempx2 = new double[residueState.atoms.length][3];
    double[] x1 = new double[this.atoms.length * 3];
    double[] x2 = new double[residueState.atoms.length * 3];
    // Create arrays of atoms for keptRotamer[k] ResidueState and newResState
    // Here, each residue state is the same residue with a different rotamer applied.
    Atom[] x1atoms = this.atoms;
    Atom[] x2atoms = residueState.atoms;
    // logger.info("Total number of atoms for current residue state: " + x1atoms.length);
    // logger.info("Total number of atoms for previously kept residue state to be tested against: "+
    // x2atoms.length);

    // Get the coordinate arrays for each atom in the current residue state using the atomMap
    for (int atomCount = 0; atomCount < x1atoms.length; atomCount++) {
      Atom exampleAtom = x1atoms[atomCount];
      tempx1[atomCount] = this.atomMap.get(exampleAtom);
    }

    // Add coordinates from tempx1[][] to x1[] for RMSD calculations
    int x1count = 0;
    for (double[] coordSet : tempx1) {
      for (int coordCount = 0; coordCount < 3; coordCount++) {
        x1[x1count] = coordSet[coordCount];
        x1count++;
      }
    }

    // Get the coorinate arrays for each atom in the previously kept residue state to be
    // tested against using the atomMap for that residue state
    for (int atomCount = 0; atomCount < x2atoms.length; atomCount++) {
      Atom exampleAtom = x2atoms[atomCount];
      tempx2[atomCount] = residueState.atomMap.get(exampleAtom);
    }

    // Add coordinate from tempx2[][] to x2[] for RMSD calculations
    int x2count = 0;
    for (double[] coordSet : tempx2) {
      for (int coordCount = 0; coordCount < 3; coordCount++) {
        x2[x2count] = coordSet[coordCount];
        x2count++;
      }
    }

    // Create a mass array of ones so the RMSD isn't mass weighted
    double[] mass = new double[x1.length];
    Arrays.fill(mass, 1);
    return ffx.potential.utils.Superpose.rmsd(x1, x2, mass);
  }

  /**
   * Getter for the field <code>atoms</code>.
   *
   * @return an array of {@link ffx.potential.bonded.Atom} objects.
   */
  public Atom[] getAtoms() {
    return atoms;
  }

  /**
   * Getter for the field <code>parent</code>.
   *
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  public Residue getParent() {
    return parent;
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return format(" ResidueState with parent residue %s, state residue " + "%s, number of atoms %d",
        parent, residue, atoms.length);
  }

  /**
   * getStateResidue.
   *
   * @return a {@link ffx.potential.bonded.Residue} object.
   */
  Residue getStateResidue() {
    return residue;
  }

  /**
   * Returns the coordinates of the selected atom; or null if atom not contained.
   *
   * @param atom An Atom.
   * @return xyz coordinates.
   */
  double[] getAtomCoords(Atom atom) {
    double[] xyz = new double[3];
    if (!atomMap.containsKey(atom)) {
      logger.info(
          format(" Illegal call to ResidueState.getAtomCoords: atom %s not found: hashcode %d", atom,
              atom.hashCode()));
      for (Atom ratom : residue.getAtomList()) {
        logger.info(format(" Atoms in residue: %s hashcode: %d", ratom, ratom.hashCode()));
      }
      for (Atom matom : atomMap.keySet()) {
        logger.info(
            format(" Atoms in ResidueState atom cache: %s hashcode: %d", matom, matom.hashCode()));
      }
      logger.log(Level.SEVERE, " Error in ResidueState.getAtomCoords.", new IllegalStateException());
    }
    System.arraycopy(atomMap.get(atom), 0, xyz, 0, 3);
    return xyz;
  }
}
