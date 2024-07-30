//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy.test

import ffx.potential.MolecularAssembly
import ffx.potential.Utilities
import ffx.potential.bonded.Angle
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.bonded.MSGroup
import ffx.potential.bonded.MSNode
import ffx.potential.bonded.Molecule
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import ffx.potential.parsers.CIFFilter
import org.apache.commons.configuration2.CompositeConfiguration
import org.openscience.cdk.AtomContainer
import org.openscience.cdk.config.AtomTypeFactory
import org.openscience.cdk.graph.rebond.RebondTool
import org.openscience.cdk.interfaces.IAtom
import org.openscience.cdk.interfaces.IBond
import org.openscience.cdk.isomorphism.AtomMatcher
import org.openscience.cdk.isomorphism.BondMatcher
import org.openscience.cdk.isomorphism.Pattern
import org.openscience.cdk.isomorphism.VentoFoggia
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import javax.vecmath.Point3d
import java.util.logging.Level

import static ffx.numerics.math.DoubleMath.dihedralAngle;
import static ffx.potential.bonded.BondedUtils.intxyz;
import static ffx.potential.utils.Superpose.applyTranslation;
import static ffx.potential.utils.Superpose.calculateTranslation;
import static ffx.potential.utils.Superpose.calculateRotation;
import static ffx.potential.utils.Superpose.applyRotation;
import static ffx.potential.utils.Superpose.rmsd;
import static ffx.potential.utils.Superpose.superpose;
import static java.lang.String.format
import static java.lang.Double.MAX_VALUE
import static org.apache.commons.io.FilenameUtils.getFullPath
import static org.apache.commons.io.FilenameUtils.getName
import static org.apache.commons.io.FilenameUtils.removeExtension
import static org.apache.commons.math3.util.FastMath.sqrt
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getSymbol

/**
 * The ReorderAtoms script sorts the atoms within an XYZ file based on atomic weight.
 *
 * @author Aaron J. Nessler
 * <br>
 * Usage:
 * <br>
 * ffxc test.ReorderAtoms &lt;filename&gt;
 */
@Command(description = " Reorder the atoms of an XYZ file based on atomic weight.",
    name = "test.ReorderAtoms")
class ReorderAtoms extends PotentialScript {

  /**
   * -a or --atomType Change atom types of first file to match second.
   */
  @Option(names = ['-a', '--atomType'], paramLabel = "false", defaultValue = "false",
          description = 'Change atom types of first file to match second.')
  private static boolean atomType

  /**
   * --sw or --swap Switch positions of all equivalent atoms (untested for >2 atoms).
   */
  @Option(names = ['--sw', '--swap'], paramLabel = "false", defaultValue = "false",
      description = 'Switch positions of equivalent atoms (untested for >2 atoms).')
  private static boolean swap

  /**
   * -m or --match Match atom indices based on structure.
   */
  @Option(names = ['-m', '--match'], paramLabel = "false", defaultValue = "false",
          description = 'Match atom indices based on structure.')
  private static boolean match

  /**
   * -r or --reorder Reorder atoms based on environment.
   */
  @Option(names = ['-r', '--reorder'], paramLabel = "false", defaultValue = "false",
          description = 'Reorder atoms based on environment.')
  private static boolean reorder

  /**
   * The final argument(s) should be two or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = 'Atomic coordinate files to reorder in XYZ/ARC format. If 2, change file 1 to match file 2.')
  List<String> filenames = null

  /**
   * Set bond tolerance (distance to determine if bond is made).
   */
  private double bondTolerance = 0.2;

  /**
   * Execute the script.
   */
  @Override
  ReorderAtoms run() {
    logger.warning(" ReorderAtoms script is still in development.\n " +
            "Use at your own risk and please verify output everytime.")

    System.setProperty("vdwterm", "false")
    // READ IN STRUCTURES
    if (!init()) {
      return null
    }
    if (reorder || swap) {
      if(match || atomType){
        // TODO allow any allowable combination of flags at the same time.
        // Atom reordering could be run every single time with the other two as auxiliary.
        logger.warning(" Currently each flag (i.e., atomType, match, or reorder must be run individually. ");
      }
      for (String filename : filenames) {
        potentialFunctions.openAll(filename)
        SystemFilter systemFilter = potentialFunctions.getFilter()
        int numModels = systemFilter.countNumModels()

        // Save new xyz file with the converted crystal.
        File saveLocation = new File(filename)
        File saveFile = potentialFunctions.versionFile(saveLocation)
        for (int i = 0; i < numModels; i++) {
          MolecularAssembly assembly = systemFilter.getActiveMolecularSystem()
          Atom[] atoms = assembly.getAtomArray()
          List<Atom> atomList0 = assembly.getAtomList()
          int nAtoms = atoms.length
          MolecularAssembly newAssembly = new MolecularAssembly(assembly.getName())
          List<Bond> bondList = assembly.getBondList()

          // OLD METHOD: based on atom environment... Breaks down with symmetry.
          // TODO handle sorting/exclusion better...
          int atomListSize = atomList0.size()
          for (int j = 0; j < atomListSize; j++) {
            for (int k = 0; k < atomListSize; k++) {
              if (j != k) {
                if (atomList0[j].getIndex() != atomList0[k].getIndex()) {
                  AtomComparator ac = new AtomComparator()
                  if (logger.isLoggable(Level.FINER)) {
                    logger.finer(" Compare atom: " + atomList0[j].toString() + " and " + atomList0[k].toString())
                  }
                  int value = ac.compare(atomList0[j], atomList0[k])
                  if (value == 0) {
                    logger.info(" This " + atomList0[j].toString() + " and " + atomList0[k].toString() +
                            " are \"equivalent\".")
                    if (swap) {
                      Atom temp = atomList0[j]
                      atomList0[j] = atomList0[k]
                      atomList0[k] = temp
                    }
                  } else if (value < 0 && reorder) {
                    if (logger.isLoggable(Level.FINER)) {
                      logger.finer(" Swapped " + atomList0[j].toString() + " and " + atomList0[k].toString() + ".")
                    }
                    Atom temp = atomList0[j]
                    atomList0[j] = atomList0[k]
                    atomList0[k] = temp
                  }
                }
              }
            }
          }
          int[] originalOrder = new int[nAtoms]
          for (int j = 0; j < nAtoms; j++) {
            originalOrder[j] = atoms[j].getIndex()
            if (logger.isLoggable(Level.FINE)) {
              logger.fine(format(" Original index for atom %3d: %3d New index: %3d", j, originalOrder[j],
                      atomList0[j].getIndex()))
            }
          }
          // Current molecule
          ArrayList<Atom> atomList = new ArrayList<>();
          int atomIndex = 1
          // Create a new set of Atoms for each molecule
          for (int j = 0; j < nAtoms; j++) {
            Atom a = atomList0.get(j)
            double[] xyz = [a.getX(), a.getY(), a.getZ()]
            Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz)
            atomList.add(atom)
          }

          // Create a new set of Bonds for each molecule
          for (Bond bond : bondList) {
            Atom a1 = bond.getAtom(0)
            Atom a2 = bond.getAtom(1)

            int index1 = -1
            int index2 = -1
            int counter = 0
            for (Atom atom : atomList0) {
              int index = atom.getIndex()
              if (index == a1.getIndex()) {
                index1 = counter
              } else if (index == a2.getIndex()) {
                index2 = counter
              }
              counter++
            }
            Atom newA1 = atomList.get(index1)
            Atom newA2 = atomList.get(index2)
            Bond b = new Bond(newA1, newA2)
            b.setBondType(bond.getBondType())
          }

          // Construct the force field for the expanded set of molecules
          newAssembly.setForceField(assembly.getForceField());
          newAssembly.setPotential(assembly.potentialEnergy)
          newAssembly.setCrystal(assembly.getCrystal())
          newAssembly.setFile(assembly.getFile())
          Utilities.biochemistry(newAssembly, atomList)

          XYZFilter filter = new XYZFilter(saveFile, newAssembly, assembly.getForceField(),
                  assembly.getProperties())
          if (numModels > 1) {
            filter.writeFile(saveFile, true)
          } else {
            filter.writeFile(saveFile, false)
          }
          systemFilter.readNext(false, false)

          if (numModels > 1) {
            logger.info(format(" Saved structure %d to " + saveFile.name, i + 1))
          } else {
            logger.info(format(" Saved to " + saveFile.name))
          }
        }
      }
    } else if (match && filenames.size() > 1) {
        potentialFunctions.openAll(filenames[0]);
        SystemFilter systemFilter1 = potentialFunctions.getFilter()
        potentialFunctions.openAll(filenames[1]);
        SystemFilter systemFilter2 = potentialFunctions.getFilter()
        logger.info(" Matching atom order between files.");
        MolecularAssembly assembly1 = systemFilter1.getActiveMolecularSystem();
        MolecularAssembly assembly2 = systemFilter2.getActiveMolecularSystem();
        int[] molNum = assembly1.getMoleculeNumbers();
        Molecule[] molecules = assembly1.getMolecules()
      ArrayList<Atom> atomList = new ArrayList<>();
      int atomIndex = 1
      File saveLocation = new File(filenames[0])
      File saveFile = potentialFunctions.versionFile(saveLocation)
      //TODO Update assembly with new order and save file.
      MolecularAssembly newAssembly = new MolecularAssembly(assembly1.getName());
      // Construct the force field for the expanded set of molecules
      newAssembly.setForceField(assembly1.getForceField());
      newAssembly.setPotential(assembly1.potentialEnergy)
      newAssembly.setCrystal(assembly1.getCrystal())
      newAssembly.setFile(assembly1.getFile())

      Molecule[] molecules2 = assembly2.getMolecules()
      int numMol2 = molecules2.length;
      boolean[] molDone = new boolean[numMol2];
      int totalNumAtoms = 0;
        for (Molecule mol : molecules) {
          List<Atom> atoms1 = mol.getAtomList();
          int numAtoms = atoms1.size();
          logger.info(format(" Molecule1 with %3d atoms: %s", numAtoms, mol.toString()));

          // Determine atoms that are unique (i.e., different environments).
          List<List<Integer>> equivalents = new ArrayList<>();
          boolean[] equivalent = new boolean[numAtoms];
          int numEquiv = 0;
          for (int i = 0; i < numAtoms; i++){
            for (int j = 1; j < numAtoms; j++){
              if(i==j){
                continue;
              }
              AtomComparator ac = new AtomComparator()
              int value = ac.compare(atoms1[i], atoms1[j])
              if (value == 0){
                equivalent[i] = true;
                equivalent[j] = true;
                // Add to object.
                int iIndex = -1;
                int jIndex = -1;
                for (int k = 0; k < numEquiv; k++){
                  List<Integer> list = equivalents.get(k)
                  if (list.contains(i)){
                    iIndex = k;
                  }
                  if(list.contains(j)){
                    jIndex = k;
                  }
                  if(iIndex >= 0 && jIndex >= 0){
                    break;
                  }
                }
                if (iIndex == -1 && jIndex == -1){
                  List<Integer> ij = new ArrayList<>();
                  ij.add(i);
                  ij.add(j);
                  equivalents.add(ij);
                  numEquiv++;
                } else if (iIndex >= 0 && jIndex == -1){
                  equivalents.get(iIndex).add(j)
                } else if (iIndex == -1 && jIndex >= 0){
                  equivalents.get(jIndex).add(i);
                }
              }
            }
          }
          int numUnique = 0;
          for( boolean value: equivalent){
            if(!value){
              numUnique++;
            }
          }
          double[] xyz1 = new double[numAtoms * 3];
          double[] mass = new double[numAtoms];

          double[] compCoords1 = new double[numUnique * 3]
          double[] compMass = new double[numUnique];
          int index = 0;
          int indexUnique = 0;
          for (int i = 0; i < numAtoms; i++) {
            Atom a = atoms1[i];
            if(!equivalent[i]){
              compCoords1[indexUnique * 3] = a.getX();
              compCoords1[indexUnique * 3 + 1] = a.getY();
              compCoords1[indexUnique * 3 + 2] = a.getZ();
              compMass[indexUnique++] = a.getMass();
            }
            xyz1[index * 3] = a.getX();
            xyz1[index * 3 + 1] = a.getY();
            xyz1[index * 3 + 2] = a.getZ();
            mass[index++] = a.getMass();
          }
          for (int l = 0; l < numMol2 ; l++) {
            Molecule mol2 = molecules2[l];
            List<Atom> atoms2 = mol2.getAtomList();
            int numAtoms2 = atoms2.size();
            logger.info(format(" Molecule2 with %3d atoms: %s", numAtoms2, mol2.toString()));
            if (numAtoms != numAtoms2 || molDone[l]) {
              continue;
            }
            double[] xyz2 = new double[numAtoms2 * 3];
            double[] compCoords2 = new double[numUnique * 3]
            index = 0;
            indexUnique = 0;
            for (int i = 0; i < numAtoms; i++) {
              Atom a = atoms2[i];
              if(!equivalent[i]){
                compCoords2[indexUnique * 3] = a.getX();
                compCoords2[indexUnique * 3 + 1] = a.getY();
                compCoords2[indexUnique++ * 3 + 2] = a.getZ();
              }
              xyz2[index * 3] = a.getX();
              xyz2[index * 3 + 1] = a.getY();
              xyz2[index++ * 3 + 2] = a.getZ();
            }
            // If no atoms are unique, try to compare with all atoms.
            if( compCoords1.length == 0 || compCoords2.length == 0){
              logger.info(format(" Comparison Coordinates invalid (1: %3d 2: %3d) %3d of %3d unique atoms.", compCoords1.length, compCoords2.length, numUnique, numAtoms))
              compCoords1 = xyz1;
              compCoords2 = xyz2;
              compMass = mass;
            }
            double[] translate1 = calculateTranslation(compCoords1, compMass);
            double[] translate2 = calculateTranslation(compCoords2, compMass);
            applyTranslation(compCoords1, translate1);
            applyTranslation(compCoords2, translate2);
            // Translate all coordinates based on unique atoms (redundant if no unique).
            applyTranslation(xyz1, translate1);
            applyTranslation(xyz2, translate2);
            double[][] rotation = calculateRotation(compCoords1, compCoords2, compMass);
            // Rotate all coordinates based on unique atoms (redundant if no unique).
            applyRotation(compCoords2, rotation);
            applyRotation(xyz2, rotation);
            double compRMSD = superpose(compCoords1, compCoords2, compMass)
            logger.info(format(" Unique atom RMSD: %9.4f A", compRMSD));
            if(compRMSD > 1.0){
              logger.warning(format(" Value of %9.3f A may lead to insufficient overlap. Double check produced ordering.", compRMSD));
            }
            // Loop through atoms that are not unique and try to find optimal match.
            for (int i = 0; i < numAtoms; i++) {
              if(equivalent[i]) {
                int ind = i * 3;
                double[] coord1 = new double[]{xyz1[ind], xyz1[ind + 1], xyz1[ind + 2]};
                double[] coord2 = new double[]{xyz2[ind], xyz2[ind + 1], xyz2[ind + 2]}
                double value = rmsd(coord1, coord2, mass[i]);
                logger.info(format(" Atom %2d distance of %9.3f A", i, value));
                // TODO update following loop to reasonably update atoms ordering.
                for (List<Integer> list : equivalents) {
                  if (list.contains(i)) {
                    int minIndex = -1;
                    double minValue = MAX_VALUE;
                    for (Integer aIndex : list) {
                      if (i != aIndex) {
                        int ind2 = aIndex * 3;
                        double[] tempCoord = new double[]{xyz1[ind2], xyz1[ind2 + 1], xyz1[ind2 + 2]};
                        double value2 = rmsd(tempCoord, coord2, mass[i]);
                        if (value2 < minValue.doubleValue()) {
                          minIndex = aIndex
                          minValue = value2
                        }
                      }
                    }
                    if (minValue < value) {
                      // Update Atom
                      Atom temp = atoms1.get(i);
                      atoms1.set(i, atoms1.get(minIndex));
                      atoms1.set(minIndex, temp);
                      // Update Coords
                      double[] tempCoord = new double[]{xyz1[ind], xyz1[ind + 1], xyz1[ind + 2]}
                      int ind2 = minIndex * 3;
                      xyz1[ind] = xyz1[ind2];
                      xyz1[ind + 1] = xyz1[ind2 + 1];
                      xyz1[ind + 2] = xyz1[ind2 + 2];
                      xyz1[ind2] = tempCoord[0];
                      xyz1[ind2 + 1] = tempCoord[1];
                      xyz1[ind2 + 2] = tempCoord[2];
                    }
                    value = rmsd(new double[]{xyz1[ind], xyz1[ind + 1], xyz1[ind + 2]},
                            new double[]{xyz2[ind], xyz2[ind + 1], xyz2[ind + 2]}, mass[i]);
                    logger.info(format(" Atom %2d distance updated to %9.3f A", i, value));
                    // Index was found. No need to check the rest.
                    break;
                  }
                }
              }
            }
            compRMSD = superpose(xyz1, xyz2, mass)
            logger.info(format(" Final RMSD: %9.4f A", compRMSD));
            // Current molecule
            // Create a new set of Atoms for each molecule
            for (int j = 0; j < numAtoms; j++) {
              Atom a = atoms1[j]
              double[] xyz = [a.getX(), a.getY(), a.getZ()]
              Atom atom = new Atom(atomIndex++, a.getName(), a.getAtomType(), xyz)
              atomList.add(atom)
            }

            // Create a new set of Bonds for each molecule
            List<Bond> bondList = mol.getBondList();
            for (Bond bond : bondList) {
              Atom a1 = bond.getAtom(0)
              Atom a2 = bond.getAtom(1)

              int index1 = -1
              int index2 = -1
              int counter = 0
              for (Atom atom : atoms1) {
                int aIndex = atom.getIndex()
                if (aIndex == a1.getIndex()) {
                  index1 = counter
                } else if (aIndex == a2.getIndex()) {
                  index2 = counter
                }
                counter++
              }
              if(index1 == -1 || index2 == -1){
                logger.severe(format(" Index1 (%3d) and Index2 (%3d) suggests error for bond: %s", index1, index2, bond))
              }
              Atom newA1 = atomList.get(index1 + totalNumAtoms)
              Atom newA2 = atomList.get(index2 + totalNumAtoms)
              logger.info(format(" index1 %3d Index2 %3d", index1, index2));
              if(!newA1.isBonded(newA2)){
                Bond b = new Bond(newA1, newA2)
                b.setBondType(bond.getBondType())
              }
            }
            molDone[l] = true;
            totalNumAtoms += numAtoms;
          }
          if(!Arrays.asList(molDone).contains(false)){
            break;
          }
        }
      Utilities.biochemistry(newAssembly, atomList)
      XYZFilter filter = new XYZFilter(saveFile, newAssembly, newAssembly.getForceField(),
              newAssembly.getProperties())
      filter.writeFile(saveFile, true)
      logger.info(format(" Saved to " + saveFile.name))

      // Merge molecules into assembly and print out.
      } else if (atomType && filenames.size() > 1){
        potentialFunctions.openAll(filenames[0]);
        SystemFilter systemFilter1 = potentialFunctions.getFilter()
        potentialFunctions.openAll(filenames[1]);
        SystemFilter systemFilter2 = potentialFunctions.getFilter()
        logger.info(" Changing atom types in file 1 to match file 2.");
        // Atom types from desired file.
        MolecularAssembly assembly1 = systemFilter1.getActiveMolecularSystem();
        do {
          MolecularAssembly assembly2 = systemFilter2.getActiveMolecularSystem();
          Atom[] atoms1 = assembly1.getAtomArray();
          Atom[] atoms2 = assembly2.getAtomArray();
          ArrayList<ArrayList<Atom>> xyzatoms = new ArrayList<>();
          MolecularAssembly outputAssembly = new MolecularAssembly(assembly1.getName());

          int numHydrogen = 0;
          for (Atom atom : atoms2) {
            if (atom.isHydrogen()) {
              numHydrogen++;
            }
          }

          List<MSNode> entitiesInput = assembly2.getAllBondedEntities();
          int[] molIndices = assembly2.getMoleculeNumbers();
          int numEntitiesInput = entitiesInput.size();
          int[] shiftIndices = new int[numEntitiesInput];
          int count = 0;
          for(int i = 0; i < numEntitiesInput; i++){
            for(int ind : molIndices){
              if(ind == i){
                count++;
              }
            }
            shiftIndices[i] = count;
          }
          if (logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Number of entities in input: %d", numEntitiesInput));
            for (MSNode entity : entitiesInput) {
              logger.fine(format(" Entity: " + entity.getName()));
              int size = entity.getAtomList().size();
              logger.fine(format("   Entity Size: %3d", size));
              if (size > 0) {
                logger.fine(format("   Entity First Atom: " + entity.getAtomList().get(0).toString()));
              } else {
                logger.warning(" Entity did not contain atoms.");
              }
            }
          }

          for (MSNode mol : entitiesInput) {
            xyzatoms.add((ArrayList<Atom>) mol.getAtomList());
          }
          int atomIndex = 1;
          // For each entity (protein, molecule, ion, etc) in the XYZ
          logger.info(" Num Entities Input: " + numEntitiesInput);
          for (int i = 0; i < numEntitiesInput; i++) {
            MSNode mol = entitiesInput.get(i);
            int numInputMolAtoms = xyzatoms.get(i).size();
            int numMolHydrogen = 0;
            for (Atom atom : xyzatoms.get(i)) {
              if (atom.isHydrogen()) {
                numMolHydrogen++;
              }
            }
            if (logger.isLoggable(Level.FINE)) {
              logger.fine(format(" Current entity number of atoms: %d (%d + %dH)", numInputMolAtoms,
                      numInputMolAtoms - numMolHydrogen, numMolHydrogen));
            }

            // Set up input (XYZ:assembly 2) file contents as CDK variable
            AtomContainer xyzCDKAtoms = new AtomContainer();
            for (Atom atom : xyzatoms.get(i)) {
              String atomName = getSymbol(atom.getAtomType().atomicNumber);
              xyzCDKAtoms.addAtom(
                      new org.openscience.cdk.Atom(atomName, new Point3d(atom.getXYZ(null))));
            }

            int zPrime;
            int nAtoms = atoms1.length;
            int nInputAtoms = atoms2.length;
            if (nAtoms % nInputAtoms == 0) {
              zPrime = (int) (nAtoms / nInputAtoms);
            } else if (nAtoms % (nInputAtoms - numHydrogen) == 0) {
              zPrime = (int) (nAtoms / (nInputAtoms - numHydrogen));
            } else {
              zPrime = 1;
            }
            logger.info(" Original ZPrime: " + zPrime);
            // Add known input bonds; a limitation is that bonds all are given a Bond order of 1.
            List<Bond> bonds = mol.getBondList();
            IBond.Order order = IBond.Order.SINGLE;
            int xyzBonds = bonds.size();
            if (xyzBonds == 0) {
              logger.warning(" No bonds detected in input structure. Please check input.\n " +
                      "If correct, separate non-bonded entities into multiple CIFs.");
            }
            int shiftIndex;
            for (Bond xyzBond : bonds) {
              int atom0Index = xyzBond.getAtom(0).getXyzIndex();
              int atom1Index = xyzBond.getAtom(1).getXyzIndex();
              shiftIndex = (molIndices[atom0Index]==0)?0:shiftIndices[molIndices[atom0Index]-1];
              if (logger.isLoggable(Level.FINER)) {
                logger.finer(format(" Bonded atom 1: %d, Bonded atom 2: %d", atom0Index, atom1Index));
                logger.finer(format(" Atom0: %3d Atom1: %3d Shift: %3d final0,1: %3d, %3d MolIndex: %3d ShiftIndex: %3d",
                        atom0Index, atom1Index, shiftIndex, atom0Index - shiftIndex - 1, atom1Index - shiftIndex - 1,
                        molIndices[atom0Index],shiftIndices[molIndices[atom0Index]]))
              }
              xyzCDKAtoms.addBond(atom0Index - shiftIndex - 1, atom1Index - shiftIndex - 1, order);
            }

            // Assign CDK atom types for the input molecule.
            AtomTypeFactory factory = AtomTypeFactory.getInstance(
                    "org/openscience/cdk/config/data/jmol_atomtypes.txt", xyzCDKAtoms.getBuilder());

            // Assign atom types to CDK object.
            for (IAtom atom : xyzCDKAtoms.atoms()) {
              CIFFilter.setAtomTypes(factory, atom);
              try {
                factory.configure(atom);
              } catch (Exception ex) {
                logger.info(" Failed to configure atoms from CIF.\n" + ex + "\n" + Arrays.toString(
                        ex.getStackTrace()));
              }
            }

            ArrayList<ArrayList<Integer>> zindices = new ArrayList<>();
            int counter = 0;
            // Bond atoms from CIF.
            int cifBonds = CIFFilter.bondAtoms(atoms1, bondTolerance);
            if (logger.isLoggable(Level.FINE)) {
              logger.fine(
                      format(" Created %d bonds between input atoms (%d in input).", cifBonds, xyzBonds));
            }
            List<Atom> atomPool = new ArrayList<>(Arrays.asList(atoms1));

            try {
              while (!atomPool.isEmpty()) {
                ArrayList<Atom> molecule = new ArrayList<>();
                CIFFilter.collectAtoms(atomPool.get(0), molecule);
                if (logger.isLoggable(Level.FINER)) {
                  logger.finer(format(" Molecule (%d) Size: %d", counter, molecule.size()));
                }
                ArrayList<Integer> indices = new ArrayList<>();
                while (molecule.size() > 0) {
                  Atom atom = molecule.remove(0);
                  indices.add(atom.getIndex());
                  atomPool.remove(atom);
                }

                if (logger.isLoggable(Level.FINER)) {
                  logger.finer(format(
                          " Molecule %d: %d atoms are ready and %d remain (%d atoms in input, %d atoms in CIF). ",
                          counter + 1, indices.size(), atomPool.size(), nInputAtoms, nAtoms));
                }
                zindices.add(indices);
                counter++;
              }
            } catch (Exception e) {
              logger.severe(" Failed to separate copies within the asymmetric unit." + e + "\n" + Utilities.stackTraceToString(e));
            }
            zPrime = zindices.size();
            logger.info(" Zprime: " + zPrime);
            // Set up CIF contents as CDK variable
            AtomContainer[] cifCDKAtomsArr = new AtomContainer[zPrime];
            AtomContainer cifCDKAtoms = new AtomContainer();
            for (int j = 0; j < zPrime; j++) {
              if (zPrime > 1) {
                logger.info(format("\n Attempting entity %d of %d", j + 1, zPrime));
              }
              ArrayList<Integer> currentList = zindices.get(j);
              int cifMolAtoms = currentList.size();
              if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" CIF atoms in current: %d", cifMolAtoms));
              }
              // Detect if CIF contains multiple copies (Z'>1)
              if (cifMolAtoms % numInputMolAtoms == 0 || cifMolAtoms % (numInputMolAtoms - numMolHydrogen) == 0) {
                cifCDKAtomsArr[j] = new AtomContainer();
                for (Integer integer : currentList) {
                  String atomName = getSymbol(atoms1[integer - 1].getAtomType().atomicNumber);
                  cifCDKAtomsArr[j].addAtom(new org.openscience.cdk.Atom(atomName, new Point3d(atoms1[integer - 1].getXYZ(null))));
                }
                AtomContainer nullCDKAtoms = new AtomContainer();

                for (IAtom atom : cifCDKAtomsArr[j].atoms()) {
                  if (atom.toString() == null) {
                    nullCDKAtoms.addAtom(atom);
                  }
                }
                cifCDKAtomsArr[j].remove(nullCDKAtoms);

                // Compute bonds for CIF molecule.
                factory = AtomTypeFactory.getInstance(
                        "org/openscience/cdk/config/data/jmol_atomtypes.txt", cifCDKAtomsArr[j].getBuilder());
                for (IAtom atom : cifCDKAtomsArr[j].atoms()) {
                  CIFFilter.setAtomTypes(factory, atom);
                  try {
                    factory.configure(atom);
                  } catch (Exception ex) {
                    logger.info(" Failed to configure CIF atoms.\n" + ex + "\n" + Utilities.stackTraceToString(ex));
                  }
                }

                RebondTool rebonder = new RebondTool(CIFFilter.MAX_COVALENT_RADIUS, CIFFilter.MIN_BOND_DISTANCE,
                        bondTolerance);
                try {
                  rebonder.rebond(cifCDKAtomsArr[j]);
                } catch (Exception ex) {
                  logger.info("Failed to rebond CIF atoms.\n" + ex + "\n" + Utilities.stackTraceToString(ex));
                }

                int cifMolBonds = cifCDKAtomsArr[j].getBondCount();
                if (logger.isLoggable(Level.FINE)) {
                  logger.fine(format(" Number of CIF bonds: %d (%d in input)", cifMolBonds, xyzBonds));
                }
                // Number of bonds matches.
                // If cifMolBonds == 0 then ion or atom with implicit hydrogens (e.g. water, methane, etc.)
                if (cifMolBonds != 0 && cifMolBonds % xyzBonds == 0) {
                  Pattern pattern = VentoFoggia.findIdentical(
                          xyzCDKAtoms, AtomMatcher.forElement(), BondMatcher.forAny());
                  int[] p = pattern.match(cifCDKAtomsArr[j]);
                  int pLength = p.length;
                  if (p != null && pLength == numInputMolAtoms) {
                    // Use matched atoms to update the positions of the input file atoms.
                    for (int k = 0; k < pLength; k++) {
                      if (logger.isLoggable(Level.FINEST)) {
                        logger.finest(
                                format(" %d input %s -> CIF %s", k, xyzCDKAtoms.getAtom(k).getSymbol(),
                                        cifCDKAtomsArr[j].getAtom(p[k]).getSymbol()));
                      }
                      Point3d point3d = cifCDKAtomsArr[j].getAtom(p[k]).getPoint3d();
                      xyzatoms.get(i).get(k).setXYZ(new double[]{point3d.x, point3d.y, point3d.z});
                    }
                  } else {
                    if (logger.isLoggable(Level.FINE)) {
                      logger.fine(
                              format(" Atoms from CIF (%d) and input (%d) structures don't match.", p.length,
                                      nAtoms));
                    }
                    continue;
                  }
                } else if ((xyzBonds - numMolHydrogen) == 0 || cifMolBonds % ((xyzBonds - numMolHydrogen)) == 0) {
                  // Hydrogen most likely missing from file. If zero, then potentially water/methane (implicit hydrogen atoms).
                  if (logger.isLoggable(Level.FINE)) {
                    logger.info(" CIF may contain implicit hydrogen -- attempting to patch.");
                  }
                  // Match heavy atoms between CIF and input
                  Pattern pattern = VentoFoggia.findSubstructure(
                          cifCDKAtomsArr[j], AtomMatcher.forElement(), BondMatcher.forAny());
                  int[] p = pattern.match(xyzCDKAtoms);
                  int pLength = p.length;
                  if (p != null && pLength == numInputMolAtoms - numMolHydrogen) {
                    // Used matched atoms to update the positions of the input file atoms.
                    for (int k = 0; k < pLength; k++) {
                      Point3d point3d = cifCDKAtomsArr[j].getAtom(k).getPoint3d();
                      xyzatoms.get(i).get(p[k]).setXYZ(new double[]{point3d.x, point3d.y, point3d.z});
                    }
                    // Add in hydrogen atoms
                    Atom lastKnownAtom1 = null;
                    int knownCount = 0;
                    for (Atom hydrogen : xyzatoms.get(i)) {
                      if (hydrogen.isHydrogen()) {
                        Bond bond0 = hydrogen.getBonds().get(0);
                        Atom atom1 = bond0.get1_2(hydrogen);
                        double angle0_2 = 0;
                        List<Angle> anglesList = hydrogen.getAngles(); // Same as bond
                        Atom atom2 = null;
                        switch (anglesList.size()) {
                          case 0 ->
                            // H-Cl No angles
                            // Place hydrogen slightly off center of bonded atom (~1Å away).
                            hydrogen.moveTo(new double[]{atom1.getX() - 0.6, atom1.getY() - 0.6,
                                    atom1.getZ() - 0.6});
                          case 1 -> {
                            // H-O=C Need torsion
                            // H-O-H
                            for (Angle angle : anglesList) {
                              atom2 = angle.get1_3(hydrogen);
                              if (atom2 != null) {
                                angle0_2 = angle.angleType.angle[0];
                              }
                              if (angle0_2 != 0) {
                                //if angle0_2 is found then no need to look for another atom.
                                break;
                              }
                            }
                            assert atom2 != null;
                            List<Bond> bonds2 = atom2.getBonds();
                            Atom atom3 =
                                    (bonds2.size() > 1 && atom1 == bonds2.get(0).get1_2(atom2)) ? bonds2.get(
                                            1).get1_2(atom2) : bonds2.get(0).get1_2(atom2);
                            double diAng = dihedralAngle(hydrogen.getXYZ(null), atom1.getXYZ(null),
                                    atom2.getXYZ(null), atom3.getXYZ(null));
                            if (atom1 != atom3) {
                              intxyz(hydrogen, atom1, bond0.bondType.distance, atom2, angle0_2, atom3,
                                      Math.toDegrees(diAng), 0);
                            } else {
                              // Likely water as Atom 3 is not unique. Since no hydrogen atoms
                              //  are present, there isn't a third atom...
                              double[] coord = new double[]{atom2.getX(), atom2.getY(), atom3.getZ()};
                              double mag = sqrt(
                                      coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
                              coord[0] /= mag;
                              coord[1] /= mag;
                              coord[2] /= mag;
                              hydrogen.moveTo(atom1.getX() - coord[0], atom1.getY() - coord[1],
                                      atom1.getZ() - coord[2]);
                            }
                          }
                          default -> {
                            // H-C(-C)(-C)-C
                            Atom atom2B = null;
                            double angle0_2B = 0;
                            Atom proposedAtom;
                            double proposedAngle = 0;
                            int chiral = 1;
                            for (Angle angle : anglesList) {
                              proposedAtom = angle.get1_3(hydrogen);
                              if (proposedAtom != null && !proposedAtom.isHydrogen()) {
                                proposedAngle = angle.angleType.angle[0];
                              }
                              if (proposedAngle != 0) {
                                // If angle1_3 is found then no need to look for another atom.
                                if (angle0_2 != 0) {
                                  atom2B = proposedAtom;
                                  angle0_2B = proposedAngle;
                                } else {
                                  atom2 = proposedAtom;
                                  angle0_2 = proposedAngle;
                                }
                                proposedAngle = 0.0;
                              }
                              if (angle0_2 != 0 && angle0_2B != 0) {
                                break;
                              }
                            }
                            if (lastKnownAtom1 == null || lastKnownAtom1 != atom1) {
                              lastKnownAtom1 = atom1;
                              knownCount = 0;
                            } else {
                              knownCount++;
                            }
                            if (angle0_2B == 0.0) {
                              // Hydrogen position depends on other hydrogen, use generic location
                              chiral = 0;
                              assert atom2 != null;
                              List<Bond> bonds2 = atom2.getBonds();
                              atom2B =
                                      (atom1 == bonds2.get(0).get1_2(atom2)) ? bonds2.get(1).get1_2(atom2)
                                              : bonds2.get(0).get1_2(atom2);
                              // Evenly space out hydrogen
                              angle0_2B = 180.0 - 120.0 * knownCount;
                            } else if (anglesList.size() == 2) {
                              // Trigonal hydrogen use equipartition between selected atoms
                              chiral = 3;
                            } else if (knownCount == 1) {
                              // Already created one hydrogen (chiral = 1), use other.
                              chiral = -1;
                            }
                            //TODO discern whether to use chiral = 1 or -1 when angles are to three heavy atoms.
                            intxyz(hydrogen, atom1, bond0.bondType.distance, atom2, angle0_2, atom2B,
                                    angle0_2B, chiral);
                          }
                        }
                      }
                    }
                  } else {
                    if (logger.isLoggable(Level.FINE)) {
                      logger.fine(" Could not match heavy atoms between CIF and input.");
                    }
                    if (p != null && logger.isLoggable(Level.FINE)) {
                      logger.fine(
                              format(" Matched %d atoms out of %d in CIF (%d in input)", pLength, nAtoms,
                                      nInputAtoms - numMolHydrogen));
                    }
                    continue;
                  }
                } else {
                  logger.info(format(" CIF (%d) and input ([%d+%dH=]%d) have a different number of bonds.",
                          cifMolBonds, xyzBonds - numMolHydrogen, numMolHydrogen, xyzBonds));
                  continue;
                }
                cifCDKAtoms.add(cifCDKAtomsArr[j]);
                MSGroup molecule;
                if (mol instanceof Polymer) {
                  molecule = new Polymer(((Polymer) mol).getChainID(), mol.getName(), true);
                } else {
                  molecule = new Molecule(mol.getName());
                }
                List<Atom> atomList = new ArrayList<>();
                for (Atom atom : xyzatoms.get(i)) {
                  if (logger.isLoggable(Level.FINER)) {
                    logger.finer(format(" Atom Residue: %s %2d", atom.getResidueName(), atom.getResidueNumber()));
                  }
                  Atom molAtom;
                  if (molecule instanceof Polymer) {
                    molAtom = new Atom(atomIndex++, atom.getName(), atom.getAltLoc(), atom.getXYZ(null),
                            atom.getResidueName(), atom.getResidueNumber(), atom.getChainID(), atom.getOccupancy(),
                            atom.getTempFactor(), atom.getSegID());
                    molAtom.setAtomType(atom.getAtomType());
                  } else {
                    molAtom = new Atom(atomIndex++, atom.getName(), atom.getAtomType(), atom.getXYZ(null));
                    if (atom.getResidueName() != null) {
                      molAtom.setResName(atom.getResidueName());
                    }
                  }
                  atomList.add(molAtom);
                }
                List<Bond> bondList = mol.getBondList();
                for (Bond bond : bondList) {
                  Atom a1 = bond.getAtom(0);
                  Atom a2 = bond.getAtom(1);
                  if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" Bonded atom 1: %d, Bonded atom 2: %d", a1.getXyzIndex(),
                            a2.getXyzIndex()));
                  }
                  Atom newA1 = atomList.get(a1.getIndex() - shiftIndex - 1);
                  Atom newA2 = atomList.get(a2.getIndex() - shiftIndex - 1);
                  Bond bond2 = new Bond(newA1, newA2);
                  bond2.setBondType(bond.getBondType());
                }
                Residue res = null;
                for (Atom atom : atomList) {
                  if (molecule instanceof Polymer) {
                    if (res == null) {
                      res = new Residue(atom.getResidueName(), atom.getResidueNumber(), ((Polymer) mol).getResidue(atom.getResidueNumber()).getResidueType());
                    } else if (res.getResidueNumber() != atom.getResidueNumber()) {
                      if (logger.isLoggable(Level.FINER)) {
                        logger.finer(" Added Residue " + res.getResidueNumber() + " " + res.getName());
                      }
                      molecule.addMSNode(res);
                      res = new Residue(atom.getResidueName(), atom.getResidueNumber(), ((Polymer) mol).getResidue(atom.getResidueNumber()).getResidueType());
                    }
                    res.addMSNode(atom);
                  } else {
                    molecule.addMSNode(atom);
                  }
                }
                if (molecule instanceof Polymer && res != null) {
                  if (logger.isLoggable(Level.FINER)) {
                    logger.finer(" Added Final Residue " + res.getResidueNumber() + " " + res.getName());
                  }
                  molecule.addMSNode(res);
                }
                outputAssembly.addMSNode(molecule);
              } else {
                if (logger.isLoggable(Level.INFO)) {
                  logger.info(
                          format(" Number of atoms in CIF (%d) molecule do not match input (%d + %dH = %d).",
                                  cifMolAtoms, nInputAtoms - numMolHydrogen, numMolHydrogen, nInputAtoms));
                }
              }
            }
          }
          // If no atoms, then conversion has failed... use active assembly
          if (logger.isLoggable(Level.FINE)) {
            logger.fine(format("\n Output Assembly Atoms: %d", outputAssembly.getAtomList().size()));
          }
          outputAssembly.setPotential(assembly2.getPotentialEnergy());
          if (assembly1.getCrystal() != null) {
            outputAssembly.setCrystal(assembly1.getCrystal());
          }
          outputAssembly.setForceField(assembly2.getForceField());
          outputAssembly.setFile(assembly1.getFile());
          outputAssembly.setName(assembly1.getName());
//          SystemFilter.setMolecularSystem(outputAssembly);

          if (outputAssembly.getAtomList().size() < 1) {
            logger.info(" Atom types could not be matched. File could not be written.");
          } else if (!writeOutputFile(outputAssembly, systemFilter1)) {
            logger.info(" Output assembly file could not be written.");
          }
          // Finsished with this block. Clean up/reset variables for next block.
          outputAssembly.destroy();
        } while (systemFilter1.readNext());
      }
    return this
  }

  /**
   * Write molecular assembly to an output file.
   * @param ma Assembly to save
   * @param file File (for path and name base)
   * @return
   */
  private static boolean writeOutputFile(MolecularAssembly ma, SystemFilter systemFilter){
    File file = systemFilter.getFiles().get(0);
    String dir = getFullPath(file.getAbsolutePath()) + File.separator;
    String fileName = removeExtension(getName(file.getAbsolutePath()));
    String spacegroup;
    if(ma.getCrystal() != null) {
      spacegroup = ma.getCrystal().getUnitCell().spaceGroup.shortName;
    }else{
      spacegroup = null;
    }
    List<MSNode> entities = ma.getAllBondedEntities();
    File saveFile;
    if(systemFilter instanceof PDBFilter){
      // Change name for different space groups. Input cannot handle multiple space groups in same file.
      if (entities.size() > 1) {
        fileName += "_z" + entities.size();
      }
      saveFile = new File(dir + fileName + ".pdb");
    }  else {
      if (systemFilter.countNumModels() > 1) {

        saveFile = new File(dir + fileName + "_reorder.arc");
      } else {
        // If only structure then create XYZ.
        saveFile = XYZFilter.version(new File(dir + fileName + ".xyz"));
      }
    }
    ForceField ff = ma.getForceField();
    CompositeConfiguration properties = ma.getProperties();
    if(systemFilter instanceof PDBFilter){
      PDBFilter pdbFilter = new PDBFilter(saveFile, ma, ff, properties);
      pdbFilter.writeFile(saveFile, true);
    }else {
      XYZFilter xyzFilter = new XYZFilter(saveFile, ma, ff, properties);
      xyzFilter.writeFile(saveFile, true);
    }
    logger.info("\n Saved output file:        " + saveFile.getAbsolutePath());
    File propertyFile = new File(removeExtension(saveFile.getName()) + ".properties");
    if (!propertyFile.exists()) {
      try {
        FileWriter fw = new FileWriter(propertyFile, false);
        BufferedWriter bw = new BufferedWriter(fw);
        ForceField forceField = ma.getForceField();
        String parameters = forceField.getString("parameters", "none");
        if (parameters != null && !parameters.equalsIgnoreCase("none")) {
          bw.write(format("parameters %s\n", parameters));
        }
        String forceFieldProperty = forceField.getString("forcefield", "none");
        if (forceFieldProperty != null && !forceFieldProperty.equalsIgnoreCase("none")) {
          bw.write(format("forcefield %s\n", forceFieldProperty));
        }
        String patch = forceField.getString("patch", "none");
        if (patch != null && !patch.equalsIgnoreCase("none")) {
          bw.write(format("patch %s\n", patch));
        }
        bw.write(format("spacegroup %s\n", spacegroup));
        if (entities.size() > 1) {
          bw.write("intermolecular-softcore true\n");
        }
        logger.info("\n Saved properties file: " + propertyFile.getAbsolutePath() + "\n");
        bw.close();
      } catch (Exception ex) {
        logger.info("Failed to write files.\n" + Utilities.stackTraceToString(ex));
      }
    } else {
      logger.info("\n Property file already exists:  " + propertyFile.getAbsolutePath() + "\n");
    }
    return true;
  }

  /**
   * Compare atoms based on connectivity and atom types to determine order.
   */
  class AtomComparator implements Comparator {
    /**
     * Compare two atoms based on their connections check environments.
     * @param obj1
     * @param obj2
     * @return 0 is uncertain, -1 if reverse, 1 if same.
     */
    int compare(Object obj1, Object obj2) {
      Atom a1 = (Atom) obj1
      Atom a2 = (Atom) obj2
      int comp = Integer.compare(a1.getMoleculeNumber(), a2.getMoleculeNumber())
      if (logger.isLoggable(Level.FINER) && comp != 0) {
        logger.finer(format(" Different Molecule Number (%d vs %d)", a1.getMoleculeNumber(),
            a2.getMoleculeNumber()))
      }
      if (comp == 0) {
        comp = Integer.compare(a1.getResidueNumber(), a2.getResidueNumber())
        if (logger.isLoggable(Level.FINER) && comp != 0) {
          logger.finer(format(" Different Residue Numbers %d vs %d", a1.getResidueNumber(),
              a2.getResidueNumber()))
        }
        if (comp == 0) {
          // Want heavier atoms first so * -1.
          comp = -Double.compare(a1.getMass(), a2.getMass());
          if (logger.isLoggable(Level.FINER) && comp != 0) {
            logger.finer(format(" Different Masses %6.3f %6.3f", a1.getMass(), a2.getMass()))
          }
          if (comp == 0) {
            comp = Integer.compare(a1.getAtomType().type, a2.getAtomType().type)
            if (logger.isLoggable(Level.FINER) && comp != 0) {
              logger.finer(
                  format(" Different Atom Types (%d vs %s)", a1.atomType.type, a2.atomType.type))
            }
            if (comp == 0) {
              comp = Integer.compare(a1.getBonds().size(), a2.getBonds().size())
              if (logger.isLoggable(Level.FINER) && comp != 0) {
                logger.finer(format(" Different Bond Sizes (%d vs %d)", a1.getBonds().size(),
                    a2.getBonds().size()))
              }
              if (comp == 0) {
                comp = Integer.compare(a1.get12List().size(), a2.get12List().size())
                if (logger.isLoggable(Level.FINER) && comp != 0) {
                  logger.finer(
                      format(" Different 1-2 Size (%d vs %d).", a1.get12List(), a2.get12List()))
                }
                if (comp == 0) {
                  for (int i = 0; i < a1.get12List().size(); i++) {
                    Atom a12 = a1.get12List()[i]
                    Atom a22 = a2.get12List()[i]
                    comp = shallowCompare(a12, a22)
                    if (comp != 0) {
                      break
                    }
                  }
                  if (logger.isLoggable(Level.FINER) && comp != 0) {
                    logger.finer(" Different 1-2 Shallow.")
                  }
                  if (comp == 0) {
                    comp = Integer.compare(a1.get13List().size(), a2.get13List().size())
                    if (logger.isLoggable(Level.FINER) && comp != 0) {
                      logger.finer(" Different 1-3 Size.")
                    }
                    if (comp == 0) {
                      for (int i = 0; i < a1.get13List().size(); i++) {
                        Atom a12 = a1.get13List()[i]
                        Atom a22 = a2.get13List()[i]
                        comp = shallowCompare(a12, a22)
                        if (comp != 0) {
                          break
                        }
                      }
                      if (logger.isLoggable(Level.FINER) && comp != 0) {
                        logger.finer(" Different 1-3 Shallow.")
                      }
                      if (comp == 0) {
                        comp = Integer.compare(a1.get14List().size(), a2.get14List().size())
                        if (logger.isLoggable(Level.FINER) && comp != 0) {
                          logger.finer(" Different 1-4 Size.")
                        }
                        if (comp == 0) {
                          for (int i = 0; i < a1.get14List().size(); i++) {
                            Atom a12 = a1.get14List()[i]
                            Atom a22 = a2.get14List()[i]
                            comp = shallowCompare(a12, a22)
                            if (comp != 0) {
                              break
                            }
                          }
                          if (logger.isLoggable(Level.FINER) && comp != 0) {
                            logger.finer(" Different 1-4 Shallow.")
                          }
                          if (comp == 0) {
                            comp = Integer.compare(a1.get15List().size(), a2.get15List().size())
                            if (logger.isLoggable(Level.FINER) && comp != 0) {
                              logger.finer(" Different 1-5 Size.")
                            }
                            if (comp == 0) {
                              for (int i = 0; i < a1.get15List().size(); i++) {
                                Atom a12 = a1.get15List()[i]
                                Atom a22 = a2.get15List()[i]
                                comp = shallowCompare(a12, a22)
                                if (comp != 0) {
                                  break
                                }
                              }
                              if (logger.isLoggable(Level.FINER) && comp != 0) {
                                logger.finer(" Different 1-5 Shallow.")
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      return comp
    }

    /**
     * Compare two atoms based on their bonded atoms' connections.
     * @param a1 Atom 1
     * @param a2 Atom 2
     * @return 0 is uncertain, -1 if reverse, 1 if same.
     */
    private int shallowCompare(Atom a1, Atom a2) {
      int comp = Integer.compare(a1.getBonds().size(), a2.getBonds().size());
      if (logger.isLoggable(Level.FINER) && comp != 0) {
        logger.finer(" Different Shallow Bond Size.")
      }
      if (comp == 0) {
        comp = Integer.compare(a1.atomType.type, a2.getAtomType().type)
        if (logger.isLoggable(Level.FINER) && comp != 0) {
          logger.finer(" Different Shallow Atom Type.")
        }
        if (comp == 0) {
          comp = Integer.compare(a1.get12List().size(), a2.get12List().size());
          if (logger.isLoggable(Level.FINER) && comp != 0) {
            logger.finer(" Different Shallow 1-2 Size.")
          }
          if (comp == 0) {
            comp = Integer.compare(a1.get13List().size(), a2.get13List().size());
            if (logger.isLoggable(Level.FINER) && comp != 0) {
              logger.finer(" Different Shallow 1-3 Size.")
            }
            if (comp == 0) {
              comp = Integer.compare(a1.get14List().size(), a2.get14List().size());
              if (logger.isLoggable(Level.FINER) && comp != 0) {
                logger.finer(" Different Shallow 1-4 Size.")
              }
              if (comp == 0) {
                comp = Integer.compare(a1.get15List().size(), a2.get15List().size());
                if (logger.isLoggable(Level.FINER) && comp != 0) {
                  logger.finer(" Different Shallow 1-5 Size.")
                }
              }
            }
          }
        }
      }
      return comp;
    }
  }
}