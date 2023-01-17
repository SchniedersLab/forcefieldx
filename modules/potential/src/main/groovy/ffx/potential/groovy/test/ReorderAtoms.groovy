//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy.test

import ffx.potential.MolecularAssembly
import ffx.potential.Utilities
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XYZFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level

import static java.lang.String.format

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
   * --sw or --swap Switch positions of all equivalent atoms (untested for >2 atoms).
   */
  @Option(names = ['--sw', '--swap'], paramLabel = "false", defaultValue = "false",
      description = 'Switch positions of equivalent atoms (untested for >2 atoms).')
  private static boolean swap

  /**
   * -o or --overwrite Attempt to replace original file with reordered coordinates.
   */
  @Option(names = ['-o', '--overwrite'], paramLabel = "false", defaultValue = "false",
      description = 'Replace original file with output (CAUTION: deletes input file before write is complete).')
  private static boolean overwrite

  /**
   * The final argument(s) should be two or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
      description = 'Atomic coordinate files to reorder in XYZ/ARC format.')
  List<String> filenames = null

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
    for (String filename : filenames) {
      potentialFunctions.openAll(filename)
      SystemFilter systemFilter = potentialFunctions.getFilter()
      int numModels = systemFilter.countNumModels()
      // Save new xyz file with the converted crystal.
      File saveLocation = new File(filename)
      File saveFile
      if (overwrite) {
        saveFile = saveLocation
      } else {
        saveFile = potentialFunctions.versionFile(saveLocation)
      }
      for (int i = 0; i < numModels; i++) {
        MolecularAssembly assembly = systemFilter.getActiveMolecularSystem()
        Atom[] atoms = assembly.getAtomArray()
        List<Atom> atomList0 = assembly.getAtomList()
        int nAtoms = atoms.length
        MolecularAssembly newAssembly = new MolecularAssembly(assembly.getName())
        List<Bond> bondList = assembly.getBondList()

        // TODO handle sorting/exclusion better...
        int alSize = atomList0.size()
        for (int j = 0; j < alSize; j++) {
          for (int k = j; k < alSize; k++) {
            if (atomList0[j].getIndex() != atomList0[k].getIndex()) {
              AtomComparator ac = new AtomComparator()
              if (logger.isLoggable(Level.FINER)) {
                logger.finer(
                    " Compare atom: " + atomList0[j].toString() + " and " + atomList0[k].toString())
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
              } else if (value < 0) {
                if (logger.isLoggable(Level.FINER)) {
                  logger.finer(
                      " Swapped " + atomList0[j].toString() + " and " + atomList0[k].toString() +
                          ".")
                }
                Atom temp = atomList0[j]
                atomList0[j] = atomList0[k]
                atomList0[k] = temp
              }
            }
          }
        }
        int[] originalOrder = new int[nAtoms]
        for (int j = 0; j < nAtoms; j++) {
          originalOrder[j] = atoms[j].getIndex()
          if (logger.isLoggable(Level.FINE)) {
            logger.fine(format(" Original index (%3d): %3d New index: %3d", j, originalOrder[j],
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

        if (overwrite && saveFile.exists()) {
          if (saveFile.delete()) {
            logger.info(" Original file deleted.")
          } else {
            logger.info(" Original file could not be deleted.")
          }
        }
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
    return this
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
          comp = Double.compare(a1.getMass(), a2.getMass())
          if (logger.isLoggable(Level.FINER) && comp != 0) {
            logger.finer(format(" Different Masses %6.3f %6.3f", a1.getMass(), a2.getMass()))
          }
          if (comp == 0) {
            comp = -Integer.compare(a1.atomType.type, a2.getAtomType().type)
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