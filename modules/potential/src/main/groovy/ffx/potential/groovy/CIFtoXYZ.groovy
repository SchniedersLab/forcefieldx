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
package ffx.potential.groovy

import ffx.crystal.Crystal
import ffx.crystal.SpaceGroup
import ffx.crystal.SpaceGroupDefinitions
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Angle
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.ForceField

import org.apache.commons.io.FilenameUtils
import org.openscience.cdk.AtomContainer
import org.openscience.cdk.config.AtomTypeFactory
import org.openscience.cdk.graph.rebond.RebondTool
import org.openscience.cdk.interfaces.IAtom
import org.openscience.cdk.interfaces.IAtomType
import org.openscience.cdk.interfaces.IBond.Order
import org.openscience.cdk.isomorphism.AtomMatcher
import org.openscience.cdk.isomorphism.BondMatcher
import org.openscience.cdk.isomorphism.Pattern
import org.openscience.cdk.isomorphism.VentoFoggia
import org.rcsb.cif.CifIO
import org.rcsb.cif.model.Column
import org.rcsb.cif.schema.StandardSchemata
import org.rcsb.cif.schema.core.*
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import javax.vecmath.Point3d
import java.nio.file.Path
import java.nio.file.Paths

import static ffx.numerics.math.DoubleMath.dihedralAngle
import static ffx.potential.bonded.BondedUtils.intxyz
import static ffx.potential.parsers.PDBFilter.toPDBAtomLine
import static ffx.potential.parsers.SystemFilter.setVersioning
import static ffx.potential.parsers.SystemFilter.Versioning.PREFIX
import static ffx.potential.parsers.SystemFilter.version
import static java.lang.String.format
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getSymbol

/**
 * The CIFtoXYZ script converts a CIF file to an XYZ file including atom types.
 * TODO: Move CIF parsing into a parsers.CIFFilter class.
 * <br>
 * Usage:
 * <br>
 * ffxc test.CIFtoXYZ &lt;filename.cif&gt; &lt;filename.pdb&gt;
 */
@Command(description = " Convert a single molecule CIF file to XYZ format.", name = "ffxc test.CIFtoXYZ")
class CIFtoXYZ extends PotentialScript {

  /**
   * --sg or --spaceGroupNumber Override the CIF space group.
   */
  @Option(names = ['--sg', '--spaceGroupNumber'], paramLabel = "-1", defaultValue = "-1",
      description = 'Override the CIF space group.')
  private int sgNum

  /**
   * --name or --spaceGroupName Override the CIF space group.
   */
  @Option(names = ['--name', '--spaceGroupName'], paramLabel = "", defaultValue = "",
      description = 'Override the CIF space group.')
  private String sgName

  /**
   * The final argument(s) should be a CIF file and an XYZ file with atom types.
   */
  @Parameters(arity = "1..2", paramLabel = "files",
      description = "A CIF file and an XYZ file (currently limited to a single molecule).")
  List<String> filenames = null

  /**
   * CIFtoXYZ Constructor.
   */
  CIFtoXYZ() {
    this(new Binding())
  }

  /**
   * CIFtoXYZ Constructor.
   * @param binding Groovy Binding to use.
   */
  CIFtoXYZ(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  CIFtoXYZ run() {

    // Turn off CDK logging.
    System.setProperty("cdk.logging.level", "fatal")

    // Turn off non-bonded interactions.
    System.setProperty("vdwterm", "false")

    if (!init()) {
      return this
    }

    //TODO support concatenated CIF files (separated by "data_######").
    CifCoreFile cifFile
    Path path
    if (filenames != null && filenames.size() == 2) {
      path = Paths.get(filenames.get(0))
      cifFile = CifIO.readFromPath(path).as(StandardSchemata.CIF_CORE)
    } else {
      logger.info(helpString())
      return this
    }

    String modelFilename = path.toAbsolutePath().toString()
    logger.info("\n Opening CIF file " + path)

//    CifCoreBlock block = cifFile.firstBlock
    for(CifCoreBlock block: cifFile.blocks) {
      logger.info(format(" Number of Blocks: %d", cifFile.blocks.size()))

      Chemical chemical = block.chemical
      Column nameCommon = chemical.nameCommon
      Column nameSystematic = chemical.nameSystematic
      int rowCount = nameCommon.rowCount
      if (rowCount > 1) {
        logger.info(" Chemical components")
        for (int i = 0; i < rowCount; i++) {
          logger.info(format("  %s", nameCommon.get(i)))
        }
      } else if (rowCount > 0) {
        logger.info(format(" Chemical component: %s", nameCommon.get(0)))
      }

      // Determine the space group.
      Symmetry symmetry = block.symmetry
      if (sgNum == -1 && sgName == "") {
        if (symmetry.intTablesNumber.rowCount > 0) {
          sgNum = symmetry.intTablesNumber.get(0)
          logger.info(format(" CIF International Tables Number: %d", sgNum))
        }
        if (symmetry.spaceGroupNameH_M.rowCount > 0) {
          sgName = symmetry.spaceGroupNameH_M.get(0)
          logger.info(format(" CIF Hermann–Mauguin Space Group: %s", sgName))
        }
      } else {
        if (sgNum != -1) {
          logger.info(format(" Command line International Tables Number: %d", sgNum))
        } else {
          logger.info(format(" Command line space group name: %s", sgName))
        }
      }

      SpaceGroup sg
      sg = SpaceGroupDefinitions.spaceGroupFactory(sgName)
      if (sg == null) {
        sg = SpaceGroupDefinitions.spaceGroupFactory(sgNum)
      }

      // Fall back to P1.
      if (sg == null) {
        logger.warning(" The space group could not be determined from the CIF file (using P1).")
        sg = SpaceGroupDefinitions.spaceGroupFactory(1)
      }

      Cell cell = block.cell
      double a = cell.lengthA.get(0)
      double b = cell.lengthB.get(0)
      double c = cell.lengthC.get(0)
      double alpha = cell.angleAlpha.get(0)
      double beta = cell.angleBeta.get(0)
      double gamma = cell.angleGamma.get(0)

      Crystal crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName)
      logger.info(" New crystal:")
      logger.info(crystal.toString())

      AtomSite atomSite = block.atomSite
      Column label = atomSite.label
      Column typeSymbol = atomSite.typeSymbol
      Column fractX = atomSite.fractX
      Column fractY = atomSite.fractY
      Column fractZ = atomSite.fractZ

      int nAtoms = label.getRowCount()
      logger.info(format("\n Number of atoms: %d", nAtoms))
      Atom[] atoms = new Atom[nAtoms]

      // Define per atom information for the PDB file.
      String resName = "CIF"
      int resID = 1
      double occupancy = 1.0
      double bfactor = 1.0
      char altLoc = ' '
      char chain = 'A'
      String segID = "A"
      String[] symbols = new String[nAtoms]

      // Loop over atoms.
      for (int i = 0; i < nAtoms; i++) {
        symbols[i] = typeSymbol.get(i)
        double x = fractX.get(i)
        double y = fractY.get(i)
        double z = fractZ.get(i)
        double[] xyz = [x, y, z]
        crystal.toCartesianCoordinates(xyz, xyz)
        atoms[i] = new Atom(i + 1, label.get(i), altLoc, xyz, resName, resID, chain, occupancy,
                bfactor, segID)
        atoms[i].setHetero(true)
      }

      // Determine where to save results.
      File saveDir = baseDir
      File cif = new File(modelFilename)
      String cifName = cif.getAbsolutePath()
      if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
        saveDir = new File(FilenameUtils.getFullPath(cifName))
      }
      String dirName = saveDir.toString() + File.separator
      String fileName = FilenameUtils.getName(cifName)

      boolean savePDB = true
      if (filenames.size() > 1) {
        savePDB = false
        // Open the XYZ file (with electrostatics turned off).
        System.setProperty("mpoleterm", "false")
        MolecularAssembly[] assemblies = potentialFunctions.openAll(filenames.get(1))
        System.clearProperty("mpoleterm")

        setActiveAssembly(assemblies[0])
        Atom[] xyzAtoms = activeAssembly.getAtomArray()
        int nXYZAtoms = xyzAtoms.length

        // Used to determine if the CIF structure is missing hydrogen later on.
        int numHydrogens = 0
        for (Atom atom in xyzAtoms) {
          if (atom.isHydrogen()) {
            numHydrogens++
          }
        }

        //TODO: Add >1 ASU
        if (nAtoms == nXYZAtoms || (nAtoms + numHydrogens) == nXYZAtoms) {
          AtomContainer cifCDKAtoms = new AtomContainer()
          AtomContainer xyzCDKAtoms = new AtomContainer()
          for (int i = 0; i < nAtoms; i++) {
            cifCDKAtoms.addAtom(new org.openscience.cdk.Atom(symbols[i],
                    new Point3d(atoms[i].getXYZ(null))))
          }
          AtomContainer nullCDKAtoms = new AtomContainer()

          for (IAtom atom in cifCDKAtoms.atoms()) {
            if (atom.toString() == null) {
              nullCDKAtoms.addAtom(atom)
            }
          }
          cifCDKAtoms.remove(nullCDKAtoms)

          for (Atom atom : xyzAtoms) {
            String atomName = getSymbol(atom.getAtomType().atomicNumber)
            xyzCDKAtoms.addAtom(new org.openscience.cdk.Atom(atomName, new Point3d(atom.getXYZ(null))))
          }

          // Add known XYZ bonds; a limitation is that bonds all are given a Bond order of 1.
          List<Bond> bonds = activeAssembly.getBondList()
          Order order = Order.SINGLE
          int xyzBonds = bonds.size()
          for (Bond bond : bonds) {
            xyzCDKAtoms.addBond(bond.getAtom(0).xyzIndex - 1, bond.getAtom(1).xyzIndex - 1, order)
          }

          // Assign CDK atom types for the XYZ molecule.
          AtomTypeFactory factory = AtomTypeFactory.getInstance(
                  "org/openscience/cdk/config/data/jmol_atomtypes.txt",
                  xyzCDKAtoms.getBuilder())

          for (IAtom atom : xyzCDKAtoms.atoms()) {
            String atomTypeName = atom.getAtomTypeName()
            if (atomTypeName == null || atomTypeName.length() == 0) {
              IAtomType[] types = factory.getAtomTypes(atom.getSymbol())
              if (types.length > 0) {
                IAtomType atomType = types[0]
                atom.setAtomTypeName(atomType.getAtomTypeName())
              } else {
                logger.info(" No atom type found for " + atom.toString())
              }
            }
            factory.configure(atom)
          }

          // Compute bonds for CIF molecule.
          factory = AtomTypeFactory.getInstance(
                  "org/openscience/cdk/config/data/jmol_atomtypes.txt",
                  cifCDKAtoms.getBuilder())

          for (IAtom atom : cifCDKAtoms.atoms()) {
            String atomTypeName = atom.getAtomTypeName()
            if (atomTypeName == null || atomTypeName.length() == 0) {
              IAtomType[] types = factory.getAtomTypes(atom.getSymbol())
              if (types.length > 0) {
                IAtomType atomType = types[0]
                atom.setAtomTypeName(atomType.getAtomTypeName())
              } else {
                logger.info(" No atom type found for " + atom.toString())
              }
            }
            factory.configure(atom)
          }

          double maxCovalentRadius = 2.0
          double minBondDistance = 0.5
          double bondTolerance = 0.5
          RebondTool rebonder = new RebondTool(maxCovalentRadius, minBondDistance, bondTolerance)
          rebonder.rebond(cifCDKAtoms)

          int cifBonds = cifCDKAtoms.bondCount

          // Number of bonds matches.
          if (cifBonds == xyzBonds) {
            Pattern pattern =
                    VentoFoggia.findIdentical(xyzCDKAtoms, AtomMatcher.forElement(), BondMatcher.forAny())
            int[] p = pattern.match(cifCDKAtoms)
            if (p != null && p.length == nAtoms) {
              // Used matched atoms to update the positions of the XYZ file atoms.
              for (int i = 0; i < p.length; i++) {
                // logger.info(format(" %d XYZ %s -> CIF %s", i, xyzCDKAtoms.getAtom(i).getSymbol(), cifCDKAtoms.getAtom(p[i]).getSymbol()))
                Point3d point3d = cifCDKAtoms.getAtom(p[i]).getPoint3d()
                xyzAtoms[i].setXYZ(point3d.x, point3d.y, point3d.z)
              }
            } else {
              logger.info(" CIF and XYZ files don't match.")
              savePDB = true
            }

          } else if ((cifBonds + numHydrogens) % xyzBonds == 0) {
            // Hydrogens most likely missing from file.
            logger.info(" CIF may contain implicit hydrogen -- attempting to patch.")
            // Match heavy atoms between CIF and XYZ
            Pattern pattern = VentoFoggia.
                    findSubstructure(cifCDKAtoms, AtomMatcher.forElement(), BondMatcher.forAny())
            int[] p = pattern.match(xyzCDKAtoms)
            if (p != null && p.length == nAtoms) {
              // Used matched atoms to update the positions of the XYZ file atoms.
              for (int i = 0; i < p.length; i++) {
                Point3d point3d = cifCDKAtoms.getAtom(i).getPoint3d()
                xyzAtoms[p[i]].setXYZ(point3d.x, point3d.y, point3d.z)
              }
              // Add in hydrogen atoms
              Atom lastKnownAtom1 = null
              int knownCount = 0
              for (Atom hydrogen in xyzAtoms) {
                if (hydrogen.isHydrogen()) {
                  Bond bond0 = hydrogen.getBonds()[0]
                  Atom atom1 = bond0.get1_2(hydrogen)
                  double angle0_2 = 0
                  List<Angle> anglesList = hydrogen.getAngles() // Same as bond
                  Atom atom2 = null
                  switch (anglesList.size()) {
                    case 0:
                      // H-Cl No angles
                      logger.info("Structure may contain only two atoms. Not supported yet.")
                      break
                    case 1:
                      // H-O=C Need torsion
                      for (Angle angle in anglesList) {
                        atom2 = angle.get1_3(hydrogen)
                        if (atom2 != null && !atom2.isHydrogen()) {
                          angle0_2 = angle.angleType.angle[0]
                        }
                        if (angle0_2 != 0) {
                          //if angle0_2 is found then no need to look for another atom.
                          break
                        }
                      }
                      List<Bond> bonds2 = atom2.getBonds()
                      Atom atom3 = (atom1 == bonds2[0].get1_2(atom2)) ? bonds2[1].get1_2(atom2) : bonds2[0].get1_2(atom2)
                      double diAng =
                              dihedralAngle(hydrogen.getXYZ(null), atom1.getXYZ(null),
                                      atom2.getXYZ(null), atom3.getXYZ(null))
                      intxyz(hydrogen, atom1, bond0.bondType.distance, atom2, angle0_2, atom3,
                              Math.toDegrees(diAng), 0)
                      break
                    default:
                      // H-C(-C)(-C)-C
                      Atom atom2B = null
                      double angle0_2B = 0
                      Atom proposedAtom
                      double proposedAngle = 0
                      int chiral = 1
                      for (Angle angle in anglesList) {
                        proposedAtom = angle.get1_3(hydrogen)
                          if (proposedAtom != null && !proposedAtom.isHydrogen()) {
                            proposedAngle = angle.angleType.angle[0]
                          }
                          if (proposedAngle != 0) {
                            // If angle1_3 is found then no need to look for another atom.
                            if (angle0_2 != 0) {
                              atom2B = proposedAtom
                              angle0_2B = proposedAngle
                            } else {
                              atom2 = proposedAtom
                              angle0_2 = proposedAngle
                            }
                            proposedAngle = 0.0
                          }
                        if (angle0_2 != 0 && angle0_2B != 0) {
                          break
                        }
                      }
                      if (lastKnownAtom1 == null || lastKnownAtom1 != atom1) {
                        lastKnownAtom1 = atom1
                        knownCount = 0
                      } else {
                        knownCount++
                      }
                      if (angle0_2B == 0.0) {
                        // Hydrogen position depends on other hydrogens, use generic location
                        chiral = 0
                        List<Bond> bonds2 = atom2.getBonds()
                        atom2B = (atom1 == bonds2[0].get1_2(atom2)) ? bonds2[1].get1_2(atom2) : bonds2[0].get1_2(atom2)
                        // Evenly space out hydrogens
                        angle0_2B = 180.0 - 120.0 * knownCount
                      } else if (anglesList.size() == 2) {
                        // Trigonal hydrogen use equipartition between selected atoms
                        chiral = 3
                      } else if (knownCount == 1) {
                        // Already created one hydrogen (chiral = 1), use other.
                        chiral = -1
                      }
                      //TODO decern whether to use chiral = 1 or -1 when angles are to three heavy atoms.
                      intxyz(hydrogen, atom1, bond0.bondType.distance, atom2, angle0_2, atom2B, angle0_2B, chiral)
                  }
                }
              }
            } else {
              logger.info(" Could not match heavy atoms between CIF and XYZ.")
              savePDB = true
            }
          } else {
            logger.info(format(" CIF (%d) and XYZ (%d) have a different number of bonds.", cifBonds,
                    xyzBonds))
            savePDB = true
          }
        } else {
          logger.info(format(" CIF (%d) and XYZ (%d) files have different numbers of atoms.", nAtoms,
                  nXYZAtoms))
          savePDB = true
        }
      }

      setVersioning(PREFIX)
      if (savePDB) {
        // Save a PDB file if there is no XYZ file supplied.
        fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
        File saveFile = version(new File(dirName + fileName))

        // Write out PDB file.
        FileWriter fw = new FileWriter(saveFile, false)
        BufferedWriter bw = new BufferedWriter(fw)
        bw.write(crystal.toCRYST1())
        for (int i = 0; i < nAtoms; i++) {
          bw.write(toPDBAtomLine(atoms[i]))
        }
        bw.write("END\n")
        bw.close()
        logger.info("\n Saved PDB file to " + saveFile.getAbsolutePath())
      } else {
        // Update the active molecular assembly to use the CIF space group and unit cell.
        activeAssembly.getPotentialEnergy().setCrystal(crystal)
        // Choose a location to save the file.
        fileName = FilenameUtils.removeExtension(fileName) + ".xyz"
        File saveFile = version(new File(dirName + fileName))

        // Write out the result.
        potentialFunctions.saveAsXYZ(activeAssembly, saveFile)
        logger.info("\n Saved XYZ file:        " + saveFile.getAbsolutePath())

        // Write out a property file.
        fileName = FilenameUtils.removeExtension(fileName) + ".properties"
        File propertyFile = version(new File(dirName + fileName))
        if (!propertyFile.exists()) {
          FileWriter fw = new FileWriter(propertyFile, false)
          BufferedWriter bw = new BufferedWriter(fw)
          ForceField forceField = activeAssembly.getForceField()
          String parameters = forceField.getString("parameters", "none")
          if (parameters != null && !parameters.equalsIgnoreCase("none")) {
            bw.write(format("parameters %s\n", parameters))
          }
          String patch = forceField.getString("patch", "none")
          if (patch != null && !patch.equalsIgnoreCase("none")) {
            bw.write(format("patch %s\n", patch))
          }
          bw.write(format("spacegroup %s\n", crystal.spaceGroup.shortName))

          bw.close()
          logger.info("\n Saved properties file: " + propertyFile.getAbsolutePath())
        } else {
          logger.info("\n Property file already exists:  " + propertyFile.getAbsolutePath())
        }
      }
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    return new ArrayList<Potential>()
  }
}