//******************************************************************************
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
//******************************************************************************
package ffx.potential.groovy.test

import javax.vecmath.Point3d

import java.nio.file.Path
import java.nio.file.Paths
import static java.lang.String.format

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
import org.rcsb.cif.schema.core.AtomSite
import org.rcsb.cif.schema.core.Cell
import org.rcsb.cif.schema.core.Chemical
import org.rcsb.cif.schema.core.CifCoreBlock
import org.rcsb.cif.schema.core.CifCoreFile
import org.rcsb.cif.schema.core.Symmetry

import ffx.crystal.Crystal
import ffx.crystal.SpaceGroup
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.ForceField
import static ffx.potential.parsers.PDBFilter.toPDBAtomLine

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The CIFtoXYZ script converts a CIF file to XYZ file with atom types.
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
    private int sgNum = -1

    /**
     * --name or --spaceGroupName Override the CIF space group.
     */
    @Option(names = ['--name', '--spaceGroupName'], paramLabel = "", defaultValue = "",
            description = 'Override the CIF space group.')
    private String sgName = ""

    /**
     * The final argument(s) should be a CIF file and an XYZ file with atom types.
     */
    @Parameters(arity = "2..*", paramLabel = "files",
            description = "A CIF file and an XYZ file (currently limited to one molecule).")
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    CIFtoXYZ run() {

        // Turn off CDK logging.
        System.setProperty("cdk.logging.level", "fatal")

        if (!init()) {
            return null
        }

        CifCoreFile cifFile
        Path path
        if (filenames != null && filenames.size() > 0) {
            path = Paths.get(filenames.get(0))
            cifFile = CifIO.readFromPath(path).as(StandardSchemata.CIF_CORE)
        } else {
            logger.info(helpString())
            return null
        }

        String modelFilename = path.toAbsolutePath().toString()
        logger.info("\n Opening CIF file " + path)

        CifCoreBlock firstBlock = cifFile.firstBlock

        Chemical chemical = firstBlock.chemical
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

        // Determine the sapce group.
        Symmetry symmetry = firstBlock.symmetry
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
        if (sgNum != -1) {
            sg = SpaceGroup.spaceGroupFactory(sgNum)
        } else {
            sg = SpaceGroup.spaceGroupFactory(sgName)
        }

        // Fall to back to P1
        if (sg == null) {
            logger.info(" The space group could not be determined from the CIF file (using P1).")
            sg = SpaceGroup.spaceGroupFactory(1)
        }

        Cell cell = firstBlock.cell
        double a = cell.lengthA.get(0)
        double b = cell.lengthB.get(0)
        double c = cell.lengthC.get(0)
        double alpha = cell.angleAlpha.get(0)
        double beta = cell.angleBeta.get(0)
        double gamma = cell.angleGamma.get(0)

        Crystal crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName)
        logger.info(crystal.toString())

        AtomSite atomSite = firstBlock.atomSite
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
            atoms[i] = new Atom(i + 1, label.get(i), altLoc, xyz, resName, resID, chain, occupancy, bfactor, segID)
            atoms[i].setHetero(true)
        }

        // Determine where to save resutls.
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

            activeAssembly = assemblies[0]
            Atom[] xyzAtoms = activeAssembly.getAtomArray()
            int nXYZAtoms = xyzAtoms.length
            if (nAtoms == nXYZAtoms) {
                AtomContainer cifCDKAtoms = new AtomContainer()
                AtomContainer xyzCDKAtoms = new AtomContainer()
                for (int i = 0; i < nAtoms; i++) {
                    cifCDKAtoms.addAtom(new org.openscience.cdk.Atom(symbols[i],
                            new Point3d(atoms[i].getXYZ(null))))
                    xyzCDKAtoms.addAtom(new org.openscience.cdk.Atom(xyzAtoms[i].getAtomType().name,
                            new Point3d(xyzAtoms[i].getXYZ(null))))
                }

                // Add known XYZ bonds. A limitation is all are given a Bond order of 1.
                ArrayList<Bond> bonds = activeAssembly.getBondList()
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
                        IAtomType[] types = factory.getAtomTypes(atom.getSymbol());
                        if (types.length > 0) {
                            IAtomType atomType = types[0]
                            atom.setAtomTypeName(atomType.getAtomTypeName())
                        } else {
                            logger.info(" No atom type found for " + atom.toString());
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
                        IAtomType[] types = factory.getAtomTypes(atom.getSymbol());
                        if (types.length > 0) {
                            IAtomType atomType = types[0]
                            atom.setAtomTypeName(atomType.getAtomTypeName())
                        } else {
                            logger.info(" No atom type found for " + atom.toString());
                        }
                    }
                    factory.configure(atom)
                }

                RebondTool rebonder = new RebondTool(2.0, 0.5, 0.5)
                rebonder.rebond(cifCDKAtoms)

                int cifBonds = cifCDKAtoms.bondCount
                logger.info(format(" CIF Bond count: %d", cifBonds))

                // Number of bonds matches.
                if (cifBonds == xyzBonds) {
                    Pattern pattern = VentoFoggia.findIdentical(xyzCDKAtoms, AtomMatcher.forElement(), BondMatcher.forAny())
                    int[] p = pattern.match(cifCDKAtoms)
                    if (p != null && p.length == nAtoms) {
                        // Used matched atoms to update the positions of the XYZ file atoms.
                        for (int i = 0; i < p.length; i++) {
                            // logger.info(format(" %d XYZ %s -> CIF %s", i, xyzCDKAtoms.getAtom(i).getSymbol(), cifCDKAtoms.getAtom(p[i]).getSymbol()))
                            Point3d point3d = cifCDKAtoms.getAtom(p[i]).getPoint3d()
                            xyzAtoms[i].setXYZ(point3d.x, point3d.y, point3d.z)
                        }

                        // Update the active molecular assembly to use the CIF space group and unit cell.
                        activeAssembly.getPotentialEnergy().setCrystal(crystal)
                        // Choose a location to save the file.
                        fileName = FilenameUtils.removeExtension(fileName) + ".xyz"
                        File modelFile = new File(dirName + fileName)
                        File saveFile = potentialFunctions.versionFile(modelFile)

                        // Write out the result.
                        potentialFunctions.saveAsXYZ(activeAssembly, saveFile)
                        logger.info("\n Saved XYZ file:        " + saveFile.getAbsolutePath())

                        // Write out a property file.
                        fileName = FilenameUtils.removeExtension(fileName) + ".properties"
                        File propertyFile = new File(dirName + fileName)
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
                    } else {
                        logger.info(" CIF and XYZ files don't match.")
                        savePDB = true
                    }
                } else {
                    logger.info(format(" CIF (%d) and XYZ (%d) have a different number of bonds.", cifBonds, xyzBonds))
                    savePDB = true
                }
            } else {
                logger.info(format(" CIF (%d) and XYZ (%d) files have different numbers of atoms.", nAtoms, nXYZAtoms))
                savePDB = true
            }
        }

        if (savePDB) {
            // Save a PDB file if there is no XYZ file supplied.
            fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
            File modelFile = new File(dirName + fileName)
            File saveFile = potentialFunctions.versionFile(modelFile)

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
        }

        return this
    }

    @Override
    List<Potential> getPotentials() {
        return new ArrayList<Potential>()
    }
}
