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
import ffx.crystal.LatticeSystem
import ffx.crystal.SpaceGroup
import ffx.crystal.SpaceGroupDefinitions
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Angle
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Bond
import ffx.potential.bonded.Molecule
import ffx.potential.cli.PotentialScript
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.XYZFilter
import ffx.potential.parsers.PDBFilter

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
import org.rcsb.cif.schema.core.Chemical
import org.rcsb.cif.schema.core.CifCoreBlock
import org.rcsb.cif.schema.core.CifCoreFile
import org.rcsb.cif.schema.core.Cell
import org.rcsb.cif.schema.core.Symmetry
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import javax.vecmath.Point3d
import java.nio.file.Path
import java.nio.file.Paths

import static ffx.crystal.SpaceGroupConversions.hrConversion
import static ffx.numerics.math.DoubleMath.dihedralAngle
import static ffx.numerics.math.DoubleMath.dist
import static ffx.potential.bonded.BondedUtils.intxyz
import static java.lang.String.format
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getSymbol
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getCovalentRadius

/**
 * The CIFtoXYZ script converts a CIF file to an XYZ file including atom types.
 * TODO: Move CIF parsing into a parsers.CIFFilter class.
 * <br>
 * Usage:
 * <br>
 * ffxc CIFtoXYZ &lt;filename.cif&gt; &lt;filename.pdb&gt;
 */
@Command(description = " Convert a single molecule CIF file to XYZ format.", name = "ffxc CIFtoXYZ")
class CIFtoXYZ extends PotentialScript {

    /**
     * --sg or --spaceGroupNumber Override the CIF space group.
     */
    @Option(names = ['--sg', '--spaceGroupNumber'], paramLabel = "-1", defaultValue = "-1",
            description = 'Override the CIF space group.')
    private int sgNum

    /**
     * --fl or --fixLattice Override CIF parameters to satisfy lattice conditions.
     */
    @Option(names = ['--fl', '--fixLattice'], paramLabel = "false", defaultValue = "false",
            description = 'Override CIF parameters to satisfy lattice conditions (Otherwise error).')
    private boolean fixLattice

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
            description = "A CIF file and an XYZ file.")
    List<String> filenames = null

    /**
     * Maximum atomic covalent radius for CDK Rebonder Tool
     */
    private static final double MAX_COVALENT_RADIUS = 2.0

    /**
     * Minimum bond distance for CDK Rebonder Tool
     */
    private static final double MIN_BOND_DISTANCE = 0.5

    /**
     * Atoms within covalent bond radii + tolerance will be bonded.
     */
    private static final double BOND_TOLERANCE = 0.2

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
        int numFailed = 0

        for (CifCoreBlock block : cifFile.blocks) {
            logger.info(" Block ID: " + block.getBlockHeader())
            Chemical chemical = block.chemical
            Column nameCommon = chemical.nameCommon
            Column nameSystematic = chemical.nameSystematic
            int rowCount = nameCommon.getRowCount()
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
                if (symmetry.getIntTablesNumber().getRowCount() > 0) {
                    sgNum = symmetry.getIntTablesNumber().get(0)
                    logger.info(format(" CIF International Tables Number: %d", sgNum))
                }
                if (symmetry.getSpaceGroupNameH_M().getRowCount() > 0) {
                    sgName = symmetry.getSpaceGroupNameH_M().get(0)
                    logger.info(format(" CIF Hermann–Mauguin Space Group: %s", sgName))
                }else if(block.getSpaceGroup().getNameH_MFull().getRowCount() > 0){
                    sgName = block.getSpaceGroup().getNameH_MFull().get(0)
                    logger.info(format(" CIF Hermann–Mauguin Space Group: %s", sgName))
                }else if(block.getSpaceGroup().getNameH_MAlt().getRowCount() > 0){
                    sgName = block.getSpaceGroup().getNameH_MAlt().get(0)
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
                logger.finer(" Space group name not found. Attempting to use space group number.")
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
            LatticeSystem latticeSystem = sg.latticeSystem
            double[] latticeParameters = [a, b, c, alpha, beta, gamma]
            Crystal crystal
            if (latticeSystem.validParameters(a, b, c, alpha, beta, gamma)) {
                crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName)
            } else {
                if (fixLattice) {
                    logger.warning(" Attempting to patch disagreement between lattice system and lattice parameters.")
                    boolean fixed = false
                    if(latticeSystem == LatticeSystem.HEXAGONAL_LATTICE || latticeSystem == LatticeSystem.RHOMBOHEDRAL_LATTICE) {
                        crystal = hrConversion(a, b, c, alpha, beta, gamma, sg)
                        latticeSystem = crystal.spaceGroup.latticeSystem
                        if (latticeSystem.validParameters(crystal.a, crystal.b, crystal.c, crystal.alpha, crystal.beta, crystal.gamma)) {
                            fixed=true
                        }
                    }
                    if(!fixed) {
                        double[] newLatticeParameters = latticeSystem.fixParameters(a, b, c, alpha, beta, gamma)
                        if (newLatticeParameters == latticeParameters) {
                            logger.warning(" Conversion Failed: The proposed lattice parameters for " + sg.pdbName
                                    + " do not satisfy the " + latticeSystem + ".")
                            numFailed++
                            continue
                        } else {
                            a = newLatticeParameters[0]
                            b = newLatticeParameters[1]
                            c = newLatticeParameters[2]
                            alpha = newLatticeParameters[3]
                            beta = newLatticeParameters[4]
                            gamma = newLatticeParameters[5]
                            crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName)
                        }
                    }
                } else {
                    logger.warning(" Conversion Failed: The proposed lattice parameters for " + sg.pdbName
                            + " do not satisfy the " + latticeSystem + ".")
                    logger.info(" Use \"--fl\" flag to attempt to fix automatically.")
                    numFailed++
                    continue
                }
            }
            logger.info(" New Crystal: " + crystal.toString())

            AtomSite atomSite = block.atomSite
            Column label = atomSite.label
            Column typeSymbol = atomSite.typeSymbol
            Column fractX = atomSite.fractX
            Column fractY = atomSite.fractY
            Column fractZ = atomSite.fractZ

            int nAtoms = label.getRowCount()
            if(nAtoms < 1){
                logger.warning(" CIF file did not contain coordinates.")
                numFailed++
                continue
            }
            logger.info(format("\n Number of Atoms: %d", nAtoms))
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
                if(typeSymbol.getRowCount() > 0){
                    symbols[i] = typeSymbol.get(i)
                }else{
                    symbols[i] = getAtomElement(label.get(i))
                }
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

            boolean savePDB = false
            // Open the XYZ file (with electrostatics turned off).
            System.setProperty("mpoleterm", "false")
            MolecularAssembly[] assemblies = potentialFunctions.openAll(filenames.get(1))
            System.clearProperty("mpoleterm")

            setActiveAssembly(assemblies[0])
            activeAssembly.setName(block.getBlockHeader())
            Atom[] xyzAtoms = activeAssembly.getAtomArray()
            int nXYZAtoms = xyzAtoms.length

            // Used to determine if the CIF structure is missing hydrogen later on.
            int numHydrogens = 0
            for (Atom atom in xyzAtoms) {
                if (atom.isHydrogen()) {
                    numHydrogens++
                }
            }

            boolean itsFine = false
            int zPrime
            // Check if there are the same number of atoms in both (CIF may contain multiple copies).
            if (nAtoms % nXYZAtoms == 0) {
                zPrime = (int) (nAtoms / nXYZAtoms)
                itsFine = true
                // Check if there are the same number of heavy atoms in both (CIF may contain multiple copies).
            } else if ((nAtoms + numHydrogens) % nXYZAtoms == 0) {
                zPrime = (int) ((nAtoms + numHydrogens) / nXYZAtoms)
                itsFine = true
            } else {
                zPrime = 1
            }
            if (zPrime > 1) {
                logger.info(format(" Detected more than one copy in asymmetric unit of CIF file (Z'=%d) -- attempting to separate.", zPrime))
            }

            MolecularAssembly outputAssembly = new MolecularAssembly(block.getBlockHeader())

            if (itsFine) {
                // Set up XYZ file contents as CDK variable
                AtomContainer xyzCDKAtoms = new AtomContainer()
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

                int nAUatoms = (int) (nAtoms / zPrime)
                int[][] zIndices = new int[zPrime][nAUatoms]
                int counter = 0
                if (zPrime == 1) {
                    for (int i = 0; i < nAtoms; i++) {
                        zIndices[0][i] = atoms[i].getIndex() - 1
                    }
                } else {
                    // Bond atoms in CIF file.
                    bondAtoms(atoms)

                    List<Atom> atomPool = new ArrayList<>()
                    for (int i = 0; i < atoms.size(); i++) {
                        atomPool.add(atoms[i])
                    }


                    try {
                        while (!atomPool.isEmpty()) {
                            List<Atom> molecule = new ArrayList<>()
                            collectAtoms(atomPool.get(0), molecule)
                            List<Integer> indices = new ArrayList<>()
                            while (molecule.size() > 0) {
                                Atom atom = molecule.remove(0)
                                indices.add(atom.getIndex())
                                atomPool.remove(atom)
                            }
                            if (counter >= zPrime) {
                                logger.info(format("ERROR: Molecule %d of %d detected, but %d atoms are ready and %d remain (%d atoms per molecule, %d atoms in CIF). ",
                                        counter + 1, zPrime, indices.size(), atomPool.size(), nAUatoms, nAtoms))
                            }
                            zIndices[counter++] = indices.stream().mapToInt(i -> i - 1).toArray()
                        }
                    } catch (Exception e) {
                        logger.info(e.toString())
                        logger.warning(" Failed to separate copies within the asymmetric unit.")
                        numFailed++
                        continue
                    }
                }

                // Set up CIF file contents as CDK variable
                AtomContainer[] cifCDKAtomsArr = new AtomContainer[zPrime]
                AtomContainer cifCDKAtoms = new AtomContainer()
                int atomIndex = 1
                for (int i = 0; i < zPrime; i++) {
                    cifCDKAtomsArr[i] = new AtomContainer()
                    for (int j = 0; j < nAUatoms; j++) {
                        cifCDKAtomsArr[i].addAtom(new org.openscience.cdk.Atom(symbols[zIndices[i][j]],
                                new Point3d(atoms[zIndices[i][j]].getXYZ(null))))
                    }
                    AtomContainer nullCDKAtoms = new AtomContainer()

                    for (IAtom atom in cifCDKAtomsArr[i].atoms()) {
                        if (atom.toString() == null) {
                            nullCDKAtoms.addAtom(atom)
                        }
                    }
                    cifCDKAtomsArr[i].remove(nullCDKAtoms)


                    // Compute bonds for CIF molecule.
                    factory = AtomTypeFactory.getInstance(
                            "org/openscience/cdk/config/data/jmol_atomtypes.txt",
                            cifCDKAtomsArr[i].getBuilder())
                    for (IAtom atom : cifCDKAtomsArr[i].atoms()) {
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

                    RebondTool rebonder = new RebondTool(MAX_COVALENT_RADIUS, MIN_BOND_DISTANCE, BOND_TOLERANCE)
                    rebonder.rebond(cifCDKAtomsArr[i])

                    int cifBonds = cifCDKAtomsArr[i].getBondCount();

                    // Number of bonds matches.
                    if (cifBonds % xyzBonds == 0) {
                        Pattern pattern =
                                VentoFoggia.findIdentical(xyzCDKAtoms, AtomMatcher.forElement(), BondMatcher.forAny())
                        int[] p = pattern.match(cifCDKAtomsArr[i])
                        if (p != null && p.length == nAUatoms) {
                            // Used matched atoms to update the positions of the XYZ file atoms.
                            for (int j = 0; j < p.length; j++) {
                                // logger.info(format(" %d XYZ %s -> CIF %s", j, xyzCDKAtoms.getAtom(j).getSymbol(), cifCDKAtoms.getAtom(p[j]).getSymbol()))
                                Point3d point3d = cifCDKAtomsArr[i].getAtom(p[j]).getPoint3d()
                                xyzAtoms[j].setXYZ(point3d.x, point3d.y, point3d.z)
                            }
                        } else {
                            logger.warning(format(" CIF (%d) and XYZ (%d) files don't match.", p.length, nAtoms))
                            savePDB = true
                        }

                    } else if ((cifBonds + numHydrogens) % xyzBonds == 0) {
                        // Hydrogens most likely missing from file.
                        logger.info(" CIF may contain implicit hydrogen -- attempting to patch.")
                        // Match heavy atoms between CIF and XYZ
                        Pattern pattern = VentoFoggia.
                                findSubstructure(cifCDKAtomsArr[i], AtomMatcher.forElement(), BondMatcher.forAny())
                        int[] p = pattern.match(xyzCDKAtoms)
                        if (p != null && p.length == nAUatoms) {
                            // Used matched atoms to update the positions of the XYZ file atoms.
                            for (int j = 0; j < p.length; j++) {
                                Point3d point3d = cifCDKAtomsArr[i].getAtom(j).getPoint3d()
                                xyzAtoms[p[j]].setXYZ(point3d.x, point3d.y, point3d.z)
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
                                            if (angle0_2B == 0.0.toDouble()) {
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
                                            //TODO discern whether to use chiral = 1 or -1 when angles are to three heavy atoms.
                                            intxyz(hydrogen, atom1, bond0.bondType.distance, atom2, angle0_2, atom2B, angle0_2B, chiral)
                                    }
                                }
                            }
                        } else {
                            logger.warning(" Could not match heavy atoms between CIF and XYZ.")
                            savePDB = true
                        }
                    } else {
                        logger.warning(format(" CIF (%d) and XYZ ([%dH+%d=]%d) have a different number of bonds.", cifBonds,
                                numHydrogens, xyzBonds-numHydrogens, xyzBonds))
                        savePDB = true
                    }
                    cifCDKAtoms.add(cifCDKAtomsArr[i])
                    Molecule molecule = new Molecule("Molecule-" + i)
                    List<Atom> atomList = new ArrayList<>()
                    for (int j = 0; j < nXYZAtoms; j++) {
                        Atom atom = xyzAtoms[j]
                        Atom molAtom = new Atom(atomIndex++, atom.getName(), atom.getAtomType(), atom.getXYZ(null))
                        atomList.add(molAtom)
                    }
                    List<Bond> bondList = activeAssembly.getBondList()
                    for (Bond bond : bondList) {
                        Atom a1 = bond.getAtom(0)
                        Atom a2 = bond.getAtom(1)
                        Atom newA1 = atomList.get(a1.getIndex() - 1)
                        Atom newA2 = atomList.get(a2.getIndex() - 1)
                        Bond bond2 = new Bond(newA1, newA2)
                        bond2.setBondType(bond.getBondType())
                    }
                    for (int j = 0; j < atomList.size(); j++) {
                        molecule.addMSNode(atomList.get(j))
                    }
                    outputAssembly.addMSNode(molecule)
                }
            } else {
                logger.warning(format(" CIF (%d) and XYZ (%d) files have different numbers of atoms.", nAtoms,
                        nXYZAtoms))
                savePDB = true
            }
            // If no atoms, then conversion has failed... use active assembly
            if (outputAssembly.getAtomList().size() < 1) {
                outputAssembly = activeAssembly
            }

            File saveFile
            if (savePDB) {
                numFailed++
                // Save a PDB file if there is no XYZ file supplied.
                fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
                saveFile = new File(dirName + fileName)

                // Write out PDB file.
                PDBFilter pdbFilter = new PDBFilter(saveFile, outputAssembly, null, null)
                //PDBs should only be written out when an error occurs... Each individually might make it easier to debug.
                pdbFilter.writeFile(saveFile, false)
                logger.info("\n Saved PDB file to " + saveFile.getAbsolutePath())
            } else {
                // Choose a location to save the file.
                fileName = FilenameUtils.removeExtension(fileName)
                if (cifFile.blocks.size() > 1) {
                    if (zPrime > 1) {
                        fileName += "_z" + zPrime.toString()
                    }
                    // Concatenated CIF files may contain more than one space group.
                    fileName = fileName + "_" + sg.shortName.replaceAll("\\/", "") + ".arc"
                    saveFile = new File(dirName + fileName)
                } else {
                    // If only structure, create new file on each run.
                    fileName = fileName + ".xyz"
                    saveFile = XYZFilter.version(new File(dirName + fileName))
                }

                // Write out the result.
                // Update the active molecular assembly to use the CIF space group and unit cell.
                outputAssembly.setPotential(activeAssembly.getPotentialEnergy())
                outputAssembly.getPotentialEnergy().setCrystal(crystal)
                XYZFilter xyzFilter = new XYZFilter(saveFile, outputAssembly, null, null)
                xyzFilter.writeFile(saveFile, true)
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
                    if (zPrime > 1) {
                        bw.write("intermolecular-softcore true")
                    }

                    bw.close()
                    logger.info("\n Saved properties file: " + propertyFile.getAbsolutePath() + "\n")
                } else {
                    logger.info("\n Property file already exists:  " + propertyFile.getAbsolutePath() + "\n")
                }
            }
            //Reset space group for next block
            sgNum = -1
            sgName = ""
            // All done with this molecule. Therefore clean up.
            outputAssembly.destroy()
        }
        if (numFailed > 0) {
            logger.info(format(" %d CIF files were not successfully converted.", numFailed))
        }
        return this
    }

    /**
     * Add bonds between atoms.
     *
     * @param atoms To potentially be bonded together.
     */
    private static void bondAtoms(Atom[] atoms) {
        for (int i = 0; i < atoms.size(); i++) {
            Atom atom1 = atoms[i]
            String atomIelement = getAtomElement(atom1)
            double radiusI = getCovalentRadius(atomIelement)
            if (radiusI == 0) {
                logger.warning(format(" Atom %d element (%s) not determined.", i + 1, atomIelement))
                return
            }
            double[] xyz = [atom1.getX(), atom1.getY(), atom1.getZ()]
            for (int j = 0; j < atoms.size(); j++) {
                if (i == j) {
                    continue
                }
                Atom atom2 = atoms[j]
                String atomJelement = getAtomElement(atom2)
                double radiusJ = getCovalentRadius(atomJelement)
                if (radiusJ == 0) {
                    logger.warning(format(" Atom %d element (%s) not determined.", j + 1, atomJelement))
                    return
                }
                double[] xyz2 = [atom2.getX(), atom2.getY(), atom2.getZ()]
                double length = dist(xyz, xyz2)
                double bondLength = radiusI + radiusJ + BOND_TOLERANCE
                if (length < bondLength) {
                    Bond bond = new Bond(atom1, atom2)
                    atom1.setBond(bond)
                    atom2.setBond(bond)
                    logger.finer(format("Bonded atom %d (%s) with atom %d (%s): bond length (%4.4f Å) < tolerance (%4.4f Å)",
                            i+1, atomIelement, j+1, atomJelement, length, bondLength))
                }
            }
        }
    }

    /**
     * Parse atom name to determine atomic element.
     *
     * @param atom Atom whose element we desire
     * @return String specifying atom element.
     */
    private static String getAtomElement(Atom atom) {
        return getAtomElement(atom.getName())
    }

    /**
     * Parse atom name to determine atomic element.
     *
     * @param atom Atom whose element we desire
     * @return String specifying atom element.
     */
    private static String getAtomElement(String name) {
        return name.replaceAll("[()]","").replaceAll("_", "").replaceAll("-", "").replaceAll(" +", "").split("[0-9]")[0]
    }

    /**
     * Finds all atoms that are bonded to one another.
     * See potential/Utilities/collectAtoms
     *
     * @param seed Starting atom.
     * @param atoms that are bonded.
     */
    private static void collectAtoms(Atom seed, List<Atom> atoms) {
        if (seed == null) {
            return
        }
        atoms.add(seed)
        for (Bond b : seed.getBonds()) {
            Atom nextAtom = b.get1_2(seed)
            if (!atoms.contains(nextAtom)) {
                collectAtoms(nextAtom, atoms)
            }
        }
    }

    @Override
    List<Potential> getPotentials() {
        return new ArrayList<Potential>()
    }
}