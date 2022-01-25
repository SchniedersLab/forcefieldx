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
package ffx.potential.parsers;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.crystal.LatticeSystem;
import ffx.crystal.SpaceGroupDefinitions;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;

import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.graph.rebond.RebondTool;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.VentoFoggia;

import org.rcsb.cif.CifIO;
import org.rcsb.cif.model.Column;
import org.rcsb.cif.model.FloatColumn;
import org.rcsb.cif.schema.core.AtomSite;
import org.rcsb.cif.schema.core.CifCoreBlock;
import org.rcsb.cif.schema.core.CifCoreFile;
import org.rcsb.cif.schema.core.Chemical;
import org.rcsb.cif.schema.core.Symmetry;
import org.rcsb.cif.schema.core.Cell;
import org.rcsb.cif.schema.StandardSchemata;

import javax.vecmath.Point3d;

import static ffx.crystal.SpaceGroupConversions.hrConversion;
import static ffx.numerics.math.DoubleMath.dihedralAngle;
import static ffx.numerics.math.DoubleMath.dist;
import static ffx.potential.bonded.BondedUtils.intxyz;
import static java.lang.String.format;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getSymbol;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getCovalentRadius;

/**
 * The CIFFilter class parses TINKER Cartesian coordinate (*.XYZ) files.
 *
 * @author Aaron J. Nessler
 * @since 1.0
 */
public class CIFFilter{

    private static final Logger logger = Logger.getLogger(CIFFilter.class.getName());

    /**
     * Maximum atomic covalent radius for CDK Rebonder Tool
     */
    private static final double MAX_COVALENT_RADIUS = 2.0;

    /**
     * Minimum bond distance for CDK Rebonder Tool
     */
    private static final double MIN_BOND_DISTANCE = 0.5;

    /**
     * List of output files created from converted CIF file.
     */
    private List<String> createdFiles = new ArrayList<>();
    /**
     * Input filenames to be converted (0: CIF file 1: XYZ file)
     */
    private final String[] fileNames;
    /**
     * Override space group based on number.
     */
    private int sgNum;
    /**
     * Override space group based on name.
     */
    private String sgName;
    /**
     * Tolerance added to covalent radii to determine bonded atoms.
     */
    private final double bondTolerance;
    /**
     * Attempt to fix inconsistencies between lattice parameters and space group.
     */
    private final boolean fixLattice;
    /**
     * Directory where conversion is taking place.
     */
    private final File baseDir;
    /**
     * XYZ assembly to base CIF file off of.
     */
    private final MolecularAssembly activeAssembly;


    public CIFFilter(String[] fileNames, MolecularAssembly activeAssembly, int sgNum,
                     String sgName, double bondTolerance, boolean fixLattice, File baseDir) {
        this.fileNames = fileNames;
        this.sgName = sgName;
        this.sgNum = sgNum;
        this.activeAssembly = activeAssembly;
        this.bondTolerance = bondTolerance;
        this.fixLattice = fixLattice;
        this.baseDir = baseDir;
    }

    /**
     * Open a CIF file and attempt to convert to XYZ.
     * @return Whether conversion was successful.
     */
    public boolean readFile() {

        CifCoreFile cifFile;
        Path path = Paths.get(fileNames[0]);
        try {
            cifFile = CifIO.readFromPath(path).as(StandardSchemata.CIF_CORE);
        } catch (Exception ex) {
            logger.info(" Failed to create CIF file object.\n" + ex + "\n" + Arrays.toString(ex.getStackTrace()));
            return false;
        }

        String modelFilename = path.toAbsolutePath().toString();
        logger.info("\n Opening CIF file " + path);
        int numFailed = 0;
        for (CifCoreBlock block : cifFile.getBlocks()) {
            logger.info(" Block ID: " + block.getBlockHeader());
            Chemical chemical = block.getChemical();
            Column<String[]> nameCommon = chemical.getNameCommon();
            Column<String[]> nameSystematic = chemical.getNameSystematic();
            int rowCount = nameCommon.getRowCount();
            if (rowCount > 1) {
                logger.info(" Chemical components");
                for (int i = 0; i < rowCount; i++) {
                    logger.info(format("  %s", nameCommon.getColumnName()));
                }
            } else if (rowCount > 0) {
                logger.info(format(" Chemical component: %s", nameCommon.getColumnName()));
            }

            // Determine the space group.
            Symmetry symmetry = block.getSymmetry();
            if (sgNum == -1 && sgName.equals("")) {
                if (symmetry.getIntTablesNumber().getRowCount() > 0) {
                    sgNum = symmetry.getIntTablesNumber().get(0);
                    logger.info(format(" CIF International Tables Number: %d", sgNum));
                }
                if (symmetry.getSpaceGroupNameH_M().getRowCount() > 0) {
                    sgName = symmetry.getSpaceGroupNameH_M().get(0);
                    logger.info(format(" CIF Hermann–Mauguin Space Group: %s", sgName));
                } else if (block.getSpaceGroup().getNameH_mFull().getRowCount() > 0) {
                    sgName = block.getSpaceGroup().getNameH_mFull().get(0);
                    logger.info(format(" CIF Hermann–Mauguin Space Group: %s", sgName));
                } else if (block.getSpaceGroup().getNameH_mAlt().getRowCount() > 0) {
                    sgName = block.getSpaceGroup().getNameH_mAlt().get(0);
                    logger.info(format(" CIF Hermann–Mauguin Space Group: %s", sgName));
                }
            } else {
                if (sgNum != -1) {
                    logger.info(format(" Command line International Tables Number: %d", sgNum));
                } else {
                    logger.info(format(" Command line space group name: %s", sgName));
                }
            }

            SpaceGroup sg;
            sg = SpaceGroupDefinitions.spaceGroupFactory(sgName);
            if (sg == null) {
                logger.finer(" Space group name not found. Attempting to use space group number.");
                sg = SpaceGroupDefinitions.spaceGroupFactory(sgNum);
            }

            // Fall back to P1.
            if (sg == null) {
                logger.warning(" The space group could not be determined from the CIF file (using P1).");
                sg = SpaceGroupDefinitions.spaceGroupFactory(1);
                numFailed++;
            }

            Cell cell = block.getCell();
            double a = cell.getLengthA().get(0);
            double b = cell.getLengthB().get(0);
            double c = cell.getLengthC().get(0);
            double alpha = cell.getAngleAlpha().get(0);
            double beta = cell.getAngleBeta().get(0);
            double gamma = cell.getAngleGamma().get(0);
            LatticeSystem latticeSystem = sg.latticeSystem;
            double[] latticeParameters = {a, b, c, alpha, beta, gamma};
            Crystal crystal = null;
            if (latticeSystem.validParameters(a, b, c, alpha, beta, gamma)) {
                crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName);
            } else {
                if (fixLattice) {
                    logger.warning(
                            " Attempting to patch disagreement between lattice system and lattice parameters.");
                    boolean fixed = false;
                    // Check if Rhombohedral lattice has been named Hexagonal
                    if (latticeSystem == LatticeSystem.HEXAGONAL_LATTICE &&
                            LatticeSystem.RHOMBOHEDRAL_LATTICE.validParameters(a, b, c, alpha, beta, gamma)) {
                        crystal = hrConversion(a, b, c, alpha, beta, gamma, sg);
                        latticeSystem = crystal.spaceGroup.latticeSystem;
                        if (
                                latticeSystem
                                        .validParameters(crystal.a, crystal.b, crystal.c, crystal.alpha, crystal.beta,
                                                crystal.gamma)) {
                            fixed = true;
                        }
                        // Check if Hexagonal lattice has been named Rhombohedral
                    } else if (latticeSystem == LatticeSystem.RHOMBOHEDRAL_LATTICE &&
                            LatticeSystem.HEXAGONAL_LATTICE.validParameters(a, b, c, alpha, beta, gamma)) {
                        crystal = hrConversion(a, b, c, alpha, beta, gamma, sg);
                        latticeSystem = crystal.spaceGroup.latticeSystem;
                        if (
                                latticeSystem
                                        .validParameters(crystal.a, crystal.b, crystal.c, crystal.alpha, crystal.beta,
                                                crystal.gamma)) {
                            fixed = true;
                        }
                    }
                    if (!fixed) {
                        double[] newLatticeParameters = latticeSystem.fixParameters(a, b, c, alpha, beta, gamma);
                        if (newLatticeParameters == latticeParameters) {
                            logger.warning(" Conversion Failed: The proposed lattice parameters for " + sg.pdbName
                                    + " do not satisfy the " + latticeSystem + ".");
                            numFailed++;
                            continue;
                        } else {
                            a = newLatticeParameters[0];
                            b = newLatticeParameters[1];
                            c = newLatticeParameters[2];
                            alpha = newLatticeParameters[3];
                            beta = newLatticeParameters[4];
                            gamma = newLatticeParameters[5];
                            crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName);
                        }
                    }
                } else {
                    logger.warning(" Conversion Failed: The proposed lattice parameters for " + sg.pdbName
                            + " do not satisfy the " + latticeSystem + ".");
                    logger.info(" Use \"--fl\" flag to attempt to fix automatically.");
                    numFailed++;
                    continue;
                }
            }
            logger.info(" New Crystal: " + crystal);

            AtomSite atomSite = block.getAtomSite();
            Column<String[]> label = atomSite.getLabel();
            Column<String[]> typeSymbol = atomSite.getTypeSymbol();
            FloatColumn fractX = atomSite.getFractX();
            FloatColumn fractY = atomSite.getFractY();
            FloatColumn fractZ = atomSite.getFractZ();

            int nAtoms = label.getRowCount();
            if (nAtoms < 1) {
                logger.warning(" CIF file did not contain coordinates.");
                numFailed++;
                continue;
            }
            if(logger.isLoggable(Level.FINE)) {
                logger.fine(format("\n Number of Atoms in CIF: %d", nAtoms));
            }
            Atom[] atoms = new Atom[nAtoms];

            // Define per atom information for the PDB file.
            String resName = "CIF";
            double occupancy = 1.0;
            double bfactor = 1.0;
            char altLoc = ' ';
            char chain = 'A';
            String segID = "A";
            String[] symbols = new String[nAtoms];

            // Loop over atoms.
            for (int i = 0; i < nAtoms; i++) {
                // Assigning each atom their own resID prevents comparator from treating them as the same atom.
                if (typeSymbol.getRowCount() > 0) {
                    symbols[i] = typeSymbol.getStringData(i);
                } else {
                    symbols[i] = getAtomElement(label.getStringData(i));
                }
                double x = fractX.get(i);
                double y = fractY.get(i);
                double z = fractZ.get(i);
                double[] xyz = {x, y, z};
                crystal.toCartesianCoordinates(xyz, xyz);
                atoms[i] = new Atom(i + 1, label.getStringData(i), altLoc, xyz, resName, i, chain, occupancy,
                        bfactor, segID);
                atoms[i].setHetero(true);
            }

            // Determine where to save results.
            File saveDir = baseDir;
            File cif = new File(modelFilename);
            String cifName = cif.getAbsolutePath();
            if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
                saveDir = new File(FilenameUtils.getFullPath(cifName));
            }
            String dirName = saveDir + File.separator;
            String fileName = FilenameUtils.getName(cifName);

            // Open the XYZ file (with electrostatics turned off).
            System.setProperty("mpoleterm", "false");
            System.clearProperty("mpoleterm");
            activeAssembly.setName(block.getBlockHeader());
            Atom[] xyzAtoms = activeAssembly.getAtomArray();
            ArrayList<ArrayList<Atom>> xyzatoms = new ArrayList<>();
            int nXYZAtoms = xyzAtoms.length;

            int numHydrogens = 0;
            for (Atom atom : xyzAtoms) {
                if (atom.isHydrogen()) {
                    numHydrogens++;
                }
            }

            List<MSNode> entitiesXYZ = activeAssembly.getAllBondedEntities();
            int numEntitiesXYZ = entitiesXYZ.size();
            if(logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Number of entities in XYZ: %d", numEntitiesXYZ));
            }
            MolecularAssembly outputAssembly = new MolecularAssembly(block.getBlockHeader());
            int zPrime = -1;
            if (nAtoms % nXYZAtoms == 0) {
                zPrime = nAtoms / nXYZAtoms;
            } else if (nAtoms % (nXYZAtoms - numHydrogens) == 0) {
                zPrime = nAtoms / (nXYZAtoms - numHydrogens);
            } else {
                zPrime = 1;
            }

            if (zPrime > 1) {
                logger.info(format(" Detected more than one copy in asymmetric unit of CIF file (Z'=%d)." +
                        " -- attempting to separate.", zPrime));
            }
            for (MSNode mol : entitiesXYZ) {
                xyzatoms.add((ArrayList<Atom>) mol.getAtomList());
            }
            int atomIndex = 1;
            for (int i = 0; i < numEntitiesXYZ; i++) {
                MSNode mol = entitiesXYZ.get(i);
                int numXYZMolAtoms = xyzatoms.get(i).size();
                int numMolHydrogens = 0;
                for (Atom atom : xyzatoms.get(i)) {
                    if (atom.isHydrogen()) {
                        numMolHydrogens++;
                    }
                }
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" Current Entity Number of Atoms: %d (%d + %dH)", numXYZMolAtoms, numXYZMolAtoms - numHydrogens, numHydrogens));
                }

                // Set up XYZ file contents as CDK variable
                AtomContainer xyzCDKAtoms = new AtomContainer();
                for (Atom atom : xyzatoms.get(i)) {
                    String atomName = getSymbol(atom.getAtomType().atomicNumber);
                    xyzCDKAtoms.addAtom(new org.openscience.cdk.Atom(atomName, new Point3d(atom.getXYZ(null))));
                }

                int lessIndex = findMaxLessIndex(xyzatoms, i);
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" Less Index: %d", lessIndex));
                }
                // Add known XYZ bonds; a limitation is that bonds all are given a Bond order of 1.
                List<Bond> bonds = mol.getBondList();
                IBond.Order order = IBond.Order.SINGLE;
                int xyzBonds = bonds.size();
                for (Bond xyzBond : bonds) {
                    if (logger.isLoggable(Level.FINE)) {
                        logger.fine(format(" Bonded atom 1: %d, Bonded atom 2: %d", xyzBond.getAtom(0).getXyzIndex(), xyzBond.getAtom(1).getXyzIndex()));
                    }
                    xyzCDKAtoms.addBond(xyzBond.getAtom(0).getXyzIndex() - lessIndex - 1,
                            xyzBond.getAtom(1).getXyzIndex() - lessIndex - 1, order);
                }

                // Assign CDK atom types for the XYZ molecule.
                AtomTypeFactory factory = AtomTypeFactory.getInstance(
                        "org/openscience/cdk/config/data/jmol_atomtypes.txt",
                        xyzCDKAtoms.getBuilder());

                for (IAtom atom : xyzCDKAtoms.atoms()) {
                    setAtomTypes(factory, atom);
                    try {
                        factory.configure(atom);
                    } catch (Exception ex) {
                        logger.info("Failed to configure atoms from CIF.\n" + ex + "\n" + Arrays.toString(ex.getStackTrace()));
                    }
                }

                ArrayList<ArrayList<Integer>> zindices = new ArrayList<>();
                int counter = 0;
                // Bond atoms in CIF file.
                int cifBonds = bondAtoms(atoms, bondTolerance);
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(format(" Created %d bonds between CIF atoms (%d in xyz).",
                            cifBonds, xyzBonds));
                }
                List<Atom> atomPool = new ArrayList<>(Arrays.asList(atoms));

                try {
                    while (!atomPool.isEmpty()) {
                        ArrayList<Atom> molecule = new ArrayList<>();
                        collectAtoms(atomPool.get(0), molecule);
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
                                    " Molecule %d: %d atoms are ready and %d remain (%d atoms in XYZ, %d atoms in CIF). ",
                                    counter + 1, indices.size(), atomPool.size(), nXYZAtoms, nAtoms));
                        }
                        zindices.add(indices);
                        counter++;
                    }
                } catch (Exception e) {
                    logger.info(e + "\n" + Arrays.toString(e.getStackTrace()));
                    logger.warning(" Failed to separate copies within the asymmetric unit.");
                    numFailed++;
                    continue;
                }
                zPrime = zindices.size();
                // Set up CIF file contents as CDK variable
                AtomContainer[] cifCDKAtomsArr = new AtomContainer[zPrime];
                AtomContainer cifCDKAtoms = new AtomContainer();
                for (int j = 0; j < zPrime; j++) {
                    ArrayList<Integer> currentList = zindices.get(j);
                    int cifMolAtoms = currentList.size();
                    if (logger.isLoggable(Level.FINE)) {
                        logger.fine(format(" CIF atoms in current: %d", cifMolAtoms));
                    }
                    if (cifMolAtoms % numXYZMolAtoms == 0 || cifMolAtoms % (numXYZMolAtoms - numMolHydrogens) == 0) {
                        cifCDKAtomsArr[j] = new AtomContainer();
                        for (Integer integer : currentList) {
                            cifCDKAtomsArr[j].addAtom(new org.openscience.cdk.Atom(symbols[integer - 1],
                                    new Point3d(atoms[integer - 1].getXYZ(null))));
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
                                "org/openscience/cdk/config/data/jmol_atomtypes.txt",
                                cifCDKAtomsArr[j].getBuilder());
                        for (IAtom atom : cifCDKAtomsArr[j].atoms()) {
                            setAtomTypes(factory, atom);
                            try {
                                factory.configure(atom);
                            } catch (Exception ex) {
                                logger.info("Failed to configure CIF atoms.\n" + ex + "\n" + Arrays.toString(ex.getStackTrace()));
                            }
                        }

                        RebondTool rebonder = new RebondTool(MAX_COVALENT_RADIUS, MIN_BOND_DISTANCE, bondTolerance);
                        try {
                            rebonder.rebond(cifCDKAtomsArr[j]);
                        } catch (Exception ex) {
                            logger.info("Failed to rebond CIF atoms.\n" + ex + "\n" + Arrays.toString(ex.getStackTrace()));
                        }

                        int cifMolBonds = cifCDKAtomsArr[j].getBondCount();
                        if (logger.isLoggable(Level.FINE)) {
                            logger.fine(format(" Number of CIF bonds: %d", cifMolBonds));
                        }
                        // Number of bonds matches.
                        // If cifMolBonds == 0 then ion or atom with implicit hydrogens (e.g. water, methane, etc.)
                        if (cifMolBonds != 0 && cifMolBonds % xyzBonds == 0) {
                            org.openscience.cdk.isomorphism.Pattern pattern =
                                    VentoFoggia
                                            .findIdentical(xyzCDKAtoms, AtomMatcher.forElement(), BondMatcher.forAny());
                            int[] p = pattern.match(cifCDKAtomsArr[j]);
                            if (p != null && p.length == numXYZMolAtoms) {
                                // Used matched atoms to update the positions of the XYZ file atoms.
                                for (int k = 0; k < p.length; k++) {
                                    if (logger.isLoggable(Level.FINEST)) {
                                        logger.finest(format(" %d XYZ %s -> CIF %s", k, xyzCDKAtoms.getAtom(k).getSymbol(), cifCDKAtoms.getAtom(p[k]).getSymbol()));
                                    }
                                    Point3d point3d = cifCDKAtomsArr[j].getAtom(p[k]).getPoint3d();
                                    xyzatoms.get(i).get(k).setXYZ(new double[]{point3d.x, point3d.y, point3d.z});
                                }
                            } else {
                                logger.info(format(" Atoms from CIF (%d) and XYZ (%d) structures don't match.", p.length, nAtoms));
                                continue;
                            }
                        } else if ((xyzBonds - numMolHydrogens) == 0 || cifMolBonds % ((xyzBonds - numMolHydrogens)) == 0) {
                            // Hydrogens most likely missing from file. If zero then potentially water.
                            logger.info(" CIF may contain implicit hydrogen -- attempting to patch.");
                            // Match heavy atoms between CIF and XYZ
                            org.openscience.cdk.isomorphism.Pattern pattern = VentoFoggia.
                                    findSubstructure(cifCDKAtomsArr[j], AtomMatcher.forElement(), BondMatcher.forAny());
                            int[] p = pattern.match(xyzCDKAtoms);
                            if (p != null && p.length == numXYZMolAtoms - numMolHydrogens) {
                                // Used matched atoms to update the positions of the XYZ file atoms.
                                for (int k = 0; k < p.length; k++) {
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
                                            case 0:
                                                // H-Cl No angles
                                                // Place hydrogen slightly off center of bonded atom (~1Å away).
                                                hydrogen.moveTo(new double[]{atom1.getX() - 0.6, atom1.getY() - 0.6, atom1.getZ() - 0.6});
                                                break;
                                            case 1:
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
                                                List<Bond> bonds2 = atom2.getBonds();
                                                Atom atom3 = (bonds2.size() > 1 && atom1 == bonds2.get(0).get1_2(atom2)) ? bonds2.get(1).get1_2(atom2) :
                                                        bonds2.get(0).get1_2(atom2);
                                                double diAng =
                                                        dihedralAngle(hydrogen.getXYZ(null), atom1.getXYZ(null),
                                                                atom2.getXYZ(null), atom3.getXYZ(null));
                                                intxyz(hydrogen, atom1, bond0.bondType.distance, atom2, angle0_2, atom3,
                                                        Math.toDegrees(diAng), 0);
                                                break;
                                            default:
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
                                                    // Hydrogen position depends on other hydrogens, use generic location
                                                    chiral = 0;
                                                    bonds2 = atom2.getBonds();
                                                    atom2B = (atom1 == bonds2.get(0).get1_2(atom2)) ? bonds2.get(1).get1_2(atom2) :
                                                            bonds2.get(0).get1_2(atom2);
                                                    // Evenly space out hydrogens
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
                            } else {
                                logger.info(" Could not match heavy atoms between CIF and XYZ.");
                                if (p != null) {
                                    logger.info(
                                            format(" Matched %d atoms out of %d in CIF (%d in XYZ)", p.length, nAtoms,
                                                    nXYZAtoms - numHydrogens));
                                }
                                continue;
                            }
                        } else {
                            logger.info(
                                    format(" CIF (%d) and XYZ ([%d+%dH=]%d) have a different number of bonds.", cifMolBonds,
                                            xyzBonds - numHydrogens, numHydrogens, xyzBonds));
                            continue;
                        }
                        cifCDKAtoms.add(cifCDKAtomsArr[j]);
                        Molecule molecule = new Molecule("Molecule-" + i * zPrime + j);
                        List<Atom> atomList = new ArrayList<>();
                        for (Atom atom : xyzatoms.get(i)) {
                            Atom molAtom = new Atom(atomIndex++, atom.getName(), atom.getAtomType(),
                                    atom.getXYZ(null));
                            atomList.add(molAtom);
                        }
                        List<Bond> bondList = mol.getBondList();
                        for (Bond bond : bondList) {
                            Atom a1 = bond.getAtom(0);
                            Atom a2 = bond.getAtom(1);
                            logger.fine(format(" Bonded atom 1: %d, Bonded atom 2: %d", a1.getXyzIndex(), a2.getXyzIndex()));
                            Atom newA1 = atomList.get(a1.getIndex() - lessIndex - 1);
                            Atom newA2 = atomList.get(a2.getIndex() - lessIndex - 1);
                            Bond bond2 = new Bond(newA1, newA2);
                            bond2.setBondType(bond.getBondType());
                        }
                        for (Atom atom : atomList) {
                            molecule.addMSNode(atom);
                        }
                        outputAssembly.addMSNode(molecule);
                    } else {
                        if(logger.isLoggable(Level.FINE)) {
                            logger.fine(format(" Number of atoms in CIF (%d) molecule do not match XYZ (%d + %dH = %d).",
                                    cifMolAtoms, nXYZAtoms - numHydrogens, numHydrogens, nXYZAtoms));
                        }
                    }
                }
            }
            // If no atoms, then conversion has failed... use active assembly
            if(logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Output Assembly Size: %d", outputAssembly.getAtomList().size()));
            }

            if (outputAssembly.getAtomList().size() < 1) {
                outputAssembly = activeAssembly;
            } else {
                if(logger.isLoggable(Level.FINEST)) {
                    for (Atom atom : outputAssembly.getAtomList()) {
                        logger.finest(atom.toString());
                    }
                }
            }

            File saveFile;
            // Choose a location to save the file.
            fileName = FilenameUtils.removeExtension(fileName);
            if (cifFile.getBlocks().size() > 1) {
                if (zPrime > 1) {
                    fileName += "_z" + zPrime;
                }
                // Concatenated CIF files may contain more than one space group.
                fileName = fileName + "_" + crystal.spaceGroup.shortName.replaceAll("\\/", "") + ".arc";
                saveFile = new File(dirName + fileName);
            } else {
                // If only structure, create new file on each run.
                fileName = fileName + ".xyz";
                saveFile = XYZFilter.version(new File(dirName + fileName));
            }

            // Write out the result.
            // Update the active molecular assembly to use the CIF space group and unit cell.
            outputAssembly.setPotential(activeAssembly.getPotentialEnergy());
            outputAssembly.getPotentialEnergy().setCrystal(crystal);
            XYZFilter xyzFilter = new XYZFilter(saveFile, outputAssembly, null, null);
            xyzFilter.writeFile(saveFile, true);
            logger.info("\n Saved XYZ file:        " + saveFile.getAbsolutePath());

            // Write out a property file.
            fileName = FilenameUtils.removeExtension(fileName) + ".properties";
            File propertyFile = new File(dirName + fileName);
            if (!propertyFile.exists()) {
                try {
                    FileWriter fw = new FileWriter(propertyFile, false);
                    BufferedWriter bw = new BufferedWriter(fw);
                    ForceField forceField = activeAssembly.getForceField();
                    String parameters = forceField.getString("parameters", "none");
                    if (parameters != null && !parameters.equalsIgnoreCase("none")) {
                        bw.write(format("parameters %s\n", parameters));
                    }
                    String forcefieldProperty = forceField.getString("forcefield", "none");
                    if (forcefieldProperty != null && !forcefieldProperty.equalsIgnoreCase("none")) {
                        bw.write(format("forcefield %s\n", forcefieldProperty));
                    }
                    String patch = forceField.getString("patch", "none");
                    if (patch != null && !patch.equalsIgnoreCase("none")) {
                        bw.write(format("patch %s\n", patch));
                    }
                    bw.write(format("spacegroup %s\n", crystal.spaceGroup.shortName));
                    if (zPrime > 1) {
                        bw.write("intermolecular-softcore true");
                    }

                    bw.close();
                } catch (Exception ex) {
                    logger.info("Failed to write files.\n" + ex + "\n" + Arrays.toString(ex.getStackTrace()));
                }
                logger.info("\n Saved properties file: " + propertyFile.getAbsolutePath() + "\n");
            } else {
                logger.info("\n Property file already exists:  " + propertyFile.getAbsolutePath() + "\n");
            }
            if (!createdFiles.contains(saveFile.getName())) {
                createdFiles.add(saveFile.getName());
            }
            //Reset space group for next block
            sgNum = -1;
            sgName = "";
            // All done with this molecule. Therefore, clean up.
            outputAssembly.destroy();
        }
        if (numFailed > 0) {
            logger.info(format(" %d CIF file(s) were not successfully converted.", numFailed));
        }
        return true;
    }

    /**
     * Specify the atom types for atoms created by factory.
     * @param factory Atom generator
     * @param atom Atom of interest
     */
    private void setAtomTypes(AtomTypeFactory factory, IAtom atom) {
        String atomTypeName = atom.getAtomTypeName();
        if (atomTypeName == null || atomTypeName.length() == 0) {
            IAtomType[] types = factory.getAtomTypes(atom.getSymbol());
            if (types.length > 0) {
                IAtomType atomType = types[0];
                atom.setAtomTypeName(atomType.getAtomTypeName());
            } else {
                logger.info(" No atom type found for " + atom);
            }
        }
    }

    /**
     * Find the maximum index that is less than those contained in current model.
     * @param xyzatoms ArrayList containing ArrayList of atoms in each entity
     * @param currentList Entity of currently being used.
     * @return Maximum index that is less than indices in current entity.
     */
    private int findMaxLessIndex(ArrayList<ArrayList<Atom>> xyzatoms, int currentList) {
        // If no less indices are found, want to return 0.
        int lessIndex = 0;
        int minIndex = Integer.MAX_VALUE;
        for(Atom atom: xyzatoms.get(currentList)){
            int temp = atom.getIndex();
            if(temp < minIndex){
                minIndex = temp;
            }
        }
        for(int i = 0; i<xyzatoms.size(); i++){
            if(i != currentList) {
                for (Atom atom : xyzatoms.get(i)) {
                    int temp = atom.getIndex();
                    if(temp < minIndex && temp > lessIndex){
                        lessIndex = temp;
                    }
                }
            }
        }
        return lessIndex;
    }

    /**
     * Obtain a list of output files written from the conversion. Used in file not committed to Git...
     * @return String[] containing output file names.
     */
    public String[] getCreatedFileNames() {
        String[] files = new String[createdFiles.size()];
        return createdFiles.toArray(files);
    }

    /**
     * Add bonds between atoms.
     *
     * @param atoms To potentially be bonded together.
     */
    private static int bondAtoms(Atom[] atoms, double bondTolerance) {
        int bondCount = 0;
        for (int i = 0; i < atoms.length; i++) {
            Atom atom1 = atoms[i];
            String atomIelement = getAtomElement(atom1);
            double radiusI = (getCovalentRadius(atomIelement) == null) ? 0: getCovalentRadius(atomIelement);
            if (radiusI == 0) {
                logger.warning(format(" Atom %d (%s) element (%s) not determined.", i + 1, atom1.getName(), atomIelement));
                return -1;
            }
            double[] xyz = {atom1.getX(), atom1.getY(), atom1.getZ()};
            for (int j = i + 1; j < atoms.length; j++) {
                Atom atom2 = atoms[j];
                String atomJelement = getAtomElement(atom2);
                double radiusJ = (getCovalentRadius(atomJelement) == null) ? 0: getCovalentRadius(atomJelement);
                if (radiusJ == 0) {
                    logger.warning(format(" Atom %d element (%s) not determined.", j + 1, atomJelement));
                    return -1;
                }
                double[] xyz2 = {atom2.getX(), atom2.getY(), atom2.getZ()};
                double length = dist(xyz, xyz2);
                double bondLength = radiusI + radiusJ + bondTolerance;
                if (length < bondLength) {
                    bondCount++;
                    Bond bond = new Bond(atom1, atom2);
                    atom1.setBond(bond);
                    atom2.setBond(bond);
                    logger.finest(format(
                            "Bonded atom %d (%s) with atom %d (%s): bond length (%4.4f Å) < tolerance (%4.4f Å)",
                            i + 1, atomIelement, j + 1, atomJelement, length, bondLength));
                }
            }
        }
        return bondCount;
    }

    /**
     * Parse atom name to determine atomic element.
     *
     * @param atom Atom whose element we desire
     * @return String specifying atom element.
     */
    private static String getAtomElement(Atom atom) {
        return getAtomElement(atom.getName());
    }

    /**
     * Parse atom name to determine atomic element.
     *
     * @param name Name of atom whose element we desire
     * @return String specifying atom element.
     */
    private static String getAtomElement(String name) {
        return name.replaceAll("[()]", "").replaceAll("_", "").replaceAll("-", "").replaceAll(" +", "")
                .split("[0-9]")[0];
    }

    /**
     * Finds all atoms that are bonded to one another.
     * See potential/Utilities/collectAtoms
     *
     * @param seed Starting atom.
     * @param atoms that are bonded.
     */
    private static void collectAtoms(Atom seed, ArrayList<Atom> atoms) {
        logger.finest(format(" Atom: %s", seed.getName()));
        atoms.add(seed);
        for (Bond b : seed.getBonds()) {
            Atom nextAtom = b.get1_2(seed);
            if (!atoms.contains(nextAtom)) {
                collectAtoms(nextAtom, atoms);
            }
        }
    }
}