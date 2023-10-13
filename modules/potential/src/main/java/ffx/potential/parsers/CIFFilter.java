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
package ffx.potential.parsers;

import ffx.crystal.SpaceGroup;
import ffx.crystal.SpaceGroupDefinitions;
import ffx.crystal.Crystal;
import ffx.crystal.LatticeSystem;
import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MSGroup;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.vecmath.Point3d;

import org.apache.commons.configuration2.CompositeConfiguration;

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.graph.rebond.RebondTool;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.isomorphism.AtomMatcher;
import org.openscience.cdk.isomorphism.BondMatcher;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.isomorphism.Pattern;

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

import static ffx.crystal.SpaceGroupConversions.hrConversion;
import static ffx.numerics.math.DoubleMath.dihedralAngle;
import static ffx.numerics.math.DoubleMath.dist;
import static ffx.potential.bonded.BondedUtils.intxyz;

import static java.lang.String.format;

import static org.apache.commons.io.FilenameUtils.getExtension;
import static org.apache.commons.io.FilenameUtils.getFullPath;
import static org.apache.commons.io.FilenameUtils.getName;
import static org.apache.commons.io.FilenameUtils.removeExtension;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getCovalentRadius;
import static org.openscience.cdk.tools.periodictable.PeriodicTable.getSymbol;

/**
 * The CIFFilter class parses CIF coordinate (*.CIF) files.
 *
 * @author Aaron J. Nessler
 * @since 1.0
 */
public class CIFFilter extends SystemFilter {
  /**
   * A logger for this class.
   */
  private static final Logger logger = Logger.getLogger(CIFFilter.class.getName());
  /**
   * Set bond tolerance (distance to determine if bond is made).
   */
  private double bondTolerance = 0.2;
  /**
   * The CIF file reader.
   */
  private final BufferedReader bufferedReader = null;
  /**
   * The CIF file object.
   */
  private CifCoreFile cifFile;
  /**
   * List of output files created from converted CIF file.
   */
  private final List<String> createdFiles = new ArrayList<>();
  /**
   * The current snapshot.
   */
  private int snapShot;
  /**
   * The space group name.
   */
  private String sgName = null;
  /**
   * The space group number.
   */
  private int sgNum = -1;
  /**
   * The crystallographic information for this system.
   */
  private Crystal crystal;
  /**
   * Whether to fix lattice parameters.
   */
  private boolean fixLattice = false;
  /**
   * Whether to use PDB format.
   */
  private boolean usePDB = false;
  /**
   * The number of entities in the asymmetric unit.
   */
  private int zPrime = -1;
  /**
   * Maximum atomic covalent radius for CDK Rebonder Tool
   */
  public static final double MAX_COVALENT_RADIUS = 2.0;

  /**
   * Minimum bond distance for CDK Rebonder Tool
   */
  public static final double MIN_BOND_DISTANCE = 0.5;

  /**
   * Constructor for CIFFilter on a single file and a single assembly.
   *
   * @param file Input file
   * @param molecularAssembly Active assembly
   * @param forceField Force field for save file.
   * @param properties Properties for save file.
   */
  public CIFFilter(File file, MolecularAssembly molecularAssembly, ForceField forceField,
      CompositeConfiguration properties, boolean saveCIF) {
    super(file, molecularAssembly, forceField, properties);
    this.fileType = Utilities.FileType.CIF;

    if (!saveCIF) {
      String dir = getFullPath(molecularAssembly.getFile().getAbsolutePath()) + File.separator;
      String cifName = removeExtension(file.getName()) + ".cif";
      Path path = Paths.get(dir + cifName);
      try {
        cifFile = CifIO.readFromPath(path).as(StandardSchemata.CIF_CORE);
      } catch (Exception ex) {
        logger.info(" Failed to create CIF file object.\n " + ex);
      }
      String ext = getExtension(file.getName());
      currentFile = new File(removeExtension(currentFile.getAbsolutePath()) + ext);
      if(ext.equalsIgnoreCase("pdb")){
        usePDB = true;
      }
    }
  }

  /**
   * Constructor for CIFFilter on a single file and multiple assemblies.
   *
   * @param file Input file
   * @param molecularAssemblies Active assemblies
   * @param forceField Force field for save file.
   * @param properties Properties for save file.
   */
  public CIFFilter(File file, List<MolecularAssembly> molecularAssemblies, ForceField forceField,
      CompositeConfiguration properties, boolean saveCIF) {
    super(file, molecularAssemblies, forceField, properties);
    this.fileType = Utilities.FileType.CIF;

    if (!saveCIF) {
      String dir =
          getFullPath(molecularAssemblies.get(0).getFile().getAbsolutePath()) + File.separator;
      String cifName = removeExtension(file.getName()) + ".cif";
      Path path = Paths.get(dir + cifName);
      try {
        cifFile = CifIO.readFromPath(path).as(StandardSchemata.CIF_CORE);
      } catch (Exception ex) {
        logger.info(" Failed to create CIF file object.\n" + ex);
      }
      currentFile = new File(removeExtension(currentFile.getAbsolutePath()) + ".xyz");
    }
  }

  /**
   * Constructor for CIFFilter on a multiple files and a single assembly.
   *
   * @param files Input files
   * @param molecularAssembly Active assembly
   * @param forceField Force field for save file.
   * @param properties Properties for save file.
   */
  public CIFFilter(List<File> files, MolecularAssembly molecularAssembly, ForceField forceField,
      CompositeConfiguration properties, boolean saveCIF) {
    super(files, molecularAssembly, forceField, properties);
    this.fileType = Utilities.FileType.CIF;

    if (!saveCIF) {
      String dir = getFullPath(molecularAssembly.getFile().getAbsolutePath()) + File.separator;
      String cifName = removeExtension(files.get(0).getName()) + ".cif";
      Path path = Paths.get(dir + cifName);
      try {
        cifFile = CifIO.readFromPath(path).as(StandardSchemata.CIF_CORE);
      } catch (Exception ex) {
        logger.info(" Failed to create CIF file object.\n" + ex);
      }
      currentFile = new File(removeExtension(currentFile.getAbsolutePath()) + ".xyz");
    }
  }

  /**
   * Add bonds between atoms.
   *
   * @param atoms To potentially be bonded together.
   */
  public static int bondAtoms(Atom[] atoms, double bondTolerance) {
    int bondCount = 0;
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      Atom atom1 = atoms[i];
      String atomIelement = getAtomElement(atom1);
      double radiusI =
          (getCovalentRadius(atomIelement) == null) ? 0 : getCovalentRadius(atomIelement);
      if (radiusI == 0) {
        logger.warning(format(" Atom %d (%s) element (%s) not determined.", i + 1, atom1.getName(),
            atomIelement));
        return -1;
      }
      double[] xyz = {atom1.getX(), atom1.getY(), atom1.getZ()};
      for (int j = i + 1; j < nAtoms; j++) {
        Atom atom2 = atoms[j];
        String atomJelement = getAtomElement(atom2);
        double radiusJ =
            (getCovalentRadius(atomJelement) == null) ? 0 : getCovalentRadius(atomJelement);
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
              " Bonded atom %d (%s) with atom %d (%s): bond length (%4.4f Å) < tolerance (%4.4f Å)",
              i + 1, atomIelement, j + 1, atomJelement, length, bondLength));
        }
      }
    }
    return bondCount;
  }

  /** {@inheritDoc} */
  @Override
  public void closeReader() {
    if (bufferedReader != null) {
      try {
        bufferedReader.close();
      } catch (IOException ex) {
        logger.warning(format(" Exception in closing CIF filter: %s", ex));
      }
    }
  }

  /**
   * Finds all atoms that are bonded to one another. See potential/Utilities/collectAtoms
   *
   * @param seed Starting atom.
   * @param atoms that are bonded.
   */
  public static void collectAtoms(Atom seed, ArrayList<Atom> atoms) {
    logger.finest(format(" Atom: %s", seed.getName()));
    atoms.add(seed);
    for (Bond b : seed.getBonds()) {
      Atom nextAtom = b.get1_2(seed);
      if (!atoms.contains(nextAtom)) {
        collectAtoms(nextAtom, atoms);
      }
    }
  }

  /**
   * Find the maximum index that is less than those contained in current model.
   *
   * @param xyzAtomLists ArrayList containing ArrayList of atoms in each entity
   * @param currentList Entity of currently being used.
   * @return Maximum index that is less than indices in current entity.
   */
  public int findMaxLessIndex(ArrayList<ArrayList<Atom>> xyzAtomLists, int currentList) {
    // If no less indices are found, want to return 0.
    int lessIndex = 0;
    int minIndex = Integer.MAX_VALUE;
    for (Atom atom : xyzAtomLists.get(currentList)) {
      int temp = atom.getIndex();
      if (temp < minIndex) {
        minIndex = temp;
      }
    }
    int xyzSize = xyzAtomLists.size();
    for (int i = 0; i < xyzSize; i++) {
      if (i != currentList) {
        for (Atom atom : xyzAtomLists.get(i)) {
          int temp = atom.getIndex();
          if (temp < minIndex && temp > lessIndex) {
            lessIndex = temp;
          }
        }
      }
    }
    return lessIndex;
  }

  /**
   * Parse atom name to determine atomic element.
   *
   * @param atom Atom whose element we desire
   * @return String specifying atom element.
   */
  public static String getAtomElement(Atom atom) {
    return getAtomElement(atom.getName());
  }

  /**
   * Parse atom name to determine atomic element.
   *
   * @param name Name of atom whose element we desire
   * @return String specifying atom element.
   */
  private static String getAtomElement(String name) {
    String value = "";
    try {
      value = name.replaceAll("[()]", "").replaceAll("_", "").replaceAll("-", "")
          .replaceAll(" +", "").split("[0-9]")[0];
    } catch (Exception e) {
      logger.severe(" Error extracting atom element. Please ensure the CIF is formatted correctly.\n" + e);
    }
    return value;
  }

  /**
   * Obtain a list of output files written from the conversion. Used in file not committed to Git...
   *
   * @return String[] containing output file names.
   */
  public String[] getCreatedFileNames() {
    String[] files = new String[createdFiles.size()];
    return createdFiles.toArray(files);
  }

  /**
   * Parse CIF block to determine crystal and atom information.
   * @param block CifCoreBlock object to be interpreted.
   * @return Atoms from CIF. (Crystal updated as global variable)
   */
  private Atom[] parseBlock(CifCoreBlock block){
    // Parse CIF file "blocks"
    logger.info("\n Block ID: " + block.getBlockHeader());
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

    // Determine the space group from CIF.
    Symmetry symmetry = block.getSymmetry();
    if (sgNum == -1 && sgName == null || sgName.equals("")) {
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
      // Parse user defined space group.
      if (sgNum != -1) {
        logger.info(format(" Command line International Tables Number: %d", sgNum));
      } else {
        logger.info(format(" Command line space group name: %s", sgName));
      }
    }

    // Set space group.
    SpaceGroup sg = null;
    if (sgName != null) {
      sg = SpaceGroupDefinitions.spaceGroupFactory(sgName);
    }
    if (sg == null) {
      logger.finer(" Space group name not found. Attempting to use space group number.");
      sg = SpaceGroupDefinitions.spaceGroupFactory(sgNum);
    }

    // Fall back to P1.
    if (sg == null) {
      logger.warning(" The space group could not be determined from the CIF file (using P1).");
      sg = SpaceGroupDefinitions.spaceGroupFactory(1);
    }

    // Generate FFX Crystal for CIF file.
    Cell cell = block.getCell();
    if(cell.isDefined()) {
      double a = cell.getLengthA().get(0);
      double b = cell.getLengthB().get(0);
      double c = cell.getLengthC().get(0);
      double alpha = cell.getAngleAlpha().get(0);
      double beta = cell.getAngleBeta().get(0);
      double gamma = cell.getAngleGamma().get(0);
      assert sg != null;
      LatticeSystem latticeSystem = sg.latticeSystem;
      double[] latticeParameters = {a, b, c, alpha, beta, gamma};
      if (latticeSystem.validParameters(a, b, c, alpha, beta, gamma)) {
        crystal = new Crystal(a, b, c, alpha, beta, gamma, sg.pdbName);
      } else {
        // Fix lattice parameters if they disobey space group.
        if (fixLattice) {
          logger.info(
                  " Attempting to patch disagreement between lattice system and lattice parameters.");
          boolean fixed = false;
          // Check if Rhombohedral lattice has been named Hexagonal
          if (latticeSystem == LatticeSystem.HEXAGONAL_LATTICE
                  && LatticeSystem.RHOMBOHEDRAL_LATTICE.validParameters(a, b, c, alpha, beta, gamma)) {
            crystal = hrConversion(a, b, c, alpha, beta, gamma, sg);
            latticeSystem = crystal.spaceGroup.latticeSystem;
            if (latticeSystem.validParameters(crystal.a, crystal.b, crystal.c, crystal.alpha,
                    crystal.beta, crystal.gamma)) {
              fixed = true;
            }
            // Check if Hexagonal lattice has been named Rhombohedral
          } else if (latticeSystem == LatticeSystem.RHOMBOHEDRAL_LATTICE
                  && LatticeSystem.HEXAGONAL_LATTICE.validParameters(a, b, c, alpha, beta, gamma)) {
            crystal = hrConversion(a, b, c, alpha, beta, gamma, sg);
            latticeSystem = crystal.spaceGroup.latticeSystem;
            if (latticeSystem.validParameters(crystal.a, crystal.b, crystal.c, crystal.alpha,
                    crystal.beta, crystal.gamma)) {
              fixed = true;
            }
          }
          if (!fixed) {
            double[] newLatticeParameters = latticeSystem.fixParameters(a, b, c, alpha, beta, gamma);
            if (newLatticeParameters == latticeParameters) {
              logger.warning(" Conversion Failed: The proposed lattice parameters for " + sg.pdbName
                      + " do not satisfy the " + latticeSystem + ".");
              return null;
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
          logger.info(" Use \"--fixLattice\" or \"--fl\" flag to attempt to fix automatically.");
          return null;
        }
      }
      activeMolecularAssembly.setName(block.getBlockHeader());
      logger.info(" New Crystal: " + crystal);
      activeMolecularAssembly.setCrystal(crystal);
    }
    // Parse CIF Atoms
    AtomSite atomSite = block.getAtomSite();
    Column<String[]> label = atomSite.getLabel();
    Column<String[]> typeSymbol = atomSite.getTypeSymbol();
    FloatColumn fractX = atomSite.getFractX();
    FloatColumn fractY = atomSite.getFractY();
    FloatColumn fractZ = atomSite.getFractZ();

    FloatColumn cartX = atomSite.getCartnX();
    FloatColumn cartY = atomSite.getCartnY();
    FloatColumn cartZ = atomSite.getCartnZ();

    int nAtoms = label.getRowCount();
    if (nAtoms < 1) {
      logger.warning(" CIF file did not contain coordinates.");
      return null;
    }
    if (logger.isLoggable(Level.FINE)) {
      logger.fine(format("\n Number of Atoms in CIF: %d", nAtoms));
    }
    Atom[] atoms = new Atom[nAtoms];

    // Define per atom information from the CIF.
    String resName = "CIF";
    double bfactor = 1.0;
    char altLoc = ' ';
    char chain = 'A';
//    String segID = "A";
    String[] symbols = new String[nAtoms];

    // Loop over atoms in CIF and generate FFX Atom(s).
    for (int i = 0; i < nAtoms; i++) {
      double occupancy = 1.0;
      if(atomSite.getOccupancy().isDefined()){
        occupancy = atomSite.getOccupancy().get(i);
      }
      // Assigning each atom their own resID prevents comparator from treating them as the same atom.
      if (typeSymbol.getRowCount() > 0) {
        symbols[i] = typeSymbol.getStringData(i);
      } else {
        symbols[i] = getAtomElement(label.getStringData(i));
      }
      double x = (fractX.isDefined()) ? fractX.get(i) : cartX.get(i);
      double y = (fractY.isDefined()) ? fractY.get(i) : cartY.get(i);
      double z = (fractZ.isDefined()) ? fractZ.get(i) : cartZ.get(i);
      double[] xyz = {x, y, z};
      if (fractX.isDefined() && crystal != null) {
        crystal.toCartesianCoordinates(xyz, xyz);
      }
      atoms[i] = new Atom(i + 1, label.getStringData(i), altLoc, xyz, resName, i, chain, occupancy,
              bfactor, symbols[i]);
      atoms[i].setHetero(true);
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format(
                " Atom (%2d) Name: " + atoms[i].getName() + " Label: " + label.getStringData(i)
                        + " Symbol: " + symbols[i], i));
      }
    }
    return atoms;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Parse the CIF File
   */
  @Override
  public boolean readFile() {
    // Open the input file (with electrostatics turned off).
    System.setProperty("mpoleterm", "false");
    Atom[] inputAtoms = activeMolecularAssembly.getAtomArray();
    int nInputAtoms = inputAtoms.length;

    int numHydrogen = 0;
    for (Atom atom : inputAtoms) {
      if (atom.isHydrogen()) {
        numHydrogen++;
      }
    }

    List<MSNode> entitiesInput = activeMolecularAssembly.getNodeList();
    int numEntitiesInput = entitiesInput.size();
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

    for (CifCoreBlock block : cifFile.getBlocks()) {
      Atom[] atoms = parseBlock(block);
      logger.info(crystal.toString());
      if(atoms == null){
        return false;
      }
      if(atoms.length == 0){
        logger.warning(" Atoms could not be obtained from CIF.");
        return false;
      }
      int nAtoms = atoms.length;

      // Define assembly object to contain CIF info.
      MolecularAssembly outputAssembly = new MolecularAssembly(block.getBlockHeader());
      int zPrime;
      if (this.zPrime > 0) {
        zPrime = this.zPrime;
      } else if (nAtoms % nInputAtoms == 0) {
        zPrime = nAtoms / nInputAtoms;
      } else if (nAtoms % (nInputAtoms - numHydrogen) == 0) {
        zPrime = nAtoms / (nInputAtoms - numHydrogen);
      } else {
        zPrime = 1;
      }

      if (zPrime > 1) {
        logger.info(format(" Detected more than one copy in asymmetric unit of CIF file (Z'=%d)."
            + " -- attempting to separate.", zPrime));
      }
      ArrayList<ArrayList<Atom>> xyzatoms = new ArrayList<>();
      int atomIndex = 1;
      int shiftIndex = 0;
      int[] moleculeNumbers = activeMolecularAssembly.getMoleculeNumbers();
      if(logger.isLoggable(Level.FINER)) {
        for (int i = 0; i < moleculeNumbers.length; i++) {
          logger.finer(format(" Molecule Number (atom %2d): %2d", i, moleculeNumbers[i]));
        }
      }
      // For each entity (protein, molecule, ion, etc) in the XYZ
      for (int i = 0; i < numEntitiesInput; i++) {
        MSNode mol = entitiesInput.get(i);
        logger.info(format(" Size of Entity: %2d", mol.getAtomList().size()));
        xyzatoms.add((ArrayList<Atom>) mol.getAtomList());
        // Set up input (XYZ) file contents as CDK variable
        AtomContainer xyzCDKAtoms = new AtomContainer();
        // Number of hydrogen in an entity from the XYZ.
        int numMolHydrogen = 0;
        int atomCount = 0;
        for (Atom atom : xyzatoms.get(i)) {
          if (atom.isHydrogen()) {
            numMolHydrogen++;
          }
          String atomName = getSymbol(atom.getAtomType().atomicNumber);
          xyzCDKAtoms.addAtom(
              new org.openscience.cdk.Atom(atomName, new Point3d(atom.getXYZ(null))));
          atomCount++;
        }
        logger.info(format(" Number of atoms in XYZ CDK: %2d", atomCount));
        int numInputMolAtoms = xyzatoms.get(i).size();
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(format(" Current entity number of atoms: %d (%d + %dH)", numInputMolAtoms,
                  numInputMolAtoms - numMolHydrogen, numMolHydrogen));
        }

        // Add known input bonds; a limitation is that bonds all are given a Bond order of 1.
        List<Bond> bonds = mol.getBondList();
        IBond.Order order = IBond.Order.SINGLE;
        int xyzBonds = bonds.size();
        if (xyzBonds == 0) {
          logger.warning(" No bonds detected in input structure. Please check input.\n " +
              "If correct, separate non-bonded entities into multiple CIFs.");
        }
        for (Bond xyzBond : bonds) {
          int atom0Index = xyzBond.getAtom(0).getXyzIndex();
          int atom1Index = xyzBond.getAtom(1).getXyzIndex();
          if (logger.isLoggable(Level.FINER)) {
            logger.finer(
                format(" Bonded atom 1: %d, Bonded atom 2: %d", atom0Index,
                    atom1Index));
          }
          logger.fine(format(" Atom 0 index: %2d Atom 1 Index: %2d lessIndex: %2d Value 0: %2d Value 1: %2d",
                  atom0Index, atom1Index, shiftIndex, atom0Index - shiftIndex - 1, atom1Index - shiftIndex - 1));
          xyzCDKAtoms.addBond(atom0Index - shiftIndex - 1,
              atom1Index - shiftIndex - 1, order);
        }

        // Assign CDK atom types for the input molecule.
        AtomTypeFactory factory = AtomTypeFactory.getInstance(
            "org/openscience/cdk/config/data/jmol_atomtypes.txt", xyzCDKAtoms.getBuilder());

        // Assign atom types to CDK object.
        for (IAtom atom : xyzCDKAtoms.atoms()) {
          setAtomTypes(factory, atom);
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
        int cifBonds = bondAtoms(atoms, bondTolerance);
        if (logger.isLoggable(Level.FINE)) {
          logger.fine(
              format(" Created %d bonds between CIF atoms (%d in input).", cifBonds, xyzBonds));
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
                  " Molecule %d: %d atoms are ready and %d remain (%d atoms in input, %d atoms in CIF). ",
                  counter + 1, indices.size(), atomPool.size(), nInputAtoms, nAtoms));
            }
            zindices.add(indices);
            counter++;
          }
        } catch (Exception e) {
          logger.info(e + "\n" + Arrays.toString(e.getStackTrace()));
          logger.warning(" Failed to separate copies within the asymmetric unit.");
          return false;
        }
        zPrime = zindices.size();
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
          if (cifMolAtoms % numInputMolAtoms == 0
              || cifMolAtoms % (numInputMolAtoms - numMolHydrogen) == 0) {
            cifCDKAtomsArr[j] = new AtomContainer();
            for (Integer integer : currentList) {
              cifCDKAtomsArr[j].addAtom(new org.openscience.cdk.Atom(atoms[integer-1].getSegID(),
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
                logger.info(" Failed to configure CIF atoms.\n" + ex + "\n" + Arrays.toString(
                    ex.getStackTrace()));
              }
            }

            RebondTool rebonder = new RebondTool(MAX_COVALENT_RADIUS, MIN_BOND_DISTANCE,
                bondTolerance);
            try {
              rebonder.rebond(cifCDKAtomsArr[j]);
            } catch (Exception ex) {
              logger.info(
                  "Failed to rebond CIF atoms.\n" + ex + "\n" + Arrays.toString(ex.getStackTrace()));
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
                // Used matched atoms to update the positions of the input file atoms.
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
            } else if ((xyzBonds - numMolHydrogen) == 0
                || cifMolBonds % ((xyzBonds - numMolHydrogen)) == 0) {
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
              if(logger.isLoggable(Level.FINER)) {
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
                if(atom.getResidueName() != null){
                  molAtom.setResName(atom.getResidueName());
                }
              }
              atomList.add(molAtom);
            }
            List<Bond> bondList = mol.getBondList();
            for (Bond bond : bondList) {
              Atom a1 = bond.getAtom(0);
              Atom a2 = bond.getAtom(1);
              if(logger.isLoggable(Level.FINE)) {
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
              if(molecule instanceof Polymer){
                if(res == null){
                  res = new Residue(atom.getResidueName(), atom.getResidueNumber(), ((Polymer) mol).getResidue(atom.getResidueNumber()).getResidueType());
                } else if(res.getResidueNumber() != atom.getResidueNumber()){
                  if(logger.isLoggable(Level.FINER)) {
                    logger.finer(" Added Residue " + res.getResidueNumber() + " " + res.getName());
                  }
                  molecule.addMSNode(res);
                  res = new Residue(atom.getResidueName(), atom.getResidueNumber(), ((Polymer) mol).getResidue(atom.getResidueNumber()).getResidueType());
                }
                res.addMSNode(atom);
              }else {
                molecule.addMSNode(atom);
              }
            }
            if(molecule instanceof Polymer && res != null){
              if(logger.isLoggable(Level.FINER)) {
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
        shiftIndex += mol.getAtomList().size();
      }
      // If no atoms, then conversion has failed... use active assembly
      if (logger.isLoggable(Level.FINE)) {
        logger.fine(format("\n Output Assembly Atoms: %d", outputAssembly.getAtomList().size()));
      }
      outputAssembly.setPotential(activeMolecularAssembly.getPotentialEnergy());
      outputAssembly.setCrystal(crystal);
      outputAssembly.setForceField(activeMolecularAssembly.getForceField());
      outputAssembly.setFile(activeMolecularAssembly.getFile());
      outputAssembly.setName(activeMolecularAssembly.getName());
      setMolecularSystem(outputAssembly);

      if (outputAssembly.getAtomList().size() < 1) {
        logger.info(" Atom types could not be matched. File could not be written.");
      } else if (!writeOutputFile()) {
        logger.info(" Input File could not be written.");
      }
      // Finsished with this block. Clean up/reset variables for next block.
      outputAssembly.destroy();
      sgNum = -1;
      sgName = null;
    }
    return true;
  }

  /**
   * Reads the next snapshot of an archive into the activeMolecularAssembly. After calling this
   * function, a BufferedReader will remain open until the <code>close</code> method is called.
   */
  @Override
  public boolean readNext() {
    return readNext(false, true);
  }

  /**
   * Reads the next snapshot of an archive into the activeMolecularAssembly. After calling this
   * function, a BufferedReader will remain open until the <code>close</code> method is called.
   */
  @Override
  public boolean readNext(boolean resetPosition) {
    return readNext(resetPosition, true);
  }

  /**
   * Reads the next snapshot of an archive into the activeMolecularAssembly. After calling this
   * function, a BufferedReader will remain open until the <code>close</code> method is called.
   */
  @Override
  public boolean readNext(boolean resetPosition, boolean print) {
    return readNext(resetPosition, print, true);
  }

  /**
   * Reads the next snapshot of an archive into the activeMolecularAssembly. After calling this
   * function, a BufferedReader will remain open until the <code>close</code> method is called.
   */
  @Override
  public boolean readNext(boolean resetPosition, boolean print, boolean parse) {
    List<CifCoreBlock> blocks = cifFile.getBlocks();
    CifCoreBlock currentBlock;
    if (!parse) {
      snapShot++;
      if (print) {
        logger.info(format(" Skipped Block: %d", snapShot));
      }
      return true;
    } else if (resetPosition) {
      currentBlock = blocks.get(0);
      snapShot = 0;
    } else if (++snapShot < blocks.size()) {
      currentBlock = blocks.get(snapShot);
    } else {
      if (print) {
        logger.info(" Reached end of available blocks in CIF file.");
      }
      return false;
    }
    if (print) {
      logger.info(" Current Block: " + currentBlock.getBlockHeader());
    }
    return true;
  }

  /**
   * Specify the atom types for atoms created by factory.
   *
   * @param factory Atom generator
   * @param atom Atom of interest
   */
  public static void setAtomTypes(AtomTypeFactory factory, IAtom atom) {
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
   * Set buffer value when bonding atoms together.
   *
   * @param bondTolerance Value added to VdW radii when determining bonding.
   */
  public void setBondTolerance(double bondTolerance) {
    this.bondTolerance = bondTolerance;
  }

  /**
   * Determine whether lattice parameters can be manipulated to follow lattice system constraints.
   *
   * @param fixLattice True if lattices should be fixed.
   */
  public void setFixLattice(boolean fixLattice) {
    this.fixLattice = fixLattice;
  }

  /**
   * Override the space group of a CIF conversion based on space group name.
   *
   * @param sgName Name of the desired space group
   */
  public void setSgName(String sgName) {
    this.sgName = sgName;
  }

  /**
   * Override the space group of a CIF conversion based on space group number.
   *
   * @param sgNum Number of the desired space group
   */
  public void setSgNum(int sgNum) {
    this.sgNum = sgNum;
  }

  /**
   * Override the number of copies in the asymmetric unit (Z').
   *
   * @param zPrime Number of the desired space group
   */
  public void setZprime(int zPrime) {
    this.zPrime = zPrime;
  }

  /**
   * Write CIF files for multiple File objects.
   *
   * @return whether conversion was successful.
   */
  public boolean writeFiles() {
    for (File file : files) {
      File saveFile = new File(removeExtension(file.getAbsolutePath()) + ".cif");
      if (!writeFile(saveFile, true, null)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Save a CIF file for the given molecular assembly.
   *
   * @return whether conversion was successful.
   */
  @Override
  public boolean writeFile(File saveFile, boolean append) {
    return writeFile(saveFile, append, null);
  }

  /**
   * Save a CIF file for the given molecular assembly.
   *
   * @return whether conversion was successful.
   */
  @Override
  public boolean writeFile(File saveFile, boolean append, String[] extraLines) {
    try {
      if (!append) {
        SystemFilter.setVersioning(Versioning.PREFIX);
        saveFile = version(saveFile);
      }
      BufferedWriter bw = new BufferedWriter(new FileWriter(saveFile, append));
      if (extraLines != null) {
        for (String line : extraLines) {
          line = line.replaceAll("\n", " ");
          bw.write("### " + line + "\n");
        }
      }
      bw.write("\ndata_" + activeMolecularAssembly.getName());
      Crystal xtal = activeMolecularAssembly.getCrystal().getUnitCell();
      bw.write("\n_symmetry_cell_setting\t" + xtal.spaceGroup.latticeSystem.name().toLowerCase()
          .replaceAll("_lattice", ""));
      bw.write("\n_symmetry_space_group_name_H-M\t" + "'" + xtal.spaceGroup.shortName + "'");
      bw.write("\n_symmetry_Int_Tables_number\t" + xtal.spaceGroup.number);
      bw.write("\nloop_\n_symmetry_equiv_pos_site_id\n_symmetry_equiv_pos_as_xyz");
      int numSymOps = xtal.spaceGroup.getNumberOfSymOps();
      for (int i = 0; i < numSymOps; i++) {
        bw.write(format(
            "\n%d " + xtal.spaceGroup.getSymOp(i).toXYZString().toLowerCase().replaceAll(" +", ""),
            i));
      }
      bw.write(format("\n_cell_length_a\t%4.4f", xtal.a));
      bw.write(format("\n_cell_length_b\t%4.4f", xtal.b));
      bw.write(format("\n_cell_length_c\t%4.4f", xtal.c));
      bw.write(format("\n_cell_angle_alpha\t%4.4f", xtal.alpha));
      bw.write(format("\n_cell_angle_beta\t%4.4f", xtal.beta));
      bw.write(format("\n_cell_angle_gamma\t%4.4f", xtal.gamma));
      bw.write(format("\n_cell_volume\t%4.4f", xtal.volume));
      int numEntities = (zPrime < 1) ? activeMolecularAssembly.getMolecules().size() : zPrime;
      if (numEntities > 1) {
        if (zPrime < 1) {
          logger.info(format(" Molecules detected, guessing a Z' of %d. Set manually using --zp.",
              numEntities));
        }
        bw.write(format("\n_cell_formula_units_Z\t%3d", numEntities));
      }
      bw.write("\nloop_");
      bw.write("\n_atom_site_label");
      bw.write("\n_atom_site_type_symbol");
      bw.write("\n_atom_site_fract_x");
      bw.write("\n_atom_site_fract_y");
      bw.write("\n_atom_site_fract_z");
      Atom[] atoms = activeMolecularAssembly.getAtomArray();
      int count = 1;
      for (Atom atom : atoms) {
        String name = atom.getName();
        String symbol = getSymbol(atom.getAtomicNumber());
        if (Objects.equals(name, symbol)) {
          name += count++;
        }
        double[] xyzC = {atom.getX(), atom.getY(), atom.getZ()};
        double[] xyzF = new double[3];
        xtal.toFractionalCoordinates(xyzC, xyzF);
        bw.write(format("\n%-3s %2s %8.6f %8.6f %8.6f", name, symbol, xyzF[0], xyzF[1], xyzF[2]));
      }
      bw.write("\n#END\n");
      bw.close();
      logger.info(format("\n Wrote CIF file: %s \n", saveFile.getAbsolutePath()));
    } catch (Exception ex) {
      logger.info(format("\n Failed to write out CIF file: %s \n" + ex, saveFile.getAbsolutePath()));
      return false;
    }
    if (!createdFiles.contains(saveFile.getAbsolutePath())) {
      createdFiles.add(saveFile.getAbsolutePath());
    }
    return true;
  }

  /**
   * Write an output file.
   *
   * @return Whether file was written successfully.
   */
  private boolean writeOutputFile() {
    String dir = getFullPath(files.get(0).getAbsolutePath()) + File.separator;
    String fileName = removeExtension(getName(files.get(0).getAbsolutePath()));
    String spacegroup;
    if(activeMolecularAssembly.getCrystal() != null) {
      spacegroup = activeMolecularAssembly.getCrystal().getUnitCell().spaceGroup.shortName;
    }else{
      spacegroup = null;
    }
    List<MSNode> entities = activeMolecularAssembly.getAllBondedEntities();

    File saveFile;
    if(usePDB){
      // Change name for different space groups. Input cannot handle multiple space groups in same file.
      if (entities.size() > 1) {
        fileName += "_z" + entities.size();
      }
      saveFile = new File(dir + fileName + ".pdb");
    } else if (cifFile.getBlocks().size() > 1) {
      // Change name for different space groups. Input cannot handle multiple space groups in same file.
      if (entities.size() > 1) {
        fileName += "_z" + entities.size();
      }
      fileName += "_" + spacegroup.replaceAll("\\/", "");
      saveFile = new File(dir + fileName + ".arc");
    } else {
      // If only structure then create XYZ.
      saveFile = XYZFilter.version(new File(dir + fileName + ".xyz"));
    }
    ForceField ff = activeMolecularAssembly.getForceField();
    CompositeConfiguration properties = activeMolecularAssembly.getProperties();
    if(usePDB){
      PDBFilter pdbFilter = new PDBFilter(saveFile, activeMolecularAssembly, ff, properties);
      pdbFilter.writeFile(saveFile, true);
    }else {
      XYZFilter xyzFilter = new XYZFilter(saveFile, activeMolecularAssembly, ff, properties);
      xyzFilter.writeFile(saveFile, true);
    }
    logger.info("\n Saved output file:        " + saveFile.getAbsolutePath());
    writePropertiesFile(dir, fileName, spacegroup, entities);
    if (!createdFiles.contains(saveFile.getAbsolutePath())) {
      createdFiles.add(saveFile.getAbsolutePath());
    }
    return true;
  }

  /**
   * Write an output file.
   *
   * @return Whether file was written successfully.
   */
  public boolean writeOutputFile(MolecularAssembly ma) {
    setMolecularSystem(ma);
    return writeOutputFile();
  }

  /**
   * Write a properties file.
   * @param dir String for directory location.
   * @param fileName String for file name.
   * @param spacegroup String for space group.
   * @param entities List of MSNodes in assembly.
   */
  public void writePropertiesFile(String dir, String fileName, String spacegroup, List<MSNode> entities) {
    // Write out a property file.
    File propertyFile = new File(dir + fileName + ".properties");
    if (!propertyFile.exists()) {
      try {
        FileWriter fw = new FileWriter(propertyFile, false);
        BufferedWriter bw = new BufferedWriter(fw);
        ForceField forceField = activeMolecularAssembly.getForceField();
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
        if(spacegroup != null) {
          bw.write(format("spacegroup %s\n", spacegroup));
        }
        if (entities.size() > 1) {
          bw.write("intermolecular-softcore true\n");
        }
        logger.info("\n Saved properties file: " + propertyFile.getAbsolutePath() + "\n");
        bw.close();
      } catch (Exception ex) {
        logger.info("Failed to write files.\n" + ex);
      }
    } else {
      logger.info("\n Property file already exists:  " + propertyFile.getAbsolutePath() + "\n");
    }
  }
}
