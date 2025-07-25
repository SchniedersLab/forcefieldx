// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import static ffx.potential.bonded.BondedUtils.numberAtoms;
import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard;
import static ffx.potential.bonded.PolymerUtils.assignAtomTypes;
import static ffx.potential.bonded.PolymerUtils.buildDisulfideBonds;
import static ffx.potential.bonded.PolymerUtils.buildMissingResidues;
import static ffx.potential.bonded.PolymerUtils.locateDisulfideBonds;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_2;
import static ffx.potential.parsers.PDBFilter.PDBFileStandard.VERSION3_3;
import static ffx.utilities.StringUtils.padLeft;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;
import static org.apache.commons.lang3.StringUtils.repeat;
import static org.apache.commons.math3.util.FastMath.min;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SpaceGroupDefinitions;
import ffx.crystal.SpaceGroupInfo;
import ffx.crystal.SymOp;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities.FileType;
import ffx.potential.bonded.*;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.parameters.ForceField;
import ffx.utilities.Hybrid36;
import ffx.utilities.StringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * The PDBFilter class parses data from a Protein DataBank (*.PDB) file. The following records are
 * recognized: ANISOU, ATOM, CONECT, CRYST1, END, HELIX, HETATM, LINK, SHEET, SSBOND, REMARK. The
 * rest are currently ignored.
 *
 * @author Michael J. Schnieders
 * @see <a href="http://www.wwpdb.org/documentation/format32/v3.2.html">PDB format 3.2</a>
 * @since 1.0
 */
public final class PDBFilter extends SystemFilter {

  private static final Logger logger = Logger.getLogger(PDBFilter.class.getName());
  private static final Set<String> backboneNames;
  private static final Set<String> constantPhBackboneNames;
  private static final Set<String> naBackboneNames;

  static {
    String[] names = {"C", "CA", "N", "O", "OXT", "OT2"};
    backboneNames = Set.of(names);

    String[] constantPhNames = {"C", "CA", "N", "O", "OXT", "OT2", "H", "HA", "H1", "H2", "H3"};
    constantPhBackboneNames = Set.of(constantPhNames);

    String[] naNames = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"};
    naBackboneNames = Set.of(naNames);
  }

  /** Map of SEQRES entries. */
  private final Map<Character, String[]> seqRes = new HashMap<>();
  /** Map of DBREF entries. */
  private final Map<Character, int[]> dbRef = new HashMap<>();
  /** List of altLoc characters seen in the PDB file. */
  private final List<Character> altLocs = new ArrayList<>();
  /**
   * List of segIDs defined for the PDB file.
   *
   * <p>The expectation is for chain naming from A-Z, then from 0-9. For large systems, chain names
   * are sometimes reused due to limitations in the PDB format.
   *
   * <p>However, we define segIDs to always be unique. For the first A-Z,0-9 series chainID ==
   * segID. Then, for second A-Z,0-9 series, the segID = 1A-1Z,10-19, and for the third series segID
   * = 2A-2Z,20-29, and so on.
   */
  private final List<String> segIDs = new ArrayList<>();

  private final Map<Character, List<String>> segidMap = new HashMap<>();
  /** Maps a chain to the number of insertion codes encountered in that chain. */
  private final Map<Character, Integer> insertionCodeCount = new HashMap<>();
  /**
   * Maps chainIDResNumInsCode to renumbered chainIDResNum. For example, residue 52A in chain C might
   * be renumbered to residue 53, and mapped as "C52A" to "C53".
   */
  private final Map<String, String> pdbToNewResMap = new HashMap<>();
  /** List of modified residues * */
  private final Map<String, String> modRes = new HashMap<>();
  /** Keep track of ATOM record serial numbers to match them with ANISOU records. */
  private final HashMap<Integer, Atom> atoms = new HashMap<>();

  private final Map<MolecularAssembly, BufferedReader> readers = new HashMap<>();
  /** The current altLoc - i.e., the one we are defining a chemical system for. */
  private Character currentAltLoc = 'A';
  /** Character for the current chain ID. */
  private Character currentChainID = null;
  /** String for the current SegID. */
  private String currentSegID = null;
  /** Flag to indicate a mutation is requested. */
  private boolean mutate = false;
  private List<Mutation> mutations = null;
  private List<Integer> resNumberList = null;
  /** Flag to indicate if missing fields should be printed (i.e. missing B-factors). */
  private boolean printMissingFields = true;
  /** Number of symmetry operators when expanding to a P1 unit cell (-1 saves as current spacegroup). */
  private int nSymOp = -1;
  /** Number of replicates in A lattice direction (-1 defaults to unit cell). */
  private int lValue = -1;
  /** Number of replicates in B lattice direction (-1 defaults to unit cell). */
  private int mValue = -1;
  /** Number of replicates in C lattice direction (-1 defaults to unit cell). */
  private int nValue = -1;
  /**
   * The serial field continues from the previous asymmetric unit when expanding to P1. This is not
   * used when saving as the current spacegroup.
   */
  private int serialP1 = 0;
  /** Assume current standard. */
  private PDBFileStandard fileStandard = VERSION3_3;
  /** If false, skip logging "Saving file". */
  private boolean logWrites = true;
  /** Keep track of the current MODEL in the file. */
  private int modelsRead = 1;
  /** Tracks output MODEL numbers. Unused if below zero. */
  private int modelsWritten = -1;
  /** Replicates vector dimensions if saving as expanded. */
  private int[] lmn = new int[]{1,1,1};
  private String versionFileName;

  private final File readFile;
  private List<String> remarkLines = Collections.emptyList();
  private double lastReadLambda = Double.NaN;

  /**
   * If true, read in titratable residues in their fully protonated form.
   */
  private boolean constantPH = false;
  /**
   * If true, read in titratable residues in their protonated form.
   */
  private boolean rotamerTitration = false;
  /**
   * List of residue to rename for constantPH simulations.
   */
  private static final HashMap<AminoAcid3, AminoAcid3> constantPHResidueMap = new HashMap<>();

  static {
    // Lysine
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.LYD, AminoAcidUtils.AminoAcid3.LYS);
    // Cysteine
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.CYD, AminoAcidUtils.AminoAcid3.CYS);
    // Histidine
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.HID, AminoAcidUtils.AminoAcid3.HIS);
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.HIE, AminoAcidUtils.AminoAcid3.HIS);
    // Aspartate
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.ASP, AminoAcidUtils.AminoAcid3.ASD);
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.ASH, AminoAcidUtils.AminoAcid3.ASD);
    // Glutamate
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.GLU, AminoAcidUtils.AminoAcid3.GLD);
    constantPHResidueMap.put(AminoAcidUtils.AminoAcid3.GLH, AminoAcidUtils.AminoAcid3.GLD);
  }

  /**
   * List of residue to rename for rotamer titration simulations.
   */
  private static final HashMap<AminoAcid3, AminoAcid3> rotamerResidueMap = new HashMap<>();

  static {
    // Lysine
    rotamerResidueMap.put(AminoAcidUtils.AminoAcid3.LYD, AminoAcidUtils.AminoAcid3.LYS);
    // Histidine
    rotamerResidueMap.put(AminoAcidUtils.AminoAcid3.HID, AminoAcidUtils.AminoAcid3.HIS);
    rotamerResidueMap.put(AminoAcidUtils.AminoAcid3.HIE, AminoAcidUtils.AminoAcid3.HIS);
    // Aspartate
    rotamerResidueMap.put(AminoAcidUtils.AminoAcid3.ASP, AminoAcidUtils.AminoAcid3.ASH);
    // Glutamate
    rotamerResidueMap.put(AminoAcidUtils.AminoAcid3.GLU, AminoAcidUtils.AminoAcid3.GLH);
    //Cysteine
    rotamerResidueMap.put(AminoAcidUtils.AminoAcid3.CYD, AminoAcidUtils.AminoAcid3.CYS);
  }

  /**
   * Constructor for PDBFilter.
   *
   * @param files a {@link java.util.List} object.
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   */
  public PDBFilter(List<File> files, MolecularAssembly molecularAssembly, ForceField forceField,
                   CompositeConfiguration properties) {
    super(files, molecularAssembly, forceField, properties);
    bondList = new ArrayList<>();
    this.fileType = FileType.PDB;
    readFile = files.get(0);
  }

  /**
   * Parse the PDB File from a URL.
   *
   * @param file a {@link java.io.File} object.
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   */
  public PDBFilter(File file, MolecularAssembly molecularAssembly, ForceField forceField,
                   CompositeConfiguration properties) {
    super(file, molecularAssembly, forceField, properties);
    bondList = new ArrayList<>();
    this.fileType = FileType.PDB;
    readFile = file;
  }

  /**
   * Parse the PDB File from a URL.
   *
   * @param file a {@link java.io.File} object.
   * @param molecularAssemblies a {@link java.util.List} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   */
  public PDBFilter(File file, List<MolecularAssembly> molecularAssemblies, ForceField forceField,
                   CompositeConfiguration properties) {
    super(file, molecularAssemblies, forceField, properties);
    bondList = new ArrayList<>();
    this.fileType = FileType.PDB;
    readFile = file;
  }

  /**
   * Constructor for PDBFilter with residue numbers.
   *
   * @param file a {@link java.util.List} object.
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param forceField a {@link ffx.potential.parameters.ForceField} object.
   * @param properties a {@link org.apache.commons.configuration2.CompositeConfiguration}
   *     object.
   * @param resNumberList a List of integer residue numbers for constant pH rotamer
   *     optimization.
   */
  public PDBFilter(File file, MolecularAssembly molecularAssembly, ForceField forceField,
                   CompositeConfiguration properties, List<Integer> resNumberList) {
    super(file, molecularAssembly, forceField, properties);
    bondList = new ArrayList<>();
    this.fileType = FileType.PDB;
    this.readFile = file;
    this.resNumberList = resNumberList;
    //this.chainList = chainList;
  }


  /**
   * Simple method useful for converting files to PDB format.
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @return Returns a PDB ATOM or HETATM record for the passed Atom.
   */
  public static String toPDBAtomLine(Atom atom) {
    StringBuilder sb;
    if (atom.isHetero()) {
      sb = new StringBuilder("HETATM");
    } else {
      sb = new StringBuilder("ATOM  ");
    }
    sb.append(repeat(" ", 74));

    String name = atom.getName();
    int nameLength = name.length();
    if (nameLength > 4) {
      name = name.substring(0, 4);
    } else if (nameLength == 1) {
      name = name + "  ";
    } else if (nameLength == 2) {
      name = name + " ";
    }
    int serial = atom.getXyzIndex();
    sb.replace(6, 16, format("%5s " + padLeft(name.toUpperCase(), 4), Hybrid36.encode(5, serial)));

    Character altLoc = atom.getAltLoc();
    if (altLoc != null) {
      sb.setCharAt(16, altLoc);
    } else {
      char blankChar = ' ';
      sb.setCharAt(16, blankChar);
    }

    String resName = atom.getResidueName();
    sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));

    char chain = atom.getChainID();
    sb.setCharAt(21, chain);

    int resID = atom.getResidueNumber();
    sb.replace(22, 26, format("%4s", Hybrid36.encode(4, resID)));

    double[] xyz = atom.getXYZ(null);
    StringBuilder decimals = new StringBuilder();
    for (int i = 0; i < 3; i++) {
      try {
        decimals.append(StringUtils.fwFpDec(xyz[i], 8, 3));
      } catch (IllegalArgumentException ex) {
        String newValue = StringUtils.fwFpTrunc(xyz[i], 8, 3);
        logger.info(
                format(" XYZ coordinate %8.3f for atom %s overflowed PDB format and is truncated to %s.",
                        xyz[i], atom, newValue));
        decimals.append(newValue);
      }
    }
    try {
      decimals.append(StringUtils.fwFpDec(atom.getOccupancy(), 6, 2));
    } catch (IllegalArgumentException ex) {
      logger.severe(
              format(" Occupancy %6.2f for atom %s must be between 0 and 1.", atom.getOccupancy(),
                      atom));
    }
    try {
      decimals.append(StringUtils.fwFpDec(atom.getTempFactor(), 6, 2));
    } catch (IllegalArgumentException ex) {
      String newValue = StringUtils.fwFpTrunc(atom.getTempFactor(), 6, 2);
      logger.info(
              format(" B-factor %6.2f for atom %s overflowed the PDB format and is truncated to %s.",
                      atom.getTempFactor(), atom, newValue));
      decimals.append(newValue);
    }

    sb.replace(30, 66, decimals.toString());
    sb.replace(78, 80, format("%2d", 0));
    sb.append("\n");
    return sb.toString();
  }

  public void setConstantPH(boolean constantPH) {
    this.constantPH = constantPH;
  }

  public void setRotamerTitration(boolean rotamerTitration) {
    this.rotamerTitration = rotamerTitration;
  }

  /** clearSegIDs */
  public void clearSegIDs() {
    segIDs.clear();
  }

  /** {@inheritDoc} */
  @Override
  public void closeReader() {
    for (MolecularAssembly system : systems) {
      BufferedReader br = readers.get(system);
      if (br != null) {
        try {
          br.close();
        } catch (IOException ex) {
          logger.warning(format(" Exception in closing system %s: %s", system.toString(), ex));
        }
      }
    }
  }

  @Override
  public int countNumModels() {
    Set<File> files = systems.stream().map(MolecularAssembly::getFile).map(File::toString).distinct()
            .map(File::new).collect(Collectors.toSet());

    // Dangers of parallelism are minimized by: unique files/filenames, read-only access.
    return files.parallelStream().mapToInt((File fi) -> {
      int nModelsLocal = 0;
      try (BufferedReader br = new BufferedReader(new FileReader(fi))) {
        String line = br.readLine();
        while (line != null) {
          if (line.startsWith("MODEL")) {
            ++nModelsLocal;
          }
          line = br.readLine();
        }
        nModelsLocal = Math.max(1, nModelsLocal);
      } catch (IOException ex) {
        logger.info(format(" Exception in parsing file %s: %s", fi, ex));
      }
      return nModelsLocal;
    }).sum();
  }

  /**
   * Get the list of alternate locations encountered.
   *
   * @return the alternate location list.
   */
  public List<Character> getAltLocs() {
    return altLocs;
  }

  /** {@inheritDoc} */
  @Override
  public OptionalDouble getLastReadLambda() {
    return Double.isNaN(lastReadLambda) ? OptionalDouble.empty() : OptionalDouble.of(lastReadLambda);
  }

  /**
   * Returns all the remark lines found by the last readFile call.
   *
   * @return Remark lines from the last readFile call.
   */
  @Override
  public String[] getRemarkLines() {
    int nRemarks = remarkLines.size();
    return remarkLines.toArray(new String[nRemarks]);
  }

  @Override
  public int getSnapshot() {
    return modelsRead;
  }

  /**
   * Mutate residue(s) as the PDB file is being parsed.
   *
   * @param mutations a {@link java.util.List} object.
   */
  public void mutate(List<Mutation> mutations) {
    mutate = true;
    if (this.mutations == null) {
      this.mutations = new ArrayList<>();
    }
    this.mutations.addAll(mutations);
  }

  /** Parse the PDB File */
  @Override
  public boolean readFile() {
    remarkLines = new ArrayList<>();
    // First atom is #1, to match xyz file format
    int xyzIndex = 1;
    setFileRead(false);
    systems.add(activeMolecularAssembly);

    List<String> conects = new ArrayList<>();
    List<String> links = new ArrayList<>();
    List<String> ssbonds = new ArrayList<>();
    List<String> structs = new ArrayList<>();
    try {
      for (File file : files) {
        currentFile = file;
        if (mutate) {
          List<Character> chainIDs = new ArrayList<>();
          try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line = br.readLine();
            while (line != null) {
              // Replace all tabs w/ 4x spaces
              line = line.replaceAll("\t", "    ");
              String identity = line;
              if (line.length() > 6) {
                identity = line.substring(0, 6);
              }
              identity = identity.trim().toUpperCase();
              Record record;
              try {
                record = Record.valueOf(identity);
              } catch (Exception e) {
                // Continue until the record is recognized.
                line = br.readLine();
                continue;
              }
              switch (record) {
                case ANISOU, HETATM, ATOM -> {
                  char c22 = line.charAt(21);
                  boolean idFound = false;
                  for (Character chainID : chainIDs) {
                    if (c22 == chainID) {
                      idFound = true;
                      break;
                    }
                  }
                  if (!idFound) {
                    chainIDs.add(c22);
                  }
                }
              }
              line = br.readLine();
            }
            for (Mutation mtn : mutations) {
              if (!chainIDs.contains(mtn.chainChar)) {
                if (chainIDs.size() == 1) {
                  logger.warning(
                          format(" Chain ID %c for mutation not found: only one chain %c found.",
                                  mtn.chainChar, chainIDs.get(0)));
                } else {
                  logger.warning(
                          format(" Chain ID %c for mutation not found: mutation will not proceed.",
                                  mtn.chainChar));
                }
              }
            }
          } catch (IOException ioException) {
            logger.fine(format(" Exception %s in parsing file to find chain IDs", ioException));
          }
        }

        // Check that the current file exists and that we can read it.
        if (currentFile == null || !currentFile.exists() || !currentFile.canRead()) {
          return false;
        }

        // Open the current file for parsing.
        try (BufferedReader br = new BufferedReader(new FileReader(currentFile))) {
          // Echo the alternate location being parsed.
          if (currentAltLoc == 'A') {
            logger.info(format(" Reading %s", currentFile.getName()));
          } else {
            logger.info(format(" Reading %s alternate location %s", currentFile.getName(), currentAltLoc));

          }
          activeMolecularAssembly.setAlternateLocation(currentAltLoc);

          // Reset the current chain and segID.
          currentChainID = null;
          currentSegID = null;
          boolean containsInsCode = false;

          // Read the first line of the file.
          String line = br.readLine();

          // Parse until END or ENDMDL is found, or to the end of the file.
          while (line != null) {
            // Replace all tabs w/ 4x spaces
            line = line.replaceAll("\t", "    ");
            String identity = line;
            if (line.length() > 6) {
              identity = line.substring(0, 6);
            }
            identity = identity.trim().toUpperCase();
            Record record;
            try {
              record = Record.valueOf(identity);
            } catch (Exception e) {
              // Continue until the record is recognized.
              line = br.readLine();
              continue;
            }

            // Switch on the known record.
            switch (record) {
              case ENDMDL:
              case END:
                // Setting "line" to null will exit the loop.
                line = null;
                continue;
              case DBREF:
// =============================================================================
//  1 -  6       Record name   "DBREF "
//  8 - 11       IDcode        idCode             ID code of this entry.
// 13            Character     chainID            Chain identifier.
// 15 - 18       Integer       seqBegin           Initial sequence number of the
//                                                PDB sequence segment.
// 19            AChar         insertBegin        Initial insertion code of the
//                                                PDB sequence segment.
// 21 - 24       Integer       seqEnd             Ending sequence number of the
//                                                PDB sequence segment.
// 25            AChar         insertEnd          Ending insertion code of the
//                                                PDB sequence segment.
// 27 - 32       LString       database           Sequence database name.
// 34 - 41       LString       dbAccession        Sequence database accession code.
// 43 - 54       LString       dbIdCode           Sequence  database identification code.
// 56 - 60       Integer       dbseqBegin         Initial sequence number of the
//                                                database seqment.
// 61            AChar         idbnsBeg           Insertion code of initial residue of the
//                                                segment, if PDB is the reference.
// 63 - 67       Integer       dbseqEnd           Ending sequence number of the
//                                                database segment.
// 68            AChar         dbinsEnd           Insertion code of the ending residue of
//                                                the segment, if PDB is the reference.
// =============================================================================
                Character chainID = line.substring(12, 13).toUpperCase().charAt(0);
                int seqBegin = parseInt(line.substring(14, 18).trim());
                int seqEnd = parseInt(line.substring(20, 24).trim());
                int[] seqRange = dbRef.computeIfAbsent(chainID, k -> new int[2]);
                seqRange[0] = seqBegin;
                seqRange[1] = seqEnd;
                break;
              case SEQRES:
// =============================================================================
//  1 -  6        Record name    "SEQRES"
//  8 - 10        Integer        serNum       Serial number of the SEQRES record for the
//                                            current  chain. Starts at 1 and increments
//                                            by one  each line. Reset to 1 for each chain.
// 12             Character      chainID      Chain identifier. This may be any single
//                                            legal  character, including a blank which is
//                                            is used if there is only one chain.
// 14 - 17        Integer        numRes       Number of residues in the chain.
//                                            This  value is repeated on every record.
// 20 - 22        Residue name   resName      Residue name.
// 24 - 26        Residue name   resName      Residue name.
// 28 - 30        Residue name   resName      Residue name.
// 32 - 34        Residue name   resName      Residue name.
// 36 - 38        Residue name   resName      Residue name.
// 40 - 42        Residue name   resName      Residue name.
// 44 - 46        Residue name   resName      Residue name.
// 48 - 50        Residue name   resName      Residue name.
// 52 - 54        Residue name   resName      Residue name.
// 56 - 58        Residue name   resName      Residue name.
// 60 - 62        Residue name   resName      Residue name.
// 64 - 66        Residue name   resName      Residue name.
// 68 - 70        Residue name   resName      Residue name.
// =============================================================================
                activeMolecularAssembly.addHeaderLine(line);
                chainID = line.substring(11, 12).toUpperCase().charAt(0);
                int serNum = parseInt(line.substring(7, 10).trim());
                String[] chain = seqRes.get(chainID);
                int numRes = parseInt(line.substring(13, 17).trim());
                if (chain == null) {
                  chain = new String[numRes];
                  seqRes.put(chainID, chain);
                }
                int resID = (serNum - 1) * 13;
                int end = line.length();
                for (int start = 19; start + 3 <= end; start += 4) {
                  String res = line.substring(start, start + 3).trim();
                  if (res.isEmpty()) {
                    break;
                  }
                  chain[resID++] = res;
                }
                break;
              case MODRES:
                String modResName = line.substring(12, 15).trim();
                String stdName = line.substring(24, 27).trim();
                modRes.put(modResName.toUpperCase(), stdName.toUpperCase());
                activeMolecularAssembly.addHeaderLine(line);
// =============================================================================
//  1 -  6        Record name     "MODRES"
//  8 - 11        IDcode          idCode         ID code of this entry.
// 13 - 15        Residue name    resName        Residue name used in this entry.
// 17             Character       chainID        Chain identifier.
// 19 - 22        Integer         seqNum         Sequence number.
// 23             AChar           iCode          Insertion code.
// 25 - 27        Residue name    stdRes         Standard residue name.
// 30 - 70        String          comment        Description of the residue modification.
// =============================================================================
                break;
              case ANISOU:
// =============================================================================
//  1 - 6        Record name   "ANISOU"
//  7 - 11       Integer       serial         Atom serial number.
// 13 - 16       Atom          name           Atom name.
// 17            Character     altLoc         Alternate location indicator
// 18 - 20       Residue name  resName        Residue name.
// 22            Character     chainID        Chain identifier.
// 23 - 26       Integer       resSeq         Residue sequence number.
// 27            AChar         iCode          Insertion code.
// 29 - 35       Integer       u[0][0]        U(1,1)
// 36 - 42       Integer       u[1][1]        U(2,2)
// 43 - 49       Integer       u[2][2]        U(3,3)
// 50 - 56       Integer       u[0][1]        U(1,2)
// 57 - 63       Integer       u[0][2]        U(1,3)
// 64 - 70       Integer       u[1][2]        U(2,3)
// 77 - 78       LString(2)    element        Element symbol, right-justified.
// 79 - 80       LString(2)    charge         Charge on the atom.
// =============================================================================
                boolean deleteAnisou = properties.getBoolean("delete-anisou", false);
                double resetBfactors = properties.getDouble("reset-bfactors", -1.0);
                if (deleteAnisou || resetBfactors >= 0.0) {
                  break;
                }
                Integer serial = Hybrid36.decode(5, line.substring(6, 11));
                Character altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                if (!altLocs.contains(altLoc)) {
                  altLocs.add(altLoc);
                }
                if (!altLoc.equals(' ') && !altLoc.equals('A') && !altLoc.equals(currentAltLoc)) {
                  break;
                }
                double[] adp = new double[6];
                adp[0] = parseInt(line.substring(28, 35).trim()) * 1.0e-4;
                adp[1] = parseInt(line.substring(35, 42).trim()) * 1.0e-4;
                adp[2] = parseInt(line.substring(42, 49).trim()) * 1.0e-4;
                adp[3] = parseInt(line.substring(49, 56).trim()) * 1.0e-4;
                adp[4] = parseInt(line.substring(56, 63).trim()) * 1.0e-4;
                adp[5] = parseInt(line.substring(63, 70).trim()) * 1.0e-4;
                if (atoms.containsKey(serial)) {
                  Atom a = atoms.get(serial);
                  a.setAltLoc(altLoc);
                  a.setAnisou(adp);
                } else {
                  logger.info(
                          format(" No ATOM record for ANISOU serial number %d has been found.", serial));
                  logger.info(format(" This ANISOU record will be ignored:\n %s", line));
                }
                break;
              case ATOM:
// =============================================================================
//  1 -  6        Record name   "ATOM  "
//  7 - 11        Integer       serial       Atom serial number.
// 13 - 16        Atom          name         Atom name.
// 17             Character     altLoc       Alternate location indicator.
// 18 - 20        Residue name  resName      Residue name.
// 22             Character     chainID      Chain identifier.
// 23 - 26        Integer       resSeq       Residue sequence number.
// 27             AChar         iCode        Code for insertion of residues.
// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
// 55 - 60        Real(6.2)     occupancy    Occupancy.
// 61 - 66        Real(6.2)     tempFactor   Temperature factor.
// 77 - 78        LString(2)    element      Element symbol, right-justified.
// 79 - 80        LString(2)    charge       Charge  on the atom.
// =============================================================================
                String name;
                String resName;
                String segID;
                int resSeq;
                boolean printAtom;
                double[] d;
                double occupancy;
                double tempFactor;
                Atom newAtom;
                Atom returnedAtom;
                // If it's a misnamed water, it will fall through to HETATM.
                if (!line.substring(17, 20).trim().equals("HOH")) {
                  serial = Hybrid36.decode(5, line.substring(6, 11));
                  name = line.substring(12, 16).trim();
                  if (name.toUpperCase().contains("1H") || name.toUpperCase().contains("2H")
                          || name.toUpperCase().contains("3H")) {
                    // VERSION3_2 is presently just a placeholder for "anything non-standard".
                    fileStandard = VERSION3_2;
                  }
                  altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                  if (!altLocs.contains(altLoc)) {
                    altLocs.add(altLoc);
                  }

                  if (!altLoc.equals(' ') && !altLoc.equals('A') && !altLoc.equals(currentAltLoc)) {
                    break;
                  }

                  resName = line.substring(17, 20).trim();
                  chainID = line.substring(21, 22).charAt(0);
                  segID = getSegID(chainID);
                  resSeq = Hybrid36.decode(4, line.substring(22, 26));

                  char insertionCode = line.charAt(26);
                  if (insertionCode != ' ' && !containsInsCode) {
                    containsInsCode = true;
                    logger.warning(
                            " FFX support for files with " + "insertion codes is experimental. "
                                    + "Residues will be renumbered to " + "eliminate insertion codes (52A "
                                    + "becomes 53, 53 becomes 54, etc)");
                  }

                  int offset = insertionCodeCount.getOrDefault(chainID, 0);
                  String pdbResNum = format("%c%d%c", chainID, resSeq, insertionCode);
                  if (!pdbToNewResMap.containsKey(pdbResNum)) {
                    if (insertionCode != ' ') {
                      ++offset;
                      insertionCodeCount.put(chainID, offset);
                    }
                    resSeq += offset;
                    if (offset != 0) {
                      logger.info(
                              format(" Chain %c " + "residue %s-%s renumbered to %c %s-%d", chainID,
                                      pdbResNum.substring(1).trim(), resName, chainID, resName, resSeq));
                    }
                    String newNum = format("%c%d", chainID, resSeq);
                    pdbToNewResMap.put(pdbResNum, newNum);
                  } else {
                    resSeq += offset;
                  }

                  printAtom = false;
                  if (mutate) {
                    boolean doBreak = false;
                    for (Mutation mtn : mutations) {
                      if (chainID == mtn.chainChar && resSeq == mtn.resID) {
                        mtn.origResName = resName;
                        resName = mtn.resName;
                        String atomName = name.toUpperCase();

                        int isAA = AminoAcidUtils.getAminoAcidNumber(resName);
                        int isNA = NucleicAcidUtils.getNucleicAcidNumber(resName);

                        if ((isNA != -1 && naBackboneNames.contains(atomName)) || (isAA != -1 && backboneNames.contains(atomName))) {
                          printAtom = true;
                        } else {
                          // grab pur-pur or pyr-pyr alchem. atoms
                          ArrayList<String> alchAtoms = mtn.getAlchemicalAtoms(false);
                          if (alchAtoms == null) {
                            // test to see if atom is involved in glycosyl torsion and if needs renaming (pyr-pur/pur-pyr)
                            String newName = mtn.isNonAlchemicalAtom(atomName);
                            if (newName != null) { // if not null -- use name it
                              printAtom = true;
                              if (newName.startsWith("~")) { // switch from purine to pyrmidine or v.v.
                                // switch name and include it as an alchemical atom
                                name = newName.substring(1);
                                logger.info(format(" DELETING atom %d %s of %s %d in chain %s", serial, atomName, resName, resSeq, chainID));
                              } else {
                                // replace name but do not include as an alchemical atom
                                name = newName;
                              }
                              doBreak = false;
                            } else if (!atomName.contains("'")) {
                              logger.info(format(" DELETING atom %d %s of %s %d in chain %s", serial, atomName, resName, resSeq, chainID));
                              doBreak = true;
                            } else {
                              printAtom = true;
                              doBreak = false;
                            }
                          } else {
                            if (alchAtoms.contains(atomName) && !atomName.contains("'")) {
                              logger.info(format(" DELETING atom %d %s of %s %d in chain %s", serial, atomName, resName, resSeq, chainID));
                              doBreak = true;
                            } else {
                              printAtom = true;
                              doBreak = false;
                            }
                          }
                          break;
                        }
                      }
                    }
                    if (doBreak) {
                      break;
                    }
                  }

                  if (constantPH) {
                    AminoAcid3 aa3 = AminoAcidUtils.getAminoAcid(resName.toUpperCase());
                    if (constantPHResidueMap.containsKey(aa3)) {
                      String atomName = name.toUpperCase();
                      AminoAcid3 aa3PH = constantPHResidueMap.get(aa3);
                      resName = aa3PH.name();
                      if (constantPhBackboneNames.contains(atomName)) {
                        logger.info(format(" %s-%d %s", resName, resSeq, atomName));
                      } else if (!atomName.startsWith("H")) {
                        logger.info(format(" %s-%d %s", resName, resSeq, atomName));
                      } else {
                        logger.info(format(" %s-%d %s skipped", resName, resSeq, atomName));
                        break;
                      }
                    }
                  } else if (rotamerTitration) {
                    AminoAcid3 aa3 = AminoAcidUtils.getAminoAcid(resName.toUpperCase());
                    if (rotamerResidueMap.containsKey(aa3) && resNumberList.contains(resSeq)) {
                      AminoAcid3 aa3rotamer = rotamerResidueMap.get(aa3);
                      resName = aa3rotamer.name();
                    }
                  }
                  d = new double[3];
                  d[0] = parseDouble(line.substring(30, 38).trim());
                  d[1] = parseDouble(line.substring(38, 46).trim());
                  d[2] = parseDouble(line.substring(46, 54).trim());
                  occupancy = 1.0;
                  tempFactor = 1.0;
                  try {
                    occupancy = parseDouble(line.substring(54, 60).trim());
                    tempFactor = parseDouble(line.substring(60, 66).trim());
                  } catch (NumberFormatException | StringIndexOutOfBoundsException e) {
                    // Use default values.
                    if (printMissingFields) {
                      logger.info(" Missing occupancy and b-factors set to 1.0.");
                      printMissingFields = false;
                    } else if (logger.isLoggable(Level.FINE)) {
                      logger.fine(" Missing occupancy and b-factors set to 1.0.");
                    }
                  }

                  double bfactor = properties.getDouble("reset-bfactors", -1.0);
                  if (bfactor >= 0.0) {
                    tempFactor = bfactor;
                  }

                  newAtom = new Atom(0, name, altLoc, d, resName, resSeq, chainID, occupancy, tempFactor, segID);

                  // Check if this is a modified residue.
                  if (modRes.containsKey(resName.toUpperCase())) {
                    newAtom.setModRes(true);
                  }
                  returnedAtom = (Atom) activeMolecularAssembly.addMSNode(newAtom);
                  if (returnedAtom != newAtom) {
                    // A previously added atom has been retained.
                    atoms.put(serial, returnedAtom);
                    if (logger.isLoggable(Level.FINE)) {
                      logger.fine(returnedAtom + " has been retained over\n" + newAtom);
                    }
                  } else {
                    // The new atom has been added.
                    atoms.put(serial, newAtom);
                    // Check if the newAtom took the xyzIndex of a previous alternate conformer.
                    if (newAtom.getIndex() == 0) {
                      newAtom.setXyzIndex(xyzIndex++);
                    }
                    if (printAtom) {
                      logger.info(newAtom.toString());
                    }
                  }
                  break;
                }
                break;
              case HETATM:
// =============================================================================
//  1 - 6        Record name    "HETATM"
//  7 - 11       Integer        serial        Atom serial number.
// 13 - 16       Atom           name          Atom name.
// 17            Character      altLoc        Alternate location indicator.
// 18 - 20       Residue name   resName       Residue name.
// 22            Character      chainID       Chain identifier.
// 23 - 26       Integer        resSeq        Residue sequence number.
// 27            AChar          iCode         Code for insertion of residues.
// 31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
// 39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
// 47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
// 55 - 60       Real(6.2)      occupancy     Occupancy.
// 61 - 66       Real(6.2)      tempFactor    Temperature factor.
// 77 - 78       LString(2)     element       Element symbol; right-justified.
// 79 - 80       LString(2)     charge        Charge on the atom.
// =============================================================================
                serial = Hybrid36.decode(5, line.substring(6, 11));
                name = line.substring(12, 16).trim();
                altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                if (!altLocs.contains(altLoc)) {
                  altLocs.add(altLoc);
                }
                if (!altLoc.equals(' ') && !altLoc.equals(currentAltLoc)) {
                  break;
                }
                // if (!altLoc.equals(' ') && !altLoc.equals('A') && !altLoc.equals(currentAltLoc)) {
                //  break;
                //}
                resName = line.substring(17, 20).trim();
                chainID = line.substring(21, 22).charAt(0);
                segID = getSegID(chainID);
                resSeq = Hybrid36.decode(4, line.substring(22, 26));

                char insertionCode = line.charAt(26);
                if (insertionCode != ' ' && !containsInsCode) {
                  containsInsCode = true;
                  logger.warning(" FFX support for files with " + "insertion codes is experimental. "
                          + "Residues will be renumbered to " + "eliminate insertion codes (52A "
                          + "becomes 53, 53 becomes 54, etc)");
                }

                int offset = insertionCodeCount.getOrDefault(chainID, 0);
                String pdbResNum = format("%c%d%c", chainID, resSeq, insertionCode);
                if (!pdbToNewResMap.containsKey(pdbResNum)) {
                  if (insertionCode != ' ') {
                    ++offset;
                    insertionCodeCount.put(chainID, offset);
                  }
                  resSeq += offset;
                  if (offset != 0) {
                    logger.info(
                            format(" Chain %c " + "molecule %s-%s renumbered to %c %s-%d", chainID,
                                    pdbResNum.substring(1).trim(), resName, chainID, resName, resSeq));
                  }
                  String newNum = format("%c%d", chainID, resSeq);
                  pdbToNewResMap.put(pdbResNum, newNum);
                } else {
                  resSeq += offset;
                }

                d = new double[3];
                d[0] = parseDouble(line.substring(30, 38).trim());
                d[1] = parseDouble(line.substring(38, 46).trim());
                d[2] = parseDouble(line.substring(46, 54).trim());
                occupancy = 1.0;
                tempFactor = 1.0;
                try {
                  occupancy = parseDouble(line.substring(54, 60).trim());
                  tempFactor = parseDouble(line.substring(60, 66).trim());
                } catch (NumberFormatException | StringIndexOutOfBoundsException e) {
                  // Use default values.
                  if (printMissingFields) {
                    logger.info(" Missing occupancy and b-factors set to 1.0.");
                    printMissingFields = false;
                  } else if (logger.isLoggable(Level.FINE)) {
                    logger.fine(" Missing occupancy and b-factors set to 1.0.");
                  }
                }

                double bfactor = properties.getDouble("reset-bfactors", -1.0);
                if (bfactor >= 0.0) {
                  tempFactor = bfactor;
                }

                newAtom = new Atom(0, name, altLoc, d, resName, resSeq, chainID, occupancy, tempFactor, segID);
                newAtom.setHetero(true);
                // Check if this is a modified residue.
                if (modRes.containsKey(resName.toUpperCase())) {
                  newAtom.setModRes(true);
                }
                returnedAtom = (Atom) activeMolecularAssembly.addMSNode(newAtom);
                if (returnedAtom != newAtom) {
                  // A previously added atom has been retained.
                  atoms.put(serial, returnedAtom);
                  if (logger.isLoggable(Level.FINE)) {
                    logger.fine(returnedAtom + " has been retained over\n" + newAtom);
                  }
                } else {
                  // The new atom has been added.
                  atoms.put(serial, newAtom);
                  newAtom.setXyzIndex(xyzIndex++);
                }
                break;
              case CRYST1:
// =============================================================================
// The CRYST1 record presents the unit cell parameters, space group, and Z
// value. If the structure was not determined by crystallographic means, CRYST1
// simply provides the unitary values, with an appropriate REMARK.
//
//  7 - 15       Real(9.3)     a              a (Angstroms).
// 16 - 24       Real(9.3)     b              b (Angstroms).
// 25 - 33       Real(9.3)     c              c (Angstroms).
// 34 - 40       Real(7.2)     alpha          alpha (degrees).
// 41 - 47       Real(7.2)     beta           beta (degrees).
// 48 - 54       Real(7.2)     gamma          gamma (degrees).
// 56 - 66       LString       sGroup         Space  group.
// 67 - 70       Integer       z              Z value.
// =============================================================================
                if (line.length() < 55) {
                  logger.severe(" CRYST1 record is improperly formatted.");
                }
                double aaxis = parseDouble(line.substring(6, 15).trim());
                double baxis = parseDouble(line.substring(15, 24).trim());
                double caxis = parseDouble(line.substring(24, 33).trim());
                double alpha = parseDouble(line.substring(33, 40).trim());
                double beta = parseDouble(line.substring(40, 47).trim());
                double gamma = parseDouble(line.substring(47, 54).trim());
                int limit = min(line.length(), 66);
                String sg = line.substring(55, limit).trim();
                properties.addProperty("a-axis", aaxis);
                properties.addProperty("b-axis", baxis);
                properties.addProperty("c-axis", caxis);
                properties.addProperty("alpha", alpha);
                properties.addProperty("beta", beta);
                properties.addProperty("gamma", gamma);
                properties.addProperty("spacegroup", SpaceGroupInfo.pdb2ShortName(sg));
                break;
              case CONECT:
// =============================================================================
//  7 - 11        Integer        serial       Atom  serial number
// 12 - 16        Integer        serial       Serial number of bonded atom
// 17 - 21        Integer        serial       Serial number of bonded atom
// 22 - 26        Integer        serial       Serial number of bonded atom
// 27 - 31        Integer        serial       Serial number of bonded atom
//
// CONECT records involving atoms for which the coordinates are not present
// in the entry (e.g., symmetry-generated) are not given.
// CONECT records involving atoms for which the coordinates are missing due
// to disorder, are also not provided.
// =============================================================================
                conects.add(line);
                break;
              case LINK:
// =============================================================================
// The LINK records specify connectivity between residues that is not implied by
// the primary structure. Connectivity is expressed in terms of the atom names.
// They also include the distance associated with each linkage following the
// symmetry operations at the end of each record.
// 13 - 16         Atom           name1           Atom name.
// 17              Character      altLoc1         Alternate location indicator.
// 18 - 20         Residue name   resName1        Residue  name.
// 22              Character      chainID1        Chain identifier.
// 23 - 26         Integer        resSeq1         Residue sequence number.
// 27              AChar          iCode1          Insertion code.
// 43 - 46         Atom           name2           Atom name.
// 47              Character      altLoc2         Alternate location indicator.
// 48 - 50         Residue name   resName2        Residue name.
// 52              Character      chainID2        Chain identifier.
// 53 - 56         Integer        resSeq2         Residue sequence number.
// 57              AChar          iCode2          Insertion code.
// 60 - 65         SymOP          sym1            Symmetry operator atom 1.
// 67 - 72         SymOP          sym2            Symmetry operator atom 2.
// 74 – 78         Real(5.2)      Length          Link distance
// =============================================================================
                char a1 = line.charAt(16);
                char a2 = line.charAt(46);
                if (a1 != a2) {
                  // logger.info(format(" Ignoring LINK record as alternate locations do not match\n
                  // %s.", line));
                  break;
                }
                if (currentAltLoc == 'A') {
                  if ((a1 == ' ' || a1 == 'A') && (a2 == ' ' || a2 == 'A')) {
                    links.add(line);
                  }
                } else if (a1 == currentAltLoc && a2 == currentAltLoc) {
                  links.add(line);
                }
                break;
              case SSBOND:
// =============================================================================
// The SSBOND record identifies each disulfide bond in protein and polypeptide
// structures by identifying the two residues involved in the bond.
// The disulfide bond distance is included after the symmetry operations at
// the end of the SSBOND record.
//
//  8 - 10        Integer         serNum       Serial number.
// 12 - 14        LString(3)      "CYS"        Residue name.
// 16             Character       chainID1     Chain identifier.
// 18 - 21        Integer         seqNum1      Residue sequence number.
// 22             AChar           icode1       Insertion code.
// 26 - 28        LString(3)      "CYS"        Residue name.
// 30             Character       chainID2     Chain identifier.
// 32 - 35        Integer         seqNum2      Residue sequence number.
// 36             AChar           icode2       Insertion code.
// 60 - 65        SymOP           sym1         Symmetry oper for 1st resid
// 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
// 74 – 78        Real(5.2)      Length        Disulfide bond distance
//
// If SG of cysteine is disordered then there are possible alternate linkages.
// wwPDB practice is to put together all possible SSBOND records. This is
// problematic because the alternate location identifier is not specified in
// the SSBOND record.
//
// Notes:
// SSBOND records may be invalid if chain IDs are reused.
// SSBOND records are applied by FFX to all conformers.
// =============================================================================
                ssbonds.add(line);
                break;
              case HELIX:
// =============================================================================
// HELIX records are used to identify the position of helices in the molecule.
// Helices are named, numbered, and classified by type. The residues where the
// helix begins and ends are noted, as well as the total length.
//
//  8 - 10        Integer        serNum        Serial number of the helix. This starts
//                                             at 1  and increases incrementally.
// 12 - 14        LString(3)     helixID       Helix  identifier. In addition to a serial
//                                             number, each helix is given an
//                                             alphanumeric character helix identifier.
// 16 - 18        Residue name   initResName   Name of the initial residue.
// 20             Character      initChainID   Chain identifier for the chain containing
//                                             this  helix.
// 22 - 25        Integer        initSeqNum    Sequence number of the initial residue.
// 26             AChar          initICode     Insertion code of the initial residue.
// 28 - 30        Residue  name  endResName    Name of the terminal residue of the helix.
// 32             Character      endChainID    Chain identifier for the chain containing
//                                             this  helix.
// 34 - 37        Integer        endSeqNum     Sequence number of the terminal residue.
// 38             AChar          endICode      Insertion code of the terminal residue.
// 39 - 40        Integer        helixClass    Helix class (see below).
// 41 - 70        String         comment       Comment about this helix.
// 72 - 76        Integer        length        Length of this helix.
//
//                                      CLASS NUMBER
// TYPE OF  HELIX                     (COLUMNS 39 - 40)
// --------------------------------------------------------------
// Right-handed alpha (default)                1
// Right-handed omega                          2
// Right-handed pi                             3
// Right-handed gamma                          4
// Right-handed 3 - 10                         5
// Left-handed alpha                           6
// Left-handed omega                           7
// Left-handed gamma                           8
// 2 - 7 ribbon/helix                          9
// Polyproline                                10
// =============================================================================
              case SHEET:
// =============================================================================
// SHEET records are used to identify the position of sheets in the molecule.
// Sheets are both named and numbered. The residues where the sheet begins and
// ends are noted.
//
//  8 - 10        Integer       strand         Strand  number which starts at 1 for each
//                                             strand within a sheet and increases by one.
// 12 - 14        LString(3)    sheetID        Sheet  identifier.
// 15 - 16        Integer       numStrands     Number  of strands in sheet.
// 18 - 20        Residue name  initResName    Residue  name of initial residue.
// 22             Character     initChainID    Chain identifier of initial residue in strand.
// 23 - 26        Integer       initSeqNum     Sequence number of initial residue in strand.
// 27             AChar         initICode      Insertion code of initial residue in  strand.
// 29 - 31        Residue name  endResName     Residue name of terminal residue.
// 33             Character     endChainID     Chain identifier of terminal residue.
// 34 - 37        Integer       endSeqNum      Sequence number of terminal residue.
// 38             AChar         endICode       Insertion code of terminal residue.
// 39 - 40        Integer       sense          Sense of strand with respect to previous
//                                             strand in the sheet. 0 if first strand,
//                                             1 if  parallel,and -1 if anti-parallel.
// 42 - 45        Atom          curAtom        Registration.  Atom name in current strand.
// 46 - 48        Residue name  curResName     Registration.  Residue name in current strand
// 50             Character     curChainId     Registration. Chain identifier in current strand.
// 51 - 54        Integer       curResSeq      Registration.  Residue sequence number
//                                             in current strand.
// 55             AChar         curICode       Registration. Insertion code in current strand.
// 57 - 60        Atom          prevAtom       Registration.  Atom name in previous strand.
// 61 - 63        Residue name  prevResName    Registration.  Residue name in previous strand.
// 65             Character     prevChainId    Registration.  Chain identifier in previous strand.
// 66 - 69        Integer       prevResSeq     Registration. Residue sequence number
//                                             in previous strand.
// 70             AChar         prevICode      Registration.  Insertion code in
//                                             previous strand.
// =============================================================================
                structs.add(line);
                break;
              case MODEL: // Currently, no handling in initial read.
                break;
              case MTRIX1:
// ================================================================================
// MTRIXn (n = 1, 2, or 3) records present transformations expressing
// non-crystallographic symmetry.
// MTRIXn will appear only when such transformations are required to generate an
// entire asymmetric unit,
// such as a large viral structure.
//
//  8 - 10        Integer       serial         Serial number.
// 11 - 20        Real(10.6)    m[n][1]        Mn1
// 21 - 30        Real(10.6)    m[n][2]        Mn2
// 31 - 40        Real(10.6)    m[n][3]        Mn3
// 21 - 30        Real(10.6)    v[n]           Vn
// 60             Integer       iGiven         1 if coordinates for the representations which are
//                                              approximately related by the transformations of the
//                                              molecule are contained in the entry. Otherwise, blank.
// =================================================================================
                StringBuilder MTRX1 = new StringBuilder(line.substring(11, 55));
                properties.addProperty("MTRIX1", MTRX1);
                break;
              case MTRIX2:
                StringBuilder MTRX2 = new StringBuilder(line.substring(11, 55));
                properties.addProperty("MTRIX2", MTRX2);
                break;
              case MTRIX3:
                StringBuilder MTRX3 = new StringBuilder(line.substring(11, 55));
                properties.addProperty("MTRIX3", MTRX3);
                break;
              case REMARK:
                remarkLines.add(line.trim());
                if (line.contains("Lambda:")) {
                  Matcher m = lambdaPattern.matcher(line);
                  if (m.find()) {
                    lastReadLambda = Double.parseDouble(m.group(1));
                  }
                }
// =================================================================================
// REMARK 350: presents all transformations, both crystallographic and non-crystallographic,
// needed to generate the biomolecule. These transformations operate on the coordinates in the
// entry. Both author and computational descriptions of assemblies are provided, if applicable.
// For strict ncs case where more than one assembly presents in asymmetric unit, only one
// chain with unit matrix will reported in REMARK 350, the other chain will be generated
// by rotation and translation.
//
// 20 - 23        Integer       serial         Serial number.
// 24 - 33        Real(10.6)    m[n][1]        Mn1
// 34 - 43        Real(10.6)    m[n][2]        Mn2
// 44 - 53        Real(10.6)    m[n][3]        Mn3
// 59 - 68        Real(10.6)    v[n]           Vn
// =================================================================================
                if (line.length() >= 68) {
                  String remarkType = line.substring(7, 10).trim();
                  if (remarkType.matches("\\d+") && parseInt(remarkType) == 350 && line.substring(13,
                          18).equalsIgnoreCase("BIOMT")) {
                    properties.addProperty("BIOMTn", new StringBuilder(line.substring(24, 68)));
                  }
                }
                break;
              default:
                break;
            }
            line = br.readLine();
          }

        } catch (FileNotFoundException fileNotFoundException) {
          logger.log(Level.SEVERE, " PDB file not found", fileNotFoundException);
        }
      }
      xyzIndex--;
      setFileRead(true);
    } catch (IOException e) {
      logger.exiting(PDBFilter.class.getName(), "readFile", e);
      return false;
    }

    // Locate disulfide bonds; bond parameters are assigned below.
    List<Bond> ssBondList = locateDisulfideBonds(ssbonds, activeMolecularAssembly, pdbToNewResMap);

    // Record the number of atoms read in from the PDB file before applying
    // algorithms that may build new atoms.
    int pdbAtoms = activeMolecularAssembly.getAtomArray().length;
    removeExcessHydrogens();
    // Build missing backbone atoms in loops.
    buildMissingResidues(xyzIndex, activeMolecularAssembly, seqRes, dbRef);

    // Assign atom types. Missing side-chains atoms and missing hydrogen will be built in.
    bondList = assignAtomTypes(activeMolecularAssembly, fileStandard);

    // Assign disulfide bonds parameters and log their creation.
    buildDisulfideBonds(ssBondList, activeMolecularAssembly, bondList);

    // Finally, re-number the atoms if missing atoms were created.
    int currentN = activeMolecularAssembly.getAtomArray().length;
    boolean renumber = forceField.getBoolean("renumber-pdb", false);
    if (pdbAtoms != currentN) {
      logger.info(format(" Renumbering PDB file due to built atoms (%d vs %d)", currentN, pdbAtoms));
      numberAtoms(activeMolecularAssembly);
    } else if (renumber) {
      logger.info(" Renumbering PDB file due to renumber-pdb flag.");
      numberAtoms(activeMolecularAssembly);
    }
    return true;
  }

  public void removeExcessHydrogens(){
    logger.info(" Removing excess Hydrogens");
    for(Residue residue: activeMolecularAssembly.getResidueList()){
      if(residue.getName().equals("ACE") || residue.getName().equals("NME") || residue.getResidueType() != Residue.ResidueType.AA){
        break;
      }
      String trueResName = residue.getAtomByName("CA", true).getResidueName();
      Atom atom;
      switch (trueResName) {
        case "HID", "GLU" -> {
          // No HE2
          atom = residue.getAtomByName("HE2", true);
        }
        case "HIE" -> {
          // No HD1
          atom = residue.getAtomByName("HD1", true);
        }
        case "ASP" -> {
          // No HD2
          atom = residue.getAtomByName("HD2", true);
        }
        case "LYD" -> {
          // No HZ3
          atom = residue.getAtomByName("HZ3", true);
        }
        case "CYD" -> {
          // No HG
          atom = residue.getAtomByName("HG", true);
        }
        default -> {
          atom = null;
        }
        // Do nothing.
      }
      if(atom != null){
        int index = activeMolecularAssembly.getResidueList().indexOf(residue);
        MSNode atoms = residue.getAtomNode();
        atoms.remove(atom);
        residue.setName(trueResName);
        activeMolecularAssembly.getResidueList().set(index, residue);
      }
    }
  }

  /** {@inheritDoc} */
  @Override
  public boolean readNext() {
    return readNext(false);
  }

  /** {@inheritDoc} */
  @Override
  public boolean readNext(boolean resetPosition) {
    return readNext(resetPosition, false);
  }

  /** {@inheritDoc} */
  @Override
  public boolean readNext(boolean resetPosition, boolean print) {
    return readNext(resetPosition, print, true);
  }

  /** {@inheritDoc} */
  @Override
  @SuppressWarnings("fallthrough")
  public boolean readNext(boolean resetPosition, boolean print, boolean parse) {
    modelsRead = resetPosition ? 1 : modelsRead + 1;
    if (!parse) {
      if (print) {
        logger.info(format(" Skipped Model %d.", modelsRead));
      }
      return true;
    }
    remarkLines = new ArrayList<>(remarkLines.size());
    // ^ is beginning of line, \\s+ means "one or more whitespace", (\\d+) means match and capture
    // one or more digits.
    Pattern modelPatt = Pattern.compile("^MODEL\\s+(\\d+)");
    boolean eof = true;
    for (MolecularAssembly system : systems) {
      try {
        BufferedReader currentReader;
        if (readers.containsKey(system)) {
          currentReader = readers.get(system);
          try {
            if (!currentReader.ready()) {
              currentReader = new BufferedReader(new FileReader(readFile));
              // Mark the start of the file.
              currentReader.mark(0);
              readers.remove(system);
              readers.put(system, currentReader);
            } else if (resetPosition) {
              // If the BufferedReader has been opened, and reset is requested, reset the position.
              currentReader.reset();
            }
          } catch (Exception exception) {
            // If all structures in the PDB file have been read, the currentReader may have closed.
            // The try block will catch this case and reset to the beginning of the file.
            currentReader = new BufferedReader(new FileReader(readFile));
            // Mark the start of the file.
            currentReader.mark(0);
            readers.remove(system);
            readers.put(system, currentReader);
          }
        } else {
          currentReader = new BufferedReader(new FileReader(readFile));
          // Mark the start of the file.
          currentReader.mark(0);
          readers.put(system, currentReader);
        }

        // Skip to appropriate model.
        String line = currentReader.readLine();
        while (line != null) {
          line = line.trim();
          Matcher m = modelPatt.matcher(line);
          if (m.find()) {
            int modelNum = parseInt(m.group(1));
            if (modelNum == modelsRead) {
              if (print) {
                logger.log(Level.INFO, format(" Reading model %d for %s", modelNum, currentFile));
              }
              eof = false;
              break;
            }
          }
          line = currentReader.readLine();
        }
        if (eof) {
          if (logger.isLoggable(Level.FINEST)) {
            logger.log(Level.FINEST, format("\n End of file reached for %s", readFile));
          }
          currentReader.close();
          return false;
        }

        // Begin parsing the model.
        currentChainID = null;
        currentSegID = null;
        boolean modelDone = false;
        line = currentReader.readLine();
        while (line != null) {
          line = line.trim();
          String recID = line.substring(0, Math.min(6, line.length())).trim();
          try {
            Record record = Record.valueOf(recID);
            boolean hetatm = true;
            switch (record) {
              // =============================================================================
              //
              //  7 - 11        Integer       serial       Atom serial number.
              // 13 - 16        Atom          name         Atom name.
              // 17             Character     altLoc       Alternate location indicator.
              // 18 - 20        Residue name  resName      Residue name.
              // 22             Character     chainID      Chain identifier.
              // 23 - 26        Integer       resSeq       Residue sequence number.
              // 27             AChar         iCode        Code for insertion of residues.
              // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in
              // Angstroms.
              // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in
              // Angstroms.
              // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in
              // Angstroms.
              // 55 - 60        Real(6.2)     occupancy    Occupancy.
              // 61 - 66        Real(6.2)     tempFactor   Temperature factor.
              // 77 - 78        LString(2)    element      Element symbol, right-justified.
              // 79 - 80        LString(2)    charge       Charge  on the atom.
              // =============================================================================
              //         1         2         3         4         5         6         7
              // 123456789012345678901234567890123456789012345678901234567890123456789012345678
              // ATOM      1  N   ILE A  16      60.614  71.140 -10.592  1.00  7.38           N
              // ATOM      2  CA  ILE A  16      60.793  72.149  -9.511  1.00  6.91           C
              case ATOM:
                hetatm = false;
              case HETATM:
                String name = line.substring(12, 16).trim();
                if (name.toUpperCase().contains("1H") || name.toUpperCase().contains("2H")
                        || name.toUpperCase().contains("3H")) {
                  // VERSION3_2 is presently just a placeholder for "anything non-standard".
                  fileStandard = VERSION3_2;
                }
                Character altLoc = line.substring(16, 17).toUpperCase().charAt(0);
                if (!altLoc.equals(' ') && !altLoc.equals(currentAltLoc)) {
                  break;
                }
                // if (!altLoc.equals(' ') && !altLoc.equals('A') && !altLoc.equals(currentAltLoc)) {
                //  break;
                // }
                String resName = line.substring(17, 20).trim();
                Character chainID = line.substring(21, 22).charAt(0);
                String segID = getExistingSegID(chainID);
                int resSeq = Hybrid36.decode(4, line.substring(22, 26));
                double[] d = new double[3];
                d[0] = parseDouble(line.substring(30, 38).trim());
                d[1] = parseDouble(line.substring(38, 46).trim());
                d[2] = parseDouble(line.substring(46, 54).trim());
                double occupancy = 1.0;
                double tempFactor = 1.0;
                Atom newAtom = new Atom(0, name, altLoc, d, resName, resSeq, chainID, occupancy,
                        tempFactor, segID);
                newAtom.setHetero(hetatm);
                // Check if this is a modified residue.
                if (modRes.containsKey(resName.toUpperCase())) {
                  newAtom.setModRes(true);
                }

                Atom returnedAtom = activeMolecularAssembly.findAtom(newAtom);
                if (returnedAtom != null) {
                  returnedAtom.setXYZ(d);
                  double[] retXYZ = new double[3];
                  returnedAtom.getXYZ(retXYZ);
                } else {
                  String message = format(" Could not find atom %s in assembly", newAtom);
                  if (dieOnMissingAtom) {
                    logger.severe(message);
                  } else {
                    logger.warning(message);
                  }
                }
                break;
              case CRYST1:
                // =============================================================================
                // The CRYST1 record presents the unit cell parameters, space group, and Z
                // value. If the structure was not determined by crystallographic means, CRYST1
                // simply provides the unitary values, with an appropriate REMARK.
                //
                //  7 - 15       Real(9.3)     a              a (Angstroms).
                // 16 - 24       Real(9.3)     b              b (Angstroms).
                // 25 - 33       Real(9.3)     c              c (Angstroms).
                // 34 - 40       Real(7.2)     alpha          alpha (degrees).
                // 41 - 47       Real(7.2)     beta           beta (degrees).
                // 48 - 54       Real(7.2)     gamma          gamma (degrees).
                // 56 - 66       LString       sGroup         Space  group.
                // 67 - 70       Integer       z              Z value.
                // =============================================================================
                logger.fine(" Crystal record found.");
                if (line.length() < 55) {
                  logger.severe(" CRYST1 record is improperly formatted.");
                }
                double aaxis = parseDouble(line.substring(6, 15).trim());
                double baxis = parseDouble(line.substring(15, 24).trim());
                double caxis = parseDouble(line.substring(24, 33).trim());
                double alpha = parseDouble(line.substring(33, 40).trim());
                double beta = parseDouble(line.substring(40, 47).trim());
                double gamma = parseDouble(line.substring(47, 54).trim());
                int limit = min(line.length(), 66);
                String sg = line.substring(55, limit).trim();
//                properties.clearProperty("a-axis");
//                properties.clearProperty("b-axis");
//                properties.clearProperty("c-axis");
//                properties.clearProperty("alpha");
//                properties.clearProperty("beta");
//                properties.clearProperty("gamma");
//                properties.clearProperty("spacegroup");
//
//                properties.addProperty("a-axis", aaxis);
//                properties.addProperty("b-axis", baxis);
//                properties.addProperty("c-axis", caxis);
//                properties.addProperty("alpha", alpha);
//                properties.addProperty("beta", beta);
//                properties.addProperty("gamma", gamma);
//                properties.addProperty("spacegroup", SpaceGroupInfo.pdb2ShortName(sg));
                Crystal crystal = activeMolecularAssembly.getCrystal();
                SpaceGroup spaceGroup = SpaceGroupDefinitions.spaceGroupFactory(sg);
                if (Objects.equals(crystal.spaceGroup.shortName, spaceGroup.shortName)) {
                  crystal.changeUnitCellParameters(aaxis, baxis, caxis, alpha, beta, gamma);
                } else {
                  // TODO: Handle changes in space groups... Means recalculating force field terms.
                  logger.warning(format(" Original space group %s could not be changed to %s",
                          crystal.spaceGroup.shortName, spaceGroup.shortName));
                }
                break;
              case ENDMDL:
              case END: // Technically speaking, END should be at the end of the file, not end of
                // the model.
                logger.log(Level.FINE, format(" Model %d successfully read", modelsRead));
                modelDone = true;
                break;
              case REMARK:
                remarkLines.add(line.trim());
                if (line.contains("Lambda:")) {
                  Matcher m = lambdaPattern.matcher(line);
                  if (m.find()) {
                    lastReadLambda = Double.parseDouble(m.group(1));
                  }
                }
                break;
              default:
                break;
            }
          } catch (Exception ex) {
            // Do nothing; it's not an ATOM/HETATM line.
          }
          if (modelDone) {
            break;
          }
          line = currentReader.readLine();
        }
        return true;
      } catch (IOException ex) {
        logger.info(
                format(" Exception in parsing frame %d of %s:" + " %s", modelsRead, system.toString(),
                        ex));
      }
    }
    return false;
  }

  /**
   * Specify the alternate location.
   *
   * @param molecularAssembly The MolecularAssembly to populate.
   * @param altLoc The alternate location to use.
   */
  public void setAltID(MolecularAssembly molecularAssembly, Character altLoc) {
    setMolecularSystem(molecularAssembly);
    currentAltLoc = altLoc;
  }

  /**
   * Sets whether this PDBFilter should log each time it saves to a file.
   *
   * @param logWrites a boolean.
   */
  public void setLogWrites(boolean logWrites) {
    this.logWrites = logWrites;
  }

  /**
   * setModelNumbering.
   *
   * @param modelsWritten the number of models written.
   */
  public void setModelNumbering(int modelsWritten) {
    this.modelsWritten = modelsWritten;
  }

  public void setLMN(int[] lmn) {
    if(lmn[0] >= 1 && lmn[1] >= 1 && lmn[2] >= 1){
      this.lmn = lmn;
    }else{
      // Provided dimensions are not handled. Revert to P1.
      this.lmn = new int[]{1,1,1};
    }
  }

  /**
   * setSymOp.
   *
   * @param symOp a int.
   */
  public void setSymOp(int symOp) {
    this.nSymOp = symOp;
  }

  /**
   * Expand the current system to P1 during the save operation.
   *
   * @param file The file to write.
   * @return Return true on a successful write.
   */
  public boolean writeFileAsP1(File file) {
    // XYZ File First Line
    final int l = lmn[0];
    final int m = lmn[1];
    final int n = lmn[2];
    final int numReplicates = l * m * n;
    Crystal crystal = activeMolecularAssembly.getCrystal();
    int nSymOps = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();

    if (nSymOps == 1 && l <= 1 && m <= 1 && n <= 1) {
      // This is a P1 system.
      if (!writeFile(file, false)) {
        logger.info(format(" Save failed for %s", activeMolecularAssembly));
        return false;
      } else {
        return true;
      }
    } else {
      Polymer[] polymers = activeMolecularAssembly.getChains();
      int chainCount = 0;
      for (Polymer polymer : polymers) {
        Character chainID = Polymer.CHAIN_IDS.charAt(chainCount++);
        polymer.setChainID(chainID);
        polymer.setSegID(chainID.toString());
      }
      nSymOp = 0;
      logWrites = false;
      boolean writeEnd = false;
      if (!writeFile(file, false, false, writeEnd)) {
        logger.info(format(" Save failed for %s", activeMolecularAssembly));
        return false;
      } else {
        for (int i = 0; i < l; i++) {
          for (int j = 0; j < m; j++) {
            for (int k = 0; k < n; k++) {
              lValue = i;
              mValue = j;
              nValue = k;
              for (int iSym = 0; iSym < nSymOps; iSym++) {
                nSymOp = iSym;
                for (Polymer polymer : polymers) {
                  Character chainID = Polymer.CHAIN_IDS.charAt(chainCount++);
                  polymer.setChainID(chainID);
                  polymer.setSegID(chainID.toString());
                }
                // If the last sym op to be written.
                writeEnd = iSym == nSymOps - 1 && i == l - 1 && j == m - 1 && k == n - 1;
                if (!writeFile(file, true, false, writeEnd)) {
                  logger.info(format(" Save failed for %s", activeMolecularAssembly));
                  return false;
                }
              }
            }
          }
        }
      }

      // Reset the chainIDs.
      chainCount = 0;
      for (Polymer polymer : polymers) {
        Character chainID = Polymer.CHAIN_IDS.charAt(chainCount++);
        polymer.setChainID(chainID);
        polymer.setSegID(chainID.toString());
      }

    }

    return true;
  }

  /**
   * writeFile
   *
   * @param saveFile a {@link java.io.File} object.
   * @param append Whether to append to saveFile (vs over-write).
   * @param printLinear Ignored (remains to present a different method signature).
   * @param writeEnd True if this is the final model.
   * @return Success of writing.
   */
  public boolean writeFile(File saveFile, boolean append, boolean printLinear, boolean writeEnd) {
    return writeFile(saveFile, append, Collections.emptySet(), writeEnd, true);
  }

  /**
   * writeFile
   *
   * @param saveFile a {@link java.io.File} object to save to.
   * @param append Whether to append to saveFile (vs over-write).
   * @param toExclude A {@link java.util.Set} of {@link ffx.potential.bonded.Atom}s to exclude
   *     from writing.
   * @param writeEnd True if this is the final model.
   * @param versioning True if the file being saved to should be versioned. False if the file
   *     being saved to should be overwritten.
   * @return Success of writing.
   */
  public boolean writeFile(File saveFile, boolean append, Set<Atom> toExclude, boolean writeEnd,
                           boolean versioning) {
    return writeFile(saveFile, append, toExclude, writeEnd, versioning, null);
  }

  /**
   * writeFile
   *
   * @param saveFile a {@link java.io.File} object to save to.
   * @param append Whether to append to saveFile (vs over-write).
   * @param toExclude A {@link java.util.Set} of {@link ffx.potential.bonded.Atom}s to exclude
   *     from writing.
   * @param writeEnd True if this is the final model.
   * @param versioning True if the file being saved to should be versioned. False if the file
   *     being saved to should be overwritten.
   * @param extraLines Extra comment/header lines to write.
   * @return Success of writing.
   */
  public boolean writeFile(File saveFile, boolean append, Set<Atom> toExclude, boolean writeEnd,
                           boolean versioning, String[] extraLines) {
    // Set standardize atom names to false in the presence of deuterium
    List<Atom> deuteriumAtoms = new ArrayList<>();
    for(Atom atom: activeMolecularAssembly.getAtomArray()){
      if(atom.getName().startsWith("D")){
        String name = atom.getName().replace("D","H");
        atom.setName(name);
        deuteriumAtoms.add(atom);
      }
    }
    if (standardizeAtomNames) {
      logger.info(" Setting atom names to PDB standard.");
      renameAtomsToPDBStandard(activeMolecularAssembly);
    }

    for(Atom atom: activeMolecularAssembly.getAtomArray()){
      if(deuteriumAtoms.contains(atom) && atom.getName().startsWith("H")){
        String name = atom.getName().replace("H","D");
        atom.setName(name);
      }
    }
    final Set<Atom> atomExclusions = toExclude == null ? Collections.emptySet() : toExclude;
    if (saveFile == null) {
      return false;
    }
    if (vdwH) {
      logger.info(" Saving hydrogen to van der Waals centers instead of nuclear locations.");
    }
    if (nSymOp > -1) {
      logger.info(format(" Saving atoms using the symmetry operator:\n%s\n",
              activeMolecularAssembly.getCrystal().getUnitCell().spaceGroup.getSymOp(nSymOp)
                      .toString()));
    }

    // Create StringBuilders for ATOM, ANISOU and TER records that can be reused.
    StringBuilder sb = new StringBuilder("ATOM  ");
    StringBuilder anisouSB = new StringBuilder("ANISOU");
    StringBuilder terSB = new StringBuilder("TER   ");
    StringBuilder model = null;
    for (int i = 6; i < 80; i++) {
      sb.append(' ');
      anisouSB.append(' ');
      terSB.append(' ');
    }

    File newFile = saveFile;
    if (!append) {
      if (versioning) {
        newFile = version(saveFile);
      }
    } else if (modelsWritten >= 0) {
      model = new StringBuilder(format("MODEL     %-4d", ++modelsWritten));
      model.append(repeat(" ", 65));
    }
    activeMolecularAssembly.setFile(newFile);
    if (activeMolecularAssembly.getName() == null) {
      activeMolecularAssembly.setName(newFile.getName());
    }
    if (logWrites) {
      logger.log(Level.INFO, " Saving {0}", newFile.getName());
    }

    try (FileWriter fw = new FileWriter(newFile, append);
         BufferedWriter bw = new BufferedWriter(fw)) {
      /*
       Will come before CRYST1 and ATOM records, but after anything
       written by writeFileWithHeader (particularly X-ray refinement statistics).
      */
      String[] headerLines = activeMolecularAssembly.getHeaderLines();
      for (String line : headerLines) {
        bw.write(format("%s\n", line));
      }
      if (extraLines != null) {
        if (rotamerTitration && extraLines[0].contains("REMARK")) {
          for (String line : extraLines) {
            bw.write(line + "\n");
          }
        } else {
          for (String line : extraLines) {
            bw.write(format("REMARK 999 %s\n", line));
          }
        }
      }
      if (model != null) {
        bw.write(model.toString());
        bw.newLine();
      }
      // =============================================================================
      // The CRYST1 record presents the unit cell parameters, space group, and Z
      // value. If the structure was not determined by crystallographic means, CRYST1
      // simply provides the unitary values, with an appropriate REMARK.
      //
      //  7 - 15       Real(9.3)     a              a (Angstroms).
      // 16 - 24       Real(9.3)     b              b (Angstroms).
      // 25 - 33       Real(9.3)     c              c (Angstroms).
      // 34 - 40       Real(7.2)     alpha          alpha (degrees).
      // 41 - 47       Real(7.2)     beta           beta (degrees).
      // 48 - 54       Real(7.2)     gamma          gamma (degrees).
      // 56 - 66       LString       sGroup         Space  group.
      // 67 - 70       Integer       z              Z value.
      // =============================================================================
      if (nSymOp < 0) {
        // Write out the unit cell.
        Crystal crystal = activeMolecularAssembly.getCrystal();
        if (crystal != null && !crystal.aperiodic()) {
          Crystal c = crystal.getUnitCell();
          if (lmn[0] > 0 || lmn[1] > 0 || lmn[2] > 0) {
            c.a = c.a * lmn[0];
            c.b = c.b * lmn[1];
            c.c = c.c * lmn[2];
          }
          bw.write(c.toCRYST1());
        }
      } else if (nSymOp == 0) {
        // Write a P1 cell.
        Crystal crystal = activeMolecularAssembly.getCrystal();
        if (crystal != null && !crystal.aperiodic()) {
          Crystal c = crystal.getUnitCell();
          Crystal p1 = new Crystal((lmn[0]>0)? c.a * lmn[0] : c.a, (lmn[1]>0)? c.b * lmn[1] : c.b, (lmn[2]>0)? c.c * lmn[2] : c.c, c.alpha, c.beta, c.gamma, "P1");
          bw.write(p1.toCRYST1());
        }
      }
      // =============================================================================
      // The SSBOND record identifies each disulfide bond in protein and polypeptide
      // structures by identifying the two residues involved in the bond.
      // The disulfide bond distance is included after the symmetry operations at
      // the end of the SSBOND record.
      //
      //  8 - 10        Integer         serNum       Serial number.
      // 12 - 14        LString(3)      "CYS"        Residue name.
      // 16             Character       chainID1     Chain identifier.
      // 18 - 21        Integer         seqNum1      Residue sequence number.
      // 22             AChar           icode1       Insertion code.
      // 26 - 28        LString(3)      "CYS"        Residue name.
      // 30             Character       chainID2     Chain identifier.
      // 32 - 35        Integer         seqNum2      Residue sequence number.
      // 36             AChar           icode2       Insertion code.
      // 60 - 65        SymOP           sym1         Symmetry oper for 1st resid
      // 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
      // 74 – 78        Real(5.2)      Length        Disulfide bond distance
      //
      // If SG of cysteine is disordered then there are possible alternate linkages.
      // wwPDB practice is to put together all possible SSBOND records. This is
      // problematic because the alternate location identifier is not specified in
      // the SSBOND record.
      // =============================================================================
      int serNum = 1;
      Polymer[] polymers = activeMolecularAssembly.getChains();
      if (polymers != null) {
        for (Polymer polymer : polymers) {
          List<Residue> residues = polymer.getResidues();
          for (Residue residue : residues) {
            if (residue.getName().equalsIgnoreCase("CYS")) {
              List<Atom> cysAtoms = residue.getAtomList().stream()
                      .filter(a -> !atomExclusions.contains(a)).toList();
              Atom SG1 = null;
              for (Atom atom : cysAtoms) {
                String atName = atom.getName().toUpperCase();
                if (atName.equals("SG") || atName.equals("SH")
                        || atom.getAtomType().atomicNumber == 16) {
                  SG1 = atom;
                  break;
                }
              }
              List<Bond> bonds = SG1.getBonds();
              for (Bond bond : bonds) {
                Atom SG2 = bond.get1_2(SG1);
                if (SG2.getAtomType().atomicNumber == 16 && !atomExclusions.contains(SG2)) {
                  if (SG1.getIndex() < SG2.getIndex()) {
                    bond.energy(false);
                    bw.write(format("SSBOND %3d CYS %1s %4s    CYS %1s %4s %36s %5.2f\n", serNum++,
                            SG1.getChainID().toString(), Hybrid36.encode(4, SG1.getResidueNumber()),
                            SG2.getChainID().toString(), Hybrid36.encode(4, SG2.getResidueNumber()), "",
                            bond.getValue()));
                  }
                }
              }
            }
          }
        }
      }
      // =============================================================================
      //
      //  7 - 11        Integer       serial       Atom serial number.
      // 13 - 16        Atom          name         Atom name.
      // 17             Character     altLoc       Alternate location indicator.
      // 18 - 20        Residue name  resName      Residue name.
      // 22             Character     chainID      Chain identifier.
      // 23 - 26        Integer       resSeq       Residue sequence number.
      // 27             AChar         iCode        Code for insertion of residues.
      // 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
      // 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
      // 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
      // 55 - 60        Real(6.2)     occupancy    Occupancy.
      // 61 - 66        Real(6.2)     tempFactor   Temperature factor.
      // 77 - 78        LString(2)    element      Element symbol, right-justified.
      // 79 - 80        LString(2)    charge       Charge  on the atom.
      // =============================================================================
      //         1         2         3         4         5         6         7
      // 123456789012345678901234567890123456789012345678901234567890123456789012345678
      // ATOM      1  N   ILE A  16      60.614  71.140 -10.592  1.00  7.38           N
      // ATOM      2  CA  ILE A  16      60.793  72.149  -9.511  1.00  6.91           C
      MolecularAssembly[] molecularAssemblies = this.getMolecularAssemblyArray();
      int serial = 1;
      if (nSymOp > 0) {
        serial = serialP1;
      }

      // Loop over biomolecular chains
      if (polymers != null) {
        for (Polymer polymer : polymers) {
          currentSegID = polymer.getName();
          currentChainID = polymer.getChainID();
          sb.setCharAt(21, currentChainID);
          // Loop over residues
          List<Residue> residues = polymer.getResidues();
          for (Residue residue : residues) {
            String resName = residue.getName();
            if (resName.length() > 3) {
              resName = resName.substring(0, 3);
            }
            int resID = residue.getResidueNumber();
            sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
            sb.replace(22, 26, format("%4s", Hybrid36.encode(4, resID)));
            // Loop over atoms
            List<Atom> residueAtoms = residue.getAtomList().stream()
                    .filter(a -> !atomExclusions.contains(a)).collect(Collectors.toList());
            boolean altLocFound = false;
            for (Atom atom : residueAtoms) {
              if (mutate) {
                for (Mutation mtn : mutations) {
                  if (resID == mtn.resID) {
                    ArrayList<String> alchAtoms = mtn.getAlchemicalAtoms(true);
                    if (alchAtoms != null) {
                      if (residue.getBackboneAtoms().contains(atom) && alchAtoms.contains(atom.getName())) {
                        logger.info(format(" MUTATION atom is %d chain %s",serial, currentChainID));
                      }
                    } else {
                      // treating pur-pyr or pyr-pur N9/N1 & C2/C4 as alchemical
                      if (residue.getBackboneAtoms().contains(atom)) {
                        logger.info(format(" MUTATION atom is %d chain %s",serial, currentChainID));
                      }
                    }
                  }
                }
              }
              writeAtom(atom, serial++, sb, anisouSB, bw);
              Character altLoc = atom.getAltLoc();
              if (altLoc != null && !altLoc.equals(' ')) {
                altLocFound = true;
              }
            }
            // Write out alternate conformers
            if (altLocFound) {
              for (int ma = 1; ma < molecularAssemblies.length; ma++) {
                MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
                Polymer altPolymer = altMolecularAssembly.getPolymer(currentChainID, currentSegID, false);
                Residue altResidue = altPolymer.getResidue(resName, resID, false, Residue.ResidueType.AA);
                if (altResidue == null) {
                  resName = AminoAcid3.UNK.name();
                  altResidue = altPolymer.getResidue(resName, resID, false, Residue.ResidueType.AA);
                }
                residueAtoms = altResidue.getAtomList().stream()
                        .filter(a -> !atomExclusions.contains(a)).collect(Collectors.toList());
                for (Atom atom : residueAtoms) {
                  if (atom.getAltLoc() != null && !atom.getAltLoc().equals(' ') && !atom.getAltLoc()
                          .equals('A')) {
                    sb.replace(17, 20, padLeft(atom.getResidueName().toUpperCase(), 3));
                    writeAtom(atom, serial++, sb, anisouSB, bw);
                  }
                }
              }
            }
          }
          terSB.replace(6, 11, format("%5s", Hybrid36.encode(5, serial++)));
          terSB.replace(12, 16, "    ");
          terSB.replace(16, 26, sb.substring(16, 26));
          bw.write(terSB.toString());
          bw.newLine();
        }
      }
      sb.replace(0, 6, "HETATM");
      sb.setCharAt(21, 'A');

      Character chainID = 'A';
      if (polymers != null) {
        chainID = polymers[0].getChainID();
      }
      activeMolecularAssembly.setChainIDAndRenumberMolecules(chainID);

      // Loop over molecules, ions and then water.
      List<MSNode> molecules = activeMolecularAssembly.getMolecules();
      int numMolecules = molecules.size();
      for (int i = 0; i < numMolecules; i++) {
        Molecule molecule = (Molecule) molecules.get(i);
        chainID = molecule.getChainID();
        sb.setCharAt(21, chainID);
        String resName = molecule.getResidueName();
        int resID = molecule.getResidueNumber();
        if (resName.length() > 3) {
          resName = resName.substring(0, 3);
        }
        sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
        sb.replace(22, 26, format("%4s", Hybrid36.encode(4, resID)));
        // List<Atom> moleculeAtoms = molecule.getAtomList();
        List<Atom> moleculeAtoms = molecule.getAtomList().stream()
                .filter(a -> !atomExclusions.contains(a)).collect(Collectors.toList());
        boolean altLocFound = false;
        for (Atom atom : moleculeAtoms) {
          writeAtom(atom, serial++, sb, anisouSB, bw);
          Character altLoc = atom.getAltLoc();
          if (altLoc != null && !altLoc.equals(' ')) {
            altLocFound = true;
          }
        }
        // Write out alternate conformers
        if (altLocFound) {
          for (int ma = 1; ma < molecularAssemblies.length; ma++) {
            MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
            MSNode altmolecule = altMolecularAssembly.getMolecules().get(i);
            moleculeAtoms = altmolecule.getAtomList();
            for (Atom atom : moleculeAtoms) {
              if (atom.getAltLoc() != null && !atom.getAltLoc().equals(' ') && !atom.getAltLoc()
                      .equals('A')) {
                writeAtom(atom, serial++, sb, anisouSB, bw);
              }
            }
          }
        }
      }

      List<MSNode> ions = activeMolecularAssembly.getIons();
      for (int i = 0; i < ions.size(); i++) {
        Molecule ion = (Molecule) ions.get(i);
        chainID = ion.getChainID();
        sb.setCharAt(21, chainID);
        String resName = ion.getResidueName();
        int resID = ion.getResidueNumber();
        if (resName.length() > 3) {
          resName = resName.substring(0, 3);
        }
        sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
        sb.replace(22, 26, format("%4s", Hybrid36.encode(4, resID)));
        // List<Atom> ionAtoms = ion.getAtomList();
        List<Atom> ionAtoms = ion.getAtomList().stream().filter(a -> !atomExclusions.contains(a))
                .collect(Collectors.toList());
        boolean altLocFound = false;
        for (Atom atom : ionAtoms) {
          writeAtom(atom, serial++, sb, anisouSB, bw);
          Character altLoc = atom.getAltLoc();
          if (altLoc != null && !altLoc.equals(' ')) {
            altLocFound = true;
          }
        }
        // Write out alternate conformers
        if (altLocFound) {
          for (int ma = 1; ma < molecularAssemblies.length; ma++) {
            MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
            MSNode altion = altMolecularAssembly.getIons().get(i);
            ionAtoms = altion.getAtomList();
            for (Atom atom : ionAtoms) {
              if (atom.getAltLoc() != null && !atom.getAltLoc().equals(' ') && !atom.getAltLoc()
                      .equals('A')) {
                writeAtom(atom, serial++, sb, anisouSB, bw);
              }
            }
          }
        }
      }

      List<MSNode> water = activeMolecularAssembly.getWater();
      for (int i = 0; i < water.size(); i++) {
        Molecule wat = (Molecule) water.get(i);
        chainID = wat.getChainID();
        sb.setCharAt(21, chainID);
        String resName = wat.getResidueName();
        int resID = wat.getResidueNumber();
        if (resName.length() > 3) {
          resName = resName.substring(0, 3);
        }
        sb.replace(17, 20, padLeft(resName.toUpperCase(), 3));
        sb.replace(22, 26, format("%4s", Hybrid36.encode(4, resID)));
        List<Atom> waterAtoms = wat.getAtomList().stream().filter(a -> !atomExclusions.contains(a))
                .collect(Collectors.toList());
        boolean altLocFound = false;
        for (Atom atom : waterAtoms) {
          writeAtom(atom, serial++, sb, anisouSB, bw);
          Character altLoc = atom.getAltLoc();
          if (altLoc != null && !altLoc.equals(' ')) {
            altLocFound = true;
          }
        }
        // Write out alternate conformers
        if (altLocFound) {
          for (int ma = 1; ma < molecularAssemblies.length; ma++) {
            MolecularAssembly altMolecularAssembly = molecularAssemblies[ma];
            MSNode altwater = altMolecularAssembly.getWater().get(i);
            waterAtoms = altwater.getAtomList();
            for (Atom atom : waterAtoms) {
              if (atom.getAltLoc() != null && !atom.getAltLoc().equals(' ') && !atom.getAltLoc()
                      .equals('A')) {
                writeAtom(atom, serial++, sb, anisouSB, bw);
              }
            }
          }
        }
      }

      if (model != null) {
        bw.write("ENDMDL");
        bw.newLine();
      }

      if (writeEnd) {
        bw.write("END");
        bw.newLine();
      }

      if (nSymOp >= 0) {
        serialP1 = serial;
      }

    } catch (Exception e) {
      String message = "Exception writing to file: " + saveFile;
      logger.log(Level.WARNING, message, e);
      return false;
    }
    return true;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Write out the Atomic information in PDB format.
   */
  @Override
  public boolean writeFile(File saveFile, boolean append) {
    return writeFile(saveFile, append, false, true);
  }

  public boolean writeFile(File saveFile, boolean append, String[] extraLines) {
    return writeFile(saveFile, append, Collections.emptySet(), false, !append, extraLines);
  }

  /**
   * Writes out the atomic information in PDB format.
   *
   * @param saveFile The file to save information to.
   * @param append True if the current data should be appended to the saveFile (as in arc
   *     files).
   * @param versioning True if the saveFile should be versioned. False if the saveFile should be
   *     overwritten.
   * @return Success of writing.
   */
  public boolean writeFile(File saveFile, boolean append, boolean versioning) {
    return writeFile(saveFile, append, Collections.emptySet(), true, versioning);
  }

  /**
   * writeFileWithHeader.
   *
   * @param saveFile a {@link java.io.File} object.
   * @param header a {@link java.lang.String} object.
   * @param append a boolean.
   * @return a boolean.
   */
  public boolean writeFileWithHeader(File saveFile, String header, boolean append) {
    // Set standardize atom names to false in the presence of deuterium
    List<Atom> deuteriumAtoms = new ArrayList<>();
    for(Atom atom: activeMolecularAssembly.getAtomArray()){
      if(atom.getName().startsWith("D")){
        String name = atom.getName().replace("D","H");
        atom.setName(name);
        deuteriumAtoms.add(atom);
      }
    }
    if (standardizeAtomNames) {
      logger.info(" Setting atom names to PDB standard.");
      renameAtomsToPDBStandard(activeMolecularAssembly);
    }

    for(Atom atom: activeMolecularAssembly.getAtomArray()){
      if(deuteriumAtoms.contains(atom) && atom.getName().startsWith("H")){
        String name = atom.getName().replace("H","D");
        atom.setName(name);
      }
    }
    activeMolecularAssembly.setFile(saveFile);
    activeMolecularAssembly.setName(saveFile.getName());

    try (FileWriter fw = new FileWriter(saveFile, append); BufferedWriter bw = new BufferedWriter(
            fw)) {
      bw.write(header);
      bw.newLine();
    } catch (Exception e) {
      String message = " Exception writing to file: " + saveFile;
      logger.log(Level.WARNING, message, e);
      return false;
    }
    if (writeFile(saveFile, true)) {
      logger.log(Level.INFO, " Wrote PDB to file {0}", saveFile.getPath());
      return true;
    } else {
      logger.log(Level.INFO, " Error writing to file {0}", saveFile.getPath());
      return false;
    }
  }

  /**
   * writeFileWithHeader.
   *
   * @param saveFile a {@link java.io.File} object.
   * @param header a {@link java.lang.String} object.
   * @return a boolean.
   */
  public boolean writeFileWithHeader(File saveFile, String header) {
    return writeFileWithHeader(saveFile, header, true);
  }

  /**
   * writeFileWithHeader.
   *
   * @param saveFile a {@link java.io.File} object.
   * @param header a {@link java.lang.StringBuilder} object.
   * @return a boolean.
   */
  public boolean writeFileWithHeader(File saveFile, StringBuilder header) {
    return writeFileWithHeader(saveFile, header.toString());
  }

  /**
   * Get unique SegID for possibly duplicate chain IDs.
   *
   * @param c chain ID just read.
   * @return a unique segID.
   */
  private String getExistingSegID(Character c) {
    if (c.equals(' ')) {
      c = 'A';
    }

    // If the chain ID has not changed, return the existing segID.
    if (c.equals(currentChainID)) {
      return currentSegID;
    } else {
      currentChainID = null;
    }

    List<String> segIDs = segidMap.get(c);
    if (segIDs != null) {
      if (segIDs.size() > 1) {
        if (currentSegID == null) {
          currentChainID = c;
          currentSegID = segIDs.get(0);
          return segIDs.get(0);
        } else if (currentSegID.length() == 1) {
          currentChainID = c;
          currentSegID = segIDs.get(1);
          return segIDs.get(1);
        } else if (currentSegID.length() == 2) {
          String s = currentSegID.substring(0,1);
          int num = -2;
          try {
            num = Integer.parseInt(s);
          } catch (NumberFormatException e) {
            logger.severe(" SegID of length 2 does not start with an integer.");
          }
          currentChainID = c;
          currentSegID = segIDs.get(num+1);
          return segIDs.get(num+1);
        } else {
          logger.info(" Too many repeated chains. Using single letter for segID.");
        }
      }
      return segIDs.get(0);
    } else {
      logger.log(Level.INFO, format(" Creating SegID for to chain %s", c));
      return getSegID(c);
    }
  }

  /**
   * Convert possibly duplicate chain IDs into unique segIDs.
   *
   * @param c chain ID just read.
   * @return a unique segID.
   */
  private String getSegID(Character c) {
    if (c.equals(' ')) {
      c = 'A';
    }

    // If the chain ID has not changed, return the existing segID.
    if (c.equals(currentChainID)) {
      return currentSegID;
    }

    // Loop through existing segIDs to find the first one that is unused.
    int count = 0;
    for (String segID : segIDs) {
      if (segID.endsWith(c.toString())) {
        count++;
      }
    }

    // If the count is greater than 0, then append it.
    String newSegID;
    if (count == 0) {
      newSegID = c.toString();
    } else {
      newSegID = count + c.toString();
    }
    segIDs.add(newSegID);
    currentChainID = c;
    currentSegID = newSegID;

    if (segidMap.containsKey(c)) {
      segidMap.get(c).add(newSegID);
    } else {
      List<String> newChainList = new ArrayList<>();
      newChainList.add(newSegID);
      segidMap.put(c, newChainList);
    }

    return newSegID;
  }

  /**
   * writeAtom
   *
   * @param atom a {@link ffx.potential.bonded.Atom} object.
   * @param serial a int.
   * @param sb a {@link java.lang.StringBuilder} object.
   * @param anisouSB a {@link java.lang.StringBuilder} object.
   * @param bw a {@link java.io.BufferedWriter} object.
   * @throws java.io.IOException if any.
   */
  private void writeAtom(Atom atom, int serial, StringBuilder sb, StringBuilder anisouSB,
                         BufferedWriter bw) throws IOException {
    String name = atom.getName();
    int nameLength = name.length();
    if (nameLength > 4) {
      name = name.substring(0, 4);
    } else if (nameLength == 1) {
      name = name + "  ";
    } else if (nameLength == 2) {
      if (atom.getAtomType().valence == 0) {
        name = name + "  ";
      } else {
        name = name + " ";
      }
    }
    double[] xyz = vdwH ? atom.getRedXYZ() : atom.getXYZ(null);
    if (nSymOp >= 0) {
      Crystal crystal = activeMolecularAssembly.getCrystal().getUnitCell();
      SymOp symOp = crystal.spaceGroup.getSymOp(nSymOp);
      double[] newXYZ = new double[xyz.length];
      crystal.applySymOp(xyz, newXYZ, symOp);
      if (lValue > 0 || mValue > 0 || nValue > 0) {
        double[] translation = new double[] {lValue, mValue, nValue};
        crystal.getUnitCell().toCartesianCoordinates(translation, translation);
        newXYZ[0] += translation[0];
        newXYZ[1] += translation[1];
        newXYZ[2] += translation[2];
      }
      xyz = newXYZ;
    }
    sb.replace(6, 16, format("%5s " + padLeft(name.toUpperCase(), 4), Hybrid36.encode(5, serial)));
    Character altLoc = atom.getAltLoc();
    sb.setCharAt(16, Objects.requireNonNullElse(altLoc, ' '));


    /*
     * On the following code:
     * #1: StringBuilder.replace will allow for longer strings, expanding the StringBuilder's length if necessary.
     *
     * #2: sb was never re-initialized, so if there was overflow,
     * sb would continue to be > 80 characters long, resulting in broken PDB files
     *
     * #3: It may be wiser to have XYZ coordinates result in shutdown, not
     * truncation of coordinates. #4: Excessive B-factors aren't much of an
     * issue; if the B-factor is past 999.99, that's the difference between
     * "density extends to Venus" and "density extends to Pluto".
     */
    StringBuilder decimals = new StringBuilder();
    for (int i = 0; i < 3; i++) {
      try {
        decimals.append(StringUtils.fwFpDec(xyz[i], 8, 3));
      } catch (IllegalArgumentException ex) {
        String newValue = StringUtils.fwFpTrunc(xyz[i], 8, 3);
        logger.info(format(" XYZ %d coordinate %8.3f for atom %s "
                + "overflowed bounds of 8.3f string specified by PDB "
                + "format; truncating value to %s", i, xyz[i], atom, newValue));
        decimals.append(newValue);
      }
    }
    try {
      decimals.append(StringUtils.fwFpDec(atom.getOccupancy(), 6, 2));
    } catch (IllegalArgumentException ex) {
      logger.severe(
              format(" Occupancy %f for atom %s is impossible; " + "value must be between 0 and 1",
                      atom.getOccupancy(), atom));
    }
    try {
      decimals.append(StringUtils.fwFpDec(atom.getTempFactor(), 6, 2));
    } catch (IllegalArgumentException ex) {
      String newValue = StringUtils.fwFpTrunc(atom.getTempFactor(), 6, 2);
      logger.info(format(" Atom temp factor %6.2f for atom %s overflowed "
                      + "bounds of 6.2f string specified by PDB format; truncating " + "value to %s",
              atom.getTempFactor(), atom, newValue));
      decimals.append(newValue);
    }
    sb.replace(30, 66, decimals.toString());

    name = Atom.ElementSymbol.values()[atom.getAtomicNumber() - 1].toString();
    name = name.toUpperCase();
    if (atom.isDeuterium()) {
      name = "D";
    }
    sb.replace(76, 78, padLeft(name, 2));
    sb.replace(78, 80, format("%2d", 0));
    bw.write(sb.toString());
    bw.newLine();
    // =============================================================================
    //  1 - 6        Record name   "ANISOU"
    //  7 - 11       Integer       serial         Atom serial number.
    // 13 - 16       Atom          name           Atom name.
    // 17            Character     altLoc         Alternate location indicator
    // 18 - 20       Residue name  resName        Residue name.
    // 22            Character     chainID        Chain identifier.
    // 23 - 26       Integer       resSeq         Residue sequence number.
    // 27            AChar         iCode          Insertion code.
    // 29 - 35       Integer       u[0][0]        U(1,1)
    // 36 - 42       Integer       u[1][1]        U(2,2)
    // 43 - 49       Integer       u[2][2]        U(3,3)
    // 50 - 56       Integer       u[0][1]        U(1,2)
    // 57 - 63       Integer       u[0][2]        U(1,3)
    // 64 - 70       Integer       u[1][2]        U(2,3)
    // 77 - 78       LString(2)    element        Element symbol, right-justified.
    // 79 - 80       LString(2)    charge         Charge on the atom.
    // =============================================================================
    double[] anisou = atom.getAnisou(null);
    if (anisou != null) {
      anisouSB.replace(6, 80, sb.substring(6, 80));
      anisouSB.replace(28, 70,
              format("%7d%7d%7d%7d%7d%7d", (int) (anisou[0] * 1e4), (int) (anisou[1] * 1e4),
                      (int) (anisou[2] * 1e4), (int) (anisou[3] * 1e4), (int) (anisou[4] * 1e4),
                      (int) (anisou[5] * 1e4)));
      bw.write(anisouSB.toString());
      bw.newLine();
    }
  }

  /** PDB records that are recognized. */
  private enum Record {
    ANISOU, ATOM, CONECT, CRYST1, DBREF, END, MODEL, ENDMDL, HELIX, HETATM, LINK, MTRIX1, MTRIX2, MTRIX3, MODRES, SEQRES, SHEET, SSBOND, REMARK
  }

  /** Presently, VERSION3_3 is default, and VERSION3_2 is anything non-standard. */
  public enum PDBFileStandard {
    VERSION2_3, VERSION3_0, VERSION3_1, VERSION3_2, VERSION3_3
  }

  public static class Mutation {

    /** Residue ID of the residue to mutate. */
    final int resID;
    /** Residue name after mutation. */
    final String resName;
    /** Character for the chain ID of the residue that will be mutated. */
    final char chainChar;
    /** Residue name before mutation. */
    String origResName;

    public Mutation(int resID, char chainChar, String newResName) {
      newResName = newResName.toUpperCase();
      if (newResName.length() != 3) {
        logger.log(Level.WARNING, format("Invalid mutation target: %s.", newResName));
      }
      int isAA = AminoAcidUtils.getAminoAcidNumber(newResName);
      int isNA = NucleicAcidUtils.getNucleicAcidNumber(newResName);
      if (isAA == -1 && isNA == -1) {
        logger.log(Level.WARNING, format("Invalid mutation target: %s.", newResName));
      }
      this.resID = resID;
      this.chainChar = chainChar;
      this.resName = newResName;
    }

    /**
     * Check to see if an atom is involved in the mutated base's glycosyl torsion. If the mutation is a switch from
     * purine to pyrimidine or vice versa, it will return '~name', meaning the name should be replaced but to include it
     * as an alchemical atom.
     * @param atomName atom name to check
     * @return new name to use if it is involved in glycosyl torsion, null otherwise
     */
    public String isNonAlchemicalAtom(String atomName) {
      if (isWtPurine()) {
        if (atomName.equals("N9")) {
          if (isMtnPyrimidine()) {
            return "~N1";
          }
          return atomName;
        } else if (atomName.equals("C4")) {
          if (isMtnPyrimidine()) {
            return "~C2";
          }
          return atomName;
        }
        return null;
      }

      if (isWtPyrimidine()) {
        if (atomName.equals("N1")) {
          if (isMtnPurine()) {
            // here
            return "~N9";
          }
          return atomName;
        } else if (atomName.equals("C2")) {
          if (isMtnPurine()) {
            // here
            return "~C4";
          }
          return atomName;
        }
        return null;
      }

      return null;
    }

    /**
     * Determines what atoms should be alchemical for a purine to purine or pyrimidine to pyrimidine mutation.
     * @param isWriting true if writing the pdb, false if reading the pdb
     * @return ArrayList of alchemical atoms, null if not a pur-pur or pyr-pyr mutation
     */
    public ArrayList<String> getAlchemicalAtoms(boolean isWriting) {
      // Log warning that the mutation input is the same residue and return nothing so prev. functionality is not changed
      if (resName.equals(origResName)) {
        logger.severe("Desired Mutation residue is the same as the original.");
        return null;
      }

      boolean purpur;
      if (isMtnPurine() && isWtPurine()) {
        purpur = true;
      } else if (isMtnPyrimidine() && isWtPyrimidine()) {
        purpur = false;
      } else {
        // Return nothing so previous functionality is not changed
        return null;
      }

      String res;

      // look at the mutation residue if writing or the original (wild type) residue if reading
      if (isWriting) {
        res = resName;
      } else {
        res = origResName;
      }

      ArrayList<String> list = new ArrayList<>();
      if (purpur) { // purine: either A to G or G to A
        if (res.equals("DAD") || res.equals("DA")) { // MTN is A
          list.add("N6");
          list.add("H61");
          list.add("H62");
          list.add("H2");
        } else { // MTN is G
          list.add("H1");
          list.add("N2");
          list.add("H21");
          list.add("H22");
          list.add("O6");
        }
      } else { // pyrimidine: either T to C or C to T
        if (res.equals("DTY") || res.equals("DT")) { // MTN is T
          list.add("H3");
          list.add("O4");
          list.add("C7");
          list.add("H71");
          list.add("H72");
          list.add("H73");
        } else { // MTN is C
          list.add("N4");
          list.add("H41");
          list.add("H42");
          list.add("H5");
        }
      }
      return list;
    }

    /**
     * Determine if the mutation residue is purine.
     * @return true if mutation is purine
     */
    public boolean isMtnPurine() {
      return resName.equals("DA") || resName.equals("DG") || resName.equals("DAD") || resName.equals("DGU");
    }

    /**
     * Determine if the mutation residue is pyrimidine.
     * @return true if mutation is pyrimidine
     */
    public boolean isMtnPyrimidine() {
      return resName.equals("DC") || resName.equals("DT") || resName.equals("DCY") || resName.equals("DTY");
    }

    /**
     * Determine if original (wild type) residue is purine.
     * @return true if original residue is purine
     */
    public boolean isWtPurine() {
      return origResName.equals("DA") || origResName.equals("DG") || origResName.equals("DAD") || origResName.equals("DGU");
    }

    /**
     * Determine if original (wild type) residue is pyrimidine.
     * @return true if original residue is pyrimidine
     */
    public boolean isWtPyrimidine() {
      return origResName.equals("DT") || origResName.equals("DC") || origResName.equals("DTY") || origResName.equals("DCY");
    }
  }
}