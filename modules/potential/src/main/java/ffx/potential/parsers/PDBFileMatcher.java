/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.potential.parsers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.PDBCrystallographicInfo;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.SSBond;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.StructurePairAligner;
import org.biojava.nbio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.nbio.structure.io.PDBFileReader;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import org.biojava.nbio.structure.xtal.CrystalCell;

/**
 * Aligns a list of files with a list of source files by RMSD, and can use the
 * source files to fix certain issues in the aligned files.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 */
public class PDBFileMatcher {

    private static final Logger logger = Logger.getLogger(PDBFileMatcher.class.getName());
    private final File[] sourceFiles;
    private final FileFilePair[] matchFilePairs;
    private String fileSuffix = "_match";
    private boolean verbose = false;
    private boolean parallel = true;
    private boolean headerLink = true;
    private boolean fixBFactors = false;
    private boolean fixSSBonds = false;
    private boolean fixAtoms = false;
    private boolean superpose = true;
    private int atomsUsed = 1;
    private ParallelTeam parallelTeam;
    private long[] iterationTimes;
    private long[] fixTimes;
    private boolean robustMatch = false;
    private boolean fixCryst = false;
    private boolean fixModel = false;

    public PDBFileMatcher(File[] sourceFiles, File[] matchFiles) {
        this.sourceFiles = sourceFiles;
        int numMatchFiles = matchFiles.length;
        matchFilePairs = new FileFilePair[numMatchFiles];
        iterationTimes = new long[numMatchFiles];
        fixTimes = new long[numMatchFiles];
        for (int i = 0; i < numMatchFiles; i++) {
            matchFilePairs[i] = new FileFilePair(matchFiles[i]);
        }
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    public void setParallel(boolean parallel) {
        this.parallel = parallel;
    }

    public void setHeaderLink(boolean headerLink) {
        this.headerLink = headerLink;
    }

    public void setFixBFactors(boolean fixBFactors) {
        this.fixBFactors = fixBFactors;
    }

    public void setFixSSBonds(boolean fixSSBonds) {
        this.fixSSBonds = fixSSBonds;
    }

    public void setSuperpose(boolean superpose) {
        this.superpose = superpose;
    }

    public void setAtomsUsed(int atomsUsed) {
        this.atomsUsed = atomsUsed;
    }

    public void setSuffix(String suffix) {
        this.fileSuffix = suffix;
    }

    public void setRobustMatch(boolean robustMatch) {
        this.robustMatch = robustMatch;
    }

    public void setFixCryst(boolean fixCryst) {
        this.fixCryst = fixCryst;
    }

    /**
     * Primary class method: matches match files to source files, and if any
     * fixer flags are set, calls fixFiles to edit the match files.
     */
    public void match() {
        try {
            if (parallel) {
                parallelTeam = new ParallelTeam();
                try {
                    parallelTeam.execute(new MatchingRegion());
                } catch (Exception ex) {
                    logger.severe(String.format(" Exception matching files in parallel: %s", ex.toString()));
                }
                for (int i = 0; i < iterationTimes.length; i++) {
                    logger.info(String.format(" Iteration %d time: %12.9f sec", i, 1.0E-9 * iterationTimes[i]));
                }
            } else {
                sequentialFileMatch();
            }
            fixAtoms = fixSSBonds || fixBFactors;
            fixModel = fixAtoms || headerLink || fixCryst;
            if (fixModel) {
                fixFiles();
            }
        } catch (Exception ex) {
            logger.severe(String.format(" Error in matching: %s", ex.toString()));
        }
    }

    /**
     * Second major class method: uses information from source files to fix
     * match files.
     */
    private void fixFiles() {
        if (parallel) {
            if (parallelTeam == null) {
                parallelTeam = new ParallelTeam();
            }
            try {
                parallelTeam.execute(new FixerRegion());
            } catch (Exception ex) {
                logger.severe(String.format(" Exception fixing files in parallel: %s", ex.toString()));
            }
        } else {
            sequentialFixer();
        }
        for (int i = 0; i < fixTimes.length; i++) {
            logger.info(String.format(" File %d fix time: %12.9f sec", i, 1.0E-9 * fixTimes[i]));
        }
    }

    private double calculateRMSD(FileFilePair matchFile, Structure sourceStructure, StructurePairAligner aligner) throws StructureException {
        Structure matchStructure = matchFile.getStructure();

        if (superpose) {
            double rmsd = Double.MAX_VALUE;
            aligner.align(matchStructure, matchStructure);
            AlternativeAlignment[] alignments = aligner.getAlignments();
            for (AlternativeAlignment alignment : alignments) {
                double alignRMSD = alignment.getRmsd();
                rmsd = alignRMSD < rmsd ? alignRMSD : rmsd;
            }
            return rmsd;
        }

        Atom[] matchArray;
        Atom[] sourceArray;
        switch (atomsUsed) {
            case 1:
                String[] atomNames = {"CA", "N1", "N9"};
                matchArray = StructureTools.getAtomArray(matchStructure, atomNames);
                sourceArray = StructureTools.getAtomArray(sourceStructure, atomNames);
                break;

            case 2:
                List<Atom> matchAtoms = new ArrayList<>();
                List<Chain> matchChains = matchStructure.getChains();
                for (Chain chain : matchChains) {
                    List<Group> matchGroups = chain.getAtomGroups(GroupType.AMINOACID);
                    matchGroups.addAll(chain.getAtomGroups(GroupType.NUCLEOTIDE));
                    for (Group group : matchGroups) {
                        matchAtoms.addAll(group.getAtoms());
                    }
                }
                matchArray = matchAtoms.toArray(new Atom[matchAtoms.size()]);

                List<Atom> sourceAtoms = new ArrayList<>();
                List<Chain> sourceChains = sourceStructure.getChains();
                for (Chain chain : sourceChains) {
                    List<Group> sourceGroups = chain.getAtomGroups(GroupType.AMINOACID);
                    sourceGroups.addAll(chain.getAtomGroups(GroupType.NUCLEOTIDE));
                    for (Group group : sourceGroups) {
                        sourceAtoms.addAll(group.getAtoms());
                    }
                }
                sourceArray = sourceAtoms.toArray(new Atom[sourceAtoms.size()]);
                break;

            case 3:
                matchAtoms = new ArrayList<>();
                matchChains = matchStructure.getChains();
                for (Chain chain : matchChains) {
                    List<Group> matchGroups = chain.getAtomGroups();
                    for (Group group : matchGroups) {
                        if (!group.isWater()) {
                            matchAtoms.addAll(group.getAtoms());
                        }
                    }
                }
                matchArray = matchAtoms.toArray(new Atom[matchAtoms.size()]);

                sourceAtoms = new ArrayList<>();
                sourceChains = sourceStructure.getChains();
                for (Chain chain : sourceChains) {
                    List<Group> matchGroups = chain.getAtomGroups();
                    for (Group group : matchGroups) {
                        if (!group.isWater()) {
                            sourceAtoms.addAll(group.getAtoms());
                        }
                    }
                }
                sourceArray = sourceAtoms.toArray(new Atom[sourceAtoms.size()]);
                break;

            case 4:
                matchArray = StructureTools.getAllAtomArray(matchStructure);
                sourceArray = StructureTools.getAllAtomArray(sourceStructure);
                break;
            default:
                throw new IllegalArgumentException("atomsUsed is not 1-4: this has not been properly checked at an earlier stage.");
        } // Returns simpleRMSD(matchArray, sourceArray) below old code.

        /*List<Atom> matchAtoms = new ArrayList<>();
         Structure matchStructure = matchFile.getStructure();
         List<Chain> matchChains = matchStructure.getChains();
         if (atomsUsed == 4) {
         matchAtoms = StructureTools.getAllAtomArray(matchStructure);
         }
         for (Chain chain : matchChains) {
         List<Group> matchGroups = chain.getSeqResGroups();
         for (Group group : matchGroups) {
         String groupType = group.getType();
         if (groupType.equalsIgnoreCase("HETATOM")) {
         if (atomsUsed < 3) {
         continue;
         } else if (atomsUsed < 4 && group.isWater()) {
         continue;
         }
         }
         if (atomsUsed != 1) {
         matchAtoms.addAll(group.getAtoms());
         } else {
         try {
         matchAtoms.add(getReferenceAtom(group));
         } catch (StructureException ex) {
         String refAtomType = (groupType.equalsIgnoreCase("AminoAcid") ? "CA" : "N1/9");
         // refAtomType should be binary between amino acid and nucleic acid, as HETATOM is eliminated.
         logger.info(String.format(" Reference atom %s could not be found for group %s",
         refAtomType, group.toString()));
         }
         }
         }
         }
         Atom[] match = new Atom[matchAtoms.size()];
         matchAtoms.toArray(match);
         matchAtoms.clear();

         List<Atom> sourceAtoms = new ArrayList<>();
         List<Chain> sourceChains = sourceStructure.getChains();
         for (Chain chain : sourceChains) {
         List<Group> sourceGroups = chain.getSeqResGroups();
         for (Group group : sourceGroups) {
         String groupType = group.getType();
         if (groupType.equalsIgnoreCase("HETATOM")) {
         if (atomsUsed < 3) {
         continue;
         } else if (atomsUsed < 4 && group.isWater()) {
         continue;
         }
         }
         if (atomsUsed != 1) {
         sourceAtoms.addAll(group.getAtoms());
         } else {
         try {
         sourceAtoms.add(getReferenceAtom(group));
         } catch (StructureException ex) {
         logger.info(String.format(" Reference atom could not be found for group %s", group.toString()));
         }
         }
         }
         }
         Atom[] sourceArray = new Atom[sourceAtoms.size()];
         sourceAtoms.toArray(sourceArray);
         sourceAtoms.clear();*/
        return simpleRMSD(sourceArray, matchArray);
    }

    private double simpleRMSD(Atom[] atomsA, Atom[] atomsB) throws IllegalArgumentException {
        double rmsd = 0.0;
        //double comp = 0.0; PDB file precision more important than addition error.
        int numMatches = 0;
        int numAtomsA = atomsA.length;
        for (int i = 0; i < numAtomsA; i++) {
            Atom atomFromA = atomsA[i];
            try {
                Atom atomFromB = getMatchingAtom(atomFromA, atomsB, i);
                ++numMatches; // Done here in case of exceptions thrown.
                rmsd += getSqDistance(atomFromA, atomFromB);
                /*double dist2 = getSqDistance(atomFromA, atomFromB) - comp; // Kahan summation algorithm.
                 double temp = rmsd + dist2;
                 comp = (temp - rmsd) - dist2;
                 rmsd = temp;*/
            } catch (IllegalArgumentException ex) {
                if (!atomFromA.getElement().isHydrogen()) {
                    // Will eventually move logger statement here.
                }
                logger.info(String.format(" Error in finding mate for atom %s", atomFromA.toString()));
            }
        }
        if (numMatches == 0) {
            throw new IllegalArgumentException(" No atomic matches found when calculating RMSD.");
        }
        rmsd = Math.sqrt(rmsd / numMatches);
        return rmsd;
    }

    private double getSqDistance(Atom a1, Atom a2) {
        double dx = a1.getX();
        dx -= a2.getX();
        double dy = a1.getY();
        dy -= a2.getY();
        double dz = a2.getZ();
        dz -= a2.getZ();
        return dx * dx + dy * dy + dz * dz;
    }

    private double getDistance(Atom a1, Atom a2) {
        return Math.sqrt(getSqDistance(a1, a2));
    }

    /**
     * Finds atom1's match in atoms2[], given an array index to search first.
     *
     * @param atom1 An Atom
     * @param atoms2 An Atom[] to search
     * @param i Index in atoms2 to check first.
     * @return atom1's match in atoms2.
     * @throws IllegalArgumentException If no match can be found.
     */
    private Atom getMatchingAtom(Atom atom1, Atom[] atoms2, int i) throws IllegalArgumentException {
        Atom atom2;
        try {
            atom2 = atoms2[i];
        } catch (ArrayIndexOutOfBoundsException ex) { // May be important if one structure lacks hydrogens.
            atom2 = atoms2[0];
        }
        if (compareAtoms(atom1, atom2)) {
            return atom2;
        }
        atom2 = getMatchingAtom(atom1, atoms2);
        return atom2;
    }

    /**
     * Finds atom1's match in structure2; uses atom1's residue number, atom
     * type, and PDB serial number as a guess; if robustMatch is set true, will
     * search all Atoms in structure2.
     *
     * @param atom1 An Atom
     * @param structure2 A Structure to search
     * @param searchAll Whether to search all Atoms in structure2.
     * @return atom1's match in structure2
     * @throws IllegalArgumentException If no match can be found.
     */
    private Atom getMatchingAtom(Atom atom1, Structure structure2, boolean searchAll) throws IllegalArgumentException {
        ResidueNumber res1 = atom1.getGroup().getResidueNumber();
        String chainID = res1.getChainId();
        Atom atom2 = null;
        try {
            Chain chain2 = structure2.getChainByPDB(chainID);
            Group group2 = chain2.getGroupByPDB(res1);
            atom2 = group2.getAtom(atom1.getName());
            if (atom1.getName().equalsIgnoreCase("H")) {
                atom2 = group2.getAtom("H1");
            } else if (atom1.getName().equalsIgnoreCase("H1")) {
                atom2 = group2.getAtom("H");
                }
                atom2 = group2.getAtom(atom1.getPDBserial());
        } catch (StructureException ex) {
            if (!searchAll) {
                throw new IllegalArgumentException("Matching atom not found.");
            }
            for (Chain chain : structure2.getChains()) {
                for (Group group : chain.getAtomGroups()) {
                    for (Atom atom : group.getAtoms()) {
                        if (compareAtoms(atom1, atom)) {
                            return atom2;
                        }
                    }
                }
            }
        }
        if (atom2 != null && compareAtoms(atom1, atom2)) {
            return atom2;
        }
        throw new IllegalArgumentException("Matching atom not found.");
    }

    /**
     * Searches for atom1's match in atoms2; first checks for an equivalent
     * atom, then performs a search over all atoms in atoms2.
     *
     * @param atom1 An Atom.
     * @param atoms2 An Atom[] to search.
     * @return atom1's match in atom2.
     * @throws IllegalArgumentException If no match could be found.
     */
    private Atom getMatchingAtom(Atom atom1, Atom[] atoms2) throws IllegalArgumentException {
        Atom atom2 = atoms2[0];
        Structure structure2 = atom2.getGroup().getChain().getStructure();
        try {
            atom2 = getMatchingAtom(atom1, structure2, false);
            return atom2;
        } catch (IllegalArgumentException ex) {
            for (Atom atom : atoms2) {
                if (compareAtoms(atom1, atom)) {
                    return atom;
                }
            }
        }
        throw new IllegalArgumentException(String.format("No matching atom for %s found", atom1.toString()));
    }

    /**
     * Mildly robust atom comparison tool which attempts to account for things
     * like missing hydrogens causing inequal PDB serial numbering. Compares on
     * element, group PDB code, atom serial (if false, also checks group
     * number), and atom name (currently only accepts H-H1 mismatches).
     *
     * @param a1
     * @param a2
     * @return if considered equivalent.
     */
    private boolean compareAtoms(Atom a1, Atom a2) {
        if (!a1.getElement().equals(a2.getElement())) {
            return false;
        }
        Group group1 = a1.getGroup();
        Group group2 = a2.getGroup();
        if (!group1.getPDBName().equalsIgnoreCase(group2.getPDBName())) {
            return false;
        }
        if (a1.getPDBserial() != a2.getPDBserial()) {
            // Second check accounts for situations where one structure may be H-trimmed.
            if (!group1.getResidueNumber().equals(group2.getResidueNumber())) {
                return false;
            }
        }
        String a1Name = a1.getName();
        String a2Name = a2.getName();
        if (!a1Name.equalsIgnoreCase(a2Name)) {
            // Courtesy of MSMBuilder, I have to keep track of Hs which should be H1.
            if (!((a1Name.equalsIgnoreCase("H") && a2Name.equalsIgnoreCase("H1"))
                    || (a2Name.equalsIgnoreCase("H") && a1Name.equalsIgnoreCase("H1")))) {
                return false;
            }
        }
        return true;
    }

    /**
     * May rework once I figure out whether BioJava likes to keep atom names
     * intact or not: returns CA (proteins), N1/9 (nucleic acids), or the first
     * returned atom (default) as in FFX's Residue.getReferenceAtom().
     *
     * @param group Residue to get reference atom for.
     * @return CA or N1/9 Atom.
     * @throws StructureException If no proper reference Atom could be found.
     */
    private Atom getReferenceAtom(Group group) throws StructureException {
        switch (group.getType()) {
            case AMINOACID:
                return group.getAtom("CA");
            case NUCLEOTIDE:
                Atom retAtom = group.getAtom("N1");
                if (retAtom == null) {
                    retAtom = group.getAtom("N9");
                }
                return retAtom;
            default:
                return group.getAtoms().get(0);
        }
    }

    private File createVersionedCopy(File origFile) throws IOException {
        String filename = origFile.getName();
        String ext = FilenameUtils.getExtension(filename);
        filename = FilenameUtils.removeExtension(origFile.getName()).concat(fileSuffix).concat(ext);
        File retFile = new File(filename);
        if (retFile.exists()) {
            for (int i = 2; i < 1000; i++) {
                retFile = new File(filename.concat("_" + i));
                if (!retFile.exists()) {
                    break;
                }
            }
            if (retFile.exists()) {
                String exSt = "Versioning failed: all filenames from " + filename
                        + " to " + filename + "_999 already exist.";
                throw new IOException(exSt);
            }
        }
        return retFile;
    }

    /**
     * Compares two SSBonds based on chain IDs, insertion codes, and residue
     * numbers; not order-sensitive.
     *
     * @param bond1
     * @param bond2
     * @return If bond1 and bond2 are equivalent.
     */
    private boolean compareSSBonds(SSBond bond1, SSBond bond2) {
        List<String> bond1ChainIDs = new ArrayList<>(2);
        bond1ChainIDs.add(bond1.getChainID1());
        bond1ChainIDs.add(bond1.getChainID2());
        List<String> bond2ChainIDs = new ArrayList<>(2);
        bond2ChainIDs.add(bond2.getChainID1());
        bond2ChainIDs.add(bond2.getChainID2());
        for (String resID : bond1ChainIDs) {
            if (!bond2ChainIDs.contains(resID)) {
                return false;
            }
        }

        List<String> bond1InsCodes = new ArrayList<>(2);
        bond1InsCodes.add(bond1.getInsCode1());
        bond1InsCodes.add(bond1.getInsCode2());
        List<String> bond2InsCodes = new ArrayList<>(2);
        bond2InsCodes.add(bond2.getInsCode1());
        bond2InsCodes.add(bond2.getInsCode2());
        for (String resID : bond1InsCodes) {
            if (!bond2InsCodes.contains(resID)) {
                return false;
            }
        }

        List<String> bond1ResNums = new ArrayList<>(2);
        bond1ResNums.add(bond1.getResnum1());
        bond1ResNums.add(bond1.getResnum2());
        List<String> bond2ResNums = new ArrayList<>(2);
        bond2ResNums.add(bond2.getResnum1());
        bond2ResNums.add(bond2.getResnum2());
        for (String resID : bond1ResNums) {
            if (!bond2ResNums.contains(resID)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Creates a duplicate PDBCrystallographicInfo object.
     *
     * @param sourceInfo
     * @return duplicate of sourceInfo
     * @throws IllegalArgumentException if sourceInfo not crystallographic.
     */
    private PDBCrystallographicInfo cloneCrystalInfo(PDBCrystallographicInfo sourceInfo)
            throws IllegalArgumentException {

        PDBCrystallographicInfo retInfo = new PDBCrystallographicInfo();
        retInfo.setSpaceGroup(sourceInfo.getSpaceGroup());

        // Newer versions of BioJava can use the CrystalCell object.
        CrystalCell cell;
        cell = new CrystalCell();
        
        cell.setA(sourceInfo.getA());
        cell.setAlpha(sourceInfo.getAlpha());
        cell.setB(sourceInfo.getB());
        cell.setBeta(sourceInfo.getBeta());
        cell.setC(sourceInfo.getC());
        cell.setGamma(sourceInfo.getGamma());
        
        retInfo.setCrystalCell(cell);

        return retInfo;
    }

    /**
     * Interior code of file matching loop: intended to ensure consistency
     * between sequential and parallel versions of code.
     *
     * @param currentFilePair
     * @param filereader
     * @param aligner
     * @throws IOException
     */
    private void matchFile(FileFilePair currentFilePair, PDBFileReader filereader,
            StructurePairAligner aligner) throws IOException {
        Structure currentStructure = filereader.getStructure(currentFilePair.getMatchedFile());
        currentFilePair.setStructure(currentStructure);
        int numSources = sourceFiles.length;
        for (int j = 0; j < numSources; j++) {
            File currentSource = sourceFiles[j];
            Structure sourceStructure = filereader.getStructure(currentSource);
            try {
                double rmsd = calculateRMSD(currentFilePair, sourceStructure, aligner);
                currentFilePair.attemptReplace(currentSource, rmsd);
            } catch (StructureException ex) {
                logger.warning(String.format(" Error in calculating RMSD for match %s and source %s : %s",
                        currentFilePair.getMatchedFile().getName(), currentSource.getName(), ex.toString()));
            }
        }
    }

    private void fixFile(FileFilePair currentPair, PDBFileReader filereader) throws IOException {
        File matchFile = currentPair.getMatchedFile();
        Structure matchStructure = currentPair.getStructure();
        if (matchStructure == null) {
            matchStructure = filereader.getStructure(matchFile);
        }

        File sourceFile = currentPair.getSourceFile();
        if (sourceFile == null) {
            throw new IOException(String.format("No source file was matched to file %s", matchFile.toString()));
        }

        Structure sourceStructure = null;
        if (fixAtoms) {
            sourceStructure = filereader.getStructure(sourceFile);
            Atom[] matchAtoms = StructureTools.getAllAtomArray(matchStructure);
            for (Atom matchAtom : matchAtoms) {
                Atom sourceAtom = getMatchingAtom(matchAtom, sourceStructure, robustMatch);
                if (fixBFactors) {
                    matchAtom.setTempFactor(sourceAtom.getTempFactor());
                }
            }
            // Other methods can go here.
        }
        if (fixSSBonds) {
            if (sourceStructure == null) {
                sourceStructure = filereader.getStructure(sourceFile);
            }
            List<SSBond> sourceBonds = sourceStructure.getSSBonds();
            List<SSBond> matchBonds = matchStructure.getSSBonds();
            for (SSBond sourceBond : sourceBonds) {
                boolean isContained = false;
                for (SSBond matchBond : matchBonds) {
                    if (compareSSBonds(matchBond, sourceBond)) {
                        isContained = true;
                        break;
                    }
                }
                if (!isContained) {
                    matchStructure.addSSBond(sourceBond.clone());
                }
            }
        }
        if (fixCryst) {
            if (sourceStructure == null) {
                sourceStructure = filereader.getStructure(sourceFile);
            }
            PDBCrystallographicInfo crystalInfo = sourceStructure.getCrystallographicInfo();
            try {
                PDBCrystallographicInfo duplicateInfo = cloneCrystalInfo(crystalInfo);
                matchStructure.setCrystallographicInfo(duplicateInfo);
            } catch (IllegalArgumentException ex) {
                logger.warning(String.format(" No crystal information for source structure "
                        + "%s: nothing attached to file %s", sourceFile.toString(), matchFile.toString()));
            }
        }
        String pdb = matchStructure.toPDB();
        if (headerLink) {
            StringBuilder pdbBuilder = new StringBuilder(pdb);
            int position = pdbBuilder.lastIndexOf("REMARK ");
            int remarkNumber = 4;
            if (position >= 0) {
                String nextLine = pdbBuilder.substring(position, position + 1000);
                int offset = nextLine.indexOf("%n");
                if (offset < 0) {
                    nextLine = pdbBuilder.substring(position);
                    offset = nextLine.indexOf("%n");
                }
                position += offset;
                String[] tok = nextLine.split(" +", 3);
                try {
                    remarkNumber = Integer.parseInt(tok[1]) + 1;
                } catch (NumberFormatException ex) {
                    // Silent.
                }
            }

            String toInsert = String.format("REMARK%4d SOURCE FILE: %s", remarkNumber, sourceFile.getName());
            toInsert = toInsert.concat(String.format("REMARK%4d RMSD:%11.6f ANGSTROMS", remarkNumber, currentPair.getRMSD()));
            pdbBuilder.insert(position, toInsert);
            pdb = pdbBuilder.toString();
        }

        File newFile = createVersionedCopy(matchFile);
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(newFile))) {
            try {
                bw.write(pdb);
            } catch (IOException ex) {
                logger.warning(String.format(" Error writing to file %s", newFile.getName()));
            }
        }
    }

    private void sequentialFileMatch() {
        PDBFileReader filereader = new PDBFileReader();
        StructurePairAligner aligner = new StructurePairAligner();
        //int numSources = sourceFiles.length;
        try {
            for (int i = 0; i < matchFilePairs.length; i++) {
                Long iterTime = -System.nanoTime();
                FileFilePair currentFilePair = matchFilePairs[i];
                matchFile(currentFilePair, filereader, aligner);
                iterationTimes[i] = iterTime + System.nanoTime();
            }
        } catch (IOException ex) {
            logger.severe(String.format(" Exception matching files %s", ex.toString()));
        }
    }

    private void sequentialFixer() {
        PDBFileReader filereader = new PDBFileReader();
        for (int i = 0; i < matchFilePairs.length; i++) {
            Long fixTime = -System.nanoTime();
            FileFilePair currentPair = matchFilePairs[i];
            try {
                fixFile(currentPair, filereader);
            } catch (IOException ex) {

            }
            fixTime += System.nanoTime();
            fixTimes[i] += fixTime;
        }
    }

    private class MatchingRegion extends ParallelRegion {

        @Override
        public void run() {
            try {
                execute(0, matchFilePairs.length, new MatchingLoop());
            } catch (Exception e) {
                String message = " Exception matching in parallel.";
                logger.log(Level.SEVERE, message, e);
            }
        }
    }

    private class MatchingLoop extends IntegerForLoop {

        @Override
        public void run(int lb, int ub) throws IOException {
            PDBFileReader filereader = new PDBFileReader();
            StructurePairAligner aligner = new StructurePairAligner();
            for (int i = lb; i <= ub; i++) {
                Long iterTime = -System.nanoTime();
                FileFilePair currentFilePair = matchFilePairs[i];
                matchFile(currentFilePair, filereader, aligner);
                iterTime += System.nanoTime();
                iterationTimes[i] = iterTime;
            }
        }
    }

    private class FixerRegion extends ParallelRegion {

        @Override
        public void run() {
            try {
                execute(0, matchFilePairs.length, new FixerLoop());
            } catch (Exception e) {
                String message = " Exception fixing files in parallel.";
                logger.log(Level.SEVERE, message, e);
            }
        }
    }

    private class FixerLoop extends IntegerForLoop {

        @Override
        public void run(int lb, int ub) throws IOException {
            PDBFileReader filereader = new PDBFileReader();
            for (int i = lb; i <= ub; i++) {
                Long fixTime = -System.nanoTime();
                FileFilePair currentPair = matchFilePairs[i];
                fixFile(currentPair, filereader);
                fixTime += System.nanoTime();
                fixTimes[i] = fixTime;
            }
        }
    }

    /**
     * Used to keep track of what is currently considered the best source file
     * candidate for an alignment file and its RMSD.
     */
    private class FileFilePair {

        private final File matchFile;
        private Structure matchStructure;
        private File sourceFile;
        private double rmsd;

        public FileFilePair(File matchFile) {
            this.matchFile = matchFile;
            rmsd = Double.MAX_VALUE;
        }

        public FileFilePair(File matchFile, File sourceFile) {
            this(matchFile);
            this.sourceFile = sourceFile;
        }

        public FileFilePair(File match, File sourceFile, double rmsd) {
            this(match, sourceFile);
            this.rmsd = rmsd;
        }

        public void setSourceFile(File sourceFile) {
            this.sourceFile = sourceFile;
        }

        public void setRMSD(double rmsd) {
            this.rmsd = rmsd;
        }

        public void setStructure(Structure structure) {
            this.matchStructure = structure;
        }

        public File getSourceFile() {
            return sourceFile;
        }

        public File getMatchedFile() {
            return matchFile;
        }

        public double getRMSD() {
            return rmsd;
        }

        public Structure getStructure() {
            return matchStructure;
        }

        /**
         * Will replace the source file if its RMSD is lower than current RMSD
         * and returns true; otherwise returns false.
         *
         * @param newSource Source file to check
         * @param rmsd A double
         * @return Whether replaced
         */
        public boolean attemptReplace(File newSource, double rmsd) {
            if (rmsd < this.rmsd) {
                this.sourceFile = newSource;
                this.rmsd = rmsd;
                return true;
            }
            return false;
        }

        /**
         * {@inheritDoc}
         *
         * Overidden equals method that return true if object is equals to this,
         * or is of the same class and has the same alignedFile and sourceFile.
         */
        @Override
        public boolean equals(Object object) {
            if (this == object) {
                return true;
            } else if (object == null || getClass() != object.getClass()) {
                return false;
            }
            FileFilePair other = (FileFilePair) object;
            return other.getMatchedFile().equals(this.matchFile) && other.getSourceFile().equals(this.sourceFile);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder("Aligned file: ").append(matchFile.getPath()).append("\nSource file: ");
            if (sourceFile != null) {
                sb.append(sourceFile.getPath());
            } else {
                sb.append("NULL");
            }
            if (rmsd != Double.MAX_VALUE) {
                sb.append(String.format("\nRMSD: %f", rmsd));
            } else {
                sb.append("\nRMSD: NAN (not set)");
            }
            return sb.toString();
        }
    }

    /*private class AtomAtomPair {
     private final Atom atom1;
     private Atom atom2;

     public AtomAtomPair(Atom atom1) {
     this.atom1 = atom1;
     }
     public AtomAtomPair(Atom atom1, Atom atom2) {
     this(atom1);
     this.atom2 = atom2;
     }

     public Atom getAtom1() {
     return atom1;
     }
     public Atom getAtom2() {
     return atom2;
     }
     public Atom[] getAtoms() {
     Atom[] ret = {atom1, atom2};
     return ret;
     }
     public void setAtom2(Atom a2) {
     this.atom2 = a2;
     }
     /**
     * Returns true and sets atom2 if a2 matches (using compareAtoms) atom1.
     * @param a2 Candidate Atom match for atom1
     * @return If a match
     */
    /*public boolean tryMatch(Atom a2) {
     if (compareAtoms(atom1, a2)) {
     this.atom2 = a2;
     return true;
     } else {
     return false;
     }
     }*/
    /**
     * {@inheritDoc}
     *
     * Overidden equals method that return true if object is equals to this, or
     * is of the same class and has the same alignedFile and sourceFile.
     */
    /*@Override
     public boolean equals(Object object) {
     if (this == object) {
     return true;
     } else if (object == null || getClass() != object.getClass()) {
     return false;
     }
     AtomAtomPair other = (AtomAtomPair) object;
     return atom1.equals(other.getAtom1()) && atom2.equals(other.getAtom2());
     }

     @Override
     public String toString() {
     return new StringBuilder("Atom 1: ").append(atom1.toString()).append("   Atom 2: ").append(atom2.toString()).toString();
     }
     }*/
}
