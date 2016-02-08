/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.potential;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.logging.Logger;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;

import static ffx.numerics.VectorMath.diff;

/**
 * The Utilities class provides methods to locate functional units of an organic
 * system.
 *
 * @author Michael J. Schnieders
 *
 */
public final class Utilities {

    private static final Logger logger = Logger.getLogger(Utilities.class.getName());

    /**
     * An enumeration of recognized file types.
     */
    public enum FileType {

        XYZ, INT, ARC, PDB, ANY, SIM, UNK
    }
    
    public enum DataType {
        BIOJAVA, UNK
    }

    /**
     * An enumeration of recognized organic polymers.
     */
    public enum PolymerType {

        AMINOACID, NUCLEICACID, UNKNOWN
    }
    /**
     * The algorithms used to split arrays of atoms into a structural hierarchy
     * use a bunch of List instances. Pooling them seems like it might be a
     * performance win - although better algorithms probably exist. This is
     * currently backed by ArrayLists.
     */
    public static List<List<Atom>> atomListPool = new ArrayList<>();
    static int count = 0;
    /**
     * Repeating atomic numbers of an amino acid chain.
     */
    public static final int AAPATTERN[] = {7, 6, 6};
    /**
     * Repeating atomic numbers of a nucleic acid chain.
     */
    public static final int NAPATTERN[] = {8, 6, 6, 6, 8};
    /**
     * Stoichiometry of side chains can be used for identification, accept for a
     * couple cases: 1.) Proline & Valine 2.) (Iso)Leucine 3.) DNA Gaunine/RNA
     * Adenine. This Hashtable returns the 3-letter name for amino acids, a
     * single character for nucleic acids, or an integer indicating a special
     * case.
     */
    private static final HashMap<String, String> sidechainStoichiometry = new HashMap<>();

    private static final double p4 = 15.236;
    private static final double p5 = 1.254;
    private static final double p5inv = 1.0 / 1.254;
    private static final double pip5 = Math.PI * p5;
    private static final double convert = -332.05382 / 2.0;
    private static final double[] x1 = new double[3];
    private static final double[] x2 = new double[3];

    static {
        // Amino Acid Side Chains
        sidechainStoichiometry.put("S1C3", "MET");
        sidechainStoichiometry.put("S1C1", "CYS");
        sidechainStoichiometry.put("O1C1", "SER");
        sidechainStoichiometry.put("O1C2", "THR");
        sidechainStoichiometry.put("O1C7", "TYR");
        sidechainStoichiometry.put("O2C2", "ASP");
        sidechainStoichiometry.put("O2C3", "GLU");
        sidechainStoichiometry.put("O1N1C2", "ASN");
        sidechainStoichiometry.put("O1N1C3", "GLN");
        sidechainStoichiometry.put("N3C4", "ARG");
        sidechainStoichiometry.put("N2C4", "HIS");
        sidechainStoichiometry.put("N1C9", "TRP");
        sidechainStoichiometry.put("N1C4", "LYS");
        sidechainStoichiometry.put("C7", "PHE");
        sidechainStoichiometry.put("H", "GLY");
        sidechainStoichiometry.put("C1", "ALA");
        // DNA
        sidechainStoichiometry.put("O2N3C6", "DC");
        sidechainStoichiometry.put("O1N5C7", "DA");
        sidechainStoichiometry.put("O3N2C7", "DT");
        // RNA
        sidechainStoichiometry.put("O3N5C7", "G");
        sidechainStoichiometry.put("O3N3C6", "C");
        sidechainStoichiometry.put("O4N2C6", "U");
        // SPECIAL CASES
        sidechainStoichiometry.put("C3", "1"); // Proline / Valine
        sidechainStoichiometry.put("C4", "2"); // (ISO)Leucine
        sidechainStoichiometry.put("O2N5C7", "3"); // DNA Gaunine / RNA Adenine
    }

    /**
     * <p>
     * addAtomListToPool</p>
     *
     * @param a a {@link java.util.List} object.
     */
    public static void addAtomListToPool(List<Atom> a) {
        a.clear();
        atomListPool.add(a);
    }

    /**
     * Collect all the atoms in a polymer opposite the end atom, and put them
     * into the residue. This assumes that the "cap" is only linked to the rest
     * of the polymer through the end atom.
     *
     * @param end Atom
     * @param seed Atom
     * @param residue Residue
     */
    public static void addCap(Atom end, Atom seed, Residue residue) {
        List<Atom> cap = getAtomListFromPool();
        cap.add(end);
        collectAtoms(seed, cap);
        // Assume the end & seed atoms are already part of the residue
        cap.remove(0);
        cap.remove(0);
        for (Atom a : cap) {
            residue.addMSNode(a);
        }
    }

    /**
     * Add a phosphate and its bonded oxygens that are not bonded to a carbon to
     * the specified residue.
     *
     * @param phosphate Atom
     * @param residue Residue
     */
    public static void addPhosphate(Atom phosphate, Residue residue) {
        if (phosphate == null) {
            return;
        }
        residue.addMSNode(phosphate);
        for (Bond b : phosphate.getBonds()) {
            Atom oxygen = b.get1_2(phosphate);
            // Add oxygens not bonded to a Carbon
            if (numberOfBondsWith(oxygen, 6) == 0) {
                residue.addMSNode(oxygen);
                // Add hydrogens atoms for protonated oxygen groups
                Atom hydrogen = findBondWith(oxygen, 1);
                if (hydrogen != null) {
                    residue.addMSNode(hydrogen);
                }
            }
        }
    }

    private static Residue assignResidue(List<Atom> backbone, int start,
            List<Atom> atoms, List<Atom> sidePolymer) {
        Atom a;
        int atomicnum;
        int bins[] = new int[5]; // 0 = S, 1 = P, 2 = O, 3 = N, 4 = C
        char chars[] = {'S', 'P', 'O', 'N', 'C'};
        for (ListIterator li = sidePolymer.listIterator(); li.hasNext();) {
            a = (Atom) li.next();
            atomicnum = a.getAtomicNumber();
            switch (atomicnum) {
                case 1:
                    // ignore hydrogens
                    break;
                case 6:
                    // Carbon
                    bins[4]++;
                    break;
                case 7:
                    // Nitrogen
                    bins[3]++;
                    break;
                case 8:
                    // Oxygen
                    bins[2]++;
                    break;
                case 15:
                    // Phosphorus
                    bins[1]++;
                    break;
                case 16:
                    // Sulfur
                    bins[0]++;
                    break;
                default:
                    return null;
            }
        }
        StringBuilder key = new StringBuilder();
        int atomCount = 0;
        for (int i = 0; i < 5; i++) {
            if (bins[i] != 0) {
                atomCount += bins[i];
                key.append(chars[i]);
                key.append(Integer.toString(bins[i]));
            }
        }
        if (atomCount == 0) {
            key.append("H"); // Glycine
        }
        String resname = sidechainStoichiometry.get(key.toString());
        if (resname == null) {
            resname = "Unknown";
        } else {
            resname = resname.intern();
        }
        if (resname.equals("1") || resname.equals("2")) {
            // Special case where atom string keys aren't unique
            Atom alpha = backbone.get(start + 1);
            Atom carbonyl = backbone.get(start + 2);
            Atom beta = null;
            List alphabonds = alpha.getBonds();
            Bond abond;
            for (ListIterator li = alphabonds.listIterator(); li.hasNext();) {
                abond = (Bond) li.next();
                beta = abond.get1_2(alpha);
                // Don't want the peptide nitrogen or alpha hydrogen or carbonyl
                // carbon
                if (beta.getAtomicNumber() != 7 && beta.getAtomicNumber() != 1 && beta != carbonyl) {
                    break;
                }
                beta = null;
            }
            if (beta == null) {
                return null;
            }
            List betabonds = beta.getBonds();
            Atom gamma;
            int carboncount = 0;
            for (ListIterator li = betabonds.listIterator(); li.hasNext();) {
                abond = (Bond) li.next();
                gamma = abond.get1_2(beta);
                if (gamma.getAtomicNumber() == 6) {
                    carboncount++;
                }
            }
            if (resname.equals("1")) {
                if (carboncount == 2) {
                    resname = "PRO";
                } else {
                    resname = "VAL";
                }
            } else {
                if (carboncount == 2) {
                    resname = "LEU";
                } else {
                    resname = "ILE";
                }
            }
        } else if (resname.equals("3")) {
            Atom c3 = backbone.get(start + 3);
            int num = countCO(c3);
            if (num == 2) {
                resname = "A";
            } else {
                resname = "DG";
            }
        }
        Residue residue = null;
        try {
            Residue.NA3.valueOf(resname.toUpperCase());
            residue = new Residue(resname, Residue.ResidueType.NA);
        } catch (Exception e) {
        }
        if (residue == null) {
            try {
                Residue.AA3.valueOf(resname.toUpperCase());
                residue = new Residue(resname, Residue.ResidueType.AA);
            } catch (Exception e) {
            }
        }
        if (residue == null) {
            residue = new Residue(resname, Residue.ResidueType.UNK);
        }
        // Create the Residue group
        for (ListIterator li = atoms.listIterator(); li.hasNext();) {
            a = (Atom) li.next();
            residue.addMSNode(a);
        }
        return residue;
    }

    /**
     * This routine sub-divides a system into groups of ions, water, hetero
     * molecules, and polynucleotides/polypeptides.
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param atoms a {@link java.util.List} object.
     */
    public static void biochemistry(MolecularAssembly molecularAssembly, List<Atom> atoms) {
        Atom atom, seed = null;
        int num = 0;
        int waterNum = 0;
        int ionNum = 0;
        int moleculeNum = 0;
        List<String> segIDs = new ArrayList<>();
        while (atoms.size() > 0) {
            /**
             * Nitrogen is used to "seed" a backbone search because carbon can
             * be separated from the backbone by a sulfur (ie. MET).
             */
            for (Atom a : atoms) {
                seed = a;
                if (seed.getAtomicNumber() == 7) {
                    break;
                }
            }
            /**
             * If no nitrogen atoms remain, there are no nucleic or amino acids.
             */
            if (seed.getAtomicNumber() != 7) {
                List<Atom> moleculeAtoms;
                while (atoms.size() > 0) {
                    atom = atoms.get(0);
                    // Check for a metal ion or noble gas
                    if (atom.getNumBonds() == 0) {
                        ionNum++;
                        Molecule ion = new Molecule(atom.getName() + "-" + ionNum);
                        ion.addMSNode(atom);
                        atoms.remove(0);
                        molecularAssembly.addMSNode(ion);
                        continue;
                    } // Check for water
                    else if (atom.getAtomicNumber() == 8 && isWaterOxygen(atom)) {
                        waterNum++;
                        Molecule water = new Molecule("Water-" + waterNum);
                        water.addMSNode(atom);
                        atoms.remove(0);
                        List<Bond> bonds = atom.getBonds();
                        for (Bond b : bonds) {
                            Atom hydrogen = b.get1_2(atom);
                            water.addMSNode(hydrogen);
                            atoms.remove(hydrogen);
                        }
                        molecularAssembly.addMSNode(water);
                        continue;
                    }
                    // Otherwise classify the molecule as a hetero
                    moleculeNum++;
                    Molecule molecule = new Molecule("Molecule-" + moleculeNum);
                    moleculeAtoms = getAtomListFromPool();
                    collectAtoms(atoms.get(0), moleculeAtoms);
                    while (moleculeAtoms.size() > 0) {
                        atom = moleculeAtoms.get(0);
                        moleculeAtoms.remove(0);
                        molecule.addMSNode(atom);
                        atoms.remove(atom);
                    }
                    molecularAssembly.addMSNode(molecule);
                }
                seed = null;
                break;
            }
            List<Atom> backbone = findPolymer(atoms, seed, null);
            if (backbone.size() > 0) {
                for (ListIterator li = backbone.listIterator(backbone.size()); li.hasPrevious();) {
                    seed = (Atom) li.previous();
                    if (seed.getAtomicNumber() == 7) {
                        break;
                    }
                }
                backbone = findPolymer(atoms, seed, null);
            }
            Character chainID = getChainID(num);
            String segID = getSegID(chainID, segIDs);
            Polymer c = new Polymer(chainID, segID, true);
            if (backbone.size() > 2 && divideBackbone(backbone, c)) {
                for (Atom a : c.getAtomList()) {
                    atoms.remove(a);
                }
                logger.fine(" Sequenced chain: " + c.getName());
                molecularAssembly.addMSNode(c);
                num++;
            } else {
                moleculeNum++;
                Molecule hetero = new Molecule("" + moleculeNum + "-Hetero");
                atom = backbone.get(0);
                List<Atom> heteroAtomList = getAtomListFromPool();
                collectAtoms(atom, heteroAtomList);
                for (Atom a : heteroAtomList) {
                    hetero.addMSNode(a);
                }
                for (Atom a : hetero.getAtomList()) {
                    atoms.remove(a);
                }
                molecularAssembly.addMSNode(hetero);
            }
        }
    }

    /**
     * Convert possibly duplicate chainID into a unique segID.
     *
     * @param c chain ID just read.
     * @return a unique segID.
     */
    private static String getSegID(Character c, List<String> segIDs) {
        if (c == null || c.equals(' ')) {
            c = 'A';
        }

        // Loop through existing segIDs to find the first one that is unused.
        int n = segIDs.size();
        int m = 0;
        for (int i = 0; i < n; i++) {
            String segID = segIDs.get(i);
            if (segID.endsWith(c.toString())) {
                m++;
            }
        }

        // If the count is greater than 0, then append it.
        String newSegID = null;
        if (m == 0) {
            newSegID = c.toString();
        } else {
            newSegID = c.toString() + Integer.toString(m);
        }

        segIDs.add(newSegID);
        return newSegID;
    }

    /**
     * Given an array of bonded atoms, this function recursively collects all
     * other connected atoms, without backtracking over atoms already in the
     * list. Disulfide bonds are not crossed. (the intent is to search along a
     * peptide backbone) Atoms preloaded into the List provide search
     * termination.
     *
     * @param seed Atom
     * @param atoms List
     */
    private static void collectAtoms(Atom seed, List<Atom> atoms) {
        if (seed == null) {
            return;
        }
        atoms.add(seed);
        for (Bond b : seed.getBonds()) {
            Atom nextAtom = b.get1_2(seed);
            if (nextAtom.getParent() != null) {
                continue;
            }
            // avoid crossing disulfides
            if ((nextAtom.getAtomicNumber() != 16 || seed.getAtomicNumber() != 16) && !atoms.contains(nextAtom)) {
                collectAtoms(nextAtom, atoms);
            }
        }
    }

    /**
     * <p>
     * countCO</p>
     *
     * @param adjacent a {@link ffx.potential.bonded.Atom} object.
     * @return a int.
     */
    public static int countCO(Atom adjacent) {
        int total = 0;
        for (Bond b : adjacent.getBonds()) {
            Atom carbonyl = b.get1_2(adjacent);
            if (carbonyl.getAtomicNumber() == 6) {
                for (Bond b2 : carbonyl.getBonds()) {
                    Atom oxygen = b2.get1_2(carbonyl);
                    if (oxygen.getAtomicNumber() == 8) {
                        total++;
                    }
                }
            }
        }
        return total;
    }

    /**
     * <p>
     * divideBackbone</p>
     *
     * @param backbone a {@link java.util.List} object.
     * @param c a {@link ffx.potential.bonded.Polymer} object.
     * @return a boolean.
     */
    public static boolean divideBackbone(List<Atom> backbone, Polymer c) {
        int length = backbone.size();
        // Try to find a Phosphorus or Nitrogen in the backbone
        int n, p;
        n = p = 0;
        for (Atom match : backbone) {
            int an = match.getAtomicNumber();
            if (an == 15) {
                p++;
            } else if (an == 7) {
                n++;
            }
        }
        PolymerType type;
        if (p >= n && p > 1) {
            type = PolymerType.NUCLEICACID;
        } else if (n > p && n > 2) {
            type = PolymerType.AMINOACID;
        } else {
            return false;
        }
        int start = -1;
        for (int i = 0; i < length; i++) {
            Residue res = patternMatch(i, backbone, type);
            if (res != null) {
                for (Atom a : res.getAtomList()) {
                    a.setParent(null);
                }
                if (!(res.getName().equals("Unknown"))) {
                    start = i;
                    // Want 5' to 3'
                    if (type == PolymerType.NUCLEICACID) {
                        Atom carbon5 = backbone.get(start + 1);
                        if (numberOfBondsWith(carbon5, 6) != 1) {
                            start = -1;
                        }
                    }
                    break;
                }
            }
        }
        if (start == -1) {
            backbone = reverseAtomList(backbone);
            for (int i = 0; i < length; i++) {
                Residue res = patternMatch(i, backbone, type);
                if (res != null) {
                    for (Atom a : res.getAtomList()) {
                        a.setParent(null);
                    }
                    if (!(res.getName().equals("Unknown"))) {
                        start = i;
                        break;
                    }
                }
            }
        }
        if (start == -1) {
            return false;
        }
        // Potential Polypeptide
        if (type == PolymerType.AMINOACID) {
            Atom nitrogen, alpha, carbonyl = null;
            Atom nitrogen2, carbonyl2;
            Residue aa = null;
            int lastRes = 0;
            int firstRes = -1;
            List<Residue> aaArray = new ArrayList<Residue>();
            while (start < length) {
                aa = patternMatch(start, backbone, PolymerType.AMINOACID);
                if (aa != null) {
                    if (firstRes == -1) {
                        firstRes = start;
                        carbonyl = findCarbonyl(backbone.get(start));
                    }
                    aaArray.add(aa);
                    lastRes = start;
                }
                start += 3;
            }
            // Make sure the fisrt residue is found
            aa = null;
            if (carbonyl != null) {
                alpha = findAlphaCarbon(carbonyl);
                if (alpha != null) {
                    nitrogen = findBondWith(alpha, 7);
                    if (nitrogen != null) {
                        nitrogen2 = findBondWith(carbonyl, 7);
                        List<Atom> firstAtoms = getAtomListFromPool();
                        firstAtoms.add(nitrogen);
                        firstAtoms.add(alpha);
                        firstAtoms.add(carbonyl);
                        firstAtoms.add(nitrogen2);
                        aa = patternMatch(0, firstAtoms, PolymerType.AMINOACID);
                        addAtomListToPool(firstAtoms);
                        if (aa != null) {
                            addCap(alpha, nitrogen, aa);
                            aaArray.add(0, aa);
                        }
                    }
                }
            }
            // Add the remaining atoms to the end of the Polymer
            if (aa == null) {
                nitrogen = backbone.get(firstRes);
                alpha = backbone.get(firstRes + 1);
                addCap(alpha, nitrogen, aaArray.get(0));
            }
            // Make sure the last residue is found
            aa = null;
            carbonyl = findCarbonyl(backbone.get(lastRes + 1));
            if (carbonyl != null) {
                nitrogen = findBondWith(carbonyl, 7);
                if (nitrogen != null) {
                    alpha = findAlphaCarbon(nitrogen);
                    if (alpha != null) {
                        carbonyl2 = findCarbonyl(alpha);
                        if (carbonyl2 != null) {
                            List<Atom> lastAtoms = getAtomListFromPool();
                            lastAtoms.add(carbonyl);
                            lastAtoms.add(nitrogen);
                            lastAtoms.add(alpha);
                            lastAtoms.add(carbonyl2);
                            aa = patternMatch(1, lastAtoms,
                                    PolymerType.AMINOACID);
                            addAtomListToPool(lastAtoms);
                            if (aa != null) {
                                addCap(alpha, carbonyl2, aa);
                                aaArray.add(aa);
                            }
                        }
                    }
                }
            }
            if (aa == null) {
                carbonyl = backbone.get(lastRes + 2);
                alpha = backbone.get(lastRes + 1);
                addCap(alpha, carbonyl, aaArray.get(aaArray.size() - 1));
            }
            int index = 1;
            for (Residue r : aaArray) {
                r.setNumber(index++);
                c.addMSNode(r);
            }
            // Potential DNA/RNA
        } else if (type == PolymerType.NUCLEICACID) {
            Residue base;
            int lastRes = 0;
            boolean firstBase = true;
            Atom phos = null;
            Atom oxygen1 = null;
            Atom phosphate1 = null;
            List<Residue> na = new ArrayList<Residue>();
            while (start < length) {
                base = patternMatch(start, backbone, PolymerType.NUCLEICACID);
                if (base != null) {
                    phos = backbone.get(start - 1);
                    if (phos != null && phos.getAtomicNumber() == 15) {
                        addPhosphate(phos, base);
                    }
                    na.add(base);
                    if (firstBase) {
                        firstBase = false;
                        phosphate1 = backbone.get(start - 1);
                        oxygen1 = backbone.get(start);
                    }
                    lastRes = start;
                }
                start += 6;
            }
            // Make sure the fisrt base is found
            Atom o2, o3;
            Atom c1, c2, c3;
            if (phosphate1 != null && oxygen1 != null) {
                o2 = findOtherOxygen(phosphate1, oxygen1);
                if (o2 != null) {
                    c1 = findBondWith(o2, 6);
                    if (c1 != null) {
                        c2 = findCO(c1);
                        if (c2 != null) {
                            c3 = findC5(c2);
                            if (c3 != null) {
                                o3 = findBondWith(c3, 8);
                                if (o3 != null) {
                                    List<Atom> firstAtoms = getAtomListFromPool();
                                    firstAtoms.add(o3);
                                    firstAtoms.add(c3);
                                    firstAtoms.add(c2);
                                    firstAtoms.add(c1);
                                    firstAtoms.add(o2);
                                    firstAtoms.add(phosphate1);
                                    base = patternMatch(0, firstAtoms, type);
                                    if (base != null) {
                                        addCap(c3, o3, base);
                                        na.add(0, base);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Make sure the last base is found
            oxygen1 = backbone.get(lastRes + 4);
            phosphate1 = backbone.get(lastRes + 5);
            if (phosphate1 != null && oxygen1 != null) {
                o2 = findOtherOxygen(phosphate1, oxygen1);
                if (o2 != null) {
                    c1 = findBondWith(o2, 6);
                    if (c1 != null) {
                        c2 = findBondWith(c1, 6);
                        if (c2 != null) {
                            c3 = findCCO(c2);
                            if (c3 != null) {
                                o3 = findBondWith(c3, 8);
                                if (o3 != null) {
                                    List<Atom> lastAtoms = getAtomListFromPool();
                                    lastAtoms.add(phosphate1);
                                    lastAtoms.add(o2);
                                    lastAtoms.add(c1);
                                    lastAtoms.add(c2);
                                    lastAtoms.add(c3);
                                    lastAtoms.add(o3);
                                    base = patternMatch(1, lastAtoms, type);
                                    if (base != null) {
                                        addPhosphate(phosphate1, base);
                                        addCap(c3, o3, base);
                                        na.add(base);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            int index = 1;
            for (Residue r : na) {
                r.setNumber(index++);
                c.addMSNode(r);
            }
        } else {
            return false;
        }
        return true;
    }

    /**
     * Returns a carbon that is bonded to the atom a, a carbonyl group, and a
     * nitrogen. O=C-alpha-N
     *
     * @param a Atom
     * @return Atom
     */
    public static Atom findAlphaCarbon(Atom a) {
        for (Bond b : a.getBonds()) {
            Atom alpha = b.get1_2(a);
            if (alpha.getAtomicNumber() == 6
                    && findCO(alpha) != null
                    && formsBondsWith(alpha, 7)) {
                return alpha;
            }
        }
        return null;
    }

    /**
     * Returns the first atom with the specified atomic number that bonds with
     * atom a, or null otherwise.
     *
     * @param a Atom
     * @param atomicNumber int
     * @return Atom
     */
    public static Atom findBondWith(Atom a, int atomicNumber) {
        for (Bond b : a.getBonds()) {
            Atom other = b.get1_2(a);
            if (other.getAtomicNumber() == atomicNumber) {
                return other;
            }
        }
        return null;
    }

    /**
     * Returns a carbon that is bonded to the adjacent atom, which bonds 1
     * carbon and 1 oxygen.
     *
     * @param adjacent Atom
     * @return Atom
     */
    public static Atom findC5(Atom adjacent) {
        for (Bond b : adjacent.getBonds()) {
            Atom carbon = b.get1_2(adjacent);
            if (carbon.getAtomicNumber() == 6 && numberOfBondsWith(carbon, 6) == 1 && numberOfBondsWith(carbon, 8) == 1) {
                return carbon;
            }
        }
        return null;
    }

    /**
     * Returns a carbon that is bonded to the adjacent atom and an oxygen.
     *
     * @param adjacent Atom
     * @return Atom
     */
    public static Atom findCarbonyl(Atom adjacent) {
        for (Bond b : adjacent.getBonds()) {
            Atom carbonyl = b.get1_2(adjacent);
            if (carbonyl.getAtomicNumber() == 6) {
                for (Bond b2 : carbonyl.getBonds()) {
                    Atom oxygen = b2.get1_2(carbonyl);
                    if (oxygen.getAtomicNumber() == 8 && oxygen.getBonds().size() == 1) {
                        return carbonyl;
                    }
                }
            }
        }
        return null;
    }

    /**
     * Returns a carbon that is bonded to the adjacent atom, which bonds 2
     * carbons and 1 oxygen.
     *
     * @param adjacent Atom
     * @return Atom
     */
    public static Atom findCCO(Atom adjacent) {
        for (Bond b : adjacent.getBonds()) {
            Atom carbon = b.get1_2(adjacent);
            if (carbon.getAtomicNumber() == 6 && numberOfBondsWith(carbon, 6) == 2 && numberOfBondsWith(carbon, 8) >= 1) {
                return carbon;
            }
        }
        return null;
    }

    /**
     * Returns the first nitrogen atom found that is bonded to the adjacent
     * atom.
     *
     * @param adjacent Atom
     * @return Atom a nitrogen atom.
     */
    public static Atom findN(Atom adjacent) {
        for (Bond b : adjacent.getBonds()) {
            Atom nitrogen = b.get1_2(adjacent);
            if (nitrogen.getAtomicNumber() == 7) {
                return nitrogen;
            }
        }
        return null;
    }

    /**
     * Returns a carbon that is bonded to the adjacent atom, which bonds at
     * least 1 oxygen.
     *
     * @param adjacent Atom
     * @return Atom
     */
    public static Atom findCO(Atom adjacent) {
        for (Bond b : adjacent.getBonds()) {
            Atom carbon = b.get1_2(adjacent);
            if (carbon.getAtomicNumber() == 6 && formsBondsWith(carbon, 8)) {
                return carbon;
            }
        }
        return null;
    }

    /**
     * Returns an oxygen atom that is bonded to atom p and a carbon, but is not
     * atom o. O-P-X-C where X is the returned atom. This is useful for
     * traversing a nucleic acid backbone.
     *
     * @param p Atom
     * @param o Atom
     * @return Atom
     */
    public static Atom findOtherOxygen(Atom p, Atom o) {
        for (Bond b : p.getBonds()) {
            Atom oxygen = b.get1_2(p);
            if (oxygen.getAtomicNumber() == 8 && oxygen != o && formsBondsWith(oxygen, 6)) {
                return oxygen;
            }
        }
        return null;
    }

    /**
     * <p>
     * findPolymer</p>
     *
     * @param atoms List
     * @param currentAtom Atom
     * @param path List
     * @return List
     */
    public static List<Atom> findPolymer(List<Atom> atoms, Atom currentAtom,
            List<Atom> path) {
        // Atom has no bonds to follow
        if (currentAtom.getBonds() == null) {
            path = getAtomListFromPool();
            path.add(currentAtom);
            return path;
        }
        // End of Recursion conditions
        if (currentAtom.getParent() != null) {
            return null;
        }
        int anum = currentAtom.getAtomicNumber();
        // Only C,N,O,P in a DNA/RNA/protein backbone
        if (anum != 6 && anum != 7 && anum != 8 && anum != 15) {
            return null;
        }
        // Allow the search to make it out of side chains, but not enter them...
        if (path != null && path.size() > 7) {
            // Oxygen is only in the backbone for Nucleic Acids in a phosphate
            // group
            if (anum == 8) {
                if (!formsBondsWith(currentAtom, 15)) {
                    return null;
                }
                // Nitrogen is only in the backbone in peptide bonds
            } else if (anum == 7) {
                Atom carbonyl = findCarbonyl(currentAtom);
                if (carbonyl == null) {
                    return null;
                }
                // Avoid more than 3 carbons in a row (phenyl groups, etc.)
            } else if (anum == 6) {
                Atom a;
                int anum2, anum3, anum4;
                int size = path.size();
                a = path.get(size - 1);
                anum2 = a.getAtomicNumber();
                if (anum2 == 6) {
                    a = path.get(size - 2);
                    anum3 = a.getAtomicNumber();
                    if (anum3 == 6) {
                        a = path.get(size - 3);
                        anum4 = a.getAtomicNumber();
                        if (anum4 == 6) {
                            return null;
                        }
                    }
                }
            }
        }
        // Atoms with only one bond are at the end of a Polymer
        Atom previousAtom = null;
        if (path != null) {
            previousAtom = path.get(path.size() - 1);
        }
        List<Bond> bonds = currentAtom.getBonds();
        if (bonds.size() == 1 && previousAtom != null) {
            return null;
        }
        // Initialization
        if (path == null) {
            path = getAtomListFromPool();
            previousAtom = null;
            // Or Continuation
        } else {
            List<Atom> pathclone = getAtomListFromPool();
            pathclone.addAll(path);
            path = pathclone;
        }
        // Add the currentAtom to the growing path
        path.add(currentAtom);
        // Continue search in each bond direction, but no backtracking over
        // previousAtom
        Atom nextAtom;
        List<Atom> newPolymer, maxPolymer = getAtomListFromPool();
        for (Bond b : bonds) {
            nextAtom = b.get1_2(currentAtom);
            // Check to avoid returning in the same direction and loops
            if (nextAtom != previousAtom && !path.contains(nextAtom)) {
                newPolymer = findPolymer(atoms, nextAtom, path);
                if (newPolymer != null) {
                    // Check to see if the Polymers contain any of the same
                    // atoms
                    // and if so, use the shorter Polymer (avoids loops)
                    if (haveCommonAtom(newPolymer, maxPolymer)) {
                        if (newPolymer.size() < maxPolymer.size()) {
                            addAtomListToPool(maxPolymer);
                            maxPolymer = newPolymer;
                        }
                    } else if (newPolymer.size() > maxPolymer.size()) {
                        addAtomListToPool(maxPolymer);
                        maxPolymer = newPolymer;
                    }
                }
            }
        }
        // Add the currentAtom to the longest discovered chain and return
        maxPolymer.add(0, currentAtom);
        return maxPolymer;
    }

    /**
     * Returns an atom bonded to the "end" atom, which is not equal to "other".
     *
     * @param end Atom
     * @param other Atom
     * @return Atom
     */
    public static Atom findSeed(Atom end, Atom other) {
        for (Bond b : end.getBonds()) {
            Atom seed = b.get1_2(end);
            if (seed != other) {
                return seed;
            }
        }
        return null;
    }

    /**
     * True if Atom a forms a bond with another atom of the specified atomic
     * Number.
     *
     * @param a Atom
     * @param atomicNumber int
     * @return boolean
     */
    public static boolean formsBondsWith(Atom a, int atomicNumber) {
        for (Bond b : a.getBonds()) {
            Atom other = b.get1_2(a);
            if (other.getAtomicNumber() == atomicNumber) {
                return true;
            }
        }
        return false;
    }

    /**
     * <p>
     * getAtomListFromPool</p>
     *
     * @return a {@link java.util.List} object.
     */
    public static List<Atom> getAtomListFromPool() {
        if (atomListPool.isEmpty()) {
            return new ArrayList<Atom>();
        }
        return atomListPool.remove(0);
    }

    /**
     * Returns true if the lists contain any atom in common.
     *
     * @param list1 List
     * @param list2 List
     * @return boolean
     */
    private static boolean haveCommonAtom(List<Atom> list1, List<Atom> list2) {
        if (list1 == null || list2 == null) {
            return false;
        }
        for (Atom a : list1) {
            if (list2.contains(a)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Returns true if Atom a is water oxygen.
     *
     * @param a Atom
     * @return boolean
     */
    public static boolean isWaterOxygen(Atom a) {
        if (a.getAtomicNumber() != 8) {
            return false;
        }
        for (Bond b : a.getBonds()) {
            Atom h = b.get1_2(a);
            if (h.getAtomicNumber() != 1) {
                return false;
            }
        }
        return true;
    }

    /**
     * This function returns the number of times atom "a" is bonded to an atom
     * of the specified atomic number.
     *
     * @param a Atom
     * @param atomicNumber int
     * @return int
     */
    public static int numberOfBondsWith(Atom a, int atomicNumber) {
        int total = 0;
        for (Bond b : a.getBonds()) {
            Atom other = b.get1_2(a);
            if (other.getAtomicNumber() == atomicNumber) {
                total++;
            }
        }
        return total;
    }

    // Check to see if a portion of the backbone matches that of a
    // biological polymer, and if so determine the respective residue
    private static Residue patternMatch(int start, List<Atom> backbone,
            PolymerType type) {
        int pattern[];
        // Initialization
        if (type == PolymerType.AMINOACID) {
            pattern = AAPATTERN;
            if (backbone.size() < start + pattern.length) {
                return null;
            }
            // Check for correct Carbonyl placement
            Atom a = backbone.get(start + 1);
            if (formsBondsWith(a, 8)) {
                return null;
            }
            a = backbone.get(start + 2);
            if (!formsBondsWith(a, 8)) {
                return null;
            }
        } else if (type == PolymerType.NUCLEICACID) {
            pattern = NAPATTERN;
            if (backbone.size() < start + pattern.length) {
                return null;
            }
        } else {
            return null;
        }
        int length = pattern.length;
        List<Atom> atoms = getAtomListFromPool();
        List<Atom> sidePolymer = getAtomListFromPool();
        for (int i = 0; i < length; i++) {
            Atom a = backbone.get(start + i);
            // add backbone atoms to terminate sidePolymer search
            sidePolymer.add(a);
            if (a.getAtomicNumber() != pattern[i]) {
                return null;
            }
        }
        // Collect all the atoms in the Residue
        // Add the atom before and after the backbone pattern to
        // terminate the search, then remove them
        if (start > 0) {
            atoms.add(backbone.get(start - 1));
        }
        if (start + length < backbone.size()) {
            atoms.add(backbone.get(start + length));
        }
        collectAtoms(backbone.get(start), atoms);
        if (start > 0) {
            atoms.remove(0);
        }
        if (start + length < backbone.size()) {
            atoms.remove(0);
            // Collect Just Side chain atoms, then remove backbone termination
        }
        if (type == PolymerType.AMINOACID) {
            collectAtoms(sidePolymer.get(1), sidePolymer);
        } else if (type == PolymerType.NUCLEICACID) {
            Atom seed = null;
            for (Atom a : atoms) {
                if (a.getAtomicNumber() == 7) {
                    seed = a;
                    break;
                }
            }
            if (seed != null && seed.getAtomicNumber() == 7) {
                sidePolymer.add(seed);
                collectAtoms(seed, sidePolymer);
            } else {
                return null;
            }
        }
        for (int i = 0; i <= length; i++) {
            sidePolymer.remove(0);
        }
        Residue res = assignResidue(backbone, start, atoms, sidePolymer);
        return res;
    }

    /**
     * Determine chainID for a given polymer number.
     *
     * @param i int
     * @return Character
     */
    public static Character getChainID(int i) {
        if (i > 35) {
            i = i % 36;
        }
        Character c = null;
        if (i < 26) {
            /**
             * 65 is 'A'. 90 is 'Z'.
             */
            c = Character.valueOf((char) (i + 65));
        } else {
            i -= 26;
            /**
             * 48 is '0'. 57 is '9'.
             */
            c = Character.valueOf((char) (i + 48));
        }
        return c;
    }

    /**
     * Returns an List with reversed ordering.
     *
     * @param atomList List
     * @return List
     */
    static private List<Atom> reverseAtomList(List<Atom> atomList) {
        List<Atom> reversedList = getAtomListFromPool();
        for (Atom a : atomList) {
            reversedList.add(0, a);
        }
        return reversedList;
    }

    /**
     * Finds the RMS deviation between the atoms of MolecularAssembly m1 and m2
     * provided they have the same number of atoms.
     *
     * @param m1 a {@link ffx.potential.MolecularAssembly} object.
     * @param m2 a {@link ffx.potential.MolecularAssembly} object.
     * @return a double.
     */
    public static double RMSCoordDev(MolecularAssembly m1, MolecularAssembly m2) {
        if (m1 == null || m2 == null) {
            return 0;
        }
        int n1 = m1.getAtomList().size();
        int n2 = m2.getAtomList().size();
        if (n1 != n2) {
            return 0;
        }
        Atom a1, a2;
        double[] d = new double[3];
        double[] da = new double[3];
        double[] db = new double[3];
        double rms = 0;
        ListIterator li, lj;
        for (li = m1.getAtomList().listIterator(), lj = m2.getAtomList().listIterator(); li.hasNext();) {
            a1 = (Atom) li.next();
            a2 = (Atom) lj.next();
            a1.getXYZ(da);
            a2.getXYZ(db);
            diff(da, db, d);
            rms += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        }
        return Math.sqrt(rms / n1);
    }
}
