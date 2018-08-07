/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.potential.bonded;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Residue.ResidueType;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import ffx.potential.parameters.AngleType;

import static ffx.potential.bonded.BondedUtils.determineIntxyz;
import static ffx.potential.bonded.BondedUtils.intxyz;

/**
 * The Rotamer Library Class manages a library of side-chain Rotamers for amino
 * acids, and a library of backbone Rotamers for nucleic acids.
 *
 * @author Ava M. Lynn
 * @author Shibo Gao
 * @author Jacob M. Litman
 */
public class RotamerLibrary {

    private static final Logger logger = Logger.getLogger(RotamerLibrary.class.getName());
    /**
     * Number of amino acid residues types currently recognized, although there
     * are not rotamer libraries for each yet.
     */
    private static final int numberOfAminoAcids = AminoAcid3.values().length;
    /**
     * Number of nucleic acid residues types currently recognized, although
     * there are not rotamer libraries for each yet.
     */
    private static final int numberOfNucleicAcids = NucleicAcid3.values().length;
    /**
     * The first time rotamers are requested for an amino acid type, they are
     * instantiated into an array, which is stored in the cache. Subsequently
     * the reference is simply returned.
     */
    private final Rotamer[][] aminoAcidRotamerCache = new Rotamer[numberOfAminoAcids][];
    /**
     * The first time rotamers are requested for a nucleic acid type, they are
     * instantiated into an array, which is stored in the cache. Subsequently
     * the reference is simply returned.
     */
    private final Rotamer[][] nucleicAcidRotamerCache = new Rotamer[numberOfNucleicAcids][];

    private ProteinLibrary proteinLibrary = ProteinLibrary.PonderAndRichards;
    
    // Can easily add an naSolvationLibrary when we want to do NA sequence optimization.
    private boolean useOrigCoordsRotamer = false;
    
    // May have to become instance variable, would break staticness of measureRotamer methods.
    private static final Map<String, NonstandardRotLibrary> nonstdRotCache = new HashMap<>();
    
    private static final RotamerLibrary defaultRotamerLibrary = new RotamerLibrary(ProteinLibrary.PonderAndRichards, false);

    private static final boolean useIdealRingGeometries;
    static {
        useIdealRingGeometries = Boolean.parseBoolean(System.getProperty("useIdealRingGeo", "true"));
    }

    private static final Map<AminoAcid3, Map<String, Double>> idealAngleGeometries;
    static {
        Map<AminoAcid3, Map<String, Double>> angleGeos = new HashMap<>();

        Map<String, Double> pheMap = new HashMap<>();
        pheMap.put("CD1-CG-CB", 120.3);
        pheMap.put("CD2-CG-CB", 120.3);
        pheMap.put("CE1-CD1-CG", 120.3);
        pheMap.put("CE2-CD2-CG", 120.3);
        pheMap.put("CZ-CE1-CD1", 120.0);
        pheMap.put("CZ-CE2-CD2", 120.0);
        angleGeos.put(AminoAcid3.PHE, Collections.unmodifiableMap(pheMap));

        Map<String, Double> tyrMap = new HashMap<>();
        tyrMap.put("CD1-CG-CB", 120.3);
        tyrMap.put("CD2-CG-CB", 120.3);
        tyrMap.put("CE1-CD1-CG", 120.3);
        tyrMap.put("CE2-CD2-CG", 120.3);
        tyrMap.put("CZ-CE1-CD1", 120.0);
        tyrMap.put("CZ-CE2-CD2", 120.0);
        angleGeos.put(AminoAcid3.TYR, Collections.unmodifiableMap(tyrMap));

        Map<String, Double> tydMap = new HashMap<>();
        tydMap.put("CD1-CG-CB", 120.5);
        tydMap.put("CD2-CG-CB", 120.5);
        tydMap.put("CE1-CD1-CG", 120.4);
        tydMap.put("CE2-CD2-CG", 120.4);
        tydMap.put("CZ-CE1-CD1", 120.8);
        tydMap.put("CZ-CE2-CD2", 120.8);
        angleGeos.put(AminoAcid3.TYD, Collections.unmodifiableMap(tydMap));
        
        Map<String, Double> hisMap = new HashMap<>();
        hisMap.put("ND1-CG-CB", 122.1);
        hisMap.put("CD2-CG-CB", 131.0);
        hisMap.put("CD2-CG-ND1", 106.8);
        hisMap.put("CE1-ND1-CG", 109.5);
        hisMap.put("NE2-CD2-CG", 107.1);
        angleGeos.put(AminoAcid3.HIS, Collections.unmodifiableMap(hisMap));

        Map<String, Double> hidMap = new HashMap<>();
        hidMap.put("ND1-CG-CB", 123.5);
        hidMap.put("CD2-CG-CB", 132.3);
        hidMap.put("CD2-CG-ND1", 104.2);
        hidMap.put("CE1-ND1-CG", 108.8);
        hidMap.put("NE2-CD2-CG", 111.2);
        angleGeos.put(AminoAcid3.HID, Collections.unmodifiableMap(hidMap));

        Map<String, Double> hieMap = new HashMap<>();
        hieMap.put("ND1-CG-CB", 120.2);
        hieMap.put("CD2-CG-CB", 129.1);
        hieMap.put("CD2-CG-ND1", 110.7);
        hieMap.put("CE1-ND1-CG", 105.1);
        hieMap.put("NE2-CD2-CG", 104.6);
        angleGeos.put(AminoAcid3.HIE, Collections.unmodifiableMap(hieMap));

        Map<String, Double> trpMap = new HashMap<>();
        trpMap.put("CD1-CG-CB", 126.4);
        trpMap.put("CD2-CG-CB", 126.5);
        trpMap.put("CD2-CG-CD1", 107.1);
        trpMap.put("NE1-CD1-CG", 106.9);
        trpMap.put("CE2-NE1-CD1", 109.4);
        trpMap.put("CE3-CD2-CE2", 121.6);
        trpMap.put("CZ2-CE2-CD2", 123.5);
        trpMap.put("CZ3-CE3-CD2", 116.7);
        trpMap.put("CH2-CZ2-CE2", 116.2);
        angleGeos.put(AminoAcid3.TRP, Collections.unmodifiableMap(trpMap));

        idealAngleGeometries = Collections.unmodifiableMap(angleGeos);
    }
    
    public RotamerLibrary(ProteinLibrary name, boolean origCoords) {
        proteinLibrary = name;
        useOrigCoordsRotamer = origCoords;
    }
    
    public static RotamerLibrary getDefaultLibrary() {
        return defaultRotamerLibrary;
    }

    /**
     * Set the protein rotamer library to use.
     *
     * @param name the ProteinLibrary to use.
     */
    public void setLibrary(ProteinLibrary name) {
        proteinLibrary = name;
        for (int i = 0; i < numberOfAminoAcids; i++) {
            aminoAcidRotamerCache[i] = null;
        }
    }

    public boolean getUsingOrigCoordsRotamer() {
        return useOrigCoordsRotamer;
    }

    public void setUseOrigCoordsRotamer(boolean set) {
        useOrigCoordsRotamer = set;
    }

    /**
     * Get the protein rotamer library.
     *
     * @return the ProteinLibrary in use.
     */
    public ProteinLibrary getLibrary() {
        return proteinLibrary;
    }

    /**
     * Return rotamer array for the given AA or NA residue.
     *
     * @param residue the Residue to examine.
     * @return Array of Rotamers for Residue's type.
     */
    Rotamer[] getRotamers(Residue residue) {
        // Package-private; intended to be accessed only by Residue and extensions
        // thereof. Otherwise, use Residue.getRotamers(RotamerLibrary library).
        if (residue == null) {
            return null;
        }
        if (isModRes(residue)) {
            if(nonstdRotCache.containsKey(residue.getName().toUpperCase())) {
                return nonstdRotCache.get(residue.getName().toUpperCase()).getRotamers();
            }
            return null;
        }
        switch (residue.getResidueType()) {
            case AA:
                AminoAcid3 aa = AminoAcid3.valueOf(residue.getName());
                return getRotamers(aa);
            case NA:
                NucleicAcid3 na = NucleicAcid3.valueOf(residue.getName());
                return getRotamers(na);
            default:
                if (nonstdRotCache.containsKey(residue.getName().toUpperCase())) {
                    return nonstdRotCache.get(residue.getName().toUpperCase()).getRotamers();
                }
                return null;
        }
    }
    
    /**      
     * Checks if a Residue has a PTM.
     * @param residue
     * @return 
     */
    private static boolean isModRes(Residue residue) {
        List<Atom> resatoms = residue.getAtomList();
        return (resatoms != null && !resatoms.isEmpty() && resatoms.get(0).isModRes());
    }

    /**
     * Return an array of Rotamers for the given amino acid.
     *
     * @param name The name of the amino acid.
     * @return An array of Rotamers.
     */
    public Rotamer[] getRotamers(AminoAcid3 name) {
        switch (proteinLibrary) {
            case PonderAndRichards:
                return getPonderAndRichardsRotamers(name);
            case Richardson:
                return getRichardsonRotamers(name);
            case None:
            default:
                return null;
        }
    }

    /**
     * Return an array of Rotamers for the given nucleic acid.
     *
     * @param name The name of the nucleic acid.
     * @return An array of Rotamers.
     */
    public Rotamer[] getRotamers(NucleicAcid3 name) {
        return getRichardsonRNARotamers(name);
    }

    /**
     * Initializes default coordinates (presently PDB coordinates) for key atoms
     * in all nucleic acid Residues. This is necessary to preserve rotamer
     * independence while still adjusting rotamers to correctly meet prior
     * residues; C4', O4', and C1' are used to build the rest of the nucleic
     * acid base, but are moved with each Rotamer to take part of the strain of
     * correctly meeting O3' of residue i-1.
     *
     * It also initializes the location of O3' in both North and South puckers;
     * while in theory this could be recalculated each time based off of C4',
     * O4', and C1', it is easier to just store these two locations and call
     * them when needed.
     *
     * This MUST be called before any applyRotamer calls are made, else invalid
     * coordinates will be stored.
     *
     * @param polymers the Polymer array to examine.
     */
    public static void initializeDefaultAtomicCoordinates(Polymer[] polymers) {
        for (Polymer polymer : polymers) {
            ArrayList<Residue> current = polymer.getResidues();
            for (Residue residuej : current) {
                switch (residuej.getResidueType()) {
                    case NA:
                        residuej.initializeDefaultAtomicCoordinates();
                        break;
                    default:
                        break;
                }
            }
        }
    }

    /**
     * TODO: Add reference to Ponder & Richard's original paper.
     *
     * @param name Type of amino acid.
     * @return Rotamer cache (double[] of torsions).
     */
    private Rotamer[] getPonderAndRichardsRotamers(AminoAcid3 name) {
        int n = name.ordinal();
        if (aminoAcidRotamerCache[n] != null) {
            return aminoAcidRotamerCache[n];
        }
        switch (name) {
            case VAL:
                aminoAcidRotamerCache[n] = new Rotamer[3];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 173.5, 9.0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -63.4, 8.1);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 69.3, 9.6);
                break;
            case LEU:
                aminoAcidRotamerCache[n] = new Rotamer[4];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -64.9, 8.2, 176.0, 9.9);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -176.4, 10.2, 63.1, 8.2);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -165.3, 10.0, 168.2, 34.2);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, 44.3, 20.0, 60.4, 18.8);
                break;
            case ILE:
                aminoAcidRotamerCache[n] = new Rotamer[5];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -60.9, 7.5, 168.7, 11.6);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -59.6, 9.6, -64.1, 14.3);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 61.7, 5.0, 163.8, 16.4);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -166.6, 10.1, 166.0, 8.9);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -174.8, 24.9, 72.1, 10.5);
                break;
            case SER:
                aminoAcidRotamerCache[n] = new Rotamer[3];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 64.7, 16.1);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -69.7, 14.6);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -176.1, 20.2);
                break;
            case THR:
                aminoAcidRotamerCache[n] = new Rotamer[3];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62.7, 8.5);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -59.7, 9.4);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -169.5, 6.6);
                break;
            case CYS:
            case CYD:
                aminoAcidRotamerCache[n] = null;
//                rotamerCache[n] = new Rotamer[3];
//                rotamerCache[n][0] = new Rotamer(name, -65.2, 10.1);
//                rotamerCache[n][1] = new Rotamer(name, -179.6, 9.5);
//                rotamerCache[n][2] = new Rotamer(name, 63.5, 9.6);
                break;
            /*
             * TODO: Figure out proline rotamers.  I have dihedrals from
             * the Richardson lab/Kinemages website (downloaded PDB, used
             * the get_dihedral function to extract the dihedrals), and they
             * conflict with these rotamers.  Plus, this library only
             * specifies one of the two necessary dihedrals.
             */
            case PRO:
                aminoAcidRotamerCache[n] = new Rotamer[3];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 24.0, 8.0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 0.0, 8.0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -24.0, 8.0);
                /*
                 aminoAcidRotamerCache[n] = new Rotamer[2];
                 // The exo- conformation
                 aminoAcidRotamerCache[n][0] = new Rotamer(name, -27.9, 0, 39.0, 0);
                 // The endo- conformation
                 aminoAcidRotamerCache[n][2] = new Rotamer(name, 25.6, 0, -35.8, 0);
                 */
                break;
            case PHE:
                aminoAcidRotamerCache[n] = new Rotamer[4];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -66.3, 10.2, 94.3, 19.5);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -179.2, 9.3, 78.9, 8.9);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 66.0, 12.0, 90.7, 9.4);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -71.9, 16.3, -0.4, 26.1);
                break;
            case TYR:
                aminoAcidRotamerCache[n] = new Rotamer[8];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -66.5, 11.4, 96.6, 21.8, 0, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -179.7, 12.6, 71.9, 13.4, 0, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 63.3, 9.4, 89.1, 13.0, 0, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -67.2, 13.2, -1.0, 20.1, 0, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -66.5, 11.4, 96.6, 21.8, 180.0, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -179.7, 12.6, 71.9, 13.4, 180.0, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, 63.3, 9.4, 89.1, 13.0, 180.0, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -67.2, 13.2, -1.0, 20.1, 180.0, 0);
                break;
            case TYD:
                aminoAcidRotamerCache[n] = new Rotamer[4];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -66.5, 11.4, 96.6, 21.8);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -179.7, 12.6, 71.9, 13.4);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 63.3, 9.4, 89.1, 13.0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -67.2, 13.2, -1.0, 20.1);
                break;
            case TRP:
                aminoAcidRotamerCache[n] = new Rotamer[6];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -70.4, 7.0, 100.5, 18.2);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 64.8, 13.0, -88.9, 5.3);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177.3, 7.9, -95.1, 7.6);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -179.5, 3.4, 87.5, 3.8);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -73.3, 6.5, -87.7, 8.1);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 62.2, 10.0, 112.5, 15.0);
                break;
            case HIS:
            case HIE:
            case HID:
                aminoAcidRotamerCache[n] = new Rotamer[6];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -62.8, 10.0, -74.3, 17.2);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -175.2, 15.4, -88.7, 43.5);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -69.8, 5.9, 96.1, 32.2);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, 67.9, 17.4, -80.5, 40.7);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -177.3, 6.3, 100.5, 14.0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 48.8, 10.0, 89.5, 30.0);
                break;
            case ASH:
                aminoAcidRotamerCache[n] = new Rotamer[6];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -68.3, 9.2, -25.7, 31.1);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -169.1, 9.5, 3.9, 38.9);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 63.7, 9.9, 2.4, 29.4);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -68.3, 9.2, 154.3, 31.1);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -169.1, 9.5, -176.1, 38.9);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 63.7, 9.9, -177.6, 29.4);
                break;
            case ASP:
                aminoAcidRotamerCache[n] = new Rotamer[3];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -68.3, 9.2, -25.7, 31.1);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -169.1, 9.5, 3.9, 38.9);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 63.7, 9.9, 2.4, 29.4);
                break;
            case ASN:
                aminoAcidRotamerCache[n] = new Rotamer[6];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -68.3, 12.3, -36.8, 25.2);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -177.1, 8.8, 1.3, 34.1);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -67.2, 10.8, 128.8, 24.2);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, 63.9, 3.7, -6.8, 13.5);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -174.9, 17.9, -156.8, 58.9);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 63.6, 6.6, 53.8, 17.1);
                break;
            case GLU:
                aminoAcidRotamerCache[n] = new Rotamer[7];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -69.6, 19.2, -177.2, 21.7, -11.4, 44.8);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -176.2, 14.9, 175.4, 10.6, -6.7, 39.0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -64.6, 13.5, -69.1, 17.3, -33.4, 27.4);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -55.6, 10.6, 77.0, 6.8, 25.3, 32.6);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, 69.8, 10.6, -179.0, 23.7, 6.6, 64.2);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -173.6, 14.6, 70.6, 8.7, 14.0, 37.1);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, 63.0, 4.3, -80.4, 13.9, 16.3, 20.8);
                break;
            case GLH:
                aminoAcidRotamerCache[n] = new Rotamer[14];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -69.6, 19.2, -177.2, 21.7, -11.4, 44.8);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -176.2, 14.9, 175.4, 10.6, -6.7, 39.0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -64.6, 13.5, -69.1, 17.3, -33.4, 27.4);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -55.6, 10.6, 77.0, 6.8, 25.3, 32.6);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, 69.8, 10.6, -179.0, 23.7, 6.6, 64.2);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -173.6, 14.6, 70.6, 8.7, 14.0, 37.1);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, 63.0, 4.3, -80.4, 13.9, 16.3, 20.8);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -69.6, 19.2, -177.2, 21.7, 168.6, 44.8);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -176.2, 14.9, 175.4, 10.6, 175.3, 39.0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -64.6, 13.5, -69.1, 17.3, 146.6, 27.4);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -55.6, 10.6, 77.0, 6.8, -154.7, 32.6);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, 69.8, 10.6, -179.0, 23.7, -173.4, 64.2);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -173.6, 14.6, 70.6, 8.7, -166.0, 37.1);
                aminoAcidRotamerCache[n][13] = new Rotamer(name, 63.0, 4.3, -80.4, 13.9, -163.7, 20.8);
                break;
            case GLN:
                aminoAcidRotamerCache[n] = new Rotamer[10];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -66.7, 14.1, -178.5, 14.9, -24.0, 38.0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -66.7, 14.1, -178.5, 14.9, 156.0, 38.0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -174.6, 11.5, -177.7, 17.2, -24.0, 38.0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -174.6, 11.5, -177.7, 17.2, 156.0, 38.0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -58.7, 11.2, -63.8, 16.1, -46.3, 27.7);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -51.3, 7.3, -90.4, 22.8, 165.0, 38.2);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -179.4, 21.5, 67.3, 7.9, 26.8, 38.4);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, 167.5, 14.8, 70.9, 3.7, 174.2, 7.1);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, 70.8, 13.0, -165.6, 9.5, -24.0, 38.0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, 70.8, 13.0, -165.6, 9.5, 156.0, 38.0);
                break;
            case MET:
                aminoAcidRotamerCache[n] = new Rotamer[13];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -64.5, 12.7, -68.5, 6.0, -75.6, 14.1);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -78.3, 5.4, -174.7, 15.7, 65.0, 20.0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -78.3, 5.4, -174.7, 15.7, 180.0, 20.0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -78.3, 5.4, -174.7, 15.7, -65.0, 20.0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, 178.9, 8.7, 179.0, 13.4, 65.0, 20.0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 178.9, 8.7, 179.0, 13.4, 65.0, 20.0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, 178.9, 8.7, 179.0, 13.4, -65.0, 20.0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 65.0, 20.0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -170.0, 24.0, 65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -170.0, 24.0, -65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -70.0, 21.0, 65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, 61.0, 21.0, 65.0, 20.0, 180.0, 20.0);
                break;
            case LYS:
            case LYD:
                aminoAcidRotamerCache[n] = new Rotamer[12];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 65.0, 20.0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -170.0, 24.0, -65.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -70.0, 21.0, 65.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, -65.0, 20.0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, 180.0, 20.0, -65.0, 20.0);
                break;
            case ARG:
                aminoAcidRotamerCache[n] = new Rotamer[14];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 61.0, 25.0, 180.0, 20.0, 65.0, 20.0, 90.0, 20.0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 61.0, 25.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 65.0, 20.0, 90.0, 20.0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 90.0, 20.0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, 180.0, 20.0, -90.0, 20.0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -170.0, 24.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 65.0, 20.0, 90.0, 20.0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, 180.0, 20.0, -90.0, 20.0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -70.0, 21.0, 180.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -170.0, 21.0, 65.0, 20.0, 65.0, 20.0, 180.0, 20.0);
                aminoAcidRotamerCache[n][13] = new Rotamer(name, -70.0, 21.0, -65.0, 20.0, -65.0, 20.0, 180.0, 20.0);
                break;
            default:
                // Handles GLY, ALA, CYX, ...
                break;
        }
        return aminoAcidRotamerCache[n];
    }

    /**
     * TODO: Add reference to Richardson et al.
     *
     * @param name Type of amino acid.
     * @return Rotamer cache (double[] of torsions).
     */
    private Rotamer[] getRichardsonRotamers(AminoAcid3 name) {
        int n = name.ordinal();
        if (aminoAcidRotamerCache[n] != null) {
            return aminoAcidRotamerCache[n];
        }
        switch (name) {
            case VAL:
                aminoAcidRotamerCache[n] = new Rotamer[3];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 64, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 175, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -60, 0);
                break;
            case LEU:
                aminoAcidRotamerCache[n] = new Rotamer[5];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 80, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -177, 0, 65, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -172, 0, 145, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -85, 0, 65, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -65, 0, 175, 0);
                break;
            case ILE:
                aminoAcidRotamerCache[n] = new Rotamer[7];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 100, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 170, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 66, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 165, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -65, 0, 100, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, 170, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -57, 0, -60, 0);
                break;
            case SER:
                aminoAcidRotamerCache[n] = new Rotamer[18];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 0, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 60, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 62, 0, 120, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, 62, 0, -60, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 62, 0, -120, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -177, 0, 0, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -177, 0, 60, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -177, 0, 120, 0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -177, 0, 180, 0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -177, 0, -60, 0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -177, 0, -120, 0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -65, 0, 0, 0);
                aminoAcidRotamerCache[n][13] = new Rotamer(name, -65, 0, 60, 0);
                aminoAcidRotamerCache[n][14] = new Rotamer(name, -65, 0, 120, 0);
                aminoAcidRotamerCache[n][15] = new Rotamer(name, -65, 0, 180, 0);
                aminoAcidRotamerCache[n][16] = new Rotamer(name, -65, 0, -60, 0);
                aminoAcidRotamerCache[n][17] = new Rotamer(name, -65, 0, -120, 0);
                break;
            case THR:
                aminoAcidRotamerCache[n] = new Rotamer[18];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 0, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 60, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 62, 0, 120, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, 62, 0, -60, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 62, 0, -120, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -175, 0, 0, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -175, 0, 60, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -175, 0, 120, 0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -175, 0, 180, 0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -175, 0, -60, 0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -175, 0, -120, 0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -65, 0, 0, 0);
                aminoAcidRotamerCache[n][13] = new Rotamer(name, -65, 0, 60, 0);
                aminoAcidRotamerCache[n][14] = new Rotamer(name, -65, 0, 120, 0);
                aminoAcidRotamerCache[n][15] = new Rotamer(name, -65, 0, 180, 0);
                aminoAcidRotamerCache[n][16] = new Rotamer(name, -65, 0, -60, 0);
                aminoAcidRotamerCache[n][17] = new Rotamer(name, -65, 0, -120, 0);
                break;
            case CYS:
            case CYD:
                aminoAcidRotamerCache[n] = null;
//                rotamerCache[n] = new Rotamer[3];
//                rotamerCache[n][0] = new Rotamer(name, 62, 0);
//                rotamerCache[n][1] = new Rotamer(name, -177, 0);
//                rotamerCache[n][2] = new Rotamer(name, -65, 0);
                break;
            case PHE:
                aminoAcidRotamerCache[n] = new Rotamer[4];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 90, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -177, 0, 80, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -65, 0, -85, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -65, 0, -30, 0);
                break;
            case TYR:
                aminoAcidRotamerCache[n] = new Rotamer[8];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 90, 0, 0, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 90, 0, 180, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 80, 0, 0, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 80, 0, 180, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -65, 0, -85, 0, 0, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, -85, 0, 180, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -65, 0, -30, 0, 0, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -65, 0, -30, 0, 180, 0);
                break;
            case TYD:
                aminoAcidRotamerCache[n] = new Rotamer[4];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 90, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, -177, 0, 80, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -65, 0, -85, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -65, 0, -30, 0);
                break;
            case TRP:
                aminoAcidRotamerCache[n] = new Rotamer[7];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, -90, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 90, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, -105, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 90, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -65, 0, -90, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, -5, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -65, 0, 95, 0);
                break;
            case HIS:
            case HIE:
            case HID:
                aminoAcidRotamerCache[n] = new Rotamer[8];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, -75, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 80, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, -165, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, -80, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -177, 0, 60, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, -70, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -65, 0, 165, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -65, 0, 80, 0);
                break;
            case ASP:
                aminoAcidRotamerCache[n] = new Rotamer[5];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 10, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 30, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 0, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 65, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -70, 0, -15, 0);
                break;
            case ASH:
                aminoAcidRotamerCache[n] = new Rotamer[10];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 10, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 30, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 0, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 65, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -70, 0, -15, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 62, 0, -170, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, 62, 0, -150, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -177, 0, -180, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -177, 0, -115, 0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -70, 0, 165, 0);
                break;
            case ASN:
                aminoAcidRotamerCache[n] = new Rotamer[7];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, -10, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 30, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -174, 0, -20, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 30, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -65, 0, -20, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, -75, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -65, 0, 120, 0);
                break;
            case GLU:
                aminoAcidRotamerCache[n] = new Rotamer[8];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, -20, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 70, 0, -80, 0, 0, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 65, 0, 10, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 180, 0, 0, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -177, 0, -80, 0, -25, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, 85, 0, 0, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -67, 0, -180, 0, -10, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -65, 0, -65, 0, -40, 0);
                break;
            case GLH:
                aminoAcidRotamerCache[n] = new Rotamer[16];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, -20, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 70, 0, -80, 0, 0, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 65, 0, 10, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 180, 0, 0, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -177, 0, -80, 0, -25, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, 85, 0, 0, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -67, 0, -180, 0, -10, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -65, 0, -65, 0, -40, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, 62, 0, 180, 0, 160, 0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, 70, 0, -80, 0, -180, 0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -177, 0, 65, 0, -170, 0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -177, 0, 180, 0, -180, 0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -177, 0, -80, 0, 155, 0);
                aminoAcidRotamerCache[n][13] = new Rotamer(name, -65, 0, 85, 0, -180, 0);
                aminoAcidRotamerCache[n][14] = new Rotamer(name, -67, 0, -180, 0, 170, 0);
                aminoAcidRotamerCache[n][15] = new Rotamer(name, -65, 0, -65, 0, 140, 0);
                break;
            case GLN:
                aminoAcidRotamerCache[n] = new Rotamer[9];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 20, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 70, 0, -75, 0, 0, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 65, 0, -100, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 65, 0, 60, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -177, 0, 180, 0, 0, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -65, 0, 85, 0, 0, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -67, 0, 180, 0, -25, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -65, 0, -65, 0, -40, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -65, 0, -65, 0, 100, 0);
                break;
            case MET:
                aminoAcidRotamerCache[n] = new Rotamer[13];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 75, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 180, 0, -75, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, -177, 0, 65, 0, 75, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, -177, 0, 65, 0, 180, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, -177, 0, 180, 0, 75, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -177, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -177, 0, 180, 0, -75, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -67, 0, 180, 0, 75, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -67, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -67, 0, 180, 0, -75, 0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -65, 0, -65, 0, 103, 0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -65, 0, -65, 0, 180, 0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -65, 0, -65, 0, -70, 0);
                break;
            case LYS:
            case LYD:
                aminoAcidRotamerCache[n] = new Rotamer[27];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 68, 0, 180, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 65.0, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0, 180, 0, -65, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, 62, 0, 180, 0, -68, 0, 180, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, -177, 0, 68, 0, 180, 0, 65, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, -177, 0, 68, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -177, 0, 68, 0, 180, 0, -65, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -177, 0, 180, 0, 68, 0, 65, 0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -177, 0, 180, 0, 68, 0, 180, 0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 65, 0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -177, 0, 180, 0, 180, 0, -65, 0);
                aminoAcidRotamerCache[n][13] = new Rotamer(name, -177, 0, 180, 0, -68, 0, 180, 0);
                aminoAcidRotamerCache[n][14] = new Rotamer(name, -177, 0, 180, 0, -68, 0, -65, 0);
                aminoAcidRotamerCache[n][15] = new Rotamer(name, -90, 0, 68, 0, 180, 0, 180);
                aminoAcidRotamerCache[n][16] = new Rotamer(name, -67, 0, 180, 0, 68, 0, -65, 0);
                aminoAcidRotamerCache[n][17] = new Rotamer(name, -67, 0, 180, 0, 68, 0, 180, 0);
                aminoAcidRotamerCache[n][18] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 65, 0);
                aminoAcidRotamerCache[n][19] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][20] = new Rotamer(name, -67, 0, 180, 0, 180, 0, -65, 0);
                aminoAcidRotamerCache[n][21] = new Rotamer(name, -67, 0, 180, 0, -68, 0, 180, 0);
                aminoAcidRotamerCache[n][22] = new Rotamer(name, -67, 0, 180, 0, -68, 0, -65, 0);
                aminoAcidRotamerCache[n][23] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 65, 0);
                aminoAcidRotamerCache[n][24] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][25] = new Rotamer(name, -62, 0, -68, 0, 180, 0, -65, 0);
                aminoAcidRotamerCache[n][26] = new Rotamer(name, -62, 0, -68, 0, -68, 0, 180, 0);
                break;
            case ARG:
                aminoAcidRotamerCache[n] = new Rotamer[34];
                aminoAcidRotamerCache[n][0] = new Rotamer(name, 62, 0, 180, 0, 65, 0, 85, 0);
                aminoAcidRotamerCache[n][1] = new Rotamer(name, 62, 0, 180, 0, 65, 0, -175, 0);
                aminoAcidRotamerCache[n][2] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 85, 0);
                aminoAcidRotamerCache[n][3] = new Rotamer(name, 62, 0, 180, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][4] = new Rotamer(name, 62, 0, 180, 0, 180, 0, -85, 0);
                aminoAcidRotamerCache[n][5] = new Rotamer(name, 62, 0, 180, 0, -65, 0, 175, 0);
                aminoAcidRotamerCache[n][6] = new Rotamer(name, 62, 0, 180, 0, -65, 0, -85, 0);
                aminoAcidRotamerCache[n][7] = new Rotamer(name, -177, 0, 65, 0, 65, 0, 85, 0);
                aminoAcidRotamerCache[n][8] = new Rotamer(name, -177, 0, 65, 0, 65, 0, -175, 0);
                aminoAcidRotamerCache[n][9] = new Rotamer(name, -177, 0, 65, 0, 180, 0, 85, 0);
                aminoAcidRotamerCache[n][10] = new Rotamer(name, -177, 0, 65, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][11] = new Rotamer(name, -177, 0, 180, 0, 65, 0, 85, 0);
                aminoAcidRotamerCache[n][12] = new Rotamer(name, -177, 0, 180, 0, 65, 0, -175, 0);
                aminoAcidRotamerCache[n][13] = new Rotamer(name, -177, 0, 180, 0, 65, 0, -105, 0);
                aminoAcidRotamerCache[n][14] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 85, 0);
                aminoAcidRotamerCache[n][15] = new Rotamer(name, -177, 0, 180, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][16] = new Rotamer(name, -177, 0, 180, 0, 180, 0, -85, 0);
                aminoAcidRotamerCache[n][17] = new Rotamer(name, -177, 0, 180, 0, -65, 0, 105, 0);
                aminoAcidRotamerCache[n][18] = new Rotamer(name, -177, 0, 180, 0, -65, 0, 175, 0);
                aminoAcidRotamerCache[n][19] = new Rotamer(name, -177, 0, 180, 0, -65, 0, -85, 0);
                aminoAcidRotamerCache[n][20] = new Rotamer(name, -67, 0, 180, 0, 65, 0, 85, 0);
                aminoAcidRotamerCache[n][21] = new Rotamer(name, -67, 0, 180, 0, 65, 0, -175, 0);
                aminoAcidRotamerCache[n][22] = new Rotamer(name, -67, 0, 180, 0, 65, 0, -105, 0);
                aminoAcidRotamerCache[n][23] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 85, 0);
                aminoAcidRotamerCache[n][24] = new Rotamer(name, -67, 0, 180, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][25] = new Rotamer(name, -67, 0, 180, 0, 180, 0, -85, 0);
                aminoAcidRotamerCache[n][26] = new Rotamer(name, -67, 0, 180, 0, -65, 0, 105, 0);
                aminoAcidRotamerCache[n][27] = new Rotamer(name, -67, 0, 180, 0, -65, 0, 175, 0);
                aminoAcidRotamerCache[n][28] = new Rotamer(name, -67, 0, -167, 0, -65, 0, -85, 0);
                aminoAcidRotamerCache[n][29] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 85, 0);
                aminoAcidRotamerCache[n][30] = new Rotamer(name, -62, 0, -68, 0, 180, 0, 180, 0);
                aminoAcidRotamerCache[n][31] = new Rotamer(name, -62, 0, -68, 0, 180, 0, -85, 0);
                aminoAcidRotamerCache[n][32] = new Rotamer(name, -62, 0, -68, 0, -65, 0, 175, 0);
                aminoAcidRotamerCache[n][33] = new Rotamer(name, -62, 0, -68, 0, -65, 0, -85, 0);
                break;
            default:
                // Handles GLY, ALA, CYX, PRO, ...
                break;
        }
        return aminoAcidRotamerCache[n];
    }

    /**
     * Returns the Rotamers for a specified nucleic acid type. Torsion angles
     * are listed from delta (i-1) to delta (i), along with standard deviations
     * calculated by Richardson et al.
     *
     * Citation: Richardson, J.S., et al., RNA backbone: Consensus all-angle
     * conformers and modular string nomenclature (an RNA Ontology Consortium
     * contribution). Rna-a Publication of the Rna Society, 2008. 14(3): p. 465-481.
     *
     * @param name Type of nucleic acid.
     * @return Rotamer cache (double[] of torsions).
     */
    private Rotamer[] getRichardsonRNARotamers(NucleicAcid3 name) {
        int n = name.ordinal();
        if (nucleicAcidRotamerCache[n] != null) {
            return nucleicAcidRotamerCache[n];
        }
        /*
         * Comments on rotamers can be found in Richardson et al, 2008.
         * Rotamers 0-45 are these rotamers in order.
         *
         * In the future, subsequent sets of 46 could be these rotamers with
         * rotations of either the base as a whole, or portions of the base
         * (such as the C6 amino group of adenine).  One suggestion is to use
         * %46 if only the backbone information is needed, or integer division
         * by 46 if the backbone information is unnecessary.
         *
         * Torsions are in order delta (i-1), epsilon (i-1), zeta (i-1), alpha,
         * beta, gamma, delta.  However, at this moment, only delta through alpha
         * (reverse order) are used to build the backbone, and delta (i-1) as
         * a binary function to determine what previous sugar pucker to expect.
         *
         * 1a, 1m, 1L, &a, 7a, 3a, 9a, 1g, 7d, 3d, 5d, 1e, 1c, 1f, 5j, 1b, 1{,
         * 3b, 1z, 5z, 7p, 1t, 5q, 1o, 7r, 2a, 4a, 0a, #a, 4g, 6g, 8d, 4d, 6d,
         * 2h, 4n, 0i, 6n, 6j, 2[, 4b, 0b, 4p, 6p, 4s, 2o
         */
        switch (name) {
            case ADE:
            case GUA:
            case CYT:
            case URI:
            case DAD:
            case DGU:
            case DCY:
            case DTY:
            case THY:
                nucleicAcidRotamerCache[n] = new Rotamer[46];
                nucleicAcidRotamerCache[n][0] = new Rotamer(name, 81, 4, -148, 10, -71, 7, -65, 8, 174, 8, 54, 6, 81, 3);
                nucleicAcidRotamerCache[n][1] = new Rotamer(name, 84, 5, -142, 16, -68, 15, -68, 16, -138, 12, 58, 10, 86, 7);
                nucleicAcidRotamerCache[n][2] = new Rotamer(name, 86, 4, -115, 6, -92, 13, -56, 8, 138, 4, 62, 10, 79, 5);
                nucleicAcidRotamerCache[n][3] = new Rotamer(name, 82, 5, -169, 7, -95, 6, -64, 9, -178, 10, 51, 7, 82, 5);
                nucleicAcidRotamerCache[n][4] = new Rotamer(name, 83, 4, -143, 23, -138, 14, -57, 9, 161, 15, 49, 6, 82, 3);
                nucleicAcidRotamerCache[n][5] = new Rotamer(name, 85, 4, -144, 24, 173, 14, -71, 12, 164, 16, 46, 7, 85, 6);//
                nucleicAcidRotamerCache[n][6] = new Rotamer(name, 83, 2, -150, 15, 121, 13, -71, 12, 157, 23, 49, 6, 81, 3);
                nucleicAcidRotamerCache[n][7] = new Rotamer(name, 81, 3, -141, 8, -69, 9, 167, 8, 160, 16, 51, 5, 85, 3);
                nucleicAcidRotamerCache[n][8] = new Rotamer(name, 84, 4, -121, 16, -103, 12, 70, 10, 170, 23, 53, 6, 85, 3);
                nucleicAcidRotamerCache[n][9] = new Rotamer(name, 85, 4, -116, 15, -156, 15, 66, 19, -179, 23, 55, 6, 86, 4);
                nucleicAcidRotamerCache[n][10] = new Rotamer(name, 80, 4, -158, 7, 63, 14, 68, 12, 143, 30, 50, 7, 83, 2);
                nucleicAcidRotamerCache[n][11] = new Rotamer(name, 81, 3, -159, 8, -79, 6, -111, 9, 83, 11, 168, 6, 86, 4);
                nucleicAcidRotamerCache[n][12] = new Rotamer(name, 80, 3, -163, 9, -69, 10, 153, 12, -166, 12, 179, 10, 84, 3);
                nucleicAcidRotamerCache[n][13] = new Rotamer(name, 81, 2, -157, 14, -66, 11, 172, 11, 139, 13, 176, 10, 84, 3);
                nucleicAcidRotamerCache[n][14] = new Rotamer(name, 87, 7, -136, 23, 80, 15, 67, 9, 109, 10, 176, 6, 84, 4);
                nucleicAcidRotamerCache[n][15] = new Rotamer(name, 84, 4, -145, 10, -71, 10, -60, 9, 177, 12, 58, 7, 145, 7);
                nucleicAcidRotamerCache[n][16] = new Rotamer(name, 83, 4, -140, 10, -71, 10, -63, 8, -138, 9, 54, 7, 144, 8);
                nucleicAcidRotamerCache[n][17] = new Rotamer(name, 85, 3, -134, 18, 168, 17, -67, 15, 178, 22, 49, 5, 148, 3);
                nucleicAcidRotamerCache[n][18] = new Rotamer(name, 83, 3, -154, 18, -82, 19, -164, 14, 162, 25, 51, 5, 145, 5);
                nucleicAcidRotamerCache[n][19] = new Rotamer(name, 83, 3, -154, 5, 53, 7, 164, 5, 148, 10, 50, 5, 148, 4);
                nucleicAcidRotamerCache[n][20] = new Rotamer(name, 84, 3, -123, 24, -140, 15, 68, 12, -160, 30, 54, 7, 146, 6);
                nucleicAcidRotamerCache[n][21] = new Rotamer(name, 81, 3, -161, 20, -71, 8, 180, 17, -165, 14, 178, 9, 147, 5);
                nucleicAcidRotamerCache[n][22] = new Rotamer(name, 82, 8, -155, 6, 69, 14, 63, 9, 115, 17, 176, 6, 146, 4);
                nucleicAcidRotamerCache[n][23] = new Rotamer(name, 84, 4, -143, 17, -73, 15, -63, 7, -135, 39, -66, 7, 151, 13);
                nucleicAcidRotamerCache[n][24] = new Rotamer(name, 85, 4, -127, 13, -112, 19, 63, 13, -178, 27, -64, 4, 150, 7);
                nucleicAcidRotamerCache[n][25] = new Rotamer(name, 145, 8, -100, 12, -71, 18, -72, 13, -167, 17, 53, 7, 84, 5);
                nucleicAcidRotamerCache[n][26] = new Rotamer(name, 146, 7, -100, 15, 170, 14, -62, 19, 170, 34, 51, 8, 84, 5);
                nucleicAcidRotamerCache[n][27] = new Rotamer(name, 149, 7, -137, 11, 139, 25, -75, 11, 158, 20, 48, 6, 84, 4);
                nucleicAcidRotamerCache[n][28] = new Rotamer(name, 148, 3, -168, 5, 146, 6, -71, 7, 151, 12, 42, 4, 85, 3);
                nucleicAcidRotamerCache[n][29] = new Rotamer(name, 148, 8, -103, 14, 165, 21, -155, 14, 165, 15, 49, 7, 83, 4);
                nucleicAcidRotamerCache[n][30] = new Rotamer(name, 145, 7, -97, 18, 80, 16, -156, 29, -170, 23, 58, 5, 85, 7);
                nucleicAcidRotamerCache[n][31] = new Rotamer(name, 149, 6, -89, 10, -119, 17, 62, 10, 176, 23, 54, 4, 87, 3);
                nucleicAcidRotamerCache[n][32] = new Rotamer(name, 150, 6, -110, 26, -172, 7, 80, 20, -162, 20, 61, 8, 89, 4);
                nucleicAcidRotamerCache[n][33] = new Rotamer(name, 147, 6, -119, 23, 89, 16, 59, 14, 161, 23, 52, 7, 83, 4);
                nucleicAcidRotamerCache[n][34] = new Rotamer(name, 148, 4, -99, 8, -70, 12, -64, 10, 177, 17, 176, 14, 87, 4);
                nucleicAcidRotamerCache[n][35] = new Rotamer(name, 144, 7, -133, 14, -156, 14, 74, 12, -143, 20, -166, 9, 81, 3);
                nucleicAcidRotamerCache[n][36] = new Rotamer(name, 149, 2, -85, 20, 100, 13, 81, 11, -112, 12, -178, 3, 83, 2);
                nucleicAcidRotamerCache[n][37] = new Rotamer(name, 150, 6, -92, 11, 85, 8, 64, 5, -169, 8, 177, 9, 86, 5);
                nucleicAcidRotamerCache[n][38] = new Rotamer(name, 142, 8, -116, 28, 66, 15, 72, 8, 122, 22, -178, 6, 84, 3);
                nucleicAcidRotamerCache[n][39] = new Rotamer(name, 146, 8, -101, 16, -69, 17, -68, 12, -150, 21, 54, 7, 148, 7);
                nucleicAcidRotamerCache[n][40] = new Rotamer(name, 145, 7, -115, 20, 163, 13, -66, 6, 172, 14, 46, 6, 146, 6);
                nucleicAcidRotamerCache[n][41] = new Rotamer(name, 148, 4, -112, 20, 112, 14, -85, 17, 165, 16, 57, 12, 146, 6);
                nucleicAcidRotamerCache[n][42] = new Rotamer(name, 150, 10, -100, 16, -146, 19, 72, 13, -152, 27, 57, 14, 148, 4);
                nucleicAcidRotamerCache[n][43] = new Rotamer(name, 146, 7, -102, 21, 90, 15, 68, 12, 173, 18, 56, 8, 148, 4);
                nucleicAcidRotamerCache[n][44] = new Rotamer(name, 150, 2, -112, 16, 170, 12, -82, 13, 84, 7, 176, 6, 148, 2);
                nucleicAcidRotamerCache[n][45] = new Rotamer(name, 147, 6, -104, 15, -64, 16, -73, 4, -165, 26, -66, 7, 150, 3);
                break;
            case MP1:
            case DP2:
            case TP3:
            case UNK:
            case M2MG:
            case H2U:
            case M2G:
            case OMC:
            case OMG:
            case PSU:
            case M5MC:
            case M7MG:
            case M5MU:
            case M1MA:
            case YYG:
            default:
                break;
        }
        return nucleicAcidRotamerCache[n];
    }

    /**
     * Measures the torsions in a list of Residues.
     *
     * @param residueList Residues to be measured.
     * @param print Verbosity flag.
     */
    public static void measureRotamers(List<Residue> residueList, boolean print) {
        double chi[] = new double[7];
        for (Residue residue : residueList) {
            chi[0] = chi[1] = chi[2] = chi[3] = chi[4] = chi[5] = chi[6] = 0.0;
            measureRotamer(residue, chi, print);
            switch (residue.getResidueType()) {
                case AA:
                    logger.info(String.format(" %c %s %8.3f %8.3f %8.3f %8.3f", residue.getChainID(), residue,
                            chi[0], chi[1], chi[2], chi[3]));
                    break;
                case NA:
                    logger.info(String.format(" %c %s %8.3f %8.3f %8.3f %8.3f %8.3f"
                            + " %8.3f %8.3f", residue.getChainID(), residue, chi[0],
                            chi[1], chi[2], chi[3], chi[4], chi[5], chi[6]));
                    break;
                default:
                    logger.info(String.format(
                            " Not recognized as a nucleic or amino acid residue"));
                    break;
            }
        }
    }

    /**
     * Measures the torsional angles of a residue's side chain.
     * @param residue
     * @param print
     * @return 
     */
    public static double[] measureRotamer(Residue residue, boolean print) {
        if (residue != null) {
            switch (residue.getResidueType()) {
                case AA:
                    double[] chi = new double[4];
                    try {
                        measureRotamer(residue, chi, print);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        String message = " Array passed to measureRotamer was not of sufficient size.";
                        logger.log(Level.WARNING, message, e);
                    }
                    return chi;
                case NA:
                    chi = new double[7];
                    try {
                        measureRotamer(residue, chi, print);
                    } catch (ArrayIndexOutOfBoundsException e) {
                        String message = " Array passed to measureRotamer was not of sufficient size.";
                        logger.log(Level.WARNING, message, e);
                    }
                    return chi;
                default:
                    return null;
            }
        }
        return null;
    }

    /**
     * Measures the torsion angles of a Residue.
     *
     * @param residue To be measured
     * @param chi Array to be filled with torsion values
     * @param print Verbosity flag
     */
    public static void measureRotamer(Residue residue, double chi[], boolean print) {
        if (residue == null) {
            return;
        }
        switch (residue.getResidueType()) {
            case AA:
                try {
                    measureAARotamer(residue, chi, print);
                } catch (ArrayIndexOutOfBoundsException e) {
                    String message = " Array passed to measureRotamer was not of sufficient size.";
                    logger.log(Level.WARNING, message, e);
                }
                break;

            case NA:
                try {
                    measureNARotamer(residue, chi, print);
                } catch (ArrayIndexOutOfBoundsException e) {
                    String message = "Array passed to measureRotamer was not of sufficient size.";
                    logger.log(Level.WARNING, message, e);
                }
                break;
            default:
                try {
                    measureUNKRotamer(residue, chi, print);
                } catch (ArrayIndexOutOfBoundsException e) {
                    String message = "Array passed to measureRotamer was not of sufficient size.";
                    logger.log(Level.WARNING, message, e);
                }
                break;
        }
    }
    
    public static void measureUNKRotamer(Residue residue, double[] chi, boolean print) {
        String resName = residue.getName().toUpperCase();
        if (nonstdRotCache.containsKey(resName)) {
            nonstdRotCache.get(resName).measureNonstdRot(residue, chi, print);
        } else {
            logger.warning(String.format(" Could not measure chi angles "
                    + "for residue %s", residue.toString()));
        }
    }

    /**
     * Measure the current torsions of a nucleic acid Residue, starting from the
     * 5'-most torsion.
     *
     * Chi[0]-chi[6] in order are: delta (i-1), epsilon (i-1), zeta (i-1), alpha
     * (i), beta (i), gamma (i) and delta (i) for residue i.
     *
     * @param residue Residue to be measured.
     * @param chi Array to be filled with torsion values.
     * @param print Verbosity flag.
     */
    public static void measureNARotamer(Residue residue, double chi[], boolean print) {
        NucleicAcid3 name = NucleicAcid3.valueOf(residue.getName());
        Residue prevResidue = residue.getPreviousResidue();
        Torsion torsion;

        Atom C5s = (Atom) residue.getAtomNode("C5\'");
        Atom C4s = (Atom) residue.getAtomNode("C4\'");
        Atom C3s = (Atom) residue.getAtomNode("C3\'");
        Atom O3s = (Atom) residue.getAtomNode("O3\'");
        Atom O5s = (Atom) residue.getAtomNode("O5\'");
        Atom P = (Atom) residue.getAtomNode("P");


        /*
         * Start by measuring delta (i-1) if available, working up to delta.  If
         * there is no prior residue, start measuring from the 5'-most torsion.
         */
        if (prevResidue == null) {
            switch (name) {
                case GUA:
                case ADE:
                case DGU:
                case DAD:
                case CYT:
                case URI:
                case THY:
                case DCY:
                case DTY:
                    /*
                     * If there is an HO5s, measure alpha based on HO5s.  Else,
                     * measure zeta (i-1) based on OP3 and alpha on P.
                     */
                    Atom HO5s = (Atom) residue.getAtomNode("HO5\'");
                    if (HO5s != null) {
                        torsion = HO5s.getTorsion(O5s, C5s, C4s);
                        chi[4] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    } else {
                        Atom OP3 = (Atom) residue.getAtomNode("OP3");
                        if (OP3 != null) {
                            torsion = OP3.getTorsion(P, O5s, C5s);
                            chi[3] = torsion.getValue();
                            if (print) {
                                logger.info(torsion.toString());
                            }
                        }

                        torsion = P.getTorsion(O5s, C5s, C4s);
                        chi[4] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    break;
                default:
                    break;
            }
        } else {
            switch (name) {
                case GUA:
                case ADE:
                case DGU:
                case DAD:
                case CYT:
                case URI:
                case THY:
                case DCY:
                case DTY:
                    Atom O3sPrev = (Atom) prevResidue.getAtomNode("O3\'");
                    Atom C3sPrev = (Atom) prevResidue.getAtomNode("C3\'");
                    Atom C4sPrev = (Atom) prevResidue.getAtomNode("C4\'");
                    Atom C5sPrev = (Atom) prevResidue.getAtomNode("C5\'");

                    torsion = C5sPrev.getTorsion(C4sPrev, C3sPrev, O3sPrev);
                    chi[0] = torsion.getValue();
                    if (print) {
                        logger.info(torsion.toString());
                    }

                    torsion = C4sPrev.getTorsion(C3sPrev, O3sPrev, P);
                    chi[1] = torsion.getValue();
                    if (print) {
                        logger.info(torsion.toString());
                    }

                    torsion = C3sPrev.getTorsion(O3sPrev, P, O5s);
                    chi[2] = torsion.getValue();
                    if (print) {
                        logger.info(torsion.toString());
                    }

                    torsion = O3sPrev.getTorsion(P, O5s, C5s);
                    chi[3] = torsion.getValue();
                    if (print) {
                        logger.info(torsion.toString());
                    }

                    torsion = P.getTorsion(O5s, C5s, C4s);
                    chi[4] = torsion.getValue();
                    if (print) {
                        logger.info(torsion.toString());
                    }
                    break;
                default:
                    break;
            }

        }
        /*
         * Measure torsions common to all nucleic acids (gamma, delta).
         */
        torsion = O5s.getTorsion(C5s, C4s, C3s);
        chi[5] = torsion.getValue();
        if (print) {
            logger.info(torsion.toString());
        }

        torsion = C5s.getTorsion(C4s, C3s, O3s);
        chi[6] = torsion.getValue();
        if (print) {
            logger.info(torsion.toString());
        }
    }

    /**
     * Measures the delta torsion (sugar pucker) of a nucleic acid Residue.
     *
     * @param residue To be measured
     * @return Delta torsion (sugar pucker angle).
     */
    public static double measureDelta(Residue residue) {
        Atom C5s = (Atom) residue.getAtomNode("C5\'");
        Atom C4s = (Atom) residue.getAtomNode("C4\'");
        Atom C3s = (Atom) residue.getAtomNode("C3\'");
        Atom O3s = (Atom) residue.getAtomNode("O3\'");
        Torsion torsion = O3s.getTorsion(C3s, C4s, C5s);
        return torsion.getValue();
    }

    /**
     * Measures the torsions of an amino acid Residue's current configuration.
     *
     * @param residue To be measured.
     * @param chi Array to be filled with torsion values.
     * @param print Verbosity flag.
     */
    public static void measureAARotamer(Residue residue, double chi[], boolean print) {
        if (residue instanceof MultiResidue) {
            residue = ((MultiResidue) residue).getActive();
        }
        AminoAcid3 name = residue.getAminoAcid3();
        /*try {
            name = AminoAcid3.valueOf(residue.getName());
        } catch (IllegalArgumentException ex) {
            logger.info(String.format("(IAE) valueOf(%s)", residue.getName()));
            throw ex;
        }*/
        switch (name) {
            case VAL: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG1)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                        break;
                    }
                }
                break;
            }
            case LEU: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case ILE: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG1)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG1, CD1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case SER: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom OG = (Atom) residue.getAtomNode("OG");
                Atom HG = (Atom) residue.getAtomNode("HG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, OG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, OG, HG)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case THR: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom OG1 = (Atom) residue.getAtomNode("OG1");
                Atom HG1 = (Atom) residue.getAtomNode("HG1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, OG1)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, OG1, HG1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case CYX: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom SG = (Atom) residue.getAtomNode("SG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, SG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                        break;
                    }
                }
                break;
            }
            case CYD: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom SG = (Atom) residue.getAtomNode("SG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, SG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                        break;
                    }
                }
                break;
            }
            case PHE: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CG = (Atom) residue.getAtomNode("CG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                        break;
                    }
                }
                break;
            }
            case PRO: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CG = (Atom) residue.getAtomNode("CG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case TYR: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom OH = (Atom) residue.getAtomNode("OH");
                Atom HH = (Atom) residue.getAtomNode("HH");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CE2, CZ, OH, HH)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case TYD: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CG = (Atom) residue.getAtomNode("CG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case TRP: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CG = (Atom) residue.getAtomNode("CG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case HIS: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, ND1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case HID: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, ND1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case HIE: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, ND1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case ASP: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                        break;
                    }
                }
                break;
            }
            case ASH: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, OD1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case ASN: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, OD1)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case GLU: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CB, CG, CD, OE1)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case GLH: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CB, CG, CD, OE1)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case GLN: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CB, CG, CD, OE1)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case MET: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom SD = (Atom) residue.getAtomNode("SD");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, SD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CB, CG, SD, CE)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case LYS: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom NZ = (Atom) residue.getAtomNode("NZ");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CB, CG, CD, CE)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CG, CD, CE, NZ)) {
                        chi[3] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case LYD: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom NZ = (Atom) residue.getAtomNode("NZ");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CB, CG, CD, CE)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CG, CD, CE, NZ)) {
                        chi[3] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }

                }
                break;
            }
            case ARG: {
                ArrayList<ROLS> torsions = residue.getTorsionList();
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom NE = (Atom) residue.getAtomNode("NE");
                for (ROLS rols : torsions) {
                    Torsion torsion = (Torsion) rols;
                    if (torsion.compare(N, CA, CB, CG)) {
                        chi[0] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CA, CB, CG, CD)) {
                        chi[1] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CB, CG, CD, NE)) {
                        chi[2] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                    if (torsion.compare(CG, CD, NE, CZ)) {
                        chi[3] = torsion.getValue();
                        if (print) {
                            logger.info(torsion.toString());
                        }
                    }
                }
                break;
            }
            case UNK:
                chi = new double[7];
                String resName = residue.getName().toUpperCase();
                if (nonstdRotCache.containsKey(resName)) {
                    nonstdRotCache.get(resName).measureNonstdRot(residue, chi, print);
                    //nonstdRotCache.get(resName).applyNonstdRotamer(residue, rotamer);
                } else {
                    throw new IllegalArgumentException(String.format("(IAE) valueOf(%s)", residue.getName()));
                }
                break;
            default: {
            }

        }
    }

    /**
     * Applies a Rotamer to a Residue by calling applyAARotamer or
     * applyNARotamer.
     *
     * @param residue the Residue whose side-chain will be moved.
     * @param rotamer the Rotamer defining the move.
     */
    public static void applyRotamer(Residue residue, Rotamer rotamer) {
        applyRotamer(residue, rotamer, false);
    }

    /**
     * Version of applyRotamer which allows for chain context-independent
     * drawing of nucleic acid Rotamers. Solely used in saveRotamers at this
     * point, although it may be useful for debugging.
     *
     * @param residue the Residue to be moved.
     * @param rotamer Rotamer to be applied.
     * @param independent Whether to draw Rotamer independent of chain context.
     */
    public static void applyRotamer(Residue residue, Rotamer rotamer, boolean independent) {
        residue.setRotamer(rotamer);
        if (rotamer.isState) {
            applyState(residue, rotamer);
        } else {
            switch (residue.getResidueType()) {
                case AA:
                    applyAARotamer(residue, rotamer);
                    break;
                case NA:
                    applyNARotamer(residue, rotamer, independent);
                    break;
                default:
                    break;
            }
        }

    }

    /**
     * Applies a coordinates-based Rotamer (defined by Cartesian coordinates
     * instead of by a set of torsion angles); intended for use with original
     * coordinates Rotamers and possibly other future cases.
     *
     * @param residue Residue to apply Rotamer for
     * @param rotamer Coordinates-based Rotamer
     */
    private static void applyState(Residue residue, Rotamer rotamer) {
        if (rotamer.isState) {
            residue.revertState(rotamer.originalState);
        } else {
            logger.warning(String.format(" Attempting to apply a ResidueState for "
                    + "a torsion-based rotamer %s for residue %s", rotamer, residue));
        }
    }

    /**
     * Obtains an idealized angle using the force field. Note: can be misleading!
     * Often, idealized angles are not equilibrium angles; for example, AMOEBA 2018
     * has phenylalanine ring interior angles that sum to > 720 degrees. This does not
     * mean that an idealized phenylalanine ring violates geometry, it means that the
     * ring is always under at least a little strain.
     *
     * @param a1 An Atom.
     * @param a2 Another Atom.
     * @param a3 Another Atom.
     * @return Force field-defined a1-a2-a3 angle in degrees.
     */
    private static double angleFromForceField(Atom a1, Atom a2, Atom a3) {
        Angle a = a1.getAngle(a2, a3);
        AngleType at = a.getAngleType();
        return at.angle[a.nh];
    }

    /**
     * Obtains an idealized angle using a lookup map of stored, idealized geometries,
     * with a fallback to using the force field. This is intended for use with ring
     * systems, where ring constraints mean that optimum-energy rings may not have
     * the values defined by the force field.
     *
     * Current lookup map is only for PHE, TYR, TYD, HIS, HID, HIE, and TRP, using
     * values obtained from a tight bonded-terms-only optimization under AMOEBA BIO 2018.
     *
     * @param resName Name of the Residue containing a1-a3.
     * @param a1 An atom.
     * @param a2 Another Atom.
     * @param a3 Another Atom.
     * @return Stored idealized a1-a2-a3 angle in degrees.
     */
    private static double idealGeometryAngle(AminoAcid3 resName, Atom a1, Atom a2, Atom a3) {
        StringBuilder sb = new StringBuilder(a1.getName());
        sb.append("-").append(a2.getName()).append("-").append(a3.getName());
        Map<String, Double> resMap = idealAngleGeometries.get(resName);
        String atomString = sb.toString();

        if (resMap.containsKey(atomString)) {
            return resMap.get(atomString);
        } else {
            sb = new StringBuilder(a3.getName());
            sb.append("-").append(a2.getName()).append("-").append(a1.getName());
            atomString = sb.toString();
            if (resMap.containsKey(atomString)) {
                return resMap.get(atomString);
            } else {
                logger.fine(String.format(" Could not find an ideal-geometry angle for %s %s-%s-%s", resName, a1, a2, a3));
                return angleFromForceField(a1, a2, a3);
            }
        }
    }

    /**
     * Gets an Angle, using the default as set by property for whether to use idealized ring
     * geometry (default), or based on force field.
     *
     * @param resName AminoAcid3 for a1-a3.
     * @param a1 An Atom.
     * @param a2 Another Atom.
     * @param a3 A third Atom.
     * @return a1-a2-a3 angle for internal geometry.
     */
    private static double getAngle(AminoAcid3 resName, Atom a1, Atom a2, Atom a3) {
        if (useIdealRingGeometries) {
            return idealGeometryAngle(resName, a1, a2, a3);
        } else {
            return angleFromForceField(a1, a2, a3);
        }
    }

    /**
     * Draws CZ of Phe/Tyr/Tyd twice (from each branch of the ring), the cuts it down the middle.
     * @param resName Residue containing CZ.
     * @param CZ CZ to be placed.
     * @param CG
     * @param CE1
     * @param CD1
     * @param CE2
     * @param CD2
     * @return Mean coordinates for CZ based on internal geometry.
     */
    private static double[] drawCZ(AminoAcid3 resName, Atom CZ, Atom CG, Atom CE1, Atom CD1, Atom CE2, Atom CD2) {
        double bondLen = CZ.getBond(CE1).bondType.distance;
        double ang = getAngle(resName, CZ, CE1, CD1);
        double[] xCG = new double[3];
        xCG = CG.getXYZ(xCG);

        double[] xCE = new double[3];
        xCE = CE1.getXYZ(xCE);
        double[] xCD = new double[3];
        xCD = CD1.getXYZ(xCD);
        double[] xyz1 = determineIntxyz(xCE, bondLen, xCD, ang, xCG, 0.0, 0);

        xCE = CE2.getXYZ(xCE);
        xCD = CD2.getXYZ(xCD);
        double[] xyz2 = determineIntxyz(xCE, bondLen, xCD, ang, xCG, 0, 0);

        for (int i = 0; i < 3; i++) {
            xyz1[i] += xyz2[i];
            xyz1[i] *= 0.5;
        }
        return xyz1;
    }

    /**
     * Moves CZ of Phe/Tyr/Tyd to a mean position determined by both branches of the ring.
     * @param resName Residue containing CZ.
     * @param CZ CZ to be placed.
     * @param CG
     * @param CE1
     * @param CD1
     * @param CE2
     * @param CD2
     */
    private static void applyCZ(AminoAcid3 resName, Atom CZ, Atom CG, Atom CE1, Atom CD1, Atom CE2, Atom CD2) {
        CZ.moveTo(drawCZ(resName, CZ, CG, CE1, CD1, CE2, CD2));
    }

    /**
     * Applies an amino acid Rotamer.
     *
     * @param residue Residue
     * @param rotamer Rotamer to be applied to Residue
     */
    private static void applyAARotamer(Residue residue, Rotamer rotamer) {
        if (residue == null || rotamer == null) {
            return;
        }
        AminoAcid3 name;
        if (residue instanceof MultiResidue) {
            name = rotamer.name;
            if (!((MultiResidue) residue).setActiveResidue(name)) {
                logger.warning(String.format(" Could not set residue %s for multi-residue %s", name, residue));
            }
        } else {
            name = residue.getAminoAcid3();
        }
        switch (name) {
            case VAL: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                Atom CG2 = (Atom) residue.getAtomNode("CG2");
                Atom HB = (Atom) residue.getAtomNode("HB");
                Atom HG11 = (Atom) residue.getAtomNode("HG11");
                Atom HG12 = (Atom) residue.getAtomNode("HG12");
                Atom HG13 = (Atom) residue.getAtomNode("HG13");
                Atom HG21 = (Atom) residue.getAtomNode("HG21");
                Atom HG22 = (Atom) residue.getAtomNode("HG22");
                Atom HG23 = (Atom) residue.getAtomNode("HG23");
                Bond CG_CB = CB.getBond(CG1);
                Bond HB_CB = CB.getBond(HB);
                Bond HG_CG = HG11.getBond(CG1);
                double dCG_CB = CG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                Angle CG_CB_CA = CG1.getAngle(CB, CA);
                Angle HB_CB_CA = HB.getAngle(CB, CA);
                Angle HG_CG_CB = HG11.getAngle(CG1, CB);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                intxyz(CG1, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CG2, CB, dCG_CB, CA, dCG_CB_CA, CG1, 109.5, -1);
                intxyz(HB, CB, dHB_CB, CA, dHB_CB_CA, CG1, 109.4, 1);
                intxyz(HG11, CG1, dHG_CG, CB, dHG_CG_CB, CA, 180.0, 0);
                intxyz(HG12, CG1, dHG_CG, CB, dHG_CG_CB, HG11, 109.4, 1);
                intxyz(HG13, CG1, dHG_CG, CB, dHG_CG_CB, HG11, 109.4, -1);
                intxyz(HG21, CG2, dHG_CG, CB, dHG_CG_CB, CA, 180.0, 0);
                intxyz(HG22, CG2, dHG_CG, CB, dHG_CG_CB, HG21, 109.4, 1);
                intxyz(HG23, CG2, dHG_CG, CB, dHG_CG_CB, HG21, 109.4, -1);
                break;
            }
            case LEU: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG = (Atom) residue.getAtomNode("HG");
                Atom HD11 = (Atom) residue.getAtomNode("HD11");
                Atom HD12 = (Atom) residue.getAtomNode("HD12");
                Atom HD13 = (Atom) residue.getAtomNode("HD13");
                Atom HD21 = (Atom) residue.getAtomNode("HD21");
                Atom HD22 = (Atom) residue.getAtomNode("HD22");
                Atom HD23 = (Atom) residue.getAtomNode("HD23");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG.getBond(CG);
                Bond HD_CD = HD11.getBond(CD1);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD1.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG.getAngle(CG, CB);
                Angle HD_CD_CG = HD11.getAngle(CD1, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CD1, 109.5, -1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG, CG, dHG_CG, CB, dHG_CG_CB, CD1, 109.4, 1);
                intxyz(HD11, CD1, dHD_CD, CG, dHD_CD_CG, CB, 180.0, 0);
                intxyz(HD12, CD1, dHD_CD, CG, dHD_CD_CG, HD11, 109.4, 1);
                intxyz(HD13, CD1, dHD_CD, CG, dHD_CD_CG, HD11, 109.4, -1);
                intxyz(HD21, CD2, dHD_CD, CG, dHD_CD_CG, CB, 180.0, 0);
                intxyz(HD22, CD2, dHD_CD, CG, dHD_CD_CG, HD21, 109.4, 1);
                intxyz(HD23, CD2, dHD_CD, CG, dHD_CD_CG, HD21, 109.4, -1);
                break;
            }
            case ILE: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG1 = (Atom) residue.getAtomNode("CG1");
                Atom CG2 = (Atom) residue.getAtomNode("CG2");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom HB = (Atom) residue.getAtomNode("HB");
                Atom HG12 = (Atom) residue.getAtomNode("HG12");
                Atom HG13 = (Atom) residue.getAtomNode("HG13");
                Atom HG21 = (Atom) residue.getAtomNode("HG21");
                Atom HG22 = (Atom) residue.getAtomNode("HG22");
                Atom HG23 = (Atom) residue.getAtomNode("HG23");
                Atom HD11 = (Atom) residue.getAtomNode("HD11");
                Atom HD12 = (Atom) residue.getAtomNode("HD12");
                Atom HD13 = (Atom) residue.getAtomNode("HD13");
                Bond CG1_CB = CG1.getBond(CB);
                Bond CG2_CB = CG2.getBond(CB);
                Bond CD1_CG1 = CD1.getBond(CG1);
                Bond HB_CB = HB.getBond(CB);
                Bond HG1_CG = HG12.getBond(CG1);
                Bond HG2_CG = HG22.getBond(CG2);
                Bond HD_CD = HD12.getBond(CD1);
                double dCG1_CB = CG1_CB.bondType.distance;
                double dCG2_CB = CG2_CB.bondType.distance;
                double dCD1_CG1 = CD1_CG1.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG1_CG = HG1_CG.bondType.distance;
                double dHG2_CG = HG2_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                Angle CG1_CB_CA = CG1.getAngle(CB, CA);
                Angle CG2_CB_CA = CG2.getAngle(CB, CA);
                Angle CD1_CG1_CB = CD1.getAngle(CG1, CB);
                Angle HB_CB_CA = HB.getAngle(CB, CA);
                Angle HG1_CG_CB = HG12.getAngle(CG1, CB);
                Angle HG2_CG_CB = HG21.getAngle(CG2, CB);
                Angle HD_CD1_CG1 = HD11.getAngle(CD1, CG1);
                double dCG1_CB_CA = CG1_CB_CA.angleType.angle[CG1_CB_CA.nh];
                double dCG2_CB_CA = CG2_CB_CA.angleType.angle[CG2_CB_CA.nh];
                double dCD1_CG1_CB = CD1_CG1_CB.angleType.angle[CD1_CG1_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG1_CG_CB = HG1_CG_CB.angleType.angle[HG1_CG_CB.nh];
                double dHG2_CG_CB = HG2_CG_CB.angleType.angle[HG2_CG_CB.nh];
                double dHD_CD1_CG1 = HD_CD1_CG1.angleType.angle[HD_CD1_CG1.nh];
                intxyz(CG1, CB, dCG1_CB, CA, dCG1_CB_CA, N, rotamer.chi1, 0);
                intxyz(CG2, CB, dCG2_CB, CA, dCG2_CB_CA, CG1, 109.5, 1);
                intxyz(CD1, CG1, dCD1_CG1, CB, dCD1_CG1_CB, CA, rotamer.chi2, 0);
                intxyz(HB, CB, dHB_CB, CA, dHB_CB_CA, CG2, 109.4, 1);
                intxyz(HG12, CG1, dHG1_CG, CB, dHG1_CG_CB, CD1, 109.4, 1);
                intxyz(HG13, CG1, dHG1_CG, CB, dHG1_CG_CB, CD1, 109.4, -1);
                intxyz(HG21, CG2, dHG2_CG, CB, dHG2_CG_CB, CG1, 180.0, 0);
                intxyz(HG22, CG2, dHG2_CG, CB, dHG2_CG_CB, HG21, 109.0, 1);
                intxyz(HG23, CG2, dHG2_CG, CB, dHG2_CG_CB, HG21, 109.0, -1);
                intxyz(HD11, CD1, dHD_CD, CG1, dHD_CD1_CG1, CB, 180.0, 0);
                intxyz(HD12, CD1, dHD_CD, CG1, dHD_CD1_CG1, HD11, 109.0, 1);
                intxyz(HD13, CD1, dHD_CD, CG1, dHD_CD1_CG1, HD11, 109.0, -1);
                break;
            }
            case SER: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom OG = (Atom) residue.getAtomNode("OG");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG = (Atom) residue.getAtomNode("HG");
                Bond OG_CB = OG.getBond(CB);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_OG = HG.getBond(OG);
                double dOG_CB = OG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_OG = HG_OG.bondType.distance;
                Angle OG_CB_CA = OG.getAngle(CB, CA);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_OG_CB = HG.getAngle(OG, CB);
                double dOG_CB_CA = OG_CB_CA.angleType.angle[OG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_OG_CB = HG_OG_CB.angleType.angle[HG_OG_CB.nh];
                intxyz(OG, CB, dOG_CB, CA, dOG_CB_CA, N, rotamer.chi1, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, OG, 106.7, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, OG, 106.7, -1);
                intxyz(HG, OG, dHG_OG, CB, dHG_OG_CB, CA, 180.0, 0);
                if (rotamer.length == 2) {
                    intxyz(HG, OG, dHG_OG, CB, dHG_OG_CB, CA, rotamer.chi2, 0);
                } else {
                    intxyz(HG, OG, dHG_OG, CB, dHG_OG_CB, CA, 180.0, 0);
                }
                break;
            }
            case THR: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom OG1 = (Atom) residue.getAtomNode("OG1");
                Atom CG2 = (Atom) residue.getAtomNode("CG2");
                Atom HB = (Atom) residue.getAtomNode("HB");
                Atom HG1 = (Atom) residue.getAtomNode("HG1");
                Atom HG21 = (Atom) residue.getAtomNode("HG21");
                Atom HG22 = (Atom) residue.getAtomNode("HG22");
                Atom HG23 = (Atom) residue.getAtomNode("HG23");
                Bond OG1_CB = OG1.getBond(CB);
                Bond CG2_CB = CG2.getBond(CB);
                Bond HB_CB = HB.getBond(CB);
                Bond HG1_OG1 = HG1.getBond(OG1);
                Bond HG2_CG2 = HG21.getBond(CG2);
                double dOG1_CB = OG1_CB.bondType.distance;
                double dCG2_CB = CG2_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG1_OG1 = HG1_OG1.bondType.distance;
                double dHG2_CG2 = HG2_CG2.bondType.distance;
                Angle OG1_CB_CA = OG1.getAngle(CB, CA);
                Angle CG2_CB_CA = CG2.getAngle(CB, CA);
                Angle HB_CB_CA = HB.getAngle(CB, CA);
                Angle HG1_OG1_CB = HG1.getAngle(OG1, CB);
                Angle HG2_CG2_CB = HG21.getAngle(CG2, CB);
                double dOG1_CB_CA = OG1_CB_CA.angleType.angle[OG1_CB_CA.nh];
                double dCG2_CB_CA = CG2_CB_CA.angleType.angle[CG2_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG1_OG1_CB = HG1_OG1_CB.angleType.angle[HG1_OG1_CB.nh];
                double dHG2_CG2_CB = HG2_CG2_CB.angleType.angle[HG2_CG2_CB.nh];
                intxyz(OG1, CB, dOG1_CB, CA, dOG1_CB_CA, N, rotamer.chi1, 0);
                intxyz(CG2, CB, dCG2_CB, CA, dCG2_CB_CA, OG1, 107.7, 1);
                intxyz(HB, CB, dHB_CB, CA, dHB_CB_CA, OG1, 106.7, -1);
                intxyz(HG1, OG1, dHG1_OG1, CB, dHG1_OG1_CB, CA, 180.0, 0);
                if (rotamer.length == 2) {
                    intxyz(HG1, OG1, dHG1_OG1, CB, dHG1_OG1_CB, CA, rotamer.chi2, 0);
                } else {
                    intxyz(HG1, OG1, dHG1_OG1, CB, dHG1_OG1_CB, CA, 180, 0);
                }
                intxyz(HG21, CG2, dHG2_CG2, CB, dHG2_CG2_CB, CA, 180.0, 0);
                intxyz(HG22, CG2, dHG2_CG2, CB, dHG2_CG2_CB, HG21, 109.0, 1);
                intxyz(HG23, CG2, dHG2_CG2, CB, dHG2_CG2_CB, HG21, 109.0, -1);
                break;
            }
            case CYS:
            case CYX: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom SG = (Atom) residue.getAtomNode("SG");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG = (Atom) residue.getAtomNode("HG");
                if (CA == null || CB == null || N == null || SG == null || HB2 == null || HB3 == null || HG == null) {
                    break;
                }
                Bond SG_CB = SG.getBond(CB);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_SG = HG.getBond(SG);
                double dSG_CB = SG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_SG = HG_SG.bondType.distance;
                Angle SG_CB_CA = SG.getAngle(CB, CA);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_SG_CB = HG.getAngle(SG, CB);
                double dSG_CB_CA = SG_CB_CA.angleType.angle[SG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_SG_CB = HG_SG_CB.angleType.angle[HG_SG_CB.nh];
                intxyz(SG, CB, dSG_CB, CA, dSG_CB_CA, N, rotamer.chi1, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, -1);
                intxyz(HG, SG, dHG_SG, CB, dHG_SG_CB, CA, 180.0, 0);
                break;
            }
            case CYD: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom SG = (Atom) residue.getAtomNode("SG");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Bond SG_CB = SG.getBond(CB);
                Bond HB_CB = HB2.getBond(CB);
                double dSG_CB = SG_CB.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                Angle SG_CB_CA = SG.getAngle(CB, CA);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                double dSG_CB_CA = SG_CB_CA.angleType.angle[SG_CB_CA.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                intxyz(SG, CB, dSG_CB, CA, dSG_CB_CA, N, rotamer.chi1, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, SG, 112.0, -1);
                break;
            }
            case PHE: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HZ = (Atom) residue.getAtomNode("HZ");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond CE_CD = CE1.getBond(CD1);
                Bond HB_CB = HB2.getBond(CB);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;

                double dHD_CD = HD1.getBond(CD1).bondType.distance;
                double dHE_CE = HE1.getBond(CE1).bondType.distance;
                double dHZ_CZ = HZ.getBond(CZ).bondType.distance;

                double dCG_CB_CA = getAngle(name, CG, CB, CA);
                double dCD_CG_CB = getAngle(name, CD1, CG, CB);
                double dCE_CD_CG = getAngle(name, CE1, CD1, CG);
                double dHB_CB_CA = getAngle(name, HB2, CB, CA);

                double dHD_CD_CG = getAngle(name, HD1, CD1, CG);
                double dHD_CD_CE = getAngle(name, HD1, CD1, CE1);
                double dHE_CE_CD = getAngle(name, HE1, CE1, CD1);
                double dHE_CE_CZ = getAngle(name, HE1, CE1, CZ);
                double dHZ_CZ_CE = getAngle(name, HZ, CZ, CE1);

                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2 + 180, 0);
                intxyz(CE1, CD1, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CE2, CD2, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                applyCZ(name, CZ, CG, CE1, CD1, CE2, CD2);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);

                intxyz(HD1, CD1, dHD_CD, CG, dHD_CD_CG, CE1, dHD_CD_CE, 3);
                intxyz(HD2, CD2, dHD_CD, CG, dHD_CD_CG, CE2, dHD_CD_CE, 3);
                intxyz(HE1, CE1, dHE_CE, CD1, dHE_CE_CD, CZ, dHE_CE_CZ, 3);
                intxyz(HE2, CE2, dHE_CE, CD2, dHE_CE_CD, CZ, dHE_CE_CZ, 3);
                intxyz(HZ, CZ, dHZ_CZ, CE1, dHZ_CZ_CE, CE2, dHZ_CZ_CE, 3);
                break;
            }
            case PRO: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, N, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, N, 109.4, -1);
                break;
            }
            case TYR: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom OH = (Atom) residue.getAtomNode("OH");
                Atom HH = (Atom) residue.getAtomNode("HH");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond CE_CD = CE1.getBond(CD1);
                Bond HB_CB = HB2.getBond(CB);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;

                double dHD_CD = HD1.getBond(CD1).bondType.distance;
                double dHE_CE = HE1.getBond(CE1).bondType.distance;
                double dOH_CZ = OH.getBond(CZ).bondType.distance;
                double dHH_OH = HH.getBond(OH).bondType.distance;

                double dCG_CB_CA = getAngle(name, CG, CB, CA);
                double dCD_CG_CB = getAngle(name, CD1, CG, CB);
                double dCE_CD_CG = getAngle(name, CE1, CD1, CG);
                double dHB_CB_CA = getAngle(name, HB2, CB, CA);

                double dHD_CD_CG = getAngle(name, HD1, CD1, CG);
                double dHD_CD_CE = getAngle(name, HD1, CD1, CE1);
                double dHE_CE_CD = getAngle(name, HE1, CE1, CD1);
                double dHE_CE_CZ = getAngle(name, HE1, CE1, CZ);
                double dOH_CZ_CE = getAngle(name, OH, CZ, CE1);
                double dHH_OH_CZ = getAngle(name, HH, OH, CZ);

                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2 + 180, 0);
                intxyz(CE1, CD1, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CE2, CD2, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                applyCZ(name, CZ, CG, CE1, CD1, CE2, CD2);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);

                intxyz(HD1, CD1, dHD_CD, CG, dHD_CD_CG, CE1, dHD_CD_CE, 3);
                intxyz(HD2, CD2, dHD_CD, CG, dHD_CD_CG, CE2, dHD_CD_CE, 3);
                intxyz(HE1, CE1, dHE_CE, CD1, dHE_CE_CD, CZ, dHE_CE_CZ, 3);
                intxyz(HE2, CE2, dHE_CE, CD2, dHE_CE_CD, CZ, dHE_CE_CZ, 3);
                intxyz(OH, CZ, dOH_CZ, CE1, dOH_CZ_CE, CE2, dOH_CZ_CE, 3);

                if (rotamer.length == 3) {
                    intxyz(HH, OH, dHH_OH, CZ, dHH_OH_CZ, CE2, rotamer.chi3, 0);
                } else {
                    intxyz(HH, OH, dHH_OH, CZ, dHH_OH_CZ, CE2, 0.0, 0);
                }
                break;
            }
            case TYD: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom OH = (Atom) residue.getAtomNode("OH");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD1.getBond(CG);
                Bond CE_CD = CE1.getBond(CD1);
                Bond HB_CB = HB2.getBond(CB);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;

                double dHD_CD = HD1.getBond(CD1).bondType.distance;
                double dHE_CE = HE1.getBond(CE1).bondType.distance;
                double dOH_CZ = OH.getBond(CZ).bondType.distance;

                double dCG_CB_CA = getAngle(name, CG, CB, CA);
                double dCD_CG_CB = getAngle(name, CD1, CG, CB);
                double dCE_CD_CG = getAngle(name, CE1, CD1, CG);
                double dHB_CB_CA = getAngle(name, HB2, CB, CA);

                double dHD_CD_CG = getAngle(name, HD1, CD1, CG);
                double dHD_CD_CE = getAngle(name, HD1, CD1, CE1);
                double dHE_CE_CD = getAngle(name, HE1, CE1, CD1);
                double dHE_CE_CZ = getAngle(name, HE1, CE1, CZ);
                double dOH_CZ_CE = getAngle(name, OH, CZ, CE1);

                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2 + 180, 0);
                intxyz(CE1, CD1, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                intxyz(CE2, CD2, dCE_CD, CG, dCE_CD_CG, CB, 180, 0);
                applyCZ(name, CZ, CG, CE1, CD1, CE2, CD2);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);

                intxyz(HD1, CD1, dHD_CD, CG, dHD_CD_CG, CE1, dHD_CD_CE, 3);
                intxyz(HD2, CD2, dHD_CD, CG, dHD_CD_CG, CE2, dHD_CD_CE, 3);
                intxyz(HE1, CE1, dHE_CE, CD1, dHE_CE_CD, CZ, dHE_CE_CZ, 3);
                intxyz(HE2, CE2, dHE_CE, CD2, dHE_CE_CD, CZ, dHE_CE_CZ, 3);
                intxyz(OH, CZ, dOH_CZ, CE1, dOH_CZ_CE, CE2, dOH_CZ_CE, 3);
                break;
            }
            case TRP: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD1 = (Atom) residue.getAtomNode("CD1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom NE1 = (Atom) residue.getAtomNode("NE1");
                Atom CE2 = (Atom) residue.getAtomNode("CE2");
                Atom CE3 = (Atom) residue.getAtomNode("CE3");
                Atom CZ2 = (Atom) residue.getAtomNode("CZ2");
                Atom CZ3 = (Atom) residue.getAtomNode("CZ3");
                Atom CH2 = (Atom) residue.getAtomNode("CH2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Atom HZ2 = (Atom) residue.getAtomNode("HZ2");
                Atom HZ3 = (Atom) residue.getAtomNode("HZ3");
                Atom HH2 = (Atom) residue.getAtomNode("HH2");
                Bond CG_CB = CG.getBond(CB);
                Bond CD1_CG = CD1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond NE1_CD1 = NE1.getBond(CD1);
                Bond CE2_NE1 = CE2.getBond(NE1);
                Bond CE3_CD2 = CE3.getBond(CD2);
                Bond CZ2_CE2 = CZ2.getBond(CE2);
                Bond CZ3_CE3 = CZ3.getBond(CE3);
                Bond CH2_CZ2 = CH2.getBond(CZ2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD1_CD1 = HD1.getBond(CD1);
                Bond HE1_NE1 = HE1.getBond(NE1);
                Bond HE3_CE3 = HE3.getBond(CE3);
                Bond HZ2_CZ2 = HZ2.getBond(CZ2);
                Bond HZ3_CZ3 = HZ3.getBond(CZ3);
                Bond HH2_CH2 = HH2.getBond(CH2);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD1_CG = CD1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dNE1_CD1 = NE1_CD1.bondType.distance;
                double dCE2_NE1 = CE2_NE1.bondType.distance;
                double dCE3_CD2 = CE3_CD2.bondType.distance;
                double dCZ2_CE2 = CZ2_CE2.bondType.distance;
                double dCZ3_CE3 = CZ3_CE3.bondType.distance;
                double dCH2_CZ2 = CH2_CZ2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD1_CD1 = HD1_CD1.bondType.distance;
                double dHE1_NE1 = HE1_NE1.bondType.distance;
                double dHE3_CE3 = HE3_CE3.bondType.distance;
                double dHZ2_CZ2 = HZ2_CZ2.bondType.distance;
                double dHZ3_CZ3 = HZ3_CZ3.bondType.distance;
                double dHH2_CH2 = HH2_CH2.bondType.distance;

                double dCG_CB_CA = getAngle(name, CG, CB, CA);
                double dCD1_CG_CB = getAngle(name, CD1, CG, CB);
                double dCD2_CG_CB = getAngle(name, CD2, CG, CB);
                double dNE1_CD1_CG = getAngle(name, NE1, CD1, CG);
                double dCE2_NE1_CD1 = getAngle(name, CE2, NE1, CD1);
                double dCE3_CD2_CE2 = getAngle(name, CE3, CD2, CE2);
                double dCZ2_CE2_CD2 = getAngle(name, CZ2, CE2, CD2);
                double dCZ3_CE3_CD2 = getAngle(name, CZ3, CE3, CD2);
                double dCH2_CZ2_CE2 = getAngle(name, CH2, CZ2, CE2);
                double dHB_CB_CA = getAngle(name, HB2, CB, CA);
                double dHD1_CD1_CG = getAngle(name, HD1, CD1, CG);
                double dHD1_CD1_NE1 = getAngle(name, HD1, CD1, NE1);
                double dHE1_NE1_CD1 = getAngle(name, HE1, NE1, CD1);
                double dHE1_NE1_CE2 = getAngle(name, HE1, NE1, CE2);
                double dHE3_CE3_CD2 = getAngle(name, HE3, CE3, CD2);
                double dHE3_CE3_CZ3 = getAngle(name, HE3, CE3, CZ3);
                double dHZ2_CZ2_CE2 = getAngle(name, HZ2, CZ2, CE2);
                double dHZ2_CZ2_CH2 = getAngle(name, HZ2, CZ2, CH2);
                double dHZ3_CZ3_CE3 = getAngle(name, HZ3, CZ3, CE3);
                double dHZ3_CZ3_CH2 = getAngle(name, HZ3, CZ3, CH2);
                double dHH2_CH2_CZ2 = getAngle(name, HH2, CH2, CZ2);
                double dHH2_CH2_CZ3 = getAngle(name, HH2, CH2, CZ3);

                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD1, CG, dCD1_CG, CB, dCD1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, CA, rotamer.chi2 + 180, 0);
                intxyz(NE1, CD1, dNE1_CD1, CG, dNE1_CD1_CG, CD2, 0.0, 0);
                intxyz(CE2, NE1, dCE2_NE1, CD1, dCE2_NE1_CD1, CG, 0.0, 0);
                intxyz(CE3, CD2, dCE3_CD2, CE2, dCE3_CD2_CE2, NE1, 180.0, 0);
                intxyz(CZ2, CE2, dCZ2_CE2, CD2, dCZ2_CE2_CD2, CE3, 0.0, 0);
                intxyz(CZ3, CE3, dCZ3_CE3, CD2, dCZ3_CE3_CD2, CE2, 0.0, 0);
                intxyz(CH2, CZ2, dCH2_CZ2, CE2, dCH2_CZ2_CE2, CD2, 0.0, 0);
                // Continue using 109.4 degrees for tetrahedral hydrogens.
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, CD1, dHD1_CD1, CG, dHD1_CD1_CG, NE1, dHD1_CD1_NE1, 3);
                intxyz(HE1, NE1, dHE1_NE1, CD1, dHE1_NE1_CD1, CE2, dHE1_NE1_CE2, 3);
                intxyz(HE3, CE3, dHE3_CE3, CD2, dHE3_CE3_CD2, CZ3, dHE3_CE3_CZ3, 3);
                intxyz(HZ2, CZ2, dHZ2_CZ2, CE2, dHZ2_CZ2_CE2, CH2, dHZ2_CZ2_CH2, 3);
                intxyz(HZ3, CZ3, dHZ3_CZ3, CE3, dHZ3_CZ3_CE3, CH2, dHZ3_CZ3_CH2, 3);
                intxyz(HH2, CH2, dHH2_CH2, CZ2, dHH2_CH2_CZ2, CZ3, dHH2_CH2_CZ3, 3);
                break;
            }
            case HIS: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Bond CG_CB = CG.getBond(CB);
                Bond ND1_CG = ND1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond CE1_ND1 = CE1.getBond(ND1);
                Bond NE2_CD2 = NE2.getBond(CD2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD1_ND1 = HD1.getBond(ND1);
                Bond HD2_CD2 = HD2.getBond(CD2);
                Bond HE1_CE1 = HE1.getBond(CE1);
                Bond HE2_NE2 = HE2.getBond(NE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dND1_CG = ND1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dCE1_ND1 = CE1_ND1.bondType.distance;
                double dNE2_CD2 = NE2_CD2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD1_ND1 = HD1_ND1.bondType.distance;
                double dHD2_CD2 = HD2_CD2.bondType.distance;
                double dHE1_CE1 = HE1_CE1.bondType.distance;
                double dHE2_NE2 = HE2_NE2.bondType.distance;

                double dCG_CB_CA = getAngle(name, CG, CB, CA);
                double dND1_CG_CB = getAngle(name, ND1, CG, CB);
                double dCD2_CG_CB = getAngle(name, CD2, CG, CB);
                double dCE1_ND1_CG = getAngle(name, CE1, ND1, CG);
                double dNE2_CD2_CG = getAngle(name, NE2, CD2, CG);
                double dHB_CB_CA = getAngle(name, HB2, CB, CA);
                double dHD1_ND1_CG = getAngle(name, HD1, ND1, CG);
                double dHD1_ND1_CE1 = getAngle(name, HD1, ND1, CE1);
                double dHD2_CD2_CG = getAngle(name, HD2, CD2, CG);
                double dHD2_CD2_NE2 = getAngle(name, HD2, CD2, NE2);
                double dHE1_CE1_ND1 = getAngle(name, HE1, CE1, ND1);
                double dHE1_CE1_NE2 = getAngle(name, HE1, CE1, NE2);
                double dHE2_NE2_CD2 = getAngle(name, HE2, NE2, CD2);
                double dHE2_NE2_CE1 = getAngle(name, HE2, NE2, CE1);

                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(ND1, CG, dND1_CG, CB, dND1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, CA, rotamer.chi2 + 180, 0);
                intxyz(CE1, ND1, dCE1_ND1, CG, dCE1_ND1_CG, CD2, 0.0, 0);
                intxyz(NE2, CD2, dNE2_CD2, CG, dNE2_CD2_CG, ND1, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, ND1, dHD1_ND1, CG, dHD1_ND1_CG, CE1, dHD1_ND1_CE1, 3);
                intxyz(HD2, CD2, dHD2_CD2, CG, dHD2_CD2_CG, NE2, dHD2_CD2_NE2, 3);
                intxyz(HE1, CE1, dHE1_CE1, ND1, dHE1_CE1_ND1, NE2, dHE1_CE1_NE2, 3);
                intxyz(HE2, NE2, dHE2_NE2, CD2, dHE2_NE2_CD2, CE1, dHE2_NE2_CE1, 3);
                break;
            }
            case HID: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD1 = (Atom) residue.getAtomNode("HD1");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Bond CG_CB = CG.getBond(CB);
                Bond ND1_CG = ND1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond CE1_ND1 = CE1.getBond(ND1);
                Bond NE2_CD2 = NE2.getBond(CD2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD1_ND1 = HD1.getBond(ND1);
                Bond HD2_CD2 = HD2.getBond(CD2);
                Bond HE1_CE1 = HE1.getBond(CE1);
                double dCG_CB = CG_CB.bondType.distance;
                double dND1_CG = ND1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dCE1_ND1 = CE1_ND1.bondType.distance;
                double dNE2_CD2 = NE2_CD2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD1_ND1 = HD1_ND1.bondType.distance;
                double dHD2_CD2 = HD2_CD2.bondType.distance;
                double dHE1_CE1 = HE1_CE1.bondType.distance;

                double dCG_CB_CA = getAngle(name, CG, CB, CA);
                double dND1_CG_CB = getAngle(name, ND1, CG, CB);
                double dCD2_CG_CB = getAngle(name, CD2, CG, CB);
                double dCE1_ND1_CG = getAngle(name, CE1, ND1, CG);
                double dNE2_CD2_CG = getAngle(name, NE2, CD2, CG);
                double dHB_CB_CA = getAngle(name, HB2, CB, CA);
                double dHD1_ND1_CG = getAngle(name, HD1, ND1, CG);
                double dHD1_ND1_CE1 = getAngle(name, HD1, ND1, CE1);
                double dHD2_CD2_CG = getAngle(name, HD2, CD2, CG);
                double dHD2_CD2_NE2 = getAngle(name, HD2, CD2, NE2);
                double dHE1_CE1_ND1 = getAngle(name, HE1, CE1, ND1);
                double dHE1_CE1_NE2 = getAngle(name, HE1, CE1, NE2);

                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(ND1, CG, dND1_CG, CB, dND1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, CA, rotamer.chi2 + 180, 0);
                intxyz(CE1, ND1, dCE1_ND1, CG, dCE1_ND1_CG, CD2, 0.0, 0);
                intxyz(NE2, CD2, dNE2_CD2, CG, dNE2_CD2_CG, ND1, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD1, ND1, dHD1_ND1, CG, dHD1_ND1_CG, CE1, dHD1_ND1_CE1, 3);
                intxyz(HD2, CD2, dHD2_CD2, CG, dHD2_CD2_CG, NE2, dHD2_CD2_NE2, 3);
                intxyz(HE1, CE1, dHE1_CE1, ND1, dHE1_CE1_ND1, NE2, dHE1_CE1_NE2, 3);
                break;
            }
            case HIE: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom ND1 = (Atom) residue.getAtomNode("ND1");
                Atom CD2 = (Atom) residue.getAtomNode("CD2");
                Atom CE1 = (Atom) residue.getAtomNode("CE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Bond CG_CB = CG.getBond(CB);
                Bond ND1_CG = ND1.getBond(CG);
                Bond CD2_CG = CD2.getBond(CG);
                Bond CE1_ND1 = CE1.getBond(ND1);
                Bond NE2_CD2 = NE2.getBond(CD2);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD2_CD2 = HD2.getBond(CD2);
                Bond HE1_CE1 = HE1.getBond(CE1);
                Bond HE2_NE2 = HE2.getBond(NE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dND1_CG = ND1_CG.bondType.distance;
                double dCD2_CG = CD2_CG.bondType.distance;
                double dCE1_ND1 = CE1_ND1.bondType.distance;
                double dNE2_CD2 = NE2_CD2.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD2_CD2 = HD2_CD2.bondType.distance;
                double dHE1_CE1 = HE1_CE1.bondType.distance;
                double dHE2_NE2 = HE2_NE2.bondType.distance;

                double dCG_CB_CA = getAngle(name, CG, CB, CA);
                double dND1_CG_CB = getAngle(name, ND1, CG, CB);
                double dCD2_CG_CB = getAngle(name, CD2, CG, CB);
                double dCE1_ND1_CG = getAngle(name, CE1, ND1, CG);
                double dNE2_CD2_CG = getAngle(name, NE2, CD2, CG);
                double dHB_CB_CA = getAngle(name, HB2, CB, CA);
                double dHD2_CD2_CG = getAngle(name, HD2, CD2, CG);
                double dHD2_CD2_NE2 = getAngle(name, HD2, CD2, NE2);
                double dHE1_CE1_ND1 = getAngle(name, HE1, CE1, ND1);
                double dHE1_CE1_NE2 = getAngle(name, HE1, CE1, NE2);
                double dHE2_NE2_CD2 = getAngle(name, HE2, NE2, CD2);
                double dHE2_NE2_CE1 = getAngle(name, HE2, NE2, CE1);

                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(ND1, CG, dND1_CG, CB, dND1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CD2, CG, dCD2_CG, CB, dCD2_CG_CB, CA, rotamer.chi2 + 180, 0);
                intxyz(CE1, ND1, dCE1_ND1, CG, dCE1_ND1_CG, CD2, 0.0, 0);
                intxyz(NE2, CD2, dNE2_CD2, CG, dNE2_CD2_CG, ND1, 0.0, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HD2, CD2, dHD2_CD2, CG, dHD2_CD2_CG, NE2, dHD2_CD2_NE2, 3);
                intxyz(HE1, CE1, dHE1_CE1, ND1, dHE1_CE1_ND1, NE2, dHE1_CE1_NE2, 3);
                intxyz(HE2, NE2, dHE2_NE2, CD2, dHE2_NE2_CD2, CE1, dHE2_NE2_CE1, 3);
                break;
            }
            case ASP: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                Atom OD2 = (Atom) residue.getAtomNode("OD2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Bond CG_CB = CG.getBond(CB);
                Bond OD1_CG = OD1.getBond(CG);
                Bond OD2_CG = OD2.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                double dCG_CB = CG_CB.bondType.distance;
                double dOD1_CG = OD1_CG.bondType.distance;
                double dOD2_CG = OD2_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle OD1_CG_CB = OD1.getAngle(CG, CB);
                Angle OD2_CG_CB = OD2.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
                double dOD2_CG_CB = OD2_CG_CB.angleType.angle[OD2_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(OD1, CG, dOD1_CG, CB, dOD1_CG_CB, CA, 0.0, 0);
                intxyz(OD2, CG, dOD2_CG, CB, dOD2_CG_CB, OD1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, -1);
                break;
            }
            case ASH: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                Atom OD2 = (Atom) residue.getAtomNode("OD2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Bond CG_CB = CG.getBond(CB);
                Bond OD1_CG = OD1.getBond(CG);
                Bond OD2_CG = OD2.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD2_OD2 = HD2.getBond(OD2);
                double dCG_CB = CG_CB.bondType.distance;
                double dOD1_CG = OD1_CG.bondType.distance;
                double dOD2_CG = OD2_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD2_OD2 = HD2_OD2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle OD1_CG_CB = OD1.getAngle(CG, CB);
                Angle OD2_CG_CB = OD2.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD2_OD2_CG = HD2.getAngle(OD2, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
                double dOD2_CG_CB = OD2_CG_CB.angleType.angle[OD2_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD2_OD2_CG = HD2_OD2_CG.angleType.angle[HD2_OD2_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(OD1, CG, dOD1_CG, CB, dOD1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OD2, CG, dOD2_CG, CB, dOD2_CG_CB, OD1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, -1);
                intxyz(HD2, OD2, dHD2_OD2, CG, dHD2_OD2_CG, OD1, 0.0, 0);
                break;
            }
            case ASN: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom OD1 = (Atom) residue.getAtomNode("OD1");
                Atom ND2 = (Atom) residue.getAtomNode("ND2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HD21 = (Atom) residue.getAtomNode("HD21");
                Atom HD22 = (Atom) residue.getAtomNode("HD22");
                Bond CG_CB = CG.getBond(CB);
                Bond OD1_CG = OD1.getBond(CG);
                Bond ND2_CG = ND2.getBond(CG);
                Bond HB_CB = HB2.getBond(CB);
                Bond HD2_ND2 = HD21.getBond(ND2);
                double dCG_CB = CG_CB.bondType.distance;
                double dOD1_CG = OD1_CG.bondType.distance;
                double dND2_CG = ND2_CG.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHD2_ND2 = HD2_ND2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle OD1_CG_CB = OD1.getAngle(CG, CB);
                Angle ND2_CG_CB = ND2.getAngle(CG, CB);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HD2_ND2_CG = HD21.getAngle(ND2, CG);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
                double dND2_CG_CB = ND2_CG_CB.angleType.angle[ND2_CG_CB.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHD2_ND2_CG = HD2_ND2_CG.angleType.angle[HD2_ND2_CG.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(OD1, CG, dOD1_CG, CB, dOD1_CG_CB, CA, rotamer.chi2, 0);
                intxyz(ND2, CG, dND2_CG, CB, dND2_CG_CB, OD1, 124.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 107.9, -1);
                intxyz(HD21, ND2, dHD2_ND2, CG, dHD2_ND2_CG, CB, 0.0, 0);
                intxyz(HD22, ND2, dHD2_ND2, CG, dHD2_ND2_CG, HD21, 120.0, 1);
                break;
            }
            case GLU: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                Atom OE2 = (Atom) residue.getAtomNode("OE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond OE1_CD = OE1.getBond(CD);
                Bond OE2_CD = OE2.getBond(CD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dOE1_CD = OE1_CD.bondType.distance;
                double dOE2_CD = OE2_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle OE1_CD_CG = OE1.getAngle(CD, CG);
                Angle OE2_CD_CG = OE2.getAngle(CD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
                double dOE2_CD_CG = OE2_CD_CG.angleType.angle[OE2_CD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, dOE1_CD, CG, dOE1_CD_CG, CB, rotamer.chi3, 0);
                intxyz(OE2, CD, dOE2_CD, CG, dOE2_CD_CG, OE1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, -1);
                break;
            }
            case GLH: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                Atom OE2 = (Atom) residue.getAtomNode("OE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond OE1_CD = OE1.getBond(CD);
                Bond OE2_CD = OE2.getBond(CD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HE2_OE2 = HE2.getBond(OE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dOE1_CD = OE1_CD.bondType.distance;
                double dOE2_CD = OE2_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHE2_OE2 = HE2_OE2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle OE1_CD_CG = OE1.getAngle(CD, CG);
                Angle OE2_CD_CG = OE2.getAngle(CD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HE2_OE2_CD = HE2.getAngle(OE2, CD);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
                double dOE2_CD_CG = OE2_CD_CG.angleType.angle[OE2_CD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHE2_OE2_CD = HE2_OE2_CD.angleType.angle[HE2_OE2_CD.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, dOE1_CD, CG, dOE1_CD_CG, CB, rotamer.chi3, 0);
                intxyz(OE2, CD, dOE2_CD, CG, dOE2_CD_CG, OE1, 126.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, -1);
                intxyz(HE2, OE2, dHE2_OE2, CD, dHE2_OE2_CD, OE1, 0.0, 0);
                break;
            }
            case GLN: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom OE1 = (Atom) residue.getAtomNode("OE1");
                Atom NE2 = (Atom) residue.getAtomNode("NE2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HE21 = (Atom) residue.getAtomNode("HE21");
                Atom HE22 = (Atom) residue.getAtomNode("HE22");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond OE1_CD = OE1.getBond(CD);
                Bond NE2_CD = NE2.getBond(CD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HE2_NE2 = HE21.getBond(NE2);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dOE1_CD = OE1_CD.bondType.distance;
                double dNE2_CD = NE2_CD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHE2_NE2 = HE2_NE2.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle OE1_CD_CG = OE1.getAngle(CD, CG);
                Angle NE2_CD_CG = NE2.getAngle(CD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HE2_NE2_CD = HE21.getAngle(NE2, CD);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
                double dNE2_CD_CG = NE2_CD_CG.angleType.angle[NE2_CD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHE2_NE2_CD = HE2_NE2_CD.angleType.angle[HE2_NE2_CD.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(OE1, CD, dOE1_CD, CG, dOE1_CD_CG, CB, rotamer.chi3, 0);
                intxyz(NE2, CD, dNE2_CD, CG, dNE2_CD_CG, OE1, 124.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 107.9, -1);
                intxyz(HE21, NE2, dHE2_NE2, CD, dHE2_NE2_CD, CG, 0.0, 0);
                intxyz(HE22, NE2, dHE2_NE2, CD, dHE2_NE2_CD, HE21, 120.0, 1);
                break;
            }
            case MET: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom SD = (Atom) residue.getAtomNode("SD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HE1 = (Atom) residue.getAtomNode("HE1");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Bond CG_CB = CG.getBond(CB);
                Bond SD_CG = SD.getBond(CG);
                Bond CE_SD = CE.getBond(SD);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HE_CE = HE1.getBond(CE);
                double dCG_CB = CG_CB.bondType.distance;
                double dSD_CG = SD_CG.bondType.distance;
                double dCE_SD = CE_SD.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle SD_CG_CB = SD.getAngle(CG, CB);
                Angle CE_SD_CG = CE.getAngle(SD, CG);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HE_CE_SD = HE1.getAngle(CE, SD);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dSD_CG_CB = SD_CG_CB.angleType.angle[SD_CG_CB.nh];
                double dCE_SD_CG = CE_SD_CG.angleType.angle[CE_SD_CG.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHE_CE_SD = HE_CE_SD.angleType.angle[HE_CE_SD.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(SD, CG, dSD_CG, CB, dSD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CE, SD, dCE_SD, CG, dCE_SD_CG, CB, rotamer.chi3, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, SD, 112.0, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, SD, 112.0, -1);
                intxyz(HE1, CE, dHE_CE, SD, dHE_CE_SD, CG, 180.0, 0);
                intxyz(HE2, CE, dHE_CE, SD, dHE_CE_SD, HE1, 109.4, 1);
                intxyz(HE3, CE, dHE_CE, SD, dHE_CE_SD, HE1, 109.4, -1);
                break;
            }
            case LYS: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom NZ = (Atom) residue.getAtomNode("NZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Atom HZ1 = (Atom) residue.getAtomNode("HZ1");
                Atom HZ2 = (Atom) residue.getAtomNode("HZ2");
                Atom HZ3 = (Atom) residue.getAtomNode("HZ3");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond CE_CD = CE.getBond(CD);
                Bond NZ_CE = NZ.getBond(CE);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                Bond HE_CE = HE2.getBond(CE);
                Bond HZ_NZ = HZ1.getBond(NZ);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dNZ_CE = NZ_CE.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                double dHZ_NZ = HZ_NZ.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle CE_CD_CG = CE.getAngle(CD, CG);
                Angle NZ_CE_CD = NZ.getAngle(CE, CD);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                Angle HE_CE_CD = HE2.getAngle(CE, CD);
                Angle HZ_NZ_CE = HZ1.getAngle(NZ, CE);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
                double dNZ_CE_CD = NZ_CE_CD.angleType.angle[NZ_CE_CD.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
                double dHZ_NZ_CE = HZ_NZ_CE.angleType.angle[HZ_NZ_CE.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CE, CD, dCE_CD, CG, dCE_CD_CG, CB, rotamer.chi3, 0);
                intxyz(NZ, CE, dNZ_CE, CD, dNZ_CE_CD, CG, rotamer.chi4, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, -1);
                intxyz(HE2, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, 1);
                intxyz(HE3, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, -1);
                intxyz(HZ1, NZ, dHZ_NZ, CE, dHZ_NZ_CE, CD, 180.0, 0);
                intxyz(HZ2, NZ, dHZ_NZ, CE, dHZ_NZ_CE, HZ1, 109.5, 1);
                intxyz(HZ3, NZ, dHZ_NZ, CE, dHZ_NZ_CE, HZ1, 109.5, -1);
                break;
            }
            case LYD: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom CE = (Atom) residue.getAtomNode("CE");
                Atom NZ = (Atom) residue.getAtomNode("NZ");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Atom HE2 = (Atom) residue.getAtomNode("HE2");
                Atom HE3 = (Atom) residue.getAtomNode("HE3");
                Atom HZ1 = (Atom) residue.getAtomNode("HZ1");
                Atom HZ2 = (Atom) residue.getAtomNode("HZ2");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond CE_CD = CE.getBond(CD);
                Bond NZ_CE = NZ.getBond(CE);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                Bond HE_CE = HE2.getBond(CE);
                Bond HZ_NZ = HZ1.getBond(NZ);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dCE_CD = CE_CD.bondType.distance;
                double dNZ_CE = NZ_CE.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_CE = HE_CE.bondType.distance;
                double dHZ_NZ = HZ_NZ.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle CE_CD_CG = CE.getAngle(CD, CG);
                Angle NZ_CE_CD = NZ.getAngle(CE, CD);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                Angle HE_CE_CD = HE2.getAngle(CE, CD);
                Angle HZ_NZ_CE = HZ1.getAngle(NZ, CE);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dCE_CD_CG = CE_CD_CG.angleType.angle[CE_CD_CG.nh];
                double dNZ_CE_CD = NZ_CE_CD.angleType.angle[NZ_CE_CD.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_CE_CD = HE_CE_CD.angleType.angle[HE_CE_CD.nh];
                double dHZ_NZ_CE = HZ_NZ_CE.angleType.angle[HZ_NZ_CE.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(CE, CD, dCE_CD, CG, dCE_CD_CG, CB, rotamer.chi3, 0);
                intxyz(NZ, CE, dNZ_CE, CD, dNZ_CE_CD, CG, rotamer.chi4, 0);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, CE, 109.4, -1);
                intxyz(HE2, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, 1);
                intxyz(HE3, CE, dHE_CE, CD, dHE_CE_CD, NZ, 108.8, -1);
                intxyz(HZ1, NZ, dHZ_NZ, CE, dHZ_NZ_CE, CD, 180.0, 0);
                intxyz(HZ2, NZ, dHZ_NZ, CE, dHZ_NZ_CE, HZ1, 109.5, 1);
                break;
            }
            case ARG: {
                Atom CA = (Atom) residue.getAtomNode("CA");
                Atom CB = (Atom) residue.getAtomNode("CB");
                Atom N = (Atom) residue.getAtomNode("N");
                Atom CG = (Atom) residue.getAtomNode("CG");
                Atom CD = (Atom) residue.getAtomNode("CD");
                Atom NE = (Atom) residue.getAtomNode("NE");
                Atom CZ = (Atom) residue.getAtomNode("CZ");
                Atom NH1 = (Atom) residue.getAtomNode("NH1");
                Atom NH2 = (Atom) residue.getAtomNode("NH2");
                Atom HB2 = (Atom) residue.getAtomNode("HB2");
                Atom HB3 = (Atom) residue.getAtomNode("HB3");
                Atom HG2 = (Atom) residue.getAtomNode("HG2");
                Atom HG3 = (Atom) residue.getAtomNode("HG3");
                Atom HD2 = (Atom) residue.getAtomNode("HD2");
                Atom HD3 = (Atom) residue.getAtomNode("HD3");
                Atom HE = (Atom) residue.getAtomNode("HE");
                Atom HH11 = (Atom) residue.getAtomNode("HH11");
                Atom HH12 = (Atom) residue.getAtomNode("HH12");
                Atom HH21 = (Atom) residue.getAtomNode("HH21");
                Atom HH22 = (Atom) residue.getAtomNode("HH22");
                Bond CG_CB = CG.getBond(CB);
                Bond CD_CG = CD.getBond(CG);
                Bond NE_CD = NE.getBond(CD);
                Bond CZ_NE = CZ.getBond(NE);
                Bond NH_CZ = NH1.getBond(CZ);
                Bond HB_CB = HB2.getBond(CB);
                Bond HG_CG = HG2.getBond(CG);
                Bond HD_CD = HD2.getBond(CD);
                Bond HE_NE = HE.getBond(NE);
                Bond HH_NH = HH11.getBond(NH1);
                double dCG_CB = CG_CB.bondType.distance;
                double dCD_CG = CD_CG.bondType.distance;
                double dNE_CD = NE_CD.bondType.distance;
                double dCZ_NE = CZ_NE.bondType.distance;
                double dNH_CZ = NH_CZ.bondType.distance;
                double dHB_CB = HB_CB.bondType.distance;
                double dHG_CG = HG_CG.bondType.distance;
                double dHD_CD = HD_CD.bondType.distance;
                double dHE_NE = HE_NE.bondType.distance;
                double dHH_NH = HH_NH.bondType.distance;
                Angle CG_CB_CA = CG.getAngle(CB, CA);
                Angle CD_CG_CB = CD.getAngle(CG, CB);
                Angle NE_CD_CG = NE.getAngle(CD, CG);
                Angle CZ_NE_CD = CZ.getAngle(NE, CD);
                Angle NH_CZ_NE = NH1.getAngle(CZ, NE);
                Angle HB_CB_CA = HB2.getAngle(CB, CA);
                Angle HG_CG_CB = HG2.getAngle(CG, CB);
                Angle HD_CD_CG = HD2.getAngle(CD, CG);
                Angle HE_NE_CD = HE.getAngle(NE, CD);
                Angle HH_NH_CZ = HH11.getAngle(NH1, CZ);
                double dCG_CB_CA = CG_CB_CA.angleType.angle[CG_CB_CA.nh];
                double dCD_CG_CB = CD_CG_CB.angleType.angle[CD_CG_CB.nh];
                double dNE_CD_CG = NE_CD_CG.angleType.angle[NE_CD_CG.nh];
                double dCZ_NE_CD = CZ_NE_CD.angleType.angle[CZ_NE_CD.nh];
                double dNH_CZ_NE = NH_CZ_NE.angleType.angle[NH_CZ_NE.nh];
                double dHB_CB_CA = HB_CB_CA.angleType.angle[HB_CB_CA.nh];
                double dHG_CG_CB = HG_CG_CB.angleType.angle[HG_CG_CB.nh];
                double dHD_CD_CG = HD_CD_CG.angleType.angle[HD_CD_CG.nh];
                double dHE_NE_CD = HE_NE_CD.angleType.angle[HE_NE_CD.nh];
                double dHH_NH_CZ = HH_NH_CZ.angleType.angle[HH_NH_CZ.nh];
                intxyz(CG, CB, dCG_CB, CA, dCG_CB_CA, N, rotamer.chi1, 0);
                intxyz(CD, CG, dCD_CG, CB, dCD_CG_CB, CA, rotamer.chi2, 0);
                intxyz(NE, CD, dNE_CD, CG, dNE_CD_CG, CB, rotamer.chi3, 0);
                intxyz(CZ, NE, dCZ_NE, CD, dCZ_NE_CD, CG, rotamer.chi4, 0);
                intxyz(NH1, CZ, dNH_CZ, NE, dNH_CZ_NE, CD, 180, 0);
                intxyz(NH2, CZ, dNH_CZ, NE, dNH_CZ_NE, NH1, 120.0, 1);
                intxyz(HB2, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, 1);
                intxyz(HB3, CB, dHB_CB, CA, dHB_CB_CA, CG, 109.4, -1);
                intxyz(HG2, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, 1);
                intxyz(HG3, CG, dHG_CG, CB, dHG_CG_CB, CD, 109.4, -1);
                intxyz(HD2, CD, dHD_CD, CG, dHD_CD_CG, NE, 109.4, 1);
                intxyz(HD3, CD, dHD_CD, CG, dHD_CD_CG, NE, 109.4, -1);
                intxyz(HE, NE, dHE_NE, CD, dHE_NE_CD, CZ, 120.0, 1);
                intxyz(HH11, NH1, dHH_NH, CZ, dHH_NH_CZ, NE, 180.0, 0);
                intxyz(HH12, NH1, dHH_NH, CZ, dHH_NH_CZ, HH11, 120.0, 1);
                intxyz(HH21, NH2, dHH_NH, CZ, dHH_NH_CZ, NE, 180.0, 0);
                intxyz(HH22, NH2, dHH_NH, CZ, dHH_NH_CZ, HH21, 120.0, 1);
                break;
            }
            case UNK:
                String resName = residue.getName().toUpperCase();
                if (nonstdRotCache.containsKey(resName)) {
                    nonstdRotCache.get(resName).applyNonstdRotamer(residue, rotamer);
                }
                break;
            default:
                break;
        }
    }

    /**
     * Applies a nucleic acid Rotamer, returning the magnitude of the correction applied to make residue i join i-1.
     *
     * Note that the independent flag is separate from DEE independence: DEE
     * independence is preserved by applying corrections based on a non-variable
     * set of coordinates, and is wholly independent of what is happening to
     * residue i-1.
     *
     * Cannot presently handle 3' phosphate caps: I do not know what they would
     * be labeled as in PDB files. A template for how to handle 3' phosphate
     * caps is written but commented out.
     *
     * @param residue Residue.
     * @param rotamer Rotamer to be applied to Residue.
     * @param independent Whether to draw NA rotamer independent of chain
     * context.
     * @return Magnitude of the correction vector.
     */
    public static double applyNARotamer(Residue residue, Rotamer rotamer, boolean independent) {
        if (rotamer.isState) {
            applyState(residue, rotamer);
            return 0;
        }
        NucleicAcid3 na = NucleicAcid3.valueOf(residue.getName());
        Residue prevResidue = residue.getPreviousResidue();
        boolean is3sTerminal = false;  // 3' terminal
        if (residue.getNextResidue() == null) {
            is3sTerminal = true;
        }

        // Check if this is a 3' phosphate being listed as its own residue.
        /* if (residue.getAtomList().size() == 1) {
         return;
         } */
        boolean isDeoxy = residue.getAtomNode("O2\'") == null;

        // Note: chi values will generally be applied from chi7 to chi1.
        // Will have to add an else-if to handle DNA C3'-exo configurations.
        NucleicSugarPucker sugarPucker = NucleicSugarPucker.checkPucker(rotamer.chi7, isDeoxy);
        NucleicSugarPucker prevSugarPucker = NucleicSugarPucker.checkPucker(rotamer.chi1, isDeoxy);

        // Revert C1', O4', and C4' coordinates to PDB defaults.
        Atom C1s = (Atom) residue.getAtomNode("C1\'");
        C1s.moveTo(residue.getC1sCoords());
        Atom O4s = (Atom) residue.getAtomNode("O4\'");
        O4s.moveTo(residue.getO4sCoords());
        Atom C4s = (Atom) residue.getAtomNode("C4\'");
        C4s.moveTo(residue.getC4sCoords());

        // Presently, the exterior method loadPriorAtomicCoordinates() directly
        // calls applySugarPucker instead of going through applyRotamer().
        applySugarPucker(residue, sugarPucker, isDeoxy, true);
        applyNABackbone(residue, rotamer, prevResidue);

        double naCorrection = 0;
        if (prevResidue != null && !independent) {
            naCorrection = applyNACorrections(residue, prevResidue, rotamer, prevSugarPucker, isDeoxy, is3sTerminal);
        } /* else if (!independent) {
         startingResidueConsistencyCheck(residue, rotamer, correctionThreshold);
         } */

        applyNASideAtoms(residue, rotamer, prevResidue, isDeoxy, is3sTerminal, prevSugarPucker);
        return naCorrection;
    }

    /**
     * If place is true, builds C2', C3', and O3' based on delta(i) and returns
     * an empty double[]; if place is false, returns a double[] filled with the
     * coordinates at which O3' would be placed by the specified pucker.
     *
     * Presently uses default locations for C1', O4', and C4' to build these
     * atoms.
     *
     * @param residue Nucleic acid Residue to which the pucker is to bea applied
     * @param pucker An int specifying pucker (1=North, 2=South).
     * @param isDeoxy Boolean
     * @param place Flag for usage case.
     * @return A double[] with O3' coordinates (place=false), or null
     * (place=true).
     */
    public static double[] applySugarPucker(Residue residue, NucleicSugarPucker pucker, boolean isDeoxy, boolean place) {
        // Torsions from http://ndb-mirror-2.rutgers.edu/ndbmodule/archives/proj/valence/table6.html
        // SP is short for Sugar Pucker (torsion).
        final double C2_SP_SOUTH_RNA = 24.2;
        final double C3_SP_SOUTH_RNA = 357.7;
        final double O3_SP_SOUTH_RNA = 268.1;
        final double C2_SP_NORTH_RNA = 324.7;
        final double C3_SP_NORTH_RNA = 20.4;
        final double O3_SP_NORTH_RNA = 201.8;

        final double C2_SP_SOUTH_DNA = 22.6;
        final double C3_SP_SOUTH_DNA = 357.7;
        final double O3_SP_SOUTH_DNA = 265.8;
        final double C2_SP_NORTH_DNA = 327.7;
        final double C3_SP_NORTH_DNA = 16.4;
        final double O3_SP_NORTH_DNA = 205.4;

        // Ret will only be filled if place is false.
        double[] ret = new double[3];

        // Constituents of the sugar
        Atom C1s = (Atom) residue.getAtomNode("C1\'");
        // C2' will only be used if place is true.
        Atom C3s = (Atom) residue.getAtomNode("C3\'");
        Atom O4s = (Atom) residue.getAtomNode("O4\'");
        Atom C4s = (Atom) residue.getAtomNode("C4\'");
        Atom O3s = (Atom) residue.getAtomNode("O3\'");

        // Bonds and angles necessary to draw C3' and O3'.
        Bond C3s_C4s = C4s.getBond(C3s);
        double dC3s_C4s = C3s_C4s.bondType.distance;
        Angle C3s_C4s_O4s = O4s.getAngle(C4s, C3s);
        double dC3s_C4s_O4s = C3s_C4s_O4s.angleType.angle[C3s_C4s_O4s.nh];

        Bond C3s_O3s = C3s.getBond(O3s);
        double dC3s_O3s = C3s_O3s.bondType.distance;
        Angle C4s_C3s_O3s = C4s.getAngle(C3s, O3s);
        double dC4s_C3s_O3s = C4s_C3s_O3s.angleType.angle[C4s_C3s_O3s.nh];

        /*
         * If place is true, place C3', C2', and O3'. Else, determine the
         * coordinates of O3' based on the invariant atoms.
         *
         * If deoxy, a different set of sugar pucker torsions is
         * applied.
         *
         * Then, if a North pucker (delta (i) is in the range of 78-90),
         * apply Cx_SP_NORTH torsion to place C3 and C2.  Else,
         * assume a South pucker (delta is between 140 and 152).
         *
         * Then, use O3_SP_NORTH or SOUTH to place O3.
         */
        if (place) {
            Atom C2s = (Atom) residue.getAtomNode("C2\'");
            Bond C3s_C2s = C3s.getBond(C2s);
            double dC3s_C2s = C3s_C2s.bondType.distance;
            Angle C4s_C3s_C2s = C4s.getAngle(C3s, C2s);
            double dC4s_C3s_C2s = C4s_C3s_C2s.angleType.angle[C4s_C3s_C2s.nh];

            if (isDeoxy) {
                if (pucker == NucleicSugarPucker.C3_ENDO) {
                    intxyz(C3s, C4s, dC3s_C4s, O4s, dC3s_C4s_O4s, C1s, C3_SP_NORTH_DNA, 0);
                    intxyz(C2s, C3s, dC3s_C2s, C4s, dC4s_C3s_C2s, O4s, C2_SP_NORTH_DNA, 0);
                    intxyz(O3s, C3s, dC3s_O3s, C4s, dC4s_C3s_O3s, O4s, O3_SP_NORTH_DNA, 0);
                } // TODO: else-if for 3'-exo configuration (DNA only)
                else {
                    intxyz(C3s, C4s, dC3s_C4s, O4s, dC3s_C4s_O4s, C1s, C3_SP_SOUTH_DNA, 0);
                    intxyz(C2s, C3s, dC3s_C2s, C4s, dC4s_C3s_C2s, O4s, C2_SP_SOUTH_DNA, 0);
                    intxyz(O3s, C3s, dC3s_O3s, C4s, dC4s_C3s_O3s, O4s, O3_SP_SOUTH_DNA, 0);
                }
            } else {
                if (pucker == NucleicSugarPucker.C3_ENDO) {
                    intxyz(C3s, C4s, dC3s_C4s, O4s, dC3s_C4s_O4s, C1s, C3_SP_NORTH_RNA, 0);
                    intxyz(C2s, C3s, dC3s_C2s, C4s, dC4s_C3s_C2s, O4s, C2_SP_NORTH_RNA, 0);
                    intxyz(O3s, C3s, dC3s_O3s, C4s, dC4s_C3s_O3s, O4s, O3_SP_NORTH_RNA, 0);
                } else {
                    intxyz(C3s, C4s, dC3s_C4s, O4s, dC3s_C4s_O4s, C1s, C3_SP_SOUTH_RNA, 0);
                    intxyz(C2s, C3s, dC3s_C2s, C4s, dC4s_C3s_C2s, O4s, C2_SP_SOUTH_RNA, 0);
                    intxyz(O3s, C3s, dC3s_O3s, C4s, dC4s_C3s_O3s, O4s, O3_SP_SOUTH_RNA, 0);
                }
            }
        } else {
            double[] C1sXYZ = new double[3];
            double[] C3sXYZ;
            double[] O4sXYZ = new double[3];
            double[] C4sXYZ = new double[3];
            C1s.getXYZ(C1sXYZ);
            O4s.getXYZ(O4sXYZ);
            C4s.getXYZ(C4sXYZ);

            // O3s coordinates will be filled into ret.
            if (isDeoxy) {
                if (pucker == NucleicSugarPucker.C3_ENDO) {
                    C3sXYZ = determineIntxyz(C4sXYZ, dC3s_C4s, O4sXYZ, dC3s_C4s_O4s, C1sXYZ, C3_SP_NORTH_DNA, 0);
                    ret = determineIntxyz(C3sXYZ, dC3s_O3s, C4sXYZ, dC4s_C3s_O3s, O4sXYZ, O3_SP_NORTH_DNA, 0);
                } // TODO: else-if for 3'-exo configuration (DNA only)
                else {
                    C3sXYZ = determineIntxyz(C4sXYZ, dC3s_C4s, O4sXYZ, dC3s_C4s_O4s, C1sXYZ, C3_SP_SOUTH_DNA, 0);
                    ret = determineIntxyz(C3sXYZ, dC3s_O3s, C4sXYZ, dC4s_C3s_O3s, O4sXYZ, O3_SP_SOUTH_DNA, 0);
                }
            } else {
                if (pucker == NucleicSugarPucker.C3_ENDO) {
                    C3sXYZ = determineIntxyz(C4sXYZ, dC3s_C4s, O4sXYZ, dC3s_C4s_O4s, C1sXYZ, C3_SP_NORTH_RNA, 0);
                    ret = determineIntxyz(C3sXYZ, dC3s_O3s, C4sXYZ, dC4s_C3s_O3s, O4sXYZ, O3_SP_NORTH_RNA, 0);
                } else {
                    C3sXYZ = determineIntxyz(C4sXYZ, dC3s_C4s, O4sXYZ, dC3s_C4s_O4s, C1sXYZ, C3_SP_SOUTH_RNA, 0);
                    ret = determineIntxyz(C3sXYZ, dC3s_O3s, C4sXYZ, dC4s_C3s_O3s, O4sXYZ, O3_SP_SOUTH_RNA, 0);
                }
            }
        }
        return ret;
    }

    /**
     * Draws nucleic acid Atoms outside the backbone. Called after corrections
     * have been applied, so that these Atoms are drawn with ideal bond lengths
     * and angles.
     *
     * @param residue Residue.
     * @param rotamer If 5' capped by HO5s, uses chi5 to draw HO5s.
     * @param prevResidue NA residue at the 5' end of residue.
     * @param isDeoxy If Residue is DNA; false means RNA.
     * @param is3sTerminal If Residue is at a 3' end.
     * @param prevSugarPucker Sugar pucker for prevResidue specified by Rotamer.
     */
    private static void applyNASideAtoms(Residue residue, Rotamer rotamer, Residue prevResidue, boolean isDeoxy, boolean is3sTerminal, NucleicSugarPucker prevSugarPucker) {
        Atom C1s = (Atom) residue.getAtomNode("C1\'");
        Atom C2s = (Atom) residue.getAtomNode("C2\'");
        Atom C3s = (Atom) residue.getAtomNode("C3\'");
        Atom C4s = (Atom) residue.getAtomNode("C4\'");
        Atom O4s = (Atom) residue.getAtomNode("O4\'");
        Atom O3s = (Atom) residue.getAtomNode("O3\'");
        // O2s will be null in DNA.
        Atom O2s = (Atom) residue.getAtomNode("O2\'");

        // Hydrogens attached to the sugar
        // Hydrogens in DNA.  Will be null in RNA.
        Atom H2ss = (Atom) residue.getAtomNode("H2\'\'");
        // Hydrogens in RNA.  Will be null in DNA.
        Atom H2s = (Atom) residue.getAtomNode("H2\'");
        Atom HO2s = (Atom) residue.getAtomNode("HO2\'");
        // Common hydrogens
        Atom H3s = (Atom) residue.getAtomNode("H3\'");
        Atom H4s = (Atom) residue.getAtomNode("H4\'");
        Atom H1s = (Atom) residue.getAtomNode("H1\'");
        Atom H5s = (Atom) residue.getAtomNode("H5\'");
        Atom H5ss = (Atom) residue.getAtomNode("H5\'\'");

        Atom C5s = (Atom) residue.getAtomNode("C5\'");
        Atom O5s = (Atom) residue.getAtomNode("O5\'");
        Atom P = (Atom) residue.getAtomNode("P");
        Atom OP1 = (Atom) residue.getAtomNode("OP1");
        Atom OP2 = (Atom) residue.getAtomNode("OP2");

        // Build atachments to C2'.
        if (isDeoxy) {
            Bond C2s_H2s = C2s.getBond(H2s);
            double dC2s_H2s = C2s_H2s.bondType.distance;
            Angle C3s_C2s_H2s = C3s.getAngle(C2s, H2s);
            double dC3s_C2s_H2s = C3s_C2s_H2s.angleType.angle[C3s_C2s_H2s.nh];
            intxyz(H2s, C2s, dC2s_H2s, C3s, dC3s_C2s_H2s, C1s, 109.4, 1);

            Bond C2s_H2ss = C2s.getBond(H2ss);
            double dC2s_H2ss = C2s_H2ss.bondType.distance;
            Angle C3s_C2s_H2ss = C3s.getAngle(C2s, H2ss);
            double dC3s_C2s_H2ss = C3s_C2s_H2ss.angleType.angle[C3s_C2s_H2ss.nh];
            intxyz(H2ss, C2s, dC2s_H2ss, C3s, dC3s_C2s_H2ss, C1s, 109.4, -1);
        } else {
            Bond C2s_H2s = C2s.getBond(H2s);
            double dC2s_H2s = C2s_H2s.bondType.distance;
            Angle C3s_C2s_H2s = C3s.getAngle(C2s, H2s);
            double dC3s_C2s_H2s = C3s_C2s_H2s.angleType.angle[C3s_C2s_H2s.nh];
            intxyz(H2s, C2s, dC2s_H2s, C3s, dC3s_C2s_H2s, C1s, 109.4, -1);

            Bond C2s_O2s = C2s.getBond(O2s);
            double dC2s_O2s = C2s_O2s.bondType.distance;
            Angle C3s_C2s_O2s = C3s.getAngle(C2s, O2s);
            double dC3s_C2s_O2s = C3s_C2s_O2s.angleType.angle[C3s_C2s_O2s.nh];
            intxyz(O2s, C2s, dC2s_O2s, C3s, dC3s_C2s_O2s, C1s, 109.4, 1);

            /*
             * The placement of HO2' may eventually become a rotameric
             * function, but is presently being defaulted to 70 degrees.
             *
             * 70 degrees was just what I got from looking at 3ZD7 in
             * PyMol
             */
            Bond O2s_HO2s = O2s.getBond(HO2s);
            double dO2s_HO2s = O2s_HO2s.bondType.distance;
            Angle C2s_O2s_HO2s = C2s.getAngle(O2s, HO2s);
            double dC2s_O2s_HO2s = C2s_O2s_HO2s.angleType.angle[C2s_O2s_HO2s.nh];
            intxyz(HO2s, O2s, dO2s_HO2s, C2s, dC2s_O2s_HO2s, C1s, 70, 0);
        }

        Bond C1s_H1s = C1s.getBond(H1s);
        double dC1s_H1s = C1s_H1s.bondType.distance;
        Angle C2s_C1s_H1s = C2s.getAngle(C1s, H1s);
        double dC2s_C1s_H1s = C2s_C1s_H1s.angleType.angle[C2s_C1s_H1s.nh];
        intxyz(H1s, C1s, dC1s_H1s, C2s, dC2s_C1s_H1s, O4s, 109.4, 1);

        Bond C3s_H3s = C3s.getBond(H3s);
        double dC3s_H3s = C3s_H3s.bondType.distance;
        Angle C2s_C3s_H3s = C2s.getAngle(C3s, H3s);
        double dC2s_C3s_H3s = C2s_C3s_H3s.angleType.angle[C2s_C3s_H3s.nh];
        intxyz(H3s, C3s, dC3s_H3s, C2s, dC2s_C3s_H3s, C4s, 109.4, 1);

        Bond C4s_H4s = C4s.getBond(H4s);
        double dC4s_H4s = C4s_H4s.bondType.distance;
        Angle C3s_C4s_H4s = C3s.getAngle(C4s, H4s);
        double dC3s_C4s_H4s = C3s_C4s_H4s.angleType.angle[C3s_C4s_H4s.nh];
        intxyz(H4s, C4s, dC4s_H4s, C3s, dC3s_C4s_H4s, O4s, 109.4, -1);

        Bond C5s_H5s = C5s.getBond(H5s);
        double dC5s_H5s = C5s_H5s.bondType.distance;
        Angle C4s_C5s_H5s = C4s.getAngle(C5s, H5s);
        double dC4s_C5s_H5s = C4s_C5s_H5s.angleType.angle[C4s_C5s_H5s.nh];
        intxyz(H5s, C5s, dC5s_H5s, C4s, dC4s_C5s_H5s, O5s, 109.4, -1);

        Bond C5s_H5ss = C5s.getBond(H5ss);
        double dC5s_H5ss = C5s_H5ss.bondType.distance;
        Angle C4s_C5s_H5ss = C4s.getAngle(C5s, H5ss);
        double dC4s_C5s_H5ss = C4s_C5s_H5ss.angleType.angle[C4s_C5s_H5ss.nh];
        intxyz(H5ss, C5s, dC5s_H5ss, C4s, dC4s_C5s_H5ss, O5s, 109.4, 1);

        if (is3sTerminal) {
            // TODO: Determine proper labels for 3' phosphate caps so they may be implemented.
            Atom HO3s = (Atom) residue.getAtomNode("HO3\'");
            // if (HO3s != null) {
            Bond O3s_HO3s = O3s.getBond(HO3s);
            double dO3s_HO3s = O3s_HO3s.bondType.distance;
            Angle C3s_O3s_HO3s = C3s.getAngle(O3s, HO3s);
            double dC3s_O3s_HO3s = C3s_O3s_HO3s.angleType.angle[C3s_O3s_HO3s.nh];
            // HO3s defaulted to the antiperiplanar value of 180 degrees.
            intxyz(HO3s, O3s, dO3s_HO3s, C3s, dC3s_O3s_HO3s, C4s, 180, 0);
        }

        if (P != null) {
            double[] PXYZ = new double[3];
            P.getXYZ(PXYZ);
            double[] O5sXYZ = new double[3];
            O5s.getXYZ(O5sXYZ);
            double[] C5sXYZ = new double[3];
            C5s.getXYZ(C5sXYZ);

            Bond P_OP1 = P.getBond(OP1);
            double dP_OP1 = P_OP1.bondType.distance;
            Angle O5s_P_OP1 = C5s.getAngle(O5s, P);
            double dO5s_P_OP1 = O5s_P_OP1.angleType.angle[O5s_P_OP1.nh];

            Bond P_OP2 = P.getBond(OP2);
            double dP_OP2 = P_OP2.bondType.distance;
            Angle O5s_P_OP2 = C5s.getAngle(O5s, P);
            double dO5s_P_OP2 = O5s_P_OP2.angleType.angle[O5s_P_OP2.nh];

            // TODO: Handle hydrogens attached to 5'-terminal phosphates.

            /*
             * If there is a prior residue, draw tetrahedrally based on O3'
             * (i-1).  Else, draw based on OP3.
             */
            if (prevResidue != null) {
                double[] O3sPriorCoords;
                if (prevSugarPucker == NucleicSugarPucker.C3_ENDO) {
                    O3sPriorCoords = prevResidue.getO3sNorth();
                } else {
                    O3sPriorCoords = prevResidue.getO3sSouth();
                }
                double[] OP1XYZ = determineIntxyz(PXYZ, dP_OP1, O5sXYZ, dO5s_P_OP1, O3sPriorCoords, 109.4, 1);
                double[] OP2XYZ = determineIntxyz(PXYZ, dP_OP2, O5sXYZ, dO5s_P_OP2, O3sPriorCoords, 109.4, -1);
                OP1.moveTo(OP1XYZ);
                OP2.moveTo(OP2XYZ);
            } else {
                Atom OP3 = (Atom) residue.getAtomNode("OP3");
                double[] OP3XYZ = new double[3];
                OP3.getXYZ(OP3XYZ);
                double[] OP1XYZ = determineIntxyz(PXYZ, dP_OP1, O5sXYZ, dO5s_P_OP1, OP3XYZ, 109.4, 1);
                double[] OP2XYZ = determineIntxyz(PXYZ, dP_OP2, O5sXYZ, dO5s_P_OP2, OP3XYZ, 109.4, -1);
                OP1.moveTo(OP1XYZ);
                OP2.moveTo(OP2XYZ);
            }
        } else {
            Atom HO5s = (Atom) residue.getAtomNode("HO5\'");
            Bond O5s_HO5s = O5s.getBond(HO5s);
            double dO5s_HO5s = O5s_HO5s.bondType.distance;
            Angle C5s_O5s_HO5s = C5s.getAngle(O5s, HO5s);
            double dC5s_O5s_HO5s = C5s_O5s_HO5s.angleType.angle[C5s_O5s_HO5s.nh];
            intxyz(HO5s, O5s, dO5s_HO5s, C5s, dC5s_O5s_HO5s, C4s, rotamer.chi5, 0);
        }
    }

    /**
     * Draws the backbone of a nucleic acid Residue.
     *
     * @param residue Residue.
     * @param rotamer Rotamer being applied to Residue.
     * @param prevResidue Residue 5' of residue.
     */
    private static void applyNABackbone(Residue residue, Rotamer rotamer, Residue prevResidue) {
        Atom O3s = (Atom) residue.getAtomNode("O3\'");
        Atom C3s = (Atom) residue.getAtomNode("C3\'");
        Atom C4s = (Atom) residue.getAtomNode("C4\'");
        Atom C5s = (Atom) residue.getAtomNode("C5\'");
        Atom O5s = (Atom) residue.getAtomNode("O5\'");
        /*
         * Two of the following atoms will be null, depending on whether
         * there is a previous residue, whether it is a 5' PO4 cap, or whether
         * it is a 5' OH cap.
         */
        Atom P = (Atom) residue.getAtomNode("P");
        Atom OP3 = (Atom) residue.getAtomNode("OP3");
        Atom HO5s = (Atom) residue.getAtomNode("HO5\'");

        Bond C4s_C5s = C4s.getBond(C5s);
        double dC4s_C5s = C4s_C5s.bondType.distance;
        Angle C3s_C4s_C5s = C3s.getAngle(C4s, C5s);
        double dC3s_C4s_C5s = C3s_C4s_C5s.angleType.angle[C3s_C4s_C5s.nh];
        intxyz(C5s, C4s, dC4s_C5s, C3s, dC3s_C4s_C5s, O3s, rotamer.chi7, 0);

        Bond C5s_O5s = C5s.getBond(O5s);
        double dC5s_O5s = C5s_O5s.bondType.distance;
        Angle C4s_C5s_O5s = C4s.getAngle(C5s, O5s);
        double dC4s_C5s_O5s = C4s_C5s_O5s.angleType.angle[C4s_C5s_O5s.nh];
        intxyz(O5s, C5s, dC5s_O5s, C4s, dC4s_C5s_O5s, C3s, rotamer.chi6, 0);

        if (prevResidue == null) {
            // If capped by HO5s, draw HO5s and return. Else, assume OP3 capped.
            if (HO5s != null) {
                Bond O5s_HO5s = O5s.getBond(HO5s);
                double dO5s_HO5s = O5s_HO5s.bondType.distance;
                Angle C5s_O5s_HO5s = C5s.getAngle(O5s, HO5s);
                double dC5s_O5s_HO5s = C5s_O5s_HO5s.angleType.angle[C5s_O5s_HO5s.nh];
                intxyz(HO5s, O5s, dO5s_HO5s, C5s, dC5s_O5s_HO5s, C4s, rotamer.chi5, 0);
            } else {
                Bond O5s_P = O5s.getBond(P);
                double dO5s_P = O5s_P.bondType.distance;
                Angle C5s_O5s_P = C5s.getAngle(O5s, P);
                double dC5s_O5s_P = C5s_O5s_P.angleType.angle[C5s_O5s_P.nh];
                intxyz(P, O5s, dO5s_P, C5s, dC5s_O5s_P, C4s, rotamer.chi5, 0);

                Bond P_OP3 = P.getBond(OP3);
                double dP_OP3 = P_OP3.bondType.distance;
                Angle O5s_P_OP3 = C5s.getAngle(O5s, P);
                double dO5s_P_OP3 = O5s_P_OP3.angleType.angle[O5s_P_OP3.nh];
                intxyz(OP3, P, dP_OP3, O5s, dO5s_P_OP3, C5s, rotamer.chi4, 0);
            }
        } else {
            Bond O5s_P = O5s.getBond(P);
            double dO5s_P = O5s_P.bondType.distance;
            Angle C5s_O5s_P = C5s.getAngle(O5s, P);
            double dC5s_O5s_P = C5s_O5s_P.angleType.angle[C5s_O5s_P.nh];
            intxyz(P, O5s, dO5s_P, C5s, dC5s_O5s_P, C4s, rotamer.chi5, 0);
        }
    }

    /**
     * Applies Cartesian translations to nucleic acid backbone atoms to allow P
     * to correctly join up with O3' of the prior Residue.
     *
     * @param residue Residue.
     * @param prevResidue Residue 5' of residue.
     * @param rotamer Rotamer being applied to residue.
     * @param prevSugarPucker Expected sugar pucker of prevResidue.
     * @return The magnitude of any applied correction.
     */
    private static double applyNACorrections(Residue residue, Residue prevResidue, Rotamer rotamer,
            NucleicSugarPucker prevSugarPucker, boolean isDeoxy, boolean is3sTerminal) {
        // Backbone atoms of this residue to be adjusted
        Atom C3s = (Atom) residue.getAtomNode("C3\'");
        Atom O4s = (Atom) residue.getAtomNode("O4\'");
        Atom C4s = (Atom) residue.getAtomNode("C4\'");
        Atom C5s = (Atom) residue.getAtomNode("C5\'");
        Atom O5s = (Atom) residue.getAtomNode("O5\'");
        Atom P = (Atom) residue.getAtomNode("P");
        Atom C1s = (Atom) residue.getAtomNode("C1\'");
        Atom C2s = (Atom) residue.getAtomNode("C2\'");

        // This reference being used solely to get ideal bond lengths & angles.
        Atom O3sPrev = (Atom) prevResidue.getAtomNode("O3\'");
        double[] O3sPriorCoords;

        // Original position of O3' (i-1). Will be used to draw the correction
        // vector.
        if (prevSugarPucker == NucleicSugarPucker.C3_ENDO) {
            O3sPriorCoords = prevResidue.getO3sNorth();
        } else {
            O3sPriorCoords = prevResidue.getO3sSouth();
        } // TODO: Else-if block for the C3'-exo configuration of DNA sugars.

        Bond P_O3sPrev = P.getBond(O3sPrev);
        double dP_O3sPrev = P_O3sPrev.bondType.distance;
        Angle O5s_P_O3sPrev = O5s.getAngle(P, O3sPrev);
        double dO5s_P_O3sPrev = O5s_P_O3sPrev.angleType.angle[O5s_P_O3sPrev.nh];
        double[] O3sHypCoords = determineIntxyz(P.getXYZ(null), dP_O3sPrev, O5s.getXYZ(null), dO5s_P_O3sPrev,
                C5s.getXYZ(null), rotamer.chi4, 0);

        // Index 5 will be full correction, and indices 0-4 will be 1/6 to 5/6
        // of the full correction in increasing order.  Index 6 is a 1/12
        // correction applied to other atoms in the sugar.
        double[][] corrections = new double[7][3];
        for (int i = 0; i < 3; i++) {
            corrections[5][i] = O3sPriorCoords[i] - O3sHypCoords[i];
            corrections[0][i] = (1.0 / 6.0) * corrections[5][i];
            corrections[1][i] = (1.0 / 3.0) * corrections[5][i];
            corrections[2][i] = (1.0 / 2.0) * corrections[5][i];
            corrections[3][i] = (2.0 / 3.0) * corrections[5][i];
            corrections[4][i] = (5.0 / 6.0) * corrections[5][i];
            corrections[6][i] = (1.0 / 12.0) * corrections[5][i];
        }

        /*
         * Move backbone atoms by an appropriate fraction of the correction
         * vector. Do this before checking the threshold, so that atoms are moved
         * in case that is needed before the exception gets thrown.
         */
        O4s.move(corrections[0]);
        C3s.move(corrections[0]);
        C4s.move(corrections[1]);
        C5s.move(corrections[2]);
        O5s.move(corrections[3]);
        P.move(corrections[4]);
        C1s.move(corrections[6]);
        C2s.move(corrections[6]);

        return ((corrections[5][0] * corrections[5][0])
                + (corrections[5][1] * corrections[5][1])
                + (corrections[5][2] * corrections[5][2]));
    }
    
    public static boolean addRotPatch(String rotFileName) {
        File rotPatchFile = new File(rotFileName);
        if (rotPatchFile.exists() && rotPatchFile.canRead()) {
            return addRotPatch(new File(rotFileName));
        }
        return false;
    }
    
    private static boolean addRotPatch(File rpatchFile) {
        try (BufferedReader br = new BufferedReader(new FileReader(rpatchFile))) {
            String resName = null;
            List<String> applyLines = new ArrayList<>();
            //List<Rotamer> rotamers = new ArrayList<>();
            List<String> rotLines = new ArrayList<>();
            ResidueType rType = ResidueType.AA;
            String line = br.readLine();
            while (line != null) {
                line = line.trim();
                if (line.startsWith("PLACE")) {
                    applyLines.add(line);
                } else if (line.startsWith("RESNAME")) {
                    String[] toks = line.split("\\s+");
                    resName = toks[1];
                } else if (line.startsWith("RESTYPE")) {
                    String[] toks = line.split("\\s+");
                    switch (toks[1]) {
                        case "AA":
                            rType = ResidueType.AA;
                            break;
                        case "NA":
                            rType = ResidueType.NA;
                            break;
                        case "UNK":
                        default:
                            rType = ResidueType.UNK;
                            break;
                    }
                } else if (line.startsWith("ROTAMER")) {
                    rotLines.add(line);
                }
                line = br.readLine();
            }
            if (resName != null) {
                List<Rotamer> rotamers = new ArrayList<>();
                for (String string : rotLines) {
                    String[] toks = string.split("\\s+");
                    int nVals = toks.length - 1;
                    double[] values = new double[nVals];
                    for (int i = 0; i < nVals; i++) {
                        values[i] = Double.parseDouble(toks[i + 1]);
                    }
                    switch (rType) {
                        case AA:
                            rotamers.add(new Rotamer(AminoAcid3.UNK, values));
                            break;
                        case NA:
                            rotamers.add(new Rotamer(NucleicAcid3.UNK, values));
                            break;
                        case UNK:
                        default:
                            rotamers.add(new Rotamer(values));
                            break;
                    }
                }
                
                if (nonstdRotCache.containsKey(resName)) {
                    logger.warning(String.format(" Rotamer library already contains "
                            + "rotamer definition for residue %s!", resName));
                } else {
                    NonstandardRotLibrary nrlib = new NonstandardRotLibrary(resName,
                            applyLines.toArray(new String[applyLines.size()]),
                            rotamers.toArray(new Rotamer[rotamers.size()]));
                    nonstdRotCache.put(resName, nrlib);
                }
                return true;
            } else {
                return false;
            }
        } catch (IOException ex) {
            logger.warning(String.format(" Exception in parsing rotamer patch "
                    + "file %s: %s", rpatchFile.getName(), ex.toString()));
            return false;
        }
    }

    public enum ProteinLibrary {

        PonderAndRichards (1), Richardson (2), None(-1);

        private final int oldIntegerConstant;

        ProteinLibrary(int oldConst) {
            this.oldIntegerConstant = oldConst;
        }

        /**
         * Parses a String input to a ProteinLibrary. Can be either a name, or an integer constant.
         *
         * The name is preferred, but the integer constant is allowed for legacy reasons.
         *
         * @param input Input to parse.
         * @return A ProteinLibrary.
         * @throws IllegalArgumentException If no matching ProteinLibrary found.
         */
        public static ProteinLibrary getProteinLibrary(String input) throws IllegalArgumentException {
            if (input.matches("^\\d+$")) {
                return int2Library(Integer.parseInt(input));
            } else {
                return Arrays.stream(ProteinLibrary.values()).
                        filter((ProteinLibrary pl) -> pl.toString().equalsIgnoreCase(input)).
                        findAny().
                        orElseThrow(() -> new IllegalArgumentException(" No protein library found that corresponds to " + input));
            }
        }

        /**
         * Converts an integer to a corresponding ProteinLibrary. Deprecated in favor of using the actual name.
         *
         * @param library Index of the library.
         * @return A ProteinLibrary.
         * @throws IllegalArgumentException If no matching ProteinLibrary found.
         */
        @Deprecated
        public static ProteinLibrary intToProteinLibrary(int library) throws IllegalArgumentException {
            return int2Library(library);
        }

        /**
         * Converts an integer to a corresponding ProteinLibrary. Wrapped by deprecated intToProteinLibrary.
         *
         * @param library Index of the library.
         * @return A ProteinLibrary.
         * @throws IllegalArgumentException If no matching ProteinLibrary found.
         */
        private static ProteinLibrary int2Library(int library) throws IllegalArgumentException {
            for (ProteinLibrary lib : ProteinLibrary.values()) {
                if (library == lib.oldIntegerConstant) {
                    return lib;
                }
            }
            throw new IllegalArgumentException(String.format(" Could not find a " +
                    "protein rotamer library to correspond with %d!", library));
        }
    }

    public enum NucleicAcidLibrary {
        RICHARDSON
    }
    
    /**
     * Class contains rotamer information for a nonstandard amino acid.
     */
    private static class NonstandardRotLibrary {
        private final String resName;
        private final String[] placeRecords; // Must be ordered correctly!
        private final List<Rotamer> stdRotamers; // "Library" rotamers for this residue.
        
        NonstandardRotLibrary(String resname, String[] placeRecords, Rotamer[] stdRotamers) {
            this.resName = resname.toUpperCase();
            this.placeRecords = Arrays.copyOf(placeRecords, placeRecords.length);
            this.stdRotamers = new ArrayList<>(Arrays.asList(stdRotamers));
        }
        
        String getResName() {
            return resName;
        }
        
        Rotamer[] getRotamers() {
            return stdRotamers.toArray(new Rotamer[stdRotamers.size()]);
        }
        
        /**
         * Intended for use when rotamers added separately from the initial
         * definition file.
         * @param newRots 
         */
        void addRotamers(Rotamer[] newRots) {
            stdRotamers.addAll(Arrays.asList(newRots));
        }
        
        void measureNonstdRot(Residue residue, double[] chi, boolean print) {
            for (String placeRec : placeRecords) {
                String[] toks = placeRec.split("\\s+");
                if (toks[0].equalsIgnoreCase("PLACECHI")) {
                    int chiNum = Integer.parseInt(toks[5]) - 1;
                    Atom at1 = (Atom) residue.getAtomNode(toks[1]);
                    Atom at2 = (Atom) residue.getAtomNode(toks[2]);
                    Atom at3 = (Atom) residue.getAtomNode(toks[3]);
                    Atom at4 = (Atom) residue.getAtomNode(toks[4]);
                    Torsion tors = at1.getTorsion(at2, at3, at4);
                    chi[chiNum] = tors.getValue();
                    if (print) {
                        logger.info(tors.toString());
                    }
                }
            }
        }
        
        /**
         * Applies a nonstandard rotamer to an amino acid.
         * @param residue
         * @param rotamer 
         */
        void applyNonstdRotamer(Residue residue, Rotamer rotamer) {
            if (!residue.getName().equalsIgnoreCase(resName)) {
                throw new IllegalArgumentException(String.format(" Residue %s is "
                        + "not of type %s", residue.toString(), resName));
            }
            for (String record : placeRecords) {
                String[] toks = record.split("\\s+");
                Atom at1 = (Atom) residue.getAtomNode(toks[1]);
                Atom at2 = (Atom) residue.getAtomNode(toks[2]);
                Atom at3 = (Atom) residue.getAtomNode(toks[3]);
                Atom at4 = (Atom) residue.getAtomNode(toks[4]);
                Bond b12 = at1.getBond(at2);
                double dbond = b12.bondType.distance;
                Angle a123 = at1.getAngle(at2, at3);
                double dang = a123.angleType.angle[a123.nh];
                double dtors;
                if (toks[0].equalsIgnoreCase("PLACECHI")) {
                    int chiNum = Integer.parseInt(toks[5]);
                    dtors = rotamer.angles[chiNum - 1];
                    /*switch (chiNum) {
                        case 1:
                            dtors = rotamer.chi1;
                            break;
                        case 2:
                            dtors = rotamer.chi2;
                            break;
                        case 3:
                            dtors = rotamer.chi3;
                            break;
                        case 4:
                            dtors = rotamer.chi4;
                            break;
                        case 5:
                            dtors = rotamer.chi5;
                            break;
                        case 6:
                            dtors = rotamer.chi6;
                            break;
                        case 7:
                            dtors = rotamer.chi7;
                            break;
                        default:
                            throw new IllegalArgumentException(" Must be chi 1-7");
                    }*/
                } else {
                    dtors = Double.parseDouble(toks[5]);
                }
                int chirality = Integer.parseInt(toks[6]);
                intxyz(at1, at2, dbond, at3, dang, at4, dtors, chirality);
            }
        }
    }

    public enum NucleicSugarPucker {
        // Used to be 2, 1, and 3 respectively.
        C2_ENDO("south"), C3_ENDO("north"), C3_EXO();
        private final List<String> alternateNames;

        NucleicSugarPucker() {
            alternateNames = Collections.emptyList();
        }

        NucleicSugarPucker(String aName) {
            alternateNames = Collections.singletonList(aName);
        }

        /**
         * Returns the sugar pucker associated with a delta torsion. Currently does
         * not support the C3'-exo DNA-only pucker.
         *
         * @param delta Delta torsion to check
         * @param isDeoxy If DNA (vs. RNA). Presently ignored.
         * @return Pucker
         */
        public static NucleicSugarPucker checkPucker(double delta, boolean isDeoxy) {
        /*
         * Midpoint between North, South is 115 degrees.
         *
         * 0-360: North is 0-115 or 295-360.
         * -180 to 180: North is -65 to 115.
         */
            delta = Crystal.mod(delta, 360.0);
            if (delta <= 115 || delta > 295) {
                return C3_ENDO;
            } else {
                return C2_ENDO;
            } // TODO: Add else-if to handle C3'-exo pucker.
        }
    }
}
