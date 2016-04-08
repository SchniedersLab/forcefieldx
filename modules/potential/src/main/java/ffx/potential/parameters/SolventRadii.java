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
package ffx.potential.parameters;

import java.util.HashMap;
import java.util.logging.Logger;

/**
 * Provides Bondi scaling factors for atomic radii in Generalized Kirkwood
 * continuum solvent.
 * @author slucore
 */
public class SolventRadii {
    
    private static final Logger logger = Logger.getLogger(ForceField.class.getName());
    
    private final String forcefield;
    private double defaultBondi;
    private double overlapScale;
    private final HashMap<Integer,Double> atomtypeToBondi;
    private final HashMap<Integer,Double> biotypeToBondi;
    
    public SolventRadii(String forcefield) {
        switch (forcefield.toUpperCase()) {
            case "AMOEBA_PROTEIN_2013":
                defaultBondi = 1.15;
                overlapScale = 0.60;
                atomtypeToBondi = amoebapro13ByAtomtype;
                biotypeToBondi = amoebapro13ByBiotype;
                break;
            case "AMBER99SB":
                defaultBondi = 1.25;
                overlapScale = 0.55;
                atomtypeToBondi = amber99sbByAtomtype;
                biotypeToBondi = amber99sbByBiotype;
                break;
            default:
                defaultBondi = 0.0;
                overlapScale = 0.0;
                atomtypeToBondi = null;
                biotypeToBondi = null;
                logger.severe("No GK solvent radii available for forcefield: " + forcefield);
        }
        this.forcefield = forcefield.toUpperCase();
        
        // Parse properties and command-line overrides.
        String shctProp = System.getProperty("gk-overlapScale");
        if (shctProp != null) {
            overlapScale = Double.parseDouble(shctProp);
        }
    }
    
    public double getDefaultBondi() {
        return defaultBondi;
    }
    
    public double getOverlapScale() {
        return overlapScale;
    }
    
    public HashMap<Integer,Double> getAtomBondiMap() {
        return atomtypeToBondi;
    }
    
    public HashMap<Integer,Double> getBioBondiMap() {
        return biotypeToBondi;
    }
    
    /**
     * Map types to bondi factors fit for AMOEBA_PROTEIN_2013.
     */
    private static final HashMap<Integer,Double> amoebapro13ByAtomtype = new HashMap<>();
    private static final HashMap<Integer,Double> amoebapro13ByBiotype = new HashMap<>();
    
    /**
     * Map types to bondi factors fit for AMBER99SB.
     */
    private static final HashMap<Integer,Double> amber99sbByAtomtype = new HashMap<>();
    private static final HashMap<Integer,Double> amber99sbByBiotype = new HashMap<>();

    static {
        // GLY
        amoebapro13ByBiotype.put( 2, 1.15);
        amoebapro13ByBiotype.put( 6, 1.15);
        // ALA
//        biotypeToBondi.put( 8, 1.60 );    // DO MANUALLY
//        biotypeToBondi.put( 12, 1.60 );   // DO MANUALLY
        amoebapro13ByBiotype.put( 13, 1.60 );
        amoebapro13ByBiotype.put( 14, 1.60 );
        //VAL
        amoebapro13ByBiotype.put( 21, 1.40 );
        amoebapro13ByBiotype.put( 22, 1.40 );
        amoebapro13ByBiotype.put( 23, 1.40 );
            amoebapro13ByBiotype.put( 25, 1.40 );//^
        amoebapro13ByBiotype.put( 24, 1.40 );
            amoebapro13ByBiotype.put( 26, 1.40 );//^
        //LEU
        amoebapro13ByBiotype.put( 33, 1.40 );
        amoebapro13ByBiotype.put( 34, 1.40 );
        amoebapro13ByBiotype.put( 35, 1.40 );
        amoebapro13ByBiotype.put( 36, 1.40 );
        amoebapro13ByBiotype.put( 37, 1.40 );
            amoebapro13ByBiotype.put( 39, 1.40 );//^
        amoebapro13ByBiotype.put( 38, 1.40 );
            amoebapro13ByBiotype.put( 40, 1.40 );//^
        //ILE
        amoebapro13ByBiotype.put( 47, 1.40 );
        amoebapro13ByBiotype.put( 48, 1.40 );
        amoebapro13ByBiotype.put( 49, 1.40 );
        amoebapro13ByBiotype.put( 50, 1.40 );
        amoebapro13ByBiotype.put( 51, 1.40 );
        amoebapro13ByBiotype.put( 52, 1.40 );
        amoebapro13ByBiotype.put( 53, 1.40 );
        amoebapro13ByBiotype.put( 54, 1.40 );
        //SER
        amoebapro13ByBiotype.put( 63, 1.0235 );
        amoebapro13ByBiotype.put( 64, 1.0235 );
        //THR
        amoebapro13ByBiotype.put( 73, 1.25 );
        amoebapro13ByBiotype.put( 74, 1.25 );
        //CYD
//        biotypeToBondi.put( 83, 1.02 );   // DO MANUALLY
//        biotypeToBondi.put( 84, 1.02 );   // DO MANUALLY
        amoebapro13ByBiotype.put( 95, 1.02 );
        amoebapro13ByBiotype.put( 104, 1.02 );
        //CYS
        amoebapro13ByBiotype.put( 85, 1.80 );
        amoebapro13ByBiotype.put( 86, 1.80 );
        //PRO
        amoebapro13ByBiotype.put( 105, 1.05 );
            amoebapro13ByBiotype.put( 644, 1.05 );//^
        amoebapro13ByBiotype.put( 106, 1.05 );
            amoebapro13ByBiotype.put( 645, 1.05 );//^
        amoebapro13ByBiotype.put( 107, 1.05 );
        amoebapro13ByBiotype.put( 108, 1.05 );
        amoebapro13ByBiotype.put( 109, 1.05 );
        amoebapro13ByBiotype.put( 110, 1.05 );
            amoebapro13ByBiotype.put( 648, 1.05 );//^
        amoebapro13ByBiotype.put( 111, 1.05 );
        amoebapro13ByBiotype.put( 112, 1.05 );
        amoebapro13ByBiotype.put( 113, 1.05 );
        amoebapro13ByBiotype.put( 114, 1.05 );
        amoebapro13ByBiotype.put( 115, 1.05 );
        //PHE
        amoebapro13ByBiotype.put( 122, 1.325 );
        amoebapro13ByBiotype.put( 123, 1.325 );
        amoebapro13ByBiotype.put( 124, 1.325 );
        amoebapro13ByBiotype.put( 125, 1.325 );
        amoebapro13ByBiotype.put( 126, 1.325 );
        amoebapro13ByBiotype.put( 127, 1.325 );
        amoebapro13ByBiotype.put( 128, 1.325 );
        amoebapro13ByBiotype.put( 129, 1.325 );
        amoebapro13ByBiotype.put( 130, 1.325 );
        //TYR
        amoebapro13ByBiotype.put( 145, 1.15 );
        amoebapro13ByBiotype.put( 146, 1.15 );
        //TYD
        amoebapro13ByBiotype.put( 161, 0.938563 );
        //TRP
        amoebapro13ByBiotype.put( 168, 1.32475 );
        amoebapro13ByBiotype.put( 169, 1.32475 );
        amoebapro13ByBiotype.put( 170, 1.32475 );
        amoebapro13ByBiotype.put( 171, 1.32475 );
        amoebapro13ByBiotype.put( 172, 1.32475 );
        amoebapro13ByBiotype.put( 173, 1.32475 );
        amoebapro13ByBiotype.put( 174, 1.32475 );
        amoebapro13ByBiotype.put( 175, 1.32475 );
        amoebapro13ByBiotype.put( 176, 1.32475 );
        amoebapro13ByBiotype.put( 177, 1.32475 );
        amoebapro13ByBiotype.put( 178, 1.32475 );
        amoebapro13ByBiotype.put( 179, 1.32475 );
        amoebapro13ByBiotype.put( 180, 1.32475 );
        amoebapro13ByBiotype.put( 181, 1.32475 );
        amoebapro13ByBiotype.put( 182, 1.32475 );
        amoebapro13ByBiotype.put( 183, 1.32475 );
        amoebapro13ByBiotype.put( 184, 1.32475 );
        //HIS
        amoebapro13ByBiotype.put( 194, 1.60 );
        amoebapro13ByBiotype.put( 195, 1.60 );
        amoebapro13ByBiotype.put( 198, 1.60 );
        amoebapro13ByBiotype.put( 199, 1.60 );
        amoebapro13ByBiotype.put( 200, 1.60 );
        amoebapro13ByBiotype.put( 201, 1.60 );
        //HID
        amoebapro13ByBiotype.put( 211, 1.1375 );
        amoebapro13ByBiotype.put( 212, 1.1375 );
        amoebapro13ByBiotype.put( 215, 1.1375 );
        amoebapro13ByBiotype.put( 216, 1.1375 );
        amoebapro13ByBiotype.put( 217, 1.1375 );
        //HIE
        amoebapro13ByBiotype.put( 227, 1.06175 );
        amoebapro13ByBiotype.put( 230, 1.06175 );
        amoebapro13ByBiotype.put( 231, 1.06175 );
        amoebapro13ByBiotype.put( 232, 1.06175 );
        amoebapro13ByBiotype.put( 233, 1.06175 );
        //ASP
        amoebapro13ByBiotype.put( 242, 1.0555 );
        amoebapro13ByBiotype.put( 243, 1.0555 );
        //ASH
        amoebapro13ByBiotype.put( 252, 1.1125 );
        amoebapro13ByBiotype.put( 253, 1.1125 );
        amoebapro13ByBiotype.put( 254, 1.1125 );
        amoebapro13ByBiotype.put( 255, 1.1125 );
        //ASN
        amoebapro13ByBiotype.put( 265, 1.118125 );
        amoebapro13ByBiotype.put( 266, 1.118125 );
        amoebapro13ByBiotype.put( 267, 1.118125 );
        //GLU
        amoebapro13ByBiotype.put( 278, 1.16 );
        amoebapro13ByBiotype.put( 279, 1.16 );
        //GLH
        amoebapro13ByBiotype.put( 290, 1.06 );
        amoebapro13ByBiotype.put( 291, 1.06 );
        amoebapro13ByBiotype.put( 292, 1.06 );
        amoebapro13ByBiotype.put( 293, 1.06 );
        //GLN
        amoebapro13ByBiotype.put( 305, 1.085 );
        amoebapro13ByBiotype.put( 306, 1.085 );
        amoebapro13ByBiotype.put( 307, 1.085 );
        //MET
        amoebapro13ByBiotype.put( 318, 1.30 );
        //LYS
        amoebapro13ByBiotype.put( 335, 1.64 );
        amoebapro13ByBiotype.put( 336, 1.64 );
        //LYD
        amoebapro13ByBiotype.put( 351, 1.562 );
        amoebapro13ByBiotype.put( 352, 1.562 );
        //ARG
        amoebapro13ByBiotype.put( 365, 1.525 );
        amoebapro13ByBiotype.put( 366, 1.525 );
        amoebapro13ByBiotype.put( 367, 1.525 );
        amoebapro13ByBiotype.put( 368, 1.525 );
        amoebapro13ByBiotype.put( 369, 1.525 );
    }
    
    static {
        //ALA
        amber99sbByBiotype.put(8, 1.55);
        amber99sbByBiotype.put(12, 1.55);
        amber99sbByBiotype.put(13, 1.55);
        amber99sbByBiotype.put(14, 1.55);
        //ARG
        amber99sbByBiotype.put(365, 1.50);
        amber99sbByBiotype.put(366, 1.50);
        amber99sbByBiotype.put(367, 1.50);
        amber99sbByBiotype.put(368, 1.50);
        amber99sbByBiotype.put(369, 1.50);
        //ASN
        amber99sbByBiotype.put(265, 1.160);
        amber99sbByBiotype.put(266, 1.160);
        amber99sbByBiotype.put(267, 1.160);
        //ASP
        amber99sbByBiotype.put(242, .934);
        amber99sbByBiotype.put(243, .934);
        //CYS
        amber99sbByBiotype.put(85, 1.125);
        amber99sbByBiotype.put(86, 1.125);
        //GLN
        amber99sbByBiotype.put(305, 1.1625);
        amber99sbByBiotype.put(306, 1.1625);
        amber99sbByBiotype.put(307, 1.1625);
        //GLU
        amber99sbByBiotype.put(278, 0.92721875);
        amber99sbByBiotype.put(279, 0.92721875);
        //GLY
        amber99sbByBiotype.put(2, 1.50);
        amber99sbByBiotype.put(6, 1.50);
        //HIS
        amber99sbByBiotype.put(194, 1.912);
        amber99sbByBiotype.put(195, 1.912);
        amber99sbByBiotype.put(198, 1.912);
        amber99sbByBiotype.put(199, 1.912);
        amber99sbByBiotype.put(200, 1.912);
        amber99sbByBiotype.put(201, 1.912);
        //ILE
        amber99sbByBiotype.put(47, 1.50);
        amber99sbByBiotype.put(48, 1.50);
        amber99sbByBiotype.put(49, 1.50);
        amber99sbByBiotype.put(50, 1.50);
        amber99sbByBiotype.put(51, 1.50);
        amber99sbByBiotype.put(52, 1.50);
        amber99sbByBiotype.put(53, 1.50);
        amber99sbByBiotype.put(54, 1.50);
        //LEU
        amber99sbByBiotype.put(33, 1.50);
        amber99sbByBiotype.put(34, 1.50);
        amber99sbByBiotype.put(35, 1.50);
        amber99sbByBiotype.put(36, 1.50);
        amber99sbByBiotype.put(37, 1.50);
        amber99sbByBiotype.put(38, 1.50);
        //LYS
        amber99sbByBiotype.put(335, 1.50);
        amber99sbByBiotype.put(336, 1.50);
        //MET
        amber99sbByBiotype.put(318, 2.05);
        //PHE
        amber99sbByBiotype.put(122, 1.1875);
        amber99sbByBiotype.put(123, 1.1875);
        amber99sbByBiotype.put(124, 1.1875);
        amber99sbByBiotype.put(125, 1.1875);
        amber99sbByBiotype.put(126, 1.1875);
        amber99sbByBiotype.put(127, 1.1875);
        amber99sbByBiotype.put(128, 1.1875);
        amber99sbByBiotype.put(129, 1.1875);
        amber99sbByBiotype.put(130, 1.1875);
        //PRO
        amber99sbByBiotype.put(105, 1.50);
        amber99sbByBiotype.put(106, 1.50);
        amber99sbByBiotype.put(107, 1.50);
        amber99sbByBiotype.put(108, 1.50);
        amber99sbByBiotype.put(109, 1.50);
        amber99sbByBiotype.put(110, 1.50);
        amber99sbByBiotype.put(111, 1.50);
        amber99sbByBiotype.put(112, 1.50);
        amber99sbByBiotype.put(113, 1.50);
        amber99sbByBiotype.put(114, 1.50);
        amber99sbByBiotype.put(115, 1.50);
        //SER
        amber99sbByBiotype.put(63, 1.174);
        amber99sbByBiotype.put(64, 1.174);
        //THR
        amber99sbByBiotype.put(73, 1.11175);
        amber99sbByBiotype.put(74, 1.11175);
        //tmp
        //TRP
        amber99sbByBiotype.put(168, 1.14375);
        amber99sbByBiotype.put(169, 1.14375);
        amber99sbByBiotype.put(170, 1.14375);
        amber99sbByBiotype.put(171, 1.14375);
        amber99sbByBiotype.put(172, 1.14375);
        amber99sbByBiotype.put(173, 1.14375);
        amber99sbByBiotype.put(174, 1.14375);
        amber99sbByBiotype.put(175, 1.14375);
        amber99sbByBiotype.put(176, 1.14375);
        amber99sbByBiotype.put(177, 1.14375);
        amber99sbByBiotype.put(178, 1.14375);
        amber99sbByBiotype.put(179, 1.14375);
        amber99sbByBiotype.put(180, 1.14375);
        amber99sbByBiotype.put(181, 1.14375);
        amber99sbByBiotype.put(182, 1.14375);
        amber99sbByBiotype.put(183, 1.14375);
        amber99sbByBiotype.put(184, 1.14375);
        //TYR
        amber99sbByBiotype.put(145, 1.12425);
        amber99sbByBiotype.put(146, 1.12425);
        //VAL
        amber99sbByBiotype.put(21, 1.15);
        amber99sbByBiotype.put(22, 1.15);
        amber99sbByBiotype.put(23, 1.15);
        amber99sbByBiotype.put(24, 1.15);
        amber99sbByBiotype.put(25, 1.15);
        amber99sbByBiotype.put(26, 1.15);
    }
    
    static {
        // GLY
//        typeToBondi.put(2,1.15);      // unnecessary and may suffer from the 
//        typeToBondi.put(6,1.15);      // problem of atomTypes 8,12
        // ALA
//        typeToBondi.put(8,1.60);    // lots of AAs use this!
//        typeToBondi.put(12,1.60);   // lots of AAs use this!
        amoebapro13ByAtomtype.put(13,1.60);
        amoebapro13ByAtomtype.put(14,1.60);
        // VAL
        amoebapro13ByAtomtype.put(15,1.40);
        amoebapro13ByAtomtype.put(16,1.40);
        amoebapro13ByAtomtype.put(17,1.40);
        amoebapro13ByAtomtype.put(18,1.40);
        // LEU
        amoebapro13ByAtomtype.put(19,1.40);
        amoebapro13ByAtomtype.put(20,1.40);
        amoebapro13ByAtomtype.put(21,1.40);
        amoebapro13ByAtomtype.put(22,1.40);
        amoebapro13ByAtomtype.put(23,1.40);
        amoebapro13ByAtomtype.put(24,1.40);
        // ILE
        amoebapro13ByAtomtype.put(25,1.40);
        amoebapro13ByAtomtype.put(26,1.40);
        amoebapro13ByAtomtype.put(27,1.40);
        amoebapro13ByAtomtype.put(28,1.40);
        amoebapro13ByAtomtype.put(29,1.40);
        amoebapro13ByAtomtype.put(30,1.40);
        amoebapro13ByAtomtype.put(31,1.40);
        amoebapro13ByAtomtype.put(32,1.40);
        // SER
        amoebapro13ByAtomtype.put(35,1.0235);
        amoebapro13ByAtomtype.put(36,1.0235);
        // THR
        amoebapro13ByAtomtype.put(39,1.25);
        amoebapro13ByAtomtype.put(40,1.25);
        // CYD
//        typeToBondi.put(43,1.02);     // shared with CYS!
//        typeToBondi.put(44,1.02);     // shared with CYS!
        amoebapro13ByAtomtype.put(48,1.02);
        amoebapro13ByAtomtype.put(49,1.02);
        // CYS
        amoebapro13ByAtomtype.put(45,1.80);
        amoebapro13ByAtomtype.put(46,1.80);
        // PRO
        amoebapro13ByAtomtype.put(50,1.05);
        amoebapro13ByAtomtype.put(51,1.05);
        amoebapro13ByAtomtype.put(52,1.05);
        amoebapro13ByAtomtype.put(53,1.05);
        amoebapro13ByAtomtype.put(54,1.05);
        amoebapro13ByAtomtype.put(55,1.05);
        amoebapro13ByAtomtype.put(56,1.05);
        amoebapro13ByAtomtype.put(57,1.05);
        amoebapro13ByAtomtype.put(58,1.05);
        amoebapro13ByAtomtype.put(59,1.05);
        amoebapro13ByAtomtype.put(60,1.05);
        // PHE
        amoebapro13ByAtomtype.put(61,1.325);
        amoebapro13ByAtomtype.put(62,1.325);
        amoebapro13ByAtomtype.put(63,1.325);
        amoebapro13ByAtomtype.put(64,1.325);
        amoebapro13ByAtomtype.put(65,1.325);
        amoebapro13ByAtomtype.put(66,1.325);
        amoebapro13ByAtomtype.put(67,1.325);
        amoebapro13ByAtomtype.put(68,1.325);
        amoebapro13ByAtomtype.put(69,1.325);
        // TYR
        amoebapro13ByAtomtype.put(78,1.15);
        amoebapro13ByAtomtype.put(79,1.15);
        // TYD
        amoebapro13ByAtomtype.put(88,0.938563);
        // TRP
        amoebapro13ByAtomtype.put(89,1.32475);
        amoebapro13ByAtomtype.put(90,1.32475);
        amoebapro13ByAtomtype.put(91,1.32475);
        amoebapro13ByAtomtype.put(92,1.32475);
        amoebapro13ByAtomtype.put(93,1.32475);
        amoebapro13ByAtomtype.put(94,1.32475);
        amoebapro13ByAtomtype.put(95,1.32475);
        amoebapro13ByAtomtype.put(96,1.32475);
        amoebapro13ByAtomtype.put(97,1.32475);
        amoebapro13ByAtomtype.put(98,1.32475);
        amoebapro13ByAtomtype.put(99,1.32475);
        amoebapro13ByAtomtype.put(100,1.32475);
        amoebapro13ByAtomtype.put(101,1.32475);
        amoebapro13ByAtomtype.put(102,1.32475);
        amoebapro13ByAtomtype.put(103,1.32475);
        amoebapro13ByAtomtype.put(104,1.32475);
        amoebapro13ByAtomtype.put(105,1.32475);
        // HIS
        amoebapro13ByAtomtype.put(109,1.60);
        amoebapro13ByAtomtype.put(110,1.60);
        amoebapro13ByAtomtype.put(113,1.60);
        amoebapro13ByAtomtype.put(114,1.60);
        amoebapro13ByAtomtype.put(115,1.60);
        amoebapro13ByAtomtype.put(116,1.60);
        // HID
        amoebapro13ByAtomtype.put(120,1.1375);
        amoebapro13ByAtomtype.put(121,1.1375);
        amoebapro13ByAtomtype.put(124,1.1375);
        amoebapro13ByAtomtype.put(125,1.1375);
        amoebapro13ByAtomtype.put(126,1.1375);
        // HIE
        amoebapro13ByAtomtype.put(130,1.06175);
        amoebapro13ByAtomtype.put(133,1.06175);
        amoebapro13ByAtomtype.put(134,1.06175);
        amoebapro13ByAtomtype.put(135,1.06175);
        amoebapro13ByAtomtype.put(136,1.06175);
        // ASP
        amoebapro13ByAtomtype.put(139,1.0555);
        amoebapro13ByAtomtype.put(140,1.0555);
        // ASH
        amoebapro13ByAtomtype.put(143,1.1125);
        amoebapro13ByAtomtype.put(144,1.1125);
        amoebapro13ByAtomtype.put(145,1.1125);
        amoebapro13ByAtomtype.put(146,1.1125);
        // ASN
        amoebapro13ByAtomtype.put(150,1.118125);
        amoebapro13ByAtomtype.put(151,1.118125);
        amoebapro13ByAtomtype.put(152,1.118125);
        // GLU
        amoebapro13ByAtomtype.put(157,1.16);
        amoebapro13ByAtomtype.put(158,1.16);
        // GLH
        amoebapro13ByAtomtype.put(163,1.06);
        amoebapro13ByAtomtype.put(164,1.06);
        amoebapro13ByAtomtype.put(165,1.06);
        amoebapro13ByAtomtype.put(166,1.06);
        // GLN
        amoebapro13ByAtomtype.put(172,1.085);
        amoebapro13ByAtomtype.put(173,1.085);
        amoebapro13ByAtomtype.put(174,1.085);
        // MET
        amoebapro13ByAtomtype.put(179,1.30);
        // LYS
        amoebapro13ByAtomtype.put(190,1.64);
        amoebapro13ByAtomtype.put(191,1.64);
        // LYD
        amoebapro13ByAtomtype.put(200,1.562);
        amoebapro13ByAtomtype.put(201,1.562);
        // ARG
        amoebapro13ByAtomtype.put(208,1.525);
        amoebapro13ByAtomtype.put(209,1.525);
        amoebapro13ByAtomtype.put(210,1.525);
        amoebapro13ByAtomtype.put(211,1.525);
        amoebapro13ByAtomtype.put(212,1.525);
    }
    
    static {
        // ALA
//        amber99sbByAtomtype.put(8, 1.55);     // lots of AAs use these
//        amber99sbByAtomtype.put(12, 1.55);
        amber99sbByAtomtype.put(13, 1.55);
        amber99sbByAtomtype.put(14, 1.55);
        // ARG
        amber99sbByAtomtype.put(299, 1.50);
        amber99sbByAtomtype.put(300, 1.50);
        amber99sbByAtomtype.put(301, 1.50);
        amber99sbByAtomtype.put(302, 1.50);
        amber99sbByAtomtype.put(303, 1.50);
        // ASN
        amber99sbByAtomtype.put(229, 1.160);
        amber99sbByAtomtype.put(230, 1.160);
        amber99sbByAtomtype.put(231, 1.160);
        // ASP
        amber99sbByAtomtype.put(218, 0.934);
        amber99sbByAtomtype.put(219, 0.934);
        // CYS
        amber99sbByAtomtype.put(85, 1.125);
        amber99sbByAtomtype.put(86, 1.125);
        // GLN
        amber99sbByAtomtype.put(255, 1.1625);
        amber99sbByAtomtype.put(256, 1.1625);
        amber99sbByAtomtype.put(257, 1.1625);
        // GLU
        amber99sbByAtomtype.put(242, 0.92721875);
        amber99sbByAtomtype.put(243, 0.92721875);
        // GLY
//        amber99sbByAtomtype.put(2, 1.50);     // lots of AAs use these
//        amber99sbByAtomtype.put(6, 1.50);
        // HIS
        amber99sbByAtomtype.put(170, 1.912);
        amber99sbByAtomtype.put(171, 1.912);
        amber99sbByAtomtype.put(174, 1.912);
        amber99sbByAtomtype.put(175, 1.912);
        amber99sbByAtomtype.put(176, 1.912);
        amber99sbByAtomtype.put(177, 1.912);
        // ILE
        amber99sbByAtomtype.put(47, 1.50);
        amber99sbByAtomtype.put(48, 1.50);
        amber99sbByAtomtype.put(49, 1.50);
        amber99sbByAtomtype.put(50, 1.50);
        amber99sbByAtomtype.put(51, 1.50);
        amber99sbByAtomtype.put(52, 1.50);
        amber99sbByAtomtype.put(53, 1.50);
        amber99sbByAtomtype.put(54, 1.50);
        // LEU
        amber99sbByAtomtype.put(33, 1.50);
        amber99sbByAtomtype.put(34, 1.50);
        amber99sbByAtomtype.put(35, 1.50);
        amber99sbByAtomtype.put(36, 1.50);
        amber99sbByAtomtype.put(37, 1.50);
        amber99sbByAtomtype.put(38, 1.50);
        // LYS
        amber99sbByAtomtype.put(285, 1.50);
        amber99sbByAtomtype.put(286, 1.50);
        // MET
        amber99sbByAtomtype.put(268, 2.05);
        // PHE
        amber99sbByAtomtype.put(113, 1.1875);
        amber99sbByAtomtype.put(114, 1.1875);
        amber99sbByAtomtype.put(115, 1.1875);
        amber99sbByAtomtype.put(116, 1.1875);
        amber99sbByAtomtype.put(117, 1.1875);
        amber99sbByAtomtype.put(118, 1.1875);
        amber99sbByAtomtype.put(119, 1.1875);
        amber99sbByAtomtype.put(120, 1.1875);
        amber99sbByAtomtype.put(121, 1.1875);
        // PRO
        amber99sbByAtomtype.put(96, 1.50);
        amber99sbByAtomtype.put(97, 1.50);
        amber99sbByAtomtype.put(98, 1.50);
        amber99sbByAtomtype.put(99, 1.50);
        amber99sbByAtomtype.put(100, 1.50);
        amber99sbByAtomtype.put(101, 1.50);
        amber99sbByAtomtype.put(102, 1.50);
        amber99sbByAtomtype.put(103, 1.50);
        amber99sbByAtomtype.put(104, 1.50);
        amber99sbByAtomtype.put(105, 1.50);
        amber99sbByAtomtype.put(106, 1.50);
        // SER
        amber99sbByAtomtype.put(63, 1.174);
        amber99sbByAtomtype.put(64, 1.174);
        // THR
        amber99sbByAtomtype.put(73, 1.11175);
        amber99sbByAtomtype.put(74, 1.11175);
        // TRP
        amber99sbByAtomtype.put(144, 1.14375);
        amber99sbByAtomtype.put(145, 1.14375);
        amber99sbByAtomtype.put(146, 1.14375);
        amber99sbByAtomtype.put(147, 1.14375);
        amber99sbByAtomtype.put(148, 1.14375);
        amber99sbByAtomtype.put(149, 1.14375);
        amber99sbByAtomtype.put(150, 1.14375);
        amber99sbByAtomtype.put(151, 1.14375);
        amber99sbByAtomtype.put(152, 1.14375);
        amber99sbByAtomtype.put(153, 1.14375);
        amber99sbByAtomtype.put(154, 1.14375);
        amber99sbByAtomtype.put(155, 1.14375);
        amber99sbByAtomtype.put(156, 1.14375);
        amber99sbByAtomtype.put(157, 1.14375);
        amber99sbByAtomtype.put(158, 1.14375);
        amber99sbByAtomtype.put(159, 1.14375);
        amber99sbByAtomtype.put(160, 1.14375);
        // TYR
        amber99sbByAtomtype.put(136, 1.12425);
        amber99sbByAtomtype.put(137, 1.12425);
        // VAL
        amber99sbByAtomtype.put(21, 1.15);
        amber99sbByAtomtype.put(22, 1.15);
        amber99sbByAtomtype.put(23, 1.15);
        amber99sbByAtomtype.put(24, 1.15);
        amber99sbByAtomtype.put(25, 1.15);
        amber99sbByAtomtype.put(26, 1.15);
    }
    
}
