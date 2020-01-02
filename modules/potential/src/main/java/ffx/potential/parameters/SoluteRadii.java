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
package ffx.potential.parameters;

import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import static java.lang.String.format;

import ffx.potential.bonded.Atom;

/**
 * Apply Generalized Kirkwood atomic radii.
 *
 * @author Michael J. Schnieders
 */
public class SoluteRadii {

    private static final Logger logger = Logger.getLogger(SoluteRadii.class.getName());

    private SoluteRadii() {
    }

    /**
     * This maps atomic number to reasonable GK base radii.
     */
    private static final Map<Integer, Double> DEFAULT_RADII = new HashMap<>();

    /**
     * This map connects AMOEBA '09 atom classes to GK base radii.
     */
    private static final Map<Integer, Double> AMOEBA_2009_GK_RADII = new HashMap<>();

    /**
     * This map connects AMOEBA '14 atom classes to GK base radii.
     */
    private static final Map<Integer, Double> AMOEBA_2014_GK_RADII = new HashMap<>();

    /**
     * This map connects AMOEBA Nucleic Acid '17 atom classes to GK base radii.
     */
    private static final Map<Integer, Double> AMOEBA_NUC_2017_GK_RADII = new HashMap<>();

    /**
     * This map connects AMOEBA Bio '18 atom classes to GK base radii.
     */
    private static final Map<Integer, Double> AMOEBA_BIO_2018_GK_RADII = new HashMap<>();


    public static void logRadiiSource(ForceField forceField) {
        String forcefieldName = forceField.getString("FORCEFIELD", ForceField.ForceFieldName.AMOEBA_BIO_2009.toString());
        forcefieldName = forcefieldName.replaceAll("_", "-");
        if (forceField.getBoolean("GK_USEFITRADII", true)) {
            if (forcefieldName.equalsIgnoreCase("AMOEBA-2009")) {
                logger.info(format("   Radii:                  %20s", forcefieldName.toUpperCase()));
            } else if (forcefieldName.equalsIgnoreCase("AMOEBA-2014")) {
                logger.info(format("   Radii:                  %20s", forcefieldName.toUpperCase()));
            } else if (forcefieldName.equalsIgnoreCase("AMOEBA-NUC-2017")) {
                logger.info(format("   Radii:                  %20s", forcefieldName.toUpperCase()));
            } else if (forcefieldName.equalsIgnoreCase("AMOEBA-BIO-2018")) {
                logger.info(format("   Radii:                  %20s", forcefieldName.toUpperCase()));
            }
        }
    }

    public static double applyGKRadii(ForceField forceField, double bondiScale,
                                      Atom[] atoms, double[] baseRadius) {
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            baseRadius[i] = DEFAULT_RADII.get(atom.getAtomicNumber()) * bondiScale;
            int key = atom.getAtomType().atomClass;
            SoluteType soluteType = forceField.getSoluteType(Integer.toString(key));
            if (soluteType != null) {
                baseRadius[i] = soluteType.diameter * 0.5;
            }
        }
        return bondiScale;
    }

    // All solvation free energies are for 1 M standard state in both vacuum and solvent.
    // No "phase potential" correction is applied for moving ions across the vacuum - liquid interface.

    /**
     * Lithium Ion Li+
     * Wang Thesis: -116.91 kca/mol
     * <p>
     * AMOEBA 2009 Class 6
     * <p>
     * Generalized Kirkwood   -119.29504846
     * Cavitation                3.34060579
     * Dispersion               -0.93488991
     * Solvation              -116.88933258
     */
    private static final double GK_AMOEBA_LITHIUM = 2.748;
    /**
     * Sodium Ion Na+
     * Wang Thesis:  -91.62 kcal/mol
     * Grossfield et al: -91.8 kcal/mol
     * <p>
     * AMOEBA 2009 Class 7
     * <p>
     * Generalized Kirkwood    -94.17489031
     * Cavitation                4.64284122
     * Dispersion               -2.05069286
     * Solvation               -91.58274195
     */
    private static final double GK_AMOEBA_SODIUM = 3.481;
    /**
     * Potassium Ion K+
     * Wang Thesis:  -74.32 kca/mol
     * Grossfield et al: -74.5 kca/mol
     * <p>
     * AMOEBA 2009 Class 8
     * <p>
     * Generalized Kirkwood    -77.49947829
     * Cavitation                6.17608374
     * Dispersion               -3.06353750
     * Solvation               -74.38693205
     */
    private static final double GK_AMOEBA_POTASSIUM = 4.230;
    /**
     * Rubidium Ion Rb+
     * Wang Thesis:  -69.05 kcal/mol
     * <p>
     * AMOEBA 2009 Class 9
     * <p>
     * Generalized Kirkwood    -72.20766370
     * Cavitation                6.70040723
     * Dispersion               -3.53557505
     * Solvation               -69.04283152
     */
    private static final double GK_AMOEBA_RUBIDIUM = 4.540;
    /**
     * Cesium Ion Cs+
     * Wang Thesis:  -63.66 kcal/mol
     * <p>
     * AMOEBA 2009 Class 10
     * <p>
     * Generalized Kirkwood    -66.79356014
     * Cavitation                7.30534462
     * Dispersion               -4.16432764
     * Solvation               -63.65254316
     */
    private static final double GK_AMOEBA_CESIUM = 4.908;
    /**
     * Magnesium Ion Mg+2
     * Ren et al. 2006: -431.1 kcal/mol
     * Expt: −435.4 kcal/mol
     * <p>
     * AMOEBA 2009 Class 11
     * <p>
     * Generalized Kirkwood   -434.20237507
     * Cavitation                4.61414342
     * Dispersion               -2.08752589
     * Solvation              -431.67575754
     */
    private static final double GK_AMOEBA_MAGNESIUM = 3.020;
    /**
     * Calcium Ion Ca+2
     * Ren et al. 2006: -354.9 kcal/mol
     * Expt: −357.2 kcal/mol
     * <p>
     * AMOEBA 2009 Class 12
     * <p>
     * Generalized Kirkwood   -357.98284813
     * Cavitation                6.06086494
     * Dispersion               -2.98489939
     * Solvation              -354.90688258
     */
    private static final double GK_AMOEBA_CALCIUM = 3.663;
    /**
     * Zinc Ion Zn+2
     * Ren et al. 2011: −458.9
     * Expt: −467.7 kcal/mol
     * <p>
     * AMOEBA 2009 Class 13
     * <p>
     * Generalized Kirkwood   -461.55972288
     * Cavitation                4.13533119
     * Dispersion               -1.73074989
     * Solvation              -459.15514158
     */
    private static final double GK_AMOEBA_ZINC = 2.841;
    /**
     * Fluoride Ion F-
     * Wang Thesis: -116.71 kcal/mol
     * <p>
     * AMOEBA 2009 Class 14
     * <p>
     * Generalized Kirkwood   -119.90592289
     * Cavitation                5.61431867
     * Dispersion               -2.41793122
     * Solvation              -116.70953544
     */
    private static final double GK_AMOEBA_FLUORIDE = 2.734;
    /**
     * Chloride Ion Cl-
     *
     * Target Data:
     * Grossfield et al: -86.5 kca/mol (from previous AMOEBA Chloride parameters).
     * Wang Thesis:  -86.12 kca/mol (from 2018).
     * Nonpolar solvation: 2.87 kcal/mol

     * <p>
     * AMOEBA 2009 Class 15
     * <p>
     * Generalized Kirkwood    -89.54460344
     * Cavitation                7.25565967
     * Dispersion               -3.86494289 (Nonpolar total: 3.39)
     * Solvation               -86.15388666
     */
    private static final double GK_AMOEBA_CHLORIDE = 3.661;
    /**
     * Bromide Ion Br-
     * Wang Thesis:  -79.66 kca/mol
     * <p>
     * AMOEBA 2009 Class 16
     * <p>
     * Generalized Kirkwood    -82.80444384
     * Cavitation                7.78219576
     * Dispersion               -4.64697404
     * Solvation               -79.66922212
     */
    private static final double GK_AMOEBA_BROMIDE = 3.959;
    /**
     * Iodide Ion I-
     * Wang Thesis:  -71.25 kca/mol
     * <p>
     * AMOEBA 2009 Class 17
     * <p>
     * Generalized Kirkwood    -74.05077777
     * Cavitation                8.59369359
     * Dispersion               -5.79174929
     * Solvation               -71.24883347
     */
    private static final double GK_AMOEBA_IODIDE = 4.427;
    /**
     * Methane CH4 Carbon
     * <p>
     * AMOEBA 2009 Class 25
     */
    private static final double GK_AMOEBA_METHANE_CH4 = 4.914;
    /**
     * Methane H4C Hydrogen
     * <p>
     * AMOEBA 2009 Class 26
     */
    private static final double GK_AMOEBA_METHANE_H4C = 3.770;
    /**
     * Ethane CH3 Carbon
     * Alkane -CH2-
     * Alkane >CH-
     * Alkane >C<
     * <p>
     * AMOEBA 2009 Classes 27, 29, 31 and 33
     */
    private static final double GK_AMOEBA_ALKANE_C = 4.584;
    /**
     * Ethane H3C Hydrogen
     * Alkane -H2C-
     * Alkane -HC<
     * <p>
     * AMOEBA 2009 Class 28, 30 and 32
     */
    private static final double GK_AMOEBA_ALKANE_H = 3.576;
    /**
     * Water Oxygen
     * <p>
     * AMOEBA 2009 Class 34
     * <p>
     * TODO: Combine with Alcohol Oxygen?
     */
    private static final double GK_AMOEBA_WATER_O = 2.9964;
    /**
     * Water Hydrogen
     * <p>
     * AMOEBA 2009 Class 35
     */
    private static final double GK_AMOEBA_WATER_H = 2.3895;
    /**
     * Methanol O, Ethanol O, isoPropanol O, Methyl Ether O
     * Phenol OH, p-Cresol OH
     * <p>
     * AMOEBA 2009 Class 36
     */
    private static final double GK_AMOEBA_ALCOHOL_O = 2.9964;
    /**
     * Methanol HO, Ethanol HO, isoPropanol HO,
     * Phenol HO, p-Cresol HO
     * <p>
     * AMOEBA 2009 Class 37
     */
    private static final double GK_AMOEBA_ALCOHOL_HO = 2.3364;
    /**
     * Methanol CH3
     * <p>
     * AMOEBA 2009 Class 38
     * <p>
     * TODO: Combine with Ethanol C?
     */
    private static final double GK_AMOEBA_METHANOL_C = 3.384;
    /**
     * Methanol H3C, Ethanol H2C, isoPropanol >CH-
     * <p>
     * AMOEBA 2009 Class 39
     * <p>
     * TODO: Combine with Ethanol HC?
     */
    private static final double GK_AMOEBA_METHANOL_HC = 2.583;
    /**
     * Ethanol Carbon, Propanol Me-CH2, Propanol CH3, isoPropanol CH3
     * Methyl Ether CH3
     * Methyl Amine CH3, Ethyl Amine CH2, Ethyl Amine CH3, Propyl Amine CH3
     * Propyl Amine Me-CH2, Dimethyl Amine CH3, Trimethyl Amine CH3
     * Pyrrolidine C-CH2-C, Pyrrolidine CH2-N
     * NMePyrrolidine CH2-N, NMePyrrolidine CH2<, NMePyrrolidine CH3-N
     * Acetamide CH3, Propamide Me-CH2, Propamide CH3
     * NMeFormamide CH3, NEtFormamide CH2-N
     * <p>
     * AMOEBA 2009 Class 40
     */
    private static final double GK_AMOEBA_POLARGROUP_C = 3.438;
    /**
     * Ethanol Hydrogen, Propanol H3C, isoPropanol H3C
     * Ethyl Amine H3C, Propyl Amine H3C
     * Propamide H3C, NEtFormamide H3C
     * Ethyl Sulfonate H3C, Propyl Sulfonate H3C
     * Ethylbenzene H2C, Ethylbenzene H3C, Toluene H3C
     * 4-Ethylimidazole H2C, 4-Ethylimidazole H3C
     * 3-Ethylindole H2C, 3-Ethylindole H3C
     * <p>
     * AMOEBA 2009 Class 41
     */
    private static final double GK_AMOEBA_POLARGROUP_HC = 2.664;
    /**
     * Propanol Me-CH2,
     * Propyl Amine Me-CH2
     * NMePyrrolidine H3C-N
     * Propamide Me-CH2
     * Propyl Sulfonate Me-CH2
     * <p>
     * AMOEBA 2009 Class 42
     * <p>
     * TODO: Combine with Ethanol HC?
     */
    private static final double GK_AMOEBA_PROPANOL_H = 2.682;
    /**
     * Isopropanol >CH-
     * <p>
     * AMOEBA 2009 Class 43
     */
    private static final double GK_AMOEBA_ISOPROPANOL_C = 3.650;
    /**
     * Methyl Ether H3C
     * <p>
     * AMOEBA 2009 Class 44
     * <p>
     * TODO: Is 3.179 too big?
     */
    private static final double GK_AMOEBA_METHYLETHER_H = 3.179;
    /**
     * Ammonia N
     * <p>
     * AMOEBA 2009 Class 45
     */
    private static final double GK_AMOEBA_AMMONIA_N = 3.710;
    /**
     * Ammonia H
     * Imidazole HN, 4-Ethylimidazole HND, 4-Ethylimidazole HNE
     * Indole HN, 3-Ethylindole HN, 3-Formylindole HN
     * <p>
     * AMOEBA 2009 Class 46
     */
    private static final double GK_AMOEBA_NITROGEN_H = 2.700;
    /**
     * Methyl Amine N, Ethyl Amine N, Dimethyl Amine N, Trimethyl Amine N
     * <p>
     * AMOEBA 2009 Class 49
     */
    private static final double GK_AMOEBA_METHYLAMINE_N = 3.1535;
    /**
     * Methyl Amine H2N, Ethyl Amine H2N, Dimethyl Amine HN
     * <p>
     * AMOEBA 2009 Class 50
     */
    private static final double GK_AMOEBA_METHYLAMINE_HN = 2.430;
    /**
     * Methyl Amine H3C, Ethyl Amine H2C
     * Dimethyl Amine H3C, Trimethyl Amine H3C
     * p-Cresol H3C
     * <p>
     * AMOEBA 2009 Class 51
     */
    private static final double GK_AMOEBA_METHYLAMINE_HC = 2.592;
    /**
     * Pyrrolidine C-CH2-C, Pyrrolidine H2C-N
     * NMePyrrolidine H2C-N, NMePyrrolidine H2C<
     * <p>
     * AMOEBA 2009 Class 52
     * <p>
     * TODO: Combine with Ethanol HC?
     */
    private static final double GK_AMOEBA_PYRROLIDINE_HC = 2.664;
    /**
     * Formamide C=O, Acetamide C=O
     * NMeFormamide C=O, NMeAcetamide C=O
     * DiMeFormamide C=O, DiMeAcetamide C=O
     * <p>
     * AMOEBA 2009 Class 53
     */
    private static final double GK_AMOEBA_AMIDE_CO = 3.438;
    /**
     * Formamide HCO, NMeFormamide HCO, DiMeFormamide HCO
     * <p>
     * AMOEBA 2009 Class 54
     */
    private static final double GK_AMOEBA_AMIDE_HCO = 2.520;
    /**
     * Formamide O, Acetamide O, NMeFormamide O
     * NMeAcetamide O, DiMeFormamide O, DiMeAcetamide O
     * Formic Acid O=C, Acetic Acid O=C, Formaldehyde O=C, Acetaldehyde O=C
     * <p>
     * AMOEBA 2009 Class 55
     */
    private static final double GK_AMOEBA_CARBONYL_O = 3.069;
    /**
     * Formamide N, Acetamide N, NMeFormamide N, NMeAcetamide N,
     * DiMeFormamide N, DiMeAcetamide N
     * <p>
     * AMOEBA 2009 Class 56
     */
    private static final double GK_AMOEBA_AMIDE_N = 3.710;
    /**
     * Formamide H2N, Acetamide H2N, NMeFormamide HN, NMeAcetamide HN
     * <p>
     * AMOEBA 2009 Class 57
     */
    private static final double GK_AMOEBA_AMIDE_HN = 2.331;
    /**
     * Acetamide H3C, NMeAcetamide H3C-C, DiMeAcetamide H3C-C
     * NMeFormamide H3C, NMeAcetamide H3C-N, DiMeFormamide H3C, DiMeAcetamide H3C-N
     * <p>
     * AMOEBA 2009 Classes 58 and 59
     * <p>
     * TODO: Combine with Ethanol HC?
     */
    private static final double GK_AMOEBA_AMIDE_H3C = 2.650;
    /**
     * Formic Acid OH, Acetic Acid OH
     * <p>
     * AMOEBA 2009 Class 60
     */
    private static final double GK_AMOEBA_CARBOXCYLIC_ACID_O = 3.30285;
    /**
     * Formic Acid HO, Acetic Acid HO
     * <p>
     * AMOEBA 2009 Class 61
     */
    private static final double GK_AMOEBA_CARBOXCYLIC_ACID_HO = 2.3895;
    /**
     * Formic Acid C=O, Acetic Acid C=O, Formaldehyde C=O, Acetaldehyde C=O
     * Methyl Sulfide CH3, Dimethyl Sulfide CH3, Dimethyl Disulfide CH3,
     * Ethyl Sulfide CH2, MeEt Sulfide CH3-S, MeEt Sulfide CH2-S
     * <p>
     * AMOEBA 2009 Class 62
     */
    private static final double GK_AMOEBA_CARBONYL_CO = 3.402;
    /**
     * Formic Acid HC=O, Formaldehyde HC=O, Acetaldehyde HC=O
     * 3-Formylindole HC=O
     * <p>
     * AMOEBA 2009 Class 63
     */
    private static final double GK_AMOEBA_FORMYL_HCO = 2.628;
    /**
     * Acetaldehyde H3C, Acetaldehyde H3C
     * <p>
     * AMOEBA 2009 Class 64
     */
    private static final double GK_AMOEBA_ACETYL_HCO = 2.682;
    /**
     * Hydrogen Sulfide S
     * <p>
     * AMOEBA 2009 Class 65
     */
    private static final double GK_AMOEBA_HSULFIDE_S = 4.4055;
    /**
     * Hydrogen Sulfide H
     * <p>
     * AMOEBA 2009 Class 66
     */
    private static final double GK_AMOEBA_HSULFIDE_H = 3.3240;
    /**
     * Methyl Sulfide S, Dimethyl Sulfide S, Dimethyl Disulfide S
     * Ethyl Sulfide S, MeEt Sulfide S
     * <p>
     * AMOEBA 2009 Class 67
     */
    private static final double GK_AMOEBA_SULFIDE_S = 5.3267;
    /**
     * Methyl Sulfide HS, Ethyl Sulfide HS
     * <p>
     * AMOEBA 2009 Class 68
     */
    private static final double GK_AMOEBA_SULFIDE_HS = 3.6841;
    /**
     * Methyl Sulfide H3C, Dimethyl Sulfide H3C, Dimethyl Disulfide H3C
     * Ethyl Sulfide H2C, MeEt Sulfide H3C-S, MeEt Sulfide CH2-S
     * <p>
     * AMOEBA 2009 Class 69
     */
    private static final double GK_AMOEBA_SULFIDE_H3C = 3.8171;
    /**
     * Ethyl Sulfide CH3, MeEt Sulfide CH3-C
     * <p>
     * AMOEBA 2009 Class 70
     */
    private static final double GK_AMOEBA_SULFIDE_CH3 = 5.0806;
    /**
     * Ethyl Sulfide H3C, MeEt Sulfide H3C-C
     * <p>
     * AMOEBA 2009 Class 71
     */
    private static final double GK_AMOEBA_ETHYL_SULFIDE_H3C = 3.9634;
    /**
     * Hydrogen Cyanide CN, Acetonitrile CN
     * <p>
     * AMOEBA 2009 Classes 79 and 82
     * <p>
     * TODO: Cobime with Alkane Carbon?
     */
    private static final double GK_AMOEBA_ETHYL_NITRILE_C = 4.584;
    /**
     * Hydrogen Cyanide CN, Acetonitrile CN
     * <p>
     * AMOEBA 2009 Class 80
     */
    private static final double GK_AMOEBA_ETHYL_NITRILE_N = 3.6516;
    /**
     * Acetonitrile H3C
     * <p>
     * AMOEBA 2009 Class 83
     */
    private static final double GK_AMOEBA_ETHYL_NITRILE_H3C = 3.4920;
    /**
     * Benzene C, Pyridinium C2, Pyridinium C3, Pyridinium C4
     * <p>
     * AMOEBA 2009 Class 87
     */
    private static final double GK_AMOEBA_BENZENE_C = 5.0540;
    /**
     * Benzene HC, Indole HC2 to Indole HC7, 3-Ethylindole HC2, 3-Ethylindole HC3 to 3-Ethylindole HC7
     * 3-Formylindole HC2, 3-Formylindole HC4 to 3-Formylindole HC7
     * Pyridinium H2 to Pyridinium H4
     * <p>
     * AMOEBA 2009 Class 88
     */
    private static final double GK_AMOEBA_BENZENE_H = 3.9634;
    /**
     * Ethylbenzene C2, Ethylbenzene C3, Ethylbenzene C4, Ethylbenzene C1-CH2
     * Phenol C1-OH, Phenol C2 to Phenol C4,
     * Toluene C1-CH3, Toluene C2 to Toluene C4
     * p-Cresol C1-CH3, p-Cresol C2, p-Cresol C3, p-Cresol C4-OH
     * Benzamidine C1-CN2, Benzamidine C2 to Benzamidine C4
     * <p>
     * AMOEBA 2009 Class 89
     */
    private static final double GK_AMOEBA_AROMATIC_C = 3.800;
    /**
     * Ethylbenzene H2 to Ethylbenzene H4
     * Phenol H2 to Phenol H4
     * Toluene H2 to Toluene H4
     * p-Cresol H2, p-Cresol H3
     * Benzamidine H2 to Benzamidine H4
     * <p>
     * AMOEBA 2009 Class 90
     */
    private static final double GK_AMOEBA_AROMATIC_H = 2.980;
    /**
     * Imidazole NH, Imidazole N=C-
     * 4-Ethylimidazole ND, 4-Ethylimidazole NE
     * Indole N
     * 3-Ethylindole N
     * 3-Formylindole N
     * <p>
     * AMOEBA 2009 Class 91
     */
    private static final double GK_AMOEBA_IMIDAZOLE_NH = 3.339;
    /**
     * Imidazole N-C-N, Imidazole C-N=C, Imidazole C-NH-
     * 4-Ethylimidazole CE, 4-Ethylimidazole CD, 4-Ethylimidazole CG
     * <p>
     * AMOEBA 2009 Class 92
     */
    private static final double GK_AMOEBA_IMIDAZOLE_C = 3.591;
    /**
     * Imidazole HC
     * 4-Ethylimidazole HCE, 4-Ethylimidazole HCD
     * <p>
     * AMOEBA 2009 Class 93
     */
    private static final double GK_AMOEBA_IMIDAZOLE_HC = 3.000;

    /**
     * TODO: Fit this so GK matches GXG trimer charging.
     */
    private static final double GK_AMOEBA_PROTEIN_CARBOXYLIC_ACID_CO = 3.100;//GK_AMOEBA_CARBONYL_CO;
    private static final double GK_AMOEBA_PROTEIN_CARBOXYLIC_ACID_O = 2.850;//GK_AMOEBA_CARBOXCYLIC_ACID_O;
    private static final double GK_AMOEBA_PROTEIN_LYSINE_NH = 3.700;//GK_AMOEBA_METHYLAMINE_N;
    private static final double GK_AMOEBA_PROTEIN_LYSINE_HN = 2.500;//GK_AMOEBA_METHYLAMINE_HN;
    private static final double GK_AMOEBA_PROTEIN_ARGININE_CZ = 5.350;//GK_AMOEBA_POLARGROUP_C;
    private static final double GK_AMOEBA_PROTEIN_TERMINAL_NH3 = 6.200;//GK_AMOEBA_IMIDAZOLE_NH;
    private static final double GK_AMOEBA_PROTEIN_TERMINAL_H3N = GK_AMOEBA_NITROGEN_H;
    private static final double GK_AMOEBA_PROTEIN_HISTIDINE_NH = 6.200;//GK_AMOEBA_IMIDAZOLE_NH;


    static {
        DEFAULT_RADII.put(0, 0.0);
        DEFAULT_RADII.put(1, 1.2);
        DEFAULT_RADII.put(2, 1.4);
        DEFAULT_RADII.put(5, 1.8);
        DEFAULT_RADII.put(6, 1.7);
        DEFAULT_RADII.put(7, 1.55);
        DEFAULT_RADII.put(8, 1.52);
        DEFAULT_RADII.put(9, 1.47);
        DEFAULT_RADII.put(10, 1.54);
        DEFAULT_RADII.put(14, 2.1);
        DEFAULT_RADII.put(15, 1.8);
        DEFAULT_RADII.put(16, 1.8);
        DEFAULT_RADII.put(17, 1.75);
        DEFAULT_RADII.put(18, 1.88);
        DEFAULT_RADII.put(34, 1.9);
        DEFAULT_RADII.put(35, 1.85);
        DEFAULT_RADII.put(36, 2.02);
        DEFAULT_RADII.put(53, 1.98);
        DEFAULT_RADII.put(54, 2.16);
        for (int i = 0; i <= 118; i++) {
            if (!DEFAULT_RADII.containsKey(i)) {
                DEFAULT_RADII.put(i, 2.0);
            }
        }

        // *************************************************
        // Initialize AMOEBA 2009 GK Radii.
//        atom          1    1    He    "Helium Atom He"               2     4.003    0
//        atom          2    2    Ne    "Neon Atom Ne"                10    20.179    0
//        atom          3    3    Ar    "Argon Atom Ar"               18    39.948    0
//        atom          4    4    Kr    "Krypton Atom Kr"             36    83.800    0
//        atom          5    5    Xe    "Xenon Atom Xe"               54   131.290    0
        AMOEBA_2009_GK_RADII.put(6, GK_AMOEBA_LITHIUM / 2.0);
        AMOEBA_2009_GK_RADII.put(7, GK_AMOEBA_SODIUM / 2.0);
        AMOEBA_2009_GK_RADII.put(8, GK_AMOEBA_POTASSIUM / 2.0);
        AMOEBA_2009_GK_RADII.put(9, GK_AMOEBA_RUBIDIUM / 2.0);
        AMOEBA_2009_GK_RADII.put(10, GK_AMOEBA_CESIUM / 2.0);
        AMOEBA_2009_GK_RADII.put(11, GK_AMOEBA_MAGNESIUM / 2.0);
        AMOEBA_2009_GK_RADII.put(12, GK_AMOEBA_CALCIUM / 2.0);
        AMOEBA_2009_GK_RADII.put(13, GK_AMOEBA_ZINC / 2.0);
        AMOEBA_2009_GK_RADII.put(14, GK_AMOEBA_FLUORIDE / 2.0);
        AMOEBA_2009_GK_RADII.put(15, GK_AMOEBA_CHLORIDE / 2.0);
        AMOEBA_2009_GK_RADII.put(16, GK_AMOEBA_BROMIDE / 2.0);
        AMOEBA_2009_GK_RADII.put(17, GK_AMOEBA_IODIDE / 2.0);
//        atom         18   18    C     "Cyanide Ion C"                6    12.011    1
//        atom         19   19    N     "Cyanide Ion N"                7    14.007    1
//        atom         20   20    B     "Tetrafluoroborate B"          5    10.810    4
//        atom         21   21    F     "Tetrafluoroborate F"          9    18.998    1
//        atom         22   22    P     "Hexafluorophosphate P"       15    30.974    6
//        atom         23   23    F     "Hexafluorophosphate F"        9    18.998    1
//        atom         24   24    N     "Dinitrogen N2"                7    14.007    1
        AMOEBA_2009_GK_RADII.put(25, GK_AMOEBA_METHANE_CH4 / 2.0);
        AMOEBA_2009_GK_RADII.put(26, GK_AMOEBA_METHANE_H4C / 2.0);
        AMOEBA_2009_GK_RADII.put(27, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(28, GK_AMOEBA_ALKANE_H / 2.0);
        AMOEBA_2009_GK_RADII.put(29, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(30, GK_AMOEBA_ALKANE_H / 2.0);
        AMOEBA_2009_GK_RADII.put(31, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(32, GK_AMOEBA_ALKANE_H / 2.0);
        AMOEBA_2009_GK_RADII.put(33, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(34, GK_AMOEBA_WATER_O / 2.0);
        AMOEBA_2009_GK_RADII.put(35, GK_AMOEBA_WATER_H / 2.0);

        AMOEBA_2009_GK_RADII.put(36, GK_AMOEBA_ALCOHOL_O / 2.0);
        AMOEBA_2009_GK_RADII.put(37, GK_AMOEBA_ALCOHOL_HO / 2.0);
        AMOEBA_2009_GK_RADII.put(38, GK_AMOEBA_METHANOL_C / 2.0);
        AMOEBA_2009_GK_RADII.put(39, GK_AMOEBA_METHANOL_HC / 2.0);
        AMOEBA_2009_GK_RADII.put(40, GK_AMOEBA_POLARGROUP_C / 2.0);
        AMOEBA_2009_GK_RADII.put(41, GK_AMOEBA_POLARGROUP_HC / 2.0);
        AMOEBA_2009_GK_RADII.put(42, GK_AMOEBA_PROPANOL_H / 2.0);
        AMOEBA_2009_GK_RADII.put(43, GK_AMOEBA_ISOPROPANOL_C / 2.0);
        AMOEBA_2009_GK_RADII.put(44, GK_AMOEBA_METHYLETHER_H / 2.0);
        AMOEBA_2009_GK_RADII.put(45, GK_AMOEBA_AMMONIA_N / 2.0);
        AMOEBA_2009_GK_RADII.put(46, GK_AMOEBA_NITROGEN_H / 2.0);
//        atom         63   47    N     "Ammonium Ion N+"              7    14.007    4
//        atom         64   48    H     "Ammonium Ion H4N+"            1     1.008    1
        AMOEBA_2009_GK_RADII.put(49, GK_AMOEBA_METHYLAMINE_N / 2.0);
        AMOEBA_2009_GK_RADII.put(50, GK_AMOEBA_METHYLAMINE_HN / 2.0);
        AMOEBA_2009_GK_RADII.put(51, GK_AMOEBA_METHYLAMINE_HC / 2.0);
        AMOEBA_2009_GK_RADII.put(52, GK_AMOEBA_PYRROLIDINE_HC / 2.0);
        AMOEBA_2009_GK_RADII.put(53, GK_AMOEBA_AMIDE_CO / 2.0);
        AMOEBA_2009_GK_RADII.put(54, GK_AMOEBA_AMIDE_HCO / 2.0);
        AMOEBA_2009_GK_RADII.put(55, GK_AMOEBA_CARBONYL_O / 2.0);
        AMOEBA_2009_GK_RADII.put(56, GK_AMOEBA_AMIDE_N / 2.0);
        AMOEBA_2009_GK_RADII.put(57, GK_AMOEBA_AMIDE_HN / 2.0);
        AMOEBA_2009_GK_RADII.put(58, GK_AMOEBA_AMIDE_H3C / 2.0);
        AMOEBA_2009_GK_RADII.put(59, GK_AMOEBA_AMIDE_H3C / 2.0);
        AMOEBA_2009_GK_RADII.put(60, GK_AMOEBA_CARBOXCYLIC_ACID_O / 2.0);
        AMOEBA_2009_GK_RADII.put(61, GK_AMOEBA_CARBOXCYLIC_ACID_HO / 2.0);
        AMOEBA_2009_GK_RADII.put(62, GK_AMOEBA_CARBONYL_CO / 2.0);
        AMOEBA_2009_GK_RADII.put(63, GK_AMOEBA_FORMYL_HCO / 2.0);
        AMOEBA_2009_GK_RADII.put(64, GK_AMOEBA_ACETYL_HCO / 2.0);
        AMOEBA_2009_GK_RADII.put(65, GK_AMOEBA_HSULFIDE_S / 2.0);
        AMOEBA_2009_GK_RADII.put(66, GK_AMOEBA_HSULFIDE_H / 2.0);
        AMOEBA_2009_GK_RADII.put(67, GK_AMOEBA_SULFIDE_S / 2.0);
        AMOEBA_2009_GK_RADII.put(68, GK_AMOEBA_SULFIDE_HS / 2.0);
        AMOEBA_2009_GK_RADII.put(69, GK_AMOEBA_SULFIDE_H3C / 2.0);
        AMOEBA_2009_GK_RADII.put(70, GK_AMOEBA_SULFIDE_CH3 / 2.0);
        AMOEBA_2009_GK_RADII.put(71, GK_AMOEBA_ETHYL_SULFIDE_H3C / 2.0);
//        atom        189   72    S     "Dimethyl Sulfoxide S=O"      16    32.066    3
//        atom        190   73    O     "Dimethyl Sulfoxide S=O"       8    15.999    1
//        atom        191   74    C     "Dimethyl Sulfoxide CH3"       6    12.011    4
//        atom        192   75    H     "Dimethyl Sulfoxide H3C"       1     1.008    1
//        atom        193   76    S     "Methyl Sulfonate SO3-"       16    32.066    4
//        atom        194   77    O     "Methyl Sulfonate SO3-"        8    15.999    1
//        atom        200   78    H     "Ethyl Sulfonate H2C"          1     1.008    1
        AMOEBA_2009_GK_RADII.put(79, GK_AMOEBA_ETHYL_NITRILE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(80, GK_AMOEBA_ETHYL_NITRILE_N / 2.0);
//        atom        209   81    H     "Hydrogen Cyanide HCN"         1     1.008    1
        AMOEBA_2009_GK_RADII.put(82, GK_AMOEBA_ETHYL_NITRILE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(83, GK_AMOEBA_ETHYL_NITRILE_H3C / 2.0);
//        atom        214   84    C     "Tricyanomethide CN"           6    12.011    2
//        atom        215   85    N     "Tricyanomethide CN"           7    14.007    1
//        atom        216   86    C     "Tricyanomethide >C-"          6    12.011    3
        AMOEBA_2009_GK_RADII.put(87, GK_AMOEBA_BENZENE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(88, GK_AMOEBA_BENZENE_H / 2.0);
        AMOEBA_2009_GK_RADII.put(89, GK_AMOEBA_AROMATIC_C / 2.0);
        AMOEBA_2009_GK_RADII.put(90, GK_AMOEBA_AROMATIC_H / 2.0);
        AMOEBA_2009_GK_RADII.put(91, GK_AMOEBA_IMIDAZOLE_NH / 2.0);
        AMOEBA_2009_GK_RADII.put(92, GK_AMOEBA_IMIDAZOLE_C / 2.0);
        AMOEBA_2009_GK_RADII.put(93, GK_AMOEBA_IMIDAZOLE_HC / 2.0);
//        atom        282   94    C     "Indole C2"                    6    12.011    3
//        atom        311   95    C     "3-Ethylindole CH2"            6    12.011    4
//        atom        331   96    O     "3-Formylindole O=C"           8    15.999    1
//        atom        333   97    N     "Benzamidine N"                7    14.007    3
//        atom        334   98    H     "Benzamidine HN"               1     1.008    1
//        atom        335   99    C     "Benzamidine N-C-N"            6    12.011    3
//        atom        343  100    N     "Pyridinium N"                 7    14.007    3
//        atom        347  101    H     "Pyridinium HN"                1     1.008    1

        // *************************************************
        // Initialize AMOEBA 2014 GK Radii.
        AMOEBA_2014_GK_RADII.put(6, GK_AMOEBA_LITHIUM / 2.0);
        AMOEBA_2014_GK_RADII.put(7, GK_AMOEBA_SODIUM / 2.0);
        AMOEBA_2014_GK_RADII.put(8, GK_AMOEBA_POTASSIUM / 2.0);
        AMOEBA_2014_GK_RADII.put(9, GK_AMOEBA_RUBIDIUM / 2.0);
        AMOEBA_2014_GK_RADII.put(10, GK_AMOEBA_CESIUM / 2.0);
        AMOEBA_2014_GK_RADII.put(11, GK_AMOEBA_MAGNESIUM / 2.0);
        AMOEBA_2014_GK_RADII.put(12, GK_AMOEBA_CALCIUM / 2.0);
        AMOEBA_2014_GK_RADII.put(13, GK_AMOEBA_ZINC / 2.0);
        AMOEBA_2014_GK_RADII.put(14, GK_AMOEBA_FLUORIDE / 2.0);
        AMOEBA_2014_GK_RADII.put(15, GK_AMOEBA_CHLORIDE / 2.0);
        AMOEBA_2014_GK_RADII.put(16, GK_AMOEBA_BROMIDE / 2.0);
        AMOEBA_2014_GK_RADII.put(17, GK_AMOEBA_IODIDE / 2.0);

        AMOEBA_2014_GK_RADII.put(25, GK_AMOEBA_METHANE_CH4 / 2.0);
        AMOEBA_2014_GK_RADII.put(26, GK_AMOEBA_METHANE_H4C / 2.0);
        AMOEBA_2014_GK_RADII.put(27, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(28, GK_AMOEBA_ALKANE_H / 2.0);
        AMOEBA_2014_GK_RADII.put(29, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(30, GK_AMOEBA_ALKANE_H / 2.0);
        AMOEBA_2014_GK_RADII.put(31, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(32, GK_AMOEBA_ALKANE_H / 2.0);
        AMOEBA_2014_GK_RADII.put(33, GK_AMOEBA_ALKANE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(34, GK_AMOEBA_WATER_O / 2.0);
        AMOEBA_2014_GK_RADII.put(35, GK_AMOEBA_WATER_H / 2.0);
        AMOEBA_2014_GK_RADII.put(36, GK_AMOEBA_ALCOHOL_O / 2.0);
        AMOEBA_2014_GK_RADII.put(37, GK_AMOEBA_ALCOHOL_HO / 2.0);
        AMOEBA_2014_GK_RADII.put(38, GK_AMOEBA_METHANOL_C / 2.0);
        AMOEBA_2014_GK_RADII.put(39, GK_AMOEBA_METHANOL_HC / 2.0);
        AMOEBA_2014_GK_RADII.put(40, GK_AMOEBA_POLARGROUP_C / 2.0);
        AMOEBA_2014_GK_RADII.put(41, GK_AMOEBA_POLARGROUP_HC / 2.0);
        AMOEBA_2014_GK_RADII.put(42, GK_AMOEBA_PROPANOL_H / 2.0);
        AMOEBA_2014_GK_RADII.put(43, GK_AMOEBA_ISOPROPANOL_C / 2.0);
        AMOEBA_2014_GK_RADII.put(44, GK_AMOEBA_METHYLETHER_H / 2.0);
        AMOEBA_2014_GK_RADII.put(45, GK_AMOEBA_AMMONIA_N / 2.0);
        AMOEBA_2014_GK_RADII.put(46, GK_AMOEBA_NITROGEN_H / 2.0);
        AMOEBA_2014_GK_RADII.put(49, GK_AMOEBA_METHYLAMINE_N / 2.0);
        AMOEBA_2014_GK_RADII.put(50, GK_AMOEBA_METHYLAMINE_HN / 2.0);
        AMOEBA_2014_GK_RADII.put(51, GK_AMOEBA_METHYLAMINE_HC / 2.0);
        AMOEBA_2014_GK_RADII.put(52, GK_AMOEBA_PYRROLIDINE_HC / 2.0);
        AMOEBA_2014_GK_RADII.put(53, GK_AMOEBA_AMIDE_CO / 2.0);
        AMOEBA_2014_GK_RADII.put(54, GK_AMOEBA_AMIDE_HCO / 2.0);
        AMOEBA_2014_GK_RADII.put(55, GK_AMOEBA_CARBONYL_O / 2.0);
        AMOEBA_2014_GK_RADII.put(56, GK_AMOEBA_AMIDE_N / 2.0);
        AMOEBA_2014_GK_RADII.put(57, GK_AMOEBA_AMIDE_HN / 2.0);
        AMOEBA_2014_GK_RADII.put(58, GK_AMOEBA_AMIDE_H3C / 2.0);
        AMOEBA_2014_GK_RADII.put(59, GK_AMOEBA_AMIDE_H3C / 2.0);
        AMOEBA_2014_GK_RADII.put(60, GK_AMOEBA_CARBOXCYLIC_ACID_O / 2.0);
        AMOEBA_2014_GK_RADII.put(61, GK_AMOEBA_CARBOXCYLIC_ACID_HO / 2.0);
        AMOEBA_2014_GK_RADII.put(62, GK_AMOEBA_CARBONYL_CO / 2.0);
        AMOEBA_2014_GK_RADII.put(63, GK_AMOEBA_FORMYL_HCO / 2.0);
        AMOEBA_2014_GK_RADII.put(64, GK_AMOEBA_ACETYL_HCO / 2.0);
        AMOEBA_2014_GK_RADII.put(65, GK_AMOEBA_HSULFIDE_S / 2.0);
        AMOEBA_2014_GK_RADII.put(66, GK_AMOEBA_HSULFIDE_H / 2.0);
        AMOEBA_2014_GK_RADII.put(67, GK_AMOEBA_SULFIDE_S / 2.0);
        AMOEBA_2014_GK_RADII.put(68, GK_AMOEBA_SULFIDE_HS / 2.0);
        AMOEBA_2014_GK_RADII.put(69, GK_AMOEBA_SULFIDE_H3C / 2.0);
        AMOEBA_2014_GK_RADII.put(70, GK_AMOEBA_SULFIDE_CH3 / 2.0);
        AMOEBA_2014_GK_RADII.put(71, GK_AMOEBA_ETHYL_SULFIDE_H3C / 2.0);

        AMOEBA_2014_GK_RADII.put(79, GK_AMOEBA_ETHYL_NITRILE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(80, GK_AMOEBA_ETHYL_NITRILE_N / 2.0);
        AMOEBA_2014_GK_RADII.put(82, GK_AMOEBA_ETHYL_NITRILE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(83, GK_AMOEBA_ETHYL_NITRILE_H3C / 2.0);

        AMOEBA_2014_GK_RADII.put(87, GK_AMOEBA_BENZENE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(88, GK_AMOEBA_BENZENE_H / 2.0);
        AMOEBA_2014_GK_RADII.put(89, GK_AMOEBA_AROMATIC_C / 2.0);
        AMOEBA_2014_GK_RADII.put(90, GK_AMOEBA_AROMATIC_H / 2.0);
        AMOEBA_2014_GK_RADII.put(91, GK_AMOEBA_IMIDAZOLE_NH / 2.0);
        AMOEBA_2014_GK_RADII.put(92, GK_AMOEBA_IMIDAZOLE_C / 2.0);
        AMOEBA_2014_GK_RADII.put(93, GK_AMOEBA_IMIDAZOLE_HC / 2.0);

        AMOEBA_2014_GK_RADII.put(403, 5.3200 / 2.0); // 403    C     "Ethene"
        AMOEBA_2014_GK_RADII.put(404, 4.1720 / 2.0); // 404    H     "Ethene"
        AMOEBA_2014_GK_RADII.put(422, 3.4132 / 2.0); // 422   N     "Pyridine"
        AMOEBA_2014_GK_RADII.put(423, 3.4776 / 2.0); // 423   C     "Pyridine C-N"
        AMOEBA_2014_GK_RADII.put(424, 2.7600 / 2.0); // 424   H     "Pyridine HCN"
        AMOEBA_2014_GK_RADII.put(425, 3.4960 / 2.0); // 425   C     "Pyridine CCN"
        AMOEBA_2014_GK_RADII.put(426, 2.7416 / 2.0); // 426   H     "Pyridine HCCN"
        AMOEBA_2014_GK_RADII.put(427, 3.4960 / 2.0); // 427   C     "Pyridine C"
        AMOEBA_2014_GK_RADII.put(428, 2.7416 / 2.0); // 428   H     "Pyridine HC"
        AMOEBA_2014_GK_RADII.put(443, 4.5360 / 2.0); // 443    C     "CH3SH"
        AMOEBA_2014_GK_RADII.put(444, 4.6058 / 2.0); // 444    S     "CH3SH"
        AMOEBA_2014_GK_RADII.put(445, 3.4440 / 2.0); // 445    H     "CH3SH"
        AMOEBA_2014_GK_RADII.put(446, 3.3240 / 2.0); // 446    H     "CH3SH"
        AMOEBA_2014_GK_RADII.put(447, 4.3470 / 2.0); // 447    C     "DMSO"
        AMOEBA_2014_GK_RADII.put(448, 4.8060 / 2.0); // 448    S     "DMSO"
        AMOEBA_2014_GK_RADII.put(449, 3.5520 / 2.0); // 449    H     "DMSO"
        AMOEBA_2014_GK_RADII.put(450, 4.2120 / 2.0); // 450    O     "DMSO"
        AMOEBA_2014_GK_RADII.put(463, 4.5840 / 2.0); // 463    C     "MeF"
        AMOEBA_2014_GK_RADII.put(464, 3.3810 / 2.0); // 464    F     "MeF"
        AMOEBA_2014_GK_RADII.put(465, 3.5520 / 2.0); // 465    H     "MeF"
        AMOEBA_2014_GK_RADII.put(466, 4.2020 / 2.0); // 466    C     "MeCl"
        AMOEBA_2014_GK_RADII.put(467, 4.2895 / 2.0); // 467    Cl    "MeCl"
        AMOEBA_2014_GK_RADII.put(468, 3.2560 / 2.0); // 468    H     "MeCl"
        AMOEBA_2014_GK_RADII.put(469, 4.5840 / 2.0); // 469    C     "MeBr"
        AMOEBA_2014_GK_RADII.put(470, 4.3780 / 2.0); // 470    Br    "MeBr"
        AMOEBA_2014_GK_RADII.put(471, 3.5520 / 2.0); // 471    H     "MeBr"
        AMOEBA_2014_GK_RADII.put(472, 4.5600 / 2.0); // 472    C     "BenF"
        AMOEBA_2014_GK_RADII.put(473, 4.5600 / 2.0); // 473    C     "BenF"
        AMOEBA_2014_GK_RADII.put(474, 4.5600 / 2.0); // 474    C     "BenF"
        AMOEBA_2014_GK_RADII.put(475, 4.5600 / 2.0); // 475    C     "BenF"
        AMOEBA_2014_GK_RADII.put(476, 3.3810 / 2.0); // 476    F     "BenF"
        AMOEBA_2014_GK_RADII.put(477, 3.5760 / 2.0); // 477    H     "BenF"
        AMOEBA_2014_GK_RADII.put(478, 3.5760 / 2.0); // 478    H     "BenF"
        AMOEBA_2014_GK_RADII.put(479, 3.5760 / 2.0); // 479    H     "BenF"
        AMOEBA_2014_GK_RADII.put(480, 4.5600 / 2.0); // 480    C     "BenCl"
        AMOEBA_2014_GK_RADII.put(481, 4.5600 / 2.0); // 481    C     "BenCl"
        AMOEBA_2014_GK_RADII.put(482, 4.5600 / 2.0); // 482    C     "BenCl"
        AMOEBA_2014_GK_RADII.put(483, 4.5600 / 2.0); // 483    C     "BenCl"
        AMOEBA_2014_GK_RADII.put(484, 4.4760 / 2.0); // 484    Cl    "BenCl"
        AMOEBA_2014_GK_RADII.put(485, 3.5760 / 2.0); // 485    H     "BenCl"
        AMOEBA_2014_GK_RADII.put(486, 3.5760 / 2.0); // 486    H     "BenCl"
        AMOEBA_2014_GK_RADII.put(487, 3.5760 / 2.0); // 487    H     "BenCl"
        AMOEBA_2014_GK_RADII.put(488, 4.5600 / 2.0); // 488    C     "BenBr"
        AMOEBA_2014_GK_RADII.put(489, 4.5600 / 2.0); // 489    C     "BenBr"
        AMOEBA_2014_GK_RADII.put(490, 4.5600 / 2.0); // 490    C     "BenBr"
        AMOEBA_2014_GK_RADII.put(491, 4.5600 / 2.0); // 491    C     "BenBr"
        AMOEBA_2014_GK_RADII.put(492, 4.3780 / 2.0); // 492    Br    "BenBr"
        AMOEBA_2014_GK_RADII.put(493, 3.5760 / 2.0); // 493    H     "BenBr"
        AMOEBA_2014_GK_RADII.put(494, 3.5760 / 2.0); // 494    H     "BenBr"
        AMOEBA_2014_GK_RADII.put(495, 3.5760 / 2.0); // 495    H     "BenBr"
        AMOEBA_2014_GK_RADII.put(499, 4.5840 / 2.0); // 499    C     "MeCl2"
        AMOEBA_2014_GK_RADII.put(500, 4.6625 / 2.0); // 500    Cl    "MeCl2"
        AMOEBA_2014_GK_RADII.put(501, 3.5760 / 2.0); // 501    H     "MeCl2"
        AMOEBA_2014_GK_RADII.put(502, 4.5840 / 2.0); // 502    C     "MeBr2"
        AMOEBA_2014_GK_RADII.put(503, 4.7760 / 2.0); // 503    Br    "MeBr2"
        AMOEBA_2014_GK_RADII.put(504, 3.5760 / 2.0); // 504    H     "MeBr2"

        AMOEBA_NUC_2017_GK_RADII.put(1, GK_AMOEBA_AROMATIC_C / 2.0);       // 1   C   Adenine C4/C5, Guanine C4/C5
        AMOEBA_NUC_2017_GK_RADII.put(2, GK_AMOEBA_AROMATIC_C / 2.0);       // 2   C   Adenine C8, Guanine C8
        AMOEBA_NUC_2017_GK_RADII.put(3, GK_AMOEBA_AROMATIC_C / 2.0);       // 3   C   Adenine C6, Cytosine C4
        AMOEBA_NUC_2017_GK_RADII.put(4, GK_AMOEBA_IMIDAZOLE_NH / 2.0);     // 4   N   Adenine N9 RNA/DNA
        AMOEBA_NUC_2017_GK_RADII.put(5, GK_AMOEBA_METHYLAMINE_N / 2.0);    // 5   N   Adenine N6, Cytosine N4, Guanine N2
        AMOEBA_NUC_2017_GK_RADII.put(6, GK_AMOEBA_AROMATIC_C / 2.0);       // 6   C   Adenine C2, Cytosine C6, Thymine C5/C6, Uracil C5/C6
        AMOEBA_NUC_2017_GK_RADII.put(7, GK_AMOEBA_IMIDAZOLE_NH / 2.0);     // 7   N   Adenine N7, Guanine N7
        AMOEBA_NUC_2017_GK_RADII.put(8, GK_AMOEBA_IMIDAZOLE_NH / 2.0);     // 8   N   Adenine N1/N3, Cytosine N3, Guanine N3
        AMOEBA_NUC_2017_GK_RADII.put(9, GK_AMOEBA_AROMATIC_H / 2.0);       // 9   H   Adenine H8, Cytosine H6, Guanine H8, Thymine H6, Uracil H6
        AMOEBA_NUC_2017_GK_RADII.put(10, GK_AMOEBA_AROMATIC_H / 2.0);      // 10   H   Adenine H2, Uracil H5
        AMOEBA_NUC_2017_GK_RADII.put(11, GK_AMOEBA_METHYLAMINE_HN / 2.0);  // 11   H   Adenine H61, Cytosine H41, Guanine H21
        AMOEBA_NUC_2017_GK_RADII.put(12, GK_AMOEBA_IMIDAZOLE_NH / 2.0);    // 12   N   Cytosine N1 RNA/DNA, Guanine N1, Thymine N3, Uracil N3
        AMOEBA_NUC_2017_GK_RADII.put(13, GK_AMOEBA_BENZENE_C / 2.0);       // 13   C   Cytosine C5
        AMOEBA_NUC_2017_GK_RADII.put(14, GK_AMOEBA_AMIDE_CO / 2.0);        // 14   C   Cytosine C2, Guanine C6, Thymine C2/C4, Uracil C2/C4
        AMOEBA_NUC_2017_GK_RADII.put(15, GK_AMOEBA_BENZENE_H / 2.0);       // 15   H   Cytosine H5
        AMOEBA_NUC_2017_GK_RADII.put(16, GK_AMOEBA_CARBONYL_O / 2.0);      // 16   O   Cytosine O2, Guanine O6, Thymine O2/O4, Uracil O2/O4
        AMOEBA_NUC_2017_GK_RADII.put(17, GK_AMOEBA_IMIDAZOLE_C / 2.0);     // 17   C   Guanine C2 (this one is tricky...)
        AMOEBA_NUC_2017_GK_RADII.put(18, GK_AMOEBA_IMIDAZOLE_NH / 2.0);    // 18   N   Guanine N9 RNA/DNA
        AMOEBA_NUC_2017_GK_RADII.put(19, GK_AMOEBA_AMIDE_HN / 2.0);        // 19   H   Guanine H1, Thymine H3, Uracil H3
        AMOEBA_NUC_2017_GK_RADII.put(20, GK_AMOEBA_POLARGROUP_C / 2.0);       // 20   C   Thymine C7 (C5-attached methyl group C)
        AMOEBA_NUC_2017_GK_RADII.put(21, GK_AMOEBA_IMIDAZOLE_NH / 2.0);    // 21   N   Thymine N1, Uracil N1
        AMOEBA_NUC_2017_GK_RADII.put(22, GK_AMOEBA_POLARGROUP_HC / 2.0);       // 22   H   Thymine H7
        //AMOEBA_NUC_2017_GK_RADII.put(23, GK_AMOEBA_CARBONYL_O / 2.0);     // 23   O   Ribose O4 (Unsure on this one)
        //AMOEBA_NUC_2017_GK_RADII.put(23, GK_AMOEBA_IMIDAZOLE_C / 2.0);    // 24   C   Ribose C4

//        atom     65     23  O     "Ribose O4'"                 8   15.9990  2
//        atom     66     24  C     "Ribose C4'"                 6   12.0110  4
//        atom     67     25  H     "Ribose H4'"                 1    1.0080  1
//        atom     68     26  C     "Ribose C1'"                 6   12.0110  4
//        atom     69     25  H     "Ribose H1'"                 1    1.0080  1
//        atom     70     26  C     "Ribose C2'"                 6   12.0110  4
//        atom     71     25  H     "Ribose H2'1"                1    1.0080  1
//        atom     72     26  C     "Ribose C3'"                 6   12.0110  4
//        atom     73     25  H     "Ribose H3'"                 1    1.0080  1
//        atom     74     27  O     "Ribose O2'"                 8   15.9990  2
//        atom     75     28  H     "Ribose HO'2"                1    1.0080  1
//        atom     76     29  O     "Ribose O3'"                 8   15.9990  2
//        atom     77     28  H     "Ribose HO'3"                1    1.0080  1
//        atom     78     30  C     "Ribose C5'"                 6   12.0110  4
//        atom     79     31  H     "Ribose H5'"                 1    1.0080  1
//        atom     80     32  O     "Ribose O5'"                 8   15.9990  2
//        atom     81     28  H     "Ribose HO'5"                1    1.0080  1

        AMOEBA_NUC_2017_GK_RADII.put(48, GK_AMOEBA_WATER_O / 2.0);
        AMOEBA_NUC_2017_GK_RADII.put(49, GK_AMOEBA_WATER_H / 2.0);
        AMOEBA_NUC_2017_GK_RADII.put(50, GK_AMOEBA_SODIUM / 2.0);
        AMOEBA_NUC_2017_GK_RADII.put(51, GK_AMOEBA_POTASSIUM / 2.0);
        AMOEBA_NUC_2017_GK_RADII.put(52, GK_AMOEBA_MAGNESIUM / 2.0);
        AMOEBA_NUC_2017_GK_RADII.put(53, GK_AMOEBA_CHLORIDE / 2.0);

        // AMOEBA BIO 2018 Protein Force Field
        AMOEBA_BIO_2018_GK_RADII.put(1, GK_AMOEBA_AMIDE_N / 2.0);       // Glycine N, Alanine N
        AMOEBA_BIO_2018_GK_RADII.put(2, GK_AMOEBA_POLARGROUP_C / 2.0);  // Glycine CA
        AMOEBA_BIO_2018_GK_RADII.put(3, GK_AMOEBA_AMIDE_CO / 2.0);      // Glycine C, Alanine C
        AMOEBA_BIO_2018_GK_RADII.put(4, GK_AMOEBA_AMIDE_HN / 2.0);      // Glycine HN, Alanine HN
        AMOEBA_BIO_2018_GK_RADII.put(5, GK_AMOEBA_CARBONYL_O / 2.0);    // Glycine O, Alanine O
        AMOEBA_BIO_2018_GK_RADII.put(6, GK_AMOEBA_POLARGROUP_HC / 2.0); // Glycine HA, Alanine HA
        AMOEBA_BIO_2018_GK_RADII.put(7, GK_AMOEBA_ISOPROPANOL_C / 2.0);  // Alanine CA, Cysteine Anion CA
        AMOEBA_BIO_2018_GK_RADII.put(8, GK_AMOEBA_ALKANE_C / 2.0); // Alanine CB
        AMOEBA_BIO_2018_GK_RADII.put(9, GK_AMOEBA_ALKANE_H / 2.0);    // Alanine HB
        AMOEBA_BIO_2018_GK_RADII.put(10, GK_AMOEBA_ALCOHOL_O / 2.0);    // Serine OG
        AMOEBA_BIO_2018_GK_RADII.put(11, GK_AMOEBA_ALCOHOL_HO / 2.0);   // Serine HG
        AMOEBA_BIO_2018_GK_RADII.put(12, GK_AMOEBA_HSULFIDE_S / 2.0);    // Cysteine SG, Cystine SG
        AMOEBA_BIO_2018_GK_RADII.put(13, GK_AMOEBA_HSULFIDE_H / 2.0);   // Cysteine HG
        AMOEBA_BIO_2018_GK_RADII.put(14, GK_AMOEBA_HSULFIDE_S / 2.0);    // TODO: Cysteine Anion S-
        AMOEBA_BIO_2018_GK_RADII.put(15, GK_AMOEBA_AMIDE_N / 2.0);      // Proline N
        AMOEBA_BIO_2018_GK_RADII.put(16, GK_AMOEBA_POLARGROUP_C / 2.0); // Proline CB, CG, CD
        AMOEBA_BIO_2018_GK_RADII.put(17, GK_AMOEBA_AROMATIC_C / 2.0);   // Phenylalanine CG, CD, CE, CZ
        AMOEBA_BIO_2018_GK_RADII.put(18, GK_AMOEBA_AROMATIC_H / 2.0);   // Phenylalanine HD, HE, HZ
        AMOEBA_BIO_2018_GK_RADII.put(19, GK_AMOEBA_ALCOHOL_O / 2.0);    // Tyrosine OH
        AMOEBA_BIO_2018_GK_RADII.put(20, GK_AMOEBA_ALCOHOL_HO / 2.0);   // Tyrosine HH
        AMOEBA_BIO_2018_GK_RADII.put(21, GK_AMOEBA_ALCOHOL_O / 2.0);    // TODO: Tyrosine O-
        AMOEBA_BIO_2018_GK_RADII.put(22, GK_AMOEBA_AROMATIC_C / 2.0);   // Tryptophan CG, CD1, CD2, CE2, CE3, CZ2, CZ3, CH2
        AMOEBA_BIO_2018_GK_RADII.put(23, GK_AMOEBA_AROMATIC_H / 2.0);   // Tryptophan HD1, HE3, HZ2, HZ3, HH2
        AMOEBA_BIO_2018_GK_RADII.put(24, GK_AMOEBA_PROTEIN_HISTIDINE_NH / 2.0); // Tryptophan NE1, Histidine (+) ND1, NE2
        AMOEBA_BIO_2018_GK_RADII.put(25, GK_AMOEBA_NITROGEN_H / 2.0);   // Tryptophan HE1
        AMOEBA_BIO_2018_GK_RADII.put(26, GK_AMOEBA_IMIDAZOLE_C / 2.0);  // Histidine CG, CD2
        AMOEBA_BIO_2018_GK_RADII.put(27, GK_AMOEBA_IMIDAZOLE_HC / 2.0); // Histidine HD2, HE1
        AMOEBA_BIO_2018_GK_RADII.put(28, GK_AMOEBA_IMIDAZOLE_C / 2.0);  // Histidine CE1
        AMOEBA_BIO_2018_GK_RADII.put(29, GK_AMOEBA_IMIDAZOLE_NH / 2.0); // Histidine (HE) ND1, (HD) NE2
        AMOEBA_BIO_2018_GK_RADII.put(30, GK_AMOEBA_PROTEIN_CARBOXYLIC_ACID_CO / 2.0);  // Aspartate CG. Glutamate CD
        AMOEBA_BIO_2018_GK_RADII.put(31, GK_AMOEBA_PROTEIN_CARBOXYLIC_ACID_O / 2.0);   // Aspartate OD, Glutamate OE
        AMOEBA_BIO_2018_GK_RADII.put(32, GK_AMOEBA_CARBONYL_CO / 2.0);  // Aspartic Acid CG, Glutamic Acid CD
        AMOEBA_BIO_2018_GK_RADII.put(33, GK_AMOEBA_CARBONYL_O / 2.0);   // Aspartic Acid OD1, Glutamic Acid OE1
        AMOEBA_BIO_2018_GK_RADII.put(34, GK_AMOEBA_CARBOXCYLIC_ACID_O / 2.0);   // Aspartic Acid OD2, Glutamic Acid OE2
        AMOEBA_BIO_2018_GK_RADII.put(35, GK_AMOEBA_CARBOXCYLIC_ACID_HO / 2.0);  // Aspartic Acid HD2, Glutamic Acid HE2
        AMOEBA_BIO_2018_GK_RADII.put(36, GK_AMOEBA_ALKANE_C / 2.0); // Lysine CG
        AMOEBA_BIO_2018_GK_RADII.put(37, GK_AMOEBA_PROTEIN_LYSINE_NH / 2.0); // Lysine NZ
        AMOEBA_BIO_2018_GK_RADII.put(38, GK_AMOEBA_PROTEIN_LYSINE_HN / 2.0); // Lysine HN
        AMOEBA_BIO_2018_GK_RADII.put(39, GK_AMOEBA_PROTEIN_ARGININE_CZ / 2.0); // Arginine CZ
        AMOEBA_BIO_2018_GK_RADII.put(40, GK_AMOEBA_POLARGROUP_C / 2.0); // Acetyl Cap CH3
        AMOEBA_BIO_2018_GK_RADII.put(41, GK_AMOEBA_PROTEIN_TERMINAL_NH3 / 2.0); // N-Terminal NH3+
        AMOEBA_BIO_2018_GK_RADII.put(42, GK_AMOEBA_PROTEIN_TERMINAL_H3N / 2.0); // N-Terminal H3N+

        // AMOEBA BIO 2018 Nucleic Acid Force Field
        AMOEBA_BIO_2018_GK_RADII.put(43, AMOEBA_NUC_2017_GK_RADII.get(1));
        AMOEBA_BIO_2018_GK_RADII.put(44, AMOEBA_NUC_2017_GK_RADII.get(2));
        AMOEBA_BIO_2018_GK_RADII.put(45, AMOEBA_NUC_2017_GK_RADII.get(3));
        AMOEBA_BIO_2018_GK_RADII.put(46, AMOEBA_NUC_2017_GK_RADII.get(4));
        AMOEBA_BIO_2018_GK_RADII.put(47, AMOEBA_NUC_2017_GK_RADII.get(5));
        AMOEBA_BIO_2018_GK_RADII.put(48, AMOEBA_NUC_2017_GK_RADII.get(6));
        AMOEBA_BIO_2018_GK_RADII.put(49, AMOEBA_NUC_2017_GK_RADII.get(7));
        AMOEBA_BIO_2018_GK_RADII.put(50, AMOEBA_NUC_2017_GK_RADII.get(8));
        AMOEBA_BIO_2018_GK_RADII.put(51, AMOEBA_NUC_2017_GK_RADII.get(9));
        AMOEBA_BIO_2018_GK_RADII.put(52, AMOEBA_NUC_2017_GK_RADII.get(10));
        AMOEBA_BIO_2018_GK_RADII.put(53, AMOEBA_NUC_2017_GK_RADII.get(11));
        AMOEBA_BIO_2018_GK_RADII.put(54, AMOEBA_NUC_2017_GK_RADII.get(12));
        AMOEBA_BIO_2018_GK_RADII.put(55, AMOEBA_NUC_2017_GK_RADII.get(13));
        AMOEBA_BIO_2018_GK_RADII.put(56, AMOEBA_NUC_2017_GK_RADII.get(14));
        AMOEBA_BIO_2018_GK_RADII.put(57, AMOEBA_NUC_2017_GK_RADII.get(15));
        AMOEBA_BIO_2018_GK_RADII.put(58, AMOEBA_NUC_2017_GK_RADII.get(16));
        AMOEBA_BIO_2018_GK_RADII.put(59, AMOEBA_NUC_2017_GK_RADII.get(17));
        AMOEBA_BIO_2018_GK_RADII.put(60, AMOEBA_NUC_2017_GK_RADII.get(18));
        AMOEBA_BIO_2018_GK_RADII.put(61, AMOEBA_NUC_2017_GK_RADII.get(19));
        AMOEBA_BIO_2018_GK_RADII.put(62, AMOEBA_NUC_2017_GK_RADII.get(20));
        AMOEBA_BIO_2018_GK_RADII.put(63, AMOEBA_NUC_2017_GK_RADII.get(21));
        AMOEBA_BIO_2018_GK_RADII.put(64, AMOEBA_NUC_2017_GK_RADII.get(22));

        AMOEBA_BIO_2018_GK_RADII.put(92, GK_AMOEBA_LITHIUM / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(93, GK_AMOEBA_SODIUM / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(94, GK_AMOEBA_POTASSIUM / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(95, GK_AMOEBA_RUBIDIUM / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(96, GK_AMOEBA_CESIUM / 2.0);
        // Class 97: Berylium Ion Be+2
        AMOEBA_BIO_2018_GK_RADII.put(98, GK_AMOEBA_MAGNESIUM / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(99, GK_AMOEBA_CALCIUM / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(100, GK_AMOEBA_ZINC / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(101, GK_AMOEBA_FLUORIDE / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(102, GK_AMOEBA_CHLORIDE / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(103, GK_AMOEBA_BROMIDE / 2.0);
        AMOEBA_BIO_2018_GK_RADII.put(104, GK_AMOEBA_IODIDE / 2.0);


    }
}
