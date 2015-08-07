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
package ffx.algorithms;

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Torsion;
import static ffx.potential.bonded.BondedUtils.intxyz;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import org.apache.commons.math3.util.FastMath;

/**
 * Checks whether asparagine and glutamine residues are in the correct orientation.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 */
public class NQFlipper {
    private static final Logger logger = Logger.getLogger(NQFlipper.class.getName());
    
    private final MolecularAssembly assembly;
    private List<Residue> residueList;
    private String chainID;
    private double[] x;
    private final Potential potential;
    /**
     * useIdealGeometry over-rides useOwnGeometry, and chooses between original 
     * coordinates exactly, original torsions with ideal bond/angle, and flipped
     * torsions with ideal bond/angle.
     */
    private boolean useIdealGeometry = true;
    /**
     * If not using ideal geometry, this flag essentially rotates the NQ, switching
     * only the torsion, while leaving the old geometry in place; if false, the 
     * atoms exactly switch coordinates.
     */
    private boolean useOwnGeometry = true;
    
    public NQFlipper(MolecularAssembly assembly, Potential potential) {
        this.assembly = assembly;
        this.potential = potential;
    }
    
    public void setResidues(List<Residue> residueList) {
        this.residueList = new ArrayList<>(residueList);
    }
    
    /**
     * Set a contiguous block of residues to optimize in a specific chain.
     *
     * @param chain
     * @param startResID
     * @param finalResID
     */
    public void setResidues(String chain, int startResID, int finalResID) {
        this.chainID = chain;
        this.setResidues(startResID, finalResID);
    }

    /**
     * Set a contiguous block of residues to optimize.
     *
     * @param startResID
     * @param finalResID
     */
    public void setResidues(int startResID, int finalResID) {
        Polymer polymer;
        Polymer[] polymers;
        
        if (chainID != null) {
            polymer = assembly.getChain(chainID);
        } else {
            polymers = assembly.getPolymers();
            polymer = polymers[0];
        }
        
        residueList = new ArrayList<>();
        for (int i = startResID; i <= finalResID; i++) {
            Residue residue = polymer.getResidue(i);
            if (residue != null) {
                Rotamer[] rotamers = residue.getRotamers(residue);
                if (rotamers != null) {
                    if (rotamers.length == 1) {
                        switch (residue.getResidueType()) {
                            case NA:
                                residue.initializeDefaultAtomicCoordinates();
                                break;
                            case AA:
                            default:
                                RotamerLibrary.applyRotamer(residue, rotamers[0]);
                                break;
                        }
                    } else {
                        residueList.add(residue);
                    }
                }
            }
        }
    }
    
    public void flipNQs() {
        String useOwnGeoStr = System.getProperty("nq-useOwnGeometry");
        if (useOwnGeoStr != null) {
            useOwnGeometry = Boolean.parseBoolean(useOwnGeoStr);
        }
        String useIdealGeoStr = System.getProperty("nq-useIdealGeometry");
        if (useIdealGeoStr != null) {
            useIdealGeometry = Boolean.parseBoolean(useIdealGeoStr);
        }
        
        logger.info(" Attempting to flip asparagines and glutamines...");
        if (useIdealGeometry) {
            flipNQsIdeal();
        } else {
            flipNQsOrig();
        }
    }
    
    private void flipNQsIdeal() {
        List<Residue> originalResidues = new ArrayList<>();
        List<Residue> origGeoResidues = new ArrayList<>();
        List<Residue> flippedResidues = new ArrayList<>();
        for (Residue residue : residueList) {
            String resName = residue.getPDBName().toUpperCase();
            if (resName.equals("ASN") || resName.equals("GLN")) {
                
                double initialEnergy = currentEnergy();
                logger.info(String.format(" Initial energy for %7s: %16.8f", residue, initialEnergy));
                ResidueState orig = residue.storeCoordinates();
                
                try {
                    applyIdealGeometry(residue, resName, false);
                } catch (IllegalArgumentException ex) {
                    logger.warning(String.format(" Could not flip residue %s due "
                            + "to null atom.", residue.toString()));
                    continue;
                }
                double origGeoEnergy = currentEnergy();
                logger.info(String.format(" Ideal-geometry energy for %7s: %16.8f", residue, origGeoEnergy));
                ResidueState origGeo = residue.storeCoordinates();
                
                applyIdealGeometry(residue, resName, true);
                double flippedEnergy = currentEnergy();
                logger.info(String.format(" Flipped energy for %7s: %16.8f", residue, flippedEnergy));
                ResidueState flipped = residue.storeCoordinates();
                
                if (flippedEnergy < initialEnergy && flippedEnergy < origGeoEnergy) {
                    logger.info(String.format(" Flipping %s with flipped energy %16.8f "
                            + "< %16.8f", residue.toString(), flippedEnergy, initialEnergy));
                    applyIdealGeometry(residue, resName, true);
                    residue.revertCoordinates(flipped);
                } else if (origGeoEnergy < initialEnergy) {
                    logger.info(String.format(" Applying ideal geometry to %s with "
                            + "energy %16.8f < %16.8f", residue, origGeoEnergy, initialEnergy));
                    applyIdealGeometry(residue, resName, false);
                    residue.revertCoordinates(origGeo);
                } else {
                    logger.info(String.format(" Retaining conformation for %s "
                            + "with energy %16.8f", residue.toString(), initialEnergy));
                    originalResidues.add(residue);
                    residue.revertCoordinates(orig);
                }
            }
        }
        logger.info(" Residues retaining original conformation: ");
        for (Residue residue : originalResidues) {
            logger.info(residue.toString());
        }
        logger.info(" Residues using idealized geometry: ");
        for (Residue residue : origGeoResidues) {
            logger.info(residue.toString());
        }
        logger.info(" Residues with flipped conformation: ");
        for (Residue residue : flippedResidues) {
            logger.info(residue.toString());
        }
    }
    
    private void applyIdealGeometry(Residue residue, String resname, boolean flip) {
        if (resname.equalsIgnoreCase("ASN")) {
            applyIdealAsn(residue, flip);
        } else if (resname.equalsIgnoreCase("GLN")) {
            applyIdealGln(residue, flip);
        } else {
            throw new IllegalArgumentException(String.format(" Non-NQ residue %s", residue.toString()));
        }
    }
    
    private void applyIdealAsn(Residue residue, boolean flip) {
        Atom OD1 = (Atom) residue.getAtomNode("OD1");
        Atom ND2 = (Atom) residue.getAtomNode("ND2");
        Atom CA = (Atom) residue.getAtomNode("CA");
        Atom CB = (Atom) residue.getAtomNode("CB");
        Atom CG = (Atom) residue.getAtomNode("CG");
        Atom HD21 = (Atom) residue.getAtomNode("HD21");
        Atom HD22 = (Atom) residue.getAtomNode("HD22");
        if (OD1 == null || ND2 == null || CB == null || CA == null || CG == null 
                || HD21 == null || HD22 == null) {
            throw new IllegalArgumentException(String.format("Null atom "
                    + "in residue %s", residue.toString()));
        }

        Bond OD1_CG = OD1.getBond(CG);
        Angle OD1_CG_CB = OD1.getAngle(CG, CB);
        Torsion OD1_CG_CB_CA = OD1.getTorsion(CG, CB, CA);
        double dOD1_CG = OD1_CG.bondType.distance;
        double dOD1_CG_CB = OD1_CG_CB.angleType.angle[OD1_CG_CB.nh];
        double dOD1_CG_CB_CA = OD1_CG_CB_CA.getValue();

        Bond ND2_CG = ND2.getBond(CG);
        Angle ND2_CG_CB = ND2.getAngle(CG, CB);
        Torsion ND2_CG_CB_CA = ND2.getTorsion(CG, CB, CA);
        double dND2_CG = ND2_CG.bondType.distance;
        double dND2_CG_CB = ND2_CG_CB.angleType.angle[ND2_CG_CB.nh];
        double dND2_CG_CB_CA = ND2_CG_CB_CA.getValue();

        Bond HD21_ND2 = HD21.getBond(ND2);
        Angle HD21_ND2_CG = HD21.getAngle(ND2, CG);
        Torsion HD21_ND2_CG_CB = HD21.getTorsion(ND2, CG, CB);
        double dHD21_ND2 = HD21_ND2.bondType.distance;
        double dHD21_ND2_CG = HD21_ND2_CG.angleType.angle[HD21_ND2_CG.nh];
        double dHD21_ND2_CG_CB = HD21_ND2_CG_CB.getValue();

        Bond HD22_ND2 = HD22.getBond(ND2);
        Angle HD22_ND2_CG = HD22.getAngle(ND2, CG);
        Torsion HD22_ND2_CG_CB = HD22.getTorsion(ND2, CG, CB);
        double dHD22_ND2 = HD22_ND2.bondType.distance;
        double dHD22_ND2_CG = HD22_ND2_CG.angleType.angle[HD22_ND2_CG.nh];
        double dHD22_ND2_CG_CB = HD22_ND2_CG_CB.getValue();

        if (flip) {
            intxyz(OD1, CG, dND2_CG, CB, dND2_CG_CB, CA, dND2_CG_CB_CA, 0);
            intxyz(ND2, CG, dOD1_CG, CB, dOD1_CG_CB, CA, dOD1_CG_CB_CA, 0);
        } else {
            intxyz(OD1, CG, dND2_CG, CB, dND2_CG_CB, CA, dOD1_CG_CB_CA, 0);
            intxyz(ND2, CG, dOD1_CG, CB, dOD1_CG_CB, CA, dND2_CG_CB_CA, 0);
        }
        intxyz(HD21, ND2, dHD21_ND2, CG, dHD21_ND2_CG, CB, dHD21_ND2_CG_CB, 0);
        intxyz(HD22, ND2, dHD22_ND2, CG, dHD22_ND2_CG, CB, dHD22_ND2_CG_CB, 0);
    }
    
    private void applyIdealGln(Residue residue, boolean flip) {
        Atom OE1 = (Atom) residue.getAtomNode("OE1");
        Atom NE2 = (Atom) residue.getAtomNode("NE2");
        Atom CB = (Atom) residue.getAtomNode("CB");
        Atom CG = (Atom) residue.getAtomNode("CG");
        Atom CD = (Atom) residue.getAtomNode("CD");
        Atom HE21 = (Atom) residue.getAtomNode("HE21");
        Atom HE22 = (Atom) residue.getAtomNode("HE22");
        if (OE1 == null || NE2 == null || CG == null || CB == null || CD == null 
                || HE21 == null || HE22 == null) {
            throw new IllegalArgumentException(String.format("Null atom "
                    + "in residue %s", residue.toString()));
        }

        Bond OE1_CD = OE1.getBond(CD);
        Angle OE1_CD_CG = OE1.getAngle(CD, CG);
        Torsion OE1_CD_CG_CB = OE1.getTorsion(CD, CG, CB);
        double dOE1_CD = OE1_CD.bondType.distance;
        double dOE1_CD_CG = OE1_CD_CG.angleType.angle[OE1_CD_CG.nh];
        double dOE1_CD_CG_CB = OE1_CD_CG_CB.getValue();

        Bond NE2_CD = NE2.getBond(CD);
        Angle NE2_CD_CG = NE2.getAngle(CD, CG);
        Torsion NE2_CD_CG_CB = NE2.getTorsion(CD, CG, CB);
        double dNE2_CD = NE2_CD.bondType.distance;
        double dNE2_CD_CG = NE2_CD_CG.angleType.angle[NE2_CD_CG.nh];
        double dNE2_CD_CG_CB = NE2_CD_CG_CB.getValue();

        Bond HE21_NE2 = HE21.getBond(NE2);
        Angle HE21_NE2_CD = HE21.getAngle(NE2, CD);
        Torsion HE21_NE2_CD_CG = HE21.getTorsion(NE2, CD, CG);
        double dHE21_NE2 = HE21_NE2.bondType.distance;
        double dHE21_NE2_CD = HE21_NE2_CD.angleType.angle[HE21_NE2_CD.nh];
        double dHE21_NE2_CD_CG = HE21_NE2_CD_CG.getValue();

        Bond HD22_NE2 = HE22.getBond(NE2);
        Angle HD22_NE2_CD = HE22.getAngle(NE2, CD);
        Torsion HD22_NE2_CD_CG = HE22.getTorsion(NE2, CD, CG);
        double dHD22_NE2 = HD22_NE2.bondType.distance;
        double dHD22_NE2_CD = HD22_NE2_CD.angleType.angle[HD22_NE2_CD.nh];
        double dHD22_NE2_CD_CG = HD22_NE2_CD_CG.getValue();

        if (flip) {
            intxyz(OE1, CD, dNE2_CD, CG, dNE2_CD_CG, CB, dNE2_CD_CG_CB, 0);
            intxyz(NE2, CD, dOE1_CD, CG, dOE1_CD_CG, CB, dOE1_CD_CG_CB, 0);
        } else {
            intxyz(OE1, CD, dNE2_CD, CG, dNE2_CD_CG, CB, dOE1_CD_CG_CB, 0);
            intxyz(NE2, CD, dOE1_CD, CG, dOE1_CD_CG, CB, dNE2_CD_CG_CB, 0);
        }
        intxyz(HE21, NE2, dHE21_NE2, CD, dHE21_NE2_CD, CG, dHE21_NE2_CD_CG, 0);
        intxyz(HE22, NE2, dHD22_NE2, CD, dHD22_NE2_CD, CG, dHD22_NE2_CD_CG, 0);
    }
    
    private void flipNQsOrig() {
        List<Residue> originalResidues = new ArrayList<>();
        List<Residue> flippedResidues = new ArrayList<>();
        for (Residue residue : residueList) {
            String resName = residue.getPDBName().toUpperCase();
            if (resName.equals("ASN") || resName.equals("GLN")) {
                double initialEnergy = currentEnergy();
                logger.info(String.format(" Initial energy for %7s: %16.8f", residue, initialEnergy));
                
                ResidueState orig = residue.storeCoordinates();
                try {
                    flipRes(residue);
                } catch (IllegalArgumentException ex) {
                    logger.warning(String.format(" Could not flip residue %s due "
                            + "to null atom.", residue.toString()));
                    continue;
                }
                
                double flippedEnergy = currentEnergy();
                if (flippedEnergy < initialEnergy) {
                    logger.info(String.format(" Flipping %s with flipped energy %16.8f "
                            + "< %16.8f", residue.toString(), flippedEnergy, initialEnergy));
                    flippedResidues.add(residue);
                } else {
                    logger.info(String.format(" Retaining conformation for %s "
                            + "with flipped energy %16.8f > %16.8f", residue.toString(), 
                            flippedEnergy, initialEnergy));
                    residue.revertCoordinates(orig);
                    originalResidues.add(residue);
                }
            }
        }
        logger.info(" Residues retaining original conformation: ");
        for (Residue residue : originalResidues) {
            logger.info(residue.toString());
        }
        logger.info(" Residues with flipped conformation: ");
        for (Residue residue : flippedResidues) {
            logger.info(residue.toString());
        }
    }
    
    private void flipRes(Residue residue) {
        if (residue.getName().equalsIgnoreCase("ASN")) {
            flipAsn(residue);
        } else if (residue.getName().equalsIgnoreCase("GLN")) {
            flipGln(residue);
        } else {
            logger.warning(String.format(" Cannot flip non-ASN/GLN residue %s", residue.toString()));
        }
    }
    
    private void flipAsn(Residue residue) {
        Atom OD1 = (Atom) residue.getAtomNode("OD1");
        Atom ND2 = (Atom) residue.getAtomNode("ND2");
        Atom CA = (Atom) residue.getAtomNode("CA");
        Atom CB = (Atom) residue.getAtomNode("CB");
        Atom CG = (Atom) residue.getAtomNode("CG");
        Atom HD21 = (Atom) residue.getAtomNode("HD21");
        Atom HD22 = (Atom) residue.getAtomNode("HD22");
        if (OD1 == null || ND2 == null || CB == null || CA == null || CG == null 
                || HD21 == null || HD22 == null) {
            throw new IllegalArgumentException(String.format("Null atom "
                    + "in residue %s", residue.toString()));
        }

        Bond OD1_CG = OD1.getBond(CG);
        Torsion OD1_CG_CB_CA = OD1.getTorsion(CG, CB, CA);
        double dOD1_CG = OD1_CG.getLength();
        double dOD1_CG_CB = calculateAngle(OD1, CG, CB);
        double dOD1_CG_CB_CA = OD1_CG_CB_CA.getValue();

        Bond ND2_CG = ND2.getBond(CG);
        Torsion ND2_CG_CB_CA = ND2.getTorsion(CG, CB, CA);
        double dND2_CG = ND2_CG.getLength();
        double dND2_CG_CB = calculateAngle(ND2, CG, CB);
        double dND2_CG_CB_CA = ND2_CG_CB_CA.getValue();

        Bond HD21_ND2 = HD21.getBond(ND2);
        Torsion HD21_ND2_CG_CB = HD21.getTorsion(ND2, CG, CB);
        double dHD21_ND2 = HD21_ND2.getLength();
        double dHD21_ND2_CG = calculateAngle(HD21, ND2, CG);
        double dHD21_ND2_CG_CB = HD21_ND2_CG_CB.getValue();

        Bond HD22_ND2 = HD22.getBond(ND2);
        Torsion HD22_ND2_CG_CB = HD22.getTorsion(ND2, CG, CB);
        double dHD22_ND2 = HD22_ND2.getLength();
        double dHD22_ND2_CG = calculateAngle(HD22, ND2, CG);
        double dHD22_ND2_CG_CB = HD22_ND2_CG_CB.getValue();

        if (useOwnGeometry) {
            intxyz(OD1, CG, dOD1_CG, CB, dOD1_CG_CB, CA, dND2_CG_CB_CA, 0);
            intxyz(ND2, CG, dND2_CG, CB, dND2_CG_CB, CA, dOD1_CG_CB_CA, 0);
        } else {
            intxyz(OD1, CG, dND2_CG, CB, dND2_CG_CB, CA, dND2_CG_CB_CA, 0);
            intxyz(ND2, CG, dOD1_CG, CB, dOD1_CG_CB, CA, dOD1_CG_CB_CA, 0);
        }
        intxyz(HD21, ND2, dHD21_ND2, CG, dHD21_ND2_CG, CB, dHD21_ND2_CG_CB, 1);
        intxyz(HD22, ND2, dHD22_ND2, CG, dHD22_ND2_CG, CB, dHD22_ND2_CG_CB, -1);
    }
    
    private void flipGln(Residue residue) {
        Atom OE1 = (Atom) residue.getAtomNode("OE1");
        Atom NE2 = (Atom) residue.getAtomNode("NE2");
        Atom CB = (Atom) residue.getAtomNode("CB");
        Atom CG = (Atom) residue.getAtomNode("CG");
        Atom CD = (Atom) residue.getAtomNode("CD");
        Atom HE21 = (Atom) residue.getAtomNode("HE21");
        Atom HE22 = (Atom) residue.getAtomNode("HE22");
        if (OE1 == null || NE2 == null || CG == null || CB == null || CD == null 
                || HE21 == null || HE22 == null) {
            throw new IllegalArgumentException(String.format("Null atom "
                    + "in residue %s", residue.toString()));
        }

        Bond OE1_CD = OE1.getBond(CD);
        Torsion OE1_CD_CG_CB = OE1.getTorsion(CD, CG, CB);
        double dOE1_CD = OE1_CD.getLength();
        double dOE1_CD_CG = calculateAngle(OE1, CD, CG);
        double dOE1_CD_CG_CB = OE1_CD_CG_CB.getValue();
        
        Bond NE2_CD = NE2.getBond(CD);
        Torsion NE2_CD_CG_CB = NE2.getTorsion(CD, CG, CB);
        double dNE2_CD = NE2_CD.getLength();
        double dNE2_CD_CG = calculateAngle(NE2, CD, CG);
        double dNE2_CD_CG_CB = NE2_CD_CG_CB.getValue();

        // Let's avoid little things like "not actually flipping" here...
        double oxygenTors = dOE1_CD_CG_CB;
        double nitrogenTors = dNE2_CD_CG_CB;

        Bond HE21_NE2 = HE21.getBond(NE2);
        Torsion HE21_NE2_CD_CG = HE21.getTorsion(NE2, CD, CG);
        double dHE21_NE2 = HE21_NE2.getLength();
        double dHE21_NE2_CD = calculateAngle(HE21, NE2, CD);
        double dHE21_NE2_CD_CG = HE21_NE2_CD_CG.getValue();

        Bond HE22_NE2 = HE22.getBond(NE2);
        Torsion HE22_NE2_CD_CG = HE22.getTorsion(NE2, CD, CG);
        double dHE22_NE2 = HE22_NE2.getLength();
        double dHE22_NE2_CD = calculateAngle(HE22, NE2, CD);
        double dHE22_NE2_CD_CG = HE22_NE2_CD_CG.getValue();

        if (useOwnGeometry) {
            intxyz(OE1, CD, dOE1_CD, CG, dOE1_CD_CG, CB, dNE2_CD_CG_CB, 0);
            intxyz(NE2, CD, dNE2_CD, CG, dNE2_CD_CG, CB, dOE1_CD_CG_CB, 0);
        } else {
            intxyz(OE1, CD, dNE2_CD, CG, dNE2_CD_CG, CB, dNE2_CD_CG_CB, 0);
            intxyz(NE2, CD, dOE1_CD, CG, dOE1_CD_CG, CB, dOE1_CD_CG_CB, 0);
        }
        intxyz(HE21, NE2, dHE21_NE2, CD, dHE21_NE2_CD, CG, dHE21_NE2_CD_CG, 0);
        intxyz(HE22, NE2, dHE22_NE2, CD, dHE22_NE2_CD, CG, dHE22_NE2_CD_CG, 0);
    }
    
    /**
     * Because Angle.getValue() does not appear to return a useable angle, this
     * method calculates the angle between three Atoms.
     * @param a1 An Atom.
     * @param a2 An Atom.
     * @param a3 An Atom.
     * @return Angle between a1, a2, and a3.
     */
    private double calculateAngle(org.biojava.nbio.structure.Atom a1, 
            org.biojava.nbio.structure.Atom a2, org.biojava.nbio.structure.Atom a3) {
        double[] x1 = a1.getCoords();
        double[] x2 = a2.getCoords();
        double[] x3 = a3.getCoords();
        double result = 0; // Starts as dot product.
        double mag1 = 0;
        double mag2 = 0;
        for (int i = 0; i < 3; i++) {
            double d12 = x1[i] - x2[i];
            double d23 = x2[i] - x3[i];
            result += (d12 * d23);
            mag1 += (d12 * d12);
            mag2 += (d23 * d23);
        }
        mag1 = FastMath.sqrt(mag1);
        mag2 = FastMath.sqrt(mag2);
        result = FastMath.toDegrees(FastMath.acos(result / (mag1 * mag2)));
        
        /* Because of how I wrote it, it returns the angle between (12) and (23),
         * when the angle of interest is between (12) and (32). Thus, subtract from 180.
        */
        return 180.0 - result;
    }
    
    /**
     * Calculates the energy at the current state. Copied from RotamerOptimizaition;
     * perhaps possible to make that a public method?
     *
     * @return Energy of the current state.
     */
    private double currentEnergy() {
        if (x == null) {
            int nVar = potential.getNumberOfVariables();
            x = new double[nVar];
        }
        try {
            potential.getCoordinates(x);
            return potential.energy(x);
        } catch (ArithmeticException ex) {
            logger.warning(ex.getMessage());
            return Double.MAX_VALUE;
        }
    }
}
