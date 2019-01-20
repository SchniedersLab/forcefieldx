/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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

import javax.swing.tree.MutableTreeNode;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.extended.TitrationUtils.TitrationType;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.TorsionType;
import static ffx.potential.bonded.BondedUtils.buildBond;
import static ffx.potential.bonded.BondedUtils.intxyz;

/**
 * The MultiResidue class allows switching between residues for uses such as
 * sequence optimization.
 *
 * @author Will Tollefson
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class MultiTerminus extends Residue {

    private static final Logger logger = Logger.getLogger(MultiResidue.class.getName());

    private final ForceField forceField;
    private final ForceFieldEnergy forceFieldEnergy;
    private final MolecularAssembly mola;
    public final END end;
    public boolean isCharged;
    /**
     * Constant <code>kB=0.83144725</code>
     */
    public static final double kB = 0.83144725;
    private final static boolean DEBUG = false;

    /**
     * <p>Constructor for MultiTerminus.</p>
     *
     * @param residue          a {@link ffx.potential.bonded.Residue} object.
     * @param forceField       a {@link ffx.potential.parameters.ForceField} object.
     * @param forceFieldEnergy a {@link ffx.potential.ForceFieldEnergy} object.
     * @param mola             a {@link ffx.potential.MolecularAssembly} object.
     */
    public MultiTerminus(Residue residue, ForceField forceField, ForceFieldEnergy forceFieldEnergy,
                         MolecularAssembly mola) {
        super(residue.getName(), residue.getResidueNumber(), residue.residueType,
                residue.getChainID(), residue.getChainID().toString());
        try {
            String ffString = forceField.getString(ForceField.ForceFieldString.FORCEFIELD);
            if (ForceField.ForceFieldName.valueOf(ffString) != ForceField.ForceFieldName.AMOEBA_PROTEIN_2013) {
                logger.severe("MultiTerminus supported only under AMOEBA_PROTEIN_2013.");
            }
        } catch (Exception ex) {
            Logger.getLogger(MultiTerminus.class.getName()).log(Level.SEVERE, null, ex);
        }
        this.forceField = forceField;
        this.forceFieldEnergy = forceFieldEnergy;
        this.mola = mola;
        if (residue.getNextResidue() == null) {
            end = END.CTERM;
            isCharged = (residue.getAtomNode("HO") == null);
        } else if (residue.getPreviousResidue() == null) {
            end = END.NTERM;
            isCharged = (residue.getAtomNode("H3") != null);
        } else {
            end = null;
            logger.severe("MultiTerminus constructed for non-terminal residue.");
        }
//        removeLeaves();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void revertState(ResidueState state) {
        if (state.getIsNeutralTerminus() && this.isCharged == true) {
            titrateTerminus_v1(298.15);     // TODO generalize
        }
        super.revertState(state);
    }

    /**
     * Useful for locating backbone atom nodes that share a name with side-chain atoms.
     */
    private Atom getBBAtom(String name) {
        ArrayList<Atom> list = this.getAtomList();
        for (Atom atom : list) {
            if (atom.getName().equals(name)) {
//                logger.info(" Found: " + atom.getName());
                return atom;
//                try {
//                    BB_TYPE bbType = BB_TYPE.valueOf(name);
//                    AtomType type = atom.getAtomType();
//                    if ((type.type == bbType.chrgType && type.atomClass == bbType.chrgClass)
//                            || (type.type == bbType.neutType && type.atomClass == bbType.neutClass)) {
//                        return atom;
//                    }
//                } catch (Exception ex) {}
            }
        }
        return null;
    }

    private void updateBondedTerms() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("Updating bonded terms: \n"));
        List<ROLS> bonds = getBondList();
        for (int i = 0; i < bonds.size(); i++) {
            Bond bond = (Bond) bonds.get(i);
            BondType oldType = bond.bondType;
            int c[] = new int[2];
            c[0] = bond.atoms[0].getAtomType().atomClass;
            c[1] = bond.atoms[1].getAtomType().atomClass;
            String key = BondType.sortKey(c);
            BondType newType = forceField.getBondType(key);
            if (oldType != newType) {
                sb.append(String.format(" Bond: %s --> %s \n", bond.bondType, newType));
                bond.setBondType(newType);
                if (newType.distance < 0.9 * oldType.distance || newType.distance > 1.1 * oldType.distance) {
                    logger.info(String.format(" Large bond distance change: %s %s,  %.2f --> %.2f ",
                            bond.atoms[0].describe(Atom.Descriptions.XyzIndex_Name), bond.atoms[1].describe(Atom.Descriptions.XyzIndex_Name),
                            oldType.distance, newType.distance));
                }
            }
        }
        List<ROLS> angles = getAngleList();
        for (int i = 0; i < angles.size(); i++) {
            Angle angle = (Angle) angles.get(i);
            AngleType oldType = angle.angleType;
            if (DEBUG) {
                logger.info(String.format(" %d ( %s %s %s ) ( %d %d %d )", i,
                        angle.atoms[0].describe(Atom.Descriptions.XyzIndex_Name),
                        angle.atoms[1].describe(Atom.Descriptions.XyzIndex_Name),
                        angle.atoms[2].describe(Atom.Descriptions.XyzIndex_Name),
                        angle.atoms[0].getAtomType().atomClass,
                        angle.atoms[1].getAtomType().atomClass,
                        angle.atoms[2].getAtomType().atomClass));
            }
            Angle dummy = Angle.angleFactory(angle.bonds[0], angle.bonds[1], forceField);
            AngleType newType = dummy.angleType;
            if (oldType != newType) {
                sb.append(String.format(" Angle: %s --> %s \n", angle.angleType, dummy.angleType));
                angle.setAngleType(dummy.angleType);
                if (newType.angle[0] < 0.9 * oldType.angle[0] || newType.angle[0] > 1.1 * oldType.angle[0]) {
                    logger.info(String.format(" Large angle change: %s %s %s,  %.2f --> %.2f ",
                            angle.atoms[0].describe(Atom.Descriptions.XyzIndex_Name), angle.atoms[1].describe(Atom.Descriptions.XyzIndex_Name), angle.atoms[2].describe(Atom.Descriptions.XyzIndex_Name),
                            oldType.angle[0], newType.angle[0]));
                }
            }
        }
        List<ROLS> torsions = getTorsionList();
        for (int i = 0; i < torsions.size(); i++) {
            Torsion tors = (Torsion) torsions.get(i);
            TorsionType oldType = tors.torsionType;
            Torsion dummy = Torsion.torsionFactory(tors.bonds[0], tors.bonds[1], tors.bonds[2], forceField);
            TorsionType newType = dummy.torsionType;
            if (oldType != newType) {
                sb.append(String.format(" Torsion: %s --> %s \n", tors.torsionType, dummy.torsionType));
                tors.torsionType = dummy.torsionType;
            }
        }
        if (DEBUG) {
            logger.info(sb.toString());
        }
    }

    private Atom uberH3;
    private Atom uberHO;
    private Bond bondH3;
    private Bond bondHO;
    private List<ROLS> rolsH3 = new ArrayList<>();
    private List<ROLS> rolsHO = new ArrayList<>();

    /**
     * Changes the charge state of this MultiTerminus.
     * Keep existing Atom objects but updates types, bonded terms, and builds new proton if necessary.
     *
     * @param temperature a double.
     * @return a {@link ffx.potential.extended.TitrationUtils.TitrationType} object.
     */
    public TitrationType titrateTerminus_v1(double temperature) {
        logger.info(String.format(" Titrating residue %s (currently %d).", this.toString(), (isCharged ? 1 : 0)));
//        StringBuilder sb = new StringBuilder();
//        sb.append(" Contents of children: ");
//        for (Atom atom : this.getAtomList()) {
//            sb.append(String.format("%s, ", atom.getName()));
//        }
//        logger.info(sb.toString());
        /**
         * Get references to the backbone atoms.
         */
        TitrationType titrationType = (isCharged) ? TitrationType.DEPROT : TitrationType.PROT;
        Atom N = getBBAtom("N");
        Atom CA = getBBAtom("CA");
        Atom C = getBBAtom("C");
        Atom O = getBBAtom("O");
        Atom OXT = getBBAtom("OXT");
        Atom H1 = getBBAtom("H1");
        Atom H2 = getBBAtom("H2");
        Atom H3 = getBBAtom("H3");
        Atom HA = getBBAtom("HA");
        Atom OH = getBBAtom("OH");
        Atom HO = getBBAtom("HO");
        String resName = C.getResidueName();
        int resSeq = C.getResidueNumber();
        Character chainID = C.getChainID();
        List<Atom> typeChanged = new ArrayList<>();
        if (DEBUG) {
            printBonds();
        }

        if (end == END.NTERM) {
            if (isCharged) {
                if (rolsH3.isEmpty()) {
                    for (Angle a : H3.getAngles()) {
                        rolsH3.add(a);
                    }
                    for (Torsion t : H3.getTorsions()) {
                        rolsH3.add(t);
                    }
                }
                N.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.N.neutType)));
                H1.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H1.neutType)));
                H2.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H2.neutType)));
//                intxyz(H1, N, 1.02, CA, 109.5, C, 180.0, 0);
//                intxyz(H2, N, 1.02, CA, 109.5, C, 0.0, 0);
                bondH3 = N.getBond(H3);
                bondH3.removeFromParent();
                for (Angle a : H3.getAngles()) {
                    a.removeFromParent();
                }
                for (Torsion t : H3.getTorsions()) {
                    t.removeFromParent();
                }
                H3.removeFromParent();
                H3.setParent(null);
                H3.setUse(false);
                uberH3 = H3;
                typeChanged.add(N);
                typeChanged.add(H1);
                typeChanged.add(H2);
//                logger.info(String.format(" Finished titration. H3 status: %b %b %b",
//                        N.getBond(H3) == null, this.getAtomNode().contains(H3) == null, H3.getParent() == null));
            } else {
                if (H3 != null) {
                    logger.severe("N-terminal found in incorrect charge state.");
                }
                if (uberH3 == null || bondH3 == null) {
                    logger.severe("Please start with (N-)termini in the charged state.");
                }
                H3 = uberH3;
                N.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.N.chrgType)));
                H1.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H1.chrgType)));
                H2.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H2.chrgType)));
                intxyz(H1, N, 1.02, CA, 109.5, C, 180.0, 0);
                intxyz(H2, N, 1.02, CA, 109.5, C, 60.0, 0);
                intxyz(H3, N, 1.02, CA, 109.5, C, -60.0, 0);
                maxwellMe(H3, temperature);
                this.getAtomNode().add(H3);
                this.add(bondH3);
//                logger.info("Number of rolsH3 terms: " + rolsH3.size());
                for (ROLS rols : rolsH3) {
                    if (rols instanceof Angle) {
                        this.add((Angle) rols);
                    } else if (rols instanceof Torsion) {
                        this.add((Torsion) rols);
                    }
                }
                H3.setParent(this.getAtomNode());
                H3.setUse(true);
                typeChanged.add(N);
                typeChanged.add(H1);
                typeChanged.add(H2);
                typeChanged.add(H3);
//                logger.info(String.format(" Finished titration. H3 statuses: "
//                        + "(They have each other: %b %b) (I have them: %b %b) (I have bond: %b)", 
//                        N.getBond(H3) != null, H3.getBond(N) != null, 
//                        this.getAtomNode().contains(H3) != null, H3.getParent() == this.getAtomNode(),
//                        this.getBondList().contains(bondH3)));
//                logger.info(String.format(" Bonds from H3: %s %s",
//                        H3.getBonds().get(0).get1_2(H3).describe(Atom.Descriptions.INDEX_NAME), 
//                        H3.getBonds().get(0).get1_2(H3).getBonds().get(0).get1_2(H3.getBonds().get(0).get1_2(H3)).describe(Atom.Descriptions.INDEX_NAME)));
            }
        } else if (end == END.CTERM) {
            if (isCharged) {
                OXT.setName("OH");
                OH = OXT;
                if (HO != null) {
                    logger.warning("C-terminal in unusual charge state.");
                }
                C.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.C.neutType)));
                O.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.O.neutType)));
                OH.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.OH.neutType)));
                if (uberHO == null) {
                    // Gotta build the HO and all its bonded terms.
                    uberHO = new Atom(mola.getAtomArray().length, "HO", OH.getAltLoc(), new double[3],
                            resName, resSeq, chainID, OH.getOccupancy(), OH.getTempFactor(), OH.getSegID(), true);
                    uberHO.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.HO.neutType)));
                    bondHO = new Bond(OH, uberHO);
                    bondHO.setBondType(forceField.getBondType(
                            String.format("%d %d", OH.getAtomType().atomClass, uberHO.getAtomType().atomClass)));
                    Angle a1 = Angle.angleFactory(bondHO, OH.getBond(C), forceField);
                    Torsion t1 = Torsion.torsionFactory(O.getBond(C), C.getBond(OH), bondHO, forceField);
                    Torsion t2 = Torsion.torsionFactory(CA.getBond(C), C.getBond(OH), bondHO, forceField);
                    this.add(a1);
                    this.add(t1);
                    this.add(t2);
                }
                HO = uberHO;
                intxyz(HO, OXT, 1.02, C, 109.5, CA, -1.7, 0);
                maxwellMe(HO, temperature);
                this.getAtomNode().add(HO);
                this.add(bondHO);
                HO.setParent(this.getAtomNode());
                HO.setUse(true);
//                logger.info("Number of rolsHO terms: " + rolsHO.size());
                for (ROLS rols : rolsHO) {
                    if (rols instanceof Angle) {
                        this.add((Angle) rols);
                    } else if (rols instanceof Torsion) {
                        this.add((Torsion) rols);
                    }
                }
                typeChanged.add(C);
                typeChanged.add(O);
                typeChanged.add(OH);
                typeChanged.add(HO);
            } else {
                if (rolsHO.isEmpty()) {
                    rolsHO = new ArrayList<>();
                    for (Angle a : HO.getAngles()) {
                        rolsHO.add(a);
                    }
                    for (Torsion t : HO.getTorsions()) {
                        rolsHO.add(t);
                    }
                }
                C.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.C.chrgType)));
                O.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.O.chrgType)));
                OH.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.OH.chrgType)));
                OH.setName("OXT");
                OXT = OH;
                bondHO = OH.getBond(HO);
                bondHO.removeFromParent();
                for (Angle a : HO.getAngles()) {
                    a.removeFromParent();
                }
                for (Torsion t : HO.getTorsions()) {
                    t.removeFromParent();
                }
                HO.removeFromParent();
                HO.setParent(null);
                HO.setUse(false);
                uberHO = HO;
                typeChanged.add(C);
                typeChanged.add(O);
                typeChanged.add(OH);
            }
        }
        updateGeometry();
        updateBondedTerms();
        isCharged = !isCharged;
        forceFieldEnergy.reInit();
        if (DEBUG) {
            printBonds();
        }
        return titrationType;
    }

    private void printBonds() {
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("BondList: %s\n", this.toString()));
        for (ROLS rols : getBondList()) {
            Bond bond = (Bond) rols;
            sb.append(String.format("  %s\n", bond.toString()));
        }
        logger.info(sb.toString());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void add(MutableTreeNode mtn) {
        super.add(mtn);
        if (DEBUG) {
            logger.info("Adding to terminus: " + mtn.toString());
        }
    }

    /**
     * Changes the charge state of this MultiTerminus.
     * Keep existing Atom objects but updates types, bonded terms, and builds new proton if necessary.
     */
    public void titrateTerminus_v0() {
        logger.info(String.format(" Titrating residue %s (currently %d).", this.toString(), (isCharged ? 1 : 0)));
        StringBuilder sb = new StringBuilder();
        sb.append(" Contents of children: ");
        for (Atom atom : this.getAtomList()) {
            sb.append(String.format("%s, ", atom.getName()));
        }
        logger.info(sb.toString());
        /**
         * Get references to the backbone atoms.
         */
        Atom N = getBBAtom("N");
        Atom CA = getBBAtom("CA");
        Atom C = getBBAtom("C");
        Atom O = getBBAtom("O");
        Atom OXT = getBBAtom("OXT");
        Atom H1 = getBBAtom("H1");
        Atom H2 = getBBAtom("H2");
        Atom H3 = getBBAtom("H3");
        Atom HA = getBBAtom("HA");
        Atom OH = getBBAtom("OH");
        Atom HO = getBBAtom("HO");
        String resName = C.getResidueName();
        int resSeq = C.getResidueNumber();
        Character chainID = C.getChainID();

        if (end == END.NTERM) {
            if (isCharged) {
                N.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.N.neutType)));
                H1.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H1.neutType)));
                H2.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H2.neutType)));
                intxyz(H1, N, 1.02, CA, 109.5, C, 180.0, 0);
                intxyz(H2, N, 1.02, CA, 109.5, C, 0.0, 0);
                Bond bondH3 = N.getBond(H3);
                bondH3.removeFromParent();
                H3.removeFromParent();
                H3.setParent(null);
//                N.getBondList().remove(bondH3);       // oughtta be in updateGeometry()
//                N.removeBond(bondH3);
//                this.getBondList().remove(bondH3);    // returns false
//                this.getAtomNode().remove(H3);        // throws "is not a child"
                updateGeometry();
                logger.info(String.format(" Finished titration. H3 status: %b %b %b",
                        N.getBond(H3) == null, this.getAtomNode().contains(H3) == null, H3.getParent() == null));
            } else {
                if (H3 != null) {
                    logger.warning("N-terminal use state toggled.");
                }
                N.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.N.chrgType)));
                H1.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H1.chrgType)));
                H2.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H2.chrgType)));
                H3 = new Atom(mola.getAtomArray().length, "H3", N.getAltLoc(), new double[3], resName, resSeq, chainID,
                        N.getOccupancy(), N.getTempFactor(), N.getSegID(), true);
                H3.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.H3.chrgType)));
                intxyz(H1, N, 1.02, CA, 109.5, C, 180.0, 0);
                intxyz(H2, N, 1.02, CA, 109.5, C, 60.0, 0);
                intxyz(H3, N, 1.02, CA, 109.5, C, -60.0, 0);
                double vel[] = new double[3];
                N.getVelocity(vel);
                H3.setVelocity(vel);
//                Bond bondH3 = buildBond(N, H3, forceField, null); // try manually
                Bond bondH3 = new Bond(N, H3);
                bondH3.setBondType(forceField.getBondType(
                        String.format("%d %d", N.getAtomType().atomClass, H3.getAtomType().atomClass)));
                this.getAtomNode().add(H3);
                this.getBondList().add(bondH3);
                this.add(bondH3);
                updateGeometry();
                logger.info(String.format(" Finished titration. H3 statuses: "
                                + "(They have each other: %b %b) (I have them: %b %b) (I have bond: %b)",
                        N.getBond(H3) != null, H3.getBond(N) != null,
                        this.getAtomNode().contains(H3) != null, H3.getParent() == this.getAtomNode(),
                        this.getBondList().contains(bondH3)));
                logger.info(String.format(" Bonds from H3: %s %s",
                        H3.getBonds().get(0).get1_2(H3).describe(Atom.Descriptions.XyzIndex_Name),
                        H3.getBonds().get(0).get1_2(H3).getBonds().get(0).get1_2(H3.getBonds().get(0).get1_2(H3)).describe(Atom.Descriptions.XyzIndex_Name)));
            }
        } else if (end == END.CTERM) {
            if (isCharged) {
                if (HO != null) {
                    logger.warning("C-terminal in unusual charge state.");
                }
                C.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.C.neutType)));
                O.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.O.neutType)));
                OXT.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.OH.neutType)));
                OXT.setName("OH");
                OH = OXT;
                HO = new Atom(mola.getAtomArray().length, "HO", OXT.getAltLoc(), new double[3], resName, resSeq, chainID,
                        OXT.getOccupancy(), OXT.getTempFactor(), OXT.getSegID(), true);
                HO.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.HO.neutType)));
                intxyz(HO, OXT, 1.02, C, 109.5, CA, -1.7, 0);
                double vel[] = new double[3];
                OH.getVelocity(vel);
                HO.setVelocity(vel);
                buildBond(OH, HO, forceField, null);
                this.getAtomNode().add(HO);
            } else {
                C.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.C.chrgType)));
                O.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.O.chrgType)));
                OH.setAtomType(forceField.getAtomType(Integer.toString(BB_TYPE.OH.chrgType)));
                OH.setName("OXT");
                Bond bondHO = OH.getBond(HO);
                bondHO.removeFromParent();
                HO.removeFromParent();
                HO.setParent(null);
                this.getBondList().remove(bondHO);
                this.getAtomNode().remove(HO);
                updateGeometry();
            }
        }
        isCharged = !isCharged;
        forceFieldEnergy.reInit();
    }

    /**
     * For testing.
     */
    private void titrateTerminusByRebuilding() {
        if (true) {
            throw new UnsupportedOperationException();
        }
        /**
         * Get references to the backbone atoms.
         */
        Atom CA = (Atom) this.getAtomNode("CA");
        Atom C = (Atom) this.getAtomNode("C");
        Atom HA = (Atom) this.getAtomNode("HA");
        Atom N = (Atom) this.getAtomNode("N");
        Atom O = (Atom) this.getAtomNode("O");
        Atom OXT = (Atom) this.getAtomNode("OXT");
        Atom NH2 = (Atom) this.getAtomNode("NH2");

        if (getNextResidue() == null) {

        } else if (getPreviousResidue() == null) {
            String resName = C.getResidueName();
            int resSeq = C.getResidueNumber();
            Character chainID = C.getChainID();
            double Cxyz[] = new double[3];
            double Oxyz[] = new double[3];
            double OXTxyz[] = new double[3];
            C.getXYZ(Cxyz);
            O.getXYZ(Oxyz);
            OXT.getXYZ(OXTxyz);

            int protCkey = 235;
            int protOkey = 236;
            int protOHkey = 237;
            int protHOkey = 238;
            Atom protC = new Atom(0, "C", C.getAltLoc(), Cxyz, resName, resSeq, chainID,
                    C.getOccupancy(), C.getTempFactor(), C.getSegID(), true);
            Atom protO = new Atom(0, "O", O.getAltLoc(), Oxyz, resName, resSeq, chainID,
                    O.getOccupancy(), O.getTempFactor(), O.getSegID(), true);
            Atom protOH = new Atom(0, "OH", OXT.getAltLoc(), OXTxyz, resName, resSeq, chainID,
                    OXT.getOccupancy(), OXT.getTempFactor(), OXT.getSegID(), true);
            Atom protHO = new Atom(0, "HO", OXT.getAltLoc(), new double[3], resName, resSeq, chainID,
                    OXT.getOccupancy(), OXT.getTempFactor(), OXT.getSegID(), true);
            intxyz(protHO, protOH, 1.02, protC, 109.5, CA, -1.7, 0);
            protC.setAtomType(forceField.getAtomType(Integer.toString(protCkey)));
            protO.setAtomType(forceField.getAtomType(Integer.toString(protOkey)));
            protOH.setAtomType(forceField.getAtomType(Integer.toString(protOHkey)));
            protHO.setAtomType(forceField.getAtomType(Integer.toString(protHOkey)));

            buildBond(CA, protO, forceField, null);
            buildBond(protC, protO, forceField, null);
            buildBond(protC, protOH, forceField, null);
            buildBond(protOH, protHO, forceField, null);
        }
    }

    /**
     * Update Atom references to local geometry.
     *
     * @param residue
     */
    private void updateGeometry() {
        /**
         * Update atom references to local geometry.
         */
        ArrayList<Atom> atoms = this.getAtomList();
        ArrayList<ROLS> bonds = this.getBondList();
        ArrayList<ROLS> angles = this.getAngleList();
        ArrayList<ROLS> torsions = this.getTorsionList();

        for (Atom atom : atoms) {
            atom.clearGeometry();
        }
        for (Atom atom : atoms) {
            for (ROLS bond : bonds) {
                Bond b = (Bond) bond;
                if (b.containsAtom(atom)) {
                    atom.setBond(b);
                }
            }
        }
        for (Atom atom : atoms) {
            for (ROLS angle : angles) {
                Angle a = (Angle) angle;
                if (a.containsAtom(atom)) {
                    atom.setAngle(a);
                }
            }
        }
        for (Atom atom : atoms) {
            for (ROLS torsion : torsions) {
                Torsion t = (Torsion) torsion;
                if (t.containsAtom(atom)) {
                    atom.setTorsion(t);
                }
            }
        }
    }

//  For reference.
//    public void addResidue(Residue newResidue) {
//        /**
//         * Add the new residue to list.
//         */
//        consideredResidues.add(newResidue);
//        /**
//         * Get references to nearby residues.
//         */
//        Residue prevResidue = activeResidue.getPreviousResidue();
//        Residue nextResidue = activeResidue.getNextResidue();
//        Residue prev2Residue = null;
//        if (prevResidue != null) {
//            prev2Residue = prevResidue.getPreviousResidue();
//        }
//        Residue next2Residue = null;
//        if (nextResidue != null) {
//            next2Residue = nextResidue.getNextResidue();
//        }
//
//        /**
//         * Move atoms from the active Residue to the new Residue.
//         */
//        moveBackBoneAtoms(activeResidue, newResidue);
//        
//        /**
//         * Pass references of the active Residues' joints to the new Residue.
//         */
//        ArrayList<Joint> joints = activeResidue.getJoints();
//        for (Joint joint : joints) {
//            newResidue.addJoint(joint);
//        }
//        /**
//         * Make the new Residue active.
//         */
//        activeResidue.removeFromParent();
//        activeResidue = newResidue;
//        add(activeResidue);
//        /**
//         * Build side-chain atoms and assign atom types for the new Residue.
//         */
//        try {
//            assignAminoAcidAtomTypes(newResidue, prevResidue, nextResidue, forceField, null);
//            if (nextResidue != null) {
//                Atom C = (Atom) newResidue.getAtomNode("C");
//                Atom nextN = (Atom) nextResidue.getAtomNode("N");
//                for (Joint joint : joints) {
//                    Bond bond = (Bond) joint.getBondList().get(0);
//                    if (bond.containsAtom(C) && bond.containsAtom(nextN)) {
//                        C.setBond(bond);
//                    }
//                }
//            }
//        } catch (MissingHeavyAtomException | MissingAtomTypeException exception) {
//            logger.severe(exception.toString());
//        }
//        newResidue.finalize(true, forceField);
//        updateGeometry(newResidue, prevResidue, nextResidue, prev2Residue, next2Residue);
//    }

    /**
     * <p>isCharged.</p>
     *
     * @return a boolean.
     */
    public boolean isCharged() {
        return isCharged;
    }

    public enum BB_TYPE {
        /**
         * [name from Biotype record appended]
         * atom    233     30  C     "C-Terminal COO-"            6   12.0110  3   C
         * atom    234     31  O     "C-Terminal COO-"            8   15.9990  1   OXT
         * atom    235     32  C     "C-Terminal COOH C=O"        6   12.0110  3   C
         * atom    236     33  O     "C-Terminal COOH O=C"        8   15.9990  1   O
         * atom    237     34  OH    "C-Terminal COOH OH"         8   15.9990  2   OH
         * atom    238     35  HO    "C-Terminal COOH HO"         1    1.0080  1   HO
         * atom    225      1  N     "Amide Cap NH2"              7   14.0070  3   N
         * atom    226      4  HN    "Amide Cap H2N"              1    1.0080  1   HN
         * atom    231     41  N     "N-Terminal NH3+"            7   14.0070  4   N
         * atom    232     42  H     "N-Terminal H3N+"            1    1.0080  1   HN
         **/
        N(231, 41, 225, 1),
        C(233, 30, 235, 32),
        O(234, 31, 236, 33),
        OXT(234, 31, 236, 33),     // On titration, this'll need renamed.
        OH(234, 31, 237, 34),      // ^
        HO(-1, -1, 238, 35),
        H1(232, 41, 226, 4),       // On titration, these'll need removed and rebuilt.
        H2(232, 41, 226, 4),
        H3(232, 41, 226, 4);

        public int chrgType, chrgClass;
        public int neutType, neutClass;

        BB_TYPE(int chrgType, int chrgClass, int neutType, int neutClass) {
            this.chrgType = chrgType;
            this.chrgClass = chrgClass;
            this.neutType = neutType;
            this.neutClass = neutClass;
        }
    }

    public enum END {
        NTERM, CTERM;
    }

    private void maxwellMe(Atom atom, double temperature) {
        double vv[] = new double[3];
        for (int i = 0; i < 3; i++) {
            vv[i] = ThreadLocalRandom.current().nextGaussian() * sqrt(kB * temperature / atom.getMass());
        }
        atom.setVelocity(vv);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        MultiTerminus other = (MultiTerminus) object;
        if (this.getResidueNumber() == other.getResidueNumber()
                && this.isCharged == other.isCharged()) {
            return true;
        } else {
            return false;
        }
    }

}
