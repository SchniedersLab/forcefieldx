/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
 */
package ffx.potential.bonded;

import ffx.potential.ResidueEnumerations.NucleicAcid3;
import ffx.potential.Rotamer;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.Hashtable;
import java.util.ListIterator;
import java.util.logging.Logger;

import javax.media.j3d.Canvas3D;
import javax.media.j3d.J3DGraphics2D;
import javax.media.j3d.Material;
import javax.media.j3d.Node;
import javax.vecmath.Color3f;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import ffx.potential.RotamerLibrary;

import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;
import java.util.logging.Level;

/**
 * The Residue class represents individual amino acids or nucleic acid bases.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 */
public class Residue extends MSGroup {

    private static final Logger logger = Logger.getLogger(ffx.potential.bonded.Residue.class.getName());
    private static final long serialVersionUID = 1L;
    private static Point3d point3d = new Point3d();
    private static Point2d point2d = new Point2d();
    /**
     * Constant <code>NA1Set</code>
     */
    public static EnumSet NA1Set = EnumSet.allOf(NA1.class);
    /**
     * Constant <code>NA3Set</code>
     */
    public static EnumSet NA3Set = EnumSet.allOf(NA3.class);
    /**
     * Constant <code>NASet</code>
     */
    public static EnumSet NASet = EnumSet.allOf(NA.class);
    /**
     * Constant <code>NA1toNA3</code>
     */
    public static Hashtable<NA1, NA3> NA1toNA3 = new Hashtable<NA1, NA3>();
    /**
     * Constant <code>NA3Color</code>
     */
    public static Hashtable<NA3, Color3f> NA3Color = new Hashtable<NA3, Color3f>();
    /**
     * Constant <code>AA1Set</code>
     */
    public static EnumSet AA1Set = EnumSet.allOf(AA1.class);
    /**
     * Constant <code>AA3Set</code>
     */
    public static EnumSet AA3Set = EnumSet.allOf(AA3.class);
    /**
     * Constant <code>AASet</code>
     */
    public static EnumSet AASet = EnumSet.allOf(AA.class);
    /**
     * Constant <code>AA1toAA3</code>
     */
    public static Hashtable<AA1, AA3> AA1toAA3 = new Hashtable<AA1, AA3>();
    /**
     * Constant <code>AA3Color</code>
     */
    public static Hashtable<AA3, Color3f> AA3Color = new Hashtable<AA3, Color3f>();
    /**
     * Constant <code>SSTypeColor</code>
     */
    public static Hashtable<SSType, Color3f> SSTypeColor = new Hashtable<SSType, Color3f>();
    /**
     * Constant <code>Ramachandran="new String[17]"</code>
     */
    public static String Ramachandran[] = new String[17];

    static {
        NA1 na1[] = NA1.values();
        NA3 na3[] = NA3.values();
        for (int i = 0; i < NA1.values().length; i++) {
            NA1toNA3.put(na1[i], na3[i]);
        }
    }

    static {
        NA3Color.put(NA3.A, RendererCache.RED);
        NA3Color.put(NA3.C, RendererCache.MAGENTA);
        NA3Color.put(NA3.G, RendererCache.BLUE);
        NA3Color.put(NA3.U, RendererCache.YELLOW);
        NA3Color.put(NA3.DA, RendererCache.RED);
        NA3Color.put(NA3.DC, RendererCache.MAGENTA);
        NA3Color.put(NA3.DG, RendererCache.BLUE);
        NA3Color.put(NA3.DT, RendererCache.ORANGE);
        NA3Color.put(NA3.MPO, RendererCache.GREEN);
        NA3Color.put(NA3.DPO, RendererCache.GREEN);
        NA3Color.put(NA3.TPO, RendererCache.GREEN);
        NA3Color.put(NA3.UNK, RendererCache.CYAN);
    }

    static {
        AA1 aa1[] = AA1.values();
        AA3 aa3[] = AA3.values();
        for (int i = 0; i < AA1.values().length; i++) {
            AA1toAA3.put(aa1[i], aa3[i]);
        }
    }

    static {
        AA3Color.put(AA3.ALA, RendererCache.GRAY);
        AA3Color.put(AA3.ARG, RendererCache.BLUE);
        AA3Color.put(AA3.ASN, RendererCache.BLUE);
        AA3Color.put(AA3.ASP, RendererCache.RED);
        AA3Color.put(AA3.CYS, RendererCache.YELLOW);
        AA3Color.put(AA3.GLN, RendererCache.BLUE);
        AA3Color.put(AA3.GLU, RendererCache.RED);
        AA3Color.put(AA3.GLY, RendererCache.GRAY);
        AA3Color.put(AA3.ILE, RendererCache.GRAY);
        AA3Color.put(AA3.LEU, RendererCache.GRAY);
        AA3Color.put(AA3.LYS, RendererCache.BLUE);
        AA3Color.put(AA3.MET, RendererCache.YELLOW);
        AA3Color.put(AA3.PHE, RendererCache.GREEN);
        AA3Color.put(AA3.PRO, RendererCache.ORANGE);
        AA3Color.put(AA3.SER, RendererCache.BLUE);
        AA3Color.put(AA3.THR, RendererCache.BLUE);
        AA3Color.put(AA3.TRP, RendererCache.GREEN);
        AA3Color.put(AA3.TYR, RendererCache.GREEN);
        AA3Color.put(AA3.VAL, RendererCache.GRAY);
        AA3Color.put(AA3.HIS, RendererCache.BLUE);
        AA3Color.put(AA3.HIE, RendererCache.BLUE);
        AA3Color.put(AA3.HID, RendererCache.BLUE);
        AA3Color.put(AA3.ORN, RendererCache.ORANGE);
        AA3Color.put(AA3.AIB, RendererCache.ORANGE);
        AA3Color.put(AA3.PCA, RendererCache.ORANGE);
        AA3Color.put(AA3.FOR, RendererCache.RED);
        AA3Color.put(AA3.ACE, RendererCache.RED);
        AA3Color.put(AA3.NH2, RendererCache.BLUE);
        AA3Color.put(AA3.NME, RendererCache.BLUE);
        AA3Color.put(AA3.UNK, RendererCache.MAGENTA);
    }

    static {
        SSTypeColor.put(SSType.NONE, RendererCache.WHITE);
        SSTypeColor.put(SSType.SHEET, RendererCache.PINK);
        SSTypeColor.put(SSType.HELIX, RendererCache.BLUE);
        SSTypeColor.put(SSType.TURN, RendererCache.YELLOW);
    }

    static {
        Ramachandran[0] = "Default (Extended)       [-135.0  135.0]";
        Ramachandran[1] = "Alpha Helix (R)          [ -57.0  -47.0]";
        Ramachandran[2] = "Alpha Helix (L)          [  57.0   47.0]";
        Ramachandran[3] = "3-10 Helix               [ -49.0  -26.0]";
        Ramachandran[4] = "Pi Helix                 [ -57.0  -70.0]";
        Ramachandran[5] = "Polyproline II Helix     [ -79.0  149.0]";
        Ramachandran[6] = "Parallel Beta Strand     [-119.0  113.0]";
        Ramachandran[7] = "Antiparallel Beta Strand [-139.0  135.0]";
        Ramachandran[8] = "Beta-Hairpin 2' (i+1)    [  90.0 -170.0]";
        Ramachandran[9] = "Beta-Hairpin 2' (i+2)    [ -80.0  -10.0]";
        Ramachandran[10] = "Beta-Hairpin 1' (i+1)    [  57.0   47.0]";
        Ramachandran[11] = "Beta-Hairpin 1' (i+2)    [  57.0   47.0]";
        Ramachandran[12] = "Beta-Hairpin 1  (i+1)    [ -57.0  -47.0]";
        Ramachandran[13] = "Beta-Hairpin 1  (i+2)    [ -57.0  -47.0]";
        Ramachandran[14] = "Beta-Hairpin 1  (i+3)    [  90.0 -170.0]";
        Ramachandran[15] = "Beta-Hairpin 3' (i+1)    [  57.0   47.0]";
        Ramachandran[16] = "Beta-Hairpin 3' (i+2)    [ -80.0  -10.0]";
    }
    /**
     * The residue number of this atom in a chain.
     */
    private int resNumber;
    /**
     * Possibly redundant PDB chain ID.
     */
    private Character chainID;
    /**
     * Unique segID.
     */
    private String segID;
    /**
     * Secondary structure type.
     */
    private SSType ssType = SSType.NONE;
    /**
     * Residue type.
     */
    protected ResidueType residueType = ResidueType.UNK;
    private AA3 aa;
    private NA3 na;
    /**
     * These arrays store default coordinates for certain atoms in nucleic acid
     * Residues. C1', O4', and C4' are the critical sugar atoms off which every
     * other atom is drawn when applyRotamer is called; the backbone
     * corrections, <<<<<<< HEAD however, move these atoms, so they must be
     * reverted to original coordinates each time applyRotamer is called.
     *
     * ======= however, move these atoms, so they must be reverted to original
     * coordinates each time applyRotamer is called.
     *
     * >>>>>>> fcefa7e8f2a22b42ca11a7021b7f01ddba9a6e2f O3' North and South
     * coordinates are technically non-essential, as they could be derived from
     * C1', O4', C4', and a given sugar pucker, however, it is much less
     * computationally expensive to calculate them once and then store them.
     * <<<<<<< HEAD
     *
     * =======
     *
     * >>>>>>> fcefa7e8f2a22b42ca11a7021b7f01ddba9a6e2f TODO: Add O3'
     * coordinates for the DNA C3'-exo configuration.
     */
    private double[] O3sNorthCoords = null;
    private double[] O3sSouthCoords = null;
    private double[] C1sCoords = null;
    private double[] O4sCoords = null;
    private double[] C4sCoords = null;

    /**
     * Default Constructor where num is this Residue's position in the Polymer.
     *
     * @param num a int.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     */
    public Residue(int num, ResidueType rt) {
        super();
        resNumber = num;
        residueType = rt;
        assignResidueType();
    }

    /**
     * <p>
     * Constructor for Residue.</p>
     *
     * @param name a {@link java.lang.String} object.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     */
    public Residue(String name, ResidueType rt) {
        super(name);
        residueType = rt;
        assignResidueType();
    }

    /**
     * Name is the residue's 3 letter abbreviation and num is its position in
     * the Polymer.
     *
     * @param name a {@link java.lang.String} object.
     * @param num a int.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     */
    public Residue(String name, int num, ResidueType rt) {
        this(name, rt);
        resNumber = num;
    }

    /**
     * Name is the residue's 3 letter abbreviation and num is its position in
     * the Polymer.
     *
     * @param name a {@link java.lang.String} object.
     * @param resNumber a int.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     * @param chainID a {@link java.lang.Character} object.
     * @param segID a {@link java.lang.String} object.
     */
    public Residue(String name, int resNumber, ResidueType rt, Character chainID,
            String segID) {
        this(name, rt);
        this.resNumber = resNumber;
        this.chainID = chainID;
        this.segID = segID;
    }

    /**
     * As above, with atoms being a FNode with this Residue's atoms as child
     * nodes
     *
     * @param name a {@link java.lang.String} object.
     * @param num a int.
     * @param atoms a {@link ffx.potential.bonded.MSNode} object.
     * @param rt a {@link ffx.potential.bonded.Residue.ResidueType} object.
     */
    public Residue(String name, int num, MSNode atoms, ResidueType rt) {
        super(name, atoms);
        resNumber = num;
        residueType = rt;
        assignResidueType();
        finalize(true);
    }

    public Rotamer[] getRotamers(Residue residue) {
        return RotamerLibrary.getRotamers(residue);
    }

    public ResidueType getResidueType() {
        return residueType;
    }

    /**
     * Returns the Residue bonded to this Residue at this Residue's 3' or 
     * C-terminal end.  Any use of this method to add Residues to a sliding
     * window or similar MUST not add that residue if that residue has no 
     * Rotamers, as several algorithms (such as the distance matrix) assume
     * that all Residues being optimized have Rotamers.
     * 
     * @return The next Residue.
     */
    public Residue getNextResidue() {
        switch (residueType) {
            case AA: {
                Atom carbon = (Atom) getAtomNode("C");
                ArrayList<Bond> bonds = carbon.getBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(carbon);
                    if (other.getName().equalsIgnoreCase("N")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            }
            case NA:
                Atom oxygen = (Atom) getAtomNode("O3\'");
                ArrayList<Bond> bonds = oxygen.getBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(oxygen);
                    if (other.getName().equalsIgnoreCase("P")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            default:
                return null;
        }
        return null;
        // Will generally indicate that you passed in a chain-terminal residue.
    }

    /**
     * Returns the Residue bonded to this Residue at this Residue's 5' or 
     * N-terminal end.  Any use of this method to add Residues to a sliding
     * window or similar MUST not add that residue if that residue has no 
     * Rotamers, as several algorithms (such as the distance matrix) assume
     * that all Residues being optimized have Rotamers.
     * 
     * @return The previous Residue.
     */
    public Residue getPreviousResidue() {
        switch (residueType) {
            case AA: {
                Atom nitrogen = (Atom) getAtomNode("N");
                ArrayList<Bond> bonds = nitrogen.getBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(nitrogen);
                    if (other.getName().equalsIgnoreCase("C")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            }
            case NA:
                Atom phosphate = (Atom) getAtomNode("P");
                if (phosphate == null) {
                    return null;
                }
                ArrayList<Bond> bonds = phosphate.getBonds();
                for (Bond b : bonds) {
                    Atom other = b.get1_2(phosphate);
                    if (other.getName().equalsIgnoreCase("O3\'")) {
                        return (Residue) other.getParent().getParent();
                    }
                }
                break;
            default:
                return null;
        }
        return null;
        // Will generally indicate that you passed in a chain-starting residue.
    }

    /**
     * {@inheritDoc}
     *
     * Allows adding Atoms to the Residue.
     */
    @Override
    public MSNode addMSNode(MSNode o) {
        Atom currentAtom = null;
        if (o instanceof Atom) {
            Atom newAtom = (Atom) o;
            Character newAlt = newAtom.getAltLoc();
            MSNode atoms = getAtomNode();
            currentAtom = (Atom) atoms.contains(newAtom);
            if (currentAtom == null) {
                currentAtom = newAtom;
                atoms.add(newAtom);
                setFinalized(false);
            } else {
                /**
                 * Allow overwriting of the root alternate conformer (' ' or
                 * 'A').
                 */
                Character currentAlt = currentAtom.getAltLoc();
                if (currentAlt.equals(' ') || currentAlt.equals('A')) {
                    if (!newAlt.equals(' ') && !newAlt.equals('A')) {
                        newAtom.xyzIndex = currentAtom.xyzIndex;
                        atoms.remove(currentAtom);
                        currentAtom = newAtom;
                        atoms.add(currentAtom);
                        setFinalized(false);
                    }
                }
            }
        } else {
            logger.warning("Only an Atom can be added to a Residue.");
        }
        return currentAtom;
    }

    /**
     * <p>
     * deleteAtom</p>
     *
     * @param atomToDelete a {@link ffx.potential.bonded.Atom} object.
     */
    public void deleteAtom(Atom atomToDelete) {
        MSNode atoms = getAtomNode();
        if (atoms.contains(atomToDelete) != null) {
            logger.info(" The following atom is being deleted from the model:\n"
                    + atomToDelete.toString());
            atoms.remove(atomToDelete);
        }
    }

    public ArrayList<Atom> getSideChainAtoms() {
        switch (residueType) {
            case NA:
                return null;
            case AA:
                ArrayList<Atom> atoms = getAtomList();
                ArrayList<Atom> ret = new ArrayList<Atom>(atoms);
                for (Atom atom : atoms) {
                    String name = atom.getName().toUpperCase();
                    if (name.equals("N") || name.equals("H") || name.equals("H1") || name.equals("H2") || name.equals("H3")
                            || name.equals("CA") || name.startsWith("HA")
                            || name.equals("C") || name.equals("O") || name.equals("OXT") || name.equals("OT2")) {
                        ret.remove(atom);
                    }
                }
                return ret;
            default:
                return null;
        }
    }

    private void assignResidueType() {
        String name = getName().toUpperCase();
        switch (residueType) {
            case AA:
                aa = null;
                try {
                    if (name.length() >= 2) {
                        aa = AA3.valueOf(name);
                    } else if (name.length() == 1) {
                        AA1 aa1 = AA1.valueOf(name);
                        aa = AA1toAA3.get(aa1);
                    }
                } catch (Exception e) {
                    aa = AA3.UNK;
                }
                break;
            case NA:
                na = null;
                try {
                    if (name.length() >= 2) {
                        na = NA3.valueOf(name);
                    } else if (name.length() == 1) {
                        NA1 na1 = NA1.valueOf(name);
                        na = NA1toNA3.get(na1);
                    }
                } catch (Exception e) {
                    na = NA3.UNK;
                }
                break;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void drawLabel(Canvas3D canvas, J3DGraphics2D g2d, Node node) {
        if (RendererCache.labelResidues) {
            double d[] = getCenter();
            point3d.x = d[0];
            point3d.y = d[1];
            point3d.z = d[2];
            RendererCache.getScreenCoordinate(canvas, node, point3d, point2d);
            g2d.drawString(getName(), (float) point2d.x, (float) point2d.y);
        }
        if (RendererCache.labelAtoms) {
            super.drawLabel(canvas, g2d, node);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Overidden equals method that return true if object is not equals to this,
     * is of the same class, has the same parent Polymer and the same sequence
     * number.
     */
    @Override
    public boolean equals(Object object) {
        if (this == object) {
            return true;
        } else if (object == null || getClass() != object.getClass()) {
            return false;
        }
        Residue other = (Residue) object;
        if (getParent() == null || other.getParent() == null) {
            return (getResidueNumber() == other.getResidueNumber());
        } else if (getParent() == other.getParent()) {
            return (getResidueNumber() == other.getResidueNumber());
        } else {
            return false;
        }
    }

    /**
     * {@inheritDoc}
     *
     * The Finalize method should be called once all atoms have been added to
     * the Residue. Geometry objects (Bonds, Angles, etc) are then formed,
     * followed by a determination of under-constrained (Dangeling) atoms.
     */
    @Override
    public void finalize(boolean finalizeGeometry) {
        setFinalized(false);
        getAtomNode().setName("Atoms (" + getAtomList().size() + ")");
        if (finalizeGeometry) {
            // constructValenceTerms();
            collectValenceTerms();
            removeLeaves();
        }
        // findDangelingAtoms();
        setCenter(getMultiScaleCenter(false));
        setFinalized(true);
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's default
     * O3' coordinates given a North pucker.
     *
     * @return a new double[] with default XYZ coordinates for O3' in a North
     * pucker.
     */
    public double[] getO3sNorth() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = O3sNorthCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's default
     * O3' coordinates given a South pucker.
     *
     * @return a new double[] with default XYZ coordinates for O3' in a South
     * pucker.
     */
    public double[] getO3sSouth() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = O3sSouthCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's original
     * C1' coordinates.
     *
     * @return a new double[] with original XYZ coordinates for C1'.
     */
    public double[] getC1sCoords() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = C1sCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's original
     * O4' coordinates.
     *
     * @return a new double[] with original XYZ coordinates for O4'.
     */
    public double[] getO4sCoords() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = O4sCoords[i];
        }
        return ret;
    }

    /**
     * Returns the position of this (presumably nucleic acid) Residue's original
     * C4' coordinates.
     *
     * @return a new double[] with original XYZ coordinates for C4'.
     */
    public double[] getC4sCoords() {
        double[] ret = new double[3];
        for (int i = 0; i < ret.length; i++) {
            ret[i] = C4sCoords[i];
        }
        return ret;
    }

    /**
     * Initializes this (presumably nucleic acid) Residue's C1s, O4s, C4s,
     * O3sNorth, and O3sSouth default coordinates based on default PDB atom
     * locations; to preserve rotamer independence, this must be called before
     * any NA rotamers are applied.
     */
    public void initializeDefaultAtomicCoordinates() {
        if (residueType != ResidueType.NA) {
            return;
        }
        boolean isDeoxy;
        try {
            switch (NucleicAcid3.valueOf(this.getName())) {
                case DAD:
                case DCY:
                case DGU:
                case DTY:
                    isDeoxy = true;
                    break;
                case CYT:
                case ADE:
                case THY:
                case URI:
                case GUA:
                default:
                    isDeoxy = false;
                    break;
            }
            C1sCoords = new double[3];
            ((Atom) getAtomNode("C1\'")).getXYZ(C1sCoords);
            O4sCoords = new double[3];
            ((Atom) getAtomNode("O4\'")).getXYZ(O4sCoords);
            C4sCoords = new double[3];
            ((Atom) getAtomNode("C4\'")).getXYZ(C4sCoords);

            /**
             * With the place flag set false, applySugarPucker returns
             * hypothetical O3' coordinates based on default atom positions and
             * the supplied sugar pucker.
             */
            O3sNorthCoords = RotamerLibrary.applySugarPucker(this, 1, isDeoxy, false);
            O3sSouthCoords = RotamerLibrary.applySugarPucker(this, 2, isDeoxy, false);
        } catch (Exception e) {
            logger.log(Level.WARNING, toString(), e);
        }
    }

    /**
     * Returns the position of this Residue's Alpha Carbon (if it is an amino
     * acid).
     *
     * @return a {@link javax.vecmath.Vector3d} object.
     */
    public Vector3d getAlpha3d() {
        Atom a = (Atom) getAtomNode("CA");
        if (a == null) {
            return null;
        }
        Vector3d temp = new Vector3d();
        a.getV3D(temp);
        return temp;
    }

    /**
     * Returns this Residues Parent Polymer name.
     *
     * @return a {@link java.lang.Character} object.
     */
    public Character getChainID() {
        return chainID;
    }

    // Public data access methods
    /**
     * Returns this Residue's sequence number.
     *
     * @return a int.
     */
    public int getResidueNumber() {
        return resNumber;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public final int hashCode() {
        int hash = hash(SEED, getParent().hashCode());
        hash = hash(hash, getResidueNumber());
        return hash(hash, getName());
    }

    /**
     * {@inheritDoc}
     *
     * Prints "Residue Number: x" to stdout.
     */
    @Override
    public void print() {
        logger.info(toString());
        super.print();
    }

    public double[] getSideChainCOM() {
        Vector3d v = new Vector3d();
        Vector3d v2 = new Vector3d();
        int count = 0;
        for (ListIterator li = getSideChainAtoms().listIterator(); li.hasNext();) {
            Atom a = (Atom) li.next();
            String id = a.getName();
            if (!id.equals("CA") && !id.equals("N") && !id.equals("C")
                    && !id.equals("O")) {
                a.getV3D(v2);
                v.add(v2);
                count++;
            } else if (id.equals("CA")) {
                a.print();
            }
        }
        v.scale(1.0 / count);
        double ret[] = new double[3];
        v.get(ret);
        return ret;
    }

    /**
     * <p>
     * printSideChainCOM</p>
     */
    public void logSideChainCOM() {
        double com[] = this.getSideChainCOM();
        logger.info(String.format(" %s %8.3f %8.3f %8.3f ", getName(), com[0], com[1], com[2]));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color,
            Material mat) {
// If Color by Residue, pass this Residue's Color
        if (newColorModel == RendererCache.ColorModel.RESIDUE) {
            switch (residueType) {
                case AA:
                    color = AA3Color.get(aa);
                    break;
                case NA:
                    color = NA3Color.get(na);
                    break;
                default:
                    color = null;
            }
            if (color == null) {
                return;
            }
            mat = RendererCache.materialFactory(color);
        } else if (newColorModel == RendererCache.ColorModel.STRUCTURE) {
            color = SSTypeColor.get(ssType);
            mat = RendererCache.materialFactory(color);
        }
        super.setColor(newColorModel, color, mat);
    }

    /**
     * <p>
     * setNumber</p>
     *
     * @param n a int.
     */
    public void setNumber(int n) {
        resNumber = n;
    }

    /**
     * <p>
     * Setter for the field <code>chainID</code>.</p>
     *
     * @param c a {@link java.lang.Character} object.
     */
    public void setChainID(Character c) {
        chainID = c;
    }

    /**
     * <p>
     * Getter for the field <code>segID</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getSegID() {
        return segID;
    }

    /**
     * <p>
     * setSSType</p>
     *
     * @param ss a {@link ffx.potential.bonded.Residue.SSType} object.
     */
    public void setSSType(SSType ss) {
        ssType = ss;
    }
    private String shortString = null;

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        if (shortString == null) {
            shortString = new String("" + resNumber + "-" + getName());
        }
        return shortString;
    }

    public enum AA {

        GLYCINE, ALANINE, VALINE, LEUCINE, ILLUECINE, SERINE, THREONINE,
        CYSTIENE, PROLINE, PHENYLALANINE, TYROSINE, TYPTOPHAN, ASPARTATE,
        ASPARTAMINE, GLUTAMATE, GLUTAMINE, METHIONINE, LYSINE, ARGININE,
        HISTIDINE;
    }

    public enum AA1 {

        G, A, V, L, I, S, T, C, P, F, Y, W, D, N, E, Q, M, K, R, H, U, Z, O, B,
        J, f, a, n, m, X;
    }

    public enum AA3 {

        GLY, ALA, VAL, LEU, ILE, SER, THR, CYS, PRO, PHE, TYR, TRP, ASP, ASN,
        GLU, GLN, MET, LYS, ARG, HIS, HID, HIE, ORN, AIB, PCA, FOR, ACE, NH2,
        NME, UNK;
    }

    public enum NA {

        ADENINE, CYTOSINE, GUANINE, URACIL, DEOXYADENINE, DEOXYCYTOSINE,
        DEOXYGUANINE, THYMINE, MONOPHOSPHATE, DIPHOSPHATE, TRIPHOSPHATE;
    }

    public enum NA1 {

        A, C, G, U, D, I, B, T, P, Q, R, X;
    }

    public enum NA3 {

        A, C, G, U, DA, DC, DG, DT, MPO, DPO, TPO, UNK;
    }

    public enum ResidueType {

        NA, AA, UNK;
    }

    public enum SSType {

        NONE, HELIX, SHEET, TURN;
    }
}
