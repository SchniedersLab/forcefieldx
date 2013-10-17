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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import static java.lang.String.format;

import javax.media.j3d.*;
import javax.vecmath.Color3f;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.VDWType;

import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

/**
 * The Atom class represents a single atom and defines its alternate
 * conformations and molecular mechanics atom type.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class Atom extends MSNode implements Comparable<Atom> {

    private static final Logger logger = Logger.getLogger(Atom.class.getName());
    private static Point3d point3d = new Point3d();
    private static Point2d point2d = new Point2d();
    /**
     * Constant <code>AtomColor</code>
     */
    public static final Map<Integer, Color3f> AtomColor;
    /**
     * Constant <code>AtomVDW</code>
     */
    public static final Map<Integer, Double> AtomVDW;
    /**
     * Constant <code>SP3=4</code>
     */
    public static final int SP = 2, SP2 = 3, SP3 = 4;
    /**
     * Constant <code>hybridTable</code>
     */
    public static final Map<String, Integer> hybridTable;

    static {
        // Initialize HashMaps
        AtomColor = new HashMap<Integer, Color3f>();
        AtomVDW = new HashMap<Integer, Double>();
        hybridTable = new HashMap<String, Integer>();
        // van der Waals
        AtomVDW.put(0, 1.0);
        AtomVDW.put(1, 1.20);
        AtomVDW.put(2, 1.22);
        AtomVDW.put(3, 0.78);
        AtomVDW.put(4, 0.34);
        AtomVDW.put(5, 2.08);
        AtomVDW.put(6, 1.85);
        AtomVDW.put(7, 1.54);
        AtomVDW.put(8, 1.40);
        AtomVDW.put(9, 1.35);
        AtomVDW.put(10, 1.60);
        AtomVDW.put(11, 0.98);
        AtomVDW.put(12, 0.78);
        AtomVDW.put(13, 0.57);
        AtomVDW.put(14, 2.00);
        AtomVDW.put(15, 1.90);
        AtomVDW.put(16, 1.85);
        AtomVDW.put(17, 1.81);
        AtomVDW.put(18, 1.91);
        AtomVDW.put(19, 1.33);
        AtomVDW.put(20, 1.06);
        AtomVDW.put(21, 0.91);
        AtomVDW.put(22, 0.83);
        AtomVDW.put(23, 0.82);
        AtomVDW.put(24, 2.00);
        AtomVDW.put(25, 2.00);
        AtomVDW.put(26, 2.00);
        AtomVDW.put(27, 2.00);
        AtomVDW.put(28, 2.00);
        AtomVDW.put(29, 2.00);
        AtomVDW.put(30, 2.00);
        AtomVDW.put(31, 2.00);
        AtomVDW.put(32, 2.00);
        AtomVDW.put(33, 2.00);
        AtomVDW.put(34, 2.00);
        AtomVDW.put(35, 1.95);
        AtomVDW.put(36, 1.89);
        for (int i = 37; i < 109; i++) {
            AtomVDW.put(i, 2.00);
        }
        // Colors
        AtomColor.put(0, RendererCache.RED);
        AtomColor.put(1, RendererCache.WHITE);
        AtomColor.put(2, RendererCache.GREEN);
        AtomColor.put(3, RendererCache.MAGENTA);
        AtomColor.put(4, RendererCache.MAGENTA);
        AtomColor.put(5, RendererCache.MAGENTA);
        AtomColor.put(6, RendererCache.GRAY);
        AtomColor.put(7, RendererCache.BLUE);
        AtomColor.put(8, RendererCache.RED);
        AtomColor.put(9, RendererCache.GREEN);
        AtomColor.put(10, RendererCache.GREEN);
        AtomColor.put(11, RendererCache.MAGENTA);
        AtomColor.put(12, RendererCache.MAGENTA);
        AtomColor.put(13, RendererCache.MAGENTA);
        AtomColor.put(14, RendererCache.GRAY);
        AtomColor.put(15, RendererCache.ORANGE);
        AtomColor.put(16, RendererCache.YELLOW);
        AtomColor.put(17, RendererCache.GREEN);
        AtomColor.put(18, RendererCache.GREEN);
        for (int i = 19; i < 109; i++) {
            if (i != 36 && i != 54 && i != 86) {
                AtomColor.put(i, RendererCache.MAGENTA);
            } else {
                AtomColor.put(i, RendererCache.GREEN);
            }
        }

        // Default hybridization
        hybridTable.put("1", 1);
        hybridTable.put("6", 4);
        hybridTable.put("7", 3);
        hybridTable.put("8", 2);
        hybridTable.put("15", 4);
        hybridTable.put("16", 2);
        hybridTable.put("19", 0);
        hybridTable.put("26", 8);
    }

    private boolean hetatm;
    /**
     * Either the PDB "name" record or the molecular mechanics atom type name.
     *
     * @since 1.0
     */
    //private String name = null;
    /**
     * Contiguous atom index ranging from 0..nAtoms.
     *
     * @since 1.0
     */
    public int xyzIndex = -1;
    /**
     * PDB "resname" record.
     *
     * @since 1.0
     */
    private String resName = null;
    /**
     * PDB "resSeq" record.
     *
     * @since 1.0
     */
    private int resSeq = -1;
    /**
     * PDB "chainID" or "segID" record.
     *
     * @since 1.0
     */
    private Character chainID = null;
    /**
     * Array of altLoc identifiers defined for this atom.
     *
     * @since 1.0
     */
    private Character altLoc;
    /**
     * Array of XYZ coordinates for each altLoc.
     *
     * @since 1.0
     */
    private double xyz[];
    /**
     * Array of XYZ gradients for each altLoc.
     *
     * @since 1.0
     */
    private double xyzGradient[];
    /**
     * Array of occupancy values for each altLoc.
     *
     * @since 1.0
     */
    private double occupancy;
    /**
     * Array of occupancy gradients for each altLoc.
     *
     * @since 1.0
     */
    private double occupancyGradient;
    /**
     * Array of tempFactor values for each altLoc.
     *
     * @since 1.0
     */
    private double tempFactor;
    /**
     * Array fo tempFactorGradients
     *
     * @since 1.0
     */
    private double tempFactorGradient;
    /**
     * List of anisou tensors for each altLoc.
     *
     * @since 1.0
     */
    private double anisou[];
    /**
     * List of anisou tensors for each altLoc.
     *
     * @since 1.0
     */
    private double anisouGradient[];
    /**
     * If the atom is active, it should be included in target functions.
     */
    private boolean active = true;
    private String segID = null;
    private double formFactorWidth = 3.5;
    private int formFactorIndex = -1;
    private ArrayList<Vector3d> trajectory;
    // Molecular Mechanics Info
    private AtomType atomType = null;
    private MultipoleType multipoleType = null;
    private PolarizeType polarizeType = null;
    private VDWType vdwType = null;
    private Atom[] multipoleReferenceSites = null;
    private double globalDipole[] = null;
    private double globalQuadrupole[][] = null;
    private boolean applyState = false;
    // solvation
    private double bornRadius;
    // Connectivity information.
    private final ArrayList<Bond> bonds = new ArrayList<Bond>();
    private final ArrayList<Angle> angles = new ArrayList<Angle>();
    private final ArrayList<Torsion> torsions = new ArrayList<Torsion>();
    private final ArrayList<Atom> one_5s = new ArrayList<Atom>();
    /**
     * *************************************************************************
     */
    // Java3D methods and variables for visualization of this Atom.
    // The current ViewModel
    private ViewModel viewModel = ViewModel.INVISIBLE;
    private ViewModel polygonType = ViewModel.FILL;
    private ColorModel colorModel = ColorModel.CPK;
    // Java3D Scenegraph Objects
    private Shape3D sphere, cylinder, cone;
    private BranchGroup branchGroup, vectorBranchGroup;
    private TransformGroup transformGroup;
    private Transform3D transform3D;
    // Appearance and Coloring
    private Appearance appearance;
    private Color3f currentCol, previousCol;
    private Color3f userColor = RendererCache.userColor;
    private int detail = RendererCache.detail;
    private double radius = RendererCache.radius;
    private double scale = 1.0;
    // "stale" is True if this Atom's J3D transforms need to be updated before
    // making it visible
    private boolean stale = false;
    private String shortString = null;
    private Vector3d vector3d = new Vector3d();

    /**
     * Default constructor.
     *
     * @param name The Atom's PDB name.
     * @since 1.0
     */
    public Atom(String name) {
        super(name, 1);
        //this.name = name;
        currentCol = previousCol = RendererCache.toAtomColor(name);
        colorModel = ColorModel.CPK;
    }

    /**
     * Constructor used when parsing XYZ files.
     *
     * @param xyzIndex Contiguous, unique atom index between 1..nAtoms.
     * @param name The Atom's molecular mechanics name.
     * @param atomType Molecular mechanics atom type.
     * @param xyz Cartesian coordinates.
     * @since 1.0
     */
    public Atom(int xyzIndex, String name, AtomType atomType, double[] xyz) {
        //this(name, atomType, xyz);
        this(name);
        this.xyzIndex = xyzIndex;
        this.atomType = atomType;
        this.xyz = xyz;
        this.xyzGradient = new double[3];
        setAllowsChildren(false);
        if (atomType != null) {
            currentCol = previousCol = AtomColor.get(atomType.atomicNumber);
        }
    }

    /**
     * Constructor used when parsing PDB files.
     *
     * @param xyzIndex Contiguous, unique atom index between 1..nAtoms.
     * @param name The Atom's molecular mechanics name.
     * @param altLoc The alternate locations (' ' or 'A' or 'B' etc.).
     * @param xyz Cartesian coordinates.
     * @param resName Residue name.
     * @param resSeq Residue sequence number.
     * @param chainID Possible redundant chain ID.
     * @param occupancy Crystallographic occupancy.
     * @param tempFactor Crystallographic B-factor.
     * @param segID Unique segment ID.
     * @since 1.0.
     */
    public Atom(int xyzIndex, String name,
            Character altLoc, double[] xyz, String resName, int resSeq,
            Character chainID, double occupancy, double tempFactor,
            String segID) {
        this(xyzIndex, name, null, xyz);
        this.resName = resName;
        this.resSeq = resSeq;
        this.chainID = chainID;
        this.altLoc = altLoc;
        this.occupancy = occupancy;
        this.tempFactor = tempFactor;
        this.segID = segID;
    }

    /**
     * {@inheritDoc}
     *
     * Implementation of the Comparable interface.
     */
    @Override
    public int compareTo(Atom a) {
        if (a == null) {
            throw new NullPointerException();
        }
        if (a == this) {
            return 0;
        }
        int a0 = a.xyzIndex;
        if (xyzIndex < a0) {
            return -1;
        }
        if (xyzIndex > a0) {
            return 1;
        }
        // There should not be duplicate, identical atom objects.
        // assert (xyzIndex != a0);
        return 0;
    }

    /**
     * <p>
     * addTrajectoryCoords</p>
     *
     * @param coords a {@link javax.vecmath.Vector3d} object.
     * @param position a int.
     */
    public void addTrajectoryCoords(Vector3d coords, int position) {
        if (trajectory == null) {
            trajectory = new ArrayList<Vector3d>();
            trajectory.add(0, new Vector3d(xyz));
        }
        trajectory.add(position, coords);
    }

    /**
     * <p>
     * setHetero</p>
     *
     * @param hetatm a boolean.
     */
    public void setHetero(boolean hetatm) {
        this.hetatm = hetatm;
    }

    /**
     * <p>
     * isHetero</p>
     *
     * @return a boolean.
     */
    public boolean isHetero() {
        return hetatm;
    }

    public Atom copy() {
        double coords[] = {xyz[0], xyz[1], xyz[2]};
        Atom atom = new Atom(getXYZIndex(), getName(), getAltLoc(), coords,
                getResidueName(), getResidueNumber(), getChainID(),
                getOccupancy(), getTempFactor(), getSegID());
        atom.setAtomType(getAtomType());
        return atom;
    }

    /**
     * <p>
     * isHydrogen</p>
     *
     * @return a boolean.
     */
    public boolean isHydrogen() {
        if (atomType == null) {
            return false;
        }
        if (atomType.atomicNumber == 1) {
            return true;
        }
        return false;
    }

    public boolean isActive() {
        return active;
    }

    public void setActive(boolean active) {
        this.active = active;
    }

    /**
     * <p>
     * isDeuterium</p>
     *
     * @return a boolean.
     */
    public boolean isDeuterium() {
        String name = getName();
        return (isHydrogen() && (name.charAt(0) == 'D'
                || name.charAt(0) == 'd'));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void drawLabel(Canvas3D canvas, J3DGraphics2D g2d, Node node) {
        if (RendererCache.labelAtoms) {
            point3d.x = getX();
            point3d.y = getY();
            point3d.z = getZ();
            RendererCache.getScreenCoordinate(canvas, node, point3d, point2d);
            g2d.drawString(toShortString(), (float) point2d.x,
                    (float) point2d.y);
        }

    }

    /**
     * {@inheritDoc}
     *
     * Overidden equals method.
     */
    @Override
    public final boolean equals(Object object) {
        if (this == object) {
            return true;
        }
        if (object == null || !(object instanceof Atom)) {
            return false;
        }
        Atom other = (Atom) object;

        return (other.resName != null && other.resName.equals(resName)
                && other.resSeq == resSeq
                && other.getName() != null && other.getName().equals(getName())
                && other.segID != null && other.segID.equals(segID));
    }

    /**
     * <p>
     * Getter for the field <code>angles</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Angle> getAngles() {
        return angles;
    }

    public Angle getAngle(Atom centralAtom, Atom endAtom) {
        for (Angle angle : angles) {
            Atom atom13 = angle.get1_3(this);
            if (atom13 != null && atom13.equals(endAtom)
                    && angle.getCentralAtom().equals(centralAtom)) {
                return angle;
            }
        }
        return null;
    }

    /**
     * <p>
     * get1_5s</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Atom> get1_5s() {
        return one_5s;
    }

    /**
     * <p>
     * getAtomAppearance</p>
     *
     * @return a {@link javax.media.j3d.Appearance} object.
     */
    public Appearance getAtomAppearance() {
        if (appearance == null) {
            appearance = RendererCache.appearanceFactory(currentCol,
                    polygonType);
        }
        return appearance;
    }

    /**
     * <p>
     * getAtomColor</p>
     *
     * @return a {@link javax.vecmath.Color3f} object.
     */
    public Color3f getAtomColor() {
        return currentCol;
    }

    /**
     * Gets the Atomic Number
     *
     * @return Atomic Number
     */
    public int getAtomicNumber() {
        return atomType.atomicNumber;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public ArrayList<Atom> getAtomList() {
        ArrayList<Atom> atoms = new ArrayList<Atom>();
        atoms.add(this);
        return atoms;
    }

    /**
     * <p>
     * Getter for the field <code>atomType</code>.</p>
     *
     * @return a {@link ffx.potential.parameters.AtomType} object.
     */
    public AtomType getAtomType() {
        return atomType;
    }

    /**
     * <p>
     * getBond</p>
     *
     * @param a a {@link ffx.potential.bonded.Atom} object.
     * @return a {@link ffx.potential.bonded.Bond} object.
     */
    public Bond getBond(Atom a) {
        if (bonds == null) {
            return null;
        }
        for (Bond bond : bonds) {
            if (bond.get1_2(a) == this) {
                return bond;
            }
        }
        return null;
    }

    /**
     * Gets the list of the Bonds <b>this</b> Atom helps to form
     *
     * @return A list of the bonds this atom helps to form
     */
    public ArrayList<Bond> getBonds() {
        return bonds;
    }

    /**
     * <p>
     * Getter for the field <code>bornRadius</code>.</p>
     *
     * @return a double.
     */
    public double getBornRadius() {
        return bornRadius;
    }

    /**
     * <p>
     * getBornVolume</p>
     *
     * @return a double.
     */
    public double getBornVolume() {
        return 2.0;
    }

    /**
     * Get the chain name
     *
     * @return String
     */
    public Character getChainID() {
        if (chainID != null) {
            return chainID;
        }
        Polymer p = (Polymer) this.getMSNode(Polymer.class);
        if (p == null) {
            return null;
        }
        chainID = p.getName().charAt(0);
        return chainID;
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
     * Gets the partial atomic charge
     *
     * @return partial atomic charge
     */
    public double getCharge() {
        return 1.0;
    }

    /**
     * Finds a Torsion which contains this atom, and atoms 2, 3, and 4.
     *
     * @param atom2
     * @param atom3
     * @param atom4
     * @return Torsion.
     */
    public Torsion getTorsion(Atom atom2, Atom atom3, Atom atom4) {
        for (Torsion torsion : torsions) {
            if (torsion.compare(this, atom2, atom3, atom4)) {
                return torsion;
            }
        }
        return null;
    }

    /**
     * <p>
     * Getter for the field <code>torsions</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Torsion> getTorsions() {
        return torsions;
    }

    /**
     * Gets the Epsilon value
     *
     * @return Epsilon value
     */
    public double getEpsilon() {
        return 1.0;
    }

    /**
     * <p>
     * applyLambda</p>
     *
     * @return a boolean.
     */
    public boolean applyLambda() {
        return applyState;
    }

    /**
     * <p>
     * setApplyLambda</p>
     *
     * @param applyState a boolean.
     */
    public void setApplyLambda(boolean applyState) {
        this.applyState = applyState;
    }

    /**
     * <p>
     * getForce</p>
     *
     * @param t an array of double.
     */
    public void getForce(double[] t) {
        if (xyzGradient == null || t == null) {
            return;
        }
        t[0] = xyzGradient[0];
        t[1] = xyzGradient[1];
        t[2] = xyzGradient[2];
    }

    /**
     * Gets the Atomic Hybridization
     *
     * @return Atomic Hybridization
     */
    public int getHybridization() {
        return atomType.valence;
    }

    /**
     * Gets the atom ID
     *
     * @return atom ID
     */
    public String getIdent() {
        return atomType.environment;
    }

    /**
     * Gets the atom Key
     *
     * @return atom Key
     */
    public String getKey() {
        if (atomType != null) {
            return atomType.getKey();
        }
        return null;
    }

    /**
     * Gets the Atomic Mass.
     *
     * @return Atomic Mass
     */
    public double getMass() {
        return atomType.atomicWeight;
    }

    /**
     * <p>
     * Getter for the field <code>multipoleReferenceSites</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public Atom[] getMultipoleReferenceSites() {
        return multipoleReferenceSites;
    }

    /**
     * <p>
     * Getter for the field <code>multipoleType</code>.</p>
     *
     * @return a {@link ffx.potential.parameters.MultipoleType} object.
     */
    public MultipoleType getMultipoleType() {
        return multipoleType;
    }

    /**
     * <p>
     * Getter for the field <code>polarizeType</code>.</p>
     *
     * @return a {@link ffx.potential.parameters.PolarizeType} object.
     */
    public PolarizeType getPolarizeType() {
        return polarizeType;
    }

    /**
     * <p>
     * getNumAngles</p>
     *
     * @return a int.
     */
    public final int getNumAngles() {
        if (angles == null) {
            return 0;
        }
        return angles.size();
    }

    /**
     * Gets the number of atoms bonded to <b>this</b> Atom
     *
     * @return Number of atoms bonded to this atom
     */
    public final int getNumBonds() {
        if (bonds == null) {
            return 0;
        }
        return bonds.size();
    }

    /**
     * <p>
     * getNumDihedrals</p>
     *
     * @return a int.
     */
    public final int getNumDihedrals() {
        if (torsions == null) {
            return 0;
        }
        return torsions.size();
    }

    /**
     * <p>
     * getRDielectric</p>
     *
     * @return a double.
     */
    public double getRDielectric() {
        return 2.0;
    }

    /**
     * Get the residue name
     *
     * @return String
     */
    public String getResidueName() {
        if (resName != null) {
            return resName;
        }
        Residue r = (Residue) getMSNode(Residue.class);
        return r.getName();
    }

    /**
     * <p>
     * getResidueNumber</p>
     *
     * @return a int.
     */
    public int getResidueNumber() {
        return resSeq;
    }

    /**
     * <p>
     * Setter for the field <code>resName</code>.</p>
     *
     * @param resName a {@link java.lang.String} object.
     */
    public void setResName(String resName) {
        this.resName = resName;
    }

    /**
     * Gets the Sigma value
     *
     * @return Sigma value
     */
    public double getSigma() {
        return 1.0;
    }

    /**
     * <p>
     * getTrajectoryCoords</p>
     *
     * @param position a int.
     * @return a {@link javax.vecmath.Vector3d} object.
     */
    public Vector3d getTrajectoryCoords(int position) {
        return trajectory.get(position);
    }

    /**
     * <p>
     * getTrajectoryLength</p>
     *
     * @return a int.
     */
    public int getTrajectoryLength() {
        return trajectory.size();
    }

    /**
     * <p>
     * getType</p>
     *
     * @return a int.
     */
    public int getType() {
        return atomType.type;
    }

    /**
     * Gets the Atom's Cartesian Coordinates return The Cartesian Coordinates
     *
     * @param temp a {@link javax.vecmath.Vector3d} object.
     */
    public void getV3D(Vector3d temp) {
        temp.set(xyz);
    }

    /**
     * Gets the van der Waals radius.
     *
     * @return a double.
     */
    public double getVDWR() {
        return 1.0;
    }

    /**
     * <p>
     * Setter for the field <code>altLoc</code>.</p>
     *
     * @param a a {@link java.lang.Character} object.
     */
    public void setAltLoc(Character a) {
        altLoc = a;
    }

    /**
     * <p>
     * Getter for the field <code>altLoc</code>.</p>
     *
     * @return a {@link java.lang.Character} object.
     */
    public Character getAltLoc() {
        return altLoc;
    }

    /**
     * Gets the x coordinate
     *
     * @return x coordinate
     */
    public double getX() {
        return xyz[0];
    }

    /**
     * <p>
     * getXYZ</p>
     *
     * @param x an array of double.
     */
    public void getXYZ(double[] x) {
        x[0] = xyz[0];
        x[1] = xyz[1];
        x[2] = xyz[2];
    }

    /**
     * <p>
     * getXYZ</p>
     *
     * @return an array of double.
     */
    public double[] getXYZ() {
        return xyz;
    }

    /**
     * Gets the XYZ Index
     *
     * @return XYZ Index
     */
    public final int getXYZIndex() {
        return xyzIndex;
    }

    /**
     * Gets the y coordinate
     *
     * @return y coordinate
     */
    public final double getY() {
        return xyz[1];
    }

    /**
     * Gets the z coordinate
     *
     * @return z coordinate
     */
    public final double getZ() {
        return xyz[2];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public final int hashCode() {
        return hash(SEED, xyzIndex);
    }

    /**
     * Create the Sphere Java3D objects.
     *
     * @param newShapes List
     */
    private void initSphere(List<BranchGroup> newShapes) {
        if (appearance == null) {
            appearance = RendererCache.appearanceFactory(currentCol,
                    ViewModel.FILL);
        }
        if (transform3D == null) {
            transform3D = RendererCache.transform3DFactory(new Vector3d(xyz),
                    scale);
        } else {
            transform3D.setTranslation(new Vector3d(xyz));
            transform3D.setScale(scale);
        }
        detail = RendererCache.detail;
        branchGroup = RendererCache.sphereFactory(appearance, detail,
                transform3D);
        transformGroup = (TransformGroup) branchGroup.getChild(0);
        sphere = (Shape3D) transformGroup.getChild(0);
        sphere.setUserData(this);
        newShapes.add(branchGroup);
    }

    /**
     * Checks to see if an Atom is bonded to <b>this</b> Atom
     *
     * @param a Atom to check
     * @return True is Atom a is bonded to <b>this</b>this atom
     */
    public boolean isBonded(Atom a) {
        for (Bond bond : bonds) {
            if (bond.get1_2(a) == this) {
                return true;
            }
        }
        return false;
    }

    /**
     * <p>
     * is_1_3</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @return a boolean.
     */
    public boolean is_1_3(Atom atom) {
        for (Angle a : angles) {
            if (a.get1_3(atom) == this) {
                return true;
            }
        }
        return false;
    }

    /**
     * <p>
     * is_12_or_13</p>
     *
     * @param a a {@link ffx.potential.bonded.Atom} object.
     * @return a boolean.
     */
    public boolean is_12_or_13(Atom a) {
        if (isBonded(a) || is_1_3(a)) {
            return true;
        }
        return false;
    }

    /**
     * Gets whether or not the Atom is under-constrained
     *
     * @return True if the atom is under-constrained (ie has can accept bonds)
     */
    public boolean isDangeling() {
        Integer hybrid = hybridTable.get("" + atomType.atomicNumber);
        if (hybrid == null) {
            return false;
        }
        int result = hybrid.compareTo(bonds.size());
        return (result > 0);
    }

    /**
     * <p>
     * isTrigonal</p>
     *
     * @return a boolean.
     */
    public boolean isTrigonal() {
        if (bonds != null && bonds.size() == 3) {
            return true;
        }
        return false;
    }

    /**
     * <p>
     * isStale</p>
     *
     * @return a boolean.
     */
    public boolean isStale() {
        return stale;
    }

    // True if this Atom's Sphere or Vector is visible
    /**
     * <p>
     * isVisible</p>
     *
     * @return a boolean.
     */
    public boolean isVisible() {
        return (viewModel != ViewModel.INVISIBLE);
    }

    /**
     * Add a vector to the Atom's current position vector
     *
     * @param d Vector to add to the current position
     */
    public void move(double[] d) {
        xyz[0] += d[0];
        xyz[1] += d[1];
        xyz[2] += d[2];
        stale = true;
    }

    /**
     * <p>
     * moveTo</p>
     *
     * @param a a double.
     * @param b a double.
     * @param c a double.
     */
    public void moveTo(double a, double b, double c) {
        xyz[0] = a;
        xyz[1] = b;
        xyz[2] = c;
        stale = true;
    }

    /**
     * <p>
     * setXYZ</p>
     *
     * @param xyz an array of double.
     */
    public void setXYZ(double xyz[]) {
        this.xyz = xyz;
    }

    /**
     * Moves the atom to the specified location
     *
     * @param d Location to move <b>this</b> Atom to
     */
    public void moveTo(double[] d) {
        moveTo(d[0], d[1], d[2]);
    }

    /**
     * <p>
     * moveTo</p>
     *
     * @param v a {@link javax.vecmath.Vector3d} object.
     */
    public void moveTo(Vector3d v) {
        moveTo(v.x, v.y, v.z);
    }

    /**
     * {@inheritDoc}
     *
     * Prints the atom identity and Cartesian coordinates to the logger.
     */
    @Override
    public final void print() {
        logger.info(toString());
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void removeFromParent() {
        super.removeFromParent();
        if (angles != null) {
            angles.clear();
        }
        if (torsions != null) {
            torsions.clear();
        }
        if (trajectory != null) {
            trajectory.clear();
        }
        trajectory = null;
        stale = false;
        viewModel = ViewModel.INVISIBLE;
        RendererCache.poolTransform3D(transform3D);
        if (branchGroup != null) {
            branchGroup.detach();
            branchGroup.setUserData(null);
            RendererCache.poolSphere(branchGroup);
            branchGroup = null;
        }
    }

    /**
     * <p>
     * Setter for the field <code>atomType</code>.</p>
     *
     * @param atomType a {@link ffx.potential.parameters.AtomType} object.
     */
    public void setAtomType(AtomType atomType) {
        this.atomType = atomType;
    }

    /**
     * <p>
     * setVDWType</p>
     *
     * @param vdwType a {@link ffx.potential.parameters.VDWType} object.
     */
    public void setVDWType(VDWType vdwType) {
        this.vdwType = vdwType;
    }

    /**
     * <p>
     * getVDWType</p>
     *
     * @return a {@link ffx.potential.parameters.VDWType} object.
     */
    public VDWType getVDWType() {
        return vdwType;
    }

    /**
     * <p>
     * Setter for the field <code>tempFactor</code>.</p>
     *
     * @param tempFactor a double.
     */
    public void setTempFactor(double tempFactor) {
        this.tempFactor = tempFactor;
    }

    /**
     * <p>
     * Getter for the field <code>tempFactor</code>.</p>
     *
     * @return a double.
     */
    public double getTempFactor() {
        return tempFactor;
    }

    /**
     * <p>
     * Setter for the field <code>tempFactorGradient</code>.</p>
     *
     * @param tempFactorGradient a double.
     */
    public void setTempFactorGradient(double tempFactorGradient) {
        this.tempFactorGradient = tempFactorGradient;
    }

    /**
     * <p>
     * addToTempFactorGradient</p>
     *
     * @param tempFactorGradient a double.
     */
    public void addToTempFactorGradient(double tempFactorGradient) {
        this.tempFactorGradient += tempFactorGradient;
    }

    /**
     * <p>
     * Getter for the field <code>tempFactorGradient</code>.</p>
     *
     * @return a double.
     */
    public double getTempFactorGradient() {
        return tempFactorGradient;
    }

    /**
     * <p>
     * Setter for the field <code>occupancy</code>.</p>
     *
     * @param occupancy a double.
     */
    public void setOccupancy(double occupancy) {
        this.occupancy = occupancy;
    }

    /**
     * <p>
     * Getter for the field <code>occupancy</code>.</p>
     *
     * @return a double.
     */
    public double getOccupancy() {
        return occupancy;
    }

    /**
     * <p>
     * Setter for the field <code>occupancyGradient</code>.</p>
     *
     * @param occupancyGradient a double.
     */
    public void setOccupancyGradient(double occupancyGradient) {
        this.occupancyGradient = occupancyGradient;
    }

    /**
     * <p>
     * addToOccupancyGradient</p>
     *
     * @param occupancyGradient a double.
     */
    public void addToOccupancyGradient(double occupancyGradient) {
        this.occupancyGradient += occupancyGradient;
    }

    /**
     * <p>
     * Getter for the field <code>occupancyGradient</code>.</p>
     *
     * @return a double.
     */
    public double getOccupancyGradient() {
        return occupancyGradient;
    }

    /**
     * <p>
     * Setter for the field <code>anisou</code>.</p>
     *
     * @param anisou an array of double.
     */
    public void setAnisou(double[] anisou) {
        this.anisou = anisou;
    }

    /**
     * <p>
     * Getter for the field <code>anisou</code>.</p>
     *
     * @return an array of double.
     */
    public double[] getAnisou() {
        return anisou;

    }

    /**
     * <p>
     * Setter for the field <code>anisouGradient</code>.</p>
     *
     * @param anisou an array of double.
     */
    public void setAnisouGradient(double[] anisou) {
        this.anisouGradient = anisou;
    }

    /**
     * <p>
     * addToAnisouGradient</p>
     *
     * @param anisouGradient an array of double.
     */
    public void addToAnisouGradient(double[] anisouGradient) {
        if (this.anisouGradient == null) {
            this.anisouGradient = new double[6];
        }
        for (int i = 0; i < 6; i++) {
            this.anisouGradient[i] += anisouGradient[i];
        }
    }

    /**
     * <p>
     * Getter for the field <code>anisouGradient</code>.</p>
     *
     * @return an array of double.
     */
    public double[] getAnisouGradient() {
        return anisouGradient;
    }

    /**
     * <p>
     * Setter for the field <code>formFactorWidth</code>.</p>
     *
     * @param width a double.
     */
    public void setFormFactorWidth(double width) {
        formFactorWidth = width;
    }

    /**
     * <p>
     * Getter for the field <code>formFactorWidth</code>.</p>
     *
     * @return a double.
     */
    public double getFormFactorWidth() {
        return formFactorWidth;
    }

    /**
     * <p>
     * Setter for the field <code>formFactorIndex</code>.</p>
     *
     * @param formFactorIndex a int.
     */
    public void setFormFactorIndex(int formFactorIndex) {
        this.formFactorIndex = formFactorIndex;
    }

    /**
     * <p>
     * Getter for the field <code>formFactorIndex</code>.</p>
     *
     * @return a int.
     */
    public int getFormFactorIndex() {
        return formFactorIndex;
    }

    /**
     * <p>
     * setAngle</p>
     *
     * @param a a {@link ffx.potential.bonded.Angle} object.
     */
    public void setAngle(Angle a) {
        angles.add(a);
    }

    /**
     * Specify that <b>this</b> Atom is part of a Bond
     *
     * @param b Bond that <b>this</b> Atom is part of
     */
    public void setBond(Bond b) {
        bonds.add(b);
    }

    /**
     * <p>
     * removeBond</p>
     *
     * @param b a {@link ffx.potential.bonded.Bond} object.
     */
    public void removeBond(Bond b) {
        bonds.remove(b);
    }

    /**
     * <p>
     * set1_5</p>
     *
     * @param a a {@link ffx.potential.bonded.Atom} object.
     */
    public void set1_5(Atom a) {
        one_5s.add(a);
    }

    /**
     * Set the effective Born Radius.
     *
     * @param bornRadius a double.
     */
    public void setBornRadius(double bornRadius) {
        this.bornRadius = bornRadius;
    }

    /**
     * {@inheritDoc}
     *
     * Polymorphic setColor method.
     */
    @Override
    public void setColor(ColorModel newColorModel, Color3f newCol, Material newMat) {
        switch (newColorModel) {
            case CPK:
                newCol = RendererCache.getColor(this, ColorModel.CPK);
                if (newCol == currentCol) {
                    return;
                }
                colorModel = newColorModel;
                currentCol = previousCol = newCol;
                break;
            case USERCOLOR:
                colorModel = newColorModel;
                currentCol = previousCol = userColor;
                break;
            case APPLYUSERCOLOR:
                userColor = RendererCache.userColor;
                currentCol = previousCol = userColor;
                break;
            case MONOCHROME:
                colorModel = newColorModel;
                currentCol = previousCol = RendererCache.WHITE;
                break;
            case SELECT:
                if (isSelected()) {
                    newCol = RendererCache.selectionColor;
                    if (newCol != currentCol) {
                        currentCol = newCol;
                    } else {
                        return;
                    }
                } else {
                    currentCol = previousCol;
                }
                break;
            case PICK:
                newCol = RendererCache.pickingColor;
                if (newCol != currentCol) {
                    currentCol = newCol;
                } else {
                    return;
                }
                break;
            case REVERT:
                if (RendererCache.highlightSelections && isSelected()) {
                    currentCol = RendererCache.selectionColor;
                } else {
                    currentCol = previousCol;
                }
                break;
            case PARTIALCHARGE:
                newCol = RendererCache.getColor(this, ColorModel.PARTIALCHARGE);
                if (newCol == currentCol) {
                    return;
                }
                colorModel = newColorModel;
                currentCol = previousCol = newCol;
                break;
            default:
                // Check for a Color Choice sent from a higher level structure
                // (residue,polymer,etc)
                if (newCol == currentCol || newCol == null) {
                    return;
                }
                colorModel = newColorModel;
                currentCol = previousCol = newCol;
        }
// Apply the Color Change
        appearance = RendererCache.appearanceFactory(currentCol, polygonType);
        if (branchGroup != null && viewModel != ViewModel.INVISIBLE) {
            sphere.setAppearance(appearance);
        }
        if (bonds != null) {
            for (Bond bond : bonds) {
                bond.setColor(this);
            }
        }
    }

    /**
     * <p>
     * setCurrentCycle</p>
     *
     * @param cycle a int.
     */
    public void setCurrentCycle(int cycle) {
        if (trajectory == null) {
            return;
        }
        if (cycle <= 0 || cycle > trajectory.size()) {
            return;
        }
        moveTo(trajectory.get(cycle - 1));
    }

    /**
     * <p>
     * setTorsion</p>
     *
     * @param torsion a {@link ffx.potential.bonded.Torsion} object.
     */
    public void setTorsion(Torsion torsion) {
        torsions.add(torsion);
        Atom a14 = torsion.get1_4(this);
        if (a14 != null) {
            // 1-5 atoms will be bonded to the 1-4 atom.
            if (a14.getBonds() != null) {
                for (Bond b : a14.getBonds()) {
                    Atom a15 = b.get1_2(a14);
                    // Do not include the 1-3 atom
                    if (!torsion.containsAtom(a15)) {
                        // Do not include the 1-5 atom more than once (rings).
                        if (!one_5s.contains(a15)) {
                            set1_5(a15);
                        }
                    }
                }
            }
        }
    }

    /**
     * <p>
     * setXYZGradient</p>
     *
     * @param x a double.
     * @param y a double.
     * @param z a double.
     */
    public void setXYZGradient(double x, double y, double z) {
        xyzGradient[0] = x;
        xyzGradient[1] = y;
        xyzGradient[2] = z;
    }

    /**
     * <p>
     * addToXYZGradient</p>
     *
     * @param x a double.
     * @param y a double.
     * @param z a double.
     */
    public void addToXYZGradient(double x, double y, double z) {
        xyzGradient[0] += x;
        xyzGradient[1] += y;
        xyzGradient[2] += z;
    }

    /**
     * <p>
     * getXYZGradient</p>
     *
     * @param x an array of double.
     */
    public void getXYZGradient(double x[]) {
        if (x == null) {
            x = new double[3];
        }
        x[0] = xyzGradient[0];
        x[1] = xyzGradient[1];
        x[2] = xyzGradient[2];
    }

    /**
     * <p>
     * setGlobalMultipole</p>
     *
     * @param dipole an array of double.
     * @param quadrupole an array of double.
     */
    public void setGlobalMultipole(double dipole[], double quadrupole[][]) {
        if (globalDipole == null) {
            globalDipole = new double[3];
        }
        if (globalQuadrupole == null) {
            globalQuadrupole = new double[3][3];
        }
        for (int i = 0; i < 3; i++) {
            globalDipole[i] = dipole[i];
            System.arraycopy(quadrupole[i], 0, globalQuadrupole[i], 0, 3);
        }
    }

    /**
     * <p>
     * Setter for the field <code>multipoleType</code>.</p>
     *
     * @param multipoleType a {@link ffx.potential.parameters.MultipoleType}
     * object.
     * @param multipoleReferenceSites an array of
     * {@link ffx.potential.bonded.Atom} objects.
     */
    public void setMultipoleType(MultipoleType multipoleType, Atom[] multipoleReferenceSites) {
        this.multipoleType = multipoleType;
        this.multipoleReferenceSites = multipoleReferenceSites;
    }

    /**
     * <p>
     * Setter for the field <code>polarizeType</code>.</p>
     *
     * @param polarizeType a {@link ffx.potential.parameters.PolarizeType}
     * object.
     */
    public void setPolarizeType(PolarizeType polarizeType) {
        this.polarizeType = polarizeType;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setSelected(boolean b) {
        selected = b;
    }

    // Vector Methods
    /**
     * <p>
     * setSphereVisible</p>
     *
     * @param sphereVisible a boolean.
     * @param newShapes a {@link java.util.List} object.
     */
    public void setSphereVisible(boolean sphereVisible, List<BranchGroup> newShapes) {
        if (!sphereVisible) {
            // Make this atom invisible.
            if (branchGroup != null) {
                sphere.setPickable(false);
                sphere.setAppearance(RendererCache.nullAp);
            }
        } else {
            // Make this atom visible.
            if (branchGroup == null) {
                initSphere(newShapes);
            }
            sphere.setAppearance(appearance);
            sphere.setPickable(true);
            updateSphere();
        }
    }

    /**
     * {@inheritDoc}
     *
     * Polymorphic setView method.
     */
    @Override
    public void setView(ViewModel newViewModel, List<BranchGroup> newShapes) {
        switch (newViewModel) {
            // case INVISIBLE through case TUBE change the "ViewModel"
            case INVISIBLE:
                viewModel = ViewModel.INVISIBLE;
                setSphereVisible(false, newShapes);
                break;
            case WIREFRAME:
                viewModel = ViewModel.INVISIBLE;
                setSphereVisible(false, newShapes);
                break;
            case SPACEFILL:
                viewModel = ViewModel.SPACEFILL;
                scale = AtomVDW.get(atomType.atomicNumber) * radius;
                setSphereVisible(true, newShapes);
                break;
            case RMIN:
                viewModel = ViewModel.RMIN;
                scale = 1.0 * radius;
                setSphereVisible(true, newShapes);
                break;
            case BALLANDSTICK:
                viewModel = ViewModel.BALLANDSTICK;
                scale = AtomVDW.get(atomType.atomicNumber) / 5.0d * radius;
                setSphereVisible(true, newShapes);
                break;
            case TUBE:
                viewModel = ViewModel.TUBE;
                scale = RendererCache.radius * 0.2d;
                setSphereVisible(true, newShapes);
                break;
            case SHOWHYDROGENS:
                if (atomType.atomicNumber == 1) {
                    return;
                }
                break;
            case HIDEHYDROGENS:
                if (atomType.atomicNumber == 1) {
                    viewModel = ViewModel.INVISIBLE;
                    setSphereVisible(false, newShapes);
                    return;
                }
                break;
            case RESTRICT:
                if (!isSelected()) {
                    viewModel = ViewModel.INVISIBLE;
                    setSphereVisible(false, newShapes);
                    return;
                }
                break;
            case DETAIL:
                int newdetail = RendererCache.detail;
                if (newdetail != detail) {
                    detail = newdetail;
                    if (sphere != null) {
                        Geometry geom = RendererCache.getSphereGeom(detail);
                        sphere.removeAllGeometries();
                        sphere.addGeometry(geom);
                    }
                }
                double newradius = RendererCache.radius;
                if (newradius != radius) {
                    radius = newradius;
                    setView(viewModel, newShapes);
                }
                break;
            // Polygon Appearance Selection
            case FILL:
            case POINTS:
            case LINES:
                polygonType = newViewModel;
                appearance = RendererCache.appearanceFactory(currentCol,
                        polygonType);
                if (viewModel != ViewModel.INVISIBLE) {
                    setSphereVisible(true, newShapes);
                }
                break;
        }
    }

    /**
     * <p>
     * setXYZIndex</p>
     *
     * @param index a int.
     */
    public void setXYZIndex(int index) {
        xyzIndex = index;
    }

    /**
     * <p>
     * toMultipoleString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String toMultipoleString() {
        {
            if (multipoleType == null || globalDipole == null || globalQuadrupole == null) {
                return null;
            }
            StringBuilder multipoleBuffer = new StringBuilder(toString());
            multipoleBuffer.append(String.format("\n%11$s % 7.5f\n" + "%11$s % 7.5f % 7.5f % 7.5f\n" + "%11$s % 7.5f\n" + "%11$s % 7.5f % 7.5f\n" + "%11$s % 7.5f % 7.5f % 7.5f",
                    multipoleType.charge, globalDipole[0], globalDipole[1],
                    globalDipole[2], globalQuadrupole[0][0],
                    globalQuadrupole[1][0], globalQuadrupole[1][1],
                    globalQuadrupole[2][0], globalQuadrupole[2][1],
                    globalQuadrupole[2][2], "                 "));
            return multipoleBuffer.toString();
        }
    }

    /**
     * <p>
     * toShortString</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String toShortString() {
        if (shortString == null) {
            shortString = format("%d-%s", xyzIndex, getName());
        }
        return shortString;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        if (altLoc != null && altLoc != ' ') {
            return String.format("%s %7d-%s %s %d (%7.2f,%7.2f,%7.2f) %s", altLoc, xyzIndex, getName(),
                    resName, resSeq, xyz[0], xyz[1], xyz[2], segID);
        }
        if (resName == null) {
            return String.format("%7d-%s (%7.2f,%7.2f,%7.2f)", xyzIndex, getName(), xyz[0], xyz[1], xyz[2]);
        }
        return String.format("%7d-%s %s %d (%7.2f,%7.2f,%7.2f) %s", xyzIndex, getName(), resName, resSeq,
                xyz[0], xyz[1], xyz[2], segID);
    }

    /**
     * {@inheritDoc}
     *
     * Java3D transforms are not updated unless they are visible. This allows
     * better performance for rendering partial structures during an interactive
     * dynamics run or during trajectory playback.
     */
    @Override
    public void update() {
        if (stale) {
            updateSphere();
            stale = false;
        }
    }

    /**
     * <p>
     * updateSphere</p>
     */
    public void updateSphere() {
        if (branchGroup != null && viewModel != ViewModel.INVISIBLE) {
            vector3d.set(xyz);
            transform3D.setTranslation(vector3d);
            transform3D.setScale(scale);
            transformGroup.setTransform(transform3D);
        }
    }

    /**
     * Element symbols for the first 109 elements.
     */
    public static enum ElementSymbol {

        H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg,
        Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga,
        Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In,
        Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho,
        Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At,
        Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr,
        Rf, Db, Sg, Bh, Hs, Mt;
    }
}
