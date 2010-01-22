/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Logger;

import javax.media.j3d.Appearance;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.Geometry;
import javax.media.j3d.J3DGraphics2D;
import javax.media.j3d.Material;
import javax.media.j3d.Node;
import javax.media.j3d.Shape3D;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;
import javax.vecmath.Color3f;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;
import java.util.Arrays;
import java.util.Vector;

/**
 * The Atom class represents a single atom and defines its alternate
 * conformations and molecular mechanics atom type.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class Atom extends MSNode implements Comparable<Atom> {

    private Logger logger = Logger.getLogger(Atom.class.getName());

    /**
     * Implementation of the Comparable interface.
     * @param a The atom to compare to.
     * @return If <code>this == a</code> return  0<br>
     *         If <code>a < this</code>  return -1<br>
     *         If <code>a > this</code>  return  1<br>
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
        assert (xyzIndex != a0);
        return 0;
    }

    public enum AtomName {

        H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr;
    }
    /**
     *
     */
    private static final long serialVersionUID = 1L;
    // Theses "max" values are used for relative vector display
    private static Point3d point3d = new Point3d();
    private static Point2d point2d = new Point2d();
    private static double[] y = {0.0d, 1.0d, 0.0d};
    public static Hashtable<Integer, Color3f> AtomColor = new Hashtable<Integer, Color3f>();
    public static Hashtable<Integer, Float> AtomVDW = new Hashtable<Integer, Float>();
    /**
     * Hybridizations
     */
    public static final int SP = 2, SP2 = 3, SP3 = 4;
    public final static Hashtable<String, Integer> hybridTable = new Hashtable<String, Integer>();

    static {
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
    }

    static {
        AtomVDW.put(0, 1.0f);
        AtomVDW.put(1, 1.20f);
        AtomVDW.put(2, 1.22f);
        AtomVDW.put(3, 0.78f);
        AtomVDW.put(4, 0.34f);
        AtomVDW.put(5, 2.08f);
        AtomVDW.put(6, 1.85f);
        AtomVDW.put(7, 1.54f);
        AtomVDW.put(8, 1.40f);
        AtomVDW.put(9, 1.35f);
        AtomVDW.put(10, 1.60f);
        AtomVDW.put(11, 0.98f);
        AtomVDW.put(12, 0.78f);
        AtomVDW.put(13, 0.57f);
        AtomVDW.put(14, 2.00f);
        AtomVDW.put(15, 1.90f);
        AtomVDW.put(16, 1.85f);
        AtomVDW.put(17, 1.81f);
        AtomVDW.put(18, 1.91f);
        AtomVDW.put(19, 1.33f);
        AtomVDW.put(20, 1.06f);
        AtomVDW.put(21, 0.91f);
        AtomVDW.put(22, 0.83f);
        AtomVDW.put(23, 0.82f);
        AtomVDW.put(24, 2.0f);
        AtomVDW.put(25, 2.0f);
        AtomVDW.put(26, 2.0f);
        AtomVDW.put(27, 2.0f);
        AtomVDW.put(28, 2.0f);
        AtomVDW.put(29, 2.0f);
        AtomVDW.put(30, 2.0f);
        AtomVDW.put(31, 2.0f);
        AtomVDW.put(32, 2.0f);
        AtomVDW.put(33, 2.0f);
        AtomVDW.put(34, 2.0f);
        AtomVDW.put(35, 1.95f);
        AtomVDW.put(36, 1.89f);
        for (int i = 37; i < 109; i++) {
            AtomVDW.put(i, 2.0f);
        }
    }

    static {
        hybridTable.put("1", 1);
        hybridTable.put("6", 4);
        hybridTable.put("7", 3);
        hybridTable.put("8", 2);
        hybridTable.put("15", 4);
        hybridTable.put("16", 2);
        hybridTable.put("19", 0);
        hybridTable.put("26", 8);
    }
    /**
     * Either the PDB "name" record or the molecular mechanics atom type name.
     *
     * @since 1.0
     */
    private String name = null;
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
    private String chainID = null;
    /**
     * Specifies which alternate conformation is currently active. Getter and
     * Setter methods respect the altID index.
     *
     * @since 1.0
     */
    private int altID = 0;
    /**
     * Array of altLoc identifiers defined for this atom.
     *
     * @since 1.0
     */
    private Character altLoc[] = new Character[1];
    /**
     * Array of XYZ coordiantes for each altLoc.
     *
     * @since 1.0
     */
    private double xyz[][] = new double[1][];
    /**
     * Array of XYZ gradients for each altLoc.
     *
     * @since 1.0
     */
    private double xyzGradient[][] = new double[1][];
    /**
     * Array of occupancy values for each altLoc.
     *
     * @since 1.0
     */
    private double occupancy[] = new double[1];
    /**
     * Array of occupancy gradients for each altLoc.
     *
     * @since 1.0
     */
    private double occupancyGradient[] = new double[1];
    /**
     * Array of tempFactor values for each altLoc.
     *
     * @since 1.0
     */
    private double tempFactor[] = new double[1];
    /**
     * Array fo tempFactorGradients
     *
     * @since 1.0
     */
    private double tempFactorGradient[] = new double[1];
    /**
     * List of anisou tensors for each altLoc.
     *
     * @since 1.0
     */
    private double anisou[][] = new double[1][];
    /**
     * List of anisou tensors for each altLoc.
     *
     * @since 1.0
     */
    private double anisouGradient[][] = new double[1][];
    private ArrayList<Vector3d> trajectory;
    // Molecular Mechanics Info
    private AtomType atomType = null;
    private MultipoleType multipoleType = null;
    private PolarizeType polarizeType = null;
    private Atom[] multipoleReferenceSites = null;
    private double globalDipole[] = null;
    private double globalQuadrupole[][] = null;
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

    /**
     * Default constructor.
     *
     * @param name The Atom's PDB name.
     *
     * @since 1.0
     */
    public Atom(String name) {
        super(name, 1);
        this.name = name;
        currentCol = previousCol = RendererCache.toAtomColor(name);
        colorModel = ColorModel.CPK;
    }

    /**
     * Constructor
     *
     * @param name The Atom's PDB name.
     * @param atomType Molecular mechanics atom type.
     * @param xyz Cartesian coordinates.
     *
     * @since 1.0
     */
    public Atom(String name, AtomType atomType, double[] xyz) {
        this(name);
        this.atomType = atomType;
        this.xyz[0] = xyz;
        this.xyzGradient[0] = new double[3];
        setAllowsChildren(false);
        if (atomType != null) {
            currentCol = previousCol = AtomColor.get(atomType.atomicNumber);
        }
    }

    /**
     * Constructor used when parsing XYZ files.
     * 
     * @param xyzIndex Contiguous, unique atom index between 0..nAtoms.
     * @param name The Atom's molecular mechanics name.
     * @param atomType Molecular mechanics atom type.
     * @param xyz Cartesian coordinates.
     * 
     * @since 1.0
     */
    public Atom(int xyzIndex, String name, AtomType atomType, double[] xyz) {
        this(name, atomType, xyz);
        this.xyzIndex = xyzIndex;
    }

    /**
     * Constructor used when parsing PDB files.
     */
    public Atom(int xyzIndex, String name,
            Character altLoc, double[] d, String resName, int resSeq,
            String chainID, double occupancy, double tempFactor) {
        this(xyzIndex, name, null, d);
        this.resName = resName;
        this.resSeq = resSeq;
        this.chainID = chainID;
        this.altLoc[0] = altLoc;
        this.occupancy[0] = occupancy;
        this.tempFactor[0] = tempFactor;
        altID = 0;
    }

    public void addTrajectoryCoords(Vector3d coords, int position) {
        if (trajectory == null) {
            trajectory = new ArrayList<Vector3d>();
            int i = 0;
            for (double[] x : xyz) {
                trajectory.add(i++, new Vector3d(x));
            }
        }
        trajectory.add(position, coords);
    }

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
     * Overidden equals method.
     *
     * @param object
     *            The Object to compare with <b>this</b>
     * @return True if <b>this</b> atom and object do not reference the same
     *         object, are of the same class, and have the same id
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

        return (other.chainID != null && other.chainID.equals(chainID)
                && other.resName != null && other.resName.equals(resName)
                && other.resSeq == resSeq
                && other.name != null && other.name.equals(name));
    }

    public ArrayList<Angle> getAngles() {
        return angles;
    }

    public ArrayList<Atom> get1_5s() {
        return one_5s;
    }

    public Appearance getAtomAppearance() {
        if (appearance == null) {
            appearance = RendererCache.appearanceFactory(currentCol,
                    polygonType);
        }
        return appearance;
    }

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

    @Override
    public ArrayList<Atom> getAtomList() {
        ArrayList<Atom> atoms = new ArrayList<Atom>();
        atoms.add(this);
        return atoms;
    }

    public AtomType getAtomType() {
        return atomType;
    }

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

    public double getBornRadius() {
        return bornRadius;
    }

    public double getBornVolume() {
        return 2.0;
    }

    /**
     * Get the chain name
     *
     * @return String
     */
    public String getChain() {
        if (chainID != null) {
            return chainID;
        }
        Polymer p = (Polymer) this.getMSNode(Polymer.class);
        if (p == null) {
            return null;
        }
        chainID = p.getName();
        return chainID;
    }

    /**
     * Gets the partial atomic charge
     *
     * @return partial atomic charge
     */
    public double getCharge() {
        return 1.0;
    }

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

    public void getForce(double[] t) {
        if (xyzGradient == null || t == null) {
            return;
        }
        double g[] = xyzGradient[altID];
        t[0] = g[0];
        t[1] = g[1];
        t[2] = g[2];
    }

    public double getGPol() {
        return 2.0;
    }

    /**
     * Gets the energy gradient
     *
     * @return energy gradient
     */
    // public Vector3d getGradient(){ return new Vector3d(gradient); }
    /**
     * Gets the Atomic Hybridization
     *
     * @return Atomic Hybridization
     */
    public int getHybridization() {
        return atomType.valence;
    }

    /**
     * Gets the ID
     *
     * @return ID
     */
    public String getID() {
        return name;
    }

    /**
     * Gets the atom ID
     *
     * @return atom ID
     */
    public String getIdent() {
        String s = new String(atomType.environment);
        return s;
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

    public Atom[] getMultipoleReferenceSites() {
        return this.multipoleReferenceSites;
    }

    public MultipoleType getMultipoleType() {
        return multipoleType;
    }

    public PolarizeType getPolarizeType() {
        return polarizeType;
    }

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

    public final int getNumDihedrals() {
        if (torsions == null) {
            return 0;
        }
        return torsions.size();
    }

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

    public int getResidueNumber() {
        return resSeq;
    }

    /**
     * Gets the Sigma value
     *
     * @return Sigma value
     */
    public double getSigma() {
        return 1.0;
    }

    public Vector3d getTrajectoryCoords(int position) {
        return trajectory.get(position);
    }

    public int getTrajectoryLength() {
        return trajectory.size();
    }

    public int getType() {
        return atomType.type;
    }

    /**
     * Gets the Atom's Cartesian Coordinates return The Cartesian Coordinates
     */
    public void getV3D(Vector3d temp) {
        temp.set(xyz[altID]);
    }

    // public final Vector3f getnewV3D(){ return new Vector3f(v3f); }
    /**
     * Gets the van der Waals radius return van der Waals radius
     */
    public double getVDWR() {
        return 1.0;
    }

    public void setAltLoc(Character a) {
        for (int i = 0; i < altLoc.length; i++) {
            if (altLoc[i] == a) {
                altID = i;
            }
        }
    }

    public Character getAltLoc() {
        return altLoc[altID];
    }

    public Character[] getAltLocs() {
        return altLoc;
    }

    public void addAltLoc(Character a, double[] xyz, double occupancy, double bFactor) {
        altID = -1;
        for (int i = 0; i < altLoc.length; i++) {
            if (altLoc[i] == a) {
                altID = i;
            }
        }
        /**
         * If we've seen this altLoc character already, overwrite the previous
         * values.
         */
        if (altID != -1) {
            this.xyz[altID] = xyz;
            this.occupancy[altID] = occupancy;
            this.tempFactor[altID] = bFactor;
        } else {
            /**
             * This is a new altLoc so we need to copy over the current arrays to
             * longer versions and store the new altLoc information.
             */
            altID = altLoc.length;
            altLoc = Arrays.copyOf(altLoc, altID + 1);
            altLoc[altID] = a;
            this.xyz = Arrays.copyOf(this.xyz, altID + 1);
            this.xyz[altID] = xyz;
            this.occupancy = Arrays.copyOf(this.occupancy, altID + 1);
            this.occupancy[altID] = occupancy;
            this.tempFactor = Arrays.copyOf(this.tempFactor, altID + 1);
            this.tempFactor[altID] = bFactor;
            xyzGradient = Arrays.copyOf(xyzGradient, altID + 1);
            xyzGradient[altID] = new double[3];
            occupancyGradient = Arrays.copyOf(occupancyGradient, altID + 1);
            tempFactorGradient = Arrays.copyOf(tempFactorGradient, altID + 1);
        }
    }

    /**
     * Gets the x coordinate
     *
     * @return x coordinate
     */
    public double getX() {
        return xyz[altID][0];
    }

    public void getXYZ(double[] x) {
        double t[] = xyz[altID];
        x[0] = t[0];
        x[1] = t[1];
        x[2] = t[2];
    }

    public double[] getXYZ() {
        return xyz[altID];
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
        return xyz[altID][1];
    }

    /**
     * Gets the z coordinate
     *
     * @return z coordinate
     */
    public final double getZ() {
        return xyz[altID][2];
    }

    @Override
    public final int hashCode() {
        int ret = HashCodeUtil.hash(HashCodeUtil.ATOMSEED, xyzIndex);
        return ret;
    }

    /**
     * Create the Sphere Java3D objects.
     *
     * @param newShapes
     *            List
     */
    private void initSphere(List<BranchGroup> newShapes) {
        if (appearance == null) {
            appearance = RendererCache.appearanceFactory(currentCol,
                    ViewModel.FILL);
        }
        if (transform3D == null) {
            transform3D = RendererCache.transform3DFactory(new Vector3d(xyz[altID]),
                    scale);
        } else {
            transform3D.setTranslation(new Vector3d(xyz[altID]));
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
     * @param a
     *            Atom to check
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

    public boolean is_1_3(Atom atom) {
        for (Angle a : angles) {
            if (a.get1_3(atom) == this) {
                return true;
            }
        }
        return false;
    }

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

    public boolean isTrigonal() {
        if (bonds != null && bonds.size() == 3) {
            return true;
        }
        return false;
    }

    public boolean isStale() {
        return stale;
    }

    // True if this Atom's Sphere or Vector is visible
    public boolean isVisible() {
        return (viewModel != ViewModel.INVISIBLE);
    }

    /**
     * Add a vector to the Atom's current position vector
     *
     * @param d Vector to add to the current position
     */
    public void move(double[] d) {
        double x[] = xyz[altID];
        x[0] += d[0];
        x[1] += d[1];
        x[2] += d[2];
        stale = true;
    }

    public void moveTo(double a, double b, double c) {
        double x[] = xyz[altID];
        x[0] = a;
        x[1] = b;
        x[2] = c;
        stale = true;
    }

    public void setXYZ(double xyz[]) {
        this.xyz[altID] = xyz;
    }

    /**
     * Moves the atom to the specified location
     *
     * @param d Location to move <b>this</b> Atom to
     */
    public void moveTo(double[] d) {
        moveTo(d[0], d[1], d[2]);
    }

    public void moveTo(Vector3d v) {
        moveTo(v.x, v.y, v.z);
    }

    /**
     * Prints the atom identity and Cartesian coordinates to the logger.
     */
    @Override
    public final void print() {
        logger.info(toString());
    }

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

    public void setAtomType(AtomType atomType) {
        this.atomType = atomType;
    }

    public void setTempFactor(double tempFactor) {
        this.tempFactor[altID] = tempFactor;
    }

    public double getTempFactor() {
        return tempFactor[altID];
    }

    public void setTempFactorGradient(double tempFactorGradient) {
        this.tempFactorGradient[altID] = tempFactorGradient;
    }

    public double getTempFactorGradient() {
        return tempFactorGradient[altID];
    }

    public void setOccupancy(double occupancy) {
        this.occupancy[altID] = occupancy;
    }

    public double getOccupancy() {
        return occupancy[altID];
    }

    public void setOccupancyGradient(double occupancyGradient) {
        this.occupancyGradient[altID] = occupancyGradient;
    }

    public double getOccupancyGradient() {
        return occupancyGradient[altID];
    }

    public void setAnisou(double[] anisou) {
        /**
         * Ensure that we have space for it.
         */
        if (this.anisou.length <= altID) {
            this.anisou = Arrays.copyOf(this.anisou, altID + 1);
            this.anisouGradient = Arrays.copyOf(this.anisouGradient, altID + 1);
            anisouGradient[altID] = new double[6];
        }
        this.anisou[altID] = anisou;
    }

    public double[] getAnisou() {
        if (anisou.length > altID) {
            return anisou[altID];
        } else {
            return null;
        }
    }

    public void setAnisouGradient(double[] anisou) {
        /**
         * Ensure that we have space for it.
         */
        assert (this.anisouGradient.length <= altID);

        this.anisouGradient[altID] = anisou;
    }

    public double[] getAnisouGradient() {
        return anisouGradient[altID];
    }

    public void setAngle(Angle a) {
        angles.add(a);
    }

    /**
     * Specify that <b>this</b> Atom is part of a Bond
     *
     * @param b
     *            Bond that <b>this</b> Atom is part of
     */
    public void setBond(Bond b) {
        bonds.add(b);
    }

    public void set1_5(Atom a) {
        one_5s.add(a);
    }

    /**
     * Set the effective Born Radius.
     *
     * @param bornRadius
     */
    public void setBornRadius(double bornRadius) {
        this.bornRadius = bornRadius;
    }

    /**
     * Polymorphic setColor method.
     *
     * @param newColorModel
     *            ColorModel
     * @param newCol
     *            Color3f
     * @param newMat
     *            Material
     */
    @Override
    public void setColor(ColorModel newColorModel, Color3f newCol,
            Material newMat) {
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

    public void setCurrentCycle(int cycle) {
        if (trajectory == null) {
            return;
        }
        if (cycle <= 0 || cycle > trajectory.size()) {
            return;
        }
        moveTo(trajectory.get(cycle - 1));
    }

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

    public void setXYZGradient(double x, double y, double z) {
        double g[] = xyzGradient[altID];
        g[0] = x;
        g[1] = y;
        g[2] = z;
    }

    public void addToXYZGradient(double x, double y, double z) {
        double g[] = xyzGradient[altID];
        g[0] += x;
        g[1] += y;
        g[2] += z;
    }

    public void getXYZGradient(double x[]) {
        if (x == null) {
            x = new double[3];
        }
        double g[] = xyzGradient[altID];
        x[0] = g[0];
        x[1] = g[1];
        x[2] = g[2];
    }

    public void setGlobalMultipole(double dipole[], double quadrupole[][]) {
        if (globalDipole == null) {
            globalDipole = new double[3];
        }
        if (globalQuadrupole == null) {
            globalQuadrupole = new double[3][3];
        }
        for (int i = 0; i < 3; i++) {
            globalDipole[i] = dipole[i];
            for (int j = 0; j < 3; j++) {
                globalQuadrupole[i][j] = quadrupole[i][j];
            }
        }
    }

    public void setMultipoleType(MultipoleType multipoleType,
            Atom[] multipoleReferenceSites) {
        this.multipoleType = multipoleType;
        this.multipoleReferenceSites = multipoleReferenceSites;
    }

    public void setPolarizeType(PolarizeType polarizeType) {
        this.polarizeType = polarizeType;
    }

    @Override
    public void setSelected(boolean b) {
        selected = b;
    }

    // Vector Methods
    public void setSphereVisible(boolean sphereVisible,
            List<BranchGroup> newShapes) {
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
     * Polymorphic setView method.
     *
     * @param newViewModel
     *            ViewModel
     * @param newShapes
     *            List
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

    public void setXYZIndex(int index) {
        xyzIndex = index;
    }

    public String toMultipoleString() {
        {
            if (multipoleType == null || globalDipole == null || globalQuadrupole == null) {
                return null;
            }
            StringBuffer multipoleBuffer = new StringBuffer(toString());
            multipoleBuffer.append(String.format("\n%11$s % 7.5f\n" + "%11$s % 7.5f % 7.5f % 7.5f\n" + "%11$s % 7.5f\n" + "%11$s % 7.5f % 7.5f\n" + "%11$s % 7.5f % 7.5f % 7.5f",
                    multipoleType.charge, globalDipole[0], globalDipole[1],
                    globalDipole[2], globalQuadrupole[0][0],
                    globalQuadrupole[1][0], globalQuadrupole[1][1],
                    globalQuadrupole[2][0], globalQuadrupole[2][1],
                    globalQuadrupole[2][2], "                 "));
            return multipoleBuffer.toString();
        }
    }
    private String shortString = null;

    public String toShortString() {
        if (shortString == null) {
            shortString = new String("" + xyzIndex + "-" + name);
        }
        return shortString;
    }

    /**
     * @return The string: "INDEX - ID (X, Y, Z)"
     */
    @Override
    public String toString() {
        double x[] = xyz[altID];
        return String.format("%7d-%s (%7.2f,%7.2f,%7.2f)", xyzIndex, name,
                x[0], x[1], x[2]);
    }

    /**
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
    private Vector3d vector3d = new Vector3d();

    public void updateSphere() {
        if (branchGroup != null && viewModel != ViewModel.INVISIBLE) {
            vector3d.set(xyz[altID]);
            transform3D.setTranslation(vector3d);
            transform3D.setScale(scale);
            transformGroup.setTransform(transform3D);
        }
    }
}
