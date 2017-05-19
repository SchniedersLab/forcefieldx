/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import static java.lang.String.format;

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

import org.apache.commons.lang.ArrayUtils;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;
import ffx.potential.extended.ExtendedVariable;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.VDWType;

import static ffx.potential.extended.SBLogger.SB;
import static ffx.utilities.HashCodeUtil.SEED;
import static ffx.utilities.HashCodeUtil.hash;

/**
 * The Atom class represents a single atom and defines its alternate
 * conformations and molecular mechanics atom type.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
@SuppressWarnings({"serial", "CloneableImplementsClone"})
public class Atom extends MSNode implements Comparable<Atom> {

    private static final Logger logger = Logger.getLogger(Atom.class.getName());

    public enum Resolution {
        FIXEDCHARGE, AMOEBA;
    }

    private Resolution resolution = Resolution.AMOEBA;

    public void setResolution(Resolution resolution) {
        this.resolution = resolution;
    }

    public Resolution getResolution() {
        return resolution;
    }

    public enum Indexing {
        XYZ, PERSIST
    }

    private boolean isBackground = false;

    public void setBackground() {
        isBackground = true;
    }

    public boolean isBackground() {
        return isBackground;
    }

    private Atom[] axisAtoms = null;
    private int[] axisAtomIndices = null;
    private MultipoleFrameDefinition frameDefinition = null;

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
        AtomColor = new HashMap<>();
        AtomVDW = new HashMap<>();
        hybridTable = new HashMap<>();
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

    /**
     * Contiguous atom index ranging from 1..nAtoms.
     *
     * @since 1.0
     */
    private int xyzIndex = -1;
    /**
     * Persistent (unmodifiable) indexing alternative to xyzIndex.
     */
    public final int persistentIndex;
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
    private final double xyz[] = new double[3];
    /**
     * Array of velocities
     */
    private final double velocity[] = new double[3];
    private final double acceleration[] = new double[3];
    private final double previousAcceleration[] = new double[3];
    /**
     * Array of XYZ coordinates for the electron (van der Waals) centers of each
     * atom: if null, methods will refer to xyz.
     *
     * @since 1.0
     */
    private double redXYZ[];
    /**
     * Array of XYZ gradient.
     *
     * @since 1.0
     */
    private final double xyzGradient[] = new double[3];
    /**
     * Array of XYZ lambda gradient.
     *
     * @since 1.0
     */
    private final double xyzLambdaGradient[] = new double[3];
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
    private double occupancyVelocity;
    private double occupancyAcceleration;
    private double occupancyPreviousAcceleration;

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
    private double tempFactorVelocity;
    private double tempFactorAcceleration;
    private double tempFactorPreviousAcceleration;

    /**
     * Anisou tensor.
     *
     * @since 1.0
     */
    private double anisou[];
    /**
     * Anisou gradient, velocity, accel and prev accel.
     *
     * @since 1.0
     */
    private double anisouGradient[];
    private double anisouVelocity[];
    private double anisouAcceleration[];
    private double anisouPreviousAcceleration[];

    /**
     * If use is true, this atom should be included in target functions.
     */
    private boolean use = true;
    /**
     * If active is true, the coordinates of this atom can be modified.
     */
    private boolean active = true;
    /**
     * If built is true, this atom was built during the parsing of a file.
     */
    private boolean built = false;
    /**
     * True if this Atom is a HETATM.
     */
    private boolean hetatm = false;
    /**
     * True if this Atom is a member of modified residue.
     */
    private boolean modres = false;
    /**
     * If electrostatics is true, include the charge, multipole and/or
     * polarizability in electrostatics calculations.
     */
    private boolean electrostatics = true;
    private String segID = null;
    private double formFactorWidth = 3.5;
    private double formFactorWidth2 = formFactorWidth * formFactorWidth;
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
    private boolean esvTerm = false;
    private ExtendedVariable esv = null;
    private Double scaledPolarizability = null;
    private Double unscaledPolarizability = null;
    private int moleculeNumber = 0;
    // solvation
    private double bornRadius;
    // Connectivity information.
    private final ArrayList<Bond> bonds = new ArrayList<>();
    private final ArrayList<Angle> angles = new ArrayList<>();
    private final ArrayList<Torsion> torsions = new ArrayList<>();
    private final ArrayList<Atom> one_5s = new ArrayList<>();
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
    private final Vector3d vector3d = new Vector3d();
    /* Extended System handling */
    private MultipoleType esvMultipole;
    private MultipoleType esvMultipoleDot;

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
        redXYZ = null;
        persistentIndex = MolecularAssembly.persistentAtomIndexer++;
        //this.atomSerial = atomSerialCount.getAndIncrement();
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
        System.arraycopy(xyz, 0, this.xyz, 0, 3);
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

    public Atom(int xyzIndex, String name,
            Character altLoc, double[] xyz, String resName, int resSeq,
            Character chainID, double occupancy, double tempFactor,
            String segID, boolean built) {
        this(xyzIndex, name, altLoc, xyz, resName, resSeq, chainID, occupancy, tempFactor, segID);
        this.built = built;
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
        return 0;
    }

    /*
    @Override
    public int compareTo(Atom other) {
        return Objects.compare(this, other, indexComparator);
    }
    public static IndexComparator indexComparator = new IndexComparator();
    public static class IndexComparator implements Comparator<Atom> {
        private IndexComparator() {}
        @Override
        public int compare(Atom me, Atom u) {
            int mi = me.getIndex();
            int ui = u.getIndex();
            // There should not be duplicate, identical atom objects.
            assert(mi != ui);
            return Integer.compare(mi, ui);
        }
    }   */
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
     * setModRes</p>
     *
     * @param modres a boolean.
     */
    public void setModRes(boolean modres) {
        this.modres = modres;
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
     * isModRes</p>
     *
     * @return a boolean.
     */
    public boolean isModRes() {
        return modres;
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
        Atom atom = new Atom(getIndex(), getName(), getAltLoc(), coords,
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
        return atomType.atomicNumber == 1;
    }

    /**
     * If true, this atom should be used in potential energy functions.
     *
     * @return true if this atom should be included in the potential energy.
     */
    public boolean getUse() {
        return use;
    }

    /**
     * If true, this atom was built during PDB file parsing.
     *
     * @return true if this atom was built during parsing of a PDB file.
     */
    public boolean getBuilt() {
        return built;
    }

    public void setBuilt(boolean built) {
        this.built = built;
    }

    /**
     * If true, this atom should be used in potential energy functions.
     *
     * @param use
     */
    public void setUse(boolean use) {
        this.use = use;
    }

    /**
     * If active, the coordinates of this atom can be modified.
     *
     * @return true if this atom's coordinates, b-factors, etc. can be modified.
     */
    public boolean isActive() {
        return active;
    }

    /**
     * If active, the coordinates of this atom can be modified.
     *
     * @param active
     */
    public void setActive(boolean active) {
        this.active = active;
    }

    public void setElectrostatics(boolean electrostatics) {
        this.electrostatics = electrostatics;
    }

    public boolean getElectrostatics() {
        return electrostatics;
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
            g2d.drawString(describe(Atom.Descriptions.XyzIndex_Name), (float) point2d.x,
                    (float) point2d.y);
        }

    }

    /**
     * This definition allows residue name to differ between two equal atoms.
     * Consider replacing with a Comparator<Atom>.
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
     * @param atom2 Atom number 2.
     * @param atom3 Atom number 3.
     * @param atom4 Atom number 4.
     * @return Torsion the Torsion if found, or null if not found.
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
     * The (single) ExtendedVariable of which this atom is a part (or null
     * otherwise).
     */
    public ExtendedVariable getEsv() {
        return esv;
    }

    public final void setEsv(ExtendedVariable esv, MultipoleType type, MultipoleType dotType, Double scaledAlpha) {
        if (esv == null || (esvTerm && esv != this.esv)) {
            SB.crash(format("Error attaching ESV to atom (multiples not supported): %s %s %s\n", this, this.esv, esv));
        }
        this.esv = esv;
        esvMultipole = type;
        esvMultipoleDot = dotType;
        unscaledPolarizability = getPolarizeType().polarizability;
        scaledPolarizability = (scaledAlpha != null) ? scaledAlpha : unscaledPolarizability;
        esvTerm = true;
    }

    public final void setEsv(ExtendedVariable esv, MultipoleType type, MultipoleType dotType) {
        setEsv(esv, type, dotType, null);
    }

    public MultipoleType getEsvMultipole() {
        return (esvTerm) ? esvMultipole : getMultipoleType();
    }

    public MultipoleType getEsvMultipoleDot() {
        return (esvTerm) ? esvMultipoleDot : getMultipoleType();
    }

    public double getScaledPolarizability() {
        return (esvTerm) ? scaledPolarizability : polarizeType.polarizability;
    }

    public double getUnscaledPolarizability() {
        return (esvTerm) ? unscaledPolarizability : polarizeType.polarizability;
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
        if (r != null) {
            return r.getName();
        }
        return null;
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
     * <p>
     * getRedXYZ</p>
     *
     * @param x an array of double.
     */
    public void getRedXYZ(double[] x) {
        if (redXYZ != null) {
            x[0] = redXYZ[0];
            x[1] = redXYZ[1];
            x[2] = redXYZ[2];
        } else {
            getXYZ(x);
        }
    }

    /**
     * <p>
     * getXYZ</p>
     *
     * @param xyz an array of double.
     *
     * @return the original array with updated coordinates, or a new array if
     * xyz was null.
     */
    public double[] getXYZ(double xyz[]) {
        if (xyz == null) {
            xyz = Arrays.copyOf(this.xyz, 3);
        } else {
            System.arraycopy(this.xyz, 0, xyz, 0, 3);
        }
        return xyz;
    }

    /**
     * <p>
     * getRedXYZ</p>
     *
     * @return an array of double.
     */
    public double[] getRedXYZ() {
        if (active) {
            return redXYZ == null ? xyz : redXYZ;
        }
        if (redXYZ != null) {
            return Arrays.copyOf(redXYZ, 3);
        }
        return Arrays.copyOf(xyz, 3);
    }

    /**
     * Gets the x coordinate
     *
     * @return x coordinate
     */
    public final double getX() {
        return xyz[0];
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
     * Gets the reduced x coordinate (van der Waals center). Will refer to xyz[]
     * if not initialized.
     *
     * @return Reduced x coordinate
     */
    public final double getRedX() {
        return redXYZ == null ? xyz[0] : redXYZ[0];
    }

    /**
     * Gets the reduced y coordinate (van der Waals center). Will refer to xyz[]
     * if not initialized.
     *
     * @return Reduced y coordinate
     */
    public final double getRedY() {
        return redXYZ == null ? xyz[1] : redXYZ[1];
    }

    /**
     * Gets the reduced z coordinate (van der Waals center). Will refer to xyz[]
     * if not initialized.
     *
     * @return Reduced z coordinate
     */
    public final double getRedZ() {
        return redXYZ == null ? xyz[2] : redXYZ[2];
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public final int hashCode() {
        //return hash(SEED, xyzIndex);
        int hash = hash(SEED, resName);
        hash = hash(hash, resSeq);
        hash = hash(hash, getName());
        return hash(hash, segID);
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
        if (active) {
            xyz[0] += d[0];
            xyz[1] += d[1];
            xyz[2] += d[2];
            stale = true;
        }
    }

    public void rotate(double[][] d) {
        int rowsInA = xyz.length;
        double columnsInA = d.length; // same as rows in d
        int columnsInB = d[0].length;
        double[][] c = new double[rowsInA][columnsInB];
        for (int i = 0; i < rowsInA; i++) {
            for (int j = 0; j < columnsInB; j++) {
                for (int k = 0; k < columnsInA; k++) {
                    c[i][j] = c[i][j] + xyz[k] * d[k][j];
                }
            }
        }
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
        if (active) {
            xyz[0] = a;
            xyz[1] = b;
            xyz[2] = c;
            stale = true;
        }
    }

    /**
     * <p>
     * setXYZ</p>
     *
     * @param xyz an array of double.
     */
    public void setXYZ(double xyz[]) {
        if (active && xyz != null) {
            System.arraycopy(xyz, 0, this.xyz, 0, 3);
        }
    }

    /**
     *
     * @param velocity an array of double.
     */
    public void setVelocity(double velocity[]) {
        if (active && velocity != null) {
            System.arraycopy(velocity, 0, this.velocity, 0, 3);
        }
    }

    /**
     *
     * @param acceleration an array of double.
     */
    public void setAcceleration(double acceleration[]) {
        if (active && acceleration != null) {
            System.arraycopy(acceleration, 0, this.acceleration, 0, 3);
        }
    }

    /**
     *
     * @param previousAcceleration an array of double.
     */
    public void setPreviousAcceleration(double previousAcceleration[]) {
        if (active && previousAcceleration != null) {
            System.arraycopy(previousAcceleration, 0, this.previousAcceleration, 0, 3);
        }
    }

    /**
     * <p>
     * setXYZ</p>
     *
     * @param redXYZ an array of double.
     */
    public void setRedXYZ(double redXYZ[]) {
        if (active) {
            if (this.redXYZ == null) {
                this.redXYZ = new double[3];
            }
            this.redXYZ[0] = redXYZ[0];
            this.redXYZ[1] = redXYZ[1];
            this.redXYZ[2] = redXYZ[2];
        }
    }

    /**
     * Moves the atom to the specified location
     *
     * @param d Location to move <b>this</b> Atom to
     */
    public void moveTo(double[] d) {
        if (active) {
            moveTo(d[0], d[1], d[2]);
        }
    }

    /**
     * <p>
     * moveTo</p>
     *
     * @param v a {@link javax.vecmath.Vector3d} object.
     */
    public void moveTo(Vector3d v) {
        if (active) {
            moveTo(v.x, v.y, v.z);
        }
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
        if (active) {
            this.tempFactor = tempFactor;
        }
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
     * Getter for the field <code>tempFactorGradient</code>.</p>
     *
     * @return a double.
     */
    public double getTempFactorGradient() {
        return tempFactorGradient;
    }

    /**
     * <p>
     * Getter for the field <code>tempFactorVelocity</code>.</p>
     *
     * @return a double.
     */
    public double getTempFactorVelocity() {
        return tempFactorVelocity;
    }

    /**
     * <p>
     * Getter for the field <code>tempFactorAcceleration</code>.</p>
     *
     * @return a double.
     */
    public double getTempFactorAcceleration() {
        return tempFactorAcceleration;
    }

    /**
     * <p>
     * Getter for the field <code>tempFactorPreviousAcceleration</code>.</p>
     *
     * @return a double.
     */
    public double getTempFactorPreviousAcceleration() {
        return tempFactorPreviousAcceleration;
    }

    /**
     * <p>
     * Setter for the field <code>tempFactorGradient</code>.</p>
     *
     * @param tempFactorGradient a double.
     */
    public void setTempFactorGradient(double tempFactorGradient) {
        if (active) {
            this.tempFactorGradient = tempFactorGradient;
        }
    }

    /**
     * <p>
     * addToTempFactorGradient</p>
     *
     * @param tempFactorGradient a double.
     */
    public void addToTempFactorGradient(double tempFactorGradient) {
        if (active) {
            this.tempFactorGradient += tempFactorGradient;
        }
    }

    /**
     * <p>
     * Setter for the field <code>tempFactorVelocity</code>.</p>
     *
     * @param tempFactorVelocity a double.
     */
    public void setTempFactorVelocity(double tempFactorVelocity) {
        if (active) {
            this.tempFactorVelocity = tempFactorVelocity;
        }
    }

    /**
     * <p>
     * Setter for the field <code>tempFactorAcceleration</code>.</p>
     *
     * @param tempFactorAcceleration a double.
     */
    public void setTempFactorAcceleration(double tempFactorAcceleration) {
        if (active) {
            this.tempFactorAcceleration = tempFactorAcceleration;
        }
    }

    /**
     * <p>
     * Setter for the field <code>tempFactorPreviousAcceleration</code>.</p>
     *
     * @param tempFactorPreviousAcceleration a double.
     */
    public void setTempFactorPreviousAcceleration(double tempFactorPreviousAcceleration) {
        if (active) {
            this.tempFactorPreviousAcceleration = tempFactorPreviousAcceleration;
        }
    }

    /**
     * <p>
     * Setter for the field <code>occupancy</code>.</p>
     *
     * @param occupancy a double.
     */
    public void setOccupancy(double occupancy) {
        if (active) {
            this.occupancy = occupancy;
        }
    }

    /**
     *
     * @param occupancyVelocity a double.
     */
    public void setOccupancyVelocity(double occupancyVelocity) {
        if (active) {
            this.occupancyVelocity = occupancyVelocity;
        }
    }

    /**
     *
     * @param occupancyAcceleration a double.
     */
    public void setOccupancyAcceleration(double occupancyAcceleration) {
        if (active) {
            this.occupancyAcceleration = occupancyAcceleration;
        }
    }

    /**
     *
     * @param occupancyPreviousAcceleration a double.
     */
    public void setOccupancyPreviousAcceleration(double occupancyPreviousAcceleration) {
        if (active) {
            this.occupancyPreviousAcceleration = occupancyPreviousAcceleration;
        }
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
        if (active) {
            this.occupancyGradient = occupancyGradient;
        }
    }

    /**
     * <p>
     * addToOccupancyGradient</p>
     *
     * @param occupancyGradient a double.
     */
    public void addToOccupancyGradient(double occupancyGradient) {
        if (active) {
            this.occupancyGradient += occupancyGradient;
        }
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
     * Getter for the field <code>occupancyVelocity</code>.</p>
     *
     * @return a double.
     */
    public double getOccupancyVelocity() {
        return occupancyVelocity;
    }

    /**
     * <p>
     * Getter for the field <code>occupancyAccelerationy</code>.</p>
     *
     * @return a double.
     */
    public double getOccupancyAcceleration() {
        return occupancyAcceleration;
    }

    /**
     * <p>
     * Getter for the field <code>occupancyPreviousAccelerationy</code>.</p>
     *
     * @return a double.
     */
    public double getOccupancyPreviousAcceleration() {
        return occupancyPreviousAcceleration;
    }

    /**
     * <p>
     * Getter for the field <code>anisou</code>.</p>
     *
     * @param anisou
     * @return an array of double.
     */
    public double[] getAnisou(double[] anisou) {
        if (this.anisou == null) {
            return null;
        } else if (anisou == null) {
            anisou = Arrays.copyOf(this.anisou, 6);
        } else {
            System.arraycopy(this.anisou, 0, anisou, 0, 6);
        }
        return anisou;
    }

    /**
     * <p>
     * Getter for the field <code>anisouGradient</code>.</p>
     *
     * @param anisouGradient
     * @return an array of double.
     */
    public double[] getAnisouGradient(double[] anisouGradient) {
        if (this.anisouGradient == null) {
            return null;
        } else if (anisouGradient == null) {
            anisouGradient = Arrays.copyOf(this.anisouGradient, 6);
        } else {
            System.arraycopy(this.anisouGradient, 0, anisouGradient, 0, 6);
        }
        return anisouGradient;
    }

    /**
     * <p>
     * Getter for the field <code>anisouVelocity</code>.</p>
     *
     * @param anisouVelocity
     * @return an array of double.
     */
    public double[] getAnisouVelocity(double[] anisouVelocity) {
        if (this.anisouVelocity == null) {
            return null;
        } else if (anisouVelocity == null) {
            anisouVelocity = Arrays.copyOf(this.anisouVelocity, 6);
        } else {
            System.arraycopy(this.anisouVelocity, 0, anisouVelocity, 0, 6);
        }
        return anisouVelocity;
    }

    /**
     * <p>
     * Getter for the field <code>anisouAcceleration</code>.</p>
     *
     * @param anisouAcceleration
     * @return an array of double.
     */
    public double[] getAnisouAcceleration(double[] anisouAcceleration) {
        if (this.anisouAcceleration == null) {
            return null;
        } else if (anisouAcceleration == null) {
            anisouAcceleration = Arrays.copyOf(this.anisouAcceleration, 6);
        } else {
            System.arraycopy(this.anisouAcceleration, 0, anisouAcceleration, 0, 6);
        }
        return anisouAcceleration;
    }

    /**
     * <p>
     * Getter for the field <code>anisouPreviousAcceleration</code>.</p>
     *
     * @param anisouPreviousAcceleration
     * @return an array of double.
     */
    public double[] getAnisouPreviousAcceleration(double[] anisouPreviousAcceleration) {
        if (this.anisouPreviousAcceleration == null) {
            return null;
        } else if (anisouPreviousAcceleration == null) {
            anisouPreviousAcceleration = Arrays.copyOf(this.anisouPreviousAcceleration, 6);
        } else {
            System.arraycopy(this.anisouPreviousAcceleration, 0, anisouPreviousAcceleration, 0, 6);
        }
        return anisouPreviousAcceleration;
    }

    /**
     * <p>
     * Setter for the field <code>anisou</code>.</p>
     *
     * @param anisou an array of double.
     */
    public void setAnisou(double[] anisou) {
        if (active) {
            if (anisou == null) {
                this.anisou = null;
            } else if (this.anisou == null) {
                this.anisou = Arrays.copyOf(anisou, 6);
            } else {
                System.arraycopy(anisou, 0, this.anisou, 0, 6);
            }
        }
    }

    /**
     * <p>
     * Setter for the field <code>anisouGradient</code>.</p>
     *
     * @param anisouGradient an array of double.
     */
    public void setAnisouGradient(double[] anisouGradient) {
        if (active) {
            if (anisouGradient == null) {
                this.anisouGradient = null;
            } else if (this.anisouGradient == null) {
                this.anisouGradient = Arrays.copyOf(anisouGradient, 6);
            } else {
                System.arraycopy(anisouGradient, 0, this.anisouGradient, 0, 6);
            }
        } else if (anisouGradient == null) {
            this.anisouGradient = null;
        } else {
            if (this.anisouGradient == null) {
                this.anisouGradient = new double[6];
            }
            Arrays.fill(anisouGradient, 0.0);
        }
    }

    /**
     * <p>
     * Setter for the field <code>anisouVelocity</code>.</p>
     *
     * @param anisouVelocity an array of double.
     */
    public void setAnisouVelocity(double[] anisouVelocity) {
        if (active) {
            if (anisouVelocity == null) {
                this.anisouVelocity = null;
            } else if (this.anisouVelocity == null) {
                this.anisouVelocity = Arrays.copyOf(anisouVelocity, 6);
            } else {
                System.arraycopy(anisouVelocity, 0, this.anisouVelocity, 0, 6);
            }
        } else if (anisouVelocity == null) {
            this.anisouVelocity = null;
        } else {
            if (this.anisouVelocity == null) {
                this.anisouVelocity = new double[6];
            }
            Arrays.fill(anisouVelocity, 0.0);
        }
    }

    /**
     * <p>
     * Setter for the field <code>anisouAcceleration</code>.</p>
     *
     * @param anisouAcceleration an array of double.
     */
    public void setAnisouAcceleration(double[] anisouAcceleration) {
        if (active) {
            if (anisouAcceleration == null) {
                this.anisouAcceleration = null;
            } else if (this.anisouAcceleration == null) {
                this.anisouAcceleration = Arrays.copyOf(anisouAcceleration, 6);
            } else {
                System.arraycopy(anisouAcceleration, 0, this.anisouAcceleration, 0, 6);
            }
        } else if (anisouAcceleration == null) {
            this.anisouAcceleration = null;
        } else {
            if (this.anisouAcceleration == null) {
                this.anisouAcceleration = new double[6];
            }
            Arrays.fill(anisouAcceleration, 0.0);
        }
    }

    /**
     * <p>
     * Setter for the field <code>anisouPreviousAcceleration</code>.</p>
     *
     * @param anisouPreviousAcceleration an array of double.
     */
    public void setAnisouPreviousAcceleration(double[] anisouPreviousAcceleration) {
        if (active) {
            if (anisouPreviousAcceleration == null) {
                this.anisouPreviousAcceleration = null;
            } else if (this.anisouPreviousAcceleration == null) {
                this.anisouPreviousAcceleration = Arrays.copyOf(anisouPreviousAcceleration, 6);
            } else {
                System.arraycopy(anisouPreviousAcceleration, 0, this.anisouPreviousAcceleration, 0, 6);
            }
        } else if (anisouPreviousAcceleration == null) {
            this.anisouPreviousAcceleration = null;
        } else {
            if (this.anisouPreviousAcceleration == null) {
                this.anisouPreviousAcceleration = new double[6];
            }
            Arrays.fill(anisouPreviousAcceleration, 0.0);
        }
    }

    /**
     * <p>
     * addToAnisouGradient</p>
     *
     * @param anisouGradient an array of double.
     */
    public void addToAnisouGradient(double[] anisouGradient) {
        if (active) {
            if (anisouGradient == null) {
            } else if (this.anisouGradient == null) {
                this.anisouGradient = Arrays.copyOf(anisouGradient, 6);
            } else {
                for (int i = 0; i < 6; i++) {
                    this.anisouGradient[i] += anisouGradient[i];
                }
            }
        }
    }

    /**
     * <p>
     * Setter for the field <code>formFactorWidth</code>.</p>
     *
     * @param width a double.
     */
    public void setFormFactorWidth(double width) {
        formFactorWidth = width;
        formFactorWidth2 = width * width;
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
     * Getter for the field <code>formFactorWidth</code>.</p>
     *
     * @return a double.
     */
    public double getFormFactorWidth2() {
        return formFactorWidth2;
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
        if (a != null && a.containsAtom(this)) {
            /**
             * for (Angle angle : angles) { if (angle == a) { return; } }
             */
            angles.add(a);
        }
    }

    /**
     * Specify that <b>this</b> Atom is part of a Bond
     *
     * @param b Bond that <b>this</b> Atom is part of
     */
    public void setBond(Bond b) {
        if (b != null && b.containsAtom(this)) {

            for (Bond bond : bonds) {
                if (bond == b) {
                    return;
                }
            }
            bonds.add(b);
        }
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
        if (a != null && !one_5s.contains(a)) {
            one_5s.add(a);
        }
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
        if (torsion == null || !torsion.containsAtom(this)) {
            return;
        }

        for (Torsion t : torsions) {
            if (torsion == t) {
                return;
            }
        }

        if (!isBackground) {
            for (Atom atom : torsion.getAtomArray()) {
                if (atom.isBackground()) {
                    return;
                }
            }
        }

        torsions.add(torsion);
        Atom a14 = torsion.get1_4(this);
        if (a14 != null) {
            // 1-5 atoms will be bonded to the 1-4 atom.
            if (a14.getBonds() != null) {
                for (Bond b : a14.getBonds()) {
                    Atom a15 = b.get1_2(a14);
                    // Do not include the 1-3 atom
                    if (!torsion.containsAtom(a15)) {
                        set1_5(a15);
                    }
                }
            }
        }
    }

    /**
     * Clear out the geometry lists for this atom.
     */
    public void clearGeometry() {
        bonds.clear();
        angles.clear();
        torsions.clear();
        one_5s.clear();
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
        if (active) {
            xyzGradient[0] = x;
            xyzGradient[1] = y;
            xyzGradient[2] = z;
        }
    }

    /**
     * <p>
     * setLambdaXYZGradient</p>
     *
     * @param x a double.
     * @param y a double.
     * @param z a double.
     */
    public void setLambdaXYZGradient(double x, double y, double z) {
        if (active) {
            xyzLambdaGradient[0] = x;
            xyzLambdaGradient[1] = y;
            xyzLambdaGradient[2] = z;
        }
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
        if (active) {
            xyzGradient[0] += x;
            xyzGradient[1] += y;
            xyzGradient[2] += z;
        }
    }

    public void addToXYZGradient(int axis, double value) {
        if (active) {
            xyzGradient[axis] += value;
        }
    }

    /**
     * <p>
     * addToLambdaXYZGradient</p>
     *
     * @param x a double.
     * @param y a double.
     * @param z a double.
     */
    public void addToLambdaXYZGradient(double x, double y, double z) {
        if (active) {
            xyzLambdaGradient[0] += x;
            xyzLambdaGradient[1] += y;
            xyzLambdaGradient[2] += z;
        }
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
     * getLambdaXYZGradient</p>
     *
     * @param x an array of double.
     */
    public void getLambdaXYZGradient(double x[]) {
        if (x == null) {
            x = new double[3];
        }
        x[0] = xyzLambdaGradient[0];
        x[1] = xyzLambdaGradient[1];
        x[2] = xyzLambdaGradient[2];
    }

    /**
     * @param velocity an array of double.
     */
    public void getVelocity(double velocity[]) {
        if (velocity == null) {
            velocity = new double[3];
        }
        velocity[0] = this.velocity[0];
        velocity[1] = this.velocity[1];
        velocity[2] = this.velocity[2];
    }

    /**
     * @param acceleration an array of double.
     */
    public void getAcceleration(double acceleration[]) {
        if (acceleration == null) {
            acceleration = new double[3];
        }
        acceleration[0] = this.acceleration[0];
        acceleration[1] = this.acceleration[1];
        acceleration[2] = this.acceleration[2];
    }

    /**
     * @param previousAcceleration an array of double.
     */
    public void getPreviousAcceleration(double previousAcceleration[]) {
        if (previousAcceleration == null) {
            previousAcceleration = new double[3];
        }
        previousAcceleration[0] = this.previousAcceleration[0];
        previousAcceleration[1] = this.previousAcceleration[1];
        previousAcceleration[2] = this.previousAcceleration[2];
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

    public void setMultipoleType(MultipoleType multipoleType) {
        this.multipoleType = multipoleType;
    }

    public void setAxisAtoms(Atom... set) {
        axisAtoms = (Atom[]) ArrayUtils.clone(set);
        if (axisAtoms == null) {
            axisAtomIndices = null;
            return;
        }
        axisAtomIndices = new int[axisAtoms.length];
        for (int i = 0; i < set.length; i++) {
            axisAtomIndices[i] = set[i].getArrayIndex();
        }
    }

    public void setAxisAtomIndices(int[] set) {
        axisAtomIndices = ArrayUtils.clone(set);
    }

    public int[] getAxisAtomIndices() {
        return ArrayUtils.clone(axisAtomIndices);
    }

    public void setFrameDefinition(MultipoleFrameDefinition frame) {
        this.frameDefinition = frame;
    }

    public void setMoleculeNumber(int molecule) {
        this.moleculeNumber = molecule;
    }

    public int getMoleculeNumber() {
        return moleculeNumber;
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
    public void setSelected(boolean selected) {
        this.selected = selected;
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

    public final int getXyzIndex() {
        return xyzIndex;
    }

    public final void setXyzIndex(int set) {
        xyzIndex = set;
    }

    /**
     * Note: the MolecularAssembly to which this Atom belongs is cached. If you
     * wanna move atoms between assemblies, un-cache it.
     */
    public final int getIndex() {
        switch (MolecularAssembly.atomIndexing) {
            case PERSIST:
                return persistentIndex;
            default:
            case XYZ:
                return xyzIndex;
        }
    }

    public final int getArrayIndex() {
        return xyzIndex - 1;
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
            multipoleBuffer.append(String.format("\n%11$s % 7.5f\n" + "%11$s % 7.5f % 7.5f % 7.5f\n" + "%11$s % 7.5f\n" + "%11$s % 7.5f % 7.5f\n" + "%11$s % 7.5f % 7.5f % 7.5f", multipoleType.getCharge(), globalDipole[0], globalDipole[1],
                    globalDipole[2], globalQuadrupole[0][0],
                    globalQuadrupole[1][0], globalQuadrupole[1][1],
                    globalQuadrupole[2][0], globalQuadrupole[2][1],
                    globalQuadrupole[2][2], "                 "));
            return multipoleBuffer.toString();
        }
    }

    public enum Descriptions {
        Default, Trim, XyzIndex_Name, ArrayIndex_Name, Resnum_Name;
    }

    public String describe(Descriptions type) {
        switch (type) {
            default:
            case Default:
                return toString();
            case Trim:
                return format("%d-%-3s %s %s%d", getIndex(), getName(), resName, segID, resSeq);
            case XyzIndex_Name:
                return format("%d-%s", getIndex(), getName());
            case ArrayIndex_Name:
                return format("%d-%s", getIndex() - 1, getName());
            case Resnum_Name:
                return format("%d-%s", resSeq, getName());
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        if (altLoc != null && altLoc != ' ') {
            return String.format("%s %7d-%s %s %d (%7.2f,%7.2f,%7.2f) %s", altLoc, getIndex(), getName(),
                    resName, resSeq, xyz[0], xyz[1], xyz[2], segID);
        }
        if (resName == null) {
            return String.format("%7d-%s (%7.2f,%7.2f,%7.2f)", getIndex(), getName(), xyz[0], xyz[1], xyz[2]);
        }
        return String.format("%7d-%s %s %d (%7.2f,%7.2f,%7.2f) %s", getIndex(), getName(), resName, resSeq,
                xyz[0], xyz[1], xyz[2], segID);
    }

    public String toShortString() {
        return this.describe(Descriptions.XyzIndex_Name);
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
