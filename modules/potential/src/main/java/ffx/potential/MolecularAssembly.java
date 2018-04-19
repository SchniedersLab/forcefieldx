/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.potential;

import javax.media.j3d.Appearance;
import javax.media.j3d.BoundingSphere;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.ColoringAttributes;
import javax.media.j3d.GeometryArray;
import javax.media.j3d.Group;
import javax.media.j3d.LineArray;
import javax.media.j3d.LineAttributes;
import javax.media.j3d.Link;
import javax.media.j3d.Material;
import javax.media.j3d.Node;
import javax.media.j3d.RenderingAttributes;
import javax.media.j3d.Shape3D;
import javax.media.j3d.SharedGroup;
import javax.media.j3d.Switch;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;
import javax.vecmath.Color3f;
import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.List;
import java.util.ListIterator;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import static java.lang.String.format;

import com.sun.j3d.utils.picking.PickTool;

import org.apache.commons.configuration.CompositeConfiguration;
import org.jdesktop.j3d.loaders.vrml97.VrmlLoader;
import org.jdesktop.j3d.loaders.vrml97.VrmlScene;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.numerics.VectorMath;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Indexing;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSGroup;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.RendererCache;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Residue.ResiduePosition;
import ffx.potential.parameters.ForceField;
import static ffx.potential.bonded.Residue.ResiduePosition.FIRST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.LAST_RESIDUE;
import static ffx.potential.bonded.Residue.ResiduePosition.MIDDLE_RESIDUE;
import static ffx.potential.extended.ExtUtils.prop;

/**
 * The MolecularAssembly class is a collection of Polymers, Hetero Molecules,
 * Ions and Water
 *
 * @author Michael J. Schnieders
 */
public class MolecularAssembly extends MSGroup {

    private static final Logger logger = Logger.getLogger(MolecularAssembly.class.getName());
    private static final long serialVersionUID = 1L;
    /**
     * Constant <code>MultiScaleLevel=4</code>
     */
    public static final int MultiScaleLevel = 4;
    /**
     * Constant <code>KCAL_TO_KJ=4.184</code>
     */
    public static final double KCAL_TO_KJ = 4.184;
    private static double[] a = new double[3];

    public static final Indexing atomIndexing = prop(Indexing.class, "sys.atomIndexing", Indexing.XYZ);
    /**
     * Persistent index parallel to xyzIndex.
     */
    public static int persistentAtomIndexer = 1;

    // MolecularAssembly member variables
    private File file;
    private final List<String> headerLines = new ArrayList<>();
    protected ForceField forceField;
    private ForceFieldEnergy potentialEnergy;
    private CompositeConfiguration properties;
    private Vector3d offset;
    private int cycles = 1;
    private int currentCycle = 1;
    private List<String> altLoc = null;
    // Data Nodes
    private final MSNode ions = new MSNode("Ions");
    private final MSNode water = new MSNode("Waters");
    private final MSNode molecules = new MSNode("Hetero Molecules");
    // 3D Graphics Nodes - There is a diagram explaining the MolecularAssembly
    // Scenegraph below
    private BranchGroup branchGroup;
    private TransformGroup originToRot;
    private Transform3D originToRotT3D;
    private Vector3d originToRotV3D;
    private TransformGroup rotToCOM;
    private Transform3D rotToCOMT3D;
    private Vector3d rotToCOMV3D;
    private BranchGroup base;
    private Switch switchGroup;
    private Shape3D wire;
    private BranchGroup vrml;
    private TransformGroup vrmlTG;
    private Transform3D vrmlTd;
    private BranchGroup childNodes;
    private Atom[] atomLookUp;
    private LineAttributes lineAttributes;
    private File vrmlFile = null;
    private URL vrmlURL = null;
    private boolean visible = false;
    private final ArrayList<BranchGroup> myNewShapes = new ArrayList<>();

    public enum FractionalMode {OFF, MOLECULE, ATOM};
    private FractionalMode fractionalMode = FractionalMode.MOLECULE;
    private double fractionalCoordinates[][];

    /**
     * <p>
     * Constructor for MolecularAssembly.</p>
     *
     * @param name a {@link java.lang.String} object.
     */
    public MolecularAssembly(String name) {
        super(name);
        MolecularAssembly.this.getAtomNode().setName("MacroMolecules");
        add(molecules);
        add(ions);
        add(water);
    }

    /**
     * <p>
     * Constructor for MolecularAssembly.</p>
     *
     * @param name     a {@link java.lang.String} object.
     * @param Polymers a {@link ffx.potential.bonded.MSNode} object.
     */
    public MolecularAssembly(String name, MSNode Polymers) {
        super(name, Polymers);
    }

    public MolecularAssembly(String name, MSNode Polymers, CompositeConfiguration properties) {
        this(name, Polymers);
        this.properties = properties;
    }

    public void setFractionalMode(FractionalMode mode) {
        fractionalMode = mode;
    }

    public FractionalMode getFractionalMode() {
        return fractionalMode;
    }

    /**
     * Count the number of fractional coordinate entities in the system. If fractionalMode is MOLECULE, then the count
     * is equal to the number of molecules. If fractionalMode is ATOM, then the count is the number of atoms. Otherwise
     * the count is zero.
     *
     * @return The number of fractional coordinate entities.
     */
    public int fractionalCount() {
        int count = 0;
        switch (fractionalMode) {
            case MOLECULE:
                // Move polymers togethers.
                Polymer polymers[] = getChains();
                if (polymers != null && polymers.length > 0) {
                    count += polymers.length;
                }
                List<Molecule> molecules = getMolecules();
                if (molecules != null) {
                    count += molecules.size();
                }

                List<MSNode> waters = getWaters();
                if (waters != null) {
                    count += waters.size();
                }

                List<MSNode> ions = getIons();
                if (ions != null) {
                    count += ions.size();
                }
                break;
            case ATOM:
                count = getAtomArray().length;
                break;
            case OFF:
                count = 0;
                break;
        }

        if (fractionalCoordinates == null || fractionalCoordinates.length != count) {
            fractionalCoordinates = new double[count][3];
        }

        return count;
    }

    /**
     * Compute fractional coordinates.
     */
    public void computeFractionalCoordinates() {

        // Count up the number of fractional coordinate entities.
        fractionalCount();

        Crystal unitCell = getCrystal().getUnitCell();
        double[] com = new double[3];

        switch (fractionalMode) {
            case MOLECULE:
                int iMolecule = 0;
                Polymer polymers[] = getChains();
                if (polymers != null && polymers.length > 0) {
                    // Find the center of mass
                    for (Polymer polymer : polymers) {
                        List<Atom> list = polymer.getAtomList();
                        com[0] = 0.0;
                        com[1] = 0.0;
                        com[2] = 0.0;
                        double totalMass = 0.0;
                        for (Atom atom : list) {
                            double m = atom.getMass();
                            com[0] += atom.getX() * m;
                            com[1] += atom.getY() * m;
                            com[2] += atom.getZ() * m;
                            totalMass += m;
                        }
                        com[0] /= totalMass;
                        com[1] /= totalMass;
                        com[2] /= totalMass;
                        unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
                    }
                }

                // Loop over each molecule
                List<Molecule> molecules = getMolecules();
                for (MSNode molecule : molecules) {
                    List<Atom> list = molecule.getAtomList();
                    // Find the center of mass
                    com[0] = 0.0;
                    com[1] = 0.0;
                    com[2] = 0.0;
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        com[0] += atom.getX() * m;
                        com[1] += atom.getY() * m;
                        com[2] += atom.getZ() * m;
                        totalMass += m;
                    }
                    com[0] /= totalMass;
                    com[1] /= totalMass;
                    com[2] /= totalMass;
                    unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
                }

                // Loop over each water
                List<MSNode> waters = getWaters();
                for (MSNode water : waters) {
                    List<Atom> list = water.getAtomList();
                    // Find the center of mass
                    com[0] = 0.0;
                    com[1] = 0.0;
                    com[2] = 0.0;
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        com[0] += atom.getX() * m;
                        com[1] += atom.getY() * m;
                        com[2] += atom.getZ() * m;
                        totalMass += m;
                    }
                    com[0] /= totalMass;
                    com[1] /= totalMass;
                    com[2] /= totalMass;
                    unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
                }

                // Loop over each ion
                List<MSNode> ions = getIons();
                for (MSNode ion : ions) {
                    List<Atom> list = ion.getAtomList();
                    // Find the center of mass
                    com[0] = 0.0;
                    com[1] = 0.0;
                    com[2] = 0.0;
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        com[0] += atom.getX() * m;
                        com[1] += atom.getY() * m;
                        com[2] += atom.getZ() * m;
                        totalMass += m;
                    }
                    com[0] /= totalMass;
                    com[1] /= totalMass;
                    com[2] /= totalMass;
                    unitCell.toFractionalCoordinates(com, fractionalCoordinates[iMolecule++]);
                }
                break;
            case ATOM:
                Atom atoms[] = getAtomArray();
                int nAtoms = atoms.length;
                for (int i = 0; i < nAtoms; i++) {
                    atoms[i].getXYZ(com);
                    unitCell.toFractionalCoordinates(com, fractionalCoordinates[i]);
                }
                break;
            case OFF:
                break;
        }
    }

    /**
     * Move to fractional coordinates.
     */
    public void moveToFractionalCoordinates() {

        if (fractionalCoordinates == null) {
            return;
        }

        Crystal unitCell = getCrystal().getUnitCell();
        double[] com = new double[3];

        switch (fractionalMode) {
            case MOLECULE:
                int iMolecule = 0;
                Polymer polymers[] = getChains();
                if (polymers != null && polymers.length > 0) {
                    // Find the center of mass
                    for (Polymer polymer : polymers) {
                        List<Atom> list = polymer.getAtomList();
                        double totalMass = 0.9;
                        com[0] = 0.0;
                        com[1] = 0.0;
                        com[2] = 0.0;
                        for (Atom atom : list) {
                            double m = atom.getMass();
                            com[0] += atom.getX() * m;
                            com[1] += atom.getY() * m;
                            com[2] += atom.getZ() * m;
                            totalMass += m;
                        }
                        com[0] /= totalMass;
                        com[1] /= totalMass;
                        com[2] /= totalMass;
                        // Find the new center of mass in fractional coordinates.
                        unitCell.toFractionalCoordinates(com, com);
                        // Find the reciprocal translation vector.
                        double[] frac = fractionalCoordinates[iMolecule++];
                        com[0] = frac[0] - com[0];
                        com[1] = frac[1] - com[1];
                        com[2] = frac[2] - com[2];
                        // Convert the fractional translation vector to Cartesian coordinates.
                        unitCell.toCartesianCoordinates(com, com);
                        // Move all atoms.
                        for (Atom atom : list) {
                            atom.move(com);
                        }
                    }
                }

                // Loop over each molecule
                List<Molecule> molecules = getMolecules();
                for (MSNode molecule : molecules) {
                    List<Atom> list = molecule.getAtomList();
                    // Find the center of mass
                    com[0] = 0.0;
                    com[1] = 0.0;
                    com[2] = 0.0;
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        com[0] += atom.getX() * m;
                        com[1] += atom.getY() * m;
                        com[2] += atom.getZ() * m;
                        totalMass += m;
                    }
                    com[0] /= totalMass;
                    com[1] /= totalMass;
                    com[2] /= totalMass;
                    // Find the new center of mass in fractional coordinates.
                    unitCell.toFractionalCoordinates(com, com);
                    // Find the reciprocal translation vector to the previous COM.
                    double[] frac = fractionalCoordinates[iMolecule++];
                    com[0] = frac[0] - com[0];
                    com[1] = frac[1] - com[1];
                    com[2] = frac[2] - com[2];
                    // Convert the fractional translation vector to Cartesian coordinates.
                    unitCell.toCartesianCoordinates(com, com);
                    // Move all atoms.
                    for (Atom atom : list) {
                        atom.move(com);
                    }
                }

                // Loop over each water
                List<MSNode> waters = getWaters();
                for (MSNode water : waters) {
                    List<Atom> list = water.getAtomList();
                    // Find the center of mass
                    com[0] = 0.0;
                    com[1] = 0.0;
                    com[2] = 0.0;
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        com[0] += atom.getX() * m;
                        com[1] += atom.getY() * m;
                        com[2] += atom.getZ() * m;
                        totalMass += m;
                    }
                    com[0] /= totalMass;
                    com[1] /= totalMass;
                    com[2] /= totalMass;
                    // Find the new center of mass in fractional coordinates.
                    unitCell.toFractionalCoordinates(com, com);
                    // Find the reciprocal translation vector to the previous COM.
                    double[] frac = fractionalCoordinates[iMolecule++];
                    com[0] = frac[0] - com[0];
                    com[1] = frac[1] - com[1];
                    com[2] = frac[2] - com[2];
                    // Convert the fractional translation vector to Cartesian coordinates.
                    unitCell.toCartesianCoordinates(com, com);

                    double r = ffx.numerics.VectorMath.r(com);
                    /**
                     * Warn if an atom is moved more than 1 Angstrom.
                     */
                    if (r > 1.0) {
                        int i = iMolecule - 1;
                        logger.info(String.format(" %d R: %16.8f", i, r));
                        logger.info(String.format(" %d FRAC %16.8f %16.8f %16.8f", i, frac[0], frac[1], frac[2]));
                        logger.info(String.format(" %d COM  %16.8f %16.8f %16.8f", i, com[0], com[1], com[2]));
                    }

                    // Move all atoms.
                    for (Atom atom : list) {
                        atom.move(com);
                    }
                }

                // Loop over each ion
                List<MSNode> ions = getIons();
                for (MSNode ion : ions) {
                    List<Atom> list = ion.getAtomList();
                    // Find the center of mass
                    com[0] = 0.0;
                    com[1] = 0.0;
                    com[2] = 0.0;
                    double totalMass = 0.0;
                    for (Atom atom : list) {
                        double m = atom.getMass();
                        com[0] += atom.getX() * m;
                        com[1] += atom.getY() * m;
                        com[2] += atom.getZ() * m;
                        totalMass += m;
                    }
                    com[0] /= totalMass;
                    com[1] /= totalMass;
                    com[2] /= totalMass;
                    // Find the new center of mass in fractional coordinates.
                    unitCell.toFractionalCoordinates(com, com);
                    // Find the reciprocal translation vector to the previous COM.
                    double[] frac = fractionalCoordinates[iMolecule++];
                    com[0] = frac[0] - com[0];
                    com[1] = frac[1] - com[1];
                    com[2] = frac[2] - com[2];
                    // Convert the fractional translation vector to Cartesian coordinates.
                    unitCell.toCartesianCoordinates(com, com);
                    // Move all atoms.
                    for (Atom atom : list) {
                        atom.move(com);
                    }
                }
                break;
            case ATOM:
                Atom atoms[] = getAtomArray();
                int nAtoms = atoms.length;
                for (int i = 0; i < nAtoms; i++) {
                    // Convert the stored factional coordinates to Cartesian coordinates in the current unitcell.
                    unitCell.toCartesianCoordinates(fractionalCoordinates[i], com);
                    atoms[i].moveTo(com);
                }
                break;
            case OFF:
                break;
        }
    }

    /**
     * <p>
     * Setter for the field <code>forceField</code>.</p>
     *
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     */
    public void setForceField(ForceField forceField) {
        this.forceField = forceField;
    }

    public void setPropertiesFromForceField() {
        this.properties = forceField.getProperties();
    }

    /**
     * <p>
     * setPotential</p>
     *
     * @param potentialEnergy a {@link ffx.potential.ForceFieldEnergy} object.
     */
    public void setPotential(ForceFieldEnergy potentialEnergy) {
        this.potentialEnergy = potentialEnergy;
    }

    /**
     * <p>
     * Getter for the field <code>potentialEnergy</code>.</p>
     *
     * @return a {@link ffx.potential.ForceFieldEnergy} object.
     */
    public ForceFieldEnergy getPotentialEnergy() {
        return potentialEnergy;
    }

    public ParallelTeam getParallelTeam() {
        if (potentialEnergy != null) {
            return potentialEnergy.getParallelTeam();
        } else {
            return null;
        }
    }

    public ResiduePosition getResiduePosition(int residueNumber) {
        ResiduePosition position;
        int numberOfResidues = 0;
        Polymer polymers[] = getChains();
        int nPolymers = polymers.length;
        for (int i = 0; i < nPolymers; i++) {
            Polymer polymer = polymers[i];
            ArrayList<Residue> residues = polymer.getResidues();
            numberOfResidues += residues.size();
        }
        if (residueNumber == 0) {
            position = FIRST_RESIDUE;
        } else if (residueNumber == numberOfResidues - 1) {
            position = LAST_RESIDUE;
        } else {
            position = MIDDLE_RESIDUE;
        }
        return position;
    }

    /**
     * Adds a header line to this MolecularAssembly (particularly for PDB
     * formats)
     *
     * @param line Line to add.
     */
    public void addHeaderLine(String line) {
        headerLines.add(line);
    }

    /**
     * Gets the header lines associated with this MolecularAssembly
     * (particularly for PDB)
     *
     * @return Header lines.
     */
    public String[] getHeaderLines() {
        String[] ret = new String[headerLines.size()];
        headerLines.toArray(ret);
        return ret;
    }

    /**
     * <p>
     * getCrystal</p>
     *
     * @return a {@link ffx.crystal.Crystal} object.
     */
    public Crystal getCrystal() {
        if (potentialEnergy == null) {
            return null;
        }
        return potentialEnergy.getCrystal();
    }

    /**
     * <p>
     * Getter for the field <code>forceField</code>.</p>
     *
     * @return a {@link ffx.potential.parameters.ForceField} object.
     */
    public ForceField getForceField() {
        return forceField;
    }

    /**
     * <p>
     * addAltLocation</p>
     *
     * @param s a {@link java.lang.String} object.
     */
    public void addAltLocation(String s) {
        if (altLoc == null) {
            altLoc = new Vector<>();
        }
        if (!altLoc.contains(s)) {
            altLoc.add(s);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public MSNode addMSNode(MSNode o) {
        ArrayList<MSNode> Polymers = getAtomNodeList();
        if (o instanceof Atom) {
            Atom atom = (Atom) o;
            if (atom.isModRes()) {
                return getResidue(atom, true, Residue.ResidueType.AA);
            } else if (!atom.isHetero()) {
                return getResidue(atom, true);
            } else {
                return getMolecule(atom, true);
            }
        } else if (o instanceof Residue) {
            Residue residue = (Residue) o;
            Character chainID = residue.getChainID();
            String segID = residue.getSegID();
            int index = Polymers.indexOf(new Polymer(chainID, segID));
            /**
             * See if the polymer already exists.
             */
            if (index != -1) {
                Polymer c = (Polymer) Polymers.get(index);
                setFinalized(false);
                return c.addMSNode(residue);
            } else {
                Polymer newc = new Polymer(chainID, segID);
                getAtomNode().add(newc);
                setFinalized(false);
                return newc.addMSNode(residue);
            }
        } else if (o instanceof Polymer) {
            Polymer c = (Polymer) o;
            int index = Polymers.indexOf(c);
            if (index == -1) {
                getAtomNode().add(c);
                setFinalized(false);
                return c;
            } else {
                return (Polymer) Polymers.get(index);
            }
        } else if (o instanceof Molecule) {
            Molecule m = (Molecule) o;
            if (m.getAtomNode().getChildCount() == 1) {
                ions.add(m);
                return m;
            } else if (Utilities.isWaterOxygen((Atom) m.getAtomNode().getChildAt(0))) {
                water.add(m);
                return m;
            } else {
                molecules.add(m);
                return m;
            }
        } else {
            String message = "Programming error in MolecularAssembly addNode";
            logger.log(Level.SEVERE, message);
            return o;
        }
    }

    /**
     * <p>
     * center</p>
     */
    public void center() {
        double center[] = getMultiScaleCenter(false);
        offset = new Vector3d(center);
        if (vrml != null) {
            vrmlTd.set(offset);
            vrmlTG.setTransform(vrmlTd);
        }
        offset.negate();
        originToRotV3D.set(offset);
        originToRotT3D.setTranslation(originToRotV3D);
        originToRot.setTransform(originToRotT3D);
        rotToCOMT3D.setIdentity();
        rotToCOM.setTransform(rotToCOMT3D);
        offset.negate();
        rotateAbout(offset);
        originToRotT3D.get(offset);
    }

    /**
     * <p>
     * centerAt</p>
     *
     * @param d an array of double.
     */
    public void centerAt(double[] d) {
        double[] Rc = {0, 0, 0};
        double[] c = new double[3];
        ListIterator<Atom> li;
        int i, num = getAtomList().size();
        for (li = getAtomList().listIterator(); li.hasNext(); ) {
            (li.next()).getXYZ(a);
            Rc[0] += a[0];
            Rc[1] += a[1];
            Rc[2] += a[2];
        }
        for (i = 0; i < 3; i++) {
            Rc[i] /= num;
        }
        VectorMath.diff(d, Rc, c);
        for (li = getAtomList().listIterator(); li.hasNext(); ) {
            (li.next()).move(c);
        }
    }

    /**
     * <p>
     * centerView</p>
     *
     * @param rot   a boolean.
     * @param trans a boolean.
     */
    public void centerView(boolean rot, boolean trans) {
        originToRot.getTransform(originToRotT3D);
        if (rot) {
            Matrix3d m3d = new Matrix3d();
            m3d.setIdentity();
            originToRotT3D.setRotation(m3d);
            // rotToCOMT3D.setRotation(m3d);
        }
        if (trans) {
            originToRotV3D.set(offset);
            originToRotT3D.set(originToRotV3D);
        }
        originToRot.setTransform(originToRotT3D);
        // rotToCOM.setTransform(rotToCOMT3D);
    }

    /**
     * <p>
     * createBox</p>
     */
    public void createBox() {
        int vertices = 8;
        LineArray la = new LineArray(4 * vertices, GeometryArray.COORDINATES | GeometryArray.COLOR_4 | GeometryArray.NORMALS);
        la.setCapability(LineArray.ALLOW_COORDINATE_WRITE);
        la.setCapability(LineArray.ALLOW_COORDINATE_READ);
        la.setCapability(LineArray.ALLOW_COLOR_WRITE);
        la.setCapability(LineArray.ALLOW_COUNT_READ);
        la.setCapability(LineArray.ALLOW_INTERSECT);
        la.setCapability(LineArray.ALLOW_FORMAT_READ);
        // Create a normal
        // for (ListIterator<MSNode> li = bondlist.listIterator(); li.hasNext(); ){
        // la.setCoordinate(i, a1);
        // la.setColor(i, col);
        // la.setNormal(i++, a1);
        // }
        ColoringAttributes cola = new ColoringAttributes(new Color3f(),
                ColoringAttributes.SHADE_GOURAUD);
        Appearance app = new Appearance();
        lineAttributes = new LineAttributes();
        lineAttributes.setLineWidth(RendererCache.bondwidth);
        lineAttributes.setCapability(LineAttributes.ALLOW_WIDTH_WRITE);
        lineAttributes.setLineAntialiasingEnable(true);
        app.setLineAttributes(lineAttributes);
        app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_READ);
        app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_WRITE);
        RenderingAttributes ra = new RenderingAttributes();
        ra.setAlphaTestValue(0.1f);
        ra.setAlphaTestFunction(RenderingAttributes.GREATER);
        ra.setDepthBufferEnable(true);
        ra.setDepthBufferWriteEnable(true);
        app.setRenderingAttributes(ra);
        app.setColoringAttributes(cola);
        Shape3D wireframe = new Shape3D(la, app);
        // PickTool.setCapabilities(wire, PickTool.INTERSECT_COORD);
        wireframe.setUserData(this);
        wireframe.setBounds(new BoundingSphere(new Point3d(0, 0, 0), 10.0));
        try {
            wireframe.setBoundsAutoCompute(false);
        } catch (Exception e) {
            e.printStackTrace();
        }
        wireframe.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        wireframe.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
        // return wire;
    }

    /**
     * The MolecularAssembly BranchGroup has two TransformGroups between it and
     * the "base" node where geometry is attached. If the point between the two
     * transformations is where user rotation occurs. For example, if rotating
     * about the center of mass of the system, the RotToCOM transformation will
     * be an identity transformation (ie. none). If rotation is about some atom
     * or group of atoms within the system, then the RotToCOM transformation
     * will be a translation from that point to the COM.
     *
     * @param zero boolean
     * @return BranchGroup
     */
    public BranchGroup createScene(boolean zero) {
        originToRotT3D = new Transform3D();
        originToRotV3D = new Vector3d();
        originToRot = new TransformGroup(originToRotT3D);
        branchGroup = new BranchGroup();
        rotToCOM = new TransformGroup();
        rotToCOMT3D = new Transform3D();
        rotToCOMV3D = new Vector3d();
        // Set capabilities needed for picking and moving the MolecularAssembly
        branchGroup.setCapability(BranchGroup.ALLOW_DETACH);
        originToRot.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
        originToRot.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        originToRot.setCapability(TransformGroup.ENABLE_PICK_REPORTING);
        rotToCOM.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        rotToCOM.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
        // Put the MolecularAssembly in the middle of the scene
        if (zero) {
            originToRotV3D.set(0.0, 0.0, 0.0);
            originToRotT3D.set(originToRotV3D);
            originToRot.setTransform(originToRotT3D);
        }
        wire = renderWire();
        switchGroup = new Switch(Switch.CHILD_NONE);
        switchGroup.setCapability(Switch.ALLOW_SWITCH_WRITE);
        base = new BranchGroup();
        base.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
        base.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
        childNodes = new BranchGroup();
        childNodes.setCapability(BranchGroup.ALLOW_DETACH);
        childNodes.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
        childNodes.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
        switchGroup.addChild(base);
        if (wire != null) {
            base.addChild(wire);
        }
        vrml = loadVRML();
        if (vrml != null) {
            vrmlTG = new TransformGroup();
            vrmlTd = new Transform3D();
            vrmlTG.setTransform(vrmlTd);
            vrmlTG.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
            vrmlTG.addChild(vrml);
            switchGroup.addChild(vrmlTG);
            setView(RendererCache.ViewModel.INVISIBLE, null);
        }
        switchGroup.setWhichChild(Switch.CHILD_ALL);
        rotToCOM.addChild(switchGroup);
        originToRot.addChild(rotToCOM);
        branchGroup.addChild(originToRot);
        branchGroup.compile();
        return branchGroup;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean destroy() {
        try {
            potentialEnergy.destroy();
        } catch (Exception ex) {
            logger.warning(String.format(" Exception in shutting down force field "
                    + "energy parallel teams: %s", ex.toString()));
        }
        detach();
        return super.destroy();
    }

    /**
     * <p>
     * detach</p>
     */
    public void detach() {
        synchronized (this) {
            if (branchGroup != null && branchGroup.isLive()) {
                branchGroup.detach();
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void finalize(boolean finalizeGroups, ForceField forceField) {
        setFinalized(false);
        if (finalizeGroups) {
            bondTime = 0;
            angleTime = 0;
            stretchBendTime = 0;
            ureyBradleyTime = 0;
            outOfPlaneBendTime = 0;
            torsionTime = 0;
            piOrbitalTorsionTime = 0;
            torsionTorsionTime = 0;
            ArrayList<MSNode> Polymers = getAtomNodeList();
            for (ListIterator<MSNode> li = Polymers.listIterator(); li.hasNext(); ) {
                MSGroup group = (MSGroup) li.next();
                if (logger.isLoggable(Level.FINE)) {
                    logger.fine(" Finalizing bonded terms for polymer " + group.toString());
                }
                try {
                    group.finalize(true, forceField);
                } catch (Exception e) {
                    String message = "Fatal exception finalizing " + group.toString();
                    logger.log(Level.SEVERE, message, e);
                    System.exit(-1);
                }
                if (logger.isLoggable(Level.FINE)) {
                    Runtime runtime = Runtime.getRuntime();
                    long occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
                    long MB = 1024 * 1024;
                    logger.fine("\n In-Use Memory   (Mb): " + occupiedMemory / MB
                            + "\n Free Memory     (Mb): " + runtime.freeMemory() / MB
                            + "\n Total Memory    (Mb): " + runtime.totalMemory() / MB);
                }
            }
            for (MSNode m : molecules.getChildList()) {
                Molecule molecule = (Molecule) m;
                molecule.finalize(true, forceField);
            }
            for (MSNode m : water.getChildList()) {
                Molecule molecule = (Molecule) m;
                molecule.finalize(true, forceField);
            }
            for (MSNode m : ions.getChildList()) {
                Molecule molecule = (Molecule) m;
                molecule.finalize(true, forceField);
            }
            if (logger.isLoggable(Level.FINE)) {
                StringBuilder sb = new StringBuilder("\n Time to create bonded energy terms\n\n");
                sb.append(String.format(" Bond Streching     %10.3f\n", bondTime * 1.0e-9));
                sb.append(String.format(" Angle Bending      %10.3f\n", angleTime * 1.0e-9));
                sb.append(String.format(" Stretch-Bend       %10.3f\n", stretchBendTime * 1.0e-9));
                sb.append(String.format(" Urey-Bradley       %10.3f\n", ureyBradleyTime * 1.0e-9));
                sb.append(String.format(" Out-of-Plane Bend  %10.3f\n", outOfPlaneBendTime * 1.0e-9));
                sb.append(String.format(" Torsionanl Angle   %10.3f\n", torsionTime * 1.0e-9));
                sb.append(String.format(" Pi-Orbital Torsion %10.3f\n", piOrbitalTorsionTime * 1.0e-9));
                sb.append(String.format(" Torsion-Torsion    %10.3f\n", torsionTorsionTime * 1.0e-9));
                logger.fine(sb.toString());
            }
        }
        if (!java.awt.GraphicsEnvironment.isHeadless()) {
            createScene(!finalizeGroups);
            center();
        }
        removeLeaves();
        setFinalized(true);
    }

    /**
     * <p>
     * findAtom</p>
     *
     * @param atom a {@link ffx.potential.bonded.Atom} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom findAtom(Atom atom) {
        if (!atom.isHetero() || atom.isModRes()) {
            Polymer polymer = getPolymer(atom.getChainID(), atom.getSegID(), false);
            if (polymer != null) {
                Residue res = polymer.getResidue(atom.getResidueName(), atom.getResidueNumber(), false);
                if (res != null) {
                    MSNode node = res.getAtomNode();
                    Atom root = (Atom) node.contains(atom);
                    return root;
                }
            }
            return null;
        } else {
            ArrayList<Molecule> list = getMolecules();
            for (MSNode node : list) {
                Molecule m = (Molecule) node;
                if (m.getSegID().equalsIgnoreCase(atom.getSegID())
                        && m.getResidueName().equalsIgnoreCase(atom.getResidueName())
                        && m.getResidueNumber() == atom.getResidueNumber()) {
                    Atom root = (Atom) node.contains(atom);
                    return root;
                }
            }
            ArrayList<MSNode> ionList = getIons();
            for (MSNode node : ionList) {
                Molecule m = (Molecule) node;
                if (m.getSegID().equalsIgnoreCase(atom.getSegID())
                        && m.getResidueName().equalsIgnoreCase(atom.getResidueName())
                        && m.getResidueNumber() == atom.getResidueNumber()) {
                    Atom root = (Atom) node.contains(atom);
                    return root;
                }
            }
            ArrayList<MSNode> waterList = getWaters();
            for (MSNode node : waterList) {
                Molecule m = (Molecule) node;
                if (m.getSegID().equalsIgnoreCase(atom.getSegID())
                        && m.getResidueName().equalsIgnoreCase(atom.getResidueName())
                        && m.getResidueNumber() == atom.getResidueNumber()) {
                    Atom root = (Atom) node.contains(atom);
                    return root;
                }
            }
            return null;
        }
    }

    /**
     * This method assigns a unique integer to every molecule in the
     * MolecularAssembly beginning at 0. An integer array with these values for
     * each atom is returned.
     *
     * @return an array of molecule numbers for each atom.
     */
    public int[] getMoleculeNumbers() {
        int moleculeNumber[] = new int[getAtomList().size()];
        int current = 0;
        // Loop over polymers together
        Polymer[] polymers = getChains();
        if (polymers != null && polymers.length > 0) {
            for (Polymer polymer : polymers) {
                List<Atom> atomList = polymer.getAtomList();
                for (Atom atom : atomList) {
                    moleculeNumber[atom.getXyzIndex() - 1] = current;
                }
                current++;
            }
        }

        // Loop over each molecule
        for (MSNode molecule : molecules.getChildList()) {
            List<Atom> atomList = molecule.getAtomList();
            for (Atom atom : atomList) {
                moleculeNumber[atom.getXyzIndex() - 1] = current;
                atom.setMoleculeNumber(current);
            }
            current++;
        }

        // Loop over each water
        for (MSNode wat : water.getChildList()) {
            List<Atom> atomList = wat.getAtomList();
            for (Atom atom : atomList) {
                moleculeNumber[atom.getXyzIndex() - 1] = current;
            }
            current++;
        }

        // Loop over each ion
        for (MSNode ion : ions.getChildList()) {
            List<Atom> atomList = ion.getAtomList();
            for (Atom atom : atomList) {
                moleculeNumber[atom.getXyzIndex() - 1] = current;
            }
            current++;
        }

        return moleculeNumber;
    }

    /**
     * <p>
     * getAltLocations</p>
     *
     * @return an array of {@link java.lang.String} objects.
     */
    public String[] getAltLocations() {
        if (altLoc == null || altLoc.isEmpty()) {
            return null;
        }

        String[] names = new String[altLoc.size()];
        int i = 0;
        for (String s : altLoc) {
            names[i++] = s;
        }

        return names;
    }

    /**
     * <p>
     * getAtomFromWireVertex</p>
     *
     * @param i a int.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom getAtomFromWireVertex(
            int i) {
        if (atomLookUp != null && atomLookUp.length > i) {
            return atomLookUp[i];
        }
        return null;
    }

    /**
     * <p>
     * getAtomArray</p>
     *
     * @return an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public Atom[] getAtomArray() {
        ArrayList<Atom> atoms = getAtomList();
        Atom[] atomArray = atoms.toArray(new Atom[atoms.size()]);

        // Arrays.sort(atomArray);

        for (int i = 0; i < atoms.size(); i++) {
            atomArray[i].setXyzIndex(i + 1);
        }

        return atomArray;
    }

    /**
     * <p>
     * getActiveAtomArray</p>
     *
     * @return an array of active {@link ffx.potential.bonded.Atom} objects.
     */
    public Atom[] getActiveAtomArray() {
        ArrayList<Atom> atoms = getAtomList();
        ArrayList<Atom> activeAtoms = new ArrayList<>();
        for (Atom a : atoms) {
            if (a.isActive()) {
                activeAtoms.add(a);
            }
        }
        Atom[] atomArray = activeAtoms.toArray(new Atom[activeAtoms.size()]);
        Arrays.sort(atomArray);
        return atomArray;
    }

    /**
     * <p>
     * getBackBoneAtoms</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Atom> getBackBoneAtoms() {
        ArrayList<Atom> backbone = new ArrayList<>();
        List<Residue> residues = getResidueList();
        for (Residue residue : residues) {
            backbone.addAll(residue.getBackboneAtoms());
        }

        /*Atom ca = new Atom("CA");
        ArrayList<ROLS> atoms = this.getList(Atom.class, new ArrayList<ROLS>());
        for (ROLS m : atoms) {
            Atom atom = (Atom) m;
            if (atom.equals(ca)) {
                backbone.add(atom);
                // else if (a.equals(new Atom("C"))) backbone.add(a);
                // else if (a.equals(new Atom("N"))) backbone.add(a);
            }
        }*/
        return backbone;
    }

    /**
     * <p>
     * Getter for the field <code>branchGroup</code>.</p>
     *
     * @return a {@link javax.media.j3d.BranchGroup} object.
     */
    public BranchGroup getBranchGroup() {
        return branchGroup;
    }

    /**
     * <p>
     * getChain</p>
     *
     * @param name a {@link java.lang.String} object.
     * @return a {@link ffx.potential.bonded.Polymer} object.
     */
    public Polymer getChain(String name) {
        for (ListIterator<MSNode> li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            MSNode node = li.next();
            if (node instanceof Polymer) {
                String chainName = node.getName();
                if (chainName.equalsIgnoreCase(name)) {
                    return (Polymer) node;
                }
            }
        }
        return null;
    }

    /**
     * <p>
     * getChainNames</p>
     *
     * @return an array of {@link java.lang.String} objects.
     */
    public String[] getChainNames() {
        ArrayList<String> temp = new ArrayList<>();
        for (ListIterator<MSNode> li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            MSNode node = li.next();
            if (node instanceof Polymer) {
                temp.add(((Polymer) node).getName());
            }

        }
        if (temp.isEmpty()) {
            return null;
        }

        String[] names = new String[temp.size()];
        for (int i = 0; i < temp.size(); i++) {
            names[i] = temp.get(i);
        }

        return names;
    }

    /**
     * <p>
     * getChains</p>
     *
     * @return an array of {@link ffx.potential.bonded.Polymer} objects.
     */
    public Polymer[] getChains() {
        ArrayList<Polymer> polymers = new ArrayList<>();
        for (ListIterator<MSNode> li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            MSNode node = li.next();
            if (node instanceof Polymer) {
                polymers.add((Polymer) node);
            }
        }
        if (polymers.isEmpty()) {
            return null;
        }
        return polymers.toArray(new Polymer[polymers.size()]);
    }

    /**
     * <p>
     * Getter for the field <code>currentCycle</code>.</p>
     *
     * @return a int.
     */
    public int getCurrentCycle() {
        return currentCycle;
    }

    /**
     * <p>
     * Getter for the field <code>cycles</code>.</p>
     *
     * @return a int.
     */
    public int getCycles() {
        return cycles;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getExtent() {
        double[] Rc = {0, 0, 0};
        int num = getAtomList().size();
        for (ListIterator<Atom> li = getAtomList().listIterator(); li.hasNext(); ) {
            li.next().getXYZ(a);
            Rc[0] += a[0];
            Rc[1] += a[1];
            Rc[2] += a[2];
        }

        for (int i = 0; i
                < 3; i++) {
            Rc[i] /= num;
        }

        double r, d = 0;
        double[] xyz = new double[3];
        for (ListIterator<Atom> li = getAtomList().listIterator(); li.hasNext(); ) {
            li.next().getXYZ(xyz);
            VectorMath.diff(xyz, Rc, xyz);
            r = VectorMath.r(xyz);
            if (d < r) {
                d = r;
            }

        }
        return d;
    }

    /**
     * <p>
     * Getter for the field <code>file</code>.</p>
     *
     * @return a {@link java.io.File} object.
     */
    public File getFile() {
        return file;
    }

    /**
     * <p>
     * Getter for the field <code>offset</code>.</p>
     *
     * @return a {@link javax.vecmath.Vector3d} object.
     */
    public Vector3d getOffset() {
        if (offset == null) {
            offset = new Vector3d(0.0, 0.0, 0.0);
        }

        return offset;
    }

    /**
     * <p>
     * Getter for the field <code>originToRot</code>.</p>
     *
     * @return a {@link javax.media.j3d.TransformGroup} object.
     */
    public TransformGroup getOriginToRot() {
        return originToRot;
    }

    public CompositeConfiguration getProperties() {
        return properties == null ? forceField.getProperties() : properties;
    }

    /**
     * <p>
     * getPolymer</p>
     *
     * @param chainID a {@link java.lang.Character} object.
     * @param segID   a {@link java.lang.String} object.
     * @param create  a boolean.
     * @return a {@link ffx.potential.bonded.Polymer} object.
     */
    public Polymer getPolymer(Character chainID, String segID, boolean create) {
        for (ListIterator<MSNode> li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            MSNode node = li.next();
            if (node instanceof Polymer) {
                Polymer polymer = (Polymer) node;
                if (polymer.getName().equals(segID)
                        && polymer.getChainID().equals(chainID)) {
                    return (Polymer) node;
                }
            }
        }
        if (create) {
            Polymer polymer = new Polymer(chainID, segID, true);
            addMSNode(polymer);
            return polymer;
        }

        return null;
    }

    private Atom getResidue(Atom atom, boolean create) {
        return getResidue(atom, create, Residue.ResidueType.UNK);
    }

    private Atom getResidue(Atom atom, boolean create, Residue.ResidueType defaultRT) {
        Character chainID = atom.getChainID();
        String resName = atom.getResidueName();
        int resNum = atom.getResidueNumber();
        String segID = atom.getSegID();
        // Find/Create the chain
        Polymer polymer = getPolymer(chainID, segID, create);
        if (polymer == null) {
            return null;
        }
        Residue res = polymer.getResidue(resName, resNum, create, defaultRT);
        if (create && res != null) {
            return (Atom) res.addMSNode(atom);
        }
        return null;
    }

    /**
     * <p>
     * Getter for the field <code>ions</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<MSNode> getIons() {
        return ions.getChildList();
    }

    /**
     * <p>
     * Getter for the field <code>molecules</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Molecule> getMolecules() {
        ArrayList<Molecule> ret = new ArrayList<>();
        for (MSNode node : molecules.getChildList()) {
            ret.add((Molecule) node);
        }
        return ret;
    }

    /**
     * <p>
     * getWaters</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<MSNode> getWaters() {
        return water.getChildList();
    }

    /**
     * <p>
     * deleteMolecule</p>
     *
     * @param molecule a {@link ffx.potential.bonded.Molecule} object.
     */
    public void deleteMolecule(Molecule molecule) {
        ArrayList<MSNode> list = ions.getChildList();
        for (MSNode node : list) {
            Molecule m = (Molecule) node;
            if (molecule == m) {
                ions.remove(m);
                return;
            }
        }
        list = water.getChildList();
        for (MSNode node : list) {
            Molecule m = (Molecule) node;
            if (molecule == m) {
                water.remove(m);
                return;
            }
        }
        list = molecules.getChildList();
        for (MSNode node : list) {
            Molecule m = (Molecule) node;
            if (molecule == m) {
                molecules.remove(m);
                return;
            }
        }
    }

    public double getMass() {
        Atom atoms[] = this.getAtomArray();
        int nAtoms = atoms.length;
        double mass = 0;
        for (int i = 0; i < nAtoms; i++) {
            mass += atoms[i].getMass();
        }
        return mass;
    }

    private Atom getMolecule(Atom atom, boolean create) {
        String resName = atom.getResidueName();
        int resNum = atom.getResidueNumber();
        Character chainID = atom.getChainID();
        String segID = atom.getSegID();
        ArrayList<MSNode> list = ions.getChildList();
        for (MSNode node : list) {
            Molecule m = (Molecule) node;
            if (m.getSegID().equalsIgnoreCase(segID) && m.getResidueName().equalsIgnoreCase(resName)
                    && m.getResidueNumber() == resNum) {
                return (Atom) m.addMSNode(atom);
            }
        }
        list = water.getChildList();
        for (MSNode node : list) {
            Molecule m = (Molecule) node;
            if (m.getSegID().equalsIgnoreCase(segID)
                    && m.getResidueName().equalsIgnoreCase(resName)
                    && m.getResidueNumber() == resNum) {
                return (Atom) m.addMSNode(atom);
            }
        }
        list = molecules.getChildList();
        for (MSNode node : list) {
            Molecule m = (Molecule) node;
            if (m.getSegID().equalsIgnoreCase(segID)
                    && m.getResidueName().equalsIgnoreCase(resName)
                    && m.getResidueNumber() == resNum) {
                return (Atom) m.addMSNode(atom);
            }
        }
        if (create) {
            Molecule m = new Molecule(resName, resNum, chainID, segID);
            m.addMSNode(atom);
            if (resName.equalsIgnoreCase("DOD")
                    || resName.equalsIgnoreCase("HOH")
                    || resName.equalsIgnoreCase("WAT")) {
                water.add(m);
                // NA, K, MG, MG2, CA, CA2, CL
            } else if (resName.equalsIgnoreCase("NA") || resName.equalsIgnoreCase("K")
                    || resName.equalsIgnoreCase("MG") || resName.equalsIgnoreCase("MG2")
                    || resName.equalsIgnoreCase("CA") || resName.equalsIgnoreCase("CA2")
                    || resName.equalsIgnoreCase("CL") || resName.equalsIgnoreCase("BR")
                    || resName.equalsIgnoreCase("ZN") || resName.equalsIgnoreCase("ZN2")) {
                ions.add(m);
            } else {
                molecules.add(m);
            }
            return atom;
        } else {
            return null;
        }
    }

    /**
     * <p>
     * getResidueList</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<Residue> getResidueList() {
        ArrayList<Residue> residues = new ArrayList<>();
        ListIterator<MSNode> li,
                lj;
        MSNode o;
        Polymer c;
        for (li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            o = li.next();
            if (o instanceof Polymer) {
                c = (Polymer) o;
                for (lj = c.getAtomNodeList().listIterator(); lj.hasNext(); ) {
                    o = lj.next();
                    if (o instanceof Residue) {
                        residues.add((Residue) o);
                    }

                }
            }
        }
        return residues;
    }

    /**
     * <p>
     * getNodeList</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<MSNode> getNodeList() {
        ArrayList<MSNode> residues = new ArrayList<>();
        ListIterator<MSNode> li,
                lj;
        MSNode o;
        Polymer c;
        for (li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            o = li.next();
            if (o instanceof Polymer) {
                c = (Polymer) o;
                for (lj = c.getAtomNodeList().listIterator(); lj.hasNext(); ) {
                    o = lj.next();
                    if (o instanceof Residue) {
                        residues.add(o);
                    }

                }
            }
        }

        ArrayList<MSNode> list = ions.getChildList();
        for (MSNode node : list) {
            residues.add(node);
        }

        list = water.getChildList();
        for (MSNode node : list) {
            residues.add(node);
        }

        list = molecules.getChildList();
        for (MSNode node : list) {
            residues.add(node);
        }
        return residues;
    }

    /**
     * Sums up charge of the system, checking nonstandard residues for
     * non-unitary charges.
     *
     * @param alwaysLog Log non-unitary charge warnings for all nodes
     * @return System charge
     */
    public double getCharge(boolean alwaysLog) {
        double totalCharge = 0;
        for (MSNode node : getNodeList()) {
            double charge = 0;
            boolean isNonstandard = false;
            for (Atom atom : node.getAtomList()) {
                charge += atom.getMultipoleType().getCharge();
                if (atom.isModRes()) {
                    isNonstandard = true;
                }
            }
            if ((alwaysLog || isNonstandard) && (Math.abs(Math.round(charge) - charge) > 1.0E-5)) {
                logger.warning(String.format(" Node %s has non-unitary charge %12.8f", node.toString(), charge));
            }
            totalCharge += charge;
        }
        return totalCharge;
    }

    /**
     * <p>
     * getTransformGroup</p>
     *
     * @return a {@link javax.media.j3d.TransformGroup} object.
     */
    public TransformGroup getTransformGroup() {
        return originToRot;
    }

    /**
     * <p>
     * getWireFrame</p>
     *
     * @return a {@link javax.media.j3d.Node} object.
     */
    public Node getWireFrame() {
        return wire;
    }

    /**
     * <p>
     * isVisible</p>
     *
     * @return a boolean.
     */
    public boolean isVisible() {
        return visible;
    }

    /**
     * <p>
     * loadVRML</p>
     *
     * @return a {@link javax.media.j3d.BranchGroup} object.
     */
    public BranchGroup loadVRML() {
        try {
            VrmlLoader loader = new VrmlLoader();
            VrmlScene scene = null;
            if (vrmlFile != null && vrmlFile.exists()) {
                scene = (VrmlScene) loader.load(vrmlFile.getAbsolutePath());
            } else if (vrmlURL != null) {
                scene = (VrmlScene) loader.load(vrmlURL);
            } else {
                return null;
            }
            BranchGroup bg = scene.getSceneGroup();
            recurseVRML(bg);
            bg.setCapability(BranchGroup.ALLOW_DETACH);
            bg.setCapability(BranchGroup.ALLOW_BOUNDS_READ);
            bg.compile();
            return bg;
        } catch (Exception e) {
            String message = "Fatal exception loading VRML.\n";
            logger.log(Level.SEVERE, message, e);
            System.exit(-1);
            return null;
        }

    }

    /**
     * <p>
     * moveCenter</p>
     *
     * @param d an array of double.
     */
    public void moveCenter(double[] d) {
        for (ListIterator<Atom> li = getAtomList().listIterator(); li.hasNext(); ) {
            (li.next()).move(d);
        }

    }

    private void recurseVRML(Node node) {
        if (node instanceof Shape3D) {
            Shape3D s3d = (Shape3D) node;
            PickTool.setCapabilities(s3d, PickTool.INTERSECT_COORD);
            return;
        } else if (node instanceof SharedGroup) {
            SharedGroup sg = (SharedGroup) node;
            for (Enumeration<Node> e = sg.getAllChildren(); e.hasMoreElements(); ) {
                recurseVRML(e.nextElement());
            }
            return;
        } else if (node instanceof BranchGroup) {
            BranchGroup bg = (BranchGroup) node;
            for (Enumeration<Node> e = bg.getAllChildren(); e.hasMoreElements(); ) {
                recurseVRML(e.nextElement());
            }
            return;
        } else if (node instanceof TransformGroup) {
            TransformGroup vrmlTG1 = (TransformGroup) node;
            for (Enumeration<Node> e = vrmlTG1.getAllChildren(); e.hasMoreElements(); ) {
                node = e.nextElement();
                recurseVRML(node);
            }
            return;
        } else if (node instanceof Link) {
            Link link = (Link) node;
            recurseVRML(link.getSharedGroup());
            return;
        } else if (node instanceof Group) {
            Group group = (Group) node;
            for (Enumeration<Node> e = group.getAllChildren(); e.hasMoreElements(); ) {
                Node n = e.nextElement();
                recurseVRML(n);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void removeLeaves() {
        super.removeLeaves();
        MSNode macroNode = getAtomNode();
        if (macroNode != null) {
            if (macroNode.getChildCount() > 0) {
                getAtomNode().setName(
                        "Macromolecules " + "(" + macroNode.getChildCount() + ")");
            } else if (macroNode.getParent() == this) {
                removeChild(macroNode);
            }
        }

        if (molecules.getChildCount() == 0) {
            removeChild(molecules);
        } else {
            molecules.setName("Hetero Molecules " + "(" + molecules.getChildCount() + ")");
        }

        if (ions.getChildCount() == 0) {
            removeChild(ions);
        } else {
            ions.setName("Ions " + "(" + ions.getChildCount() + ")");
        }

        if (water.getChildCount() == 0) {
            removeChild(water);
        } else {
            water.setName("Water " + "(" + water.getChildCount() + ")");
        }

    }

    private Shape3D renderWire() {
        ArrayList<ROLS> bonds = getBondList();
        int numbonds = bonds.size();
        if (numbonds < 1) {
            return null;
        }

        Vector3d bondmidpoint = new Vector3d();
        double[] mid = {0, 0, 0};
        Vector3d v1 = new Vector3d();
        Vector3d v2 = new Vector3d();
        float[] a1 = {0, 0, 0};
        float[] a2 = {0, 0, 0};
        float[] col = new float[4];
        Bond bond;

        Atom atom1,
                atom2;
        LineArray la = new LineArray(4 * numbonds, GeometryArray.COORDINATES | GeometryArray.COLOR_4 | GeometryArray.NORMALS);
        la.setCapability(LineArray.ALLOW_COORDINATE_WRITE);
        la.setCapability(LineArray.ALLOW_COORDINATE_READ);
        la.setCapability(LineArray.ALLOW_COLOR_WRITE);
        la.setCapability(LineArray.ALLOW_COUNT_READ);
        la.setCapability(LineArray.ALLOW_INTERSECT);
        la.setCapability(LineArray.ALLOW_FORMAT_READ);
        atomLookUp
                = new Atom[4 * numbonds];
        int i = 0;
        col[3] = 0.9f;
        for (ListIterator<ROLS> li = bonds.listIterator(); li.hasNext(); ) {
            bond = (Bond) li.next();
            bond.setWire(la, i);
            atom1 = bond.getAtom(0);
            atom2 = bond.getAtom(1);
            atom1.getV3D(v1);
            atom2.getV3D(v2);
            a1[0] = (float) v1.x;
            a1[1] = (float) v1.y;
            a1[2] = (float) v1.z;
            a2[0] = (float) v2.x;
            a2[1] = (float) v2.y;
            a2[2] = (float) v2.z;
            // Find the bond center
            bondmidpoint.add(v1, v2);
            bondmidpoint.scale(0.5d);
            bondmidpoint.get(mid);

            // Atom #1
            Atom.AtomColor.get(atom1.getAtomicNumber()).get(col);
            atomLookUp[i] = atom1;
            la.setCoordinate(i, a1);
            la.setColor(i, col);
            la.setNormal(i, a2);
            i++;

            atomLookUp[i] = atom1;
            la.setCoordinate(i, mid);
            la.setColor(i, col);
            la.setNormal(i, a2);
            i++;

            // Atom #2
            Atom.AtomColor.get(atom2.getAtomicNumber()).get(col);
            atomLookUp[i] = atom2;
            la.setCoordinate(i, a2);
            la.setColor(i, col);
            la.setNormal(i, a1);
            i++;

            atomLookUp[i] = atom2;
            la.setCoordinate(i, mid);
            la.setColor(i, col);
            la.setNormal(i, a1);
            i++;

        }

        ColoringAttributes cola = new ColoringAttributes(new Color3f(),
                ColoringAttributes.SHADE_GOURAUD);
        Appearance app = new Appearance();
        lineAttributes
                = new LineAttributes();
        lineAttributes.setLineWidth(RendererCache.bondwidth);
        lineAttributes.setCapability(LineAttributes.ALLOW_WIDTH_WRITE);
        lineAttributes.setLineAntialiasingEnable(true);
        app.setLineAttributes(lineAttributes);
        app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_READ);
        app.setCapability(Appearance.ALLOW_LINE_ATTRIBUTES_WRITE);
        RenderingAttributes ra = new RenderingAttributes();
        ra.setAlphaTestValue(0.1f);
        ra.setAlphaTestFunction(RenderingAttributes.GREATER);
        ra.setDepthBufferEnable(true);
        ra.setDepthBufferWriteEnable(true);
        app.setRenderingAttributes(ra);
        app.setColoringAttributes(cola);
        Shape3D wireframe = new Shape3D(la, app);
        // PickTool.setCapabilities(wire, PickTool.INTERSECT_COORD);
        wireframe.setUserData(this);
        wireframe.setBounds(new BoundingSphere(new Point3d(0, 0, 0), 1000.0));
        try {
            wireframe.setBoundsAutoCompute(false);
        } catch (Exception e) {
            e.printStackTrace();
        }

        wireframe.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        wireframe.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
        wireframe.setCapability(Shape3D.ALLOW_LOCAL_TO_VWORLD_READ);
        return wireframe;
    }

    /**
     * Rotate about a point in given in the System's Local Coordinates
     *
     * @param v Vector3d
     */
    public void rotateAbout(Vector3d v) {
        Vector3d newRotPoint = new Vector3d(v);
        originToRot.getTransform(originToRotT3D);
        originToRotT3D.get(originToRotV3D);
        originToRotT3D.setTranslation(new Vector3d(0, 0, 0));
        rotToCOM.getTransform(rotToCOMT3D);
        rotToCOMT3D.get(rotToCOMV3D);
        newRotPoint.add(rotToCOMV3D);
        originToRotT3D.transform(newRotPoint);
        newRotPoint.add(originToRotV3D);
        originToRotT3D.setTranslation(newRotPoint);
        rotToCOMV3D.set(v);
        rotToCOMV3D.negate();
        rotToCOMT3D.setTranslation(rotToCOMV3D);
        originToRot.setTransform(originToRotT3D);
        rotToCOM.setTransform(rotToCOMT3D);
    }

    /**
     * <p>
     * sceneGraphChange</p>
     *
     * @param newShapes a {@link java.util.List} object.
     */
    public void sceneGraphChange(List<BranchGroup> newShapes) {
        if (newShapes == null) {
            newShapes = myNewShapes;
        }

        if (newShapes.isEmpty()) {
            return;
        }

        boolean reCompile = false;
        // Check for nodes (new and/or recycled) being added to this
        // MolecularAssembly
        for (ListIterator<BranchGroup> li = newShapes.listIterator(); li.hasNext(); ) {
            BranchGroup group = li.next();
            li.remove();
            // This is code for cycling between two MolecularAssemblies
            if (group.getUserData() != null) {
                logger.info(format("%s %s", group.toString(), group.getUserData().toString()));
                /*
                 * Object userData = group.getUserData(); if (userData!=this) {
                 * // The appearance has already been set during a recursive
                 * call to setView, // although we need to turn back on
                 * Pickablility. TransformGroup tg = (TransformGroup)
                 * group.getChild(0); Shape3D shape = (Shape3D) tg.getChild(0);
                 * shape.setPickable(true); group.setUserData(this); if
                 * (!reCompile) { if (childNodes.isLive()) {
                 * childNodes.detach(); } reCompile = true; }
                 * childNodes.moveTo(group);
                 */
            } else {
                // This is a new group since it has no userData.
                // We can not query for the identity of its parent later, so
                // we will store it as a userData reference.
                group.setUserData(this);
                if (!reCompile) {
                    if (childNodes.isLive()) {
                        childNodes.detach();
                    }

                    reCompile = true;
                }

                childNodes.addChild(group);
            }

        }
        if (reCompile) {
            childNodes.compile();
            base.addChild(childNodes);
        }

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color,
                         Material mat) {
        for (ListIterator<MSNode> li = getAtomNodeList().listIterator(); li.hasNext(); ) {
            MSGroup group = (MSGroup) li.next();
            group.setColor(newColorModel, color, mat);
        }

        for (MSNode m : molecules.getChildList()) {
            m.setColor(newColorModel, color, mat);
        }

        for (MSNode m : water.getChildList()) {
            m.setColor(newColorModel, color, mat);
        }

        for (MSNode m : ions.getChildList()) {
            m.setColor(newColorModel, color, mat);
        }

    }

    /**
     * <p>
     * Setter for the field <code>currentCycle</code>.</p>
     *
     * @param c a int.
     */
    public void setCurrentCycle(int c) {
        if (c <= cycles && c > 0) {
            currentCycle = c;
            for (ListIterator<Atom> li = getAtomList().listIterator(); li.hasNext(); ) {
                (li.next()).setCurrentCycle(currentCycle);
            }

        }
    }

    /**
     * <p>
     * Setter for the field <code>cycles</code>.</p>
     *
     * @param c a int.
     */
    public void setCycles(int c) {
        cycles = c;
    }

    /**
     * <p>
     * Setter for the field <code>file</code>.</p>
     *
     * @param f a {@link java.io.File} object.
     */
    public void setFile(File f) {
        if (f == null) {
            return;
        }
        file = f;
    }

    /**
     * <p>
     * Setter for the field <code>offset</code>.</p>
     *
     * @param o a {@link javax.vecmath.Vector3d} object.
     */
    public void setOffset(Vector3d o) {
        offset = o;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setView(RendererCache.ViewModel newViewModel,
                        List<BranchGroup> newShapes) {
        // Just Detach the whole system branch group
        if (newViewModel == RendererCache.ViewModel.DESTROY) {
            if (switchGroup != null) {
                switchGroup.setWhichChild(Switch.CHILD_NONE);
            }
            visible = false;
        } else if (newViewModel == RendererCache.ViewModel.SHOWVRML) {
            switchGroup.setWhichChild(Switch.CHILD_ALL);
        } else if (newViewModel == RendererCache.ViewModel.HIDEVRML) {
            switchGroup.setWhichChild(0);
        } else {
            setWireWidth(RendererCache.bondwidth);
            if (newViewModel == RendererCache.ViewModel.DETAIL && childNodes.isLive()) {
                childNodes.detach();
            }
            /**
             * We'll collect new Scenegraph Shapes in our newShapeNode This is
             * to avoid the case where setView is called from the root node and
             * all new shapes for every MolecularAssembly would then be put into
             * the same ArrayList.
             */
            super.setView(newViewModel, myNewShapes);
            ArrayList<ROLS> moleculeList = getList(Molecule.class, new ArrayList<>());
            for (ROLS m : moleculeList) {
                m.setView(newViewModel, myNewShapes);
            }
            for (MSNode m : molecules.getChildList()) {
                m.setView(newViewModel, myNewShapes);
            }
            for (MSNode m : water.getChildList()) {
                m.setView(newViewModel, myNewShapes);
            }
            for (MSNode m : ions.getChildList()) {
                m.setView(newViewModel, myNewShapes);
            }
            if (newViewModel == RendererCache.ViewModel.INVISIBLE) {
                switchGroup.setWhichChild(0);
            }
            if (newViewModel == RendererCache.ViewModel.DETAIL) {
                childNodes.compile();
                base.addChild(childNodes);
            }
        }
    }

    /**
     * <p>
     * setVRML</p>
     *
     * @param v a {@link javax.media.j3d.BranchGroup} object.
     */
    public void setVRML(BranchGroup v) {
        vrmlURL = null;
        vrmlFile = null;
        vrml = v;
    }

    /**
     * <p>
     * setVRML</p>
     *
     * @param file a {@link java.io.File} object.
     */
    public void setVRML(File file) {
        vrmlFile = file;
        vrmlURL
                = null;
    }

    /**
     * <p>
     * setVRML</p>
     *
     * @param url a {@link java.net.URL} object.
     */
    public void setVRML(URL url) {
        vrmlURL = url;
        vrmlFile
                = null;
    }

    /**
     * <p>
     * setWireWidth</p>
     *
     * @param f a float.
     */
    public void setWireWidth(float f) {
        if (wire == null) {
            return;
        }

        lineAttributes.setLineWidth(f);
    }

    /**
     * <p>
     * sidePolymerCOM</p>
     */
    public void sidePolymerCOM() {
        ArrayList<Residue> residues = getResidueList();
        Residue r;

        ListIterator<Residue> li;

        for (li = residues.listIterator(); li.hasNext(); ) {
            r = li.next();
            r.logSideChainCOM();
        }
    }

    /**
     * Moves the center of all chemical entities into the primary
     * unit cell. Somewhat experimental feature; use with caution.
     */
    public void moveAllIntoUnitCell() {
        moveIntoUnitCell(getChains());
        moveIntoUnitCell(getWaters());
        moveIntoUnitCell(getIons());

        // Unsure why this mapping is flagged as necessary
        List<MSNode> molNodes = getMolecules().stream().
                map((Molecule m) -> (MSNode) m).
                collect(Collectors.toList());
        moveIntoUnitCell(molNodes);
    }

    private void moveIntoUnitCell(MSNode[] groups) {
        if (groups != null && groups.length > 0) {
            moveIntoUnitCell(Arrays.asList(groups));
        }
    }

    /**
     * Move the center of each listed chemical entity into the primary unit cell.
     *
     * @param groups
     */
    private void moveIntoUnitCell(List<MSNode> groups) {
        Crystal cryst = getCrystal();
        if (cryst.aperiodic()) {
            return;
        }

        for (MSNode group : groups) {
            double[] com = new double[3];
            double[] xyz = new double[3];
            double[] translate = new double[3];

            List<Atom> atoms = group.getAtomList();
            double totMass = 0;

            for (Atom atom : atoms) {
                double mass = atom.getMass();
                totMass += mass;
                xyz = atom.getXYZ(xyz);
                for (int i = 0; i < 3; i++) {
                    double diff = xyz[i] - com[i];
                    diff /= totMass;
                    com[i] += diff;
                }
            }

            cryst.toPrimaryCell(com, translate);
            for (int i = 0; i < 3; i++) {
                translate[i] -= com[i];
            }

            for (Atom atom : atoms) {
                atom.move(translate);
            }
        }
    }
}
