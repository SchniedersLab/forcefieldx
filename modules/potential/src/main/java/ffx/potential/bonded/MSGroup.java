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
package ffx.potential.bonded;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import java.util.logging.Logger;

import static java.lang.String.format;

import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Material;
import org.jogamp.vecmath.Color3f;

import ffx.numerics.math.VectorMath;
import ffx.potential.parameters.ForceField;

/**
 * The MSGroup class has one subnode containing atoms, and one that contains
 * molecular mechanics/geometry terms.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings({"serial", "CloneableImplementsClone"})
public abstract class MSGroup extends MSNode {

    private static final Logger logger = Logger.getLogger(MSGroup.class.getName());
    /**
     * Constructs the Geometry of this MultiScaleGroup and stores the terms as children in the Term node.
     */
    // Atoms Node
    private MSNode atomNode = new MSNode("Atoms");
    // Valence Energy Terms
    private MSNode termNode = new MSNode("Valence Terms");
    private MSNode bondNode = new MSNode("Bonds");
    private MSNode angleNode = new MSNode("Angles");
    private MSNode stretchBendNode = new MSNode("Stretch-Bends");
    private MSNode ureyBradleyNode = new MSNode("Urey-Bradleys");
    private MSNode outOfPlaneBendNode = new MSNode("Out-of-Plane Bends");
    private MSNode torsionNode = new MSNode("Torsions");
    private MSNode stretchTorsionNode = new MSNode("Stretch-Torsions");
    private MSNode angleTorsionNode = new MSNode("Angle-Torsions");
    private MSNode piOrbitalTorsionNode = new MSNode("Pi-Orbital Torsions");
    private MSNode torsionTorsionNode = new MSNode("Torsion-Torsions");
    private MSNode improperTorsionNode = new MSNode("Improper Torsions");
    /**
     * List of Joints.
     */
    private ArrayList<Joint> joints = new ArrayList<>();
    /**
     * Whether the terms are current.
     */
    private boolean finalized;
    /**
     * Center of the MultiScaleGroup
     */
    private double[] center;
    /**
     * List of underconstrained Atoms
     */
    private ArrayList<Atom> dangelingatoms;
    /**
     * Constant <code>bondTime=0</code>
     */
    protected static long bondTime = 0;
    /**
     * Constant <code>angleTime=0</code>
     */
    protected static long angleTime = 0;
    /**
     * Constant <code>stretchBendTime=0</code>
     */
    protected static long stretchBendTime = 0;
    /**
     * Constant <code>ureyBradleyTime=0</code>
     */
    protected static long ureyBradleyTime = 0;
    /**
     * Constant <code>outOfPlaneBendTime=0</code>
     */
    protected static long outOfPlaneBendTime = 0;
    /**
     * Constant <code>torsionTime=0</code>
     */
    protected static long torsionTime = 0;
    /**
     * Constant <code>stretchTorsionTime=0</code>
     */
    protected static long stretchTorsionTime = 0;
    /**
     * Constant <code>angleTorsionTime=0</code>
     */
    protected static long angleTorsionTime = 0;
    /**
     * Constant <code>piOrbitalTorsionTime=0</code>
     */
    protected static long piOrbitalTorsionTime = 0;
    /**
     * Constant <code>torsionTorsionTime=0</code>
     */
    protected static long torsionTorsionTime = 0;
    /**
     * Constant <code>torsionTorsionTime=0</code>
     */
    protected static long improperTorsionTime = 0;

    /**
     * Default Constructor initializes a MultiScaleGroup and a few of its
     * sub-nodes.
     */
    public MSGroup() {
        super("", 2);
        finalized = false;
        termNode.add(bondNode);
        termNode.add(angleNode);
        termNode.add(stretchBendNode);
        termNode.add(ureyBradleyNode);
        termNode.add(outOfPlaneBendNode);
        termNode.add(torsionNode);
        termNode.add(stretchTorsionNode);
        termNode.add(angleTorsionNode);
        termNode.add(piOrbitalTorsionNode);
        termNode.add(torsionTorsionNode);
        termNode.add(improperTorsionNode);
        add(atomNode);
        add(termNode);
    }

    /**
     * Constructs a MultiScaleGroup object with name n.
     *
     * @param n a {@link java.lang.String} object.
     */
    public MSGroup(String n) {
        this();
        setName(n);
    }

    /**
     * Constructs a MultiScaleGroup object with name n and sets its AtomGroup
     * node equals to node.
     *
     * @param n    a {@link java.lang.String} object.
     * @param node a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSGroup(String n, MSNode node) {
        this(n);
        atomNode = node;
        add(atomNode);
    }

    /**
     * Abstract method that should specify how to add various MSNodes subclasses
     * (such as Atoms, Residues and Polymers) to the MSGroup
     *
     * @param m a {@link ffx.potential.bonded.MSNode} object.
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public abstract MSNode addMSNode(MSNode m);

    /**
     * <p>
     * reOrderAtoms</p>
     */
    public void reOrderAtoms() {
        ArrayList<Atom> atomList = getAtomList();
        int nAtoms = atomList.size();
        Atom[] atoms = new Atom[nAtoms];
        atoms = atomList.toArray(atoms);

        boolean sorted;
        do {
            sorted = true;
            for (int i = 1; i < nAtoms; i++) {
                Atom a1 = atoms[i - 1];
                Atom a2 = atoms[i];
                if (a1.getName().compareToIgnoreCase(a2.getName()) > 0) {
                    int i1 = a1.getIndex();
                    int i2 = a2.getIndex();
                    atoms[i - 1] = a2;
                    atoms[i] = a1;
                    a1.setXyzIndex(i2);
                    a2.setXyzIndex(i1);
                    sorted = false;
                }
            }
        } while (!sorted);
    }

    /**
     * Special-case check for if a Bond connects atoms in a hetero-group polymer linkage.
     *
     * @param polyLinks polymerlink records from properties
     * @param bond      Bond to check
     * @return If bond should be special-cased
     */
    private boolean polymerLinked(String[] polyLinks, Bond bond) {
        Atom a1 = bond.getAtom(0);
        Atom a2 = bond.get1_2(a1);

        String resName = a1.getResidueName();
        if (resName == null) {
            logger.warning(format(" Residue name of atom 1 (%s) is null!", a1));
            return false;
        } else {
            String r2name = a2.getResidueName();
            if (r2name == null) {
                logger.warning(format(" Residue name of atom 2 (%s) is null!", a2));
                return false;
            } else if (!resName.equals(r2name)) {
                return false;
            }
        }
        String name1 = a1.getName();
        String name2 = a2.getName();

        for (String pLink : polyLinks) {
            String[] toks = pLink.split("\\s+");
            if (resName.equalsIgnoreCase(toks[0])) {
                if (name1.equalsIgnoreCase(toks[1]) && name2.equalsIgnoreCase(toks[2])) {
                    return true;
                } else if (name1.equalsIgnoreCase(toks[2]) && name2.equalsIgnoreCase(toks[1])) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * <p>
     * assignBondedTerms</p>
     *
     * @param forceField the ForceField to use when creating bonded terms.
     */
    public void assignBondedTerms(ForceField forceField) {
        MSNode newBondNode = new MSNode("Bonds");
        MSNode newAngleNode = new MSNode("Angles");
        MSNode newStretchBendNode = new MSNode("Stretch-Bends");
        MSNode newUreyBradleyNode = new MSNode("Urey-Bradleys");
        MSNode newOutOfPlaneBendNode = new MSNode("Out-of-Plane Bends");
        MSNode newTorsionNode = new MSNode("Torsions");
        MSNode newStretchTorsionNode = new MSNode("Stretch-Torsions");
        MSNode newAngleTorsionNode = new MSNode("Angle-Torsions");
        MSNode newPiOrbitalTorsionNode = new MSNode("Pi-Orbital Torsions");
        MSNode newTorsionTorsionNode = new MSNode("Torsion-Torsions");
        MSNode newImproperTorsionNode = new MSNode("Improper Torsions");

        String[] polyLinks = forceField.getProperties().getStringArray("polymerlink");

        // Collect all bonds for which both atoms are in this Group
        long time = System.nanoTime();
        ArrayList<Bond> bonds = new ArrayList<>();
        for (Atom atom : getAtomList()) {
            // This check is probably unnecessary.
            if (atom.getNumBonds() != 0) {
                for (Bond bond : atom.getBonds()) {
                    if ((bond.sameGroup() || polymerLinked(polyLinks, bond)) && bond.getParent() == null) {
                        newBondNode.add(bond);
                        bonds.add(bond);
                    }
                }
            }
        }
        newBondNode.setName("Bonds (" + newBondNode.getChildCount() + ")");
        setBonds(newBondNode);
        bondTime += System.nanoTime() - time;

        // Find intra-group angles.
        time = System.nanoTime();
        ArrayList<Angle> angles = new ArrayList<>();
        for (Atom atom : getAtomList()) {
            if (atom.getNumBonds() != 0) {
                int index = 0;
                for (Bond bond : atom.getBonds()) {
                    index++;
                    if (bond.sameGroup() || polymerLinked(polyLinks, bond)) {
                        for (ListIterator<Bond> li = atom.getBonds().listIterator(index); li.hasNext(); ) {
                            Bond bond2 = li.next();
                            if (bond2.sameGroup() || polymerLinked(polyLinks, bond2)) {
                                Angle newAngle = Angle.angleFactory(bond, bond2, forceField);
                                newAngleNode.insert(newAngle, 0);
                                angles.add(newAngle);
                            }
                        }
                    }
                }
            }
        }
        newAngleNode.setName("Angles (" + newAngleNode.getChildCount() + ")");
        setAngles(newAngleNode);
        angleTime += System.nanoTime() - time;

        // Find Stretch-Bends.
        time = System.nanoTime();
        for (Angle angle : angles) {
            StretchBend newStretchBend = StretchBend.stretchBendFactory(angle, forceField);
            if (newStretchBend != null) {
                newStretchBendNode.insert(newStretchBend, 0);
            }
        }
        newStretchBendNode.setName("Stretch-Bends (" + newStretchBendNode.getChildCount() + ")");
        setStretchBends(newStretchBendNode);
        stretchBendTime += System.nanoTime() - time;

        // Urey-Bradleys.
        time = System.nanoTime();
        for (Angle angle : angles) {
            UreyBradley newUreyBradley = UreyBradley.ureyBradlyFactory(angle, forceField);
            if (newUreyBradley != null) {
                newUreyBradleyNode.insert(newUreyBradley, 0);
            }
        }
        newUreyBradleyNode.setName("Urey-Bradleys (" + newUreyBradleyNode.getChildCount() + ")");
        setUreyBradleys(newUreyBradleyNode);
        ureyBradleyTime += System.nanoTime() - time;

        // Out-of-Plane Bends.
        time = System.nanoTime();
        for (Angle angle : angles) {
            OutOfPlaneBend opBend = OutOfPlaneBend.outOfPlaneBendFactory(angle, forceField);
            if (opBend != null) {
                newOutOfPlaneBendNode.insert(opBend, 0);
            }
        }
        newOutOfPlaneBendNode.setName("Out-of-Plane Bends (" + newOutOfPlaneBendNode.getChildCount() + ")");
        setOutOfPlaneBends(newOutOfPlaneBendNode);
        outOfPlaneBendTime += System.nanoTime() - time;

        // Find Intra-Group Torsions.
        time = System.nanoTime();
        ArrayList<Torsion> torsions = new ArrayList<>();
        for (Bond middleBond : bonds) {
            Atom atom1 = middleBond.getAtom(0);
            Atom atom2 = middleBond.getAtom(1);
            if (atom1.getNumBonds() != 0 && atom2.getNumBonds() != 0) {
                for (Bond bond1 : atom1.getBonds()) {
                    if (bond1 != middleBond) {
                        for (Bond bond3 : atom2.getBonds()) {
                            if (bond3 != middleBond) {
                                Torsion torsion = Torsion.torsionFactory(bond1, middleBond, bond3, forceField);
                                if (torsion != null) {
                                    newTorsionNode.add(torsion);
                                    torsions.add(torsion);
                                }
                            }
                        }
                    }
                }
            }
        }
        newTorsionNode.setName("Torsions (" + newTorsionNode.getChildCount() + ")");
        setTorsions(newTorsionNode);
        torsionTime += System.nanoTime() - time;

        // Find Intra-Group Stretch-Torsions.
        time = System.nanoTime();
        for (Torsion torsion : torsions) {
            StretchTorsion stretchTorsion = StretchTorsion.stretchTorsionFactory(torsion, forceField);
            if (stretchTorsion != null) {
                newStretchTorsionNode.add(stretchTorsion);
            }
        }
        newStretchTorsionNode.setName("Stretch-Torsions (" + newStretchTorsionNode.getChildCount() + ")");
        setStretchTorsions(newStretchTorsionNode);
        stretchTorsionTime += System.nanoTime() - time;

        // Find Intra-Group Angle-Torsions.
        time = System.nanoTime();
        for (Torsion torsion : torsions) {
            AngleTorsion angleTorsion = AngleTorsion.angleTorsionFactory(torsion, forceField);
            if (angleTorsion != null) {
                newAngleTorsionNode.add(angleTorsion);
            }
        }
        newAngleTorsionNode.setName("Angle-Torsions (" + newAngleTorsionNode.getChildCount() + ")");
        setAngleTorsions(newAngleTorsionNode);
        angleTorsionTime += System.nanoTime() - time;

        // Find Pi-Orbital Torsions.
        time = System.nanoTime();
        for (Bond bond : bonds) {
            PiOrbitalTorsion piOrbitalTorsion = PiOrbitalTorsion.piOrbitalTorsionFactory(bond, forceField);
            if (piOrbitalTorsion != null) {
                newPiOrbitalTorsionNode.add(piOrbitalTorsion);
            }
        }
        newPiOrbitalTorsionNode.setName("Pi-Orbital Torsions (" + newPiOrbitalTorsionNode.getChildCount() + ")");
        setPiOrbitalTorsions(newPiOrbitalTorsionNode);
        piOrbitalTorsionTime += System.nanoTime() - time;

        // Find Improper-Torsions.
        time = System.nanoTime();
        ArrayList<Atom> atoms = getAtomList();
        for (Atom atom : atoms) {
            if (atom.isTrigonal()) {
                ArrayList<ImproperTorsion> improperTorsions = ImproperTorsion.improperTorsionFactory(atom, forceField);
                if (improperTorsions != null) {
                    for (ImproperTorsion improperTorsion : improperTorsions) {
                        newImproperTorsionNode.add(improperTorsion);
                    }
                }
            }
        }
        newImproperTorsionNode.setName("Improper Torsions (" + newImproperTorsionNode.getChildCount() + ")");
        setImproperTorsions(newImproperTorsionNode);
        improperTorsionTime += System.nanoTime() - time;

        // Find Torsion-Torsions.
        time = System.nanoTime();
        for (Angle angle : angles) {
            Atom atom1 = angle.atoms[0];
            Atom atom2 = angle.atoms[1];
            Atom atom3 = angle.atoms[2];
            for (Bond firstBond : atom1.getBonds()) {
                Atom atom0 = firstBond.get1_2(atom1);
                if (atom0 != atom2 && atom0 != atom3) {
                    for (Bond lastBond : atom3.getBonds()) {
                        Atom atom4 = lastBond.get1_2(atom3);
                        if (atom4 != atom0 && atom4 != atom1 && atom4 != atom2) {
                            TorsionTorsion torsionTorsion = TorsionTorsion.
                                    torsionTorsionFactory(firstBond, angle, lastBond, forceField);
                            if (torsionTorsion != null) {
                                newTorsionTorsionNode.insert(torsionTorsion, 0);
                            }
                        }
                    }
                }
            }
        }
        newTorsionTorsionNode.setName(
                "Torsion-Torsions (" + newTorsionTorsionNode.getChildCount() + ")");
        setTorsionTorsions(newTorsionTorsionNode);
        torsionTorsionTime += System.nanoTime() - time;
        int numberOfValenceTerms = newBondNode.getChildCount() + newAngleNode.getChildCount()
                + newStretchBendNode.getChildCount() + newUreyBradleyNode.getChildCount()
                + newOutOfPlaneBendNode.getChildCount() + newTorsionNode.getChildCount()
                + newStretchTorsionNode.getChildCount() + newAngleTorsionNode.getChildCount()
                + newPiOrbitalTorsionNode.getChildCount() + newTorsionTorsionNode.getChildCount()
                + newImproperTorsionNode.getChildCount();
        termNode.setName(
                "Valence Terms (" + numberOfValenceTerms + ")");
    }

    /**
     * <p>
     * constructValenceTerms</p>
     */
    public void constructValenceTerms() {
        MSNode b = new MSNode("Bonds");
        MSNode a = new MSNode("Angles");
        MSNode d = new MSNode("Dihedrals");
        int index = 0;
        double[] da = new double[3];
        double[] db = new double[3];
        ArrayList<Atom> atomList = getAtomList();
        for (Atom a1 : atomList) {
            index++;
            for (ListIterator li = atomList.listIterator(index); li.hasNext(); ) {
                Atom a2 = (Atom) li.next();
                a1.getXYZ(da);
                a2.getXYZ(db);
                double d1 = VectorMath.dist(da, db);
                double d2 = Bond.BUFF + a1.getVDWR() / 2 + a2.getVDWR() / 2;
                if (d1 < d2) {
                    ArrayList<Bond> adjunctBonds = new ArrayList<>();
                    if (a1.getNumBonds() > 0) {
                        adjunctBonds.addAll(a1.getBonds());
                    }
                    if (a2.getNumBonds() > 0) {
                        adjunctBonds.addAll(a2.getBonds());
                    }
                    Bond newbond = new Bond(a1, a2);
                    b.add(newbond);
                    for (Bond adjunctBond : adjunctBonds) {
                        if (newbond == adjunctBond) {
                            logger.info("New Bond = Adjunct Bond");
                        } else {
                            Angle newangle = new Angle(newbond, adjunctBond);
                            a.add(newangle);
                            Atom atom13 = adjunctBond.getOtherAtom(newbond);
                            for (Bond bond14 : atom13.getBonds()) {
                                if (bond14 != adjunctBond) {
                                    d.add(new Torsion(newangle, bond14));
                                }
                            }
                        }
                    }
                }
            }
        }
        setBonds(b);
        setAngles(a);
        setTorsions(d);
    }

    /**
     * Create a joint between two chemical groups.
     *
     * @param bond       Bond
     * @param group1     a {@link ffx.potential.bonded.MSGroup} object.
     * @param group2     a {@link ffx.potential.bonded.MSGroup} object.
     * @param forceField the ForceField parameters to use when creating the
     *                   joint.
     * @return Joint the created Joint.
     */
    public Joint createJoint(Bond bond, MSGroup group1, MSGroup group2, ForceField forceField) {
        if (bond == null || group1 == null || group2 == null) {
            return null;
        }

        MSNode newBondNode = new MSNode("Bonds");
        MSNode newAngleNode = new MSNode("Angles");
        MSNode newStretchBendNode = new MSNode("Stretch-Bends");
        MSNode newUreyBradleyNode = new MSNode("Urey-Bradleys");
        MSNode newOutOfPlaneNode = new MSNode("Out-of-Plane Bends");
        MSNode newTorsionNode = new MSNode("Torsions");
        MSNode newStretchTorsionNode = new MSNode("Stretch-Torsions");
        MSNode newAngleTorsionNode = new MSNode("Angle-Torsions");
        MSNode newPiOrbitalTorsionNode = new MSNode("Pi-Orbital Torsions");
        MSNode newTorsionTorsionNode = new MSNode("Torsion-Torsions");
        //MSNode newImproperTorsionNode = new MSNode("Improper Torsions");
        newBondNode.add(bond);
        newBondNode.setName("Bonds (" + newBondNode.getChildCount() + ")");

        // Collect Angles that include the joining bond(s)
        ArrayList<Angle> angles = new ArrayList<>();
        // Chemical Group #1
        Atom atom1 = bond.getAtom(0);
        for (Bond bond2 : atom1.getBonds()) {
            if (bond != bond2 && bond.getOtherAtom(bond2) != null) {
                Angle newAngle = Angle.angleFactory(bond, bond2, forceField);
                if (newAngle != null) {
                    newAngleNode.add(newAngle);
                    angles.add(newAngle);
                }
            }
        }
        // Chemical Group #2
        Atom atom2 = bond.getAtom(1);
        for (Bond bond2 : atom2.getBonds()) {
            if (bond != bond2 && bond.getOtherAtom(bond2) != null) {
                Angle newAngle = Angle.angleFactory(bond, bond2, forceField);
                if (newAngle != null) {
                    newAngleNode.add(newAngle);
                    angles.add(newAngle);
                }
            }
        }
        newAngleNode.setName("Angles (" + newAngleNode.getChildCount() + ")");

        for (Angle angle : angles) {
            StretchBend stretchBend = StretchBend.stretchBendFactory(angle, forceField);
            if (stretchBend != null) {
                newStretchBendNode.add(stretchBend);
            }
            UreyBradley ureyBradley = UreyBradley.ureyBradlyFactory(angle, forceField);
            if (ureyBradley != null) {
                newUreyBradleyNode.add(ureyBradley);
            }
            OutOfPlaneBend outOfPlaneBend = OutOfPlaneBend.outOfPlaneBendFactory(angle, forceField);
            if (outOfPlaneBend != null) {
                newOutOfPlaneNode.add(outOfPlaneBend);
            }
        }
        newStretchBendNode.setName(
                "Stretch-Bends (" + newStretchBendNode.getChildCount() + ")");
        newUreyBradleyNode.setName(
                "Urey-Bradleys (" + newUreyBradleyNode.getChildCount() + ")");
        newOutOfPlaneNode.setName(
                "Out-of-Plane Bends (" + newOutOfPlaneNode.getChildCount() + ")");

        // Find torsions across the joint.
        atom1 = bond.getAtom(0);
        atom2 = bond.getAtom(1);
        ArrayList<Torsion> torsions = new ArrayList<>();
        if (atom1.getNumBonds() != 0 && atom2.getNumBonds() != 0) {
            for (Bond firstBond : atom1.getBonds()) {
                if (firstBond != bond) {
                    for (Bond lastBond : atom2.getBonds()) {
                        if (lastBond != bond) {
                            Torsion torsion = Torsion.torsionFactory(firstBond, bond, lastBond, forceField);
                            if (torsion != null) {
                                newTorsionNode.add(torsion);
                                torsions.add(torsion);
                            }
                        }
                    }
                }
            }
        }
        newTorsionNode.setName(
                "Torsional Angles (" + newTorsionNode.getChildCount() + ")");

        for (Torsion torsion : torsions) {
            StretchTorsion stretchTorsion = StretchTorsion.stretchTorsionFactory(torsion, forceField);
            if (stretchTorsion != null) {
                newStretchTorsionNode.add(stretchTorsion);
            }
            AngleTorsion angleTorsion = AngleTorsion.angleTorsionFactory(torsion, forceField);
            if (angleTorsion != null) {
                newAngleTorsionNode.add(angleTorsion);
            }
        }
        newStretchTorsionNode.setName("Stretch-Torsions (" + stretchTorsionNode.getChildCount() + ")");
        newAngleTorsionNode.setName("Angle-Torsions (" + angleTorsionNode.getChildCount() + ")");

        // Find Pi-Orbital Torsions across the joint.
        PiOrbitalTorsion piOrbitalTorsion = PiOrbitalTorsion.piOrbitalTorsionFactory(bond, forceField);
        if (piOrbitalTorsion != null) {
            newPiOrbitalTorsionNode.add(piOrbitalTorsion);
        }
        newPiOrbitalTorsionNode.setName(
                "Pi-Orbital Torsions (" + newPiOrbitalTorsionNode.getChildCount() + ")");

        // Find Torsion-Torsions across the joint.
        for (Angle angle : angles) {
            atom1 = angle.atoms[0];
            atom2 = angle.atoms[1];
            Atom atom3 = angle.atoms[2];
            for (Bond firstBond : atom1.getBonds()) {
                Atom atom0 = firstBond.get1_2(atom1);
                if (atom0 != atom2 && atom0 != atom3) {
                    for (Bond lastBond : atom3.getBonds()) {
                        Atom atom4 = lastBond.get1_2(atom3);
                        if (atom4 != atom0 && atom4 != atom1 && atom4 != atom2) {
                            TorsionTorsion torsionTorsion = TorsionTorsion.
                                    torsionTorsionFactory(firstBond, angle, lastBond, forceField);
                            if (torsionTorsion != null) {
                                newTorsionTorsionNode.add(torsionTorsion);
                            }
                        }
                    }
                }
            }
        }

        newTorsionTorsionNode.setName(
                "Torsion-Torsions (" + newTorsionTorsionNode.getChildCount() + ")");

        Joint newJoint = new Joint(group1, group2, newBondNode, newAngleNode,
                newStretchBendNode, newUreyBradleyNode, newOutOfPlaneNode,
                newTorsionNode, newStretchTorsionNode, newAngleTorsionNode,
                newPiOrbitalTorsionNode, newTorsionTorsionNode);

        group1.addJoint(newJoint);
        group2.addJoint(newJoint);

        return newJoint;
    }

    /**
     * <p>addJoint.</p>
     *
     * @param newJoint a {@link ffx.potential.bonded.Joint} object.
     */
    void addJoint(Joint newJoint) {
        joints.add(newJoint);
    }

    /**
     * <p>clearJoints.</p>
     */
    public void clearJoints() {
        joints.clear();
    }

    /**
     * <p>Getter for the field <code>joints</code>.</p>
     *
     * @return a {@link java.util.ArrayList} object.
     */
    ArrayList<Joint> getJoints() {
        return joints;
    }

    /**
     * Joiner joins Moieties m1 and m2 and returns the Geometry objects formed
     * in a Joint.
     *
     * @param group1     a {@link ffx.potential.bonded.MSGroup} object.
     * @param group2     a {@link ffx.potential.bonded.MSGroup} object.
     * @param forceField the ForceField parameters to use when creating the
     *                   joint.
     * @return a {@link ffx.potential.bonded.Joint} object.
     */
    public Joint createJoint(MSGroup group1, MSGroup group2, ForceField forceField) {
        Joint joint = null;
        double[] da = new double[3];
        double[] db = new double[3];
        for (Atom a1 : group1.getAtomList()) {
            a1.getXYZ(da);
            for (Atom a2 : group2.getAtomList()) {
                a2.getXYZ(db);
                double d1 = VectorMath.dist(da, db);
                double d2 = Bond.BUFF + a1.getVDWR() / 2 + a2.getVDWR() / 2;
                if (d1 < d2) {
                    Bond b = new Bond(a1, a2);
                    Joint newJoint = createJoint(b, group1, group2, forceField);
                    if (joint != null) {
                        joint.merge(newJoint);
                    } else {
                        joint = newJoint;
                    }
                }
            }
        }
        return joint;
    }

    /**
     * Abstract method that should specify how to finalize a MSGroup
     *
     * @param finalizeGroups a boolean.
     * @param forceField     the ForceField parameters to use when finalizing the
     *                       MSGroup.
     */
    public abstract void finalize(boolean finalizeGroups, ForceField forceField);

    /**
     * This method constructs an ArrayList of atoms which are under-constrained.
     * (i.e. They can except more bonds)
     */
    public void findDangelingAtoms() {
        ArrayList<Atom> d = new ArrayList<>();
        for (Atom a : getAtomList()) {
            if (a.isDangeling()) {
                d.add(a);
            }
        }
        setDangelingAtoms(d);
    }

    /**
     * Returns the MultiScaleGroup's angles FNode.
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getAngles() {
        return angleNode;
    }

    /**
     * Returns the AtomNode.
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getAtomNode() {
        return atomNode;
    }

    /**
     * Returns the MSNode at the given index.
     *
     * @param index a int.
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getAtomNode(int index) {
        return getAtomNodeList().get(index);
    }

    /**
     * Returns the AtomNode specified by the String n.
     *
     * @param n a {@link java.lang.String} object.
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getAtomNode(String n) {
        ArrayList<MSNode> list = getAtomNodeList();
        for (MSNode msNode : list) {
            if (msNode.getName().compareTo(n) == 0) {
                return msNode;
            }
        }
        return null;
    }

    /**
     * <p>getAtomByName.</p>
     *
     * @param n               a {@link java.lang.String} object.
     * @param caseInsensitive a boolean.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom getAtomByName(String n, boolean caseInsensitive) {
        ArrayList<MSNode> list = getAtomNodeList();
        for (MSNode msNode : list) {
            if (msNode.getName().equals(n)
                    || (caseInsensitive && msNode.getName().equalsIgnoreCase(n))) {
                return (Atom) msNode;
            }
        }
        return null;
    }

    /**
     * Returns an ArrayList of the AtomNode's children.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList<MSNode> getAtomNodeList() {
        return atomNode.getChildList();
    }

    /**
     * Returns the Bond at the supplied index.
     *
     * @param index a int.
     * @return a {@link ffx.potential.bonded.Bond} object.
     */
    public Bond getBond(int index) {
        return (Bond) bondNode.getChildAt(index);
    }

    /**
     * Returns the Bond with the given id.
     *
     * @param id a {@link java.lang.String} object.
     * @return a {@link ffx.potential.bonded.Bond} object.
     */
    public Bond getBond(String id) {
        int i = bondNode.getIndex(new Bond(id));
        if (i == -1) {
            return null;
        }
        return (Bond) bondNode.getChildAt(i);
    }

    /**
     * Returns the MultiScaleGroup's bonds FNode.
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getBonds() {
        return bondNode;
    }

    /**
     * Returns the MultiScaleGroup's center as a double[3].
     *
     * @return an array of double.
     */
    public double[] getCenter() {
        return center;
    }

    /**
     * Returns the MultiScaleGroup's dangelingatoms list.
     *
     * @return a {@link java.util.ArrayList} object.
     */
    public ArrayList getDangelingAtoms() {
        return dangelingatoms;
    }

    /**
     * Returns the MultiScaleGroup's Torsion MSNode.
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getTorsions() {
        return torsionNode;
    }

    /**
     * This method finds the Geometrical center of this MultiScaleGroup, or the
     * atomicWeight-weighted center if w is set to true, and returns it as a
     * double[3].
     *
     * @param w a boolean.
     * @return an array of double.
     */
    public double[] getMultiScaleCenter(boolean w) {
        // Find the center of atomicWeight if w == true, the center of geometry
        // if w ==
        // false
        double[] Rc = {0.0d, 0.0d, 0.0d};
        ArrayList<Atom> atoms = getAtomList();
        if (atoms == null) {
            return Rc;
        }
        double sum = 0.0d;
        if (w) {
            for (Atom a : atoms) {
                double mass = a.getMass();
                Rc[0] += mass * a.getX();
                Rc[1] += mass * a.getY();
                Rc[2] += mass * a.getZ();
                sum += mass;
            }
        } else {
            for (Atom a : atoms) {
                Rc[0] += a.getX();
                Rc[1] += a.getY();
                Rc[2] += a.getZ();
            }
            sum = atoms.size();
        }
        Rc[0] /= sum;
        Rc[1] /= sum;
        Rc[2] /= sum;
        return Rc;
    }

    /**
     * Returns the MultiScaleGroup's terms FNode.
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    public MSNode getTermNode() {
        return termNode;
    }

    /**
     * Returns true if the MultiScaleGroup is finalized.
     *
     * @return a boolean.
     */
    public boolean isFinalized() {
        return finalized;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Prints the MultiScaleGroup's Atoms and Bonds.
     */
    @Override
    public void print() {
        super.print();
        for (Atom a : atomNode.getAtomList()) {
            a.print();
        }
        for (ROLS m : bondNode.getBondList()) {
            Bond b = (Bond) m;
            b.print();
        }
    }

    /**
     * <p>
     * removeLeaves</p>
     */
    protected void removeLeaves() {
        if (termNode.getParent() == null) {
            return;
        }
        if (bondNode.getChildCount() == 0 && !(bondNode.getParent() == null)) {
            termNode.removeChild(bondNode);
        }
        if (angleNode.getChildCount() == 0 && !(angleNode.getParent() == null)) {
            termNode.removeChild(angleNode);
        }
        if (stretchBendNode.getChildCount() == 0 && !(stretchBendNode.getParent() == null)) {
            termNode.removeChild(stretchBendNode);
        }
        if (ureyBradleyNode.getChildCount() == 0 && !(ureyBradleyNode.getParent() == null)) {
            termNode.removeChild(ureyBradleyNode);
        }
        if (outOfPlaneBendNode.getChildCount() == 0 && (outOfPlaneBendNode.getParent() != null)) {
            termNode.removeChild(outOfPlaneBendNode);
        }
        if (torsionNode.getChildCount() == 0 && !(torsionNode.getParent() == null)) {
            termNode.removeChild(torsionNode);
        }
        if (stretchTorsionNode.getChildCount() == 0 && !(stretchTorsionNode.getParent() == null)) {
            termNode.removeChild(stretchTorsionNode);
        }
        if (angleTorsionNode.getChildCount() == 0 && !(angleTorsionNode.getParent() == null)) {
            termNode.removeChild(angleTorsionNode);
        }
        if (piOrbitalTorsionNode.getChildCount() == 0 && !(piOrbitalTorsionNode.getParent() == null)) {
            termNode.removeChild(piOrbitalTorsionNode);
        }
        if (torsionTorsionNode.getChildCount() == 0 && !(torsionTorsionNode.getParent() == null)) {
            termNode.removeChild(torsionTorsionNode);
        }
        if (improperTorsionNode.getChildCount() == 0 && !(improperTorsionNode.getParent() == null)) {
            termNode.removeChild(improperTorsionNode);
        }
        if (termNode.getChildCount() == 0) {
            removeChild(termNode);
        }
        if (atomNode.getChildCount() == 0) {
            removeChild(atomNode);
        }
    }

    /**
     * Sets the Angles node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setAngles(MSNode t) {
        termNode.removeChild(angleNode);
        angleNode = t;
        termNode.add(angleNode);
    }

    /**
     * Sets the ImproperTorsion node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setImproperTorsions(MSNode t) {
        termNode.removeChild(improperTorsionNode);
        improperTorsionNode = t;
        termNode.add(improperTorsionNode);
    }

    /**
     * Sets the Moieties node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setAtomNode(MSNode t) {
        remove(atomNode);
        atomNode = t;
        add(atomNode);
    }

    /**
     * Sets the Bonds node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setBonds(MSNode t) {
        termNode.removeChild(bondNode);
        bondNode.removeAllChildren();
        bondNode = t;
        termNode.add(bondNode);
    }

    /**
     * Sets the Stretch-Bends node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setStretchBends(MSNode t) {
        termNode.removeChild(stretchBendNode);
        stretchBendNode.removeAllChildren();
        stretchBendNode = t;
        termNode.add(stretchBendNode);
    }

    /**
     * Sets the Urey-Bradley node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setUreyBradleys(MSNode t) {
        termNode.removeChild(ureyBradleyNode);
        ureyBradleyNode.removeAllChildren();
        ureyBradleyNode = t;
        termNode.add(ureyBradleyNode);
    }

    /**
     * Sets the Out-of-Plane Bend node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setOutOfPlaneBends(MSNode t) {
        termNode.removeChild(outOfPlaneBendNode);
        outOfPlaneBendNode.removeAllChildren();
        outOfPlaneBendNode = t;
        termNode.add(outOfPlaneBendNode);
    }

    /**
     * Set the value of Center to d.
     *
     * @param d an array of double.
     */
    public void setCenter(double[] d) {
        center = d;
    }

    /**
     * Sets the MultiScaleGroup's dangelingatoms member to a.
     *
     * @param a a {@link java.util.ArrayList} object.
     */
    public void setDangelingAtoms(ArrayList<Atom> a) {
        dangelingatoms = a;
    }

    /**
     * Sets the MultiScaleGroup's torsion node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setTorsions(MSNode t) {
        termNode.removeChild(torsionNode);
        torsionNode = t;
        termNode.add(torsionNode);
    }

    /**
     * Sets the MultiScaleGroup's stretch-torsion node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setStretchTorsions(MSNode t) {
        termNode.removeChild(stretchTorsionNode);
        stretchTorsionNode = t;
        termNode.add(stretchTorsionNode);
    }

    /**
     * Sets the MultiScaleGroup's angle-torsion node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setAngleTorsions(MSNode t) {
        termNode.removeChild(angleTorsionNode);
        angleTorsionNode = t;
        termNode.add(angleTorsionNode);
    }

    /**
     * Sets the MultiScaleGroup's Pi-Orbital Torsion node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setPiOrbitalTorsions(MSNode t) {
        termNode.removeChild(piOrbitalTorsionNode);
        piOrbitalTorsionNode = t;
        termNode.add(piOrbitalTorsionNode);
    }

    /**
     * Sets the MultiScaleGroup's Torsion-Torsion node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setTorsionTorsions(MSNode t) {
        termNode.removeChild(torsionTorsionNode);
        torsionTorsionNode = t;
        termNode.add(torsionTorsionNode);
    }

    /**
     * Specifies whether the MultiScaleGroup has been finalized.
     *
     * @param t a boolean.
     */
    public void setFinalized(boolean t) {
        finalized = t;
    }

    /**
     * Sets the MultiScaleGroup's terms node to t.
     *
     * @param t a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setTerms(MSNode t) {
        remove(termNode);
        termNode = t;
        add(termNode);
    }

    /**
     * <p>
     * updateAtoms</p>
     */
    public void updateAtoms() {
        for (Atom a : getAtomList()) {
            a.update();
        }
    }

    /**
     * <p>
     * updateBonds</p>
     */
    public void updateBonds() {
        for (ROLS b : getBondList()) {
            b.update();
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setColor(RendererCache.ColorModel newColorModel, Color3f color, Material mat) {
        if (newColorModel == RendererCache.ColorModel.MOLECULE && (color == null || mat == null)) {
            return;
        }
        atomNode.setColor(newColorModel, color, mat);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setView(RendererCache.ViewModel newViewModel, List<BranchGroup> newShapes) {
        atomNode.setView(newViewModel, newShapes);
        bondNode.setView(newViewModel, newShapes);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns the MSGroup's name.
     */
    @Override
    public String toString() {
        return getName();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void update() {
        updateAtoms();
        updateBonds();
    }

}
