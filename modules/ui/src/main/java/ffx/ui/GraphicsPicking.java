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
package ffx.ui;

import javax.swing.tree.TreePath;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.jogamp.java3d.Bounds;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Node;
import org.jogamp.java3d.SceneGraphPath;
import org.jogamp.java3d.Shape3D;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.utils.picking.PickCanvas;
import org.jogamp.java3d.utils.picking.PickIntersection;
import org.jogamp.java3d.utils.picking.PickResult;
import org.jogamp.vecmath.Vector3d;
import static org.apache.commons.math3.util.FastMath.toDegrees;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.RendererCache;
import ffx.potential.bonded.Residue;
import ffx.ui.behaviors.PickMouseBehavior;
import static ffx.numerics.math.VectorMath.bondAngle;
import static ffx.numerics.math.VectorMath.dihedralAngle;
import static ffx.numerics.math.VectorMath.dist;

/**
 * The GraphicsPicking class is used to make selections and measurements.
 *
 * @author Michael J. Schnieders
 */
public class GraphicsPicking extends PickMouseBehavior {

    private static final Logger logger = Logger.getLogger(GraphicsPicking.class.getName());

    public enum PickLevel {

        PICKATOM, PICKBOND, PICKANGLE, PICKDIHEDRAL, PICKRESIDUE, PICKMOLECULE,
        PICKPOLYMER, PICKSYSTEM, MEASUREDISTANCE, MEASUREANGLE, MEASUREDIHEDRAL
    }

    /**
     * Constant <code>pickLevelHash</code>
     */
    static final Hashtable<String, PickLevel> pickLevelHash = new Hashtable<>();

    static {
        PickLevel[] values = PickLevel.values();
        for (PickLevel value : values) {
            pickLevelHash.put(value.toString(), value);
        }
    }

    private MainPanel mainPanel;
    // Turn On/Off picking
    private boolean picking = false;
    // Picking Level
    private PickLevel pickLevel = PickLevel.PICKATOM;
    private PickLevel newPickLevel = PickLevel.PICKATOM;
    // Previously picked Atom
    private Atom previousAtom = null;
    // Number of times the previousAtom has been picked consecutively
    private int pickNumber = 0;
    // Previously picked MSNode
    private MSNode previousPick = null;
    // Selected Atoms for Measuring
    private Vector<Atom> atomCache = new Vector<>(4);
    private int count = 0;
    // A few static variables for reuse
    private double[] a = new double[3];
    private double[] b = new double[3];
    private double[] c = new double[3];
    private double[] d = new double[3];
    private Transform3D systemTransform3D = new Transform3D();
    private Vector3d syspos = new Vector3d();
    private Vector3d atpos = new Vector3d();

    /**
     * Constructor
     *
     * @param base           Base of the Scenegraph
     * @param bounds         Behavior bounds
     * @param graphicsCanvas Scene Canvas3D
     * @param mainPanel      MainPanel
     */
    public GraphicsPicking(BranchGroup base, Bounds bounds, GraphicsCanvas graphicsCanvas,
                           MainPanel mainPanel) {
        super(graphicsCanvas, base, bounds);
        this.mainPanel = mainPanel;
        pickCanvas.setMode(PickCanvas.GEOMETRY);
        pickCanvas.setTolerance(3.0f);
    }

    /**
     * Clear currently selected nodes
     */
    public void clear() {
        if (previousPick != null) {
            mainPanel.getHierarchy().collapsePath(new TreePath(previousPick.getPath()));
            previousPick.setSelected(false);
            previousPick.setColor(RendererCache.ColorModel.SELECT, null, null);
            previousPick = null;
            pickNumber = 0;
        }
        for (Atom a : atomCache) {
            a.setSelected(false);
            a.setColor(RendererCache.ColorModel.SELECT, null, null);
        }
        atomCache.clear();
    }

    private void distance(Atom atom, double[] pos) {
        MolecularAssembly m = (MolecularAssembly) atom
                .getMSNode(MolecularAssembly.class);
        m.getTransformGroup().getTransform(systemTransform3D);
        systemTransform3D.get(syspos);
        systemTransform3D.setScale(1.0d);
        systemTransform3D.setTranslation(new Vector3d(0, 0, 0));
        atom.getV3D(atpos);
        systemTransform3D.transform(atpos);
        atpos.add(syspos);
        atpos.get(pos);
    }

    /**
     * <p>
     * getPick</p>
     *
     * @return a {@link ffx.potential.bonded.MSNode} object.
     */
    MSNode getPick() {
        return previousPick;
    }

    /**
     * <p>
     * Getter for the field <code>picking</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getPicking() {
        return picking;
    }

    /**
     * <p>
     * Getter for the field <code>pickLevel</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    String getPickLevel() {
        return pickLevel.toString();
    }


    private void measure() {
        String measurement;
        double value;
        Atom a1, a2, a3, a4;
        switch (pickLevel) {
            case MEASUREDISTANCE:
                if (atomCache.size() < 2) {
                    return;
                }
                a1 = atomCache.get(0);
                a2 = atomCache.get(1);
                distance(a1, a);
                distance(a2, b);
                value = dist(a, b);
                measurement = "\nDistance\t" + a1.getIndex() + ", " + a2.getIndex()
                        + ":   \t" + format("%10.5f", value);
                break;
            case MEASUREANGLE:
                if (atomCache.size() < 3) {
                    return;
                }
                a1 = atomCache.get(0);
                a2 = atomCache.get(1);
                a3 = atomCache.get(2);
                distance(a1, a);
                distance(a2, b);
                distance(a3, c);
                value = bondAngle(a, b, c);
                value = toDegrees(value);
                measurement = "\nAngle\t" + a1.getIndex() + ", " + a2.getIndex() + ", "
                        + a3.getIndex() + ":   \t" + format("%10.5f", value);
                break;
            case MEASUREDIHEDRAL:
                if (atomCache.size() < 4) {
                    return;
                }
                a1 = atomCache.get(0);
                a2 = atomCache.get(1);
                a3 = atomCache.get(2);
                a4 = atomCache.get(3);
                distance(a1, a);
                distance(a2, b);
                distance(a3, c);
                distance(a4, d);
                value = dihedralAngle(a, b, c, d);
                value = toDegrees(value);
                measurement = "\nDihedral\t" + a1.getIndex() + ", " + a2.getIndex()
                        + ", " + a3.getIndex() + ", " + a4.getIndex() + ":\t"
                        + format("%10.5f", value);
                break;
            default:
                return;
        }
        logger.info(measurement);
        ModelingShell modelingShell = mainPanel.getModelingShell();
        modelingShell.setMeasurement(measurement, value);
        count = 0;
    }

    /**
     * <p>
     * resetCount</p>
     */
    void resetCount() {
        count = 0;
    }

    /**
     * <p>
     * Setter for the field <code>picking</code>.</p>
     *
     * @param m a boolean.
     */
    public void setPicking(boolean m) {
        picking = m;
        if (!picking) {
            clear();
        }
    }

    /**
     * <p>
     * Setter for the field <code>pickLevel</code>.</p>
     *
     * @param newPick a {@link java.lang.String} object.
     */
    void setPickLevel(String newPick) {
        if (pickLevelHash.containsKey(newPick.toUpperCase())) {
            newPickLevel = pickLevelHash.get(newPick.toUpperCase());
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Called by Java3D when an atom is picked
     */
    public void updateScene(int xpos, int ypos) {
        if (!picking) {
            return;
        }
        // Determine what FNode was picked
        pickCanvas.setShapeLocation(xpos, ypos);
        PickResult result = pickCanvas.pickClosest();
        if (result != null) {
            SceneGraphPath sceneGraphPath = result.getSceneGraphPath();
            Node node = sceneGraphPath.getObject();
            if (!(node instanceof Shape3D)) {
                return;
            }
            Shape3D pickedShape3D = (Shape3D) node;
            Object userData = pickedShape3D.getUserData();
            if (userData instanceof MolecularAssembly) {
                FFXSystem sys = (FFXSystem) userData;
                if (result.numIntersections() > 0) {
                    PickIntersection pickIntersection = result.getIntersection(0);
                    int[] coords = pickIntersection.getPrimitiveCoordinateIndices();
                    userData = sys.getAtomFromWireVertex(coords[0]);
                } else {
                    return;
                }
            }
            if (userData instanceof Atom) {
                Atom a = (Atom) userData;
                // Check to see if the pickLevel has changed
                if (!(pickLevel == newPickLevel)) {
                    pickLevel = newPickLevel;
                    pickNumber = 0;
                }
                // Clear selections between measurements
                String pickLevelString = pickLevel.toString();
                boolean measure = pickLevelString.startsWith("MEASURE");
                if (!measure || count == 0) {
                    for (Atom matom : atomCache) {
                        matom.setSelected(false);
                        matom.setColor(RendererCache.ColorModel.SELECT, null, null);
                    }
                    atomCache.clear();
                    count = 0;
                }
                // If measuring, select the current atom and add it to the cache
                if (measure && !atomCache.contains(a)) {
                    atomCache.add(0, a);
                    a.setSelected(true);
                    a.setColor(RendererCache.ColorModel.PICK, null, null);
                    count++;
                    measure();
                }
                if (!measure) {
                    // Check to see if the same Atom has been selected twice in a row.
                    // This allows iteration through the atom's terms.
                    if (a == previousAtom) {
                        pickNumber++;
                    } else {
                        previousAtom = a;
                        pickNumber = 0;
                    }
                    MSNode currentPick = null;
                    switch (pickLevel) {
                        case PICKATOM:
                            currentPick = a;
                            break;
                        case PICKBOND:
                        case PICKANGLE:
                        case PICKDIHEDRAL:
                            ArrayList terms;
                            if (pickLevel == PickLevel.PICKBOND) {
                                terms = a.getBonds();
                            } else if (pickLevel == PickLevel.PICKANGLE) {
                                terms = a.getAngles();
                            } else {
                                terms = a.getTorsions();
                            }
                            if (terms == null) {
                                return;
                            }
                            int num = terms.size();
                            if (pickNumber >= num) {
                                pickNumber = 0;
                            }
                            currentPick = (BondedTerm) terms.get(pickNumber);
                            break;
                        case PICKRESIDUE:
                        case PICKPOLYMER:
                        case PICKMOLECULE:
                        case PICKSYSTEM:
                            MSNode dataNode;
                            if (pickLevel == PickLevel.PICKRESIDUE) {
                                dataNode = (MSNode) a.getMSNode(Residue.class);
                            } else if (pickLevel == PickLevel.PICKPOLYMER) {
                                dataNode = (MSNode) a.getMSNode(Polymer.class);
                            } else if (pickLevel == PickLevel.PICKSYSTEM) {
                                dataNode = (MSNode) a.getMSNode(MolecularAssembly.class);
                            } else {
                                dataNode = (MSNode) a.getMSNode(Molecule.class);
                                if (dataNode == null) {
                                    dataNode = (MSNode) a.getMSNode(Polymer.class);
                                }
                            }
                            currentPick = dataNode;
                            break;
                        case MEASUREANGLE:
                        case MEASUREDIHEDRAL:
                        case MEASUREDISTANCE:
                            break;
                    }
                    // Add the selected node to the Tree View
                    if (currentPick != null) {
                        if (controlButton) {
                            mainPanel.getHierarchy().toggleSelection(currentPick);
                        } else if (currentPick != previousPick) {
                            mainPanel.getHierarchy().onlySelection(currentPick);
                        }
                        // Color the Current Pick by Picking Color
                        mainPanel.getGraphics3D().updateScene(currentPick,
                                false, false, null, true,
                                RendererCache.ColorModel.PICK);
                    }
                    // Remove picking color from the previousPick
                    if (previousPick != null && previousPick != currentPick) {
                        previousPick.setColor(RendererCache.ColorModel.REVERT, null, null);
                    }
                    previousPick = currentPick;
                }
            }
        }
    }
}
