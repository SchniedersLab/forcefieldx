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

import java.awt.AWTEvent;
import java.awt.event.MouseEvent;
import java.util.Enumeration;

import org.jogamp.java3d.Behavior;
import org.jogamp.java3d.Bounds;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Node;
import org.jogamp.java3d.SceneGraphPath;
import org.jogamp.java3d.Shape3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.WakeupCriterion;
import org.jogamp.java3d.WakeupOnAWTEvent;
import org.jogamp.java3d.WakeupOnBehaviorPost;
import org.jogamp.java3d.WakeupOr;
import org.jogamp.java3d.utils.picking.PickCanvas;
import org.jogamp.java3d.utils.picking.PickIntersection;
import org.jogamp.java3d.utils.picking.PickResult;
import org.jogamp.java3d.utils.universe.SimpleUniverse;
import org.jogamp.vecmath.Point3d;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.ui.GraphicsCanvas.LeftButtonMode;
import ffx.ui.behaviors.GlobalBehavior;
import ffx.ui.behaviors.MouseRotate;
import ffx.ui.behaviors.MouseTranslate;
import ffx.ui.behaviors.MouseZoom;

/**
 * The GraphicsEvents class listens for mouse events over the Java3D
 * GraphicsCanvas, dispatching work to more specialized System Rotation and
 * Translation Behaviors or to the GlobalOrbitBehavior.
 *
 * @author Michael J. Schnieders
 */
public class GraphicsEvents extends Behavior {

    // Behavior Post IDs
    /**
     * Constant <code>ROTATEPOST=1</code>
     */
    private static int ROTATEPOST = 1;
    /**
     * Constant <code>TRANSLATEPOST=2</code>
     */
    private static int TRANSLATEPOST = 2;
    /**
     * Constant <code>ZOOMPOST=3</code>
     */
    private static int ZOOMPOST = 3;
    /**
     * Constant <code>BEHAVIORDONEPOST=4</code>
     */
    private static int BEHAVIORDONEPOST = 4;

    // GUI Panels
    private MainPanel mainPanel;
    private GraphicsCanvas graphicsCanvas;
    // Wake up conditions
    private WakeupOr mouseCriterion;
    private WakeupOr postCriterion;
    // Mouse/Pick state upon wake up event
    private boolean buttonPress;
    private boolean leftButton;
    private boolean rightButton;
    private boolean middleButton;
    private PickCanvas pickCanvas;
    private PickResult pickResult;
    private Atom atom;
    private boolean axisSelected;
    // Behaviors
    private MouseRotate systemRotate;
    private MouseTranslate systemTranslate;
    private MouseZoom globalZoom;
    private GlobalBehavior viewOrbitBehavior;

    /**
     * <p>
     * Constructor for GraphicsEvents.</p>
     *
     * @param mainPanel      a {@link ffx.ui.MainPanel} object.
     * @param graphicsCanvas a {@link ffx.ui.GraphicsCanvas} object.
     * @param graphicsAxis   a {@link ffx.ui.GraphicsAxis} object.
     * @param simpleUniverse a {@link org.jogamp.java3d.utils.universe.SimpleUniverse} object.
     * @param bounds         a {@link org.jogamp.java3d.Bounds} object.
     * @param root           a {@link org.jogamp.java3d.BranchGroup} object.
     * @param transformGroup a {@link org.jogamp.java3d.TransformGroup} object.
     */
    public GraphicsEvents(MainPanel mainPanel, GraphicsCanvas graphicsCanvas, GraphicsAxis graphicsAxis,
                          SimpleUniverse simpleUniverse, Bounds bounds, BranchGroup root, TransformGroup transformGroup) {
        this.mainPanel = mainPanel;
        this.graphicsCanvas = graphicsCanvas;
        // Scenegraph Nodes
        TransformGroup viewTransformGroup = simpleUniverse.getViewingPlatform().getViewPlatformTransform();
        setSchedulingBounds(bounds);
        // Initialize the System Rotate Behavior
        systemRotate = new MouseRotate(MouseRotate.MANUAL_WAKEUP,
                viewTransformGroup, this, ROTATEPOST, BEHAVIORDONEPOST);
        systemRotate.setFactor(0.025);
        systemRotate.setSchedulingBounds(bounds);
        root.addChild(systemRotate);
        // Initialize the System Translate Behavior
        systemTranslate = new MouseTranslate(MouseTranslate.MANUAL_WAKEUP,
                viewTransformGroup, this, TRANSLATEPOST, BEHAVIORDONEPOST);
        systemTranslate.setFactor(0.5);
        systemTranslate.setSchedulingBounds(bounds);
        root.addChild(systemTranslate);
        // Initialize the globalZoom Behavior
        globalZoom = new MouseZoom(MouseZoom.MANUAL_WAKEUP, viewTransformGroup,
                this, ZOOMPOST, BEHAVIORDONEPOST);
        globalZoom.setFactor(0.0005);
        globalZoom.setSchedulingBounds(bounds);
        globalZoom.setTransformGroup(transformGroup);
        root.addChild(globalZoom);
        // Initialize the viewOrbitBehavior
        viewOrbitBehavior = new GlobalBehavior(this.graphicsCanvas);
        viewOrbitBehavior.setUpCallback(graphicsAxis);
        viewOrbitBehavior.setSchedulingBounds(bounds);
        simpleUniverse.getViewingPlatform().setViewPlatformBehavior(viewOrbitBehavior);
        // Initialize the PickCanvas
        pickCanvas = new PickCanvas(this.graphicsCanvas, simpleUniverse.getLocale());
        pickCanvas.setMode(PickCanvas.GEOMETRY);
        pickCanvas.setTolerance(20.0f);
    }

    /**
     * <p>
     * centerView</p>
     *
     * @param resetRotation    a boolean.
     * @param resetTranslation a boolean.
     * @param resetZoom        a boolean.
     */
    void centerView(boolean resetRotation, boolean resetTranslation, boolean resetZoom) {
        viewOrbitBehavior.centerView(resetRotation, resetTranslation);
    }

    private boolean globalZoom() {
        postId(ZOOMPOST);
        return true;
    }

    /**
     * <p>
     * initialize</p>
     */
    public void initialize() {
        WakeupCriterion[] behaviorPost = new WakeupCriterion[3];
        behaviorPost[0] = new WakeupOnBehaviorPost(systemRotate, BEHAVIORDONEPOST);
        behaviorPost[1] = new WakeupOnBehaviorPost(systemTranslate, BEHAVIORDONEPOST);
        behaviorPost[2] = new WakeupOnBehaviorPost(globalZoom, BEHAVIORDONEPOST);
        postCriterion = new WakeupOr(behaviorPost);
        WakeupCriterion[] awtCriterion = new WakeupCriterion[1];
        awtCriterion[0] = new WakeupOnAWTEvent(java.awt.AWTEvent.MOUSE_EVENT_MASK);
        mouseCriterion = new WakeupOr(awtCriterion);
        wakeupOn(mouseCriterion);
    }

    /**
     * <p>
     * processMouseEvent</p>
     *
     * @param evt a {@link java.awt.event.MouseEvent} object.
     */
    private void processMouseEvent(MouseEvent evt) {
        buttonPress = false;
        leftButton = false;
        middleButton = false;
        rightButton = false;
        int mod = evt.getModifiersEx();
        if (evt.getID() == MouseEvent.MOUSE_PRESSED) {
            buttonPress = true;
        }
        // Left Button
        if ((mod & MouseEvent.BUTTON1_DOWN_MASK) == MouseEvent.BUTTON1_DOWN_MASK) {
            leftButton = true;
        }
        // Middle Button
        if ((mod & MouseEvent.BUTTON2_DOWN_MASK) == MouseEvent.BUTTON2_DOWN_MASK) {
            middleButton = true;
        }
        // Alternatively, map "alt + button1" to the middle button
        if ((mod & MouseEvent.ALT_DOWN_MASK) == MouseEvent.ALT_DOWN_MASK) {
            if (leftButton) {
                middleButton = true;
                leftButton = false;
            }
        }
        // Right Button
        if ((mod & MouseEvent.BUTTON3_DOWN_MASK) == MouseEvent.BUTTON3_DOWN_MASK) {
            rightButton = true;
        }
        // Alternatively, map "shift + button1" to the right button
        if ((mod & MouseEvent.SHIFT_DOWN_MASK) == MouseEvent.SHIFT_DOWN_MASK) {
            if (leftButton) {
                rightButton = true;
                leftButton = false;
            }
        }
        int x = evt.getX();
        int y = evt.getY();
        atom = null;
        axisSelected = false;
        if (buttonPress) {
            // Picking Results
            pickCanvas.setShapeLocation(x, y);
            // Once in a while "pickClosest" throws an exception due to
            // not being able to invert a matrix??
            // Catch and ignore this until a fix is determined...
            try {
                pickResult = pickCanvas.pickClosest();
            } catch (Exception e) {
                pickResult = null;
            }
            if (pickResult != null) {
                SceneGraphPath sgp = pickResult.getSceneGraphPath();
                Node node = sgp.getObject();
                if (node instanceof Shape3D) {
                    Shape3D s = (Shape3D) node;
                    Object o = s.getUserData();
                    if (o instanceof MolecularAssembly) {
                        MolecularAssembly sys = (MolecularAssembly) o;
                        if (pickResult.numIntersections() > 0) {
                            PickIntersection pi = pickResult.getIntersection(0);
                            int[] coords = pi.getPrimitiveCoordinateIndices();
                            atom = sys.getAtomFromWireVertex(coords[0]);
                        }
                    } else if (o instanceof Atom) {
                        atom = (Atom) o;
                    } else if (o instanceof GraphicsAxis) {
                        axisSelected = true;
                    }
                }
            }
        }
    }

    /**
     * Most of the logic for mouse interaction with the Scenegraph is here.
     * <p>
     * {@inheritDoc}
     */
    public void processStimulus(Enumeration criteria) {
        viewOrbitBehavior.setEnable(false);
        AWTEvent[] awtEvents;
        while (criteria.hasMoreElements()) {
            WakeupCriterion wakeup = (WakeupCriterion) criteria.nextElement();
            if (wakeup instanceof WakeupOnAWTEvent) {
                awtEvents = ((WakeupOnAWTEvent) wakeup).getAWTEvent();
                if (awtEvents == null) {
                    continue;
                }
                for (AWTEvent awtEvent : awtEvents) {
                    MouseEvent mouseEvent;
                    if (awtEvent instanceof MouseEvent) {
                        mouseEvent = (MouseEvent) awtEvent;
                        processMouseEvent(mouseEvent);
                    } else {
                        continue;
                    }
                    if (!axisSelected) {
                        // Wake Up System Translate Behavior
                        if (rightButton && buttonPress) {
                            systemTranslate.setMouseButton(MouseEvent.BUTTON3_DOWN_MASK);
                            if (systemTranslate()) {
                                wakeupOn(postCriterion);
                                return;
                            }
                        }
                        // Wake Up Left Button Mode
                        if (leftButton && buttonPress) {
                            LeftButtonMode leftButtonMode = graphicsCanvas.getLeftButtonMode();
                            switch (leftButtonMode) {
                                case ROTATE:
                                    if (systemRotate()) {
                                        wakeupOn(postCriterion);
                                        return;
                                    }
                                    break;
                                case TRANSLATE:
                                    systemTranslate.setMouseButton(MouseEvent.BUTTON1_DOWN_MASK);
                                    if (systemTranslate()) {
                                        wakeupOn(postCriterion);
                                        return;
                                    }
                                    break;
                                case ZOOM:
                                    globalZoom.setMouseButton(MouseEvent.BUTTON1_DOWN_MASK);
                                    if (globalZoom()) {
                                        wakeupOn(postCriterion);
                                        return;
                                    }
                            }
                        }
                        // Wake up Global Zoom Behavior
                        if (middleButton && buttonPress) {
                            globalZoom.setMouseButton(MouseEvent.BUTTON2_DOWN_MASK);
                            if (globalZoom()) {
                                wakeupOn(postCriterion);
                                return;
                            }
                        }
                    } else {
                        viewOrbitBehavior.setEnable(true);
                        wakeupOn(mouseCriterion);
                        return;
                    }
                }
            }
        }
        wakeupOn(mouseCriterion);
    }

    /**
     * <p>
     * setGlobalCenter</p>
     *
     * @param d an array of double.
     */
    public void setGlobalCenter(double[] d) {
        Point3d point = new Point3d(d);
        viewOrbitBehavior.setRotationCenter(point);
    }

    private boolean systemRotate() {
        TransformGroup transformGroup = getTransformGroup();
        if (transformGroup != null) {
            systemRotate.setTransformGroup(transformGroup);
            postId(ROTATEPOST);
            return true;
        }
        return false;
    }

    private boolean systemTranslate() {
        TransformGroup transformGroup = getTransformGroup();
        if (transformGroup != null) {
            systemTranslate.setTransformGroup(transformGroup);
            postId(TRANSLATEPOST);
            return true;
        }
        return false;
    }

    private TransformGroup getTransformGroup() {
        TransformGroup transformGroup = null;
        GraphicsCanvas.MouseMode mouseMode = graphicsCanvas.getMouseMode();
        if ((mouseMode == GraphicsCanvas.MouseMode.SYSTEMBELOWMOUSE) && atom != null) {
            transformGroup = (TransformGroup) pickResult.getNode(PickResult.TRANSFORM_GROUP);
        } else if (mouseMode == GraphicsCanvas.MouseMode.ACTIVESYSTEM) {
            if (mainPanel.getHierarchy().getActive() != null) {
                transformGroup = mainPanel.getHierarchy().getActive().getTransformGroup();
            }
        }
        if (transformGroup != null) {
            // Make sure we can read and write the TransformGroup.
            if (!transformGroup.getCapability(TransformGroup.ALLOW_TRANSFORM_READ)
                    || !transformGroup.getCapability(TransformGroup.ALLOW_TRANSFORM_WRITE)) {
                transformGroup = null;
            }
        }
        return transformGroup;
    }
}
