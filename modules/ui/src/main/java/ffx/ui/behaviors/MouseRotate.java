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
package ffx.ui.behaviors;

import java.awt.AWTEvent;
import java.awt.event.MouseEvent;
import java.util.Iterator;

import org.jogamp.java3d.Behavior;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.WakeupCriterion;
import org.jogamp.java3d.WakeupOnAWTEvent;
import org.jogamp.vecmath.Matrix4d;
import org.jogamp.vecmath.Vector3d;

// import ffe.panels.*;

/**
 * The MouseRotate class implements a mouse rotation behavior.
 *
 * @author Michael J. Schnieders
 */
public class MouseRotate extends MouseBehavior {

    private static final Vector3d zero3d = new Vector3d(0.0, 0.0, 0.0);
    private static Transform3D VPTG_T3D = new Transform3D();
    private static Vector3d translation = new Vector3d();
    private static Matrix4d mat = new Matrix4d();
    double xAngle, yAngle;
    double xFactor = 0.001;
    double yFactor = 0.001;
    int doneID = 0;
    private MouseBehaviorCallback callback = null;

    /**
     * <p>
     * Constructor for MouseRotate.</p>
     *
     * @param flags a int.
     * @param VPTG  a {@link org.jogamp.java3d.TransformGroup} object.
     */
    public MouseRotate(int flags, TransformGroup VPTG) {
        super(flags, VPTG);
    }

    /**
     * <p>
     * Constructor for MouseRotate.</p>
     *
     * @param flags    a int.
     * @param VPTG     a {@link org.jogamp.java3d.TransformGroup} object.
     * @param behavior a {@link org.jogamp.java3d.Behavior} object.
     * @param postID   a int.
     * @param dID      a int.
     */
    public MouseRotate(int flags, TransformGroup VPTG, Behavior behavior, int postID, int dID) {
        super(flags, VPTG, behavior, postID);
        doneID = dID;
    }

    /**
     * Return the x-axis movement multipler.
     *
     * @return a double.
     */
    public double getXFactor() {
        return xFactor;
    }

    /**
     * Return the y-axis movement multipler.
     *
     * @return a double.
     */
    public double getYFactor() {
        return yFactor;
    }

    /**
     * <p>
     * initialize</p>
     */
    public void initialize() {
        super.initialize();
        xAngle = 0;
        yAngle = 0;
        if ((flags & INVERT_INPUT) == INVERT_INPUT) {
            invert = true;
            xFactor *= -1;
            yFactor *= -1;
        }
    }

    /**
     * {@inheritDoc}
     */
    public void processStimulus(Iterator<WakeupCriterion> criteria) {
        while (criteria.hasNext()) {
            WakeupCriterion wakeup = criteria.next();
            if (wakeup instanceof WakeupOnAWTEvent) {
                AWTEvent[] event = ((WakeupOnAWTEvent) wakeup).getAWTEvent();
                for (AWTEvent awtEvent : event) {
                    MouseEvent mevent = (MouseEvent) awtEvent;
                    processMouseEvent(mevent);
                    int id = awtEvent.getID();
                    // Drag and Button 1 down
                    if ((id == MouseEvent.MOUSE_DRAGGED)
                            && ((mevent.getModifiersEx() & MouseEvent.BUTTON1_DOWN_MASK) == MouseEvent.BUTTON1_DOWN_MASK)
                            && transformGroup != null) {
                        x = ((MouseEvent) awtEvent).getX();
                        y = ((MouseEvent) awtEvent).getY();
                        int dx = x - xLast;
                        int dy = y - yLast;
                        if (!reset) {
                            xAngle = dy * yFactor;
                            yAngle = dx * xFactor;
                            transformX.rotX(xAngle);
                            transformY.rotY(yAngle);
                            transformGroup.getTransform(currXform);
                            currXform.get(mat);
                            currXform.setTranslation(zero3d);
                            if (ViewerTG != null) {
                                ViewerTG.getTransform(VPTG_T3D);
                                VPTG_T3D.setTranslation(zero3d);
                                VPTG_T3D.invert();
                                currXform.mul(VPTG_T3D, currXform);
                            }
                            if (invert) {
                                currXform.mul(currXform, transformX);
                                currXform.mul(currXform, transformY);
                            } else {
                                currXform.mul(transformX, currXform);
                                currXform.mul(transformY, currXform);
                            }
                            if (ViewerTG != null) {
                                VPTG_T3D.invert();
                                currXform.mul(VPTG_T3D, currXform);
                            }
                            translation.set(mat.m03, mat.m13, mat.m23);
                            currXform.setTranslation(translation);
                            transformGroup.setTransform(currXform);
                            transformChanged(currXform);
                            if (callback != null) {
                                callback.transformChanged(MouseBehaviorCallback.TRANSLATE, currXform);
                            }
                        } else {
                            reset = false;
                        }
                        xLast = x;
                        yLast = y;
                    } else if (id == MouseEvent.MOUSE_PRESSED) {
                        xLast = ((MouseEvent) awtEvent).getX();
                        yLast = ((MouseEvent) awtEvent).getY();
                    } else if (id == MouseEvent.MOUSE_RELEASED) {
                        setTransformGroup(null);
                    }
                }
            }
        }
        if (transformGroup != null || postCriterion == null) {
            wakeupOn(mouseCriterion);
        } else {
            postId(doneID);
            wakeupOn(postCriterion);
        }
    }

    /**
     * Set the x-axis amd y-axis movement multipler with factor.
     *
     * @param factor a double.
     */
    public void setFactor(double factor) {
        xFactor = yFactor = factor;
    }

    /**
     * Set the x-axis amd y-axis movement multipler with xFactor and yFactor
     * respectively.
     *
     * @param xFactor a double.
     * @param yFactor a double.
     */
    public void setFactor(double xFactor, double yFactor) {
        this.xFactor = xFactor;
        this.yFactor = yFactor;
    }

    /**
     * <p>
     * setupCallback</p>
     *
     * @param c a {@link ffx.ui.behaviors.MouseBehaviorCallback} object.
     */
    public void setupCallback(MouseBehaviorCallback c) {
        callback = c;
    }

    /**
     * <p>
     * transformChanged</p>
     *
     * @param transform a {@link org.jogamp.java3d.Transform3D} object.
     */
    public void transformChanged(Transform3D transform) {
    }
}
