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

/**
 * The MouseZoom class implements a Mouse Zoom behavior.
 *
 * @author Michael J. Schnieders
 */
public class MouseZoom extends MouseBehavior {

    private double zFactor = 0.0002;
    private MouseBehaviorCallback callback = null;
    private int mouseButton = MouseEvent.BUTTON2_DOWN_MASK;
    private int doneID = 0;

    /**
     * <p>
     * Constructor for MouseZoom.</p>
     *
     * @param flags a int.
     * @param VPTG  a {@link org.jogamp.java3d.TransformGroup} object.
     */
    public MouseZoom(int flags, TransformGroup VPTG) {
        super(flags, VPTG);
    }

    /**
     * <p>
     * Constructor for MouseZoom.</p>
     *
     * @param flags    a int.
     * @param VPTG     a {@link org.jogamp.java3d.TransformGroup} object.
     * @param behavior a {@link org.jogamp.java3d.Behavior} object.
     * @param postID   a int.
     * @param dID      a int.
     */
    public MouseZoom(int flags, TransformGroup VPTG, Behavior behavior, int postID, int dID) {
        super(flags, VPTG, behavior, postID);
        doneID = dID;
    }

    /**
     * Return the y-axis movement multipler.
     *
     * <p>
     * getFactor</p>
     *
     * @return a double.
     */
    public double getFactor() {
        return zFactor;
    }

    /**
     * <p>
     * initialize</p>
     */
    public void initialize() {
        super.initialize();
        if ((flags & INVERT_INPUT) == INVERT_INPUT) {
            zFactor *= -1;
            invert = true;
        }
    }

    /**
     * <p>
     * Setter for the field <code>mouseButton</code>.</p>
     *
     * @param button a int.
     */
    public void setMouseButton(int button) {
        mouseButton = button;
    }

    /**
     * {@inheritDoc}
     */
    public void processStimulus(Iterator<WakeupCriterion> criteria) {
        boolean done = false;
        while (criteria.hasNext()) {
            WakeupCriterion wakeup = criteria.next();
            if (wakeup instanceof WakeupOnAWTEvent) {
                AWTEvent[] event = ((WakeupOnAWTEvent) wakeup).getAWTEvent();
                for (AWTEvent awtEvent : event) {
                    processMouseEvent((MouseEvent) awtEvent);
                    int id = awtEvent.getID();
                    MouseEvent mevent = (MouseEvent) awtEvent;
                    int mod = mevent.getModifiersEx();
                    boolean middleButton = ((mod & mouseButton) == mouseButton);
                    if (!middleButton) {
                        middleButton = ((mod & MouseEvent.ALT_DOWN_MASK) == MouseEvent.ALT_DOWN_MASK);
                    }
                    if ((id == MouseEvent.MOUSE_DRAGGED) && middleButton) {
                        y = ((MouseEvent) awtEvent).getY();
                        int dy = y - yLast;
                        if (!reset) {
                            transformGroup.getTransform(currXform);
                            double z = (-1.0) * dy * zFactor;
                            double scale = currXform.getScale() + z;
                            if (scale > 0) {
                                currXform.setScale(scale);
                                transformGroup.setTransform(currXform);
                                transformChanged(currXform);
                            }
                            if (callback != null) {
                                callback.transformChanged(MouseBehaviorCallback.ZOOM, currXform);
                            }
                        } else {
                            reset = false;
                        }
                        xLast = x;
                        yLast = y;
                    }
                    if (id == MouseEvent.MOUSE_PRESSED) {
                        xLast = ((MouseEvent) awtEvent).getX();
                        yLast = ((MouseEvent) awtEvent).getY();
                    } else if (id == MouseEvent.MOUSE_RELEASED) {
                        done = true;
                    }
                }
            }
        }
        if (!done) {
            wakeupOn(mouseCriterion);
        } else {
            reset = true;
            mouseButton = MouseEvent.BUTTON2_DOWN_MASK;
            postId(doneID);
            wakeupOn(postCriterion);
        }
    }

    /**
     * Set the y-axis movement multipler with factor.
     *
     * <p>
     * setFactor</p>
     *
     * @param factor a double.
     */
    public void setFactor(double factor) {
        zFactor = factor;
    }

    /**
     * The transformChanged method in the callback class will be called every
     * time the transform is updated.
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
