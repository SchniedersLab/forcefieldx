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

import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.WakeupCriterion;
import org.jogamp.java3d.WakeupOnAWTEvent;

/**
 * The MouseProperties class is simple extension of MouseBehavior.
 *
 * @author Michael J. Schnieders
 */
public class MouseProperties extends MouseBehavior {

    double xAngle, yAngle;
    double xFactor = .03;
    double yFactor = .03;
    // private MouseBehaviorCallback callback = null;

    /**
     * <p>
     * Constructor for MouseProperties.</p>
     *
     * @param flags a int.
     * @param VPTG  a {@link org.jogamp.java3d.TransformGroup} object.
     */
    public MouseProperties(int flags, TransformGroup VPTG) {
        super(flags, VPTG);
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
                    processMouseEvent((MouseEvent) awtEvent);
                    if (((buttonPress) && ((flags & MANUAL_WAKEUP) == 0))
                            || ((wakeUp) && ((flags & MANUAL_WAKEUP) != 0))) {
                        int id = awtEvent.getID();
                        if ((id == MouseEvent.MOUSE_DRAGGED)) {
                            xLast = ((MouseEvent) awtEvent).getX();
                            yLast = ((MouseEvent) awtEvent).getY();
                        } else if (id == MouseEvent.MOUSE_PRESSED) {
                            xLast = ((MouseEvent) awtEvent).getX();
                            yLast = ((MouseEvent) awtEvent).getY();
                        }
                    }
                }
            }
        }
        wakeupOn(mouseCriterion);
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
     * The transformChanged method in the callback class will be called every
     * time the transform is updated
     *
     * @param c a {@link ffx.ui.behaviors.MouseBehaviorCallback} object.
     */
    public void setupCallback(MouseBehaviorCallback c) {
        // callback = c;
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
