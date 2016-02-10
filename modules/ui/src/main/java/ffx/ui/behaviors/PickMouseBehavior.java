/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.ui.behaviors;

import java.awt.AWTEvent;
import java.awt.Event;
import java.awt.event.MouseEvent;
import java.util.Enumeration;

import javax.media.j3d.Behavior;
import javax.media.j3d.Bounds;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.TransformGroup;
import javax.media.j3d.WakeupCriterion;
import javax.media.j3d.WakeupOnAWTEvent;
import javax.media.j3d.WakeupOr;

import com.sun.j3d.utils.picking.PickCanvas;

/**
 * The PickMouseBehavior class is the base class for mouse picking behaviors.
 *
 * @author Michael J. Schnieders
 *
 */
public abstract class PickMouseBehavior extends Behavior {

    static int count = 0;
    protected PickCanvas pickCanvas;
    protected WakeupCriterion[] conditions;
    protected WakeupOr wakeupCondition;
    protected boolean buttonPress = false;
    protected boolean shiftButton = false;
    protected boolean controlButton = false;
    protected TransformGroup currGrp;
    protected MouseEvent mevent;

    /*
     * Creates a PickMouseBehavior given current canvas, root of the tree to
     * operate on, and the bounds.
     */
    /**
     * <p>
     * Constructor for PickMouseBehavior.</p>
     *
     * @param canvas a {@link javax.media.j3d.Canvas3D} object.
     * @param root a {@link javax.media.j3d.BranchGroup} object.
     * @param bounds a {@link javax.media.j3d.Bounds} object.
     */
    public PickMouseBehavior(Canvas3D canvas, BranchGroup root, Bounds bounds) {
        super();
        setSchedulingBounds(bounds);
        currGrp = new TransformGroup();
        currGrp.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
        currGrp.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        root.addChild(currGrp);
        pickCanvas = new PickCanvas(canvas, root);
        pickCanvas.setMode(PickCanvas.GEOMETRY);
        pickCanvas.setTolerance(10.0f);
    }

    /**
     * <p>
     * initialize</p>
     */
    public void initialize() {
        conditions = new WakeupCriterion[1];
        conditions[0] = new WakeupOnAWTEvent(Event.MOUSE_DOWN);
        wakeupCondition = new WakeupOr(conditions);
        wakeupOn(wakeupCondition);
    }

    private void processMouseEvent(MouseEvent evt) {
        buttonPress = false;
        if (evt.getID() == MouseEvent.MOUSE_PRESSED
                | evt.getID() == MouseEvent.MOUSE_CLICKED
                | evt.getID() == MouseEvent.MOUSE_RELEASED) {
            buttonPress = true;
        }
        if (evt.isControlDown()) {
            controlButton = true;
        } else {
            controlButton = false;
        }
        if (evt.isShiftDown()) {
            shiftButton = true;
        } else {
            shiftButton = false;
        }
    }

    /**
     * {@inheritDoc}
     */
    public void processStimulus(Enumeration criteria) {
        WakeupCriterion wakeup;
        AWTEvent[] evt = null;
        int xpos = 0, ypos = 0;
        while (criteria.hasMoreElements()) {
            wakeup = (WakeupCriterion) criteria.nextElement();
            if (wakeup instanceof WakeupOnAWTEvent) {
                evt = ((WakeupOnAWTEvent) wakeup).getAWTEvent();
            }
        }
        if (evt[0] instanceof MouseEvent) {
            mevent = (MouseEvent) evt[0];
            processMouseEvent((MouseEvent) evt[0]);
            xpos = mevent.getPoint().x;
            ypos = mevent.getPoint().y;
        }
        if (buttonPress) {
            updateScene(xpos, ypos);
        }
        wakeupOn(wakeupCondition);
    }

    /**
     * <p>
     * setTolerance</p>
     *
     * @param tol a float.
     */
    public void setTolerance(float tol) {
        if (pickCanvas != null) {
            pickCanvas.setTolerance(tol);
        }
    }

    /*
     * Subclasses shall implement this update function
     */
    /**
     * <p>
     * updateScene</p>
     *
     * @param xpos a int.
     * @param ypos a int.
     */
    public abstract void updateScene(int xpos, int ypos);
}
