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

import java.awt.event.MouseEvent;
import java.util.Enumeration;

import org.jogamp.java3d.Behavior;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.WakeupCriterion;
import org.jogamp.java3d.WakeupOnAWTEvent;
import org.jogamp.java3d.WakeupOnBehaviorPost;
import org.jogamp.java3d.WakeupOr;

/**
 * The MouseBehavior class is the Base class for all mouse manipulators.
 *
 * @author Michael J. Schnieders
 */
public abstract class MouseBehavior extends Behavior {
    /*
     * Set this flag if you want to manually wakeup the behavior.
     */

    /**
     * Constant <code>MANUAL_WAKEUP=0x1</code>
     */
    public static final int MANUAL_WAKEUP = 0x1;
    /*
     * Set this flag if you want to invert the inputs. This is useful when the
     * transform for the view platform is being changed instead of the transform
     * for the object.
     */
    /**
     * Constant <code>INVERT_INPUT=0x2</code>
     */
    public static final int INVERT_INPUT = 0x2;
    protected WakeupCriterion[] mouseEvents;
    protected WakeupOr mouseCriterion;
    protected Behavior poster;
    protected int id;
    protected WakeupOnBehaviorPost postCriterion;
    protected int x, y;
    protected int x_last, y_last;
    protected TransformGroup transformGroup;
    protected Transform3D transformX;
    protected Transform3D transformY;
    protected Transform3D currXform;
    protected boolean buttonPress = false;
    protected boolean reset = false;
    protected boolean invert = false;
    protected boolean wakeUp = false;
    protected int flags = 0;
    protected TransformGroup ViewerTG;
    /*
     * Swap a new transformGroup replacing the old one. This allows manipulators
     * to operate on different nodes.
     * @param transformGroup The *new* transform group to be manipulated.
     */
    Transform3D t3d = new Transform3D();

    /**
     * <p>
     * Constructor for MouseBehavior.</p>
     *
     * @param format a int.
     * @param VPTG   a {@link org.jogamp.java3d.TransformGroup} object.
     */
    public MouseBehavior(int format, TransformGroup VPTG) {
        super();
        flags = format;
        ViewerTG = VPTG;
        currXform = new Transform3D();
        transformX = new Transform3D();
        transformY = new Transform3D();
        reset = true;
    }

    /**
     * <p>
     * Constructor for MouseBehavior.</p>
     *
     * @param format a int.
     * @param VPTG   a {@link org.jogamp.java3d.TransformGroup} object.
     * @param b      a {@link org.jogamp.java3d.Behavior} object.
     * @param i      a int.
     */
    public MouseBehavior(int format, TransformGroup VPTG, Behavior b, int i) {
        this(format, VPTG);
        poster = b;
        id = i;
    }

    /*
     * Initializes the behavior.
     */

    /**
     * <p>
     * initialize</p>
     */
    public void initialize() {
        mouseEvents = new WakeupCriterion[3];
        mouseEvents[0] = new WakeupOnAWTEvent(MouseEvent.MOUSE_DRAGGED);
        mouseEvents[1] = new WakeupOnAWTEvent(MouseEvent.MOUSE_PRESSED);
        mouseEvents[2] = new WakeupOnAWTEvent(MouseEvent.MOUSE_RELEASED);
        mouseCriterion = new WakeupOr(mouseEvents);
        if (poster == null) {
            wakeupOn(mouseCriterion);
        } else {
            postCriterion = new WakeupOnBehaviorPost(poster, id);
            wakeupOn(postCriterion);
        }
        x = 0;
        y = 0;
        x_last = 0;
        y_last = 0;
    }

    /*
     * Handles mouse events
     */

    /**
     * <p>
     * processMouseEvent</p>
     *
     * @param evt a {@link java.awt.event.MouseEvent} object.
     */
    public void processMouseEvent(MouseEvent evt) {
        if (evt.getID() == MouseEvent.MOUSE_PRESSED) {
            buttonPress = true;
        } else if (evt.getID() == MouseEvent.MOUSE_RELEASED) {
            buttonPress = false;
            wakeUp = false;
        }
    }

    /**
     * {@inheritDoc}
     */
    public abstract void processStimulus(Enumeration criteria);

    /**
     * <p>
     * Setter for the field <code>transformGroup</code>.</p>
     *
     * @param t a {@link org.jogamp.java3d.TransformGroup} object.
     */
    public void setTransformGroup(TransformGroup t) {
        transformGroup = t;
        currXform = new Transform3D();
        transformX = new Transform3D();
        transformY = new Transform3D();
        reset = true;
    }

    /*
     * Manually wake up the behavior. If MANUAL_WAKEUP flag was set upon
     * creation, you must wake up this behavior each time it is handled.
     */

    /**
     * <p>
     * wakeup</p>
     */
    public void wakeup() {
        wakeUp = true;
    }
}
