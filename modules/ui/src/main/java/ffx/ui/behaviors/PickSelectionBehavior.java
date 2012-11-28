/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */
package ffx.ui.behaviors;

import java.awt.event.MouseEvent;

import javax.media.j3d.*;

import com.sun.j3d.utils.picking.PickResult;

/**
 * The PickSelectionBehavior class implements a mouse based selections behavior.
 *
 * @author Michael J. Schnieders
 *
 */
public class PickSelectionBehavior extends PickMouseBehavior implements
        MouseBehaviorCallback {

    MouseSelection drag;
    private PickingCallback callback = null;
    private TransformGroup currentTG;

    /**
     * <p>Constructor for PickSelectionBehavior.</p>
     *
     * @param root a {@link javax.media.j3d.BranchGroup} object.
     * @param canvas a {@link javax.media.j3d.Canvas3D} object.
     * @param bounds a {@link javax.media.j3d.Bounds} object.
     * @param VPTG a {@link javax.media.j3d.TransformGroup} object.
     * @param pickMode a int.
     */
    public PickSelectionBehavior(BranchGroup root, Canvas3D canvas,
            Bounds bounds, TransformGroup VPTG, int pickMode) {
        super(canvas, root, bounds);
        drag = new MouseSelection(MouseSelection.MANUAL_WAKEUP, VPTG);
        drag.setTransformGroup(currGrp);
        currGrp.addChild(drag);
        drag.setSchedulingBounds(bounds);
        setSchedulingBounds(bounds);
        pickCanvas.setMode(pickMode);
    }

    /*
     * Return the pickMode component of this PickRotateBehavior.
     */
    /**
     * <p>getPickMode</p>
     *
     * @return a int.
     */
    public int getPickMode() {
        return pickCanvas.getMode();
    }

    /*
     * Sets the pickMode component of this PickRotateBehavior to the value of
     * the passed pickMode. @param pickMode the pickMode to be copied.
     */
    /**
     * <p>setPickMode</p>
     *
     * @param pickMode a int.
     */
    public void setPickMode(int pickMode) {
        pickCanvas.setMode(pickMode);
    }

    /*
     * Register the class @param callback to be called each time the picked
     * object moves
     */
    /**
     * <p>setupCallback</p>
     *
     * @param c a {@link ffx.ui.behaviors.PickingCallback} object.
     */
    public void setupCallback(PickingCallback c) {
        callback = c;
        if (callback == null) {
            drag.setupCallback(null);
        } else {
            drag.setupCallback(this);
        }
    }

    /**
     * {@inheritDoc}
     */
    public void transformChanged(int type, Transform3D transform) {
        callback.transformChanged(PickingCallback.SELECTION, currentTG);
    }

    /**
     * {@inheritDoc}
     */
    public void transformClicked(int type, Transform3D transform) {
        callback.transformClicked(PickingCallback.SELECTION, currentTG);
    }

    /**
     * {@inheritDoc}
     */
    public void transformDoubleClicked(int type, Transform3D transform) {
        callback.transformDoubleClicked(PickingCallback.SELECTION, currentTG);
    }

    /*
     * Update the scene to manipulate any nodes.
     * @param xpos Current mouse X pos. @param ypos Current mouse Y pos.
     */
    /**
     * {@inheritDoc}
     */
    public void updateScene(int xpos, int ypos) {
        TransformGroup tg = null;
        if ((mevent.getModifiersEx() & MouseEvent.BUTTON1) == MouseEvent.BUTTON1) {
            pickCanvas.setShapeLocation(xpos, ypos);
            PickResult r = pickCanvas.pickClosest();
            if (r != null) {
                tg = (TransformGroup) r.getNode(PickResult.TRANSFORM_GROUP);
                if ((tg != null)
                        && (tg
                        .getCapability(TransformGroup.ALLOW_TRANSFORM_READ))
                        && (tg
                        .getCapability(TransformGroup.ALLOW_TRANSFORM_WRITE))) {
                    drag.wakeup();
                    currentTG = tg;
                    if (callback != null) {
                        callback.transformClicked(PickingCallback.SELECTION,
                                currentTG);
                    }
                }
            } else if (callback != null) {
                callback.transformClicked(PickingCallback.NO_PICK, null);
            }
        }
    }
}
