//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

import org.jogamp.java3d.Bounds;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Canvas3D;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.utils.picking.PickResult;

/**
 * The PickZoomBehavior class implements a scaling behavior on a picked
 * scenegraph object.
 *
 * @author Michael J. Schnieders
 */
public class PickZoomBehavior extends PickMouseBehavior implements
        MouseBehaviorCallback {

    public MouseZoom zoom;
    private PickingCallback callback = null;
    private TransformGroup currentTG;

    /**
     * <p>
     * Constructor for PickZoomBehavior.</p>
     *
     * @param root     a {@link org.jogamp.java3d.BranchGroup} object.
     * @param canvas   a {@link org.jogamp.java3d.Canvas3D} object.
     * @param bounds   a {@link org.jogamp.java3d.Bounds} object.
     * @param VPTG     a {@link org.jogamp.java3d.TransformGroup} object.
     * @param pickMode a int.
     */
    public PickZoomBehavior(BranchGroup root, Canvas3D canvas, Bounds bounds,
                            TransformGroup VPTG, int pickMode) {
        super(canvas, root, bounds);
        zoom = new MouseZoom(MouseBehavior.MANUAL_WAKEUP, VPTG);
        currGrp.addChild(zoom);
        setSchedulingBounds(bounds);
        pickCanvas.setMode(pickMode);
    }

    /*
     * Return the pickMode component of this PickTranslateBehavior.
     */

    /**
     * <p>
     * getPickMode</p>
     *
     * @return a int.
     */
    public int getPickMode() {
        return pickCanvas.getMode();
    }

    /*
     * Sets the pickMode component of this PickTranslateBehavior to the value of
     * the passed pickMode. @param pickMode the pickMode to be copied.
     */

    /**
     * <p>
     * setPickMode</p>
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
     * <p>
     * setupCallback</p>
     *
     * @param c a {@link ffx.ui.behaviors.PickingCallback} object.
     */
    public void setupCallback(PickingCallback c) {
        callback = c;
        if (callback == null) {
            zoom.setupCallback(null);
        } else {
            zoom.setupCallback(this);
        }
    }

    /**
     * {@inheritDoc}
     */
    public void transformChanged(int type, Transform3D transform) {
        callback.transformChanged(PickingCallback.ZOOM, currentTG);
    }

    /**
     * {@inheritDoc}
     */
    public void transformClicked(int type, Transform3D transform) {
        callback.transformClicked(PickingCallback.ZOOM, currentTG);
    }

    /**
     * {@inheritDoc}
     */
    public void transformDoubleClicked(int type, Transform3D transform) {
        callback.transformDoubleClicked(PickingCallback.ZOOM, currentTG);
    }

    /*
     * Update the scene to manipulate any nodes. This is not meant to be called
     * by users. Behavior automatically calls this. You can call this only if
     * you know what you are doing.
     * @param xpos Current mouse X pos. @param ypos Current mouse Y pos.
     */

    /**
     * {@inheritDoc}
     */
    public void updateScene(int xpos, int ypos) {
        TransformGroup tg = null;
        if ((mevent.getModifiersEx() & MouseEvent.BUTTON3) == MouseEvent.BUTTON3) {
            pickCanvas.setShapeLocation(xpos, ypos);
            PickResult r = pickCanvas.pickClosest();
            if (r != null) {
                tg = (TransformGroup) r.getNode(PickResult.TRANSFORM_GROUP);
                if ((tg != null)
                        && (tg
                        .getCapability(TransformGroup.ALLOW_TRANSFORM_READ))
                        && (tg
                        .getCapability(TransformGroup.ALLOW_TRANSFORM_WRITE))) {
                    zoom.setTransformGroup(tg);
                    zoom.wakeup();
                    currentTG = tg;
                    if (callback != null) {
                        callback.transformClicked(PickingCallback.ZOOM,
                                currentTG);
                    }
                } else if (callback != null) {
                    callback.transformChanged(PickingCallback.NO_PICK, null);
                }
            }
        }
    }
}
