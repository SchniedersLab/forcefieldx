/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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

import javax.media.j3d.Bounds;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;

import java.awt.event.MouseEvent;

/**
 * The PickRotateBehavior class implements a mouse rotate behavior on a picked
 * object.
 *
 * @author Michael J. Schnieders
 *
 */
public class PickRotateBehavior extends PickMouseBehavior implements
        MouseBehaviorCallback {

    public MouseRotate drag;
    private PickingCallback callback = null;
    private TransformGroup currentTG;

    /**
     * <p>
     * Constructor for PickRotateBehavior.</p>
     *
     * @param bg a {@link javax.media.j3d.BranchGroup} object.
     * @param canvas a {@link javax.media.j3d.Canvas3D} object.
     * @param bounds a {@link javax.media.j3d.Bounds} object.
     * @param VPTG a {@link javax.media.j3d.TransformGroup} object.
     * @param pickMode a int.
     */
    public PickRotateBehavior(BranchGroup bg, Canvas3D canvas, Bounds bounds,
            TransformGroup VPTG, int pickMode) {
        super(canvas, bg, bounds);
        drag = new MouseRotate(MouseRotate.MANUAL_WAKEUP, VPTG);
        drag.setTransformGroup(currGrp);
        currGrp.addChild(drag);
        drag.setFactor(0.025);
        setSchedulingBounds(bounds);
        drag.setSchedulingBounds(bounds);
        pickCanvas.setMode(pickMode);
    }

    /*
     * Return the pickMode component of this PickRotateBehavior.
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
     * Sets the pickMode component of this PickRotateBehavior to the value of
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
            drag.setupCallback(null);
        } else {
            drag.setupCallback(this);
        }
    }

    /*
     * Callback method from MouseRotate This is used when the Picking callback
     * is enabled
     */
    /**
     * {@inheritDoc}
     */
    public void transformChanged(int type, Transform3D transform) {
        callback.transformChanged(PickingCallback.ROTATE, currentTG);
    }

    /**
     * {@inheritDoc}
     */
    public void transformClicked(int type, Transform3D transform) {
        callback.transformClicked(PickingCallback.ROTATE, currentTG);
    }

    /**
     * {@inheritDoc}
     */
    public void transformDoubleClicked(int type, Transform3D transform) {
        callback.transformDoubleClicked(PickingCallback.ROTATE, currentTG);
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
        if ((mevent.getModifiersEx() & MouseEvent.BUTTON1_DOWN_MASK) == MouseEvent.BUTTON1_DOWN_MASK) {
            pickCanvas.setShapeLocation(xpos, ypos);
            // PickResult r = pickCanvas.pickClosest();
            if (callback != null) {
                callback.transformChanged(PickingCallback.NO_PICK, null);
            }
        }
    }
}
/*
 * Copyright (c) 1996-1998 Sun MicroFSystems, Inc. All Rights Reserved. Sun
 * grants you ("Licensee") a non-exclusive, royalty free, license to use, modify
 * and redistribute this software in source and binary code form, provided that
 * i) this copyright notice and license appear on all copies of the software;
 * and ii) Licensee does not utilize the software in a manner which is
 * disparaging to Sun. This software is provided "AS IS," without a warranty of
 * any kind. ALL EXPRESS OR IMPLIED CONDITIONS, REPRESENTATIONS AND WARRANTIES,
 * INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
 * PURPOSE OR NON-INFRINGEMENT, ARE HEREBY EXCLUDED. SUN AND ITS LICENSORS SHALL
 * NOT BE LIABLE FOR ANY DAMAGES SUFFXRED BY LICENSEE AS A RESULT OF USING,
 * MODIFYING OR DISTRIBUTING THE SOFTWARE OR ITS DERIVATIVES. IN NO EVENT WILL
 * SUN OR ITS LICENSORS BE LIABLE FOR ANY LOST REVENUE, PROFIT OR DATA, OR FOR
 * DIRECT, INDIRECT, SPECIAL, CONSEQUENTIAL, INCIDENTAL OR PUNITIVE DAMAGES,
 * HOWEVER CAUSED AND REGARDLESS OF THE THEORY OF LIABILITY, ARISING OUT OF THE
 * USE OF OR INABILITY TO USE SOFTWARE, EVEN IF SUN HAS BEEN ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGES. This software is not designed or intended for
 * use in on-line control of aircraft, air traffic, aircraft navigation or
 * aircraft communications; or in the design, construction, operation or
 * maintenance of any nuclear facility. Licensee represents and warrants that it
 * will not use or redistribute the Software for such purposes.
 */
