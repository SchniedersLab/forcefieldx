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

import javax.media.j3d.Canvas3D;
import javax.media.j3d.Transform3D;
import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;

import com.sun.j3d.utils.behaviors.vp.OrbitBehavior;

/**
 * The GlobalBehavior class allows mouse control over camera position, adding a
 * few functions to the OrbitBehavior class.
 *
 * @author Michael J. Schnieders
 *
 */
public class GlobalBehavior extends OrbitBehavior {

    private Transform3D trans3D = new Transform3D();
    private Matrix3d homeQuat = new Matrix3d();
    private Vector3d homeTrans = new Vector3d(0.0, 0.0, 2.0);
    private Vector3d trans = new Vector3d();
    private MouseBehaviorCallback navigation = null;
    private boolean first = false;

    /**
     * <p>Constructor for GlobalBehavior.</p>
     */
    public GlobalBehavior() {
        super();
    }

    /**
     * <p>Constructor for GlobalBehavior.</p>
     *
     * @param canvas a {@link javax.media.j3d.Canvas3D} object.
     */
    public GlobalBehavior(Canvas3D canvas) {
        super(canvas, OrbitBehavior.MOUSE_MOTION_LISTENER);
        trans3D.setTranslation(homeTrans);
        setHomeTransform(trans3D);
        setReverseRotate(true);
        setReverseTranslate(true);
        setRotFactors(2.0, 2.0);
        setProportionalZoom(true);
        setEnable(false);
        homeQuat.setIdentity();
    }

    /**
     * <p>centerView</p>
     *
     * @param resetRotation a boolean.
     * @param resetTranslation a boolean.
     */
    public void centerView(boolean resetRotation, boolean resetTranslation) {
        if (!resetRotation && !resetTranslation) {
            return;
        }
        vp.getViewPlatformTransform().getTransform(trans3D);
        trans3D.get(trans);
        if (resetRotation) {
            trans3D.set(homeQuat);
            trans3D.setTranslation(homeTrans);
        }
        if (resetTranslation) {
            trans3D.setTranslation(homeTrans);
        }
        setHomeTransform(trans3D);
        goHome();
        if (resetRotation) {
            navigation.transformChanged(MouseBehaviorCallback.ORBIT, trans3D);
        }
    }

    /**
     * <p>integrateTransforms</p>
     */
    public void integrateTransforms() {
        // The "first" flag allows the mouse motion to be reset
        // (ie. dx = x - x_last where x_last is wrong initially)
        if (first) {
            vp.getViewPlatformTransform().getTransform(trans3D);
        }
        super.integrateTransforms();
        if (first) {
            first = false;
            vp.getViewPlatformTransform().setTransform(trans3D);
        }
        vp.getViewPlatformTransform().getTransform(trans3D);
        navigation.transformChanged(MouseBehaviorCallback.ORBIT, trans3D);
    }

    /**
     * {@inheritDoc}
     */
    public void setEnable(boolean b) {
        super.setEnable(b);
        if (b) {
            first = true;
        }
    }

    /**
     * <p>setUpCallback</p>
     *
     * @param m a {@link ffx.ui.behaviors.MouseBehaviorCallback} object.
     */
    public void setUpCallback(MouseBehaviorCallback m) {
        navigation = m;
    }
}
