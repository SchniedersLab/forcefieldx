/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.ui;

import javax.media.j3d.Appearance;
import javax.media.j3d.Bounds;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.Font3D;
import javax.media.j3d.FontExtrusion;
import javax.media.j3d.Geometry;
import javax.media.j3d.Group;
import javax.media.j3d.Material;
import javax.media.j3d.Shape3D;
import javax.media.j3d.Text3D;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;
import javax.vecmath.AxisAngle4f;
import javax.vecmath.Color3f;
import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;

import java.awt.Color;
import java.awt.Font;
import java.util.Enumeration;

import ffx.ui.behaviors.MouseBehaviorCallback;

import com.sun.j3d.utils.geometry.Cone;
import com.sun.j3d.utils.geometry.Cylinder;
import com.sun.j3d.utils.geometry.Sphere;
import com.sun.j3d.utils.picking.PickTool;
import com.sun.j3d.utils.universe.ViewingPlatform;

/**
 * The GraphicsAxis class encapsulates the 3D Axis that is used to display and
 * control rotation/translation in the global frame.
 *
 * @author Michael J. Schnieders
 *
 */
public final class GraphicsAxis extends Group implements MouseBehaviorCallback {

    Canvas3D canvas;
    ViewingPlatform viewingPlatform;
    Bounds axisBounds;
    BranchGroup axisBranchGroup = new BranchGroup();
    TransformGroup axisTransformGroup = new TransformGroup();
    Transform3D axisTransform3D = new Transform3D();
    Vector3d axisVector3d = new Vector3d(-0.7, -0.6, -1.25);
    public Matrix3d matrix = new Matrix3d();

    /**
     * <p>
     * Constructor for GraphicsAxis.</p>
     *
     * @param v a {@link com.sun.j3d.utils.universe.ViewingPlatform} object.
     * @param c a {@link javax.media.j3d.Canvas3D} object.
     * @param b a {@link javax.media.j3d.Bounds} object.
     */
    public GraphicsAxis(ViewingPlatform v, Canvas3D c, Bounds b) {
        viewingPlatform = v;
        canvas = c;
        axisBounds = b;
        setCapability(Group.ENABLE_PICK_REPORTING);
        createAxis();
        setBounds(axisBounds);
        axisTransform3D.setTranslation(axisVector3d);
        axisTransform3D.setScale(0.015);
        axisTransformGroup.setTransform(axisTransform3D);
        axisTransformGroup.addChild(this);
        axisTransformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
        axisTransformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        axisTransformGroup.setCapability(TransformGroup.ENABLE_PICK_REPORTING);
        axisBranchGroup.setCapability(BranchGroup.ALLOW_DETACH);
        axisBranchGroup.setCapability(BranchGroup.ENABLE_PICK_REPORTING);
        axisBranchGroup.addChild(axisTransformGroup);
        axisBranchGroup.compile();
        viewingPlatform.getViewPlatformTransform().addChild(axisBranchGroup);
    }

    /**
     * <p>
     * center</p>
     */
    public void center() {
        axisTransform3D.setIdentity();
        axisTransform3D.setScale(0.015);
        axisTransform3D.setTranslation(axisVector3d);
        axisTransformGroup.setTransform(axisTransform3D);
    }

    /**
     * <p>
     * createAxis</p>
     */
    public void createAxis() {
        Appearance ap = new Appearance();
        Color3f col = new Color3f(Color.lightGray);
        Color3f black = new Color3f(Color.black);
        Color3f white = new Color3f(Color.white);
        Material mat = new Material(col, black, col, white, 50.0f);
        mat.setLightingEnable(true);
        ap.setMaterial(mat);
        // X-Axis
        Cone xcone = new Cone(2.0f, 3.0f, ap);
        xcone.setUserData(this);
        Transform3D xconeT3d = new Transform3D();
        xconeT3d.setTranslation(new Vector3d(10.0f, 0.0f, 0.0f));
        xconeT3d.setRotation(new AxisAngle4f(0.0f, 0.0f, 1.0f, (float) Math.PI
                / -2.0f));
        TransformGroup xconeTG = new TransformGroup(xconeT3d);
        xconeTG.addChild(xcone);
        Cylinder xcylinder = new Cylinder(1.0f, 9.0f, ap);
        xcylinder.setUserData(this);
        Transform3D xcyT3d = new Transform3D();
        xcyT3d.setTranslation(new Vector3d(4.5, 0.0, 0.0));
        xcyT3d.setRotation(new AxisAngle4f(0.0f, 0.0f, 1.0f,
                (float) Math.PI / 2.0f));
        TransformGroup xcyTG = new TransformGroup(xcyT3d);
        xcyTG.addChild(xcylinder);
        setCapabilities(xcone, xcylinder);
        addChild(xconeTG);
        addChild(xcyTG);
        // Y-Axis
        Cone ycone = new Cone(2.0f, 3.0f, ap);
        ycone.setUserData(this);
        Transform3D yconeT3d = new Transform3D();
        yconeT3d.setTranslation(new Vector3d(0.0f, 10.0f, 0.0f));
        TransformGroup yconeTG = new TransformGroup(yconeT3d);
        yconeTG.addChild(ycone);
        Cylinder ycylinder = new Cylinder(1.0f, 9.0f, ap);
        ycylinder.setUserData(this);
        Transform3D ycyT3d = new Transform3D();
        ycyT3d.setTranslation(new Vector3d(0.0, 4.5, 0.0));
        TransformGroup ycyTG = new TransformGroup(ycyT3d);
        ycyTG.addChild(ycylinder);
        setCapabilities(ycone, ycylinder);
        addChild(yconeTG);
        addChild(ycyTG);
        // Z-Axis
        Cone zcone = new Cone(2.0f, 3.0f, ap);
        zcone.setUserData(this);
        Transform3D zconeT3d = new Transform3D();
        zconeT3d.setTranslation(new Vector3d(0.0f, 0.0f, 10.0f));
        zconeT3d.setRotation(new AxisAngle4f(1.0f, 0.0f, 0.0f,
                (float) Math.PI / 2.0f));
        TransformGroup zconeTG = new TransformGroup(zconeT3d);
        zconeTG.addChild(zcone);
        Cylinder zcylinder = new Cylinder(1.0f, 9.0f, ap);
        zcylinder.setUserData(this);
        Transform3D zcyT3d = new Transform3D();
        zcyT3d.setTranslation(new Vector3d(0.0, 0.0, 4.5));
        zcyT3d.setRotation(new AxisAngle4f(1.0f, 0.0f, 0.0f,
                (float) Math.PI / 2.0f));
        TransformGroup zcyTG = new TransformGroup(zcyT3d);
        zcyTG.addChild(zcylinder);
        setCapabilities(zcone, zcylinder);
        addChild(zconeTG);
        addChild(zcyTG);
        Sphere sphere = new Sphere(1.0f, ap);
        if (!sphere.getShape().getGeometry(0).isCompiled()
                && !sphere.getShape().getGeometry(0).isLive()) {
            PickTool.setCapabilities(sphere.getShape(),
                    PickTool.INTERSECT_COORD);
        }
        addChild(sphere);
        // Labels
        ap = new Appearance();
        col = new Color3f(Color.green);
        mat = new Material(col, black, col, white, 50.0f);
        mat.setLightingEnable(true);
        ap.setMaterial(mat);
        Font font = new Font("Arial", Font.PLAIN, 4);
        Font3D font3d = new Font3D(font, new FontExtrusion());
        addChild(createAxisLabel("X", font3d, ap, 11.0, 0.0, 0.0));
        addChild(createAxisLabel("Y", font3d, ap, 0.0, 11.0, 0.0));
        addChild(createAxisLabel("Z", font3d, ap, 0.0, 0.0, 11.0));
    }

    @SuppressWarnings("unchecked")
    private TransformGroup createAxisLabel(String letter, Font3D font3d,
            Appearance ap, double x, double y, double z) {
        Text3D text = new Text3D(font3d, letter);
        text.setUserData(this);
        Transform3D t3D = new Transform3D();
        t3D.setTranslation(new Vector3d(x, y, z));
        TransformGroup tg = new TransformGroup(t3D);
        Shape3D text3d = new Shape3D(text, ap);
        text3d.setUserData(this);
        for (Enumeration<Geometry> e = text3d.getAllGeometries(); e.hasMoreElements();) {
            Geometry g = e.nextElement();
            g.setCapability(Geometry.ALLOW_INTERSECT);
        }
        text3d.setCapability(Shape3D.ENABLE_PICK_REPORTING);
        text3d.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        tg.addChild(text3d);
        return tg;
    }

    /**
     * <p>
     * Getter for the field <code>axisTransformGroup</code>.</p>
     *
     * @return a {@link javax.media.j3d.TransformGroup} object.
     */
    public TransformGroup getAxisTransformGroup() {
        return axisTransformGroup;
    }

    private void setCapabilities(Cone cone, Cylinder cylinder) {
        if (!cone.getShape(0).getGeometry(0).isLive()
                && !cone.getShape(0).getGeometry(0).isCompiled()) {
            PickTool.setCapabilities(cone.getShape(0), PickTool.INTERSECT_COORD);
            cone.getShape(0).setUserData(this);
            PickTool.setCapabilities(cone.getShape(1), PickTool.INTERSECT_COORD);
            cone.getShape(1).setUserData(this);
            PickTool.setCapabilities(cylinder.getShape(0),
                    PickTool.INTERSECT_COORD);
            cylinder.getShape(0).setUserData(this);
            PickTool.setCapabilities(cylinder.getShape(1),
                    PickTool.INTERSECT_COORD);
            cylinder.getShape(1).setUserData(this);
            PickTool.setCapabilities(cylinder.getShape(2),
                    PickTool.INTERSECT_COORD);
            cylinder.getShape(2).setUserData(this);
        }
    }

    /**
     * <p>
     * showAxis</p>
     *
     * @param b a boolean.
     */
    public void showAxis(boolean b) {
        if (b == true && !axisBranchGroup.isLive()) {
            viewingPlatform.getViewPlatformTransform().addChild(axisBranchGroup);
        } else if (b == false && axisBranchGroup.isLive()) {
            axisBranchGroup.detach();
        }
    }

    /**
     * {@inheritDoc}
     */
    public void transformChanged(int type, Transform3D viewTransform) {
        viewTransform.get(matrix);
        matrix.invert();
        viewTransform.set(matrix);
        axisTransform3D.set(viewTransform);
        axisTransform3D.setTranslation(axisVector3d);
        axisTransform3D.setScale(0.015);
        axisTransformGroup.setTransform(axisTransform3D);
    }

    /**
     * {@inheritDoc}
     */
    public void transformClicked(int type, Transform3D transform) {
        transformChanged(type, transform);
    }

    /**
     * {@inheritDoc}
     */
    public void transformDoubleClicked(int type, Transform3D transform) {
        transformChanged(type, transform);
    }
}
