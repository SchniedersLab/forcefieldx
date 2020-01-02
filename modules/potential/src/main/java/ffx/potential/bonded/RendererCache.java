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
package ffx.potential.bonded;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;

import org.jogamp.java3d.Appearance;
import org.jogamp.java3d.BranchGroup;
import org.jogamp.java3d.Canvas3D;
import org.jogamp.java3d.ColoringAttributes;
import org.jogamp.java3d.Geometry;
import org.jogamp.java3d.GeometryArray;
import org.jogamp.java3d.LineAttributes;
import org.jogamp.java3d.Material;
import org.jogamp.java3d.Node;
import org.jogamp.java3d.PointAttributes;
import org.jogamp.java3d.PolygonAttributes;
import org.jogamp.java3d.RenderingAttributes;
import org.jogamp.java3d.ShaderAppearance;
import org.jogamp.java3d.ShaderProgram;
import org.jogamp.java3d.Shape3D;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.java3d.TransparencyAttributes;
import org.jogamp.java3d.utils.geometry.Cone;
import org.jogamp.java3d.utils.geometry.Cylinder;
import org.jogamp.java3d.utils.geometry.Sphere;
import org.jogamp.vecmath.Color3f;
import org.jogamp.vecmath.Point2d;
import org.jogamp.vecmath.Point3d;
import org.jogamp.vecmath.Vector3d;

/**
 * The RendererCache class defines constants related to rendering modes and
 * caches primitives for the Renderer.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class RendererCache {

    public enum ColorModel {

        CPK, GROUP, RESIDUE, POLYMER, MOLECULE, MONOCHROME, USERCOLOR, PARTIALCHARGE, PICK, SELECT, REVERT, STRUCTURE, APPLYUSERCOLOR
    }

    public enum ViewModel {

        WIREFRAME, BALLANDSTICK, SPACEFILL, RMIN, TUBE, INVISIBLE, RESTRICT, SHOWHYDROGENS, HIDEHYDROGENS, DETAIL, RIBBON, SHOWVRML, HIDEVRML, INDUCEDDIPOLE, FORCE, VELOCITY, ACCELERATION, HIDEVECTORS, UNIT, RELATIVE, ABSOLUTE, POINTS, LINES, FILL, DESTROY
    }

    /**
     * Constant <code>ORANGE</code>
     */
    static final Color3f ORANGE = new Color3f(Color.orange);
    /**
     * Constant <code>RED</code>
     */
    static final Color3f RED = new Color3f(Color.red);
    /**
     * Constant <code>BLUE</code>
     */
    static final Color3f BLUE = new Color3f(Color.blue);
    /**
     * Constant <code>GRAY</code>
     */
    static final Color3f GRAY = new Color3f(Color.lightGray);
    /**
     * Constant <code>YELLOW</code>
     */
    static final Color3f YELLOW = new Color3f(Color.yellow);
    /**
     * Constant <code>CYAN</code>
     */
    static final Color3f CYAN = new Color3f(Color.cyan);
    /**
     * Constant <code>GREEN</code>
     */
    static final Color3f GREEN = new Color3f(Color.green);
    /**
     * Constant <code>WHITE</code>
     */
    static final Color3f WHITE = new Color3f(Color.white);
    /**
     * Constant <code>PINK</code>
     */
    static final Color3f PINK = new Color3f(Color.pink);
    /**
     * Constant <code>MAGENTA</code>
     */
    static final Color3f MAGENTA = new Color3f(Color.magenta);
    /**
     * Constant <code>nullAp</code>
     */
    static final Appearance nullAp;
    /**
     * Constant <code>BLACK</code>
     */
    public static final Color3f BLACK = new Color3f(Color.black);
    /**
     * Constant <code>detail=3</code>
     */
    public static int detail = 3;
    /**
     * Constant <code>radius=1.0d</code>
     */
    public static double radius = 1.0d;
    /**
     * Constant <code>bondwidth=3</code>
     */
    public static int bondwidth = 3;
    /**
     * Constant <code>highlightSelections=false</code>
     */
    public static boolean highlightSelections = false;
    /**
     * Constant <code>labelAtoms=false</code>
     */
    public static boolean labelAtoms = false;
    /**
     * Constant <code>labelResidues=false</code>
     */
    public static boolean labelResidues = false;
    /**
     * Constant <code>pickingColor</code>
     */
    public static Color3f pickingColor = MAGENTA;
    /**
     * Constant <code>selectionColor</code>
     */
    public static Color3f selectionColor = YELLOW;
    /**
     * Constant <code>userColor</code>
     */
    public static Color3f userColor = WHITE;
    /**
     * Constant <code>viewModelHash</code>
     */
    public static final Hashtable<String, ViewModel> viewModelHash = new Hashtable<>();
    /**
     * Constant <code>colorModelHash</code>
     */
    public static final Hashtable<String, ColorModel> colorModelHash = new Hashtable<>();

    /**
     * <p>
     * Constructor for RendererCache.</p>
     */
    public RendererCache() {
    }

    /**
     * <p>
     * appearanceFactory</p>
     *
     * @param col         a {@link org.jogamp.vecmath.Color3f} object.
     * @param polygonType a {@link ffx.potential.bonded.RendererCache.ViewModel}
     *                    object.
     * @return a {@link org.jogamp.java3d.Appearance} object.
     */
    static Appearance appearanceFactory(Color3f col, ViewModel polygonType) {
        if (col == null) {
            return null;
        }
        Appearance ap;
        if (polygonType == RendererCache.ViewModel.FILL) {
            ap = fillAppearances.get(col);
        } else if (polygonType == RendererCache.ViewModel.POINTS) {
            ap = pointAppearances.get(col);
        } else {
            ap = lineAppearances.get(col);
        }
        if (ap == null) {
            ap = createAppearance(col, polygonType);
        }
        return ap;
    }

    /**
     * <p>
     * coneFactory</p>
     *
     * @param ap  a {@link org.jogamp.java3d.Appearance} object.
     * @param res a int.
     * @return a {@link org.jogamp.java3d.Shape3D} object.
     */
    protected static Shape3D coneFactory(Appearance ap, int res) {
        if (res > 3) {
            res = 3;
        }
        Shape3D cone = new Shape3D();
        cone.setAppearance(ap);
        cone.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        cone.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
        cone.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);
        cone.setCapability(Shape3D.ALLOW_PICKABLE_WRITE);
        cone.addGeometry(getConeGeom(0, res));
        cone.addGeometry(getConeGeom(1, res));
        return cone;
    }

    /**
     * <p>
     * getPolarGeom</p>
     *
     * @param res a int.
     * @return a {@link org.jogamp.java3d.Geometry} object.
     */
    protected static Geometry getPolarGeom(int res) {
        return getSphereGeom(res);
    }

    /**
     * <p>
     * doubleCylinderFactory</p>
     *
     * @param a1  a {@link ffx.potential.bonded.Atom} object.
     * @param a2  a {@link ffx.potential.bonded.Atom} object.
     * @param div a int.
     * @return a {@link org.jogamp.java3d.BranchGroup} object.
     */
    static BranchGroup doubleCylinderFactory(Atom a1, Atom a2, int div) {
        BranchGroup branchGroup;
        if (doubleCylinderPool.size() > 0) {
            branchGroup = doubleCylinderPool.remove(0);
            if (branchGroup != null) {
                TransformGroup cy1tg = (TransformGroup) branchGroup.getChild(0);
                Shape3D cy1 = (Shape3D) cy1tg.getChild(0);
                cy1.setAppearance(a1.getAtomAppearance());
                cy1.setUserData(a1);
                TransformGroup cy2tg = (TransformGroup) branchGroup.getChild(1);
                Shape3D cy2 = (Shape3D) cy2tg.getChild(0);
                cy2.setUserData(a2);
                cy2.setAppearance(a2.getAtomAppearance());
                return branchGroup;
            }
        }
        branchGroup = new BranchGroup();
        branchGroup.setCapability(BranchGroup.ALLOW_DETACH);
        branchGroup.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
        branchGroup.setCapability(BranchGroup.ENABLE_PICK_REPORTING);
        TransformGroup cy1tg = createTransformGroup(null);
        Shape3D cy1 = createCylinder(a1.getAtomAppearance(), div);
        cy1.setUserData(a1);
        cy1tg.addChild(cy1);
        branchGroup.addChild(cy1tg);
        TransformGroup cy2tg = createTransformGroup(null);
        Shape3D cy2 = createCylinder(a2.getAtomAppearance(), div);
        cy2.setUserData(a2);
        cy2tg.addChild(cy2);
        branchGroup.addChild(cy2tg);
        branchGroup.compile();
        return branchGroup;
    }

    /**
     * <p>
     * getColor</p>
     *
     * @param a    a {@link ffx.potential.bonded.Atom} object.
     * @param mode a {@link ffx.potential.bonded.RendererCache.ColorModel}
     *             object.
     * @return a {@link org.jogamp.vecmath.Color3f} object.
     */
    static Color3f getColor(Atom a, ColorModel mode) {
        switch (mode) {
            case CPK:
                return Atom.AtomColor.get(a.getAtomicNumber());
            case PICK:
                return pickingColor;
            case SELECT:
                return selectionColor;
            case PARTIALCHARGE:
                int index;
                double charge = a.getCharge();
                if (charge < 0.0d) {
                    float c = (float) (charge * 1000.0);
                    index = -1 * Math.round(c);
                    if (index > 999) {
                        index = 999;
                    }
                    if (negCharge[index] == null) {
                        float value = index * 0.001f;
                        negCharge[index] = new Color3f(1.0f, 1.0f - value,
                                1.0f - value);
                    }
                    return negCharge[index];
                } else if (charge == 0) {
                    return WHITE;
                } else {
                    float c = (float) (charge * 1000.0);
                    index = Math.round(c);
                    if (index > 999) {
                        index = 999;
                    }
                    if (posCharge[index] == null) {
                        float value = index * 0.001f;
                        posCharge[index] = new Color3f(1.0f - value, 1.0f - value,
                                1.0f);
                    }
                    return posCharge[index];
                }
            default:
                return NULLColor;
        }
    }

    /**
     * <p>
     * getCylinderGeom</p>
     *
     * @param num a int.
     * @param res a int.
     * @return a {@link org.jogamp.java3d.Geometry} object.
     */
    static Geometry getCylinderGeom(int num, int res) {
        if (res < 0) {
            res = 0;
        }
        if (res > 10) {
            res = 10;
        }
        if (cylgeom[num][res] == null) {
            initCylinderGeom(res);
        }
        return cylgeom[num][res];
    }

    /**
     * <p>
     * getScreenCoordinate</p>
     *
     * @param canvas  a {@link org.jogamp.java3d.Canvas3D} object.
     * @param node    a {@link org.jogamp.java3d.Node} object.
     * @param point3d a {@link org.jogamp.vecmath.Point3d} object.
     * @param point   a {@link org.jogamp.vecmath.Point2d} object.
     */
    static void getScreenCoordinate(Canvas3D canvas, Node node, Point3d point3d, final Point2d point) {
        if (point == null) {
            return;
        }
        // Get the transform to put the node in the world coordinate system
        node.getLocalToVworld(localToVworld);
        // Get the image plate transform
        canvas.getVworldToImagePlate(worldToImagePlate);
        // Transform into world coordinates
        localToVworld.transform(point3d);
        // Transform into imageplate coordinates
        worldToImagePlate.transform(point3d);
        // Final step to the 2D Screen.
        canvas.getPixelLocationFromImagePlate(point3d, point);
        /*
          Now we have the location where the point will be rendered on the
          screen depending on resize, placement, size, and eye point policies.
          This should only be called on points that reside within the clipping
          planes.
         */
    }

    /**
     * <p>
     * Getter for the field <code>sphereGeom</code>.</p>
     *
     * @param res a int.
     * @return a {@link org.jogamp.java3d.Geometry} object.
     */
    static Geometry getSphereGeom(int res) {
        if (res < 0) {
            res = 0;
        }
        if (res > 10) {
            res = 10;
        }
        if (sphereGeom[res] == null) {
            initSphereGeom(res);
        }
        return sphereGeom[res];
    }

    /**
     * <p>
     * materialFactory</p>
     *
     * @param col a {@link org.jogamp.vecmath.Color3f} object.
     * @return a {@link org.jogamp.java3d.Material} object.
     */
    static Material materialFactory(Color3f col) {
        if (col == null) {
            return null;
        }
        Material mat = materials.get(col);
        if (mat == null) {
            mat = new Material(col, BLACK, col, WHITE, 75.0f);
            mat.setLightingEnable(true);
            materials.put(col, mat);
        }
        return mat;
    }

    /**
     * <p>
     * poolDoubleCylinder</p>
     *
     * @param branchGroup a {@link org.jogamp.java3d.BranchGroup} object.
     */
    static void poolDoubleCylinder(BranchGroup branchGroup) {
        if (branchGroup != null) {
            doubleCylinderPool.add(branchGroup);
        }
    }

    /**
     * <p>
     * poolSphere</p>
     *
     * @param tg a {@link org.jogamp.java3d.BranchGroup} object.
     */
    static void poolSphere(BranchGroup tg) {
        if (tg != null) {
            spherePool.add(tg);
        }
    }

    /**
     * <p>
     * poolTransform3D</p>
     *
     * @param transform3D a {@link org.jogamp.java3d.Transform3D} object.
     */
    static void poolTransform3D(Transform3D transform3D) {
        if (transform3D != null) {
            transform3DPool.add(transform3D);
        }
    }

    /**
     * <p>
     * sphereFactory</p>
     *
     * @param ap          a {@link org.jogamp.java3d.Appearance} object.
     * @param div         a int.
     * @param transform3D a {@link org.jogamp.java3d.Transform3D} object.
     * @return a {@link org.jogamp.java3d.BranchGroup} object.
     */
    static BranchGroup sphereFactory(Appearance ap, int div, Transform3D transform3D) {
        BranchGroup branchGroup;
        if (spherePool.size() > 0) {
            branchGroup = spherePool.remove(0);
            if (branchGroup != null) {
                TransformGroup transformGroup = (TransformGroup) branchGroup.getChild(0);
                transformGroup.setTransform(transform3D);
                Shape3D sphere = (Shape3D) transformGroup.getChild(0);
                sphere.setAppearance(ap);
                return branchGroup;
            }
        }
        branchGroup = new BranchGroup();
        branchGroup.setCapability(BranchGroup.ALLOW_DETACH);
        branchGroup.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
        branchGroup.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
        TransformGroup transformGroup = createTransformGroup(transform3D);
        Shape3D sphere = createSphere(ap, div);
        transformGroup.addChild(sphere);
        branchGroup.addChild(transformGroup);
        branchGroup.compile();
        return branchGroup;
    }

    /**
     * <p>
     * toAtomColor</p>
     *
     * @param s a {@link java.lang.String} object.
     * @return a {@link org.jogamp.vecmath.Color3f} object.
     */
    static Color3f toAtomColor(String s) {
        String c = s.toLowerCase();
        if (c.startsWith("h")) {
            return WHITE;
        }
        if (c.startsWith("c")) {
            return GRAY;
        }
        if (c.startsWith("n")) {
            return BLUE;
        }
        if (c.startsWith("o")) {
            return RED;
        }
        if (c.startsWith("p")) {
            return GREEN;
        }
        if (c.startsWith("s")) {
            return YELLOW;
        }
        return NULLColor;
    }

    /**
     * <p>
     * transform3DFactory</p>
     *
     * @return a {@link org.jogamp.java3d.Transform3D} object.
     */
    static Transform3D transform3DFactory() {
        Transform3D transform3D;
        if (transform3DPool.size() > 0) {
            transform3D = transform3DPool.get(0);
            if (transform3D != null) {
                return transform3D;
            }
        }
        transform3D = new Transform3D();
        return transform3D;
    }

    /**
     * <p>
     * transform3DFactory</p>
     *
     * @param position a {@link org.jogamp.vecmath.Vector3d} object.
     * @param scale    a double.
     * @return a {@link org.jogamp.java3d.Transform3D} object.
     */
    static Transform3D transform3DFactory(Vector3d position, double scale) {
        Transform3D transform3D;
        if (transform3DPool.size() > 0) {
            transform3D = transform3DPool.get(0);
            if (transform3D != null) {
                transform3D.setTranslation(position);
                transform3D.setScale(scale);
                return transform3D;
            }
        }
        transform3D = new Transform3D();
        transform3D.setTranslation(position);
        transform3D.setScale(scale);
        return transform3D;
    }

    /**
     * <p>
     * getConeGeom</p>
     *
     * @param num a int.
     * @param res a int.
     * @return a {@link org.jogamp.java3d.Geometry} object.
     */
    private static Geometry getConeGeom(int num, int res) {
        if (res > 3) {
            res = 3;
        }
        if (conegeom[num][res] == null) {
            initConeGeom(res);
        }
        return conegeom[num][res];
    }

    private static void initConeGeom(int res) {
        Cone cone = new Cone(1.0f, 1.0f, Cone.GENERATE_NORMALS
                | Cone.ENABLE_GEOMETRY_PICKING | Cone.ENABLE_APPEARANCE_MODIFY,
                (res + 1) * 4, 1, nullAp);
        for (int i = 0; i < 2; i++) {
            conegeom[i][res] = cone.getShape(i).getGeometry();
            /*
             * conegeom[i][res].setCapability(Geometry.ALLOW_INTERSECT);
             * conegeom[i][res].setCapability(GeometryArray.ALLOW_FORMAT_READ);
             * conegeom[i][res].setCapability(GeometryArray.ALLOW_COUNT_READ);
             * conegeom
             * [i][res].setCapability(GeometryArray.ALLOW_COORDINATE_READ);
             */
        }
    }

    private static void initCylinderGeom(int res) {
        Appearance ap = new Appearance();
        Cylinder cyl = new Cylinder(1.0f, 1.0f, Cylinder.GENERATE_NORMALS
                | Cylinder.ENABLE_APPEARANCE_MODIFY
                | Cylinder.ENABLE_GEOMETRY_PICKING, 2 + res, 1, ap);
        for (int i = 0; i < 3; i++) {
            cylgeom[i][res] = cyl.getShape(i).getGeometry();
            try {
                cylgeom[i][res].setCapability(Geometry.ALLOW_INTERSECT);
                cylgeom[i][res].setCapability(GeometryArray.ALLOW_FORMAT_READ);
                cylgeom[i][res].setCapability(GeometryArray.ALLOW_COUNT_READ);
                cylgeom[i][res].setCapability(GeometryArray.ALLOW_COORDINATE_READ);
            } catch (Exception e) {
                return;
            }
        }
    }

    private static void initSphereGeom(int res) {
        Appearance ap = new Appearance();
        Sphere sphere;
        sphere = new Sphere(1.0f, Sphere.GENERATE_NORMALS
                | Sphere.ENABLE_APPEARANCE_MODIFY
                | Sphere.ENABLE_GEOMETRY_PICKING, 4 + 3 * res, ap);
        sphereGeom[res] = sphere.getShape().getGeometry();
        // GeometryArray g = (GeometryArray) sphereGeom[res];
        /*
         * if (!g.isLive()) { g.setCapability(g.ALLOW_FORMAT_READ);
         * g.setCapability(g.ALLOW_COUNT_READ);
         * g.setCapability(g.ALLOW_COORDINATE_READ); }
         */
    }

    /**
     * This method creates a Cylinder
     *
     * @param ap  a {@link org.jogamp.java3d.Appearance} object.
     * @param res a int.
     * @return a {@link org.jogamp.java3d.Shape3D} object.
     */
    private static Shape3D createCylinder(Appearance ap, int res) {
        if (res < 0) {
            res = 0;
        }
        if (res > 10) {
            res = 10;
        }
        final Shape3D cyl = new Shape3D();
        cyl.setAppearance(ap);
        cyl.addGeometry(getCylinderGeom(0, res));
        cyl.addGeometry(getCylinderGeom(1, res));
        cyl.addGeometry(getCylinderGeom(2, res));
        cyl.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        cyl.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
        cyl.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);
        cyl.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
        cyl.setCapability(Shape3D.ENABLE_PICK_REPORTING);
        cyl.setCapability(Shape3D.ALLOW_PICKABLE_WRITE);
        return cyl;
    }

    /**
     * This method creates a single Sphere from the given appearance
     */
    private static Shape3D createSphere(Appearance ap, int div) {
        Shape3D shape3d = new Shape3D();
        shape3d.setAppearance(ap);
        shape3d.addGeometry(getSphereGeom(div));
        shape3d.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
        shape3d.setCapability(Shape3D.ALLOW_GEOMETRY_WRITE);
        shape3d.setCapability(Shape3D.ALLOW_APPEARANCE_WRITE);
        shape3d.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
        shape3d.setCapability(Shape3D.ENABLE_PICK_REPORTING);
        shape3d.setCapability(Shape3D.ALLOW_PICKABLE_WRITE);
        return shape3d;
    }

    private static TransformGroup createTransformGroup(Transform3D transform3D) {
        TransformGroup transformGroup;
        if (transform3D == null) {
            transformGroup = new TransformGroup();
        } else {
            transformGroup = new TransformGroup(transform3D);
        }
        transformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        transformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
        transformGroup.setCapability(TransformGroup.ALLOW_CHILDREN_READ);
        return transformGroup;
    }

    private static Appearance createAppearance(Color3f col, ViewModel polygonType) {
        Appearance ap = null;
        if (shaderProgram != null) {
            ShaderAppearance sap = new ShaderAppearance();
            sap.setShaderProgram(shaderProgram);
            ap = sap;
        }
        if (ap == null) {
            ap = new Appearance();
        }
        Material mat = materialFactory(col);
        ap.setMaterial(mat);
        ap.setRenderingAttributes(renderingAttributes);
        ap.setColoringAttributes(coloringAttributes);
        ap.setLineAttributes(lineAttributes);
        ap.setPointAttributes(pointAttributes);
        if (polygonType == RendererCache.ViewModel.FILL) {
            ap.setPolygonAttributes(fillPolygonAttributes);
            fillAppearances.put(col, ap);
        } else if (polygonType == RendererCache.ViewModel.POINTS) {
            ap.setPolygonAttributes(pointPolygonAttributes);
            pointAppearances.put(col, ap);
        } else {
            ap.setPolygonAttributes(linePolygonAttributes);
            lineAppearances.put(col, ap);
        }
        return ap;
    }

    private static final Transform3D localToVworld = new Transform3D();
    private static final Transform3D worldToImagePlate = new Transform3D();
    private static final Hashtable<Color3f, Material> materials = new Hashtable<>();
    private static List<Transform3D> transform3DPool = Collections.synchronizedList(new ArrayList<>());
    private static List<BranchGroup> spherePool = Collections.synchronizedList(new ArrayList<>());
    private static List<BranchGroup> doubleCylinderPool = Collections.synchronizedList(new ArrayList<>());
    private static final Geometry[] sphereGeom = new Geometry[11];
    private static final Geometry[][] cylgeom = new Geometry[3][11];
    private static final Geometry[][] conegeom = new Geometry[2][4];
    private static final Hashtable<Color3f, Appearance> pointAppearances = new Hashtable<>();
    private static final Hashtable<Color3f, Appearance> lineAppearances = new Hashtable<>();
    private static final Hashtable<Color3f, Appearance> fillAppearances = new Hashtable<>();
    private static ShaderProgram shaderProgram = null;
    private static final Color3f[] negCharge = new Color3f[1000];
    private static final Color3f[] posCharge = new Color3f[1000];
    /**
     * Constant <code>NULLColor</code>
     */
    private static final Color3f NULLColor = new Color3f(Color.darkGray);
    /**
     * Constant <code>lineAttributes</code>
     */
    private static final LineAttributes lineAttributes = new LineAttributes();
    /**
     * Constant <code>pointAttributes</code>
     */
    private static final PointAttributes pointAttributes = new PointAttributes();
    /**
     * Constant <code>coloringAttributes</code>
     */
    private static final ColoringAttributes coloringAttributes = new ColoringAttributes();
    /**
     * Constant <code>renderingAttributes</code>
     */
    private static final RenderingAttributes renderingAttributes = new RenderingAttributes();
    /**
     * Constant <code>transparencyAttributes</code>
     */
    private static final TransparencyAttributes transparencyAttributes = new TransparencyAttributes();
    /**
     * Constant <code>fillPolygonAttributes</code>
     */
    private static final PolygonAttributes fillPolygonAttributes = new PolygonAttributes();
    /**
     * Constant <code>pointPolygonAttributes</code>
     */
    private static final PolygonAttributes pointPolygonAttributes = new PolygonAttributes();
    /**
     * Constant <code>linePolygonAttributes</code>
     */
    private static final PolygonAttributes linePolygonAttributes = new PolygonAttributes();

    static {
        coloringAttributes.setShadeModel(ColoringAttributes.NICEST);
        coloringAttributes.setColor(new Color3f(0, 0, 0));
        lineAttributes.setLineAntialiasingEnable(true);
        lineAttributes.setLinePattern(LineAttributes.PATTERN_SOLID);
        lineAttributes.setLineWidth(1.0f);
        pointAttributes.setPointAntialiasingEnable(true);
        pointAttributes.setPointSize(1.0f);
        fillPolygonAttributes.setPolygonMode(PolygonAttributes.POLYGON_FILL);
        fillPolygonAttributes.setCullFace(PolygonAttributes.CULL_BACK);
        linePolygonAttributes.setPolygonMode(PolygonAttributes.POLYGON_LINE);
        pointPolygonAttributes.setPolygonMode(PolygonAttributes.POLYGON_POINT);
        renderingAttributes.setVisible(true);
        renderingAttributes.setDepthBufferEnable(true);
        renderingAttributes.setDepthBufferWriteEnable(true);
        renderingAttributes.setIgnoreVertexColors(true);
        transparencyAttributes.setTransparencyMode(TransparencyAttributes.NONE);

        ViewModel[] values = ViewModel.values();
        for (ViewModel value : values) {
            viewModelHash.put(value.toString(), value);
        }

        ColorModel[] colorModelValues = ColorModel.values();
        for (ColorModel value : colorModelValues) {
            colorModelHash.put(value.toString(), value);
        }

        nullAp = new Appearance();
        RenderingAttributes ra = new RenderingAttributes();
        ra.setVisible(false);
        nullAp.setRenderingAttributes(ra);
    }
}
