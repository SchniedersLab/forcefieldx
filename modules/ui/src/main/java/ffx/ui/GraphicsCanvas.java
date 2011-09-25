/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.ui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsEnvironment;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.logging.Logger;
import java.util.prefs.Preferences;

import javax.imageio.ImageIO;
import javax.media.j3d.AmbientLight;
import javax.media.j3d.Background;
import javax.media.j3d.BoundingSphere;
import javax.media.j3d.Bounds;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.DirectionalLight;
import javax.media.j3d.GraphicsConfigTemplate3D;
import javax.media.j3d.GraphicsContext3D;
import javax.media.j3d.ImageComponent;
import javax.media.j3d.ImageComponent2D;
import javax.media.j3d.J3DGraphics2D;
import javax.media.j3d.Raster;
import javax.media.j3d.Transform3D;
import javax.media.j3d.TransformGroup;
import javax.media.j3d.View;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JColorChooser;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.vecmath.Color3f;
import javax.vecmath.Point3d;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3d;
import javax.vecmath.Vector3f;

import com.sun.j3d.utils.universe.SimpleUniverse;

import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.RendererCache;
import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;
import java.util.logging.Level;

/**
 * The GraphicsCanvas class provides a Canvas on which to render 3D Graphics.
 * The following display types are currently supported: Wireframe, Ball & Stick,
 * Spacefill/CPK, RMIN and Tube.
 *
 * @author schnied
 * @version $Id: $
 */
@SuppressWarnings("serial")
public class GraphicsCanvas extends Canvas3D implements ActionListener {

    private static final Logger logger = Logger.getLogger(GraphicsCanvas.class.getName());

    /**
     * The ImageFormat enum lists supported image formats.
     */
    public enum ImageFormat {

        BMP, GIF, JPEG, PNG, WBMP;
    }

    /**
     * The MouseMode enum describes what system is affected by mouse drags.
     */
    public enum MouseMode {

        SYSTEMBELOWMOUSE, ACTIVESYSTEM;
    }

    /**
     * The LeftButtonMode enum describes what the left mouse button does.
     */
    public enum LeftButtonMode {

        ROTATE, TRANSLATE, ZOOM;
    }
    /** Constant <code>imageFormatHash</code> */
    public static final Hashtable<String, ImageFormat> imageFormatHash = new Hashtable<String, ImageFormat>();

    static {
        ImageFormat values[] = ImageFormat.values();
        for (ImageFormat value : values) {
            imageFormatHash.put(value.toString(), value);
        }
    }
    // Controller Classes
    private ffx.potential.bonded.Renderer renderer;
    private GraphicsEvents graphicsEvents;
    private GraphicsPicking rendererPicking;
    private MainPanel mainPanel;
    private GraphicsAxis graphicsAxis;
    private GraphicsFullScreen fullScreenWindow;
    // 3D Universe Variables
    private SimpleUniverse universe;
    private Background background;
    private BranchGroup baseBranchGroup;
    private TransformGroup baseTransformGroup;
    private Transform3D baseTransform3D = new Transform3D();
    private Bounds bounds;
    // State Variables
    private GraphicsPrefs graphics3DPrefs = null;
    private MouseMode mouseMode = MouseMode.ACTIVESYSTEM;
    private ImageFormat imageFormat = ImageFormat.PNG;
    private LeftButtonMode leftButtonMode = LeftButtonMode.ROTATE;
    private boolean imageCapture = false;
    private File imageName;

    /**
     * The GraphicsCanvas constructor initializes the Java3D Universe and
     * Behaviors.
     *
     * @param config a {@link java.awt.GraphicsConfiguration} object.
     * @param f a {@link ffx.ui.MainPanel} object.
     */
    public GraphicsCanvas(GraphicsConfiguration config, MainPanel f) {
        super(config);
        mainPanel = f;
        // logger.info("CG: " +
        // isShadingLanguageSupported(Shader.SHADING_LANGUAGE_CG));
        // logger.info("GLSL: " +
        // isShadingLanguageSupported(Shader.SHADING_LANGUAGE_GLSL));
        initialize();
    }

    /**
     * <p>Constructor for GraphicsCanvas.</p>
     *
     * @param f a {@link ffx.ui.MainPanel} object.
     */
    public GraphicsCanvas(MainPanel f) {
        this(GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getBestConfiguration(
                new GraphicsConfigTemplate3D()), f);
    }

    /**
     * {@inheritDoc}
     *
     * Handles ActionEvents from the Selection, Display, Color, Options, and
     * Picking Menus.
     */
    @Override
    public void actionPerformed(ActionEvent evt) {
        String arg = evt.getActionCommand();
        // Selection Menu
        if (arg.equals("LabelSelectedAtoms")) {
            labelSelectedAtoms();
        } else if (arg.equals("LabelSelectedResidues")) {
            labelSelectedResidues();
        } else if (arg.equals("SetLabelFontColor")) {
            setLabelFontColor();
        } else if (arg.equals("SetLabelFontSize")) {
            setLabelFontSize();
            // Display Menu
        } else if (RendererCache.viewModelHash.containsKey(arg.toUpperCase())) {
            setViewModel(arg);
        } else if (arg.equals("Preferences")) {
            preferences();
            // Color Menu
        } else if (RendererCache.colorModelHash.containsKey(arg.toUpperCase())) {
            setColorModel(arg);
        } else if (arg.equals("SetSelectionColor")) {
            setSelectionColor();
        } else if (arg.equals("SetUserColor")) {
            setUserColor();
            // Options
        } else if (arg.equals("SystemBelowMouse")) {
            mouseMode = MouseMode.SYSTEMBELOWMOUSE;
        } else if (arg.equals("ActiveSystem")) {
            mouseMode = MouseMode.ACTIVESYSTEM;
        } else if (arg.equals("Rotate")) {
            leftButtonMode = LeftButtonMode.ROTATE;
        } else if (arg.equals("Translate")) {
            leftButtonMode = LeftButtonMode.TRANSLATE;
        } else if (arg.equals("Zoom")) {
            leftButtonMode = LeftButtonMode.ZOOM;
        } else if (arg.equals("ResetRotation")) {
            resetRotation();
        } else if (arg.equals("ResetTranslation")) {
            resetTranslation();
        } else if (arg.equals("ResetRotationAndTranslation")) {
            resetRotationAndTranslation();
        } else if (arg.equalsIgnoreCase("RotateAboutPick")) {
            rotateAboutPick();
        } else if (arg.equalsIgnoreCase("RotateAboutCenter")) {
            rotateAboutCenter();
        } else if (arg.equalsIgnoreCase("ResetGlobalTranslation")) {
            resetGlobalTranslation();
        } else if (arg.equalsIgnoreCase("ResetGlobalRotation")) {
            resetGlobalRotation();
        } else if (arg.equalsIgnoreCase("ResetGlobalZoom")) {
            resetGlobalZoom();
        } else if (arg.equalsIgnoreCase("ResetGlobalView")) {
            resetGlobalView();
        } else if (arg.equals("FullScreen")) {
            fullScreen();
        } else if (arg.equals("SetBackgroundColor")) {
            setBackgroundColor();
        } else if (arg.equalsIgnoreCase("ZoomIn")) {
            zoomIn();
        } else if (arg.equalsIgnoreCase("ZoomOut")) {
            zoomOut();
            // Picking Menu
        } else if (arg.equalsIgnoreCase("GraphicsPicking")) {
            graphicsPicking(evt);
        } else if (imageFormatHash.containsKey(arg.toUpperCase())) {
            setImageFormat(arg);
        } else if (arg.equals("CaptureGraphics")) {
            captureGraphics();
        } else if (GraphicsPicking.pickLevelHash.containsKey(arg.toUpperCase())) {
            setPickingLevel(arg);
        } else if (arg.equals("SetGraphicsPickingColor")) {
            setGraphicsPickingColor();
        } else {
            logger.warning("Graphics Menu command not found: "
                    + arg.toString());
        }
    }

    /**
     * This attaches a MolecularAssembly to the Scene BranchGroup.
     *
     * @param s
     *            MolecularAssembly to attach.
     */
    public void attachModel(MolecularAssembly s) {
        if (s == null) {
            return;
        }
        synchronized (this) {
            BranchGroup bg = s.getBranchGroup();
            resetGlobalView();
            baseBranchGroup.addChild(bg);
        }
    }

    private void captureGraphics() {
        MolecularAssembly active = mainPanel.getHierarchy().getActive();
        if (active == null) {
            return;
        }
        imageName = null;
        String name = active.getName();
        JFileChooser fileChooser = MainPanel.resetFileChooser();
        fileChooser.setAcceptAllFileFilterUsed(true);
        if (mainPanel.getHierarchy().getActive() != null) {
            imageName = mainPanel.getHierarchy().getActive().getFile();
        } else {
            imageName = null;
        }
        if (imageName != null) {
            if (name.indexOf(".") > 0) {
                name = name.substring(0, name.indexOf("."));
            }
            imageName = new File(imageName.getParentFile() + File.separator
                    + name + "." + imageFormat);
            fileChooser.setSelectedFile(imageName);
        }
        fileChooser.setDialogTitle("Select Name for Screen Capture " + "("
                + imageFormat + ")");
        fileChooser.setCurrentDirectory(MainPanel.getPWD());
        int result = fileChooser.showSaveDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
            imageName = fileChooser.getSelectedFile();
            mainPanel.setCWD(fileChooser.getCurrentDirectory());
            imageCapture = true;
            repaint();
        }
    }

    /**
     * <p>fullScreen</p>
     */
    public void fullScreen() {
        if (fullScreenWindow == null) {
            fullScreenWindow = new GraphicsFullScreen(mainPanel.getFrame(),
                    this);
        }
        fullScreenWindow.enterFullScreen();
    }

    /**
     * <p>Getter for the field <code>mouseMode</code>.</p>
     *
     * @return a {@link ffx.ui.GraphicsCanvas.MouseMode} object.
     */
    public MouseMode getMouseMode() {
        return mouseMode;
    }

    /**
     * <p>Getter for the field <code>leftButtonMode</code>.</p>
     *
     * @return a {@link ffx.ui.GraphicsCanvas.LeftButtonMode} object.
     */
    public LeftButtonMode getLeftButtonMode() {
        return leftButtonMode;
    }

    /**
     * <p>getNavigation</p>
     *
     * @return a {@link ffx.ui.GraphicsAxis} object.
     */
    public GraphicsAxis getNavigation() {
        return graphicsAxis;
    }

    /**
     * <p>getStatusBar</p>
     *
     * @return a {@link javax.swing.JLabel} object.
     */
    public JLabel getStatusBar() {
        return mainPanel.getStatusBar();
    }

    /**
     * <p>graphicsPicking</p>
     *
     * @param evt a {@link java.awt.event.ActionEvent} object.
     */
    public void graphicsPicking(ActionEvent evt) {
        if (evt.getSource() instanceof JButton) {
            MainMenu m = mainPanel.getMainMenu();
            boolean picking = m.getPicking();
            if (picking) {
                rendererPicking.clear();
                rendererPicking.setPicking(false);
                m.setPickBehavior(false);
            } else {
                rendererPicking.setPicking(true);
                m.setPickBehavior(true);
            }
        } else if (evt.getSource() instanceof JCheckBoxMenuItem) {
            JCheckBoxMenuItem jcbmi = (JCheckBoxMenuItem) evt.getSource();
            if (jcbmi.isSelected()) {
                rendererPicking.setPicking(true);
            } else {
                rendererPicking.setPicking(false);
            }
        }
    }

    /**
     * Initialization of the GraphisCanvas.
     */
    private void initialize() {
        setBackground(Color.black);
        universe = new SimpleUniverse(this);
        SimpleUniverse.setJ3DThreadPriority(Thread.MAX_PRIORITY);
        universe.getViewingPlatform().setNominalViewingTransform();
        // Create the Scene Root BranchGroup
        BranchGroup objRoot = new BranchGroup();
        baseTransformGroup = new TransformGroup();
        Transform3D t3d = new Transform3D();
        t3d.setScale(0.1d);
        baseTransformGroup.setTransform(t3d);
        baseTransformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
        baseTransformGroup.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
        // Set the Background
        Color3f bgColor = new Color3f(RendererCache.BLACK);
        background = new Background(bgColor);
        background.setCapability(Background.ALLOW_COLOR_READ);
        background.setCapability(Background.ALLOW_COLOR_WRITE);
        bounds = new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 2000.0);
        background.setApplicationBounds(bounds);
        // Create lights
        AmbientLight aLgt = new AmbientLight(new Color3f(Color.darkGray));
        aLgt.setInfluencingBounds(bounds);
        Vector3f dir = new Vector3f(0.0f, -1.0f, -1.0f);
        Color3f dLgtColor = new Color3f(Color.lightGray);
        DirectionalLight dLgt = new DirectionalLight(dLgtColor, dir);
        dLgt.setInfluencingBounds(bounds);
        dir = new Vector3f(0.0f, 1.0f, -1.0f);
        dLgtColor = new Color3f(0.1f, 0.1f, 0.1f);
        DirectionalLight dLgt2 = new DirectionalLight(dLgtColor, dir);
        dLgt2.setInfluencingBounds(bounds);
        // Create the Base of the Molecular Scene
        baseBranchGroup = new BranchGroup();
        baseBranchGroup.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
        baseBranchGroup.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
        baseBranchGroup.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
        baseBranchGroup.setCapability(BranchGroup.ALLOW_BOUNDS_READ);
        // Add children created above to the base TransformGroup
        baseTransformGroup.addChild(background);
        baseTransformGroup.addChild(baseBranchGroup);
        objRoot.addChild(baseTransformGroup);
        // Position the view platmform and add lights
        View v = universe.getViewer().getView();
        v.setProjectionPolicy(View.PARALLEL_PROJECTION);
        v.setFrontClipPolicy(View.VIRTUAL_EYE);
        v.setFrontClipDistance(1.0);
        v.setBackClipPolicy(View.VIRTUAL_EYE);
        v.setBackClipDistance(10.0);
        v.setTransparencySortingPolicy(View.TRANSPARENCY_SORT_NONE);
        Transform3D trans = new Transform3D();
        trans.set(new Vector3d(0.0d, 0.0d, 2.0d));
        TransformGroup vptg = universe.getViewingPlatform().getViewPlatformTransform();
        vptg.setTransform(trans);
        BranchGroup viewBranch = new BranchGroup();
        viewBranch.addChild(aLgt);
        viewBranch.addChild(dLgt);
        viewBranch.addChild(dLgt2);
        vptg.addChild(viewBranch);
        // Initialize Behaviors
        graphicsAxis = new GraphicsAxis(universe.getViewingPlatform(), this,
                bounds);
        graphicsEvents = new GraphicsEvents(mainPanel, this, graphicsAxis,
                universe, bounds, baseBranchGroup, baseTransformGroup);
        baseBranchGroup.addChild(graphicsEvents);
        rendererPicking = new GraphicsPicking(baseBranchGroup, bounds, this,
                mainPanel);
        baseBranchGroup.addChild(rendererPicking);
        renderer = new ffx.potential.bonded.Renderer(bounds, mainPanel.getStatusBar());
        baseBranchGroup.addChild(renderer);
        // Compile the Root BranchGroup and add it to the Universe
        objRoot.compile();
        universe.addBranchGraph(objRoot);
    }

    // protected void processMouseEvent(java.awt.event.MouseEvent event){
    // super.processMouseEvent(event);
    // logger.info("MouseEvent: " + event.toString());
    // }
    /**
     * <p>isCacheFull</p>
     *
     * @return a boolean.
     */
    public boolean isCacheFull() {
        return renderer.isCacheFull();
    }

    /**
     * <p>isSceneRendering</p>
     *
     * @return a boolean.
     */
    public boolean isSceneRendering() {
        return renderer.isArmed();
    }

    /**
     * <p>labelSelectedAtoms</p>
     */
    public void labelSelectedAtoms() {
        if (RendererCache.labelAtoms == true) {
            RendererCache.labelAtoms = false;
            getStatusBar().setText("  Atom Labeling Turned Off");
        } else {
            RendererCache.labelAtoms = true;
            getStatusBar().setText("  Atom Labeling Turned On");
        }
        repaint();
    }

    /**
     **********************************************************************
     */
    // Selection Commands
    public void labelSelectedResidues() {
        if (RendererCache.labelResidues == true) {
            RendererCache.labelResidues = false;
            getStatusBar().setText("  Residue Labeling Turned Off");
        } else {
            RendererCache.labelResidues = true;
            getStatusBar().setText("  Residue Labeling Turned On");
        }
        repaint();
    }

    /**
     * Load preferences from the user node.
     */
    public void loadPrefs() {
        String c = GraphicsCanvas.class.getName();
        RendererCache.bondwidth = prefs.getInt(c + ".bondwidth", 3);
        RendererCache.detail = prefs.getInt(c + ".detail", 3);
        RendererCache.radius = prefs.getDouble(c + ".radius", 1.0d);
        String s = prefs.get(c + ".mouse", MouseMode.ACTIVESYSTEM.name());
        if (s.equalsIgnoreCase("ACTIVESYSTEM") || s.equalsIgnoreCase("SYSTEMBELOWMOUSE")) {
            mouseMode = MouseMode.valueOf(s);
        } else {
            mouseMode = MouseMode.ACTIVESYSTEM;
        }
        mainPanel.getMainMenu().setMouseMode(mouseMode);
        RendererCache.highlightSelections = prefs.getBoolean(c + ".highlight", false);
        mainPanel.getMainMenu().setHighlighting(
                RendererCache.highlightSelections);
        String[] hlColor = prefs.get(c + ".highlightColor", "153 153 255").trim().split(" +");
        Color newColor = new Color(Integer.parseInt(hlColor[0]), Integer.parseInt(hlColor[1]), Integer.parseInt(hlColor[2]));
        RendererCache.selectionColor = new Color3f(newColor);
        RendererCache.labelAtoms = prefs.getBoolean(c + ".labelAtoms", false);
        mainPanel.getMainMenu().setAtomLabels(RendererCache.labelAtoms);
        RendererCache.labelResidues = prefs.getBoolean(c + ".labelResidues", false);
        mainPanel.getMainMenu().setResidueLabels(RendererCache.labelResidues);
        /*
         * int fontSize = prefs.getInt("Graphics_labelSize", 12); J3DGraphics2D
         * j2D = getOffscreenCanvas3D().getGraphics2D(); Font currentFont =
         * j2D.getFont(); Font newFont = new Font(currentFont.getName(),
         * currentFont.getStyle(), fontSize); j2D.setFont(newFont); String[]
         * fontColor = prefs.get("Graphics_labelColor", "255 255 255")
         * .trim().split(" +"); newColor = new
         * Color(Integer.parseInt(fontColor[0]), Integer
         * .parseInt(fontColor[1]), Integer.parseInt(fontColor[2]));
         * j2D.setPaint(newColor);
         */
        String[] pickColor = prefs.get(c + ".pickColor", "102 255 102").trim().split(" +");
        newColor = new Color(Integer.parseInt(pickColor[0]), Integer.parseInt(pickColor[1]), Integer.parseInt(pickColor[2]));
        RendererCache.pickingColor = new Color3f(newColor);
        String pickLevel = prefs.get(c + ".pickLevel", "PickAtom");
        mainPanel.getMainMenu().setPickLevel(pickLevel);
        boolean pickMode = prefs.getBoolean(c + ".picking", false);
        if (pickMode) {
            rendererPicking.setPicking(true);
        }
        mainPanel.getMainMenu().setPickBehavior(pickMode);
        String[] userColor = prefs.get(c + ".userColor", "255 255 255").trim().split(" +");
        newColor = new Color(Integer.parseInt(userColor[0]), Integer.parseInt(userColor[1]), Integer.parseInt(userColor[2]));
        RendererCache.userColor = new Color3f(newColor);
        /*
         * String[] bgColor = prefs.get("Graphics_backgroundColor", "0 0 0")
         * .trim().split(" +"); newColor = new
         * Color(Integer.parseInt(bgColor[0]), Integer .parseInt(bgColor[1]),
         * Integer.parseInt(bgColor[2])); background.setColor(new
         * Color3f(newColor));
         */
    }

    /** {@inheritDoc} */
    // The following three methods modify default Canvas3D methods.
    @Override
    public void paint(java.awt.Graphics g) {
        super.paint(g);
        Toolkit.getDefaultToolkit().sync();
    }

    /**
     * {@inheritDoc}
     *
     * Labels are drawn in postRender.
     */
    @Override
    public void postRender() {
        if (RendererCache.labelAtoms || RendererCache.labelResidues) {
            J3DGraphics2D g2D = getGraphics2D();
            synchronized (mainPanel.getHierarchy()) {
                ArrayList<MSNode> nodes = mainPanel.getHierarchy().getActiveNodes();
                if (nodes != null && nodes.size() > 0) {
                    for (MSNode node : nodes) {
                        MolecularAssembly sys = (MolecularAssembly) node.getMSNode(MolecularAssembly.class);
                        if (sys != null) {
                            node.drawLabel(this, g2D, sys.getWireFrame());
                        }
                    }
                } else {
                    return;
                }
            }
            g2D.flush(true);
        }
    }

    /**
     * {@inheritDoc}
     *
     * Image capture from the 3D Canvas is done in postSwap.
     */
    @Override
    public void postSwap() {
        if (!imageCapture || mainPanel.getHierarchy().getActive() == null) {
            return;
        }
        GraphicsContext3D ctx = getGraphicsContext3D();
        Rectangle rect = getBounds();
        BufferedImage img = new BufferedImage(rect.width, rect.height,
                BufferedImage.TYPE_INT_RGB);
        ImageComponent2D comp = new ImageComponent2D(ImageComponent.FORMAT_RGB,
                img);
        Raster ras = new Raster(new Point3f(-1.0f, -1.0f, -1.0f),
                Raster.RASTER_COLOR, 0, 0, rect.width, rect.height, comp, null);
        ctx.readRaster(ras);
        img = ras.getImage().getImage();
        try {
            if (!ImageIO.write(img, imageFormat.toString(), imageName)) {
                String os = System.getProperty("os.name");
                logger.severe("No image writer was found for "
                        + imageFormat.toString()
                        + " on "
                        + os
                        + ".\nThis Java Runtime Environment (JRE) fails to meet the javax.imageio specification.\n"
                        + "Shame on you " + os + "!\n"
                        + "Please try a different image format.\n");
                imageName.delete();
            } else {
                logger.info("" + imageName + " was captured.");
            }
        } catch (IOException e) {
            logger.warning("Image capture failed.\n" + e);
        }
        imageCapture = false;
    }

    /**
     * <p>preferences</p>
     */
    public void preferences() {
        if (graphics3DPrefs == null) {
            graphics3DPrefs = new GraphicsPrefs(mainPanel.getFrame(), mainPanel.getDataRoot());
            graphics3DPrefs.setModal(false);
        }
        graphics3DPrefs.setVisible(true);
        graphics3DPrefs.toFront();
    }

    /**
     * <p>resetGlobalRotation</p>
     */
    public void resetGlobalRotation() {
        graphicsEvents.centerView(true, false, false);
    }

    /**
     * <p>resetGlobalTranslation</p>
     */
    public void resetGlobalTranslation() {
        graphicsEvents.centerView(false, true, false);
    }

    /**
     * This functions centers the scene.
     */
    public void resetGlobalView() {
        double radius = mainPanel.getDataRoot().getExtent();
        Transform3D t3d = new Transform3D();
        t3d.setScale(1.0d / (1.2d * radius));
        baseTransformGroup.setTransform(t3d);
        graphicsEvents.centerView(true, true, true);
    }

    /**
     * <p>resetGlobalZoom</p>
     */
    public void resetGlobalZoom() {
        double radius = mainPanel.getDataRoot().getExtent();
        baseTransformGroup.getTransform(baseTransform3D);
        baseTransform3D.setScale(1.0d / (1.2d * radius));
        baseTransformGroup.setTransform(baseTransform3D);
    }

    /**
     **********************************************************************
     */
    // Options Commands
    public void resetRotation() {
        MolecularAssembly sys = mainPanel.getHierarchy().getActive();
        if (sys != null) {
            sys.centerView(true, false);
        }
    }

    /**
     * <p>resetRotationAndTranslation</p>
     */
    public void resetRotationAndTranslation() {
        MolecularAssembly sys = mainPanel.getHierarchy().getActive();
        if (sys != null) {
            sys.centerView(true, true);
        }
    }

    /**
     * <p>resetTranslation</p>
     */
    public void resetTranslation() {
        MolecularAssembly sys = mainPanel.getHierarchy().getActive();
        if (sys != null) {
            sys.centerView(false, true);
        }
    }

    /**
     * <p>rotateAboutCenter</p>
     */
    public void rotateAboutCenter() {
        MolecularAssembly sys = mainPanel.getHierarchy().getActive();
        double center[] = sys.getMultiScaleCenter(false);
        sys.rotateAbout(new Vector3d(center));
    }

    /**
     * <p>rotateAboutPick</p>
     */
    public void rotateAboutPick() {
        MSNode node = rendererPicking.getPick();
        if (node != null) {
            double center[] = node.getCenter(false);
            MolecularAssembly m = (MolecularAssembly) node.getMSNode(MolecularAssembly.class);
            m.rotateAbout(new Vector3d(center));
        }
    }
    /** ******************************************************************** */
    // The following three methods modidfy default Canvas3D methods.
    /**
     * Save preferences to the user node.
     */
    private static final Preferences prefs = Preferences.userNodeForPackage(GraphicsCanvas.class);

    /**
     * <p>savePrefs</p>
     */
    public void savePrefs() {
        String c = GraphicsCanvas.class.getName();
        prefs.putInt(c + ".bondwidth", RendererCache.bondwidth);
        prefs.putInt(c + ".detail", RendererCache.detail);
        prefs.putDouble(c + ".radius", RendererCache.radius);
        prefs.put(c + ".mouse", mouseMode.name());
        prefs.putBoolean(c + ".highlight", RendererCache.highlightSelections);
        Color col = RendererCache.selectionColor.get();
        prefs.put(c + ".highlightColor", "" + col.getRed() + " " + col.getGreen() + " " + col.getBlue());
        prefs.putBoolean(c + ".labelAtoms", RendererCache.labelAtoms);
        prefs.putBoolean(c + ".labelResidues", RendererCache.labelResidues);
        prefs.putInt(c + ".labelSize", getGraphics2D().getFont().getSize());
        Color fontColor = getGraphics2D().getColor();
        prefs.put(c + ".labelColor", "" + fontColor.getRed() + " "
                + fontColor.getGreen() + " " + fontColor.getBlue());
        col = RendererCache.pickingColor.get();
        prefs.put(c + ".pickColor", "" + col.getRed() + " " + col.getGreen() + " " + col.getBlue());
        prefs.putBoolean(c + ".picking", rendererPicking.getPicking());
        prefs.put(c + ".pickLevel", rendererPicking.getPickLevel());
        col = RendererCache.userColor.get();
        prefs.put(c + ".userColor", "" + col.getRed() + " " + col.getGreen() + " " + col.getBlue());
        Color3f temp = new Color3f();
        background.getColor(temp);
        col = temp.get();
        prefs.put(c + ".backgroundColor", "" + col.getRed() + " " + col.getGreen() + " " + col.getBlue());
    }

    /**
     * <p>selected</p>
     */
    public void selected() {
        validate();
        repaint();
    }

    /**
     * <p>setAxisShowing</p>
     *
     * @param b a boolean.
     */
    public void setAxisShowing(boolean b) {
        if (b && graphicsAxis == null) {
            graphicsAxis = new GraphicsAxis(universe.getViewingPlatform(),
                    this, bounds);
        } else if (graphicsAxis != null) {
            graphicsAxis.showAxis(b);
        }
    }

    /**
     * <p>setBackgroundColor</p>
     */
    public void setBackgroundColor() {
        Color3f col = new Color3f();
        background.getColor(col);
        Color newcolor = JColorChooser.showDialog(this,
                "Choose Background Color", col.get());
        if (newcolor != null && newcolor != col.get()) {
            background.setColor(new Color3f(newcolor));
        }
    }

    /**
     * <p>setCaptures</p>
     *
     * @param c a boolean.
     */
    public void setCaptures(boolean c) {
        imageCapture = c;
    }

    /**
     * <p>setColor</p>
     *
     * @param model a {@link java.lang.String} object.
     */
    public void setColor(String model) {
        setColorModel(model);
    }

    /** ********************************************************************** */
    // Color Commands
    /**
     * Operates on the Active nodes.
     *
     * @param model
     *            String
     */
    public void setColorModel(String model) {
        if (!RendererCache.colorModelHash.containsKey(model.toUpperCase())) {
            return;
        }
        RendererCache.ColorModel colorModel = RendererCache.colorModelHash.get(model.toUpperCase());
        ArrayList<MSNode> active = mainPanel.getHierarchy().getActiveNodes();
        if (active == null) {
            return;
        }
        renderer.arm(active, false, false, null, true, colorModel);
    }

    /**
     * Operates on the passed node.
     *
     * @param model
     *            String
     * @param node a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setColorModel(String model, MSNode node) {
        if (node == null) {
            return;
        }
        if (!RendererCache.colorModelHash.containsKey(model.toUpperCase())) {
            return;
        }
        RendererCache.ColorModel colorModel = RendererCache.colorModelHash.get(model.toUpperCase());
        renderer.arm(node, false, false, null, true, colorModel);
    }

    /**
     * <p>setGraphicsPickingColor</p>
     */
    public void setGraphicsPickingColor() {
        Color3f col = new Color3f();
        RendererCache.pickingColor.get(col);
        Color newcolor = JColorChooser.showDialog(this, "Choose Picking Color",
                col.get());
        if (newcolor != null) {
            RendererCache.pickingColor = new Color3f(newcolor);
        }
    }

    /**
     **********************************************************************
     *
     * @param format a {@link java.lang.String} object.
     */
    // Export Commands
    public void setImageFormat(String format) {
        if (format == null) {
            return;
        }
        format = format.toUpperCase();
        if (!imageFormatHash.containsKey(format)) {
            return;
        }
        imageFormat = imageFormatHash.get(format);
    }

    /**
     * <p>setLabelFontColor</p>
     */
    public void setLabelFontColor() {
        Color color = getGraphics2D().getColor();
        Color newColor = JColorChooser.showDialog(this, "Choose Font Color",
                color);
        if (newColor != null && newColor != color) {
            getGraphics2D().setPaint(newColor);
            if (RendererCache.labelAtoms || RendererCache.labelResidues) {
                repaint();
            }
            getStatusBar().setText(
                    "  Label Font Color Changed to (" + newColor.getRed() + ","
                    + newColor.getGreen() + "," + newColor.getBlue()
                    + ")");
        }
    }

    /**
     * <p>setLabelFontSize</p>
     */
    public void setLabelFontSize() {
        Font currentFont = getGraphics2D().getFont();
        int currentSize = currentFont.getSize();
        String size = new String("" + currentSize);
        size = JOptionPane.showInputDialog("Set the Font Size (8 to 64)", size);
        try {
            int f = Integer.parseInt(size);
            if (f < 8 || f > 64 || f == currentSize) {
                return;
            }
            Font newFont = new Font(currentFont.getName(), currentFont.getStyle(), f);
            getGraphics2D().setFont(newFont);
            if (RendererCache.labelAtoms || RendererCache.labelResidues) {
                repaint();
            }
            getStatusBar().setText(
                    "  Label Font Size Changed to " + newFont.getSize());
        } catch (NumberFormatException e) {
            return;
        }
    }

    /**
     **********************************************************************
     */
    // Misc. get and set methods.
    public void setLabelsUpdated() {
        repaint();
    }

    /**
     **********************************************************************
     *
     * @param level a {@link java.lang.String} object.
     */
    // Picking Commands
    public void setPickingLevel(String level) {
        if (level == null) {
            return;
        }
        level = level.toUpperCase();
        if (GraphicsPicking.pickLevelHash.containsKey(level)) {
            GraphicsPicking.PickLevel pickLevel = GraphicsPicking.pickLevelHash.get(level);
            switch (pickLevel) {
                case PICKATOM:
                case PICKBOND:
                case PICKANGLE:
                case PICKDIHEDRAL:
                case PICKRESIDUE:
                case PICKPOLYMER:
                case PICKSYSTEM:
                    rendererPicking.setPickLevel(level);
                    break;
                case MEASUREDISTANCE:
                case MEASUREANGLE:
                case MEASUREDIHEDRAL:
                    rendererPicking.setPickLevel(level);
                    rendererPicking.setPicking(true);
                    mainPanel.getMainMenu().setPickBehavior(true);
                    rendererPicking.clear();
                    rendererPicking.resetCount();
                    break;
                default:
                    logger.warning("Unexpected PickingLevel");
            }
        }
    }

    /**
     * <p>setPosition</p>
     */
    public void setPosition() {
        setPosition(mainPanel.getHierarchy().getActive());
    }

    /**
     * <p>setPosition</p>
     *
     * @param node a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setPosition(MSNode node) {
        updateScene(node, true, false, null, true, null);
    }

    /**
     * <p>setSelectionColor</p>
     */
    public void setSelectionColor() {
        Color3f col = new Color3f();
        RendererCache.selectionColor.get(col);
        Color newcolor = JColorChooser.showDialog(this,
                "Choose Selection Color", col.get());
        if (newcolor != null) {
            RendererCache.selectionColor = new Color3f(newcolor);
        }
        if (RendererCache.highlightSelections) {
            this.updateScene(mainPanel.getDataRoot(), false, false, null, true,
                    RendererCache.ColorModel.SELECT);
        }
    }

    /**
     * <p>setUserColor</p>
     */
    public void setUserColor() {
        Color3f col = new Color3f();
        RendererCache.selectionColor.get(col);
        Color newcolor = JColorChooser.showDialog(this, "Choose User Color",
                col.get());
        if (newcolor != null) {
            RendererCache.userColor = new Color3f(newcolor);
        }
    }

    /**
     **********************************************************************
     *
     * @param model a {@link java.lang.String} object.
     */
    // Display Commands
    public void setView(String model) {
        setViewModel(model);
    }

    /**
     * Operates on the active nodes.
     *
     * @param model
     *            String
     */
    public void setViewModel(String model) {
        if (!RendererCache.viewModelHash.containsKey(model.toUpperCase())) {
            return;
        }
        RendererCache.ViewModel viewModel = RendererCache.viewModelHash.get(model.toUpperCase());
        if (viewModel == RendererCache.ViewModel.RESTRICT) {
            renderer.arm(mainPanel.getDataRoot(), false, true, viewModel,
                    false, null);
            return;
        }
        ArrayList<MSNode> active = mainPanel.getHierarchy().getActiveNodes();
        if (active == null) {
            return;
        }
        renderer.arm(active, false, true, viewModel, false, null);
    }

    /**
     * Operates on the supplied node.
     *
     * @param model
     *            String
     * @param node a {@link ffx.potential.bonded.MSNode} object.
     */
    public void setViewModel(String model, MSNode node) {
        if (node == null) {
            return;
        }
        if (!RendererCache.viewModelHash.containsKey(model.toUpperCase())) {
            return;
        }
        RendererCache.ViewModel viewModel = RendererCache.viewModelHash.get(model.toUpperCase());
        renderer.arm(node, false, true, viewModel, false, null);
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
        return "3D Graphics";
    }

    /**
     * <p>updateScene</p>
     *
     * @param n a {@link java.util.ArrayList} object.
     * @param t a boolean.
     * @param v a boolean.
     * @param newViewModel a {@link ffx.potential.bonded.RendererCache.ViewModel} object.
     * @param c a boolean.
     * @param newColorModel a {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void updateScene(ArrayList<MSNode> n, boolean t, boolean v,
            RendererCache.ViewModel newViewModel, boolean c,
            RendererCache.ColorModel newColorModel) {
        if (n != null) {
            renderer.arm(n, t, v, newViewModel, c, newColorModel);
        }
    }

    /**
     * <p>updateScene</p>
     *
     * @param n a {@link ffx.potential.bonded.MSNode} object.
     * @param t a boolean.
     * @param v a boolean.
     * @param newViewModel a {@link ffx.potential.bonded.RendererCache.ViewModel} object.
     * @param c a boolean.
     * @param newColorModel a {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void updateScene(MSNode n, boolean t, boolean v,
            RendererCache.ViewModel newViewModel, boolean c,
            RendererCache.ColorModel newColorModel) {
        if (n != null) {
            renderer.arm(n, t, v, newViewModel, c, newColorModel);
        }
    }

    /**
     * <p>updateSceneWait</p>
     *
     * @param n a {@link java.util.ArrayList} object.
     * @param t a boolean.
     * @param v a boolean.
     * @param newViewModel a {@link ffx.potential.bonded.RendererCache.ViewModel} object.
     * @param c a boolean.
     * @param newColorModel a {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void updateSceneWait(ArrayList<MSNode> n, boolean t, boolean v,
            RendererCache.ViewModel newViewModel, boolean c,
            RendererCache.ColorModel newColorModel) {
        if (n != null) {
            renderer.arm(n, t, v, newViewModel, c, newColorModel);
        }
        while (isSceneRendering() || isCacheFull()) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    /**
     * <p>viewWait</p>
     *
     * @param viewModel a {@link ffx.potential.bonded.RendererCache.ViewModel} object.
     */
    public void viewWait(ViewModel viewModel) {
        if (viewModel == null) {
            logger.info("Null view.");
            return;
        }
        ArrayList<MSNode> nodes = mainPanel.getHierarchy().getActiveNodes();
        if (nodes == null) {
            logger.info("No active nodes.");
            return;
        }
        updateSceneWait(nodes, false, true, viewModel, false, null);
    }

    /**
     * <p>viewWait</p>
     *
     * @param viewMode a {@link java.lang.String} object.
     */
    public void viewWait(String viewMode) {
        if (viewMode == null) {
            logger.info("Null view.");
            return;
        }
        try {
            ViewModel viewModel = ViewModel.valueOf(viewMode.toUpperCase());
            viewWait(viewModel);
        } catch (Exception e) {
            logger.info("Unknown view command.");
        }
    }

    /**
     * <p>colorWait</p>
     *
     * @param colorMode a {@link java.lang.String} object.
     */
    public void colorWait(String colorMode) {
        if (colorMode == null) {
            logger.info("Null color.");
            return;
        }
        try {
            ColorModel colorModel = ColorModel.valueOf(colorMode.toUpperCase());
            colorWait(colorModel);
        } catch (Exception e) {
            logger.info("Unknown color command.");
        }

    }

    /**
     * <p>colorWait</p>
     *
     * @param colorModel a {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void colorWait(ColorModel colorModel) {
        if (colorModel == null) {
            logger.info("Null color.");
            return;
        }
        ArrayList<MSNode> nodes = mainPanel.getHierarchy().getActiveNodes();
        if (nodes == null) {
            logger.info("No active nodes.");
            return;
        }
        updateSceneWait(nodes, false, false, null, true, colorModel);
    }

    /**
     * <p>updateSceneWait</p>
     *
     * @param n a {@link ffx.potential.bonded.MSNode} object.
     * @param t a boolean.
     * @param v a boolean.
     * @param newViewModel a {@link ffx.potential.bonded.RendererCache.ViewModel} object.
     * @param c a boolean.
     * @param newColorModel a {@link ffx.potential.bonded.RendererCache.ColorModel} object.
     */
    public void updateSceneWait(MSNode n, boolean t, boolean v,
            RendererCache.ViewModel newViewModel, boolean c,
            RendererCache.ColorModel newColorModel) {
        if (n != null) {
            renderer.arm(n, t, v, newViewModel, c, newColorModel);
        }
        while (isSceneRendering() || isCacheFull()) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    String message = "Exception waiting for a Graphics operation.";
                    logger.log(Level.WARNING, message, e);
                }
            }
        }
    }

    /**
     * <p>zoomIn</p>
     */
    public void zoomIn() {
        baseTransformGroup.getTransform(baseTransform3D);
        double scale = baseTransform3D.getScale() + 0.01;
        baseTransform3D.setScale(scale);
        baseTransformGroup.setTransform(baseTransform3D);
    }

    /**
     * <p>zoomOut</p>
     */
    public void zoomOut() {
        baseTransformGroup.getTransform(baseTransform3D);
        double scale = baseTransform3D.getScale() - 0.01;
        if (scale > 0.0) {
            baseTransform3D.setScale(scale);
            baseTransformGroup.setTransform(baseTransform3D);
        }
    }
}
