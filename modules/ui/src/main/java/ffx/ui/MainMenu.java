/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2014.
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
package ffx.ui;

import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.*;

import static javax.swing.Action.*;

import org.apache.commons.lang3.SystemUtils;

import ffx.ui.properties.FFXLocale;

/**
 * The MainMenu class creates the Force Field X Menu Bar
 *
 * @author Michael J. Schnieders
 *
 */
public class MainMenu extends JMenuBar {

    private static final long serialVersionUID = 1L;
    private static final int keyMask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();
    // Controller References
    private final MainPanel mainPanel;
    private final GraphicsCanvas graphics;
    // Locale and ClassLoader
    private FFXLocale locale;
    private ClassLoader loader;
    // Toolbar
    private JToolBar toolBar;
    private Insets insets;
    private ImageIcon blankIcon;
    // Selection Menu
    private JCheckBoxMenuItem highlightCBMI;
    private JCheckBoxMenuItem labelResiduesMI;
    private JCheckBoxMenuItem labelAtomsMI;
    // Options Menu
    private JRadioButtonMenuItem activeRBMI;
    private JRadioButtonMenuItem mouseRBMI;
    private ButtonGroup dragModeButtonGroup;
    private ButtonGroup leftMouseButtonGroup;
    private JRadioButtonMenuItem rotateRBMI;
    // Picking Menu
    private JCheckBoxMenuItem pickingCBMI;
    private JRadioButtonMenuItem atomRBMI;
    private JRadioButtonMenuItem bondRBMI;
    private JRadioButtonMenuItem angleRBMI;
    private JRadioButtonMenuItem dihedralRBMI;
    private JRadioButtonMenuItem residueRBMI;
    private JRadioButtonMenuItem polymerRBMI;
    private JRadioButtonMenuItem moleculeRBMI;
    private JRadioButtonMenuItem systemRBMI;
    private JRadioButtonMenuItem measureDistanceRBMI;
    private JRadioButtonMenuItem measureAngleRBMI;
    private JRadioButtonMenuItem measureDihedralRBMI;
    private ButtonGroup levelBG;
    // Trajectory Menu
    private JCheckBoxMenuItem oscillateCBMI;
    // Simulation Menu
    private JMenuItem localMI;
    private JMenuItem remoteMI;
    private JMenuItem releaseMI;
    // Export Menu
    private ButtonGroup captureFormatButtonGroup;
    // Window Menu
    private JCheckBoxMenuItem systemsCBMI;
    private JCheckBoxMenuItem toolBarCBMI;
    private JCheckBoxMenuItem globalAxisCBMI;
    private final String icons = "ffx/ui/icons/";

    /*
     * Constructor
     * @param f Main application controller.
     */
    /**
     * <p>Constructor for MainMenu.</p>
     *
     * @param f a {@link ffx.ui.MainPanel} object.
     */
    public MainMenu(MainPanel f) {

        // Create the Tool Bar
        toolBar = new JToolBar("Toolbar", JToolBar.HORIZONTAL);
        toolBar.setBorderPainted(true);
        toolBar.setRollover(true);
        JButton temp = new JButton();
        insets = temp.getInsets();
        insets.set(2, 2, 2, 2);

        mainPanel = f;
        graphics = mainPanel.getGraphics3D();
        locale = mainPanel.getFFXLocale();
        loader = getClass().getClassLoader();
        blankIcon = new ImageIcon(loader.getResource(icons + "blank.gif"));

        /**
         * Main Menubar
         */
        JMenu fileMenu = addMenu("File", 'F');
        JMenu selectionMenu = addMenu("Selection", 'E');
        JMenu displayMenu = addMenu("Display", 'D');
        JMenu colorMenu = addMenu("Color", 'C');
        JMenu optionsMenu = addMenu("Options", 'O');
        JMenu pickingMenu = addMenu("Picking", 'P');
        JMenu trajectoryMenu = addMenu("Trajectory", 'T');
        JMenu exportMenu = addMenu("Export", 'X');
        JMenu windowMenu = addMenu("Window", 'W');
        JMenu helpMenu = addMenu("Help", 'H');

        /**
         * File Menu - Events Handled by the MainPanel Class.
         */
        addMenuItem(fileMenu, icons + "folder_page", "Open", 'O', KeyEvent.VK_O, mainPanel);
        addMenuItem(fileMenu, icons + "disk", "SaveAs", 'S', KeyEvent.VK_S, mainPanel);
        addMenuItem(fileMenu, icons + "cancel", "Close", 'C', -1, mainPanel);
        addMenuItem(fileMenu, "BLANK", "CloseAll", 'A', -1, mainPanel);
        fileMenu.addSeparator();
        addMenuItem(fileMenu, icons + "drive_web", "DownloadFromPDB", 'D', KeyEvent.VK_D, mainPanel);
        fileMenu.addSeparator();
        addMenuItem(fileMenu, "BLANK", "ChooseKeyFile", 'R', -1, mainPanel);
        addMenuItem(fileMenu, "BLANK", "ChooseLogFile", 'I', -1, mainPanel);
        addMenuItem(fileMenu, "BLANK", "LoadRestartData", 'R', -1, mainPanel);
        if (!SystemUtils.IS_OS_MAC_OSX) {
            fileMenu.addSeparator();
            addMenuItem(fileMenu, "BLANK", "Exit", 'E', KeyEvent.VK_Q, mainPanel);
        }
        toolBar.addSeparator();

        /**
         * Selection Menu - Events Handled by the MainPanel and GraphicsCanvas.
         */
        addMenuItem(selectionMenu, icons + "add", "SelectAll", 'A', KeyEvent.VK_A, mainPanel);
        addMenuItem(selectionMenu, "BLANK", "RestrictToSelections", 'R', -1, graphics);
        addMenuItem(selectionMenu, icons + "arrow_merge", "MergeSelections", 'M', -1, mainPanel);
        selectionMenu.addSeparator();
        highlightCBMI = addCBMenuItem(selectionMenu, icons + "asterisk_yellow", "HighlightSelections", 'H', KeyEvent.VK_H, mainPanel);
        addMenuItem(selectionMenu, "BLANK", "SetSelectionColor", 'S', -1, graphics);
        selectionMenu.addSeparator();
        labelAtomsMI = addCBMenuItem(selectionMenu, "BLANK", "LabelSelectedAtoms", 'O', -1, graphics);
        labelResiduesMI = addCBMenuItem(selectionMenu, "BLANK", "LabelSelectedResidues", 'R', -1, graphics);
        addMenuItem(selectionMenu, "BLANK", "SetLabelFontSize", 'Z', -1, graphics);
        addMenuItem(selectionMenu, "BLANK", "SetLabelFontColor", 'C', -1, graphics);
        highlightCBMI.setSelected(false);
        labelAtomsMI.setSelected(false);
        labelResiduesMI.setSelected(false);
        toolBar.addSeparator();

        /**
         * Display Menu - Events handled by the GraphicsCanvas.
         */
        addMenuItem(displayMenu, "BLANK", "Wireframe", 'W', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "Tube", 'T', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "Spacefill", 'S', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "BallAndStick", 'B', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "Invisible", 'I', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "RMIN", 'R', -1, graphics);
        displayMenu.addSeparator();
        addMenuItem(displayMenu, "BLANK", "ShowHydrogens", 'H', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "HideHydrogens", 'Y', -1, graphics);
        displayMenu.addSeparator();
        addMenuItem(displayMenu, "BLANK", "Fill", 'F', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "Points", 'P', -1, graphics);
        addMenuItem(displayMenu, "BLANK", "Lines", 'I', -1, graphics);
        displayMenu.addSeparator();
        addMenuItem(displayMenu, "BLANK", "Preferences", 'P', -1, graphics);

        /**
         * Color Menu - Events handled by the GraphicsCanvas.
         */
        addMenuItem(colorMenu, "BLANK", "Monochrome", 'M', -1, graphics);
        addMenuItem(colorMenu, "BLANK", "CPK", 'C', -1, graphics);
        addMenuItem(colorMenu, "BLANK", "Residue", 'R', -1, graphics);
        addMenuItem(colorMenu, "BLANK", "Structure", 'S', -1, graphics);
        addMenuItem(colorMenu, "BLANK", "Polymer", 'M', -1, graphics);
        addMenuItem(colorMenu, "BLANK", "PartialCharge", 'P', -1, graphics);
        addMenuItem(colorMenu, "BLANK", "UserColor", 'U', -1, graphics);
        colorMenu.addSeparator();
        addMenuItem(colorMenu, "BLANK", "ApplyUserColor", 'A', -1, graphics);
        addMenuItem(colorMenu, "BLANK", "SetUserColor", 'C', -1, graphics);

        /**
         * Options Menu - Events handled by the GraphicsCanvas.
         */
        dragModeButtonGroup = new ButtonGroup();
        activeRBMI = addBGMI(dragModeButtonGroup, optionsMenu, "BLANK", "ActiveSystem", 'A', KeyEvent.VK_A, graphics);
        mouseRBMI = addBGMI(dragModeButtonGroup, optionsMenu, "BLANK", "SystemBelowMouse", 'S', KeyEvent.VK_M, graphics);
        activeRBMI.setSelected(true);
        optionsMenu.addSeparator();
        JMenu leftMouseMenu = addSubMenu(optionsMenu, "LeftMouseButton", 'M');
        leftMouseButtonGroup = new ButtonGroup();
        rotateRBMI = addBGMI(leftMouseButtonGroup, leftMouseMenu, "BLANK", "Rotate", 'R', KeyEvent.VK_R, graphics);
        addBGMI(leftMouseButtonGroup, leftMouseMenu, "BLANK", "Translate", 'T', KeyEvent.VK_T, graphics);
        addBGMI(leftMouseButtonGroup, leftMouseMenu, "BLANK", "Zoom", 'Z', KeyEvent.VK_Z, graphics);
        rotateRBMI.setSelected(true);
        optionsMenu.addSeparator();
        addMenuItem(optionsMenu, "BLANK", "RotateAboutCenter", 'C', KeyEvent.VK_C, graphics);
        addMenuItem(optionsMenu, "BLANK", "RotateAboutPick", 'P', KeyEvent.VK_P, graphics);
        addMenuItem(optionsMenu, "BLANK", "ResetRotation", 'R', -1, graphics);
        addMenuItem(optionsMenu, "BLANK", "ResetTranslation", 'T', -1, graphics);
        addMenuItem(optionsMenu, icons + "arrow_refresh", "ResetRotationAndTranslation", 'E', -1, graphics);
        optionsMenu.addSeparator();
        addMenuItem(optionsMenu, icons + "magnifier_zoom_in", "ZoomIn", 'I', -1, graphics);
        addMenuItem(optionsMenu, icons + "magnifier_zoom_out", "ZoomOut", 'O', -1, graphics);
        addMenuItem(optionsMenu, "BLANK", "ResetGlobalZoom", 'Z', -1, graphics);
        addMenuItem(optionsMenu, "BLANK", "ResetGlobalRotation", 'N', -1, graphics);
        addMenuItem(optionsMenu, "BLANK", "ResetGlobalTranslation", 'O', -1, graphics);
        addMenuItem(optionsMenu, icons + "house", "ResetGlobalView", 'V', -1, graphics);
        optionsMenu.addSeparator();
        addMenuItem(optionsMenu, "BLANK", "SetBackgroundColor", 'B', -1, graphics);
        toolBar.addSeparator();

        /**
         * Picking Menu - Events handled by the GraphicsCanvas.
         */
        levelBG = new ButtonGroup();
        pickingCBMI = addCBMenuItem(pickingMenu, icons + "wand", "GraphicsPicking", 'G', KeyEvent.VK_0, graphics);
        pickingCBMI.setSelected(false);
        pickingMenu.addSeparator();
        atomRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickAtom", 'A', KeyEvent.VK_1, graphics);
        bondRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickBond", 'B', KeyEvent.VK_2, graphics);
        angleRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickAngle", 'N', KeyEvent.VK_3, graphics);
        dihedralRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickDihedral", 'D', KeyEvent.VK_4, graphics);
        residueRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickResidue", 'R', KeyEvent.VK_5, graphics);
        polymerRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickPolymer", 'P', KeyEvent.VK_6, graphics);
        moleculeRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickMolecule", 'M', KeyEvent.VK_7, graphics);
        systemRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickSystem", 'S', KeyEvent.VK_8, graphics);
        pickingMenu.addSeparator();
        measureDistanceRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "MeasureDistance", 'I', -1, graphics);
        measureAngleRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "MeasureAngle", 'L', -1, graphics);
        measureDihedralRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "MeasureDihedral", 'H', -1, graphics);
        atomRBMI.setSelected(true);
        pickingMenu.addSeparator();
        addMenuItem(pickingMenu, "BLANK", "SetGraphicsPickingColor", 'S', -1, graphics);
        toolBar.addSeparator();

        /**
         * Trajectory Menu - Events handled by the MainPanel.
         */
        oscillateCBMI = addCBMenuItem(trajectoryMenu, icons + "control_repeat_blue", "Oscillate", 'O', -1, mainPanel);
        oscillateCBMI.setSelected(false);
        addMenuItem(trajectoryMenu, "BLANK", "Frame", 'A', -1, mainPanel);
        addMenuItem(trajectoryMenu, "BLANK", "Speed", 'E', -1, mainPanel);
        addMenuItem(trajectoryMenu, "BLANK", "Skip", 'K', -1, mainPanel);
        trajectoryMenu.addSeparator();
        addMenuItem(trajectoryMenu, icons + "control_play_blue", "Play", 'P', -1, mainPanel);
        addMenuItem(trajectoryMenu, icons + "control_stop_blue", "Stop", 'S', -1, mainPanel);
        addMenuItem(trajectoryMenu, icons + "control_fastforward_blue", "StepForward", 'F', -1, mainPanel);
        addMenuItem(trajectoryMenu, icons + "control_rewind_blue", "StepBack", 'B', -1, mainPanel);
        addMenuItem(trajectoryMenu, icons + "control_start_blue", "Reset", 'R', -1, mainPanel);
        toolBar.addSeparator();

        /**
         * Export Menu - Events handled by the GraphicsCanvas.
         */
        addMenuItem(exportMenu, icons + "camera", "CaptureGraphics", 'C', KeyEvent.VK_G, graphics);
        exportMenu.addSeparator();
        captureFormatButtonGroup = new ButtonGroup();
        addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "PNG", 'P', -1, graphics).setSelected(true);
        addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "JPEG", 'J', -1, graphics);
        addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "BMP", 'B', -1, graphics);
        addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "WBMP", 'W', -1, graphics);
        addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "GIF", 'G', -1, graphics);
        toolBar.addSeparator();

        /**
         * Window Menu - Events handled by the GraphicsCanvas.
         */
        addMenuItem(windowMenu, icons + "application_home", "ResetPanes", 'R', -1, mainPanel);
        addMenuItem(windowMenu, icons + "application_osx_terminal", "ResetConsole", 'L', -1, mainPanel);
        windowMenu.addSeparator();
        addMenuItem(windowMenu, icons + "application_side_contract", "ExpandGraphicsWindow", 'E', -1, mainPanel);
        addMenuItem(windowMenu, icons + "application_side_expand", "ShrinkGraphicsWindow", 'L', -1, mainPanel);
        windowMenu.addSeparator();
        systemsCBMI = addCBMenuItem(windowMenu, "BLANK", "ShowTree", 'T', -1, mainPanel);
        toolBarCBMI = addCBMenuItem(windowMenu, "BLANK", "ShowToolBar", 'B', -1, mainPanel);
        globalAxisCBMI = addCBMenuItem(windowMenu, "BLANK", "ShowGlobalAxes", 'C', -1, mainPanel);
        globalAxisCBMI.setSelected(true);
        toolBar.addSeparator();

        /**
         * Help Menu - Events handled by the MainPanel.
         */
        Action a = addMenuItem(helpMenu, icons + "help", "HelpContents", 'H', KeyEvent.VK_HELP, mainPanel);
        /**
         * Fix the ACCELERATOR_KEY for the Help menu item; no modifiers will be
         * used.
         */
        a.putValue(ACCELERATOR_KEY, KeyStroke.getKeyStroke(KeyEvent.VK_HELP, 0));

        addMenuItem(helpMenu, "BLANK", "About", 'A', -1, mainPanel);

    }

    private JMenu addMenu(String name, char mnemonic) {
        JMenu menu = new JMenu(locale.getValue(name));
        add(menu);
        if (mnemonic != '.') {
            menu.setMnemonic(mnemonic);
        }
        return menu;
    }

    private Action addMenuItem(JMenu menu, String icon,
            String actionCommand, int mnemonic, int accelerator,
            final ActionListener actionListener) {
        Action a = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                actionListener.actionPerformed(e);
            }
        };
        configureAction(a, icon, actionCommand, mnemonic, accelerator);
        JMenuItem menuItem = new JMenuItem(a);
        menu.add(menuItem);
        return a;
    }

    private JRadioButtonMenuItem addBGMI(ButtonGroup buttonGroup, JMenu menu,
            String icon, String actionCommand, int mnemonic, int accelerator,
            final ActionListener actionListener) {

        final JRadioButtonMenuItem menuItem = new JRadioButtonMenuItem();

        Action a = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                /**
                 * If the ActionEvent is from a ToolBar button, pass it through
                 * the JRadioButtonMenuItem.
                 */
                if (e.getSource() != menuItem) {
                    menuItem.doClick();
                    return;
                }
                actionListener.actionPerformed(e);
            }
        };
        this.configureAction(a, icon, actionCommand, mnemonic, accelerator);
        menuItem.setAction(a);
        buttonGroup.add(menuItem);
        menu.add(menuItem);
        return menuItem;
    }

    private JCheckBoxMenuItem addCBMenuItem(JMenu menu, String icon,
            String actionCommand, int mnemonic, int accelerator,
            final ActionListener actionListener) {

        final JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();

        Action a = new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                /**
                 * If the ActionEvent is from a ToolBar button, pass it through
                 * the JCheckBoxMenuItem.
                 */
                if (e.getSource() != menuItem) {
                    menuItem.doClick();
                    return;
                }
                actionListener.actionPerformed(e);
            }
        };
        configureAction(a, icon, actionCommand, mnemonic, accelerator);
        menuItem.setAction(a);
        menu.add(menuItem);
        return menuItem;
    }

    private void configureAction(Action a, String icon, String actionCommand,
            int mnemonic, int accelerator) {
        String name = locale.getValue(actionCommand);
        ImageIcon imageIcon = getIcon(icon);
        KeyStroke keyStroke = KeyStroke.getKeyStroke(accelerator, keyMask);
        a.putValue(ACCELERATOR_KEY, keyStroke);
        a.putValue(ACTION_COMMAND_KEY, actionCommand);
        a.putValue(MNEMONIC_KEY, mnemonic);
        a.putValue(NAME, name);
        a.putValue(SHORT_DESCRIPTION, name);
        a.putValue(LONG_DESCRIPTION, name);
        a.putValue(SMALL_ICON, imageIcon);
        a.putValue(LARGE_ICON_KEY, imageIcon);
        if (!icon.equalsIgnoreCase("BLANK")) {
            toolBar.add(a);
        }
    }

    private ImageIcon getIcon(String icon) {
        ImageIcon imageIcon = null;
        if (!icon.equalsIgnoreCase("BLANK")) {
            imageIcon = new ImageIcon(loader.getResource(icon + ".png"));
        } else {
            imageIcon = blankIcon;
        }
        return imageIcon;
    }

    private JMenu addSubMenu(JMenu parent, String name, char mnemonic) {
        JMenu menu = new JMenu(locale.getValue(name));
        parent.add(menu);
        if (mnemonic != '.') {
            menu.setMnemonic(mnemonic);
        }
        menu.setIcon(blankIcon);
        return menu;
    }

    /**
     * <p>getHighlighting</p>
     *
     * @return a boolean.
     */
    public boolean getHighlighting() {
        return highlightCBMI.isSelected();
    }

    /**
     * <p>getMouseMode</p>
     *
     * @return a {@link ffx.ui.GraphicsCanvas.MouseMode} object.
     */
    public GraphicsCanvas.MouseMode getMouseMode() {
        if (activeRBMI.isSelected()) {
            return GraphicsCanvas.MouseMode.ACTIVESYSTEM;
        }
        return GraphicsCanvas.MouseMode.SYSTEMBELOWMOUSE;
    }

    /**
     * <p>getPicking</p>
     *
     * @return a boolean.
     */
    public boolean getPicking() {
        return pickingCBMI.isSelected();
    }

    /**
     * Get a reference the tool bar
     *
     * @return Force Field X ToolBar
     */
    public JToolBar getToolBar() {
        return toolBar;
    }

    /**
     * <p>isAxisShowing</p>
     *
     * @return a boolean.
     */
    public boolean isAxisShowing() {
        return globalAxisCBMI.isSelected();
    }

    /**
     * <p>isMenuShowing</p>
     *
     * @return a boolean.
     */
    public boolean isMenuShowing() {
        return toolBarCBMI.isSelected();
    }

    /**
     * <p>isPickingActive</p>
     *
     * @return a boolean.
     */
    public boolean isPickingActive() {
        return pickingCBMI.isSelected();
    }

    /**
     * <p>toggleSystemShowing</p>
     */
    public void toggleSystemShowing() {
        systemsCBMI.doClick();
    }

    /**
     * <p>isSystemShowing</p>
     *
     * @return a boolean.
     */
    public boolean isSystemShowing() {
        return systemsCBMI.isSelected();
    }

    /**
     * <p>toggleToolBarShowing</p>
     */
    public void toggleToolBarShowing() {
        toolBarCBMI.doClick();
    }

    /**
     * <p>setAtomLabels</p>
     *
     * @param b a boolean.
     */
    public void setAtomLabels(boolean b) {
        labelAtomsMI.setSelected(b);
    }

    /**
     * <p>setAxisShowing</p>
     *
     * @param b a boolean.
     */
    public void setAxisShowing(boolean b) {
        globalAxisCBMI.setSelected(b);
    }

    /**
     * Toggle connection status
     *
     * @param b a boolean.
     */
    public void setConnect(boolean b) {
        localMI.setEnabled(b);
        remoteMI.setEnabled(b);
        releaseMI.setEnabled(!b);
    }

    /**
     * <p>setHighlighting</p>
     *
     * @param h a boolean.
     */
    public void setHighlighting(boolean h) {
        highlightCBMI.setSelected(h);
    }

    /**
     * <p>setMenuShowing</p>
     *
     * @param b a boolean.
     */
    public void setMenuShowing(boolean b) {
        toolBarCBMI.setSelected(b);
    }

    /**
     * <p>setMouseMode</p>
     *
     * @param m a {@link ffx.ui.GraphicsCanvas.MouseMode} object.
     */
    public void setMouseMode(GraphicsCanvas.MouseMode m) {
        if (m == GraphicsCanvas.MouseMode.ACTIVESYSTEM) {
            activeRBMI.doClick();
        } else {
            mouseRBMI.doClick();
        }
    }

    /**
     * <p>setPickBehavior</p>
     *
     * @param pick a boolean.
     */
    public void setPickBehavior(boolean pick) {
        pickingCBMI.setSelected(pick);
    }

    /**
     * <p>setPickLevel</p>
     *
     * @param arg a {@link java.lang.String} object.
     */
    public void setPickLevel(String arg) {
        if (arg.equals("PickAtom")) {
            atomRBMI.doClick();
        } else if (arg.equals("PickBond")) {
            bondRBMI.doClick();
        } else if (arg.equals("PickAngle")) {
            angleRBMI.doClick();
        } else if (arg.equals("PickDihedral")) {
            dihedralRBMI.doClick();
        } else if (arg.equals("PickResidue")) {
            residueRBMI.doClick();
        } else if (arg.equals("PickPolymer")) {
            polymerRBMI.doClick();
        } else if (arg.equals("PickMolecule")) {
            moleculeRBMI.doClick();
        } else if (arg.equals("PickSystem")) {
            systemRBMI.doClick();
        } else if (arg.equals("MeasureDistance")) {
            measureDistanceRBMI.doClick();
        } else if (arg.equals("MeasureAngle")) {
            measureAngleRBMI.doClick();
        } else if (arg.equals("MeasureDihedral")) {
            measureDihedralRBMI.doClick();
        }
    }

    /**
     * <p>setResidueLabels</p>
     *
     * @param b a boolean.
     */
    public void setResidueLabels(boolean b) {
        labelResiduesMI.setSelected(b);
    }

    /**
     * <p>setSystemShowing</p>
     *
     * @param b a boolean.
     */
    public void setSystemShowing(boolean b) {
        systemsCBMI.setSelected(b);
    }

    /**
     * <p>systemClick</p>
     */
    public void systemClick() {
        systemsCBMI.doClick();
    }
}
