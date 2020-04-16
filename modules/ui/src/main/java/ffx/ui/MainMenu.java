// ******************************************************************************
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
// ******************************************************************************
package ffx.ui;

import static javax.swing.Action.ACCELERATOR_KEY;
import static javax.swing.Action.ACTION_COMMAND_KEY;
import static javax.swing.Action.LARGE_ICON_KEY;
import static javax.swing.Action.LONG_DESCRIPTION;
import static javax.swing.Action.MNEMONIC_KEY;
import static javax.swing.Action.NAME;
import static javax.swing.Action.SHORT_DESCRIPTION;
import static javax.swing.Action.SMALL_ICON;

import ffx.ui.properties.FFXLocale;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JToolBar;
import javax.swing.KeyStroke;
import org.apache.commons.lang3.SystemUtils;

/**
 * The MainMenu class creates the Force Field X Menu Bar
 *
 * @author Michael J. Schnieders
 */
public class MainMenu extends JMenuBar {

  /**
   * Note that the getMenuShortcutKeyMask() is deprecated in JDK 10 and replaced with
   * getMenuShortcutKeyMaskEx(). However, this later method is not present in JDK 8/9.
   */
  @SuppressWarnings("deprecation")
  private static final int keyMask = java.awt.Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();

  // Locale and ClassLoader
  private FFXLocale locale;
  private ClassLoader loader;
  // Toolbar
  private JToolBar toolBar;
  private ImageIcon blankIcon;

  // Selection Menu
  private JCheckBoxMenuItem highlightCBMI;
  private JCheckBoxMenuItem labelResiduesMI;
  private JCheckBoxMenuItem labelAtomsMI;
  // Options Menu
  private JRadioButtonMenuItem activeRBMI;
  private JRadioButtonMenuItem mouseRBMI;
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
  // Simulation Menu
  private JMenuItem localMI;
  private JMenuItem remoteMI;
  private JMenuItem releaseMI;
  // Window Menu
  private JCheckBoxMenuItem systemsCBMI;
  private JCheckBoxMenuItem toolBarCBMI;
  private JCheckBoxMenuItem globalAxisCBMI;

  /*
   * Constructor
   * @param f Main application controller.
   */

  /**
   * Constructor for MainMenu.
   *
   * @param f a {@link ffx.ui.MainPanel} object.
   */
  public MainMenu(MainPanel f) {

    // Create the Tool Bar
    toolBar = new JToolBar("Toolbar", JToolBar.HORIZONTAL);
    toolBar.setBorderPainted(true);
    toolBar.setRollover(true);
    JButton temp = new JButton();
    Insets insets = temp.getInsets();
    insets.set(2, 2, 2, 2);

    // Controller References
    GraphicsCanvas graphics = f.getGraphics3D();
    locale = f.getFFXLocale();
    loader = getClass().getClassLoader();
    String icons = "ffx/ui/icons/";
    blankIcon = new ImageIcon(loader.getResource(icons + "blank.gif"));

    String value = System.getProperty("structures", "false").trim();
    // Structure Menu
    boolean includeStructureMenu;
    try {
      includeStructureMenu = Boolean.parseBoolean(value);
    } catch (Exception e) {
      includeStructureMenu = false;
    }

    // Main Menubar
    JMenu fileMenu = addMenu("File", 'F');
    JMenu selectionMenu = addMenu("Selection", 'E');
    JMenu structureMenu = null;
    if (includeStructureMenu) {
      structureMenu = addMenu("Structure", 'S');
    }
    JMenu displayMenu = addMenu("Display", 'D');
    JMenu colorMenu = addMenu("Color", 'C');
    JMenu optionsMenu = addMenu("Options", 'O');
    JMenu pickingMenu = addMenu("Picking", 'P');
    JMenu trajectoryMenu = addMenu("Trajectory", 'T');
    JMenu exportMenu = addMenu("Export", 'X');
    JMenu windowMenu = addMenu("Window", 'W');
    JMenu helpMenu = addMenu("Help", 'H');

    // File Menu - Events Handled by the MainPanel Class.
    addMenuItem(fileMenu, icons + "folder_page", "Open", 'O', KeyEvent.VK_O, f);
    addMenuItem(fileMenu, icons + "disk", "SaveAs", 'S', KeyEvent.VK_S, f);
    addMenuItem(fileMenu, icons + "cancel", "Close", 'C', -1, f);
    addMenuItem(fileMenu, "BLANK", "CloseAll", 'A', -1, f);
    fileMenu.addSeparator();
    addMenuItem(fileMenu, icons + "drive_web", "DownloadFromPDB", 'D', KeyEvent.VK_D, f);
    fileMenu.addSeparator();
    addMenuItem(fileMenu, "BLANK", "ChooseKeyFile", 'R', -1, f);
    addMenuItem(fileMenu, "BLANK", "ChooseLogFile", 'I', -1, f);
    addMenuItem(fileMenu, "BLANK", "LoadRestartData", 'R', -1, f);
    if (!SystemUtils.IS_OS_MAC_OSX) {
      fileMenu.addSeparator();
      addMenuItem(fileMenu, "BLANK", "Exit", 'E', KeyEvent.VK_Q, f);
    }
    toolBar.addSeparator();

    // Selection Menu - Events Handled by the MainPanel and GraphicsCanvas.
    addMenuItem(selectionMenu, icons + "add", "SelectAll", 'A', KeyEvent.VK_A, f);
    addMenuItem(selectionMenu, "BLANK", "RestrictToSelections", 'R', -1, graphics);
    addMenuItem(selectionMenu, icons + "arrow_merge", "MergeSelections", 'M', -1, f);
    selectionMenu.addSeparator();
    highlightCBMI =
        addCBMenuItem(
            selectionMenu, icons + "asterisk_yellow", "HighlightSelections", 'H', KeyEvent.VK_H, f);
    addMenuItem(selectionMenu, "BLANK", "SetSelectionColor", 'S', -1, graphics);
    selectionMenu.addSeparator();
    labelAtomsMI = addCBMenuItem(selectionMenu, "BLANK", "LabelSelectedAtoms", 'O', -1, graphics);
    labelResiduesMI =
        addCBMenuItem(selectionMenu, "BLANK", "LabelSelectedResidues", 'R', -1, graphics);
    addMenuItem(selectionMenu, "BLANK", "SetLabelFontSize", 'Z', -1, graphics);
    addMenuItem(selectionMenu, "BLANK", "SetLabelFontColor", 'C', -1, graphics);
    highlightCBMI.setSelected(false);
    labelAtomsMI.setSelected(false);
    labelResiduesMI.setSelected(false);
    toolBar.addSeparator();

    // Structure Menu - Events Handled by the MainPanel.
    if (includeStructureMenu) {
      // Locate a jar file that has PDB Structures.
      String file = "ffx/xray/structures/1N7S.pdb";
      addMenuItem(structureMenu, "BLANK", file, '.', -1, f);
    }

    // Display Menu - Events handled by the GraphicsCanvas.
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

    // Color Menu - Events handled by the GraphicsCanvas.
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

    // Options Menu - Events handled by the GraphicsCanvas.
    ButtonGroup dragModeButtonGroup = new ButtonGroup();
    activeRBMI =
        addBGMI(
            dragModeButtonGroup,
            optionsMenu,
            "BLANK",
            "ActiveSystem",
            'A',
            KeyEvent.VK_A,
            graphics);
    mouseRBMI =
        addBGMI(
            dragModeButtonGroup,
            optionsMenu,
            "BLANK",
            "SystemBelowMouse",
            'S',
            KeyEvent.VK_M,
            graphics);
    activeRBMI.setSelected(true);
    optionsMenu.addSeparator();
    JMenu leftMouseMenu = addSubMenu(optionsMenu, "LeftMouseButton", 'M');
    ButtonGroup leftMouseButtonGroup = new ButtonGroup();
    JRadioButtonMenuItem rotateRBMI =
        addBGMI(
            leftMouseButtonGroup, leftMouseMenu, "BLANK", "Rotate", 'R', KeyEvent.VK_R, graphics);
    addBGMI(
        leftMouseButtonGroup, leftMouseMenu, "BLANK", "Translate", 'T', KeyEvent.VK_T, graphics);
    addBGMI(leftMouseButtonGroup, leftMouseMenu, "BLANK", "Zoom", 'Z', KeyEvent.VK_Z, graphics);
    rotateRBMI.setSelected(true);
    optionsMenu.addSeparator();
    addMenuItem(optionsMenu, "BLANK", "RotateAboutCenter", 'C', KeyEvent.VK_C, graphics);
    addMenuItem(optionsMenu, "BLANK", "RotateAboutPick", 'P', KeyEvent.VK_P, graphics);
    addMenuItem(optionsMenu, "BLANK", "ResetRotation", 'R', -1, graphics);
    addMenuItem(optionsMenu, "BLANK", "ResetTranslation", 'T', -1, graphics);
    addMenuItem(
        optionsMenu, icons + "arrow_refresh", "ResetRotationAndTranslation", 'E', -1, graphics);
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

    // Picking Menu - Events handled by the GraphicsCanvas.
    ButtonGroup levelBG = new ButtonGroup();
    pickingCBMI =
        addCBMenuItem(pickingMenu, icons + "wand", "GraphicsPicking", 'G', KeyEvent.VK_0, graphics);
    pickingCBMI.setSelected(false);
    pickingMenu.addSeparator();
    atomRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickAtom", 'A', KeyEvent.VK_1, graphics);
    bondRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickBond", 'B', KeyEvent.VK_2, graphics);
    angleRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickAngle", 'N', KeyEvent.VK_3, graphics);
    dihedralRBMI =
        addBGMI(levelBG, pickingMenu, "BLANK", "PickDihedral", 'D', KeyEvent.VK_4, graphics);
    residueRBMI =
        addBGMI(levelBG, pickingMenu, "BLANK", "PickResidue", 'R', KeyEvent.VK_5, graphics);
    polymerRBMI =
        addBGMI(levelBG, pickingMenu, "BLANK", "PickPolymer", 'P', KeyEvent.VK_6, graphics);
    moleculeRBMI =
        addBGMI(levelBG, pickingMenu, "BLANK", "PickMolecule", 'M', KeyEvent.VK_7, graphics);
    systemRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "PickSystem", 'S', KeyEvent.VK_8, graphics);
    pickingMenu.addSeparator();
    measureDistanceRBMI =
        addBGMI(levelBG, pickingMenu, "BLANK", "MeasureDistance", 'I', -1, graphics);
    measureAngleRBMI = addBGMI(levelBG, pickingMenu, "BLANK", "MeasureAngle", 'L', -1, graphics);
    measureDihedralRBMI =
        addBGMI(levelBG, pickingMenu, "BLANK", "MeasureDihedral", 'H', -1, graphics);
    atomRBMI.setSelected(true);
    pickingMenu.addSeparator();
    addMenuItem(pickingMenu, "BLANK", "SetGraphicsPickingColor", 'S', -1, graphics);
    toolBar.addSeparator();

    // Trajectory Menu - Events handled by the MainPanel.
    JCheckBoxMenuItem oscillateCBMI =
        addCBMenuItem(trajectoryMenu, icons + "control_repeat_blue", "Oscillate", 'O', -1, f);
    oscillateCBMI.setSelected(false);
    addMenuItem(trajectoryMenu, "BLANK", "Frame", 'A', -1, f);
    addMenuItem(trajectoryMenu, "BLANK", "Speed", 'E', -1, f);
    addMenuItem(trajectoryMenu, "BLANK", "Skip", 'K', -1, f);
    trajectoryMenu.addSeparator();
    addMenuItem(trajectoryMenu, icons + "control_play_blue", "Play", 'P', -1, f);
    addMenuItem(trajectoryMenu, icons + "control_stop_blue", "Stop", 'S', -1, f);
    addMenuItem(trajectoryMenu, icons + "control_fastforward_blue", "StepForward", 'F', -1, f);
    addMenuItem(trajectoryMenu, icons + "control_rewind_blue", "StepBack", 'B', -1, f);
    addMenuItem(trajectoryMenu, icons + "control_start_blue", "Reset", 'R', -1, f);
    toolBar.addSeparator();

    // Export Menu - Events handled by the GraphicsCanvas.
    addMenuItem(exportMenu, icons + "camera", "CaptureGraphics", 'C', KeyEvent.VK_G, graphics);
    exportMenu.addSeparator();
    ButtonGroup captureFormatButtonGroup = new ButtonGroup();
    addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "PNG", 'P', -1, graphics)
        .setSelected(true);
    addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "JPEG", 'J', -1, graphics);
    addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "BMP", 'B', -1, graphics);
    addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "WBMP", 'W', -1, graphics);
    addBGMI(captureFormatButtonGroup, exportMenu, "BLANK", "GIF", 'G', -1, graphics);
    toolBar.addSeparator();

    // Window Menu - Events handled by the GraphicsCanvas.
    addMenuItem(windowMenu, icons + "application_home", "ResetPanes", 'R', -1, f);
    addMenuItem(windowMenu, icons + "application_osx_terminal", "ResetConsole", 'L', -1, f);
    windowMenu.addSeparator();
    addMenuItem(
        windowMenu, icons + "application_side_contract", "ExpandGraphicsWindow", 'E', -1, f);
    addMenuItem(windowMenu, icons + "application_side_expand", "ShrinkGraphicsWindow", 'L', -1, f);
    windowMenu.addSeparator();
    systemsCBMI = addCBMenuItem(windowMenu, "BLANK", "ShowTree", 'T', -1, f);
    toolBarCBMI = addCBMenuItem(windowMenu, "BLANK", "ShowToolBar", 'B', -1, f);
    globalAxisCBMI = addCBMenuItem(windowMenu, "BLANK", "ShowGlobalAxes", 'C', -1, f);
    globalAxisCBMI.setSelected(true);
    toolBar.addSeparator();

    addMenuItem(helpMenu, "BLANK", "About", 'A', -1, f);
  }

  /**
   * getHighlighting
   *
   * @return a boolean.
   */
  public boolean getHighlighting() {
    return highlightCBMI.isSelected();
  }

  /**
   * setHighlighting
   *
   * @param h a boolean.
   */
  void setHighlighting(boolean h) {
    highlightCBMI.setSelected(h);
  }

  /**
   * getMouseMode
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
   * setMouseMode
   *
   * @param m a {@link ffx.ui.GraphicsCanvas.MouseMode} object.
   */
  void setMouseMode(GraphicsCanvas.MouseMode m) {
    if (m == GraphicsCanvas.MouseMode.ACTIVESYSTEM) {
      activeRBMI.doClick();
    } else {
      mouseRBMI.doClick();
    }
  }

  /**
   * getPicking
   *
   * @return a boolean.
   */
  public boolean getPicking() {
    return pickingCBMI.isSelected();
  }

  /**
   * isPickingActive
   *
   * @return a boolean.
   */
  public boolean isPickingActive() {
    return pickingCBMI.isSelected();
  }

  /** toggleSystemShowing */
  public void toggleSystemShowing() {
    systemsCBMI.doClick();
  }

  private JMenu addMenu(String name, char mnemonic) {
    JMenu menu = new JMenu(locale.getValue(name));
    add(menu);
    if (mnemonic != '.') {
      menu.setMnemonic(mnemonic);
    }
    return menu;
  }

  private Action addMenuItem(
      JMenu menu,
      String icon,
      String actionCommand,
      int mnemonic,
      int accelerator,
      final ActionListener actionListener) {
    Action a =
        new AbstractAction() {
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

  private JRadioButtonMenuItem addBGMI(
      ButtonGroup buttonGroup,
      JMenu menu,
      String icon,
      String actionCommand,
      int mnemonic,
      int accelerator,
      final ActionListener actionListener) {

    final JRadioButtonMenuItem menuItem = new JRadioButtonMenuItem();

    Action a =
        new AbstractAction() {
          @Override
          public void actionPerformed(ActionEvent e) {

            // If the ActionEvent is from a ToolBar button, pass it through the
            // JRadioButtonMenuItem.
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

  private JCheckBoxMenuItem addCBMenuItem(
      JMenu menu,
      String icon,
      String actionCommand,
      int mnemonic,
      int accelerator,
      final ActionListener actionListener) {

    final JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem();

    Action a =
        new AbstractAction() {
          @Override
          public void actionPerformed(ActionEvent e) {
            // If the ActionEvent is from a ToolBar button, pass it through the JCheckBoxMenuItem.
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

  private void configureAction(
      Action a, String icon, String actionCommand, int mnemonic, int accelerator) {
    String name;
    try {
      name = locale.getValue(actionCommand);
    } catch (Exception e) {
      name = actionCommand;
    }
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
   * Get a reference the tool bar
   *
   * @return Force Field X ToolBar
   */
  JToolBar getToolBar() {
    return toolBar;
  }

  /**
   * isAxisShowing
   *
   * @return a boolean.
   */
  boolean isAxisShowing() {
    return globalAxisCBMI.isSelected();
  }

  /**
   * setAxisShowing
   *
   * @param b a boolean.
   */
  public void setAxisShowing(boolean b) {
    globalAxisCBMI.setSelected(b);
  }

  /**
   * isMenuShowing
   *
   * @return a boolean.
   */
  boolean isMenuShowing() {
    return toolBarCBMI.isSelected();
  }

  /**
   * setMenuShowing
   *
   * @param b a boolean.
   */
  void setMenuShowing(boolean b) {
    toolBarCBMI.setSelected(b);
  }

  /**
   * isSystemShowing
   *
   * @return a boolean.
   */
  boolean isSystemShowing() {
    return systemsCBMI.isSelected();
  }

  /**
   * setSystemShowing
   *
   * @param b a boolean.
   */
  void setSystemShowing(boolean b) {
    systemsCBMI.setSelected(b);
  }

  /** toggleToolBarShowing */
  void toggleToolBarShowing() {
    toolBarCBMI.doClick();
  }

  /**
   * setAtomLabels
   *
   * @param b a boolean.
   */
  void setAtomLabels(boolean b) {
    labelAtomsMI.setSelected(b);
  }

  /**
   * Toggle connection status
   *
   * @param b a boolean.
   */
  void setConnect(boolean b) {
    localMI.setEnabled(b);
    remoteMI.setEnabled(b);
    releaseMI.setEnabled(!b);
  }

  /**
   * setPickBehavior
   *
   * @param pick a boolean.
   */
  void setPickBehavior(boolean pick) {
    pickingCBMI.setSelected(pick);
  }

  /**
   * setPickLevel
   *
   * @param arg a {@link java.lang.String} object.
   */
  void setPickLevel(String arg) {
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
   * setResidueLabels
   *
   * @param b a boolean.
   */
  void setResidueLabels(boolean b) {
    labelResiduesMI.setSelected(b);
  }

  /** systemClick */
  void systemClick() {
    systemsCBMI.doClick();
  }
}
