// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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

import static ffx.utilities.FileUtils.copyInputStreamToTmpFile;
import static ffx.utilities.StringUtils.pdbForID;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MSRoot;
import ffx.potential.bonded.RendererCache;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ARCFileFilter;
import ffx.potential.parsers.FFXFileFilter;
import ffx.potential.parsers.ForceFieldFileFilter;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.INTFileFilter;
import ffx.potential.parsers.INTFilter;
import ffx.potential.parsers.InducedFileFilter;
import ffx.potential.parsers.InducedFilter;
import ffx.potential.parsers.KeyFileFilter;
import ffx.potential.parsers.KeyFilter;
import ffx.potential.parsers.MergeFilter;
import ffx.potential.parsers.PDBFileFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFileFilter;
import ffx.potential.parsers.XYZFilter;
import ffx.ui.properties.FFXLocale;
import ffx.utilities.Keyword;
import ffx.utilities.Resources;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsEnvironment;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.Preferences;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileSystemView;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.SystemUtils;
import org.jogamp.java3d.GraphicsConfigTemplate3D;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.vecmath.Vector3d;

/**
 * The MainPanel class is the main container for Force Field X, handles file input/output and is used
 * to pass references among the various sub-Panels.
 *
 * @author Michael J. Schnieders
 */
public final class MainPanel extends JPanel implements ActionListener, ChangeListener {

  @Serial
  private static final long serialVersionUID = 1L;

  /**
   * Constant <code>version="1.0.0"</code>
   */
  public static final String version = "1.0.0";
  /**
   * Constant <code>date="January 2025"</code>
   */
  public static final String date = "January 2025";
  /**
   * Constant
   */
  public static final String border =
      " _________________________________________________________________________\n";
  /**
   * Constant
   */
  public static final String title = "        FORCE FIELD X - Polyglot Software for Molecular Biophysics \n";

  public static final String aboutString;
  /**
   * Constant <code>KEYWORDS=1</code>
   */
  static final int KEYWORDS = 1;
  /**
   * Constant <code>MODELING=2</code>
   */
  static final int MODELING = 2;
  /**
   * Constant <code>forceFieldFileFilter</code>
   */
  static final ForceFieldFileFilter forceFieldFileFilter = new ForceFieldFileFilter();
  /**
   * Constant <code>keyFileFilter</code>
   */
  static final KeyFileFilter keyFileFilter = new KeyFileFilter();

  private static final Logger logger = Logger.getLogger(MainPanel.class.getName());
  /**
   * Constant <code>GRAPHICS=0</code>
   */
  private static final int GRAPHICS = 0;
  /**
   * Constant <code>xyzFileFilter</code>
   */
  private static final XYZFileFilter xyzFileFilter = new XYZFileFilter();
  /**
   * Constant <code>arcFileFilter</code>
   */
  private static final ARCFileFilter arcFileFilter = new ARCFileFilter();
  /**
   * Constant <code>intFileFilter</code>
   */
  private static final INTFileFilter intFileFilter = new INTFileFilter();
  /**
   * Constant <code>indFileFilter</code>
   */
  private static final InducedFileFilter indFileFilter = new InducedFileFilter();
  /**
   * Constant <code>pdbFileFilter</code>
   */
  private static final PDBFileFilter pdbFileFilter = new PDBFileFilter();
  /**
   * Constant <code>ffxFileFilter</code>
   */
  private static final FFXFileFilter ffxFileFilter = new FFXFileFilter();

  private static final Preferences preferences = Preferences.userNodeForPackage(MainPanel.class);
  /**
   * Constant <code>classpath=""</code>
   */
  static String classpath;
  /**
   * Constant <code>ffxDir</code>
   */
  static File ffxDir;
  /**
   * Present working directory.
   */
  private static File pwd;
  /**
   * JFileChooser for choosing a file.
   */
  private static JFileChooser fileChooser = null;

  static {
    var basedir = System.getProperty("basedir");
    var mvnProps = new File(basedir + "/bin/build.properties");
    var commitVersion = version + "-unknown";
    var commitDate = date;
    var commitSCM = "";
    if (mvnProps.exists()) {
      try (BufferedReader br = new BufferedReader(new FileReader(mvnProps))) {
        var ffxVersion = "1.0.0";
        var ffxVersionProp = "ffx.version=";
        var gitCommitsCount = "";
        var gitCommitsCountProp = "git.total.commit.count=";
        var timestampProp = "timestamp=";
        var gitRevisionProp = "git.commit.id.full=";
        var line = br.readLine();
        while (line != null) {
          line = line.trim();
          if (line.startsWith(ffxVersionProp)) {
            ffxVersion = line.replaceFirst(ffxVersionProp, "");
          } else if (line.startsWith(gitCommitsCountProp)) {
            gitCommitsCount = line.replaceFirst(gitCommitsCountProp, "");
          } else if (line.startsWith(timestampProp)) {
            var timeStr = line.replaceFirst(timestampProp, "");
            // Expected to be MM-dd-yyyy
            var timeToks = timeStr.split("-");
            try {
              var year = timeToks[2];
              int mon = Integer.parseInt(timeToks[0]);
              Month month = Month.of(mon);
              var mstr = month.toString();
              commitDate = format("%c%s %s", mstr.charAt(0), mstr.substring(1).toLowerCase(), year);
            } catch (Exception ex) {
              commitDate = date;
            }
          } else if (line.startsWith(gitRevisionProp) && !line.contains("UNKNOWN_REVISION")) {
            var scm = line.replaceFirst(gitRevisionProp, "");
            commitSCM = format("        %s %s \n", "Git Revision ", scm);
          }
          line = br.readLine();
        }
        var sb = new StringBuilder(ffxVersion).append("-");
        if (!gitCommitsCount.isEmpty()) {
          sb.append(gitCommitsCount);
        } else {
          sb.append("unknown");
        }
        commitVersion = sb.toString();
      } catch (Exception ex) {
        //
      }
    }
    aboutString = "        Version "
        + commitVersion
        + "  "
        + commitDate
        + " \n"
        + commitSCM // Will contain its own spacing/newline, or be empty.
        + " \n"
        + """                
                Please cite the following reference when using Force Field X:
        
                RA Gogal, AJ Nessler, AC Thiel, HV Bernabe, RA Corrigan Grove,
                LM Cousineau, JM Litman, JM Miller, G Qi, MJ Speranza,
                MR Tollefson, TD Fenn, JJ Michaelson, O Okada, JP Piquemal,
                JW Ponder, J Shen, RJH Smith, W Yang, P Ren and MJ Schnieders,
                2024, Journal of Chemical Physics, 161 (1).
        
                Copyright (c)  Michael J. Schnieders  2001-2025
                All Rights Reserved
        
                Force Field X is distributed under the GPL v. 3 license
                with linking exception and is hosted by the Schnieders Lab
                at The University of Iowa.
        
                User Manual:   https://ffx.biochem.uiowa.edu/manual.html
                Publications:  https://ffx.biochem.uiowa.edu/publications.html
                License:       https://ffx.biochem.uiowa.edu/licenses.html
        """;

    try {
      String ffxString = System.getProperty("ffx.dir", ".");
      ffxDir = new File(ffxString);
      classpath = System.getProperty("java.class.path");
      pwd = MainPanel.getPWD();
    } catch (Exception e) {
      logger.severe(" FFX directory could not be found.\n" + e);
    }
  }

  /**
   * Main FFX JFrame.
   */
  private final JFrame frame;
  /**
   * Number of File Opener Threads.
   */
  private int fileOpenerThreads = -1;
  /**
   * Root of the structural hierarchy.
   */
  private MSRoot dataRoot;
  /**
   * The structural hierarchy.
   */
  private Hierarchy hierarchy;
  /**
   * The FFX Main Menu.
   */
  private MainMenu mainMenu;
  /**
   * The Graphics Panel.
   */
  private GraphicsPanel graphicsPanel;
  /**
   * The Modeling Panel.
   */
  private ModelingPanel modelingPanel;
  /**
   * The Keyword Panel.
   */
  private KeywordPanel keywordPanel;
  /**
   * A reference to the Modeling Shell.
   */
  private ModelingShell modelingShell = null;
  /**
   * The Java3D Graphics Canvas.
   */
  private GraphicsCanvas graphicsCanvas;
  /**
   * The SplitPane holds the Hierarchy and JTabbedPane.
   */
  private JSplitPane splitPane;
  /**
   * The value fo the Split Pane Divider.
   */
  private int splitPaneDivider;
  /**
   * Status Label.
   */
  private JLabel statusLabel;
  /**
   * Filter to open a force field file.
   */
  private ForceFieldFilter forceFieldFilter;
  /**
   * The FFX Locale.
   */
  private FFXLocale locale = null;
  /**
   * The FFX About Dialog.
   */
  private JDialog aboutDialog = null;
  /**
   * The FFX About Text Area.
   */
  private JTextArea aboutTextArea = null;
  /**
   * Thread to open systems.
   */
  private Thread openThread = null;
  /**
   * The active system filter.
   */
  private SystemFilter activeFilter = null;
  /**
   * Flag to indicate oscillation.
   */
  private boolean oscillate = false;
  /**
   * Reference to a Simulation Loader.
   */
  private SimulationLoader simulation;
  /**
   * IP of the simulation.
   */
  private String ip = "";
  /**
   * Simulation port.
   */
  private int port = 2000;
  /**
   * InetAddress of the simulation.
   */
  private InetAddress address = null;
  /**
   * InetSocketAddress of the simulation.
   */
  private InetSocketAddress socketAddress = new InetSocketAddress(port);
  /**
   * Initialize all the sub-Panels and put them together
   */
  private boolean init = false;
  /**
   * Exit status to describe how FFX is terminating.
   */
  private ExitStatus exitType = ExitStatus.NORMAL;

  /**
   * MainPanel Constructor
   *
   * @param f Application Frame
   */
  public MainPanel(JFrame f) {
    frame = f;
  }

  /**
   * Constructor for MainPanel.
   */
  public MainPanel() {
    frame = null;
  }

  /**
   * getPWD
   *
   * @return a {@link java.io.File} object.
   */
  static File getPWD() {
    if (pwd == null) {
      pwd =
          new File(
              System.getProperty(
                  "user.dir",
                  FileSystemView.getFileSystemView().getDefaultDirectory().getAbsolutePath()));
    }
    return pwd;
  }

  /**
   * JFileChooser
   *
   * @return a {@link javax.swing.JFileChooser} object.
   */
  static JFileChooser resetFileChooser() {
    if (fileChooser == null) {
      fileChooser = new JFileChooser();
    }
    fileChooser.resetChoosableFileFilters();
    fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
    fileChooser.setFileHidingEnabled(false);
    fileChooser.setAcceptAllFileFilterUsed(true);
    fileChooser.setCurrentDirectory(getPWD());
    fileChooser.setMultiSelectionEnabled(false);
    fileChooser.setSelectedFile(null);
    return fileChooser;
  }

  /**
   * about
   */
  public void about() {
    if (aboutDialog == null) {
      aboutDialog = new JDialog(frame, "About... ", true);
      URL ffxURL = getClass().getClassLoader().getResource("ffx/ui/icons/splash.png");
      ImageIcon logoIcon = new ImageIcon(ffxURL);
      JLabel logoLabel = new JLabel(logoIcon);
      logoLabel.setSize(600, 600);
      logoLabel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
      Container contentpane = aboutDialog.getContentPane();
      contentpane.setLayout(new BorderLayout());
      initAbout();
      contentpane.add(aboutTextArea, BorderLayout.SOUTH);
      contentpane.add(logoLabel, BorderLayout.CENTER);
      aboutDialog.pack();
      Dimension dim = getToolkit().getScreenSize();
      Dimension ddim = aboutDialog.getSize();
      aboutDialog.setLocation((dim.width - ddim.width) / 2, (dim.height - ddim.height) / 2);
      aboutDialog.setResizable(false);
    }
    aboutDialog.setVisible(true);
  }

  /**
   * {@inheritDoc}
   *
   * <p>Handle most File, Selection, Trajectory, Simulation, Window and Help Menu Commands This
   * should probably be partitioned between a few different handlers
   */
  @Override
  public void actionPerformed(ActionEvent evt) {
    String arg = evt.getActionCommand();
    if (logger.isLoggable(Level.FINEST)) {
      logger.finest(" Action: " + arg);
    }
    // File Commands
    switch (arg) {
      case "Open" -> open();
      case "DownloadFromPDB" -> openFromPDB();
      case "SaveAs" -> saveAsXYZ(null);
      case "Close" -> close();
      case "CloseAll" -> closeAll();
      case "ChooseKeyFile" -> chooseKey();
      case "ChooseLogFile" -> chooseLog();
      case "LoadInducedData" -> openInduced();

      // Selection Commands
      case "SelectAll" -> selectAll();
      case "MergeSelections" -> merge();
      case "HighlightSelections" -> highlightSelections(evt);

      // Trajectory
      case "Play" -> play();
      case "Stop" -> stop();
      case "StepForward" -> stepForward();
      case "StepBack" -> stepBack();
      case "Reset" -> reset();
      case "Oscillate" -> oscillate(evt);
      case "Frame" -> frame();
      case "Speed" -> speed();
      case "Skip" -> skip();

      // Simulation
      case "ConnectToLocalJob" -> connectToTINKER(null, null);
      case "ConnectToRemoteJob" -> connect();
      case "ReleaseJob" -> release();
      case "SetPort" -> setPort();
      case "SetRemoteJobAddress" -> setRemoteJobAddress();

      // Window
      case "ShowToolBar" -> showToolBar(evt);
      case "ShowTree" -> showTree(evt);
      case "ShowGlobalAxes" -> showGlobalAxes(evt);
      case "ResetPanes" -> resetPanes();
      case "ResetConsole" -> resetShell();
      case "OceanLookAndFeel" -> oceanLookAndFeel();
      case "WindowsLookAndFeel", "MacOSXLookAndFeel", "MotifLookAndFeel" -> platformLookAndFeel();
      case "ShrinkGraphicsWindow" -> resizePanes(20);
      case "ExpandGraphicsWindow" -> resizePanes(-20);
      case "About" -> about();

      // Others
      case "GarbageCollect" -> Runtime.getRuntime().gc();
      case "Exit" -> exit();
      default -> {
        try {
          ClassLoader cl = MainPanel.class.getClassLoader();
          URL url = cl.getResource(arg);
          logger.info(url.toString());
          File structureFile = new File(url.getFile());
          logger.info(structureFile.toString());
          String tempFile =
              copyInputStreamToTmpFile(url.openStream(), "ffx", structureFile.getName(), "pdb");
          open(tempFile);
        } catch (Exception e) {
          System.err.println("MainPanel - Menu command not found: " + arg);
        }
      }
    }
  }

  /**
   * Detach the active FSystem's BranchGroup from the Scene and clear that FSystem's data
   *
   * @return a {@link java.lang.Thread} object.
   */
  public Thread close() {
    FFXSystem m = hierarchy.getActive();
    return close(m);
  }

  /**
   * close
   *
   * @param closedModel a {@link ffx.ui.FFXSystem} object.
   * @return a {@link java.lang.Thread} object.
   */
  public Thread close(FFXSystem closedModel) {
    if (closedModel.getParent() != null) {
      Trajectory traj = closedModel.getTrajectory();
      if (traj != null) {
        traj.stop();
      }
      if (simulation != null && simulation.getFSystem() == closedModel) {
        release();
      }
      hierarchy.removeTreeNode(closedModel);
      closedModel.setView(RendererCache.ViewModel.DESTROY, null);
      Thread thread = new Thread(new UIFileCloser(closedModel));
      thread.start();
      return thread;
    }
    return null;
  }

  /**
   * exit with current exit code (default: 0 (ExitStatus.NORMAL))
   */
  public void exit() {
    exit(exitType);
  }

  /**
   * frame
   */
  public void frame() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    String frameNumber = "" + trajectory.getFrame();
    frameNumber = JOptionPane.showInputDialog("Enter the Frame Number", frameNumber);
    try {
      int f = Integer.parseInt(frameNumber);
      trajectory.setFrame(f);
    } catch (NumberFormatException e) {
      //
    }
  }

  /**
   * Return the active SystemFilter.
   *
   * @return the active SystemFilter.
   */
  public SystemFilter getFilter() {
    return activeFilter;
  }

  /**
   * Getter for the field <code>frame</code>.
   *
   * @return a {@link java.awt.Frame} object.
   */
  public Frame getFrame() {
    return frame;
  }

  /**
   * Getter for the field <code>hierarchy</code>.
   *
   * @return a {@link ffx.ui.Hierarchy} object.
   */
  public Hierarchy getHierarchy() {
    return hierarchy;
  }

  /**
   * Getter for the field <code>mainMenu</code>.
   *
   * @return a {@link ffx.ui.MainMenu} object.
   */
  public MainMenu getMainMenu() {
    return mainMenu;
  }

  /**
   * Getter for the field <code>modelingShell</code>.
   *
   * @return a {@link ffx.ui.ModelingShell} object.
   */
  public ModelingShell getModelingShell() {
    if (modelingShell == null) {
      modelingShell = new ModelingShell(this);
      modelingShell.run();
    }
    return modelingShell;
  }

  /**
   * initialize
   */
  public void initialize() {
    if (init) {
      return;
    }
    init = true;

    String dir =
        System.getProperty(
            "user.dir", FileSystemView.getFileSystemView().getDefaultDirectory().getAbsolutePath());
    setCWD(new File(dir));
    locale = new FFXLocale("en", "US");

    // Create the Root Node
    dataRoot = new MSRoot();
    Border bb = BorderFactory.createEtchedBorder(EtchedBorder.RAISED);
    statusLabel = new JLabel("  ");
    JLabel stepLabel = new JLabel("  ");
    stepLabel.setHorizontalAlignment(JLabel.RIGHT);
    JLabel energyLabel = new JLabel("  ");
    energyLabel.setHorizontalAlignment(JLabel.RIGHT);
    JPanel statusPanel = new JPanel(new GridLayout(1, 3));
    statusPanel.setBorder(bb);
    statusPanel.add(statusLabel);
    statusPanel.add(stepLabel);
    statusPanel.add(energyLabel);

    // Initialize various Panels
    setLayout(new BorderLayout());
    hierarchy = new Hierarchy(this);
    hierarchy.setStatus(statusLabel, stepLabel, energyLabel);
    keywordPanel = new KeywordPanel(this);
    modelingPanel = new ModelingPanel(this);

    JPanel treePane = new JPanel(new BorderLayout());
    JScrollPane scrollPane = new JScrollPane(hierarchy, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
            JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
    treePane.add(scrollPane, BorderLayout.CENTER);

    if (!GraphicsEnvironment.isHeadless()) {
      GraphicsConfigTemplate3D template3D = new GraphicsConfigTemplate3D();
      // template3D.setDoubleBuffer(GraphicsConfigTemplate.PREFERRED);
      GraphicsConfiguration gc = null;
      try {
        gc = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getBestConfiguration(template3D);
      } catch (Exception e) {
        logger.log(Level.SEVERE, " Exception encountered when trying to get the best GraphicsConfiguration", e);
      }
      try {
        graphicsCanvas = new GraphicsCanvas(gc, this);
      } catch (Exception e) {
        logger.log(Level.SEVERE, " Exception encountered when trying to create the GraphicsCanvas", e);
      }
      graphicsPanel = new GraphicsPanel(graphicsCanvas, statusPanel);
    }

    // Holds 3D Graphics, Keyword Editor, Modeling Commands and Log Panels
    // JTabbedPane tabbedPane = new JTabbedPane();
    // ClassLoader loader = getClass().getClassLoader();
    // ImageIcon graphicsIcon = new ImageIcon(loader.getResource("ffx/ui/icons/monitor.png"));
    // ImageIcon keywordIcon = new ImageIcon(loader.getResource("ffx/ui/icons/key.png"));
    // ImageIcon modelingIcon = new ImageIcon(loader.getResource("ffx/ui/icons/cog.png"));
    // tabbedPane.addTab(locale.getValue("Graphics"), graphicsIcon, graphicsPanel);
    // tabbedPane.addTab(locale.getValue("KeywordEditor"), keywordIcon, keywordPanel);
    // tabbedPane.addTab(locale.getValue("ModelingCommands"), modelingIcon, modelingPanel);
    // tabbedPane.addChangeListener(this);
    // splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, false, treePane, tabbedPane);

    splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, false, treePane, graphicsPanel);
    splitPane.setResizeWeight(0.25);
    splitPane.setOneTouchExpandable(true);
    add(splitPane, BorderLayout.CENTER);

    if (!GraphicsEnvironment.isHeadless()) {
      mainMenu = new MainMenu(this);
      add(mainMenu.getToolBar(), BorderLayout.NORTH);
      getModelingShell();
      loadPrefs();
    }
  }

  /**
   * merge
   *
   * @param nodesToMerge an array of {@link ffx.potential.bonded.MSNode} objects.
   */
  public void merge(MSNode[] nodesToMerge) {
    ArrayList<MSNode> activeNodes = new ArrayList<>();
    for (MSNode node : nodesToMerge) {
      if (node != null) {
        activeNodes.add(node);
      }
    }
    if (activeNodes.size() > 1) {
      merge(activeNodes);
    }
  }

  public Thread open(File file, String commandDescription) {
    UIFileOpener opener = openInit(file, commandDescription);
    openThread = new Thread(opener);
    openThread.start();
    setPanel(GRAPHICS);
    return openThread;
  }

  public Thread open(List<File> files, String commandDescription) {
    UIFileOpener openFile = openInit(files, commandDescription);
    openThread = new Thread(openFile);
    openThread.start();
    setPanel(GRAPHICS);
    return openThread;
  }

  /**
   * open
   *
   * @param name a {@link java.lang.String} object.
   * @return a {@link java.lang.Thread} object.
   */
  public Thread open(String name) {
    File file = resolveName(name);
    if (file == null) {
      logger.log(Level.WARNING, "{0}: could not be found.", name);
      return null;
    }
    return open(file, null);
  }

  /**
   * open
   *
   * @param names an array of {@link java.lang.String} objects.
   * @return a {@link java.lang.Thread} object.
   */
  public Thread open(String[] names) {
    if (names == null) {
      return null;
    }
    List<File> files = new ArrayList<>();
    // Resolve all file names.
    for (String name : names) {
      File file = resolveName(name);
      if (file == null || !file.exists()) {
        return null;
      }
      files.add(file);
    }
    return open(files, null);
  }

  public synchronized MolecularAssembly[] openWaitUtils(String file) {
    UIFileOpener opener = openFromUtils(file);
    Thread thread = new Thread(opener);
    while (thread.isAlive()) {
      try {
        wait(1);
      } catch (InterruptedException e) {
        String message = "Exception waiting for " + file + " to open.";
        logger.log(Level.WARNING, message, e);
        return null;
      }
    }

    MolecularAssembly[] systems = activeFilter.getMolecularAssemblyArray();
    if (systems != null) {
      int n = systems.length;
      FFXSystem[] ffxSystems = new FFXSystem[n];
      FFXSystem[] allSystems = getHierarchy().getSystems();
      int total = allSystems.length;
      System.arraycopy(allSystems, total - n, ffxSystems, 0, n);
      return ffxSystems;
    } else {
      return null;
    }
  }

  /**
   * reset
   */
  public void reset() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    trajectory.stop();
    trajectory.rewind();
  }

  public void saveKeywordFile(File file) {
    keywordPanel.keySave(file);
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void stateChanged(ChangeEvent evt) {
    JTabbedPane jtp = (JTabbedPane) evt.getSource();
    int index = jtp.getSelectedIndex();
    if (index == 0) {
      graphicsCanvas.selected();
    } else if (index == 1) {
      keywordPanel.selected();
    } else if (index == 2) {
      modelingPanel.selected();
    }
  }

  /**
   * stop
   */
  public void stop() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    trajectory.stop();
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return "Program Control";
  }

  /**
   * Prompt the user to select an alternate key file.
   */
  private void chooseKey() {
    JFileChooser d = MainPanel.resetFileChooser();
    d.setFileFilter(new KeyFileFilter());
    d.setAcceptAllFileFilterUsed(false);
    FFXSystem sys = getHierarchy().getActive();
    if (sys != null) {
      File newCWD = sys.getFile();
      if (newCWD != null && newCWD.getParentFile() != null) {
        d.setCurrentDirectory(newCWD.getParentFile());
      }
    } else {
      return;
    }
    int result = d.showOpenDialog(this);
    if (result == JFileChooser.APPROVE_OPTION) {
      File f = d.getSelectedFile();
      sys.setKeyFile(f);
      sys.setKeywords(KeyFilter.open(f));
      getKeywordPanel().loadActive(sys);
    }
  }

  /**
   * Prompt the user to select an alternate log file
   */
  private void chooseLog() {
    JFileChooser d = resetFileChooser();
    FFXSystem sys = getHierarchy().getActive();
    if (sys != null) {
      File newCWD = sys.getFile();
      if (newCWD != null && newCWD.getParentFile() != null) {
        d.setCurrentDirectory(newCWD.getParentFile());
      }
    } else {
      return;
    }
    d.setDialogTitle("Select a log file");
    d.setAcceptAllFileFilterUsed(true);
    int result = d.showOpenDialog(this);
    if (result == JFileChooser.APPROVE_OPTION) {
      File f = d.getSelectedFile();
      if (f != null) {
        sys.setLogFile(f);
        setCWD(d.getCurrentDirectory());
        // getModelingPanel().selected();
      }
    }
  }

  /**
   * closeWait
   */
  synchronized void closeWait() {
    FFXSystem active = hierarchy.getActive();
    if (active == null) {
      logger.log(Level.INFO, " No active system to close.");
      return;
    }
    Thread thread = close(active);
    while (thread != null && thread.isAlive()) {
      try {
        wait(1);
      } catch (InterruptedException e) {
        String message = "Exception waiting for " + active + " to close.";
        logger.log(Level.WARNING, message, e);
      }
    }
  }

  /**
   * Close all open systems.
   */
  synchronized void closeAll() {
    while (hierarchy.getActive() != null) {
      closeWait();
    }
  }

  /**
   * Attempt to connect to a TINKER Simulation
   */
  private void connect() {
    if (simulation == null || simulation.isFinished()) {
      if (simulation != null) {
        simulation.release();
      }
      simulation = new SimulationLoader(null, null, this, socketAddress);
      simulation.connect();
      mainMenu.setConnect(false);
      setPanel(GRAPHICS);
    }
  }

  /**
   * connectToTINKER
   *
   * @param system         a {@link ffx.ui.FFXSystem} object.
   * @param modelingThread a {@link java.lang.Thread} object.
   */
  void connectToTINKER(FFXSystem system, Thread modelingThread) {
    if (simulation == null || simulation.isFinished()) {
      if (simulation != null) {
        simulation.release();
      }
      InetSocketAddress tempAddress;
      try {
        tempAddress = new InetSocketAddress(InetAddress.getLocalHost(), port);
      } catch (Exception e) {
        try {
          tempAddress = new InetSocketAddress(InetAddress.getByName(null), port);
        } catch (Exception ex) {
          System.err.println("Could not determine Local Host: " + ex);
          return;
        }
      }
      simulation = new SimulationLoader(system, modelingThread, this, tempAddress);
      if (modelingThread != null) {
        modelingThread.start();
      }
      simulation.connect();
      mainMenu.setConnect(false);
      setPanel(GRAPHICS);
    }
  }

  /**
   * createKeyFile
   *
   * @param system a {@link ffx.ui.FFXSystem} object.
   * @return a boolean.
   */
  boolean createKeyFile(FFXSystem system) {
    String message = "Please select a parameter file " + "and a TINKER Key file will be created.";
    String params = (String) JOptionPane.showInputDialog(this, message, "Parameter File",
        JOptionPane.QUESTION_MESSAGE, null, keywordPanel.getParamFiles(), null);
    if (params != null) {
      if (params.equalsIgnoreCase("Use an existing TINKER Key file")) {
        JFileChooser fc = resetFileChooser();
        fc.setDialogTitle("Choose a KEY File");
        fc.setCurrentDirectory(pwd);
        fc.setSelectedFile(null);
        fc.setFileFilter(keyFileFilter);
        int result = fc.showOpenDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
          File keyfile = fc.getSelectedFile();
          if (keyfile.exists()) {
            Hashtable<String, Keyword> keywordHash = KeyFilter.open(keyfile);
            if (keywordHash != null) {
              system.setKeywords(keywordHash);
            } else {
              return false;
            }
            system.setKeyFile(keyfile);
            system.setForceField(null);
            return true;
          }
        }
      } else {
        File tempFile = system.getFile();
        if (tempFile.getParentFile().canWrite()) {
          String path = system.getFile().getParent() + File.separatorChar;
          String keyFileName = system.getName() + ".key";
          File keyfile = new File(path + keyFileName);
          try {
            FileWriter fw = new FileWriter(keyfile);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("\n");
            bw.write("# Force Field Selection\n");
            String tempParm = keywordPanel.getParamPath(params);
            if (tempParm.indexOf(" ") > 0) {
              tempParm = "\"" + keywordPanel.getParamPath(params) + "\"";
            }
            bw.write("PARAMETERS        " + tempParm + "\n");
            bw.close();
            fw.close();
            Hashtable<String, Keyword> keywordHash = KeyFilter.open(keyfile);
            if (keywordHash != null) {
              system.setKeywords(keywordHash);
            } else {
              return false;
            }
            system.setKeyFile(keyfile);
            system.setForceField(null);
            return true;
          } catch (Exception e) {
            logger.warning("" + e);
            message = "There was an error creating " + keyfile.getAbsolutePath();
            JOptionPane.showMessageDialog(this, message);
          }
        } else {
          message =
              "Could not create a Key file because " + pwd.getAbsolutePath() + " is not writable";
          JOptionPane.showMessageDialog(this, message);
        }
      }
    }
    return false;
  }

  /**
   * exit with a target ExitStatus
   *
   * @param exitStatus How FFX has closed.
   */
  void exit(ExitStatus exitStatus) {
    // Package-private out of conservatism; may be safe to make public.
    savePrefs();

    Resources.logResources();

    System.exit(exitStatus.getExitCode());
  }

  /**
   * Set the current exit code.
   *
   * @param exitType Enumerated type for exit codes.
   */
  void setExitType(ExitStatus exitType) {
    // Package-private out of conservatism; may be safe to make public.
    this.exitType = exitType;
  }

  /**
   * Getter for the field <code>dataRoot</code>.
   *
   * @return a {@link ffx.potential.bonded.MSRoot} object.
   */
  MSRoot getDataRoot() {
    return dataRoot;
  }

  /**
   * getFFXLocale
   *
   * @return a {@link ffx.ui.properties.FFXLocale} object.
   */
  FFXLocale getFFXLocale() {
    return locale;
  }

  /**
   * getGraphics3D
   *
   * @return a {@link ffx.ui.GraphicsCanvas} object.
   */
  GraphicsCanvas getGraphics3D() {
    return graphicsCanvas;
  }

  /**
   * Getter for the field <code>keywordPanel</code>.
   *
   * @return a {@link ffx.ui.KeywordPanel} object.
   */
  KeywordPanel getKeywordPanel() {
    return keywordPanel;
  }

  ModelingPanel getModelingPanel() {
    return modelingPanel;
  }

  /**
   * getStatusBar
   *
   * @return a {@link javax.swing.JLabel} object.
   */
  JLabel getStatusBar() {
    return statusLabel;
  }

  /**
   * Get the Trajectory wrapper for the active system
   *
   * @return trajectory
   */
  private Trajectory getTrajectory() {
    FFXSystem system = hierarchy.getActive();
    if (system == null) {
      return null;
    }
    Trajectory trajectory = system.getTrajectory();
    if (trajectory != null) {
      return trajectory;
    }
    trajectory = new Trajectory(system, this);
    trajectory.setOscillate(oscillate);
    system.setTrajectory(trajectory);
    return trajectory;
  }

  /**
   * highlightSelections
   *
   * @param evt a {@link java.awt.event.ActionEvent} object.
   */
  private void highlightSelections(ActionEvent evt) {
    if (evt.getSource() instanceof JCheckBoxMenuItem jcb) {
      hierarchy.setHighlighting(jcb.isSelected());
    } else {
      boolean highlighting = RendererCache.highlightSelections;
      if (highlighting) {
        hierarchy.setHighlighting(false);
        mainMenu.setHighlighting(false);
      } else {
        hierarchy.setHighlighting(true);
        mainMenu.setHighlighting(true);
      }
    }
  }

  private void initAbout() {
    aboutTextArea = new JTextArea();
    Font font = Font.decode(Font.MONOSPACED);
    aboutTextArea.setFont(font);
    aboutTextArea.setText(aboutString);
    aboutTextArea.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
    aboutTextArea.setEditable(false);
  }

  /**
   * isOpening
   *
   * @return a boolean.
   */
  boolean isOpening() {
    return (openThread != null && openThread.isAlive());
  }

  /**
   * Load preferences from the user node
   */
  private void loadPrefs() {
    String c = MainPanel.class.getName();
    JFrame frame1 = (JFrame) SwingUtilities.getRoot(this);
    Toolkit toolkit = getToolkit();
    Dimension screenSize = toolkit.getScreenSize();
    int x = preferences.getInt(c + ".x", screenSize.width / 8);
    int y = preferences.getInt(c + ".y", screenSize.height / 8);
    int width = preferences.getInt(c + ".width", screenSize.width * 3 / 4);
    int height = preferences.getInt(c + ".height", screenSize.height * 3 / 4);
    if (width > screenSize.width * 0.4
        && width < screenSize.width * 0.8
        && height > screenSize.height * 0.4
        && height < screenSize.height * 0.8) {
      frame1.setSize(width, height);
    } else {
      frame1.setSize(screenSize.width * 4 / 5, screenSize.height * 4 / 5);
    }
    if (x > 0 && x < screenSize.width / 2 && y > 0 && y < screenSize.height / 2) {
      frame1.setLocation(x, y);
    } else {
      frame1.setLocation(screenSize.width / 8, screenSize.height / 8);
    }
    splitPaneDivider = preferences.getInt(c + ".divider", 200);
    if (splitPaneDivider < frame1.getWidth() * (1.0f / 4.0f)) {
      splitPaneDivider = (int) (frame1.getWidth() * (1.0f / 4.0f));
    }
    splitPane.setDividerLocation(splitPaneDivider);
    if (!preferences.getBoolean(c + ".system", true)) {
      mainMenu.setSystemShowing(false);
      splitPane.setDividerLocation(0);
    } else {
      mainMenu.setSystemShowing(true);
    }
    if (!preferences.getBoolean(c + ".menu", true)) {
      remove(mainMenu.getToolBar());
      mainMenu.setMenuShowing(false);
      validate();
    } else {
      mainMenu.setMenuShowing(true);
    }
    try {
      port = preferences.getInt(c + ".port", 2000);
      ip = preferences.get(c + ".ip", InetAddress.getLocalHost().getHostAddress());
      if (ip != null) {
        address = InetAddress.getByName(ip);
        socketAddress = new InetSocketAddress(address, port);
      } else {
        socketAddress = new InetSocketAddress(port);
      }
    } catch (Exception e) {
      logger.log(Level.WARNING, e.toString());
    }
    if (graphicsCanvas != null) {
      graphicsCanvas.loadPrefs();
    }
  }

  /**
   * merge
   */
  private void merge() {
    ArrayList<MSNode> activeNodes = hierarchy.getActiveNodes();
    if (activeNodes.size() >= 2) {
      merge(activeNodes);
    }
  }

  /**
   * Merge two or more selected FSystem Nodes into one FSystem node. There are a few gotchas that
   * need to be fixed
   *
   * @param nodesToMerge a {@link java.util.ArrayList} object.
   */
  private void merge(ArrayList<MSNode> nodesToMerge) {
    ArrayList<MSNode> activeNodes = new ArrayList<>();
    for (MSNode node : nodesToMerge) {
      if (node != null && !(node instanceof MSRoot)) {
        activeNodes.add(node);
      }
    }
    if (activeNodes.size() <= 1) {
      return;
    }
    // Set up a structure to hold the new system
    FFXSystem active = hierarchy.getActive();
    File file = SystemFilter.version(hierarchy.getActive().getFile());
    FFXSystem system = new FFXSystem(file, "Merge Result", active.getProperties());
    system.setKeyFile(active.getKeyFile());
    system.setKeywords(KeyFilter.open(active.getKeyFile()));
    // Fill arrays with the atoms and bonds from the systems to be combined
    ArrayList<Atom> mergedAtoms = new ArrayList<>();
    ArrayList<Bond> mergedBonds = new ArrayList<>();
    ArrayList<FFXSystem> systems = new ArrayList<>();
    TransformGroup parentTransformGroup;
    FFXSystem parentSystem;
    Transform3D parentTransform3D = new Transform3D();
    Vector3d parentPosition = new Vector3d();
    Vector3d atomPosition = new Vector3d();
    // TINKER Atom Numbers start at 1
    int atomNum = 1;
    Vector3d zero = new Vector3d(0.0, 0.0, 0.0);
    for (MSNode m : activeNodes) {
      parentSystem = m.getMSNode(FFXSystem.class);
      if (parentSystem == null) {
        return;
      }
      if (!systems.contains(parentSystem)) {
        graphicsCanvas.updateSceneWait(
            parentSystem, false, true, RendererCache.ViewModel.WIREFRAME, false, null);
        systems.add(parentSystem);
      }
      // Move each atom into the global frame by applying the System
      // Transform to
      // relative atomic position
      parentTransformGroup = parentSystem.getOriginToRot();
      parentTransformGroup.getTransform(parentTransform3D);
      parentTransform3D.get(parentPosition);
      parentTransform3D.setTranslation(zero);
      // parentTransform3D.setScale(1.0d);
      for (Atom atom : m.getAtomList()) {
        atom.removeFromParent();
        atom.setXyzIndex(atomNum++);
        mergedAtoms.add(atom);
        atom.getV3D(atomPosition);
        parentTransform3D.transform(atomPosition);
        atomPosition.add(parentPosition);
        atom.moveTo(atomPosition);
      }
      for (Bond bond : m.getBondList()) {
        bond.removeFromParent();
        mergedBonds.add(bond);
      }
    }
    for (FFXSystem sys : systems) {
      close(sys);
    }
    MergeFilter mergeFilter = new MergeFilter(system, mergedAtoms, mergedBonds);
    UIFileOpener fileOpener = new UIFileOpener(mergeFilter, this);
    if (fileOpenerThreads > 0) {
      fileOpener.setNThreads(fileOpenerThreads);
    }
    Thread thread = new Thread(fileOpener);
    thread.start();
  }

  /**
   * oceanLookAndFeel
   */
  private void oceanLookAndFeel() {
    try {
      UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
      SwingUtilities.updateComponentTreeUI(SwingUtilities.getRoot(this));
    } catch (Exception e) {
      //
    }
  }

  /**
   * Trys to convert a file picked from a JFileChooser
   */
  private Thread open() {
    if (openThread != null && openThread.isAlive()) {
      return null;
    }
    JFileChooser fc = resetFileChooser();
    fc.setDialogTitle("Choose FFX, PDB, XYZ or ARC");
    fc.addChoosableFileFilter(xyzFileFilter);
    fc.addChoosableFileFilter(pdbFileFilter);
    fc.addChoosableFileFilter(intFileFilter);
    fc.addChoosableFileFilter(arcFileFilter);
    fc.addChoosableFileFilter(ffxFileFilter);
    fc.setAcceptAllFileFilterUsed(true);
    int result = fc.showOpenDialog(this);
    if (result == JFileChooser.APPROVE_OPTION) {
      File file = fc.getSelectedFile();
      return open(file, null);
    }
    return null;
  }

  private UIFileOpener openFromUtils(File file, String commandDescription) {
    UIFileOpener opener = openInit(file, commandDescription);
    openThread = new Thread(opener);
    openThread.start();
    setPanel(GRAPHICS);
    return opener;
  }

  /**
   * Attempts to load the supplied file
   *
   * @param file               File to open
   * @param commandDescription Description of the command that created this file.
   * @return a {@link java.lang.Thread} object.
   */
  private UIFileOpener openInit(File file, String commandDescription) {
    if (file == null || !file.isFile() || !file.canRead()) {
      return null;
    }
    file = new File(FilenameUtils.normalize(file.getAbsolutePath()));
    // Set the Current Working Directory based on this file.
    setCWD(file.getParentFile());

    // Create the CompositeConfiguration properties.
    CompositeConfiguration properties = Keyword.loadProperties(file);
    // Create an FFXSystem for this file.
    FFXSystem newSystem = new FFXSystem(file, commandDescription, properties);
    // Create a Force Field.
    forceFieldFilter = new ForceFieldFilter(properties);
    ForceField forceField = forceFieldFilter.parse();
    String[] patches = properties.getStringArray("patch");
    for (String patch : patches) {
      logger.info(" Attempting to read force field patch from " + patch + ".");
      CompositeConfiguration patchConfiguration = new CompositeConfiguration();
      patchConfiguration.addProperty("parameters", patch);
      forceFieldFilter = new ForceFieldFilter(patchConfiguration);
      ForceField patchForceField = forceFieldFilter.parse();
      forceField.append(patchForceField);
      if (RotamerLibrary.addRotPatch(patch)) {
        logger.info(format(" Loaded rotamer definitions from patch %s.", patch));
      }
    }
    newSystem.setForceField(forceField);
    SystemFilter systemFilter;

    // Decide what parser to use.
    if (xyzFileFilter.acceptDeep(file)) {
      // Use the TINKER Cartesian Coordinate File Parser.
      systemFilter = new XYZFilter(file, newSystem, forceField, properties);
    } else if (intFileFilter.acceptDeep(file)) {
      // Use the TINKER Internal Coordinate File Parser.
      systemFilter = new INTFilter(file, newSystem, forceField, properties);
    } else {
      // Use the PDB File Parser.
      systemFilter = new PDBFilter(file, newSystem, forceField, properties);
    }

    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    activeFilter = systemFilter;
    UIFileOpener fileOpener = new UIFileOpener(systemFilter, this);
    if (fileOpenerThreads > 0) {
      fileOpener.setNThreads(fileOpenerThreads);
    }
    return fileOpener;
  }

  private UIFileOpener openFromUtils(List<File> files, String commandDescription) {
    UIFileOpener openFile = openInit(files, commandDescription);
    openThread = new Thread(openFile);
    openThread.start();
    setPanel(GRAPHICS);
    return openFile;
  }

  /**
   * Attempts to load the supplied file
   *
   * @param files              Files to open
   * @param commandDescription Description of the command that created this file.
   * @return a {@link java.lang.Thread} object.
   */
  private UIFileOpener openInit(List<File> files, String commandDescription) {
    if (files == null) {
      return null;
    }
    File file = new File(FilenameUtils.normalize(files.get(0).getAbsolutePath()));
    // Set the Current Working Directory based on this file.
    setCWD(file.getParentFile());

    // Create the CompositeConfiguration properties.
    CompositeConfiguration properties = Keyword.loadProperties(file);
    forceFieldFilter = new ForceFieldFilter(properties);
    ForceField forceField = forceFieldFilter.parse();

    // Create an FFXSystem for this file.
    FFXSystem newSystem = new FFXSystem(file, commandDescription, properties);
    String[] patches = properties.getStringArray("patch");
    for (String patch : patches) {
      logger.info(" Attempting to read force field patch from " + patch + ".");
      CompositeConfiguration patchConfiguration = new CompositeConfiguration();
      patchConfiguration.addProperty("parameters", patch);
      forceFieldFilter = new ForceFieldFilter(patchConfiguration);
      ForceField patchForceField = forceFieldFilter.parse();
      forceField.append(patchForceField);
      if (RotamerLibrary.addRotPatch(patch)) {
        logger.info(format(" Loaded rotamer definitions from patch %s.", patch));
      }
    }
    newSystem.setForceField(forceField);
    // Decide what parser to use.
    SystemFilter systemFilter;
    if (xyzFileFilter.acceptDeep(file)) {
      // Use the TINKER Cartesian Coordinate File Parser.
      systemFilter = new XYZFilter(files, newSystem, forceField, properties);
    } else if (intFileFilter.acceptDeep(file)) {
      // Use the TINKER Internal Coordinate File Parser.
      systemFilter = new INTFilter(files, newSystem, forceField, properties);
    } else {
      // Use the PDB File Parser.
      systemFilter = new PDBFilter(files, newSystem, forceField, properties);
    }

    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    activeFilter = systemFilter;
    // return new UIFileOpener(systemFilter, this);
    UIFileOpener fileOpener = new UIFileOpener(systemFilter, this);
    if (fileOpenerThreads > 0) {
      fileOpener.setNThreads(fileOpenerThreads);
    }
    return fileOpener;
  }

  /**
   * openWait
   *
   * @param file a {@link java.lang.String} object.
   * @return an array of {@link ffx.ui.FFXSystem} objects.
   */
  synchronized FFXSystem[] openWait(String file) {
    Thread thread = open(file);
    while (thread != null && thread.isAlive()) {
      try {
        wait(1);
      } catch (InterruptedException e) {
        String message = "Exception waiting for " + file + " to open.";
        logger.log(Level.WARNING, message, e);
        return null;
      }
    }

    MolecularAssembly[] systems = activeFilter.getMolecularAssemblyArray();
    if (systems != null) {
      int n = systems.length;
      FFXSystem[] ffxSystems = new FFXSystem[n];
      FFXSystem[] allSystems = getHierarchy().getSystems();
      int total = allSystems.length;
      arraycopy(allSystems, total - n, ffxSystems, 0, n);
      return ffxSystems;
    } else {
      return null;
    }
  }

  /**
   * openWait
   *
   * @param file     a {@link java.lang.String} object.
   * @param nThreads the number of threads.
   * @return an array of {@link ffx.ui.FFXSystem} objects.
   */
  synchronized FFXSystem[] openWait(String file, int nThreads) {
    fileOpenerThreads = nThreads;
    FFXSystem[] systs = openWait(file);
    fileOpenerThreads = -1;
    return systs;
  }

  /**
   * openWait
   *
   * @param files an array of {@link java.lang.String} objects.
   * @return an array of {@link ffx.ui.FFXSystem} objects.
   */
  synchronized FFXSystem[] openWait(String[] files) {
    Thread thread = open(files);
    while (thread != null && thread.isAlive()) {
      try {
        wait(1);
      } catch (InterruptedException e) {
        String message = "Exception waiting for " + files[0] + " to open.";
        logger.log(Level.WARNING, message, e);
        return null;
      }
    }

    MolecularAssembly[] systems = activeFilter.getMolecularAssemblyArray();
    if (systems != null) {
      int n = systems.length;
      FFXSystem[] ffxSystems = new FFXSystem[n];
      FFXSystem[] allSystems = getHierarchy().getSystems();
      int total = allSystems.length;
      arraycopy(allSystems, total - n, ffxSystems, 0, n);
      return ffxSystems;
    } else {
      return null;
    }
  }

  /**
   * openWait
   *
   * @param files    an array of {@link java.lang.String} objects.
   * @param nThreads the number of threads.
   * @return an array of {@link ffx.ui.FFXSystem} objects.
   */
  synchronized FFXSystem[] openWait(String[] files, int nThreads) {
    fileOpenerThreads = nThreads;
    FFXSystem[] systs = openWait(files);
    fileOpenerThreads = -1;
    return systs;
  }

  private UIFileOpener openFromUtils(String name) {
    File file = resolveName(name);
    if (file == null) {
      logger.log(Level.WARNING, "{0}: could not be found.", name);
      return null;
    }
    return openFromUtils(file, null);
  }

  private File resolveName(String name) {
    // Return null if name == null.
    if (name == null) {
      return null;
    }
    File file = new File(name);
    // If the file exists, return it.
    if (file.exists()) {
      return file;
    }
    // Check for a file in the CWD.
    file = new File(pwd + File.separator + name);
    if (file.exists()) {
      return file;
    }
    // Check for an HTTP address
    if (name.startsWith("http://")) {
      String fileName = FilenameUtils.getName(name);
      if (fileName == null) {
        return null;
      }
      return downloadURL(name);
    }
    // Check for a PDB ID.
    if (name.length() == 4) {
      String fileName = name + ".pdb";
      String path = getPWD().getAbsolutePath();
      File pdbFile = new File(path + File.separatorChar + fileName);
      if (!pdbFile.exists()) {
        String fromURL = pdbForID(name);
        return downloadURL(fromURL);
      } else {
        return pdbFile;
      }
    }
    return null;
  }

  private UIFileOpener openFromUtils(String[] names) {
    if (names == null) {
      return null;
    }
    List<File> files = new ArrayList<>();
    // Resolve all file names.
    for (String name : names) {
      File file = resolveName(name);
      if (file == null || !file.exists()) {
        return null;
      }
      files.add(file);
    }
    return openFromUtils(files, null);
  }

  /**
   * Opens a file from the PDB
   */
  private void openFromPDB() {
    if (openThread != null && openThread.isAlive()) {
      return;
    }
    String code = JOptionPane.showInputDialog("Enter the PDB Identifier (4 characters)", "");
    if (code == null) {
      return;
    }
    code = code.toLowerCase().trim();
    if (code.length() != 4) {
      return;
    }
    String fileName = code + ".pdb";
    String path = getPWD().getAbsolutePath();
    File pdbFile = new File(path + File.separatorChar + fileName);
    CompositeConfiguration properties = Keyword.loadProperties(pdbFile);
    forceFieldFilter = new ForceFieldFilter(properties);
    ForceField forceField = forceFieldFilter.parse();
    FFXSystem newSystem = new FFXSystem(pdbFile, "PDB", properties);
    newSystem.setForceField(forceField);
    if (!pdbFile.exists()) {
      String fromURL = pdbForID(code);
      pdbFile = downloadURL(fromURL);
      if (pdbFile == null || !pdbFile.exists()) {
        return;
      }
    } else {
      String message = format(" Reading the local copy of the PDB file %s.", pdbFile);
      logger.info(message);
    }
    PDBFilter pdbFilter = new PDBFilter(pdbFile, newSystem, forceField, properties);
    setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    UIFileOpener openFile = new UIFileOpener(pdbFilter, this);
    if (fileOpenerThreads > 0) {
      openFile.setNThreads(fileOpenerThreads);
    }
    openThread = new Thread(openFile);
    openThread.start();
    setPanel(GRAPHICS);
  }

  private File downloadURL(String fromString) {
    // Check for null input.
    if (fromString == null) {
      return null;
    }

    // Convert the string to a URL instance.
    URL fromURL;
    try {
      URI uri = new URI(fromString);
      fromURL = uri.toURL();
    } catch (MalformedURLException | URISyntaxException e) {
      String message = format(" URL incorrectly formatted %s.", fromString);
      logger.log(Level.INFO, message, e);
      return null;
    }

    // Download the URL to a local file.
    logger.info(format(" Downloading %s", fromString));
    try {
      File toFile = new File(FilenameUtils.getName(fromURL.getPath()));
      FileUtils.copyURLToFile(fromURL, toFile, 1000, 1000);
      logger.info(format(" Saved to %s\n", toFile.getPath()));
      return toFile;
    } catch (IOException ex) {
      logger.log(Level.INFO, " Failed to read URL " + fromURL.getPath(), ex);
      return null;
    }
  }

  private void openInduced() {
    FFXSystem active = hierarchy.getActive();
    resetFileChooser();
    fileChooser.setCurrentDirectory(pwd);
    fileChooser.setSelectedFile(active.getFile());
    fileChooser.setDialogTitle("Choose Induced Dipole File");
    fileChooser.addChoosableFileFilter(indFileFilter);
    fileChooser.setAcceptAllFileFilterUsed(true);
    fileChooser.setFileFilter(indFileFilter);
    int result = fileChooser.showOpenDialog(this);
    if (result == JFileChooser.APPROVE_OPTION) {
      File f = fileChooser.getSelectedFile();
      InducedFilter indFilter = new InducedFilter(active, f);
      indFilter.read();
    }
  }

  /**
   * Attempt to convert a TINKER *.key file
   *
   * @param newSystem FFXSystem that needs an associated Key File
   * @param createKey flag to create a key file be created
   * @return Key file that was found, or null if nothing could be found
   */
  boolean openKey(FFXSystem newSystem, boolean createKey) {
    String keyFileName;
    String temp = newSystem.getFile().getName();
    int dot = temp.lastIndexOf(".");
    if (dot > 0) {
      keyFileName = temp.substring(0, dot) + ".key";
    } else {
      keyFileName = temp + ".key";
    }
    String path = newSystem.getFile().getParent() + File.separator;
    File keyfile = new File(path + keyFileName);
    if (keyfile.exists()) {
      Hashtable<String, Keyword> keywordHash = KeyFilter.open(keyfile);
      if (keywordHash != null) {
        newSystem.setKeywords(keywordHash);
      } else {
        return false;
      }
      newSystem.setKeyFile(keyfile);
      newSystem.setForceField(null);
      return true;
    }
    keyfile = new File(path + "tinker.key");
    if (keyfile.exists()) {
      logger.info("Using tinker.key: " + keyfile);
      Hashtable<String, Keyword> keywordHash = KeyFilter.open(keyfile);
      if (keywordHash != null) {
        newSystem.setKeywords(keywordHash);
      } else {
        return false;
      }
      newSystem.setKeyFile(keyfile);
      newSystem.setForceField(null);
      return true;
    }
    if (createKey) {
      return createKeyFile(newSystem);
    }
    return false;
  }

  /**
   * openOn
   *
   * @param f         a {@link java.io.File} object.
   * @param oldSystem a {@link ffx.ui.FFXSystem} object.
   * @param command   a {@link java.lang.String} object.
   */
  void openOn(File f, FFXSystem oldSystem, String command) {
    XYZFilter.readOnto(f, oldSystem);
    oldSystem.setCommandDescription(command);
    graphicsCanvas.updateScene(oldSystem, true, false, null, false, null);
    getHierarchy().updateStatus();
    getHierarchy().repaint();
  }

  /**
   * oscillate
   *
   * @param evt a {@link java.awt.event.ActionEvent} object.
   */
  private void oscillate(ActionEvent evt) {
    oscillate = ((JCheckBoxMenuItem) evt.getSource()).isSelected();
    FFXSystem[] systems = getHierarchy().getSystems();

    if (systems == null) {
      return;
    }

    for (FFXSystem system : systems) {
      Trajectory trajectory = system.getTrajectory();
      if (trajectory != null) {
        trajectory.setOscillate(oscillate);
      }
    }
  }

  /**
   * platformLookAndFeel
   */
  private void platformLookAndFeel() {
    try {
      if (SystemUtils.IS_OS_LINUX) {
        UIManager.setLookAndFeel("com.sun.java.swing.plaf.motif.MotifLookAndFeel");
      } else {
        UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
      }
      SwingUtilities.updateComponentTreeUI(SwingUtilities.getRoot(this));
    } catch (ClassNotFoundException
             | InstantiationException
             | IllegalAccessException
             | UnsupportedLookAndFeelException e) {
      logger.log(Level.WARNING, "Can''t set look and feel: {0}", e);
    }
  }

  /**
   * play
   */
  private void play() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    trajectory.start();
  }

  /**
   * Close the connection to a running simulation
   */
  private void release() {
    if (simulation != null) {
      simulation.release();
      simulation = null;
      mainMenu.setConnect(true);
    }
  }

  /**
   * resetPanes
   */
  private void resetPanes() {
    resizePanes(0);
  }

  /**
   * resetShell
   */
  void resetShell() {
    if (!GraphicsEnvironment.isHeadless()) {
      modelingShell = getModelingShell();
      try {
        modelingShell.exit();
      } catch (NullPointerException e) {
        //
      } finally {
        modelingShell = null;
      }
      modelingShell = getModelingShell();
    }
  }

  /**
   * Set the split panes to their default proportions
   *
   * @param move a int.
   */
  private void resizePanes(int move) {
    if (move == 0) {
      splitPaneDivider = 0;
      mainMenu.setMenuShowing(false);
      mainMenu.toggleToolBarShowing();
      mainMenu.setSystemShowing(false);
      mainMenu.systemClick();
    } else {
      splitPane.setDividerLocation(splitPane.getDividerLocation() + move);
    }
  }

  /**
   * Save the currently selected FFXSystem to disk.
   *
   * @param file File to save the system to.
   * @since 1.0
   */
  void saveAsXYZ(File file) {
    FFXSystem system = hierarchy.getActive();
    if (system != null && !system.isClosing()) {
      File saveFile = file;
      if (saveFile == null) {
        resetFileChooser();
        fileChooser.setCurrentDirectory(pwd);
        fileChooser.setFileFilter(xyzFileFilter);
        fileChooser.setAcceptAllFileFilterUsed(false);
        int result = fileChooser.showSaveDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
          saveFile = fileChooser.getSelectedFile();
          pwd = saveFile.getParentFile();
        }
      }
      if (saveFile != null) {
        SystemFilter filter = new XYZFilter(saveFile, system, null, null);
        if (filter.writeFile(saveFile, false)) {
          // Refresh Panels with the new System name
          hierarchy.setActive(system);
          activeFilter = filter;
        }
      }
    }
  }

  /**
   * Save the currently selected FFXSystem to disk.
   *
   * @param file File to save the system to.
   * @since 1.0
   */
  void saveAsP1(File file) {
    FFXSystem system = hierarchy.getActive();
    if (system != null && !system.isClosing()) {
      File saveFile = file;
      if (saveFile == null) {
        resetFileChooser();
        fileChooser.setCurrentDirectory(pwd);
        fileChooser.setFileFilter(xyzFileFilter);
        fileChooser.setAcceptAllFileFilterUsed(false);
        int result = fileChooser.showSaveDialog(this);
        if (result == JFileChooser.APPROVE_OPTION) {
          saveFile = fileChooser.getSelectedFile();
          pwd = saveFile.getParentFile();
        }
      }
      if (saveFile != null) {
        XYZFilter filter = new XYZFilter(saveFile, system, null, null);
        Crystal crystal = system.getCrystal().getUnitCell();
        if (crystal instanceof ReplicatesCrystal) {
          crystal = crystal.getUnitCell();
        }
        if (filter.writeFileAsP1(saveFile, false, crystal)) {
          // Refresh Panels with the new System name
          hierarchy.setActive(system);
        }
        activeFilter = filter;
      }
    }
  }

  /**
   * Save the currently selected FFXSystem to a PDB file.
   *
   * @param file File to save the system to.
   * @since 1.0
   */
  void saveAsPDB(File file) {
    saveAsPDB(file, true);
  }

  /**
   * Save the currently selected FFXSystem to a PDB file.
   *
   * @param file     File to save the system to.
   * @param writeEnd Flag to indicate this is the last snapshot.
   * @since 1.0
   */
  private void saveAsPDB(File file, boolean writeEnd) {
    saveAsPDB(file, writeEnd, false);
  }

  /**
   * Save the currently selected FFXSystem to a PDB file.
   *
   * @param file     File to save the system to.
   * @param writeEnd Flag to indicate this is the last snapshot.
   * @param append   Flag to indicate appending this snapshot.
   * @since 1.0
   */
  void saveAsPDB(File file, boolean writeEnd, boolean append) {
    FFXSystem system = hierarchy.getActive();
    if (system == null) {
      logger.log(Level.INFO, " No active system to save.");
      return;
    }
    if (system.isClosing()) {
      logger.log(Level.INFO, " {0} is being closed and can no longer be saved.", system);
      return;
    }
    File saveFile = file;
    if (saveFile == null) {
      resetFileChooser();
      fileChooser.setCurrentDirectory(pwd);
      fileChooser.setFileFilter(pdbFileFilter);
      fileChooser.setAcceptAllFileFilterUsed(false);
      int result = fileChooser.showSaveDialog(this);
      if (result == JFileChooser.APPROVE_OPTION) {
        saveFile = fileChooser.getSelectedFile();
        pwd = saveFile.getParentFile();
      }
    }
    if (saveFile == null) {
      logger.log(Level.INFO, " No filename is defined for {0}.", system);
      return;
    }
    PDBFilter pdbFilter = new PDBFilter(saveFile, system, null, null);
    if (pdbFilter.writeFile(saveFile, append, false, writeEnd)) {
      // Refresh Panels with the new System name
      hierarchy.setActive(system);
      activeFilter = pdbFilter;
    } else {
      logger.log(Level.INFO, " Save failed for: {0}", system);
    }
  }

  void savePDBasP1(File file) {
    FFXSystem system = hierarchy.getActive();
    if (system == null) {
      logger.log(Level.INFO, " No active system to save.");
      return;
    }
    if (system.isClosing()) {
      logger.log(Level.INFO, " {0} is being closed and can no longer be saved.", system);
      return;
    }
    File saveFile = file;
    if (saveFile == null) {
      resetFileChooser();
      fileChooser.setCurrentDirectory(pwd);
      fileChooser.setFileFilter(pdbFileFilter);
      fileChooser.setAcceptAllFileFilterUsed(false);
      int result = fileChooser.showSaveDialog(this);
      if (result == JFileChooser.APPROVE_OPTION) {
        saveFile = fileChooser.getSelectedFile();
        pwd = saveFile.getParentFile();
      }
    }
    if (saveFile == null) {
      logger.log(Level.INFO, " No filename is defined for {0}.", system);
      return;
    }

    PDBFilter pdbFilter = new PDBFilter(saveFile, system, null, null);
    pdbFilter.writeFileAsP1(saveFile);
  }

  /**
   * saveAsPDB
   *
   * @param activeSystems an array of {@link ffx.potential.MolecularAssembly} objects.
   * @param file          a {@link java.io.File} object.
   */
  void saveAsPDB(MolecularAssembly[] activeSystems, File file) {
    File saveFile = file;
    if (saveFile == null) {
      resetFileChooser();
      fileChooser.setCurrentDirectory(pwd);
      fileChooser.setFileFilter(pdbFileFilter);
      fileChooser.setAcceptAllFileFilterUsed(false);
      int result = fileChooser.showSaveDialog(this);
      if (result == JFileChooser.APPROVE_OPTION) {
        saveFile = fileChooser.getSelectedFile();
        pwd = saveFile.getParentFile();
      }
    }
    if (saveFile != null) {
      PDBFilter pdbFilter = new PDBFilter(saveFile, Arrays.asList(activeSystems), null, null);
      pdbFilter.writeFile(saveFile, false);
      activeFilter = pdbFilter;
    }
  }

  /**
   * Save preferences to the user node
   */
  private void savePrefs() {
    String c = MainPanel.class.getName();
    if (!GraphicsEnvironment.isHeadless()) {
      preferences.putInt(c + ".x", frame.getLocation().x);
      preferences.putInt(c + ".y", frame.getLocation().y);
      preferences.putInt(c + ".width", frame.getWidth());
      preferences.putInt(c + ".height", frame.getHeight());
      preferences.putBoolean(c + ".system", mainMenu.isSystemShowing());
      preferences.putInt(c + ".divider", splitPane.getDividerLocation());
      preferences.putBoolean(c + ".menu", mainMenu.isMenuShowing());
      preferences.putBoolean(c + ".axis", mainMenu.isAxisShowing());
    }
    if (ip == null) {
      ip = "";
    }
    if (address != null) {
      String s = address.getHostAddress();
      if (s != null) {
        preferences.put(c + ".ip", s);
      }
      preferences.putInt(c + ".port", socketAddress.getPort());
    }
    preferences.put(c + ".cwd", pwd.toString());

    if (modelingPanel != null) {
      modelingPanel.savePrefs();
    }
    if (keywordPanel != null) {
      keywordPanel.savePrefs();
    }
    if (graphicsCanvas != null) {
      graphicsCanvas.savePrefs();
    }
  }

  /**
   * selectAll
   */
  private void selectAll() {
    if (dataRoot.getChildCount() == 0) {
      return;
    }
    hierarchy.selectAll();
  }

  /**
   * setCWD
   *
   * @param file a {@link java.io.File} object.
   */
  void setCWD(File file) {
    if ((file == null) || (!file.exists())) {
      return;
    }
    pwd = file;
  }

  /**
   * setPanel
   *
   * @param panel a int.
   */
  void setPanel(int panel) {
    // tabbedPane.setSelectedIndex(panel);
  }

  /**
   * Setter for the field <code>port</code>.
   */
  private void setPort() {
    String s = "" + port;
    s = JOptionPane.showInputDialog("Enter a port number", s);
    if (s == null) {
      return;
    }
    int temp;
    try {
      temp = Integer.parseInt(s);
    } catch (NumberFormatException e) {
      return;
    }
    port = temp;
    socketAddress = new InetSocketAddress(address, port);
  }

  /**
   * setRemoteJobAddress
   */
  private void setRemoteJobAddress() {
    if (address == null) {
      try {
        address = InetAddress.getLocalHost();
      } catch (Exception e) {
        try {
          address = InetAddress.getByName(null);
        } catch (Exception ex) {
          return;
        }
      }
    }
    String s = address.getHostAddress();
    s = JOptionPane.showInputDialog("Enter an IP Address (XXX.XXX.XXX.XXX)", s);
    if (s == null) {
      return;
    }
    InetAddress newAddress;
    InetSocketAddress newSocketAddress;
    try {
      newAddress = InetAddress.getByName(s);
      newSocketAddress = new InetSocketAddress(newAddress, port);
    } catch (Exception e) {
      return;
    }
    address = newAddress;
    socketAddress = newSocketAddress;
  }

  /**
   * showGlobalAxes
   *
   * @param evt a {@link java.awt.event.ActionEvent} object.
   */
  private void showGlobalAxes(ActionEvent evt) {
    JCheckBoxMenuItem showAxesCheckBox = (JCheckBoxMenuItem) evt.getSource();
    graphicsCanvas.setAxisShowing(showAxesCheckBox.isSelected());
  }

  /**
   * showToolBar
   *
   * @param evt a {@link java.awt.event.ActionEvent} object.
   */
  private void showToolBar(ActionEvent evt) {
    JCheckBoxMenuItem toolBarCheckBox = (JCheckBoxMenuItem) evt.getSource();
    if (toolBarCheckBox.isSelected()) {
      add(mainMenu.getToolBar(), BorderLayout.NORTH);
      frame.validate();
    } else {
      remove(mainMenu.getToolBar());
      frame.validate();
    }
  }

  /**
   * showTree
   *
   * @param evt a {@link java.awt.event.ActionEvent} object.
   */
  private void showTree(ActionEvent evt) {
    JCheckBoxMenuItem treeCheckBox = (JCheckBoxMenuItem) evt.getSource();
    if (treeCheckBox.isSelected()) {
      if (splitPaneDivider < frame.getWidth() * (1.0f / 4.0f)) {
        splitPaneDivider = (int) (frame.getWidth() * (1.0f / 4.0f));
      }
      splitPane.setDividerLocation(splitPaneDivider);
    } else {
      splitPaneDivider = splitPane.getDividerLocation();
      splitPane.setDividerLocation(0.0);
    }
  }

  /**
   * skip
   */
  private void skip() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    String skip = "" + trajectory.getSkip();
    skip = JOptionPane.showInputDialog("Enter the Number of Frames to Skip", skip);
    try {
      int f = Integer.parseInt(skip);
      trajectory.setSkip(f);
    } catch (NumberFormatException e) {
      //
    }
  }

  /**
   * speed
   */
  private void speed() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    String rate = "" + trajectory.getRate();
    rate = JOptionPane.showInputDialog("Enter the Frame Rate (1-100)", rate);
    try {
      int f = Integer.parseInt(rate);
      trajectory.setRate(f);
    } catch (NumberFormatException e) {
      //
    }
  }

  /**
   * stepBack
   */
  private void stepBack() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    trajectory.stop();
    trajectory.back();
  }

  /**
   * stepForward
   */
  private void stepForward() {
    Trajectory trajectory = getTrajectory();
    if (trajectory == null) {
      return;
    }
    trajectory.stop();
    trajectory.forward();
  }

  /**
   * Enumerates the exit status codes FFX may terminate with; 0 will indicate normal execution, 1-99
   * will indicate fatal errors, 100-199 non-fatal errors, and 200-254 other exit statuses.
   * Presently, only 0, 1, 3, 100, and 200 have been defined for FFX.
   *
   * <p>When adding to this enumeration, avoid the ranges 2, 64-78, 126-128, 130, 137, and 255 or
   * greater (see https://tldp.org/LDP/abs/html/exitcodes.html and the C/C++ standard
   * /usr/include/sysexits.h).
   */
  enum ExitStatus {
    // Normal termination.
    NORMAL(0),
    // Indicates some uncaught Exception, Error, or Throwable. As of now,
    // this enum value is unused, and we rely on the JVM automatically exiting
    // with a system code of 1 under these circumstances.
    EXCEPTION(1),
    // A call to logger.severe() resulted in program termination.
    SEVERE(3),
    // Algorithm did not complete properly (a minimization had a bad
    // interpolation, etc).
    ALGORITHM_FAILURE(100),
    // Some issue not listed here.
    OTHER(200);

    // This gets sent to System.exit().
    private final int exitCode;

    ExitStatus(int exitCode) {
      this.exitCode = exitCode;
    }

    /**
     * Gets the exit code associated with this exit status.
     *
     * @return JVM exit code.
     */
    int getExitCode() {
      return exitCode;
    }
  }
}
