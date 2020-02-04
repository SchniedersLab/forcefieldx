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
package ffx.ui;

import javax.help.HelpSet;
import javax.help.JHelp;
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
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.MalformedURLException;
import java.net.URL;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.Preferences;
import static java.lang.System.arraycopy;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.SystemUtils;
import org.biojava.nbio.structure.Structure;
import org.jogamp.java3d.GraphicsConfigTemplate3D;
import org.jogamp.java3d.Transform3D;
import org.jogamp.java3d.TransformGroup;
import org.jogamp.vecmath.Vector3d;

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
import ffx.potential.parsers.BiojavaDataFilter;
import ffx.potential.parsers.BioJavaFilter;
import ffx.potential.parsers.ConversionFilter;
import ffx.potential.parsers.DYNFileFilter;
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
import ffx.potential.utils.PotentialsDataConverter;
import ffx.ui.properties.FFXLocale;
import ffx.utilities.Keyword;
import ffx.utilities.Resources;
import static ffx.utilities.FileUtils.copyInputStreamToTmpFile;
import static ffx.utilities.StringUtils.pdbForID;

/**
 * The MainPanel class is the main container for Force Field X, handles file
 * input/output and is used to pass references among the various sub-Panels.
 *
 * @author Michael J. Schnieders
 */
public final class MainPanel extends JPanel implements ActionListener,
        ChangeListener {

    private static final Logger logger = Logger.getLogger(MainPanel.class.getName());
    /**
     * Constant <code>GRAPHICS=0</code>
     */
    private static final int GRAPHICS = 0;
    /**
     * Constant <code>KEYWORDS=1</code>
     */
    static final int KEYWORDS = 1;
    /**
     * Constant <code>MODELING=2</code>
     */
    static final int MODELING = 2;
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
     * Constant <code>dynFileFilter</code>
     */
    public static final DYNFileFilter dynFileFilter = new DYNFileFilter();
    /**
     * Constant <code>biojavaDataFilter</code>
     */
    private static final BiojavaDataFilter biojavaDataFilter = new BiojavaDataFilter();
    /**
     * Constant <code>indFileFilter</code>
     */
    private static final InducedFileFilter indFileFilter = new InducedFileFilter();
    /**
     * Constant <code>forceFieldFileFilter</code>
     */
    static final ForceFieldFileFilter forceFieldFileFilter = new ForceFieldFileFilter();
    /**
     * Constant <code>pdbFileFilter</code>
     */
    private static final PDBFileFilter pdbFileFilter = new PDBFileFilter();
    /**
     * Constant <code>keyFileFilter</code>
     */
    static final KeyFileFilter keyFileFilter = new KeyFileFilter();
    /**
     * Constant <code>ffxFileFilter</code>
     */
    private static final FFXFileFilter ffxFileFilter = new FFXFileFilter();
    /**
     * Number of File Opener Threads.
     */
    private int fileOpenerThreads = -1;

    static {
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
    private JFrame frame;
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
    // private LogPanel logPanel;
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
     * Filter to for system conversion.
     */
    private ConversionFilter activeConvFilter = null;
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
     * <p>
     * Constructor for MainPanel.</p>
     */
    public MainPanel() {
        frame = null;
    }

    /**
     * <p>
     * getPWD</p>
     *
     * @return a {@link java.io.File} object.
     */
    static File getPWD() {
        if (pwd == null) {
            pwd = new File(System.getProperty("user.dir",
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
     * <p>
     * help</p>
     */
    private void help() {
        String helpHS = "ffx/help/jhelpset.hs";
        ClassLoader cl = getClass().getClassLoader();
        HelpSet hs;
        try {
            URL hsURL = HelpSet.findHelpSet(cl, helpHS);
            hs = new HelpSet(null, hsURL);
        } catch (Exception e) {
            logger.warning("HelpSet not found: " + e);
            return;
        }
        JHelp jhelp = new JHelp(hs);
        JFrame helpFrame = new JFrame();
        JFrame.setDefaultLookAndFeelDecorated(true);
        helpFrame.add(jhelp);
        helpFrame.setTitle("Force Field X Help");
        helpFrame.setSize(jhelp.getPreferredSize());
        helpFrame.setVisible(true);
        jhelp.setCurrentID("ForceFieldXBook");
        helpFrame.toFront();
    }

    /**
     * <p>
     * about</p>
     */
    public void about() {
        if (aboutDialog == null) {
            aboutDialog = new JDialog(frame, "About... ", true);
            URL ffxURL = getClass().getClassLoader().getResource(
                    "ffx/ui/icons/splash.png");
            ImageIcon logoIcon = new ImageIcon(ffxURL);
            JLabel logoLabel = new JLabel(logoIcon);
            logoLabel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
            Container contentpane = aboutDialog.getContentPane();
            contentpane.setLayout(new BorderLayout());
            initAbout();
            contentpane.add(aboutTextArea, BorderLayout.SOUTH);
            contentpane.add(logoLabel, BorderLayout.CENTER);
            aboutDialog.pack();
            Dimension dim = getToolkit().getScreenSize();
            Dimension ddim = aboutDialog.getSize();
            aboutDialog.setLocation((dim.width - ddim.width) / 2,
                    (dim.height - ddim.height) / 2);
            aboutDialog.setResizable(false);
        }
        aboutDialog.setVisible(true);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Handle most File, Selection, Trajectory, Simulation, Window and Help Menu
     * Commands This should probably be partitioned between a few different
     * handlers
     */
    @Override
    public void actionPerformed(ActionEvent evt) {
        String arg = evt.getActionCommand();
        if (logger.isLoggable(Level.FINEST)) {
            logger.finest(" Action: " + arg);
        }
        // File Commands
        if (arg.equals("Open")) {
            open();
        } else if (arg.equals("DownloadFromPDB")) {
            openFromPDB();
        } else if (arg.equals("SaveAs")) {
            saveAsXYZ(null);
        } else if (arg.equals("Close")) {
            close();
        } else if (arg.equals("CloseAll")) {
            closeAll();
        } else if (arg.equals("ChooseKeyFile")) {
            chooseKey();
        } else if (arg.equals("ChooseLogFile")) {
            chooseLog();
        } else if (arg.equals("LoadInducedData")) {
            openInduced();
            // Selection Commands
        } else if (arg.equals("SelectAll")) {
            selectAll();
        } else if (arg.equals("MergeSelections")) {
            merge();
        } else if (arg.equals("HighlightSelections")) {
            highlightSelections(evt);
            // Trajectory
        } else if (arg.equals("Play")) {
            play();
        } else if (arg.equals("Stop")) {
            stop();
        } else if (arg.equals("StepForward")) {
            stepForward();
        } else if (arg.equals("StepBack")) {
            stepBack();
        } else if (arg.equals("Reset")) {
            reset();
        } else if (arg.equals("Oscillate")) {
            oscillate(evt);
        } else if (arg.equals("Frame")) {
            frame();
        } else if (arg.equals("Speed")) {
            speed();
        } else if (arg.equals("Skip")) {
            skip();
            // Simulation
        } else if (arg.equals("ConnectToLocalJob")) {
            connectToTINKER(null, null);
        } else if (arg.equals("ConnectToRemoteJob")) {
            connect();
        } else if (arg.equals("ReleaseJob")) {
            release();
        } else if (arg.equals("SetPort")) {
            setPort();
        } else if (arg.equals("SetRemoteJobAddress")) {
            setRemoteJobAddress();
            // Window
        } else if (arg.equals("ShowToolBar")) {
            showToolBar(evt);
        } else if (arg.equals("ShowTree")) {
            showTree(evt);
        } else if (arg.equals("ShowGlobalAxes")) {
            showGlobalAxes(evt);
        } else if (arg.equals("ResetPanes")) {
            resetPanes();
        } else if (arg.equals("ResetConsole")) {
            resetShell();
        } else if (arg.equals("OceanLookAndFeel")) {
            oceanLookAndFeel();
        } else if (arg.equals("WindowsLookAndFeel") || arg.equals("MacOSXLookAndFeel") || arg.equals("MotifLookAndFeel")) {
            platformLookAndFeel();
        } else if (arg.equals("ShrinkGraphicsWindow")) {
            resizePanes(20);
        } else if (arg.equals("ExpandGraphicsWindow")) {
            resizePanes(-20);
            // Help
        } else if (arg.equals("HelpContents")) {
            help();
        } else if (arg.equals("About")) {
            about();
            // Others
        } else if (arg.equals("GarbageCollect")) {
            Runtime.getRuntime().runFinalization();
            Runtime.getRuntime().gc();
        } else if (arg.equals("Exit")) {
            exit();
        } else {
            try {
                ClassLoader cl = MainPanel.class.getClassLoader();
                URL url = cl.getResource(arg);
                logger.info(url.toString());
                File structureFile = new File(url.getFile());
                logger.info(structureFile.toString());
                String tempFile = copyInputStreamToTmpFile(url.openStream(), "ffx", structureFile.getName(), "pdb");
                open(tempFile);
            } catch (Exception e) {
                System.err.println("MainPanel - Menu command not found: " + arg);
            }

        }
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
                //getModelingPanel().selected();
            }
        }
    }

    /**
     * Detach the active FSystem's BranchGroup from the Scene and clear that
     * FSystem's data
     *
     * @return a {@link java.lang.Thread} object.
     */
    public Thread close() {
        FFXSystem m = hierarchy.getActive();
        return close(m);
    }

    /**
     * <p>
     * closeWait</p>
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
     * <p>
     * close</p>
     *
     * @param closedModel a {@link ffx.ui.FFXSystem} object.
     * @return a {@link java.lang.Thread} object.
     */
    public Thread close(FFXSystem closedModel) {
        if (closedModel != null && closedModel.getParent() != null) {
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
     * <p>
     * connectToTINKER</p>
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
                tempAddress = new InetSocketAddress(InetAddress.getLocalHost(),
                        port);
            } catch (Exception e) {
                try {
                    tempAddress = new InetSocketAddress(InetAddress.getByName(null), port);
                } catch (Exception ex) {
                    System.err.println("Could not determine Local Host: " + ex);
                    return;
                }
            }
            simulation = new SimulationLoader(system, modelingThread, this,
                    tempAddress);
            if (modelingThread != null) {
                modelingThread.start();
            }
            simulation.connect();
            mainMenu.setConnect(false);
            setPanel(GRAPHICS);
        }
    }

    /**
     * <p>
     * createKeyFile</p>
     *
     * @param system a {@link ffx.ui.FFXSystem} object.
     * @return a boolean.
     */
    boolean createKeyFile(FFXSystem system) {
        String message = new String("Please select a parameter file " + "and a TINKER Key file will be created.");
        String params = (String) JOptionPane.showInputDialog(this, message,
                "Parameter File", JOptionPane.QUESTION_MESSAGE, null,
                keywordPanel.getParamFiles(), null);
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
                    message = "Could not create a Key file because " + pwd.getAbsolutePath() + " is not writable";
                    JOptionPane.showMessageDialog(this, message);
                }
            }
        }
        return false;
    }

    /**
     * <p>
     * exit with current exit code (default: 0 (ExitStatus.NORMAL))</p>
     */
    public void exit() {
        exit(exitType);
    }

    /**
     * <p>
     * exit with a target ExitStatus</p>
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
     * <p>
     * frame</p>
     */
    public void frame() {
        Trajectory trajectory = getTrajectory();
        if (trajectory == null) {
            return;
        }
        String frameNumber = "" + trajectory.getFrame();
        frameNumber = JOptionPane.showInputDialog("Enter the Frame Number",
                frameNumber);
        try {
            int f = Integer.parseInt(frameNumber);
            trajectory.setFrame(f);
        } catch (NumberFormatException e) {
            //
        }
    }

    /**
     * <p>
     * Getter for the field <code>dataRoot</code>.</p>
     *
     * @return a {@link ffx.potential.bonded.MSRoot} object.
     */
    MSRoot getDataRoot() {
        return dataRoot;
    }

    /**
     * <p>
     * getFFXLocale</p>
     *
     * @return a {@link ffx.ui.properties.FFXLocale} object.
     */
    FFXLocale getFFXLocale() {
        return locale;
    }

    /**
     * <p>
     * getGraphics3D</p>
     *
     * @return a {@link ffx.ui.GraphicsCanvas} object.
     */
    GraphicsCanvas getGraphics3D() {
        return graphicsCanvas;
    }

    /**
     * <p>
     * Getter for the field <code>hierarchy</code>.</p>
     *
     * @return a {@link ffx.ui.Hierarchy} object.
     */
    public Hierarchy getHierarchy() {
        return hierarchy;
    }

    /**
     * <p>
     * Getter for the field <code>keywordPanel</code>.</p>
     *
     * @return a {@link ffx.ui.KeywordPanel} object.
     */
    KeywordPanel getKeywordPanel() {
        return keywordPanel;
    }

    /**
     * <p>
     * Getter for the field <code>mainMenu</code>.</p>
     *
     * @return a {@link ffx.ui.MainMenu} object.
     */
    public MainMenu getMainMenu() {
        return mainMenu;
    }

    /**
     * <p>
     * Getter for the field <code>frame</code>.</p>
     *
     * @return a {@link java.awt.Frame} object.
     */
    public Frame getFrame() {
        return frame;
    }

    ModelingPanel getModelingPanel() {
        return modelingPanel;
    }

    /**
     * <p>
     * Getter for the field <code>modelingShell</code>.</p>
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
     * <p>
     * getStatusBar</p>
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
     * Return the active SystemFilter.
     *
     * @return the active SystemFilter.
     */
    public SystemFilter getFilter() {
        return activeFilter;
    }

    /**
     * <p>
     * highlightSelections</p>
     *
     * @param evt a {@link java.awt.event.ActionEvent} object.
     */
    private void highlightSelections(ActionEvent evt) {
        if (evt.getSource() instanceof JCheckBoxMenuItem) {
            JCheckBoxMenuItem jcb = (JCheckBoxMenuItem) evt.getSource();
            if (jcb.isSelected()) {
                hierarchy.setHighlighting(true);
            } else {
                hierarchy.setHighlighting(false);
            }
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

    /**
     * Constant <code>version="1.0.0-BETA"</code>
     */
    public static final String version = "1.0.0-BETA";
    /**
     * Constant <code>date="January 2020"</code>
     */
    public static final String date = "October 2019";

    private static final String commitVersion;
    private static final String commitDate;
    private static final String commitSCM; // Presently using Git.

    /*
      Attempts to initialize version, date, and SCM versioning from
      target/ffx-mvn.properties. Fallback is to use defaults from above.
     */
    static {
        String basedir = System.getProperty("basedir");
        File mvnProps = new File(basedir + "/target/ffx-mvn.properties");
        String cVers = version + "-unknown";
        String cDate = date;
        String cSCM = "";
        if (mvnProps.exists()) {
            try (BufferedReader br = new BufferedReader(new FileReader(mvnProps))) {
                String ffxVers = "1.0.0-BETA";
                String ffxvProp = "ffx.version=";

                String commitsCount = "";
                String ccProp = "git.commitsCount=";

                String timeProp = "timestamp=";

                String scmProp = "git.revision=";

                String line = br.readLine();
                while (line != null) {
                    line = line.trim();
                    if (line.startsWith(ffxvProp)) {
                        ffxVers = line.replaceFirst(ffxvProp, "");
                    }/* else if (line.startsWith(gtProp)) {
                        gitTag = line.replaceFirst(gtProp, "");
                    } */ else if (line.startsWith(ccProp)) {
                        commitsCount = line.replaceFirst(ccProp, "");
                    } else if (line.startsWith(timeProp)) {
                        String timeStr = line.replaceFirst(timeProp, "");
                        // Expected to be MM-dd-yyyy
                        String[] timeToks = timeStr.split("-");
                        try {
                            String year = timeToks[2];
                            int mon = Integer.parseInt(timeToks[0]);
                            Month month = Month.of(mon);
                            String mstr = month.toString();
                            cDate = String.format("%c%s %s", mstr.charAt(0), mstr.substring(1).toLowerCase(), year);
                        } catch (Exception ex) {
                            cDate = date;
                        }
                    } else if (line.startsWith(scmProp) && !line.contains("UNKNOWN_REVISION")) {
                        String scm = line.replaceFirst(scmProp, "");
                        cSCM = String.format("        %s %s \n", "SCM version ", scm);
                    }
                    line = br.readLine();
                }

                StringBuilder sb = new StringBuilder(ffxVers).append("-");

                if (!commitsCount.isEmpty()) {
                    sb.append(commitsCount);
                } else {
                    sb.append("unknown");
                }
                cVers = sb.toString();
            } catch (Exception ex) {
                //
            }
        }
        commitVersion = cVers;
        commitDate = cDate;
        commitSCM = cSCM;
    }

    /**
     * Constant
     */
    public static final String border
            = " ______________________________________________________________________________\n";
    /**
     * Constant
     */
    public static final String title
            = "        FORCE FIELD X - Software for Molecular Biophysics \n";
    /**
     * Constant
     */
    public static final String aboutString
            = "        Version " + commitVersion + "  " + commitDate + " \n"
            + commitSCM // Will contain its own spacing/newline, or be empty.
            + "\n        Copyright (c)  Michael J. Schnieders    2001-2020 \n"
            + "        Portions Copyright (c)  Timothy D. Fenn 2009-2020 \n"
            + "        Portions Copyright (c)  Jacob M. Litman 2015-2020 \n"
            + "\n"
            + "        All Rights Reserved \n"
            + "\n"
            + "        Force Field X is distributed under the GPL v. 3 license\n"
            + "        with linking exception and is hosted by the Schnieders Lab \n"
            + "        at The University of Iowa.\n"
            + "\n"
            + "        For publications see    http://ffx.biochem.uiowa.edu/publications.html \n"
            + "        For the GPL v.3 license http://ffx.biochem.uiowa.edu/license.html \n";

    private void initAbout() {
        aboutTextArea = new JTextArea();
        Font font = Font.decode(Font.MONOSPACED);
        aboutTextArea.setFont(font);
        aboutTextArea.setText(aboutString);
        aboutTextArea.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
        aboutTextArea.setEditable(false);
    }

    /**
     * <p>
     * initialize</p>
     */
    public void initialize() {
        if (init) {
            return;
        }
        init = true;

        String dir = System.getProperty("user.dir", FileSystemView.getFileSystemView().getDefaultDirectory().getAbsolutePath());
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
        JScrollPane scrollPane = new JScrollPane(hierarchy,
                JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
                JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        treePane.add(scrollPane, BorderLayout.CENTER);

        if (!GraphicsEnvironment.isHeadless()) {
            GraphicsConfigTemplate3D template3D = new GraphicsConfigTemplate3D();
            // template3D.setDoubleBuffer(GraphicsConfigTemplate.PREFERRED);
            GraphicsConfiguration gc = null;
            try {
                gc = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getBestConfiguration(template3D);
                // gc = GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getDefaultConfiguration();
            } catch (Exception e) {
                e.fillInStackTrace();
                e.printStackTrace();
                logger.log(Level.SEVERE,
                        " Exception encountered when trying to get the best GraphicsConfiguration", e);
            }
            try {
                graphicsCanvas = new GraphicsCanvas(gc, this);
            } catch (Exception e) {
                e.fillInStackTrace();
                e.printStackTrace();
                logger.log(Level.SEVERE,
                        " Exception encountered when trying to create the GraphicsCanvas", e);
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
     * <p>
     * isOpening</p>
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
        if (width > screenSize.width * 0.4 && width < screenSize.width * 0.8 && height > screenSize.height * 0.4 && height < screenSize.height * 0.8) {
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
     * <p>
     * merge</p>
     */
    private void merge() {
        ArrayList<MSNode> activeNodes = hierarchy.getActiveNodes();
        if (activeNodes.size() >= 2) {
            merge(activeNodes);
        }
    }

    /**
     * Merge two or more selected FSystem Nodes into one FSystem node. There are
     * a few gotchas that need to be fixed
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
            parentSystem = (FFXSystem) m.getMSNode(FFXSystem.class);
            if (parentSystem == null) {
                return;
            }
            if (!systems.contains(parentSystem)) {
                graphicsCanvas.updateSceneWait(parentSystem, false, true,
                        RendererCache.ViewModel.WIREFRAME, false, null);
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
        MergeFilter mergeFilter = new MergeFilter(system, mergedAtoms,
                mergedBonds);
        UIFileOpener fileOpener = new UIFileOpener(mergeFilter, this);
        if (fileOpenerThreads > 0) {
            fileOpener.setNThreads(fileOpenerThreads);
        }
        Thread thread = new Thread(fileOpener);
        thread.start();
    }

    /**
     * <p>
     * merge</p>
     *
     * @param nodesToMerge an array of {@link ffx.potential.bonded.MSNode}
     *                     objects.
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

    /**
     * <p>
     * oceanLookAndFeel</p>
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

    public Thread open(File file, String commandDescription) {
        UIFileOpener opener = openInit(file, commandDescription);
        openThread = new Thread(opener);
        openThread.start();
        setPanel(GRAPHICS);
        return openThread;
    }

    public Thread convert(Object data, File file, String commandDescription) {
        UIDataConverter converter = convertInit(data, file, commandDescription);
        openThread = new Thread(converter);
        openThread.start();
        setPanel(GRAPHICS);
        return openThread;
    }

    /**
     * Attempts to load the supplied file
     *
     * @param file               File to open
     * @param commandDescription Description of the command that created this
     *                           file.
     * @return a {@link java.lang.Thread} object.
     */
    private UIFileOpener openInit(File file, String commandDescription) {
        if (file == null || !file.isFile() || !file.canRead()) {
            return null;
        }
        file = new File(FilenameUtils.normalize(file.getAbsolutePath()));
        // Set the Current Working Directory based on this file.
        setCWD(file.getParentFile());

        // Get "filename" from "filename.extension".
        String name = file.getName();
        String extension = FilenameUtils.getExtension(name);

        // Run a Force Field X script.
        if (extension.equalsIgnoreCase("ffx") || extension.equalsIgnoreCase("groovy")) {
            ModelingShell shell = getModelingShell();
            shell.runFFXScript(file);
            boolean shutDown = Boolean.parseBoolean(System.getProperty("ffx.shutDown", "true"));
            if (java.awt.GraphicsEnvironment.isHeadless() && shutDown) {
                exit();
            } else {
                return null;
            }
        }

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
                logger.info(String.format(" Loaded rotamer definitions from patch %s.", patch));
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
        //return new UIFileOpener(systemFilter, this);
    }

    /**
     * Attempts to load from the supplied data structure
     *
     * @param data               Data structure to load from
     * @param file               Source file
     * @param commandDescription Description of the command that created this
     *                           file.
     * @return A thread-based UIDataConverter
     */
    private UIDataConverter convertInit(Object data, File file, String commandDescription) {
        if (data == null) {
            return null;
        }

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
                logger.info(String.format(" Loaded rotamer definitions from patch %s.", patch));
            }
        }
        newSystem.setForceField(forceField);
        ConversionFilter convFilter;

        // Decide what parser to use.
        if (biojavaDataFilter.accept(data)) {
            convFilter = new BioJavaFilter((Structure) data, newSystem, forceField, properties);
        } else {
            throw new IllegalArgumentException("Not a recognized data structure.");
        }

        setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
        activeConvFilter = convFilter;
        return new UIDataConverter(data, file, convFilter, this);
    }

    private UIFileOpener openFromUtils(List<File> files, String commandDescription) {
        UIFileOpener openFile = openInit(files, commandDescription);
        openThread = new Thread(openFile);
        openThread.start();
        setPanel(GRAPHICS);
        return openFile;
    }

    public Thread open(List<File> files, String commandDescription) {
        UIFileOpener openFile = openInit(files, commandDescription);
        openThread = new Thread(openFile);
        openThread.start();
        setPanel(GRAPHICS);
        return openThread;
    }

    /**
     * Attempts to load the supplied file
     *
     * @param files              Files to open
     * @param commandDescription Description of the command that created this
     *                           file.
     * @return a {@link java.lang.Thread} object.
     */
    private UIFileOpener openInit(List<File> files, String commandDescription) {
        if (files == null) {
            return null;
        }
        File file = new File(FilenameUtils.normalize(files.get(0).getAbsolutePath()));
        // Set the Current Working Directory based on this file.
        setCWD(file.getParentFile());

        // Get "filename" from "filename.extension".
        String name = file.getName();

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
                logger.info(String.format(" Loaded rotamer definitions from patch %s.", patch));
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
        //return new UIFileOpener(systemFilter, this);
        UIFileOpener fileOpener = new UIFileOpener(systemFilter, this);
        if (fileOpenerThreads > 0) {
            fileOpener.setNThreads(fileOpenerThreads);
        }
        return fileOpener;
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

        MolecularAssembly[] systems = activeFilter.getMolecularAssemblys();
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
     * <p>
     * openWait</p>
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

        MolecularAssembly[] systems = activeFilter.getMolecularAssemblys();
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
     * <p>
     * openWait</p>
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
     * <p>
     * openWait</p>
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

        MolecularAssembly[] systems = activeFilter.getMolecularAssemblys();
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
     * <p>
     * openWait</p>
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

    /**
     * Converts a non-Force Field X data structure into an array of FFXSystem[].
     * Presently does not yet have support for array- or list-based data
     * structures, only singular objects.
     *
     * @param data Outside data structure
     * @param file Source file
     * @return An array of FFXSystem
     */
    synchronized FFXSystem[] convertWait(Object data, File file) {
        if (file == null) {
            try {
                file = PotentialsDataConverter.getDefaultFile(data);
            } catch (FileNotFoundException | IllegalArgumentException ex) {
                logger.warning(String.format(" Exception in finding file for data structure: %s", ex.toString()));
                return null;
            }
        }
        String name = file.getName();

        Thread thread = convert(data, file, null);
        while (thread != null && thread.isAlive()) {
            try {
                wait(1);
            } catch (InterruptedException e) {
                String message = "Exception waiting for " + name + " to open.";
                logger.log(Level.WARNING, message, e);
                return null;
            }
        }

        MolecularAssembly[] systems = activeConvFilter.getMolecularAssemblys();
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
     * <p>
     * open</p>
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

    /**
     * <p>
     * open</p>
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
        String code = JOptionPane.showInputDialog(
                "Enter the PDB Identifier (4 characters)", "");
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
            String message = String.format(" Reading the local copy of the PDB file %s.", pdbFile);
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
            fromURL = new URL(fromString);
        } catch (MalformedURLException e) {
            String message = String.format(" URL incorrectly formatted %s.", fromString);
            logger.log(Level.INFO, message, e);
            return null;
        }

        // Download the URL to a local file.
        logger.info(String.format(" Downloading %s", fromString));
        try {
            File toFile = new File(FilenameUtils.getName(fromURL.getPath()));
            FileUtils.copyURLToFile(fromURL, toFile, 1000, 1000);
            logger.info(String.format(" Saved to %s\n", toFile.getPath()));
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
     * <p>
     * openOn</p>
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
     * <p>
     * oscillate</p>
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
     * <p>
     * platformLookAndFeel</p>
     */
    private void platformLookAndFeel() {
        try {
            if (SystemUtils.IS_OS_LINUX) {
                UIManager.setLookAndFeel("com.sun.java.swing.plaf.motif.MotifLookAndFeel");
            } else {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            }
            SwingUtilities.updateComponentTreeUI(SwingUtilities.getRoot(this));
        } catch (ClassNotFoundException | InstantiationException
                | IllegalAccessException | UnsupportedLookAndFeelException e) {
            logger.log(Level.WARNING, "Can''t set look and feel: {0}", e);
        }
    }

    /**
     * <p>
     * play</p>
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
     * <p>
     * reset</p>
     */
    public void reset() {
        Trajectory trajectory = getTrajectory();
        if (trajectory == null) {
            return;
        }
        trajectory.stop();
        trajectory.rewind();
    }

    /**
     * <p>
     * resetPanes</p>
     */
    private void resetPanes() {
        resizePanes(0);
    }

    /**
     * <p>
     * resetShell</p>
     */
    void resetShell() {
        if (!GraphicsEnvironment.isHeadless()) {
            modelingShell = getModelingShell();
            modelingShell.savePrefs();
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

    public void saveKeywordFile(File file) {
        keywordPanel.keySave(file);
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

    void savePDBSymMates(File file, String suffix) {
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
        String filename = FilenameUtils.removeExtension(file.getName());
        PDBFilter pdbFilter = new PDBFilter(saveFile, system, null, null);
        if (pdbFilter.writeFile(saveFile, false)) {
            // Refresh Panels with the new System name
            hierarchy.setActive(system);
            activeFilter = pdbFilter;
        } else {
            logger.log(Level.INFO, " Save failed for: {0}", system);
        }

        Crystal crystal = system.getCrystal();
        int nSymOps = crystal.spaceGroup.getNumberOfSymOps();
        logger.info(String.format(" Writing %d symmetry mates for %s", nSymOps, system.toString()));
        for (int i = 1; i < nSymOps; i++) {
            pdbFilter.setSymOp(i);
            String saveFileName = filename + suffix + "_" + i + ".pdb";
            saveFile = new File(saveFileName);
            for (int j = 1; j < 1000; j++) {
                if (!saveFile.exists()) {
                    break;
                }
                saveFile = new File(saveFileName + "_" + j);
            }

            StringBuilder symSb = new StringBuilder();
            String[] symopLines = crystal.spaceGroup.getSymOp(i).toString().split("\\r?\\n");
            int nLines = symopLines.length;
            symSb.append("REMARK 350\nREMARK 350 SYMMETRY OPERATORS");
            for (int j = 0; j < nLines; j++) {
                symSb.append("\nREMARK 350 ").append(symopLines[j]);
            }

            symopLines = crystal.spaceGroup.getSymOp(i).toXYZString().split("\\r?\\n");
            nLines = symopLines.length;
            symSb.append("\nREMARK 350\nREMARK 350 SYMMETRY OPERATORS XYZ FORM");
            for (int j = 0; j < nLines; j++) {
                symSb.append("\nREMARK 350 ").append(symopLines[j]);
            }

            if (saveFile.exists()) {
                logger.warning(String.format(" Could not successfully version file "
                        + "%s: appending to file %s", saveFileName, saveFile.getName()));
                if (!pdbFilter.writeFileWithHeader(saveFile, symSb.toString(), true)) {
                    logger.log(Level.INFO, " Save failed for: {0}", system);
                }
            } else {
                if (!pdbFilter.writeFileWithHeader(saveFile, symSb.toString(), false)) {
                    logger.log(Level.INFO, " Save failed for: {0}", system);
                }
            }
        }
    }

    /**
     * <p>
     * saveAsPDB</p>
     *
     * @param activeSystems an array of {@link ffx.potential.MolecularAssembly}
     *                      objects.
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
            PDBFilter pdbFilter = new PDBFilter(saveFile,
                    Arrays.asList(activeSystems), null, null);
            pdbFilter.writeFile(saveFile, false);
            activeFilter = pdbFilter;
        }
    }

    private static final Preferences preferences = Preferences.userNodeForPackage(MainPanel.class);

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
        if (modelingShell != null) {
            modelingShell.savePrefs();
        }
        if (graphicsCanvas != null) {
            graphicsCanvas.savePrefs();
        }
    }

    /**
     * <p>
     * selectAll</p>
     */
    private void selectAll() {
        if (dataRoot.getChildCount() == 0) {
            return;
        }
        hierarchy.selectAll();
    }

    /**
     * <p>
     * setCWD</p>
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
     * <p>
     * setPanel</p>
     *
     * @param panel a int.
     */
    void setPanel(int panel) {
        // tabbedPane.setSelectedIndex(panel);
    }

    /**
     * <p>
     * Setter for the field <code>port</code>.</p>
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
     * <p>
     * setRemoteJobAddress</p>
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
        String s = "" + address.getHostAddress();
        s = JOptionPane.showInputDialog(
                "Enter an IP Address (XXX.XXX.XXX.XXX)", s);
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
     * <p>
     * showGlobalAxes</p>
     *
     * @param evt a {@link java.awt.event.ActionEvent} object.
     */
    private void showGlobalAxes(ActionEvent evt) {
        JCheckBoxMenuItem showAxesCheckBox = (JCheckBoxMenuItem) evt.getSource();
        graphicsCanvas.setAxisShowing(showAxesCheckBox.isSelected());
    }

    /**
     * <p>
     * showToolBar</p>
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
     * <p>
     * showTree</p>
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
     * <p>
     * skip</p>
     */
    private void skip() {
        Trajectory trajectory = getTrajectory();
        if (trajectory == null) {
            return;
        }
        String skip = "" + trajectory.getSkip();
        skip = JOptionPane.showInputDialog(
                "Enter the Number of Frames to Skip", skip);
        try {
            int f = Integer.parseInt(skip);
            trajectory.setSkip(f);
        } catch (NumberFormatException e) {
            //
        }
    }

    /**
     * <p>
     * speed</p>
     */
    private void speed() {
        Trajectory trajectory = getTrajectory();
        if (trajectory == null) {
            return;
        }
        String rate = "" + trajectory.getRate();
        rate = JOptionPane.showInputDialog("Enter the Frame Rate (1-100)",
                rate);
        try {
            int f = Integer.parseInt(rate);
            trajectory.setRate(f);
        } catch (NumberFormatException e) {
            //
        }
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
     * <p>
     * stepBack</p>
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
     * <p>
     * stepForward</p>
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
     * <p>
     * stop</p>
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
     * Enumerates the exit status codes FFX may terminate with; 0 will indicate
     * normal execution, 1-99 will indicate fatal errors, 100-199 non-fatal
     * errors, and 200-254 other exit statuses. Presently, only 0, 1, 3, 100,
     * and 200 have been defined for FFX.
     * <p>
     * When adding to this enumeration, avoid the ranges 2, 64-78, 126-128, 130,
     * 137, and 255 or greater (see http://tldp.org/LDP/abs/html/exitcodes.html
     * and the C/C++ standard /usr/include/sysexits.h).
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
