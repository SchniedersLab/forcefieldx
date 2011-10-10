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
package ffx;

import java.awt.GraphicsEnvironment;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.net.InetAddress;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import static java.lang.String.format;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.UIManager;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.SystemUtils;
import org.apache.commons.lang.builder.ToStringBuilder;
import org.apache.commons.lang.time.StopWatch;

import edu.rit.pj.Comm;

import ffx.ui.LogHandler;
import ffx.ui.MainPanel;
import ffx.ui.macosx.OSXAdapter;
import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Date;

/**
 * The Main class is the entry point to the graphical user interface version of
 * Force Field X.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 * @version $Id: $
 */
public class Main extends JFrame {

    /**
     * Process any "-D" command line flags.
     */
    private static String[] processProperties(String args[]) {
        List newArgs = new ArrayList<String>();
        for (int i = 0; i < args.length; i++) {
            String arg = args[i].trim();
            if (arg.startsWith("-D")) {
                // Remove -D from the front of String.
                arg = arg.substring(2);
                // Split at the equals if it exists.
                if (arg.contains("=")) {
                    int equalsPosition = arg.indexOf("=");
                    String key = arg.substring(0, equalsPosition);
                    String value = arg.substring(equalsPosition + 1);
                    System.setProperty(key, value);
                } else {
                    if (arg.length() > 0) {
                        System.setProperty(arg, "");
                    }
                }
            } else {
                newArgs.add(arg);
            }
        }
        args = new String[newArgs.size()];
        newArgs.toArray(args);

        return args;
    }

    /**
     * Print out credits.
     */
    private static void header(String args[]) {
        StringBuilder sb = new StringBuilder();
        sb.append(MainPanel.border).append("\n");
        sb.append(MainPanel.title).append("\n");
        sb.append(MainPanel.aboutString).append("\n");
        sb.append(MainPanel.border);
        sb.append("\n ").append(new Date());
        sb.append("\n Command line arguments:\n ");
        sb.append(Arrays.toString(args));
        logger.info(sb.toString());
    }

    /**
     * Print out help for the command line version of Force Field X.
     */
    private static void commandLineInteraceHelp() {
        logger.info(" usage: ffxc [-D<property=value>] <command> [-options] <PDB|XYZ>");
        logger.info("\n where commands include:\n");
        ClassLoader classLoader = ClassLoader.getSystemClassLoader();
        classLoader.getResource("List all scripts");
        logger.info("\n For help on a spcific command: ffxc command -h\n");
        System.exit(0);
    }

    /**
     * Determine the host name, process ID, and FFX base directory.
     */
    private static void environment() {
        try {
            InetAddress addr = InetAddress.getLocalHost();
            hostName = addr.getHostName();
        } catch (UnknownHostException e) {
            // Do nothing.
        }
        procID = Integer.parseInt(System.getProperty("app.pid"));
        ffxDirectory = new File(System.getProperty("basedir"));
        try {
            logger.fine(String.format(" Force Field X directory is %s", ffxDirectory.getCanonicalPath()));
        } catch (Exception e) {
            // Do Nothing.
        }

    }

    /**
     * Start up the Parallel Java communication layer.
     */
    private static void startParallelJava(String args[]) {

        /**
         * Capture Parallel Java Output.
         */
        PrintStream origOut = System.out;
        PrintStream origErr = System.err;

        OutputStream redirectOut = new ByteArrayOutputStream();
        PrintStream out = new PrintStream(redirectOut);
        System.setOut(out);

        OutputStream redirectErr = new ByteArrayOutputStream();
        PrintStream err = new PrintStream(redirectErr);
        System.setErr(err);

        try {
            Comm.init(args);
            Comm world = Comm.world();
            rank = world.rank();
            processes = world.size();
        } catch (Exception e) {
            String message = " Exception starting up Parallel Java communication layer.";
            logger.log(Level.WARNING, message, e.toString());
        }

        /**
         * Revert standard output/error.
         */
        System.setOut(origOut);
        System.setErr(origErr);

        
        /**
         * If more than one process is requested, log any Parallel Java output.
         */
        String processesRequested = System.getProperty("pj.nn", "1");
        if (!processesRequested.equals("1")) {
            out.flush();
            err.flush();
            logger.info(redirectOut.toString());
            logger.info(redirectErr.toString());
        }

        /**
         * Close the temporary streams.
         */
        out.close();
        err.close();
    }

    /**
     * Start the Force Field X command line interface.
     * @param commandLineFile
     * @param argList
     */
    private static void startCommandLineInterface(File commandLineFile, List<String> argList) {
        if (processes == 1) {
            logger.info(String.format(" Process ID %d on %s.", procID, hostName));
            logger.info(String.format(" Starting up the command line interface."));
        } else {
            logger.info(String.format(" Starting up process ID %d (rank %d) on %s.\n", procID, rank, hostName));
        }
        HeadlessMain m = new HeadlessMain(commandLineFile, argList, logHandler);
    }

    /**
     * Start the Force Field X graphical user interface.
     * @param commandLineFile
     * @param argList
     */
    private static void startGraphicalUserInterface(File commandLineFile, List<String> argList) {
        logger.info(String.format(" Process ID %d on %s.", procID, hostName));
        logger.info(String.format(" Starting up the graphical user interface."));

        /**
         * Some Mac OS X specific features that help FFX look native.
         * These need to be set before the MainPanel is created.
         */
        if (SystemUtils.IS_OS_MAC_OSX) {
            OSXAdapter.setOSXProperties();
        }
        /**
         * Set some Swing Constants.
         */
        UIManager.put("swing.boldMetal", Boolean.FALSE);
        setDefaultLookAndFeelDecorated(false);
        /**
         * Initialize the Main frame and Force Field X MainPanel.
         */
        Main m = new Main(commandLineFile, argList);
    }

    /**
     * Create an instance of Force Field X
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.lang.Exception if any.
     */
    public static void main(String[] args) throws Exception {
        /**
         * Print out the header.
         */
        header(args);

        /**
         * Print out help for the command line interface.
         */
        if (GraphicsEnvironment.isHeadless() && args.length < 2) {
            commandLineInteraceHelp();
        }

        /**
         * Process any "-D" command line flags.
         */
        args = processProperties(args);

        /**
         * Start up the Parallel Java communication layer.
         */
        startParallelJava(args);

        /**
         * Determine host name and process ID.
         */
        environment();

        File commandLineFile = null;
        int nArgs = args.length;
        if (nArgs > 0) {
            commandLineFile = new File(args[0]);
            // Resolve a relavtive path
            if (commandLineFile.exists()) {
                commandLineFile = new File(FilenameUtils.normalize(
                        commandLineFile.getAbsolutePath()));
            }
        }

        /**
         * Convert the args to a List<String>.
         */
        List<String> argList = new ArrayList<String>(nArgs);
        if (nArgs > 1) {
            for (int i = 1; i < nArgs; i++) {
                argList.add(args[i]);
            }
        }

        /**
         * Start up the GUI or CLI version of Force Field X.
         */
        if (!GraphicsEnvironment.isHeadless()) {
            startGraphicalUserInterface(commandLineFile, argList);
        } else {
            startCommandLineInterface(commandLineFile, argList);
        }
    }

    /**
     * Main does some window initializations.
     *
     * @param commandLineFile a {@link java.io.File} object.
     * @param argList a {@link java.util.List} object.
     */
    public Main(File commandLineFile, List<String> argList) {
        super("Force Field X");
        // Start the clock.
        stopWatch.start();
        setVisible(false);

        // Create the MainPanel and MainMenu, then add them to the JFrame
        java.awt.Toolkit.getDefaultToolkit().setDynamicLayout(true);
        mainPanel = new MainPanel(this);
        logHandler.setMainPanel(mainPanel);
        add(mainPanel);
        mainPanel.initialize();
        setJMenuBar(mainPanel.getMainMenu());
        // Set the Title and Icon
        setTitle("Force Field X");
        URL iconURL = getClass().getClassLoader().getResource(
                "ffx/ui/icons/icon64.png");
        ImageIcon icon = new ImageIcon(iconURL);
        setIconImage(icon.getImage());
        addWindowListener(new WindowAdapter() {

            @Override
            public void windowClosing(WindowEvent e) {
                if (mainPanel != null) {
                    mainPanel.exit();
                }
                System.exit(0);
            }
        });
        // This is a hack to get GraphicsCanvis to initialize on some
        // platform/Java3D combinations.
        mainPanel.setPanel(MainPanel.KEYWORDS);
        setVisible(true);
        mainPanel.setPanel(MainPanel.GRAPHICS);
        // Mac OS X specific features that help Force Field X look native
        // on Macs. This needs to be done after the MainPanel is created.
        if (SystemUtils.IS_OS_MAC_OSX) {
            OSXAdapter.macOSXRegistration(mainPanel);
        }
        // Finally, open the supplied file if necessary.
        if (commandLineFile != null && !commandLineFile.exists()) {
            /**
             * See if the commandLineFile is an embedded script.
             */
            String name = commandLineFile.getName();
            name = name.replace('.', File.separatorChar);
            ClassLoader loader = getClass().getClassLoader();
            URL embeddedScript = loader.getResource("ffx/scripts/" + name + ".ffx");
            if (embeddedScript == null) {
                embeddedScript = loader.getResource("ffx/scripts/" + name + ".groovy");
            }
            if (embeddedScript != null) {
                try {
                    commandLineFile = new File(
                            FFXClassLoader.copyInputStreamToTmpFile(
                            embeddedScript.openStream(), commandLineFile.getName(), ".ffx"));
                } catch (Exception e) {
                    logger.warning("Exception extracting embedded script "
                                   + embeddedScript.toString() + "\n" + e.toString());
                }
            }
        }
        if (commandLineFile != null) {
            if (commandLineFile.exists()) {
                mainPanel.getModelingShell().setArgList(argList);
                mainPanel.open(commandLineFile, null);
            } else {
                logger.warning(format("%s was not found.", commandLineFile.toString()));
            }
        }
        if (System.getProperty("ffx.timer") != null) {
            StringBuilder sb = new StringBuilder();
            sb.append(format("\n Start-up Time (msec): %s.", stopWatch.getTime()));
            Runtime runtime = Runtime.getRuntime();
            runtime.runFinalization();
            runtime.gc();
            long occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
            long KB = 1024;
            sb.append(format("\n In-Use Memory   (Kb): %d", occupiedMemory / KB));
            sb.append(format("\n Free Memory     (Kb): %d", runtime.freeMemory() / KB));
            sb.append(format("\n Total Memory    (Kb): %d", runtime.totalMemory() / KB));
            logger.info(sb.toString());
        }
    }

    /**
     * {@inheritDoc}
     *
     * Commons.Lang Style toString.
     */
    @Override
    public String toString() {
        ToStringBuilder toStringBuilder = new ToStringBuilder(this).append(
                "Up Time: " + stopWatch).append("Logger: " + logger.getName());
        return toStringBuilder.toString();
    }
    private static final Logger logger = Logger.getLogger(Main.class.getName());
    private static final Level level;
    private static final LogHandler logHandler;
    /** Constant <code>stopWatch</code> */
    public static StopWatch stopWatch = new StopWatch();
    /**
     * This is the main application wrapper.
     */
    public MainPanel mainPanel;
    /**
     * Rank of this process for a multi-process Parallel Java FFX job.
     */
    private static int rank = 0;
    /**
     * Number of processes for a multi-process Parallel Java FFX job.
     */
    private static int processes = 1;
    /**
     * Name of the machine FFX is running on.
     */
    private static String hostName = null;
    /**
     * Force Field X base directory.
     */
    private static File ffxDirectory = null;
    /**
     * Force Field X process ID.
     */
    private static int procID = -1;

    /**
     * Replace the default console handler with our custom FFX handler.
     */
    static {
        try {
            Logger defaultLogger = LogManager.getLogManager().getLogger("");
            Handler defaultHandlers[] = defaultLogger.getHandlers();
            for (Handler h : defaultHandlers) {
                defaultLogger.removeHandler(h);
            }
        } catch (Exception e) {
            System.err.println(e.toString());
        }

        String logLevel = System.getProperty("ffx.log", "info");
        Level tempLevel;
        try {
            tempLevel = Level.parse(logLevel.toUpperCase());
        } catch (Exception e) {
            tempLevel = Level.INFO;
        }

        level = tempLevel;
        logHandler = new LogHandler();
        logHandler.setLevel(level);
        Logger ffxLogger = Logger.getLogger("ffx");
        ffxLogger.addHandler(logHandler);
        ffxLogger.setLevel(level);
    }
}
