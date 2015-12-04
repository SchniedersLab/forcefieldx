/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx;

import java.awt.GraphicsEnvironment;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.net.InetAddress;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
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
import org.apache.commons.lang3.SystemUtils;
import org.apache.commons.lang3.builder.ToStringBuilder;
import org.apache.commons.lang3.time.StopWatch;

import edu.rit.pj.Comm;

import ffx.ui.LogHandler;
import ffx.ui.MainPanel;
import ffx.ui.OSXAdapter;

/**
 * The Main class is the entry point to the graphical user interface version of
 * Force Field X.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class Main extends JFrame {

    private static final Logger logger = Logger.getLogger(Main.class.getName());
    private static Level level;
    private static LogHandler logHandler;

    /**
     * Process any "-D" command line flags.
     */
    private static String[] processProperties(String args[]) {
        List<String> newArgs = new ArrayList<>();
        for (String arg : args) {
            arg = arg.trim();
            if (arg.startsWith("-D")) {
                // Remove -D from the front of String.
                arg = arg.substring(2);
                // Split at the first equals if it exists.
                if (arg.contains("=")) {
                    int equalsPosition = arg.indexOf("=");
                    String key = arg.substring(0, equalsPosition);
                    String value = arg.substring(equalsPosition + 1);
                    // Set the system property.
                    System.setProperty(key, value);
                } else {
                    if (arg.length() > 0) {
                        System.setProperty(arg, "");
                    }
                }
            } else {
                // Collect non "-D" arguments.
                newArgs.add(arg);
            }
        }
        // Return the remaining arguments.
        args = new String[newArgs.size()];
        newArgs.toArray(args);
        return args;
    }

    /**
     * Replace the default console handler with our custom FFX handler.
     */
    private static void startLogging() {
        // Remove all log handlers from the default logger.
        try {
            Logger defaultLogger = LogManager.getLogManager().getLogger("");
            Handler defaultHandlers[] = defaultLogger.getHandlers();
            for (Handler h : defaultHandlers) {
                defaultLogger.removeHandler(h);
            }
        } catch (Exception e) {
            System.err.println(e.toString());
        }

        // Retrieve the log level from the ffx.log system property.
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

    /**
     * Print out a promo.
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
    private static void commandLineInterfaceHelp() {
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

        String procString = System.getProperty("app.pid");
        if (procString != null) {
            procID = Integer.parseInt(procString);
        } else {
            procID = 0;
        }

        String dirString = System.getProperty("basedir");
        if (dirString != null) {
            ffxDirectory = new File(dirString);
        } else {
            ffxDirectory = new File(".");
        }

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
        try {
            Comm.init(args);
            Comm world = Comm.world();
            rank = world.rank();
            processes = world.size();
        } catch (Exception e) {
            String message = String.format(" Exception starting up the Parallel Java communication layer.");
            logger.log(Level.WARNING, message, e.toString());
        }
    }

    /**
     * Start the Force Field X command line interface.
     *
     * @param commandLineFile
     * @param argList
     */
    private static void startCommandLineInterface(File commandLineFile, List<String> argList) {
        if (processes == 1) {
            logger.info(String.format(" Process ID %d on %s.", procID, hostName));
            logger.info(String.format(" Starting up the command line interface."));
        } else {
            logger.info(String.format(" Starting up rank %d on %s.\n", rank, hostName));
        }
        HeadlessMain m = new HeadlessMain(commandLineFile, argList, logHandler);
    }

    /**
     * Start the Force Field X graphical user interface.
     *
     * @param commandLineFile
     * @param argList
     */
    private static void startGraphicalUserInterface(File commandLineFile, List<String> argList) {
        logger.info(String.format(" Process ID %d on %s.", procID, hostName));
        logger.info(String.format(" Starting up the graphical user interface."));

        /**
         * Some Mac OS X specific features that help FFX look native. These need
         * to be set before the MainPanel is created.
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
         * Process any "-D" command line flags.
         */
        args = processProperties(args);

        /**
         * Configure our logging.
         */
        startLogging();

        /**
         * Print out help for the command line interface.
         */
        if (GraphicsEnvironment.isHeadless() && args.length < 2) {
            commandLineInterfaceHelp();
        }

        /**
         * Determine host name and process ID.
         */
        environment();

        /**
         * Start up the Parallel Java communication layer.
         */
        startParallelJava(args);

        /**
         * Run the pKa input GUI if requested. Halts execution until GUI exits.
         */
        /**
         * if (System.getProperty("pKaCalc") != null) { if
         * (System.getProperty("pKaCalc").equals("true")) { ffx.pka.pKaRun
         * runnable = new ffx.pka.pKaRun(); Thread t = new Thread(runnable,"pKa
         * Thread"); t.start(); t.join(); final int NUM_PKA_ARGS = 25; String[]
         * newArgs = new String[NUM_PKA_ARGS]; int currentArg = 0; for (int i=0;
         * i < newArgs.length; i++) { newArgs[currentArg] = runnable.getArg(i);
         * if (runnable.getArg(i) == null) { String temp = runnable.getArg(i -
         * 1); if (temp.startsWith("-s") || temp.startsWith("-f")) {
         * currentArg--; } } else { currentArg++; } } args = newArgs; } }
         */
        // Print the header.
        // Moved this here so I could see the args being supplied by pKaRun.
        header(args);

        /**
         * Parse the specified command or structure file.
         */
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
        List<String> argList = new ArrayList<>(nArgs);
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
            name = name.replace('.', '/');
            String pathName = "ffx/scripts/" + name;
            ClassLoader loader = getClass().getClassLoader();
            URL embeddedScript = loader.getResource(pathName + ".ffx");
            if (embeddedScript == null) {
                embeddedScript = loader.getResource(pathName + ".groovy");
            }
            if (embeddedScript != null) {
                try {
                    commandLineFile = new File(
                            FFXClassLoader.copyInputStreamToTmpFile(
                                    embeddedScript.openStream(), commandLineFile.getName(), ".ffx"));
                } catch (Exception e) {
                    logger.info(String.format(" The embedded script %s could not be extracted.", embeddedScript));
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

        if (logger.isLoggable(Level.FINE)) {
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
            logger.fine(sb.toString());
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
    /**
     * Constant <code>stopWatch</code>
     */
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
}
