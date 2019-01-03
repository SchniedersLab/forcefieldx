/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import java.awt.GraphicsEnvironment;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Field;
import java.net.InetAddress;
import java.net.URL;
import java.net.URLDecoder;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.Logger;
import java.util.zip.ZipEntry;
import static java.lang.String.format;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.SystemUtils;
import org.apache.commons.lang3.builder.ToStringBuilder;
import org.apache.commons.lang3.time.StopWatch;

import edu.rit.pj.Comm;

import ffx.ui.LogHandler;
import ffx.ui.MainPanel;
import ffx.ui.OSXAdapter;

import sun.misc.Unsafe;

/**
 * The Main class is the entry point to the graphical user interface version of
 * Force Field X.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public final class Main extends JFrame {

    private static final Logger logger = Logger.getLogger(Main.class.getName());
    /**
     * Constant <code>stopWatch</code>
     */
    public static StopWatch stopWatch = new StopWatch();
    /**
     * This is the main application wrapper.
     */
    public static MainPanel mainPanel;
    /**
     * An Adapter class to integrate with OSX.
     */
    public static OSXAdapter osxAdapter;
    private static Level level;
    private static LogHandler logHandler;
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
     * Main does some window initializations.
     *
     * @param commandLineFile a {@link java.io.File} object.
     * @param argList         a {@link java.util.List} object.
     */
    public Main(File commandLineFile, List<String> argList) {
        super("Force Field X");
        // Start the clock.
        stopWatch.start();

        // Init the GUI.
        try {
            SwingUtilities.invokeAndWait(initGUI);
        } catch (Exception e) {
            e.printStackTrace();
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
                    commandLineFile = new File(copyInputStreamToTmpFile(
                                    embeddedScript.openStream(), commandLineFile.getName(), ".groovy"));
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
     * This Runnable is used to init the GUI using SwingUtilities.
     */
    final private Runnable initGUI = new Runnable() {
        /**
         * Create the MainPanel and MainMenu, then add them to the JFrame.
         */
        public void run() {

            Toolkit.getDefaultToolkit().setDynamicLayout(true);

            // Set the Title and Icon
            Main.this.setTitle("Force Field X");
            URL iconURL = getClass().getClassLoader().getResource(
                    "ffx/ui/icons/icon64.png");
            ImageIcon icon = new ImageIcon(iconURL);
            Main.this.setIconImage(icon.getImage());
            Main.this.addWindowListener(new WindowAdapter() {
                @Override
                public void windowClosing(WindowEvent e) {
                    if (mainPanel != null) {
                        mainPanel.exit();
                    }
                    System.exit(0);
                }
            });

            mainPanel = new MainPanel(Main.this);
            mainPanel.setVisible(false);
            logHandler.setMainPanel(mainPanel);
            Main.this.add(mainPanel);
            mainPanel.initialize();
            Main.this.setVisible(true);
            mainPanel.setVisible(true);
            Main.this.setJMenuBar(mainPanel.getMainMenu());

            // This is a hack to get GraphicsCanvis to initialize on some platform/Java3D combinations.
            // mainPanel.setPanel(MainPanel.KEYWORDS);
            // mainPanel.setVisible(true);
            // mainPanel.setPanel(MainPanel.GRAPHICS);

            // Mac OS X specific features that help Force Field X look native
            // on Macs. This needs to be done after the MainPanel is created.
            if (SystemUtils.IS_OS_MAC_OSX) {
                osxAdapter = new OSXAdapter(mainPanel);
            }

            SwingUtilities.updateComponentTreeUI(SwingUtilities.getRoot(Main.this));
        }
    };

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

        /**
         * This removes logger warnings about Illegal Access.
         */
        try {
            Field theUnsafe = Unsafe.class.getDeclaredField("theUnsafe");
            theUnsafe.setAccessible(true);
            Unsafe u = (Unsafe) theUnsafe.get(null);

            Class cls = Class.forName("jdk.internal.module.IllegalAccessLogger");
            Field logger = cls.getDeclaredField("logger");
            u.putObjectVolatile(cls, u.staticFieldOffset(logger), null);
        } catch (Exception e) {
            // ignore
        }
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
        // Print out command line arguments if the array is not null.
        if (args != null && args.length > 0) {
            sb.append("\n Command line arguments:\n ");
            sb.append(Arrays.toString(args));
        }
        logger.info(sb.toString());
    }

    /**
     * Print out help for the command line version of Force Field X.
     */
    private static void commandLineInterfaceHelp(boolean listTestScripts) {
        logger.info("\n usage: ffxc [-D<property=value>] <command> [-options] <PDB|XYZ>");
        logger.info("\n where commands include:\n");
        if (listTestScripts) {
            try {
                Main.listGroovyScripts(false, true);
            } catch (IOException e) {
                logger.log(Level.INFO, " Exception listing scripts.", e);
            }
            logger.info("\n For help on an experimental or test command use:  ffxc <command> -h\n");
        } else {
            try {
                Main.listGroovyScripts(true, false);
            } catch (IOException e) {
                logger.log(Level.INFO, " Exception listing scripts.", e);
            }
            logger.info("\n To list experimental & test scripts: ffxc --test");
            logger.info(" For help on a specific command use:  ffxc <command> -h\n");
        }
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
            logger.info(String.format(" Starting up the command line interface.\n"));
        } else {
            logger.info(String.format(" Starting up rank %d on %s.\n", rank, hostName));
        }
        HeadlessMain m = new HeadlessMain(commandLineFile, argList, logHandler);
        mainPanel = m.mainPanel;
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
         * Set some Swing Constants.
         */
        UIManager.put("swing.boldMetal", Boolean.FALSE);
        setDefaultLookAndFeelDecorated(false);

        /**
         * Some Mac OS X specific features that help FFX look native. These need
         * to be set before the MainPanel is created.
         */
        if (SystemUtils.IS_OS_MAC_OSX) {
            OSXAdapter.setOSXProperties();
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception e) {
                //
            }
        }

        /**
         * Initialize the Main frame and Force Field X MainPanel.
         */
        Main m = new Main(commandLineFile, argList);
    }

    /**
     * Create an instance of Force Field X
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String[] args) {
        /**
         * Process any "-D" command line flags.
         */
        args = processProperties(args);

        /**
         * Configure our logging.
         */
        startLogging();

        /**
         * Determine host name and process ID.
         */
        environment();

        /**
         * Start up the Parallel Java communication layer.
         */
        startParallelJava(args);

        // Print the header.
        header(args);

        /**
         * Print out help for the command line interface.
         */
        if (GraphicsEnvironment.isHeadless() && args.length < 2) {
            if (args.length == 1 && args[0].toUpperCase().contains("TEST")) {
                commandLineInterfaceHelp(true);
            }
            commandLineInterfaceHelp(false);
        }

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
     * {@inheritDoc}
     * <p>
     * Commons.Lang Style toString.
     */
    @Override
    public String toString() {
        ToStringBuilder toStringBuilder = new ToStringBuilder(this).append(
                "Up Time: " + stopWatch).append("Logger: " + logger.getName());
        return toStringBuilder.toString();
    }

    public static Map<String, List<String>> listGroovyScripts(
            boolean logScripts, boolean logTestScripts) throws IOException {

        // Find the location of ffx-all.jar
        ClassLoader classLoader = ClassLoader.getSystemClassLoader();
        URL scriptURL = classLoader.getResource("ffx/scripts");
        String ffx = scriptURL.getPath().substring(5, scriptURL.getPath().indexOf("!"));

        JarFile jar = new JarFile(URLDecoder.decode(ffx, "UTF-8"));

        Map<String, List<String>> classes = new HashMap<>();
        List<String> scripts = new ArrayList<>();
        List<String> testScripts = new ArrayList<>();
        classes.put("scripts", scripts);
        classes.put("testScripts", testScripts);

        // Getting the files into the jar
        Enumeration<? extends JarEntry> enumeration = jar.entries();

        // Iterates into the files in the jar file
        while (enumeration.hasMoreElements()) {
            ZipEntry zipEntry = enumeration.nextElement();
            String className = zipEntry.getName();

            // Find FFX groovy scripts.
            if (className.startsWith("ffx") && className.endsWith(".groovy")) {
                className = className.replace(".groovy", "").replace("/", ".");
                className = className.replace("ffx.scripts.", "");
                if (className.toUpperCase().contains("TEST")) {
                    testScripts.add(className);
                } else {
                    scripts.add(className);
                }
            }
        }

        // Sort the scripts alphabetically.
        Collections.sort(scripts);
        Collections.sort(testScripts);

        // Log the script names.
        if (logTestScripts) {
            for (String script : testScripts) {
                logger.info(" " + script);
            }
        }
        if (logScripts) {
            for (String script : scripts) {
                logger.info(" " + script);
            }
        }

        return classes;
    }

    /**
     * Returns the file name of a temporary copy of <code>input</code> content.
     *
     * @param input  a {@link java.io.InputStream} object.
     * @param name   a {@link java.lang.String} object.
     * @param suffix a {@link java.lang.String} object.
     * @return a {@link java.lang.String} object.
     * @throws java.io.IOException if any.
     */
    public static String copyInputStreamToTmpFile(final InputStream input,
                                                  String name, final String suffix) throws IOException {
        File tmpFile = null;

        try {
            name = "ffx." + name + ".";
            tmpFile = File.createTempFile(name, suffix);
        } catch (IOException e) {
            System.out.println(" Could not extract a Force Field X library.");
            System.err.println(e.toString());
            System.exit(-1);
        }

        tmpFile.deleteOnExit();

        OutputStream output = null;
        try {
            output = new BufferedOutputStream(new FileOutputStream(tmpFile));
            byte[] buffer = new byte[8192];
            int size;
            while ((size = input.read(buffer)) != -1) {
                output.write(buffer, 0, size);
            }
        } finally {
            if (input != null) {
                input.close();
            }
            if (output != null) {
                output.close();
            }
        }

        return tmpFile.toString();
    }
}

