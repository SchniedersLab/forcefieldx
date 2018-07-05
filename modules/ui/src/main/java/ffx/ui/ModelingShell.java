/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
package ffx.ui;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextPane;
import javax.swing.JViewport;
import javax.swing.SwingUtilities;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Menu;
import java.awt.MenuBar;
import java.awt.Point;
import java.awt.Rectangle;
import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.prefs.Preferences;

import org.codehaus.groovy.runtime.MethodClosure;

import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Minimize;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.Terminatable;
import ffx.algorithms.integrators.IntegratorEnum;
import ffx.algorithms.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;
import ffx.potential.utils.PotentialsFunctions;

import groovy.ui.Console;

/**
 * The ModelingShell is used to script Multiscale Modeling Routines via the
 * Groovy scripting language. Functionality available through the modeling shell
 * includes the Force Field X API, Java API and Groovy extensions.
 *
 * @author Michael J. Schnieders
 *
 */
public final class ModelingShell extends Console implements AlgorithmListener {

    /**
     * The logger for this class.
     */
    private static final Logger logger = Logger.getLogger(ModelingShell.class.getName());
    /**
     * A reference to the main application container.
     */
    private final MainPanel mainPanel;
    /**
     * The flag interrupted is true if a script is running and the user requests
     * it be canceled.
     */
    private boolean interrupted;
    /**
     * The flag headless is true for the CLI and false for the GUI.
     */
    private final boolean headless;
    /**
     * An algorithm that implements the Terminatable interface can be cleanly
     * terminated before completion. For example, after completion of an
     * optimization step or MD step.
     */
    private Terminatable terminatableAlgorithm = null;
    /**
     * Flag to indicate if a script is running.
     */
    public boolean scriptRunning;
    /**
     * Timing.
     */
    private long time;
    private long subTime;
    private static final double toSeconds = 1.0e-9;

    private List<String> args;

    /**
     * <p>
     * Constructor for ModelingShell.</p>
     *
     * @param mainPanel a reference to the MainPanel.
     */
    public ModelingShell(MainPanel mainPanel) {
        super();
        this.mainPanel = mainPanel;
        headless = java.awt.GraphicsEnvironment.isHeadless();
        initContext();
        loadPrefs();
    }

    /**
     * Initialize access to Force Field X variables and methods from with the
     * Shell.
     */
    private void initContext() {
        setVariable("dat", mainPanel.getHierarchy());
        setVariable("cmd", mainPanel);
        setVariable("vis", mainPanel.getGraphics3D());
        setVariable("sh", this);
        setVariable("active", mainPanel.getHierarchy().getActive());
        setVariable("logger", logger);

        // Timer
        setVariable("time", new MethodClosure(this, "time"));

        // File
        setVariable("open", new MethodClosure(mainPanel, "openWait"));
        setVariable("convertWait", new MethodClosure(mainPanel, "convertWait"));
        setVariable("save", new MethodClosure(mainPanel, "saveAsXYZ"));
        setVariable("saveAsXYZ", new MethodClosure(mainPanel, "saveAsXYZ"));
        setVariable("saveAsP1", new MethodClosure(mainPanel, "saveAsP1"));
        setVariable("saveAsPDB", new MethodClosure(mainPanel, "saveAsPDB"));
        setVariable("close", new MethodClosure(mainPanel, "closeWait"));
        setVariable("closeAll", new MethodClosure(mainPanel, "closeAll"));

        // Select
        setVariable("select", new MethodClosure(this, "select"));

        // Display and View menus.
        if (!headless) {
            GraphicsCanvas graphics = mainPanel.getGraphics3D();

            // Display
            int index = 0;
            for (ViewModel view : ViewModel.values()) {
                setVariable(view.name(), view);
                index++;
                if (index > 8) {
                    break;
                }
            }
            setVariable("view", new MethodClosure(graphics, "viewWait"));

            // Color
            index = 0;
            for (ColorModel color : ColorModel.values()) {
                setVariable(color.name(), color);
                index++;
                if (index > 6) {
                    break;
                }
            }
            setVariable("color", new MethodClosure(graphics, "colorWait"));
        }

        // Algorithms
        setVariable("returnEnergy", new MethodClosure(this, "returnEnergy"));
        setVariable("energy", new MethodClosure(this, "energy"));
        setVariable("analyze", new MethodClosure(this, "analyze"));
        setVariable("minimize", new MethodClosure(this, "minimize"));
        setVariable("minimize_2", new MethodClosure(this, "minimize_2"));
        setVariable("md", new MethodClosure(this, "md"));
        setVariable("potential", new MethodClosure(this, "potential"));
        setVariable("poledit", new MethodClosure(this, "poledit"));
        setVariable("superpose", new MethodClosure(this, "superpose"));

        // Obtain UIUtils object
        setVariable("functions", new UIUtils(this, mainPanel));

        // Define a listener variable to send updates back to the GUI.
        setVariable("listener", this);

    }

    /**
     * Update the shell menu items.
     */
    private void initMenus() {
        JFrame frame = (JFrame) this.getFrame();
        MenuBar menuBar = frame.getMenuBar();
        /**
         * Remove "Capture Std. Out", "Capture Std. Error" & "Detached Output"
         * from the View menu.
         */
        Menu menu = menuBar.getMenu(2);
        menu.remove(5);
        menu.remove(5);
        menu.remove(9);
        /**
         * Edit the Script menu.
         */
        menu = menuBar.getMenu(4);
        menu.remove(4);
        menu.remove(4);
        menu.remove(4);
        menu.remove(5);
        menu.remove(7);
    }

    /**
     * Configure the Swing GUI for the shell.
     */
    @Override
    public void run() {
        if (!headless) {
            try {
                super.run();
                // Output JTextPane
                JTextPane output = getOutputArea();
                output.setBackground(Color.BLACK);
                output.setForeground(Color.WHITE);
                // Input JTextPane
                JTextPane input = getInputArea();
                input.setBackground(Color.WHITE);
                input.setForeground(Color.BLACK);
                // Output StyledDocument Styles
                StyledDocument doc = output.getStyledDocument();
                Style defStyle = StyleContext.getDefaultStyleContext().getStyle(StyleContext.DEFAULT_STYLE);
                Style regular = doc.addStyle("regular", defStyle);
                Style prompt = doc.addStyle("prompt", regular);
                Style command = doc.addStyle("command", regular);
                Style result = doc.addStyle("result", regular);
                StyleConstants.setFontFamily(regular, "Monospaced");
                setPromptStyle(prompt);
                setCommandStyle(command);
                setResultStyle(result);
                StyleConstants.setForeground(prompt, Color.ORANGE);
                StyleConstants.setForeground(command, Color.GREEN);
                StyleConstants.setForeground(result, Color.GREEN);
                StyleConstants.setBackground(result, Color.BLACK);
                clearOutput();
                //initMenus();

                // Set labels and icon for Force Field X.
                getStatusLabel().setText("Welcome to the Force Field X Shell.");
                JFrame frame = (JFrame) this.getFrame();
                frame.setTitle("Force Field X Shell");
                URL iconURL = getClass().getClassLoader().getResource(
                        "ffx/ui/icons/icon64.png");
                ImageIcon icon = new ImageIcon(iconURL);
                frame.setIconImage(icon.getImage());
                frame.setSize(600, 600);
            } catch (Exception e) {
                System.out.println(e.toString());
                e.printStackTrace();
                logger.warning(e.toString());
            }
        }
    }

    /**
     * <p>
     * setArgList</p>
     *
     * @param argList a {@link java.util.List} object.
     */
    public void setArgList(List<String> argList) {
        args = new ArrayList<>(argList);
        setVariable("args", argList);
    }

    List<String> getArgs() {
        return new ArrayList<>(args);
    }

    /**
     * <p>
     * runFFXScript</p>
     *
     * @param file a {@link java.io.File} object.
     */
    public void runFFXScript(File file) {
        try {
            before();
            try {
                getShell().run(file, args);
            } catch (Exception ex) {
                logger.log(Level.SEVERE, " Uncaught error: FFX is shutting down.\n", ex);
            }
            after();
        } catch (Exception e) { // Replacing this with a "Multi-Catch" leads to specific Exceptions not present in some versions of Groovy.
            String message = "Error evaluating script.";
            logger.log(Level.WARNING, message, e);
        }
    }

    /**
     * <p>
     * time</p>
     *
     * @return a {@link java.lang.Double} object.
     */
    public Double time() {
        long current = System.nanoTime();
        double timer = (current - subTime) * toSeconds;
        subTime = current;
        appendOutput(String.format("\n Intermediate time: %8.3f (sec)\n", timer), getPromptStyle());
        return timer;
    }

    /**
     * <p>
     * select</p>
     *
     * @param node a {@link ffx.potential.bonded.MSNode} object.
     */
    public void select(MSNode node) {
        if (node != null) {
            mainPanel.getHierarchy().onlySelection(node);
            sync();
            logger.info(String.format(" Selected: %s.", node));
        }
    }

    /**
     * <p>
     * energy</p>
     *
     * @return a {@link ffx.potential.ForceFieldEnergy} object.
     */
    public ForceFieldEnergy energy() {
        if (interrupted) {
            logger.info(" Algorithm interrupted - skipping energy.");
            return null;
        }
        if (terminatableAlgorithm != null) {
            logger.info(" Algorithm already running - skipping energy.");
            return null;
        }

        MolecularAssembly active = mainPanel.getHierarchy().getActive();
        if (active != null) {
            ForceFieldEnergy energy = active.getPotentialEnergy();
            if (energy == null) {
                energy = ForceFieldEnergy.energyFactory(active);
                active.setPotential(energy);
            }
            energy.energy(false, true);
            return energy;
        }
        return null;
    }

    /**
     * <p>
     * returnEnergy</p>
     *
     * @return Current system energy (a double).
     */
    public double returnEnergy() {
        if (interrupted) {
            logger.info(" Algorithm interrupted - skipping energy.");
            return 0.0;
        }
        if (terminatableAlgorithm != null) {
            logger.info(" Algorithm already running - skipping energy.");
            return 0.0;
        }

        MolecularAssembly active = mainPanel.getHierarchy().getActive();
        if (active != null) {
            ForceFieldEnergy energy = active.getPotentialEnergy();
            if (energy == null) {
                energy = ForceFieldEnergy.energyFactory(active);
                active.setPotential(energy);
            }
            return energy.energy(false, true);
        }
        logger.warning(" Energy could not be calculated");
        return 0.0;
    }

    public AlgorithmFunctions getUIAlgorithmUtils() {
        return new UIUtils(this, mainPanel);
    }

    public PotentialsFunctions getUIPotentialsUtils() {
        return new UIUtils(this, mainPanel);
    }

    /**
     * <p>
     * minimize</p>
     *
     * @param eps a double.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps) {
        if (interrupted) {
            logger.info(" Algorithm interrupted - skipping minimization.");
            return null;
        }
        if (terminatableAlgorithm != null) {
            logger.info(" Algorithm already running - skipping minimization.");
            return null;
        }
        MolecularAssembly active = mainPanel.getHierarchy().getActive();
        if (active != null) {
            Minimize minimize = new Minimize(active, this);
            terminatableAlgorithm = minimize;
            Potential potential = minimize.minimize(eps);
            terminatableAlgorithm = null;
            return potential;
        } else {
            logger.info(" No active system to minimize.");
        }
        return null;
    }

    /**
     * <p>
     * md</p>
     *
     * @param nStep a int.
     * @param timeStep a double.
     * @param printInterval a double.
     * @param saveInterval a double.
     * @param temperature a double.
     * @param initVelocities a boolean.
     * @param dyn a {@link java.io.File} object.
     */
    public void md(int nStep, double timeStep, double printInterval,
                   double saveInterval, double temperature, boolean initVelocities,
                   File dyn) {
        if (interrupted || terminatableAlgorithm != null) {
            return;
        }
        FFXSystem active = mainPanel.getHierarchy().getActive();
        if (active != null) {
            MolecularDynamics molecularDynamics = new MolecularDynamics(active,
                    active.getPotentialEnergy(),
                    active.getProperties(),
                    this, ThermostatEnum.BUSSI, IntegratorEnum.BEEMAN);
            terminatableAlgorithm = molecularDynamics;
            molecularDynamics.dynamic(nStep, timeStep, printInterval,
                    saveInterval, temperature, initVelocities, dyn);
            terminatableAlgorithm = null;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Return false if user elects to cancel.
     */
    @Override
    public boolean askToSaveFile() {
        File file = (File) getScriptFile();
        if (file == null || !getDirty()) {
            return true;
        }
        switch (JOptionPane.showConfirmDialog((Component) getFrame(),
                "Save changes to " + file.getName() + "?",
                "Force Field X Shell", JOptionPane.YES_NO_CANCEL_OPTION)) {
            case JOptionPane.YES_OPTION:
                return fileSave();
            case JOptionPane.NO_OPTION:
                return true;
            default:
                return false;
        }
    }

    /**
     * {@inheritDoc}
     *
     * Print out the Force Field X promo.
     */
    @Override
    final public void clearOutput() {
        if (!java.awt.GraphicsEnvironment.isHeadless()) {
            JTextPane output = getOutputArea();
            output.setText("");
            appendOutput(MainPanel.border + "\n" + MainPanel.title
                    + MainPanel.aboutString + "\n" + MainPanel.border + "\n", getCommandStyle());
        }
    }

    @Override
    public void clearOutput(EventObject evt) {
        clearOutput();
    }

    @Override
    public void fileNewWindow() {
        mainPanel.resetShell();
    }

    @Override
    public void fileNewWindow(EventObject evt) {
        fileNewWindow();
    }

    @Override
    public void showAbout() {
        mainPanel.about();
    }

    @Override
    public void showAbout(EventObject evt) {
        showAbout();
    }

    /**
     * Clear output text from any previous script and then output a message
     * about the new script.
     */
    private void scriptStartup() {
        clearOutput();

        /**
         * Attempt to get the script's name.
         */
        Object name = getScriptFile();
        if (name != null && name instanceof File) {
            name = ((File) name).getName();
        }
        /**
         * A short message about the script to be evaluated.
         */
        String message;
        if (name == null || name.toString().equalsIgnoreCase("null")) {
            message = String.format("\n Evaluating...\n\n");
        } else {
            message = String.format("\n Evaluating " + name + "...\n\n");
        }
        appendOutput(message, getPromptStyle());
    }

    /**
     * <p>
     * before</p>
     */
    public void before() {
        interrupted = false;
        terminatableAlgorithm = null;
        time = System.nanoTime();
        subTime = time;
        mainPanel.getModelingPanel().enableLaunch(false);
    }

    /**
     * <p>
     * after</p>
     */
    public void after() {
        time = System.nanoTime() - time;
        if (!interrupted) {
            appendOutput(String.format("\n Script wall clock time: %6.3f (sec)", time * toSeconds), getPromptStyle());
        } else {
            appendOutput(String.format("\n Script interrupted after: %6.3f (sec)", time * toSeconds), getPromptStyle());
        }
        mainPanel.getModelingPanel().enableLaunch(true);
    }

    /**
     * Fix up the "Result: " message, then call the original method.
     *
     * @param string String to ouput.
     * @param style Style to use.
     */
    public void appendOutputNl(String string, Style style) {
        if (interrupted) {
            return;
        }
        if (headless) {
            logger.info(string);
        }
        if (string.equals("Result: ")) {
            string = " Script result: \n";
        } else if (string.equals("groovy> ")) {
            string = " ffx> ";
        }

        super.appendOutputNl(string, style);

        if (EventQueue.isDispatchThread()) {
            scroll();
        } else {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    scroll();
                }
            });
        }
    }

    /**
     * <p>
     * scroll</p>
     */
    public void scroll() {
        JTextPane output = getOutputArea();
        JSplitPane splitPane = getSplitPane();
        JScrollPane scrollPane = (JScrollPane) splitPane.getBottomComponent();
        JViewport viewport = scrollPane.getViewport();
        Rectangle visibleSize = viewport.getVisibleRect();
        Dimension totalSize = output.getSize();
        Point point = new Point(0, totalSize.height - visibleSize.height);
        viewport.setViewPosition(point);
    }

    /**
     * <p>
     * appendOutput</p>
     *
     * @param string a {@link java.lang.String} object.
     * @param style a {@link javax.swing.text.Style} object.
     */
    public void appendOutput(String string, Style style) {
        if (interrupted) {
            return;
        }

        if (headless) {
            logger.info(string);
            return;
        }

        super.appendOutput(string, style);
        if (EventQueue.isDispatchThread()) {
            scroll();
        } else {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    scroll();
                }
            });
        }
    }

    public void interruptScript() {
        if (!scriptRunning) {
            return;
        }
        doInterrupt();
    }

    /**
     * If at exit time, a script is running, the user is given an option to
     * interrupt it first
     *
     * {@inheritDoc}
     */
    @Override
    public Object askToInterruptScript() {
        if (!scriptRunning) {
            return true;
        }
        int rc = JOptionPane.showConfirmDialog(getScrollArea(),
                "Script executing. Press 'OK' to attempt to interrupt it before exiting.",
                "Force Field X Shell", JOptionPane.OK_CANCEL_OPTION);
        if (rc == JOptionPane.OK_OPTION) {
            doInterrupt();
            return true;
        } else {
            return false;
        }
    }

    private static final Preferences preferences = Preferences.userNodeForPackage(ModelingShell.class);

    /**
     * <p>
     * loadPrefs</p>
     */
    final public void loadPrefs() {
    }

    /**
     * <p>
     * savePrefs</p>
     */
    public void savePrefs() {
    }

    /**
     * <p>
     * setMeasurement</p>
     *
     * @param measurement a {@link java.lang.String} object.
     * @param d a double.
     */
    public void setMeasurement(String measurement, double d) {
        try {
            appendOutput(measurement, getOutputStyle());
        } catch (Exception e) {
            String message = "Exception appending measurement to Shell.\n";
            logger.log(Level.WARNING, message, e);
        }
    }

    /**
     * <p>
     * sync</p>
     */
    public void sync() {
        try {
            setVariable("active", mainPanel.getHierarchy().getActive());
        } catch (Exception e) {
            String message = " Exception syncing shell variables.\n";
            logger.log(Level.WARNING, message, e);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        return "Force Field X Shell";
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean algorithmUpdate(MolecularAssembly active) {
        if (interrupted) {
            return false;
        }

        GraphicsCanvas graphics = mainPanel.getGraphics3D();
        if (graphics != null) {
            /**
             * Use the blocking graphics update method so that only
             * self-consistent coordinate sets are displayed.
             */
            //if (SwingUtilities.isEventDispatchThread()) {
            graphics.updateSceneWait(active, true, false, null, false, null);
            //}
        }

        if (interrupted) {
            /**
             * The algorithm could have been interrupted during the graphics
             * update.
             */
            return false;
        } else {
            return true;
        }
    }

    @Override
    public void clearContext() {
        super.clearContext();
        initContext();
    }

    @Override
    public void clearContext(EventObject evt) {
        super.clearContext(evt);
        initContext();
    }

    @Override
    public void updateTitle() {
        JFrame frame = (JFrame) this.getFrame();
        File file = (File) getScriptFile();
        if (file != null) {
            String name = file.getName();
            frame.setTitle(name + " - Force Field X Shell");
        } else {
            frame.setTitle("Force Field X Shell");
        }
    }
}
