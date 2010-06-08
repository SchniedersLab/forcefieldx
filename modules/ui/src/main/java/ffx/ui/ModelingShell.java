/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
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
import java.awt.Component;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Rectangle;
import java.io.File;
import java.net.URL;
import java.util.EventObject;
import java.util.prefs.Preferences;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.Box;
import javax.swing.ImageIcon;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextPane;
import javax.swing.SwingUtilities;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;
import javax.swing.text.StyledDocument;

import groovy.ui.Console;

import org.codehaus.groovy.runtime.MethodClosure;

import ffx.algorithms.Minimize;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.MolecularDynamics;
import ffx.algorithms.Terminatable;
import ffx.algorithms.Thermostat.Thermostats;
import ffx.numerics.Optimizable;

import ffx.potential.PotentialEnergy;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;

/**
 * The ModelingShell is used to script Multiscale Modeling Routines via the
 * Groovy scripting language. Functionality available through the modeling shell
 * includes the Force Field X API, Java API and Groovy extensions.
 */
public class ModelingShell extends Console implements AlgorithmListener {

    private static final Logger logger = Logger.getLogger(ModelingShell.class.getName());
    private static final long serialVersionUID = 1L;
    private MainPanel mainPanel;
    private boolean interrupted;
    private boolean headless;
    private Terminatable terminatableAlgorithm = null;
    private Rectangle outputSize = new Rectangle(0, 0, 0, 0);
    private Rectangle scrollTo = new Rectangle(0, 0, 0, 0);
    private long time;
    private long subTime;
    public boolean scriptRunning;
    private static final double toSeconds = 1.0e-9;

    public ModelingShell(MainPanel m) {
        mainPanel = m;
        headless = java.awt.GraphicsEnvironment.isHeadless();

        if (!headless) {
            run();
        }

        setVariable("dat", mainPanel.getHierarchy());
        setVariable("cmd", mainPanel);
        setVariable("vis", mainPanel.getGraphics3D());
        setVariable("sh", this);
        setVariable("active", mainPanel.getHierarchy().getActive());

        // Timer
        setVariable("time", new MethodClosure(this, "time"));

        // File
        setVariable("open", new MethodClosure(mainPanel, "openWait"));
        setVariable("save", new MethodClosure(mainPanel, "saveAsXYZ"));
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
        setVariable("energy", new MethodClosure(this, "energy"));
        setVariable("minimize", new MethodClosure(this, "minimize"));
        setVariable("md", new MethodClosure(this, "md"));

        /**
         * Configure the Swing GUI for the shell.
         */
        if (!headless) {

            JTextPane output = getOutputArea();
            outputSize = output.getVisibleRect();

            output.setBackground(Color.BLACK);
            output.setForeground(Color.WHITE);
            JTextPane input = getInputArea();
            input.setBackground(Color.WHITE);
            input.setForeground(Color.BLACK);

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

            getStatusLabel().setText("Welcome to the Force Field X Shell.");

            JFrame frame = (JFrame) this.getFrame();
            frame.setTitle("Force Field X Shell");
            URL iconURL = getClass().getClassLoader().getResource(
                    "ffx/ui/icons/icon64.png");
            ImageIcon icon = new ImageIcon(iconURL);
            frame.setIconImage(icon.getImage());

            clearOutput();
        }

        loadPrefs();
    }

    public void headlessRun(File file) {
        try {
            before();
            getShell().evaluate(file);
            after();
        } catch (Exception e) {
            String message = "Error evaluating script.";
            logger.log(Level.WARNING, message, e);
        }
    }

    public Double time() {
        long current = System.nanoTime();
        double timer = (current - subTime) * toSeconds;
        subTime = current;
        appendOutput(String.format("\n Intermediate time: %8.3f (sec)\n", timer), getPromptStyle());
        return timer;
    }

    public void select(MSNode node) {
        if (node != null) {
            mainPanel.getHierarchy().onlySelection(node);
            sync();
            logger.info(" Selected: " + node);
        }
    }

    public PotentialEnergy energy() {
        if (interrupted) {
            return null;
        }
        MolecularAssembly active = mainPanel.getHierarchy().getActive();
        if (active != null) {
            PotentialEnergy energy = active.getPotentialEnergy();
            if (energy == null) {
                energy = new PotentialEnergy(active);
                active.setPotential(energy);
            }
            energy.energy(false, true);
            return energy;
        }
        return null;
    }

    public Optimizable minimize(double eps) {
        if (interrupted || terminatableAlgorithm != null) {
            return null;
        }
        MolecularAssembly active = mainPanel.getHierarchy().getActive();
        if (active != null) {
            Minimize minimize = new Minimize(active, this);
            terminatableAlgorithm = minimize;
            Optimizable energy = minimize.minimize(eps);
            terminatableAlgorithm = null;
            return energy;
        }
        return null;
    }

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
                                                                        this, Thermostats.BUSSI);
            terminatableAlgorithm = molecularDynamics;
            molecularDynamics.dynamic(nStep, timeStep, printInterval,
                                      saveInterval, temperature, initVelocities, dyn);
            terminatableAlgorithm = null;
        }
    }

    /*
    @Override
    public void runScript() {
    runScript(null);
    }

    @Override
    public void runScript(EventObject evt) {
    scriptStartup();
    super.runScript(evt);
    }

    @Override
    public void runSelectedScript() {
    runSelectedScript(null);
    }

    @Override
    public void runSelectedScript(EventObject evt) {
    scriptStartup();
    super.runSelectedScript(evt);
    } */
    /**
     * Return false if user elects to cancel.
     *
     * @return The result of the <code>fileSave</code> method.
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
     * Print out the Force Field X promo.
     */
    @Override
    final public void clearOutput() {
        if (!java.awt.GraphicsEnvironment.isHeadless()) {
            JTextPane output = getOutputArea();
            output.setText("");
            appendOutput(MainPanel.border + MainPanel.title
                         + MainPanel.aboutString + "\n" + MainPanel.border, getCommandStyle());
        }
    }

    /*
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

    @Override
    public void updateTitle() {
    JFrame frame = (JFrame) getFrame();
    frame.setTitle("Force Field X Shell");
    } */
    private void scriptStartup() {
        clearOutput();
        Object name = getScriptFile();
        if (name != null && name instanceof File) {
            name = ((File) name).getName();
        }
        if (name == null || name.toString().equalsIgnoreCase("null")) {
            appendOutput(String.format("\n Evaluating...\n\n"), getPromptStyle());
        } else {
            appendOutput(String.format("\n Evaluating " + name + "...\n\n"), getPromptStyle());
        }
    }

    public void before() {
        interrupted = false;
        terminatableAlgorithm = null;
        time = System.nanoTime();
        subTime = time;
        if (!headless) {
            outputSize = getOutputArea().getVisibleRect();
            scrollTo.x = 0;
            scrollTo.height = outputSize.height;
            scrollTo.width = outputSize.width;
        }
    }

    public void after() {
        time = System.nanoTime() - time;
        if (!interrupted) {
            appendOutput(String.format("\n Total script time: %8.3f (sec)", time * toSeconds), getPromptStyle());
        }
    }

    /**
     * Fix up the "Result: " message, then call the original method.
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

    public void scroll() {
        JTextPane output = getOutputArea();
        Dimension dim = output.getSize();
        scrollTo.y = dim.height - outputSize.y;
        output.scrollRectToVisible(scrollTo);
    }

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
    /**
     * Confirm whether to interrupt the running thread.
     * @param evt
     */
    /* @Override
    public void confirmRunInterrupt(EventObject evt) {
    int rc = JOptionPane.showConfirmDialog((Component) getFrame(), "Attempt to interrupt script?",
    "Force Field X Shell", JOptionPane.YES_NO_OPTION);
    if (rc == JOptionPane.YES_OPTION) {
    interrupted = true;
    if (terminatableAlgorithm != null) {
    terminatableAlgorithm.terminate();
    }
    Thread thread = getRunThread();
    if (thread != null) {
    thread.interrupt();
    }
    }
    } */
    private static final Preferences preferences = Preferences.userNodeForPackage(ModelingShell.class);

    final public void loadPrefs() {
    }

    public void savePrefs() {
    }

    public void setMeasurement(String measurement, double d) {
        try {
            appendOutput(measurement, getOutputStyle());
        } catch (Exception e) {
            String message = "Exception appending measurement to Shell.\n";
            logger.log(Level.WARNING, message, e);
        }
    }

    public void sync() {
        try {
            setVariable("active", mainPanel.getHierarchy().getActive());
        } catch (Exception e) {
            String message = "Exception syncing shell variables.\n";
            logger.log(Level.WARNING, message, e);
        }
    }

    @Override
    public String toString() {
        return "Force Field X Shell";
    }

    @Override
    public boolean algorithmUpdate(MolecularAssembly active) {
        if (interrupted) {
            return !interrupted;
        }

        GraphicsCanvas graphics = mainPanel.getGraphics3D();
        if (graphics != null) {
            /**
             * Use the blocking graphics update method so that only self-consistent
             * coordinate sets are displayed.
             */
            graphics.updateSceneWait(active, true, false, null, false, null);
        }

        return !interrupted;
    }
}
