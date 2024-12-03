// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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

import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.algorithms.cli.AlgorithmsScript;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.algorithms.optimize.Minimize;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.RendererCache.ColorModel;
import ffx.potential.bonded.RendererCache.ViewModel;
import ffx.potential.cli.PotentialScript;
import ffx.potential.utils.PotentialsFunctions;
import ffx.utilities.Console;
import ffx.utilities.GroovyFileFilter;
import groovy.lang.Binding;
import groovy.lang.Script;

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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
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

import org.codehaus.groovy.runtime.MethodClosure;

import org.graalvm.polyglot.Context;
import org.graalvm.polyglot.Engine;
import org.graalvm.polyglot.Source;
import org.graalvm.polyglot.Value;
import org.graalvm.polyglot.Language;
import org.graalvm.polyglot.proxy.ProxyArray;

/**
 * The ModelingShell is used to script Multiscale Modeling Routines via the Groovy scripting
 * language. Functionality available through the modeling shell includes the Force Field X API, Java
 * API and Groovy extensions.
 *
 * @author Michael J. Schnieders
 */
public class ModelingShell extends Console implements AlgorithmListener {

  /**
   * The logger for this class.
   */
  private static final Logger logger = Logger.getLogger(ModelingShell.class.getName());
  private static final double toSeconds = 1.0e-9;

  /**
   * A reference to the main application container.
   */
  private final MainPanel mainPanel;
  /**
   * The flag headless is true for the CLI and false for the GUI.
   */
  private final boolean headless;
  /**
   * The flag interrupted is true if a script is running and the user requests it be canceled.
   */
  private boolean interrupted;
  /**
   * An algorithm that implements the Terminatable interface can be cleanly terminated before
   * completion. For example, after completion of an optimization step or MD step.
   */
  private Terminatable terminatableAlgorithm = null;
  /**
   * Flag to indicate if a script is running.
   */
  private boolean scriptRunning;
  /**
   * Timing.
   */
  private long time;

  private long subTime;
  private List<String> args;

  /**
   * Constructor for ModelingShell.
   *
   * @param mainPanel a reference to the MainPanel.
   */
  public ModelingShell(MainPanel mainPanel) {
    super();
    this.mainPanel = mainPanel;
    headless = java.awt.GraphicsEnvironment.isHeadless();
    initContext(getShell().getContext());
  }

  /**
   * after
   */
  public void after() {
    time = System.nanoTime() - time;
    scriptRunning = false;
    if (!interrupted) {
      appendOutput(String.format("\n Script wall clock time: %6.3f (sec)", time * toSeconds),
          getPromptStyle());
    } else {
      appendOutput(String.format("\n Script interrupted after: %6.3f (sec)", time * toSeconds),
          getPromptStyle());
    }
    mainPanel.getModelingPanel().enableLaunch(true);
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
      /*
       Use the blocking graphics update method so that only
       self-consistent coordinate sets are displayed.
      */
      // if (SwingUtilities.isEventDispatchThread()) {
      graphics.updateSceneWait(active, true, false, null, false, null);
      // }
    }

    // The algorithm could have been interrupted during the graphics update.
    return !interrupted;
  }

  /**
   * If at exit time, a script is running, the user is given an option to interrupt it first
   *
   * <p>{@inheritDoc}
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

  /**
   * {@inheritDoc}
   *
   * <p>Return false if user elects to cancel.
   */
  @Override
  public boolean askToSaveFile() {
    File file = (File) getScriptFile();
    if (file == null || !getDirty()) {
      return true;
    }
    return switch (JOptionPane.showConfirmDialog((Component) getFrame(),
        "Save changes to " + file.getName() + "?", "Force Field X Shell",
        JOptionPane.YES_NO_CANCEL_OPTION)) {
      case JOptionPane.YES_OPTION -> fileSave();
      case JOptionPane.NO_OPTION -> true;
      default -> false;
    };
  }

  /**
   * before
   */
  public void before() {
    interrupted = false;
    terminatableAlgorithm = null;
    time = System.nanoTime();
    subTime = time;
    mainPanel.getModelingPanel().enableLaunch(false);
    scriptRunning = true;
  }

  @Override
  public void clearContext() {
    super.clearContext();
    initContext(getShell().getContext());
  }

  @Override
  public void clearContext(EventObject evt) {
    super.clearContext(evt);
    initContext(getShell().getContext());
  }

  /**
   * {@inheritDoc}
   *
   * <p>Print out the Force Field X promo.
   */
  @Override
  public void clearOutput() {
    if (!java.awt.GraphicsEnvironment.isHeadless()) {
      JTextPane output = getOutputArea();
      output.setText("");
      appendOutput(
          MainPanel.border + "\n" + MainPanel.title + MainPanel.aboutString + "\n" + MainPanel.border
              + "\n", getCommandStyle());
    }
  }

  @Override
  public void clearOutput(EventObject evt) {
    clearOutput();
  }

  /**
   * energy
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

  @Override
  public void fileNewWindow() {
    mainPanel.resetShell();
  }

  @Override
  public void fileNewWindow(EventObject evt) {
    fileNewWindow();
  }

  public AlgorithmFunctions getUIAlgorithmUtils() {
    return new UIUtils(this, mainPanel);
  }

  public PotentialsFunctions getUIPotentialsUtils() {
    return new UIUtils(this, mainPanel);
  }

  /**
   * md
   *
   * @param nStep          The number of MD steps.
   * @param timeStep       a double.
   * @param printInterval  a double.
   * @param saveInterval   a double.
   * @param temperature    a double.
   * @param initVelocities a boolean.
   * @param dyn            a {@link java.io.File} object.
   */
  public void md(int nStep, double timeStep, double printInterval, double saveInterval,
                 double temperature, boolean initVelocities, File dyn) {
    if (interrupted || terminatableAlgorithm != null) {
      return;
    }
    FFXSystem active = mainPanel.getHierarchy().getActive();
    if (active != null) {
      MolecularDynamics molecularDynamics = new MolecularDynamics(active,
          active.getPotentialEnergy(), this, ThermostatEnum.BUSSI, IntegratorEnum.BEEMAN);
      terminatableAlgorithm = molecularDynamics;
      molecularDynamics.dynamic(nStep, timeStep, printInterval, saveInterval, temperature,
          initVelocities, dyn);
      terminatableAlgorithm = null;
    }
  }

  /**
   * Configure the Swing GUI for the shell.
   */
  @Override
  public void run() {
    if (!headless) {
      if (SwingUtilities.isEventDispatchThread()) {
        init();
      } else {
        try {
          SwingUtilities.invokeAndWait(this::init);
        } catch (Exception e) {
          //
        }
      }
    }
  }

  /**
   * select
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
   * setArgList
   *
   * @param argList a {@link java.util.List} object.
   */
  public void setArgList(List<String> argList) {
    args = new ArrayList<>(argList);
    setVariable("args", argList);
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
   * time
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
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    return "Force Field X Shell";
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

  /**
   * Initialize access to Force Field X variables and methods from with the Shell.
   */
  private void initContext(Binding binding) {
    binding.setVariable("dat", mainPanel.getHierarchy());
    binding.setVariable("cmd", mainPanel);
    binding.setVariable("vis", mainPanel.getGraphics3D());
    binding.setVariable("sh", this);
    binding.setVariable("active", mainPanel.getHierarchy().getActive());
    binding.setVariable("logger", logger);

    // Timer
    binding.setVariable("time", new MethodClosure(this, "time"));

    // File
    binding.setVariable("open", new MethodClosure(mainPanel, "openWait"));
    binding.setVariable("convertWait", new MethodClosure(mainPanel, "convertWait"));
    binding.setVariable("save", new MethodClosure(mainPanel, "saveAsXYZ"));
    binding.setVariable("saveAsXYZ", new MethodClosure(mainPanel, "saveAsXYZ"));
    binding.setVariable("saveAsP1", new MethodClosure(mainPanel, "saveAsP1"));
    binding.setVariable("saveAsPDB", new MethodClosure(mainPanel, "saveAsPDB"));
    binding.setVariable("close", new MethodClosure(mainPanel, "closeWait"));
    binding.setVariable("closeAll", new MethodClosure(mainPanel, "closeAll"));

    // Select
    binding.setVariable("select", new MethodClosure(this, "select"));

    // Display and View menus.
    if (!headless) {
      GraphicsCanvas graphics = mainPanel.getGraphics3D();

      // Display
      int index = 0;
      for (ViewModel view : ViewModel.values()) {
        binding.setVariable(view.name(), view);
        index++;
        if (index > 8) {
          break;
        }
      }
      binding.setVariable("view", new MethodClosure(graphics, "viewWait"));

      // Color
      index = 0;
      for (ColorModel color : ColorModel.values()) {
        binding.setVariable(color.name(), color);
        index++;
        if (index > 6) {
          break;
        }
      }
      binding.setVariable("color", new MethodClosure(graphics, "colorWait"));
    }

    // Algorithms
    binding.setVariable("returnEnergy", new MethodClosure(this, "returnEnergy"));
    binding.setVariable("energy", new MethodClosure(this, "energy"));
    binding.setVariable("analyze", new MethodClosure(this, "analyze"));
    binding.setVariable("minimize", new MethodClosure(this, "minimize"));
    binding.setVariable("minimize_2", new MethodClosure(this, "minimize_2"));
    binding.setVariable("md", new MethodClosure(this, "md"));
    binding.setVariable("potential", new MethodClosure(this, "potential"));
    binding.setVariable("poledit", new MethodClosure(this, "poledit"));
    binding.setVariable("superpose", new MethodClosure(this, "superpose"));

    // Obtain UIUtils object
    binding.setVariable("functions", new UIUtils(this, mainPanel));

    // Define a listener variable to send updates back to the GUI.
    binding.setVariable("listener", this);
  }

  /**
   * Update the shell menu items.
   */
  private void initMenus() {
    JFrame frame = (JFrame) this.getFrame();
    MenuBar menuBar = frame.getMenuBar();
    // Remove "Capture Std. Out", "Capture Std. Error" & "Detached Output" from the View menu.
    Menu menu = menuBar.getMenu(2);
    menu.remove(5);
    menu.remove(5);
    menu.remove(9);

    // Edit the Script menu.
    menu = menuBar.getMenu(4);
    menu.remove(4);
    menu.remove(4);
    menu.remove(4);
    menu.remove(5);
    menu.remove(7);
  }

  /**
   * Initialize the Shell.
   */
  private void init() {
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
      // initMenus();

      // Set labels and icon for Force Field X.
      getStatusLabel().setText("Welcome to the Force Field X Shell.");
      JFrame frame = (JFrame) this.getFrame();
      frame.setTitle("Force Field X Shell");
      URL iconURL = getClass().getClassLoader().getResource("ffx/ui/icons/icon64.png");
      ImageIcon icon = new ImageIcon(iconURL);
      frame.setIconImage(icon.getImage());
      frame.setSize(600, 600);
    } catch (Exception e) {
      logger.warning(" Exception starting up the FFX console.\n" + e);
    }
  }

  List<String> getArgs() {
    return new ArrayList<>(args);
  }

  /**
   * runFFXScript - Execute a FFX script.
   *
   * @param file    a {@link java.io.File} object.
   * @param argList List of String inputs to the script.
   * @return Returns a reference to the executed script.
   */
  public Script runFFXScript(File file, List<String> argList) {
    GroovyFileFilter groovyFileFilter = new GroovyFileFilter();
    // Check that the file is a Groovy script.
    if (!groovyFileFilter.accept(file)) {
      // If not, assume its a Python script.
      return runNonGroovyScript(file, argList);
    }

    logger.info(" Executing external Groovy script: " + file.getAbsolutePath() + "\n");
    try {
      before();
      Script script = null;
      try {
        // Run the file using the current Shell and its Binding.
        Object o = getShell().run(file, argList);

        if (o instanceof Script) {
          script = (Script) o;
        }

        // Do not destroy the system when using the GUI.
        if (headless) {
          if (o instanceof PotentialScript) {
            ((PotentialScript) o).destroyPotentials();
          } else if (o instanceof AlgorithmsScript) {
            ((AlgorithmsScript) o).destroyPotentials();
          }
        }
      } catch (Exception ex) {
        logger.log(Level.SEVERE, " Uncaught error: FFX is shutting down.\n", ex);
      }
      after();
      return script;
    } catch (Exception e) {
      // Replacing this with a "Multi-Catch" leads to specific Exceptions not present in
      // some versions of Groovy.
      String message = "Error evaluating script.";
      logger.log(Level.WARNING, message, e);
    }

    return null;
  }

  /**
   * runPythonScript - Execute a Python script.
   *
   * @param file    a {@link java.io.File} object.
   * @param argList List of String inputs to the script.
   * @return Returns a reference to the executed script.
   */
  public Script runNonGroovyScript(File file, List<String> argList) {
    logger.info(" Attempting to execute Polyglot script:\n  " + file.getAbsolutePath() + "\n");

    if (logger.isLoggable(Level.FINE)) {
      logger.fine(" Available languages: ");
      try (Engine engine = Engine.create()) {
        Map<String, Language> map = engine.getLanguages();
        for (Map.Entry<String, Language> entry : map.entrySet()) {
          logger.fine("  " + entry.getKey());
        }
      } catch (Exception ex) {
        logger.log(Level.SEVERE, " Uncaught error: FFX is shutting down.\n", ex);
      }
    }

    try {
      before();
      Script script = null;
      String language = Source.findLanguage(file);
      logger.info(" Detected script language: " + language);
      Source source = Source.newBuilder(language, file).build();
      try (Context context = getContext(language)) {

        // Get the bindings for the language.
        Value bindings = context.getBindings(language);
        // Use the Polyglot ProxyArray to pass the command line arguments to the script.
        List<Object> objList = new ArrayList<>(argList);
        ProxyArray argArray = ProxyArray.fromList(objList);
        bindings.putMember("args", argArray);
        // Run the file using the current Shell and its Binding.
        Value result = context.eval(source);
        logger.info(" Execution of Polyglot script completed.");
      } catch (Exception ex) {
        logger.log(Level.SEVERE, " Uncaught error: FFX is shutting down.\n", ex);
      }
      after();
      return script;
    } catch (Exception e) {
      // Replacing this with a "Multi-Catch" leads to specific Exceptions not present in some versions of Groovy.
      String message = "Error evaluating script.";
      logger.log(Level.WARNING, message, e);
    }
    return null;
  }

  /**
   * Create a Polyglot Context for the specified language.
   * @param language a String specifying the language.
   * @return a Polyglot Context.
   */
  private Context getContext(String language) {
    if (language.equalsIgnoreCase("python")) {
      // For Python, try to locate the Graal Python executable.

      // The default location is $FFX_HOME/ffx_venv/bin/graalpy.
      String FFX_HOME = System.getProperty("basedir");
      Path graalpy = Paths.get(FFX_HOME, "ffx_venv", "bin", "graalpy");

      // Override the default location with the -Dgraalpy=path.to.graalpy option.
      String graalpyString = System.getProperty("graalpy", graalpy.toString());
      graalpy = Paths.get(graalpyString);
      if (graalpy.toFile().exists()) {
        logger.info(" graalpy (-Dgraalpy=path.to.graalpy):             " + graalpy);
        return Context.newBuilder(language).allowAllAccess(true).
            option("python.Executable", graalpyString).
            option("python.PythonPath", ".").build();
      }
    }
    // Fall through to default.
    return Context.newBuilder(language).allowAllAccess(true).build();
  }

  /**
   * runFFXScript - Execute a compiled FFX script.
   *
   * @param script  a compiled FFX script.
   * @param argList List of String inputs to the script.
   * @return Returns a reference to the executed script.
   */
  public Script runFFXScript(Class<? extends Script> script, List<String> argList) {
    logger.info(" Executing internal script: " + script.getCanonicalName() + "\n");
    try {
      before();
      Script groovyScript = null;
      try {
        // Create a Binding for command line arguments and FFX User Interface variables.
        Binding binding = new Binding();
        binding.setVariable("args", argList);
        initContext(binding);

        // Create a new instance of the script and run it.
        groovyScript = script.getDeclaredConstructor().newInstance();
        groovyScript.setBinding(binding);
        groovyScript.run();

        // Do not destroy the system when using the GUI.
        if (headless) {
          if (groovyScript instanceof PotentialScript) {
            ((PotentialScript) groovyScript).destroyPotentials();
          } else if (groovyScript instanceof AlgorithmsScript) {
            ((AlgorithmsScript) groovyScript).destroyPotentials();
          }
        }
      } catch (Exception ex) {
        logger.log(Level.SEVERE, " Uncaught error: FFX is shutting down.\n", ex);
      }
      after();

      return groovyScript;
    } catch (Exception e) {
      // Replacing this with a "Multi-Catch" leads to specific Exceptions not present in
      // some versions of Groovy.
      String message = "Error evaluating script.";
      logger.log(Level.WARNING, message, e);
    }
    return null;
  }

  /**
   * returnEnergy
   *
   * @return Current system energy (a double).
   */
  double returnEnergy() {
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

  /**
   * minimize
   *
   * @param eps a double.
   * @return a {@link ffx.numerics.Potential} object.
   */
  Potential minimize(double eps) {
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
   * Fix up the "Result: " message, then call the original method.
   *
   * @param string String to ouput.
   * @param style  Style to use.
   */
  void appendOutputNl(String string, Style style) {
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
      SwingUtilities.invokeLater(this::scroll);
    }
  }

  /**
   * scroll
   */
  private void scroll() {
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
   * appendOutput
   *
   * @param string a {@link java.lang.String} object.
   * @param style  a {@link javax.swing.text.Style} object.
   */
  private void appendOutput(String string, Style style) {
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
      SwingUtilities.invokeLater(this::scroll);
    }
  }

  /**
   * setMeasurement
   *
   * @param measurement a {@link java.lang.String} object.
   * @param d           a double.
   */
  void setMeasurement(String measurement, double d) {
    try {
      appendOutput(measurement, getOutputStyle());
    } catch (Exception e) {
      String message = "Exception appending measurement to Shell.\n";
      logger.log(Level.WARNING, message, e);
    }
  }

  /**
   * sync
   */
  void sync() {
    try {
      setVariable("active", mainPanel.getHierarchy().getActive());
    } catch (Exception e) {
      String message = " Exception syncing shell variables.\n";
      logger.log(Level.WARNING, message, e);
    }
  }
}
