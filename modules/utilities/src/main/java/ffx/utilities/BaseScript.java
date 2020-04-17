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
package ffx.utilities;

import static java.lang.String.format;
import static java.util.Collections.sort;
import static picocli.CommandLine.usage;

import groovy.lang.Binding;
import groovy.lang.Script;
import java.awt.GraphicsEnvironment;
import java.io.ByteArrayOutputStream;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.ZipEntry;
import picocli.CommandLine;
import picocli.CommandLine.Help.Ansi;
import picocli.CommandLine.Option;
import picocli.CommandLine.ParseResult;

/**
 * BaseScript class.
 *
 * @author Michael J. Schnieders
 */
public abstract class BaseScript extends Script {

  /** The logger for this class. */
  protected static final Logger logger = Logger.getLogger(BaseScript.class.getName());

  /**
   * Unix shells are able to evaluate PicoCLI ANSI color codes, but right now the FFX GUI Shell does
   * not.
   *
   * <p>In a headless environment, color will be ON for command line help, but OFF for the GUI.
   */
  public final Ansi color;

  /** The Groovy Binding contains defined variables, closures, etc. */
  public Binding context;

  /** The array of args passed into the Script. */
  public String[] args;

  /** Parse Result. */
  public ParseResult parseResult = null;

  /** -V or --version Prints the FFX version and exits. */
  @Option(
      names = {"-V", "--version"},
      versionHelp = true,
      description = "Print the Force Field X version and exit.")
  public boolean version = false;

  /** -h or --help Prints a help message. */
  @Option(
      names = {"-h", "--help"},
      usageHelp = true,
      description = "Print command help and exit.")
  public boolean help = false;

  /** Default constructor for an FFX Script. */
  public BaseScript() {
    if (GraphicsEnvironment.isHeadless()) {
      color = Ansi.ON;
    } else {
      color = Ansi.OFF;
    }
  }

  /**
   * Use the System ClassLoader to find the requested script.
   *
   * @param name Name of the script to load (e.g. Energy).
   * @return The Script, if found, or null.
   */
  public static Class<? extends BaseScript> getScript(String name) {
    ClassLoader loader = ClassLoader.getSystemClassLoader();
    String pathName = name;
    Class<?> script;
    try {
      // First try to load the class directly.
      script = loader.loadClass(pathName);
    } catch (ClassNotFoundException e) {
      // Next, try to load a script from the potential package.
      pathName = "ffx.potential.groovy." + name;
      try {
        script = loader.loadClass(pathName);
      } catch (ClassNotFoundException e2) {
        // Next, try to load a script from the algorithm package.
        pathName = "ffx.algorithms.groovy." + name;
        try {
          script = loader.loadClass(pathName);
        } catch (ClassNotFoundException e3) {
          // Finally, try to load a script from the xray package.
          pathName = "ffx.xray.groovy." + name;
          try {
            script = loader.loadClass(pathName);
          } catch (ClassNotFoundException e4) {
            logger.warning(format(" %s was not found.", name));
            return null;
          }
        }
      }
    }
    return script.asSubclass(BaseScript.class);
  }

  /**
   * List the embedded FFX Groovy Scripts.
   *
   * @param logScripts List Scripts.
   * @param logTestScripts List Test Scripts.
   */
  public static void listGroovyScripts(boolean logScripts, boolean logTestScripts) {

    ClassLoader classLoader = ClassLoader.getSystemClassLoader();
    String location = "ffx";
    URL scriptURL = classLoader.getResource(location);
    if (scriptURL == null) {
      logger.info(format(" The %s resource could not be found by the classloader.", location));
      return;
    }
    String scriptPath = scriptURL.getPath();
    String ffx = scriptPath.substring(5, scriptURL.getPath().indexOf("!"));
    JarFile jar;
    try {
      jar = new JarFile(URLDecoder.decode(ffx, StandardCharsets.UTF_8));
    } catch (Exception e) {
      logger.info(format(" The %s resource could not be decoded.", location));
      return;
    }

    List<String> scripts = new ArrayList<>();
    List<String> testScripts = new ArrayList<>();
    // Iterates over Jar entries.
    Enumeration<JarEntry> enumeration = jar.entries();
    while (enumeration.hasMoreElements()) {
      ZipEntry zipEntry = enumeration.nextElement();
      String className = zipEntry.getName();
      if (className.startsWith("ffx")
          && className.contains("groovy")
          && className.endsWith(".class")
          && !className.contains("$")) {
        className = className.replace("/", ".");
        className = className.replace(".class", "");
        // Present the classes using "short-cut" names.
        className = className.replace("ffx.potential.groovy.", "");
        className = className.replace("ffx.algorithms.groovy.", "");
        className = className.replace("ffx.realspace.groovy", "realspace");
        className = className.replace("ffx.xray.groovy", "xray");
        if (className.toUpperCase().contains("TEST")) {
          testScripts.add(className);
        } else {
          scripts.add(className);
        }
      }
    }

    // Sort the scripts alphabetically.
    sort(scripts);
    sort(testScripts);

    // Log the script names.
    if (logTestScripts) {
      for (String script : testScripts) {
        logger.info("   " + script);
      }
    }
    if (logScripts) {
      for (String script : scripts) {
        logger.info("   " + script);
      }
    }
  }

  /**
   * Default help information.
   *
   * @return String describing how to use this command.
   */
  public String helpString() {
    try {
      StringOutputStream sos = new StringOutputStream(new ByteArrayOutputStream());
      usage(this, sos, color);
      return " " + sos.toString();
    } catch (UnsupportedEncodingException e) {
      logger.log(Level.WARNING, e.toString());
      return null;
    }
  }

  /**
   * Initialize this Script based on the specified command line arguments.
   *
   * @return boolean Returns true if the script should continue and false to exit.
   */
  public boolean init() {
    context = getBinding();

    // The args property could either be a list or an array of String arguments.
    Object arguments = context.getProperty("args");
    if (arguments instanceof List<?>) {
      List<?> list = (List<?>) arguments;
      int numArgs = list.size();
      args = new String[numArgs];
      for (int i = 0; i < numArgs; i++) {
        args[i] = (String) list.get(i);
      }
    } else {
      args = (String[]) context.getProperty("args");
    }

    CommandLine commandLine = new CommandLine(this);
    try {
      parseResult = commandLine.parseArgs(args);
    } catch (CommandLine.UnmatchedArgumentException uae) {
      logger.warning(
          " The usual source of this exception is when long-form arguments (such as --uaA) are only preceded by one dash (such as -uaA, which is an error).");
      throw uae;
    }

    // Print help info exit.
    if (help) {
      logger.info(helpString());
      return false;
    }

    // Version info is printed by default.
    if (version) {
      // This should not be reached, due to the FFX Main class handling the "-V, --version" flag and
      // exiting.
      return false;
    }

    return true;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Execute the script.
   */
  @Override
  public BaseScript run() {
    logger.info(helpString());
    return this;
  }
}
