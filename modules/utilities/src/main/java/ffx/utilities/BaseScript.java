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
package ffx.utilities;

import java.awt.GraphicsEnvironment;
import java.io.ByteArrayOutputStream;
import java.io.UnsupportedEncodingException;
import java.util.logging.Level;
import java.util.logging.Logger;

import groovy.lang.Binding;
import groovy.lang.Script;
import picocli.CommandLine;
import picocli.CommandLine.Help.Ansi;
import picocli.CommandLine.Option;
import picocli.CommandLine.ParseResult;

/**
 * <p>BaseScript class.</p>
 *
 * @author Michael J. Schnieders
 */
public class BaseScript extends Script {

    /**
     * The logger for this class.
     */
    public static final Logger logger = Logger.getLogger(BaseScript.class.getName());

    /**
     * Unix shells are able to evaluate PicoCLI ANSI color codes, but right now the FFX GUI Shell does not.
     * <p>
     * In a headless environment, color will be ON for command line help, but OFF for the GUI.
     */
    public final Ansi color;

    /**
     * The Groovy Binding contains defined variables, closures, etc.
     */
    public Binding context;

    /**
     * The array of args passed into the Script.
     */
    public String[] args;

    /**
     * Parse Result.
     */
    public ParseResult parseResult = null;

    /**
     * -h or --help to print a help message
     */
    @Option(names = {"-h", "--help"}, usageHelp = true, description = "Print this help message.")
    public boolean help = false;

    /**
     * Default constructor for an FFX Script.
     */
    public BaseScript() {
        if (GraphicsEnvironment.isHeadless()) {
            color = Ansi.ON;
        } else {
            color = Ansi.OFF;
        }
    }

    /**
     * Initialize this Script based on the specified command line arguments.
     *
     * @return boolean Returns true if the script should continue and false to exit.
     */
    public boolean init() {
        context = getBinding();
        args = (String[]) context.getProperty("args");

        CommandLine commandLine = new CommandLine(this);
        try {
            parseResult = commandLine.parseArgs(args);
        } catch (CommandLine.UnmatchedArgumentException uae) {
            logger.warning(" The usual source of this exception is when long-form arguments (such as --uaA) are only preceded by one dash (such as -uaA, which is an error).");
            throw uae;
        }

        if (help) {
            logger.info(helpString());
            return false;
        }

        return true;
    }

    /**
     * Default help information.
     *
     * @return String describing how to use this command.
     */
    public String helpString() {
        try {
            StringOutputStream sos = new StringOutputStream(new ByteArrayOutputStream());
            picocli.CommandLine.usage(this, sos, color);
            return " " + sos.toString();
        } catch (UnsupportedEncodingException e) {
            logger.log(Level.WARNING, e.toString());
            return null;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Execute the script.
     */
    @Override
    public BaseScript run() {
        logger.info(helpString());
        return this;
    }

}
