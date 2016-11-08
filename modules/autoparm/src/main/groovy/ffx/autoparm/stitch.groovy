
/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

// STITCH
// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.util.CliBuilder;
import groovy.lang.MissingPropertyException;

// FFX Imports
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.PatchCombiner
import ffx.potential.MolecularAssembly;

//Java Imports
import java.io.FileReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.Object;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.io.CharArrayReader;
import static java.lang.System.out;
import ffx.potential.ForceFieldEnergy
import java.lang.ArrayIndexOutOfBoundsException;
import java.lang.IndexOutOfBoundsException;
import java.lang.String;
import java.util.regex.Pattern;
import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * PolType parameter stitching code
 * 
 * @author Rae Ann Corrigan
 */

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc stitch [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();

if (options.h || arguments == null) {
    return cli.usage();
} 

List<File> textFiles = new ArrayList<>();
List<File> fragmentParameters = new ArrayList<>();
// Read in command line.
//while there is a "next argument"
for(int i = 0; i < (arguments.size()); i++){
    if(arguments.get(i).contains(".txt")){
        //add to "textFiles" array list
        textFiles.add(arguments.get(i));
    }
    if(arguments.get(i).contains(".key_5")){
        //add to "fragmentParameters" array list (PolType outputs)
        fragmentParameters.add(arguments.get(i));
    }
    
    
}//end read command line dynamically "for"

Stitch(fragmentParameters, textFiles);

return;

