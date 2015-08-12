
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


// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc stitch [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
//if (options.h || arguments == null || arguments.size() != 1){ original code
if (options.h || arguments == null) {
    return cli.usage();
} 
// Read in command line.
String mapname = arguments.get(0);
String patch1 = arguments.get(1);
String patch2 = arguments.get(2);
String patch3 = arguments.get(3);
String patch4 = arguments.get(4);
String patch5 = arguments.get(5);
String patch6 = arguments.get(6);
//String fullPatch = arguments.get(7);

//systems = open(xyzname);
//energy();

//define file opener class, then read map file, and print to screen
//remember to build mvn in ffx dir. in a separate terminal window, not the build in Netbeans, before running stitch
/*try {
    File file = new File(mapname);
    FileReader fileReader = new FileReader(file);
    BufferedReader bufferedReader = new BufferedReader(fileReader);
    StringBuffer stringBuffer = new StringBuffer();
    String line;
    while ((line = bufferedReader.readLine()) != null) {
        String[] parts = (line.split(Pattern.quote(" : ")));
        String AtomName = parts[0];
        String FragCode = parts[1];
        String fAtomName = parts[2];
        String TypeNum = parts[3];
    }
    fileReader.close();
    //System.out.println(Arrays.toString(parts));
} catch (IOException e) {
    e.printStackTrace();
}*/



HashMap<String,String> myMap = new HashMap<>();
myMap.put("StringKey", "StringValue");
if (!myMap.get("StringKey").equals("StringValue")) {
    System.out.println("Error!");
}

/*HashMap<Integer,Double> vdwMap = new HashMap<>();
vdwMap.put(509, 2.58);

HashMap<List<Integer>,Integer> typeMap = new HashMap<>();
List<Integer> chlorides = new ArrayList<>();
chlorides.add(409);
chlorides.add(411);
typeMap.put(chlorides, 509);*/

open("CMB.pdb");
MolecularAssembly CMB = (MolecularAssembly) active;
open("CPP.pdb");
MolecularAssembly CPP = (MolecularAssembly) active;
open("FBN.pdb");
MolecularAssembly FBN = (MolecularAssembly) active;
open("NDF.pdb");
MolecularAssembly NDF = (MolecularAssembly) active;
open("POS.pdb");
MolecularAssembly POS = (MolecularAssembly) active;
open("AGN.pdb");
MolecularAssembly AGN = (MolecularAssembly) active;

ForceField cmbFF = CMB.getForceField();
ForceField cppFF = CPP.getForceField();
ForceField fbnFF = FBN.getForceField();
ForceField ndfFF = NDF.getForceField();
ForceField posFF = POS.getForceField();
ForceField agnFF = AGN.getForceField();

List<MolecularAssembly> molecularAssemblies = new ArrayList<>();
molecularAssemblies.add(CMB);
molecularAssemblies.add(CPP);
molecularAssemblies.add(FBN);
molecularAssemblies.add(NDF);
molecularAssemblies.add(POS);
molecularAssemblies.add(AGN);


PatchCombiner pc = new PatchCombiner(molecularAssemblies, mapname, patch1, patch2, patch3, patch4, patch5, patch6);
pc.myMethod();

return;

