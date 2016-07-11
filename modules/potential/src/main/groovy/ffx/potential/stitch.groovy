
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
//if (options.h || arguments == null || arguments.size() != 1){ original code
if (options.h || arguments == null) {
    return cli.usage();
} 

//TODO: Generalize Stitch 
List<MolecularAssembly> molecularAssemblies = new ArrayList<>();
// Read in command line.
//General
String mapname = arguments.get(0);
//while there is a "next argument", make it a string and name it "patch#"
for(int i = 1; i < (arguments.size() - 1 ); i++){
    String strBegin = "patch";
    String strEnd = Integer.toString(i);
    String strName = strBegin.concat(strEnd);
    
    strName = arguments.get(i);
    
    System.out.println("Args name test: ");
    System.out.println(arguments.getAt(i));
    
    //open MolecularAssemblies for each pdb; named the same as the patches
    //name PDBs to open
    String pdbname = arguments.getAt(i);
    pdbname = pdbname.replaceAll("patch","pdb");
    String patchname = arguments.getAt(i);
    patchname = patchname.replaceAll(".patch", "\0");
    String ffname = patchname.concat("FF");
    String ffename = patchname.concat("energy");
    
    System.out.println("PDB name: "+pdbname);
    System.out.println("Patchname: "+patchname);
    System.out.println("ForceField name: "+ffname);
    System.out.println("ForceField Energy name: "+ffename+"\n");
    
    //open Molecular assembiles accordingly
    open(pdbname);
    //MolecularAssembly [patchname] = (MolecularAssembly) active;
    //molecularAssemblies.add([patchname]);
    
    //Forcefields and forcefield energies
    //ForceField [ffname] = [patchname].getForceField();
    //ForceFieldEnergy [ffename] = [patchname].getPotentialEnergy();

    // Loop over force field terms from ForceFieldEnergy instances

    // Bonds

    //Bond bonds[] = [ffename].getBonds();
}

System.out.println();

//Specific to Vemurafenib
/*String mapname = arguments.get(0);
String patch1 = arguments.get(1);
String patch2 = arguments.get(2);
String patch3 = arguments.get(3);
String patch4 = arguments.get(4);
String patch5 = arguments.get(5);
String patch6 = arguments.get(6);
//String fullPatch = arguments.get(7);*/

/*HashMap<String,String> myMap = new HashMap<>();
myMap.put("StringKey", "StringValue");
if (!myMap.get("StringKey").equals("StringValue")) {
System.out.println("Error!");
}*/

//Specific to Vemurafenib
/*open("CMB.pdb");
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
ForceField agnFF = AGN.getForceField();*/

//Specific to Vemurafenib
/*molecularAssemblies.add(CMB);
molecularAssemblies.add(CPP);
molecularAssemblies.add(FBN);
molecularAssemblies.add(NDF);
molecularAssemblies.add(POS);
molecularAssemblies.add(AGN);*/



/*PatchCombiner pc = new PatchCombiner(molecularAssemblies, mapname, patch1, patch2, patch3, patch4, patch5, patch6);
pc.myMethod();
 */
return;

