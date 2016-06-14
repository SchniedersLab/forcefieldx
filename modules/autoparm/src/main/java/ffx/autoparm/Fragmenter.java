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
package ffx.autoparm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fragment.ExhaustiveFragmenter;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.modeling.builder3d.ModelBuilder3D;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

/**
 * Splits large molecules into fragments for PolType
 *
 * @author Rae Corrigan
 */
public class Fragmenter {

    protected String sdffile;

    public Fragmenter(String sdffile) {
        this.sdffile = sdffile;
    }

    private enum FIELD {
        XLogP,
        ALogP,
        ALogp2,
        AMR,
        SMILES_Kekule,
        SMILES_Aromatic
    }

    private static final int SIZE = 30;

    @SuppressWarnings({"null", "CallToPrintStackTrace"})
    public void readSDF() throws FileNotFoundException, IOException {
        File file = new File(sdffile);
        BufferedReader read = null;

        try {
            FileReader fileReader = new FileReader(file);
            read = new BufferedReader(fileReader);
        } catch (IOException e) {
            e.printStackTrace();
        }

        IteratingSDFReader reader = null;

        try {

            reader = new IteratingSDFReader(read, SilentChemObjectBuilder.getInstance());
            while (reader.hasNext()) {

                IAtomContainer molecule = reader.next();

                try {
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                    CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(molecule);

                    //Fragmentation call
                    fragment(molecule);

                } catch (Exception x) {
                    System.err.println("*");
                }
            }
        } catch (Exception x) {
            x.printStackTrace();
        } finally {
            try {
                reader.close();
            } catch (Exception x) {
            }
        }
    }

    String[] fArray = new String[]{"For SMILES"};
    String[] frame = new String[]{"For Frameworks"};
    String[] rings = new String[]{"For Ring Systems"};
    String[] eArray = new String[]{"For Exh Fragments"};
    String[] eatArray = new String[]{"For testing for eaten fragments"};
    String[] rhArray = new String[]{"For removed-hydrogen SMILES"};
    int numSubstructures = 0;
    List<String> toRemove = new ArrayList<>();
    List<String> finalList = new ArrayList<>();
    List<IAtomContainer> rhList = new ArrayList<>();
    List<Integer> indicelist = new ArrayList<>();
    List<String> finallist = new ArrayList<>();

    @SuppressWarnings("CallToPrintStackTrace")
    protected void fragment(IAtomContainer molecule) throws Exception {
        //MurckoFragmenter implimentation
        MurckoFragmenter murk = new MurckoFragmenter();
        murk.generateFragments(molecule);
        System.out.println("MURCKO FRAGMENTS");
        fArray = murk.getFragments();

        System.out.println(Arrays.toString(fArray));

        //ExhaustiveFragmenter implimentation
        ExhaustiveFragmenter exh = new ExhaustiveFragmenter();
        exh.setMinimumFragmentSize(20);
        exh.generateFragments(molecule);
        System.out.println("\nEXHAUSTIVE FRAGMENTS");
        eArray = exh.getFragments();
        eatArray = exh.getFragments();
        rhArray = exh.getFragments();
        int orig = eatArray.length;

        System.out.println(Arrays.toString(eArray) + "\n");

        //checking for "eaten" fragments
        //remove hydrogens for more accurate substructure checking
        //for (int s = 0; s < rhArray.length; s++)
        for (String rhArray1 : rhArray) {
            //convert each array entry (SMILES) to IAtomContainer
            IAtomContainer molec = null;
            try {
                SmilesParser smp = new SmilesParser(SilentChemObjectBuilder.getInstance());
                molec = smp.parseSmiles(rhArray1);
            } catch (InvalidSmilesException ise) {
                System.out.println(ise.toString());
            }
            //remove hydrogens using CDK AtomContainerManipulator.removeHydrogens(IAtomContainer)
            try {
                AtomContainerManipulator.removeHydrogens(molec);
            } catch (Exception e) {
                e.printStackTrace();
            }
            rhList.add(molec);
        }

        //check for substructures and collect indicies to take out entires from full-H array
        for (int t = 0; t < rhList.size(); t++) {
            IAtomContainer query = rhList.get(t);
            Pattern pattern = VentoFoggia.findSubstructure(query);

            for (int u = 0; u < rhList.size(); u++) {
                IAtomContainer tester = rhList.get(u);

                //is "Query" is a substructure of "Tester"
                //makes sure query and tester aren't the same molecule and that 
                //     query is smaller than tester (substructures have to be 
                //     smaller than the main structure)
                if (pattern.matches(tester) && (tester != query) && (tester.getAtomCount() > query.getAtomCount()) && (tester.getAtomCount() < SIZE)) {
                    indicelist.add(t);
                }
            }
        }

        for (int v = 0; v < rhList.size(); v++) {
            if (!indicelist.contains(v)) {
                finallist.add(eatArray[v]);
            } else {
                numSubstructures++;
            }
        }

        //Make final list of non-substructures into an array to pass on
        String[] finalArray = new String[finallist.size()];

        for (int n = 0; n < finallist.size(); n++) {
            finalArray[n] = finallist.get(n);
        }

        System.out.print("Substructures removed: ");
        System.out.println(numSubstructures);

        System.out.println("Orig length: " + orig);
        System.out.println("Final length: " + finalArray.length + "\n");

        //write SMILES file
        //pass fArray to smilesToObject to convert Murcko fragments to SDF
        //pass eArray to smilesToObject to convert Exhaustive fragments to SDF
        //pass finalArray to smilesToObject to convert non-substructure fragments to SDF
        smilesToObject(finalArray);

    }

    protected void smilesToObject(String[] smilesArr) throws Exception {

        for (int i = 0; i < smilesArr.length; i++) {
            String content = smilesArr[i];
            //entry number in SMILES array to be used later in SDF writing
            int num = i;
            //System.out.println("\nPassing SMILES string to doConversion\n");
            doConversion(content, num);
        }

    }

    protected void doConversion(String smi, int num) throws Exception {
        IAtomContainer mol = null;
        //entry number in SMILES array to be used later in SDF writing
        int number = num;

        try {
            SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
            mol = sp.parseSmiles(smi);
        } catch (InvalidSmilesException ise) {
            System.out.println(ise.toString());
        }

        try {
            @SuppressWarnings("null")
            CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(mol.getBuilder());
            for (IAtom atom : mol.atoms()) {
                IAtomType type = matcher.findMatchingAtomType(mol, atom);
                AtomTypeManipulator.configure(atom, type);
            }
        } catch (Exception e) {
            System.out.println(e);
        }

        try {
            @SuppressWarnings("null")
            CDKHydrogenAdder ha = CDKHydrogenAdder.getInstance(mol.getBuilder());
            ha.addImplicitHydrogens(mol);
            AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
        } catch (CDKException e) {
            System.out.println(e);
        }

        ModelBuilder3D mb3d;
        mb3d = ModelBuilder3D.getInstance(SilentChemObjectBuilder.getInstance());
        @SuppressWarnings("UnusedAssignment")
        IAtomContainer molecule = null;
        molecule = mb3d.generate3DCoordinates(mol, false);

        //Fragmenting checks
        //30 atoms or less
        if (molecule.getAtomCount() < SIZE) {

            //"eaten" fragments checked for already
            //writeSDF
            writeSDF(molecule, number);
        }

    }

    @SuppressWarnings("CallToPrintStackTrace")
    protected void writeSDF(IAtomContainer sm, int n) throws Exception {
        //convert SMILES to sdf's
        //entry number in SMILES array, to aid in file naming
        int j = n;
        String fileBegin = "fragment.sdf_";
        String fileEnd = Integer.toString(j);
        String filename = fileBegin.concat(fileEnd);
        System.out.print("filename: ");
        System.out.println(filename);
        System.out.println();

        SDFWriter swrite = new SDFWriter();

        try {
            File sdfFromSMILES = new File(filename);
            FileWriter filew = new FileWriter(sdfFromSMILES.getAbsoluteFile());
            BufferedWriter out = new BufferedWriter(filew);

            swrite.setWriter(out);
            swrite.write(sm);
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            swrite.close();
        }

    }

}


/*SMIlES pattern matching substructure check

int contains = 0;
        int len = 0;

        //checking "eaten" fragments
        System.out.println("\nChecking for substructures\n");
        for (int l = 0; l < eArray.length; l++) {
            String str1 = eatArray[l];

            for (int m = 0; m < (eArray.length); m++) {
                String str2 = eatArray[m];
                if (str2.contains(str1)) {
                    contains = 1;
                } else {
                    contains = 0;
                }
                if (str2.length() != str1.length()) {
                    len = 1;
                } else {
                    len = 0;
                }

                if (str2.contains(str1) && (str2.length() != str1.length())) {
                    toRemove.add(str1);
                    numSubstructures++;
                }
            } 
        }//System.out.println("Done looking for substructures\n");

        System.out.println("toRemove: ");
        System.out.println(Arrays.toString(toRemove.toArray()));

        int alreadyThere = 0;
        int keepCount = 0;
        System.out.println("Creating final array\n");
        for (int f = 0; f < eArray.length; f++) {
            String test = eatArray[f];
            alreadyThere = 0;
            keepCount = 0;
            for (int g = 0; g < toRemove.size(); g++) {
                if ((test != toRemove.get(g))) {
                    keepCount++;
                    alreadyThere = 1;
                } else if (test == toRemove.get(g)) {
                    //System.out.println("Found frag to leave out");
                }
            }

            if (keepCount == toRemove.size()) {
                finalList.add(test);
            }
        }*/
