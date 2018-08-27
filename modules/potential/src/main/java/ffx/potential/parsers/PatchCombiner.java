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
package ffx.potential.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import ffx.potential.MolecularAssembly;

/**
 * Patch Combiner merges multiple patch files.
 *
 * @author Rae Ann Corrigan
 *
 */
public class PatchCombiner {

    private static final Logger logger = Logger.getLogger(ForceFieldFilter.class.getName());
    private final List<MolecularAssembly> molecularAssemblies;
    private final String mapname;
    private final String patches[];

    List<String> atomNameList;
    List<String> fragCodeList;
    List<String> fAtomNameList;
    List<String> typeNumList;
    List<String> vTypeNumList;

    String[] atomName;
    String[] fragCode;
    String[] fAtomName;
    String[] typeNum;
    String[] vTypeNum;
    int max;

    HashMap<String, String> typeNumMap;
    HashMap<Integer, String> patchMap;
    HashMap<String, String> Map4to51;
    HashMap<String, String> Map4to52;
    HashMap<String, String> Map4to53;
    HashMap<String, String> Map4to54;
    HashMap<String, String> Map4to55;
    HashMap<String, String> Map4to56;

    /**
     * <p>Constructor for PatchCombiner.</p>
     *
     * @param molecularAssemblies a {@link java.util.List} object.
     * @param mapname a {@link java.lang.String} object.
     * @param patch a {@link java.lang.String} object.
     */
    public PatchCombiner(List<MolecularAssembly> molecularAssemblies, String mapname,
                         String... patch) {
        this.molecularAssemblies = molecularAssemblies;
        this.mapname = mapname;
        this.patches = patch;
    }

    /**
     * <p>combinePatches.</p>
     *
     * @throws java.io.FileNotFoundException if any.
     */
    public void combinePatches() throws FileNotFoundException {

        atomNameList = new ArrayList<>();
        fragCodeList = new ArrayList<>();
        fAtomNameList = new ArrayList<>();
        typeNumList = new ArrayList<>();
        vTypeNumList = new ArrayList<>();
        typeNumMap = new HashMap<>();

        try {
            File file = new File(mapname);
            FileReader fileReader = new FileReader(file);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String line;

            while ((line = bufferedReader.readLine()) != null) {
                String[] pieces = line.split(" : ");
                if (pieces.length == 5) {
                    atomNameList.add(pieces[0]);
                    fragCodeList.add(pieces[1]);
                    fAtomNameList.add(pieces[2]);
                    typeNumList.add(pieces[3]);
                    vTypeNumList.add(pieces[4]);
                    typeNumMap.put(pieces[3], pieces[4]);
                } else {
                    System.out.println("Skipped Line: " + line);
                }

            }

            fileReader.close();
            atomNameList.remove("Atom Name");
            fragCodeList.remove("Fragment Code");
            fAtomNameList.remove("fAtom Name");
            typeNumList.remove("fAtom Type Number");
            vTypeNumList.remove("vAtom Type Number");
        } catch (IOException e) {
        }

        patchMap = new HashMap<>();
        patchMap.put(1, patches[0]);
        patchMap.put(2, patches[1]);
        patchMap.put(3, patches[2]);
        patchMap.put(4, patches[3]);
        patchMap.put(5, patches[4]);
        patchMap.put(6, patches[5]);

        atomName = new String[atomNameList.size()];
        atomNameList.toArray(atomName);
        fragCode = new String[fragCodeList.size()];
        fragCodeList.toArray(fragCode);
        fAtomName = new String[fAtomNameList.size()];
        fAtomNameList.toArray(fAtomName);
        typeNum = new String[typeNumList.size()];
        typeNumList.toArray(typeNum);
        vTypeNum = new String[vTypeNumList.size()];
        vTypeNumList.toArray(vTypeNum);

        Map4to51 = new HashMap<>();
        Map4to52 = new HashMap<>();
        Map4to53 = new HashMap<>();
        Map4to54 = new HashMap<>();
        Map4to55 = new HashMap<>();
        Map4to56 = new HashMap<>();

        for (int l = 0; l < fragCode.length; ++l) {
            if (Integer.parseInt(fragCode[l]) == 1) {
                Map4to51.put(typeNum[l], vTypeNum[l]);
            } else if (Integer.parseInt(fragCode[l]) == 2) {
                Map4to52.put(typeNum[l], vTypeNum[l]);
            } else if (Integer.parseInt(fragCode[l]) == 3) {
                Map4to53.put(typeNum[l], vTypeNum[l]);
            } else if (Integer.parseInt(fragCode[l]) == 4) {
                Map4to54.put(typeNum[l], vTypeNum[l]);
            } else if (Integer.parseInt(fragCode[l]) == 5) {
                Map4to55.put(typeNum[l], vTypeNum[l]);
            } else if (Integer.parseInt(fragCode[l]) == 6) {
                Map4to56.put(typeNum[l], vTypeNum[l]);
            }
        }

        max = atomName.length;

        /**
         * Atom writing code
         */
        writeAtoms();

        /**
         * *
         * Multipole averaging/writing code
         */
        writeMultipoles();

        /**
         * Polarize averaging/writing code
         */
        writePolarize();

        /**
         * VDW averaging/writing code
         */
        writeVDW();

        /**
         * Bond averaging/writing code
         */
    }

    /**
     * <p>writeAtoms.</p>
     *
     * @throws java.io.FileNotFoundException if any.
     */
    public void writeAtoms() throws FileNotFoundException {
        PrintWriter atomWriter1 = new PrintWriter("Oatom.patch");

        for (int count = 1; count < 7; ++count) {

            String code = Integer.toString(count);
            String patchname = patchMap.get(Integer.parseInt(code));

            try {
                File file = new File(patchname);
                FileReader atomfr = new FileReader(file);
                BufferedReader atombr = new BufferedReader(atomfr);
                StringBuilder stringBuilder = new StringBuilder();
                String line;
                while ((line = atombr.readLine()) != null) {
                    stringBuilder.append(line);
                    stringBuilder.append("\n");

                    String patch = line;

                    //find atom values
                    if (patch.contains("atom")) {

                        String[] atomSplit1 = line.split(" ");

                        //remove spaces
                        List<String> list = new ArrayList<>(Arrays.asList(atomSplit1));
                        list.removeAll(Arrays.asList("", null));

                        String[] atomSplit2 = new String[list.size()];
                        list.toArray(atomSplit2);

                        //change 4** typeNum to 5** vTypeNum for entry into patch file
                        if (Integer.parseInt(code) == 1) {
                            atomSplit2[1] = Map4to51.get(atomSplit2[1]);
                            atomSplit2[2] = Map4to51.get(atomSplit2[2]);
                        } else if (Integer.parseInt(code) == 2) {
                            atomSplit2[1] = Map4to52.get(atomSplit2[1]);
                            atomSplit2[2] = Map4to52.get(atomSplit2[2]);
                        } else if (Integer.parseInt(code) == 3) {
                            atomSplit2[1] = Map4to53.get(atomSplit2[1]);
                            atomSplit2[2] = Map4to53.get(atomSplit2[2]);
                        } else if (Integer.parseInt(code) == 4) {
                            atomSplit2[1] = Map4to54.get(atomSplit2[1]);
                            atomSplit2[2] = Map4to54.get(atomSplit2[2]);
                        } else if (Integer.parseInt(code) == 5) {
                            atomSplit2[1] = Map4to55.get(atomSplit2[1]);
                            atomSplit2[2] = Map4to55.get(atomSplit2[2]);
                        } else if (Integer.parseInt(code) == 6) {
                            atomSplit2[1] = Map4to56.get(atomSplit2[1]);
                            atomSplit2[2] = Map4to56.get(atomSplit2[2]);
                        } else {
                            System.out.println("Not in a patch");
                        }

                        if (atomSplit2[1] != null) {
                            atomWriter1.write("atom        " + atomSplit2[1] + "   " + atomSplit2[2] + "    " + atomSplit2[3] + "     " + "'Vemurafenib'" + "         " + atomSplit2[5] + "    " + atomSplit2[6] + "    " + atomSplit2[7]);
                            atomWriter1.write("\n");
                        }

                    }

                }
            } catch (IOException e) {
            }
        }
        atomWriter1.close();

        //remove duplicate lines
        String result;
        StringBuilder resultBuilder = new StringBuilder();
        Set<String> alreadyPresent1 = new HashSet<>();

        try {
            FileReader fr = new FileReader("Oatom.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            String line1;
            while ((line1 = bufferedReader.readLine()) != null) {

                String[] aParts = line1.split(" ");

                boolean first = true;
                if (!alreadyPresent1.contains(aParts[8])) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilder.append("\n");
                    }

                    if (!alreadyPresent1.contains(aParts[8])) {
                        resultBuilder.append(line1);
                        resultBuilder.append("\n");
                    }

                    alreadyPresent1.add(aParts[8]);

                }
                result = resultBuilder.toString();

                try {

                    PrintWriter atomWriter2 = new PrintWriter("atom.patch");
                    atomWriter2.write(result);
                    atomWriter2.write("\n");
                    atomWriter2.close();

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Atom Patch: COMPLETE");
    }

    /**
     * <p>writeMultipoles.</p>
     *
     * @throws java.io.FileNotFoundException if any.
     */
    public void writeMultipoles() throws FileNotFoundException {
        List<String> multipoles = new ArrayList<>();
        List<String> mCol1List = new ArrayList<>();
        List<String> mCol2List = new ArrayList<>();
        List<String> mCol3List = new ArrayList<>();
        List<String> mValue1List = new ArrayList<>();

        String mpLine1 = "                                        0.";
        String mpLine2 = "                                       -0.";
        String mpLine3 = "                                       -1.";
        String mpLine4 = "                                        1.";
        int lineCount = 0;

        List<Double> l2c1List = new ArrayList<>();
        List<Double> l2c2List = new ArrayList<>();
        List<Double> l2c3List = new ArrayList<>();
        List<Double> l3c1List = new ArrayList<>();
        List<Double> l4c1List = new ArrayList<>();
        List<Double> l4c2List = new ArrayList<>();
        List<Double> l5c1List = new ArrayList<>();
        List<Double> l5c2List = new ArrayList<>();
        List<Double> l5c3List = new ArrayList<>();

        for (int j = 0; j < max; ++j) {

            String code = fragCode[j];
            String patchname = patchMap.get(Integer.parseInt(code));
            String multipole = null;

            try {
                File file = new File(patchname);
                FileReader fileReader = new FileReader(file);
                BufferedReader bufferedReader = new BufferedReader(fileReader);
                StringBuilder stringBuffer = new StringBuilder();
                String line;

                while ((line = bufferedReader.readLine()) != null) {
                    stringBuffer.append(line);
                    stringBuffer.append("\n");

                    String patch = line;

                    if (patch.contains("multipole")) {
                        multipole = line;
                        lineCount = 0;

                        String[] multipole1 = multipole.split(" ");
                        String[] multipole2 = new String[multipole1.length];

                        for (int k = 0; k < multipole1.length; ++k) {
                            multipole2[k] = multipole1[k];
                        }

                        //need to eliminate negative signs from mapping and put them back after mapping
                        String[] chars = multipole.split("(?!^)");
                        String first = (chars[12] + chars[13] + chars[14]);
                        String second = (chars[17] + chars[18] + chars[19]);
                        String third = (chars[22] + chars[23] + chars[24]);

                        //replaces fragment atom type numbers with vem. atom type numbers
                        if (Integer.parseInt(code) == 1) {

                            multipole2[3] = Map4to51.get(first);
                            multipole2[5] = Map4to51.get(second);
                            multipole2[7] = Map4to51.get(third);

                        } else if (Integer.parseInt(code) == 2) {

                            multipole2[3] = Map4to52.get(first);
                            multipole2[5] = Map4to52.get(second);
                            multipole2[7] = Map4to52.get(third);

                        } else if (Integer.parseInt(code) == 3) {

                            multipole2[3] = Map4to53.get(first);
                            multipole2[5] = Map4to53.get(second);
                            multipole2[7] = Map4to53.get(third);

                        } else if (Integer.parseInt(code) == 4) {

                            multipole2[3] = Map4to54.get(first);
                            multipole2[5] = Map4to54.get(second);
                            multipole2[7] = Map4to54.get(third);

                        } else if (Integer.parseInt(code) == 5) {

                            multipole2[3] = Map4to55.get(first);
                            multipole2[5] = Map4to55.get(second);
                            multipole2[7] = Map4to55.get(third);

                        } else if (Integer.parseInt(code) == 6) {

                            multipole2[3] = Map4to56.get(first);
                            multipole2[5] = Map4to56.get(second);
                            multipole2[7] = Map4to56.get(third);

                        } else {
                            System.out.println("Not in a patch");
                        }

                        //makes sure negative signs stay with negative type numbers
                        if ("-".equals(chars[11])) {
                            multipole2[3] = chars[11].concat(multipole2[3]);
                        }
                        if ("-".equals(chars[16])) {
                            multipole2[5] = chars[16].concat(multipole2[5]);
                        }
                        if ("-".equals(chars[21])) {
                            multipole2[7] = chars[21].concat(multipole2[7]);
                        }

                        mCol1List.add(multipole2[3]);
                        mCol2List.add(multipole2[5]);
                        mCol3List.add(multipole2[7]);
                        mValue1List.add(multipole2[multipole2.length - 1]);

                        //adds to multipoles list
                        multipoles.add("multipole   " + multipole2[3] + " " + multipole2[5] + " " + multipole2[7] + "              " + multipole2[multipole2.length - 1]);

                    } //ends if(patch.contains("multipole"))
                    //enters else if loop for lines 2-5 that don't begin with the word multipole
                    else if (patch.contains(mpLine1) || patch.contains(mpLine2) || patch.contains(mpLine3) || patch.contains(mpLine4)) {
                        if (lineCount == 1) {
                            String line2 = line;
                            String[] line2a = line2.split(" ");

                            //remove spaces
                            List<String> list = new ArrayList<>(Arrays.asList(line2a));
                            list.removeAll(Arrays.asList("", null));
                            String[] line2b = new String[list.size()];
                            list.toArray(line2b);

                            l2c1List.add(Double.parseDouble(line2b[0]));
                            l2c2List.add(Double.parseDouble(line2b[1]));
                            l2c3List.add(Double.parseDouble(line2b[2]));

                        } else if (lineCount == 2) {
                            String line3 = line;
                            String[] line3a = line3.split(" ");

                            //remove spaces
                            List<String> list = new ArrayList<>(Arrays.asList(line3a));
                            list.removeAll(Arrays.asList("", null));
                            String[] line3b = new String[list.size()];
                            list.toArray(line3b);

                            l3c1List.add(Double.parseDouble(line3b[0]));

                        } else if (lineCount == 3) {
                            String line4 = line;
                            String[] line4a = line4.split(" ");

                            //remove spaces
                            List<String> list = new ArrayList<>(Arrays.asList(line4a));
                            list.removeAll(Arrays.asList("", null));
                            String[] line4b = new String[list.size()];
                            list.toArray(line4b);

                            l4c1List.add(Double.parseDouble(line4b[0]));
                            l4c2List.add(Double.parseDouble(line4b[1]));

                        } else if (lineCount == 4) {
                            String line5 = line;
                            String[] line5a = line5.split(" ");

                            //remove spaces
                            List<String> list = new ArrayList<>(Arrays.asList(line5a));
                            list.removeAll(Arrays.asList("", null));
                            String[] line5b = new String[list.size()];
                            list.toArray(line5b);

                            l5c1List.add(Double.parseDouble(line5b[0]));
                            l5c2List.add(Double.parseDouble(line5b[1]));
                            l5c3List.add(Double.parseDouble(line5b[2]));

                        }
                    }

                    ++lineCount;
                }//ends reader while
            }//end try
            catch (IOException e) {

            }
        }//ends j loop

        //define new arrays to enter list data into
        String[] mCol1 = new String[mCol1List.size()];
        String[] mCol2 = new String[mCol2List.size()];
        String[] mCol3 = new String[mCol3List.size()];
        String[] mValue1 = new String[mValue1List.size()];
        Double[] l2c1 = new Double[l2c1List.size()];
        Double[] l2c2 = new Double[l2c2List.size()];
        Double[] l2c3 = new Double[l2c3List.size()];
        Double[] l3c1 = new Double[l3c1List.size()];
        Double[] l4c1 = new Double[l4c1List.size()];
        Double[] l4c2 = new Double[l4c2List.size()];
        Double[] l5c1 = new Double[l5c1List.size()];
        Double[] l5c2 = new Double[l5c2List.size()];
        Double[] l5c3 = new Double[l5c3List.size()];

        //convert string lists to string arrays for further use
        mCol1List.toArray(mCol1);
        mCol2List.toArray(mCol2);
        mCol3List.toArray(mCol3);
        mValue1List.toArray(mValue1);

        l2c1List.toArray(l2c1);
        l2c2List.toArray(l2c2);
        l2c3List.toArray(l2c3);
        l3c1List.toArray(l3c1);
        l4c1List.toArray(l4c1);
        l4c2List.toArray(l4c2);
        l5c1List.toArray(l5c1);
        l5c2List.toArray(l5c2);
        l5c3List.toArray(l5c3);

        PrintWriter multipoleWriter1 = new PrintWriter("Omultipole.patch");

        List<String> finall2c1L = new ArrayList<>();
        List<String> finall2c2L = new ArrayList<>();
        List<String> finall2c3L = new ArrayList<>();
        List<String> finall3c1L = new ArrayList<>();
        List<String> finall4c1L = new ArrayList<>();
        List<String> finall4c2L = new ArrayList<>();
        List<String> finall5c1L = new ArrayList<>();
        List<String> finall5c2L = new ArrayList<>();
        List<String> finall5c3L = new ArrayList<>();

        for (int i = 0; i < multipoles.size(); ++i) {

            List<String> col1 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> l2c1List1 = new ArrayList<>();
            List<Double> l2c2List1 = new ArrayList<>();
            List<Double> l2c3List1 = new ArrayList<>();
            List<Double> l3c1List1 = new ArrayList<>();
            List<Double> l4c1List1 = new ArrayList<>();
            List<Double> l4c2List1 = new ArrayList<>();
            List<Double> l5c1List1 = new ArrayList<>();
            List<Double> l5c2List1 = new ArrayList<>();
            List<Double> l5c3List1 = new ArrayList<>();

            for (int j = 0; j < multipoles.size(); ++j) {

                List<Double> col1List = new ArrayList<>();

                List<Double> l2c1List2 = new ArrayList<>();
                List<Double> l2c2List2 = new ArrayList<>();
                List<Double> l2c3List2 = new ArrayList<>();
                List<Double> l3c1List2 = new ArrayList<>();
                List<Double> l4c1List2 = new ArrayList<>();
                List<Double> l4c2List2 = new ArrayList<>();
                List<Double> l5c1List2 = new ArrayList<>();
                List<Double> l5c2List2 = new ArrayList<>();
                List<Double> l5c3List2 = new ArrayList<>();

                if ((mCol1[i].equals(mCol1[j]) && mCol2[i].equals(mCol2[j]) && mCol3[i].equals(mCol3[j]))
                        || (mCol1[i].equals(mCol1[j]) && mCol2[i].equals(mCol3[j]) && mCol3[i].equals(mCol2[j]))
                        || (mCol1[i].equals(mCol2[j]) && mCol2[i].equals(mCol1[j]) && mCol3[i].equals(mCol3[j]))
                        || (mCol1[i].equals(mCol2[j]) && mCol2[i].equals(mCol3[j]) && mCol3[i].equals(mCol1[j]))
                        || (mCol1[i].equals(mCol3[j]) && mCol2[i].equals(mCol1[j]) && mCol3[i].equals(mCol2[j]))
                        || (mCol1[i].equals(mCol3[j]) && mCol2[i].equals(mCol2[j]) && mCol3[i].equals(mCol1[j]))) {

                    //enter multipole values into lists for storage before averaging
                    col1.add(mValue1[j]);

                    //convert string list to double list for averaging
                    for (String multipoleN1 : col1) {

                        col1List.add(Double.parseDouble(multipoleN1));
                    }

                    l2c1List2.add(l2c1[j]);
                    l2c2List2.add(l2c2[j]);
                    l2c3List2.add(l2c3[j]);
                    l3c1List2.add(l3c1[j]);
                    l4c1List2.add(l4c1[j]);
                    l4c2List2.add(l4c2[j]);
                    l5c1List2.add(l5c1[j]);
                    l5c2List2.add(l5c2[j]);
                    l5c3List2.add(l5c3[j]);

                    col1List1 = col1List;

                    l2c1List1 = l2c1List2;
                    l2c2List1 = l2c2List2;
                    l2c3List1 = l2c3List2;
                    l3c1List1 = l3c1List2;
                    l4c1List1 = l4c1List2;
                    l4c2List1 = l4c2List2;
                    l5c1List1 = l5c1List2;
                    l5c2List1 = l5c2List2;
                    l5c3List1 = l5c3List2;

                }

            } //ends j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] l2c1Array = new Double[l2c1List1.size()];
            Double[] l2c2Array = new Double[l2c2List1.size()];
            Double[] l2c3Array = new Double[l2c3List1.size()];
            Double[] l3c1Array = new Double[l3c1List1.size()];
            Double[] l4c1Array = new Double[l4c1List1.size()];
            Double[] l4c2Array = new Double[l4c2List1.size()];
            Double[] l5c1Array = new Double[l5c1List1.size()];
            Double[] l5c2Array = new Double[l5c2List1.size()];
            Double[] l5c3Array = new Double[l5c3List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            l2c1List1.toArray(l2c1Array);
            l2c2List1.toArray(l2c2Array);
            l2c3List1.toArray(l2c3Array);
            l3c1List1.toArray(l3c1Array);
            l4c1List1.toArray(l4c1Array);
            l4c2List1.toArray(l4c2Array);
            l5c1List1.toArray(l5c1Array);
            l5c2List1.toArray(l5c2Array);
            l5c3List1.toArray(l5c3Array);

            int size1m = col1Array.length;
            int size2m = l2c1Array.length;
            int size3m = l2c2Array.length;
            int size4m = l2c3Array.length;
            int size5m = l3c1Array.length;
            int size6m = l4c1Array.length;
            int size7m = l4c2Array.length;
            int size8m = l5c1Array.length;
            int size9m = l5c2Array.length;
            int size10m = l5c3Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;
            double sum4 = 0;
            double sum5 = 0;
            double sum6 = 0;
            double sum7 = 0;
            double sum8 = 0;
            double sum9 = 0;
            double sum10 = 0;

            //average
            for (int k = 0; k < size1m; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1m;

            for (int k = 0; k < size2m; ++k) {
                sum2 = sum2 + l2c1Array[k];
            }
            double average2 = sum2 / size2m;

            for (int k = 0; k < size3m; ++k) {
                sum3 = sum3 + l2c2Array[k];
            }
            double average3 = sum3 / size3m;

            for (int k = 0; k < size4m; ++k) {
                sum4 = sum4 + l2c3Array[k];
            }
            double average4 = sum4 / size4m;

            for (int k = 0; k < size5m; ++k) {
                sum5 = sum5 + l3c1Array[k];
            }
            double average5 = sum5 / size5m;

            for (int k = 0; k < size6m; ++k) {
                sum6 = sum6 + l4c1Array[k];
            }
            double average6 = sum6 / size6m;

            for (int k = 0; k < size7m; ++k) {
                sum7 = sum7 + l4c2Array[k];
            }
            double average7 = sum7 / size7m;

            for (int k = 0; k < size8m; ++k) {
                sum8 = sum8 + l5c1Array[k];
            }
            double average8 = sum8 / size8m;

            for (int k = 0; k < size9m; ++k) {
                sum9 = sum9 + l5c2Array[k];
            }
            double average9 = sum9 / size9m;

            for (int k = 0; k < size10m; ++k) {
                sum10 = sum10 + l5c3Array[k];
            }
            double average10 = sum10 / size10m;

            DecimalFormat avgFormat = new DecimalFormat("0.00000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);
            String avg3 = avgFormat.format(average3);
            String avg4 = avgFormat.format(average4);
            String avg5 = avgFormat.format(average5);
            String avg6 = avgFormat.format(average6);
            String avg7 = avgFormat.format(average7);
            String avg8 = avgFormat.format(average8);
            String avg9 = avgFormat.format(average9);
            String avg10 = avgFormat.format(average10);

            String leadSpace = "                                        ";

            //ensure trace = 0
            Double trace;
            Double finall3c1D;
            Double finall4c2D;
            Double finall5c3D;

            String finall3c1 = null;
            String finall4c2 = null;
            String finall5c3 = null;

            trace = Double.parseDouble(avg5) + Double.parseDouble(avg7) + Double.parseDouble(avg10);

            if (trace == 0) {
                finall3c1D = Double.parseDouble(avg5);
                finall4c2D = Double.parseDouble(avg7);
                finall5c3D = Double.parseDouble(avg10);

                finall3c1 = avgFormat.format(finall3c1D);
                finall4c2 = avgFormat.format(finall4c2D);
                finall5c3 = avgFormat.format(finall5c3D);

            } else {
                Double third = trace / 3;
                finall3c1D = Double.parseDouble(avg5) - third;
                finall4c2D = Double.parseDouble(avg7) - third;
                finall5c3D = Double.parseDouble(avg10) - third;

                finall3c1 = avgFormat.format(finall3c1D);
                finall4c2 = avgFormat.format(finall4c2D);
                finall5c3 = avgFormat.format(finall5c3D);
            }

            //create line 2-5 values lists for final printing
            finall2c1L.add(avg2);
            finall2c2L.add(avg3);
            finall2c3L.add(avg4);
            finall3c1L.add(finall3c1);
            finall4c1L.add(avg6);
            finall4c2L.add(finall4c2);
            finall5c1L.add(avg8);
            finall5c2L.add(avg9);
            finall5c3L.add(finall5c3);

            //ensures proper spacing for lines with and without negative type numbers
            if (!" ".equals(mCol1[i]) && !" ".equals(mCol2[i]) && !" ".equals(mCol3[i])) {
                if (mCol1[i].contains("-") && mCol2[i].contains("-") && mCol3[i].contains("-")) {
                    multipoleWriter1.write("multipole  " + mCol1[i] + " " + mCol2[i] + " " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                } else if (mCol1[i].contains("-") && mCol2[i].contains("-")) {
                    multipoleWriter1.write("multipole  " + mCol1[i] + " " + mCol2[i] + "  " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                } else if (mCol1[i].contains("-") && mCol3[i].contains("-")) {
                    multipoleWriter1.write("multipole  " + mCol1[i] + "  " + mCol2[i] + " " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                } else if (mCol2[i].contains("-") && mCol3[i].contains("-")) {
                    multipoleWriter1.write("multipole   " + mCol1[i] + " " + mCol2[i] + " " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                } else if (mCol1[i].contains("-")) {
                    multipoleWriter1.write("multipole  " + mCol1[i] + "  " + mCol2[i] + "  " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                } else if (mCol2[i].contains("-")) {
                    multipoleWriter1.write("multipole   " + mCol1[i] + " " + mCol2[i] + "  " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                } else if (mCol3[i].contains("-")) {
                    multipoleWriter1.write("multipole   " + mCol1[i] + "  " + mCol2[i] + " " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                } else {
                    multipoleWriter1.write("multipole   " + mCol1[i] + "  " + mCol2[i] + "  " + mCol3[i] + "               " + avg1 + "\n"
                            + leadSpace + avg2 + "    " + avg3 + "    " + avg4 + "\n"
                            + leadSpace + finall3c1 + "\n"
                            + leadSpace + avg6 + "    " + finall4c2 + "\n"
                            + leadSpace + avg8 + "    " + avg9 + "    " + finall5c3 + "\n");
                }
            }

        }

        multipoleWriter1.close();

        //String result = null;
        //StringBuilder resultBuilder = new StringBuilder();
        Set<String> alreadyPresentm = new HashSet<>();
        int lineCounter = 0;

        //convert line 2-5 lists to arrays for final printing
        String[] l2c1A = new String[finall2c1L.size()];
        finall2c1L.toArray(l2c1A);
        String[] l2c2A = new String[finall2c2L.size()];
        finall2c2L.toArray(l2c2A);
        String[] l2c3A = new String[finall2c3L.size()];
        finall2c3L.toArray(l2c3A);
        String[] l3c1A = new String[finall3c1L.size()];
        finall3c1L.toArray(l3c1A);
        String[] l4c1A = new String[finall4c1L.size()];
        finall4c1L.toArray(l4c1A);
        String[] l4c2A = new String[finall4c2L.size()];
        finall4c2L.toArray(l4c2A);
        String[] l5c1A = new String[finall5c1L.size()];
        finall5c1L.toArray(l5c1A);
        String[] l5c2A = new String[finall5c2L.size()];
        finall5c2L.toArray(l5c2A);
        String[] l5c3A = new String[finall5c3L.size()];
        finall5c3L.toArray(l5c3A);

        //make sure there is only one multipole parameter (all five lines) for each atom combination
        try {
            FileReader fr = new FileReader("Omultipole.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder resultBuilder = new StringBuilder();
            String line;
            String leadSpace = "                                        ";
            while ((line = bufferedReader.readLine()) != null) {

                if (line.contains("multipole")) {

                    //set variables for comparison/individuality testing testing
                    String[] chars = line.split("(?!^)");
                    String[] test1 = new String[25];
                    String[] test2 = new String[25];
                    String[] test3 = new String[25];
                    String[] test4 = new String[25];
                    String[] test5 = new String[25];
                    String[] test6 = new String[25];

                    for (int k = 0; k < 25; ++k) {
                        test1[k] = chars[k];
                    }
                    for (int k = 0; k < 25; ++k) {
                        test2[k] = chars[k];
                    }
                    for (int k = 0; k < 25; ++k) {
                        test3[k] = chars[k];
                    }
                    for (int k = 0; k < 25; ++k) {
                        test4[k] = chars[k];
                    }
                    for (int k = 0; k < 25; ++k) {
                        test5[k] = chars[k];
                    }
                    for (int k = 0; k < 25; ++k) {
                        test6[k] = chars[k];
                    }

                    //switch the multipole atom terms for checking purposes
                    //if the original line (test1/testS) is: "multipole     ABC   DEF   GHI"
                    //test2/testS2: "multipole     ABC   GHI   DEF"
                    test2[16] = test1[21];
                    test2[17] = test1[22];
                    test2[18] = test1[23];
                    test2[19] = test1[24];

                    test2[21] = test1[16];
                    test2[22] = test1[17];
                    test2[23] = test1[18];
                    test2[24] = test1[19];

                    //test3/testS3: "multipole     DEF   ABC   GHI"
                    test3[11] = test1[16];
                    test3[12] = test1[17];
                    test3[13] = test1[18];
                    test3[14] = test1[19];

                    test3[16] = test1[11];
                    test3[17] = test1[12];
                    test3[18] = test1[13];
                    test3[19] = test1[14];

                    //test4/testS4: "multipole     DEF   GHI   ABC"
                    test4[11] = test1[16];
                    test4[12] = test1[17];
                    test4[13] = test1[18];
                    test4[14] = test1[19];

                    test4[16] = test1[21];
                    test4[17] = test1[22];
                    test4[18] = test1[23];
                    test4[19] = test1[24];

                    test4[21] = test1[11];
                    test4[22] = test1[12];
                    test4[23] = test1[13];
                    test4[24] = test1[14];

                    //test5/testS5: "multipole     GHI   ABC   DEF"
                    test5[11] = test1[21];
                    test5[12] = test1[22];
                    test5[13] = test1[23];
                    test5[14] = test1[24];

                    test5[16] = test1[11];
                    test5[17] = test1[12];
                    test5[18] = test1[13];
                    test5[19] = test1[14];

                    test5[21] = test1[16];
                    test5[22] = test1[17];
                    test5[23] = test1[18];
                    test5[24] = test1[19];

                    //test6/test6S: "multipole     GHI   DEF   ABC"
                    test6[11] = test1[21];
                    test6[12] = test1[22];
                    test6[13] = test1[23];
                    test6[14] = test1[24];

                    test6[21] = test1[11];
                    test6[22] = test1[12];
                    test6[23] = test1[13];
                    test6[24] = test1[14];

                    //converts test1-6 arrays to strings for comparison testing
                    String testS = Arrays.toString(test1);
                    String testS2 = Arrays.toString(test2);
                    String testS3 = Arrays.toString(test3);
                    String testS4 = Arrays.toString(test4);
                    String testS5 = Arrays.toString(test5);
                    String testS6 = Arrays.toString(test6);

                    boolean first = true;
                    if (!alreadyPresentm.contains(testS) || !alreadyPresentm.contains(testS2) || !alreadyPresentm.contains(testS3)
                            || !alreadyPresentm.contains(testS4) || !alreadyPresentm.contains(testS5) || !alreadyPresentm.contains(testS6)) {
                        if (first) {
                            first = false;
                        } else {
                            resultBuilder.append("\n");
                        }

                        if (!alreadyPresentm.contains(testS) || !alreadyPresentm.contains(testS2) || !alreadyPresentm.contains(testS3)
                                || !alreadyPresentm.contains(testS4) || !alreadyPresentm.contains(testS5) || !alreadyPresentm.contains(testS6)) {
                            resultBuilder.append(line).append("\n"); //l1
                            resultBuilder.append(leadSpace).append(l2c1A[lineCounter]).append("    ").append(l2c2[lineCounter]).append("    ").append(l2c3A[lineCounter]).append("\n"); //l2
                            resultBuilder.append(leadSpace).append(l3c1A[lineCounter]).append("\n");//l3
                            resultBuilder.append(leadSpace).append(l4c1A[lineCounter]).append("    ").append(l4c2A[lineCounter]).append("\n");//l4
                            resultBuilder.append(leadSpace).append(l5c1A[lineCounter]).append("    ").append(l5c2A[lineCounter]).append("    ").append(l5c3A[lineCounter]);//l5
                            resultBuilder.append("\n");
                        }

                        alreadyPresentm.add(testS);
                        alreadyPresentm.add(testS2);
                        alreadyPresentm.add(testS3);
                        alreadyPresentm.add(testS4);
                        alreadyPresentm.add(testS5);
                        alreadyPresentm.add(testS6);

                    }
                    String result = resultBuilder.toString();

                    //write final multipole values to patch file
                    try {

                        try (PrintWriter multipoleWriter2 = new PrintWriter("multipole.patch")) {
                            multipoleWriter2.write(result);
                            multipoleWriter2.write("\n");
                            multipoleWriter2.close();
                        }

                    } catch (IOException e) {
                    }

                    ++lineCounter;
                } //ends if(line.contains("multipole"))
            } //ends reader while

        } catch (IOException e) {
        }

        System.out.println("Multipole Patch: COMPLETE");
    }

    /**
     * <p>writePolarize.</p>
     *
     * @throws java.io.FileNotFoundException if any.
     */
    public void writePolarize() throws FileNotFoundException {
        String polarAtom = null;

        PrintWriter polarWriter1 = new PrintWriter("Opolarize.patch");

        String type = null;

        List<String> polar4final = new ArrayList<>();
        List<String> polar5final = new ArrayList<>();
        List<String> polar6final = new ArrayList<>();

        List<String> pPresent = new ArrayList<>();

        for (int i = 0; i < max; ++i) {

            List<String> polarTypes4 = new ArrayList<>();
            List<String> polarTypes5 = new ArrayList<>();
            List<String> polarTypes6 = new ArrayList<>();

            List<String> col1 = new ArrayList<>();
            List<String> col2 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> col2List1 = new ArrayList<>();

            int ptcount = 0;
            String check4 = "500";
            String check5 = "500";
            String check6 = "500";

            for (int j = i; j < max; ++j) {

                List<Double> col1List = new ArrayList<>();
                List<Double> col2List = new ArrayList<>();

                if (atomName[i].equals(atomName[j])) {

                    //find corresponding fragCode value for index j
                    //from this fragCode, using the fragCodeMap, find matching patch file
                    String code = fragCode[j];
                    String num = typeNum[j];
                    type = vTypeNum[j];

                    String patchname = patchMap.get(Integer.parseInt(code));
                    String polar;

                    //read in patch file according to mapped fragCode
                    try {
                        File file = new File(patchname);
                        try (FileReader fileReader = new FileReader(file)) {
                            BufferedReader bufferedReader = new BufferedReader(fileReader);
                            StringBuilder stringBuffer = new StringBuilder();
                            String line;
                            while ((line = bufferedReader.readLine()) != null) {
                                stringBuffer.append(line);
                                stringBuffer.append("\n");

                                String patch = line;

                                //find polarize values
                                //find proper typeNum in polarize values
                                if (patch.contains("polarize")) {
                                    polar = line;
                                    String[] polarNums2 = polar.split(" ");

                                    List<String> list1 = new ArrayList<>(Arrays.asList(polarNums2));
                                    list1.removeAll(Arrays.asList(""));
                                    String[] polarNums3 = new String[list1.size()];
                                    list1.toArray(polarNums3);

                                    if (num.equals(polarNums3[1])) {
                                        polarAtom = line;

                                        String[] polarNums1 = polarAtom.split(" ");

                                        //remove white space
                                        List<String> list = new ArrayList<>(Arrays.asList(polarNums1));
                                        list.removeAll(Arrays.asList(""));
                                        String[] polarNums = new String[list.size()];
                                        list.toArray(polarNums);

                                        //write polarize values for proper typeNum into a file to work with them
                                        col1.add(polarNums[2]);
                                        col2.add(polarNums[3]);

                                        //write the trailing values for each line to respective arrays for further processing
                                        if (polarNums.length < 5) {
                                            polarTypes4.add(" ");
                                            polarTypes5.add(" ");
                                            polarTypes6.add(" ");
                                        }
                                        if (polarNums.length == 5) {
                                            polarTypes4.add(polarNums[4]);
                                            polarTypes5.add(" ");
                                            polarTypes6.add(" ");
                                        }
                                        if (polarNums.length == 6) {
                                            if (polarNums[5] != null) {
                                                polarTypes4.add(polarNums[4]);
                                                polarTypes5.add(polarNums[5]);
                                                polarTypes6.add(" ");
                                            }
                                        }
                                        if (polarNums.length == 7) {
                                            if (polarNums[6] != null) {
                                                polarTypes4.add(polarNums[4]);
                                                polarTypes5.add(polarNums[5]);
                                                polarTypes6.add(polarNums[6]);
                                            }
                                        }

                                        //convert string list to double list for averaging
                                        for (String polarN1 : col1) {
                                            col1List.add(Double.parseDouble(polarN1));
                                        }
                                        for (String polarN2 : col2) {
                                            col2List.add(Double.parseDouble(polarN2));
                                        }

                                        col1List1 = col1List;
                                        col2List1 = col2List;

                                        String polarTypes4S2 = null;
                                        String polarTypes5S2 = null;
                                        String polarTypes6S2 = null;
                                        String polarTypes4S3 = null;
                                        String polarTypes5S3 = null;
                                        String polarTypes6S3 = null;

                                        //polarTypes 4/5/6S are strings that contain 4** type numbers
                                        String[] polarTypes4A = new String[polarTypes4.size()];
                                        polarTypes4.toArray(polarTypes4A);
                                        String polarTypes4S = polarTypes4A[0];
                                        String[] polarTypes5A = new String[polarTypes5.size()];
                                        polarTypes5.toArray(polarTypes5A);
                                        String polarTypes5S = polarTypes5A[0];
                                        String[] polarTypes6A = new String[polarTypes6.size()];
                                        polarTypes6.toArray(polarTypes6A);
                                        String polarTypes6S = polarTypes6A[0];

                                        //for subsequent additions to polarTypes4/5/6 lists and polarTypes4/5/6A arrays
                                        if (polarTypes4A.length == 2) {
                                            polarTypes4S2 = polarTypes4A[1];
                                            polarTypes5S2 = polarTypes5A[1];
                                            polarTypes6S2 = polarTypes6A[1];
                                        }
                                        if (polarTypes4A.length == 3) {
                                            polarTypes4S2 = polarTypes4A[1];
                                            polarTypes5S2 = polarTypes5A[1];
                                            polarTypes6S2 = polarTypes6A[1];

                                            polarTypes4S3 = polarTypes4A[2];
                                            polarTypes5S3 = polarTypes5A[2];
                                            polarTypes6S3 = polarTypes6A[2];
                                        }

                                        //polarTypes4/5/6A2 will contain 5** type numbers
                                        String polarTypes4A2;
                                        String polarTypes5A2;
                                        String polarTypes6A2;

                                        //converts read-in type number (4**) to final type numbers (5**)
                                        //sets first column of 4** numbers to their corresponding 5** numbers
                                        if (polarTypes4S != null && !" ".equals(polarTypes4S)) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes4A2 = Map4to51.get(polarTypes4S);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes4A2 = Map4to52.get(polarTypes4S);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes4A2 = Map4to53.get(polarTypes4S);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes4A2 = Map4to54.get(polarTypes4S);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes4A2 = Map4to55.get(polarTypes4S);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes4A2 = Map4to56.get(polarTypes4S);
                                            } else {
                                                polarTypes4A2 = " ";
                                            }
                                        } else {
                                            polarTypes4A2 = " ";
                                        }

                                        //if the polarize line has 2 ending typeNum values
                                        //sets second column of 4** numbers to their corresponding 5** numbers
                                        if (!" ".equals(polarTypes5S) && polarTypes5S != null) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes5A2 = Map4to51.get(polarTypes5S);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes5A2 = Map4to52.get(polarTypes5S);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes5A2 = Map4to53.get(polarTypes5S);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes5A2 = Map4to54.get(polarTypes5S);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes5A2 = Map4to55.get(polarTypes5S);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes5A2 = Map4to56.get(polarTypes5S);
                                            } else {
                                                polarTypes5A2 = " ";
                                            }
                                        } else {
                                            polarTypes5A2 = " ";
                                        }

                                        //if the polarize line has 3 ending typeNum values
                                        //sets third column of 4** numbers to their corresponding 5** numbers
                                        if (!" ".equals(polarTypes6S) && polarTypes6S != null) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes6A2 = Map4to51.get(polarTypes6S);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes6A2 = Map4to52.get(polarTypes6S);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes6A2 = Map4to53.get(polarTypes6S);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes6A2 = Map4to54.get(polarTypes6S);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes6A2 = Map4to55.get(polarTypes6S);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes6A2 = Map4to56.get(polarTypes6S);
                                            } else {
                                                polarTypes6A2 = " ";
                                            }
                                        } else {
                                            polarTypes6A2 = " ";
                                        }

                                        String polarTypes4A3 = null;
                                        String polarTypes5A3 = null;
                                        String polarTypes6A3 = null;

                                        //for 4/5/6S2 values: second read-in; only happens if their are 2+ similar atoms
                                        //convert ending typeNum values to vTypeNum values
                                        //sets first column of 4** numbers to their corresponding 5** numbers
                                        if (polarTypes4S2 != null && !" ".equals(polarTypes4S2)) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes4A3 = Map4to51.get(polarTypes4S2);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes4A3 = Map4to52.get(polarTypes4S2);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes4A3 = Map4to53.get(polarTypes4S2);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes4A3 = Map4to54.get(polarTypes4S2);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes4A3 = Map4to55.get(polarTypes4S2);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes4A3 = Map4to56.get(polarTypes4S2);
                                            } else {
                                                polarTypes4A3 = " ";
                                            }
                                        } else {
                                            polarTypes4A3 = " ";
                                        }

                                        //if the polarize line has 2 ending typeNum values
                                        //sets second column of 4** numbers to their corresponding 5** numbers
                                        if (!" ".equals(polarTypes5S2) && polarTypes5S2 != null) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes5A3 = Map4to51.get(polarTypes5S2);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes5A3 = Map4to52.get(polarTypes5S2);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes5A3 = Map4to53.get(polarTypes5S2);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes5A3 = Map4to54.get(polarTypes5S2);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes5A3 = Map4to55.get(polarTypes5S2);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes5A3 = Map4to56.get(polarTypes5S2);
                                            } else {
                                                polarTypes5A3 = " ";
                                            }
                                        } else {
                                            polarTypes5A3 = " ";
                                        }

                                        //if the polarize line has 3 ending typeNum values
                                        //sets third column of 4** numbers to their corresponding 5** numbers
                                        if (!" ".equals(polarTypes6S2) && polarTypes6S2 != null) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes6A3 = Map4to51.get(polarTypes6S2);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes6A3 = Map4to52.get(polarTypes6S2);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes6A3 = Map4to53.get(polarTypes6S2);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes6A3 = Map4to54.get(polarTypes6S2);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes6A3 = Map4to55.get(polarTypes6S2);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes6A3 = Map4to56.get(polarTypes6S2);
                                            } else {
                                                polarTypes6A3 = " ";
                                            }
                                        } else {
                                            polarTypes6A3 = " ";
                                        }

                                        String polarTypes4A4 = null;
                                        String polarTypes5A4 = null;
                                        String polarTypes6A4 = null;

                                        //for 4/5/6S3 values: third read-in; only happens if there are 3+ similar atoms
                                        //convert ending typeNum values to vTypeNum values
                                        //sets first column of 4** numbers to their corresponding 5** numbers
                                        if (polarTypes4S3 != null && !" ".equals(polarTypes4S3)) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes4A4 = Map4to51.get(polarTypes4S3);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes4A4 = Map4to52.get(polarTypes4S3);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes4A4 = Map4to53.get(polarTypes4S3);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes4A4 = Map4to54.get(polarTypes4S3);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes4A4 = Map4to55.get(polarTypes4S3);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes4A4 = Map4to56.get(polarTypes4S3);
                                            } else {
                                                polarTypes4A4 = " ";
                                            }
                                        } else {
                                            polarTypes4A4 = " ";
                                        }

                                        //if the polarize line has 2 ending typeNum values
                                        //sets second column of 4** numbers to their corresponding 5** numbers
                                        if (!" ".equals(polarTypes5S3) && polarTypes5S3 != null) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes5A4 = Map4to51.get(polarTypes5S3);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes5A4 = Map4to52.get(polarTypes5S3);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes5A4 = Map4to53.get(polarTypes5S3);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes5A4 = Map4to54.get(polarTypes5S3);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes5A4 = Map4to55.get(polarTypes5S3);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes5A4 = Map4to56.get(polarTypes5S3);
                                            } else {
                                                polarTypes5A4 = " ";
                                            }
                                        } else {
                                            polarTypes5A4 = " ";
                                        }

                                        //if the polarize line has 3 ending typeNum values
                                        //sets third column of 4** numbers to their corresponding 5** numbers
                                        if (!" ".equals(polarTypes6S3) && polarTypes6S3 != null) {
                                            if (Integer.parseInt(code) == 1) {
                                                polarTypes6A4 = Map4to51.get(polarTypes6S3);
                                            } else if (Integer.parseInt(code) == 2) {
                                                polarTypes6A4 = Map4to52.get(polarTypes6S3);
                                            } else if (Integer.parseInt(code) == 3) {
                                                polarTypes6A4 = Map4to53.get(polarTypes6S3);
                                            } else if (Integer.parseInt(code) == 4) {
                                                polarTypes6A4 = Map4to54.get(polarTypes6S3);
                                            } else if (Integer.parseInt(code) == 5) {
                                                polarTypes6A4 = Map4to55.get(polarTypes6S3);
                                            } else if (Integer.parseInt(code) == 6) {
                                                polarTypes6A4 = Map4to56.get(polarTypes6S3);
                                            } else {
                                                polarTypes6A4 = " ";
                                            }
                                        } else {
                                            polarTypes6A4 = " ";
                                        }

                                        //similar atom incrimentor
                                        ++ptcount;

                                        //for similar atoms in many patches, certain polarizing atoms may not be listed
                                        //in every patch.  These "extra" values ensure that all atoms that contribute
                                        //to the polarization of the atom in question are listed in the "polarize" section
                                        //of the final, vem patch file
                                        String extra1 = null;
                                        String extra2 = null;

                                        //check and see if all recorded values are already listed in the final "polarize"
                                        //line for the atom in question.
                                        //if not, they are entered into "extra1" or "extra2" variables for adding to
                                        //final "polarize" types array for printing in final patch file
                                        if (ptcount == 2) {
                                            if (!" ".equals(polarTypes4A3) && polarTypes4A3 != null && !check4.equals(polarTypes4A3) && !check5.equals(polarTypes4A3)
                                                    && !check6.equals(polarTypes4A3) && extra1 == null) {
                                                extra1 = polarTypes4A3;
                                            }
                                            if (!" ".equals(polarTypes5A3) && polarTypes5A3 != null && !check4.equals(polarTypes5A3) && !check5.equals(polarTypes5A3)
                                                    && !check6.equals(polarTypes5A3) && extra1 == null) {
                                                extra1 = polarTypes5A3;
                                            } else if (!" ".equals(polarTypes5A3) && polarTypes5A3 != null && !check4.equals(polarTypes5A3) && !check5.equals(polarTypes5A3)
                                                    && !check6.equals(polarTypes5A3) && extra2 == null) {
                                                extra2 = polarTypes5A3;
                                            }
                                            if (!" ".equals(polarTypes6A3) && polarTypes6A3 != null && !check4.equals(polarTypes6A3) && !check5.equals(polarTypes6A3)
                                                    && !check6.equals(polarTypes6A3) && extra1 == null) {
                                                extra1 = polarTypes6A3;
                                            } else if (!" ".equals(polarTypes6A3) && polarTypes6A3 != null && !check4.equals(polarTypes6A3) && !check5.equals(polarTypes6A3)
                                                    && !check6.equals(polarTypes6A3) && extra2 == null) {
                                                extra2 = polarTypes6A3;
                                            }

                                        }

                                        //identifies extra atoms, not listed in the first atom's "polarize" line, and stores
                                        //them for addition to the final "polarize" line for the atom in question
                                        if (ptcount == 3) {
                                            if (!" ".equals(polarTypes4A4) && polarTypes4A4 != null && !check4.equals(polarTypes4A4) && !check5.equals(polarTypes4A4)
                                                    && !check6.equals(polarTypes4A4)) {
                                                if (extra1 == null) {
                                                    extra1 = polarTypes4A4;
                                                } else if (extra2 == null) {
                                                    extra2 = polarTypes4A4;
                                                } else {
                                                    System.out.println("Need another extra variable");
                                                }
                                            }
                                            if (!" ".equals(polarTypes5A4) && polarTypes5A4 != null && !check4.equals(polarTypes5A4) && !check5.equals(polarTypes5A4)
                                                    && !check6.equals(polarTypes5A4)) {
                                                if (extra1 == null) {
                                                    extra1 = polarTypes5A4;
                                                } else if (extra2 == null) {
                                                    extra2 = polarTypes5A4;
                                                } else {
                                                    System.out.println("Need another extra variable");
                                                }
                                            }
                                            if (!" ".equals(polarTypes6A4) && polarTypes6A4 != null && !check4.equals(polarTypes6A4) && !check5.equals(polarTypes6A4)
                                                    && !check6.equals(polarTypes6A4)) {
                                                if (extra1 == null) {
                                                    extra1 = polarTypes6A4;
                                                } else if (extra2 == null) {
                                                    extra2 = polarTypes6A4;
                                                } else {
                                                    System.out.println("Need another extra value");
                                                }
                                            }

                                        }
                                        //if the extra values are the same, then only one is needed
                                        //so extra2 is set back to null
                                        if (extra1 == null ? extra2 == null : extra1.equals(extra2)) {
                                            extra2 = null;
                                        }

                                        //add updated polarize types to a final list for printing
                                        List<String> pt4Aa4 = new ArrayList<>();
                                        List<String> pt5Aa4 = new ArrayList<>();
                                        List<String> pt6Aa4 = new ArrayList<>();

                                        for (int n = 0; n < 1; ++n) {

                                            boolean first = true;
                                            if (!pPresent.contains(atomName[i]) || atomName[i].contains("X")) {
                                                if (first) {
                                                    first = false;
                                                } else {
                                                    pt4Aa4.add("\n");
                                                    pt5Aa4.add("\n");
                                                    pt6Aa4.add("\n");
                                                }

                                                if (!pPresent.contains(atomName[i]) || atomName[i].contains("X")) {
                                                    pt4Aa4.add(polarTypes4A2);
                                                    pt5Aa4.add(polarTypes5A2);
                                                    pt6Aa4.add(polarTypes6A2);
                                                }

                                                pPresent.add(atomName[i]);

                                            }
                                        }

                                        String[] polarTypes4Aa4 = new String[1];
                                        pt4Aa4.toArray(polarTypes4Aa4);
                                        String[] polarTypes5Aa4 = new String[1];
                                        pt5Aa4.toArray(polarTypes5Aa4);
                                        String[] polarTypes6Aa4 = new String[1];
                                        pt6Aa4.toArray(polarTypes6Aa4);

                                        //ensures proper addition of the identified "extra" atoms to the final "polarize" line
                                        //for the atom in question
                                        if (!check4.equals(polarTypes6Aa4[0]) && !check5.equals(polarTypes6Aa4[0])
                                                && !check6.equals(polarTypes6Aa4[0]) && extra1 != null) {
                                            polarTypes6Aa4[0] = extra1;
                                        } else if (!check4.equals(polarTypes5Aa4[0]) && !check5.equals(polarTypes5Aa4[0])
                                                && !check6.equals(polarTypes5Aa4[0]) && extra1 != null) {
                                            polarTypes5Aa4[0] = extra1;
                                        }

                                        if (!check4.equals(polarTypes6Aa4[0]) && !check5.equals(polarTypes6Aa4[0])
                                                && !check6.equals(polarTypes6Aa4[0]) && extra2 != null) {
                                            polarTypes6Aa4[0] = extra2;
                                        }

                                        //updates the values used for checking for extras
                                        //sets check4/5/6 to the polarize 5** values for the first occurance
                                        //of the atom in question.
                                        if (ptcount == 1 && polarTypes4Aa4[0] != null && polarTypes5Aa4[0] != null && polarTypes6Aa4 != null) {
                                            check4 = polarTypes4Aa4[0];
                                            check5 = polarTypes5Aa4[0];
                                            check6 = polarTypes6Aa4[0];
                                        }

                                        //adds "polarize" line trailing values to final arrays for
                                        //printing as part of the final "polarize" line for the atom in question
                                        //in the final patch file.
                                        for (int count = 0; count < polarTypes4Aa4.length; ++count) {
                                            polar4final.add(polarTypes4Aa4[count]);
                                        }
                                        for (int count = 0; count < polarTypes5Aa4.length; ++count) {
                                            polar5final.add(polarTypes5Aa4[count]);
                                        }
                                        for (int count = 0; count < polarTypes6Aa4.length; ++count) {
                                            polar6final.add(polarTypes6Aa4[count]);
                                        }
                                    }
                                } //end if(patch.contains("polarize");
                            } //end reader while
                        } //end try

                    } catch (IOException e) {
                    }

                } //end if(atomName[i].equals(atomName[j]))

            } //end j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] col2Array = new Double[col2List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            col2List1.toArray(col2Array);

            int size1p = col1Array.length;
            int size2p = col2Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;

            for (int k = 0; k < size1p; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1p;

            for (int k = 0; k < size2p; ++k) {
                sum2 = sum2 + col2Array[k];
            }
            double average2 = sum2 / size2p;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);

            //convert lists to arrays for printing
            String[] polarize4null = new String[polar4final.size()];
            polar4final.toArray(polarize4null);
            String[] polarize5null = new String[polar5final.size()];
            polar5final.toArray(polarize5null);
            String[] polarize6null = new String[polar6final.size()];
            polar6final.toArray(polarize6null);

            //remove nulls
            List<String> list4 = new ArrayList<>(Arrays.asList(polarize4null));
            list4.removeAll(Collections.singleton(null));
            String[] polarize4final = new String[list4.size()];
            list4.toArray(polarize4final);
            List<String> list5 = new ArrayList<>(Arrays.asList(polarize5null));
            list5.removeAll(Collections.singleton(null));
            String[] polarize5final = new String[list5.size()];
            list5.toArray(polarize5final);
            List<String> list6 = new ArrayList<>(Arrays.asList(polarize6null));
            list6.removeAll(Collections.singleton(null));
            String[] polarize6final = new String[list6.size()];
            list6.toArray(polarize6final);

            //set incrimentor variables for writer1
            int pf4 = polarize4final.length;
            int pf5 = polarize5final.length;
            int pf6 = polarize6final.length;

            //write values to Opolarize.patch (original polarize values, hence "Opolarize")
            polarWriter1.write("polarize     " + type + "          " + avg1 + "     " + avg2 + " " + polarize4final[pf4 - 1] + " " + polarize5final[pf5 - 1] + " " + polarize6final[pf6 - 1] + "\n");

        }

        polarWriter1.close();

        //String result;
        //StringBuilder resultBuilder = new StringBuilder();
        Set<String> alreadyPresent = new HashSet<>();

        //make sure there is only one polarize parameter line for each atom type
        try {
            FileReader fr = new FileReader("Opolarize.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder resultBuilder = new StringBuilder();
            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] parts = line.split(" ");

                if (parts[5].contains("5")) {
                    boolean first = true;
                    if (!alreadyPresent.contains(parts[5])) {
                        if (first) {
                            first = false;
                        } else {
                            resultBuilder.append("\n");
                        }

                        if (!alreadyPresent.contains(parts[5])) {
                            resultBuilder.append(line);
                            resultBuilder.append("\n");
                        }

                        alreadyPresent.add(parts[5]);

                    }
                    String result = resultBuilder.toString();

                    //write final polarize values to patch file
                    try {

                        try (PrintWriter polarWriter2 = new PrintWriter("polarize.patch")) {
                            polarWriter2.write(result);
                            polarWriter2.write("\n");
                            polarWriter2.close();
                        }

                    } catch (IOException e) {
                    }
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Polarize Patch: COMPLETE");
    }

    /**
     * <p>writeVDW.</p>
     *
     * @throws java.io.FileNotFoundException if any.
     */
    public void writeVDW() throws FileNotFoundException {
        String vdwAtom = null;

        PrintWriter vdwWriter1 = new PrintWriter("Ovdw.patch");

        String type2 = null;

        for (int i = 0; i < max; ++i) {

            List<String> col1 = new ArrayList<>();
            List<String> col2 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> col2List1 = new ArrayList<>();

            for (int j = i; j < max; ++j) {

                List<Double> col1List = new ArrayList<>();
                List<Double> col2List = new ArrayList<>();

                if (atomName[i].equals(atomName[j])) {
                    //find corresponding fragCode value for index i
                    //from this fragCode, using the fragCodeMap, find matching patch file
                    String code = fragCode[j];
                    String num = typeNum[j];
                    type2 = vTypeNum[j];

                    String patchname = patchMap.get(Integer.parseInt(code));
                    String vdw;

                    //read in patch file according to mapped fragCode
                    try {
                        File file = new File(patchname);
                        try (FileReader fileReader = new FileReader(file)) {
                            BufferedReader bufferedReader = new BufferedReader(fileReader);
                            StringBuilder stringBuffer = new StringBuilder();
                            String line;
                            while ((line = bufferedReader.readLine()) != null) {
                                stringBuffer.append(line);
                                stringBuffer.append("\n");

                                String patch = line;

                                //find vdw values
                                //find proper typeNum in vdw values
                                if (patch.contains("vdw")) {
                                    vdw = line;
                                    if (vdw.contains(num)) {
                                        vdwAtom = line;

                                        String[] vdwNums = line.split(" ");

                                        //write vdw values for proper typeNum into a file to work with them
                                        col1.add(vdwNums[9]);
                                        col2.add(vdwNums[12]);

                                        //convert string list to double list for averaging
                                        for (String vdwN1 : col1) {
                                            col1List.add(Double.parseDouble(vdwN1));
                                        }
                                        for (String vdwN2 : col2) {
                                            col2List.add(Double.parseDouble(vdwN2));
                                        }

                                        col1List1 = col1List;
                                        col2List1 = col2List;

                                    }
                                } //end if(patch.contains("vdw"))
                            } //end reader while
                        } //end try

                    } catch (IOException e) {
                    }
                } //end if(atomName[i].equals(atomName[j]))
            } //end j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] col2Array = new Double[col2List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            col2List1.toArray(col2Array);

            int size1v = col1Array.length;
            int size2v = col2Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;

            for (int k = 0; k < size1v; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1v;

            for (int k = 0; k < size2v; ++k) {
                sum2 = sum2 + col2Array[k];
            }
            double average2 = sum2 / size2v;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);

            vdwWriter1.write("vdw       " + type2 + "  " + avg1 + "   " + avg2 + "\n");

        }

        vdwWriter1.close();

        //String resultV = null;
        //StringBuilder resultBuilderV = new StringBuilder();
        Set<String> alreadyPresentV = new HashSet<>();

        //make sure there is only one vdw parameter line for each atom type
        try {
            FileReader fr = new FileReader("Ovdw.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder resultBuilder = new StringBuilder();
            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] parts = line.split(" ");

                boolean first = true;
                if (!alreadyPresentV.contains(parts[7])) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilder.append("\n");
                    }

                    if (!alreadyPresentV.contains(parts[7])) {
                        resultBuilder.append(line);
                        resultBuilder.append("\n");
                    }

                    alreadyPresentV.add(parts[7]);

                }
                String result = resultBuilder.toString();

                //write final vdw values to patch file
                try {

                    PrintWriter vdwWriter2 = new PrintWriter("vdw.patch");
                    vdwWriter2.write(result);
                    vdwWriter2.close();

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("van der Waals Patch: COMPLETE");
    }

    /**
     * <p>writeBonded.</p>
     *
     * @throws java.io.FileNotFoundException if any.
     */
    public void writeBonded() throws FileNotFoundException {
        List<String> bonds = new ArrayList<>();
        List<String> firstNList = new ArrayList<>();
        List<String> secondNList = new ArrayList<>();
        List<String> value1List = new ArrayList<>();
        List<String> value2List = new ArrayList<>();

        for (int j = 0; j < max; ++j) {

            String code = fragCode[j];
            String patchname = patchMap.get(Integer.parseInt(code));
            String bond = null;

            try {
                File file = new File(patchname);
                FileReader fileReader = new FileReader(file);
                BufferedReader bufferedReader = new BufferedReader(fileReader);
                StringBuilder stringBuffer = new StringBuilder();
                String line;

                while ((line = bufferedReader.readLine()) != null) {
                    stringBuffer.append(line);
                    stringBuffer.append("\n");

                    String patch = line;

                    if (patch.contains("bond")) {
                        bond = line;

                        String[] bond1 = (bond.split(" "));
                        String[] bond2 = new String[bond1.length];

                        for (int k = 0; k < bond1.length; ++k) {
                            bond2[k] = bond1[k];
                        }

                        //replaces fragment atom type numbers with vem. atom type numbers
                        if (Integer.parseInt(code) == 1) {

                            bond2[6] = Map4to51.get(bond1[6]);
                            bond2[9] = Map4to51.get(bond1[9]);

                        } else if (Integer.parseInt(code) == 2) {

                            bond2[6] = Map4to52.get(bond1[6]);
                            bond2[9] = Map4to52.get(bond1[9]);

                        } else if (Integer.parseInt(code) == 3) {

                            bond2[6] = Map4to53.get(bond1[6]);
                            bond2[9] = Map4to53.get(bond1[9]);

                        } else if (Integer.parseInt(code) == 4) {

                            bond2[6] = Map4to54.get(bond1[6]);
                            bond2[9] = Map4to54.get(bond1[9]);

                        } else if (Integer.parseInt(code) == 5) {

                            bond2[6] = Map4to55.get(bond1[6]);
                            bond2[9] = Map4to55.get(bond1[9]);

                        } else if (Integer.parseInt(code) == 6) {

                            bond2[6] = Map4to56.get(bond1[6]);
                            bond2[9] = Map4to56.get(bond1[9]);

                        } else {
                            System.out.println("Not in a patch");
                        }

                        firstNList.add(bond2[6]);
                        secondNList.add(bond2[9]);
                        value1List.add(bond2[12]);
                        value2List.add(bond2[16]);

                        //add to bonds list
                        bonds.add("bond      " + bond2[6] + "   " + bond2[9] + "   " + bond2[12] + "    " + bond2[16] + "\n");

                    } //ends if(patch.contains("bond"))
                } //ends reader while

            } catch (IOException e) {
            }
        } //ends j loop

        //define new arrays to enter list data into
        String[] firstN = new String[firstNList.size()];
        String[] secondN = new String[secondNList.size()];
        String[] value1 = new String[value1List.size()];
        String[] value2 = new String[value2List.size()];

        //convert string lists to string arrays for further use
        firstNList.toArray(firstN);
        secondNList.toArray(secondN);
        value1List.toArray(value1);
        value2List.toArray(value2);

        PrintWriter bondWriter1 = new PrintWriter("Obond.patch");

        for (int i = 0; i < bonds.size(); ++i) {

            List<String> col1 = new ArrayList<>();
            List<String> col2 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> col2List1 = new ArrayList<>();

            for (int j = i; j < bonds.size(); ++j) {

                List<Double> col1List = new ArrayList<>();
                List<Double> col2List = new ArrayList<>();

                if (firstN[i].equals(firstN[j]) && secondN[i].equals(secondN[j])) {

                    col1.add(value1[j]);
                    col2.add(value2[j]);

                    //convert string list to double list for averaging
                    for (String bondN1 : col1) {

                        col1List.add(Double.parseDouble(bondN1));
                    }
                    for (String bondN2 : col2) {
                        col2List.add(Double.parseDouble(bondN2));
                    }

                    col1List1 = col1List;
                    col2List1 = col2List;

                }
            } //ends j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] col2Array = new Double[col2List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            col2List1.toArray(col2Array);

            int size1b = col1Array.length;
            int size2b = col2Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;

            for (int k = 0; k < size1b; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1b;

            for (int k = 0; k < size2b; ++k) {
                sum2 = sum2 + col2Array[k];
            }
            double average2 = sum2 / size2b;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);

            if (secondN[i] != null) {
                bondWriter1.write("bond       " + firstN[i] + "   " + secondN[i] + "   " + avg1 + "   " + avg2 + "\n");
            }
        } //ends i "for" loop

        bondWriter1.close();

        String resultB;
        StringBuilder resultBuilderb = new StringBuilder();
        Set<String> alreadyPresentb = new HashSet<>();

        //make sure there is only one bond parameter line for each atom combination
        try {
            FileReader fr = new FileReader("Obond.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] chars = line.split("(?!^)");
                String[] test1 = new String[20];
                String[] test2 = new String[20];

                for (int k = 0; k < 20; ++k) {
                    test1[k] = chars[k];
                }

                for (int k = 0; k < 20; ++k) {
                    test2[k] = chars[k];
                }

                //switch the bonded atom terms for checking purposes
                test2[11] = test1[17];
                test2[12] = test1[18];
                test2[13] = test1[19];
                test2[17] = test1[11];
                test2[18] = test1[12];
                test2[19] = test1[13];

                String testS = Arrays.toString(test1);
                String testS2 = Arrays.toString(test2);

                boolean first = true;
                if (!alreadyPresentb.contains(testS) || !alreadyPresentb.contains(testS2)) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilderb.append("\n");
                    }

                    if (!alreadyPresentb.contains(testS) || !alreadyPresentb.contains(testS2)) {
                        resultBuilderb.append(line);
                        resultBuilderb.append("\n");
                    }

                    alreadyPresentb.add(testS);
                    alreadyPresentb.add(testS2);

                }
                resultB = resultBuilderb.toString();

                //write final bond values to patch file
                try {

                    try (PrintWriter bondWriter2 = new PrintWriter("bond.patch")) {
                        bondWriter2.write(resultB);
                        bondWriter2.close();
                    }

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Bond Patch: COMPLETE");

        /**
         * *
         * End bond averaging/writing code
         */
        /**
         * *
         * Angle averaging/writing code
         */
        List<String> angles = new ArrayList<>();
        List<String> aCol1List = new ArrayList<>();
        List<String> aCol2List = new ArrayList<>();
        List<String> aCol3List = new ArrayList<>();
        List<String> aValue1List = new ArrayList<>();
        List<String> aValue2List = new ArrayList<>();

        for (int j = 0; j < max; ++j) {

            String code = fragCode[j];
            String patchname = patchMap.get(Integer.parseInt(code));
            String angle = null;

            try {
                File file = new File(patchname);
                FileReader fileReader = new FileReader(file);
                BufferedReader bufferedReader = new BufferedReader(fileReader);
                StringBuilder stringBuffer = new StringBuilder();
                String line;

                while ((line = bufferedReader.readLine()) != null) {
                    stringBuffer.append(line);
                    stringBuffer.append("\n");

                    String patch = line;

                    if (patch.contains("angle")) {
                        angle = line;

                        String[] angle1 = angle.split(" ");
                        String[] angle2 = new String[angle1.length];

                        for (int k = 0; k < angle1.length; ++k) {
                            angle2[k] = angle1[k];
                        }

                        //replaces fragment atom type numbers with vem. atom type numbers
                        if (Integer.parseInt(code) == 1) {

                            angle2[5] = Map4to51.get(angle1[5]);
                            angle2[8] = Map4to51.get(angle1[8]);
                            angle2[11] = Map4to51.get(angle1[11]);

                        } else if (Integer.parseInt(code) == 2) {

                            angle2[5] = Map4to52.get(angle1[5]);
                            angle2[8] = Map4to52.get(angle1[8]);
                            angle2[11] = Map4to52.get(angle1[11]);

                        } else if (Integer.parseInt(code) == 3) {

                            angle2[5] = Map4to53.get(angle1[5]);
                            angle2[8] = Map4to53.get(angle1[8]);
                            angle2[11] = Map4to53.get(angle1[11]);

                        } else if (Integer.parseInt(code) == 4) {

                            angle2[5] = Map4to54.get(angle1[5]);
                            angle2[8] = Map4to54.get(angle1[8]);
                            angle2[11] = Map4to54.get(angle1[11]);

                        } else if (Integer.parseInt(code) == 5) {

                            angle2[5] = Map4to55.get(angle1[5]);
                            angle2[8] = Map4to55.get(angle1[8]);
                            angle2[11] = Map4to55.get(angle1[11]);

                        } else if (Integer.parseInt(code) == 6) {

                            angle2[5] = Map4to56.get(angle1[5]);
                            angle2[8] = Map4to56.get(angle1[8]);
                            angle2[11] = Map4to56.get(angle1[11]);

                        } else {
                            System.out.println("Not in a patch");
                        }

                        aCol1List.add(angle2[5]);
                        aCol2List.add(angle2[8]);
                        aCol3List.add(angle2[11]);
                        aValue1List.add(angle2[15]);
                        aValue2List.add(angle2[17]);

                        //adds to angles list
                        angles.add("angle     " + angle2[5] + "   " + angle2[8] + "   " + angle2[11] + "    " + angle2[15] + "  " + angle2[17]);

                    } //ends if(patch.contains("angle"))

                }//ends reader while
            }//end try
            catch (IOException e) {

            }
        }//ends j loop

        //define new arrays to enter list data into
        String[] aCol1 = new String[aCol1List.size()];
        String[] aCol2 = new String[aCol2List.size()];
        String[] aCol3 = new String[aCol3List.size()];
        String[] aValue1 = new String[aValue1List.size()];
        String[] aValue2 = new String[aValue2List.size()];

        //convert string lists to string arrays for further use
        aCol1List.toArray(aCol1);
        aCol2List.toArray(aCol2);
        aCol3List.toArray(aCol3);
        aValue1List.toArray(aValue1);
        aValue2List.toArray(aValue2);

        PrintWriter angleWriter1 = new PrintWriter("Oangle.patch");

        for (int i = 0; i < angles.size(); ++i) {

            List<String> col1 = new ArrayList<>();
            List<String> col2 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> col2List1 = new ArrayList<>();

            for (int j = 0; j < angles.size(); ++j) {

                List<Double> col1List = new ArrayList<>();
                List<Double> col2List = new ArrayList<>();

                if (aCol1[i].equals(aCol1[j]) && aCol2[i].equals(aCol2[j]) && aCol3[i].equals(aCol3[j])) {

                    //enter angle values into lists for storage before averaging
                    col1.add(aValue1[j]);
                    col2.add(aValue2[j]);

                    //convert string list to double list for averaging
                    for (String angleN1 : col1) {

                        col1List.add(Double.parseDouble(angleN1));
                    }
                    for (String angleN2 : col2) {
                        col2List.add(Double.parseDouble(angleN2));
                    }

                    col1List1 = col1List;
                    col2List1 = col2List;

                }

            } //ends j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] col2Array = new Double[col2List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            col2List1.toArray(col2Array);

            int size1a = col1Array.length;
            int size2a = col2Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;

            for (int k = 0; k < size1a; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1a;

            for (int k = 0; k < size2a; ++k) {
                sum2 = sum2 + col2Array[k];
            }
            double average2 = sum2 / size2a;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);

            if (!" ".equals(aCol3[i])) {
                angleWriter1.write("angle     " + aCol1[i] + "   " + aCol2[i] + "   " + aCol3[i] + "    " + avg1 + "  " + avg2 + "\n");
            }

        }

        angleWriter1.close();

        //String resultAn = null;
        //StringBuilder resultBuilderan = new StringBuilder();
        Set<String> alreadyPresentan = new HashSet<>();

        //make sure there is only one angle parameter line for each atom combination
        try {
            FileReader fr = new FileReader("Oangle.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder resultBuilder = new StringBuilder();
            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] chars = line.split("(?!^)");
                String[] test1 = new String[25];
                String[] test2 = new String[25];
                String[] test3 = new String[25];
                String[] test4 = new String[25];
                String[] test5 = new String[25];
                String[] test6 = new String[25];

                for (int k = 0; k < 25; ++k) {
                    test1[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test2[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test3[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test4[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test5[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test6[k] = chars[k];
                }

                //switch the angle atom terms for checking purposes
                //if the original line (test1/testS) is: "strbnd     ABC   DEF   GHI"
                //test2/testS2: "strbnd     ABC   GHI   DEF"
                test2[16] = test1[22];
                test2[17] = test1[23];
                test2[18] = test1[24];
                test2[22] = test1[16];
                test2[23] = test1[17];
                test2[24] = test1[18];

                //test3/testS3: "strbnd     DEF   ABC   GHI"
                test3[10] = test1[16];
                test3[11] = test1[17];
                test3[12] = test1[18];
                test3[16] = test1[10];
                test3[17] = test1[11];
                test3[18] = test1[12];

                //test4/testS4: "strbnd     DEF   GHI   ABC"
                test4[10] = test1[16];
                test4[11] = test1[17];
                test4[12] = test1[18];
                test4[16] = test1[22];
                test4[17] = test1[23];
                test4[18] = test1[24];
                test4[22] = test1[10];
                test4[23] = test1[11];
                test4[24] = test1[12];

                //test5/testS5: "strbnd     GHI   ABC   DEF"
                test5[10] = test1[22];
                test5[11] = test1[23];
                test5[12] = test1[24];
                test5[16] = test1[10];
                test5[17] = test1[11];
                test5[18] = test1[12];
                test5[22] = test1[16];
                test5[23] = test1[17];
                test5[24] = test1[18];

                //test6/test6S: "strbnd     GHI   DEF   ABC"
                test6[10] = test1[22];
                test6[11] = test1[23];
                test6[12] = test1[24];
                test6[22] = test1[10];
                test6[23] = test1[11];
                test6[24] = test1[12];

                //converts test1-6 arrays to strings for comparison testing
                String testS = Arrays.toString(test1);
                String testS2 = Arrays.toString(test2);
                String testS3 = Arrays.toString(test3);
                String testS4 = Arrays.toString(test4);
                String testS5 = Arrays.toString(test5);
                String testS6 = Arrays.toString(test6);

                boolean first = true;
                if (!alreadyPresentan.contains(testS) || !alreadyPresentan.contains(testS2) || !alreadyPresentan.contains(testS3)
                        || !alreadyPresentan.contains(testS4) || !alreadyPresentan.contains(testS5) || !alreadyPresentan.contains(testS6)) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilder.append("\n");
                    }

                    if (!alreadyPresentan.contains(testS) || !alreadyPresentan.contains(testS2) || !alreadyPresentan.contains(testS3)
                            || !alreadyPresentan.contains(testS4) || !alreadyPresentan.contains(testS5) || !alreadyPresentan.contains(testS6)) {
                        resultBuilder.append(line);
                        resultBuilder.append("\n");
                    }

                    alreadyPresentan.add(testS);
                    alreadyPresentan.add(testS2);
                    alreadyPresentan.add(testS3);
                    alreadyPresentan.add(testS4);
                    alreadyPresentan.add(testS5);
                    alreadyPresentan.add(testS6);

                }
                String result = resultBuilder.toString();

                //write final angle values to patch file
                try {

                    try (PrintWriter angleWriter2 = new PrintWriter("angle.patch")) {
                        angleWriter2.write(result);
                        angleWriter2.close();
                    }

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Angle Patch: COMPLETE");

        /**
         * *
         * End angle averaging/writing code
         */
        /**
         * *
         * Stretch bend (strbnd) averaging/writing code
         */
        List<String> strbnds = new ArrayList<>();
        List<String> stCol1List = new ArrayList<>();
        List<String> stCol2List = new ArrayList<>();
        List<String> stCol3List = new ArrayList<>();
        List<String> stValue1List = new ArrayList<>();
        List<String> stValue2List = new ArrayList<>();

        for (int j = 0; j < max; ++j) {

            String code = fragCode[j];
            String patchname = patchMap.get(Integer.parseInt(code));
            String strbnd = null;

            try {
                File file = new File(patchname);
                FileReader fileReader = new FileReader(file);
                BufferedReader bufferedReader = new BufferedReader(fileReader);
                StringBuilder stringBuffer = new StringBuilder();
                String line;

                while ((line = bufferedReader.readLine()) != null) {
                    stringBuffer.append(line);
                    stringBuffer.append("\n");

                    String patch = line;

                    if (patch.contains("strbnd")) {
                        strbnd = line;

                        String[] strbnd1 = strbnd.split(" ");
                        String[] strbnd2 = new String[strbnd1.length];

                        for (int k = 0; k < strbnd1.length; ++k) {
                            strbnd2[k] = strbnd1[k];
                        }

                        //replaces fragment atom type numbers with vem. atom type numbers
                        if (Integer.parseInt(code) == 1) {

                            strbnd2[4] = Map4to51.get(strbnd1[4]);
                            strbnd2[7] = Map4to51.get(strbnd1[7]);
                            strbnd2[10] = Map4to51.get(strbnd1[10]);

                        } else if (Integer.parseInt(code) == 2) {

                            strbnd2[4] = Map4to52.get(strbnd1[4]);
                            strbnd2[7] = Map4to52.get(strbnd1[7]);
                            strbnd2[10] = Map4to52.get(strbnd1[10]);

                        } else if (Integer.parseInt(code) == 3) {

                            strbnd2[4] = Map4to53.get(strbnd1[4]);
                            strbnd2[7] = Map4to53.get(strbnd1[7]);
                            strbnd2[10] = Map4to53.get(strbnd1[10]);

                        } else if (Integer.parseInt(code) == 4) {

                            strbnd2[4] = Map4to54.get(strbnd1[4]);
                            strbnd2[7] = Map4to54.get(strbnd1[7]);
                            strbnd2[10] = Map4to54.get(strbnd1[10]);

                        } else if (Integer.parseInt(code) == 5) {

                            strbnd2[4] = Map4to55.get(strbnd1[4]);
                            strbnd2[7] = Map4to55.get(strbnd1[7]);
                            strbnd2[10] = Map4to55.get(strbnd1[10]);

                        } else if (Integer.parseInt(code) == 6) {

                            strbnd2[4] = Map4to56.get(strbnd1[4]);
                            strbnd2[7] = Map4to56.get(strbnd1[7]);
                            strbnd2[10] = Map4to56.get(strbnd1[10]);

                        } else {
                            System.out.println("Not in a patch");
                        }

                        stCol1List.add(strbnd2[4]);
                        stCol2List.add(strbnd2[7]);
                        stCol3List.add(strbnd2[10]);
                        stValue1List.add(strbnd2[14]);
                        stValue2List.add(strbnd2[17]);

                        //adds to strbnds list
                        strbnds.add("strbnd    " + strbnd2[4] + "   " + strbnd2[7] + "   " + strbnd2[10] + "    " + strbnd2[14] + "   " + strbnd2[17]);

                    } //ends if(patch.contains("strbnd"))

                }//ends reader while
            }//end try
            catch (IOException e) {

            }
        }//ends j loop

        //define new arrays to enter list data into
        String[] stCol1 = new String[stCol1List.size()];
        String[] stCol2 = new String[stCol2List.size()];
        String[] stCol3 = new String[stCol3List.size()];
        String[] stValue1 = new String[stValue1List.size()];
        String[] stValue2 = new String[stValue2List.size()];

        //convert string lists to string arrays for further use
        stCol1List.toArray(stCol1);
        stCol2List.toArray(stCol2);
        stCol3List.toArray(stCol3);
        stValue1List.toArray(stValue1);
        stValue2List.toArray(stValue2);

        PrintWriter strbndWriter1 = new PrintWriter("Ostrbnd.patch");

        for (int i = 0; i < strbnds.size(); ++i) {

            List<String> col1 = new ArrayList<>();
            List<String> col2 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> col2List1 = new ArrayList<>();

            for (int j = 0; j < strbnds.size(); ++j) {

                List<Double> col1List = new ArrayList<>();
                List<Double> col2List = new ArrayList<>();

                if (stCol1[i].equals(stCol1[j]) && stCol2[i].equals(stCol2[j]) && stCol3[i].equals(stCol3[j])) {

                    //enter angle values into lists for storage before averaging
                    col1.add(stValue1[j]);
                    col2.add(stValue2[j]);

                    //convert string list to double list for averaging
                    for (String strbndN1 : col1) {

                        col1List.add(Double.parseDouble(strbndN1));
                    }
                    for (String strbndN2 : col2) {
                        col2List.add(Double.parseDouble(strbndN2));
                    }

                    col1List1 = col1List;
                    col2List1 = col2List;

                }

            } //ends j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] col2Array = new Double[col2List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            col2List1.toArray(col2Array);

            int size1st = col1Array.length;
            int size2st = col2Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;

            for (int k = 0; k < size1st; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1st;

            for (int k = 0; k < size2st; ++k) {
                sum2 = sum2 + col2Array[k];
            }
            double average2 = sum2 / size2st;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);

            //writes out strbnd initial values
            if (!" ".equals(stCol3[i])) {
                strbndWriter1.write("strbnd    " + stCol1[i] + "   " + stCol2[i] + "   " + stCol3[i] + "    " + avg1 + "   " + avg2 + "\n");
            }

        }

        strbndWriter1.close();
        Set<String> alreadyPresentst = new HashSet<>();

        //make sure there is only one strbnd parameter line for each atom combination
        try {
            FileReader fr = new FileReader("Ostrbnd.patch");
            StringBuilder resultBuilder = new StringBuilder();
            BufferedReader bufferedReader = new BufferedReader(fr);
            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] chars = line.split("(?!^)");
                String[] test1 = new String[25];
                String[] test2 = new String[25];
                String[] test3 = new String[25];
                String[] test4 = new String[25];
                String[] test5 = new String[25];
                String[] test6 = new String[25];

                for (int k = 0; k < 25; ++k) {
                    test1[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test2[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test3[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test4[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test5[k] = chars[k];
                }
                for (int k = 0; k < 25; ++k) {
                    test6[k] = chars[k];
                }

                //switch the strbnd atom terms for checking purposes
                //if the original line (test1/testS) is: "angle     ABC   DEF   GHI"
                //test2/testS2: "angle     ABC   GHI   DEF"
                test2[16] = test1[22];
                test2[17] = test1[23];
                test2[18] = test1[24];
                test2[22] = test1[16];
                test2[23] = test1[17];
                test2[24] = test1[18];

                //test3/testS3: "angle     DEF   ABC   GHI"
                test3[10] = test1[16];
                test3[11] = test1[17];
                test3[12] = test1[18];
                test3[16] = test1[10];
                test3[17] = test1[11];
                test3[18] = test1[12];

                //test4/testS4: "angle     DEF   GHI   ABC"
                test4[10] = test1[16];
                test4[11] = test1[17];
                test4[12] = test1[18];
                test4[16] = test1[22];
                test4[17] = test1[23];
                test4[18] = test1[24];
                test4[22] = test1[10];
                test4[23] = test1[11];
                test4[24] = test1[12];

                //test5/testS5: "angle     GHI   ABC   DEF"
                test5[10] = test1[22];
                test5[11] = test1[23];
                test5[12] = test1[24];
                test5[16] = test1[10];
                test5[17] = test1[11];
                test5[18] = test1[12];
                test5[22] = test1[16];
                test5[23] = test1[17];
                test5[24] = test1[18];

                //test6/test6S: "angle     GHI   DEF   ABC"
                test6[10] = test1[22];
                test6[11] = test1[23];
                test6[12] = test1[24];
                test6[22] = test1[10];
                test6[23] = test1[11];
                test6[24] = test1[12];

                //converts test1-6 arrays to strings for comparison testing
                String testS = Arrays.toString(test1);
                String testS2 = Arrays.toString(test2);
                String testS3 = Arrays.toString(test3);
                String testS4 = Arrays.toString(test4);
                String testS5 = Arrays.toString(test5);
                String testS6 = Arrays.toString(test6);

                boolean first = true;
                if (!alreadyPresentst.contains(testS) || !alreadyPresentst.contains(testS2) || !alreadyPresentst.contains(testS3)
                        || !alreadyPresentst.contains(testS4) || !alreadyPresentst.contains(testS5) || !alreadyPresentst.contains(testS6)) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilder.append("\n");
                    }

                    if (!alreadyPresentst.contains(testS) || !alreadyPresentst.contains(testS2) || !alreadyPresentst.contains(testS3)
                            || !alreadyPresentst.contains(testS4) || !alreadyPresentst.contains(testS5) || !alreadyPresentst.contains(testS6)) {
                        resultBuilder.append(line);
                        resultBuilder.append("\n");
                    }

                    alreadyPresentst.add(testS);
                    alreadyPresentst.add(testS2);
                    alreadyPresentst.add(testS3);
                    alreadyPresentst.add(testS4);
                    alreadyPresentst.add(testS5);
                    alreadyPresentst.add(testS6);

                }
                String result = resultBuilder.toString();

                //write final strbnd values to patch file
                try {

                    try (PrintWriter strbndWriter2 = new PrintWriter("strbnd.patch")) {
                        strbndWriter2.write(result);
                        strbndWriter2.close();
                    }

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Stretch Bend Patch: COMPLETE");

        /**
         * *
         * End strbnd averaging/writing code
         */
        /**
         * *
         * Out-of-plane bend (opbend) averaging/writing code
         */
        List<String> opbends = new ArrayList<>();
        List<String> opCol1List = new ArrayList<>();
        List<String> opCol2List = new ArrayList<>();
        List<String> opValue1List = new ArrayList<>();
        List<String> opValue2List = new ArrayList<>();

        for (int j = 0; j < max; ++j) {

            String code = fragCode[j];
            String patchname = patchMap.get(Integer.parseInt(code));
            String opbend = null;

            try {
                File file = new File(patchname);
                FileReader fileReader = new FileReader(file);
                BufferedReader bufferedReader = new BufferedReader(fileReader);
                StringBuilder stringBuffer = new StringBuilder();
                String line;

                while ((line = bufferedReader.readLine()) != null) {
                    stringBuffer.append(line);
                    stringBuffer.append("\n");

                    String patch = line;

                    if (patch.contains("opbend")) {
                        opbend = line;

                        String[] opbend1 = (opbend.split(" "));
                        String[] opbend2 = new String[opbend1.length];

                        for (int k = 0; k < opbend1.length; ++k) {
                            opbend2[k] = opbend1[k];
                        }

                        //replaces fragment atom type numbers with vem. atom type numbers
                        if (Integer.parseInt(code) == 1) {

                            opbend2[2] = Map4to51.get(opbend1[2]);
                            opbend2[3] = Map4to51.get(opbend1[3]);

                        } else if (Integer.parseInt(code) == 2) {

                            opbend2[2] = Map4to52.get(opbend1[2]);
                            opbend2[3] = Map4to52.get(opbend1[3]);

                        } else if (Integer.parseInt(code) == 3) {

                            opbend2[2] = Map4to53.get(opbend1[2]);
                            opbend2[3] = Map4to53.get(opbend1[3]);

                        } else if (Integer.parseInt(code) == 4) {

                            opbend2[2] = Map4to54.get(opbend1[2]);
                            opbend2[3] = Map4to54.get(opbend1[3]);

                        } else if (Integer.parseInt(code) == 5) {

                            opbend2[2] = Map4to55.get(opbend1[2]);
                            opbend2[3] = Map4to55.get(opbend1[3]);

                        } else if (Integer.parseInt(code) == 6) {

                            opbend2[2] = Map4to56.get(opbend1[2]);
                            opbend2[3] = Map4to56.get(opbend1[3]);

                        } else {
                            System.out.println("Not in a patch");
                        }

                        opCol1List.add(opbend2[2]);
                        opCol2List.add(opbend2[3]);
                        opValue1List.add(opbend2[10]);

                        //add to opbends list
                        opbends.add("opbend  " + opbend2[2] + "  " + opbend2[3] + " 0   0   " + opbend2[10] + "\n");

                    } //ends if(patch.contains("opbend"))
                } //ends reader while loop

            } catch (IOException e) {
            }
        } //ends j loop

        //define new arrays to enter list data into
        String[] opCol1 = new String[opCol1List.size()];
        String[] opCol2 = new String[opCol2List.size()];
        String[] opValue1 = new String[opValue1List.size()];

        //convert string lists to string arrays for further use
        opCol1List.toArray(opCol1);
        opCol2List.toArray(opCol2);
        opValue1List.toArray(opValue1);

        PrintWriter opWriter1 = new PrintWriter("Oopbend.patch");

        for (int i = 0; i < opbends.size(); ++i) {

            List<String> col1 = new ArrayList<>();
            List<String> col2 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> col2List1 = new ArrayList<>();

            for (int j = i; j < opbends.size(); ++j) {

                List<Double> col1List = new ArrayList<>();
                List<Double> col2List = new ArrayList<>();

                if (opCol1[i].equals(opCol1[j]) && opCol2[i].equals(opCol2[j])) {

                    col1.add(opValue1[j]);

                    //convert string list to double list for averaging
                    for (String opbendN1 : col1) {

                        col1List.add(Double.parseDouble(opbendN1));
                    }

                    col1List1 = col1List;
                    col2List1 = col2List;

                }
            } //ends j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] col2Array = new Double[col2List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            col2List1.toArray(col2Array);

            int size1op = col1Array.length;
            int size2op = col2Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;

            for (int k = 0; k < size1op; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1op;

            for (int k = 0; k < size2op; ++k) {
                sum2 = sum2 + col2Array[k];
            }
            double average2 = sum2 / size2op;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);

            if (!" ".equals(opCol1[i]) && !" ".equals(opCol2[i])) {
                opWriter1.write("opbend  " + opCol1[i] + " " + opCol2[i] + " 0  0  " + avg1 + "   " + "\n");
            }
        } //ends i loop

        opWriter1.close();

        Set<String> alreadyPresentop = new HashSet<>();
        //make sure there is only one opbend parameter line for each atom combination
        try {
            FileReader fr = new FileReader("Oopbend.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder resultBuilder = new StringBuilder();

            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] chars = line.split("(?!^)");
                String[] test1 = new String[15];
                String[] test2 = new String[15];

                for (int k = 0; k < 15; ++k) {
                    test1[k] = chars[k];
                }

                for (int k = 0; k < 15; ++k) {
                    test2[k] = chars[k];
                }

                //switch the opbend atom terms for checking purposes
                test2[8] = test1[12];
                test2[9] = test1[13];
                test2[10] = test1[14];
                test2[12] = test1[8];
                test2[13] = test1[9];
                test2[14] = test1[10];

                String testS = Arrays.toString(test1);
                String testS2 = Arrays.toString(test2);

                boolean first = true;
                if (!alreadyPresentop.contains(testS) || !alreadyPresentop.contains(testS2)) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilder.append("\n");
                    }

                    if (!alreadyPresentop.contains(testS) || !alreadyPresentop.contains(testS2)) {
                        resultBuilder.append(line);
                        resultBuilder.append("\n");
                    }

                    alreadyPresentop.add(testS);
                    alreadyPresentop.add(testS2);

                }
                String result = resultBuilder.toString();

                //write final opbend values to patch file
                try {

                    try (PrintWriter opWriter2 = new PrintWriter("opbend.patch")) {
                        opWriter2.write(result);
                        opWriter2.close();
                    }

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Out-of-Plane Bend Patch: COMPLETE");

        /**
         * *
         * End opbend averaging/writing code
         */
        /**
         * *
         * Torsion averaging/writing code
         */
        List<String> torsions = new ArrayList<>();
        List<String> tCol1List = new ArrayList<>();
        List<String> tCol2List = new ArrayList<>();
        List<String> tCol3List = new ArrayList<>();
        List<String> tCol4List = new ArrayList<>();
        List<String> tValue1List = new ArrayList<>();
        List<String> tValue2List = new ArrayList<>();
        List<String> tValue3List = new ArrayList<>();

        for (int j = 0; j < max; ++j) {

            String code = fragCode[j];
            String patchname = patchMap.get(Integer.parseInt(code));
            String torsion = null;

            try {
                File file = new File(patchname);
                FileReader fileReader = new FileReader(file);
                BufferedReader bufferedReader = new BufferedReader(fileReader);
                StringBuilder stringBuffer = new StringBuilder();
                String line;

                while ((line = bufferedReader.readLine()) != null) {
                    stringBuffer.append(line);
                    stringBuffer.append("\n");

                    String patch = line;

                    if (patch.contains("torsion")) {
                        torsion = line;

                        String[] torsionSplit = torsion.split(" ");

                        //remove spaces to avoid (-) issues
                        List<String> list = new ArrayList<>(Arrays.asList(torsionSplit));
                        list.removeAll(Arrays.asList("", null));

                        String[] torsion1 = new String[list.size()];
                        list.toArray(torsion1);

                        String[] torsion2 = new String[torsion1.length];

                        for (int k = 0; k < torsion1.length; ++k) {
                            torsion2[k] = torsion1[k];
                        }

                        //replaces fragment atom type numbers with vem. atom type numbers
                        if (Integer.parseInt(code) == 1) {

                            torsion2[1] = Map4to51.get(torsion1[1]);
                            torsion2[2] = Map4to51.get(torsion1[2]);
                            torsion2[3] = Map4to51.get(torsion1[3]);
                            torsion2[4] = Map4to51.get(torsion1[4]);

                        } else if (Integer.parseInt(code) == 2) {

                            torsion2[1] = Map4to52.get(torsion1[1]);
                            torsion2[2] = Map4to52.get(torsion1[2]);
                            torsion2[3] = Map4to52.get(torsion1[3]);
                            torsion2[4] = Map4to52.get(torsion1[4]);

                        } else if (Integer.parseInt(code) == 3) {

                            torsion2[1] = Map4to53.get(torsion1[1]);
                            torsion2[2] = Map4to53.get(torsion1[2]);
                            torsion2[3] = Map4to53.get(torsion1[3]);
                            torsion2[4] = Map4to53.get(torsion1[4]);

                        } else if (Integer.parseInt(code) == 4) {

                            torsion2[1] = Map4to54.get(torsion1[1]);
                            torsion2[2] = Map4to54.get(torsion1[2]);
                            torsion2[3] = Map4to54.get(torsion1[3]);
                            torsion2[4] = Map4to54.get(torsion1[4]);

                        } else if (Integer.parseInt(code) == 5) {

                            torsion2[1] = Map4to55.get(torsion1[1]);
                            torsion2[2] = Map4to55.get(torsion1[2]);
                            torsion2[3] = Map4to55.get(torsion1[3]);
                            torsion2[4] = Map4to55.get(torsion1[4]);

                        } else if (Integer.parseInt(code) == 6) {

                            torsion2[1] = Map4to56.get(torsion1[1]);
                            torsion2[2] = Map4to56.get(torsion1[2]);
                            torsion2[3] = Map4to56.get(torsion1[3]);
                            torsion2[4] = Map4to56.get(torsion1[4]);

                        } else {
                            System.out.println("Not in a patch");
                        }

                        tCol1List.add(torsion2[1]);
                        tCol2List.add(torsion2[2]);
                        tCol3List.add(torsion2[3]);
                        tCol4List.add(torsion2[4]);
                        tValue1List.add(torsion2[5]);
                        tValue2List.add(torsion2[8]);
                        tValue3List.add(torsion2[11]);

                        //adds to torsions list
                        torsions.add("torsion   " + torsion2[1] + "   " + torsion2[2] + "   " + torsion2[3] + "   " + torsion2[4] + "     "
                                + torsion2[5] + " 0.0 1     " + torsion2[8] + " 180.0 2     " + torsion2[11] + " 0.0 3");

                    } //ends if(patch.contains("torsion"))

                }//ends reader while
            }//end try
            catch (IOException e) {

            }
        }//ends j loop

        //define new arrays to enter list data into
        String[] tCol1 = new String[tCol1List.size()];
        String[] tCol2 = new String[tCol2List.size()];
        String[] tCol3 = new String[tCol3List.size()];
        String[] tCol4 = new String[tCol4List.size()];
        String[] tValue1 = new String[tValue1List.size()];
        String[] tValue2 = new String[tValue2List.size()];
        String[] tValue3 = new String[tValue3List.size()];

        //convert string lists to string arrays for further use
        tCol1List.toArray(tCol1);
        tCol2List.toArray(tCol2);
        tCol3List.toArray(tCol3);
        tCol4List.toArray(tCol4);
        tValue1List.toArray(tValue1);
        tValue2List.toArray(tValue2);
        tValue3List.toArray(tValue3);

        PrintWriter torsionWriter1 = new PrintWriter("Otorsion.patch");

        for (int i = 0; i < torsions.size(); ++i) {

            List<String> col1 = new ArrayList<>();
            List<String> col2 = new ArrayList<>();
            List<String> col3 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();
            List<Double> col2List1 = new ArrayList<>();
            List<Double> col3List1 = new ArrayList<>();

            for (int j = 0; j < torsions.size(); ++j) {

                List<Double> col1List = new ArrayList<>();
                List<Double> col2List = new ArrayList<>();
                List<Double> col3List = new ArrayList<>();

                if (tCol1[i].equals(tCol1[j]) && tCol2[i].equals(tCol2[j]) && tCol3[i].equals(tCol3[j]) && tCol4[i].equals(tCol4[j])) {

                    //enter torsion values into lists for storage before averaging
                    col1.add(tValue1[j]);
                    col2.add(tValue2[j]);
                    col3.add(tValue3[j]);

                    //convert string list to double list for averaging
                    for (String torsionN1 : col1) {

                        col1List.add(Double.parseDouble(torsionN1));
                    }
                    for (String torsionN2 : col2) {
                        col2List.add(Double.parseDouble(torsionN2));
                    }
                    for (String torsionN3 : col3) {
                        col3List.add(Double.parseDouble(torsionN3));
                    }

                    col1List1 = col1List;
                    col2List1 = col2List;
                    col3List1 = col3List;

                }

            } //ends j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];
            Double[] col2Array = new Double[col2List1.size()];
            Double[] col3Array = new Double[col3List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);
            col2List1.toArray(col2Array);
            col3List1.toArray(col3Array);

            int size1t = col1Array.length;
            int size2t = col2Array.length;
            int size3t = col3Array.length;

            //average variables and function
            double sum1 = 0;
            double sum2 = 0;
            double sum3 = 0;

            for (int k = 0; k < size1t; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1t;

            for (int k = 0; k < size2t; ++k) {
                sum2 = sum2 + col2Array[k];
            }
            double average2 = sum2 / size2t;

            for (int k = 0; k < size3t; ++k) {
                sum3 = sum3 + col3Array[k];
            }
            double average3 = sum3 / size3t;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);
            String avg2 = avgFormat.format(average2);
            String avg3 = avgFormat.format(average3);

            //writes out torsion initial values
            if (!" ".equals(tCol1[i]) && !" ".equals(tCol2[i]) && !" ".equals(tCol3[i]) && !" ".equals(tCol4[i])) {

                torsionWriter1.write("torsion   " + tCol1[i] + "   " + tCol2[i] + "   " + tCol3[i] + "   " + tCol4[i] + "     "
                        + avg1 + " 0.0 1     " + avg2 + " 180.0 2     " + avg3 + " 0.0 3" + "\n");
            }

        }

        torsionWriter1.close();

        Set<String> alreadyPresentt = new HashSet<>();
        //make sure there is only one torsion parameter line for each atom combination
        try {
            FileReader fr = new FileReader("Otorsion.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder resultBuilder = new StringBuilder();
            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] chars = line.split(" ");

                //remove spaces to avoid (-) issues
                List<String> listC = new ArrayList<>(Arrays.asList(chars));
                listC.removeAll(Arrays.asList("", null));

                String[] testC = new String[listC.size()];
                listC.toArray(testC);

                //create comparison/uniqueness testing strings
                String testS = ("torsion   " + testC[1] + "   " + testC[2] + "   " + testC[3] + "   " + testC[4]);
                String testS2 = ("torsion   " + testC[1] + "   " + testC[2] + "   " + testC[4] + "   " + testC[3]);
                String testS3 = ("torsion   " + testC[1] + "   " + testC[3] + "   " + testC[2] + "   " + testC[4]);
                String testS4 = ("torsion   " + testC[1] + "   " + testC[3] + "   " + testC[4] + "   " + testC[2]);
                String testS5 = ("torsion   " + testC[1] + "   " + testC[4] + "   " + testC[2] + "   " + testC[3]);
                String testS6 = ("torsion   " + testC[1] + "   " + testC[4] + "   " + testC[3] + "   " + testC[2]);

                String testS7 = ("torsion   " + testC[2] + "   " + testC[1] + "   " + testC[3] + "   " + testC[4]);
                String testS8 = ("torsion   " + testC[2] + "   " + testC[1] + "   " + testC[4] + "   " + testC[3]);
                String testS9 = ("torsion   " + testC[2] + "   " + testC[3] + "   " + testC[1] + "   " + testC[4]);
                String testS10 = ("torsion   " + testC[2] + "   " + testC[3] + "   " + testC[4] + "   " + testC[1]);
                String testS11 = ("torsion   " + testC[2] + "   " + testC[4] + "   " + testC[1] + "   " + testC[3]);
                String testS12 = ("torsion   " + testC[2] + "   " + testC[4] + "   " + testC[3] + "   " + testC[1]);

                String testS13 = ("torsion   " + testC[3] + "   " + testC[1] + "   " + testC[2] + "   " + testC[4]);
                String testS14 = ("torsion   " + testC[3] + "   " + testC[1] + "   " + testC[4] + "   " + testC[2]);
                String testS15 = ("torsion   " + testC[3] + "   " + testC[2] + "   " + testC[1] + "   " + testC[4]);
                String testS16 = ("torsion   " + testC[3] + "   " + testC[2] + "   " + testC[4] + "   " + testC[1]);
                String testS17 = ("torsion   " + testC[3] + "   " + testC[4] + "   " + testC[1] + "   " + testC[2]);
                String testS18 = ("torsion   " + testC[3] + "   " + testC[4] + "   " + testC[2] + "   " + testC[1]);

                String testS19 = ("torsion   " + testC[4] + "   " + testC[1] + "   " + testC[2] + "   " + testC[3]);
                String testS20 = ("torsion   " + testC[4] + "   " + testC[1] + "   " + testC[3] + "   " + testC[2]);
                String testS21 = ("torsion   " + testC[4] + "   " + testC[2] + "   " + testC[1] + "   " + testC[3]);
                String testS22 = ("torsion   " + testC[4] + "   " + testC[2] + "   " + testC[3] + "   " + testC[1]);
                String testS23 = ("torsion   " + testC[4] + "   " + testC[3] + "   " + testC[1] + "   " + testC[2]);
                String testS24 = ("torsion   " + testC[4] + "   " + testC[3] + "   " + testC[2] + "   " + testC[1]);

                boolean first = true;
                if (!alreadyPresentt.contains(testS) || !alreadyPresentt.contains(testS2) || !alreadyPresentt.contains(testS3)
                        || !alreadyPresentt.contains(testS4) || !alreadyPresentt.contains(testS5) || !alreadyPresentt.contains(testS6)
                        || !alreadyPresentt.contains(testS7) || !alreadyPresentt.contains(testS8) || !alreadyPresentt.contains(testS9)
                        || !alreadyPresentt.contains(testS10) || !alreadyPresentt.contains(testS11) || !alreadyPresentt.contains(testS12)
                        || !alreadyPresentt.contains(testS13) || !alreadyPresentt.contains(testS14) || !alreadyPresentt.contains(testS15)
                        || !alreadyPresentt.contains(testS16) || !alreadyPresentt.contains(testS17) || !alreadyPresentt.contains(testS18)
                        || !alreadyPresentt.contains(testS19) || !alreadyPresentt.contains(testS20) || !alreadyPresentt.contains(testS21)
                        || !alreadyPresentt.contains(testS22) || !alreadyPresentt.contains(testS23) || !alreadyPresentt.contains(testS24)) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilder.append("\n");
                    }

                    if (!alreadyPresentt.contains(testS) || !alreadyPresentt.contains(testS2) || !alreadyPresentt.contains(testS3)
                            || !alreadyPresentt.contains(testS4) || !alreadyPresentt.contains(testS5) || !alreadyPresentt.contains(testS6)
                            || !alreadyPresentt.contains(testS7) || !alreadyPresentt.contains(testS8) || !alreadyPresentt.contains(testS9)
                            || !alreadyPresentt.contains(testS10) || !alreadyPresentt.contains(testS11) || !alreadyPresentt.contains(testS12)
                            || !alreadyPresentt.contains(testS13) || !alreadyPresentt.contains(testS14) || !alreadyPresentt.contains(testS15)
                            || !alreadyPresentt.contains(testS16) || !alreadyPresentt.contains(testS17) || !alreadyPresentt.contains(testS18)
                            || !alreadyPresentt.contains(testS19) || !alreadyPresentt.contains(testS20) || !alreadyPresentt.contains(testS21)
                            || !alreadyPresentt.contains(testS22) || !alreadyPresentt.contains(testS23) || !alreadyPresentt.contains(testS24)) {
                        resultBuilder.append(line);
                        resultBuilder.append("\n");
                    }

                    alreadyPresentt.add(testS);
                    alreadyPresentt.add(testS2);
                    alreadyPresentt.add(testS3);
                    alreadyPresentt.add(testS4);
                    alreadyPresentt.add(testS5);
                    alreadyPresentt.add(testS6);
                    alreadyPresentt.add(testS7);
                    alreadyPresentt.add(testS8);
                    alreadyPresentt.add(testS9);
                    alreadyPresentt.add(testS10);
                    alreadyPresentt.add(testS11);
                    alreadyPresentt.add(testS12);
                    alreadyPresentt.add(testS13);
                    alreadyPresentt.add(testS14);
                    alreadyPresentt.add(testS15);
                    alreadyPresentt.add(testS16);
                    alreadyPresentt.add(testS17);
                    alreadyPresentt.add(testS18);
                    alreadyPresentt.add(testS19);
                    alreadyPresentt.add(testS20);
                    alreadyPresentt.add(testS21);
                    alreadyPresentt.add(testS22);
                    alreadyPresentt.add(testS23);
                    alreadyPresentt.add(testS24);

                }
                String result = resultBuilder.toString();

                //write final torsion values to patch file
                try {

                    try (PrintWriter torsionWriter2 = new PrintWriter("torsion.patch")) {
                        torsionWriter2.write(result);
                        torsionWriter2.close();
                    }

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Torsion Patch: COMPLETE");

        /**
         * *
         * End torsion averaging/writing code
         */
        /**
         * *
         * Pitors averaging/writing code
         */
        List<String> pitorss = new ArrayList<>();
        List<String> piCol1List = new ArrayList<>();
        List<String> piCol2List = new ArrayList<>();
        List<String> piValue1List = new ArrayList<>();
        List<String> piValue2List = new ArrayList<>();

        for (int j = 0; j < max; ++j) {

            String code = fragCode[j];
            String patchname = patchMap.get(Integer.parseInt(code));
            String pitors = null;

            try {
                File file = new File(patchname);
                FileReader fileReader = new FileReader(file);
                BufferedReader bufferedReader = new BufferedReader(fileReader);
                StringBuilder stringBuffer = new StringBuilder();
                String line;

                while ((line = bufferedReader.readLine()) != null) {
                    stringBuffer.append(line);
                    stringBuffer.append("\n");

                    String patch = line;

                    if (patch.contains("pitors")) {
                        pitors = line;

                        String[] pitors1 = (pitors.split(" "));
                        String[] pitors2 = new String[pitors1.length];

                        for (int k = 0; k < pitors1.length; ++k) {
                            pitors2[k] = pitors1[k];
                        }

                        //replaces fragment atom type numbers with vem. atom type numbers
                        if (Integer.parseInt(code) == 1) {

                            pitors2[4] = Map4to51.get(pitors1[4]);
                            pitors2[7] = Map4to51.get(pitors1[7]);

                        } else if (Integer.parseInt(code) == 2) {

                            pitors2[4] = Map4to52.get(pitors1[4]);
                            pitors2[7] = Map4to52.get(pitors1[7]);

                        } else if (Integer.parseInt(code) == 3) {

                            pitors2[4] = Map4to53.get(pitors1[4]);
                            pitors2[7] = Map4to53.get(pitors1[7]);

                        } else if (Integer.parseInt(code) == 4) {

                            pitors2[4] = Map4to54.get(pitors1[4]);
                            pitors2[7] = Map4to54.get(pitors1[7]);

                        } else if (Integer.parseInt(code) == 5) {

                            pitors2[4] = Map4to55.get(pitors1[4]);
                            pitors2[7] = Map4to55.get(pitors1[7]);

                        } else if (Integer.parseInt(code) == 6) {

                            pitors2[4] = Map4to56.get(pitors1[4]);
                            pitors2[7] = Map4to56.get(pitors1[7]);

                        } else {
                            System.out.println("Not in a patch");
                        }

                        piCol1List.add(pitors2[4]);
                        piCol2List.add(pitors2[7]);
                        piValue1List.add(pitors2[12]);

                        //add to pitorss list
                        pitorss.add("pitors    " + pitors2[4] + "   " + pitors2[7] + "     " + pitors2[12] + "\n");

                    } //ends if(patch.contains("pitors"))
                } //ends reader while loop

            } catch (IOException e) {
            }
        } //ends j loop

        //define new arrays to enter list data into
        String[] piCol1 = new String[piCol1List.size()];
        String[] piCol2 = new String[piCol2List.size()];
        String[] piValue1 = new String[piValue1List.size()];

        //convert string lists to string arrays for further use
        piCol1List.toArray(piCol1);
        piCol2List.toArray(piCol2);
        piValue1List.toArray(piValue1);

        PrintWriter piWriter1 = new PrintWriter("Opitors.patch");

        for (int i = 0; i < pitorss.size(); ++i) {

            List<String> col1 = new ArrayList<>();

            List<Double> col1List1 = new ArrayList<>();

            for (int j = i; j < pitorss.size(); ++j) {

                List<Double> col1List = new ArrayList<>();

                if (piCol1[i].equals(piCol1[j]) && piCol2[i].equals(piCol2[j])) {

                    col1.add(piValue1[j]);

                    //convert string list to double list for averaging
                    for (String pitorsN1 : col1) {

                        col1List.add(Double.parseDouble(pitorsN1));
                    }

                    col1List1 = col1List;

                }
            } //ends j loop

            //Define new array for averaging
            Double[] col1Array = new Double[col1List1.size()];

            //convert double list to double array for averaging
            col1List1.toArray(col1Array);

            int size1pi = col1Array.length;

            //average variables and function
            double sum1 = 0;

            for (int k = 0; k < size1pi; ++k) {
                sum1 = sum1 + col1Array[k];
            }
            double average1 = sum1 / size1pi;

            DecimalFormat avgFormat = new DecimalFormat("0.0000");

            String avg1 = avgFormat.format(average1);

            if (!" ".equals(piCol1[i]) && !" ".equals(piCol2[i])) {
                piWriter1.write("pitors    " + piCol1[i] + "   " + piCol2[i] + "     " + avg1 + "\n");
            }
        } //ends i "for" loop

        piWriter1.close();

        Set<String> alreadyPresentpi = new HashSet<>();

        //make sure there is only one pitors parameter line for each atom combination
        try {
            FileReader fr = new FileReader("Opitors.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder resultBuilder = new StringBuilder();
            String line;
            while ((line = bufferedReader.readLine()) != null) {

                String[] chars = line.split("(?!^)");
                String[] test1 = new String[19];
                String[] test2 = new String[19];

                for (int k = 0; k < 19; ++k) {
                    test1[k] = chars[k];
                }

                for (int k = 0; k < 19; ++k) {
                    test2[k] = chars[k];
                }

                //switch the pitors atom type terms for checking purposes
                test2[10] = test1[16];
                test2[11] = test1[17];
                test2[12] = test1[18];
                test2[16] = test1[10];
                test2[17] = test1[11];
                test2[18] = test1[12];

                String testS = Arrays.toString(test1);
                String testS2 = Arrays.toString(test2);

                boolean first = true;
                if (!alreadyPresentpi.contains(testS) || !alreadyPresentpi.contains(testS2)) {
                    if (first) {
                        first = false;
                    } else {
                        resultBuilder.append("\n");
                    }

                    if (!alreadyPresentpi.contains(testS) || !alreadyPresentpi.contains(testS2)) {
                        resultBuilder.append(line);
                        resultBuilder.append("\n");
                    }

                    alreadyPresentpi.add(testS);
                    alreadyPresentpi.add(testS2);

                }
                String result = resultBuilder.toString();

                //write final pitors values to patch file
                try {

                    try (PrintWriter piWriter2 = new PrintWriter("pitors.patch")) {
                        piWriter2.write(result);
                        piWriter2.close();
                    }

                } catch (IOException e) {
                }

            }

        } catch (IOException e) {
        }

        System.out.println("Pitors Patch: COMPLETE");

        /**
         * *
         * End pitors averaging/writing code
         */
        /**
         * *
         * Final writing code
         */
        //read each file
        //print contents to final.patch
        PrintWriter finalWriter = new PrintWriter("final.patch");

        try {
            FileReader fr = new FileReader("pitors.patch");
            BufferedReader bufferedReader = new BufferedReader(fr);
            StringBuilder stringBuffer = new StringBuilder();
            String line;
            while ((line = bufferedReader.readLine()) != null) {
                finalWriter.write(line + "\n");
            }

            finalWriter.close();

        } catch (IOException e) {

        }

        System.out.println("Final Patch: COMPLETE");
    }

}
