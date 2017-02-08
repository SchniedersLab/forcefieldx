/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;

import ffx.autoparm.PME_2.Polarization;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;

import static ffx.autoparm.PME_2.atoms;
import static ffx.autoparm.PME_2.axisAtom;
import static ffx.autoparm.PME_2.coordinates;
import static ffx.autoparm.PME_2.directDipole;
import static ffx.autoparm.PME_2.directDipolep;
import static ffx.autoparm.PME_2.field1;
import static ffx.autoparm.PME_2.field2;
import static ffx.autoparm.PME_2.frame;
import static ffx.autoparm.PME_2.globalMultipole;
import static ffx.autoparm.PME_2.induce_pedit;
import static ffx.autoparm.PME_2.inducedDipole;
import static ffx.autoparm.PME_2.inducedDipolep;
import static ffx.autoparm.PME_2.ip11;
import static ffx.autoparm.PME_2.ip12;
import static ffx.autoparm.PME_2.ip13;
import static ffx.autoparm.PME_2.ip14;
import static ffx.autoparm.PME_2.localMultipole;
import static ffx.autoparm.PME_2.maxThreads;
import static ffx.autoparm.PME_2.nAtoms;
import static ffx.autoparm.PME_2.nSymm;
import static ffx.autoparm.PME_2.neighborLists;
import static ffx.autoparm.PME_2.off2;
import static ffx.autoparm.PME_2.parallelTeam;
import static ffx.autoparm.PME_2.pdamp;
import static ffx.autoparm.PME_2.pedit;
import static ffx.autoparm.PME_2.polargrp;
import static ffx.autoparm.PME_2.polarizability;
import static ffx.autoparm.PME_2.polarization;
import static ffx.autoparm.PME_2.poleps;
import static ffx.autoparm.PME_2.polsor;
import static ffx.autoparm.PME_2.rotateMultipolesRegion;
import static ffx.autoparm.PME_2.thole;
import static ffx.potential.parameters.MultipoleType.BOHR;

/**
 * Poledit Provides Multipole Parameters from GDMA Output.
 *
 * Possible issue: xyzIndex and Type are the same always.
 *
 * @author Gaurav Chattree
 * @since 1.0
 *
 */
public class Poledit {

    private ArrayList<Atom> atomslist = new ArrayList<Atom>();
    private ArrayList<Double> polaritylist = new ArrayList<Double>();
    private ArrayList<Double> pdamplist = new ArrayList<Double>();
    private int[] xaxis;
    private static final Logger logger = Logger.getLogger(Poledit.class.getName());
    private boolean remove_symmetry = false;

    /**
     * This method takes data from GDMA and prints out multipole parameters
     *
     * @param gdmaoutfname File location of multipole params output by GDMA
     * @param peditinfname File location of molecular polarization group
     * information
     */
    public Poledit(String gdmaoutfname, String peditinfname) {
        pedit = true;
        nSymm = 1;
        readGDMA(gdmaoutfname);
        int index = gdmaoutfname.lastIndexOf(".");
        String name = gdmaoutfname.substring(0, index);
        setup_print_xyz(name);
        parallelTeam = new ParallelTeam();
        maxThreads = parallelTeam.getThreadCount();

        nAtoms = atomslist.size();
        ip11 = new int[nAtoms][];
        ip12 = new int[nAtoms][];
        ip13 = new int[nAtoms][];
        ip14 = new int[nAtoms][];
        field1 = new double[nAtoms][3];
        field2 = new double[nAtoms][3];
        inducedDipole = new double[1][nAtoms][3];
        inducedDipolep = new double[1][nAtoms][3];
        directDipole = new double[nAtoms][3];
        directDipolep = new double[nAtoms][3];
        neighborLists = new int[1][nAtoms][3];
        poleps = 1e-6;
        polsor = .7;
        off2 = 49;

        /*set local frames.*/
        setframe(peditinfname);
        printglobalmpoles();
        /*Go from global to local since: inverse = true*/
        rotateMultipolesRegion = new PME_2.RotateMultipolesRegion(maxThreads, true);

        try {
            parallelTeam.execute(rotateMultipolesRegion);
        } catch (Exception ex) {
            Logger.getLogger(Poledit.class.getName()).log(Level.SEVERE, null, ex);
        }
        printlocalmpoles();
        /*create polarization groups*/
        polargrp();
        polarization = Polarization.MUTUAL;

        /*find induced values and remove them from the multipoles*/
        induce_pedit();
        printlocalmpoles();

        removeInducedFromGlobal();

        /*Go from global to local since: inverse = true*/
        try {
            parallelTeam.execute(rotateMultipolesRegion);
        } catch (Exception ex) {
            Logger.getLogger(Poledit.class.getName()).log(Level.SEVERE, null, ex);
        }

        printlocalmpoles();
        fixpolar();
        printlocalmpoles();

        prtpolar(name);
    }

    /**
     * <p>
     * posinarray</p>
     *
     * @param find a {@link java.lang.String} object.
     * @param ln a {@link java.lang.String} object.
     * @return a int.
     */
    public int posinarray(String find, String ln) {
        int pos = -1;

        for (int i = 0; i < ln.split(" +").length; i++) {
            if (ln.split(" +")[i].equals(find)) {
                pos = i;
                break;
            }
        }

        return pos;
    }

    /**
     * <p>
     * readGDMA</p>
     *
     * @param gdmaoutfname a {@link java.lang.String} object.
     */
    public void readGDMA(String gdmaoutfname) {
        ArrayList<double[]> sphmultipole = new ArrayList<double[]>();
        File gdmaoutf = new File(gdmaoutfname);
        try {
            if (gdmaoutf != null && gdmaoutf.exists() && gdmaoutf.canRead()) {
                BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(gdmaoutf)));
                String line;
                int i = 1;
                String element;
                while ((line = br.readLine()) != null) {
                    if (line.contains("x =") && !line.contains("origin")) {
                        element = line.split(" +")[0];
                        double xyz[] = new double[3];
                        xyz[0] = Double.parseDouble(line.split(" +")[3]);
                        xyz[1] = Double.parseDouble(line.split(" +")[6]);
                        xyz[2] = Double.parseDouble(line.split(" +")[9]);
                        line = br.readLine();
                        //radius = Double.parseDouble(line.split(" +")[8]);
                        double mp[] = new double[9];
                        for (int t = 0; t < mp.length; t++) {
                            mp[t] = 0;
                        }
                        line = br.readLine();
                        line = line + " " + br.readLine();
                        line = line + " " + br.readLine();
                        line = line + " " + br.readLine();

                        if (line.contains("Q00")) {
                            int x = posinarray("Q00", line);
                            mp[0] = Double.parseDouble(line.split(" +")[x + 2]); //Q00
                        }
                        if (line.contains("Q10")) {
                            int x = posinarray("Q10", line);
                            mp[1] = Double.parseDouble(line.split(" +")[x + 2]); //Q10
                        }
                        if (line.contains("Q11c")) {
                            int x = posinarray("Q11c", line);
                            mp[2] = Double.parseDouble(line.split(" +")[x + 2]); //Q11c
                        }
                        if (line.contains("Q11s")) {
                            int x = posinarray("Q11s", line);
                            mp[3] = Double.parseDouble(line.split(" +")[x + 2]); //Q11s
                        }
                        if (line.contains("Q20")) {
                            int x = posinarray("Q20", line);
                            mp[4] = Double.parseDouble(line.split(" +")[x + 2]); //Q20
                        }
                        if (line.contains("Q21c")) {
                            int x = posinarray("Q21c", line);
                            mp[5] = Double.parseDouble(line.split(" +")[x + 2]); //Q21c
                        }
                        if (line.contains("Q21s")) {
                            int x = posinarray("Q21s", line);
                            mp[6] = Double.parseDouble(line.split(" +")[x + 2]); //Q21s
                        }
                        if (line.contains("Q22c")) {
                            int x = posinarray("Q22c", line);
                            mp[7] = Double.parseDouble(line.split(" +")[x + 2]); //Q22c
                        }
                        if (line.contains("Q22s")) {
                            int x = posinarray("Q22s", line);
                            mp[8] = Double.parseDouble(line.split(" +")[x + 2]); //Q22s
                        }

//                        if(line.contains("Q00")){
//                            mp[0] = Double.parseDouble(line.split(" +")[3]); //Q00
//                        }
//                        line = br.readLine();
//                        if(line.contains("Q10") && line.contains("Q11c") && line.contains("Q11s")){
//                            mp[1] = Double.parseDouble(line.split(" +")[5]); //Q10
//                            mp[2] = Double.parseDouble(line.split(" +")[8]); //Q11c
//                            mp[3] = Double.parseDouble(line.split(" +")[11]); //Q11s
//                        }
//                        else if(line.contains("Q10") && !line.contains("Q11c") && !line.contains("Q11s")){
//                            mp[1] = Double.parseDouble(line.split(" +")[5]); //Q10
//                        }
//                        else if(!line.contains("Q10") && line.contains("Q11c") && !line.contains("Q11s")){
//                            mp[2] = Double.parseDouble(line.split(" +")[5]); //Q10
//                        }
//                        else if(!line.contains("Q10") && !line.contains("Q11c") && line.contains("Q11s")){
//                            mp[3] = Double.parseDouble(line.split(" +")[5]); //Q10
//                        }
//                        else if(!line.contains("Q10") && line.contains("Q11c") && line.contains("Q11s")){
//                            mp[2] = Double.parseDouble(line.split(" +")[5]); //Q11c
//                            mp[3] = Double.parseDouble(line.split(" +")[8]); //Q11s
//                        }
//                        else if(line.contains("Q10") && line.contains("Q11c") && !line.contains("Q11s")){
//                            mp[1] = Double.parseDouble(line.split(" +")[5]); //Q10
//                            mp[2] = Double.parseDouble(line.split(" +")[8]); //Q11c
//                        }
//                        else if(line.contains("Q10") && !line.contains("Q11c") && line.contains("Q11s")){
//                            mp[1] = Double.parseDouble(line.split(" +")[5]); //Q10
//                            mp[3] = Double.parseDouble(line.split(" +")[8]); //Q11s
//                        }
//                        line = br.readLine();
//                        if(line.contains("Q20") && line.contains("Q21c") && line.contains("Q21s")){
//                            mp[4] = Double.parseDouble(line.split(" +")[5]); //Q20
//                            mp[5] = Double.parseDouble(line.split(" +")[8]); //Q21c
//                            mp[6] = Double.parseDouble(line.split(" +")[11]); //Q21s
//                        }
//                        else if(line.contains("Q20") && !line.contains("Q21c") && !line.contains("Q21s")){
//                            mp[4] = Double.parseDouble(line.split(" +")[5]); //Q20
//                        }
//                        else if(!line.contains("Q20") && line.contains("Q21c") && !line.contains("Q21s")){
//                            mp[5] = Double.parseDouble(line.split(" +")[5]); //Q21c
//                        }
//                        else if(!line.contains("Q20") && !line.contains("Q21c") && line.contains("Q21s")){
//                            mp[6] = Double.parseDouble(line.split(" +")[5]); //Q21s
//                        }
//                        else if(!line.contains("Q20") && line.contains("Q21c") && line.contains("Q21s")){
//                            mp[5] = Double.parseDouble(line.split(" +")[5]); //Q21c
//                            mp[6] = Double.parseDouble(line.split(" +")[8]); //Q21s
//                        }
//                        else if(line.contains("Q20") && line.contains("Q21c") && !line.contains("Q21s")){
//                            mp[4] = Double.parseDouble(line.split(" +")[5]); //Q20
//                            mp[5] = Double.parseDouble(line.split(" +")[8]); //Q21c
//                        }
//                        else if(line.contains("Q20") && !line.contains("Q21c") && line.contains("Q21s")){
//                            mp[4] = Double.parseDouble(line.split(" +")[5]); //Q20
//                            mp[6] = Double.parseDouble(line.split(" +")[8]); //Q21s
//                        }
//                        else if(line.contains("Q20") && !line.contains("Q21c") && line.contains("Q21s")){
//                            mp[4] = Double.parseDouble(line.split(" +")[5]); //Q20
//                            mp[6] = Double.parseDouble(line.split(" +")[8]); //Q21s
//                        }
//                        else if(line.contains("Q20") && !line.contains("Q21c") && line.contains("Q21s")){
//                            mp[4] = Double.parseDouble(line.split(" +")[5]); //Q20
//                            mp[6] = Double.parseDouble(line.split(" +")[8]); //Q21s
//                        }
//                        line = br.readLine();
//                        if(line.contains("Q22c") && line.contains("Q22s")){
//                            mp[7] = Double.parseDouble(line.split(" +")[3]); //Q22c
//                            mp[8] = Double.parseDouble(line.split(" +")[6]); //Q22s
//                        }
//                        else if(line.contains("Q22c") && !line.contains("Q22s")){
//                            mp[7] = Double.parseDouble(line.split(" +")[3]); //Q22s
//                        }
//                        else if(!line.contains("Q22c") && line.contains("Q22s")){
//                            mp[8] = Double.parseDouble(line.split(" +")[3]); //Q22s
//                        }
                        double radius[] = {0};
                        double polart[] = {0};
                        double pd[] = {0};
                        AtomType at = get_atom_type(i, element, radius, polart, pd);

                        Atom a = new Atom(Integer.toString(i));
                        a.setXYZ(xyz);
                        a.setAtomType(at);
                        a.setBornRadius(radius[0]);
                        a.setXyzIndex(i);

                        //What is this? CHECK
                        polaritylist.add(polart[0]);
                        /////

                        pdamplist.add(pd[0]);
                        atomslist.add(a);
                        sphmultipole.add(mp);
                        i++;
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Error Reading File");
            System.exit(1);
        }

        //Initialize various arrays0
        atoms = new Atom[atomslist.size()];
        for (int i = 0; i < atomslist.size(); i++) {
            atoms[i] = atomslist.get(i);
        }
        coordinates = new double[nSymm][3][atomslist.size()];
        for (int i = 0; i < atomslist.size(); i++) {
            coordinates[0][0][i] = atoms[i].getX();
            coordinates[0][1][i] = atoms[i].getY();
            coordinates[0][2][i] = atoms[i].getZ();
        }
        pdamp = new double[atomslist.size()];
        for (int i = 0; i < atomslist.size(); i++) {
            pdamp[i] = pdamplist.get(i);
        }

        thole = new double[atomslist.size()];
        for (int i = 0; i < atomslist.size(); i++) {
            thole[i] = .39;
        }

        polarizability = new double[atomslist.size()];
        for (int i = 0; i < atomslist.size(); i++) {
            polarizability[i] = polaritylist.get(i);
        }

        //Convert multipoles from spherical harmonics to cartesian coordinates
        //should this be global or local?//
        globalMultipole = new double[1][atoms.length][13];
        localMultipole = new double[atoms.length][13];
        double term = Math.sqrt(.75);
        for (int i = 0; i < atoms.length; i++) {
            globalMultipole[0][i][0] = sphmultipole.get(i)[0]; //Q00
            globalMultipole[0][i][1] = sphmultipole.get(i)[2]; //Q11c
            globalMultipole[0][i][2] = sphmultipole.get(i)[3]; //Q11s
            globalMultipole[0][i][3] = sphmultipole.get(i)[1]; //Q10
            globalMultipole[0][i][4] = -.5 * sphmultipole.get(i)[4] + term * sphmultipole.get(i)[7]; //-.5 * Q20 + term * Q22c
            globalMultipole[0][i][5] = term * sphmultipole.get(i)[8]; //term * Q22s
            globalMultipole[0][i][6] = term * sphmultipole.get(i)[5]; //term * Q21c
            globalMultipole[0][i][7] = term * sphmultipole.get(i)[8]; //term * Q22s (same as 5)
            globalMultipole[0][i][8] = -.5 * sphmultipole.get(i)[4] - term * sphmultipole.get(i)[7]; //-.5 * Q20 - term * Q22c
            globalMultipole[0][i][9] = term * sphmultipole.get(i)[6]; // term * Q21s
            globalMultipole[0][i][10] = term * sphmultipole.get(i)[5]; //term * Q21c (same as 6)
            globalMultipole[0][i][11] = term * sphmultipole.get(i)[6]; // term * Q21s (same as 9)
            globalMultipole[0][i][12] = sphmultipole.get(i)[4]; //Q20
        }

        //Maintain traceless values and convert units?
        for (int i = 0; i < atoms.length; i++) {
            for (int j = 1; j < 4; j++) {
                globalMultipole[0][i][j] = globalMultipole[0][i][j] * BOHR;
            }
            for (int j = 4; j < 13; j++) {
                globalMultipole[0][i][j] = globalMultipole[0][i][j] * Math.pow(BOHR, 2) / 3;
            }
        }
    }

    /**
     * Creates and returns an atomType object from an element name
     *
     * To work on: What is 'thole' and is it necessary in this function.
     *
     * @param type a int.
     * @param element a {@link java.lang.String} object.
     * @param radius an array of double.
     * @param pdamp an array of double.
     * @return atomType
     * @param polarity an array of double.
     */
    public AtomType get_atom_type(int type, String element, double[] radius,
            double[] polarity, double[] pdamp) {
        int atmnum = 0;
        double atmmass;
        double rad, pol, pd;

        // Get atomic number from element name
        if (element.equalsIgnoreCase("SI")) {
            atmnum = 14;
        } else if (element.equalsIgnoreCase("CL")) {
            atmnum = 17;
        } else if (element.equalsIgnoreCase("BR")) {
            atmnum = 35;
        } else if (element.equalsIgnoreCase("H")) {
            atmnum = 1;
        } else if (element.equalsIgnoreCase("B")) {
            atmnum = 5;
        } else if (element.equalsIgnoreCase("C")) {
            atmnum = 6;
        } else if (element.equalsIgnoreCase("N")) {
            atmnum = 7;
        } else if (element.equalsIgnoreCase("O")) {
            atmnum = 8;
        } else if (element.equalsIgnoreCase("F")) {
            atmnum = 9;
        } else if (element.equalsIgnoreCase("P")) {
            atmnum = 15;
        } else if (element.equalsIgnoreCase("S")) {
            atmnum = 16;
        } else if (element.equalsIgnoreCase("I")) {
            atmnum = 53;
        }
        // Set atomic radius
        radius[0] = 0.77;
        if (atmnum == 0) {
            radius[0] = 0.00;
        }
        if (atmnum == 1) {
            radius[0] = 0.37;
        }
        if (atmnum == 2) {
            radius[0] = 0.32;
        }
        if (atmnum == 6) {
            radius[0] = 0.77;
        }
        if (atmnum == 7) {
            radius[0] = 0.75;
        }
        if (atmnum == 8) {
            radius[0] = 0.73;
        }
        if (atmnum == 9) {
            radius[0] = 0.71;
        }
        if (atmnum == 10) {
            radius[0] = 0.69;
        }
        if (atmnum == 14) {
            radius[0] = 1.11;
        }
        if (atmnum == 15) {
            radius[0] = 1.06;
        }
        if (atmnum == 16) {
            radius[0] = 1.02;
        }
        if (atmnum == 17) {
            radius[0] = 0.99;
        }
        if (atmnum == 18) {
            radius[0] = 0.97;
        }
        if (atmnum == 35) {
            radius[0] = 1.14;
        }
        if (atmnum == 36) {
            radius[0] = 1.10;
        }
        if (atmnum == 53) {
            radius[0] = 1.33;
        }
        if (atmnum == 54) {
            radius[0] = 1.30;
        }
        radius[0] = 1.1 * radius[0];

        // Set polarity
        atmmass = 1;
        polarity[0] = 0.0;
        if (atmnum == 1) {
            atmmass = 1.008;
            polarity[0] = 0.496;
        } else if (atmnum == 5) {
            atmmass = 10.810;
            polarity[0] = 1.600;
        } else if (atmnum == 6) {
            atmmass = 12.011;
            polarity[0] = 1.334;
        } else if (atmnum == 7) {
            atmmass = 14.007;
            polarity[0] = 1.073;
        } else if (atmnum == 8) {
            atmmass = 15.999;
            polarity[0] = 0.837;
        } else if (atmnum == 9) {
            atmmass = 18.998;
        } else if (atmnum == 14) {
            atmmass = 28.086;
        } else if (atmnum == 15) {
            atmmass = 30.974;
            polarity[0] = 1.828;
        } else if (atmnum == 16) {
            atmmass = 32.066;
            polarity[0] = 2.800;
        } else if (atmnum == 17) {
            atmmass = 35.453;
            polarity[0] = 4.000;
        } else if (atmnum == 35) {
            atmmass = 79.904;
            polarity[0] = 5.650;
        } else if (atmnum == 53) {
            atmmass = 126.904;
            polarity[0] = 7.250;
        }
        double sixth = 1.0 / 6.0;
        pdamp[0] = Math.pow(polarity[0], sixth);

        AtomType at = new AtomType(type, type, element, null, atmnum, atmmass, 0);
        return at;
    }

    /**
     * Sets up connectivities based on radii and then prints this information
     * out to an xyz file
     *
     * @param name a {@link java.lang.String} object.
     */
    public void setup_print_xyz(String name) {
        //Set up connectivities based on radii
        for (int i = 0; i < atoms.length - 1; i++) {
            Atom ai = atoms[i];
            for (int j = i + 1; j < atoms.length; j++) {
                Atom aj = atoms[j];
                double xr = aj.getX() - ai.getX();
                double yr = aj.getY() - ai.getY();
                double zr = aj.getZ() - ai.getZ();
                double rij = ai.getBornRadius() + aj.getBornRadius();
                double dij = Math.sqrt(xr * xr + yr * yr + zr * zr);
                if (dij < rij) {;
                    Bond b = new Bond(ai, aj);
                }
            }
        }

        File outf = new File(name + ".xyz");
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outf)));
            DecimalFormat myFormatter = new DecimalFormat(" ##########0.00000;-##########0.00000");

            bw.write(String.format("%5d\n", atoms.length));
            for (Atom a : atoms) {
                String output = String.format("%6d", a.getAtomType().type) + "  " + a.getAtomType().name + " " + String.format("%12s %12s %12s", myFormatter.format(a.getX()), myFormatter.format(a.getY()), myFormatter.format(a.getZ())) + " " + String.format("%6d", a.getAtomType().atomClass);
                for (int i = 0; i < a.getBonds().size(); i++) {
                    output += String.format("%6d", a.getBonds().get(i).get1_2(a).getAtomType().type);
                }
                bw.write(output + "\n");
            }
            bw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * An example MultipoleFrameTypes would be: {401,401,404}*
     */
    public void setframeAutomatic() {
        for (int i = 0; i < atoms.length; i++) {
            double charge = globalMultipole[0][i][0];
            double dipole[] = {globalMultipole[0][i][1], globalMultipole[0][i][2], globalMultipole[0][i][3]};
            double quadrupole[][] = {{globalMultipole[0][i][4], globalMultipole[0][i][5], globalMultipole[0][i][6]}, {globalMultipole[0][i][7], globalMultipole[0][i][8], globalMultipole[0][i][9]}, {globalMultipole[0][i][10], globalMultipole[0][i][11], globalMultipole[0][i][12]}};
            int[] multipoleFrameTypes = new int[3];
            MultipoleFrameDefinition frameDefinition = null;
            Atom a = atoms[i];
            int j = a.getBonds().size();
            if (j == 0) {
                multipoleFrameTypes[0] = 0;
                multipoleFrameTypes[1] = 0;
                multipoleFrameTypes[2] = 0;
            } else if (j == 1) {
                Atom ia = a.getBonds().get(0).get1_2(a);
                multipoleFrameTypes[2] = ia.getType();
                if (ia.getBonds().size() == 1) {
                    multipoleFrameTypes[0] = 0;
                } else {
                    int m = 0;
                    for (int k = 0; k < ia.getBonds().size(); k++) {
                        Atom kb = ia.getBonds().get(k).get1_2(ia);
                        if (kb.getAtomicNumber() > m && kb.getIndex() != a.getIndex()) {
                            multipoleFrameTypes[0] = kb.getType();
                            m = kb.getAtomicNumber();
                        }
                    }
                }
                multipoleFrameTypes[1] = 0;
            } else if (j == 2) {
                Atom ia = a.getBonds().get(0).get1_2(a);
                Atom ib = a.getBonds().get(1).get1_2(a);
                Atom kab = priority(a, ia, ib);
                if (kab.getIndex() == ia.getIndex()) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                } else if (kab.getIndex() == ib.getIndex()) {
                    multipoleFrameTypes[2] = ib.getType();
                    multipoleFrameTypes[0] = ia.getType();
                } else {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                }
                multipoleFrameTypes[1] = 0;
            } else if (j == 3) {
                Atom ia = a.getBonds().get(0).get1_2(a);
                Atom ib = a.getBonds().get(1).get1_2(a);
                Atom ic = a.getBonds().get(2).get1_2(a);
                Atom kab = priority(a, ia, ib);
                Atom kac = priority(a, ia, ic);
                Atom kbc = priority(a, ib, ic);
                if (kab == null && kac == null) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                } else if (kab == null) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kac == null) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kbc == null) {
                    multipoleFrameTypes[2] = ib.getType();
                    multipoleFrameTypes[0] = ic.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kab.getIndex() == ia.getIndex() && kac.getIndex() == ia.getIndex()) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ((kbc.getIndex() == ic.getIndex()) ? ic.getType() : ib.getType());
                } else if (kab.getIndex() == ib.getIndex() && kbc.getIndex() == ib.getIndex()) {
                    multipoleFrameTypes[2] = ib.getType();
                    multipoleFrameTypes[0] = ((kac.getIndex() == ic.getIndex()) ? ic.getType() : ia.getType());
                } else if (kac.getIndex() == ic.getIndex() && kbc.getIndex() == ic.getIndex()) {
                    multipoleFrameTypes[2] = ic.getType();
                    multipoleFrameTypes[0] = ((kab.getIndex() == ib.getIndex()) ? ib.getType() : ia.getType());
                }

            } else if (j == 4) {
                Atom ia = a.getBonds().get(0).get1_2(a);
                Atom ib = a.getBonds().get(1).get1_2(a);
                Atom ic = a.getBonds().get(2).get1_2(a);
                Atom id = a.getBonds().get(3).get1_2(a);
                Atom kab = priority(a, ia, ib);
                Atom kac = priority(a, ia, ic);
                Atom kad = priority(a, ia, id);
                Atom kbc = priority(a, ib, ic);
                Atom kbd = priority(a, ib, id);
                Atom kcd = priority(a, ic, id);
                if (kab == null && kac == null && kad == null) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                } else if (kab.getIndex() == ia.getIndex() && kac.getIndex() == ia.getIndex() && kad.getIndex() == ia.getIndex()) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                    if (kbc.getIndex() == ic.getIndex() && (kcd == null || kcd.getIndex() == ic.getIndex())) {
                        multipoleFrameTypes[0] = ic.getType();
                    } else if (kbd.getIndex() == id.getIndex() && kcd.getIndex() == id.getIndex()) {
                        multipoleFrameTypes[0] = id.getType();
                    }
                } else if (kab.getIndex() == ib.getIndex() && kbc.getIndex() == ib.getIndex() && kbd.getIndex() == ib.getIndex()) {
                    multipoleFrameTypes[2] = ib.getType();
                    multipoleFrameTypes[0] = ia.getType();
                    if (kac.getIndex() == ic.getIndex() && (kcd == null || kcd.getIndex() == ic.getIndex())) {
                        multipoleFrameTypes[0] = ic.getType();
                    } else if (kad.getIndex() == id.getIndex() && kcd.getIndex() == id.getIndex()) {
                        multipoleFrameTypes[0] = id.getType();
                    }
                } else if (kac.getIndex() == ic.getIndex() && kbc.getIndex() == ic.getIndex() && kcd.getIndex() == ic.getIndex()) {
                    multipoleFrameTypes[2] = ic.getType();
                    multipoleFrameTypes[0] = ia.getType();
                    if (kab.getIndex() == ib.getIndex() && (kbd == null || kbd.getIndex() == ib.getIndex())) {
                        multipoleFrameTypes[0] = ib.getType();
                    }
                    if (kad.getIndex() == id.getIndex() && kbd.getIndex() == id.getIndex()) {
                        multipoleFrameTypes[0] = id.getType();
                    }
                } else if (kad.getIndex() == id.getIndex() && kbd.getIndex() == id.getIndex()) {
                    multipoleFrameTypes[2] = id.getType();
                    multipoleFrameTypes[0] = ia.getType();
                    if (kab.getIndex() == ib.getIndex() && (kbc == null || kbc.getIndex() == ib.getIndex())) {
                        multipoleFrameTypes[0] = ib.getType();
                    } else if (kac.getIndex() == ic.getIndex() && kbc.getIndex() == ic.getIndex()) {
                        multipoleFrameTypes[0] = ic.getType();
                    }
                } else if (kab == null && kac.getIndex() == ia.getIndex() && kad.getIndex() == ia.getIndex()) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ib.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kac == null && kab.getIndex() == ia.getIndex() && kad.getIndex() == ia.getIndex()) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = ic.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kad == null && kab.getIndex() == ia.getIndex() && kac.getIndex() == ia.getIndex()) {
                    multipoleFrameTypes[2] = ia.getType();
                    multipoleFrameTypes[0] = id.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kbc == null && kab.getIndex() == ib.getIndex() && kbd.getIndex() == ib.getIndex()) {
                    multipoleFrameTypes[2] = ib.getType();
                    multipoleFrameTypes[0] = ic.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kbd == null && kab.getIndex() == ib.getIndex() && kbc.getIndex() == ib.getIndex()) {
                    multipoleFrameTypes[2] = ib.getType();
                    multipoleFrameTypes[0] = id.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                } else if (kcd == null && kac.getIndex() == ic.getIndex() && kbc.getIndex() == ic.getIndex()) {
                    multipoleFrameTypes[2] = ic.getType();
                    multipoleFrameTypes[0] = id.getType();
                    frameDefinition = MultipoleFrameDefinition.BISECTOR;
                }
            }
            //Set up stuff. Charge, Dipole, Quadrupole, FrameTypes, FrameDef
            MultipoleType m = new MultipoleType(charge, dipole.clone(), quadrupole.clone(), multipoleFrameTypes.clone(), frameDefinition);
            a.setMultipoleType(m);
        }
    }

    /**
     * <p>
     * setframe</p>
     *
     * @param peditinfname a {@link java.lang.String} object.
     */
    public void setframe(String peditinfname) {
        axisAtom = new int[atoms.length][3];
        xaxis = new int[atoms.length];
        frame = new MultipoleFrameDefinition[atoms.length];
        File peditin = new File(peditinfname);
        try {
            if (peditin != null && peditin.exists() && peditin.canRead()) {
                BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(peditin)));
                String line;
                int ia, ib, ic;
                while (!(line = br.readLine()).equals("")) {
                    int index = Integer.parseInt(line.split(" +")[0]) - 1;
                    axisAtom[index][2] = -1;
                    ia = Integer.parseInt(line.split(" +")[1]);
                    axisAtom[index][0] = Math.abs(ia) - 1;
                    ib = Integer.parseInt(line.split(" +")[2]);
                    axisAtom[index][1] = Math.abs(ib) - 1;
                    xaxis[index] = axisAtom[index][1];
                    ic = -1;
                    if (line.split(" +").length > 3) {
                        ic = Integer.parseInt(line.split(" +")[3]);
                        axisAtom[index][2] = Math.abs(ic) - 1;
                    }
//                    if(ia > 0 && ib == 0){
//                    	for(int r = 0; r < atoms[index].getNumBonds(); r++){
//                        	Atom b = atoms[index].getBonds().get(r).get1_2(atoms[index]);
//                        	if(b.getIndex() != Math.abs(ia)){
//                        		ib = b.getIndex();
//                        		axisAtom[index][1] = ib - 1;
//                        		break;
//                        	}
//                    	}
//                    	axisAtom[index][1] = ib;
//                    }
                    if (ia > 0 && ib >= 0) {
                        frame[index] = MultipoleFrameDefinition.ZTHENX;
                    }
                    if (ia < 0 || ib < 0) {
                        frame[index] = MultipoleFrameDefinition.BISECTOR;
                        if (ib < 0) {
                            xaxis[index] = -1;
                        }
                    }
                    if (line.split(" +").length > 3 && ib < 0 && ic < 0) {
                        frame[index] = MultipoleFrameDefinition.ZTHENBISECTOR;
                    }
                    if (line.split(" +").length > 3 && ia < 0 && ib < 0 && ic < 0) {
                        frame[index] = MultipoleFrameDefinition.TRISECTOR;
                    }
                }

                for (int i = 0; i < nAtoms; i++) {
                    Atom a = atoms[i];
                    int[] polgrp = new int[a.getNumBonds()];
                    for (int j = 0; j < a.getNumBonds(); j++) {
                        polgrp[j] = a.getBonds().get(j).get1_2(a).getIndex();
                    }
                    PolarizeType p = new PolarizeType(a.getType(), polarizability[i], thole[i], polgrp.clone());
                    a.setPolarizeType(p);
                }

                while ((line = br.readLine()) != null) {
                    if (line.contains("Y")) {
                        remove_symmetry = true;
                        break;
                    } else if (line.equals("")) {
                        continue;
                    }
                    if (isInteger(line.split(" +")[0]) && isInteger(line.split(" +")[1])) {
                        ia = Integer.parseInt(line.split(" +")[0]) - 1;
                        ib = Integer.parseInt(line.split(" +")[1]) - 1;
                        for (int i = 0; i < atoms[ia].getNumBonds(); i++) {
                            if (atoms[ia].getPolarizeType().polarizationGroup[i] == ib + 1) {
                                for (int j = i + 1; j < atoms[ia].getNumBonds(); j++) {
                                    atoms[ia].getPolarizeType().polarizationGroup[j - 1] = atoms[ia].getPolarizeType().polarizationGroup[j];
                                }
                                atoms[ia].getPolarizeType().polarizationGroup[atoms[ia].getNumBonds() - 1] = 0;
                            }
                        }
                        for (int i = 0; i < atoms[ib].getNumBonds(); i++) {
                            if (atoms[ib].getPolarizeType().polarizationGroup[i] == ia + 1) {
                                for (int j = i + 1; j < atoms[ib].getNumBonds(); j++) {
                                    atoms[ib].getPolarizeType().polarizationGroup[j - 1] = atoms[ib].getPolarizeType().polarizationGroup[j];
                                }
                                atoms[ib].getPolarizeType().polarizationGroup[atoms[ib].getNumBonds() - 1] = 0;
                            }
                        }
                    } else if (isInteger(line.split(" +")[0]) && !isInteger(line.split(" +")[1])) {
                        ia = Integer.parseInt(line.split(" +")[0]) - 1;
                        polarizability[ia] = Double.parseDouble(line.split(" +")[1]);

                    }

                }
            }
        } catch (Exception e) {
            Logger.getLogger(Poledit.class.getName()).log(Level.SEVERE, null, e);
        }

    }

    /**
     * <p>
     * isInteger</p>
     *
     * @param string a {@link java.lang.String} object.
     * @return a boolean.
     */
    public boolean isInteger(String string) {
        try {
            Integer.valueOf(string);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    /**
     * <p>
     * priority</p>
     *
     * @param o a {@link ffx.potential.bonded.Atom} object.
     * @param a a {@link ffx.potential.bonded.Atom} object.
     * @param b a {@link ffx.potential.bonded.Atom} object.
     * @return a {@link ffx.potential.bonded.Atom} object.
     */
    public Atom priority(Atom o, Atom a, Atom b) {
        int ka = a.getAtomicNumber();
        int kb = b.getAtomicNumber();
        if (ka > kb) {
            return a;
        } else if (kb > ka) {
            return b;
        } else {
            ka = 0;
            for (int i = 0; i < a.getBonds().size(); i++) {
                Atom m = a.getBonds().get(i).get1_2(a);
                if (o.getIndex() != m.getIndex()) {
                    if (m.getAtomicNumber() > ka) {
                        ka = m.getAtomicNumber();
                    }
                }
            }
            kb = 0;
            for (int i = 0; i < b.getBonds().size(); i++) {
                Atom m = b.getBonds().get(i).get1_2(b);
                if (o.getIndex() != m.getIndex()) {
                    if (m.getAtomicNumber() > kb) {
                        kb = m.getAtomicNumber();
                    }
                }
            }
            if (a.getBonds().size() > b.getBonds().size()) {
                return a;
            } else if (b.getBonds().size() > a.getBonds().size()) {
                return b;
            } else if (ka > kb) {
                return a;
            } else if (kb > kb) {
                return b;
            } else {
                return null;
            }
        }
    }

    /**
     * <p>
     * printlocalmpoles</p>
     */
    public void printlocalmpoles() {
        for (int i = 0; i < nAtoms; i++) {
            System.out.println("\n\nSite: " + (i + 1) + " Name: " + atoms[i].getAtomType().name + " Atomic Number: " + atoms[i].getAtomicNumber());
            System.out.println("\n\nLocal Frame: " + frame[i] + " " + (axisAtom[i][0] + 1) + " " + (xaxis[i] >= 0 ? axisAtom[i][1] + 1 : 0) + " " + (axisAtom[i][2] >= 0 ? (axisAtom[i][2] + 1) : 0));
            System.out.println("\nCharge: " + localMultipole[i][0]);
            System.out.println("\nDipole: " + localMultipole[i][1] / BOHR + " " + localMultipole[i][2] / BOHR + " " + localMultipole[i][3] / BOHR);
            System.out.println("\nQuadrupole: " + 3 * localMultipole[i][4] / (BOHR * BOHR));
            System.out.println("\n            " + 3 * localMultipole[i][7] / (BOHR * BOHR) + " " + 3 * localMultipole[i][8] / (BOHR * BOHR));
            System.out.println("\n            " + 3 * localMultipole[i][10] / (BOHR * BOHR) + " " + 3 * localMultipole[i][11] / (BOHR * BOHR) + " " + 3 * localMultipole[i][12] / (BOHR * BOHR));
        }
    }

    /**
     * <p>
     * printglobalmpoles</p>
     */
    public void printglobalmpoles() {
        for (int i = 0; i < nAtoms; i++) {
            System.out.println("\n\nSite: " + (i + 1) + " Name: " + atoms[i].getAtomType().name + " Atomic Number: " + atoms[i].getAtomicNumber());
            System.out.println("\n\nLocal Frame: " + frame[i] + " " + (axisAtom[i][0] + 1) + " " + (xaxis[i] >= 0 ? axisAtom[i][1] + 1 : 0) + " " + (axisAtom[i][2] >= 0 ? (axisAtom[i][2] + 1) : 0));
            System.out.println("\nCharge: " + globalMultipole[0][i][0]);
            System.out.println("\nDipole: " + globalMultipole[0][i][1] / BOHR + " " + globalMultipole[0][i][2] / BOHR + " " + globalMultipole[0][i][3] / BOHR);
            System.out.println("\nQuadrupole: " + 3 * globalMultipole[0][i][4] / (BOHR * BOHR));
            System.out.println("\n            " + 3 * globalMultipole[0][i][7] / (BOHR * BOHR) + " " + 3 * globalMultipole[0][i][8] / (BOHR * BOHR));
            System.out.println("\n            " + 3 * globalMultipole[0][i][10] / (BOHR * BOHR) + " " + 3 * globalMultipole[0][i][11] / (BOHR * BOHR) + " " + 3 * globalMultipole[0][i][12] / (BOHR * BOHR));
        }
    }

    /**
     * <p>
     * fixpolar</p>
     */
    public void fixpolar() {
        if (remove_symmetry) {
            for (int i = 0; i < nAtoms; i++) {
                boolean yzero = false;
                if (axisAtom[i][2] == -1) {
                    yzero = true;
                }
                if (frame[i] == MultipoleType.MultipoleFrameDefinition.BISECTOR) {
                    yzero = true;
                }
                if (frame[i] == MultipoleType.MultipoleFrameDefinition.ZTHENBISECTOR) {
                    yzero = true;
                }
                //check, what represents the zaxis atom? axisAtom[i][1] or axisAtom[i][0]
                if (axisAtom[i][0] == -1) {
                    localMultipole[i][12] = 0;
                }
                if (xaxis[i] == -1) {
                    localMultipole[i][1] = 0;
                    localMultipole[i][4] = -0.5 * localMultipole[i][12];
                    localMultipole[i][6] = 0;
                    localMultipole[i][8] = localMultipole[i][5];
                    localMultipole[i][10] = 0;
                }
                if (yzero) {
                    localMultipole[i][2] = 0;
                    localMultipole[i][5] = 0;
                    localMultipole[i][7] = 0;
                    localMultipole[i][9] = 0;
                    localMultipole[i][11] = 0;
                }

            }
        }
        //maintain integer net charge for whole system
        int k = 0;
        double big = 0.0;
        double sum = 0.0;
        for (int i = 0; i < nAtoms; i++) {
            sum = sum + localMultipole[i][0];
            double ci = Math.abs(localMultipole[i][0]);
            if (ci > big) {
                k = i;
                big = ci;
            }
        }
        sum = sum - Math.round(sum);
        if (k != 0) {
            localMultipole[k][0] = localMultipole[k][0] - sum;
        }

        //maintain traceless quadrupole at each multipole site
        for (int i = 0; i < nAtoms; i++) {
            sum = localMultipole[i][4] + localMultipole[i][8] + localMultipole[i][12];
            big = Math.max(localMultipole[i][4], Math.max(localMultipole[i][8], localMultipole[i][12]));
            k = 0;
            if (big == Math.abs(localMultipole[i][4])) {
                k = 4;
            } else if (big == Math.abs(localMultipole[i][8])) {
                k = 8;
            } else if (big == Math.abs(localMultipole[i][12])) {
                k = 12;
            }
            if (k != 0) {
                localMultipole[i][k] = localMultipole[i][k] - sum;
            }
        }

    }

    /**
     * <p>
     * removeInducedFromGlobal</p>
     */
    public void removeInducedFromGlobal() {
        for (int i = 0; i < nAtoms; i++) {
            for (int j = 0; j < 3; j++) {
                globalMultipole[0][i][j + 1] = globalMultipole[0][i][j + 1] - inducedDipole[0][i][j];
            }
        }
    }

    /**
     * <p>
     * prtpolar</p>
     *
     * @param name a {@link java.lang.String} object.
     */
    public void prtpolar(String name) {
        File outf = new File(name + ".key");
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outf)));
            DecimalFormat myFormatter = new DecimalFormat(" ##########0.00000;-##########0.00000");
            String output = " ";
            int index = name.lastIndexOf("/");
            String molname = name.substring(index + 1, name.length());
            for (int i = 0; i < atoms.length; i++) {
                output = String.format("atom%5d%5d%3s \"%-20s\" %5d%6s%5d", atoms[i].getIndex(), atoms[i].getType(), atoms[i].getAtomType().name, molname, atoms[i].getAtomicNumber(), myFormatter.format(atoms[i].getMass()), atoms[i].getNumBonds());
                //output = "atom\t"+atoms[i].getIndex()+" "+atoms[i].getType()+" "+atoms[i].getName()+" "+molname+" "+atoms[i].getAtomicNumber()+" "+atoms[i].getMass()+" "+atoms[i].getNumBonds();
                bw.write(output + "\n");
            }
            bw.write("\n");
            for (int i = 0; i < atoms.length; i++) {
                if (frame[i] == MultipoleType.MultipoleFrameDefinition.ZTHENX) {
                    if (axisAtom[i][2] == -1) {
                        //output = String.format("multipole %5d %5d %5d%9s%s",atoms[i].getType(),axisAtom[i][0]+1,(xaxis[i] > 0 ? axisAtom[i][1]+1 : 0)," ",myFormatter.format(localMultipole[i][0]));
                        output = String.format("multipole %5d %5d %5d%9s%s", atoms[i].getType(), axisAtom[i][0] + 1, axisAtom[i][1] + 1, " ", myFormatter.format(localMultipole[i][0]));

                    } else {
                        //output = String.format("multipole %5d %5d %5d %5d%9s%s",atoms[i].getType(),axisAtom[i][0]+1,(xaxis[i] > 0 ? axisAtom[i][1]+1 : 0),axisAtom[i][2]+1,myFormatter.format(localMultipole[i][0]));
                        output = String.format("multipole %5d %5d %5d %5d%9s%s", atoms[i].getType(), axisAtom[i][0] + 1, axisAtom[i][1] + 1, axisAtom[i][2] + 1, myFormatter.format(localMultipole[i][0]));

                    }
                } else if (frame[i] == MultipoleType.MultipoleFrameDefinition.BISECTOR) {
                    if (axisAtom[i][2] == -1) {

                        //output = String.format("multipole %5d %5d %5d%9s%s",atoms[i].getType(),axisAtom[i][0]+1,(xaxis[i] > 0 ? -(axisAtom[i][1]+1) : 0)," ",myFormatter.format(localMultipole[i][0]));
                        output = String.format("multipole %5d %5d %5d%9s%s", atoms[i].getType(), axisAtom[i][0] + 1, -(axisAtom[i][1] + 1), " ", myFormatter.format(localMultipole[i][0]));

                    } else {
                        //output = String.format("multipole %5d %5d %5d %5d%9s%s",atoms[i].getType(),axisAtom[i][0]+1,(xaxis[i] > 0 ? -(axisAtom[i][1]+1) : 0),axisAtom[i][2]+1,myFormatter.format(localMultipole[i][0]));
                        output = String.format("multipole %5d %5d %5d %5d%9s%s", atoms[i].getType(), axisAtom[i][0] + 1, -(axisAtom[i][1] + 1), axisAtom[i][2] + 1, myFormatter.format(localMultipole[i][0]));
                    }
                } else if (frame[i] == MultipoleType.MultipoleFrameDefinition.ZTHENBISECTOR) {
                    //output = String.format("multipole %5d %5d %5d %5d%9s%s",atoms[i].getType(),axisAtom[i][0]+1,(xaxis[i] > 0 ? -(axisAtom[i][1]+1) : 0),-(axisAtom[i][2]+1),myFormatter.format(localMultipole[i][0]));
                    output = String.format("multipole %5d %5d %5d %5d%9s%s", atoms[i].getType(), axisAtom[i][0] + 1, -(axisAtom[i][1] + 1), -(axisAtom[i][2] + 1), myFormatter.format(localMultipole[i][0]));
                } else if (frame[i] == MultipoleType.MultipoleFrameDefinition.TRISECTOR) {
                    //output = String.format("multipole %5d %5d %5d %5d%9s%11.5f",atoms[i].getType(),-axisAtom[i][0]+1,(xaxis[i] > 0 ? -(axisAtom[i][1]+1) : 0),-(axisAtom[i][2]+1),localMultipole[i][0]);
                    output = String.format("multipole %5d %5d %5d %5d%9s%11.5f", atoms[i].getType(), -axisAtom[i][0] + 1, -(axisAtom[i][1] + 1), -(axisAtom[i][2] + 1), localMultipole[i][0]);

                }
                bw.write(output + "\n");
                output = String.format("%35s", "") + " " + String.format("%s %s %s", myFormatter.format(localMultipole[i][1] / BOHR), myFormatter.format(localMultipole[i][2] / BOHR), myFormatter.format(localMultipole[i][3] / BOHR));
                bw.write(output + "\n");
                output = String.format("%35s", "") + " " + String.format("%s", myFormatter.format(3 * localMultipole[i][4] / (BOHR * BOHR)));
                bw.write(output + "\n");
                output = String.format("%35s", "") + " " + String.format("%s %s", myFormatter.format(3 * localMultipole[i][7] / (BOHR * BOHR)), myFormatter.format(3 * localMultipole[i][8] / (BOHR * BOHR)));
                bw.write(output + "\n");
                output = String.format("%35s", "") + " " + String.format("%s %s %s", myFormatter.format(3 * localMultipole[i][10] / (BOHR * BOHR)), myFormatter.format(3 * localMultipole[i][11] / (BOHR * BOHR)), myFormatter.format(3 * localMultipole[i][12] / (BOHR * BOHR)));
                bw.write(output + "\n");

            }
            bw.write("\n");
            for (int i = 0; i < nAtoms; i++) {
                output = String.format("polarize %5d %29s %7s", atoms[i].getIndex(), myFormatter.format(polarizability[i]), myFormatter.format(thole[i]));
                for (int k = 0; k < atoms[i].getPolarizeType().polarizationGroup.length; k++) {
                    if (atoms[i].getPolarizeType().polarizationGroup[k] != 0) {
                        output = output + "   " + atoms[i].getPolarizeType().polarizationGroup[k];
                    }
                }
                bw.write(output + "\n");
            }
            bw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * <p>
     * main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String args[]) {
        //Poledit p = new Poledit("/users/gchattree/Research/Compounds/test_compounds/12-ethanediol-test/12-ethanediol.gdmaout", "/users/gchattree/Research/Compounds/test_compounds/12-ethanediol-test/12-ethanediol-peditin.txt");
        //Poledit p2 = new Poledit("/users/gchattree/Research/Compounds/test_compounds/phenobarbital-tinker-goal/phenobarbital.gdmaout","/users/gchattree/Research/Compounds/test_compounds/phenobarbital-test/phenobarbital-peditin.txt");
        //Poledit p3 = new Poledit("/users/gchattree/Research/Compounds/poltypeffx-2/di-n-propyl_sulfide-test/di-n-propyl_sulfide.gdmaout", "/users/gchattree/Research/Compounds/poltypeffx-2/di-n-propyl_sulfide-test/di-n-propyl_sulfide-peditin.txt");
        //Poledit p4 = new Poledit("/users/gchattree/Research/Compounds/easycompounds/2-ethoxyethanol/2-ethoxyethanol.gdmaout","/users/gchattree/Research/Compounds/easycompounds/2-ethoxyethanol/2-ethoxyethanol-peditin.txt");
    }
}
