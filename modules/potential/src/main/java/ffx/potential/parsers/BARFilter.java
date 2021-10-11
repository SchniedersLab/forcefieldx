//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
//******************************************************************************
package ffx.potential.parsers;

import ffx.potential.Utilities;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BondType;

import static ffx.potential.bonded.Bond.logNoBondType;
import static java.lang.Double.parseDouble;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import java.io.*;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;
/**
 * The BARFilter class parses TINKER bar(*.BAR) files.
 *
 * @author Rose A. Gogal
 * @since 1.0
 */

public class BARFilter {

    private BufferedReader bufferedReader = null;
    private String remarkLine;
    private static final Logger logger = Logger.getLogger(XYZFilter.class.getName());
    private File givenFile;
    private double[] e1l1;
    private double[] e1l2;
    private double[] e2l1;
    private double[] e2l2;
    private double[] volume1;
    private double[] volume2;
    private int snaps;
    private int snaps2;
    private int startingSnap = 0;
    private int endingSnap = 0;
    private int count = 0;
    private double temp;
    private boolean getLambda = false;



    /**
     * BARFilter constructor
     * @param barFile a {@link java.util.List} object.
     *
     */
    public BARFilter(File barFile) {
        this.givenFile = barFile;

    }

    /**
     * BARFilter constructor
     * @param barFile a {@link java.util.List} object.
     * @param startingSnap a {@link java.util.List} object.
     * @param endingSnap a {@link java.util.List} object.
     */
    public BARFilter(File barFile, int startingSnap, int endingSnap) {

        this.givenFile = barFile;
        this.startingSnap = startingSnap;
        this.endingSnap = endingSnap;

    }

    /**
     * BARFilter constructor
     * @param xyzFile a {@link java.util.List} object.
     * @param e1l1 energy in ensemble 1 at lambda 1
     * @param e1l2 energy in ensemble 1 at lambda 2
     * @param e2l1 energy in ensemble 2 at lambda 1
     * @param e2l2 energy in ensemble 2 at lambda 2
     * @param volume1 volume in ensemble 1
     * @param volume2 volume in ensemble 2
     * @param temp temperature
     */
    public BARFilter(File xyzFile,
                     double[] e1l1, double[] e1l2,
                     double[] e2l1,
                     double[] e2l2,
                     double[] volume1,
                     double[] volume2, double temp) {
        this.givenFile = xyzFile;
        this.e1l1 = e1l1;
        this.e1l2 = e1l2;
        this.e2l1 = e2l1;
        this.e2l2 = e2l2;
        this.volume1 = volume1;
        this.volume2 = volume2;
        this.temp = temp;

    }



    /**
     * Read TINKER bar files and parse the snapshots into energy arrays
     * @return
     */
    public boolean readFile() {
        ArrayList<Double> ens1lam1 = new ArrayList<Double>();
        ArrayList<Double> ens1lam2 = new ArrayList<Double>();
        ArrayList<Double> ens2lam1 = new ArrayList<Double>();
        ArrayList<Double> ens2lam2 = new ArrayList<Double>();
        ArrayList<Double> vol1 = new ArrayList<Double>();
        ArrayList<Double> vol2 = new ArrayList<Double>();
        int snapshots = 0;
        int xyzCount=0;
        try (BufferedReader br = new BufferedReader(new FileReader(givenFile))) {

            String str = "";
            String data = "";

            while ((data = br.readLine()) != null) {
                count += 1;
                String[] tokens = data.trim().split(" +");
                if(startingSnap != 0 || endingSnap !=0){
                    if(startingSnap != 0 && endingSnap == 0){
                        endingSnap = snaps;
                    }
                    if (data.contains(".xyz") || tokens.length < 3) {
                        snaps = Integer.parseInt(tokens[0]);
                    }
                    snapshots = (endingSnap - startingSnap) + 1;
                    if(count >= startingSnap + 1 && count <= endingSnap + 1){
                        if (tokens.length == 4) {
                            vol1.add(Double.parseDouble(tokens[3]));
                        }
                        ens1lam1.add(Double.parseDouble(tokens[1]));
                        ens1lam2.add(Double.parseDouble(tokens[2]));
                    } else if (count >= snaps+startingSnap+2 && count <= snaps + endingSnap + 2 ){
                        if (tokens.length == 4) {
                            vol2.add(Double.parseDouble(tokens[3]));
                        }
                        ens2lam1.add(Double.parseDouble(tokens[1]));
                        ens2lam2.add(Double.parseDouble(tokens[2]));
                    }

                } else {

                    if (data.contains(".xyz") || tokens.length < 3) {
                        xyzCount += 1;
                        if(xyzCount == 1) {
                            snaps = Integer.parseInt(tokens[0]);
                        } else if (xyzCount ==2){
                            snaps2 = Integer.parseInt(tokens[0]);
                        }

                    } else if (count <= snaps + 1 && count != 1) {
                        if (tokens.length == 4) {
                            vol1.add(Double.parseDouble(tokens[3]));
                        }

                        ens1lam1.add(Double.parseDouble(tokens[1]));
                        ens1lam2.add(Double.parseDouble(tokens[2]));
                    } else if (count > snaps + 2) {
                        if (tokens.length == 4) {
                            vol2.add(Double.parseDouble(tokens[3]));
                        }
                        ens2lam1.add(Double.parseDouble(tokens[1]));
                        ens2lam2.add(Double.parseDouble(tokens[2]));
                    }
                }

            }
            if(snapshots != 0){
                e1l1 = new double[snapshots];
                e1l2 = new double[snapshots];
                e2l1 = new double[snapshots];
                e2l2 = new double[snapshots];
                volume1 = new double[snapshots];
                volume2 = new double[snapshots];
                snaps = snapshots;
            } else {
                e1l1 = new double[snaps];
                e1l2 = new double[snaps];
                e2l1 = new double[snaps2];
                e2l2 = new double[snaps2];
                volume1 = new double[snaps];
                volume2 = new double[snaps2];
            }


            for (int i = 0; i < ens1lam1.size(); i++) {
                e1l1[i] = ens1lam1.get(i);
                e1l2[i] = ens1lam2.get(i);
                if (!vol1.isEmpty()) {
                    volume1[i] = vol1.get(i);
                }

            }

            for (int i = 0; i < ens2lam1.size(); i++) {
                e2l1[i] = ens2lam1.get(i);
                e2l2[i] = ens2lam2.get(i);
                if (!vol1.isEmpty()) {
                    volume2[i] = vol2.get(i);
                }

            }



            // Read blank lines at the top of the file

            if (data == null) {
                return false;
            }


        } catch (FileNotFoundException fileNotFoundException) {
            fileNotFoundException.printStackTrace();
        } catch (IOException ioException) {
            ioException.printStackTrace();
        }

        return true;


    }


    /**
     * Write TINKER bar files
     * @param saveFile filename
     * @param isPBC include volume
     * @return
     * @throws IOException
     */
    public boolean writeFile(String saveFile, boolean isPBC) throws IOException {
        int snaps = e1l1.length;
        int snaps2 = e2l1.length;
        String name = givenFile.getName();


        File newFile = new File(saveFile);
        logger.info(format("\n Writing Tinker-compatible BAR file to %s.", newFile));
        try (FileWriter fw = new FileWriter(newFile, newFile.exists());
             BufferedWriter bw = new BufferedWriter(fw)) {
            StringBuilder fileName = new StringBuilder();
            fileName.append("  ").append(name);
            bw.write(format("%8d %9.3f %s\n", snaps, temp, name));
            for (int i = 0; i < snaps; i++) {
                if (isPBC) {
                    bw.write(format("%8d %20.10f %20.10f %20.10f\n", i + 1, e1l1[i], e1l2[i], volume2[i]));
                } else {
                    bw.write(format("%8d %20.10f %20.10f\n", i + 1, e1l1[i], e1l2[i]));
                }
            }
            StringBuilder fileName2 = new StringBuilder();

            fileName2.append("  ").append(name);

            bw.write(format("%8d %9.3f  %s\n", snaps2, temp, name));
            for (int i = 0; i < snaps2; i++) {
                if (isPBC) {
                    bw.write(format("%8d %20.10f %20.10f %20.10f\n", i + 1, e2l1[i], e2l2[i], volume2[i]));
                } else {
                    bw.write(format("%8d %20.10f %20.10f\n", i + 1, e2l1[i], e2l2[i]));
                }
            }
        }


        return false;
    }

    public double[] getE1l1() {
        return e1l1;
    }

    public double[] getE2l1() {
        return e2l1;
    }

    public double[] getE2l2() {
        return e2l2;
    }

    public double[] getE1l2() {
        return e1l2;
    }

    public double[] getVolume1() {
        return volume1;
    }

    public double[] getVolume2() {
        return volume2;
    }

    public int getSnaps() {
        return snaps;
    }


}
