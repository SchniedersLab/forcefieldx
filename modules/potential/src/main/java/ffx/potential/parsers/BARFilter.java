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
    private int count = 0;
    private double temp;


    //pass in four 1 dimensional arrays
    public BARFilter(File barFile) {
        this.givenFile = barFile;
    }


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

    public void closeReader() {

    }

    public boolean readFile() {
        ArrayList<Double> ens1lam1 = new ArrayList<Double>();
        ArrayList<Double> ens1lam2 = new ArrayList<Double>();
        ArrayList<Double> ens2lam1 = new ArrayList<Double>();
        ArrayList<Double> ens2lam2 = new ArrayList<Double>();
        ArrayList<Double> vol1 = new ArrayList<Double>();
        ArrayList<Double> vol2 = new ArrayList<Double>();
        try (BufferedReader br = new BufferedReader(new FileReader(givenFile))) {

            String str = "";
            String data = "";

            while ((data = br.readLine()) != null) {
                count += 1;

                String[] tokens = data.trim().split(" +");
                if (data.contains(".xyz") || tokens.length < 3) {
                    snaps = Integer.parseInt(tokens[0]);
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
            e1l1 = new double[snaps];
            e1l2 = new double[snaps];
            e2l1 = new double[snaps];
            e2l2 = new double[snaps];
            volume1 = new double[snaps];
            volume2 = new double[snaps];

            for (int i = 0; i < ens1lam1.size(); i++) {
                e1l1[i] = ens1lam1.get(i);
                e1l2[i] = ens1lam2.get(i);
                e2l1[i] = ens2lam1.get(i);
                e2l2[i] = ens2lam2.get(i);


                if (!vol1.isEmpty()) {
                    volume1[i] = vol1.get(i);
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


    //one file at a time
    public boolean writeFile(String saveFile, boolean isPBC) throws IOException {
        int snaps = e1l1.length;
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

            bw.write(format("%8d %9.3f  %s\n", snaps, temp, name));
            for (int i = 0; i < snaps; i++) {
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
