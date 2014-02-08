/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
 */
package ffx.autoparm;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.logging.Logger;

import static java.lang.Math.max;
import static java.lang.String.format;

import org.apache.commons.configuration.CompositeConfiguration;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaals.VDW_FORM;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parsers.XYZFilter;

import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.r;

/**
 * Compute the potential energy and derivatives of an AMOEBA system.
 *
 * @author Gaurav Chattree and Michael J. Schnieders
 * @since 1.0
 *
 */
public class Energy {

    private static final Logger logger = Logger.getLogger(Energy.class.getName());
    private Atom[] atoms;
    private Crystal crystal;
    private ParallelTeam parallelTeam;
    private VanDerWaals vanderWaals;
    private PME_2 pme2;
    protected int nAtoms;
    protected int nVanDerWaals, nPME;
    private File structure_key;
    private File structure_xyz;
    InputStreamReader stdinput = new InputStreamReader(System.in);
    BufferedReader stdreader = new BufferedReader(stdinput);
    private MolecularAssembly molecularAssembly;
    private ArrayList<String> key = new ArrayList<String>();
    private ForceField forceField;
    private boolean do_propyze = false;
    private boolean do_detail = false;

    /**
     * <p>Constructor for Energy.</p>
     *
     * @param xyz_filename a {@link java.lang.String} object.
     * @param keyfname a {@link java.lang.String} object.
     * @param options a {@link java.lang.String} object.
     * @throws java.io.IOException if any.
     */
    public Energy(String xyz_filename, String keyfname, String options) throws IOException {

        parallelTeam = new ParallelTeam();
        logger.info(format(" Constructing Force Field"));
        logger.info(format("\n SMP threads:                        %10d", parallelTeam.getThreadCount()));

        structure_xyz = new File(xyz_filename);
        if (!(structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead())) {
            System.out.println("Couldn't find xyz file");
            System.exit(1);
        }
        int n = 1;
        String oxyzfname = null;
        String old = xyz_filename;
        while (structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead()) {
            oxyzfname = xyz_filename;
            n++;
            xyz_filename = old;
            xyz_filename = xyz_filename + "_" + Integer.toString(n);
            structure_xyz = new File(xyz_filename);
        }
        structure_xyz = new File(oxyzfname);
        int index = oxyzfname.lastIndexOf(".");
        String name = oxyzfname.substring(0, index);

        if (keyfname != null) {
            structure_key = new File(keyfname);
            if (!(structure_key != null && structure_key.exists() && structure_key.canRead())) {
                System.out.println("Couldn't find key file");
                System.exit(1);
            }
        } else {
            keyfname = name + ".key";
            structure_key = new File(keyfname);
            if (!(structure_key != null && structure_key.exists() && structure_key.canRead())) {
                System.out.println("Couldn't find key file");
                System.exit(1);
            }
            n = 1;
            String okeyfname = null;
            old = keyfname;
            while (structure_key != null && structure_key.exists() && structure_key.canRead() && structure_key.length() != 0) {
                okeyfname = keyfname;
                n++;
                keyfname = old;
                keyfname = keyfname + "_" + Integer.toString(n);
                structure_key = new File(keyfname);
            }
            structure_key = new File(okeyfname);
        }

        molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure_xyz);
        CompositeConfiguration properties = Keyword_poltype.loadProperties(structure_key);
        ForceFieldFilter_2 forceFieldFilter = new ForceFieldFilter_2(properties, structure_key);
        forceField = forceFieldFilter.parse();
        molecularAssembly.setForceField(forceField);
        XYZFilter xyzFilter = new XYZFilter(structure_xyz, molecularAssembly, forceField, properties);
        xyzFilter.readFile();
        Utilities.biochemistry(molecularAssembly, xyzFilter.getAtomList());
        molecularAssembly.finalize(true);


        //Read options
        if (options != null) {
            if (options.toLowerCase().contains("p")) {
                do_propyze = true;
            }
//            if(options.toLowerCase().contains("d")){
//            	do_detail = true;
//            }
        }
        if (do_propyze) {
            // Get a reference to the sorted atom array.
            atoms = molecularAssembly.getAtomArray();
            nAtoms = atoms.length;


            // Define the cutoff lengths.
            double vdwOff = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
            double ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
            double buff = 2.0;
            double cutOff2 = 2.0 * (max(vdwOff, ewaldOff) + buff);

            // Determine the unit cell dimensions and Spacegroup
            String spacegroup;
            double a, b, c, alpha, beta, gamma;
            boolean aperiodic;
            try {
                a = forceField.getDouble(ForceFieldDouble.A_AXIS);
                aperiodic = false;
                b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
                c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
                alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
                beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
                gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
                spacegroup = forceField.getString(ForceFieldString.SPACEGROUP, "P1");
            } catch (Exception e) {
                logger.info(" The system will be treated as aperiodic.");
                aperiodic = true;
                spacegroup = "P1";
                /**
                 * Search all atom pairs to find the largest pair-wise distance.
                 */
                double xr[] = new double[3];
                double maxr = 0.0;
                for (int i = 0; i < nAtoms - 1; i++) {
                    double[] xi = atoms[i].getXYZ();
                    for (int j = i + 1; j < nAtoms; j++) {
                        double[] xj = atoms[j].getXYZ();
                        diff(xi, xj, xr);
                        double r = r(xr);
                        if (r > maxr) {
                            maxr = r;
                        }
                    }
                }
                a = 2.0 * (maxr + ewaldOff);
                b = a;
                c = a;
                alpha = 90.0;
                beta = 90.0;
                gamma = 90.0;
            }
            Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
            unitCell.setAperiodic(aperiodic);

            /**
             * If necessary, create a ReplicatesCrystal.
             */
            if (!aperiodic) {
                this.crystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, cutOff2);
            } else {
                this.crystal = unitCell;
            }

            vanderWaals = new VanDerWaals(molecularAssembly, crystal, parallelTeam, VDW_FORM.BUFFERED_14_7);
            pme2 = new PME_2(forceField, atoms, crystal, parallelTeam, vanderWaals.getNeighborLists(), key);
            pme2.propyze = true;
            pme2.init_prms();
        }
    }

    /**
     * <p>energy</p>
     *
     * @param gradient a boolean.
     * @param print a boolean.
     */
    public void energy(boolean gradient, boolean print) {
        ForceFieldEnergy energy = new ForceFieldEnergy(molecularAssembly);
        molecularAssembly.setPotential(energy);
//        if(do_detail){
//        	energy.tor_verbose = true;
//        }
        energy.energy(gradient, print);
        if (do_propyze) {
            system_mpoles();
        }
    }

    /**
     * <p>torsional_angles</p>
     */
    public void torsional_angles() {
    }

    /**
     * <p>system_mpoles</p>
     */
    public void system_mpoles() {
        //Find center of mass.
        double weigh = 0;
        double xyzmid[] = {0, 0, 0};
        double xyzcm[][] = new double[nAtoms][3];
        for (int i = 0; i < nAtoms; i++) {
            weigh = weigh + atoms[i].getMass();
            for (int j = 0; j < 3; j++) {
                xyzmid[j] = xyzmid[j] + atoms[i].getXYZ()[j] * atoms[i].getMass();
            }
        }
        if (weigh != 0) {
            for (int j = 0; j < 3; j++) {
                xyzmid[j] = xyzmid[j] / weigh;
            }
        }

        for (int i = 0; i < nAtoms; i++) {
            for (int j = 0; j < 3; j++) {
                xyzcm[i][j] = atoms[i].getXYZ()[j] - xyzmid[j];
            }
        }
        addInducedToGlobal();
        double netchg = 0, xdpl = 0, ydpl = 0, zdpl = 0, xxqdp = 0, xyqdp = 0, xzqdp = 0, yxqdp = 0, yyqdp = 0, yzqdp = 0, zxqdp = 0, zyqdp = 0, zzqdp = 0;
        for (int i = 0; i < nAtoms; i++) {
            double charge = atoms[i].getMultipoleType().charge;
            double[] dipole = {pme2.globalMultipole[0][i][1], pme2.globalMultipole[0][i][2], pme2.globalMultipole[0][i][3]};
            //double[] dipole = atoms[i].getMultipoleType().dipole;
            netchg = netchg + charge;
            xdpl = xdpl + xyzcm[i][0] * charge + dipole[0];
            ydpl = ydpl + xyzcm[i][1] * charge + dipole[1];
            zdpl = zdpl + xyzcm[i][2] * charge + dipole[2];
            xxqdp = xxqdp + xyzcm[i][0] * xyzcm[i][0] * charge + 2 * xyzcm[i][0] * dipole[0];
            xyqdp = xyqdp + xyzcm[i][0] * xyzcm[i][1] * charge + xyzcm[i][0] * dipole[1] + xyzcm[i][1] * dipole[0];
            xzqdp = xzqdp + xyzcm[i][0] * xyzcm[i][2] * charge + xyzcm[i][0] * dipole[2] + xyzcm[i][2] * dipole[0];

            yxqdp = yxqdp + xyzcm[i][0] * xyzcm[i][1] * charge + xyzcm[i][0] * dipole[1] + xyzcm[i][1] * dipole[0];
            yyqdp = yyqdp + xyzcm[i][1] * xyzcm[i][1] * charge + 2 * xyzcm[i][1] * dipole[1];
            yzqdp = yzqdp + xyzcm[i][1] * xyzcm[i][2] * charge + xyzcm[i][1] * dipole[2] + xyzcm[i][2] * dipole[1];

            //zxqdp = zxqdp + xyzcm[i][2] * xyzcm[i][0] * charge + xyzcm[i][2] * dipole[0] + xyzcm[i][0] * dipole[2];
            zxqdp = zxqdp + xyzcm[i][0] * xyzcm[i][2] * charge + xyzcm[i][0] * dipole[2] + xyzcm[i][2] * dipole[0];
            //zyqdp = zyqdp + xyzcm[i][2] * xyzcm[i][1] * charge + xyzcm[i][2] * dipole[1] + xyzcm[i][1] * dipole[2];
            zyqdp = zyqdp + xyzcm[i][1] * xyzcm[i][2] * charge + xyzcm[i][1] * dipole[2] + xyzcm[i][2] * dipole[1];
            zzqdp = zzqdp + xyzcm[i][2] * xyzcm[i][2] * charge + 2 * xyzcm[i][2] * dipole[2];
        }



        double qave = (xxqdp + yyqdp + zzqdp) / 3;
        xxqdp = 1.5 * (xxqdp - qave);
        xyqdp = 1.5 * xyqdp;
        xzqdp = 1.5 * xzqdp;
        yxqdp = 1.5 * yxqdp;
        yyqdp = 1.5 * (yyqdp - qave);
        yzqdp = 1.5 * yzqdp;
        zxqdp = 1.5 * zxqdp;
        zyqdp = 1.5 * zyqdp;
        zzqdp = 1.5 * (zzqdp - qave);

        for (int i = 0; i < nAtoms; i++) {
            double[][] quadrupole = {{pme2.globalMultipole[0][i][4], pme2.globalMultipole[0][i][7], pme2.globalMultipole[0][i][8]}, {pme2.globalMultipole[0][i][7], pme2.globalMultipole[0][i][5], pme2.globalMultipole[0][i][9]}, {pme2.globalMultipole[0][i][8], pme2.globalMultipole[0][i][9], pme2.globalMultipole[0][i][6]}};
            //double[][] quadrupole = atoms[i].getMultipoleType().quadrupole;
            xxqdp = xxqdp + 3 * quadrupole[0][0];
            xyqdp = xyqdp + 3 * quadrupole[0][1];
            xzqdp = xzqdp + 3 * quadrupole[0][2];
            yxqdp = yxqdp + 3 * quadrupole[1][0];
            yyqdp = yyqdp + 3 * quadrupole[1][1];
            yzqdp = yzqdp + 3 * quadrupole[1][2];
            zxqdp = zxqdp + 3 * quadrupole[2][0];
            zyqdp = zyqdp + 3 * quadrupole[2][1];
            zzqdp = zzqdp + 3 * quadrupole[2][2];
        }

        xdpl = MultipoleType.DEBYE * xdpl;
        ydpl = MultipoleType.DEBYE * ydpl;
        zdpl = MultipoleType.DEBYE * zdpl;

        xxqdp = MultipoleType.DEBYE * xxqdp;
        xyqdp = MultipoleType.DEBYE * xyqdp;
        xzqdp = MultipoleType.DEBYE * xzqdp;
        yxqdp = MultipoleType.DEBYE * yxqdp;
        yyqdp = MultipoleType.DEBYE * yyqdp;
        yzqdp = MultipoleType.DEBYE * yzqdp;
        zxqdp = MultipoleType.DEBYE * zxqdp;
        zyqdp = MultipoleType.DEBYE * zyqdp;
        zzqdp = MultipoleType.DEBYE * zzqdp;

        double netdpl = Math.sqrt(xdpl * xdpl + ydpl * ydpl + zdpl * zdpl);

        RealMatrix a = new Array2DRowRealMatrix(new double[][]{{xxqdp, xyqdp, xzqdp}, {yxqdp, yyqdp, yzqdp}, {zxqdp, zyqdp, zzqdp}});

        EigenDecompositionImpl e = new EigenDecompositionImpl(a, 1);
        a = e.getD();
        double[] netqdp = {a.getColumn(0)[0], a.getColumn(1)[1], a.getColumn(2)[2]};

        DecimalFormat myFormatter = new DecimalFormat(" ##########0.00000;-##########0.00000");
        String output;
        output = String.format(" Total Electric Charge:   %13s %s Electrons\n", " ", myFormatter.format(netchg));
        System.out.println(output);
        output = String.format(" Dipole Moment Magnitude: %13s %s Debyes\n", " ", myFormatter.format(netdpl));
        System.out.println(output);
        output = String.format(" Dipole X,Y,Z-Components: %13s %s %s %s\n", " ", myFormatter.format(xdpl), myFormatter.format(ydpl), myFormatter.format(zdpl));
        System.out.println(output);
        output = String.format(" Quadrupole Moment Tensor:%13s %s %s %s", " ", myFormatter.format(xxqdp), myFormatter.format(xyqdp), myFormatter.format(xzqdp));
        System.out.println(output);
        output = String.format("      (Buckinghams)       %13s %s %s %s", " ", myFormatter.format(yxqdp), myFormatter.format(yyqdp), myFormatter.format(yzqdp));
        System.out.println(output);
        output = String.format("                          %13s %s %s %s\n", " ", myFormatter.format(zxqdp), myFormatter.format(zyqdp), myFormatter.format(zzqdp));
        System.out.println(output);
        output = String.format("Principle Axes Quadrupole:%13s %s %s %s", " ", myFormatter.format(netqdp[2]), myFormatter.format(netqdp[1]), myFormatter.format(netqdp[0]));
        System.out.println(output);

    }

    /**
     * <p>addInducedToGlobal</p>
     */
    public void addInducedToGlobal() {
        for (int i = 0; i < nAtoms; i++) {
            for (int j = 0; j < 3; j++) {
                pme2.globalMultipole[0][i][j + 1] = pme2.globalMultipole[0][i][j + 1] + pme2.inducedDipole[0][i][j];
            }
        }
    }

    /**
     * <p>main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     * @throws java.io.IOException if any.
     */
    public static void main(String args[]) throws IOException {
        Energy e = new Energy("/users/gchattree/Research/Compounds/s_test3_compounds/famotidine/ttt.xyz", null, "d");
        e.energy(false, true);
    }
}
