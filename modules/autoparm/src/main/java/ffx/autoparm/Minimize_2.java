/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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

import java.io.*;
import java.text.DecimalFormat;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration.CompositeConfiguration;

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.Terminatable;
import ffx.numerics.LBFGS;
import ffx.numerics.LineSearch.LineSearchResult;
import ffx.numerics.OptimizationListener;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.XYZFilter;

/**
 * Minimize the potential energy of a system to an RMS gradient per atom
 * convergence criteria.
 *
 * @author Gaurav Chattree and Michael J. Schnieders
 * @since 1.0
 *
 */
public class Minimize_2 implements OptimizationListener, Terminatable {

    private static final Logger logger = Logger.getLogger(Minimize_2.class.getName());
    private int n;
    private final double[] x;
    private final double[] grad;
    private final double[] scaling;
    private MolecularAssembly molecularAssembly;
    private final Potential potential;
    private AlgorithmListener algorithmListener;
    private boolean done = false;
    private boolean terminate = false;
    private long time;
    private double grms;
    private int nSteps;
    private File structure_key;
    private File structure_xyz;
    InputStreamReader stdinput = new InputStreamReader(System.in);
    private Atom atoms[];
    private String name;
    private int add;
    BufferedReader stdreader = new BufferedReader(stdinput);
    private ForceField forceField;

    /**
     * <p>Constructor for Minimize_2.</p>
     *
     * @param xyz_filename a {@link java.lang.String} object.
     * @param keyfname a {@link java.lang.String} object.
     * @throws java.io.IOException if any.
     */
    public Minimize_2(String xyz_filename, String keyfname) throws IOException {

        structure_xyz = new File(xyz_filename);
        if (!(structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead())) {
            System.out.println("Couldn't find xyz file");
            System.exit(1);
        }
        int r = 1;
        String oxyzfname = null;
        String old = xyz_filename;
        while (structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead()) {
            oxyzfname = xyz_filename;
            r++;
            xyz_filename = old;
            xyz_filename = xyz_filename + "_" + Integer.toString(r);
            structure_xyz = new File(xyz_filename);
        }
        add = r;
        structure_xyz = new File(oxyzfname);
        int index = oxyzfname.lastIndexOf(".");
        name = oxyzfname.substring(0, index);

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
        atoms = molecularAssembly.getAtomArray();
        molecularAssembly.finalize(true);
        //algorithmListener = this;
        if (molecularAssembly.getPotentialEnergy() == null) {
            molecularAssembly.setPotential(new ForceFieldEnergy(molecularAssembly));
        }
        potential = molecularAssembly.getPotentialEnergy();
        n = potential.getNumberOfVariables();
        x = new double[n];
        grad = new double[n];
        scaling = new double[n];
        for (int i = 0; i < n; i++) {
            scaling[i] = 12.0;
        }
        potential.setScaling(scaling);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void terminate() {
        terminate = true;
        while (!done) {
            synchronized (this) {
                try {
                    wait(1);
                } catch (Exception e) {
                    logger.log(Level.WARNING, "Exception terminating minimization.\n", e);
                }
            }
        }
    }

    /**
     * <p>minimize</p>
     *
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize() {
        return minimize(1.0);
    }

    /**
     * <p>minimize</p>
     *
     * @param eps a double.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(double eps) {
        return minimize(7, eps);
    }

    /**
     * <p>minimize</p>
     *
     * @param m a int.
     * @param eps a double.
     * @return a {@link ffx.numerics.Potential} object.
     */
    public Potential minimize(int m, double eps) {
        time = System.nanoTime();
        potential.getCoordinates(x);
        /**
         * Scale coordinates.
         */
        for (int i = 0; i < n; i++) {
            x[i] *= scaling[i];
        }

        done = false;
        int status = 2;
        //print();
        double e = potential.energyAndGradient(x, grad);
        //print();
        status = LBFGS.minimize(n, m, x, e, grad, eps, potential, this);
        print();
        done = true;

        switch (status) {
            case 0:
                logger.info(String.format("\n Optimization achieved convergence criteria: %8.5f\n", grms));
                break;
            case 1:
                logger.info(String.format("\n Optimization terminated at step %d.\n", nSteps));
                break;
            default:
                logger.warning("\n Optimization failed.\n");
        }
        return potential;
    }

    /**
     * <p>print</p>
     */
    public void print() {
//    	for(int i = 0; i < x.length; i+=3){
//    		System.out.println(x[i]/12+" "+x[i+1]/12+" "+x[i+2]/12);
//    	}
//    	double norm = 0;
//    	for(int i = 0; i < grad.length; i++){
//    		norm += Math.pow(grad[i],2);
//    	}
//    	norm = Math.sqrt(norm);
//    	System.out.println("gradnorm: "+norm);
//    	System.out.println(potential.energyAndGradient(x, grad));
        File outf = new File(name + ".xyz_" + add);
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outf)));
            DecimalFormat myFormatter = new DecimalFormat(" ##########0.00000;-##########0.00000");

            bw.write(String.format("%6d\n", atoms.length));
            for (Atom a : atoms) {
                String output = String.format("%6d", a.xyzIndex) + "  " + a.getAtomType().name + " " + String.format("%12s %12s %12s", myFormatter.format(a.getX()), myFormatter.format(a.getY()), myFormatter.format(a.getZ())) + " " + String.format("%6d", a.getAtomType().atomClass);
                for (int i = 0; i < a.getBonds().size(); i++) {
                    output += String.format("%6d", a.getBonds().get(i).get1_2(a).xyzIndex);
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
     * {@inheritDoc}
     *
     * Implement the OptimizationListener interface.
     *
     * @since 1.0
     */
    @Override
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df, double angle, LineSearchResult info) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        if (iter == 0) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time\n");
        }
        if (info == null) {
            logger.info(String.format("%6d%13.4f%11.4f", iter, f, grms));
        } else {
            if (info == LineSearchResult.Success) {
                logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8.3f",
                        iter, f, grms, df, xrms, angle, nfun, seconds));
            } else {
                logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%9.2f%7d %8s",
                        iter, f, grms, df, xrms, angle, nfun, info.toString()));
            }
        }
        // Update the listener and check for an termination request.
//        if (algorithmListener != null) {
//            algorithmListener.algorithmUpdate(molecularAssembly);
//        }
        if (terminate) {
            logger.info(" The optimization recieved a termination request.");
            // Tell the L-BFGS optimizer to terminate.
            return false;
        }
        return true;
    }

    /**
     * Implement the OptimizationListener interface.
     *
     * @param iter Number of iterations so far.
     * @param nfun Number of function evaluations so far.
     * @param grms Gradient RMS at current solution.
     * @param xrms Coordinate change RMS at current solution.
     * @param f Function value at current solution.
     * @param df Change in the function value compared to the previous solution.
     * @since 1.0
     * @return a boolean.
     */
    public boolean optimizationUpdate(int iter, int nfun, double grms, double xrms, double f, double df) {
        long currentTime = System.nanoTime();
        Double seconds = (currentTime - time) * 1.0e-9;
        time = currentTime;
        this.grms = grms;
        this.nSteps = iter;

        if (iter == 1) {
            logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
            logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Evals     Time\n");
        }
        logger.info(String.format("%6d%13.4f%11.4f%11.4f%10.4f%7d %8.3f",
                iter, f, grms, df, xrms, nfun, seconds));
        // Update the listener and check for an termination request.
        if (algorithmListener != null) {
            algorithmListener.algorithmUpdate(molecularAssembly);
        }
        if (terminate) {
            logger.info(" The optimization recieved a termination request.");
            // Tell the optimizer to terminate.
            return false;
        }
        return true;
    }

    /**
     * <p>main</p>
     *
     * @param args an array of {@link java.lang.String} objects.
     */
    public static void main(String args[]) {
        try {
            Minimize_2 m = new Minimize_2("/home/gchattree/Research/Compounds/s_test3_compounds/famotidine/ttt.xyz", null);
            m.minimize(0.1);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
