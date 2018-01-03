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
package ffx.algorithms;

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Math.abs;

import ffx.crystal.Crystal;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;

/**
 * The Rattle classes implements the RATTLE distance constraint method.
 *
 * @author Jill J. Hauer
 */
public class Rattle {

    private static final Logger logger = Logger.getLogger(Rattle.class.getName());
    private MolecularAssembly molAss;
    private double[] xold;
    private double[] yold;
    private double[] zold;
    private double[][] vel;
    private int nVariables;
    private Crystal cryst = molAss.getCrystal();

    //set the velocity and original coordinates of global variables in Rattle Class

    public Rattle(int nVariables, MolecularAssembly molAss, double v[]) {

        int i, j = 0;
        Atom[] tempArray = molAss.getAtomArray();
        int size = tempArray.length;
        this.nVariables = nVariables;

        //get velocities prior to integrator
        for (i = 0; i < v.length; i = i + 3) {
            this.vel[j][0] = v[i];
            this.vel[j][1] = v[i + 1];
            this.vel[j][2] = v[i + 2];
            j++;
        }

        //get coordinates prior to integrator
        for (i = 0; i < size; i++) {
            this.xold[i] = tempArray[i].getX();
            this.yold[i] = tempArray[i].getY();
            this.zold[i] = tempArray[i].getZ();
        }

    }

    // Rattle algorithm applied to distances and half-step velocity values
    public void postForce1(MolecularAssembly molAss, double dt) {
        int i, n;
        int niter, maxiter;
        int ia, ib;
        int nrat; // number of bonds Rattle is applied to (all bonds)
        double eps, sor;
        double rateps;
        double xr, yr, zr;
        double xo, yo, zo;
        double dot, rma, rmb;
        double dist2;
        double delta = 0; // needs initializing to eliminate compile error on line 159
        double term;
        double xterm, yterm, zterm;
        double[] xyzr = new double[3];
        double[] xyzo = new double[3];
        boolean done;

        n = molAss.getAtomArray().length;
        nrat = molAss.getBondList().size();

        // Perform dynamic allocation of some local arrays
        boolean[] moved = new boolean[n];
        boolean[] update = new boolean[n];
        Atom[] atomArray = new Atom[n];
        double[] x = new double[n];
        double[] y = new double[n];
        double[] z = new double[n];
        double[] krat = new double[nrat]; //equilibrium bond length distances
        int[][] bondAtmNum = new int[nrat][2]; // rows correspond with bondArray index

        //initialize a list of all bonds in the system 
        for (i = 0; i < nrat; i++) {
            bondAtmNum[i][0] = molAss.getBond(i).getAtomArray()[0].getIndex();
            bondAtmNum[i][1] = molAss.getBond(i).getAtomArray()[1].getIndex();
            krat[i] = molAss.getBond(i).bondType.distance; // equilibrium distance

        }

        // Initialize the lists of atoms previously corrected
        for (i = 0; i < n; i++) {
            if (atomArray[i].isActive()) {
                moved[i] = true;
            } else {
                moved[i] = false;
            }

            update[i] = false;

        }

        // Set the iteration counter, termination and tolerance
        maxiter = 500;
        sor = 1.250;   
        rateps = 0.0000010; // default in Tinker
        eps = rateps;

        // Apply RATTLE to distances and half-step velocity values
        niter = 0;
        done = false;

        while ((!done) && (niter < maxiter)) {
            niter++;
            done = true;

            for (i = 0; i < nrat; i++) {
                ia = bondAtmNum[i][0];
                ib = bondAtmNum[i][1];
                if (moved[ia] || moved[ib]) {
                    xr = x[ib] - x[ia];
                    yr = y[ib] - y[ia];
                    zr = z[ib] - z[ia];
                    xyzr[0] = xr;
                    xyzr[1] = yr;
                    xyzr[2] = zr;

                    dist2 = cryst.image(xyzr); 
            
                    //update xr, yr, and zr
                    xr = xyzr[0];
                    yr = xyzr[1];
                    zr = xyzr[2];
                    delta = (krat[i] * krat[i]) - dist2;
                    
                    if (abs(delta) > eps) {
                        done = false;
                        update[ia] = true;
                        update[ib] = true;

                        xo = this.xold[ib] - this.xold[ia];
                        yo = this.yold[ib] - this.yold[ia];
                        zo = this.zold[ib] - this.zold[ia];
                        xyzo[0] = xo;
                        xyzo[1] = yo;
                        xyzo[2] = zo;

                        cryst.image(xyzo);
                        xo = xyzo[0];
                        yo = xyzo[1];
                        zo = xyzo[2];
                        
                        dot = (xr * xo) + (yr * yo) + (zr * zo);
                        rma = 1.00 / atomArray[ia].getMass();
                        rmb = 1.00 / atomArray[ib].getMass();
                        term = sor * delta / (2 * (rma + rmb) * dot);
                        xterm = xo * term;
                        yterm = yo * term;
                        zterm = zo * term;
                        //This is where global coordinates should be adjusted
                        x[ia] = x[ia] - (xterm * rma);
                        y[ia] = y[ia] - (yterm * rma);
                        z[ia] = z[ia] - (zterm * rma);
                        x[ib] = x[ib] - (xterm * rmb);
                        y[ib] = y[ib] - (yterm * rmb);
                        z[ib] = z[ib] - (zterm * rmb);
                        rma = rma / dt;
                        rmb = rmb / dt;
                        // This is where global velocity should be adjusted
                        this.vel[ia][0] = this.vel[ia][0] - (xterm * rma);
                        this.vel[ia][1] = this.vel[ia][1] - (yterm * rma);
                        this.vel[ia][2] = this.vel[ia][2] - (zterm * rma);
                        this.vel[ib][0] = this.vel[ib][0] - (xterm * rmb);
                        this.vel[ib][1] = this.vel[ib][1] - (yterm * rmb);
                        this.vel[ib][2] = this.vel[ib][2] - (zterm * rmb);

                    }
                }

            }

            for (i = 0; i < n; i++) {
                moved[i] = update[i];
                update[i] = false;
            }
        }

        // Write information on the number of iterations needed
        if (niter == maxiter) {
            logger.log(Level.SEVERE, "Rattle: Distance constraints not satisfied.\n");
        }
    }

    public void postForce2(double dt) {
        int i, n;
        int ia, ib;
        int niter, maxiter;
        int nrat;
        double eps, sor;
        double xr, yr, zr;
        double xv, yv, zv;
        double dot, rma, rmb;
        double term;
        double xterm, yterm, zterm;
        boolean done;
        double[] xyzr = new double[3];

        // Perform dynamic allocation of some local arrays
        n = molAss.getAtomArray().length;
        boolean[] moved = new boolean[n];
        boolean[] update = new boolean[n];
        boolean[] ratimage = new boolean[molAss.getBondList().size()]; // flag to use minimum image for holonomic constraint
        Atom[] atomArray = new Atom[n];
        int[][] bondAtmNum = new int[molAss.getBondList().size()][2];
        double[] krat = new double[molAss.getBondList().size()];
        double[] x = new double[molAss.getAtomArray().length];
        double[] y = new double[molAss.getAtomArray().length];
        double[] z = new double[molAss.getAtomArray().length];

        //initialize a list of all atom coordinates in the system
        n = molAss.getAtomArray().length;
        for (i = 0; i < n; i++) {
            x[i] = atomArray[i].getX();
            y[i] = atomArray[i].getY();
            z[i] = atomArray[i].getZ();
        }

        //initialize a list of all bonds in the system 
        nrat = molAss.getBondList().size();
        for (i = 0; i < nrat; i++) {
            bondAtmNum[i][0] = molAss.getBond(i).getAtomArray()[0].getIndex();
            bondAtmNum[i][1] = molAss.getBond(i).getAtomArray()[1].getIndex();
            krat[i] = molAss.getBond(i).bondType.distance; // equilibrium distance

        }

        //initialize ratimage for minimum image convention
        for (i = 0; i < nrat; i++) 
        {
            ratimage[i] = true;
        }

        //initialize lists of atoms previously corrected
        atomArray = molAss.getAtomArray();
        for (i = 0; i < n; i++) {
            if (atomArray[i].isActive()) {
                moved[i] = true;
            } else {
                moved[i] = false;
            }

            update[i] = false;
        }

        //set iteration counter, termination, and tolerance
        maxiter = 500;
        niter = 0;
        done = false;
        sor = 1.25;
        eps = 0.000001 / dt;

        //apply the RATTLE algorithm to correct velocities
        while ((!done) && (niter < maxiter)) {
            niter++;
            done = true;
            for (i = 1; i < nrat; i++) 
            {
                ia = bondAtmNum[i][0];
                ib = bondAtmNum[i][1];

                if (moved[ia] || moved[ib]) {
                    xr = x[ib] - x[ia];
                    yr = y[ib] - y[ia];
                    zr = z[ib] - z[ia];
                    
                    xyzr[0] = xr;
                    xyzr[1] = yr;
                    xyzr[2] = zr;
                    
                    cryst.image(xyzr);
                    xr = xyzr[0];
                    yr = xyzr[1];
                    zr = xyzr[2];

                    xv = this.vel[ib][0] - this.vel[ia][0];
                    yv = this.vel[ib][1] - this.vel[ia][1];
                    zv = this.vel[ib][2] - this.vel[ia][2];
                    dot = (xr * xv) + (yr * yv) + (zr * zv);
                    rma = 1.00 / atomArray[ia].getMass();
                    rmb = 1.00 / atomArray[ib].getMass();
                    term = -dot / ((rma + rmb) * (krat[i] * krat[i]));

                    if (abs(term) > eps) {
                        done = false;
                        update[ia] = true;
                        update[ib] = true;
                        term = sor * term;
                        xterm = xr * term;
                        yterm = yr * term;
                        zterm = zr * term;
                        this.vel[ia][0] = this.vel[ia][0] - (xterm * rma);
                        this.vel[ia][1] = this.vel[ia][1] - (yterm * rma);
                        this.vel[ia][2] = this.vel[ia][2] - (zterm * rma);
                        this.vel[ib][0] = this.vel[ib][0] + (xterm * rma);
                        this.vel[ib][1] = this.vel[ib][1] + (yterm * rma);
                        this.vel[ib][2] = this.vel[ib][2] + (zterm * rma);

                    }

                }

            }
            for(i = 0; i < n; i++)
            {
                moved[i] = update[i];
                update[i] = false;
            }
        }

        // Write information on the number of iterations needed
        if (niter == maxiter) {
            logger.log(Level.SEVERE, "Rattle: Velocity constraints not satisfied.\n");
        }

    }
}
