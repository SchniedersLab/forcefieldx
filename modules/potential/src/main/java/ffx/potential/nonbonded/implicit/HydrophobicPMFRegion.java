//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
package ffx.potential.nonbonded.implicit;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.tanh;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;

import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Torsion;

/**
 * Initial implementation of a Hydrophobic PMF.
 *
 * @since 1.0
 */
public class HydrophobicPMFRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(HydrophobicPMFRegion.class.getName());

    private Atom[] atoms;
    private int nAtoms;
    private double[] x, y, z;
    private boolean[] use;
    private double[][][] grad;

    // Radius of a carbon atom.
    private final double rCarbon = 1.7;
    // Radius of a water molecule.
    private final double rWater = 1.4;
    // Constant for calculation of atomic surface area.
    private final double safact = 0.3516;
    // Surface area of a hydrophobic carbon atom.
    private final double acSurf = 120.7628;
    // tanh slope (set very steep).
    private final double tSlope = 100.0;
    // Shift the tanh plot along the x-axis.
    private final double tOffset = 6.0;
    // Cutoff distance for pairwise HPMF interactions.
    private final double hpmfCut = 11.0;
    // Cutoff squared
    private final double hpmfCut2 = hpmfCut * hpmfCut;
    // Hydrophobic PMF well depth parameter.
    private final double h1 = -0.7308004860404441194;
    private final double h2 = 0.2001645051578760659;
    private final double h3 = -0.0905499953418473502;
    // Hydrophobic PMF well center point.
    private final double c1 = 3.8167879266271396155;
    private final double c2 = 5.4669162286016419472;
    private final double c3 = 7.1167694861385353278;
    // Reciprocal of the hydrophobic PMF well width.
    private final double w1 = 1.6858993102248638341;
    private final double w2 = 1.3906405621629980285;
    private final double w3 = 1.5741657341338335385;
    private final double rSurf = rCarbon + 2.0 * rWater;
    private final double piSurf = PI * (rCarbon + rWater);
    // Radius of each atom for use with hydrophobic PMF.
    private final double[] rPMF;
    // Number of hydrophobic carbon atoms in the system.
    private final int nCarbon;
    // Number of the atom for each HPMF carbon atom site.
    private final int[] iCarbon;
    // SASA value for each hydrophobic PMF carbon atom
    private final double[] carbonSASA;
    // SASA value
    private final double[] tanhSA;
    // Loop to find the SASA value of each hydrophobic carbon.
    private final CarbonSASALoop[] carbonSASALoop;
    // Loop to find the hydrophobic energy.
    private final HydrophobicPMFLoop[] hydrophobicPMFLoop;
    // Loop to find the SASA chain rule derivatives.
    private final CarbonSASACRLoop[] carbonSASACRLoop;
    // Shared energy variable.
    private final SharedDouble sharedEnergy;
    private boolean gradient;
    private final double[] dtanhSA;
    private final double[] sasa;
    private final double[] carbonSASACR;

    public HydrophobicPMFRegion(Atom[] atoms, double[] x, double[] y, double[] z,
                                boolean[] use, double[][][] grad, int nt) {
        logger.info(format(" Hydrophobic PMF cut-off:              %8.2f (A)", hpmfCut));

        this.atoms = atoms;
        nAtoms = atoms.length;
        this.x = x;
        this.y = y;
        this.z = z;
        this.use = use;
        this.grad = grad;

        /**
         * Count hydrophobic carbons.
         */
        int count = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            int atomicNumber = atom.getAtomicNumber();
            if (atomicNumber == 6) {
                List<Bond> bonds = atom.getBonds();
                int bondCount = bonds.size();
                if (bondCount <= 2) {
                    continue;
                }
                boolean keep = true;
                for (int j = 0; j < bondCount; j++) {
                    Atom atom2 = bonds.get(j).get1_2(atom);
                    int atomicNumber2 = atom2.getAtomicNumber();
                    if (bondCount == 3 && atomicNumber2 == 8) {
                        keep = false;
                        break;
                    }
                }
                if (keep) {
                    count++;
                }
            }
        }
        /**
         * Allocate arrays.
         */
        rPMF = new double[nAtoms];
        nCarbon = count;
        iCarbon = new int[nCarbon];
        carbonSASA = new double[nCarbon];
        carbonSASACR = new double[nCarbon];
        tanhSA = new double[nCarbon];
        dtanhSA = new double[nCarbon];
        sasa = new double[nCarbon];

        carbonSASALoop = new CarbonSASALoop[nt];
        hydrophobicPMFLoop = new HydrophobicPMFLoop[nt];
        carbonSASACRLoop = new CarbonSASACRLoop[nt];
        for (int i = 0; i < nt; i++) {
            carbonSASALoop[i] = new CarbonSASALoop();
            hydrophobicPMFLoop[i] = new HydrophobicPMFLoop();
            carbonSASACRLoop[i] = new CarbonSASACRLoop();
        }
        sharedEnergy = new SharedDouble();

        /**
         * Assign hydrophobic carbon values.
         */
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            int atomicNumber = atom.getAtomicNumber();
            if (atomicNumber == 6) {
                int nh = 0;
                List<Bond> bonds = atom.getBonds();
                int bondCount = bonds.size();
                if (bondCount <= 2) {
                    continue;
                }
                boolean keep = true;
                for (int j = 0; j < bondCount; j++) {
                    Atom atom2 = bonds.get(j).get1_2(atom);
                    int atomicNumber2 = atom2.getAtomicNumber();
                    if (atomicNumber2 == 1) {
                        nh++;
                    }
                    if (bondCount == 3 && atomicNumber2 == 8) {
                        keep = false;
                    }
                }
                if (keep) {
                    iCarbon[index] = i;
                    carbonSASA[index] = 1.0;
                    if (bondCount == 3 && nh == 0) {
                        carbonSASA[index] = 1.554;
                    } else if (bondCount == 3 && nh == 1) {
                        carbonSASA[index] = 1.073;
                    } else if (bondCount == 4 && nh == 1) {
                        carbonSASA[index] = 1.276;
                    } else if (bondCount == 4 && nh == 2) {
                        carbonSASA[index] = 1.045;
                    } else if (bondCount == 4 && nh == 3) {
                        carbonSASA[index] = 0.880;
                    }
                    carbonSASA[index] = carbonSASA[index] * safact / acSurf;
                    if (logger.isLoggable(Level.FINEST)) {
                        logger.finest(format(" %d Base HPMF SASA for atom %d: %10.8f",
                                index + 1, i + 1, carbonSASA[index]));
                    }
                    index++;
                }
            }
        }

        /**
         * Assign HPMF atomic radii from traditional Bondi values
         */
        for (int i = 0; i < nAtoms; i++) {
            rPMF[i] = 2.0;
            int atmnum = atoms[i].getAtomicNumber();
            switch (atmnum) {
                case 0:
                    rPMF[i] = 0.0;
                    break;
                case 1:
                    rPMF[i] = 1.20;
                    break;
                case 2:
                    rPMF[i] = 1.40;
                    break;
                case 5:
                    rPMF[i] = 1.80;
                    break;
                case 6:
                    rPMF[i] = 1.70;
                    break;
                case 7:
                    rPMF[i] = 1.55;
                    break;
                case 8:
                    rPMF[i] = 1.50;
                    break;
                case 9:
                    rPMF[i] = 1.47;
                    break;
                case 10:
                    rPMF[i] = 1.54;
                    break;
                case 14:
                    rPMF[i] = 2.10;
                    break;
                case 15:
                    rPMF[i] = 1.80;
                    break;
                case 16:
                    rPMF[i] = 1.80;
                    break;
                case 17:
                    rPMF[i] = 1.75;
                    break;
                case 18:
                    rPMF[i] = 1.88;
                    break;
                case 34:
                    rPMF[i] = 1.90;
                    break;
                case 35:
                    rPMF[i] = 1.85;
                    break;
                case 36:
                    rPMF[i] = 2.02;
                    break;
                case 53:
                    rPMF[i] = 1.98;
                    break;
                case 54:
                    rPMF[i] = 2.16;
                    break;
            }
        }
    }

    public void setGradient(boolean gradient) {
        this.gradient = gradient;
    }

    public double getEnergy() {
        return sharedEnergy.get();
    }

    @Override
    public void start() {
        sharedEnergy.set(0);
    }

    @Override
    public void run() {
        int ti = getThreadIndex();
        try {
            execute(0, nCarbon - 1, carbonSASALoop[ti]);
            execute(0, nCarbon - 2, hydrophobicPMFLoop[ti]);
            if (gradient) {
                execute(0, nCarbon - 1, carbonSASACRLoop[ti]);
            }
        } catch (Exception e) {
            String message = "Fatal exception computing Born radii in thread " + ti + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * Compute Hydrophobic PMF radii.
     *
     * @since 1.0
     */
    private class CarbonSASALoop extends IntegerForLoop {

        @Override
        public void run(int lb, int ub) {
            /**
             * Get the surface area for each hydrophobic carbon atom.
             */
            for (int ii = lb; ii <= ub; ii++) {
                final int i = iCarbon[ii];
                if (!use[i]) {
                    continue;
                }
                final double carbonSA = carbonSASA[ii];
                double sa = acSurf;
                final double xi = x[i];
                final double yi = y[i];
                final double zi = z[i];
                int count = 0;
                for (int k = 0; k < nAtoms; k++) {
                    if (i != k && use[k]) {
                        final double xr = x[k] - xi;
                        final double yr = y[k] - yi;
                        final double zr = z[k] - zi;
                        final double r2 = xr * xr + yr * yr + zr * zr;
                        double rk = rPMF[k];
                        double rBig = rk + rSurf;
                        if (r2 < rBig * rBig) {
                            final double r = sqrt(r2);
                            final double rSmall = rk - rCarbon;
                            final double part = piSurf * (rBig - r) * (1.0 + rSmall / r);
                            sa *= (1.0 - carbonSA * part);
                            count++;
                        }
                    }
                }
                sasa[ii] = sa;
                //sasa[ii] = carbonSA;
                double tSA = tanh(tSlope * (sa - tOffset));
                tanhSA[ii] = 0.5 * (1.0 + tSA);
                dtanhSA[ii] = 0.5 * tSlope * (1.0 - tSA * tSA);
            }
        }
    }

    /**
     * Compute Born radii for a range of atoms via the Grycuk method.
     *
     * @since 1.0
     */
    private class HydrophobicPMFLoop extends IntegerForLoop {

        private double energy;
        // Omit
        private final int[] omit;
        private double[] gX;
        private double[] gY;
        private double[] gZ;

        public HydrophobicPMFLoop() {
            omit = new int[nAtoms];
        }

        @Override
        public void start() {
            energy = 0.0;
            for (int i = 0; i < nAtoms; i++) {
                omit[i] = -1;
            }
            int threadID = getThreadIndex();
            gX = grad[threadID][0];
            gY = grad[threadID][1];
            gZ = grad[threadID][2];
            if (gradient) {
                for (int i = 0; i < nCarbon; i++) {
                    carbonSASACR[i] = 0.0;
                }
            }
        }

        @Override
        public void run(int lb, int ub) {
            /**
             * Hydrophobic PME energy.
             */
            for (int ii = lb; ii <= ub; ii++) {
                final int i = iCarbon[ii];
                if (!use[i]) {
                    continue;
                }
                double tanhSAi = tanhSA[ii];
                Atom ssAtom = null;
                Atom atom = atoms[i];
                List<Bond> bonds = atom.getBonds();
                for (Bond bond : bonds) {
                    Atom atom2 = bond.get1_2(atom);
                    int k = atom2.getIndex() - 1;
                    if (!use[k]) {
                        continue;
                    }
                    omit[k] = i;
                    if (atom2.getAtomicNumber() == 16) {
                        ssAtom = atom2;
                    }
                }
                List<Angle> angles = atom.getAngles();
                for (Angle angle : angles) {
                    Atom atom2 = angle.get1_3(atom);
                    if (atom2 != null) {
                        int k = atom2.getIndex() - 1;
                        if (!use[k]) {
                            continue;
                        }
                        omit[k] = i;
                    }
                }
                List<Torsion> torsions = atom.getTorsions();
                for (Torsion torsion : torsions) {
                    Atom atom2 = torsion.get1_4(atom);
                    if (atom2 != null) {
                        int k = atom2.getIndex() - 1;
                        if (!use[k]) {
                            continue;
                        }
                        omit[k] = i;
                        if (ssAtom != null) {
                            List<Bond> bonds2 = atom2.getBonds();
                            for (Bond bond : bonds2) {
                                Atom s = bond.get1_2(atom2);
                                if (s.getAtomicNumber() == 16) {
                                    List<Bond> sBonds = s.getBonds();
                                    for (Bond sBond : sBonds) {
                                        Atom s2 = sBond.get1_2(s);
                                        if (s2 == ssAtom) {
                                            omit[k] = -1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                final double xi = x[i];
                final double yi = y[i];
                final double zi = z[i];
                double e = 0.0;
                for (int kk = ii + 1; kk < nCarbon; kk++) {
                    int k = iCarbon[kk];
                    if (!use[k]) {
                        continue;
                    }
                    if (omit[k] != i) {
                        final double xr = xi - x[k];
                        final double yr = yi - y[k];
                        final double zr = zi - z[k];
                        final double r2 = xr * xr + yr * yr + zr * zr;
                        if (r2 < hpmfCut2) {
                            final double r = sqrt(r2);
                            final double a1 = (r - c1) * w1;
                            final double a2 = (r - c2) * w2;
                            final double a3 = (r - c3) * w3;
                            final double e1 = h1 * exp(-a1 * a1);
                            final double e2 = h2 * exp(-a2 * a2);
                            final double e3 = h3 * exp(-a3 * a3);
                            final double t1t2 = tanhSAi * tanhSA[kk];
                            final double sum = (e1 + e2 + e3);
                            e += sum * t1t2;
                            if (gradient) {
                                /**
                                 * First part of hydrophobic PMF derivative
                                 * calculation.
                                 */
                                double de1 = -2.0 * e1 * a1 * w1;
                                double de2 = -2.0 * e2 * a2 * w2;
                                double de3 = -2.0 * e3 * a3 * w3;
                                double dsum = (de1 + de2 + de3) * t1t2 / r;
                                double dedx = dsum * xr;
                                double dedy = dsum * yr;
                                double dedz = dsum * zr;
                                gX[i] += dedx;
                                gY[i] += dedy;
                                gZ[i] += dedz;
                                gX[k] -= dedx;
                                gY[k] -= dedy;
                                gZ[k] -= dedz;
                                /**
                                 * Chain Rule Term.
                                 */
                                carbonSASACR[ii] += sum * tanhSA[kk] * dtanhSA[ii];
                                carbonSASACR[kk] += sum * tanhSA[ii] * dtanhSA[kk];
                            }
                        }
                    }
                }
                energy += e;
            }
        }

        @Override
        public void finish() {
            sharedEnergy.addAndGet(energy);
        }
    }

    /**
     * Compute Hydrophobic PMF chain rule term.
     *
     * @since 1.0
     */
    private class CarbonSASACRLoop extends IntegerForLoop {

        private double[] gX;
        private double[] gY;
        private double[] gZ;

        public CarbonSASACRLoop() {
        }

        @Override
        public void start() {
            int threadID = getThreadIndex();
            gX = grad[threadID][0];
            gY = grad[threadID][1];
            gZ = grad[threadID][2];
        }

        @Override
        public void run(int lb, int ub) {
            for (int ii = lb; ii <= ub; ii++) {
                final int i = iCarbon[ii];
                if (!use[i]) {
                    continue;
                }
                final double carbonSA = carbonSASA[ii];
                final double xi = x[i];
                final double yi = y[i];
                final double zi = z[i];
                for (int k = 0; k < nAtoms; k++) {
                    if (i != k && use[k]) {
                        final double xr = xi - x[k];
                        final double yr = yi - y[k];
                        final double zr = zi - z[k];
                        final double r2 = xr * xr + yr * yr + zr * zr;
                        double rk = rPMF[k];
                        double rBig = rk + rSurf;
                        if (r2 <= rBig * rBig) {
                            final double r = sqrt(r2);
                            final double rSmall = rk - rCarbon;
                            final double rr = 1.0 / r;
                            final double rr2 = rr * rr;
                            final double part = piSurf * (rBig - r) * (1.0 + rSmall * rr);
                            double t1b = -piSurf * (1.0 + rBig * rSmall * rr2);
                            double t1a = -sasa[ii] / (1.0 / carbonSA - part);
                            double de = t1a * t1b * rr * carbonSASACR[ii];
                            double dedx = de * xr;
                            double dedy = de * yr;
                            double dedz = de * zr;
                            gX[i] += dedx;
                            gY[i] += dedy;
                            gZ[i] += dedz;
                            gX[k] -= dedx;
                            gY[k] -= dedy;
                            gZ[k] -= dedz;
                        }
                    }
                }
            }
        }
    }
}
