/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.potential.bonded;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.exp;

import ffx.potential.MolecularAssembly;
import ffx.potential.parsers.PDBFilter;

/**
 * @author Mallory R. Tollefson
 */
public class SturmMethod {

    private static final Logger logger = Logger.getLogger(SturmMethod.class.getName());

    final int PRINT_LEVEL = 0;
    final int MAX_ORDER = 16;
    final int MAXPOW = 32;
    final double SMALL_ENOUGH = 1.0e-18;
    double RELERROR;
    int MAXIT, MAX_ITER_SECANT;
    double[] roots;
    double[][] xyz_o = new double[5][3];

    /**
     * Class that holds polynomial information.
     *
     * Poly objects have one position that holds the order of a polynomial in
     * poly.ord.
     *
     * Poly objects have up to 16 positions that hold the coefficients of a 16th
     * degree polynomial.
     *
     * The coefficients are located in poly.coef[position].
     */
    private class Poly {

        int ord;
        double coef[];

        private Poly() {
            coef = new double[MAX_ORDER + 1];
        }
    }

    /**
     * Sets termination criteria for polynomial solver.
     *
     * @param tol_secant
     * @param max_iter_sturm
     * @param max_iter_secant
     */
    void initializeSturm(double[][] tol_secant, int[][] max_iter_sturm, int[][] max_iter_secant) {
        RELERROR = tol_secant[0][0];
        MAXIT = max_iter_sturm[0][0];
        MAX_ITER_SECANT = max_iter_secant[0][0];
    }

    void solveSturm(int[] p_order, int[] n_root, double[] poly_coeffs, double[] roots) {

        Poly[] sseq = new Poly[MAX_ORDER * 2];
        double min, max;
        int order, nroots, nchanges, np;
        int[] atmin = new int[1];
        int[] atmax = new int[1];
        this.roots = roots;

        for (int i = 0; i < MAX_ORDER * 2; i++) {
            sseq[i] = new Poly();
        }

        order = p_order[0];

        for (int i = order; i >= 0; i--) {
            sseq[0].coef[i] = poly_coeffs[i];
        }

        if (logger.isLoggable(Level.FINE)) {
            StringBuilder string = new StringBuilder();
            for (int i = order; i >= 0; i--) {
                string.append(String.format(" Coefficients in Sturm solver\n"));
                string.append(String.format("%d %f\n", i, sseq[0].coef[i]));
            }
            logger.fine(string.toString());
        }

        np = buildSturm(order, sseq);

        if (logger.isLoggable(Level.FINE)) {
            StringBuilder string1 = new StringBuilder();
            string1.append(String.format(" Sturm sequence for:\n"));
            for (int i = order; i >= 0; i--) {
                string1.append(String.format("%f ", sseq[0].coef[i]));
                string1.append("\n");
            }
            for (int i = 0; i <= np; i++) {
                for (int j = sseq[i].ord; j >= 0; j--) {
                    string1.append(String.format("%f ", sseq[i].coef[j]));
                    string1.append("\n");
                }
            }
            logger.fine(string1.toString());
        }

        //get the number of real roots
        nroots = numRoots(np, sseq, atmin, atmax);

        if (nroots == 0) {
            n_root[0] = nroots;
        }

        if (logger.isLoggable(Level.FINE)) {
            logger.fine(String.format(" Number of real roots: %d\n", nroots));
        }

        //calculate the bracket that the roots live in
        min = -1.0;
        nchanges = numChanges(np, sseq, min);

        for (int i = 0; nchanges != atmin[0] && i != MAXPOW; i++) {
            min *= 10.0;
            nchanges = numChanges(np, sseq, min);
        }

        if (nchanges != atmin[0]) {
            logger.fine(String.format(" Solve: unable to bracket all negative roots\n"));
            atmin[0] = nchanges;
        }

        max = 1.0;
        nchanges = numChanges(np, sseq, max);

        for (int i = 0; nchanges != atmax[0] && i != MAXPOW; i++) {
            max *= 10.0;
            nchanges = numChanges(np, sseq, max);
        }
        if (nchanges != atmax[0]) {
            logger.fine(String.format(" Solve: unable to bracket all positive roots\n"));
            atmax[0] = nchanges;
        }

        // Perform the bisection
        nroots = atmin[0] - atmax[0];
        sbisect(np, sseq, min, max, atmin[0], atmax[0], this.roots);
        n_root[0] = nroots;

        //write out the roots
        if (logger.isLoggable(Level.FINE)) {
            if (nroots == 1) {
                logger.fine(String.format("\n One distinct real root at x = %f\n", this.roots[0]));
            } else {
                StringBuilder string2 = new StringBuilder();
                string2.append(String.format("\n %d distinct real roots for x: \n", nroots));
                for (int i = 0; i != nroots; i++) {
                    string2.append(String.format("%f\n", this.roots[i]));
                }
                logger.fine(string2.toString());
            }
        }
    }

    double hyperTan(double a, double x) {
        double ax = a * x;
        if (ax > 100.0) {
            return (1.0);
        } else if (ax < -100.0) {
            return (-1.0);
        } else {
            double exp_x1 = exp(ax);
            double exp_x2 = exp(-ax);
            return (exp_x1 - exp_x2) / (exp_x1 + exp_x2);
        }
    }

    /**
     * Calculates the modulus of u(x)/v(x) leaving it in r, it returns 0 if r(x)
     * is constant. This function assumes the leading coefficient of v is is 1
     * or -1.
     *
     * Modp was originally a static int and returned r.ord as an int. The
     * buildSturm function requires that a boolean is returned from the modp
     * function, so modp is set to return a boolean. This new boolean return
     * could be problematic elsewhere in the code if modp is used to return a
     * value. This should be checked.
     *
     * @param u
     * @param v
     * @param r
     * @return
     */
    static boolean modp(Poly u, Poly v, Poly r) {

        double nr[] = r.coef;
        int end = u.ord;
        double uc[] = u.coef;

        for (int i = 0; i <= end; i++) {
            nr[i] = uc[i];

        }
        if (v.coef[v.ord] < 0.0) {
            for (int k = u.ord - v.ord - 1; k >= 0; k -= 2) {
                r.coef[k] = -r.coef[k];
            }

            for (int k = u.ord - v.ord; k >= 0; k--) {
                for (int j = v.ord + k - 1; j >= k; j--) {
                    r.coef[j] = -r.coef[j] - r.coef[v.ord + k] * v.coef[j - k];
                }
            }
        } else {
            for (int k = u.ord - v.ord; k >= 0; k--) {
                for (int j = v.ord + k - 1; j >= k; j--) {
                    r.coef[j] -= r.coef[v.ord + k] * v.coef[j - k];
                }
            }
        }

        int k = v.ord - 1;
        while (k >= 0 && abs(r.coef[k]) < 1.0e-18) {
            r.coef[k] = 0.0;
            k--;
        }

        if (k < 0) {
            r.ord = 0;
        } else {
            r.ord = k;
        }

        if (r.ord > 0) {
            return true;
        } else if (r.ord == 0) {
            return false;
        } else {
            return false;
        }
    }

    /**
     * Build up a sturm sequence for a polynomial, and return the number of
     * polynomials in the sequence.
     *
     * @param ord
     * @param sseq
     * @return
     */
    int buildSturm(int ord, Poly[] sseq) {
        double f;
        double[] fp;
        double[] fc;

        sseq[0].ord = ord;
        sseq[1].ord = ord - 1;

        // Calculate the derivative and normalise the leading coefficient
        f = abs(sseq[0].coef[ord] * ord);
        fp = sseq[1].coef;
        fc = Arrays.copyOfRange(sseq[0].coef, 1, sseq[0].coef.length);

        int j = 0;
        for (int i = 1; i <= ord; i++) {
            fp[j] = fc[j] * i / f;
            j++;
        }

        // Construct the rest of the Sturm sequence. (Double check this... )
        int i;
        for (i = 0; i < sseq[0].coef.length - 2 && modp(sseq[i], sseq[i + 1], sseq[i + 2]); i++) {
            //reverse the sign and normalise
            f = -Math.abs(sseq[i + 2].coef[sseq[i + 2].ord]);
            for (j = sseq[i + 2].ord; j >= 0; j--) {
                sseq[i + 2].coef[j] /= f;
            }
        }

        // Reverse the sign.
        sseq[i + 2].coef[0] = -sseq[i + 2].coef[0];

        return (sseq[0].ord - sseq[i + 2].ord);
    }

    /**
     * Return the number of distinct real roots of the polynomial described in
     * sseq.
     *
     * @param np
     * @param sseq
     * @param atneg
     * @param atpos
     * @return
     */
    int numRoots(int np, Poly[] sseq, int[] atneg, int[] atpos) {

        int atposinf = 0;
        int atneginf = 0;

        // Changes at positive infinity.
        double lf = sseq[0].coef[sseq[0].ord];

        for (int i = 1; i <= np; i++) //for (s = sseq + 1; s <= sseq + np; s++)
        {
            double f = sseq[i].coef[sseq[i].ord];
            if (lf == 0.0 || lf * f < 0) {
                atposinf++;
            }
            lf = f;
        }

        // Changes at negative infinity.
        if ((sseq[0].ord & 1) != 0) {
            lf = -sseq[0].coef[sseq[0].ord];
        } else {
            lf = sseq[0].coef[sseq[0].ord];
        }

        for (int i = 1; i <= np; i++) {
            double f;
            if ((sseq[i].ord & 1) != 0) {
                f = -sseq[i].coef[sseq[i].ord];
            } else {
                f = sseq[i].coef[sseq[i].ord];
            }
            if (lf == 0.0 || lf * f < 0) {
                atneginf++;
            }
            lf = f;
        }

        atneg[0] = atneginf;
        atpos[0] = atposinf;

        return (atneginf - atposinf);
    }

    /**
     * Return the number of sign changesin the Sturm sequence in sseq at the
     * value a.
     *
     * @param np
     * @param sseq
     * @param a
     * @return
     */
    int numChanges(int np, Poly[] sseq, double a) {

        int changes = 0;
        double lf = evalpoly(sseq[0].ord, sseq[0].coef, a);
        for (int i = 1; i <= np; i++) {
            double f = evalpoly(sseq[i].ord, sseq[i].coef, a);
            if (lf == 0.0 || lf * f < 0) {
                changes++;
            }

            lf = f;
        }

        return changes;
    }

    /**
     * Uses a bisection based on the Sturm sequence for the polynomial described
     * in sseq to isolate intervals in which roots occur, the roots are returned
     * in the roots array in order of magnitude.
     *
     * @param np
     * @param sseq
     * @param min
     * @param max
     * @param atmin
     * @param atmax
     * @param roots
     */
    void sbisect(int np, Poly[] sseq, double min, double max, int atmin, int atmax, double[] roots) {
        double mid = 0.;
        int n1 = 0;
        int n2 = 0;
        int its, atmid, nroot;

        int remainder = 0;
        int count = 0;

        nroot = atmin - atmax;

        if (nroot == 1) {
            //first try a less expensive technique
            if (modrf(sseq[0].ord, sseq[0].coef, min, max, roots)) {
                remainder = this.roots.length - roots.length;
                count = 0;
                for (int q = remainder; q < this.roots.length; q++) {
                    this.roots[q] = roots[count];
                    count++;
                }
                return;

            }

            //When code reaches this point, the root must be evaluated
            //      using the Sturm sequence.
            for (its = 0; its < MAXIT; its++) {
                mid = (min + max) / 2;
                atmid = numChanges(np, sseq, mid);
                if (abs(mid) > RELERROR) {
                    if (abs((max - min) / mid) < RELERROR) {
                        roots[0] = mid;
                        remainder = this.roots.length - roots.length;
                        count = 0;
                        for (int q = remainder; q < this.roots.length; q++) {
                            this.roots[q] = roots[count];
                            count++;
                        }
                        return;
                    }
                } else if (abs(max - min) < RELERROR) {
                    roots[0] = mid;
                    remainder = this.roots.length - roots.length;
                    count = 0;
                    for (int q = remainder; q < this.roots.length; q++) {
                        this.roots[q] = roots[count];
                        count++;
                    }
                    return;
                }
                if ((atmin - atmid) == 0) {
                    min = mid;
                } else {
                    max = mid;
                }
            }
            if (its == MAXIT) {
                logger.info(String.format(" sbisect: overflow min %f max %f diff %f nroot %d n1 %d n2 %d\n",
                        min, max, max - min, nroot, n1, n2));
                roots[0] = mid;
            }
            remainder = this.roots.length - roots.length;
            count = 0;
            for (int q = remainder; q < this.roots.length; q++) {
                this.roots[q] = roots[count];
                count++;
            }
            return;

        }

        // More than one root in the interval, must bisect.
        for (its = 0; its < MAXIT; its++) {
            mid = (min + max) / 2;
            atmid = numChanges(np, sseq, mid);
            n1 = atmin - atmid;
            n2 = atmid - atmax;
            if (n1 != 0 && n2 != 0) {
                sbisect(np, sseq, min, mid, atmin, atmid, roots);
                sbisect(np, sseq, mid, max, atmid, atmax, Arrays.copyOfRange(roots, n1, roots.length));
                break;
            }
            if (n1 == 0) {
                min = mid;
            } else {
                max = mid;
            }
        }

        if (its == MAXIT) {
            for (n1 = atmax; n1 < atmin; n1++) {
                roots[n1 - atmax] = mid;
            }
        }

        //may not need this~~~~
        remainder = this.roots.length - roots.length;
        count = 0;
        for (int q = remainder; q < this.roots.length; q++) {
            this.roots[q] = roots[count];
            count++;
        }
    }

    /**
     * Evaluate polynomial defined in coef returning its value.
     *
     * @param ord
     * @param coef
     * @param x
     * @return
     */
    double evalpoly(int ord, double[] coef, double x) {

        double fp[] = coef;
        double f = fp[ord];

        int i = ord;
        for (i--; i >= 0; i--) {
            f = x * f + fp[i];
        }
        return f;
    }

    /**
     * Uses the modified regula-falsi method to evaluate the root in interval
     * [a,b] of the polynomial described in coef. The root is returned in val.
     * The routine returns zero if it can't converge.
     *
     * @param ord
     * @param coef
     * @param a
     * @param b
     * @param val
     * @return
     */
    boolean modrf(int ord, double[] coef, double a, double b, double[] val) {
        double fa = coef[ord];
        double fb = fa;

        for (int i = ord - 1; i >= 0; i--) {
            fa = a * fa + coef[i];
            fb = b * fb + coef[i];
        }
        if (fa * fb > 0.0) {
            return false;
        }
        double lfx = fa;
        for (int its = 0; its < MAX_ITER_SECANT; its++) {
            double x = (fb * a - fa * b) / (fb - fa);
            // constrain that x stays in the bounds
            if (x < a || x > b) {
                x = 0.5 * (a + b);
            }
            double fx = coef[ord];
            for (int i = ord - 1; i >= 0; i--) {
                fx = x * fx + coef[i];
            }
            if (abs(x) > RELERROR) {
                if (abs(fx / x) < RELERROR) {
                    val[0] = x;
                    return true;
                }
            } else if (abs(fx) < RELERROR) {
                val[0] = x;
                return true;
            }
            if ((fa * fx) < 0) {
                b = x;
                fb = fx;
                if ((lfx * fx) > 0) {
                    fa /= 2;
                }
            } else {
                a = x;
                fa = fx;
                if ((lfx * fx) > 0) {
                    fb /= 2;
                }
            }

            lfx = fx;
        }

        return false;
    }

    /**
     * Write out loop coordinates and determine oxygen placement.
     *
     * @param pdb_name
     * @param res_name
     * @param r_n
     * @param r_a
     * @param r_c
     * @param stt_res
     * @param end_res
     * @param chain_n
     * @param chain_a
     * @param chain_c
     * @param molAss
     * @param counter
     * @param writeFile
     * @return
     */
    File writePDBBackbone(char[] pdb_name, String[] res_name, double[][] r_n, double[][] r_a, double[][] r_c, int stt_res, int end_res, char[] chain_n, char[] chain_a, char[] chain_c, MolecularAssembly molAss, int counter, boolean writeFile) {

        Polymer[] newChain = molAss.getChains();
        ArrayList<Atom> backBoneAtoms;
        double[] xyz_n = new double[3];
        double[] xyz_a = new double[3];
        double[] xyz_c = new double[3];

        xyz_o[0][0] = 0.0;
        xyz_o[0][1] = 0.0;
        xyz_o[0][2] = 0.0;
        xyz_o[4][0] = 0.0;
        xyz_o[4][1] = 0.0;
        xyz_o[4][2] = 0.0;

        ArrayList<Atom> OAtoms = new ArrayList<>();

        for (int i = stt_res + 1; i < end_res; i++) {
            Residue newResidue = newChain[0].getResidue(i);
            backBoneAtoms = newResidue.getBackboneAtoms();
            for (Atom backBoneAtom : backBoneAtoms) {
                switch (backBoneAtom.getAtomType().name) {
                    case "C":
                        xyz_c[0] = r_c[i - stt_res][0];
                        xyz_c[1] = r_c[i - stt_res][1];
                        xyz_c[2] = r_c[i - stt_res][2];
                        backBoneAtom.moveTo(xyz_c);
                        break;
                    case "N":
                        xyz_n[0] = r_n[i - stt_res][0];
                        xyz_n[1] = r_n[i - stt_res][1];
                        xyz_n[2] = r_n[i - stt_res][2];
                        backBoneAtom.moveTo(xyz_n);
                        break;
                    case "CA":
                        xyz_a[0] = r_a[i - stt_res][0];
                        xyz_a[1] = r_a[i - stt_res][1];
                        xyz_a[2] = r_a[i - stt_res][2];
                        backBoneAtom.moveTo(xyz_a);
                        break;
                    case "HA":
                        newResidue.deleteAtom(backBoneAtom);
                        break;
                    case "H":
                        newResidue.deleteAtom(backBoneAtom);
                        break;
                    case "O":
                        OAtoms.add(backBoneAtom);
                        break;
                    default:
                        newResidue.deleteAtom(backBoneAtom);
                        break;
                }
            }

            ArrayList<Atom> sideChainAtoms = newResidue.getSideChainAtoms();
            for (Atom sideChainAtom : sideChainAtoms) {
                newResidue.deleteAtom(sideChainAtom);
            }
        }

        int oCount = 0;
        for (int i = stt_res + 1; i < end_res; i++) {
            Residue newResidue = newChain[0].getResidue(i);
            backBoneAtoms = newResidue.getBackboneAtoms();
            Atom CA = new Atom("CA");
            Atom N = new Atom("N");
            Atom C = new Atom("C");
            Atom O = OAtoms.get(oCount);

            for (Atom backBoneAtom : backBoneAtoms) {
                switch (backBoneAtom.getAtomType().name) {
                    case "C":
                        C = backBoneAtom;
                        break;
                    case "N":
                        N = backBoneAtom;
                        break;
                    case "CA":
                        CA = backBoneAtom;
                        break;
                    default:
                        break;
                }
            }
            BondedUtils.intxyz(O, C, 1.2255, CA, 122.4, N, 180, 0);
            xyz_o[i - stt_res][0] = O.getX();
            xyz_o[i - stt_res][1] = O.getY();
            xyz_o[i - stt_res][2] = O.getZ();

            oCount++;
        }

        File file = molAss.getFile();

        /**
         * for (int i = stt_res; i <= end_res; i++) { double[] xyz =
         * molAss.getBackBoneAtoms().get(i).getXYZ(); }
         */
        String filename = FilenameUtils.removeExtension(file.getAbsolutePath());
        if (!filename.contains("_loop")) {
            filename = filename + "_loop";
        }

        File modifiedFile = new File(filename + ".pdb_" + counter);
        PDBFilter modFilter = new PDBFilter(modifiedFile, molAss, null, null);

        if (writeFile) {
            modFilter.writeFile(modifiedFile, true);
        }

        return (modifiedFile);
    }

    /**
     * Used only in JUnit testing.
     *
     * @return
     */
    public double[][] getr_o() {
        return xyz_o;
    }
}
