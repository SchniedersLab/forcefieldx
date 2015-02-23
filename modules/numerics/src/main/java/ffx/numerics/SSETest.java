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
package ffx.numerics;

import java.util.Random;

/**
 * @author M. J. Schnieders
 */
public class SSETest {

    public static void main(String args[]) {

        int m = 500;
        int n = 500;
        if (args != null && args.length > 0) {
            m = Integer.parseInt(args[0]);
            n = Integer.parseInt(args[1]);
        }

        SSETest test = new SSETest(m, n);

        int nLoops = 100000;

        test.init(m, n);
        Random r = new Random(0);
        double temp = 0;
        long time = 0;
        for (int i = 1; i <= nLoops; i++) {
            time -= System.nanoTime();
            double y[] = test.matVec(test.A, test.x, m, n);
            time += System.nanoTime();
            int j = (int) Math.floor(m * r.nextDouble());
            temp += y[j];
            if (i % 10000 == 0) {
                System.out.println(" Nested: " + temp + " " + time * 1.0e-9);
                time = 0;
            }
        }

        test.init(m, n);
        r = new Random(0);
        temp = 0;
        time = 0;
        for (int i = 1; i <= nLoops; i++) {
            time -= System.nanoTime();
            double y[] = test.matVec(test.flatA, test.x, m, n);
            time += System.nanoTime();
            int j = (int) Math.floor(m * r.nextDouble());
            temp += y[j];
            if (i % 10000 == 0) {
                System.out.println(" Flat: " + temp + " " + time * 1.0e-9);
                time = 0;
            }
        }
    }

    public final double A[][];
    public final double flatA[];
    public final double x[];

    public SSETest(int m, int n) {
        A = new double[m][n];
        flatA = new double[m * n];
        x = new double[n];
    }

    private void init(int n, int m) {
        Random r = new Random(0);

        // Initialize A and flatA.
        for (int i = 0; i < m; i++) {
            int idx = i * n;
            for (int j = 0; j < n; j++) {
                A[i][j] = r.nextDouble();
                flatA[idx + j] = A[i][j];
            }
        }
        // Initialize X.
        for (int j = 0; j < n; j++) {
            x[j] = r.nextDouble();
        }
    }

    public final double[] matVec(final double[][] A, final double[] x,
            final int m, final int n) {
        double[] y = new double[m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                y[i] += A[i][j] * x[j];
            }
        }
        return y;
    }

    public final double[] matVec(final double A[], final double[] x,
            final int m, final int n) {
        double[] y = new double[m];
        final int extra = n - n % 8;
        final int ub = ((n / 8) * 8) - 1;
        int idx = 0;
        for (int i = 0; i < m; i++) {
            double acc = 0;
            for (int j = 0; j < ub; j += 8) {
                int ptr = idx + j;
                y[i] += A[ptr] * x[j]
                        + A[ptr + 1] * x[j + 1]
                        + A[ptr + 2] * x[j + 2]
                        + A[ptr + 3] * x[j + 3];
                acc += A[ptr + 4] * x[j + 4]
                        + A[ptr + 5] * x[j + 5]
                        + A[ptr + 6] * x[j + 6]
                        + A[ptr + 7] * x[j + 7];
            }
            y[i] += acc;
            for (int j = extra; j < n; j++) {
                y[i] += A[idx + j] * x[j];
            }
            idx += n;
        }
        return y;
    }

}
