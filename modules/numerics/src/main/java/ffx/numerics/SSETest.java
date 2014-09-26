package ffx.numerics;

import java.util.Random;

/**
 * @author M. J. Schnieders
 */
public class SSETest {

    public static void main(String args[]) {

        int m = Integer.parseInt(args[0]);
        int n = Integer.parseInt(args[1]);

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
