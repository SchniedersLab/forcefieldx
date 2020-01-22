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

import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.atan2;
import static org.apache.commons.math3.util.FastMath.cbrt;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.reduction.SharedDouble;

import ffx.numerics.switching.MultiplicativeSwitch;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.utils.EnergyException;
import static ffx.numerics.math.VectorMath.cross;
import static ffx.numerics.math.VectorMath.dist;
import static ffx.numerics.math.VectorMath.dist2;
import static ffx.numerics.math.VectorMath.dot;
import static ffx.numerics.math.VectorMath.norm;
import static ffx.numerics.math.VectorMath.r;

/**
 * Initial port of the TINKER volume area code.
 */
public class VolumeRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(VolumeRegion.class.getName());

    private final int nAtoms;
    private final Atom[] atoms;
    private final double[] x, y, z;
    /**
     * Base atomic radii without the "exclude" buffer.
     */
    private final double[] baseRadius;
    /**
     * Atomic radii including the "exclude" buffer.
     */
    private final double[] radius;
    /**
     * If true, atom is not used.
     */
    private final boolean[] skip;

    private final VolumeLoop[] volumeLoop;
    private final SharedDouble sharedVolume;
    private final SharedDouble sharedArea;
    /**
     * Surface area (Ang^2).
     */
    private double surfaceArea;
    /**
     * Surface area energy (kcal/mol).
     */
    private double surfaceAreaEnergy;
    /**
     * Volume (Ang^3).
     */
    private double volume;
    /**
     * Volume energy (kcal/mol).
     */
    private double volumeEnergy;
    /**
     * Cavitation energy, which is a function of volume and/or surface area.
     */
    private double cavitationEnergy;
    /**
     * Probe is used for molecular (contact/reentrant) volume and surface area.
     */
    private double probe = 0.0;
    /**
     * Exclude is used for excluded volume and accessible surface area.
     */
    private double exclude = 1.4;
    /**
     * Solvent pressure for small solutes.
     */
    private double solventPressure = GeneralizedKirkwood.DEFAULT_SOLVENT_PRESSURE;
    /**
     * Chandler cross-over point between small solutes and large solutes.
     */
    private double crossOver = GeneralizedKirkwood.DEFAULT_CROSSOVER;
    /**
     * Surface tension for large solutes.
     */
    private double surfaceTension = GeneralizedKirkwood.DEFAULT_CAVDISP_SURFACE_TENSION;

    /**
     * Effective radius probe.
     * <p>
     * In cavitation volume scaling regime, approximate solvent excluded volume
     * and effective radius are computed as follow.
     * <p>
     * 1) GaussVol vdW volume is computed from defined radii.
     * 2) Effective radius is computed as Reff = cbrt(3.0 * volume / (4.0 * PI)) + effectiveRadiusProbe.
     * 3) Solvent Excluded Volume = 4/3 * Pi * Reff^3
     */
    private double effectiveRadius;
    private double switchRange = 3.5;
    private double saSwitchRangeOff = 3.9;
    /**
     * Begin turning off the Volume term.
     */
    private double volumeOff = crossOver - switchRange;
    /**
     * Volume term is zero at the cut-off.
     */
    private double volumeCut = crossOver + switchRange;
    /**
     * Begin turning off the SA term.
     */
    private double surfaceAreaOff = crossOver + saSwitchRangeOff;
    /**
     * SA term is zero at the cut-off.
     */
    private double surfaceAreaCut = crossOver - switchRange;
    /**
     * Volume multiplicative switch.
     */
    private MultiplicativeSwitch volumeSwitch = new MultiplicativeSwitch(volumeCut, volumeOff);
    /**
     * Surface area multiplicative switch.
     */
    private MultiplicativeSwitch surfaceAreaSwitch = new MultiplicativeSwitch(surfaceAreaCut, surfaceAreaOff);

    private final int[] itab;
    private final static int MAXCUBE = 40;
    private final static int MAXARC = 1000;
    private final static int MAXMNB = 500;
    /**
     * maximum number of cycle convex edges.
     */
    private final static int MAXCYEP = 30;
    /**
     * maximum number of convex face cycles.
     */
    private final static int MAXFPCY = 10;
    /**
     * maximum number of saddle faces.
     */
    private final int maxfs;
    /**
     * maximum number of convex edges.
     */
    private final int maxep;
    /**
     * maximum number of neighboring atom pairs.
     */
    private final int maxcls;
    /**
     * maximum number of circles.
     */
    private final int maxc;
    /**
     * maximum number of total tori.
     */
    private final int maxt;
    /**
     * maximum number of temporary tori.
     */
    private final int maxtt;
    /**
     * maximum number of concave edges.
     */
    private final int maxen;
    /**
     * maximum number of vertices.
     */
    private final int maxv;
    /**
     * maximum number of probe positions.
     */
    private final int maxp;
    /**
     * maximum number of concave faces.
     */
    private final int maxfn;
    /**
     * maximum number of convex faces.
     */
    private final int maxfp;
    /**
     * maximum number of cycles.
     */
    private final int maxcy;
    /**
     * Copy of the atomic coordinates. [X,Y,Z][Atom Number]
     */
    private final double[][] a;
    /**
     * If true, atom has no free surface.
     */
    private final boolean[] nosurf;
    /**
     * Atom free of neighbors.
     */
    private final boolean[] afree;
    /**
     * Atom buried.
     */
    private final boolean[] abur;
    /**
     * Begin and end pointers for atoms neighbors.
     */
    private final int[][] acls;
    /**
     * Atom numbers of neighbors.
     */
    private final int[] cls;
    /**
     * Pointer from neighbor to torus.
     */
    private final int[] clst;
    /**
     * Number of temporary tori.
     */
    private int ntt;
    /**
     * Temporary torus atom numbers.
     */
    private final int[][] tta;
    /**
     * First edge of each temporary torus.
     */
    private final int[] ttfe;
    /**
     * Last edge of each temporary torus.
     */
    private final int[] ttle;
    /**
     * Pointer to next edge of temporary torus.
     */
    private final int[] enext;
    /**
     * Temporary torus buried.
     */
    private final boolean[] ttbur;
    /**
     * Temporary torus free.
     */
    private final boolean[] ttfree;
    /**
     * Torus center.
     */
    private final double[][] t;
    /**
     * Torus radius.
     */
    private final double[] tr;
    /**
     * Torus axis.
     */
    private final double[][] tax;
    /**
     * Number of tori.
     */
    private int nt;
    /**
     * Torus atom number.
     */
    private final int[][] ta;
    /**
     * Torus first edge.
     */
    private final int[] tfe;
    /**
     * Torus free edge of neighbor.
     */
    private final boolean[] tfree;
    /**
     * Probe position coordinates.
     */
    private final double[][] p;
    /**
     * Number of probe positions.
     */
    private int np;
    /**
     * Probe position atom numbers.
     */
    private final int[][] pa;
    /**
     * Vertex coordinates.
     */
    private final double[][] v;
    /**
     * Number of vertices.
     */
    private int nv;
    /**
     * Vertex atom number.
     */
    private final int[] va;
    /**
     * Vertex probe number.
     */
    private final int[] vp;
    /**
     * Number of concave edges.
     */
    private int nen;
    /**
     * Vertex numbers for each concave edge.
     */
    private final int[][] env;
    /**
     * Concave face concave edge numbers.
     */
    private final int[][] fnen;
    /**
     * Circle center.
     */
    private final double[][] c;
    /**
     * Circle radius.
     */
    private final double[] cr;
    /**
     * Number of circles.
     */
    private int nc;
    /**
     * Circle atom number.
     */
    private final int[] ca;
    /**
     * Circle torus number.
     */
    private final int[] ct;
    /**
     * Number of convex edges.
     */
    private int nep;
    /**
     * Convex edge circle number.
     */
    private final int[] epc;
    /**
     * Convex edge vertex number.
     */
    private final int[][] epv;
    /**
     * First convex edge of each atom.
     */
    private final int[] afe;
    /**
     * Last convex edge of each atom.
     */
    private final int[] ale;
    /**
     * Pointer to next convex edge of atom.
     */
    private final int[] epnext;
    /**
     * Number of saddle faces.
     */
    private int nfs;
    /**
     * Saddle face concave edge numbers.
     */
    private final int[][] fsen;
    /**
     * Saddle face convex edge numbers.
     */
    private final int[][] fsep;
    /**
     * Number of cycles.
     */
    private int ncy;
    /**
     * Number of convex edges in cycle.
     */
    private final int[] cynep;
    /**
     * Cycle convex edge numbers.
     */
    private final int[][] cyep;
    /**
     * Number of convex faces.
     */
    private int nfp;
    /**
     * Atom number of convex face.
     */
    private final int[] fpa;
    /**
     * Convex face cycle numbers
     */
    private final int[][] fpcy;
    /**
     * Number of cycles bounding convex face.
     */
    private final int[] fpncy;

    // These are from the "nearby" method.
    /**
     * True if cube contains active atoms.
     */
    private final boolean[] activeCube;
    /**
     * True if cube or adjacent cubes have active atoms.
     */
    private final boolean[] activeAdjacentCube;
    /**
     * Pointer to first atom in list for cube.
     */
    private final int[] firstAtomPointer;
    /**
     * Integer cube coordinates.
     */
    private final int[][] cubeCoordinates;
    /**
     * Pointer to next atom in cube.
     */
    private final int[] nextAtomPointer;

    private final static double SIZE = 0.000001;
    private final double[] vector = new double[3];

    /**
     * VolumeRegion constructor.
     *
     * @param atoms Array of atom instances.
     * @param x     X-coordinates.
     * @param y     Y-coordinates.
     * @param z     Z-cooddinates.
     * @param nt    Number of threads.
     */
    public VolumeRegion(Atom[] atoms, double[] x, double[] y, double[] z,
                        double[] baseRadius, int nt) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.x = x;
        this.y = y;
        this.z = z;
        this.baseRadius = baseRadius;

        volumeLoop = new VolumeLoop[nt];
        for (int i = 0; i < nt; i++) {
            volumeLoop[i] = new VolumeLoop();
        }
        sharedVolume = new SharedDouble();
        sharedArea = new SharedDouble();
        itab = new int[nAtoms];

        radius = new double[nAtoms];
        skip = new boolean[nAtoms];
        maxfs = 6 * nAtoms;
        maxep = 12 * nAtoms;
        maxcls = 240 * nAtoms;
        maxc = 8 * nAtoms;
        maxt = 4 * nAtoms;
        maxtt = 120 * nAtoms;
        maxen = 12 * nAtoms;
        maxv = 12 * nAtoms;
        maxp = 4 * nAtoms;
        maxfn = 4 * nAtoms;
        maxfp = 2 * nAtoms;
        maxcy = 3 * nAtoms;

        a = new double[3][nAtoms];
        nosurf = new boolean[nAtoms];
        afree = new boolean[nAtoms];
        abur = new boolean[nAtoms];
        acls = new int[2][nAtoms];
        cls = new int[maxcls];
        clst = new int[maxcls];
        tta = new int[2][maxtt];
        ttfe = new int[maxtt];
        ttle = new int[maxtt];
        enext = new int[maxen];
        ttbur = new boolean[maxtt];
        ttfree = new boolean[maxtt];
        t = new double[3][maxt];
        tr = new double[maxt];
        tax = new double[3][maxt];

        ta = new int[2][maxt];
        tfe = new int[maxt];
        tfree = new boolean[maxt];
        p = new double[3][maxp];
        pa = new int[3][maxp];
        v = new double[3][maxv];
        va = new int[maxv];
        vp = new int[maxv];
        env = new int[2][maxen];
        fnen = new int[3][maxfn];
        c = new double[3][maxc];
        cr = new double[maxc];
        ca = new int[maxc];
        ct = new int[maxc];
        epc = new int[maxep];
        epv = new int[2][maxep];
        afe = new int[nAtoms];
        ale = new int[nAtoms];
        epnext = new int[maxep];
        fsen = new int[2][maxfs];
        fsep = new int[2][maxfs];
        cynep = new int[maxcy];
        cyep = new int[MAXCYEP][maxcy];
        fpa = new int[maxfp];
        fpcy = new int[MAXFPCY][maxfp];
        fpncy = new int[maxfp];

        activeCube = new boolean[MAXCUBE * MAXCUBE * MAXCUBE];
        activeAdjacentCube = new boolean[MAXCUBE * MAXCUBE * MAXCUBE];
        firstAtomPointer = new int[MAXCUBE * MAXCUBE * MAXCUBE];
        cubeCoordinates = new int[3][nAtoms];
        nextAtomPointer = new int[nAtoms];
    }

    public void setSolventPressure(double solventPressure) {
        this.solventPressure = solventPressure;
    }

    public void setSurfaceTension(double surfaceTension) {
        this.surfaceTension = surfaceTension;
    }

    public void setCrossOver(double crossOver) {
        this.crossOver = crossOver;
        volumeOff = crossOver - switchRange;
        volumeCut = crossOver + switchRange;
        surfaceAreaOff = crossOver + saSwitchRangeOff;
        surfaceAreaCut = crossOver - switchRange;
        volumeSwitch = new MultiplicativeSwitch(volumeCut, volumeOff);
        surfaceAreaSwitch = new MultiplicativeSwitch(surfaceAreaCut, surfaceAreaOff);
    }

    public void setExclude(double exclude) {
        this.exclude = exclude;
    }

    public void setProbe(double probe) {
        this.probe = probe;
    }

    public double getSurfaceArea() {
        return surfaceArea;
    }

    public double getVolume() {
        return volume;
    }

    public double getSolventPressure() {
        return solventPressure;
    }

    public double getSurfaceTension() {
        return surfaceTension;
    }

    public double getCrossOver() { return crossOver; }

    /**
     * Return Volume based cavitation energy.
     *
     * @return Volume based cavitation energy.
     */
    public double getVolumeEnergy() {
        return volumeEnergy;
    }

    /**
     * Return Surface Area based cavitation energy.
     *
     * @return Surface Area based cavitation energy.
     */
    public double getSurfaceAreaEnergy() {
        return surfaceAreaEnergy;
    }

    public double getEnergy() {
        return cavitationEnergy;
    }

    public double getEffectiveRadius() {
        return effectiveRadius;
    }

    public double getProbe() {
        return probe;
    }

    public double getExclude() {
        return exclude;
    }

    @Override
    public void start() {
        sharedVolume.set(0.0);
        sharedArea.set(0.0);
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            a[0][i] = atom.getX();
            a[1][i] = atom.getY();
            a[2][i] = atom.getZ();
        }
        wiggle();
    }

    @Override
    public void run() {
        try {
            execute(0, nAtoms - 1, volumeLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = "Fatal exception computing Volume energy in thread " + getThreadIndex() + "\n";
            logger.log(Level.SEVERE, message, e);
        }
    }

    @Override
    public void finish() {
        energy();
    }

    private void wiggle() {
        // Apply a small perturbation of fixed magnitude to each atom.
        for (int i = 0; i < nAtoms; i++) {
            getRandomVector(vector);
            a[0][i] = x[i] + (SIZE * vector[0]);
            a[1][i] = y[i] + (SIZE * vector[1]);
            a[2][i] = z[i] + (SIZE * vector[2]);
        }
    }

    private void getRandomVector(double[] vector) {
        double x, y, s;
        x = 0;
        y = 0;

        // Get a pair of appropriate components in the plane.
        s = 2.0;
        while (s >= 1.0) {
            x = (2.0 * Math.random()) - 1.0;
            y = (2.0 * Math.random()) - 1.0;
            s = (x * x) + (y * y);
        }

        // Construct the 3-dimensional random unit vector.
        vector[2] = 1.0 - 2.0 * s;
        s = 2.0 * sqrt(1.0 - s);
        vector[1] = s * y;
        vector[0] = s * x;
    }

    private void energy() {

        // Calculate a purely surface area based cavitation energy.
        surfaceArea = sharedArea.get();
        surfaceAreaEnergy = surfaceArea * surfaceTension;

        // Calculate a purely volume based cavitation energy.
        volume = sharedVolume.get();
        volumeEnergy = volume * solventPressure;

        // Use Volume to find an effective cavity radius.
        effectiveRadius = cbrt(3.0 * volume / (4.0 * PI));
        double reff = effectiveRadius;
        double reff2 = reff * reff;
        double reff3 = reff2 * reff;
        double reff4 = reff3 * reff;
        double reff5 = reff4 * reff;
        double vdWVolPI23 = pow(volume / PI, 2.0 / 3.0);
        // double dReffdvdW = 1.0 / (pow(6.0, 2.0 / 3.0) * PI * vdWVolPI23);

        // Find the cavitation energy using a combination of volume and surface area dependence.
        if (reff < volumeOff) {
            // Find cavity energy from only the molecular volume.
            cavitationEnergy = volumeEnergy;
            // addVolumeGradient(dVolSEVdVolvdW * solventPressure, gradient);
        } else if (reff <= volumeCut) {
            // Include a tapered molecular volume.
            double taper = volumeSwitch.taper(reff, reff2, reff3, reff4, reff5);
            cavitationEnergy = taper * volumeEnergy;
            // double dtaper = volumeSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
            // double factor = dtaper * volumeEnergy + taper * solventPressure;
            // addVolumeGradient(factor, gradient);
        }

        if (reff > surfaceAreaOff) {
            // Find cavity energy from only SA.
            cavitationEnergy = surfaceAreaEnergy;
            // addSurfaceAreaGradient(0.0, surfaceTension, gradient);
        } else if (reff >= surfaceAreaCut) {
            // Include a tapered surface area term.
            double taperSA = surfaceAreaSwitch.taper(reff, reff2, reff3, reff4, reff5);
            cavitationEnergy += taperSA * surfaceAreaEnergy;
            // double dtaperSA = surfaceAreaSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
            // addSurfaceAreaGradient(dtaperSA * surfaceAreaEnergy, taperSA * surfaceTension, gradient);
        }

        if (logger.isLoggable(Level.FINE)) {
            logger.fine(format("\n Volume:              %8.3f (Ang^3)", volume));
            logger.fine(format(" Volume Energy:       %8.3f (kcal/mol)", volumeEnergy));
            logger.fine(format(" Surface Area:        %8.3f (Ang^2)", surfaceArea));
            logger.fine(format(" Surface Area Energy: %8.3f (kcal/mol)", surfaceAreaEnergy));
            logger.fine(format(" Volume + SA Energy:  %8.3f (kcal/mol)", cavitationEnergy));
            logger.fine(format(" Effective Radius:    %8.3f (Ang)", reff));
        }
    }

    /**
     * Row major indexing (the last dimension is contiguous in memory).
     *
     * @param i x index.
     * @param j y index.
     * @param k z index.
     * @return the row major index.
     */
    private static int index(int i, int j, int k) {
        return k + MAXCUBE * (j + MAXCUBE * i);
    }

    /**
     * Compute Volume energy for a range of atoms.
     *
     * @since 1.0
     */
    private class VolumeLoop extends IntegerForLoop {

        private int nfn;
        private int iv;
        private int itemp;
        private final int mxcube = 15;
        private final int[] inov = new int[MAXARC];
        private final int[][] cube = new int[2][mxcube * mxcube * mxcube];
        private double xmin, ymin, zmin;
        private double xmax, ymax, zmax;
        private double theta1;
        private double theta2;
        private double pix2;
        private double rmax;
        private final double[] arci = new double[MAXARC];
        private final double[] arcf = new double[MAXARC];
        private final double[] dx = new double[MAXARC];
        private final double[] dy = new double[MAXARC];
        private final double[] dsq = new double[MAXARC];
        private final double[] d = new double[MAXARC];
        private boolean ttok;
        private final double[] vdwrad = new double[nAtoms];
        private final double[][] dex = new double[3][nAtoms];
        private double localVolume;
        private double localSurfaceArea;

        /**
         * Extra padding to avert cache interface.
         */
        private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        private long pad8, pad9, pada, padb, padc, padd, pade, padf;

        /**
         * Row major indexing (the last dimension is contiguous in memory).
         *
         * @param i x index.
         * @param j y index.
         * @param k z index.
         * @return the row major index.
         */
        private int index2(int i, int j, int k) {
            return k + mxcube * (j + mxcube * i);
        }

        @Override
        public void start() {
            fill(dex[0], 0.0);
            fill(dex[1], 0.0);
            fill(dex[2], 0.0);
            localVolume = 0.0;
            localSurfaceArea = 0.0;
        }

        @Override
        public void finish() {
            sharedVolume.addAndGet(localVolume);
            sharedArea.addAndGet(localSurfaceArea);
        }

        void setRadius() {
            // Initialize minimum and maximum range of atoms.
            pix2 = 2.0 * PI;
            rmax = 0.0;
            xmin = x[0];
            xmax = x[0];
            ymin = y[0];
            ymax = y[0];
            zmin = z[0];
            zmax = z[0];

            /*
              Assign van der Waals radii to the atoms; note that the radii
              are incremented by the size of the probe; then get the
              maximum and minimum ranges of atoms.
             */
            for (int i = 0; i < nAtoms; i++) {
                if (baseRadius[i] == 0.0) {
                    radius[i] = 0.0;
                    vdwrad[i] = 0.0;
                    skip[i] = true;
                } else {
                    skip[i] = false;
                    vdwrad[i] = baseRadius[i] + exclude;
                    radius[i] = baseRadius[i] + exclude;
                    if (vdwrad[i] > rmax) {
                        rmax = vdwrad[i];
                    }
                    if (x[i] < xmin) {
                        xmin = x[i];
                    }
                    if (x[i] > xmax) {
                        xmax = x[i];
                    }
                    if (y[i] < ymin) {
                        ymin = y[i];
                    }
                    if (y[i] > ymax) {
                        ymax = y[i];
                    }
                    if (z[i] < zmin) {
                        zmin = z[i];
                    }
                    if (z[i] > zmax) {
                        zmax = z[i];
                    }
                }
            }
        }

        /**
         * Find the analytical volume and surface area.
         */
        void calcVolume() {
            nearby();
            torus();
            place();
            compress();
            saddles();
            contact();
            vam();
        }

        void getVector(double[] ai, double[][] temp, int index) {
            ai[0] = temp[0][index];
            ai[1] = temp[1][index];
            ai[2] = temp[2][index];
        }

        void getVector(double[] ai, double[][][] temp, int index1, int index2) {
            ai[0] = temp[0][index1][index2];
            ai[1] = temp[1][index1][index2];
            ai[2] = temp[2][index1][index2];
        }

        /**
         * The gettor method tests for a possible torus position at the
         * interface between two atoms, and finds the torus radius, center
         * and axis.
         */
        boolean gettor(int ia, int ja, double[] torcen, double[] torad, double[] torax) {

            double dij, temp;
            double temp1, temp2;
            double[] vij = new double[3];
            double[] uij = new double[3];
            double[] bij = new double[3];
            double[] ai = new double[3];
            double[] aj = new double[3];

            // Get the distance between the two atoms.
            ttok = false;
            getVector(ai, a, ia);
            getVector(aj, a, ja);
            dij = dist(ai, aj);

            // Find a unit vector along interatomic (torus) axis.
            for (int k = 0; k < 3; k++) {
                vij[k] = a[k][ja] - a[k][ia];
                uij[k] = vij[k] / dij;
            }

            // Find coordinates of the center of the torus.
            temp = 1.0 + ((radius[ia] + probe) * (radius[ia] + probe) - (radius[ja] + probe) * (radius[ja] + probe)) / (dij * dij);
            for (int k = 0; k < 3; k++) {
                bij[k] = a[k][ia] + 0.5 * vij[k] * temp;
            }

            // Skip if atoms too far apart (should not happen).
            temp1 = (radius[ia] + radius[ja] + 2.0 * probe) * (radius[ia] + radius[ja] + 2.0 * probe) - dij * dij;
            if (temp1 >= 0.0) {

                // Skip if one atom is inside the other.
                temp2 = dij * dij - (radius[ia] - radius[ja]) * (radius[ia] - radius[ja]);
                if (temp2 >= 0.0) {

                    // Store the torus radius, center, and axis.
                    ttok = true;
                    torad[0] = sqrt(temp1 * temp2) / (2.0 * dij);
                    for (int k = 0; k < 3; k++) {
                        torcen[k] = bij[k];
                        torax[k] = uij[k];
                    }
                }
            }
            return ttok;
        }

        /**
         * The nearby method finds all of the through-space neighbors of
         * each atom for use in surface area and volume calculations.
         */
        public void nearby() {
            int maxclsa = 1000;
            int iptr, juse;
            int i1, j1, k1;
            int iatom, jatom;
            int ici, icj, ick;
            int jci, jcj, jck;
            int jcls, jmin;
            int jmincls = 0;
            int jmold;
            int ncls, nclsa;
            int[] clsa = new int[maxclsa];

            // Temporary neighbor list, before sorting.
            int[] tempNeighborList = new int[maxclsa];
            double radmax, width;
            double sum, sumi;
            double d2, r2;
            double vect1, vect2, vect3;

            // Minimum atomic coordinates (cube corner).
            double[] minAtomicCoordinates = new double[3];
            double[] ai = new double[3];
            double[] aj = new double[3];

            /*
             * Ignore all atoms that are completely inside another atom;
             * may give nonsense results if this step is not taken.
             */
            for (int i = 0; i < nAtoms - 1; i++) {
                if (!skip[i]) {
                    getVector(ai, a, i);
                    for (int j = i + 1; j < nAtoms; j++) {
                        getVector(aj, a, j);
                        d2 = dist2(ai, aj);
                        r2 = (radius[i] - radius[j]) * (radius[i] - radius[j]);
                        if (!skip[j] && d2 < r2) {
                            if (radius[i] < radius[j]) {
                                skip[i] = true;
                            } else {
                                skip[j] = true;
                            }
                        }
                    }
                }
            }

            // Check for new coordinate minima and radii maxima.
            radmax = 0.0;
            for (int k = 0; k < 3; k++) {
                minAtomicCoordinates[k] = a[k][0];
            }
            for (int i = 0; i < nAtoms; i++) {
                for (int k = 0; k < 3; k++) {
                    if (a[k][i] < minAtomicCoordinates[k]) {
                        minAtomicCoordinates[k] = a[k][i];
                    }
                }
                if (radius[i] > radmax) {
                    radmax = radius[i];
                }
            }

            // Calculate width of cube from maximum atom radius and probe radius.
            width = 2.0 * (radmax + probe);

            // Set up cube arrays; first the integer coordinate arrays.
            for (int i = 0; i < nAtoms; i++) {
                for (int k = 0; k < 3; k++) {
                    cubeCoordinates[k][i] = (int) ((a[k][i] - minAtomicCoordinates[k]) / width);
                    if (cubeCoordinates[k][i] < 0) {
                        //logger.severe("Cube Coordinate Too Small");
                        throw new EnergyException("Cube Coordinate Too Small", false);
                    } else if (cubeCoordinates[k][i] > MAXCUBE) {
                        //logger.severe("Cube Coordinate Too Large");
                        throw new EnergyException("Cube Coordinate Too Large", false);
                    }
                }
            }

            // Initialize head pointer and srn=2 arrays.
//            for (int i = 0; i < MAXCUBE; i++) {
//                for (int j = 0; j < MAXCUBE; j++) {
//                    for (int k = 0; k < MAXCUBE; k++) {
//                        firstAtomPointer[i][j][k] = -1;
//                        activeCube[i][j][k] = false;
//                        activeAdjacentCube[i][j][k] = false;
//                    }
//                }
//            }
            fill(firstAtomPointer, -1);
            fill(activeCube, false);
            fill(activeAdjacentCube, false);

            // Initialize linked list pointers.
            fill(nextAtomPointer, -1);

            // Set up head and later pointers for each atom.
            outerloop:
            for (iatom = 0; iatom < nAtoms; iatom++) {

                // Skip atoms with surface request numbers of zero.
                if (skip[iatom]) {
                    continue;
                }
                getVector(ai, a, iatom);
                int i = cubeCoordinates[0][iatom];
                int j = cubeCoordinates[1][iatom];
                int k = cubeCoordinates[2][iatom];
                if (firstAtomPointer[index(i, j, k)] <= -1) {
                    // First atom in this cube.
                    firstAtomPointer[index(i, j, k)] = iatom;
                } else {
                    int counter = 1;
                    // Add to end of linked list.
                    iptr = firstAtomPointer[index(i, j, k)];
                    innerloop:
                    for (i = 0; i < counter; i++) {
                        getVector(aj, a, iptr);

                        // Check for duplicate atoms, turn off one of them.
                        if (dist2(ai, aj) <= 0.0) {
                            skip[iatom] = true;
                            continue outerloop;
                        }

                        // Move on down the list.
                        if (nextAtomPointer[iptr] <= -1.0) {
                            continue;
                        }
                        iptr = nextAtomPointer[iptr];
                        counter++;
                    }

                    // Store atom number.
                    nextAtomPointer[iptr] = iatom;
                }

                // Check for surfaced atom.
                if (!skip[iatom]) {
                    activeCube[index(i, j, k)] = true;
                }
            }

            // Check if this cube or any adjacent cube has active atoms.
            for (int k = 0; k < MAXCUBE; k++) {
                for (int j = 0; j < MAXCUBE; j++) {
                    for (int i = 0; i < MAXCUBE; i++) {
                        if (firstAtomPointer[index(i, j, k)] != -1) {
                            for (k1 = Math.max(k - 1, 0); k1 < Math.min(k + 2, MAXCUBE); k1++) {
                                for (j1 = Math.max(j - 1, 0); j1 < Math.min(j + 2, MAXCUBE); j1++) {
                                    for (i1 = Math.max(i - 1, 0); i1 < Math.min(i + 2, MAXCUBE); i1++) {
                                        if (activeCube[index(i1, j1, k1)]) {
                                            activeAdjacentCube[index(i, j, k)] = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ncls = -1;

            // Zero pointers for atom and find its cube.
            for (int i = 0; i < nAtoms; i++) {
                nclsa = -1;
                nosurf[i] = skip[i];
                acls[0][i] = -1;
                acls[1][i] = -1;
                if (skip[i]) {
                    continue;
                }
                ici = cubeCoordinates[0][i];
                icj = cubeCoordinates[1][i];
                ick = cubeCoordinates[2][i];

                // Skip iatom if its cube and adjoining cubes contain only blockers.
                if (!activeAdjacentCube[index(ici, icj, ick)]) {
                    continue;
                }
                sumi = 2.0 * probe + radius[i];

                // Check iatom cube and adjacent cubes for neighboring atoms.
                for (jck = max(ick - 1, 0); jck < min(ick + 2, MAXCUBE); jck++) {
                    for (jcj = max(icj - 1, 0); jcj < min(icj + 2, MAXCUBE); jcj++) {
                        for4:
                        for (jci = max(ici - 1, 0); jci < min(ici + 2, MAXCUBE); jci++) {
                            int j = firstAtomPointer[index(jci, jcj, jck)];
                            int counter = 1;
                            for (int q = 0; q < counter; q++) {
                                for (int z = 0; z < 1; z++) {
                                    // Check for end of linked list for this cube.
                                    if (j <= -1) {
                                        continue for4;
                                    } else if (i == j) {
                                        continue;
                                    } else if (skip[j]) {
                                        continue;
                                    }

                                    // Distance check.
                                    sum = sumi + radius[j];
                                    vect1 = abs(a[0][j] - a[0][i]);
                                    if (vect1 >= sum) {
                                        continue;
                                    }
                                    vect2 = abs(a[1][j] - a[1][i]);
                                    if (vect2 >= sum) {
                                        continue;
                                    }
                                    vect3 = abs(a[2][j] - a[2][i]);
                                    if (vect3 >= sum) {
                                        continue;
                                    }
                                    d2 = (vect1 * vect1) + (vect2 * vect2) + (vect3 * vect3);
                                    if (d2 >= sum * sum) {
                                        continue;
                                    }

                                    // Atoms are neighbors, save atom number in temporary array.
                                    if (!skip[j]) {
                                        nosurf[i] = false;
                                    }
                                    nclsa++;
                                    if (nclsa > maxclsa) {
                                        //logger.severe("Too many Neighbors for Atom");
                                        throw new EnergyException("Too many Neighbors for Atom", false);
                                    }
                                    tempNeighborList[nclsa] = j;
                                }

                                // Get number of next atom in cube.
                                j = nextAtomPointer[j];
                                counter++;
                            }
                        }
                    }
                }
                if (nosurf[i]) {
                    continue;
                }

                // Set up neighbors arrays with jatom in increasing order.
                jmold = -1;
                for (juse = 0; juse < nclsa + 1; juse++) {
                    jmin = nAtoms;
                    for (jcls = 0; jcls < nclsa + 1; jcls++) {

                        // Don't use ones already sorted.
                        if (tempNeighborList[jcls] > jmold) {
                            if (tempNeighborList[jcls] < jmin) {
                                jmin = tempNeighborList[jcls];
                                jmincls = jcls;
                            }
                        }
                    }
                    jmold = jmin;
                    jcls = jmincls;
                    jatom = tempNeighborList[jcls];
                    clsa[juse] = jatom;
                }

                // Set up pointers to first and last neighbors of atom.
                if (nclsa > -1) {
                    acls[0][i] = ncls + 1;
                    for (int m = 0; m < nclsa + 1; m++) {
                        ncls++;
                        if (ncls > maxcls) {
                            //logger.severe("Too many Neighboring Atom Pairs");
                            throw new EnergyException("Too many Neighboring Atom Pairs", false);
                        }
                        cls[ncls] = clsa[m];
                    }
                    acls[1][i] = ncls;
                }

            }
        }

        /**
         * The torus method sets a list of all of the temporary torus
         * positions by testing for a torus between each atom and its
         * neighbors
         */
        void torus() {
            int ia, ja, jn;
            int ibeg = -1;
            int iend = -1;
            double[] tt = new double[3];
            double[] ttax = new double[3];

            // No torus is possible if there is only one atom.
            ntt = -1;
            for (ia = 0; ia < nAtoms; ia++) {
                afree[ia] = true;
            }
            if (nAtoms < 1) {
                return;
            }

            // Get beginning and end pointers to neighbors of this atom.
            for (ia = 0; ia < nAtoms; ia++) {
                if (!nosurf[ia]) {
                    ibeg = acls[0][ia];
                    iend = acls[1][ia];
                }

                // Check for no neighbors.
                if (ibeg > -1) {
                    for (jn = ibeg; jn < iend + 1; jn++) {

                        // Clear pointer from neighbor to torus.
                        clst[jn] = -1;

                        // Get atom number of neighbor.
                        ja = cls[jn];

                        // Don't create torus twice.
                        if (ja >= ia) {

                            // Do some solid geometry.
                            double[] ttr = {0.0};

                            ttok = gettor(ia, ja, tt, ttr, ttax);
                            if (ttok) {
                                // We have temporary torus; set up variables.
                                ntt++;
                                if (ntt > maxtt) {
                                    //logger.severe("Too many Temporary Tori");
                                    throw new EnergyException("Too many Temporary Tori", false);
                                }

                                // Mark both atoms not free.
                                afree[ia] = false;
                                afree[ja] = false;
                                tta[0][ntt] = ia;
                                tta[1][ntt] = ja;

                                // Pointer from neighbor to torus.
                                clst[jn] = ntt;

                                // Initialize torus as both free and buried.
                                ttfree[ntt] = true;
                                ttbur[ntt] = true;

                                // Clear pointers from torus to first and last concave edges.
                                ttfe[ntt] = -1;
                                ttle[ntt] = -1;
                            }
                        }
                    }
                }
            }
        }

        /**
         * The place method finds the probe sites by putting the probe
         * sphere tangent to each triple of neighboring atoms.
         */
        public void place() {
            int[] mnb = new int[MAXMNB];
            int[] ikt = new int[MAXMNB];
            int[] jkt = new int[MAXMNB];
            int[] lkcls = new int[MAXMNB];
            double[] tik = new double[3];
            double[] tij = new double[3];
            double[] uij = new double[3];
            double[] uik = new double[3];
            double[] uijk = new double[3];
            double[] bij = new double[3];
            double[] bijk = new double[3];
            double[] aijk = new double[3];
            double[] pijk = new double[3];
            double[] tijik = new double[3];
            double[] tempv = new double[3];
            double[] utb = new double[3];
            double[] ai = new double[3];
            double[] ak = new double[3];
            double[] discls = new double[MAXMNB];
            double[] sumcls = new double[MAXMNB];
            boolean tb, tok, prbok;

            // No possible placement if there are no temporary tori.
            if (ntt <= -1) {
                return;
            }

            np = -1;
            nfn = -1;
            nen = -1;
            nv = -1;

            // Consider each torus in turn.
            for (int itt = 0; itt < ntt + 1; itt++) {
                // Get atom numbers.
                int ia = tta[0][itt];
                int ja = tta[1][itt];

                // Form mutual neighbor list; clear number of mutual neighbors of atoms ia and ja.
                int nmnb = -1;

                // Get beginning and end pointers for each atom's neighbor list.
                int iptr = acls[0][ia];
                int jptr = acls[0][ja];

                // For loops like this help replace go to statements from the Tinker code.
                outerloop:
                for (int i = 0; i < 1; i++) {

                    if (iptr <= -1 || jptr <= -1) {
                        continue;
                    }
                    int iend = acls[1][ia];
                    int jend = acls[1][ja];

                    // Collect mutual neighbors.
                    int counter = 1;
                    int counter2 = 1;
                    for (int t = 0; t < counter2; t++) {
                        for (int q = 0; q < counter; q++) {

                            // Check for end of loop.
                            if (iptr > iend || jptr > jend) {
                                continue outerloop;
                            }

                            // Go move the lagging pointer.
                            if (cls[iptr] < cls[jptr]) {
                                iptr++;
                                counter++;
                                continue;
                            }
                            if (cls[jptr] < cls[iptr]) {
                                jptr++;
                                counter++;
                            }
                        }
                            /*
                              Both point at same neighbor; one more mutual
                              neighbor save atom number of mutual neighbor.
                             */
                        nmnb++;
                        if (nmnb > MAXMNB) {
                            //logger.severe("Too many Mutual Neighbors");
                            throw new EnergyException("Too many Mutual Neighbors", false);
                        }
                        mnb[nmnb] = cls[iptr];

                        // Save pointers to second and third tori.
                        ikt[nmnb] = clst[iptr];
                        jkt[nmnb] = clst[jptr];
                        iptr++;
                        counter2++;
                    }
                }
                // We have all the mutual neighbors of ia and ja if no mutual neighbors, skip to end of loop.
                if (nmnb <= -1) {
                    ttbur[itt] = false;
                    continue;
                }
                double[] hij = {0.0};
                ttok = gettor(ia, ja, bij, hij, uij);
                for (int km = 0; km < nmnb + 1; km++) {
                    int ka = mnb[km];
                    getVector(ak, a, ka);
                    discls[km] = dist2(bij, ak);
                    sumcls[km] = (probe + radius[ka]) * (probe + radius[ka]);
                    // Initialize link to next farthest out neighbor.
                    lkcls[km] = -1;
                }

                // Set up a linked list of neighbors in order of increasing distance from ia-ja torus center.
                int lkf = 0;
                for (int z = 0; z < 1; z++) {
                    if (nmnb <= 0) {
                        continue;
                    }
                    // Put remaining neighbors in linked list at proper position.
                    int l;
                    for (l = 1; l < nmnb + 1; l++) {

                        int l1 = -1;
                        int l2 = lkf;
                        int counter2 = 1;
                        for (int w = 0; w < counter2; w++) {
                            if (discls[l] < discls[l2]) {
                                continue;
                            }
                            l1 = l2;
                            l2 = lkcls[l2];
                            if (l2 != -1) {
                                counter2++;
                            }
                        }

                        // Add to list.
                        if (l1 == -1) {
                            lkf = l;
                            lkcls[l] = l2;
                        } else {
                            lkcls[l1] = l;
                            lkcls[l] = l2;
                        }
                    }
                }
                // Loop through mutual neighbors.
                for (int km = 0; km < nmnb + 1; km++) {

                    // Get atom number of neighbors.
                    int ka = mnb[km];
                    if (skip[ia] && skip[ja] && skip[ka]) {
                        continue;
                    }

                    // Get tori numbers for neighbor.
                    int ik = ikt[km];
                    int jk = jkt[km];

                    // Possible new triple, do some geometry to retrieve saddle center, axis and radius.
                    prbok = false;
                    tb = false;
                    double[] rij = {0.0};
                    double hijk = 0.0;
                    tok = gettor(ia, ja, tij, rij, uij);
                    for (int w = 0; w < 1; w++) {
                        if (tok) {
                            getVector(ai, a, ka);
                            double dat2 = dist2(ai, tij);
                            double rad2 = (radius[ka] + probe) * (radius[ka] + probe) - rij[0] * rij[0];

                            // If "ka" less than "ja", then all we care about is whether the torus is buried.
                            if (ka < ja) {
                                if (rad2 <= 0.0 || dat2 > rad2) {
                                    continue;
                                }
                            }

                            double[] rik = {0.0};
                            tok = gettor(ia, ka, tik, rik, uik);
                            if (!tok) {
                                continue;
                            }
                            double dotijk = dot(uij, uik);
                            dotijk = check(dotijk);
                            double wijk = acos(dotijk);
                            double swijk = sin(wijk);

                            // If the three atoms are colinear, then there is no probe placement; but we still care
                            // whether the torus is buried by atom "k".
                            if (swijk == 0.0) {
                                tb = (rad2 > 0.0 && dat2 <= rad2);
                                continue;
                            }
                            cross(uij, uik, uijk);
                            for (int k = 0; k < 3; k++) {
                                uijk[k] = uijk[k] / swijk;
                            }
                            cross(uijk, uij, utb);
                            for (int k = 0; k < 3; k++) {
                                tijik[k] = tik[k] - tij[k];
                            }
                            double dotut = dot(uik, tijik);
                            double fact = dotut / swijk;
                            for (int k = 0; k < 3; k++) {
                                bijk[k] = tij[k] + utb[k] * fact;
                            }
                            getVector(ai, a, ia);
                            double dba = dist2(ai, bijk);
                            double rip2 = (radius[ia] + probe) * (radius[ia] + probe);
                            double rad = rip2 - dba;
                            if (rad < 0.0) {
                                tb = (rad2 > 0.0 && dat2 <= rad2);
                            } else {
                                prbok = true;
                                hijk = sqrt(rad);
                            }
                        }
                    }
                    if (tb) {
                        ttbur[itt] = true;
                        ttfree[itt] = false;
                        continue;
                    }

                    // Check for duplicate triples or any possible probe positions.
                    if (ka < ja || !prbok) {
                        continue;
                    }

                    // Altitude vector.
                    for (int k = 0; k < 3; k++) {
                        aijk[k] = hijk * uijk[k];
                    }

                    // We try two probe placements.
                    for3:
                    for (int ip = 0; ip < 2; ip++) {
                        for (int k = 0; k < 3; k++) {
                            if (ip == 0) {
                                pijk[k] = bijk[k] + aijk[k];
                            } else {
                                pijk[k] = bijk[k] - aijk[k];
                            }
                        }

                        // Mark three tori not free.
                        ttfree[itt] = false;
                        ttfree[ik] = false;
                        ttfree[jk] = false;

                        // Check for collisions.
                        int lm = lkf;
                        int counter3 = 1;
                        for4:
                        for (int t = 0; t < counter3; t++) {
                            for (int p = 0; p < 1; p++) {
                                if (lm < 0) {
                                    continue for4;
                                }

                                // Get atom number of mutual neighbor.
                                int la = mnb[lm];
                                // Must not equal third atom.
                                if (la == ka) {
                                    continue;
                                }
                                getVector(ak, a, la);
                                // Compare distance to sum of radii.
                                if (dist2(pijk, ak) <= sumcls[lm]) {
                                    continue for3;
                                }
                            }
                            lm = lkcls[lm];
                            counter3++;
                        }
                        // We have a new probe position.
                        np++;
                        if (np > maxp) {
                            //logger.severe("Too many Probe Positions");
                            throw new EnergyException("Too many Probe Positions", false);
                        }
                        // Mark three tori not buried.
                        ttbur[itt] = false;
                        ttbur[ik] = false;
                        ttbur[jk] = false;
                        // Store probe center.
                        for (int k = 0; k < 3; k++) {
                            p[k][np] = pijk[k];
                        }

                        // Calculate vectors from probe to atom centers.
                        if (nv + 3 > maxv) {
                            //logger.severe("Too many Vertices");
                            throw new EnergyException("Too many Vertices", false);
                        }
                        for (int k = 0; k < 3; k++) {
                            v[k][nv + 1] = a[k][ia] - p[k][np];
                            v[k][nv + 2] = a[k][ja] - p[k][np];
                            v[k][nv + 3] = a[k][ka] - p[k][np];
                        }
                        //double matrix[] = new double[9];
                        //int a = 0;
                        //for (int b = 0; b < 3; b++) {
                        //    for (int c = 0; c < 3; c++) {
                        //        matrix[a++] = v[b][nv + c + 1];
                        //    }
                        //}

                        // Calculate determinant of vectors defining triangle.
                        double det = v[0][nv + 1] * v[1][nv + 2] * v[2][nv + 3]
                                + v[0][nv + 2] * v[1][nv + 3] * v[2][nv + 1]
                                + v[0][nv + 3] * v[1][nv + 1] * v[2][nv + 2]
                                - v[0][nv + 3] * v[1][nv + 2] * v[2][nv + 1]
                                - v[0][nv + 2] * v[1][nv + 1] * v[2][nv + 3]
                                - v[0][nv + 1] * v[1][nv + 3] * v[2][nv + 2];
                        // Now add probe coordinates to vertices.
                        for (int k = 0; k < 3; k++) {
                            v[k][nv + 1] = p[k][np] + (v[k][nv + 1] * probe / (radius[ia] + probe));
                            v[k][nv + 2] = p[k][np] + (v[k][nv + 2] * probe / (radius[ja] + probe));
                            v[k][nv + 3] = p[k][np] + (v[k][nv + 3] * probe / (radius[ka] + probe));
                        }
                        // Want the concave face to have counter-clockwise orientation.
                        if (det > 0.0) {
                            // Swap second and third vertices.
                            for (int k = 0; k < 3; k++) {
                                tempv[k] = v[k][nv + 2];
                                v[k][nv + 2] = v[k][nv + 3];
                                v[k][nv + 3] = tempv[k];
                            }
                            // Set up pointers from probe to atoms.
                            pa[0][np] = ia;
                            pa[1][np] = ka;
                            pa[2][np] = ja;
                            // Set pointers from vertices to atoms.
                            va[nv + 1] = ia;
                            va[nv + 2] = ka;
                            va[nv + 3] = ja;
                            // Insert concave edges into linked lists for appropriate tori.
                            inedge(nen + 1, ik);
                            inedge(nen + 2, jk);
                            inedge(nen + 3, itt);
                        } else {
                            // Similarly, if face already counter-clockwise.
                            pa[0][np] = ia;
                            pa[1][np] = ja;
                            pa[2][np] = ka;
                            va[nv + 1] = ia;
                            va[nv + 2] = ja;
                            va[nv + 3] = ka;
                            inedge(nen + 1, itt);
                            inedge(nen + 2, jk);
                            inedge(nen + 3, ik);
                        }
                        // Set up pointers from vertices to probe.
                        for (int kv = 1; kv < 4; kv++) {
                            vp[nv + kv] = np;
                        }
                        // Set up concave edges and concave face.
                        if (nen + 3 > maxen) {
                            //logger.severe("Too many Concave Edges");
                            throw new EnergyException("Too many Concave Edges", false);
                        }
                        // Edges point to vertices.
                        env[0][nen + 1] = nv + 1;
                        env[1][nen + 1] = nv + 2;
                        env[0][nen + 2] = nv + 2;
                        env[1][nen + 2] = nv + 3;
                        env[0][nen + 3] = nv + 3;
                        env[1][nen + 3] = nv + 1;
                        if (nfn + 1 > maxfn) {
                            //logger.severe("Too many Concave Faces");
                            throw new EnergyException("Too many Concave Faces", false);
                        }
                        // Face points to edges.
                        for (int ke = 0; ke < 3; ke++) {
                            fnen[ke][nfn + 1] = nen + ke + 1;
                        }
                        // Increment counters for number of faces, edges, and vertices.
                        nfn++;
                        nen += 3;
                        nv += 3;
                    }
                }
            }
        }

        /**
         * The inedge method inserts a concave edge into the linked list for
         * its temporary torus.
         */
        void inedge(int edgeNumber, int torusNumber) {
            // Check for a serious error in the calling arguments.
            if (edgeNumber <= -1) {
                //logger.severe("Bad Edge Number in INEDGE");
                throw new EnergyException("Bad Edge Number in INEDGE", false);
            }
            if (torusNumber <= -1) {
                //logger.severe("Bad Torus Number in INEDGE");
                throw new EnergyException("Bad Torus Number in INEDGE", false);
            }
            // Set beginning of list or add to end.
            if (ttfe[torusNumber] == -1) {
                ttfe[torusNumber] = edgeNumber;
                enext[edgeNumber] = -1;
                ttle[torusNumber] = edgeNumber;
            } else {
                enext[ttle[torusNumber]] = edgeNumber;
                enext[edgeNumber] = -1;
                ttle[torusNumber] = edgeNumber;
            }
        }

        /**
         * The compress method transfers only the non-buried tori from the
         * temporary tori arrays to the final tori arrays.
         */
        void compress() {
            double[] ai = new double[3];
            double[] aj = new double[3];
            // Initialize the number of non-buried tori.
            nt = -1;
            if (ntt <= -1) {
                return;
            }
            // If torus is free, then it is not buried; skip to end of loop if buried torus.
            double[] trtemp = {0};
            for (int itt = 0; itt < ntt + 1; itt++) {
                if (ttfree[itt]) {
                    ttbur[itt] = false;
                }
                if (!ttbur[itt]) {
                    // First, transfer information.
                    nt++;
                    if (nt > maxt) {
                        throw new EnergyException("Too many non-buried tori.", false);
                    }
                    int ia = tta[0][itt];
                    int ja = tta[1][itt];
                    getVector(ai, t, nt);
                    getVector(aj, tax, nt);
                    ttok = gettor(ia, ja, ai, trtemp, aj);
                    t[0][nt] = ai[0];
                    t[1][nt] = ai[1];
                    t[2][nt] = ai[2];
                    tax[0][nt] = aj[0];
                    tax[1][nt] = aj[1];
                    tax[2][nt] = aj[2];
                    tr[nt] = trtemp[0];
                    ta[0][nt] = ia;
                    ta[1][nt] = ja;
                    tfree[nt] = ttfree[itt];
                    tfe[nt] = ttfe[itt];
                    // Special check for inconsistent probes.
                    int iptr = tfe[nt];
                    int ned = -1;
                    while (iptr != -1) {
                        ned++;
                        iptr = enext[iptr];
                    }
                    if ((ned % 2) == 0) {
                        iptr = tfe[nt];
                        while (iptr != -1) {
                            int iv1 = env[0][iptr];
                            int iv2 = env[1][iptr];
                            int ip1 = vp[iv1];
                            int ip2 = vp[iv2];
                            logger.warning(format("Odd Torus for Probes IP1 %d and IP2 %d", ip1, ip2));
                            iptr = enext[iptr];
                        }
                    }
                }
            }
        }

        /**
         * The saddles method constructs circles, convex edges, and saddle faces.
         */
        void saddles() {
            final int maxent = 500;
            int k, ia, in, ip;
            int it, iv, itwo;
            int ien, ient, nent;
            int l1, l2, m1, n1;
            int[] ten = new int[maxent];
            int[] nxtang = new int[maxent];
            double triple, factor;
            double dtev, dt;
            double[] ai = new double[3];
            double[] aj = new double[3];
            double[] ak = new double[3];
            double[] atvect = new double[3];
            double[] teang = new double[maxent];
            double[][] tev = new double[3][maxent];
            boolean[] sdstrt = new boolean[maxent];

            // Zero the number of circles, convex edges, and saddle faces.
            nc = -1;
            nep = -1;
            nfs = -1;
            fill(afe, -1);
            fill(ale, -1);
            fill(abur, true);

            // No saddle faces if no tori.
            if (nt < 0) {
                return;
            }

            // Cycle through tori.
            for (it = 0; it <= nt; it++) {
                if (skip[ta[0][it]] && skip[ta[1][it]]) {
                    continue;
                }

                // Set up two circles.
                for (in = 0; in < 2; in++) {
                    ia = ta[in][it];

                    // Mark atom not buried.
                    abur[ia] = false;

                    // Vector from atom to torus center.
                    for (k = 0; k < 3; k++) {
                        atvect[k] = t[k][it] - a[k][ia];
                    }
                    factor = radius[ia] / (radius[ia] + probe);

                    // One more circle.
                    nc++;
                    if (nc > maxc) {
                        throw new EnergyException("Too many Circles", false);
                    }

                    // Circle center.
                    for (k = 0; k < 3; k++) {
                        c[k][nc] = a[k][ia] + factor * atvect[k];
                    }

                    // Pointer from circle to atom.
                    ca[nc] = ia;

                    // Pointer from circle to torus.
                    ct[nc] = it;

                    // Circle radius.
                    cr[nc] = factor * tr[it];
                }

                // Skip to special "free torus" code.
                if (!tfree[it]) {
                    /*
                      Now we collect all the concave edges for this
                      torus; for each concave edge, calculate vector
                      from torus center through probe center and the
                      angle relative to first such vector.
                     */

                    // Clear the number of concave edges for torus.
                    nent = -1;

                    // Pointer to start of linked list.
                    ien = tfe[it];
                    // 10 Continue
                    while (ien >= 0) {
                        // One more concave edge.
                        nent++;
                        if (nent > maxent) {
                            //logger.severe("Too many Edges for Torus");
                            throw new EnergyException("Too many Edges for Torus", false);
                        }
                        // First vertex of edge.
                        iv = env[0][ien];

                        // Probe number of vertex.
                        ip = vp[iv];
                        for (k = 0; k < 3; k++) {
                            tev[k][nent] = p[k][ip] - t[k][it];
                        }
                        dtev = 0.0;
                        for (k = 0; k < 3; k++) {
                            dtev += tev[k][nent] * tev[k][nent];
                        }
                        if (dtev <= 0.0) {
                            throw new EnergyException("Probe on Torus Axis", false);
                        }
                        dtev = sqrt(dtev);
                        for (k = 0; k < 3; k++) {
                            tev[k][nent] = tev[k][nent] / dtev;
                        }

                        // Store concave edge number.
                        ten[nent] = ien;
                        if (nent > 0) {
                            // Calculate angle between this vector and first vector.
                            dt = 0.0;
                            for (k = 0; k < 3; k++) {
                                dt += tev[k][0] * tev[k][nent];
                            }
                            if (dt > 1.0) {
                                dt = 1.0;
                            }
                            if (dt < -1.0) {
                                dt = -1.0;
                            }

                            // Store angle.
                            teang[nent] = Math.acos(dt);

                            ai[0] = tev[0][0];
                            ai[1] = tev[1][0];
                            ai[2] = tev[2][0];
                            aj[0] = tev[0][nent];
                            aj[1] = tev[1][nent];
                            aj[2] = tev[2][nent];
                            ak[0] = tax[0][it];
                            ak[1] = tax[1][it];
                            ak[2] = tax[2][it];

                            // Get the sign right.
                            triple = triple(ai, aj, ak);
                            if (triple < 0.0) {
                                teang[nent] = 2.0 * Math.PI - teang[nent];
                            }
                        } else {
                            teang[0] = 0.0;
                        }
                        // Saddle face starts with this edge if it points parallel to torus axis vector
                        // (which goes from first to second atom).
                        sdstrt[nent] = (va[iv] == ta[0][it]);
                        // Next edge in list.
                        ien = enext[ien];
                    }

                    if (nent <= -1) {
                        logger.severe("No Edges for Non-free Torus");
                    }

                    itwo = 2;
                    if ((nent % itwo) == 0) {
                        throw new EnergyException("Odd Number of Edges for Torus", false);
                    }

                    // Set up linked list of concave edges in order of increasing angle
                    // around the torus axis;
                    // clear second linked (angle-ordered) list pointers.
                    for (ient = 0; ient <= nent; ient++) {
                        nxtang[ient] = -1;
                    }
                    for (ient = 1; ient <= nent; ient++) {
                        // We have an entry to put into linked list search for place to put it.
                        l1 = -1;
                        l2 = 0;
                        while (l2 != -1 && teang[ient] >= teang[l2]) {
                            l1 = l2;
                            l2 = nxtang[l2];
                        }
                        // We are at end of linked list or between l1 and l2; insert edge.
                        if (l1 <= -1) {
                            throw new EnergyException("Logic Error in SADDLES", true);
                        }
                        nxtang[l1] = ient;
                        nxtang[ient] = l2;
                    }

                    // Collect pairs of concave edges into saddles.
                    // Create convex edges while you're at it.
                    l1 = 0;
                    while (l1 > -1) {
                        // Check for start of saddle.
                        if (sdstrt[l1]) {
                            // One more saddle face.
                            nfs++;
                            if (nfs > maxfs) {
                                throw new EnergyException("Too many Saddle Faces", false);
                            }
                            // Get edge number.
                            ien = ten[l1];
                            // First concave edge of saddle.
                            fsen[0][nfs] = ien;
                            // One more convex edge.
                            nep++;
                            if (nep > maxep) {
                                throw new EnergyException("Too many Convex Edges", false);
                            }
                            // First convex edge points to second circle.
                            epc[nep] = nc;
                            // Atom circle lies on.
                            ia = ca[nc];
                            // Insert convex edge into linked list for atom.
                            ipedge(nep, ia);
                            // First vertex of convex edge is second vertex of concave edge.
                            epv[0][nep] = env[1][ien];
                            // First convex edge of saddle.
                            fsep[0][nfs] = nep;
                            // One more convex edge.
                            nep++;
                            if (nep > maxep) {
                                throw new EnergyException("Too many Convex Edges", false);
                            }
                            // Second convex edge points to fist circle.
                            epc[nep] = nc - 1;
                            ia = ca[nc - 1];
                            // Insert convex edge into linked list for atom.
                            ipedge(nep, ia);
                            // Second vertex of second convex edge is first vertex of first concave edge.
                            epv[1][nep] = env[0][ien];
                            l1 = nxtang[l1];
                            // Wrap around.
                            if (l1 <= -1) {
                                l1 = 0;
                            }
                            if (sdstrt[l1]) {
                                m1 = nxtang[l1];
                                if (m1 <= -1) {
                                    m1 = 0;
                                }
                                if (sdstrt[m1]) {
                                    throw new EnergyException("Three Starts in a Row", false);
                                }
                                n1 = nxtang[m1];
                                // The old switcheroo.
                                nxtang[l1] = n1;
                                nxtang[m1] = l1;
                                l1 = m1;
                            }
                            ien = ten[l1];
                            // Second concave edge for saddle face.
                            fsen[1][nfs] = ien;
                            // Second vertex of first convex edge is first vertex of second concave edge.
                            epv[1][nep - 1] = env[0][ien];
                            // First vertex of second convex edge is second vertex of second concave edge.
                            epv[0][nep] = env[1][ien];
                            fsep[1][nfs] = nep;
                            // Quit if we have wrapped around to first edge.
                            if (l1 == 0) {
                                break;
                            }
                        }
                        // Next concave edge.
                        l1 = nxtang[l1];
                    }
                } else {
                    // Free torus
                    // Set up entire circles as convex edges for new saddle surface; one more saddle face.
                    nfs++;
                    if (nfs > maxfs) {
                        throw new EnergyException("Too many Saddle Faces", false);
                    }
                    // No concave edge for saddle.
                    fsen[0][nfs] = -1;
                    fsen[1][nfs] = -1;
                    // One more convex edge.
                    nep++;
                    ia = ca[nc];
                    // Insert convex edge into linked list of atom.
                    ipedge(nep, ia);
                    // No vertices for convex edge.
                    epv[0][nep] = -1;
                    epv[1][nep] = -1;
                    // Pointer from convex edge to second circle.
                    epc[nep] = nc;
                    // First convex edge for saddle face.
                    fsep[0][nfs] = nep;
                    // One more convex edge.
                    nep++;
                    ia = ca[nc - 1];
                    // Insert second convex edge into linked list.
                    ipedge(nep, ia);
                    // No vertices for convex edge.
                    epv[0][nep] = -1;
                    epv[1][nep] = -1;
                    // Convex edge points to first circle.
                    epc[nep] = nc - 1;
                    // Second convex edge for saddle face.
                    fsep[1][nfs] = nep;
                    // Buried torus; do nothing with it.
                }
            }
        }

        /**
         * The triple method finds the triple product of three vectors; used
         * as a service routine by the Connolly surface area and volume
         * computation.
         */
        public double triple(double[] x, double[] y, double[] z) {
            double triple;
            double[] xy = new double[3];
            cross(x, y, xy);
            triple = dot(xy, z);
            return triple;
        }

        /**
         * The ipedge method inserts a convex edge into linked list for
         * atom.
         *
         * @param edgeNumber
         * @param atomNumber
         */
        void ipedge(int edgeNumber, int atomNumber) {
            // First, check for an error condition.
            if (edgeNumber <= -1) {
                //logger.severe("Bad Edge Number in IPEDGE");
                throw new EnergyException("Bad Edge Number in IPEDGE", true);
            }
            if (atomNumber <= -1) {
                //logger.severe("Bad Atom Number in IPEDGE");
                throw new EnergyException("Bad Atom Number in IPEDGE", true);
            }
            // Set beginning of list or add to end.
            if (afe[atomNumber] == -1) {
                afe[atomNumber] = edgeNumber;
                epnext[edgeNumber] = -1;
                ale[atomNumber] = edgeNumber;
            } else {
                epnext[ale[atomNumber]] = edgeNumber;
                epnext[edgeNumber] = -1;
                ale[atomNumber] = edgeNumber;
            }
        }

        /**
         * The contact method constructs the contact surface, cycles and
         * convex faces.
         */
        public void contact() {
            final int maxepa = 300;
            final int maxcypa = 100;
            int jepa;
            int[] aic = new int[maxepa];
            int[] aia = new int[maxepa];
            int[] aep = new int[maxepa];
            int[][] av = new int[2][maxepa];
            int[] ncyepa = new int[maxcypa];
            int[][] cyepa = new int[MAXCYEP][maxcypa];
            double[] ai = new double[3];
            double[][] acvect = new double[3][maxepa];
            double[][] aavect = new double[3][maxepa];
            double[] pole = new double[3];
            double[] acr = new double[maxepa];
            double[] unvect = new double[3];
            boolean[] epused = new boolean[maxepa];
            boolean[][] cycy = new boolean[maxcypa][maxcypa];
            boolean[] cyused = new boolean[maxcypa];
            boolean[][] samef = new boolean[maxcypa][maxcypa];

            // Zero out the number of cycles and convex faces.
            ncy = -1;
            nfp = -1;

            // Mark all free atoms not buried.
            for (int ia = 0; ia < nAtoms; ia++) {
                if (afree[ia]) {
                    abur[ia] = false;
                }
            }

            // Go through all atoms.
            firstloop:
            for (int ia = 0; ia < nAtoms; ia++) {
                if (skip[ia] || abur[ia]) {
                    continue;
                }

                // Special code for completely solvent-accessible atom.
                for (int o = 0; o < 1; o++) {
                    if (afree[ia]) {
                        continue;
                    }
                    // Gather convex edges for atom Clear number of convex edges for atom.
                    int nepa = -1;
                    // Pointer to first edge.
                    int iep = afe[ia];
                    int counter = 1;
                    for (int j = 0; j < counter; j++) {
                        // Check whether finished gathering.
                        if (iep <= -1) {
                            continue;
                        }
                        // One more edge.
                        nepa++;
                        if (nepa > maxepa) {
                            //logger.severe("Too many Convex Edges for Atom");
                            throw new EnergyException("Too many Convex Edges for Atom", false);
                        }
                        // Store vertices of edge.
                        av[0][nepa] = epv[0][iep];
                        av[1][nepa] = epv[1][iep];
                        // Store convex edge number.
                        aep[nepa] = iep;
                        int ic = epc[iep];
                        // Store circle number.
                        aic[nepa] = ic;
                        // Get neighboring atom.
                        int it = ct[ic];
                        int ia2;
                        if (ta[0][it] == ia) {
                            ia2 = ta[1][it];
                        } else {
                            ia2 = ta[0][it];
                        }
                        // Store other atom number, we might need it sometime.
                        aia[nepa] = ia2;
                            /*
                              Vector from atom to circle center; also vector
                              from atom to center of neighboring atom sometimes
                              we use one vector, sometimes the other.
                             */
                        for (int k = 0; k < 3; k++) {
                            acvect[k][nepa] = c[k][ic] - a[k][ia];
                            aavect[k][nepa] = a[k][ia2] - a[k][ia];
                        }
                        // Circle radius.
                        acr[nepa] = cr[ic];
                        // Pointer to next edge.
                        iep = epnext[iep];
                        counter++;
                    }
                    if (nepa <= -1) {
                        //logger.severe("No Edges for Non-buried, Non-free Atom");
                        throw new EnergyException("No Edges for Non-buried, Non-free Atom", false);
                    }
                    // Form cycles; initialize all the convex edges as not used in cycle.
                    for (int iepa = 0; iepa < nepa + 1; iepa++) {
                        epused[iepa] = false;
                    }
                    // Save old number of cycles.
                    int ncyold = ncy;
                    int nused = -1;
                    int ncypa = -1;
                    // Look for starting edge.
                    int iepa = 0;
                    int counter3 = 1;
                    outerloop:
                    for (int z = 0; z < counter3; z++) {
                        innerloop:
                        for (int w = 0; w < 1; w++) {
                            for (iepa = 0; iepa < nepa + 1; iepa++) {
                                if (!epused[iepa]) {
                                    continue innerloop;
                                }
                            }
                            continue outerloop;
                        }
                        // Cannot find starting edge; finished.
                        // Pointer to edge.
                        iep = aep[iepa];
                        // One edge so far on this cycle.
                        int ncyep = 0;
                        // One more cycle for atom.
                        ncypa++;
                        if (ncypa > maxcypa) {
                            //logger.severe("Too many Cycles per Atom");
                            throw new EnergyException("Too many Cycles per Atom", false);
                        }
                        // Mark edge used in cycle.
                        epused[iepa] = true;
                        nused++;
                        // One more cycle for molecule.
                        ncy++;
                        if (ncy > maxcy) {
                            //logger.severe("Too many Cycles");
                            throw new EnergyException("Too many Cycles", false);
                        }
                        // Index of edge in atom cycle array.
                        cyepa[ncyep][ncypa] = iepa;
                        // Store in molecule cycle array a pointer to edge.
                        cyep[ncyep][ncy] = iep;
                        // Second vertex of this edge is the vertex to look for
                        // next as the first vertex of another edge.
                        int lookv = av[1][iepa];
                        // If no vertex, this cycle is finished.
                        for (int t = 0; t < 1; t++) {
                            if (lookv <= -1) {
                                continue;
                            }
                            int counter2 = 1;

                            for (int r = 0; r < counter2; r++) {
                                // Look for next connected edge.
                                outer:
                                for (jepa = 0; jepa < nepa + 1; jepa++) {
                                    for (int q = 0; q < 1; q++) {
                                        if (epused[jepa]) {
                                            continue;
                                        }
                                        // Check second vertex of iepa versus first vertex of jepa.
                                        if (av[0][jepa] != lookv) {
                                            continue;
                                        }
                                        // Edges are connected pointer to edge.
                                        iep = aep[jepa];
                                        // One more edge for this cycle.
                                        ncyep++;
                                        if (ncyep > MAXCYEP) {
                                            //logger.severe("Too many Edges per Cycle");
                                            throw new EnergyException("Too many Edges per Cycle", false);
                                        }
                                        epused[jepa] = true;
                                        nused++;
                                        // Store index in local edge array.
                                        cyepa[ncyep][ncypa] = jepa;
                                        // Store pointer to edge.
                                        cyep[ncyep][ncy] = iep;
                                        // New vertex to look for.
                                        lookv = av[1][jepa];
                                        // If no vertex, this cycle is in trouble.
                                        if (lookv <= -1) {
                                            //logger.severe("Pointer Error in Cycle");
                                            throw new EnergyException("Pointer Error in Cycle", true);
                                        }
                                        counter2++;
                                        break outer;
                                    }
                                }
                            }
                            // It better connect to first edge of cycle.
                            if (lookv != av[0][iepa]) {
                                throw new EnergyException("Cycle does not Close", true);
                            }
                        }
                        // This cycle is finished, store number of edges in cycle.
                        ncyepa[ncypa] = ncyep;
                        cynep[ncy] = ncyep;
                        // Look for more cycles.
                        if (nused < nepa) {
                            counter3++;
                        }
                    }
                    // Compare cycles for inside/outside relation; check to
                    // see if cycle i is inside cycle j.
                    for (int icya = 0; icya < ncypa + 1; icya++) {
                        innerloop:
                        for (int jcya = 0; jcya < ncypa + 1; jcya++) {
                            int jcy = ncyold + jcya + 1;
                            // Initialize.
                            cycy[icya][jcya] = true;
                            // Check for same cycle.
                            if (icya == jcya || ncyepa[jcya] < 2) {
                                continue;
                            }
                            // If cycles i and j have a pair of edges belonging to the same circle,
                            // then they are outside each other.
                            for (int icyep = 0; icyep < ncyepa[icya] + 1; icyep++) {
                                iepa = cyepa[icyep][icya];
                                int ic = aic[iepa];
                                //for (int jcyep = 0; jcyep < ncyepa[jcya]; jcyep++) {
                                for (int jcyep = 0; jcyep < ncyepa[jcya] + 1; jcyep++) {
                                    jepa = cyepa[jcyep][jcya];
                                    int jc = aic[jepa];
                                    if (ic == jc) {
                                        cycy[icya][jcya] = false;
                                        continue innerloop;
                                    }
                                }
                            }
                            iepa = cyepa[0][icya];
                            ai[0] = aavect[0][iepa];
                            ai[1] = aavect[1][iepa];
                            ai[2] = aavect[2][iepa];
                            double anaa = r(ai);
                            double factor = radius[ia] / anaa;
                            // North pole and unit vector pointer south.
                            for (int k = 0; k < 3; k++) {
                                pole[k] = factor * aavect[k][iepa] + a[k][ia];
                                unvect[k] = -aavect[k][iepa] / anaa;
                            }
                            cycy[icya][jcya] = ptincy(pole, unvect, jcy);
                        }
                    }
                    // Group cycles into faces; direct comparison for i and j.
                    for (int icya = 0; icya < ncypa + 1; icya++) {
                        for (int jcya = 0; jcya < ncypa + 1; jcya++) {
                            // Tentatively say that cycles i and j bound the
                            // same face if they are inside each other.
                            samef[icya][jcya] = (cycy[icya][jcya] && cycy[jcya][icya]);
                        }
                    }
                    // If i is in exterior of k, and k is in interior of i
                    // and j, then i and j do not bound the same face.
                    for (int icya = 0; icya < ncypa + 1; icya++) {
                        for (int jcya = 0; jcya < ncypa + 1; jcya++) {
                            if (icya != jcya) {
                                for (int kcya = 0; kcya < ncypa + 1; kcya++) {
                                    if (kcya != icya && kcya != jcya) {
                                        if (cycy[kcya][icya] && cycy[kcya][jcya] && !cycy[icya][kcya]) {
                                            samef[icya][jcya] = false;
                                            samef[jcya][icya] = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // Fill gaps so that "samef" falls into complete blocks.
                    for (int icya = 0; icya < ncypa - 1; icya++) {
                        for (int jcya = icya + 1; jcya < ncypa; jcya++) {
                            if (samef[icya][jcya]) {
                                for (int kcya = jcya + 1; kcya < ncypa + 1; kcya++) {
                                    if (samef[jcya][kcya]) {
                                        samef[icya][kcya] = true;
                                        samef[kcya][icya] = true;
                                    }
                                }
                            }
                        }
                    }
                    // Group cycles belonging to the same face.
                    for (int icya = 0; icya < ncypa + 1; icya++) {
                        cyused[icya] = false;
                    }
                    // Clear number of cycles used in bounding faces.
                    nused = -1;
                    for (int icya = 0; icya < ncypa + 1; icya++) {
                        // Check for already used.
                        if (cyused[icya]) {
                            continue;
                        }
                        // One more convex face.
                        nfp++;
                        if (nfp > maxfp) {
                            //logger.severe("Too many Convex Faces");
                            throw new EnergyException("Too many Convex Faces", false);
                        }
                        // Clear number of cycles for face.
                        fpncy[nfp] = -1;
                        // Pointer from face to atom.
                        fpa[nfp] = ia;
                        // Look for all other cycles belonging to same face.
                        for (int jcya = 0; jcya < ncypa + 1; jcya++) {
                            // Check for cycle alraDY used in another face.
                            if (cyused[jcya] || !samef[icya][jcya]) {
                                continue;
                            }
                            // Mark cycle used.
                            cyused[jcya] = true;
                            nused++;
                            // One more cycle for face.
                            fpncy[nfp]++;
                            if (fpncy[nfp] > MAXFPCY) {
                                //logger.severe("Too many Cycles bounding Convex Face");
                                throw new EnergyException("Too many Cycles bounding Convex Face", false);
                            }
                            int i = fpncy[nfp];
                            // Store cycle number.
                            fpcy[i][nfp] = ncyold + jcya + 1;
                            // Check for finished.
                            if (nused >= ncypa) {
                                continue firstloop;
                            }
                        }
                    }
                    // Should not fall though end of for loops.
                    //logger.severe("Not all Cycles grouped into Convex Faces");
                    throw new EnergyException("Not all Cycles grouped into Convex Faces", true);
                }
                // Once face for free atom; no cycles.
                nfp++;
                if (nfp > maxfp) {
                    //logger.severe("Too many Convex Faces");
                    throw new EnergyException("Too many Convex Faces", false);
                }
                fpa[nfp] = ia;
                fpncy[nfp] = -1;
            }
        }

        boolean ptincy(double[] pnt, double[] unvect, int icy) {
            double[] acvect = new double[3];
            double[] cpvect = new double[3];
            double[] polev = new double[3];
            double[][] spv = new double[3][MAXCYEP];
            double[][] epu = new double[3][MAXCYEP];
            // Check for being eaten by neighbor.
            int iatom = 0;
            for (int ke = 0; ke < cynep[icy]; ke++) {
                int iep = cyep[ke][icy];
                int ic = epc[iep];
                int it = ct[ic];
                iatom = ca[ic];
                int iaoth;
                if (ta[0][it] == iatom) {
                    iaoth = ta[1][it];
                } else {
                    iaoth = ta[0][it];
                }
                for (int k = 0; k < 3; k++) {
                    acvect[k] = a[k][iaoth] - a[k][iatom];
                    cpvect[k] = pnt[k] - c[k][ic];
                }
                if (dot(acvect, cpvect) >= 0.0) {
                    return false;
                }
            }
            if (cynep[icy] <= 1) {
                return true;
            }
            int nedge = cynep[icy];
            for (int ke = 0; ke < cynep[icy]; ke++) {

                // Vertex number (use first vertex of edge).
                int iep = cyep[ke][icy];
                iv = epv[0][iep];
                if (iv != 0) {
                    // Vector from north pole to vertex.
                    for (int k = 0; k < 3; k++) {
                        polev[k] = v[k][iv] - pnt[k];
                    }
                    // Calculate multiplication factor.
                    double dt = dot(polev, unvect);
                    if (dt == 0.0) {
                        return true;
                    }
                    double f = (radius[iatom] * radius[iatom]) / dt;
                    if (f < 1.0) {
                        return true;
                    }
                    // Projected vertex for this convex edge.
                    for (int k = 0; k < 3; k++) {
                        spv[k][ke] = pnt[k] + f * polev[k];
                    }
                }
            }
            epuclc(spv, nedge, epu);
            double totang = rotang(epu, nedge, unvect);
            return (totang > 0.0);
        }

        double rotang(double[][] epu, int nedge, double[] unvect) {
            double[] crs = new double[3];
            double[] ai = new double[3];
            double[] aj = new double[3];
            double[] ak = new double[3];

            double totang = 0.0;

            // Sum angles at vertices of cycle.
            for (int ke = 0; ke < nedge + 1; ke++) {
                double dt;
                if (ke < nedge) {
                    ai[0] = epu[0][ke];
                    ai[1] = epu[1][ke];
                    ai[2] = epu[2][ke];
                    aj[0] = epu[0][ke + 1];
                    aj[1] = epu[1][ke + 1];
                    aj[2] = epu[2][ke + 1];
                    dt = dot(ai, aj);
                    cross(ai, ai, crs);
                } else {
                    ak[0] = epu[0][0];
                    ak[1] = epu[1][0];
                    ak[2] = epu[2][0];
                    // Closing edge of cycle.
                    dt = dot(ai, ak);
                    cross(ai, ak, crs);
                }
                if (dt < -1.0) {
                    dt = -1.0;
                } else if (dt > 1.0) {
                    dt = 1.0;
                }
                double ang = acos(dt);
                if (dot(crs, unvect) > 0.0) {
                    ang = -ang;
                }
                // Add to total for cycle.
                totang += ang;
            }
            return totang;
        }

        void epuclc(double[][] spv, int nedge, double[][] epu) {
            double[] ai = new double[3];
            // Calculate unit vectors along edges.
            for (int ke = 0; ke < nedge; ke++) {
                // Get index of second edge of corner.
                int ke2;
                if (ke < nedge) {
                    ke2 = ke + 1;
                } else {
                    ke2 = 0;
                }
                // Unit vector along edge of cycle.
                for (int k = 0; k < 3; k++) {
                    epu[k][ke] = spv[k][ke2] - spv[k][ke];
                }
                getVector(ai, epu, ke);
                double epun = r(ai);
                if (epun <= 0.0) {
                    //logger.severe("Null Edge in Cycle");
                    throw new EnergyException("Null Edge in Cycle", true);
                }
                // Normalize.
                if (epun > 0.0) {
                    for (int k = 0; k < 3; k++) {
                        epu[k][ke] = epu[k][ke] / epun;
                    }
                } else {
                    for (int k = 0; k < 3; k++) {
                        epu[k][ke] = 0.0;
                    }
                }
            }
            // Vectors for null edges come from following or preceding edges.
            for (int ke = 0; ke < nedge; ke++) {
                getVector(ai, epu, ke);
                if (r(ai) <= 0.0) {
                    int le = ke - 1;
                    if (le <= 0) {
                        le = nedge;
                    }
                    for (int k = 0; k < 3; k++) {
                        epu[k][ke] = epu[k][le];
                    }
                }
            }
        }

        /**
         * The measpm method computes the volume of a single prism section
         * of the full interior polyhedron.
         */
        double measpm(int ifn) {
            double[][] pav = new double[3][3];
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] vect3 = new double[3];
            double height = 0.0;
            for (int ke = 0; ke < 3; ke++) {
                int ien = fnen[ke][ifn];
                iv = env[0][ien];
                int ia = va[iv];
                height += a[2][ia];
                int ip = vp[iv];
                for (int k = 0; k < 3; k++) {
                    pav[k][ke] = a[k][ia] - p[k][ip];
                }
            }
            height *= 1 / 3.0;
            for (int k = 0; k < 3; k++) {
                vect1[k] = pav[k][1] - pav[k][0];
                vect2[k] = pav[k][2] - pav[k][0];
            }
            cross(vect1, vect2, vect3);
            return height * vect3[2] / 2.0;
        }

        void measfp(int ifp, double[] av) {
            double angle;
            double[] ai = new double[3];
            double[] aj = new double[3];
            double[] ak = new double[3];
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] acvect = new double[3];
            double[] aavect = new double[3];
            double[][][] tanv = new double[3][2][MAXCYEP];
            double[][] radial = new double[3][MAXCYEP];
            double pcurve = 0.0;
            double gcurve = 0.0;
            int ia = fpa[ifp];
            int ncycle = fpncy[ifp];
            int ieuler;
            if (ncycle > -1) {
                ieuler = 1 - ncycle;
            } else {
                ieuler = 2;
            }
            for (int icyptr = 0; icyptr < ncycle + 1; icyptr++) {
                int icy = fpcy[icyptr][ifp];
                int nedge = cynep[icy];
                for (int ke = 0; ke < nedge + 1; ke++) {
                    int iep = cyep[ke][icy];
                    int ic = epc[iep];
                    int it = ct[ic];
                    int ia2;
                    if (ia == ta[0][it]) {
                        ia2 = ta[1][it];
                    } else {
                        ia2 = ta[0][it];
                    }
                    for (int k = 0; k < 3; k++) {
                        acvect[k] = c[k][ic] - a[k][ia];
                        aavect[k] = a[k][ia2] - a[k][ia];
                    }
                    norm(aavect, aavect);
                    double dt = dot(acvect, aavect);
                    double geo = -dt / (radius[ia] * cr[ic]);
                    int iv1 = epv[0][iep];
                    int iv2 = epv[1][iep];
                    if (iv1 == -1 || iv2 == -1) {
                        angle = 2.0 * PI;
                    } else {
                        for (int k = 0; k < 3; k++) {
                            vect1[k] = v[k][iv1] - c[k][ic];
                            vect2[k] = v[k][iv2] - c[k][ic];
                            radial[k][ke] = v[k][iv1] - a[k][ia];
                        }
                        getVector(ai, radial, ke);
                        norm(ai, ai);
                        radial[0][ke] = ai[0];
                        radial[1][ke] = ai[1];
                        radial[2][ke] = ai[2];
                        getVector(ai, tanv, 0, ke);
                        cross(vect1, aavect, aj);
                        norm(aj, aj);
                        tanv[0][0][ke] = aj[0];
                        tanv[1][0][ke] = aj[1];
                        tanv[2][0][ke] = aj[2];
                        getVector(ak, tanv, 1, ke);
                        cross(vect2, aavect, ak);
                        norm(ak, ak);
                        tanv[0][1][ke] = ak[0];
                        tanv[1][1][ke] = ak[1];
                        tanv[2][1][ke] = ak[2];
                        angle = vecang(vect1, vect2, aavect, -1.0);
                    }
                    gcurve += cr[ic] * angle * geo;
                    if (nedge != 0) {
                        if (ke > 0) {
                            getVector(ai, tanv, 1, ke - 1);
                            getVector(aj, tanv, 0, ke);
                            getVector(ak, radial, ke);
                            angle = vecang(ai, aj, ak, 1.0);
                            if (angle < 0.0) {
                                //logger.severe("Negative Angle in MEASFP");
                                throw new EnergyException("Negative Angle in MEASFP", true);
                            }
                            pcurve += angle;
                        }
                    }
                }
                if (nedge > 0) {
                    getVector(ai, tanv, 1, nedge);
                    getVector(aj, tanv, 0, 0);
                    getVector(ak, radial, 0);
                    angle = vecang(ai, aj, ak, 1.0);
                    if (angle < 0.0) {
                        //logger.severe("Negative Angle in MEASFP");
                        throw new EnergyException("Negative Angle in MEASFP", true);
                    }
                    pcurve += angle;
                }
            }
            double gauss = 2.0 * PI * ieuler - pcurve - gcurve;
            double areap = gauss * (radius[ia] * radius[ia]);
            double volp = areap * radius[ia] / 3.0;
            av[0] = areap;
            av[1] = volp;
        }

        void measfs(int ifs, double[] saddle) {
            double areas = 0.0;
            double vols = 0.0;
            double areasp = 0.0;
            double volsp = 0.0;
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] aavect = new double[3];

            int iep = fsep[0][ifs];
            int ic = epc[iep];
            int it = ct[ic];
            int ia1 = ta[0][it];
            int ia2 = ta[1][it];
            for (int k = 0; k < 3; k++) {
                aavect[k] = a[k][ia2] - a[k][ia1];
            }
            norm(aavect, aavect);
            int iv1 = epv[0][iep];
            int iv2 = epv[1][iep];
            double phi;
            if (iv1 == -1 || iv2 == -1) {
                phi = 2.0 * PI;
            } else {
                for (int k = 0; k < 3; k++) {
                    vect1[k] = v[k][iv1] - c[k][ic];
                    vect2[k] = v[k][iv2] - c[k][ic];
                }
                phi = vecang(vect1, vect2, aavect, 1.0);
            }
            for (int k = 0; k < 3; k++) {
                vect1[k] = a[k][ia2] - t[k][it];
                vect2[k] = a[k][ia2] - t[k][it];
            }
            double d1 = -1.0 * dot(vect1, aavect);
            double d2 = dot(vect2, aavect);
            theta1 = atan2(d1, tr[it]);
            theta2 = atan2(d2, tr[it]);

            // Check for cusps.
            double thetaq;
            boolean cusp;
            if (tr[it] < probe && theta1 > 0.0 && theta2 > 0.0) {
                cusp = true;
                double rat = tr[it] / probe;
                if (rat > 1.0) {
                    rat = 1.0;
                } else if (rat < -1.0) {
                    rat = -1.0;
                }
                thetaq = acos(rat);
            } else {
                cusp = false;
                thetaq = 0.0;
                areasp = 0.0;
                volsp = 0.0;
            }
            double term1 = tr[it] * probe * (theta1 + theta2);
            double term2 = (probe * probe) * (sin(theta1) + sin(theta2));
            areas = phi * (term1 - term2);
            if (cusp) {
                double spin = tr[it] * probe * thetaq - probe * probe * sin(thetaq);
                areasp = 2.0 * phi * spin;
            }

            iep = fsep[0][ifs];
            int ic2 = epc[iep];
            iep = fsep[1][ifs];
            int ic1 = epc[iep];
            if (ca[ic1] != ia1) {
                throw new EnergyException("IA1 Inconsistency in MEASFS", true);
            }
            for (int k = 0; k < 3; k++) {
                vect1[k] = c[k][ic1] - a[k][ia1];
                vect2[k] = c[k][ic2] - a[k][ia2];
            }
            double w1 = dot(vect1, aavect);
            double w2 = -1.0 * dot(vect2, aavect);
            double cone1 = phi * ((w1 * cr[ic1] * cr[ic1])) / 6.0;
            double cone2 = phi * ((w2 * cr[ic2] * cr[ic2])) / 6.0;
            term1 = (tr[it] * tr[it]) * probe * (sin(theta1) + sin(theta2));
            term2 = sin(theta1) * cos(theta1) + theta1 + sin(theta2) * cos(theta2) + theta2;
            term2 = tr[it] * (probe * probe) * term2;
            double term3 = sin(theta1) * cos(theta1) * cos(theta1)
                    + 2.0 * sin(theta1) + sin(theta2) * cos(theta2) * cos(theta2)
                    + 2.0 * sin(theta2);
            term3 = (probe * probe * probe / 3.0) * term3;
            double volt = (phi / 2.0) * (term1 - term2 + term3);
            vols = volt + cone1 + cone2;
            if (cusp) {
                term1 = (tr[it] * tr[it]) * probe * sin(thetaq);
                term2 = sin(thetaq) * cos(thetaq) + thetaq;
                term2 = tr[it] * (probe * probe) * term2;
                term3 = sin(thetaq) * cos(thetaq) * cos(thetaq) + 2.0 * sin(thetaq);
                term3 = (probe * probe * probe / 3.0) * term3;
                volsp = phi * (term1 - term2 + term3);
            }
            saddle[0] = areas;
            saddle[1] = vols;
            saddle[2] = areasp;
            saddle[3] = volsp;
        }

        void measfn(int ifn, double[] av) {
            double[] ai = new double[3];
            double[] aj = new double[3];
            double[] ak = new double[3];
            double[] angle = new double[3];
            double[][] pvv = new double[3][3];
            double[][] pav = new double[3][3];
            double[][] planev = new double[3][3];
            double arean;
            double voln;
            for (int ke = 0; ke < 3; ke++) {
                int ien = fnen[ke][ifn];
                int iv = env[0][ien];
                int ia = va[iv];
                int ip = vp[iv];
                for (int k = 0; k < 3; k++) {
                    pvv[k][ke] = v[k][iv] - p[k][ip];
                    pav[k][ke] = a[k][ia] - p[k][ip];
                }
                if (probe > 0.0) {
                    getVector(ai, pvv, ke);
                    norm(ai, ai);
                }
            }
            if (probe <= 0.0) {
                arean = 0.0;
            } else {
                for (int ke = 0; ke < 3; ke++) {
                    int je = ke + 1;
                    if (je > 2) {
                        je = 0;
                    }
                    getVector(ai, pvv, ke);
                    getVector(aj, pvv, je);
                    getVector(ak, planev, ke);
                    cross(ai, aj, ak);
                    norm(ak, ak);
                }
                for (int ke = 0; ke < 3; ke++) {
                    int je = ke - 1;
                    if (je < 0) {
                        je = 2;
                    }
                    getVector(ai, planev, je);
                    getVector(aj, planev, ke);
                    getVector(ak, pvv, ke);
                    angle[ke] = vecang(ai, aj, ak, -1.0);
                    if (angle[ke] < 0.0) {
                        throw new EnergyException("Negative Angle in MEASFN", true);
                    }
                }
                double defect = 2.0 * PI - (angle[0] + angle[1] + angle[2]);
                arean = (probe * probe) * defect;
            }
            getVector(ai, pav, 0);
            getVector(aj, pav, 1);
            getVector(ak, pav, 2);
            double simplx = -triple(ai, aj, ak) / 6.0;
            voln = simplx - arean * probe / 3.0;
            av[0] = arean;
            av[1] = voln;

        }

        /**
         * The vam method takes the analytical molecular surface defined as
         * a collection of spherical and toroidal polygons and uses it to
         * compute the volume and surface area
         */
        void vam() {
            if (nfn < 0) {
                nfn = 0;
            }
            final int maxdot = 1000;
            final int maxop = 100;
            final int nscale = 20;
            int[] ivs = new int[3];
            int[] ispind = new int[3];
            int[] ispnd2 = new int[3];
            int[] ifnop = new int[maxop];
            int[] nlap = new int[nfn];
            int[] enfs = new int[20 * nAtoms];
            int[][] fnt = new int[3][nfn];
            int[][] nspt = new int[3][nfn];
            double[][] cenop = new double[3][maxop];
            double[] sdot = new double[3];
            double[] dotv = new double[nscale];
            double[] tau = new double[3];
            double[] ppm = new double[3];
            double[] xpnt1 = new double[3];
            double[] xpnt2 = new double[3];
            double[] qij = new double[3];
            double[] qji = new double[3];
            double[][] vects = new double[3][3];
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] vect3 = new double[3];
            double[] vect4 = new double[3];
            double[] vect5 = new double[3];
            double[] vect6 = new double[3];
            double[] vect7 = new double[3];
            double[] vect8 = new double[3];
            double[] upp = new double[3];
            double[] thetaq = new double[3];
            double[] sigmaq = new double[3];
            double[] umq = new double[3];
            double[] upq = new double[3];
            double[] uc = new double[3];
            double[] uq = new double[3];
            double[] uij = new double[3];
            double[] ai = new double[3];
            double[] aj = new double[3];
            double[] ak = new double[3];
            double[] al = new double[3];
            double[][] dots = new double[3][maxdot];
            double[][] tdots = new double[3][maxdot];
            double[] atmarea = new double[nAtoms];
            double[] depths = new double[nfn];
            double[] cora = new double[nfn];
            double[] corv = new double[nfn];
            double[][] alts = new double[3][nfn];
            double[][] fncen = new double[3][nfn];
            double[][][] fnvect = new double[3][3][nfn];
            boolean[] ate = new boolean[maxop];
            boolean[] badav = new boolean[nfn];
            boolean[] badt = new boolean[nfn];
            boolean[][] fcins = new boolean[3][nfn];
            boolean[][] fcint = new boolean[3][nfn];
            boolean[][] fntrev = new boolean[3][nfn];
            boolean[] vip = new boolean[3];
            boolean[] cintp = new boolean[1];
            boolean[] cinsp = new boolean[1];

            // Compute the volume of the interior polyhedron.
            double polyhedronVolume = 0.0;
            for (int ifn = 0; ifn < nfn + 1; ifn++) {
                polyhedronVolume += measpm(ifn);
            }

            // Compute the area and volume due to convex faces as well as
            // the area partitioned among the atoms.
            double totap = 0.0;
            double totvp = 0.0;
            fill(atmarea, 0.0);
            double[] convexFaces = {0.0, 0.0};
            for (int ifp = 0; ifp < nfp + 1; ifp++) {
                measfp(ifp, convexFaces);
                int ia = fpa[ifp];
                atmarea[ia] += convexFaces[0];
                totap += convexFaces[0];
                totvp += convexFaces[1];
            }

            // Compute the area and volume due to saddle faces
            // as well as the spindle correction value.
            double totas = 0.0;
            double totvs = 0.0;
            double totasp = 0.0;
            double totvsp = 0.0;
            double[] saddle = {0.0, 0.0, 0.0, 0.0};
            for (int ifs = 0; ifs < nfs + 1; ifs++) {
                for (int k = 0; k < 2; k++) {
                    int ien = fsen[k][ifs];
                    if (ien > -1) {
                        enfs[ien] = ifs;
                    }
                }
                measfs(ifs, saddle);
                double areas = saddle[0];
                double vols = saddle[1];
                double areasp = saddle[2];
                double volsp = saddle[3];
                totas += areas;
                totvs += vols;
                totasp += areasp;
                totvsp += volsp;
                if (areas - areasp < 0.0) {
                    //logger.severe("Negative Area for Saddle Face");
                    throw new EnergyException("Negative Area for Saddle Face", true);
                }
            }

            // Compute the area and volume due to concave faces.
            double totan = 0.0;
            double totvn = 0.0;
            double[] concaveFaces = {0.0, 0.0};
            for (int ifn = 0; ifn < nfn + 1; ifn++) {
                measfn(ifn, concaveFaces);
                double arean = concaveFaces[0];
                double voln = concaveFaces[1];
                totan += arean;
                totvn += voln;
            }

            // Compute the area and volume lens correction values.
            double alenst = 0.0;
            double alensn = 0.0;
            double vlenst = 0.0;
            double vlensn = 0.0;
            for (int j = 0; j < 1; j++) {
                if (probe <= 0.0) {
                    continue;
                }
                int[] ndots = {maxdot};
                gendot(ndots, dots, probe, 0.0, 0.0, 0.0);
                double dota = (4.0 * PI * probe * probe) / ndots[0];
                for (int ifn = 0; ifn < nfn + 1; ifn++) {
                    nlap[ifn] = -1;
                    cora[ifn] = 0.0;
                    corv[ifn] = 0.0;
                    badav[ifn] = false;
                    badt[ifn] = false;
                    for (int k = 0; k < 3; k++) {
                        nspt[k][ifn] = -1;
                    }
                    int ien = fnen[0][ifn];
                    iv = env[0][ien];
                    int ip = vp[iv];
                    getVector(ai, alts, ifn);
                    depths[ifn] = depth(ip, ai);
                    for (int k = 0; k < 3; k++) {
                        fncen[k][ifn] = p[k][ip];
                    }
                    // This assigned value was never used?
                    int ia = va[iv];
                    // Get vertices and vectors.
                    for (int ke = 0; ke < 3; ke++) {
                        ien = fnen[ke][ifn];
                        ivs[ke] = env[0][ien];
                        ia = va[ivs[ke]];
                        int ifs = enfs[ien];
                        int iep = fsep[0][ifs];
                        int ic = epc[iep];
                        int it = ct[ic];
                        fnt[ke][ifn] = it;
                        fntrev[ke][ifn] = (ta[0][it] != ia);
                    }
                    for (int ke = 0; ke < 3; ke++) {
                        for (int k = 0; k < 3; k++) {
                            vects[k][ke] = v[k][ivs[ke]] - p[k][ip];
                        }
                    }
                    // Calculate normal vectors for the three planes that
                    // cut out the geodesic triangle.
                    getVector(ai, vects, 0);
                    getVector(aj, vects, 1);
                    getVector(ak, fnvect, 0, ifn);
                    cross(ai, aj, ak);
                    norm(ak, ak);
                    getVector(ai, vects, 2);
                    getVector(ak, fnvect, 1, ifn);
                    cross(aj, ai, ak);
                    norm(ak, ak);
                    getVector(aj, vects, 0);
                    getVector(ak, fnvect, 2, ifn);
                    cross(ai, aj, ak);
                    norm(ak, ak);
                }
                for (int ifn = 0; ifn < nfn; ifn++) {
                    for (int jfn = ifn + 1; jfn < nfn + 1; jfn++) {
                        getVector(ai, fncen, ifn);
                        getVector(aj, fncen, jfn);
                        double dij2 = dist2(ai, aj);
                        if (dij2 > 4.0 * probe * probe) {
                            continue;
                        }
                        if (depths[ifn] > probe && depths[jfn] > probe) {
                            continue;
                        }
                        // These two probes may have intersecting surfaces.
                        double dpp = dist(ai, aj);
                        // Compute the midpoint.
                        for (int k = 0; k < 3; k++) {
                            ppm[k] = (fncen[k][ifn] + fncen[k][jfn]) / 2.0;
                            upp[k] = (fncen[k][jfn] - fncen[k][ifn]) / dpp;
                        }
                        double rm = probe * probe - (dpp / 2.0) * (dpp / 2.0);
                        if (rm < 0.0) {
                            rm = 0.0;
                        }
                        rm = sqrt(rm);
                        double rat = dpp / (2.0 * probe);
                        check(rat);
                        double rho = asin(rat);
                        // Use circle-place intersection routine.
                        boolean alli = true;
                        boolean anyi = false;
                        boolean spindl = false;
                        for (int k = 0; k < 3; k++) {
                            ispind[k] = -1;
                            ispnd2[k] = -1;
                        }
                        for (int ke = 0; ke < 3; ke++) {
                            thetaq[ke] = 0.0;
                            sigmaq[ke] = 0.0;
                            tau[ke] = 0.0;
                            getVector(ai, fncen, ifn);
                            getVector(aj, fnvect, ke, ifn);
                            cirpln(ppm, rm, upp, ai, aj, cintp, cinsp, xpnt1, xpnt2);
                            fcins[ke][ifn] = cinsp[0];
                            fcint[ke][ifn] = cintp[0];
                            if (!cinsp[0]) {
                                alli = false;
                            }
                            if (cintp[0]) {
                                anyi = true;
                            }
                            if (!cintp[0]) {
                                continue;
                            }
                            int it = fnt[ke][ifn];
                            if (tr[it] > probe) {
                                continue;
                            }
                            for (int ke2 = 0; ke2 < 3; ke2++) {
                                if (it == fnt[ke2][jfn]) {
                                    ispind[ke] = it;
                                    nspt[ke][ifn]++;
                                    ispnd2[ke2] = it;
                                    nspt[ke][jfn]++;
                                    spindl = true;
                                }
                            }
                            if (ispind[ke] == -1) {
                                continue;
                            }

                            // Check that the two ways of calculating intersection points match.
                            rat = tr[it] / probe;
                            check(rat);
                            thetaq[ke] = acos(rat);
                            double stq = sin(thetaq[ke]);
                            if (fntrev[ke][ifn]) {
                                for (int k = 0; k < 3; k++) {
                                    uij[k] = -tax[k][it];
                                }
                            } else {
                                for (int k = 0; k < 3; k++) {
                                    uij[k] = tax[k][it];
                                }
                            }
                            for (int k = 0; k < 3; k++) {
                                qij[k] = t[k][it] - stq * probe * uij[k];
                                qji[k] = t[k][it] + stq * probe * uij[k];
                            }
                            for (int k = 0; k < 3; k++) {
                                umq[k] = (qij[k] - ppm[k]) / rm;
                                upq[k] = (qij[k] - fncen[k][ifn]) / probe;
                            }
                            cross(uij, upp, vect1);
                            double dt = dot(umq, vect1);
                            check(dt);
                            sigmaq[ke] = acos(dt);
                            getVector(ai, fnvect, ke, ifn);
                            cross(upq, ai, vect1);
                            norm(vect1, uc);
                            cross(upp, upq, vect1);
                            norm(vect1, uq);
                            dt = dot(uc, uq);
                            check(dt);
                            tau[ke] = PI - acos(dt);
                        }
                        boolean allj = true;
                        boolean anyj = false;
                        for (int ke = 0; ke < 3; ke++) {
                            getVector(ai, fncen, jfn);
                            getVector(aj, fnvect, ke, jfn);
                            cirpln(ppm, rm, upp, ai, aj, cinsp, cintp, xpnt1, xpnt2);
                            fcins[ke][jfn] = cinsp[0];
                            fcint[ke][jfn] = cintp[0];
                            if (!cinsp[0]) {
                                allj = false;
                            }
                            if (cintp[0]) {
                                anyj = true;
                            }
                        }
                        boolean case1 = (alli && allj && !anyi && !anyj);
                        boolean case2 = (anyi && anyj && spindl);
                        if (!case1 && !case2) {
                            continue;
                        }
                        // This kind of overlap can be handled.
                        nlap[ifn]++;
                        nlap[jfn]++;
                        for (int ke = 0; ke < 3; ke++) {
                            int ien = fnen[ke][ifn];
                            int iv1 = env[0][ien];
                            int iv2 = env[1][ien];
                            for (int k = 0; k < 3; k++) {
                                vect3[k] = v[k][iv1] - fncen[k][ifn];
                                vect4[k] = v[k][iv2] - fncen[k][ifn];
                            }
                            for (int ke2 = 0; ke2 < 3; ke2++) {
                                if (ispind[ke] == ispnd2[ke2] || ispind[ke] == -1) {
                                    continue;
                                }
                                getVector(ai, fncen, ifn);
                                getVector(aj, fnvect, ke, ifn);
                                getVector(ak, fncen, jfn);
                                getVector(al, fnvect, ke2, jfn);
                                cirpln(ai, probe, aj, ak, al, cinsp, cintp, xpnt1, xpnt2);
                                if (!cintp[0]) {
                                    continue;
                                }
                                ien = fnen[ke2][jfn];
                                iv1 = env[0][ien];
                                iv2 = env[1][ien];
                                for (int k = 0; k < 3; k++) {
                                    vect7[k] = v[k][iv1] - fncen[k][jfn];
                                    vect8[k] = v[k][iv2] - fncen[k][jfn];
                                }
                                // Check whether point lies on spindle arc.
                                for (int k = 0; k < 3; k++) {
                                    vect1[k] = xpnt1[k] - fncen[k][ifn];
                                    vect2[k] = xpnt2[k] - fncen[k][ifn];
                                    vect5[k] = xpnt1[k] - fncen[k][jfn];
                                    vect6[k] = xpnt2[k] - fncen[k][jfn];
                                }
                                // Continue to next if statement if any of the following are true.
                                getVector(ai, fnvect, ke, ifn);
                                getVector(aj, fnvect, ke2, jfn);
                                endloop:
                                for (int u = 0; u < 1; u++) {
                                    outerloop:
                                    for (int z = 0; z < 1; z++) {
                                        for (int t = 0; t < 1; t++) {
                                            if (triple(vect3, vect1, ai) < 0.0
                                                    || triple(vect1, vect4, ai) < 0.0
                                                    || triple(vect7, vect5, aj) < 0.0
                                                    || triple(vect5, vect8, aj) < 0.0) {
                                                continue;
                                            }
                                            continue outerloop;
                                        }
                                        if (triple(vect3, vect2, ai) < 0.0
                                                || triple(vect2, vect4, ai) < 0.0
                                                || triple(vect7, vect6, aj) < 0.0
                                                || triple(vect6, vect8, aj) < 0.0) {
                                            continue endloop;
                                        }
                                    }
                                    badav[ifn] = true;
                                }
                            }
                        }
                        for (int ke = 0; ke < 3; ke++) {
                            int ien = fnen[ke][ifn];
                            int iv1 = env[0][ien];
                            int iv2 = env[1][ien];
                            for (int k = 0; k < 3; k++) {
                                vect3[k] = v[k][iv1] - fncen[k][ifn];
                                vect4[k] = v[k][iv2] - fncen[k][ifn];
                            }
                            for (int ke2 = 0; ke2 < 3; ke2++) {
                                if (ispind[ke] == ispnd2[ke2] || ispnd2[ke2] == -1) {
                                    continue;
                                }
                                getVector(ai, fncen, ifn);
                                getVector(aj, fnvect, ke, ifn);
                                getVector(ak, fncen, jfn);
                                getVector(al, fnvect, ke2, jfn);
                                cirpln(ak, probe, al, ai, aj, cinsp, cintp, xpnt1, xpnt2);
                                if (!cintp[0]) {
                                    continue;
                                }
                                ien = fnen[ke2][jfn];
                                iv1 = env[0][ien];
                                iv2 = env[1][ien];
                                for (int k = 0; k < 3; k++) {
                                    vect7[k] = v[k][iv1] - fncen[k][jfn];
                                    vect8[k] = v[k][iv2] - fncen[k][jfn];
                                }
                                // Check whether point lies on spindle arc.
                                for (int k = 0; k < 3; k++) {
                                    vect1[k] = xpnt1[k] - fncen[k][ifn];
                                    vect2[k] = xpnt2[k] - fncen[k][ifn];
                                    vect5[k] = xpnt1[k] - fncen[k][jfn];
                                    vect6[k] = xpnt2[k] - fncen[k][jfn];
                                }
                                // Continue to next if statement if any of the following are true.
                                getVector(ai, fnvect, ke, ifn);
                                getVector(aj, fnvect, ke2, jfn);
                                if (triple(vect3, vect1, ai) < 0.0
                                        || triple(vect1, vect4, ai) < 0.0
                                        || triple(vect7, vect5, aj) < 0.0
                                        || triple(vect5, vect8, aj) < 0.0) {
                                    if (!(triple(vect3, vect2, ai) < 0.0
                                            || triple(vect2, vect4, ai) < 0.0
                                            || triple(vect7, vect6, aj) < 0.0
                                            || triple(vect6, vect8, aj) < 0.0)) {
                                        badav[jfn] = true;
                                    }
                                } else {
                                    badav[jfn] = true;
                                }
                            }

                            double sumlam = 0.0;
                            double sumsig = 0.0;
                            double sumsc = 0.0;
                            for (int k = 0; k < 3; k++) {
                                if (ispind[ke] != 0) {
                                    sumlam += PI - tau[ke];
                                    sumsig += sigmaq[ke] - PI;
                                    sumsc += sin(sigmaq[ke]) * cos(sigmaq[ke]);
                                }
                            }
                            double alens = 2.0 * probe * probe
                                    * (PI - sumlam - sin(rho) * (PI + sumsig));
                            double vint = alens * probe / 3.0;
                            double vcone = probe * rm * rm * sin(rho) * (PI + sumsig) / 3.0;
                            double vpyr = probe * rm * rm * sin(rho) * sumsc / 3.0;
                            double vlens = vint - vcone + vpyr;
                            cora[ifn] += alens;
                            cora[jfn] += alens;
                            corv[ifn] += vlens;
                            corv[jfn] += vlens;
                        }
                        // Check for vertex on opposing probe in face.
                        outerloop:
                        for (int kv = 0; kv < 3; kv++) {
                            vip[kv] = false;
                            int ien = fnen[kv][jfn];
                            iv = env[0][ien];
                            for (int k = 0; k < 3; k++) {
                                vect1[k] = v[k][iv] - fncen[k][ifn];
                            }
                            norm(vect1, vect1);
                            for (int ke = 0; ke < 3; ke++) {
                                getVector(ai, fnvect, ke, ifn);
                                getVector(aj, v, iv);
                                double dt = dot(ai, aj);
                                if (dt > 0.0) {
                                    continue outerloop;
                                }
                                vip[kv] = true;
                            }
                        }
                    }
                }
                for (int ifn = 0; ifn < nfn + 1; ifn++) {
                    for (int ke = 0; ke < 3; ke++) {
                        if (nspt[ke][ifn] > 0) {
                            badt[ifn] = true;
                        }
                    }
                }
                for (int ifn = 0; ifn < nfn + 1; ifn++) {
                    if (nlap[ifn] <= -1) {
                        continue;
                    }
                    // Gather all overlapping probes.
                    int nop = -1;
                    for (int jfn = 0; jfn < nfn + 1; jfn++) {
                        if (ifn != jfn) {
                            getVector(ai, fncen, ifn);
                            getVector(aj, fncen, jfn);
                            double dij2 = dist2(ai, aj);
                            if (dij2 <= 4.0 * probe * probe) {
                                if (depths[jfn] <= probe) {
                                    nop++;
                                    if (nop > maxop) {
                                        //logger.severe("NOP Overflow in VAM");
                                        throw new EnergyException("NOP Overflow in VAM", false);
                                    }
                                    ifnop[nop] = jfn;
                                    for (int k = 0; k < 3; k++) {
                                        cenop[k][nop] = fncen[k][jfn];
                                    }
                                }
                            }
                        }
                        // Numerical calculation of the correction.
                        double areado = 0.0;
                        double voldo = 0.0;
                        double scinc = 1.0 / nscale;
                        for (int isc = 0; isc < nscale; isc++) {
                            double rsc = isc - 0.5;
                            dotv[isc] = probe * dota * rsc * rsc * scinc * scinc * scinc;
                        }
                        for (int iop = 0; iop < nop + 1; iop++) {
                            ate[iop] = false;
                        }
                        int neatmx = 0;
                        for (int idot = 0; idot < ndots[0]; idot++) {
                            boolean move = false;
                            for (int ke = 0; ke < 3; ke++) {
                                getVector(ai, fnvect, ke, ifn);
                                getVector(aj, dots, idot);
                                double dt = dot(ai, aj);
                                if (dt > 0.0) {
                                    move = true;
                                    break;
                                }
                            }
                            if (move) {
                                continue;
                            }
                            for (int k = 0; k < 3; k++) {
                                tdots[k][idot] = fncen[k][ifn] + dots[k][idot];
                            }
                            for (int iop = 0; iop < nop + 1; iop++) {
                                jfn = ifnop[iop];
                                getVector(ai, dots, idot);
                                getVector(aj, fncen, jfn);
                                double ds2 = dist2(ai, aj);
                                if (ds2 > probe * probe) {
                                    areado += dota;
                                    break;
                                }
                            }
                            for (int isc = 0; isc < nscale; isc++) {
                                double rsc = isc - 0.5;
                                for (int k = 0; k < 3; k++) {
                                    sdot[k] = fncen[k][ifn] + rsc * scinc * dots[k][idot];
                                }
                                int neat = 0;
                                for (int iop = 0; iop < nop + 1; iop++) {
                                    jfn = ifnop[iop];
                                    getVector(ai, fncen, jfn);
                                    double ds2 = dist2(sdot, ai);
                                    if (ds2 > probe * probe) {
                                        for (int k = 0; k < 3; k++) {
                                            vect1[k] = sdot[k] - fncen[k][jfn];
                                        }
                                        for (int ke = 0; ke < 3; ke++) {
                                            getVector(ai, fnvect, ke, jfn);
                                            double dt = dot(ai, vect1);
                                            if (dt > 0.0) {
                                                move = true;
                                                break;
                                            }
                                        }
                                        if (move) {
                                            break;
                                        }
                                        neat++;
                                        ate[iop] = true;
                                    }
                                }
                                if (neat > neatmx) {
                                    neatmx = neat;
                                }
                                if (neat > 0) {
                                    voldo += dotv[isc] * (neat / (1.0 + neat));
                                }
                            }
                        }
                        double coran = areado;
                        double corvn = voldo;
                        int nate = 0;
                        for (int iop = 0; iop < nop + 1; iop++) {
                            if (ate[iop]) {
                                nate++;
                            }
                        }
                        // Use either the analytical or numerical correction.
                        boolean usenum = (nate > nlap[ifn] || neatmx > 1 || badt[ifn]);
                        if (usenum) {
                            cora[ifn] = coran;
                            corv[ifn] = corvn;
                            alensn += cora[ifn];
                            vlensn += corv[ifn];
                        } else if (badav[ifn]) {
                            corv[ifn] = corvn;
                            vlensn += corv[ifn];
                        }
                        alenst += cora[ifn];
                        vlenst += corv[ifn];
                    }
                }
            }
            // Finally, compute the total area and total volume.
            double area = totap + totas + totan - totasp - alenst;
            double volume = totvp + totvs + totvn + polyhedronVolume - totvsp + vlenst;
            //logger.info(format(" Volume = %16.8f, Area = %16.8f", volume, area));
            //logger.info(format(" Total Volume        %16.8f", volume));
            //logger.info(format(" Total Area          %16.8f", area));

            localVolume += volume;
            localSurfaceArea += area;
        }

        /**
         * The gendot method finds the coordinates of a specified number of
         * surface points for a sphere with the input radius and coordinate
         * center.
         */
        void gendot(int[] ndots, double[][] dots, double radius,
                    double xcenter, double ycenter, double zcenter) {
            int nequat = (int) sqrt(PI * ((double) ndots[0]));
            int nvert = nequat / 2;
            if (nvert < 0) {
                nvert = 0;
            }
            int k = 0;
            outerloop:
            for (int i = -1; i < nvert + 1; i++) {
                double fi = (PI * ((double) i)) / ((double) nvert);
                double z = cos(fi);
                double xy = sin(fi);
                int nhoriz = (int) (nequat * xy);
                if (nhoriz < 0) {
                    nhoriz = 0;
                }
                for (int j = -1; j < nhoriz; j++) {
                    double fj = (2.0 * PI * ((double) (j))) / ((double) (nhoriz));
                    double x = cos(fj) * xy;
                    double y = sin(fj) * xy;
                    k++;
                    dots[0][k] = x * radius + xcenter;
                    dots[1][k] = y * radius + ycenter;
                    dots[2][k] = z * radius + zcenter;
                    if (k >= ndots[0]) {
                        break outerloop;
                    }
                }
            }
            ndots[0] = k;
        }

        /**
         * The cirpln method determines the points of intersection between a
         * specified circle and plane.
         */
        boolean cirpln(double[] circen, double cirrad, double[] cirvec,
                       double[] plncen, double[] plnvec, boolean[] cinsp,
                       boolean[] cintp, double[] xpnt1, double[] xpnt2) {
            double[] cpvect = new double[3];
            double[] pnt1 = new double[3];
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] uvect1 = new double[3];
            double[] uvect2 = new double[3];
            for (int k = 0; k < 3; k++) {
                cpvect[k] = plncen[k] - circen[k];
            }
            double dcp = dot(cpvect, plnvec);
            cinsp[0] = (dcp > 0.0);
            cross(plnvec, cirvec, vect1);
            if (r(vect1) > 0.0) {
                norm(vect1, uvect1);
                cross(cirvec, uvect1, vect2);
                if (r(vect2) > 0.0) {
                    norm(vect2, uvect2);
                    double dir = dot(uvect2, plnvec);
                    if (dir != 0.0) {
                        double ratio = dcp / dir;
                        if (abs(ratio) <= cirrad) {
                            for (int k = 0; k < 3; k++) {
                                pnt1[k] = circen[k] + ratio * uvect2[k];
                            }
                            double rlen = cirrad * cirrad - ratio * ratio;
                            if (rlen < 0.0) {
                                rlen = 0.0;
                            }
                            rlen = sqrt(rlen);
                            for (int k = 0; k < 3; k++) {
                                xpnt1[k] = pnt1[k] - rlen * uvect1[k];
                                xpnt2[k] = pnt1[k] + rlen * uvect1[k];
                            }
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        /**
         * The vecang method finds the angle between two vectors handed with
         * respect to a coordinate axis; returns an angle in the range
         * [0,2*PI].
         */
        double vecang(double[] v1, double[] v2, double[] axis, double hand) {
            double a1 = r(v1);
            double a2 = r(v2);
            double dt = dot(v1, v2);
            double a12 = a1 * a2;
            if (abs(a12) != 0.0) {
                dt = dt / a12;
            }
            dt = check(dt);
            double angle = acos(dt);
            double vecang;
            if (hand * triple(v1, v2, axis) < 0.0) {
                vecang = 2.0 * PI - angle;
            } else {
                vecang = angle;
            }
            return vecang;
        }

        public double depth(int ip, double[] alt) {
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] vect3 = new double[3];
            double[] vect4 = new double[3];
            int ia1 = pa[0][ip];
            int ia2 = pa[1][ip];
            int ia3 = pa[2][ip];
            for (int k = 0; k < 3; k++) {
                vect1[k] = a[k][ia1] - a[k][ia3];
                vect2[k] = a[k][ia2] - a[k][ia3];
                vect3[k] = p[k][ip] - a[k][ia3];
            }
            cross(vect1, vect2, vect4);
            norm(vect4, vect4);
            double dot = dot(vect4, vect3);
            for (int k = 0; k < 3; k++) {
                alt[k] = vect4[k];
            }
            return dot;
        }

        public double check(double angle) {
            if (angle > 1.0) {
                angle = 1.0;
            } else if (angle < -1.0) {
                angle = -1.0;
            }
            return angle;
        }

        void calcDerivative(int lb, int ub) {
                /*
                  Fix the step size in the z-direction; this value sets the
                  accuracy of the numerical derivatives; zstep=0.06 is a good
                  balance between compute time and accuracy.
                 */
            double zstep = 0.0601;

            // Load the cubes based on coarse lattice; first of all set edge
            // length to the maximum diameter of any atom.
            double edge = 2.0 * rmax;
            int nx = (int) ((xmax - xmin) / edge);
            int ny = (int) ((ymax - ymin) / edge);
            int nz = (int) ((zmax - zmin) / edge);
            if (max(max(nx, ny), nz) > mxcube) {
                //logger.severe(" VOLUME1  --  Increase the Value of MAXCUBE");
                throw new EnergyException(" VOLUME1  --  Increase the Value of MAXCUBE", false);
            }

            // Initialize the coarse lattice of cubes.
//            for (int i = 0; i <= nx; i++) {
//                for (int j = 0; j <= ny; j++) {
//                    for (int k = 0; k <= nz; k++) {
//                        cube[0][index2(i,j,k)] = -1;
//                        cube[1][index2(i,j,k)] = -1;
//                    }
//                }
//            }
            fill(cube[0], -1);
            fill(cube[1], -1);

            // Find the number of atoms in each cube.
            for (int m = 0; m < nAtoms; m++) {
                if (!skip[m]) {
                    int i = (int) ((x[m] - xmin) / edge);
                    int j = (int) ((y[m] - ymin) / edge);
                    int k = (int) ((z[m] - zmin) / edge);
                    cube[0][index2(i, j, k)]++;
                }
            }

                /*
                  Determine the highest index in the array "itab" for the atoms
                  that fall into each cube; the first cube that has atoms
                  defines the first index for "itab"; the final index for the
                  atoms in the present cube is the final index of the last cube
                  plus the number of atoms in the present cube.
                 */
            int isum = 0;
            int tcube;
            for (int i = 0; i <= nx; i++) {
                for (int j = 0; j <= ny; j++) {
                    for (int k = 0; k <= nz; k++) {
                        tcube = cube[0][index2(i, j, k)];
                        if (tcube != -1) {
                            isum += tcube;
                            cube[1][index2(i, j, k)] = isum;
                        }
                    }
                }
            }

            // "cube(1,,,)" now contains a pointer to the array "itab"
            // giving the position of the last entry for the list of atoms
            // in that cube of total number equal to "cube(0,,,)".
            for (int m = 0; m < nAtoms; m++) {
                if (!skip[m]) {
                    int i = (int) ((x[m] - xmin) / edge);
                    int j = (int) ((y[m] - ymin) / edge);
                    int k = (int) ((z[m] - zmin) / edge);
                    tcube = cube[1][index2(i, j, k)];
                    itab[tcube] = m;
                    cube[1][index2(i, j, k)]--;
                }
            }

            // Set "cube(1,,,)" to be the starting index in "itab" for atom
            // list of that cube; and "cube(0,,,)" to be the stop index.
            isum = 0;
            for (int i = 0; i <= nx; i++) {
                for (int j = 0; j <= ny; j++) {
                    for (int k = 0; k <= nz; k++) {
                        tcube = cube[0][index2(i, j, k)];
                        if (tcube != -1) {
                            isum += tcube;
                            cube[0][index2(i, j, k)] = isum;
                            cube[1][index2(i, j, k)]++;
                        }
                    }
                }
            }

            // Process in turn each atom from the coordinate list; first
            // select the potential intersecting atoms.
            int ir;
            for (ir = 0; ir < nAtoms; ir++) {
                double pre_dx = 0.0;
                double pre_dy = 0.0;
                double pre_dz = 0.0;
                if (skip[ir]) {
                    continue;
                }
                double rr = vdwrad[ir];
                double rrx2 = 2.0 * rr;
                double rrsq = rr * rr;
                double xr = x[ir];
                double yr = y[ir];
                double zr = z[ir];
                // Find cubes to search for overlaps for current atom.
                int istart = (int) ((xr - xmin) / edge);
                int istop = min(istart + 2, nx + 1);
                istart = max(istart, 1);
                int jstart = (int) ((yr - ymin) / edge);
                int jstop = min(jstart + 2, ny + 1);
                jstart = max(jstart, 1);
                int kstart = (int) ((zr - zmin) / edge);
                int kstop = min(kstart + 2, nz + 1);
                kstart = max(kstart, 1);
                // Load all overlapping atoms into "inov".
                int io = -1;
                int in;
                for (int i = istart - 1; i < istop; i++) {
                    for (int j = jstart - 1; j < jstop; j++) {
                        for (int k = kstart - 1; k < kstop; k++) {
                            int mstart = cube[1][index2(i, j, k)];
                            if (mstart != -1) {
                                int mstop = cube[0][index2(i, j, k)];
                                for (int m = mstart; m <= mstop; m++) {
                                    in = itab[m];
                                    if (in != ir) {
                                        io++;
                                        if (io > MAXARC) {
                                            //logger.severe(" VOLUME1  --  Increase the Value of MAXARC");
                                        }
                                        dx[io] = x[in] - xr;
                                        dy[io] = y[in] - yr;
                                        dsq[io] = (dx[io] * dx[io]) + (dy[io] * dy[io]);
                                        double dist2 = dsq[io] + ((z[in] - zr) * (z[in] - zr));
                                        double vdwsum = (rr + vdwrad[in]) * (rr + vdwrad[in]);
                                        if (dist2 > vdwsum || dist2 == 0.0) {
                                            io--;
                                        } else {
                                            d[io] = sqrt(dsq[io]);
                                            inov[io] = in;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Determine resolution along the z-axis.
                if (io != -1) {
                    double ztop = zr + rr;
                    double ztopshave = ztop - zstep;
                    double zgrid = zr - rr;
                    // Half of the part not covered by the planes.
                    zgrid += 0.5 * (rrx2 - (((int) (rrx2 / zstep)) * zstep));
                    double zstart = zgrid;
                    // Section atom spheres perpendicular to the z-axis.
                    while (zgrid <= ztop) {
                        // "rsecr" is radius of circle of intersection of "ir" sphere on the current sphere.
                        double rsec2r = rrsq - ((zgrid - zr) * (zgrid - zr));
                        if (rsec2r < 0.0) {
                            rsec2r = 0.000001;
                        }
                        double rsecr = sqrt(rsec2r);
                        double phi1;
                        double cos_phi1;
                        if (zgrid >= ztopshave) {
                            cos_phi1 = 1.0;
                            phi1 = 0.0;
                        } else {
                            cos_phi1 = (zgrid + (0.5 * zstep) - zr) / rr;
                            phi1 = acos(cos_phi1);
                        }
                        double phi2;
                        double cos_phi2;
                        if (zgrid == zstart) {
                            cos_phi2 = -1.0;
                            phi2 = PI;
                        } else {
                            cos_phi2 = (zgrid - (0.5 * zstep) - zr) / rr;
                            phi2 = acos(cos_phi2);
                        }
                        // Check intersection of neighbor circles.
                        int narc = -1;
                        for (int k = 0; k <= io; k++) {
                            in = inov[k];
                            double rinsq = vdwrad[in] * vdwrad[in];
                            double rsec2n = rinsq - ((zgrid - z[in]) * (zgrid - z[in]));
                            if (rsec2n > 0.0) {
                                double rsecn = sqrt(rsec2n);
                                if (d[k] < (rsecr + rsecn)) {
                                    double rdiff = rsecr - rsecn;
                                    if (d[k] <= abs(rdiff)) {
                                        if (rdiff < 0.0) {
                                            narc = 0;
                                            arci[narc] = 0.0;
                                            arcf[narc] = pix2;
                                        }
                                        continue;
                                    }
                                    narc++;
                                    if (narc > MAXARC) {
                                        //logger.info("VOLUME1 -- Increase the Value of MAXARC");
                                    }
                                        /*
                                          Initial and final arc endpoints are
                                          found for intersection of "ir" circle
                                          with another circle contained in same
                                          plane; the initial endpoint of the
                                          enclosed arc is stored in "arci", the
                                          final endpoint in "arcf"; get
                                          "cosine" via law of cosines.
                                         */
                                    double cosine = (dsq[k] + rsec2r - rsec2n) / (2.0 * d[k] * rsecr);
                                    cosine = min(1.0, max(-1.0, cosine));
                                        /*
                                          "alpha" is the angle between a line
                                          containing either point of
                                          intersection and the reference circle
                                          center and the line containing both
                                          circle centers; "beta" is the angle
                                          between the line containing both
                                          circle centers and x-axis.
                                         */
                                    double alpha = acos(cosine);
                                    double beta = atan2(dy[k], dx[k]);
                                    if (dy[k] < 0.0) {
                                        beta += pix2;
                                    }
                                    double ti = beta - alpha;
                                    double tf = beta + alpha;
                                    if (ti < 0.0) {
                                        ti += pix2;
                                    }
                                    if (tf > pix2) {
                                        tf -= pix2;
                                    }
                                    arci[narc] = ti;
                                    // If the arc crosses zero, then it is broken into two segments; the first
                                    // ends at two pi and the second begins at zero.
                                    if (tf < ti) {
                                        arcf[narc] = pix2;
                                        narc++;
                                        arci[narc] = 0.0;
                                    }
                                    arcf[narc] = tf;
                                }
                            }
                        }

                        // Find the pre-area and pre-forces on this section
                        // (band), "pre-" means a multiplicative factor is
                        // yet to be applied.
                        double seg_dz;
                        if (narc == -1) {
                            seg_dz = pix2 * ((cos_phi1 * cos_phi1) - (cos_phi2 * cos_phi2));
                            pre_dz += seg_dz;
                        } else {
                            // Sort the arc endpoint arrays, each with
                            // "narc" entries, in order of increasing values
                            // of the arguments in "arci".
                            double temp;
                            for (int k = 0; k < narc; k++) {
                                double aa = arci[k];
                                double bb = arcf[k];
                                temp = 1000000.0;
                                for (int i = k; i <= narc; i++) {
                                    if (arci[i] <= temp) {
                                        temp = arci[i];
                                        itemp = i;
                                    }
                                }
                                arci[k] = arci[itemp];
                                arcf[k] = arcf[itemp];
                                arci[itemp] = aa;
                                arcf[itemp] = bb;
                            }
                            // Consolidate arcs by removing overlapping arc endpoints.
                            temp = arcf[0];
                            int j = 0;
                            for (int k = 1; k <= narc; k++) {
                                if (temp < arci[k]) {
                                    arcf[j] = temp;
                                    j++;
                                    arci[j] = arci[k];
                                    temp = arcf[k];
                                } else if (temp < arcf[k]) {
                                    temp = arcf[k];
                                }
                            }
                            arcf[j] = temp;
                            narc = j;
                            if (narc == 0) {
                                narc = 1;
                                arcf[1] = pix2;
                                arci[1] = arcf[0];
                                arcf[0] = arci[0];
                                arci[0] = 0.0;
                            } else {
                                temp = arci[0];
                                for (int k = 0; k < narc; k++) {
                                    arci[k] = arcf[k];
                                    arcf[k] = arci[k + 1];
                                }

                                if (temp == 0.0 && arcf[narc] == pix2) {
                                    narc--;
                                } else {
                                    arci[narc] = arcf[narc];
                                    arcf[narc] = temp;
                                }
                            }

                            // Compute the numerical pre-derivative values.
                            for (int k = 0; k <= narc; k++) {
                                theta1 = arci[k];
                                theta2 = arcf[k];
                                double dtheta;
                                if (theta2 >= theta1) {
                                    dtheta = theta2 - theta1;
                                } else {
                                    dtheta = (theta2 + pix2) - theta1;
                                }
                                double phi_term = phi2 - phi1 - 0.5 * (sin(2.0 * phi2) - sin(2.0 * phi1));
                                double seg_dx = (sin(theta2) - sin(theta1)) * phi_term;
                                double seg_dy = (cos(theta1) - cos(theta2)) * phi_term;
                                seg_dz = dtheta * ((cos_phi1 * cos_phi1) - (cos_phi2 * cos_phi2));
                                pre_dx += seg_dx;
                                pre_dy += seg_dy;
                                pre_dz += seg_dz;
                            }
                        }
                        zgrid += zstep;
                    }
                }
                dex[0][ir] = 0.5 * rrsq * pre_dx;
                dex[1][ir] = 0.5 * rrsq * pre_dy;
                dex[2][ir] = 0.5 * rrsq * pre_dz;
            }
        }

        @Override
        public void run(int lb, int ub) {
            setRadius();
            calcVolume();
            // calcDerivative(lb, ub);
        }
    }
}