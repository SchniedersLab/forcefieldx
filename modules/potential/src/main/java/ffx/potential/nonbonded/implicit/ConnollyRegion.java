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
import static java.lang.System.arraycopy;
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
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;

import ffx.numerics.atomic.AtomicDoubleArray3D;
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
 * ConnollyRegion uses the algorithms from the AMS/VAM programs of
 * Michael Connolly to compute the analytical molecular surface
 * area and volume of a collection of spherical atoms; thus
 * it implements Fred Richards' molecular surface definition as
 * a set of analytically defined spherical and toroidal polygons.
 * <p>
 * Numerical derivatives of the volume are available.
 * <p>
 * Literature references:
 * M. L. Connolly, "Analytical Molecular Surface Calculation",
 * Journal of Applied Crystallography, 16, 548-558 (1983)
 * <p>
 * M. L. Connolly, "Computation of Molecular Volume", Journal
 * of the American Chemical Society, 107, 1118-1124 (1985)
 * <p>
 * C. E. Kundrot, J. W. Ponder and F. M. Richards, "Algorithms for
 * Calculating Excluded Volume and Its Derivatives as a Function
 * of Molecular Conformation and Their Use in Energy Minimization",
 * Journal of Computational Chemistry, 12, 402-409 (1991)
 *
 * @author Guowei Qi
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ConnollyRegion extends ParallelRegion {

    private static final Logger logger = Logger.getLogger(ConnollyRegion.class.getName());

    /**
     * Number of atoms in the system.
     */
    private final int nAtoms;
    /**
     * Array of atoms.
     */
    private Atom[] atoms;
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
    /**
     * If true, compute the gradient
     */
    private boolean gradient = false;
    /**
     * Array to store volume gradient
     **/
    private final double[][] dex;
    private final int[] itab;
    /**
     * A 3D array to store the gradient.
     */
    private AtomicDoubleArray3D grad;

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
     * Default size of a vector to randomly perturb coordinates.
     */
    public final static double DEFAULT_WIGGLE = 1.0e-8;
    /**
     * Size of a vector to randomly perturb coordinates.
     */
    private double wiggle = DEFAULT_WIGGLE;
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

    /**
     * 3D grid for finding atom neighbors.
     */
    private final static int MAXCUBE = 40;
    /**
     * maximum number of cycle convex edges.
     */
    private final static int MAXCYEP = 30;
    /**
     * maximum number of convex face cycles.
     */
    private final static int MAXFPCY = 10;
    /**
     * Maximum number of saddle faces.
     */
    private final int maxSaddleFace;
    /**
     * maximum number of convex edges.
     */
    private final int maxConvexEdges;
    /**
     * maximum number of neighboring atom pairs.
     */
    private final int maxAtomPairs;
    /**
     * maximum number of circles.
     */
    private final int maxCircles;
    /**
     * maximum number of total tori.
     */
    private final int maxTori;
    /**
     * maximum number of temporary tori.
     */
    private final int maxTempTori;
    /**
     * maximum number of concave edges.
     */
    private final int maxConcaveEdges;
    /**
     * maximum number of vertices.
     */
    private final int maxVertices;
    /**
     * maximum number of probe positions.
     */
    private final int maxProbePositions;
    /**
     * maximum number of concave faces.
     */
    private final int maxConcaveFaces;
    /**
     * maximum number of convex faces.
     */
    private final int maxConvexFaces;
    /**
     * maximum number of cycles.
     */
    private final int maxCycles;
    /**
     * Copy of the atomic coordinates. [X,Y,Z][Atom Number]
     */
    private final double[][] atomCoords;
    private final double[] x, y, z;
    /**
     * If true, atom has no free surface.
     */
    private final boolean[] noFreeSurface;
    /**
     * Atom free of neighbors.
     */
    private final boolean[] atomFreeOfNeighbors;
    /**
     * Atom buried.
     */
    private final boolean[] atomBuried;
    /**
     * Begin and end pointers for atoms neighbors.
     */
    private final int[][] atomNeighborPointers;
    /**
     * Atom numbers of neighbors.
     */
    private final int[] neighborAtomNumbers;
    /**
     * Pointer from neighbor to torus.
     */
    private final int[] neighborToTorus;
    /**
     * Number of temporary tori.
     */
    private int nTempTori;
    /**
     * Temporary torus atom numbers.
     */
    private final int[][] tempToriAtomNumbers;
    /**
     * First edge of each temporary torus.
     */
    private final int[] tempToriFirstEdge;
    /**
     * Last edge of each temporary torus.
     */
    private final int[] tempToriLastEdge;
    /**
     * Pointer to next edge of temporary torus.
     */
    private final int[] tempToriNextEdge;
    /**
     * Temporary torus buried.
     */
    private final boolean[] tempToriBuried;
    /**
     * Temporary torus free.
     */
    private final boolean[] tempToriFree;
    /**
     * Torus center.
     */
    private final double[][] torusCenter;
    /**
     * Torus radius.
     */
    private final double[] torusRadius;
    /**
     * Torus axis.
     */
    private final double[][] torusAxis;
    /**
     * Number of tori.
     */
    private int nTori;
    /**
     * Torus atom number.
     */
    private final int[][] torusAtomNumber;
    /**
     * Torus first edge.
     */
    private final int[] torusFirstEdge;
    /**
     * Torus free edge of neighbor.
     */
    private final boolean[] torusNeighborFreeEdge;
    /**
     * Probe position coordinates.
     */
    private final double[][] probeCoords;
    /**
     * Number of probe positions.
     */
    private int nProbePositions;
    /**
     * Probe position atom numbers.
     */
    private final int[][] probeAtomNumbers;
    /**
     * Vertex coordinates.
     */
    private final double[][] vertexCoords;
    /**
     * Number of concave faces.
     */
    private int nConcaveFaces;
    /**
     * Number of vertices.
     */
    private int nConcaveVerts;
    /**
     * Vertex atom number.
     */
    private final int[] vertexAtomNumbers;
    /**
     * Vertex probe number.
     */
    private final int[] vertexProbeNumber;
    /**
     * Number of concave edges.
     */
    private int nConcaveEdges;
    /**
     * Vertex numbers for each concave edge.
     */
    private final int[][] concaveEdgeVertexNumbers;
    /**
     * Concave face concave edge numbers.
     */
    private final int[][] concaveFaceEdgeNumbers;
    /**
     * Circle center.
     */
    private final double[][] circleCenter;
    /**
     * Circle radius.
     */
    private final double[] circleRadius;
    /**
     * Number of circles.
     */
    private int nCircles;
    /**
     * Circle atom number.
     */
    private final int[] circleAtomNumber;
    /**
     * Circle torus number.
     */
    private final int[] circleTorusNumber;
    /**
     * Number of convex edges.
     */
    private int nConvexEdges;
    /**
     * Convex edge circle number.
     */
    private final int[] convexEdgeCircleNumber;
    /**
     * Convex edge vertex number.
     */
    private final int[][] convexEdgeVertexNumber;
    /**
     * First convex edge of each atom.
     */
    private final int[] convexEdgeFirst;
    /**
     * Last convex edge of each atom.
     */
    private final int[] convexEdgeLast;
    /**
     * Pointer to next convex edge of atom.
     */
    private final int[] convexEdgeNext;
    /**
     * Number of saddle faces.
     */
    private int nSaddleFaces;
    /**
     * Saddle face concave edge numbers.
     */
    private final int[][] saddleConcaveEdgeNumbers;
    /**
     * Saddle face convex edge numbers.
     */
    private final int[][] saddleConvexEdgeNumbers;
    /**
     * Number of cycles.
     */
    private int nCycles;
    /**
     * Number of convex edges in cycle.
     */
    private final int[] convexEdgeCycleNum;
    /**
     * Cycle convex edge numbers.
     */
    private final int[][] convexEdgeCycleNumbers;
    /**
     * Number of convex faces.
     */
    private int nConvexFaces;
    /**
     * Atom number of convex face.
     */
    private final int[] convexFaceAtomNumber;
    /**
     * Convex face cycle numbers
     */
    private final int[][] convexFaceCycleNumbers;
    /**
     * Number of cycles bounding convex face.
     */
    private final int[] convexFaceNumCycles;

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

    private final ParallelTeam parallelTeam;

    /**
     * VolumeRegion constructor.
     *
     * @param atoms      Array of atom instances.
     * @param baseRadius Base radius for each atom (no added probe).
     * @param nThreads   Number of threads.
     */
    public ConnollyRegion(Atom[] atoms, double[] baseRadius, int nThreads) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.baseRadius = baseRadius;

        // Parallelization variables.
        parallelTeam = new ParallelTeam(1);
        volumeLoop = new VolumeLoop[nThreads];
        for (int i = 0; i < nThreads; i++) {
            volumeLoop[i] = new VolumeLoop();
        }
        sharedVolume = new SharedDouble();
        sharedArea = new SharedDouble();

        // Atom variables.
        radius = new double[nAtoms];
        skip = new boolean[nAtoms];
        atomCoords = new double[3][nAtoms];
        x = atomCoords[0];
        y = atomCoords[1];
        z = atomCoords[2];
        noFreeSurface = new boolean[nAtoms];
        atomFreeOfNeighbors = new boolean[nAtoms];
        atomBuried = new boolean[nAtoms];
        atomNeighborPointers = new int[2][nAtoms];

        // Atom pair variables.
        maxAtomPairs = 240 * nAtoms;
        neighborAtomNumbers = new int[maxAtomPairs];
        neighborToTorus = new int[maxAtomPairs];

        // Temporary Torus variables.
        maxTempTori = 120 * nAtoms;
        maxConcaveEdges = 12 * nAtoms;
        tempToriAtomNumbers = new int[2][maxTempTori];
        tempToriFirstEdge = new int[maxTempTori];
        tempToriLastEdge = new int[maxTempTori];
        tempToriBuried = new boolean[maxTempTori];
        tempToriFree = new boolean[maxTempTori];
        tempToriNextEdge = new int[maxConcaveEdges];

        // Torus variables.
        maxTori = 4 * nAtoms;
        torusCenter = new double[3][maxTori];
        torusRadius = new double[maxTori];
        torusAxis = new double[3][maxTori];
        torusAtomNumber = new int[2][maxTori];
        torusFirstEdge = new int[maxTori];
        torusNeighborFreeEdge = new boolean[maxTori];

        // Probe variables.
        maxProbePositions = 4 * nAtoms;
        probeCoords = new double[3][maxProbePositions];
        probeAtomNumbers = new int[3][maxProbePositions];

        // Vertex variables.
        maxVertices = 12 * nAtoms;
        vertexCoords = new double[3][maxVertices];
        vertexAtomNumbers = new int[maxVertices];
        vertexProbeNumber = new int[maxVertices];

        // Concave face variables.
        maxConcaveFaces = 4 * nAtoms;
        concaveFaceEdgeNumbers = new int[3][maxConcaveFaces];
        concaveEdgeVertexNumbers = new int[2][maxConcaveEdges];

        // Circle variables.
        maxCircles = 8 * nAtoms;
        circleCenter = new double[3][maxCircles];
        circleRadius = new double[maxCircles];
        circleAtomNumber = new int[maxCircles];
        circleTorusNumber = new int[maxCircles];

        // Convex face variables.
        maxConvexEdges = 12 * nAtoms;
        maxConvexFaces = 2 * nAtoms;
        maxCycles = 3 * nAtoms;
        convexEdgeCircleNumber = new int[maxConvexEdges];
        convexEdgeVertexNumber = new int[2][maxConvexEdges];
        convexEdgeFirst = new int[nAtoms];
        convexEdgeLast = new int[nAtoms];
        convexEdgeNext = new int[maxConvexEdges];
        convexEdgeCycleNum = new int[maxCycles];
        convexEdgeCycleNumbers = new int[MAXCYEP][maxCycles];
        convexFaceAtomNumber = new int[maxConvexFaces];
        convexFaceNumCycles = new int[maxConvexFaces];
        convexFaceCycleNumbers = new int[MAXFPCY][maxConvexFaces];

        // Saddle face variables.
        maxSaddleFace = 6 * nAtoms;
        saddleConcaveEdgeNumbers = new int[2][maxSaddleFace];
        saddleConvexEdgeNumbers = new int[2][maxSaddleFace];

        // Domain decomposition
        activeCube = new boolean[MAXCUBE * MAXCUBE * MAXCUBE];
        activeAdjacentCube = new boolean[MAXCUBE * MAXCUBE * MAXCUBE];
        firstAtomPointer = new int[MAXCUBE * MAXCUBE * MAXCUBE];
        cubeCoordinates = new int[3][nAtoms];
        nextAtomPointer = new int[nAtoms];

        // Volume derivative variables.
        dex = new double[3][nAtoms];
        itab = new int[nAtoms];
    }

    /**
     * Initialize this VolumeRegion instance for an energy evaluation.
     * <p>
     * Currently the number of atoms cannot change.
     *
     * @param atoms    Array of atoms.
     * @param gradient Compute the atomic coordinate gradient.
     * @param grad     Array to accumulate the gradient.
     */
    public void init(Atom[] atoms, boolean gradient, AtomicDoubleArray3D grad) {
        this.atoms = atoms;
        this.gradient = gradient;
        this.grad = grad;
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

    /**
     * Apply a random perturbation to the atomic coordinates
     * to avoid numerical instabilities for various linear, planar and
     * symmetric structures.
     *
     * @param wiggle Size of the random vector move in Angstroms.
     */
    public void setWiggle(double wiggle) {
        this.wiggle = wiggle;
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

    public double getCrossOver() {
        return crossOver;
    }

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

    /**
     * Execute the VolumeRegion with a private, single threaded ParallelTeam.
     */
    public void runVolume() {
        try {
            parallelTeam.execute(this);
        } catch (Exception e) {
            String message = " Fatal exception computing the Connolly surface area and volume.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    @Override
    public void start() {
        sharedVolume.set(0.0);
        sharedArea.set(0.0);
        double[] vector = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            getRandomVector(vector);
            Atom atom = atoms[i];
            atomCoords[0][i] = atom.getX() + (wiggle * vector[0]);
            atomCoords[1][i] = atom.getY() + (wiggle * vector[1]);
            atomCoords[2][i] = atom.getZ() + (wiggle * vector[2]);
        }
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

    /**
     * Construct the 3-dimensional random unit vector.
     *
     * @param vector The generated vector.
     */
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

    /**
     * Compute the cavitation energy.
     */
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
        double dReffdvdW = 1.0 / (pow(6.0, 2.0 / 3.0) * PI * vdWVolPI23);

        if (gradient && logger.isLoggable(Level.FINE)) {
            for (int i = 0; i < nAtoms; i++) {
                logger.fine(format(" Gradient %d (%16.8f, %16.8f, %16.8f)", i, dex[0][i], dex[1][i], dex[2][i]));
            }
        }

        // Find the cavitation energy using a combination of volume and surface area dependence.
        if (reff < volumeOff) {
            // Find cavity energy from only the molecular volume.
            cavitationEnergy = volumeEnergy;
            if (gradient) {
                for (int i = 0; i < nAtoms; i++) {
                    double dx = solventPressure * dex[0][i];
                    double dy = solventPressure * dex[1][i];
                    double dz = solventPressure * dex[2][i];
                    grad.add(0, i, dx, dy, dz);
                }
            }
        } else if (reff <= volumeCut) {
            // Include a tapered molecular volume.
            double taper = volumeSwitch.taper(reff, reff2, reff3, reff4, reff5);
            cavitationEnergy = taper * volumeEnergy;

            if (gradient) {
                double dtaper = volumeSwitch.dtaper(reff, reff2, reff3, reff4) * dReffdvdW;
                double factor = dtaper * volumeEnergy + taper * solventPressure;
                for (int i = 0; i < nAtoms; i++) {
                    double dx = factor * dex[0][i];
                    double dy = factor * dex[1][i];
                    double dz = factor * dex[2][i];
                    grad.add(0, i, dx, dy, dz);
                }
            }
        }

        // TODO: We need to include the surface area contribution to the gradient.
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
     * Compute Volume energy for a range of atoms.
     *
     * @since 1.0
     */
    private class VolumeLoop extends IntegerForLoop {
        /**
         * 3D grid for numeric volume derivatives.
         */
        private final static int MXCUBE = 15;
        /**
         * Maximum arcs for volume derivatives.
         */
        private final static int MAXARC = 1000;
        /**
         * Used by the place method to find probe sites.
         */
        private final static int MAXMNB = 500;

        private double localVolume;
        private double localSurfaceArea;

        @Override
        public void start() {
            localVolume = 0.0;
            localSurfaceArea = 0.0;
            if (gradient) {
                fill(dex[0], 0.0);
                fill(dex[1], 0.0);
                fill(dex[2], 0.0);
            }
        }

        @Override
        public void run(int lb, int ub) {
            setRadius();
            computeVolumeAndArea();
            if (gradient) {
                computeVolumeGradient();
            }
        }

        @Override
        public void finish() {
            sharedVolume.addAndGet(localVolume);
            sharedArea.addAndGet(localSurfaceArea);
        }

        /**
         * Assign van der Waals radii to the atoms; note that the radii
         * are incremented by the size of the probe.
         */
        private void setRadius() {
            for (int i = 0; i < nAtoms; i++) {
                if (baseRadius[i] == 0.0) {
                    radius[i] = 0.0;
                    skip[i] = true;
                } else {
                    skip[i] = false;
                    radius[i] = baseRadius[i] + exclude;
                }
            }
        }

        /**
         * Find the analytical volume and surface area.
         */
        private void computeVolumeAndArea() {
            nearby();
            torus();
            place();
            compress();
            saddles();
            contact();
            vam();
        }

        /**
         * The getTorus method tests for a possible torus position at the
         * interface between two atoms, and finds the torus radius, center
         * and axis.
         *
         * @param ia          First atom.
         * @param ja          Second atom.
         * @param torusCenter Torus center.
         * @param torusRadius Torus radius.
         * @param torusAxis   Torus axis.
         * @return
         */
        private boolean getTorus(int ia, int ja, double[] torusCenter, double[] torusRadius, double[] torusAxis) {
            boolean foundTorus = false;

            // Get the distance between the two atoms.
            double[] ai = new double[3];
            double[] aj = new double[3];
            getVector(ai, atomCoords, ia);
            getVector(aj, atomCoords, ja);
            double dij = dist(ai, aj);

            // Find a unit vector along interatomic (torus) axis.
            double[] vij = new double[3];
            double[] uij = new double[3];
            for (int k = 0; k < 3; k++) {
                vij[k] = atomCoords[k][ja] - atomCoords[k][ia];
                uij[k] = vij[k] / dij;
            }

            // Find coordinates of the center of the torus.
            double temp = 1.0 + ((radius[ia] + probe) * (radius[ia] + probe) - (radius[ja] + probe) * (radius[ja] + probe)) / (dij * dij);
            double[] bij = new double[3];
            for (int k = 0; k < 3; k++) {
                bij[k] = atomCoords[k][ia] + 0.5 * vij[k] * temp;
            }

            // Skip if atoms too far apart (should not happen).
            double temp1 = (radius[ia] + radius[ja] + 2.0 * probe) * (radius[ia] + radius[ja] + 2.0 * probe) - dij * dij;
            if (temp1 >= 0.0) {

                // Skip if one atom is inside the other.
                double temp2 = dij * dij - (radius[ia] - radius[ja]) * (radius[ia] - radius[ja]);
                if (temp2 >= 0.0) {

                    // Store the torus radius, center, and axis.
                    foundTorus = true;
                    torusRadius[0] = sqrt(temp1 * temp2) / (2.0 * dij);
                    for (int k = 0; k < 3; k++) {
                        torusCenter[k] = bij[k];
                        torusAxis[k] = uij[k];
                    }
                }
            }
            return foundTorus;
        }

        /**
         * The nearby method finds all of the through-space neighbors of
         * each atom for use in surface area and volume calculations.
         */
        private void nearby() {
            int maxclsa = 1000;
            int[] clsa = new int[maxclsa];
            // Temporary neighbor list, before sorting.
            int[] tempNeighborList = new int[maxclsa];
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
                    getVector(ai, atomCoords, i);
                    for (int j = i + 1; j < nAtoms; j++) {
                        getVector(aj, atomCoords, j);
                        double d2 = dist2(ai, aj);
                        double r2 = (radius[i] - radius[j]) * (radius[i] - radius[j]);
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
            double radmax = 0.0;
            for (int k = 0; k < 3; k++) {
                minAtomicCoordinates[k] = atomCoords[k][0];
            }
            for (int i = 0; i < nAtoms; i++) {
                for (int k = 0; k < 3; k++) {
                    if (atomCoords[k][i] < minAtomicCoordinates[k]) {
                        minAtomicCoordinates[k] = atomCoords[k][i];
                    }
                }
                if (radius[i] > radmax) {
                    radmax = radius[i];
                }
            }

            // Calculate width of cube from maximum atom radius and probe radius.
            double width = 2.0 * (radmax + probe);

            // Set up cube arrays; first the integer coordinate arrays.
            for (int i = 0; i < nAtoms; i++) {
                for (int k = 0; k < 3; k++) {
                    cubeCoordinates[k][i] = (int) ((atomCoords[k][i] - minAtomicCoordinates[k]) / width);
                    if (cubeCoordinates[k][i] < 0) {
                        throw new EnergyException(" Cube Coordinate Too Small", false);
                    } else if (cubeCoordinates[k][i] > MAXCUBE) {
                        throw new EnergyException(" Cube Coordinate Too Large", false);
                    }
                }
            }

            // Initialize head pointer and srn=2 arrays.
            fill(firstAtomPointer, -1);
            fill(activeCube, false);
            fill(activeAdjacentCube, false);

            // Initialize linked list pointers.
            fill(nextAtomPointer, -1);

            // Set up head and later pointers for each atom.
            atomLoop:
            for (int iatom = 0; iatom < nAtoms; iatom++) {
                // Skip atoms with surface request numbers of zero.
                if (skip[iatom]) {
                    continue;
                }

                getVector(ai, atomCoords, iatom);
                int i = cubeCoordinates[0][iatom];
                int j = cubeCoordinates[1][iatom];
                int k = cubeCoordinates[2][iatom];

                if (firstAtomPointer[index(i, j, k)] <= -1) {
                    // First atom in this cube.
                    firstAtomPointer[index(i, j, k)] = iatom;
                } else {
                    // Add to end of linked list.
                    int iptr = firstAtomPointer[index(i, j, k)];
                    boolean done = false;
                    while (!done) {
                        getVector(aj, atomCoords, iptr);
                        // Check for duplicate atoms, turn off one of them.
                        if (dist2(ai, aj) <= 0.0) {
                            skip[iatom] = true;
                            continue atomLoop;
                        }
                        // Move on down the list.
                        if (nextAtomPointer[iptr] <= -1.0) {
                            done = true;
                        } else {
                            iptr = nextAtomPointer[iptr];
                        }
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
                            checkCube(i, j, k);
                        }
                    }
                }
            }

            int ncls = -1;

            // Zero pointers for atom and find its cube.
            for (int i = 0; i < nAtoms; i++) {
                int nclsa = -1;
                noFreeSurface[i] = skip[i];
                atomNeighborPointers[0][i] = -1;
                atomNeighborPointers[1][i] = -1;
                if (skip[i]) {
                    continue;
                }
                int ici = cubeCoordinates[0][i];
                int icj = cubeCoordinates[1][i];
                int ick = cubeCoordinates[2][i];

                // Skip iatom if its cube and adjoining cubes contain only blockers.
                if (!activeAdjacentCube[index(ici, icj, ick)]) {
                    continue;
                }
                double sumi = 2.0 * probe + radius[i];

                // Check iatom cube and adjacent cubes for neighboring atoms.
                for (int jck = max(ick - 1, 0); jck < min(ick + 2, MAXCUBE); jck++) {
                    for (int jcj = max(icj - 1, 0); jcj < min(icj + 2, MAXCUBE); jcj++) {
                        // TODO: Try to remove this labeled for loop.
                        for4:
                        for (int jci = max(ici - 1, 0); jci < min(ici + 2, MAXCUBE); jci++) {
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
                                    double sum = sumi + radius[j];
                                    double vect1 = abs(atomCoords[0][j] - atomCoords[0][i]);
                                    if (vect1 >= sum) {
                                        continue;
                                    }
                                    double vect2 = abs(atomCoords[1][j] - atomCoords[1][i]);
                                    if (vect2 >= sum) {
                                        continue;
                                    }
                                    double vect3 = abs(atomCoords[2][j] - atomCoords[2][i]);
                                    if (vect3 >= sum) {
                                        continue;
                                    }
                                    double d2 = (vect1 * vect1) + (vect2 * vect2) + (vect3 * vect3);
                                    if (d2 >= sum * sum) {
                                        continue;
                                    }

                                    // Atoms are neighbors, save atom number in temporary array.
                                    if (!skip[j]) {
                                        noFreeSurface[i] = false;
                                    }
                                    nclsa++;
                                    if (nclsa >= maxclsa) {
                                        throw new EnergyException(" Too many Neighbors for Atom", false);
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

                if (noFreeSurface[i]) {
                    continue;
                }

                // Set up neighbors arrays with jatom in increasing order.
                int jmold = -1;
                for (int juse = 0; juse <= nclsa; juse++) {
                    int jmin = nAtoms;
                    int jmincls = 0;
                    for (int jcls = 0; jcls < nclsa + 1; jcls++) {
                        // Don't use ones already sorted.
                        if (tempNeighborList[jcls] > jmold) {
                            if (tempNeighborList[jcls] < jmin) {
                                jmin = tempNeighborList[jcls];
                                jmincls = jcls;
                            }
                        }
                    }
                    jmold = jmin;
                    int jcls = jmincls;
                    int jatom = tempNeighborList[jcls];
                    clsa[juse] = jatom;
                }

                // Set up pointers to first and last neighbors of atom.
                if (nclsa > -1) {
                    atomNeighborPointers[0][i] = ncls + 1;
                    for (int m = 0; m < nclsa + 1; m++) {
                        ncls++;
                        if (ncls >= maxAtomPairs) {
                            throw new EnergyException(" Too many Neighboring Atom Pairs", false);
                        }
                        neighborAtomNumbers[ncls] = clsa[m];
                    }
                    atomNeighborPointers[1][i] = ncls;
                }
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
        private int index(int i, int j, int k) {
            return k + MAXCUBE * (j + MAXCUBE * i);
        }

        /**
         * Check if this cube or any adjacent cube has active atoms.
         *
         * @param i X-index for activeAdjacentCube.
         * @param j Y-index for activeAdjacentCube.
         * @param k Z-index for activeAdjacentCube.
         */
        private void checkCube(int i, int j, int k) {
            for (int k1 = max(k - 1, 0); k1 < min(k + 2, MAXCUBE); k1++) {
                for (int j1 = max(j - 1, 0); j1 < min(j + 2, MAXCUBE); j1++) {
                    for (int i1 = max(i - 1, 0); i1 < min(i + 2, MAXCUBE); i1++) {
                        if (activeCube[index(i1, j1, k1)]) {
                            activeAdjacentCube[index(i, j, k)] = true;
                            return;
                        }
                    }
                }
            }
        }

        /**
         * The torus method sets a list of all of the temporary torus
         * positions by testing for a torus between each atom and its
         * neighbors
         */
        private void torus() {
            // No torus is possible if there is only one atom.
            nTempTori = -1;
            fill(atomFreeOfNeighbors, true);
            if (nAtoms <= 1) {
                return;
            }

            double[] tt = new double[3];
            double[] ttax = new double[3];
            // Get beginning and end pointers to neighbors of this atom.
            for (int ia = 0; ia < nAtoms; ia++) {
                if (!noFreeSurface[ia]) {
                    int ibeg = atomNeighborPointers[0][ia];
                    int iend = atomNeighborPointers[1][ia];
                    if (ibeg > -1) {
                        // Check for no neighbors.
                        for (int jn = ibeg; jn <= iend; jn++) {
                            // Clear pointer from neighbor to torus.
                            neighborToTorus[jn] = -1;
                            // Get atom number of neighbor.
                            int ja = neighborAtomNumbers[jn];
                            // Don't create torus twice.
                            if (ja >= ia) {
                                // Do some solid geometry.
                                double[] ttr = {0.0};
                                boolean ttok = getTorus(ia, ja, tt, ttr, ttax);
                                if (ttok) {
                                    // We have temporary torus; set up variables.
                                    nTempTori++;
                                    if (nTempTori >= maxTempTori) {
                                        throw new EnergyException(" Too many Temporary Tori", false);
                                    }
                                    // Mark both atoms not free.
                                    atomFreeOfNeighbors[ia] = false;
                                    atomFreeOfNeighbors[ja] = false;
                                    tempToriAtomNumbers[0][nTempTori] = ia;
                                    tempToriAtomNumbers[1][nTempTori] = ja;
                                    // Pointer from neighbor to torus.
                                    neighborToTorus[jn] = nTempTori;
                                    // Initialize torus as both free and buried.
                                    tempToriFree[nTempTori] = true;
                                    tempToriBuried[nTempTori] = true;
                                    // Clear pointers from torus to first and last concave edges.
                                    tempToriFirstEdge[nTempTori] = -1;
                                    tempToriLastEdge[nTempTori] = -1;
                                }
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
        private void place() {
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

            nProbePositions = -1;
            nConcaveFaces = -1;
            nConcaveEdges = -1;
            nConcaveVerts = -1;

            // No possible placement if there are no temporary tori.
            if (nTempTori <= -1) {
                return;
            }

            // Consider each torus in turn.
            for (int itt = 0; itt < nTempTori + 1; itt++) {
                // Get atom numbers.
                int ia = tempToriAtomNumbers[0][itt];
                int ja = tempToriAtomNumbers[1][itt];

                // Form mutual neighbor list; clear number of mutual neighbors of atoms ia and ja.
                int nmnb = -1;

                // Get beginning and end pointers for each atom's neighbor list.
                int iptr = atomNeighborPointers[0][ia];
                int jptr = atomNeighborPointers[0][ja];

                // For loops like this help replace go to statements from the Tinker code.
                outerloop:
                for (int i = 0; i < 1; i++) {
                    if (iptr <= -1 || jptr <= -1) {
                        continue;
                    }
                    int iend = atomNeighborPointers[1][ia];
                    int jend = atomNeighborPointers[1][ja];

                    // Collect mutual neighbors.
                    int counter = 1;
                    int counter2 = 1;
                    // TODO: Try to turn this into a while loop.
                    for (int t = 0; t < counter2; t++) {
                        // TODO: Try to turn this into a while loop.
                        for (int q = 0; q < counter; q++) {
                            // Check for end of loop.
                            if (iptr > iend || jptr > jend) {
                                continue outerloop;
                            }
                            // Go move the lagging pointer.
                            if (neighborAtomNumbers[iptr] < neighborAtomNumbers[jptr]) {
                                iptr++;
                                counter++;
                                continue;
                            }
                            if (neighborAtomNumbers[jptr] < neighborAtomNumbers[iptr]) {
                                jptr++;
                                counter++;
                            }
                        }

                        /*
                          Both point at same neighbor; one more mutual
                          neighbor save atom number of mutual neighbor.
                         */
                        nmnb++;
                        if (nmnb >= MAXMNB) {
                            throw new EnergyException(" Too many Mutual Neighbors", false);
                        }
                        mnb[nmnb] = neighborAtomNumbers[iptr];

                        // Save pointers to second and third tori.
                        ikt[nmnb] = neighborToTorus[iptr];
                        jkt[nmnb] = neighborToTorus[jptr];
                        iptr++;
                        counter2++;
                    }
                }
                // We have all the mutual neighbors of ia and ja if no mutual neighbors, skip to end of loop.
                if (nmnb <= -1) {
                    tempToriBuried[itt] = false;
                    continue;
                }
                double[] hij = {0.0};
                boolean ttok = getTorus(ia, ja, bij, hij, uij);
                for (int km = 0; km < nmnb + 1; km++) {
                    int ka = mnb[km];
                    getVector(ak, atomCoords, ka);
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
                    tok = getTorus(ia, ja, tij, rij, uij);
                    for (int w = 0; w < 1; w++) {
                        if (tok) {
                            getVector(ai, atomCoords, ka);
                            double dat2 = dist2(ai, tij);
                            double rad2 = (radius[ka] + probe) * (radius[ka] + probe) - rij[0] * rij[0];

                            // If "ka" less than "ja", then all we care about is whether the torus is buried.
                            if (ka < ja) {
                                if (rad2 <= 0.0 || dat2 > rad2) {
                                    continue;
                                }
                            }

                            double[] rik = {0.0};
                            tok = getTorus(ia, ka, tik, rik, uik);
                            if (!tok) {
                                continue;
                            }
                            double dotijk = check(dot(uij, uik));
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
                            getVector(ai, atomCoords, ia);
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
                        tempToriBuried[itt] = true;
                        tempToriFree[itt] = false;
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
                        tempToriFree[itt] = false;
                        tempToriFree[ik] = false;
                        tempToriFree[jk] = false;

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
                                getVector(ak, atomCoords, la);
                                // Compare distance to sum of radii.
                                if (dist2(pijk, ak) <= sumcls[lm]) {
                                    continue for3;
                                }
                            }
                            lm = lkcls[lm];
                            counter3++;
                        }
                        // We have a new probe position.
                        nProbePositions++;
                        if (nProbePositions >= maxProbePositions) {
                            throw new EnergyException(" Too many Probe Positions", false);
                        }
                        // Mark three tori not buried.
                        tempToriBuried[itt] = false;
                        tempToriBuried[ik] = false;
                        tempToriBuried[jk] = false;
                        // Store probe center.
                        for (int k = 0; k < 3; k++) {
                            probeCoords[k][nProbePositions] = pijk[k];
                        }

                        // Calculate vectors from probe to atom centers.
                        if (nConcaveVerts + 3 >= maxVertices) {
                            throw new EnergyException(" Too many Vertices", false);
                        }
                        for (int k = 0; k < 3; k++) {
                            vertexCoords[k][nConcaveVerts + 1] = atomCoords[k][ia] - probeCoords[k][nProbePositions];
                            vertexCoords[k][nConcaveVerts + 2] = atomCoords[k][ja] - probeCoords[k][nProbePositions];
                            vertexCoords[k][nConcaveVerts + 3] = atomCoords[k][ka] - probeCoords[k][nProbePositions];
                        }
                        //double matrix[] = new double[9];
                        //int a = 0;
                        //for (int b = 0; b < 3; b++) {
                        //    for (int c = 0; c < 3; c++) {
                        //        matrix[a++] = v[b][nv + c + 1];
                        //    }
                        //}

                        // Calculate determinant of vectors defining triangle.
                        double det = vertexCoords[0][nConcaveVerts + 1] * vertexCoords[1][nConcaveVerts + 2] * vertexCoords[2][nConcaveVerts + 3]
                                + vertexCoords[0][nConcaveVerts + 2] * vertexCoords[1][nConcaveVerts + 3] * vertexCoords[2][nConcaveVerts + 1]
                                + vertexCoords[0][nConcaveVerts + 3] * vertexCoords[1][nConcaveVerts + 1] * vertexCoords[2][nConcaveVerts + 2]
                                - vertexCoords[0][nConcaveVerts + 3] * vertexCoords[1][nConcaveVerts + 2] * vertexCoords[2][nConcaveVerts + 1]
                                - vertexCoords[0][nConcaveVerts + 2] * vertexCoords[1][nConcaveVerts + 1] * vertexCoords[2][nConcaveVerts + 3]
                                - vertexCoords[0][nConcaveVerts + 1] * vertexCoords[1][nConcaveVerts + 3] * vertexCoords[2][nConcaveVerts + 2];
                        // Now add probe coordinates to vertices.
                        for (int k = 0; k < 3; k++) {
                            vertexCoords[k][nConcaveVerts + 1] = probeCoords[k][nProbePositions] + (vertexCoords[k][nConcaveVerts + 1] * probe / (radius[ia] + probe));
                            vertexCoords[k][nConcaveVerts + 2] = probeCoords[k][nProbePositions] + (vertexCoords[k][nConcaveVerts + 2] * probe / (radius[ja] + probe));
                            vertexCoords[k][nConcaveVerts + 3] = probeCoords[k][nProbePositions] + (vertexCoords[k][nConcaveVerts + 3] * probe / (radius[ka] + probe));
                        }
                        // Want the concave face to have counter-clockwise orientation.
                        if (det > 0.0) {
                            // Swap second and third vertices.
                            for (int k = 0; k < 3; k++) {
                                tempv[k] = vertexCoords[k][nConcaveVerts + 2];
                                vertexCoords[k][nConcaveVerts + 2] = vertexCoords[k][nConcaveVerts + 3];
                                vertexCoords[k][nConcaveVerts + 3] = tempv[k];
                            }
                            // Set up pointers from probe to atoms.
                            probeAtomNumbers[0][nProbePositions] = ia;
                            probeAtomNumbers[1][nProbePositions] = ka;
                            probeAtomNumbers[2][nProbePositions] = ja;
                            // Set pointers from vertices to atoms.
                            vertexAtomNumbers[nConcaveVerts + 1] = ia;
                            vertexAtomNumbers[nConcaveVerts + 2] = ka;
                            vertexAtomNumbers[nConcaveVerts + 3] = ja;
                            // Insert concave edges into linked lists for appropriate tori.
                            insertEdge(nConcaveEdges + 1, ik);
                            insertEdge(nConcaveEdges + 2, jk);
                            insertEdge(nConcaveEdges + 3, itt);
                        } else {
                            // Similarly, if face already counter-clockwise.
                            probeAtomNumbers[0][nProbePositions] = ia;
                            probeAtomNumbers[1][nProbePositions] = ja;
                            probeAtomNumbers[2][nProbePositions] = ka;
                            vertexAtomNumbers[nConcaveVerts + 1] = ia;
                            vertexAtomNumbers[nConcaveVerts + 2] = ja;
                            vertexAtomNumbers[nConcaveVerts + 3] = ka;
                            insertEdge(nConcaveEdges + 1, itt);
                            insertEdge(nConcaveEdges + 2, jk);
                            insertEdge(nConcaveEdges + 3, ik);
                        }
                        // Set up pointers from vertices to probe.
                        for (int kv = 1; kv < 4; kv++) {
                            vertexProbeNumber[nConcaveVerts + kv] = nProbePositions;
                        }
                        // Set up concave edges and concave face.
                        if (nConcaveEdges + 3 >= maxConcaveEdges) {
                            throw new EnergyException(" Too many Concave Edges", false);
                        }
                        // Edges point to vertices.
                        concaveEdgeVertexNumbers[0][nConcaveEdges + 1] = nConcaveVerts + 1;
                        concaveEdgeVertexNumbers[1][nConcaveEdges + 1] = nConcaveVerts + 2;
                        concaveEdgeVertexNumbers[0][nConcaveEdges + 2] = nConcaveVerts + 2;
                        concaveEdgeVertexNumbers[1][nConcaveEdges + 2] = nConcaveVerts + 3;
                        concaveEdgeVertexNumbers[0][nConcaveEdges + 3] = nConcaveVerts + 3;
                        concaveEdgeVertexNumbers[1][nConcaveEdges + 3] = nConcaveVerts + 1;
                        if (nConcaveFaces + 1 >= maxConcaveFaces) {
                            throw new EnergyException(" Too many Concave Faces", false);
                        }
                        // Face points to edges.
                        for (int ke = 0; ke < 3; ke++) {
                            concaveFaceEdgeNumbers[ke][nConcaveFaces + 1] = nConcaveEdges + ke + 1;
                        }
                        // Increment counters for number of faces, edges, and vertices.
                        nConcaveFaces++;
                        nConcaveEdges += 3;
                        nConcaveVerts += 3;
                    }
                }
            }
        }

        /**
         * The insertEdge method inserts a concave edge into the linked list for
         * its temporary torus.
         *
         * @param edgeNumber  The concave edge number.
         * @param torusNumber The temporary torus.
         */
        private void insertEdge(int edgeNumber, int torusNumber) {
            // Check for a serious error in the calling arguments.
            if (edgeNumber <= -1) {
                throw new EnergyException(" Bad Edge Number in INEDGE", false);
            }
            if (torusNumber <= -1) {
                throw new EnergyException(" Bad Torus Number in INEDGE", false);
            }
            // Set beginning of list or add to end.
            if (tempToriFirstEdge[torusNumber] == -1) {
                tempToriFirstEdge[torusNumber] = edgeNumber;
                tempToriNextEdge[edgeNumber] = -1;
                tempToriLastEdge[torusNumber] = edgeNumber;
            } else {
                tempToriNextEdge[tempToriLastEdge[torusNumber]] = edgeNumber;
                tempToriNextEdge[edgeNumber] = -1;
                tempToriLastEdge[torusNumber] = edgeNumber;
            }
        }

        /**
         * The compress method transfers only the non-buried tori from the
         * temporary tori arrays to the final tori arrays.
         */
        private void compress() {
            double[] torCenter = new double[3];
            double[] torAxis = new double[3];
            // Initialize the number of non-buried tori.
            nTori = -1;
            if (nTempTori <= -1) {
                return;
            }
            // If torus is free, then it is not buried; skip to end of loop if buried torus.
            double[] torRad = {0};
            for (int itt = 0; itt <= nTempTori; itt++) {
                if (tempToriFree[itt]) {
                    tempToriBuried[itt] = false;
                }
                if (!tempToriBuried[itt]) {
                    // First, transfer information.
                    nTori++;
                    if (nTori >= maxTori) {
                        throw new EnergyException(" Too many non-buried tori.", false);
                    }
                    int ia = tempToriAtomNumbers[0][itt];
                    int ja = tempToriAtomNumbers[1][itt];
                    getVector(torCenter, torusCenter, nTori);
                    getVector(torAxis, torusAxis, nTori);
                    getTorus(ia, ja, torCenter, torRad, torAxis);
                    torusCenter[0][nTori] = torCenter[0];
                    torusCenter[1][nTori] = torCenter[1];
                    torusCenter[2][nTori] = torCenter[2];
                    torusAxis[0][nTori] = torAxis[0];
                    torusAxis[1][nTori] = torAxis[1];
                    torusAxis[2][nTori] = torAxis[2];
                    torusRadius[nTori] = torRad[0];
                    torusAtomNumber[0][nTori] = ia;
                    torusAtomNumber[1][nTori] = ja;
                    torusNeighborFreeEdge[nTori] = tempToriFree[itt];
                    torusFirstEdge[nTori] = tempToriFirstEdge[itt];
                    // Special check for inconsistent probes.
                    int iptr = torusFirstEdge[nTori];
                    int ned = -1;
                    while (iptr != -1) {
                        ned++;
                        iptr = tempToriNextEdge[iptr];
                    }
                    if ((ned % 2) == 0) {
                        iptr = torusFirstEdge[nTori];
                        while (iptr != -1) {
                            int iv1 = concaveEdgeVertexNumbers[0][iptr];
                            int iv2 = concaveEdgeVertexNumbers[1][iptr];
                            int ip1 = vertexProbeNumber[iv1];
                            int ip2 = vertexProbeNumber[iv2];
                            logger.warning(format(" Odd Torus for Probes IP1 %d and IP2 %d", ip1, ip2));
                            iptr = tempToriNextEdge[iptr];
                        }
                    }
                }
            }
        }

        /**
         * The saddles method constructs circles, convex edges, and saddle faces.
         */
        private void saddles() {
            final int maxent = 500;
            int[] ten = new int[maxent];
            int[] nxtang = new int[maxent];
            double[] ai = new double[3];
            double[] aj = new double[3];
            double[] ak = new double[3];
            double[] atvect = new double[3];
            double[] teang = new double[maxent];
            double[][] tev = new double[3][maxent];
            boolean[] sdstrt = new boolean[maxent];

            // Zero the number of circles, convex edges, and saddle faces.
            nCircles = -1;
            nConvexEdges = -1;
            nSaddleFaces = -1;
            fill(convexEdgeFirst, -1);
            fill(convexEdgeLast, -1);
            fill(atomBuried, true);

            // No saddle faces if no tori.
            if (nTori < 0) {
                return;
            }

            // Cycle through tori.
            for (int it = 0; it <= nTori; it++) {
                if (skip[torusAtomNumber[0][it]] && skip[torusAtomNumber[1][it]]) {
                    continue;
                }

                // Set up two circles.
                for (int in = 0; in < 2; in++) {
                    int ia = torusAtomNumber[in][it];

                    // Mark atom not buried.
                    atomBuried[ia] = false;

                    // Vector from atom to torus center.
                    for (int k = 0; k < 3; k++) {
                        atvect[k] = torusCenter[k][it] - atomCoords[k][ia];
                    }
                    double factor = radius[ia] / (radius[ia] + probe);

                    // One more circle.
                    nCircles++;
                    if (nCircles >= maxCircles) {
                        throw new EnergyException(" Too many Circles", false);
                    }

                    // Circle center.
                    for (int k = 0; k < 3; k++) {
                        circleCenter[k][nCircles] = atomCoords[k][ia] + factor * atvect[k];
                    }

                    // Pointer from circle to atom.
                    circleAtomNumber[nCircles] = ia;

                    // Pointer from circle to torus.
                    circleTorusNumber[nCircles] = it;

                    // Circle radius.
                    circleRadius[nCircles] = factor * torusRadius[it];
                }

                // Skip to special "free torus" code.
                if (!torusNeighborFreeEdge[it]) {
                    /*
                      Now we collect all the concave edges for this
                      torus; for each concave edge, calculate vector
                      from torus center through probe center and the
                      angle relative to first such vector.
                     */

                    // Clear the number of concave edges for torus.
                    int nent = -1;

                    // Pointer to start of linked list.
                    int ien = torusFirstEdge[it];
                    // 10 Continue
                    while (ien >= 0) {
                        // One more concave edge.
                        nent++;
                        if (nent >= maxent) {
                            throw new EnergyException(" Too many Edges for Torus", false);
                        }
                        // First vertex of edge.
                        int iv = concaveEdgeVertexNumbers[0][ien];

                        // Probe number of vertex.
                        int ip = vertexProbeNumber[iv];
                        for (int k = 0; k < 3; k++) {
                            tev[k][nent] = probeCoords[k][ip] - torusCenter[k][it];
                        }
                        double dtev = 0.0;
                        for (int k = 0; k < 3; k++) {
                            dtev += tev[k][nent] * tev[k][nent];
                        }
                        if (dtev <= 0.0) {
                            throw new EnergyException(" Probe on Torus Axis", false);
                        }
                        dtev = sqrt(dtev);
                        for (int k = 0; k < 3; k++) {
                            tev[k][nent] = tev[k][nent] / dtev;
                        }

                        // Store concave edge number.
                        ten[nent] = ien;
                        if (nent > 0) {
                            // Calculate angle between this vector and first vector.
                            double dt = 0.0;
                            for (int k = 0; k < 3; k++) {
                                dt += tev[k][0] * tev[k][nent];
                            }
                            if (dt > 1.0) {
                                dt = 1.0;
                            }
                            if (dt < -1.0) {
                                dt = -1.0;
                            }

                            // Store angle.
                            teang[nent] = acos(dt);

                            ai[0] = tev[0][0];
                            ai[1] = tev[1][0];
                            ai[2] = tev[2][0];
                            aj[0] = tev[0][nent];
                            aj[1] = tev[1][nent];
                            aj[2] = tev[2][nent];
                            ak[0] = torusAxis[0][it];
                            ak[1] = torusAxis[1][it];
                            ak[2] = torusAxis[2][it];

                            // Get the sign right.
                            double triple = tripleProduct(ai, aj, ak);
                            if (triple < 0.0) {
                                teang[nent] = 2.0 * Math.PI - teang[nent];
                            }
                        } else {
                            teang[0] = 0.0;
                        }
                        // Saddle face starts with this edge if it points parallel to torus axis vector
                        // (which goes from first to second atom).
                        sdstrt[nent] = (vertexAtomNumbers[iv] == torusAtomNumber[0][it]);
                        // Next edge in list.
                        ien = tempToriNextEdge[ien];
                    }

                    if (nent <= -1) {
                        logger.severe(" No Edges for Non-free Torus");
                    }

                    if ((nent % 2) == 0) {
                        throw new EnergyException(" Odd Number of Edges for Torus", false);
                    }

                    // Set up linked list of concave edges in order of increasing angle
                    // around the torus axis;
                    // clear second linked (angle-ordered) list pointers.
                    for (int ient = 0; ient <= nent; ient++) {
                        nxtang[ient] = -1;
                    }
                    for (int ient = 1; ient <= nent; ient++) {
                        // We have an entry to put into linked list search for place to put it.
                        int l1 = -1;
                        int l2 = 0;
                        while (l2 != -1 && teang[ient] >= teang[l2]) {
                            l1 = l2;
                            l2 = nxtang[l2];
                        }
                        // We are at end of linked list or between l1 and l2; insert edge.
                        if (l1 <= -1) {
                            throw new EnergyException(" Logic Error in SADDLES", true);
                        }
                        nxtang[l1] = ient;
                        nxtang[ient] = l2;
                    }

                    // Collect pairs of concave edges into saddles.
                    // Create convex edges while you're at it.
                    int l1 = 0;
                    while (l1 > -1) {
                        // Check for start of saddle.
                        if (sdstrt[l1]) {
                            // One more saddle face.
                            nSaddleFaces++;
                            if (nSaddleFaces >= maxSaddleFace) {
                                throw new EnergyException(" Too many Saddle Faces", false);
                            }
                            // Get edge number.
                            ien = ten[l1];
                            // First concave edge of saddle.
                            saddleConcaveEdgeNumbers[0][nSaddleFaces] = ien;
                            // One more convex edge.
                            nConvexEdges++;
                            if (nConvexEdges >= maxConvexEdges) {
                                throw new EnergyException(" Too many Convex Edges", false);
                            }
                            // First convex edge points to second circle.
                            convexEdgeCircleNumber[nConvexEdges] = nCircles;
                            // Atom circle lies on.
                            int ia = circleAtomNumber[nCircles];
                            // Insert convex edge into linked list for atom.
                            insertConvexEdgeForAtom(nConvexEdges, ia);
                            // First vertex of convex edge is second vertex of concave edge.
                            convexEdgeVertexNumber[0][nConvexEdges] = concaveEdgeVertexNumbers[1][ien];
                            // First convex edge of saddle.
                            saddleConvexEdgeNumbers[0][nSaddleFaces] = nConvexEdges;
                            // One more convex edge.
                            nConvexEdges++;
                            if (nConvexEdges >= maxConvexEdges) {
                                throw new EnergyException(" Too many Convex Edges", false);
                            }
                            // Second convex edge points to fist circle.
                            convexEdgeCircleNumber[nConvexEdges] = nCircles - 1;
                            ia = circleAtomNumber[nCircles - 1];
                            // Insert convex edge into linked list for atom.
                            insertConvexEdgeForAtom(nConvexEdges, ia);
                            // Second vertex of second convex edge is first vertex of first concave edge.
                            convexEdgeVertexNumber[1][nConvexEdges] = concaveEdgeVertexNumbers[0][ien];
                            l1 = nxtang[l1];
                            // Wrap around.
                            if (l1 <= -1) {
                                l1 = 0;
                            }
                            if (sdstrt[l1]) {
                                int m1 = nxtang[l1];
                                if (m1 <= -1) {
                                    m1 = 0;
                                }
                                if (sdstrt[m1]) {
                                    throw new EnergyException(" Three Starts in a Row", false);
                                }
                                int n1 = nxtang[m1];
                                // The old switcheroo.
                                nxtang[l1] = n1;
                                nxtang[m1] = l1;
                                l1 = m1;
                            }
                            ien = ten[l1];
                            // Second concave edge for saddle face.
                            saddleConcaveEdgeNumbers[1][nSaddleFaces] = ien;
                            // Second vertex of first convex edge is first vertex of second concave edge.
                            convexEdgeVertexNumber[1][nConvexEdges - 1] = concaveEdgeVertexNumbers[0][ien];
                            // First vertex of second convex edge is second vertex of second concave edge.
                            convexEdgeVertexNumber[0][nConvexEdges] = concaveEdgeVertexNumbers[1][ien];
                            saddleConvexEdgeNumbers[1][nSaddleFaces] = nConvexEdges;
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
                    nSaddleFaces++;
                    if (nSaddleFaces >= maxSaddleFace) {
                        throw new EnergyException(" Too many Saddle Faces", false);
                    }
                    // No concave edge for saddle.
                    saddleConcaveEdgeNumbers[0][nSaddleFaces] = -1;
                    saddleConcaveEdgeNumbers[1][nSaddleFaces] = -1;
                    // One more convex edge.
                    nConvexEdges++;
                    int ia = circleAtomNumber[nCircles];
                    // Insert convex edge into linked list of atom.
                    insertConvexEdgeForAtom(nConvexEdges, ia);
                    // No vertices for convex edge.
                    convexEdgeVertexNumber[0][nConvexEdges] = -1;
                    convexEdgeVertexNumber[1][nConvexEdges] = -1;
                    // Pointer from convex edge to second circle.
                    convexEdgeCircleNumber[nConvexEdges] = nCircles;
                    // First convex edge for saddle face.
                    saddleConvexEdgeNumbers[0][nSaddleFaces] = nConvexEdges;
                    // One more convex edge.
                    nConvexEdges++;
                    ia = circleAtomNumber[nCircles - 1];
                    // Insert second convex edge into linked list.
                    insertConvexEdgeForAtom(nConvexEdges, ia);
                    // No vertices for convex edge.
                    convexEdgeVertexNumber[0][nConvexEdges] = -1;
                    convexEdgeVertexNumber[1][nConvexEdges] = -1;
                    // Convex edge points to first circle.
                    convexEdgeCircleNumber[nConvexEdges] = nCircles - 1;
                    // Second convex edge for saddle face.
                    saddleConvexEdgeNumbers[1][nSaddleFaces] = nConvexEdges;
                    // Buried torus; do nothing with it.
                }
            }

            logger.fine(format(" Number of circles      %d", nCircles + 1));
            logger.fine(format(" Number of convex edges %d", nConvexEdges + 1));
            logger.fine(format(" Number of saddle faces %d", nSaddleFaces + 1));
        }

        /**
         * The contact method constructs the contact surface, cycles and convex faces.
         */
        private void contact() {
            final int maxepa = 300;
            final int maxcypa = 100;
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
            nCycles = -1;
            nConvexFaces = -1;

            // Mark all free atoms not buried.
            for (int ia = 0; ia < nAtoms; ia++) {
                if (atomFreeOfNeighbors[ia]) {
                    atomBuried[ia] = false;
                }
            }

            // Go through all atoms.
            atomLoop:
            for (int ia = 0; ia < nAtoms; ia++) {
                if (skip[ia] || atomBuried[ia]) {
                    continue;
                }

                // Special code for completely solvent-accessible atom.
                for (int o = 0; o < 1; o++) {
                    if (atomFreeOfNeighbors[ia]) {
                        continue;
                    }
                    // Gather convex edges for atom Clear number of convex edges for atom.
                    int nepa = -1;
                    // Pointer to first edge.
                    int iep = convexEdgeFirst[ia];
                    // Check whether finished gathering.
                    while (iep > -1) {
                        // One more edge.
                        nepa++;
                        if (nepa >= maxepa) {
                            throw new EnergyException(" Too many Convex Edges for Atom", false);
                        }
                        // Store vertices of edge.
                        av[0][nepa] = convexEdgeVertexNumber[0][iep];
                        av[1][nepa] = convexEdgeVertexNumber[1][iep];
                        // Store convex edge number.
                        aep[nepa] = iep;
                        int ic = convexEdgeCircleNumber[iep];
                        // Store circle number.
                        aic[nepa] = ic;
                        // Get neighboring atom.
                        int it = circleTorusNumber[ic];
                        int ia2;
                        if (torusAtomNumber[0][it] == ia) {
                            ia2 = torusAtomNumber[1][it];
                        } else {
                            ia2 = torusAtomNumber[0][it];
                        }
                        // Store other atom number, we might need it sometime.
                        aia[nepa] = ia2;
                        /*
                          Vector from atom to circle center; also vector
                          from atom to center of neighboring atom sometimes
                          we use one vector, sometimes the other.
                         */
                        for (int k = 0; k < 3; k++) {
                            acvect[k][nepa] = circleCenter[k][ic] - atomCoords[k][ia];
                            aavect[k][nepa] = atomCoords[k][ia2] - atomCoords[k][ia];
                        }
                        // Circle radius.
                        acr[nepa] = circleRadius[ic];
                        // Pointer to next edge.
                        iep = convexEdgeNext[iep];
                    }
                    if (nepa <= -1) {
                        throw new EnergyException(" No Edges for Non-buried, Non-free Atom", false);
                    }
                    // Form cycles; initialize all the convex edges as not used in cycle.
                    for (int iepa = 0; iepa < nepa + 1; iepa++) {
                        epused[iepa] = false;
                    }
                    // Save old number of cycles.
                    int ncyold = nCycles;
                    int nused = -1;
                    int ncypa = -1;
                    while (nused < nepa) {
                        // Look for starting edge.
                        int iepa;
                        for (iepa = 0; iepa <= nepa; iepa++) {
                            if (!epused[iepa]) {
                                break;
                            }
                        }
                        if (iepa > nepa) {
                            // Cannot find starting edge; finished.
                            break;
                        }
                        // Pointer to edge.
                        iep = aep[iepa];
                        // One edge so far on this cycle.
                        int ncyep = 0;
                        // One more cycle for atom.
                        ncypa++;
                        if (ncypa >= maxcypa) {
                            throw new EnergyException(" Too many Cycles per Atom", false);
                        }
                        // Mark edge used in cycle.
                        epused[iepa] = true;
                        nused++;
                        // One more cycle for molecule.
                        nCycles++;
                        if (nCycles >= maxCycles) {
                            throw new EnergyException(" Too many Cycles", false);
                        }
                        // Index of edge in atom cycle array.
                        cyepa[ncyep][ncypa] = iepa;
                        // Store in molecule cycle array a pointer to edge.
                        convexEdgeCycleNumbers[ncyep][nCycles] = iep;
                        // Second vertex of this edge is the vertex to look for
                        // next as the first vertex of another edge.
                        int lookv = av[1][iepa];
                        // If no vertex, this cycle is finished.
                        if (lookv > -1) {
                            for (int jepa = 0; jepa < nepa + 1; jepa++) {
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
                                if (ncyep >= MAXCYEP) {
                                    throw new EnergyException(" Too many Edges per Cycle", false);
                                }
                                epused[jepa] = true;
                                nused++;
                                // Store index in local edge array.
                                cyepa[ncyep][ncypa] = jepa;
                                // Store pointer to edge.
                                convexEdgeCycleNumbers[ncyep][nCycles] = iep;
                                // New vertex to look for.
                                lookv = av[1][jepa];
                                // If no vertex, this cycle is in trouble.
                                if (lookv <= -1) {
                                    throw new EnergyException(" Pointer Error in Cycle", true);
                                }
                                jepa = 0;
                            }
                            // It better connect to first edge of cycle.
                            if (lookv != av[0][iepa]) {
                                throw new EnergyException(" Cycle does not Close", true);
                            }
                        }
                        // This cycle is finished, store number of edges in cycle.
                        ncyepa[ncypa] = ncyep;
                        convexEdgeCycleNum[nCycles] = ncyep;
                    }
                    // Compare cycles for inside/outside relation; check to
                    // see if cycle i is inside cycle j.
                    for (int icya = 0; icya <= ncypa; icya++) {
                        innerLoop:
                        for (int jcya = 0; jcya <= ncypa; jcya++) {
                            int jcy = ncyold + jcya + 1;
                            // Initialize.
                            cycy[icya][jcya] = true;
                            // Check for same cycle.
                            if (icya == jcya || ncyepa[jcya] < 2) {
                                continue;
                            }
                            // If cycles i and j have a pair of edges belonging to the same circle,
                            // then they are outside each other.
                            for (int icyep = 0; icyep <= ncyepa[icya]; icyep++) {
                                int iepa = cyepa[icyep][icya];
                                int ic = aic[iepa];
                                for (int jcyep = 0; jcyep <= ncyepa[jcya]; jcyep++) {
                                    int jepa = cyepa[jcyep][jcya];
                                    int jc = aic[jepa];
                                    if (ic == jc) {
                                        cycy[icya][jcya] = false;
                                        continue innerLoop;
                                    }
                                }
                            }
                            int iepa = cyepa[0][icya];
                            ai[0] = aavect[0][iepa];
                            ai[1] = aavect[1][iepa];
                            ai[2] = aavect[2][iepa];
                            double anaa = r(ai);
                            double factor = radius[ia] / anaa;
                            // North pole and unit vector pointer south.
                            for (int k = 0; k < 3; k++) {
                                pole[k] = factor * aavect[k][iepa] + atomCoords[k][ia];
                                unvect[k] = -aavect[k][iepa] / anaa;
                            }
                            cycy[icya][jcya] = ptincy(pole, unvect, jcy);
                        }
                    }
                    // Group cycles into faces; direct comparison for i and j.
                    for (int icya = 0; icya <= ncypa; icya++) {
                        for (int jcya = 0; jcya <= ncypa; jcya++) {
                            // Tentatively say that cycles i and j bound the
                            // same face if they are inside each other.
                            samef[icya][jcya] = (cycy[icya][jcya] && cycy[jcya][icya]);
                        }
                    }
                    // If i is in exterior of k, and k is in interior of i
                    // and j, then i and j do not bound the same face.
                    for (int icya = 0; icya <= ncypa; icya++) {
                        for (int jcya = 0; jcya <= ncypa; jcya++) {
                            if (icya != jcya) {
                                for (int kcya = 0; kcya <= ncypa; kcya++) {
                                    if (kcya != icya && kcya != jcya) {
                                        if (cycy[kcya][icya] && cycy[kcya][jcya] &&
                                                (!cycy[icya][kcya])) {
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
                                for (int kcya = jcya + 1; kcya <= ncypa; kcya++) {
                                    if (samef[jcya][kcya]) {
                                        samef[icya][kcya] = true;
                                        samef[kcya][icya] = true;
                                    }
                                }
                            }
                        }
                    }
                    // Group cycles belonging to the same face.
                    for (int icya = 0; icya <= ncypa; icya++) {
                        cyused[icya] = false;
                    }
                    // Clear number of cycles used in bounding faces.
                    nused = -1;
                    for (int icya = 0; icya <= ncypa; icya++) {
                        // Check for already used.
                        if (cyused[icya]) {
                            continue;
                        }
                        // One more convex face.
                        nConvexFaces++;
                        if (nConvexFaces >= maxConvexFaces) {
                            throw new EnergyException(" Too many Convex Faces", false);
                        }
                        // Clear number of cycles for face.
                        convexFaceNumCycles[nConvexFaces] = -1;
                        // Pointer from face to atom.
                        convexFaceAtomNumber[nConvexFaces] = ia;
                        // Look for all other cycles belonging to same face.
                        for (int jcya = 0; jcya <= ncypa; jcya++) {
                            // Check for cycle already used in another face.
                            if (cyused[jcya]) {
                                continue;

                            }
                            // Cycles i and j belonging to same face
                            if (!samef[icya][jcya]) {
                                continue;
                            }
                            // Mark cycle used.
                            cyused[jcya] = true;
                            nused++;
                            // One more cycle for face.
                            convexFaceNumCycles[nConvexFaces]++;
                            if (convexFaceNumCycles[nConvexFaces] >= MAXFPCY) {
                                throw new EnergyException(" Too many Cycles bounding Convex Face", false);
                            }
                            int i = convexFaceNumCycles[nConvexFaces];
                            // Store cycle number.
                            convexFaceCycleNumbers[i][nConvexFaces] = ncyold + jcya + 1;
                            // Check for finished.
                            if (nused >= ncypa) {
                                continue atomLoop;
                            }
                        }
                    }
                    // Should not fall though end of for loops.
                    throw new EnergyException(" Not all Cycles grouped into Convex Faces", true);
                }
                // Once face for free atom; no cycles.
                nConvexFaces++;
                if (nConvexFaces >= maxConvexFaces) {
                    throw new EnergyException(" Too many Convex Faces", false);
                }
                convexFaceAtomNumber[nConvexFaces] = ia;
                convexFaceNumCycles[nConvexFaces] = -1;
            }
        }

        /**
         * The vam method takes the analytical molecular surface defined as
         * a collection of spherical and toroidal polygons and uses it to
         * compute the volume and surface area
         */
        private void vam() {
            final int maxdot = 1000;
            final int maxop = 100;
            final int nscale = 20;
            int[] ivs = new int[3];
            int[] ispind = new int[3];
            int[] ispnd2 = new int[3];
            int[] ifnop = new int[maxop];
            int[] enfs = new int[20 * nAtoms];
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
            double[] atomArea = new double[nAtoms];
            boolean[] ate = new boolean[maxop];
            boolean[] vip = new boolean[3];
            boolean[] cinsp = new boolean[1];

            int[] nlap = null;
            int[][] fnt = null;
            int[][] nspt = null;
            double[] depths = null;
            double[] cora = null;
            double[] corv = null;
            double[][] alts = null;
            double[][] fncen = null;
            double[][][] fnvect = null;
            boolean[] badav = null;
            boolean[] badt = null;
            boolean[][] fcins = null;
            boolean[][] fcint = null;
            boolean[][] fntrev = null;
            if (nConcaveFaces >= 0) {
                nlap = new int[nConcaveFaces + 1];
                fnt = new int[3][nConcaveFaces + 1];
                nspt = new int[3][nConcaveFaces + 1];
                depths = new double[nConcaveFaces + 1];
                cora = new double[nConcaveFaces + 1];
                corv = new double[nConcaveFaces + 1];
                alts = new double[3][nConcaveFaces + 1];
                fncen = new double[3][nConcaveFaces + 1];
                fnvect = new double[3][3][nConcaveFaces + 1];
                badav = new boolean[nConcaveFaces + 1];
                badt = new boolean[nConcaveFaces + 1];
                fcins = new boolean[3][nConcaveFaces + 1];
                fcint = new boolean[3][nConcaveFaces + 1];
                fntrev = new boolean[3][nConcaveFaces + 1];
            }

            // Compute the volume of the interior polyhedron.
            double polyhedronVolume = 0.0;
            for (int i = 0; i <= nConcaveFaces; i++) {
                polyhedronVolume += measurePrism(i);
            }

            // Compute the area and volume due to convex faces as well as
            // the area partitioned among the atoms.
            double convexFaceArea = 0.0;
            double convexFaceVolume = 0.0;
            fill(atomArea, 0.0);
            double[] convexFaces = {0.0, 0.0};
            for (int i = 0; i <= nConvexFaces; i++) {
                measureConvexFace(i, convexFaces);
                int ia = convexFaceAtomNumber[i];
                atomArea[ia] += convexFaces[0];
                convexFaceArea += convexFaces[0];
                convexFaceVolume += convexFaces[1];
            }

            // Compute the area and volume due to saddle faces
            // as well as the spindle correction value.
            double saddleFaceArea = 0.0;
            double saddleFaceVolume = 0.0;
            double spindleArea = 0.0;
            double spindleVolume = 0.0;
            double[] saddle = {0.0, 0.0, 0.0, 0.0};
            for (int i = 0; i <= nSaddleFaces; i++) {
                for (int k = 0; k < 2; k++) {
                    int ien = saddleConcaveEdgeNumbers[k][i];
                    if (ien > -1) {
                        enfs[ien] = i;
                    }
                }
                measureSaddleFace(i, saddle);
                double areas = saddle[0];
                double vols = saddle[1];
                double areasp = saddle[2];
                double volsp = saddle[3];
                saddleFaceArea += areas;
                saddleFaceVolume += vols;
                spindleArea += areasp;
                spindleVolume += volsp;
                if (areas - areasp < 0.0) {
                    throw new EnergyException(" Negative Area for Saddle Face", true);
                }
            }

            // Compute the area and volume due to concave faces.
            double concaveFaceArea = 0.0;
            double concaveFaceVolume = 0.0;
            double[] concaveFaces = {0.0, 0.0};
            for (int i = 0; i <= nConcaveFaces; i++) {
                measureConcaveFace(i, concaveFaces);
                double arean = concaveFaces[0];
                double voln = concaveFaces[1];
                concaveFaceArea += arean;
                concaveFaceVolume += voln;
            }

            // Compute the area and volume lens correction values.
            double lensArea = 0.0;
            double lensAreaNumeric = 0.0;
            double lensVolume = 0.0;
            double lensVolumeNumeric = 0.0;

            if (probe > 0.0) {
                int ndots = genDots(maxdot, dots, probe);
                double dota = (4.0 * PI * probe * probe) / ndots;
                for (int ifn = 0; ifn <= nConcaveFaces; ifn++) {
                    nlap[ifn] = -1;
                    cora[ifn] = 0.0;
                    corv[ifn] = 0.0;
                    badav[ifn] = false;
                    badt[ifn] = false;
                    for (int k = 0; k < 3; k++) {
                        nspt[k][ifn] = -1;
                    }
                    int ien = concaveFaceEdgeNumbers[0][ifn];
                    int iv = concaveEdgeVertexNumbers[0][ien];
                    int ip = vertexProbeNumber[iv];
                    getVector(ai, alts, ifn);
                    depths[ifn] = depth(ip, ai);
                    for (int k = 0; k < 3; k++) {
                        fncen[k][ifn] = probeCoords[k][ip];
                    }
                    // Get vertices and vectors.
                    for (int ke = 0; ke < 3; ke++) {
                        ien = concaveFaceEdgeNumbers[ke][ifn];
                        ivs[ke] = concaveEdgeVertexNumbers[0][ien];
                        int ia = vertexAtomNumbers[ivs[ke]];
                        int ifs = enfs[ien];
                        int iep = saddleConvexEdgeNumbers[0][ifs];
                        int ic = convexEdgeCircleNumber[iep];
                        int it = circleTorusNumber[ic];
                        fnt[ke][ifn] = it;
                        fntrev[ke][ifn] = (torusAtomNumber[0][it] != ia);
                    }
                    for (int ke = 0; ke < 3; ke++) {
                        for (int k = 0; k < 3; k++) {
                            vects[k][ke] = vertexCoords[k][ivs[ke]] - probeCoords[k][ip];
                        }
                    }
                    // Calculate normal vectors for the three planes that
                    // cut out the geodesic triangle.
                    getVector(ai, vects, 0);
                    getVector(aj, vects, 1);
                    cross(ai, aj, ak);
                    norm(ak, ak);
                    putVector(ak, fnvect, 0, ifn);
                    getVector(ak, vects, 2);
                    cross(aj, ak, ai);
                    norm(ai, ai);
                    putVector(ai, fnvect, 1, ifn);
                    getVector(ai, vects, 0);
                    cross(ak, ai, aj);
                    norm(aj, aj);
                    putVector(aj, fnvect, 2, ifn);
                }
                for (int ifn = 0; ifn <= nConcaveFaces - 1; ifn++) {
                    for (int jfn = ifn + 1; jfn <= nConcaveFaces; jfn++) {
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
                        double rat = check(dpp / (2.0 * probe));
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
                            boolean cintp = circlePlane(ppm, rm, upp, ai, aj, cinsp, xpnt1, xpnt2);
                            fcins[ke][ifn] = cinsp[0];
                            fcint[ke][ifn] = cintp;
                            if (!cinsp[0]) {
                                alli = false;
                            }
                            if (cintp) {
                                anyi = true;
                            }
                            if (!cintp) {
                                continue;
                            }
                            int it = fnt[ke][ifn];
                            if (torusRadius[it] > probe) {
                                continue;
                            }
                            for (int ke2 = 0; ke2 < 3; ke2++) {
                                if (it == fnt[ke2][jfn]) {
                                    ispind[ke] = it;
                                    nspt[ke][ifn]++;
                                    ispnd2[ke2] = it;
                                    nspt[ke2][jfn]++;
                                    spindl = true;
                                }
                            }
                            if (ispind[ke] == -1) {
                                continue;
                            }

                            // Check that the two ways of calculating intersection points match.
                            rat = check(torusRadius[it] / probe);
                            thetaq[ke] = acos(rat);
                            double stq = sin(thetaq[ke]);
                            if (fntrev[ke][ifn]) {
                                for (int k = 0; k < 3; k++) {
                                    uij[k] = -torusAxis[k][it];
                                }
                            } else {
                                for (int k = 0; k < 3; k++) {
                                    uij[k] = torusAxis[k][it];
                                }
                            }
                            for (int k = 0; k < 3; k++) {
                                qij[k] = torusCenter[k][it] - stq * probe * uij[k];
                                qji[k] = torusCenter[k][it] + stq * probe * uij[k];
                            }
                            for (int k = 0; k < 3; k++) {
                                umq[k] = (qij[k] - ppm[k]) / rm;
                                upq[k] = (qij[k] - fncen[k][ifn]) / probe;
                            }
                            cross(uij, upp, vect1);
                            double dt = check(dot(umq, vect1));
                            sigmaq[ke] = acos(dt);
                            getVector(ai, fnvect, ke, ifn);
                            cross(upq, ai, vect1);
                            norm(vect1, uc);
                            cross(upp, upq, vect1);
                            norm(vect1, uq);
                            dt = check(dot(uc, uq));
                            tau[ke] = PI - acos(dt);
                        }
                        boolean allj = true;
                        boolean anyj = false;
                        for (int ke = 0; ke < 3; ke++) {
                            getVector(ai, fncen, jfn);
                            getVector(aj, fnvect, ke, jfn);
                            boolean cintp = circlePlane(ppm, rm, upp, ai, aj, cinsp, xpnt1, xpnt2);
                            fcins[ke][jfn] = cinsp[0];
                            fcint[ke][jfn] = cintp;
                            if (!cinsp[0]) {
                                allj = false;
                            }
                            if (cintp) {
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
                            int ien = concaveFaceEdgeNumbers[ke][ifn];
                            int iv1 = concaveEdgeVertexNumbers[0][ien];
                            int iv2 = concaveEdgeVertexNumbers[1][ien];
                            for (int k = 0; k < 3; k++) {
                                vect3[k] = vertexCoords[k][iv1] - fncen[k][ifn];
                                vect4[k] = vertexCoords[k][iv2] - fncen[k][ifn];
                            }
                            for (int ke2 = 0; ke2 < 3; ke2++) {
                                if (ispind[ke] == ispnd2[ke2] || ispind[ke] == -1) {
                                    continue;
                                }
                                getVector(ai, fncen, ifn);
                                getVector(aj, fnvect, ke, ifn);
                                getVector(ak, fncen, jfn);
                                getVector(al, fnvect, ke2, jfn);
                                boolean cintp = circlePlane(ai, probe, aj, ak, al, cinsp, xpnt1, xpnt2);
                                if (!cintp) {
                                    continue;
                                }
                                ien = concaveFaceEdgeNumbers[ke2][jfn];
                                iv1 = concaveEdgeVertexNumbers[0][ien];
                                iv2 = concaveEdgeVertexNumbers[1][ien];
                                for (int k = 0; k < 3; k++) {
                                    vect7[k] = vertexCoords[k][iv1] - fncen[k][jfn];
                                    vect8[k] = vertexCoords[k][iv2] - fncen[k][jfn];
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
                                if (tripleProduct(vect3, vect1, ai) < 0.0
                                        || tripleProduct(vect1, vect4, ai) < 0.0
                                        || tripleProduct(vect7, vect5, aj) < 0.0
                                        || tripleProduct(vect5, vect8, aj) < 0.0) {
                                    if (!(tripleProduct(vect3, vect2, ai) < 0.0
                                            || tripleProduct(vect2, vect4, ai) < 0.0
                                            || tripleProduct(vect7, vect6, aj) < 0.0
                                            || tripleProduct(vect6, vect8, aj) < 0.0)) {
                                        badav[ifn] = true;
                                    }
                                } else {
                                    badav[ifn] = true;
                                }
                            }
                        }
                        for (int ke = 0; ke < 3; ke++) {
                            int ien = concaveFaceEdgeNumbers[ke][ifn];
                            int iv1 = concaveEdgeVertexNumbers[0][ien];
                            int iv2 = concaveEdgeVertexNumbers[1][ien];
                            for (int k = 0; k < 3; k++) {
                                vect3[k] = vertexCoords[k][iv1] - fncen[k][ifn];
                                vect4[k] = vertexCoords[k][iv2] - fncen[k][ifn];
                            }
                            for (int ke2 = 0; ke2 < 3; ke2++) {
                                if (ispind[ke] == ispnd2[ke2] || ispnd2[ke2] == -1) {
                                    continue;
                                }
                                getVector(ai, fncen, ifn);
                                getVector(aj, fnvect, ke, ifn);
                                getVector(ak, fncen, jfn);
                                getVector(al, fnvect, ke2, jfn);
                                boolean cintp = circlePlane(ak, probe, al, ai, aj, cinsp, xpnt1, xpnt2);
                                if (!cintp) {
                                    continue;
                                }
                                ien = concaveFaceEdgeNumbers[ke2][jfn];
                                iv1 = concaveEdgeVertexNumbers[0][ien];
                                iv2 = concaveEdgeVertexNumbers[1][ien];
                                for (int k = 0; k < 3; k++) {
                                    vect7[k] = vertexCoords[k][iv1] - fncen[k][jfn];
                                    vect8[k] = vertexCoords[k][iv2] - fncen[k][jfn];
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
                                if (tripleProduct(vect3, vect1, ai) < 0.0
                                        || tripleProduct(vect1, vect4, ai) < 0.0
                                        || tripleProduct(vect7, vect5, aj) < 0.0
                                        || tripleProduct(vect5, vect8, aj) < 0.0) {
                                    if (!(tripleProduct(vect3, vect2, ai) < 0.0
                                            || tripleProduct(vect2, vect4, ai) < 0.0
                                            || tripleProduct(vect7, vect6, aj) < 0.0
                                            || tripleProduct(vect6, vect8, aj) < 0.0)) {
                                        badav[jfn] = true;
                                    }
                                } else {
                                    badav[jfn] = true;
                                }
                            }
                        }

                        double sumlam = 0.0;
                        double sumsig = 0.0;
                        double sumsc = 0.0;
                        for (int ke = 0; ke < 3; ke++) {
                            if (ispind[ke] != -1) {
                                sumlam += (PI - tau[ke]);
                                sumsig += (sigmaq[ke] - PI);
                                sumsc += (sin(sigmaq[ke]) * cos(sigmaq[ke]));
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

                        // Check for vertex on opposing probe in face.
                        outerloop:
                        for (int kv = 0; kv < 3; kv++) {
                            vip[kv] = false;
                            int ien = concaveFaceEdgeNumbers[kv][jfn];
                            int iv = concaveEdgeVertexNumbers[0][ien];
                            for (int k = 0; k < 3; k++) {
                                vect1[k] = vertexCoords[k][iv] - fncen[k][ifn];
                            }
                            norm(vect1, vect1);
                            for (int ke = 0; ke < 3; ke++) {
                                getVector(ai, fnvect, ke, ifn);
                                getVector(aj, vertexCoords, iv);
                                double dt = dot(ai, aj);
                                if (dt > 0.0) {
                                    continue outerloop;
                                }
                                vip[kv] = true;
                            }
                        }
                    }
                }

                outerLoop:
                for (int ifn = 0; ifn <= nConcaveFaces; ifn++) {
                    for (int ke = 0; ke < 3; ke++) {
                        if (nspt[ke][ifn] > 0) {
                            badt[ifn] = true;
                            break outerLoop;
                        }
                    }
                }

                double fourProbe2 = 4.0 * probe * probe;
                for (int ifn = 0; ifn <= nConcaveFaces; ifn++) {
                    if (nlap[ifn] <= -1) {
                        continue;
                    }
                    // Gather all overlapping probes.
                    int nop = -1;
                    for (int jfn = 0; jfn <= nConcaveFaces; jfn++) {
                        if (ifn != jfn) {
                            getVector(ai, fncen, ifn);
                            getVector(aj, fncen, jfn);
                            double dij2 = dist2(ai, aj);
                            if (dij2 <= fourProbe2 && depths[jfn] <= probe) {
                                nop++;
                                if (nop >= maxop) {
                                    throw new EnergyException(" NOP Overflow in VAM", false);
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
                        double rsc = (isc + 1) - 0.5;
                        dotv[isc] = probe * dota * rsc * rsc * scinc * scinc * scinc;
                    }
                    for (int iop = 0; iop <= nop; iop++) {
                        ate[iop] = false;
                    }
                    int neatmx = 0;
                    dotLoop:
                    for (int idot = 0; idot < ndots; idot++) {
                        for (int ke = 0; ke < 3; ke++) {
                            getVector(ai, fnvect, ke, ifn);
                            getVector(aj, dots, idot);
                            double dt = dot(ai, aj);
                            if (dt > 0.0) {
                                continue dotLoop;
                            }
                        }
                        for (int k = 0; k < 3; k++) {
                            tdots[k][idot] = fncen[k][ifn] + dots[k][idot];
                        }
                        for (int iop = 0; iop <= nop; iop++) {
                            int jfn = ifnop[iop];
                            getVector(ai, tdots, idot);
                            getVector(aj, fncen, jfn);
                            double ds2 = dist2(ai, aj);
                            if (ds2 < probe * probe) {
                                areado += dota;
                                break;
                            }
                        }
                        for (int isc = 0; isc < nscale; isc++) {
                            double rsc = (isc + 1) - 0.5;
                            for (int k = 0; k < 3; k++) {
                                sdot[k] = fncen[k][ifn] + rsc * scinc * dots[k][idot];
                            }
                            int neat = 0;
                            optLoop:
                            for (int iop = 0; iop <= nop; iop++) {
                                int jfn = ifnop[iop];
                                getVector(ai, fncen, jfn);
                                double ds2 = dist2(sdot, ai);
                                if (ds2 < probe * probe) {
                                    for (int k = 0; k < 3; k++) {
                                        vect1[k] = sdot[k] - fncen[k][jfn];
                                    }
                                    for (int ke = 0; ke < 3; ke++) {
                                        getVector(ai, fnvect, ke, jfn);
                                        double dt = dot(ai, vect1);
                                        if (dt > 0.0) {
                                            continue optLoop;
                                        }
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
                    for (int iop = 0; iop <= nop; iop++) {
                        if (ate[iop]) {
                            nate++;
                        }
                    }

                    // Use either the analytical or numerical correction.
                    boolean usenum = (nate > nlap[ifn] + 1 || neatmx > 1 || badt[ifn]);
                    if (usenum) {
                        cora[ifn] = coran;
                        corv[ifn] = corvn;
                        lensAreaNumeric += cora[ifn];
                        lensVolumeNumeric += corv[ifn];
                    } else if (badav[ifn]) {
                        corv[ifn] = corvn;
                        lensVolumeNumeric += corv[ifn];
                    }
                    lensArea += cora[ifn];
                    lensVolume += corv[ifn];
                }
            }

            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(" Convex Faces:     %6d Area: %12.6f Volume: %12.6f", nConvexFaces + 1, convexFaceArea, convexFaceVolume));
                logger.fine(format(" Saddle Faces:     %6d Area: %12.6f Volume: %12.6f", nSaddleFaces + 1, saddleFaceArea, saddleFaceVolume));
                logger.fine(format(" Concave Faces:    %6d Area: %12.6f Volume: %12.6f", nConcaveFaces + 1, concaveFaceArea, concaveFaceVolume));
                logger.fine(format(" Buried Polyhedron:                          Volume: %12.6f", polyhedronVolume));
                if (probe > 0.0) {
                    logger.fine(format("\n Spindle Correction:            Area: %12.6f Volume: %12.6f", -spindleArea, -spindleVolume));
                    logger.fine(format(" Lens Analytic Correction:      Area: %12.6f Volume: %12.6f",
                            -lensArea - lensAreaNumeric, lensVolume - lensVolumeNumeric));
                    logger.fine(format(" Lens Numerical Correction:     Area: %12.6f Volume: %12.6f",
                            lensAreaNumeric, lensVolumeNumeric));
                }
            }

            // Finally, compute the total area and total volume.
            double area = convexFaceArea + saddleFaceArea + concaveFaceArea - spindleArea - lensArea;
            double volume = convexFaceVolume + saddleFaceVolume + concaveFaceVolume + polyhedronVolume - spindleVolume + lensVolume;

            localVolume += volume;
            localSurfaceArea += area;
        }

        /**
         * The triple method finds the triple product of three vectors; used
         * as a service routine by the Connolly surface area and volume
         * computation.
         *
         * @param x Vector 1.
         * @param y Vector 2.
         * @param z Vector 3.
         * @return The triple product.
         */
        private double tripleProduct(double[] x, double[] y, double[] z) {
            double[] xy = new double[3];
            cross(x, y, xy);
            return dot(xy, z);
        }

        /**
         * The insertConvexEdgeForAtom method inserts a convex edge into the linked list for an atom.
         *
         * @param edgeNumber The edge number.
         * @param atomNumber The atom number.
         */
        private void insertConvexEdgeForAtom(int edgeNumber, int atomNumber) {
            if (edgeNumber <= -1) {
                throw new EnergyException(" Bad Edge Number in IPEDGE", true);
            }
            if (atomNumber <= -1) {
                throw new EnergyException(" Bad Atom Number in IPEDGE", true);
            }
            // Set beginning of list or add to end.
            if (convexEdgeFirst[atomNumber] == -1) {
                convexEdgeFirst[atomNumber] = edgeNumber;
                convexEdgeNext[edgeNumber] = -1;
                convexEdgeLast[atomNumber] = edgeNumber;
            } else {
                convexEdgeNext[convexEdgeLast[atomNumber]] = edgeNumber;
                convexEdgeNext[edgeNumber] = -1;
                convexEdgeLast[atomNumber] = edgeNumber;
            }
        }

        /**
         * Check for eaten by neighbor.
         *
         * @param northPole  North pole vector.
         * @param unitVector Unit vector.
         * @param icy        Cycle index.
         * @return True if the angle sum is greater than 0.
         */
        private boolean ptincy(double[] northPole, double[] unitVector, int icy) {
            double[] acvect = new double[3];
            double[] cpvect = new double[3];
            double[][] projectedVertsForConvexEdges = new double[3][MAXCYEP];
            double[][] unitVectorsAlongEdges = new double[3][MAXCYEP];
            // Check for being eaten by neighbor.
            int iatom = -1;
            for (int ke = 0; ke <= convexEdgeCycleNum[icy]; ke++) {
                int iep = convexEdgeCycleNumbers[ke][icy];
                int ic = convexEdgeCircleNumber[iep];
                int it = circleTorusNumber[ic];
                iatom = circleAtomNumber[ic];
                int iaoth;
                if (torusAtomNumber[0][it] == iatom) {
                    iaoth = torusAtomNumber[1][it];
                } else {
                    iaoth = torusAtomNumber[0][it];
                }
                for (int k = 0; k < 3; k++) {
                    acvect[k] = atomCoords[k][iaoth] - atomCoords[k][iatom];
                    cpvect[k] = northPole[k] - circleCenter[k][ic];
                }
                if (dot(acvect, cpvect) >= 0.0) {
                    return false;
                }
            }
            if (convexEdgeCycleNum[icy] <= 1) {
                return true;
            }
            if (projectedVertexForEachConvexEdge(northPole, unitVector, icy, iatom, projectedVertsForConvexEdges)) {
                return true;
            }
            int nEdges = convexEdgeCycleNum[icy];
            unitVectorsAlongEdges(projectedVertsForConvexEdges, nEdges, unitVectorsAlongEdges);
            double angleSum = sumOfAnglesAtVerticesForCycle(unitVectorsAlongEdges, nEdges, unitVector);
            return (angleSum > 0.0);
        }

        /**
         * Projected vertex for each convex edge.
         *
         * @param northPole                    North pole vector.
         * @param unitVector                   Unit vector.
         * @param icy                          Cycle number.
         * @param ia                           Atom number.
         * @param projectedVertsForConvexEdges Projected vertex for each convex edge.
         * @return Returns true if the projection terminates before completion.
         */
        private boolean projectedVertexForEachConvexEdge(double[] northPole, double[] unitVector,
                                                         int icy, int ia, double[][] projectedVertsForConvexEdges) {
            double[] polev = new double[3];
            for (int ke = 0; ke <= convexEdgeCycleNum[icy]; ke++) {
                // Vertex number (use first vertex of edge).
                int iep = convexEdgeCycleNumbers[ke][icy];
                int iv = convexEdgeVertexNumber[0][iep];
                if (iv != -1) {
                    // Vector from north pole to vertex.
                    for (int k = 0; k < 3; k++) {
                        polev[k] = vertexCoords[k][iv] - northPole[k];
                    }
                    // Calculate multiplication factor.
                    double dt = dot(polev, unitVector);
                    if (dt == 0.0) {
                        return true;
                    }
                    double f = (2 * radius[ia]) / dt;
                    if (f < 1.0) {
                        return true;
                    }
                    // Projected vertex for this convex edge.
                    for (int k = 0; k < 3; k++) {
                        projectedVertsForConvexEdges[k][ke] = northPole[k] + f * polev[k];
                    }
                }
            }
            return false;
        }

        /**
         * Sum of angles at vertices of cycle.
         *
         * @param unitVectorsAlongEdges Input vertices.
         * @param nEdges                Number of edges.
         * @param unitVector            Input vector.
         * @return Sum of angles at vertices of cycle
         */
        private double sumOfAnglesAtVerticesForCycle(double[][] unitVectorsAlongEdges, int nEdges, double[] unitVector) {
            double[] crs = new double[3];
            double[] ai = new double[3];
            double[] aj = new double[3];
            double sumOfAngles = 0.0;
            // Sum of angles at vertices of cycle.
            for (int ke = 0; ke <= nEdges; ke++) {
                double dt;
                if (ke < nEdges) {
                    getVector(ai, unitVectorsAlongEdges, ke);
                    getVector(aj, unitVectorsAlongEdges, ke + 1);
                    dt = dot(ai, aj);
                    cross(ai, aj, crs);
                } else {
                    // Closing edge of cycle.
                    getVector(ai, unitVectorsAlongEdges, ke);
                    getVector(aj, unitVectorsAlongEdges, 0);
                    dt = dot(ai, aj);
                    cross(ai, aj, crs);
                }
                dt = check(dt);
                double ang = acos(dt);
                if (dot(crs, unitVector) > 0.0) {
                    ang = -ang;
                }
                // Add to total for cycle.
                sumOfAngles += ang;
            }
            return sumOfAngles;
        }

        /**
         * Calculate unit vectors along edges.
         *
         * @param projectedVertsForConvexEdges Projected verts for convex edges.
         * @param nEdges                       Number of edges.
         * @param unitVectorsAlongEdges        Unit vectors along edges.
         */
        private void unitVectorsAlongEdges(double[][] projectedVertsForConvexEdges, int nEdges, double[][] unitVectorsAlongEdges) {
            double[] ai = new double[3];
            // Calculate unit vectors along edges.
            for (int ke = 0; ke <= nEdges; ke++) {
                // Get index of second edge of corner.
                int ke2;
                if (ke < nEdges) {
                    ke2 = ke + 1;
                } else {
                    ke2 = 0;
                }
                // Unit vector along edge of cycle.
                for (int k = 0; k < 3; k++) {
                    unitVectorsAlongEdges[k][ke] = projectedVertsForConvexEdges[k][ke2] - projectedVertsForConvexEdges[k][ke];
                }
                getVector(ai, unitVectorsAlongEdges, ke);
                double epun = r(ai);
                if (epun <= 0.0) {
                    throw new EnergyException(" Null Edge in Cycle", true);
                }
                // Normalize.
                if (epun > 0.0) {
                    for (int k = 0; k < 3; k++) {
                        unitVectorsAlongEdges[k][ke] = unitVectorsAlongEdges[k][ke] / epun;
                    }
                } else {
                    for (int k = 0; k < 3; k++) {
                        unitVectorsAlongEdges[k][ke] = 0.0;
                    }
                }
            }
            // Vectors for null edges come from following or preceding edges.
            for (int ke = 0; ke <= nEdges; ke++) {
                getVector(ai, unitVectorsAlongEdges, ke);
                if (r(ai) <= 0.0) {
                    int le = ke - 1;
                    if (le < 0) {
                        le = nEdges;
                    }
                    for (int k = 0; k < 3; k++) {
                        unitVectorsAlongEdges[k][ke] = unitVectorsAlongEdges[k][le];
                    }
                }
            }
        }

        /**
         * The measpm method computes the volume of a single prism section
         * of the full interior polyhedron.
         */
        private double measurePrism(int ifn) {
            double[][] pav = new double[3][3];
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] vect3 = new double[3];
            double height = 0.0;
            for (int ke = 0; ke < 3; ke++) {
                int ien = concaveFaceEdgeNumbers[ke][ifn];
                int iv = concaveEdgeVertexNumbers[0][ien];
                int ia = vertexAtomNumbers[iv];
                height += atomCoords[2][ia];
                int ip = vertexProbeNumber[iv];
                for (int k = 0; k < 3; k++) {
                    pav[k][ke] = atomCoords[k][ia] - probeCoords[k][ip];
                }
            }
            height *= 1.0 / 3.0;
            for (int k = 0; k < 3; k++) {
                vect1[k] = pav[k][1] - pav[k][0];
                vect2[k] = pav[k][2] - pav[k][0];
            }
            cross(vect1, vect2, vect3);
            return height * vect3[2] / 2.0;
        }

        /**
         * Compute the area and volume of a convex face.
         *
         * @param iConvexFace The convex face index.
         * @param areaVolume  The measured area and volume.
         */
        private void measureConvexFace(int iConvexFace, double[] areaVolume) {
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
            int ia = convexFaceAtomNumber[iConvexFace];
            int ncycle = convexFaceNumCycles[iConvexFace];
            int ieuler;
            if (ncycle > -1) {
                ieuler = 1 - ncycle;
            } else {
                ieuler = 2;
            }
            for (int icyptr = 0; icyptr < ncycle + 1; icyptr++) {
                int icy = convexFaceCycleNumbers[icyptr][iConvexFace];
                int nedge = convexEdgeCycleNum[icy];
                for (int ke = 0; ke < nedge + 1; ke++) {
                    int iep = convexEdgeCycleNumbers[ke][icy];
                    int ic = convexEdgeCircleNumber[iep];
                    int it = circleTorusNumber[ic];
                    int ia2;
                    if (ia == torusAtomNumber[0][it]) {
                        ia2 = torusAtomNumber[1][it];
                    } else {
                        ia2 = torusAtomNumber[0][it];
                    }
                    for (int k = 0; k < 3; k++) {
                        acvect[k] = circleCenter[k][ic] - atomCoords[k][ia];
                        aavect[k] = atomCoords[k][ia2] - atomCoords[k][ia];
                    }
                    norm(aavect, aavect);
                    double dt = dot(acvect, aavect);
                    double geo = -dt / (radius[ia] * circleRadius[ic]);
                    int iv1 = convexEdgeVertexNumber[0][iep];
                    int iv2 = convexEdgeVertexNumber[1][iep];
                    if (iv1 == -1 || iv2 == -1) {
                        angle = 2.0 * PI;
                    } else {
                        for (int k = 0; k < 3; k++) {
                            vect1[k] = vertexCoords[k][iv1] - circleCenter[k][ic];
                            vect2[k] = vertexCoords[k][iv2] - circleCenter[k][ic];
                            radial[k][ke] = vertexCoords[k][iv1] - atomCoords[k][ia];
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
                        angle = vectorAngle(vect1, vect2, aavect, -1.0);
                    }
                    gcurve += circleRadius[ic] * angle * geo;
                    if (nedge != 0) {
                        if (ke > 0) {
                            getVector(ai, tanv, 1, ke - 1);
                            getVector(aj, tanv, 0, ke);
                            getVector(ak, radial, ke);
                            angle = vectorAngle(ai, aj, ak, 1.0);
                            if (angle < 0.0) {
                                throw new EnergyException(" Negative Angle in measureConvexFace", true);
                            }
                            pcurve += angle;
                        }
                    }
                }
                if (nedge > 0) {
                    getVector(ai, tanv, 1, nedge);
                    getVector(aj, tanv, 0, 0);
                    getVector(ak, radial, 0);
                    angle = vectorAngle(ai, aj, ak, 1.0);
                    if (angle < 0.0) {
                        throw new EnergyException(" Negative Angle in measureConvexFace", true);
                    }
                    pcurve += angle;
                }
            }
            double gauss = 2.0 * PI * ieuler - pcurve - gcurve;
            double areap = gauss * (radius[ia] * radius[ia]);
            double volp = areap * radius[ia] / 3.0;
            areaVolume[0] = areap;
            areaVolume[1] = volp;
        }

        /**
         * Compute the area and volume of a saddle face.
         *
         * @param iSaddleFace The saddle face index.
         * @param areaVolume  The measured area and volume.
         */
        private void measureSaddleFace(int iSaddleFace, double[] areaVolume) {
            double areas = 0.0;
            double vols = 0.0;
            double areasp = 0.0;
            double volsp = 0.0;
            int iep = saddleConvexEdgeNumbers[0][iSaddleFace];
            int ic = convexEdgeCircleNumber[iep];
            int it = circleTorusNumber[ic];
            int ia1 = torusAtomNumber[0][it];
            int ia2 = torusAtomNumber[1][it];
            double[] aavect = new double[3];
            for (int k = 0; k < 3; k++) {
                aavect[k] = atomCoords[k][ia2] - atomCoords[k][ia1];
            }
            norm(aavect, aavect);
            int iv1 = convexEdgeVertexNumber[0][iep];
            int iv2 = convexEdgeVertexNumber[1][iep];
            double phi;
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            if (iv1 == -1 || iv2 == -1) {
                phi = 2.0 * PI;
            } else {
                for (int k = 0; k < 3; k++) {
                    vect1[k] = vertexCoords[k][iv1] - circleCenter[k][ic];
                    vect2[k] = vertexCoords[k][iv2] - circleCenter[k][ic];
                }
                phi = vectorAngle(vect1, vect2, aavect, 1.0);
            }
            for (int k = 0; k < 3; k++) {
                vect1[k] = atomCoords[k][ia1] - torusCenter[k][it];
                vect2[k] = atomCoords[k][ia2] - torusCenter[k][it];
            }
            double d1 = -1.0 * dot(vect1, aavect);
            double d2 = dot(vect2, aavect);
            double theta1 = atan2(d1, torusRadius[it]);
            double theta2 = atan2(d2, torusRadius[it]);

            // Check for cusps.
            double thetaq;
            boolean cusp;
            if (torusRadius[it] < probe && theta1 > 0.0 && theta2 > 0.0) {
                cusp = true;
                double rat = torusRadius[it] / probe;
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
            double term1 = torusRadius[it] * probe * (theta1 + theta2);
            double term2 = (probe * probe) * (sin(theta1) + sin(theta2));
            areas = phi * (term1 - term2);
            if (cusp) {
                double spin = torusRadius[it] * probe * thetaq - probe * probe * sin(thetaq);
                areasp = 2.0 * phi * spin;
            }

            iep = saddleConvexEdgeNumbers[0][iSaddleFace];
            int ic2 = convexEdgeCircleNumber[iep];
            iep = saddleConvexEdgeNumbers[1][iSaddleFace];
            int ic1 = convexEdgeCircleNumber[iep];
            if (circleAtomNumber[ic1] != ia1) {
                throw new EnergyException(" Atom number inconsistency in measureSaddleFace", true);
            }
            for (int k = 0; k < 3; k++) {
                vect1[k] = circleCenter[k][ic1] - atomCoords[k][ia1];
                vect2[k] = circleCenter[k][ic2] - atomCoords[k][ia2];
            }
            double w1 = dot(vect1, aavect);
            double w2 = -1.0 * dot(vect2, aavect);
            double cone1 = phi * ((w1 * circleRadius[ic1] * circleRadius[ic1])) / 6.0;
            double cone2 = phi * ((w2 * circleRadius[ic2] * circleRadius[ic2])) / 6.0;
            term1 = (torusRadius[it] * torusRadius[it]) * probe * (sin(theta1) + sin(theta2));
            term2 = sin(theta1) * cos(theta1) + theta1 + sin(theta2) * cos(theta2) + theta2;
            term2 = torusRadius[it] * (probe * probe) * term2;
            double term3 = sin(theta1) * cos(theta1) * cos(theta1)
                    + 2.0 * sin(theta1) + sin(theta2) * cos(theta2) * cos(theta2)
                    + 2.0 * sin(theta2);
            term3 = (probe * probe * probe / 3.0) * term3;
            double volt = (phi / 2.0) * (term1 - term2 + term3);
            vols = volt + cone1 + cone2;
            if (cusp) {
                term1 = (torusRadius[it] * torusRadius[it]) * probe * sin(thetaq);
                term2 = sin(thetaq) * cos(thetaq) + thetaq;
                term2 = torusRadius[it] * (probe * probe) * term2;
                term3 = sin(thetaq) * cos(thetaq) * cos(thetaq) + 2.0 * sin(thetaq);
                term3 = (probe * probe * probe / 3.0) * term3;
                volsp = phi * (term1 - term2 + term3);
            }
            areaVolume[0] = areas;
            areaVolume[1] = vols;
            areaVolume[2] = areasp;
            areaVolume[3] = volsp;
        }

        /**
         * Compute the area and volume of a concave face.
         *
         * @param iConcaveFace The concave face index.
         * @param areaVolume   The measured area and volume.
         */
        private void measureConcaveFace(int iConcaveFace, double[] areaVolume) {
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
                int ien = concaveFaceEdgeNumbers[ke][iConcaveFace];
                int iv = concaveEdgeVertexNumbers[0][ien];
                int ia = vertexAtomNumbers[iv];
                int ip = vertexProbeNumber[iv];
                for (int k = 0; k < 3; k++) {
                    pvv[k][ke] = vertexCoords[k][iv] - probeCoords[k][ip];
                    pav[k][ke] = atomCoords[k][ia] - probeCoords[k][ip];
                }
                if (probe > 0.0) {
                    getVector(ai, pvv, ke);
                    norm(ai, ai);
                    putVector(ai, pvv, ke);
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
                    cross(ai, aj, ak);
                    norm(ak, ak);
                    putVector(ak, planev, ke);
                }
                for (int ke = 0; ke < 3; ke++) {
                    int je = ke - 1;
                    if (je < 0) {
                        je = 2;
                    }
                    getVector(ai, planev, je);
                    getVector(aj, planev, ke);
                    getVector(ak, pvv, ke);
                    angle[ke] = vectorAngle(ai, aj, ak, -1.0);
                    if (angle[ke] < 0.0) {
                        throw new EnergyException(" Negative angle in measureConcaveFace", true);
                    }
                }
                double defect = 2.0 * PI - (angle[0] + angle[1] + angle[2]);
                arean = (probe * probe) * defect;
            }
            getVector(ai, pav, 0);
            getVector(aj, pav, 1);
            getVector(ak, pav, 2);
            double simplx = -tripleProduct(ai, aj, ak) / 6.0;
            voln = simplx - arean * probe / 3.0;
            areaVolume[0] = arean;
            areaVolume[1] = voln;
        }

        /**
         * The gendot method finds the coordinates of a specified number of
         * surface points for a sphere with the input radius.
         *
         * @param ndots  Number of dots to generate.
         * @param dots   Coordinates of generaged dots.
         * @param radius Sphere radius.
         */
        private int genDots(int ndots, double[][] dots, double radius) {
            int nequat = (int) sqrt(PI * ((double) ndots));
            int nvert = (int) (0.5 * (double) nequat);
            if (nvert < 1) {
                nvert = 1;
            }
            int k = -1;
            outerloop:
            for (int i = 0; i <= nvert; i++) {
                double fi = (PI * ((double) i)) / ((double) nvert);
                double z = cos(fi);
                double xy = sin(fi);
                int nhoriz = (int) (nequat * xy);
                if (nhoriz < 1) {
                    nhoriz = 1;
                }
                for (int j = 0; j < nhoriz; j++) {
                    double fj = (2.0 * PI * ((double) j)) / ((double) nhoriz);
                    double x = cos(fj) * xy;
                    double y = sin(fj) * xy;
                    k++;
                    dots[0][k] = x * radius;
                    dots[1][k] = y * radius;
                    dots[2][k] = z * radius;
                    if (k >= ndots) {
                        break outerloop;
                    }
                }
            }
            return k;
        }

        /**
         * The cirpln method determines the points of intersection between a
         * specified circle and plane.
         */
        private boolean circlePlane(double[] circen, double cirrad, double[] cirvec,
                                    double[] plncen, double[] plnvec, boolean[] cinsp,
                                    double[] xpnt1, double[] xpnt2) {
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
         * respect to a coordinate axis.
         *
         * @param v1   First vector.
         * @param v2   Second vector.
         * @param axis Axis.
         * @param hand Hand.
         * @return An angle in the range [0, 2*PI].
         */
        private double vectorAngle(double[] v1, double[] v2, double[] axis, double hand) {
            double a1 = r(v1);
            double a2 = r(v2);
            double dt = dot(v1, v2);
            double a12 = a1 * a2;
            if (abs(a12) != 0.0) {
                dt = dt / a12;
            }
            dt = check(dt);
            double angle = acos(dt);
            if (hand * tripleProduct(v1, v2, axis) < 0.0) {
                return 2.0 * PI - angle;
            } else {
                return angle;
            }
        }

        /**
         * Compute the probe depth.
         *
         * @param ip  The probe number.
         * @param alt The depth vector.
         * @return The dot product.
         */
        private double depth(int ip, double[] alt) {
            double[] vect1 = new double[3];
            double[] vect2 = new double[3];
            double[] vect3 = new double[3];
            double[] vect4 = new double[3];
            int ia1 = probeAtomNumbers[0][ip];
            int ia2 = probeAtomNumbers[1][ip];
            int ia3 = probeAtomNumbers[2][ip];
            for (int k = 0; k < 3; k++) {
                vect1[k] = atomCoords[k][ia1] - atomCoords[k][ia3];
                vect2[k] = atomCoords[k][ia2] - atomCoords[k][ia3];
                vect3[k] = probeCoords[k][ip] - atomCoords[k][ia3];
            }
            cross(vect1, vect2, vect4);
            norm(vect4, vect4);
            double dot = dot(vect4, vect3);
            arraycopy(vect4, 0, alt, 0, 3);
            return dot;
        }

        /**
         * Constrain angle to be in the range -1.0 <= angle <= 1.0.
         *
         * @param angle The value to check.
         * @return the constrained value.
         */
        private double check(double angle) {
            return max(min(angle, 1.0), -1.0);
        }

        /**
         * Row major indexing (the last dimension is contiguous in memory).
         *
         * @param i x index.
         * @param j y index.
         * @param k z index.
         * @return the row major index.
         */
        private int index2(int i, int j, int k) {
            return k + MXCUBE * (j + MXCUBE * i);
        }

        private void getVector(double[] ai, double[][] temp, int index) {
            ai[0] = temp[0][index];
            ai[1] = temp[1][index];
            ai[2] = temp[2][index];
        }

        private void putVector(double[] ai, double[][] temp, int index) {
            temp[0][index] = ai[0];
            temp[1][index] = ai[1];
            temp[2][index] = ai[2];
        }

        private void getVector(double[] ai, double[][][] temp, int index1, int index2) {
            ai[0] = temp[0][index1][index2];
            ai[1] = temp[1][index1][index2];
            ai[2] = temp[2][index1][index2];
        }

        private void putVector(double[] ai, double[][][] temp, int index1, int index2) {
            temp[0][index1][index2] = ai[0];
            temp[1][index1][index2] = ai[1];
            temp[2][index1][index2] = ai[2];
        }

        private void computeVolumeGradient() {
            int[][] cube = new int[2][MXCUBE * MXCUBE * MXCUBE];
            int[] inov = new int[MAXARC];
            double[] arci = new double[MAXARC];
            double[] arcf = new double[MAXARC];
            double[] dx = new double[MAXARC];
            double[] dy = new double[MAXARC];
            double[] dsq = new double[MAXARC];
            double[] d = new double[MAXARC];

            // Initialize minimum and maximum range of atoms.
            double rmax = 0.0;
            double xmin = x[0];
            double xmax = x[0];
            double ymin = y[0];
            double ymax = y[0];
            double zmin = z[0];
            double zmax = z[0];
            for (int i = 0; i < nAtoms; i++) {
                if (radius[i] > 0.0) {
                    if (radius[i] > rmax) {
                        rmax = radius[i];
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
            if (max(max(nx, ny), nz) > MXCUBE) {
                throw new EnergyException(" VolumeRegion derivative  --  Increase the Value of mxcube", false);
            }

            // Initialize the coarse lattice of cubes.
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
                double rr = radius[ir];
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
                                            logger.severe(" VolumeRegion gradient  --  Increase the Value of MAXARC");
                                        }
                                        dx[io] = x[in] - xr;
                                        dy[io] = y[in] - yr;
                                        dsq[io] = (dx[io] * dx[io]) + (dy[io] * dy[io]);
                                        double dist2 = dsq[io] + ((z[in] - zr) * (z[in] - zr));
                                        double vdwsum = (rr + radius[in]) * (rr + radius[in]);
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
                        double twoPI = 2.0 * PI;
                        for (int k = 0; k <= io; k++) {
                            in = inov[k];
                            double rinsq = radius[in] * radius[in];
                            double rsec2n = rinsq - ((zgrid - z[in]) * (zgrid - z[in]));
                            if (rsec2n > 0.0) {
                                double rsecn = sqrt(rsec2n);
                                if (d[k] < (rsecr + rsecn)) {
                                    double rdiff = rsecr - rsecn;
                                    if (d[k] <= abs(rdiff)) {
                                        if (rdiff < 0.0) {
                                            narc = 0;
                                            arci[narc] = 0.0;
                                            arcf[narc] = twoPI;
                                        }
                                        continue;
                                    }
                                    narc++;
                                    if (narc > MAXARC) {
                                        logger.severe(" VolumeRegion gradient  --  Increase the Value of MAXARC");
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
                                        beta += twoPI;
                                    }
                                    double ti = beta - alpha;
                                    double tf = beta + alpha;
                                    if (ti < 0.0) {
                                        ti += twoPI;
                                    }
                                    if (tf > twoPI) {
                                        tf -= twoPI;
                                    }
                                    arci[narc] = ti;
                                    // If the arc crosses zero, then it is broken into two segments; the first
                                    // ends at two pi and the second begins at zero.
                                    if (tf < ti) {
                                        arcf[narc] = twoPI;
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
                            seg_dz = twoPI * ((cos_phi1 * cos_phi1) - (cos_phi2 * cos_phi2));
                            pre_dz += seg_dz;
                        } else {
                            // Sort the arc endpoint arrays, each with
                            // "narc" entries, in order of increasing values
                            // of the arguments in "arci".
                            for (int k = 0; k < narc; k++) {
                                double aa = arci[k];
                                double bb = arcf[k];
                                double temp = 1000000.0;
                                int itemp = 0;
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
                            double temp = arcf[0];
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
                                arcf[1] = twoPI;
                                arci[1] = arcf[0];
                                arcf[0] = arci[0];
                                arci[0] = 0.0;
                            } else {
                                temp = arci[0];
                                for (int k = 0; k < narc; k++) {
                                    arci[k] = arcf[k];
                                    arcf[k] = arci[k + 1];
                                }

                                if (temp == 0.0 && arcf[narc] == twoPI) {
                                    narc--;
                                } else {
                                    arci[narc] = arcf[narc];
                                    arcf[narc] = temp;
                                }
                            }

                            // Compute the numerical pre-derivative values.
                            for (int k = 0; k <= narc; k++) {
                                double theta1 = arci[k];
                                double theta2 = arcf[k];
                                double dtheta;
                                if (theta2 >= theta1) {
                                    dtheta = theta2 - theta1;
                                } else {
                                    dtheta = (theta2 + twoPI) - theta1;
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
    }
}