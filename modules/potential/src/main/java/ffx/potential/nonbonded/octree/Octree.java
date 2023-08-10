//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.potential.nonbonded.octree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import static ffx.utilities.Constants.dWater;
import static java.lang.String.format;
import static ffx.potential.parameters.MultipoleType.t000;
import static ffx.potential.parameters.MultipoleType.t001;
import static ffx.potential.parameters.MultipoleType.t002;
import static ffx.potential.parameters.MultipoleType.t010;
import static ffx.potential.parameters.MultipoleType.t011;
import static ffx.potential.parameters.MultipoleType.t020;
import static ffx.potential.parameters.MultipoleType.t100;
import static ffx.potential.parameters.MultipoleType.t101;
import static ffx.potential.parameters.MultipoleType.t110;
import static ffx.potential.parameters.MultipoleType.t200;
import static ffx.utilities.Constants.ELEC_ANG_TO_DEBYE;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.numerics.atomic.AtomicDoubleArray;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.numerics.multipole.GKEnergyQI;
import ffx.numerics.multipole.PolarizableMultipole;
import ffx.numerics.multipole.QIFrame;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.nonbonded.pme.Polarization;
import ffx.potential.parameters.ForceField;
import ffx.utilities.Constants;
import org.apache.commons.lang3.BooleanUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

/**
 * Octree method presented in the Fast Multipole Method (FMM) tutorial from the Barba Group:
 * https://github.com/barbagroup/FMM_tutorial
 */
public class Octree {

    private static final Logger logger = Logger.getLogger(Octree.class.getName());
    /**
     * List of all leaf cells
     */
    private final ArrayList<OctreeCell> leaves = new ArrayList<>();
    /**
     * Critical (maximum allowed) number of points allowed in any one cell: If a cell already contains
     * nCritical points, it needs to be split
     */
    private final int nCritical;
    /**
     * List of particles (atoms)
     */
    private final Atom[] particles;
    /**
     * Tolerance parameter
     */
    private final double theta;

    private final double[][][] globalMultipole;
    /**
     * List of cells
     */
    private ArrayList<OctreeCell> cells = new ArrayList<>();

    private MolecularAssembly molecularAssembly;

    private double gkPermanentEnergy = 0.0;
    private double gkPolarizationEnergy = 0.0;
    private double gkEnergy = 0.0;

    protected AtomicDoubleArray3D torque;
    protected AtomicDoubleArray3D grad;
    protected AtomicDoubleArray sharedBornGrad;
    protected AtomicDoubleArray selfEnergy;
    protected AtomicDoubleArray crossEnergy;

    private boolean gradient;

    private double[] selfEnergies;

    private int[] interactionCount;

    private double[] phi;

    private double[] smallestSide = new double[2];
    private int[][] interactionTracker;

    /**
     * Default constructor: only need to pass in a list of particles nCritical and theta set to
     * defaults
     *
     * @param particles ArrayList of type OctreeParticles of all particles to be used in tree
     */
//  public Octree(Atom[] particles) {
//    this.particles = particles;
//    this.nCritical = 10;
//    this.theta = 0.5;
//  }

    /**
     * Constructor allowing the specification of nCritical, default theta value
     *
     * @param nCritical Critical number of particles; cells must split when they reach nCritical
     * @param particles ArrayList of type OctreeParticles of all particles to be used in tree
     */
//  public Octree(int nCritical, Atom[] particles) {
//    this.nCritical = nCritical;
//    this.particles = particles;
//    this.theta = 0.5;
//  }

    /**
     * Constructor allowing the specification of nCritical and theta
     *
     * @param nCritical Critical number of particles; cells must split when they reach nCritical
     * @param particles ArrayList of type OctreeParticles of all particles to be used in tree
     * @param theta     Specifies near field vs far field
     */
    public Octree(int nCritical, Atom[] particles, double theta, double[][][] globalMultipole, MolecularAssembly molecularAssembly, boolean gradient, double[] selfEnergies) {
        this.nCritical = nCritical;
        this.particles = particles;
        this.theta = theta;
        this.globalMultipole = globalMultipole;
        this.phi = new double[particles.length];
        Arrays.fill(phi, 0.0);
        this.interactionCount = new int[particles.length];
        Arrays.fill(interactionCount, 0);
        this.molecularAssembly = molecularAssembly;
        this.gradient = gradient;
        this.selfEnergies = selfEnergies;
        this.interactionTracker = new int[particles.length][particles.length];
//        Arrays.fill(interactionTracker, 0);
    }

    /**
     * Add a child.
     *
     * @param octant The octant.
     * @param p      Cell index p.
     */
    public void addChild(int octant, int p) {
        OctreeCell tempCell = new OctreeCell(nCritical);

        // Create new cell
        // TODO: should the cells array be passed in or the "this" value from the particular Octree?
        cells.add(tempCell);

        // Last element of cells list is new child, c
        int c = cells.size() - 1;

        // Geometric reference between parent and child
        cells.get(c).setR((cells.get(p).getR()) * 0.5);
//    logger.info(format("oct & 1 %d",(octant & 1)));
//    logger.info(format("oct & 2 %d",(octant & 2)));
//    logger.info(format("oct & 4 %d",(octant & 4)));
        cells.get(c).setX((cells.get(p).getX()) + cells.get(c).getR() * ((octant & 1) * 2 - 1));
        cells.get(c).setY((cells.get(p).getY()) + cells.get(c).getR() * ((octant & 2) - 1));
        cells.get(c).setZ((cells.get(p).getZ()) + cells.get(c).getR() * ((octant & 4) * 0.5 - 1));

        // Establish mutual reference in cells list
        cells.get(c).setParentIndex(p);
        cells.get(p).setChildren(octant, c);
        // Bitwise operation to determine empty/non-empty child cells
//    logger.info(format("nChild before bitwise operation %d",cells.get(p).getnChild()));
        cells.get(p).setnChild(cells.get(p).getnChild() | (1 << octant));
        logger.info(format("+++cell %d is created as a child of cell %d", c, p));
//    logger.info(format("** Cell %d sidelength = %4.3f",c,cells.get(c).getR()*2));
        if (smallestSide[1] > cells.get(c).getR() * 2) {
            smallestSide[0] = c;
            smallestSide[1] = cells.get(c).getR() * 2;
        }
//    logger.info(format("nChild after  bitwise operation %d",cells.get(p).getnChild()));
    }

    /**
     * Build the tree.
     *
     * @param root The Octree root cell.
     */
    public void buildTree(OctreeCell root) {
        // (Re)Initialize cells to an empty array list
        cells = new ArrayList<>();

        // Set root cell - add it to the cells array list
        cells.add(root);
        smallestSide[0] = 0.0;
        smallestSide[1] = cells.get(0).getR() * 2;

        // Build tree
        int n = particles.length;

        for (int i = 0; i < n; i++) {
            int current = 0;

//      logger.info(format("nLeaf at i %d : %d",i,cells.get(current).getNumLeaves()));

            while (cells.get(current).getNumLeaves() >= nCritical) {
                cells.get(current).setNumLeaves(cells.get(current).getNumLeaves() + 1);
                int octX = 0;
                int octY = 0;
                int octZ = 0;

                if (particles[i].getX() > cells.get(current).getX()) {
                    octX = 1;
                }
                if (particles[i].getY() > cells.get(current).getY()) {
                    octY = 1;
                }
                if (particles[i].getZ() > cells.get(current).getZ()) {
                    octZ = 1;
                }

                // Find particle's octant - should be an integer from 0 to 7
                int octant = octX + (octY << 1) + (octZ << 2);

                // If there's not a child cell in the particle's octant, create one
                boolean noChildInOctant =
                        BooleanUtils.toBoolean(cells.get(current).getnChild() & (1 << octant));
                if (!noChildInOctant) {
                    addChild(octant, current);
                }

                current = cells.get(current).getChildAtIndex(octant);
            }

            // Allocate the particle in the leaf cell
            cells.get(current).setLeaf(cells.get(current).getNumLeaves(), i);
            cells.get(current).setNumLeaves(cells.get(current).getNumLeaves() + 1);
            logger.info(format("particle %d stored in cell %d", i, current));

//      logger.info(format("nLeaf at i %d after allocating particle: %d",i,cells.get(current).getNumLeaves()));

            // Check whether to split cell
            if (cells.get(current).getNumLeaves() >= nCritical) {
//        logger.info("Splitting cells in buildTree");
                splitCell(current);
            }
        }

        // return cells;
    }

    /**
     * Compute multipole moments for an array of atoms using the geometric center of a domain decomposition box.
     * Atoms passed in should be from the associated domain decomposition Cell
     *
     * @param cellAtoms Atom array to consider.
     * @param cell      NeighborList Cell containing cell index information and atom index information
     */
    public double[] computeMomentsGeometric(Atom[] cellAtoms, OctreeCell cell) {
//    logger.info("** Enters computeMomentsGeometric");
        // Zero out total charge, dipole and quadrupole components.
        var netchg = 0.0;
        var netdpl = 0.0;
        var xdpl = 0.0;
        var ydpl = 0.0;
        var zdpl = 0.0;
        var xxqdp = 0.0;
        var xyqdp = 0.0;
        var xzqdp = 0.0;
        var yxqdp = 0.0;
        var yyqdp = 0.0;
        var yzqdp = 0.0;
        var zxqdp = 0.0;
        var zyqdp = 0.0;
        var zzqdp = 0.0;

        // Find the geometric center of the cell.
        double xmid = cell.getX();
        double ymid = cell.getY();
        double zmid = cell.getZ();

//    logger.info(format("Geometric Center of Cell : %4.4f , %4.4f , %4.4f",xmid,ymid,zmid));

        // Get atoms within OctreeCell
        int n = cellAtoms.length;

        double[] xcm = new double[n];
        double[] ycm = new double[n];
        double[] zcm = new double[n];
        int k = 0;
        for (Atom atom : cellAtoms) {
            xcm[k] = atom.getX() - xmid;
            ycm[k] = atom.getY() - ymid;
            zcm[k] = atom.getZ() - zmid;
            k++;
        }

        // Account for charge, dipoles and induced dipoles.
        k = 0;
        for (Atom atom : cellAtoms) {
            int i = atom.getIndex() - 1;
            double[] globalMultipolei = globalMultipole[0][i];
//      double[] inducedDipolei = inducedDipole[0][i];

            var ci = globalMultipolei[t000];
            var dix = globalMultipolei[t100];
            var diy = globalMultipolei[t010];
            var diz = globalMultipolei[t001];
//      var uix = inducedDipolei[0];
//      var uiy = inducedDipolei[1];
//      var uiz = inducedDipolei[2];
//      logger.info(format("Global Multipole for atom %d: %4.3f %4.3f %4.3f %4.3f",i,ci,dix,diy,diz));

            netchg += ci;
            xdpl += xcm[k] * ci + dix;
            ydpl += ycm[k] * ci + diy;
            zdpl += zcm[k] * ci + diz;
            xxqdp += xcm[k] * xcm[k] * ci + 2.0 * xcm[k] * (dix);
            xyqdp += xcm[k] * ycm[k] * ci + xcm[k] * (diy) + ycm[k] * (dix);
            xzqdp += xcm[k] * zcm[k] * ci + xcm[k] * (diz) + zcm[k] * (dix);
            yxqdp += ycm[k] * xcm[k] * ci + ycm[k] * (dix) + xcm[k] * (diy);
            yyqdp += ycm[k] * ycm[k] * ci + 2.0 * ycm[k] * (diy);
            yzqdp += ycm[k] * zcm[k] * ci + ycm[k] * (diz) + zcm[k] * (diy);
            zxqdp += zcm[k] * xcm[k] * ci + zcm[k] * (dix) + xcm[k] * (diz);
            zyqdp += zcm[k] * ycm[k] * ci + zcm[k] * (diy) + ycm[k] * (diz);
            zzqdp += zcm[k] * zcm[k] * ci + 2.0 * zcm[k] * (diz);

//      logger.info(format("Cell %d   Charge %4.4f    Dipole %4.4f %4.4f %4.4f",cell.getParentIndex(),netchg,xdpl,ydpl,zdpl));
//      logger.info(format("Cell ##   Charge %4.4f    Dipole %4.4f %4.4f %4.4f",netchg,xdpl,ydpl,zdpl));

//      xdpl += xcm[k] * ci + dix + uix;
//      ydpl += ycm[k] * ci + diy + uiy;
//      zdpl += zcm[k] * ci + diz + uiz;
//      xxqdp += xcm[k] * xcm[k] * ci + 2.0 * xcm[k] * (dix + uix);
//      xyqdp += xcm[k] * ycm[k] * ci + xcm[k] * (diy + uiy) + ycm[k] * (dix + uix);
//      xzqdp += xcm[k] * zcm[k] * ci + xcm[k] * (diz + uiz) + zcm[k] * (dix + uix);
//      yxqdp += ycm[k] * xcm[k] * ci + ycm[k] * (dix + uix) + xcm[k] * (diy + uiy);
//      yyqdp += ycm[k] * ycm[k] * ci + 2.0 * ycm[k] * (diy + uiy);
//      yzqdp += ycm[k] * zcm[k] * ci + ycm[k] * (diz + uiz) + zcm[k] * (diy + uiy);
//      zxqdp += zcm[k] * xcm[k] * ci + zcm[k] * (dix + uix) + xcm[k] * (diz + uiz);
//      zyqdp += zcm[k] * ycm[k] * ci + zcm[k] * (diy + uiy) + ycm[k] * (diz + uiz);
//      zzqdp += zcm[k] * zcm[k] * ci + 2.0 * zcm[k] * (diz + uiz);
//      logger.info(format("Multipole for atom %d: %f4.3 %f4.3 %f4.3 %f4.3 %f4.3 %f4.3 %f4.3 %f4.3 %f4.3 %f4.3",
//              i,netchg,xdpl,ydpl,zdpl,xxqdp,xyqdp,xzqdp,yxqdp,yyqdp,yzqdp,zxqdp,zyqdp,zzqdp));
            k++;
        }

        // Convert  the quadrupole from traced to traceless form.
        //TODO: add this
        var qave = (xxqdp + yyqdp + zzqdp) / 3.0;
        var xx = 1.5 * (xxqdp - qave);
        var xy = 1.5 * xyqdp;
        var xz = 1.5 * xzqdp;
        var yx = 1.5 * yxqdp;
        var yy = 1.5 * (yyqdp - qave);
        var yz = 1.5 * yzqdp;
        var zx = 1.5 * zxqdp;
        var zy = 1.5 * zyqdp;
        var zz = 1.5 * (zzqdp - qave);
//
        // Add the traceless atomic quadrupoles to total quadrupole.
        //TODO: add this back in for testing
        for (Atom atom : cellAtoms) {
          int i = atom.getIndex() - 1;
          double[] globalMultipolei = globalMultipole[0][i];
          var qixx = globalMultipolei[t200];
          var qiyy = globalMultipolei[t020];
          var qizz = globalMultipolei[t002];
          var qixy = globalMultipolei[t110];
          var qixz = globalMultipolei[t101];
          var qiyz = globalMultipolei[t011];
          xx += qixx;
          xy += qixy;
          xz += qixz;
          yx += qixy;
          yy += qiyy;
          yz += qiyz;
          zx += qixz;
          zy += qiyz;
          zz += qizz;
        }
        double[] calculatedTracelessQDP = {xx, xy, xz, yx, yy, yz, zx, zy, zz};
        cell.addToTracelessQuadrapole(calculatedTracelessQDP);
//
//    // Convert dipole to Debye and quadrupole to Buckingham.
//    xdpl = xdpl * ELEC_ANG_TO_DEBYE;
//    ydpl = ydpl * ELEC_ANG_TO_DEBYE;
//    zdpl = zdpl * ELEC_ANG_TO_DEBYE;
//    xxqdp = xxqdp * ELEC_ANG_TO_DEBYE;
//    xyqdp = xyqdp * ELEC_ANG_TO_DEBYE;
//    xzqdp = xzqdp * ELEC_ANG_TO_DEBYE;
//    yxqdp = yxqdp * ELEC_ANG_TO_DEBYE;
//    yyqdp = yyqdp * ELEC_ANG_TO_DEBYE;
//    yzqdp = yzqdp * ELEC_ANG_TO_DEBYE;
//    zxqdp = zxqdp * ELEC_ANG_TO_DEBYE;
//    zyqdp = zyqdp * ELEC_ANG_TO_DEBYE;
//    zzqdp = zzqdp * ELEC_ANG_TO_DEBYE;
//
//    // Get dipole magnitude and diagonalize quadrupole tensor.
//    netdpl = sqrt(xdpl * xdpl + ydpl * ydpl + zdpl * zdpl);
//    double[][] a = new double[3][3];
//    a[0][0] = xxqdp;
//    a[0][1] = xyqdp;
//    a[0][2] = xzqdp;
//    a[1][0] = yxqdp;
//    a[1][1] = yyqdp;
//    a[1][2] = yzqdp;
//    a[2][0] = zxqdp;
//    a[2][1] = zyqdp;
//    a[2][2] = zzqdp;
//    EigenDecomposition e = new EigenDecomposition(new Array2DRowRealMatrix(a));
//    // Eigenvalues are returned in descending order, but logged below in ascending order.
//    var netqdp = e.getRealEigenvalues();

//    logger.info(format("\n Electric Moments for Cell %d\n",cell.getParentIndex()));
//    logger.info(format("\n Electric Moments for Cell ##\n"));
//    logger.info(format("  Total Electric Charge:    %13.5f Electrons\n", netchg));
//    logger.info(format("  Dipole Moment Magnitude:  %13.5f Debye\n", netdpl));
//    logger.info(format("  Dipole X,Y,Z-Components:  %13.5f %13.5f %13.5f\n", xdpl, ydpl, zdpl));
//    logger.info(format("  Quadrupole Moment Tensor: %13.5f %13.5f %13.5f", xxqdp, xyqdp, xzqdp));
//    logger.info(format("       (Buckinghams)        %13.5f %13.5f %13.5f", yxqdp, yyqdp, yzqdp));
//    logger.info(format("                            %13.5f %13.5f %13.5f\n", zxqdp, zyqdp, zzqdp));
//    logger.info(
//            format(
//                    "  Principal Axes Quadrupole %13.5f %13.5f %13.5f\n", netqdp[2], netqdp[1], netqdp[0]));

        return new double[]{netchg, xdpl, ydpl, zdpl, xxqdp, yyqdp, zzqdp, xyqdp, xzqdp, yzqdp};
    }

//  /** Direct summation. */
//  public void directSum() {
//    for (int i = 0; i < particles.length; i++) {
//      for (int j = 0; j < particles.length; j++) {
//        if (j != i) {
////          double r = particles[i].distance(particles[j]);
//          double r = distance(particles[i].getXYZ(null),particles[j].getXYZ(null));
////          particles[j].addToPhi(particles[j].getCharge() / r);
//        }
//      }
//    }
    // Reset potential for all particles
//    for (Atom particle : particles) {
//      particle.addToPhi(0);
//    }
//  }

    /**
     * Compute a distance between a position and the coordinates of another atom or Octree Cell.
     *
     * @param array Position.
     * @param other Position of another cell/atom.
     * @return Returns the distance.
     */
    public double distance(double[] array, double[] other) {
        return Math.sqrt(
                Math.pow((array[0] - other[0]), 2)
                        + Math.pow((array[1] - other[1]), 2)
                        + Math.pow((array[2] - other[2]), 2));
    }
    //  public double distance(double[] array, OctreeCell other) {
//    return Math.sqrt(
//        Math.pow((array[0] - other.getX()), 2)
//            + Math.pow((array[1] - other.getY()), 2)
//            + Math.pow((array[2] - other.getZ()), 2));
//  }

    /**
     * Evaluate potential at all target points
     */
    public void evalPotential() {
        for (int i = 0; i < particles.length; i++) {
            evalAtTarget(0, i);
//            logger.info(format("Number of interactions for atom %d = %d", i, interactionCount[i]));
            logger.info(format("GK energy = %4.4f", gkEnergy));
//      logger.info(format("Atom %d potential = %4.3f",i,phi[i]));
        }
    }

    /**
     * Compute interactions between two atoms.
     *
     * @param i Atom index of atom of interest.
     * @param k Atom index of reaction atom.
     */
    private void interaction(int i, int k) {
        double[] dx_local = new double[3];
        double trqxi = 0.0;
        double trqyi = 0.0;
        double trqzi = 0.0;
        double dborni = 0.0;
        double dedxi = 0.0;
        double dedyi = 0.0;
        double dedzi = 0.0;
        final double[] gradI = new double[3];
        final double[] torqueI = new double[3];
        final double[] torqueK = new double[3];
        final double[] gI = new double[3];
        final double[] tI = new double[3];
        final double[] tK = new double[3];
        double[][] transOp = new double[3][3];
        transOp[0][0] = 1.0;
        transOp[1][1] = 1.0;
        transOp[2][2] = 1.0;
//        logger.info(format("transOP: %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f",
//                transOp[0][0],transOp[0][1],transOp[0][2],transOp[1][0],transOp[1][1],transOp[1][2],transOp[2][0],transOp[2][1],transOp[2][2]));

        dx_local[0] = particles[k].getX() - particles[i].getX();
        dx_local[1] = particles[k].getY() - particles[i].getY();
        dx_local[2] = particles[k].getZ() - particles[i].getZ();
        double r2 = dx_local[0] * dx_local[0] + dx_local[1] * dx_local[1] + dx_local[2] * dx_local[2];

        // Set the multipole moments for site I.
        PolarizableMultipole mI = new PolarizableMultipole();
        PolarizableMultipole mK = new PolarizableMultipole();
        mI.setPermanentMultipole(globalMultipole[0][i]);
        // Set the multipole moments for site K.
        mK.setPermanentMultipole(globalMultipole[0][k]);

        double selfScale = 1.0;
        QIFrame qiFrame = new QIFrame();
        if (i == k) {
            selfScale = 0.5;
        } else {
            qiFrame.setAndRotate(dx_local, mI, mK);
        }

        // Set the GK energy tensor.
        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        GeneralizedKirkwood generalizedKirkwood = forceFieldEnergy.getGK();
        torque = generalizedKirkwood.getTorque();
        grad = generalizedKirkwood.getGrad();
        sharedBornGrad = AtomicDoubleArray.atomicDoubleArrayFactory(
                AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI, 1, particles.length);
        selfEnergy = AtomicDoubleArray.atomicDoubleArrayFactory(
                AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI, 1, particles.length);
        crossEnergy = AtomicDoubleArray.atomicDoubleArrayFactory(
                AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI, 1, particles.length);

        double[] born = generalizedKirkwood.getBorn();
        double rbi = born[i];
        double rbk = born[k];
        ForceField forceField = molecularAssembly.getForceField();
        double soluteDielectric = forceField.getDouble("SOLUTE_DIELECTRIC", 1.0);
        double solventDielectric = forceField.getDouble("SOLVENT_DIELECTRIC", dWater);
        double gkc = forceField.getDouble("GKC", 2.455);
        GKEnergyQI gkEnergyQI = new GKEnergyQI(soluteDielectric, solventDielectric, gkc, gradient);
        gkEnergyQI.initPotential(dx_local, r2, rbi, rbk);

        // Compute the GK permanent multipole interaction.
        double eik;

        String polar = forceField.getString("POLARIZATION", "MUTUAL");
        Polarization polarization;
        boolean polarizationTerm = forceField.getBoolean("POLARIZETERM", true);
        if (!polarizationTerm || polar.equalsIgnoreCase("NONE")) {
            polarization = Polarization.NONE;
        } else if (polar.equalsIgnoreCase("DIRECT")) {
            polarization = Polarization.DIRECT;
        } else {
            polarization = Polarization.MUTUAL;
        }
        double electric = forceField.getDouble("ELECTRIC", Constants.DEFAULT_ELECTRIC);

        if (!gradient) {
            // Compute the GK permanent interaction energy.
            eik = electric * selfScale * gkEnergyQI.multipoleEnergy(mI, mK);
            gkPermanentEnergy += eik;

            if (polarization != Polarization.NONE) {
                // Compute the GK polarization interaction energy.
                double ep = electric * selfScale * gkEnergyQI.polarizationEnergy(mI, mK);
                gkPolarizationEnergy += ep;
                eik += ep;
            }
        } else {
            if (i == k) {
                // Compute the GK permanent interaction energy.
                eik = electric * selfScale * gkEnergyQI.multipoleEnergy(mI, mK);
                gkPermanentEnergy += eik;
                // There is no dE/dR or torque for the i == k permanent self-energy.
                if (polarization != Polarization.NONE) {
                    // Compute the GK polarization interaction energy.
                    // There is no dE/dR, but there can be a torque for the i == k polarization self-energy
                    // if the permanent dipole and induced dipole are not aligned.
                    double mutualMask = 1.0;
                    if (polarization == Polarization.DIRECT) {
                        mutualMask = 0.0;
                    }
                    double ep =
                            electric * selfScale * gkEnergyQI.polarizationEnergyAndGradient(mI, mK, mutualMask,
                                    gradI, torqueI, torqueK);
                    gkPolarizationEnergy += ep;
                    eik += ep;
                    // The torque for the i == k polarization self-energy.
                    trqxi += electric * selfScale * torqueI[0];
                    trqyi += electric * selfScale * torqueI[1];
                    trqzi += electric * selfScale * torqueI[2];
                    double tkx = electric * selfScale * torqueK[0];
                    double tky = electric * selfScale * torqueK[1];
                    double tkz = electric * selfScale * torqueK[2];
                    final double rtkx = tkx * transOp[0][0] + tky * transOp[1][0] + tkz * transOp[2][0];
                    final double rtky = tkx * transOp[0][1] + tky * transOp[1][1] + tkz * transOp[2][1];
                    final double rtkz = tkx * transOp[0][2] + tky * transOp[1][2] + tkz * transOp[2][2];
                    torque.add(0, k, rtkx, rtky, rtkz);
                }
                // Compute the GK permanent Born chain-rule term.
                gkEnergyQI.initBorn(dx_local, r2, rbi, rbk);
                double db = gkEnergyQI.multipoleEnergyBornGrad(mI, mK);
                if (polarization != Polarization.NONE) {
                    db += gkEnergyQI.polarizationEnergyBornGrad(mI, mK, polarization == Polarization.MUTUAL);
                }
                dborni += electric * rbi * db;
            } else {
                // Sum the GK permanent interaction energy.
                eik = electric * gkEnergyQI.multipoleEnergyAndGradient(mI, mK, gradI, torqueI, torqueK);
                gkPermanentEnergy += eik;
                if (polarization != Polarization.NONE) {
                    double mutualMask = 1.0;
                    if (polarization == Polarization.DIRECT) {
                        mutualMask = 0.0;
                    }
                    // Sum the GK polarization interaction energy.
                    double ep =
                            electric * gkEnergyQI.polarizationEnergyAndGradient(mI, mK, mutualMask, gI, tI, tK);
                    eik += ep;
                    gkPolarizationEnergy += ep;
                    for (int j = 0; j < 3; j++) {
                        gradI[j] += gI[j];
                        torqueI[j] += tI[j];
                        torqueK[j] += tK[j];
                    }
                }

                // Convert to units of kcal/mol
                for (int j = 0; j < 3; j++) {
                    gradI[j] *= electric;
                    torqueI[j] *= electric;
                    torqueK[j] *= electric;
                }

                // Rotate gradient and torques from the QI frame into the Global frame.
                qiFrame.toGlobal(gradI);
                qiFrame.toGlobal(torqueI);
                qiFrame.toGlobal(torqueK);

                // Accumulate partial derivatives on site I.
                dedxi += gradI[0];
                dedyi += gradI[1];
                dedzi += gradI[2];
                trqxi += torqueI[0];
                trqyi += torqueI[1];
                trqzi += torqueI[2];

                // Accumulate partial derivatives on site K.
                double dedx = -gradI[0];
                double dedy = -gradI[1];
                double dedz = -gradI[2];
                // Handle symmetry mate rotation.
                final double dedxk = dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
                final double dedyk = dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
                final double dedzk = dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];
                double tkx = torqueK[0];
                double tky = torqueK[1];
                double tkz = torqueK[2];
                final double rtkx = tkx * transOp[0][0] + tky * transOp[1][0] + tkz * transOp[2][0];
                final double rtky = tkx * transOp[0][1] + tky * transOp[1][1] + tkz * transOp[2][1];
                final double rtkz = tkx * transOp[0][2] + tky * transOp[1][2] + tkz * transOp[2][2];
                grad.add(0, k, dedxk, dedyk, dedzk);
                torque.add(0, k, rtkx, rtky, rtkz);

                // Compute the Born-chain rule term.
                gkEnergyQI.initBorn(dx_local, r2, rbi, rbk);
                double db = gkEnergyQI.multipoleEnergyBornGrad(mI, mK);
                if (polarization != Polarization.NONE) {
                    db += gkEnergyQI.polarizationEnergyBornGrad(mI, mK, polarization == Polarization.MUTUAL);
                }
                db *= electric;
                dborni += rbk * db;
                sharedBornGrad.add(0, k, rbi * db);
            }
        }

        if (i == k) {
            selfEnergy.add(0, i, eik);
        } else {
            double half = 0.5 * eik;
            crossEnergy.add(0, i, half);
            crossEnergy.add(0, k, half);
        }

//        logger.info(format("Atom %d and %d interaction eik = %4.4f",i,k,eik));
        gkEnergy += eik;
//    count++;
    }

    /**
     * Compute interactions between an atom and an Octree cell.
     *
     * @param i Atom index.
     * @param k Cell index.
     */
    private void interactionWithCell(int i, int k) {
        double[] dx_local = new double[3];
        double trqxi = 0.0;
        double trqyi = 0.0;
        double trqzi = 0.0;
        double dborni = 0.0;
        double dedxi = 0.0;
        double dedyi = 0.0;
        double dedzi = 0.0;
        final double[] gradI = new double[3];
        final double[] torqueI = new double[3];
        final double[] torqueK = new double[3];
        final double[] gI = new double[3];
        final double[] tI = new double[3];
        final double[] tK = new double[3];
        double[][] transOp = new double[3][3];
        transOp[0][0] = 1.0;
        transOp[1][1] = 1.0;
        transOp[2][2] = 1.0;

        dx_local[0] = cells.get(k).getX() - particles[i].getX();
        dx_local[1] = cells.get(k).getY() - particles[i].getY();
        dx_local[2] = cells.get(k).getZ() - particles[i].getZ();
        double r2 = dx_local[0] * dx_local[0] + dx_local[1] * dx_local[1] + dx_local[2] * dx_local[2];

        // Set the multipole moments for site I.
        PolarizableMultipole mI = new PolarizableMultipole();
        PolarizableMultipole mK = new PolarizableMultipole();
        mI.setPermanentMultipole(globalMultipole[0][i]);
//        mK.setPermanentMultipole(cells.get(k).getMultipole());

        double[] multipoleCell = cells.get(k).getMultipole();
        double[] tracelessQDBCell = cells.get(k).getTracelessQDP();
        multipoleCell[4] = tracelessQDBCell[0] * 3.0; // xx
        multipoleCell[5] = tracelessQDBCell[4] * 3.0; // yy
        multipoleCell[6] = tracelessQDBCell[8] * 3.0; // zz
        multipoleCell[7] = tracelessQDBCell[1] * 3.0/2.0; // xy
        multipoleCell[8] = tracelessQDBCell[2] * 3.0/2.0; // xz
        multipoleCell[9] = tracelessQDBCell[5] * 3.0/2.0; // yz
        //TODO: Multiplying by 3 or 3/2 is done to try to cancel the effect in PolarizableMultipole.setMultipole()
        mK.setPermanentMultipole(multipoleCell);

        // Set the multipole moments for site K.
//        double[] multipoleCell = cells.get(k).getMultipole();
//        double r = Math.sqrt(r2);
//        double r3 = Math.pow(r,3);
//        double r5 = Math.pow(r,5);
//        double dx = dx_local[0];
//        double dy = dx_local[1];
//        double dz = dx_local[2];
//        double[] weight = new double[10];
//        weight[0] = 1 / r;
//        weight[1] = -dx / r3;
//        weight[2] = -dy / r3;
//        weight[3] = -dz / r3;
//        weight[4] = (3 * Math.pow(dx, 2)) / r5 - (1 / r3);
//        weight[5] = (3 * Math.pow(dy, 2)) / r5 - (1 / r3);
//        weight[6] = (3 * Math.pow(dz, 2)) / r5 - (1 / r3);
//        weight[7] = 3 * dx * dy / r5;
//        weight[8] = 3 * dx * dz / r5;
//        weight[9] = 3 * dy * dz / r5;
//        for (int q = 0; q < 10; q++) {
//            multipoleCell[q] = multipoleCell[q] * weight[q];
//        }
//        mK.setPermanentMultipole(multipoleCell);

        double selfScale = 1.0;
        QIFrame qiFrame = new QIFrame();

        qiFrame.setAndRotate(dx_local, mI, mK);


        // Set the GK energy tensor.
        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        GeneralizedKirkwood generalizedKirkwood = forceFieldEnergy.getGK();
        torque = generalizedKirkwood.getTorque();
        grad = generalizedKirkwood.getGrad();
        sharedBornGrad = AtomicDoubleArray.atomicDoubleArrayFactory(
                AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI, 1, particles.length);
        selfEnergy = AtomicDoubleArray.atomicDoubleArrayFactory(
                AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI, 1, particles.length);
        crossEnergy = AtomicDoubleArray.atomicDoubleArrayFactory(
                AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI, 1, particles.length);

        double[] born = generalizedKirkwood.getBorn();
        double rbi = born[i];
//        logger.info(format("Atom %d atom born radius = %4.8f",i,born[i]));

//        AtomicDoubleArray selfEnergies = generalizedKirkwood.getSelfEnergy();
        cells.get(k).setTotalBornToZero();
        calculateTotalSelfEnergy(k,k);
        cells.get(k).setCellBornRadius();
//        determineMaxBornRadius(k,k);
//        double rbk = cells.get(k).getMaxBornRadius(); // radius = max born radius of the atoms in the cell
        double rbk = cells.get(k).getCellBornRadius();
//        rbk = rbk + cells.get(k).getR(); // radius = born radius estimate + cell radius (sidelength/2)
//        double rbk = cells.get(k).getR(); // radius = cell radius (sidelength/2)
//        double rbk = Math.sqrt(r2); // radius = distance between
        logger.info(format("Cell %d cell born radius = %4.8f",k,rbk));
        logger.info(format("Cell %d radius = %4.8f",k,cells.get(k).getR()));
        logger.info(format("Cell %d distance from atom %d = %4.8f",k,i,Math.sqrt(r2)));

        ForceField forceField = molecularAssembly.getForceField();
        double soluteDielectric = forceField.getDouble("SOLUTE_DIELECTRIC", 1.0);
        double solventDielectric = forceField.getDouble("SOLVENT_DIELECTRIC", dWater);
        double gkc = forceField.getDouble("GKC", 2.455);
        GKEnergyQI gkEnergyQI = new GKEnergyQI(soluteDielectric, solventDielectric, gkc, gradient);
        gkEnergyQI.initPotential(dx_local, r2, rbi, rbk);

        // Compute the GK permanent multipole interaction.
        double eik;

        String polar = forceField.getString("POLARIZATION", "MUTUAL");
        Polarization polarization;
        boolean polarizationTerm = forceField.getBoolean("POLARIZETERM", true);
        if (!polarizationTerm || polar.equalsIgnoreCase("NONE")) {
            polarization = Polarization.NONE;
        } else if (polar.equalsIgnoreCase("DIRECT")) {
            polarization = Polarization.DIRECT;
        } else {
            polarization = Polarization.MUTUAL;
        }
        double electric = forceField.getDouble("ELECTRIC", Constants.DEFAULT_ELECTRIC);

        if (!gradient) {
            // Compute the GK permanent interaction energy.
            eik = electric * selfScale * gkEnergyQI.multipoleEnergy(mI, mK);
            gkPermanentEnergy += eik;

            if (polarization != Polarization.NONE) {
                // Compute the GK polarization interaction energy.
                double ep = electric * selfScale * gkEnergyQI.polarizationEnergy(mI, mK);
                gkPolarizationEnergy += ep;
                eik += ep;
            }
        } else {
            // Sum the GK permanent interaction energy.
            eik = electric * gkEnergyQI.multipoleEnergyAndGradient(mI, mK, gradI, torqueI, torqueK);
            gkPermanentEnergy += eik;
            if (polarization != Polarization.NONE) {
                double mutualMask = 1.0;
                if (polarization == Polarization.DIRECT) {
                    mutualMask = 0.0;
                }
                // Sum the GK polarization interaction energy.
                double ep =
                        electric * gkEnergyQI.polarizationEnergyAndGradient(mI, mK, mutualMask, gI, tI, tK);
                eik += ep;
                gkPolarizationEnergy += ep;
                for (int j = 0; j < 3; j++) {
                    gradI[j] += gI[j];
                    torqueI[j] += tI[j];
                    torqueK[j] += tK[j];
                }
            }

            // Convert to units of kcal/mol
            for (int j = 0; j < 3; j++) {
                gradI[j] *= electric;
                torqueI[j] *= electric;
                torqueK[j] *= electric;
            }

            // Rotate gradient and torques from the QI frame into the Global frame.
            qiFrame.toGlobal(gradI);
            qiFrame.toGlobal(torqueI);
            qiFrame.toGlobal(torqueK);

            // Accumulate partial derivatives on site I.
            dedxi += gradI[0];
            dedyi += gradI[1];
            dedzi += gradI[2];
            trqxi += torqueI[0];
            trqyi += torqueI[1];
            trqzi += torqueI[2];

            // Accumulate partial derivatives on site K.
            double dedx = -gradI[0];
            double dedy = -gradI[1];
            double dedz = -gradI[2];
            // Handle symmetry mate rotation.
            final double dedxk = dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
            final double dedyk = dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
            final double dedzk = dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];
            double tkx = torqueK[0];
            double tky = torqueK[1];
            double tkz = torqueK[2];
            final double rtkx = tkx * transOp[0][0] + tky * transOp[1][0] + tkz * transOp[2][0];
            final double rtky = tkx * transOp[0][1] + tky * transOp[1][1] + tkz * transOp[2][1];
            final double rtkz = tkx * transOp[0][2] + tky * transOp[1][2] + tkz * transOp[2][2];
//            grad.add(0, k, dedxk, dedyk, dedzk);
//            torque.add(0, k, rtkx, rtky, rtkz);

            // Compute the Born-chain rule term.
            gkEnergyQI.initBorn(dx_local, r2, rbi, rbk);
            double db = gkEnergyQI.multipoleEnergyBornGrad(mI, mK);
            if (polarization != Polarization.NONE) {
                db += gkEnergyQI.polarizationEnergyBornGrad(mI, mK, polarization == Polarization.MUTUAL);
            }
            db *= electric;
            dborni += rbk * db;
//            sharedBornGrad.add(0, k, rbi * db);

        }

        double half = 0.5 * eik;
        crossEnergy.add(0, i, half);
//        crossEnergy.add(0, k, half);

//        logger.info(format("Atom %d and cell %d interaction eik = %4.4f",i,k,eik));
        gkEnergy += eik;// * cells.get(k).getNumLeaves();//2;
//    count++;
    }

    /**
     * Compute the L2 error.
     *
     * @param phiDirect Potential from direct summation.
     * @param phiTree   Potential from the tree.
     */
    public void l2Error(double[] phiDirect, double[] phiTree) {
        double errorSumNum = 0.0;
        double errorSumDenom = 0.0;
        for (int i = 0; i < phiDirect.length; i++) {
            errorSumNum = errorSumNum + Math.pow((phiDirect[i] - phiTree[i]), 2);
            errorSumDenom = errorSumDenom + Math.pow(phiDirect[i], 2);
        }
        double error = Math.sqrt(errorSumNum / errorSumDenom);
        logger.info("L2 Norm Error: " + error);
    }

    /**
     * Update sweep.
     */
    public void upwardSweep() {
        logger.info(format("******************Start of the Upward Sweep**********************"));
        for (int c = (cells.size() - 1); c > 0; c--) {
            int p = cells.get(c).getParentIndex();
            M2M(p, c);
//      double[] calculatedMultipole = cells.get(p).getMultipole();
//      logger.info(format("Geometric center of cell p%d: %4.3f %4.3f %4.3f",p,cells.get(p).getX(),cells.get(p).getY(),cells.get(p).getZ()));
//      logger.info(format("Multipole for parent cell %d : ch %4.4f ",p,calculatedMultipole[0]));
//      logger.info(format("                               dx %4.4f dy %4.4f dz %4.4f",calculatedMultipole[1],calculatedMultipole[2],calculatedMultipole[3]));
//      logger.info(format("                               dxx %4.4f dyy %4.4f dzz %4.4f",calculatedMultipole[4],calculatedMultipole[5],calculatedMultipole[6]));
//      logger.info(format("                               dxy %4.4f dxz %4.4f dyz %4.4f",calculatedMultipole[7],calculatedMultipole[8],calculatedMultipole[9]));
            if (c == 1) {
                double[] calculatedMultipole = cells.get(p).getMultipole();
                double[] calculatedTracelessQDP = cells.get(p).getTracelessQDP();
                logger.info(format("Geometric center of cell p%d: %4.3f %4.3f %4.3f", p, cells.get(p).getX(), cells.get(p).getY(), cells.get(p).getZ()));
                logger.info(format("Multipole for parent cell %d : ch %4.4f ", p, calculatedMultipole[0]));
                logger.info(format("                               dx %4.4f dy %4.4f dz %4.4f", calculatedMultipole[1], calculatedMultipole[2], calculatedMultipole[3]));
                logger.info(format("                               dxx %4.4f dyy %4.4f dzz %4.4f", calculatedMultipole[4], calculatedMultipole[5], calculatedMultipole[6]));
                logger.info(format("                               dxy %4.4f dxz %4.4f dyz %4.4f", calculatedMultipole[7], calculatedMultipole[8], calculatedMultipole[9]));
                logger.info(format("Traceless quadrapole for parent cell %d : ", p));
                logger.info(format("                                          xx %4.4f xy %4.4f xz %4.4f", calculatedTracelessQDP[0], calculatedTracelessQDP[1], calculatedTracelessQDP[2]));
                logger.info(format("                                          yx %4.4f yy %4.4f yz %4.4f", calculatedTracelessQDP[3], calculatedTracelessQDP[4], calculatedTracelessQDP[5]));
                logger.info(format("                                          zx %4.4f zy %4.4f zz %4.4f", calculatedTracelessQDP[6], calculatedTracelessQDP[7], calculatedTracelessQDP[8]));
            }
        }
        if (cells.size() == 1) {
            double[] calculatedMultipole = cells.get(0).getMultipole();
//      logger.info(format("Geometric center of cell p%d c%d : %4.3f %4.3f %4.3f",p,c,cells.get(c).getX(),cells.get(c).getY(),cells.get(c).getZ()));
            logger.info(format("Multipole for parent cell 0 : ch %4.4f ", calculatedMultipole[0]));
            logger.info(format("                               dx %4.4f dy %4.4f dz %4.4f", calculatedMultipole[1], calculatedMultipole[2], calculatedMultipole[3]));
            logger.info(format("                               dxx %4.4f dyy %4.4f dzz %4.4f", calculatedMultipole[4], calculatedMultipole[5], calculatedMultipole[6]));
            logger.info(format("                               dxy %4.4f dxz %4.4f dyz %4.4f", calculatedMultipole[7], calculatedMultipole[8], calculatedMultipole[9]));
        }
    }

    /**
     * Spit a cell.
     *
     * @param p Cell index.
     */
    private void splitCell(int p) {
        logger.info(format("================= start of split cell %d =====================", p));
        for (int i = 0; i < nCritical; i++) {
//      logger.info(format("parent index %d nCrit %d",p,nCritical));
            int octX = 0;
            int octY = 0;
            int octZ = 0;

            if (particles[i].getX() > cells.get(p).getX()) {
                octX = 1;
            }
            if (particles[i].getY() > cells.get(p).getY()) {
                octY = 1;
            }
            if (particles[i].getZ() > cells.get(p).getZ()) {
                octZ = 1;
            }

            // Find particle's octant - should be an integer from 0 to 7
            int octant = octX + (octY << 1) + (octZ << 2);
//      logger.info(format("Child octant %d",octant));

            // If there's not a child cell in the particle's octant, create one
            boolean noChildInOctant = BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << octant));
            if (!noChildInOctant) {
//        logger.info(format("No child at octant %d for parent %d",octant,p));
                addChild(octant, p);
            }

            // Reallocate the particle in the child cell
            int c = cells.get(p).getChildAtIndex(octant);
//      logger.info(format("c %d c leaves %d",c, cells.get(c).getNumLeaves()));
            cells.get(c).setLeaf(cells.get(c).getNumLeaves(), cells.get(p).getLeavesValueAtIndex(i));
            cells.get(c).setNumLeaves(cells.get(c).getNumLeaves() + 1);
            logger.info(format(">>>>particle %d is reallocated to cell %d", cells.get(p).getLeavesValueAtIndex(i), c));

            // Check if child cell reaches nCritical - split recursively if so
            if (cells.get(c).getNumLeaves() >= nCritical) {
                splitCell(c);
            }
        }
        logger.info(format("================= end of split cell %d =====================", p));
    }

    /**
     * Get multipoles for all leaf cells.
     *
     * @param p      Cell index (should be 0, or root, traverse down).
     * @param leaves The list of leaves (cells).
     */
//  public void getMultipole(int p, ArrayList<Integer> leaves) {
    public void getMultipole(int p) {
        // If the current cell is not a leaf cell, traverse down
        if (cells.get(p).getNumLeaves() >= nCritical) {
            for (int c = 0; c < 8; c++) {
                if (BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << c))) {
//          getMultipole(cells.get(p).getChildAtIndex(c), leaves);
                    getMultipole(cells.get(p).getChildAtIndex(c));
                }
            }
        } else { // Otherwise, cell p is a leaf cell
            // Loop in leaf particles, do P2M
            for (int i = 0; i < cells.get(p).getNumLeaves(); i++) {
                int l = cells.get(p).getLeavesValueAtIndex(i);
                Atom[] atom = new Atom[]{particles[l]};

                // Calculate Multipole and fill array
                double[] calculatedMultipole = computeMomentsGeometric(atom, cells.get(p));
//        logger.info(format("Multipole for cell %d : ch %4.4f ",p,calculatedMultipole[0]));
//        logger.info(format("                        dx %4.4f dy %4.4f dz %4.4f",calculatedMultipole[1],calculatedMultipole[2],calculatedMultipole[3]));
//        logger.info(format("                        dxx %4.4f dyy %4.4f dzz %4.4f",calculatedMultipole[4],calculatedMultipole[5],calculatedMultipole[6]));
//        logger.info(format("                        dxy %4.4f dxz %4.4f dyz %4.4f",calculatedMultipole[7],calculatedMultipole[8],calculatedMultipole[9]));

                // Set Multipole
                cells.get(p).addToMultipole(calculatedMultipole);
            }
        }
    }

    /**
     * M2M.
     *
     * @param p Cell index p.
     * @param c Cell index c.
     */
    private void M2M(int p, int c) {
//    logger.info(format("M2M for p = %d and c = %d",p,c));
        double dx = cells.get(c).getX() - cells.get(p).getX();
        double dy = cells.get(c).getY() - cells.get(p).getY();
        double dz = cells.get(c).getZ() - cells.get(p).getZ();

        double[] Dxyz = new double[]{dx, dy, dz};
        double[] Dyzx = new double[]{dy, dz, dx};

        double[] currentChildMultipole = cells.get(c).getMultipole();
        double[] currentChildTracelessQDP = cells.get(c).getTracelessQDP();

        cells.get(p).addToMultipole(currentChildMultipole);
        cells.get(p).addToTracelessQuadrapole(currentChildTracelessQDP);

        // Additional Multipole Terms
        double[] additionalMultipoleTerms = new double[10];
        // Added to charge
        additionalMultipoleTerms[0] = 0;
        // Added to Dipole
        additionalMultipoleTerms[1] = currentChildMultipole[0] * Dxyz[0];
        additionalMultipoleTerms[2] = currentChildMultipole[0] * Dxyz[1];
        additionalMultipoleTerms[3] = currentChildMultipole[0] * Dxyz[2];

        // Added to Quadropole
        additionalMultipoleTerms[4] =
                2.0 * currentChildMultipole[1] * Dxyz[0] + currentChildMultipole[0] * Math.pow(Dxyz[0], 2);
        additionalMultipoleTerms[5] =
                2.0 * currentChildMultipole[2] * Dxyz[1] + currentChildMultipole[0] * Math.pow(Dxyz[1], 2);
        additionalMultipoleTerms[6] =
                2.0 * currentChildMultipole[3] * Dxyz[2] + currentChildMultipole[0] * Math.pow(Dxyz[2], 2);

        additionalMultipoleTerms[7] =
                currentChildMultipole[2] * Dxyz[0]
                        + currentChildMultipole[1] * Dyzx[0]
                        + currentChildMultipole[0] * Dxyz[0] * Dyzx[0];
        additionalMultipoleTerms[8] =
                currentChildMultipole[1] * Dxyz[2]
                        + currentChildMultipole[3] * Dyzx[2]
                        + currentChildMultipole[0] * Dxyz[2] * Dyzx[2];
        additionalMultipoleTerms[9] =
                currentChildMultipole[3] * Dxyz[1]
                        + currentChildMultipole[2] * Dyzx[1]
                        + currentChildMultipole[0] * Dxyz[1] * Dyzx[1];

        var qave = (additionalMultipoleTerms[4] + additionalMultipoleTerms[5] + additionalMultipoleTerms[6]) / 3.0;
        var xx = 1.5 * (additionalMultipoleTerms[4] - qave);
        var xy = 1.5 * additionalMultipoleTerms[7];
        var xz = 1.5 * additionalMultipoleTerms[8];
        var yx = 1.5 * additionalMultipoleTerms[7];
        var yy = 1.5 * (additionalMultipoleTerms[5] - qave);
        var yz = 1.5 * additionalMultipoleTerms[9];
        var zx = 1.5 * additionalMultipoleTerms[8];
        var zy = 1.5 * additionalMultipoleTerms[9];
        var zz = 1.5 * (additionalMultipoleTerms[6] - qave);

        // Add the traceless atomic quadrupoles to total quadrupole.
//        for (Atom atom : cellAtoms) {
//            int i = atom.getIndex() - 1;
//            double[] globalMultipolei = globalMultipole[0][i];
//            var qixx = globalMultipolei[t200];
//            var qiyy = globalMultipolei[t020];
//            var qizz = globalMultipolei[t002];
//            var qixy = globalMultipolei[t110];
//            var qixz = globalMultipolei[t101];
//            var qiyz = globalMultipolei[t011];
//            xx += qixx;
//            xy += qixy;
//            xz += qixz;
//            yx += qixy;
//            yy += qiyy;
//            yz += qiyz;
//            zx += qixz;
//            zy += qiyz;
//            zz += qizz;
//        }
        double[] additionalQDPTerms = {xx, xy, xz, yx, yy, yz, zx, zy, zz};

        cells.get(p).addToMultipole(additionalMultipoleTerms);
        cells.get(p).addToTracelessQuadrapole(additionalQDPTerms);
    }

    /**
     * Evaluate potential at one target
     *
     * @param p Index of parent cell
     * @param i Index of target particle
     */
    private void evalAtTarget(int p, int i) {
        if (p == 0) {
            logger.info(format("******Evaluating potential for atom %d********", i));
//      logger.info(format("Parent cell for atom %d is %d",i)); atoms/particles do not store this info
        }

        // Non-leaf cell
//        logger.info(format("Cell %d number of leaves = %d", p, cells.get(p).getNumLeaves()));
        if (cells.get(p).getNumLeaves() >= nCritical) {

            // Loop through p's child cells (8 octants)
            for (int oct = 0; oct < 8; oct++) {
                if (BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << oct))) {
                    int c = cells.get(p).getChildAtIndex(oct);
//          double r = distance(particles[i].getXYZ(null),cells.get(c));
                    double r = distance(particles[i].getXYZ(null), new double[]{cells.get(c).getX(), cells.get(c).getY(), cells.get(c).getZ()});

                    // Near field child cell
                    if (cells.get(c).getR() > theta * r) {
//                        logger.info(format("Cell %d is NEAR field for atom %d", c, i));
                        evalAtTarget(c, i);
                    } else { // Far field child cell
//                        logger.info(format("Cell %d is FAR field for atom %d", c, i));
                        interactionCount[i] += 1;
                        interactionWithCell(i, c);

                        double dx = particles[i].getX() - cells.get(c).getX();
                        double dy = particles[i].getY() - cells.get(c).getY();
                        double dz = particles[i].getZ() - cells.get(c).getZ();
                        double r3 = Math.pow(r, 3);
                        double r5 = r3 * Math.pow(r, 2);

                        // Calculate the weight from each multipole
                        double[] weight = new double[10];
                        weight[0] = 1 / r;
                        weight[1] = -dx / r3;
                        weight[2] = -dy / r3;
                        weight[3] = -dz / r3;
                        weight[4] = (3 * Math.pow(dx, 2)) / r5 - (1 / r3);
                        weight[5] = (3 * Math.pow(dy, 2)) / r5 - (1 / r3);
                        weight[6] = (3 * Math.pow(dz, 2)) / r5 - (1 / r3);
                        weight[7] = 3 * dx * dy / r5;
                        weight[8] = 3 * dx * dz / r5;
                        weight[9] = 3 * dy * dz / r5;

                        // Calculate dot product of multipole array and weight array
                        double dotProduct = 0.0;
                        double[] multipoleArray = cells.get(c).getMultipole();
                        for (int d = 0; d < weight.length; d++) {
                            dotProduct = dotProduct + multipoleArray[d] * weight[d];
                        }
                        phi[i] += dotProduct;
//            logger.info(format("Atom %d far field dot product = %4.3f. phi = %4.3f",i,dotProduct,phi[i]));
//            particles[i].addToPhi(dotProduct);
                    }
                }
            }
        } else { // Leaf Cell
//            logger.info(format("Cell %d is a LEAF cell for atom %d", p, i));
            // Loop in twig cell's particles
            for (int j = 0; j < cells.get(p).getNumLeaves(); j++) {
//            OctreeParticle source = particles.get(cells.get(p).getLeavesValueAtIndex(j));
                interactionCount[i] += 1;
                Atom source = particles[cells.get(p).getLeavesValueAtIndex(j)];
                double r = distance(particles[i].getXYZ(null), source.getXYZ(null));
//                if (r != 0) {
//          logger.info(format("Atom %d to atom %d r = %4.3f",i,cells.get(p).getLeavesValueAtIndex(j),r));
//              particles.get(i).addToPhi(source.getCharge() / r);
                phi[i] += source.getCharge() / r;
                //TODO: fix logic
//                if (i <= cells.get(p).getLeavesValueAtIndex(j)) { //if statment to restrict atom interactions happening x2
////                    logger.info(format("Interacting atom %d with %d", i, cells.get(p).getLeavesValueAtIndex(j)));
//                    interaction(i, cells.get(p).getLeavesValueAtIndex(j));
//                }
                if (interactionTracker[i][cells.get(p).getLeavesValueAtIndex(j)] == 0) {
                    interaction(i, cells.get(p).getLeavesValueAtIndex(j));
                    interactionTracker[i][cells.get(p).getLeavesValueAtIndex(j)] = 1;
                    interactionTracker[cells.get(p).getLeavesValueAtIndex(j)][i] = 1;
                }
//          logger.info(format("Atom %d twig particle %d contribution = %4.3f. phi = %4.3f",i,source.getIndex()-1,source.getCharge() / r,phi[i]));
//                }
            }
        }
    }

//    private void evalAtTargetFF(int p, int i) {
////        if (p == 0) {
////            logger.info(format("******Evaluating potential for atom %d********", i));
//////      logger.info(format("Parent cell for atom %d is %d",i)); atoms/particles do not store this info
////        }
//
//        // Non-leaf cell
////        logger.info(format("Cell %d number of leaves = %d", p, cells.get(p).getNumLeaves()));
//        if (cells.get(p).getNumLeaves() >= nCritical) {
//
//            // Loop through p's child cells (8 octants)
//            for (int oct = 0; oct < 8; oct++) {
//                if (BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << oct))) {
//                    int c = cells.get(p).getChildAtIndex(oct);
////          double r = distance(particles[i].getXYZ(null),cells.get(c));
//                    double r = distance(particles[i].getXYZ(null), new double[]{cells.get(c).getX(), cells.get(c).getY(), cells.get(c).getZ()});
//
//                    // Near field child cell
//                    if (cells.get(c).getR() > theta * r) {
////                        logger.info(format("Cell %d is NEAR field for atom %d", c, i));
//                        evalAtTargetFF(c, i);
//                    } else { // Far field child cell
////                        logger.info(format("Cell %d is FAR field for atom %d", c, i));
//                        interactionCount[i] += 1;
//                        interactionWithCell(i, c);
//
//                        double dx = particles[i].getX() - cells.get(c).getX();
//                        double dy = particles[i].getY() - cells.get(c).getY();
//                        double dz = particles[i].getZ() - cells.get(c).getZ();
//
//                        double r3 = Math.pow(r, 3);
//                        double r5 = r3 * Math.pow(r, 2);
//
//                        // Calculate the weight from each multipole
//                        double[] weight = new double[10];
//                        weight[0] = 1 / r;
//                        weight[1] = -dx / r3;
//                        weight[2] = -dy / r3;
//                        weight[3] = -dz / r3;
//                        weight[4] = (3 * Math.pow(dx, 2)) / r5 - (1 / r3);
//                        weight[5] = (3 * Math.pow(dy, 2)) / r5 - (1 / r3);
//                        weight[6] = (3 * Math.pow(dz, 2)) / r5 - (1 / r3);
//                        weight[7] = 3 * dx * dy / r5;
//                        weight[8] = 3 * dx * dz / r5;
//                        weight[9] = 3 * dy * dz / r5;
//
//                        // Calculate dot product of multipole array and weight array
//                        double dotProduct = 0.0;
//                        double[] multipoleArray = cells.get(c).getMultipole();
//                        for (int d = 0; d < weight.length; d++) {
//                            dotProduct = dotProduct + multipoleArray[d] * weight[d];
//                        }
//                        phi[i] += dotProduct;
////            logger.info(format("Atom %d far field dot product = %4.3f. phi = %4.3f",i,dotProduct,phi[i]));
////            particles[i].addToPhi(dotProduct);
//                    }
//                }
//            }
//        } else { // Leaf Cell
////            logger.info(format("Cell %d is a LEAF cell for atom %d", p, i));
//            // Loop in twig cell's particles
//            for (int j = 0; j < cells.get(p).getNumLeaves(); j++) {
////            OctreeParticle source = particles.get(cells.get(p).getLeavesValueAtIndex(j));
////                interactionCount[i] += 1;
//                Atom source = particles[cells.get(p).getLeavesValueAtIndex(j)];
//                double r = distance(particles[i].getXYZ(null), source.getXYZ(null));
////                if (r != 0) {
////          logger.info(format("Atom %d to atom %d r = %4.3f",i,cells.get(p).getLeavesValueAtIndex(j),r));
////              particles.get(i).addToPhi(source.getCharge() / r);
//                phi[i] += source.getCharge() / r;
////                if (i <= cells.get(p).getLeavesValueAtIndex(j)) {
//////                    logger.info(format("Interacting atom %d with %d", i, cells.get(p).getLeavesValueAtIndex(j)));
////                    interaction(i, cells.get(p).getLeavesValueAtIndex(j));
////                }
////          logger.info(format("Atom %d twig particle %d contribution = %4.3f. phi = %4.3f",i,source.getIndex()-1,source.getCharge() / r,phi[i]));
////                }
//            }
//        }
//    }


    private void calculateTotalSelfEnergy(int p, int pPermanent) {
//        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
//        GeneralizedKirkwood generalizedKirkwood = forceFieldEnergy.getGK();
//        AtomicDoubleArray selfEnergies = generalizedKirkwood.getSelfEnergy();

        if (cells.get(p).getNumLeaves() >= nCritical) {
            for (int oct = 0; oct < 8; oct++) {
                if (BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << oct))) {
                    int c = cells.get(p).getChildAtIndex(oct);
                    calculateTotalSelfEnergy(c, pPermanent);
                }
            }
        } else {
            for (int j = 0; j < cells.get(p).getNumLeaves(); j++) {
                Atom currentAtom = particles[cells.get(p).getLeavesValueAtIndex(j)];
                double currentCharge = currentAtom.getCharge();

                cells.get(pPermanent).addToTotalSelfEnergy(Math.pow(currentCharge, 2) / selfEnergies[cells.get(p).getLeavesValueAtIndex(j)]);
//                cells.get(pPermanent).addToTotalSelfEnergy(selfEnergies[cells.get(p).getLeavesValueAtIndex(j)]);
//                logger.info(format("Atom %d self energy = %4.4f and charge = %4.4f",
//                        cells.get(p).getLeavesValueAtIndex(j),selfEnergies[cells.get(p).getLeavesValueAtIndex(j)],currentCharge));
            }
        }
    }

    private void determineMaxBornRadius(int p, int pPermanent) {
        ForceFieldEnergy forceFieldEnergy = molecularAssembly.getPotentialEnergy();
        GeneralizedKirkwood generalizedKirkwood = forceFieldEnergy.getGK();
        double[] bornRadii = generalizedKirkwood.getBorn();

        if (cells.get(p).getNumLeaves() >= nCritical) {
            for (int oct = 0; oct < 8; oct++) {
                if (BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << oct))) {
                    int c = cells.get(p).getChildAtIndex(oct);
                    determineMaxBornRadius(c, pPermanent);
                }
            }
        } else {
            for (int j = 0; j < cells.get(p).getNumLeaves(); j++) {
                Atom currentAtom = particles[cells.get(p).getLeavesValueAtIndex(j)];
                if (cells.get(pPermanent).getMaxBornRadius() < bornRadii[currentAtom.getIndex()-1]) {
                    cells.get(pPermanent).setMaxBornRadius(bornRadii[currentAtom.getIndex()-1]);
                }
            }
        }
    }

    /**
     * Calculation of Near-Field and Far-Field boxes in each hierarchy level
     *
     * @param ws      Index of parent cell
     * @param nLevels Index of target particle
     */

    public void neighborCount(int ws) {
        double H = Math.ceil((Math.log(Math.ceil(cells.get(0).getR() * 2 / smallestSide[1]))) / Math.log(2.0));

        int nLevels = (int) H;
        logger.info(format("Total number of hierarchy levels in the system (H) = %d", nLevels));

        // Initialize counters
        double nNFtotal = 0.0;
        double nFFtotal = 0.0;

        logger.info(format("ws = %d", ws));
        logger.info("hLevel    nNFlevel    nNFtotal    nFFlevel    nFFtotal");

        // Loop over hierarchy levels
        for (int hLevel = 1; hLevel <= nLevels; hLevel++) {
            int nBoxes = (int) Math.pow(2, hLevel);
            int nScan = nBoxes / 2;

            // Loop over each cube x,y,z in the hierarchy
            // work on a single octant and then multiply the total count by 8
            double nNFlevel = 0.0;
            double nFFlevel = 0.0;
            for (int x = 1; x <= nScan; x++) {
                for (int y = 1; y <= nScan; y++) {
                    for (int z = 1; z <= nScan; z++) {
                        // NF count
                        double superBoxNF = (Math.min(x + ws, nBoxes) - Math.max(x - ws, 1) + 1)
                                * (Math.min(y + ws, nBoxes) - Math.max(y - ws, 1) + 1)
                                * (Math.min(z + ws, nBoxes) - Math.max(z - ws, 1) + 1);
                        nNFlevel = nNFlevel + superBoxNF - 1.0;

                        // FF count
                        double superBoxFF = (Math.min((Math.floor((x + 1) / 2) + ws) * 2, nBoxes) - Math.max((Math.floor((x - 1) / 2) - ws) * 2, 0))
                                * (Math.min((Math.floor((y + 1) / 2) + ws) * 2, nBoxes) - Math.max((Math.floor((y - 1) / 2) - ws) * 2, 0))
                                * (Math.min((Math.floor((z + 1) / 2) + ws) * 2, nBoxes) - Math.max((Math.floor((z - 1) / 2) - ws) * 2, 0));
                        nFFlevel = nFFlevel + superBoxFF - superBoxNF;
                    }
                }
            }
            nNFtotal = nNFtotal + nNFlevel;
            nFFtotal = nFFtotal + nFFlevel;

            // Print NF and FF counts for each hierarchy level
            logger.info(format("%d    %4.1f    %4.1f    %4.1f    %4.1f", hLevel, nNFlevel * 8, nNFtotal * 8, nFFlevel * 8, nFFtotal * 8));
        }
    }

    public void printCells() {
        for (OctreeCell cell : cells) {
            logger.info(format("Cell parent index %d nLeaves %d children %d %d %d %d %d %d %d %d \n",
                    cell.getParentIndex(), cell.getNumLeaves(), cell.getChildAtIndex(0), cell.getChildAtIndex(1),
                    cell.getChildAtIndex(2), cell.getChildAtIndex(3), cell.getChildAtIndex(4),
                    cell.getChildAtIndex(5), cell.getChildAtIndex(6), cell.getChildAtIndex(7)));
        }
    }

    public void printPhi() {
        double totalPhi = 0.0;
        for (int i = 0; i < phi.length; i++) {
            totalPhi += phi[i];
        }
        logger.info(format("Total potential = %4.3f", totalPhi));
    }
}