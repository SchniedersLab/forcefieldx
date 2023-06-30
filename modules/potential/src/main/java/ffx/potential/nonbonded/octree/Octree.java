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
import java.util.List;
import java.util.logging.Logger;
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

import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.NeighborList;
import org.apache.commons.lang3.BooleanUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

/**
 * Octree method presented in the Fast Multipole Method (FMM) tutorial from the Barba Group:
 * https://github.com/barbagroup/FMM_tutorial
 */
public class Octree {

  private static final Logger logger = Logger.getLogger(Octree.class.getName());
  /** List of all leaf cells */
  private final ArrayList<OctreeCell> leaves = new ArrayList<>();
  /**
   * Critical (maximum allowed) number of points allowed in any one cell: If a cell already contains
   * nCritical points, it needs to be split
   */
  private final int nCritical;
  /** List of particles (atoms) */
  private final Atom[] particles;
  /** Tolerance parameter */
  private final double theta;

  private final double[][][] globalMultipole;
  /** List of cells */
  private ArrayList<OctreeCell> cells = new ArrayList<>();

  private double[] phi = new double[304];

  private double[] smallestSide = new double[2];

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
   * @param theta Specifies near field vs far field
   */
  public Octree(int nCritical, Atom[] particles, double theta, double[][][] globalMultipole) {
    this.nCritical = nCritical;
    this.particles = particles;
    this.theta = theta;
    this.globalMultipole = globalMultipole;
  }

  /**
   * Add a child.
   *
   * @param octant The octant.
   * @param p Cell index p.
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
    logger.info(format("+++cell %d is created as a child of cell %d",c,p));
    logger.info(format("** Cell %d sidelength = %4.3f",c,cells.get(c).getR()*2));
    if (smallestSide[1] > cells.get(c).getR()*2) {
      smallestSide[0] = c;
      smallestSide[1] = cells.get(c).getR()*2;
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
    smallestSide[1] = cells.get(0).getR()*2;

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
      logger.info(format("particle %d stored in cell %d",i,current));

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
   * @param cell NeighborList Cell containing cell index information and atom index information
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

//    // Convert the quadrupole from traced to traceless form.
//    var qave = (xxqdp + yyqdp + zzqdp) / 3.0;
//    xxqdp = 1.5 * (xxqdp - qave);
//    xyqdp = 1.5 * xyqdp;
//    xzqdp = 1.5 * xzqdp;
//    yxqdp = 1.5 * yxqdp;
//    yyqdp = 1.5 * (yyqdp - qave);
//    yzqdp = 1.5 * yzqdp;
//    zxqdp = 1.5 * zxqdp;
//    zyqdp = 1.5 * zyqdp;
//    zzqdp = 1.5 * (zzqdp - qave);
//
//    // Add the traceless atomic quadrupoles to total quadrupole.
//    for (Atom atom : cellAtoms) {
//      int i = atom.getIndex() - 1;
//      double[] globalMultipolei = globalMultipole[0][i];
//      var qixx = globalMultipolei[t200];
//      var qiyy = globalMultipolei[t020];
//      var qizz = globalMultipolei[t002];
//      var qixy = globalMultipolei[t110];
//      var qixz = globalMultipolei[t101];
//      var qiyz = globalMultipolei[t011];
//      xxqdp += qixx;
//      xyqdp += qixy;
//      xzqdp += qixz;
//      yxqdp += qixy;
//      yyqdp += qiyy;
//      yzqdp += qiyz;
//      zxqdp += qixz;
//      zyqdp += qiyz;
//      zzqdp += qizz;
//    }
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

    return new double[]{netchg,xdpl,ydpl,zdpl,xxqdp,yyqdp,zzqdp,xyqdp,xzqdp,yzqdp};
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

  /** Evaluate potential at all target points */
  public void evalPotential() {
    for (int i = 0; i < particles.length; i++) {
      evalAtTarget(0, i);
//      logger.info(format("Atom %d potential = %4.3f",i,phi[i]));
    }
  }

  /**
   * Compute the L2 error.
   *
   * @param phiDirect Potential from direct summation.
   * @param phiTree Potential from the tree.
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

  /** Update sweep. */
  public void upwardSweep() {
    logger.info(format("******************Start of the Upward Sweep**********************"));
    for (int c = (cells.size() - 1); c > 0; c--) {
      int p = cells.get(c).getParentIndex();
      M2M(p, c);
      double[] calculatedMultipole = cells.get(p).getMultipole();
      logger.info(format("Geometric center of cell p%d: %4.3f %4.3f %4.3f",p,cells.get(p).getX(),cells.get(p).getY(),cells.get(p).getZ()));
      logger.info(format("Multipole for parent cell %d : ch %4.4f ",p,calculatedMultipole[0]));
      logger.info(format("                               dx %4.4f dy %4.4f dz %4.4f",calculatedMultipole[1],calculatedMultipole[2],calculatedMultipole[3]));
      logger.info(format("                               dxx %4.4f dyy %4.4f dzz %4.4f",calculatedMultipole[4],calculatedMultipole[5],calculatedMultipole[6]));
      logger.info(format("                               dxy %4.4f dxz %4.4f dyz %4.4f",calculatedMultipole[7],calculatedMultipole[8],calculatedMultipole[9]));
    }
    if (cells.size() == 1){
      double[] calculatedMultipole = cells.get(0).getMultipole();
//      logger.info(format("Geometric center of cell p%d c%d : %4.3f %4.3f %4.3f",p,c,cells.get(c).getX(),cells.get(c).getY(),cells.get(c).getZ()));
      logger.info(format("Multipole for parent cell 0 : ch %4.4f ",calculatedMultipole[0]));
      logger.info(format("                               dx %4.4f dy %4.4f dz %4.4f",calculatedMultipole[1],calculatedMultipole[2],calculatedMultipole[3]));
      logger.info(format("                               dxx %4.4f dyy %4.4f dzz %4.4f",calculatedMultipole[4],calculatedMultipole[5],calculatedMultipole[6]));
      logger.info(format("                               dxy %4.4f dxz %4.4f dyz %4.4f",calculatedMultipole[7],calculatedMultipole[8],calculatedMultipole[9]));
    }
  }

  /**
   * Spit a cell.
   *
   * @param p Cell index.
   */
  private void splitCell(int p) {
    logger.info(format("================= start of split cell %d =====================",p));
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
      logger.info(format(">>>>particle %d is reallocated to cell %d",cells.get(p).getLeavesValueAtIndex(i),c));

      // Check if child cell reaches nCritical - split recursively if so
      if (cells.get(c).getNumLeaves() >= nCritical) {
        splitCell(c);
      }
    }
    logger.info(format("================= end of split cell %d =====================",p));
  }

  /**
   * Get multipoles for all leaf cells.
   *
   * @param p Cell index (should be 0, or root, traverse down).
   * @param leaves The list of leaves (cells).
   */
//  public void getMultipole(int p, ArrayList<Integer> leaves) {
  public void getMultipole(int p){
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
        double[] calculatedMultipole = computeMomentsGeometric(atom,cells.get(p));
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

    double[] Dxyz = new double[] {dx, dy, dz};
    double[] Dyzx = new double[] {dy, dz, dx};

    double[] currentChildMultipole = cells.get(c).getMultipole();

    cells.get(p).addToMultipole(currentChildMultipole);

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

    cells.get(p).addToMultipole(additionalMultipoleTerms);
  }

  /**
   * Evaluate potential at one target
   *
   * @param p Index of parent cell
   * @param i Index of target particle
   */
  private void evalAtTarget(int p, int i) {

    // Non-leaf cell
    if (cells.get(p).getNumLeaves() >= nCritical) {

      // Loop through p's child cells (8 octants)
      for (int oct = 0; oct < 8; oct++) {
        if (BooleanUtils.toBoolean(cells.get(p).getnChild() & (1 << oct))) {
          int c = cells.get(p).getChildAtIndex(oct);
//          double r = distance(particles[i].getXYZ(null),cells.get(c));
          double r = distance(particles[i].getXYZ(null), new double[]{cells.get(c).getX(), cells.get(c).getY(), cells.get(c).getZ()});

          // Near field child cell
          if (cells.get(c).getR() > theta * r) {
            evalAtTarget(c, i);
          } else { // Far field child cell
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
            phi[particles[i].getIndex()-1] += dotProduct;
//            logger.info(format("Atom %d far field dot product = %4.3f. phi = %4.3f",i,dotProduct,phi[i]));
//            particles[i].addToPhi(dotProduct);
          }
        }
      }
    } else { // Leaf Cell
          // Loop in twig cell's particles
      for (int j = 0; j < cells.get(p).getNumLeaves(); j++) {
//            OctreeParticle source = particles.get(cells.get(p).getLeavesValueAtIndex(j));
        Atom source = particles[cells.get(p).getLeavesValueAtIndex(j)];
        double r = distance(particles[i].getXYZ(null),source.getXYZ(null));
        if (r != 0) {
//          logger.info(format("Atom %d to atom %d r = %4.3f",i,cells.get(p).getLeavesValueAtIndex(j),r));
//              particles.get(i).addToPhi(source.getCharge() / r);
          phi[particles[i].getIndex()-1] += source.getCharge() / r;
//          logger.info(format("Atom %d twig particle %d contribution = %4.3f. phi = %4.3f",i,source.getIndex()-1,source.getCharge() / r,phi[i]));
        }
      }
    }
  }

  /**
   * Calculation of Near-Field and Far-Field boxes in each hierarchy level
   *
   * @param ws Index of parent cell
   * @param nLevels Index of target particle
   */

  public void neighborCount(int ws){
    // TODO: Add H
    double H = Math.ceil((Math.log(Math.ceil(cells.get(0).getR()*2/smallestSide[1])))/Math.log(2.0));
    logger.info(format("Cell count %4.1f Smallest sidelength %4.3f",smallestSide[0],smallestSide[1]));

    int nLevels = 2;

    // Initialize counters
    double nNFtotal = 0;
    double nFFtotal = 0;

    logger.info(format("ws = %d",ws));
    logger.info("hLevel    nNFlevel    nNFtotal    nFFlevel    nFFtotal");

    // Loop over hierarchy levels
    for (int hLevel = 1; hLevel <= nLevels; hLevel++) {
      double nBoxes = Math.pow(2,hLevel);
      double nScan = nBoxes / 2;

      // Loop over each cube x,y,z in the hierarchy
      // work on a single octant and then multiply the total count by 8
      double nNFlevel = 0;
      double nFFlevel = 0;
      for (int x = 0; x < nScan; x++) {
        for (int y = 0; y < nScan; y++) {
          for (int z = 0; z < nScan; z++) {
            // NF count
            double superBoxNF = (Math.min(x+ws,nBoxes) - Math.max(x-ws,1) + 1)
                    * (Math.min(y+ws,nBoxes) - Math.max(y-ws,1) + 1)
                    * (Math.min(z+ws,nBoxes) - Math.max(z-ws,1) + 1);
            nNFlevel = nNFlevel + superBoxNF - 1;

            // FF count
            double superBoxFF = (Math.min((Math.floor((x+1)/2)+ws)*2, nBoxes) - Math.max((Math.floor((x-1)/2)-ws)*2, 0))
                    * (Math.min((Math.floor((y+1)/2)+ws)*2, nBoxes) - Math.max((Math.floor((y-1)/2)-ws)*2, 0))
                    * (Math.min((Math.floor((z+1)/2)+ws)*2, nBoxes) - Math.max((Math.floor((z-1)/2)-ws)*2, 0));
            nFFlevel = nFFlevel + superBoxFF - superBoxNF;
          }
        }
      }
      nNFtotal = nNFtotal + nNFlevel;
      nFFtotal = nFFtotal + nFFlevel;

      // Print NF and FF counts for each hierarchy level
      logger.info(format("%d    %4.1f    %4.1f    %4.1f    %4.1f",hLevel, nNFlevel*8, nNFtotal*8, nFFlevel*8, nFFtotal*8));
    }
  }

  public void printCells(){
    for(OctreeCell cell : cells){
      logger.info(format("Cell parent index %d nLeaves %d children %d %d %d %d %d %d %d %d \n",
              cell.getParentIndex(), cell.getNumLeaves(), cell.getChildAtIndex(0), cell.getChildAtIndex(1),
              cell.getChildAtIndex(2), cell.getChildAtIndex(3), cell.getChildAtIndex(4),
              cell.getChildAtIndex(5), cell.getChildAtIndex(6), cell.getChildAtIndex(7)));
    }
  }
}
