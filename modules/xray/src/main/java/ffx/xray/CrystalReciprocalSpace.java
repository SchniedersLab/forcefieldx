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
package ffx.xray;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static org.apache.commons.math.util.FastMath.exp;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;

import ffx.crystal.*;
import ffx.numerics.ComplexNumber;
import ffx.numerics.fft.Complex;
import ffx.numerics.fft.Complex3DParallel;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.SpatialDensityLoop;
import ffx.potential.nonbonded.SpatialDensityRegion;
import ffx.xray.RefinementMinimize.RefinementMode;

import static ffx.numerics.fft.Complex3D.iComplex3D;

/**
 * Structure factor calculation (including bulk solvent structure factors)
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @see <a href="http://dx.doi.org/10.1107/S0567739473000458" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1973). A29, 183-191.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0567739477001211" target="_blank">
 * L. F. Ten Eyck, Acta Cryst. (1977). A33, 486-492.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0365110X55001862" target="_blank">
 * J. Waser, Acta Cryst. (1955). 8, 595.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0108767388009183" target="_blank">
 * A. T. Brunger, Acta Cryst. (1989). A45, 42-50.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/97809553602060000551"
 * target="_blank"> G. Bricogne, Int. Tables Cryst. (2006). Vol. B, ch. 1.3, pp.
 * 25-98.</a>
 *
 * @see <a href="http://dx.doi.org/10.1002/jcc.1032" target="_blank"> J. A.
 * Grant, B. T. Pickup, A. Nicholls, J. Comp. Chem. (2001). 22, 608-640</a>
 *
 * @see <a href="http://dx.doi.org/10.1006/jmbi.1994.1633" target="_blank"> J.
 * S. Jiang, A. T. Brunger, JMB (1994) 243, 100-115.</a>
 *
 * @see <a href="http://dx.doi.org/10.1107/S0907444910031045" target="_blank">
 * T.D. Fenn, M. J. Schnieders, A. T. Brunger, Acta Cryst. (2010). D66,
 * 1024-1031.</a>
 * @author Michael J. Schnieders
 *
 */
public class CrystalReciprocalSpace {

    /**
     * the possible solvent model methods
     */
    public static interface SolventModel {

        /**
         * do not model solvent scattering
         */
        public static final int NONE = 1;
        /**
         * the classical binary (0, 1) model
         */
        public static final int BINARY = 2;
        /**
         * a Gaussian switch model
         */
        public static final int GAUSSIAN = 3;
        /**
         * a cubic polynomial switch model (default)
         */
        public static final int POLYNOMIAL = 4;
    }
    private static final Logger logger = Logger.getLogger(CrystalReciprocalSpace.class.getName());
    private static double toSeconds = 0.000000001;
    private final double badd;
    private double arad;
    private int aradgrid;
    private final Crystal crystal;
    private final Resolution resolution;
    private final ReflectionList reflectionlist;
    private boolean solvent = false;
    protected int solventmodel;
    private final boolean neutron;
    private boolean use_3g = true;
    private double weight = 1.0;
    protected double solvent_a;
    protected double solvent_b;
    private final int nSymm;
    private final Atom atoms[];
    private final int nAtoms;
    // not final for purposes of finite differences
    private double coordinates[][][];
    private final FormFactor atomffactors[];
    private final FormFactor solventffactors[];
    private final int fftX, fftY, fftZ;
    private final double fftscale;
    private final int complexFFT3DSpace;
    private final int halfFFTX, halfFFTY, halfFFTZ;
    protected final double densityGrid[];
    protected final double solventGrid[];
    private final ParallelTeam parallelTeam;
    private final int threadCount;
    private final SpatialDensityRegion spatialDensityRegion;
    private final SpatialDensityRegion bulkSpatialDensityRegion;
    private final AtomicDensityLoop atomicDensityLoops[];
    private final AtomicGradientRegion atomicGradientRegion;
    private final SolventDensityLoop solventDensityLoops[];
    private final SolventDensityLoop bulkSolventDensityLoops[];
    private final SolventGradientRegion solventGradientRegion;
    private final ParallelTeam fftTeam;
    private final Complex3DParallel complexFFT3D;
    protected boolean lambdaTerm = false;
    protected double lambda = 1.0;

    /**
     * Crystal Reciprocal Space constructor, assumes this is not a bulk solvent
     * mask and is not a neutron data set
     *
     * @param reflectionlist the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam) {
        this(reflectionlist, atoms, fftTeam, parallelTeam, false, false,
                SolventModel.POLYNOMIAL);
    }

    /**
     * Crystal Reciprocal Space constructor, assumes this is not a neutron data
     * set and implements a polynomial bulk solvent mask if needed
     *
     * @param reflectionlist the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventmask true if this is a bulk solvent mask
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventmask) {
        this(reflectionlist, atoms, fftTeam, parallelTeam, solventmask, false,
                SolventModel.POLYNOMIAL);
    }

    /**
     * Crystal Reciprocal Space constructor, assumes a polynomial bulk solvent
     * mask if needed
     *
     * @param reflectionlist the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventmask true if this is a bulk solvent mask
     * @param neutron true if this is a neutron structure
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventmask, boolean neutron) {
        this(reflectionlist, atoms, fftTeam, parallelTeam, solventmask, neutron,
                SolventModel.POLYNOMIAL);
    }

    /**
     * Crystal Reciprocal Space constructor, all parameters provided
     *
     * @param reflectionlist the {@link ReflectionList} to fill with structure
     * factors
     * @param atoms array of {@link ffx.potential.bonded.Atom atoms} for
     * structure factor computation
     * @param fftTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param parallelTeam {@link edu.rit.pj.ParallelTeam} for parallelization
     * @param solventmask true if this is a bulk solvent mask
     * @param neutron true if this is a neutron structure
     * @param solventmodel bulk solvent model type
     * @see CrystalReciprocalSpace.SolventModel
     */
    public CrystalReciprocalSpace(ReflectionList reflectionlist,
            Atom atoms[],
            ParallelTeam fftTeam, ParallelTeam parallelTeam,
            boolean solventmask, boolean neutron, int solventmodel) {
        this.reflectionlist = reflectionlist;
        this.crystal = reflectionlist.crystal;
        this.resolution = reflectionlist.resolution;
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.fftTeam = fftTeam;
        this.parallelTeam = parallelTeam;
        this.solvent = solventmask;
        this.solventmodel = solventmodel;
        this.neutron = neutron;
        this.nSymm = 1;
        threadCount = parallelTeam.getThreadCount();
        // necssary for the bulksolvent expansion!
        int bulknsym = crystal.spaceGroup.symOps.size();

        double density = 2.0 / resolution.sampling_limit();
        double res = resolution.res_limit();

        // Set default FFT grid size from unit cell dimensions.
        int nX = (int) Math.floor(crystal.a * density / res) + 1;
        int nY = (int) Math.floor(crystal.b * density / res) + 1;
        int nZ = (int) Math.floor(crystal.c * density / res) + 1;
        if (nX % 2 != 0) {
            nX += 1;
        }
        // Use preferred dimensions.
        while (!Complex.preferredDimension(nX)) {
            nX += 2;
        }
        while (!Complex.preferredDimension(nY)) {
            nY += 1;
        }
        while (!Complex.preferredDimension(nZ)) {
            nZ += 1;
        }

        fftX = nX;
        fftY = nY;
        fftZ = nZ;
        halfFFTX = nX / 2;
        halfFFTY = nY / 2;
        halfFFTZ = nZ / 2;
        fftscale = crystal.volume;
        complexFFT3DSpace = fftX * fftY * fftZ * 2;
        densityGrid = new double[complexFFT3DSpace];

        String solventname = "none";
        if (solvent) {
            badd = 0.0;
            atomffactors = null;
            double vdwr;
            switch (solventmodel) {
                case SolventModel.BINARY:
                    solventname = "binary";
                    solvent_a = 1.0;
                    solvent_b = 1.0;
                    solventGrid = new double[complexFFT3DSpace];
                    solventffactors = new SolventBinaryFormFactor[nAtoms];
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventffactors[i] = new SolventBinaryFormFactor(atoms[i], vdwr + solvent_a);
                    }
                    break;
                case SolventModel.GAUSSIAN:
                    solventname = "Gaussian";
                    solvent_a = 11.5;
                    solvent_b = 0.55;
                    solventGrid = new double[complexFFT3DSpace];
                    solventffactors = new SolventGaussFormFactor[nAtoms];
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventffactors[i] = new SolventGaussFormFactor(atoms[i], vdwr * solvent_b);
                    }
                    break;
                case SolventModel.POLYNOMIAL:
                    solventname = "polynomial switch";
                    solvent_a = 0.0;
                    solvent_b = 0.8;
                    solventGrid = new double[complexFFT3DSpace];
                    solventffactors = new SolventPolyFormFactor[nAtoms];
                    for (int i = 0; i < nAtoms; i++) {
                        vdwr = atoms[i].getVDWType().radius * 0.5;
                        solventffactors[i] = new SolventPolyFormFactor(atoms[i], vdwr + solvent_a, solvent_b);
                    }
                    break;
                default:
                    solventGrid = null;
                    solventffactors = null;
                    break;
            }
        } else {
            badd = 2.0;
            if (neutron) {
                atomffactors = new NeutronFormFactor[nAtoms];
                for (int i = 0; i < nAtoms; i++) {
                    atomffactors[i] = new NeutronFormFactor(atoms[i], badd);
                }
            } else {
                atomffactors = new XRayFormFactor[nAtoms];
                for (int i = 0; i < nAtoms; i++) {
                    atomffactors[i] = new XRayFormFactor(atoms[i], use_3g, badd);
                }
            }
            solventGrid = null;
            solventffactors = null;
        }

        // determine number of grid points to sample density on
        arad = -1.0;
        for (Atom a : atoms) {
            double vdwr = a.getVDWType().radius * 0.5;
            if (!solvent) {
                arad = Math.max(arad, a.getFormFactorWidth());
            } else {
                switch (solventmodel) {
                    case SolventModel.BINARY:
                        arad = Math.max(arad, vdwr + solvent_a + 0.2);
                        break;
                    case SolventModel.GAUSSIAN:
                        arad = Math.max(arad, vdwr * solvent_b + 2.0);
                        break;
                    case SolventModel.POLYNOMIAL:
                        arad = Math.max(arad, vdwr + solvent_b + 0.2);
                        break;
                }
            }
        }
        aradgrid = (int) Math.floor(arad * nX / crystal.a) + 1;

        // local copy of coordinates - note use of bulknsym
        coordinates = new double[bulknsym][3][nAtoms];
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        for (int i = 0; i < bulknsym; i++) {
            double x[] = coordinates[i][0];
            double y[] = coordinates[i][1];
            double z[] = coordinates[i][2];
            for (int j = 0; j < nAtoms; j++) {
                Atom aj = atoms[j];
                crystal.applySymOp(aj.getXYZ(), xyz, symops.get(i));
                x[j] = xyz[0];
                y[j] = xyz[1];
                z[j] = xyz[2];
            }
        }

        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder();
            if (solvent) {
                sb.append(String.format("  Bulk Solvent Grid\n"));
                sb.append(String.format("  Bulk solvent model type:       %s\n",
                        solventname));
            } else {
                sb.append(String.format("  Atomic Grid\n"));
            }
            sb.append(String.format("  Form factor grid radius (radius):  %d (%8.3f)\n",
                    aradgrid, arad));
            sb.append(String.format("  Grid density:               %8.3f\n",
                    density));
            sb.append(String.format("  Grid dimensions:           (%d,%d,%d)\n",
                    fftX, fftY, fftZ));
            logger.info(sb.toString());
        }

        if (solvent) {
            int minWork = nSymm;

            if (solventmodel != SolventModel.NONE) {
                spatialDensityRegion =
                        new SpatialDensityRegion(fftX, fftY, fftZ,
                        densityGrid, (aradgrid + 2) * 2, nSymm, minWork,
                        threadCount, crystal, atoms, coordinates);
                bulkSpatialDensityRegion =
                        new BulkSolventDensityRegion(fftX, fftY, fftZ,
                        solventGrid, (aradgrid + 2) * 2, bulknsym, minWork,
                        threadCount, crystal, atoms, coordinates, 4.0, parallelTeam);
                if (solventmodel == SolventModel.GAUSSIAN) {
                    spatialDensityRegion.setInitValue(0.0);
                } else {
                    spatialDensityRegion.setInitValue(1.0);
                }

                atomicDensityLoops = null;
                solventDensityLoops = new SolventDensityLoop[threadCount];
                for (int i = 0; i < threadCount; i++) {
                    solventDensityLoops[i] =
                            new SolventDensityLoop(spatialDensityRegion);
                }
                spatialDensityRegion.setDensityLoop(solventDensityLoops);

                bulkSolventDensityLoops = new SolventDensityLoop[threadCount];
                for (int i = 0; i < threadCount; i++) {
                    bulkSolventDensityLoops[i] =
                            new SolventDensityLoop(bulkSpatialDensityRegion);
                }
                bulkSpatialDensityRegion.setDensityLoop(bulkSolventDensityLoops);

                atomicGradientRegion = null;
                solventGradientRegion = new SolventGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES);
            } else {
                spatialDensityRegion = null;
                bulkSpatialDensityRegion = null;
                atomicDensityLoops = null;
                solventDensityLoops = null;
                bulkSolventDensityLoops = null;
                atomicGradientRegion = null;
                solventGradientRegion = null;
            }

        } else {
            /**
             * Create nSymm pieces of work per thread; the empty pieces will be
             * removed leaving 1 piece of work per thread.
             */
            int minWork = nSymm;
            spatialDensityRegion = new SpatialDensityRegion(fftX, fftY, fftZ,
                    densityGrid, (aradgrid + 2) * 2, nSymm, minWork,
                    threadCount, crystal, atoms, coordinates);
            bulkSpatialDensityRegion = null;
            spatialDensityRegion.setInitValue(0.0);

            atomicDensityLoops = new AtomicDensityLoop[threadCount];
            solventDensityLoops = null;
            bulkSolventDensityLoops = null;
            for (int i = 0; i < threadCount; i++) {
                atomicDensityLoops[i] =
                        new AtomicDensityLoop(spatialDensityRegion);
            }
            spatialDensityRegion.setDensityLoop(atomicDensityLoops);

            atomicGradientRegion = new AtomicGradientRegion(RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES);
            solventGradientRegion = null;
        }
        complexFFT3D = new Complex3DParallel(fftX, fftY, fftZ, fftTeam);
    }

    /**
     * Set the current value of the state variable.
     *
     * @param lambda a double.
     */
    protected void setLambda(double lambda) {
        this.lambda = lambda;
    }

    /**
     * set the bulk solvent parameters
     *
     * @param a added atom width (binary, polynomial only)
     * @param b falloff of the switch (Gaussian, polynomial only)
     */
    public void setSolventAB(double a, double b) {
        double vdwr;
        this.solvent_a = a;
        this.solvent_b = b;
        if (!solvent) {
            return;
        }
        switch (solventmodel) {
            case SolventModel.BINARY:
                for (int i = 0; i < nAtoms; i++) {
                    vdwr = atoms[i].getVDWType().radius * 0.5;
                    solventffactors[i] = new SolventBinaryFormFactor(atoms[i], vdwr + solvent_a);
                }
                break;
            case SolventModel.GAUSSIAN:
                for (int i = 0; i < nAtoms; i++) {
                    vdwr = atoms[i].getVDWType().radius * 0.5;
                    solventffactors[i] = new SolventGaussFormFactor(atoms[i], vdwr * solvent_b);
                }
                break;
            case SolventModel.POLYNOMIAL:
                for (int i = 0; i < nAtoms; i++) {
                    vdwr = atoms[i].getVDWType().radius * 0.5;
                    solventffactors[i] = new SolventPolyFormFactor(atoms[i], vdwr + solvent_a, solvent_b);
                }
                break;
        }
    }

    /**
     * should the structure factor computation use 3 Gaussians or 6 for atoms?
     *
     * @param use_3g if true, use 3 Gaussian
     */
    public void setUse3G(boolean use_3g) {
        this.use_3g = use_3g;
        if (solvent || neutron) {
            return;
        }
        for (int i = 0; i < nAtoms; i++) {
            atomffactors[i] = new XRayFormFactor(atoms[i], use_3g, badd);
        }
    }

    /**
     * return dataset weight
     *
     * @return the weight wA
     */
    public double getWeight() {
        return weight;
    }

    /**
     * set the dataset weight
     *
     * @param weight desired weight wA
     */
    public void setWeight(double weight) {
        this.weight = weight;
    }

    /**
     * offset X coordinates (mostly for finite difference checks)
     *
     * @param n {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    public void deltaX(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[0] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][0][n] = symxyz[0];
        }
    }

    /**
     * offset Y coordinates (mostly for finite difference checks)
     *
     * @param n {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    public void deltaY(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[1] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][1][n] = symxyz[1];
        }
    }

    /**
     * offset Z coordinates (mostly for finite difference checks)
     *
     * @param n {@link ffx.potential.bonded.Atom} to apply delta to
     * @param delta amount to shift atom by
     */
    public void deltaZ(int n, double delta) {
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];
        atoms[n].getXYZ(xyz);
        xyz[2] += delta;
        for (int i = 0; i < nSymm; i++) {
            crystal.applySymOp(xyz, symxyz, symops.get(i));
            coordinates[i][2][n] = symxyz[2];
        }
    }

    /**
     * set atom coordinates
     *
     * @param coords new coordinate positions (3 params per atom)
     */
    public void setCoordinates(double coords[]) {
        assert (coords != null);
        List<SymOp> symops = crystal.spaceGroup.symOps;
        double xyz[] = new double[3];
        double symxyz[] = new double[3];

        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            xyz[0] = coords[index++];
            xyz[1] = coords[index++];
            xyz[2] = coords[index++];
            for (int j = 0; j < nSymm; j++) {
                crystal.applySymOp(xyz, symxyz, symops.get(j));
                coordinates[j][0][i] = symxyz[0];
                coordinates[j][1][i] = symxyz[1];
                coordinates[j][2][i] = symxyz[2];
            }
        }
    }

    /**
     * parallelized computation of structure factors
     *
     * @param hkldata structure factor list to fill in
     * @see DiffractionRefinementData
     */
    public void computeDensity(double hkldata[][]) {
        computeDensity(hkldata, false);
    }

    /**
     * parallelized computation of structure factors
     *
     * @param hkldata structure factor list to fill in
     * @param print if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    public void computeDensity(double hkldata[][], boolean print) {
        if (solvent) {
            if (solventmodel != SolventModel.NONE) {
                computeSolventDensity(hkldata, print);
            }
        } else {
            computeAtomicDensity(hkldata, print);
        }
    }

    /**
     * compute inverse FFT to determine atomic gradients
     *
     * @param hkldata structure factors to apply inverse FFT
     * @param freer array of free r flags corresponding to hkldata
     * @param flag Rfree flag value
     * @param refinementmode
     * {@link RefinementMinimize.RefinementMode refinement mode}
     * @see RefinementMinimize.RefinementMode
     * @see DiffractionRefinementData
     */
    public void computeAtomicGradients(double hkldata[][],
            int freer[], int flag, RefinementMode refinementmode) {
        computeAtomicGradients(hkldata, freer, flag, refinementmode, false);
    }

    /**
     * compute inverse FFT to determine atomic gradients
     *
     * @param hkldata structure factors to apply inverse FFT
     * @param freer array of free r flags corresponding to hkldata
     * @param flag Rfree flag value
     * @param refinementmode
     * {@link RefinementMinimize.RefinementMode refinement mode}
     * @param print if true, print information on timings during the calculation
     * @see RefinementMinimize.RefinementMode
     * @see DiffractionRefinementData
     */
    public void computeAtomicGradients(double hkldata[][],
            int freer[], int flag, RefinementMode refinementmode,
            boolean print) {

        if (solvent && solventmodel == SolventModel.NONE) {
            return;
        }

        // zero out the density
        for (int i = 0; i < complexFFT3DSpace; i++) {
            densityGrid[i] = 0.0;
        }

        int nfree = 0;
        StringBuilder sb = new StringBuilder();
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        // int nsym = 1;
        List<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        ComplexNumber cj = new ComplexNumber();
        HKL ij = new HKL();
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = hkldata[ih.index()];
            if (Double.isNaN(fc[0])) {
                continue;
            }

            // cross validation check!!!
            if (freer != null) {
                if (freer[ih.index()] == flag) {
                    nfree++;
                    continue;
                }
            }

            c.re(fc[0]);
            c.im(fc[1]);
            // scale
            c.times_ip(2.0 / fftscale);

            // apply symmetry
            for (int j = 0; j < nsym; j++) {
                cj.copy(c);
                crystal.applyTransSymRot(ih, ij, symops.get(j));
                double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                int h = Crystal.mod(ij.h(), fftX);
                int k = Crystal.mod(ij.k(), fftY);
                int l = Crystal.mod(ij.l(), fftZ);

                if (h < halfFFTX + 1) {
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    cj.phase_shift_ip(shift);
                    densityGrid[ii] += cj.re();
                    densityGrid[ii + 1] += -cj.im();
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    cj.phase_shift_ip(shift);
                    densityGrid[ii] += cj.re();
                    densityGrid[ii + 1] += cj.im();
                }
            }
        }
        symtime += System.nanoTime();

        long startTime = System.nanoTime();
        complexFFT3D.ifft(densityGrid);
        long fftTime = System.nanoTime() - startTime;


        /*
         CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
         mapout.write(densityGrid);
         */

        startTime = System.nanoTime();
        long permanentDensityTime = 0;
        try {
            if (solvent) {
                solventGradientRegion.setRefinementMode(refinementmode);
                parallelTeam.execute(solventGradientRegion);
            } else {
                atomicGradientRegion.setRefinementMode(refinementmode);
                parallelTeam.execute(atomicGradientRegion);
            }
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Exception computing atomic gradients.";
            logger.log(Level.SEVERE, message, e);
        }

        if (solvent) {
            sb.append(String.format(" Solvent symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(String.format(" Solvent inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(String.format(" Grid solvent gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(String.format(" %d reflections ignored (cross validation set)\n", nfree));
        } else {
            sb.append(String.format(" Atomic symmetry: %8.3f\n", symtime * toSeconds));
            sb.append(String.format(" Atomic inverse FFT: %8.3f\n", fftTime * toSeconds));
            sb.append(String.format(" Grid atomic gradients: %8.3f\n", permanentDensityTime * toSeconds));
            sb.append(String.format(" %d reflections ignored (cross validation set)\n", nfree));
        }

        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    /**
     * parallelized computation of structure factors (identical to
     * compuateDensity)
     *
     * @param hkldata structure factor list to fill in
     * @see DiffractionRefinementData
     */
    public void computeAtomicDensity(double hkldata[][]) {
        computeAtomicDensity(hkldata, false);
    }

    /**
     * parallelized computation of structure factors (identical to
     * computeDensity)
     *
     * @param hkldata structure factor list to fill in
     * @param print if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    public void computeAtomicDensity(double hkldata[][], boolean print) {
        StringBuilder sb = new StringBuilder();
        // clear out the reflection data
        int n = reflectionlist.hkllist.size();
        for (int i = 0; i < n; i++) {
            hkldata[i][0] = hkldata[i][1] = 0.0;
        }

        spatialDensityRegion.assignAtomsToCells();
        long permanentDensityTime = 0;
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Fatal exception evaluating atomic electron density.";
            logger.log(Level.SEVERE, message, e);
        }

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(densityGrid);

        long startTime = System.nanoTime();
        complexFFT3D.fft(densityGrid);
        long fftTime = System.nanoTime() - startTime;

        // extract structure factors
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        List<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        HKL ij = new HKL();
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = hkldata[ih.index()];

            // apply symmetry
            for (int j = 0; j < nsym; j++) {
                crystal.applyTransSymRot(ih, ij, symops.get(j));
                double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                int h = Crystal.mod(ij.h(), fftX);
                int k = Crystal.mod(ij.k(), fftY);
                int l = Crystal.mod(ij.l(), fftZ);

                if (h < halfFFTX + 1) {
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    c.re(densityGrid[ii]);
                    c.im(densityGrid[ii + 1]);
                    c.phase_shift_ip(shift);
                    fc[0] += c.re();
                    fc[1] += c.im();
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    c.re(densityGrid[ii]);
                    c.im(-densityGrid[ii + 1]);
                    c.phase_shift_ip(shift);
                    fc[0] += c.re();
                    fc[1] += c.im();
                }
            }
        }

        // scale
        double scale = (fftscale) / (fftX * fftY * fftZ);
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = hkldata[ih.index()];
            c.re(fc[0]);
            c.im(fc[1]);
            // remove Badd
            double s = Crystal.invressq(crystal, ih);
            c.times_ip(scale * exp(0.25 * badd * s));
            c.conjugate_ip();

            fc[0] = c.re();
            fc[1] = c.im();
        }
        symtime += System.nanoTime();

        sb.append(String.format(" Atomic grid density: %8.3f\n", permanentDensityTime * toSeconds));
        sb.append(String.format(" Atomic FFT: %8.3f\n", fftTime * toSeconds));
        sb.append(String.format(" Atomic symmetry/scaling: %8.3f\n", symtime * toSeconds));

        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    /**
     * parallelized computation of bulk solvent structure factors
     *
     * @param hkldata structure factor list to fill in
     * @see DiffractionRefinementData
     */
    public void computeSolventDensity(double hkldata[][]) {
        computeSolventDensity(hkldata, false);
    }

    /**
     * parallelized computation of bulk solvent structure factors
     *
     * @param hkldata structure factor list to fill in
     * @param print if true, print information on timings during the calculation
     * @see DiffractionRefinementData
     */
    public void computeSolventDensity(double hkldata[][], boolean print) {
        StringBuilder sb = new StringBuilder();
        // clear out the reflection data
        int n = reflectionlist.hkllist.size();
        for (int i = 0; i < n; i++) {
            hkldata[i][0] = hkldata[i][1] = 0.0;
        }

        spatialDensityRegion.assignAtomsToCells();
        long permanentDensityTime = 0;
        try {
            long startTime = System.nanoTime();
            parallelTeam.execute(spatialDensityRegion);
            permanentDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }

        long startTime = System.nanoTime();
        if (solventmodel == SolventModel.BINARY) {
            /*
             * need to shrink mask
             * first, copy densitygrid to df_map
             * temporarily to store original map
             */
            int nmap = densityGrid.length;
            System.arraycopy(densityGrid, 0, solventGrid, 0, nmap);

            // Logic to loop within the cutoff box.
            double xyz[] = {solvent_b, solvent_b, solvent_b};
            double uvw[] = new double[3];
            crystal.toFractionalCoordinates(xyz, uvw);
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;

            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY);

                        if (densityGrid[ii] == 1.0) {
                            continue;
                        }

                        for (int iz = k - ifrz; iz <= k + ifrz; iz++) {
                            int giz = Crystal.mod(iz, fftZ);
                            for (int iy = j - ifry; iy <= j + ifry; iy++) {
                                int giy = Crystal.mod(iy, fftY);
                                for (int ix = i - ifrx; ix <= i + ifrx; ix++) {
                                    int gix = Crystal.mod(ix, fftX);

                                    final int jj = iComplex3D(gix, giy, giz, fftX, fftY);

                                    if (densityGrid[jj] == 1.0) {
                                        solventGrid[ii] = 1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // copy back
            System.arraycopy(solventGrid, 0, densityGrid, 0, nmap);
        }

        // copy stuff over for derivatives
        System.arraycopy(densityGrid, 0, solventGrid, 0, densityGrid.length);

        // babinet principle
        for (int k = 0; k < fftZ; k++) {
            for (int j = 0; j < fftY; j++) {
                for (int i = 0; i < fftX; i++) {
                    final int ii = iComplex3D(i, j, k, fftX, fftY);
                    if (solventmodel == SolventModel.BINARY
                            || solventmodel == SolventModel.POLYNOMIAL) {
                        densityGrid[ii] = 1.0 - densityGrid[ii];
                    } else if (solventmodel == SolventModel.GAUSSIAN) {
                        densityGrid[ii] = 1.0 - exp(-solvent_a * densityGrid[ii]);
                    }
                }
            }
        }
        long expTime = System.nanoTime() - startTime;

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(densityGrid);

        startTime = System.nanoTime();
        complexFFT3D.fft(densityGrid);
        long fftTime = System.nanoTime() - startTime;

        // extract structure factors
        long symtime = -System.nanoTime();
        int nsym = crystal.spaceGroup.symOps.size();
        List<SymOp> symops = crystal.spaceGroup.symOps;
        ComplexNumber c = new ComplexNumber();
        HKL ij = new HKL();
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = hkldata[ih.index()];

            // apply symmetry
            for (int j = 0; j < nsym; j++) {
                crystal.applyTransSymRot(ih, ij, symops.get(j));
                double shift = Crystal.sym_phase_shift(ih, symops.get(j));

                int h = Crystal.mod(ij.h(), fftX);
                int k = Crystal.mod(ij.k(), fftY);
                int l = Crystal.mod(ij.l(), fftZ);

                if (h < halfFFTX + 1) {
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    c.re(densityGrid[ii]);
                    c.im(densityGrid[ii + 1]);
                    c.phase_shift_ip(shift);
                    fc[0] += c.re();
                    fc[1] += c.im();
                } else {
                    h = (fftX - h) % fftX;
                    k = (fftY - k) % fftY;
                    l = (fftZ - l) % fftZ;
                    final int ii = iComplex3D(h, k, l, fftX, fftY);
                    c.re(densityGrid[ii]);
                    c.im(-densityGrid[ii + 1]);
                    c.phase_shift_ip(shift);
                    fc[0] += c.re();
                    fc[1] += c.im();
                }
            }
        }

        // scale
        double scale = (fftscale) / (fftX * fftY * fftZ);
        for (HKL ih : reflectionlist.hkllist) {
            double fc[] = hkldata[ih.index()];
            c.re(fc[0]);
            c.im(fc[1]);
            c.times_ip(scale);
            c.conjugate_ip();

            // negative: babinet
            fc[0] = -c.re();
            fc[1] = -c.im();
        }
        symtime += System.nanoTime();

        // now expand derivative mask so it includes nearby local density
        bulkSpatialDensityRegion.assignAtomsToCells();
        long solventDensityTime = 0;
        try {
            startTime = System.nanoTime();
            parallelTeam.execute(bulkSpatialDensityRegion);
            solventDensityTime = System.nanoTime() - startTime;
        } catch (Exception e) {
            String message = "Fatal exception evaluating solvent electron density.";
            logger.log(Level.SEVERE, message, e);
        }

        long solventExpTime = 0;
        if (solventmodel == SolventModel.GAUSSIAN) {
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        final int ii = iComplex3D(i, j, k, fftX, fftY);
                        solventGrid[ii] = exp(-solvent_a * solventGrid[ii]);
                    }
                }
            }
        }
        solventExpTime = System.nanoTime() - startTime;

        // CCP4MapWriter mapout = new CCP4MapWriter(fftX, fftY, fftZ, crystal, "/tmp/foo.map");
        // mapout.write(solventGrid);

        sb.append(String.format(" Solvent grid density: %8.3f\n", permanentDensityTime * toSeconds));
        sb.append(String.format(" bulk solvent exponentiation: %8.3f\n", expTime * toSeconds));
        sb.append(String.format(" Solvent FFT: %8.3f\n", fftTime * toSeconds));
        sb.append(String.format(" Solvent symmetry/scaling: %8.3f\n", symtime * toSeconds));
        sb.append(String.format(" Solvent grid expansion: %8.3f\n", solventDensityTime * toSeconds));
        sb.append(String.format(" Solvent grid expansion exponentiation: %8.3f\n", solventExpTime * toSeconds));

        if (logger.isLoggable(Level.INFO) && print) {
            logger.info(sb.toString());
        }
    }

    private class AtomicDensityLoop extends SpatialDensityLoop {

        final double xyz[] = new double[3];
        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public AtomicDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            final double lambdai = atoms[n].applyLambda() ? lambda : 1.0;
            xyz[0] = coordinates[iSymm][0][n];
            xyz[1] = coordinates[iSymm][1][n];
            xyz[2] = coordinates[iSymm][2][n];
            FormFactor atomff = atomffactors[n];
            atomff.update(xyz, badd);
            crystal.toFractionalCoordinates(xyz, uvw);
            final int frad = Math.min(aradgrid,
                    (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                int gix = Crystal.mod(ix, fftX);
                xf[0] = ix / (double) fftX;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy / (double) fftY;
                    for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                        int giz = Crystal.mod(iz, fftZ);
                        xf[2] = iz / (double) fftZ;
                        crystal.toCartesianCoordinates(xf, xc);

                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = atomff.rho(grid[ii], lambdai, xc);
                    }
                }
            }
        }
    }

    private class SolventDensityLoop extends SpatialDensityLoop {

        final double xyz[] = new double[3];
        final double uvw[] = new double[3];
        final double xc[] = new double[3];
        final double xf[] = new double[3];
        final double grid[];

        public SolventDensityLoop(SpatialDensityRegion region) {
            super(region, region.getNsymm(), region.actualCount);
            grid = region.getGrid();
        }

        @Override
        public void gridDensity(int iSymm, int n) {
            final double lambdai = atoms[n].applyLambda() ? lambda : 1.0;
            xyz[0] = coordinates[iSymm][0][n];
            xyz[1] = coordinates[iSymm][1][n];
            xyz[2] = coordinates[iSymm][2][n];
            FormFactor solventff = solventffactors[n];
            solventff.update(xyz);
            crystal.toFractionalCoordinates(xyz, uvw);
            double vdwr = atoms[n].getVDWType().radius * 0.5;
            int frad = aradgrid;
            switch (solventmodel) {
                case SolventModel.BINARY:
                    frad = Math.min(aradgrid,
                            (int) Math.floor((vdwr + solvent_a + 0.2) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.GAUSSIAN:
                    frad = Math.min(aradgrid,
                            (int) Math.floor((vdwr * solvent_b + 2.0) * fftX / crystal.a) + 1);
                    break;
                case SolventModel.POLYNOMIAL:
                    frad = Math.min(aradgrid,
                            (int) Math.floor((vdwr + solvent_b + 0.2) * fftX / crystal.a) + 1);
                    break;
            }

            // Logic to loop within the cutoff box.
            final double frx = fftX * uvw[0];
            final int ifrx = (int) frx;
            final int ifrxu = ifrx + frad;

            final double fry = fftY * uvw[1];
            final int ifry = (int) fry;
            final int ifryu = ifry + frad;

            final double frz = fftZ * uvw[2];
            final int ifrz = (int) frz;
            final int ifrzu = ifrz + frad;

            for (int ix = ifrx - frad; ix <= ifrxu; ix++) {
                int gix = Crystal.mod(ix, fftX);
                xf[0] = ix / (double) fftX;
                for (int iy = ifry - frad; iy <= ifryu; iy++) {
                    int giy = Crystal.mod(iy, fftY);
                    xf[1] = iy / (double) fftY;
                    for (int iz = ifrz - frad; iz <= ifrzu; iz++) {
                        int giz = Crystal.mod(iz, fftZ);
                        xf[2] = iz / (double) fftZ;
                        crystal.toCartesianCoordinates(xf, xc);

                        final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                        grid[ii] = solventff.rho(grid[ii], lambdai, xc);
                    }
                }
            }
        }
    }

    private class AtomicGradientRegion extends ParallelRegion {

        private final AtomicGradientLoop atomicGradientLoop[];
        private RefinementMode refinementmode;

        public AtomicGradientRegion(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
            atomicGradientLoop = new AtomicGradientLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                atomicGradientLoop[i] = new AtomicGradientLoop();
            }
        }

        public RefinementMode getRefinementMode() {
            return refinementmode;
        }

        public void setRefinementMode(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, atomicGradientLoop[getThreadIndex()]);
            } catch (Exception e) {
                e.printStackTrace();
                logger.severe(e.toString());
            }
        }

        private class AtomicGradientLoop extends IntegerForLoop {

            final double xyz[] = new double[3];
            final double uvw[] = new double[3];
            final double xc[] = new double[3];
            final double xf[] = new double[3];
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    if (lambdaTerm) {
                        if (!atoms[n].applyLambda()) {
                            continue;
                        }
                    }
                    xyz[0] = coordinates[0][0][n];
                    xyz[1] = coordinates[0][1][n];
                    xyz[2] = coordinates[0][2][n];
                    FormFactor atomff = atomffactors[n];
                    atomff.update(xyz, 0.0);
                    crystal.toFractionalCoordinates(xyz, uvw);
                    final int dfrad = Math.min(aradgrid,
                            (int) Math.floor(atoms[n].getFormFactorWidth() * fftX / crystal.a) + 1);

                    // Logic to loop within the cutoff box.
                    final double frx = fftX * uvw[0];
                    final int ifrx = (int) frx;

                    final double fry = fftY * uvw[1];
                    final int ifry = (int) fry;

                    final double frz = fftZ * uvw[2];
                    final int ifrz = (int) frz;

                    for (int ix = ifrx - dfrad; ix <= ifrx + dfrad; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        for (int iy = ifry - dfrad; iy <= ifry + dfrad; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy / (double) fftY;
                            for (int iz = ifrz - dfrad; iz <= ifrz + dfrad; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz / (double) fftZ;
                                crystal.toCartesianCoordinates(xf, xc);

                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                atomff.rho_grad(xc, weight * densityGrid[ii], refinementmode);
                            }
                        }
                    }
                }
            }
        }
    }

    private class SolventGradientRegion extends ParallelRegion {

        private final SolventGradientLoop solventGradientLoop[];
        private RefinementMode refinementmode;

        public SolventGradientRegion(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
            solventGradientLoop = new SolventGradientLoop[threadCount];
            for (int i = 0; i < threadCount; i++) {
                solventGradientLoop[i] = new SolventGradientLoop();
            }
        }

        public RefinementMode getRefinementMode() {
            return refinementmode;
        }

        public void setRefinementMode(RefinementMode refinementmode) {
            this.refinementmode = refinementmode;
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, solventGradientLoop[getThreadIndex()]);
            } catch (Exception e) {
                logger.severe(e.toString());
            }
        }

        private class SolventGradientLoop extends IntegerForLoop {

            double xyz[] = new double[3];
            double uvw[] = new double[3];
            double xc[] = new double[3];
            double xf[] = new double[3];
            // Extra padding to avert cache interference.
            long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(final int lb, final int ub) {
                for (int n = lb; n <= ub; n++) {
                    if (lambdaTerm) {
                        if (!atoms[n].applyLambda()) {
                            continue;
                        }
                    }
                    xyz[0] = coordinates[0][0][n];
                    xyz[1] = coordinates[0][1][n];
                    xyz[2] = coordinates[0][2][n];
                    FormFactor solventff = solventffactors[n];
                    solventff.update(xyz);
                    crystal.toFractionalCoordinates(xyz, uvw);
                    double vdwr = atoms[n].getVDWType().radius * 0.5;
                    double dfcmult = 1.0;
                    int dfrad = aradgrid;
                    switch (solventmodel) {
                        case SolventModel.BINARY:
                            dfrad = Math.min(aradgrid,
                                    (int) Math.floor((vdwr + solvent_a + 0.2) * fftX / crystal.a) + 1);
                            break;
                        case SolventModel.GAUSSIAN:
                            dfrad = Math.min(aradgrid,
                                    (int) Math.floor((vdwr * solvent_b + 2.0) * fftX / crystal.a) + 1);
                            dfcmult = solvent_a;
                            break;
                        case SolventModel.POLYNOMIAL:
                            dfrad = Math.min(aradgrid,
                                    (int) Math.floor((vdwr + solvent_b + 0.2) * fftX / crystal.a) + 1);
                            break;
                    }

                    // Logic to loop within the cutoff box.
                    final double frx = fftX * uvw[0];
                    final int ifrx = (int) frx;
                    final int ifrxu = ifrx + dfrad;

                    final double fry = fftY * uvw[1];
                    final int ifry = (int) fry;
                    final int ifryu = ifry + dfrad;

                    final double frz = fftZ * uvw[2];
                    final int ifrz = (int) frz;
                    final int ifrzu = ifrz + dfrad;

                    for (int ix = ifrx - dfrad; ix <= ifrxu; ix++) {
                        int gix = Crystal.mod(ix, fftX);
                        xf[0] = ix / (double) fftX;
                        for (int iy = ifry - dfrad; iy <= ifryu; iy++) {
                            int giy = Crystal.mod(iy, fftY);
                            xf[1] = iy / (double) fftY;
                            for (int iz = ifrz - dfrad; iz <= ifrzu; iz++) {
                                int giz = Crystal.mod(iz, fftZ);
                                xf[2] = iz / (double) fftZ;
                                crystal.toCartesianCoordinates(xf, xc);

                                final int ii = iComplex3D(gix, giy, giz, fftX, fftY);
                                solventff.rho_grad(xc,
                                        weight * densityGrid[ii] * dfcmult * solventGrid[ii],
                                        refinementmode);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * <p>getXDim</p>
     *
     * @return a double.
     */
    public double getXDim() {
        return fftX;
    }

    /**
     * <p>getYDim</p>
     *
     * @return a double.
     */
    public double getYDim() {
        return fftY;
    }

    /**
     * <p>getZDim</p>
     *
     * @return a double.
     */
    public double getZDim() {
        return fftZ;
    }

    /**
     * <p>densityNorm</p>
     *
     * @param data an array of double.
     * @param meansd an array of double.
     * @param norm a boolean.
     */
    public void densityNorm(double data[], double meansd[], boolean norm) {
        double mean, sd;

        mean = sd = 0.0;
        int n = 0;
        for (int k = 0; k < fftZ; k++) {
            for (int j = 0; j < fftY; j++) {
                for (int i = 0; i < fftX; i++) {
                    int index = iComplex3D(i, j, k, fftX, fftY);
                    n++;
                    mean += (data[index] - mean) / n;
                }
            }
        }

        n = 0;
        for (int k = 0; k < fftZ; k++) {
            for (int j = 0; j < fftY; j++) {
                for (int i = 0; i < fftX; i++) {
                    int index = iComplex3D(i, j, k, fftX, fftY);
                    sd += Math.pow(data[index] - mean, 2.0);
                    n++;
                }
            }
        }
        sd = Math.sqrt(sd / n);

        if (meansd != null) {
            meansd[0] = mean;
            meansd[1] = sd;
        }

        if (norm) {
            for (int k = 0; k < fftZ; k++) {
                for (int j = 0; j < fftY; j++) {
                    for (int i = 0; i < fftX; i++) {
                        int index = iComplex3D(i, j, k, fftX, fftY);
                        data[index] = (data[index] - mean) / sd;
                    }
                }
            }
        }
    }
}
