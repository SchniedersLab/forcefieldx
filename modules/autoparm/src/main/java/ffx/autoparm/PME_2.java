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
package ffx.autoparm;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelSection;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.numerics.Potential;
import ffx.numerics.TensorRecursion;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Torsion;
import ffx.potential.nonbonded.NeighborList;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.ReciprocalSpace.FFTMethod;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.parameters.ForceField.ForceFieldType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;

import static ffx.numerics.Erf.erfc;
import static ffx.numerics.VectorMath.cross;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.dot;
import static ffx.numerics.VectorMath.norm;
import static ffx.numerics.VectorMath.r;
import static ffx.numerics.VectorMath.scalar;
import static ffx.numerics.VectorMath.sum;
import static ffx.numerics.VectorMath.transpose3;
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

/**
 * This Particle Mesh Ewald class implements PME for the AMOEBA polarizable
 * mutlipole force field in parallel using a {@link NeighborList} for any
 * {@link Crystal}. The real space contribution is contained within this class,
 * but the reciprocal space contribution is delegated to the
 * {@link ReciprocalSpace} class.
 *
 * @author Michael J. Schnieders<br> derived from:<br> TINKER code by Jay
 * Ponder, Pengyu Ren and Tom Darden.<br>
 * @see <a href="http://dx.doi.org/10.1063/1.1630791" target="_blank"> C. Sagui,
 * L. G. Pedersen, and T. A. Darden, Journal of Chemical Physics 120 (1), 73
 * (2004)</a><br> <a href="http://link.aip.org/link/?JCPSA6/98/10089/1"
 * target="_blank"> T. Darden, D. York, and L. Pedersen, Journal of Chemical
 * Physics 98 (12), 10089 (1993)</a><br> <a href="http://www.ccp5.org"
 * target="_blank"> W. Smith, "Point Multipoles in the Ewald Summation
 * (Revisited)", CCP5 Newsletter, 46, 18-30, 1998</a><br>
 *
 */
public class PME_2 implements Potential {

    /**
     * Constant <code>pedit=false</code>
     */
    public static boolean pedit = false;
    /**
     * Constant <code>propyze=false</code>
     */
    public static boolean propyze = false;
    private Boolean use_pme = false;
    private double target_grid[][][];//nSymm, nAtoms, 4
    private double scaling[];
    private ArrayList<String> key;
    public boolean fitmpl = true, fitdpl = true, fitqdpl = true;
    private int nvars;
    /**
     * Constant <code>logger</code>
     */
    public static final Logger logger = Logger.getLogger(PME_2.class.getName());
    private double totalEnergy;

    @Override
    public double energy(double[] x) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Polarization modes include "direct", in which induced dipoles do not
     * interact, and "mutual" that converges the self-consistent field to a
     * tolerance specified by the "polar-eps" keyword.
     */
    public enum Polarization {

        MUTUAL, DIRECT, NONE
    }
    private int interactions;
    private double multipoleEnergy;
    private double polarizationEnergy;
    /**
     * Constant <code>lambda=1.0</code>
     */
    public static double lambda = 1.0;
    /**
     * Reference to the force field being used.
     */
    private final ForceField forceField;
    /**
     * Unit cell and spacegroup information.
     */
    public static Crystal crystal;
    /**
     * Number of symmetry operators.
     */
    public static int nSymm;
    /**
     * An ordered array of atoms in the system.
     */
    public static Atom atoms[];
    /**
     * The number of atoms in the system.
     */
    public static int nAtoms;
    /**
     * Dimensions of [nsymm][3][nAtoms].
     */
    public static double coordinates[][][];
    /**
     * Constant <code>neighborLists=</code>
     */
    public static int neighborLists[][][];
    private final int[][][] ewaldLists;
    private final int[][] ewaldCounts;
    /**
     * *************************************************************************
     * Permanent multipole variables.
     */
    /**
     * Permanent multipoles in their local frame.
     */
    public static double localMultipole[][];
    /**
     * Constant <code>frame</code>
     */
    public static MultipoleType.MultipoleFrameDefinition frame[];
    /**
     * Constant <code>axisAtom=</code>
     */
    public static int axisAtom[][];
    /**
     * Constant <code>lol=</code>
     */
    public static int lol[];
    /**
     * Dimensions of [nsymm][nAtoms][10]
     */
    public static double globalMultipole[][][];
    private final double cartesianMultipolePhi[][];
    /**
     * The interaction energy between 1-2 multipoles is scaled by m12scale.
     */
    private final double m12scale;
    /**
     * The interaction energy between 1-3 multipoles is scaled by m13scale.
     */
    private final double m13scale;
    /**
     * The interaction energy between 1-4 multipoles is scaled by m14scale.
     */
    private final double m14scale;
    /**
     * The interaction energy between 1-5 multipoles is scaled by m15scale.
     */
    private final double m15scale;
    /**
     * *************************************************************************
     * Induced dipole variables.
     */
    /**
     * Polarization mode.
     */
    public static Polarization polarization;
    /**
     * Constant <code>polsor=</code>
     */
    public static double polsor;
    /**
     * Constant <code>poleps=</code>
     */
    public static double poleps;
    /**
     * Direct polarization field due to permanent multipoles at polarizable
     * sites within their group are scaled. The scaling is 0.0 in AMOEBA.
     */
    public static double d11scale;
    /**
     * The interaction energy between a permanent multipole and polarizable site
     * that are 1-2 is scaled by p12scale.
     */
    public static double p12scale;
    /**
     * The interaction energy between a permanent multipole and polarizable site
     * that are 1-3 is scaled by p13scale.
     */
    public static double p13scale;
    /**
     * Constant <code>pdamp=</code>
     */
    public static double pdamp[];
    /**
     * Constant <code>thole=</code>
     */
    public static double thole[];
    /**
     * Constant <code>polarizability=</code>
     */
    public static double polarizability[];
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public static double inducedDipole[][][];
    /**
     * Constant <code>inducedDipolep=</code>
     */
    public static double inducedDipolep[][][];
    /**
     * Constant <code>directDipole=</code>
     */
    public static double directDipole[][];
    /**
     * Constant <code>directDipolep=</code>
     */
    public static double directDipolep[][];
    private final double cartesianDipolePhi[][];
    private final double cartesianDipolepPhi[][];
    /**
     * Dimensions of [nsymm][nAtoms][3]
     */
    public static double field1[][];
    /**
     * Constant <code>field2=</code>
     */
    public static double field2[][];
    /**
     * Constant <code>ip11=</code>
     */
    public static int ip11[][];
    /**
     * Constant <code>ip12=</code>
     */
    public static int ip12[][];
    /**
     * Constant <code>ip13=</code>
     */
    public static int ip13[][];
    /**
     * Constant <code>ip14=</code>
     */
    public static int ip14[][]; //added by gchattree
    /**
     * *************************************************************************
     * Mutable Particle Mesh Ewald constants.
     */
    private double aewald;
    private double alsq2;
    private double piEwald;
    private double aewald3;
    private double off;
    /**
     * Constant <code>off2=</code>
     */
    public static double off2;
    private double permanentSelfEnergy;
    /**
     * *************************************************************************
     * Parallel variables.
     */
    /**
     * By default, maxThreads is set to the number of available SMP cores.
     */
    public static int maxThreads;
    /**
     * The default ParallelTeam encapsulates the maximum number of theads used
     * to parallelize the electrostatics calculation.
     */
    public static ParallelTeam parallelTeam;
    /**
     * Either 1 or 2; see description below.
     */
    private final int sectionThreads;
    /**
     * The sectionTeam encapsulates 1 or 2 threads.
     *
     * If it contains 1 thread, the real and reciprocal space calculations are
     * done sequentially.
     *
     * If it contains 2 threads, the real and reciprocal space calculations will
     * be done concurrently.
     */
    private final ParallelTeam sectionTeam;
    /**
     * If the real and reciprocal space parts of PME are done sequentially, then
     * the realSpaceTeam is equal parallalTeam.
     *
     * If the real and reciprocal space parts of PME are done concurrently, then
     * the realSpaceTeam will have fewer threads than the default parallelTeam.
     */
    private final ParallelTeam realSpaceTeam;
    /**
     * If real and reciprocal space are done sequentially or OpenCL is used,
     * then realSpaceThreads == maxThreads.
     *
     * Otherwise the number of realSpaceThreads is set to ffx.realSpaceThreads.
     */
    private final int realSpaceThreads;
    /**
     * If the real and reciprocal space parts of PME are done sequentially, then
     * the reciprocalSpaceTeam is equal parallalTeam.
     *
     * If the real and reciprocal space parts of PME are done concurrently, then
     * the reciprocalSpaceTeam will have fewer threads than the default
     * parallelTeam.
     */
    private final ParallelTeam fftTeam;
    /**
     * If real and reciprocal space are done sequentially then
     * reciprocalSpaceThreads == maxThreads
     *
     * If OpenCL is used, reciprocalSpaceThreads == 1
     *
     * Otherwise, reciprocalSpaceThreads = maxThreads - realSpaceThreads
     */
    private final int fftThreads;
    /**
     * Constant <code>rotateMultipolesRegion</code>
     */
    public static RotateMultipolesRegion rotateMultipolesRegion;
    private final ExpandCoordinatesRegion expandCoordinatesRegion;
    /**
     * Constant <code>expandInducedDipolesRegion</code>
     */
    public static ExpandInducedDipolesRegion expandInducedDipolesRegion;
    //private ReciprocalSpace reciprocalSpace;
    //private PermanentFieldRegion permanentFieldRegion;
    //private InducedDipoleFieldRegion inducedDipoleFieldRegion;
    private RealSpaceEnergyRegion realSpaceEnergyRegion;
    private TorqueRegion torqueRegion;
    /**
     * Constant <code>pairWiseSchedule</code>
     */
    public static IntegerSchedule pairWiseSchedule;
    private final SharedDoubleArray sharedGrad[];
    private final SharedDoubleArray sharedTorque[];
    private final boolean gpuFFT;
    private long realSpaceTime;
    //private long reciprocalSpaceTime;
    //private long bsplineTime, densityTime, realAndFFTTime, phiTime;
    /**
     * Constant <code>toSeconds=1.0e-9</code>
     */
    public static double toSeconds = 1.0e-9;
    /**
     * Conversion from electron^2/Ang to Kcal/mole
     */
    public static final double electric = 332.063709;
    /**
     * The sqrt of PI.
     */
    public static final double sqrtPi = sqrt(Math.PI);
    //private ArrayList<Double> put;

    /**
     * <p>
     * Constructor for PME_2.</p>
     *
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
     * @param neighborLists an array of int.
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     * @param key a {@link java.util.ArrayList} object.
     */
    public PME_2(ForceField forceField, Atom[] atoms,
            Crystal crystal, ParallelTeam parallelTeam, int neighborLists[][][], ArrayList<String> key) {
        this.forceField = forceField;
        this.atoms = atoms;
        this.crystal = crystal;
        this.parallelTeam = parallelTeam;
        this.neighborLists = neighborLists;
        this.key = key;
        nAtoms = atoms.length;
        nSymm = crystal.spaceGroup.getNumberOfSymOps();
        maxThreads = parallelTeam.getThreadCount();

        coordinates = new double[nSymm][3][nAtoms];
        inducedDipole = new double[nSymm][nAtoms][3];
        inducedDipolep = new double[nSymm][nAtoms][3];
        neighborLists = new int[nSymm][][];
        double x[] = coordinates[0][0];
        double y[] = coordinates[0][1];
        double z[] = coordinates[0][2];
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            double xyz[] = ai.getXYZ(null);
            x[i] = xyz[0];
            y[i] = xyz[1];
            z[i] = xyz[2];
        }

        /**
         * The size of reduced neighbor list and cache depend on the size of the
         * real space cutoff.
         */
        ewaldLists = new int[nSymm][nAtoms][];
        ewaldCounts = new int[nSymm][nAtoms];

        polsor = forceField.getDouble(ForceFieldDouble.POLAR_SOR, 0.70);
        m12scale = forceField.getDouble(ForceFieldDouble.MPOLE_12_SCALE, 0.0);
        m13scale = forceField.getDouble(ForceFieldDouble.MPOLE_13_SCALE, 0.0);
        m14scale = forceField.getDouble(ForceFieldDouble.MPOLE_14_SCALE, 0.4);
        m15scale = forceField.getDouble(ForceFieldDouble.MPOLE_15_SCALE, 0.8);
        d11scale = forceField.getDouble(ForceFieldDouble.DIRECT_11_SCALE, 0.0);
        p12scale = forceField.getDouble(ForceFieldDouble.POLAR_12_SCALE, 0.0);
        p13scale = forceField.getDouble(ForceFieldDouble.POLAR_13_SCALE, 0.0);
        String polar = forceField.getString(ForceFieldString.POLARIZATION, "DEFAULT");
        boolean polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);

        Polarization initPolarization;
        try {
            initPolarization = Polarization.valueOf(polar.toUpperCase());
        } catch (Exception e) {
            initPolarization = Polarization.MUTUAL;
        }

        if (!polarizationTerm) {
            initPolarization = Polarization.NONE;
        }

        polarization = initPolarization;
        poleps = forceField.getDouble(ForceFieldDouble.POLAR_EPS, 1e-6);

        String temp = forceField.getString(ForceField.ForceFieldString.FFT_METHOD, "PJ");
        FFTMethod method;
        try {
            method = ReciprocalSpace.FFTMethod.valueOf(temp.toUpperCase().trim());
        } catch (Exception e) {
            method = ReciprocalSpace.FFTMethod.PJ;
        }
        gpuFFT = method != ReciprocalSpace.FFTMethod.PJ;

        localMultipole = new double[nAtoms][10];
        frame = new MultipoleType.MultipoleFrameDefinition[nAtoms];
        axisAtom = new int[nAtoms][];
        assignMultipoles();
        globalMultipole = new double[nSymm][nAtoms][10];
        cartesianMultipolePhi = new double[nAtoms][tensorCount];
        directDipole = new double[nAtoms][3];
        directDipolep = new double[nAtoms][3];
        cartesianDipolePhi = new double[nAtoms][tensorCount];
        cartesianDipolepPhi = new double[nAtoms][tensorCount];
        field1 = new double[nAtoms][3];
        field2 = new double[nAtoms][3];
        ip11 = new int[nAtoms][];
        ip12 = new int[nAtoms][];
        ip13 = new int[nAtoms][];
        ip14 = new int[nAtoms][];//added by chattree
        if (!use_pme) {
            polargrp();
        } else {
            assignPolarizationGroups();
        }
        thole = new double[nAtoms];
        pdamp = new double[nAtoms];
        polarizability = new double[nAtoms];
        for (Atom ai : atoms) {
            PolarizeType polarizeType = ai.getPolarizeType();
            int index = ai.xyzIndex - 1;
            thole[index] = polarizeType.thole;
            pdamp[index] = polarizeType.pdamp;
            polarizability[index] = polarizeType.polarizability;
        }
        sharedGrad = new SharedDoubleArray[3];
        sharedTorque = new SharedDoubleArray[3];
        for (int i = 0; i < 3; i++) {
            sharedGrad[i] = new SharedDoubleArray(nAtoms);
            sharedTorque[i] = new SharedDoubleArray(nAtoms);
        }
        setEwaldParameters();
        if (logger.isLoggable(Level.INFO)) {
            StringBuilder sb = new StringBuilder("\n Electrostatics\n");
            sb.append(format(" Polarization:                         %8s\n", polarization.toString()));
            if (polarization == Polarization.MUTUAL) {
                sb.append(format(" SCF convergence criteria:            %8.3e\n", poleps));
                sb.append(format(" SOR parameter                         %8.3f\n", polsor));
            }
            sb.append(format(" Real space cut-off:                   %8.3f (A)\n", off));
            sb.append(format(" Ewald coefficient:                    %8.3f", aewald));
            logger.info(sb.toString());
        }

        if (gpuFFT) {
            sectionThreads = 2;
            realSpaceThreads = parallelTeam.getThreadCount();
            fftThreads = 1;
            sectionTeam = new ParallelTeam(sectionThreads);
            realSpaceTeam = parallelTeam;
            fftTeam = new ParallelTeam(fftThreads);
        } else {
            boolean concurrent;
            int realThreads = 1;
            try {
                realThreads = forceField.getInteger(ForceField.ForceFieldInteger.PME_REAL_THREADS);
                if (realThreads >= maxThreads || realThreads < 1) {
                    throw new Exception("pme-real-threads must be < ffx.nt and greater than 0");
                }
                concurrent = true;
            } catch (Exception e) {
                concurrent = false;
            }
            if (concurrent) {
                sectionThreads = 2;
                realSpaceThreads = realThreads;
                fftThreads = maxThreads - realThreads;
                sectionTeam = new ParallelTeam(sectionThreads);
                realSpaceTeam = new ParallelTeam(realSpaceThreads);
                fftTeam = new ParallelTeam(fftThreads);
            } else {
                /**
                 * If pme-real-threads is not defined, then do real and
                 * reciprocal space parts sequentially.
                 */
                sectionThreads = 1;
                realSpaceThreads = maxThreads;
                fftThreads = maxThreads;
                sectionTeam = new ParallelTeam(sectionThreads);
                realSpaceTeam = parallelTeam;
                fftTeam = parallelTeam;
            }
        }

        boolean available = false;
        String pairWiseStrategy = null;
        try {
            pairWiseStrategy = forceField.getString(ForceField.ForceFieldString.REAL_SCHEDULE);
            IntegerSchedule.parse(pairWiseStrategy);
            available = true;
        } catch (Exception e) {
            available = false;
        }
        if (available) {
            pairWiseSchedule = IntegerSchedule.parse(pairWiseStrategy);
            logger.info(format(" Electrostatics pairwise schedule: %s", pairWiseStrategy));
        } else {
            pairWiseSchedule = IntegerSchedule.fixed();
        }

        rotateMultipolesRegion = new RotateMultipolesRegion(maxThreads, false);
        expandCoordinatesRegion = new ExpandCoordinatesRegion(maxThreads);
        expandInducedDipolesRegion = new ExpandInducedDipolesRegion(maxThreads);
//        if(use_pme){
//        	/**
//        	 * Note that we always pass on the unit cell crystal to ReciprocalSpace
//        	 * instance even if the real space calculations require
//        	 * a ReplicatesCrystal.
//        	 */
//        	reciprocalSpace = new ReciprocalSpace(crystal.getUnitCell(), forceField,
//        			coordinates, atoms, aewald, fftTeam, parallelTeam);
//        	permanentFieldRegion = new PermanentFieldRegion(realSpaceTeam);
//        	inducedDipoleFieldRegion = new InducedDipoleFieldRegion(realSpaceTeam);
//        	realSpaceEnergyRegion = new RealSpaceEnergyRegion(maxThreads);
//        	torqueRegion = new TorqueRegion(maxThreads);
//        }

        Boolean gradient = true;
//        if(use_pme){
//        	multipoleEnergy = 0.0;
//            polarizationEnergy = 0.0;
//            interactions = 0;
//            realSpaceTime = 0;
//            reciprocalSpaceTime = 0;
//        }
        /**
         * Initialize the coordinates.
         */
//        double x[] = coordinates[0][0];
//        double y[] = coordinates[0][1];
//        double z[] = coordinates[0][2];
        for (int i = 0; i < nAtoms; i++) {
            double xyz2[] = atoms[i].getXYZ(null);
            x[i] = xyz2[0];
            y[i] = xyz2[1];
            z[i] = xyz2[2];
        }

        checkCacheSize();

        /**
         * Initialize the gradient accumulation arrays.
         */
        if (gradient) {
            for (int j = 0; j < nAtoms; j++) {
                sharedGrad[0].set(j, 0.0);
                sharedGrad[1].set(j, 0.0);
                sharedGrad[2].set(j, 0.0);
                sharedTorque[0].set(j, 0.0);
                sharedTorque[1].set(j, 0.0);
                sharedTorque[2].set(j, 0.0);
            }
        }
        /**
         * Expand coordinates and rotate multipoles for all atoms in the unit
         * cell.
         */
//        expandCoordinates();
//        rotateMulitpoles();
        if (!use_pme) {
            init_prms();
        }
//        else{
//            /**
//             * Find the permanent multipole potential and its gradients.
//             */
//        	try {
//                parallelTeam.execute(expandCoordinatesRegion);
//                parallelTeam.execute(rotateMultipolesRegion);
//
//                bsplineTime = -System.nanoTime();
//                reciprocalSpace.computeBSplines();
//                bsplineTime += System.nanoTime();
//
//                densityTime = -System.nanoTime();
//                reciprocalSpace.splinePermanentMultipoles(globalMultipole, null);
//                densityTime += System.nanoTime();
//
//                /**
//                 * Here the real space contribution to the field is calculated at
//                 * the same time the reciprocal space convolution is being done.
//                 * This is useful since the reciprocal space convolution
//                 * (the 3D FFT and inverse FFT) do not parallelize well.
//                 */
//                realAndFFTTime = -System.nanoTime();
//                sectionTeam.execute(permanentFieldRegion);
//                realAndFFTTime += System.nanoTime();
//
//                phiTime = -System.nanoTime();
//                reciprocalSpace.computePermanentPhi(cartesianMultipolePhi);
//                phiTime += System.nanoTime();
//            } catch (Exception e) {
//                String message = "Fatal exception computing the permanent multipole field.\n";
//                logger.log(Level.SEVERE, message, e);
//            }
//
//            /**
//             * Do the self consistent field calculation.
//             */
//            selfConsistentField(logger.isLoggable(Level.FINE));
//
//            if (logger.isLoggable(Level.FINE)) {
//                StringBuilder sb = new StringBuilder();
//                sb.append(format("\n b-Spline:   %8.3f (sec)\n", bsplineTime * toSeconds));
//                sb.append(format(" Density:    %8.3f (sec)\n", densityTime * toSeconds));
//                sb.append(format(" Real + FFT: %8.3f (sec)\n", realAndFFTTime * toSeconds));
//                sb.append(format(" Phi:        %8.3f (sec)\n", phiTime * toSeconds));
//                logger.fine(sb.toString());
//            }
//        }

        //added by gchattree, sets scaling for LGBFS and read in keywords
        String ln;
        for (int i = 0; i < key.size(); i++) {
            ln = key.get(i);
            if (ln.toUpperCase().contains("FIX-MONOPOLE")) {
                fitmpl = false;
            } else if (ln.toUpperCase().contains("FIX-DIPOLE")) {
                fitdpl = false;
            } else if (ln.toUpperCase().contains("FIX-QUADRUPOLE")) {
                fitqdpl = false;
            }
        }
        nvars = getCoordinates(null).length;
        scaling = new double[nvars];
        for (int i = 0; i < scaling.length; i++) {
            scaling[i] = 1;
        }
    }

//    @Override
//    public double getdEdL() {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    @Override
//    public double getd2EdL2() {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
//
//    @Override
//    public void getdEdXdL(double[] gradients) {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }
    /**
     * <p>
     * set_target_grid</p>
     *
     * @param target_grid an array of double.
     */
    public void set_target_grid(double target_grid[][][]) {
        this.target_grid = target_grid;
    }

    /**
     * <p>
     * varprm</p>
     *
     * @param x an array of double.
     * @param ivar a int.
     * @param eps a double.
     */
    public void varprm(double x[], int ivar, double eps) {
        int n = 0;
        ArrayList<Integer> types = new ArrayList<Integer>();
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            int itype = ai.getType();
            if (!types.contains(itype)) {
                types.add(itype);
                if (fitmpl) {
                    if (localMultipole[i][t000] != 0) {
                        localMultipole[i][t000] = x[n] + ((n == ivar) ? eps : 0);
                        n++;
                    }
                }
                if (fitdpl) {
                    if (localMultipole[i][t100] != 0) {
                        localMultipole[i][t100] = x[n] + ((n == ivar) ? eps : 0);
                        n++;
                    }
                    if (localMultipole[i][t010] != 0) {
                        localMultipole[i][t010] = x[n] + ((n == ivar) ? eps : 0);
                        n++;
                    }
                    if (localMultipole[i][t001] != 0) {
                        localMultipole[i][t001] = x[n] + ((n == ivar) ? eps : 0);
                        n++;
                    }
                }
                if (fitqdpl) {
                    if (localMultipole[i][t200] != 0) {
                        localMultipole[i][t200] = (x[n] + ((n == ivar) ? eps : 0)) * 3;
                        n++;
                    }
                    if (localMultipole[i][t020] != 0) {
                        localMultipole[i][t020] = (x[n] + ((n == ivar) ? eps : 0)) * 3;
                        n++;
                    }
                    //Keep ZZ-Quad Fixed. Improves Optimization
//        			if(localMultipole[i][t002] != 0){
//        				localMultipole[i][t002] = x[n] + ((n == ivar) ? 0 : 0);
//        				n++;
//        			}
                    if (localMultipole[i][t110] != 0) {
                        localMultipole[i][t110] = (x[n] + ((n == ivar) ? eps : 0)) * 3;
                        n++;
                    }
                    if (localMultipole[i][t101] != 0) {
                        localMultipole[i][t101] = (x[n] + ((n == ivar) ? eps : 0)) * 3;
                        n++;
                    }
                    if (localMultipole[i][t011] != 0) {
                        localMultipole[i][t011] = (x[n] + ((n == ivar) ? eps : 0)) * 3;
                        n++;
                    }
                    localMultipole[i][t002] = -localMultipole[i][t200] - localMultipole[i][t020];
                }
                for (int k = 0; k < nAtoms; k++) {
                    if (k != i) {
                        Atom ak = atoms[k];
                        if (ak.getType() == itype) {
                            for (int j = 0; j < 10; j++) {
                                localMultipole[k][j] = localMultipole[i][j];
                            }
                        }
                    }
                }
            }
        }
    }

    /*Methods needed to implement Potential interface*/
    /**
     * <p>
     * energyAndGradient</p>
     *
     * @param x an array of double.
     * @param g an array of double.
     * @return a double.
     */
    public double energyAndGradient(double x[], double g[]) {
        //long currenttime = System.nanoTime();
        double pot_grid[][][] = new double[nSymm][][];
        Double xyz[] = new Double[3];
        //Change parameters
        varprm(x, -1, 0);
//    	for(int i = 0; i < x.length; i++){
//    		System.out.println(x[i]);
//    	}
        //rotate multipoles and induce
        init_prms();
        //calc new energy grid
        for (int i = 0; i < nSymm; i++) {
            pot_grid[i] = new double[target_grid[i].length][1];
            for (int j = 0; j < target_grid[i].length; j++) {
                xyz[0] = target_grid[i][j][0];
                xyz[1] = target_grid[i][j][1];
                xyz[2] = target_grid[i][j][2];
                pot_grid[i][j][0] = potpoint(xyz);
            }
        }
        //calc error
        double total_error = 0;
        double cscale = 100000000.0 / nSymm;
        double tscale = 10000.0 / nSymm;
        double eps = 0.000001;
        //double eps = 1;
        double netchg;
        int npoints = target_grid[0].length;
        double er = 0;
        double ec = 0;
        double et = 0;
        for (int i = 0; i < nSymm; i++) {
            netchg = 0;
            for (int j = 0; j < npoints; j++) {
                er += Math.pow(pot_grid[i][j][0] - target_grid[i][j][3], 2);
            }
            for (int j = 0; j < nAtoms; j++) {
                netchg += globalMultipole[i][j][0];
            }
            ec += cscale * Math.pow(netchg - Math.round(netchg), 2);
        }
        er = Math.sqrt(er / npoints);
        total_error = er + ec + et;
        //System.out.println(total_error+" "+er+" "+ec+" "+et+" "+npoints);
        //set up gradient array
        double e0 = 0;
        double e = 0;
        //int nvars = getNumberOfVariables();
        for (int k = 0; k < nvars; k++) {
            varprm(x, k, -.5 * eps);
            init_prms();
            //calc new energy grid
            for (int i = 0; i < nSymm; i++) {
                pot_grid[i] = new double[target_grid[i].length][1];
                for (int j = 0; j < target_grid[i].length; j++) {
                    xyz[0] = target_grid[i][j][0];
                    xyz[1] = target_grid[i][j][1];
                    xyz[2] = target_grid[i][j][2];
                    pot_grid[i][j][0] = potpoint(xyz);
                }
            }
            er = 0;
            ec = 0;
            et = 0;
            for (int i = 0; i < nSymm; i++) {
                netchg = 0;
                for (int j = 0; j < npoints; j++) {
                    er += Math.pow(pot_grid[i][j][0] - target_grid[i][j][3], 2);
                }
                //System.out.println("er2 = "+er);
                for (int j = 0; j < nAtoms; j++) {
                    netchg += globalMultipole[i][j][0];
                }
                ec += cscale * Math.pow(netchg - Math.round(netchg), 2);
                //System.out.println("ec2 = "+ec);
            }
            er = Math.sqrt(er / npoints);
            e0 = er + ec + et;
            varprm(x, k, .5 * eps);
            init_prms();
            for (int i = 0; i < nSymm; i++) {
                pot_grid[i] = new double[target_grid[i].length][1];
                for (int j = 0; j < target_grid[i].length; j++) {
                    xyz[0] = target_grid[i][j][0];
                    xyz[1] = target_grid[i][j][1];
                    xyz[2] = target_grid[i][j][2];
                    pot_grid[i][j][0] = potpoint(xyz);
                }
            }
            er = 0;
            ec = 0;
            et = 0;
            for (int i = 0; i < nSymm; i++) {
                netchg = 0;
                for (int j = 0; j < npoints; j++) {
                    er += Math.pow(pot_grid[i][j][0] - target_grid[i][j][3], 2);
                }
                for (int j = 0; j < nAtoms; j++) {
                    netchg += globalMultipole[i][j][0];
                }
                ec += cscale * Math.pow(netchg - Math.round(netchg), 2);
            }
            er = Math.sqrt(er / npoints);
            e = er + ec + et;
            g[k] = (e - e0) / eps;
//        	System.out.println(k+" g[k] = "+g[k]+" "+e0+" "+e);
        }
        //long endtime = System.nanoTime();
        //System.out.println("TIME "+(endtime - currenttime)*1e-9);
        totalEnergy = total_error;
        return total_error;
    }

    /**
     * <p>
     * Setter for the field <code>scaling</code>.</p>
     *
     * @param scaling an array of double.
     */
    public void setScaling(double scaling[]) {
        this.scaling = scaling;
    }

    /**
     * <p>
     * Getter for the field <code>scaling</code>.</p>
     *
     * @return an array of double.
     */
    public double[] getScaling() {
        return scaling;
    }

    /**
     * <p>
     * Getter for the field <code>coordinates</code>.</p>
     *
     * @param parameters an array of double.
     * @return an array of double.
     */
    public double[] getCoordinates(double parameters[]) {
        ArrayList<Double> p = new ArrayList<Double>();
        ArrayList<Integer> types = new ArrayList<Integer>();

        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            if (!types.contains(ai.getType())) {
                types.add(ai.getType());
                if (fitmpl) {
                    if (localMultipole[i][t000] != 0) {
                        p.add(localMultipole[i][t000]);
                    }
                }
                if (fitdpl) {
                    if (localMultipole[i][t100] != 0) {
                        p.add(localMultipole[i][t100]);
                    }
                    if (localMultipole[i][t010] != 0) {
                        p.add(localMultipole[i][t010]);
                    }
                    if (localMultipole[i][t001] != 0) {
                        p.add(localMultipole[i][t001]);
                    }
                }
                if (fitqdpl) {
                    if (localMultipole[i][t200] != 0) {
                        p.add(localMultipole[i][t200] / 3);
                    }
                    if (localMultipole[i][t020] != 0) {
                        p.add(localMultipole[i][t020] / 3);
                    }
                    //Keep ZZ-Quad Fixed. Improves Optimization
                    //if(localMultipole[i][t002] != 0) p.add(localMultipole[i][t002]*3);
                    //System.out.println(localMultipole[i][t002]);
                    if (localMultipole[i][t110] != 0) {
                        p.add(localMultipole[i][t110] / 3);
                    }
                    if (localMultipole[i][t101] != 0) {
                        p.add(localMultipole[i][t101] / 3);
                    }
                    if (localMultipole[i][t011] != 0) {
                        p.add(localMultipole[i][t011] / 3);
                    }
                }
            }
        }
        if (parameters == null) {
            parameters = new double[p.size()];
        }
        //System.out.println(p.size());
        for (int i = 0; i < p.size(); i++) {
            parameters[i] = p.get(i).doubleValue();
            //System.out.println("parameters[i] "+i+" "+parameters[i]);
        }
        return parameters;
    }

    /**
     * <p>
     * getallmpoles</p>
     *
     * @param x an array of double.
     * @return an array of double.
     */
    public double[] getallmpoles(double x[]) {
        ArrayList<Double> p = new ArrayList<Double>();
        ArrayList<Integer> types = new ArrayList<Integer>();
        int r = 0;
//    	for(int i = 0; i < nAtoms; i++){
//    		System.out.println(frame[i]);
//    		Atom ai = atoms[i];
//    		if(!types.contains(ai.getType())){
//    			types.add(ai.getType());
//    		}
//    	}
//    	Collections.sort(types);
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            if (!types.contains(ai.getType())) {
                types.add(ai.getType());
                if (fitmpl) {
                    if (localMultipole[i][t000] == 0) {
                        p.add(localMultipole[i][t000]);
                    } else {
                        p.add(x[r]);
                        r = r + 1;
                    }
                }
                if (fitdpl) {
                    if (localMultipole[i][t100] == 0) {
                        p.add(localMultipole[i][t100]);
                    } else {
                        p.add(x[r]);
                        r = r + 1;
                    }
                    if (localMultipole[i][t010] == 0) {
                        p.add(localMultipole[i][t010]);
                    } else {
                        p.add(x[r]);
                        r = r + 1;
                    }
                    if (localMultipole[i][t001] == 0) {
                        p.add(localMultipole[i][t001]);
                    } else {
                        p.add(x[r]);
                        r = r + 1;
                    }
                }
                if (fitqdpl) {
                    if (localMultipole[i][t200] == 0) {
                        p.add(localMultipole[i][t200]);
                    } else {
                        p.add(x[r] * 3);
                        r = r + 1;
                    }
                    if (localMultipole[i][t020] == 0) {
                        p.add(localMultipole[i][t020]);
                    } else {
                        p.add(x[r] * 3);
                        r = r + 1;
                    }
                    //Keep ZZ-Quad Fixed. Improves Optimization
                    p.add(localMultipole[i][t002]);
                    //System.out.println(localMultipole[i][t002]);
                    //if(localMultipole[i][t002] == 0){
                    //	p.add(localMultipole[i][t002]);
                    //}
                    //else{
                    //	p.add(x[r]);
                    //	r = r + 1;
                    //}
                    if (localMultipole[i][t110] == 0) {
                        p.add(localMultipole[i][t110]);
                    } else {
                        p.add(x[r] * 3);
                        r = r + 1;
                    }
                    if (localMultipole[i][t101] == 0) {
                        p.add(localMultipole[i][t101]);
                    } else {
                        p.add(x[r] * 3);
                        r = r + 1;
                    }
                    if (localMultipole[i][t011] == 0) {
                        p.add(localMultipole[i][t011]);
                    } else {
                        p.add(x[r] * 3);
                        r = r + 1;
                    }
                }
            }
        }
        double parameters[] = new double[p.size()];
        for (int i = 0; i < p.size(); i++) {
            parameters[i] = p.get(i).doubleValue();
        }
        return parameters;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getMass() {
        return null;
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        return nvars;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    /**
     * <p>
     * Setter for the field <code>key</code>.</p>
     *
     * @param key a {@link java.util.ArrayList} object.
     */
    public void setkey(ArrayList<String> key) {
        this.key = key;
    }

    //added by gchattree
    /**
     * <p>
     * polargrp</p>
     */
    public static void polargrp() {
        int index;

        ArrayList<ArrayList<Integer>> temp_groups1 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> temp_groups2 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> temp_groups3 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> temp_groups4 = new ArrayList<ArrayList<Integer>>();

        ArrayList<Integer> polarizationGroup = new ArrayList<Integer>();

        ArrayList<Integer> list = new ArrayList<Integer>();
        int nlist = 0;
        ArrayList<Integer> keep = new ArrayList<Integer>();
        //int nkeep = 0;
        ArrayList<Integer> mask = new ArrayList<Integer>();

        ArrayList<Integer> jg;
        ArrayList<Integer> ig;
        int kk;
        int jj;
        int start;
        int stop;
        boolean done;

        for (Atom ai : atoms) {
            ArrayList<Integer> group = new ArrayList<Integer>();
            polarizationGroup.clear();
            index = ai.getXYZIndex() - 1;
            group.add(index);
            //polarizationGroup.add(ai.getType());
            PolarizeType polarizeType = ai.getPolarizeType();
            if (polarizeType != null) {
                if (polarizeType.polarizationGroup != null) {
                    for (int i : polarizeType.polarizationGroup) {
                        if (!polarizationGroup.contains(i)) {
                            polarizationGroup.add(i);
                        }
                    }
                }
            }
            for (Bond bi : ai.getFFXBonds()) {
                Atom aj = bi.get1_2(ai);
                int tj = aj.getType();
                for (int g : polarizationGroup) {
                    if (g == tj) {
                        Integer index2 = aj.getXYZIndex() - 1;
                        group.add(index2);
                    }
                }
            }
            Collections.sort(group);
            temp_groups1.add(group);
        }

        //Next part of ip11 creation
        for (int n = 0; n < nAtoms; n++) {
            list.add(n, -1);
        }

        for (int i = 0; i < nAtoms; i++) {
            ig = temp_groups1.get(i);
            done = false;
            start = 1;
            stop = ig.size();
            for (int j = start - 1; j < stop; j++) {
                jj = ig.get(j);
                if (jj < i) {
                    done = true;
                    jg = temp_groups1.get(jj);
                    for (int k = 0; k < jg.size(); k++) {
                        if (k > ig.size() - 1) {
                            for (int s = ig.size(); s < k + 1; s++) {
                                ig.add(0);
                            }
                            ig.set(k, jg.get(k));
                        } else {
                            ig.set(k, jg.get(k));
                        }
                    }
                } else {
                    list.set(jj, i);
                }
            }
            while (!done) {
                done = true;
                for (int j = start - 1; j < stop; j++) {
                    jj = ig.get(j);
                    jg = temp_groups1.get(jj);
                    for (int k = 0; k < jg.size(); k++) {
                        kk = jg.get(k);
                        if (list.get(kk) != i) {
                            ig.add(kk);
                            list.set(kk, i);
                        }
                    }
                }
                if (ig.size() != stop) {
                    done = false;
                    start = stop + 1;
                    stop = ig.size();
                }
            }
            Collections.sort(ig);
        }

        //final part of ip11 array creation
        for (int n = 0; n < nAtoms; n++) {
            ArrayList<Integer> group = temp_groups1.get(n);
            Collections.sort(group);
            //System.out.println(group);
            ip11[n] = new int[group.size()];
            int j = 0;
            for (int k : group) {
                ip11[n][j++] = k;
            }
        }

        //start ip12 creation
        for (int n = 0; n < nAtoms; n++) {
            mask.add(n, -1);
        }
        for (int i = 0; i < nAtoms; i++) {
            list.clear();
            keep.clear();
            ArrayList<Integer> group = new ArrayList<Integer>();
            ig = temp_groups1.get(i);
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                list.add(jj);
                mask.set(jj, i);
            }
            for (int j = 0; j < list.size(); j++) {
                jj = list.get(j);
                Atom ajj = atoms[jj];
                for (int k = 0; k < ajj.getFFXBonds().size(); k++) {
                    kk = ajj.getFFXBonds().get(k).get1_2(ajj).getXYZIndex() - 1;
                    //System.out.println(mask.get(kk)+" "+i);
                    if (mask.get(kk) != i) {
                        keep.add(kk);
                    }
                }
            }
            nlist = 0;
            list.clear();
            for (int j = 0; j < keep.size(); j++) {
                jj = keep.get(j);
                jg = temp_groups1.get(jj);
                for (int k = 0; k < jg.size(); k++) {
                    kk = jg.get(k);
                    //System.out.println((j+1)+" "+(jj+1)+" "+(k+1)+" "+(kk+1));
                    nlist++;
                    //list.set(nlist, kk);
                    if (nlist - 1 < list.size()) {
                        list.set(nlist - 1, kk);
                    } else {
                        list.add(kk);
                    }
                }
            }
            Collections.sort(list);
            for (int j = 0; j < list.size(); j++) {
                group.add(j, list.get(j));
            }
            temp_groups2.add(group);
        }
        //final part of ip12 array creation
        for (int n = 0; n < nAtoms; n++) {
            ArrayList<Integer> group = temp_groups2.get(n);
            Collections.sort(group);
            //System.out.println(group);
            ip12[n] = new int[group.size()];
            int j = 0;
            for (int k : group) {
                ip12[n][j++] = k;
            }
        }

        //start ip13 creation
        mask.clear();
        for (int n = 0; n < nAtoms; n++) {
            mask.add(n, -1);
        }
        for (int i = 0; i < nAtoms; i++) {
            list.clear();
            ArrayList<Integer> group = new ArrayList<Integer>();
            ig = temp_groups1.get(i);
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                mask.set(jj, i);
            }
            ig = temp_groups2.get(i);
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                mask.set(jj, i);
            }
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                jg = temp_groups2.get(jj);
                for (int k = 0; k < jg.size(); k++) {
                    kk = jg.get(k);
                    if (mask.get(kk) != i) {
                        list.add(kk);
                    }
                }
            }
            Collections.sort(list);
            for (int j = 0; j < list.size(); j++) {
                group.add(j, list.get(j));
            }
            temp_groups3.add(group);
        }
        //final part of ip13 array creation
        for (int n = 0; n < nAtoms; n++) {
            ArrayList<Integer> group = temp_groups3.get(n);
            Collections.sort(group);
            //System.out.println(group);
            ip13[n] = new int[group.size()];
            int j = 0;
            for (int k : group) {
                ip13[n][j++] = k;
            }
        }

        //start ip14 creation
        mask.clear();
        for (int n = 0; n < nAtoms; n++) {
            mask.add(n, -1);
        }
        for (int i = 0; i < nAtoms; i++) {
            list.clear();
            ArrayList<Integer> group = new ArrayList<Integer>();
            ig = temp_groups1.get(i);
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                mask.set(jj, i);
            }
            ig = temp_groups2.get(i);
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                mask.set(jj, i);
            }
            ig = temp_groups3.get(i);
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                mask.set(jj, i);
            }
            for (int j = 0; j < ig.size(); j++) {
                jj = ig.get(j);
                jg = temp_groups2.get(jj);
                for (int k = 0; k < jg.size(); k++) {
                    kk = jg.get(k);
                    if (mask.get(kk) != i) {
                        list.add(kk);
                    }
                }
            }
            Collections.sort(list);
            for (int j = 0; j < list.size(); j++) {
                group.add(j, list.get(j));
            }
            temp_groups4.add(group);
        }
        //final part of ip14 array creation
        for (int n = 0; n < nAtoms; n++) {
            ArrayList<Integer> group = temp_groups3.get(n);
            Collections.sort(group);
            //System.out.println(group);
            ip14[n] = new int[group.size()];
            int j = 0;
            for (int k : group) {
                ip14[n][j++] = k;
            }
        }
    }

    /**
     * <p>
     * induce_pedit</p>
     */
    public static void induce_pedit() {
        double d12scale = 1;
        double d13scale = 1;
        double d14scale = 1;
        double mask_local[];
        double maskp_local[];
        double fx_local[];
        double fy_local[];
        double fz_local[];
        double fxp_local[];
        double fyp_local[];
        double fzp_local[];
        double dx_local[];
        mask_local = new double[nAtoms];
        maskp_local = new double[nAtoms];
        fx_local = new double[nAtoms];
        fy_local = new double[nAtoms];
        fz_local = new double[nAtoms];
        fxp_local = new double[nAtoms];
        fyp_local = new double[nAtoms];
        fzp_local = new double[nAtoms];
        dx_local = new double[3];

        for (int i = 0; i < nAtoms; i++) {
            mask_local[i] = 1.0;
            maskp_local[i] = 1.0;
        }

        for (int i = 0; i < nAtoms; i++) {
            fx_local[i] = 0.0;
            fy_local[i] = 0.0;
            fz_local[i] = 0.0;
            fxp_local[i] = 0.0;
            fyp_local[i] = 0.0;
            fzp_local[i] = 0.0;
        }
        int lists[][] = neighborLists[0];
        final double x[] = coordinates[0][0];
        final double y[] = coordinates[0][1];
        final double z[] = coordinates[0][2];
        final double mpole[][] = globalMultipole[0];
        /**
         * System.out.println((i+1)+" "+ind[0]+" "+polar+" "+fieldi[0]); Loop
         * over atoms
         *
         */
        for (int i = 0; i < (nAtoms - 1); i++) {
            final double pdi = pdamp[i];
            final double pti = thole[i];
            final double xi = x[i];
            final double yi = y[i];
            final double zi = z[i];
            final double globalMultipolei[] = mpole[i];
            final double ci = globalMultipolei[0];
            final double dix = globalMultipolei[1];
            final double diy = globalMultipolei[2];
            final double diz = globalMultipolei[3];
            final double qixx = globalMultipolei[4];
            final double qiyy = globalMultipolei[8];
            final double qizz = globalMultipolei[12];
            final double qixy = globalMultipolei[5];
            final double qixz = globalMultipolei[6];
            final double qiyz = globalMultipolei[9];
            /**
             * Apply energy masking rules. (not working for maskp)
             */
            Atom ai = atoms[i];
            for (Torsion torsion : ai.getTorsions()) {
                Atom ak = torsion.get1_4(ai);
                if (ak != null) {
                    int index = ak.xyzIndex - 1;
                    for (int k : ip11[i]) {
                        if (k == index) {
                            maskp_local[index] = 0.5;
                        }
                    }
                }
            }
            for (Angle angle : ai.getAngles()) {
                Atom ak = angle.get1_3(ai);
                if (ak != null) {
                    int index = ak.xyzIndex - 1;
                    maskp_local[index] = p13scale;
                }
            }
            for (Bond bond : ai.getFFXBonds()) {
                int index = bond.get1_2(ai).xyzIndex - 1;
                maskp_local[index] = p12scale;
            }
            /**
             * Apply group based polarization masking rule.1.334^(1/6)
             *
             */
            for (int index : ip11[i]) {
                mask_local[index] = d11scale;

            }
            for (int index : ip12[i]) {
                mask_local[index] = d12scale;

            }
            for (int index : ip13[i]) {
                mask_local[index] = d13scale;
            }
            for (int index : ip14[i]) {
                mask_local[index] = d14scale;
            }
            for (int k = i; k < nAtoms; k++) {
                if (i != k) {
                    final double xk = x[k];
                    final double yk = y[k];
                    final double zk = z[k];
                    dx_local[0] = xk - xi;
                    dx_local[1] = yk - yi;
                    dx_local[2] = zk - zi;
                    //final double r2 = crystal.image(dx_local);
                    final double r2 = dx_local[0] * dx_local[0] + dx_local[1] * dx_local[1] + dx_local[2] * dx_local[2];
                    if (r2 <= off2) {
                        final double xr = dx_local[0];
                        final double yr = dx_local[1];
                        final double zr = dx_local[2];
                        final double pdk = pdamp[k];
                        final double ptk = thole[k];
                        final double globalMultipolek[] = mpole[k];
                        final double ck = globalMultipolek[0];
                        final double dkx = globalMultipolek[1];
                        final double dky = globalMultipolek[2];
                        final double dkz = globalMultipolek[3];
                        final double qkxx = globalMultipolek[4];
                        final double qkyy = globalMultipolek[8];
                        final double qkzz = globalMultipolek[12];
                        final double qkxy = globalMultipolek[5];
                        final double qkxz = globalMultipolek[6];
                        final double qkyz = globalMultipolek[9];
                        final double r = sqrt(r2);
                        /**
                         * Compute the error function scaled and unscaled terms.
                         */
                        double scale3 = mask_local[k];
                        double scale5 = mask_local[k];
                        double scale7 = mask_local[k];
                        double damp = pdi * pdk;
                        double expdamp = 0.0;
                        if (damp != 0.0) {
                            final double pgamma = min(pti, ptk);
                            final double rdamp = r / damp;
                            damp = -pgamma * rdamp * rdamp * rdamp;
                            if (damp > -50.0) {
                                expdamp = exp(damp);
                                scale3 = scale3 * (1.0 - expdamp);
                                scale5 = scale5 * (1.0 - expdamp * (1.0 - damp));
                                scale7 = scale7 * (1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp));
                            }
                        }
                        final double rr3 = scale3 * 1.0 / (r * r2);
                        final double rr5 = (scale3 != 0) ? scale5 / scale3 * 3.0 * rr3 / r2 : 0;
                        final double rr7 = (scale5 != 0) ? scale7 / scale5 * 5.0 * rr5 / r2 : 0;
                        final double dir = dix * xr + diy * yr + diz * zr;
                        final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                        final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                        final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                        final double qir = (qix * xr + qiy * yr + qiz * zr) / 2.0;
                        final double dkr = dkx * xr + dky * yr + dkz * zr;
                        final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                        final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                        final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                        final double qkr = (qkx * xr + qky * yr + qkz * zr) / 2.0;
                        final double fidx = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkx + rr5 * qkx;
                        final double fidy = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dky + rr5 * qky;
                        final double fidz = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkz + rr5 * qkz;
                        final double fkdx = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * dix - rr5 * qix;
                        final double fkdy = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diy - rr5 * qiy;
                        final double fkdz = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diz - rr5 * qiz;
                        //System.out.println((i+1)+" "+(k+1)+" "+xr+" "+rr3+" "+ck+" "+rr5+" "+dkr+" "+rr7+" "+qkr);
                        fx_local[i] += fidx;
                        fy_local[i] += fidy;
                        fz_local[i] += fidz;
                        fx_local[k] += fkdx;
                        fy_local[k] += fkdy;
                        fz_local[k] += fkdz;
                    }
                }
            }
        }
        //set up fields
        for (int i = 0; i < nAtoms; i++) {
            double fieldi[] = field1[i];
            double field2i[] = field2[i];
            fieldi[0] = fx_local[i];
            fieldi[1] = fy_local[i];
            fieldi[2] = fz_local[i];
            //System.out.println(fx_local[i]+" "+fy_local[i]+" "+fz_local[i]);
            field2i[0] = fxp_local[i];
            field2i[1] = fyp_local[i];
            field2i[2] = fzp_local[i];
        }
        setDipoleMoments(false);
    }

    //added by gchattree
    /**
     * <p>
     * induce0a</p>
     */
    public static void induce0a() {
        double d12scale = 1;
        double d13scale = 1;
        double d14scale = 1;
        double mask_local[];
        double maskp_local[];
        double fx_local[];
        double fy_local[];
        double fz_local[];
        double fxp_local[];
        double fyp_local[];
        double fzp_local[];
        double dx_local[];
        mask_local = new double[nAtoms];
        maskp_local = new double[nAtoms];
        fx_local = new double[nAtoms];
        fy_local = new double[nAtoms];
        fz_local = new double[nAtoms];
        fxp_local = new double[nAtoms];
        fyp_local = new double[nAtoms];
        fzp_local = new double[nAtoms];
        dx_local = new double[3];

        for (int i = 0; i < nAtoms; i++) {
            mask_local[i] = 1.0;
            maskp_local[i] = 1.0;
        }

        for (int i = 0; i < nAtoms; i++) {
            fx_local[i] = 0.0;
            fy_local[i] = 0.0;
            fz_local[i] = 0.0;
            fxp_local[i] = 0.0;
            fyp_local[i] = 0.0;
            fzp_local[i] = 0.0;
        }
        int lists[][] = neighborLists[0];
        final double x[] = coordinates[0][0];
        final double y[] = coordinates[0][1];
        final double z[] = coordinates[0][2];
        final double mpole[][] = globalMultipole[0];
        /**
         * Loop over atoms
         *
         */
        for (int i = 0; i < (nAtoms - 1); i++) {
            final double pdi = pdamp[i];
            final double pti = thole[i];
            final double xi = x[i];
            final double yi = y[i];
            final double zi = z[i];
            final double globalMultipolei[] = mpole[i];
            final double ci = globalMultipolei[0];
            final double dix = globalMultipolei[t100];
            final double diy = globalMultipolei[t010];
            final double diz = globalMultipolei[t001];
            final double qixx = globalMultipolei[t200] / 3.0;
            final double qiyy = globalMultipolei[t020] / 3.0;
            final double qizz = globalMultipolei[t002] / 3.0;
            final double qixy = globalMultipolei[t110] / 3.0;
            final double qixz = globalMultipolei[t101] / 3.0;
            final double qiyz = globalMultipolei[t011] / 3.0;
            /**
             * Apply energy masking rules. (not working for maskp)
             */
            Atom ai = atoms[i];
            for (Torsion torsion : ai.getTorsions()) {
                Atom ak = torsion.get1_4(ai);
                if (ak != null) {
                    int index = ak.xyzIndex - 1;
                    for (int k : ip11[i]) {
                        if (k == index) {
                            maskp_local[index] = 0.5;
                        }
                    }
                }
            }
            for (Angle angle : ai.getAngles()) {
                Atom ak = angle.get1_3(ai);
                if (ak != null) {
                    int index = ak.xyzIndex - 1;
                    maskp_local[index] = p13scale;
                }
            }
            for (Bond bond : ai.getFFXBonds()) {
                int index = bond.get1_2(ai).xyzIndex - 1;
                maskp_local[index] = p12scale;
            }
            /**
             * Apply group based polarization masking rule.
             *
             */
            for (int index : ip11[i]) {
                mask_local[index] = d11scale;

            }
            for (int index : ip12[i]) {
                mask_local[index] = d12scale;

            }
            for (int index : ip13[i]) {
                mask_local[index] = d13scale;
            }
            for (int index : ip14[i]) {
                mask_local[index] = d14scale;
            }
//            if(i == 3){
//            	for (int index : ip11[i]) System.out.println("1 "+index);
//            	for (int index : ip12[i]) System.out.println("2 "+index);
//            	for (int index : ip13[i]) System.out.println("3 "+index);
//            	for (int index : ip14[i]) System.out.println("4 "+index);
//            }
            for (int k = i; k < nAtoms; k++) {
                if (i != k) {
                    final double xk = x[k];
                    final double yk = y[k];
                    final double zk = z[k];
                    dx_local[0] = xk - xi;
                    dx_local[1] = yk - yi;
                    dx_local[2] = zk - zi;
                    //final double r2 = crystal.image(dx_local);
                    final double r2 = dx_local[0] * dx_local[0] + dx_local[1] * dx_local[1] + dx_local[2] * dx_local[2];
                    if (r2 <= off2) {
                        final double xr = dx_local[0];
                        final double yr = dx_local[1];
                        final double zr = dx_local[2];
                        final double pdk = pdamp[k];
                        final double ptk = thole[k];
                        final double globalMultipolek[] = mpole[k];
                        final double ck = globalMultipolek[t000];
                        final double dkx = globalMultipolek[t100];
                        final double dky = globalMultipolek[t010];
                        final double dkz = globalMultipolek[t001];
                        final double qkxx = globalMultipolek[t200] / 3.0;
                        final double qkyy = globalMultipolek[t020] / 3.0;
                        final double qkzz = globalMultipolek[t002] / 3.0;
                        final double qkxy = globalMultipolek[t110] / 3.0;
                        final double qkxz = globalMultipolek[t101] / 3.0;
                        final double qkyz = globalMultipolek[t011] / 3.0;
                        final double r = sqrt(r2);
                        /**
                         * Compute the error function scaled and unscaled terms.
                         */
                        double scale3 = 1.0;
                        double scale5 = 1.0;
                        double scale7 = 1.0;
                        double damp = pdi * pdk;
                        double expdamp = 0.0;
                        if (damp != 0.0) {
                            final double pgamma = min(pti, ptk);
                            final double rdamp = r / damp;
                            damp = -pgamma * rdamp * rdamp * rdamp;
                            if (damp > -50.0) {
                                expdamp = exp(damp);
                                scale3 = 1.0 - expdamp;
                                scale5 = 1.0 - expdamp * (1.0 - damp);
                                scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                            }
                        }
                        final double rr3 = scale3 * 1.0 / (r * r2);
                        final double rr5 = (scale3 != 0) ? scale5 / scale3 * 3.0 * rr3 / r2 : 0;
                        final double rr7 = (scale5 != 0) ? scale7 / scale5 * 5.0 * rr5 / r2 : 0;
                        final double dir = dix * xr + diy * yr + diz * zr;
                        final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                        final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                        final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                        final double qir = (qix * xr + qiy * yr + qiz * zr) / 2.0;
                        final double dkr = dkx * xr + dky * yr + dkz * zr;
                        final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                        final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                        final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                        final double qkr = (qkx * xr + qky * yr + qkz * zr) / 2.0;
                        final double fidx = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkx + rr5 * qkx;
                        final double fidy = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dky + rr5 * qky;
                        final double fidz = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkz + rr5 * qkz;
                        final double fkdx = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * dix - rr5 * qix;
                        final double fkdy = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diy - rr5 * qiy;
                        final double fkdz = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diz - rr5 * qiz;
                        final double fipx = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkx + rr5 * qkx;
                        final double fipy = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dky + rr5 * qky;
                        final double fipz = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkz + rr5 * qkz;
                        final double fkpx = xr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * dix - rr5 * qix;
                        final double fkpy = yr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diy - rr5 * qiy;
                        final double fkpz = zr * (rr3 * ci + rr5 * dir + rr7 * qir) - rr3 * diz - rr5 * qiz;

                        fx_local[i] += fidx * mask_local[k];
                        fy_local[i] += fidy * mask_local[k];
                        fz_local[i] += fidz * mask_local[k];
                        fx_local[k] += fkdx * mask_local[k];
                        fy_local[k] += fkdy * mask_local[k];
                        fz_local[k] += fkdz * mask_local[k];
                        //System.out.println(mask_local[k]);
                        //not working
                        fxp_local[i] += fipx * maskp_local[k];
                        //System.out.println((i+1)+" "+(k+1)+" "+xr+" "+rr3+" "+ck+" "+rr5+" "+dkr+" "+rr7+" "+qkr);
                        fyp_local[i] += fipy * maskp_local[k];
                        fzp_local[i] += fipz * maskp_local[k];
                        fxp_local[k] += fkpx * maskp_local[k];
                        fyp_local[k] += fkpy * maskp_local[k];
                        fzp_local[k] += fkpz * maskp_local[k];
                    }
                }
            }
        }
        /**
         * Loop over symmetry mates. THIS IS PROBABLY NOT RIGHT, CHECK!
         */
//        for (int iSymm = 1; iSymm < nSymm; iSymm++) {
//            lists = neighborLists[iSymm];
//            double xs[] = coordinates[iSymm][0];
//            double ys[] = coordinates[iSymm][1];
//            double zs[] = coordinates[iSymm][2];
//            double mpoles[][] = globalMultipole[iSymm];
//            /**
//             * Loop over atoms in a chunk of the asymmetric unit.
//             */
//            for (int i = 0; i < nAtoms; i++) {
//                /**
//                 * Apply energy masking rules. (not working for maskp)
//                 */
//                Atom ai = atoms[i];
//                for (Torsion torsion : ai.getTorsions()) {
//                    Atom ak = torsion.get1_4(ai);
//                    if (ak != null) {
//                        int index = ak.xyzIndex - 1;
//                        for (int k : ip11[i]) {
//                            if (k == index) {
//                                maskp_local[index] = 0.5;
//                            }
//                        }
//                    }
//                }
//                for (Angle angle : ai.getAngles()) {
//                    Atom ak = angle.get1_3(ai);
//                    if (ak != null) {
//                        int index = ak.xyzIndex - 1;
//                        maskp_local[index] = p13scale;
//                    }
//                }
//                for (Bond bond : ai.getFFXBonds()) {
//                    int index = bond.get1_2(ai).xyzIndex - 1;
//                    maskp_local[index] = p12scale;
//                }
//                /**
//                 * Apply group based polarization masking rule.
//                 * FIX
//                 */
//                for (int index : ip11[i]) {
//                    mask_local[index] = d11scale;
//                }
//                for (int index : ip12[i]) {
//                    mask_local[index] = d12scale;
//                }
//                for (int index : ip13[i]) {
//                    mask_local[index] = d13scale;
//                }
//                for (int index : ip14[i]) {
//                    mask_local[index] = d14scale;
//                }
//
//                final double pdi = pdamp[i];
//                final double pti = thole[i];
//                final double xi = x[i];
//                final double yi = y[i];
//                final double zi = z[i];
//                /**
//                 * Loop over the neighbor list.
//                 */
//                final int list[] = lists[i];
//                final int npair = list.length;
//                for (int j = 0; j < npair; j++) {
//                    int k = list[j];
//                    final double xk = xs[k];
//                    final double yk = ys[k];
//                    final double zk = zs[k];
//                    dx_local[0] = xk - xi;
//                    dx_local[1] = yk - yi;
//                    dx_local[2] = zk - zi;
//                    final double r2 = crystal.image(dx_local);
//                    if (r2 <= off2) {
//                        //ewald[counts[i]++] = k;
//                        final double xr = dx_local[0];
//                        final double yr = dx_local[1];
//                        final double zr = dx_local[2];
//                        final double pdk = pdamp[k];
//                        final double ptk = thole[k];
//                        final double multipolek[] = mpoles[k];
//                        final double ck = multipolek[t000];
//                        final double dkx = multipolek[t100];
//                        final double dky = multipolek[t010];
//                        final double dkz = multipolek[t001];
//                        final double qkxx = multipolek[t200] / 3.0;
//                        final double qkyy = multipolek[t020] / 3.0;
//                        final double qkzz = multipolek[t002] / 3.0;
//                        final double qkxy = multipolek[t110] / 3.0;
//                        final double qkxz = multipolek[t101] / 3.0;
//                        final double qkyz = multipolek[t011] / 3.0;
//                        final double r = sqrt(r2);
//                        /**
//                         * Compute the error function scaled and
//                         * unscaled terms.
//                         */
//                        double scale3 = 1.0;
//                        double scale5 = 1.0;
//                        double scale7 = 1.0;
//                        double damp = pdi * pdk;
//                        double expdamp = 0.0;
//                        if (damp != 0.0) {
//                            final double pgamma = min(pti, ptk);
//                            final double rdamp = r / damp;
//                            damp = -pgamma * rdamp * rdamp * rdamp;
//                            if (damp > -50.0) {
//                                expdamp = exp(damp);
//                                scale3 = 1.0 - expdamp;
//                                scale5 = 1.0 - expdamp * (1.0 - damp);
//                                scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
//                            }
//                        }
//                        final double rr3 = scale3 * 1.0 / (r * r2);
//                        final double rr5 = scale5 * 3.0 * rr3 / r2;
//                        final double rr7 = scale7 * 5.0 * rr5 / r2;
//                        final double dkr = dkx * xr + dky * yr + dkz * zr;
//                        final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
//                        final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
//                        final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
//                        final double qkr = (qkx * xr + qky * yr + qkz * zr) / 2.0;
//                        final double fidx = -xr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkx + rr5 * qkx;
//                        final double fidy = -yr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dky + rr5 * qky;
//                        final double fidz = -zr * (rr3 * ck - rr5 * dkr + rr7 * qkr) - rr3 * dkz + rr5 * qkz;
//                        fx_local[i] += fidx;
//                        fy_local[i] += fidy;
//                        fz_local[i] += fidz;
//                        fxp_local[i] += fidx;
//                        fyp_local[i] += fidy;
//                        fzp_local[i] += fidz;
//                    }
//                }
//            }
//        }
        //set up fields
        for (int i = 0; i < nAtoms; i++) {
            double fieldi[] = field1[i];
            double field2i[] = field2[i];
            fieldi[0] = fx_local[i];
            fieldi[1] = fy_local[i];
            fieldi[2] = fz_local[i];
            field2i[0] = fxp_local[i];
            field2i[1] = fyp_local[i];
            field2i[2] = fzp_local[i];
        }
        setDipoleMoments(false);
    }

    /**
     * <p>
     * setDipoleMoments</p>
     *
     * @param print a boolean.
     */
    public static void setDipoleMoments(boolean print) {
        long startTime = System.nanoTime();
        //set the induced dipoles
        if (polarization == Polarization.NONE) {
            for (int i = 0; i < nAtoms; i++) {
                inducedDipole[0][i][0] = 0.0;
                inducedDipole[0][i][1] = 0.0;
                inducedDipole[0][i][2] = 0.0;
                inducedDipolep[0][i][0] = 0.0;
                inducedDipolep[0][i][1] = 0.0;
                inducedDipolep[0][i][2] = 0.0;
            }
            try {
                //CHECK! Do I need this?
                //parallelTeam.execute(expandInducedDipolesRegion);
            } catch (Exception e) {
                e.printStackTrace();
            }
            return;
        }
        /**
         * Set the induced dipoles to the polarizability times the direct field.
         */
        final double induced0[][] = inducedDipole[0];
        final double inducedp0[][] = inducedDipolep[0];
        for (int i = 0; i < nAtoms; i++) {
            final double polar = polarizability[i];
            final double fieldi[] = field1[i];
            final double ind[] = induced0[i];
            final double directi[] = directDipole[i];
            ind[0] = polar * fieldi[0];
            ind[1] = polar * fieldi[1];
            ind[2] = polar * fieldi[2];
            //System.out.println((i+1)+" "+ind[0]+" "+polar+" "+fieldi[0]);
            directi[0] = ind[0];
            directi[1] = ind[1];
            directi[2] = ind[2];
            final double field2i[] = field2[i];
            final double inp[] = inducedp0[i];
            final double directpi[] = directDipolep[i];
            inp[0] = polar * field2i[0];
            inp[1] = polar * field2i[1];
            inp[2] = polar * field2i[2];
            directpi[0] = inp[0];
            directpi[1] = inp[1];
            directpi[2] = inp[2];
        }
        //Expands dipole moments to other symmetries?
        //CHECK. Need this expandInducedDipolesRegion?
//        try {
//            parallelTeam.execute(expandInducedDipolesRegion);
//        } catch (Exception e) {
//            e.printStackTrace();
//        }

        if (polarization == Polarization.MUTUAL) {
            calcMutualDipoleMoments(print, startTime);
        }
    }

    /**
     * <p>
     * calcMutualDipoleMoments</p>
     *
     * @param print a boolean.
     * @param startTime a long.
     */
    public static void calcMutualDipoleMoments(boolean print, long startTime) {
        double mask_local[] = new double[nAtoms];
        double dx_local[] = new double[3];
        final double induced0[][] = inducedDipole[0];
        final double inducedp0[][] = inducedDipolep[0];
        double d12scale = 1;
        double d13scale = 1;
        double d14scale = 1;
        StringBuffer sb = null;
        long directTime = System.nanoTime() - startTime;
        if (print) {
            sb = new StringBuffer(
                    "\n SELF-CONSISTENT FIELD\n"
                    + " Iter     RMS Change (Debyes)   Time\n");
        }

        boolean done = false;
        int maxiter = 1000;
        int iter = 0;
        double eps = 100.0;
        double epsold;

        while (!done) {
            //reset fields
            for (int i = 0; i < nAtoms; i++) {
                double fieldi[] = field1[i];
                double field2i[] = field2[i];
                fieldi[0] = 0;
                fieldi[1] = 0;
                fieldi[2] = 0;
                field2i[0] = 0;
                field2i[1] = 0;
                field2i[2] = 0;
            }

            double fx_local[] = new double[nAtoms];
            double fy_local[] = new double[nAtoms];
            double fz_local[] = new double[nAtoms];
            double fxp_local[] = new double[nAtoms];
            double fyp_local[] = new double[nAtoms];
            double fzp_local[] = new double[nAtoms];

            long cycleTime = -System.nanoTime();
            for (int i = 0; i < (nAtoms - 1); i++) {
                double ind[] = induced0[i];
                double indp[] = inducedp0[i];
                double xyzi[] = atoms[i].getXYZ(null);
                double pdi = pdamp[i];
                double pti = thole[i];
                double duix = ind[0];
                double duiy = ind[1];
                double duiz = ind[2];
                double puix = indp[0];
                double puiy = indp[1];
                double puiz = indp[2];
                for (int j = i + 1; j < nAtoms; j++) {
                    //If it is not pedit, mask_local is set to 1 or the 'dXscale'
                    mask_local[j] = !pedit ? 1 : 0;
                }
                for (int index : ip11[i]) {
                    double u1scale = 1;
                    if (pedit) {
                        u1scale = u1scale - d11scale;
                    }
                    mask_local[index] = u1scale;
                }
                for (int index : ip12[i]) {
                    double u2scale = 1;
                    if (pedit) {
                        u2scale = u2scale - d12scale;
                    }
                    mask_local[index] = u2scale;
                }
                for (int index : ip13[i]) {
                    double u3scale = 1;
                    if (pedit) {
                        u3scale = u3scale - d13scale;
                    }
                    mask_local[index] = u3scale;
                }
                for (int index : ip14[i]) {
                    double u4scale = 1;
                    if (pedit) {
                        u4scale = u4scale - d14scale;
                    }
                    mask_local[index] = u4scale;
                }

                for (int k = i + 1; k < nAtoms; k++) {
                    double indk[] = induced0[k];
                    double indkp[] = inducedp0[k];
                    double xyzk[] = atoms[k].getXYZ(null);
                    dx_local[0] = xyzk[0] - xyzi[0];
                    dx_local[1] = xyzk[1] - xyzi[1];
                    dx_local[2] = xyzk[2] - xyzi[2];
                    //final double r2 = crystal.image(dx_local);
                    final double r2 = dx_local[0] * dx_local[0] + dx_local[1] * dx_local[1] + dx_local[2] * dx_local[2];
                    if (r2 <= off2) {
                        double xr = dx_local[0];
                        double yr = dx_local[1];
                        double zr = dx_local[2];
                        double r = sqrt(r2);
                        double dukx = indk[0];
                        double duky = indk[1];
                        double dukz = indk[2];
                        double pukx = indkp[0];
                        double puky = indkp[1];
                        double pukz = indkp[2];
                        double scale3 = mask_local[k];
                        double scale5 = mask_local[k];
                        double damp = pdi * pdamp[k];
                        if (damp != 0) {
                            final double pgamma = min(pti, thole[k]);
                            damp = -pgamma * pow((r / damp), 3);
                            if (damp > -50) {
                                double expdamp = exp(damp);
                                scale3 = scale3 * (1 - expdamp);
                                scale5 = scale5 * (1 - expdamp * (1 - damp));
                            }
                        }
                        double rr3 = scale3 / (r * r2);
                        double rr5 = 3 * scale5 / (r * r2 * r2);
                        double duir = xr * duix + yr * duiy + zr * duiz;
                        double dukr = xr * dukx + yr * duky + zr * dukz;
                        double puir = xr * puix + yr * puiy + zr * puiz;
                        double pukr = xr * pukx + yr * puky + zr * pukz;
                        double fidx = -rr3 * dukx + rr5 * dukr * xr;
                        double fidy = -rr3 * duky + rr5 * dukr * yr;
                        double fidz = -rr3 * dukz + rr5 * dukr * zr;
                        double fkdx = -rr3 * duix + rr5 * duir * xr;
                        double fkdy = -rr3 * duiy + rr5 * duir * yr;
                        double fkdz = -rr3 * duiz + rr5 * duir * zr;
                        double fipx = -rr3 * pukx + rr5 * pukr * xr;
                        double fipy = -rr3 * puky + rr5 * pukr * yr;
                        double fipz = -rr3 * pukz + rr5 * pukr * zr;
                        double fkpx = -rr3 * puix + rr5 * puir * xr;
                        double fkpy = -rr3 * puiy + rr5 * puir * yr;
                        double fkpz = -rr3 * puiz + rr5 * puir * zr;
                        fx_local[i] += fidx;
                        fy_local[i] += fidy;
                        fz_local[i] += fidz;
                        fx_local[k] += fkdx;
                        fy_local[k] += fkdy;
                        fz_local[k] += fkdz;
                        fxp_local[i] += fipx;
                        fyp_local[i] += fipy;
                        fzp_local[i] += fipz;
                        fxp_local[k] += fkpx;
                        fyp_local[k] += fkpy;
                        fzp_local[k] += fkpz;
//                        System.out.println(i+" "+k+" "+fidx+" "+fkdx+" "+duix+" "+dukx);
                    }
                }
            }

            //add to the fields
            for (int i = 0; i < nAtoms; i++) {
                double fieldi[] = field1[i];
                double field2i[] = field2[i];
                fieldi[0] += fx_local[i];
                fieldi[1] += fy_local[i];
                fieldi[2] += fz_local[i];
                field2i[0] += fxp_local[i];
                field2i[1] += fyp_local[i];
                field2i[2] += fzp_local[i];
            }
            //Fix this later
            for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                for (int i = 0; i < (nAtoms - 1); i++) {
                    double ind[] = induced0[i];
                    double indp[] = inducedp0[i];
                    double xyzi[] = atoms[i].getXYZ(null);
                    double pdi = pdamp[i];
                    double pti = thole[i];
                    double duix = ind[0];
                    double duiy = ind[1];
                    double duiz = ind[2];
                    double puix = indp[0];
                    double puiy = indp[1];
                    double puiz = indp[2];
                    for (int index : ip11[i]) {
                        final double u1scale = 1;
                        mask_local[index] = u1scale;
                    }
                    for (int index : ip12[i]) {
                        final double u2scale = 1;
                        mask_local[index] = u2scale;
                    }
                    for (int index : ip13[i]) {
                        final double u3scale = 1;
                        mask_local[index] = u3scale;
                    }
                    for (int index : ip14[i]) {
                        final double u4scale = 1;
                        mask_local[index] = u4scale;
                    }
                    for (int k = i + 1; k < nAtoms; k++) {
                        double indk[] = induced0[k];
                        double indkp[] = inducedp0[k];
                        double xyzk[] = atoms[k].getXYZ(null);
                        dx_local[0] = xyzk[0] - xyzi[0];
                        dx_local[1] = xyzk[1] - xyzi[1];
                        dx_local[2] = xyzk[2] - xyzi[2];
                        final double r2 = crystal.image(dx_local);
                        if (r2 <= off2) {
                            double xr = dx_local[0];
                            double yr = dx_local[1];
                            double zr = dx_local[2];
                            double r = sqrt(r2);
                            double dukx = indk[0];
                            double duky = indk[1];
                            double dukz = indk[2];
                            double pukx = indkp[0];
                            double puky = indkp[1];
                            double pukz = indkp[2];
                            double scale3 = mask_local[k];
                            double scale5 = mask_local[k];
                            double damp = pdi * pdamp[k];
                            if (damp != 0) {
                                final double pgamma = min(pti, thole[k]);
                                damp = -pgamma * pow((r / damp), 3);
                                if (damp > -50) {
                                    double expdamp = exp(damp);
                                    scale3 = scale3 * (1 - expdamp);
                                    scale5 = scale5 * (1 - expdamp * (1 - damp));
                                }
                            }
                            double rr3 = scale3 / (r * r2);
                            double rr5 = 3 * scale5 / (r * r2 * r2);
                            double duir = xr * duix + yr * duiy + zr * duiz;
                            double dukr = xr * dukx + yr * duky + zr * dukz;
                            double puir = xr * puix + yr * puiy + zr * puiz;
                            double pukr = xr * pukx + yr * puky + zr * pukz;
                            double fidx = -rr3 * dukx + rr5 * dukr * xr;
                            double fidy = -rr3 * duky + rr5 * dukr * yr;
                            double fidz = -rr3 * dukz + rr5 * dukr * zr;
                            double fkdx = -rr3 * duix + rr5 * duir * xr;
                            double fkdy = -rr3 * duiy + rr5 * duir * yr;
                            double fkdz = -rr3 * duiz + rr5 * duir * zr;
                            double fipx = -rr3 * pukx + rr5 * pukr * xr;
                            double fipy = -rr3 * puky + rr5 * pukr * yr;
                            double fipz = -rr3 * pukz + rr5 * pukr * zr;
                            double fkpx = -rr3 * puix + rr5 * puir * xr;
                            double fkpy = -rr3 * puiy + rr5 * puir * yr;
                            double fkpz = -rr3 * puiz + rr5 * puir * zr;
                            fx_local[i] += fidx;
                            fy_local[i] += fidy;
                            fz_local[i] += fidz;
                            fx_local[k] += fkdx;
                            fy_local[k] += fkdy;
                            fz_local[k] += fkdz;
                            fxp_local[i] += fipx;
                            fyp_local[i] += fipy;
                            fzp_local[i] += fipz;
                            fxp_local[k] += fkpx;
                            fyp_local[k] += fkpy;
                            fzp_local[k] += fkpz;
                        }
                    }
                }
                //add to the fields
                for (int i = 0; i < nAtoms; i++) {
                    double fieldi[] = field1[i];
                    double field2i[] = field2[i];
                    fieldi[0] += fx_local[i];
                    fieldi[1] += fy_local[i];
                    fieldi[2] += fz_local[i];
                    field2i[0] += fxp_local[i];
                    field2i[1] += fyp_local[i];
                    field2i[2] += fzp_local[i];
                }
            }
            iter++;
            epsold = eps;
            eps = 0.0;
            double epsp = 0.0;
            for (int i = 0; i < nAtoms; i++) {
                final double ind[] = induced0[i];
                final double indp[] = inducedp0[i];
                final double direct[] = directDipole[i];
                final double directp[] = directDipolep[i];
                final double fieldi[] = field1[i];
                final double field2i[] = field2[i];
                final double polar = polarizability[i];
                for (int j = 0; j < 3; j++) {
                    double previous = ind[j];
                    double mutual = polar * fieldi[j];
                    ind[j] = direct[j] + mutual;
                    double delta = polsor * (ind[j] - previous);
                    ind[j] = previous + delta;
                    eps += delta * delta;
                    previous = indp[j];
                    mutual = polar * field2i[j];
                    indp[j] = directp[j] + mutual;
                    delta = polsor * (indp[j] - previous);
                    indp[j] = previous + delta;
                    epsp += delta * delta;
                }
            }
            //CHECK. NEED THIS?
//            try {
//    			parallelTeam.execute(expandInducedDipolesRegion);
//    		} catch (Exception e) {
//    			e.printStackTrace();
//    		}
            eps = max(eps, epsp);
            eps = MultipoleType.DEBYE * sqrt(eps / (double) nAtoms);
            cycleTime += System.nanoTime();

            if (print) {
                sb.append(String.format(
                        " %4d  %15.10f      %8.3f\n", iter, eps, cycleTime * toSeconds));
            }
            if (eps < poleps) {
                done = true;
            }
            if (eps > epsold) {
                done = true;
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = String.format("Fatal convergence failure: (%10.5f > %10.5f)\n", eps, epsold);
                logger.severe(message);
            }
            if (iter >= maxiter) {
                done = true;
                if (sb != null) {
                    logger.warning(sb.toString());
                }
                String message = String.format("Maximum iterations reached: (%d)\n", iter);
                logger.severe(message);
            }
        }

        if (print) {
            sb.append(String.format("\n Direct:                    %8.3f\n",
                    toSeconds * directTime));
            startTime = System.nanoTime() - startTime;
            sb.append(String.format(" SCF Total:                 %8.3f\n",
                    startTime * toSeconds));
            logger.info(sb.toString());
        }
        if (false) {
            sb = new StringBuffer();
            for (int i = 0; i < 100; i++) {
                sb.append(String.format(
                        "Induced Dipole  %d %15.8f %15.8f %15.8f\n", i + 1,
                        MultipoleType.DEBYE * inducedDipole[0][i][0],
                        MultipoleType.DEBYE * inducedDipole[0][i][1],
                        MultipoleType.DEBYE * inducedDipole[0][i][2]));
            }
            logger.info(sb.toString());
        }
    }

    /**
     * <p>
     * init_prms</p>
     */
    public void init_prms() {
        if (propyze) {
            dividempoles3();
        }
        //rotate multipoles into new frame and induce
        try {
            //Don't need to expand CoordinatesRegion. Decide later.
            parallelTeam.execute(expandCoordinatesRegion);
            parallelTeam.execute(rotateMultipolesRegion);
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (!propyze) {
            induce0a();
        }
    }

    /**
     * <p>
     * dividempoles3</p>
     */
    public void dividempoles3() {
        for (int i = 0; i < nAtoms; i++) {
            for (int j = 4; j < 10; j++) {
                localMultipole[i][j] = localMultipole[i][j] / 3;
            }
        }
    }

    //added by gchattree
    /**
     * <p>
     * potpoint</p>
     *
     * @param xyz an array of {@link java.lang.Double} objects.
     * @return a double.
     */
    public double potpoint(Double xyz[]) {
        double epot = 0;
        double xi = xyz[0];
        double yi = xyz[1];
        double zi = xyz[2];
        double mpole[][] = globalMultipole[0];
        double ind[][] = inducedDipole[0];
        double indp[][] = inducedDipolep[0];
        double ci;
        double dix;
        double diy;
        double diz;
        double qixx;
        double qiyy;
        double qizz;
        double qixy;
        double qixz;
        double qiyz;
        double uix;
        double uiy;
        double uiz;
        double qkx;
        double qky;
        double qkz;
        double scq;
        double scd;
        double scu;
        double globalMultipolei[];
        double inducedDipolei[];
        double inducedDipolepi[];
        double xr;
        double yr;
        double zr;
        double r2;
        double r;
        double rr1;
        double rr3;
        double rr5;
        double ei;
        double e;
        double em = 0;
        double ep = 0;
        double[] xyza;

        for (int i = 0; i < mpole.length; i++) {
            xyza = atoms[i].getXYZ(null);
            xr = xyza[0] - xi;
            yr = xyza[1] - yi;
            zr = xyza[2] - zi;
            r2 = xr * xr + yr * yr + zr * zr;
            r = sqrt(r2);
            globalMultipolei = mpole[i];
            inducedDipolei = ind[i];
            inducedDipolepi = indp[i];
            ci = globalMultipolei[t000];
            dix = globalMultipolei[t100];
            diy = globalMultipolei[t010];
            diz = globalMultipolei[t001];
            qixx = globalMultipolei[t200] / 3.0;
            qiyy = globalMultipolei[t020] / 3.0;
            qizz = globalMultipolei[t002] / 3.0;
            qixy = globalMultipolei[t110] / 3.0;
            qixz = globalMultipolei[t101] / 3.0;
            qiyz = globalMultipolei[t011] / 3.0;
            uix = inducedDipolei[0];
            uiy = inducedDipolei[1];
            uiz = inducedDipolei[2];

            qkx = qixx * xr + qixy * yr + qixz * zr;
            qky = qixy * xr + qiyy * yr + qiyz * zr;
            qkz = qixz * xr + qiyz * yr + qizz * zr;

            scd = dix * xr + diy * yr + diz * zr;
            scq = qkx * xr + qky * yr + qkz * zr;
            scu = uix * xr + uiy * yr + uiz * zr;

            rr1 = 1 / r;
            rr3 = rr1 / r2;
            rr5 = 3 * rr3 / r2;
            e = ci * rr1 - scd * rr3 + scq * rr5;
            ei = -scu * rr3;
            e = electric * e;
            ei = electric * ei;
            em = em + e;
            ep = ep + ei;

        }
        epot = ep + em;
        return epot;
    }

    private void setEwaldParameters() {
        off = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        aewald = forceField.getDouble(ForceFieldDouble.EWALD_ALPHA, ewaldCoefficient(off));
        off2 = off * off;
        alsq2 = 2.0 * aewald * aewald;
        piEwald = 1.0 / (sqrtPi * aewald);
        aewald3 = 4.0 / 3.0 * pow(aewald, 3.0) / sqrtPi;
        permanentSelfEnergy = permanentSelfEnergy();
    }

    private void checkCacheSize() {
        for (int i = 0; i < nSymm; i++) {
            for (int j = 0; j < nAtoms; j++) {
                int size = neighborLists[i][j].length;
                if (ewaldLists[i][j] == null || ewaldLists[i][j].length < size) {
                    ewaldLists[i][j] = new int[size];
                }
            }
        }
    }

    /**
     * Set the electrostatic lambda scaling factor. Atoms whose permanent
     * multipole and polarizability will be scaled are echoed to the logger.
     *
     * @param lambda Must satisfy greater than or equal to 0.0 and less than or
     * equal to 1.0.
     */
//    @Override
//    public void setLambda(double lambda) {
//        assert (lambda >= 0.0 && lambda <= 1.0);
//        this.lambda = lambda;
//        /**
//         * Log the new lambda value.
//         */
//        StringBuilder sb = new StringBuilder(" Multipoles and polarizabilities will be scaled for:\n");
//        boolean scaleAtoms = false;
//        for (int i = 0; i < nAtoms; i++) {
//            if (atoms[i].applyLambda()) {
//                scaleAtoms = true;
//                sb.append(" " + atoms[i].toShortString() + "\n");
//            }
//        }
//        if (scaleAtoms) {
//            logger.info(format(" Electrostatic lambda value is set to %8.3f", lambda));
//            logger.info(sb.toString());
//        } else {
//            logger.warning(" No atoms are selected for multipole and polarizability scaling.\n");
//        }
//    }
    /**
     * Get the current lambda scale value.
     *
     * @return lambda
     */
//    @Override
//    public double getLambda() {
//        return lambda;
//    }
    /**
     * Calculate the PME electrostatic energy.
     *
     * @param gradient If <code>true</code>, the gradient will be calculated.
     * @param print If <code>true</code>, extra logging is enabled.
     * @return return the total electrostatic energy (permanent + polarization).
     */
//    public double energy(boolean gradient, boolean print) {
//        /**
//         * Initialize the energy components.
//         */
//        double eself = 0.0;
//        double erecip = 0.0;
//        double ereal = 0.0;
//        double eselfi = 0.0;
//        double erecipi = 0.0;
//        double ereali = 0.0;
//        multipoleEnergy = 0.0;
//        polarizationEnergy = 0.0;
//        interactions = 0;
//        realSpaceTime = 0;
//        reciprocalSpaceTime = 0;
//
//        /**
//         * Initialize the coordinates.
//         */
//        double x[] = coordinates[0][0];
//        double y[] = coordinates[0][1];
//        double z[] = coordinates[0][2];
//        for (int i = 0; i < nAtoms; i++) {
//            double xyz[] = atoms[i].getXYZ();
//            x[i] = xyz[0];
//            y[i] = xyz[1];
//            z[i] = xyz[2];
//        }
//
//        checkCacheSize();
//
//        /**
//         * Initialize the gradient accumulation arrays.
//         */
//        if (gradient) {
//            for (int j = 0; j < nAtoms; j++) {
//                sharedGrad[0].set(j, 0.0);
//                sharedGrad[1].set(j, 0.0);
//                sharedGrad[2].set(j, 0.0);
//                sharedTorque[0].set(j, 0.0);
//                sharedTorque[1].set(j, 0.0);
//                sharedTorque[2].set(j, 0.0);
//            }
//        }
//
//        /**
//         * Find the permanent multipole potential and its gradients.
//         */
//        try {
//            parallelTeam.execute(expandCoordinatesRegion);
//            parallelTeam.execute(rotateMultipolesRegion);
//
//            bsplineTime = -System.nanoTime();
//            reciprocalSpace.computeBSplines();
//            bsplineTime += System.nanoTime();
//
//            densityTime = -System.nanoTime();
//            reciprocalSpace.splinePermanentMultipoles(globalMultipole, null);
//            densityTime += System.nanoTime();
//            /**
//             * Here the real space contribution to the field is calculated at
//             * the same time the reciprocal space convolution is being done.
//             * This is useful since the reciprocal space convolution
//             * (the 3D FFT and inverse FFT) do not parallelize well.
//             */
//            realAndFFTTime = -System.nanoTime();
//            sectionTeam.execute(permanentFieldRegion);
//            realAndFFTTime += System.nanoTime();
//
//            phiTime = -System.nanoTime();
//            reciprocalSpace.computePermanentPhi(cartesianMultipolePhi);
//            phiTime += System.nanoTime();
//        } catch (Exception e) {
//            String message = "Fatal exception computing the permanent multipole field.\n";
//            logger.log(Level.SEVERE, message, e);
//        }
//
//        /**
//         * Do the self consistent field calculation.
//         */
//        selfConsistentField(logger.isLoggable(Level.FINE));
//
//        if (logger.isLoggable(Level.FINE)) {
//            StringBuilder sb = new StringBuilder();
//            sb.append(format("\n b-Spline:   %8.3f (sec)\n", bsplineTime * toSeconds));
//            sb.append(format(" Density:    %8.3f (sec)\n", densityTime * toSeconds));
//            sb.append(format(" Real + FFT: %8.3f (sec)\n", realAndFFTTime * toSeconds));
//            sb.append(format(" Phi:        %8.3f (sec)\n", phiTime * toSeconds));
//            logger.fine(sb.toString());
//        }
//
//        /**
//         * The self energy of the permanent multipoles is constant.
//         */
//        eself = permanentSelfEnergy;
//        interactions = nAtoms;
//
//        /**
//         * The energy of the permanent multipoles in their reciprocal space
//         * potential. This potential was computed prior to the SCF.
//         */
//        erecip = permanentReciprocalSpaceEnergy(gradient);
//
//        /**
//         * Find the total real space energy. This includes the permanent
//         * multipoles in their own real space potential and the interaction of
//         * permanent multipoles with induced dipoles.
//         */
//        long time = System.nanoTime();
//        try {
//            realSpaceEnergyRegion.setGradient(gradient);
//            parallelTeam.execute(realSpaceEnergyRegion);
//            ereal = realSpaceEnergyRegion.getPermanentEnergy();
//            ereali = realSpaceEnergyRegion.getPolarizationEnergy();
//            interactions += realSpaceEnergyRegion.getInteractions();
//        } catch (Exception e) {
//            String message = "Exception computing the real space energy.\n";
//            logger.log(Level.SEVERE, message, e);
//        }
//        time = System.nanoTime() - time;
//        realSpaceTime += time;
//
//        if (polarization != Polarization.NONE) {
//            /**
//             * The induced dipole self energy due to interaction with the
//             * permanent dipole.
//             */
//            eselfi = inducedDipoleSelfEnergy(gradient);
    /**
     * Calculate the PME electrostatic energy.
     *
     * @return return the total electrostatic energy (permanent + polarization).
     */
//    public double energy(boolean gradient, boolean print) {
//        /**
//         * Initialize the energy components.
//         */
//        double eself = 0.0;
//        double erecip = 0.0;
//        double ereal = 0.0;
//        double eselfi = 0.0;
//        double erecipi = 0.0;
//        double ereali = 0.0;
//        multipoleEnergy = 0.0;
//        polarizationEnergy = 0.0;
//        interactions = 0;
//        realSpaceTime = 0;
//        reciprocalSpaceTime = 0;
//
//        /**
//         * Initialize the coordinates.
//         */
//        double x[] = coordinates[0][0];
//        double y[] = coordinates[0][1];
//        double z[] = coordinates[0][2];
//        for (int i = 0; i < nAtoms; i++) {
//            double xyz[] = atoms[i].getXYZ();
//            x[i] = xyz[0];
//            y[i] = xyz[1];
//            z[i] = xyz[2];
//        }
//
//        checkCacheSize();
//
//        /**
//         * Initialize the gradient accumulation arrays.
//         */
//        if (gradient) {
//            for (int j = 0; j < nAtoms; j++) {
//                sharedGrad[0].set(j, 0.0);
//                sharedGrad[1].set(j, 0.0);
//                sharedGrad[2].set(j, 0.0);
//                sharedTorque[0].set(j, 0.0);
//                sharedTorque[1].set(j, 0.0);
//                sharedTorque[2].set(j, 0.0);
//            }
//        }
//
//        /**
//         * Find the permanent multipole potential and its gradients.
//         */
//        try {
//            parallelTeam.execute(expandCoordinatesRegion);
//            parallelTeam.execute(rotateMultipolesRegion);
//
//            bsplineTime = -System.nanoTime();
//            reciprocalSpace.computeBSplines();
//            bsplineTime += System.nanoTime();
//
//            densityTime = -System.nanoTime();
//            reciprocalSpace.splinePermanentMultipoles(globalMultipole, null);
//            densityTime += System.nanoTime();
//            /**
//             * Here the real space contribution to the field is calculated at
//             * the same time the reciprocal space convolution is being done.
//             * This is useful since the reciprocal space convolution
//             * (the 3D FFT and inverse FFT) do not parallelize well.
//             */
//            realAndFFTTime = -System.nanoTime();
//            sectionTeam.execute(permanentFieldRegion);
//            realAndFFTTime += System.nanoTime();
//
//            phiTime = -System.nanoTime();
//            reciprocalSpace.computePermanentPhi(cartesianMultipolePhi);
//            phiTime += System.nanoTime();
//        } catch (Exception e) {
//            String message = "Fatal exception computing the permanent multipole field.\n";
//            logger.log(Level.SEVERE, message, e);
//        }
//
//        /**
//         * Do the self consistent field calculation.
//         */
//        selfConsistentField(logger.isLoggable(Level.FINE));
//
//        if (logger.isLoggable(Level.FINE)) {
//            StringBuilder sb = new StringBuilder();
//            sb.append(format("\n b-Spline:   %8.3f (sec)\n", bsplineTime * toSeconds));
//            sb.append(format(" Density:    %8.3f (sec)\n", densityTime * toSeconds));
//            sb.append(format(" Real + FFT: %8.3f (sec)\n", realAndFFTTime * toSeconds));
//            sb.append(format(" Phi:        %8.3f (sec)\n", phiTime * toSeconds));
//            logger.fine(sb.toString());
//        }
//
//        /**
//         * The self energy of the permanent multipoles is constant.
//         */
//        eself = permanentSelfEnergy;
//        interactions = nAtoms;
//
//        /**
//         * The energy of the permanent multipoles in their reciprocal space
//         * potential. This potential was computed prior to the SCF.
//         */
//        erecip = permanentReciprocalSpaceEnergy(gradient);
//
//        /**
//         * Find the total real space energy. This includes the permanent
//         * multipoles in their own real space potential and the interaction of
//         * permanent multipoles with induced dipoles.
//         */
//        long time = System.nanoTime();
//        try {
//            realSpaceEnergyRegion.setGradient(gradient);
//            parallelTeam.execute(realSpaceEnergyRegion);
//            ereal = realSpaceEnergyRegion.getPermanentEnergy();
//            ereali = realSpaceEnergyRegion.getPolarizationEnergy();
//            interactions += realSpaceEnergyRegion.getInteractions();
//        } catch (Exception e) {
//            String message = "Exception computing the real space energy.\n";
//            logger.log(Level.SEVERE, message, e);
//        }
//        time = System.nanoTime() - time;
//        realSpaceTime += time;
//
//        if (polarization != Polarization.NONE) {
//            /**
//             * The induced dipole self energy due to interaction with the
//             * permanent dipole.
//             */
//            eselfi = inducedDipoleSelfEnergy(gradient);
//            /**
//             * The energy of the permanent multipoles in the induced dipole
//             * reciprocal potential.
//             */
//            time = System.nanoTime();
//            erecipi = inducedDipoleReciprocalSpaceEnergy(gradient);
//            time = System.nanoTime() - time;
//            reciprocalSpaceTime += time;
//        }
//        if (logger.isLoggable(Level.FINE)) {
//            StringBuilder sb = new StringBuilder();
//            sb.append(format("\n Total Time =    Real +   Recip (sec)\n"));
//            sb.append(format("   %8.3f =%8.3f +%8.3f\n", toSeconds * (realSpaceTime + reciprocalSpaceTime), toSeconds * realSpaceTime,
//                             toSeconds * reciprocalSpaceTime));
//            sb.append(format(" Multipole Self-Energy:   %16.8f\n", eself));
//            sb.append(format(" Multipole Reciprocal:    %16.8f\n", erecip));
//            sb.append(format(" Multipole Real Space:    %16.8f\n", ereal));
//            sb.append(format(" Polarization Self-Energy:%16.8f\n", eselfi));
//            sb.append(format(" Polarization Reciprocal: %16.8f\n", erecipi));
//            sb.append(format(" Polarization Real Space: %16.8f\n", ereali));
//            logger.fine(sb.toString());
//        }
//        // Collect energy terms.
//        multipoleEnergy = eself + erecip + ereal;
//        polarizationEnergy = eselfi + erecipi + ereali;
//        // Add electrostatic gradient to total atomic gradient.
//        if (gradient) {
//            // Convert torques to forces.
//            try {
//                parallelTeam.execute(torqueRegion);
//            } catch (Exception e) {
//                String message = "Exception calculating torques.";
//                logger.log(Level.SEVERE, message, e);
//            }
//            for (int i = 0; i < nAtoms; i++) {
//                atoms[i].addToXYZGradient(sharedGrad[0].get(i), sharedGrad[1].get(i), sharedGrad[2].get(i));
//            }
//        }
//        return multipoleEnergy + polarizationEnergy;
//    }
    public int getInteractions() {
        return interactions;
    }

    /**
     * <p>
     * getPermanentEnergy</p>
     *
     * @return a double.
     */
    public double getPermanentEnergy() {
        return multipoleEnergy;
    }

    /**
     * <p>
     * Getter for the field <code>polarizationEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    /**
     * <p>
     * getGradients</p>
     *
     * @param grad an array of double.
     */
    public void getGradients(double grad[][]) {
        for (int i = 0; i < nAtoms; i++) {
            grad[0][i] = sharedGrad[0].get(i);
            grad[1][i] = sharedGrad[1].get(i);
            grad[2][i] = sharedGrad[2].get(i);
        }
    }

//    private void selfConsistentField(boolean print) {
//        long startTime = System.nanoTime();
//        /**
//         * Initialize the electric field to the direct field.
//         */
//        for (int i = 0; i < nAtoms; i++) {
//            double fieldi[] = field1[i];
//            double fieldpi[] = field2[i];
//            double mpolei[] = globalMultipole[0][i];
//            double phii[] = cartesianMultipolePhi[i];
//            double fx = aewald3 * mpolei[t100] - phii[t100];
//            double fy = aewald3 * mpolei[t010] - phii[t010];
//            double fz = aewald3 * mpolei[t001] - phii[t001];
//            fieldi[0] += fx;
//            fieldi[1] += fy;
//            fieldi[2] += fz;
//            fieldpi[0] += fx;
//            fieldpi[1] += fy;
//            fieldpi[2] += fz;
//        }
//
//        if (polarization == Polarization.NONE) {
//            for (int i = 0; i < nAtoms; i++) {
//                inducedDipole[0][i][0] = 0.0;
//                inducedDipole[0][i][1] = 0.0;
//                inducedDipole[0][i][2] = 0.0;
//                inducedDipolep[0][i][0] = 0.0;
//                inducedDipolep[0][i][1] = 0.0;
//                inducedDipolep[0][i][2] = 0.0;
//            }
//            try {
//                parallelTeam.execute(expandInducedDipolesRegion);
//            } catch (Exception e) {
//                String message = "Exception expanding induced dipoles.";
//                logger.log(Level.SEVERE, message, e);
//            }
//            return;
//        }
//        /**
//         * Set the induced dipoles to the polarizability times the direct field.
//         */
//        final double induced0[][] = inducedDipole[0];
//        final double inducedp0[][] = inducedDipolep[0];
//        for (int i = 0; i < nAtoms; i++) {
//            final double polar = polarizability[i];
//            final double fieldi[] = field1[i];
//            final double ind[] = induced0[i];
//            final double directi[] = directDipole[i];
//            ind[0] = polar * fieldi[0];
//            ind[1] = polar * fieldi[1];
//            ind[2] = polar * fieldi[2];
//            directi[0] = ind[0];
//            directi[1] = ind[1];
//            directi[2] = ind[2];
//            final double field2i[] = field2[i];
//            final double inp[] = inducedp0[i];
//            final double directpi[] = directDipolep[i];
//            inp[0] = polar * field2i[0];
//            inp[1] = polar * field2i[1];
//            inp[2] = polar * field2i[2];
//            directpi[0] = inp[0];
//            directpi[1] = inp[1];
//            directpi[2] = inp[2];
//        }
//
//        try {
//            parallelTeam.execute(expandInducedDipolesRegion);
//        } catch (Exception e) {
//            String message = "Exception expanding induced dipoles.";
//            logger.log(Level.SEVERE, message, e);
//        }
//
//        if (polarization == Polarization.MUTUAL) {
//            StringBuilder sb = null;
//            long directTime = System.nanoTime() - startTime;
//            if (print) {
//                sb = new StringBuilder(
//                        "\n SELF-CONSISTENT FIELD\n"
//                        + " Iter     RMS Change (Debyes)   Time\n");
//            }
//            boolean done = false;
//            int maxiter = 1000;
//            int iter = 0;
//            double eps = 100.0;
//            double epsold;
//            while (!done) {
//                long cycleTime = -System.nanoTime();
//                /**
//                 * Find the induced dipole field.
//                 */
//                try {
//                    densityTime -= System.nanoTime();
//                    reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipolep, null);
//                    densityTime += System.nanoTime();
//
//                    realAndFFTTime -= System.nanoTime();
//                    sectionTeam.execute(inducedDipoleFieldRegion);
//                    realAndFFTTime += System.nanoTime();
//
//                    phiTime -= System.nanoTime();
//                    reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolepPhi);
//                    phiTime += System.nanoTime();
//                } catch (Exception e) {
//                    String message = "Fatal exception computing the induced dipole field.\n";
//                    logger.log(Level.SEVERE, message, e);
//                }
//
//                /**
//                 * Add the self and reciprocal space fields to the
//                 * real space field.
//                 */
//                for (int i = 0; i < nAtoms; i++) {
//                    double fieldi[] = field1[i];
//                    double fieldpi[] = field2[i];
//                    double dipolei[] = induced0[i];
//                    double dipolepi[] = inducedp0[i];
//                    final double phii[] = cartesianDipolePhi[i];
//                    final double phipi[] = cartesianDipolepPhi[i];
//                    fieldi[0] += aewald3 * dipolei[0];
//                    fieldi[1] += aewald3 * dipolei[1];
//                    fieldi[2] += aewald3 * dipolei[2];
//                    fieldi[0] -= phii[t100];
//                    fieldi[1] -= phii[t010];
//                    fieldi[2] -= phii[t001];
//                    fieldpi[0] += aewald3 * dipolepi[0];
//                    fieldpi[1] += aewald3 * dipolepi[1];
//                    fieldpi[2] += aewald3 * dipolepi[2];
//                    fieldpi[0] -= phipi[t100];
//                    fieldpi[1] -= phipi[t010];
//                    fieldpi[2] -= phipi[t001];
//                }
//
//                /**
//                 * Apply Successive Over-Relaxation (SOR) and
//                 * check for convergence of the SCF.
//                 */
//                iter++;
//                epsold = eps;
//                eps = 0.0;
//                double epsp = 0.0;
//                for (int i = 0; i < nAtoms; i++) {
//                    final double ind[] = induced0[i];
//                    final double indp[] = inducedp0[i];
//                    final double direct[] = directDipole[i];
//                    final double directp[] = directDipolep[i];
//                    final double fieldi[] = field1[i];
//                    final double field2i[] = field2[i];
//                    final double polar = polarizability[i];
//                    for (int j = 0; j < 3; j++) {
//                        double previous = ind[j];
//                        double mutual = polar * fieldi[j];
//                        ind[j] = direct[j] + mutual;
//                        double delta = polsor * (ind[j] - previous);
//                        ind[j] = previous + delta;
//                        eps += delta * delta;
//                        previous = indp[j];
//                        mutual = polar * field2i[j];
//                        indp[j] = directp[j] + mutual;
//                        delta = polsor * (indp[j] - previous);
//                        indp[j] = previous + delta;
//                        epsp += delta * delta;
//                    }
//                }
//                try {
//                    parallelTeam.execute(expandInducedDipolesRegion);
//                } catch (Exception e) {
//                    String message = "Exception expanding induced dipoles.";
//                    logger.log(Level.SEVERE, message, e);
//                }
//                eps = max(eps, epsp);
//                eps = MultipoleType.DEBYE * sqrt(eps / (double) nAtoms);
//                cycleTime += System.nanoTime();
//
//                if (print) {
//                    sb.append(format(
//                            " %4d  %15.10f      %8.3f\n", iter, eps, cycleTime * toSeconds));
//                }
//                if (eps < poleps) {
//                    done = true;
//                }
//                if (eps > epsold) {
//                    if (sb != null) {
//                        logger.warning(sb.toString());
//                    }
//                    String message = format("Fatal convergence failure: (%10.5f > %10.5f)\n", eps, epsold);
//                    logger.severe(message);
//                }
//                if (iter >= maxiter) {
//                    if (sb != null) {
//                        logger.warning(sb.toString());
//                    }
//                    String message = format("Maximum iterations reached: (%d)\n", iter);
//                    logger.severe(message);
//                }
//            }
//            if (print) {
//                sb.append(format("\n Direct:                    %8.3f\n",
//                                 toSeconds * directTime));
//                startTime = System.nanoTime() - startTime;
//                sb.append(format(" SCF Total:                 %8.3f\n",
//                                 startTime * toSeconds));
//                logger.info(sb.toString());
//            }
//            if (false) {
//                sb = new StringBuilder();
//                for (int i = 0; i < 100; i++) {
//                    sb.append(format(
//                            "Induced Dipole  %d %15.8f %15.8f %15.8f\n", i + 1,
//                            MultipoleType.DEBYE * inducedDipole[0][i][0],
//                            MultipoleType.DEBYE * inducedDipole[0][i][1],
//                            MultipoleType.DEBYE * inducedDipole[0][i][2]));
//                }
//                logger.info(sb.toString());
//            }
//        }
//    }
    private double permanentSelfEnergy() {
        double e = 0.0;
        double term = 2.0 * aewald * aewald;
        double fterm = -electric * aewald / sqrtPi;
        for (int i = 0; i < nAtoms; i++) {
            double in[] = localMultipole[i];
            double cii = in[t000] * in[t000];
            double dii = in[t100] * in[t100] + in[t010] * in[t010] + in[t001] * in[t001];
            double qii = in[t200] * in[t200] + in[t020] * in[t020] + in[t002] * in[t002] + 2.0 * (in[t110] * in[t110] + in[t101] * in[t101] + in[t011] * in[t011]);
            e += fterm * (cii + term * (dii / 3.0 + 2.0 * term * qii / 45.0));
        }
        return e;
    }

    /**
     * The Permanent Field Region should be executed by a ParallelTeam with
     * exactly 2 threads. The Real Space and Reciprocal Space Sections will be
     * run concurrently, each with the number of threads defined by their
     * respective ParallelTeam instances.
     */
//    private class PermanentFieldRegion extends ParallelRegion {
//
//        private PermanentRealSpaceFieldSection permanentRealSpaceFieldSection;
//        private PermanentReciprocalSection permanentReciprocalSection;
//
//        public PermanentFieldRegion(ParallelTeam pt) {
//            permanentRealSpaceFieldSection = new PermanentRealSpaceFieldSection(pt);
//            permanentReciprocalSection = new PermanentReciprocalSection();
//        }
//
//        @Override
//        public void run() {
//            try {
//                execute(permanentRealSpaceFieldSection, permanentReciprocalSection);
//            } catch (Exception e) {
//                String message = "Fatal exception computing the permanent multipole field.\n";
//                logger.log(Level.SEVERE, message, e);
//            }
//        }
//    }
    /**
     * Computes the Permanent Multipole Real Space Field.
     */
    private class PermanentRealSpaceFieldSection extends ParallelSection {

        private final PermanentRealSpaceFieldRegion permanentRealSpaceFieldRegion;
        private final ParallelTeam pt;

        public PermanentRealSpaceFieldSection(ParallelTeam pt) {
            this.pt = pt;
            int nt = pt.getThreadCount();
            permanentRealSpaceFieldRegion = new PermanentRealSpaceFieldRegion(nt);
        }

        @Override
        public void run() {
            try {
                long time = -System.nanoTime();
                pt.execute(permanentRealSpaceFieldRegion);
                permanentRealSpaceFieldRegion.setField(field1, field2);
                time += System.nanoTime();
                realSpaceTime += time;
            } catch (Exception e) {
                String message = "Fatal exception computing the real space field.\n";
                logger.log(Level.SEVERE, message, e);
            }
        }
    }

    /**
     * Compute the permanent multipole reciprocal space contribution to the
     * electric potential, field, etc. using the number of threads specified by
     * the ParallelTeam used to construct the ReciprocalSpace instance.
     */
//    private class PermanentReciprocalSection extends ParallelSection {
//
//        @Override
//        public void run() {
//            reciprocalSpace.permanentMultipoleConvolution();
//        }
//    }
    /**
     * The Induced Dipole Field Region should be executed by a ParallelTeam with
     * exactly 2 threads. The Real Space and Reciprocal Space Sections will be
     * run concurrently, each with the number of threads defined by their
     * respective ParallelTeam instances.
     */
//    private class InducedDipoleFieldRegion extends ParallelRegion {
//
//        private InducedDipoleRealSpaceFieldSection inducedRealSpaceFieldSection;
//        private InducedDipoleReciprocalFieldSection inducedReciprocalFieldSection;
//
//        public InducedDipoleFieldRegion(ParallelTeam pt) {
//            inducedRealSpaceFieldSection = new InducedDipoleRealSpaceFieldSection(pt);
//            inducedReciprocalFieldSection = new InducedDipoleReciprocalFieldSection();
//        }
//
//        @Override
//        public void run() {
//            try {
//                execute(inducedRealSpaceFieldSection, inducedReciprocalFieldSection);
//            } catch (Exception e) {
//                String message = "Fatal exception computing the induced dipole field.\n";
//                logger.log(Level.SEVERE, message, e);
//            }
//        }
//    }
//
//    private class InducedDipoleRealSpaceFieldSection extends ParallelSection {
//
//        private final InducedDipoleRealSpaceFieldRegion polarizationRealSpaceFieldRegion;
//        private final ParallelTeam pt;
//
//        public InducedDipoleRealSpaceFieldSection(ParallelTeam pt) {
//            this.pt = pt;
//            int nt = pt.getThreadCount();
//            polarizationRealSpaceFieldRegion = new InducedDipoleRealSpaceFieldRegion(nt);
//        }
//
//        @Override
//        public void run() {
//            try {
//                long time = System.nanoTime();
//                pt.execute(polarizationRealSpaceFieldRegion);
//                polarizationRealSpaceFieldRegion.setField(field1, field2);
//                time = System.nanoTime() - time;
//                realSpaceTime += time;
//            } catch (Exception e) {
//                String message = "Fatal exception computing the real space field.\n";
//                logger.log(Level.SEVERE, message, e);
//            }
//        }
//    }
//    private class InducedDipoleReciprocalFieldSection extends ParallelSection {
//
//        @Override
//        public void run() {
//            reciprocalSpace.inducedDipoleConvolution();
//        }
//    }
//    private double permanentReciprocalSpaceEnergy(boolean gradient) {
//        double erecip = 0.0;
//        final double pole[][] = globalMultipole[0];
//        final double fpole[][] = reciprocalSpace.getFracMultipoles();
//        final double fractionalMultipolePhi[][] = reciprocalSpace.getFracMultipolePhi();
//        final double nfftX = reciprocalSpace.getXDim();
//        final double nfftY = reciprocalSpace.getYDim();
//        final double nfftZ = reciprocalSpace.getZDim();
//        for (int i = 0; i < nAtoms; i++) {
//            final double phi[] = cartesianMultipolePhi[i];
//            final double fPhi[] = fractionalMultipolePhi[i];
//            final double mpole[] = pole[i];
//            final double fmpole[] = fpole[i];
//            double e = fmpole[t000] * fPhi[t000] + fmpole[t100] * fPhi[t100]
//                       + fmpole[t010] * fPhi[t010] + fmpole[t001] * fPhi[t001]
//                       + fmpole[t200] * fPhi[t200] + fmpole[t020] * fPhi[t020]
//                       + fmpole[t002] * fPhi[t002] + fmpole[t110] * fPhi[t110]
//                       + fmpole[t101] * fPhi[t101] + fmpole[t011] * fPhi[t011];
//            erecip += e;
//            if (gradient) {
//                double gx = fmpole[t000] * fPhi[t100] + fmpole[t100] * fPhi[t200] + fmpole[t010] * fPhi[t110] + fmpole[t001] * fPhi[t101] + fmpole[t200] * fPhi[t300] + fmpole[t020] * fPhi[t120] + fmpole[t002] * fPhi[t102] + fmpole[t110] * fPhi[t210] + fmpole[t101] * fPhi[t201] + fmpole[t011] * fPhi[t111];
//                double gy = fmpole[t000] * fPhi[t010] + fmpole[t100] * fPhi[t110] + fmpole[t010] * fPhi[t020] + fmpole[t001] * fPhi[t011] + fmpole[t200] * fPhi[t210] + fmpole[t020] * fPhi[t030] + fmpole[t002] * fPhi[t012] + fmpole[t110] * fPhi[t120] + fmpole[t101] * fPhi[t111] + fmpole[t011] * fPhi[t021];
//                double gz = fmpole[t000] * fPhi[t001] + fmpole[t100] * fPhi[t101] + fmpole[t010] * fPhi[t011] + fmpole[t001] * fPhi[t002] + fmpole[t200] * fPhi[t201] + fmpole[t020] * fPhi[t021] + fmpole[t002] * fPhi[t003] + fmpole[t110] * fPhi[t111] + fmpole[t101] * fPhi[t102] + fmpole[t011] * fPhi[t012];
//                gx *= nfftX;
//                gy *= nfftY;
//                gz *= nfftZ;
//                final double recip[][] = crystal.getUnitCell().A;
//                final double dfx = recip[0][0] * gx + recip[0][1] * gy + recip[0][2] * gz;
//                final double dfy = recip[1][0] * gx + recip[1][1] * gy + recip[1][2] * gz;
//                final double dfz = recip[2][0] * gx + recip[2][1] * gy + recip[2][2] * gz;
//                // Compute dipole torques
//                double tqx = -mpole[t010] * phi[t001] + mpole[t001] * phi[t010];
//                double tqy = -mpole[t001] * phi[t100] + mpole[t100] * phi[t001];
//                double tqz = -mpole[t100] * phi[t010] + mpole[t010] * phi[t100];
//                // Compute quadrupole torques
//                tqx -= 2.0 / 3.0 * (mpole[t110] * phi[t101] + mpole[t020] * phi[t011] + mpole[t011] * phi[t002] - mpole[t101] * phi[t110] - mpole[t011] * phi[t020] - mpole[t002] * phi[t011]);
//                tqy -= 2.0 / 3.0 * (mpole[t101] * phi[t200] + mpole[t011] * phi[t110] + mpole[t002] * phi[t101] - mpole[t200] * phi[t101] - mpole[t110] * phi[t011] - mpole[t101] * phi[t002]);
//                tqz -= 2.0 / 3.0 * (mpole[t200] * phi[t110] + mpole[t110] * phi[t020] + mpole[t101] * phi[t011] - mpole[t110] * phi[t200] - mpole[t020] * phi[t110] - mpole[t011] * phi[t101]);
//                sharedGrad[0].addAndGet(i, electric * dfx);
//                sharedGrad[1].addAndGet(i, electric * dfy);
//                sharedGrad[2].addAndGet(i, electric * dfz);
//                sharedTorque[0].addAndGet(i, electric * tqx);
//                sharedTorque[1].addAndGet(i, electric * tqy);
//                sharedTorque[2].addAndGet(i, electric * tqz);
//            }
//        }
//        erecip = 0.5 * electric * erecip;
//        return erecip;
//    }
//    private class InducedDipoleFieldRegion extends ParallelRegion {
//
//        private InducedDipoleRealSpaceFieldSection inducedRealSpaceFieldSection;
//        private InducedDipoleReciprocalFieldSection inducedReciprocalFieldSection;
//
//        public InducedDipoleFieldRegion(ParallelTeam pt) {
//            inducedRealSpaceFieldSection = new InducedDipoleRealSpaceFieldSection(pt);
//            inducedReciprocalFieldSection = new InducedDipoleReciprocalFieldSection();
//        }
//
//        @Override
//        public void run() {
//            try {
//                execute(inducedRealSpaceFieldSection, inducedReciprocalFieldSection);
//            } catch (Exception e) {
//                String message = "Fatal exception computing the induced dipole field.\n";
//                logger.log(Level.SEVERE, message, e);
//            }
//        }
//    }
//
//    private class InducedDipoleRealSpaceFieldSection extends ParallelSection {
//
//        private final InducedDipoleRealSpaceFieldRegion polarizationRealSpaceFieldRegion;
//        private final ParallelTeam pt;
//
//        public InducedDipoleRealSpaceFieldSection(ParallelTeam pt) {
//            this.pt = pt;
//            int nt = pt.getThreadCount();
//            polarizationRealSpaceFieldRegion = new InducedDipoleRealSpaceFieldRegion(nt);
//        }
//
//        @Override
//        public void run() {
//            try {
//                long time = System.nanoTime();
//                pt.execute(polarizationRealSpaceFieldRegion);
//                polarizationRealSpaceFieldRegion.setField(field1, field2);
//                time = System.nanoTime() - time;
//                realSpaceTime += time;
//            } catch (Exception e) {
//                String message = "Fatal exception computing the real space field.\n";
//                logger.log(Level.SEVERE, message, e);
//            }
//        }
//    }
//    private class InducedDipoleReciprocalFieldSection extends ParallelSection {
//
//        @Override
//        public void run() {
//            reciprocalSpace.inducedDipoleConvolution();
//        }
//    }
//    private double permanentReciprocalSpaceEnergy(boolean gradient) {
//        double erecip = 0.0;
//        final double pole[][] = globalMultipole[0];
//        final double fpole[][] = reciprocalSpace.getFracMultipoles();
//        final double fractionalMultipolePhi[][] = reciprocalSpace.getFracMultipolePhi();
//        final double nfftX = reciprocalSpace.getXDim();
//        final double nfftY = reciprocalSpace.getYDim();
//        final double nfftZ = reciprocalSpace.getZDim();
//        for (int i = 0; i < nAtoms; i++) {
//            final double phi[] = cartesianMultipolePhi[i];
//            final double fPhi[] = fractionalMultipolePhi[i];
//            final double mpole[] = pole[i];
//            final double fmpole[] = fpole[i];
//            double e = fmpole[t000] * fPhi[t000] + fmpole[t100] * fPhi[t100]
//                       + fmpole[t010] * fPhi[t010] + fmpole[t001] * fPhi[t001]
//                       + fmpole[t200] * fPhi[t200] + fmpole[t020] * fPhi[t020]
//                       + fmpole[t002] * fPhi[t002] + fmpole[t110] * fPhi[t110]
//                       + fmpole[t101] * fPhi[t101] + fmpole[t011] * fPhi[t011];
//            erecip += e;
//            if (gradient) {
//                double gx = fmpole[t000] * fPhi[t100] + fmpole[t100] * fPhi[t200] + fmpole[t010] * fPhi[t110] + fmpole[t001] * fPhi[t101] + fmpole[t200] * fPhi[t300] + fmpole[t020] * fPhi[t120] + fmpole[t002] * fPhi[t102] + fmpole[t110] * fPhi[t210] + fmpole[t101] * fPhi[t201] + fmpole[t011] * fPhi[t111];
//                double gy = fmpole[t000] * fPhi[t010] + fmpole[t100] * fPhi[t110] + fmpole[t010] * fPhi[t020] + fmpole[t001] * fPhi[t011] + fmpole[t200] * fPhi[t210] + fmpole[t020] * fPhi[t030] + fmpole[t002] * fPhi[t012] + fmpole[t110] * fPhi[t120] + fmpole[t101] * fPhi[t111] + fmpole[t011] * fPhi[t021];
//                double gz = fmpole[t000] * fPhi[t001] + fmpole[t100] * fPhi[t101] + fmpole[t010] * fPhi[t011] + fmpole[t001] * fPhi[t002] + fmpole[t200] * fPhi[t201] + fmpole[t020] * fPhi[t021] + fmpole[t002] * fPhi[t003] + fmpole[t110] * fPhi[t111] + fmpole[t101] * fPhi[t102] + fmpole[t011] * fPhi[t012];
//                gx *= nfftX;
//                gy *= nfftY;
//                gz *= nfftZ;
//                final double recip[][] = crystal.getUnitCell().A;
//                final double dfx = recip[0][0] * gx + recip[0][1] * gy + recip[0][2] * gz;
//                final double dfy = recip[1][0] * gx + recip[1][1] * gy + recip[1][2] * gz;
//                final double dfz = recip[2][0] * gx + recip[2][1] * gy + recip[2][2] * gz;
//                // Compute dipole torques
//                double tqx = -mpole[t010] * phi[t001] + mpole[t001] * phi[t010];
//                double tqy = -mpole[t001] * phi[t100] + mpole[t100] * phi[t001];
//                double tqz = -mpole[t100] * phi[t010] + mpole[t010] * phi[t100];
//                // Compute quadrupole torques
//                tqx -= 2.0 / 3.0 * (mpole[t110] * phi[t101] + mpole[t020] * phi[t011] + mpole[t011] * phi[t002] - mpole[t101] * phi[t110] - mpole[t011] * phi[t020] - mpole[t002] * phi[t011]);
//                tqy -= 2.0 / 3.0 * (mpole[t101] * phi[t200] + mpole[t011] * phi[t110] + mpole[t002] * phi[t101] - mpole[t200] * phi[t101] - mpole[t110] * phi[t011] - mpole[t101] * phi[t002]);
//                tqz -= 2.0 / 3.0 * (mpole[t200] * phi[t110] + mpole[t110] * phi[t020] + mpole[t101] * phi[t011] - mpole[t110] * phi[t200] - mpole[t020] * phi[t110] - mpole[t011] * phi[t101]);
//                sharedGrad[0].addAndGet(i, electric * dfx);
//                sharedGrad[1].addAndGet(i, electric * dfy);
//                sharedGrad[2].addAndGet(i, electric * dfz);
//                sharedTorque[0].addAndGet(i, electric * tqx);
//                sharedTorque[1].addAndGet(i, electric * tqy);
//                sharedTorque[2].addAndGet(i, electric * tqz);
//            }
//        }
//        erecip = 0.5 * electric * erecip;
//        return erecip;
//    }
    private double inducedDipoleSelfEnergy(boolean gradient) {
        double e = 0.0;
        final double term = -2.0 / 3.0 * electric * aewald * aewald * aewald / sqrtPi;
        final double ind[][] = inducedDipole[0];
        final double indp[][] = inducedDipolep[0];
        final double mpole[][] = globalMultipole[0];
        for (int i = 0; i < nAtoms; i++) {
            final double indi[] = ind[i];
            final double multipolei[] = mpole[i];
            final double dix = multipolei[t100];
            final double diy = multipolei[t010];
            final double diz = multipolei[t001];
            final double dii = indi[0] * dix + indi[1] * diy + indi[2] * diz;
            e += term * dii;
        }
        if (gradient) {
            final double fterm = -2.0 * term;
            for (int i = 0; i < nAtoms; i++) {
                final double indi[] = ind[i];
                final double indpi[] = indp[i];
                final double multipolei[] = mpole[i];
                final double dix = multipolei[t100];
                final double diy = multipolei[t010];
                final double diz = multipolei[t001];
                final double uix = 0.5 * (indi[0] + indpi[0]);
                final double uiy = 0.5 * (indi[1] + indpi[1]);
                final double uiz = 0.5 * (indi[2] + indpi[2]);
                sharedTorque[0].addAndGet(i, fterm * (diy * uiz - diz * uiy));
                sharedTorque[1].addAndGet(i, fterm * (diz * uix - dix * uiz));
                sharedTorque[2].addAndGet(i, fterm * (dix * uiy - diy * uix));
            }
        }
        return e;
    }

//    private double inducedDipoleReciprocalSpaceEnergy(boolean gradient) {
//        double e = 0.0;
//        if (gradient && polarization == Polarization.DIRECT) {
//            try {
//                reciprocalSpace.splineInducedDipoles(inducedDipole, inducedDipolep, null);
//                sectionTeam.execute(inducedDipoleFieldRegion);
//                reciprocalSpace.computeInducedPhi(cartesianDipolePhi, cartesianDipolepPhi);
//            } catch (Exception ex) {
//                String message = "Fatal exception computing the induced reciprocal space field.\n";
//                logger.log(Level.SEVERE, message, ex);
//            }
//            for (int i = 0; i < nAtoms; i++) {
//                final double fieldi[] = field1[i];
//                final double phii[] = cartesianDipolePhi[i];
//                fieldi[0] -= phii[t100];
//                fieldi[1] -= phii[t010];
//                fieldi[2] -= phii[t001];
//                final double fieldpi[] = field2[i];
//                final double phipi[] = cartesianDipolepPhi[i];
//                fieldpi[0] -= phipi[t100];
//                fieldpi[1] -= phipi[t010];
//                fieldpi[2] -= phipi[t001];
//            }
//        } else {
//            reciprocalSpace.cartToFracInducedDipoles(inducedDipole, inducedDipolep);
//        }
//        final double nfftX = reciprocalSpace.getXDim();
//        final double nfftY = reciprocalSpace.getYDim();
//        final double nfftZ = reciprocalSpace.getZDim();
//        final double mpole[][] = globalMultipole[0];
//        final double fractionalMultipolePhi[][] = reciprocalSpace.getFracMultipolePhi();
//        final double fractionalInducedDipolePhi[][] = reciprocalSpace.getFracInducedDipolePhi();
//        final double fractionalInducedDipolepPhi[][] = reciprocalSpace.getFracInducedDipoleCRPhi();
//        final double fmpole[][] = reciprocalSpace.getFracMultipoles();
//        final double find[][] = reciprocalSpace.getFracInducedDipoles();
//        final double finp[][] = reciprocalSpace.getFracInducedDipolesCR();
//        for (int i = 0; i < nAtoms; i++) {
//            final double fPhi[] = fractionalMultipolePhi[i];
//            final double findi[] = find[i];
//            final double indx = findi[0];
//            final double indy = findi[1];
//            final double indz = findi[2];
//            e += indx * fPhi[t100] + indy * fPhi[t010] + indz * fPhi[t001];
//            if (gradient) {
//                final double iPhi[] = cartesianDipolePhi[i];
//                final double ipPhi[] = cartesianDipolepPhi[i];
//                final double fiPhi[] = fractionalInducedDipolePhi[i];
//                final double fipPhi[] = fractionalInducedDipolepPhi[i];
//                final double mpolei[] = mpole[i];
//                final double fmpolei[] = fmpole[i];
//                final double finpi[] = finp[i];
//                final double inpx = finpi[0];
//                final double inpy = finpi[1];
//                final double inpz = finpi[2];
//                final double insx = indx + inpx;
//                final double insy = indy + inpy;
//                final double insz = indz + inpz;
//                for (int t = 0; t < tensorCount; t++) {
//                    sPhi[t] = 0.5 * (iPhi[t] + ipPhi[t]);
//                    sfPhi[t] = fiPhi[t] + fipPhi[t];
//                }
//                double gx = insx * fPhi[t200] + insy * fPhi[t110] + insz * fPhi[t101];
//                double gy = insx * fPhi[t110] + insy * fPhi[t020] + insz * fPhi[t011];
//                double gz = insx * fPhi[t101] + insy * fPhi[t011] + insz * fPhi[t002];
//                if (polarization == Polarization.MUTUAL) {
//                    gx += indx * fipPhi[t200] + inpx * fiPhi[t200] + indy * fipPhi[t110] + inpy * fiPhi[t110] + indz * fipPhi[t101] + inpz * fiPhi[t101];
//                    gy += indx * fipPhi[t110] + inpx * fiPhi[t110] + indy * fipPhi[t020] + inpy * fiPhi[t020] + indz * fipPhi[t011] + inpz * fiPhi[t011];
//                    gz += indx * fipPhi[t101] + inpx * fiPhi[t101] + indy * fipPhi[t011] + inpy * fiPhi[t011] + indz * fipPhi[t002] + inpz * fiPhi[t002];
//                }
//                gx += fmpolei[t000] * sfPhi[t100] + fmpolei[t100] * sfPhi[t200] + fmpolei[t010] * sfPhi[t110] + fmpolei[t001] * sfPhi[t101] + fmpolei[t200] * sfPhi[t300] + fmpolei[t020] * sfPhi[t120] + fmpolei[t002] * sfPhi[t102] + fmpolei[t110] * sfPhi[t210] + fmpolei[t101] * sfPhi[t201] + fmpolei[t011] * sfPhi[t111];
//                gy += fmpolei[t000] * sfPhi[t010] + fmpolei[t100] * sfPhi[t110] + fmpolei[t010] * sfPhi[t020] + fmpolei[t001] * sfPhi[t011] + fmpolei[t200] * sfPhi[t210] + fmpolei[t020] * sfPhi[t030] + fmpolei[t002] * sfPhi[t012] + fmpolei[t110] * sfPhi[t120] + fmpolei[t101] * sfPhi[t111] + fmpolei[t011] * sfPhi[t021];
//                gz += fmpolei[t000] * sfPhi[t001] + fmpolei[t100] * sfPhi[t101] + fmpolei[t010] * sfPhi[t011] + fmpolei[t001] * sfPhi[t002] + fmpolei[t200] * sfPhi[t201] + fmpolei[t020] * sfPhi[t021] + fmpolei[t002] * sfPhi[t003] + fmpolei[t110] * sfPhi[t111] + fmpolei[t101] * sfPhi[t102] + fmpolei[t011] * sfPhi[t012];
//                gx *= nfftX;
//                gy *= nfftY;
//                gz *= nfftZ;
//                double recip[][] = crystal.getUnitCell().A;
//                final double dfx = recip[0][0] * gx + recip[0][1] * gy + recip[0][2] * gz;
//                final double dfy = recip[1][0] * gx + recip[1][1] * gy + recip[1][2] * gz;
//                final double dfz = recip[2][0] * gx + recip[2][1] * gy + recip[2][2] * gz;
//                // Compute dipole torques
//                double tqx = -mpolei[t010] * sPhi[t001] + mpolei[t001] * sPhi[t010];
//                double tqy = -mpolei[t001] * sPhi[t100] + mpolei[t100] * sPhi[t001];
//                double tqz = -mpolei[t100] * sPhi[t010] + mpolei[t010] * sPhi[t100];
//                // Compute quadrupole torques
//                tqx -= 2.0 / 3.0 * (mpolei[t110] * sPhi[t101] + mpolei[t020] * sPhi[t011] + mpolei[t011] * sPhi[t002] - mpolei[t101] * sPhi[t110] - mpolei[t011] * sPhi[t020] - mpolei[t002] * sPhi[t011]);
//                tqy -= 2.0 / 3.0 * (mpolei[t101] * sPhi[t200] + mpolei[t011] * sPhi[t110] + mpolei[t002] * sPhi[t101] - mpolei[t200] * sPhi[t101] - mpolei[t110] * sPhi[t011] - mpolei[t101] * sPhi[t002]);
//                tqz -= 2.0 / 3.0 * (mpolei[t200] * sPhi[t110] + mpolei[t110] * sPhi[t020] + mpolei[t101] * sPhi[t011] - mpolei[t110] * sPhi[t200] - mpolei[t020] * sPhi[t110] - mpolei[t011] * sPhi[t101]);
//                sharedGrad[0].addAndGet(i, 0.5 * electric * dfx);
//                sharedGrad[1].addAndGet(i, 0.5 * electric * dfy);
//                sharedGrad[2].addAndGet(i, 0.5 * electric * dfz);
//                sharedTorque[0].addAndGet(i, electric * tqx);
//                sharedTorque[1].addAndGet(i, electric * tqy);
//                sharedTorque[2].addAndGet(i, electric * tqz);
//            }
//        }
//        e *= 0.5 * electric;
//        return e;
//    }
    private class PermanentRealSpaceFieldRegion extends ParallelRegion {

        private final PermanentRealSpaceFieldLoop permanentRealSpaceFieldLoop[];
        private final SharedDoubleArray sharedField[];
        private final SharedDoubleArray sharedFieldp[];

        public PermanentRealSpaceFieldRegion(int nt) {
            super();
            sharedField = new SharedDoubleArray[3];
            sharedField[0] = new SharedDoubleArray(nAtoms);
            sharedField[1] = new SharedDoubleArray(nAtoms);
            sharedField[2] = new SharedDoubleArray(nAtoms);
            sharedFieldp = new SharedDoubleArray[3];
            sharedFieldp[0] = new SharedDoubleArray(nAtoms);
            sharedFieldp[1] = new SharedDoubleArray(nAtoms);
            sharedFieldp[2] = new SharedDoubleArray(nAtoms);
            permanentRealSpaceFieldLoop = new PermanentRealSpaceFieldLoop[nt];
            for (int i = 0; i < nt; i++) {
                permanentRealSpaceFieldLoop[i] = new PermanentRealSpaceFieldLoop();
            }
        }

        public void setField(double fld[][], double fldp[][]) {
            for (int i = 0; i < nAtoms; i++) {
                fld[i][0] = sharedField[0].get(i);
                fld[i][1] = sharedField[1].get(i);
                fld[i][2] = sharedField[2].get(i);
                fldp[i][0] = sharedFieldp[0].get(i);
                fldp[i][1] = sharedFieldp[1].get(i);
                fldp[i][2] = sharedFieldp[2].get(i);
            }
        }

        @Override
        public void start() {
            for (int i = 0; i < nAtoms; i++) {
                sharedField[0].set(i, 0.0);
                sharedField[1].set(i, 0.0);
                sharedField[2].set(i, 0.0);
                sharedFieldp[0].set(i, 0.0);
                sharedFieldp[1].set(i, 0.0);
                sharedFieldp[2].set(i, 0.0);
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1,
                        permanentRealSpaceFieldLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing the real space field in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
                System.exit(-1);
            }
        }

        @Override
        public void finish() {
            // logger.info(format("\nPermanent Real Space Field: %10.3f seconds\n",
            // (System.nanoTime() - time) * 0.000000001));
        }

        private class PermanentRealSpaceFieldLoop extends IntegerForLoop {

            private final double mask_local[];
            private final double maskp_local[];
            private final double fx_local[];
            private final double fy_local[];
            private final double fz_local[];
            private final double fxp_local[];
            private final double fyp_local[];
            private final double fzp_local[];
            private final double dx_local[];

            public PermanentRealSpaceFieldLoop() {
                super();
                mask_local = new double[nAtoms];
                maskp_local = new double[nAtoms];
                fx_local = new double[nAtoms];
                fy_local = new double[nAtoms];
                fz_local = new double[nAtoms];
                fxp_local = new double[nAtoms];
                fyp_local = new double[nAtoms];
                fzp_local = new double[nAtoms];
                dx_local = new double[3];
                for (int i = 0; i < nAtoms; i++) {
                    mask_local[i] = 1.0;
                    maskp_local[i] = 1.0;
                }
            }

            @Override
            public IntegerSchedule schedule() {
                return pairWiseSchedule;
            }

            @Override
            public void start() {
                for (int i = 0; i < nAtoms; i++) {
                    fx_local[i] = 0.0;
                    fy_local[i] = 0.0;
                    fz_local[i] = 0.0;
                    fxp_local[i] = 0.0;
                    fyp_local[i] = 0.0;
                    fzp_local[i] = 0.0;
                }
            }

            @Override
            public void run(int lb, int ub) {
                int lists[][] = neighborLists[0];
                int ewalds[][] = ewaldLists[0];
                int counts[] = ewaldCounts[0];
                final double x[] = coordinates[0][0];
                final double y[] = coordinates[0][1];
                final double z[] = coordinates[0][2];
                final double mpole[][] = globalMultipole[0];
                /**
                 * Loop over atom chunk.
                 */
                for (int i = lb; i <= ub; i++) {
                    final double pdi = pdamp[i];
                    final double pti = thole[i];
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double globalMultipolei[] = mpole[i];
                    final double ci = globalMultipolei[0];
                    final double dix = globalMultipolei[t100];
                    final double diy = globalMultipolei[t010];
                    final double diz = globalMultipolei[t001];
                    final double qixx = globalMultipolei[t200] / 3.0;
                    final double qiyy = globalMultipolei[t020] / 3.0;
                    final double qizz = globalMultipolei[t002] / 3.0;
                    final double qixy = globalMultipolei[t110] / 3.0;
                    final double qixz = globalMultipolei[t101] / 3.0;
                    final double qiyz = globalMultipolei[t011] / 3.0;
                    /**
                     * Apply energy masking rules.
                     */
                    Atom ai = atoms[i];
                    for (Torsion torsion : ai.getTorsions()) {
                        Atom ak = torsion.get1_4(ai);
                        if (ak != null) {
                            int index = ak.xyzIndex - 1;
                            for (int k : ip11[i]) {
                                if (k == index) {
                                    maskp_local[index] = 0.5;
                                }
                            }
                        }
                    }
                    for (Angle angle : ai.getAngles()) {
                        Atom ak = angle.get1_3(ai);
                        if (ak != null) {
                            int index = ak.xyzIndex - 1;
                            maskp_local[index] = p13scale;
                        }
                    }
                    for (Bond bond : ai.getFFXBonds()) {
                        int index = bond.get1_2(ai).xyzIndex - 1;
                        maskp_local[index] = p12scale;
                    }
                    /**
                     * Apply group based polarization masking rule.
                     */
                    for (int index : ip11[i]) {
                        mask_local[index] = d11scale;
                    }
                    /**
                     * Loop over the neighbor list.
                     */
                    final int list[] = lists[i];
                    int npair = list.length;
                    counts[i] = 0;
                    final int ewald[] = ewalds[i];
                    for (int j = 0; j < npair; j++) {
                        int k = list[j];
                        final double xk = x[k];
                        final double yk = y[k];
                        final double zk = z[k];
                        dx_local[0] = xk - xi;
                        dx_local[1] = yk - yi;
                        dx_local[2] = zk - zi;
                        final double r2 = crystal.image(dx_local);
                        if (r2 <= off2) {
                            ewald[counts[i]++] = k;
                            final double xr = dx_local[0];
                            final double yr = dx_local[1];
                            final double zr = dx_local[2];
                            final double pdk = pdamp[k];
                            final double ptk = thole[k];
                            final double globalMultipolek[] = mpole[k];
                            final double ck = globalMultipolek[t000];
                            final double dkx = globalMultipolek[t100];
                            final double dky = globalMultipolek[t010];
                            final double dkz = globalMultipolek[t001];
                            final double qkxx = globalMultipolek[t200] / 3.0;
                            final double qkyy = globalMultipolek[t020] / 3.0;
                            final double qkzz = globalMultipolek[t002] / 3.0;
                            final double qkxy = globalMultipolek[t110] / 3.0;
                            final double qkxz = globalMultipolek[t101] / 3.0;
                            final double qkyz = globalMultipolek[t011] / 3.0;
                            final double r = sqrt(r2);
                            /**
                             * Calculate the error function damping terms.
                             */
                            final double ralpha = aewald * r;
                            final double bn0 = erfc(ralpha) / r;
                            double alsq2n = piEwald;
                            final double exp2a = exp(-ralpha * ralpha);
                            alsq2n = alsq2 * alsq2n;
                            final double bn1 = (bn0 + alsq2n * exp2a) / r2;
                            alsq2n = alsq2 * alsq2n;
                            final double bn2 = (3.0 * bn1 + alsq2n * exp2a) / r2;
                            alsq2n = alsq2 * alsq2n;
                            final double bn3 = (5.0 * bn2 + alsq2n * exp2a) / r2;
                            /**
                             * Compute the error function scaled and unscaled
                             * terms.
                             */
                            double scale3 = 1.0;
                            double scale5 = 1.0;
                            double scale7 = 1.0;
                            double damp = pdi * pdk;
                            double expdamp = 0.0;
                            if (damp != 0.0) {
                                final double pgamma = min(pti, ptk);
                                final double rdamp = r / damp;
                                damp = -pgamma * rdamp * rdamp * rdamp;
                                if (damp > -50.0) {
                                    expdamp = exp(damp);
                                    scale3 = 1.0 - expdamp;
                                    scale5 = 1.0 - expdamp * (1.0 - damp);
                                    scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                                }
                            }
                            final double scale = mask_local[k];
                            final double scalep = maskp_local[k];
                            final double dsc3 = scale3 * scale;
                            final double dsc5 = scale5 * scale;
                            final double dsc7 = scale7 * scale;
                            final double psc3 = scale3 * scalep;
                            final double psc5 = scale5 * scalep;
                            final double psc7 = scale7 * scalep;
                            final double rr3 = 1.0 / (r * r2);
                            final double rr5 = 3.0 * rr3 / r2;
                            final double rr7 = 5.0 * rr5 / r2;
                            final double drr3 = (1.0 - dsc3) * rr3;
                            final double drr5 = (1.0 - dsc5) * rr5;
                            final double drr7 = (1.0 - dsc7) * rr7;
                            final double prr3 = (1.0 - psc3) * rr3;
                            final double prr5 = (1.0 - psc5) * rr5;
                            final double prr7 = (1.0 - psc7) * rr7;
                            final double dir = dix * xr + diy * yr + diz * zr;
                            final double qix = 2.0 * (qixx * xr + qixy * yr + qixz * zr);
                            final double qiy = 2.0 * (qixy * xr + qiyy * yr + qiyz * zr);
                            final double qiz = 2.0 * (qixz * xr + qiyz * yr + qizz * zr);
                            final double qir = (qix * xr + qiy * yr + qiz * zr) / 2.0;
                            final double dkr = dkx * xr + dky * yr + dkz * zr;
                            final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                            final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                            final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                            final double qkr = (qkx * xr + qky * yr + qkz * zr) / 2.0;
                            final double fimx = -xr * (bn1 * ck - bn2 * dkr + bn3 * qkr) - bn1 * dkx + bn2 * qkx;
                            final double fimy = -yr * (bn1 * ck - bn2 * dkr + bn3 * qkr) - bn1 * dky + bn2 * qky;
                            final double fimz = -zr * (bn1 * ck - bn2 * dkr + bn3 * qkr) - bn1 * dkz + bn2 * qkz;
                            final double fkmx = xr * (bn1 * ci + bn2 * dir + bn3 * qir) - bn1 * dix - bn2 * qix;
                            final double fkmy = yr * (bn1 * ci + bn2 * dir + bn3 * qir) - bn1 * diy - bn2 * qiy;
                            final double fkmz = zr * (bn1 * ci + bn2 * dir + bn3 * qir) - bn1 * diz - bn2 * qiz;
                            final double fidx = -xr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * dkx + drr5 * qkx;
                            final double fidy = -yr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * dky + drr5 * qky;
                            final double fidz = -zr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * dkz + drr5 * qkz;
                            final double fkdx = xr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * dix - drr5 * qix;
                            final double fkdy = yr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * diy - drr5 * qiy;
                            final double fkdz = zr * (drr3 * ci + drr5 * dir + drr7 * qir) - drr3 * diz - drr5 * qiz;
                            final double fipx = -xr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * dkx + prr5 * qkx;
                            final double fipy = -yr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * dky + prr5 * qky;
                            final double fipz = -zr * (prr3 * ck - prr5 * dkr + prr7 * qkr) - prr3 * dkz + prr5 * qkz;
                            final double fkpx = xr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * dix - prr5 * qix;
                            final double fkpy = yr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * diy - prr5 * qiy;
                            final double fkpz = zr * (prr3 * ci + prr5 * dir + prr7 * qir) - prr3 * diz - prr5 * qiz;
                            fx_local[i] += fimx - fidx;
                            fy_local[i] += fimy - fidy;
                            fz_local[i] += fimz - fidz;
                            fx_local[k] += fkmx - fkdx;
                            fy_local[k] += fkmy - fkdy;
                            fz_local[k] += fkmz - fkdz;
                            fxp_local[i] += fimx - fipx;
                            fyp_local[i] += fimy - fipy;
                            fzp_local[i] += fimz - fipz;
                            fxp_local[k] += fkmx - fkpx;
                            fyp_local[k] += fkmy - fkpy;
                            fzp_local[k] += fkmz - fkpz;
                        }
                    }
                    for (Torsion torsion : ai.getTorsions()) {
                        Atom ak = torsion.get1_4(ai);
                        if (ak != null) {
                            int index = ak.xyzIndex - 1;
                            maskp_local[index] = 1.0;
                        }
                    }
                    for (Angle angle : ai.getAngles()) {
                        Atom ak = angle.get1_3(ai);
                        if (ak != null) {
                            int index = ak.xyzIndex - 1;
                            maskp_local[index] = 1.0;
                        }
                    }
                    for (Bond bond : ai.getFFXBonds()) {
                        int index = bond.get1_2(ai).xyzIndex - 1;
                        maskp_local[index] = 1.0;
                    }
                    for (int index : ip11[i]) {
                        mask_local[index] = 1.0;
                    }
                }
                /**
                 * Loop over symmetry mates.
                 */
                for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                    lists = neighborLists[iSymm];
                    ewalds = ewaldLists[iSymm];
                    counts = ewaldCounts[iSymm];
                    double xs[] = coordinates[iSymm][0];
                    double ys[] = coordinates[iSymm][1];
                    double zs[] = coordinates[iSymm][2];
                    double mpoles[][] = globalMultipole[iSymm];
                    /**
                     * Loop over atoms in a chunk of the asymmetric unit.
                     */
                    for (int i = lb; i <= ub; i++) {
                        final double pdi = pdamp[i];
                        final double pti = thole[i];
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        /**
                         * Loop over the neighbor list.
                         */
                        final int list[] = lists[i];
                        final int npair = list.length;
                        counts[i] = 0;
                        final int ewald[] = ewalds[i];
                        for (int j = 0; j < npair; j++) {
                            int k = list[j];
                            final double xk = xs[k];
                            final double yk = ys[k];
                            final double zk = zs[k];
                            dx_local[0] = xk - xi;
                            dx_local[1] = yk - yi;
                            dx_local[2] = zk - zi;
                            final double r2 = crystal.image(dx_local);
                            if (r2 <= off2) {
                                ewald[counts[i]++] = k;
                                final double xr = dx_local[0];
                                final double yr = dx_local[1];
                                final double zr = dx_local[2];
                                final double pdk = pdamp[k];
                                final double ptk = thole[k];
                                final double multipolek[] = mpoles[k];
                                final double ck = multipolek[t000];
                                final double dkx = multipolek[t100];
                                final double dky = multipolek[t010];
                                final double dkz = multipolek[t001];
                                final double qkxx = multipolek[t200] / 3.0;
                                final double qkyy = multipolek[t020] / 3.0;
                                final double qkzz = multipolek[t002] / 3.0;
                                final double qkxy = multipolek[t110] / 3.0;
                                final double qkxz = multipolek[t101] / 3.0;
                                final double qkyz = multipolek[t011] / 3.0;
                                final double r = sqrt(r2);
                                /**
                                 * Calculate the error function damping terms.
                                 */
                                final double ralpha = aewald * r;
                                final double bn0 = erfc(ralpha) / r;
                                double alsq2n = piEwald;
                                final double exp2a = exp(-ralpha * ralpha);
                                alsq2n = alsq2 * alsq2n;
                                final double bn1 = (bn0 + alsq2n * exp2a) / r2;
                                alsq2n = alsq2 * alsq2n;
                                final double bn2 = (3.0 * bn1 + alsq2n * exp2a) / r2;
                                alsq2n = alsq2 * alsq2n;
                                final double bn3 = (5.0 * bn2 + alsq2n * exp2a) / r2;
                                /**
                                 * Compute the error function scaled and
                                 * unscaled terms.
                                 */
                                double scale3 = 1.0;
                                double scale5 = 1.0;
                                double scale7 = 1.0;
                                double damp = pdi * pdk;
                                double expdamp = 0.0;
                                if (damp != 0.0) {
                                    final double pgamma = min(pti, ptk);
                                    final double rdamp = r / damp;
                                    damp = -pgamma * rdamp * rdamp * rdamp;
                                    if (damp > -50.0) {
                                        expdamp = exp(damp);
                                        scale3 = 1.0 - expdamp;
                                        scale5 = 1.0 - expdamp * (1.0 - damp);
                                        scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                                    }
                                }
                                final double dsc3 = scale3;
                                final double dsc5 = scale5;
                                final double dsc7 = scale7;
                                final double rr3 = 1.0 / (r * r2);
                                final double rr5 = 3.0 * rr3 / r2;
                                final double rr7 = 5.0 * rr5 / r2;
                                final double drr3 = (1.0 - dsc3) * rr3;
                                final double drr5 = (1.0 - dsc5) * rr5;
                                final double drr7 = (1.0 - dsc7) * rr7;
                                final double dkr = dkx * xr + dky * yr + dkz * zr;
                                final double qkx = 2.0 * (qkxx * xr + qkxy * yr + qkxz * zr);
                                final double qky = 2.0 * (qkxy * xr + qkyy * yr + qkyz * zr);
                                final double qkz = 2.0 * (qkxz * xr + qkyz * yr + qkzz * zr);
                                final double qkr = (qkx * xr + qky * yr + qkz * zr) / 2.0;
                                final double fimx = -xr * (bn1 * ck - bn2 * dkr + bn3 * qkr) - bn1 * dkx + bn2 * qkx;
                                final double fimy = -yr * (bn1 * ck - bn2 * dkr + bn3 * qkr) - bn1 * dky + bn2 * qky;
                                final double fimz = -zr * (bn1 * ck - bn2 * dkr + bn3 * qkr) - bn1 * dkz + bn2 * qkz;
                                final double fidx = -xr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * dkx + drr5 * qkx;
                                final double fidy = -yr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * dky + drr5 * qky;
                                final double fidz = -zr * (drr3 * ck - drr5 * dkr + drr7 * qkr) - drr3 * dkz + drr5 * qkz;
                                fx_local[i] += fimx - fidx;
                                fy_local[i] += fimy - fidy;
                                fz_local[i] += fimz - fidz;
                                fxp_local[i] += fimx - fidx;
                                fyp_local[i] += fimy - fidy;
                                fzp_local[i] += fimz - fidz;
                            }
                        }
                    }
                }
            }

            @Override
            public void finish() {
                sharedField[0].reduce(fx_local, DoubleOp.SUM);
                sharedField[1].reduce(fy_local, DoubleOp.SUM);
                sharedField[2].reduce(fz_local, DoubleOp.SUM);
                sharedFieldp[0].reduce(fxp_local, DoubleOp.SUM);
                sharedFieldp[1].reduce(fyp_local, DoubleOp.SUM);
                sharedFieldp[2].reduce(fzp_local, DoubleOp.SUM);
            }
        }
    }

    private class InducedDipoleRealSpaceFieldRegion extends ParallelRegion {

        private final PolarizationRealSpaceFieldLoop polarizationRealSpaceFieldLoop[];
        private final SharedDoubleArray sharedField[];
        private final SharedDoubleArray sharedFieldp[];

        public InducedDipoleRealSpaceFieldRegion(int nt) {
            super();
            sharedField = new SharedDoubleArray[3];
            sharedField[0] = new SharedDoubleArray(nAtoms);
            sharedField[1] = new SharedDoubleArray(nAtoms);
            sharedField[2] = new SharedDoubleArray(nAtoms);
            sharedFieldp = new SharedDoubleArray[3];
            sharedFieldp[0] = new SharedDoubleArray(nAtoms);
            sharedFieldp[1] = new SharedDoubleArray(nAtoms);
            sharedFieldp[2] = new SharedDoubleArray(nAtoms);
            polarizationRealSpaceFieldLoop = new PolarizationRealSpaceFieldLoop[nt];
            for (int i = 0; i < nt; i++) {
                polarizationRealSpaceFieldLoop[i] = new PolarizationRealSpaceFieldLoop();
            }
        }

        public void setField(double fld[][], double fldp[][]) {
            for (int i = 0; i < nAtoms; i++) {
                fld[i][0] = sharedField[0].get(i);
                fld[i][1] = sharedField[1].get(i);
                fld[i][2] = sharedField[2].get(i);
                fldp[i][0] = sharedFieldp[0].get(i);
                fldp[i][1] = sharedFieldp[1].get(i);
                fldp[i][2] = sharedFieldp[2].get(i);
            }
        }

        /**
         * Initialize the shared field arrays to 0.
         */
        @Override
        public void start() {
            for (int i = 0; i < nAtoms; i++) {
                sharedField[0].set(i, 0.0);
                sharedField[1].set(i, 0.0);
                sharedField[2].set(i, 0.0);
                sharedFieldp[0].set(i, 0.0);
                sharedFieldp[1].set(i, 0.0);
                sharedFieldp[2].set(i, 0.0);
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1,
                        polarizationRealSpaceFieldLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing the induced real space field in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class PolarizationRealSpaceFieldLoop extends IntegerForLoop {

            private int list[], lists[][], counts[];
            private int npair, i, j, k, iSymm;
            private double fx, fy, fz;
            private double px, py, pz;
            private double xi, yi, zi;
            private double pdi, pdk, pti, ptk;
            private double dipolei[], dipolepi[];
            private double uix, uiy, uiz;
            private double pix, piy, piz;
            private double xr, yr, zr;
            private double dipolek[], dipolepk[];
            private double ukx, uky, ukz;
            private double pkx, pky, pkz;
            private double bn0, bn1, bn2;
            private double scale3, scale5, damp, expdamp, pgamma, rdamp;
            private double r, ralpha, alsq2n, exp2a, r2, rr3, rr5;
            private double uir, ukr, pir, pkr;
            private double bn2ukr, bn2uir, bn2pkr, bn2pir;
            private double rr5ukr, rr5uir, rr5pkr, rr5pir;
            private double fimx, fimy, fimz;
            private double fkmx, fkmy, fkmz;
            private double fidx, fidy, fidz;
            private double fkdx, fkdy, fkdz;
            private double pimx, pimy, pimz;
            private double pkmx, pkmy, pkmz;
            private double pidx, pidy, pidz;
            private double pkdx, pkdy, pkdz;
            private double xs[], ys[], zs[];
            private double inds[][], indps[][];
            private final double fx_local[];
            private final double fy_local[];
            private final double fz_local[];
            private final double fxp_local[];
            private final double fyp_local[];
            private final double fzp_local[];
            private final double dx_local[];
            private final double x[] = coordinates[0][0];
            private final double y[] = coordinates[0][1];
            private final double z[] = coordinates[0][2];
            private final double ind[][] = inducedDipole[0];
            private final double inp[][] = inducedDipolep[0];

            public PolarizationRealSpaceFieldLoop() {
                super();
                fx_local = new double[nAtoms];
                fy_local = new double[nAtoms];
                fz_local = new double[nAtoms];
                fxp_local = new double[nAtoms];
                fyp_local = new double[nAtoms];
                fzp_local = new double[nAtoms];
                dx_local = new double[3];
            }

            @Override
            public IntegerSchedule schedule() {
                return pairWiseSchedule;
            }

            @Override
            public void start() {
                for (i = 0; i < nAtoms; i++) {
                    fx_local[i] = 0.0;
                    fy_local[i] = 0.0;
                    fz_local[i] = 0.0;
                    fxp_local[i] = 0.0;
                    fyp_local[i] = 0.0;
                    fzp_local[i] = 0.0;
                }
            }

            @Override
            public void run(int lb, int ub) {
                /**
                 * Loop over a chunk of atoms.
                 */
                lists = ewaldLists[0];
                counts = ewaldCounts[0];
                for (i = lb; i <= ub; i++) {
                    fx = 0.0;
                    fy = 0.0;
                    fz = 0.0;
                    px = 0.0;
                    py = 0.0;
                    pz = 0.0;
                    xi = x[i];
                    yi = y[i];
                    zi = z[i];
                    dipolei = ind[i];
                    uix = dipolei[0];
                    uiy = dipolei[1];
                    uiz = dipolei[2];
                    dipolepi = inp[i];
                    pix = dipolepi[0];
                    piy = dipolepi[1];
                    piz = dipolepi[2];
                    pdi = pdamp[i];
                    pti = thole[i];
                    /**
                     * Loop over the neighbor list.
                     */
                    list = lists[i];
                    npair = counts[i];
                    for (j = 0; j < npair; j++) {
                        k = list[j];
                        pdk = pdamp[k];
                        ptk = thole[k];
                        dx_local[0] = x[k] - xi;
                        dx_local[1] = y[k] - yi;
                        dx_local[2] = z[k] - zi;
                        r2 = crystal.image(dx_local);
                        xr = dx_local[0];
                        yr = dx_local[1];
                        zr = dx_local[2];
                        dipolek = ind[k];
                        ukx = dipolek[0];
                        uky = dipolek[1];
                        ukz = dipolek[2];
                        dipolepk = inp[k];
                        pkx = dipolepk[0];
                        pky = dipolepk[1];
                        pkz = dipolepk[2];
                        uir = uix * xr + uiy * yr + uiz * zr;
                        ukr = ukx * xr + uky * yr + ukz * zr;
                        pir = pix * xr + piy * yr + piz * zr;
                        pkr = pkx * xr + pky * yr + pkz * zr;
                        /**
                         * Calculate the error function damping terms.
                         */
                        r = sqrt(r2);
                        ralpha = aewald * r;
                        bn0 = erfc(ralpha) / r;
                        alsq2n = piEwald;
                        exp2a = exp(-ralpha * ralpha);
                        alsq2n = alsq2 * alsq2n;
                        bn1 = (bn0 + alsq2n * exp2a) / r2;
                        alsq2n = alsq2 * alsq2n;
                        bn2 = (3.0 * bn1 + alsq2n * exp2a) / r2;
                        scale3 = 1.0;
                        scale5 = 1.0;
                        damp = pdi * pdk;
                        expdamp = 0.0;
                        if (damp != 0.0) {
                            pgamma = min(pti, ptk);
                            rdamp = r / damp;
                            damp = -pgamma * rdamp * rdamp * rdamp;
                            if (damp > -50.0) {
                                expdamp = exp(damp);
                                scale3 = 1.0 - expdamp;
                                scale5 = 1.0 - expdamp * (1.0 - damp);
                            }
                        }
                        rr3 = 1.0 / (r * r2);
                        rr5 = 3.0 * rr3 / r2;
                        rr3 *= (1.0 - scale3);
                        rr5 *= (1.0 - scale5);
                        bn2ukr = bn2 * ukr;
                        bn2uir = bn2 * uir;
                        bn2pkr = bn2 * pkr;
                        bn2pir = bn2 * pir;
                        rr5ukr = rr5 * ukr;
                        rr5uir = rr5 * uir;
                        rr5pkr = rr5 * pkr;
                        rr5pir = rr5 * pir;
                        fimx = -bn1 * ukx + bn2ukr * xr;
                        fimy = -bn1 * uky + bn2ukr * yr;
                        fimz = -bn1 * ukz + bn2ukr * zr;
                        fkmx = -bn1 * uix + bn2uir * xr;
                        fkmy = -bn1 * uiy + bn2uir * yr;
                        fkmz = -bn1 * uiz + bn2uir * zr;
                        fidx = -rr3 * ukx + rr5ukr * xr;
                        fidy = -rr3 * uky + rr5ukr * yr;
                        fidz = -rr3 * ukz + rr5ukr * zr;
                        fkdx = -rr3 * uix + rr5uir * xr;
                        fkdy = -rr3 * uiy + rr5uir * yr;
                        fkdz = -rr3 * uiz + rr5uir * zr;
                        pimx = -bn1 * pkx + bn2pkr * xr;
                        pimy = -bn1 * pky + bn2pkr * yr;
                        pimz = -bn1 * pkz + bn2pkr * zr;
                        pkmx = -bn1 * pix + bn2pir * xr;
                        pkmy = -bn1 * piy + bn2pir * yr;
                        pkmz = -bn1 * piz + bn2pir * zr;
                        pidx = -rr3 * pkx + rr5pkr * xr;
                        pidy = -rr3 * pky + rr5pkr * yr;
                        pidz = -rr3 * pkz + rr5pkr * zr;
                        pkdx = -rr3 * pix + rr5pir * xr;
                        pkdy = -rr3 * piy + rr5pir * yr;
                        pkdz = -rr3 * piz + rr5pir * zr;
                        fx += fimx - fidx;
                        fy += fimy - fidy;
                        fz += fimz - fidz;
                        px += pimx - pidx;
                        py += pimy - pidy;
                        pz += pimz - pidz;
                        fx_local[k] += fkmx - fkdx;
                        fy_local[k] += fkmy - fkdy;
                        fz_local[k] += fkmz - fkdz;
                        fxp_local[k] += pkmx - pkdx;
                        fyp_local[k] += pkmy - pkdy;
                        fzp_local[k] += pkmz - pkdz;
                    }
                    fx_local[i] += fx;
                    fy_local[i] += fy;
                    fz_local[i] += fz;
                    fxp_local[i] += px;
                    fyp_local[i] += py;
                    fzp_local[i] += pz;
                }
                /**
                 * Loop over symmetry mates.
                 */
                for (iSymm = 1; iSymm < nSymm; iSymm++) {
                    lists = ewaldLists[iSymm];
                    counts = ewaldCounts[iSymm];
                    xs = coordinates[iSymm][0];
                    ys = coordinates[iSymm][1];
                    zs = coordinates[iSymm][2];
                    inds = inducedDipole[iSymm];
                    indps = inducedDipolep[iSymm];
                    /**
                     * Loop over a chunk of atoms.
                     */
                    for (i = lb; i <= ub; i++) {
                        fx = 0.0;
                        fy = 0.0;
                        fz = 0.0;
                        px = 0.0;
                        py = 0.0;
                        pz = 0.0;
                        xi = x[i];
                        yi = y[i];
                        zi = z[i];
                        pdi = pdamp[i];
                        pti = thole[i];
                        /**
                         * Loop over the neighbor list.
                         */
                        list = lists[i];
                        npair = counts[i];
                        for (j = 0; j < npair; j++) {
                            k = list[j];
                            pdk = pdamp[k];
                            ptk = thole[k];
                            dx_local[0] = xs[k] - xi;
                            dx_local[1] = ys[k] - yi;
                            dx_local[2] = zs[k] - zi;
                            r2 = crystal.image(dx_local);
                            xr = dx_local[0];
                            yr = dx_local[1];
                            zr = dx_local[2];
                            dipolek = inds[k];
                            dipolepk = indps[k];
                            ukx = dipolek[0];
                            uky = dipolek[1];
                            ukz = dipolek[2];
                            pkx = dipolepk[0];
                            pky = dipolepk[1];
                            pkz = dipolepk[2];
                            ukr = ukx * xr + uky * yr + ukz * zr;
                            pkr = pkx * xr + pky * yr + pkz * zr;
                            /**
                             * Calculate the error function damping terms.
                             */
                            r = sqrt(r2);
                            ralpha = aewald * r;
                            bn0 = erfc(ralpha) / r;
                            alsq2n = piEwald;
                            exp2a = exp(-ralpha * ralpha);
                            alsq2n = alsq2 * alsq2n;
                            bn1 = (bn0 + alsq2n * exp2a) / r2;
                            alsq2n = alsq2 * alsq2n;
                            bn2 = (3.0 * bn1 + alsq2n * exp2a) / r2;
                            scale3 = 1.0;
                            scale5 = 1.0;
                            damp = pdi * pdk;
                            expdamp = 0.0;
                            if (damp != 0.0) {
                                pgamma = min(pti, ptk);
                                rdamp = r / damp;
                                damp = -pgamma * rdamp * rdamp * rdamp;
                                if (damp > -50.0) {
                                    expdamp = exp(damp);
                                    scale3 = 1.0 - expdamp;
                                    scale5 = 1.0 - expdamp * (1.0 - damp);
                                }
                            }
                            rr3 = 1.0 / (r * r2);
                            rr5 = 3.0 * rr3 / r2;
                            rr3 *= (1.0 - scale3);
                            rr5 *= (1.0 - scale5);
                            bn2ukr = bn2 * ukr;
                            bn2pkr = bn2 * pkr;
                            rr5ukr = rr5 * ukr;
                            rr5pkr = rr5 * pkr;
                            fimx = -bn1 * ukx + bn2ukr * xr;
                            fimy = -bn1 * uky + bn2ukr * yr;
                            fimz = -bn1 * ukz + bn2ukr * zr;
                            fidx = -rr3 * ukx + rr5ukr * xr;
                            fidy = -rr3 * uky + rr5ukr * yr;
                            fidz = -rr3 * ukz + rr5ukr * zr;
                            pimx = -bn1 * pkx + bn2pkr * xr;
                            pimy = -bn1 * pky + bn2pkr * yr;
                            pimz = -bn1 * pkz + bn2pkr * zr;
                            pidx = -rr3 * pkx + rr5pkr * xr;
                            pidy = -rr3 * pky + rr5pkr * yr;
                            pidz = -rr3 * pkz + rr5pkr * zr;
                            fx += fimx - fidx;
                            fy += fimy - fidy;
                            fz += fimz - fidz;
                            px += pimx - pidx;
                            py += pimy - pidy;
                            pz += pimz - pidz;
                        }
                        fx_local[i] += fx;
                        fy_local[i] += fy;
                        fz_local[i] += fz;
                        fxp_local[i] += px;
                        fyp_local[i] += py;
                        fzp_local[i] += pz;
                    }
                }
            }

            /**
             * Reduce this thread's field contribution into the shared arrays.
             */
            @Override
            public void finish() {
                sharedField[0].reduce(fx_local, DoubleOp.SUM);
                sharedField[1].reduce(fy_local, DoubleOp.SUM);
                sharedField[2].reduce(fz_local, DoubleOp.SUM);
                sharedFieldp[0].reduce(fxp_local, DoubleOp.SUM);
                sharedFieldp[1].reduce(fyp_local, DoubleOp.SUM);
                sharedFieldp[2].reduce(fzp_local, DoubleOp.SUM);
            }
        }
    }

    /**
     * The Real Space Gradient Region class parallelizes evaluation of the real
     * space energy and gradients using an array of Real Space Gradient Loops.
     */
    private class RealSpaceEnergyRegion extends ParallelRegion {

        private final SharedDouble sharedPermanentEnergy;
        private final SharedDouble sharedPolarizationEnergy;
        private final SharedInteger sharedInteractions;
        private final RealSpaceEnergyLoop realSpaceEnergyLoop[];
        private boolean gradient;
        private long overheadTime;

        public RealSpaceEnergyRegion(int nt) {
            super();
            sharedPermanentEnergy = new SharedDouble();
            sharedPolarizationEnergy = new SharedDouble();
            sharedInteractions = new SharedInteger();
            realSpaceEnergyLoop = new RealSpaceEnergyLoop[nt];
            for (int i = 0; i < nt; i++) {
                realSpaceEnergyLoop[i] = new RealSpaceEnergyLoop();
            }
        }

        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }

        public double getPermanentEnergy() {
            return sharedPermanentEnergy.get();
        }

        public double getPolarizationEnergy() {
            return sharedPolarizationEnergy.get();
        }

        public int getInteractions() {
            return sharedInteractions.get();
        }

        @Override
        public void start() {
            overheadTime = System.nanoTime();
            sharedPermanentEnergy.set(0.0);
            sharedPolarizationEnergy.set(0.0);
            sharedInteractions.set(0);
        }

        @Override
        public void run() {
            try {
                int threadIndex = getThreadIndex();
                realSpaceEnergyLoop[threadIndex].setGradient(gradient);
                execute(0, nAtoms - 1, realSpaceEnergyLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing the real space energy in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        @Override
        public void finish() {
            long computeTime = 0;
            for (int i = 0; i < maxThreads; i++) {
                computeTime += realSpaceEnergyLoop[i].getComputeTime();
            }
            overheadTime = System.nanoTime() - overheadTime;
            overheadTime = overheadTime - computeTime / maxThreads;
            //double compute = (double) computeTime / threadCount * toSeconds;
            //double overhead = (double) overheadTime * toSeconds;
            //double efficiency = compute / (compute + overhead) * 100;
            /*
             * logger.info(format("Real Space Energy Parallel Performance\n"
             * + "Avg. Compute Time  %10.3f (sec)\n" +
             * "Overhead Time      %10.3f (sec)\n" +
             * "Efficiency         %10.3f\n", compute, overhead, efficiency));
             */
        }

        /**
         * The Real Space Gradient Loop class contains methods and thread local
         * variables to parallelize the evaluation of the real space permanent
         * and polarization energies and gradients.
         */
        private class RealSpaceEnergyLoop extends IntegerForLoop {

            private long computeTime;
            private boolean gradient;
            private double ci;
            private double dix, diy, diz;
            private double qixx, qiyy, qizz, qixy, qixz, qiyz;
            private double ck;
            private double dkx, dky, dkz;
            private double qkxx, qkyy, qkzz, qkxy, qkxz, qkyz;
            private double uix, uiy, uiz;
            private double pix, piy, piz;
            private double xr, yr, zr;
            private double ukx, uky, ukz;
            private double pkx, pky, pkz;
            private double bn0, bn1, bn2, bn3, bn4, bn5;
            private double rr1, rr3, rr5, rr7, rr9, rr11;
            private double scale, scale3, scale5, scale7;
            private double scalep, scaled;
            private double ddsc3x, ddsc3y, ddsc3z;
            private double ddsc5x, ddsc5y, ddsc5z;
            private double ddsc7x, ddsc7y, ddsc7z;
            private double permanentEnergy;
            private double inducedEnergy;
            private int i, k, iSymm, count;
            private final double dx_local[];
            private final double gxi_local[];
            private final double gyi_local[];
            private final double gzi_local[];
            private final double gxk_local[];
            private final double gyk_local[];
            private final double gzk_local[];
            private final double txi_local[];
            private final double tyi_local[];
            private final double tzi_local[];
            private final double txk_local[];
            private final double tyk_local[];
            private final double tzk_local[];
            private final double masking_local[];
            private final double maskingp_local[];
            private final double maskingd_local[];

            public RealSpaceEnergyLoop() {
                super();
                gxi_local = new double[nAtoms];
                gyi_local = new double[nAtoms];
                gzi_local = new double[nAtoms];
                gxk_local = new double[nAtoms];
                gyk_local = new double[nAtoms];
                gzk_local = new double[nAtoms];
                txi_local = new double[nAtoms];
                tyi_local = new double[nAtoms];
                tzi_local = new double[nAtoms];
                txk_local = new double[nAtoms];
                tyk_local = new double[nAtoms];
                tzk_local = new double[nAtoms];
                masking_local = new double[nAtoms];
                maskingp_local = new double[nAtoms];
                maskingd_local = new double[nAtoms];
                dx_local = new double[3];
            }

            public long getComputeTime() {
                return computeTime;
            }

            public void setGradient(boolean gradient) {
                this.gradient = gradient;
            }

            @Override
            public IntegerSchedule schedule() {
                return pairWiseSchedule;
            }

            @Override
            public void start() {
                permanentEnergy = 0.0;
                inducedEnergy = 0.0;
                count = 0;
                for (int j = 0; j < nAtoms; j++) {
                    masking_local[j] = 1.0;
                    maskingp_local[j] = 1.0;
                    maskingd_local[j] = 1.0;
                }
                if (gradient) {
                    for (int j = 0; j < nAtoms; j++) {
                        gxi_local[j] = 0.0;
                        gyi_local[j] = 0.0;
                        gzi_local[j] = 0.0;
                        txi_local[j] = 0.0;
                        tyi_local[j] = 0.0;
                        tzi_local[j] = 0.0;
                    }
                }
                computeTime = 0;
            }

            @Override
            public void run(int lb, int ub) {
                long startTime = System.nanoTime();
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                for (iSymm = 0; iSymm < nSymm; iSymm++) {
                    if (gradient && iSymm > 0) {
                        for (int j = 0; j < nAtoms; j++) {
                            gxk_local[j] = 0.0;
                            gyk_local[j] = 0.0;
                            gzk_local[j] = 0.0;
                            txk_local[j] = 0.0;
                            tyk_local[j] = 0.0;
                            tzk_local[j] = 0.0;
                        }
                    }
                    realSpaceChunk(lb, ub);
                    if (gradient && iSymm > 0) {
                        /**
                         * Apply the symmetry rotation to the force and torque
                         * of atoms in a symmetry mate.
                         */
                        SymOp symOp = symOps.get(iSymm);
                        crystal.applySymRot(nAtoms, gxk_local, gyk_local,
                                gzk_local, gxk_local, gyk_local, gzk_local,
                                symOp);
                        crystal.applySymRot(nAtoms, txk_local, tyk_local,
                                tzk_local, txk_local, tyk_local, tzk_local,
                                symOp);
                        /**
                         * The two force and torque arrays can now be condensed
                         * into single arrays.
                         */
                        for (int j = 0; j < nAtoms; j++) {
                            gxi_local[j] += gxk_local[j];
                            gyi_local[j] += gyk_local[j];
                            gzi_local[j] += gzk_local[j];
                            txi_local[j] += txk_local[j];
                            tyi_local[j] += tyk_local[j];
                            tzi_local[j] += tzk_local[j];
                        }
                    }
                }
                computeTime += System.nanoTime() - startTime;
            }

            @Override
            public void finish() {
                sharedInteractions.addAndGet(count);
                sharedPermanentEnergy.addAndGet(permanentEnergy * electric);
                sharedPolarizationEnergy.addAndGet(inducedEnergy * electric);
                if (gradient) {
                    for (int j = 0; j < nAtoms; j++) {
                        gxi_local[j] *= electric;
                        gyi_local[j] *= electric;
                        gzi_local[j] *= electric;
                        txi_local[j] *= electric;
                        tyi_local[j] *= electric;
                        tzi_local[j] *= electric;
                    }
                    /**
                     * Reduce the force and torque contributions computed by the
                     * current thread into the shared arrays.
                     */
                    sharedGrad[0].reduce(gxi_local, DoubleOp.SUM);
                    sharedGrad[1].reduce(gyi_local, DoubleOp.SUM);
                    sharedGrad[2].reduce(gzi_local, DoubleOp.SUM);
                    sharedTorque[0].reduce(txi_local, DoubleOp.SUM);
                    sharedTorque[1].reduce(tyi_local, DoubleOp.SUM);
                    sharedTorque[2].reduce(tzi_local, DoubleOp.SUM);
                }
            }

            /**
             * Evaluate the real space permanent energy and polarization energy
             * for a chunk of atoms.
             *
             * @param lb The lower bound of the chunk.
             * @param ub The upper bound of the chunk.
             */
            private void realSpaceChunk(final int lb, final int ub) {
                final double x[] = coordinates[0][0];
                final double y[] = coordinates[0][1];
                final double z[] = coordinates[0][2];
                final double mpole[][] = globalMultipole[0];
                final double ind[][] = inducedDipole[0];
                final double indp[][] = inducedDipolep[0];
                final int lists[][] = ewaldLists[iSymm];
                final double neighborX[] = coordinates[iSymm][0];
                final double neighborY[] = coordinates[iSymm][1];
                final double neighborZ[] = coordinates[iSymm][2];
                final double neighborMultipole[][] = globalMultipole[iSymm];
                final double neighborInducedDipole[][] = inducedDipole[iSymm];
                final double neighborInducedDipolep[][] = inducedDipolep[iSymm];
                double asymmetric = 1.0;
                if (iSymm > 0) {
                    asymmetric = 0.5;
                }
                for (i = lb; i <= ub; i++) {
                    final Atom ai = atoms[i];
                    if (iSymm == 0) {
                        for (Atom ak : ai.get1_5s()) {
                            /*
                             dx_local[0] = ak.getX() - ai.getX();
                             dx_local[1] = ak.getY() - ai.getY();
                             dx_local[2] = ak.getZ() - ai.getZ();
                             final double r2 = crystal.image(dx_local);
                             if (r2 > off2) {
                             logger.warning("Needed Ewald interaction is outside the cutoff: " + ai + ak);
                             } */
                            masking_local[ak.xyzIndex - 1] = m15scale;
                        }
                        for (Torsion torsion : ai.getTorsions()) {
                            Atom ak = torsion.get1_4(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = m14scale;
                                for (int j : ip11[i]) {
                                    if (j == index) {
                                        maskingp_local[index] = 0.5;
                                    }
                                }
                            }
                        }
                        for (Angle angle : ai.getAngles()) {
                            Atom ak = angle.get1_3(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = m13scale;
                                maskingp_local[index] = p13scale;
                            }
                        }
                        for (Bond bond : ai.getFFXBonds()) {
                            int index = bond.get1_2(ai).xyzIndex - 1;
                            masking_local[index] = m12scale;
                            maskingp_local[index] = p12scale;
                        }
                        for (int j : ip11[i]) {
                            maskingd_local[j] = d11scale;
                        }
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double globalMultipolei[] = mpole[i];
                    final double inducedDipolei[] = ind[i];
                    final double inducedDipolepi[] = indp[i];
                    ci = globalMultipolei[t000];
                    dix = globalMultipolei[t100];
                    diy = globalMultipolei[t010];
                    diz = globalMultipolei[t001];
                    qixx = globalMultipolei[t200] / 3.0;
                    qiyy = globalMultipolei[t020] / 3.0;
                    qizz = globalMultipolei[t002] / 3.0;
                    qixy = globalMultipolei[t110] / 3.0;
                    qixz = globalMultipolei[t101] / 3.0;
                    qiyz = globalMultipolei[t011] / 3.0;
                    uix = inducedDipolei[0];
                    uiy = inducedDipolei[1];
                    uiz = inducedDipolei[2];
                    pix = inducedDipolepi[0];
                    piy = inducedDipolepi[1];
                    piz = inducedDipolepi[2];
                    final double pdi = pdamp[i];
                    final double pti = thole[i];
                    final int list[] = lists[i];
                    final int npair = ewaldCounts[iSymm][i];
                    for (int j = 0; j < npair; j++) {
                        k = list[j];
                        final double xk = neighborX[k];
                        final double yk = neighborY[k];
                        final double zk = neighborZ[k];
                        dx_local[0] = xk - xi;
                        dx_local[1] = yk - yi;
                        dx_local[2] = zk - zi;
                        final double r2 = crystal.image(dx_local);
                        xr = dx_local[0];
                        yr = dx_local[1];
                        zr = dx_local[2];
                        final double globalMultipolek[] = neighborMultipole[k];
                        ck = globalMultipolek[t000];
                        dkx = globalMultipolek[t100];
                        dky = globalMultipolek[t010];
                        dkz = globalMultipolek[t001];
                        qkxx = globalMultipolek[t200] / 3.0;
                        qkyy = globalMultipolek[t020] / 3.0;
                        qkzz = globalMultipolek[t002] / 3.0;
                        qkxy = globalMultipolek[t110] / 3.0;
                        qkxz = globalMultipolek[t101] / 3.0;
                        qkyz = globalMultipolek[t011] / 3.0;
                        final double inducedDipolek[] = neighborInducedDipole[k];
                        ukx = inducedDipolek[0];
                        uky = inducedDipolek[1];
                        ukz = inducedDipolek[2];
                        final double inducedDipolepk[] = neighborInducedDipolep[k];
                        pkx = inducedDipolepk[0];
                        pky = inducedDipolepk[1];
                        pkz = inducedDipolepk[2];
                        final double pdk = pdamp[k];
                        final double ptk = thole[k];
                        scale = masking_local[k];
                        scalep = maskingp_local[k];
                        scaled = maskingd_local[k];
                        scale3 = 1.0;
                        scale5 = 1.0;
                        scale7 = 1.0;
                        final double r = sqrt(r2);
                        final double ralpha = aewald * r;
                        bn0 = erfc(ralpha) / r;
                        double alsq2n = piEwald;
                        final double exp2a = exp(-ralpha * ralpha);
                        alsq2n = alsq2 * alsq2n;
                        bn1 = (bn0 + alsq2n * exp2a) / r2;
                        alsq2n = alsq2 * alsq2n;
                        bn2 = (3.0 * bn1 + alsq2n * exp2a) / r2;
                        alsq2n = alsq2 * alsq2n;
                        bn3 = (5.0 * bn2 + alsq2n * exp2a) / r2;
                        alsq2n = alsq2 * alsq2n;
                        bn4 = (7.0 * bn3 + alsq2n * exp2a) / r2;
                        alsq2n = alsq2 * alsq2n;
                        bn5 = (9.0 * bn4 + alsq2n * exp2a) / r2;
                        rr1 = 1.0 / r;
                        rr3 = rr1 / r2;
                        rr5 = 3.0 * rr3 / r2;
                        rr7 = 5.0 * rr5 / r2;
                        rr9 = 7.0 * rr7 / r2;
                        rr11 = 9.0 * rr9 / r2;
                        ddsc3x = 0.0;
                        ddsc3y = 0.0;
                        ddsc3z = 0.0;
                        ddsc5x = 0.0;
                        ddsc5y = 0.0;
                        ddsc5z = 0.0;
                        ddsc7x = 0.0;
                        ddsc7y = 0.0;
                        ddsc7z = 0.0;
                        double damp = pdi * pdk;
                        if (damp != 0.0) {
                            final double pgamma = min(pti, ptk);
                            final double rdamp = r / damp;
                            damp = -pgamma * rdamp * rdamp * rdamp;
                            if (damp > -50.0) {
                                final double expdamp = exp(damp);
                                scale3 = 1.0 - expdamp;
                                scale5 = 1.0 - expdamp * (1.0 - damp);
                                scale7 = 1.0 - expdamp * (1.0 - damp + 0.6 * damp * damp);
                                final double temp3 = -3.0 * damp * expdamp / r2;
                                final double temp5 = -damp;
                                final double temp7 = -0.2 - 0.6 * damp;
                                ddsc3x = temp3 * xr;
                                ddsc3y = temp3 * yr;
                                ddsc3z = temp3 * zr;
                                ddsc5x = temp5 * ddsc3x;
                                ddsc5y = temp5 * ddsc3y;
                                ddsc5z = temp5 * ddsc3z;
                                ddsc7x = temp7 * ddsc5x;
                                ddsc7y = temp7 * ddsc5y;
                                ddsc7z = temp7 * ddsc5z;
                            }
                        }
                        double e = asymmetric * realSpacePermanentPair();
                        permanentEnergy += e;
                        if (polarization != Polarization.NONE) {
                            inducedEnergy += asymmetric * realSpacePolarizationPair();
                        }
                        /*
                         if (i == 0) {
                         System.out.println(format("%5d %10.5f %10.5f %5d %10.5f %10.5f %10.5f", k + 1, e * electric / asymmetric, r, iSymm + 1, xk, yk, zk));
                         }
                         if (k == 0) {
                         System.out.println(format("%5d %10.5f %10.5f %5d %10.5f %10.5f %10.5f", i + 1, e * electric / asymmetric, r, iSymm + 1, xi, yi, zi));
                         }
                         */
                        count++;
                    }
                    if (iSymm == 0) {
                        for (Atom ak : ai.get1_5s()) {
                            int index = ak.xyzIndex - 1;
                            masking_local[index] = 1.0;
                        }
                        for (Torsion torsion : ai.getTorsions()) {
                            Atom ak = torsion.get1_4(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = 1.0;
                                maskingp_local[index] = 1.0;
                            }
                        }
                        for (Angle angle : ai.getAngles()) {
                            Atom ak = angle.get1_3(ai);
                            if (ak != null) {
                                int index = ak.xyzIndex - 1;
                                masking_local[index] = 1.0;
                                maskingp_local[index] = 1.0;
                            }
                        }
                        for (Bond bond : ai.getFFXBonds()) {
                            int index = bond.get1_2(ai).xyzIndex - 1;
                            masking_local[index] = 1.0;
                            maskingp_local[index] = 1.0;
                        }
                        for (int j : ip11[i]) {
                            maskingd_local[j] = 1.0;
                        }
                    }
                }
            }

            /**
             * Evaluate the real space permanent energy for a pair of multipole
             * sites.
             *
             * @return the permanent multipole energy.
             */
            private double realSpacePermanentPair() {
                final double dixdkx = diy * dkz - diz * dky;
                final double dixdky = diz * dkx - dix * dkz;
                final double dixdkz = dix * dky - diy * dkx;
                final double dixrx = diy * zr - diz * yr;
                final double dixry = diz * xr - dix * zr;
                final double dixrz = dix * yr - diy * xr;
                final double dkxrx = dky * zr - dkz * yr;
                final double dkxry = dkz * xr - dkx * zr;
                final double dkxrz = dkx * yr - dky * xr;
                final double qirx = qixx * xr + qixy * yr + qixz * zr;
                final double qiry = qixy * xr + qiyy * yr + qiyz * zr;
                final double qirz = qixz * xr + qiyz * yr + qizz * zr;
                final double qkrx = qkxx * xr + qkxy * yr + qkxz * zr;
                final double qkry = qkxy * xr + qkyy * yr + qkyz * zr;
                final double qkrz = qkxz * xr + qkyz * yr + qkzz * zr;
                final double qiqkrx = qixx * qkrx + qixy * qkry + qixz * qkrz;
                final double qiqkry = qixy * qkrx + qiyy * qkry + qiyz * qkrz;
                final double qiqkrz = qixz * qkrx + qiyz * qkry + qizz * qkrz;
                final double qkqirx = qkxx * qirx + qkxy * qiry + qkxz * qirz;
                final double qkqiry = qkxy * qirx + qkyy * qiry + qkyz * qirz;
                final double qkqirz = qkxz * qirx + qkyz * qiry + qkzz * qirz;
                final double qixqkx = qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy - qiyz * qkyy - qizz * qkyz;
                final double qixqky = qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz - qixy * qkyz - qixz * qkzz;
                final double qixqkz = qixx * qkxy + qixy * qkyy + qixz * qkyz - qixy * qkxx - qiyy * qkxy - qiyz * qkxz;
                final double rxqirx = yr * qirz - zr * qiry;
                final double rxqiry = zr * qirx - xr * qirz;
                final double rxqirz = xr * qiry - yr * qirx;
                final double rxqkrx = yr * qkrz - zr * qkry;
                final double rxqkry = zr * qkrx - xr * qkrz;
                final double rxqkrz = xr * qkry - yr * qkrx;
                final double rxqikrx = yr * qiqkrz - zr * qiqkry;
                final double rxqikry = zr * qiqkrx - xr * qiqkrz;
                final double rxqikrz = xr * qiqkry - yr * qiqkrx;
                final double rxqkirx = yr * qkqirz - zr * qkqiry;
                final double rxqkiry = zr * qkqirx - xr * qkqirz;
                final double rxqkirz = xr * qkqiry - yr * qkqirx;
                final double qkrxqirx = qkry * qirz - qkrz * qiry;
                final double qkrxqiry = qkrz * qirx - qkrx * qirz;
                final double qkrxqirz = qkrx * qiry - qkry * qirx;
                final double qidkx = qixx * dkx + qixy * dky + qixz * dkz;
                final double qidky = qixy * dkx + qiyy * dky + qiyz * dkz;
                final double qidkz = qixz * dkx + qiyz * dky + qizz * dkz;
                final double qkdix = qkxx * dix + qkxy * diy + qkxz * diz;
                final double qkdiy = qkxy * dix + qkyy * diy + qkyz * diz;
                final double qkdiz = qkxz * dix + qkyz * diy + qkzz * diz;
                final double dixqkrx = diy * qkrz - diz * qkry;
                final double dixqkry = diz * qkrx - dix * qkrz;
                final double dixqkrz = dix * qkry - diy * qkrx;
                final double dkxqirx = dky * qirz - dkz * qiry;
                final double dkxqiry = dkz * qirx - dkx * qirz;
                final double dkxqirz = dkx * qiry - dky * qirx;
                final double rxqidkx = yr * qidkz - zr * qidky;
                final double rxqidky = zr * qidkx - xr * qidkz;
                final double rxqidkz = xr * qidky - yr * qidkx;
                final double rxqkdix = yr * qkdiz - zr * qkdiy;
                final double rxqkdiy = zr * qkdix - xr * qkdiz;
                final double rxqkdiz = xr * qkdiy - yr * qkdix;
                /**
                 * Calculate the scalar products for permanent multipoles.
                 */
                final double sc2 = dix * dkx + diy * dky + diz * dkz;
                final double sc3 = dix * xr + diy * yr + diz * zr;
                final double sc4 = dkx * xr + dky * yr + dkz * zr;
                final double sc5 = qirx * xr + qiry * yr + qirz * zr;
                final double sc6 = qkrx * xr + qkry * yr + qkrz * zr;
                final double sc7 = qirx * dkx + qiry * dky + qirz * dkz;
                final double sc8 = qkrx * dix + qkry * diy + qkrz * diz;
                final double sc9 = qirx * qkrx + qiry * qkry + qirz * qkrz;
                final double sc10 = 2.0 * (qixy * qkxy + qixz * qkxz + qiyz * qkyz) + qixx * qkxx + qiyy * qkyy + qizz * qkzz;
                /**
                 * Calculate the gl functions for permanent multipoles.
                 */
                final double gl0 = ci * ck;
                final double gl1 = ck * sc3 - ci * sc4;
                final double gl2 = ci * sc6 + ck * sc5 - sc3 * sc4;
                final double gl3 = sc3 * sc6 - sc4 * sc5;
                final double gl4 = sc5 * sc6;
                final double gl5 = -4.0 * sc9;
                final double gl6 = sc2;
                final double gl7 = 2.0 * (sc7 - sc8);
                final double gl8 = 2.0 * sc10;
                /**
                 * Compute the energy contributions for this interaction.
                 */
                double e = gl0 * bn0 + (gl1 + gl6) * bn1 + (gl2 + gl7 + gl8) * bn2 + (gl3 + gl5) * bn3 + gl4 * bn4;
                final double efix = gl0 * rr1 + (gl1 + gl6) * rr3 + (gl2 + gl7 + gl8) * rr5 + (gl3 + gl5) * rr7 + gl4 * rr9;
                e = e - efix * (1.0 - scale);
                if (!gradient) {
                    return e;
                }
                boolean dorl = false;
                if (scale != 1.0) {
                    dorl = true;
                }
                final double gf1 = bn1 * gl0 + bn2 * (gl1 + gl6) + bn3 * (gl2 + gl7 + gl8) + bn4 * (gl3 + gl5) + bn5 * gl4;
                final double gf2 = -ck * bn1 + sc4 * bn2 - sc6 * bn3;
                final double gf3 = ci * bn1 + sc3 * bn2 + sc5 * bn3;
                final double gf4 = 2.0 * bn2;
                final double gf5 = 2.0 * (-ck * bn2 + sc4 * bn3 - sc6 * bn4);
                final double gf6 = 2.0 * (-ci * bn2 - sc3 * bn3 - sc5 * bn4);
                final double gf7 = 4.0 * bn3;
                /*
                 * Get the permanent force with screening.
                 */
                double ftm2x = gf1 * xr + gf2 * dix + gf3 * dkx + gf4 * (qkdix - qidkx) + gf5 * qirx + gf6 * qkrx + gf7 * (qiqkrx + qkqirx);
                double ftm2y = gf1 * yr + gf2 * diy + gf3 * dky + gf4 * (qkdiy - qidky) + gf5 * qiry + gf6 * qkry + gf7 * (qiqkry + qkqiry);
                double ftm2z = gf1 * zr + gf2 * diz + gf3 * dkz + gf4 * (qkdiz - qidkz) + gf5 * qirz + gf6 * qkrz + gf7 * (qiqkrz + qkqirz);
                /*
                 * Get the permanent torque with screening.
                 */
                double ttm2x = -bn1 * dixdkx + gf2 * dixrx + gf4 * (dixqkrx + dkxqirx + rxqidkx - 2.0 * qixqkx) - gf5 * rxqirx - gf7 * (rxqikrx + qkrxqirx);
                double ttm2y = -bn1 * dixdky + gf2 * dixry + gf4 * (dixqkry + dkxqiry + rxqidky - 2.0 * qixqky) - gf5 * rxqiry - gf7 * (rxqikry + qkrxqiry);
                double ttm2z = -bn1 * dixdkz + gf2 * dixrz + gf4 * (dixqkrz + dkxqirz + rxqidkz - 2.0 * qixqkz) - gf5 * rxqirz - gf7 * (rxqikrz + qkrxqirz);
                double ttm3x = bn1 * dixdkx + gf3 * dkxrx - gf4 * (dixqkrx + dkxqirx + rxqkdix - 2.0 * qixqkx) - gf6 * rxqkrx - gf7 * (rxqkirx - qkrxqirx);
                double ttm3y = bn1 * dixdky + gf3 * dkxry - gf4 * (dixqkry + dkxqiry + rxqkdiy - 2.0 * qixqky) - gf6 * rxqkry - gf7 * (rxqkiry - qkrxqiry);
                double ttm3z = bn1 * dixdkz + gf3 * dkxrz - gf4 * (dixqkrz + dkxqirz + rxqkdiz - 2.0 * qixqkz) - gf6 * rxqkrz - gf7 * (rxqkirz - qkrxqirz);
                double ftm2rx = 0.0;
                double ftm2ry = 0.0;
                double ftm2rz = 0.0;
                double ttm2rx = 0.0;
                double ttm2ry = 0.0;
                double ttm2rz = 0.0;
                double ttm3rx = 0.0;
                double ttm3ry = 0.0;
                double ttm3rz = 0.0;
                if (dorl) {
                    final double gfr1 = rr3 * gl0 + rr5 * (gl1 + gl6) + rr7 * (gl2 + gl7 + gl8) + rr9 * (gl3 + gl5) + rr11 * gl4;
                    final double gfr2 = -ck * rr3 + sc4 * rr5 - sc6 * rr7;
                    final double gfr3 = ci * rr3 + sc3 * rr5 + sc5 * rr7;
                    final double gfr4 = 2.0 * rr5;
                    final double gfr5 = 2.0 * (-ck * rr5 + sc4 * rr7 - sc6 * rr9);
                    final double gfr6 = 2.0 * (-ci * rr5 - sc3 * rr7 - sc5 * rr9);
                    final double gfr7 = 4.0 * rr7;
                    /*
                     * Get the permanent force without screening.
                     */
                    ftm2rx = gfr1 * xr + gfr2 * dix + gfr3 * dkx + gfr4 * (qkdix - qidkx) + gfr5 * qirx + gfr6 * qkrx + gfr7 * (qiqkrx + qkqirx);
                    ftm2ry = gfr1 * yr + gfr2 * diy + gfr3 * dky + gfr4 * (qkdiy - qidky) + gfr5 * qiry + gfr6 * qkry + gfr7 * (qiqkry + qkqiry);
                    ftm2rz = gfr1 * zr + gfr2 * diz + gfr3 * dkz + gfr4 * (qkdiz - qidkz) + gfr5 * qirz + gfr6 * qkrz + gfr7 * (qiqkrz + qkqirz);
                    /*
                     * Get the permanent torque without screening.
                     */
                    ttm2rx = -rr3 * dixdkx + gfr2 * dixrx + gfr4 * (dixqkrx + dkxqirx + rxqidkx - 2.0 * qixqkx) - gfr5 * rxqirx - gfr7 * (rxqikrx + qkrxqirx);
                    ttm2ry = -rr3 * dixdky + gfr2 * dixry + gfr4 * (dixqkry + dkxqiry + rxqidky - 2.0 * qixqky) - gfr5 * rxqiry - gfr7 * (rxqikry + qkrxqiry);
                    ttm2rz = -rr3 * dixdkz + gfr2 * dixrz + gfr4 * (dixqkrz + dkxqirz + rxqidkz - 2.0 * qixqkz) - gfr5 * rxqirz - gfr7 * (rxqikrz + qkrxqirz);
                    ttm3rx = rr3 * dixdkx + gfr3 * dkxrx - gfr4 * (dixqkrx + dkxqirx + rxqkdix - 2.0 * qixqkx) - gfr6 * rxqkrx - gfr7 * (rxqkirx - qkrxqirx);
                    ttm3ry = rr3 * dixdky + gfr3 * dkxry - gfr4 * (dixqkry + dkxqiry + rxqkdiy - 2.0 * qixqky) - gfr6 * rxqkry - gfr7 * (rxqkiry - qkrxqiry);
                    ttm3rz = rr3 * dixdkz + gfr3 * dkxrz - gfr4 * (dixqkrz + dkxqirz + rxqkdiz - 2.0 * qixqkz) - gfr6 * rxqkrz - gfr7 * (rxqkirz - qkrxqirz);
                }
                /**
                 * Handle the case where scaling is used.
                 */
                final double scale1 = 1.0 - scale;
                ftm2x = ftm2x - scale1 * ftm2rx;
                ftm2y = ftm2y - scale1 * ftm2ry;
                ftm2z = ftm2z - scale1 * ftm2rz;
                ttm2x = ttm2x - scale1 * ttm2rx;
                ttm2y = ttm2y - scale1 * ttm2ry;
                ttm2z = ttm2z - scale1 * ttm2rz;
                ttm3x = ttm3x - scale1 * ttm3rx;
                ttm3y = ttm3y - scale1 * ttm3ry;
                ttm3z = ttm3z - scale1 * ttm3rz;
                if (iSymm == 0) {
                    gxi_local[i] += ftm2x;
                    gyi_local[i] += ftm2y;
                    gzi_local[i] += ftm2z;
                    gxi_local[k] -= ftm2x;
                    gyi_local[k] -= ftm2y;
                    gzi_local[k] -= ftm2z;
                    txi_local[i] += ttm2x;
                    tyi_local[i] += ttm2y;
                    tzi_local[i] += ttm2z;
                    txi_local[k] += ttm3x;
                    tyi_local[k] += ttm3y;
                    tzi_local[k] += ttm3z;
                } else {
                    gxi_local[i] += 0.5 * ftm2x;
                    gyi_local[i] += 0.5 * ftm2y;
                    gzi_local[i] += 0.5 * ftm2z;
                    gxk_local[k] -= 0.5 * ftm2x;
                    gyk_local[k] -= 0.5 * ftm2y;
                    gzk_local[k] -= 0.5 * ftm2z;
                    txi_local[i] += 0.5 * ttm2x;
                    tyi_local[i] += 0.5 * ttm2y;
                    tzi_local[i] += 0.5 * ttm2z;
                    txk_local[k] += 0.5 * ttm3x;
                    tyk_local[k] += 0.5 * ttm3y;
                    tzk_local[k] += 0.5 * ttm3z;
                }
                return e;
            }

            /**
             * Evaluate the polarization energy for a pair of polarizable
             * multipole sites.
             *
             * @return the polarization energy.
             */
            private double realSpacePolarizationPair() {
                final double dsc3 = 1.0 - scale3 * scaled;
                final double dsc5 = 1.0 - scale5 * scaled;
                final double dsc7 = 1.0 - scale7 * scaled;
                final double psc3 = 1.0 - scale3 * scalep;
                final double psc5 = 1.0 - scale5 * scalep;
                final double psc7 = 1.0 - scale7 * scalep;
                final double usc3 = 1.0 - scale3;
                final double usc5 = 1.0 - scale5;
                final double dixukx = diy * ukz - diz * uky;
                final double dixuky = diz * ukx - dix * ukz;
                final double dixukz = dix * uky - diy * ukx;
                final double dkxuix = dky * uiz - dkz * uiy;
                final double dkxuiy = dkz * uix - dkx * uiz;
                final double dkxuiz = dkx * uiy - dky * uix;
                final double dixukpx = diy * pkz - diz * pky;
                final double dixukpy = diz * pkx - dix * pkz;
                final double dixukpz = dix * pky - diy * pkx;
                final double dkxuipx = dky * piz - dkz * piy;
                final double dkxuipy = dkz * pix - dkx * piz;
                final double dkxuipz = dkx * piy - dky * pix;
                final double dixrx = diy * zr - diz * yr;
                final double dixry = diz * xr - dix * zr;
                final double dixrz = dix * yr - diy * xr;
                final double dkxrx = dky * zr - dkz * yr;
                final double dkxry = dkz * xr - dkx * zr;
                final double dkxrz = dkx * yr - dky * xr;
                final double qirx = qixx * xr + qixy * yr + qixz * zr;
                final double qiry = qixy * xr + qiyy * yr + qiyz * zr;
                final double qirz = qixz * xr + qiyz * yr + qizz * zr;
                final double qkrx = qkxx * xr + qkxy * yr + qkxz * zr;
                final double qkry = qkxy * xr + qkyy * yr + qkyz * zr;
                final double qkrz = qkxz * xr + qkyz * yr + qkzz * zr;
                final double rxqirx = yr * qirz - zr * qiry;
                final double rxqiry = zr * qirx - xr * qirz;
                final double rxqirz = xr * qiry - yr * qirx;
                final double rxqkrx = yr * qkrz - zr * qkry;
                final double rxqkry = zr * qkrx - xr * qkrz;
                final double rxqkrz = xr * qkry - yr * qkrx;
                final double qiukx = qixx * ukx + qixy * uky + qixz * ukz;
                final double qiuky = qixy * ukx + qiyy * uky + qiyz * ukz;
                final double qiukz = qixz * ukx + qiyz * uky + qizz * ukz;
                final double qkuix = qkxx * uix + qkxy * uiy + qkxz * uiz;
                final double qkuiy = qkxy * uix + qkyy * uiy + qkyz * uiz;
                final double qkuiz = qkxz * uix + qkyz * uiy + qkzz * uiz;
                final double qiukpx = qixx * pkx + qixy * pky + qixz * pkz;
                final double qiukpy = qixy * pkx + qiyy * pky + qiyz * pkz;
                final double qiukpz = qixz * pkx + qiyz * pky + qizz * pkz;
                final double qkuipx = qkxx * pix + qkxy * piy + qkxz * piz;
                final double qkuipy = qkxy * pix + qkyy * piy + qkyz * piz;
                final double qkuipz = qkxz * pix + qkyz * piy + qkzz * piz;
                final double uixqkrx = uiy * qkrz - uiz * qkry;
                final double uixqkry = uiz * qkrx - uix * qkrz;
                final double uixqkrz = uix * qkry - uiy * qkrx;
                final double ukxqirx = uky * qirz - ukz * qiry;
                final double ukxqiry = ukz * qirx - ukx * qirz;
                final double ukxqirz = ukx * qiry - uky * qirx;
                final double uixqkrpx = piy * qkrz - piz * qkry;
                final double uixqkrpy = piz * qkrx - pix * qkrz;
                final double uixqkrpz = pix * qkry - piy * qkrx;
                final double ukxqirpx = pky * qirz - pkz * qiry;
                final double ukxqirpy = pkz * qirx - pkx * qirz;
                final double ukxqirpz = pkx * qiry - pky * qirx;
                final double rxqiukx = yr * qiukz - zr * qiuky;
                final double rxqiuky = zr * qiukx - xr * qiukz;
                final double rxqiukz = xr * qiuky - yr * qiukx;
                final double rxqkuix = yr * qkuiz - zr * qkuiy;
                final double rxqkuiy = zr * qkuix - xr * qkuiz;
                final double rxqkuiz = xr * qkuiy - yr * qkuix;
                final double rxqiukpx = yr * qiukpz - zr * qiukpy;
                final double rxqiukpy = zr * qiukpx - xr * qiukpz;
                final double rxqiukpz = xr * qiukpy - yr * qiukpx;
                final double rxqkuipx = yr * qkuipz - zr * qkuipy;
                final double rxqkuipy = zr * qkuipx - xr * qkuipz;
                final double rxqkuipz = xr * qkuipy - yr * qkuipx;
                /**
                 * Calculate the scalar products for permanent multipoles.
                 */
                final double sc3 = dix * xr + diy * yr + diz * zr;
                final double sc4 = dkx * xr + dky * yr + dkz * zr;
                final double sc5 = qirx * xr + qiry * yr + qirz * zr;
                final double sc6 = qkrx * xr + qkry * yr + qkrz * zr;
                /**
                 * Calculate the scalar products for polarization components.
                 */
                final double sci1 = uix * dkx + uiy * dky + uiz * dkz + dix * ukx + diy * uky + diz * ukz;
                final double sci3 = uix * xr + uiy * yr + uiz * zr;
                final double sci4 = ukx * xr + uky * yr + ukz * zr;
                final double sci7 = qirx * ukx + qiry * uky + qirz * ukz;
                final double sci8 = qkrx * uix + qkry * uiy + qkrz * uiz;
                final double scip1 = pix * dkx + piy * dky + piz * dkz + dix * pkx + diy * pky + diz * pkz;
                final double scip2 = uix * pkx + uiy * pky + uiz * pkz + pix * ukx + piy * uky + piz * ukz;
                final double scip3 = pix * xr + piy * yr + piz * zr;
                final double scip4 = pkx * xr + pky * yr + pkz * zr;
                final double scip7 = qirx * pkx + qiry * pky + qirz * pkz;
                final double scip8 = qkrx * pix + qkry * piy + qkrz * piz;
                /**
                 * Calculate the gl functions for polarization components.
                 */
                final double gli1 = ck * sci3 - ci * sci4;
                final double gli2 = -sc3 * sci4 - sci3 * sc4;
                final double gli3 = sci3 * sc6 - sci4 * sc5;
                final double gli6 = sci1;
                final double gli7 = 2.0 * (sci7 - sci8);
                final double glip1 = ck * scip3 - ci * scip4;
                final double glip2 = -sc3 * scip4 - scip3 * sc4;
                final double glip3 = scip3 * sc6 - scip4 * sc5;
                final double glip6 = scip1;
                final double glip7 = 2.0 * (scip7 - scip8);
                /**
                 * Compute the energy contributions for this interaction.
                 */
                double ei = (gli1 + gli6) * bn1 + (gli2 + gli7) * bn2 + gli3 * bn3;
                final double eifix = (gli1 + gli6) * rr3 * psc3 + (gli2 + gli7) * rr5 * psc5 + gli3 * rr7 * psc7;
                ei = ei - eifix;
                ei = 0.5 * ei;
                if (!gradient) {
                    return ei;
                }
                boolean dorli = false;
                if (psc3 != 0.0 || dsc3 != 0.0 || usc3 != 0.0) {
                    dorli = true;
                }
                /*
                 * Get the induced force with screening.
                 */
                final double gfi1 = 0.5 * bn2 * (gli1 + glip1 + gli6 + glip6) + 0.5 * bn2 * scip2 + 0.5 * bn3 * (gli2 + glip2 + gli7 + glip7) - 0.5 * bn3 * (sci3 * scip4 + scip3 * sci4) + 0.5 * bn4 * (gli3 + glip3);
                final double gfi2 = -ck * bn1 + sc4 * bn2 - sc6 * bn3;
                final double gfi3 = ci * bn1 + sc3 * bn2 + sc5 * bn3;
                final double gfi4 = 2.0 * bn2;
                final double gfi5 = bn3 * (sci4 + scip4);
                final double gfi6 = -bn3 * (sci3 + scip3);
                double ftm2ix = gfi1 * xr + 0.5 * (gfi2 * (uix + pix) + bn2 * (sci4 * pix + scip4 * uix) + gfi3 * (ukx + pkx) + bn2 * (sci3 * pkx + scip3 * ukx) + (sci4 + scip4) * bn2 * dix + (sci3 + scip3) * bn2 * dkx + gfi4 * (qkuix + qkuipx - qiukx - qiukpx)) + gfi5 * qirx + gfi6 * qkrx;
                double ftm2iy = gfi1 * yr + 0.5 * (gfi2 * (uiy + piy) + bn2 * (sci4 * piy + scip4 * uiy) + gfi3 * (uky + pky) + bn2 * (sci3 * pky + scip3 * uky) + (sci4 + scip4) * bn2 * diy + (sci3 + scip3) * bn2 * dky + gfi4 * (qkuiy + qkuipy - qiuky - qiukpy)) + gfi5 * qiry + gfi6 * qkry;
                double ftm2iz = gfi1 * zr + 0.5 * (gfi2 * (uiz + piz) + bn2 * (sci4 * piz + scip4 * uiz) + gfi3 * (ukz + pkz) + bn2 * (sci3 * pkz + scip3 * ukz) + (sci4 + scip4) * bn2 * diz + (sci3 + scip3) * bn2 * dkz + gfi4 * (qkuiz + qkuipz - qiukz - qiukpz)) + gfi5 * qirz + gfi6 * qkrz;
                /*
                 * Get the induced torque with screening.
                 */
                final double gti2 = 0.5 * bn2 * (sci4 + scip4);
                final double gti3 = 0.5 * bn2 * (sci3 + scip3);
                final double gti4 = gfi4;
                final double gti5 = gfi5;
                final double gti6 = gfi6;
                double ttm2ix = -0.5 * bn1 * (dixukx + dixukpx) + gti2 * dixrx - gti5 * rxqirx + 0.5 * gti4 * (ukxqirx + rxqiukx + ukxqirpx + rxqiukpx);
                double ttm2iy = -0.5 * bn1 * (dixuky + dixukpy) + gti2 * dixry - gti5 * rxqiry + 0.5 * gti4 * (ukxqiry + rxqiuky + ukxqirpy + rxqiukpy);
                double ttm2iz = -0.5 * bn1 * (dixukz + dixukpz) + gti2 * dixrz - gti5 * rxqirz + 0.5 * gti4 * (ukxqirz + rxqiukz + ukxqirpz + rxqiukpz);
                double ttm3ix = -0.5 * bn1 * (dkxuix + dkxuipx) + gti3 * dkxrx - gti6 * rxqkrx - 0.5 * gti4 * (uixqkrx + rxqkuix + uixqkrpx + rxqkuipx);
                double ttm3iy = -0.5 * bn1 * (dkxuiy + dkxuipy) + gti3 * dkxry - gti6 * rxqkry - 0.5 * gti4 * (uixqkry + rxqkuiy + uixqkrpy + rxqkuipy);
                double ttm3iz = -0.5 * bn1 * (dkxuiz + dkxuipz) + gti3 * dkxrz - gti6 * rxqkrz - 0.5 * gti4 * (uixqkrz + rxqkuiz + uixqkrpz + rxqkuipz);
                double ftm2rix = 0.0;
                double ftm2riy = 0.0;
                double ftm2riz = 0.0;
                double ttm2rix = 0.0;
                double ttm2riy = 0.0;
                double ttm2riz = 0.0;
                double ttm3rix = 0.0;
                double ttm3riy = 0.0;
                double ttm3riz = 0.0;
                if (dorli) {
                    /*
                     * Get the induced force without screening.
                     */
                    final double gfri1 = 0.5 * rr5 * ((gli1 + gli6) * psc3 + (glip1 + glip6) * dsc3 + scip2 * usc3) + 0.5 * rr7 * ((gli7 + gli2) * psc5 + (glip7 + glip2) * dsc5 - (sci3 * scip4 + scip3 * sci4) * usc5) + 0.5 * rr9 * (gli3 * psc7 + glip3 * dsc7);
                    final double gfri4 = 2.0 * rr5;
                    final double gfri5 = rr7 * (sci4 * psc7 + scip4 * dsc7);
                    final double gfri6 = -rr7 * (sci3 * psc7 + scip3 * dsc7);
                    ftm2rix = gfri1 * xr + 0.5 * (-rr3 * ck * (uix * psc3 + pix * dsc3) + rr5 * sc4 * (uix * psc5 + pix * dsc5) - rr7 * sc6 * (uix * psc7 + pix * dsc7)) + (rr3 * ci * (ukx * psc3 + pkx * dsc3) + rr5 * sc3 * (ukx * psc5 + pkx * dsc5) + rr7 * sc5 * (ukx * psc7 + pkx * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * pix + scip4 * uix + sci3 * pkx + scip3 * ukx) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * dix + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dkx + 0.5 * gfri4 * ((qkuix - qiukx) * psc5 + (qkuipx - qiukpx) * dsc5) + gfri5 * qirx + gfri6 * qkrx;
                    ftm2riy = gfri1 * yr + 0.5 * (-rr3 * ck * (uiy * psc3 + piy * dsc3) + rr5 * sc4 * (uiy * psc5 + piy * dsc5) - rr7 * sc6 * (uiy * psc7 + piy * dsc7)) + (rr3 * ci * (uky * psc3 + pky * dsc3) + rr5 * sc3 * (uky * psc5 + pky * dsc5) + rr7 * sc5 * (uky * psc7 + pky * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * piy + scip4 * uiy + sci3 * pky + scip3 * uky) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * diy + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dky + 0.5 * gfri4 * ((qkuiy - qiuky) * psc5 + (qkuipy - qiukpy) * dsc5) + gfri5 * qiry + gfri6 * qkry;
                    ftm2riz = gfri1 * zr + 0.5 * (-rr3 * ck * (uiz * psc3 + piz * dsc3) + rr5 * sc4 * (uiz * psc5 + piz * dsc5) - rr7 * sc6 * (uiz * psc7 + piz * dsc7)) + (rr3 * ci * (ukz * psc3 + pkz * dsc3) + rr5 * sc3 * (ukz * psc5 + pkz * dsc5) + rr7 * sc5 * (ukz * psc7 + pkz * dsc7)) * 0.5 + rr5 * usc5 * (sci4 * piz + scip4 * uiz + sci3 * pkz + scip3 * ukz) * 0.5 + 0.5 * (sci4 * psc5 + scip4 * dsc5) * rr5 * diz + 0.5 * (sci3 * psc5 + scip3 * dsc5) * rr5 * dkz + 0.5 * gfri4 * ((qkuiz - qiukz) * psc5 + (qkuipz - qiukpz) * dsc5) + gfri5 * qirz + gfri6 * qkrz;
                    /*
                     * Get the induced torque without screening.
                     */
                    final double gtri2 = 0.5 * rr5 * (sci4 * psc5 + scip4 * dsc5);
                    final double gtri3 = 0.5 * rr5 * (sci3 * psc5 + scip3 * dsc5);
                    final double gtri4 = gfri4;
                    final double gtri5 = gfri5;
                    final double gtri6 = gfri6;
                    ttm2rix = -rr3 * (dixukx * psc3 + dixukpx * dsc3) * 0.5 + gtri2 * dixrx - gtri5 * rxqirx + gtri4 * ((ukxqirx + rxqiukx) * psc5 + (ukxqirpx + rxqiukpx) * dsc5) * 0.5;
                    ttm2riy = -rr3 * (dixuky * psc3 + dixukpy * dsc3) * 0.5 + gtri2 * dixry - gtri5 * rxqiry + gtri4 * ((ukxqiry + rxqiuky) * psc5 + (ukxqirpy + rxqiukpy) * dsc5) * 0.5;
                    ttm2riz = -rr3 * (dixukz * psc3 + dixukpz * dsc3) * 0.5 + gtri2 * dixrz - gtri5 * rxqirz + gtri4 * ((ukxqirz + rxqiukz) * psc5 + (ukxqirpz + rxqiukpz) * dsc5) * 0.5;
                    ttm3rix = -rr3 * (dkxuix * psc3 + dkxuipx * dsc3) * 0.5 + gtri3 * dkxrx - gtri6 * rxqkrx - gtri4 * ((uixqkrx + rxqkuix) * psc5 + (uixqkrpx + rxqkuipx) * dsc5) * 0.5;
                    ttm3riy = -rr3 * (dkxuiy * psc3 + dkxuipy * dsc3) * 0.5 + gtri3 * dkxry - gtri6 * rxqkry - gtri4 * ((uixqkry + rxqkuiy) * psc5 + (uixqkrpy + rxqkuipy) * dsc5) * 0.5;
                    ttm3riz = -rr3 * (dkxuiz * psc3 + dkxuipz * dsc3) * 0.5 + gtri3 * dkxrz - gtri6 * rxqkrz - gtri4 * ((uixqkrz + rxqkuiz) * psc5 + (uixqkrpz + rxqkuipz) * dsc5) * 0.5;
                }
                /*
                 * Account for partially excluded induced interactions.
                 */
                double temp3 = 0.5 * rr3 * ((gli1 + gli6) * scalep + (glip1 + glip6) * scaled);
                double temp5 = 0.5 * rr5 * ((gli2 + gli7) * scalep + (glip2 + glip7) * scaled);
                final double temp7 = 0.5 * rr7 * (gli3 * scalep + glip3 * scaled);
                final double fridmpx = temp3 * ddsc3x + temp5 * ddsc5x + temp7 * ddsc7x;
                final double fridmpy = temp3 * ddsc3y + temp5 * ddsc5y + temp7 * ddsc7y;
                final double fridmpz = temp3 * ddsc3z + temp5 * ddsc5z + temp7 * ddsc7z;
                /*
                 * Find some scaling terms for induced-induced force.
                 */
                temp3 = 0.5 * rr3 * scip2;
                temp5 = -0.5 * rr5 * (sci3 * scip4 + scip3 * sci4);
                final double findmpx = temp3 * ddsc3x + temp5 * ddsc5x;
                final double findmpy = temp3 * ddsc3y + temp5 * ddsc5y;
                final double findmpz = temp3 * ddsc3z + temp5 * ddsc5z;
                /*
                 * Modify the forces for partially excluded interactions.
                 */
                ftm2ix = ftm2ix - fridmpx - findmpx;
                ftm2iy = ftm2iy - fridmpy - findmpy;
                ftm2iz = ftm2iz - fridmpz - findmpz;
                /*
                 * Correction to convert mutual to direct polarization force.
                 */
                if (polarization == Polarization.DIRECT) {
                    final double gfd = 0.5 * (bn2 * scip2 - bn3 * (scip3 * sci4 + sci3 * scip4));
                    final double gfdr = 0.5 * (rr5 * scip2 * usc3 - rr7 * (scip3 * sci4 + sci3 * scip4) * usc5);
                    ftm2ix = ftm2ix - gfd * xr - 0.5 * bn2 * (sci4 * pix + scip4 * uix + sci3 * pkx + scip3 * ukx);
                    ftm2iy = ftm2iy - gfd * yr - 0.5 * bn2 * (sci4 * piy + scip4 * uiy + sci3 * pky + scip3 * uky);
                    ftm2iz = ftm2iz - gfd * zr - 0.5 * bn2 * (sci4 * piz + scip4 * uiz + sci3 * pkz + scip3 * ukz);
                    final double fdirx = gfdr * xr + 0.5 * usc5 * rr5 * (sci4 * pix + scip4 * uix + sci3 * pkx + scip3 * ukx);
                    final double fdiry = gfdr * yr + 0.5 * usc5 * rr5 * (sci4 * piy + scip4 * uiy + sci3 * pky + scip3 * uky);
                    final double fdirz = gfdr * zr + 0.5 * usc5 * rr5 * (sci4 * piz + scip4 * uiz + sci3 * pkz + scip3 * ukz);
                    ftm2ix = ftm2ix + fdirx + findmpx;
                    ftm2iy = ftm2iy + fdiry + findmpy;
                    ftm2iz = ftm2iz + fdirz + findmpz;
                }
                /**
                 * Handle the case where scaling is used.
                 */
                ftm2ix = ftm2ix - ftm2rix;
                ftm2iy = ftm2iy - ftm2riy;
                ftm2iz = ftm2iz - ftm2riz;
                ttm2ix = ttm2ix - ttm2rix;
                ttm2iy = ttm2iy - ttm2riy;
                ttm2iz = ttm2iz - ttm2riz;
                ttm3ix = ttm3ix - ttm3rix;
                ttm3iy = ttm3iy - ttm3riy;
                ttm3iz = ttm3iz - ttm3riz;
                if (iSymm == 0) {
                    gxi_local[i] += ftm2ix;
                    gyi_local[i] += ftm2iy;
                    gzi_local[i] += ftm2iz;
                    gxi_local[k] -= ftm2ix;
                    gyi_local[k] -= ftm2iy;
                    gzi_local[k] -= ftm2iz;
                    txi_local[i] += ttm2ix;
                    tyi_local[i] += ttm2iy;
                    tzi_local[i] += ttm2iz;
                    txi_local[k] += ttm3ix;
                    tyi_local[k] += ttm3iy;
                    tzi_local[k] += ttm3iz;
                } else {
                    gxi_local[i] += 0.5 * ftm2ix;
                    gyi_local[i] += 0.5 * ftm2iy;
                    gzi_local[i] += 0.5 * ftm2iz;
                    txi_local[i] += 0.5 * ttm2ix;
                    tyi_local[i] += 0.5 * ttm2iy;
                    tzi_local[i] += 0.5 * ttm2iz;
                    gxk_local[k] -= 0.5 * ftm2ix;
                    gyk_local[k] -= 0.5 * ftm2iy;
                    gzk_local[k] -= 0.5 * ftm2iz;
                    txk_local[k] += 0.5 * ttm3ix;
                    tyk_local[k] += 0.5 * ttm3iy;
                    tzk_local[k] += 0.5 * ttm3iz;
                }
                return ei;
            }
        }
    }

    private class ExpandCoordinatesRegion extends ParallelRegion {

        private final ExpandCoordinatesLoop expandCoordinatesLoop[];

        public ExpandCoordinatesRegion(int maxThreads) {
            expandCoordinatesLoop = new ExpandCoordinatesLoop[maxThreads];
            for (int i = 0; i < maxThreads; i++) {
                expandCoordinatesLoop[i] = new ExpandCoordinatesLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, expandCoordinatesLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception expanding coordinates in thread: " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class ExpandCoordinatesLoop extends IntegerForLoop {

            private final double in[] = new double[3];
            private final double out[] = new double[3];
            private final double x[] = coordinates[0][0];
            private final double y[] = coordinates[0][1];
            private final double z[] = coordinates[0][2];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public IntegerSchedule schedule() {
                return pairWiseSchedule;
            }

            @Override
            public void run(int lb, int ub) {
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                for (int iSymm = 1; iSymm < nSymm; iSymm++) {
                    SymOp symOp = symOps.get(iSymm);
                    double xs[] = coordinates[iSymm][0];
                    double ys[] = coordinates[iSymm][1];
                    double zs[] = coordinates[iSymm][2];
                    for (int i = lb; i <= ub; i++) {
                        in[0] = x[i];
                        in[1] = y[i];
                        in[2] = z[i];
                        crystal.applySymOp(in, out, symOp);
                        xs[i] = out[0];
                        ys[i] = out[1];
                        zs[i] = out[2];
                    }
                }
            }
        }
    }

    public static class RotateMultipolesRegion extends ParallelRegion {

        private final RotateMultipolesLoop rotateMultipolesLoop[];

        public RotateMultipolesRegion(int nt, boolean inverse) {
            rotateMultipolesLoop = new RotateMultipolesLoop[nt];
            for (int i = 0; i < nt; i++) {
                rotateMultipolesLoop[i] = new RotateMultipolesLoop(inverse);
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, rotateMultipolesLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception rotating multipoles in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        public class RotateMultipolesLoop extends IntegerForLoop {
            // Local variables

            private final double localOrigin[] = new double[3];
            private final double xAxis[] = new double[3];
            private final double yAxis[] = new double[3];
            private final double zAxis[] = new double[3];
            private double rotmat[][] = new double[3][3];
            private final double tempDipole[] = new double[3];
            private final double tempQuadrupole[][] = new double[3][3];
            private final double dipole[] = new double[3];
            private final double quadrupole[][] = new double[3][3];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;
            private boolean inverse = false;

            public RotateMultipolesLoop(boolean inverse) {
                this.inverse = inverse;
            }

            @Override
            public void run(int lb, int ub) {
                for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                    final double x[] = coordinates[iSymm][0];
                    final double y[] = coordinates[iSymm][1];
                    final double z[] = coordinates[iSymm][2];
                    for (int ii = lb; ii <= ub; ii++) {
                        final double in[] = inverse ? globalMultipole[iSymm][ii] : localMultipole[ii];
                        final double out[] = inverse ? localMultipole[ii] : globalMultipole[iSymm][ii];
                        localOrigin[0] = x[ii];
                        localOrigin[1] = y[ii];
                        localOrigin[2] = z[ii];
                        int referenceSites[] = axisAtom[ii];
                        for (int i = 0; i < 3; i++) {
                            zAxis[i] = 0.0;
                            xAxis[i] = 0.0;
                            dipole[i] = 0.0;
                            for (int j = 0; j < 3; j++) {
                                quadrupole[i][j] = 0.0;
                            }
                        }
                        if (referenceSites == null || referenceSites.length < 2) {
                            out[t000] = in[0];
                            out[t100] = 0.0;
                            out[t010] = 0.0;
                            out[t001] = 0.0;
                            out[t200] = 0.0;
                            out[t020] = 0.0;
                            out[t002] = 0.0;
                            out[t110] = 0.0;
                            out[t101] = 0.0;
                            out[t011] = 0.0;
                            continue;
                        }
                        switch (frame[ii]) {
                            case BISECTOR:
                                int index = referenceSites[0];
                                zAxis[0] = x[index];
                                zAxis[1] = y[index];
                                zAxis[2] = z[index];
                                index = referenceSites[1];
                                xAxis[0] = x[index];
                                xAxis[1] = y[index];
                                xAxis[2] = z[index];
                                diff(zAxis, localOrigin, zAxis);
                                norm(zAxis, zAxis);
                                diff(xAxis, localOrigin, xAxis);
                                norm(xAxis, xAxis);
                                sum(xAxis, zAxis, zAxis);
                                norm(zAxis, zAxis);
                                rotmat[0][2] = zAxis[0];
                                rotmat[1][2] = zAxis[1];
                                rotmat[2][2] = zAxis[2];
                                double dot = dot(xAxis, zAxis);
                                scalar(zAxis, dot, zAxis);
                                diff(xAxis, zAxis, xAxis);
                                norm(xAxis, xAxis);
                                rotmat[0][0] = xAxis[0];
                                rotmat[1][0] = xAxis[1];
                                rotmat[2][0] = xAxis[2];
                                break;
                            case ZTHENBISECTOR:
                                index = referenceSites[0];
                                zAxis[0] = x[index];
                                zAxis[1] = y[index];
                                zAxis[2] = z[index];
                                index = referenceSites[1];
                                xAxis[0] = x[index];
                                xAxis[1] = y[index];
                                xAxis[2] = z[index];
                                index = referenceSites[2];
                                yAxis[0] = x[index];
                                yAxis[1] = y[index];
                                yAxis[2] = z[index];
                                diff(zAxis, localOrigin, zAxis);
                                norm(zAxis, zAxis);
                                rotmat[0][2] = zAxis[0];
                                rotmat[1][2] = zAxis[1];
                                rotmat[2][2] = zAxis[2];
                                diff(xAxis, localOrigin, xAxis);
                                norm(xAxis, xAxis);
                                diff(yAxis, localOrigin, yAxis);
                                norm(yAxis, yAxis);
                                sum(xAxis, yAxis, xAxis);
                                norm(xAxis, xAxis);
                                dot = dot(xAxis, zAxis);
                                scalar(zAxis, dot, zAxis);
                                diff(xAxis, zAxis, xAxis);
                                norm(xAxis, xAxis);
                                rotmat[0][0] = xAxis[0];
                                rotmat[1][0] = xAxis[1];
                                rotmat[2][0] = xAxis[2];
                                break;
                            case ZTHENX:
                            default:
                                index = referenceSites[0];
                                zAxis[0] = x[index];
                                zAxis[1] = y[index];
                                zAxis[2] = z[index];
                                index = referenceSites[1];
                                if (index != -1) {
                                    xAxis[0] = x[index];
                                    xAxis[1] = y[index];
                                    xAxis[2] = z[index];
                                } else if (index == -1) {
                                    xAxis[0] = 0;
                                    xAxis[1] = 0;
                                    xAxis[2] = 0;
                                }
                                diff(zAxis, localOrigin, zAxis);
                                norm(zAxis, zAxis);
                                rotmat[0][2] = zAxis[0];
                                rotmat[1][2] = zAxis[1];
                                rotmat[2][2] = zAxis[2];
                                diff(xAxis, localOrigin, xAxis);
                                dot = dot(xAxis, zAxis);
                                scalar(zAxis, dot, zAxis);
                                diff(xAxis, zAxis, xAxis);
                                norm(xAxis, xAxis);
                                rotmat[0][0] = xAxis[0];
                                rotmat[1][0] = xAxis[1];
                                rotmat[2][0] = xAxis[2];
                        }
                        // Finally the Y elements.
                        rotmat[0][1] = rotmat[2][0] * rotmat[1][2] - rotmat[1][0] * rotmat[2][2];
                        rotmat[1][1] = rotmat[0][0] * rotmat[2][2] - rotmat[2][0] * rotmat[0][2];
                        rotmat[2][1] = rotmat[1][0] * rotmat[0][2] - rotmat[0][0] * rotmat[1][2];
                        // Do the rotation.
                        tempDipole[0] = in[t100];
                        tempDipole[1] = in[t010];
                        tempDipole[2] = in[t001];
                        if (pedit) {
                            tempQuadrupole[0][0] = in[4];
                            tempQuadrupole[1][1] = in[8];
                            tempQuadrupole[2][2] = in[12];
                            tempQuadrupole[0][1] = in[5];
                            tempQuadrupole[0][2] = in[6];
                            tempQuadrupole[1][2] = in[9];
                            tempQuadrupole[1][0] = in[7];
                            tempQuadrupole[2][0] = in[10];
                            tempQuadrupole[2][1] = in[11];
                        } else {
                            tempQuadrupole[0][0] = in[t200];
                            tempQuadrupole[1][1] = in[t020];
                            tempQuadrupole[2][2] = in[t002];
                            tempQuadrupole[0][1] = in[t110];
                            tempQuadrupole[0][2] = in[t101];
                            tempQuadrupole[1][2] = in[t011];
                            tempQuadrupole[1][0] = in[t110];
                            tempQuadrupole[2][0] = in[t101];
                            tempQuadrupole[2][1] = in[t011];
                        }

                        // Check for chiral flipping.
                        if (frame[ii] == MultipoleType.MultipoleFrameDefinition.ZTHENX
                                && referenceSites.length == 3 && !pedit) {
                            localOrigin[0] = x[ii];
                            localOrigin[1] = y[ii];
                            localOrigin[2] = z[ii];
                            int index = referenceSites[0];
                            zAxis[0] = x[index];
                            zAxis[1] = y[index];
                            zAxis[2] = z[index];
                            index = referenceSites[1];
                            xAxis[0] = x[index];
                            xAxis[1] = y[index];
                            xAxis[2] = z[index];
                            index = referenceSites[2];
                            yAxis[0] = x[index];
                            yAxis[1] = y[index];
                            yAxis[2] = z[index];
                            diff(localOrigin, yAxis, localOrigin);
                            diff(zAxis, yAxis, zAxis);
                            diff(xAxis, yAxis, xAxis);
                            double c1 = zAxis[1] * xAxis[2] - zAxis[2] * xAxis[1];
                            double c2 = xAxis[1] * localOrigin[2] - xAxis[2] * localOrigin[1];
                            double c3 = localOrigin[1] * zAxis[2] - localOrigin[2] * zAxis[1];
                            double vol = localOrigin[0] * c1 + zAxis[0] * c2 + xAxis[0] * c3;
                            if (vol < 0.0) {
                                tempDipole[1] = -tempDipole[1];
                                tempQuadrupole[0][1] = -tempQuadrupole[0][1];
                                tempQuadrupole[1][0] = -tempQuadrupole[1][0];
                                tempQuadrupole[1][2] = -tempQuadrupole[1][2];
                                tempQuadrupole[2][1] = -tempQuadrupole[2][1];
                                //logger.info("Chiral flip of atom " + atoms[ii]);
                            }
                        }
                        if (inverse) {
                            rotmat = transpose3(rotmat);
                        }
                        for (int i = 0; i < 3; i++) {
                            double[] rotmati = rotmat[i];
                            double[] quadrupolei = quadrupole[i];
                            for (int j = 0; j < 3; j++) {
                                double[] rotmatj = rotmat[j];
                                dipole[i] += rotmati[j] * tempDipole[j];
                                if (j < i) {
                                    quadrupolei[j] = quadrupole[j][i];
                                } else {
                                    for (int k = 0; k < 3; k++) {
                                        double[] localQuadrupolek = tempQuadrupole[k];
                                        quadrupolei[j] += rotmati[k]
                                                * (rotmatj[0] * localQuadrupolek[0]
                                                + rotmatj[1] * localQuadrupolek[1]
                                                + rotmatj[2] * localQuadrupolek[2]);
                                    }
                                }
                            }
                        }
                        if (frame[ii] != MultipoleType.MultipoleFrameDefinition.ZTHENX
                                && referenceSites.length == 3 && pedit && axisAtom[ii][2] != 0) {
                            localOrigin[0] = x[ii];
                            localOrigin[1] = y[ii];
                            localOrigin[2] = z[ii];
                            int index = referenceSites[0];
                            zAxis[0] = x[index];
                            zAxis[1] = y[index];
                            zAxis[2] = z[index];
                            index = referenceSites[1];
                            xAxis[0] = x[index];
                            xAxis[1] = y[index];
                            xAxis[2] = z[index];
                            index = referenceSites[2];
                            yAxis[0] = x[index];
                            yAxis[1] = y[index];
                            yAxis[2] = z[index];
                            diff(localOrigin, yAxis, localOrigin);
                            diff(zAxis, yAxis, zAxis);
                            diff(xAxis, yAxis, xAxis);
                            double c1 = zAxis[1] * xAxis[2] - zAxis[2] * xAxis[1];
                            double c2 = xAxis[1] * localOrigin[2] - xAxis[2] * localOrigin[1];
                            double c3 = localOrigin[1] * zAxis[2] - localOrigin[2] * zAxis[1];
                            double vol = localOrigin[0] * c1 + zAxis[0] * c2 + xAxis[0] * c3;
                            if ((axisAtom[ii][2] < 0 && vol > 0) || (axisAtom[ii][2] > 0 && vol < 0.0)) {
                                axisAtom[ii][2] = -axisAtom[ii][2];
                                dipole[1] = -dipole[1];
                                quadrupole[0][1] = -quadrupole[0][1];
                                quadrupole[1][0] = -quadrupole[1][0];
                                quadrupole[1][2] = -quadrupole[1][2];
                                quadrupole[2][1] = -quadrupole[2][1];
                                //logger.info("Chiral flip of atom " + atoms[ii]);
                            }
                        }
                        /**
                         * R.Q.R^t is correct. RealMatrix quad = new
                         * Array2DRowRealMatrix(tempQuadrupole); RealMatrix rot
                         * = new Array2DRowRealMatrix(rotmat); quad =
                         * rot.multiply(quad).multiply(rot.transpose());
                         * logger.info(format("%d %8.3f %8.3f", ii,
                         * quad.getEntry(0,0), quadrupole[0][0]));
                         * logger.info(format("%d %8.3f %8.3f\n", ii,
                         * quad.getEntry(0,1), quadrupole[0][1]));
                         */
                        double scale = 1.0;
                        Atom a = atoms[ii];
                        if (a.applyLambda()) {
                            scale = lambda;
                        }
                        double dipoleScale = 1.0;
                        double quadrupoleScale = 1.0;
                        out[t000] = scale * in[0];
                        out[t100] = scale * dipole[0] * dipoleScale;
                        out[t010] = scale * dipole[1] * dipoleScale;
                        out[t001] = scale * dipole[2] * dipoleScale;
                        if (pedit) {
                            out[4] = scale * quadrupole[0][0] * quadrupoleScale;
                            out[8] = scale * quadrupole[1][1] * quadrupoleScale;
                            out[12] = scale * quadrupole[2][2] * quadrupoleScale;
                            out[7] = scale * quadrupole[0][1] * quadrupoleScale;
                            out[10] = scale * quadrupole[0][2] * quadrupoleScale;
                            out[11] = scale * quadrupole[1][2] * quadrupoleScale;
                        } else {
                            out[t200] = scale * quadrupole[0][0] * quadrupoleScale;
                            out[t020] = scale * quadrupole[1][1] * quadrupoleScale;
                            out[t002] = scale * quadrupole[2][2] * quadrupoleScale;
                            out[t110] = scale * quadrupole[0][1] * quadrupoleScale;
                            out[t101] = scale * quadrupole[0][2] * quadrupoleScale;
                            out[t011] = scale * quadrupole[1][2] * quadrupoleScale;
                        }

                        /* No rotation (useful for checking non-torque parts
                         * of the multipole gradient
                         out[t000] = scale * in[t000];
                         out[t100] = scale * in[t100] * dipoleScale;
                         out[t010] = scale * in[t010] * dipoleScale;
                         out[t001] = scale * in[t001] * dipoleScale;
                         out[t200] = scale * in[t200] * quadrupoleScale;
                         out[t020] = scale * in[t020] * quadrupoleScale;
                         out[t002] = scale * in[t002] * quadrupoleScale;
                         out[t110] = scale * in[t110] * quadrupoleScale;
                         out[t101] = scale * in[t101] * quadrupoleScale;
                         out[t011] = scale * in[t011] * quadrupoleScale;
                         */
                        if (!pedit) {
                            PolarizeType polarizeType = a.getPolarizeType();
                            polarizability[ii] = scale * polarizeType.polarizability;
                        }
                    }
                }
            }
        }
    }

    public static class ExpandInducedDipolesRegion extends ParallelRegion {

        private final ExpandInducedDipoleLoop expandInducedDipoleLoop[];

        public ExpandInducedDipolesRegion(int maxThreads) {
            expandInducedDipoleLoop = new ExpandInducedDipoleLoop[maxThreads];
            for (int i = 0; i < maxThreads; i++) {
                expandInducedDipoleLoop[i] = new ExpandInducedDipoleLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, expandInducedDipoleLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception expanding coordinates in thread: " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        public class ExpandInducedDipoleLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return pairWiseSchedule;
            }

            @Override
            public void run(int lb, int ub) {
                for (int s = 1; s < nSymm; s++) {
                    SymOp symOp = crystal.spaceGroup.symOps.get(s);
                    for (int ii = lb; ii <= ub; ii++) {
                        crystal.applySymRot(inducedDipole[0][ii], inducedDipole[s][ii], symOp);
                        crystal.applySymRot(inducedDipolep[0][ii], inducedDipolep[s][ii], symOp);
                    }
                }
            }
        }
    }

    private class TorqueRegion extends ParallelRegion {

        private final TorqueLoop torqueLoop[];

        public TorqueRegion(int maxThreads) {
            torqueLoop = new TorqueLoop[maxThreads];
            for (int i = 0; i < maxThreads; i++) {
                torqueLoop[i] = new TorqueLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, torqueLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing torque in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class TorqueLoop extends IntegerForLoop {

            private final double trq[] = new double[3];
            private final double u[] = new double[3];
            private final double v[] = new double[3];
            private final double w[] = new double[3];
            private final double r[] = new double[3];
            private final double s[] = new double[3];
            private final double uv[] = new double[3];
            private final double uw[] = new double[3];
            private final double vw[] = new double[3];
            private final double ur[] = new double[3];
            private final double us[] = new double[3];
            private final double vs[] = new double[3];
            private final double ws[] = new double[3];
            private final double t1[] = new double[3];
            private final double t2[] = new double[3];
            private final double localOrigin[] = new double[3];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    final int ax[] = axisAtom[i];
                    // Ions, for example, have no torque.
                    if (ax == null || ax.length < 2) {
                        continue;
                    }
                    final int ia = ax[0];
                    final int ib = i;
                    final int ic = ax[1];
                    int id = 0;
                    trq[0] = sharedTorque[0].get(i);
                    trq[1] = sharedTorque[1].get(i);
                    trq[2] = sharedTorque[2].get(i);
                    double x[] = coordinates[0][0];
                    double y[] = coordinates[0][1];
                    double z[] = coordinates[0][2];
                    localOrigin[0] = x[ib];
                    localOrigin[1] = y[ib];
                    localOrigin[2] = z[ib];
                    u[0] = x[ia];
                    u[1] = y[ia];
                    u[2] = z[ia];
                    v[0] = x[ic];
                    v[1] = y[ic];
                    v[2] = z[ic];
                    // Construct the three rotation axes for the local frame
                    diff(u, localOrigin, u);
                    diff(v, localOrigin, v);
                    switch (frame[i]) {
                        default:
                        case ZTHENX:
                        case BISECTOR:
                            cross(u, v, w);
                            break;
                        case TRISECTOR:
                        case ZTHENBISECTOR:
                            id = ax[2];
                            w[0] = x[id];
                            w[1] = y[id];
                            w[2] = z[id];
                            diff(w, localOrigin, w);
                    }

                    double ru = r(u);
                    double rv = r(v);
                    double rw = r(w);
                    scalar(u, 1.0 / ru, u);
                    scalar(v, 1.0 / rv, v);
                    scalar(w, 1.0 / rw, w);
                    // Find the perpendicular and angle for each pair of axes.
                    cross(v, u, uv);
                    cross(w, u, uw);
                    cross(w, v, vw);
                    double ruv = r(uv);
                    double ruw = r(uw);
                    double rvw = r(vw);
                    scalar(uv, 1.0 / ruv, uv);
                    scalar(uw, 1.0 / ruw, uw);
                    scalar(vw, 1.0 / rvw, vw);
                    // Compute the sine of the angle between the rotation axes.
                    double uvcos = dot(u, v);
                    double uvsin = sqrt(1.0 - uvcos * uvcos);
                    //double uwcos = dot(u, w);
                    //double uwsin = sqrt(1.0 - uwcos * uwcos);
                    //double vwcos = dot(v, w);
                    //double vwsin = sqrt(1.0 - vwcos * vwcos);
                    /*
                     * Negative of dot product of torque with unit vectors gives result of
                     * infinitesimal rotation along these vectors.
                     */
                    double dphidu = -(trq[0] * u[0] + trq[1] * u[1] + trq[2] * u[2]);
                    double dphidv = -(trq[0] * v[0] + trq[1] * v[1] + trq[2] * v[2]);
                    double dphidw = -(trq[0] * w[0] + trq[1] * w[1] + trq[2] * w[2]);
                    switch (frame[i]) {
                        case ZTHENBISECTOR:
                            // Build some additional axes needed for the Z-then-Bisector method
                            sum(v, w, r);
                            cross(u, r, s);
                            double rr = r(r);
                            double rs = r(s);
                            scalar(r, 1.0 / rr, r);
                            scalar(s, 1.0 / rs, s);
                            // Find the perpendicular and angle for each pair of axes.
                            cross(r, u, ur);
                            cross(s, u, us);
                            cross(s, v, vs);
                            cross(s, w, ws);
                            double rur = r(ur);
                            double rus = r(us);
                            double rvs = r(vs);
                            double rws = r(ws);
                            scalar(ur, 1.0 / rur, ur);
                            scalar(us, 1.0 / rus, us);
                            scalar(vs, 1.0 / rvs, vs);
                            scalar(ws, 1.0 / rws, ws);
                            // Compute the sine of the angle between the rotation axes
                            double urcos = dot(u, r);
                            double ursin = sqrt(1.0 - urcos * urcos);
                            //double uscos = dot(u, s);
                            //double ussin = sqrt(1.0 - uscos * uscos);
                            double vscos = dot(v, s);
                            double vssin = sqrt(1.0 - vscos * vscos);
                            double wscos = dot(w, s);
                            double wssin = sqrt(1.0 - wscos * wscos);
                            // Compute the projection of v and w onto the ru-plane
                            scalar(s, -vscos, t1);
                            scalar(s, -wscos, t2);
                            sum(v, t1, t1);
                            sum(w, t2, t2);
                            double rt1 = r(t1);
                            double rt2 = r(t2);
                            scalar(t1, 1.0 / rt1, t1);
                            scalar(t2, 1.0 / rt2, t2);
                            double ut1cos = dot(u, t1);
                            double ut1sin = sqrt(1.0 - ut1cos * ut1cos);
                            double ut2cos = dot(u, t2);
                            double ut2sin = sqrt(1.0 - ut2cos * ut2cos);
                            double dphidr = -(trq[0] * r[0] + trq[1] * r[1] + trq[2] * r[2]);
                            double dphids = -(trq[0] * s[0] + trq[1] * s[1] + trq[2] * s[2]);
                            for (int j = 0; j < 3; j++) {
                                double du = ur[j] * dphidr / (ru * ursin) + us[j] * dphids / ru;
                                double dv = (vssin * s[j] - vscos * t1[j]) * dphidu / (rv * (ut1sin + ut2sin));
                                double dw = (wssin * s[j] - wscos * t2[j]) * dphidu / (rw * (ut1sin + ut2sin));
                                sharedGrad[j].addAndGet(ia, du);
                                sharedGrad[j].addAndGet(ic, dv);
                                sharedGrad[j].addAndGet(id, dw);
                                sharedGrad[j].addAndGet(ib, -du - dv - dw);
                            }
                            break;
                        case ZTHENX:
                            for (int j = 0; j < 3; j++) {
                                double du = uv[j] * dphidv / (ru * uvsin) + uw[j] * dphidw / ru;
                                double dv = -uv[j] * dphidu / (rv * uvsin);
                                sharedGrad[j].addAndGet(ia, du);
                                sharedGrad[j].addAndGet(ic, dv);
                                sharedGrad[j].addAndGet(ib, -du - dv);
                            }
                            break;
                        case BISECTOR:
                            for (int j = 0; j < 3; j++) {
                                double du = uv[j] * dphidv / (ru * uvsin) + 0.5 * uw[j] * dphidw / ru;
                                double dv = -uv[j] * dphidu / (rv * uvsin) + 0.5 * vw[j] * dphidw / rv;
                                sharedGrad[j].addAndGet(ia, du);
                                sharedGrad[j].addAndGet(ic, dv);
                                sharedGrad[j].addAndGet(ib, -du - dv);
                            }
                            break;
                        default:
                            String message = "Fatal exception: Unknown frame definition: " + frame[i] + "\n";
                            logger.log(Level.SEVERE, message);
                    }
                }
            }
        }
    }

    private double ewaldCoefficient(double cutoff) {
        /*
         * Set the tolerance value; use of 1.0d-8 results in large Ewald
         * coefficients that ensure continuity in the gradient
         */
        double eps = 1.0e-8;
        /*
         * Get an approximate value from cutoff and tolerance.
         */
        double ratio = eps + 1.0;
        double x = 0.5;
        int i = 0;
        // Larger values lead to a more "delta-function-like" Gaussian
        while (ratio >= eps) {
            i++;
            x *= 2.0;
            ratio = erfc(x * cutoff) / cutoff;
        }
        /*
         * Use a binary search to refine the coefficient.
         */
        int k = i + 60;
        double xlo = 0.0;
        double xhi = x;
        for (int j = 0; j < k; j++) {
            x = (xlo + xhi) / 2.0;
            ratio = erfc(x * cutoff) / cutoff;
            if (ratio >= eps) {
                xlo = x;
            } else {
                xhi = x;
            }
        }
        return x;
    }

    /**
     * <p>
     * ewaldCutoff</p>
     *
     * @param coeff a double.
     * @param maxCutoff a double.
     * @param eps a double.
     * @return a double.
     */
    public static double ewaldCutoff(double coeff, double maxCutoff, double eps) {
        /*
         * Set the tolerance value; use of 1.0d-8 requires strict convergence
         * of the real Space sum.
         */
        double ratio = erfc(coeff * maxCutoff) / maxCutoff;

        if (ratio > eps) {
            return maxCutoff;
        }

        /*
         * Use a binary search to refine the coefficient.
         */
        double xlo = 0.0;
        double xhi = maxCutoff;
        double cutoff = 0.0;
        for (int j = 0; j < 100; j++) {
            cutoff = (xlo + xhi) / 2.0;
            ratio = erfc(coeff * cutoff) / cutoff;
            if (ratio >= eps) {
                xlo = cutoff;
            } else {
                xhi = cutoff;
            }
        }
        return cutoff;
    }

    /**
     * Given an array of atoms (with atom types), assign multipole types and
     * reference sites.
     *
     * @param atoms List
     * @param forceField ForceField
     */
    private void assignMultipoles() {
        if (forceField == null) {
            String message = "No force field is defined.\n";
            logger.log(Level.SEVERE, message);
        }
        if (forceField.getForceFieldTypeCount(ForceFieldType.MULTIPOLE) < 1) {
            String message = "Force field has no multipole types.\n";
            logger.log(Level.SEVERE, message);
            return;
        }
        if (nAtoms < 1) {
            String message = "No atoms are defined.\n";
            logger.log(Level.SEVERE, message);
            return;
        }
        for (int i = 0; i < nAtoms; i++) {
            if (!assignMultipole(i)) {
                Atom atom = atoms[i];
                String message = "No multipole could be assigned to atom:\n"
                        + atom + "\nof type:\n" + atom.getAtomType();
                logger.log(Level.SEVERE, message);
            }
        }
        /**
         * Check for multipoles that were not assigned correctly.
         */
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < nAtoms; i++) {
            boolean flag = false;
            for (int j = 0; j < 10; j++) {
                if (Double.isNaN(localMultipole[i][j])) {
                    flag = true;
                    break;
                }
            }
            if (flag) {
                sb.append("\n" + atoms[i].toString() + "\n");
                sb.append(format("%d", i + 1));
                for (int j = 0; j < 10; j++) {
                    sb.append(format(" %8.3f", localMultipole[i][j]));
                }
                sb.append("\n");
            }
        }
        if (sb.length() > 0) {
            String message = "Fatal exception: Error assigning multipoles. " + sb.toString();
            logger.log(Level.SEVERE, message);
            System.exit(-1);
        }
    }

    private boolean assignMultipole(int i) {
        Atom atom = atoms[i];
        AtomType atomType = atoms[i].getAtomType();
        if (atomType == null) {
            String message = "Fatal exception: Multipoles can only be assigned to atoms that have been typed.";
            logger.severe(message);
        }
        PolarizeType polarizeType = forceField.getPolarizeType(atomType.getKey());
        if (polarizeType != null) {
            atom.setPolarizeType(polarizeType);
        } else {
            String message = "Fatal Exception: No polarization type was found for " + atom.toString();
            logger.severe(message);
        }
        MultipoleType multipoleType = null;
        String key = null;
        // No reference atoms.
        key = atomType.getKey() + " 0 0";
        multipoleType = forceField.getMultipoleType(key);
        if (multipoleType != null) {
            atom.setMultipoleType(multipoleType, null);
            localMultipole[i][t000] = multipoleType.charge;
            localMultipole[i][t100] = multipoleType.dipole[0];
            localMultipole[i][t010] = multipoleType.dipole[1];
            localMultipole[i][t001] = multipoleType.dipole[2];
            localMultipole[i][t200] = multipoleType.quadrupole[0][0];
            localMultipole[i][t020] = multipoleType.quadrupole[1][1];
            localMultipole[i][t002] = multipoleType.quadrupole[2][2];
            localMultipole[i][t110] = multipoleType.quadrupole[0][1];
            localMultipole[i][t101] = multipoleType.quadrupole[0][2];
            localMultipole[i][t011] = multipoleType.quadrupole[1][2];
            axisAtom[i] = null;
            frame[i] = multipoleType.frameDefinition;
            return true;
        }
        // No bonds.
        List<Bond> bonds = atom.getFFXBonds();
        if (bonds == null || bonds.size() < 1) {
            String message = "Multipoles can only be assigned after bonded relationships are defined.\n";
            logger.severe(message);
        }
        // 1 reference atom.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            key = atomType.getKey() + " " + atom2.getAtomType().getKey() + " 0";
            multipoleType = forceField.getMultipoleType(key);
            if (multipoleType != null) {
                int multipoleReferenceAtoms[] = new int[1];
                multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                atom.setMultipoleType(multipoleType, null);
                localMultipole[i][0] = multipoleType.charge;
                localMultipole[i][1] = multipoleType.dipole[0];
                localMultipole[i][2] = multipoleType.dipole[1];
                localMultipole[i][3] = multipoleType.dipole[2];
                localMultipole[i][4] = multipoleType.quadrupole[0][0];
                localMultipole[i][5] = multipoleType.quadrupole[1][1];
                localMultipole[i][6] = multipoleType.quadrupole[2][2];
                localMultipole[i][7] = multipoleType.quadrupole[0][1];
                localMultipole[i][8] = multipoleType.quadrupole[0][2];
                localMultipole[i][9] = multipoleType.quadrupole[1][2];
                axisAtom[i] = multipoleReferenceAtoms;
                frame[i] = multipoleType.frameDefinition;
                return true;
            }
        }
        // 2 reference atoms.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                key = atomType.getKey() + " " + key2 + " " + key3;
                multipoleType = forceField.getMultipoleType(key);
                if (multipoleType != null) {
                    int multipoleReferenceAtoms[] = new int[2];
                    multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                    multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                    atom.setMultipoleType(multipoleType, null);
                    localMultipole[i][0] = multipoleType.charge;
                    localMultipole[i][1] = multipoleType.dipole[0];
                    localMultipole[i][2] = multipoleType.dipole[1];
                    localMultipole[i][3] = multipoleType.dipole[2];
                    localMultipole[i][4] = multipoleType.quadrupole[0][0];
                    localMultipole[i][5] = multipoleType.quadrupole[1][1];
                    localMultipole[i][6] = multipoleType.quadrupole[2][2];
                    localMultipole[i][7] = multipoleType.quadrupole[0][1];
                    localMultipole[i][8] = multipoleType.quadrupole[0][2];
                    localMultipole[i][9] = multipoleType.quadrupole[1][2];
                    axisAtom[i] = multipoleReferenceAtoms;
                    frame[i] = multipoleType.frameDefinition;
                    return true;
                }
            }
        }
        // 3 reference atoms.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                for (Bond b3 : bonds) {
                    if (b == b3 || b2 == b3) {
                        continue;
                    }
                    Atom atom4 = b3.get1_2(atom);
                    String key4 = atom4.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        int multipoleReferenceAtoms[] = new int[3];
                        multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                        multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                        multipoleReferenceAtoms[2] = atom4.xyzIndex - 1;
                        atom.setMultipoleType(multipoleType, null);
                        localMultipole[i][0] = multipoleType.charge;
                        localMultipole[i][1] = multipoleType.dipole[0];
                        localMultipole[i][2] = multipoleType.dipole[1];
                        localMultipole[i][3] = multipoleType.dipole[2];
                        localMultipole[i][4] = multipoleType.quadrupole[0][0];
                        localMultipole[i][5] = multipoleType.quadrupole[1][1];
                        localMultipole[i][6] = multipoleType.quadrupole[2][2];
                        localMultipole[i][7] = multipoleType.quadrupole[0][1];
                        localMultipole[i][8] = multipoleType.quadrupole[0][2];
                        localMultipole[i][9] = multipoleType.quadrupole[1][2];
                        axisAtom[i] = multipoleReferenceAtoms;
                        frame[i] = multipoleType.frameDefinition;
                        return true;
                    }
                }
                List<Angle> angles = atom.getAngles();
                for (Angle angle : angles) {
                    Atom atom4 = angle.get1_3(atom);
                    if (atom4 != null) {
                        String key4 = atom4.getAtomType().getKey();
                        key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                        multipoleType = forceField.getMultipoleType(key);
                        if (multipoleType != null) {
                            int multipoleReferenceAtoms[] = new int[3];
                            multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                            multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                            multipoleReferenceAtoms[2] = atom4.xyzIndex - 1;
                            atom.setMultipoleType(multipoleType, null);
                            localMultipole[i][0] = multipoleType.charge;
                            localMultipole[i][1] = multipoleType.dipole[0];
                            localMultipole[i][2] = multipoleType.dipole[1];
                            localMultipole[i][3] = multipoleType.dipole[2];
                            localMultipole[i][4] = multipoleType.quadrupole[0][0];
                            localMultipole[i][5] = multipoleType.quadrupole[1][1];
                            localMultipole[i][6] = multipoleType.quadrupole[2][2];
                            localMultipole[i][7] = multipoleType.quadrupole[0][1];
                            localMultipole[i][8] = multipoleType.quadrupole[0][2];
                            localMultipole[i][9] = multipoleType.quadrupole[1][2];
                            axisAtom[i] = multipoleReferenceAtoms;
                            frame[i] = multipoleType.frameDefinition;
                            return true;
                        }
                    }
                }
            }
        }
        // Revert to a 2 reference atom definition that may include a 1-3 site.
        // For example a hydrogen on water.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            List<Angle> angles = atom.getAngles();
            for (Angle angle : angles) {
                Atom atom3 = angle.get1_3(atom);
                if (atom3 != null) {
                    String key3 = atom3.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        int multipoleReferenceAtoms[] = new int[2];
                        multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                        multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                        atom.setMultipoleType(multipoleType, null);
                        localMultipole[i][0] = multipoleType.charge;
                        localMultipole[i][1] = multipoleType.dipole[0];
                        localMultipole[i][2] = multipoleType.dipole[1];
                        localMultipole[i][3] = multipoleType.dipole[2];
                        localMultipole[i][4] = multipoleType.quadrupole[0][0];
                        localMultipole[i][5] = multipoleType.quadrupole[1][1];
                        localMultipole[i][6] = multipoleType.quadrupole[2][2];
                        localMultipole[i][7] = multipoleType.quadrupole[0][1];
                        localMultipole[i][8] = multipoleType.quadrupole[0][2];
                        localMultipole[i][9] = multipoleType.quadrupole[1][2];
                        axisAtom[i] = multipoleReferenceAtoms;
                        frame[i] = multipoleType.frameDefinition;
                        return true;
                    }
                    for (Angle angle2 : angles) {
                        Atom atom4 = angle2.get1_3(atom);
                        if (atom4 != null && atom4 != atom3) {
                            String key4 = atom4.getAtomType().getKey();
                            key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                            multipoleType = forceField.getMultipoleType(key);
                            if (multipoleType != null) {
                                int multipoleReferenceAtoms[] = new int[3];
                                multipoleReferenceAtoms[0] = atom2.xyzIndex - 1;
                                multipoleReferenceAtoms[1] = atom3.xyzIndex - 1;
                                multipoleReferenceAtoms[2] = atom4.xyzIndex - 1;
                                atom.setMultipoleType(multipoleType, null);
                                localMultipole[i][0] = multipoleType.charge;
                                localMultipole[i][1] = multipoleType.dipole[0];
                                localMultipole[i][2] = multipoleType.dipole[1];
                                localMultipole[i][3] = multipoleType.dipole[2];
                                localMultipole[i][4] = multipoleType.quadrupole[0][0];
                                localMultipole[i][5] = multipoleType.quadrupole[1][1];
                                localMultipole[i][6] = multipoleType.quadrupole[2][2];
                                localMultipole[i][7] = multipoleType.quadrupole[0][1];
                                localMultipole[i][8] = multipoleType.quadrupole[0][2];
                                localMultipole[i][9] = multipoleType.quadrupole[1][2];
                                axisAtom[i] = multipoleReferenceAtoms;
                                frame[i] = multipoleType.frameDefinition;
                                return true;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }

    private void assignPolarizationGroups() {
        /**
         * Find directly connected group members for each atom.
         */
        List<Integer> group = new ArrayList<Integer>();
        List<Integer> polarizationGroup = new ArrayList<Integer>();
        //int g11 = 0;
        for (Atom ai : atoms) {
            group.clear();
            polarizationGroup.clear();
            Integer index = ai.getXYZIndex() - 1;
            group.add(index);
            polarizationGroup.add(ai.getType());
            PolarizeType polarizeType = ai.getPolarizeType();
            if (polarizeType != null) {
                if (polarizeType.polarizationGroup != null) {
                    for (int i : polarizeType.polarizationGroup) {
                        if (!polarizationGroup.contains(i)) {
                            polarizationGroup.add(i);
                        }
                    }
                    growGroup(polarizationGroup, group, ai);
                    Collections.sort(group);
                    ip11[index] = new int[group.size()];
                    int j = 0;
                    for (int k : group) {
                        ip11[index][j++] = k;
                    }
                } else {
                    ip11[index] = new int[group.size()];
                    int j = 0;
                    for (int k : group) {
                        ip11[index][j++] = k;
                    }
                }
                //g11 += ip11[index].length;
                //System.out.println(format("%d %d", index + 1, g11));
            } else {
                String message = "The polarize keyword was not found for atom "
                        + (index + 1) + " with type " + ai.getType();
                logger.severe(message);
            }
        }
        /**
         * Find 1-2 group relationships.
         */
        int mask[] = new int[nAtoms];
        List<Integer> list = new ArrayList<Integer>();
        List<Integer> keep = new ArrayList<Integer>();
        for (int i = 0; i < nAtoms; i++) {
            mask[i] = -1;
        }
        for (int i = 0; i < nAtoms; i++) {
            list.clear();
            for (int j : ip11[i]) {
                list.add(j);
                mask[j] = i;
            }
            keep.clear();
            for (int j : list) {
                Atom aj = atoms[j];
                ArrayList<Bond> bonds = aj.getFFXBonds();
                for (Bond b : bonds) {
                    Atom ak = b.get1_2(aj);
                    int k = ak.getXYZIndex() - 1;
                    if (mask[k] != i) {
                        keep.add(k);
                    }
                }
            }
            list.clear();
            for (int j : keep) {
                for (int k : ip11[j]) {
                    list.add(k);
                }
            }
            Collections.sort(list);
            ip12[i] = new int[list.size()];
            int j = 0;
            for (int k : list) {
                ip12[i][j++] = k;
            }
        }
        /**
         * Find 1-3 group relationships.
         */
        for (int i = 0; i < nAtoms; i++) {
            mask[i] = -1;
        }
        for (int i = 0; i < nAtoms; i++) {
            for (int j : ip11[i]) {
                mask[j] = i;
            }
            for (int j : ip12[i]) {
                mask[j] = i;
            }
            list.clear();
            for (int j : ip12[i]) {
                for (int k : ip12[j]) {
                    if (mask[k] != i) {
                        if (!list.contains(k)) {
                            list.add(k);
                        }
                    }
                }
            }
            ip13[i] = new int[list.size()];
            Collections.sort(list);
            int j = 0;
            for (int k : list) {
                ip13[i][j++] = k;
            }
        }
    }

    /**
     * A recursive method that checks all atoms bonded to the seed atom for
     * inclusion in the polarization group. The method is called on each newly
     * found group member.
     *
     * @param polarizationGroup Atom types that should be included in the group.
     * @param group XYZ indeces of current group members.
     * @param seed The bonds of the seed atom are queried for inclusion in the
     * group.
     */
    private void growGroup(List<Integer> polarizationGroup,
            List<Integer> group, Atom seed) {
        List<Bond> bonds = seed.getFFXBonds();
        for (Bond bi : bonds) {
            Atom aj = bi.get1_2(seed);
            int tj = aj.getType();
            boolean added = false;
            for (int g : polarizationGroup) {
                if (g == tj) {
                    Integer index = aj.getXYZIndex() - 1;
                    if (!group.contains(index)) {
                        group.add(index);
                        added = true;
                        break;
                    }
                }
            }
            if (added) {
                PolarizeType polarizeType = aj.getPolarizeType();
                for (int i : polarizeType.polarizationGroup) {
                    if (!polarizationGroup.contains(i)) {
                        polarizationGroup.add(i);
                    }
                }
                growGroup(polarizationGroup, group, aj);
            }
        }
    }
    /**
     * Number of unique tensors for given order.
     */
    public static final int tensorCount = TensorRecursion.tensorCount(3);
    private final double sfPhi[] = new double[tensorCount];
    private final double sPhi[] = new double[tensorCount];

    public void setEnergyTermState(STATE state) {
    }
    
    @Override
    public void reInit() {
        throw new UnsupportedOperationException(String.format(" No reInit method defined for %s", PME_2.class.toString()));
    }
}
