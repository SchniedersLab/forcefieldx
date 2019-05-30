//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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
package ffx.potential.nonbonded;

import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.SolventRadii;
import ffx.potential.parameters.VDWType;
import ffx.potential.utils.EnergyException;
import static ffx.potential.parameters.ForceField.toEnumForm;
import static ffx.potential.parameters.MultipoleType.ELECTRIC;
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
 * This Generalized Kirkwood class implements GK for the AMOEBA polarizable
 * atomic multipole force field in parallel using a {@link ffx.potential.nonbonded.NeighborList}.
 *
 * @author Michael J. Schnieders<br> derived from:<br> TINKER code by Michael J.
 * Schnieders and Jay W. Ponder<br>
 * @see <a href="http://dx.doi.org/10.1021/ct7001336" target="_blank">M. J.
 * Schnieders and J. W. Ponder, Polarizable atomic multipole solutes in a
 * generalized Kirkwood continuum, Journal of Chemical Theory and Computation
 * 2007, 3, (6), 2083-2097.</a><br>
 */
public class GeneralizedKirkwood implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(GeneralizedKirkwood.class.getName());

    /**
     * Empirical constant that controls the GK cross-term.
     */
    private static final double gkc = 2.455;
    /**
     * Kirkwood monopole reaction field constant.
     */
    private final double fc;
    /**
     * Kirkwood dipole reaction field constant.
     */
    private final double fd;
    /**
     * Kirkwood quadrupole reaction field constant.
     */
    private final double fq;
    /**
     * Permittivity of water at STP.
     */
    private static final double dWater = 78.3;
    /**
     * Default Bondi scale factor.
     */
    private static final double DEFAULT_BONDI_SCALE = 1.03;
    /**
     * Default overlap scale factor for the Hawkins, Cramer & Truhlar pairwise descreening algorithm.
     */
    private static final double DEFAULT_OVERLAP_SCALE = 0.69;
    /**
     * Default surface tension for apolar models with an explicit dispersion term.
     */
    private static final double DEFAULT_CAVDISP_SURFACE_TENSION = 0.08;
    /**
     * Default solvent pressure for apolar models with an explicit volume term.
     */
    private static final double DEFAULT_SOLVENT_PRESSURE = 0.0327;
    /**
     * Default probe radius for use with Gaussian Volumes.
     */
    private static final double DEFAULT_GAUSSVOL_PROBE = 1.92;
    /**
     * Default surface tension for apolar models without an explicit dispersion
     * term. This is lower than CAVDISP, since the favorable dispersion term is
     * implicitly included.
     */
    private static final double DEFAULT_CAV_SURFACE_TENSION = 0.0049;
    /**
     * Empirical scaling of the Bondi radii.
     */
    private double bondiScale;
    /**
     * Cavitation surface tension coefficient (kcal/mol/A^2).
     */
    private final double surfaceTension;
    /**
     * The requested permittivity.
     */
    private double epsilon;
    /**
     * Water probe radius.
     */
    public double probe;
    /**
     * Dielectric offset from:
     * <p>
     * W. C. Still, A. Tempczyk, R. C. Hawley and T. Hendrickson, "A Semianalytical Treatment of Solvation for Molecular
     * Mechanics and Dynamics", J. Amer. Chem. Soc., 112, 6127-6129 (1990)
     */
    private final double dOffset = 0.09;
    /**
     * Force field in use.
     */
    private final ForceField forceField;
    /**
     * Treatment of polarization.
     */
    private final Polarization polarization;
    /**
     * Treatment of non-polar interactions.
     */
    private final NonPolar nonPolar;
    /**
     * Array of Atoms being considered.
     */
    private Atom[] atoms;
    /**
     * Number of atoms.
     */
    private int nAtoms;
    /**
     * Cartesian coordinates of each atom.
     */
    private double[] x, y, z;
    /**
     * Base radius of each atom.
     */
    private double[] baseRadius;
    /**
     * Overlap scale factor for each atom, when using the Hawkins, Cramer & Truhlar pairwise descreening algorithm.
     * <p>
     * G. D. Hawkins, C. J. Cramer and D. G. Truhlar, "Parametrized Models of Aqueous Free Energies of Solvation Based on Pairwise
     * Descreening of Solute Atomic Charges from a Dielectric Medium", J. Phys. Chem., 100, 19824-19839 (1996).
     */
    private double[] overlapScale;
    /**
     * Over-ride the overlap scale factor for hydrogen atoms.
     */
    private final double heavyAtomOverlapScale;
    /**
     * Over-ride the overlap scale factor for hydrogen atoms.
     */
    private final double hydrogenOverlapScale;
    /**
     * Radius of each atom for calculation of dispersion energy.
     */
    private double[] rDisp;
    /**
     * Born radius of each atom.
     */
    private double[] born;
    /**
     * Flag to indicate if an atom should be included.
     */
    private boolean[] use = null;
    /**
     * Periodic boundary conditions and symmetry.
     */
    private Crystal crystal;
    /**
     * Particle mesh Ewald instance, which contains variables such as expanded coordinates and multipoles
     * in the global frame that GK uses.
     */
    private final ParticleMeshEwald particleMeshEwald;
    /**
     * Atomic coordinates for each symmetry operator.
     */
    private double[][][] sXYZ;
    /**
     * Multipole moments for each symmetry operator.
     */
    private double[][][] globalMultipole;
    /**
     * Induced dipoles for each symmetry operator.
     */
    private double[][][] inducedDipole;
    /**
     * Induced dipole chain rule terms for each symmetry operator.
     */
    private double[][][] inducedDipoleCR;
    /**
     * Gradient array for each thread.
     */
    private double[][][] grad;
    /**
     * Torque array for each thread.
     */
    private double[][][] torque;
    /**
     * Lambda gradient array for each thread (dU/dX/dL)
     */
    private double[][][] lambdaGrad;
    /**
     * Lambda torque array for each thread.
     */
    private double[][][] lambdaTorque;
    /**
     * Neighbor lists for each atom and symmetry operator.
     */
    private int[][][] neighborLists;
    /**
     * This field is because re-initializing the force field resizes some arrays
     * but not others; that second category must, when called on, be resized not
     * to the current number of atoms but to the maximum number of atoms (and
     * thus to the size of the first category of arrays).
     */
    private int maxNumAtoms;
    /**
     * Parallel team object for shared memory parallelization.
     */
    private final ParallelTeam parallelTeam;
    /**
     * Parallel computation of Born Radii.
     */
    private final BornRadiiRegion bornRadiiRegion;
    /**
     * Parallel computation of the Permanent GK Field.
     */
    private final PermanentGKFieldRegion permanentGKFieldRegion;
    /**
     * Parallel computation of the Induced GK Field.
     */
    private final InducedGKFieldRegion inducedGKFieldRegion;
    /**
     * Parallel computation of the GK continuum electrostatics energy.
     */
    private final GKEnergyRegion gkEnergyRegion;
    /**
     * Parallel computation of Born radii chain rule term.
     */
    private final BornCRRegion bornGradRegion;
    /**
     * Parallel computation of Dispersion energy.
     */
    private final DispersionRegion dispersionRegion;
    /**
     * Parallel computation of Cavitation.
     */
    private final CavitationRegion cavitationRegion;
    /**
     * Parallel computation of Volume.
     */
    private final VolumeRegion volumeRegion;
    /**
     * Gaussian Based Volume and Surface Area
     */
    private final GaussVol gaussVol;
    /**
     * Shared array for computation of Born radii gradient.
     */
    private SharedDoubleArray sharedBornGrad;
    /**
     * Shared arrays for computation of the GK field for each symmetry operator.
     */
    protected SharedDoubleArray[] sharedGKField;
    /**
     * Shared arrays for computation of the GK field chain-rule term for each symmetry operator.
     */
    protected SharedDoubleArray[] sharedGKFieldCR;
    /**
     * GK cut-off distance.
     */
    private double cutoff;
    /**
     * GK cut-off distance squared.
     */
    private double cut2;
    /**
     * Boolean flag to indicate GK will be scaled by the lambda state variable.
     */
    private boolean lambdaTerm;
    /**
     * The current value of the lambda state variable.
     */
    private double lambda = 1.0;
    /**
     * lPow equals lambda^polarizationLambdaExponent, where polarizationLambdaExponent also used by PME.
     */
    private double lPow = 1.0;
    /**
     * First derivative of lPow with respect to l.
     */
    private double dlPow = 0.0;
    /**
     * Second derivative of lPow with respect to l.
     */
    private double dl2Pow = 0.0;
    /**
     * Electrostatic Solvation Energy.
     */
    private double solvationEnergy = 0.0;
    /**
     * Dispersion Solvation Energy.
     */
    private double dispersionEnergy = 0.0;
    /**
     * Cavitation Solvation Energy.
     */
    private double cavitationEnergy = 0.0;
    /**
     * Time to compute GK electrostatics.
     */
    private long gkTime = 0;
    /**
     * Time to compute Dispersion energy.
     */
    private long dispersionTime = 0;
    /**
     * Time to compute Cavitation energy.
     */
    private long cavitationTime = 0;
    /**
     * Use base radii defined by AtomType rather than by atomic number.
     */
    private boolean verboseRadii;
    /**
     * If true, prevents Born radii from updating.
     */
    private boolean fixedRadii = false;
    /**
     * Forces all atoms to be considered during Born radius updates.
     */
    private boolean nativeEnvironmentApproximation;
    /**
     * Maps radii overrides (by AtomType) specified from the command line. e.g.
     * -DradiiOverride=134r1.20,135r1.20 sets atom types 134,135 to Bondi=1.20
     */
    private final HashMap<Integer, Double> radiiOverride = new HashMap<>();
    /**
     * Control GK logging.
     */
    private final Level GK_WARN_LEVEL;
    /**
     * Flag to indicate calculation of molecular volume.
     */
    private boolean doVolume;

    /**
     * <p>
     * Constructor for GeneralizedKirkwood.</p>
     *
     * @param forceField        a {@link ffx.potential.parameters.ForceField} object.
     * @param atoms             an array of {@link ffx.potential.bonded.Atom} objects.
     * @param particleMeshEwald a {@link ffx.potential.nonbonded.ParticleMeshEwald} object.
     * @param crystal           a {@link ffx.crystal.Crystal} object.
     * @param parallelTeam      a {@link edu.rit.pj.ParallelTeam} object.
     */
    public GeneralizedKirkwood(ForceField forceField, Atom[] atoms,
                               ParticleMeshEwald particleMeshEwald, Crystal crystal,
                               ParallelTeam parallelTeam) {

        this.forceField = forceField;
        this.atoms = atoms;
        this.particleMeshEwald = particleMeshEwald;
        this.crystal = crystal;
        this.parallelTeam = parallelTeam;
        nAtoms = atoms.length;
        maxNumAtoms = nAtoms;
        polarization = particleMeshEwald.polarization;

        String suppressGKwarnings = System.getProperty("gk-suppressWarnings");
        if (suppressGKwarnings != null && Boolean.parseBoolean(suppressGKwarnings)) {
            GK_WARN_LEVEL = Level.FINE;
        } else {
            GK_WARN_LEVEL = Level.WARNING;
        }

        // Set the Kirkwood multipolar reaction field constants.
        epsilon = forceField.getDouble(ForceField.ForceFieldDouble.GK_EPSILON, dWater);
        fc = 1.0 * (1.0 - epsilon) / (0.0 + 1.0 * epsilon);
        fd = 2.0 * (1.0 - epsilon) / (1.0 + 2.0 * epsilon);
        fq = 3.0 * (1.0 - epsilon) / (2.0 + 3.0 * epsilon);

        // Define default Bondi scale factor, and HCT overlap scale factors.
        bondiScale = forceField.getDouble(ForceField.ForceFieldDouble.GK_BONDIOVERRIDE, DEFAULT_BONDI_SCALE);
        heavyAtomOverlapScale = forceField.getDouble(ForceField.ForceFieldDouble.GK_OVERLAPSCALE, DEFAULT_OVERLAP_SCALE);
        hydrogenOverlapScale = forceField.getDouble(ForceField.ForceFieldDouble.GK_HYDROGEN_OVERLAPSCALE, heavyAtomOverlapScale);

        // Process any radii override values.
        String radiiProp = forceField.getString(ForceField.ForceFieldString.GK_RADIIOVERRIDE, null);
        if (radiiProp != null) {
            String[] tokens = radiiProp.split("A");
            for (String token : tokens) {
                if (!token.contains("r")) {
                    logger.severe("Invalid radius override.");
                }
                int separator = token.indexOf("r");
                int type = Integer.parseInt(token.substring(0, separator));
                double factor = Double.parseDouble(token.substring(separator + 1));
                logger.info(format(" (GK) Scaling AtomType %d with bondi factor %.2f", type, factor));
                radiiOverride.put(type, factor);
            }
        }

        NonPolar nonpolarModel;
        try {
            String cavModel = forceField.getString(ForceField.ForceFieldString.CAVMODEL, "CAV_DISP").toUpperCase();
            nonpolarModel = getNonPolarModel(cavModel);
        } catch (Exception ex) {
            nonpolarModel = NonPolar.NONE;
            logger.warning(format(" Error parsing non-polar model (set to NONE) %s", ex.toString()));
        }
        nonPolar = nonpolarModel;

        sharedGKField = new SharedDoubleArray[3];
        sharedGKFieldCR = new SharedDoubleArray[3];

        nativeEnvironmentApproximation = forceField.getBoolean(
                ForceField.ForceFieldBoolean.NATIVE_ENVIRONMENT_APPROXIMATION, false);
        probe = forceField.getDouble(ForceField.ForceFieldDouble.PROBE_RADIUS, 1.4);
        cutoff = forceField.getDouble(ForceField.ForceFieldDouble.GK_CUTOFF, particleMeshEwald.getEwaldCutoff());
        cut2 = cutoff * cutoff;
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.GK_LAMBDATERM,
                forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false));

        /*
          If polarization lambda exponent is set to 0.0, then we're running
          Dual-Topology and the GK energy will be scaled with the overall system lambda value.
         */
        double polLambdaExp = forceField.getDouble(ForceField.ForceFieldDouble.POLARIZATION_LAMBDA_EXPONENT, 3.0);
        if (polLambdaExp == 0.0) {
            lambdaTerm = false;
            logger.info(" GK lambda term set to false.");
        }

        // If PME includes polarization and is a function of lambda, GK must also.
        if (!lambdaTerm && particleMeshEwald.getPolarizationType() != Polarization.NONE) {
            if (forceField.getBoolean(ForceField.ForceFieldBoolean.ELEC_LAMBDATERM,
                    forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false))) {
                logger.info(" If PME includes polarization and is a function of lambda, GK must also.");
                lambdaTerm = true;
            }
        }

        int threadCount = parallelTeam.getThreadCount();
        bornRadiiRegion = new BornRadiiRegion(threadCount);
        permanentGKFieldRegion = new PermanentGKFieldRegion(threadCount);
        inducedGKFieldRegion = new InducedGKFieldRegion(threadCount);
        gkEnergyRegion = new GKEnergyRegion(threadCount);
        bornGradRegion = new BornCRRegion(threadCount);

        initAtomArrays();

        doVolume = forceField.getBoolean(ForceField.ForceFieldBoolean.VOLUME, false);

        double tensionDefault;
        switch (nonPolar) {
            case CAV:
                tensionDefault = DEFAULT_CAV_SURFACE_TENSION;
                cavitationRegion = new CavitationRegion(atoms, x, y, z, use, neighborLists, grad, lambdaGrad,
                        threadCount, probe, lambdaTerm, tensionDefault);
                volumeRegion = new VolumeRegion(atoms, x, y, z, tensionDefault, threadCount);
                dispersionRegion = null;
                gaussVol = null;
                break;
            case CAV_DISP:
                tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
                cavitationRegion = new CavitationRegion(atoms, x, y, z, use, neighborLists, grad, lambdaGrad,
                        threadCount, probe, lambdaTerm, tensionDefault);
                volumeRegion = new VolumeRegion(atoms, x, y, z, tensionDefault, threadCount);
                dispersionRegion = new DispersionRegion(threadCount);
                gaussVol = null;
                break;
            case GAUSS_DISP:
                dispersionRegion = new DispersionRegion(threadCount);
                cavitationRegion = null;
                volumeRegion = null;

                gaussVol = new GaussVol(nAtoms, null);
                tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
                boolean[] isHydrogen = new boolean[nAtoms];
                double[] radii = new double[nAtoms];
                double[] volume = new double[nAtoms];
                double[] gamma = new double[nAtoms];

                double fourThirdsPI = 4.0 / 3.0 * PI;
                double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0);

                probe = forceField.getDouble(ForceField.ForceFieldDouble.PROBE_RADIUS, DEFAULT_GAUSSVOL_PROBE);
                int index = 0;
                for (Atom atom : atoms) {
                    isHydrogen[index] = atom.isHydrogen();
                    radii[index] = atom.getVDWType().radius / 2.0 * rminToSigma;
                    radii[index] += probe;
                    volume[index] = fourThirdsPI * pow(radii[index], 3);
                    gamma[index] = 1.0;
                    index++;
                }

                try {
                    gaussVol.setGammas(gamma);
                    gaussVol.setRadii(radii);
                    gaussVol.setVolumes(volume);
                    gaussVol.setIsHydrogen(isHydrogen);
                } catch (Exception e) {
                    logger.severe(" Exception creating GaussVol: " + e.toString());
                }
                break;
            case BORN_CAV_DISP:
                tensionDefault = DEFAULT_CAVDISP_SURFACE_TENSION;
                cavitationRegion = null;
                volumeRegion = null;
                gaussVol = null;
                dispersionRegion = new DispersionRegion(threadCount);
                break;
            case HYDROPHOBIC_PMF:
            case BORN_SOLV:
            case NONE:
            default:
                tensionDefault = DEFAULT_CAV_SURFACE_TENSION;
                cavitationRegion = null;
                volumeRegion = null;
                dispersionRegion = null;
                gaussVol = null;
                break;
        }

        surfaceTension = forceField.getDouble(ForceField.ForceFieldDouble.SURFACE_TENSION, tensionDefault);

        logger.info("  Continuum Solvation ");
        logger.info(format("   Generalized Kirkwood Cut-Off:       %8.3f (A)", cutoff));
        logger.info(format("   Solvent Dielectric:                 %8.3f", epsilon));
        SolventRadii.logRadiiSource(forceField);
        logger.info(format("   Non-Polar Model:                  %10s",
                nonPolar.toString().replace('_', '-')));

        if (cavitationRegion != null) {
            logger.info(format("   Cavitation Probe Radius:            %8.3f (A)", probe));
            logger.info(format("   Cavitation Surface Tension:         %8.3f (Kcal/mol/A^2)", surfaceTension));
        } else if (gaussVol != null) {
            logger.info(format("   Cavitation Probe Radius:            %8.3f (A)", probe));
            logger.info(format("   Cavitation Solvent Pressure:        %8.3f (Kcal/mol/A^3)", gaussVol.getSolventPressure()));
            logger.info(format("   Cavitation Surface Tension:         %8.3f (Kcal/mol/A^2)", gaussVol.getSurfaceTension()));
        }


        // Print out all Base Radii
        if (logger.isLoggable(Level.FINE)) {
            logger.fine("   GK Base Radii");
            for (int i = 0; i < nAtoms; i++) {
                logger.info(format("   %s %8.6f %5.3f", atoms[i].toString(), baseRadius[i], overlapScale[i]));
            }
        }

    }

    /**
     * <p>Getter for the field <code>overlapScale</code>.</p>
     *
     * @return an array of {@link double} objects.
     */
    public double[] getOverlapScale() {
        return overlapScale;
    }

    /**
     * <p>getBaseRadii.</p>
     *
     * @return an array of {@link double} objects.
     */
    public double[] getBaseRadii() {
        return baseRadius;
    }

    /**
     * <p>Getter for the field <code>surfaceTension</code>.</p>
     *
     * @return a double.
     */
    public double getSurfaceTension() {
        return surfaceTension;
    }

    /**
     * <p>getNonPolarModel.</p>
     *
     * @return a {@link ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar} object.
     */
    public NonPolar getNonPolarModel() {
        return nonPolar;
    }

    /**
     * Returns the dielectric offset (in Angstroms).
     *
     * @return Currently: 0.09 Angstroms.
     */
    public double getDielecOffset() {
        return dOffset;
    }

    /**
     * Returns the solvent relative permittivity (typically 78.3).
     *
     * @return Relative permittivity of solvent.
     */
    public double getSolventPermittivity() {
        return epsilon;
    }

    /**
     * Returns the probe radius (typically 1.4 Angstroms).
     *
     * @return Radius of the solvent probe.
     */
    public double getProbeRadius() {
        return probe;
    }

    /**
     * <p>Setter for the field <code>cutoff</code>.</p>
     *
     * @param cutoff a double.
     */
    public void setCutoff(double cutoff) {
        this.cutoff = cutoff;
        this.cut2 = cutoff * cutoff;
    }

    /**
     * <p>Getter for the field <code>cutoff</code>.</p>
     *
     * @return a double.
     */
    public double getCutoff() {
        return cutoff;
    }

    /**
     * <p>Setter for the field <code>crystal</code>.</p>
     *
     * @param crystal a {@link ffx.crystal.Crystal} object.
     */
    public void setCrystal(Crystal crystal) {
        this.crystal = crystal;
    }

    /**
     * <p>setNeighborList.</p>
     *
     * @param neighbors an array of {@link int} objects.
     */
    public void setNeighborList(int[][][] neighbors) {
        this.neighborLists = neighbors;
    }

    /**
     * <p>Setter for the field <code>atoms</code>.</p>
     *
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public void setAtoms(Atom[] atoms) {
        this.atoms = atoms;
        nAtoms = atoms.length;
        maxNumAtoms = nAtoms > maxNumAtoms ? nAtoms : maxNumAtoms;
        initAtomArrays();
    }

    /**
     * <p>Setter for the field <code>fixedRadii</code>.</p>
     *
     * @param fixedRadii a boolean.
     */
    public void setFixedRadii(boolean fixedRadii) {
        this.fixedRadii = fixedRadii;
    }

    /**
     * <p>Getter for the field <code>fixedRadii</code>.</p>
     *
     * @return a boolean.
     */
    public boolean getFixedRadii() {
        return fixedRadii;
    }

    /**
     * <p>
     * Sets whether GK uses the Native Environment Approximation.
     * </p>
     * <p>
     * This (previously known as born-use-all) is useful for rotamer
     * optimization under continuum solvent. If a large number of sidechain
     * atoms are completely removed from the GK/GB calculation, the remaining
     * sidechains are overly solvated. The NEA says "we will keep all sidechains
     * not under optimization in some default state and let them contribute to
     * Born radii calculations, but still exclude their solvation energy
     * components."
     * </p>
     *
     * @param nativeEnvironmentApproximation Whether to use the NEA.
     */
    public void setNativeEnvironmentApproximation(boolean nativeEnvironmentApproximation) {
        this.nativeEnvironmentApproximation = nativeEnvironmentApproximation;
    }

    /**
     * <p>
     * Checks whether GK uses the Native Environment Approximation.
     * </p>
     * <p>
     * This (previously known as born-use-all) is useful for rotamer
     * optimization under continuum solvent. If a large number of sidechain
     * atoms are completely removed from the GK/GB calculation, the remaining
     * sidechains are overly solvated. The NEA says "we will keep all sidechains
     * not under optimization in some default state and let them contribute to
     * Born radii calculations, but still exclude their solvation energy
     * components."
     * </p>
     *
     * @return Whether the NEA is in use.
     */
    public boolean getNativeEnvironmentApproximation() {
        return nativeEnvironmentApproximation;
    }

    private void initAtomArrays() {
        if (fixedRadii) {
            fixedRadii = false;
        }

        sXYZ = particleMeshEwald.coordinates;

        x = sXYZ[0][0];
        y = sXYZ[0][1];
        z = sXYZ[0][2];

        globalMultipole = particleMeshEwald.globalMultipole;
        inducedDipole = particleMeshEwald.inducedDipole;
        inducedDipoleCR = particleMeshEwald.inducedDipoleCR;
        neighborLists = particleMeshEwald.neighborLists;
        grad = particleMeshEwald.getGradient();
        torque = particleMeshEwald.getTorque();
        lambdaGrad = particleMeshEwald.getLambdaGradient();
        lambdaTorque = particleMeshEwald.getLambdaTorque();

        if (sharedGKField[0] == null || sharedGKField[0].length() < nAtoms) {
            sharedGKField[0] = new SharedDoubleArray(nAtoms);
            sharedGKField[1] = new SharedDoubleArray(nAtoms);
            sharedGKField[2] = new SharedDoubleArray(nAtoms);
            sharedGKFieldCR[0] = new SharedDoubleArray(nAtoms);
            sharedGKFieldCR[1] = new SharedDoubleArray(nAtoms);
            sharedGKFieldCR[2] = new SharedDoubleArray(nAtoms);
            sharedBornGrad = new SharedDoubleArray(nAtoms);
            baseRadius = new double[nAtoms];
            overlapScale = new double[nAtoms];
            rDisp = new double[nAtoms];
            born = new double[nAtoms];
            use = new boolean[nAtoms];
        }

        fill(use, true);

        bondiScale = SolventRadii.applyGKRadii(forceField, bondiScale, atoms, baseRadius);

        // Set up HCT overlap scale factors and any requested radii overrides.
        for (int i = 0; i < nAtoms; i++) {
            overlapScale[i] = heavyAtomOverlapScale;
            int atomicNumber = atoms[i].getAtomicNumber();
            if (atomicNumber == 1) {
                overlapScale[i] = hydrogenOverlapScale;
            }

            AtomType atomType = atoms[i].getAtomType();
            if (radiiOverride.containsKey(atomType.type)) {
                double override = radiiOverride.get(atomType.type);
                // Remove default bondiFactor, and apply override.
                baseRadius[i] = baseRadius[i] * override / bondiScale;
                logger.info(format(" Scaling %s (atom type %d) to %7.4f (Bondi factor %7.4f)",
                        atoms[i], atomType.type, baseRadius[i], override));
            }
        }

        // Resets verboseRadii; reduces logging messages when mutating MultiResidues.
        if (dispersionRegion != null) {
            dispersionRegion.init();
        }

        if (cavitationRegion != null) {
            cavitationRegion.init();
        }

        if (gaussVol != null) {
            boolean[] isHydrogen = new boolean[nAtoms];
            double[] radii = new double[nAtoms];
            double[] volume = new double[nAtoms];
            double[] gamma = new double[nAtoms];

            double fourThirdsPI = 4.0 / 3.0 * PI;
            double rminToSigma = 1.0 / pow(2.0, 1.0 / 6.0);

            probe = forceField.getDouble(ForceField.ForceFieldDouble.PROBE_RADIUS, DEFAULT_GAUSSVOL_PROBE);

            int index = 0;
            for (Atom atom : atoms) {
                isHydrogen[index] = atom.isHydrogen();
                radii[index] = atom.getVDWType().radius / 2.0 * rminToSigma;
                radii[index] += probe;
                volume[index] = fourThirdsPI * pow(radii[index], 3);
                gamma[index] = 1.0;
                index++;
            }

            try {
                gaussVol.setGammas(gamma);
                gaussVol.setRadii(radii);
                gaussVol.setVolumes(volume);
                gaussVol.setIsHydrogen(isHydrogen);
            } catch (Exception e) {
                logger.severe(" Exception creating GaussVol: " + e.toString());
            }

        }

    }

    /**
     * <p>Setter for the field <code>use</code>.</p>
     *
     * @param use an array of {@link boolean} objects.
     */
    public void setUse(boolean[] use) {
        this.use = use;
    }

    /**
     * <p>
     * computeBornRadii</p>
     */
    public void computeBornRadii() {

        // Born radii are fixed.
        if (fixedRadii) {
            return;
        }

        try {
            parallelTeam.execute(bornRadiiRegion);
        } catch (Exception e) {
            String message = "Fatal exception computing Born radii.";
            logger.log(Level.SEVERE, message, e);
        }

        for (int i = 0; i < nAtoms; i++) {
            if (use[i]) {
                double borni = born[i];
                if (isInfinite(borni) || isNaN(borni)) {
                    throw new EnergyException(format(" %s\n Born radii %d %8.3f", atoms[i], i, born[i]), true);
                }
            }
            // logger.info(format(" %s Born radii %d %8.3f", atoms[i], i, born[i]));
        }
    }

    /**
     * <p>
     * computePermanentGKField</p>
     */
    public void computePermanentGKField() {
        try {
            parallelTeam.execute(permanentGKFieldRegion);
        } catch (Exception e) {
            String message = "Fatal exception computing permanent GK field.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * <p>
     * computeInducedGKField</p>
     */
    public void computeInducedGKField() {
        try {
            parallelTeam.execute(inducedGKFieldRegion);
        } catch (Exception e) {
            String message = "Fatal exception computing induced GK field.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * <p>
     * solvationEnergy</p>
     *
     * @param gradient a boolean.
     * @param print    a boolean.
     * @return a double.
     */
    public double solvationEnergy(boolean gradient, boolean print) {

        // Initialize the gradient accumulation arrays.
        if (gradient) {
            for (int j = 0; j < nAtoms; j++) {
                sharedBornGrad.set(j, 0.0);
            }
        }

        try {
            // Find the GK energy.
            gkTime = -System.nanoTime();
            gkEnergyRegion.setGradient(gradient);
            parallelTeam.execute(gkEnergyRegion);
            gkTime += System.nanoTime();

            // Find the nonpolar energy.
            switch (nonPolar) {
                case CAV:
                    cavitationTime = -System.nanoTime();
                    parallelTeam.execute(cavitationRegion);
                    if (doVolume) {
                        ParallelTeam serialTeam = new ParallelTeam(1);
                        serialTeam.execute(volumeRegion);
                    }
                    cavitationTime += System.nanoTime();
                    break;
                case CAV_DISP:
                    dispersionTime = -System.nanoTime();
                    dispersionRegion.setGradient(gradient);
                    parallelTeam.execute(dispersionRegion);
                    dispersionTime += System.nanoTime();
                    cavitationTime = -System.nanoTime();
                    parallelTeam.execute(cavitationRegion);
                    if (doVolume) {
                        ParallelTeam serialTeam = new ParallelTeam(1);
                        serialTeam.execute(volumeRegion);
                    }
                    cavitationTime += System.nanoTime();
                    break;
                case GAUSS_DISP:
                    dispersionTime = -System.nanoTime();
                    dispersionRegion.setGradient(gradient);
                    parallelTeam.execute(dispersionRegion);
                    dispersionTime += System.nanoTime();

                    cavitationTime = -System.nanoTime();
                    double[][] positions = new double[nAtoms][3];
                    for (int i = 0; i < nAtoms; i++) {
                        positions[i][0] = atoms[i].getX();
                        positions[i][1] = atoms[i].getY();
                        positions[i][2] = atoms[i].getZ();
                    }
                    gaussVol.energyAndGradient(positions, grad[0]);
                    cavitationTime += System.nanoTime();

                    break;
                case BORN_CAV_DISP:
                    dispersionTime = -System.nanoTime();
                    dispersionRegion.setGradient(gradient);
                    parallelTeam.execute(dispersionRegion);
                    dispersionTime += System.nanoTime();
                    break;
                case HYDROPHOBIC_PMF:
                case BORN_SOLV:
                case NONE:
                default:
                    break;
            }
        } catch (Exception e) {
            String message = "Fatal exception computing the continuum solvation energy.";
            logger.log(Level.SEVERE, message, e);
        }

        // Compute the Born radii chain rule term.
        if (gradient) {
            try {
                gkTime -= System.nanoTime();
                parallelTeam.execute(bornGradRegion);
                gkTime += System.nanoTime();
            } catch (Exception e) {
                String message = "Fatal exception computing Born radii chain rule term.";
                logger.log(Level.SEVERE, message, e);
            }

            for (int i = 0; i < nAtoms; i++) {
                double borni = sharedBornGrad.get(i);
                if (use[i]) {
                    if (isInfinite(borni) || isNaN(borni)) {
                        throw new EnergyException(format(" %s\n Born radii %d %8.3f", atoms[i], i, borni), true);
                    }
                }
                // logger.info(format(" %s Born radius grad %d %8.3f", atoms[i], i, borni));
            }

        }

        if (print) {
            logger.info(format(" Generalized Kirkwood%16.8f %10.3f",
                    gkEnergyRegion.getEnergy(), gkTime * 1e-9));
            switch (nonPolar) {
                case CAV:
                    if (doVolume) {
                        cavitationEnergy = volumeRegion.getEnergy();
                    } else {
                        cavitationEnergy = cavitationRegion.getEnergy();
                    }
                    logger.info(format(" Cavitation          %16.8f %10.3f",
                            cavitationEnergy, cavitationTime * 1e-9));
                    break;
                case CAV_DISP:
                    if (doVolume) {
                        cavitationEnergy = volumeRegion.getEnergy();
                    } else {
                        cavitationEnergy = cavitationRegion.getEnergy();
                    }
                    dispersionEnergy = dispersionRegion.getEnergy();
                    logger.info(format(" Cavitation          %16.8f %10.3f",
                            cavitationEnergy, cavitationTime * 1e-9));
                    logger.info(format(" Dispersion          %16.8f %10.3f",
                            dispersionEnergy, dispersionTime * 1e-9));
                    break;
                case GAUSS_DISP:
                    cavitationEnergy = gaussVol.getEnergy();
                    dispersionEnergy = dispersionRegion.getEnergy();
                    logger.info(format(" Cavitation          %16.8f %10.3f",
                            cavitationEnergy, cavitationTime * 1e-9));
                    logger.info(format(" Dispersion          %16.8f %10.3f",
                            dispersionEnergy, dispersionTime * 1e-9));
                    break;
                case BORN_CAV_DISP:
                    dispersionEnergy = dispersionRegion.getEnergy();
                    logger.info(format(" Dispersion          %16.8f %10.3f",
                            dispersionEnergy, dispersionTime * 1e-9));
                    break;
                case HYDROPHOBIC_PMF:
                case BORN_SOLV:
                case NONE:
                default:
                    break;
            }
        }

        switch (nonPolar) {
            case CAV:
                if (doVolume) {
                    solvationEnergy = gkEnergyRegion.getEnergy() + volumeRegion.getEnergy();
                } else {
                    solvationEnergy = gkEnergyRegion.getEnergy() + cavitationRegion.getEnergy();
                }
                break;
            case CAV_DISP:
                if (doVolume) {
                    solvationEnergy = gkEnergyRegion.getEnergy() + dispersionRegion.getEnergy()
                            + volumeRegion.getEnergy();
                } else {
                    solvationEnergy = gkEnergyRegion.getEnergy() + dispersionRegion.getEnergy()
                            + cavitationRegion.getEnergy();
                }
                break;
            case GAUSS_DISP:
                solvationEnergy = gkEnergyRegion.getEnergy() + dispersionRegion.getEnergy()
                        + gaussVol.getEnergy();
                break;
            case BORN_CAV_DISP:
                solvationEnergy = gkEnergyRegion.getEnergy() + dispersionRegion.getEnergy();
                break;
            case HYDROPHOBIC_PMF:
            case BORN_SOLV:
            case NONE:
            default:
                solvationEnergy = gkEnergyRegion.getEnergy();
                break;
        }

        if (lambdaTerm) {
            return lPow * solvationEnergy;
        } else {
            return solvationEnergy;
        }
    }

    /**
     * Returns the cavitation component (if applicable) of GK energy. If this GK
     * is operating without a cavitation term, it either returns 0, or throws an
     * error if throwError is true.
     *
     * @param throwError a boolean.
     * @return Cavitation energy
     */
    public double getCavitationEnergy(boolean throwError) {
        switch (nonPolar) {
            case CAV:
            case CAV_DISP:
                return cavitationEnergy;
            default:
                if (throwError) {
                    throw new IllegalArgumentException(" GK is operating without a cavitation term");
                } else {
                    return 0.0;
                }
        }
    }

    /**
     * Returns the dispersion component (if applicable) of GK energy. If this GK
     * is operating without a dispersion term, it either returns 0, or throws an
     * error if throwError is true.
     *
     * @param throwError a boolean.
     * @return Cavitation energy
     */
    public double getDispersionEnergy(boolean throwError) {
        switch (nonPolar) {
            case CAV_DISP:
            case BORN_CAV_DISP:
                return dispersionEnergy;
            default:
                if (throwError) {
                    throw new IllegalArgumentException(" GK is operating without a dispersion term");
                } else {
                    return 0.0;
                }
        }
    }

    /**
     * <p>
     * getInteractions</p>
     *
     * @return a int.
     */
    public int getInteractions() {
        return gkEnergyRegion.getInteractions();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Updates the value of lPow.
     */
    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            this.lambda = lambda;
        } else {

            // If the lambdaTerm flag is false, lambda must be set to one.
            this.lambda = 1.0;
            lPow = 1.0;
            dlPow = 0.0;
            dl2Pow = 0.0;
        }
    }

    void setLambdaFunction(double lPow, double dlPow, double dl2Pow) {
        if (lambdaTerm) {
            this.lPow = lPow;
            this.dlPow = dlPow;
            this.dl2Pow = dl2Pow;
        } else {
            // If the lambdaTerm flag is false, lambda must be set to one.
            this.lambda = 1.0;
            this.lPow = 1.0;
            this.dlPow = 0.0;
            this.dl2Pow = 0.0;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getLambda() {
        return lambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        if (lambdaTerm) {
            return dlPow * solvationEnergy;
        }
        return 0.0;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The 2nd derivative is 0.0. (U=Lambda*Egk, dU/dL=Egk, d2U/dL2=0.0)
     */
    @Override
    public double getd2EdL2() {
        if (lambdaTerm) {
            return dl2Pow * solvationEnergy;
        } else {
            return 0.0;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * These contributions are already aggregated into the arrays used by PME.
     */
    @Override
    public void getdEdXdL(double[] gradient) {
    }

    /**
     * <p>getNonPolarModel.</p>
     *
     * @param nonpolarModel a {@link java.lang.String} object.
     * @return a {@link ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar} object.
     */
    public static NonPolar getNonPolarModel(String nonpolarModel) {
        try {
            return NonPolar.valueOf(toEnumForm(nonpolarModel));
        } catch (IllegalArgumentException ex) {
            logger.warning(" Unrecognized nonpolar model requested; defaulting to NONE.");
            return NonPolar.NONE;
        }
    }

    /**
     * Compute Born radii in parallel via the Grycuk method.
     *
     * @since 1.0
     */
    private class BornRadiiRegion extends ParallelRegion {

        private final BornRadiiLoop[] bornRadiiLoop;
        private SharedDoubleArray sharedBorn;
        private SharedDouble ecavTot;

        BornRadiiRegion(int nt) {
            bornRadiiLoop = new BornRadiiLoop[nt];
            for (int i = 0; i < nt; i++) {
                bornRadiiLoop[i] = new BornRadiiLoop();
            }
            ecavTot = new SharedDouble(0.0);
        }

        @Override
        public void start() {
            if (sharedBorn == null || sharedBorn.length() < nAtoms) {
                sharedBorn = new SharedDoubleArray(nAtoms);
            }
            for (int i = 0; i < nAtoms; i++) {
                sharedBorn.set(i, 0.0);
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, bornRadiiLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing Born radii in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        @Override
        public void finish() {
            for (int i = 0; i < nAtoms; i++) {
                final double baseRi = baseRadius[i];
                if (!use[i]) {
                    born[i] = baseRi;
                } else {
                    double sum = sharedBorn.get(i);
                    if (sum <= 0.0) {
                        sum = PI4_3 * 1.0e-9;
                        born[i] = 1.0 / pow(sum / PI4_3, oneThird);
                        //logger.info(format(" I < 0; Resetting %d to %12.6f", i, born[i]));
                        logger.log(GK_WARN_LEVEL, format(" I < 0; Resetting %d to %12.6f", i, born[i]));
                        continue;
                    }
                    born[i] = 1.0 / pow(sum / PI4_3, oneThird);
                    if (verboseRadii) {
                        logger.info(format(" Atom %s: Base radius %10.6f, Born radius %10.6f", atoms[i], baseRi, born[i]));
                    }
                    if (born[i] < baseRi) {
                        // logger.info(format(" Less than base radii; resetting to %d %12.6f", i, baseRi));
                        born[i] = baseRi;
                        continue;
                    }
                    if (isInfinite(born[i]) || isNaN(born[i])) {
                        logger.info(format(" NaN / Infinite: Resetting Base Radii %d %12.6f", i, baseRi));
                        born[i] = baseRi;
                    }
                }
            }
            if (verboseRadii) {
                // This could get very verbose if printed at each step.
                logger.info(" Disabling verbose radii printing.");
                verboseRadii = false;
            }
        }

        /**
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class BornRadiiLoop extends IntegerForLoop {

            private double[] localBorn;
            private double ecav;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            BornRadiiLoop() {
                ecav = 0.0;
            }

            @Override
            public void start() {
                if (localBorn == null || localBorn.length < nAtoms) {
                    localBorn = new double[nAtoms];
                }
                fill(localBorn, 0.0);
            }

            @Override
            public void finish() {
                sharedBorn.reduce(localBorn, DoubleOp.SUM);
                ecavTot.addAndGet(ecav);
            }

            /**
             * Use pairwise descreening to compute integral of 1/r^6.
             *
             * @param r            atomic separation.
             * @param r2           atomic separation squared.
             * @param radius       base radius of the atom being descreened.
             * @param scaledRadius scaled raduis of the atom doing the descreening.
             * @return this contribution to the descreening integral.
             */
            private double integral(double r, double r2, double radius, double scaledRadius) {
                double integral = 0.0;

                // Descreen only if atom I does not engulf atom K.
                if (radius < r + scaledRadius) {
                    // Atom i is engulfed by atom k.
                    if (radius + r < scaledRadius) {
                        final double lower = radius;
                        final double upper = scaledRadius - r;
                        integral = (PI4_3 * (1.0 / (upper * upper * upper) - 1.0 / (lower * lower * lower)));
                    }

                    // Upper integration bound is always the same.
                    double upper = r + scaledRadius;

                    // Lower integration bound depends on atoms sizes and separation.
                    double lower;
                    if (radius + r < scaledRadius) {
                        // Atom i is engulfed by atom k.
                        lower = scaledRadius - r;
                    } else if (r < radius + scaledRadius) {
                        // Atoms are overlapped, begin integration from ri.
                        lower = radius;
                    } else {
                        // No overlap between atoms.
                        lower = r - scaledRadius;
                    }

                    double l2 = lower * lower;
                    double l4 = l2 * l2;
                    double lr = lower * r;
                    double l4r = l4 * r;
                    double u2 = upper * upper;
                    double u4 = u2 * u2;
                    double ur = upper * r;
                    double u4r = u4 * r;
                    double scaledRk2 = scaledRadius * scaledRadius;
                    double term = (3.0 * (r2 - scaledRk2) + 6.0 * u2 - 8.0 * ur) / u4r
                            - (3.0 * (r2 - scaledRk2) + 6.0 * l2 - 8.0 * lr) / l4r;
                    integral -= PI_12 * term;
                }

                return integral;
            }

            @Override
            public void run(int lb, int ub) {
                // The descreening integral is initialized to the limit of the atom alone in solvent.
                for (int i = lb; i <= ub; i++) {
                    final double baseRi = baseRadius[i];
                    localBorn[i] = PI4_3 / (baseRi * baseRi * baseRi);
                }

                int nSymm = crystal.spaceGroup.symOps.size();
                if (nSymm == 0) {
                    nSymm = 1;
                }

                for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                    double[][] xyz = sXYZ[iSymOp];
                    for (int i = lb; i <= ub; i++) {
                        if (!nativeEnvironmentApproximation && !use[i]) {
                            continue;
                        }
                        final double baseRi = baseRadius[i];
                        assert (baseRi > 0.0);
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        int[] list = neighborLists[iSymOp][i];
                        for (int k : list) {
                            final double baseRk = baseRadius[k];
                            assert (baseRk > 0.0);
                            if (!nativeEnvironmentApproximation && !use[k]) {
                                continue;
                            }
                            if (i != k) {
                                final double xr = xyz[0][k] - xi;
                                final double yr = xyz[1][k] - yi;
                                final double zr = xyz[2][k] - zi;
                                final double r2 = crystal.image(xr, yr, zr);
                                if (r2 > cut2) {
                                    continue;
                                }
                                final double r = sqrt(r2);

                                // Atom i being descreeened by atom k.
                                double scaledRk = baseRk * overlapScale[k];
                                localBorn[i] += integral(r, r2, baseRi, scaledRk);

                                // Atom k being descreeened by atom i.
                                double scaledRi = baseRi * overlapScale[i];
                                localBorn[k] += integral(r, r2, baseRk, scaledRi);
                            } else if (iSymOp > 0) {
                                final double xr = xyz[0][k] - xi;
                                final double yr = xyz[1][k] - yi;
                                final double zr = xyz[2][k] - zi;
                                final double r2 = crystal.image(xr, yr, zr);
                                if (r2 > cut2) {
                                    continue;
                                }
                                final double r = sqrt(r2);

                                // Atom i being descreeened by atom k.
                                double scaledRk = baseRk * overlapScale[k];
                                localBorn[i] += integral(r, r2, baseRi, scaledRk);

                                // For symmetry mates, atom k is not descreeened by atom i.
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Compute the Generalized Kirkwood permanent reaction field in parallel.
     *
     * @since 1.0
     */
    private class PermanentGKFieldRegion extends ParallelRegion {

        private final PermanentGKFieldLoop[] permanentGKFieldLoop;

        private PermanentGKFieldRegion(int nt) {
            permanentGKFieldLoop = new PermanentGKFieldLoop[nt];
            for (int i = 0; i < nt; i++) {
                permanentGKFieldLoop[i] = new PermanentGKFieldLoop();
            }
        }

        @Override
        public void start() {
            for (int i = 0; i < nAtoms; i++) {
                sharedGKField[0].set(i, 0.0);
                sharedGKField[1].set(i, 0.0);
                sharedGKField[2].set(i, 0.0);
            }
        }

        @Override
        public void run() {
            try {
                int threadIndex = getThreadIndex();
                execute(0, nAtoms - 1, permanentGKFieldLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing GK Energy in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute the Generalized Kirkwood permanent reaction field.
         *
         * @since 1.0
         */
        private class PermanentGKFieldLoop extends IntegerForLoop {

            private final double[][] a;
            private final double[] gc;
            private final double[] gux, guy, guz;
            private final double[] gqxx, gqyy, gqzz;
            private final double[] gqxy, gqxz, gqyz;
            private double[] fx_local;
            private double[] fy_local;
            private double[] fz_local;
            private final double[] dx_local;
            private double xi, yi, zi;
            private double ci, uxi, uyi, uzi, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi;
            private double rbi;
            private int iSymm;
            private double[][] transOp;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public PermanentGKFieldLoop() {
                a = new double[4][3];
                gc = new double[11];
                gux = new double[11];
                guy = new double[11];
                guz = new double[11];
                gqxx = new double[11];
                gqyy = new double[11];
                gqzz = new double[11];
                gqxy = new double[11];
                gqxz = new double[11];
                gqyz = new double[11];
                dx_local = new double[3];
                transOp = new double[3][3];
            }

            @Override
            public void start() {
                if (fx_local == null || fx_local.length < maxNumAtoms) {
                    fx_local = new double[maxNumAtoms];
                    fy_local = new double[maxNumAtoms];
                    fz_local = new double[maxNumAtoms];
                }
                fill(fx_local, 0.0);
                fill(fy_local, 0.0);
                fill(fz_local, 0.0);
            }

            @Override
            public void finish() {
                // Reduce the field contributions computed by the current thread into the shared arrays.
                sharedGKField[0].reduce(fx_local, DoubleOp.SUM);
                sharedGKField[1].reduce(fy_local, DoubleOp.SUM);
                sharedGKField[2].reduce(fz_local, DoubleOp.SUM);
            }

            @Override
            public void run(int lb, int ub) {
                int nSymm = crystal.spaceGroup.symOps.size();
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                for (iSymm = 0; iSymm < nSymm; iSymm++) {
                    SymOp symOp = symOps.get(iSymm);
                    crystal.getTransformationOperator(symOp, transOp);
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        xi = x[i];
                        yi = y[i];
                        zi = z[i];
                        double[] multipolei = globalMultipole[0][i];
                        ci = multipolei[t000];
                        uxi = multipolei[t100];
                        uyi = multipolei[t010];
                        uzi = multipolei[t001];
                        qxxi = multipolei[t200] * oneThird;
                        qxyi = multipolei[t110] * oneThird;
                        qxzi = multipolei[t101] * oneThird;
                        qyyi = multipolei[t020] * oneThird;
                        qyzi = multipolei[t011] * oneThird;
                        qzzi = multipolei[t002] * oneThird;
                        rbi = born[i];
                        int[] list = neighborLists[iSymm][i];
                        for (int k : list) {
                            if (!use[k]) {
                                continue;
                            }
                            permanentGKField(i, k);
                        }

                        // Include the self permanent reaction field, which is not in the neighbor list.
                        if (iSymm == 0) {
                            permanentGKField(i, i);
                        }
                    }
                }
            }

            private void permanentGKField(int i, int k) {
                dx_local[0] = sXYZ[iSymm][0][k] - xi;
                dx_local[1] = sXYZ[iSymm][1][k] - yi;
                dx_local[2] = sXYZ[iSymm][2][k] - zi;
                final double r2 = crystal.image(dx_local);
                if (r2 > cut2) {
                    return;
                }
                double xr = dx_local[0];
                double yr = dx_local[1];
                double zr = dx_local[2];
                double xr2 = xr * xr;
                double yr2 = yr * yr;
                double zr2 = zr * zr;
                final double rbk = born[k];
                final double[] multipolek = globalMultipole[iSymm][k];
                final double ck = multipolek[t000];
                final double uxk = multipolek[t100];
                final double uyk = multipolek[t010];
                final double uzk = multipolek[t001];
                final double qxxk = multipolek[t200] * oneThird;
                final double qxyk = multipolek[t110] * oneThird;
                final double qxzk = multipolek[t101] * oneThird;
                final double qyyk = multipolek[t020] * oneThird;
                final double qyzk = multipolek[t011] * oneThird;
                final double qzzk = multipolek[t002] * oneThird;
                final double rb2 = rbi * rbk;
                final double expterm = exp(-r2 / (gkc * rb2));
                final double expc = expterm / gkc;
                final double expc1 = 1.0 - expc;
                final double dexpc = -2.0 / (gkc * rb2);
                final double expcdexpc = -expc * dexpc;
                final double gf2 = 1.0 / (r2 + rb2 * expterm);
                final double gf = sqrt(gf2);
                final double gf3 = gf2 * gf;
                final double gf5 = gf3 * gf2;
                final double gf7 = gf5 * gf2;

                // Reaction potential auxiliary terms.
                a[0][0] = gf;
                a[1][0] = -gf3;
                a[2][0] = 3.0 * gf5;
                a[3][0] = -15.0 * gf7;

                // Reaction potential gradient auxiliary terms.
                a[0][1] = expc1 * a[1][0];
                a[1][1] = expc1 * a[2][0];
                a[2][1] = expc1 * a[3][0];
                // 2nd reaction potential gradient auxiliary terms.
                a[1][2] = expc1 * a[2][1] + expcdexpc * a[2][0];

                // Multiply the potential auxiliary terms by their dielectric functions.
                a[0][1] = fc * a[0][1];
                a[1][0] = fd * a[1][0];
                a[1][1] = fd * a[1][1];
                a[1][2] = fd * a[1][2];
                a[2][0] = fq * a[2][0];
                a[2][1] = fq * a[2][1];

                // Unweighted reaction potential tensor.
                gux[1] = xr * a[1][0];
                guy[1] = yr * a[1][0];
                guz[1] = zr * a[1][0];

                // Unweighted reaction potential gradient tensor.
                gc[2] = xr * a[0][1];
                gc[3] = yr * a[0][1];
                gc[4] = zr * a[0][1];
                gux[2] = a[1][0] + xr2 * a[1][1];
                gux[3] = xr * yr * a[1][1];
                gux[4] = xr * zr * a[1][1];
                guy[2] = gux[3];
                guy[3] = a[1][0] + yr2 * a[1][1];
                guy[4] = yr * zr * a[1][1];
                guz[2] = gux[4];
                guz[3] = guy[4];
                guz[4] = a[1][0] + zr2 * a[1][1];
                gqxx[2] = xr * (2.0 * a[2][0] + xr2 * a[2][1]);
                gqxx[3] = yr * xr2 * a[2][1];
                gqxx[4] = zr * xr2 * a[2][1];
                gqyy[2] = xr * yr2 * a[2][1];
                gqyy[3] = yr * (2.0 * a[2][0] + yr2 * a[2][1]);
                gqyy[4] = zr * yr2 * a[2][1];
                gqzz[2] = xr * zr2 * a[2][1];
                gqzz[3] = yr * zr2 * a[2][1];
                gqzz[4] = zr * (2.0 * a[2][0] + zr2 * a[2][1]);
                gqxy[2] = yr * (a[2][0] + xr2 * a[2][1]);
                gqxy[3] = xr * (a[2][0] + yr2 * a[2][1]);
                gqxy[4] = zr * xr * yr * a[2][1];
                gqxz[2] = zr * (a[2][0] + xr2 * a[2][1]);
                gqxz[3] = gqxy[4];
                gqxz[4] = xr * (a[2][0] + zr2 * a[2][1]);
                gqyz[2] = gqxy[4];
                gqyz[3] = zr * (a[2][0] + yr2 * a[2][1]);
                gqyz[4] = yr * (a[2][0] + zr2 * a[2][1]);

                // Unweighted 2nd reaction potential gradient tensor.
                gux[5] = xr * (3.0 * a[1][1] + xr2 * a[1][2]);
                gux[6] = yr * (a[1][1] + xr2 * a[1][2]);
                gux[7] = zr * (a[1][1] + xr2 * a[1][2]);
                gux[8] = xr * (a[1][1] + yr2 * a[1][2]);
                gux[9] = zr * xr * yr * a[1][2];
                gux[10] = xr * (a[1][1] + zr2 * a[1][2]);
                guy[5] = yr * (a[1][1] + xr2 * a[1][2]);
                guy[6] = xr * (a[1][1] + yr2 * a[1][2]);
                guy[7] = gux[9];
                guy[8] = yr * (3.0 * a[1][1] + yr2 * a[1][2]);
                guy[9] = zr * (a[1][1] + yr2 * a[1][2]);
                guy[10] = yr * (a[1][1] + zr2 * a[1][2]);
                guz[5] = zr * (a[1][1] + xr2 * a[1][2]);
                guz[6] = gux[9];
                guz[7] = xr * (a[1][1] + zr2 * a[1][2]);
                guz[8] = zr * (a[1][1] + yr2 * a[1][2]);
                guz[9] = yr * (a[1][1] + zr2 * a[1][2]);
                guz[10] = zr * (3.0 * a[1][1] + zr2 * a[1][2]);

                // Generalized Kirkwood permanent reaction field.
                double fix = uxk * gux[2] + uyk * gux[3] + uzk * gux[4]
                        + 0.5 * (ck * gux[1] + qxxk * gux[5]
                        + qyyk * gux[8] + qzzk * gux[10]
                        + 2.0 * (qxyk * gux[6] + qxzk * gux[7]
                        + qyzk * gux[9]))
                        + 0.5 * (ck * gc[2] + qxxk * gqxx[2]
                        + qyyk * gqyy[2] + qzzk * gqzz[2]
                        + 2.0 * (qxyk * gqxy[2] + qxzk * gqxz[2]
                        + qyzk * gqyz[2]));
                double fiy = uxk * guy[2] + uyk * guy[3] + uzk * guy[4]
                        + 0.5 * (ck * guy[1] + qxxk * guy[5]
                        + qyyk * guy[8] + qzzk * guy[10]
                        + 2.0 * (qxyk * guy[6] + qxzk * guy[7]
                        + qyzk * guy[9]))
                        + 0.5 * (ck * gc[3] + qxxk * gqxx[3]
                        + qyyk * gqyy[3] + qzzk * gqzz[3]
                        + 2.0 * (qxyk * gqxy[3] + qxzk * gqxz[3]
                        + qyzk * gqyz[3]));
                double fiz = uxk * guz[2] + uyk * guz[3] + uzk * guz[4]
                        + 0.5 * (ck * guz[1] + qxxk * guz[5]
                        + qyyk * guz[8] + qzzk * guz[10]
                        + 2.0 * (qxyk * guz[6] + qxzk * guz[7]
                        + qyzk * guz[9]))
                        + 0.5 * (ck * gc[4] + qxxk * gqxx[4]
                        + qyyk * gqyy[4] + qzzk * gqzz[4]
                        + 2.0 * (qxyk * gqxy[4] + qxzk * gqxz[4]
                        + qyzk * gqyz[4]));
                double fkx = uxi * gux[2] + uyi * gux[3] + uzi * gux[4]
                        - 0.5 * (ci * gux[1] + qxxi * gux[5]
                        + qyyi * gux[8] + qzzi * gux[10]
                        + 2.0 * (qxyi * gux[6] + qxzi * gux[7]
                        + qyzi * gux[9]))
                        - 0.5 * (ci * gc[2] + qxxi * gqxx[2]
                        + qyyi * gqyy[2] + qzzi * gqzz[2]
                        + 2.0 * (qxyi * gqxy[2] + qxzi * gqxz[2]
                        + qyzi * gqyz[2]));
                double fky = uxi * guy[2] + uyi * guy[3] + uzi * guy[4]
                        - 0.5 * (ci * guy[1] + qxxi * guy[5]
                        + qyyi * guy[8] + qzzi * guy[10]
                        + 2.0 * (qxyi * guy[6] + qxzi * guy[7]
                        + qyzi * guy[9]))
                        - 0.5 * (ci * gc[3] + qxxi * gqxx[3]
                        + qyyi * gqyy[3] + qzzi * gqzz[3]
                        + 2.0 * (qxyi * gqxy[3] + qxzi * gqxz[3]
                        + qyzi * gqyz[3]));
                double fkz = uxi * guz[2] + uyi * guz[3] + uzi * guz[4]
                        - 0.5 * (ci * guz[1] + qxxi * guz[5]
                        + qyyi * guz[8] + qzzi * guz[10]
                        + 2.0 * (qxyi * guz[6] + qxzi * guz[7]
                        + qyzi * guz[9]))
                        - 0.5 * (ci * gc[4] + qxxi * gqxx[4]
                        + qyyi * gqyy[4] + qzzi * gqzz[4]
                        + 2.0 * (qxyi * gqxy[4] + qxzi * gqxz[4]
                        + qyzi * gqyz[4]));

                // Scale the self-field by half, such that it sums to one below.
                if (i == k) {
                    fix *= 0.5;
                    fiy *= 0.5;
                    fiz *= 0.5;
                    fkx *= 0.5;
                    fky *= 0.5;
                    fkz *= 0.5;
                }
                fx_local[i] += fix;
                fy_local[i] += fiy;
                fz_local[i] += fiz;

                double xc = fkx;
                double yc = fky;
                double zc = fkz;
                fkx = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
                fky = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
                fkz = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
                fx_local[k] += fkx;
                fy_local[k] += fky;
                fz_local[k] += fkz;
            }
        }
    }

    /**
     * Compute the Generalized Kirkwood induced reaction field in parallel.
     *
     * @since 1.0
     */
    private class InducedGKFieldRegion extends ParallelRegion {

        private final InducedGKFieldLoop[] inducedGKFieldLoop;

        InducedGKFieldRegion(int nt) {
            inducedGKFieldLoop = new InducedGKFieldLoop[nt];
            for (int i = 0; i < nt; i++) {
                inducedGKFieldLoop[i] = new InducedGKFieldLoop();
            }
        }

        @Override
        public void start() {
            for (int i = 0; i < nAtoms; i++) {
                sharedGKField[0].set(i, 0.0);
                sharedGKField[1].set(i, 0.0);
                sharedGKField[2].set(i, 0.0);
                sharedGKFieldCR[0].set(i, 0.0);
                sharedGKFieldCR[1].set(i, 0.0);
                sharedGKFieldCR[2].set(i, 0.0);
            }
        }

        @Override
        public void run() {
            try {
                int threadIndex = getThreadIndex();
                execute(0, nAtoms - 1, inducedGKFieldLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing GK field in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute the Generalized Kirkwood induced reaction field.
         *
         * @since 1.0
         */
        private class InducedGKFieldLoop extends IntegerForLoop {

            private final double[][] a;
            private final double[] gux;
            private final double[] guy;
            private final double[] guz;
            private double[] fx_local;
            private double[] fy_local;
            private double[] fz_local;
            private double[] fxCR_local;
            private double[] fyCR_local;
            private double[] fzCR_local;
            private final double[] dx_local;
            private double xi, yi, zi;
            private double uix, uiy, uiz;
            private double uixCR, uiyCR, uizCR;
            private double rbi;
            private int iSymm;
            private double[][] transOp;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            InducedGKFieldLoop() {
                a = new double[3][2];
                gux = new double[5];
                guy = new double[5];
                guz = new double[5];
                dx_local = new double[3];
                transOp = new double[3][3];
            }

            @Override
            public void start() {
                if (fx_local == null || fx_local.length < maxNumAtoms) {
                    fx_local = new double[maxNumAtoms];
                    fy_local = new double[maxNumAtoms];
                    fz_local = new double[maxNumAtoms];
                    fxCR_local = new double[maxNumAtoms];
                    fyCR_local = new double[maxNumAtoms];
                    fzCR_local = new double[maxNumAtoms];
                }
                fill(fx_local, 0.0);
                fill(fy_local, 0.0);
                fill(fz_local, 0.0);
                fill(fxCR_local, 0.0);
                fill(fyCR_local, 0.0);
                fill(fzCR_local, 0.0);
            }

            @Override
            public void finish() {

                // Reduce the field contributions computed by the current thread into the shared arrays.
                sharedGKField[0].reduce(fx_local, DoubleOp.SUM);
                sharedGKField[1].reduce(fy_local, DoubleOp.SUM);
                sharedGKField[2].reduce(fz_local, DoubleOp.SUM);
                sharedGKFieldCR[0].reduce(fxCR_local, DoubleOp.SUM);
                sharedGKFieldCR[1].reduce(fyCR_local, DoubleOp.SUM);
                sharedGKFieldCR[2].reduce(fzCR_local, DoubleOp.SUM);
            }

            @Override
            public void run(int lb, int ub) {
                int nSymm = crystal.spaceGroup.symOps.size();
                List<SymOp> symOps = crystal.spaceGroup.symOps;
                for (iSymm = 0; iSymm < nSymm; iSymm++) {
                    SymOp symOp = symOps.get(iSymm);
                    crystal.getTransformationOperator(symOp, transOp);
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        xi = x[i];
                        yi = y[i];
                        zi = z[i];
                        uix = inducedDipole[0][i][0];
                        uiy = inducedDipole[0][i][1];
                        uiz = inducedDipole[0][i][2];
                        uixCR = inducedDipoleCR[0][i][0];
                        uiyCR = inducedDipoleCR[0][i][1];
                        uizCR = inducedDipoleCR[0][i][2];
                        rbi = born[i];
                        int[] list = neighborLists[iSymm][i];
                        for (int k : list) {
                            if (!use[k]) {
                                continue;
                            }
                            inducedGKField(i, k);
                        }

                        // Include the self induced reaction field, which is not in the neighbor list.
                        if (iSymm == 0) {
                            inducedGKField(i, i);
                        }
                    }
                }
            }

            private void inducedGKField(int i, int k) {
                dx_local[0] = sXYZ[iSymm][0][k] - xi;
                dx_local[1] = sXYZ[iSymm][1][k] - yi;
                dx_local[2] = sXYZ[iSymm][2][k] - zi;
                final double r2 = crystal.image(dx_local);
                if (r2 > cut2) {
                    return;
                }
                double xr = dx_local[0];
                double yr = dx_local[1];
                double zr = dx_local[2];
                double xr2 = xr * xr;
                double yr2 = yr * yr;
                double zr2 = zr * zr;
                final double ukx = inducedDipole[iSymm][k][0];
                final double uky = inducedDipole[iSymm][k][1];
                final double ukz = inducedDipole[iSymm][k][2];
                final double ukxCR = inducedDipoleCR[iSymm][k][0];
                final double ukyCR = inducedDipoleCR[iSymm][k][1];
                final double ukzCR = inducedDipoleCR[iSymm][k][2];
                final double rbk = born[k];
                final double rb2 = rbi * rbk;
                final double expterm = exp(-r2 / (gkc * rb2));
                final double expc = expterm / gkc;
                final double expc1 = 1.0 - expc;
                final double gf2 = 1.0 / (r2 + rb2 * expterm);
                final double gf = sqrt(gf2);
                final double gf3 = gf2 * gf;
                final double gf5 = gf3 * gf2;

                // Reaction potential auxiliary terms.
                a[1][0] = -gf3;
                a[2][0] = 3.0 * gf5;

                // Reaction potential gradient auxiliary term.
                a[1][1] = expc1 * a[2][0];

                // Multiply the potential auxiliary terms by their dielectric functions.
                a[1][0] = fd * a[1][0];
                a[1][1] = fd * a[1][1];
                a[2][0] = fq * a[2][0];

                // Unweighted reaction potential gradient tensor.
                gux[2] = a[1][0] + xr2 * a[1][1];
                gux[3] = xr * yr * a[1][1];
                gux[4] = xr * zr * a[1][1];
                guy[2] = gux[3];
                guy[3] = a[1][0] + yr2 * a[1][1];
                guy[4] = yr * zr * a[1][1];
                guz[2] = gux[4];
                guz[3] = guy[4];
                guz[4] = a[1][0] + zr2 * a[1][1];

                // Compute the reaction field due to induced dipoles.
                double fix = ukx * gux[2] + uky * guy[2] + ukz * guz[2];
                double fiy = ukx * gux[3] + uky * guy[3] + ukz * guz[3];
                double fiz = ukx * gux[4] + uky * guy[4] + ukz * guz[4];
                double fkx = uix * gux[2] + uiy * guy[2] + uiz * guz[2];
                double fky = uix * gux[3] + uiy * guy[3] + uiz * guz[3];
                double fkz = uix * gux[4] + uiy * guy[4] + uiz * guz[4];
                double fixCR = ukxCR * gux[2] + ukyCR * guy[2] + ukzCR * guz[2];
                double fiyCR = ukxCR * gux[3] + ukyCR * guy[3] + ukzCR * guz[3];
                double fizCR = ukxCR * gux[4] + ukyCR * guy[4] + ukzCR * guz[4];
                double fkxCR = uixCR * gux[2] + uiyCR * guy[2] + uizCR * guz[2];
                double fkyCR = uixCR * gux[3] + uiyCR * guy[3] + uizCR * guz[3];
                double fkzCR = uixCR * gux[4] + uiyCR * guy[4] + uizCR * guz[4];

                // Scale the self-field by half, such that it sums to one below.
                if (i == k) {
                    fix *= 0.5;
                    fiy *= 0.5;
                    fiz *= 0.5;
                    fkx *= 0.5;
                    fky *= 0.5;
                    fkz *= 0.5;
                    fixCR *= 0.5;
                    fiyCR *= 0.5;
                    fizCR *= 0.5;
                    fkxCR *= 0.5;
                    fkyCR *= 0.5;
                    fkzCR *= 0.5;
                }

                fx_local[i] += fix;
                fy_local[i] += fiy;
                fz_local[i] += fiz;
                double xc = fkx;
                double yc = fky;
                double zc = fkz;
                fkx = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
                fky = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
                fkz = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
                fx_local[k] += fkx;
                fy_local[k] += fky;
                fz_local[k] += fkz;

                fxCR_local[i] += fixCR;
                fyCR_local[i] += fiyCR;
                fzCR_local[i] += fizCR;
                xc = fkxCR;
                yc = fkyCR;
                zc = fkzCR;
                fkxCR = xc * transOp[0][0] + yc * transOp[1][0] + zc * transOp[2][0];
                fkyCR = xc * transOp[0][1] + yc * transOp[1][1] + zc * transOp[2][1];
                fkzCR = xc * transOp[0][2] + yc * transOp[1][2] + zc * transOp[2][2];
                fxCR_local[k] += fkxCR;
                fyCR_local[k] += fkyCR;
                fzCR_local[k] += fkzCR;
            }
        }
    }

    /**
     * Compute the Generalized Kirkwood reaction field energy.
     *
     * @since 1.0
     */
    private class GKEnergyRegion extends ParallelRegion {

        private boolean gradient = false;
        private final SharedDouble sharedGKEnergy;
        private final SharedInteger sharedInteractions;
        private final GKEnergyLoop[] gkEnergyLoop;

        GKEnergyRegion(int nt) {
            gkEnergyLoop = new GKEnergyLoop[nt];
            for (int i = 0; i < nt; i++) {
                gkEnergyLoop[i] = new GKEnergyLoop();
            }
            sharedGKEnergy = new SharedDouble();
            sharedInteractions = new SharedInteger();
        }

        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }

        public double getEnergy() {
            return sharedGKEnergy.get();
        }

        public int getInteractions() {
            return sharedInteractions.get();
        }

        @Override
        public void start() {
            sharedGKEnergy.set(0.0);
            sharedInteractions.set(0);
        }

        @Override
        public void run() {
            try {
                int threadIndex = getThreadIndex();
                gkEnergyLoop[threadIndex].setGradient(gradient);
                execute(0, nAtoms - 1, gkEnergyLoop[threadIndex]);
            } catch (Exception e) {
                String message = "Fatal exception computing GK Energy in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class GKEnergyLoop extends IntegerForLoop {

            private final double[][] a;
            private final double[][] b;
            private final double[] gc;
            private final double[] gux;
            private final double[] guy;
            private final double[] guz;
            private final double[] gqxx;
            private final double[] gqyy;
            private final double[] gqzz;
            private final double[] gqxy;
            private final double[] gqxz;
            private final double[] gqyz;
            private double[] gb_local;
            private double[] gbi_local;
            private final double[] dx_local;
            private double[] gX;
            private double[] gY;
            private double[] gZ;
            private double[] tX;
            private double[] tY;
            private double[] tZ;
            private double[] lgX;
            private double[] lgY;
            private double[] lgZ;
            private double[] ltX;
            private double[] ltY;
            private double[] ltZ;
            private double ci, uxi, uyi, uzi, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi;
            private double ck, uxk, uyk, uzk, qxxk, qxyk, qxzk, qyyk, qyzk, qzzk;
            private double dxi, dyi, dzi, pxi, pyi, pzi, sxi, syi, szi;
            private double dxk, dyk, dzk, pxk, pyk, pzk, sxk, syk, szk;
            private double xr, yr, zr, xr2, yr2, zr2, rbi, rbk;
            private double xi, yi, zi;
            private boolean gradient = false;
            private int count;
            private int iSymm;
            private double[][] transOp;
            private double gkEnergy;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            GKEnergyLoop() {
                a = new double[6][4];
                b = new double[5][3];
                gc = new double[31];
                gux = new double[31];
                guy = new double[31];
                guz = new double[31];
                gqxx = new double[31];
                gqyy = new double[31];
                gqzz = new double[31];
                gqxy = new double[31];
                gqxz = new double[31];
                gqyz = new double[31];
                dx_local = new double[3];
                transOp = new double[3][3];
            }

            public void setGradient(boolean gradient) {
                this.gradient = gradient;
            }

            @Override
            public void start() {

                if (gb_local == null || gb_local.length < nAtoms) {
                    gb_local = new double[nAtoms];
                    gbi_local = new double[nAtoms];
                }

                gkEnergy = 0.0;
                count = 0;
                int threadID = getThreadIndex();
                gX = grad[threadID][0];
                gY = grad[threadID][1];
                gZ = grad[threadID][2];
                tX = torque[threadID][0];
                tY = torque[threadID][1];
                tZ = torque[threadID][2];
                if (gradient) {
                    fill(gb_local, 0.0);
                    fill(gbi_local, 0.0);
                }
                if (lambdaTerm) {
                    lgX = lambdaGrad[threadID][0];
                    lgY = lambdaGrad[threadID][1];
                    lgZ = lambdaGrad[threadID][2];
                    ltX = lambdaTorque[threadID][0];
                    ltY = lambdaTorque[threadID][1];
                    ltZ = lambdaTorque[threadID][2];
                }
            }

            @Override
            public void run(int lb, int ub) {
                int nSymm = crystal.spaceGroup.symOps.size();
                List<SymOp> symOps = crystal.spaceGroup.symOps;

                for (iSymm = 0; iSymm < nSymm; iSymm++) {
                    SymOp symOp = symOps.get(iSymm);
                    crystal.getTransformationOperator(symOp, transOp);
                    for (int i = lb; i <= ub; i++) {
                        if (!use[i]) {
                            continue;
                        }
                        xi = x[i];
                        yi = y[i];
                        zi = z[i];
                        final double[] multipolei = globalMultipole[0][i];
                        ci = multipolei[t000];
                        uxi = multipolei[t100];
                        uyi = multipolei[t010];
                        uzi = multipolei[t001];
                        qxxi = multipolei[t200] * oneThird;
                        qxyi = multipolei[t110] * oneThird;
                        qxzi = multipolei[t101] * oneThird;
                        qyyi = multipolei[t020] * oneThird;
                        qyzi = multipolei[t011] * oneThird;
                        qzzi = multipolei[t002] * oneThird;
                        dxi = inducedDipole[0][i][0];
                        dyi = inducedDipole[0][i][1];
                        dzi = inducedDipole[0][i][2];
                        pxi = inducedDipoleCR[0][i][0];
                        pyi = inducedDipoleCR[0][i][1];
                        pzi = inducedDipoleCR[0][i][2];
                        sxi = dxi + pxi;
                        syi = dyi + pyi;
                        szi = dzi + pzi;
                        rbi = born[i];
                        int[] list = neighborLists[iSymm][i];
                        for (int k : list) {
                            if (!use[k]) {
                                continue;
                            }
                            interaction(i, k);

                        }
                        if (iSymm == 0) {
                            // Include self-interactions for the asymmetric unit atoms.
                            interaction(i, i);

                            /*
                              Formula for Born energy approximation for cavitation energy is:
                              e = surfaceTension / 6 * (ri + probe)^2 * (ri/rb)^6.
                              ri is the base atomic radius the atom.
                              rb is Born radius of the atom.
                             */
                            switch (nonPolar) {
                                case BORN_SOLV:
                                case BORN_CAV_DISP:
                                    double r = baseRadius[i] + dOffset + probe;
                                    double ratio = (baseRadius[i] + dOffset) / born[i];
                                    ratio *= ratio;
                                    ratio *= (ratio * ratio);
                                    double saTerm = surfaceTension * r * r * ratio / 6.0;
                                    gkEnergy += saTerm;
                                    gb_local[i] -= 6.0 * saTerm / born[i];
                                    break;
                                default:
                                    break;
                            }
                        }
                    }
                }
            }

            @Override
            public void finish() {
                sharedInteractions.addAndGet(count);
                sharedGKEnergy.addAndGet(gkEnergy);
                if (gradient) {

                    // Reduce the torque contributions computed by the current thread into the shared array.
                    sharedBornGrad.reduce(gb_local, DoubleOp.SUM);
                }
            }

            private void interaction(int i, int k) {
                dx_local[0] = sXYZ[iSymm][0][k] - xi;
                dx_local[1] = sXYZ[iSymm][1][k] - yi;
                dx_local[2] = sXYZ[iSymm][2][k] - zi;
                double r2 = crystal.image(dx_local);
                if (r2 > cut2) {
                    return;
                }
                xr = dx_local[0];
                yr = dx_local[1];
                zr = dx_local[2];
                xr2 = xr * xr;
                yr2 = yr * yr;
                zr2 = zr * zr;
                rbk = born[k];

                final double[] multipolek = globalMultipole[iSymm][k];
                ck = multipolek[t000];
                uxk = multipolek[t100];
                uyk = multipolek[t010];
                uzk = multipolek[t001];
                qxxk = multipolek[t200] * oneThird;
                qxyk = multipolek[t110] * oneThird;
                qxzk = multipolek[t101] * oneThird;
                qyyk = multipolek[t020] * oneThird;
                qyzk = multipolek[t011] * oneThird;
                qzzk = multipolek[t002] * oneThird;
                dxk = inducedDipole[iSymm][k][0];
                dyk = inducedDipole[iSymm][k][1];
                dzk = inducedDipole[iSymm][k][2];
                pxk = inducedDipoleCR[iSymm][k][0];
                pyk = inducedDipoleCR[iSymm][k][1];
                pzk = inducedDipoleCR[iSymm][k][2];
                sxk = dxk + pxk;
                syk = dyk + pyk;
                szk = dzk + pzk;
                final double rb2 = rbi * rbk;
                final double expterm = exp(-r2 / (gkc * rb2));
                final double expc = expterm / gkc;
                final double expc1 = 1.0 - expc;
                final double expcr = r2 * expterm / (gkc * gkc * rb2 * rb2);
                final double dexpc = -2.0 / (gkc * rb2);
                double expcdexpc = -expc * dexpc;
                final double dexpcr = 2.0 / (gkc * rb2 * rb2);
                final double dgfdr = 0.5 * expterm * (1.0 + r2 / (rb2 * gkc));
                final double gf2 = 1.0 / (r2 + rb2 * expterm);
                final double gf = sqrt(gf2);
                final double gf3 = gf2 * gf;
                final double gf5 = gf3 * gf2;
                final double gf7 = gf5 * gf2;
                final double gf9 = gf7 * gf2;
                final double gf11 = gf9 * gf2;

                // Reaction potential auxiliary terms.
                a[0][0] = gf;
                a[1][0] = -gf3;
                a[2][0] = 3.0 * gf5;
                a[3][0] = -15.0 * gf7;
                a[4][0] = 105.0 * gf9;
                a[5][0] = -945.0 * gf11;

                // Reaction potential gradient auxiliary terms.
                a[0][1] = expc1 * a[1][0];
                a[1][1] = expc1 * a[2][0];
                a[2][1] = expc1 * a[3][0];
                a[3][1] = expc1 * a[4][0];
                a[4][1] = expc1 * a[5][0];

                // 2nd reaction potential gradient auxiliary terms.
                a[0][2] = expc1 * a[1][1] + expcdexpc * a[1][0];
                a[1][2] = expc1 * a[2][1] + expcdexpc * a[2][0];
                a[2][2] = expc1 * a[3][1] + expcdexpc * a[3][0];
                a[3][2] = expc1 * a[4][1] + expcdexpc * a[4][0];

                if (gradient) {

                    // 3rd reaction potential gradient auxiliary terms.
                    expcdexpc = 2.0 * expcdexpc;
                    a[0][3] = expc1 * a[1][2] + expcdexpc * a[1][1];
                    a[1][3] = expc1 * a[2][2] + expcdexpc * a[2][1];
                    a[2][3] = expc1 * a[3][2] + expcdexpc * a[3][1];
                    expcdexpc = -expc * dexpc * dexpc;
                    a[0][3] = a[0][3] + expcdexpc * a[1][0];
                    a[1][3] = a[1][3] + expcdexpc * a[2][0];
                    a[2][3] = a[2][3] + expcdexpc * a[3][0];

                    // Born radii derivatives of reaction potential auxiliary terms.
                    b[0][0] = dgfdr * a[1][0];
                    b[1][0] = dgfdr * a[2][0];
                    b[2][0] = dgfdr * a[3][0];
                    b[3][0] = dgfdr * a[4][0];
                    b[4][0] = dgfdr * a[5][0];

                    // Born radii gradients of reaction potential gradient auxiliary terms.
                    b[0][1] = b[1][0] - expcr * a[1][0] - expc * b[1][0];
                    b[1][1] = b[2][0] - expcr * a[2][0] - expc * b[2][0];
                    b[2][1] = b[3][0] - expcr * a[3][0] - expc * b[3][0];
                    b[3][1] = b[4][0] - expcr * a[4][0] - expc * b[4][0];

                    // Born radii derivatives of the 2nd reaction potential gradient auxiliary terms.
                    b[0][2] = b[1][1] - (expcr * (a[1][1] + dexpc * a[1][0])
                            + expc * (b[1][1] + dexpcr * a[1][0] + dexpc * b[1][0]));
                    b[1][2] = b[2][1] - (expcr * (a[2][1] + dexpc * a[2][0])
                            + expc * (b[2][1] + dexpcr * a[2][0] + dexpc * b[2][0]));
                    b[2][2] = b[3][1] - (expcr * (a[3][1] + dexpc * a[3][0])
                            + expc * (b[3][1] + dexpcr * a[3][0] + dexpc * b[3][0]));

                    // Multiply the Born radii auxiliary terms by their dielectric functions.
                    b[0][0] = ELECTRIC * fc * b[0][0];
                    b[0][1] = ELECTRIC * fc * b[0][1];
                    b[0][2] = ELECTRIC * fc * b[0][2];
                    b[1][0] = ELECTRIC * fd * b[1][0];
                    b[1][1] = ELECTRIC * fd * b[1][1];
                    b[1][2] = ELECTRIC * fd * b[1][2];
                    b[2][0] = ELECTRIC * fq * b[2][0];
                    b[2][1] = ELECTRIC * fq * b[2][1];
                    b[2][2] = ELECTRIC * fq * b[2][2];
                }

                // Multiply the potential auxiliary terms by their dielectric functions.
                a[0][0] = ELECTRIC * fc * a[0][0];
                a[0][1] = ELECTRIC * fc * a[0][1];
                a[0][2] = ELECTRIC * fc * a[0][2];
                a[0][3] = ELECTRIC * fc * a[0][3];
                a[1][0] = ELECTRIC * fd * a[1][0];
                a[1][1] = ELECTRIC * fd * a[1][1];
                a[1][2] = ELECTRIC * fd * a[1][2];
                a[1][3] = ELECTRIC * fd * a[1][3];
                a[2][0] = ELECTRIC * fq * a[2][0];
                a[2][1] = ELECTRIC * fq * a[2][1];
                a[2][2] = ELECTRIC * fq * a[2][2];
                a[2][3] = ELECTRIC * fq * a[2][3];

                // Compute the GK tensors required to compute the energy.
                energyTensors();

                // Compute the GK interaction energy.
                double eik = energy(i, k);

                gkEnergy += eik;
                count++;

                if (gradient || lambdaTerm) {
                    // Compute the additional GK tensors required to compute the energy gradient.
                    gradientTensors();

                    // Compute the permanent GK energy gradient.
                    permanentEnergyGradient(i, k);
                    if (polarization != Polarization.NONE) {
                        // Compute the induced GK energy gradient.
                        polarizationEnergyGradient(i, k);
                    }
                }
            }

            private void energyTensors() {
                // Unweighted reaction potential tensor.
                gc[1] = a[0][0];
                gux[1] = xr * a[1][0];
                guy[1] = yr * a[1][0];
                guz[1] = zr * a[1][0];
                gqxx[1] = xr2 * a[2][0];
                gqyy[1] = yr2 * a[2][0];
                gqzz[1] = zr2 * a[2][0];
                gqxy[1] = xr * yr * a[2][0];
                gqxz[1] = xr * zr * a[2][0];
                gqyz[1] = yr * zr * a[2][0];

                // Unweighted reaction potential gradient tensor.
                gc[2] = xr * a[0][1];
                gc[3] = yr * a[0][1];
                gc[4] = zr * a[0][1];
                gux[2] = a[1][0] + xr2 * a[1][1];
                gux[3] = xr * yr * a[1][1];
                gux[4] = xr * zr * a[1][1];
                guy[2] = gux[3];
                guy[3] = a[1][0] + yr2 * a[1][1];
                guy[4] = yr * zr * a[1][1];
                guz[2] = gux[4];
                guz[3] = guy[4];
                guz[4] = a[1][0] + zr2 * a[1][1];
                gqxx[2] = xr * (2.0 * a[2][0] + xr2 * a[2][1]);
                gqxx[3] = yr * xr2 * a[2][1];
                gqxx[4] = zr * xr2 * a[2][1];
                gqyy[2] = xr * yr2 * a[2][1];
                gqyy[3] = yr * (2.0 * a[2][0] + yr2 * a[2][1]);
                gqyy[4] = zr * yr2 * a[2][1];
                gqzz[2] = xr * zr2 * a[2][1];
                gqzz[3] = yr * zr2 * a[2][1];
                gqzz[4] = zr * (2.0 * a[2][0] + zr2 * a[2][1]);
                gqxy[2] = yr * (a[2][0] + xr2 * a[2][1]);
                gqxy[3] = xr * (a[2][0] + yr2 * a[2][1]);
                gqxy[4] = zr * xr * yr * a[2][1];
                gqxz[2] = zr * (a[2][0] + xr2 * a[2][1]);
                gqxz[3] = gqxy[4];
                gqxz[4] = xr * (a[2][0] + zr2 * a[2][1]);
                gqyz[2] = gqxy[4];
                gqyz[3] = zr * (a[2][0] + yr2 * a[2][1]);
                gqyz[4] = yr * (a[2][0] + zr2 * a[2][1]);

                // Unweighted 2nd reaction potential gradient tensor.
                gc[5] = a[0][1] + xr2 * a[0][2];
                gc[6] = xr * yr * a[0][2];
                gc[7] = xr * zr * a[0][2];
                gc[8] = a[0][1] + yr2 * a[0][2];
                gc[9] = yr * zr * a[0][2];
                gc[10] = a[0][1] + zr2 * a[0][2];
                gux[5] = xr * (3.0 * a[1][1] + xr2 * a[1][2]);
                gux[6] = yr * (a[1][1] + xr2 * a[1][2]);
                gux[7] = zr * (a[1][1] + xr2 * a[1][2]);
                gux[8] = xr * (a[1][1] + yr2 * a[1][2]);
                gux[9] = zr * xr * yr * a[1][2];
                gux[10] = xr * (a[1][1] + zr2 * a[1][2]);
                guy[5] = yr * (a[1][1] + xr2 * a[1][2]);
                guy[6] = xr * (a[1][1] + yr2 * a[1][2]);
                guy[7] = gux[9];
                guy[8] = yr * (3.0 * a[1][1] + yr2 * a[1][2]);
                guy[9] = zr * (a[1][1] + yr2 * a[1][2]);
                guy[10] = yr * (a[1][1] + zr2 * a[1][2]);
                guz[5] = zr * (a[1][1] + xr2 * a[1][2]);
                guz[6] = gux[9];
                guz[7] = xr * (a[1][1] + zr2 * a[1][2]);
                guz[8] = zr * (a[1][1] + yr2 * a[1][2]);
                guz[9] = yr * (a[1][1] + zr2 * a[1][2]);
                guz[10] = zr * (3.0 * a[1][1] + zr2 * a[1][2]);
                gqxx[5] = 2.0 * a[2][0] + xr2 * (5.0 * a[2][1] + xr2 * a[2][2]);
                gqxx[6] = yr * xr * (2.0 * a[2][1] + xr2 * a[2][2]);
                gqxx[7] = zr * xr * (2.0 * a[2][1] + xr2 * a[2][2]);
                gqxx[8] = xr2 * (a[2][1] + yr2 * a[2][2]);
                gqxx[9] = zr * yr * xr2 * a[2][2];
                gqxx[10] = xr2 * (a[2][1] + zr2 * a[2][2]);
                gqyy[5] = yr2 * (a[2][1] + xr2 * a[2][2]);
                gqyy[6] = xr * yr * (2.0 * a[2][1] + yr2 * a[2][2]);
                gqyy[7] = xr * zr * yr2 * a[2][2];
                gqyy[8] = 2.0 * a[2][0] + yr2 * (5.0 * a[2][1] + yr2 * a[2][2]);
                gqyy[9] = yr * zr * (2.0 * a[2][1] + yr2 * a[2][2]);
                gqyy[10] = yr2 * (a[2][1] + zr2 * a[2][2]);
                gqzz[5] = zr2 * (a[2][1] + xr2 * a[2][2]);
                gqzz[6] = xr * yr * zr2 * a[2][2];
                gqzz[7] = xr * zr * (2.0 * a[2][1] + zr2 * a[2][2]);
                gqzz[8] = zr2 * (a[2][1] + yr2 * a[2][2]);
                gqzz[9] = yr * zr * (2.0 * a[2][1] + zr2 * a[2][2]);
                gqzz[10] = 2.0 * a[2][0] + zr2 * (5.0 * a[2][1] + zr2 * a[2][2]);
                gqxy[5] = xr * yr * (3.0 * a[2][1] + xr2 * a[2][2]);
                gqxy[6] = a[2][0] + (xr2 + yr2) * a[2][1] + xr2 * yr2 * a[2][2];
                gqxy[7] = zr * yr * (a[2][1] + xr2 * a[2][2]);
                gqxy[8] = xr * yr * (3.0 * a[2][1] + yr2 * a[2][2]);
                gqxy[9] = zr * xr * (a[2][1] + yr2 * a[2][2]);
                gqxy[10] = xr * yr * (a[2][1] + zr2 * a[2][2]);
                gqxz[5] = xr * zr * (3.0 * a[2][1] + xr2 * a[2][2]);
                gqxz[6] = yr * zr * (a[2][1] + xr2 * a[2][2]);
                gqxz[7] = a[2][0] + (xr2 + zr2) * a[2][1] + xr2 * zr2 * a[2][2];
                gqxz[8] = xr * zr * (a[2][1] + yr2 * a[2][2]);
                gqxz[9] = xr * yr * (a[2][1] + zr2 * a[2][2]);
                gqxz[10] = xr * zr * (3.0 * a[2][1] + zr2 * a[2][2]);
                gqyz[5] = zr * yr * (a[2][1] + xr2 * a[2][2]);
                gqyz[6] = xr * zr * (a[2][1] + yr2 * a[2][2]);
                gqyz[7] = xr * yr * (a[2][1] + zr2 * a[2][2]);
                gqyz[8] = yr * zr * (3.0 * a[2][1] + yr2 * a[2][2]);
                gqyz[9] = a[2][0] + (yr2 + zr2) * a[2][1] + yr2 * zr2 * a[2][2];
                gqyz[10] = yr * zr * (3.0 * a[2][1] + zr2 * a[2][2]);
            }

            private double energy(int i, int k) {
                // Electrostatic solvation energy of the permanent multipoles in their own GK reaction potential.
                double esym = ci * ck * gc[1]
                        - (uxi * (uxk * gux[2] + uyk * guy[2] + uzk * guz[2])
                        + uyi * (uxk * gux[3] + uyk * guy[3] + uzk * guz[3])
                        + uzi * (uxk * gux[4] + uyk * guy[4] + uzk * guz[4]));
                double ewi = ci * (uxk * gc[2] + uyk * gc[3] + uzk * gc[4])
                        - ck * (uxi * gux[1] + uyi * guy[1] + uzi * guz[1])
                        + ci * (qxxk * gc[5] + qyyk * gc[8] + qzzk * gc[10]
                        + 2.0 * (qxyk * gc[6] + qxzk * gc[7] + qyzk * gc[9]))
                        + ck * (qxxi * gqxx[1] + qyyi * gqyy[1] + qzzi * gqzz[1]
                        + 2.0 * (qxyi * gqxy[1] + qxzi * gqxz[1] + qyzi * gqyz[1]))
                        - uxi * (qxxk * gux[5] + qyyk * gux[8] + qzzk * gux[10]
                        + 2.0 * (qxyk * gux[6] + qxzk * gux[7] + qyzk * gux[9]))
                        - uyi * (qxxk * guy[5] + qyyk * guy[8] + qzzk * guy[10]
                        + 2.0 * (qxyk * guy[6] + qxzk * guy[7] + qyzk * guy[9]))
                        - uzi * (qxxk * guz[5] + qyyk * guz[8] + qzzk * guz[10]
                        + 2.0 * (qxyk * guz[6] + qxzk * guz[7] + qyzk * guz[9]))
                        + uxk * (qxxi * gqxx[2] + qyyi * gqyy[2] + qzzi * gqzz[2]
                        + 2.0 * (qxyi * gqxy[2] + qxzi * gqxz[2] + qyzi * gqyz[2]))
                        + uyk * (qxxi * gqxx[3] + qyyi * gqyy[3] + qzzi * gqzz[3]
                        + 2.0 * (qxyi * gqxy[3] + qxzi * gqxz[3] + qyzi * gqyz[3]))
                        + uzk * (qxxi * gqxx[4] + qyyi * gqyy[4] + qzzi * gqzz[4]
                        + 2.0 * (qxyi * gqxy[4] + qxzi * gqxz[4] + qyzi * gqyz[4]))
                        + qxxi * (qxxk * gqxx[5] + qyyk * gqxx[8] + qzzk * gqxx[10]
                        + 2.0 * (qxyk * gqxx[6] + qxzk * gqxx[7] + qyzk * gqxx[9]))
                        + qyyi * (qxxk * gqyy[5] + qyyk * gqyy[8] + qzzk * gqyy[10]
                        + 2.0 * (qxyk * gqyy[6] + qxzk * gqyy[7] + qyzk * gqyy[9]))
                        + qzzi * (qxxk * gqzz[5] + qyyk * gqzz[8] + qzzk * gqzz[10]
                        + 2.0 * (qxyk * gqzz[6] + qxzk * gqzz[7] + qyzk * gqzz[9]))
                        + 2.0 * (qxyi * (qxxk * gqxy[5] + qyyk * gqxy[8] + qzzk * gqxy[10]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxy[7] + qyzk * gqxy[9]))
                        + qxzi * (qxxk * gqxz[5] + qyyk * gqxz[8] + qzzk * gqxz[10]
                        + 2.0 * (qxyk * gqxz[6] + qxzk * gqxz[7] + qyzk * gqxz[9]))
                        + qyzi * (qxxk * gqyz[5] + qyyk * gqyz[8] + qzzk * gqyz[10]
                        + 2.0 * (qxyk * gqyz[6] + qxzk * gqyz[7] + qyzk * gqyz[9])));
                double ewk = ci * (uxk * gux[1] + uyk * guy[1] + uzk * guz[1])
                        - ck * (uxi * gc[2] + uyi * gc[3] + uzi * gc[4])
                        + ci * (qxxk * gqxx[1] + qyyk * gqyy[1] + qzzk * gqzz[1]
                        + 2.0 * (qxyk * gqxy[1] + qxzk * gqxz[1] + qyzk * gqyz[1]))
                        + ck * (qxxi * gc[5] + qyyi * gc[8] + qzzi * gc[10]
                        + 2.0 * (qxyi * gc[6] + qxzi * gc[7] + qyzi * gc[9]))
                        - uxi * (qxxk * gqxx[2] + qyyk * gqyy[2] + qzzk * gqzz[2]
                        + 2.0 * (qxyk * gqxy[2] + qxzk * gqxz[2] + qyzk * gqyz[2]))
                        - uyi * (qxxk * gqxx[3] + qyyk * gqyy[3] + qzzk * gqzz[3]
                        + 2.0 * (qxyk * gqxy[3] + qxzk * gqxz[3] + qyzk * gqyz[3]))
                        - uzi * (qxxk * gqxx[4] + qyyk * gqyy[4] + qzzk * gqzz[4]
                        + 2.0 * (qxyk * gqxy[4] + qxzk * gqxz[4] + qyzk * gqyz[4]))
                        + uxk * (qxxi * gux[5] + qyyi * gux[8] + qzzi * gux[10]
                        + 2.0 * (qxyi * gux[6] + qxzi * gux[7] + qyzi * gux[9]))
                        + uyk * (qxxi * guy[5] + qyyi * guy[8] + qzzi * guy[10]
                        + 2.0 * (qxyi * guy[6] + qxzi * guy[7] + qyzi * guy[9]))
                        + uzk * (qxxi * guz[5] + qyyi * guz[8] + qzzi * guz[10]
                        + 2.0 * (qxyi * guz[6] + qxzi * guz[7] + qyzi * guz[9]))
                        + qxxi * (qxxk * gqxx[5] + qyyk * gqyy[5] + qzzk * gqzz[5]
                        + 2.0 * (qxyk * gqxy[5] + qxzk * gqxz[5] + qyzk * gqyz[5]))
                        + qyyi * (qxxk * gqxx[8] + qyyk * gqyy[8] + qzzk * gqzz[8]
                        + 2.0 * (qxyk * gqxy[8] + qxzk * gqxz[8] + qyzk * gqyz[8]))
                        + qzzi * (qxxk * gqxx[10] + qyyk * gqyy[10] + qzzk * gqzz[10]
                        + 2.0 * (qxyk * gqxy[10] + qxzk * gqxz[10] + qyzk * gqyz[10]))
                        + 2.0 * (qxyi * (qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]))
                        + qxzi * (qxxk * gqxx[7] + qyyk * gqyy[7] + qzzk * gqzz[7]
                        + 2.0 * (qxyk * gqxy[7] + qxzk * gqxz[7] + qyzk * gqyz[7]))
                        + qyzi * (qxxk * gqxx[9] + qyyk * gqyy[9] + qzzk * gqzz[9]
                        + 2.0 * (qxyk * gqxy[9] + qxzk * gqxz[9] + qyzk * gqyz[9])));
                double e = esym + 0.5 * (ewi + ewk);
                double ei = 0.0;
                if (polarization != Polarization.NONE) {

                    // Electrostatic solvation energy of the permanent multipoles in the 
                    // GK reaction potential of the induced dipoles.
                    double esymi = -uxi * (dxk * gux[2] + dyk * guy[2] + dzk * guz[2])
                            - uyi * (dxk * gux[3] + dyk * guy[3] + dzk * guz[3])
                            - uzi * (dxk * gux[4] + dyk * guy[4] + dzk * guz[4])
                            - uxk * (dxi * gux[2] + dyi * guy[2] + dzi * guz[2])
                            - uyk * (dxi * gux[3] + dyi * guy[3] + dzi * guz[3])
                            - uzk * (dxi * gux[4] + dyi * guy[4] + dzi * guz[4]);
                    double ewii = ci * (dxk * gc[2] + dyk * gc[3] + dzk * gc[4])
                            - ck * (dxi * gux[1] + dyi * guy[1] + dzi * guz[1])
                            - dxi * (qxxk * gux[5] + qyyk * gux[8] + qzzk * gux[10]
                            + 2.0 * (qxyk * gux[6] + qxzk * gux[7] + qyzk * gux[9]))
                            - dyi * (qxxk * guy[5] + qyyk * guy[8] + qzzk * guy[10]
                            + 2.0 * (qxyk * guy[6] + qxzk * guy[7] + qyzk * guy[9]))
                            - dzi * (qxxk * guz[5] + qyyk * guz[8] + qzzk * guz[10]
                            + 2.0 * (qxyk * guz[6] + qxzk * guz[7] + qyzk * guz[9]))
                            + dxk * (qxxi * gqxx[2] + qyyi * gqyy[2] + qzzi * gqzz[2]
                            + 2.0 * (qxyi * gqxy[2] + qxzi * gqxz[2] + qyzi * gqyz[2]))
                            + dyk * (qxxi * gqxx[3] + qyyi * gqyy[3] + qzzi * gqzz[3]
                            + 2.0 * (qxyi * gqxy[3] + qxzi * gqxz[3] + qyzi * gqyz[3]))
                            + dzk * (qxxi * gqxx[4] + qyyi * gqyy[4] + qzzi * gqzz[4]
                            + 2.0 * (qxyi * gqxy[4] + qxzi * gqxz[4] + qyzi * gqyz[4]));
                    double ewki = ci * (dxk * gux[1] + dyk * guy[1] + dzk * guz[1])
                            - ck * (dxi * gc[2] + dyi * gc[3] + dzi * gc[4])
                            - dxi * (qxxk * gqxx[2] + qyyk * gqyy[2] + qzzk * gqzz[2]
                            + 2.0 * (qxyk * gqxy[2] + qxzk * gqxz[2] + qyzk * gqyz[2]))
                            - dyi * (qxxk * gqxx[3] + qyyk * gqyy[3] + qzzk * gqzz[3]
                            + 2.0 * (qxyk * gqxy[3] + qxzk * gqxz[3] + qyzk * gqyz[3]))
                            - dzi * (qxxk * gqxx[4] + qyyk * gqyy[4] + qzzk * gqzz[4]
                            + 2.0 * (qxyk * gqxy[4] + qxzk * gqxz[4] + qyzk * gqyz[4]))
                            + dxk * (qxxi * gux[5] + qyyi * gux[8] + qzzi * gux[10]
                            + 2.0 * (qxyi * gux[6] + qxzi * gux[7] + qyzi * gux[9]))
                            + dyk * (qxxi * guy[5] + qyyi * guy[8] + qzzi * guy[10]
                            + 2.0 * (qxyi * guy[6] + qxzi * guy[7] + qyzi * guy[9]))
                            + dzk * (qxxi * guz[5] + qyyi * guz[8] + qzzi * guz[10]
                            + 2.0 * (qxyi * guz[6] + qxzi * guz[7] + qyzi * guz[9]));
                    ei = 0.5 * (esymi + 0.5 * (ewii + ewki));
                }
                if (i == k) {
                    e *= 0.5;
                    ei *= 0.5;
                }
                return e + ei;
            }

            private void gradientTensors() {
                // Born radii gradients of unweighted reaction potential tensor.
                gc[21] = b[0][0];
                gux[21] = xr * b[1][0];
                guy[21] = yr * b[1][0];
                guz[21] = zr * b[1][0];
                gqxx[21] = xr2 * b[2][0];
                gqyy[21] = yr2 * b[2][0];
                gqzz[21] = zr2 * b[2][0];
                gqxy[21] = xr * yr * b[2][0];
                gqxz[21] = xr * zr * b[2][0];
                gqyz[21] = yr * zr * b[2][0];

                // Born gradients of the unweighted reaction potential gradient tensor
                gc[22] = xr * b[0][1];
                gc[23] = yr * b[0][1];
                gc[24] = zr * b[0][1];
                gux[22] = b[1][0] + xr2 * b[1][1];
                gux[23] = xr * yr * b[1][1];
                gux[24] = xr * zr * b[1][1];
                guy[22] = gux[23];
                guy[23] = b[1][0] + yr2 * b[1][1];
                guy[24] = yr * zr * b[1][1];
                guz[22] = gux[24];
                guz[23] = guy[24];
                guz[24] = b[1][0] + zr2 * b[1][1];
                gqxx[22] = xr * (2.0 * b[2][0] + xr2 * b[2][1]);
                gqxx[23] = yr * xr2 * b[2][1];
                gqxx[24] = zr * xr2 * b[2][1];
                gqyy[22] = xr * yr2 * b[2][1];
                gqyy[23] = yr * (2.0 * b[2][0] + yr2 * b[2][1]);
                gqyy[24] = zr * yr2 * b[2][1];
                gqzz[22] = xr * zr2 * b[2][1];
                gqzz[23] = yr * zr2 * b[2][1];
                gqzz[24] = zr * (2.0 * b[2][0] + zr2 * b[2][1]);
                gqxy[22] = yr * (b[2][0] + xr2 * b[2][1]);
                gqxy[23] = xr * (b[2][0] + yr2 * b[2][1]);
                gqxy[24] = zr * xr * yr * b[2][1];
                gqxz[22] = zr * (b[2][0] + xr2 * b[2][1]);
                gqxz[23] = gqxy[24];
                gqxz[24] = xr * (b[2][0] + zr2 * b[2][1]);
                gqyz[22] = gqxy[24];
                gqyz[23] = zr * (b[2][0] + yr2 * b[2][1]);
                gqyz[24] = yr * (b[2][0] + zr2 * b[2][1]);

                // Born radii derivatives of the unweighted 2nd reaction potential gradient tensor.
                gc[25] = b[0][1] + xr2 * b[0][2];
                gc[26] = xr * yr * b[0][2];
                gc[27] = xr * zr * b[0][2];
                gc[28] = b[0][1] + yr2 * b[0][2];
                gc[29] = yr * zr * b[0][2];
                gc[30] = b[0][1] + zr2 * b[0][2];
                gux[25] = xr * (3.0 * b[1][1] + xr2 * b[1][2]);
                gux[26] = yr * (b[1][1] + xr2 * b[1][2]);
                gux[27] = zr * (b[1][1] + xr2 * b[1][2]);
                gux[28] = xr * (b[1][1] + yr2 * b[1][2]);
                gux[29] = zr * xr * yr * b[1][2];
                gux[30] = xr * (b[1][1] + zr2 * b[1][2]);
                guy[25] = yr * (b[1][1] + xr2 * b[1][2]);
                guy[26] = xr * (b[1][1] + yr2 * b[1][2]);
                guy[27] = gux[29];
                guy[28] = yr * (3.0 * b[1][1] + yr2 * b[1][2]);
                guy[29] = zr * (b[1][1] + yr2 * b[1][2]);
                guy[30] = yr * (b[1][1] + zr2 * b[1][2]);
                guz[25] = zr * (b[1][1] + xr2 * b[1][2]);
                guz[26] = gux[29];
                guz[27] = xr * (b[1][1] + zr2 * b[1][2]);
                guz[28] = zr * (b[1][1] + yr2 * b[1][2]);
                guz[29] = yr * (b[1][1] + zr2 * b[1][2]);
                guz[30] = zr * (3.0 * b[1][1] + zr2 * b[1][2]);
                gqxx[25] = 2.0 * b[2][0] + xr2 * (5.0 * b[2][1] + xr2 * b[2][2]);
                gqxx[26] = yr * xr * (2.0 * b[2][1] + xr2 * b[2][2]);
                gqxx[27] = zr * xr * (2.0 * b[2][1] + xr2 * b[2][2]);
                gqxx[28] = xr2 * (b[2][1] + yr2 * b[2][2]);
                gqxx[29] = zr * yr * xr2 * b[2][2];
                gqxx[30] = xr2 * (b[2][1] + zr2 * b[2][2]);
                gqyy[25] = yr2 * (b[2][1] + xr2 * b[2][2]);
                gqyy[26] = xr * yr * (2.0 * b[2][1] + yr2 * b[2][2]);
                gqyy[27] = xr * zr * yr2 * b[2][2];
                gqyy[28] = 2.0 * b[2][0] + yr2 * (5.0 * b[2][1] + yr2 * b[2][2]);
                gqyy[29] = yr * zr * (2.0 * b[2][1] + yr2 * b[2][2]);
                gqyy[30] = yr2 * (b[2][1] + zr2 * b[2][2]);
                gqzz[25] = zr2 * (b[2][1] + xr2 * b[2][2]);
                gqzz[26] = xr * yr * zr2 * b[2][2];
                gqzz[27] = xr * zr * (2.0 * b[2][1] + zr2 * b[2][2]);
                gqzz[28] = zr2 * (b[2][1] + yr2 * b[2][2]);
                gqzz[29] = yr * zr * (2.0 * b[2][1] + zr2 * b[2][2]);
                gqzz[30] = 2.0 * b[2][0] + zr2 * (5.0 * b[2][1] + zr2 * b[2][2]);
                gqxy[25] = xr * yr * (3.0 * b[2][1] + xr2 * b[2][2]);
                gqxy[26] = b[2][0] + (xr2 + yr2) * b[2][1] + xr2 * yr2 * b[2][2];
                gqxy[27] = zr * yr * (b[2][1] + xr2 * b[2][2]);
                gqxy[28] = xr * yr * (3.0 * b[2][1] + yr2 * b[2][2]);
                gqxy[29] = zr * xr * (b[2][1] + yr2 * b[2][2]);
                gqxy[30] = xr * yr * (b[2][1] + zr2 * b[2][2]);
                gqxz[25] = xr * zr * (3.0 * b[2][1] + xr2 * b[2][2]);
                gqxz[26] = yr * zr * (b[2][1] + xr2 * b[2][2]);
                gqxz[27] = b[2][0] + (xr2 + zr2) * b[2][1] + xr2 * zr2 * b[2][2];
                gqxz[28] = xr * zr * (b[2][1] + yr2 * b[2][2]);
                gqxz[29] = xr * yr * (b[2][1] + zr2 * b[2][2]);
                gqxz[30] = xr * zr * (3.0 * b[2][1] + zr2 * b[2][2]);
                gqyz[25] = zr * yr * (b[2][1] + xr2 * b[2][2]);
                gqyz[26] = xr * zr * (b[2][1] + yr2 * b[2][2]);
                gqyz[27] = xr * yr * (b[2][1] + zr2 * b[2][2]);
                gqyz[28] = yr * zr * (3.0 * b[2][1] + yr2 * b[2][2]);
                gqyz[29] = b[2][0] + (yr2 + zr2) * b[2][1] + yr2 * zr2 * b[2][2];
                gqyz[30] = yr * zr * (3.0 * b[2][1] + zr2 * b[2][2]);

                // Unweighted 3rd reaction potential gradient tensor.
                gc[11] = xr * (3.0 * a[0][2] + xr2 * a[0][3]);
                gc[12] = yr * (a[0][2] + xr2 * a[0][3]);
                gc[13] = zr * (a[0][2] + xr2 * a[0][3]);
                gc[14] = xr * (a[0][2] + yr2 * a[0][3]);
                gc[15] = xr * yr * zr * a[0][3];
                gc[16] = xr * (a[0][2] + zr2 * a[0][3]);
                gc[17] = yr * (3.0 * a[0][2] + yr2 * a[0][3]);
                gc[18] = zr * (a[0][2] + yr2 * a[0][3]);
                gc[19] = yr * (a[0][2] + zr2 * a[0][3]);
                gc[20] = zr * (3.0 * a[0][2] + zr2 * a[0][3]);

                gux[11] = 3.0 * a[1][1] + xr2 * (6.0 * a[1][2] + xr2 * a[1][3]);
                gux[12] = xr * yr * (3.0 * a[1][2] + xr2 * a[1][3]);
                gux[13] = xr * zr * (3.0 * a[1][2] + xr2 * a[1][3]);
                gux[14] = a[1][1] + (xr2 + yr2) * a[1][2] + xr2 * yr2 * a[1][3];
                gux[15] = yr * zr * (a[1][2] + xr2 * a[1][3]);
                gux[16] = a[1][1] + (xr2 + zr2) * a[1][2] + xr2 * zr2 * a[1][3];
                gux[17] = xr * yr * (3.0 * a[1][2] + yr2 * a[1][3]);
                gux[18] = xr * zr * (a[1][2] + yr2 * a[1][3]);
                gux[19] = xr * yr * (a[1][2] + zr2 * a[1][3]);
                gux[20] = xr * zr * (3.0 * a[1][2] + zr2 * a[1][3]);

                guy[11] = gux[12];
                guy[12] = gux[14];
                guy[13] = gux[15];
                guy[14] = gux[17];
                guy[15] = gux[18];
                guy[16] = gux[19];
                guy[17] = 3.0 * a[1][1] + yr2 * (6.0 * a[1][2] + yr2 * a[1][3]);
                guy[18] = yr * zr * (3.0 * a[1][2] + yr2 * a[1][3]);
                guy[19] = a[1][1] + (yr2 + zr2) * a[1][2] + yr2 * zr2 * a[1][3];
                guy[20] = yr * zr * (3.0 * a[1][2] + zr2 * a[1][3]);

                guz[11] = gux[13];
                guz[12] = gux[15];
                guz[13] = gux[16];
                guz[14] = gux[18];
                guz[15] = gux[19];
                guz[16] = gux[20];
                guz[17] = guy[18];
                guz[18] = guy[19];
                guz[19] = guy[20];
                guz[20] = 3.0 * a[1][1] + zr2 * (6.0 * a[1][2] + zr2 * a[1][3]);

                gqxx[11] = xr * (12.0 * a[2][1] + xr2 * (9.0 * a[2][2] + xr2 * a[2][3]));
                gqxx[12] = yr * (2.0 * a[2][1] + xr2 * (5.0 * a[2][2] + xr2 * a[2][3]));
                gqxx[13] = zr * (2.0 * a[2][1] + xr2 * (5.0 * a[2][2] + xr2 * a[2][3]));
                gqxx[14] = xr * (2.0 * a[2][1] + yr2 * 2.0 * a[2][2] + xr2 * (a[2][2] + yr2 * a[2][3]));
                gqxx[15] = xr * yr * zr * (2.0 * a[2][2] + xr2 * a[2][3]);
                gqxx[16] = xr * (2.0 * a[2][1] + zr2 * 2.0 * a[2][2] + xr2 * (a[2][2] + zr2 * a[2][3]));
                gqxx[17] = yr * xr2 * (3.0 * a[2][2] + yr2 * a[2][3]);
                gqxx[18] = zr * xr2 * (a[2][2] + yr2 * a[2][3]);
                gqxx[19] = yr * xr2 * (a[2][2] + zr2 * a[2][3]);
                gqxx[20] = zr * xr2 * (3.0 * a[2][2] + zr2 * a[2][3]);

                gqxy[11] = yr * (3.0 * a[2][1] + xr2 * (6.0 * a[2][2] + xr2 * a[2][3]));
                gqxy[12] = xr * (3.0 * (a[2][1] + yr2 * a[2][2]) + xr2 * (a[2][2] + yr2 * a[2][3]));
                gqxy[13] = xr * yr * zr * (3.0 * a[2][2] + xr2 * a[2][3]);
                gqxy[14] = yr * (3.0 * (a[2][1] + xr2 * a[2][2]) + yr2 * (a[2][2] + xr2 * a[2][3]));
                gqxy[15] = zr * (a[2][1] + (yr2 + xr2) * a[2][2] + yr2 * xr2 * a[2][3]);
                gqxy[16] = yr * (a[2][1] + (xr2 + zr2) * a[2][2] + xr2 * zr2 * a[2][3]);
                gqxy[17] = xr * (3.0 * (a[2][1] + yr2 * a[2][2]) + yr2 * (3.0 * a[2][2] + yr2 * a[2][3]));
                gqxy[18] = xr * yr * zr * (3.0 * a[2][2] + yr2 * a[2][3]);
                gqxy[19] = xr * (a[2][1] + (yr2 + zr2) * a[2][2] + yr2 * zr2 * a[2][3]);
                gqxy[20] = xr * yr * zr * (3.0 * a[2][2] + zr2 * a[2][3]);

                gqxz[11] = zr * (3.0 * a[2][1] + xr2 * (6.0 * a[2][2] + xr2 * a[2][3]));
                gqxz[12] = xr * yr * zr * (3.0 * a[2][2] + xr2 * a[2][3]);
                gqxz[13] = xr * (3.0 * (a[2][1] + zr2 * a[2][2]) + xr2 * (a[2][2] + zr2 * a[2][3]));
                gqxz[14] = zr * (a[2][1] + (xr2 + yr2) * a[2][2] + xr2 * yr2 * a[2][3]);
                gqxz[15] = yr * (a[2][1] + (xr2 + zr2) * a[2][2] + zr2 * xr2 * a[2][3]);
                gqxz[16] = zr * (3.0 * (a[2][1] + xr2 * a[2][2]) + zr2 * (a[2][2] + xr2 * a[2][3]));
                gqxz[17] = xr * yr * zr * (3.0 * a[2][2] + yr2 * a[2][3]);
                gqxz[18] = xr * (a[2][1] + (zr2 + yr2) * a[2][2] + zr2 * yr2 * a[2][3]);
                gqxz[19] = xr * yr * zr * (3.0 * a[2][2] + zr2 * a[2][3]);
                gqxz[20] = xr * (3.0 * a[2][1] + zr2 * (6.0 * a[2][2] + zr2 * a[2][3]));

                gqyy[11] = xr * yr2 * (3.0 * a[2][2] + xr2 * a[2][3]);
                gqyy[12] = yr * (2.0 * a[2][1] + xr2 * 2.0 * a[2][2] + yr2 * (a[2][2] + xr2 * a[2][3]));
                gqyy[13] = zr * yr2 * (a[2][2] + xr2 * a[2][3]);
                gqyy[14] = xr * (2.0 * a[2][1] + yr2 * (5.0 * a[2][2] + yr2 * a[2][3]));
                gqyy[15] = xr * yr * zr * (2.0 * a[2][2] + yr2 * a[2][3]);
                gqyy[16] = xr * yr2 * (a[2][2] + zr2 * a[2][3]);
                gqyy[17] = yr * (12.0 * a[2][1] + yr2 * (9.0 * a[2][2] + yr2 * a[2][3]));
                gqyy[18] = zr * (2.0 * a[2][1] + yr2 * (5.0 * a[2][2] + yr2 * a[2][3]));
                gqyy[19] = yr * (2.0 * a[2][1] + zr2 * 2.0 * a[2][2] + yr2 * (a[2][2] + zr2 * a[2][3]));
                gqyy[20] = zr * yr2 * (3.0 * a[2][2] + zr2 * a[2][3]);

                gqyz[11] = xr * yr * zr * (3.0 * a[2][2] + xr2 * a[2][3]);
                gqyz[12] = zr * (a[2][1] + (xr2 + yr2) * a[2][2] + xr2 * yr2 * a[2][3]);
                gqyz[13] = yr * (a[2][1] + (xr2 + zr2) * a[2][2] + xr2 * zr2 * a[2][3]);
                gqyz[14] = xr * yr * zr * (3.0 * a[2][2] + yr2 * a[2][3]);
                gqyz[15] = xr * (a[2][1] + (yr2 + zr2) * a[2][2] + yr2 * zr2 * a[2][3]);
                gqyz[16] = xr * yr * zr * (3.0 * a[2][2] + zr2 * a[2][3]);
                gqyz[17] = zr * (3.0 * a[2][1] + yr2 * (6.0 * a[2][2] + yr2 * a[2][3]));
                gqyz[18] = yr * (3.0 * (a[2][1] + zr2 * a[2][2]) + yr2 * (a[2][2] + zr2 * a[2][3]));
                gqyz[19] = zr * (3.0 * (a[2][1] + yr2 * a[2][2]) + zr2 * (a[2][2] + yr2 * a[2][3]));
                gqyz[20] = yr * (3.0 * a[2][1] + zr2 * (6.0 * a[2][2] + zr2 * a[2][3]));

                gqzz[11] = xr * zr2 * (3.0 * a[2][2] + xr2 * a[2][3]);
                gqzz[12] = yr * (zr2 * a[2][2] + xr2 * (zr2 * a[2][3]));
                gqzz[13] = zr * (2.0 * a[2][1] + xr2 * 2.0 * a[2][2] + zr2 * (a[2][2] + xr2 * a[2][3]));
                gqzz[14] = xr * zr2 * (a[2][2] + yr2 * a[2][3]);
                gqzz[15] = xr * yr * zr * (2.0 * a[2][2] + zr2 * a[2][3]);
                gqzz[16] = xr * (2.0 * a[2][1] + zr2 * (5.0 * a[2][2] + zr2 * a[2][3]));
                gqzz[17] = yr * zr2 * (3.0 * a[2][2] + yr2 * a[2][3]);
                gqzz[18] = zr * (2.0 * a[2][1] + yr2 * 2.0 * a[2][2] + zr2 * (a[2][2] + yr2 * a[2][3]));
                gqzz[19] = yr * (2.0 * a[2][1] + zr2 * (5.0 * a[2][2] + zr2 * a[2][3]));
                gqzz[20] = zr * (12.0 * a[2][1] + zr2 * (9.0 * a[2][2] + zr2 * a[2][3]));
            }

            private void permanentEnergyGradient(int i, int k) {
                final double desymdr = ci * ck * gc[21]
                        - (uxi * (uxk * gux[22] + uyk * guy[22] + uzk * guz[22])
                        + uyi * (uxk * gux[23] + uyk * guy[23] + uzk * guz[23])
                        + uzi * (uxk * gux[24] + uyk * guy[24] + uzk * guz[24]));
                final double dewidr = ci * (uxk * gc[22] + uyk * gc[23] + uzk * gc[24])
                        - ck * (uxi * gux[21] + uyi * guy[21] + uzi * guz[21])
                        + ci * (qxxk * gc[25] + qyyk * gc[28] + qzzk * gc[30]
                        + 2.0 * (qxyk * gc[26] + qxzk * gc[27] + qyzk * gc[29]))
                        + ck * (qxxi * gqxx[21] + qyyi * gqyy[21] + qzzi * gqzz[21]
                        + 2.0 * (qxyi * gqxy[21] + qxzi * gqxz[21] + qyzi * gqyz[21]))
                        - uxi * (qxxk * gux[25] + qyyk * gux[28] + qzzk * gux[30]
                        + 2.0 * (qxyk * gux[26] + qxzk * gux[27] + qyzk * gux[29]))
                        - uyi * (qxxk * guy[25] + qyyk * guy[28] + qzzk * guy[30]
                        + 2.0 * (qxyk * guy[26] + qxzk * guy[27] + qyzk * guy[29]))
                        - uzi * (qxxk * guz[25] + qyyk * guz[28] + qzzk * guz[30]
                        + 2.0 * (qxyk * guz[26] + qxzk * guz[27] + qyzk * guz[29]))
                        + uxk * (qxxi * gqxx[22] + qyyi * gqyy[22] + qzzi * gqzz[22]
                        + 2.0 * (qxyi * gqxy[22] + qxzi * gqxz[22] + qyzi * gqyz[22]))
                        + uyk * (qxxi * gqxx[23] + qyyi * gqyy[23] + qzzi * gqzz[23]
                        + 2.0 * (qxyi * gqxy[23] + qxzi * gqxz[23] + qyzi * gqyz[23]))
                        + uzk * (qxxi * gqxx[24] + qyyi * gqyy[24] + qzzi * gqzz[24]
                        + 2.0 * (qxyi * gqxy[24] + qxzi * gqxz[24] + qyzi * gqyz[24]))
                        + qxxi * (qxxk * gqxx[25] + qyyk * gqxx[28] + qzzk * gqxx[30]
                        + 2.0 * (qxyk * gqxx[26] + qxzk * gqxx[27] + qyzk * gqxx[29]))
                        + qyyi * (qxxk * gqyy[25] + qyyk * gqyy[28] + qzzk * gqyy[30]
                        + 2.0 * (qxyk * gqyy[26] + qxzk * gqyy[27] + qyzk * gqyy[29]))
                        + qzzi * (qxxk * gqzz[25] + qyyk * gqzz[28] + qzzk * gqzz[30]
                        + 2.0 * (qxyk * gqzz[26] + qxzk * gqzz[27] + qyzk * gqzz[29]))
                        + 2.0 * (qxyi * (qxxk * gqxy[25] + qyyk * gqxy[28] + qzzk * gqxy[30]
                        + 2.0 * (qxyk * gqxy[26] + qxzk * gqxy[27] + qyzk * gqxy[29]))
                        + qxzi * (qxxk * gqxz[25] + qyyk * gqxz[28] + qzzk * gqxz[30]
                        + 2.0 * (qxyk * gqxz[26] + qxzk * gqxz[27] + qyzk * gqxz[29]))
                        + qyzi * (qxxk * gqyz[25] + qyyk * gqyz[28] + qzzk * gqyz[30]
                        + 2.0 * (qxyk * gqyz[26] + qxzk * gqyz[27] + qyzk * gqyz[29])));
                final double dewkdr = ci * (uxk * gux[21] + uyk * guy[21] + uzk * guz[21])
                        - ck * (uxi * gc[22] + uyi * gc[23] + uzi * gc[24])
                        + ci * (qxxk * gqxx[21] + qyyk * gqyy[21] + qzzk * gqzz[21]
                        + 2.0 * (qxyk * gqxy[21] + qxzk * gqxz[21] + qyzk * gqyz[21]))
                        + ck * (qxxi * gc[25] + qyyi * gc[28] + qzzi * gc[30]
                        + 2.0 * (qxyi * gc[26] + qxzi * gc[27] + qyzi * gc[29]))
                        - uxi * (qxxk * gqxx[22] + qyyk * gqyy[22] + qzzk * gqzz[22]
                        + 2.0 * (qxyk * gqxy[22] + qxzk * gqxz[22] + qyzk * gqyz[22]))
                        - uyi * (qxxk * gqxx[23] + qyyk * gqyy[23] + qzzk * gqzz[23]
                        + 2.0 * (qxyk * gqxy[23] + qxzk * gqxz[23] + qyzk * gqyz[23]))
                        - uzi * (qxxk * gqxx[24] + qyyk * gqyy[24] + qzzk * gqzz[24]
                        + 2.0 * (qxyk * gqxy[24] + qxzk * gqxz[24] + qyzk * gqyz[24]))
                        + uxk * (qxxi * gux[25] + qyyi * gux[28] + qzzi * gux[30]
                        + 2.0 * (qxyi * gux[26] + qxzi * gux[27] + qyzi * gux[29]))
                        + uyk * (qxxi * guy[25] + qyyi * guy[28] + qzzi * guy[30]
                        + 2.0 * (qxyi * guy[26] + qxzi * guy[27] + qyzi * guy[29]))
                        + uzk * (qxxi * guz[25] + qyyi * guz[28] + qzzi * guz[30]
                        + 2.0 * (qxyi * guz[26] + qxzi * guz[27] + qyzi * guz[29]))
                        + qxxi * (qxxk * gqxx[25] + qyyk * gqyy[25] + qzzk * gqzz[25]
                        + 2.0 * (qxyk * gqxy[25] + qxzk * gqxz[25] + qyzk * gqyz[25]))
                        + qyyi * (qxxk * gqxx[28] + qyyk * gqyy[28] + qzzk * gqzz[28]
                        + 2.0 * (qxyk * gqxy[28] + qxzk * gqxz[28] + qyzk * gqyz[28]))
                        + qzzi * (qxxk * gqxx[30] + qyyk * gqyy[30] + qzzk * gqzz[30]
                        + 2.0 * (qxyk * gqxy[30] + qxzk * gqxz[30] + qyzk * gqyz[30]))
                        + 2.0 * (qxyi * (qxxk * gqxx[26] + qyyk * gqyy[26] + qzzk * gqzz[26]
                        + 2.0 * (qxyk * gqxy[26] + qxzk * gqxz[26] + qyzk * gqyz[26]))
                        + qxzi * (qxxk * gqxx[27] + qyyk * gqyy[27] + qzzk * gqzz[27]
                        + 2.0 * (qxyk * gqxy[27] + qxzk * gqxz[27] + qyzk * gqyz[27]))
                        + qyzi * (qxxk * gqxx[29] + qyyk * gqyy[29] + qzzk * gqzz[29]
                        + 2.0 * (qxyk * gqxy[29] + qxzk * gqxz[29] + qyzk * gqyz[29])));
                final double dsumdr = desymdr + 0.5 * (dewidr + dewkdr);
                final double drbi = rbk * dsumdr;
                final double drbk = rbi * dsumdr;

                double selfScale = 1.0;
                if (i == k) {
                    if (iSymm == 0) {
                        gb_local[i] += drbi;
                        return;
                    } else {
                        selfScale = 0.5;
                    }
                }

                // Increment the gradients and Born chain rule term.
                final double dedx = selfScale * dEdX();
                final double dedy = selfScale * dEdY();
                final double dedz = selfScale * dEdZ();

                gX[i] -= lPow * dedx;
                gY[i] -= lPow * dedy;
                gZ[i] -= lPow * dedz;
                gb_local[i] += selfScale * drbi;

                final double dedxk = dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
                final double dedyk = dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
                final double dedzk = dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];

                gX[k] += lPow * dedxk;
                gY[k] += lPow * dedyk;
                gZ[k] += lPow * dedzk;
                gb_local[k] += selfScale * drbk;
                if (lambdaTerm) {
                    lgX[i] -= dlPow * dedx;
                    lgY[i] -= dlPow * dedy;
                    lgZ[i] -= dlPow * dedz;
                    lgX[k] += dlPow * dedxk;
                    lgY[k] += dlPow * dedyk;
                    lgZ[k] += dlPow * dedzk;
                }
                permanentEnergyTorque(i, k);
            }

            private double dEdZ() {
                final double desymdz = ci * ck * gc[4]
                        - (uxi * (uxk * gux[7] + uyk * guy[7] + uzk * guz[7])
                        + uyi * (uxk * gux[9] + uyk * guy[9] + uzk * guz[9])
                        + uzi * (uxk * gux[10] + uyk * guy[10] + uzk * guz[10]));
                final double dewidz = ci * (uxk * gc[7] + uyk * gc[9] + uzk * gc[10])
                        - ck * (uxi * gux[4] + uyi * guy[4] + uzi * guz[4])
                        + ci * (qxxk * gc[13] + qyyk * gc[18] + qzzk * gc[20]
                        + 2.0 * (qxyk * gc[15] + qxzk * gc[16] + qyzk * gc[19]))
                        + ck * (qxxi * gqxx[4] + qyyi * gqyy[4] + qzzi * gqzz[4]
                        + 2.0 * (qxyi * gqxy[4] + qxzi * gqxz[4] + qyzi * gqyz[4]))
                        - uxi * (qxxk * gux[13] + qyyk * gux[18] + qzzk * gux[20]
                        + 2.0 * (qxyk * gux[15] + qxzk * gux[16] + qyzk * gux[19]))
                        - uyi * (qxxk * guy[13] + qyyk * guy[18] + qzzk * guy[20]
                        + 2.0 * (qxyk * guy[15] + qxzk * guy[16] + qyzk * guy[19]))
                        - uzi * (qxxk * guz[13] + qyyk * guz[18] + qzzk * guz[20]
                        + 2.0 * (qxyk * guz[15] + qxzk * guz[16] + qyzk * guz[19]))
                        + uxk * (qxxi * gqxx[7] + qyyi * gqyy[7] + qzzi * gqzz[7]
                        + 2.0 * (qxyi * gqxy[7] + qxzi * gqxz[7] + qyzi * gqyz[7]))
                        + uyk * (qxxi * gqxx[9] + qyyi * gqyy[9] + qzzi * gqzz[9]
                        + 2.0 * (qxyi * gqxy[9] + qxzi * gqxz[9] + qyzi * gqyz[9]))
                        + uzk * (qxxi * gqxx[10] + qyyi * gqyy[10] + qzzi * gqzz[10]
                        + 2.0 * (qxyi * gqxy[10] + qxzi * gqxz[10] + qyzi * gqyz[10]))
                        + qxxi * (qxxk * gqxx[13] + qyyk * gqxx[18] + qzzk * gqxx[20]
                        + 2.0 * (qxyk * gqxx[15] + qxzk * gqxx[16] + qyzk * gqxx[19]))
                        + qyyi * (qxxk * gqyy[13] + qyyk * gqyy[18] + qzzk * gqyy[20]
                        + 2.0 * (qxyk * gqyy[15] + qxzk * gqyy[16] + qyzk * gqyy[19]))
                        + qzzi * (qxxk * gqzz[13] + qyyk * gqzz[18] + qzzk * gqzz[20]
                        + 2.0 * (qxyk * gqzz[15] + qxzk * gqzz[16] + qyzk * gqzz[19]))
                        + 2.0 * (qxyi * (qxxk * gqxy[13] + qyyk * gqxy[18] + qzzk * gqxy[20]
                        + 2.0 * (qxyk * gqxy[15] + qxzk * gqxy[16] + qyzk * gqxy[19]))
                        + qxzi * (qxxk * gqxz[13] + qyyk * gqxz[18] + qzzk * gqxz[20]
                        + 2.0 * (qxyk * gqxz[15] + qxzk * gqxz[16] + qyzk * gqxz[19]))
                        + qyzi * (qxxk * gqyz[13] + qyyk * gqyz[18] + qzzk * gqyz[20]
                        + 2.0 * (qxyk * gqyz[15] + qxzk * gqyz[16] + qyzk * gqyz[19])));
                final double dewkdz = ci * (uxk * gux[4] + uyk * guy[4] + uzk * guz[4])
                        - ck * (uxi * gc[7] + uyi * gc[9] + uzi * gc[10])
                        + ci * (qxxk * gqxx[4] + qyyk * gqyy[4] + qzzk * gqzz[4]
                        + 2.0 * (qxyk * gqxy[4] + qxzk * gqxz[4] + qyzk * gqyz[4]))
                        + ck * (qxxi * gc[13] + qyyi * gc[18] + qzzi * gc[20]
                        + 2.0 * (qxyi * gc[15] + qxzi * gc[16] + qyzi * gc[19]))
                        - uxi * (qxxk * gqxx[7] + qyyk * gqyy[7] + qzzk * gqzz[7]
                        + 2.0 * (qxyk * gqxy[7] + qxzk * gqxz[7] + qyzk * gqyz[7]))
                        - uyi * (qxxk * gqxx[9] + qyyk * gqyy[9] + qzzk * gqzz[9]
                        + 2.0 * (qxyk * gqxy[9] + qxzk * gqxz[9] + qyzk * gqyz[9]))
                        - uzi * (qxxk * gqxx[10] + qyyk * gqyy[10] + qzzk * gqzz[10]
                        + 2.0 * (qxyk * gqxy[10] + qxzk * gqxz[10] + qyzk * gqyz[10]))
                        + uxk * (qxxi * gux[13] + qyyi * gux[18] + qzzi * gux[20]
                        + 2.0 * (qxyi * gux[15] + qxzi * gux[16] + qyzi * gux[19]))
                        + uyk * (qxxi * guy[13] + qyyi * guy[18] + qzzi * guy[20]
                        + 2.0 * (qxyi * guy[15] + qxzi * guy[16] + qyzi * guy[19]))
                        + uzk * (qxxi * guz[13] + qyyi * guz[18] + qzzi * guz[20]
                        + 2.0 * (qxyi * guz[15] + qxzi * guz[16] + qyzi * guz[19]))
                        + qxxi * (qxxk * gqxx[13] + qyyk * gqyy[13] + qzzk * gqzz[13]
                        + 2.0 * (qxyk * gqxy[13] + qxzk * gqxz[13] + qyzk * gqyz[13]))
                        + qyyi * (qxxk * gqxx[18] + qyyk * gqyy[18] + qzzk * gqzz[18]
                        + 2.0 * (qxyk * gqxy[18] + qxzk * gqxz[18] + qyzk * gqyz[18]))
                        + qzzi * (qxxk * gqxx[20] + qyyk * gqyy[20] + qzzk * gqzz[20]
                        + 2.0 * (qxyk * gqxy[20] + qxzk * gqxz[20] + qyzk * gqyz[20]))
                        + 2.0 * (qxyi * (qxxk * gqxx[15] + qyyk * gqyy[15] + qzzk * gqzz[15]
                        + 2.0 * (qxyk * gqxy[15] + qxzk * gqxz[15] + qyzk * gqyz[15]))
                        + qxzi * (qxxk * gqxx[16] + qyyk * gqyy[16] + qzzk * gqzz[16]
                        + 2.0 * (qxyk * gqxy[16] + qxzk * gqxz[16] + qyzk * gqyz[16]))
                        + qyzi * (qxxk * gqxx[19] + qyyk * gqyy[19] + qzzk * gqzz[19]
                        + 2.0 * (qxyk * gqxy[19] + qxzk * gqxz[19] + qyzk * gqyz[19])));
                return desymdz + 0.5 * (dewidz + dewkdz);
            }

            private double dEdY() {
                final double desymdy = ci * ck * gc[3]
                        - (uxi * (uxk * gux[6] + uyk * guy[6] + uzk * guz[6])
                        + uyi * (uxk * gux[8] + uyk * guy[8] + uzk * guz[8])
                        + uzi * (uxk * gux[9] + uyk * guy[9] + uzk * guz[9]));
                final double dewidy = ci * (uxk * gc[6] + uyk * gc[8] + uzk * gc[9])
                        - ck * (uxi * gux[3] + uyi * guy[3] + uzi * guz[3])
                        + ci * (qxxk * gc[12] + qyyk * gc[17] + qzzk * gc[19]
                        + 2.0 * (qxyk * gc[14] + qxzk * gc[15] + qyzk * gc[18]))
                        + ck * (qxxi * gqxx[3] + qyyi * gqyy[3] + qzzi * gqzz[3]
                        + 2.0 * (qxyi * gqxy[3] + qxzi * gqxz[3] + qyzi * gqyz[3]))
                        - uxi * (qxxk * gux[12] + qyyk * gux[17] + qzzk * gux[19]
                        + 2.0 * (qxyk * gux[14] + qxzk * gux[15] + qyzk * gux[18]))
                        - uyi * (qxxk * guy[12] + qyyk * guy[17] + qzzk * guy[19]
                        + 2.0 * (qxyk * guy[14] + qxzk * guy[15] + qyzk * guy[18]))
                        - uzi * (qxxk * guz[12] + qyyk * guz[17] + qzzk * guz[19]
                        + 2.0 * (qxyk * guz[14] + qxzk * guz[15] + qyzk * guz[18]))
                        + uxk * (qxxi * gqxx[6] + qyyi * gqyy[6] + qzzi * gqzz[6]
                        + 2.0 * (qxyi * gqxy[6] + qxzi * gqxz[6] + qyzi * gqyz[6]))
                        + uyk * (qxxi * gqxx[8] + qyyi * gqyy[8] + qzzi * gqzz[8]
                        + 2.0 * (qxyi * gqxy[8] + qxzi * gqxz[8] + qyzi * gqyz[8]))
                        + uzk * (qxxi * gqxx[9] + qyyi * gqyy[9] + qzzi * gqzz[9]
                        + 2.0 * (qxyi * gqxy[9] + qxzi * gqxz[9] + qyzi * gqyz[9]))
                        + qxxi * (qxxk * gqxx[12] + qyyk * gqxx[17] + qzzk * gqxx[19]
                        + 2.0 * (qxyk * gqxx[14] + qxzk * gqxx[15] + qyzk * gqxx[18]))
                        + qyyi * (qxxk * gqyy[12] + qyyk * gqyy[17] + qzzk * gqyy[19]
                        + 2.0 * (qxyk * gqyy[14] + qxzk * gqyy[15] + qyzk * gqyy[18]))
                        + qzzi * (qxxk * gqzz[12] + qyyk * gqzz[17] + qzzk * gqzz[19]
                        + 2.0 * (qxyk * gqzz[14] + qxzk * gqzz[15] + qyzk * gqzz[18]))
                        + 2.0 * (qxyi * (qxxk * gqxy[12] + qyyk * gqxy[17] + qzzk * gqxy[19]
                        + 2.0 * (qxyk * gqxy[14] + qxzk * gqxy[15] + qyzk * gqxy[18]))
                        + qxzi * (qxxk * gqxz[12] + qyyk * gqxz[17] + qzzk * gqxz[19]
                        + 2.0 * (qxyk * gqxz[14] + qxzk * gqxz[15] + qyzk * gqxz[18]))
                        + qyzi * (qxxk * gqyz[12] + qyyk * gqyz[17] + qzzk * gqyz[19]
                        + 2.0 * (qxyk * gqyz[14] + qxzk * gqyz[15] + qyzk * gqyz[18])));
                final double dewkdy = ci * (uxk * gux[3] + uyk * guy[3] + uzk * guz[3])
                        - ck * (uxi * gc[6] + uyi * gc[8] + uzi * gc[9])
                        + ci * (qxxk * gqxx[3] + qyyk * gqyy[3] + qzzk * gqzz[3]
                        + 2.0 * (qxyk * gqxy[3] + qxzk * gqxz[3] + qyzk * gqyz[3]))
                        + ck * (qxxi * gc[12] + qyyi * gc[17] + qzzi * gc[19]
                        + 2.0 * (qxyi * gc[14] + qxzi * gc[15] + qyzi * gc[18]))
                        - uxi * (qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]))
                        - uyi * (qxxk * gqxx[8] + qyyk * gqyy[8] + qzzk * gqzz[8]
                        + 2.0 * (qxyk * gqxy[8] + qxzk * gqxz[8] + qyzk * gqyz[8]))
                        - uzi * (qxxk * gqxx[9] + qyyk * gqyy[9] + qzzk * gqzz[9]
                        + 2.0 * (qxyk * gqxy[9] + qxzk * gqxz[9] + qyzk * gqyz[9]))
                        + uxk * (qxxi * gux[12] + qyyi * gux[17] + qzzi * gux[19]
                        + 2.0 * (qxyi * gux[14] + qxzi * gux[15] + qyzi * gux[18]))
                        + uyk * (qxxi * guy[12] + qyyi * guy[17] + qzzi * guy[19]
                        + 2.0 * (qxyi * guy[14] + qxzi * guy[15] + qyzi * guy[18]))
                        + uzk * (qxxi * guz[12] + qyyi * guz[17] + qzzi * guz[19]
                        + 2.0 * (qxyi * guz[14] + qxzi * guz[15] + qyzi * guz[18]))
                        + qxxi * (qxxk * gqxx[12] + qyyk * gqyy[12] + qzzk * gqzz[12]
                        + 2.0 * (qxyk * gqxy[12] + qxzk * gqxz[12] + qyzk * gqyz[12]))
                        + qyyi * (qxxk * gqxx[17] + qyyk * gqyy[17] + qzzk * gqzz[17]
                        + 2.0 * (qxyk * gqxy[17] + qxzk * gqxz[17] + qyzk * gqyz[17]))
                        + qzzi * (qxxk * gqxx[19] + qyyk * gqyy[19] + qzzk * gqzz[19]
                        + 2.0 * (qxyk * gqxy[19] + qxzk * gqxz[19] + qyzk * gqyz[19]))
                        + 2.0 * (qxyi * (qxxk * gqxx[14] + qyyk * gqyy[14] + qzzk * gqzz[14]
                        + 2.0 * (qxyk * gqxy[14] + qxzk * gqxz[14] + qyzk * gqyz[14]))
                        + qxzi * (qxxk * gqxx[15] + qyyk * gqyy[15] + qzzk * gqzz[15]
                        + 2.0 * (qxyk * gqxy[15] + qxzk * gqxz[15] + qyzk * gqyz[15]))
                        + qyzi * (qxxk * gqxx[18] + qyyk * gqyy[18] + qzzk * gqzz[18]
                        + 2.0 * (qxyk * gqxy[18] + qxzk * gqxz[18] + qyzk * gqyz[18])));

                return desymdy + 0.5 * (dewidy + dewkdy);
            }

            private double dEdX() {
                final double desymdx = ci * ck * gc[2]
                        - (uxi * (uxk * gux[5] + uyk * guy[5] + uzk * guz[5])
                        + uyi * (uxk * gux[6] + uyk * guy[6] + uzk * guz[6])
                        + uzi * (uxk * gux[7] + uyk * guy[7] + uzk * guz[7]));
                final double dewidx = ci * (uxk * gc[5] + uyk * gc[6] + uzk * gc[7])
                        - ck * (uxi * gux[2] + uyi * guy[2] + uzi * guz[2])
                        + ci * (qxxk * gc[11] + qyyk * gc[14] + qzzk * gc[16]
                        + 2.0 * (qxyk * gc[12] + qxzk * gc[13] + qyzk * gc[15]))
                        + ck * (qxxi * gqxx[2] + qyyi * gqyy[2] + qzzi * gqzz[2]
                        + 2.0 * (qxyi * gqxy[2] + qxzi * gqxz[2] + qyzi * gqyz[2]))
                        - uxi * (qxxk * gux[11] + qyyk * gux[14] + qzzk * gux[16]
                        + 2.0 * (qxyk * gux[12] + qxzk * gux[13] + qyzk * gux[15]))
                        - uyi * (qxxk * guy[11] + qyyk * guy[14] + qzzk * guy[16]
                        + 2.0 * (qxyk * guy[12] + qxzk * guy[13] + qyzk * guy[15]))
                        - uzi * (qxxk * guz[11] + qyyk * guz[14] + qzzk * guz[16]
                        + 2.0 * (qxyk * guz[12] + qxzk * guz[13] + qyzk * guz[15]))
                        + uxk * (qxxi * gqxx[5] + qyyi * gqyy[5] + qzzi * gqzz[5]
                        + 2.0 * (qxyi * gqxy[5] + qxzi * gqxz[5] + qyzi * gqyz[5]))
                        + uyk * (qxxi * gqxx[6] + qyyi * gqyy[6] + qzzi * gqzz[6]
                        + 2.0 * (qxyi * gqxy[6] + qxzi * gqxz[6] + qyzi * gqyz[6]))
                        + uzk * (qxxi * gqxx[7] + qyyi * gqyy[7] + qzzi * gqzz[7]
                        + 2.0 * (qxyi * gqxy[7] + qxzi * gqxz[7] + qyzi * gqyz[7]))
                        + qxxi * (qxxk * gqxx[11] + qyyk * gqxx[14] + qzzk * gqxx[16]
                        + 2.0 * (qxyk * gqxx[12] + qxzk * gqxx[13] + qyzk * gqxx[15]))
                        + qyyi * (qxxk * gqyy[11] + qyyk * gqyy[14] + qzzk * gqyy[16]
                        + 2.0 * (qxyk * gqyy[12] + qxzk * gqyy[13] + qyzk * gqyy[15]))
                        + qzzi * (qxxk * gqzz[11] + qyyk * gqzz[14] + qzzk * gqzz[16]
                        + 2.0 * (qxyk * gqzz[12] + qxzk * gqzz[13] + qyzk * gqzz[15]))
                        + 2.0 * (qxyi * (qxxk * gqxy[11] + qyyk * gqxy[14] + qzzk * gqxy[16]
                        + 2.0 * (qxyk * gqxy[12] + qxzk * gqxy[13] + qyzk * gqxy[15]))
                        + qxzi * (qxxk * gqxz[11] + qyyk * gqxz[14] + qzzk * gqxz[16]
                        + 2.0 * (qxyk * gqxz[12] + qxzk * gqxz[13] + qyzk * gqxz[15]))
                        + qyzi * (qxxk * gqyz[11] + qyyk * gqyz[14] + qzzk * gqyz[16]
                        + 2.0 * (qxyk * gqyz[12] + qxzk * gqyz[13] + qyzk * gqyz[15])));
                final double dewkdx = ci * (uxk * gux[2] + uyk * guy[2] + uzk * guz[2])
                        - ck * (uxi * gc[5] + uyi * gc[6] + uzi * gc[7])
                        + ci * (qxxk * gqxx[2] + qyyk * gqyy[2] + qzzk * gqzz[2]
                        + 2.0 * (qxyk * gqxy[2] + qxzk * gqxz[2] + qyzk * gqyz[2]))
                        + ck * (qxxi * gc[11] + qyyi * gc[14] + qzzi * gc[16]
                        + 2.0 * (qxyi * gc[12] + qxzi * gc[13] + qyzi * gc[15]))
                        - uxi * (qxxk * gqxx[5] + qyyk * gqyy[5] + qzzk * gqzz[5]
                        + 2.0 * (qxyk * gqxy[5] + qxzk * gqxz[5] + qyzk * gqyz[5]))
                        - uyi * (qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]))
                        - uzi * (qxxk * gqxx[7] + qyyk * gqyy[7] + qzzk * gqzz[7]
                        + 2.0 * (qxyk * gqxy[7] + qxzk * gqxz[7] + qyzk * gqyz[7]))
                        + uxk * (qxxi * gux[11] + qyyi * gux[14] + qzzi * gux[16]
                        + 2.0 * (qxyi * gux[12] + qxzi * gux[13] + qyzi * gux[15]))
                        + uyk * (qxxi * guy[11] + qyyi * guy[14] + qzzi * guy[16]
                        + 2.0 * (qxyi * guy[12] + qxzi * guy[13] + qyzi * guy[15]))
                        + uzk * (qxxi * guz[11] + qyyi * guz[14] + qzzi * guz[16]
                        + 2.0 * (qxyi * guz[12] + qxzi * guz[13] + qyzi * guz[15]))
                        + qxxi * (qxxk * gqxx[11] + qyyk * gqyy[11] + qzzk * gqzz[11]
                        + 2.0 * (qxyk * gqxy[11] + qxzk * gqxz[11] + qyzk * gqyz[11]))
                        + qyyi * (qxxk * gqxx[14] + qyyk * gqyy[14] + qzzk * gqzz[14]
                        + 2.0 * (qxyk * gqxy[14] + qxzk * gqxz[14] + qyzk * gqyz[14]))
                        + qzzi * (qxxk * gqxx[16] + qyyk * gqyy[16] + qzzk * gqzz[16]
                        + 2.0 * (qxyk * gqxy[16] + qxzk * gqxz[16] + qyzk * gqyz[16]))
                        + 2.0 * (qxyi * (qxxk * gqxx[12] + qyyk * gqyy[12] + qzzk * gqzz[12]
                        + 2.0 * (qxyk * gqxy[12] + qxzk * gqxz[12] + qyzk * gqyz[12]))
                        + qxzi * (qxxk * gqxx[13] + qyyk * gqyy[13] + qzzk * gqzz[13]
                        + 2.0 * (qxyk * gqxy[13] + qxzk * gqxz[13] + qyzk * gqyz[13]))
                        + qyzi * (qxxk * gqxx[15] + qyyk * gqyy[15] + qzzk * gqzz[15]
                        + 2.0 * (qxyk * gqxy[15] + qxzk * gqxz[15] + qyzk * gqyz[15])));

                return desymdx + 0.5 * (dewidx + dewkdx);
            }

            private void permanentEnergyTorque(int i, int k) {

                // Torque on permanent dipoles due to permanent reaction field.
                final double ix = uxk * gux[2] + uyk * gux[3] + uzk * gux[4]
                        + 0.5 * (ck * gux[1] + qxxk * gux[5] + qyyk * gux[8] + qzzk * gux[10]
                        + 2.0 * (qxyk * gux[6] + qxzk * gux[7] + qyzk * gux[9])
                        + ck * gc[2] + qxxk * gqxx[2] + qyyk * gqyy[2] + qzzk * gqzz[2]
                        + 2.0 * (qxyk * gqxy[2] + qxzk * gqxz[2] + qyzk * gqyz[2]));
                final double iy = uxk * guy[2] + uyk * guy[3] + uzk * guy[4]
                        + 0.5 * (ck * guy[1] + qxxk * guy[5] + qyyk * guy[8] + qzzk * guy[10]
                        + 2.0 * (qxyk * guy[6] + qxzk * guy[7] + qyzk * guy[9])
                        + ck * gc[3] + qxxk * gqxx[3] + qyyk * gqyy[3] + qzzk * gqzz[3]
                        + 2.0 * (qxyk * gqxy[3] + qxzk * gqxz[3] + qyzk * gqyz[3]));
                final double iz = uxk * guz[2] + uyk * guz[3] + uzk * guz[4]
                        + 0.5 * (ck * guz[1] + qxxk * guz[5] + qyyk * guz[8] + qzzk * guz[10]
                        + 2.0 * (qxyk * guz[6] + qxzk * guz[7] + qyzk * guz[9])
                        + ck * gc[4] + qxxk * gqxx[4] + qyyk * gqyy[4] + qzzk * gqzz[4]
                        + 2.0 * (qxyk * gqxy[4] + qxzk * gqxz[4] + qyzk * gqyz[4]));
                final double kx = uxi * gux[2] + uyi * gux[3] + uzi * gux[4]
                        - 0.5 * (ci * gux[1] + qxxi * gux[5] + qyyi * gux[8] + qzzi * gux[10]
                        + 2.0 * (qxyi * gux[6] + qxzi * gux[7] + qyzi * gux[9])
                        + ci * gc[2] + qxxi * gqxx[2] + qyyi * gqyy[2] + qzzi * gqzz[2]
                        + 2.0 * (qxyi * gqxy[2] + qxzi * gqxz[2] + qyzi * gqyz[2]));
                final double ky = uxi * guy[2] + uyi * guy[3] + uzi * guy[4]
                        - 0.5 * (ci * guy[1] + qxxi * guy[5] + qyyi * guy[8] + qzzi * guy[10]
                        + 2.0 * (qxyi * guy[6] + qxzi * guy[7] + qyzi * guy[9])
                        + ci * gc[3] + qxxi * gqxx[3] + qyyi * gqyy[3] + qzzi * gqzz[3]
                        + 2.0 * (qxyi * gqxy[3] + qxzi * gqxz[3] + qyzi * gqyz[3]));
                final double kz = uxi * guz[2] + uyi * guz[3] + uzi * guz[4]
                        - 0.5 * (ci * guz[1] + qxxi * guz[5] + qyyi * guz[8] + qzzi * guz[10]
                        + 2.0 * (qxyi * guz[6] + qxzi * guz[7] + qyzi * guz[9])
                        + ci * gc[4] + qxxi * gqxx[4] + qyyi * gqyy[4] + qzzi * gqzz[4]
                        + 2.0 * (qxyi * gqxy[4] + qxzi * gqxz[4] + qyzi * gqyz[4]));
                double tix = uyi * iz - uzi * iy;
                double tiy = uzi * ix - uxi * iz;
                double tiz = uxi * iy - uyi * ix;
                double tkx = uyk * kz - uzk * ky;
                double tky = uzk * kx - uxk * kz;
                double tkz = uxk * ky - uyk * kx;

                // Torque on quadrupoles due to permanent reaction field gradient.
                final double ixx
                        = -0.5 * (ck * gqxx[1] + uxk * gqxx[2] + uyk * gqxx[3] + uzk * gqxx[4]
                        + qxxk * gqxx[5] + qyyk * gqxx[8] + qzzk * gqxx[10]
                        + 2.0 * (qxyk * gqxx[6] + qxzk * gqxx[7] + qyzk * gqxx[9])
                        + ck * gc[5] + uxk * gux[5] + uyk * guy[5] + uzk * guz[5]
                        + qxxk * gqxx[5] + qyyk * gqyy[5] + qzzk * gqzz[5]
                        + 2.0 * (qxyk * gqxy[5] + qxzk * gqxz[5] + qyzk * gqyz[5]));
                final double ixy
                        = -0.5 * (ck * gqxy[1] + uxk * gqxy[2] + uyk * gqxy[3] + uzk * gqxy[4]
                        + qxxk * gqxy[5] + qyyk * gqxy[8] + qzzk * gqxy[10]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxy[7] + qyzk * gqxy[9])
                        + ck * gc[6] + uxk * gux[6] + uyk * guy[6] + uzk * guz[6]
                        + qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]));
                final double ixz
                        = -0.5 * (ck * gqxz[1] + uxk * gqxz[2] + uyk * gqxz[3] + uzk * gqxz[4]
                        + qxxk * gqxz[5] + qyyk * gqxz[8] + qzzk * gqxz[10]
                        + 2.0 * (qxyk * gqxz[6] + qxzk * gqxz[7] + qyzk * gqxz[9])
                        + ck * gc[7] + uxk * gux[7] + uyk * guy[7] + uzk * guz[7]
                        + qxxk * gqxx[7] + qyyk * gqyy[7] + qzzk * gqzz[7]
                        + 2.0 * (qxyk * gqxy[7] + qxzk * gqxz[7] + qyzk * gqyz[7]));
                final double iyy
                        = -0.5 * (ck * gqyy[1] + uxk * gqyy[2] + uyk * gqyy[3] + uzk * gqyy[4]
                        + qxxk * gqyy[5] + qyyk * gqyy[8] + qzzk * gqyy[10]
                        + 2.0 * (qxyk * gqyy[6] + qxzk * gqyy[7] + qyzk * gqyy[9])
                        + ck * gc[8] + uxk * gux[8] + uyk * guy[8] + uzk * guz[8]
                        + qxxk * gqxx[8] + qyyk * gqyy[8] + qzzk * gqzz[8]
                        + 2.0 * (qxyk * gqxy[8] + qxzk * gqxz[8] + qyzk * gqyz[8]));
                final double iyz
                        = -0.5 * (ck * gqyz[1] + uxk * gqyz[2] + uyk * gqyz[3] + uzk * gqyz[4]
                        + qxxk * gqyz[5] + qyyk * gqyz[8] + qzzk * gqyz[10]
                        + 2.0 * (qxyk * gqyz[6] + qxzk * gqyz[7] + qyzk * gqyz[9])
                        + ck * gc[9] + uxk * gux[9] + uyk * guy[9] + uzk * guz[9]
                        + qxxk * gqxx[9] + qyyk * gqyy[9] + qzzk * gqzz[9]
                        + 2.0 * (qxyk * gqxy[9] + qxzk * gqxz[9] + qyzk * gqyz[9]));
                final double izz
                        = -0.5 * (ck * gqzz[1] + uxk * gqzz[2] + uyk * gqzz[3] + uzk * gqzz[4]
                        + qxxk * gqzz[5] + qyyk * gqzz[8] + qzzk * gqzz[10]
                        + 2.0 * (qxyk * gqzz[6] + qxzk * gqzz[7] + qyzk * gqzz[9])
                        + ck * gc[10] + uxk * gux[10] + uyk * guy[10] + uzk * guz[10]
                        + qxxk * gqxx[10] + qyyk * gqyy[10] + qzzk * gqzz[10]
                        + 2.0 * (qxyk * gqxy[10] + qxzk * gqxz[10] + qyzk * gqyz[10]));
                final double iyx = ixy;
                final double izx = ixz;
                final double izy = iyz;
                final double kxx
                        = -0.5 * (ci * gqxx[1] - uxi * gqxx[2] - uyi * gqxx[3] - uzi * gqxx[4]
                        + qxxi * gqxx[5] + qyyi * gqxx[8] + qzzi * gqxx[10]
                        + 2.0 * (qxyi * gqxx[6] + qxzi * gqxx[7] + qyzi * gqxx[9])
                        + ci * gc[5] - uxi * gux[5] - uyi * guy[5] - uzi * guz[5]
                        + qxxi * gqxx[5] + qyyi * gqyy[5] + qzzi * gqzz[5]
                        + 2.0 * (qxyi * gqxy[5] + qxzi * gqxz[5] + qyzi * gqyz[5]));
                double kxy
                        = -0.5 * (ci * gqxy[1] - uxi * gqxy[2] - uyi * gqxy[3] - uzi * gqxy[4]
                        + qxxi * gqxy[5] + qyyi * gqxy[8] + qzzi * gqxy[10]
                        + 2.0 * (qxyi * gqxy[6] + qxzi * gqxy[7] + qyzi * gqxy[9])
                        + ci * gc[6] - uxi * gux[6] - uyi * guy[6] - uzi * guz[6]
                        + qxxi * gqxx[6] + qyyi * gqyy[6] + qzzi * gqzz[6]
                        + 2.0 * (qxyi * gqxy[6] + qxzi * gqxz[6] + qyzi * gqyz[6]));
                double kxz
                        = -0.5 * (ci * gqxz[1] - uxi * gqxz[2] - uyi * gqxz[3] - uzi * gqxz[4]
                        + qxxi * gqxz[5] + qyyi * gqxz[8] + qzzi * gqxz[10]
                        + 2.0 * (qxyi * gqxz[6] + qxzi * gqxz[7] + qyzi * gqxz[9])
                        + ci * gc[7] - uxi * gux[7] - uyi * guy[7] - uzi * guz[7]
                        + qxxi * gqxx[7] + qyyi * gqyy[7] + qzzi * gqzz[7]
                        + 2.0 * (qxyi * gqxy[7] + qxzi * gqxz[7] + qyzi * gqyz[7]));
                double kyy
                        = -0.5 * (ci * gqyy[1] - uxi * gqyy[2] - uyi * gqyy[3] - uzi * gqyy[4]
                        + qxxi * gqyy[5] + qyyi * gqyy[8] + qzzi * gqyy[10]
                        + 2.0 * (qxyi * gqyy[6] + qxzi * gqyy[7] + qyzi * gqyy[9])
                        + ci * gc[8] - uxi * gux[8] - uyi * guy[8] - uzi * guz[8]
                        + qxxi * gqxx[8] + qyyi * gqyy[8] + qzzi * gqzz[8]
                        + 2.0 * (qxyi * gqxy[8] + qxzi * gqxz[8] + qyzi * gqyz[8]));
                double kyz
                        = -0.5 * (ci * gqyz[1] - uxi * gqyz[2] - uyi * gqyz[3] - uzi * gqyz[4]
                        + qxxi * gqyz[5] + qyyi * gqyz[8] + qzzi * gqyz[10]
                        + 2.0 * (qxyi * gqyz[6] + qxzi * gqyz[7] + qyzi * gqyz[9])
                        + ci * gc[9] - uxi * gux[9] - uyi * guy[9] - uzi * guz[9]
                        + qxxi * gqxx[9] + qyyi * gqyy[9] + qzzi * gqzz[9]
                        + 2.0 * (qxyi * gqxy[9] + qxzi * gqxz[9] + qyzi * gqyz[9]));
                double kzz
                        = -0.5 * (ci * gqzz[1] - uxi * gqzz[2] - uyi * gqzz[3] - uzi * gqzz[4]
                        + qxxi * gqzz[5] + qyyi * gqzz[8] + qzzi * gqzz[10]
                        + 2.0 * (qxyi * gqzz[6] + qxzi * gqzz[7] + qyzi * gqzz[9])
                        + ci * gc[10] - uxi * gux[10] - uyi * guy[10] - uzi * guz[10]
                        + qxxi * gqxx[10] + qyyi * gqyy[10] + qzzi * gqzz[10]
                        + 2.0 * (qxyi * gqxy[10] + qxzi * gqxz[10] + qyzi * gqyz[10]));
                double kyx = kxy;
                double kzx = kxz;
                double kzy = kyz;
                tix += 2.0 * (qxyi * ixz + qyyi * iyz + qyzi * izz - qxzi * ixy - qyzi * iyy - qzzi * izy);
                tiy += 2.0 * (qxzi * ixx + qyzi * iyx + qzzi * izx - qxxi * ixz - qxyi * iyz - qxzi * izz);
                tiz += 2.0 * (qxxi * ixy + qxyi * iyy + qxzi * izy - qxyi * ixx - qyyi * iyx - qyzi * izx);
                tkx += 2.0 * (qxyk * kxz + qyyk * kyz + qyzk * kzz - qxzk * kxy - qyzk * kyy - qzzk * kzy);
                tky += 2.0 * (qxzk * kxx + qyzk * kyx + qzzk * kzx - qxxk * kxz - qxyk * kyz - qxzk * kzz);
                tkz += 2.0 * (qxxk * kxy + qxyk * kyy + qxzk * kzy - qxyk * kxx - qyyk * kyx - qyzk * kzx);

                if (i == k) {
                    double selfScale = 0.5;
                    tix *= selfScale;
                    tiy *= selfScale;
                    tiz *= selfScale;
                    tkx *= selfScale;
                    tky *= selfScale;
                    tkz *= selfScale;
                }

                tX[i] += lPow * tix;
                tY[i] += lPow * tiy;
                tZ[i] += lPow * tiz;

                final double rtkx = tkx * transOp[0][0] + tky * transOp[1][0] + tkz * transOp[2][0];
                final double rtky = tkx * transOp[0][1] + tky * transOp[1][1] + tkz * transOp[2][1];
                final double rtkz = tkx * transOp[0][2] + tky * transOp[1][2] + tkz * transOp[2][2];
                tX[k] += lPow * rtkx;
                tY[k] += lPow * rtky;
                tZ[k] += lPow * rtkz;
                if (lambdaTerm) {
                    ltX[i] += dlPow * tix;
                    ltY[i] += dlPow * tiy;
                    ltZ[i] += dlPow * tiz;
                    ltX[k] += dlPow * rtkx;
                    ltY[k] += dlPow * rtky;
                    ltZ[k] += dlPow * rtkz;
                }
            }

            private void polarizationEnergyGradient(int i, int k) {
                // Electrostatic solvation free energy gradient of the permanent
                // multipoles in the reaction potential of the induced dipoles.
                final double dpsymdx = -uxi * (sxk * gux[5] + syk * guy[5] + szk * guz[5])
                        - uyi * (sxk * gux[6] + syk * guy[6] + szk * guz[6])
                        - uzi * (sxk * gux[7] + syk * guy[7] + szk * guz[7])
                        - uxk * (sxi * gux[5] + syi * guy[5] + szi * guz[5])
                        - uyk * (sxi * gux[6] + syi * guy[6] + szi * guz[6])
                        - uzk * (sxi * gux[7] + syi * guy[7] + szi * guz[7]);
                final double dpwidx = ci * (sxk * gc[5] + syk * gc[6] + szk * gc[7])
                        - ck * (sxi * gux[2] + syi * guy[2] + szi * guz[2])
                        - sxi * (qxxk * gux[11] + qyyk * gux[14] + qzzk * gux[16]
                        + 2.0 * (qxyk * gux[12] + qxzk * gux[13] + qyzk * gux[15]))
                        - syi * (qxxk * guy[11] + qyyk * guy[14] + qzzk * guy[16]
                        + 2.0 * (qxyk * guy[12] + qxzk * guy[13] + qyzk * guy[15]))
                        - szi * (qxxk * guz[11] + qyyk * guz[14] + qzzk * guz[16]
                        + 2.0 * (qxyk * guz[12] + qxzk * guz[13] + qyzk * guz[15]))
                        + sxk * (qxxi * gqxx[5] + qyyi * gqyy[5] + qzzi * gqzz[5]
                        + 2.0 * (qxyi * gqxy[5] + qxzi * gqxz[5] + qyzi * gqyz[5]))
                        + syk * (qxxi * gqxx[6] + qyyi * gqyy[6] + qzzi * gqzz[6]
                        + 2.0 * (qxyi * gqxy[6] + qxzi * gqxz[6] + qyzi * gqyz[6]))
                        + szk * (qxxi * gqxx[7] + qyyi * gqyy[7] + qzzi * gqzz[7]
                        + 2.0 * (qxyi * gqxy[7] + qxzi * gqxz[7] + qyzi * gqyz[7]));
                final double dpwkdx = ci * (sxk * gux[2] + syk * guy[2] + szk * guz[2])
                        - ck * (sxi * gc[5] + syi * gc[6] + szi * gc[7])
                        - sxi * (qxxk * gqxx[5] + qyyk * gqyy[5] + qzzk * gqzz[5]
                        + 2.0 * (qxyk * gqxy[5] + qxzk * gqxz[5] + qyzk * gqyz[5]))
                        - syi * (qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]))
                        - szi * (qxxk * gqxx[7] + qyyk * gqyy[7] + qzzk * gqzz[7]
                        + 2.0 * (qxyk * gqxy[7] + qxzk * gqxz[7] + qyzk * gqyz[7]))
                        + sxk * (qxxi * gux[11] + qyyi * gux[14] + qzzi * gux[16]
                        + 2.0 * (qxyi * gux[12] + qxzi * gux[13] + qyzi * gux[15]))
                        + syk * (qxxi * guy[11] + qyyi * guy[14] + qzzi * guy[16]
                        + 2.0 * (qxyi * guy[12] + qxzi * guy[13] + qyzi * guy[15]))
                        + szk * (qxxi * guz[11] + qyyi * guz[14] + qzzi * guz[16]
                        + 2.0 * (qxyi * guz[12] + qxzi * guz[13] + qyzi * guz[15]));
                final double dpsymdy = -uxi * (sxk * gux[6] + syk * guy[6] + szk * guz[6])
                        - uyi * (sxk * gux[8] + syk * guy[8] + szk * guz[8])
                        - uzi * (sxk * gux[9] + syk * guy[9] + szk * guz[9])
                        - uxk * (sxi * gux[6] + syi * guy[6] + szi * guz[6])
                        - uyk * (sxi * gux[8] + syi * guy[8] + szi * guz[8])
                        - uzk * (sxi * gux[9] + syi * guy[9] + szi * guz[9]);
                final double dpwidy = ci * (sxk * gc[6] + syk * gc[8] + szk * gc[9])
                        - ck * (sxi * gux[3] + syi * guy[3] + szi * guz[3])
                        - sxi * (qxxk * gux[12] + qyyk * gux[17] + qzzk * gux[19]
                        + 2.0 * (qxyk * gux[14] + qxzk * gux[15] + qyzk * gux[18]))
                        - syi * (qxxk * guy[12] + qyyk * guy[17] + qzzk * guy[19]
                        + 2.0 * (qxyk * guy[14] + qxzk * guy[15] + qyzk * guy[18]))
                        - szi * (qxxk * guz[12] + qyyk * guz[17] + qzzk * guz[19]
                        + 2.0 * (qxyk * guz[14] + qxzk * guz[15] + qyzk * guz[18]))
                        + sxk * (qxxi * gqxx[6] + qyyi * gqyy[6] + qzzi * gqzz[6]
                        + 2.0 * (qxyi * gqxy[6] + qxzi * gqxz[6] + qyzi * gqyz[6]))
                        + syk * (qxxi * gqxx[8] + qyyi * gqyy[8] + qzzi * gqzz[8]
                        + 2.0 * (qxyi * gqxy[8] + qxzi * gqxz[8] + qyzi * gqyz[8]))
                        + szk * (qxxi * gqxx[9] + qyyi * gqyy[9] + qzzi * gqzz[9]
                        + 2.0 * (qxyi * gqxy[9] + qxzi * gqxz[9] + qyzi * gqyz[9]));
                final double dpwkdy = ci * (sxk * gux[3] + syk * guy[3] + szk * guz[3])
                        - ck * (sxi * gc[6] + syi * gc[8] + szi * gc[9])
                        - sxi * (qxxk * gqxx[6] + qyyk * gqyy[6] + qzzk * gqzz[6]
                        + 2.0 * (qxyk * gqxy[6] + qxzk * gqxz[6] + qyzk * gqyz[6]))
                        - syi * (qxxk * gqxx[8] + qyyk * gqyy[8] + qzzk * gqzz[8]
                        + 2.0 * (qxyk * gqxy[8] + qxzk * gqxz[8] + qyzk * gqyz[8]))
                        - szi * (qxxk * gqxx[9] + qyyk * gqyy[9] + qzzk * gqzz[9]
                        + 2.0 * (qxyk * gqxy[9] + qxzk * gqxz[9] + qyzk * gqyz[9]))
                        + sxk * (qxxi * gux[12] + qyyi * gux[17] + qzzi * gux[19]
                        + 2.0 * (qxyi * gux[14] + qxzi * gux[15] + qyzi * gux[18]))
                        + syk * (qxxi * guy[12] + qyyi * guy[17] + qzzi * guy[19]
                        + 2.0 * (qxyi * guy[14] + qxzi * guy[15] + qyzi * guy[18]))
                        + szk * (qxxi * guz[12] + qyyi * guz[17] + qzzi * guz[19]
                        + 2.0 * (qxyi * guz[14] + qxzi * guz[15] + qyzi * guz[18]));
                final double dpsymdz = -uxi * (sxk * gux[7] + syk * guy[7] + szk * guz[7])
                        - uyi * (sxk * gux[9] + syk * guy[9] + szk * guz[9])
                        - uzi * (sxk * gux[10] + syk * guy[10] + szk * guz[10])
                        - uxk * (sxi * gux[7] + syi * guy[7] + szi * guz[7])
                        - uyk * (sxi * gux[9] + syi * guy[9] + szi * guz[9])
                        - uzk * (sxi * gux[10] + syi * guy[10] + szi * guz[10]);
                final double dpwidz = ci * (sxk * gc[7] + syk * gc[9] + szk * gc[10])
                        - ck * (sxi * gux[4] + syi * guy[4] + szi * guz[4])
                        - sxi * (qxxk * gux[13] + qyyk * gux[18] + qzzk * gux[20]
                        + 2.0 * (qxyk * gux[15] + qxzk * gux[16] + qyzk * gux[19]))
                        - syi * (qxxk * guy[13] + qyyk * guy[18] + qzzk * guy[20]
                        + 2.0 * (qxyk * guy[15] + qxzk * guy[16] + qyzk * guy[19]))
                        - szi * (qxxk * guz[13] + qyyk * guz[18] + qzzk * guz[20]
                        + 2.0 * (qxyk * guz[15] + qxzk * guz[16] + qyzk * guz[19]))
                        + sxk * (qxxi * gqxx[7] + qyyi * gqyy[7] + qzzi * gqzz[7]
                        + 2.0 * (qxyi * gqxy[7] + qxzi * gqxz[7] + qyzi * gqyz[7]))
                        + syk * (qxxi * gqxx[9] + qyyi * gqyy[9] + qzzi * gqzz[9]
                        + 2.0 * (qxyi * gqxy[9] + qxzi * gqxz[9] + qyzi * gqyz[9]))
                        + szk * (qxxi * gqxx[10] + qyyi * gqyy[10] + qzzi * gqzz[10]
                        + 2.0 * (qxyi * gqxy[10] + qxzi * gqxz[10] + qyzi * gqyz[10]));
                final double dpwkdz = ci * (sxk * gux[4] + syk * guy[4] + szk * guz[4])
                        - ck * (sxi * gc[7] + syi * gc[9] + szi * gc[10])
                        - sxi * (qxxk * gqxx[7] + qyyk * gqyy[7] + qzzk * gqzz[7]
                        + 2.0 * (qxyk * gqxy[7] + qxzk * gqxz[7] + qyzk * gqyz[7]))
                        - syi * (qxxk * gqxx[9] + qyyk * gqyy[9] + qzzk * gqzz[9]
                        + 2.0 * (qxyk * gqxy[9] + qxzk * gqxz[9] + qyzk * gqyz[9]))
                        - szi * (qxxk * gqxx[10] + qyyk * gqyy[10] + qzzk * gqzz[10]
                        + 2.0 * (qxyk * gqxy[10] + qxzk * gqxz[10] + qyzk * gqyz[10]))
                        + sxk * (qxxi * gux[13] + qyyi * gux[18] + qzzi * gux[20]
                        + 2.0 * (qxyi * gux[15] + qxzi * gux[16] + qyzi * gux[19]))
                        + syk * (qxxi * guy[13] + qyyi * guy[18] + qzzi * guy[20]
                        + 2.0 * (qxyi * guy[15] + qxzi * guy[16] + qyzi * guy[19]))
                        + szk * (qxxi * guz[13] + qyyi * guz[18] + qzzi * guz[20]
                        + 2.0 * (qxyi * guz[15] + qxzi * guz[16] + qyzi * guz[19]));

                // Effective radii chain rule terms for the electrostatic solvation free energy
                // gradient of the permanent multipoles in the reaction potential of the induced dipoles.
                final double dsymdr = -uxi * (sxk * gux[22] + syk * guy[22] + szk * guz[22])
                        - uyi * (sxk * gux[23] + syk * guy[23] + szk * guz[23])
                        - uzi * (sxk * gux[24] + syk * guy[24] + szk * guz[24])
                        - uxk * (sxi * gux[22] + syi * guy[22] + szi * guz[22])
                        - uyk * (sxi * gux[23] + syi * guy[23] + szi * guz[23])
                        - uzk * (sxi * gux[24] + syi * guy[24] + szi * guz[24]);
                final double dwipdr = ci * (sxk * gc[22] + syk * gc[23] + szk * gc[24])
                        - ck * (sxi * gux[21] + syi * guy[21] + szi * guz[21])
                        - sxi * (qxxk * gux[25] + qyyk * gux[28] + qzzk * gux[30]
                        + 2.0 * (qxyk * gux[26] + qxzk * gux[27] + qyzk * gux[29]))
                        - syi * (qxxk * guy[25] + qyyk * guy[28] + qzzk * guy[30]
                        + 2.0 * (qxyk * guy[26] + qxzk * guy[27] + qyzk * guy[29]))
                        - szi * (qxxk * guz[25] + qyyk * guz[28] + qzzk * guz[30]
                        + 2.0 * (qxyk * guz[26] + qxzk * guz[27] + qyzk * guz[29]))
                        + sxk * (qxxi * gqxx[22] + qyyi * gqyy[22] + qzzi * gqzz[22]
                        + 2.0 * (qxyi * gqxy[22] + qxzi * gqxz[22] + qyzi * gqyz[22]))
                        + syk * (qxxi * gqxx[23] + qyyi * gqyy[23] + qzzi * gqzz[23]
                        + 2.0 * (qxyi * gqxy[23] + qxzi * gqxz[23] + qyzi * gqyz[23]))
                        + szk * (qxxi * gqxx[24] + qyyi * gqyy[24] + qzzi * gqzz[24]
                        + 2.0 * (qxyi * gqxy[24] + qxzi * gqxz[24] + qyzi * gqyz[24]));
                final double dwkpdr = ci * (sxk * gux[21] + syk * guy[21] + szk * guz[21])
                        - ck * (sxi * gc[22] + syi * gc[23] + szi * gc[24])
                        - sxi * (qxxk * gqxx[22] + qyyk * gqyy[22] + qzzk * gqzz[22]
                        + 2.0 * (qxyk * gqxy[22] + qxzk * gqxz[22] + qyzk * gqyz[22]))
                        - syi * (qxxk * gqxx[23] + qyyk * gqyy[23] + qzzk * gqzz[23]
                        + 2.0 * (qxyk * gqxy[23] + qxzk * gqxz[23] + qyzk * gqyz[23]))
                        - szi * (qxxk * gqxx[24] + qyyk * gqyy[24] + qzzk * gqzz[24]
                        + 2.0 * (qxyk * gqxy[24] + qxzk * gqxz[24] + qyzk * gqyz[24]))
                        + sxk * (qxxi * gux[25] + qyyi * gux[28] + qzzi * gux[30]
                        + 2.0 * (qxyi * gux[26] + qxzi * gux[27] + qyzi * gux[29]))
                        + syk * (qxxi * guy[25] + qyyi * guy[28] + qzzi * guy[30]
                        + 2.0 * (qxyi * guy[26] + qxzi * guy[27] + qyzi * guy[29]))
                        + szk * (qxxi * guz[25] + qyyi * guz[28] + qzzi * guz[30]
                        + 2.0 * (qxyi * guz[26] + qxzi * guz[27] + qyzi * guz[29]));
                double dpdx = 0.5 * (dpsymdx + 0.5 * (dpwidx + dpwkdx));
                double dpdy = 0.5 * (dpsymdy + 0.5 * (dpwidy + dpwkdy));
                double dpdz = 0.5 * (dpsymdz + 0.5 * (dpwidz + dpwkdz));
                double dsumdri = dsymdr + 0.5 * (dwipdr + dwkpdr);
                double dbi = 0.5 * rbk * dsumdri;
                double dbk = 0.5 * rbi * dsumdri;
                if (polarization == Polarization.MUTUAL) {
                    dpdx -= 0.5 * (dxi * (pxk * gux[5] + pyk * gux[6] + pzk * gux[7])
                            + dyi * (pxk * guy[5] + pyk * guy[6] + pzk * guy[7])
                            + dzi * (pxk * guz[5] + pyk * guz[6] + pzk * guz[7])
                            + dxk * (pxi * gux[5] + pyi * gux[6] + pzi * gux[7])
                            + dyk * (pxi * guy[5] + pyi * guy[6] + pzi * guy[7])
                            + dzk * (pxi * guz[5] + pyi * guz[6] + pzi * guz[7]));
                    dpdy -= 0.5 * (dxi * (pxk * gux[6] + pyk * gux[8] + pzk * gux[9])
                            + dyi * (pxk * guy[6] + pyk * guy[8] + pzk * guy[9])
                            + dzi * (pxk * guz[6] + pyk * guz[8] + pzk * guz[9])
                            + dxk * (pxi * gux[6] + pyi * gux[8] + pzi * gux[9])
                            + dyk * (pxi * guy[6] + pyi * guy[8] + pzi * guy[9])
                            + dzk * (pxi * guz[6] + pyi * guz[8] + pzi * guz[9]));
                    dpdz -= 0.5 * (dxi * (pxk * gux[7] + pyk * gux[9] + pzk * gux[10])
                            + dyi * (pxk * guy[7] + pyk * guy[9] + pzk * guy[10])
                            + dzi * (pxk * guz[7] + pyk * guz[9] + pzk * guz[10])
                            + dxk * (pxi * gux[7] + pyi * gux[9] + pzi * gux[10])
                            + dyk * (pxi * guy[7] + pyi * guy[9] + pzi * guy[10])
                            + dzk * (pxi * guz[7] + pyi * guz[9] + pzi * guz[10]));
                    final double duvdr = dxi * (pxk * gux[22] + pyk * gux[23] + pzk * gux[24])
                            + dyi * (pxk * guy[22] + pyk * guy[23] + pzk * guy[24])
                            + dzi * (pxk * guz[22] + pyk * guz[23] + pzk * guz[24])
                            + dxk * (pxi * gux[22] + pyi * gux[23] + pzi * gux[24])
                            + dyk * (pxi * guy[22] + pyi * guy[23] + pzi * guy[24])
                            + dzk * (pxi * guz[22] + pyi * guz[23] + pzi * guz[24]);
                    dbi -= 0.5 * rbk * duvdr;
                    dbk -= 0.5 * rbi * duvdr;
                }

                // Increment the gradients and Born chain rule term.
                if (i == k && iSymm == 0) {
                    gb_local[i] += dbi;
                } else {
                    if (i == k) {
                        dpdx *= 0.5;
                        dpdy *= 0.5;
                        dpdz *= 0.5;
                        dbi *= 0.5;
                        dbk *= 0.5;
                    }
                    gX[i] -= lPow * dpdx;
                    gY[i] -= lPow * dpdy;
                    gZ[i] -= lPow * dpdz;
                    gb_local[i] += dbi;

                    final double rdpdx = dpdx * transOp[0][0] + dpdy * transOp[1][0] + dpdz * transOp[2][0];
                    final double rdpdy = dpdx * transOp[0][1] + dpdy * transOp[1][1] + dpdz * transOp[2][1];
                    final double rdpdz = dpdx * transOp[0][2] + dpdy * transOp[1][2] + dpdz * transOp[2][2];
                    gX[k] += lPow * rdpdx;
                    gY[k] += lPow * rdpdy;
                    gZ[k] += lPow * rdpdz;
                    gb_local[k] += dbk;

                    if (lambdaTerm) {
                        lgX[i] -= dlPow * dpdx;
                        lgY[i] -= dlPow * dpdy;
                        lgZ[i] -= dlPow * dpdz;
                        lgX[k] += dlPow * rdpdx;
                        lgY[k] += dlPow * rdpdy;
                        lgZ[k] += dlPow * rdpdz;
                    }
                }
                polarizationEnergyTorque(i, k);
            }

            private void polarizationEnergyTorque(int i, int k) {
                double fix = 0.5 * (sxk * gux[2] + syk * guy[2] + szk * guz[2]);
                double fiy = 0.5 * (sxk * gux[3] + syk * guy[3] + szk * guz[3]);
                double fiz = 0.5 * (sxk * gux[4] + syk * guy[4] + szk * guz[4]);
                double fkx = 0.5 * (sxi * gux[2] + syi * guy[2] + szi * guz[2]);
                double fky = 0.5 * (sxi * gux[3] + syi * guy[3] + szi * guz[3]);
                double fkz = 0.5 * (sxi * gux[4] + syi * guy[4] + szi * guz[4]);
                double fixx = -0.25
                        * ((sxk * gqxx[2] + syk * gqxx[3] + szk * gqxx[4])
                        + (sxk * gux[5] + syk * guy[5] + szk * guz[5]));
                double fixy = -0.25
                        * ((sxk * gqxy[2] + syk * gqxy[3] + szk * gqxy[4])
                        + (sxk * gux[6] + syk * guy[6] + szk * guz[6]));
                double fixz = -0.25
                        * ((sxk * gqxz[2] + syk * gqxz[3] + szk * gqxz[4])
                        + (sxk * gux[7] + syk * guy[7] + szk * guz[7]));
                double fiyy = -0.25
                        * ((sxk * gqyy[2] + syk * gqyy[3] + szk * gqyy[4])
                        + (sxk * gux[8] + syk * guy[8] + szk * guz[8]));
                double fiyz = -0.25
                        * ((sxk * gqyz[2] + syk * gqyz[3] + szk * gqyz[4])
                        + (sxk * gux[9] + syk * guy[9] + szk * guz[9]));
                double fizz = -0.25
                        * ((sxk * gqzz[2] + syk * gqzz[3] + szk * gqzz[4])
                        + (sxk * gux[10] + syk * guy[10] + szk * guz[10]));
                double fiyx = fixy;
                double fizx = fixz;
                double fizy = fiyz;
                double fkxx = 0.25
                        * ((sxi * gqxx[2] + syi * gqxx[3] + szi * gqxx[4])
                        + (sxi * gux[5] + syi * guy[5] + szi * guz[5]));
                double fkxy = 0.25
                        * ((sxi * gqxy[2] + syi * gqxy[3] + szi * gqxy[4])
                        + (sxi * gux[6] + syi * guy[6] + szi * guz[6]));
                double fkxz = 0.25
                        * ((sxi * gqxz[2] + syi * gqxz[3] + szi * gqxz[4])
                        + (sxi * gux[7] + syi * guy[7] + szi * guz[7]));
                double fkyy = 0.25
                        * ((sxi * gqyy[2] + syi * gqyy[3] + szi * gqyy[4])
                        + (sxi * gux[8] + syi * guy[8] + szi * guz[8]));
                double fkyz = 0.25
                        * ((sxi * gqyz[2] + syi * gqyz[3] + szi * gqyz[4])
                        + (sxi * gux[9] + syi * guy[9] + szi * guz[9]));
                double fkzz = 0.25
                        * ((sxi * gqzz[2] + syi * gqzz[3] + szi * gqzz[4])
                        + (sxi * gux[10] + syi * guy[10] + szi * guz[10]));
                double fkyx = fkxy;
                double fkzx = fkxz;
                double fkzy = fkyz;
                if (i == k) {
                    fix *= 0.5;
                    fiy *= 0.5;
                    fiz *= 0.5;
                    fkx *= 0.5;
                    fky *= 0.5;
                    fkz *= 0.5;
                    fixx *= 0.5;
                    fixy *= 0.5;
                    fixz *= 0.5;
                    fiyx *= 0.5;
                    fiyy *= 0.5;
                    fiyz *= 0.5;
                    fizx *= 0.5;
                    fizy *= 0.5;
                    fizz *= 0.5;
                    fkxx *= 0.5;
                    fkxy *= 0.5;
                    fkxz *= 0.5;
                    fkyx *= 0.5;
                    fkyy *= 0.5;
                    fkyz *= 0.5;
                    fkzx *= 0.5;
                    fkzy *= 0.5;
                    fkzz *= 0.5;
                }

                // Torque due to induced reaction field on permanent dipoles.
                double tix = uyi * fiz - uzi * fiy;
                double tiy = uzi * fix - uxi * fiz;
                double tiz = uxi * fiy - uyi * fix;
                double tkx = uyk * fkz - uzk * fky;
                double tky = uzk * fkx - uxk * fkz;
                double tkz = uxk * fky - uyk * fkx;

                // Torque due to induced reaction field gradient on quadrupoles.
                tix += 2.0 * (qxyi * fixz + qyyi * fiyz + qyzi * fizz - qxzi * fixy - qyzi * fiyy - qzzi * fizy);
                tiy += 2.0 * (qxzi * fixx + qyzi * fiyx + qzzi * fizx - qxxi * fixz - qxyi * fiyz - qxzi * fizz);
                tiz += 2.0 * (qxxi * fixy + qxyi * fiyy + qxzi * fizy - qxyi * fixx - qyyi * fiyx - qyzi * fizx);
                tkx += 2.0 * (qxyk * fkxz + qyyk * fkyz + qyzk * fkzz - qxzk * fkxy - qyzk * fkyy - qzzk * fkzy);
                tky += 2.0 * (qxzk * fkxx + qyzk * fkyx + qzzk * fkzx - qxxk * fkxz - qxyk * fkyz - qxzk * fkzz);
                tkz += 2.0 * (qxxk * fkxy + qxyk * fkyy + qxzk * fkzy - qxyk * fkxx - qyyk * fkyx - qyzk * fkzx);
                tX[i] += lPow * tix;
                tY[i] += lPow * tiy;
                tZ[i] += lPow * tiz;

                final double rx = tkx;
                final double ry = tky;
                final double rz = tkz;
                tkx = rx * transOp[0][0] + ry * transOp[1][0] + rz * transOp[2][0];
                tky = rx * transOp[0][1] + ry * transOp[1][1] + rz * transOp[2][1];
                tkz = rx * transOp[0][2] + ry * transOp[1][2] + rz * transOp[2][2];

                tX[k] += lPow * tkx;
                tY[k] += lPow * tky;
                tZ[k] += lPow * tkz;
                if (lambdaTerm) {
                    ltX[i] += dlPow * tix;
                    ltY[i] += dlPow * tiy;
                    ltZ[i] += dlPow * tiz;
                    ltX[k] += dlPow * tkx;
                    ltY[k] += dlPow * tky;
                    ltZ[k] += dlPow * tkz;
                }
            }
        }
    }

    /**
     * Compute Born radii chain rule terms in parallel via the Grycuk method.
     *
     * @since 1.0
     */
    private class BornCRRegion extends ParallelRegion {

        private final BornCRLoop[] bornCRLoop;

        BornCRRegion(int nt) {
            bornCRLoop = new BornCRLoop[nt];
            for (int i = 0; i < nt; i++) {
                bornCRLoop[i] = new BornCRLoop();
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, bornCRLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing Born radii chain rule term in thread "
                        + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute Born radii chain rule terms for a range of atoms via the
         * Grycuk method.
         *
         * @since 1.0
         */
        private class BornCRLoop extends IntegerForLoop {

            private final double factor = -pow(PI, oneThird) * pow(6.0, (2.0 * oneThird)) / 9.0;
            private final double[] dx_local;
            private double[] gX;
            private double[] gY;
            private double[] gZ;
            private double[] lgX;
            private double[] lgY;
            private double[] lgZ;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            BornCRLoop() {
                dx_local = new double[3];
                /**
                 * Surface area based cavitation energy.
                 */
                double ecav = 0.0;
            }

            @Override
            public void start() {
                int threadID = getThreadIndex();
                gX = grad[threadID][0];
                gY = grad[threadID][1];
                gZ = grad[threadID][2];
                if (lambdaTerm) {
                    lgX = lambdaGrad[threadID][0];
                    lgY = lambdaGrad[threadID][1];
                    lgZ = lambdaGrad[threadID][2];
                }
            }

            /**
             * Use pairwise descreening to compute derivative of the integral of
             * 1/r^6 with respect to r.
             *
             * @param r            separation distance.
             * @param r2           separation distance squared.
             * @param radius       base radius of descreened atom.
             * @param scaledRadius scaled radius descreening atom.
             * @return the derivative.
             */
            private double integralDerivative(double r, double r2, double radius, double scaledRadius) {
                double de = 0.0;
                // Descreen only if the descreened atom does not engulf the descreener.
                if (radius < r + scaledRadius) {
                    // Atom i is engulfed by atom k.
                    if (radius + r < scaledRadius) {
                        double uik = scaledRadius - r;
                        double uik2 = uik * uik;
                        double uik4 = uik2 * uik2;
                        de = -4.0 * PI / uik4;
                    }

                    // Lower integration bound depends on atoms sizes and separation.
                    double sk2 = scaledRadius * scaledRadius;
                    if (radius + r < scaledRadius) {
                        // Atom i is engulfed by atom k.
                        double lik = scaledRadius - r;
                        double lik2 = lik * lik;
                        double lik4 = lik2 * lik2;
                        de = de + 0.25 * PI * (sk2 - 4.0 * scaledRadius * r + 17.0 * r2) / (r2 * lik4);
                    } else if (r < radius + scaledRadius) {
                        // Atoms are overlapped, begin integration from ri.
                        double lik = radius;
                        double lik2 = lik * lik;
                        double lik4 = lik2 * lik2;
                        de = de + 0.25 * PI * (2.0 * radius * radius - sk2 - r2) / (r2 * lik4);
                    } else {
                        // No overlap between atoms.
                        double lik = r - scaledRadius;
                        double lik2 = lik * lik;
                        double lik4 = lik2 * lik2;
                        de = de + 0.25 * PI * (sk2 - 4.0 * scaledRadius * r + r2) / (r2 * lik4);
                    }
                    // Upper integration bound is always the same.
                    double uik = r + scaledRadius;
                    double uik2 = uik * uik;
                    double uik4 = uik2 * uik2;
                    de = de - 0.25 * PI * (sk2 + 4.0 * scaledRadius * r + r2) / (r2 * uik4);
                }

                return de;
            }

            /**
             * Accumulate a contribution to the gradient and dU/dX/dL.
             *
             * @param i  index of atom i.
             * @param k  index of atom k.
             * @param dE partial derivative of the energy with respect to R.
             * @param xr x-component of the separation vector.
             * @param yr y-component of the separation vector.
             * @param zr z-component of the separation vector.
             */
            private void incrementGradient(int i, int k, double dE, double xr, double yr, double zr, double[][] transOp) {
                double dedx = dE * xr;
                double dedy = dE * yr;
                double dedz = dE * zr;
                gX[i] += lPow * dedx;
                gY[i] += lPow * dedy;
                gZ[i] += lPow * dedz;

                final double dedxk = dedx * transOp[0][0] + dedy * transOp[1][0] + dedz * transOp[2][0];
                final double dedyk = dedx * transOp[0][1] + dedy * transOp[1][1] + dedz * transOp[2][1];
                final double dedzk = dedx * transOp[0][2] + dedy * transOp[1][2] + dedz * transOp[2][2];

                gX[k] -= lPow * dedxk;
                gY[k] -= lPow * dedyk;
                gZ[k] -= lPow * dedzk;
                if (lambdaTerm) {
                    lgX[i] += dlPow * dedx;
                    lgY[i] += dlPow * dedy;
                    lgZ[i] += dlPow * dedz;
                    lgX[k] -= dlPow * dedxk;
                    lgY[k] -= dlPow * dedyk;
                    lgZ[k] -= dlPow * dedzk;
                }
            }

            @Override
            public void run(int lb, int ub) {
                int nSymm = crystal.spaceGroup.symOps.size();
                for (int iSymOp = 0; iSymOp < nSymm; iSymOp++) {
                    SymOp symOp = crystal.spaceGroup.symOps.get(iSymOp);
                    double[][] transOp = new double[3][3];
                    double[][] xyz = sXYZ[iSymOp];
                    crystal.getTransformationOperator(symOp, transOp);
                    for (int i = lb; i <= ub; i++) {
                        if (!nativeEnvironmentApproximation && !use[i]) {
                            continue;
                        }
                        final double ri = baseRadius[i];
                        assert (ri > 0.0);
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final double rbi = born[i];
                        double termi = PI4_3 / (rbi * rbi * rbi);
                        termi = factor / pow(termi, (4.0 * oneThird));

                        int[] list = neighborLists[iSymOp][i];
                        for (int k : list) {
                            if (!nativeEnvironmentApproximation && !use[k]) {
                                continue;
                            }
                            final double rk = baseRadius[k];
                            assert (rk > 0.0);
                            if (k != i) {
                                dx_local[0] = xyz[0][k] - xi;
                                dx_local[1] = xyz[1][k] - yi;
                                dx_local[2] = xyz[2][k] - zi;
                                double r2 = crystal.image(dx_local);
                                if (r2 > cut2) {
                                    continue;
                                }
                                final double xr = dx_local[0];
                                final double yr = dx_local[1];
                                final double zr = dx_local[2];
                                final double r = sqrt(r2);

                                // Atom i being descreeened by atom k.
                                final double sk = rk * overlapScale[k];
                                double de = integralDerivative(r, r2, ri, sk);
                                double dbr = termi * de / r;
                                de = dbr * sharedBornGrad.get(i);
                                incrementGradient(i, k, de, xr, yr, zr, transOp);

                                // Atom k being descreeened by atom i.
                                double rbk = born[k];
                                double termk = PI4_3 / (rbk * rbk * rbk);
                                termk = factor / pow(termk, (4.0 * oneThird));

                                final double si = ri * overlapScale[i];
                                de = integralDerivative(r, r2, rk, si);
                                dbr = termk * de / r;
                                de = dbr * sharedBornGrad.get(k);
                                incrementGradient(i, k, de, xr, yr, zr, transOp);

                            } else if (iSymOp > 0) {
                                dx_local[0] = xyz[0][k] - xi;
                                dx_local[1] = xyz[1][k] - yi;
                                dx_local[2] = xyz[2][k] - zi;
                                double r2 = crystal.image(dx_local);
                                if (r2 > cut2) {
                                    continue;
                                }
                                final double xr = dx_local[0];
                                final double yr = dx_local[1];
                                final double zr = dx_local[2];
                                final double r = sqrt(r2);

                                // Atom i being descreeened by atom k.
                                final double sk = rk * overlapScale[k];
                                double de = integralDerivative(r, r2, ri, sk);
                                double dbr = termi * de / r;
                                de = dbr * sharedBornGrad.get(i);
                                incrementGradient(i, k, de, xr, yr, zr, transOp);

                                // For symmetry mates, atom k is not descreeened by atom i.
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Compute Dispersion energy in parallel via pairwise descreening.
     *
     * @since 1.0
     */
    private class DispersionRegion extends ParallelRegion {

        private final DispersionLoop[] dispersionLoop;
        private final SharedDouble sharedDispersion;
        private boolean gradient = false;
        private double[] cdisp;
        private static final double DISP_OVERLAP_SCALE_FACTOR = 0.81;
        private static final double SLEVY = 1.0;
        private static final double AWATER = 0.033428;
        private static final double EPSO = 0.1100;
        private static final double EPSH = 0.0135;
        private static final double RMINO = 1.7025;
        private static final double RMINH = 1.3275;

        DispersionRegion(int nt) {
            dispersionLoop = new DispersionLoop[nt];
            for (int i = 0; i < nt; i++) {
                dispersionLoop[i] = new DispersionLoop();
            }
            sharedDispersion = new SharedDouble();
            cdisp = new double[nAtoms];

            for (int i = 0; i < nAtoms; i++) {
                VDWType type = atoms[i].getVDWType();
                double rmini = type.radius;
                rDisp[i] = rmini / 2.0;
            }
            maxDispersionEnergy();
        }

        public void init() {
            cdisp = new double[nAtoms];

            for (int i = 0; i < nAtoms; i++) {
                VDWType type = atoms[i].getVDWType();
                double rmini = type.radius;
                rDisp[i] = rmini / 2.0;
            }
            maxDispersionEnergy();
        }

        public double getEnergy() {
            return sharedDispersion.get();
        }

        /**
         * Compute the maximum Dispersion energy for each atom in isolation. The
         * loss of dispersion energy due to descreening of other atoms is then
         * calculated in the DispersionLoop.
         */
        private void maxDispersionEnergy() {
            for (int i = 0; i < nAtoms; i++) {
                VDWType type = atoms[i].getVDWType();
                double epsi = type.wellDepth;
                double rmini = type.radius / 2.0;
                if (rDisp[i] > 0.0 && epsi > 0.0) {
                    double sqEpsoEpsi = sqrt(EPSO) + sqrt(epsi);
                    double sqEpshEpsi = sqrt(EPSH) + sqrt(epsi);
                    double emixo = 4.0 * EPSO * epsi / (pow(sqEpsoEpsi, 2));
                    double rmixo = 2.0 * (pow(RMINO, 3) + pow(rmini, 3)) / (pow(RMINO, 2) + pow(rmini, 2));
                    double rmixo3 = pow(rmixo, 3);
                    double rmixo7 = pow(rmixo, 7);
                    double ao = emixo * rmixo7;
                    double emixh = 4.0 * EPSH * epsi / (pow(sqEpshEpsi, 2));
                    double rmixh = 2.0 * (pow(RMINH, 3) + pow(rmini, 3)) / (pow(RMINH, 2) + pow(rmini, 2));
                    double rmixh3 = pow(rmixh, 3);
                    double rmixh7 = pow(rmixh, 7);
                    double ah = emixh * rmixh7;
                    double ri = rDisp[i] + 0.26;
                    double ri3 = pow(ri, 3);
                    double ri7 = pow(ri, 7);
                    double ri11 = pow(ri, 11);
                    if (ri < rmixh) {
                        cdisp[i] = -4.0 * PI * emixh * (rmixh3 - ri3) / 3.0;
                        cdisp[i] = cdisp[i] - emixh * 18.0 / 11.0 * rmixh3 * PI;
                    } else {
                        cdisp[i] = 2.0 * PI * (2.0 * rmixh7 - 11.0 * ri7) * ah;
                        cdisp[i] = cdisp[i] / (11.0 * ri11);
                    }
                    cdisp[i] = 2.0 * cdisp[i];
                    if (ri < rmixo) {
                        cdisp[i] = cdisp[i] - 4.0 * PI * emixo * (rmixo3 - ri3) / 3.0;
                        cdisp[i] = cdisp[i] - emixo * 18.0 / 11.0 * rmixo3 * PI;
                    } else {
                        cdisp[i] = cdisp[i] + 2.0 * PI * (2.0 * rmixo7 - 11.0 * ri7) * ao / (11.0 * ri11);
                    }
                }
                cdisp[i] = SLEVY * AWATER * cdisp[i];
            }
        }

        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }

        @Override
        public void start() {
            sharedDispersion.set(0.0);
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, dispersionLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing Dispersion energy in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute Dispersion energy for a range of atoms via pairwise
         * descreening.
         *
         * @since 1.0
         */
        private class DispersionLoop extends IntegerForLoop {

            private double[] gX;
            private double[] gY;
            private double[] gZ;
            private double[] lgX;
            private double[] lgY;
            private double[] lgZ;
            private double edisp;
            private final double[] dx_local;
            private double r, r2, r3;
            private double xr, yr, zr;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            DispersionLoop() {
                dx_local = new double[3];
            }

            @Override
            public void start() {
                int threadID = getThreadIndex();
                gX = grad[threadID][0];
                gY = grad[threadID][1];
                gZ = grad[threadID][2];
                if (lambdaTerm) {
                    lgX = lambdaGrad[threadID][0];
                    lgY = lambdaGrad[threadID][1];
                    lgZ = lambdaGrad[threadID][2];
                }
                edisp = 0;
            }

            @Override
            public void finish() {
                sharedDispersion.addAndGet(edisp);
            }

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }

                    // Begin with the limit of atom alone in solvent.
                    edisp += cdisp[i];

                    // Now descreen over neighbors.
                    double sum = 0.0;
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    int[] list = neighborLists[0][i];
                    for (int k : list) {
                        final double rk = rDisp[k];
                        if (i != k && rk > 0.0 && use[k]) {
                            dx_local[0] = xi - x[k];
                            dx_local[1] = yi - y[k];
                            dx_local[2] = zi - z[k];
                            r2 = crystal.image(dx_local);
                            if (r2 > cut2) {
                                continue;
                            }
                            xr = dx_local[0];
                            yr = dx_local[1];
                            zr = dx_local[2];
                            r = sqrt(r2);
                            r3 = r * r2;

                            // Atom i descreened by atom k.
                            sum += descreen(i, k);

                            // Flip the sign on {xr, yr, zr};
                            xr = -xr;
                            yr = -yr;
                            zr = -zr;

                            // Atom k descreened by atom i.
                            sum += descreen(k, i);
                        }
                    }

                    // Subtract descreening.
                    edisp -= SLEVY * AWATER * sum;
                }
            }

            private double descreen(int i, int k) {
                double sum = 0.0;
                VDWType type = atoms[i].getVDWType();
                double epsi = type.wellDepth;
                double rmini = type.radius / 2.0;
                double emixo = (4.0 * EPSO * epsi) / (pow(sqrt(EPSO) + sqrt(epsi), 2));
                double rmixo = 2.0 * (pow(RMINO, 3) + pow(rmini, 3)) / (pow(RMINO, 2) + pow(rmini, 2));
                double rmixo7 = pow(rmixo, 7);
                double ao = emixo * rmixo7;
                double emixh = 4.0 * EPSH * epsi / (pow(sqrt(EPSH) + sqrt(epsi), 2));
                double rmixh = 2.0 * (pow(RMINH, 3) + pow(rmini, 3)) / (pow(RMINH, 2) + pow(rmini, 2));
                double rmixh7 = pow(rmixh, 7);
                double ah = emixh * rmixh7;
                double ri = rDisp[i];
                double rk = rDisp[k];
                double sk = rk * DISP_OVERLAP_SCALE_FACTOR;
                double sk2 = sk * sk;
                if (ri < r + sk) {
                    double de = 0.0;
                    double rmax = max(ri, r - sk);
                    double lik = rmax;
                    double lik2 = lik * lik;
                    double lik3 = lik2 * lik;
                    double lik4 = lik3 * lik;
                    if (lik < rmixo) {
                        double uik = min(r + sk, rmixo);
                        double uik2 = uik * uik;
                        double uik3 = uik2 * uik;
                        double uik4 = uik3 * uik;
                        double term = 4.0 * PI / (48.0 * r) * (3.0 * (lik4 - uik4)
                                - 8.0 * r * (lik3 - uik3) + 6.0 * (r2 - sk2) * (lik2 - uik2));
                        double iwca = -emixo * term;
                        sum = sum + iwca;
                        if (gradient) {
                            double dl;
                            if (ri > r - sk) {
                                dl = -lik2 + 2.0 * r2 + 2.0 * sk2;
                                dl = dl * lik2;
                            } else {
                                dl = -lik3 + 4.0 * lik2 * r
                                        - 6.0 * lik * r2
                                        + 2.0 * lik * sk2 + 4.0 * r3
                                        - 4.0 * r * sk2;
                                dl = dl * lik;
                            }
                            double du;
                            if (r + sk > rmixo) {
                                du = -uik2 + 2.0 * r2 + 2.0 * sk2;
                                du = -du * uik2;
                            } else {
                                du = -uik3 + 4.0 * uik2 * r
                                        - 6.0 * uik * r2
                                        + 2.0 * uik * sk2 + 4.0 * r3
                                        - 4.0 * r * sk2;
                                du = -du * uik;
                            }
                            de = de - emixo * PI * (dl + du) / (4.0 * r2);
                        }
                    }
                    if (lik < rmixh) {
                        double uik = min(r + sk, rmixh);
                        double uik2 = uik * uik;
                        double uik3 = uik2 * uik;
                        double uik4 = uik3 * uik;
                        double term = 4.0 * PI / (48.0 * r) * (3.0 * (lik4 - uik4)
                                - 8.0 * r * (lik3 - uik3) + 6.0 * (r2 - sk2) * (lik2 - uik2));
                        double iwca = -2.0 * emixh * term;
                        sum = sum + iwca;
                        if (gradient) {
                            double dl;
                            if (ri > r - sk) {
                                dl = -lik2 + 2.0 * r2 + 2.0 * sk2;
                                dl = dl * lik2;
                            } else {
                                dl = -lik3 + 4.0 * lik2 * r - 6.0 * lik * r2
                                        + 2.0 * lik * sk2 + 4.0 * r3 - 4.0 * r * sk2;
                                dl = dl * lik;
                            }
                            double du;
                            if (r + sk > rmixh) {
                                du = -uik2 + 2.0 * r2 + 2.0 * sk2;
                                du = -du * uik2;
                            } else {
                                du = -uik3 + 4.0 * uik2 * r - 6.0 * uik * r2
                                        + 2.0 * uik * sk2 + 4.0 * r3 - 4.0 * r * sk2;
                                du = -du * uik;
                            }
                            de = de - 2.0 * emixh * PI * (dl + du) / (4.0 * r2);
                        }
                    }
                    double uik = r + sk;
                    double uik2 = uik * uik;
                    double uik3 = uik2 * uik;
                    double uik4 = uik3 * uik;
                    double uik5 = uik4 * uik;
                    double uik6 = uik5 * uik;
                    double uik10 = uik5 * uik5;
                    double uik11 = uik10 * uik;
                    double uik12 = uik11 * uik;
                    double uik13 = uik12 * uik;
                    if (uik > rmixo) {
                        lik = max(rmax, rmixo);
                        lik2 = lik * lik;
                        lik3 = lik2 * lik;
                        lik4 = lik3 * lik;
                        double lik5 = lik4 * lik;
                        double lik6 = lik5 * lik;
                        double lik10 = lik5 * lik5;
                        double lik11 = lik10 * lik;
                        double lik12 = lik11 * lik;
                        double lik13 = lik12 * lik;
                        double term = 4.0 * PI / (120.0 * r * lik5 * uik5) * (15.0 * uik * lik * r * (uik4 - lik4)
                                - 10.0 * uik2 * lik2 * (uik3 - lik3) + 6.0 * (sk2 - r2) * (uik5 - lik5));
                        double term2 = 4.0 * PI / (2640.0 * r * lik12 * uik12) * (120.0 * uik * lik * r * (uik11 - lik11)
                                - 66.0 * uik2 * lik2 * (uik10 - lik10) + 55.0 * (sk2 - r2) * (uik12 - lik12));
                        double idisp = -2.0 * ao * term;
                        double irep = ao * rmixo7 * term2;
                        sum = sum + irep + idisp;
                        if (gradient) {
                            double dl;
                            if (ri > r - sk || rmax < rmixo) {
                                dl = -5.0 * lik2 + 3.0 * r2 + 3.0 * sk2;
                                dl = -dl / lik5;
                            } else {
                                dl = 5.0 * lik3 - 33.0 * lik * r2 - 3.0 * lik * sk2
                                        + 15.0 * (lik2 * r + r3 - r * sk2);
                                dl = dl / lik6;
                            }
                            double du;
                            du = 5.0 * uik3 - 33.0 * uik * r2 - 3.0 * uik * sk2
                                    + 15.0 * (uik2 * r + r3 - r * sk2);
                            du = -du / uik6;
                            de = de - 2.0 * ao * PI * (dl + du) / (15.0 * r2);

                            if (ri > r - sk || rmax < rmixo) {
                                dl = -6.0 * lik2 + 5.0 * r2 + 5.0 * sk2;
                                dl = -dl / lik12;
                            } else {
                                dl = 6.0 * lik3 - 125.0 * lik * r2 - 5.0 * lik * sk2 + 60.0 * (lik2 * r + r3 - r * sk2);
                                dl = dl / lik13;
                            }
                            du = 6.0 * uik3 - 125.0 * uik * r2 - 5.0 * uik * sk2 + 60.0 * (uik2 * r + r3 - r * sk2);
                            du = -du / uik13;
                            de = de + ao * rmixo7 * PI * (dl + du) / (60.0 * r2);
                        }

                    }
                    if (uik > rmixh) {
                        lik = max(rmax, rmixh);
                        lik2 = lik * lik;
                        lik3 = lik2 * lik;
                        lik4 = lik3 * lik;
                        double lik5 = lik4 * lik;
                        double lik6 = lik5 * lik;
                        double lik10 = lik5 * lik5;
                        double lik11 = lik10 * lik;
                        double lik12 = lik11 * lik;
                        double lik13 = lik12 * lik;
                        double term = 4.0 * PI / (120.0 * r * lik5 * uik5) * (15.0 * uik * lik * r * (uik4 - lik4)
                                - 10.0 * uik2 * lik2 * (uik3 - lik3) + 6.0 * (sk2 - r2) * (uik5 - lik5));
                        double term2 = 4.0 * PI / (2640.0 * r * lik12 * uik12) * (120.0 * uik * lik * r * (uik11 - lik11)
                                - 66.0 * uik2 * lik2 * (uik10 - lik10) + 55.0 * (sk2 - r2) * (uik12 - lik12));
                        double idisp = -4.0 * ah * term;
                        double irep = 2.0 * ah * rmixh7 * term2;
                        sum = sum + irep + idisp;
                        if (gradient) {
                            double dl;
                            if (ri > r - sk || rmax < rmixh) {
                                dl = -5.0 * lik2 + 3.0 * r2 + 3.0 * sk2;
                                dl = -dl / lik5;
                            } else {
                                dl = 5.0 * lik3 - 33.0 * lik * r2
                                        - 3.0 * lik * sk2 + 15.0 * (lik2 * r + r3 - r * sk2);
                                dl = dl / lik6;
                            }
                            double du;
                            du = 5.0 * uik3 - 33.0 * uik * r2
                                    - 3.0 * uik * sk2 + 15.0 * (uik2 * r + r3 - r * sk2);
                            du = -du / uik6;
                            de = de - 4.0 * ah * PI * (dl + du) / (15.0 * r2);
                            if (ri > r - sk || rmax < rmixh) {
                                dl = -6.0 * lik2 + 5.0 * r2 + 5.0 * sk2;
                                dl = -dl / lik12;
                            } else {
                                dl = 6.0 * lik3 - 125.0 * lik * r2 - 5.0 * lik * sk2 + 60.0 * (lik2 * r + r3 - r * sk2);
                                dl = dl / lik13;
                            }
                            du = 6.0 * uik3 - 125.0 * uik * r2 - 5.0 * uik * sk2 + 60.0 * (uik2 * r + r3 - r * sk2);
                            du = -du / uik13;
                            de = de + ah * rmixh7 * PI * (dl + du) / (30.0 * r2);
                        }

                    }

                    // Increment the individual dispersion gradient components.
                    if (gradient) {
                        de = -de / r * SLEVY * AWATER;
                        double dedx = de * xr;
                        double dedy = de * yr;
                        double dedz = de * zr;
                        gX[i] += lPow * dedx;
                        gY[i] += lPow * dedy;
                        gZ[i] += lPow * dedz;
                        gX[k] -= lPow * dedx;
                        gY[k] -= lPow * dedy;
                        gZ[k] -= lPow * dedz;
                        if (lambdaTerm) {
                            lgX[i] += dlPow * dedx;
                            lgY[i] += dlPow * dedy;
                            lgZ[i] += dlPow * dedz;
                            lgX[k] -= dlPow * dedx;
                            lgY[k] -= dlPow * dedy;
                            lgZ[k] -= dlPow * dedz;
                        }
                    }
                }
                return sum;
            }
        }
    }

    /**
     * Some static constants.
     */
    private static final double PI4_3 = 4.0 / 3.0 * PI;
    private static final double PI_12 = PI / 12.0;
    /**
     * Constant factor used with quadrupoles.
     */
    private static final double oneThird = 1.0 / 3.0;

    public enum NonPolar {

        CAV, CAV_DISP, GAUSS_DISP, HYDROPHOBIC_PMF, BORN_CAV_DISP, BORN_SOLV, NONE
    }

}
