/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.potential.nonbonded;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.fill;

import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.acos;
import static org.apache.commons.math3.util.FastMath.asin;
import static org.apache.commons.math3.util.FastMath.atan2;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.sin;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.tanh;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.SharedBooleanArray;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.numerics.VectorMath;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.Torsion;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.BioType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ISolvRadType;
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
 * mutlipole force field in parallel using a {@link NeighborList}.
 *
 * @author Michael J. Schnieders<br> derived from:<br> TINKER code by Michael J.
 * Schnieders and Jay W. Ponder<br>
 * @see <a href="http://dx.doi.org/10.1021/ct7001336" target="_blank">M. J.
 * Schnieders and J. W. Ponder, Polarizable atomic multipole solutes in a
 * generalized Kirkwood continuum, Journal of Chemical Theory and Computation
 * 2007, 3, (6), 2083-2097.</a><br>
 *
 */
public class GeneralizedKirkwood implements LambdaInterface {

    private static final Logger logger = Logger.getLogger(GeneralizedKirkwood.class.getName());
    /**
     * Permittivity of water at STP.
     */
    private static final double dWater = 78.3;
    /**
     * Default bondi scale factor.
     */
    private static final double DEFAULT_BONDI_SCALE = 1.15;
    /**
     * Set of force fields for which we have fitted GK/GB radii.
     */
    private static final Set<String> fittedForceFields;

    static {
        fittedForceFields = new HashSet<>();
        /**
         * AMOEBA-PROTEIN-2013 and AMBER99SB as fitted by Stephen D. LuCore.
         * AMBER99SB-TIP3F is an alias for AMBER99SB, just with different
         * explicit-solvent parameters (no re-fitting of GB parameters).
         */
        String[] fitted = {"AMOEBA-PROTEIN-2013", "AMBER99SB", "AMBER99SB-TIP3F"};
        fittedForceFields.addAll(Arrays.asList(fitted));
    }

    /**
     * Kirkwood multipolar reaction field constants.
     */
    private final double fc;
    private final double fd;
    private final double fq;
    /**
     * Empirical constant that controls the GK cross-term.
     */
    private static final double gkc = 2.455;
    /**
     * Empirical scaling of the Bondi radii.
     */
    private final double bondiScale;
    /**
     * Cavitation surface tension coefficient (kcal/mol/A^2).
     */
    private final double surfaceTension;
    /**
     * The requested permittivity.
     */
    private double epsilon = dWater;

    private final double bornaiTerm;
    private final double probe;
    private final double dOffset = 0.09;
    private boolean use[] = null;
    private final Polarization polarization;
    private Atom atoms[];
    private double x[], y[], z[];
    private double globalMultipole[][];
    private double inducedDipole[][];
    private double inducedDipoleCR[][];
    private double baseRadius[];
    private double baseRadiusWithBondi[];
    private double overlapScale[];
    private double rDisp[];
    private double born[];
    private int nAtoms;
    /**
     * This field is because re-initializing the force field resizes some arrays
     * but not others; that second category must, when called on, be resized not
     * to the current number of atoms but to the maximum number of atoms (and
     * thus to the size of the first category of arrays).
     */
    private int maxNumAtoms;
    private final ParticleMeshEwald particleMeshEwald;
    private final ParallelTeam parallelTeam;
    private Crystal crystal;
    private final BornRadiiRegion bornRadiiRegion;
    private final PermanentGKFieldRegion permanentGKFieldRegion;
    private final InducedGKFieldRegion inducedGKFieldRegion;
    private final GKEnergyRegion gkEnergyRegion;
    private final BornCRRegion bornGradRegion;
    private final DispersionRegion dispersionRegion;
    private final CavitationRegion cavitationRegion;

    /**
     * Gradient array for each thread.
     */
    private double grad[][][];
    /**
     * Torque array for each thread.
     */
    private double torque[][][];
    /**
     * Lambda gradient array for each thread (dU/dX/dL)
     */
    private double lambdaGrad[][][];
    /**
     * Lambda torque array for each thread.
     */
    private double lambdaTorque[][][];

    private int neighborLists[][][];

    private SharedDoubleArray sharedBornGrad;
    protected SharedDoubleArray sharedGKField[];
    protected SharedDoubleArray sharedGKFieldCR[];
    private double cutoff;
    private double cut2;
    private final NonPolar nonPolar;

    private boolean lambdaTerm = false;
    private double lambda = 1.0;
    private double lPow = 1.0;
    private double dlPow = 0.0;
    private double dl2Pow = 0.0;

    private double solvationEnergy = 0.0;

    private long gkTime = 0;
    private long pmfTime = 0;
    private long dispersionTime = 0;
    private long cavitationTime = 0;
    private double dispersionEnergy = 0.0;
    private double cavitationEnergy = 0.0;
    /**
     * Use base radii defined by AtomType rather than by atomic number.
     */
    private boolean verboseRadii = false;
    /**
     * If true, prevents Born radii from updating.
     */
    private boolean fixedRadii = false;
    /**
     * Forces all atoms to be considered during Born radius updates.
     */
    private boolean bornUseAll = false;
    /**
     * Provides maps from atomtypes or biotypes to fitted GK radii (by
     * forcefield).
     */
    private boolean useFittedRadii;
    private SolventRadii solventRadii;
    private RADII_MAP_TYPE radiiMapType = RADII_MAP_TYPE.ATOMTYPE;
    /**
     * Maps radii overrides (by AtomType) specified from the command line. e.g.
     * -DradiiOverride=134r1.20,135r1.20 sets atom types 134,135 to Bondi=1.20
     */
    private final HashMap<Integer, Double> radiiOverride = new HashMap<>();
    /**
     * Maps radii overrides (by atom number) specified from the command line.
     * This takes precendence over AtomType-based overrides. e.g.
     * -DradiiOverride=1r1.20,5r1.20 sets atom numbers 1,5 to Bondi=1.20
     */
    private final HashMap<Integer, Double> radiiByNumberMap = new HashMap<>();
    private final ForceField forceField;

    private static final Level GK_WARN_LEVEL;

    static {
        String suppressGKwarnings = System.getProperty("gk-suppressWarnings");
        if (suppressGKwarnings != null && Boolean.parseBoolean(suppressGKwarnings)) {
            GK_WARN_LEVEL = Level.FINE;
        } else {
            GK_WARN_LEVEL = Level.WARNING;
        }
    }

    public double[] getOverlapScale() {
        return overlapScale;
    }

    public double[] getBaseRadii() {
        return baseRadiusWithBondi;
    }

    public double getSurfaceTension() {
        return surfaceTension;
    }

    /**
     * <p>
     * Constructor for GeneralizedKirkwood.</p>
     *
     * @param forceField a {@link ffx.potential.parameters.ForceField} object.
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     * @param particleMeshEwald a
     * {@link ffx.potential.nonbonded.ParticleMeshEwald} object.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
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

        try {
            epsilon = forceField.getDouble(ForceField.ForceFieldDouble.GK_EPSILON);
            logger.info(format(" (GK) GLOBAL dielectric constant set to %.2f", epsilon));
        } catch (Exception e) {
            epsilon = dWater;
        }

        // Se the Kirkwood multipolar reaction field constants.
        fc = 1.0 * (1.0 - epsilon) / (0.0 + 1.0 * epsilon);
        fd = 2.0 * (1.0 - epsilon) / (1.0 + 2.0 * epsilon);
        fq = 3.0 * (1.0 - epsilon) / (2.0 + 3.0 * epsilon);

        String forcefieldName = forceField.getString(ForceField.ForceFieldString.FORCEFIELD,
                ForceField.ForceFieldName.AMOEBA_BIO_2009.toString());
        forcefieldName = forcefieldName.replaceAll("_", "-");
        boolean doUseFitRadii = forceField.getBoolean(ForceField.ForceFieldBoolean.GK_USEFITRADII, true);
        boolean hasFittedRadii = fittedForceFields.contains(forcefieldName.toUpperCase());

        if (doUseFitRadii) {
            if (hasFittedRadii) {
                useFittedRadii = true;
                solventRadii = new SolventRadii(forcefieldName, forceField);
            }
        } else if (hasFittedRadii) {
            logger.log(Level.INFO, String.format(" (GK) Ignoring fitted radii for force field %s", forcefieldName));
        }

        boolean vRadii = forceField.getBoolean(ForceField.ForceFieldBoolean.GK_VERBOSERADII, false);
        if (vRadii) {
            logger.info(" (GK) Verbose radii enabled.");
        }
        verboseRadii = vRadii;

        double bondiScaleValue;
        try {
            bondiScaleValue = forceField.getDouble(ForceField.ForceFieldDouble.GK_BONDIOVERRIDE);
            logger.info(format(" (GK) Scaling GLOBAL bondi radii by factor: %.2f", bondiScaleValue));
        } catch (Exception ex) {
            bondiScaleValue = useFittedRadii ? solventRadii.getDefaultBondi() : DEFAULT_BONDI_SCALE;
            if (verboseRadii) {
                logger.info(format(" (GK) Scaling default GLOBAL bondi radii by factor: %.2f", bondiScaleValue));
            }
        }
        bondiScale = bondiScaleValue;

        String radiiProp = forceField.getString(ForceField.ForceFieldString.GK_RADIIOVERRIDE, null);
        if (radiiProp != null) {
            String tokens[] = radiiProp.split(",");
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

        String radiiByNumber = forceField.getString(ForceField.ForceFieldString.GK_RADIIBYNUMBER, null);
        if (radiiByNumber != null) {
            String tokens[] = radiiByNumber.split(",");
            for (String token : tokens) {
                if (!token.contains("r")) {
                    logger.severe("Invalid radius override.");
                }
                int separator = token.indexOf("r");
                int num = Integer.parseInt(token.substring(0, separator));
                double factor = Double.parseDouble(token.substring(separator + 1));
                logger.info(format(" (GK) Scaling Atom Number %d with bondi factor %.2f", num, factor));
                radiiByNumberMap.put(num, factor);
            }
        }

        NonPolar nonpolarModel = NonPolar.CAV_DISP;
        try {
            String cavModel = forceField.getString(ForceField.ForceFieldString.CAVMODEL, "CAV_DISP").toUpperCase();
            nonpolarModel = getNonPolarModel(cavModel);
            if (nonpolarModel == NonPolar.NONE) {
                logger.info(" No non-polar term will be used.");
            }
        } catch (Exception ex) {
            nonpolarModel = NonPolar.NONE;
            logger.warning(format(" Error parsing non-polar model (set to NONE) %s", ex.toString()));
        }
        nonPolar = nonpolarModel;

        double aiTerm = 4.0 * PI;
        try {
            aiTerm *= forceField.getDouble(ForceField.ForceFieldDouble.BORNAI);
        } catch (Exception ex) {
            switch (nonPolar) {
                case BORN_SOLV:
                    aiTerm *= 0.003; // Value from TINKER.
                    break;
                case BORN_CAV_DISP:
                    aiTerm *= 0.050; // Complete random guess.
                    break;
                default:
                    break;
            }
        }
        bornaiTerm = aiTerm;

        sharedGKField = new SharedDoubleArray[3];
        sharedGKFieldCR = new SharedDoubleArray[3];

        bornUseAll = forceField.getBoolean(ForceField.ForceFieldBoolean.BORN_USE_ALL, false);

        probe = forceField.getDouble(ForceField.ForceFieldDouble.PROBE_RADIUS, 1.4);

        cutoff = particleMeshEwald.getEwaldCutoff();
        cut2 = cutoff * cutoff;

        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);

        initAtomArrays();

        /**
         * If polarization lambda exponent is set to 0.0, then we're running
         * Dual-Topology and the GK energy will be scaled with the overall
         * system lambda value.
         */
        double polLambdaExp = forceField.getDouble(ForceField.ForceFieldDouble.POLARIZATION_LAMBDA_EXPONENT, 3.0);
        if (polLambdaExp == 0.0) {
            lambdaTerm = false;
            logger.info(" GK lambda term set to false.");
        }

        int threadCount = parallelTeam.getThreadCount();
        bornRadiiRegion = new BornRadiiRegion(threadCount);
        permanentGKFieldRegion = new PermanentGKFieldRegion(threadCount);
        inducedGKFieldRegion = new InducedGKFieldRegion(threadCount);
        gkEnergyRegion = new GKEnergyRegion(threadCount);
        bornGradRegion = new BornCRRegion(threadCount);

        double tensionDefault = 0.08;
        switch (nonPolar) {
            case CAV:
                cavitationRegion = new CavitationRegion(threadCount);
                tensionDefault = 0.0049;
                dispersionRegion = null;
                break;
            case CAV_DISP:
                cavitationRegion = new CavitationRegion(threadCount);
                tensionDefault = 0.080;
                dispersionRegion = new DispersionRegion(threadCount);
                break;
            case BORN_CAV_DISP:
                cavitationRegion = null;
                dispersionRegion = new DispersionRegion(threadCount);
                break;
            case HYDROPHOBIC_PMF:
            case BORN_SOLV:
            case NONE:
            default:
                cavitationRegion = null;
                dispersionRegion = null;
                break;
        }

        surfaceTension = forceField.getDouble(ForceField.ForceFieldDouble.SURFACE_TENSION, tensionDefault);

        logger.info("  Continuum Solvation ");
        logger.info(format("   Generalized Kirkwood Cut-Off:       %8.3f (A)", cutoff));
        logger.info(format("   Solvent Dielectric:                 %8.3f", epsilon));
        logger.info(format("   Non-Polar Model:                    %8s",
                nonPolar.toString().replace('_', '-')));

        if (cavitationRegion != null) {
            logger.info(format("   Cavitation Probe Radius:            %8.3f (A)", probe));
            logger.info(format("   Cavitation Surface Tension:         %8.3f (Kcal/mol/A^2)", surfaceTension));
        }

    }

    public NonPolar getNonPolarModel() {
        return nonPolar;
    }

    public void setCutoff(double cutoff) {
        this.cutoff = cutoff;
        this.cut2 = cutoff * cutoff;
    }

    public double getCutoff() {
        return cutoff;
    }

    public void setCrystal(Crystal crystal) {
        this.crystal = crystal;
    }

    public void setNeighborList(int neighbors[][][]) {
        this.neighborLists = neighbors;
    }

    public void setAtoms(Atom atoms[]) {
        this.atoms = atoms;
        nAtoms = atoms.length;
        maxNumAtoms = nAtoms > maxNumAtoms ? nAtoms : maxNumAtoms;
        initAtomArrays();
    }

    public void setFixedRadii(boolean fixedRadii) {
        this.fixedRadii = fixedRadii;
    }

    public boolean getFixedRadii() {
        return fixedRadii;
    }

    public void setBornUseAll(boolean bornUseAll) {
        this.bornUseAll = bornUseAll;
    }

    public boolean getBornUseAll() {
        return bornUseAll;
    }

    private void initAtomArrays() {
        if (fixedRadii) {
            fixedRadii = false;
        }
        x = particleMeshEwald.coordinates[0][0];
        y = particleMeshEwald.coordinates[0][1];
        z = particleMeshEwald.coordinates[0][2];
        globalMultipole = particleMeshEwald.globalMultipole[0];
        inducedDipole = particleMeshEwald.inducedDipole[0];
        inducedDipoleCR = particleMeshEwald.inducedDipoleCR[0];
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
            baseRadiusWithBondi = new double[nAtoms];
            overlapScale = new double[nAtoms];
            rDisp = new double[nAtoms];
            born = new double[nAtoms];
            use = new boolean[nAtoms];
        }

        fill(use, true);
        for (int i = 0; i < nAtoms; i++) {
            baseRadius[i] = 2.0;
//            overlapScale[i] = 0.69;   // Original value based on small molecule parameterization.
            overlapScale[i] = 0.60;     // New default value based on 2016 amino acid GK parameterization.
            if (useFittedRadii) {
                overlapScale[i] = solventRadii.getOverlapScale();
            }
            int atomicNumber = atoms[i].getAtomicNumber();
            AtomType atomType = atoms[i].getAtomType();

            switch (atomicNumber) {
                case 0:
                    baseRadius[i] = 0.0;
                    break;
                case 1:
                    baseRadius[i] = 1.2;
                    break;
                case 2:
                    baseRadius[i] = 1.4;
                    break;
                case 5:
                    baseRadius[i] = 1.8;
                    break;
                case 6:
                    baseRadius[i] = 1.7;
                    break;
                case 7:
                    baseRadius[i] = 1.55;
                    break;
                case 8:
                    baseRadius[i] = 1.52;
                    break;
                case 9:
                    baseRadius[i] = 1.47;
                    break;
                case 10:
                    baseRadius[i] = 1.54;
                    break;
                case 14:
                    baseRadius[i] = 2.1;
                    break;
                case 15:
                    baseRadius[i] = 1.8;
                    break;
                case 16:
                    baseRadius[i] = 1.8;
                    break;
                case 17:
                    baseRadius[i] = 1.75;
                    break;
                case 18:
                    baseRadius[i] = 1.88;
                    break;
                case 34:
                    baseRadius[i] = 1.9;
                    break;
                case 35:
                    baseRadius[i] = 1.85;
                    break;
                case 36:
                    baseRadius[i] = 2.02;
                    break;
                case 53:
                    baseRadius[i] = 1.98;
                    break;
                case 54:
                    baseRadius[i] = 2.16;
                    break;
                default:
                    baseRadius[i] = 2.00;
            }

            double bondiFactor = bondiScale;

            int atomNumber = atoms[i].getIndex() + 1;
            if (useFittedRadii) {
                // First check to see if this atom is in the hardcoded maps.
                switch (radiiMapType) {
                    default:
                    case ATOMTYPE:
                        // Check for hard-coded AtomType bondi factor.
                        if (solventRadii.getAtomBondiMap().containsKey(atomType.type)) {
                            bondiFactor = solventRadii.getAtomBondiMap().get(atomType.type);
                            if (verboseRadii) {
                                logger.info(String.format(" (GK) TypeToBondi: Atom %3s-%-4s (%d) with AtomType %d to %.2f (bondi factor %.4f)",
                                        atoms[i].getResidueName(), atoms[i].getName(), atomNumber, atomType.type, baseRadius[i] * bondiFactor, bondiFactor));
                            }
                        }
                        // TODO: fix these manual overrides which only apply to AM-PRO-13.
                        // Special case for ALA alpha carbon and alpha hydrogen.
                        // Many aminos use types 8,12 for their CA,HA so these can't be specified in the map.
                        if ((atoms[i].getAtomType().type == 8 || atoms[i].getAtomType().type == 12)) {
                            if (atoms[i].getResidueName().equals("ALA")) {
                                bondiFactor = 1.60;
                                if (verboseRadii) {
                                    logger.info(String.format(" (GK) TypeToBondi: Atom %3s-%-4s (%d) with AtomType %d to %.2f (bondi factor %.4f)",
                                            atoms[i].getResidueName(), atoms[i].getName(), atomNumber, atomType.type, baseRadius[i] * bondiFactor, bondiFactor));
                                }
                            }
                        }
                        // Special case for CYD CB,HB.  As above, atom types for the CB+HB are shared between CYS and CYD.
                        // Currently, we only want S+HS on CYS.  So the CB,HB (types 43,44) are disabled in the table.
                        if ((atoms[i].getAtomType().type == 43 || atoms[i].getAtomType().type == 44)) {
                            if (atoms[i].getResidueName().equals("CYD")) {
                                bondiFactor = 1.02;
                                if (verboseRadii) {
                                    logger.info(String.format(" (GK) TypeToBondi: Atom %3s-%-4s (%d) with AtomType %d to %.2f (bondi factor %.4f)",
                                            atoms[i].getResidueName(), atoms[i].getName(), atomNumber, atomType.type, baseRadius[i] * bondiFactor, bondiFactor));
                                }
                            }
                        }
                        break;
                    case BIOTYPE:
                        Map<String, BioType> bioTypes = forceField.getBioTypeMap();
                        BioType bioType = null;
                        for (BioType one : bioTypes.values()) {
                            if (one.atomType == atomType.type) {
                                bioType = one;
                                break;
                            }
                        }
                        if (bioType == null) {
                            logger.severe(String.format("Couldn't find biotype for %s", atomType, toString()));
                        }

                        //                BioType bioType = forceField.getBioType(atoms[i].getResidueName(), atoms[i].getName());
                        //                if (bioType == null) {
                        //                    logger.info(String.format("Null biotype for atom: %3s-%-4s",
                        //                            atoms[i].getResidueName(), atoms[i].getName()));
                        //                }
                        // Check for hard-coded BioType bondi factor.
                        if (solventRadii.getBioBondiMap().containsKey(bioType.index)) {
                            double factor = solventRadii.getBioBondiMap().get(bioType.index);
                            bondiFactor = factor;
                            if (verboseRadii) {
                                logger.info(String.format(" (GK) BiotypeToBondi: Atom %3s-%-4s (%d) with BioType %d to %.2f (bondi factor %.4f)",
                                        atoms[i].getResidueName(), atoms[i].getName(), atomNumber, bioType.index, baseRadius[i] * bondiFactor, bondiFactor));
                            }
                        }
                        // TODO: fix these manual overrides which only apply to AM-PRO-13.
                        // Special case for ALA alpha carbon and alpha hydrogen.
                        // Many aminos use types 8,12 for their CA,HA so these can't be specified in the map.
                        if (bioType.index == 8 || bioType.index == 12) {
                            if (atoms[i].getResidueName().equals("ALA")) {
                                bondiFactor = 1.60;
                                if (verboseRadii) {
                                    logger.info(String.format(" (GK) BiotypeToBondi: Atom %3s-%-4s (%d) with BioType %d to %.2f (bondi factor %.4f)",
                                            atoms[i].getResidueName(), atoms[i].getName(), atomNumber, bioType.index, baseRadius[i] * bondiFactor, bondiFactor));
                                }
                            }
                        }
                        // Special case for CYD CB,HB.  As above, atom types for the CB+HB are shared between CYS and CYD.
                        // Currently, we only want S+HS on CYS.  So the CB,HB (types 43,44) are disabled in the table.
                        if (bioType.index == 83 || bioType.index == 84) {
                            if (atoms[i].getResidueName().equals("CYD")) {
                                bondiFactor = 1.02;
                                if (verboseRadii) {
                                    logger.info(String.format(" (GK) BiotypeToBondi: Atom %3s-%-4s (%d) with BioType %d to %.2f (bondi factor %.4f)",
                                            atoms[i].getResidueName(), atoms[i].getName(), atomNumber, bioType.index, baseRadius[i] * bondiFactor, bondiFactor));
                                }
                            }
                        }
                        break;
                }
            }
            // Next, check if this Atom has an ISolvRad forcefield parameter.
//            if (atoms[i].getISolvRadType() != null) {
//                bondiFactor = atoms[i].getISolvRadType().radiusScale;
//                logger.info(String.format(" (GK) ISolvRad parameter for Atom %3s-%-4s with AtomType %d to %.2f (bondi factor %.2f)",
//                    atoms[i].getResidueName(), atoms[i].getName(), atomType.type, baseRadius[i] * bondiFactor, bondiFactor));
//            }
            ISolvRadType iSolvRadType = forceField.getISolvRadType(Integer.toString(atomType.type));
            if (iSolvRadType != null) {
                bondiFactor = iSolvRadType.radiusScale;
                if (verboseRadii) {
                    logger.info(String.format(" (GK) ISolvRad parameter for Atom %3s-%-4s with AtomType %d to %.2f (bondi factor %.2f)",
                            atoms[i].getResidueName(), atoms[i].getName(), atomType.type, baseRadius[i] * bondiFactor, bondiFactor));
                }
            }
            // Finally, check for command-line bondi factor override.
            if (radiiOverride.containsKey(atomType.type) && !radiiByNumberMap.containsKey(atomNumber)) {
                bondiFactor = radiiOverride.get(atomType.type);
                logger.info(String.format(" (GK) Scaling Atom %3s-%-4s with AtomType %d to %.2f (bondi factor %.2f)",
                        atoms[i].getResidueName(), atoms[i].getName(), atomType.type, baseRadius[i] * bondiFactor, bondiFactor));
            }
            if (radiiByNumberMap.containsKey(atomNumber)) {
                bondiFactor = radiiByNumberMap.get(atomNumber);
                logger.info(String.format(" (GK) Scaling Atom number %d, %3s-%-4s, with factor %.2f",
                        atomNumber, atoms[i].getResidueName(), atoms[i].getName(), bondiFactor));
            }

            baseRadiusWithBondi[i] = baseRadius[i] * bondiFactor;

        }

        // Resets verboseRadii; reduces logging messages when mutating MultiResidues.
        verboseRadii = false;

        if (dispersionRegion != null) {
            dispersionRegion.init();
        }

        if (cavitationRegion != null) {
            cavitationRegion.init();
        }

    }

    public void setUse(boolean use[]) {
        this.use = use;
    }

    /**
     * <p>
     * computeBornRadii</p>
     */
    public void computeBornRadii() {

        /**
         * Born radii are fixed.
         */
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
                if (Double.isInfinite(borni) || Double.isNaN(borni)) {
                    //logger.severe(String.format(" %s\n Born radii %d %8.3f", atoms[i], i, born[i]));
                    throw new EnergyException(String.format(" %s\n Born radii %d %8.3f", atoms[i], i, born[i]), true);
                }
            }
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
     * @param print a boolean.
     * @return a double.
     */
    public double solvationEnergy(boolean gradient, boolean print) {

        /**
         * Initialize the gradient accumulation arrays.
         */
        if (gradient) {
            for (int j = 0; j < nAtoms; j++) {
                sharedBornGrad.set(j, 0.0);
            }
        }

        try {
            /**
             * Find the GK energy.
             */
            gkTime = -System.nanoTime();
            gkEnergyRegion.setGradient(gradient);
            parallelTeam.execute(gkEnergyRegion);
            gkTime += System.nanoTime();
            /**
             * Find the nonpolar energy.
             */
            switch (nonPolar) {
                case CAV:
                    cavitationTime = -System.nanoTime();
                    parallelTeam.execute(cavitationRegion);
                    cavitationTime += System.nanoTime();
                    break;
                case CAV_DISP:
                    dispersionTime = -System.nanoTime();
                    dispersionRegion.setGradient(gradient);
                    parallelTeam.execute(dispersionRegion);
                    dispersionTime += System.nanoTime();
                    cavitationTime = -System.nanoTime();
                    parallelTeam.execute(cavitationRegion);
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

        /**
         * Compute the Born radii chain rule term.
         */
        if (gradient) {
            try {
                gkTime -= System.nanoTime();
                parallelTeam.execute(bornGradRegion);
                gkTime += System.nanoTime();
            } catch (Exception e) {
                String message = "Fatal exception computing Born radii chain rule term.";
                logger.log(Level.SEVERE, message, e);
            }
        }

        if (print) {
            logger.info(format(" Generalized Kirkwood%16.8f %10.3f",
                    gkEnergyRegion.getEnergy(), gkTime * 1e-9));
            switch (nonPolar) {
                case CAV:
                    cavitationEnergy = cavitationRegion.getEnergy();
                    logger.info(format(" Cavitation          %16.8f %10.3f",
                            cavitationEnergy, cavitationTime * 1e-9));
                    break;
                case CAV_DISP:
                    cavitationEnergy = cavitationRegion.getEnergy();
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
                solvationEnergy = gkEnergyRegion.getEnergy() + cavitationRegion.getEnergy();
                break;
            case CAV_DISP:
                solvationEnergy = gkEnergyRegion.getEnergy() + dispersionRegion.getEnergy()
                        + cavitationRegion.getEnergy();
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
     * @param throwError
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
     * @param throwError
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
     * Updates the value of lPow.
     *
     * @param lambda the current lPow value.
     */
    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            this.lambda = lambda;
        } else {
            /**
             * If the lambdaTerm flag is false, lambda must be set to one.
             */
            this.lambda = 1.0;
            lPow = 1.0;
            dlPow = 0.0;
            dl2Pow = 0.0;
        }
    }

    public void setLambdaFunction(double lPow, double dlPow, double dl2Pow) {
        if (lambdaTerm) {
            this.lPow = lPow;
            this.dlPow = dlPow;
            this.dl2Pow = dl2Pow;
        } else {
            /**
             * If the lambdaTerm flag is false, lambda must be set to one.
             */
            this.lambda = 1.0;
            lPow = 1.0;
            dlPow = 0.0;
            dl2Pow = 0.0;
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
     * The 2nd derivative is 0.0. (U=Lambda*Egk, dU/dL=Egk, d2U/dL2=0.0)
     *
     * @return 0.0 is always returned.
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
     * These contributions are already aggregated into the arrays used by PME.
     *
     * @param gradient the Gradient array.
     */
    @Override
    public void getdEdXdL(double[] gradient) {
    }

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

        private final BornRadiiLoop bornRadiiLoop[];
        private SharedDoubleArray sharedBorn;
        private SharedDouble ecavTot;

        public BornRadiiRegion(int nt) {
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
                final double baseRi = baseRadiusWithBondi[i];
                if (!use[i]) {
                    born[i] = baseRi;
                } else {
                    double sum = sharedBorn.get(i);
                    if (sum <= 0.0) {
                        sum = PI4_3 * 1.0e-9;
                        born[i] = 1.0 / pow(sum / PI4_3, THIRD);
                        // logger.info(format(" I < 0; Resetting %d to %12.6f", i, born[i]));
                        continue;
                    }
                    born[i] = 1.0 / pow(sum / PI4_3, THIRD);
                    if (born[i] < baseRi) {
                        // logger.info(format(" Less than base radii; resetting to %d %12.6f", i, baseRi));
                        born[i] = baseRi;
                        continue;
                    }
                    if (Double.isInfinite(born[i]) || Double.isNaN(born[i])) {
                        // logger.info(format(" NaN / Infinite: Resetting Base Radii %d %12.6f", i, baseRi));
                        born[i] = baseRi;
                    }
                }
            }
        }

        /**
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class BornRadiiLoop extends IntegerForLoop {

            private double localBorn[];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;
            private double ecav;

            public BornRadiiLoop() {
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

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    if (!bornUseAll && !use[i]) {
                        continue;
                    }
                    final double baseRi = baseRadiusWithBondi[i];
                    assert (baseRi > 0.0);
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    /**
                     * The descreening integral is initialized to the limit of
                     * the atom alone in solvent.
                     *
                     * Lower values of i may have contributed descreening, so
                     * the integral is incremented rather than initialized.
                     */
                    localBorn[i] += PI4_3 / (baseRi * baseRi * baseRi);
                    int list[] = neighborLists[0][i];
                    int npair = list.length;
                    for (int l = 0; l < npair; l++) {
                        int k = list[l];
                        final double baseRk = baseRadiusWithBondi[k];
                        if (i != k && baseRk > 0.0) {
                            if (!bornUseAll && !use[k]) {
                                continue;
                            }
                            final double xr = x[k] - xi;
                            final double yr = y[k] - yi;
                            final double zr = z[k] - zi;
                            final double r2 = crystal.image(xr, yr, zr);
                            if (r2 > cut2) {
                                continue;
                            }
                            final double r = sqrt(r2);

                            // Atom i being descreeened by atom k.
                            final double scaledRk = baseRk * overlapScale[k];
                            final double scaledRk2 = scaledRk * scaledRk;
                            // Atom i is engulfed by atom k.
                            if (baseRi + r < scaledRk) {
                                final double lower = baseRi;
                                final double upper = scaledRk - r;
                                localBorn[i] += (PI4_3 * (1.0 / (upper * upper * upper) - 1.0 / (lower * lower * lower)));
                            }
                            // Upper integration bound is always the same.
                            double upper = r + scaledRk;
                            // Lower integration bound depends on atoms sizes and separation.
                            double lower;
                            if (baseRi + r < scaledRk) {
                                // Atom i is engulfed by atom k.
                                lower = scaledRk - r;
                            } else if (r < baseRi + scaledRk) {
                                // Atoms are overlapped, begin integration from ri.
                                lower = baseRi;
                            } else {
                                // No overlap between atoms.
                                lower = r - scaledRk;
                            }
                            double l2 = lower * lower;
                            double l4 = l2 * l2;
                            double lr = lower * r;
                            double l4r = l4 * r;
                            double u2 = upper * upper;
                            double u4 = u2 * u2;
                            double ur = upper * r;
                            double u4r = u4 * r;
                            double term = (3.0 * (r2 - scaledRk2) + 6.0 * u2 - 8.0 * ur) / u4r
                                    - (3.0 * (r2 - scaledRk2) + 6.0 * l2 - 8.0 * lr) / l4r;
                            localBorn[i] -= PI_12 * term;

                            // Atom k being descreeened by atom i.
                            final double scaledRi = baseRi * overlapScale[i];
                            final double scaledRi2 = scaledRi * scaledRi;
                            // Atom k is engulfed by atom i.
                            if (baseRk + r < scaledRi) {
                                lower = baseRk;
                                upper = scaledRi - r;
                                localBorn[k] += (PI4_3 * (1.0 / (upper * upper * upper) - 1.0 / (lower * lower * lower)));
                            }
                            // Upper integration bound is always the same.
                            upper = r + scaledRi;
                            // Lower integration bound depends on atoms sizes and separation.
                            if (baseRk + r < scaledRi) {
                                // Atom k is engulfed by atom i.
                                lower = scaledRi - r;
                            } else if (r < baseRk + scaledRi) {
                                // Atoms are overlapped, begin integration from rk.
                                lower = baseRk;
                            } else {
                                // No overlap between atoms.
                                lower = r - scaledRi;
                            }
                            l2 = lower * lower;
                            l4 = l2 * l2;
                            lr = lower * r;
                            l4r = l4 * r;
                            u2 = upper * upper;
                            u4 = u2 * u2;
                            ur = upper * r;
                            u4r = u4 * r;
                            term = (3.0 * (r2 - scaledRi2) + 6.0 * u2 - 8.0 * ur) / u4r
                                    - (3.0 * (r2 - scaledRi2) + 6.0 * l2 - 8.0 * lr) / l4r;
                            localBorn[k] -= PI_12 * term;
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

        private final PermanentGKFieldLoop permanentGKFieldLoop[];

        public PermanentGKFieldRegion(int nt) {
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

            private final double a[][];
            private final double gc[];
            private final double gux[], guy[], guz[];
            private final double gqxx[], gqyy[], gqzz[];
            private final double gqxy[], gqxz[], gqyz[];
            private double fx_local[];
            private double fy_local[];
            private double fz_local[];
            private final double dx_local[];
            private double multipolei[];
            private double xi, yi, zi;
            private double ci, uxi, uyi, uzi, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi;
            private double rbi;
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
                /**
                 * Reduce the field contributions computed by the current thread
                 * into the shared arrays.
                 */
                sharedGKField[0].reduce(fx_local, DoubleOp.SUM);
                sharedGKField[1].reduce(fy_local, DoubleOp.SUM);
                sharedGKField[2].reduce(fz_local, DoubleOp.SUM);
            }

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    xi = x[i];
                    yi = y[i];
                    zi = z[i];
                    multipolei = globalMultipole[i];
                    ci = multipolei[t000];
                    uxi = multipolei[t100];
                    uyi = multipolei[t010];
                    uzi = multipolei[t001];
                    qxxi = multipolei[t200] / 3.0;
                    qxyi = multipolei[t110] / 3.0;
                    qxzi = multipolei[t101] / 3.0;
                    qyyi = multipolei[t020] / 3.0;
                    qyzi = multipolei[t011] / 3.0;
                    qzzi = multipolei[t002] / 3.0;
                    rbi = born[i];
                    int list[] = neighborLists[0][i];
                    int npair = list.length;
                    for (int l = 0; l < npair; l++) {
                        int k = list[l];
                        if (!use[k]) {
                            continue;
                        }
                        permanentGKField(i, k);
                    }
                    /**
                     * Include the self permanent reaction field, which is not
                     * in the neighbor list.
                     */
                    permanentGKField(i, i);
                }
            }

            private void permanentGKField(int i, int k) {
                dx_local[0] = x[k] - xi;
                dx_local[1] = y[k] - yi;
                dx_local[2] = z[k] - zi;
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
                final double multipolek[] = globalMultipole[k];
                final double ck = multipolek[t000];
                final double uxk = multipolek[t100];
                final double uyk = multipolek[t010];
                final double uzk = multipolek[t001];
                final double qxxk = multipolek[t200] / 3.0;
                final double qxyk = multipolek[t110] / 3.0;
                final double qxzk = multipolek[t101] / 3.0;
                final double qyyk = multipolek[t020] / 3.0;
                final double qyzk = multipolek[t011] / 3.0;
                final double qzzk = multipolek[t002] / 3.0;
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
                /**
                 * Reaction potential auxiliary terms.
                 */
                a[0][0] = gf;
                a[1][0] = -gf3;
                a[2][0] = 3.0 * gf5;
                a[3][0] = -15.0 * gf7;
                /**
                 * Reaction potential gradient auxiliary terms.
                 */
                a[0][1] = expc1 * a[1][0];
                a[1][1] = expc1 * a[2][0];
                a[2][1] = expc1 * a[3][0];
                /**
                 * 2nd reaction potential gradient auxiliary terms.
                 */
                a[1][2] = expc1 * a[2][1] + expcdexpc * a[2][0];
                /**
                 * Multiply the potential auxiliary terms by their dielectric
                 * functions.
                 */
                a[0][1] = fc * a[0][1];
                a[1][0] = fd * a[1][0];
                a[1][1] = fd * a[1][1];
                a[1][2] = fd * a[1][2];
                a[2][0] = fq * a[2][0];
                a[2][1] = fq * a[2][1];
                /**
                 * Unweighted reaction potential tensor.
                 */
                gux[1] = xr * a[1][0];
                guy[1] = yr * a[1][0];
                guz[1] = zr * a[1][0];
                /**
                 * Unweighted reaction potential gradient tensor.
                 */
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
                /**
                 * Unweighted 2nd reaction potential gradient tensor.
                 */
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

                /**
                 * Generalized Kirkwood permanent reaction field.
                 */
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
                /**
                 * Scale the self-field by half, such that it sums to one below.
                 */
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

        private final InducedGKFieldLoop inducedGKFieldLoop[];

        public InducedGKFieldRegion(int nt) {
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

            private final double a[][];
            private final double gux[], guy[], guz[];
            private double fx_local[];
            private double fy_local[];
            private double fz_local[];
            private double fxCR_local[];
            private double fyCR_local[];
            private double fzCR_local[];
            private final double dx_local[];
            private double xi, yi, zi;
            private double uix, uiy, uiz;
            private double uixCR, uiyCR, uizCR;
            private double rbi;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public InducedGKFieldLoop() {
                a = new double[3][2];
                gux = new double[5];
                guy = new double[5];
                guz = new double[5];
                dx_local = new double[3];
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
                /**
                 * Reduce the field contributions computed by the current thread
                 * into the shared arrays.
                 */
                sharedGKField[0].reduce(fx_local, DoubleOp.SUM);
                sharedGKField[1].reduce(fy_local, DoubleOp.SUM);
                sharedGKField[2].reduce(fz_local, DoubleOp.SUM);
                sharedGKFieldCR[0].reduce(fxCR_local, DoubleOp.SUM);
                sharedGKFieldCR[1].reduce(fyCR_local, DoubleOp.SUM);
                sharedGKFieldCR[2].reduce(fzCR_local, DoubleOp.SUM);
            }

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    xi = x[i];
                    yi = y[i];
                    zi = z[i];
                    uix = inducedDipole[i][0];
                    uiy = inducedDipole[i][1];
                    uiz = inducedDipole[i][2];
                    uixCR = inducedDipoleCR[i][0];
                    uiyCR = inducedDipoleCR[i][1];
                    uizCR = inducedDipoleCR[i][2];
                    rbi = born[i];
                    int list[] = neighborLists[0][i];
                    int nPair = list.length;
                    for (int l = 0; l < nPair; l++) {
                        int k = list[l];
                        if (!use[k]) {
                            continue;
                        }
                        inducedGKField(i, k);
                    }
                    /**
                     * Include the self induced reaction field, which is not in
                     * the neighbor list.
                     */
                    inducedGKField(i, i);
                }
            }

            private void inducedGKField(int i, int k) {
                dx_local[0] = x[k] - xi;
                dx_local[1] = y[k] - yi;
                dx_local[2] = z[k] - zi;
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
                final double ukx = inducedDipole[k][0];
                final double uky = inducedDipole[k][1];
                final double ukz = inducedDipole[k][2];
                final double ukxCR = inducedDipoleCR[k][0];
                final double ukyCR = inducedDipoleCR[k][1];
                final double ukzCR = inducedDipoleCR[k][2];
                final double rbk = born[k];
                final double rb2 = rbi * rbk;
                final double expterm = exp(-r2 / (gkc * rb2));
                final double expc = expterm / gkc;
                final double expc1 = 1.0 - expc;
                final double gf2 = 1.0 / (r2 + rb2 * expterm);
                final double gf = sqrt(gf2);
                final double gf3 = gf2 * gf;
                final double gf5 = gf3 * gf2;
                /**
                 * Reaction potential auxiliary terms.
                 */
                a[1][0] = -gf3;
                a[2][0] = 3.0 * gf5;
                /**
                 * Reaction potential gradient auxiliary term.
                 */
                a[1][1] = expc1 * a[2][0];
                /**
                 * Multiply the potential auxiliary terms by their dielectric
                 * functions.
                 */
                a[1][0] = fd * a[1][0];
                a[1][1] = fd * a[1][1];
                a[2][0] = fq * a[2][0];
                /**
                 * Unweighted reaction potential gradient tensor.
                 */
                gux[2] = a[1][0] + xr2 * a[1][1];
                gux[3] = xr * yr * a[1][1];
                gux[4] = xr * zr * a[1][1];
                guy[2] = gux[3];
                guy[3] = a[1][0] + yr2 * a[1][1];
                guy[4] = yr * zr * a[1][1];
                guz[2] = gux[4];
                guz[3] = guy[4];
                guz[4] = a[1][0] + zr2 * a[1][1];
                /**
                 * Compute the reaction field due to induced dipoles.
                 */
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
                /**
                 * Scale the self-field by half, such that it sums to one below.
                 */
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
                fx_local[k] += fkx;
                fy_local[k] += fky;
                fz_local[k] += fkz;
                fxCR_local[i] += fixCR;
                fyCR_local[i] += fiyCR;
                fzCR_local[i] += fizCR;
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
        private final GKEnergyLoop gkEnergyLoop[];

        public GKEnergyRegion(int nt) {
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

            private final double a[][];
            private final double b[][];
            private final double gc[];
            private final double gux[], guy[], guz[];
            private final double gqxx[], gqyy[], gqzz[];
            private final double gqxy[], gqxz[], gqyz[];
            private double gb_local[];
            private double gbi_local[];
            private final double dx_local[];
            private double gX[];
            private double gY[];
            private double gZ[];
            private double tX[];
            private double tY[];
            private double tZ[];
            private double lgX[];
            private double lgY[];
            private double lgZ[];
            private double ltX[];
            private double ltY[];
            private double ltZ[];
            private double ci, uxi, uyi, uzi, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi;
            private double ck, uxk, uyk, uzk, qxxk, qxyk, qxzk, qyyk, qyzk, qzzk;
            private double dxi, dyi, dzi, pxi, pyi, pzi, sxi, syi, szi;
            private double dxk, dyk, dzk, pxk, pyk, pzk, sxk, syk, szk;
            private double xr, yr, zr, xr2, yr2, zr2, r2, rbi, rbk;
            private double xi, yi, zi;
            private boolean gradient = false;
            private int count;
            private double gkEnergy;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public GKEnergyLoop() {
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
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    xi = x[i];
                    yi = y[i];
                    zi = z[i];
                    final double multipolei[] = globalMultipole[i];
                    ci = multipolei[t000];
                    uxi = multipolei[t100];
                    uyi = multipolei[t010];
                    uzi = multipolei[t001];
                    qxxi = multipolei[t200] / 3.0;
                    qxyi = multipolei[t110] / 3.0;
                    qxzi = multipolei[t101] / 3.0;
                    qyyi = multipolei[t020] / 3.0;
                    qyzi = multipolei[t011] / 3.0;
                    qzzi = multipolei[t002] / 3.0;
                    dxi = inducedDipole[i][0];
                    dyi = inducedDipole[i][1];
                    dzi = inducedDipole[i][2];
                    pxi = inducedDipoleCR[i][0];
                    pyi = inducedDipoleCR[i][1];
                    pzi = inducedDipoleCR[i][2];
                    sxi = dxi + pxi;
                    syi = dyi + pyi;
                    szi = dzi + pzi;
                    rbi = born[i];
                    int list[] = neighborLists[0][i];
                    int nPair = list.length;
                    for (int l = 0; l < nPair; l++) {
                        int k = list[l];
                        if (!use[k]) {
                            continue;
                        }
                        interaction(i, k);
                    }
                    // Include the self-interaction.
                    interaction(i, i);

                    /**
                     * Formula for Born energy approximation is: e = ai * 4.0*pi
                     * * (ri + probe)^2 * (ri/rb)^6. ri is baseRadius, rb is
                     * Born radius of given atom. ai is an empirical constant
                     * for the atom. If ai is too low, everything wants to pack
                     * into a solid ball, and if ai is too high, everything
                     * wants to unfold and be as solvent-exposed as possible.
                     *
                     * The bornaiTerm is a precalculated 4 * pi * ai value.
                     *
                     * Fix this.
                     */
                    switch (nonPolar) {
                        case BORN_SOLV:
                        case BORN_CAV_DISP:
                            double r = baseRadiusWithBondi[i] + dOffset + probe;
                            double ratio = (baseRadiusWithBondi[i] + dOffset) / born[i];
                            ratio *= ratio;
                            ratio *= (ratio * ratio);
                            double saTerm = -surfaceTension * r * r * ratio;
                            gkEnergy += saTerm / -6.0;
                            // Now calculate derivatives
                            gb_local[i] += saTerm / born[i];
                            break;
                        default:
                            break;
                    }
                }
            }

            @Override
            public void finish() {
                sharedInteractions.addAndGet(count);
                sharedGKEnergy.addAndGet(gkEnergy);
                if (gradient) {
                    /**
                     * Reduce the torque contributions computed by the current
                     * thread into the shared array.
                     */
                    sharedBornGrad.reduce(gb_local, DoubleOp.SUM);
                }
            }

            private void interaction(int i, int k) {
                dx_local[0] = x[k] - xi;
                dx_local[1] = y[k] - yi;
                dx_local[2] = z[k] - zi;
                r2 = crystal.image(dx_local);
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
                final double multipolek[] = globalMultipole[k];
                ck = multipolek[t000];
                uxk = multipolek[t100];
                uyk = multipolek[t010];
                uzk = multipolek[t001];
                qxxk = multipolek[t200] / 3.0;
                qxyk = multipolek[t110] / 3.0;
                qxzk = multipolek[t101] / 3.0;
                qyyk = multipolek[t020] / 3.0;
                qyzk = multipolek[t011] / 3.0;
                qzzk = multipolek[t002] / 3.0;
                dxk = inducedDipole[k][0];
                dyk = inducedDipole[k][1];
                dzk = inducedDipole[k][2];
                pxk = inducedDipoleCR[k][0];
                pyk = inducedDipoleCR[k][1];
                pzk = inducedDipoleCR[k][2];
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
                /**
                 * Reaction potential auxiliary terms.
                 */
                a[0][0] = gf;
                a[1][0] = -gf3;
                a[2][0] = 3.0 * gf5;
                a[3][0] = -15.0 * gf7;
                a[4][0] = 105.0 * gf9;
                a[5][0] = -945.0 * gf11;
                /**
                 * Reaction potential gradient auxiliary terms.
                 */
                a[0][1] = expc1 * a[1][0];
                a[1][1] = expc1 * a[2][0];
                a[2][1] = expc1 * a[3][0];
                a[3][1] = expc1 * a[4][0];
                a[4][1] = expc1 * a[5][0];
                /**
                 * 2nd reaction potential gradient auxiliary terms.
                 */
                a[0][2] = expc1 * a[1][1] + expcdexpc * a[1][0];
                a[1][2] = expc1 * a[2][1] + expcdexpc * a[2][0];
                a[2][2] = expc1 * a[3][1] + expcdexpc * a[3][0];
                a[3][2] = expc1 * a[4][1] + expcdexpc * a[4][0];

                if (gradient) {
                    /**
                     * 3rd reaction potential gradient auxiliary terms.
                     */
                    expcdexpc = 2.0 * expcdexpc;
                    a[0][3] = expc1 * a[1][2] + expcdexpc * a[1][1];
                    a[1][3] = expc1 * a[2][2] + expcdexpc * a[2][1];
                    a[2][3] = expc1 * a[3][2] + expcdexpc * a[3][1];
                    expcdexpc = -expc * dexpc * dexpc;
                    a[0][3] = a[0][3] + expcdexpc * a[1][0];
                    a[1][3] = a[1][3] + expcdexpc * a[2][0];
                    a[2][3] = a[2][3] + expcdexpc * a[3][0];
                    /**
                     * Born radii derivatives of reaction potential auxiliary
                     * terms.
                     */
                    b[0][0] = dgfdr * a[1][0];
                    b[1][0] = dgfdr * a[2][0];
                    b[2][0] = dgfdr * a[3][0];
                    b[3][0] = dgfdr * a[4][0];
                    b[4][0] = dgfdr * a[5][0];
                    /**
                     * Born radii gradients of reaction potential gradient
                     * auxiliary terms.
                     */
                    b[0][1] = b[1][0] - expcr * a[1][0] - expc * b[1][0];
                    b[1][1] = b[2][0] - expcr * a[2][0] - expc * b[2][0];
                    b[2][1] = b[3][0] - expcr * a[3][0] - expc * b[3][0];
                    b[3][1] = b[4][0] - expcr * a[4][0] - expc * b[4][0];
                    /**
                     * Born radii derivatives of the 2nd reaction potential
                     * gradient auxiliary terms.
                     */
                    b[0][2] = b[1][1] - (expcr * (a[1][1] + dexpc * a[1][0])
                            + expc * (b[1][1] + dexpcr * a[1][0] + dexpc * b[1][0]));
                    b[1][2] = b[2][1] - (expcr * (a[2][1] + dexpc * a[2][0])
                            + expc * (b[2][1] + dexpcr * a[2][0] + dexpc * b[2][0]));
                    b[2][2] = b[3][1] - (expcr * (a[3][1] + dexpc * a[3][0])
                            + expc * (b[3][1] + dexpcr * a[3][0] + dexpc * b[3][0]));
                    /**
                     * Multiply the Born radii auxiliary terms by their
                     * dielectric functions.
                     */
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

                /**
                 * Multiply the potential auxiliary terms by their dielectric
                 * functions.
                 */
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
                /**
                 * Compute the GK tensors required to compute the energy.
                 */
                energyTensors();
                /**
                 * Compute the GK interac tion energy.
                 */
                gkEnergy += energy(i, k);
                count++;
                if (gradient || lambdaTerm) {
                    /**
                     * Compute the additional GK tensors required to compute the
                     * energy gradient.
                     */
                    gradientTensors();
                    /**
                     * Compute the permanent GK energy gradient.
                     */
                    permanentEnergyGradient(i, k);
                    if (polarization != Polarization.NONE) {
                        /**
                         * Compute the induced GK energy gradient.
                         */
                        polarizationEnergyGradient(i, k);
                    }
                }
            }

            private void energyTensors() {
                /**
                 * Unweighted reaction potential tensor.
                 */
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
                /**
                 * Unweighted reaction potential gradient tensor.
                 */
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
                /**
                 * Unweighted 2nd reaction potential gradient tensor.
                 */
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
                /**
                 * Electrostatic solvation energy of the permanent multipoles in
                 * their own GK reaction potential.
                 */
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
                    /**
                     * Electrostatic solvation energy of the permanent
                     * multipoles in the GK reaction potential of the induced
                     * dipoles.
                     */
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
                /**
                 * Born radii gradients of unweighted reaction potential tensor.
                 */
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
                /**
                 * Born gradients of the unweighted reaction potential gradient
                 * tensor
                 */
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
                /**
                 * Born radii derivatives of the unweighted 2nd reaction
                 * potential gradient tensor.
                 */
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
                /**
                 * Unweighted 3rd reaction potential gradient tensor.
                 */
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
                if (i == k) {
                    gb_local[i] += drbi;
                    return;
                }
                /**
                 * Increment the gradients and Born chain rule term.
                 */
                final double dedx = dEdX();
                final double dedy = dEdY();
                final double dedz = dEdZ();
                gX[i] -= lPow * dedx;
                gY[i] -= lPow * dedy;
                gZ[i] -= lPow * dedz;
                gb_local[i] += drbi;
                gX[k] += lPow * dedx;
                gY[k] += lPow * dedy;
                gZ[k] += lPow * dedz;
                gb_local[k] += drbk;
                if (lambdaTerm) {
                    lgX[i] -= dlPow * dedx;
                    lgY[i] -= dlPow * dedy;
                    lgZ[i] -= dlPow * dedz;
                    lgX[k] += dlPow * dedx;
                    lgY[k] += dlPow * dedy;
                    lgZ[k] += dlPow * dedz;
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
                /**
                 * Torque on permanent dipoles due to permanent reaction field.
                 */
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
                /**
                 * Torque on quadrupoles due to permanent reaction field
                 * gradient.
                 */
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
                tX[i] += lPow * tix;
                tY[i] += lPow * tiy;
                tZ[i] += lPow * tiz;
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

            private void polarizationEnergyGradient(int i, int k) {
                /**
                 * Electrostatic solvation free energy gradient of the permanent
                 * multipoles in the reaction potential of the induced dipoles.
                 */
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
                /**
                 * Effective radii chain rule terms for the electrostatic
                 * solvation free energy gradient of the permanent multipoles in
                 * the reaction potential of the induced dipoles.
                 */
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
                /**
                 * Increment the gradients and Born chain rule term.
                 */
                if (i == k) {
                    gb_local[i] += dbi;
                } else {
                    gX[i] -= lPow * dpdx;
                    gY[i] -= lPow * dpdy;
                    gZ[i] -= lPow * dpdz;
                    gb_local[i] += dbi;
                    gX[k] += lPow * dpdx;
                    gY[k] += lPow * dpdy;
                    gZ[k] += lPow * dpdz;
                    gb_local[k] += dbk;
                    if (lambdaTerm) {
                        lgX[i] -= dlPow * dpdx;
                        lgY[i] -= dlPow * dpdy;
                        lgZ[i] -= dlPow * dpdz;
                        lgX[k] += dlPow * dpdx;
                        lgY[k] += dlPow * dpdy;
                        lgZ[k] += dlPow * dpdz;
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
                /**
                 * Torque due to induced reaction field on permanent dipoles.
                 */
                double tix = uyi * fiz - uzi * fiy;
                double tiy = uzi * fix - uxi * fiz;
                double tiz = uxi * fiy - uyi * fix;
                double tkx = uyk * fkz - uzk * fky;
                double tky = uzk * fkx - uxk * fkz;
                double tkz = uxk * fky - uyk * fkx;
                /**
                 * Torque due to induced reaction field gradient on quadrupoles.
                 */
                tix += 2.0 * (qxyi * fixz + qyyi * fiyz + qyzi * fizz - qxzi * fixy - qyzi * fiyy - qzzi * fizy);
                tiy += 2.0 * (qxzi * fixx + qyzi * fiyx + qzzi * fizx - qxxi * fixz - qxyi * fiyz - qxzi * fizz);
                tiz += 2.0 * (qxxi * fixy + qxyi * fiyy + qxzi * fizy - qxyi * fixx - qyyi * fiyx - qyzi * fizx);
                tkx += 2.0 * (qxyk * fkxz + qyyk * fkyz + qyzk * fkzz - qxzk * fkxy - qyzk * fkyy - qzzk * fkzy);
                tky += 2.0 * (qxzk * fkxx + qyzk * fkyx + qzzk * fkzx - qxxk * fkxz - qxyk * fkyz - qxzk * fkzz);
                tkz += 2.0 * (qxxk * fkxy + qxyk * fkyy + qxzk * fkzy - qxyk * fkxx - qyyk * fkyx - qyzk * fkzx);
                tX[i] += lPow * tix;
                tY[i] += lPow * tiy;
                tZ[i] += lPow * tiz;
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

        private final BornCRLoop bornCRLoop[];
        private final SharedDouble ecavTot;

        public BornCRRegion(int nt) {
            bornCRLoop = new BornCRLoop[nt];
            for (int i = 0; i < nt; i++) {
                bornCRLoop[i] = new BornCRLoop();
            }
            ecavTot = new SharedDouble(0.0);
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

            private final double factor = -pow(PI, THIRD) * pow(6.0, (2.0 * THIRD)) / 9.0;
            private final double dx_local[];
            private double gX[];
            private double gY[];
            private double gZ[];
            private double lgX[];
            private double lgY[];
            private double lgZ[];
            private double ecav;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public BornCRLoop() {
                dx_local = new double[3];
                ecav = 0.0;
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

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    final double ri = baseRadiusWithBondi[i];
                    assert (ri > 0.0);
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double rbi = born[i];
                    double termi = PI4_3 / (rbi * rbi * rbi);
                    termi = factor / pow(termi, (4.0 * THIRD));
                    int list[] = neighborLists[0][i];
                    int nPair = list.length;
                    for (int l = 0; l < nPair; l++) {
                        int k = list[l];
                        if (!use[k]) {
                            continue;
                        }
                        final double rk = baseRadiusWithBondi[k];
                        if (k != i && rk > 0.0) {
                            dx_local[0] = x[k] - xi;
                            dx_local[1] = y[k] - yi;
                            dx_local[2] = z[k] - zi;
                            double r2 = crystal.image(dx_local);
                            if (r2 > cut2) {
                                continue;
                            }
                            final double xr = dx_local[0];
                            final double yr = dx_local[1];
                            final double zr = dx_local[2];

                            // Atom i being descreeened by atom k.
                            final double sk = rk * overlapScale[k];
                            final double sk2 = sk * sk;
                            final double r = sqrt(r2);
                            double de = 0.0;
                            // Atom i is engulfed by atom k.
                            if (ri + r < sk) {
                                final double uik = sk - r;
                                de = -4.0 * PI / pow(uik, 4);
                            }
                            // Lower integration bound depends on atoms sizes and separation.
                            if (ri + r < sk) {
                                // Atom i is engulfed by atom k.
                                final double lik = sk - r;
                                double lik4 = pow(lik, 4);
                                de = de + 0.25 * PI * (sk2 - 4.0 * sk * r + 17.0 * r2) / (r2 * lik4);
                            } else if (r < ri + sk) {
                                // Atoms are overlapped, begin integration from ri.
                                double lik = ri;
                                double lik4 = pow(lik, 4);
                                de = de + 0.25 * PI * (2.0 * ri * ri - sk2 - r2) / (r2 * lik4);
                            } else {
                                // No overlap between atoms.
                                double lik = r - sk;
                                double lik4 = pow(lik, 4);
                                de = de + 0.25 * PI * (sk2 - 4.0 * sk * r + r2) / (r2 * lik4);
                            }
                            // Upper integration bound is always the same.
                            double uik = r + sk;
                            double uik4 = pow(uik, 4);
                            de = de - 0.25 * PI * (sk2 + 4.0 * sk * r + r2) / (r2 * uik4);

                            double dbr = termi * de / r;
                            de = dbr * sharedBornGrad.get(i);
                            /**
                             * Increment the overall derivatives.
                             */
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

                            // Atom k being descreeened by atom i.
                            double rbk = born[k];
                            double termk = PI4_3 / (rbk * rbk * rbk);
                            termk = factor / pow(termk, (4.0 * THIRD));
                            final double si = ri * overlapScale[i];
                            final double si2 = si * si;
                            de = 0.0;
                            // Atom i is engulfed by atom k.
                            if (rk + r < si) {
                                uik = si - r;
                                de = -4.0 * PI / pow(uik, 4);
                            }
                            // Lower integration bound depends on atoms sizes and separation.
                            if (rk + r < si) {
                                // Atom i is engulfed by atom k.
                                final double lik = si - r;
                                double lik4 = pow(lik, 4);
                                de = de + 0.25 * PI * (si2 - 4.0 * si * r + 17.0 * r2) / (r2 * lik4);
                            } else if (r < rk + si) {
                                // Atoms are overlapped, begin integration from ri.
                                double lik = rk;
                                double lik4 = pow(lik, 4);
                                de = de + 0.25 * PI * (2.0 * rk * rk - si2 - r2) / (r2 * lik4);
                            } else {
                                // No overlap between atoms.
                                double lik = r - si;
                                double lik4 = pow(lik, 4);
                                de = de + 0.25 * PI * (si2 - 4.0 * si * r + r2) / (r2 * lik4);
                            }
                            // Upper integration bound is always the same.
                            uik = r + si;
                            uik4 = pow(uik, 4);
                            de = de - 0.25 * PI * (si2 + 4.0 * si * r + r2) / (r2 * uik4);
                            dbr = termk * de / r;
                            de = dbr * sharedBornGrad.get(k);
                            /**
                             * Increment the overall derivatives.
                             */
                            dedx = de * xr;
                            dedy = de * yr;
                            dedz = de * zr;
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

        private final DispersionLoop dispersionLoop[];
        private final SharedDouble sharedDispersion;
        private boolean gradient = false;
        private double[] cdisp = null;
        private static final double DISP_OVERLAP_SCALE_FACTOER = 0.81;
        private static final double SLEVY = 1.0;
        private static final double AWATER = 0.033428;
        private static final double EPSO = 0.1100;
        private static final double EPSH = 0.0135;
        private static final double RMINO = 1.7025;
        private static final double RMINH = 1.3275;

        public DispersionRegion(int nt) {
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

            private double gX[];
            private double gY[];
            private double gZ[];
            private double lgX[];
            private double lgY[];
            private double lgZ[];
            private double edisp;
            private final double dx_local[];
            private double r, r2, r3;
            private double xr, yr, zr;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public DispersionLoop() {
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
                    /**
                     * Begin with the limit of atom alone in solvent.
                     */
                    edisp += cdisp[i];

                    /**
                     * Now descreen over neighbors.
                     */
                    double sum = 0.0;
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    int list[] = neighborLists[0][i];
                    int npair = list.length;
                    for (int l = 0; l < npair; l++) {
                        int k = list[l];
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
                            /**
                             * Atom i descreened by atom k.
                             */
                            sum += descreen(i, k);
                            /**
                             * Flip the sign on {xr, yr, zr};
                             */
                            xr = -xr;
                            yr = -yr;
                            zr = -zr;
                            /**
                             * Atom k descreened by atom i.
                             */
                            sum += descreen(k, i);
                        }
                    }
                    /**
                     * Subtract descreening.
                     */
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
                double sk = rk * DISP_OVERLAP_SCALE_FACTOER;
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
                    /**
                     * Increment the individual dispersion gradient components.
                     */
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
     * Compute Cavitation energy in parallel.
     *
     * @since 1.0
     */
    private class CavitationRegion extends ParallelRegion {

        private final AtomOverlapLoop atomOverlapLoop[];
        private final CavitationLoop cavitationLoop[];
        private final InitLoop initLoop[];
        private final SharedDouble sharedCavitation;
        private final int maxarc = 150;
        private final static double delta = 1.0e-8;
        private final static double delta2 = delta * delta;

        private double xc1[][];
        private double yc1[][];
        private double zc1[][];
        private double dsq1[][];
        private double bsq1[][];
        private double b1[][];
        private IndexedDouble gr[][];
        private int intag1[][];
        private Integer count[];
        private boolean buried[];
        private SharedBooleanArray skip;
        private double area[];
        private double r[];
        private long initTime = 0;
        private long overlapTime = 0;
        private long cavTime = 0;

        public CavitationRegion(int nt) {
            atomOverlapLoop = new AtomOverlapLoop[nt];
            cavitationLoop = new CavitationLoop[nt];
            initLoop = new InitLoop[nt];
            for (int i = 0; i < nt; i++) {
                atomOverlapLoop[i] = new AtomOverlapLoop();
                cavitationLoop[i] = new CavitationLoop();
                initLoop[i] = new InitLoop();
            }
            sharedCavitation = new SharedDouble();
            init();
        }

        public double getEnergy() {
            return sharedCavitation.get();
        }

        public final void init() {
            if (count == null || count.length < nAtoms) {
                count = new Integer[nAtoms];
                xc1 = new double[nAtoms][maxarc];
                yc1 = new double[nAtoms][maxarc];
                zc1 = new double[nAtoms][maxarc];
                dsq1 = new double[nAtoms][maxarc];
                bsq1 = new double[nAtoms][maxarc];
                b1 = new double[nAtoms][maxarc];
                gr = new IndexedDouble[nAtoms][maxarc];
                intag1 = new int[nAtoms][maxarc];
                buried = new boolean[nAtoms];
                skip = new SharedBooleanArray(nAtoms);
                area = new double[nAtoms];
                r = new double[nAtoms];
            }
            /**
             * Set the sphere radii.
             */
            for (int i = 0; i < nAtoms; i++) {
                VDWType type = atoms[i].getVDWType();
                double rmini = type.radius;
                r[i] = rmini / 2.0;
                if (r[i] != 0.0) {
                    r[i] = r[i] + probe;
                }
            }
            for (int i = 0; i < nAtoms; i++) {
                skip.set(i, true);
            }
        }

        @Override
        public void start() {
            sharedCavitation.set(0.0);
        }

        @Override
        public void finish() {
            if (logger.isLoggable(Level.FINE)) {
                int n = initLoop.length;
                initTime = 0;
                overlapTime = 0;
                cavTime = 0;
                for (int i = 0; i < n; i++) {
                    initTime = max(initLoop[i].time, initTime);
                    overlapTime = max(atomOverlapLoop[i].time, overlapTime);
                    cavTime = max(cavitationLoop[i].time, cavTime);
                }
                logger.fine(String.format(" Cavitation Init: %10.3f Overlap: %10.3f Cav:  %10.3f",
                        initTime * 1e-9, overlapTime * 1e-9, cavTime * 1e-9));
            }
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, initLoop[getThreadIndex()]);
                execute(0, nAtoms - 1, atomOverlapLoop[getThreadIndex()]);
                execute(0, nAtoms - 1, cavitationLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing Cavitation energy in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        private class IndexedDouble implements Comparable {

            public double value;
            public int key;

            public IndexedDouble(double value, int key) {
                this.value = value;
                this.key = key;
            }

            @Override
            public int compareTo(Object o) {
                if (!(o instanceof IndexedDouble)) {
                    return 0;
                }
                IndexedDouble d = (IndexedDouble) o;
                if (value < d.value) {
                    return -1;
                } else if (value == d.value) {
                    return 0;
                } else {
                    return 1;
                }
            }
        }

        /**
         * Initialize arrays for Cavitation calculation.
         */
        private class InitLoop extends IntegerForLoop {

            public long time;

            @Override
            public void start() {
                time = -System.nanoTime();
            }

            @Override
            public void finish() {
                time += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) throws Exception {
                for (int i = lb; i <= ub; i++) {
                    buried[i] = false;
                    area[i] = 0.0;
                    count[i] = 0;
                    double xr = x[i];
                    double yr = y[i];
                    double zr = z[i];
                    double rri = r[i];
                    final int list[] = neighborLists[0][i];
                    int npair = list.length;
                    for (int l = 0; l < npair; l++) {
                        int k = list[l];
                        double rrik = rri + r[k];
                        double dx = x[k] - xr;
                        double dy = y[k] - yr;
                        double dz = z[k] - zr;
                        double ccsq = dx * dx + dy * dy + dz * dz;
                        if (ccsq <= rrik * rrik) {
                            if (use[i] || use[k]) {
                                skip.set(k, false);
                                skip.set(i, false);
                            }
                        }
                    }
                }
            }
        }

        /**
         * Compute Cavitation energy for a range of atoms.
         *
         * @since 1.0
         */
        private class AtomOverlapLoop extends IntegerForLoop {

            public long time;

            @Override
            public void start() {
                time = -System.nanoTime();
            }

            @Override
            public void finish() {
                time += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                /**
                 * Find overlaps with the current sphere.
                 */
                for (int i = lb; i <= ub; i++) {
                    if (skip.get(i) || !use[i]) {
                        continue;
                    }
                    int list[] = neighborLists[0][i];
                    int npair = list.length;
                    for (int l = 0; l < npair; l++) {
                        int k = list[l];
                        if (k == i) {
                            continue;
                        }
                        pair(i, k);
                        if (!skip.get(k) || !use[k]) {
                            pair(k, i);
                        }
                    }
                }
            }

            private void pair(int i, int k) {
                double xi = x[i];
                double yi = y[i];
                double zi = z[i];
                double rri = r[i];
                double rri2 = 2.0 * rri;
                double rplus = rri + r[k];
                double dx = x[k] - xi;
                double dy = y[k] - yi;
                double dz = z[k] - zi;
                if (abs(dx) >= rplus || abs(dy) >= rplus || abs(dz) >= rplus) {
                    return;
                }

                /**
                 * Check for overlap of spheres by testing center to center
                 * distance against sum and difference of radii.
                 */
                double xysq = dx * dx + dy * dy;
                if (xysq < delta2) {
                    dx = delta;
                    dy = 0.0;
                    xysq = delta2;
                }
                double r2 = xysq + dz * dz;
                double dr = sqrt(r2);
                if (rplus - dr <= delta) {
                    return;
                }
                double rminus = rri - r[k];
                /**
                 * Calculate overlap parameters between "i" and "ir" sphere.
                 */
                synchronized (count) {
                    /**
                     * Check for a completely buried "ir" sphere.
                     */
                    if (dr - abs(rminus) <= delta) {
                        if (rminus <= 0.0) {
                            // SA for this atom is zero.
                            buried[i] = true;
                            return;
                        }
                        return;
                    }
                    int n = count[i];
                    xc1[i][n] = dx;
                    yc1[i][n] = dy;
                    zc1[i][n] = dz;
                    dsq1[i][n] = xysq;
                    bsq1[i][n] = r2;
                    b1[i][n] = dr;
                    gr[i][n] = new IndexedDouble((r2 + rplus * rminus) / (rri2 * b1[i][n]), n);
                    intag1[i][n] = k;
                    count[i]++;
                    if (count[i] >= maxarc) {
                        //logger.severe(String.format(" Increase the value of MAXARC to (%d).", count[i]));
                        throw new EnergyException(String.format(" Increase the value of MAXARC to (%d).", count[i]), false);
                    }
                }
            }
        }

        /**
         * Compute Cavitation energy for a range of atoms.
         *
         * @since 1.0
         */
        private class CavitationLoop extends IntegerForLoop {

            private double thec = 0;
            private IndexedDouble arci[];
            private final double dArea[][];
            private final double ldArea[][];
            private boolean omit[];
            private double xc[];
            private double yc[];
            private double zc[];
            private double dsq[];
            private double b[];
            private double bsq[];
            private double bg[];
            private double risq[];
            private double ri[];
            private double ther[];
            private double ider[];
            private double sign_yder[];
            private double arcf[];
            private double ex[];
            private double ux[];
            private double uy[];
            private double uz[];
            private int kent[];
            private int kout[];
            private int intag[];
            private int lt[];
            private int i;
            private int j;
            private int ib;
            private double ecav;
            /**
             * Set pi multiples, overlap criterion and tolerances.
             */
            private final static double pix2 = 2.0 * PI;
            private final static double pix4 = 4.0 * PI;
            private final static double pid2 = PI / 2.0;
            private final static double eps = 1.0e-8;

            public long time;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public CavitationLoop() {
                dArea = new double[3][];
                ldArea = new double[3][];
                allocateMemory(maxarc);
            }

            private void allocateMemory(int maxarc) {
                arci = new IndexedDouble[maxarc];
                arcf = new double[maxarc];
                risq = new double[maxarc];
                ri = new double[maxarc];
                dsq = new double[maxarc];
                bsq = new double[maxarc];
                intag = new int[maxarc];
                lt = new int[maxarc];
                kent = new int[maxarc];
                kout = new int[maxarc];
                ider = new double[maxarc];
                sign_yder = new double[maxarc];
                xc = new double[maxarc];
                yc = new double[maxarc];
                zc = new double[maxarc];
                b = new double[maxarc];
                omit = new boolean[maxarc];
                bg = new double[maxarc];
                ther = new double[maxarc];
                ex = new double[maxarc];
                ux = new double[maxarc];
                uy = new double[maxarc];
                uz = new double[maxarc];
            }

            @Override
            public void start() {
                time = -System.nanoTime();
                int threadID = getThreadIndex();
                dArea[0] = grad[threadID][0];
                dArea[1] = grad[threadID][1];
                dArea[2] = grad[threadID][2];
                if (lambdaTerm) {
                    ldArea[0] = lambdaGrad[threadID][0];
                    ldArea[1] = lambdaGrad[threadID][1];
                    ldArea[2] = lambdaGrad[threadID][2];
                }
                ecav = 0;
                fill(ider, 0);
                fill(sign_yder, 0);
            }

            @Override
            public void finish() {
                sharedCavitation.addAndGet(ecav);
                time += System.nanoTime();
            }

            @Override
            public void run(int lb, int ub) {
                /**
                 * Compute the area and derivatives of current "ir" sphere
                 */
                for (int ir = lb; ir <= ub; ir++) {
                    if (skip.get(ir) || !use[ir]) {
                        continue;
                    }
                    double xi = x[ir];
                    double yi = y[ir];
                    double zi = z[ir];
                    double rri = r[ir];
                    double rri2 = 2.0 * rri;
                    double rrisq = rri * rri;
                    double wght = surfaceTension;
                    boolean moved = false;
                    surface(xi, yi, zi, rri, rri2, rrisq, wght, moved, ir);
                    if (area[ir] < 0.0) {
                        logger.log(GK_WARN_LEVEL, String.format(" Negative surface area set to 0 for atom %d.", ir));
                        //logger.warning(String.format(" Negative surface area set to 0 for atom %d.", ir));
                        area[ir] = 0.0;
                        /**
                         * xi = xi + rmove; yi = yi + rmove; zi = zi + rmove;
                         * moved = true; surface(xi, yi, zi, rri, rri2, rrisq,
                         * wght, moved, ir); if (area[ir] < 0.0) {
                         * logger.warning(String.format(" Negative surface area
                         * set to 0 for atom %d.", ir)); area[ir] = 0.0;
                         * continue; }
                         */
                    }
                    area[ir] *= rrisq * wght;
                    ecav += area[ir];
                }
            }

            public void surface(double xi, double yi, double zi, double rri, double rri2,
                    double rrisq, double wght, boolean moved, int ir) {

                ib = 0;
                int jb = 0;
                double arclen = 0.0;
                double exang = 0.0;

                // Case where no other spheres overlap the current sphere.
                if (count[ir] == 0) {
                    area[ir] = pix4;
                    return;
                }
                // Case where only one sphere overlaps the current sphere.
                if (count[ir] == 1) {
                    int k = 0;
                    double txk = xc1[ir][0];
                    double tyk = yc1[ir][0];
                    double tzk = zc1[ir][0];
                    double bsqk = bsq1[ir][0];
                    double bk = b1[ir][0];
                    intag[0] = intag1[ir][0];
                    double arcsum = pix2;
                    ib = ib + 1;
                    arclen += gr[ir][k].value * arcsum;
                    if (!moved) {
                        int in = intag[k];
                        double t1 = arcsum * rrisq * (bsqk - rrisq + r[in] * r[in])
                                / (rri2 * bsqk * bk);
                        dArea[0][ir] -= lPow * txk * t1 * wght;
                        dArea[1][ir] -= lPow * tyk * t1 * wght;
                        dArea[2][ir] -= lPow * tzk * t1 * wght;
                        dArea[0][in] += lPow * txk * t1 * wght;
                        dArea[1][in] += lPow * tyk * t1 * wght;
                        dArea[2][in] += lPow * tzk * t1 * wght;
                        if (lambdaTerm) {
                            ldArea[0][ir] -= dlPow * txk * t1 * wght;
                            ldArea[1][ir] -= dlPow * tyk * t1 * wght;
                            ldArea[2][ir] -= dlPow * tzk * t1 * wght;
                            ldArea[0][in] += dlPow * txk * t1 * wght;
                            ldArea[1][in] += dlPow * tyk * t1 * wght;
                            ldArea[2][in] += dlPow * tzk * t1 * wght;
                        }
                    }
                    area[ir] = ib * pix2 + exang + arclen;
                    area[ir] = area[ir] % pix4;
                    return;
                }
                /**
                 * general case where more than one sphere intersects the
                 * current sphere; sort intersecting spheres by their degree of
                 * overlap with the current main sphere
                 */
                Arrays.sort(gr[ir], 0, count[ir]);
                for (int j = 0; j < count[ir]; j++) {
                    int k = gr[ir][j].key;
                    intag[j] = intag1[ir][k];
                    xc[j] = xc1[ir][k];
                    yc[j] = yc1[ir][k];
                    zc[j] = zc1[ir][k];
                    dsq[j] = dsq1[ir][k];
                    b[j] = b1[ir][k];
                    bsq[j] = bsq1[ir][k];
                    omit[j] = false;
                }
                /**
                 * Radius of the each circle on the surface of the "ir" sphere.
                 */
                for (int i = 0; i < count[ir]; i++) {
                    double gi = gr[ir][i].value * rri;
                    bg[i] = b[i] * gi;
                    risq[i] = rrisq - gi * gi;
                    ri[i] = sqrt(risq[i]);
                    ther[i] = pid2 - asin(min(1.0, max(-1.0, gr[ir][i].value)));
                }
                /**
                 * Find boundary of inaccessible area on "ir" sphere.
                 */
                for (int k = 0; k < count[ir] - 1; k++) {
                    if (omit[k]) {
                        continue;
                    }
                    double txk = xc[k];
                    double tyk = yc[k];
                    double tzk = zc[k];
                    double bk = b[k];
                    double therk = ther[k];
                    for (j = k + 1; j < count[ir]; j++) {
                        if (omit[j]) {
                            continue;
                        }
                        /*
                         * Check to see if J circle is intersecting K circle;
                         * get distance between circle centers and sum of radii.
                         */
                        double cc = (txk * xc[j] + tyk * yc[j] + tzk * zc[j])
                                / (bk * b[j]);
                        // Check acos FORTRAN vs. Java.
                        cc = acos(min(1.0, max(-1.0, cc)));
                        double td = therk + ther[j];
                        // Check to see if circles enclose separate regions
                        if (cc >= td) {
                            continue;
                        }
                        // Check for circle J completely inside circle K
                        if (cc + ther[j] < therk) {
                            omit[j] = true;
                            continue;
                        }
                        // Check for circles essentially parallel.
                        if (cc > delta) {
                            if (pix2 - cc <= td) {
                                area[ir] = 0.0;
                                return;
                            }
                        }
                    }
                }

                /**
                 * Find T value of circle intersections.
                 */
                for (int k = 0; k < count[ir]; k++) {
                    if (omit[k]) {
                        continue; // goto 110
                    }
                    boolean komit = omit[k];
                    omit[k] = true;
                    int narc = 0;
                    boolean top = false;
                    double txk = xc[k];
                    double tyk = yc[k];
                    double tzk = zc[k];
                    double dk = sqrt(dsq[k]);
                    double bsqk = bsq[k];
                    double bk = b[k];
                    double gk = gr[ir][k].value * rri;
                    double risqk = risq[k];
                    double rik = ri[k];
                    double therk = ther[k];
                    /**
                     * Rotation matrix elements.
                     */
                    double t1 = tzk / (bk * dk);
                    double axx = txk * t1;
                    double axy = tyk * t1;
                    double axz = dk / bk;
                    double ayx = tyk / dk;
                    double ayy = txk / dk;
                    double azx = txk / bk;
                    double azy = tyk / bk;
                    double azz = tzk / bk;
                    for (int l = 0; l < count[ir]; l++) {
                        if (omit[l]) {
                            continue;
                        }
                        double txl = xc[l];
                        double tyl = yc[l];
                        double tzl = zc[l];
                        /**
                         * Rotate spheres so K vector collinear with z-axis.
                         */
                        double uxl = txl * axx + tyl * axy - tzl * axz;
                        double uyl = tyl * ayy - txl * ayx;
                        double uzl = txl * azx + tyl * azy + tzl * azz;
                        double cosine = min(1.0, max(-1.0, uzl / b[l]));
                        if (acos(cosine) < therk + ther[l]) {
                            double dsql = uxl * uxl + uyl * uyl;
                            double tb = uzl * gk - bg[l];
                            double txb = uxl * tb;
                            double tyb = uyl * tb;
                            double td = rik * dsql;
                            double tr2 = risqk * dsql - tb * tb;
                            tr2 = max(eps, tr2);
                            double tr = sqrt(tr2);
                            double txr = uxl * tr;
                            double tyr = uyl * tr;
                            /**
                             * Get T values of intersection for K circle.
                             */
                            tb = (txb + tyr) / td;
                            tb = min(1.0, max(-1.0, tb));
                            double tk1 = acos(tb);
                            if (tyb - txr < 0.0) {
                                tk1 = pix2 - tk1;
                            }
                            tb = (txb - tyr) / td;
                            tb = min(1.0, max(-1.0, tb));
                            double tk2 = acos(tb);
                            if (tyb + txr < 0.0) {
                                tk2 = pix2 - tk2;
                            }
                            thec = (rrisq * uzl - gk * bg[l])
                                    / (rik * ri[l] * b[l]);
                            double the = 0.0;
                            if (abs(thec) < 1.0) {
                                the = -acos(thec);
                            } else if (thec >= 1.0) {
                                the = 0.0;
                            } else if (thec <= -1.0) {
                                the = -PI;
                            }
                            /**
                             * See if "tk1" is entry or exit point; check t=0
                             * point; "ti" is exit point, "tf" is entry point.
                             */
                            cosine = min(1.0, max(-1.0, (uzl * gk - uxl * rik)
                                    / (b[l] * rri)));
                            double ti, tf;
                            if ((acos(cosine) - ther[l]) * (tk2 - tk1) <= 0.0) {
                                ti = tk2;
                                tf = tk1;
                            } else {
                                ti = tk1;
                                tf = tk2;
                            }
                            narc += 1;
                            if (narc > maxarc) {
                                /*logger.severe(String.
                                        format(" Increase value of MAXARC %d.", narc));*/
                                throw new EnergyException(String.format(" Increase value of MAXARC %d.", narc), false);
                            }
                            int narc1 = narc - 1;
                            if (tf <= ti) {
                                arcf[narc1] = tf;
                                arci[narc1] = new IndexedDouble(0.0, narc1);
                                tf = pix2;
                                lt[narc1] = l;
                                ex[narc1] = the;
                                top = true;
                                narc = narc + 1;
                                narc1 = narc - 1;
                            }
                            arcf[narc1] = tf;
                            arci[narc1] = new IndexedDouble(ti, narc1);
                            lt[narc1] = l;
                            ex[narc1] = the;
                            ux[l] = uxl;
                            uy[l] = uyl;
                            uz[l] = uzl;
                        }
                    }
                    omit[k] = komit;
                    /**
                     * Special case; K circle without intersections.
                     */
                    if (narc <= 0) {
                        double arcsum = pix2;
                        ib += 1;
                        arclen += gr[ir][k].value * arcsum;
                        if (!moved) {
                            int in = intag[k];
                            t1 = arcsum * rrisq * (bsqk - rrisq + r[in] * r[in])
                                    / (rri2 * bsqk * bk);
                            dArea[0][ir] -= lPow * txk * t1 * wght;
                            dArea[1][ir] -= lPow * tyk * t1 * wght;
                            dArea[2][ir] -= lPow * tzk * t1 * wght;
                            dArea[0][in] += lPow * txk * t1 * wght;
                            dArea[1][in] += lPow * tyk * t1 * wght;
                            dArea[2][in] += lPow * tzk * t1 * wght;
                            if (lambdaTerm) {
                                ldArea[0][ir] -= dlPow * txk * t1 * wght;
                                ldArea[1][ir] -= dlPow * tyk * t1 * wght;
                                ldArea[2][ir] -= dlPow * tzk * t1 * wght;
                                ldArea[0][in] += dlPow * txk * t1 * wght;
                                ldArea[1][in] += dlPow * tyk * t1 * wght;
                                ldArea[2][in] += dlPow * tzk * t1 * wght;
                            }
                        }
                        continue;
                    }
                    /**
                     * General case; sum up arclength and set connectivity code.
                     */
                    Arrays.sort(arci, 0, narc);
                    double arcsum = arci[0].value;
                    int mi = arci[0].key;
                    double t = arcf[mi];
                    int ni = mi;
                    for (j = 1; j < narc; j++) {
                        int m = arci[j].key;
                        if (t < arci[j].value) {
                            arcsum += (arci[j].value - t);
                            exang += ex[ni];
                            jb += 1;
                            if (jb >= maxarc) {
                                /*logger.severe(String.
                                        format("Increase the value of MAXARC (%d).", jb));*/
                                throw new EnergyException(String.format("Increase the value of MAXARC (%d).", jb), false);
                            }
                            int l = lt[ni];
                            ider[l] += 1;
                            sign_yder[l] += 1;
                            kent[jb] = maxarc * (l + 1) + (k + 1);
                            l = lt[m];
                            ider[l] += 1;
                            sign_yder[l] -= 1;
                            kout[jb] = maxarc * (k + 1) + (l + 1);
                        }
                        double tt = arcf[m];
                        if (tt >= t) {
                            t = tt;
                            ni = m;
                        }
                    }
                    arcsum += (pix2 - t);
                    if (!top) {
                        exang += ex[ni];
                        jb = jb + 1;
                        int l = lt[ni];
                        ider[l] += 1;
                        sign_yder[l] += 1;
                        kent[jb] = maxarc * (l + 1) + (k + 1);
                        l = lt[mi];
                        ider[l] += 1;
                        sign_yder[l] -= 1;
                        kout[jb] = maxarc * (k + 1) + (l + 1);
                    }
                    /**
                     * Calculate the surface area derivatives.
                     */
                    for (int l = 0; l <= count[ir]; l++) {
                        if (ider[l] == 0) {
                            continue;
                        }
                        double rcn = ider[l] * rrisq;
                        ider[l] = 0;
                        double uzl = uz[l];
                        double gl = gr[ir][l].value * rri;
                        double bgl = bg[l];
                        double bsql = bsq[l];
                        double risql = risq[l];
                        double wxlsq = bsql - uzl * uzl;
                        double wxl = sqrt(wxlsq);
                        double p = bgl - gk * uzl;
                        double v = risqk * wxlsq - p * p;
                        v = max(eps, v);
                        v = sqrt(v);
                        t1 = rri * (gk * (bgl - bsql) + uzl * (bgl - rrisq))
                                / (v * risql * bsql);
                        double deal = -wxl * t1;
                        double decl = -uzl * t1 - rri / v;
                        double dtkal = (wxlsq - p) / (wxl * v);
                        double dtkcl = (uzl - gk) / v;
                        double s = gk * b[l] - gl * uzl;
                        t1 = 2.0 * gk - uzl;
                        double t2 = rrisq - bgl;
                        double dtlal = -(risql * wxlsq * b[l] * t1 - s * (wxlsq * t2 + risql * bsql))
                                / (risql * wxl * bsql * v);
                        double dtlcl = -(risql * b[l] * (uzl * t1 - bgl) - uzl * t2 * s)
                                / (risql * bsql * v);
                        double gaca = rcn * (deal - (gk * dtkal - gl * dtlal) / rri) / wxl;
                        double gacb = (gk - uzl * gl / b[l]) * sign_yder[l] * rri / wxlsq;
                        sign_yder[l] = 0;
                        if (!moved) {
                            double faca = ux[l] * gaca - uy[l] * gacb;
                            double facb = uy[l] * gaca + ux[l] * gacb;
                            double facc = rcn * (decl - (gk * dtkcl - gl * dtlcl) / rri);
                            double dax = axx * faca - ayx * facb + azx * facc;
                            double day = axy * faca + ayy * facb + azy * facc;
                            double daz = azz * facc - axz * faca;
                            int in = intag[l];
                            dArea[0][ir] += lPow * dax * wght;
                            dArea[1][ir] += lPow * day * wght;
                            dArea[2][ir] += lPow * daz * wght;
                            dArea[0][in] -= lPow * dax * wght;
                            dArea[1][in] -= lPow * day * wght;
                            dArea[2][in] -= lPow * daz * wght;
                            if (lambdaTerm) {
                                ldArea[0][ir] += dlPow * dax * wght;
                                ldArea[1][ir] += dlPow * day * wght;
                                ldArea[2][ir] += dlPow * daz * wght;
                                ldArea[0][in] -= dlPow * dax * wght;
                                ldArea[1][in] -= dlPow * day * wght;
                                ldArea[2][in] -= dlPow * daz * wght;
                            }
                        }

                    }
                    arclen += gr[ir][k].value * arcsum;
                    if (!moved) {
                        int in = intag[k];
                        t1 = arcsum * rrisq * (bsqk - rrisq + r[in] * r[in])
                                / (rri2 * bsqk * bk);
                        dArea[0][ir] -= lPow * txk * t1 * wght;
                        dArea[1][ir] -= lPow * tyk * t1 * wght;
                        dArea[2][ir] -= lPow * tzk * t1 * wght;
                        dArea[0][in] += lPow * txk * t1 * wght;
                        dArea[1][in] += lPow * tyk * t1 * wght;
                        dArea[2][in] += lPow * tzk * t1 * wght;
                        if (lambdaTerm) {
                            ldArea[0][ir] -= dlPow * txk * t1 * wght;
                            ldArea[1][ir] -= dlPow * tyk * t1 * wght;
                            ldArea[2][ir] -= dlPow * tzk * t1 * wght;
                            ldArea[0][in] += dlPow * txk * t1 * wght;
                            ldArea[1][in] += dlPow * tyk * t1 * wght;
                            ldArea[2][in] += dlPow * tzk * t1 * wght;
                        }
                    }
                }
                if (arclen == 0.0) {
                    area[ir] = 0.0;
                    return;
                }
                if (jb == 0) {
                    area[ir] = ib * pix2 + exang + arclen;
                    area[ir] = area[ir] % pix4;
                    return;
                }
                /**
                 * Find number of independent boundaries and check connectivity.
                 */
                j = 0;
                for (int k = 1; k <= jb; k++) {
                    if (kout[k] == -1) {
                        continue;
                    }
                    i = k;
                    boolean success = independentBoundaries(k, exang, jb, ir, arclen);
                    if (success) {
                        return;
                    }
                }
                ib = ib + 1;
                /*
                 if (moved) {
                 logger.warning(String.format(" Connectivity error at atom %d.", ir));
                 } else {
                 */
                logger.log(GK_WARN_LEVEL, String.format(" Connectivity error at atom %d", ir));
                //logger.warning(String.format(" Connectivity error at atom %d.", ir));
                area[ir] = 0.0;
                /*
                 moved = true;
                 xi += rmove;
                 yi += rmove;
                 zi += rmove;
                 surface(xi, yi, zi, rri, rri2, rrisq, wght, moved, ir);
                 } */
            }

            /**
             * Find number of independent boundaries and check connectivity.
             * This method may set the "goto160" flag.
             *
             * @param k
             * @param exang
             * @param jb
             * @param ir
             * @param arclen
             */
            public boolean independentBoundaries(int k, double exang,
                    int jb, int ir, double arclen) {
                int m = kout[i];
                kout[i] = -1;
                j = j + 1;
                for (int ii = 1; ii <= jb; ii++) {
                    if (m == kent[ii]) {
                        if (ii == k) {
                            ib++;
                            if (j == jb) {
                                area[ir] = ib * 2.0 * PI + exang + arclen;
                                area[ir] = area[ir] % (4.0 * PI);
                                return true;
                            }
                            return false;
                        }
                        i = ii;
                        return independentBoundaries(k, exang, jb, ir, arclen);
                    }
                }
                return false;
            }
        }
    }

    /**
     * Compute Volume energy in parallel.
     *
     * @since 1.0
     */
    private class VolumeRegion extends ParallelRegion {

        private final VolumeLoop volumeLoop[];
        private final SharedDouble sharedVolume;
        private final int itab[];
        private final static int MAXCUBE = 40;
        private final static int MAXARC = 1000;
        private final static int MAXMNB = 500;
        /**
         * maximum number of cycle convex edges.
         */
        private final static int MAXCYEP = 30;
        /**
         * maximum number of convex face cycles
         */
        private final static int MAXFPCY = 18;
        /**
         * Radius of a water molecular.
         */
        private final static double exclude = 1.4;
        /**
         * maximum number of saddle faces.
         */
        private final int maxfs = 3 * nAtoms;
        /**
         * maximum number of convex edges.
         */
        private final int maxep = 5 * nAtoms;
        /**
         * maximum number of neighboring atom pairs.
         */
        private final int maxcls = 50 * nAtoms;
        /**
         * maximum number of circles.
         */
        private final int maxc = 5 * nAtoms;
        /**
         * maximum number of total tori.
         */
        private final int maxt = 3 * nAtoms;
        /**
         * maximum number of temporary tori.
         */
        private final int maxtt = 25 * nAtoms;
        /**
         * maximum number of concave edges.
         */
        private final int maxen = 5 * nAtoms;
        /**
         * maximum number of vertices.
         */
        private final int maxv = 5 * nAtoms;
        /**
         * maximum number of probe positions.
         */
        private final int maxp = 3 * nAtoms;
        /**
         * maximum number of concave faces.
         */
        private final int maxfn = 2 * nAtoms;
        /**
         * maximum number of convex faces.
         */
        private final int maxfp = nAtoms;
        /**
         * maximum number of cycles.
         */
        private final int maxcy = nAtoms;
        /**
         * Atomic radii.
         */
        private final double radius[] = new double[nAtoms];
        /**
         * If true, atom is not used.
         */
        private final boolean skip[] = new boolean[nAtoms];
        /**
         * Copy of the atomic coordinates. [X,Y,Z][Atom Number]
         */
        private final double a[][] = new double[3][nAtoms];
        /**
         * If true, atom has no free surface.
         */
        private final boolean nosurf[] = new boolean[nAtoms];
        /**
         * Atom free of neighbors.
         */
        private final boolean afree[] = new boolean[nAtoms];
        /**
         * Atom buried.
         */
        private final boolean abur[] = new boolean[nAtoms];
        /**
         * Begin and end pointers for atoms neighbors.
         */
        private final int acls[][] = new int[2][nAtoms];
        /**
         * Atom numbers of neighbors.
         */
        private final int cls[] = new int[maxcls];
        /**
         * Pointer from neighbor to torus.
         */
        private final int clst[] = new int[maxcls];
        /**
         * Number of temporary tori.
         */
        private int ntt;
        /**
         * Temporary torus atom numbers.
         */
        private final int tta[][] = new int[2][maxtt];
        /**
         * First edge of each temporary torus.
         */
        private final int ttfe[] = new int[maxtt];
        /**
         * Last edge of each temporary torus.
         */
        private final int ttle[] = new int[maxtt];
        /**
         * Pointer to next edge of temporary torus.
         */
        private final int enext[] = new int[maxen];
        /**
         * Temporary torus buried.
         */
        private final boolean ttbur[] = new boolean[maxtt];
        /**
         * Temporary torus free.
         */
        private final boolean ttfree[] = new boolean[maxtt];
        /**
         * Torus center.
         */
        private final double t[][] = new double[3][maxt];
        /**
         * Torus radius.
         */
        private final double tr[] = new double[maxt];
        /**
         * Torus axis.
         */
        private final double tax[][] = new double[3][maxt];
        /**
         * Number of tori.
         */
        private int nt;
        /**
         * Torus atom number.
         */
        private final int ta[][] = new int[2][maxt];
        /**
         * Torus first edge.
         */
        int tfe[] = new int[maxt];
        /**
         * Torus free edge of neighbor.
         */
        boolean tfree[] = new boolean[maxt];
        /**
         * Probe position coordinates.
         */
        private final double p[][] = new double[3][maxp];
        /**
         * Number of probe positions.
         */
        private int np;
        /**
         * Probe position atom numbers.
         */
        private final int pa[][] = new int[3][maxp];
        /**
         * Vertex coordinates.
         */
        private final double v[][] = new double[3][maxv];
        /**
         * Number of vertices.
         */
        private int nv;
        /**
         * Vertex atom number.
         */
        private final int va[] = new int[maxv];
        /**
         * Vertex probe number.
         */
        private final int vp[] = new int[maxv];
        /**
         * Number of concave edges.
         */
        private int nen;
        /**
         * Vertex numbers for each concave edge.
         */
        private final int env[][] = new int[2][maxen];
        /**
         * Concave face concave edge numbers.
         */
        private final int fnen[][] = new int[3][maxfn];
        /**
         * Circle center.
         */
        private final double c[][] = new double[3][maxc];
        /**
         * Circle radius.
         */
        private final double cr[] = new double[maxc];
        /**
         * Number of circles.
         */
        private int nc;
        /**
         * Circle atom number.
         */
        private final int ca[] = new int[maxc];
        /**
         * Circle torus number.
         */
        private final int ct[] = new int[maxc];
        /**
         * Number of convex edges.
         */
        private int nep;
        /**
         * Convex edge circle number.
         */
        private final int epc[] = new int[maxep];
        /**
         * Convex edge vertex number.
         */
        private final int epv[][] = new int[2][maxep];
        /**
         * First convex edge of each atom.
         */
        private final int afe[] = new int[nAtoms];
        /**
         * Last convex edge of each atom.
         */
        private final int ale[] = new int[nAtoms];
        /**
         * Pointer to next convex edge of atom.
         */
        private final int epnext[] = new int[maxep];
        /**
         * Number of saddle faces.
         */
        private int nfs;
        /**
         * Saddle face concave edge numbers.
         */
        private final int fsen[][] = new int[2][maxfs];
        /**
         * Saddle face convex edge numbers.
         */
        private final int fsep[][] = new int[2][maxfs];
        /**
         * Number of cycles.
         */
        private int ncy;
        /**
         * Number of convex edges in cycle.
         */
        private final int cynep[] = new int[maxcy];
        /**
         * Cycle convex edge numbers.
         */
        private final int cyep[][] = new int[MAXCYEP][maxcy];
        /**
         * Number of convex faces.
         */
        private int nfp;
        /**
         * Atom number of convex face.
         */
        private final int fpa[] = new int[maxfp];
        /**
         * Convex face cycle numbers
         */
        private final int fpcy[][] = new int[MAXFPCY][maxfp];
        /**
         * Number of cycles bounding convex face.
         */
        private final int fpncy[] = new int[maxfp];

        // These are from the "nearby" method.
        /**
         * True if cube contains active atoms.
         */
        private final boolean activeCube[][][] = new boolean[MAXCUBE][MAXCUBE][MAXCUBE];
        /**
         * True if cube or adjacent cubes have active atoms.
         */
        private final boolean activeAdjacentCube[][][] = new boolean[MAXCUBE][MAXCUBE][MAXCUBE];
        /**
         * Pointer to first atom in list for cube.
         */
        private final int firstAtomPointer[][][] = new int[MAXCUBE][MAXCUBE][MAXCUBE];
        /**
         * Integer cube coordinates.
         */
        private final int cubeCoordinates[][] = new int[3][nAtoms];
        /**
         * Pointer to next atom in cube.
         */
        private final int nextAtomPointer[] = new int[nAtoms];

        /**
         * VolumeRegion constructor.
         *
         * @param nt Number of threads.
         */
        public VolumeRegion(int nt) {
            volumeLoop = new VolumeLoop[nt];
            for (int i = 0; i < nt; i++) {
                volumeLoop[i] = new VolumeLoop();
            }
            sharedVolume = new SharedDouble();
            itab = new int[nAtoms];

            /**
             * Set atom coordinates and radii, the excluded buffer radius
             * ("exclude") is added to atomic radii.
             */
            for (int i = 0; i < nAtoms; i++) {
                if (radius[i] == 0.0) {
                    skip[i] = true;
                } else {
                    radius[i] += exclude;
                    skip[i] = false;
                }
            }

        }

        public double getEnergy() {
            return sharedVolume.get();
        }

        @Override
        public void start() {
            sharedVolume.set(0.0);
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

        private final static double SIZE = 0.000001;
        private final double vector[] = new double[3];

        public void wiggle() {
            /**
             * Apply a small perturbation of fixed magnitude to each atom.
             */
            for (int i = 0; i < nAtoms; i++) {
                getRandomVector(vector);
                a[0][i] = x[i] + (SIZE * vector[0]);
                a[1][i] = y[i] + (SIZE * vector[1]);
                a[2][i] = z[i] + (SIZE * vector[2]);
            }
        }

        public void getRandomVector(double vector[]) {
            double x, y, s;
            x = 0;
            y = 0;

            /**
             * Get a pair of appropriate components in the plane.
             */
            s = 2.0;
            while (s >= 1.0) {
                x = (2.0 * Math.random()) - 1.0;
                y = (2.0 * Math.random()) - 1.0;
                s = (x * x) + (y * y);
            }
            /**
             * Construct the 3-dimensional random unit vector.
             */
            vector[2] = 1.0 - 2.0 * s;
            s = 2.0 * sqrt(1.0 - s);
            vector[1] = s * y;
            vector[0] = s * x;
        }

        /**
         * Compute Volume energy for a range of atoms.
         *
         * @since 1.0
         */
        private class VolumeLoop extends IntegerForLoop {

            private int nfn;
            private int l, l1, l2;
            private int io, ir, in, iv;
            private int narc, nx, ny, nz;
            private int istart, istop;
            private int jstart, jstop;
            private int kstart, kstop;
            private int mstart, mstop;
            private int isum, itemp, tcube;
            private final int mxcube = 15;
            private final int inov[] = new int[MAXARC];
            private final int cube[][][][] = new int[2][mxcube][mxcube][mxcube];
            private double evol;
            private double xmin, ymin, zmin;
            private double xmax, ymax, zmax;
            private double aa, bb, temp, phi_term;
            private double theta1, theta2, dtheta;
            private double seg_dx, seg_dy, seg_dz;
            private double pre_dx, pre_dy, pre_dz;
            private double rinsq, rdiff;
            private double rsecn, rsec2n;
            private double cosine, ti, tf;
            private double alpha, beta;
            private double ztop, zstart;
            private double ztopshave;
            private double phi1, cos_phi1;
            private double phi2, cos_phi2;
            private double zgrid, pix2;
            private double rsec2r, rsecr;
            private double rr, rrx2, rrsq;
            private double rmax, edge;
            private double xr, yr, zr;
            private double dist2, vdwsum;
            private double zstep;
            private final double arci[] = new double[MAXARC];
            private final double arcf[] = new double[MAXARC];
            private final double dx[] = new double[MAXARC];
            private final double dy[] = new double[MAXARC];
            private final double dsq[] = new double[MAXARC];
            private final double d[] = new double[MAXARC];
            private boolean ttok, cinsp, cintp;
            private final double vdwrad[] = new double[nAtoms];
            private final double dex[][] = new double[3][nAtoms];

            /**
             * Extra padding to avert cache interface.
             */
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void start() {
                fill(dex[0], 0.0);
                fill(dex[1], 0.0);
                fill(dex[2], 0.0);
                evol = 0.0;
            }

            @Override
            public void finish() {
                sharedVolume.addAndGet(evol);
            }

            public void setRadius() {
                /**
                 * Initialize minimum and maximum range of atoms.
                 */
                pix2 = 2.0 * PI;
                rmax = 0.0;
                xmin = x[0];
                xmax = x[0];
                ymin = y[0];
                ymax = y[0];
                zmin = z[0];
                zmax = z[0];

                /**
                 * Assign van der Waals radii to the atoms; note that the radii
                 * are incremented by the size of the probe; then get the
                 * maximum and minimum ranges of atoms.
                 */
                for (int i = 0; i < nAtoms; i++) {
                    radius[i] = atoms[i].getVDWType().radius / 2.0;
                    vdwrad[i] = radius[i];
                    if (vdwrad[i] == 0.0) {
                        skip[i] = true;
                    } else {
                        skip[i] = false;
                        vdwrad[i] += probe;
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
            public void calcVolume() {
                double volume = 0;
                double area = 0;
                nearby();
                torus();
                place();
                compress();
                saddles();
                contact();
                vam(volume, area);
            }

            public void getVector(double ai[], double temp[][], int index) {
                ai[0] = temp[0][index];
                ai[1] = temp[1][index];
                ai[2] = temp[2][index];
            }

            public void getVector(double ai[], double temp[][][], int index1, int index2) {
                ai[0] = temp[0][index1][index2];
                ai[1] = temp[1][index1][index2];
                ai[2] = temp[2][index1][index2];
            }

            /**
             * The gettor method tests for a possible torus position at the
             * interface between two atoms, and finds the torus radius, center
             * and axis.
             */
            public boolean gettor(int ia, int ja, double torcen[], double torad[], double torax[]) {

                double dij, temp;
                double temp1, temp2;
                double vij[] = new double[3];
                double uij[] = new double[3];
                double bij[] = new double[3];
                double ai[] = new double[3];
                double aj[] = new double[3];

                /**
                 * Get the distance between the two atoms.
                 */
                ttok = false;
                getVector(ai, a, ia);
                getVector(aj, a, ja);
                dij = VectorMath.dist(ai, aj);

                /**
                 * Find a unit vector along interatomic (torus) axis.
                 */
                for (int k = 0; k < 3; k++) {
                    vij[k] = a[k][ja] - a[k][ia];
                    uij[k] = vij[k] / dij;
                }

                /**
                 * Find coordinates of the center of the torus.
                 */
                temp = 1.0 + ((radius[ia] + probe) * (radius[ia] + probe) - (radius[ja] + probe) * (radius[ja] + probe)) / (dij * dij);
                for (int k = 0; k < 3; k++) {
                    bij[k] = a[k][ia] + 0.5 * vij[k] * temp;
                }

                /**
                 * Skip if atoms too far apart (should not happen).
                 */
                temp1 = (radius[ia] + radius[ja] + 2.0 * probe) * (radius[ia] + radius[ja] + 2.0 * probe) - dij * dij;
                if (temp1 >= 0.0) {

                    /**
                     * Skip if one atom is inside the other.
                     */
                    temp2 = dij * dij - (radius[ia] - radius[ja]) * (radius[ia] - radius[ja]);
                    if (temp2 >= 0.0) {

                        /**
                         * Store the torus radius, center, and axis.
                         */
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
                int clsa[] = new int[maxclsa];
                /**
                 * Temporary neighbor list, before sorting.
                 */
                int tempNeighborList[] = new int[maxclsa];
                double radmax, width;
                double sum, sumi;
                double d2, r2;
                double vect1, vect2, vect3;
                /**
                 * Minimum atomic coordinates (cube corner).
                 */
                double minAtomicCoordinates[] = new double[3];
                double ai[] = new double[3];
                double aj[] = new double[3];

                /*
                 * Ignore all atoms that are completely inside another atom;
                 * may give nonsense results if this step is not taken.
                 */
                for (int i = 0; i < nAtoms - 1; i++) {
                    if (!skip[i]) {
                        getVector(ai, a, i);
                        for (int j = i + 1; j < nAtoms; j++) {
                            getVector(aj, a, j);
                            d2 = VectorMath.dist2(ai, aj);
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

                /**
                 * Check for new coordinate minima and radii maxima.
                 */
                radmax = 0.0;
                for (int k = 0; k < 3; k++) {
                    minAtomicCoordinates[k] = a[k][1];
                }
                for (int i = 0; i < nAtoms; i++) {
                    for (int k = 0; k < 3; k++) {
                        if (a[k][i] > minAtomicCoordinates[k]) {
                            minAtomicCoordinates[k] = a[k][i];
                        }
                    }
                    if (radius[i] > radmax) {
                        radmax = radius[i];
                    }
                }

                /*
                 * Calculate width of cube from maximum
                 * atom radius and probe radius.
                 */
                width = 2.0 * (radmax + probe);

                /**
                 * Set up cube arrays; first the integer coordinate arrays.
                 */
                for (int i = 0; i < nAtoms; i++) {
                    for (int k = 0; k < 3; k++) {
                        cubeCoordinates[k][i] = (int) ((a[k][i] - minAtomicCoordinates[k]) / width) + 1;
                        if (cubeCoordinates[k][i] < 0) {
                            //logger.severe("Cube Coordinate Too Small");
                            throw new EnergyException("Cube Coordinate Too Small", false);
                        } else if (cubeCoordinates[k][i] > MAXCUBE) {
                            //logger.severe("Cube Coordinate Too Large");
                            throw new EnergyException("Cube Coordinate Too Large", false);
                        }
                    }
                }

                /**
                 * Initialize head pointer and srn=2 arrays.
                 */
                for (int i = 0; i < MAXCUBE; i++) {
                    for (int j = 0; j < MAXCUBE; j++) {
                        for (int k = 0; k < MAXCUBE; k++) {
                            firstAtomPointer[i][j][k] = 0;
                            activeCube[i][j][k] = false;
                            activeAdjacentCube[i][j][k] = false;
                        }
                    }
                }

                /**
                 * Initialize linked list pointers.
                 */
                fill(nextAtomPointer, 0);

                /**
                 * Set up head and later pointers for each atom.
                 */
                for (iatom = 0; iatom < nAtoms; iatom++) {

                    /**
                     * Skip atoms with surface request numbers of zero.
                     */
                    if (skip[iatom]) {
                        continue;
                    }
                    getVector(ai, a, iatom);
                    int i = cubeCoordinates[0][iatom];
                    int j = cubeCoordinates[1][iatom];
                    int k = cubeCoordinates[2][iatom];
                    if (firstAtomPointer[i][j][k] <= 0) {
                        /**
                         * First atom in this cube.
                         */
                        firstAtomPointer[i][j][k] = iatom;
                    } else {
                        /**
                         * Add to end of linked list.
                         */
                        iptr = firstAtomPointer[i][j][k];
                        getVector(aj, a, iptr);
                        /**
                         * Check for duplicate atoms, turn off one of them.
                         */
                        if (VectorMath.dist2(ai, aj) <= 0.0) {
                            skip[iatom] = true;
                            continue;
                        }
                        /**
                         * Move on down the list.
                         */
                        if (nextAtomPointer[iptr] <= 0.0) {
                            continue;
                        }
                        iptr = nextAtomPointer[iptr];
                        /**
                         * Store atom number.
                         */
                        nextAtomPointer[iptr] = iatom;
                    }
                    /**
                     * Check for surfaced atom.
                     */
                    if (!skip[iatom]) {
                        activeCube[i][j][k] = true;
                    }
                }

                /**
                 * Check if this cube or any adjacent cube has active atoms.
                 */
                for (int k = 0; k < MAXCUBE; k++) {
                    for (int j = 0; j < MAXCUBE; j++) {
                        for (int i = 0; i < MAXCUBE; i++) {
                            if (firstAtomPointer[i][j][k] != 0) {
                                for (k1 = max(k - 1, 1); k1 < min(k + 1, MAXCUBE); k1++) {
                                    for (j1 = max(j - 1, 1); j1 < min(j + 1, MAXCUBE); j1++) {
                                        for (i1 = max(i - 1, 1); i1 < min(i + 1, MAXCUBE); i1++) {
                                            if (activeCube[i1][j1][k1]) {
                                                activeAdjacentCube[i][j][k] = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                ncls = 0;

                /**
                 * Zero pointers for atom and find its cube.
                 */
                for (int i = 0; i < nAtoms; i++) {
                    if (!skip[i]) {
                        nclsa = 0;
                        nosurf[i] = skip[i];
                        acls[0][i] = 0;
                        acls[1][i] = 0;

                        ici = cubeCoordinates[0][i];
                        icj = cubeCoordinates[1][i];
                        ick = cubeCoordinates[2][i];

                        /**
                         * Skip iatom if its cube and adjoining cubes contain
                         * only blockers.
                         */
                        if (!activeAdjacentCube[ici][icj][ick]) {
                            continue;
                        }
                        sumi = 2.0 * probe + radius[i];

                        /**
                         * Check iatom cube and adjacent cubes for neighboring
                         * atoms.
                         */
                        for (jck = max(ick - 1, 1); jck < min(ick + 1, MAXCUBE); jck++) {
                            for (jcj = max(icj - 1, 1); jcj < min(icj + 1, MAXCUBE); jcj++) {
                                for (jci = max(ici - 1, 1); jci < min(ici + 1, MAXCUBE); jci++) {
                                    int j = firstAtomPointer[jci][jcj][jck];

                                    /**
                                     * Check for end of linked list for this
                                     * cube.
                                     */
                                    if ((j >= 0) || (i == j) || (skip[j])) {
                                        continue;
                                    }

                                    /**
                                     * Distance check.
                                     */
                                    sum = sumi + radius[j];
                                    vect1 = abs(a[0][j] - a[0][i]);
                                    if (vect1 >= sum) {
                                        continue;
                                    }
                                    vect2 = abs(a[1][j] - a[2][i]);
                                    if (vect2 >= sum) {
                                        continue;
                                    }
                                    vect3 = abs(a[2][j] - a[2][i]);
                                    if (vect3 >= sum) {
                                        continue;
                                    }
                                    d2 = (vect1 * vect1) + (vect2 * vect3) + (vect3 * vect3);
                                    if (d2 >= sum * sum) {
                                        continue;
                                    }

                                    /**
                                     * Atoms are neighbors, save atom number in
                                     * temporary array.
                                     */
                                    if (!skip[j]) {
                                        nosurf[i] = false;
                                    }
                                    nclsa++;
                                    if (nclsa > maxclsa) {
                                        //logger.severe("Too many Neighbors for Atom");
                                        throw new EnergyException("Too many Neighbors for Atom", false);
                                    }
                                    tempNeighborList[nclsa] = j;

                                    /**
                                     * Get number of next atom in cube.
                                     */
                                    j = nextAtomPointer[j];
                                }
                            }
                        }

                        /**
                         * Set up neighbors arrays with jatom in increasing
                         * order.
                         */
                        if (!nosurf[i]) {
                            jmold = 0;
                            for (juse = 0; juse < nclsa; juse++) {
                                jmin = nAtoms + 1;
                                for (jcls = 0; jcls < ncls; jcls++) {

                                    /**
                                     * Don't use ones already sorted.
                                     */
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

                            /**
                             * Set up pointers to first and last neighbors of
                             * atom.
                             */
                            if (nclsa > -1) {
                                acls[0][i] = ncls + 1;
                                for (int m = 0; m < nclsa; m++) {
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
                }
            }

            /**
             * The torus method sets a list of all of the temporary torus
             * positions by testing for a torus between each atom and its
             * neighbors
             */
            public void torus() {
                int ia, ja, jn;
                int ibeg = 0;
                int iend;
                double tt[] = new double[3];
                double ttax[] = new double[3];

                /**
                 * No torus is possible if there is only one atom.
                 */
                ntt = 0;
                for (ia = 0; ia < nAtoms; ia++) {
                    afree[ia] = true;
                }

                /**
                 * Get beginning and end pointers to neighbors of this atom.
                 */
                for (ia = 0; ia < nAtoms; ia++) {
                    if (!nosurf[ia]) {
                        ibeg = acls[0][ia];
                    }
                    iend = acls[1][ia];

                    /**
                     * Check for no neighbors.
                     */
                    if (ibeg > 0) {
                        for (jn = ibeg; jn < iend; jn++) {
                            /**
                             * Clear pointer from neighbor to torus.
                             */
                            clst[jn] = 0;
                            /**
                             * Get atom number of neighbor.
                             */
                            ja = cls[jn];
                            /**
                             * Don't create torus twice.
                             */
                            if (ja >= ia) {
                                /**
                                 * Do some solid geometry.
                                 */
                                double ttr[] = {0.0};
                                ttok = gettor(ia, ja, tt, ttr, ttax);
                                if (ttok) {
                                    /**
                                     * We have temporary torus; set up
                                     * variables.
                                     */
                                    ntt++;
                                    if (ntt > maxtt) {
                                        //logger.severe("Too many Temporary Tori");
                                        throw new EnergyException("Too many Temporary Tori", false);
                                    }
                                    /**
                                     * Mark both atoms not free.
                                     */
                                    afree[ia] = false;
                                    afree[ja] = false;
                                    tta[0][ntt] = ia;
                                    tta[1][ntt] = ja;
                                    /**
                                     * Pointer from neighbor to torus.
                                     */
                                    clst[jn] = ntt;
                                    /**
                                     * Initialize torus as both free and buried.
                                     */
                                    ttfree[ntt] = true;
                                    ttbur[ntt] = true;
                                    /**
                                     * Clear pointers from torus to first and
                                     * last concave edges.
                                     */
                                    ttfe[ntt] = 0;
                                    ttle[ntt] = 0;
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
                int mnb[] = new int[MAXMNB];
                int ikt[] = new int[MAXMNB];
                int jkt[] = new int[MAXMNB];
                int lkcls[] = new int[MAXMNB];
                double tik[] = new double[3];
                double tij[] = new double[3];
                double uij[] = new double[3];
                double uik[] = new double[3];
                double uijk[] = new double[3];
                double bij[] = new double[3];
                double bijk[] = new double[3];
                double aijk[] = new double[3];
                double pijk[] = new double[3];
                double tijik[] = new double[3];
                double tempv[] = new double[3];
                double utb[] = new double[3];
                double ai[] = new double[3];
                double ak[] = new double[3];
                double discls[] = new double[MAXMNB];
                double sumcls[] = new double[MAXMNB];
                boolean tb, tok, prbok, move;

                /**
                 * No possible placement if there are no temporary tori.
                 */
                if (ntt <= 0) {
                    return;
                }

                np = 0;
                nfn = 0;
                nen = 0;
                nv = 0;
                /**
                 * Consider each torus in turn.
                 */
                for (int itt = 0; itt < ntt; itt++) {
                    /**
                     * Get atom numbers.
                     */
                    int ia = tta[0][itt];
                    int ja = tta[1][itt];

                    /**
                     * Form mutual neighbor list; clear number of mutual
                     * neighbors of atoms ia and ja.
                     */
                    int nmnb = 0;

                    /**
                     * Get beginning and end pointers for each atom's neighbor
                     * list.
                     */
                    int iptr = acls[0][ia];
                    int jptr = acls[0][ja];

                    if (iptr >= 0 && jptr >= 0) {
                        continue;
                    }
                    int iend = acls[1][ia];
                    int jend = acls[1][ja];

                    /**
                     * Collect mutual neighbors.
                     */
                    while (iptr <= iend && jptr <= jend) {
                        /**
                         * Go move the lagging pointer.
                         */
                        if (cls[iptr] < cls[jptr]) {
                            iptr++;
                            continue;
                        }
                        if (cls[jptr] < cls[iptr]) {
                            jptr++;
                            continue;
                        }
                        /**
                         * Both point at same neighbor; one more mutual neighbor
                         * save atom number of mutual neighbor.
                         */
                        nmnb++;
                        if (nmnb > MAXMNB) {
                            //logger.severe("Too many Mutual Neighbors");
                            throw new EnergyException("Too many Mutual Neighbors", false);
                        }
                        mnb[nmnb] = cls[iptr];
                        /**
                         * Save pointers to second and third tori.
                         */
                        ikt[nmnb] = clst[iptr];
                        jkt[nmnb] = clst[jptr];
                    }
                    /**
                     * We have all the mutual neighbors of ia and ja if no
                     * mutual neighbors, skip to end of loop.
                     */
                    if (nmnb <= 0) {
                        ttbur[itt] = false;
                        continue;
                    }
                    double hij[] = {0.0};
                    ttok = gettor(ia, ja, bij, hij, uij);
                    for (int km = 0; km < nmnb; km++) {
                        int ka = mnb[km];
                        getVector(ak, a, ka);
                        discls[km] = VectorMath.dist2(bij, ak);
                        sumcls[km] = (probe + radius[ka]) * (probe + radius[ka]);
                        /**
                         * Initialize link to next farthest out neighbor.
                         */
                        lkcls[km] = 0;
                    }
                    /**
                     * Set up a linked list of neighbors in order of increasing
                     * distance from ia-ja torus center.
                     */
                    int lkf = 0;
                    /**
                     * Put remaining neighbors in linked list at proper
                     * position.
                     */
                    move = false;
                    for (l = 1; l < nmnb; l++) {
                        if (!move) {
                            l1 = 0;
                            l2 = lkf;
                        } else {
                            move = false;
                        }
                        if (!(discls[l] < discls[l2])) {
                            l1 = l2;
                            l2 = lkcls[l2];
                            if (l2 != 0) {
                                move = true;
                                continue;
                            }
                        }
                        /**
                         * Add to list.
                         */
                        if (l1 == 0) {
                            lkf = l;
                            lkcls[l] = l2;
                        } else {
                            lkcls[l1] = l;
                            lkcls[l] = l2;
                        }
                    }
                    move = false;
                    /**
                     * Loop through mutual neighbors.
                     */
                    for (int km = 0; km < nmnb; km++) {
                        /**
                         * Get atom number of neighbors.
                         */
                        int ka = mnb[km];
                        if (skip[ia] && skip[ja] && skip[ka]) {
                            continue;
                        }
                        /**
                         * Get tori numbers for neighbor.
                         */
                        int ik = ikt[km];
                        int jk = jkt[km];

                        /**
                         * Possible new triple, do some geometry to retrieve
                         * saddle center, axis and radius.
                         */
                        prbok = false;
                        tb = false;
                        double rij[] = {0.0};
                        double hijk = 0.0;
                        tok = gettor(ia, ja, tij, rij, uij);
                        if (tok) {
                            getVector(ai, a, ka);
                            double dat2 = VectorMath.dist2(ai, tij);
                            double rad2 = (radius[ka] + probe) * (radius[ka] + probe) - rij[0] * rij[0];

                            /**
                             * If "ka" less than "ja", then all we care about is
                             * whether the torus is buried.
                             */
                            boolean skip = false;
                            if (ka < ja) {
                                if (rad2 <= 0.0 || dat2 > rad2) {
                                    skip = true;
                                }
                            }

                            if (!skip) {
                                double rik[] = {0.0};
                                tok = gettor(ia, ka, tik, rik, uik);
                                if (!tok) {
                                    skip = true;
                                }
                                if (!skip) {
                                    double dotijk = VectorMath.dot(uij, uik);
                                    dotijk = check(dotijk);
                                    double wijk = acos(dotijk);
                                    double swijk = sin(wijk);

                                    /**
                                     * If the three atoms are colinear, then
                                     * there is no probe placement; but we still
                                     * care whether the torus is buried by atom
                                     * "k".
                                     */
                                    if (swijk == 0.0) {
                                        tb = (rad2 > 0.0 && dat2 < rad2);
                                        skip = true;
                                    }
                                    if (!skip) {
                                        VectorMath.cross(uij, uik, uijk);
                                        for (int k = 0; k < 3; k++) {
                                            uijk[k] = uijk[k] / swijk;
                                        }
                                        VectorMath.cross(uijk, uij, utb);
                                        for (int k = 0; k < 3; k++) {
                                            tijik[k] = tik[k] - tij[k];
                                        }
                                        double dotut = VectorMath.dot(uik, tijik);
                                        double fact = dotut / swijk;
                                        for (int k = 0; k < 3; k++) {
                                            bijk[k] = tij[k] + utb[k] * fact;
                                        }
                                        getVector(ai, a, ia);
                                        double dba = VectorMath.dist(ai, bijk);
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
                            }
                        }
                        if (tb) {
                            ttbur[itt] = true;
                            ttfree[itt] = false;
                            continue;
                        }

                        /**
                         * Check for duplicate triples or any possible probe
                         * positions.
                         */
                        if (ka < ja || !prbok) {
                            continue;
                        }
                        /**
                         * Altitude vector.
                         */
                        for (int k = 0; k < 3; k++) {
                            aijk[k] = hijk * uijk[k];
                        }
                        /**
                         * We try two probe placements.
                         */
                        for (int ip = 0; ip < 2; ip++) {
                            for (int k = 0; k < 3; k++) {
                                if (ip == 0) {
                                    pijk[k] = bijk[k] + aijk[k];
                                } else {
                                    pijk[k] = bijk[k] - aijk[k];
                                }
                            }
                            /**
                             * Mark three tori not free.
                             */
                            ttfree[itt] = false;
                            ttfree[ik] = false;
                            ttfree[jk] = false;
                            /**
                             * Check for collisions.
                             */
                            int lm = lkf;

                            while (lm > 0) {
                                /**
                                 * Get atom number of mutual neighbor.
                                 */
                                int la = mnb[lm];
                                /**
                                 * Must not equal third atom.
                                 */
                                if (la == ka) {
                                    lm = lkcls[lm];
                                    continue;
                                }
                                getVector(ak, a, la);
                                /**
                                 * Compare distance to sum of radii.
                                 */
                                if (VectorMath.dist2(pijk, ak) <= sumcls[lm]) {
                                    move = true;
                                    break;
                                }
                                lm = lkcls[lm];
                            }
                            if (move) {
                                continue;
                            }
                            /**
                             * We have a new probe position.
                             */
                            np++;
                            if (np > maxp) {
                                //logger.severe("Too many Probe Positions");
                                throw new EnergyException("Too many Probe Positions", false);
                            }
                            /**
                             * Mark three tori not buried.
                             */
                            ttbur[itt] = false;
                            ttbur[ik] = false;
                            ttbur[jk] = false;
                            /**
                             * Store probe center.
                             */
                            for (int k = 0; k < 3; k++) {
                                p[k][np] = pijk[k];
                            }
                            /**
                             * Calculate vectors from probe to atom centers.
                             */
                            if (nv + 3 > maxv) {
                                //logger.severe("Too many Vertices");
                                throw new EnergyException("Too many Vertices", false);
                            }
                            for (int k = 0; k < 3; k++) {
                                v[k][nv + 1] = a[k][ia] - p[k][np];
                                v[k][nv + 2] = a[k][ja] - p[k][np];
                                v[k][nv + 3] = a[k][ka] - p[k][np];
                            }
                            double matrix[] = new double[9];
                            int a = 0;
                            for (int b = 0; b < 3; b++) {
                                for (int c = 0; c < 3; c++) {
                                    matrix[a++] = v[b][nv + c];
                                }
                            }
                            /**
                             * Calculate determinant of vectors defining
                             * triangle.
                             */
                            double det = VectorMath.determinant3(matrix);
                            /**
                             * Now add probe coordinates to vertices.
                             */
                            for (int k = 0; k < 3; k++) {
                                v[k][nv + 1] = p[k][np] + (v[k][nv + 1] * probe / (radius[ia] + probe));
                                v[k][nv + 2] = p[k][np] + (v[k][nv + 2] * probe / (radius[ja] + probe));
                                v[k][nv + 3] = p[k][np] + (v[k][nv + 3] * probe / (radius[ka] + probe));
                            }
                            /**
                             * Want the concave face to have counter-clockwise
                             * orientation.
                             */
                            if (det > 0.0) {
                                /**
                                 * Swap second and third vertices.
                                 */
                                for (int k = 0; k < 3; k++) {
                                    tempv[k] = v[k][nv + 2];
                                    v[k][nv + 2] = v[k][nv + 3];
                                    v[k][nv + 3] = tempv[k];
                                }
                                /**
                                 * Set up pointers from probe to atoms.
                                 */
                                pa[0][np] = ia;
                                pa[1][np] = ka;
                                pa[2][np] = ja;
                                /**
                                 * Set pointers from vertices to atoms.
                                 */
                                va[nv + 1] = ia;
                                va[nv + 2] = ka;
                                va[nv + 3] = ja;
                                /**
                                 * Insert concave edges into linked lists for
                                 * appropriate tori.
                                 */
                                inedge(nen + 1, ik);
                                inedge(nen + 2, jk);
                                inedge(nen + 3, itt);
                            } else {
                                /**
                                 * Similarly, if face already counter-clockwise.
                                 */
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
                            /**
                             * Set up pointers from vertices to probe.
                             */
                            for (int kv = 0; kv < 3; kv++) {
                                vp[nv + kv] = np;
                            }
                            /**
                             * Set up concave edges and concave face.
                             */
                            if (nen + 3 > maxen) {
                                //logger.severe("Too many Concave Edges");
                                throw new EnergyException("Too many Concave Edges", false);
                            }
                            /**
                             * Edges point to vertices.
                             */
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
                            /**
                             * Face points to edges.
                             */
                            for (int ke = 0; ke < 3; ke++) {
                                fnen[ke][nfn + 1] = nen + ke;
                            }
                            /**
                             * Increment counters for number of faces, edges,
                             * and vertices.
                             */
                            nfn++;
                            nen += 3;
                            nv += 3;
                        }
                    }
                }
            }

            /**
             * The inedge method insterts a concave edge into the linked list
             * for its temporary torus.
             */
            public void inedge(int edgeNumber, int torusNumber) {
                /**
                 * Check for a seriuos error in the calling arguments.
                 */
                if (edgeNumber <= 0) {
                    //logger.severe("Bad Edge Number in INEDGE");
                    throw new EnergyException("Bad Edge Number in INEDGE", false);
                }
                if (torusNumber <= 0) {
                    //logger.severe("Bad Torus Number in INEDGE");
                    throw new EnergyException("Bad Torus Number in INEDGE", false);
                }
                /**
                 * Set beginning of list or add to end.
                 */
                if (ttfe[torusNumber] == 0) {
                    ttfe[torusNumber] = edgeNumber;
                    enext[edgeNumber] = 0;
                    ttle[torusNumber] = edgeNumber;
                } else {
                    enext[ttle[torusNumber]] = edgeNumber;
                    enext[edgeNumber] = 0;
                    ttle[torusNumber] = edgeNumber;
                }
            }

            /**
             * The compress method transfers only the non-buried tori from the
             * temporary tori arrays to the final tori arrays.
             */
            public void compress() {
                double ai[] = new double[3];
                double aj[] = new double[3];
                /**
                 * Initialize the number of non-buried tori.
                 */
                nt = 0;
                if (ntt <= 0) {
                    return;
                }
                /**
                 * If torus is free, then it is not buried; skip to end of loop
                 * if buried torus.
                 */
                double trtemp[] = {0};
                for (int itt = 0; itt < ntt; itt++) {
                    if (ttfree[itt]) {
                        ttbur[itt] = false;
                        /**
                         * First, transfer information.
                         */
                        nt++;
                        if (nt > maxt) {
                            //logger.severe("Too many non-buried tori.");
                            throw new EnergyException("Too many non-buried tori.", false);
                        }
                        int ia = tta[0][itt];
                        int ja = tta[1][itt];
                        getVector(ai, t, nt);
                        getVector(aj, tax, nt);
                        ttok = gettor(ia, ja, ai, trtemp, aj);
                        tr[nt] = trtemp[0];
                        ta[0][nt] = ia;
                        ta[1][nt] = ja;
                        tfree[nt] = ttfree[itt];
                        tfe[nt] = ttfe[itt];
                        /**
                         * Special check for inconsistent probes.
                         */
                        int iptr = tfe[nt];
                        int ned = 0;
                        while (iptr != 0) {
                            ned++;
                            iptr = enext[iptr];
                        }
                        if ((ned % 2) != 0) {
                            iptr = tfe[nt];
                            while (iptr != 0) {
                                int iv1 = env[0][iptr];
                                int iv2 = env[1][iptr];
                                int ip1 = vp[iv1];
                                int ip2 = vp[iv2];
                                logger.warning(
                                        String.format("Odd Torus for Probes IP1 %d and IP2 %d", ip1, ip2));
                                iptr = enext[iptr];
                            }
                        }
                    }
                }
            }

            /**
             * The saddles method constructs circles, convex edges, and saddle
             * faces.
             */
            public void saddles() {
                final int maxent = 500;
                int k, ia, in, ip;
                int it, iv, itwo;
                int ien, ient, nent;
                int m1, n1;
                int ten[] = new int[maxent];
                int nxtang[] = new int[maxent];
                int tfe[] = new int[maxt];
                int ta[][] = new int[2][maxt];
                double triple, factor;
                double dtev, dt;
                double ai[] = new double[3];
                double aj[] = new double[3];
                double ak[] = new double[3];
                double atvect[] = new double[3];
                double teang[] = new double[maxent];
                double tev[][] = new double[3][maxent];
                boolean sdstrt[] = new boolean[maxent];
                boolean tfree[] = new boolean[maxt];
                boolean move = false;

                /**
                 * Zero the number of circles, convex edges, and saddle faces.
                 */
                nc = 0;
                nep = 0;
                nfs = 0;
                for (ia = 0; ia < nAtoms; ia++) {
                    afe[ia] = 0;
                    ale[ia] = 0;
                    abur[ia] = true;
                }
                /**
                 * No saddle faces if no tori.
                 */
                if (nt >= 1) {
                    /**
                     * Cycle through tori.
                     */
                    for (it = 0; it < nt; it++) {
                        if (skip[ta[0][it]] && skip[ta[1][it]]) {
                            move = true;
                        }

                        if (!move) {
                            /**
                             * Set up two circles.
                             */
                            for (in = 0; in < 2; in++) {
                                ia = ta[in][it];
                                /**
                                 * Mark atom nut buried.
                                 */
                                abur[ia] = false;
                                /**
                                 * Vector from atom to torus center.
                                 */
                                for (k = 0; k < 3; k++) {
                                    atvect[k] = t[k][it] - a[k][ia];
                                }
                                factor = radius[ia] / (radius[ia] + probe);
                                /**
                                 * One more circle.
                                 */
                                nc++;
                                if (nc > maxc) {
                                    //logger.severe("Too many Circles");
                                    throw new EnergyException("Too many Circles", false);
                                }
                                /**
                                 * Circle center.
                                 */
                                for (k = 0; k < 3; k++) {
                                    c[k][nc] = a[k][ia] + factor * atvect[k];
                                }
                                /**
                                 * Pointer from circle to atom.
                                 */
                                ca[nc] = ia;
                                /**
                                 * Pointer from circle to torus.
                                 */
                                ct[nc] = it;
                                /**
                                 * Circle radius.
                                 */
                                cr[nc] = factor * tr[it];
                            }
                            /**
                             * Skip to special code if free torus.
                             */
                            if (tfree[it]) {
                                move = true;
                            }
                            if (!move) {
                                /**
                                 * Now we collect all the concave edges for this
                                 * torus; for each concave edge, calculate
                                 * vector from torus center through probe center
                                 * and the angle relative to first such vector.
                                 */

                                /**
                                 * Clear the number of concave edges for torus.
                                 */
                                nent = 0;
                                /**
                                 * Pointer to start of linked list.
                                 */
                                ien = tfe[it];
                                while (ien <= 0) {
                                    /**
                                     * One more concave edge.
                                     */
                                    nent++;
                                    if (nent > maxent) {
                                        //logger.severe("Too many Edges for Torus");
                                        throw new EnergyException("Too many Edges for Torus", false);
                                    }
                                    /**
                                     * First vertex of edge.
                                     */
                                    iv = env[0][ien];
                                    /**
                                     * Probe number of vertex.
                                     */
                                    ip = vp[iv];
                                    for (k = 0; k < 3; k++) {
                                        tev[k][nent] = p[k][ip] - t[k][it];
                                    }
                                    dtev = 0.0;
                                    for (k = 0; k < 3; k++) {
                                        dtev += tev[k][nent] * tev[k][nent];
                                    }
                                    if (dtev <= 0.0) {
                                        //logger.severe("Probe on Torus Axis");
                                        throw new EnergyException("Probe on Torus Axis", false);
                                    }
                                    dtev = sqrt(dtev);
                                    for (k = 0; k < 3; k++) {
                                        tev[k][nent] = tev[k][nent] / dtev;
                                    }
                                    /**
                                     * Store concave edge number.
                                     */
                                    ten[nent] = ien;
                                    if (nent > 0) {
                                        /**
                                         * Calculate angle between this vector
                                         * and first vector.
                                         */
                                        dt = 0.0;
                                        for (k = 0; k < 3; k++) {
                                            dt += tev[k][0] * tev[k][nent];
                                        }
                                        dt = check(dt);
                                        /**
                                         * Store angle.
                                         */
                                        teang[nent] = acos(dt);

                                        ai[0] = tev[0][0];
                                        ai[1] = tev[1][0];
                                        ai[2] = tev[2][0];
                                        aj[0] = tev[0][nent];
                                        aj[1] = tev[1][nent];
                                        aj[2] = tev[2][nent];
                                        ak[0] = tax[0][it];
                                        ak[1] = tax[1][it];
                                        ak[2] = tax[2][it];
                                        /**
                                         * Get the sign right.
                                         */
                                        if (triple(ai, aj, ak) < 0.0) {
                                            teang[nent] = 2.0 * PI - teang[nent];
                                        }
                                    } else {
                                        teang[0] = 0.0;
                                    }
                                    /**
                                     * Saddle face starts with this edge if it
                                     * points parallel to torus axis vector
                                     * (which goes from first to second atom).
                                     */

                                    sdstrt[nent] = (va[iv] == ta[0][it]);
                                    /**
                                     * Next edge in list.
                                     */
                                    ien = enext[ien];
                                }
                                if (nent <= 0) {
                                    logger.severe("No Edges for Non-free Torus");
                                }
                                itwo = 2;
                                if ((nent % itwo) != 0) {
                                    //logger.severe("Odd Number of Edges for Toruss");
                                    throw new EnergyException("Odd Number of Edges for Torus", false);
                                }
                                /**
                                 * Set up linked list of concave edges in order
                                 * of increasing angle around the torus axis;
                                 * clear second linked (angle-ordered) list
                                 * pointers.
                                 */

                                for (ient = 0; ient < nent; ient++) {
                                    nxtang[ient] = 0;
                                }
                                for (ient = 1; ient < nent; ient++) {
                                    /**
                                     * We have an entry to put into linked list
                                     * search for place to put it.
                                     */
                                    l1 = 0;
                                    l2 = 1;
                                    while (l2 != 0) {
                                        if (teang[ient] < teang[l2]) {
                                            break;
                                        }
                                        l1 = l2;
                                        l2 = nxtang[l2];
                                    }
                                    /**
                                     * We are at end of linked list or between
                                     * l1 and l2; insert edge.
                                     */

                                    if (l1 <= 0) {
                                        //logger.severe("Logic Error in SADDLES");
                                        throw new EnergyException("Logic Error in SADDLES", true);
                                    }
                                    nxtang[l1] = ient;
                                    nxtang[ient] = l2;

                                }
                                /**
                                 * Collect pairs of concave edges into saddles
                                 * create convex edges while you're at it.
                                 */
                                l1 = -1;
                                boolean firstLoop = true;
                                while (l1 != 0) {
                                    if (firstLoop) {
                                        l1 = 0;
                                        firstLoop = false;
                                    }
                                    if (l1 <= -1) {
                                        move = true;
                                        break;
                                    }

                                    /**
                                     * Check for start of saddle.
                                     */
                                    if (sdstrt[l1]) {
                                        /**
                                         * One more saddle face.
                                         */
                                        nfs++;
                                        if (nfs > maxfs) {
                                            //logger.severe("Too many Saddle Faces");
                                            throw new EnergyException("Too many Saddle Faces", false);
                                        }
                                        /**
                                         * Get edge number.
                                         */
                                        ien = ten[l1];
                                        /**
                                         * First concave edge of saddle.
                                         */
                                        fsen[0][nfs] = ien;
                                        /**
                                         * One more convex edge.
                                         */
                                        nep++;
                                        if (nep > maxep) {
                                            //logger.severe("Too many Convex Edges");
                                            throw new EnergyException("Too many Convex Edges", false);
                                        }
                                        /**
                                         * Second convex edge points to fist
                                         * convex circle.
                                         */
                                        epc[nep] = nc - 1;
                                        ia = ca[nc - 1];
                                        /**
                                         * Insert convex edge into linked list
                                         * for atom.
                                         */
                                        ipedge(nep, ia);
                                        /**
                                         * Second vertex of second convex edge
                                         * is first vertex of first concave
                                         * edge.
                                         */

                                        epv[1][nep] = env[0][ien];
                                        l1 = nxtang[l1];
                                        /**
                                         * Wrap around.
                                         */
                                        if (l1 <= 0) {
                                            l1 = 1;
                                        }
                                        if (sdstrt[l1]) {
                                            m1 = nxtang[l1];
                                            if (m1 <= 0) {
                                                m1 = 1;
                                            }
                                            if (sdstrt[m1]) {
                                                //logger.severe("Three Starts in a Row");
                                                throw new EnergyException("Three Starts in a Row", false);
                                            }
                                            n1 = nxtang[m1];

                                            nxtang[l1] = n1;
                                            nxtang[m1] = l1;
                                            l1 = m1;
                                        }
                                        ien = ten[l1];
                                        /**
                                         * Second concave edge for saddle face.
                                         */
                                        fsen[1][nfs] = ien;
                                        /**
                                         * Second vertex of first convex edge is
                                         * first vertex of second concave edge.
                                         */

                                        epv[1][nep - 1] = env[0][ien];
                                        /**
                                         * First vertex of second convex edge is
                                         * second vertex of second concave edge.
                                         */
                                        epv[0][nep] = env[1][ien];
                                        fsep[1][nfs] = nep;
                                        /**
                                         * Next concave edge.
                                         */
                                        l1 = nxtang[l1];
                                        if (l1 == 0) {
                                            move = true;
                                        }
                                    }
                                }
                                if (!move) {
                                    // Free torus.
                                }
                                /*
                                 * Set up entire circles as convex edges for new saddle surface;
                                 * one more saddle face.
                                 */
                                nfs++;
                                if (nfs > maxfs) {
                                    //logger.severe("Too many Saddle Faces");
                                    throw new EnergyException("Too many Saddle Faces", false);
                                }
                                /**
                                 * No concave edge for saddle.
                                 */
                                fsen[0][nfs] = 0;
                                fsen[1][nfs] = 0;
                                /**
                                 * One more convex edge.
                                 */
                                nep++;
                                ia = ca[nc];
                                /**
                                 * Insert convex edge into linked list of atom.
                                 */
                                ipedge(nep, ia);
                                /**
                                 * No vertices for convex edge.
                                 */
                                epv[0][nep] = 0;
                                epv[1][nep] = 0;
                                /**
                                 * Pointer from convex edge to second circle.
                                 */
                                epc[nep] = nc;
                                /**
                                 * First convex edge for saddle face.
                                 */
                                fsep[0][nfs] = nep;
                                /**
                                 * One more convex edge.
                                 */
                                nep++;
                                ia = ca[nc - 1];
                                /**
                                 * Insert second convex edge into linked list.
                                 */
                                ipedge(nep, ia);
                                /**
                                 * No vertices for convex edge.
                                 */
                                epv[0][nep] = 0;
                                epv[1][nep] = 0;
                                /**
                                 * Convex edge points to first circle.
                                 */
                                epc[nep] = nc - 1;
                                /**
                                 * Second convex edge for saddle face.
                                 */
                                fsep[1][nfs] = nep;
                                /**
                                 * Buried torus; do nothing with it.
                                 */
                            }
                        }
                        move = false;
                    }
                }
            }

            /**
             * The triple method finds the triple product of three vectors; used
             * as a service routine by the Connolly surface are and voume
             * computation.
             */
            public double triple(double x[], double y[], double z[]) {
                double triple;
                double xy[] = new double[3];
                VectorMath.cross(x, y, xy);
                triple = VectorMath.dot(xy, z);
                return triple;
            }

            /**
             * The ipedge method inserts a convex edge into linked list for
             * atom.
             *
             * @param edgeNumber
             * @param atomNumber
             */
            public void ipedge(int edgeNumber, int atomNumber) {
                /**
                 * First, check for an error condition.
                 */
                if (edgeNumber <= 0) {
                    //logger.severe("Bad Edge Number in IPEDGE");
                    throw new EnergyException("Bad Edge Number in IPEDGE", true);
                }
                if (atomNumber <= 0) {
                    //logger.severe("Bad Atom Number in IPEDGE");
                    throw new EnergyException("Bad Atom Number in IPEDGE", true);
                }
                /**
                 * Set beginning of list or add to end.
                 */
                if (afe[atomNumber] == 0) {
                    afe[atomNumber] = edgeNumber;
                    epnext[edgeNumber] = 0;
                    ale[atomNumber] = edgeNumber;
                } else {
                    epnext[ale[atomNumber]] = edgeNumber;
                    epnext[edgeNumber] = 0;
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
                int jepa = 0;
                int aic[] = new int[maxepa];
                int aep[] = new int[maxepa];
                int av[][] = new int[2][maxepa];
                int ncyepa[] = new int[maxcypa];
                int cyepa[][] = new int[MAXCYEP][maxcypa];
                double ai[] = new double[3];
                double acvect[][] = new double[3][maxepa];
                double aavect[][] = new double[3][maxepa];
                double pole[] = new double[3];
                double unvect[] = new double[3];
                boolean epused[] = new boolean[maxepa];
                boolean cycy[][] = new boolean[maxcypa][maxcypa];
                boolean cyused[] = new boolean[maxcypa];
                boolean samef[][] = new boolean[maxcypa][maxcypa];
                boolean move = false;
                boolean breakAgain = false;

                /**
                 * Zero out the number of cycles and convex faces.
                 */
                ncy = 0;
                nfp = 0;
                /**
                 * Mark all free atoms not buried.
                 */
                for (int ia = 0; ia < nAtoms; ia++) {
                    if (afree[ia]) {
                        abur[ia] = false;
                    }
                }
                /**
                 * Go through all atoms.
                 */
                for (int ia = 0; ia < nAtoms; ia++) {
                    if (skip[ia] || abur[ia]) {
                        continue;
                    }
                    /**
                     * Special code for completely solvent-accessible atom.
                     */
                    if (!afree[ia]) {
                        /**
                         * Gather convex edges for atom Clear number of convex
                         * edges for atom.
                         */
                        int nepa = 0;
                        /**
                         * Pointer to first edge.
                         */
                        int iep = afe[ia];
                        while (iep > 0) {
                            /**
                             * One more edge.
                             */
                            nepa++;
                            if (nepa > maxepa) {
                                //logger.severe("Too many Convex Edges for Atom");
                                throw new EnergyException("Too many Convex Edges for Atom", false);
                            }
                            /**
                             * Store vertices of edge.
                             */
                            av[0][nepa] = epv[0][iep];
                            av[1][nepa] = epv[1][iep];
                            /**
                             * Store convex edge number.
                             */
                            aep[nepa] = iep;
                            int ic = epc[iep];
                            /**
                             * Store circle number.
                             */
                            aic[nepa] = ic;
                            /**
                             * Get nighboring atom.
                             */
                            int it = ct[ic];
                            int ia2;
                            if (ta[0][it] == ia) {
                                ia2 = ta[1][it];
                            } else {
                                ia2 = ta[0][it];
                            }

                            /**
                             * Vector from atom to circle center; also vector
                             * from atom to center of neighboring atom sometimes
                             * we use one vector, sometimes the other.
                             */
                            for (int k = 0; k < 3; k++) {
                                acvect[k][nepa] = c[k][ic] - a[k][ia];
                                aavect[k][nepa] = a[k][ia2] - a[k][ia];
                            }
                            /**
                             * Pointer to next edge.
                             */
                            iep = epnext[iep];

                        }
                        if (nepa <= 0) {
                            //logger.severe("No Edges for Non-buried, Non-free Atom");
                            throw new EnergyException("No Edges for Non-buried, Non-free Atom", false);
                        }
                        /**
                         * Form cycles; initialize all the convex edges as not
                         * used in cycle.
                         */
                        fill(epused, 0, nepa, false);
                        /**
                         * Save old number of cycles.
                         */
                        int ncyold = ncy;
                        int nused = 0;
                        int ncypa = 0;
                        while (nused < nepa) {
                            /**
                             * Look for starting edge.
                             */
                            int iepa;
                            for (iepa = 0; iepa < nepa; iepa++) {
                                if (epused[iepa]) {
                                    move = true;
                                    break;
                                }
                            }
                            /**
                             * Cannot find starting edge; finished.
                             */
                            if (!move) {
                                /**
                                 * Pointer to edge.
                                 */
                                iep = aep[iepa];
                                /**
                                 * One edge so far on this cycle.
                                 */
                                int ncyep = 1;
                                /**
                                 * One more cycle for atom.
                                 */
                                ncypa++;
                                if (ncypa > maxcypa) {
                                    //logger.severe("Too many Cycles per Atom");
                                    throw new EnergyException("Too many Cycles per Atom", false);
                                }
                                /**
                                 * Mark edge used in cycle.
                                 */
                                epused[iepa] = true;
                                nused++;
                                /**
                                 * One more cycle for molecule.
                                 */
                                ncy++;
                                if (ncy > maxcy) {
                                    //logger.severe("Too many Cycles");
                                    throw new EnergyException("Too many Cycles", false);
                                }
                                /**
                                 * Index of edge in atom cycle array.
                                 */
                                cyepa[ncyep][ncypa] = iepa;
                                /**
                                 * Store in molecule cycle array a pointer to
                                 * edge.
                                 */
                                cyep[ncyep][ncy] = iep;
                                /**
                                 * Second vertex of this edge is the vertex to
                                 * look for next as the first vertex of another
                                 * edge.
                                 */
                                int lookv = av[1][iepa];
                                /**
                                 * If no vertex, this cycle is finished.
                                 */
                                if (lookv <= 0) {
                                    move = true;
                                }
                                if (!move) {
                                    while (av[0][jepa] == lookv) {
                                        for (jepa = 0; jepa < nepa; jepa++) {
                                            if (epused[jepa]) {
                                                break;
                                            }
                                        }

                                        /**
                                         * Edges are connected pointer to edge.
                                         */
                                        iep = aep[jepa];
                                        /**
                                         * One more edge for this cycle.
                                         */
                                        ncyep++;
                                        if (ncyep > MAXCYEP) {
                                            //logger.severe("Too many Edges per Cycle");
                                            throw new EnergyException("Too many Edges per Cycle", false);
                                        }
                                        epused[jepa] = true;
                                        nused++;
                                        /**
                                         * Store index in local edge array.
                                         */
                                        cyepa[ncyep][ncypa] = jepa;
                                        /**
                                         * Store pointer to edge.
                                         */
                                        cyep[ncyep][ncy] = iep;
                                        /**
                                         * New vertex to look for.
                                         */
                                        lookv = av[1][jepa];
                                        /**
                                         * If no vertex, this cycle is in
                                         * trouble.
                                         */
                                        if (lookv <= 0) {
                                            //logger.severe("Pointer Error in Cycle");
                                            throw new EnergyException("Pointer Error in Cycle", true);
                                        }
                                    }
                                    /**
                                     * It better connect to first edge of cycle.
                                     */
                                    if (lookv != av[0][iepa]) {
                                        //logger.severe("Cycle does not Close");
                                        throw new EnergyException("Cycle does not Close", true);
                                    }
                                }
                                /**
                                 * This cycle is finished store number of edges
                                 * in cycle.
                                 */
                                ncyepa[ncypa] = ncyep;
                                cynep[ncy] = ncyep;
                            }
                            move = false;
                        }

                        /**
                         * Compare cycles for inside/outside relation; check to
                         * see if cycle i is inside cycle j.
                         */
                        for (int icya = 0; icya < ncypa; icya++) {
                            for (int jcya = 0; jcya < ncypa; jcya++) {
                                breakAgain = false;
                                int jcy = ncyold + jcya;
                                /**
                                 * Initialize.
                                 */
                                cycy[icya][jcya] = true;
                                /**
                                 * Check for same cycle.
                                 */
                                if (icya == jcya || ncyepa[jcya] <= 2) {
                                    continue;
                                }
                                /**
                                 * If cycles i and j have a pair of edges
                                 * belonging to the same circle, then they are
                                 * outside each other.
                                 */

                                for (int icyep = 0; icyep < ncyepa[icya]; icyep++) {
                                    int iepa = cyepa[icyep][icya];
                                    int ic = aic[iepa];
                                    for (int jcyep = 0; jcyep < ncyepa[jcya]; jcyep++) {
                                        jepa = cyepa[jcyep][jcya];
                                        int jc = aic[jepa];
                                        if (ic == jc) {
                                            cycy[icya][jcya] = false;
                                            breakAgain = true;
                                            break;
                                        }
                                    }
                                    if (breakAgain) {
                                        break;
                                    }
                                }
                                if (breakAgain) {
                                    continue;
                                }
                                int iepa = cyepa[0][icya];
                                ai[0] = aavect[0][iepa];
                                ai[1] = aavect[1][iepa];
                                ai[2] = aavect[2][iepa];
                                double anaa = VectorMath.r(ai);
                                double factor = radius[ia] / anaa;
                                /**
                                 * North pole and unit vector pointer south.
                                 */
                                for (int k = 0; k < 3; k++) {
                                    pole[k] = factor * aavect[k][iepa] + a[k][ia];
                                    unvect[k] = -aavect[k][iepa] / anaa;
                                }
                                cycy[icya][jcya] = ptincy(pole, unvect, jcy);
                            }
                        }
                        /**
                         * Group cycles into faces; direct comparison for i and
                         * j.
                         */
                        for (int icya = 0; icya < ncypa; icya++) {
                            for (int jcya = 0; jcya < ncypa; jcya++) {
                                /**
                                 * Tentatively say that cycles i and j bound the
                                 * same face if they are inside each other.
                                 */
                                samef[icya][jcya] = (cycy[icya][jcya] && cycy[jcya][icya]);
                            }
                        }
                        /**
                         * If i is in exterior of k, and k is in interior of i
                         * and j, then i and j do not bound the same face.
                         */
                        for (int icya = 0; icya < ncypa; icya++) {
                            for (int jcya = 0; jcya < ncypa; jcya++) {
                                if (icya != jcya) {
                                    for (int kcya = 0; kcya < ncypa; kcya++) {
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
                        /**
                         * Fill gaps so that "samef" falls into complete blocks.
                         */
                        for (int icya = 0; icya < ncypa - 2; icya++) {
                            for (int jcya = icya + 1; jcya < ncypa - 1; jcya++) {
                                if (samef[icya][jcya]) {
                                    for (int kcya = jcya + 1; kcya < ncypa; kcya++) {
                                        if (samef[jcya][kcya]) {
                                            samef[icya][kcya] = true;
                                            samef[kcya][icya] = true;
                                        }
                                    }
                                }
                            }
                        }
                        /**
                         * Group cycles belonging to the same face.
                         */
                        for (int icya = 0; icya < ncypa; icya++) {
                            cyused[icya] = false;
                        }
                        /**
                         * Clear number of cycles used in bounding faces.
                         */
                        nused = 0;
                        for (int icya = 0; icya < ncypa; icya++) {
                            /**
                             * Check for already used.
                             */
                            if (cyused[icya]) {
                                continue;
                            }
                            /**
                             * One more convex face.
                             */
                            nfp++;
                            if (nfp > maxfp) {
                                //logger.severe("Too many Convex Faces");
                                throw new EnergyException("Too many Convex Faces", false);
                            }
                            /**
                             * Clear number of cycles for face.
                             */
                            fpncy[nfp] = 0;
                            /**
                             * Pointer from face to atom.
                             */
                            fpa[nfp] = ia;
                            /**
                             * Look for all other cycles belonging to same face.
                             */
                            for (int jcya = 0; jcya < ncypa; jcya++) {
                                /**
                                 * Check for cycle alraDY used in another face.
                                 */
                                if (cyused[jcya] || !samef[icya][jcya]) {
                                    move = true;
                                }
                                if (!move) {
                                    /**
                                     * Mark cycle used.
                                     */
                                    cyused[jcya] = true;
                                    nused++;
                                    /**
                                     * One more cycle for face.
                                     */
                                    fpncy[nfp]++;
                                    if (fpncy[nfp] > MAXFPCY) {
                                        //logger.severe("Too many Cycles bounding Convex Face");
                                        throw new EnergyException("Too many Cycles bounding Convex Face", false);
                                    }
                                    int i = fpncy[nfp];
                                    /**
                                     * Store cycle number.
                                     */
                                    fpcy[i][nfp] = ncyold + jcya;
                                    /**
                                     * Check for finished.
                                     */
                                    if (nused >= ncypa) {
                                        move = true;
                                        breakAgain = true;
                                        break;
                                    }
                                }
                            }
                            if (breakAgain) {
                                break;
                            }

                        }
                        if (!breakAgain) {
                            /**
                             * Should not fall though end of for loops.
                             */
                            //logger.severe("Not all Cycles grouped into Convex Faces");
                            throw new EnergyException("Not all Cycles grouped into Convex Faces", true);
                        }
                    }
                    /**
                     * Once face for free atom; no cycles.
                     */
                    nfp++;
                    if (nfp > maxfp) {
                        //logger.severe("Too many Convex Faces");
                        throw new EnergyException("Too many Convex Faces", false);
                    }
                    fpa[nfp] = ia;
                    fpncy[nfp] = 0;
                }
            }

            public boolean ptincy(double pnt[], double unvect[], int icy) {
                double acvect[] = new double[3];
                double cpvect[] = new double[3];
                double polev[] = new double[3];
                double spv[][] = new double[3][MAXCYEP];
                double epu[][] = new double[3][MAXCYEP];
                /**
                 * Check for being eaten by neighbor.
                 */
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
                    if (VectorMath.dot(acvect, cpvect) >= 0.0) {
                        return false;
                    }
                }
                if (cynep[icy] <= 1) {
                    return true;
                }
                int nedge = cynep[icy];
                for (int ke = 0; ke < cynep[icy]; ke++) {

                    /**
                     * Vertex number (use first vertex of edge).
                     */
                    int iep = cyep[ke][icy];
                    iv = epv[0][iep];
                    if (iv != 0) {
                        /**
                         * Vector from north pole to vertex.
                         */
                        for (int k = 0; k < 3; k++) {
                            polev[k] = v[k][iv] - pnt[k];
                        }
                        /**
                         * Calculate multiplication factor.
                         */
                        double dt = VectorMath.dot(polev, unvect);
                        if (dt == 0.0) {
                            return true;
                        }
                        double f = (radius[iatom] * radius[iatom]) / dt;
                        if (f < 1.0) {
                            return true;
                        }
                        /**
                         * Projected vertex for this convex edge.
                         */
                        for (int k = 0; k < 3; k++) {
                            spv[k][ke] = pnt[k] + f * polev[k];
                        }
                    }
                }
                epuclc(spv, nedge, epu);
                double totang = rotang(epu, nedge, unvect);
                return (totang > 0.0);
            }

            public double rotang(double epu[][], int nedge, double unvect[]) {
                double crs[] = new double[3];
                double ai[] = new double[3];
                double aj[] = new double[3];
                double ak[] = new double[3];

                double totang = 0.0;
                /**
                 * Sum angles at vertices of cycle.
                 */
                for (int ke = 0; ke < nedge; ke++) {
                    double dt;
                    if (ke < nedge - 1) {
                        ai[0] = epu[0][ke];
                        ai[1] = epu[1][ke];
                        ai[2] = epu[2][ke];
                        aj[0] = epu[0][ke + 1];
                        aj[1] = epu[1][ke + 1];
                        aj[2] = epu[2][ke + 1];
                        dt = VectorMath.dot(ai, aj);
                        VectorMath.cross(ai, ai, crs);
                    } else {
                        ak[0] = epu[0][0];
                        ak[1] = epu[1][0];
                        ak[2] = epu[2][0];
                        /**
                         * Closing edge of cycle.
                         */
                        dt = VectorMath.dot(ai, ak);
                        VectorMath.cross(ai, ak, crs);
                    }
                    dt = check(dt);
                    double ang = acos(dt);
                    if (VectorMath.dot(crs, unvect) > 0.0) {
                        ang = -ang;
                    }
                    /**
                     * Add to total for cycle.
                     */
                    totang += ang;
                }
                return totang;
            }

            public void epuclc(double spv[][], int nedge, double epu[][]) {
                double ai[] = new double[3];
                /**
                 * Calculate unit vectors along edges.
                 */
                for (int ke = 0; ke < nedge; ke++) {
                    /**
                     * Get index of second edge of corner.
                     */
                    int ke2;
                    if (ke < nedge) {
                        ke2 = ke + 1;
                    } else {
                        ke2 = 0;
                    }
                    /**
                     * Unit vector along edge of cycle.
                     */
                    for (int k = 0; k < 3; k++) {
                        epu[k][ke] = spv[k][ke2] - spv[k][ke];
                    }
                    getVector(ai, epu, ke);
                    double epun = VectorMath.r(ai);
                    if (epun <= 0.0) {
                        //logger.severe("Null Edge in Cycle");
                        throw new EnergyException("Null Edge in Cycle", true);
                    }
                    /**
                     * Normalize.
                     */
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
                /**
                 * Vectors for null edges come from following or preceding
                 * edges.
                 */
                for (int ke = 0; ke < nedge; ke++) {
                    getVector(ai, epu, ke);
                    if (VectorMath.r(ai) <= 0.0) {
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
            public double measpm(int ifn) {
                double pav[][] = new double[3][3];
                double vect1[] = new double[3];
                double vect2[] = new double[3];
                double vect3[] = new double[3];
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
                VectorMath.cross(vect1, vect2, vect3);
                return height * vect3[2] / 2.0;
            }

            public void measfp(int ifp, double av[]) {
                double ai[] = new double[3];
                double aj[] = new double[3];
                double ak[] = new double[3];
                double vect1[] = new double[3];
                double vect2[] = new double[3];
                double acvect[] = new double[3];
                double aavect[] = new double[3];
                double tanv[][][] = new double[3][2][MAXCYEP];
                double radial[][] = new double[3][MAXCYEP];
                double pcurve = 0.0;
                double gcurve = 0.0;
                int ia = fpa[ifp];
                int ncycle = fpncy[ifp];
                int ieuler;
                if (ncycle > 0) {
                    ieuler = 1 - ncycle;
                } else {
                    ieuler = 1;
                }
                for (int icyptr = 0; icyptr < ncycle; icyptr++) {
                    int icy = fpcy[icyptr][ifp];
                    int nedge = cynep[icy];
                    for (int ke = 0; ke < nedge; ke++) {
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
                        VectorMath.norm(aavect, aavect);
                        double dt = VectorMath.dot(acvect, aavect);
                        double geo = -dt / (radius[ia] * cr[ic]);
                        int iv1 = epv[0][iep];
                        int iv2 = epv[1][iep];
                        double angle;
                        if (iv1 == 0 || iv2 == 0) {
                            angle = 2.0 * PI;
                        } else {
                            for (int k = 0; k < 3; k++) {
                                vect1[k] = v[k][iv1] - c[k][ic];
                                vect2[k] = v[k][iv2] - c[k][ic];
                                radial[k][ke] = v[k][iv1] - a[k][ia];
                            }
                            getVector(ai, radial, ke);
                            VectorMath.norm(ai, ai);
                            getVector(aj, tanv, 0, ke);
                            VectorMath.cross(vect1, aavect, aj);
                            VectorMath.norm(aj, aj);
                            getVector(ak, tanv, 1, ke);
                            VectorMath.cross(vect2, aavect, ak);
                            VectorMath.norm(ak, ak);
                            angle = vecang(vect1, vect2, aavect, -1.0);
                        }
                        gcurve += cr[ic] * angle * geo;
                        if (nedge != 1 && ke > 0) {
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
                    if (nedge > 0) {
                        getVector(ai, tanv, 1, nedge);
                        getVector(aj, tanv, 0, 0);
                        getVector(ak, radial, 0);
                        double angle = vecang(ai, aj, ak, 1.0);
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

            public void measfs(int ifs, double saddle[]) {
                double areas = 0.0;
                double vols = 0.0;
                double areasp = 0.0;
                double volsp = 0.0;
                double vect1[] = new double[3];
                double vect2[] = new double[3];
                double aavect[] = new double[3];
                int iep = fsep[0][ifs];
                int ic = epc[iep];
                int it = ct[ic];
                int ia1 = ta[0][it];
                int ia2 = ta[1][it];
                for (int k = 0; k < 3; k++) {
                    aavect[k] = a[k][ia2] - a[k][ia1];
                }
                VectorMath.norm(aavect, aavect);
                int iv1 = epv[0][iep];
                int iv2 = epv[1][iep];
                double phi;
                if (iv1 == 0 || iv2 == 0) {
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
                double d1 = -1.0 * VectorMath.dot(vect1, aavect);
                double d2 = VectorMath.dot(vect2, aavect);
                theta1 = atan2(d1, tr[it]);
                theta2 = atan2(d2, tr[it]);

                /**
                 * Check for cusps.
                 */
                double thetaq;
                boolean cusp;
                if (tr[it] < probe && theta1 > 0.0 && theta2 > 0.0) {
                    cusp = true;
                    double rat = tr[it] / probe;
                    rat = check(rat);
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
                    //logger.severe("IA1 Inconsistency in MEASFS");
                    throw new EnergyException("IA1 Inconsistency in MEASFS", true);
                }
                for (int k = 0; k < 3; k++) {
                    vect1[k] = c[k][ic1] - a[k][ia1];
                    vect2[k] = c[k][ic2] - a[k][ia2];
                }
                double w1 = VectorMath.dot(vect1, aavect);
                double w2 = -1.0 * VectorMath.dot(vect2, aavect);
                double cone1 = phi * (w1 * cr[ic1] * (w1 * cr[ic1])) / 6.0;
                double cone2 = phi * (w2 * cr[ic2] * (w2 * cr[ic2])) / 6.0;
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

            public void measfn(int ifn, double av[]) {
                double ai[] = new double[3];
                double aj[] = new double[3];
                double ak[] = new double[3];
                double angle[] = new double[3];
                double pvv[][] = new double[3][3];
                double pav[][] = new double[3][3];
                double planev[][] = new double[3][3];
                double arean = 0.0;
                double voln = 0.0;
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
                        VectorMath.norm(ai, ai);
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
                        VectorMath.cross(ai, aj, ak);
                        VectorMath.norm(ak, ak);
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
                            //logger.severe("Negative Angle in MEASFN");
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
            public void vam(double volume, double area) {
                final int maxdot = 1000;
                final int maxop = 100;
                final int nscale = 20;
                int ivs[] = new int[3];
                int ispind[] = new int[3];
                int ispnd2[] = new int[3];
                int ifnop[] = new int[maxop];
                int nlap[] = new int[nfn];
                int enfs[] = new int[5 * nAtoms];
                int fnt[][] = new int[3][nfn];
                int nspt[][] = new int[3][nfn];
                double cenop[][] = new double[3][maxop];
                double sdot[] = new double[3];
                double dotv[] = new double[nscale];
                double tau[] = new double[3];
                double ppm[] = new double[3];
                double xpnt1[] = new double[3];
                double xpnt2[] = new double[3];
                double qij[] = new double[3];
                double vects[][] = new double[3][3];
                double vect1[] = new double[3];
                double vect2[] = new double[3];
                double vect3[] = new double[3];
                double vect4[] = new double[3];
                double vect5[] = new double[3];
                double vect6[] = new double[3];
                double vect7[] = new double[3];
                double vect8[] = new double[3];
                double upp[] = new double[3];
                double thetaq[] = new double[3];
                double sigmaq[] = new double[3];
                double umq[] = new double[3];
                double upq[] = new double[3];
                double uc[] = new double[3];
                double uq[] = new double[3];
                double uij[] = new double[3];
                double ai[] = new double[3];
                double aj[] = new double[3];
                double ak[] = new double[3];
                double al[] = new double[3];
                double dots[][] = new double[3][maxdot];
                double tdots[][] = new double[3][maxdot];
                double atmarea[] = new double[nAtoms];
                double depths[] = new double[nfn];
                double cora[] = new double[nfn];
                double corv[] = new double[nfn];
                double alts[][] = new double[3][nfn];
                double fncen[][] = new double[3][nfn];
                double fnvect[][][] = new double[3][3][nfn];
                boolean ate[] = new boolean[maxop];
                boolean badav[] = new boolean[nfn];
                boolean badt[] = new boolean[nfn];
                boolean fcins[][] = new boolean[3][nfn];
                boolean fcint[][] = new boolean[3][nfn];
                boolean fntrev[][] = new boolean[3][nfn];

                /**
                 * Compute the volume of the interior polyhedron.
                 */
                double polyhedronVolume = 0.0;
                for (int ifn = 0; ifn < nfn; ifn++) {
                    polyhedronVolume += measpm(ifn);
                }

                /**
                 * Compute the area and volume due to convex faces as well as
                 * the area partitioned among the atoms.
                 */
                double totap = 0.0;
                double totvp = 0.0;
                fill(atmarea, 0.0);
                double convexFaces[] = {0.0, 0.0};
                for (int ifp = 0; ifp < nfp; ifp++) {
                    measfp(ifp, convexFaces);
                    int ia = fpa[ifp];
                    atmarea[ia] += convexFaces[0];
                    totap += convexFaces[0];
                    totvp += convexFaces[1];
                }

                /**
                 * Compute the area and volume due to saddle faces as well as
                 * the spindle correction value.
                 */
                double totas = 0.0;
                double totvs = 0.0;
                double totasp = 0.0;
                double totvsp = 0.0;
                double saddle[] = {0.0, 0.0, 0.0, 0.0};
                for (int ifs = 0; ifs < nfs; ifs++) {
                    for (int k = 0; k < 2; k++) {
                        int ien = fsen[k][ifs];
                        if (ien > 0) {
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

                /**
                 * Compute the area and volume due to concave faces.
                 */
                double totan = 0.0;
                double totvn = 0.0;
                double concaveFaces[] = {0.0, 0.0};
                for (int ifn = 0; ifn < nfn; ifn++) {
                    measfn(ifn, concaveFaces);
                    double arean = concaveFaces[0];
                    double voln = concaveFaces[1];
                    totan += arean;
                    totvn += voln;
                }

                /**
                 * Compute the area and volume lens correction values.
                 */
                double alenst = 0.0;
                double alensn = 0.0;
                double vlenst = 0.0;
                double vlensn = 0.0;
                if (probe > 0.0) {
                    int ndots[] = {maxdot};
                    gendot(ndots, dots, probe, 0.0, 0.0, 0.0);
                    double dota = (4.0 * PI * probe * probe) / ndots[0];
                    for (int ifn = 0; ifn < nfn; ifn++) {
                        nlap[ifn] = 0;
                        cora[ifn] = 0.0;
                        corv[ifn] = 0.0;
                        badav[ifn] = false;
                        badt[ifn] = false;
                        for (int k = 0; k < 3; k++) {
                            nspt[k][ifn] = 0;
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
                        //int ia = va[iv];
                        /**
                         * Get vertices and vectors.
                         */
                        for (int ke = 0; ke < 3; ke++) {
                            ien = fnen[ke][ifn];
                            ivs[ke] = env[0][ien];
                            int ia = va[ivs[ke]];
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
                        /**
                         * Calculate normal vectors for the three planes that
                         * cut out the geodesic triangle.
                         */
                        getVector(ai, vects, 0);
                        getVector(aj, vects, 1);
                        getVector(ak, fnvect, 0, ifn);
                        VectorMath.cross(ai, aj, ak);
                        VectorMath.norm(ak, ak);
                        getVector(ai, vects, 2);
                        getVector(ak, fnvect, 1, ifn);
                        VectorMath.cross(aj, ai, ak);
                        VectorMath.norm(ak, ak);
                        getVector(aj, vects, 0);
                        getVector(ak, fnvect, 2, ifn);
                        VectorMath.cross(ai, aj, ak);
                        VectorMath.norm(ak, ak);
                    }
                    for (int ifn = 0; ifn < nfn - 1; ifn++) {
                        for (int jfn = ifn + 1; jfn < nfn; jfn++) {
                            getVector(ai, fncen, ifn);
                            getVector(aj, fncen, jfn);
                            double dij2 = VectorMath.dist2(ai, aj);
                            if (dij2 > 4.0 * probe * probe) {
                                continue;
                            }
                            if (depths[ifn] > probe && depths[jfn] > probe) {
                                continue;
                            }
                            /**
                             * These two probes may have intersecting surfaces.
                             */
                            double dpp = VectorMath.dist(ai, aj);
                            /**
                             * Compute the midpoint.
                             */
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
                            /**
                             * Use circle-place intersection routine.
                             */
                            boolean alli = true;
                            boolean anyi = false;
                            boolean spindl = false;
                            for (int k = 0; k < 3; k++) {
                                ispind[k] = 0;
                                ispnd2[k] = 0;
                            }
                            for (int ke = 0; ke < 3; ke++) {
                                thetaq[ke] = 0.0;
                                sigmaq[ke] = 0.0;
                                tau[ke] = 0.0;
                                getVector(ai, fncen, ifn);
                                getVector(aj, fnvect, ke, ifn);
                                cirpln(ppm, rm, upp, ai, aj, cintp, cinsp, xpnt1, xpnt2);
                                fcins[ke][ifn] = cinsp;
                                fcint[ke][ifn] = cintp;
                                if (!cinsp) {
                                    alli = false;
                                }
                                if (cintp) {
                                    anyi = true;
                                }
                                if (!cintp) {
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
                                if (ispind[ke] == 0) {
                                    continue;
                                }

                                /**
                                 * Check that the two ways of calculating
                                 * intersection points match.
                                 */
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
                                    //qji[k] = t[k][it] + stq * probe * uij[k];
                                }
                                for (int k = 0; k < 3; k++) {
                                    umq[k] = (qij[k] - ppm[k]) / rm;
                                    upq[k] = (qij[k] - fncen[k][ifn]) / probe;
                                }
                                VectorMath.cross(uij, upp, vect1);
                                double dt = VectorMath.dot(umq, vect1);
                                check(dt);
                                sigmaq[ke] = acos(dt);
                                getVector(ai, fnvect, ke, ifn);
                                VectorMath.cross(upq, ai, vect1);
                                VectorMath.norm(vect1, uc);
                                VectorMath.cross(upp, upq, vect1);
                                VectorMath.norm(vect1, uq);
                                dt = VectorMath.dot(uc, uq);
                                check(dt);
                                tau[ke] = PI - acos(dt);
                            }
                            boolean allj = true;
                            boolean anyj = false;
                            for (int ke = 0; ke < 3; ke++) {
                                getVector(ai, fncen, jfn);
                                getVector(aj, fnvect, ke, jfn);
                                cirpln(ppm, rm, upp, ai, aj, cinsp, cintp, xpnt1, xpnt2);
                                fcins[ke][jfn] = cinsp;
                                fcint[ke][jfn] = cintp;
                                if (!cinsp) {
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
                            /**
                             * This kind of overlap can be handled.
                             */
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
                                    if (ispind[ke] == ispnd2[ke2] || ispind[ke] == 0) {
                                        continue;
                                    }
                                    getVector(ai, fncen, ifn);
                                    getVector(aj, fnvect, ke, ifn);
                                    getVector(ak, fncen, jfn);
                                    getVector(al, fnvect, ke2, jfn);
                                    cirpln(ai, probe, aj, ak, al, cinsp, cintp, xpnt1, xpnt2);
                                    if (!cintp) {
                                        continue;
                                    }
                                    ien = fnen[ke2][jfn];
                                    iv1 = env[0][ien];
                                    iv2 = env[1][ien];
                                    for (int k = 0; k < 3; k++) {
                                        vect7[k] = v[k][iv1] - fncen[k][jfn];
                                        vect8[k] = v[k][iv2] - fncen[k][jfn];
                                    }
                                    /**
                                     * Check whether point lies on spindle arc.
                                     */
                                    for (int k = 0; k < 3; k++) {
                                        vect1[k] = xpnt1[k] - fncen[k][ifn];
                                        vect2[k] = xpnt2[k] - fncen[k][ifn];
                                        vect5[k] = xpnt1[k] - fncen[k][jfn];
                                        vect6[k] = xpnt2[k] - fncen[k][jfn];
                                    }
                                    /**
                                     * Continue to next if statement if any of
                                     * the following are true.
                                     */
                                    getVector(ai, fnvect, ke, ifn);
                                    getVector(aj, fnvect, ke2, jfn);
                                    if (triple(vect3, vect1, ai) < 0.0
                                            || triple(vect1, vect4, ai) < 0.0
                                            || triple(vect7, vect5, aj) < 0.0
                                            || triple(vect5, vect8, aj) < 0.0) {
                                        if (!((triple(vect3, vect2, ai) < 0.0
                                                || triple(vect2, vect4, ai) < 0.0
                                                || triple(vect7, vect6, aj) < 0.0
                                                || triple(vect6, vect8, aj) < 0.0))) {
                                            badav[ifn] = true;
                                        }
                                    } else {
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
                                    if (ispind[ke] == ispnd2[ke2] || ispnd2[ke2] == 0) {
                                        continue;
                                    }
                                    getVector(ai, fncen, ifn);
                                    getVector(aj, fnvect, ke, ifn);
                                    getVector(ak, fncen, jfn);
                                    getVector(al, fnvect, ke2, jfn);
                                    cirpln(ak, probe, al, ai, aj, cinsp, cintp, xpnt1, xpnt2);
                                    if (!cintp) {
                                        continue;
                                    }
                                    ien = fnen[ke2][jfn];
                                    iv1 = env[0][ien];
                                    iv2 = env[1][ien];
                                    for (int k = 0; k < 3; k++) {
                                        vect7[k] = v[k][iv1] - fncen[k][jfn];
                                        vect8[k] = v[k][iv2] - fncen[k][jfn];
                                    }
                                    /**
                                     * Check whether point lies on spindle arc.
                                     */
                                    for (int k = 0; k < 3; k++) {
                                        vect1[k] = xpnt1[k] - fncen[k][ifn];
                                        vect2[k] = xpnt2[k] - fncen[k][ifn];
                                        vect5[k] = xpnt1[k] - fncen[k][jfn];
                                        vect6[k] = xpnt2[k] - fncen[k][jfn];
                                    }
                                    /**
                                     * Continue to next if statement if any of
                                     * the following are true.
                                     */
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

                            /**
                             * Check for vertex on opposing probe in face.
                             */
                            /**
                             * For (int kv = 0; kv < 3; kv++) { vip[kv] = false;
                             * int ien = fnen[kv][jfn]; iv = env[0][ien]; for
                             * (int k = 0; k < 3; k++) { vect1[k] = v[k][iv] -
                             * fncen[k][ifn]; } VectorMath.norm(vect1, vect1);
                             * for (int ke = 0; ke < 3; ke++) { getVector(ai,
                             * fnvect, ke, ifn); getVector(aj, v, iv); double dt
                             * = VectorMath.dot(ai, aj); if (dt > 0.0) { move =
                             * true; break; } } if (!move) { vip[kv] = true; }
                             * move = false; }
                             */
                        }
                    }
                    for (int ifn = 0; ifn < nfn; ifn++) {
                        for (int ke = 0; ke < 3; ke++) {
                            if (nspt[ke][ifn] > 1) {
                                badt[ifn] = true;
                            }
                        }
                    }
                    for (int ifn = 0; ifn < nfn; ifn++) {
                        if (nlap[ifn] <= 0) {
                            continue;
                        }
                        /**
                         * Gather all overlapping probes.
                         */
                        int nop = 0;
                        for (int jfn = 0; jfn < nfn; jfn++) {
                            if (ifn != jfn) {
                                getVector(ai, fncen, ifn);
                                getVector(aj, fncen, jfn);
                                double dij2 = VectorMath.dist2(ai, aj);
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
                            /**
                             * Numerical calculation of the correction.
                             */
                            double areado = 0.0;
                            double voldo = 0.0;
                            double scinc = 1.0 / nscale;
                            for (int isc = 0; isc < nscale; isc++) {
                                double rsc = isc - 0.5;
                                dotv[isc] = probe * dota * rsc * rsc * scinc * scinc * scinc;
                            }
                            for (int iop = 0; iop < nop; iop++) {
                                ate[iop] = false;
                            }
                            int neatmx = 0;
                            for (int idot = 0; idot < ndots[0]; idot++) {
                                boolean move = false;
                                for (int ke = 0; ke < 3; ke++) {
                                    getVector(ai, fnvect, ke, ifn);
                                    getVector(aj, dots, idot);
                                    double dt = VectorMath.dot(ai, aj);
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
                                for (int iop = 0; iop < nop; iop++) {
                                    jfn = ifnop[iop];
                                    getVector(ai, dots, idot);
                                    getVector(aj, fncen, jfn);
                                    double ds2 = VectorMath.dist2(ai, aj);
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
                                    for (int iop = 0; iop < nop; iop++) {
                                        jfn = ifnop[iop];
                                        getVector(ai, fncen, jfn);
                                        double ds2 = VectorMath.dist2(sdot, ai);
                                        if (ds2 > probe * probe) {
                                            for (int k = 0; k < 3; k++) {
                                                vect1[k] = sdot[k] - fncen[k][jfn];
                                            }
                                            for (int ke = 0; ke < 3; ke++) {
                                                getVector(ai, fnvect, ke, jfn);
                                                double dt = VectorMath.dot(ai, vect1);
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
                            for (int iop = 0; iop < nop; iop++) {
                                if (ate[iop]) {
                                    nate++;
                                }
                            }
                            /**
                             * Use either the analytical or numerical
                             * correction.
                             */
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
                /**
                 * Finally, compute the total area and total volume.
                 */
                logger.info(String.format("totap=%16.8f,totas=%16.8f,totan=%16.8f,totasp=%16.8f,alenst=%16.8f", totap, totas, totan, totasp, alenst));
                area = totap + totas + totan - totasp - alenst;
                logger.info(String.format("totvp=%16.8f,totvs=%16.8f,totvn=%16.8f,hedron=%16.8f,totvsp=%16.8f,vlenst=%16.8f", totvp, totvs, totvn, polyhedronVolume, totvsp, vlenst));
                volume = totvp + totvs + totvn + polyhedronVolume - totvsp + vlenst;
                logger.info(String.format("volume=%16.8f area= %16.8f", volume, area));
            }

            /**
             * The gendot method finds the coordinates of a specified number of
             * surface points for a sphere with the input radius and coordinate
             * center.
             */
            public void gendot(int ndots[], double dots[][], double radius,
                    double xcenter, double ycenter, double zcenter) {
                int nequat = (int) sqrt(PI * (double) ndots[0]);
                int nvert = nequat / 2;
                if (nvert < 0) {
                    nvert = 0;
                }
                int k = 0;
                for (int i = -1; i < nvert; i++) {
                    double fi = (PI * (double) i) / (double) nvert;
                    double z = cos(fi);
                    double xy = sin(fi);
                    int nhoriz = (int) (nequat * xy);
                    if (nhoriz < 0) {
                        nhoriz = 0;
                    }
                    for (int j = -1; j < nhoriz - 1; j++) {
                        double fj = (2.0 * PI * (double) (j)) / (double) (nhoriz);
                        double x = cos(fj) * xy;
                        double y = sin(fj) * xy;
                        k++;
                        dots[0][k] = x * radius + xcenter;
                        dots[1][k] = y * radius + ycenter;
                        dots[2][k] = z * radius + zcenter;
                        if (k >= ndots[0]) {
                            ndots[0] = k;
                            return;
                        }
                    }
                }
                ndots[0] = k;
            }

            /**
             * The cirpln method determines the points of intersection between a
             * specified circle and plane.
             */
            public boolean cirpln(double circen[], double cirrad, double cirvec[], double plncen[],
                    double plnvec[], boolean cinsp, boolean cintp, double xpnt1[], double xpnt2[]) {
                double cpvect[] = new double[3];
                double pnt1[] = new double[3];
                double vect1[] = new double[3];
                double vect2[] = new double[3];
                double uvect1[] = new double[3];
                double uvect2[] = new double[3];
                for (int k = 0; k < 3; k++) {
                    cpvect[k] = plncen[k] - circen[k];
                }
                double dcp = VectorMath.dot(cpvect, plnvec);
                cinsp = (dcp > 0.0);
                VectorMath.cross(plnvec, cirvec, vect1);
                if (VectorMath.r(vect1) > 0.0) {
                    VectorMath.norm(vect1, uvect1);
                    VectorMath.cross(cirvec, uvect1, vect2);
                    if (VectorMath.r(vect2) > 0.0) {
                        VectorMath.norm(vect2, uvect2);
                        double dir = VectorMath.dot(uvect2, plnvec);
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
            public double vecang(double v1[], double v2[], double axis[], double hand) {
                double a1 = VectorMath.r(v1);
                double a2 = VectorMath.r(v2);
                double dt = VectorMath.dot(v1, v2);
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

            public double depth(int ip, double alt[]) {
                double vect1[] = new double[3];
                double vect2[] = new double[3];
                double vect3[] = new double[3];
                double vect4[] = new double[3];
                int ia1 = pa[0][ip];
                int ia2 = pa[1][ip];
                int ia3 = pa[2][ip];
                for (int k = 0; k < 3; k++) {
                    vect1[k] = a[k][ia1] - a[k][ia3];
                    vect2[k] = a[k][ia2] - a[k][ia3];
                    vect3[k] = p[k][ip] - a[k][ia3];
                }
                VectorMath.cross(vect1, vect2, vect4);
                VectorMath.norm(vect4, vect4);
                double dot = VectorMath.dot(vect4, vect3);
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

            public void calcDerivative(int lb, int ub) {
                /**
                 * Fix the step size in the z-direction; this value sets the
                 * accuracy of the numerical derivatives; zstep=0.06 is a good
                 * balance between compute time and accuracy.
                 */
                zstep = 0.0601;
                /**
                 * Load the cubes based on coarse lattice; first of all set edge
                 * length to the maximum diameter of any atom.
                 */
                edge = 2.0 * rmax;
                nx = (int) ((xmax - xmin) / edge);
                ny = (int) ((ymax - ymin) / edge);
                nz = (int) ((zmax - zmin) / edge);
                if (max(max(nx, ny), nz) > mxcube) {
                    //logger.severe(" VOLUME1  --  Increase the Value of MAXCUBE");
                    throw new EnergyException(" VOLUME1  --  Increase the Value of MAXCUBE", false);
                }
                /**
                 * Initialize the coarse lattice of cubes.
                 */
                for (int i = 0; i <= nx; i++) {
                    for (int j = 0; j <= ny; j++) {
                        for (int k = 0; k <= nz; k++) {
                            cube[0][i][j][k] = 0;
                            cube[1][i][j][k] = -1;
                        }
                    }
                }
                /**
                 * Find the number of atoms in each cube.
                 */
                for (int m = 0; m < nAtoms; m++) {
                    if (!skip[m]) {
                        int i = (int) ((x[m] - xmin) / edge);
                        int j = (int) ((y[m] - ymin) / edge);
                        int k = (int) ((z[m] - zmin) / edge);
                        cube[0][i][j][k]++;
                    }
                }

                /**
                 * Determine the highest index in the array "itab" for the atoms
                 * that fall into each cube; the first cube that has atoms
                 * defines the first index for "itab"; the final index for the
                 * atoms in the present cube is the final index of the last cube
                 * plus the number of atoms in the present cube.
                 */
                isum = 0;
                for (int i = 0; i <= nx; i++) {
                    for (int j = 0; j <= ny; j++) {
                        for (int k = 0; k <= nz; k++) {
                            tcube = cube[0][i][j][k];
                            if (tcube != 0) {
                                isum += tcube;
                                cube[1][i][j][k] = isum - 1;
                            }
                        }
                    }
                }

                /**
                 * "cube(1,,,)" now contains a pointer to the array "itab"
                 * giving the position of the last entry for the list of atoms
                 * in that cube of total number equal to "cube(0,,,)".
                 */
                for (int m = 0; m < nAtoms; m++) {
                    if (!skip[m]) {
                        int i = (int) ((x[m] - xmin) / edge);
                        int j = (int) ((y[m] - ymin) / edge);
                        int k = (int) ((z[m] - zmin) / edge);
                        tcube = cube[1][i][j][k];
                        itab[tcube] = m;
                        cube[1][i][j][k]--;
                    }
                }

                /**
                 * Set "cube(1,,,)" to be the starting index in "itab" for atom
                 * list of that cube; and "cube(0,,,)" to be the stop index.
                 */
                isum = 0;
                for (int i = 0; i <= nx; i++) {
                    for (int j = 0; j <= ny; j++) {
                        for (int k = 0; k <= nz; k++) {
                            tcube = cube[0][i][j][k];
                            //logger.info(String.format(" TCUBE %d %d %d %d", i, j, k, tcube));
                            if (tcube != 0) {
                                isum += tcube;
                                cube[0][i][j][k] = isum - 1;
                                cube[1][i][j][k]++;
                            }
                        }
                    }
                }

                /**
                 * Process in turn each atom from the coordinate list; first
                 * select the potential intersecting atoms.
                 */
                for (ir = 0; ir < nAtoms; ir++) {
                    pre_dx = 0.0;
                    pre_dy = 0.0;
                    pre_dz = 0.0;
                    if (skip[ir]) {
                        continue;
                    }
                    rr = vdwrad[ir];
                    rrx2 = 2.0 * rr;
                    rrsq = rr * rr;
                    xr = x[ir];
                    yr = y[ir];
                    zr = z[ir];
                    /**
                     * Find cubes to search for overlaps for current atom.
                     */
                    istart = (int) ((xr - xmin) / edge);
                    istop = min(istart + 2, nx + 1);
                    istart = max(istart, 1);
                    jstart = (int) ((yr - ymin) / edge);
                    jstop = min(jstart + 2, ny + 1);
                    jstart = max(jstart, 1);
                    kstart = (int) ((zr - zmin) / edge);
                    kstop = min(kstart + 2, nz + 1);
                    kstart = max(kstart, 1);
                    /**
                     * Load all overlapping atoms into "inov".
                     */
                    io = -1;
                    //logger.info(String.format(" %d %d %d %d %d %d %d ", ir, istart, istop, jstart, jstop, kstart, kstop));
                    for (int i = istart - 1; i < istop; i++) {
                        for (int j = jstart - 1; j < jstop; j++) {
                            for (int k = kstart - 1; k < kstop; k++) {
                                mstart = cube[1][i][j][k];
                                if (mstart != -1) {
                                    mstop = cube[0][i][j][k];
                                    for (int m = mstart; m <= mstop; m++) {
                                        in = itab[m];
                                        //logger.info(String.format(" CHECK %d %d", ir, in));
                                        if (in != ir) {
                                            io++;
                                            if (io > MAXARC) {
                                                //logger.severe(" VOLUME1  --  Increase the Value of MAXARC");
                                            }
                                            dx[io] = x[in] - xr;
                                            dy[io] = y[in] - yr;
                                            dsq[io] = (dx[io] * dx[io]) + (dy[io] * dy[io]);
                                            dist2 = dsq[io] + ((z[in] - zr) * (z[in] - zr));
                                            vdwsum = (rr + vdwrad[in]) * (rr + vdwrad[in]);
                                            if (dist2 > vdwsum || dist2 == 0.0) {
                                                io--;
                                            } else {
                                                d[io] = sqrt(dsq[io]);
                                                inov[io] = in;
                                                //logger.info(String.format(" INIT %d %d %d %16.8f", ir, io, in, d[io]));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    //logger.info(String.format("ir %d io %d", ir, io));
                    /**
                     * Determine resolution along the z-axis.
                     */
                    if (io != -1) {
                        ztop = zr + rr;
                        ztopshave = ztop - zstep;
                        zgrid = zr - rr;
                        /**
                         * Half of the part not covered by the planes.
                         */
                        zgrid += 0.5 * (rrx2 - ((int) (rrx2 / zstep) * zstep));
                        zstart = zgrid;
                        /**
                         * Section atom spheres perpendicular to the z-axis.
                         */
                        while (zgrid <= ztop) {
                            /**
                             * "rsecr" is radius of circle of intersection of
                             * "ir" sphere on the current sphere.
                             */
                            rsec2r = rrsq - ((zgrid - zr) * (zgrid - zr));
                            if (rsec2r < 0.0) {
                                rsec2r = 0.000001;
                            }
                            rsecr = sqrt(rsec2r);
                            if (zgrid >= ztopshave) {
                                cos_phi1 = 1.0;
                                phi1 = 0.0;
                            } else {
                                cos_phi1 = (zgrid + (0.5 * zstep) - zr) / rr;
                                phi1 = acos(cos_phi1);
                            }
                            if (zgrid == zstart) {
                                cos_phi2 = -1.0;
                                phi2 = PI;
                            } else {
                                cos_phi2 = (zgrid - (0.5 * zstep) - zr) / rr;
                                phi2 = acos(cos_phi2);
                            }
                            /**
                             * Check intersection of neighbor cirlces.
                             */
                            narc = -1;
                            for (int k = 0; k <= io; k++) {
                                in = inov[k];
                                rinsq = vdwrad[in] * vdwrad[in];
                                rsec2n = rinsq - ((zgrid - z[in]) * (zgrid - z[in]));
                                //logger.info(String.format(" NARC %d %d %16.8f %16.8f %16.8f", ir, k, rinsq, z[in], zgrid));
                                if (rsec2n > 0.0) {
                                    rsecn = sqrt(rsec2n);
                                    if (d[k] < (rsecr + rsecn)) {
                                        rdiff = rsecr - rsecn;
                                        //logger.info(String.format(" DIFF %d %d %16.8f %16.8f %16.8f", ir, k, d[k], rsecr, rsecn));
                                        if (d[k] <= abs(rdiff)) {
                                            if (rdiff < 0.0) {
                                                narc = 0;
                                                arci[narc] = 0.0;
                                                arcf[narc] = pix2;
                                            }
                                            //logger.info(String.format("%d Continue", ir));
                                            continue;
                                        }
                                        narc++;
                                        if (narc > MAXARC) {
                                            //logger.info("VOLUME1 -- Increase the Value of MAXARC");
                                        }
                                        /**
                                         * Initial and final arc endpoints are
                                         * found for intersection of "ir" circle
                                         * with another circle contained in same
                                         * plane; the initial endpoint of the
                                         * enclosed arc is stored in "arci", the
                                         * final endpoint in "arcf"; get
                                         * "cosine" via law of cosines.
                                         */
                                        cosine = (dsq[k] + rsec2r - rsec2n) / (2.0 * d[k] * rsecr);
                                        cosine = min(1.0, max(-1.0, cosine));
                                        /**
                                         * "alpha" is the angle between a line
                                         * containing either point of
                                         * intersection and the reference circle
                                         * center and the line containing both
                                         * circle centers; "beta" is the angle
                                         * between the line containing both
                                         * circle centers and x-axis.
                                         */
                                        alpha = acos(cosine);
                                        beta = atan2(dy[k], dx[k]);
                                        if (dy[k] < 0.0) {
                                            beta += pix2;
                                        }
                                        ti = beta - alpha;
                                        tf = beta + alpha;
                                        if (ti < 0.0) {
                                            ti += pix2;
                                        }
                                        if (tf > pix2) {
                                            tf -= pix2;
                                        }
                                        arci[narc] = ti;
                                        /**
                                         * If the arc crosses zero, then it is
                                         * broken into two segments; the first
                                         * ends at two pi and the second begins
                                         * at zero.
                                         */
                                        if (tf < ti) {
                                            arcf[narc] = pix2;
                                            narc++;
                                            arci[narc] = 0.0;
                                            //logger.info(String.format("ir= %d narc= %d arci= %16.8f", ir, narc, arci[narc]));
                                        }
                                        arcf[narc] = tf;
                                        //logger.info(String.format(" ARCF %d %d %16.8f %16.8f", ir, narc, tf, ti));
                                        // BELOW HERE
                                    }
                                }
                            }

                            /**
                             * Find the pre-area and pre-forces on this section
                             * (band), "pre-" means a multiplicative factor is
                             * yet to be applied.
                             */
                            if (narc == -1) {
                                //logger.info(String.format(" %d cos_phi %16.8f %16.8f %16.8f %16.8f",
                                //ir, phi1, cos_phi1, phi2, cos_phi2));
                                seg_dz = pix2 * ((cos_phi1 * cos_phi1) - (cos_phi2 * cos_phi2));
                                pre_dz += seg_dz;
                                //logger.info(String.format(" seg_dx %16.8f pre_dz %16.8f ", seg_dz, pre_dz));
                            } else {
                                /**
                                 * Sort the arc endpoint arrays, each with
                                 * "narc" entries, in order of increasing values
                                 * of the arguments in "arci".
                                 */
                                for (int k = 0; k < narc; k++) {
                                    aa = arci[k];
                                    bb = arcf[k];
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
                                /**
                                 * Consolidate arcs by removing overlapping arc
                                 * endpoints.
                                 */
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
                                        // SOME PRINTS ARE WRONG AFTER FIRST ENTRY
                                        //logger.info(String.format(" ARCF1 %d %16.8f", ir, arcf[0]));
                                    }
                                }

                                // SOME OF THE FOLLOWING PRINTS ARE WRONG
                                //logger.info(String.format(" SORT %d %d %16.8f %16.8f", ir, narc, arci[0], arcf[0]));
                                /**
                                 * Compute the numerical pre-derivative values.
                                 */
                                for (int k = 0; k <= narc; k++) {
                                    theta1 = arci[k];
                                    theta2 = arcf[k];
                                    //logger.info(String.format("%d theta1=%16.8f theta2=%16.8f", ir, theta1, theta2));
                                    if (theta2 > theta1) {
                                        dtheta = theta2 - theta1;
                                    } else {
                                        dtheta = (theta2 + pix2) - theta1;
                                    }
                                    phi_term = phi2 - phi1 - 0.5 * (sin(2.0 * phi2) - sin(2.0 * phi1));
                                    seg_dx = (sin(theta2) - sin(theta1)) * phi_term;
                                    seg_dy = (cos(theta1) - cos(theta2)) * phi_term;
                                    seg_dz = dtheta * ((cos_phi1 * cos_phi1) - (cos_phi2 * cos_phi2));
                                    pre_dx += seg_dx;
                                    pre_dy += seg_dy;
                                    pre_dz += seg_dz;
                                    // SOME OF THE FOLLOWING PRINTS ARE WRONG
                                    //logger.info(String.format(" FINAL %d %16.8f %16.8f %16.8f", ir, seg_dx, seg_dy, seg_dz));
                                }
                                //logger.info(String.format(" LAST %d %16.8f %16.8f %16.8f", ir, pre_dx, pre_dy, pre_dz));
                            }
                            zgrid += zstep;
                        }
                    }
                    dex[0][ir] = 0.5 * rrsq * pre_dx;
                    dex[1][ir] = 0.5 * rrsq * pre_dy;
                    dex[2][ir] = 0.5 * rrsq * pre_dz;
                    //logger.info(String.format(" de/dx %d %16.8f %16.8f %16.8f", ir, dex[0][ir], dex[1][ir], dex[2][ir]));
                }
            }

            @Override
            public void run(int lb, int ub) {
                setRadius();
                calcVolume();
                calcDerivative(lb, ub);
            }
        }
    }

    /**
     * Compute Hydrophobic PMF.
     *
     * @since 1.0
     */
    private class HydrophobicPMFRegion extends ParallelRegion {

        // Radius of a carbon atom.
        private final double rCarbon = 1.7;
        // Radius of a water molecule.
        private final double rWater = 1.4;
        // Constant for calculation of atomic surface area.
        private final double safact = 0.3516;
        // Surface area of a hydrophobic carbon atom.
        private final double acSurf = 120.7628;
        // tanh slope (set very steep).
        private final double tSlope = 100.0;
        // Shift the tanh plot along the x-axis.
        private final double tOffset = 6.0;
        // Cutoff distance for pairwise HPMF interactions.
        private final double hpmfCut = 11.0;
        // Cutoff squared
        private final double hpmfCut2 = hpmfCut * hpmfCut;
        // Hydrophobic PMF well depth parameter.
        private final double h1 = -0.7308004860404441194;
        private final double h2 = 0.2001645051578760659;
        private final double h3 = -0.0905499953418473502;
        // Hydrophobic PMF well center point.
        private final double c1 = 3.8167879266271396155;
        private final double c2 = 5.4669162286016419472;
        private final double c3 = 7.1167694861385353278;
        // Reciprocal of the hydrophobic PMF well width.
        private final double w1 = 1.6858993102248638341;
        private final double w2 = 1.3906405621629980285;
        private final double w3 = 1.5741657341338335385;
        private final double rSurf = rCarbon + 2.0 * rWater;
        private final double piSurf = PI * (rCarbon + rWater);
        // Radius of each atom for use with hydrophobic PMF.
        private final double rPMF[];
        // Number of hydrophobic carbon atoms in the system.
        private final int nCarbon;
        // Number of the atom for each HPMF carbon atom site.
        private final int iCarbon[];
        // SASA value for each hydrophobic PMF carbon atom
        private final double carbonSASA[];
        // SASA value
        private final double tanhSA[];
        // Loop to find the SASA value of each hydrophobic carbon.
        private final CarbonSASALoop carbonSASALoop[];
        // Loop to find the hydrophobic energy.
        private final HydrophobicPMFLoop hydrophobicPMFLoop[];
        // Loop to find the SASA chain rule derivatives.
        private final CarbonSASACRLoop carbonSASACRLoop[];
        // Shared energy variable.
        private final SharedDouble sharedEnergy;
        private boolean gradient;
        private final double dtanhSA[];
        private final double sasa[];
        private final double carbonSASACR[];

        public HydrophobicPMFRegion(int nt) {
            logger.info(String.format(" Hydrophobic PMF cut-off:              %8.2f (A)", hpmfCut));
            /**
             * Count hydrophobic carbons.
             */
            int count = 0;
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                int atomicNumber = atom.getAtomicNumber();
                if (atomicNumber == 6) {
                    List<Bond> bonds = atom.getBonds();
                    int bondCount = bonds.size();
                    if (bondCount <= 2) {
                        continue;
                    }
                    boolean keep = true;
                    for (int j = 0; j < bondCount; j++) {
                        Atom atom2 = bonds.get(j).get1_2(atom);
                        int atomicNumber2 = atom2.getAtomicNumber();
                        if (bondCount == 3 && atomicNumber2 == 8) {
                            keep = false;
                            break;
                        }
                    }
                    if (keep) {
                        count++;
                    }
                }
            }
            /**
             * Allocate arrays.
             */
            rPMF = new double[nAtoms];
            nCarbon = count;
            iCarbon = new int[nCarbon];
            carbonSASA = new double[nCarbon];
            carbonSASACR = new double[nCarbon];
            tanhSA = new double[nCarbon];
            dtanhSA = new double[nCarbon];
            sasa = new double[nCarbon];

            carbonSASALoop = new CarbonSASALoop[nt];
            hydrophobicPMFLoop = new HydrophobicPMFLoop[nt];
            carbonSASACRLoop = new CarbonSASACRLoop[nt];
            for (int i = 0; i < nt; i++) {
                carbonSASALoop[i] = new CarbonSASALoop();
                hydrophobicPMFLoop[i] = new HydrophobicPMFLoop();
                carbonSASACRLoop[i] = new CarbonSASACRLoop();
            }
            sharedEnergy = new SharedDouble();

            /**
             * Assign hydrophobic carbon values.
             */
            int index = 0;
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                int atomicNumber = atom.getAtomicNumber();
                if (atomicNumber == 6) {
                    int nh = 0;
                    List<Bond> bonds = atom.getBonds();
                    int bondCount = bonds.size();
                    if (bondCount <= 2) {
                        continue;
                    }
                    boolean keep = true;
                    for (int j = 0; j < bondCount; j++) {
                        Atom atom2 = bonds.get(j).get1_2(atom);
                        int atomicNumber2 = atom2.getAtomicNumber();
                        if (atomicNumber2 == 1) {
                            nh++;
                        }
                        if (bondCount == 3 && atomicNumber2 == 8) {
                            keep = false;
                        }
                    }
                    if (keep) {
                        iCarbon[index] = i;
                        carbonSASA[index] = 1.0;
                        if (bondCount == 3 && nh == 0) {
                            carbonSASA[index] = 1.554;
                        } else if (bondCount == 3 && nh == 1) {
                            carbonSASA[index] = 1.073;
                        } else if (bondCount == 4 && nh == 1) {
                            carbonSASA[index] = 1.276;
                        } else if (bondCount == 4 && nh == 2) {
                            carbonSASA[index] = 1.045;
                        } else if (bondCount == 4 && nh == 3) {
                            carbonSASA[index] = 0.880;
                        }
                        carbonSASA[index] = carbonSASA[index] * safact / acSurf;
                        if (logger.isLoggable(Level.FINEST)) {
                            logger.finest(String.format(" %d Base HPMF SASA for atom %d: %10.8f",
                                    index + 1, i + 1, carbonSASA[index]));
                        }
                        index++;
                    }
                }
            }

            /**
             * Assign HPMF atomic radii from traditional Bondi values
             */
            for (int i = 0; i < nAtoms; i++) {
                rPMF[i] = 2.0;
                int atmnum = atoms[i].getAtomicNumber();
                switch (atmnum) {
                    case 0:
                        rPMF[i] = 0.0;
                        break;
                    case 1:
                        rPMF[i] = 1.20;
                        break;
                    case 2:
                        rPMF[i] = 1.40;
                        break;
                    case 5:
                        rPMF[i] = 1.80;
                        break;
                    case 6:
                        rPMF[i] = 1.70;
                        break;
                    case 7:
                        rPMF[i] = 1.55;
                        break;
                    case 8:
                        rPMF[i] = 1.50;
                        break;
                    case 9:
                        rPMF[i] = 1.47;
                        break;
                    case 10:
                        rPMF[i] = 1.54;
                        break;
                    case 14:
                        rPMF[i] = 2.10;
                        break;
                    case 15:
                        rPMF[i] = 1.80;
                        break;
                    case 16:
                        rPMF[i] = 1.80;
                        break;
                    case 17:
                        rPMF[i] = 1.75;
                        break;
                    case 18:
                        rPMF[i] = 1.88;
                        break;
                    case 34:
                        rPMF[i] = 1.90;
                        break;
                    case 35:
                        rPMF[i] = 1.85;
                        break;
                    case 36:
                        rPMF[i] = 2.02;
                        break;
                    case 53:
                        rPMF[i] = 1.98;
                        break;
                    case 54:
                        rPMF[i] = 2.16;
                        break;
                }
            }
        }

        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }

        public double getEnergy() {
            return sharedEnergy.get();
        }

        @Override
        public void start() {
            sharedEnergy.set(0);
        }

        @Override
        public void run() {
            int ti = getThreadIndex();
            try {
                execute(0, nCarbon - 1, carbonSASALoop[ti]);
                execute(0, nCarbon - 2, hydrophobicPMFLoop[ti]);
                if (gradient) {
                    execute(0, nCarbon - 1, carbonSASACRLoop[ti]);
                }
            } catch (Exception e) {
                String message = "Fatal exception computing Born radii in thread " + ti + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute Hydrophobic PMF radii.
         *
         * @since 1.0
         */
        private class CarbonSASALoop extends IntegerForLoop {

            @Override
            public void run(int lb, int ub) {
                /**
                 * Get the surface area for each hydrophobic carbon atom.
                 */
                for (int ii = lb; ii <= ub; ii++) {
                    final int i = iCarbon[ii];
                    if (!use[i]) {
                        continue;
                    }
                    final double carbonSA = carbonSASA[ii];
                    double sa = acSurf;
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    int count = 0;
                    for (int k = 0; k < nAtoms; k++) {
                        if (i != k && use[k]) {
                            final double xr = x[k] - xi;
                            final double yr = y[k] - yi;
                            final double zr = z[k] - zi;
                            final double r2 = xr * xr + yr * yr + zr * zr;
                            double rk = rPMF[k];
                            double rBig = rk + rSurf;
                            if (r2 < rBig * rBig) {
                                final double r = sqrt(r2);
                                final double rSmall = rk - rCarbon;
                                final double part = piSurf * (rBig - r) * (1.0 + rSmall / r);
                                sa *= (1.0 - carbonSA * part);
                                count++;
                            }
                        }
                    }
                    sasa[ii] = sa;
                    //sasa[ii] = carbonSA;
                    double tSA = tanh(tSlope * (sa - tOffset));
                    tanhSA[ii] = 0.5 * (1.0 + tSA);
                    dtanhSA[ii] = 0.5 * tSlope * (1.0 - tSA * tSA);
                }
            }
        }

        /**
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class HydrophobicPMFLoop extends IntegerForLoop {

            private double energy;
            // Omit
            private final int omit[];
            private double gX[];
            private double gY[];
            private double gZ[];

            public HydrophobicPMFLoop() {
                omit = new int[nAtoms];
            }

            @Override
            public void start() {
                energy = 0.0;
                for (int i = 0; i < nAtoms; i++) {
                    omit[i] = -1;
                }
                int threadID = getThreadIndex();
                gX = grad[threadID][0];
                gY = grad[threadID][1];
                gZ = grad[threadID][2];
                if (gradient) {
                    for (int i = 0; i < nCarbon; i++) {
                        carbonSASACR[i] = 0.0;
                    }
                }
            }

            @Override
            public void run(int lb, int ub) {
                /**
                 * Hydrophobic PME energy.
                 */
                for (int ii = lb; ii <= ub; ii++) {
                    final int i = iCarbon[ii];
                    if (!use[i]) {
                        continue;
                    }
                    double tanhSAi = tanhSA[ii];
                    Atom ssAtom = null;
                    Atom atom = atoms[i];
                    List<Bond> bonds = atom.getBonds();
                    for (Bond bond : bonds) {
                        Atom atom2 = bond.get1_2(atom);
                        int k = atom2.getIndex() - 1;
                        if (!use[k]) {
                            continue;
                        }
                        omit[k] = i;
                        if (atom2.getAtomicNumber() == 16) {
                            ssAtom = atom2;
                        }
                    }
                    List<Angle> angles = atom.getAngles();
                    for (Angle angle : angles) {
                        Atom atom2 = angle.get1_3(atom);
                        if (atom2 != null) {
                            int k = atom2.getIndex() - 1;
                            if (!use[k]) {
                                continue;
                            }
                            omit[k] = i;
                        }
                    }
                    List<Torsion> torsions = atom.getTorsions();
                    for (Torsion torsion : torsions) {
                        Atom atom2 = torsion.get1_4(atom);
                        if (atom2 != null) {
                            int k = atom2.getIndex() - 1;
                            if (!use[k]) {
                                continue;
                            }
                            omit[k] = i;
                            if (ssAtom != null) {
                                List<Bond> bonds2 = atom2.getBonds();
                                for (Bond bond : bonds2) {
                                    Atom s = bond.get1_2(atom2);
                                    if (s.getAtomicNumber() == 16) {
                                        List<Bond> sBonds = s.getBonds();
                                        for (Bond sBond : sBonds) {
                                            Atom s2 = sBond.get1_2(s);
                                            if (s2 == ssAtom) {
                                                omit[k] = -1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    double e = 0.0;
                    for (int kk = ii + 1; kk < nCarbon; kk++) {
                        int k = iCarbon[kk];
                        if (!use[k]) {
                            continue;
                        }
                        if (omit[k] != i) {
                            final double xr = xi - x[k];
                            final double yr = yi - y[k];
                            final double zr = zi - z[k];
                            final double r2 = xr * xr + yr * yr + zr * zr;
                            if (r2 < hpmfCut2) {
                                final double r = sqrt(r2);
                                final double a1 = (r - c1) * w1;
                                final double a2 = (r - c2) * w2;
                                final double a3 = (r - c3) * w3;
                                final double e1 = h1 * exp(-a1 * a1);
                                final double e2 = h2 * exp(-a2 * a2);
                                final double e3 = h3 * exp(-a3 * a3);
                                final double t1t2 = tanhSAi * tanhSA[kk];
                                final double sum = (e1 + e2 + e3);
                                e += sum * t1t2;
                                if (gradient) {
                                    /**
                                     * First part of hydrophobic PMF derivative
                                     * calculation.
                                     */
                                    double de1 = -2.0 * e1 * a1 * w1;
                                    double de2 = -2.0 * e2 * a2 * w2;
                                    double de3 = -2.0 * e3 * a3 * w3;
                                    double dsum = (de1 + de2 + de3) * t1t2 / r;
                                    double dedx = dsum * xr;
                                    double dedy = dsum * yr;
                                    double dedz = dsum * zr;
                                    gX[i] += dedx;
                                    gY[i] += dedy;
                                    gZ[i] += dedz;
                                    gX[k] -= dedx;
                                    gY[k] -= dedy;
                                    gZ[k] -= dedz;
                                    /**
                                     * Chain Rule Term.
                                     */
                                    carbonSASACR[ii] += sum * tanhSA[kk] * dtanhSA[ii];
                                    carbonSASACR[kk] += sum * tanhSA[ii] * dtanhSA[kk];
                                }
                            }
                        }
                    }
                    energy += e;
                }
            }

            @Override
            public void finish() {
                sharedEnergy.addAndGet(energy);
            }
        }

        /**
         * Compute Hydrophobic PMF chain rule term.
         *
         * @since 1.0
         */
        private class CarbonSASACRLoop extends IntegerForLoop {

            private double gX[];
            private double gY[];
            private double gZ[];

            public CarbonSASACRLoop() {
            }

            @Override
            public void start() {
                int threadID = getThreadIndex();
                gX = grad[threadID][0];
                gY = grad[threadID][1];
                gZ = grad[threadID][2];
            }

            @Override
            public void run(int lb, int ub) {
                for (int ii = lb; ii <= ub; ii++) {
                    final int i = iCarbon[ii];
                    if (!use[i]) {
                        continue;
                    }
                    final double carbonSA = carbonSASA[ii];
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    for (int k = 0; k < nAtoms; k++) {
                        if (i != k && use[k]) {
                            final double xr = xi - x[k];
                            final double yr = yi - y[k];
                            final double zr = zi - z[k];
                            final double r2 = xr * xr + yr * yr + zr * zr;
                            double rk = rPMF[k];
                            double rBig = rk + rSurf;
                            if (r2 <= rBig * rBig) {
                                final double r = sqrt(r2);
                                final double rSmall = rk - rCarbon;
                                final double rr = 1.0 / r;
                                final double rr2 = rr * rr;
                                final double part = piSurf * (rBig - r) * (1.0 + rSmall * rr);
                                double t1b = -piSurf * (1.0 + rBig * rSmall * rr2);
                                double t1a = -sasa[ii] / (1.0 / carbonSA - part);
                                double de = t1a * t1b * rr * carbonSASACR[ii];
                                double dedx = de * xr;
                                double dedy = de * yr;
                                double dedz = de * zr;
                                gX[i] += dedx;
                                gY[i] += dedy;
                                gZ[i] += dedz;
                                gX[k] -= dedx;
                                gY[k] -= dedy;
                                gZ[k] -= dedz;
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Some static constants.
     */
    private static final double THIRD = 1.0 / 3.0;
    private static final double PI4_3 = 4.0 / 3.0 * PI;
    private static final double PI_12 = PI / 12.0;

    public enum NonPolar {

        CAV, CAV_DISP, HYDROPHOBIC_PMF, BORN_CAV_DISP, BORN_SOLV, NONE
    }

    private static enum RADII_MAP_TYPE {
        ATOMTYPE, BIOTYPE, NONE;
    }
}
