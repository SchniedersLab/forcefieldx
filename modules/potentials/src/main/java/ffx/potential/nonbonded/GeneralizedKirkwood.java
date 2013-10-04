/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2013.
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
package ffx.potential.nonbonded;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.Math.tanh;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.DoubleOp;
import edu.rit.pj.reduction.SharedDouble;
import edu.rit.pj.reduction.SharedDoubleArray;
import edu.rit.pj.reduction.SharedInteger;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Torsion;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.VDWType;

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
public class GeneralizedKirkwood {

    private static final Logger logger = Logger.getLogger(GeneralizedKirkwood.class.getName());
    /**
     * Some static constants.
     */
    private static final double third = 1.0 / 3.0;
    private static final double pi43 = 4.0 / 3.0 * PI;
    private static final double pi12 = PI / 12.0;
    /**
     * Permittivity of water at STP.
     */
    private static final double dWater = 78.3;
    /**
     * Kirkwood multipolar reaction field constants.
     */
    private static final double fc = 1.0 * (1.0 - dWater) / (0.0 + 1.0 * dWater);
    private static final double fd = 2.0 * (1.0 - dWater) / (1.0 + 2.0 * dWater);
    private static final double fq = 3.0 * (1.0 - dWater) / (2.0 + 3.0 * dWater);
    /**
     * Empirical constant that controls the GK cross-term.
     */
    private static final double gkc = 2.455;
    /**
     * Empirical scaling of the Bondi radii.
     */
    private static final double bondiScale = 1.03;
    private static final double dispersionOverlapScaleFactor = 0.81;
    private static final double slevy = 1.0;
    private static final double awater = 0.033428;
    private static final double epso = 0.1100;
    private static final double epsh = 0.0135;
    private static final double rmino = 1.7025;
    private static final double rminh = 1.3275;
    private static final double probe = 1.4;
    private boolean use[] = null;
    private final Polarization polarization;
    private final Atom atoms[];
    private final double x[], y[], z[];
    private final double globalMultipole[][];
    private final double inducedDipole[][];
    private final double inducedDipoleCR[][];
    private final double baseRadius[];
    private final double overlapScale[];
    private final double rDisp[];
    private final double born[];
    private final int nAtoms;
    private final ParallelTeam parallelTeam;
    private final Crystal crystal;
    private final BornRadiiRegion bornRadiiRegion;
    private final PermanentGKFieldRegion permanentGKFieldRegion;
    private final InducedGKFieldRegion inducedGKFieldRegion;
    private final GKEnergyRegion gkEnergyRegion;
    private final BornCRRegion bornGradRegion;
    private final HydrophobicPMFRegion hydrophobicPMFRegion;
    private final DispersionRegion dispersionRegion;
    private final CavitationRegion cavitationRegion;
    private final double grad[][][];
    private final double torque[][][];
    private final int neighborLists[][][];
    private final SharedDoubleArray sharedBornGrad;
    protected final SharedDoubleArray sharedGKField[];
    protected final SharedDoubleArray sharedGKFieldCR[];
    private final double cutoff;
    private final double cut2;
    private NonPolar nonPolar = NonPolar.CAV_DISP;

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

        this.parallelTeam = parallelTeam;
        nAtoms = atoms.length;
        polarization = particleMeshEwald.polarization;
        neighborLists = particleMeshEwald.neighborLists;
        this.crystal = crystal;
        this.atoms = atoms;
        x = particleMeshEwald.coordinates[0][0];
        y = particleMeshEwald.coordinates[0][1];
        z = particleMeshEwald.coordinates[0][2];
        globalMultipole = particleMeshEwald.globalMultipole[0];
        grad = particleMeshEwald.getGradient();
        torque = particleMeshEwald.getTorque();

        sharedGKField = new SharedDoubleArray[3];
        sharedGKField[0] = new SharedDoubleArray(nAtoms);
        sharedGKField[1] = new SharedDoubleArray(nAtoms);
        sharedGKField[2] = new SharedDoubleArray(nAtoms);

        sharedGKFieldCR = new SharedDoubleArray[3];
        sharedGKFieldCR[0] = new SharedDoubleArray(nAtoms);
        sharedGKFieldCR[1] = new SharedDoubleArray(nAtoms);
        sharedGKFieldCR[2] = new SharedDoubleArray(nAtoms);

        sharedBornGrad = new SharedDoubleArray(nAtoms);
        baseRadius = new double[nAtoms];
        overlapScale = new double[nAtoms];
        rDisp = new double[nAtoms];
        born = new double[nAtoms];

        use = new boolean[nAtoms];
        Arrays.fill(use, true);

        cutoff = forceField.getDouble(ForceField.ForceFieldDouble.VDW_CUTOFF, 12.0);
        cut2 = cutoff * cutoff;

        for (int i = 0; i < nAtoms; i++) {
            baseRadius[i] = 2.0;
            overlapScale[i] = 0.69;
            int atomicNumber = atoms[i].getAtomicNumber();
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
            baseRadius[i] *= bondiScale;
        }
        int threadCount = parallelTeam.getThreadCount();
        bornRadiiRegion = new BornRadiiRegion(threadCount);

        permanentGKFieldRegion = new PermanentGKFieldRegion(threadCount);
        inducedGKFieldRegion = new InducedGKFieldRegion(threadCount);
        inducedDipole = particleMeshEwald.inducedDipole[0];
        inducedDipoleCR = particleMeshEwald.inducedDipoleCR[0];

        gkEnergyRegion = new GKEnergyRegion(threadCount);
        bornGradRegion = new BornCRRegion(threadCount);

        logger.info(" Continuum Solvation ");
        logger.info(String.format("  Generalized Kirkwood cut-off:        %8.2f (A)", cutoff));
        logger.info(String.format("  Non-Polar Model:                     %s",
                nonPolar.toString().replace('_', '-')));

        if (nonPolar == NonPolar.CAV_DISP) {
            dispersionRegion = new DispersionRegion(threadCount);
            cavitationRegion = new CavitationRegion(threadCount);
            hydrophobicPMFRegion = null;
        } else {
            hydrophobicPMFRegion = new HydrophobicPMFRegion(threadCount);
            dispersionRegion = null;
            cavitationRegion = null;
        }

        logger.info("");
    }

    /**
     * <p>
     * computeBornRadii</p>
     */
    public void computeBornRadii() {
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
                    Atom atom = atoms[i];
                    logger.severe(String.format(" %s\n Born radii %d %8.3f", atom, i, born[i]));
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

        for (int j = 0; j < nAtoms; j++) {
            use[j] = atoms[j].isActive();
        }

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
            gkEnergyRegion.setGradient(gradient);
            parallelTeam.execute(gkEnergyRegion);

            /**
             * Find the nonpolar energy.
             */
            if (nonPolar == NonPolar.HYDROPHOBIC_PMF) {
                hydrophobicPMFRegion.setGradient(gradient);
                parallelTeam.execute(hydrophobicPMFRegion);
            } else {
                dispersionRegion.setGradient(gradient);
                cavitationRegion.setGradient(gradient);
                parallelTeam.execute(dispersionRegion);
                parallelTeam.execute(cavitationRegion);
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
                parallelTeam.execute(bornGradRegion);
            } catch (Exception e) {
                String message = "Fatal exception computing Born radii chain rule term.";
                logger.log(Level.SEVERE, message, e);
            }
        }

        if (print) {
            logger.info(String.format(" Generalized Kirkwood %16.8f",
                    gkEnergyRegion.getEnergy()));
            if (nonPolar == NonPolar.HYDROPHOBIC_PMF) {
                logger.info(String.format(" Hydrophibic PMF      %16.8f", hydrophobicPMFRegion.getEnergy()));
            } else {
                logger.info(String.format(" Dispersion           %16.8f", dispersionRegion.getEnergy()));
                logger.info(String.format(" Cavitation           %16.8f", cavitationRegion.getEnergy()));
            }
        }

        if (nonPolar == NonPolar.HYDROPHOBIC_PMF) {
            return gkEnergyRegion.getEnergy() + hydrophobicPMFRegion.getEnergy();
        } else {
            return gkEnergyRegion.getEnergy() + dispersionRegion.getEnergy() + cavitationRegion.getEnergy();
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

    private enum NonPolar {

        CAV_DISP, HYDROPHOBIC_PMF
    }

    /**
     * Compute Born radii in parallel via the Grycuk method.
     *
     * @since 1.0
     */
    private class BornRadiiRegion extends ParallelRegion {

        private final BornRadiiLoop bornRadiiLoop[];

        public BornRadiiRegion(int nt) {
            bornRadiiLoop = new BornRadiiLoop[nt];
            for (int i = 0; i < nt; i++) {
                bornRadiiLoop[i] = new BornRadiiLoop();
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

        /**
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class BornRadiiLoop extends IntegerForLoop {

            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    born[i] = 0.0;
                    final double baseRi = baseRadius[i];
                    if (!use[i]) {
                        born[i] = baseRi;
                        continue;
                    }
                    if (baseRi > 0.0) {
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        // Integral initialized to the limit of atom alone in solvent.
                        double sum = pi43 / (baseRi * baseRi * baseRi);
                        int list[] = neighborLists[0][i];
                        int npair = list.length;
                        for (int l = 0; l < npair; l++) {
                            int k = list[l];
                            final double baseRk = baseRadius[k];
                            if (i != k && baseRk > 0.0 && use[k]) {
                                final double xr = x[k] - xi;
                                final double yr = y[k] - yi;
                                final double zr = z[k] - zi;
                                final double r2 = crystal.image(xr, yr, zr);
                                if (r2 > cut2) {
                                    continue;
                                }
                                final double r = sqrt(r2);
                                final double scaledRk = baseRk * overlapScale[k];
                                final double scaledRk2 = scaledRk * scaledRk;
                                // Atom i is engulfed atom by atom k.
                                if (baseRi + r < scaledRk) {
                                    final double lower = baseRi;
                                    final double upper = scaledRk - r;
                                    sum = sum + pi43 * (1.0 / (upper * upper * upper)
                                            - 1.0 / (lower * lower * lower));
                                }
                                // Upper integration bound is always the same.
                                final double upper = r + scaledRk;
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
                                final double l2 = lower * lower;
                                final double l4 = l2 * l2;
                                final double lr = lower * r;
                                final double l4r = l4 * r;
                                final double u2 = upper * upper;
                                final double u4 = u2 * u2;
                                final double ur = upper * r;
                                final double u4r = u4 * r;
                                final double term = (3.0 * (r2 - scaledRk2) + 6.0 * u2 - 8.0 * ur) / u4r
                                        - (3.0 * (r2 - scaledRk2) + 6.0 * l2 - 8.0 * lr) / l4r;
                                sum = sum - pi12 * term;
                            }
                        }
                        if (sum <= 0.0) {
                            sum = 0.001;
                        }
                        born[i] = pow(sum / pi43, third);
                        born[i] = 1.0 / born[i];
                    }
                }
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
                        int k = atom2.xyzIndex - 1;
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
                            int k = atom2.xyzIndex - 1;
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
                            int k = atom2.xyzIndex - 1;
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
     * Compute the Generalized Kirkwood reaction field energy.
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
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class PermanentGKFieldLoop extends IntegerForLoop {

            private final double a[][];
            private final double gc[];
            private final double gux[], guy[], guz[];
            private final double gqxx[], gqyy[], gqzz[];
            private final double gqxy[], gqxz[], gqyz[];
            private final double fx_local[];
            private final double fy_local[];
            private final double fz_local[];
            private final double dx_local[];
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
                fx_local = new double[nAtoms];
                fy_local = new double[nAtoms];
                fz_local = new double[nAtoms];
                dx_local = new double[3];
            }

            @Override
            public void start() {
                for (int j = 0; j < nAtoms; j++) {
                    fx_local[j] = 0.0;
                    fy_local[j] = 0.0;
                    fz_local[j] = 0.0;
                }
            }

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double multipolei[] = globalMultipole[i];
                    final double ci = multipolei[t000];
                    final double uxi = multipolei[t100];
                    final double uyi = multipolei[t010];
                    final double uzi = multipolei[t001];
                    final double qxxi = multipolei[t200] / 3.0;
                    final double qxyi = multipolei[t110] / 3.0;
                    final double qxzi = multipolei[t101] / 3.0;
                    final double qyyi = multipolei[t020] / 3.0;
                    final double qyzi = multipolei[t011] / 3.0;
                    final double qzzi = multipolei[t002] / 3.0;
                    final double rbi = born[i];
                    int list[] = neighborLists[0][i];
                    int npair = list.length;
                    for (int l = 0; l < npair; l++) {
                        int k = list[l];
                        if (!use[k]) {
                            continue;
                        }
                        dx_local[0] = x[k] - xi;
                        dx_local[1] = y[k] - yi;
                        dx_local[2] = z[k] - zi;
                        final double r2 = crystal.image(dx_local);
                        if (r2 > cut2) {
                            continue;
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
                         * Multiply the potential auxiliary terms by their
                         * dielectric functions.
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
                         * Scale the self-field by half, such that it sums to
                         * one below.
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
        }
    }

    /**
     * Compute the Generalized Kirkwood reaction field energy.
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
         * Compute Born radii for a range of atoms via the Grycuk method.
         *
         * @since 1.0
         */
        private class InducedGKFieldLoop extends IntegerForLoop {

            private final double a[][];
            private final double gux[], guy[], guz[];
            private final double fx_local[];
            private final double fy_local[];
            private final double fz_local[];
            private final double fxCR_local[];
            private final double fyCR_local[];
            private final double fzCR_local[];
            private final double dx_local[];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public InducedGKFieldLoop() {
                a = new double[3][2];
                gux = new double[5];
                guy = new double[5];
                guz = new double[5];
                fx_local = new double[nAtoms];
                fy_local = new double[nAtoms];
                fz_local = new double[nAtoms];
                fxCR_local = new double[nAtoms];
                fyCR_local = new double[nAtoms];
                fzCR_local = new double[nAtoms];
                dx_local = new double[3];
            }

            @Override
            public void start() {
                for (int j = 0; j < nAtoms; j++) {
                    fx_local[j] = 0.0;
                    fy_local[j] = 0.0;
                    fz_local[j] = 0.0;
                    fxCR_local[j] = 0.0;
                    fyCR_local[j] = 0.0;
                    fzCR_local[j] = 0.0;
                }
            }

            @Override
            public void run(int lb, int ub) {
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    final double uix = inducedDipole[i][0];
                    final double uiy = inducedDipole[i][1];
                    final double uiz = inducedDipole[i][2];
                    final double uixCR = inducedDipoleCR[i][0];
                    final double uiyCR = inducedDipoleCR[i][1];
                    final double uizCR = inducedDipoleCR[i][2];
                    final double rbi = born[i];
                    int list[] = neighborLists[0][i];
                    int nPair = list.length;
                    for (int l = 0; l < nPair; l++) {
                        int k = list[l];
                        if (!use[k]) {
                            continue;
                        }
                        dx_local[0] = x[k] - xi;
                        dx_local[1] = y[k] - yi;
                        dx_local[2] = z[k] - zi;
                        final double r2 = crystal.image(dx_local);
                        if (r2 > cut2) {
                            continue;
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
                         * Multiply the potential auxiliary terms by their
                         * dielectric functions.
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
                         * Scale the self-field by half, such that it sums to
                         * one below.
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

            private boolean gradient = false;
            private int count;
            private double gkEnergy;
            private final double a[][];
            private final double b[][];
            private final double gc[];
            private final double gux[], guy[], guz[];
            private final double gqxx[], gqyy[], gqzz[];
            private final double gqxy[], gqxz[], gqyz[];
            private double gX[];
            private double gY[];
            private double gZ[];
            private double tX[];
            private double tY[];
            private double tZ[];
            private final double gb_local[];
            private final double gbi_local[];
            private final double dx_local[];
            private double ci, uxi, uyi, uzi, qxxi, qxyi, qxzi, qyyi, qyzi, qzzi;
            private double ck, uxk, uyk, uzk, qxxk, qxyk, qxzk, qyyk, qyzk, qzzk;
            private double dxi, dyi, dzi, pxi, pyi, pzi, sxi, syi, szi;
            private double dxk, dyk, dzk, pxk, pyk, pzk, sxk, syk, szk;
            private double xr, yr, zr, xr2, yr2, zr2, r2, rbi, rbk;
            private int i, k;
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
                gb_local = new double[nAtoms];
                gbi_local = new double[nAtoms];
                dx_local = new double[3];
            }

            public void setGradient(boolean gradient) {
                this.gradient = gradient;
            }

            @Override
            public void start() {
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
                    Arrays.fill(gb_local, 0.0);
                    Arrays.fill(gbi_local, 0.0);
                }
            }

            @Override
            public void run(int lb, int ub) {
                for (i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
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
                        k = list[l];
                        if (!use[k]) {
                            continue;
                        }
                        dx_local[0] = x[k] - xi;
                        dx_local[1] = y[k] - yi;
                        dx_local[2] = z[k] - zi;
                        r2 = crystal.image(dx_local);
                        if (r2 > cut2) {
                            continue;
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
                             * Born radii derivatives of reaction potential
                             * auxiliary terms.
                             */
                            b[0][0] = dgfdr * a[1][0];
                            b[1][0] = dgfdr * a[2][0];
                            b[2][0] = dgfdr * a[3][0];
                            b[3][0] = dgfdr * a[4][0];
                            b[4][0] = dgfdr * a[5][0];
                            /**
                             * Born radii gradients of reaction potential
                             * gradient auxiliary terms.
                             */
                            b[0][1] = b[1][0] - expcr * a[1][0] - expc * b[1][0];
                            b[1][1] = b[2][0] - expcr * a[2][0] - expc * b[2][0];
                            b[2][1] = b[3][0] - expcr * a[3][0] - expc * b[3][0];
                            b[3][1] = b[4][0] - expcr * a[4][0] - expc * b[4][0];
                            /**
                             * Born radii derivatives of the 2nd reaction
                             * potential gradient auxiliary terms.
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
                         * Multiply the potential auxiliary terms by their
                         * dielectric functions.
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
                         * Compute the GK tensors required to compute the
                         * energy.
                         */
                        energyTensors();
                        double gkPairEnergy = energy();
                        gkEnergy += gkPairEnergy;
                        count++;
                        if (gradient) {
                            /**
                             * Compute the additional GK tensors required to
                             * compute the energy gradient.
                             */
                            gradientTensors();
                            permanentEnergyGradient();
                            if (polarization != Polarization.NONE) {
                                polarizationEnergyGradient();
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
                    /**
                     * Reduce the force and torque contributions computed by the
                     * current thread into the shared arrays.
                     */
                    sharedBornGrad.reduce(gb_local, DoubleOp.SUM);
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

            private double energy() {
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

            private void permanentEnergyGradient() {
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

                gX[i] -= dedx;
                gY[i] -= dedy;
                gZ[i] -= dedz;
                gb_local[i] += drbi;

                gX[k] += dedx;
                gY[k] += dedy;
                gZ[k] += dedz;
                gb_local[k] += drbk;
                permanentEnergyTorque();

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

            private void permanentEnergyTorque() {
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
                tX[i] += tix;
                tY[i] += tiy;
                tZ[i] += tiz;
                tX[k] += tkx;
                tY[k] += tky;
                tZ[k] += tkz;
            }

            private void polarizationEnergyGradient() {
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
                    gX[i] -= dpdx;
                    gY[i] -= dpdy;
                    gZ[i] -= dpdz;
                    gb_local[i] += dbi;

                    gX[k] += dpdx;
                    gY[k] += dpdy;
                    gZ[k] += dpdz;
                    gb_local[k] += dbk;
                }
                polarizationEnergyTorque();
            }

            private void polarizationEnergyTorque() {
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
                tX[i] += tix;
                tY[i] += tiy;
                tZ[i] += tiz;
                tX[k] += tkx;
                tY[k] += tky;
                tZ[k] += tkz;
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

        public BornCRRegion(int nt) {
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

            private final double factor = -pow(PI, third) * pow(6.0, (2.0 * third)) / 9.0;
            private final double dx_local[];
            private double gX[];
            private double gY[];
            private double gZ[];
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public BornCRLoop() {
                dx_local = new double[3];
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
                for (int i = lb; i <= ub; i++) {
                    if (!use[i]) {
                        continue;
                    }
                    final double ri = baseRadius[i];
                    if (ri > 0.0) {
                        final double xi = x[i];
                        final double yi = y[i];
                        final double zi = z[i];
                        final double rbi = born[i];
                        double term = pi43 / (rbi * rbi * rbi);
                        term = factor / pow(term, (4.0 * third));
                        int list[] = neighborLists[0][i];
                        int nPair = list.length;
                        for (int l = 0; l < nPair; l++) {
                            int k = list[l];
                            if (!use[k]) {
                                continue;
                            }
                            final double rk = baseRadius[k];
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
                                final double sk = rk * overlapScale[k];
                                final double sk2 = sk * sk;
                                final double r = sqrt(r2);
                                double de = 0.0;
                                if (ri + r < sk) {
                                    final double uik = sk - r;
                                    de = -4.0 * PI / pow(uik, 4);
                                }
                                if (ri + r < sk) {
                                    final double lik = sk - r;
                                    double lik4 = pow(lik, 4);
                                    de = de + 0.25 * PI * (sk2 - 4.0 * sk * r + 17.0 * r2) / (r2 * lik4);
                                } else if (r < ri + sk) {
                                    double lik = ri;
                                    double lik4 = pow(lik, 4);
                                    de = de + 0.25 * PI * (2.0 * ri * ri - sk2 - r2) / (r2 * lik4);
                                } else {
                                    double lik = r - sk;
                                    double lik4 = pow(lik, 4);
                                    de = de + 0.25 * PI * (sk2 - 4.0 * sk * r + r2) / (r2 * lik4);
                                }
                                double uik = r + sk;
                                double uik4 = pow(uik, 4);
                                de = de - 0.25 * PI * (sk2 + 4.0 * sk * r + r2) / (r2 * uik4);
                                double dbr = term * de / r;
                                de = dbr * sharedBornGrad.get(i);
                                /**
                                 * Increment the overall derivatives.
                                 */
                                final double dedx = de * xr;
                                final double dedy = de * yr;
                                final double dedz = de * zr;
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
     * Compute Dispersion energy in parallel via pairwise descreening.
     *
     * @since 1.0
     */
    private class DispersionRegion extends ParallelRegion {

        private final DispersionLoop dispersionLoop[];
        private final SharedDouble sharedDispersion;
        private boolean gradient = false;
        private double[] cdisp = null;

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
                    double sqEpsoEpsi = sqrt(epso) + sqrt(epsi);
                    double sqEpshEpsi = sqrt(epsh) + sqrt(epsi);
                    double emixo = 4.0 * epso * epsi / (pow(sqEpsoEpsi, 2));
                    double rmixo = 2.0 * (pow(rmino, 3) + pow(rmini, 3)) / (pow(rmino, 2) + pow(rmini, 2));
                    double rmixo3 = pow(rmixo, 3);
                    double rmixo7 = pow(rmixo, 7);
                    double ao = emixo * rmixo7;
                    double emixh = 4.0 * epsh * epsi / (pow(sqEpshEpsi, 2));
                    double rmixh = 2.0 * (pow(rminh, 3) + pow(rmini, 3)) / (pow(rminh, 2) + pow(rmini, 2));
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
                cdisp[i] = slevy * awater * cdisp[i];
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
            private double edisp;
            private final double dx_local[];
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
                    VDWType type = atoms[i].getVDWType();
                    double epsi = type.wellDepth;
                    double rmini = type.radius / 2.0;
                    double emixo = (4.0 * epso * epsi) / (pow(sqrt(epso) + sqrt(epsi), 2));
                    double rmixo = 2.0 * (pow(rmino, 3) + pow(rmini, 3)) / (pow(rmino, 2) + pow(rmini, 2));
                    double rmixo7 = pow(rmixo, 7);
                    double ao = emixo * rmixo7;
                    double emixh = 4.0 * epsh * epsi / (pow(sqrt(epsh) + sqrt(epsi), 2));
                    double rmixh = 2.0 * (pow(rminh, 3) + pow(rmini, 3)) / (pow(rminh, 2) + pow(rmini, 2));
                    double rmixh7 = pow(rmixh, 7);
                    double ah = emixh * rmixh7;
                    final double ri = rDisp[i];
                    final double xi = x[i];
                    final double yi = y[i];
                    final double zi = z[i];
                    // Integral initialized to the limit of atom alone in solvent.
                    double sum = 0.0;

                    int list[] = neighborLists[0][i];
                    int npair = list.length;
                    for (int l = 0; l < npair; l++) {
                        int k = list[l];
                        final double rk = rDisp[k];
                        if (i != k && rk > 0.0 && use[k]) {
                            dx_local[0] = xi - x[k];
                            dx_local[1] = yi - y[k];
                            dx_local[2] = zi - z[k];
                            double r2 = crystal.image(dx_local);
                            if (r2 > cut2) {
                                continue;
                            }
                            double xr = dx_local[0];
                            double yr = dx_local[1];
                            double zr = dx_local[2];
                            final double r = sqrt(r2);
                            final double r3 = r * r2;
                            double sk = rk * dispersionOverlapScaleFactor;
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
                                //   increment the individual dispersion gradient components
                                if (gradient) {
                                    de = -de / r * slevy * awater;
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
                    // increment the overall dispersion energy component
                    double e = cdisp[i] - slevy * awater * sum;
                    edisp += e;
                }
            }
        }
    }

    /**
     * Compute Cavitation energy in parallel.
     *
     * @since 1.0
     */
    private class CavitationRegion extends ParallelRegion {

        private final CavitationLoop cavitationLoop[];
        private final SharedDouble sharedCavitation;
        private boolean gradient = false;
        private final int maxarc = 1000;

        public CavitationRegion(int nt) {
            cavitationLoop = new CavitationLoop[nt];
            for (int i = 0; i < nt; i++) {
                cavitationLoop[i] = new CavitationLoop();
            }
            sharedCavitation = new SharedDouble();
        }

        public double getEnergy() {
            return sharedCavitation.get();
        }

        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }

        @Override
        public void start() {
            sharedCavitation.set(0.0);
        }

        @Override
        public void run() {
            try {
                execute(0, nAtoms - 1, cavitationLoop[getThreadIndex()]);
            } catch (Exception e) {
                String message = "Fatal exception computing Cavitation energy in thread " + getThreadIndex() + "\n";
                logger.log(Level.SEVERE, message, e);
            }
        }

        /**
         * Compute Cavitation energy for a range of atoms.
         *
         * @since 1.0
         */
        private class CavitationLoop extends IntegerForLoop {

            private double thec = 0;
            private IndexedDouble gr[];
            private IndexedDouble arci[];
            private double area[];
            private double darea[][];
            private double r[];
            private boolean skip[];
            private boolean omit[];
            private double xc1[];
            private double yc1[];
            private double zc1[];
            private double dsq1[];
            private double bsq1[];
            private double b1[];
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
            private int intag1[];
            private int lt[];
            private int i;
            private int j;
            private int ib;
            private double ecav;
            // Set pi multiples, overlap criterion and tolerances.
            private final static double pix2 = 2.0 * PI;
            private final static double pix4 = 4.0 * PI;
            private final static double pid2 = PI / 2.0;
            private final static double eps = 1.0e-8;
            private final static double delta = 1.0e-8;
            private final static double delta2 = delta * delta;
            private final static double rmove = 1.0e-8;
            private final static double surfaceTension = 0.08;
            // Extra padding to avert cache interference.
            private long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
            private long pad8, pad9, pada, padb, padc, padd, pade, padf;

            public CavitationLoop() {
                area = new double[nAtoms];
                r = new double[nAtoms];
                skip = new boolean[nAtoms];
                darea = new double[3][];
                allocateMemory(maxarc);
            }

            private void allocateMemory(int maxarc) {
                gr = new IndexedDouble[maxarc];
                arci = new IndexedDouble[maxarc];
                arcf = new double[maxarc];
                risq = new double[maxarc];
                ri = new double[maxarc];
                dsq1 = new double[maxarc];
                bsq1 = new double[maxarc];
                dsq = new double[maxarc];
                bsq = new double[maxarc];
                intag = new int[maxarc];
                intag1 = new int[maxarc];
                lt = new int[maxarc];
                kent = new int[maxarc];
                kout = new int[maxarc];
                ider = new double[maxarc];
                sign_yder = new double[maxarc];
                xc1 = new double[maxarc];
                yc1 = new double[maxarc];
                zc1 = new double[maxarc];
                b1 = new double[maxarc];
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
                int threadID = getThreadIndex();
                darea[0] = grad[threadID][0];
                darea[1] = grad[threadID][1];
                darea[2] = grad[threadID][2];
                ecav = 0;
                // Initialize the area.
                Arrays.fill(area, 0.0);
                Arrays.fill(ider, 0);
                Arrays.fill(sign_yder, 0);
                Arrays.fill(skip, true);
                /**
                 * Set the sphere radii.
                 */
                for (int i = 0; i < nAtoms; i++) {
                    r[i] = rDisp[i];
                    if (r[i] != 0.0) {
                        r[i] = r[i] + probe;
                    }
                }
                /**
                 * Set the "skip" array to exclude all inactive atoms that do
                 * not overlap any of the current active atoms.
                 */
                for (int i = 0; i < nAtoms; i++) {
                    if (use[i]) {
                        double xr = x[i];
                        double yr = y[i];
                        double zr = z[i];
                        double rr = r[i];
                        for (int k = 0; k < nAtoms; k++) {
                            double rrk = rr + r[k];
                            double rplus = rrk * rrk;
                            double dx = x[k] - xr;
                            double dy = y[k] - yr;
                            double dz = z[k] - zr;
                            double ccsq = dx * dx + dy * dy + dz * dz;
                            if (ccsq <= rplus) {
                                skip[k] = false;
                            }
                        }
                    }
                }
            }

            @Override
            public void finish() {
                sharedCavitation.addAndGet(ecav);
            }

            @Override
            public void run(int lb, int ub) {
                /**
                 * Compute the area and derivatives of current "ir" sphere
                 */
                for (int ir = lb; ir <= ub; ir++) {
                    if (!use[ir] || skip[ir]) {
                        continue;
                    }
                    double xr = x[ir];
                    double yr = y[ir];
                    double zr = z[ir];
                    double rr = r[ir];
                    double rrx2 = 2.0 * rr;
                    double rrsq = rr * rr;
                    double wght = surfaceTension;
                    boolean moved = false;
                    surface(xr, yr, zr, rr, rrx2, rrsq, wght, moved, ir);
                    if (area[ir] < 0.0) {
                        xr = xr + rmove;
                        yr = yr + rmove;
                        zr = zr + rmove;
                        moved = true;
                        surface(xr, yr, zr, rr, rrx2, rrsq, wght, moved, ir);
                        if (area[ir] < 0.0) {
                            logger.warning(String.format(" Negative surface area set to 0 for atom %d.", ir));
                            area[ir] = 0.0;
                            continue;
                        }
                    }
                    area[ir] *= rrsq * wght;
                    ecav += area[ir];
                }
                /**
                 * Zero out the area derivatives for the inactive atoms.
                 */
                for (int i = 0; i < nAtoms; i++) {
                    if (!use[i]) {
                        darea[0][i] = 0.0;
                        darea[1][i] = 0.0;
                        darea[2][i] = 0.0;
                    }
                }
            }

            public void surface(double xr, double yr, double zr, double rr, double rrx2,
                    double rrsq, double wght, boolean moved, int ir) {
                int io = 0;
                int jb = 0;
                ib = 0;
                double arclen = 0.0;
                double exang = 0.0;
                /**
                 * Test each sphere to see if it overlaps the "ir" sphere. TODO:
                 * Use the neighbor list & loop over symmetry mates!
                 */
                for (int i = 0; i < nAtoms; i++) {
                    if (i == ir || !use[i]) {
                        continue;
                    }
                    double rplus = rr + r[i];
                    double tx = x[i] - xr;
                    double ty = y[i] - yr;
                    double tz = z[i] - zr;
                    if (abs(tx) >= rplus || abs(ty) >= rplus || abs(tz) >= rplus) {
                        continue;
                    }
                    /**
                     * Check for overlap of spheres by testing center to center
                     * distance against sum and difference of radii.
                     */
                    double xysq = tx * tx + ty * ty;
                    if (xysq < delta2) {
                        tx = delta;
                        ty = 0.0;
                        xysq = delta2;
                    }
                    double ccsq = xysq + tz * tz;
                    double cc = sqrt(ccsq);
                    if (rplus - cc <= delta) {
                        continue;
                    }
                    double rminus = rr - r[i];
                    /**
                     * Check for a completely buried "ir" sphere.
                     */
                    if (cc - abs(rminus) <= delta) {
                        if (rminus <= 0.0) {
                            // SA for this atom is zero.
                            area[ir] = 0.0;
                            return;
                        }
                        continue;
                    }
                    /**
                     * Calculate overlap parameters between "i" and "ir" sphere.
                     */
                    io = io + 1;
                    int io1 = io - 1;
                    xc1[io1] = tx;
                    yc1[io1] = ty;
                    zc1[io1] = tz;
                    dsq1[io1] = xysq;
                    bsq1[io1] = ccsq;
                    b1[io1] = cc;
                    gr[io1] = new IndexedDouble((ccsq + rplus * rminus) / (rrx2 * b1[io1]), io1);
                    intag1[io1] = i;
                    if (io > maxarc) {
                        logger.severe(String.
                                format(" Increase the value of MAXARC to (%d).", io));
                    }
                }
                // Case where no other spheres overlap the current sphere.
                if (io == 0) {
                    area[ir] = pix4;
                    return;
                }
                // Case where only one sphere overlaps the current sphere.
                if (io == 0) {
                    int k = 0;
                    double txk = xc1[0];
                    double tyk = yc1[0];
                    double tzk = zc1[0];
                    double bsqk = bsq1[0];
                    double bk = b1[0];
                    intag[0] = intag1[0];
                    double arcsum = pix2;
                    ib = ib + 1;
                    arclen += gr[k].value * arcsum;
                    if (!moved) {
                        int in = intag[k];
                        double t1 = arcsum * rrsq * (bsqk - rrsq + r[in] * r[in])
                                / (rrx2 * bsqk * bk);
                        darea[0][ir] -= txk * t1 * wght;
                        darea[1][ir] -= tyk * t1 * wght;
                        darea[2][ir] -= tzk * t1 * wght;
                        darea[0][in] += txk * t1 * wght;
                        darea[1][in] += tyk * t1 * wght;
                        darea[2][in] += tzk * t1 * wght;
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
                Arrays.sort(gr, 0, io);
                for (int i = 0; i < io; i++) {
                    int k = gr[i].key;
                    intag[i] = intag1[k];
                    xc[i] = xc1[k];
                    yc[i] = yc1[k];
                    zc[i] = zc1[k];
                    dsq[i] = dsq1[k];
                    b[i] = b1[k];
                    bsq[i] = bsq1[k];
                    omit[i] = false;
                }
                /**
                 * Radius of the each circle on the surface of the "ir" sphere.
                 */
                for (int i = 0; i < io; i++) {
                    double gi = gr[i].value * rr;
                    bg[i] = b[i] * gi;
                    risq[i] = rrsq - gi * gi;
                    ri[i] = sqrt(risq[i]);
                    ther[i] = pid2 - asin(min(1.0, max(-1.0, gr[i].value)));
                }
                /**
                 * Find boundary of inaccessible area on "ir" sphere.
                 */
                for (int k = 0; k < io - 1; k++) {
                    if (omit[k]) {
                        continue;
                    }
                    double txk = xc[k];
                    double tyk = yc[k];
                    double tzk = zc[k];
                    double bk = b[k];
                    double therk = ther[k];
                    for (j = k + 1; j < io; j++) {
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
                for (int k = 0; k < io; k++) {
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
                    double gk = gr[k].value * rr;
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
                    for (int l = 0; l < io; l++) {
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
                            thec = (rrsq * uzl - gk * bg[l])
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
                                    / (b[l] * rr)));
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
                                logger.severe(String.
                                        format(" Increase value of MAXARC %d.", narc));
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
                        arclen += gr[k].value * arcsum;
                        if (!moved) {
                            int in = intag[k];
                            t1 = arcsum * rrsq * (bsqk - rrsq + r[in] * r[in])
                                    / (rrx2 * bsqk * bk);
                            //logger.info(String.format(" %d %d t1 %16.8f %16.8f %16.8f %16.8f",
                            //        ir, in, txk, tyk, tzk, t1));
                            darea[0][ir] -= txk * t1 * wght;
                            darea[1][ir] -= tyk * t1 * wght;
                            darea[2][ir] -= tzk * t1 * wght;
                            darea[0][in] += txk * t1 * wght;
                            darea[1][in] += tyk * t1 * wght;
                            darea[2][in] += tzk * t1 * wght;
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
                                logger.severe(String.
                                        format("Increase the value of MAXARC (%d).", jb));
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
                    for (int l = 0; l <= io; l++) {
                        if (ider[l] == 0) {
                            continue;
                        }
                        double rcn = ider[l] * rrsq;
                        ider[l] = 0;
                        double uzl = uz[l];
                        double gl = gr[l].value * rr;
                        double bgl = bg[l];
                        double bsql = bsq[l];
                        double risql = risq[l];
                        double wxlsq = bsql - uzl * uzl;
                        double wxl = sqrt(wxlsq);
                        double p = bgl - gk * uzl;
                        double v = risqk * wxlsq - p * p;
                        v = max(eps, v);
                        v = sqrt(v);
                        t1 = rr * (gk * (bgl - bsql) + uzl * (bgl - rrsq))
                                / (v * risql * bsql);
                        double deal = -wxl * t1;
                        double decl = -uzl * t1 - rr / v;
                        double dtkal = (wxlsq - p) / (wxl * v);
                        double dtkcl = (uzl - gk) / v;
                        double s = gk * b[l] - gl * uzl;
                        t1 = 2.0 * gk - uzl;
                        double t2 = rrsq - bgl;
                        double dtlal = -(risql * wxlsq * b[l] * t1 - s * (wxlsq * t2 + risql * bsql))
                                / (risql * wxl * bsql * v);
                        double dtlcl = -(risql * b[l] * (uzl * t1 - bgl) - uzl * t2 * s)
                                / (risql * bsql * v);
                        double gaca = rcn * (deal - (gk * dtkal - gl * dtlal) / rr) / wxl;
                        double gacb = (gk - uzl * gl / b[l]) * sign_yder[l] * rr / wxlsq;
                        sign_yder[l] = 0;
                        if (!moved) {
                            double faca = ux[l] * gaca - uy[l] * gacb;
                            double facb = uy[l] * gaca + ux[l] * gacb;
                            double facc = rcn * (decl - (gk * dtkcl - gl * dtlcl) / rr);
                            double dax = axx * faca - ayx * facb + azx * facc;
                            double day = axy * faca + ayy * facb + azy * facc;
                            double daz = azz * facc - axz * faca;
                            int in = intag[l];
                            //logger.info(String.format(" %d %d dax,y,z %16.8f %16.8f %16.8f", ir, in, dax,day,daz));
                            darea[0][ir] += dax * wght;
                            darea[1][ir] += day * wght;
                            darea[2][ir] += daz * wght;
                            darea[0][in] -= dax * wght;
                            darea[1][in] -= day * wght;
                            darea[2][in] -= daz * wght;
                        }

                    }
                    arclen = arclen + gr[k].value * arcsum;
                    if (!moved) {
                        int in = intag[k];
                        t1 = arcsum * rrsq * (bsqk - rrsq + r[in] * r[in])
                                / (rrx2 * bsqk * bk);
                        //logger.info(String.format(" %d %d t1 %16.8f %16.8f %16.8f %16.8f",
                        //        ir, in, txk, tyk, tzk, t1));
                        darea[0][ir] -= txk * t1 * wght;
                        darea[1][ir] -= tyk * t1 * wght;
                        darea[2][ir] -= tzk * t1 * wght;
                        darea[0][in] += txk * t1 * wght;
                        darea[1][in] += tyk * t1 * wght;
                        darea[2][in] += tzk * t1 * wght;
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
                if (moved) {
                    logger.warning(String.format(" Connectivity error at atom %d.", ir));
                } else {
                    moved = true;
                    xr += rmove;
                    yr += rmove;
                    zr += rmove;
                    surface(xr, yr, zr, rr, rrx2, rrsq, wght, moved, ir);
                }
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
                            ib = ib + 1;
                            if (j == jb) {
                                area[ir] = ib * 2 * PI + exang + arclen;
                                area[ir] = area[ir] % (4 * PI);
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
        }
    }
}
