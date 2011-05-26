/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2011
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.potential;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static java.lang.Math.exp;
import static java.lang.String.format;

import static ffx.numerics.VectorMath.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.numerics.Potential;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;

/**
 * Compute the potential energy and derivatives of an AMOEBA system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ForceFieldEnergy implements Potential, LambdaInterface {

    private static final Logger logger = Logger.getLogger(ForceFieldEnergy.class.getName());
    private final Atom[] atoms;
    private final Crystal crystal;
    private final ParallelTeam parallelTeam;
    private final Bond bonds[];
    private final Angle angles[];
    private final StretchBend stretchBends[];
    private final UreyBradley ureyBradleys[];
    private final OutOfPlaneBend outOfPlaneBends[];
    private final Torsion torsions[];
    private final PiOrbitalTorsion piOrbitalTorsions[];
    private final TorsionTorsion torsionTorsions[];
    private final VanDerWaals vanderWaals;
    private final ParticleMeshEwald particleMeshEwald;
    protected final int nAtoms;
    protected final int nBonds;
    protected final int nAngles;
    protected final int nStretchBends;
    protected final int nUreyBradleys;
    protected final int nOutOfPlaneBends;
    protected final int nTorsions;
    protected final int nPiOrbitalTorsions;
    protected final int nTorsionTorsions;
    protected int nVanDerWaals, nPME, nGK;
    protected final boolean bondTerm;
    protected final boolean angleTerm;
    protected final boolean stretchBendTerm;
    protected final boolean ureyBradleyTerm;
    protected final boolean outOfPlaneBendTerm;
    protected final boolean torsionTerm;
    protected final boolean piOrbitalTorsionTerm;
    protected final boolean torsionTorsionTerm;
    protected final boolean vanderWaalsTerm;
    protected final boolean multipoleTerm;
    protected final boolean polarizationTerm;
    protected final boolean generalizedKirkwoodTerm;
    protected boolean lambdaTerm;
    protected double bondEnergy, bondRMSD;
    protected double angleEnergy, angleRMSD;
    protected double stretchBendEnergy;
    protected double ureyBradleyEnergy;
    protected double outOfPlaneBendEnergy;
    protected double torsionEnergy;
    protected double piOrbitalTorsionEnergy;
    protected double torsionTorsionEnergy;
    protected double totalBondedEnergy;
    protected double vanDerWaalsEnergy;
    protected double permanentMultipoleEnergy;
    protected double polarizationEnergy;
    protected double totalElectrostaticEnergy;
    protected double totalNonBondedEnergy;
    protected double solvationEnergy;
    protected double biasEnergy;
    protected double totalEnergy;
    protected long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    protected long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime;
    protected long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    protected long recursionTime;
    protected long totalTime;
    protected double[] optimizationScaling = null;
    protected double lambda;
    protected boolean lambdaGradient = false;
    protected boolean lambdaElecOff = false;
    protected boolean doCounting = true;
    protected int energyCount;
    protected int λBins = 100;
    protected int FλBins = 401;
    protected int recursionKernel[][];
    protected double dλ = 1.0 / λBins;
    protected double dλ_2 = dλ / 2.0;
    protected double dFλ = 2.0;
    protected double dFλ_2 = dFλ / 2.0;
    protected double minFλ = -(dFλ * FλBins) / 2.0;
    protected double maxFλ = minFλ + FλBins * dFλ;
    protected double dEdλ = 0.0;
    protected double d2Edλ2 = 0.0;
    protected double dUdXdL[] = null;
    protected double gaussianMag = 0.01;
    //protected double gaussianMag = 0.005;
    protected double Fλ[];
    private static final double toSeconds = 0.000000001;

    public ForceFieldEnergy(MolecularAssembly molecularAssembly) {
        parallelTeam = new ParallelTeam();
        logger.info(format(" Constructing Force Field"));
        logger.info(format("\n SMP threads:                        %10d", parallelTeam.getThreadCount()));

        // Get a reference to the sorted atom array.
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;

        ForceField forceField = molecularAssembly.getForceField();
        bondTerm = forceField.getBoolean(ForceFieldBoolean.BONDTERM, true);
        angleTerm = forceField.getBoolean(ForceFieldBoolean.ANGLETERM, true);
        stretchBendTerm = forceField.getBoolean(ForceFieldBoolean.STRBNDTERM, true);
        ureyBradleyTerm = forceField.getBoolean(ForceFieldBoolean.UREYTERM, true);
        outOfPlaneBendTerm = forceField.getBoolean(ForceFieldBoolean.OPBENDTERM, true);
        torsionTerm = forceField.getBoolean(ForceFieldBoolean.TORSIONTERM, true);
        piOrbitalTorsionTerm = forceField.getBoolean(ForceFieldBoolean.PITORSTERM, true);
        torsionTorsionTerm = forceField.getBoolean(ForceFieldBoolean.TORTORTERM, true);
        vanderWaalsTerm = forceField.getBoolean(ForceFieldBoolean.VDWTERM, true);
        multipoleTerm = forceField.getBoolean(ForceFieldBoolean.MPOLETERM, true);
        polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);
        generalizedKirkwoodTerm = forceField.getBoolean(ForceFieldBoolean.GKTERM, false);
        lambdaTerm = forceField.getBoolean(ForceFieldBoolean.LAMBDATERM, false);

        // Define the cutoff lengths.
        double vdwOff = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
        double ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        double buff = 2.0;
        double cutOff2 = 2.0 * (max(vdwOff, ewaldOff) + buff);

        // Determine the unit cell dimensions and Spacegroup
        String spacegroup;
        double a, b, c, alpha, beta, gamma;
        boolean aperiodic;
        try {
            a = forceField.getDouble(ForceFieldDouble.A_AXIS);
            aperiodic = false;
            b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
            c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
            alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
            beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
            gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
            spacegroup = forceField.getString(ForceFieldString.SPACEGROUP, "P1");
        } catch (Exception e) {
            logger.info(" The system will be treated as aperiodic.");
            aperiodic = true;
            spacegroup = "P1";
            /**
             * Search all atom pairs to find the largest pair-wise distance.
             */
            double xr[] = new double[3];
            double maxr = 0.0;
            for (int i = 0; i < nAtoms - 1; i++) {
                double[] xi = atoms[i].getXYZ();
                for (int j = i + 1; j < nAtoms; j++) {
                    double[] xj = atoms[j].getXYZ();
                    diff(xi, xj, xr);
                    double r = r(xr);
                    if (r > maxr) {
                        maxr = r;
                    }
                }
            }
            a = 2.0 * (maxr + ewaldOff);
            b = a;
            c = a;
            alpha = 90.0;
            beta = 90.0;
            gamma = 90.0;
        }
        Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
        unitCell.setAperiodic(aperiodic);

        /**
         * Do we need a ReplicatesCrystal?
         */
        int l = 1;
        int m = 1;
        int n = 1;
        while (unitCell.a * l < cutOff2) {
            l++;
        }
        while (unitCell.b * m < cutOff2) {
            m++;
        }
        while (unitCell.c * n < cutOff2) {
            n++;
        }

        if (l * m * n > 1 && !aperiodic) {
            this.crystal = new ReplicatesCrystal(unitCell, l, m, n);
        } else {
            this.crystal = unitCell;
        }

        logger.info(crystal.toString());

        boolean rigidHydrogens = forceField.getBoolean(ForceFieldBoolean.RIGID_HYDROGENS, false);
        double rigidScale = forceField.getDouble(ForceFieldDouble.RIGID_SCALE, 10.0);

        if (rigidScale <= 1.0) {
            rigidScale = 1.0;
        }

        logger.info("\n Bonded Terms\n");
        if (rigidHydrogens && rigidScale > 1.0) {
            logger.info(format(" Rigid hydrogens:                    %10.2f", rigidScale));
        }

        // Collect, count, pack and sort bonds.
        if (bondTerm) {
            ArrayList<ROLS> bond = molecularAssembly.getBondList();
            nBonds = bond.size();
            bonds = bond.toArray(new Bond[nBonds]);
            Arrays.sort(bonds);
            logger.info(format(" Bonds:                              %10d", nBonds));
        } else {
            nBonds = 0;
            bonds = null;
        }

        // Collect, count, pack and sort angles.
        if (angleTerm) {
            ArrayList<ROLS> angle = molecularAssembly.getAngleList();
            nAngles = angle.size();
            angles = angle.toArray(new Angle[nAngles]);
            Arrays.sort(angles);
            logger.info(format(" Angles:                             %10d", nAngles));
        } else {
            nAngles = 0;
            angles = null;
        }

        // Collect, count, pack and sort stretch-bends.
        if (stretchBendTerm) {
            ArrayList<ROLS> stretchBend = molecularAssembly.getStretchBendList();
            nStretchBends = stretchBend.size();
            stretchBends = stretchBend.toArray(new StretchBend[nStretchBends]);
            Arrays.sort(stretchBends);
            logger.info(format(" Stretch-Bends:                      %10d", nStretchBends));
        } else {
            nStretchBends = 0;
            stretchBends = null;
        }

        // Collect, count, pack and sort Urey-Bradleys.
        if (ureyBradleyTerm) {
            ArrayList<ROLS> ureyBradley = molecularAssembly.getUreyBradleyList();
            nUreyBradleys = ureyBradley.size();
            ureyBradleys = ureyBradley.toArray(new UreyBradley[nUreyBradleys]);
            Arrays.sort(ureyBradleys);
            logger.info(format(" Urey-Bradleys:                      %10d", nUreyBradleys));
        } else {
            nUreyBradleys = 0;
            ureyBradleys = null;
        }

        /**
         * Set a multiplier on the force constants of bonded terms containing
         * hydrogens.
         */
        if (rigidHydrogens) {
            if (bonds != null) {
                for (Bond bond : bonds) {
                    if (bond.containsHydrogen()) {
                        bond.setRigidScale(rigidScale);
                    }
                }
            }
            if (angles != null) {
                for (Angle angle : angles) {
                    if (angle.containsHydrogen()) {
                        angle.setRigidScale(rigidScale);
                    }
                }
            }
            if (stretchBends != null) {
                for (StretchBend stretchBend : stretchBends) {
                    if (stretchBend.containsHydrogen()) {
                        stretchBend.setRigidScale(rigidScale);
                    }
                }
            }
            if (ureyBradleys != null) {
                for (UreyBradley ureyBradley : ureyBradleys) {
                    if (ureyBradley.containsHydrogen()) {
                        ureyBradley.setRigidScale(rigidScale);
                    }
                }
            }
        }

        // Collect, count, pack and sort out-of-plane bends.
        if (outOfPlaneBendTerm) {
            ArrayList<ROLS> outOfPlaneBend = molecularAssembly.getOutOfPlaneBendList();
            nOutOfPlaneBends = outOfPlaneBend.size();
            outOfPlaneBends = outOfPlaneBend.toArray(new OutOfPlaneBend[nOutOfPlaneBends]);
            Arrays.sort(outOfPlaneBends);
            logger.info(format(" Out-of-Plane Bends:                 %10d", nOutOfPlaneBends));
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<ROLS> torsion = molecularAssembly.getTorsionList();
            nTorsions = torsion.size();
            torsions = torsion.toArray(new Torsion[nTorsions]);
            logger.info(format(" Torsions:                           %10d", nTorsions));
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<ROLS> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = piOrbitalTorsion.size();
            piOrbitalTorsions = piOrbitalTorsion.toArray(new PiOrbitalTorsion[nPiOrbitalTorsions]);
            logger.info(format(" Pi-Orbital Torsions:                %10d", nPiOrbitalTorsions));
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<ROLS> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = torsionTorsion.size();
            torsionTorsions = torsionTorsion.toArray(new TorsionTorsion[nTorsionTorsions]);
            logger.info(format(" Torsion-Torsions:                   %10d", nTorsionTorsions));
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }

        logger.info("\n Non-Bonded Terms");

        if (vanderWaalsTerm) {
            vanderWaals = new VanDerWaals(forceField, atoms, crystal, parallelTeam);
        } else {
            vanderWaals = null;
        }

        if (multipoleTerm) {
            particleMeshEwald = new ParticleMeshEwald(forceField, atoms, crystal, parallelTeam,
                                                      vanderWaals.getNeighborLists());
        } else {
            particleMeshEwald = null;
        }

        molecularAssembly.setPotential(this);

    }

    public double energy(boolean gradient, boolean print) {
        bondTime = 0;
        angleTime = 0;
        stretchBendTime = 0;
        ureyBradleyTime = 0;
        outOfPlaneBendTime = 0;
        torsionTime = 0;
        piOrbitalTorsionTime = 0;
        torsionTorsionTime = 0;
        vanDerWaalsTime = 0;
        electrostaticTime = 0;
        recursionTime = 0;
        totalTime = System.nanoTime();

        // Zero out the potential energy of each bonded term.
        bondEnergy = 0.0;
        angleEnergy = 0.0;
        stretchBendEnergy = 0.0;
        ureyBradleyEnergy = 0.0;
        outOfPlaneBendEnergy = 0.0;
        torsionEnergy = 0.0;
        piOrbitalTorsionEnergy = 0.0;
        torsionTorsionEnergy = 0.0;
        totalBondedEnergy = 0.0;

        // Zero out bond and angle RMSDs.
        bondRMSD = 0.0;
        angleRMSD = 0.0;

        // Zero out the potential energy of each non-bonded term.
        vanDerWaalsEnergy = 0.0;
        permanentMultipoleEnergy = 0.0;
        polarizationEnergy = 0.0;
        totalElectrostaticEnergy = 0.0;
        totalNonBondedEnergy = 0.0;

        // Zero out the solvation energy.
        solvationEnergy = 0.0;

        // Zero out the recusion kernel energy.
        biasEnergy = 0.0;

        // Zero out the total potential energy.
        totalEnergy = 0.0;

        // Zero out the Cartesian coordinate gradient for each atom.
        if (gradient) {
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                atom.setXYZGradient(0.0, 0.0, 0.0);
            }
        }

        if (bondTerm) {
            bondTime = System.nanoTime();
            for (int i = 0; i < nBonds; i++) {
                Bond b = bonds[i];
                bondEnergy += b.energy(gradient);
                double value = b.getValue();
                bondRMSD += value * value;
            }
            bondRMSD = sqrt(bondRMSD / bonds.length);
            bondTime = System.nanoTime() - bondTime;
        }


        if (angleTerm) {
            angleTime = System.nanoTime();
            for (int i = 0; i < nAngles; i++) {
                Angle a = angles[i];
                angleEnergy += a.energy(gradient);
                double value = a.getValue();
                angleRMSD += value * value;
            }
            angleRMSD = sqrt(angleRMSD / angles.length);
            angleTime = System.nanoTime() - angleTime;
        }

        if (stretchBendTerm) {
            stretchBendTime = System.nanoTime();
            for (int i = 0; i < nStretchBends; i++) {
                stretchBendEnergy += stretchBends[i].energy(gradient);
            }
            stretchBendTime = System.nanoTime() - stretchBendTime;
        }

        if (ureyBradleyTerm) {
            ureyBradleyTime = System.nanoTime();
            for (int i = 0; i < nUreyBradleys; i++) {
                ureyBradleyEnergy += ureyBradleys[i].energy(gradient);
            }
            ureyBradleyTime = System.nanoTime() - ureyBradleyTime;
        }

        if (outOfPlaneBendTerm) {
            outOfPlaneBendTime = System.nanoTime();
            for (int i = 0; i < nOutOfPlaneBends; i++) {
                outOfPlaneBendEnergy += outOfPlaneBends[i].energy(gradient);
            }
            outOfPlaneBendTime = System.nanoTime() - outOfPlaneBendTime;
        }

        if (torsionTerm) {
            torsionTime = System.nanoTime();
            for (int i = 0; i < nTorsions; i++) {
                torsionEnergy += torsions[i].energy(gradient);
            }
            torsionTime = System.nanoTime() - torsionTime;
        }

        if (piOrbitalTorsionTerm) {
            piOrbitalTorsionTime = System.nanoTime();
            for (int i = 0; i < nPiOrbitalTorsions; i++) {
                piOrbitalTorsionEnergy += piOrbitalTorsions[i].energy(gradient);
            }
            piOrbitalTorsionTime = System.nanoTime() - piOrbitalTorsionTime;
            torsionTorsionTime = System.nanoTime();
        }

        if (torsionTorsionTerm) {
            for (int i = 0; i < nTorsionTorsions; i++) {
                torsionTorsionEnergy += torsionTorsions[i].energy(gradient);
            }
            torsionTorsionTime = System.nanoTime() - torsionTorsionTime;
        }

        if (vanderWaalsTerm) {
            vanDerWaalsTime = System.nanoTime();
            vanDerWaalsEnergy = vanderWaals.energy(gradient, print);
            nVanDerWaals = this.vanderWaals.getInteractions();
            vanDerWaalsTime = System.nanoTime() - vanDerWaalsTime;
        }

        if (multipoleTerm) {
            electrostaticTime = System.nanoTime();
            totalElectrostaticEnergy = particleMeshEwald.energy(gradient, print);
            permanentMultipoleEnergy = particleMeshEwald.getPermanentEnergy();
            polarizationEnergy = particleMeshEwald.getPolarizationEnergy();
            nPME = particleMeshEwald.getInteractions();

            solvationEnergy = particleMeshEwald.getGKEnergy();
            nGK = particleMeshEwald.getGKInteractions();

            electrostaticTime = System.nanoTime() - electrostaticTime;
        }

        if (doCounting) {
            energyCount++;
        }

        if (lambdaGradient) {
            dEdλ = 0.0;
            d2Edλ2 = 0.0;
            if (vanderWaalsTerm) {
                dEdλ = vanderWaals.getdEdLambda();
                d2Edλ2 = vanderWaals.getd2EdLambda2();
            }
            if (multipoleTerm) {
                if (lambdaElecOff) {
                    if (lambda >= 0.5) {
                        dEdλ += particleMeshEwald.getdEdLambda() * 2.0;
                        d2Edλ2 += particleMeshEwald.getd2EdLambda2() * 4.0;
                    }
                } else {
                    dEdλ += particleMeshEwald.getdEdLambda();
                    d2Edλ2 += particleMeshEwald.getd2EdLambda2();
                }
            }

            int λBin = (int) Math.floor(lambda * λBins);
            if (λBin >= λBins) {
                λBin = λBins - 1;
            }

            /**
             * If necessary, allocate more space for dEdλ.
             */
            if (dEdλ > maxFλ) {
                int newFλBins = FλBins;
                while (minFλ + newFλBins * dFλ < dEdλ) {
                    newFλBins += 100;
                }
                int newRecursionKernel[][] = new int[λBins][newFλBins];
                /**
                 * We have added bins above the indeces of the current counts
                 * just copy them into the new array.
                 */
                for (int i = 0; i < λBins; i++) {
                    for (int j = 0; j < FλBins; j++) {
                        newRecursionKernel[i][j] = recursionKernel[i][j];
                    }
                }
                recursionKernel = newRecursionKernel;
                FλBins = newFλBins;
                maxFλ = minFλ + dFλ * FλBins;
            }
            if (dEdλ < minFλ) {
                int offset = 100;
                while (dEdλ < minFλ - offset * dFλ) {
                    offset += 100;
                }
                int newFλBins = FλBins + offset;
                int newRecursionKernel[][] = new int[λBins][newFλBins];
                /**
                 * We have added bins below the current counts,
                 * so their indeces must be increased by: 
                 * offset = newFλBins - FλBins 
                 */
                for (int i = 0; i < λBins; i++) {
                    for (int j = 0; j < FλBins; j++) {
                        newRecursionKernel[i][j + offset] = recursionKernel[i][j];
                    }
                }
                recursionKernel = newRecursionKernel;
                minFλ = minFλ - offset * dFλ;
                FλBins = newFλBins;
            }

            int FλBin = (int) Math.floor((dEdλ - minFλ) / dFλ);
            if (FλBin == FλBins) {
                FλBin = FλBins - 1;
            }
            assert (FλBin < FλBins);
            assert (FλBin >= 0);

            /**
             * Calculate recursion kernel G(λ, dEdλ) and gradient.
             */
            double dGdλ = 0.0;
            double dGdFλ = 0.0;
            double λs2 = 2.0 / λBins * 2.0 / λBins;
            double Fλs2 = dFλ * 2.0 * dFλ * 2.0;
            for (int iλ = -5; iλ <= 5; iλ++) {
                int lcenter = λBin + iλ;
                double deltaλ = lambda - (lcenter * dλ + dλ_2);
                double deltaλ2 = deltaλ * deltaλ;
                // Mirror conditions for recursion kernel counts.
                int lcount = lcenter;
                if (lcount < 0) {
                    lcount = -lcount;
                }
                if (lcount >= λBins) {
                    lcount = λBins - (lcount - λBins) - 1;
                }
                for (int iFλ = -5; iFλ <= 5; iFλ++) {
                    int Fλcenter = FλBin + iFλ;
                    /**
                     * If either of the following edge conditions are true, 
                     * then there are no counts and we continue. 
                     */
                    if (Fλcenter < 0 || Fλcenter >= FλBins) {
                        continue;
                    }
                    double deltaFλ = dEdλ - (minFλ + Fλcenter * dFλ + dFλ_2);
                    double deltaFλ2 = deltaFλ * deltaFλ;
                    int weight = recursionKernel[lcount][Fλcenter];
                    double e = weight * gaussianMag
                               * exp(-deltaλ2 / (2.0 * λs2))
                               * exp(-deltaFλ2 / (2.0 * Fλs2));
                    biasEnergy += e;
                    dGdλ -= deltaλ / λs2 * e;
                    dGdFλ -= deltaFλ / Fλs2 * e;
                }
            }

            /**
             * λ gradient due to recursion kernel G(λ, dEdλ).
             */
            dEdλ += dGdλ + dGdFλ * d2Edλ2;

            /**
             * Atomic gradient due to recursion kernel G(λ, dEdλ).
             */
            for (int i = 0; i < 3 * nAtoms; i++) {
                dUdXdL[i] = 0.0;
            }
            getdEdLambdaGradient(dUdXdL);
            double grad[] = new double[3];
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                atom.getXYZGradient(grad);
                grad[0] += dGdFλ * dUdXdL[i * 3];
                grad[1] += dGdFλ * dUdXdL[i * 3 + 1];
                grad[2] += dGdFλ * dUdXdL[i * 3 + 2];
                atom.setXYZGradient(grad[0], grad[1], grad[2]);
            }

            // Update free energy F(λ) every ~100 steps.
            if (energyCount % 100 == 0 && doCounting) {
                double freeEnergy = 0.0;
                logger.info(" Lambda       Count [   F_L Range    ] <  F_L  > Cumulative");
                for (int iλ = 0; iλ < λBins; iλ++) {
                    int ulFλ = -1;
                    int llFλ = -1;
                    // Find the smallest Fλ bin.
                    for (int jFλ = 0; jFλ < FλBins; jFλ++) {
                        int count = recursionKernel[iλ][jFλ];
                        if (count > 0) {
                            llFλ = jFλ;
                            break;
                        }
                    }
                    // Find the largest Fλ bin.
                    for (int jFλ = FλBins - 1; jFλ >= 0; jFλ--) {
                        int count = recursionKernel[iλ][jFλ];
                        if (count > 0) {
                            ulFλ = jFλ;
                            break;
                        }
                    }

                    int λcount = 0;
                    // The Fλ range that has been sampled for iλ*dλ to (iλ+1)*dλ
                    double lla = minFλ + llFλ * dFλ;
                    double ula = minFλ + ulFλ * dFλ + dFλ;
                    if (ulFλ == -1) {
                        Fλ[iλ] = 0.0;
                        lla = 0.0;
                        ula = 0.0;
                    } else {
                        double sumFλ = 0.0;
                        double partitionFunction = 0.0;
                        for (int jFλ = llFλ; jFλ <= ulFλ; jFλ++) {
                            double a = minFλ + jFλ * dFλ + dFλ_2;
                            double e = exp(gKernel(iλ, jFλ) / (R * 300.0));
                            sumFλ += a * e;
                            partitionFunction += e;
                            λcount += recursionKernel[iλ][jFλ];
                        }
                        Fλ[iλ] = sumFλ / partitionFunction;
                    }
                    freeEnergy += Fλ[iλ] * dλ;
                    logger.info(String.format(" %5.3f..%5.3f %5d [%8.3f %8.3f] <%8.3f> %8.3f",
                                              iλ * dλ, (iλ+1) * dλ, λcount, lla, ula, 
                                              Fλ[iλ], freeEnergy));
                }
                logger.info(" Lambda       Count [   F_L Range    ] <  F_L  > Cumulative\n\n");
            }

            /**
             * Compute the current value of the recursion slave at F(λ)
             * using interpolation.
             */
            for (int i = -1; i < λBins; i++) {

                int iλ0 = i;
                /**
                 * Handle extrapolation from 0.0 to dλ_2.
                 */
                if (iλ0 == -1) {
                    iλ0++;
                }
                /**
                 * Handle extrapolation from 1.0-dλ_2 to 1.0.
                 */
                if (iλ0 == λBins - 1) {
                    iλ0--;
                }
                int iλ1 = iλ0 + 1;
                /**
                 * Find bin centers and values for 
                 * interpolation / extrapolation points. 
                 */
                double λ0 = iλ0 * dλ + dλ_2;
                double λ1 = iλ1 * dλ + dλ_2;
                double Fλ0 = Fλ[iλ0];
                double Fλ1 = Fλ[iλ1];
                double deltaFλ = Fλ1 - Fλ0;
                /**
                 * The integration range is usually λ0 .. λ1.
                 */
                double λll = λ0;
                double λul = λ1;
                /**
                 * Adjust the range for the limiting cases.
                 */
                if (i == -1) {
                    λll = 0.0;
                    λul = dλ_2;
                } else if (i == λBins - 1) {
                    λll = 1.0 - dλ_2;
                    λul = 1.0;
                }
                /**
                 * If the lambda is less than or equal to the upper limit,
                 * this is the final interval. Set the upper limit to λ,
                 * compute the partial derivative and break.
                 */
                boolean done = false;
                if (lambda <= λul) {
                    done = true;
                    λul = lambda;
                }
                /**
                 * Upper limit - lower limit of the integral of the 
                 * extrapolation / interpolation.
                 */
                biasEnergy += Fλ0 * λul + deltaFλ * λul * (0.5 * λul - λ0) / dλ;
                biasEnergy -= Fλ0 * λll + deltaFλ * λll * (0.5 * λll - λ0) / dλ;
                if (done) {
                    /**
                     * Compute the gradient d F(λ) / dλ at λ.  
                     */
                    dEdλ += Fλ0 + (lambda - λ0) * deltaFλ / dλ;
                    /**
                     * Exit the loop.
                     */
                    break;
                }
            }
            logger.info(String.format(" Lambda %8.6f, Bin %d, G %10.4f, dE/dLambda %10.4f",
                                      lambda, λBin, biasEnergy, dEdλ));
            /**
             * Meta-dynamics grid counts (every ~10 steps).
             */
            if (energyCount % 10 == 0 && doCounting) {
                recursionKernel[λBin][FλBin]++;
            }
        }

        totalTime = System.nanoTime() - totalTime;

        totalBondedEnergy = bondEnergy + angleEnergy + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy;
        totalNonBondedEnergy = vanDerWaalsEnergy + totalElectrostaticEnergy;
        totalEnergy = totalBondedEnergy + totalNonBondedEnergy + solvationEnergy + biasEnergy;

        if (print) {
            StringBuilder sb = new StringBuilder("\n");
            if (gradient) {
                sb.append(" Computed Potential Energy and Atomic Coordinate Gradients\n");
            } else {
                sb.append(" Computed Potential Energy\n");
            }
            sb.append(this);
            logger.info(sb.toString());
        }
        return totalEnergy;
    }
    /**
     * Gas constant (in Kcal/mole/Kelvin).
     */
    public static final double R = 1.9872066e-3;

    public void doLambdaCounting(boolean doCounting) {
        this.doCounting = doCounting;
    }

    public double gKernel(int cλ, int cFλ) {
        /**
         * Compute the value of λ and Fλ for the
         * center of the current bin. 
         */
        double vλ = cλ * dλ + dλ_2;
        double vFλ = minFλ + cFλ * dFλ + dFλ_2;
        /**
         * Set the variances for the Gaussian bias.
         */
        double λs2 = 2.0 * dλ * 2.0 * dλ;
        double Fλs2 = 2.0 * dFλ * 2.0 * dFλ;
        double sum = 0.0;
        for (int iλ = -5; iλ <= 5; iλ++) {
            int λcenter = cλ + iλ;
            double deltaλ = vλ - (λcenter * dλ + dλ_2);
            double deltaλ2 = deltaλ * deltaλ;

            // Mirror condition for Lambda counts.
            int lcount = λcenter;
            if (lcount < 0) {
                lcount = -lcount;
            }
            if (lcount >= λBins) {
                lcount = λBins - (lcount - λBins) - 1;
            }

            for (int jFλ = -5; jFλ <= 5; jFλ++) {
                int Fλcenter = cFλ + jFλ;
                /**
                 * For Fλ outside the count matrix the weight is
                 * 0 so we continue.
                 */
                if (Fλcenter < 0 || Fλcenter >= FλBins) {
                    continue;
                }

                double deltaFλ = vFλ - (minFλ + Fλcenter * dFλ + dFλ_2);
                double deltaFλ2 = deltaFλ * deltaFλ;
                int weight = recursionKernel[lcount][Fλcenter];
                if (weight > 0) {
                    double e = weight * gaussianMag * exp(-deltaλ2 / (2.0 * λs2))
                               * exp(-deltaFλ2 / (2.0 * Fλs2));
                    sum += e;
                }
            }
        }
        return sum;
    }

    public double getTotal() {
        return totalEnergy;
    }

    public String getPDBHeaderString() {
        energy(false, false);
        StringBuilder sb = new StringBuilder();
        sb.append("REMARK   3  CALCULATED POTENTIAL ENERGY\n");
        if (bondTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "BOND STRETCHING            : ", bondEnergy, bonds.length));
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "BOND RMSD                  : ", bondRMSD));
        }
        if (angleTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "ANGLE BENDING              : ", angleEnergy, angles.length));
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "ANGLE RMSD                 : ", angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "STRETCH-BEND               : ", stretchBendEnergy, stretchBends.length));
        }
        if (ureyBradleyTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "UREY-BRADLEY               : ", ureyBradleyEnergy, ureyBradleys.length));
        }
        if (outOfPlaneBendTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "OUT-OF-PLANE BEND          : ", outOfPlaneBendEnergy, outOfPlaneBends.length));
        }
        if (torsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "TORSIONAL ANGLE            : ", torsionEnergy, torsions.length));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "PI-ORBITAL TORSION         : ", piOrbitalTorsionEnergy, piOrbitalTorsions.length));
        }
        if (torsionTorsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "TORSION-TORSION            : ", torsionTorsionEnergy, torsionTorsions.length));
        }
        if (vanderWaalsTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "VAN DER WAALS              : ", vanDerWaalsEnergy, nVanDerWaals));
        }
        if (multipoleTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "ATOMIC MULTIPOLES          : ", permanentMultipoleEnergy, nPME));
        }
        if (polarizationTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "POLARIZATION               : ", polarizationEnergy, nPME));
        }
        sb.append(String.format("REMARK   3   %s %g\n",
                                "TOTAL POTENTIAL (KCAL/MOL) : ", totalEnergy));
        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "UNIT CELL POTENTIAL        : ", totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(String.format("REMARK   3   %s %g\n",
                                    "REPLICATES CELL POTENTIAL  : ", totalEnergy * nsymm));
        }
        sb.append("REMARK   3\n");

        return sb.toString();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("\n");
        if (bondTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f (%8.5f)\n",
                                    "Bond Streching    ", bondEnergy, bonds.length,
                                    bondTime * toSeconds, bondRMSD));
        }
        if (angleTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f (%8.5f)\n",
                                    "Angle Bending     ", angleEnergy, angles.length,
                                    angleTime * toSeconds, angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Stretch-Bend      ", stretchBendEnergy,
                                    stretchBends.length, stretchBendTime * toSeconds));
        }
        if (ureyBradleyTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Urey-Bradley      ", ureyBradleyEnergy,
                                    ureyBradleys.length, ureyBradleyTime * toSeconds));
        }
        if (outOfPlaneBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Out-of-Plane Bend ", outOfPlaneBendEnergy,
                                    outOfPlaneBends.length, outOfPlaneBendTime * toSeconds));
        }
        if (torsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsional Angle   ", torsionEnergy, torsions.length,
                                    torsionTime * toSeconds));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Pi-Orbital Torsion", piOrbitalTorsionEnergy,
                                    piOrbitalTorsions.length, piOrbitalTorsionTime * toSeconds));
        }
        if (torsionTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsion-Torsion   ", torsionTorsionEnergy,
                                    torsionTorsions.length, torsionTorsionTime * toSeconds));
        }
        if (vanderWaalsTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Van der Waals     ", vanDerWaalsEnergy,
                                    nVanDerWaals, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm) {
            sb.append(String.format(" %s %16.8f %12d\n",
                                    "Atomic Multipoles ", permanentMultipoleEnergy, nPME));
        }
        if (polarizationTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Polarization      ", polarizationEnergy,
                                    nPME, electrostaticTime * toSeconds));
        }

        if (generalizedKirkwoodTerm) {
            sb.append(String.format(" %s %16.8f %12d\n",
                                    "Solvation         ", solvationEnergy, nGK));
        }

        if (lambdaTerm) {
            sb.append(String.format(" %s %16.8f\n",
                                    "Recursion Kernel  ", biasEnergy));
        }

        sb.append(String.format("\n %s %16.8f  %s %12.3f (sec)\n",
                                "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));
        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(String.format(" %s %16.8f\n", "Unit Cell         ",
                                    totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(String.format(" %s %16.8f\n", "Replicates Cell   ",
                                    totalEnergy * nsymm));
        }

        return sb.toString();
    }

    public Crystal getCrystal() {
        return crystal;
    }

    @Override
    public void setLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
            if (vanderWaalsTerm) {
                vanderWaals.setLambda(lambda);
            }
            if (multipoleTerm) {
                if (lambdaElecOff) {
                    if (lambda >= 0.5) {
                        /**
                         * For lambda above 0.5, turn off electrostatics.
                         */
                        double elecLambda = (lambda - 0.5) * 2.0;
                        if (multipoleTerm) {
                            particleMeshEwald.setLambda(elecLambda);
                            particleMeshEwald.lambdaGradient(true);
                        }
                    } else {
                        particleMeshEwald.setLambda(0.0);
                        particleMeshEwald.lambdaGradient(false);
                    }
                } else {
                    particleMeshEwald.setLambda(lambda);
                }
            }
        } else {
            String message = String.format("Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    @Override
    public void setScaling(double scaling[]) {
        if (scaling != null) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    @Override
    public double energyAndGradient(double x[], double g[]) {
        /**
         * Unscale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] /= optimizationScaling[i];
            }
        }
        setCoordinates(x);
        double e = energy(true, false);
        getGradients(g);
        /**
         * Scale the coordinates and gradients.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
                g[i] /= optimizationScaling[i];
            }
        }
        return e;
    }

    public void getGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZGradient(grad);
            double gx = grad[0];
            double gy = grad[1];
            double gz = grad[2];
            if (gx == Double.NaN || gx == Double.NEGATIVE_INFINITY || gx == Double.POSITIVE_INFINITY
                || gy == Double.NaN || gy == Double.NEGATIVE_INFINITY || gy == Double.POSITIVE_INFINITY
                || gz == Double.NaN || gz == Double.NEGATIVE_INFINITY || gz == Double.POSITIVE_INFINITY) {
                String message = format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                                        a.toString(), gx, gy, gz);
                logger.warning(message);
            }
            g[index++] = gx;
            g[index++] = gy;
            g[index++] = gz;
        }
    }

    private void setCoordinates(double coords[]) {
        assert (coords != null);
        int index = 0;
        for (Atom a : atoms) {
            double x = coords[index++];
            double y = coords[index++];
            double z = coords[index++];
            a.moveTo(x, y, z);
        }
    }

    @Override
    public double[] getCoordinates(double x[]) {
        int n = getNumberOfVariables();
        if (x == null || x.length < n) {
            x = new double[n];
        }
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZ(xyz);
            x[index++] = xyz[0];
            x[index++] = xyz[1];
            x[index++] = xyz[2];
        }
        return x;
    }

    @Override
    public double[] getMass() {
        int n = getNumberOfVariables();
        double mass[] = new double[n];
        int i = 0;
        for (Atom a : atoms) {
            double m = a.getMass();
            mass[i++] = m;
            mass[i++] = m;
            mass[i++] = m;
        }
        return mass;
    }

    @Override
    public int getNumberOfVariables() {
        return nAtoms * 3;
    }

    @Override
    public double getdEdLambda() {
        return dEdλ;
    }

    @Override
    public void lambdaGradient(boolean lambdaGradient) {
        this.lambdaGradient = lambdaGradient;
        if (vanderWaalsTerm) {
            vanderWaals.lambdaGradient(lambdaGradient);
        }
        if (multipoleTerm) {
            particleMeshEwald.lambdaGradient(lambdaGradient);
        }
        energyCount = 0;
        if (recursionKernel == null) {
            /**
             * Start with a 401 bins
             */
            FλBins = 401;
            minFλ = -dFλ * FλBins / 2.0;
            maxFλ = minFλ + FλBins * dFλ;
            recursionKernel = new int[λBins][FλBins];
            Fλ = new double[λBins];
        }
        if (dUdXdL == null) {
            dUdXdL = new double[nAtoms * 3];
        }
    }

    @Override
    public void getdEdLambdaGradient(double gradients[]) {
        if (multipoleTerm) {
            if (lambdaElecOff) {
                if (lambda >= 0.5) {
                    particleMeshEwald.getdEdLambdaGradient(gradients);
                    int nAtoms3 = nAtoms * 3;
                    for (int i = 0; i < nAtoms3; i++) {
                        gradients[i] *= 2.0;
                    }
                }
            } else {
                particleMeshEwald.getdEdLambdaGradient(gradients);
            }
        }
        if (vanderWaalsTerm) {
            vanderWaals.getdEdLambdaGradient(gradients);
        }
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getd2EdLambda2() {
        return d2Edλ2;
    }
}
