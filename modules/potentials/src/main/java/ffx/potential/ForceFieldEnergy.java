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
    protected double recursionKernelEnergy;
    protected double totalEnergy;
    protected long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    protected long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime;
    protected long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    protected long recursionTime;
    protected long totalTime;
    protected double[] optimizationScaling = null;
    protected double lambda;
    protected int energyCount;
    protected int lambdaBins = 100;
    protected double λBinWidth = 1.0 / lambdaBins;
    protected double λBinHalfWidth = λBinWidth / 2.0;
    protected double mindEdλ = -500.0, maxdEdλ = 500.0;
    protected double dEdλSpan = maxdEdλ - mindEdλ;
    protected double dEdλWidth = 2.0;
    protected int dEdλBins = (int) Math.floor(dEdλSpan / dEdλWidth);
    protected int recursionKernel[][];
    protected double dEdλ = 0.0;
    protected double d2Edλ2 = 0.0;
    protected double gaussianMag = 0.002;
    protected double dAdλ[];
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
            //vanderWaals = new CellCellVanDerWaals(forceField, atoms, crystal, parallelTeam);
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
        recursionKernelEnergy = 0.0;

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

        energyCount++;

        if (lambdaTerm) {

            dEdλ = vanderWaals.getdEdLambda();
            d2Edλ2 = vanderWaals.getd2EdLambda2();
            
            logger.info(String.format("Lambda %20.8f, dE/dLambda %20.8f",lambda, dEdλ));

            int lambdaBin = (int) Math.floor(lambda * lambdaBins);
            if (lambda == 1.0) {
                lambdaBin = lambdaBins - 1;
            }

            int dUdLBin = (int) Math.floor((dEdλ - mindEdλ) / dEdλWidth);
            if (dEdλ == maxdEdλ) {
                dUdLBin = dUdLBin - 1;
            }
            
            /**
             * Allocate more space for dEdλ.
             */
            if (dUdLBin >= dEdλBins) {
                double newMaxdEdλ = maxdEdλ + 100.0;
                int newdEdλBins = (int) Math.floor((newMaxdEdλ - mindEdλ) / dEdλWidth);                
                int newRecursionKernel[][] = new int[lambdaBins][newdEdλBins];
                for (int i=0; i<lambdaBins; i++) {
                    for (int j=0; j<dEdλBins; j++) {
                        newRecursionKernel[i][j] = recursionKernel[i][j];
                    }
                }
                recursionKernel = newRecursionKernel;
                maxdEdλ = newMaxdEdλ;
                dEdλBins = newdEdλBins;
            }

            if (dUdLBin < 0) {
                double newMindEdλ = mindEdλ - 100.0;
                int newdEdλBins = (int) Math.floor((maxdEdλ - newMindEdλ) / dEdλWidth);                
                int newRecursionKernel[][] = new int[lambdaBins][newdEdλBins];
                for (int i=0; i<lambdaBins; i++) {
                    for (int j=0; j<dEdλBins; j++) {
                        newRecursionKernel[i][j] = recursionKernel[i][j];
                    }
                }
                recursionKernel = newRecursionKernel;
                mindEdλ = newMindEdλ;
                dEdλBins = newdEdλBins;
            }
            
            /**
             * Calculate recursion kernel G(λ, dEdλ) and gradient.
             */
            double dgdL = 0.0;
            double dgdEdL = 0.0;
            double ls2 = 2.0 / lambdaBins * 2.0 / lambdaBins;
            double dEdLs2 = dEdλWidth * 2.0 * dEdλWidth * 2.0;
            for (int l = -5; l <= 5; l++) {
                int lcenter = lambdaBin + 1;
                double lv = lcenter / lambdaBins + 0.5 / lambdaBins;
                double lv2 = (lambda - lv) * (lambda - lv);

                // Mirror conditions for recursion kernel counts.
                int lcount = lcenter;
                if (lcount < 0) {
                    lcount = -lcount;
                }
                if (lcount >= lambdaBins) {
                    lcount = lambdaBins - (lcount - lambdaBins) - 1;
                }

                for (int dl = -5; dl <= 5; dl++) {
                    int dlcenter = dUdLBin + dl;
                    double dlv = dlcenter * dEdλWidth + dEdλWidth / 2.0;
                    double dlv2 = (dEdλ - dlv) * (dEdλ - dlv);
                    int weight = recursionKernel[lcount][dlcenter];
                    double e = weight * gaussianMag * Math.exp(-lv2 / (2.0 * ls2)) * Math.exp(-dlv2 / (2.0 * dEdLs2));
                    recursionKernelEnergy += e;
                    dgdL += -lv / ls2 * e;
                    dgdEdL += -dlv / dEdLs2 * e;
                }
            }

            /**
             * λ gradient due recursion kernel G(λ, dEdλ).
             */
            dEdλ += dgdL + dgdEdL * d2Edλ2;

            /**
             * Atomic gradient due to recursion kernel G(λ, dEdλ).
             */
            double dUdXdL[] = new double[nAtoms * 3];
            vanderWaals.getdEdLambdaGradient(dUdXdL);
            double grad[] = new double[3];
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                atom.getXYZGradient(grad);
                grad[0] += dgdEdL * dUdXdL[i * 3];
                grad[1] += dgdEdL * dUdXdL[i * 3 + 1];
                grad[2] += dgdEdL * dUdXdL[i * 3 + 2];
                atom.setXYZGradient(grad[0], grad[1], grad[2]);
            }

            // Update free energy F(λ) every ~100 steps.
            if (energyCount % 100 == 0) {
                double freeEnergy = 0.0;
                double binHalf = 0.5 / lambdaBins;
                logger.info("λ              dAdλ           Cumulative");
                for (int i = 0; i < lambdaBins; i++) {
                    double lc = binHalf + i * λBinWidth;

                    int ul = -1;
                    int ll = -1;
                    // Find the smallest dUdL bin.
                    for (int j = 0; j < dEdλBins; j++) {
                        int count = recursionKernel[i][j];
                        if (count != 0 && ll == -1) {
                            ll = j;
                            break;
                        }
                    }
                    // Find the largest dUdL bin.
                    for (int j = dEdλBins - 1; j >= 0; j--) {
                        int count = recursionKernel[i][j];
                        if (count != 0 && ul == -1) {
                            ul = j;
                            break;
                        }
                    }

                    if (ul == -1) {
                        dAdλ[i] = 0.0;
                    } else {
                        double wdUdL = 0.0;
                        double part = 0.0;
                        for (int j = ll; j <= ul; j++) {
                            double dUdLc = mindEdλ + (j + 0.5) * dEdλWidth;
                            double e = Math.exp(gKernel(lc, dUdLc) / (R * 300.0));
                            wdUdL += dUdLc * e;
                            part += e;
                        }
                        dAdλ[i] = wdUdL / part;
                    }
                    freeEnergy += dAdλ[i] * λBinWidth;
                    logger.info(String.format("%15.8f %15.8f %15.8f",
                                              i * λBinWidth + λBinHalfWidth, dAdλ[i], freeEnergy));
                }
            }

            /**
             * Force interpolation due to recursion slave F(λ).
             */
            if (lambda > λBinHalfWidth && lambda < 1.0 - λBinHalfWidth) {
                double binCenter = lambdaBin * λBinWidth + λBinHalfWidth;
                int lb;
                int ub;
                if (lambda > binCenter) {
                    lb = lambdaBin;
                } else {
                    lb = lambdaBin - 1;
                }
                ub = lb + 1;
                double dAdLm = dAdλ[lb];
                double dAdLp = dAdλ[ub];
                double m1c = lb * λBinWidth + λBinHalfWidth;
                double p1c = ub * λBinWidth + λBinHalfWidth;
                dEdλ -= ((lambda - m1c) * dAdLp + (p1c - lambda) * dAdLm) / λBinWidth;
            } else if (lambda <= λBinHalfWidth) {
                double mlc = λBinHalfWidth;
                double plc = λBinWidth + λBinHalfWidth;
                dEdλ -= ((lambda - mlc) * dAdλ[1] + (plc - lambda) * dAdλ[0]) / λBinWidth;
            } else {
                double mlc = 1.0 - 1.5 * λBinWidth;
                double plc = 1.0 - λBinHalfWidth;
                dEdλ -= ((lambda - mlc) * dAdλ[lambdaBins - 1] + (plc - lambda) * dAdλ[lambdaBins - 2]) / λBinWidth;
            }

            /**
             * Meta-dynamic grid counts (every ~10 steps).
             */
            if (energyCount % 10 == 0) {
                recursionKernel[lambdaBin][dUdLBin]++;
            }
        }

        totalTime = System.nanoTime() - totalTime;

        totalBondedEnergy = bondEnergy + angleEnergy + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy;
        totalNonBondedEnergy = vanDerWaalsEnergy + totalElectrostaticEnergy;
        totalEnergy = totalBondedEnergy + totalNonBondedEnergy + solvationEnergy + recursionKernelEnergy;
        
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

    public double gKernel(double lambda, double dUdL) {

        int lambdaBin = (int) Math.ceil(lambda * lambdaBins) - 1;
        int dUdLBin = (int) Math.ceil((dUdL - mindEdλ) / dEdλSpan * dEdλBins) - 1;

        double sum = 0.0;
        double ls2 = 2.0 / lambdaBins * 2.0 / lambdaBins;
        double dUdLs2 = dEdλWidth * 2.0 * dEdλWidth * 2.0;

        for (int l = -5; l <= 5; l++) {
            int lcenter = lambdaBin + 1;
            double lv = lcenter / lambdaBins + 0.5 / lambdaBins;
            double lv2 = (lambda - lv) * (lambda - lv);

            // Mirror condition for Lambda counts.
            int lcount = lcenter;
            if (lcount < 0) {
                lcount = -lcount;
            }
            if (lcount >= lambdaBins) {
                lcount = lambdaBins - (lcount - lambdaBins) - 1;
            }

            for (int dl = -5; dl <= 5; dl++) {
                int dlcenter = dUdLBin + dl;
                double dlv = dlcenter * dEdλWidth + dEdλWidth / 2.0;
                double dlv2 = (dUdL - dlv) * (dUdL - dlv);
                int weight = recursionKernel[lcount][dlcenter];
                double e = weight * gaussianMag * Math.exp(-lv2 / (2.0 * ls2)) * Math.exp(-dlv2 / (2.0 * dUdLs2));
                sum += e;
            }
        }
        return sum;
    }

    public double getTotal() {
        return totalEnergy;
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
                                    "Recursion Kernel  ", solvationEnergy));
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

    public void setSoftCoreLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            vanderWaals.setLambda(lambda);
        } else {
            String message = String.format("Softcore lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    public void setElectrostaticsLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            particleMeshEwald.setLambda(lambda);
        } else {
            String message = String.format("Electrostatics lambda value %8.3f is not in the range [0..1].", lambda);
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
        this.lambdaTerm = lambdaGradient;
        if (vanderWaalsTerm) {
            vanderWaals.lambdaGradient(lambdaGradient);
        }
        energyCount = 0;
        if (recursionKernel == null) {
            double dUdLestimate1 = 0.0;
            double dUdLestimate0 = 0.0;
            maxdEdλ = dUdLestimate1 + 500.0;
            mindEdλ = dUdLestimate1 - 500.0;
            dEdλBins = (int) Math.floor((maxdEdλ - mindEdλ) / dEdλWidth);
            recursionKernel = new int[lambdaBins][dEdλBins];
            dAdλ = new double[lambdaBins];
        }
    }

    @Override
    public void getdEdLambdaGradient(double gradients[]) {
        if (vanderWaalsTerm) {
            vanderWaals.getdEdLambdaGradient(gradients);
        }
    }

    @Override
    public void setLambda(double lambda) {
        this.lambda = lambda;
        if (vanderWaalsTerm) {
            vanderWaals.setLambda(lambda);
        }
        if (multipoleTerm) {
            particleMeshEwald.setLambda(lambda);
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
