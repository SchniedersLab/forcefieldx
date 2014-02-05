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
package ffx.potential;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.potential.nonbonded.NCSRestraint;
import ffx.crystal.ReplicatesCrystal;
import ffx.numerics.Potential;
import ffx.potential.bonded.*;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;

/**
 * Compute the potential energy and derivatives of an AMOEBA system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 *
 */
public class ForceFieldEnergy implements Potential, LambdaInterface {

    private static final Logger logger = Logger.getLogger(ForceFieldEnergy.class.getName());
    private static final double toSeconds = 0.000000001;
    private final Atom[] atoms;
    private Crystal crystal;
    private final ParallelTeam parallelTeam;
    private final Bond bonds[];
    private final Angle angles[];
    private final StretchBend stretchBends[];
    private final UreyBradley ureyBradleys[];
    private final OutOfPlaneBend outOfPlaneBends[];
    private final Torsion torsions[];
    private final PiOrbitalTorsion piOrbitalTorsions[];
    private final TorsionTorsion torsionTorsions[];
    private RestraintBond restraintBonds[];
    private final VanDerWaals vanderWaals;
    private final ParticleMeshEwald particleMeshEwald;
    private final NCSRestraint ncsRestraint;
    private final CoordRestraint coordRestraint;
    protected final int nAtoms;
    protected final int nBonds;
    protected final int nAngles;
    protected final int nStretchBends;
    protected final int nUreyBradleys;
    protected final int nOutOfPlaneBends;
    protected final int nTorsions;
    protected final int nPiOrbitalTorsions;
    protected final int nTorsionTorsions;
    protected int nRestraintBonds;
    protected int nVanDerWaals, nPME, nGK;
    protected boolean bondTerm;
    protected boolean angleTerm;
    protected boolean stretchBendTerm;
    protected boolean ureyBradleyTerm;
    protected boolean outOfPlaneBendTerm;
    protected boolean torsionTerm;
    protected boolean piOrbitalTorsionTerm;
    protected boolean torsionTorsionTerm;
    protected boolean restraintBondTerm;
    protected boolean vanderWaalsTerm;
    protected boolean multipoleTerm;
    protected boolean polarizationTerm;
    protected boolean generalizedKirkwoodTerm;
    protected boolean ncsTerm;
    protected boolean restrainTerm;
    protected boolean lambdaBondedTerms = false;
    protected boolean bondTermOrig;
    protected boolean angleTermOrig;
    protected boolean stretchBendTermOrig;
    protected boolean ureyBradleyTermOrig;
    protected boolean outOfPlaneBendTermOrig;
    protected boolean torsionTermOrig;
    protected boolean piOrbitalTorsionTermOrig;
    protected boolean torsionTorsionTermOrig;
    protected boolean restraintBondTermOrig;
    protected boolean vanderWaalsTermOrig;
    protected boolean multipoleTermOrig;
    protected boolean polarizationTermOrig;
    protected boolean generalizedKirkwoodTermOrig;
    protected boolean ncsTermOrig;
    protected boolean restrainTermOrig;
    protected double bondEnergy, bondRMSD;
    protected double angleEnergy, angleRMSD;
    protected double stretchBendEnergy;
    protected double ureyBradleyEnergy;
    protected double outOfPlaneBendEnergy;
    protected double torsionEnergy;
    protected double piOrbitalTorsionEnergy;
    protected double torsionTorsionEnergy;
    protected double restraintBondEnergy;
    protected double totalBondedEnergy;
    protected double vanDerWaalsEnergy;
    protected double permanentMultipoleEnergy;
    protected double polarizationEnergy;
    protected double totalElectrostaticEnergy;
    protected double totalNonBondedEnergy;
    protected double solvationEnergy;
    protected double ncsEnergy;
    protected double restrainEnergy;
    protected double totalEnergy;
    protected long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    protected long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime;
    protected long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    protected long restraintBondTime, ncsTime, coordRestraintTime;
    protected long totalTime;
    protected double lambda = 1.0;
    protected double[] optimizationScaling = null;
    protected VARIABLE_TYPE[] variableTypes = null;
    private double xyz[] = null;

    /**
     * <p>
     * Constructor for ForceFieldEnergy.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.bonded.MolecularAssembly}
     * object.
     */
    public ForceFieldEnergy(MolecularAssembly molecularAssembly) {
        parallelTeam = new ParallelTeam();
        logger.info(format(" Constructing Force Field"));
        logger.info(format("\n SMP threads:                        %10d", parallelTeam.getThreadCount()));

        // Get a reference to the sorted atom array.
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        xyz = new double[nAtoms * 3];

        // Check that atom ordering is correct.
        for (int i = 0; i < nAtoms; i++) {
            int index = atoms[i].xyzIndex - 1;
            assert (i == index);
        }

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
        if (vanderWaalsTerm) {
            multipoleTerm = forceField.getBoolean(ForceFieldBoolean.MPOLETERM, true);
            if (multipoleTerm) {
                polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);
                generalizedKirkwoodTerm = forceField.getBoolean(ForceFieldBoolean.GKTERM, false);
            } else {
                /**
                 * If multipole electrostatics is turned off, turn off all
                 * electrostatics.
                 */
                polarizationTerm = false;
                generalizedKirkwoodTerm = false;
            }
        } else {
            /**
             * If van der Waals is turned off, turn off all non-bonded terms.
             */
            multipoleTerm = false;
            polarizationTerm = false;
            generalizedKirkwoodTerm = false;
        }
        restraintBondTerm = false;
        restrainTerm = forceField.getBoolean(ForceFieldBoolean.RESTRAINTERM, false);

        //For respa
        bondTermOrig = bondTerm;
        angleTermOrig = angleTerm;
        stretchBendTermOrig = stretchBendTerm;
        ureyBradleyTermOrig = ureyBradleyTerm;
        outOfPlaneBendTermOrig = outOfPlaneBendTerm;
        torsionTermOrig = torsionTerm;
        piOrbitalTorsionTermOrig = piOrbitalTorsionTerm;
        torsionTorsionTermOrig = torsionTorsionTerm;
        restraintBondTermOrig = restraintBondTerm;
        vanderWaalsTermOrig = vanderWaalsTerm;
        multipoleTermOrig = multipoleTerm;
        polarizationTermOrig = polarizationTerm;
        generalizedKirkwoodTermOrig = generalizedKirkwoodTerm;
        ncsTermOrig = ncsTerm;
        restrainTermOrig = restrainTerm;

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
            spacegroup = forceField.getString(ForceFieldString.SPACEGROUP, "P 1");
        } catch (Exception e) {
            logger.info(" The system will be treated as aperiodic.");
            aperiodic = true;

            double maxr = 10.0;
            for (int i = 0; i < nAtoms - 1; i++) {
                Atom ai = atoms[i];
                for (int j = 1; j < nAtoms; j++) {
                    Atom aj = atoms[j];
                    double dx = ai.getX() - aj.getX();
                    double dy = ai.getY() - aj.getY();
                    double dz = ai.getZ() - aj.getZ();
                    double r = sqrt(dx * dx + dy * dy + dz * dz);
                    maxr = max(r, maxr);
                }
            }

            /**
             * Turn off reciprocal space calculations.
             */
            forceField.addForceFieldDouble(ForceFieldDouble.EWALD_ALPHA, 0.0);
            // Specify some dummy values for the crystal.
            spacegroup = "P1";
            a = 4.0 * maxr;
            b = 4.0 * maxr;
            c = 4.0 * maxr;
            alpha = 90.0;
            beta = 90.0;
            gamma = 90.0;
        }
        Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
        unitCell.setAperiodic(aperiodic);

        /**
         * If necessary, create a ReplicatesCrystal.
         */
        if (!aperiodic) {
            this.crystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, cutOff2);
            logger.info(crystal.toString());
        } else {
            this.crystal = unitCell;
        }

        boolean rigidHydrogens = forceField.getBoolean(ForceFieldBoolean.RIGID_HYDROGENS, false);
        double rigidScale = forceField.getDouble(ForceFieldDouble.RIGID_SCALE, 10.0);

        if (rigidScale <= 1.0) {
            rigidScale = 1.0;
        }

        logger.info("\n Bonded Terms");
        if (rigidHydrogens && rigidScale > 1.0) {
            logger.info(format("  Rigid hydrogens:                   %10.2f", rigidScale));
        }

        // Collect, count, pack and sort bonds.
        if (bondTerm) {
            ArrayList<ROLS> bond = molecularAssembly.getBondList();
            nBonds = bond.size();
            bonds = bond.toArray(new Bond[nBonds]);
            Arrays.sort(bonds);
            logger.info(format("  Bonds:                             %10d", nBonds));
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
            logger.info(format("  Angles:                            %10d", nAngles));
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
            logger.info(format("  Stretch-Bends:                     %10d", nStretchBends));
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
            logger.info(format("  Urey-Bradleys:                     %10d", nUreyBradleys));
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
            logger.info(format("  Out-of-Plane Bends:                %10d", nOutOfPlaneBends));
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<ROLS> torsion = molecularAssembly.getTorsionList();
            nTorsions = torsion.size();
            torsions = torsion.toArray(new Torsion[nTorsions]);
            logger.info(format("  Torsions:                          %10d", nTorsions));
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<ROLS> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = piOrbitalTorsion.size();
            piOrbitalTorsions = piOrbitalTorsion.toArray(new PiOrbitalTorsion[nPiOrbitalTorsions]);
            logger.info(format("  Pi-Orbital Torsions:               %10d", nPiOrbitalTorsions));
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<ROLS> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = torsionTorsion.size();
            torsionTorsions = torsionTorsion.toArray(new TorsionTorsion[nTorsionTorsions]);
            logger.info(format("  Torsion-Torsions:                  %10d", nTorsionTorsions));
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }

        logger.info("\n Non-Bonded Terms");

        if (vanderWaalsTerm) {
            String name = forceField.toString().toUpperCase();
            if (name.contains("OPLS")) {
                vanderWaals = new VanDerWaals(molecularAssembly, crystal, parallelTeam, 12.0, 6.0, 0.0, 0.0);
            } else {
                vanderWaals = new VanDerWaals(molecularAssembly, crystal, parallelTeam, 14.0, 7.0, 0.07, 0.12);
            }
        } else {
            vanderWaals = null;
        }

        if (multipoleTerm) {
            particleMeshEwald = new ParticleMeshEwald(molecularAssembly, crystal,
                    vanderWaals.getNeighborList(), parallelTeam);
        } else {
            particleMeshEwald = null;
        }

        if (unitCell.spaceGroup.number == 1) {
            ncsTerm = forceField.getBoolean(ForceFieldBoolean.NCSTERM, false);
            ncsTermOrig = ncsTerm;
            if (ncsTerm) {
                String sg = forceField.getString(ForceFieldString.NCSGROUP, "P 1");
                Crystal ncsCrystal = new Crystal(a, b, c, alpha, beta, gamma, sg);
                ncsRestraint = new NCSRestraint(molecularAssembly, ncsCrystal);
            } else {
                ncsRestraint = null;
            }
        } else {
            ncsTerm = false;
            ncsTermOrig = false;
            ncsRestraint = null;
        }

        if (restrainTerm) {
            this.coordRestraint = new CoordRestraint(molecularAssembly, crystal);
        } else {
            coordRestraint = null;
        }

        molecularAssembly.setPotential(this);
    }

    public void setLambdaBondedTerms(boolean lambdaBondedTerms) {
        this.lambdaBondedTerms = lambdaBondedTerms;
    }

    public void setNoSoftCoreElectrostatics(boolean noElec) {
        if (particleMeshEwald != null) {
            particleMeshEwald.setNoSoftCoreElectrostatics(noElec);
        }
    }

    /**
     * <p>
     * energy</p>
     *
     * @param gradient a boolean.
     * @param print a boolean.
     * @return a double.
     */
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
        restraintBondTime = 0;
        ncsTime = 0;
        coordRestraintTime = 0;
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

        // Zero out potential energy of restraint terms
        restraintBondEnergy = 0.0;
        ncsEnergy = 0.0;
        restrainEnergy = 0.0;

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
                if (lambdaBondedTerms && !b.applyLambda()) {
                    continue;
                }
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
                if (lambdaBondedTerms && !a.applyLambda()) {
                    continue;
                }
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
                StretchBend sb = stretchBends[i];
                if (lambdaBondedTerms && !sb.applyLambda()) {
                    continue;
                }
                stretchBendEnergy += sb.energy(gradient);
            }
            stretchBendTime = System.nanoTime() - stretchBendTime;
        }

        if (ureyBradleyTerm) {
            ureyBradleyTime = System.nanoTime();
            for (int i = 0; i < nUreyBradleys; i++) {
                UreyBradley ub = ureyBradleys[i];
                if (lambdaBondedTerms && !ub.applyLambda()) {
                    continue;
                }
                ureyBradleyEnergy += ub.energy(gradient);
            }
            ureyBradleyTime = System.nanoTime() - ureyBradleyTime;
        }

        if (outOfPlaneBendTerm) {
            outOfPlaneBendTime = System.nanoTime();
            for (int i = 0; i < nOutOfPlaneBends; i++) {
                OutOfPlaneBend ob = outOfPlaneBends[i];
                if (lambdaBondedTerms && !ob.applyLambda()) {
                    continue;
                }
                outOfPlaneBendEnergy += ob.energy(gradient);
            }
            outOfPlaneBendTime = System.nanoTime() - outOfPlaneBendTime;
        }

        if (torsionTerm) {
            torsionTime = System.nanoTime();
            for (int i = 0; i < nTorsions; i++) {
                Torsion t = torsions[i];
                if (lambdaBondedTerms && !t.applyLambda()) {
                    continue;
                }
                torsionEnergy += t.energy(gradient);
            }
            torsionTime = System.nanoTime() - torsionTime;
        }

        if (piOrbitalTorsionTerm) {
            piOrbitalTorsionTime = System.nanoTime();
            for (int i = 0; i < nPiOrbitalTorsions; i++) {
                PiOrbitalTorsion po = piOrbitalTorsions[i];
                if (lambdaBondedTerms && !po.applyLambda()) {
                    continue;
                }
                piOrbitalTorsionEnergy += po.energy(gradient);
            }
            piOrbitalTorsionTime = System.nanoTime() - piOrbitalTorsionTime;
            torsionTorsionTime = System.nanoTime();
        }

        if (torsionTorsionTerm) {
            for (int i = 0; i < nTorsionTorsions; i++) {
                TorsionTorsion tt = torsionTorsions[i];
                if (lambdaBondedTerms && !tt.applyLambda()) {
                    continue;
                }
                torsionTorsionEnergy += torsionTorsions[i].energy(gradient);
            }
            torsionTorsionTime = System.nanoTime() - torsionTorsionTime;
        }

        if (restraintBondTerm) {
            restraintBondTime = System.nanoTime();
            for (int i = 0; i < nRestraintBonds; i++) {
                RestraintBond rb = restraintBonds[i];
                if (lambdaBondedTerms && !rb.applyLambda()) {
                    continue;
                }
                restraintBondEnergy += rb.energy(gradient);
            }
            restraintBondTime = System.nanoTime() - restraintBondTime;
        }

        if (restrainTerm && lambdaBondedTerms) {
            coordRestraintTime = -System.nanoTime();
            restrainEnergy = coordRestraint.residual(gradient, print);
            coordRestraintTime += System.nanoTime();
        }

        if (!lambdaBondedTerms) {

            if (ncsTerm) {
                ncsTime = -System.nanoTime();
                ncsEnergy = ncsRestraint.residual(gradient, print);
                ncsTime += System.nanoTime();
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
        }

        totalTime = System.nanoTime() - totalTime;

        totalBondedEnergy = bondEnergy + restraintBondEnergy + angleEnergy
                + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy
                + torsionEnergy + piOrbitalTorsionEnergy
                + torsionTorsionEnergy + ncsEnergy + restrainEnergy;
        totalNonBondedEnergy = vanDerWaalsEnergy + totalElectrostaticEnergy;
        totalEnergy = totalBondedEnergy + totalNonBondedEnergy + solvationEnergy;

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
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return totalEnergy;
    }

    /**
     * <p>
     * getPDBHeaderString</p>
     *
     * @return a {@link java.lang.String} object.
     */
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
        if (restraintBondTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "RESTRAINT BOND STRETCHING            : ", restraintBondEnergy, restraintBonds.length));
        }

        if (ncsTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "NCS RESRAINT               : ", ncsEnergy, nAtoms));
        }

        if (restrainTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "COORDINATE RESRAINT        : ", restrainEnergy, nAtoms));
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

    /**
     * {@inheritDoc}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder("\n");
        if (bondTerm && nBonds > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f (%8.5f)\n",
                    "Bond Streching    ", bondEnergy, nBonds,
                    bondTime * toSeconds, bondRMSD));
        }
        if (angleTerm && nAngles > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f (%8.5f)\n",
                    "Angle Bending     ", angleEnergy, nAngles,
                    angleTime * toSeconds, angleRMSD));
        }
        if (stretchBendTerm && nStretchBends > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Stretch-Bend      ", stretchBendEnergy,
                    nStretchBends, stretchBendTime * toSeconds));
        }
        if (ureyBradleyTerm && nUreyBradleys > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Urey-Bradley      ", ureyBradleyEnergy,
                    nUreyBradleys, ureyBradleyTime * toSeconds));
        }
        if (outOfPlaneBendTerm && nOutOfPlaneBends > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Out-of-Plane Bend ", outOfPlaneBendEnergy,
                    nOutOfPlaneBends, outOfPlaneBendTime * toSeconds));
        }
        if (torsionTerm && nTorsions > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Torsional Angle   ", torsionEnergy, nTorsions,
                    torsionTime * toSeconds));
        }
        if (piOrbitalTorsionTerm && nPiOrbitalTorsions > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Pi-Orbital Torsion", piOrbitalTorsionEnergy,
                    nPiOrbitalTorsions, piOrbitalTorsionTime * toSeconds));
        }
        if (torsionTorsionTerm && nTorsionTorsions > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Torsion-Torsion   ", torsionTorsionEnergy,
                    nTorsionTorsions, torsionTorsionTime * toSeconds));
        }
        if (restraintBondTerm && nRestraintBonds > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Bond Restraint    ", restraintBondEnergy, nRestraintBonds,
                    restraintBondTime * toSeconds));
        }
        if (ncsTerm) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "NCS Restraint     ", ncsEnergy, nAtoms,
                    ncsTime * toSeconds));
        }
        if (restrainTerm) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Coord. Restraint  ", restrainEnergy, nAtoms,
                    coordRestraintTime * toSeconds));
        }
        if (vanderWaalsTerm && nVanDerWaals > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Van der Waals     ", vanDerWaalsEnergy,
                    nVanDerWaals, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm && nPME > 0) {
            if (polarizationTerm) {
                sb.append(String.format("  %s %16.8f %12d\n",
                        "Atomic Multipoles ", permanentMultipoleEnergy, nPME));
            } else {
                sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                        "Atomic Multipoles ", permanentMultipoleEnergy, nPME, electrostaticTime * toSeconds));
            }
        }
        if (polarizationTerm && nPME > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Polarization      ", polarizationEnergy,
                    nPME, electrostaticTime * toSeconds));
        }
        if (generalizedKirkwoodTerm && nGK > 0) {
            sb.append(String.format("  %s %16.8f %12d\n",
                    "Solvation         ", solvationEnergy, nGK));
        }

        sb.append(String.format("  %s %16.8f  %s %12.3f (sec)\n",
                "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));

        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(String.format("  %s %16.8f\n", "Unit Cell         ",
                    totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(String.format("  %s %16.8f\n", "Replicates Cell   ",
                    totalEnergy * nsymm));
        }

        return sb.toString();
    }

    public ParallelTeam getParallelTeam() {
        return parallelTeam;
    }

    /**
     * <p>
     * Getter for the field <code>crystal</code>.</p>
     *
     * @return a {@link ffx.crystal.Crystal} object.
     */
    public Crystal getCrystal() {
        return crystal;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (lambda <= 1.0 && lambda >= 0.0) {
            this.lambda = lambda;
            if (vanderWaalsTerm) {
                vanderWaals.setLambda(lambda);
            }
            if (multipoleTerm) {
                particleMeshEwald.setLambda(lambda);
            }
            if (restraintBondTerm && restraintBonds != null) {
                for (int i = 0; i < restraintBonds.length; i++) {
                    restraintBonds[i].setLambda(lambda);
                }
            }
            if (ncsTerm && ncsRestraint != null) {
                ncsRestraint.setLambda(lambda);
            }
            if (restrainTerm && coordRestraint != null) {
                coordRestraint.setLambda(lambda);
            }
        } else {
            String message = String.format("Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setScaling(double scaling[]) {
        optimizationScaling = scaling;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getScaling() {
        return optimizationScaling;
    }

    /**
     * Return a reference to each variables type.
     *
     * @return the type of each variable.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        int n = getNumberOfVariables();
        VARIABLE_TYPE type[] = new VARIABLE_TYPE[n];
        int i = 0;
        for (Atom a : atoms) {
            type[i++] = VARIABLE_TYPE.X;
            type[i++] = VARIABLE_TYPE.Y;
            type[i++] = VARIABLE_TYPE.Z;
        }
        return type;
    }

    @Override
    public double energy(double[] x) {
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
        double e = energy(false, false);
        /**
         * Rescale the coordinates.
         */
        if (optimizationScaling != null) {
            int len = x.length;
            for (int i = 0; i < len; i++) {
                x[i] *= optimizationScaling[i];
            }
        }
        return e;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double x[], double g[]) {
        /**
         * Un-scale the coordinates.
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

    /**
     * <p>
     * getGradients</p>
     *
     * @param g an array of double.
     */
    public void getGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int index = 0;
        int len = atoms.length;
        for (int i = 0; i < len; i++) {
            Atom a = atoms[i];
            a.getXYZGradient(grad);
            double gx = grad[0];
            double gy = grad[1];
            double gz = grad[2];
            if (Double.isNaN(gx) || Double.isInfinite(gx)
                    || Double.isNaN(gy) || Double.isInfinite(gy)
                    || Double.isNaN(gz) || Double.isInfinite(gz)) {
                String message = format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                        a.toString(), gx, gy, gz);
                logger.severe(message);
            }
            g[index++] = gx;
            g[index++] = gy;
            g[index++] = gz;
        }
    }

    private void setCoordinates(double coords[]) {
        assert (coords != null);
        int index = 0;
        int len = atoms.length;
        for (int i = 0; i < len; i++) {
            Atom a = atoms[i];
            double x = coords[index++];
            double y = coords[index++];
            double z = coords[index++];
            a.moveTo(x, y, z);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double x[]) {
        int n = getNumberOfVariables();
        if (x == null || x.length < n) {
            x = new double[n];
        }
        double xi[] = new double[3];
        int index = 0;
        int len = atoms.length;
        for (int i = 0; i < len; i++) {
            Atom a = atoms[i];
            a.getXYZ(xi);
            x[index++] = xi[0];
            x[index++] = xi[1];
            x[index++] = xi[2];
        }
        return x;
    }

    /**
     * {@inheritDoc}
     */
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

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        return nAtoms * 3;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        double dEdLambda = 0.0;
        if (!lambdaBondedTerms) {
            if (vanderWaalsTerm) {
                dEdLambda = vanderWaals.getdEdL();
            }
            if (multipoleTerm) {
                dEdLambda += particleMeshEwald.getdEdL();
            }
            if (restraintBondTerm) {
                for (int i = 0; i < nRestraintBonds; i++) {
                    dEdLambda += restraintBonds[i].getdEdL();
                }
            }
            if (ncsTerm && ncsRestraint != null) {
                dEdLambda += ncsRestraint.getdEdL();
            }
        } else {
            if (restrainTerm && coordRestraint != null) {
                dEdLambda += coordRestraint.getdEdL();
            }
        }
        return dEdLambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double gradients[]) {
        if (!lambdaBondedTerms) {
            if (vanderWaalsTerm) {
                vanderWaals.getdEdXdL(gradients);
            }
            if (multipoleTerm) {
                particleMeshEwald.getdEdXdL(gradients);
            }
            if (restraintBondTerm) {
                for (int i = 0; i < nRestraintBonds; i++) {
                    restraintBonds[i].getdEdXdL(gradients);
                }
            }
            if (ncsTerm && ncsRestraint != null) {
                ncsRestraint.getdEdXdL(gradients);
            }
        } else {
            if (restrainTerm && coordRestraint != null) {
                coordRestraint.getdEdXdL(gradients);
            }
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
    public double getd2EdL2() {
        double d2EdLambda2 = 0.0;
        if (!lambdaBondedTerms) {
            if (vanderWaalsTerm) {
                d2EdLambda2 = vanderWaals.getd2EdL2();
            }
            if (multipoleTerm) {
                d2EdLambda2 += particleMeshEwald.getd2EdL2();
            }
            if (restraintBondTerm) {
                for (int i = 0; i < nRestraintBonds; i++) {
                    d2EdLambda2 += restraintBonds[i].getd2EdL2();
                }
            }
            if (ncsTerm && ncsRestraint != null) {
                d2EdLambda2 += ncsRestraint.getd2EdL2();
            }
        } else {
            if (restrainTerm && coordRestraint != null) {
                d2EdLambda2 += coordRestraint.getd2EdL2();
            }
        }
        return d2EdLambda2;
    }

    /**
     * <p>
     * setRestraintBond</p>
     *
     * @param a1 a {@link ffx.potential.bonded.Atom} object.
     * @param a2 a {@link ffx.potential.bonded.Atom} object.
     * @param distance a double.
     * @param forceConstant the force constant in kcal/mole
     */
    public void setRestraintBond(Atom a1, Atom a2, double distance, double forceConstant) {
        restraintBondTerm = true;
        RestraintBond rb = new RestraintBond(a1, a2, crystal);
        int classes[] = {a1.getAtomType().atomClass, a2.getAtomType().atomClass};
        rb.setBondType((new BondType(classes, forceConstant, distance, BondType.BondFunction.HARMONIC)));
        nRestraintBonds = 1;
        restraintBonds = new RestraintBond[nRestraintBonds];
        restraintBonds[0] = rb;
        rb.energy(false);
        rb.log();
    }

    /**
     * This method is for Respa integrator only.
     *
     * @param state The STATE is FAST, SLOW or BOTH.
     */
    @Override
    public void setEnergyTermState(STATE state) {
        switch (state) {
            case FAST:
                bondTerm = bondTermOrig;
                angleTerm = angleTermOrig;
                stretchBendTerm = stretchBendTermOrig;
                ureyBradleyTerm = ureyBradleyTermOrig;
                outOfPlaneBendTerm = outOfPlaneBendTermOrig;
                torsionTerm = torsionTermOrig;
                piOrbitalTorsionTerm = piOrbitalTorsionTermOrig;
                torsionTorsionTerm = torsionTorsionTermOrig;
                restraintBondTerm = restraintBondTermOrig;
                restrainTerm = restrainTermOrig;
                vanderWaalsTerm = false;
                multipoleTerm = false;
                polarizationTerm = false;
                generalizedKirkwoodTerm = false;
                break;

            case SLOW:
                vanderWaalsTerm = vanderWaalsTermOrig;
                multipoleTerm = multipoleTermOrig;
                polarizationTerm = polarizationTermOrig;
                generalizedKirkwoodTerm = generalizedKirkwoodTermOrig;
                bondTerm = false;
                angleTerm = false;
                stretchBendTerm = false;
                ureyBradleyTerm = false;
                outOfPlaneBendTerm = false;
                torsionTerm = false;
                piOrbitalTorsionTerm = false;
                torsionTorsionTerm = false;
                restraintBondTerm = false;
                restrainTerm = false;
                break;

            default:
                bondTerm = bondTermOrig;
                angleTerm = angleTermOrig;
                stretchBendTerm = stretchBendTermOrig;
                ureyBradleyTerm = ureyBradleyTermOrig;
                outOfPlaneBendTerm = outOfPlaneBendTermOrig;
                torsionTerm = torsionTermOrig;
                piOrbitalTorsionTerm = piOrbitalTorsionTermOrig;
                torsionTorsionTerm = torsionTorsionTermOrig;
                restraintBondTerm = restraintBondTermOrig;
                restrainTermOrig = restrainTerm;
                vanderWaalsTerm = vanderWaalsTermOrig;
                multipoleTerm = multipoleTermOrig;
                polarizationTerm = polarizationTermOrig;
                generalizedKirkwoodTerm = generalizedKirkwoodTermOrig;
        }
    }

    /**
     * Set the boundary conditions for this calculation.
     *
     * @param crystal
     */
    public void setCrystal(Crystal crystal) {
        this.crystal = crystal;
        /**
         * Update VanDerWaals first, in case the NeighborList needs to be
         * re-allocated to include a larger number of replicated cells.
         */
        if (vanderWaalsTerm == true) {
            vanderWaals.setCrystal(crystal);
        }
        if (multipoleTerm == true) {
            particleMeshEwald.setCrystal(crystal);
        }
        /**
         * TODO: update GeneralizedKirkwood to include support for symmetry
         * operators and periodic boundary conditions.
         */
    }
}
