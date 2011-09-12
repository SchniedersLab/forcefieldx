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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Logger;
import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.numerics.Potential;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
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
    private RestraintBond restraintBonds[];
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
    protected int nRestraintBonds;
    protected int nVanDerWaals, nPME, nGK;
    protected final boolean bondTerm;
    protected final boolean angleTerm;
    protected final boolean stretchBendTerm;
    protected final boolean ureyBradleyTerm;
    protected final boolean outOfPlaneBendTerm;
    protected final boolean torsionTerm;
    protected final boolean piOrbitalTorsionTerm;
    protected final boolean torsionTorsionTerm;
    protected boolean restraintBondTerm;
    protected final boolean vanderWaalsTerm;
    protected final boolean multipoleTerm;
    protected final boolean polarizationTerm;
    protected final boolean generalizedKirkwoodTerm;
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
    protected double totalEnergy;
    protected long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    protected long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime;
    protected long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    protected long restraintBondTime;
    protected long totalTime;
    protected double lambda = 1.0;
    protected double[] optimizationScaling = null;
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
        restraintBondTerm = false;


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
         * Do we need a ReplicatesCrystal?
         */
        int l = 1;
        int m = 1;
        int n = 1;
        if (!aperiodic) {
            while (unitCell.a * l < cutOff2) {
                l++;
            }
            while (unitCell.b * m < cutOff2) {
                m++;
            }
            while (unitCell.c * n < cutOff2) {
                n++;
            }
        }

        if (l * m * n > 1) {
            this.crystal = new ReplicatesCrystal(unitCell, l, m, n);
        } else {
            this.crystal = unitCell;
        }

        if (!aperiodic) {
            logger.info(crystal.toString());
        }

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
                                                      vanderWaals.getNeighborLists(), vanderWaals.getPairwiseSchedule());
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

        if (restraintBondTerm) {
            restraintBondTime = System.nanoTime();
            for (int i = 0; i < restraintBonds.length; i++) {
                RestraintBond rb = restraintBonds[i];
                restraintBondEnergy += rb.energy(gradient);
            }
            restraintBondTime = System.nanoTime() - restraintBondTime;
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

        totalTime = System.nanoTime() - totalTime;

        totalBondedEnergy = bondEnergy + restraintBondEnergy + angleEnergy + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy;
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

    @Override
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
        if (restraintBondTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                                    "RESTRAINT BOND STRETCHING            : ", restraintBondEnergy, restraintBonds.length));
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
                                    "Bond Streching    ", bondEnergy, nBonds,
                                    bondTime * toSeconds, bondRMSD));
        }
        if (angleTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f (%8.5f)\n",
                                    "Angle Bending     ", angleEnergy, nAngles,
                                    angleTime * toSeconds, angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Stretch-Bend      ", stretchBendEnergy,
                                    nStretchBends, stretchBendTime * toSeconds));
        }
        if (ureyBradleyTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Urey-Bradley      ", ureyBradleyEnergy,
                                    nUreyBradleys, ureyBradleyTime * toSeconds));
        }
        if (outOfPlaneBendTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Out-of-Plane Bend ", outOfPlaneBendEnergy,
                                    nOutOfPlaneBends, outOfPlaneBendTime * toSeconds));
        }
        if (torsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsional Angle   ", torsionEnergy, nTorsions,
                                    torsionTime * toSeconds));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Pi-Orbital Torsion", piOrbitalTorsionEnergy,
                                    nPiOrbitalTorsions, piOrbitalTorsionTime * toSeconds));
        }
        if (torsionTorsionTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Torsion-Torsion   ", torsionTorsionEnergy,
                                    nTorsionTorsions, torsionTorsionTime * toSeconds));
        }
        if (restraintBondTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Bond Restraint    ", restraintBondEnergy, nRestraintBonds,
                                    restraintBondTime));
        }
        if (vanderWaalsTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                    "Van der Waals     ", vanDerWaalsEnergy,
                                    nVanDerWaals, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm) {
            if (polarizationTerm) {
                sb.append(String.format(" %s %16.8f %12d\n",
                                        "Atomic Multipoles ", permanentMultipoleEnergy, nPME));
            } else {
                sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                                        "Atomic Multipoles ", permanentMultipoleEnergy, nPME, electrostaticTime * toSeconds));
            }
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
                particleMeshEwald.setLambda(lambda);
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
            if (Double.isNaN(gx) || Double.isInfinite(gx)
                || Double.isNaN(gy) || Double.isInfinite(gy)
                || Double.isNaN(gz) || Double.isInfinite(gz)) {
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
    public double getdEdL() {
        double dEdLambda = 0.0;
        if (vanderWaalsTerm) {
            dEdLambda = vanderWaals.getdEdL();
        }
        if (multipoleTerm) {
            dEdLambda += particleMeshEwald.getdEdL();
        }
        return dEdLambda;
    }

    @Override
    public void getdEdXdL(double gradients[]) {
        if (vanderWaalsTerm) {
            vanderWaals.getdEdXdL(gradients);
        }
        if (multipoleTerm) {
            particleMeshEwald.getdEdXdL(gradients);
        }
    }

    @Override
    public double getLambda() {
        return lambda;
    }

    @Override
    public double getd2EdL2() {
        double d2EdLambda2 = 0.0;
        if (vanderWaalsTerm) {
            d2EdLambda2 = vanderWaals.getd2EdL2();
        }
        if (multipoleTerm) {
            d2EdLambda2 += particleMeshEwald.getd2EdL2();
        }
        return d2EdLambda2;
    }

    public void setRestraintBond(Atom a1, Atom a2, double distance) {
        double forceConstant = 1.0;
        restraintBondTerm = true;
        RestraintBond rb = new RestraintBond(a1, a2);
        int classes[] = {a1.getAtomType().atomClass, a2.getAtomType().atomClass};
        rb.setBondType((new BondType(classes, forceConstant, distance)));
        RestraintBond restraint[] = new RestraintBond[nRestraintBonds + 1];
        System.arraycopy(restraintBonds, 0, restraint, 0, nRestraintBonds++);
    }
}
