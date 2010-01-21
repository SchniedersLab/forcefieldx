/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2009
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

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.numerics.Optimizable;
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
public class PotentialEnergy implements Optimizable {

    private static final Logger logger = Logger.getLogger(PotentialEnergy.class.getName());
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
    protected int nVanDerWaals;
    protected final boolean bondTerm;
    protected final boolean angleTerm;
    protected final boolean stretchBendTerm;
    protected final boolean ureyBradleyTerm;
    protected final boolean outOfPlaneBendTerm;
    protected final boolean torsionTerm;
    protected final boolean piOrbitalTorsionTerm;
    protected final boolean torsionTorsionTerm;
    protected final boolean vanDerWaalsTerm;
    protected final boolean multipoleTerm;
    protected final boolean polarizationTerm;
    protected double bondEnergy;
    protected double angleEnergy;
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
    protected double totalEnergy;
    protected long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    protected long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime;
    protected long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    protected long totalTime;
    protected double[] optimizationScaling = null;
    private static final double toSeconds = 0.000000001;

    public PotentialEnergy(MolecularAssembly molecularAssembly) {
        parallelTeam = new ParallelTeam();
        logger.info(" Number of Threads: " + parallelTeam.getThreadCount());

        ForceField forceField = molecularAssembly.getForceField();

        bondTerm = forceField.getBoolean(ForceFieldBoolean.BONDTERM, true);
        angleTerm = forceField.getBoolean(ForceFieldBoolean.ANGLETERM, true);
        stretchBendTerm = forceField.getBoolean(ForceFieldBoolean.STRBNDTERM, true);
        ureyBradleyTerm = forceField.getBoolean(ForceFieldBoolean.UREYTERM, true);
        outOfPlaneBendTerm = forceField.getBoolean(ForceFieldBoolean.OPBENDTERM, true);
        torsionTerm = forceField.getBoolean(ForceFieldBoolean.TORSIONTERM, true);
        piOrbitalTorsionTerm = forceField.getBoolean(ForceFieldBoolean.PITORSTERM, true);
        torsionTorsionTerm = forceField.getBoolean(ForceFieldBoolean.TORTORTERM, true);
        vanDerWaalsTerm = forceField.getBoolean(ForceFieldBoolean.VDWTERM, true);
        multipoleTerm = forceField.getBoolean(ForceFieldBoolean.MPOLETERM, true);
        polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);

        // Collect, count, pack and sort atoms.
        ArrayList<Atom> atomList = molecularAssembly.getAtomList();
        nAtoms = atomList.size();
        atoms = atomList.toArray(new Atom[nAtoms]);
        Arrays.sort(atoms);

        // Collect, count, pack and sort bonds.
        if (bondTerm) {
            ArrayList<ROLS> bond = molecularAssembly.getBondList();
            nBonds = bond.size();
            bonds = bond.toArray(new Bond[nBonds]);
            Arrays.sort(bonds);
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
        } else {
            nStretchBends = 0;
            stretchBends = null;
        }

        // Collect, count, pack and sort Urey-Brandleys.
        if (ureyBradleyTerm) {
            ArrayList<ROLS> ureyBradley = molecularAssembly.getUreyBradleyList();
            nUreyBradleys = ureyBradley.size();
            ureyBradleys = ureyBradley.toArray(new UreyBradley[nUreyBradleys]);
            Arrays.sort(ureyBradleys);
        } else {
            nUreyBradleys = 0;
            ureyBradleys = null;
        }


        // Collect, count, pack and sort out-of-plane-bends.
        if (outOfPlaneBendTerm) {
            ArrayList<ROLS> outOfPlaneBend = molecularAssembly.getOutOfPlaneBendList();
            nOutOfPlaneBends = outOfPlaneBend.size();
            outOfPlaneBends = outOfPlaneBend.toArray(new OutOfPlaneBend[nOutOfPlaneBends]);
            Arrays.sort(outOfPlaneBends);
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<ROLS> torsion = molecularAssembly.getTorsionList();
            nTorsions = torsion.size();
            torsions = torsion.toArray(new Torsion[nTorsions]);
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<ROLS> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = piOrbitalTorsion.size();
            piOrbitalTorsions = piOrbitalTorsion.toArray(new PiOrbitalTorsion[nPiOrbitalTorsions]);
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<ROLS> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = torsionTorsion.size();
            torsionTorsions = torsionTorsion.toArray(new TorsionTorsion[nTorsionTorsions]);
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }


        // Determine the unit cell dimensions and Spacegroup
        final double a = forceField.getDouble(ForceFieldDouble.A_AXIS, 10.0);
        final double b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
        final double c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
        final double alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
        final double beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
        final double gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
        final String spacegroup = forceField.getString(
                ForceFieldString.SPACEGROUP, "P1");
        crystal = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
        logger.info(crystal.toString());

        if (vanDerWaalsTerm) {
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

        // Zero out the potential energy of each non-bonded term.
        vanDerWaalsEnergy = 0.0;
        permanentMultipoleEnergy = 0.0;
        polarizationEnergy = 0.0;
        totalElectrostaticEnergy = 0.0;
        totalNonBondedEnergy = 0.0;

        // Zero out the total potential energy.
        totalEnergy = 0.0;

        // Zero out the Cartesian coordinate gradient for each atom.
        if (gradient) {
            for (Atom atom : atoms) {
                atom.setXYZGradient(0.0, 0.0, 0.0);
            }
        }

        if (bondTerm) {
            bondTime = System.nanoTime();
            for (int i = 0; i < nBonds; i++) {
                bondEnergy += bonds[i].energy(gradient);
            }
            bondTime = System.nanoTime() - bondTime;
        }

        if (angleTerm) {
            angleTime = System.nanoTime();
            for (int i = 0; i < nAngles; i++) {
                angleEnergy += angles[i].energy(gradient);
            }
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

        if (vanDerWaalsTerm) {
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
            electrostaticTime = System.nanoTime() - electrostaticTime;
        }
        totalTime = System.nanoTime() - totalTime;

        totalBondedEnergy = bondEnergy + angleEnergy + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy + piOrbitalTorsionEnergy + torsionTorsionEnergy;
        totalNonBondedEnergy = vanDerWaalsEnergy + totalElectrostaticEnergy;
        totalEnergy = totalBondedEnergy + totalNonBondedEnergy;
        if (print) {
            StringBuffer sb = new StringBuffer("\n");
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

    public double getTotal() {
        return totalEnergy;
    }

    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer("\n");
        if (bondTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                    "Bond Streching    ", bondEnergy, bonds.length,
                    bondTime * toSeconds));
        }
        if (angleTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                    "Angle Bending     ", angleEnergy, angles.length,
                    angleTime * toSeconds));
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
        if (vanDerWaalsTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                    "Van der Waals     ", vanDerWaalsEnergy,
                    nVanDerWaals, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm) {
            sb.append(String.format(" %s %16.8f %12d\n",
                    "Atomic Multipoles ", permanentMultipoleEnergy,
                    particleMeshEwald.getInteractions()));
        }
        if (polarizationTerm) {
            sb.append(String.format(" %s %16.8f %12d %12.3f\n",
                    "Polarization      ", polarizationEnergy,
                    particleMeshEwald.getInteractions(), electrostaticTime * toSeconds));
        }
        sb.append(String.format("\n %s %16.8f  %s %12.3f (sec)\n",
                "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));
        int nsymm = crystal.spaceGroup.symOps.size();
        if (nsymm > 1) {
            sb.append(String.format(" %s %16.8f\n", "Unit Cell         ",
                    totalEnergy * nsymm));
        }
        return sb.toString();
    }

    public Crystal getCrystal() {
        return crystal;
    }

    @Override
    public void setOptimizationScaling(double scaling[]) {
        if (scaling != null && scaling.length == nAtoms * 3) {
            optimizationScaling = scaling;
        } else {
            optimizationScaling = null;
        }
    }

    @Override
    public double[] getOptimizationScaling() {
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

    private void getGradients(double g[]) {
        assert (g != null && g.length == nAtoms * 3);
        double grad[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZGradient(grad);
            g[index++] = grad[0];
            g[index++] = grad[1];
            g[index++] = grad[2];
        }
    }

    private void setCoordinates(double coords[]) {
        assert (coords != null && coords.length == nAtoms * 3);
        int index = 0;
        for (Atom a : atoms) {
            double x = coords[index++];
            double y = coords[index++];
            double z = coords[index++];
            a.moveTo(x, y, z);
        }
    }

    public void getCoordinates(double x[]) {
        assert (x != null && x.length == nAtoms * 3);
        double xyz[] = new double[3];
        int index = 0;
        for (Atom a : atoms) {
            a.getXYZ(xyz);
            x[index++] = xyz[0];
            x[index++] = xyz[1];
            x[index++] = xyz[2];
        }
    }
}
