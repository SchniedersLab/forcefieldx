/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.algorithms;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.random;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.crystal.SpaceGroup;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;

import static ffx.crystal.SpaceGroup.CrystalSystem.CUBIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.HEXAGONAL;
import static ffx.crystal.SpaceGroup.CrystalSystem.MONOCLINIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.ORTHORHOMBIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.TETRAGONAL;
import static ffx.crystal.SpaceGroup.CrystalSystem.TRICLINIC;
import static ffx.crystal.SpaceGroup.CrystalSystem.TRIGONAL;

/**
 * The Barostat class maintains constant pressure using random trial moves in
 * lattice parameters, which are consistent with the space group.
 *
 * @see D. Frenkel and B. Smit, "Understanding Molecular Simulation, 2nd
 * Edition", Academic Press, San Diego, CA, 2002; Section 5.4
 *
 * @author Michael J. Schnieders
 */
public class Barostat implements Potential, LambdaInterface {

    private static final Logger logger = Logger.getLogger(Barostat.class.getName());
    /**
     * Ideal gas constant in kcal/mol/K
     */
    private static final double kB = 1.98720415e-3;
    /**
     * Conversion from kcal/mol/Ang^3 to Atm.
     */
    private static final double PRESCON = 6.85684112e4;
    /**
     * Avogadro's number.
     */
    private static final double AVOGADRO = 6.02214129e23;
    /**
     * Sampling temperature (K).
     */
    private final double temperature = 298.15;
    /**
     * Sampling pressure (atm)
     */
    private double pressure = 1.0;
    /**
     * Ideal gas constant * temperature (kcal/mol).
     */
    private final double kT = temperature * kB;

    /**
     * Flag to turn the Barostat on or off. If false, MC moves will not be
     * tried.
     */
    private boolean active = true;

    /**
     * Default edge length move (A).
     */
    private double maxSideMove = 0.5;
    /**
     * Default angular move (degrees).
     */
    private double maxAngleMove = 1.0;
    /**
     * A carbon atom cannot fit into a unit cell without
     * interfacial radii greater than ~1.2 Angstroms.
     */
    private double minInterfacialRadius = 1.2;
    /**
     * MolecularAssembly being simulated.
     */
    private final MolecularAssembly molecularAssembly;
    /**
     * Atoms in the system.
     */
    private final Atom atoms[];
    /**
     * Number of atoms.
     */
    private final int nAtoms;
    /**
     * Mass of the system.
     */
    private double mass;
    /**
     * ForceFieldEnergy that describes the system.
     */
    private final ForceFieldEnergy potential;
    /**
     * Boundary conditions and symmetry operators (may be a ReplicatedCrystal).
     */
    private Crystal crystal;
    /**
     * The unit cell.
     */
    private Crystal unitCell;
    /**
     * The number of space group symmetry operators.
     */
    private final int nSymm;
    /**
     * Number of independent molecules in the simulation cell.
     */
    private final SpaceGroup spaceGroup;
    /**
     * Number of molecules in the systems.
     */
    private final int nMolecules;
    /**
     * Center of mass of each molecule in fractional coordinates.
     */
    private final double fractionalCOM[][];
    /**
     * Number of energy evaluations between application of MC moves.
     */
    private int meanBarostatInterval = 1;
    /**
     * A counter for the number of barostat calls.
     */
    private int barostatCount = 0;
    /**
     * Number of Monte Carlo moves attempted.
     */
    private int sideMovesAttempted = 0;
    /**
     * Number of Monte Carlo moves accepted.
     */
    private int sideMovesAccepted = 0;
    /**
     * Number of Monte Carlo moves attempted.
     */
    private int angleMovesAttempted = 0;
    /**
     * Number of Monte Carlo moves accepted.
     */
    private int angleMovesAccepted = 0;
    /**
     * Energy STATE.
     */
    private STATE state = STATE.BOTH;
    /**
     * True when a Monte Carlo move is accepted.
     */
    private boolean moveAccepted = false;
    /**
     * Current density value.
     */
    private double currentDensity = 0;
    /**
     * Mean density value.
     */
    private double densityMean = 0;
    /**
     * Mean density squared for STD calculation.
     */
    private double densityMean2 = 0;
    /**
     * Standard deviation for density.
     */
    private double densitySD = 0;
    /**
     * Current lattice parameters.
     */
    private double a = 0;
    private double b = 0;
    private double c = 0;
    private double alpha = 0;
    private double beta = 0;
    private double gamma = 0;
    /**
     * Mean values of lattice parameters.
     */
    private double aMean = 0;
    private double bMean = 0;
    private double cMean = 0;
    private double alphaMean = 0;
    private double betaMean = 0;
    private double gammaMean = 0;
    /**
     * Mean values squared of lattice parameters for STD calculation.
     */
    private double aMean2 = 0;
    private double bMean2 = 0;
    private double cMean2 = 0;
    private double alphaMean2 = 0;
    private double betaMean2 = 0;
    private double gammaMean2 = 0;
    /**
     * Standard deviation of lattice parameters.
     */
    private double aSD = 0;
    private double bSD = 0;
    private double cSD = 0;
    private double alphaSD = 0;
    private double betaSD = 0;
    private double gammaSD = 0;
    private final int printFrequency = 1000;
    private double minDensity = 0.75;
    private double maxDensity = 1.50;
    private MoveType moveType = MoveType.SIDE;

    /**
     * Initialize the Barostat.
     *
     * @param molecularAssembly The molecular assembly to apply the MC barostat
     * to.
     */
    public Barostat(MolecularAssembly molecularAssembly) {

        this.molecularAssembly = molecularAssembly;
        potential = molecularAssembly.getPotentialEnergy();
        crystal = potential.getCrystal();
        unitCell = crystal.getUnitCell();
        spaceGroup = unitCell.spaceGroup;
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        nSymm = spaceGroup.getNumberOfSymOps();

        mass = 0.0;
        for (int i = 0; i < nAtoms; i++) {
            mass += atoms[i].getMass();
        }

        nMolecules = countMolecules();
        logger.info(String.format(" There are %d molecules.", nMolecules));
        fractionalCOM = new double[nMolecules][3];
    }

    @Override
    public double energy(double[] x) {
        // Do not apply the Barostat for energy only evaluations.
        return potential.energy(x);
    }

    public void setMeanBarostatInterval(int meanBarostatInterval) {
        this.meanBarostatInterval = meanBarostatInterval;
    }

    public void setMinDensity(double minDensity) {
        this.minDensity = minDensity;
    }

    public void setMaxDensity(double maxDensity) {
        this.maxDensity = maxDensity;
    }

    public void setMaxAngleMove(double maxAngleMove) {
        this.maxAngleMove = maxAngleMove;
    }

    public void setMaxSideMove(double maxSideMove) {
        this.maxSideMove = maxSideMove;
    }

    public void setPressure(double pressure) {
        this.pressure = pressure;
    }

    private double mcStep(double currentE, double currentV) {

        switch (moveType) {
            case SIDE:
                sideMovesAttempted++;
                break;
            case ANGLE:
                angleMovesAttempted++;
                break;
            default:
                break;
        }

        /**
         * Enforce minimum & maximum density constraints.
         */
        double den = density();
        if (den < minDensity || den > maxDensity) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(
                        " MC Density %10.6f is outside the range %10.6f - %10.6f.",
                        den, minDensity, maxDensity));
            }
            // Fail moves outside the specified density range.
            crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
            return currentE;
        }

        /**
         * Enforce minimum interfacial radii of 1.2 Angstroms.
         */
        if (unitCell.interfacialRadiusA < minInterfacialRadius ||
                unitCell.interfacialRadiusB < minInterfacialRadius ||
                unitCell.interfacialRadiusC < minInterfacialRadius) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(
                        " MC An interfacial radius (%10.6f,%10.6f,%10.6f) is below the minimium %10.6f",
                        unitCell.interfacialRadiusA,
                        unitCell.interfacialRadiusB,
                        unitCell.interfacialRadiusC,
                        minInterfacialRadius));
            }
            // Fail small axis length trial moves.
            crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);
            return currentE;
        }

        /**
         * Apply the boundary condition for the proposed move. If the move is
         * rejected, then the previous boundary conditions must be restored.
         */
        potential.setCrystal(crystal);

        // Update atomic coordinates to maintain molecular fractional centers of mass.
        moveToFractionalCOM();

        // Save the new volume
        double newV = unitCell.volume / nSymm;

        // Compute the new energy
        double newE = potential.energy(false, false);

        // Compute the change in potential energy
        double dE = newE - currentE;

        // Compute the pressure-volume work for the asymmetric unit.
        double dV = pressure * (newV - currentV) / PRESCON;
        // double dV = 0.0;

        // Compute the volume entropy
        double dS = -nMolecules * kT * Math.log(newV / currentV);
        //double dS = 0.0;

        // Add up the contributions
        double dT = dE + dV + dS;

        if (dT < 0.0) {
            // Energy decreased; accept the move.
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(format(
                        " MC Energy: %12.6f (dE) + %12.6f (dV) + %12.6f (dS) = %12.6f with (V=%12.6f, E=%12.6f)",
                        dE, dV, dS, dT, newV, newE));
            }
            moveAccepted = true;
            switch (moveType) {
                case SIDE:
                    sideMovesAccepted++;
                    break;
                case ANGLE:
                    angleMovesAccepted++;
                    break;
                default:
                    break;
            }
            return newE;
        }

        // Apply the Metropolis criteria.
        double boltzmann = exp(-dT / kT);
        double metropolis = random();

        // Energy increase without Metropolis criteria satisified.
        if (metropolis > boltzmann) {

            rejectMove();

            if (logger.isLoggable(Level.FINE)) {
                logger.fine(
                        format(" MC Reject: %12.6f (dE) + %12.6f (dV) + %12.6f (dS) = %12.6f with (V=%12.6f, E=%12.6f)",
                                dE, dV, dS, dT, currentV, currentE));
            }
            return currentE;
        }

        // Energy increase with Metropolis criteria satisified.
        if (logger.isLoggable(Level.FINE)) {
            logger.fine(
                    format(" MC Accept: %12.6f (dE) + %12.6f (dV) + %12.6f (dS) = %12.6f with (V=%12.6f, E=%12.6f)",
                            dE, dV, dS, dT, newV, newE));
        }

        moveAccepted = true;
        switch (moveType) {
            case SIDE:
                sideMovesAccepted++;
                break;
            case ANGLE:
                angleMovesAccepted++;
                break;
            default:
                break;
        }

        return newE;
    }

    /**
     * Reset the to state prior to trial move.
     */
    private void rejectMove() {
        // Reset the unit cell parameters
        crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma);

        // Reset the potential PBC.
        potential.setCrystal(crystal);

        // Reset the atomic coordinates to maintain molecular fractional centers of mass.
        moveToFractionalCOM();
    }

    public double density() {
        return (mass * nSymm / AVOGADRO) * (1.0e24 / unitCell.volume);
    }

    private double mcA(double currentE) {
        moveType = MoveType.SIDE;
        double currentV = unitCell.volume / nSymm;
        double move = maxSideMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a + move, b, c, alpha, beta, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the a-axis (%6.3f) of %6.3f A", a, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcB(double currentE) {
        moveType = MoveType.SIDE;
        double currentV = unitCell.volume / nSymm;
        double move = maxSideMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b + move, c, alpha, beta, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the b-axis (%6.3f) of %6.3f A", b, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcC(double currentE) {
        moveType = MoveType.SIDE;
        double currentV = unitCell.volume / nSymm;
        double move = maxSideMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b, c + move, alpha, beta, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the c-axis (%6.3f) of %6.3f A", c, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcAB(double currentE) {
        moveType = MoveType.SIDE;
        double currentV = unitCell.volume / nSymm;
        double move = maxSideMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a + move, b + move, c, alpha, beta, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the a,b-axis (%6.3f) of %6.3f A", a, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcABC(double currentE) {
        moveType = MoveType.SIDE;
        double currentV = unitCell.volume / nSymm;
        double move = maxSideMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a + move, b + move, c + move, alpha, beta, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the a,b,c-axis (%6.3f) of %6.3f A", a, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcAlpha(double currentE) {
        moveType = MoveType.ANGLE;
        double currentV = unitCell.volume / nSymm;
        double move = maxAngleMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b, c, alpha + move, beta, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the alpha angle (%6.3f) of %6.3f (degrees)", alpha, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcBeta(double currentE) {
        moveType = MoveType.ANGLE;
        double currentV = unitCell.volume / nSymm;
        double move = maxAngleMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b, c, alpha, beta + move, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the beta angle (%6.3f) of %6.3f (degrees)", beta, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcGamma(double currentE) {
        moveType = MoveType.ANGLE;
        double currentV = unitCell.volume / nSymm;
        double move = maxAngleMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b, c, alpha, beta, gamma + move);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the gamma angle (%6.3f) of %6.3f (degrees)", gamma, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcABeta(double currentE) {
        moveType = MoveType.ANGLE;
        double currentV = unitCell.volume / nSymm;
        double move = maxAngleMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b, c, alpha + move, beta + move, gamma);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the alpha/beta angles (%6.3f) of %6.3f (degrees)", alpha, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcAG(double currentE) {
        moveType = MoveType.ANGLE;
        double currentV = unitCell.volume / nSymm;
        double move = maxAngleMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b, c, alpha + move, beta, gamma + move);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the alpha/gamma angles (%6.3f) of %6.3f (degrees)", alpha, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private double mcABG(double currentE) {
        moveType = MoveType.ANGLE;
        double currentV = unitCell.volume / nSymm;
        double move = maxAngleMove * (2.0 * Math.random() - 1.0);
        boolean succeed = crystal.changeUnitCellParameters(a, b, c, alpha + move, beta + move, gamma + move);
        if (succeed) {
            if (logger.isLoggable(Level.FINE)) {
                logger.fine(String.format(" Propsing MC change to the alpha/beta/gamma angles (%6.3f) of %6.3f (degrees)", alpha, move));
            }
            return mcStep(currentE, currentV);
        }
        return currentE;
    }

    private int countMolecules() {
        int count = 0;
        // Move polymers togethers.
        Polymer polymers[] = molecularAssembly.getChains();
        if (polymers != null && polymers.length > 0) {
            count += polymers.length;
        }
        List<Molecule> molecules = molecularAssembly.getMolecules();
        if (molecules != null) {
            count += molecules.size();
        }
        List<MSNode> waters = molecularAssembly.getWaters();
        if (waters != null) {
            count += waters.size();
        }
        List<MSNode> ions = molecularAssembly.getIons();
        if (ions != null) {
            count += ions.size();
        }
        return count;
    }

    public void setDensity(double density) {
        computeFractionalCOM();

        crystal.setDensity(density, mass);

        potential.setCrystal(crystal);

        moveToFractionalCOM();
    }

    private void computeFractionalCOM() {

        int iMolecule = 0;
        double[] com = new double[3];

        Polymer polymers[] = molecularAssembly.getChains();
        if (polymers != null && polymers.length > 0) {
            // Find the center of mass
            for (Polymer polymer : polymers) {
                List<Atom> list = polymer.getAtomList();
                com[0] = 0.0;
                com[1] = 0.0;
                com[2] = 0.0;
                double totalMass = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    com[0] += atom.getX() * m;
                    com[1] += atom.getY() * m;
                    com[2] += atom.getZ() * m;
                    totalMass += m;
                }
                com[0] /= totalMass;
                com[1] /= totalMass;
                com[2] /= totalMass;
                unitCell.toFractionalCoordinates(com, fractionalCOM[iMolecule++]);
            }
        }

        // Loop over each molecule
        List<Molecule> molecules = molecularAssembly.getMolecules();
        for (MSNode molecule : molecules) {
            List<Atom> list = molecule.getAtomList();
            // Find the center of mass
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            double totalMass = 0.0;
            for (Atom atom : list) {
                double m = atom.getMass();
                com[0] += atom.getX() * m;
                com[1] += atom.getY() * m;
                com[2] += atom.getZ() * m;
                totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            unitCell.toFractionalCoordinates(com, fractionalCOM[iMolecule++]);
        }

        // Loop over each water
        List<MSNode> waters = molecularAssembly.getWaters();
        for (MSNode water : waters) {
            List<Atom> list = water.getAtomList();
            // Find the center of mass
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            double totalMass = 0.0;
            for (Atom atom : list) {
                double m = atom.getMass();
                com[0] += atom.getX() * m;
                com[1] += atom.getY() * m;
                com[2] += atom.getZ() * m;
                totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            unitCell.toFractionalCoordinates(com, fractionalCOM[iMolecule++]);
        }

        // Loop over each ion
        List<MSNode> ions = molecularAssembly.getIons();
        for (MSNode ion : ions) {
            List<Atom> list = ion.getAtomList();
            // Find the center of mass
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            double totalMass = 0.0;
            for (Atom atom : list) {
                double m = atom.getMass();
                com[0] += atom.getX() * m;
                com[1] += atom.getY() * m;
                com[2] += atom.getZ() * m;
                totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            unitCell.toFractionalCoordinates(com, fractionalCOM[iMolecule++]);
        }
    }

    private void moveToFractionalCOM() {
        int iMolecule = 0;
        double[] com = new double[3];

        Polymer polymers[] = molecularAssembly.getChains();
        if (polymers != null && polymers.length > 0) {
            // Find the center of mass
            for (Polymer polymer : polymers) {
                List<Atom> list = polymer.getAtomList();
                double totalMass = 0.9;
                com[0] = 0.0;
                com[1] = 0.0;
                com[2] = 0.0;
                for (Atom atom : list) {
                    double m = atom.getMass();
                    com[0] += atom.getX() * m;
                    com[1] += atom.getY() * m;
                    com[2] += atom.getZ() * m;
                    totalMass += m;
                }
                com[0] /= totalMass;
                com[1] /= totalMass;
                com[2] /= totalMass;
                // Find the new center of mass in fractional coordinates.
                unitCell.toFractionalCoordinates(com, com);
                // Find the reciprocal translation vector.
                double[] frac = fractionalCOM[iMolecule++];
                com[0] = frac[0] - com[0];
                com[1] = frac[1] - com[1];
                com[2] = frac[2] - com[2];
                // Convert the fractional translation vector to Cartesian coordinates.
                unitCell.toCartesianCoordinates(com, com);
                // Move all atoms.
//            for (Polymer polymer : polymers) {
//                List<Atom> list = polymer.getAtomList();
                for (Atom atom : list) {
                    atom.move(com);
                }
            }
        }

        // Loop over each molecule
        List<Molecule> molecules = molecularAssembly.getMolecules();
        for (MSNode molecule : molecules) {
            List<Atom> list = molecule.getAtomList();
            // Find the center of mass
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            double totalMass = 0.0;
            for (Atom atom : list) {
                double m = atom.getMass();
                com[0] += atom.getX() * m;
                com[1] += atom.getY() * m;
                com[2] += atom.getZ() * m;
                totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            // Find the new center of mass in fractional coordinates.
            unitCell.toFractionalCoordinates(com, com);
            // Find the reciprocal translation vector to the previous COM.
            double[] frac = fractionalCOM[iMolecule++];
            com[0] = frac[0] - com[0];
            com[1] = frac[1] - com[1];
            com[2] = frac[2] - com[2];
            // Convert the fractional translation vector to Cartesian coordinates.
            unitCell.toCartesianCoordinates(com, com);
            // Move all atoms.
            for (Atom atom : list) {
                atom.move(com);
            }
        }

        // Loop over each water
        List<MSNode> waters = molecularAssembly.getWaters();
        for (MSNode water : waters) {
            List<Atom> list = water.getAtomList();
            // Find the center of mass
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            double totalMass = 0.0;
            for (Atom atom : list) {
                double m = atom.getMass();
                com[0] += atom.getX() * m;
                com[1] += atom.getY() * m;
                com[2] += atom.getZ() * m;
                totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            // Find the new center of mass in fractional coordinates.
            unitCell.toFractionalCoordinates(com, com);
            // Find the reciprocal translation vector to the previous COM.
            double[] frac = fractionalCOM[iMolecule++];
            com[0] = frac[0] - com[0];
            com[1] = frac[1] - com[1];
            com[2] = frac[2] - com[2];
            // Convert the fractional translation vector to Cartesian coordinates.
            unitCell.toCartesianCoordinates(com, com);

            double r = ffx.numerics.VectorMath.r(com);
            /**
             * Warn if an atom is moved more than 1 Angstrom.
             */
            if (r > 1.0) {
                int i = iMolecule - 1;
                logger.info(String.format(" %d R: %16.8f", i, r));
                logger.info(String.format(" %d FRAC %16.8f %16.8f %16.8f", i, frac[0], frac[1], frac[2]));
                logger.info(String.format(" %d COM  %16.8f %16.8f %16.8f", i, com[0], com[1], com[2]));
            }

            // Move all atoms.
            for (Atom atom : list) {
                atom.move(com);
            }
        }

        // Loop over each ion
        List<MSNode> ions = molecularAssembly.getIons();
        for (MSNode ion : ions) {
            List<Atom> list = ion.getAtomList();
            // Find the center of mass
            com[0] = 0.0;
            com[1] = 0.0;
            com[2] = 0.0;
            double totalMass = 0.0;
            for (Atom atom : list) {
                double m = atom.getMass();
                com[0] += atom.getX() * m;
                com[1] += atom.getY() * m;
                com[2] += atom.getZ() * m;
                totalMass += m;
            }
            com[0] /= totalMass;
            com[1] /= totalMass;
            com[2] /= totalMass;
            // Find the new center of mass in fractional coordinates.
            unitCell.toFractionalCoordinates(com, com);
            // Find the reciprocal translation vector to the previous COM.
            double[] frac = fractionalCOM[iMolecule++];
            com[0] = frac[0] - com[0];
            com[1] = frac[1] - com[1];
            com[2] = frac[2] - com[2];
            // Convert the fractional translation vector to Cartesian coordinates.
            unitCell.toCartesianCoordinates(com, com);
            // Move all atoms.
            for (Atom atom : list) {
                atom.move(com);
            }
        }
    }

    private double applyBarostat(double currentE) {
        // Determine the current molecular centers of mass in fractional coordinates.
        computeFractionalCOM();

        // Collect the current unit cell parameters.
        crystal = potential.getCrystal();
        unitCell = crystal.getUnitCell();
        a = unitCell.a;
        b = unitCell.b;
        c = unitCell.c;
        alpha = unitCell.alpha;
        beta = unitCell.beta;
        gamma = unitCell.gamma;

        switch (spaceGroup.crystalSystem) {
            case MONOCLINIC: {
                int move = (int) floor(random() * 4.0);
                switch (move) {
                    case 0:
                        currentE = mcA(currentE);
                        break;
                    case 1:
                        currentE = mcB(currentE);
                        break;
                    case 2:
                        currentE = mcC(currentE);
                        break;
                    case 3:
                        currentE = mcBeta(currentE);
                        break;
                    default:
                        logger.severe(" Barostat programming error.");
                }
                break;
            }
            case ORTHORHOMBIC: {
                // alpha == beta == gamma == 90.0
                int move = (int) floor(random() * 3.0);
                switch (move) {
                    case 0:
                        currentE = mcA(currentE);
                        break;
                    case 1:
                        currentE = mcB(currentE);
                        break;
                    case 2:
                        currentE = mcC(currentE);
                        break;
                    default:
                        logger.severe(" Barostat programming error.");
                }
                break;
            }
            case TETRAGONAL: {
                // (a == b, alpha == beta == gamma == 90.0
                int move = (int) floor(random() * 2.0);
                switch (move) {
                    case 0:
                        currentE = mcAB(currentE);
                        break;
                    case 1:
                        currentE = mcC(currentE);
                        break;
                    default:
                        logger.severe(" Barostat programming error.");
                }
                break;
            }
            case TRIGONAL: {
                if (a == b && b == c && alpha == beta && beta == gamma) {
                    // Rombohedral axes, primitive cell.
                    int move = (int) floor(random() * 2.0);
                    switch (move) {
                        case 0:
                            currentE = mcABC(currentE);
                            break;
                        case 1:
                            currentE = mcABG(currentE);
                            break;
                        default:
                            logger.severe(" Barostat programming error.");
                    }
                } else if (a == b && alpha == 90.0 && beta == 90.0 && gamma == 120.0) {
                    // Hexagonal axes, triple obverse cell.
                    int move = (int) floor(random() * 2.0);
                    switch (move) {
                        case 0:
                            currentE = mcAB(currentE);
                            break;
                        case 1:
                            currentE = mcC(currentE);
                            break;
                        default:
                            logger.severe(" Barostat programming error.");
                    }
                } else {
                    logger.warning(" Trigonal constraints not satisfied.");
                }
                break;
            }
            case HEXAGONAL: {
                // a == b, alpha == beta == 90.0, gamma == 120.0
                int move = (int) floor(random() * 2.0);
                switch (move) {
                    case 0:
                        currentE = mcAB(currentE);
                        break;
                    case 1:
                        currentE = mcC(currentE);
                        break;
                    default:
                        logger.severe(" Barostat programming error.");
                }
                break;
            }
            case CUBIC:
                // a == b == c, alpha == beta == gamma == 90.0
                currentE = mcABC(currentE);
                break;
            case TRICLINIC:
            default: {
                if (a == b && b == c && alpha == 90.0 && beta == 90.0 && gamma == 90.0) {
                    currentE = mcABC(currentE);
                } else {
                    int move = (int) floor(random() * 6.0);
                    switch (move) {
                        case 0:
                            currentE = mcA(currentE);
                            break;
                        case 1:
                            currentE = mcB(currentE);
                            break;
                        case 2:
                            currentE = mcC(currentE);
                            break;
                        case 3:
                            currentE = mcAlpha(currentE);
                            break;
                        case 4:
                            currentE = mcBeta(currentE);
                            break;
                        case 5:
                            currentE = mcGamma(currentE);
                            break;
                        default:
                            logger.severe(" Barostat programming error.");
                    }
                }
            }
        }

        currentDensity = density();
        if (moveAccepted) {
            if (angleMovesAttempted > 0) {
                logger.info(format(" Density: %5.3f UC: %s MCS: %5.1f MCA: %5.1f", currentDensity, unitCell.toShortString(),
                        (double) sideMovesAccepted / sideMovesAttempted * 100.0,
                        (double) angleMovesAccepted / angleMovesAttempted * 100.0));
            } else {
                logger.info(format(" Density: %5.3f UC: %s MCS: %5.1f", currentDensity, unitCell.toShortString(),
                        (double) sideMovesAccepted / sideMovesAttempted * 100.0));
            }
        } else {
            // Check that the unit cell parameters have not changed.
            if (unitCell.a != a || unitCell.b != b || unitCell.c != c
                    || unitCell.alpha != alpha || unitCell.beta != beta || unitCell.gamma != gamma) {
                logger.severe(" Reversion of unit cell parameters did not succeed after failed Barostat MC move.");
            }

        }

        return currentE;
    }

    public void setActive(boolean active) {
        this.active = active;
    }

    /**
     * {@inheritDoc}
     *
     * @param x
     * @param g
     * @return
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        /**
         * Calculate the energy and gradient as usual.
         */
        double energy = potential.energyAndGradient(x, g);
        /**
         * Apply the barostat during computation of slowly varying forces.
         */
        if (active && state != STATE.FAST) {
            if (Math.random() < (1.0 / meanBarostatInterval)) {

                // Attempt to change the unit cell parameters.
                moveAccepted = false;
                applyBarostat(energy);

                // Collect Statistics.
                collectStats();

                /**
                 * If a move was accepted, then re-calculate the gradient so
                 * that it's consistent with the current unit cell parameters.
                 */
                if (moveAccepted) {
                    energy = potential.energyAndGradient(x, g);
                }
            }
        }
        return energy;
    }

    private void collectStats() {
        // Collect statistics.
        barostatCount++;
        densityMean += currentDensity;
        densityMean2 += currentDensity * currentDensity;
        aMean += a;
        bMean += b;
        cMean += c;
        alphaMean += alpha;
        betaMean += beta;
        gammaMean += gamma;
        aMean2 += a * a;
        bMean2 += b * b;
        cMean2 += c * c;
        alphaMean2 += alpha * alpha;
        betaMean2 += beta * beta;
        gammaMean2 += gamma * gamma;
        if (barostatCount % printFrequency == 0) {
            densityMean = densityMean / printFrequency;
            densityMean2 = densityMean2 / printFrequency;
            densitySD = sqrt(densityMean2 - densityMean * densityMean);
            logger.info(String.format(" Density: %5.3f +/- %5.3f", densityMean, densitySD));
            aMean = aMean / printFrequency;
            bMean = bMean / printFrequency;
            cMean = cMean / printFrequency;
            alphaMean = alphaMean / printFrequency;
            betaMean = betaMean / printFrequency;
            gammaMean = gammaMean / printFrequency;
            aMean2 = aMean2 / printFrequency;
            bMean2 = bMean2 / printFrequency;
            cMean2 = cMean2 / printFrequency;
            alphaMean2 = alphaMean2 / printFrequency;
            betaMean2 = betaMean2 / printFrequency;
            gammaMean2 = gammaMean2 / printFrequency;
            aSD = sqrt(aMean2 - aMean * aMean);
            bSD = sqrt(bMean2 - bMean * bMean);
            cSD = sqrt(cMean2 - cMean * cMean);
            alphaSD = sqrt(alphaMean2 - alphaMean * alphaMean);
            betaSD = sqrt(betaMean2 - betaMean * betaMean);
            gammaSD = sqrt(gammaMean2 - gammaMean * gammaMean);
            logger.info(String.format(
                    " Lattice a: %4.2f +/-%4.2f b: %4.2f +/-%4.2f c: %4.2f +/-%4.2f alpha: %5.2f +/-%3.2f beta: %5.2f +/-%3.2f gamma: %5.2f +/-%3.2f",
                    aMean, aSD, bMean, bSD, cMean, cSD, alphaMean, alphaSD, betaMean, betaSD, gammaMean, gammaSD));
            densityMean = 0;
            densityMean2 = 0;
            aMean = 0;
            bMean = 0;
            cMean = 0;
            alphaMean = 0;
            betaMean = 0;
            gammaMean = 0;
            aMean2 = 0;
            bMean2 = 0;
            cMean2 = 0;
            alphaMean2 = 0;
            betaMean2 = 0;
            gammaMean2 = 0;
        }
    }

    /**
     * {@inheritDoc}
     *
     * @param scaling
     */
    @Override
    public void setScaling(double[] scaling) {
        potential.setScaling(scaling);
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public double[] getScaling() {
        return potential.getScaling();
    }

    /**
     * {@inheritDoc}
     *
     * @param parameters
     * @return
     */
    @Override
    public double[] getCoordinates(double[] parameters) {
        return potential.getCoordinates(parameters);
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public double[] getMass() {
        return potential.getMass();
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public double getTotalEnergy() {
        return potential.getTotalEnergy();
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public int getNumberOfVariables() {
        return potential.getNumberOfVariables();
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return potential.getVariableTypes();
    }

    /**
     * {@inheritDoc}
     *
     * @param state
     */
    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
        potential.setEnergyTermState(state);

    }

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    @Override
    public void setVelocity(double[] velocity) {
        potential.setVelocity(velocity);
    }

    @Override
    public void setAcceleration(double[] acceleration) {
        potential.setAcceleration(acceleration);
    }

    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        potential.setPreviousAcceleration(previousAcceleration);
    }

    @Override
    public double[] getVelocity(double[] velocity) {
        return potential.getVelocity(velocity);
    }

    @Override
    public double[] getAcceleration(double[] acceleration) {
        return potential.getAcceleration(acceleration);
    }

    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        return potential.getPreviousAcceleration(previousAcceleration);
    }

    @Override
    public void setLambda(double lambda) {
        potential.setLambda(lambda);
    }

    @Override
    public double getLambda() {
        return potential.getLambda();
    }

    @Override
    public double getdEdL() {
        return potential.getdEdL();
    }

    @Override
    public double getd2EdL2() {
        return potential.getd2EdL2();
    }

    @Override
    public void getdEdXdL(double[] gradient) {
        potential.getdEdXdL(gradient);
    }

    private enum MoveType {

        SIDE, ANGLE
    }
}
