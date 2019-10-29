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
package ffx.algorithms.dynamics;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.floor;
import static org.apache.commons.math3.util.FastMath.log;
import static org.apache.commons.math3.util.FastMath.random;
import static org.apache.commons.math3.util.FastMath.sqrt;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.crystal.SpaceGroup;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import static ffx.utilities.Constants.AVOGADRO;
import static ffx.utilities.Constants.PRESCON;

/**
 * The Barostat class maintains constant pressure using random trial moves in
 * lattice parameters, which are consistent with the space group.
 * <p>
 * Frenkel and Smit, "Understanding Molecular Simulation, 2nd
 * Edition", Academic Press, San Diego, CA, 2002; Section 5.4
 *
 * @author Michael J. Schnieders
 */
public class Barostat implements CrystalPotential {

    private static final Logger logger = Logger.getLogger(Barostat.class.getName());
    /**
     * Ideal gas constant in kcal/mol/K
     */
    private static final double kB = 1.98720415e-3;
    /**
     * Ideal gas constant * temperature (kcal/mol).
     */
    private double kT;
    /**
     * Sampling pressure (atm)
     */
    private double pressure = 1.0;
    /**
     * Flag to turn the Barostat on or off. If false, MC moves will not be tried.
     */
    private boolean active = true;
    /**
     * Default edge length move (A).
     */
    private double maxSideMove = 0.25;
    /**
     * Default angle move (degrees).
     */
    private double maxAngleMove = 0.5;
    /**
     * MolecularAssembly being simulated.
     */
    private final MolecularAssembly molecularAssembly;
    /**
     * Mass of the system.
     */
    private final double mass;
    /**
     * ForceFieldEnergy that describes the system.
     */
    private final CrystalPotential potential;
    /**
     * Atomic coordinates.
     */
    private final double[] x;
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
    private double minDensity = 0.75;
    // NNQQ Peptide has a density of 1.4.
    private double maxDensity = 1.60;
    private MoveType moveType = MoveType.SIDE;

    /**
     * Initialize the Barostat.
     *
     * @param molecularAssembly The molecular assembly to apply the MC barostat to.
     * @param potential         a {@link ffx.crystal.CrystalPotential} object.
     */
    public Barostat(MolecularAssembly molecularAssembly, CrystalPotential potential) {
        this(molecularAssembly, potential, 298.15);
    }

    /**
     * Initialize the Barostat.
     *
     * @param molecularAssembly The molecular assembly to apply the MC barostat to.
     * @param potential         a {@link ffx.crystal.CrystalPotential} object.
     * @param temperature       The Metropolis Monte Carlo temperature (Kelvin).
     */
    public Barostat(MolecularAssembly molecularAssembly, CrystalPotential potential, double temperature) {

        this.molecularAssembly = molecularAssembly;
        this.potential = potential;
        this.kT = temperature * kB;

        crystal = potential.getCrystal();
        unitCell = crystal.getUnitCell();
        spaceGroup = unitCell.spaceGroup;
        // Atoms in the system.
        Atom[] atoms = molecularAssembly.getAtomArray();
        // Number of atoms.
        int nAtoms = atoms.length;
        nSymm = spaceGroup.getNumberOfSymOps();
        mass = molecularAssembly.getMass();
        x = new double[3 * nAtoms];
        nMolecules = molecularAssembly.fractionalCount();
        logger.info(String.format(" There are %d molecules.", nMolecules));
    }

    /**
     * Set the Metropolis Monte Carlo temperature.
     *
     * @param temperature Temperature (Kelvin).
     */
    public void setTemperature(double temperature) {
        this.kT = temperature * kB;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x) {
        // Do not apply the Barostat for energy only evaluations.
        return potential.energy(x);
    }

    /**
     * <p>Setter for the field <code>meanBarostatInterval</code>.</p>
     *
     * @param meanBarostatInterval a int.
     */
    public void setMeanBarostatInterval(int meanBarostatInterval) {
        this.meanBarostatInterval = meanBarostatInterval;
    }

    /**
     * Returns the mean number of steps between barostat applications.
     *
     * @return Mean steps between barostat applications.
     */
    public int getMeanBarostatInterval() {
        return meanBarostatInterval;
    }

    /**
     * <p>Setter for the field <code>minDensity</code>.</p>
     *
     * @param minDensity a double.
     */
    public void setMinDensity(double minDensity) {
        this.minDensity = minDensity;
    }

    /**
     * <p>Setter for the field <code>maxDensity</code>.</p>
     *
     * @param maxDensity a double.
     */
    public void setMaxDensity(double maxDensity) {
        this.maxDensity = maxDensity;
    }

    /**
     * <p>Setter for the field <code>maxAngleMove</code>.</p>
     *
     * @param maxAngleMove a double.
     */
    public void setMaxAngleMove(double maxAngleMove) {
        this.maxAngleMove = maxAngleMove;
    }

    /**
     * <p>Setter for the field <code>maxSideMove</code>.</p>
     *
     * @param maxSideMove a double.
     */
    public void setMaxSideMove(double maxSideMove) {
        this.maxSideMove = maxSideMove;
    }

    /**
     * <p>Setter for the field <code>pressure</code>.</p>
     *
     * @param pressure a double.
     */
    public void setPressure(double pressure) {
        this.pressure = pressure;
    }

    /**
     * Gets the pressure of this Barostat in atm.
     * @return Pressure in atm.
     */
    public double getPressure() {
        return pressure;
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

        // Enforce minimum & maximum density constraints.
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


        // A carbon atom cannot fit into a unit cell without interfacial radii greater than ~1.2 Angstroms.
        double minInterfacialRadius = 1.2;
        // Enforce minimum interfacial radii of 1.2 Angstroms.
        if (unitCell.interfacialRadiusA < minInterfacialRadius
                || unitCell.interfacialRadiusB < minInterfacialRadius
                || unitCell.interfacialRadiusC < minInterfacialRadius) {
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

        // Apply the boundary condition for the proposed move. If the move is
        // rejected, then the previous boundary conditions must be restored.
        potential.setCrystal(crystal);

        // Update atomic coordinates to maintain molecular fractional centers of mass.
        molecularAssembly.moveToFractionalCoordinates();

        // Save the new volume
        double newV = unitCell.volume / nSymm;

        potential.getCoordinates(x);

        // Compute the new energy
        double newE = potential.energy(x);

        // Compute the change in potential energy
        double dE = newE - currentE;

        // Compute the pressure-volume work for the asymmetric unit.
        double dV = pressure * (newV - currentV) / PRESCON;
        // double dV = 0.0;

        // Compute the volume entropy
        double dS = -nMolecules * kT * log(newV / currentV);
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
        molecularAssembly.moveToFractionalCoordinates();
    }

    /**
     * <p>density.</p>
     *
     * @return a double.
     */
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

    /**
     * <p>setDensity.</p>
     *
     * @param density a double.
     */
    public void setDensity(double density) {
        molecularAssembly.computeFractionalCoordinates();
        crystal.setDensity(density, mass);
        potential.setCrystal(crystal);
        molecularAssembly.moveToFractionalCoordinates();
    }

    private double applyBarostat(double currentE) {
        // Determine the current molecular centers of mass in fractional coordinates.
        molecularAssembly.computeFractionalCoordinates();

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

    /**
     * <p>Setter for the field <code>active</code>.</p>
     *
     * @param active a boolean.
     */
    public void setActive(boolean active) {
        this.active = active;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {

        // Calculate the energy and gradient as usual.
        double energy = potential.energyAndGradient(x, g);

        // Apply the barostat during computation of slowly varying forces.
        if (active && state != STATE.FAST) {
            if (random() < (1.0 / meanBarostatInterval)) {

                // Attempt to change the unit cell parameters.
                moveAccepted = false;

                applyBarostat(energy);

                // Collect Statistics.
                collectStats();

                // If a move was accepted, then re-calculate the gradient so
                // that it's consistent with the current unit cell parameters.
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
        int printFrequency = 1000;
        if (barostatCount % printFrequency == 0) {
            densityMean = densityMean / printFrequency;
            densityMean2 = densityMean2 / printFrequency;
            double densitySD = sqrt(densityMean2 - densityMean * densityMean);
            logger.info(format(" Density: %5.3f +/- %5.3f", densityMean, densitySD));
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
            double aSD = sqrt(aMean2 - aMean * aMean);
            double bSD = sqrt(bMean2 - bMean * bMean);
            double cSD = sqrt(cMean2 - cMean * cMean);
            double alphaSD = sqrt(alphaMean2 - alphaMean * alphaMean);
            double betaSD = sqrt(betaMean2 - betaMean * betaMean);
            double gammaSD = sqrt(gammaMean2 - gammaMean * gammaMean);
            logger.info(format(" Lattice a: %4.2f +/-%4.2f b: %4.2f +/-%4.2f c: %4.2f +/-%4.2f alpha: %5.2f +/-%3.2f beta: %5.2f +/-%3.2f gamma: %5.2f +/-%3.2f",
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
     */
    @Override
    public void setScaling(double[] scaling) {
        potential.setScaling(scaling);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getScaling() {
        return potential.getScaling();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double[] parameters) {
        return potential.getCoordinates(parameters);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getMass() {
        return potential.getMass();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getTotalEnergy() {
        return potential.getTotalEnergy();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        return potential.getNumberOfVariables();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        return potential.getVariableTypes();
    }

    @Override
    public List<Potential> getUnderlyingPotentials() {
        List<Potential> underlying = new ArrayList<>();
        underlying.add(potential);
        underlying.addAll(potential.getUnderlyingPotentials());
        return underlying;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
        potential.setEnergyTermState(state);

    }

    /**
     * {@inheritDoc}
     */
    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setVelocity(double[] velocity) {
        potential.setVelocity(velocity);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setAcceleration(double[] acceleration) {
        potential.setAcceleration(acceleration);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        potential.setPreviousAcceleration(previousAcceleration);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getVelocity(double[] velocity) {
        return potential.getVelocity(velocity);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getAcceleration(double[] acceleration) {
        return potential.getAcceleration(acceleration);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        return potential.getPreviousAcceleration(previousAcceleration);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Crystal getCrystal() {
        return potential.getCrystal();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setCrystal(Crystal crystal) {
        potential.setCrystal(crystal);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean destroy() {
        // Nothing at this level to destroy.
        return potential.destroy();
    }

    private enum MoveType {

        SIDE, ANGLE
    }
}
