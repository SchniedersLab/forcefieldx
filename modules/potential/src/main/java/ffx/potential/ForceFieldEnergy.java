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
package ffx.potential;

import java.io.File;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static java.util.Arrays.fill;
import static java.util.Arrays.sort;

import ffx.numerics.Constraint;
import ffx.potential.constraint.SettleConstraint;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import edu.rit.pj.ParallelTeam;
import edu.rit.pj.reduction.SharedDouble;

import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.crystal.NCSCrystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SymOp;
import ffx.numerics.atomic.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.atomic.AtomicDoubleArray3D;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.RelativeSolvation;
import ffx.potential.bonded.RelativeSolvation.SolvationLibrary;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.StretchTorsion;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.COMRestraint;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.NCSRestraint;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.ELEC_FORM;
import ffx.potential.nonbonded.ParticleMeshEwaldCart;
import ffx.potential.nonbonded.ParticleMeshEwaldQI;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import static ffx.potential.parameters.ForceField.toEnumForm;

/**
 * Compute the potential energy and derivatives of an AMOEBA system.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class ForceFieldEnergy implements CrystalPotential, LambdaInterface {

    /**
     * A Logger for the ForceFieldEnergy class.
     */
    private static final Logger logger = Logger.getLogger(ForceFieldEnergy.class.getName());
    /**
     * Convert from nanoseconds to seconds.
     */
    private static final double toSeconds = 1.0e-9;
    /**
     * Default tolerance for numerical methods of solving constraints.
     */
    public static final double DEFAULT_CONSTRAINT_TOLERANCE = 1E-4;
    /**
     * The MolecularAssembly associated with this force field energy.
     */
    protected final MolecularAssembly molecularAssembly;
    /**
     * The array of Atoms being evaluated.
     */
    private Atom[] atoms;
    /**
     * The boundary conditions used when evaluating the force field energy.
     */
    private Crystal crystal;
    /**
     * The non-bonded cut-off plus buffer distance (Angstroms).
     */
    private double cutoffPlusBuffer;
    /**
     * The Parallel Java ParallelTeam instance.
     */
    private final ParallelTeam parallelTeam;
    /**
     * A Parallel Java Region used to evaluate Bonded energy values.
     */
    private BondedRegion bondedRegion;
    /**
     * An instance of the STATE enumeration to specify calculation of slowly varying energy terms, fast varying or both.
     */
    private STATE state = STATE.BOTH;
    /**
     * An array of Bond terms.
     */
    private Bond[] bonds;
    /**
     * An array of Angle terms.
     */
    private Angle[] angles;
    /**
     * An array of Stretch-Bend terms.
     */
    private StretchBend[] stretchBends;
    /**
     * An array of Urey-Bradley terms.
     */
    private UreyBradley[] ureyBradleys;
    /**
     * An array of Out of Plane Bend terms.
     */
    private OutOfPlaneBend[] outOfPlaneBends;
    /**
     * An array of Torsion terms.
     */
    private Torsion[] torsions;
    /**
     * An array of Stretch-Torsion terms.
     */
    private StretchTorsion[] stretchTorsions;
    /**
     * An array of Angle-Torsion terms.
     */
    private AngleTorsion[] angleTorsions;
    /**
     * An array of Pi-Orbital Torsion terms.
     */
    private PiOrbitalTorsion[] piOrbitalTorsions;
    /**
     * An array of Torsion-Torsion terms.
     */
    private TorsionTorsion[] torsionTorsions;
    /**
     * An array of Improper Torsion terms.
     */
    private ImproperTorsion[] improperTorsions;
    /**
     * An array of Bond Restraint terms.
     */
    private RestraintBond[] restraintBonds;
    /**
     * An array of Coordinate Restraint terms.
     */
    private final List<CoordRestraint> coordRestraints;
    /**
     * An NCS restraint term.
     */
    private final NCSRestraint ncsRestraint;
    /**
     * A coordinate restraint term.
     */
    private final CoordRestraint autoCoordRestraint;
    /**
     * A Center-of-Mass restraint term.
     */
    private final COMRestraint comRestraint;
    /**
     * Non-Bonded van der Waals energy.
     */
    private final VanDerWaals vanderWaals;
    /**
     * Particle-Mesh Ewald electrostatic energy.
     */
    private ParticleMeshEwald particleMeshEwald;
    /**
     * Number of atoms in the system.
     */
    private int nAtoms;
    /**
     * Number of bond terms in the system.
     */
    private int nBonds;
    /**
     * Number of angle terms in the system.
     */
    private int nAngles;
    /**
     * Number of stretch-bend terms in the system.
     */
    private int nStretchBends;
    /**
     * Number of Urey-Bradley terms in the system.
     */
    private int nUreyBradleys;
    /**
     * Number of Out of Plane Bend terms in the system.
     */
    private int nOutOfPlaneBends;
    /**
     * Number of Torsion terms in the system.
     */
    private int nTorsions;
    /**
     * Number of Angle-Torsion terms in the system.
     */
    private int nAngleTorsions;
    /**
     * Number of Stretch-Torsion terms in the system.
     */
    private int nStretchTorsions;
    /**
     * Number of Improper Torsion terms in the system.
     */
    private int nImproperTorsions;
    /**
     * Number of Pi-Orbital Torsion terms in the system.
     */
    private int nPiOrbitalTorsions;
    /**
     * Number of Torsion-Torsion terms in the system.
     */
    private int nTorsionTorsions;
    /**
     * Number of Restraint Bond terms in the system.
     */
    private int nRestraintBonds = 0;
    /**
     * Number of van der Waals interactions evaluated.
     */
    private int nVanDerWaalInteractions;
    /**
     * Number of electrostatic interactions evaluated.
     */
    private int nPermanentInteractions;
    /**
     * Number of implicit solvent interactions evaluated.
     */
    private int nGKInteractions;
    /**
     * Evaluate Bond energy terms.
     */
    private boolean bondTerm;
    /**
     * Original state of the Bond energy term flag.
     */
    private boolean bondTermOrig;
    /**
     * Evaluate Angle energy terms.
     */
    private boolean angleTerm;
    /**
     * Original state of the Angle energy term flag.
     */
    private boolean angleTermOrig;
    /**
     * Evaluate Stretch-Bend energy terms.
     */
    private boolean stretchBendTerm;
    /**
     * Original state of the Stretch-Bend energy term flag.
     */
    private boolean stretchBendTermOrig;
    /**
     * Evaluate Urey-Bradley energy terms.
     */
    private boolean ureyBradleyTerm;
    /**
     * Original state of the Urey-Bradley energy term flag.
     */
    private boolean ureyBradleyTermOrig;
    /**
     * Evaluate Out of Plane Bend energy terms.
     */
    private boolean outOfPlaneBendTerm;
    /**
     * Original state of the Out-of-Plane Bend energy term flag.
     */
    private boolean outOfPlaneBendTermOrig;
    /**
     * Evaluate Torsion energy terms.
     */
    private boolean torsionTerm;
    /**
     * Original state of the Torsion energy term flag.
     */
    private boolean torsionTermOrig;
    /**
     * Evaluate Stretch-Torsion energy terms.
     */
    private boolean stretchTorsionTerm;
    /**
     * Original state of the Stretch-Torsion energy term flag.
     */
    private boolean stretchTorsionTermOrig;
    /**
     * Evaluate Angle-Torsion energy terms.
     */
    private boolean angleTorsionTerm;
    /**
     * Original state of the Angle-Torsion energy term flag.
     */
    private boolean angleTorsionTermOrig;
    /**
     * Evaluate Improper Torsion energy terms.
     */
    private boolean improperTorsionTerm;
    /**
     * Original state of the Improper Torsion energy term flag.
     */
    private boolean improperTorsionTermOrig;
    /**
     * Evaluate Pi-Orbital Torsion energy terms.
     */
    private boolean piOrbitalTorsionTerm;
    /**
     * Original state of the Pi-Orbital Torsion energy term flag.
     */
    private boolean piOrbitalTorsionTermOrig;
    /**
     * Evaluate Torsion-Torsion energy terms.
     */
    private boolean torsionTorsionTerm;
    /**
     * Original state of the Torsion-Torsion energy term flag.
     */
    private boolean torsionTorsionTermOrig;
    /**
     * Evaluate Restraint Bond energy terms.
     */
    private boolean restraintBondTerm;
    /**
     * Original state of the Restraint Bond energy term flag.
     */
    private boolean restraintBondTermOrig;
    /**
     * Evaluate van der Waals energy term.
     */
    private boolean vanderWaalsTerm;
    /**
     * Original state of the van der Waals energy term flag.
     */
    private boolean vanderWaalsTermOrig;
    /**
     * Evaluate permanent multipole electrostatics energy term.
     */
    private boolean multipoleTerm;
    /**
     * Original state of the multipole energy term flag.
     */
    private boolean multipoleTermOrig;
    /**
     * Evaluate polarization energy term.
     */
    private boolean polarizationTerm;
    /**
     * Original state of the polarization energy term flag.
     */
    private boolean polarizationTermOrig;
    /**
     * Evaluate generalized Kirkwood energy term.
     */
    private boolean generalizedKirkwoodTerm;
    /**
     * Original state of the GK energy term flag.
     */
    private boolean generalizedKirkwoodTermOrig;
    /**
     * Evaluate NCS energy term.
     */
    private boolean ncsTerm;
    /**
     * Original state of the NCS energy term flag.
     */
    private boolean ncsTermOrig;
    /**
     * Evaluate Restrain energy term.
     */
    private boolean restrainTerm;
    /**
     * Original state of the Restrain energy term flag.
     */
    private boolean restrainTermOrig;
    /**
     * Evaluate COM energy term.
     */
    private boolean comTerm;
    /**
     * Original state of the COM energy term flag.
     */
    private boolean comTermOrig;
    /**
     * Flag to indicate hydrogen bonded terms should be scaled up.
     */
    private boolean rigidHydrogens;
    /**
     * Scale factor for increasing the strength of bonded terms involving hydrogen atoms.
     */
    private double rigidScale;
    /**
     * The total Bond term energy.
     */
    private double bondEnergy;
    /**
     * The RMSD of all Bond distances relative to their equilibrium values.
     */
    private double bondRMSD;
    /**
     * The total Angle term energy.
     */
    private double angleEnergy;
    /**
     * The RMSD of all Angle bends relative to their equilibrium values.
     */
    private double angleRMSD;
    /**
     * The total Stretch-Bend term energy.
     */
    private double stretchBendEnergy;
    /**
     * The total Urey-Bradley term energy.
     */
    private double ureyBradleyEnergy;
    /**
     * The total Out-of-Plane term energy.
     */
    private double outOfPlaneBendEnergy;
    /**
     * The total Torsion term energy.
     */
    private double torsionEnergy;
    /**
     * The total Stretch-Torsion term energy.
     */
    private double stretchTorsionEnergy;
    /**
     * The total Angle-Torsion term energy.
     */
    private double angleTorsionEnergy;
    /**
     * The total Improper Torsion term energy.
     */
    private double improperTorsionEnergy;
    /**
     * The total Pi-Orbital Torsion term energy.
     */
    private double piOrbitalTorsionEnergy;
    /**
     * The total Torsion-Torsion term energy.
     */
    private double torsionTorsionEnergy;
    /**
     * The total Restraint Bond term energy.
     */
    private double restraintBondEnergy;
    /**
     * The total energy for all Bonded Energy terms.
     */
    private double totalBondedEnergy;
    /**
     * The total van der Waals energy.
     */
    private double vanDerWaalsEnergy;
    /**
     * The permanent Multipole energy.
     */
    private double permanentMultipoleEnergy;
    /**
     * The permanent multipole real space energy.
     */
    private double permanentRealSpaceEnergy;
    /**
     * The total polarization energy.
     */
    private double polarizationEnergy;
    /**
     * The total electrostatic energy.
     */
    private double totalMultipoleEnergy;
    /**
     * The total non-bonded energy (van der Waals plus electrostatic).
     */
    private double totalNonBondedEnergy;
    /**
     * The solvation energy (GK).
     */
    private double solvationEnergy;
    /**
     * The total NCS Energy.
     */
    private double ncsEnergy;
    /**
     * The total Restrain Energy.
     */
    private double restrainEnergy;
    /**
     * The total COM Restraint Energy.
     */
    private double comRestraintEnergy;
    /**
     * The total system energy.
     */
    private double totalEnergy;
    /**
     * Time to evaluate Bond terms.
     */
    private long bondTime;
    /**
     * Time to evaluate Angle terms.
     */
    private long angleTime;
    /**
     * Time to evaluate Stretch-Bend terms.
     */
    private long stretchBendTime;
    /**
     * Time to evaluate Urey-Bradley terms.
     */
    private long ureyBradleyTime;
    /**
     * Time to evaluate Out-Of-Plane Bend terms.
     */
    private long outOfPlaneBendTime;
    /**
     * Time to evaluate Torsion terms.
     */
    private long torsionTime;
    /**
     * Time to evaluate Angle-Torsion terms.
     */
    private long angleTorsionTime;
    /**
     * Time to evaluate Stretch-Torsion terms.
     */
    private long stretchTorsionTime;
    /**
     * Time to evaluate Pi-Orbital Torsion terms.
     */
    private long piOrbitalTorsionTime;
    /**
     * Time to evaluate Improper Torsion terms.
     */
    private long improperTorsionTime;
    /**
     * Time to evaluate Torsion-Torsion terms.
     */
    private long torsionTorsionTime;
    /**
     * Time to evaluate van der Waals term.
     */
    private long vanDerWaalsTime;
    /**
     * Time to evaluate electrostatics term.
     */
    private long electrostaticTime;
    /**
     * Time to evaluate Restraint Bond term.
     */
    private long restraintBondTime;
    /**
     * Time to evaluate NCS term.
     */
    private long ncsTime;
    /**
     * Time to evaluate coordinate restraint term.
     */
    private long coordRestraintTime;
    /**
     * Time to evaluate Center of Mass restraint term.
     */
    private long comRestraintTime;
    /**
     * Time to evaluate all energy terms.
     */
    private long totalTime;
    /**
     * Current value of the Lambda state variable.
     */
    protected double lambda = 1.0;
    /**
     * Indicates use of the Lambda state variable.
     */
    protected boolean lambdaTerm;
    /**
     * Indicates only bonded energy terms effected by Lambda should be evaluated.
     */
    boolean lambdaBondedTerms = false;
    /**
     * Indicates application of lambda scaling to all Torsion based energy terms.
     */
    private boolean lambdaTorsions;
    /**
     * Optimization scaling value to use for each degree of freedom.
     */
    protected double[] optimizationScaling = null;
    /**
     * Value of each degree of freedom.
     */
    private double[] xyz;
    /**
     * Enable verbose printing if large energy gradient components are observed.
     */
    private boolean printOnFailure;
    /**
     * If the absolute value of a gradient component is greater than "maxDebugGradient", verbose logging results.
     */
    final double maxDebugGradient;
    /**
     * Flag to indicate proper shutdown of the ForceFieldEnergy.
     */
    boolean destroyed = false;

    /**
     * Indicate resolution of this ForceFieldEnergy (TODO: needs further testing).
     */
    private Resolution resolution = Resolution.AMOEBA;

    /**
     * Constant pH extended system (TODO: needs further testing).
     */
    private ExtendedSystem esvSystem = null;
    private boolean esvTerm;
    private double esvBias;

    /**
     * Relative solvation term (TODO: needs further testing).
     */
    private RelativeSolvation relativeSolvation;
    private int nRelativeSolvations;
    private boolean relativeSolvationTerm;
    private double relativeSolvationEnergy;
    private Platform platform = Platform.FFX;

    private final List<Constraint> constraints;


    /**
     * <p>
     * Constructor for ForceFieldEnergy.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     */
    protected ForceFieldEnergy(MolecularAssembly molecularAssembly) {
        this(molecularAssembly, null);
    }

    /**
     * <p>
     * Constructor for ForceFieldEnergy.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param restraints        list of {@link ffx.potential.nonbonded.CoordRestraint} objects.
     */
    protected ForceFieldEnergy(MolecularAssembly molecularAssembly, List<CoordRestraint> restraints) {
        this(molecularAssembly, restraints, ParallelTeam.getDefaultThreadCount());
    }

    /**
     * <p>Constructor for ForceFieldEnergy.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
     * @param restraints        a {@link java.util.List} object.
     * @param numThreads        a int.
     */
    protected ForceFieldEnergy(MolecularAssembly molecularAssembly, List<CoordRestraint> restraints, int numThreads) {
        // Get a reference to the sorted atom array.
        this.molecularAssembly = molecularAssembly;
        atoms = molecularAssembly.getAtomArray();
        nAtoms = atoms.length;
        xyz = new double[nAtoms * 3];

        // Check that atom ordering is correct and count the number of active atoms.
        for (int i = 0; i < nAtoms; i++) {
            int index = atoms[i].getXyzIndex() - 1;
            assert (i == index);
        }

        // Enforce that the number of threads be less than or equal to the number of atoms.
        /*int nThreads = ParallelTeam.getDefaultThreadCount();
        nThreads = nAtoms < nThreads ? nAtoms : nThreads;*/
        int nThreads = nAtoms < numThreads ? nAtoms : numThreads;
        parallelTeam = new ParallelTeam(nThreads);

        ForceField forceField = molecularAssembly.getForceField();
        String name = forceField.toString().toUpperCase();

        logger.info(format(" Constructing Force Field %s", name));
        logger.info(format("\n SMP threads:                        %10d", nThreads));

        bondTerm = forceField.getBoolean(ForceFieldBoolean.BONDTERM, true);
        angleTerm = forceField.getBoolean(ForceFieldBoolean.ANGLETERM, true);
        stretchBendTerm = forceField.getBoolean(ForceFieldBoolean.STRBNDTERM, true);
        ureyBradleyTerm = forceField.getBoolean(ForceFieldBoolean.UREYTERM, true);
        outOfPlaneBendTerm = forceField.getBoolean(ForceFieldBoolean.OPBENDTERM, true);
        torsionTerm = forceField.getBoolean(ForceFieldBoolean.TORSIONTERM, true);
        stretchTorsionTerm = forceField.getBoolean(ForceFieldBoolean.STRTORSTERM, true);
        angleTorsionTerm = forceField.getBoolean(ForceFieldBoolean.ANGTORSTERM, true);
        piOrbitalTorsionTerm = forceField.getBoolean(ForceFieldBoolean.PITORSTERM, true);
        torsionTorsionTerm = forceField.getBoolean(ForceFieldBoolean.TORTORTERM, true);
        improperTorsionTerm = forceField.getBoolean(ForceFieldBoolean.IMPROPERTERM, true);
        vanderWaalsTerm = forceField.getBoolean(ForceFieldBoolean.VDWTERM, true);
        if (vanderWaalsTerm) {
            multipoleTerm = forceField.getBoolean(ForceFieldBoolean.MPOLETERM, true);
            if (multipoleTerm) {
                polarizationTerm = forceField.getBoolean(ForceFieldBoolean.POLARIZETERM, true);
                generalizedKirkwoodTerm = forceField.getBoolean(ForceFieldBoolean.GKTERM, false);
            } else {
                // If multipole electrostatics is turned off, turn off all electrostatics.
                polarizationTerm = false;
                generalizedKirkwoodTerm = false;
            }
        } else {
            // If van der Waals is turned off, turn off all non-bonded terms.
            multipoleTerm = false;
            polarizationTerm = false;
            generalizedKirkwoodTerm = false;
        }
        restraintBondTerm = false;
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);
        restrainTerm = forceField.getBoolean(ForceFieldBoolean.RESTRAINTERM, false);
        comTerm = forceField.getBoolean(ForceFieldBoolean.COMRESTRAINTERM, false);
        lambdaTorsions = forceField.getBoolean(ForceFieldBoolean.TORSION_LAMBDATERM, false);
        printOnFailure = forceField.getBoolean(ForceFieldBoolean.PRINT_ON_FAILURE, false);

        // For RESPA
        bondTermOrig = bondTerm;
        angleTermOrig = angleTerm;
        stretchBendTermOrig = stretchBendTerm;
        ureyBradleyTermOrig = ureyBradleyTerm;
        outOfPlaneBendTermOrig = outOfPlaneBendTerm;
        torsionTermOrig = torsionTerm;
        angleTorsionTermOrig = angleTorsionTerm;
        stretchTorsionTermOrig = stretchTorsionTerm;
        improperTorsionTermOrig = improperTorsionTerm;
        piOrbitalTorsionTermOrig = piOrbitalTorsionTerm;
        torsionTorsionTermOrig = torsionTorsionTerm;
        restraintBondTermOrig = restraintBondTerm;
        vanderWaalsTermOrig = vanderWaalsTerm;
        multipoleTermOrig = multipoleTerm;
        polarizationTermOrig = polarizationTerm;
        generalizedKirkwoodTermOrig = generalizedKirkwoodTerm;
        ncsTermOrig = ncsTerm;
        restrainTermOrig = restrainTerm;
        comTermOrig = comTerm;

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

            if (a == 1.0 && b == 1.0 && c == 1.0) {
                String message = " A-, B-, and C-axis values equal to 1.0.";
                logger.info(message);
                throw new Exception(message);
            }

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

            // Turn off reciprocal space calculations.
            forceField.addForceFieldDouble(ForceFieldDouble.EWALD_ALPHA, 0.0);

            // Specify some dummy values for the crystal.
            spacegroup = "P1";
            a = 2.0 * maxr;
            b = 2.0 * maxr;
            c = 2.0 * maxr;
            alpha = 90.0;
            beta = 90.0;
            gamma = 90.0;
        }
        Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
        unitCell.setAperiodic(aperiodic);

        // Define the cutoff lengths.
        double buff = 2.0;
        double defaultVdwCut = 12.0;
        if (unitCell.aperiodic()) {
            double maxDim = max(max(unitCell.a, unitCell.b), unitCell.c);
            defaultVdwCut = (maxDim * 0.5) - (buff + 1.0);
        }
        double vdwOff = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, defaultVdwCut);
        double ewaldOff = aperiodic ? ParticleMeshEwald.APERIODIC_DEFAULT_EWALD_CUTOFF : ParticleMeshEwald.PERIODIC_DEFAULT_EWALD_CUTOFF;
        ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, ewaldOff);
        if (ewaldOff > vdwOff) {
            vdwOff = ewaldOff;
        }

        /*
          Neighbor list cutoff is at least max(PME cutoff, vdW cutoff).
          For a non-frozen neighbor list, it is max(PME, vdW, GK).
          Then, if a larger neighbor-list cutoff is specified, we use that.
         
          GK cutoff may be > neighbor-list cutoff for frozen neighbor lists.
          This indicates we are running full real-space GK on a frozen neighbor list.
         
          Error is indicated for a frozen neighbor list if a specified NL cutoff < PME,vdW.
          Error is indicated for a non-frozen neighbor list if the NL cutoff was specified by the user.
         */

        double nlistCutoff = vdwOff;
        // Check for a frozen neighbor list.
        boolean disabledNeighborUpdates = forceField.getBoolean(ForceField.ForceFieldBoolean.DISABLE_NEIGHBOR_UPDATES, false);
        if (disabledNeighborUpdates) {
            nlistCutoff = forceField.getDouble(ForceFieldDouble.NEIGHBOR_LIST_CUTOFF, vdwOff);
            logger.info(format(" Neighbor list updates disabled; interactions will " +
                    "only be calculated between atoms that started the simulation " +
                    "within a radius of %9.3g Angstroms of each other", nlistCutoff));
            if (nlistCutoff < vdwOff) {
                logger.severe(format(" Specified a neighbor-list cutoff %10.4g < max(PME, vdW) cutoff %10.4g !", nlistCutoff, vdwOff));
            }
        } else {
            nlistCutoff = Math.max(forceField.getDouble(ForceFieldDouble.GK_CUTOFF, 0), nlistCutoff);
            if (forceField.hasDouble(ForceFieldDouble.NEIGHBOR_LIST_CUTOFF)) {
                logger.severe(" Specified a neighbor list cutoff for a non-frozen neighbor list!");
            }
        }

        cutoffPlusBuffer = nlistCutoff + buff;
        unitCell = configureNCS(forceField, unitCell);

        // If necessary, create a ReplicatesCrystal.
        if (!aperiodic) {
            double cutOff2 = 2.0 * cutoffPlusBuffer;
            this.crystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, cutOff2);
            logger.info(format("\n Density:                                %6.3f (g/cc)", crystal.getDensity(molecularAssembly.getMass())));
            logger.info(crystal.toString());
        } else {
            this.crystal = unitCell;
        }

        if (!unitCell.aperiodic() && unitCell.spaceGroup.number == 1) {
            ncsTerm = forceField.getBoolean(ForceFieldBoolean.NCSTERM, false);
            ncsTermOrig = ncsTerm;
        } else {
            ncsTerm = false;
            ncsTermOrig = false;
        }

        rigidHydrogens = forceField.getBoolean(ForceFieldBoolean.RIGID_HYDROGENS, false);
        rigidScale = forceField.getDouble(ForceFieldDouble.RIGID_SCALE, 10.0);

        nRelativeSolvations = 0;
        String relSolvLibrary = forceField.getString(ForceFieldString.RELATIVE_SOLVATION, "NONE").toUpperCase();
        SolvationLibrary library = SolvationLibrary.valueOf(relSolvLibrary);
        if (library == SolvationLibrary.AUTO) {
            if (generalizedKirkwoodTerm && name.toUpperCase().contains("OPLS")) {
                library = SolvationLibrary.MACCALLUM_TIP4P;
                // Change when we have good Generalized Born numbers of our own for OPLS.
            } else if (generalizedKirkwoodTerm) {
                library = SolvationLibrary.GK;
            } else {
                library = SolvationLibrary.NONE;
            }
        }
        if (library != SolvationLibrary.NONE) {
            relativeSolvationTerm = true;
            relativeSolvation = new RelativeSolvation(library, forceField);
        } else {
            relativeSolvationTerm = false;
            relativeSolvation = null;
        }

        boolean checkAllNodeCharges = forceField.getBoolean(ForceFieldBoolean.CHECK_ALL_NODE_CHARGES, false);

        if (rigidScale <= 1.0) {
            rigidScale = 1.0;
        }

        logger.info("\n Bonded Terms");
        if (rigidHydrogens && rigidScale > 1.0) {
            logger.info(format("  Rigid hydrogens:                   %10.2f", rigidScale));
        }

        // Collect, count, pack and sort bonds.
        if (bondTerm) {
            ArrayList<Bond> bond = molecularAssembly.getBondList();
            nBonds = bond.size();
            bonds = bond.toArray(new Bond[0]);
            sort(bonds);
            if (nBonds > 0) {
                logger.info(format("  Bonds:                             %10d", nBonds));
            }
        } else {
            nBonds = 0;
            bonds = null;
        }

        // Collect, count, pack and sort angles.
        if (angleTerm) {
            ArrayList<Angle> angle = molecularAssembly.getAngleList();
            nAngles = angle.size();
            angles = angle.toArray(new Angle[0]);
            sort(angles);
            if (nAngles > 0) {
                logger.info(format("  Angles:                            %10d", nAngles));
            }
        } else {
            nAngles = 0;
            angles = null;
        }

        // Collect, count, pack and sort stretch-bends.
        if (stretchBendTerm) {
            ArrayList<StretchBend> stretchBend = molecularAssembly.getStretchBendList();
            nStretchBends = stretchBend.size();
            stretchBends = stretchBend.toArray(new StretchBend[0]);
            sort(stretchBends);
            if (nStretchBends > 0) {
                logger.info(format("  Stretch-Bends:                     %10d", nStretchBends));
            }
        } else {
            nStretchBends = 0;
            stretchBends = null;
        }

        // Collect, count, pack and sort Urey-Bradleys.
        if (ureyBradleyTerm) {
            ArrayList<UreyBradley> ureyBradley = molecularAssembly.getUreyBradleyList();
            nUreyBradleys = ureyBradley.size();
            ureyBradleys = ureyBradley.toArray(new UreyBradley[0]);
            sort(ureyBradleys);
            if (nUreyBradleys > 0) {
                logger.info(format("  Urey-Bradleys:                     %10d", nUreyBradleys));
            }
        } else {
            nUreyBradleys = 0;
            ureyBradleys = null;
        }

        // Set a multiplier on the force constants of bonded terms containing hydrogens.
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
            ArrayList<OutOfPlaneBend> outOfPlaneBend = molecularAssembly.getOutOfPlaneBendList();
            nOutOfPlaneBends = outOfPlaneBend.size();
            outOfPlaneBends = outOfPlaneBend.toArray(new OutOfPlaneBend[0]);
            sort(outOfPlaneBends);
            if (nOutOfPlaneBends > 0) {
                logger.info(format("  Out-of-Plane Bends:                %10d", nOutOfPlaneBends));
            }
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        double torsionScale = forceField.getDouble(ForceFieldDouble.TORSION_SCALE, 1.0);
        if (torsionScale != 1.0) {
            forceField.setTorsionScale(torsionScale);
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<Torsion> torsion = molecularAssembly.getTorsionList();
            nTorsions = torsion.size();
            torsions = torsion.toArray(new Torsion[0]);
            if (nTorsions > 0) {
                if (torsionScale == 1.0) {
                    logger.info(format("  Torsions:                          %10d", nTorsions));
                } else {
                    logger.info(format("  Torsions (%5.2f):                  %10d", torsionScale, nTorsions));
                }
            }
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort stretch torsions.
        if (stretchTorsionTerm) {
            ArrayList<StretchTorsion> stretchTorsion = molecularAssembly.getStretchTorsionList();
            nStretchTorsions = stretchTorsion.size();
            stretchTorsions = stretchTorsion.toArray(new StretchTorsion[0]);
            if (nStretchTorsions > 0) {
                if (torsionScale == 1.0) {
                    logger.info(format("  Stretch-Torsions:                  %10d", nStretchTorsions));
                } else {
                    logger.info(format("  Stretch-Torsions (%5.2f):          %10d", torsionScale, nStretchTorsions));
                }
            }
        } else {
            nStretchTorsions = 0;
            stretchTorsions = null;
        }

        // Collect, count, pack and sort angle torsions.
        if (angleTorsionTerm) {
            ArrayList<AngleTorsion> angleTorsion = molecularAssembly.getAngleTorsionList();
            nAngleTorsions = angleTorsion.size();
            angleTorsions = angleTorsion.toArray(new AngleTorsion[0]);
            if (nAngleTorsions > 0) {
                if (torsionScale == 1.0) {
                    logger.info(format("  Angle-Torsions:                    %10d", nAngleTorsions));
                } else {
                    logger.info(format("  Angle-Torsions (%5.2f):            %10d", torsionScale, nAngleTorsions));
                }
            }
        } else {
            nAngleTorsions = 0;
            angleTorsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<PiOrbitalTorsion> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = piOrbitalTorsion.size();
            piOrbitalTorsions = piOrbitalTorsion.toArray(new PiOrbitalTorsion[0]);
            if (nPiOrbitalTorsions > 0) {
                if (torsionScale == 1.0) {
                    logger.info(format("  Pi-Orbital Torsions:               %10d", nPiOrbitalTorsions));
                } else {
                    logger.info(format("  Pi-Orbital Torsions (%5.2f):       %10d", torsionScale, nPiOrbitalTorsions));
                }
            }
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<TorsionTorsion> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = torsionTorsion.size();
            torsionTorsions = torsionTorsion.toArray(new TorsionTorsion[0]);
            if (nTorsionTorsions > 0) {
                logger.info(format("  Torsion-Torsions:                  %10d", nTorsionTorsions));
            }
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }

        // Collect, count, pack and sort improper torsions.
        if (improperTorsionTerm) {
            ArrayList<ImproperTorsion> improperTorsion = molecularAssembly.getImproperTorsionList();
            nImproperTorsions = improperTorsion.size();
            improperTorsions = improperTorsion.toArray(new ImproperTorsion[0]);
            if (nImproperTorsions > 0) {
                logger.info(format("  Improper Torsions:                 %10d", nImproperTorsions));
            }
        } else {
            nImproperTorsions = 0;
            improperTorsions = null;
        }

        logger.info("\n Non-Bonded Terms");

        int[] molecule = molecularAssembly.getMoleculeNumbers();
        if (vanderWaalsTerm) {
            vanderWaals = new VanDerWaals(atoms, molecule, crystal, forceField, parallelTeam, vdwOff, nlistCutoff);
        } else {
            vanderWaals = null;
        }

        if (multipoleTerm) {
            ELEC_FORM form;
            if (name.contains("OPLS") || name.contains("AMBER") || name.contains("CHARMM")) {
                form = ELEC_FORM.FIXED_CHARGE;
            } else {
                form = ELEC_FORM.PAM;
            }

            boolean pmeQI = forceField.getBoolean(ForceFieldBoolean.PME_QI, false);

            if (pmeQI) {
                particleMeshEwald = new ParticleMeshEwaldQI(atoms, molecule, forceField, crystal,
                        vanderWaals.getNeighborList(), form, parallelTeam);
            } else {
                particleMeshEwald = new ParticleMeshEwaldCart(atoms, molecule, forceField, crystal,
                        vanderWaals.getNeighborList(), form, parallelTeam);
            }
            double charge = molecularAssembly.getCharge(checkAllNodeCharges);
            logger.info(format("\n  Overall system charge:             %10.3f", charge));
        } else {
            particleMeshEwald = null;
        }

        if (ncsTerm) {
            String sg = forceField.getString(ForceFieldString.NCSGROUP, "P 1");
            Crystal ncsCrystal = new Crystal(a, b, c, alpha, beta, gamma, sg);
            ncsRestraint = new NCSRestraint(atoms, forceField, ncsCrystal);
        } else {
            ncsRestraint = null;
        }

        coordRestraints = new ArrayList<>();
        if (restraints != null) {
            coordRestraints.addAll(restraints);
        }
        if (restrainTerm) {
            this.autoCoordRestraint = new CoordRestraint(atoms, forceField);
            coordRestraints.add(autoCoordRestraint);
        } else {
            autoCoordRestraint = null;
        }
        if (!coordRestraints.isEmpty()) {
            restrainTerm = true;
            restrainTermOrig = restrainTerm;
            logger.log(Level.FINE, " restrainTerm set true");
        }

        if (comTerm) {
            Polymer[] polymers = molecularAssembly.getChains();
            List<MSNode> molecules = molecularAssembly.getMolecules();
            List<MSNode> waters = molecularAssembly.getWaters();
            List<MSNode> ions = molecularAssembly.getIons();
            comRestraint = new COMRestraint(atoms, polymers, molecules, waters, ions, forceField);
        } else {
            comRestraint = null;
        }

        bondedRegion = new BondedRegion();

        maxDebugGradient = forceField.getDouble(ForceFieldDouble.MAX_DEBUG_GRADIENT, Double.POSITIVE_INFINITY);

        molecularAssembly.setPotential(this);

        // Add restrain-bond records. If no restrain-distance records exist, the empty array will be returned.
        String[] bondRestraints = molecularAssembly.getProperties().getStringArray("restrain-distance");
        for (String bondRest : bondRestraints) {
            try {
                String[] toks = bondRest.split("\\s+");
                if (toks.length < 2) {
                    throw new IllegalArgumentException(format(" restrain-distance value %s could not be parsed!", bondRest));
                }
                int at1 = Integer.parseInt(toks[0]) - 1;
                int at2 = Integer.parseInt(toks[1]) - 1;

                double forceConst = 100.0;
                double flatBottomRadius = 0;
                Atom a1 = atoms[at1];
                Atom a2 = atoms[at2];

                if (toks.length > 2) {
                    forceConst = Double.parseDouble(toks[2]);
                }
                double dist;
                switch (toks.length) {
                    case 3:
                        double[] xyz1 = new double[3];
                        xyz1 = a1.getXYZ(xyz1);
                        double[] xyz2 = new double[3];
                        xyz2 = a2.getXYZ(xyz2);
                        dist = crystal.minDistOverSymOps(xyz1, xyz2);
                        break;
                    case 4:
                        dist = Double.parseDouble(toks[3]);
                        break;
                    case 5:
                        double minDist = Double.parseDouble(toks[3]);
                        double maxDist = Double.parseDouble(toks[4]);
                        dist = 0.5 * (minDist + maxDist);
                        flatBottomRadius = 0.5 * Math.abs(maxDist - minDist);
                        break;
                    default:
                        throw new IllegalArgumentException(format(" restrain-distance value %s could not be parsed!", bondRest));
                }

                setRestraintBond(a1, a2, dist, forceConst, flatBottomRadius);
            } catch (Exception ex) {
                logger.info(format(" Exception in parsing restrain-distance: %s", ex.toString()));
            }
        }

        String constraintStrings = forceField.getString(ForceFieldString.CONSTRAIN, forceField.getString(ForceFieldString.RATTLE, null));
        if (constraintStrings != null) {
            constraints = new ArrayList<>();

            logger.info(format(" Experimental: parsing constraints option %s", constraintStrings));
            if (constraintStrings.isEmpty() || constraintStrings.matches("^\\s*$")) {
                // Assume constraining only X-H bonds (i.e. RIGID-HYDROGEN).
                logger.info(" Constraining X-H bonds.");
                logger.severe(" TODO: Implement this.");
            } else {
                String[] constraintToks = constraintStrings.split("\\s+");
                for (String tok : constraintToks) {
                    if (tok.equalsIgnoreCase("WATER")) {
                        logger.info(" Constraining waters to be rigid based on angle & bonds.");
                        // XYZ files, in particular, have waters mislabeled as generic Molecules.
                        // First, find any such mislabeled water.
                        Stream<MSNode> settleStream = molecularAssembly.getMolecules().stream().
                                filter((MSNode m) -> m.getAtomList().size() == 3).
                                filter((MSNode m) -> {
                                    List<Atom> atoms = m.getAtomList();
                                    Atom O = null;
                                    List<Atom> H = new ArrayList<>(2);
                                    for (Atom at : atoms) {
                                        int atN = at.getAtomicNumber();
                                        if (atN == 8) {
                                            O = at;
                                        } else if (atN == 1) {
                                            H.add(at);
                                        }
                                    }
                                    return O != null && H.size() == 2;
                                });
                        // Now concatenate the stream with the properly labeled waters.
                        settleStream = Stream.concat(settleStream, molecularAssembly.getWaters().stream());
                        // Map them into new Settle constraints and collect.
                        List<SettleConstraint> settleConstraints = settleStream.map((MSNode m) -> m.getAngleList().get(0)).
                                map(SettleConstraint::new).
                                collect(Collectors.toList());
                        constraints.addAll(settleConstraints);

                    } else {
                        logger.severe(" Implement constraints that aren't SETTLE constraints.");
                    }
                }
            }

            logger.info(format(" Added %d constraints.", constraints.size()));
        } else {
            constraints = Collections.emptyList();
        }

        if (lambdaTerm) {
            this.setLambda(1.0);
        }
    }

    /**
     * Static factory method to create a ForceFieldEnergy, possibly via FFX or
     * OpenMM implementations.
     *
     * @param assembly To create FFE over
     * @return a {@link ffx.potential.ForceFieldEnergy} object.
     */
    public static ForceFieldEnergy energyFactory(MolecularAssembly assembly) {
        return energyFactory(assembly, null);
    }

    /**
     * Static factory method to create a ForceFieldEnergy, possibly via FFX or
     * OpenMM implementations.
     *
     * @param assembly   To create FFE over
     * @param restraints Harmonic restraints
     * @return a {@link ffx.potential.ForceFieldEnergy} object.
     */
    public static ForceFieldEnergy energyFactory(MolecularAssembly assembly, List<CoordRestraint> restraints) {
        return energyFactory(assembly, restraints, ParallelTeam.getDefaultThreadCount());
    }

    /**
     * Static factory method to create a ForceFieldEnergy, possibly via FFX or
     * OpenMM implementations.
     *
     * @param assembly   To create FFE over
     * @param restraints Harmonic restraints
     * @param numThreads Number of threads to use for FFX energy
     * @return A ForceFieldEnergy on some Platform
     */
    public static ForceFieldEnergy energyFactory(MolecularAssembly assembly, List<CoordRestraint> restraints, int numThreads) {
        ForceField ffield = assembly.getForceField();
        String platformString = toEnumForm(ffield.getString(ForceFieldString.PLATFORM, "FFX"));
        try {
            Platform platform = Platform.valueOf(platformString);
            switch (platform) {
                case OMM:
                case OMM_REF: // Should be split from the code once we figure out how to specify a kernel.
                case OMM_CUDA:
                    try {
                        return new ForceFieldEnergyOpenMM(assembly, platform, restraints, numThreads);
                    } catch (Exception ex) {
                        logger.warning(format(" Exception creating ForceFieldEnergyOpenMM: %s", ex));
                        ex.printStackTrace();

                        ForceFieldEnergy ffxEnergy = assembly.getPotentialEnergy();
                        if (ffxEnergy == null) {
                            ffxEnergy = new ForceFieldEnergy(assembly, restraints, numThreads);
                            assembly.setPotential(ffxEnergy);
                        }
                        return ffxEnergy;
                    }
                case OMM_OPENCL:
                case OMM_OPTCPU:
                    logger.warning(format(" Platform %s not supported; defaulting to FFX", platform));
                case FFX:
                default:
                    ForceFieldEnergy ffxEnergy = assembly.getPotentialEnergy();
                    if (ffxEnergy == null) {
                        ffxEnergy = new ForceFieldEnergy(assembly, restraints, numThreads);
                        assembly.setPotential(ffxEnergy);
                    }
                    return ffxEnergy;
            }
        } catch (IllegalArgumentException | NullPointerException ex) {
            logger.warning(format(" String %s did not match a known energy implementation", platformString));
            ForceFieldEnergy ffxEnergy = assembly.getPotentialEnergy();
            if (ffxEnergy == null) {
                ffxEnergy = new ForceFieldEnergy(assembly, restraints, numThreads);
                assembly.setPotential(ffxEnergy);
            }
            return ffxEnergy;
        }
    }

    /**
     * Overwrites current esvSystem if present. Multiple ExtendedSystems is
     * possible but unnecessary; add all ESVs to one system (per FFE, at least).
     *
     * @param system a {@link ffx.potential.extended.ExtendedSystem} object.
     */
    public void attachExtendedSystem(ExtendedSystem system) {
        if (system == null) {
            throw new IllegalArgumentException();
        }
        esvTerm = true;
        esvSystem = system;
        if (vanderWaalsTerm) {
            if (vanderWaals == null) {
                logger.warning("Null VdW during ESV setup.");
            } else {
                vanderWaals.attachExtendedSystem(system);
            }
        }
        if (multipoleTerm) {
            if (particleMeshEwald == null) {
                logger.warning("Null PME during ESV setup.");
            }
            if (!(particleMeshEwald instanceof ParticleMeshEwaldQI)) {
                logger.severe("Extended systems can attach only to Quasi-Internal PME. Try -Dpme-qi=true.");
            }
            ((ParticleMeshEwaldQI) particleMeshEwald).attachExtendedSystem(system);
        }
        if (crystal != null) {
            crystal.setSpecialPositionCutoff(0.0);
        }
        reInit();
    }

    /**
     * <p>detachExtendedSystem.</p>
     */
    public void detachExtendedSystem() {
        esvTerm = false;
        esvSystem = null;
        if (vanderWaalsTerm && vanderWaals != null) {
            vanderWaals.detachExtendedSystem();
        }
        if (multipoleTerm && particleMeshEwald != null) {
            if (particleMeshEwald instanceof ParticleMeshEwaldQI) {
                ((ParticleMeshEwaldQI) particleMeshEwald).detachExtendedSystem();
            }
        }
        reInit();
    }

    /**
     * <p>Setter for the field <code>resolution</code>.</p>
     *
     * @param resolution a {@link ffx.potential.bonded.Atom.Resolution} object.
     */
    public void setResolution(Resolution resolution) {
        this.resolution = resolution;

        if (vanderWaals != null) {
            vanderWaals.setResolution(resolution);
        }

        if (resolution == Resolution.FIXEDCHARGE) {
            multipoleTerm = false;
            polarizationTerm = false;
            generalizedKirkwoodTerm = false;
        }

    }

    private boolean keep(BondedTerm term) {
        switch (resolution) {
            case AMOEBA:
                return term.isResolution(Resolution.AMOEBA);
            case FIXEDCHARGE:
                return term.containsResolution(Resolution.FIXEDCHARGE);
            default:
                return true;
        }
    }

    /**
     * Need to remove degrees of freedom that are lost to prevent heating.
     */
    public void reInit() {
        int[] molecule;
        if (esvTerm) {
            atoms = esvSystem.getExtendedAndBackgroundAtoms();
            molecule = esvSystem.getExtendedAndBackgroundMolecule();
        } else {
            atoms = molecularAssembly.getAtomArray();
            molecule = molecularAssembly.getMoleculeNumbers();
        }
        nAtoms = atoms.length;

        xyz = new double[nAtoms * 3];
        getCoordinates(xyz);

        // Check that atom ordering is correct and count number of Active atoms.
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            int index = atom.getXyzIndex() - 1;
            if (index != i) {
                atom.setXyzIndex(i + 1);
            }
        }

        // Collect, count, pack and sort bonds.
        if (bondTerm) {
            ArrayList<Bond> bond = molecularAssembly.getBondList();
            nBonds = 0;
            for (Bond r : bond) {
                if (keep(r)) {
                    nBonds++;
                }
            }
            if (nBonds > bonds.length) {
                bonds = new Bond[nBonds];
            }
            fill(bonds, null);
            nBonds = 0;
            for (Bond r : bond) {
                if (keep(r)) {
                    bonds[nBonds++] = r;
                }
            }
            sort(bonds, 0, nBonds);
            if (nBonds > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Bonds:                             %10d", nBonds));
            }
        } else {
            nBonds = 0;
            bonds = null;
        }

        // Collect, count, pack and sort angles.
        if (angleTerm) {
            ArrayList<Angle> angle = molecularAssembly.getAngleList();
            nAngles = 0;
            for (Angle r : angle) {
                if (keep(r)) {
                    nAngles++;
                }
            }
            if (nAngles > angles.length) {
                angles = new Angle[nAngles];
            }
            fill(angles, null);
            nAngles = 0;
            for (Angle r : angle) {
                if (keep(r)) {
                    angles[nAngles++] = r;
                }
            }

            sort(angles, 0, nAngles);
            if (nAngles > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Angles:                            %10d", nAngles));
            }
        } else {
            nAngles = 0;
            angles = null;
        }

        // Collect, count, pack and sort stretch-bends.
        if (stretchBendTerm) {
            ArrayList<StretchBend> stretchBend = molecularAssembly.getStretchBendList();
            nStretchBends = 0;
            for (StretchBend r : stretchBend) {
                if (keep(r)) {
                    nStretchBends++;
                }
            }
            if (nStretchBends > stretchBends.length) {
                stretchBends = new StretchBend[nStretchBends];
            }
            fill(stretchBends, null);
            nStretchBends = 0;
            for (StretchBend r : stretchBend) {
                if (keep(r)) {
                    stretchBends[nStretchBends++] = r;
                }
            }
            sort(stretchBends, 0, nStretchBends);
            if (nStretchBends > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Stretch-Bends:                     %10d", nStretchBends));
            }
        } else {
            nStretchBends = 0;
            stretchBends = null;
        }

        // Collect, count, pack and sort Urey-Bradleys.
        if (ureyBradleyTerm) {
            ArrayList<UreyBradley> ureyBradley = molecularAssembly.getUreyBradleyList();
            nUreyBradleys = 0;
            for (UreyBradley r : ureyBradley) {
                if (keep(r)) {
                    nUreyBradleys++;
                }
            }
            if (nUreyBradleys > ureyBradleys.length) {
                ureyBradleys = new UreyBradley[nUreyBradleys];
            }
            fill(ureyBradleys, null);
            nUreyBradleys = 0;
            for (UreyBradley r : ureyBradley) {
                if (keep(r)) {
                    ureyBradleys[nUreyBradleys++] = (UreyBradley) r;
                }
            }
            sort(ureyBradleys, 0, nUreyBradleys);
            if (nUreyBradleys > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Urey-Bradleys:                     %10d", nUreyBradleys));
            }
        } else {
            nUreyBradleys = 0;
            ureyBradleys = null;
        }

        // Set a multiplier on the force constants of bonded terms containing hydrogens.
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
            ArrayList<OutOfPlaneBend> outOfPlaneBend = molecularAssembly.getOutOfPlaneBendList();
            nOutOfPlaneBends = 0;
            for (OutOfPlaneBend r : outOfPlaneBend) {
                if (keep(r)) {
                    nOutOfPlaneBends++;
                }
            }
            if (nOutOfPlaneBends > outOfPlaneBends.length) {
                outOfPlaneBends = new OutOfPlaneBend[nOutOfPlaneBends];
            }
            fill(outOfPlaneBends, null);
            nOutOfPlaneBends = 0;
            for (OutOfPlaneBend r : outOfPlaneBend) {
                if (keep(r)) {
                    outOfPlaneBends[nOutOfPlaneBends++] = r;
                }
            }
            sort(outOfPlaneBends, 0, nOutOfPlaneBends);
            if (nOutOfPlaneBends > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Out-of-Plane Bends:                %10d", nOutOfPlaneBends));
            }
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<Torsion> torsion = molecularAssembly.getTorsionList();
            nTorsions = 0;
            for (Torsion r : torsion) {
                if (keep(r)) {
                    nTorsions++;
                }
            }
            if (nTorsions >= torsions.length) {
                torsions = new Torsion[nTorsions];
            }
            fill(torsions, null);
            nTorsions = 0;
            for (Torsion r : torsion) {
                if (keep(r)) {
                    torsions[nTorsions++] = r;
                }
            }
            // Arrays.sort(torsions);
            if (nTorsions > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Torsions:                          %10d", nTorsions));
            }
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<PiOrbitalTorsion> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = 0;
            for (PiOrbitalTorsion r : piOrbitalTorsion) {
                if (keep(r)) {
                    nPiOrbitalTorsions++;
                }
            }
            if (nPiOrbitalTorsions >= piOrbitalTorsions.length) {
                piOrbitalTorsions = new PiOrbitalTorsion[nPiOrbitalTorsions];
            }
            fill(piOrbitalTorsions, null);
            nPiOrbitalTorsions = 0;
            for (PiOrbitalTorsion r : piOrbitalTorsion) {
                if (keep(r)) {
                    piOrbitalTorsions[nPiOrbitalTorsions++] = r;
                }
            }
            if (nPiOrbitalTorsions > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Pi-Orbital Torsions:               %10d", nPiOrbitalTorsions));
            }
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<TorsionTorsion> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = 0;
            for (TorsionTorsion r : torsionTorsion) {
                if (keep(r)) {
                    nTorsionTorsions++;
                }
            }
            if (nTorsionTorsions >= torsionTorsions.length) {
                torsionTorsions = new TorsionTorsion[nTorsionTorsions];
            }
            fill(torsionTorsions, null);
            nTorsionTorsions = 0;
            for (TorsionTorsion r : torsionTorsion) {
                if (keep(r)) {
                    torsionTorsions[nTorsionTorsions++] = r;
                }
            }
            if (nTorsionTorsions > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Torsion-Torsions:                  %10d", nTorsionTorsions));
            }
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }

        // Collect, count, pack and sort improper torsions.
        if (improperTorsionTerm) {
            ArrayList<ImproperTorsion> improperTorsion = molecularAssembly.getImproperTorsionList();
            nImproperTorsions = 0;
            for (ImproperTorsion r : improperTorsion) {
                if (keep(r)) {
                    nImproperTorsions++;
                }
            }
            if (nImproperTorsions >= improperTorsions.length) {
                improperTorsions = new ImproperTorsion[nImproperTorsions];
            }
            fill(improperTorsions, null);
            nImproperTorsions = 0;
            for (ImproperTorsion r : improperTorsion) {
                if (keep(r)) {
                    improperTorsions[nImproperTorsions++] = r;
                }
            }
            if (nImproperTorsions > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Improper Torsions:                 %10d", nImproperTorsions));
            }
        } else {
            nImproperTorsions = 0;
            improperTorsions = null;
        }

        if (vanderWaalsTerm) {
            if (esvTerm) {
                vanderWaals.setAtoms(esvSystem.getExtendedAtoms(), esvSystem.getExtendedMolecule());
            } else {
                vanderWaals.setAtoms(atoms, molecule);
            }
        }

        if (multipoleTerm) {
            if (esvTerm) {
                particleMeshEwald.setAtoms(esvSystem.getExtendedAtoms(), esvSystem.getExtendedMolecule());
            } else {
                particleMeshEwald.setAtoms(atoms, molecule);
            }
        }

        if (ncsTerm) {
            logger.severe(" NCS energy term cannot be used with variable systems sizes.");
        }

        if (restrainTerm) {
            logger.severe(" Restrain energy term cannot be used with variable systems sizes.");
        }

        if (comTerm) {
            logger.severe(" COM restrain energy term cannot be used with variable systems sizes.");
        }

        bondedRegion = new BondedRegion();
    }

    /**
     * <p>setFixedCharges.</p>
     *
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     */
    public void setFixedCharges(Atom[] atoms) {
        if (particleMeshEwald != null) {
            particleMeshEwald.setFixedCharges(atoms);
        }
    }

    /**
     * <p>Setter for the field <code>lambdaBondedTerms</code>.</p>
     *
     * @param lambdaBondedTerms a boolean.
     */
    void setLambdaBondedTerms(boolean lambdaBondedTerms) {
        this.lambdaBondedTerms = lambdaBondedTerms;
    }

    /**
     * <p>energy.</p>
     *
     * @return a double.
     */
    public double energy() {
        return energy(false, false);
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * energy</p>
     */
    public double energy(boolean gradient, boolean print) {

        try {
            bondTime = 0;
            angleTime = 0;
            stretchBendTime = 0;
            ureyBradleyTime = 0;
            outOfPlaneBendTime = 0;
            torsionTime = 0;
            stretchTorsionTime = 0;
            angleTorsionTime = 0;
            piOrbitalTorsionTime = 0;
            torsionTorsionTime = 0;
            improperTorsionTime = 0;
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
            angleTorsionEnergy = 0.0;
            stretchTorsionEnergy = 0.0;
            piOrbitalTorsionEnergy = 0.0;
            torsionTorsionEnergy = 0.0;
            improperTorsionEnergy = 0.0;
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
            permanentRealSpaceEnergy = 0.0;
            polarizationEnergy = 0.0;
            totalMultipoleEnergy = 0.0;
            totalNonBondedEnergy = 0.0;

            // Zero out the solvation energy.
            solvationEnergy = 0.0;

            // Zero out the relative solvation energy (sequence optimization)
            relativeSolvationEnergy = 0.0;
            nRelativeSolvations = 0;

            esvBias = 0.0;

            // Zero out the total potential energy.
            totalEnergy = 0.0;

            // Zero out the Cartesian coordinate gradient for each atom.
            if (gradient) {
                for (int i = 0; i < nAtoms; i++) {
                    atoms[i].setXYZGradient(0.0, 0.0, 0.0);
                    atoms[i].setLambdaXYZGradient(0.0, 0.0, 0.0);
                }
            }

            // Computed the bonded energy terms in parallel.
            try {
                bondedRegion.setGradient(gradient);
                parallelTeam.execute(bondedRegion);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception during bonded term calculation.");
                throw ex;
            } catch (Exception ex) {
                logger.info(Utilities.stackTraceToString(ex));
                logger.severe(ex.toString());
            }

            if (!lambdaBondedTerms) {
                // Compute restraint terms.
                if (ncsTerm) {
                    ncsTime = -System.nanoTime();
                    ncsEnergy = ncsRestraint.residual(gradient, print);
                    ncsTime += System.nanoTime();
                }
                if (restrainTerm && !coordRestraints.isEmpty()) {
                    coordRestraintTime = -System.nanoTime();
                    for (CoordRestraint restraint : coordRestraints) {
                        restrainEnergy += restraint.residual(gradient, print);
                    }
                    coordRestraintTime += System.nanoTime();
                }
                if (comTerm) {
                    comRestraintTime = -System.nanoTime();
                    comRestraintEnergy = comRestraint.residual(gradient, print);
                    comRestraintTime += System.nanoTime();
                }
                // Compute non-bonded terms.
                if (vanderWaalsTerm) {
                    vanDerWaalsTime = -System.nanoTime();
                    vanDerWaalsEnergy = vanderWaals.energy(gradient, print);
                    nVanDerWaalInteractions = this.vanderWaals.getInteractions();
                    vanDerWaalsTime += System.nanoTime();
                }
                if (multipoleTerm) {
                    electrostaticTime = -System.nanoTime();
                    totalMultipoleEnergy = particleMeshEwald.energy(gradient, print);
                    permanentMultipoleEnergy = particleMeshEwald.getPermanentEnergy();
                    permanentRealSpaceEnergy = particleMeshEwald.getPermRealEnergy();
                    polarizationEnergy = particleMeshEwald.getPolarizationEnergy();
                    nPermanentInteractions = particleMeshEwald.getInteractions();
                    solvationEnergy = particleMeshEwald.getGKEnergy();
                    nGKInteractions = particleMeshEwald.getGKInteractions();
                    electrostaticTime += System.nanoTime();
                }
            }

            if (relativeSolvationTerm) {
                List<Residue> residuesList = molecularAssembly.getResidueList();
                for (Residue residue : residuesList) {
                    if (residue instanceof MultiResidue) {
                        Atom refAtom = residue.getSideChainAtoms().get(0);
                        if (refAtom != null && refAtom.getUse()) {
                            // Reasonably confident that it should be -=,
                            // as we are trying to penalize residues with strong solvation energy.
                            double thisSolvation = relativeSolvation.getSolvationEnergy(residue, false);
                            relativeSolvationEnergy -= thisSolvation;
                            if (thisSolvation != 0) {
                                nRelativeSolvations++;
                            }
                        }
                    }
                }
            }

            totalTime = System.nanoTime() - totalTime;

            totalBondedEnergy = bondEnergy + restraintBondEnergy + angleEnergy
                    + stretchBendEnergy + ureyBradleyEnergy + outOfPlaneBendEnergy
                    + torsionEnergy + angleTorsionEnergy + stretchTorsionEnergy
                    + piOrbitalTorsionEnergy + improperTorsionEnergy
                    + torsionTorsionEnergy + ncsEnergy + restrainEnergy;
            totalNonBondedEnergy = vanDerWaalsEnergy + totalMultipoleEnergy + relativeSolvationEnergy;
            totalEnergy = totalBondedEnergy + totalNonBondedEnergy + solvationEnergy;
            if (esvTerm) {
                esvBias = esvSystem.getBiasEnergy();
                totalEnergy += esvBias;
            }
        } catch (EnergyException ex) {
            if (printOnFailure) {
                File origFile = molecularAssembly.getFile();
                String timeString = LocalDateTime.now().format(DateTimeFormatter.
                        ofPattern("yyyy_MM_dd-HH_mm_ss"));

                String filename = format("%s-ERROR-%s.pdb",
                        FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                        timeString);

                PotentialsFunctions ef = new PotentialsUtils();
                filename = ef.versionFile(filename);
                logger.info(format(" Writing on-error snapshot to file %s", filename));
                ef.saveAsPDB(molecularAssembly, new File(filename));
                molecularAssembly.setFile(origFile);
            }

            if (ex.doCauseSevere()) {
                logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
                return 0.0;
            } else {
                logger.log(Level.INFO, format(" Exception in energy calculation: %s", ex.toString()));
                throw ex; // Rethrow exception
            }
        }

        if (print) {
            StringBuilder sb = new StringBuilder();
            if (gradient) {
                sb.append("\n Computed Potential Energy and Atomic Coordinate Gradients\n");
            } else {
                sb.append("\n Computed Potential Energy\n");
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
     * Return the non-bonded components of energy (vdW, electrostatics).
     *
     * @param includeSolv Include solvation energy
     * @return Nonbonded energy
     */
    private double getNonbondedEnergy(boolean includeSolv) {
        return (includeSolv ? (totalNonBondedEnergy + solvationEnergy) : totalNonBondedEnergy);
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
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "BOND STRETCHING            : ", bondEnergy, nBonds));
            sb.append(format("REMARK   3   %s %g\n",
                    "BOND RMSD                  : ", bondRMSD));
        }
        if (angleTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "ANGLE BENDING              : ", angleEnergy, nAngles));
            sb.append(format("REMARK   3   %s %g\n",
                    "ANGLE RMSD                 : ", angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "STRETCH-BEND               : ", stretchBendEnergy, nStretchBends));
        }
        if (ureyBradleyTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "UREY-BRADLEY               : ", ureyBradleyEnergy, nUreyBradleys));
        }
        if (outOfPlaneBendTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "OUT-OF-PLANE BEND          : ", outOfPlaneBendEnergy, nOutOfPlaneBends));
        }
        if (torsionTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "TORSIONAL ANGLE            : ", torsionEnergy, nTorsions));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "PI-ORBITAL TORSION         : ", piOrbitalTorsionEnergy, nPiOrbitalTorsions));
        }
        if (torsionTorsionTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "TORSION-TORSION            : ", torsionTorsionEnergy, nTorsionTorsions));
        }
        if (improperTorsionTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "IMPROPER TORSION           : ", improperTorsionEnergy, nImproperTorsions));
        }
        if (restraintBondTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "RESTRAINT BOND STRETCHING            : ", restraintBondEnergy, nRestraintBonds));
        }

        if (ncsTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "NCS RESTRAINT              : ", ncsEnergy, nAtoms));
        }

        if (restrainTerm && !coordRestraints.isEmpty()) {
            int nRests = 0;
            for (CoordRestraint restraint : coordRestraints) {
                nRests += restraint.getNumAtoms();
            }
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "COORDINATE RESTRAINTS      : ", restrainEnergy, nRests));
        }

        if (comTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "COM RESTRAINT              : ", comRestraintEnergy, nAtoms));
        }

        if (vanderWaalsTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "VAN DER WAALS              : ", vanDerWaalsEnergy, nVanDerWaalInteractions));
        }
        if (multipoleTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "ATOMIC MULTIPOLES          : ", permanentMultipoleEnergy, nPermanentInteractions));
        }
        if (polarizationTerm) {
            sb.append(format("REMARK   3   %s %g (%d)\n",
                    "POLARIZATION               : ", polarizationEnergy, nPermanentInteractions));
        }
        sb.append(format("REMARK   3   %s %g\n",
                "TOTAL POTENTIAL (KCAL/MOL) : ", totalEnergy));
        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(format("REMARK   3   %s %g\n",
                    "UNIT CELL POTENTIAL        : ", totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(format("REMARK   3   %s %g\n",
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
        StringBuilder sb = new StringBuilder();

        if (bondTerm && nBonds > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f (%8.5f)\n",
                    "Bond Stretching   ", bondEnergy, nBonds,
                    bondTime * toSeconds, bondRMSD));
        }
        if (angleTerm && nAngles > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f (%8.5f)\n",
                    "Angle Bending     ", angleEnergy, nAngles,
                    angleTime * toSeconds, angleRMSD));
        }
        if (stretchBendTerm && nStretchBends > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Stretch-Bend      ", stretchBendEnergy,
                    nStretchBends, stretchBendTime * toSeconds));
        }
        if (ureyBradleyTerm && nUreyBradleys > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Urey-Bradley      ", ureyBradleyEnergy,
                    nUreyBradleys, ureyBradleyTime * toSeconds));
        }
        if (outOfPlaneBendTerm && nOutOfPlaneBends > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Out-of-Plane Bend ", outOfPlaneBendEnergy,
                    nOutOfPlaneBends, outOfPlaneBendTime * toSeconds));
        }
        if (torsionTerm && nTorsions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Torsional Angle   ", torsionEnergy, nTorsions,
                    torsionTime * toSeconds));
        }
        if (piOrbitalTorsionTerm && nPiOrbitalTorsions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Pi-Orbital Torsion", piOrbitalTorsionEnergy,
                    nPiOrbitalTorsions, piOrbitalTorsionTime * toSeconds));
        }
        if (stretchTorsionTerm && nStretchTorsions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Stretch-Torsion   ", stretchTorsionEnergy, nStretchTorsions,
                    stretchTorsionTime * toSeconds));
        }
        if (angleTorsionTerm && nAngleTorsions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Angle-Torsion     ", angleTorsionEnergy, nAngleTorsions,
                    angleTorsionTime * toSeconds));
        }
        if (torsionTorsionTerm && nTorsionTorsions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Torsion-Torsion   ", torsionTorsionEnergy,
                    nTorsionTorsions, torsionTorsionTime * toSeconds));
        }
        if (improperTorsionTerm && nImproperTorsions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Improper Torsion  ", improperTorsionEnergy,
                    nImproperTorsions, improperTorsionTime * toSeconds));
        }

        if (restraintBondTerm && nRestraintBonds > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Bond Restraint    ", restraintBondEnergy, nRestraintBonds,
                    restraintBondTime * toSeconds));
        }
        if (ncsTerm) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "NCS Restraint     ", ncsEnergy, nAtoms,
                    ncsTime * toSeconds));
        }
        if (restrainTerm && !coordRestraints.isEmpty()) {
            int nRests = 0;
            for (CoordRestraint restraint : coordRestraints) {
                nRests += restraint.getNumAtoms();
            }
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Coord. Restraints ", restrainEnergy, nRests,
                    coordRestraintTime * toSeconds));
        }
        if (comTerm) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "COM Restraint     ", comRestraintEnergy, nAtoms,
                    comRestraintTime * toSeconds));
        }
        if (vanderWaalsTerm && nVanDerWaalInteractions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Van der Waals     ", vanDerWaalsEnergy,
                    nVanDerWaalInteractions, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm && nPermanentInteractions > 0) {
            String pmeTitle = (particleMeshEwald instanceof ParticleMeshEwaldQI)
                    ? "Q.Int. Multipoles "
                    : "Atomic Multipoles ";
            if (polarizationTerm) {
                sb.append(format("  %s %16.8f %12d\n",
                        pmeTitle, permanentMultipoleEnergy, nPermanentInteractions));
            } else {
                sb.append(format("  %s %16.8f %12d %12.3f\n",
                        pmeTitle, permanentMultipoleEnergy, nPermanentInteractions, electrostaticTime * toSeconds));
            }
        }
        if (polarizationTerm && nPermanentInteractions > 0) {
            sb.append(format("  %s %16.8f %12d %12.3f\n",
                    "Polarization      ", polarizationEnergy,
                    nPermanentInteractions, electrostaticTime * toSeconds));
        }
        if (generalizedKirkwoodTerm && nGKInteractions > 0) {
            sb.append(format("  %s %16.8f %12d\n",
                    "Solvation         ", solvationEnergy, nGKInteractions));
        }

        if (relativeSolvationTerm) {
            sb.append(format("  %s %16.8f %12d\n",
                    "Relative Solvation", relativeSolvationEnergy, nRelativeSolvations));
        }

        if (esvTerm) {
            sb.append(format("  %s %16.8f  %s\n",
                    "ExtendedSystemBias", esvBias, esvSystem.getLambdaList()));
            sb.append(esvSystem.getBiasDecomposition());
        }

        sb.append(format("  %s %16.8f  %s %12.3f (sec)",
                "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));

        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(format("\n  %s %16.8f", "Unit Cell         ",
                    totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(format("\n  %s %16.8f", "Replicates Cell   ",
                    totalEnergy * nsymm));
        }
        return sb.toString();
    }

    /**
     * <p>Getter for the field <code>parallelTeam</code>.</p>
     *
     * @return a {@link edu.rit.pj.ParallelTeam} object.
     */
    public ParallelTeam getParallelTeam() {
        return parallelTeam;
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * Getter for the field <code>crystal</code>.</p>
     */
    @Override
    public Crystal getCrystal() {
        return crystal;
    }

    /**
     * <p>Getter for the field <code>cutoffPlusBuffer</code>.</p>
     *
     * @return a double.
     */
    public double getCutoffPlusBuffer() {
        return cutoffPlusBuffer;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            if (lambda <= 1.0 && lambda >= 0.0) {
                this.lambda = lambda;
                if (vanderWaalsTerm) {
                    vanderWaals.setLambda(lambda);
                }
                if (multipoleTerm) {
                    particleMeshEwald.setLambda(lambda);
                }
                if (restraintBondTerm && restraintBonds != null) {
                    for (RestraintBond restraintBond : restraintBonds) {
                        restraintBond.setLambda(lambda);
                    }
                }
                if (ncsTerm && ncsRestraint != null) {
                    ncsRestraint.setLambda(lambda);
                }
                if (restrainTerm && !coordRestraints.isEmpty()) {
                    for (CoordRestraint restraint : coordRestraints) {
                        restraint.setLambda(lambda);
                    }
                }
                if (comTerm && comRestraint != null) {
                    comRestraint.setLambda(lambda);
                }
                if (lambdaTorsions) {
                    for (int i = 0; i < nTorsions; i++) {
                        torsions[i].setLambda(lambda);
                    }
                    for (int i = 0; i < nPiOrbitalTorsions; i++) {
                        piOrbitalTorsions[i].setLambda(lambda);
                    }
                    for (int i = 0; i < nTorsionTorsions; i++) {
                        torsionTorsions[i].setLambda(lambda);
                    }
                }
            } else {
                String message = format("Lambda value %8.3f is not in the range [0..1].", lambda);
                logger.warning(message);
            }
        } else {
            logger.fine(" Attempting to set a lambda value on a ForceFieldEnergy with lambdaterm false.");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setScaling(double[] scaling) {
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
     * {@inheritDoc}
     * <p>
     * Return a reference to each variables type.
     */
    @Override
    public VARIABLE_TYPE[] getVariableTypes() {
        int n = getNumberOfVariables();
        VARIABLE_TYPE[] type = new VARIABLE_TYPE[n];
        int i = 0;
        for (int j = 0; j < nAtoms; j++) {
            if (atoms[j].isActive()) {
                type[i++] = VARIABLE_TYPE.X;
                type[i++] = VARIABLE_TYPE.Y;
                type[i++] = VARIABLE_TYPE.Z;
            }
        }
        return type;
    }

    /**
     * Returns a copy of the list of constraints this ForceFieldEnergy has.
     *
     * @return Copied list of constraints.
     */
    @Override
    public List<Constraint> getConstraints() {
        return constraints.isEmpty() ? Collections.emptyList() : new ArrayList<>(constraints);
    }

    /**
     * Applies constraints to positions
     * @param xPrior
     * @param xNew
     */
    public void applyAllConstraintPositions(double[] xPrior, double[] xNew) {
        applyAllConstraintPositions(xPrior, xNew, DEFAULT_CONSTRAINT_TOLERANCE);
    }

    public void applyAllConstraintPositions(double[] xPrior, double[] xNew, double tol) {
        if (xPrior == null) {
            xPrior = Arrays.copyOf(xNew, xNew.length);
        }
        for (Constraint constraint : constraints) {
            constraint.applyConstraintToStep(xPrior, xNew, getMass(), tol);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x) {
        return energy(x, false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x, boolean verbose) {
        assert Arrays.stream(x).allMatch(Double::isFinite);

        // Unscale the coordinates.
        unscaleCoordinates(x);

        // Set coordinates.
        setCoordinates(x);

        double e = this.energy(false, verbose);

        // Rescale the coordinates.
        scaleCoordinates(x);

        return e;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double[] x, double[] g) {
        return energyAndGradient(x, g, false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double[] x, double[] g, boolean verbose) {
        assert Arrays.stream(x).allMatch(Double::isFinite);

        // Un-scale the coordinates.
        unscaleCoordinates(x);

        // Set coordinates.
        setCoordinates(x);
        double e = energy(true, verbose);

        // Try block already exists inside energy(boolean, boolean), so only
        // need to try-catch fillGradient.
        try {
            fillGradient(g);

            // Scale the coordinates and gradients.
            scaleCoordinatesAndGradient(x, g);

            if (maxDebugGradient < Double.MAX_VALUE) {
                boolean extremeGrad = Arrays.stream(g).anyMatch((double gi) ->
                        (gi > maxDebugGradient || gi < -maxDebugGradient));
                if (extremeGrad) {
                    File origFile = molecularAssembly.getFile();
                    String timeString = LocalDateTime.now().format(DateTimeFormatter.
                            ofPattern("yyyy_MM_dd-HH_mm_ss"));

                    String filename = format("%s-LARGEGRAD-%s.pdb",
                            FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                            timeString);
                    PotentialsFunctions ef = new PotentialsUtils();
                    filename = ef.versionFile(filename);

                    logger.warning(format(" Excessively large gradient detected; printing snapshot to file %s", filename));
                    ef.saveAsPDB(molecularAssembly, new File(filename));
                    molecularAssembly.setFile(origFile);
                }
            }
            return e;
        } catch (EnergyException ex) {
            if (printOnFailure) {
                logger.info(Utilities.stackTraceToString(ex));
                String timeString = LocalDateTime.now().format(DateTimeFormatter.
                        ofPattern("yyyy_MM_dd-HH_mm_ss"));

                String filename = format("%s-ERROR-%s.pdb",
                        FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                        timeString);

                PotentialsFunctions ef = new PotentialsUtils();
                filename = ef.versionFile(filename);
                logger.info(format(" Writing on-error snapshot to file %s", filename));
                ef.saveAsPDB(molecularAssembly, new File(filename));
            }

            if (ex.doCauseSevere()) {
                logger.info(Utilities.stackTraceToString(ex));
                logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
            } else {
                logger.log(Level.INFO, format(" Exception in energy calculation: %s", ex.toString()));
            }

            throw ex;
        }
    }

    /**
     * <p>
     * getGradient</p>
     *
     * @param g an array of double.
     * @return an array of {@link double} objects.
     */
    public double[] getGradient(double[] g) {
        return fillGradient(g);
    }

    /**
     * Private method for internal use, so we don't have subclasses calling super.energy, and this class delegating to
     * the subclass's getGradient method.
     *
     * @param g Gradient array to fill.
     * @return Gradient array.
     */
    private double[] fillGradient(double[] g) {
        assert (g != null);
        double[] grad = new double[3];
        int n = getNumberOfVariables();
        if (g==null || g.length < n) {
            g = new double[n];
        }
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                a.getXYZGradient(grad);
                double gx = grad[0];
                double gy = grad[1];
                double gz = grad[2];
                if (isNaN(gx) || isInfinite(gx) || isNaN(gy) || isInfinite(gy) || isNaN(gz) || isInfinite(gz)) {
                    StringBuilder sb = new StringBuilder(format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                            a.toString(), gx, gy, gz));
                    double[] vals = new double[3];
                    a.getVelocity(vals);
                    sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getAcceleration(vals);
                    sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getPreviousAcceleration(vals);
                    sb.append(format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));

                    throw new EnergyException(sb.toString());
                }
                g[index++] = gx;
                g[index++] = gy;
                g[index++] = gz;
            }
        }
        return g;
    }

    /**
     * The coordinate array should only contain active atoms.
     *
     * @param coords an array of {@link double} objects.
     */
    public void setCoordinates(double[] coords) {
        if (coords == null) {
            return;
        }
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive() && !a.isBackground()) {
                double x = coords[index++];
                double y = coords[index++];
                double z = coords[index++];
                a.moveTo(x, y, z);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getCoordinates(double[] x) {
        int n = getNumberOfVariables();
        if (x == null || x.length < n) {
            x = new double[n];
        }
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive() && !a.isBackground()) {
                x[index++] = a.getX();
                x[index++] = a.getY();
                x[index++] = a.getZ();
            }
        }
        return x;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double[] getMass() {
        int n = getNumberOfVariables();
        double[] mass = new double[n];
        int index = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                double m = a.getMass();
                mass[index++] = m;
                mass[index++] = m;
                mass[index++] = m;
            }
        }
        return mass;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNumberOfVariables() {
        int nActive = 0;
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive() && !atoms[i].isBackground()) {
                nActive++;
            }
        }
        return nActive * 3;
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
            if (restrainTerm && !coordRestraints.isEmpty()) {
                for (CoordRestraint restraint : coordRestraints) {
                    dEdLambda += restraint.getdEdL();
                }
            }
            if (comTerm && comRestraint != null) {
                dEdLambda += comRestraint.getdEdL();
            }
            if (lambdaTorsions) {
                for (int i = 0; i < nTorsions; i++) {
                    dEdLambda += torsions[i].getdEdL();
                }
                for (int i = 0; i < nPiOrbitalTorsions; i++) {
                    dEdLambda += piOrbitalTorsions[i].getdEdL();
                }
                for (int i = 0; i < nTorsionTorsions; i++) {
                    dEdLambda += torsionTorsions[i].getdEdL();
                }
            }
        }
        return dEdLambda;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns true if lambda term is not enabled for this ForceFieldEnergy.
     */
    @Override
    public boolean dEdLZeroAtEnds() {
        // This may actually be true even with softcored atoms.
        // For now, serves the purpose of reporting true when nothing is softcored.
        return !lambdaTerm;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradients) {
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
            if (restrainTerm && !coordRestraints.isEmpty()) {
                for (CoordRestraint restraint : coordRestraints) {
                    restraint.getdEdXdL(gradients);
                }
            }
            if (comTerm && comRestraint != null) {
                comRestraint.getdEdXdL(gradients);
            }
            if (lambdaTorsions) {
                double[] grad = new double[3];
                int index = 0;
                for (int i = 0; i < nAtoms; i++) {
                    Atom a = atoms[i];
                    if (a.isActive()) {
                        a.getLambdaXYZGradient(grad);
                        gradients[index++] += grad[0];
                        gradients[index++] += grad[1];
                        gradients[index++] += grad[2];
                    }
                }
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
            if (restrainTerm && !coordRestraints.isEmpty()) {
                for (CoordRestraint restraint : coordRestraints) {
                    d2EdLambda2 += restraint.getd2EdL2();
                }
            }
            if (comTerm && comRestraint != null) {
                d2EdLambda2 += comRestraint.getd2EdL2();
            }
            if (lambdaTorsions) {
                for (int i = 0; i < nTorsions; i++) {
                    d2EdLambda2 += torsions[i].getd2EdL2();
                }
                for (int i = 0; i < nPiOrbitalTorsions; i++) {
                    d2EdLambda2 += piOrbitalTorsions[i].getd2EdL2();
                }
                for (int i = 0; i < nTorsionTorsions; i++) {
                    d2EdLambda2 += torsionTorsions[i].getd2EdL2();
                }
            }
        }
        return d2EdLambda2;
    }

    /**
     * Sets the printOnFailure flag; if override is true, over-rides any
     * existing property. Essentially sets the default value of printOnFailure
     * for an algorithm. For example, rotamer optimization will generally run
     * into force field issues in the normal course of execution as it tries
     * unphysical self and pair configurations, so the algorithm should not
     * print out a large number of error PDBs.
     *
     * @param onFail   To set
     * @param override Override properties
     */
    public void setPrintOnFailure(boolean onFail, boolean override) {
        if (override) {
            // Ignore any pre-existing value
            printOnFailure = onFail;
        } else {
            try {
                molecularAssembly.getForceField().getBoolean(ForceFieldBoolean.PRINT_ON_FAILURE);
                /*
                  If the call was successful, the property was explicitly set
                  somewhere and should be kept. If an exception was thrown, the
                  property was never set explicitly, so over-write.
                 */
            } catch (Exception ex) {
                printOnFailure = onFail;
            }
        }
    }

    /**
     * Gets the Platform associated with this force field energy. For the reference platform, always returns FFX.
     *
     * @return A Platform.
     */
    public Platform getPlatform() {
        return platform;
    }

    /**
     * <p>
     * setRestraintBond</p>
     *
     * @param a1            a {@link ffx.potential.bonded.Atom} object.
     * @param a2            a {@link ffx.potential.bonded.Atom} object.
     * @param distance      a double.
     * @param forceConstant the force constant in kcal/mole
     */
    public void setRestraintBond(Atom a1, Atom a2, double distance, double forceConstant) {
        setRestraintBond(a1, a2, distance, forceConstant, 0);
    }

    /**
     * <p>
     * setRestraintBond</p>
     *
     * @param a1            a {@link ffx.potential.bonded.Atom} object.
     * @param a2            a {@link ffx.potential.bonded.Atom} object.
     * @param distance      a double.
     * @param forceConstant the force constant in kcal/mole.
     * @param flatBottom    Radius of a flat-bottom potential in Angstroms.
     */
    private void setRestraintBond(Atom a1, Atom a2, double distance, double forceConstant, double flatBottom) {
        restraintBondTerm = true;
        RestraintBond rb = new RestraintBond(a1, a2, crystal);
        int[] classes = {a1.getAtomType().atomClass, a2.getAtomType().atomClass};
        if (flatBottom != 0) {
            rb.setBondType(new BondType(classes, forceConstant, distance, BondType.BondFunction.FLAT_BOTTOM_HARMONIC, flatBottom));
        } else {
            rb.setBondType((new BondType(classes, forceConstant, distance, BondType.BondFunction.HARMONIC)));
        }
        // As long as we continue to add elements one-at-a-time to an array, this code will continue to be ugly.
        RestraintBond[] newRbs = new RestraintBond[++nRestraintBonds];
        if (restraintBonds != null && restraintBonds.length != 0) {
            System.arraycopy(restraintBonds, 0, newRbs, 0, (nRestraintBonds - 1));
        }
        newRbs[nRestraintBonds - 1] = rb;
        restraintBonds = newRbs;
        rb.energy(false);
        rb.log();
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
     * <p>
     * This method is for the RESPA integrator only.
     */
    @Override
    public void setEnergyTermState(STATE state) {
        this.state = state;
        switch (state) {
            case FAST:
                bondTerm = bondTermOrig;
                angleTerm = angleTermOrig;
                stretchBendTerm = stretchBendTermOrig;
                ureyBradleyTerm = ureyBradleyTermOrig;
                outOfPlaneBendTerm = outOfPlaneBendTermOrig;
                torsionTerm = torsionTermOrig;
                stretchTorsionTerm = stretchTorsionTermOrig;
                angleTorsionTerm = angleTorsionTermOrig;
                piOrbitalTorsionTerm = piOrbitalTorsionTermOrig;
                torsionTorsionTerm = torsionTorsionTermOrig;
                improperTorsionTerm = improperTorsionTermOrig;
                restraintBondTerm = restraintBondTermOrig;
                ncsTerm = ncsTermOrig;
                restrainTerm = restrainTermOrig;
                comTerm = comTermOrig;
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
                stretchTorsionTerm = false;
                angleTorsionTerm = false;
                piOrbitalTorsionTerm = false;
                torsionTorsionTerm = false;
                improperTorsionTerm = false;
                restraintBondTerm = false;
                ncsTerm = false;
                restrainTerm = false;
                comTerm = false;
                break;
            default:
                bondTerm = bondTermOrig;
                angleTerm = angleTermOrig;
                stretchBendTerm = stretchBendTermOrig;
                ureyBradleyTerm = ureyBradleyTermOrig;
                outOfPlaneBendTerm = outOfPlaneBendTermOrig;
                torsionTerm = torsionTermOrig;
                stretchTorsionTerm = stretchTorsionTermOrig;
                angleTorsionTerm = angleTorsionTermOrig;
                piOrbitalTorsionTerm = piOrbitalTorsionTermOrig;
                torsionTorsionTerm = torsionTorsionTermOrig;
                improperTorsionTerm = improperTorsionTermOrig;
                restraintBondTerm = restraintBondTermOrig;
                ncsTerm = ncsTermOrig;
                restrainTermOrig = restrainTerm;
                comTermOrig = comTerm;
                vanderWaalsTerm = vanderWaalsTermOrig;
                multipoleTerm = multipoleTermOrig;
                polarizationTerm = polarizationTermOrig;
                generalizedKirkwoodTerm = generalizedKirkwoodTermOrig;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Set the boundary conditions for this calculation.
     */
    @Override
    public void setCrystal(Crystal crystal) {
        this.crystal = crystal;
        /*
          Update VanDerWaals first, in case the NeighborList needs to be
          re-allocated to include a larger number of replicated cells.
         */
        if (vanderWaalsTerm) {
            vanderWaals.setCrystal(crystal);
        }
        if (multipoleTerm) {
            particleMeshEwald.setCrystal(crystal);
        }
    }


    private Crystal configureNCS(ForceField forceField, Crystal unitCell) {
        // MTRIXn to be permuted with standard space group in NCSCrystal.java for experimental refinement.
        if (forceField.getProperties().containsKey("MTRIX1") && forceField.getProperties().containsKey("MTRIX2") && forceField.getProperties().containsKey("MTRIX3")) {
            Crystal unitCell2 = new Crystal(unitCell.a, unitCell.b, unitCell.c, unitCell.alpha,
                    unitCell.beta, unitCell.gamma, unitCell.spaceGroup.pdbName);
            SpaceGroup spaceGroup = unitCell2.spaceGroup;
            // Separate string list MTRIXn into Double matricies then pass into symops
            CompositeConfiguration properties = forceField.getProperties();
            String[] MTRX1List = properties.getStringArray("MTRIX1");
            String[] MTRX2List = properties.getStringArray("MTRIX2");
            String[] MTRX3List = properties.getStringArray("MTRIX3");
            spaceGroup.symOps.clear();
            double number1;
            double number2;
            double number3;
            for (int i = 0; i < MTRX1List.length; i++) {
                double[][] Rot_MTRX = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                double[] Tr_MTRX = {0, 0, 0};
                String[] tokens1 = MTRX1List[i].trim().split(" +"); // 4 items: rot [0][1-3] * trans[0]
                String[] tokens2 = MTRX2List[i].trim().split(" +"); // 4 items: rot [1][1-3] * trans[1]
                String[] tokens3 = MTRX3List[i].trim().split(" +"); // 4 items: rot [2][1-3] * trans[2]
                for (int k = 0; k < 4; k++) {
                    number1 = Double.parseDouble(tokens1[k]);
                    number2 = Double.parseDouble(tokens2[k]);
                    number3 = Double.parseDouble(tokens3[k]);
                    if (k != 3) {
                        Rot_MTRX[0][k] = number1;
                        Rot_MTRX[1][k] = number2;
                        Rot_MTRX[2][k] = number3;
                    } else {
                        Tr_MTRX[0] = number1;
                        Tr_MTRX[1] = number2;
                        Tr_MTRX[2] = number3;
                    }
                }
                SymOp symOp = new SymOp(Rot_MTRX, Tr_MTRX);
                if (logger.isLoggable(Level.FINEST)) {
                    logger.info(format(" MTRIXn SymOp: %d of %d\n" + symOp.toString(), i + 1, MTRX1List.length));
                }
                spaceGroup.symOps.add(symOp);
            }
            unitCell = NCSCrystal.NCSCrystalFactory(unitCell, spaceGroup.symOps);
            unitCell.updateCrystal();
        }
        return unitCell;
    }

    /**
     * Frees up assets associated with this ForceFieldEnergy, such as worker Threads.
     *
     * @return If successful in freeing up assets.
     */
    public boolean destroy() {
        if (destroyed) {
            logger.info(format(" This ForceFieldEnergy is already destroyed: %s", this.toString()));
            return true;
        } else {
            try {
                if (parallelTeam != null) {
                    parallelTeam.shutdown();
                }
                if (vanderWaals != null) {
                    vanderWaals.destroy();
                }
                if (particleMeshEwald != null) {
                    particleMeshEwald.destroy();
                }
                molecularAssembly.finishDestruction();
                destroyed = true;
                return true;
            } catch (Exception ex) {
                logger.warning(format(" Exception in shutting down a ForceFieldEnergy: %s", ex));
                logger.info(Utilities.stackTraceToString(ex));
                return false;
            }
        }
    }

    /**
     * <p>Getter for the field <code>bondEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getBondEnergy() {
        return bondEnergy;
    }

    /**
     * <p>getNumberofBonds.</p>
     *
     * @return a int.
     */
    public int getNumberofBonds() {
        return nBonds;
    }

    /**
     * <p>Getter for the field <code>angleEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getAngleEnergy() {
        return angleEnergy;
    }

    /**
     * <p>getNumberofAngles.</p>
     *
     * @return a int.
     */
    public int getNumberofAngles() {
        return nAngles;
    }

    /**
     * <p>getStrenchBendEnergy.</p>
     *
     * @return a double.
     */
    public double getStrenchBendEnergy() {
        return stretchBendEnergy;
    }

    /**
     * <p>getNumberofStretchBends.</p>
     *
     * @return a int.
     */
    public int getNumberofStretchBends() {
        return nStretchBends;
    }

    /**
     * <p>Getter for the field <code>ureyBradleyEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getUreyBradleyEnergy() {
        return ureyBradleyEnergy;
    }

    /**
     * <p>getNumberofUreyBradleys.</p>
     *
     * @return a int.
     */
    public int getNumberofUreyBradleys() {
        return nUreyBradleys;
    }

    /**
     * <p>Getter for the field <code>outOfPlaneBendEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getOutOfPlaneBendEnergy() {
        return outOfPlaneBendEnergy;
    }

    /**
     * <p>getNumberofOutOfPlaneBends.</p>
     *
     * @return a int.
     */
    public int getNumberofOutOfPlaneBends() {
        return nOutOfPlaneBends;
    }

    /**
     * <p>Getter for the field <code>torsionEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getTorsionEnergy() {
        return torsionEnergy;
    }

    /**
     * <p>getNumberofTorsions.</p>
     *
     * @return a int.
     */
    public int getNumberofTorsions() {
        return nTorsions;
    }

    /**
     * <p>Getter for the field <code>stretchTorsionEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getStretchTorsionEnergy() {
        return stretchTorsionEnergy;
    }

    /**
     * <p>getNumberofStretchTorsions.</p>
     *
     * @return a int.
     */
    public int getNumberofStretchTorsions() {
        return nStretchTorsions;
    }

    /**
     * <p>Getter for the field <code>angleTorsionEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getAngleTorsionEnergy() {
        return angleTorsionEnergy;
    }

    /**
     * <p>getNumberofAngleTorsions.</p>
     *
     * @return a int.
     */
    public int getNumberofAngleTorsions() {
        return nAngleTorsions;
    }

    /**
     * <p>Getter for the field <code>improperTorsionEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getImproperTorsionEnergy() {
        return improperTorsionEnergy;
    }

    /**
     * <p>getNumberofImproperTorsions.</p>
     *
     * @return a int.
     */
    public int getNumberofImproperTorsions() {
        return nImproperTorsions;
    }

    /**
     * <p>Getter for the field <code>piOrbitalTorsionEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getPiOrbitalTorsionEnergy() {
        return piOrbitalTorsionEnergy;
    }

    /**
     * <p>getNumberofPiOrbitalTorsions.</p>
     *
     * @return a int.
     */
    public int getNumberofPiOrbitalTorsions() {
        return nPiOrbitalTorsions;
    }

    /**
     * <p>Getter for the field <code>torsionTorsionEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getTorsionTorsionEnergy() {
        return torsionTorsionEnergy;
    }

    /**
     * <p>getNumberofTorsionTorsions.</p>
     *
     * @return a int.
     */
    public int getNumberofTorsionTorsions() {
        return nTorsionTorsions;
    }

    /**
     * <p>Getter for the field <code>vanDerWaalsEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getVanDerWaalsEnergy() {
        return vanDerWaalsEnergy;
    }

    /**
     * <p>getVanDerWaalsInteractions.</p>
     *
     * @return a int.
     */
    public int getVanDerWaalsInteractions() {
        return nVanDerWaalInteractions;
    }

    /**
     * <p>setLambdaMultipoleScale.</p>
     *
     * @param scale a double.
     */
    public void setLambdaMultipoleScale(double scale) {
        if (particleMeshEwald != null) {
            particleMeshEwald.setLambdaMultipoleScale(scale);
        }
    }

    /**
     * <p>getEnergyComponent.</p>
     *
     * @param component a {@link ffx.potential.PotentialComponent} object.
     * @return a double.
     */
    public double getEnergyComponent(PotentialComponent component) {
        switch (component) {
            case Topology:
            case ForceFieldEnergy:
                return totalEnergy;
            case VanDerWaals:
                return vanDerWaalsEnergy;
            case Multipoles:
                return particleMeshEwald.getTotalMultipoleEnergy();
            case Permanent:
                return particleMeshEwald.getPermanentEnergy();
            case Induced:
                return particleMeshEwald.getPolarizationEnergy();
            case PermanentRealSpace:
                return particleMeshEwald.getPermRealEnergy();
            case PermanentSelf:
                return particleMeshEwald.getPermSelfEnergy();
            case PermanentReciprocal:
                return particleMeshEwald.getPermRecipEnergy();
            case InducedRealSpace:
                return particleMeshEwald.getIndRealEnergy();
            case InducedSelf:
                return particleMeshEwald.getIndSelfEnergy();
            case InducedReciprocal:
                return particleMeshEwald.getIndRecipEnergy();
            case GeneralizedKirkwood:
                return particleMeshEwald.getGKEnergy();
            case Bonded:
                return totalBondedEnergy;
            case Bond:
                return bondEnergy;
            case Angle:
                return angleEnergy;
            case Torsion:
                return torsionEnergy;
            case StretchBend:
                return stretchBendEnergy;
            case OutOfPlaneBend:
                return outOfPlaneBendEnergy;
            case PiOrbitalTorsion:
                return piOrbitalTorsionEnergy;
            case TorsionTorsion:
                return torsionTorsionEnergy;
            case UreyBradley:
                return ureyBradleyEnergy;
            case RestraintBond:
                return restraintBondEnergy;
            case ImproperTorsion:
                return improperTorsionEnergy;
            case NCS:
                return ncsEnergy;
            case Restrain:
                return restrainEnergy;
            case pHMD:
            case Bias:
            case Discretizer:
            case Acidostat:
                return (esvTerm) ? esvSystem.getEnergyComponent(component) : 0.0;
            case XRay:
            default:
                throw new AssertionError(component.name());
        }
    }

    /**
     * <p>Getter for the field <code>permanentMultipoleEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getPermanentMultipoleEnergy() {
        return permanentMultipoleEnergy;
    }

    /**
     * <p>Getter for the field <code>permanentRealSpaceEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getPermanentRealSpaceEnergy() {
        return permanentRealSpaceEnergy;
    }

    /**
     * <p>getPermanentReciprocalSelfEnergy.</p>
     *
     * @return a double.
     */
    public double getPermanentReciprocalSelfEnergy() {
        return particleMeshEwald.getPermSelfEnergy();
    }

    /**
     * <p>getPermanentReciprocalMpoleEnergy.</p>
     *
     * @return a double.
     */
    public double getPermanentReciprocalMpoleEnergy() {
        return particleMeshEwald.getPermRecipEnergy();
    }

    /**
     * <p>getPermanentInteractions.</p>
     *
     * @return a int.
     */
    public int getPermanentInteractions() {
        return nPermanentInteractions;
    }

    /**
     * <p>Getter for the field <code>polarizationEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    /**
     * <p>getTotalElectrostaticEnergy.</p>
     *
     * @return a double.
     */
    public double getTotalElectrostaticEnergy() {
        return totalMultipoleEnergy + solvationEnergy;
    }

    /**
     * <p>getElectrostaticEnergy.</p>
     *
     * @return a double.
     */
    public double getElectrostaticEnergy() {
        return totalMultipoleEnergy;
    }

    /**
     * <p>Getter for the field <code>solvationEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getSolvationEnergy() {
        return solvationEnergy;
    }

    /**
     * <p>getCavitationEnergy.</p>
     *
     * @return a double.
     */
    public double getCavitationEnergy() {
        return particleMeshEwald.getCavitationEnergy(false);
    }

    /**
     * <p>getDispersionEnergy.</p>
     *
     * @return a double.
     */
    public double getDispersionEnergy() {
        return particleMeshEwald.getDispersionEnergy(false);
    }

    /**
     * <p>getEsvBiasEnergy.</p>
     *
     * @return a double.
     */
    public double getEsvBiasEnergy() {
        return esvBias;
    }

    /**
     * <p>getExtendedSystem.</p>
     *
     * @return a {@link ffx.potential.extended.ExtendedSystem} object.
     */
    public ExtendedSystem getExtendedSystem() {
        return esvSystem;
    }

    /**
     * <p>getGK.</p>
     *
     * @return a {@link ffx.potential.nonbonded.GeneralizedKirkwood} object.
     */
    public GeneralizedKirkwood getGK() {
        if (particleMeshEwald != null) {
            return particleMeshEwald.getGK();
        } else {
            return null;
        }
    }

    /**
     * <p>getVdwNode.</p>
     *
     * @return a {@link ffx.potential.nonbonded.VanDerWaals} object.
     */
    public VanDerWaals getVdwNode() {
        return vanderWaals;
    }

    /**
     * <p>getPmeNode.</p>
     *
     * @return a {@link ffx.potential.nonbonded.ParticleMeshEwald} object.
     */
    public ParticleMeshEwald getPmeNode() {
        return particleMeshEwald;
    }

    /**
     * <p>getPmeCartNode.</p>
     *
     * @return a {@link ffx.potential.nonbonded.ParticleMeshEwaldCart} object.
     */
    public ParticleMeshEwaldCart getPmeCartNode() {
        return (particleMeshEwald instanceof ParticleMeshEwaldCart)
                ? (ParticleMeshEwaldCart) particleMeshEwald : null;
    }

    /**
     * <p>getPmeQiNode.</p>
     *
     * @return a {@link ffx.potential.nonbonded.ParticleMeshEwaldQI} object.
     */
    public ParticleMeshEwaldQI getPmeQiNode() {
        return (particleMeshEwald instanceof ParticleMeshEwaldQI)
                ? (ParticleMeshEwaldQI) particleMeshEwald : null;
    }

    /**
     * <p>Getter for the field <code>coordRestraints</code>.</p>
     *
     * @return a {@link java.util.List} object.
     */
    public List<CoordRestraint> getCoordRestraints() {
        return new ArrayList<>(coordRestraints);
    }

    /**
     * <p>Getter for the field <code>restraintBonds</code>.</p>
     *
     * @return a {@link java.util.List} object.
     */
    protected List<RestraintBond> getRestraintBonds() {
        if (restraintBonds != null && restraintBonds.length > 0) {
            return Arrays.asList(restraintBonds);
        } else {
            return null;
        }
    }

    /**
     * <p>getSolvationInteractions.</p>
     *
     * @return a int.
     */
    public int getSolvationInteractions() {
        return nGKInteractions;
    }

    /**
     * <p>Getter for the field <code>relativeSolvationEnergy</code>.</p>
     *
     * @return a double.
     */
    public double getRelativeSolvationEnergy() {
        return relativeSolvationEnergy;
    }

    /**
     * {@inheritDoc}
     * <p>
     * The velocity array should only contain velocity data for active atoms.
     */
    @Override
    public void setVelocity(double[] velocity) {
        if (velocity == null) {
            return;
        }
        int index = 0;
        double[] vel = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive()) {
                vel[0] = velocity[index++];
                vel[1] = velocity[index++];
                vel[2] = velocity[index++];
                atoms[i].setVelocity(vel);
            }
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * The acceleration array should only contain acceleration data for active
     * atoms.
     */
    @Override
    public void setAcceleration(double[] acceleration) {
        if (acceleration == null) {
            return;
        }
        int index = 0;
        double[] accel = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive()) {
                accel[0] = acceleration[index++];
                accel[1] = acceleration[index++];
                accel[2] = acceleration[index++];
                atoms[i].setAcceleration(accel);
            }
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * The previousAcceleration array should only contain previous acceleration
     * data for active atoms.
     */
    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        if (previousAcceleration == null) {
            return;
        }
        int index = 0;
        double[] prev = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive()) {
                prev[0] = previousAcceleration[index++];
                prev[1] = previousAcceleration[index++];
                prev[2] = previousAcceleration[index++];
                atoms[i].setPreviousAcceleration(prev);
            }
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns an array of velocity values for active atoms.
     */
    @Override
    public double[] getVelocity(double[] velocity) {
        int n = getNumberOfVariables();
        if (velocity == null || velocity.length < n) {
            velocity = new double[n];
        }
        int index = 0;
        double[] v = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                a.getVelocity(v);
                velocity[index++] = v[0];
                velocity[index++] = v[1];
                velocity[index++] = v[2];
            }
        }
        return velocity;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns an array of acceleration values for active atoms.
     */
    @Override
    public double[] getAcceleration(double[] acceleration) {
        int n = getNumberOfVariables();
        if (acceleration == null || acceleration.length < n) {
            acceleration = new double[n];
        }
        int index = 0;
        double[] a = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive()) {
                atoms[i].getAcceleration(a);
                acceleration[index++] = a[0];
                acceleration[index++] = a[1];
                acceleration[index++] = a[2];
            }
        }
        return acceleration;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns an array of previous acceleration values for active atoms.
     */
    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        int n = getNumberOfVariables();
        if (previousAcceleration == null || previousAcceleration.length < n) {
            previousAcceleration = new double[n];
        }
        int index = 0;
        double[] a = new double[3];
        for (int i = 0; i < nAtoms; i++) {
            if (atoms[i].isActive()) {
                atoms[i].getPreviousAcceleration(a);
                previousAcceleration[index++] = a[0];
                previousAcceleration[index++] = a[1];
                previousAcceleration[index++] = a[2];
            }
        }
        return previousAcceleration;
    }

    /**
     * <p>Getter for the field <code>bonds</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.Bond} objects.
     */
    public Bond[] getBonds() {
        return bonds;
    }

    /**
     * <p>Getter for the field <code>angles</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.Angle} objects.
     */
    public Angle[] getAngles() {
        return angles;
    }

    /**
     * <p>Getter for the field <code>improperTorsions</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.ImproperTorsion} objects.
     */
    public ImproperTorsion[] getImproperTorsions() {
        return improperTorsions;
    }

    /**
     * <p>Getter for the field <code>ureyBradleys</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.UreyBradley} objects.
     */
    public UreyBradley[] getUreyBradleys() {
        return ureyBradleys;
    }

    /**
     * <p>Getter for the field <code>outOfPlaneBends</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.OutOfPlaneBend} objects.
     */
    public OutOfPlaneBend[] getOutOfPlaneBends() {
        return outOfPlaneBends;
    }

    /**
     * <p>Getter for the field <code>stretchBends</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.StretchBend} objects.
     */
    public StretchBend[] getStretchBends() {
        return stretchBends;
    }

    /**
     * <p>Getter for the field <code>torsions</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.Torsion} objects.
     */
    public Torsion[] getTorsions() {
        return torsions;
    }

    /**
     * <p>Getter for the field <code>piOrbitalTorsions</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.PiOrbitalTorsion} objects.
     */
    public PiOrbitalTorsion[] getPiOrbitalTorsions() {
        return piOrbitalTorsions;
    }

    /**
     * <p>Getter for the field <code>torsionTorsions</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.TorsionTorsion} objects.
     */
    public TorsionTorsion[] getTorsionTorsions() {
        return torsionTorsions;
    }

    /**
     * <p>Getter for the field <code>stretchTorsions</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.StretchTorsion} objects.
     */
    public StretchTorsion[] getStretchTorsions() {
        return stretchTorsions;
    }

    /**
     * <p>Getter for the field <code>angleTorsions</code>.</p>
     *
     * @return an array of {@link ffx.potential.bonded.AngleTorsion} objects.
     */
    public AngleTorsion[] getAngleTorsions() {
        return angleTorsions;
    }

    private class BondedRegion extends ParallelRegion {

        // Flag to indicate gradient computation.
        private boolean gradient = false;

        private AtomicDoubleArrayImpl atomicDoubleArrayImpl;
        private AtomicDoubleArray3D grad;
        private AtomicDoubleArray3D lambdaGrad;

        // Shared RMSD variables.
        private final SharedDouble sharedBondRMSD;
        private final SharedDouble sharedAngleRMSD;

        // Shared force field bonded energy accumulation variables.
        private final SharedDouble sharedBondEnergy;
        private final SharedDouble sharedAngleEnergy;
        private final SharedDouble sharedOutOfPlaneBendEnergy;
        private final SharedDouble sharedPiOrbitalTorsionEnergy;
        private final SharedDouble sharedStretchBendEnergy;
        private final SharedDouble sharedUreyBradleyEnergy;
        private final SharedDouble sharedImproperTorsionEnergy;
        private final SharedDouble sharedTorsionEnergy;
        private final SharedDouble sharedStretchTorsionEnergy;
        private final SharedDouble sharedAngleTorsionEnergy;
        private final SharedDouble sharedTorsionTorsionEnergy;
        // Shared restraint terms.
        private final SharedDouble sharedRestraintBondEnergy;

        // Number of threads.
        private final int nThreads;

        // Gradient loops.
        private final GradInitLoop[] gradInitLoops;
        private final GradReduceLoop[] gradReduceLoops;

        // Force field bonded energy parallel loops.
        private final BondedTermLoop[] bondLoops;
        private final BondedTermLoop[] angleLoops;
        private final BondedTermLoop[] outOfPlaneBendLoops;
        private final BondedTermLoop[] improperTorsionLoops;
        private final BondedTermLoop[] piOrbitalTorsionLoops;
        private final BondedTermLoop[] stretchBendLoops;
        private final BondedTermLoop[] torsionLoops;
        private final BondedTermLoop[] stretchTorsionLoops;
        private final BondedTermLoop[] angleTorsionLoops;
        private final BondedTermLoop[] torsionTorsionLoops;
        private final BondedTermLoop[] ureyBradleyLoops;

        // Retraint energy parallel loops.
        private final BondedTermLoop[] restraintBondLoops;

        BondedRegion() {

            // Allocate shared RMSD variables.
            sharedAngleRMSD = new SharedDouble();
            sharedBondRMSD = new SharedDouble();

            // Allocate shared force field bonded energy variables.
            sharedBondEnergy = new SharedDouble();
            sharedAngleEnergy = new SharedDouble();
            sharedImproperTorsionEnergy = new SharedDouble();
            sharedOutOfPlaneBendEnergy = new SharedDouble();
            sharedPiOrbitalTorsionEnergy = new SharedDouble();
            sharedStretchBendEnergy = new SharedDouble();
            sharedTorsionEnergy = new SharedDouble();
            sharedStretchTorsionEnergy = new SharedDouble();
            sharedAngleTorsionEnergy = new SharedDouble();
            sharedTorsionTorsionEnergy = new SharedDouble();
            sharedUreyBradleyEnergy = new SharedDouble();

            // Allocate shared restraint variables.
            sharedRestraintBondEnergy = new SharedDouble();

            nThreads = parallelTeam.getThreadCount();

            // Gradient initialization loops.
            gradInitLoops = new GradInitLoop[nThreads];
            gradReduceLoops = new GradReduceLoop[nThreads];

            // Allocate memory for force field bonded energy loop arrays.
            angleLoops = new BondedTermLoop[nThreads];
            bondLoops = new BondedTermLoop[nThreads];
            improperTorsionLoops = new BondedTermLoop[nThreads];
            outOfPlaneBendLoops = new BondedTermLoop[nThreads];
            piOrbitalTorsionLoops = new BondedTermLoop[nThreads];
            stretchBendLoops = new BondedTermLoop[nThreads];
            torsionLoops = new BondedTermLoop[nThreads];
            stretchTorsionLoops = new BondedTermLoop[nThreads];
            angleTorsionLoops = new BondedTermLoop[nThreads];
            torsionTorsionLoops = new BondedTermLoop[nThreads];
            ureyBradleyLoops = new BondedTermLoop[nThreads];

            // Allocate memory for restrain energy terms.
            restraintBondLoops = new BondedTermLoop[nThreads];

            // Define how the gradient will be accumulated.
            atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
            ForceField forceField = molecularAssembly.getForceField();
            String value = forceField.getString(ForceFieldString.ARRAY_REDUCTION, "MULTI");
            try {
                atomicDoubleArrayImpl = AtomicDoubleArrayImpl.valueOf(toEnumForm(value));
            } catch (Exception e) {
                logger.info(format(" Unrecognized ARRAY-REDUCTION %s; defaulting to %s", value, atomicDoubleArrayImpl));
            }

            grad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, nThreads);
            if (lambdaTerm) {
                lambdaGrad = new AtomicDoubleArray3D(atomicDoubleArrayImpl, nAtoms, nThreads);
            }
        }

        public void setGradient(boolean gradient) {
            this.gradient = gradient;
        }

        @Override
        public void start() {
            // Zero out shared RMSD values.
            sharedAngleRMSD.set(0.0);
            sharedBondRMSD.set(0.0);

            // Zero out shared energy values.
            sharedAngleEnergy.set(0.0);
            sharedBondEnergy.set(0.0);
            sharedImproperTorsionEnergy.set(0.0);
            sharedOutOfPlaneBendEnergy.set(0.0);
            sharedPiOrbitalTorsionEnergy.set(0.0);
            sharedStretchBendEnergy.set(0.0);
            sharedTorsionEnergy.set(0.0);
            sharedStretchTorsionEnergy.set(0.0);
            sharedAngleTorsionEnergy.set(0.0);
            sharedTorsionTorsionEnergy.set(0.0);
            sharedUreyBradleyEnergy.set(0.0);

            // Zero out shared restraint energy values.
            sharedRestraintBondEnergy.set(0.0);

            // Assure capacity of the gradient arrays.
            if (gradient) {
                grad.alloc(nAtoms);
            }
            if (lambdaTerm) {
                lambdaGrad.alloc(nAtoms);
            }
        }

        @Override
        public void finish() {
            // Finalize bond and angle RMSD values.
            if (bondTerm) {
                bondRMSD = sqrt(sharedBondRMSD.get() / bonds.length);
            }
            if (angleTerm) {
                angleRMSD = sqrt(sharedAngleRMSD.get() / angles.length);
            }

            // Load shared energy values into their respective fields.
            angleEnergy = sharedAngleEnergy.get();
            bondEnergy = sharedBondEnergy.get();
            improperTorsionEnergy = sharedImproperTorsionEnergy.get();
            outOfPlaneBendEnergy = sharedOutOfPlaneBendEnergy.get();
            piOrbitalTorsionEnergy = sharedPiOrbitalTorsionEnergy.get();
            stretchBendEnergy = sharedStretchBendEnergy.get();
            ureyBradleyEnergy = sharedUreyBradleyEnergy.get();
            torsionEnergy = sharedTorsionEnergy.get();
            stretchTorsionEnergy = sharedStretchTorsionEnergy.get();
            angleTorsionEnergy = sharedAngleTorsionEnergy.get();
            torsionTorsionEnergy = sharedTorsionTorsionEnergy.get();
            ureyBradleyEnergy = sharedUreyBradleyEnergy.get();

            // Load shared restraint energy values.
            restraintBondEnergy = sharedRestraintBondEnergy.get();

            if (esvTerm) {
                if (angleTerm) {
                    for (BondedTerm term : angles) {
                        term.reduceEsvDeriv();
                    }
                }
                if (bondTerm) {
                    for (BondedTerm term : bonds) {
                        term.reduceEsvDeriv();
                    }
                }
                if (improperTorsionTerm) {
                    for (BondedTerm term : improperTorsions) {
                        term.reduceEsvDeriv();
                    }
                }
                if (outOfPlaneBendTerm) {
                    for (BondedTerm term : outOfPlaneBends) {
                        term.reduceEsvDeriv();
                    }
                }
                if (piOrbitalTorsionTerm) {
                    for (BondedTerm term : piOrbitalTorsions) {
                        term.reduceEsvDeriv();
                    }
                }
                if (stretchBendTerm) {
                    for (BondedTerm term : stretchBends) {
                        term.reduceEsvDeriv();
                    }
                }
                if (torsionTerm) {
                    for (BondedTerm term : torsions) {
                        term.reduceEsvDeriv();
                    }
                }
                if (stretchTorsionTerm) {
                    for (BondedTerm term : stretchTorsions) {
                        term.reduceEsvDeriv();
                    }
                }
                if (angleTorsionTerm) {
                    for (BondedTerm term : angleTorsions) {
                        term.reduceEsvDeriv();
                    }
                }
                if (torsionTorsionTerm) {
                    for (BondedTerm term : torsionTorsions) {
                        term.reduceEsvDeriv();
                    }
                }
                if (ureyBradleyTerm) {
                    for (BondedTerm term : ureyBradleys) {
                        term.reduceEsvDeriv();
                    }
                }
                if (restraintBondTerm) {
                    for (BondedTerm term : restraintBonds) {
                        term.reduceEsvDeriv();
                    }
                }
            }
        }

        @Override
        public void run() throws Exception {
            int threadID = getThreadIndex();

            // Initialize the Gradient to zero.
            if (gradient) {
                if (gradInitLoops[threadID] == null) {
                    gradInitLoops[threadID] = new GradInitLoop();
                }
                execute(0, nAtoms - 1, gradInitLoops[threadID]);
            }

            // Evaluate force field bonded energy terms in parallel.
            if (angleTerm) {
                if (angleLoops[threadID] == null) {
                    angleLoops[threadID] = new BondedTermLoop(angles, sharedAngleEnergy, sharedAngleRMSD);
                }
                if (threadID == 0) {
                    angleTime = -System.nanoTime();
                }
                execute(0, nAngles - 1, angleLoops[threadID]);
                if (threadID == 0) {
                    angleTime += System.nanoTime();
                }
            }

            if (bondTerm) {
                if (bondLoops[threadID] == null) {
                    bondLoops[threadID] = new BondedTermLoop(bonds, sharedBondEnergy, sharedBondRMSD);
                }
                if (threadID == 0) {
                    bondTime = -System.nanoTime();
                }
                execute(0, nBonds - 1, bondLoops[threadID]);
                if (threadID == 0) {
                    bondTime += System.nanoTime();
                }
            }

            if (improperTorsionTerm) {
                if (improperTorsionLoops[threadID] == null) {
                    improperTorsionLoops[threadID] = new BondedTermLoop(improperTorsions, sharedImproperTorsionEnergy);
                }
                if (threadID == 0) {
                    improperTorsionTime = -System.nanoTime();
                }
                execute(0, nImproperTorsions - 1, improperTorsionLoops[threadID]);
                if (threadID == 0) {
                    improperTorsionTime += System.nanoTime();
                }
            }

            if (outOfPlaneBendTerm) {
                if (outOfPlaneBendLoops[threadID] == null) {
                    outOfPlaneBendLoops[threadID] = new BondedTermLoop(outOfPlaneBends, sharedOutOfPlaneBendEnergy);
                }
                if (threadID == 0) {
                    outOfPlaneBendTime = -System.nanoTime();
                }
                execute(0, nOutOfPlaneBends - 1, outOfPlaneBendLoops[threadID]);
                if (threadID == 0) {
                    outOfPlaneBendTime += System.nanoTime();
                }
            }

            if (piOrbitalTorsionTerm) {
                if (piOrbitalTorsionLoops[threadID] == null) {
                    piOrbitalTorsionLoops[threadID] = new BondedTermLoop(piOrbitalTorsions, sharedPiOrbitalTorsionEnergy);
                }
                if (threadID == 0) {
                    piOrbitalTorsionTime = -System.nanoTime();
                }
                execute(0, nPiOrbitalTorsions - 1, piOrbitalTorsionLoops[threadID]);
                if (threadID == 0) {
                    piOrbitalTorsionTime += System.nanoTime();
                }
            }

            if (stretchBendTerm) {
                if (stretchBendLoops[threadID] == null) {
                    stretchBendLoops[threadID] = new BondedTermLoop(stretchBends, sharedStretchBendEnergy);
                }
                if (threadID == 0) {
                    stretchBendTime = -System.nanoTime();
                }
                execute(0, nStretchBends - 1, stretchBendLoops[threadID]);
                if (threadID == 0) {
                    stretchBendTime += System.nanoTime();
                }
            }

            if (torsionTerm) {
                if (torsionLoops[threadID] == null) {
                    torsionLoops[threadID] = new BondedTermLoop(torsions, sharedTorsionEnergy);
                }
                if (threadID == 0) {
                    torsionTime = -System.nanoTime();
                }
                execute(0, nTorsions - 1, torsionLoops[threadID]);
                if (threadID == 0) {
                    torsionTime += System.nanoTime();
                }
            }

            if (stretchTorsionTerm) {
                if (stretchTorsionLoops[threadID] == null) {
                    stretchTorsionLoops[threadID] = new BondedTermLoop(stretchTorsions, sharedStretchTorsionEnergy);
                }
                if (threadID == 0) {
                    stretchTorsionTime = -System.nanoTime();
                }
                execute(0, nStretchTorsions - 1, stretchTorsionLoops[threadID]);
                if (threadID == 0) {
                    stretchTorsionTime += System.nanoTime();
                }
            }

            if (angleTorsionTerm) {
                if (angleTorsionLoops[threadID] == null) {
                    angleTorsionLoops[threadID] = new BondedTermLoop(angleTorsions, sharedAngleTorsionEnergy);
                }
                if (threadID == 0) {
                    angleTorsionTime = -System.nanoTime();
                }
                execute(0, nAngleTorsions - 1, angleTorsionLoops[threadID]);
                if (threadID == 0) {
                    angleTorsionTime += System.nanoTime();
                }
            }

            if (torsionTorsionTerm) {
                if (torsionTorsionLoops[threadID] == null) {
                    torsionTorsionLoops[threadID] = new BondedTermLoop(torsionTorsions, sharedTorsionTorsionEnergy);
                }
                if (threadID == 0) {
                    torsionTorsionTime = -System.nanoTime();
                }
                execute(0, nTorsionTorsions - 1, torsionTorsionLoops[threadID]);
                if (threadID == 0) {
                    torsionTorsionTime += System.nanoTime();
                }
            }

            if (ureyBradleyTerm) {
                if (ureyBradleyLoops[threadID] == null) {
                    ureyBradleyLoops[threadID] = new BondedTermLoop(ureyBradleys, sharedUreyBradleyEnergy);
                }
                if (threadID == 0) {
                    ureyBradleyTime = -System.nanoTime();
                }
                execute(0, nUreyBradleys - 1, ureyBradleyLoops[threadID]);
                if (threadID == 0) {
                    ureyBradleyTime += System.nanoTime();
                }
            }

            // Evaluate restraint terms in parallel.
            if (restraintBondTerm) {
                if (restraintBondLoops[threadID] == null) {
                    restraintBondLoops[threadID] = new BondedTermLoop(restraintBonds, sharedRestraintBondEnergy);
                }
                if (threadID == 0) {
                    restraintBondTime = -System.nanoTime();
                }
                execute(0, nRestraintBonds - 1, restraintBondLoops[threadID]);
                if (threadID == 0) {
                    restraintBondTime += System.nanoTime();
                }
            }

            // Reduce the Gradient and load it into Atom instances.
            if (gradient) {
                if (gradReduceLoops[threadID] == null) {
                    gradReduceLoops[threadID] = new GradReduceLoop();
                }
                execute(0, nAtoms - 1, gradReduceLoops[threadID]);
            }
        }

        private class GradInitLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int first, int last) throws Exception {
                int threadID = getThreadIndex();
                if (gradient) {
                    grad.reset(threadID, first, last);
                }
                if (lambdaTerm) {
                    lambdaGrad.reset(threadID, first, last);
                }
            }
        }

        private class GradReduceLoop extends IntegerForLoop {

            @Override
            public IntegerSchedule schedule() {
                return IntegerSchedule.fixed();
            }

            @Override
            public void run(int first, int last) throws Exception {
                if (gradient) {
                    grad.reduce(first, last);
                    for (int i = first; i <= last; i++) {
                        Atom a = atoms[i];
                        a.setXYZGradient(grad.getX(i), grad.getY(i), grad.getZ(i));
                    }
                }
                if (lambdaTerm) {
                    lambdaGrad.reduce(first, last);
                    for (int i = first; i <= last; i++) {
                        Atom a = atoms[i];
                        a.setLambdaXYZGradient(lambdaGrad.getX(i), lambdaGrad.getY(i), lambdaGrad.getZ(i));
                    }
                }
            }
        }

        private class BondedTermLoop extends IntegerForLoop {

            private final BondedTerm[] terms;
            private final SharedDouble sharedEnergy;
            private final SharedDouble sharedRMSD;
            private final boolean computeRMSD;
            private double localEnergy;
            private double localRMSD;
            private int threadID;

            BondedTermLoop(BondedTerm[] terms, SharedDouble sharedEnergy) {
                this(terms, sharedEnergy, null);
            }

            BondedTermLoop(BondedTerm[] terms, SharedDouble sharedEnergy, SharedDouble sharedRMSD) {
                this.terms = terms;
                this.sharedEnergy = sharedEnergy;
                this.sharedRMSD = sharedRMSD;
                computeRMSD = (sharedRMSD != null);
            }

            @Override
            public void start() {
                localEnergy = 0.0;
                localRMSD = 0.0;
                threadID = getThreadIndex();
            }

            @Override
            public void finish() {
                sharedEnergy.addAndGet(localEnergy);
                if (computeRMSD) {
                    sharedRMSD.addAndGet(localRMSD);
                }
            }

            @Override
            public void run(int first, int last) throws Exception {
                for (int i = first; i <= last; i++) {
                    BondedTerm term = terms[i];
                    if (!lambdaBondedTerms || term.applyLambda()) {
                        localEnergy += term.energy(gradient, threadID, grad, lambdaGrad);
                        if (computeRMSD) {
                            double value = term.getValue();
                            localRMSD += value * value;
                        }
                    }
                }
            }
        }
    }

    /**
     * Platform describes a set of force field implementations that include a
     * pure Java reference implementation (FFX), and two OpenMM implementations
     * (OMM_CUDA and OMM_REF are supported)
     * <p>
     * FFX: reference FFX implementation
     * <p>
     * OMM: Currently an alias for OMM_CUDA, may eventually become "try to find
     * best OpenMM implementation" OMM_CUDA:
     * <p>
     * OpenMM CUDA implementation OMM_REF: OpenMM reference implementation
     * <p>
     * OMM_OPTCPU: Optimized OpenMM CPU implementation (no AMOEBA)
     * <p>
     * OMM_OPENCL: OpenMM OpenCL implementation (no AMOEBA)
     */
    public enum Platform {
        FFX, OMM, OMM_CUDA, OMM_REF, OMM_OPTCPU, OMM_OPENCL;
    }
}
