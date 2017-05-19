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
package ffx.potential;

import java.io.File;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;
import static java.util.Arrays.fill;

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
import ffx.crystal.ReplicatesCrystal;
import ffx.numerics.AdderDoubleArray;
import ffx.numerics.AtomicDoubleArray;
import ffx.numerics.AtomicDoubleArray.AtomicDoubleArrayImpl;
import ffx.numerics.MultiDoubleArray;
import ffx.numerics.PJDoubleArray;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Atom.Resolution;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.BondedTerm;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.MultiResidue;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.ROLS;
import ffx.potential.bonded.RelativeSolvation;
import ffx.potential.bonded.RelativeSolvation.SolvationLibrary;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.bonded.StretchBend;
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

import static ffx.numerics.AtomicDoubleArray.AtomicDoubleArrayImpl.MULTI;
import static ffx.potential.PotentialComponent.*;
import static ffx.potential.bonded.Atom.Resolution.AMOEBA;
import static ffx.potential.bonded.Atom.Resolution.FIXEDCHARGE;
import static ffx.potential.extended.ExtUtils.prop;
import static ffx.potential.parameters.ForceField.ForceFieldString.ARRAY_REDUCTION;
import static ffx.potential.parameters.ForceField.toEnumForm;

/**
 * Compute the potential energy and derivatives of an AMOEBA system.
 *
 * @author Michael J. Schnieders
 *
 * @since 1.0
 */
public class ForceFieldEnergy implements CrystalPotential, LambdaInterface {

    private static final Logger logger = Logger.getLogger(ForceFieldEnergy.class.getName());

    private static final double toSeconds = 0.000000001;

    protected final MolecularAssembly molecularAssembly;
    private Atom[] atoms;
    /**
     * Contains ALL atoms, both foreground and background. Background atoms need
     * to be present to be included in bonded terms
     */
    private Atom[] atomsExtended;
    private Crystal crystal;
    private final ParallelTeam parallelTeam;
    private BondedRegion bondedRegion;
    private STATE state = STATE.BOTH;
    private Bond bonds[];
    private Angle angles[];
    private StretchBend stretchBends[];
    private UreyBradley ureyBradleys[];
    private OutOfPlaneBend outOfPlaneBends[];
    private Torsion torsions[];
    private PiOrbitalTorsion piOrbitalTorsions[];
    private TorsionTorsion torsionTorsions[];
    private ImproperTorsion improperTorsions[];
    private RestraintBond restraintBonds[];
    private RelativeSolvation relativeSolvation;
    private final VanDerWaals vanderWaals;
    private ParticleMeshEwald particleMeshEwald;
    private final NCSRestraint ncsRestraint;
    private final List<CoordRestraint> coordRestraints;
    private final CoordRestraint autoCoordRestraint;
    private final COMRestraint comRestraint;
    private int nAtoms;
    private int nBonds;
    private int nAngles;
    private int nStretchBends;
    private int nUreyBradleys;
    private int nOutOfPlaneBends;
    private int nTorsions;
    private int nImproperTorsions;
    private int nPiOrbitalTorsions;
    private int nTorsionTorsions;
    private int nRestraintBonds;
    private int nVanDerWaalInteractions;
    private int nPermanentInteractions;
    private int nGKInteractions;
    private int nRelativeSolvations;
    private boolean bondTerm;
    private boolean angleTerm;
    private boolean stretchBendTerm;
    private boolean ureyBradleyTerm;
    private boolean outOfPlaneBendTerm;
    private boolean torsionTerm;
    private boolean improperTorsionTerm;
    private boolean piOrbitalTorsionTerm;
    private boolean torsionTorsionTerm;
    private boolean restraintBondTerm;
    private boolean vanderWaalsTerm;
    private boolean multipoleTerm;
    private boolean polarizationTerm;
    private boolean generalizedKirkwoodTerm;
    private boolean ncsTerm;
    private boolean restrainTerm;
    private boolean restrainWithLambda;
    private boolean comTerm;
    private boolean esvTerm;
    private boolean lambdaTerm;
    protected boolean lambdaBondedTerms = false;
    private boolean relativeSolvationTerm;
    private boolean rigidHydrogens = false;
    private boolean lambdaTorsions = false;
    private double rigidScale = 1.0;
    private boolean bondTermOrig;
    private boolean angleTermOrig;
    private boolean stretchBendTermOrig;
    private boolean ureyBradleyTermOrig;
    private boolean outOfPlaneBendTermOrig;
    private boolean torsionTermOrig;
    private boolean improperTorsionTermOrig;
    private boolean piOrbitalTorsionTermOrig;
    private boolean torsionTorsionTermOrig;
    private boolean restraintBondTermOrig;
    private boolean vanderWaalsTermOrig;
    private boolean multipoleTermOrig;
    private boolean polarizationTermOrig;
    private boolean generalizedKirkwoodTermOrig;
    private boolean ncsTermOrig;
    private boolean restrainTermOrig;
    private boolean comTermOrig;
    private double bondEnergy, bondRMSD;
    private double angleEnergy, angleRMSD;
    private double stretchBendEnergy;
    private double ureyBradleyEnergy;
    private double outOfPlaneBendEnergy;
    private double torsionEnergy;
    private double improperTorsionEnergy;
    private double piOrbitalTorsionEnergy;
    private double torsionTorsionEnergy;
    private double restraintBondEnergy;
    private double totalBondedEnergy;
    private double vanDerWaalsEnergy;
    private double permanentMultipoleEnergy;
    private double permanentRealSpaceEnergy;
    private double permanentSelfEnergy;
    private double permanentReciprocalEnergy;
    private double polarizationEnergy;
    private double inducedRealSpaceEnergy;
    private double inducedSelfEnergy;
    private double inducedReciprocalEnergy;
    private double totalMultipoleEnergy;
    private double totalNonBondedEnergy;
    private double solvationEnergy;
    private double relativeSolvationEnergy;
    private double ncsEnergy;
    private double restrainEnergy;
    private double comRestraintEnergy;
    private double esvBias;
    private double totalEnergy;
    private long bondTime, angleTime, stretchBendTime, ureyBradleyTime;
    private long outOfPlaneBendTime, torsionTime, piOrbitalTorsionTime, improperTorsionTime;
    private long torsionTorsionTime, vanDerWaalsTime, electrostaticTime;
    private long restraintBondTime, ncsTime, coordRestraintTime, comRestraintTime;
    private long totalTime;
    private double lambda = 1.0;
    protected double[] optimizationScaling = null;
    private VARIABLE_TYPE[] variableTypes = null;
    private double xyz[] = null;
    private boolean printOnFailure;
    private boolean printCompact = prop("ffe.combineBonded", false);
    private boolean printOverride = prop("ffe.printOverride", false);
    private final boolean noHeader = prop("ffe.noHeader", false);
    private final boolean decomposePme = prop("pme.decompose", false);
    /**
     * *************************************
     */
    /*      Extended System Variables       */
    private ExtendedSystem esvSystem = null;
    private final boolean pmeQI = prop("pme.qi", false);
    /**
     * *************************************
     */

    private Resolution resolution = Resolution.AMOEBA;

    /**
     * <p>
     * Constructor for ForceFieldEnergy.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     */
    public ForceFieldEnergy(MolecularAssembly molecularAssembly) {
        this(molecularAssembly, null);
    }

    /**
     * <p>
     * Constructor for ForceFieldEnergy.</p>
     *
     * @param molecularAssembly a {@link ffx.potential.MolecularAssembly}
     * object.
     * @param restraints list of {@link ffx.potential.nonbonded.CoordRestraint}
     * objects.
     */
    public ForceFieldEnergy(MolecularAssembly molecularAssembly, List<CoordRestraint> restraints) {
        this(molecularAssembly, restraints, ParallelTeam.getDefaultThreadCount());
    }

    public ForceFieldEnergy(MolecularAssembly molecularAssembly, List<CoordRestraint> restraints, int numThreads) {
        if (noHeader) {
            logger.setLevel(Level.WARNING);
        }

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
        lambdaTerm = forceField.getBoolean(ForceField.ForceFieldBoolean.LAMBDATERM, false);
        restrainTerm = forceField.getBoolean(ForceFieldBoolean.RESTRAINTERM, false);
        restrainWithLambda = forceField.getBoolean(ForceField.ForceFieldBoolean.RESTRAIN_WITH_LAMBDA, true);
        comTerm = forceField.getBoolean(ForceFieldBoolean.COMRESTRAINTERM, false);
        lambdaTorsions = forceField.getBoolean(ForceFieldBoolean.LAMBDA_TORSIONS, false);
        printOnFailure = forceField.getBoolean(ForceFieldBoolean.PRINT_ON_FAILURE, false);

        // For RESPA
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
        comTermOrig = comTerm;

        // Define the cutoff lengths.
        double vdwOff = forceField.getDouble(ForceFieldDouble.VDW_CUTOFF, 9.0);
        double ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
        if (ewaldOff > vdwOff) {
            vdwOff = ewaldOff;
        }
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
            a = 2.0 * maxr;
            b = 2.0 * maxr;
            c = 2.0 * maxr;
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
            ArrayList<ROLS> bond = molecularAssembly.getBondList();
            nBonds = bond.size();
            bonds = bond.toArray(new Bond[nBonds]);
            Arrays.sort(bonds);
            if (nBonds > 0) {
                logger.info(format("  Bonds:                             %10d", nBonds));
            }
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
            if (nAngles > 0) {
                logger.info(format("  Angles:                            %10d", nAngles));
            }
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
            if (nStretchBends > 0) {
                logger.info(format("  Stretch-Bends:                     %10d", nStretchBends));
            }
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
            if (nUreyBradleys > 0) {
                logger.info(format("  Urey-Bradleys:                     %10d", nUreyBradleys));
            }
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
            ArrayList<ROLS> torsion = molecularAssembly.getTorsionList();
            nTorsions = torsion.size();
            torsions = torsion.toArray(new Torsion[nTorsions]);
            if (nTorsions > 0 && torsionScale == 1.0) {
                logger.info(format("  Torsions:                          %10d", nTorsions));
            } else {
                logger.info(format("  Torsions (%5.2f):                  %10d", torsionScale, nTorsions));
            }
        } else {
            nTorsions = 0;
            torsions = null;
        }

        // Collect, count, pack and sort pi-orbital torsions.
        if (piOrbitalTorsionTerm) {
            ArrayList<ROLS> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = piOrbitalTorsion.size();
            piOrbitalTorsions = piOrbitalTorsion.toArray(new PiOrbitalTorsion[nPiOrbitalTorsions]);
            if (nPiOrbitalTorsions > 0 && torsionScale == 1.0) {
                logger.info(format("  Pi-Orbital Torsions:               %10d", nPiOrbitalTorsions));
            } else {
                logger.info(format("  Pi-Orbital Torsions (%5.2f):       %10d", torsionScale, nPiOrbitalTorsions));
            }
        } else {
            nPiOrbitalTorsions = 0;
            piOrbitalTorsions = null;
        }

        // Collect, count, pack and sort torsion-torsions.
        if (torsionTorsionTerm) {
            ArrayList<ROLS> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = torsionTorsion.size();
            torsionTorsions = torsionTorsion.toArray(new TorsionTorsion[nTorsionTorsions]);
            if (nTorsionTorsions > 0) {
                logger.info(format("  Torsion-Torsions:                  %10d", nTorsionTorsions));
            }
        } else {
            nTorsionTorsions = 0;
            torsionTorsions = null;
        }

        // Collect, count, pack and sort improper torsions.
        if (improperTorsionTerm) {
            ArrayList<ROLS> improperTorsion = molecularAssembly.getImproperTorsionList();
            nImproperTorsions = improperTorsion.size();
            improperTorsions = improperTorsion.toArray(new ImproperTorsion[nImproperTorsions]);
            if (nImproperTorsions > 0) {
                logger.info(format("  Improper Torsions:                 %10d", nImproperTorsions));
            }

        } else {
            nImproperTorsions = 0;
            improperTorsions = null;
        }

        logger.info("\n Non-Bonded Terms");

        int molecule[] = molecularAssembly.getMoleculeNumbers();
        if (vanderWaalsTerm) {
            vanderWaals = new VanDerWaals(atoms, molecule, crystal, forceField, parallelTeam);
        } else {
            vanderWaals = null;
        }

        if (multipoleTerm) {
            ELEC_FORM form;
            if (name.contains("OPLS") || name.contains("AMBER")) {
                form = ELEC_FORM.FIXED_CHARGE;
            } else {
                form = ELEC_FORM.PAM;
            }
            if (pmeQI) {
                particleMeshEwald = new ParticleMeshEwaldQI(atoms, molecule, forceField, crystal,
                        vanderWaals.getNeighborList(), form, parallelTeam);
            } else {
                particleMeshEwald = new ParticleMeshEwaldCart(atoms, molecule, forceField, crystal,
                        vanderWaals.getNeighborList(), form, parallelTeam);
            }
            double charge = molecularAssembly.getCharge(checkAllNodeCharges);
            logger.info(String.format(" Overall system charge:            %10.3f", charge));
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
            Polymer polymers[] = molecularAssembly.getChains();
            List<Molecule> molecules = molecularAssembly.getMolecules();
            List<MSNode> waters = molecularAssembly.getWaters();
            List<MSNode> ions = molecularAssembly.getIons();
            comRestraint = new COMRestraint(atoms, polymers, molecules, waters, ions, forceField);
        } else {
            comRestraint = null;
        }

        bondedRegion = new BondedRegion();

        molecularAssembly.setPotential(this);
        if (noHeader) {
            logger.setLevel(logger.getParent().getLevel());
        }
    }

    /**
     * Static factory method to create a ForceFieldEnergy, possibly via FFX or
     * OpenMM implementations.
     *
     * @param assembly To create FFE over
     * @return
     */
    public static ForceFieldEnergy energyFactory(MolecularAssembly assembly) {
        return energyFactory(assembly, null);
    }

    /**
     * Static factory method to create a ForceFieldEnergy, possibly via FFX or
     * OpenMM implementations.
     *
     * @param assembly To create FFE over
     * @param restraints Harmonic restraints
     * @return
     */
    public static ForceFieldEnergy energyFactory(MolecularAssembly assembly, List<CoordRestraint> restraints) {
        return energyFactory(assembly, restraints, ParallelTeam.getDefaultThreadCount());
    }

    /**
     * Static factory method to create a ForceFieldEnergy, possibly via FFX or
     * OpenMM implementations.
     *
     * @param assembly To create FFE over
     * @param restraints Harmonic restraints
     * @param numThreads Number of threads to use for FFX energy
     * @return
     */
    public static ForceFieldEnergy energyFactory(MolecularAssembly assembly, List<CoordRestraint> restraints, int numThreads) {

        ForceField ffield = assembly.getForceField();
        String eImString = ffield.getString(ForceFieldString.ENERGY_IMPLEMENTATION, "FFX").toUpperCase().replaceAll("-", "_");

        try {
            EnergyImplementation eImpl = EnergyImplementation.valueOf(eImString);
            switch (eImpl) {
                case FFX:
                    ForceFieldEnergy ffxEnergy = new ForceFieldEnergy(assembly, restraints, numThreads);
                    return ffxEnergy;
                case OMM:
                case OMM_REF: // Should be split from the code once we figure out how to specify a kernel.
                case OMM_CUDA:
                    try {
                        OpenMMForceFieldEnergy oEnergy = new OpenMMForceFieldEnergy(assembly);
                        return oEnergy;
                    } catch (Exception ex) {
                        logger.warning(String.format(" Exception in creating OpenMM wrapper over force field energy: %s", ex));
                        ffxEnergy = new ForceFieldEnergy(assembly, restraints, numThreads);
                        return ffxEnergy;
                    }
                case OMM_OPENCL:
                case OMM_OPTCPU:
                default:
                    logger.warning(String.format(" Energy implementation type %s not actually implemented at this time", eImpl));
                    ffxEnergy = new ForceFieldEnergy(assembly, restraints, numThreads);
                    return ffxEnergy;
            }
        } catch (IllegalArgumentException | NullPointerException ex) {
            logger.warning(String.format(" String %s did not match a known energy implementation", eImString));
            ForceFieldEnergy ffxEnergy = new ForceFieldEnergy(assembly, restraints, numThreads);
            return ffxEnergy;
        }
    }

    /**
     * Overwrites current esvSystem if present. Multiple ExtendedSystems is
     * possible but unnecessary; add all ESVs to one system (per FFE, at least).
     *
     * @param system
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
            }
            vanderWaals.attachExtendedSystem(system);
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
     * Need to remove degree of freedom that are lost to prevent heating.
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

        /* TODO Decide on only growing vs. always modifying xyz.
        if (xyz.length < 3 * nAtoms) {
            xyz = new double[nAtoms * 3];
        }   */
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
            ArrayList<ROLS> bond = molecularAssembly.getBondList();
            nBonds = 0;
            for (ROLS r : bond) {
                if (keep((Bond) r)) {
                    nBonds++;
                }
            }
            if (nBonds > bonds.length) {
                bonds = new Bond[nBonds];
            }
            Arrays.fill(bonds, null);
            nBonds = 0;
            for (ROLS r : bond) {
                if (keep((Bond) r)) {
                    bonds[nBonds++] = (Bond) r;
                }
            }
            Arrays.sort(bonds, 0, nBonds);
            if (nBonds > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Bonds:                             %10d", nBonds));
            }
        } else {
            nBonds = 0;
            bonds = null;
        }

        // Collect, count, pack and sort angles.
        if (angleTerm) {
            ArrayList<ROLS> angle = molecularAssembly.getAngleList();
            nAngles = 0;
            for (ROLS r : angle) {
                if (keep((Angle) r)) {
                    nAngles++;
                }
            }
            if (nAngles > angles.length) {
                angles = new Angle[nAngles];
            }
            Arrays.fill(angles, null);
            nAngles = 0;
            for (ROLS r : angle) {
                if (keep((Angle) r)) {
                    angles[nAngles++] = (Angle) r;
                }
            }

            Arrays.sort(angles, 0, nAngles);
            if (nAngles > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Angles:                            %10d", nAngles));
            }
        } else {
            nAngles = 0;
            angles = null;
        }

        // Collect, count, pack and sort stretch-bends.
        if (stretchBendTerm) {
            ArrayList<ROLS> stretchBend = molecularAssembly.getStretchBendList();
            nStretchBends = 0;
            for (ROLS r : stretchBend) {
                if (keep((StretchBend) r)) {
                    nStretchBends++;
                }
            }
            if (nStretchBends > stretchBends.length) {
                stretchBends = new StretchBend[nStretchBends];
            }
            Arrays.fill(stretchBends, null);
            nStretchBends = 0;
            for (ROLS r : stretchBend) {
                if (keep((StretchBend) r)) {
                    stretchBends[nStretchBends++] = (StretchBend) r;
                }
            }
            Arrays.sort(stretchBends, 0, nStretchBends);
            if (nStretchBends > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Stretch-Bends:                     %10d", nStretchBends));
            }
        } else {
            nStretchBends = 0;
            stretchBends = null;
        }

        // Collect, count, pack and sort Urey-Bradleys.
        if (ureyBradleyTerm) {
            ArrayList<ROLS> ureyBradley = molecularAssembly.getUreyBradleyList();
            nUreyBradleys = 0;
            for (ROLS r : ureyBradley) {
                if (keep((UreyBradley) r)) {
                    nUreyBradleys++;
                }
            }
            if (nUreyBradleys > ureyBradleys.length) {
                ureyBradleys = new UreyBradley[nUreyBradleys];
            }
            fill(ureyBradleys, null);
            nUreyBradleys = 0;
            for (ROLS r : ureyBradley) {
                if (keep((UreyBradley) r)) {
                    ureyBradleys[nUreyBradleys++] = (UreyBradley) r;
                }
            }
            Arrays.sort(ureyBradleys, 0, nUreyBradleys);
            if (nUreyBradleys > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Urey-Bradleys:                     %10d", nUreyBradleys));
            }
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
            nOutOfPlaneBends = 0;
            for (ROLS r : outOfPlaneBend) {
                if (keep((OutOfPlaneBend) r)) {
                    nOutOfPlaneBends++;
                }
            }
            if (nOutOfPlaneBends > outOfPlaneBends.length) {
                outOfPlaneBends = new OutOfPlaneBend[nOutOfPlaneBends];
            }
            fill(outOfPlaneBends, null);
            nOutOfPlaneBends = 0;
            for (ROLS r : outOfPlaneBend) {
                if (keep((OutOfPlaneBend) r)) {
                    outOfPlaneBends[nOutOfPlaneBends++] = (OutOfPlaneBend) r;
                }
            }
            Arrays.sort(outOfPlaneBends, 0, nOutOfPlaneBends);
            if (nOutOfPlaneBends > 0 && logger.isLoggable(Level.FINEST)) {
                logger.finest(format("  Out-of-Plane Bends:                %10d", nOutOfPlaneBends));
            }
        } else {
            nOutOfPlaneBends = 0;
            outOfPlaneBends = null;
        }

        // Collect, count, pack and sort torsions.
        if (torsionTerm) {
            ArrayList<ROLS> torsion = molecularAssembly.getTorsionList();
            nTorsions = 0;
            for (ROLS r : torsion) {
                if (keep((Torsion) r)) {
                    nTorsions++;
                }
            }
            if (nTorsions >= torsions.length) {
                torsions = new Torsion[nTorsions];
            }
            fill(torsions, null);
            nTorsions = 0;
            for (ROLS r : torsion) {
                if (keep((Torsion) r)) {
                    torsions[nTorsions++] = (Torsion) r;
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
            ArrayList<ROLS> piOrbitalTorsion = molecularAssembly.getPiOrbitalTorsionList();
            nPiOrbitalTorsions = 0;
            for (ROLS r : piOrbitalTorsion) {
                if (keep((PiOrbitalTorsion) r)) {
                    nPiOrbitalTorsions++;
                }
            }
            if (nPiOrbitalTorsions >= piOrbitalTorsions.length) {
                piOrbitalTorsions = new PiOrbitalTorsion[nPiOrbitalTorsions];
            }
            fill(piOrbitalTorsions, null);
            nPiOrbitalTorsions = 0;
            for (ROLS r : piOrbitalTorsion) {
                if (keep((PiOrbitalTorsion) r)) {
                    piOrbitalTorsions[nPiOrbitalTorsions++] = (PiOrbitalTorsion) r;
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
            ArrayList<ROLS> torsionTorsion = molecularAssembly.getTorsionTorsionList();
            nTorsionTorsions = 0;
            for (ROLS r : torsionTorsion) {
                if (keep((TorsionTorsion) r)) {
                    nTorsionTorsions++;
                }
            }
            if (nTorsionTorsions >= torsionTorsions.length) {
                torsionTorsions = new TorsionTorsion[nTorsionTorsions];
            }
            fill(torsionTorsions, null);
            nTorsionTorsions = 0;
            for (ROLS r : torsionTorsion) {
                if (keep((TorsionTorsion) r)) {
                    torsionTorsions[nTorsionTorsions++] = (TorsionTorsion) r;
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
            ArrayList<ROLS> improperTorsion = molecularAssembly.getImproperTorsionList();
            nImproperTorsions = 0;
            for (ROLS r : improperTorsion) {
                if (keep((ImproperTorsion) r)) {
                    nImproperTorsions++;
                }
            }
            if (nImproperTorsions >= improperTorsions.length) {
                improperTorsions = new ImproperTorsion[nImproperTorsions];
            }
            fill(improperTorsions, null);
            nImproperTorsions = 0;
            for (ROLS r : improperTorsion) {
                if (keep((ImproperTorsion) r)) {
                    improperTorsions[nImproperTorsions++] = (ImproperTorsion) r;
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

    public void setFixedCharges(Atom atoms[]) {
        if (particleMeshEwald != null) {
            particleMeshEwald.setFixedCharges(atoms);
        }
    }

    public void setLambdaBondedTerms(boolean lambdaBondedTerms) {
        this.lambdaBondedTerms = lambdaBondedTerms;
    }

    /**
     * @see ForceFieldEnergy::energy(boolean, boolean)
     */
    public double energy() {
        return energy(false, false);
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
        try {
            bondTime = 0;
            angleTime = 0;
            stretchBendTime = 0;
            ureyBradleyTime = 0;
            outOfPlaneBendTime = 0;
            torsionTime = 0;
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
            permanentSelfEnergy = 0.0;
            permanentReciprocalEnergy = 0.0;
            polarizationEnergy = 0.0;
            inducedRealSpaceEnergy = 0.0;
            inducedSelfEnergy = 0.0;
            inducedReciprocalEnergy = 0.0;
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

            /**
             * Computed the bonded energy terms in parallel.
             */
            try {
                bondedRegion.setGradient(gradient);
                parallelTeam.execute(bondedRegion);
            } catch (RuntimeException ex) {
                logger.warning("Runtime exception during bonded term calculation.");
                throw ex;
            } catch (Exception ex) {
                Utilities.printStackTrace(ex);
                logger.severe(ex.toString());
            }

            if (!lambdaBondedTerms) {
                /**
                 * Compute restraint terms.
                 */
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
                /**
                 * Compute non-bonded terms.
                 */
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
                    permanentSelfEnergy = particleMeshEwald.getPermSelfEnergy();
                    permanentReciprocalEnergy = particleMeshEwald.getPermRecipEnergy();

                    polarizationEnergy = particleMeshEwald.getPolarizationEnergy();
                    inducedRealSpaceEnergy = particleMeshEwald.getIndRealEnergy();
                    inducedSelfEnergy = particleMeshEwald.getIndSelfEnergy();
                    inducedReciprocalEnergy = particleMeshEwald.getIndRecipEnergy();
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
                            /**
                             * Reasonably confident that it should be -=, as we
                             * are trying to penalize residues with strong
                             * solvation energy.
                             */
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
                    + torsionEnergy + piOrbitalTorsionEnergy + improperTorsionEnergy
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

                String filename = String.format("%s-ERROR-%s.pdb",
                        FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                        timeString);

                PotentialsFunctions ef = new PotentialsUtils();
                filename = ef.versionFile(filename);
                logger.info(String.format(" Writing on-error snapshot to file %s", filename));
                ef.saveAsPDB(molecularAssembly, new File(filename));
                molecularAssembly.setFile(origFile);
            }

            if (ex.doCauseSevere()) {
                logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
                return 0.0;
            } else {
                throw ex; // Rethrow exception
            }
        }

        if (print || printOverride) {
            if (printCompact) {
                logger.info(this.toString());
            } else {
                StringBuilder sb = new StringBuilder();
                if (gradient) {
                    sb.append("\n Computed Potential Energy and Atomic Coordinate Gradients\n");
                } else {
                    sb.append("\n Computed Potential Energy\n");
                }
                sb.append(this);
                logger.info(sb.toString());
            }
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
     * Return the non-bonded components of energy (vdW, electrostatics,
     * solvation).
     *
     * @return Nonbonded energy
     */
    public double getNonbondedEnergy() {
        return getNonbondedEnergy(true);
    }

    /**
     * Return the non-bonded components of energy (vdW, electrostatics).
     *
     * @param includeSolv Include solvation energy
     * @return Nonbonded energy
     */
    public double getNonbondedEnergy(boolean includeSolv) {
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
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "BOND STRETCHING            : ", bondEnergy, nBonds));
            sb.append(String.format("REMARK   3   %s %g\n",
                    "BOND RMSD                  : ", bondRMSD));
        }
        if (angleTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "ANGLE BENDING              : ", angleEnergy, nAngles));
            sb.append(String.format("REMARK   3   %s %g\n",
                    "ANGLE RMSD                 : ", angleRMSD));
        }
        if (stretchBendTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "STRETCH-BEND               : ", stretchBendEnergy, nStretchBends));
        }
        if (ureyBradleyTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "UREY-BRADLEY               : ", ureyBradleyEnergy, nUreyBradleys));
        }
        if (outOfPlaneBendTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "OUT-OF-PLANE BEND          : ", outOfPlaneBendEnergy, nOutOfPlaneBends));
        }
        if (torsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "TORSIONAL ANGLE            : ", torsionEnergy, nTorsions));
        }
        if (piOrbitalTorsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "PI-ORBITAL TORSION         : ", piOrbitalTorsionEnergy, nPiOrbitalTorsions));
        }
        if (torsionTorsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "TORSION-TORSION            : ", torsionTorsionEnergy, nTorsionTorsions));
        }
        if (improperTorsionTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "IMPROPER TORSION           : ", improperTorsionEnergy, nImproperTorsions));
        }
        if (restraintBondTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "RESTRAINT BOND STRETCHING            : ", restraintBondEnergy, nRestraintBonds));
        }

        if (ncsTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "NCS RESTRAINT              : ", ncsEnergy, nAtoms));
        }

        if (restrainTerm && !coordRestraints.isEmpty()) {
            int nRests = 0;
            for (CoordRestraint restraint : coordRestraints) {
                nRests += restraint.getNumAtoms();
            }
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "COORDINATE RESTRAINTS      : ", restrainEnergy, nRests));
        }

        if (comTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "COM RESTRAINT              : ", comRestraintEnergy, nAtoms));
        }

        if (vanderWaalsTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "VAN DER WAALS              : ", vanDerWaalsEnergy, nVanDerWaalInteractions));
        }
        if (multipoleTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "ATOMIC MULTIPOLES          : ", permanentMultipoleEnergy, nPermanentInteractions));
        }
        if (polarizationTerm) {
            sb.append(String.format("REMARK   3   %s %g (%d)\n",
                    "POLARIZATION               : ", polarizationEnergy, nPermanentInteractions));
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
        StringBuilder sb = new StringBuilder();
        if (printCompact) {
            double totalBondedEnergy = bondEnergy + angleEnergy + stretchBendEnergy
                    + ureyBradleyEnergy + outOfPlaneBendEnergy + torsionEnergy
                    + piOrbitalTorsionEnergy + torsionTorsionEnergy + improperTorsionEnergy;
            int totalBondedInteractions = nBonds + nAngles + nStretchBends
                    + nUreyBradleys + nOutOfPlaneBends + nTorsions
                    + nPiOrbitalTorsions + nTorsionTorsions + nImproperTorsions;
            double totalBondedTime = (bondTime + angleTime) * toSeconds;
            sb.append(String.format("  %s %16.8f %12d %12.3f (%6.4f, %6.4f)\n",
                    "Bonded Terms      ", totalBondedEnergy, totalBondedInteractions,
                    totalBondedTime, bondRMSD, angleRMSD));
        } else {
            if (bondTerm && nBonds > 0) {
                sb.append(String.format("  %s %16.8f %12d %12.3f (%8.5f)\n",
                        "Bond Stretching   ", bondEnergy, nBonds,
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
            if (improperTorsionTerm && nImproperTorsions > 0) {
                sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                        "Improper Torsion  ", improperTorsionEnergy,
                        nImproperTorsions, improperTorsionTime * toSeconds));
            }
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
        if (restrainTerm && !coordRestraints.isEmpty()) {
            int nRests = 0;
            for (CoordRestraint restraint : coordRestraints) {
                nRests += restraint.getNumAtoms();
            }
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Coord. Restraints ", restrainEnergy, nRests,
                    coordRestraintTime * toSeconds));
        }
        if (comTerm) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "COM Restraint     ", comRestraintEnergy, nAtoms,
                    comRestraintTime * toSeconds));
        }
        if (vanderWaalsTerm && nVanDerWaalInteractions > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Van der Waals     ", vanDerWaalsEnergy,
                    nVanDerWaalInteractions, vanDerWaalsTime * toSeconds));
        }
        if (multipoleTerm && nPermanentInteractions > 0) {
            String pmeTitle = (particleMeshEwald instanceof ParticleMeshEwaldQI)
                    ? "Q.Int. Multipoles "
                    : "Atomic Multipoles ";
            if (polarizationTerm) {
                sb.append(String.format("  %s %16.8f %12d\n",
                        pmeTitle, permanentMultipoleEnergy, nPermanentInteractions));
            } else {
                sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                        pmeTitle, permanentMultipoleEnergy, nPermanentInteractions, electrostaticTime * toSeconds));
            }
        }
        if (polarizationTerm && nPermanentInteractions > 0) {
            sb.append(String.format("  %s %16.8f %12d %12.3f\n",
                    "Polarization      ", polarizationEnergy,
                    nPermanentInteractions, electrostaticTime * toSeconds));
        }
        if (generalizedKirkwoodTerm && nGKInteractions > 0) {
            sb.append(String.format("  %s %16.8f %12d\n",
                    "Solvation         ", solvationEnergy, nGKInteractions));
        }

        if (relativeSolvationTerm) {
            sb.append(String.format("  %s %16.8f %12d\n",
                    "Relative Solvation", relativeSolvationEnergy, nRelativeSolvations));
        }

        if (esvTerm) {
            sb.append(String.format("  %s %16.8f  %s\n",
                    "ExtendedSystemBias", esvBias, esvSystem.getLambdaList()));
            sb.append(esvSystem.getBiasDecomposition());
        }

        sb.append(String.format("  %s %16.8f  %s %12.3f (sec)",
                "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));

        int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
        if (nsymm > 1) {
            sb.append(String.format("\n  %s %16.8f", "Unit Cell         ",
                    totalEnergy * nsymm));
        }
        if (crystal.getUnitCell() != crystal) {
            nsymm = crystal.spaceGroup.getNumberOfSymOps();
            sb.append(String.format("\n  %s %16.8f", "Replicates Cell   ",
                    totalEnergy * nsymm));
        }
        sb.append("\n");
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
    @Override
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
            if (restrainTerm && !coordRestraints.isEmpty()) {
                //autoCoordRestraint.setLambda(lambda);
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
            String message = String.format("Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
        }
    }

    public void setPrintOverride(boolean set) {
        this.printOverride = set;
    }

    public void setBondedCombined(boolean set) {
        this.printCompact = set;
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
     * {@inheritDoc}
     */
    @Override
    public double energy(double[] x) {
        return energy(x, false);
    }

    @Override
    public double energy(double[] x, boolean verbose) {
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
        double e = energy(false, verbose);
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
        return energyAndGradient(x, g, false);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double energyAndGradient(double x[], double g[], boolean verbose) {
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
        double e = energy(true, verbose);

        // Try block already exists inside energy(boolean, boolean), so only
        // need to try-catch getGradients.
        try {
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
        } catch (EnergyException ex) {
            ex.printStackTrace();
            if (printOnFailure) {
                String timeString = LocalDateTime.now().format(DateTimeFormatter.
                        ofPattern("yyyy_MM_dd-HH_mm_ss"));

                String filename = String.format("%s-ERROR-%s.pdb",
                        FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                        timeString);

                PotentialsFunctions ef = new PotentialsUtils();
                filename = ef.versionFile(filename);
                logger.info(String.format(" Writing on-error snapshot to file %s", filename));
                ef.saveAsPDB(molecularAssembly, new File(filename));
            }
            if (ex.doCauseSevere()) {
                Utilities.printStackTrace(ex);
                logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
            } else {
                ex.printStackTrace();
                throw ex; // Rethrow exception
            }

            throw ex; // Should ordinarily be unreachable.
        }
    }

    /**
     * <p>
     * getGradients</p>
     *
     * @param g an array of double.
     */
    public double[] getGradients(double g[]) {
        assert (g != null);
        double grad[] = new double[3];
        int n = getNumberOfVariables();
        if (g.length < n) {
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
                if (Double.isNaN(gx) || Double.isInfinite(gx)
                        || Double.isNaN(gy) || Double.isInfinite(gy)
                        || Double.isNaN(gz) || Double.isInfinite(gz)) {
                    /*String message = format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                            a.toString(), gx, gy, gz);*/
                    StringBuilder sb = new StringBuilder(format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                            a.toString(), gx, gy, gz));
                    double[] vals = new double[3];
                    a.getVelocity(vals);
                    sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getAcceleration(vals);
                    sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getPreviousAcceleration(vals);
                    sb.append(format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));

                    //logger.severe(message);
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
     * @param coords
     */
    protected void setCoordinates(double coords[]) {
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

    public void checkAtoms() {
        double vel[] = new double[3];
        double accel[] = new double[3];
        double grad[] = new double[3];
        StringBuilder sb;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                sb = new StringBuilder();
                sb.append("Atom: " + a + "\n");
                try {
                    sb.append("   XYZ :  " + a.getX() + ", " + a.getY() + ", " + a.getZ() + "\n");
                    a.getVelocity(vel);
                    sb.append("   Vel:   " + vel[0] + ", " + vel[1] + ", " + vel[2] + "\n");
                    a.getVelocity(accel);
                    sb.append("   Accel: " + accel[0] + ", " + accel[1] + ", " + accel[2] + "\n");
                    a.getXYZGradient(grad);
                    sb.append("   Grad:  " + grad[0] + ", " + grad[1] + ", " + grad[2] + "\n");
                    sb.append("   Mass:  " + a.getMass() + "\n");
                    if (atoms[i].getEsv() != null) {
                        sb.append(atoms[i].getEsv().toString());
                    }
                } catch (Exception e) {
                }
                logger.info(sb.toString());
            }
        }
    }

    /**
     * {@inheritDoc}
     *
     * @param x
     */
    @Override
    public double[] getCoordinates(double x[]) {
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
        double mass[] = new double[n];
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
                //dEdLambda += autoCoordRestraint.getdEdL();
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
     *
     * @param gradients
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
            if (restrainTerm && !coordRestraints.isEmpty()) {
                for (CoordRestraint restraint : coordRestraints) {
                    restraint.getdEdXdL(gradients);
                }
                //autoCoordRestraint.getdEdXdL(gradients);
            }
            if (comTerm && comRestraint != null) {
                comRestraint.getdEdXdL(gradients);
            }
            if (lambdaTorsions) {
                double grad[] = new double[3];
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
                //d2EdLambda2 += autoCoordRestraint.getd2EdL2();
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
     * @param onFail To set
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
                 * If the call was successful, the property was explicitly set
                 * somewhere and should be kept. If an exception was thrown, the
                 * property was never set explicitly, so over-write.
                 */
            } catch (Exception ex) {
                printOnFailure = onFail;
            }
        }
    }

    public boolean getPrintOnFailure() {
        return printOnFailure;
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

    @Override
    public STATE getEnergyTermState() {
        return state;
    }

    /**
     * This method is for the RESPA integrator only.
     *
     * @param state The STATE is FAST, SLOW or BOTH.
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
     * Set the boundary conditions for this calculation.
     *
     * @param crystal the Crystal contains symmetry and PBC conditions.
     */
    @Override
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

    public void destroy() throws Exception {
        if (parallelTeam != null) {
            try {
                parallelTeam.shutdown();
            } catch (Exception ex) {
                String message = " Error in shutting down the ParallelTeam.";
                logger.log(Level.WARNING, message, ex);
            }
        }
        if (vanderWaals != null) {
            vanderWaals.destroy();
        }
        if (particleMeshEwald != null) {
            particleMeshEwald.destroy();
        }
    }

    public int getNumberofAtoms() {
        return nAtoms;
    }

    public double getBondEnergy() {
        return bondEnergy;
    }

    public int getNumberofBonds() {
        return nBonds;
    }

    public double getAngleEnergy() {
        return angleEnergy;
    }

    public int getNumberofAngles() {
        return nAngles;
    }

    public double getStrenchBendEnergy() {
        return stretchBendEnergy;
    }

    public int getNumberofStretchBends() {
        return nStretchBends;
    }

    public double getUreyBradleyEnergy() {
        return ureyBradleyEnergy;
    }

    public int getNumberofUreyBradleys() {
        return nUreyBradleys;
    }

    public double getOutOfPlaneBendEnergy() {
        return outOfPlaneBendEnergy;
    }

    public int getNumberofOutOfPlaneBends() {
        return nOutOfPlaneBends;
    }

    public double getTorsionEnergy() {
        return torsionEnergy;
    }

    public int getNumberofImproperTorsions() {
        return nImproperTorsions;
    }

    public double getImproperTorsionEnergy() {
        return improperTorsionEnergy;
    }

    public int getNumberofTorsions() {
        return nTorsions;
    }

    public double getPiOrbitalTorsionEnergy() {
        return piOrbitalTorsionEnergy;
    }

    public int getNumberofPiOrbitalTorsions() {
        return nPiOrbitalTorsions;
    }

    public double getTorsionTorsionEnergy() {
        return torsionTorsionEnergy;
    }

    public int getNumberofTorsionTorsions() {
        return nTorsionTorsions;
    }

    public double getVanDerWaalsEnergy() {
        return vanDerWaalsEnergy;
    }

    public int getVanDerWaalsInteractions() {
        return nVanDerWaalInteractions;
    }

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

    public double getPermanentMultipoleEnergy() {
        return permanentMultipoleEnergy;
    }

    public double getPermanentRealSpaceEnergy() {
        return permanentRealSpaceEnergy;
    }

    public double getPermanentReciprocalSelfEnergy() {
        return particleMeshEwald.getPermSelfEnergy();
    }

    public double getPermanentReciprocalMpoleEnergy() {
        return particleMeshEwald.getPermRecipEnergy();
    }

    public int getPermanentInteractions() {
        return nPermanentInteractions;
    }

    public double getPolarizationEnergy() {
        return polarizationEnergy;
    }

    /**
     * Returns total electrostatic energy
     *
     * @param includeSolvation Whether to include solvation energy
     * @return Electrostatic energy
     */
    public double getTotalElectrostaticEnergy(boolean includeSolvation) {
        return (includeSolvation ? getTotalElectrostaticEnergy() : totalMultipoleEnergy);
    }

    public double getTotalElectrostaticEnergy() {
        return totalMultipoleEnergy + solvationEnergy;
    }

    public double getElectrostaticEnergy() {
        return totalMultipoleEnergy;
    }

    public double getSolvationEnergy() {
        return solvationEnergy;
    }

    public double getCavitationEnergy() {
        return particleMeshEwald.getCavitationEnergy(false);
    }

    public double getDispersionEnergy() {
        return particleMeshEwald.getDispersionEnergy(false);
    }

    public double getEsvBiasEnergy() {
        return esvBias;
    }

    public ExtendedSystem getExtendedSystem() {
        return esvSystem;
    }

    public GeneralizedKirkwood getGK() {
        if (particleMeshEwald != null) {
            return particleMeshEwald.getGK();
        } else {
            return null;
        }
    }

    public VanDerWaals getVdwNode() {
        return vanderWaals;
    }

    public ParticleMeshEwald getPmeNode() {
        return particleMeshEwald;
    }

    public ParticleMeshEwaldQI getPmeQiNode() {
        if (!(particleMeshEwald instanceof ParticleMeshEwaldQI)) {
            throw new IllegalStateException();
        }
        return (ParticleMeshEwaldQI) particleMeshEwald;
    }

    public int getSolvationInteractions() {
        return nGKInteractions;
    }

    public double getRelativeSolvationEnergy() {
        return relativeSolvationEnergy;
    }

    /**
     * The velocity array should only contain velocity data for active atoms.
     *
     * @param velocity
     */
    @Override
    public void setVelocity(double[] velocity) {
        if (velocity == null) {
            return;
        }
        int index = 0;
        double vel[] = new double[3];
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
     * The acceleration array should only contain acceleration data for active
     * atoms.
     *
     * @param acceleration
     */
    @Override
    public void setAcceleration(double[] acceleration) {
        if (acceleration == null) {
            return;
        }
        int index = 0;
        double accel[] = new double[3];
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
     * The previousAcceleration array should only contain previous acceleration
     * data for active atoms.
     *
     * @param previousAcceleration
     */
    @Override
    public void setPreviousAcceleration(double[] previousAcceleration) {
        if (previousAcceleration == null) {
            return;
        }
        int index = 0;
        double prev[] = new double[3];
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
     * Returns an array of velocity values for active atoms.
     *
     * @param velocity if the velocity array is null, it will be allocated.
     * @return the velocity array is returned.
     */
    @Override
    public double[] getVelocity(double[] velocity) {
        int n = getNumberOfVariables();
        if (velocity == null || velocity.length < n) {
            velocity = new double[n];
        }
        int index = 0;
        double v[] = new double[3];
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
     * Returns an array of acceleration values for active atoms.
     *
     * @param acceleration if the acceleration array is null, it will be
     * allocated.
     * @return the acceleration array is returned.
     */
    @Override
    public double[] getAcceleration(double[] acceleration) {
        int n = getNumberOfVariables();
        if (acceleration == null || acceleration.length < n) {
            acceleration = new double[n];
        }
        int index = 0;
        double a[] = new double[3];
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
     * Returns an array of previous acceleration values for active atoms.
     *
     * @param previousAcceleration if the previousAcceleration array is null, it
     * will be allocated.
     * @return the previousAcceleration array is returned.
     */
    @Override
    public double[] getPreviousAcceleration(double[] previousAcceleration) {
        int n = getNumberOfVariables();
        if (previousAcceleration == null || previousAcceleration.length < n) {
            previousAcceleration = new double[n];
        }
        int index = 0;
        double a[] = new double[3];
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

    public Bond[] getBonds() {
        return bonds;
    }

    public Angle[] getAngles() {
        return angles;
    }

    public ImproperTorsion[] getImproperTorsions() {
        return improperTorsions;
    }

    public UreyBradley[] getUreyBradleys() {
        return ureyBradleys;
    }

    public OutOfPlaneBend[] getOutOfPlaneBends() {
        return outOfPlaneBends;
    }

    public StretchBend[] getStretchBends() {
        return stretchBends;
    }

    public Torsion[] getTorsions() {
        return torsions;
    }

    public PiOrbitalTorsion[] getPiOrbitalTorsions() {
        return piOrbitalTorsions;
    }

    public TorsionTorsion[] getTorsionTorsions() {
        return torsionTorsions;
    }

    private class BondedRegion extends ParallelRegion {

        // Flag to indicate gradient computation.
        private boolean gradient = false;

        private AtomicDoubleArrayImpl atomicDoubleArrayImpl;
        /**
         * X-component of the Cartesian coordinate gradient.
         */
        private final AtomicDoubleArray gradX;
        /**
         * Y-component of the Cartesian coordinate gradient.
         */
        private final AtomicDoubleArray gradY;
        /**
         * Z-component of the Cartesian coordinate gradient.
         */
        private final AtomicDoubleArray gradZ;
        /**
         * X-component of the dU/dX/dL coordinate gradient.
         */
        private final AtomicDoubleArray lambdaGradX;
        /**
         * Y-component of the dU/dX/dL coordinate gradient.
         */
        private final AtomicDoubleArray lambdaGradY;
        /**
         * Z-component of the dU/dX/dL coordinate gradient.
         */
        private final AtomicDoubleArray lambdaGradZ;

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
        private final SharedDouble sharedTorsionTorsionEnergy;
        // Shared restraint terms.
        private final SharedDouble sharedRestraintBondEnergy;

        // Number of threads.
        private final int nThreads;

        // Gradient loops.
        private final GradInitLoop gradInitLoops[];
        private final GradReduceLoop gradReduceLoops[];

        // Force field bonded energy parallel loops.
        private final BondedTermLoop[] bondLoops;
        private final BondedTermLoop[] angleLoops;
        private final BondedTermLoop[] outOfPlaneBendLoops;
        private final BondedTermLoop[] improperTorsionLoops;
        private final BondedTermLoop[] piOrbitalTorsionLoops;
        private final BondedTermLoop[] stretchBendLoops;
        private final BondedTermLoop[] torsionLoops;
        private final BondedTermLoop[] torsionTorsionLoops;
        private final BondedTermLoop[] ureyBradleyLoops;

        // Retraint energy parallel loops.
        private final BondedTermLoop[] restraintBondLoops;

        public BondedRegion() {

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
            torsionTorsionLoops = new BondedTermLoop[nThreads];
            ureyBradleyLoops = new BondedTermLoop[nThreads];

            // Allocate memory for restrain energy terms.
            restraintBondLoops = new BondedTermLoop[nThreads];

            /**
             * Define how the gradient will be accumulated.
             */
            atomicDoubleArrayImpl = AtomicDoubleArrayImpl.MULTI;
            ForceField forceField = molecularAssembly.getForceField();
            String value = forceField.getString(ARRAY_REDUCTION, "MULTI");
            try {
                atomicDoubleArrayImpl = AtomicDoubleArrayImpl.valueOf(toEnumForm(value));
            } catch (Exception e) {
                logger.info(format(" Unrecognized ARRAY-REDUCTION %s; defaulting to %s", value, atomicDoubleArrayImpl));
            }
            switch (atomicDoubleArrayImpl) {
                case MULTI:
                default:
                    gradX = new MultiDoubleArray(nThreads, nAtoms);
                    gradY = new MultiDoubleArray(nThreads, nAtoms);
                    gradZ = new MultiDoubleArray(nThreads, nAtoms);
                    if (lambdaTerm) {
                        lambdaGradX = new MultiDoubleArray(nThreads, nAtoms);
                        lambdaGradY = new MultiDoubleArray(nThreads, nAtoms);
                        lambdaGradZ = new MultiDoubleArray(nThreads, nAtoms);
                    } else {
                        lambdaGradX = null;
                        lambdaGradY = null;
                        lambdaGradZ = null;
                    }
                    break;
                case PJ:
                    gradX = new PJDoubleArray(nThreads, nAtoms);
                    gradY = new PJDoubleArray(nThreads, nAtoms);
                    gradZ = new PJDoubleArray(nThreads, nAtoms);
                    if (lambdaTerm) {
                        lambdaGradX = new PJDoubleArray(nThreads, nAtoms);
                        lambdaGradY = new PJDoubleArray(nThreads, nAtoms);
                        lambdaGradZ = new PJDoubleArray(nThreads, nAtoms);
                    } else {
                        lambdaGradX = null;
                        lambdaGradY = null;
                        lambdaGradZ = null;
                    }
                    break;
                case ADDER:
                    gradX = new AdderDoubleArray(nThreads, nAtoms);
                    gradY = new AdderDoubleArray(nThreads, nAtoms);
                    gradZ = new AdderDoubleArray(nThreads, nAtoms);
                    if (lambdaTerm) {
                        lambdaGradX = new AdderDoubleArray(nThreads, nAtoms);
                        lambdaGradY = new AdderDoubleArray(nThreads, nAtoms);
                        lambdaGradZ = new AdderDoubleArray(nThreads, nAtoms);
                    } else {
                        lambdaGradX = null;
                        lambdaGradY = null;
                        lambdaGradZ = null;
                    }
                    break;
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
            sharedTorsionTorsionEnergy.set(0.0);
            sharedUreyBradleyEnergy.set(0.0);

            // Zero out shared restraint energy values.
            sharedRestraintBondEnergy.set(0.0);

            // Assure capacity of the gradient arrays.
            if (gradient) {
                gradX.alloc(nAtoms);
                gradY.alloc(nAtoms);
                gradZ.alloc(nAtoms);
            }
            if (lambdaTerm) {
                lambdaGradX.alloc(nAtoms);
                lambdaGradY.alloc(nAtoms);
                lambdaGradZ.alloc(nAtoms);
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
                    angleLoops[threadID]
                            = new BondedTermLoop(angles, sharedAngleEnergy, sharedAngleRMSD);
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
                    bondLoops[threadID]
                            = new BondedTermLoop(bonds, sharedBondEnergy, sharedBondRMSD);
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
                    improperTorsionLoops[threadID]
                            = new BondedTermLoop(improperTorsions, sharedImproperTorsionEnergy);
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
                    outOfPlaneBendLoops[threadID]
                            = new BondedTermLoop(outOfPlaneBends, sharedOutOfPlaneBendEnergy);
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
                    piOrbitalTorsionLoops[threadID]
                            = new BondedTermLoop(piOrbitalTorsions, sharedPiOrbitalTorsionEnergy);
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
                    stretchBendLoops[threadID]
                            = new BondedTermLoop(stretchBends, sharedStretchBendEnergy);
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
                    torsionLoops[threadID]
                            = new BondedTermLoop(torsions, sharedTorsionEnergy);
                }
                if (threadID == 0) {
                    torsionTime = -System.nanoTime();
                }
                execute(0, nTorsions - 1, torsionLoops[threadID]);
                if (threadID == 0) {
                    torsionTime += System.nanoTime();
                }
            }

            if (torsionTorsionTerm) {
                if (torsionTorsionLoops[threadID] == null) {
                    torsionTorsionLoops[threadID]
                            = new BondedTermLoop(torsionTorsions, sharedTorsionTorsionEnergy);
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
                    ureyBradleyLoops[threadID]
                            = new BondedTermLoop(ureyBradleys, sharedUreyBradleyEnergy);
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
                    restraintBondLoops[threadID]
                            = new BondedTermLoop(restraintBonds, sharedRestraintBondEnergy);
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
                    gradX.reset(threadID, first, last);
                    gradY.reset(threadID, first, last);
                    gradZ.reset(threadID, first, last);
                }
                if (lambdaTerm) {
                    lambdaGradX.reset(threadID, first, last);
                    lambdaGradY.reset(threadID, first, last);
                    lambdaGradZ.reset(threadID, first, last);
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
                    gradX.reduce(first, last);
                    gradY.reduce(first, last);
                    gradZ.reduce(first, last);
                    for (int i = first; i <= last; i++) {
                        Atom a = atoms[i];
                        a.setXYZGradient(gradX.get(i), gradY.get(i), gradZ.get(i));
                    }
                }
                if (lambdaTerm) {
                    lambdaGradX.reduce(first, last);
                    lambdaGradY.reduce(first, last);
                    lambdaGradZ.reduce(first, last);
                    for (int i = first; i <= last; i++) {
                        Atom a = atoms[i];
                        a.setLambdaXYZGradient(lambdaGradX.get(i),
                                lambdaGradY.get(i), lambdaGradZ.get(i));
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

            public BondedTermLoop(BondedTerm[] terms, SharedDouble sharedEnergy) {
                this(terms, sharedEnergy, null);
            }

            public BondedTermLoop(BondedTerm[] terms,
                    SharedDouble sharedEnergy, SharedDouble sharedRMSD) {
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
                        localEnergy += term.energy(gradient, threadID,
                                gradX, gradY, gradZ,
                                lambdaGradX, lambdaGradY, lambdaGradZ);
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
     * An EnergyImplementation describes a set of force field implementations;
     * currently FFX for the pure Java reference implementation, and various
     * OpenMM implementations, of which OMM_CUDA is preferred. Note that AMOEBA
     * only runs under the FFX, OMM_CUDA, and OMM_REF implementations.
     *
     * FFX: reference FFX implementation OMM: Currently an alias for OMM_CUDA,
     * may eventually become "try to find best OpenMM implementation" OMM_CUDA:
     * OpenMM CUDA implementation OMM_REF: OpenMM reference implementation
     * OMM_OPTCPU: Optimized OpenMM CPU implementation (no AMOEBA) OMM_OPENCL:
     * OpenMM OpenCL implementation (no AMOEBA)
     */
    public static enum EnergyImplementation {
        FFX, OMM, OMM_CUDA, OMM_REF, OMM_OPTCPU, OMM_OPENCL;
    }
}
