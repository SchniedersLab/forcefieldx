// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
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
// ******************************************************************************
package ffx.potential;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.crystal.CrystalPotential;
import ffx.crystal.LatticeSystem;
import ffx.crystal.NCSCrystal;
import ffx.crystal.ReplicatesCrystal;
import ffx.crystal.SpaceGroup;
import ffx.crystal.SpaceGroupDefinitions;
import ffx.crystal.SymOp;
import ffx.numerics.Constraint;
import ffx.numerics.math.Double3;
import ffx.numerics.switching.ConstantSwitch;
import ffx.numerics.switching.UnivariateFunctionFactory;
import ffx.numerics.switching.UnivariateSwitchingFunction;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.RestrainDistance;
import ffx.potential.bonded.RestrainPosition;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.StretchTorsion;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.constraint.CcmaConstraint;
import ffx.potential.constraint.SettleConstraint;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.COMRestraint;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.NCSRestraint;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.RestrainGroups;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsTornado;
import ffx.potential.openmm.OpenMMEnergy;
import ffx.potential.parameters.AngleType.AngleMode;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ELEC_FORM;
import ffx.potential.parameters.TorsionType;
import ffx.potential.terms.AnglePotentialEnergy;
import ffx.potential.terms.AngleTorsionPotentialEnergy;
import ffx.potential.terms.BondPotentialEnergy;
import ffx.potential.terms.EnergyTermRegion;
import ffx.potential.terms.ImproperTorsionPotentialEnergy;
import ffx.potential.terms.OutOfPlaneBendPotentialEnergy;
import ffx.potential.terms.PiOrbitalTorsionPotentialEnergy;
import ffx.potential.terms.RestrainDistancePotentialEnergy;
import ffx.potential.terms.RestrainPositionPotentialEnergy;
import ffx.potential.terms.RestrainTorsionPotentialEnergy;
import ffx.potential.terms.StretchBendPotentialEnergy;
import ffx.potential.terms.StretchTorsionPotentialEnergy;
import ffx.potential.terms.TorsionPotentialEnergy;
import ffx.potential.terms.TorsionTorsionPotentialEnergy;
import ffx.potential.terms.UreyBradleyPotentialEnergy;
import ffx.potential.utils.ConvexHullOps;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.FFXProperty;
import org.apache.commons.configuration2.CompositeConfiguration;

import javax.annotation.Nullable;
import java.io.File;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static ffx.potential.bonded.BondedTerm.removeNeuralNetworkTerms;
import static ffx.potential.bonded.RestrainPosition.parseRestrainPositions;
import static ffx.potential.nonbonded.VanDerWaalsForm.DEFAULT_VDW_CUTOFF;
import static ffx.potential.nonbonded.pme.EwaldParameters.DEFAULT_EWALD_CUTOFF;
import static ffx.potential.parameters.ForceField.toEnumForm;
import static ffx.potential.parsers.XYZFileFilter.isXYZ;
import static ffx.utilities.PropertyGroup.NonBondedCutoff;
import static ffx.utilities.PropertyGroup.PotentialFunctionSelection;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static org.apache.commons.io.FilenameUtils.removeExtension;
import static org.apache.commons.math3.util.FastMath.max;
import static org.apache.commons.math3.util.FastMath.min;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * Compute the potential energy and derivatives of a molecular system described by a force field.
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
   * RestrainMode specifies how restrain terms are applied for dual topology calculations.
   */
  public enum RestrainMode {
    /**
     * Restrain terms are applied with the energy.
     */
    ENERGY,
    /**
     * Restrain terms are applied with the alchemical restraints (Lambda Bonded Terms).
     */
    ALCHEMICAL
  }

  /**
   * The ForceFieldEnergy class implements the FFX platform.
   */
  private final Platform platform = Platform.FFX;

  /**
   * The MolecularAssembly associated with this force field energy.
   */
  protected final MolecularAssembly molecularAssembly;
  /**
   * The array of Atoms being evaluated.
   */
  private final Atom[] atoms;
  /**
   * The boundary conditions used when evaluating the force field energy.
   */
  private Crystal crystal;
  /**
   * The Parallel Java ParallelTeam instance.
   */
  private final ParallelTeam parallelTeam;
  /**
   * An NCS restraint term.
   */
  private final NCSRestraint ncsRestraint;
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
  private final ParticleMeshEwald particleMeshEwald;
  /**
   * ANI2x Neutral Network Potential.
   */
  private final ANIEnergy aniEnergy;

  /**
   * List of constraints to be applied to the system.
   */
  private final List<Constraint> constraints;
  /**
   * Default tolerance for numerical methods of solving constraints.
   */
  public static final double DEFAULT_CONSTRAINT_TOLERANCE = 1E-4;

  @FFXProperty(name = "vdw-cutoff", propertyGroup = NonBondedCutoff, defaultValue = "12.0", description = """
      Sets the cutoff distance value in Angstroms for van der Waals potential energy interactions.
      The energy for any pair of van der Waals sites beyond the cutoff distance will be set to zero.
      Other properties can be used to define the smoothing scheme near the cutoff distance.
      The default cutoff distance in the absence of the vdw-cutoff keyword is infinite for nonperiodic
      systems and 12.0 for periodic systems.
      """)
  private double vdwCutoff;

  @FFXProperty(name = "ewald-cutoff", propertyGroup = NonBondedCutoff, defaultValue = "7.0", description = """
      Sets the value in Angstroms of the real-space distance cutoff for use during Ewald summation.
      By default, in the absence of the ewald-cutoff property, a value of 7.0 is used.
      """)
  private double ewaldCutoff;

  @FFXProperty(name = "gk-cutoff", propertyGroup = NonBondedCutoff, defaultValue = "12.0", description = """
      Sets the value in Angstroms of the generalized Kirkwood distance cutoff for use
      during implicit solvent simulations. By default, in the absence of the gk-cutoff property,
      no cutoff is used under aperiodic boundary conditions and the vdw-cutoff is used under PBC.
      """)
  private double gkCutoff;

  private static final double DEFAULT_LIST_BUFFER = 2.0;

  @FFXProperty(name = "list-buffer", propertyGroup = NonBondedCutoff, defaultValue = "2.0", description = """
      Sets the size of the neighbor list buffer in Angstroms for potential energy functions.
      This value is added to the actual cutoff distance to determine which pairs will be kept on the neighbor list.
      This buffer value is used for all potential function neighbor lists.
      The default value in the absence of the list-buffer keyword is 2.0 Angstroms.
      """)
  private double listBuffer;

  /**
   * 2.0 times the neighbor list cutoff.
   */
  private final double cutOff2;
  /**
   * The non-bonded cut-off plus buffer distance (Angstroms).
   */
  private final double cutoffPlusBuffer;
  /**
   * The Electrostatic functional form.
   */
  private final ELEC_FORM elecForm;

  /**
   * Indicates use of the Lambda state variable.
   */
  @FFXProperty(name = "lambdaterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "false", description = "Specifies use of the Lambda state variable.")
  protected boolean lambdaTerm;
  /**
   * Current value of the Lambda state variable.
   */
  private double lambda = 1.0;
  /**
   * Indicates only bonded energy terms effected by Lambda should be evaluated.
   */
  public boolean lambdaBondedTerms = false;
  /**
   * Indicates all bonded energy terms should be evaluated if lambdaBondedTerms is true.
   */
  boolean lambdaAllBondedTerms = false;
  /**
   * Indicates application of lambda scaling to all Torsion based energy terms.
   */
  private final boolean lambdaTorsions;

  /**
   * An instance of the STATE enumeration to specify calculation of slowly varying energy terms, fast
   * varying or both.
   */
  private STATE state = STATE.BOTH;
  /**
   * Optimization scaling value to use for each degree of freedom.
   */
  protected double[] optimizationScaling = null;

  /**
   * Evaluate bonded energy terms in parallel.
   */
  private final EnergyTermRegion forceFieldBondedEnergyRegion;
  /**
   * Implements the Bond potential.
   */
  private final BondPotentialEnergy bondPotentialEnergy;
  /**
   * Implements the Angle potential.
   */
  private final AnglePotentialEnergy anglePotentialEnergy;
  /**
   * Implements the Stretch-Bend potential.
   */
  private final StretchBendPotentialEnergy stretchBendPotentialEnergy;
  /**
   * Implements the Urey-Bradley potential.
   */
  private final UreyBradleyPotentialEnergy ureyBradleyPotentialEnergy;
  /**
   * Implements the Out-of-Plane Bend potential.
   */
  private final OutOfPlaneBendPotentialEnergy outOfPlaneBendPotentialEnergy;
  /**
   * Implements the Torsion potential.
   */
  private final TorsionPotentialEnergy torsionPotentialEnergy;
  /**
   * Implements the Stretch-Torsion potential.
   */
  private final StretchTorsionPotentialEnergy stretchTorsionPotentialEnergy;
  /**
   * Implements the Angle-Torsion potential.
   */
  private final AngleTorsionPotentialEnergy angleTorsionPotentialEnergy;
  /**
   * Implements the Improper Torsion potential.
   */
  private final ImproperTorsionPotentialEnergy improperTorsionPotentialEnergy;
  /**
   * Implements the Pi-Orbital Torsion potential.
   */
  private final PiOrbitalTorsionPotentialEnergy piOrbitalTorsionPotentialEnergy;
  /**
   * Implements the Torsion-Torsion potential.
   */
  private final TorsionTorsionPotentialEnergy torsionTorsionPotentialEnergy;

  /**
   * Evaluate bonded energy terms in parallel.
   */
  private final EnergyTermRegion restrainEnergyRegion;
  /**
   * RestrainMode specifies how restrain terms are applied for dual topology calculations.
   */
  private final RestrainMode restrainMode;
  /**
   * Implements the Restrain-Position potential.
   */
  private final RestrainPositionPotentialEnergy restrainPositionPotentialEnergy;
  /**
   * Implements the Restrain-Distance potential.
   */
  private final RestrainDistancePotentialEnergy restrainDistancePotentialEnergy;
  /**
   * Implements the Restrain-Torsion potential.
   */
  private final RestrainTorsionPotentialEnergy restrainTorsionPotentialEnergy;
  /**
   * The Restrain Groups energy term.
   */
  private final RestrainGroups restrainGroups;

  /**
   * Number of atoms in the system.
   */
  private final int nAtoms;
  /**
   * Number of Restrain Groups in the system.
   */
  private int nRestrainGroups;
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
   * Specifies use of an intramolecular neural network.
   */
  private boolean nnTerm;

  /**
   * Specifies use of the bond stretch potential.
   */
  @FFXProperty(name = "bondterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the bond stretch potential.")
  private final boolean bondTerm;

  /**
   * Specifies use of the angle bend potential.
   */
  @FFXProperty(name = "angleterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the angle bend potential.")
  private final boolean angleTerm;

  /**
   * Evaluate Stretch-Bend energy terms.
   */
  @FFXProperty(name = "strbndterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the stretch-bend potential.")
  private final boolean stretchBendTerm;

  /**
   * Evaluate Urey-Bradley energy terms.
   */
  @FFXProperty(name = "ureyterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the Urey-Bradley potential.")
  private final boolean ureyBradleyTerm;

  /**
   * Evaluate Out of Plane Bend energy terms.
   */
  @FFXProperty(name = "opbendterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the out-of-plane potential.")
  private final boolean outOfPlaneBendTerm;

  /**
   * Evaluate Torsion energy terms.
   */
  @FFXProperty(name = "torsionterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the torsional potential.")
  private final boolean torsionTerm;

  /**
   * Evaluate Stretch-Torsion energy terms.
   */
  @FFXProperty(name = "strtorterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the stretch-torsion potential.")
  private final boolean stretchTorsionTerm;

  /**
   * Evaluate Angle-Torsion energy terms.
   */
  @FFXProperty(name = "angtorterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the angle-torsion potential.")
  private final boolean angleTorsionTerm;

  /**
   * Evaluate Improper Torsion energy terms.
   */
  @FFXProperty(name = "imptorterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the improper torsion potential.")
  private final boolean improperTorsionTerm;

  /**
   * Evaluate Pi-Orbital Torsion energy terms.
   */
  @FFXProperty(name = "pitorsterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the pi-system torsion potential.")
  private final boolean piOrbitalTorsionTerm;

  /**
   * Evaluate Torsion-Torsion energy terms.
   */
  @FFXProperty(name = "tortorterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = "Specifies use of the pi-system torsion potential.")
  private final boolean torsionTorsionTerm;

  /**
   * Evaluate RestrainPosition term.
   */
  private final boolean restrainPositionTerm;
  /**
   * Evaluate RestrainDistance terms.
   */
  private final boolean restrainDistanceTerm;
  /**
   * Evaluate RestrainTorsion terms.
   */
  private final boolean restrainTorsionTerm;
  /**
   * Evaluate RestrainGroup terms.
   */
  private boolean restrainGroupTerm;

  /**
   * Evaluate van der Waals energy term.
   */
  @FFXProperty(name = "vdwterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = """
      Specifies use of the vdw der Waals potential.
      If set to false, all non-bonded terms are turned off.
      """)
  private boolean vanderWaalsTerm;

  /**
   * Evaluate permanent multipole electrostatics energy term.
   */
  @FFXProperty(name = "mpoleterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = """
      Specifies use of the fixed charge electrostatic potential.
      Setting mpoleterm to false also turns off polarization and generalized Kirkwood,
      overriding the polarizeterm and gkterm properties.
      """)
  private boolean multipoleTerm;

  /**
   * Evaluate polarization energy term.
   */
  @FFXProperty(name = "polarizeterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "true", description = """
      Specifies use of the polarizable electrostatic potential.
      Setting polarizeterm to false overrides the polarization property.
      """)
  private boolean polarizationTerm;

  /**
   * Evaluate generalized Kirkwood energy term.
   */
  @FFXProperty(name = "gkterm", clazz = Boolean.class, propertyGroup = PotentialFunctionSelection,
      defaultValue = "false", description = "Specifies use of generalized Kirkwood electrostatics.")
  private boolean generalizedKirkwoodTerm;

  /**
   * Evaluate COM energy term.
   */
  private boolean comTerm;
  /**
   * Evaluate NCS energy term.
   */
  private boolean ncsTerm;

  /**
   * Original state of the neural network energy term flag.
   */
  private final boolean nnTermOrig;
  /**
   * Original state of the RestrainGroup term flag.
   */
  private final boolean restrainGroupTermOrig;
  /**
   * Original state of the van der Waals energy term flag.
   */
  private final boolean vanderWaalsTermOrig;
  /**
   * Original state of the multipole energy term flag.
   */
  private final boolean multipoleTermOrig;
  /**
   * Original state of the polarization energy term flag.
   */
  private final boolean polarizationTermOrig;
  /**
   * Original state of the GK energy term flag.
   */
  private final boolean generalizedKirkwoodTermOrig;
  /**
   * Original state of the ESV energy term flag.
   */
  private boolean esvTermOrig;
  /**
   * Original state of the NCS energy term flag.
   */
  private boolean ncsTermOrig;
  /**
   * Original state of the COM energy term flag.
   */
  private boolean comTermOrig;

  /**
   * The total neutral network energy.
   */
  private double nnEnergy;
  /**
   * The RestrainGroup Energy.
   */
  private double restrainGroupEnergy;
  /**
   * The total van der Waals energy.
   */
  private double vanDerWaalsEnergy;
  /**
   * The permanent Multipole energy.
   */
  private double permanentMultipoleEnergy;
  /**
   * The total polarization energy.
   */
  private double polarizationEnergy;
  /**
   * The solvation energy (GK).
   */
  private double solvationEnergy;
  /**
   * The total NCS Energy.
   */
  private double ncsEnergy;
  /**
   * The total COM Restraint Energy.
   */
  private double comRestraintEnergy;
  /**
   * The total system energy.
   */
  private double totalEnergy;

  /**
   * Time to evaluate NCS term.
   */
  private long ncsTime;
  /**
   * Time to evaluate the neutral network.
   */
  private long nnTime;
  /**
   * Time to evaluate restrain group term.
   */
  private long restrainGroupTime;
  /**
   * Time to evaluate Center of Mass restraint term.
   */
  private long comRestraintTime;
  /**
   * Time to evaluate van der Waals term.
   */
  private long vanDerWaalsTime;
  /**
   * Time to evaluate electrostatics term.
   */
  private long electrostaticTime;
  /**
   * Time to evaluate all energy terms.
   */
  private long totalTime;

  /**
   * Constant pH extended system (TODO: needs further testing).
   */
  private ExtendedSystem esvSystem = null;
  private boolean esvTerm;
  private double esvBias;

  /**
   * Flag to indicate proper shutdown of the ForceFieldEnergy.
   */
  boolean destroyed = false;
  /**
   * If the absolute value of a gradient component is greater than "maxDebugGradient", verbose
   * logging results.
   */
  public final double maxDebugGradient;
  /**
   * Enable verbose printing if large energy gradient components are observed.
   */
  private boolean printOnFailure;

  /**
   * Constructor for ForceFieldEnergy.
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   */
  protected ForceFieldEnergy(MolecularAssembly molecularAssembly) {
    this(molecularAssembly, ParallelTeam.getDefaultThreadCount());
  }

  /**
   * Constructor for ForceFieldEnergy.
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param numThreads        a int.
   */
  protected ForceFieldEnergy(MolecularAssembly molecularAssembly, int numThreads) {
    // Get a reference to the sorted atom array.
    this.molecularAssembly = molecularAssembly;
    atoms = molecularAssembly.getAtomArray();
    nAtoms = atoms.length;

    // Check that atom ordering is correct and count the number of active atoms.
    for (int i = 0; i < nAtoms; i++) {
      int index = atoms[i].getXyzIndex() - 1;
      assert (i == index);
    }

    // Enforce that the number of threads be less than or equal to the number of atoms.
    int nThreads = min(nAtoms, numThreads);
    parallelTeam = new ParallelTeam(nThreads);

    ForceField forceField = molecularAssembly.getForceField();
    CompositeConfiguration properties = forceField.getProperties();
    String name = forceField.toString().toUpperCase();

    logger.info(format(" Constructing Force Field %s", name));
    logger.info(format("\n SMP threads:                        %10d", nThreads));

    // Rename water protons if requested.
    boolean standardizeAtomNames = properties.getBoolean("standardizeAtomNames", true);
    if (standardizeAtomNames) {
      molecularAssembly.renameWaterProtons();
    }

    nnTerm = forceField.getBoolean("NNTERM", false);
    // For now, all atoms are neural network atoms, or none are.
    if (nnTerm) {
      for (Atom atom : atoms) {
        atom.setNeuralNetwork(true);
      }
    }
    bondTerm = forceField.getBoolean("BONDTERM", true);
    angleTerm = forceField.getBoolean("ANGLETERM", true);
    stretchBendTerm = forceField.getBoolean("STRBNDTERM", true);
    ureyBradleyTerm = forceField.getBoolean("UREYTERM", true);
    outOfPlaneBendTerm = forceField.getBoolean("OPBENDTERM", true);
    torsionTerm = forceField.getBoolean("TORSIONTERM", true);
    stretchTorsionTerm = forceField.getBoolean("STRTORTERM", true);
    angleTorsionTerm = forceField.getBoolean("ANGTORTERM", true);
    piOrbitalTorsionTerm = forceField.getBoolean("PITORSTERM", true);
    torsionTorsionTerm = forceField.getBoolean("TORTORTERM", true);
    improperTorsionTerm = forceField.getBoolean("IMPROPERTERM", true);
    vanderWaalsTerm = forceField.getBoolean("VDWTERM", true);

    boolean tornadoVM = forceField.getBoolean("tornado", false);
    if (vanderWaalsTerm && !tornadoVM) {
      multipoleTerm = forceField.getBoolean("MPOLETERM", true);
      if (multipoleTerm) {
        String polarizeString = forceField.getString("POLARIZATION", "NONE");
        boolean defaultPolarizeTerm = !polarizeString.equalsIgnoreCase("NONE");
        polarizationTerm = forceField.getBoolean("POLARIZETERM", defaultPolarizeTerm);
        generalizedKirkwoodTerm = forceField.getBoolean("GKTERM", false);
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

    lambdaTerm = forceField.getBoolean("LAMBDATERM", false);
    comTerm = forceField.getBoolean("COMRESTRAINTERM", false);
    lambdaTorsions = forceField.getBoolean("TORSION_LAMBDATERM", false);
    printOnFailure = forceField.getBoolean("PRINT_ON_FAILURE", false);

    String mode = forceField.getString("RESTRAIN_MODE", "ENERGY");
    if (mode.equalsIgnoreCase("ALCHEMICAL")) {
      restrainMode = RestrainMode.ALCHEMICAL;
    } else {
      restrainMode = RestrainMode.ENERGY;
    }

    // Detect "restrain-groups" property.
    if (properties.containsKey("restrain-groups")) {
      restrainGroupTerm = true;
    } else {
      restrainGroupTerm = false;
    }

    // For RESPA
    nnTermOrig = nnTerm;
    restrainGroupTermOrig = restrainGroupTerm;
    vanderWaalsTermOrig = vanderWaalsTerm;
    multipoleTermOrig = multipoleTerm;
    polarizationTermOrig = polarizationTerm;
    generalizedKirkwoodTermOrig = generalizedKirkwoodTerm;
    ncsTermOrig = ncsTerm;
    comTermOrig = comTerm;

    // Determine the unit cell dimensions and space group.
    String spacegroup;
    double a, b, c, alpha, beta, gamma;
    boolean aperiodic;
    try {
      spacegroup = forceField.getString("SPACEGROUP", "P 1");
      SpaceGroup sg = SpaceGroupDefinitions.spaceGroupFactory(spacegroup);
      LatticeSystem latticeSystem = sg.latticeSystem;
      a = forceField.getDouble("A_AXIS");
      aperiodic = false;
      b = forceField.getDouble("B_AXIS", latticeSystem.getDefaultBAxis(a));
      c = forceField.getDouble("C_AXIS", latticeSystem.getDefaultCAxis(a, b));
      alpha = forceField.getDouble("ALPHA", latticeSystem.getDefaultAlpha());
      beta = forceField.getDouble("BETA", latticeSystem.getDefaultBeta());
      gamma = forceField.getDouble("GAMMA", latticeSystem.getDefaultGamma());
      if (!sg.latticeSystem.validParameters(a, b, c, alpha, beta, gamma)) {
        logger.severe(" Check lattice parameters.");
      }
      if (a == 1.0 && b == 1.0 && c == 1.0) {
        String message = " A-, B-, and C-axis values equal to 1.0.";
        logger.info(message);
        throw new Exception(message);
      }
      if (a <= 0.0 || b <= 0.0 || c <= 0.0 || alpha <= 0.0 || beta <= 0.0 || gamma <= 0.0) {
        // Parameters are not valid. Switch to aperiodic.
        String message = " Crystal parameters are not valid due to negative or zero value.";
        logger.warning(message);
        throw new Exception(message);
      }
    } catch (Exception e) {
      aperiodic = true;

      // Determine the maximum separation between atoms.
      double maxR = 0.0;
      if (nAtoms < 10) {
        for (int i = 0; i < nAtoms - 1; i++) {
          Double3 xi = atoms[i].getXYZ();
          for (int j = 1; j < nAtoms; j++) {
            double r = atoms[j].getXYZ().dist(xi);
            maxR = max(r, maxR);
          }
        }
      } else {
        try {
          maxR = ConvexHullOps.maxDist(ConvexHullOps.constructHull(atoms));
        } catch (Exception ex) {
          // If Convex Hull approach fails (e.g., coplanar input, brute force...
          logger.info(" Convex Hull operation failed with message " + ex + "\n Trying brute force approach...");
          if (logger.isLoggable(Level.FINE)) {
            logger.fine(Utilities.stackTraceToString(ex));
          }
          for (int i = 0; i < nAtoms - 1; i++) {
            Double3 xi = atoms[i].getXYZ();
            for (int j = 1; j < nAtoms; j++) {
              double r = atoms[j].getXYZ().dist(xi);
              maxR = max(r, maxR);
            }
          }
        }
      }
      maxR = max(10.0, maxR);

      logger.info(format(" The system will be treated as aperiodic (max separation = %6.1f A).", maxR));

      // Turn off reciprocal space calculations.
      forceField.addProperty("EWALD_ALPHA", "0.0");

      // Specify some dummy values for the crystal.
      spacegroup = "P1";
      a = 2.0 * maxR;
      b = 2.0 * maxR;
      c = 2.0 * maxR;
      alpha = 90.0;
      beta = 90.0;
      gamma = 90.0;
    }
    Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
    unitCell.setAperiodic(aperiodic);

    // First set cutoffs assuming PBC.
    vdwCutoff = DEFAULT_VDW_CUTOFF;
    ewaldCutoff = DEFAULT_EWALD_CUTOFF;
    gkCutoff = vdwCutoff;
    if (unitCell.aperiodic()) {
      // Default to no cutoffs for aperiodic systems.
      vdwCutoff = Double.POSITIVE_INFINITY;
      ewaldCutoff = Double.POSITIVE_INFINITY;
      gkCutoff = Double.POSITIVE_INFINITY;
    }
    // Load user specified cutoffs.
    vdwCutoff = forceField.getDouble("VDW_CUTOFF", vdwCutoff);
    ewaldCutoff = forceField.getDouble("EWALD_CUTOFF", ewaldCutoff);
    gkCutoff = forceField.getDouble("GK_CUTOFF", gkCutoff);

    // Neighbor list cutoff is the maximum of the vdW, Ewald and GK cutoffs (i.e. to use a single list).
    double nlistCutoff = vdwCutoff;
    if (multipoleTerm) {
      nlistCutoff = max(vdwCutoff, ewaldCutoff);
    }
    if (generalizedKirkwoodTerm) {
      nlistCutoff = max(nlistCutoff, gkCutoff);
    }

    // Check for a frozen neighbor list.
    boolean disabledNeighborUpdates = forceField.getBoolean("DISABLE_NEIGHBOR_UPDATES", false);
    if (disabledNeighborUpdates) {
      logger.info(format(" Neighbor list updates disabled; interactions will "
          + "only be calculated between atoms that started the simulation "
          + "within a radius of %9.3g Angstroms of each other.", nlistCutoff));
      // The individual cutoffs are infinity (i.e. they are limited by interactions in the list).
      vdwCutoff = Double.POSITIVE_INFINITY;
      ewaldCutoff = Double.POSITIVE_INFINITY;
      gkCutoff = Double.POSITIVE_INFINITY;
    }

    listBuffer = forceField.getDouble("LIST_BUFFER", DEFAULT_LIST_BUFFER);
    cutoffPlusBuffer = nlistCutoff + listBuffer;
    cutOff2 = 2.0 * cutoffPlusBuffer;
    unitCell = configureNCS(forceField, unitCell);

    // If necessary, create a ReplicatesCrystal.
    if (!aperiodic) {
      this.crystal = ReplicatesCrystal.replicatesCrystalFactory(unitCell, cutOff2);
      logger.info(format("\n Density:                                %6.3f (g/cc)",
          crystal.getDensity(molecularAssembly.getMass())));
      logger.info(crystal.toString());
    } else {
      this.crystal = unitCell;
    }

    if (!unitCell.aperiodic() && unitCell.spaceGroup.number == 1) {
      ncsTerm = forceField.getBoolean("NCSTERM", false);
      ncsTermOrig = ncsTerm;
    } else {
      ncsTerm = false;
      ncsTermOrig = false;
    }

    // Now that the crystal is set, check for special positions.
    int nSpecial = checkForSpecialPositions(forceField);

    // Flag to indicate bonded terms with hydrogen should be scaled up.
    boolean rigidHydrogens = forceField.getBoolean("RIGID_HYDROGENS", false);
    // Scale factor for increasing the strength of bonded terms involving hydrogen atoms.
    double rigidScale = forceField.getDouble("RIGID_SCALE", 10.0);

    boolean checkAllNodeCharges = forceField.getBoolean("CHECK_ALL_NODE_CHARGES", false);

    if (rigidScale <= 1.0) {
      rigidScale = 1.0;
    }

    logger.info("\n Bonded Terms");
    if (rigidHydrogens && rigidScale > 1.0) {
      logger.info(format("  Rigid hydrogens:                   %10.2f", rigidScale));
    }

    forceFieldBondedEnergyRegion = new EnergyTermRegion(parallelTeam, molecularAssembly, lambdaTerm);
    restrainEnergyRegion = new EnergyTermRegion(parallelTeam, molecularAssembly, lambdaTerm);

    // Collect bonds.
    if (bondTerm) {
      List<Bond> bondList = molecularAssembly.getBondList();
      if (nnTerm) {
        removeNeuralNetworkTerms(bondList);
      }
      if (!bondList.isEmpty()) {
        int forceGroup = forceField.getInteger("BOND_FORCE_GROUP", 0);
        bondPotentialEnergy = new BondPotentialEnergy("Bonds", forceGroup, bondList);
        forceFieldBondedEnergyRegion.addEnergyTerm(bondPotentialEnergy);
      } else {
        bondPotentialEnergy = null;
      }
    } else {
      bondPotentialEnergy = null;
    }

    // Collect angles.
    if (angleTerm) {
      List<Angle> angleList = molecularAssembly.getAngleList();
      if (nnTerm) {
        removeNeuralNetworkTerms(angleList);
      }
      if (!angleList.isEmpty()) {
        int forceGroup = forceField.getInteger("ANGLE_FORCE_GROUP", 0);
        anglePotentialEnergy = new AnglePotentialEnergy("Angles", forceGroup, angleList);
        forceFieldBondedEnergyRegion.addEnergyTerm(anglePotentialEnergy);
      } else {
        anglePotentialEnergy = null;
      }
    } else {
      anglePotentialEnergy = null;
    }

    // Collect stretch-bends.
    if (stretchBendTerm) {
      List<StretchBend> stretchBendList = molecularAssembly.getStretchBendList();
      if (nnTerm) {
        removeNeuralNetworkTerms(stretchBendList);
      }
      if (!stretchBendList.isEmpty()) {
        int forceGroup = forceField.getInteger("STRETCH_BEND_FORCE_GROUP", 0);
        stretchBendPotentialEnergy = new StretchBendPotentialEnergy("Stretch-Bends", forceGroup, stretchBendList);
        forceFieldBondedEnergyRegion.addEnergyTerm(stretchBendPotentialEnergy);
      } else {
        stretchBendPotentialEnergy = null;
      }
    } else {
      stretchBendPotentialEnergy = null;
    }

    // Collect Urey-Bradleys.
    if (ureyBradleyTerm) {
      List<UreyBradley> ureyBradleyList = molecularAssembly.getUreyBradleyList();
      if (nnTerm) {
        removeNeuralNetworkTerms(ureyBradleyList);
      }
      if (!ureyBradleyList.isEmpty()) {
        int forceGroup = forceField.getInteger("UREY_BRADLEY_FORCE_GROUP", 0);
        ureyBradleyPotentialEnergy = new UreyBradleyPotentialEnergy("Urey-Bradleys", forceGroup, ureyBradleyList);
        forceFieldBondedEnergyRegion.addEnergyTerm(ureyBradleyPotentialEnergy);
      } else {
        ureyBradleyPotentialEnergy = null;
      }
    } else {
      ureyBradleyPotentialEnergy = null;
    }

    // Set a multiplier on the force constants of bonded terms containing hydrogen atoms.
    if (rigidHydrogens) {
      if (bondPotentialEnergy != null) {
        for (Bond bond : bondPotentialEnergy.getBonds()) {
          if (bond.containsHydrogen()) {
            bond.setRigidScale(rigidScale);
          }
        }
      }
      if (anglePotentialEnergy != null) {
        for (Angle angle : anglePotentialEnergy.getAngles()) {
          if (angle.containsHydrogen()) {
            angle.setRigidScale(rigidScale);
          }
        }
      }
      if (stretchBendPotentialEnergy != null) {
        for (StretchBend stretchBend : stretchBendPotentialEnergy.getStretchBends()) {
          if (stretchBend.containsHydrogen()) {
            stretchBend.setRigidScale(rigidScale);
          }
        }
      }
      if (ureyBradleyPotentialEnergy != null) {
        for (UreyBradley ureyBradley : ureyBradleyPotentialEnergy.getUreyBradleys()) {
          if (ureyBradley.containsHydrogen()) {
            ureyBradley.setRigidScale(rigidScale);
          }
        }
      }
    }

    // Collect out-of-plane bends.
    if (outOfPlaneBendTerm) {
      List<OutOfPlaneBend> outOfPlaneBendList = molecularAssembly.getOutOfPlaneBendList();
      if (nnTerm) {
        removeNeuralNetworkTerms(outOfPlaneBendList);
      }
      if (!outOfPlaneBendList.isEmpty()) {
        int forceGroup = forceField.getInteger("OUT_OF_PLANE_BEND_FORCE_GROUP", 0);
        outOfPlaneBendPotentialEnergy = new OutOfPlaneBendPotentialEnergy("Out-0f-Plane Bends", forceGroup, outOfPlaneBendList);
        forceFieldBondedEnergyRegion.addEnergyTerm(outOfPlaneBendPotentialEnergy);
      } else {
        outOfPlaneBendPotentialEnergy = null;
      }
    } else {
      outOfPlaneBendPotentialEnergy = null;
    }

    // Collect torsions.
    if (torsionTerm) {
      List<Torsion> torsionList = molecularAssembly.getTorsionList();
      if (nnTerm) {
        removeNeuralNetworkTerms(torsionList);
      }
      if (!torsionList.isEmpty()) {
        int forceGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);
        torsionPotentialEnergy = new TorsionPotentialEnergy("Torsions", forceGroup, torsionList);
        forceFieldBondedEnergyRegion.addEnergyTerm(torsionPotentialEnergy);
      } else
        torsionPotentialEnergy = null;
    } else {
      torsionPotentialEnergy = null;
    }

    // Collect stretch-torsions.
    if (stretchTorsionTerm) {
      List<StretchTorsion> stretchTorsionList = molecularAssembly.getStretchTorsionList();
      if (nnTerm) {
        removeNeuralNetworkTerms(stretchTorsionList);
      }
      if (!stretchTorsionList.isEmpty()) {
        int forceGroup = forceField.getInteger("STRETCH_TORSION_FORCE_GROUP", 0);
        stretchTorsionPotentialEnergy = new StretchTorsionPotentialEnergy("Stretch-Torsions", forceGroup, stretchTorsionList);
        forceFieldBondedEnergyRegion.addEnergyTerm(stretchTorsionPotentialEnergy);
      } else {
        stretchTorsionPotentialEnergy = null;
      }
    } else {
      stretchTorsionPotentialEnergy = null;
    }

    // Collect angle-torsions.
    if (angleTorsionTerm) {
      List<AngleTorsion> angleTorsionList = molecularAssembly.getAngleTorsionList();
      if (nnTerm) {
        removeNeuralNetworkTerms(angleTorsionList);
      }
      if (!angleTorsionList.isEmpty()) {
        int forceGroup = forceField.getInteger("ANGLE_TORSION_FORCE_GROUP", 0);
        angleTorsionPotentialEnergy = new AngleTorsionPotentialEnergy("Angle-Torsions", forceGroup, angleTorsionList);
        forceFieldBondedEnergyRegion.addEnergyTerm(angleTorsionPotentialEnergy);
      } else {
        angleTorsionPotentialEnergy = null;
      }
    } else {
      angleTorsionPotentialEnergy = null;
    }

    // Collect improper torsions.
    if (improperTorsionTerm) {
      List<ImproperTorsion> improperTorsionList = molecularAssembly.getImproperTorsionList();
      if (nnTerm) {
        removeNeuralNetworkTerms(improperTorsionList);
      }
      if (!improperTorsionList.isEmpty()) {
        int forceGroup = forceField.getInteger("IMPROPER_TORSION_FORCE_GROUP", 0);
        improperTorsionPotentialEnergy = new ImproperTorsionPotentialEnergy("Improper Torsions", forceGroup, improperTorsionList);
        forceFieldBondedEnergyRegion.addEnergyTerm(improperTorsionPotentialEnergy);
      } else {
        improperTorsionPotentialEnergy = null;
      }
    } else {
      improperTorsionPotentialEnergy = null;
    }

    // Collect pi-orbital torsions.
    if (piOrbitalTorsionTerm) {
      List<PiOrbitalTorsion> piOrbitalTorsionList = molecularAssembly.getPiOrbitalTorsionList();
      if (nnTerm) {
        removeNeuralNetworkTerms(piOrbitalTorsionList);
      }
      if (!piOrbitalTorsionList.isEmpty()) {
        int forceGroup = forceField.getInteger("PI_ORBITAL_TORSION_FORCE_GROUP", 0);
        piOrbitalTorsionPotentialEnergy = new PiOrbitalTorsionPotentialEnergy("Pi-Orbital Torsions", forceGroup, piOrbitalTorsionList);
        forceFieldBondedEnergyRegion.addEnergyTerm(piOrbitalTorsionPotentialEnergy);
      } else {
        piOrbitalTorsionPotentialEnergy = null;
      }
    } else {
      piOrbitalTorsionPotentialEnergy = null;
    }

    // Collect torsion-torsions.
    if (torsionTorsionTerm) {
      List<TorsionTorsion> torsionTorsionList = molecularAssembly.getTorsionTorsionList();
      if (nnTerm) {
        removeNeuralNetworkTerms(torsionTorsionList);
      }
      if (!torsionTorsionList.isEmpty()) {
        int forceGroup = forceField.getInteger("TORSION_TORSION_FORCE_GROUP", 0);
        torsionTorsionPotentialEnergy = new TorsionTorsionPotentialEnergy("Torsion-Torsions", forceGroup, torsionTorsionList);
        forceFieldBondedEnergyRegion.addEnergyTerm(torsionTorsionPotentialEnergy);
      } else {
        torsionTorsionPotentialEnergy = null;
      }
    } else {
      torsionTorsionPotentialEnergy = null;
    }

    // Collect restrain positions.
    RestrainPosition[] restrainPositions = parseRestrainPositions(molecularAssembly);
    if (restrainPositions != null && restrainPositions.length > 0) {
      restrainPositionTerm = true;
    } else {
      restrainPositionTerm = false;
    }
    if (restrainPositionTerm) {
      int forceGroup = forceField.getInteger("RESTRAIN_POSITION_FORCE_GROUP", 0);
      restrainPositionPotentialEnergy = new RestrainPositionPotentialEnergy("Restrain Positions", forceGroup, java.util.Arrays.asList(restrainPositions));
      restrainEnergyRegion.addEnergyTerm(restrainPositionPotentialEnergy);
    } else {
      restrainPositionPotentialEnergy = null;
    }

    // Collect restrain distances.
    RestrainDistance[] restrainDistances = configureRestrainDistances(properties);
    if (restrainDistances != null && restrainDistances.length > 0) {
      restrainDistanceTerm = true;
    } else {
      restrainDistanceTerm = false;
    }
    if (restrainDistanceTerm) {
      int forceGroup = forceField.getInteger("RESTRAIN_DISTANCE_FORCE_GROUP", 0);
      restrainDistancePotentialEnergy = new RestrainDistancePotentialEnergy("Restrain Distances", forceGroup, java.util.Arrays.asList(restrainDistances));
      restrainEnergyRegion.addEnergyTerm(restrainDistancePotentialEnergy);
    } else {
      restrainDistancePotentialEnergy = null;
    }

    // Collect restrain torsions.
    Torsion[] restrainTorsions = configureRestrainTorsions(properties, forceField);
    if (restrainTorsions != null && restrainTorsions.length > 0) {
      restrainTorsionTerm = true;
    } else {
      restrainTorsionTerm = false;
    }
    if (restrainTorsionTerm) {
      logger.info(" Restrain Torsion Mode: " + restrainMode);
      int forceGroup = forceField.getInteger("RESTRAIN_TORSION_FORCE_GROUP", 0);
      restrainTorsionPotentialEnergy = new RestrainTorsionPotentialEnergy("Restrain Torsions", forceGroup, java.util.Arrays.asList(restrainTorsions));
      restrainEnergyRegion.addEnergyTerm(restrainTorsionPotentialEnergy);
    } else {
      restrainTorsionPotentialEnergy = null;
    }

    int[] molecule = molecularAssembly.getMoleculeNumbers();
    if (vanderWaalsTerm) {
      logger.info("\n Non-Bonded Terms");
      boolean[] nn = null;
      if (nnTerm) {
        nn = molecularAssembly.getNeuralNetworkIdentity();
      } else {
        nn = new boolean[nAtoms];
        Arrays.fill(nn, false);
      }

      if (!tornadoVM) {
        vanderWaals = new VanDerWaals(atoms, molecule, nn, crystal, forceField, parallelTeam,
            vdwCutoff, nlistCutoff);
      } else {
        vanderWaals = new VanDerWaalsTornado(atoms, crystal, forceField, vdwCutoff);
      }
    } else {
      vanderWaals = null;
    }

    if (nnTerm) {
      aniEnergy = new ANIEnergy(molecularAssembly);
    } else {
      aniEnergy = null;
    }

    if (multipoleTerm) {
      if (name.contains("OPLS") || name.contains("AMBER") || name.contains("CHARMM")) {
        elecForm = ELEC_FORM.FIXED_CHARGE;
      } else {
        elecForm = ELEC_FORM.PAM;
      }
      particleMeshEwald = new ParticleMeshEwald(atoms, molecule, forceField, crystal,
          vanderWaals.getNeighborList(), elecForm, ewaldCutoff, gkCutoff, parallelTeam);
      double charge = molecularAssembly.getCharge(checkAllNodeCharges);
      logger.info(format("\n  Overall system charge:             %10.3f", charge));
    } else {
      elecForm = null;
      particleMeshEwald = null;
    }

    if (ncsTerm) {
      String sg = forceField.getString("NCSGROUP", "P 1");
      Crystal ncsCrystal = new Crystal(a, b, c, alpha, beta, gamma, sg);
      ncsRestraint = new NCSRestraint(atoms, forceField, ncsCrystal);
    } else {
      ncsRestraint = null;
    }

    if (restrainGroupTerm) {
      restrainGroups = new RestrainGroups(molecularAssembly);
      nRestrainGroups = restrainGroups.getNumberOfRestraints();
    } else {
      restrainGroups = null;
    }

    if (comTerm) {
      comRestraint = new COMRestraint(molecularAssembly);
    } else {
      comRestraint = null;
    }

    // bondedRegion = new BondedRegion();
    maxDebugGradient = forceField.getDouble("MAX_DEBUG_GRADIENT", Double.POSITIVE_INFINITY);

    molecularAssembly.setPotential(this);

    // Configure constraints.
    constraints = configureConstraints(forceField);

    if (lambdaTerm) {
      this.setLambda(1.0);
      if (nSpecial == 0) {
        // For lambda calculations (e.g., polymorph searches), turn-off special position checks.
        // In this case, as the crystal lattice is alchemically turned on, special position overlaps are expected.
        crystal.setSpecialPositionCutoff(0.0);
      } else {
        // If special positions have already been identified, leave checks on.
        logger.info(" Special positions checking will be performed during a lambda simulation.");
      }
    }
  }

  /**
   * Get the MolecularAssembly associated with this ForceFieldEnergy.
   *
   * @return a {@link ffx.potential.MolecularAssembly} object.
   */
  public MolecularAssembly getMolecularAssembly() {
    return molecularAssembly;
  }

  /**
   * Set the lambdaTerm flag.
   *
   * @param lambdaTerm The value to set.
   */
  public void setLambdaTerm(boolean lambdaTerm) {
    this.lambdaTerm = lambdaTerm;
  }

  /**
   * Get the lambdaTerm flag.
   *
   * @return lambdaTerm.
   */
  public boolean getLambdaTerm() {
    return lambdaTerm;
  }

  private int checkForSpecialPositions(ForceField forceField) {
    // Check for atoms at special positions. These should normally be set to inactive.
    boolean specialPositionsInactive = forceField.getBoolean("SPECIAL_POSITIONS_INACTIVE", true);
    int nSpecial = 0;
    if (specialPositionsInactive) {
      int nSymm = crystal.getNumSymOps();
      if (nSymm > 1) {
        SpaceGroup spaceGroup = crystal.spaceGroup;
        double sp2 = crystal.getSpecialPositionCutoff2();
        double[] mate = new double[3];
        StringBuilder sb = new StringBuilder("\n Atoms at Special Positions set to Inactive:\n");
        for (int i = 0; i < nAtoms; i++) {
          Atom atom = atoms[i];
          double[] xyz = atom.getXYZ(null);
          for (int iSymm = 1; iSymm < nSymm; iSymm++) {
            SymOp symOp = spaceGroup.getSymOp(iSymm);
            crystal.applySymOp(xyz, mate, symOp);
            double dr2 = crystal.image(xyz, mate);
            if (dr2 < sp2) {
              sb.append(
                  format("  %s separation with SymOp %d at %8.6f A.\n", atoms[i], iSymm, sqrt(dr2)));
              atom.setActive(false);
              nSpecial++;
              break;
            }
          }
        }
        if (nSpecial > 0) {
          logger.info(sb.toString());
        }
      }
    }
    return nSpecial;
  }

  /**
   * Method to parse the restrain-torsion-cos records.
   *
   * @param properties Configuration properties.
   * @param forceField Force field properties.
   * @return An array of restrain torsions, or null if none were found.
   */
  private Torsion[] configureRestrainTorsions(CompositeConfiguration properties, ForceField forceField) {
    StringBuilder restrainLog = new StringBuilder("\n  Restrain-Torsions");

    String[] restrainTorsions = properties.getStringArray("restrain-torsion");
    double torsionUnits = forceField.getDouble("TORSIONUNIT", TorsionType.DEFAULT_TORSION_UNIT);
    List<Torsion> restrainTorsionList = new ArrayList<>(restrainTorsions.length);
    for (String restrainString : restrainTorsions) {
      // Add the key back to the input line.
      restrainString = "restrain-torsion " + restrainString;
      // Split the line on the pound symbol to remove comments.
      String input = restrainString.split("#+")[0];
      // Split the line on whitespace.
      String[] tokens = input.trim().split(" +");
      // Restrain torsion records have a similar form as torsion records.
      // The first four tokens are atom indices instead of atom classes.
      TorsionType torsionType = TorsionType.parse(input, tokens);
      torsionType.torsionUnit = torsionUnits;

      // Collect the atom indices.
      int[] atomIndices = torsionType.atomClasses;
      int ai1 = atomIndices[0] - 1;
      int ai2 = atomIndices[1] - 1;
      int ai3 = atomIndices[2] - 1;
      int ai4 = atomIndices[3] - 1;
      Atom a1 = atoms[ai1];
      Atom a2 = atoms[ai2];
      Atom a3 = atoms[ai3];
      Atom a4 = atoms[ai4];

      // Collect the bonds between the atoms making up the restrain torsion.
      Bond firstBond = a1.getBond(a2);
      Bond middleBond = a2.getBond(a3);
      Bond lastBond = a3.getBond(a4);
      Torsion torsion = new Torsion(firstBond, middleBond, lastBond);
      torsion.setTorsionType(torsionType);
      restrainTorsionList.add(torsion);
      restrainLog.append("\n   ").append(torsion);
    }

    if (!restrainTorsionList.isEmpty()) {
      logger.info(restrainLog.toString());
      return restrainTorsionList.toArray(new Torsion[0]);
    } else {
      return null;
    }
  }

  /**
   * Method to parse the restrain-distance records.
   *
   * @param properties Configuration properties.
   */
  private RestrainDistance[] configureRestrainDistances(CompositeConfiguration properties) {
    List<RestrainDistance> restrainDistanceList = new ArrayList<>();
    String[] bondRestraints = properties.getStringArray("restrain-distance");
    for (String bondRest : bondRestraints) {
      try {
        String[] toks = bondRest.split("\\s+");
        if (toks.length < 2) {
          throw new IllegalArgumentException(
              format(" restrain-distance value %s could not be parsed!", bondRest));
        }
        // Internally, everything starts with 0, but restrain distance starts at 1, so that 1 has to
        // be subtracted
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
          case 2:
          case 3:
            double[] xyz1 = new double[3];
            xyz1 = a1.getXYZ(xyz1);
            double[] xyz2 = new double[3];
            xyz2 = a2.getXYZ(xyz2);
            // Current distance between restrained atoms
            dist = crystal.minDistOverSymOps(xyz1, xyz2);
            break;
          case 4:
            dist = Double.parseDouble(toks[3]);
            break;
          case 5:
          default:
            double minDist = Double.parseDouble(toks[3]);
            double maxDist = Double.parseDouble(toks[4]);
            dist = 0.5 * (minDist + maxDist);
            flatBottomRadius = 0.5 * Math.abs(maxDist - minDist);
            break;
        }

        UnivariateSwitchingFunction switchF;
        double lamStart = RestrainDistance.DEFAULT_RB_LAM_START;
        double lamEnd = RestrainDistance.DEFAULT_RB_LAM_END;
        if (toks.length > 5) {
          int offset = 5;
          if (toks[5].matches("^[01](?:\\.[0-9]*)?")) {
            offset = 6;
            lamStart = Double.parseDouble(toks[5]);
            if (toks[6].matches("^[01](?:\\.[0-9]*)?")) {
              offset = 7;
              lamEnd = Double.parseDouble(toks[6]);
            }
          }
          switchF = UnivariateFunctionFactory.parseUSF(toks, offset);
        } else {
          switchF = new ConstantSwitch();
        }

        RestrainDistance restrainDistance = createRestrainDistance(a1, a2, dist,
            forceConst, flatBottomRadius, lamStart, lamEnd, switchF);
        restrainDistanceList.add(restrainDistance);
      } catch (Exception ex) {
        logger.info(format(" Exception in parsing restrain-distance: %s", ex));
      }
    }
    if (!restrainDistanceList.isEmpty()) {
      return restrainDistanceList.toArray(new RestrainDistance[0]);
    } else {
      return null;
    }
  }

  /**
   * Configure bonded constraints.
   *
   * @param forceField Force field properties.
   * @return List of constraints.
   */
  private List<Constraint> configureConstraints(ForceField forceField) {
    String constraintStrings = forceField.getString("CONSTRAIN", forceField.getString("RATTLE", null));
    if (constraintStrings == null) {
      return Collections.emptyList();
    }

    ArrayList<Constraint> constraints = new ArrayList<>();
    logger.info(format(" Experimental: parsing constraints option %s", constraintStrings));

    Set<Bond> numericBonds = new HashSet<>(1);
    Set<Angle> numericAngles = new HashSet<>(1);

    // Totally empty constrain option: constrain only X-H bonds. No other options applied.
    if (constraintStrings.isEmpty() || constraintStrings.matches("^\\s*$")) {
      // Assume constraining only X-H bonds (i.e. RIGID-HYDROGEN).
      logger.info(" Constraining X-H bonds.");
      Bond[] bonds = bondPotentialEnergy.getBondArray();
      numericBonds = Arrays.stream(bonds).filter(
          (Bond bond) -> bond.getAtom(0).getAtomicNumber() == 1
              || bond.getAtom(1).getAtomicNumber() == 1).collect(Collectors.toSet());
    } else {
      String[] constraintToks = constraintStrings.split("\\s+");

      // First, accumulate SETTLE constraints.
      for (String tok : constraintToks) {
        if (tok.equalsIgnoreCase("WATER")) {
          logger.info(" Constraining waters to be rigid based on angle & bonds.");
          // XYZ files, in particular, have waters mislabeled as generic Molecules.
          // First, find any such mislabeled water.
          Stream<MSNode> settleStream = molecularAssembly.getMolecules().stream()
              .filter((MSNode m) -> m.getAtomList().size() == 3).filter((MSNode m) -> {
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
          settleStream = Stream.concat(settleStream, molecularAssembly.getWater().stream());
          // Map them into new Settle constraints and collect.
          List<SettleConstraint> settleConstraints = settleStream.map(
              (MSNode m) -> m.getAngleList().getFirst()).map(SettleConstraint::settleFactory).toList();
          constraints.addAll(settleConstraints);

        } else if (tok.equalsIgnoreCase("DIATOMIC")) {
          logger.severe(" Diatomic distance constraints not yet implemented properly.");
        } else if (tok.equalsIgnoreCase("TRIATOMIC")) {
          logger.severe(
              " Triatomic SETTLE constraints for non-water molecules not yet implemented properly.");
        }
      }

      // Second, accumulate bond/angle constraints.
      for (String tok : constraintToks) {
        if (tok.equalsIgnoreCase("BONDS") && bondPotentialEnergy != null) {
          Bond[] bonds = bondPotentialEnergy.getBondArray();
          numericBonds = new HashSet<>(Arrays.asList(bonds));
        } else if (tok.equalsIgnoreCase("ANGLES") && anglePotentialEnergy != null) {
          Angle[] angles = anglePotentialEnergy.getAngleArray();
          numericAngles = new HashSet<>(Arrays.asList(angles));
        }
      }
    }

    // Remove bonds that are already dealt with via angles.
    for (Angle angle : numericAngles) {
      angle.getBondList().forEach(numericBonds::remove);
    }

    // Remove already-constrained angles and bonds (e.g. SETTLE-constrained ones).
    List<Angle> ccmaAngles = numericAngles.stream().filter((Angle ang) -> !ang.isConstrained()).collect(Collectors.toList());
    List<Bond> ccmaBonds = numericBonds.stream().filter((Bond bond) -> !bond.isConstrained()).collect(Collectors.toList());

    CcmaConstraint ccmaConstraint = CcmaConstraint.ccmaFactory(ccmaBonds, ccmaAngles, atoms,
        getMass(), CcmaConstraint.DEFAULT_CCMA_NONZERO_CUTOFF);
    constraints.add(ccmaConstraint);
    logger.info(format(" Added %d constraints.", constraints.size()));

    return constraints;
  }

  /**
   * Static factory method to create a ForceFieldEnergy, possibly via FFX or OpenMM implementations.
   *
   * @param assembly To create FFE over
   * @return a {@link ffx.potential.ForceFieldEnergy} object.
   */
  public static ForceFieldEnergy energyFactory(MolecularAssembly assembly) {
    return energyFactory(assembly, ParallelTeam.getDefaultThreadCount());
  }

  /**
   * Static factory method to create a ForceFieldEnergy, possibly via FFX or OpenMM implementations.
   *
   * @param assembly   To create FFE over
   * @param numThreads Number of threads to use for FFX energy
   * @return A ForceFieldEnergy on some Platform
   */
  @SuppressWarnings("fallthrough")
  public static ForceFieldEnergy energyFactory(MolecularAssembly assembly, int numThreads) {
    ForceField forceField = assembly.getForceField();
    String platformString = toEnumForm(forceField.getString("PLATFORM", "FFX"));
    Platform platform;
    try {
      platform = Platform.valueOf(platformString);
    } catch (IllegalArgumentException e) {
      logger.warning(format(" String %s did not match a known energy implementation", platformString));
      platform = Platform.FFX;
    }

    // Check if the dual-topology platform flag is set.
    String dtPlatformString = toEnumForm(forceField.getString("PLATFORM_DT", "FFX"));
    try {
      Platform dtPlatform = Platform.valueOf(dtPlatformString);
      // If the dtPlatform uses OpenMM, then single topology energy will use FFX.
      if (dtPlatform != Platform.FFX) {
        // If the dtPlatform platform uses OpenMM, then use FFX for single topology.
        platform = Platform.FFX;
      }
    } catch (IllegalArgumentException e) {
      // If the dtPlatform is not recognized, ignore it.
    }

    switch (platform) {
      case OMM, OMM_REF, OMM_CUDA, OMM_OPENCL:
        try {
          return new OpenMMEnergy(assembly, platform, numThreads);
        } catch (Exception ex) {
          logger.warning(format(" Exception creating OpenMMEnergy: %s", ex));
          ForceFieldEnergy ffxEnergy = assembly.getPotentialEnergy();
          if (ffxEnergy == null) {
            ffxEnergy = new ForceFieldEnergy(assembly, numThreads);
            assembly.setPotential(ffxEnergy);
          }
          return ffxEnergy;
        }
      case OMM_CPU:
        logger.warning(format(" Platform %s not supported; defaulting to FFX", platform));
      default:
        ForceFieldEnergy ffxEnergy = assembly.getPotentialEnergy();
        if (ffxEnergy == null) {
          ffxEnergy = new ForceFieldEnergy(assembly, numThreads);
          assembly.setPotential(ffxEnergy);
        }
        return ffxEnergy;
    }
  }

  /**
   * Get all atoms that make up this ForceFieldEnergy.
   *
   * @return An array of Atoms.
   */
  public Atom[] getAtomArray() {
    return atoms;
  }

  /**
   * Applies constraints to positions
   *
   * @param xPrior Prior coodinates.
   * @param xNew   New coordinates.
   */
  public void applyAllConstraintPositions(double[] xPrior, double[] xNew) {
    applyAllConstraintPositions(xPrior, xNew, DEFAULT_CONSTRAINT_TOLERANCE);
  }

  /**
   * Applies constraints to positions
   *
   * @param xPrior Prior coodinates.
   * @param xNew   New coordinates.
   * @param tol    Constraint Tolerance.
   */
  public void applyAllConstraintPositions(double[] xPrior, double[] xNew, double tol) {
    if (xPrior == null) {
      xPrior = Arrays.copyOf(xNew, xNew.length);
    }
    for (Constraint constraint : constraints) {
      constraint.applyConstraintToStep(xPrior, xNew, getMass(), tol);
    }
  }

  /**
   * Overwrites current esvSystem if present. Multiple ExtendedSystems is possible but unnecessary;
   * add all ESVs to one system (per FFE, at least).
   *
   * @param system a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public void attachExtendedSystem(ExtendedSystem system) {
    if (system == null) {
      throw new IllegalArgumentException();
    }
    esvTerm = true;
    esvTermOrig = esvTerm;
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
      } else {
        particleMeshEwald.attachExtendedSystem(system);
      }
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns true if lambda term is not enabled for this ForceFieldEnergy.
   */
  @Override
  public boolean dEdLZeroAtEnds() {
    // This may actually be true even with softcored atoms.
    // For now, serves the purpose of reporting true when nothing is softcored.
    return !lambdaTerm;
  }

  /**
   * Frees up assets associated with this ForceFieldEnergy, such as worker Threads.
   *
   * @return If successful in freeing up assets.
   */
  public boolean destroy() {
    if (destroyed) {
      // This regularly occurs with Repex OST, as multiple OrthogonalSpaceTempering objects wrap a
      // single FFE.
      logger.fine(format(" This ForceFieldEnergy is already destroyed: %s", this));
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
   * energy.
   *
   * @return a double.
   */
  public double energy() {
    return energy(false, false);
  }

  /**
   * Compute the potential energy of the system.
   *
   * @param gradient If true, compute the Cartesian coordinate gradient.
   * @param print    If true, print the energy terms.
   * @return the energy in kcal/mol.
   */
  public double energy(boolean gradient, boolean print) {
    try {
      totalTime = System.nanoTime();
      nnTime = 0;
      restrainGroupTime = 0;
      vanDerWaalsTime = 0;
      electrostaticTime = 0;
      ncsTime = 0;

      // Zero out the potential energy of each bonded term.
      double forceFieldBondedEnergy = 0.0;
      nnEnergy = 0.0;

      // Zero out potential energy of restraint terms
      double restrainEnergy = 0.0;
      restrainGroupEnergy = 0.0;
      ncsEnergy = 0.0;

      // Zero out the potential energy of each non-bonded term.
      double totalNonBondedEnergy = 0.0;
      double totalMultipoleEnergy = 0.0;
      vanDerWaalsEnergy = 0.0;
      permanentMultipoleEnergy = 0.0;
      polarizationEnergy = 0.0;

      // Zero out the solvation energy.
      solvationEnergy = 0.0;
      esvBias = 0.0;

      // Zero out the total potential energy.
      totalEnergy = 0.0;

      // BondedEnergyRegion to compute bonded terms in parallel.
      try {
        forceFieldBondedEnergyRegion.setGradient(gradient);
        forceFieldBondedEnergyRegion.setLambdaBondedTerms(lambdaBondedTerms);
        forceFieldBondedEnergyRegion.setLambdaAllBondedTerms(lambdaAllBondedTerms);
        forceFieldBondedEnergyRegion.setInitAtomGradients(true);
        parallelTeam.execute(forceFieldBondedEnergyRegion);
        forceFieldBondedEnergy = forceFieldBondedEnergyRegion.getEnergy();
      } catch (RuntimeException ex) {
        logger.warning("Runtime exception during bonded term calculation.");
        throw ex;
      } catch (Exception ex) {
        logger.info(Utilities.stackTraceToString(ex));
        logger.severe(ex.toString());
      }

      // Restrain Energy Region to compute restrain terms in parallel.
      try {
        if ((restrainMode == RestrainMode.ENERGY) ||
            (restrainMode == RestrainMode.ALCHEMICAL && lambdaBondedTerms)) {
          boolean checkAlchemicalAtoms = (restrainMode == RestrainMode.ENERGY);
          restrainEnergyRegion.setGradient(gradient);
          restrainEnergyRegion.setLambdaBondedTerms(lambdaBondedTerms);
          restrainEnergyRegion.setLambdaAllBondedTerms(lambdaAllBondedTerms);
          restrainEnergyRegion.setCheckAlchemicalAtoms(checkAlchemicalAtoms);
          restrainEnergyRegion.setInitAtomGradients(false);
          parallelTeam.execute(restrainEnergyRegion);
          restrainEnergy = restrainEnergyRegion.getEnergy();
        }
      } catch (RuntimeException ex) {
        logger.warning("Runtime exception during restrain term calculation.");
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

        /*
        if (restrainPositionTerm) {
          restrainPositionTime = -System.nanoTime();
          for (RestrainPosition restrainPosition : restrainPositions) {
            restrainPositionEnergy += restrainPosition.residual(gradient);
          }
          restrainPositionTime += System.nanoTime();
        }
        */

        if (restrainGroupTerm) {
          restrainGroupTime = -System.nanoTime();
          restrainGroupEnergy = restrainGroups.energy(gradient);
          restrainGroupTime += System.nanoTime();
        }

        if (comTerm) {
          comRestraintTime = -System.nanoTime();
          comRestraintEnergy = comRestraint.residual(gradient, print);
          comRestraintTime += System.nanoTime();
        }


        // Compute the neural network term.
        if (nnTerm) {
          nnTime = -System.nanoTime();
          nnEnergy = aniEnergy.energy(gradient, print);
          nnTime += System.nanoTime();
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
          polarizationEnergy = particleMeshEwald.getPolarizationEnergy();
          nPermanentInteractions = particleMeshEwald.getInteractions();
          solvationEnergy = particleMeshEwald.getSolvationEnergy();
          nGKInteractions = particleMeshEwald.getGKInteractions();
          electrostaticTime += System.nanoTime();
        }
      }

      totalTime = System.nanoTime() - totalTime;

      double totalBondedEnergy = forceFieldBondedEnergy + restrainEnergy
          + nnEnergy + ncsEnergy + comRestraintEnergy + restrainGroupEnergy;

      totalNonBondedEnergy = vanDerWaalsEnergy + totalMultipoleEnergy;
      totalEnergy = totalBondedEnergy + totalNonBondedEnergy + solvationEnergy;
      if (esvTerm) {
        esvBias = esvSystem.getBiasEnergy();
        totalEnergy += esvBias;
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
    } catch (EnergyException ex) {
      if (printOnFailure) {
        printFailure();
      }
      if (ex.doCauseSevere()) {
        logger.info(Utilities.stackTraceToString(ex));
        logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
      } else {
        logger.log(Level.INFO, format(" Exception in energy calculation:\n %s", ex));
      }
      throw ex;
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
        boolean extremeGrad = Arrays.stream(g)
            .anyMatch((double gi) -> (gi > maxDebugGradient || gi < -maxDebugGradient));
        if (extremeGrad) {
          File origFile = molecularAssembly.getFile();
          String timeString = LocalDateTime.now()
              .format(DateTimeFormatter.ofPattern("yyyy_MM_dd-HH_mm_ss"));

          String filename = format("%s-LARGEGRAD-%s.pdb",
              removeExtension(molecularAssembly.getFile().getName()), timeString);
          PotentialsFunctions ef = new PotentialsUtils();
          filename = ef.versionFile(filename);

          logger.warning(format(" Excessively large gradient detected; printing snapshot to file %s",
              filename));
          ef.saveAsPDB(molecularAssembly, new File(filename));
          molecularAssembly.setFile(origFile);
        }
      }
      return e;
    } catch (EnergyException ex) {
      if (printOnFailure) {
        printFailure();
      }
      if (ex.doCauseSevere()) {
        logger.info(Utilities.stackTraceToString(ex));
        logger.log(Level.SEVERE, " Error in calculating energies or gradients", ex);
      } else {
        logger.log(Level.INFO, format(" Exception in energy calculation:\n %s", ex));
      }
      throw ex;
    }
  }

  /**
   * Save coordinates when an EnergyException is caught.
   */
  private void printFailure() {
    String timeString = LocalDateTime.now()
        .format(DateTimeFormatter.ofPattern("yyyy_MM_dd-HH_mm_ss"));
    File file = molecularAssembly.getFile();
    String ext = "pdb";
    if (isXYZ(file)) {
      ext = "xyz";
    }
    String filename = format("%s-ERROR-%s.%s", removeExtension(file.getName()), timeString, ext);
    PotentialsFunctions ef = new PotentialsUtils();
    filename = ef.versionFile(filename);
    logger.info(format(" Writing on-error snapshot to file %s", filename));
    ef.save(molecularAssembly, new File(filename));
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns an array of acceleration values for active atoms.
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
   * Getter for the field <code>angles</code> with only <code>AngleMode</code> angles.
   *
   * @param angleMode Only angles of this mode will be returned.
   * @return an array of {@link ffx.potential.bonded.Angle} objects.
   */
  public Angle[] getAngles(AngleMode angleMode) {
    if (anglePotentialEnergy == null) {
      return null;
    }
    Angle[] angles = anglePotentialEnergy.getAngleArray();
    int nAngles = angles.length;
    List<Angle> angleList = new ArrayList<>();
    // Sort all normal angles from in-plane angles
    for (int i = 0; i < nAngles; i++) {
      if (angles[i].getAngleMode() == angleMode) {
        angleList.add(angles[i]);
      }
    }
    nAngles = angleList.size();
    if (nAngles < 1) {
      return null;
    }
    return angleList.toArray(new Angle[0]);
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
      if (a.isActive()) {
        x[index++] = a.getX();
        x[index++] = a.getY();
        x[index++] = a.getZ();
      }
    }
    return x;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Getter for the field <code>crystal</code>.
   */
  @Override
  public Crystal getCrystal() {
    return crystal;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Set the boundary conditions for this calculation.
   */
  @Override
  public void setCrystal(Crystal crystal) {
    setCrystal(crystal, false);
  }

  /**
   * Getter for the field <code>cutoffPlusBuffer</code>.
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
  public STATE getEnergyTermState() {
    return state;
  }

  /**
   * {@inheritDoc}
   *
   * <p>This method is for the RESPA integrator only.
   */
  @Override
  public void setEnergyTermState(STATE state) {
    this.state = state;
    if (forceFieldBondedEnergyRegion != null) {
      forceFieldBondedEnergyRegion.setState(state);
    }
    if (restrainEnergyRegion != null) {
      restrainEnergyRegion.setState(state);
    }
    switch (state) {
      case FAST:
        nnTerm = nnTermOrig;
        restrainGroupTerm = restrainGroupTermOrig;
        ncsTerm = ncsTermOrig;
        comTerm = comTermOrig;
        vanderWaalsTerm = false;
        multipoleTerm = false;
        polarizationTerm = false;
        generalizedKirkwoodTerm = false;
        esvTerm = false;
        break;
      case SLOW:
        vanderWaalsTerm = vanderWaalsTermOrig;
        multipoleTerm = multipoleTermOrig;
        polarizationTerm = polarizationTermOrig;
        generalizedKirkwoodTerm = generalizedKirkwoodTermOrig;
        esvTerm = esvTermOrig;
        nnTerm = false;
        restrainGroupTerm = false;
        ncsTerm = false;
        comTerm = false;
        break;
      default:
        nnTerm = nnTermOrig;
        restrainGroupTerm = restrainGroupTermOrig;
        ncsTerm = ncsTermOrig;
        comTermOrig = comTerm;
        vanderWaalsTerm = vanderWaalsTermOrig;
        multipoleTerm = multipoleTermOrig;
        polarizationTerm = polarizationTermOrig;
        generalizedKirkwoodTerm = generalizedKirkwoodTermOrig;
    }
  }

  /**
   * getEsvBiasEnergy.
   *
   * @return a double.
   */
  public double getEsvBiasEnergy() {
    return esvBias;
  }

  /**
   * getExtendedSystem.
   *
   * @return a {@link ffx.potential.extended.ExtendedSystem} object.
   */
  public ExtendedSystem getExtendedSystem() {
    return esvSystem;
  }

  /**
   * getGK.
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
   * Returns the gradient array for this ForceFieldEnergy.
   *
   * @param g an array of double.
   * @return the gradient array.
   */
  public double[] getGradient(double[] g) {
    return fillGradient(g);
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

        /*
        if (restrainPositionTerm && nRestrainPositions > 0) {
          for (RestrainPosition restrainPosition : restrainPositions) {
            restrainPosition.setLambda(lambda);
          }
        } */

        if (restrainDistanceTerm && restrainDistancePotentialEnergy != null) {
          for (RestrainDistance restrainDistance : restrainDistancePotentialEnergy.getRestrainDistances()) {
            restrainDistance.setLambda(lambda);
          }
        }
        if (restrainTorsionTerm && restrainTorsionPotentialEnergy != null) {
          for (Torsion restrainTorsion : restrainTorsionPotentialEnergy.getRestrainTorsions()) {
            restrainTorsion.setLambda(lambda);
          }
        }
        if (ncsTerm && ncsRestraint != null) {
          ncsRestraint.setLambda(lambda);
        }
        if (comTerm && comRestraint != null) {
          comRestraint.setLambda(lambda);
        }
        if (lambdaTorsions) {
          if (torsionPotentialEnergy != null) {
            torsionPotentialEnergy.setLambda(lambda);
          }
          if (piOrbitalTorsionPotentialEnergy != null) {
            piOrbitalTorsionPotentialEnergy.setLambda(lambda);
          }
          if (torsionTorsionPotentialEnergy != null) {
            torsionTorsionPotentialEnergy.setLambda(lambda);
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
      if (atoms[i].isActive()) {
        nActive++;
      }
    }
    return nActive * 3;
  }

  /**
   * Create a PDB REMARK 3 string containing the potential energy terms.
   *
   * @return a String containing the PDB REMARK 3 formatted energy terms.
   */
  public String getPDBHeaderString() {
    energy(false, false);
    StringBuilder sb = new StringBuilder();
    sb.append("REMARK   3  CALCULATED POTENTIAL ENERGY\n");
    sb.append(forceFieldBondedEnergyRegion.toPDBString());
    sb.append(restrainEnergyRegion.toPDBString());
    if (nnTerm) {
      sb.append(
          format("REMARK   3   %s %g (%d)\n", "NEUTRAL NETWORK            : ", nnEnergy, nAtoms));
    }
    if (ncsTerm) {
      sb.append(
          format("REMARK   3   %s %g (%d)\n", "NCS RESTRAINT              : ", ncsEnergy, nAtoms));
    }
    if (comTerm) {
      sb.append(
          format("REMARK   3   %s %g (%d)\n", "COM RESTRAINT              : ", comRestraintEnergy,
              nAtoms));
    }
    if (vanderWaalsTerm) {
      sb.append(
          format("REMARK   3   %s %g (%d)\n", "VAN DER WAALS              : ", vanDerWaalsEnergy,
              nVanDerWaalInteractions));
    }
    if (multipoleTerm) {
      sb.append(format("REMARK   3   %s %g (%d)\n", "ATOMIC MULTIPOLES          : ",
          permanentMultipoleEnergy, nPermanentInteractions));
    }
    if (polarizationTerm) {
      sb.append(format("REMARK   3   %s %g (%d)\n", "POLARIZATION               : ",
          polarizationEnergy, nPermanentInteractions));
    }
    sb.append(format("REMARK   3   %s %g\n", "TOTAL POTENTIAL (KCAL/MOL) : ", totalEnergy));
    int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
    if (nsymm > 1) {
      sb.append(format("REMARK   3   %s %g\n", "UNIT CELL POTENTIAL        : ", totalEnergy * nsymm));
    }
    if (crystal.getUnitCell() != crystal) {
      nsymm = crystal.spaceGroup.getNumberOfSymOps();
      if (nsymm > 1) {
        sb.append(format("REMARK   3   %s %g\n", "REPLICATES CELL POTENTIAL  : ", totalEnergy * nsymm));
      }
    }
    sb.append("REMARK   3\n");

    return sb.toString();
  }

  /**
   * Getter for the field <code>parallelTeam</code>.
   *
   * @return a {@link edu.rit.pj.ParallelTeam} object.
   */
  public ParallelTeam getParallelTeam() {
    return parallelTeam;
  }

  /**
   * getPermanentInteractions.
   *
   * @return a int.
   */
  public int getPermanentInteractions() {
    return nPermanentInteractions;
  }

  /**
   * Getter for the field <code>permanentMultipoleEnergy</code>.
   *
   * @return a double.
   */
  public double getPermanentMultipoleEnergy() {
    return permanentMultipoleEnergy;
  }

  /**
   * Gets the Platform associated with this force field energy. For the reference platform, always
   * returns FFX.
   *
   * @return A Platform.
   */
  public Platform getPlatform() {
    return platform;
  }

  /**
   * getPmeNode.
   *
   * @return a {@link ParticleMeshEwald} object.
   */
  public ParticleMeshEwald getPmeNode() {
    return particleMeshEwald;
  }

  /**
   * Getter for the field <code>polarizationEnergy</code>.
   *
   * @return a double.
   */
  public double getPolarizationEnergy() {
    return polarizationEnergy;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns an array of previous acceleration values for active atoms.
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
   * {@inheritDoc}
   */
  @Override
  public double[] getScaling() {
    return optimizationScaling;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public void setScaling(double[] scaling) {
    optimizationScaling = scaling;
  }

  /**
   * Getter for the field <code>solvationEnergy</code>.
   *
   * @return a double.
   */
  public double getSolvationEnergy() {
    return solvationEnergy;
  }

  /**
   * getSolvationInteractions.
   *
   * @return a int.
   */
  public int getSolvationInteractions() {
    return nGKInteractions;
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public double getTotalEnergy() {
    return totalEnergy;
  }

  /**
   * Getter for the field <code>vanDerWaalsEnergy</code>.
   *
   * @return a double.
   */
  public double getVanDerWaalsEnergy() {
    return vanDerWaalsEnergy;
  }

  /**
   * getVanDerWaalsInteractions.
   *
   * @return a int.
   */
  public int getVanDerWaalsInteractions() {
    return nVanDerWaalInteractions;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Return a reference to each variables type.
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
   * getVdwNode.
   *
   * @return a {@link ffx.potential.nonbonded.VanDerWaals} object.
   */
  public VanDerWaals getVdwNode() {
    return vanderWaals;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Returns an array of velocity values for active atoms.
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
      /*
      if (restrainPositionTerm && nRestrainPositions > 0) {
        for (RestrainPosition restrainPosition : restrainPositions) {
          d2EdLambda2 += restrainPosition.getd2EdL2();
        }
      }
      */
      if (restrainDistanceTerm && restrainDistancePotentialEnergy != null) {
        for (RestrainDistance restrainDistance : restrainDistancePotentialEnergy.getRestrainDistanceArray()) {
          d2EdLambda2 += restrainDistance.getd2EdL2();
        }
      }
      if (ncsTerm && ncsRestraint != null) {
        d2EdLambda2 += ncsRestraint.getd2EdL2();
      }
      if (comTerm && comRestraint != null) {
        d2EdLambda2 += comRestraint.getd2EdL2();
      }
      if (lambdaTorsions) {
        if (torsionPotentialEnergy != null) {
          d2EdLambda2 += torsionPotentialEnergy.getd2EdL2();
        }
        if (piOrbitalTorsionPotentialEnergy != null) {
          d2EdLambda2 += piOrbitalTorsionPotentialEnergy.getd2EdL2();
        }
        if (torsionTorsionPotentialEnergy != null) {
          d2EdLambda2 += torsionTorsionPotentialEnergy.getd2EdL2();
        }
      }
    }
    return d2EdLambda2;
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
      if (restrainDistanceTerm && restrainDistancePotentialEnergy != null) {
        for (RestrainDistance restrainDistance : restrainDistancePotentialEnergy.getRestrainDistanceArray()) {
          dEdLambda += restrainDistance.getdEdL();
        }
      }
      /*
      if (restrainPositionTerm && nRestrainPositions > 0) {
        for (RestrainPosition restrainPosition : restrainPositions) {
          dEdLambda += restrainPosition.getdEdL();
        }
      }
      */
      if (ncsTerm && ncsRestraint != null) {
        dEdLambda += ncsRestraint.getdEdL();
      }
      if (comTerm && comRestraint != null) {
        dEdLambda += comRestraint.getdEdL();
      }
      if (lambdaTorsions) {
        if (torsionPotentialEnergy != null) {
          dEdLambda += torsionPotentialEnergy.getdEdL();
        }
        if (piOrbitalTorsionPotentialEnergy != null) {
          dEdLambda += piOrbitalTorsionPotentialEnergy.getdEdL();
        }
        if (torsionTorsionPotentialEnergy != null) {
          dEdLambda += torsionTorsionPotentialEnergy.getdEdL();
        }
      }
    }
    return dEdLambda;
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
      if (restrainDistanceTerm && restrainDistancePotentialEnergy != null) {
        for (RestrainDistance restrainDistance : restrainDistancePotentialEnergy.getRestrainDistanceArray()) {
          restrainDistance.getdEdXdL(gradients);
        }
      }
      /*
      if (restrainPositionTerm && nRestrainPositions > 0) {
        for (RestrainPosition restrainPosition : restrainPositions) {
          restrainPosition.getdEdXdL(gradients);
        }
      }
      */
      if (ncsTerm && ncsRestraint != null) {
        ncsRestraint.getdEdXdL(gradients);
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
   *
   * <p>The acceleration array should only contain acceleration data for active atoms.
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
   * The coordinate array should only contain active atoms.
   *
   * @param coords the coordinates to set.
   */
  public void setCoordinates(@Nullable double[] coords) {
    if (coords == null) {
      return;
    }
    int index = 0;
    for (Atom a : atoms) {
      if (a.isActive()) {
        double x = coords[index++];
        double y = coords[index++];
        double z = coords[index++];
        a.moveTo(x, y, z);
      }
    }
  }

  /**
   * Set the boundary conditions for this calculation.
   *
   * @param crystal             Crystal to set.
   * @param checkReplicatesCell Check if a replicates cell must be created.
   */
  public void setCrystal(Crystal crystal, boolean checkReplicatesCell) {
    if (checkReplicatesCell) {
      this.crystal = ReplicatesCrystal.replicatesCrystalFactory(crystal.getUnitCell(), cutOff2);
    } else {
      this.crystal = crystal;
    }
    /*
     Update VanDerWaals first, in case the NeighborList needs to be
     re-allocated to include a larger number of replicated cells.
    */
    if (vanderWaalsTerm) {
      vanderWaals.setCrystal(this.crystal);
    }
    if (multipoleTerm) {
      particleMeshEwald.setCrystal(this.crystal);
    }
  }

  /**
   * {@inheritDoc}
   *
   * <p>The previousAcceleration array should only contain previous acceleration data for active
   * atoms.
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
   * Sets the printOnFailure flag; if override is true, over-rides any existing property. Essentially
   * sets the default value of printOnFailure for an algorithm. For example, rotamer optimization
   * will generally run into force field issues in the normal course of execution as it tries
   * unphysical self and pair configurations, so the algorithm should not print out a large number of
   * error PDBs.
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
        molecularAssembly.getForceField().getBoolean("PRINT_ON_FAILURE");
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
   * {@inheritDoc}
   *
   * <p>The velocity array should only contain velocity data for active atoms.
   */
  @Override
  public void setVelocity(double[] velocity) {
    if (velocity == null) {
      return;
    }
    int index = 0;
    double[] vel = new double[3];
    for (Atom a : atoms) {
      if (a.isActive()) {
        vel[0] = velocity[index++];
        vel[1] = velocity[index++];
        vel[2] = velocity[index++];
      } else {
        // If the atom is not active, set the velocity to zero.
        vel[0] = 0.0;
        vel[1] = 0.0;
        vel[2] = 0.0;
      }
      a.setVelocity(vel);
    }
  }

  /**
   * {@inheritDoc}
   */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append(forceFieldBondedEnergyRegion.toString());
    sb.append(restrainEnergyRegion.toString());
    if (restrainGroupTerm && nRestrainGroups > 0) {
      sb.append(format("  %s %20.8f %12d %12.3f\n", "Restrain Groups   ", restrainGroupEnergy,
          nRestrainGroups, restrainGroupTime * toSeconds));
    }
    if (ncsTerm) {
      sb.append(format("  %s %20.8f %12d %12.3f\n", "NCS Restraint     ", ncsEnergy, nAtoms,
          ncsTime * toSeconds));
    }
    if (comTerm) {
      sb.append(format("  %s %20.8f %12d %12.3f\n", "COM Restraint     ", comRestraintEnergy, nAtoms,
          comRestraintTime * toSeconds));
    }
    if (vanderWaalsTerm && nVanDerWaalInteractions > 0) {
      sb.append(format("  %s %20.8f %12d %12.3f\n", "Van der Waals     ", vanDerWaalsEnergy,
          nVanDerWaalInteractions, vanDerWaalsTime * toSeconds));
    }
    if (multipoleTerm && nPermanentInteractions > 0) {
      if (polarizationTerm) {
        sb.append(format("  %s %20.8f %12d\n", "Atomic Multipoles ", permanentMultipoleEnergy,
            nPermanentInteractions));
      } else {
        if (elecForm == ELEC_FORM.FIXED_CHARGE) {
          sb.append(format("  %s %20.8f %12d %12.3f\n", "Atomic Charges    ", permanentMultipoleEnergy,
              nPermanentInteractions, electrostaticTime * toSeconds));
        } else
          sb.append(format("  %s %20.8f %12d %12.3f\n", "Atomic Multipoles ", permanentMultipoleEnergy,
              nPermanentInteractions, electrostaticTime * toSeconds));
      }
    }
    if (polarizationTerm && nPermanentInteractions > 0) {
      sb.append(format("  %s %20.8f %12d %12.3f\n", "Polarization      ", polarizationEnergy,
          nPermanentInteractions, electrostaticTime * toSeconds));
    }
    if (generalizedKirkwoodTerm && nGKInteractions > 0) {
      sb.append(
          format("  %s %20.8f %12d\n", "Solvation         ", solvationEnergy, nGKInteractions));
    }
    if (esvTerm) {
      sb.append(
          format("  %s %20.8f  %s\n", "ExtendedSystemBias", esvBias, esvSystem.getLambdaList()));
      sb.append(esvSystem.getBiasDecomposition());
    }
    if (nnTerm) {
      sb.append(format("  %s %20.8f %12d %12.3f\n", "Neural Network    ", nnEnergy, nAtoms,
          nnTime * toSeconds));
    }
    sb.append(format("  %s %20.8f  %s %12.3f (sec)",
        "Total Potential   ", totalEnergy, "(Kcal/mole)", totalTime * toSeconds));
    int nsymm = crystal.getUnitCell().spaceGroup.getNumberOfSymOps();
    if (nsymm > 1) {
      sb.append(format("\n  %s %20.8f", "Unit Cell         ", totalEnergy * nsymm));
    }
    if (crystal.getUnitCell() != crystal) {
      nsymm = crystal.spaceGroup.getNumberOfSymOps();
      sb.append(format("\n  %s %20.8f", "Replicates Cell   ", totalEnergy * nsymm));
    }

    return sb.toString();
  }

  /**
   * Log bonded energy terms and restraints.
   */
  public void logBondedTermsAndRestraints() {
    if (forceFieldBondedEnergyRegion != null) {
      forceFieldBondedEnergyRegion.log();
    }
    if (restrainEnergyRegion != null) {
      restrainEnergyRegion.log();
    }
    if (restrainGroupTerm && nRestrainGroups > 0) {
      logger.info("\n Restrain Group Interactions:");
      logger.info(restrainGroups.toString());
    }
  }

  /**
   * Setter for the field <code>lambdaBondedTerms</code>.
   *
   * @param lambdaBondedTerms    a boolean.
   * @param lambdaAllBondedTerms a boolean.
   */
  void setLambdaBondedTerms(boolean lambdaBondedTerms, boolean lambdaAllBondedTerms) {
    this.lambdaBondedTerms = lambdaBondedTerms;
    this.lambdaAllBondedTerms = lambdaAllBondedTerms;
  }

  /**
   * Private method for internal use, so we don't have subclasses calling super.energy, and this
   * class delegating to the subclass's getGradient method.
   *
   * @param g Gradient array to fill.
   * @return Gradient array.
   */
  private double[] fillGradient(double[] g) {
    assert (g != null);
    double[] grad = new double[3];
    int n = getNumberOfVariables();
    if (g == null || g.length < n) {
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
        if (isNaN(gx) || isInfinite(gx) || isNaN(gy) || isInfinite(gy) || isNaN(gz) || isInfinite(
            gz)) {
          StringBuilder sb = new StringBuilder(
              format("The gradient of atom %s is (%8.3f,%8.3f,%8.3f).", a, gx, gy, gz));
          double[] vals = new double[3];
          a.getVelocity(vals);
          sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          a.getAcceleration(vals);
          sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
          a.getPreviousAcceleration(vals);
          sb.append(
              format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));

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
   * setRestrainDistance
   *
   * @param a1                a {@link ffx.potential.bonded.Atom} object.
   * @param a2                a {@link ffx.potential.bonded.Atom} object.
   * @param distance          a double.
   * @param forceConstant     the force constant in kcal/mole.
   * @param flatBottom        Radius of a flat-bottom potential in Angstroms.
   * @param lamStart          At what lambda does the restraint begin to take effect?
   * @param lamEnd            At what lambda does the restraint hit full strength?
   * @param switchingFunction Switching function to use as a lambda dependence.
   */
  private RestrainDistance createRestrainDistance(Atom a1, Atom a2, double distance, double forceConstant,
                                                  double flatBottom, double lamStart, double lamEnd,
                                                  UnivariateSwitchingFunction switchingFunction) {
    boolean rbLambda = !(switchingFunction instanceof ConstantSwitch) && lambdaTerm;
    RestrainDistance restrainDistance = new RestrainDistance(a1, a2, crystal, rbLambda, lamStart, lamEnd, switchingFunction);
    int[] classes = {a1.getAtomType().atomClass, a2.getAtomType().atomClass};
    if (flatBottom != 0) {
      BondType bondType = new BondType(classes, forceConstant, distance,
          BondType.BondFunction.FLAT_BOTTOM_HARMONIC, flatBottom);
      restrainDistance.setBondType(bondType);
    } else {
      BondType bondType = new BondType(classes, forceConstant, distance,
          BondType.BondFunction.HARMONIC);
      restrainDistance.setBondType(bondType);
    }
    return restrainDistance;
  }

  private Crystal configureNCS(ForceField forceField, Crystal unitCell) {
    // MTRIXn to be permuted with standard space group in NCSCrystal.java for experimental
    // refinement.
    if (forceField.getProperties().containsKey("MTRIX1")
        && forceField.getProperties().containsKey("MTRIX2") && forceField.getProperties()
        .containsKey("MTRIX3")) {
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
          logger.info(
              format(" MTRIXn SymOp: %d of %d\n" + symOp, i + 1, MTRX1List.length));
        }
        spaceGroup.symOps.add(symOp);
      }
      unitCell = NCSCrystal.NCSCrystalFactory(unitCell, spaceGroup.symOps);
      unitCell.updateCrystal();
    }
    return unitCell;
  }

  public RestrainMode getRestrainMode() {
    return restrainMode;
  }

  /**
   * Getter for the RestrainGroup field.
   *
   * @return Returns RestrainGroup.
   */
  public RestrainGroups getRestrainGroups() {
    return restrainGroups;
  }

  /**
   * Get the TorsionTorsionPotentialEnergy.
   *
   * @return The TorsionTorsionPotentialEnergy.
   */
  public TorsionTorsionPotentialEnergy getTorsionTorsionPotentialEnergy() {
    return torsionTorsionPotentialEnergy;
  }

  /**
   * Get the AnglePotentialEnergy.
   *
   * @return The AnglePotentialEnergy.
   */
  public AnglePotentialEnergy getAnglePotentialEnergy() {
    return anglePotentialEnergy;
  }

  /**
   * Get the BondPotentialEnergy.
   *
   * @return The BondPotentialEnergy.
   */
  public BondPotentialEnergy getBondPotentialEnergy() {
    return bondPotentialEnergy;
  }

  /**
   * Get the StretchBendPotentialEnergy.
   *
   * @return The StretchBendPotentialEnergy.
   */
  public StretchBendPotentialEnergy getStretchBendPotentialEnergy() {
    return stretchBendPotentialEnergy;
  }

  /**
   * Get the UreyBradleyPotentialEnergy.
   *
   * @return The UreyBradleyPotentialEnergy.
   */
  public UreyBradleyPotentialEnergy getUreyBradleyPotentialEnergy() {
    return ureyBradleyPotentialEnergy;
  }

  /**
   * Get the OutOfPlaneBendPotentialEnergy.
   *
   * @return The OutOfPlaneBendPotentialEnergy.
   */
  public OutOfPlaneBendPotentialEnergy getOutOfPlaneBendPotentialEnergy() {
    return outOfPlaneBendPotentialEnergy;
  }

  /**
   * Get the TorsionPotentialEnergy.
   *
   * @return The TorsionPotentialEnergy.
   */
  public TorsionPotentialEnergy getTorsionPotentialEnergy() {
    return torsionPotentialEnergy;
  }

  /**
   * Get the StretchTorsionPotentialEnergy.
   *
   * @return The StretchTorsionPotentialEnergy.
   */
  public StretchTorsionPotentialEnergy getStretchTorsionPotentialEnergy() {
    return stretchTorsionPotentialEnergy;
  }

  /**
   * Get the AngleTorsionPotentialEnergy.
   *
   * @return The AngleTorsionPotentialEnergy.
   */
  public AngleTorsionPotentialEnergy getAngleTorsionPotentialEnergy() {
    return angleTorsionPotentialEnergy;
  }

  /**
   * Get the ImproperTorsionPotentialEnergy.
   *
   * @return The ImproperTorsionPotentialEnergy.
   */
  public ImproperTorsionPotentialEnergy getImproperTorsionPotentialEnergy() {
    return improperTorsionPotentialEnergy;
  }

  /**
   * Get the PiOrbitalTorsionPotentialEnergy.
   *
   * @return The PiOrbitalTorsionPotentialEnergy.
   */
  public PiOrbitalTorsionPotentialEnergy getPiOrbitalTorsionPotentialEnergy() {
    return piOrbitalTorsionPotentialEnergy;
  }

  /**
   * Get the RestrainPositionPotentialEnergy.
   *
   * @return The RestrainPositionPotentialEnergy.
   */
  public RestrainPositionPotentialEnergy getRestrainPositionPotentialEnergy() {
    return restrainPositionPotentialEnergy;
  }

  /**
   * Get the RestrainDistancePotentialEnergy.
   *
   * @return The RestrainDistancePotentialEnergy.
   */
  public RestrainDistancePotentialEnergy getRestrainDistancePotentialEnergy() {
    return restrainDistancePotentialEnergy;
  }

  /**
   * Get the RestrainTorsionPotentialEnergy.
   *
   * @return The RestrainTorsionPotentialEnergy.
   */
  public RestrainTorsionPotentialEnergy getRestrainTorsionPotentialEnergy() {
    return restrainTorsionPotentialEnergy;
  }

}
