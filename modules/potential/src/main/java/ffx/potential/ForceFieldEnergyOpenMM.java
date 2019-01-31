/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import static java.lang.String.format;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.sqrt;

import edu.rit.mp.CharacterBuf;
import edu.rit.pj.Comm;
import edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod;
import edu.uiowa.jopenmm.OpenMMLibrary.*;
import edu.uiowa.jopenmm.OpenMMUtils;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_3D_DoubleArray_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_3D_DoubleArray_destroy;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_3D_DoubleArray_set;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaAngleForce_addAngle;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaAngleForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaBondForce_addBond;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaBondForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaInPlaneAngleForce_addAngle;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaInPlaneAngleForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent12;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent13;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent14;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent15;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent11;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_Bisector;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_NoAxisType;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ThreeFold;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZBisect;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZOnly;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZThenX;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_NoCutoff;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_PME;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Direct;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Extrapolated;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Mutual;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_addMultipole;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setAEwald;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setCovalentMap;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setCutoffDistance;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setMultipoleParameters;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setPmeGridDimensions;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_setPolarizationType;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaMultipoleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaOutOfPlaneBendForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaPiTorsionForce_addPiTorsion;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaPiTorsionForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaStretchBendForce_addStretchBend;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaStretchBendForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaTorsionTorsionForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_addParticle;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setCutoffDistance;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setParticleExclusions;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setParticleParameters;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setUseDispersionCorrection;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_updateParametersInContext;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_addParticle;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_create;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setAwater;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setDispoff;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setEpsh;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setEpso;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setParticleParameters;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setRminh;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setRmino;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setShctd;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_setSlevy;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaWcaDispersionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AngstromsPerNm;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_RadiansPerDegree;
import static edu.uiowa.jopenmm.OpenMMLibrary.*;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePair;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePairNoExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_SingleParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_ParameterDerivatives;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.StretchTorsion;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.AngleType.AngleFunction;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.BondType.BondFunction;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldBoolean;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ImproperTorsionType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.PiTorsionType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWType;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import static ffx.potential.nonbonded.VanDerWaalsForm.EPSILON_RULE.GEOMETRIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_RULE.ARITHMETIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_SIZE.RADIUS;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_TYPE.R_MIN;
import static ffx.potential.nonbonded.VanDerWaalsForm.VDW_TYPE.LENNARD_JONES;

/**
 * Compute the potential energy and derivatives using OpenMM.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("deprecation")
public class ForceFieldEnergyOpenMM extends ForceFieldEnergy {

    private static final Logger logger = Logger.getLogger(ForceFieldEnergyOpenMM.class.getName());

    /**
     * Requested Platform (i.e. Java or an OpenMM platform).
     */
    private final Platform ffxPlatform;
    /**
     * OpenMM Platform.
     */
    private PointerByReference platform = null;
    /**
     * OpenMM System.
     */
    private PointerByReference system = null;
    /**
     * OpenMM Context.
     */
    private PointerByReference context = null;
    /**
     * OpenMM Integrator.
     */
    private PointerByReference integrator = null;
    /**
     * Integrator string (default = VERLET).
     */
    private String integratorString = "VERLET";
    /**
     * Time step (default = 1.0 fsec).
     */
    private double timeStep = 1.0;
    /**
     * Temperature (default = 298.15).
     */
    private double temperature = 298.15;
    /**
     * OpenMM thermostat. Currently an Andersen thermostat is supported.
     */
    private PointerByReference ommThermostat = null;
    /**
     * Barostat to be added if NPT (isothermal-isobaric) dynamics is requested.
     */
    private PointerByReference ommBarostat = null;
    /**
     * OpenMM State.
     */
    private PointerByReference state;
    /**
     * OpenMM Forces.
     */
    private PointerByReference forces;
    /**
     * OpenMM Positions.
     */
    private PointerByReference positions = null;
    /**
     * OpenMM Velocities.
     */
    private PointerByReference velocities = null;
    /**
     * OpenMM center-of-mass motion remover.
     */
    private PointerByReference commRemover = null;
    /**
     * Number of particles.
     */
    private int numParticles = 0;
    /**
     * OpenMM AMOEBA Bond Force.
     */
    private PointerByReference amoebaBondForce = null;
    /**
     * OpenMM AMOEBA Angle Force.
     */
    private PointerByReference amoebaAngleForce = null;
    /**
     * OpenMM AMOEBA In-Plane Angle Force.
     */
    private PointerByReference amoebaInPlaneAngleForce = null;
    /**
     * OpenMM AMOEBA Urey Bradley Force.
     */
    private PointerByReference amoebaUreyBradleyForce = null;
    /**
     * OpenMM AMOEBA Out-of-Plane Bend Force.
     */
    private PointerByReference amoebaOutOfPlaneBendForce = null;
    /**
     * OpenMM AMOEBA Stretch Bend Force.
     */
    private PointerByReference amoebaStretchBendForce = null;
    /**
     * OpenMM AMOEBA Torsion Force.
     */
    private PointerByReference amoebaTorsionForce = null;
    /**
     * OpenMM Improper Torsion Force.
     */
    private PointerByReference amoebaImproperTorsionForce = null;
    /**
     * OpenMM AMOEBA Pi Torsion Force.
     */
    private PointerByReference amoebaPiTorsionForce = null;
    /**
     * OpenMM AMOEBA Torsion Torsion Force.
     */
    private PointerByReference amoebaTorsionTorsionForce = null;
    /**
     * OpenMM AMOEBA van der Waals Force.
     */
    private PointerByReference amoebaVDWForce = null;
    /**
     * OpenMM AMOEBA Multipole Force.
     */
    private PointerByReference amoebaMultipoleForce = null;
    /**
     * OpenMM Generalized Kirkwood Force.
     */
    private PointerByReference amoebaGeneralizedKirkwoodForce = null;
    /**
     * OpenMM AMOEBA WCA Dispersion Force.
     */
    private PointerByReference amoebaWcaDispersionForce = null;
    /**
     * OpenMM Fixed Charge Non-Bonded Force.
     */
    private PointerByReference fixedChargeNonBondedForce = null;
    /**
     * OpenMM Stretch-Torsion couplings (as found in phosphate/nucleic acid
     * AMOEBA force fields).
     */
    private PointerByReference stretchTorsionForce = null;
    /**
     * OpenMM Angle-Torsion couplings (as found in phosphate/nucleic acid AMOEBA
     * force fields).
     */
    private PointerByReference angleTorsionForce = null;
    /**
     * Fixed charge softcore vdW force boolean.
     */
    boolean softcoreCreated = false;
    /**
     * Fixed charge softcore force.
     */
    private PointerByReference fixedChargeSoftcore = null;
    /**
     * Sterics force between alchemical atoms.
     */
    private PointerByReference alchemicalAlchemicalStericsForce = null;
    /**
     * Sterics force between alchemical and non alchemical atoms.
     */
    private PointerByReference nonAlchemicalAlchemicalStericsForce = null;
    /**
     * Boolean array, holds charge exclusion list.
     */
    private boolean[] chargeExclusion;
    /**
     * Boolean array, holds van Der Waals exclusion list.
     */
    private boolean[] vdWExclusion;
    /**
     * Double array, holds charge quantity value for exceptions.
     */
    private double[] exceptionChargeProd;
    /**
     * Double array, holds epsilon quantity value for exceptions.
     */
    private double[] exceptionEps;
    /**
     * OpenMM Custom GB Force.
     */
    private PointerByReference customGBForce = null;
    /**
     * Map from bond functional forms to the restraint-bonds using that
     * functional form. The PointerByReference should point to a
     * CustomBondForce.
     */
    private final Map<BondType.BondFunction, PointerByReference> restraintForces = new HashMap<>();
    /**
     * Langevin friction coefficient.
     */
    private double frictionCoeff;
    /**
     * Andersen thermostat collision frequency.
     */
    private double collisionFreq;
    /**
     * Lambda flag to indicate control of electrostatic scaling. If both elec and vdW are being scaled, then vdW
     * is scaled first, followed by elec.
     */
    private boolean elecLambdaTerm;
    /**
     * Lambda flag to indicate control of vdW scaling. If both elec and vdW are being scaled, then vdW
     * is scaled first, followed by elec.
     */
    private boolean vdwLambdaTerm;
    /**
     * Lambda flag to indicate control of torsional force constants (L=0 corresponds to torsions being off, and L=1 to
     * torsions at full strength.
     */
    private boolean torsionLambdaTerm;
    /**
     * Value of the lambda state variable.
     */
    private double lambda = 1.0;
    /**
     * Value of the van der Waals lambda state variable.
     */
    private double lambdaVDW = 1.0;
    /**
     * Value of the electrostatics lambda state variable.
     */
    private double lambdaElec = 1.0;
    /**
     * Value of the electrostatics lambda state variable.
     */
    private double lambdaTorsion = 1.0;
    /**
     * Derivative of van Der Waals contribution to the potential energy with respect to lambda.
     */
    private double vdwdUdL = 0.0;
    /**
     * Lambda step size for finite difference dU/dL.
     */
    private double fdDLambda;
    /**
     * Flag to set water molecule bonds as rigid.
     */
    private boolean rigidHydrogen;

    /**
     * Whether to enforce periodic boundary conditions when obtaining new
     * States.
     */
    public final int enforcePBC;


    /**
     * Boolean to control logging statements to the screen, typically used when an integrator other than the default is chosen for dynamics.
     */
    private boolean quiet = true;
    /**
     * Integer that controls the value of the quiet boolean.
     */
    private int quietInt = 0;

    // private boolean doOpenMMdEdL = false;
    // private boolean doFFXdEdL = true;
    private boolean testdEdL = true;

    /**
     * ForceFieldEnergyOpenMM constructor; offloads heavy-duty computation to an
     * OpenMM Platform while keeping track of information locally.
     *
     * @param molecularAssembly Assembly to construct energy for.
     * @param requestedPlatform requested OpenMM platform to be used.
     * @param restraints        Harmonic coordinate restraints.
     * @param nThreads          Number of threads to use in the super class
     *                          ForceFieldEnergy instance.
     */
    protected ForceFieldEnergyOpenMM(MolecularAssembly molecularAssembly, Platform requestedPlatform,
                                     List<CoordRestraint> restraints, int nThreads) {
        super(molecularAssembly, restraints, nThreads);

        Crystal crystal = getCrystal();
        int symOps = crystal.spaceGroup.getNumberOfSymOps();
        if (symOps > 1) {
            logger.info("");
            logger.severe(" OpenMM does not support symmetry operators.");
        }

        //super.energy(false, true);
        logger.info("\n Initializing OpenMM");

        loadPlatform(requestedPlatform);
        ffxPlatform = requestedPlatform;

        // Create the OpenMM System
        system = OpenMM_System_create();
        logger.info(" System created.");

        ForceField forceField = molecularAssembly.getForceField();

        // Load atoms.
        try {
            addAtoms();
        } catch (Exception e) {
            logger.severe(" Atom without mass encountered.");
        }

        rigidHydrogen = forceField.getBoolean(ForceField.ForceFieldBoolean.RIGID_HYDROGEN, false);

        if (rigidHydrogen) {
            setUpHydrogenConstraints(system);
        }

        // Add Bond Force.
        addBondForce(forceField);

        // Add Angle Force.
        addAngleForce(forceField);
        addInPlaneAngleForce(forceField);

        // Add Stretch-Bend Force.
        addStretchBendForce(forceField);

        // Add Urey-Bradley Force.
        addUreyBradleyForce(forceField);

        // Out-of Plane Bend Force.
        addOutOfPlaneBendForce(forceField);

        // Add Torsion Force.
        addTorsionForce(forceField);

        // Add Improper Torsion Force.
        addImproperTorsionForce(forceField);

        // Add Pi-Torsion Force.
        addPiTorsionForce(forceField);

        // Add Torsion-Torsion Force.
        addTorsionTorsionForce(forceField);

        // Add coordinate restraints.
        addHarmonicRestraintForce(forceField);

        // Add bond restraints.
        addRestraintBonds(forceField);

        // Add stretch-torsion coupling terms.
        addStretchTorsionCoupling(forceField);

        // Add angle-torsion coupling terms.
        addAngleTorsionCoupling(forceField);

        VanDerWaals vdW = super.getVdwNode();
        if (vdW != null) {
            VanDerWaalsForm vdwForm = vdW.getVDWForm();
            if (vdwForm.vdwType == LENNARD_JONES) {
                addFixedChargeNonBondedForce(forceField);
            } else {
                // Add vdW Force.
                addAmoebaVDWForce(forceField);

                // Add Multipole Force.
                addAmoebaMultipoleForce(forceField);
            }
        }

        // Set periodic box vectors.
        setDefaultPeriodicBoxVectors();

        frictionCoeff = forceField.getDouble(ForceFieldDouble.FRICTION_COEFF, 91.0);
        collisionFreq = forceField.getDouble(ForceFieldDouble.COLLISION_FREQ, 0.01);

        createContext(integratorString, timeStep, temperature);

        // Set initial positions.
        double x[] = new double[numParticles * 3];
        int index = 0;
        Atom atoms[] = molecularAssembly.getAtomArray();
        for (int i = 0; i < numParticles; i++) {
            Atom atom = atoms[i];
            x[index] = atom.getX();
            x[index + 1] = atom.getY();
            x[index + 2] = atom.getZ();
            index += 3;
        }
        setOpenMMPositions(x, numParticles * 3);

        int infoMask = OpenMM_State_Positions;
        infoMask += OpenMM_State_Forces;
        infoMask += OpenMM_State_Energy;

        boolean aperiodic = super.getCrystal().aperiodic();
        boolean pbcEnforced = forceField.getBoolean(ForceField.ForceFieldBoolean.ENFORCE_PBC, !aperiodic);
        enforcePBC = pbcEnforced ? OpenMM_True : OpenMM_False;

        state = OpenMM_Context_getState(context, infoMask, enforcePBC);
        forces = OpenMM_State_getForces(state);
        double openMMPotentialEnergy = OpenMM_State_getPotentialEnergy(state) / OpenMM_KJPerKcal;

        logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", openMMPotentialEnergy));
        fdDLambda = forceField.getDouble(ForceFieldDouble.FD_DLAMBDA, 0.001);

        OpenMM_State_destroy(state);

        elecLambdaTerm = forceField.getBoolean(ForceFieldBoolean.ELEC_LAMBDATERM,false);
        vdwLambdaTerm = forceField.getBoolean(ForceFieldBoolean.VDW_LAMBDATERM, false);
        torsionLambdaTerm = forceField.getBoolean(ForceFieldBoolean.TORSION_LAMBDATERM, false);

        if (elecLambdaTerm || vdwLambdaTerm || torsionLambdaTerm) {
            lambdaTerm = true;
        } else {
            lambdaTerm = false;
        }

        if (lambdaTerm) {
            logger.info(format(" Lambda scales vdW interactions: %s", vdwLambdaTerm));
            logger.info(format(" Lambda scales electrostatics:   %s", elecLambdaTerm));
            logger.info(format(" Lambda scales torsions:   %s", torsionLambdaTerm));
        }

//        CompositeConfiguration properties = molecularAssembly.getProperties();
//        if (properties.containsKey("openMMdEdL")) {
//            doOpenMMdEdL = true;
//            doFFXdEdL = false;
//        }

    }

    /**
     * Load an OpenMM Platform
     */
    private void loadPlatform(Platform requestedPlatform) {

        OpenMMUtils.init();

        // Print out the OpenMM Version.
        Pointer version = OpenMM_Platform_getOpenMMVersion();
        logger.log(Level.INFO, " Version: {0}", version.getString(0));

        // Print out the OpenMM plugin directory.
        logger.log(Level.INFO, " Plugin Dir: {0}", OpenMMUtils.OPENMM_PLUGIN_DIR);

        /**
         * Load plugins and print out plugins.
         *
         * Call the method twice to avoid a bug where not all platforms are list
         * after the first call.
         */
        PointerByReference plugins = OpenMM_Platform_loadPluginsFromDirectory(OpenMMUtils.OPENMM_PLUGIN_DIR);
        OpenMM_StringArray_destroy(plugins);

        plugins = OpenMM_Platform_loadPluginsFromDirectory(OpenMMUtils.OPENMM_PLUGIN_DIR);
        int numPlugins = OpenMM_StringArray_getSize(plugins);

        logger.log(Level.INFO, " Number of Plugins: {0}", numPlugins);
        boolean cuda = false;
        for (int i = 0; i < numPlugins; i++) {
            String pluginString = stringFromArray(plugins, i);
            logger.log(Level.INFO, "  Plugin: {0}", pluginString);
            if (pluginString.toUpperCase().contains("AMOEBACUDA")) {
                cuda = true;
            }
        }
        OpenMM_StringArray_destroy(plugins);

        /**
         * Extra logging to print out plugins that failed to load.
         */
        if (logger.isLoggable(Level.FINE)) {
            PointerByReference pluginFailers = OpenMM_Platform_getPluginLoadFailures();
            int numFailures = OpenMM_StringArray_getSize(pluginFailers);
            for (int i = 0; i < numFailures; i++) {
                String pluginString = stringFromArray(pluginFailers, i);
                logger.log(Level.FINE, " Plugin load failure: {0}", pluginString);
            }
            OpenMM_StringArray_destroy(pluginFailers);
        }

        int numPlatforms = OpenMM_Platform_getNumPlatforms();
        logger.log(Level.INFO, " Number of Platforms: {0}", numPlatforms);

        String defaultPrecision = "mixed";
        String precision = molecularAssembly.getForceField().getString(ForceField.ForceFieldString.PRECISION, defaultPrecision).toLowerCase();
        precision = precision.replace("-precision", "");
        switch (precision) {
            case "double":
            case "mixed":
            case "single":
                logger.info(String.format(" Precision level: %s", precision));
                break;
            default:
                logger.info(String.format(" Could not interpret precision level %s, defaulting to %s", precision, defaultPrecision));
                precision = defaultPrecision;
                break;
        }

        if (cuda && requestedPlatform != Platform.OMM_REF) {
            int defaultDevice = getDefaultDevice(molecularAssembly.getProperties());
            platform = OpenMM_Platform_getPlatformByName("CUDA");
            int deviceID = molecularAssembly.getForceField().getInteger(ForceField.ForceFieldInteger.CUDA_DEVICE, defaultDevice);
            String deviceIDString = Integer.toString(deviceID);

            OpenMM_Platform_setPropertyDefaultValue(platform, pointerForString("CudaDeviceIndex"), pointerForString(deviceIDString));
            OpenMM_Platform_setPropertyDefaultValue(platform, pointerForString("Precision"), pointerForString(precision));
            logger.info(String.format(" Platform: AMOEBA CUDA (Device ID %d)", deviceID));
            try {
                Comm world = Comm.world();
                if (world != null) {
                    logger.info(String.format(" Running on host %s, rank %d", world.host(), world.rank()));
                }
            } catch (IllegalStateException ise) {
                logger.fine(" Could not find the world communicator!");
            }
        } else {
            platform = OpenMM_Platform_getPlatformByName("Reference");
            logger.info(" Platform: AMOEBA CPU Reference");
        }
    }

    /**
     * Gets the default co-processor device, ignoring any CUDA_DEVICE over-ride.
     * This is either determined by process rank and the
     * availableDevices/CUDA_DEVICES property, or just 0 if neither property is
     * sets.
     *
     * @param props Properties in use.
     * @return Pre-override device index.
     */
    private static int getDefaultDevice(CompositeConfiguration props) {
        String availDeviceProp = props.getString("availableDevices",
                props.getString("CUDA_DEVICES"));
        if (availDeviceProp == null) {
            int nDevs = props.getInt("numCudaDevices", 1);
            availDeviceProp = IntStream.range(0, nDevs).
                    mapToObj(Integer::toString).
                    collect(Collectors.joining(" "));
        }
        availDeviceProp = availDeviceProp.trim();
        String[] availDevices = availDeviceProp.split("\\s+");
        int nDevs = availDevices.length;
        int[] devs = new int[nDevs];
        for (int i = 0; i < nDevs; i++) {
            devs[i] = Integer.parseInt(availDevices[i]);
        }
        int index = 0;
        try {
            Comm world = Comm.world();
            if (world != null) {
                int size = world.size();

                // Format the host as a CharacterBuf of length 100.
                int messageLen = 100;
                String host = world.host();
                // Truncate to max 100 characters.
                host = host.substring(0, Math.min(messageLen, host.length()));
                // Pad to 100 characters.
                host = String.format("%-100s", host);
                char[] messageOut = host.toCharArray();
                CharacterBuf out = CharacterBuf.buffer(messageOut);

                // Now create CharacterBuf array for all incoming messages.
                char[][] incoming = new char[size][messageLen];
                CharacterBuf[] in = new CharacterBuf[size];
                for (int i = 0; i < size; i++) {
                    in[i] = CharacterBuf.buffer(incoming[i]);
                }

                try {
                    world.allGather(out, in);
                } catch (IOException ex) {
                    logger.severe(String.format(" Failure at the allGather step for determining rank: %s\n%s", ex, Utilities.stackTraceToString(ex)));
                }
                int ownIndex = -1;
                int rank = world.rank();
                boolean selfFound = false;

                for (int i = 0; i < size; i++) {
                    String hostI = new String(incoming[i]);
                    if (hostI.equalsIgnoreCase(host)) {
                        ++ownIndex;
                        if (i == rank) {
                            selfFound = true;
                            break;
                        }
                    }
                }
                if (!selfFound) {
                    logger.severe(String.format(" Rank %d: Could not find any incoming host messages matching self %s!", rank, host.trim()));
                } else {
                    index = ownIndex % nDevs;
                }
            }
        } catch (IllegalStateException ise) {
            // Behavior is just to keep index = 0.
        }
        return devs[index];
    }

    /**
     * createContext takes in a parameters to determine which integrator the
     * user requested during the start up of the simulation. A switch statement
     * is used with Strings as the variable to determine between Lengevin,
     * Brownian, Custom, Compound and Verlet integrator
     *
     * @param integratorString a {@link java.lang.String} object.
     * @param timeStep         a double.
     * @param temperature      a double.
     */
    public void createContext(String integratorString, double timeStep, double temperature) {

        this.integratorString = integratorString;
        this.timeStep = timeStep;
        this.temperature = temperature;

        CompositeConfiguration properties = molecularAssembly.getProperties();

        //logger.info(String.format(" quietInt is %d", quietInt));
        if (quietInt == 1) {
            quiet = false;
        }

        OpenMM_Context_destroy(context);

        double dt = timeStep * 1.0e-3;
        switch (integratorString) {
            case "LANGEVIN":
                if (!quiet) {
                    logger.log(Level.INFO, String.format(" Created Langevin integrator with time step %6.3e (psec).", dt));
                }
                integrator = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, dt);

                if (properties.containsKey("randomseed")) {
                    int randomSeed = properties.getInt("randomseed", 0);
                    logger.info(String.format(" Setting random seed %d for Langevin dynamics", randomSeed));
                    OpenMM_LangevinIntegrator_setRandomNumberSeed(integrator, randomSeed);
                }
                break;
            case "RESPA":
                // Read in the inner time step in fsec, then convert to psec.
                int in = molecularAssembly.getProperties().getInt("respa-dt", 4);
                //double in = molecularAssembly.getProperties().getDouble("respa-dt",0.1);
                if (in < 2) {
                    in = 2;
                }
                double inner = dt / in;
                if (!quiet) {
                    logger.log(Level.INFO, String.format(" Created a RESPA integrator with outer %6.3e and inner %6.3e time steps (psec).", dt, inner));
                }
                integrator = addRESPA(inner, dt);
                break;
            /*
            case "BROWNIAN":
                integrator = OpenMM_BrownianIntegrator_create(temperature, frictionCoeff, dt);
                break;
            case "CUSTOM":
                integrator = OpenMM_CustomIntegrator_create(dt);
                break;
            case "COMPOUND":
                integrator = OpenMM_CompoundIntegrator_create();
                break;
             */
            case "VERLET":
            default:
                if (!quiet) {
                    logger.log(Level.INFO, String.format(" Created Verlet integrator with time step %6.3e (psec).", dt));
                }
                integrator = OpenMM_VerletIntegrator_create(dt);
        }

        // Create a context.
        context = OpenMM_Context_create_2(system, integrator, platform);

        // Set initial positions.
        double x[] = new double[numParticles * 3];
        int index = 0;
        Atom atoms[] = molecularAssembly.getAtomArray();
        for (int i = 0; i < numParticles; i++) {
            Atom atom = atoms[i];
            x[index] = atom.getX();
            x[index + 1] = atom.getY();
            x[index + 2] = atom.getZ();
            index += 3;
        }

        setOpenMMPositions(x, numParticles * 3);

        quietInt++;
        //logger.info(String.format(" quietInt is %d", quietInt));
    }

    /**
     * <p>Getter for the field <code>integratorString</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getIntegratorString() {
        return integratorString;
    }

    /**
     * <p>Getter for the field <code>temperature</code>.</p>
     *
     * @return a double.
     */
    public double getTemperature() {
        return temperature;
    }

    /**
     * <p>Getter for the field <code>timeStep</code>.</p>
     *
     * @return a double.
     */
    public double getTimeStep() {
        return timeStep;
    }

    /**
     * <p>setCoeffOfFriction.</p>
     *
     * @param coeffOfFriction a double.
     */
    public void setCoeffOfFriction(double coeffOfFriction) {
        this.frictionCoeff = coeffOfFriction;
    }

    /**
     * <p>Setter for the field <code>collisionFreq</code>.</p>
     *
     * @param collisionFreq a double.
     */
    public void setCollisionFreq(double collisionFreq) {
        this.collisionFreq = collisionFreq;
    }

    /**
     * Create a JNA Pointer to a String.
     *
     * @param string WARNING: assumes ascii-only string
     * @return pointer.
     */
    private Pointer pointerForString(String string) {
        Pointer pointer = new Memory(string.length() + 1);
        pointer.setString(0, string);
        return pointer;
    }

    /**
     * Returns the platform array as a String
     *
     * @param stringArray
     * @param i
     * @return String
     */
    private String stringFromArray(PointerByReference stringArray, int i) {
        Pointer platformPtr = OpenMM_StringArray_get(stringArray, i);
        if (platformPtr == null) {
            return null;
        }
        return platformPtr.getString(0);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean destroy() {
        boolean ffxFFEDestroy = super.destroy();
        freeOpenMM();
        logger.fine(" Destroyed the Context, Integrator, and OpenMMSystem.");
        return ffxFFEDestroy;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void finalize() throws Throwable {
        // Safer to leave super.finalize() in, even though right now that calls Object.finalize().
        logger.info(" ForceFieldEnergyOpenMM instance is being finalized.");
        super.finalize();
        if (destroyed) {
            logger.info(String.format(" Finalize called on a destroyed OpenMM ForceFieldEnergy %s", this.toString()));
        } else {
            destroy();
        }
    }

    /**
     * Destroys pointer references to Context, Integrator and System to free up
     * memory.
     */
    private void freeOpenMM() {
        if (context != null) {
            OpenMM_Context_destroy(context);
            context = null;
        }
        if (integrator != null) {
            OpenMM_Integrator_destroy(integrator);
            integrator = null;
        }
        if (system != null) {
            OpenMM_System_destroy(system);
            system = null;
        }
    }

    /**
     * Adds atoms from the molecular assembly to the OpenMM System and reports
     * to the user the number of particles added.
     */
    private void addAtoms() throws Exception {
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        numParticles = 0;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            OpenMM_System_addParticle(system, atom.getMass());
            if (atom.getMass() <= 0.0) {
                throw new Exception(" Atom without mass greater than 0.");
            }
            numParticles++;
        }
        logger.log(Level.INFO, "  Atoms {0}", nAtoms);
    }

    /**
     * Adds a center-of-mass motion remover to the Potential. Not advised for
     * anything not running MD using the OpenMM library (i.e.
     * OpenMMMolecularDynamics). Has caused bugs with the FFX MD class.
     */
    public void addCOMMRemover() {
        addCOMMRemover(false);
    }

    /**
     * Adds a center-of-mass motion remover to the Potential. Not advised for
     * anything not running MD using the OpenMM library (i.e.
     * OpenMMMolecularDynamics). Has caused bugs with the FFX MD class.
     *
     * @param addIfDuplicate Add a CCOM remover even if it already exists
     */
    public void addCOMMRemover(boolean addIfDuplicate) {
        if (commRemover == null || addIfDuplicate) {
            if (commRemover != null) {
                logger.warning(" Adding a second center-of-mass remover; this is probably incorrect!");
            }
            int frequency = 100;
            commRemover = OpenMM_CMMotionRemover_create(frequency);
            OpenMM_System_addForce(system, commRemover);
            logger.log(Level.INFO, " Added center of mass motion remover (frequency: {0})", frequency);
        } else {
            logger.warning(" Attempted to add a second center-of-mass motion remover when one already exists!");
        }
    }

    /**
     * Add an Andersen thermostat to the system.
     *
     * @param targetTemp Target temperature in Kelvins.
     */
    public void addAndersenThermostat(double targetTemp) {
        addAndersenThermostat(targetTemp, collisionFreq);
    }

    /**
     * Add an Andersen thermostat to the system.
     *
     * @param targetTemp    Target temperature in Kelvins.
     * @param collisionFreq Collision frequency in 1/psec
     */
    public void addAndersenThermostat(double targetTemp, double collisionFreq) {
        if (ommThermostat == null) {
            ommThermostat = OpenMM_AndersenThermostat_create(targetTemp, collisionFreq);
            OpenMM_System_addForce(system, ommThermostat);
            logger.info(format(" Added an Andersen thermostat at %10.6fK and collison frequency %10.6f.", targetTemp, collisionFreq));
        } else {
            logger.info(" Attempted to add a second thermostat to an OpenMM force field!");
        }
    }

    public void addMonteCarloBarostat(double targetPressure, double targetTemp, int frequency) {
        if (ommBarostat == null) {
            ommBarostat = OpenMM_MonteCarloBarostat_create(targetPressure, targetTemp, frequency);
            OpenMM_System_addForce(system, ommBarostat);
            logger.info(format(" Added a Monte Carlo barostat at target pressure %10.6f bar, target temperature %10.6fK and MC move frequency %d.", targetPressure, targetTemp, frequency));

            CompositeConfiguration properties = molecularAssembly.getProperties();

            if (properties.containsKey("randomseed")) {
                int randomSeed = properties.getInt("randomseed", 0);
                logger.info(String.format(" Setting random seed %d for Monte Carlo Barostat", randomSeed));
                OpenMM_MonteCarloBarostat_setRandomNumberSeed(integrator, randomSeed);
            }

        } else {
            logger.info(" Attempted to add a second barostat to an OpenMM force field!");
        }
    }

    /**
     * Adds the 1-2 bond force to the OpenMM System.
     *
     * @param forceField ForceField in use.
     */
    private void addBondForce(ForceField forceField) {
        Bond bonds[] = super.getBonds();
        if (bonds == null || bonds.length < 1) {
            return;
        }
        int nBonds = bonds.length;
        amoebaBondForce = OpenMM_AmoebaBondForce_create();

        double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

        for (int i = 0; i < nBonds; i++) {
            Bond bond = bonds[i];
            int i1 = bond.getAtom(0).getXyzIndex() - 1;
            int i2 = bond.getAtom(1).getXyzIndex() - 1;
            BondType bondType = bond.bondType;
            OpenMM_AmoebaBondForce_addBond(amoebaBondForce, i1, i2,
                    bond.bondType.distance * OpenMM_NmPerAngstrom,
                    kParameterConversion * bondType.forceConstant * BondType.units);

        }

        if (bonds[0].bondType.bondFunction == BondFunction.QUARTIC) {
            OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(amoebaBondForce,
                    BondType.cubic / OpenMM_NmPerAngstrom);
            OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(amoebaBondForce,
                    BondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
        }

        ForceField.ForceFieldInteger bondFGroup = ForceField.ForceFieldInteger.BOND_FORCE_GROUP;
        int fGroup = forceField.getInteger(bondFGroup, bondFGroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaBondForce, fGroup);
        OpenMM_System_addForce(system, amoebaBondForce);
        logger.log(Level.INFO, "  Bonds {0}", nBonds);
    }

    private void addAngleForce(ForceField forceField) {
        Angle angles[] = super.getAngles();
        if (angles == null || angles.length < 1) {
            return;
        }
        int nAngles = angles.length;
        List<Angle> normalAngles = new ArrayList<>();
        // Sort all normal angles from in-plane angles
        for (int i = 0; i < nAngles; i++) {
            if (angles[i].getAngleMode() == Angle.AngleMode.NORMAL) {
                normalAngles.add(angles[i]);
            }
        }
        nAngles = normalAngles.size();
        if (nAngles < 1) {
            return;
        }
        amoebaAngleForce = OpenMM_AmoebaAngleForce_create();
        for (int i = 0; i < nAngles; i++) {
            Angle angle = normalAngles.get(i);
            int i1 = angle.getAtom(0).getXyzIndex() - 1;
            int i2 = angle.getAtom(1).getXyzIndex() - 1;
            int i3 = angle.getAtom(2).getXyzIndex() - 1;
            int nh = angle.nh;
            OpenMM_AmoebaAngleForce_addAngle(amoebaAngleForce, i1, i2, i3,
                    angle.angleType.angle[nh], OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant);
        }

        if (angles[0].angleType.angleFunction == AngleFunction.SEXTIC) {
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(amoebaAngleForce, AngleType.cubic);
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(amoebaAngleForce, AngleType.quartic);
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(amoebaAngleForce, AngleType.quintic);
            OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(amoebaAngleForce, AngleType.sextic);
        }

        ForceField.ForceFieldInteger angleFgroup = ForceField.ForceFieldInteger.ANGLE_FORCE_GROUP;
        int fGroup = forceField.getInteger(angleFgroup, angleFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaAngleForce, fGroup);
        OpenMM_System_addForce(system, amoebaAngleForce);
        logger.log(Level.INFO, "  Angles {0}", nAngles);
    }

    private void addInPlaneAngleForce(ForceField forceField) {
        Angle angles[] = super.getAngles();
        if (angles == null || angles.length < 1) {
            return;
        }
        int nAngles = angles.length;
        List<Angle> inPlaneAngles = new ArrayList<>();
        //Sort all in-plane angles from normal angles
        for (int i = 0; i < nAngles; i++) {
            if (angles[i].getAngleMode() == Angle.AngleMode.IN_PLANE) {
                inPlaneAngles.add(angles[i]);
            }
        }
        nAngles = inPlaneAngles.size();
        if (nAngles < 1) {
            return;
        }
        amoebaInPlaneAngleForce = OpenMM_AmoebaInPlaneAngleForce_create();
        for (int i = 0; i < nAngles; i++) {
            Angle angle = inPlaneAngles.get(i);
            int i1 = angle.getAtom(0).getXyzIndex() - 1;
            int i2 = angle.getAtom(1).getXyzIndex() - 1;
            int i3 = angle.getAtom(2).getXyzIndex() - 1;
            int i4 = angle.getAtom4().getXyzIndex() - 1;
            int nh = angle.nh;
            OpenMM_AmoebaInPlaneAngleForce_addAngle(amoebaInPlaneAngleForce, i1, i2, i3, i4,
                    angle.angleType.angle[nh], OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant);
        }
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleCubic(amoebaInPlaneAngleForce, AngleType.cubic);
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleQuartic(amoebaInPlaneAngleForce, AngleType.quartic);
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAnglePentic(amoebaInPlaneAngleForce, AngleType.quintic);
        OpenMM_AmoebaInPlaneAngleForce_setAmoebaGlobalInPlaneAngleSextic(amoebaInPlaneAngleForce, AngleType.sextic);

        ForceField.ForceFieldInteger inplaneAngleGroup = ForceField.ForceFieldInteger.IN_PLANE_ANGLE_FORCE_GROUP;
        int fGroup = forceField.getInteger(inplaneAngleGroup, inplaneAngleGroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaInPlaneAngleForce, fGroup);
        OpenMM_System_addForce(system, amoebaInPlaneAngleForce);
        logger.log(Level.INFO, "  In-plane Angles {0}", nAngles);
    }

    private void addUreyBradleyForce(ForceField forceField) {
        UreyBradley ureyBradleys[] = super.getUreyBradleys();
        if (ureyBradleys == null || ureyBradleys.length < 1) {
            return;
        }
        amoebaUreyBradleyForce = OpenMM_AmoebaBondForce_create();
        double kParameterConversion = UreyBradleyType.units * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
        int nUreys = ureyBradleys.length;
        for (int i = 0; i < nUreys; i++) {
            UreyBradley ureyBradley = ureyBradleys[i];
            int i1 = ureyBradley.getAtom(0).getXyzIndex() - 1;
            int i2 = ureyBradley.getAtom(2).getXyzIndex() - 1;
            UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
            OpenMM_AmoebaBondForce_addBond(amoebaUreyBradleyForce, i1, i2,
                    ureyBradleyType.distance * OpenMM_NmPerAngstrom,
                    ureyBradleyType.forceConstant * kParameterConversion);
        }

        OpenMM_AmoebaBondForce_setAmoebaGlobalBondCubic(amoebaUreyBradleyForce,
                UreyBradleyType.cubic / OpenMM_NmPerAngstrom);
        OpenMM_AmoebaBondForce_setAmoebaGlobalBondQuartic(amoebaUreyBradleyForce,
                UreyBradleyType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));

        ForceField.ForceFieldInteger ubradFgroup = ForceField.ForceFieldInteger.UREY_BRADLEY_FORCE_GROUP;
        int fGroup = forceField.getInteger(ubradFgroup, ubradFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaUreyBradleyForce, fGroup);
        OpenMM_System_addForce(system, amoebaUreyBradleyForce);
        logger.log(Level.INFO, "  Urey-Bradleys {0}", nUreys);
    }

    private void addOutOfPlaneBendForce(ForceField forceField) {
        OutOfPlaneBend outOfPlaneBends[] = super.getOutOfPlaneBends();
        if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
            return;
        }
        amoebaOutOfPlaneBendForce = OpenMM_AmoebaOutOfPlaneBendForce_create();
        int nOutOfPlaneBends = outOfPlaneBends.length;
        for (int i = 0; i < nOutOfPlaneBends; i++) {
            OutOfPlaneBend outOfPlaneBend = outOfPlaneBends[i];
            int i1 = outOfPlaneBend.getAtom(0).getXyzIndex() - 1;
            int i2 = outOfPlaneBend.getAtom(1).getXyzIndex() - 1;
            int i3 = outOfPlaneBend.getAtom(2).getXyzIndex() - 1;
            int i4 = outOfPlaneBend.getAtom(3).getXyzIndex() - 1;
            OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;

            OpenMM_AmoebaOutOfPlaneBendForce_addOutOfPlaneBend(amoebaOutOfPlaneBendForce, i1, i2, i3, i4,
                    OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * OutOfPlaneBendType.units);
        }
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendCubic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.cubic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendQuartic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.quartic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendPentic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.quintic);
        OpenMM_AmoebaOutOfPlaneBendForce_setAmoebaGlobalOutOfPlaneBendSextic(amoebaOutOfPlaneBendForce, OutOfPlaneBendType.sextic);

        ForceField.ForceFieldInteger oopBendFgroup = ForceField.ForceFieldInteger.OUT_OF_PLANE_BEND_FORCE_GROUP;
        int fGroup = forceField.getInteger(oopBendFgroup, oopBendFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaOutOfPlaneBendForce, fGroup);
        OpenMM_System_addForce(system, amoebaOutOfPlaneBendForce);
        logger.log(Level.INFO, "  Out-of-Plane Bends {0}", nOutOfPlaneBends);
    }

    private void addStretchBendForce(ForceField forceField) {
        StretchBend stretchBends[] = super.getStretchBends();
        if (stretchBends == null || stretchBends.length < 1) {
            return;
        }
        int nStretchBends = stretchBends.length;
        amoebaStretchBendForce = OpenMM_AmoebaStretchBendForce_create();
        for (int i = 0; i < nStretchBends; i++) {
            StretchBend stretchBend = stretchBends[i];
            int i1 = stretchBend.getAtom(0).getXyzIndex() - 1;
            int i2 = stretchBend.getAtom(1).getXyzIndex() - 1;
            int i3 = stretchBend.getAtom(2).getXyzIndex() - 1;
            double angle = stretchBend.angleEq;
            double beq0 = stretchBend.bond0Eq;
            double beq1 = stretchBend.bond1Eq;
            double fc0 = stretchBend.force0;
            double fc1 = stretchBend.force1;
            OpenMM_AmoebaStretchBendForce_addStretchBend(amoebaStretchBendForce, i1, i2, i3,
                    beq0 * OpenMM_NmPerAngstrom, beq1 * OpenMM_NmPerAngstrom, OpenMM_RadiansPerDegree * angle,
                    (OpenMM_KJPerKcal / OpenMM_NmPerAngstrom) * fc0, (OpenMM_KJPerKcal / OpenMM_NmPerAngstrom) * fc1);

        }

        ForceField.ForceFieldInteger sbendFgroup = ForceField.ForceFieldInteger.STRETCH_BEND_FORCE_GROUP;
        int fGroup = forceField.getInteger(sbendFgroup, sbendFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaStretchBendForce, fGroup);
        OpenMM_System_addForce(system, amoebaStretchBendForce);
        logger.log(Level.INFO, "  Stretch-Bends {0}", nStretchBends);
    }

    private void addTorsionForce(ForceField forceField) {
        Torsion torsions[] = super.getTorsions();
        if (torsions == null || torsions.length < 1) {
            return;
        }
        int nTorsions = torsions.length;
        amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create();
        for (int i = 0; i < nTorsions; i++) {
            Torsion torsion = torsions[i];
            int a1 = torsion.getAtom(0).getXyzIndex() - 1;
            int a2 = torsion.getAtom(1).getXyzIndex() - 1;
            int a3 = torsion.getAtom(2).getXyzIndex() - 1;
            int a4 = torsion.getAtom(3).getXyzIndex() - 1;
            TorsionType torsionType = torsion.torsionType;
            int nTerms = torsionType.phase.length;
            for (int j = 0; j < nTerms; j++) {
                OpenMM_PeriodicTorsionForce_addTorsion(amoebaTorsionForce,
                        a1, a2, a3, a4, j + 1,
                        torsionType.phase[j] * OpenMM_RadiansPerDegree,
                        OpenMM_KJPerKcal * torsion.units * torsionType.amplitude[j]);
            }
        }

        ForceField.ForceFieldInteger torsFgroup = ForceField.ForceFieldInteger.TORSION_FORCE_GROUP;
        int fGroup = forceField.getInteger(torsFgroup, torsFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaTorsionForce, fGroup);
        OpenMM_System_addForce(system, amoebaTorsionForce);
        logger.log(Level.INFO, "  Torsions {0}", nTorsions);
    }

    private void addImproperTorsionForce(ForceField forceField) {
        ImproperTorsion impropers[] = super.getImproperTorsions();
        if (impropers == null || impropers.length < 1) {
            return;
        }
        int nImpropers = impropers.length;
        amoebaImproperTorsionForce = OpenMM_PeriodicTorsionForce_create();

        for (int i = 0; i < nImpropers; i++) {
            ImproperTorsion improperTorsion = impropers[i];
            int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
            int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
            int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
            int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
            ImproperTorsionType improperTorsionType = improperTorsion.improperType;
            OpenMM_PeriodicTorsionForce_addTorsion(amoebaImproperTorsionForce,
                    a1, a2, a3, a4, improperTorsionType.periodicity,
                    improperTorsionType.phase * OpenMM_RadiansPerDegree,
                    OpenMM_KJPerKcal * improperTorsion.units
                            * improperTorsion.scaleFactor * improperTorsionType.k);
        }

        ForceField.ForceFieldInteger imptorsFgroup = ForceField.ForceFieldInteger.IMPROPER_TORSION_FORCE_GROUP;
        int fGroup = forceField.getInteger(imptorsFgroup, imptorsFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaImproperTorsionForce, fGroup);
        OpenMM_System_addForce(system, amoebaImproperTorsionForce);
        logger.log(Level.INFO, "  Improper Torsions {0} ", nImpropers);
    }

    private void addPiTorsionForce(ForceField forceField) {
        PiOrbitalTorsion piOrbitalTorsions[] = super.getPiOrbitalTorsions();
        if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
            return;
        }
        int nPiOrbitalTorsions = piOrbitalTorsions.length;
        amoebaPiTorsionForce = OpenMM_AmoebaPiTorsionForce_create();
        double units = PiTorsionType.units;
        for (int i = 0; i < nPiOrbitalTorsions; i++) {
            PiOrbitalTorsion piOrbitalTorsion = piOrbitalTorsions[i];
            int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
            int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
            int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
            int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
            int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
            int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
            PiTorsionType type = piOrbitalTorsion.piTorsionType;
            OpenMM_AmoebaPiTorsionForce_addPiTorsion(amoebaPiTorsionForce,
                    a1, a2, a3, a4, a5, a6,
                    OpenMM_KJPerKcal * type.forceConstant * units);
        }

        ForceField.ForceFieldInteger piOrbTorsFgroup = ForceField.ForceFieldInteger.PI_ORBITAL_TORSION_FORCE_GROUP;
        int fGroup = forceField.getInteger(piOrbTorsFgroup, piOrbTorsFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaPiTorsionForce, fGroup);
        OpenMM_System_addForce(system, amoebaPiTorsionForce);
        logger.log(Level.INFO, "  Pi-Orbital Torsions {0}", nPiOrbitalTorsions);
    }

    private void addTorsionTorsionForce(ForceField forceField) {
        TorsionTorsion torsionTorsions[] = super.getTorsionTorsions();
        if (torsionTorsions == null || torsionTorsions.length < 1) {
            return;
        }
        /**
         * Load the torsion-torsions.
         */

        int nTypes = 0;
        LinkedHashMap<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();

        int nTorsionTorsions = torsionTorsions.length;
        amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create();
        for (int i = 0; i < nTorsionTorsions; i++) {
            TorsionTorsion torsionTorsion = torsionTorsions[i];
            int ia = torsionTorsion.getAtom(0).getXyzIndex() - 1;
            int ib = torsionTorsion.getAtom(1).getXyzIndex() - 1;
            int ic = torsionTorsion.getAtom(2).getXyzIndex() - 1;
            int id = torsionTorsion.getAtom(3).getXyzIndex() - 1;
            int ie = torsionTorsion.getAtom(4).getXyzIndex() - 1;

            TorsionTorsionType torsionTorsionType = torsionTorsion.torsionTorsionType;
            String key = torsionTorsionType.getKey();
            /**
             * Check if the TorTor parameters have already been added to the
             * Hash.
             */
            int gridIndex = 0;
            if (torTorTypes.containsKey(key)) {
                /**
                 * If the TorTor has been added, get its (ordered) index in the
                 * Hash.
                 */
                int index = 0;
                for (String entry : torTorTypes.keySet()) {
                    if (entry.equalsIgnoreCase(key)) {
                        gridIndex = index;
                        break;
                    } else {
                        index++;
                    }
                }
            } else {
                /**
                 * Add the new TorTor.
                 */
                torTorTypes.put(key, torsionTorsionType);
                gridIndex = nTypes;
                nTypes++;
            }

            Atom atom = torsionTorsion.getChiralAtom();
            int iChiral = -1;
            if (atom != null) {
                iChiral = atom.getXyzIndex() - 1;
            }
            OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(amoebaTorsionTorsionForce,
                    ia, ib, ic, id, ie, iChiral, gridIndex);
        }
        /**
         * Load the Torsion-Torsion parameters.
         */
        PointerByReference values = OpenMM_DoubleArray_create(6);
        int gridIndex = 0;
        for (String key : torTorTypes.keySet()) {
            TorsionTorsionType torTorType = torTorTypes.get(key);
            int nx = torTorType.nx;
            int ny = torTorType.ny;
            double tx[] = torTorType.tx;
            double ty[] = torTorType.ty;
            double f[] = torTorType.energy;
            double dx[] = torTorType.dx;
            double dy[] = torTorType.dy;
            double dxy[] = torTorType.dxy;
            /**
             * Create the 3D grid.
             */
            PointerByReference grid3D = OpenMM_3D_DoubleArray_create(nx, ny, 6);
            int xIndex = 0;
            int yIndex = 0;
            for (int j = 0; j < nx * ny; j++) {
                int addIndex = 0;
                OpenMM_DoubleArray_set(values, addIndex++, tx[xIndex]);
                OpenMM_DoubleArray_set(values, addIndex++, ty[yIndex]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * f[j]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dx[j]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dy[j]);
                OpenMM_DoubleArray_set(values, addIndex++, OpenMM_KJPerKcal * dxy[j]);
                OpenMM_3D_DoubleArray_set(grid3D, yIndex, xIndex, values);
                xIndex++;
                if (xIndex == nx) {
                    xIndex = 0;
                    yIndex++;
                }
            }
            OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(amoebaTorsionTorsionForce, gridIndex++, grid3D);
            OpenMM_3D_DoubleArray_destroy(grid3D);
        }
        OpenMM_DoubleArray_destroy(values);

        ForceField.ForceFieldInteger torTorFgroup = ForceField.ForceFieldInteger.TORSION_TORSION_FORCE_GROUP;
        int fGroup = forceField.getInteger(torTorFgroup, torTorFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaTorsionTorsionForce, fGroup);
        OpenMM_System_addForce(system, amoebaTorsionTorsionForce);
        logger.log(Level.INFO, "  Torsion-Torsions {0}", nTorsionTorsions);
    }

    /**
     * Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic
     * acid model).
     */
    private void addStretchTorsionCoupling(ForceField forceField) {
        StretchTorsion[] strTorsions = super.getStretchTorsions();
        if (strTorsions != null && strTorsions.length > 0) {
            stretchTorsionForce = OpenMM_CustomCompoundBondForce_create(4, StretchTorsion.stretchTorsionForm());
            OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi1", 0);
            OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi2", Math.PI);
            OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi3", 0);

            for (int m = 1; m < 4; m++) {
                for (int n = 1; n < 4; n++) {
                    OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchTorsionForce, String.format("k%d%d", m, n));
                }
            }

            for (int m = 1; m < 4; m++) {
                OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchTorsionForce, String.format("b%d", m));
            }

            final double unitConv = OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;

            for (StretchTorsion strTors : strTorsions) {
                double[] constants = strTors.getConstants();
                PointerByReference strTorsParams = OpenMM_DoubleArray_create(0);
                for (int m = 0; m < 3; m++) {
                    for (int n = 0; n < 3; n++) {
                        int index = (3 * m) + n;
                        double kmn = constants[index] * unitConv;
                        OpenMM_DoubleArray_append(strTorsParams, kmn);
                    }
                }

                OpenMM_DoubleArray_append(strTorsParams, strTors.bondType1.distance * OpenMM_NmPerAngstrom);
                OpenMM_DoubleArray_append(strTorsParams, strTors.bondType2.distance * OpenMM_NmPerAngstrom);
                OpenMM_DoubleArray_append(strTorsParams, strTors.bondType3.distance * OpenMM_NmPerAngstrom);

                PointerByReference strTorsParticles = OpenMM_IntArray_create(0);
                Atom[] atoms = strTors.getAtomArray(true);
                for (int i = 0; i < 4; i++) {
                    OpenMM_IntArray_append(strTorsParticles, atoms[i].getXyzIndex() - 1);
                }

                OpenMM_CustomCompoundBondForce_addBond(stretchTorsionForce, strTorsParticles, strTorsParams);
                OpenMM_DoubleArray_destroy(strTorsParams);
                OpenMM_IntArray_destroy(strTorsParticles);
            }

            ForceField.ForceFieldInteger sTorFgroup = ForceField.ForceFieldInteger.STRETCH_TORSION_FORCE_GROUP;
            int fGroup = forceField.getInteger(sTorFgroup, sTorFgroup.getDefaultValue());

            OpenMM_Force_setForceGroup(stretchTorsionForce, fGroup);
            OpenMM_System_addForce(system, stretchTorsionForce);
        }
    }

    /**
     * Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic
     * acid model).
     */
    private void addAngleTorsionCoupling(ForceField forceField) {
        AngleTorsion[] angleTorsions = super.getAngleTorsions();
        if (angleTorsions != null && angleTorsions.length > 0) {
            angleTorsionForce = OpenMM_CustomCompoundBondForce_create(4, AngleTorsion.angleTorsionForm());
            OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi1", 0);
            OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi2", Math.PI);
            OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi3", 0);

            for (int m = 1; m < 3; m++) {
                for (int n = 1; n < 4; n++) {
                    OpenMM_CustomCompoundBondForce_addPerBondParameter(angleTorsionForce, String.format("k%d%d", m, n));
                }
            }

            for (int m = 1; m < 3; m++) {
                OpenMM_CustomCompoundBondForce_addPerBondParameter(angleTorsionForce, String.format("a%d", m));
            }

            final double unitConv = OpenMM_KJPerKcal / OpenMM_RadiansPerDegree;

            for (AngleTorsion ators : angleTorsions) {
                double[] constants = ators.getConstants();
                PointerByReference atorsParams = OpenMM_DoubleArray_create(0);
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 3; n++) {
                        int index = (3 * m) + n;
                        double kmn = constants[index] * unitConv;
                        OpenMM_DoubleArray_append(atorsParams, kmn);
                    }
                }

                Atom[] atoms = ators.getAtomArray(true);

                // One thing that concerns me is whether it's correct to get angle[0] instead of angle[num hydrogens].
                // This is the way it is in FFX, but that may be a bug.
                OpenMM_DoubleArray_append(atorsParams, ators.angleType1.angle[0] * OpenMM_RadiansPerDegree);
                OpenMM_DoubleArray_append(atorsParams, ators.angleType2.angle[0] * OpenMM_RadiansPerDegree);

                PointerByReference atorsParticles = OpenMM_IntArray_create(0);
                for (int i = 0; i < 4; i++) {
                    OpenMM_IntArray_append(atorsParticles, atoms[i].getXyzIndex() - 1);
                }

                OpenMM_CustomCompoundBondForce_addBond(angleTorsionForce, atorsParticles, atorsParams);
                OpenMM_DoubleArray_destroy(atorsParams);
                OpenMM_IntArray_destroy(atorsParticles);
            }

            ForceField.ForceFieldInteger aTorFgroup = ForceField.ForceFieldInteger.ANGLE_TORSION_FORCE_GROUP;
            int fGroup = forceField.getInteger(aTorFgroup, aTorFgroup.getDefaultValue());

            OpenMM_Force_setForceGroup(angleTorsionForce, fGroup);
            OpenMM_System_addForce(system, angleTorsionForce);
        }
    }

    /**
     * Uses arithmetic mean to define sigma and geometric mean for epsilon.
     */
    private void addFixedChargeNonBondedForce(ForceField forceField) {
        VanDerWaals vdW = super.getVdwNode();
        if (vdW == null) {
            return;
        }
        /**
         * Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
         * for epsilon is supported.
         */
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        if (vdwForm.vdwType != LENNARD_JONES
                || vdwForm.radiusRule != ARITHMETIC
                || vdwForm.epsilonRule != GEOMETRIC) {
            logger.info(format(" VDW Type:         %s", vdwForm.vdwType));
            logger.info(format(" VDW Radius Rule:  %s", vdwForm.radiusRule));
            logger.info(format(" VDW Epsilon Rule: %s", vdwForm.epsilonRule));
            logger.log(Level.SEVERE, " Unsupported van der Waals functional form.");
            return;
        }

        fixedChargeNonBondedForce = OpenMM_NonbondedForce_create();

        /**
         * OpenMM vdW force requires a diameter (i.e. not radius).
         */
        double radScale = 1.0;
        if (vdwForm.radiusSize == RADIUS) {
            radScale = 2.0;
        }
        /**
         * OpenMM vdw force requires atomic sigma values (i.e. not r-min).
         */
        if (vdwForm.radiusType == R_MIN) {
            radScale /= 1.122462048309372981;
        }

        /**
         * Add particles.
         */
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];

            VDWType vdwType = atom.getVDWType();
            double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
            double eps = OpenMM_KJPerKcal * vdwType.wellDepth;

            double charge = 0.0;
            MultipoleType multipoleType = atom.getMultipoleType();
            if (multipoleType != null && atoms[i].getElectrostatics()) {
                charge = multipoleType.charge;
            }

            OpenMM_NonbondedForce_addParticle(fixedChargeNonBondedForce, charge, sigma, eps);
        }

        /**
         * Define 1-4 scale factors.
         */
        double lj14Scale = vdwForm.getScale14();
        double coulomb14Scale = 1.0 / 1.2;

        ParticleMeshEwald pme = super.getPmeNode();
        Bond bonds[] = super.getBonds();
        if (bonds != null && bonds.length > 0) {
            int nBonds = bonds.length;
            PointerByReference bondArray;
            bondArray = OpenMM_BondArray_create(0);
            for (int i = 0; i < nBonds; i++) {
                Bond bond = bonds[i];
                int i1 = bond.getAtom(0).getXyzIndex() - 1;
                int i2 = bond.getAtom(1).getXyzIndex() - 1;
                OpenMM_BondArray_append(bondArray, i1, i2);
            }
            if (pme != null) {
                coulomb14Scale = pme.getScale14();
            }
            OpenMM_NonbondedForce_createExceptionsFromBonds(fixedChargeNonBondedForce, bondArray, coulomb14Scale, lj14Scale);

            OpenMM_BondArray_destroy(bondArray);

            int num = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);
            chargeExclusion = new boolean[num];
            vdWExclusion = new boolean[num];
            exceptionChargeProd = new double[num];
            exceptionEps = new double[num];

            IntByReference particle1 = new IntByReference();
            IntByReference particle2 = new IntByReference();
            DoubleByReference chargeProd = new DoubleByReference();
            DoubleByReference sigma = new DoubleByReference();
            DoubleByReference eps = new DoubleByReference();

            for (int i = 0; i < num; i++) {
                OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i,
                        particle1, particle2, chargeProd, sigma, eps);
                if (abs(chargeProd.getValue()) > 0.0) {
                    chargeExclusion[i] = false;
                    exceptionChargeProd[i] = chargeProd.getValue();
                } else {
                    exceptionChargeProd[i] = 0.0;
                    chargeExclusion[i] = true;
                }
                if (abs(eps.getValue()) > 0.0) {
                    vdWExclusion[i] = false;
                    exceptionEps[i] = eps.getValue();
                } else {
                    vdWExclusion[i] = true;
                    exceptionEps[i] = 0.0;
                }
            }
        }

        Crystal crystal = super.getCrystal();
        if (crystal.aperiodic()) {
            OpenMM_NonbondedForce_setNonbondedMethod(fixedChargeNonBondedForce,
                    OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_NoCutoff);
        } else {
            OpenMM_NonbondedForce_setNonbondedMethod(fixedChargeNonBondedForce,
                    OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_PME);

            if (pme != null) {
                // Units of the Ewald coefficient are A^-1; Multiply by AngstromsPerNM to convert to (Nm^-1).
                double aEwald = OpenMM_AngstromsPerNm * pme.getEwaldCoefficient();
                int nx = pme.getReciprocalSpace().getXDim();
                int ny = pme.getReciprocalSpace().getYDim();
                int nz = pme.getReciprocalSpace().getZDim();
                OpenMM_NonbondedForce_setPMEParameters(fixedChargeNonBondedForce, aEwald, nx, ny, nz);
            }

            NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
            double off = nonbondedCutoff.off;
            double cut = nonbondedCutoff.cut;
            OpenMM_NonbondedForce_setCutoffDistance(fixedChargeNonBondedForce, OpenMM_NmPerAngstrom * off);
            OpenMM_NonbondedForce_setUseSwitchingFunction(fixedChargeNonBondedForce, OpenMM_True);
            if (cut == off) {
                logger.warning(" OpenMM does not properly handle cutoffs where cut == off!");
                if (cut == Double.MAX_VALUE || cut == Double.POSITIVE_INFINITY) {
                    logger.info(" Detected infinite or max-value cutoff; setting cut to 1E+40 for OpenMM.");
                    cut = 1E40;
                } else {
                    logger.info(String.format(" Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut, off));
                    cut *= 0.99;
                }
            }
            OpenMM_NonbondedForce_setSwitchingDistance(fixedChargeNonBondedForce, OpenMM_NmPerAngstrom * cut);
        }

        OpenMM_NonbondedForce_setUseDispersionCorrection(fixedChargeNonBondedForce, OpenMM_False);

        ForceField.ForceFieldInteger vdwFgroup = ForceField.ForceFieldInteger.VDW_FORCE_GROUP;
        int fGroup = forceField.getInteger(vdwFgroup, vdwFgroup.getDefaultValue());

        ForceField.ForceFieldInteger electroFgroup = ForceField.ForceFieldInteger.PME_FORCE_GROUP;
        int pmeGroup = forceField.getInteger(electroFgroup, electroFgroup.getDefaultValue());

        if (fGroup != pmeGroup) {
            logger.severe(String.format(" ERROR: VDW-FORCE-GROUP is %d while PME-FORCE-GROUP is %d. "
                    + "This is invalid for fixed-charge force fields with combined nonbonded forces.", fGroup, pmeGroup));
        }

        OpenMM_Force_setForceGroup(fixedChargeNonBondedForce, fGroup);
        OpenMM_System_addForce(system, fixedChargeNonBondedForce);
        logger.log(Level.INFO, String.format("  Fixed charge non-bonded force"));

        GeneralizedKirkwood gk = super.getGK();
        if (gk != null) {
            addCustomGBForce(forceField);
        }
    }

    /**
     * 1.) Handle interactions between non-alchemical atoms with our default
     * OpenMM NonBondedForce. Note that alchemical atoms must have eps=0 to turn
     * them off in this force.
     * <p>
     * 2.) Handle interactions between alchemical atoms and mixed non-alchemical
     * <-> alchemical interactions with an OpenMM CustomNonBondedForce.
     */
    private void addCustomNonbondedSoftcoreForce(ForceField forceField) {

        VanDerWaals vdW = super.getVdwNode();
        if (vdW == null) {
            return;
        }

        /**
         * Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
         * for epsilon is supported.
         */
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        if (vdwForm.vdwType != LENNARD_JONES
                || vdwForm.radiusRule != ARITHMETIC
                || vdwForm.epsilonRule != GEOMETRIC) {
            logger.info(format(" VDW Type:         %s", vdwForm.vdwType));
            logger.info(format(" VDW Radius Rule:  %s", vdwForm.radiusRule));
            logger.info(format(" VDW Epsilon Rule: %s", vdwForm.epsilonRule));
            logger.log(Level.SEVERE, String.format(" Unsupported van der Waals functional form."));
            return;
        }

        // Sterics mixing rules.
        String stericsMixingRules = " epsilon = sqrt(epsilon1*epsilon2);";
        stericsMixingRules += " rmin = 0.5 * (sigma1 + sigma2) * 1.122462048309372981;";

        // Softcore Lennard-Jones, with a form equivalent to that used in FFX VanDerWaals class.
        String stericsEnergyExpression = "(vdw_lambda^beta)*epsilon*x*(x-2.0);";
        // Effective softcore distance for sterics.
        stericsEnergyExpression += " x = 1.0 / (alpha*(1.0-vdw_lambda)^2.0 + (r/rmin)^6.0);";
        // Define energy expression for sterics.
        String energyExpression = stericsEnergyExpression + stericsMixingRules;

        fixedChargeSoftcore = OpenMM_CustomNonbondedForce_create(energyExpression);

        // Get the Alpha and Beta constants from the VanDerWaals instance.
        double alpha = vdW.getAlpha();
        double beta = vdW.getBeta();

        logger.info(format(" Custom non-bonded force with alpha = %8.6f and beta = %8.6f", alpha, beta));

        OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "vdw_lambda", 1.0);
        OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "alpha", alpha);
        OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "beta", beta);
        OpenMM_CustomNonbondedForce_addPerParticleParameter(fixedChargeSoftcore, "sigma");
        OpenMM_CustomNonbondedForce_addPerParticleParameter(fixedChargeSoftcore, "epsilon");

        /**
         * Add particles.
         */
        PointerByReference alchemicalGroup = OpenMM_IntSet_create();
        PointerByReference nonAlchemicalGroup = OpenMM_IntSet_create();
        DoubleByReference charge = new DoubleByReference();
        DoubleByReference sigma = new DoubleByReference();
        DoubleByReference eps = new DoubleByReference();

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            OpenMM_NonbondedForce_getParticleParameters(fixedChargeNonBondedForce, i, charge, sigma, eps);
            if (atom.applyLambda()) {
                OpenMM_IntSet_insert(alchemicalGroup, i);
                logger.info(format(" Adding alchemical atom %s.", atom));
            } else {
                OpenMM_IntSet_insert(nonAlchemicalGroup, i);
            }
            PointerByReference particleParameters = OpenMM_DoubleArray_create(0);
            OpenMM_DoubleArray_append(particleParameters, sigma.getValue());
            OpenMM_DoubleArray_append(particleParameters, eps.getValue());
            OpenMM_CustomNonbondedForce_addParticle(fixedChargeSoftcore, particleParameters);
            OpenMM_DoubleArray_destroy(particleParameters);
        }

        OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup, alchemicalGroup);
        OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup, nonAlchemicalGroup);
        OpenMM_IntSet_destroy(alchemicalGroup);
        OpenMM_IntSet_destroy(nonAlchemicalGroup);

        Crystal crystal = super.getCrystal();
        if (crystal.aperiodic()) {
            OpenMM_CustomNonbondedForce_setNonbondedMethod(fixedChargeSoftcore,
                    OpenMM_CustomNonbondedForce_NonbondedMethod.OpenMM_CustomNonbondedForce_NoCutoff);
        } else {
            OpenMM_CustomNonbondedForce_setNonbondedMethod(fixedChargeSoftcore,
                    OpenMM_CustomNonbondedForce_NonbondedMethod.OpenMM_CustomNonbondedForce_CutoffPeriodic);
        }

        NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
        double off = nonbondedCutoff.off;
        double cut = nonbondedCutoff.cut;
        if (cut == off) {
            logger.warning(" OpenMM does not properly handle cutoffs where cut == off!");
            if (cut == Double.MAX_VALUE || cut == Double.POSITIVE_INFINITY) {
                logger.info(" Detected infinite or max-value cutoff; setting cut to 1E+40 for OpenMM.");
                cut = 1E40;
            } else {
                logger.info(String.format(" Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut, off));
                cut *= 0.99;
            }
        }

        OpenMM_CustomNonbondedForce_setCutoffDistance(fixedChargeSoftcore, OpenMM_NmPerAngstrom * off);
        OpenMM_CustomNonbondedForce_setUseSwitchingFunction(fixedChargeSoftcore, OpenMM_True);
        OpenMM_CustomNonbondedForce_setSwitchingDistance(fixedChargeSoftcore, OpenMM_NmPerAngstrom * cut);

        // Add energy parameter derivative
        OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(fixedChargeSoftcore, "vdw_lambda");

        ForceField.ForceFieldInteger vdwFgroup = ForceField.ForceFieldInteger.VDW_FORCE_GROUP;
        int fGroup = forceField.getInteger(vdwFgroup, vdwFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(fixedChargeSoftcore, fGroup);
        OpenMM_System_addForce(system, fixedChargeSoftcore);
        logger.log(Level.INFO, String.format(" Added fixed charge softcore sterics force."));

        // Alchemical with Alchemical could be either softcore or normal interactions (softcore here).
        alchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(stericsEnergyExpression);

        // Non-Alchemical with Alchemical is essentially always softcore.
        nonAlchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(stericsEnergyExpression);

        // Currently both are treated the same (so we could condense the code below).
        OpenMM_CustomBondForce_addPerBondParameter(alchemicalAlchemicalStericsForce, "rmin");
        OpenMM_CustomBondForce_addPerBondParameter(alchemicalAlchemicalStericsForce, "epsilon");
        OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "vdw_lambda", 1.0);
        OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "alpha", alpha);
        OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "beta", beta);

        OpenMM_CustomBondForce_addPerBondParameter(nonAlchemicalAlchemicalStericsForce, "rmin");
        OpenMM_CustomBondForce_addPerBondParameter(nonAlchemicalAlchemicalStericsForce, "epsilon");
        OpenMM_CustomBondForce_addGlobalParameter(nonAlchemicalAlchemicalStericsForce, "vdw_lambda", 1.0);
        OpenMM_CustomBondForce_addGlobalParameter(nonAlchemicalAlchemicalStericsForce, "alpha", alpha);
        OpenMM_CustomBondForce_addGlobalParameter(nonAlchemicalAlchemicalStericsForce, "beta", beta);

        int range = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);

        IntByReference atomi = new IntByReference();
        IntByReference atomj = new IntByReference();
        int torsionMask[][] = vdW.getTorsionMask();

        for (int i = 0; i < range; i++) {
            OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i, atomi, atomj, charge, sigma, eps);

            // Omit both Exclusions (1-2, 1-3) and Exceptions (scaled 1-4) from the CustomNonbondedForce.
            OpenMM_CustomNonbondedForce_addExclusion(fixedChargeSoftcore, atomi.getValue(), atomj.getValue());

            // Deal with scaled 1-4 torsions using the CustomBondForce
            int maskI[] = torsionMask[atomi.getValue()];
            int jID = atomj.getValue();
            boolean epsException = false;
            for (int j = 0; j < maskI.length; j++) {
                if (maskI[j] == jID) {
                    epsException = true;
                    break;
                }
            }

            if (epsException) {
                Atom atom1 = atoms[atomi.getValue()];
                Atom atom2 = atoms[atomj.getValue()];

                boolean bothAlchemical = false;
                boolean oneAlchemical = false;

                if (atom1.applyLambda() && atom2.applyLambda()) {
                    bothAlchemical = true;
                } else if ((atom1.applyLambda() && !atom2.applyLambda()) || (!atom1.applyLambda() && atom2.applyLambda())) {
                    oneAlchemical = true;
                }

                if (bothAlchemical) {
                    PointerByReference bondParameters = OpenMM_DoubleArray_create(0);
                    OpenMM_DoubleArray_append(bondParameters, sigma.getValue() * 1.122462048309372981);
                    OpenMM_DoubleArray_append(bondParameters, eps.getValue());
                    OpenMM_CustomBondForce_addBond(alchemicalAlchemicalStericsForce, atomi.getValue(), atomj.getValue(), bondParameters);
                    OpenMM_DoubleArray_destroy(bondParameters);
                } else if (oneAlchemical) {
                    PointerByReference bondParameters = OpenMM_DoubleArray_create(0);
                    OpenMM_DoubleArray_append(bondParameters, sigma.getValue() * 1.122462048309372981);
                    OpenMM_DoubleArray_append(bondParameters, eps.getValue());
                    OpenMM_CustomBondForce_addBond(nonAlchemicalAlchemicalStericsForce, atomi.getValue(), atomj.getValue(), bondParameters);
                    OpenMM_DoubleArray_destroy(bondParameters);
                }
            }
        }

//        for (int i = 0; i < range; i++){
//            OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i, atomi, atomj, charge, sigma, eps);
//
//            Atom atom1 = atoms[atomi.getValue()];
//            Atom atom2 = atoms[atomj.getValue()];
//
//            if (atom1.applyLambda() || atom2.applyLambda()){
//                OpenMM_NonbondedForce_setExceptionParameters(fixedChargeNonBondedForce, i, atomi.getValue(), atomj.getValue(), abs(0.0*charge.getValue()), sigma.getValue(), abs(0.0*eps.getValue()));
//            }
//        }

        OpenMM_CustomBondForce_addEnergyParameterDerivative(alchemicalAlchemicalStericsForce, "vdw_lambda");
        OpenMM_System_addForce(system, alchemicalAlchemicalStericsForce);
        OpenMM_CustomBondForce_addEnergyParameterDerivative(nonAlchemicalAlchemicalStericsForce, "vdw_lambda");
        OpenMM_System_addForce(system, nonAlchemicalAlchemicalStericsForce);

    }

    private void addCustomGBForce(ForceField forceField) {
        GeneralizedKirkwood gk = super.getGK();
        if (gk == null) {
            return;
        }

        double sTens = 0.0;
        if (gk.getNonPolarModel() == NonPolar.BORN_SOLV || gk.getNonPolarModel() == NonPolar.BORN_CAV_DISP) {
            sTens = gk.getSurfaceTension();
            sTens *= OpenMM_KJPerKcal;
            sTens *= 100.0; // 100 square Angstroms per square nanometer.
            // logger.info(String.format(" FFX surface tension: %9.5g kcal/mol/Ang^2", sTens));
            // logger.info(String.format(" OpenMM surface tension: %9.5g kJ/mol/nm^2", sTens));
        }

        customGBForce = OpenMM_CustomGBForce_create();
        OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "q");
        OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "radius");
        OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "scale");
        OpenMM_CustomGBForce_addPerParticleParameter(customGBForce, "surfaceTension");

        OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "solventDielectric", gk.getSolventPermittivity());
        OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "soluteDielectric", 1.0);
        OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "dOffset", gk.getDielecOffset() * OpenMM_NmPerAngstrom); // Factor of 0.1 for Ang to nm.
        OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "probeRadius", gk.getProbeRadius() * OpenMM_NmPerAngstrom);

        OpenMM_CustomGBForce_addComputedValue(customGBForce, "I",
                // "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                // "step(r+sr2-or1)*0.5*((1/L^3-1/U^3)/3+(1/U^4-1/L^4)/8*(r-sr2*sr2/r)+0.25*(1/U^2-1/L^2)/r+C);"
                "0.5*((1/L^3-1/U^3)/3.0+(1/U^4-1/L^4)/8.0*(r-sr2*sr2/r)+0.25*(1/U^2-1/L^2)/r+C);"
                        + "U=r+sr2;"
                        // + "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                        + "C=2/3*(1/or1^3-1/L^3)*step(sr2-r-or1);"
                        // + "L=step(or1-D)*or1 + (1-step(or1-D))*D;"
                        // + "D=step(r-sr2)*(r-sr2) + (1-step(r-sr2))*(sr2-r);"
                        + "L = step(sr2 - r1r)*sr2mr + (1 - step(sr2 - r1r))*L;"
                        + "sr2mr = sr2 - r;"
                        + "r1r = radius1 + r;"
                        + "L = step(r1sr2 - r)*radius1 + (1 - step(r1sr2 - r))*L;"
                        + "r1sr2 = radius1 + sr2;"
                        + "L = r - sr2;"
                        + "sr2 = scale2 * radius2;"
                        + "or1 = radius1; or2 = radius2",
                OpenMM_CustomGBForce_ParticlePairNoExclusions);

        OpenMM_CustomGBForce_addComputedValue(customGBForce, "B",
                // "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                // "psi=I*or; or=radius-0.009"
                "step(BB-radius)*BB + (1 - step(BB-radius))*radius;"
                        + "BB = 1 / ( (3.0*III)^(1.0/3.0) );"
                        + "III = step(II)*II + (1 - step(II))*1.0e-9/3.0;"
                        + "II = maxI - I;"
                        + "maxI = 1/(3.0*radius^3)",
                OpenMM_CustomGBForce_SingleParticle);

        OpenMM_CustomGBForce_addEnergyTerm(customGBForce,
                "surfaceTension*(radius+probeRadius+dOffset)^2*((radius+dOffset)/B)^6/6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B",
                OpenMM_CustomGBForce_SingleParticle);

        /**
         * Particle pair term is the generalized Born cross term.
         */
        OpenMM_CustomGBForce_addEnergyTerm(customGBForce,
                "-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                        + "f=sqrt(r^2+B1*B2*exp(-r^2/(2.455*B1*B2)))",
                OpenMM_CustomGBForce_ParticlePair);

        double baseRadii[] = gk.getBaseRadii();
        double overlapScale[] = gk.getOverlapScale();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        PointerByReference doubleArray = OpenMM_DoubleArray_create(0);
        for (int i = 0; i < nAtoms; i++) {
            MultipoleType multipoleType = atoms[i].getMultipoleType();
            OpenMM_DoubleArray_append(doubleArray, multipoleType.charge);
            OpenMM_DoubleArray_append(doubleArray, OpenMM_NmPerAngstrom * baseRadii[i]);
            OpenMM_DoubleArray_append(doubleArray, overlapScale[i]);
            OpenMM_DoubleArray_append(doubleArray, sTens);

            OpenMM_CustomGBForce_addParticle(customGBForce, doubleArray);
            OpenMM_DoubleArray_resize(doubleArray, 0);
        }
        OpenMM_DoubleArray_destroy(doubleArray);

        double cut = gk.getCutoff();
        OpenMM_CustomGBForce_setCutoffDistance(customGBForce, cut);

        ForceField.ForceFieldInteger gbForceGroup = ForceField.ForceFieldInteger.GK_FORCE_GROUP;
        int fGroup = forceField.getInteger(gbForceGroup, gbForceGroup.getDefaultValue());

        OpenMM_Force_setForceGroup(customGBForce, fGroup);
        OpenMM_System_addForce(system, customGBForce);

        logger.log(Level.INFO, "  Generalized Born force");
    }

    private void addAmoebaVDWForce(ForceField forceField) {
        VanDerWaals vdW = super.getVdwNode();
        if (vdW == null) {
            return;
        }

        amoebaVDWForce = OpenMM_AmoebaVdwForce_create();
        OpenMM_System_addForce(system, amoebaVDWForce);

        ForceField.ForceFieldInteger vdwFgroup = ForceField.ForceFieldInteger.VDW_FORCE_GROUP;
        int fGroup = forceField.getInteger(vdwFgroup, vdwFgroup.getDefaultValue());
        OpenMM_Force_setForceGroup(amoebaVDWForce, fGroup);

        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
        Crystal crystal = super.getCrystal();

        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }

        /**
         * Note that the API says it wants a SIGMA value.
         */
        if (vdwForm.radiusType == VanDerWaalsForm.RADIUS_TYPE.R_MIN) {
            //radScale *= 1.122462048309372981;
        }

        int ired[] = vdW.getReductionIndex();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            OpenMM_AmoebaVdwForce_addParticle(amoebaVDWForce,
                    ired[i], OpenMM_NmPerAngstrom * vdwType.radius * radScale,
                    OpenMM_KJPerKcal * vdwType.wellDepth,
                    vdwType.reductionFactor);
        }

        // OpenMM_AmoebaVdwForce_setSigmaCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.radiusRule.name()));
        // OpenMM_AmoebaVdwForce_setEpsilonCombiningRule(amoebaVdwForce, toPropertyForm(vdwForm.epsilonRule.name()));
        OpenMM_AmoebaVdwForce_setCutoffDistance(amoebaVDWForce, nonbondedCutoff.off * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_Boolean.OpenMM_False);

        if (crystal.aperiodic()) {
            OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
                    OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_NoCutoff);
        } else {
            OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
                    OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_CutoffPeriodic);
        }

        /**
         * Create exclusion lists.
         */
        PointerByReference exclusions = OpenMM_IntArray_create(0);
        double mask[] = new double[nAtoms];
        Arrays.fill(mask, 1.0);
        for (int i = 0; i < nAtoms; i++) {
            OpenMM_IntArray_append(exclusions, i);
            vdW.applyMask(mask, i);
            for (int j = 0; j < nAtoms; j++) {
                if (mask[j] == 0.0) {
                    OpenMM_IntArray_append(exclusions, j);
                }
            }
            vdW.removeMask(mask, i);
            OpenMM_AmoebaVdwForce_setParticleExclusions(amoebaVDWForce, i, exclusions);
            OpenMM_IntArray_resize(exclusions, 0);
        }
        OpenMM_IntArray_destroy(exclusions);
        logger.log(Level.INFO, "  van der Waals force");
    }

    /**
     * Experimental. Virtual hydrogen sites require creation of new particles,
     * which then need to be handled (ignored?) for the multiple force.
     */
    private void createVirtualHydrogenSites() {

        VanDerWaals vdW = super.getVdwNode();
        if (vdW == null) {
            return;
        }
        int ired[] = vdW.getReductionIndex();

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            if (vdwType.reductionFactor < 1.0) {
                double factor = vdwType.reductionFactor;
                // Create the virtual site.
                PointerByReference virtualSite = OpenMM_TwoParticleAverageSite_create(i, ired[i], factor, 1.0 - factor);
                // Create a massless particle for the hydrogen vdW site.
                int id = OpenMM_System_addParticle(system, 0.0);
                // Denote the massless particle is a virtual site
                OpenMM_System_setVirtualSite(system, id, virtualSite);
            }
        }
    }

    private void addAmoebaMultipoleForce(ForceField forceField) {
        ParticleMeshEwald pme = super.getPmeNode();
        if (pme == null) {
            return;
        }
        int axisAtom[][] = pme.getAxisAtoms();
        double dipoleConversion = OpenMM_NmPerAngstrom;
        double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
        double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom
                * OpenMM_NmPerAngstrom;
        double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

        amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();

        double polarScale = 1.0;
        if (pme.getPolarizationType() != Polarization.MUTUAL) {
            OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Direct);
            if (pme.getPolarizationType() == Polarization.NONE) {
                polarScale = 0.0;
            }
        } else {
            String algorithm = forceField.getString(ForceField.ForceFieldString.SCF_ALGORITHM, "CG");
            ParticleMeshEwald.SCFAlgorithm scfAlgorithm;

            try {
                algorithm = algorithm.replaceAll("-", "_").toUpperCase();
                scfAlgorithm = ParticleMeshEwald.SCFAlgorithm.valueOf(algorithm);
            } catch (Exception e) {
                scfAlgorithm = ParticleMeshEwald.SCFAlgorithm.CG;
            }

            switch (scfAlgorithm) {
                case EPT:
                    logger.info(" Using extrapolated perturbation theory approximation instead of full SCF calculations. Not supported in FFX reference implementation.");
                    OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Extrapolated);
                    PointerByReference exptCoefficients = OpenMM_DoubleArray_create(4);
                    OpenMM_DoubleArray_set(exptCoefficients, 0, -0.154);
                    OpenMM_DoubleArray_set(exptCoefficients, 1, 0.017);
                    OpenMM_DoubleArray_set(exptCoefficients, 2, 0.657);
                    OpenMM_DoubleArray_set(exptCoefficients, 3, 0.475);
                    OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(amoebaMultipoleForce, exptCoefficients);
                    OpenMM_DoubleArray_destroy(exptCoefficients);
                    break;
                case CG:
                case SOR:
                default:
                    OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Mutual);
                    break;
            }
        }

        PointerByReference dipoles = OpenMM_DoubleArray_create(3);
        PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            MultipoleType multipoleType = atom.getMultipoleType();
            PolarizeType polarType = atom.getPolarizeType();

            /**
             * Define the frame definition.
             */
            int axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
            switch (multipoleType.frameDefinition) {
                case ZONLY:
                    axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
                    break;
                case ZTHENX:
                    axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
                    break;
                case BISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_Bisector;
                    break;
                case ZTHENBISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
                    break;
                case TRISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
                    break;
                default:
                    break;
            }

            double useFactor = 1.0;
            if (!atoms[i].getUse() || !atoms[i].getElectrostatics()) {
                //if (!atoms[i].getUse()) {
                useFactor = 0.0;
            }

            double lambdaScale = lambdaElec; // Should be 1.0 at this point.
            if (!atom.applyLambda()) {
                lambdaScale = 1.0;
            }

            useFactor *= lambdaScale;

            /**
             * Load local multipole coefficients.
             */
            for (int j = 0; j < 3; j++) {
                OpenMM_DoubleArray_set(dipoles, j, multipoleType.dipole[j] * dipoleConversion * useFactor);

            }
            int l = 0;
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    OpenMM_DoubleArray_set(quadrupoles, l++, multipoleType.quadrupole[j][k] * quadrupoleConversion * useFactor / 3.0);
                }
            }

            int zaxis = 1;
            int xaxis = 1;
            int yaxis = 1;
            int refAtoms[] = axisAtom[i];
            if (refAtoms != null) {
                zaxis = refAtoms[0];
                if (refAtoms.length > 1) {
                    xaxis = refAtoms[1];
                    if (refAtoms.length > 2) {
                        yaxis = refAtoms[2];
                    }
                }
            } else {
                axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
                logger.info(String.format(" Atom type %s", atom.getAtomType().toString()));
            }

            double charge = multipoleType.charge * useFactor;

            /**
             * Add the multipole.
             */
            OpenMM_AmoebaMultipoleForce_addMultipole(amoebaMultipoleForce,
                    charge, dipoles, quadrupoles,
                    axisType, zaxis, xaxis, yaxis,
                    polarType.thole,
                    polarType.pdamp * dampingFactorConversion,
                    polarType.polarizability * polarityConversion * polarScale);
        }
        OpenMM_DoubleArray_destroy(dipoles);
        OpenMM_DoubleArray_destroy(quadrupoles);

        Crystal crystal = super.getCrystal();
        if (!crystal.aperiodic()) {
            OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_PME);
            OpenMM_AmoebaMultipoleForce_setCutoffDistance(amoebaMultipoleForce,
                    pme.getEwaldCutoff() * OpenMM_NmPerAngstrom);
            OpenMM_AmoebaMultipoleForce_setAEwald(amoebaMultipoleForce,
                    pme.getEwaldCoefficient() / OpenMM_NmPerAngstrom);

            double ewaldTolerance = 1.0e-04;
            OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(amoebaMultipoleForce, ewaldTolerance);

            PointerByReference gridDimensions = OpenMM_IntArray_create(3);
            ReciprocalSpace recip = pme.getReciprocalSpace();
            OpenMM_IntArray_set(gridDimensions, 0, recip.getXDim());
            OpenMM_IntArray_set(gridDimensions, 1, recip.getYDim());
            OpenMM_IntArray_set(gridDimensions, 2, recip.getZDim());
            OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(amoebaMultipoleForce, gridDimensions);
            OpenMM_IntArray_destroy(gridDimensions);
        } else {
            OpenMM_AmoebaMultipoleForce_setNonbondedMethod(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_NoCutoff);
        }

        OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(amoebaMultipoleForce, 500);
        OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(amoebaMultipoleForce, pme.getPolarEps());

        int ip11[][] = pme.getPolarization11();
        int ip12[][] = pme.getPolarization12();
        int ip13[][] = pme.getPolarization13();

        ArrayList<Integer> list12 = new ArrayList<>();
        ArrayList<Integer> list13 = new ArrayList<>();
        ArrayList<Integer> list14 = new ArrayList<>();

        PointerByReference covalentMap = OpenMM_IntArray_create(0);
        for (int i = 0; i < nAtoms; i++) {
            Atom ai = atoms[i];
            list12.clear();
            list13.clear();
            list14.clear();

            for (Bond bond : ai.getBonds()) {
                int index = bond.get1_2(ai).getIndex() - 1;
                OpenMM_IntArray_append(covalentMap, index);
                list12.add(index);
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent12, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Angle angle : ai.getAngles()) {
                Atom ak = angle.get1_3(ai);
                if (ak != null) {
                    int index = ak.getIndex() - 1;
                    if (!list12.contains(index)) {
                        list13.add(index);
                        OpenMM_IntArray_append(covalentMap, index);
                    }
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent13, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Torsion torsion : ai.getTorsions()) {
                Atom ak = torsion.get1_4(ai);
                if (ak != null) {
                    int index = ak.getIndex() - 1;
                    if (!list12.contains(index)
                            && !list13.contains(index)) {
                        list14.add(index);
                        OpenMM_IntArray_append(covalentMap, index);
                    }
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent14, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (Atom ak : ai.get1_5s()) {
                int index = ak.getIndex() - 1;
                if (!list12.contains(index)
                        && !list13.contains(index)
                        && !list14.contains(index)) {
                    OpenMM_IntArray_append(covalentMap, index);
                }
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_Covalent15, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

            for (int j = 0; j < ip11[i].length; j++) {
                OpenMM_IntArray_append(covalentMap, ip11[i][j]);
            }
            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent11, covalentMap);
            OpenMM_IntArray_resize(covalentMap, 0);

//            for (int j = 0; j < ip12[i].length; j++) {
//                OpenMM_IntArray_append(covalentMap, ip12[i][j]);
//            }
//            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
//                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent12, covalentMap);
//            OpenMM_IntArray_resize(covalentMap, 0);
//
//            for (int j = 0; j < ip13[i].length; j++) {
//                OpenMM_IntArray_append(covalentMap, ip13[i][j]);
//            }
//            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
//                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent13, covalentMap);
//            OpenMM_IntArray_resize(covalentMap, 0);
//
//            OpenMM_AmoebaMultipoleForce_setCovalentMap(amoebaMultipoleForce, i,
//                    OpenMM_AmoebaMultipoleForce_PolarizationCovalent14, covalentMap);
        }

        OpenMM_IntArray_destroy(covalentMap);

        ForceField.ForceFieldInteger pmeFgroup = ForceField.ForceFieldInteger.PME_FORCE_GROUP;
        int fGroup = forceField.getInteger(pmeFgroup, pmeFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaMultipoleForce, fGroup);
        OpenMM_System_addForce(system, amoebaMultipoleForce);
        logger.log(Level.INFO, "  Polarizable multipole force");

        GeneralizedKirkwood gk = super.getGK();
        if (gk != null) {
            addGKForce(forceField);
        }

    }

    private void addGKForce(ForceField forceField) {

        GeneralizedKirkwood gk = super.getGK();

        amoebaGeneralizedKirkwoodForce = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
        OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(amoebaGeneralizedKirkwoodForce, gk.getSolventPermittivity());
        OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(amoebaGeneralizedKirkwoodForce, 1.0);

        double overlapScale[] = gk.getOverlapScale();
        double baseRadii[] = gk.getBaseRadii();
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            MultipoleType multipoleType = atoms[i].getMultipoleType();
            OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(amoebaGeneralizedKirkwoodForce,
                    multipoleType.charge, OpenMM_NmPerAngstrom * baseRadii[i], overlapScale[i]);
        }

        OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(amoebaGeneralizedKirkwoodForce, gk.getProbeRadius() * OpenMM_NmPerAngstrom);

        NonPolar nonpolar = gk.getNonPolarModel();
        switch (nonpolar) {
            case BORN_SOLV:
            case BORN_CAV_DISP:
            default:
                // Configure a Born Radii based surface area term.
                double surfaceTension = gk.getSurfaceTension() * OpenMM_KJPerKcal
                        * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm;
                OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(amoebaGeneralizedKirkwoodForce, OpenMM_True);
                OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(amoebaGeneralizedKirkwoodForce, -surfaceTension);
                break;
            case CAV:
            case CAV_DISP:
            case HYDROPHOBIC_PMF:
            case NONE:
                // This NonPolar model does not use a Born Radii based surface area term.
                OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(amoebaGeneralizedKirkwoodForce, OpenMM_False);
                break;
        }

        ForceField.ForceFieldInteger gkFgroup = ForceField.ForceFieldInteger.GK_FORCE_GROUP;
        int fGroup = forceField.getInteger(gkFgroup, gkFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaGeneralizedKirkwoodForce, fGroup);
        OpenMM_System_addForce(system, amoebaGeneralizedKirkwoodForce);

        switch (nonpolar) {
            case CAV_DISP:
            case BORN_CAV_DISP:
                addWCAForce(forceField);
                break;
            case CAV:
            case HYDROPHOBIC_PMF:
            case BORN_SOLV:
            case NONE:
            default:
                // WCA force is not being used.
        }

        logger.log(Level.INFO, "  Generalized Kirkwood force");
    }

    private void addWCAForce(ForceField forceField) {

        double epso = 0.1100;
        double epsh = 0.0135;
        double rmino = 1.7025;
        double rminh = 1.3275;
        double awater = 0.033428;
        double slevy = 1.0;
        double dispoff = 0.26;
        double shctd = 0.81;

        VanDerWaals vdW = super.getVdwNode();
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }

        amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create();

        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;

        for (int i = 0; i < nAtoms; i++) {
            // cdispTotal += nonpol__.cdisp[ii];
            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            double radius = vdwType.radius;
            double eps = vdwType.wellDepth;
            OpenMM_AmoebaWcaDispersionForce_addParticle(amoebaWcaDispersionForce,
                    OpenMM_NmPerAngstrom * radius * radScale,
                    OpenMM_KJPerKcal * eps);
        }

        OpenMM_AmoebaWcaDispersionForce_setEpso(amoebaWcaDispersionForce, epso * OpenMM_KJPerKcal);
        OpenMM_AmoebaWcaDispersionForce_setEpsh(amoebaWcaDispersionForce, epsh * OpenMM_KJPerKcal);
        OpenMM_AmoebaWcaDispersionForce_setRmino(amoebaWcaDispersionForce, rmino * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaWcaDispersionForce_setRminh(amoebaWcaDispersionForce, rminh * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaWcaDispersionForce_setDispoff(amoebaWcaDispersionForce, dispoff * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaWcaDispersionForce_setAwater(amoebaWcaDispersionForce,
                awater / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
        OpenMM_AmoebaWcaDispersionForce_setSlevy(amoebaWcaDispersionForce, slevy);
        OpenMM_AmoebaWcaDispersionForce_setShctd(amoebaWcaDispersionForce, shctd);

        ForceField.ForceFieldInteger gkFgroup = ForceField.ForceFieldInteger.GK_FORCE_GROUP;
        int fGroup = forceField.getInteger(gkFgroup, gkFgroup.getDefaultValue());

        OpenMM_Force_setForceGroup(amoebaWcaDispersionForce, fGroup);
        OpenMM_System_addForce(system, amoebaWcaDispersionForce);
        logger.log(Level.INFO, "  WCA dispersion force.");

    }

    /**
     * Adds harmonic restraints (CoordRestraint objects) to OpenMM as a custom
     * external force.
     *
     * TODO: Make robust to flat-bottom restraints.
     */
    private void addHarmonicRestraintForce(ForceField forceField) {
        for (CoordRestraint restraint : super.getCoordRestraints()) {
            double forceConst = restraint.getForceConstant();
            forceConst *= OpenMM_KJPerKcal;
            forceConst *= (OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm);
            Atom[] restAtoms = restraint.getAtoms();
            int nRestAts = restraint.getNumAtoms();
            double[][] oCoords = restraint.getOriginalCoordinates();
            for (int i = 0; i < nRestAts; i++) {
                oCoords[i][0] *= OpenMM_NmPerAngstrom;
                oCoords[i][1] *= OpenMM_NmPerAngstrom;
                oCoords[i][2] *= OpenMM_NmPerAngstrom;
            }

            PointerByReference theRestraint = OpenMM_CustomExternalForce_create("k*periodicdistance(x,y,z,x0,y0,z0)^2");
            OpenMM_CustomExternalForce_addGlobalParameter(theRestraint, "k", forceConst);
            OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "x0");
            OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "y0");
            OpenMM_CustomExternalForce_addPerParticleParameter(theRestraint, "z0");

            PointerByReference xyzOrigArray = OpenMM_DoubleArray_create(3);
            for (int i = 0; i < nRestAts; i++) {
                int ommIndex = restAtoms[i].getXyzIndex() - 1;
                for (int j = 0; j < 3; j++) {
                    OpenMM_DoubleArray_set(xyzOrigArray, j, oCoords[i][j]);
                }
                OpenMM_CustomExternalForce_addParticle(theRestraint, ommIndex, xyzOrigArray);
            }
            OpenMM_DoubleArray_destroy(xyzOrigArray);

            // TODO: Be able to set this on a per-restraint basis.
            ForceField.ForceFieldInteger xyzRestFgroup = ForceField.ForceFieldInteger.COORD_RESTRAINT_FORCE_GROUP;
            int fGroup = forceField.getInteger(xyzRestFgroup, xyzRestFgroup.getDefaultValue());

            OpenMM_Force_setForceGroup(theRestraint, fGroup);
            OpenMM_System_addForce(system, theRestraint);
        }
    }

    /**
     * Adds restraint bonds, if any.
     */
    private void addRestraintBonds(ForceField forceField) {
        List<RestraintBond> restraintBonds = super.getRestraintBonds();

        if (restraintBonds != null && !restraintBonds.isEmpty()) {
            // OpenMM's HarmonicBondForce class uses k, not 1/2*k as does FFX.
            double kParameterConversion = BondType.units * 2.0 * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

            for (RestraintBond rbond : super.getRestraintBonds()) {
                PointerByReference theForce; // Is not a valid substitute for a targeting computer.
                BondType bType = rbond.bondType;
                BondType.BondFunction funct = bType.bondFunction;
                if (restraintForces.containsKey(funct)) {
                    theForce = restraintForces.get(funct);
                } else {
                    theForce = OpenMM_CustomBondForce_create(funct.toMathematicalForm());
                    OpenMM_CustomBondForce_addPerBondParameter(theForce, "k");
                    OpenMM_CustomBondForce_addPerBondParameter(theForce, "r0");
                    if (funct.hasFlatBottom()) {
                        OpenMM_CustomBondForce_addPerBondParameter(theForce, "fb");
                    }

                    // Wholly untested code.
                    switch (funct) {
                        case QUARTIC:
                        case FLAT_BOTTOM_QUARTIC:
                            OpenMM_CustomBondForce_addGlobalParameter(theForce, "cubic",
                                    BondType.cubic / OpenMM_NmPerAngstrom);
                            OpenMM_CustomBondForce_addGlobalParameter(theForce, "quartic",
                                    BondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
                            break;
                        default:
                            break;
                    }

                    // TODO: Set this on a per-restraint basis.
                    ForceField.ForceFieldInteger bondRestFgroup = ForceField.ForceFieldInteger.BOND_RESTRAINT_FORCE_GROUP;
                    int fGroup = forceField.getInteger(bondRestFgroup, bondRestFgroup.getDefaultValue());

                    OpenMM_Force_setForceGroup(theForce, fGroup);
                    OpenMM_System_addForce(system, theForce);
                }
                double forceConst = bType.forceConstant * kParameterConversion;
                double equilDist = bType.distance * OpenMM_NmPerAngstrom;
                Atom[] ats = rbond.getAtomArray();
                int at1 = ats[0].getXyzIndex() - 1;
                int at2 = ats[1].getXyzIndex() - 1;

                PointerByReference bondParams = OpenMM_DoubleArray_create(0);
                OpenMM_DoubleArray_append(bondParams, forceConst);
                OpenMM_DoubleArray_append(bondParams, equilDist);
                if (funct.hasFlatBottom()) {
                    OpenMM_DoubleArray_append(bondParams, bType.flatBottomRadius * OpenMM_NmPerAngstrom);
                }
                OpenMM_CustomBondForce_addBond(theForce, at1, at2, bondParams);
                OpenMM_DoubleArray_destroy(bondParams);
            }
        }
    }

    /**
     * Update parameters if the Use flags changed.
     */
    private void updateParameters(double x[]) {

        Atom[] atoms = molecularAssembly.getAtomArray();

        if (vdwLambdaTerm) {
            if (amoebaVDWForce != null) {
                logger.severe(" Softcore vdW is not yet supported for AMOEBA.");
            }
            if (!softcoreCreated) {
                ForceField forceField = molecularAssembly.getForceField();
                addCustomNonbondedSoftcoreForce(forceField);
                // Reset the context.
                createContext(integratorString, timeStep, temperature);
                OpenMM_Context_setParameter(context, "vdw_lambda", lambdaVDW);
                softcoreCreated = true;
                if (x != null) {
                    double energy = energy(x);
                    logger.info(format(" OpenMM Energy (L=%6.3f): %16.8f", lambdaVDW, energy));
                }
            } else {
                OpenMM_Context_setParameter(context, "vdw_lambda", lambdaVDW);
            }
        }

        if (torsionLambdaTerm && amoebaTorsionForce != null) {
            updateTorsionForce(atoms);
        }

        // Update fixed charge non-bonded parameters.
        if (fixedChargeNonBondedForce != null) {
            updateFixedChargeNonBondedForce(atoms);
        }

        // Update fixed charge GB parameters.
        if (customGBForce != null) {
            updateCustomGBForce(atoms);
        }

        // Update AMOEBA vdW parameters.
        if (amoebaVDWForce != null) {
            updateAmoebaVDWForce(atoms);
        }

        // Update AMOEBA polarizable multipole parameters.
        if (amoebaMultipoleForce != null) {
            updateAmoebaMultipoleForce(atoms);
        }

        // Update GK force.
        if (amoebaGeneralizedKirkwoodForce != null) {
            updateAmoebaGeneralizedKirkwoodForce(atoms);
        }

        // Update WCA Force.
        if (amoebaWcaDispersionForce != null) {
            updateWCAForce(atoms);
        }
    }

    /**
     * Updates the AMOEBA van der Waals force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     */
    private void updateAmoebaVDWForce(Atom[] atoms) {
        VanDerWaals vdW = super.getVdwNode();
        VanDerWaalsForm vdwForm = vdW.getVDWForm();

        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }

        /**
         * Note that the API says it wants a SIGMA value.
         */
        if (vdwForm.radiusType == VanDerWaalsForm.RADIUS_TYPE.R_MIN) {
            //radScale *= 1.122462048309372981;
        }

        int ired[] = vdW.getReductionIndex();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            double useFactor = 1.0;
            if (!atoms[i].getUse()) {
                useFactor = 0.0;
            }
            double eps = OpenMM_KJPerKcal * vdwType.wellDepth * useFactor;
            OpenMM_AmoebaVdwForce_setParticleParameters(amoebaVDWForce,
                    i, ired[i], OpenMM_NmPerAngstrom * vdwType.radius * radScale,
                    eps, vdwType.reductionFactor);
        }
        OpenMM_AmoebaVdwForce_updateParametersInContext(amoebaVDWForce, context);
    }

    /**
     * Updates the fixed-charge non-bonded force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     */
    private void updateFixedChargeNonBondedForce(Atom[] atoms) {
        VanDerWaals vdW = super.getVdwNode();
        /**
         * Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
         * for epsilon is supported.
         */
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        if (vdwForm.vdwType != LENNARD_JONES
                || vdwForm.radiusRule != ARITHMETIC
                || vdwForm.epsilonRule != GEOMETRIC) {
            logger.log(Level.SEVERE, String.format(" Unsupported van der Waals functional form."));
            return;
        }

        /**
         * OpenMM vdW force requires a diameter (i.e. not radius).
         */
        double radScale = 1.0;
        if (vdwForm.radiusSize == RADIUS) {
            radScale = 2.0;
        }
        /**
         * OpenMM vdw force requires atomic sigma values (i.e. not r-min).
         */
        if (vdwForm.radiusType == R_MIN) {
            radScale /= 1.122462048309372981;
        }

        /**
         * Update parameters.
         */
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            boolean applyLambda = atom.applyLambda();

            double charge = Double.MIN_VALUE;
            MultipoleType multipoleType = atom.getMultipoleType();
            if (multipoleType != null && atom.getElectrostatics()) {
                charge = multipoleType.charge;
            }

            VDWType vdwType = atom.getVDWType();
            double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
            double eps = OpenMM_KJPerKcal * vdwType.wellDepth;

            if (applyLambda) {
                // If we're using vdwLambdaTerm, this atom's vdW interactions are handled by the Custom Non-Bonded force.
                if (vdwLambdaTerm) {
                    eps = 0.0;
                }
                // Always scale the charge by lambdaElec
                charge *= lambdaElec;
                // logger.info(format(" %s using electrostatics %b, with charge %8.3f", atom.toString(), atom.getElectrostatics(), charge));
            }

            if (!atoms[i].getUse()) {
                eps = 0.0;
                charge = 0.0;
            }

            OpenMM_NonbondedForce_setParticleParameters(fixedChargeNonBondedForce, i, charge, sigma, eps);
        }

        /**
         * Update Exceptions.
         */
        IntByReference particle1 = new IntByReference();
        IntByReference particle2 = new IntByReference();
        DoubleByReference chargeProd = new DoubleByReference();
        DoubleByReference sigma = new DoubleByReference();
        DoubleByReference eps = new DoubleByReference();

        int numExceptions = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);

        for (int i = 0; i < numExceptions; i++) {

            /**
             * Only update exceptions.
             */
            if (chargeExclusion[i] && vdWExclusion[i]) {
                continue;
            }

            OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i,
                    particle1, particle2, chargeProd, sigma, eps);

            int i1 = particle1.getValue();
            int i2 = particle2.getValue();

            double qq = exceptionChargeProd[i];
            double epsilon = exceptionEps[i];

            Atom atom1 = atoms[i1];
            Atom atom2 = atoms[i2];

            double lambdaValue = lambdaElec;
            if (lambdaValue == 0.0) {
                lambdaValue = 1.0e-6;
            }

            if (atom1.applyLambda()) {
                qq *= lambdaValue;
                if (vdwLambdaTerm) {
                    epsilon = 1.0e-6;
                    // qq = 1.0e-6;
                }
            }

            if (atom2.applyLambda()) {
                qq *= lambdaValue;
                if (vdwLambdaTerm) {
                    epsilon = 1.0e-6;
                    // qq = 1.0e-6;
                }
            }

            if (!atom1.getUse() || !atom2.getUse()) {
                qq = 1.0e-6;
                epsilon = 1.0e-6;
            }

            OpenMM_NonbondedForce_setExceptionParameters(fixedChargeNonBondedForce, i,
                    i1, i2, qq, sigma.getValue(), epsilon);

            /**
             * logger.info(format(" B Exception %d %d %d q=%10.8f s=%10.8f
             * e=%10.8f.", i, i1, i2, chargeProd.getValue(), sigma.getValue(),
             * eps.getValue()));
             *
             * logger.info(format(" E Exception %d %d %d q=%10.8f s=%10.8f
             * e=%10.8f.", i, i1, i2, qq, sigma.getValue(), epsilon));
             */
        }

        OpenMM_NonbondedForce_updateParametersInContext(fixedChargeNonBondedForce, context);
    }

    /**
     * Updates the custom GB force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     */
    private void updateCustomGBForce(Atom[] atoms) {
        GeneralizedKirkwood gk = super.getGK();
        double[] baseRadii = gk.getBaseRadii();
        double[] overlapScale = gk.getOverlapScale();
        PointerByReference doubleArray = OpenMM_DoubleArray_create(0);
        boolean nea = gk.getNativeEnvironmentApproximation();

        double sTens = 0.0;
        if (gk.getNonPolarModel() == NonPolar.BORN_SOLV || gk.getNonPolarModel() == NonPolar.BORN_CAV_DISP) {
            sTens = gk.getSurfaceTension();
            sTens *= OpenMM_KJPerKcal;
            sTens *= 100.0; // 100 square Angstroms per square nanometer.
        }

        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            double chargeUseFactor = 1.0;
            if (!atoms[i].getUse() || !atoms[i].getElectrostatics()) {
                chargeUseFactor = 0.0;
            }
            double lambdaScale = lambdaElec;
            if (!atom.applyLambda()) {
                lambdaScale = 1.0;
            }

            chargeUseFactor *= lambdaScale;
            double overlapScaleUseFactor = nea ? 1.0 : chargeUseFactor;
            double oScale = overlapScale[i] * overlapScaleUseFactor;

            MultipoleType multipoleType = atom.getMultipoleType();
            double charge = multipoleType.charge * chargeUseFactor;
            double surfaceTension = sTens * chargeUseFactor;

            double baseRadius = baseRadii[i];

            OpenMM_DoubleArray_append(doubleArray, charge);
            OpenMM_DoubleArray_append(doubleArray, OpenMM_NmPerAngstrom * baseRadius);
            OpenMM_DoubleArray_append(doubleArray, oScale);
            OpenMM_DoubleArray_append(doubleArray, surfaceTension);

            OpenMM_CustomGBForce_setParticleParameters(customGBForce, i, doubleArray);
            OpenMM_DoubleArray_resize(doubleArray, 0);
        }
        OpenMM_DoubleArray_destroy(doubleArray);
        OpenMM_CustomGBForce_updateParametersInContext(customGBForce, context);
    }

    /**
     * Updates the Amoeba electrostatic multipolar force for change in Use
     * flags.
     *
     * @param atoms Array of all Atoms in the system
     */
    private void updateAmoebaMultipoleForce(Atom[] atoms) {
        ParticleMeshEwald pme = super.getPmeNode();
        int axisAtom[][] = pme.getAxisAtoms();
        double dipoleConversion = OpenMM_NmPerAngstrom;
        double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
        double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom
                * OpenMM_NmPerAngstrom;
        double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

        double polarScale = 1.0;
        if (pme.getPolarizationType() == Polarization.NONE) {
            polarScale = 0.0;
        }

        PointerByReference dipoles = OpenMM_DoubleArray_create(3);
        PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom atom = atoms[i];
            MultipoleType multipoleType = atom.getMultipoleType();
            PolarizeType polarType = atom.getPolarizeType();
            double useFactor = 1.0;

            if (!atoms[i].getUse() || !atoms[i].getElectrostatics()) {
                //if (!atoms[i].getUse()) {
                useFactor = 0.0;
            }

            double lambdaScale = lambdaElec;
            if (!atom.applyLambda()) {
                lambdaScale = 1.0;
            }

            useFactor *= lambdaScale;

            /**
             * Define the frame definition.
             */
            int axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
            switch (multipoleType.frameDefinition) {
                case ZONLY:
                    axisType = OpenMM_AmoebaMultipoleForce_ZOnly;
                    break;
                case ZTHENX:
                    axisType = OpenMM_AmoebaMultipoleForce_ZThenX;
                    break;
                case BISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_Bisector;
                    break;
                case ZTHENBISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ZBisect;
                    break;
                case TRISECTOR:
                    axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
                    break;
                default:
                    break;
            }

            /**
             * Load local multipole coefficients.
             */
            for (int j = 0; j < 3; j++) {
                OpenMM_DoubleArray_set(dipoles, j, multipoleType.dipole[j] * dipoleConversion * useFactor);

            }
            int l = 0;
            for (int j = 0; j < 3; j++) {
                for (int k = 0; k < 3; k++) {
                    OpenMM_DoubleArray_set(quadrupoles, l++, multipoleType.quadrupole[j][k]
                            * quadrupoleConversion / 3.0 * useFactor);
                }
            }

            //int zaxis = 0;
            int zaxis = 1;
            //int xaxis = 0;
            int xaxis = 1;
            //int yaxis = 0;
            int yaxis = 1;
            int refAtoms[] = axisAtom[i];
            if (refAtoms != null) {
                zaxis = refAtoms[0];
                if (refAtoms.length > 1) {
                    xaxis = refAtoms[1];
                    if (refAtoms.length > 2) {
                        yaxis = refAtoms[2];
                    }
                }
            } else {
                axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
            }

            /**
             * Add the multipole.
             */
            OpenMM_AmoebaMultipoleForce_setMultipoleParameters(amoebaMultipoleForce, i,
                    multipoleType.charge * useFactor, dipoles, quadrupoles,
                    axisType, zaxis, xaxis, yaxis,
                    polarType.thole,
                    polarType.pdamp * dampingFactorConversion,
                    polarType.polarizability * polarityConversion * polarScale * useFactor);
        }
        OpenMM_DoubleArray_destroy(dipoles);
        OpenMM_DoubleArray_destroy(quadrupoles);

        OpenMM_AmoebaMultipoleForce_updateParametersInContext(amoebaMultipoleForce, context);
    }

    /**
     * Updates the AMOEBA Generalized Kirkwood force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system
     */
    private void updateAmoebaGeneralizedKirkwoodForce(Atom[] atoms) {
        GeneralizedKirkwood gk = super.getGK();
        double overlapScale[] = gk.getOverlapScale();
        double baseRadii[] = gk.getBaseRadii();
        int nAtoms = atoms.length;
        boolean nea = gk.getNativeEnvironmentApproximation();

        for (int i = 0; i < nAtoms; i++) {
            double chargeUseFactor = 1.0;
            if (!atoms[i].getUse() || !atoms[i].getElectrostatics()) {
                chargeUseFactor = 0.0;
            }

            double lambdaScale = lambdaElec;
            if (!atoms[i].applyLambda()) {
                lambdaScale = 1.0;
            }

            chargeUseFactor *= lambdaScale;
            double overlapScaleUseFactor = nea ? 1.0 : chargeUseFactor;

            MultipoleType multipoleType = atoms[i].getMultipoleType();
            OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(amoebaGeneralizedKirkwoodForce, i,
                    multipoleType.charge * chargeUseFactor,
                    OpenMM_NmPerAngstrom * baseRadii[i], overlapScale[i] * overlapScaleUseFactor);
        }
        OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(amoebaGeneralizedKirkwoodForce, context);
    }

    /**
     * Updates the WCA force for change in Use flags.
     *
     * @param atoms Array of all Atoms in the system.
     */
    private void updateWCAForce(Atom[] atoms) {
        VanDerWaals vdW = super.getVdwNode();
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        double radScale = 1.0;
        if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
            radScale = 0.5;
        }
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            double useFactor = 1.0;
            if (!atoms[i].getUse()) {
                useFactor = 0.0;
            }

            // Scale all implicit solvent terms with the lambda for electrostatics.
            double lambdaScale = lambdaElec;
            if (!atoms[i].applyLambda()) {
                lambdaScale = 1.0;
            }
            useFactor *= lambdaScale;

            Atom atom = atoms[i];
            VDWType vdwType = atom.getVDWType();
            double radius = vdwType.radius;
            double eps = vdwType.wellDepth;
            OpenMM_AmoebaWcaDispersionForce_setParticleParameters(amoebaWcaDispersionForce, i,
                    OpenMM_NmPerAngstrom * radius * radScale,
                    OpenMM_KJPerKcal * eps * useFactor);
        }
        OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(amoebaWcaDispersionForce, context);
    }

    /**
     * Updates the Torsion force for application of lambda scaling.
     *
     * @param atoms Array of all Atoms in the system.
     */
    private void updateTorsionForce(Atom[] atoms) {

        // Only update parameters if torsions are being scaled by lambda.
        if (!torsionLambdaTerm) {
            return;
        }

        // Check if this system has torsions.
        Torsion torsions[] = super.getTorsions();
        if (torsions == null || torsions.length < 1) {
            return;
        }

        int nTorsions = torsions.length;
        int index = 0;
        for (int i = 0; i < nTorsions; i++) {
            Torsion torsion = torsions[i];
            int a1 = torsion.getAtom(0).getXyzIndex() - 1;
            int a2 = torsion.getAtom(1).getXyzIndex() - 1;
            int a3 = torsion.getAtom(2).getXyzIndex() - 1;
            int a4 = torsion.getAtom(3).getXyzIndex() - 1;
            TorsionType torsionType = torsion.torsionType;
            int nTerms = torsionType.phase.length;
            for (int j = 0; j < nTerms; j++) {

                double forceConstant = OpenMM_KJPerKcal * torsion.units * torsionType.amplitude[j];
                forceConstant *= lambdaTorsion;

                OpenMM_PeriodicTorsionForce_setTorsionParameters(
                        amoebaTorsionForce, index,
                        a1, a2, a3, a4, j + 1,
                        torsionType.phase[j] * OpenMM_RadiansPerDegree,
                        forceConstant);
                index++;
            }
        }

        OpenMM_PeriodicTorsionForce_updateParametersInContext(amoebaTorsionForce, context);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (lambdaTerm) {
            if (lambda >= 0.0 && lambda <= 1.0) {
                this.lambda = lambda;
                super.setLambda(lambda);

                // Initially set all lambda values to 1.0.
                lambdaTorsion = 1.0;
                lambdaVDW = 1.0;
                lambdaElec = 1.0;

                if (torsionLambdaTerm) {
                    lambdaTorsion = lambda;
                }


                if (elecLambdaTerm && vdwLambdaTerm) {
                    // Lambda effects both vdW and electrostatics.
                    lambdaVDW = lambda;
                    if (lambda < 0.5) {
                        // Begin turning vdW on with electrostatics off.
                        lambdaElec = 0.0;
                    } else {
                        // Turn electrostatics on during the latter part of the path.
                        lambdaElec = 2.0 * (lambda - 0.5);
                    }
                } else if (vdwLambdaTerm) {
                    // Lambda effects vdW, with electrostatics turned off.
                    lambdaVDW = lambda;
                    lambdaElec = 0.0;
                } else if (elecLambdaTerm) {
                    // Lambda effects electrostatics, but not vdW.
                    lambdaElec = lambda;
                }

                updateParameters(null);
            } else {
                String message = format(" Lambda value %8.3f is not in the range [0..1].", lambda);
                logger.warning(message);
            }
        } else {
            logger.fine(" Attempting to set a lambda value on a ForceFieldEnergyOpenMM with lambdaterm false.");
        }
    }

    /**
     * Evaluates energy both with OpenMM and reference potential, and returns
     * the difference FFX-OpenMM.
     *
     * @param x       Coordinate array
     * @param verbose a boolean.
     * @return Energy discrepancy
     */
    public double energyVsFFX(double[] x, boolean verbose) {
        double ffxE = super.energy(x, verbose);
        double thisE = energy(x, verbose);
        return ffxE - thisE;
    }

    /**
     * Evaluates energy explicitly using the Java implementation backing this
     * ForceFieldEnergyOpenMM.
     *
     * @param x       Coordinate array
     * @param verbose Verbosity of energy call
     * @return Total energy calculated by reference FFX implementation.
     */
    public double ffxEnergy(double[] x, boolean verbose) {
        return super.energy(x, verbose);
    }

    /**
     * Evaluates energy and gradients both with OpenMM and reference potential,
     * and returns the difference FFX-OpenMM.
     *
     * @param x       Coordinate array
     * @param gFFX    Array for FFX gradients to be stored in
     * @param gOMM    Array for OpenMM gradients to be stored in
     * @param verbose a boolean.
     * @return Energy discrepancy
     */
    public double energyAndGradVsFFX(double[] x, double[] gFFX, double[] gOMM, boolean verbose) {
        double ffxE = super.energyAndGradient(x, gFFX, verbose);
        double thisE = energyAndGradient(x, gOMM, verbose);
        return ffxE - thisE;
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

        if (lambdaBondedTerms) {
            return 0.0;
        }

        updateParameters(x);

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
        setOpenMMPositions(x, x.length);

        int infoMask = OpenMM_State_Energy;
        state = OpenMM_Context_getState(context, infoMask, enforcePBC);
        double e = OpenMM_State_getPotentialEnergy(state) / OpenMM_KJPerKcal;
        if (!Double.isFinite(e)) {
            String message = String.format(" Energy from OpenMM was a non-finite %8g", e);
            logger.warning(message);
            throw new EnergyException(message);
        }
        OpenMM_State_destroy(state);

        if (verbose) {
            logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", e));
        }

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
        if (lambdaBondedTerms) {
            return 0.0;
        }

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
        setOpenMMPositions(x, x.length);

        int infoMask = OpenMM_State_Energy;
        infoMask += OpenMM_State_Forces;

        if (vdwLambdaTerm) {
            infoMask += OpenMM_State_ParameterDerivatives;
        }

        state = OpenMM_Context_getState(context, infoMask, enforcePBC);
        double e = OpenMM_State_getPotentialEnergy(state) / OpenMM_KJPerKcal;
        if (!Double.isFinite(e)) {
            String message = String.format(" Energy from OpenMM was a non-finite %8g", e);
            logger.warning(message);
            throw new EnergyException(message);
        }

        if (vdwLambdaTerm) {
            PointerByReference parameterArray = OpenMM_State_getEnergyParameterDerivatives(state);
            int numDerives = OpenMM_ParameterArray_getSize(parameterArray);
            if (numDerives > 0) {
                vdwdUdL = OpenMM_ParameterArray_get(parameterArray, pointerForString("vdw_lambda")) / OpenMM_KJPerKcal;
            }
        }

        if (maxDebugGradient < Double.MAX_VALUE) {
            boolean extremeGrad = Arrays.stream(g).anyMatch((double gi) -> {
                return (gi > maxDebugGradient || gi < -maxDebugGradient);
            });
            if (extremeGrad) {
                File origFile = molecularAssembly.getFile();
                String timeString = LocalDateTime.now().format(DateTimeFormatter.
                        ofPattern("yyyy_MM_dd-HH_mm_ss"));

                String filename = String.format("%s-LARGEGRAD-%s.pdb",
                        FilenameUtils.removeExtension(molecularAssembly.getFile().getName()),
                        timeString);
                PotentialsFunctions ef = new PotentialsUtils();
                filename = ef.versionFile(filename);

                logger.warning(String.format(" Excessively large gradients detected; printing snapshot to file %s", filename));
                ef.saveAsPDB(molecularAssembly, new File(filename));
                molecularAssembly.setFile(origFile);
            }
        }

        if (verbose) {
            logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", e));
        }

        forces = OpenMM_State_getForces(state);

        fillGradients(g);
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

        OpenMM_State_destroy(state);
        return e;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setCrystal(Crystal crystal) {
        super.setCrystal(crystal);
        setDefaultPeriodicBoxVectors();
        //loadFFXPositionToOpenMM();
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * getGradients</p>
     */
    @Override
    public double[] getGradients(double g[]) {
        return fillGradients(g);
    }

    /**
     * Private method for internal use, so we don't have subclasses calling
     * super.energy, and this class delegating to the subclass's getGradients
     * method.
     *
     * @param g Gradient array to fill.
     * @return Gradient array.
     */
    public double[] fillGradients(double[] g) {
        assert (g != null);
        int n = getNumberOfVariables();
        if (g.length < n) {
            g = new double[n];
        }
        int index = 0;
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                OpenMM_Vec3 posInNm = OpenMM_Vec3Array_get(forces, i);
                /**
                 * Convert OpenMM Forces in KJ/Nm into an FFX gradient in
                 * Kcal/A.
                 */
                double gx = -posInNm.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double gy = -posInNm.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double gz = -posInNm.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
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
                a.setXYZGradient(gx, gy, gz);
                g[index++] = gx;
                g[index++] = gy;
                g[index++] = gz;
            }
        }
        return g;
    }

    /**
     * setOpenMMPositions takes in an array of doubles generated by the DYN
     * reader method and appends these values to a Vec3Array. Finally this
     * method sets the created Vec3Array as the positions of the context.
     *
     * @param x                 an array of {@link double} objects.
     * @param numberOfVariables a int.
     */
    public void setOpenMMPositions(double x[], int numberOfVariables) {
        assert numberOfVariables == getNumberOfVariables();
        if (positions == null) {
            positions = OpenMM_Vec3Array_create(0);
        } else {
            OpenMM_Vec3Array_resize(positions, 0);
        }
        OpenMM_Vec3.ByValue pos = new OpenMM_Vec3.ByValue();
        for (int i = 0; i < numberOfVariables; i = i + 3) {
            pos.x = x[i] * OpenMM_NmPerAngstrom;
            pos.y = x[i + 1] * OpenMM_NmPerAngstrom;
            pos.z = x[i + 2] * OpenMM_NmPerAngstrom;
            OpenMM_Vec3Array_append(positions, pos);
        }
        OpenMM_Context_setPositions(context, positions);
    }

    /**
     * setOpenMMVelocities takes in an array of doubles generated by the DYN
     * reader method and appends these values to a Vec3Array. Finally this
     * method sets the created Vec3Arrat as the velocities of the context.
     *
     * @param v                 an array of {@link double} objects.
     * @param numberOfVariables a int.
     */
    public void setOpenMMVelocities(double v[], int numberOfVariables) {
        assert numberOfVariables == getNumberOfVariables();
        if (velocities == null) {
            velocities = OpenMM_Vec3Array_create(0);
        } else {
            OpenMM_Vec3Array_resize(velocities, 0);
        }
        OpenMM_Vec3.ByValue vel = new OpenMM_Vec3.ByValue();
        for (int i = 0; i < numberOfVariables; i = i + 3) {
            vel.x = v[i] * OpenMM_NmPerAngstrom;
            vel.y = v[i + 1] * OpenMM_NmPerAngstrom;
            vel.z = v[i + 2] * OpenMM_NmPerAngstrom;
            OpenMM_Vec3Array_append(velocities, vel);
        }
        OpenMM_Context_setVelocities(context, velocities);
    }

    /**
     * getOpenMMPositions takes in a PointerByReference containing the position
     * information of the context. This method creates a Vec3Array that contains
     * the three dimensional information of the positions of the atoms. The
     * method then adds these values to a new double array x and returns it to
     * the method call
     *
     * @param positions         a {@link com.sun.jna.ptr.PointerByReference} object.
     * @param numberOfVariables a int.
     * @param x                 an array of {@link double} objects.
     * @return x
     */
    public double[] getOpenMMPositions(PointerByReference positions, int numberOfVariables, double x[]) {
        assert numberOfVariables == getNumberOfVariables();
        if (x == null || x.length < numberOfVariables) {
            x = new double[numberOfVariables];
        }
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            int offset = i * 3;
            OpenMM_Vec3 pos = OpenMM_Vec3Array_get(positions, i);
            x[offset] = pos.x * OpenMM_AngstromsPerNm;
            x[offset + 1] = pos.y * OpenMM_AngstromsPerNm;
            x[offset + 2] = pos.z * OpenMM_AngstromsPerNm;
            Atom atom = atoms[i];
            atom.moveTo(x[offset], x[offset + 1], x[offset + 2]);
        }
        return x;
    }

    /**
     * getOpenMMVelocities takes in a PointerByReference containing the velocity
     * information of the context. This method creates a Vec3Array that contains
     * the three dimensional information of the velocities of the atoms. This
     * method then adds these values to a new double array v and returns it to
     * the method call
     *
     * @param velocities        a {@link com.sun.jna.ptr.PointerByReference} object.
     * @param numberOfVariables a int.
     * @param v                 an array of {@link double} objects.
     * @return an array of {@link double} objects.
     */
    public double[] getOpenMMVelocities(PointerByReference velocities, int numberOfVariables, double v[]) {
        assert numberOfVariables == getNumberOfVariables();
        if (v == null || v.length < numberOfVariables) {
            v = new double[numberOfVariables];
        }
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            int offset = i * 3;
            OpenMM_Vec3 vel = OpenMM_Vec3Array_get(velocities, i);
            v[offset] = vel.x * OpenMM_AngstromsPerNm;
            v[offset + 1] = vel.y * OpenMM_AngstromsPerNm;
            v[offset + 2] = vel.z * OpenMM_AngstromsPerNm;
            Atom atom = atoms[i];
            double velocity[] = {v[offset], v[offset + 1], v[offset + 2]};
            atom.setVelocity(velocity);
        }
        return v;
    }

    /**
     * getOpenMMAccelerations takes in a PointerByReference containing the force
     * information of the context. This method creates a Vec3Array that contains
     * the three dimensional information of the forces on the atoms. This method
     * then adds these values (divided by mass, effectively turning them into
     * accelerations) to a new double array a and returns it to the method call
     *
     * @param forces            a {@link com.sun.jna.ptr.PointerByReference} object.
     * @param numberOfVariables a int.
     * @param mass              an array of {@link double} objects.
     * @param a                 an array of {@link double} objects.
     * @return an array of {@link double} objects.
     */
    public double[] getOpenMMAccelerations(PointerByReference forces, int numberOfVariables, double[] mass, double[] a) {
        assert numberOfVariables == getNumberOfVariables();
        if (a == null || a.length < numberOfVariables) {
            a = new double[numberOfVariables];
        }
        Atom[] atoms = molecularAssembly.getAtomArray();
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            int offset = i * 3;
            OpenMM_Vec3 acc = OpenMM_Vec3Array_get(forces, i);
            a[offset] = (acc.x * 10.0) / mass[i];
            a[offset + 1] = (acc.y * 10.0) / mass[i + 1];
            a[offset + 2] = (acc.z * 10.0) / mass[i + 2];
            Atom atom = atoms[i];
            double acceleration[] = {a[offset], a[offset + 1], a[offset + 2]};
            atom.setAcceleration(acceleration);
        }
        return a;
    }

    private void setDefaultPeriodicBoxVectors() {

        OpenMM_Vec3 a = new OpenMM_Vec3();
        OpenMM_Vec3 b = new OpenMM_Vec3();
        OpenMM_Vec3 c = new OpenMM_Vec3();

        Crystal crystal = super.getCrystal();

        if (!crystal.aperiodic()) {
            a.x = crystal.a * OpenMM_NmPerAngstrom;
            a.y = 0.0 * OpenMM_NmPerAngstrom;
            a.z = 0.0 * OpenMM_NmPerAngstrom;
            b.x = 0.0 * OpenMM_NmPerAngstrom;
            b.y = crystal.b * OpenMM_NmPerAngstrom;
            b.z = 0.0 * OpenMM_NmPerAngstrom;
            c.x = 0.0 * OpenMM_NmPerAngstrom;
            c.y = 0.0 * OpenMM_NmPerAngstrom;
            c.z = crystal.c * OpenMM_NmPerAngstrom;
            OpenMM_System_setDefaultPeriodicBoxVectors(system, a, b, c);
        }
    }

    private void getPeriodicBoxVectors(PointerByReference state) {

        OpenMM_Vec3 a = new OpenMM_Vec3();
        OpenMM_Vec3 b = new OpenMM_Vec3();
        OpenMM_Vec3 c = new OpenMM_Vec3();

        OpenMM_State_getPeriodicBoxVectors(state, a, b, c);
    }

    /**
     * getIntegrator returns the integrator used for the context
     *
     * @return integrator
     */
    public PointerByReference getIntegrator() {
        return integrator;
    }

    /**
     * Returns the context created for the ForceFieldEnergyOpenMM object
     *
     * @return context
     */
    public PointerByReference getContext() {
        return context;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Platform getPlatform() {
        return ffxPlatform;
    }

    /**
     * Sets the finite-difference step size used for getdEdL.
     *
     * @param fdDLambda FD step size.
     */
    public void setFdDLambda(double fdDLambda) {
        this.fdDLambda = fdDLambda;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getdEdL() {
        // No lambda dependence.
        if (!lambdaTerm) {
            return 0.0;
        }

        // Only vdW has a lambda dependence.
        // if (vdwLambdaTerm && !elecLambdaTerm) {
        //     return vdwdUdL;
        // }

        if (testdEdL) {
            testLambda();
            testdEdL = false;
        }

        double currentLambda = lambda;
        double width;
        double ePlus;
        double eMinus;
        double dEdL;
        double openMMdEdL;

        // Small optimization to only create the x array once.
        double[] x = new double[getNumberOfVariables()];
        getCoordinates(x);

//        if (doOpenMMdEdL) {
        // logger.info(String.format(" Calculating lambda with OpenMM"));
        width = fdDLambda;
        if (currentLambda + fdDLambda > 1.0) {
            logger.fine(" Could not test the upper point, as current lambda + fdDL > 1");
            ePlus = energy(x);
            setLambda(currentLambda - fdDLambda);
            eMinus = energy(x);
        } else if (currentLambda - fdDLambda < 0.0) {
            logger.fine(" Could not test the lower point, as current lambda - fdDL < 1");
            eMinus = energy(x);
            setLambda(currentLambda + fdDLambda);
            ePlus = energy(x);
        } else {
            setLambda(currentLambda + fdDLambda);
            ePlus = energy(x);
            setLambda(currentLambda - fdDLambda);
            eMinus = energy(x);
            //logger.info(String.format(" OpenMM %12.8f %12.8f", ePlus, eMinus));
            width *= 2.0;
        }
        // Reset Lambda.
        setLambda(currentLambda);
        openMMdEdL = (ePlus - eMinus) / width;
        //logger.info(String.format(" Step: %16.10f OpenMM: %16.10f", fdDLambda, openMMdEdL));
        dEdL = openMMdEdL;
//        }

//        double ffxdEdL = 0.0;
//        if (doFFXdEdL || !doOpenMMdEdL) {
//            // logger.info(String.format(" Calculating lambda with FFX"));
//            width = fdDLambda;
//            // This section technically not robust to the case that fdDLambda > 0.5.
//            // However, that should be an error case checked when fdDLambda is set.
//            super.setLambda(1.0);
//            if (currentLambda + fdDLambda > 1.0) {
//                logger.fine(" Could not test the upper point, as current lambda + fdDL > 1");
//                super.setLambdaMultipoleScale(currentLambda);
//                ePlus = super.energy(x, false);
//                // ePlus = super.getTotalElectrostaticEnergy();
//                super.setLambdaMultipoleScale(currentLambda - fdDLambda);
//                eMinus = super.energy(x, false);
//                // eMinus = super.getTotalElectrostaticEnergy();
//            } else if (currentLambda - fdDLambda < 0.0) {
//                logger.fine(" Could not test the lower point, as current lambda - fdDL < 1");
//                super.setLambdaMultipoleScale(currentLambda);
//                eMinus = super.energy(x, false);
//                // eMinus = super.getTotalElectrostaticEnergy();
//                super.setLambdaMultipoleScale(currentLambda + fdDLambda);
//                ePlus = super.energy(x, false);
//                // ePlus = super.getTotalElectrostaticEnergy();
//            } else {
//                super.setLambdaMultipoleScale(currentLambda + fdDLambda);
//                ePlus = super.energy(x, false);
//                // ePlus = super.getTotalElectrostaticEnergy();
//                super.setLambdaMultipoleScale(currentLambda - fdDLambda);
//                eMinus = super.energy(x, false);
//                // eMinus = super.getTotalElectrostaticEnergy();
//                //logger.info(String.format(" FFX    %12.8f %12.8f", ePlus, eMinus));
//                width *= 2.0;
//            }
//            super.setLambdaMultipoleScale(currentLambda);
//            ffxdEdL = (ePlus - eMinus) / width;
//            //logger.info(String.format(" Step: %16.10f FFX:    %16.10f", fdDLambda, ffxdEdL));
//            dEdL = ffxdEdL;
//        }

        return dEdL;
    }

    /**
     * Test the OpenMM and FFX energy for Lambda = 0 and Lambda = 1.
     */
    public void testLambda() {
        // Save the current value of Lambda.
        double currentLambda = lambda;

        // Small optimization to only create the x array once.
        double[] x = new double[getNumberOfVariables()];
        getCoordinates(x);

        // Test OpenMM at L=0.0 and L=1.
        setLambda(0.0);
        double openMMEnergyZero = energy(x);
        setLambda(1.0);
        double openMMEnergyOne = energy(x);

        // Test FFX at L=0 and L=1.
        super.setLambda(1.0);
        super.setLambdaMultipoleScale(0.0);
        double ffxEnergyZero = super.energy(x, false);
        super.setLambdaMultipoleScale(1.0);
        double ffxEnergyOne = super.energy(x, false);
        super.setLambdaMultipoleScale(currentLambda);

        setLambda(currentLambda);

        logger.info(format(" OpenMM Energy at L=0.0: %16.8f and L=1: %16.8f", openMMEnergyZero, openMMEnergyOne));
        logger.info(format(" FFX Energy    at L=0.0: %16.8f and L=1: %16.8f", ffxEnergyZero, ffxEnergyOne));
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double gradients[]) {
        // Note for ForceFieldEnergyOpenMM this method is not implemented.
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public double getd2EdL2() {
        return 0.0;
    }

    /**
     * <p>getVDWRadius.</p>
     *
     * @return a double.
     */
    public double getVDWRadius() {
        VanDerWaals vdW = super.getVdwNode();
        if (vdW == null) {
            return -1.0;
        }

        /**
         * Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
         * for epsilon is supported.
         */
        VanDerWaalsForm vdwForm = vdW.getVDWForm();
        if (vdwForm.vdwType != LENNARD_JONES
                || vdwForm.radiusRule != ARITHMETIC
                || vdwForm.epsilonRule != GEOMETRIC) {
            logger.log(Level.SEVERE, String.format(" Unsupported van der Waals functional form."));
            return -1.0;
        }

        /**
         * OpenMM vdW force requires a diameter (i.e. not radius).
         */
        double radScale = 1.0;
        if (vdwForm.radiusSize == RADIUS) {
            radScale = 2.0;
        }
        /**
         * OpenMM vdw force requires atomic sigma values (i.e. not r-min).
         */
        if (vdwForm.radiusType == R_MIN) {
            radScale /= 1.122462048309372981;
        }

        return radScale;
    }

    /**
     * <p>setUpHydrogenConstraints.</p>
     *
     * @param system a {@link com.sun.jna.ptr.PointerByReference} object.
     */
    public void setUpHydrogenConstraints(PointerByReference system) {
        int i;
        int iAtom1;
        int iAtom2;

        //Atom[] atoms = molecularAssembly.getAtomArray();
        Bond[] bonds = super.getBonds();

        logger.info(String.format(" Setting up Hydrogen constraints"));

        if (bonds == null || bonds.length < 1) {
            return;
        }
        int nBonds = bonds.length;
        Atom atom1;
        Atom atom2;
        Atom parentAtom;
        Bond bondForBondLength;
        BondType bondType;

        for (i = 0; i < nBonds; i++) {
            Bond bond = bonds[i];
            atom1 = bond.getAtom(0);
            atom2 = bond.getAtom(1);
            if (atom1.isHydrogen()) {
                parentAtom = atom1.getBonds().get(0).get1_2(atom1);
                bondForBondLength = atom1.getBonds().get(0);
                bondType = bondForBondLength.bondType;
                iAtom1 = atom1.getXyzIndex() - 1;
                iAtom2 = parentAtom.getXyzIndex() - 1;
                OpenMM_System_addConstraint(system, iAtom1, iAtom2, bondForBondLength.bondType.distance * OpenMM_NmPerAngstrom);
            } else if (atom2.isHydrogen()) {
                parentAtom = atom2.getBonds().get(0).get1_2(atom2);
                bondForBondLength = atom2.getBonds().get(0);
                bondType = bondForBondLength.bondType;
                iAtom1 = atom2.getXyzIndex() - 1;
                iAtom2 = parentAtom.getXyzIndex() - 1;
                OpenMM_System_addConstraint(system, iAtom1, iAtom2, bondForBondLength.bondType.distance * OpenMM_NmPerAngstrom);
            }
        }
    }

    /**
     * <p>calculateDegreesOfFreedom.</p>
     *
     * @return a int.
     */
    public int calculateDegreesOfFreedom() {
        int dof = numParticles * 3;
        dof = dof - OpenMM_System_getNumConstraints(system);
        if (commRemover != null) {
            dof -= 3;
        }
        return dof;
    }

    /**
     * <p>Getter for the field <code>numParticles</code>.</p>
     *
     * @return a int.
     */
    public int getNumParticles() {
        return numParticles;
    }

    /**
     * <p>addRESPA.</p>
     *
     * @param inner a double.
     * @param dt    a double.
     * @return a {@link com.sun.jna.ptr.PointerByReference} object.
     */
    public PointerByReference addRESPA(double inner, double dt) {

        integrator = OpenMM_CustomIntegrator_create(dt);

        // Update the Context with the new integrator (i.e. without creating a new context).
        // OpenMM_CustomIntegrator_addUpdateContextState (integrator);
        int n = (int) (Math.round(dt / inner));
        StringBuffer e1 = new StringBuffer("v+0.5*(dt/" + n + ")*f0/m");
        StringBuffer e11 = new StringBuffer(n + "*(x-x1)/dt+" + e1);
        StringBuffer e2 = new StringBuffer("x+(dt/" + n + ")*v");

        OpenMM_CustomIntegrator_addPerDofVariable(integrator, "x1", 0.0);
        OpenMM_CustomIntegrator_addComputePerDof(integrator, "v", "v+0.5*dt*f1/m");
        for (int i = 0; i < n; i++) {
            OpenMM_CustomIntegrator_addComputePerDof(integrator, "v", e1.toString());
            OpenMM_CustomIntegrator_addComputePerDof(integrator, "x", e2.toString());
            OpenMM_CustomIntegrator_addComputePerDof(integrator, "x1", "x");
            OpenMM_CustomIntegrator_addConstrainPositions(integrator);
            OpenMM_CustomIntegrator_addComputePerDof(integrator, "v", e11.toString());
            OpenMM_CustomIntegrator_addConstrainVelocities(integrator);
        }
        OpenMM_CustomIntegrator_addComputePerDof(integrator, "v", "v+0.5*dt*f1/m");
        OpenMM_CustomIntegrator_addConstrainVelocities(integrator);

        return integrator;
    }

}
