//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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
import static java.lang.Double.isFinite;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.round;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toRadians;

import edu.rit.mp.CharacterBuf;
import edu.rit.pj.Comm;
import edu.uiowa.jopenmm.AmoebaOpenMMLibrary;
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
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setAlchemicalMethod;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setCutoffDistance;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setParticleExclusions;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setParticleParameters;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setSoftcoreAlpha;
import static edu.uiowa.jopenmm.AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_setSoftcorePower;
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
import ffx.utilities.Constants;
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
     * OpenMM Context.
     */
    private OpenMMContext openMMContext;
    /**
     * OpenMM System.
     */
    private OpenMMSystem openMMSystem;
    /**
     * The atoms this ForceFieldEnergyOpenMM operates on.
     */
    private final Atom[] atoms;
    /**
     * Number of particles.
     */
    private int numParticles = 0;
    /**
     * Whether to enforce periodic boundary conditions when obtaining new States.
     */
    public final int enforcePBC;
    /**
     * Truncate the normal OpenMM Lambda Path from 0..1 to Lambda_Start..1. This is useful for conformational
     * optimization if full removal of vdW interactions is not desired (i.e. lambdaStart = ~0.2).
     */
    private double lambdaStart = 0.0;
    /**
     * Lambda step size for finite difference dU/dL.
     */
    private double finiteDifferenceStepSize;
    /**
     * Use two-sided finite difference dU/dL.
     */
    private boolean twoSidedFiniteDifference = true;

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

        logger.info("\n Initializing OpenMM");

        ForceField forceField = molecularAssembly.getForceField();
        openMMContext = new OpenMMContext(forceField, requestedPlatform);

        atoms = molecularAssembly.getAtomArray();
        openMMSystem = new OpenMMSystem(molecularAssembly);

        // Set periodic box vectors.
        setDefaultPeriodicBoxVectors();

        boolean aperiodic = super.getCrystal().aperiodic();
        boolean pbcEnforced = forceField.getBoolean("ENFORCE_PBC", !aperiodic);
        enforcePBC = pbcEnforced ? OpenMM_True : OpenMM_False;

        finiteDifferenceStepSize = forceField.getDouble("FD_DLAMBDA", 0.001);
        twoSidedFiniteDifference = forceField.getBoolean("FD_TWO_SIDED", twoSidedFiniteDifference);
    }

    /**
     * Create an OpenMM Context.
     *
     * @param integratorString Integrator to use.
     * @param timeStep         Time step.
     * @param temperature      Temperature (K).
     */
    public void createContext(String integratorString, double timeStep, double temperature) {
        openMMContext.createContext(integratorString, timeStep, temperature);
    }

    /**
     * <p>
     * Getter for the field <code>integratorString</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getIntegratorString() {
        return openMMContext.integratorString;
    }

    /**
     * <p>
     * Getter for the field <code>temperature</code>.</p>
     *
     * @return a double.
     */
    public double getTemperature() {
        return openMMContext.temperature;
    }

    /**
     * <p>
     * Getter for the field <code>timeStep</code>.</p>
     *
     * @return a double.
     */
    public double getTimeStep() {
        return openMMContext.timeStep;
    }

    /**
     * <p>
     * setCoeffOfFriction.</p>
     *
     * @param coeffOfFriction a double.
     */
    public void setCoeffOfFriction(double coeffOfFriction) {
        openMMContext.openMMIntegrator.frictionCoeff = coeffOfFriction;
    }

    /**
     * <p>
     * Setter for the field <code>collisionFreq</code>.</p>
     *
     * @param collisionFreq a double.
     */
    public void setCollisionFreq(double collisionFreq) {
        openMMContext.openMMIntegrator.collisionFreq = collisionFreq;
    }

    /**
     * Adds a center-of-mass motion remover to the Potential. Not advised for
     * anything not running MD using the OpenMM library (i.e. OpenMMMolecularDynamics).
     * Has caused bugs with the FFX MD class.
     *
     * @param addIfDuplicate Add a CCOM remover even if it already exists
     */
    public void addCOMMRemover(boolean addIfDuplicate) {
        openMMSystem.addCOMMRemover(addIfDuplicate);
    }

    /**
     * Add an Andersen thermostat to the system.
     *
     * @param targetTemp Target temperature in Kelvins.
     */
    public void addAndersenThermostat(double targetTemp) {
        /*
         * Citation:
         * Andersen, H. C., Molecular dynamics simulations at constant pressure and/or temperature. The Journal of Chemical Physics 1980, 72 (4), 2384-2393.
         */
        addAndersenThermostat(targetTemp, openMMContext.openMMIntegrator.collisionFreq);
    }

    /**
     * Add a Monte Carlo Barostat to the system.
     *
     * @param targetPressure The target pressure.
     * @param targetTemp     The target temperature.
     * @param frequency      The frequency to apply the barostat.
     */
    public void addMonteCarloBarostat(double targetPressure, double targetTemp, int frequency) {
        openMMSystem.addMonteCarloBarostat(targetPressure, targetTemp, frequency);
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
     * Evaluates energy and gradients with both OpenMM and reference potential,
     * and returns the difference in energies.
     *
     * @param x       Coordinate array.
     * @param gOMM    To be filled with OpenMM-calculated gradients.
     * @param gFFX    To be filled with FFX-calculated gradients.
     * @param verbose Whether to be verbose.
     * @return Difference in energies.
     */
    public double energyAndGradientVsFFX(double[] x, double[] gOMM, double[] gFFX, boolean verbose) {
        double ffxE = super.energyAndGradient(x, gFFX, verbose);
        double thisE = energyAndGradient(x, gOMM, verbose);
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
     * setOpenMMPositions takes in an array of doubles generated by the DYN
     * reader method and appends these values to a Vec3Array. Finally this
     * method sets the created Vec3Array as the positions of the context.
     *
     * @param x                 an array of {@link double} objects.
     * @param numberOfVariables a int.
     */
    public void setOpenMMPositions(double[] x, int numberOfVariables) {
        int nOMMVar = atoms.length * 3;
        int xOffset = 0;
        if (numberOfVariables != nOMMVar) {
            double[] newX = new double[nOMMVar];
            for (int i = 0; i < atoms.length; i++) {
                int i3 = 3 * i;
                if (atoms[i].isActive()) {
                    System.arraycopy(x, xOffset, newX, i3, 3);
                    xOffset += 3;
                } else {
                    double[] currXYZ = new double[3];
                    currXYZ = atoms[i].getXYZ(currXYZ);
                    System.arraycopy(currXYZ, 0, newX, i3, 3);
                }
            }
            x = newX;
        }
        openMMContext.setOpenMMCoordinates(x, nOMMVar);
    }

    /**
     * setOpenMMVelocities takes in an array of doubles generated by the DYN
     * reader method and appends these values to a Vec3Array. Finally this
     * method sets the created Vec3Arrat as the velocities of the context.
     *
     * @param v                 an array of {@link double} objects.
     * @param numberOfVariables a int.
     */
    public void setOpenMMVelocities(double[] v, int numberOfVariables) {
        openMMContext.setOpenMMVelocities(v, numberOfVariables);
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
    public double[] getOpenMMPositions(PointerByReference positions, int numberOfVariables, double[] x) {
        return openMMContext.getOpenMMPositions(positions, numberOfVariables, x);
    }

    public void getOpenMMVelocities(PointerByReference velocities, int numberOfVariables, double[] v) {
        openMMContext.getOpenMMVelocities(velocities, numberOfVariables, v);
    }

    public void getOpenMMAccelerations(PointerByReference accelerations, int numberOfVariables, double[] m, double[] a) {
        openMMContext.getOpenMMAccelerations(accelerations, numberOfVariables, m, a);
    }

    /**
     * Update active atoms.
     */
    public void setActiveAtoms() {
        openMMSystem.updateAtomMass();
        openMMContext.reinitContext();
    }

    /**
     * getIntegrator returns the integrator used for the context
     *
     * @return integrator
     */
    public PointerByReference getIntegrator() {
        return openMMContext.getIntegrator();
    }

    /**
     * Returns the context created for the ForceFieldEnergyOpenMM object
     *
     * @return context
     */
    public PointerByReference getContext() {
        return openMMContext.getContext();
    }

    /**
     * Describes the Context this is running with.
     *
     * @return String describing the context.
     */
    public String describeContext() {
        return openMMContext.toString();
    }

    /**
     * Sets the finite-difference step size used for getdEdL.
     *
     * @param finiteDifferenceStepSize FD step size.
     */
    public void setFiniteDifferenceStepSize(double finiteDifferenceStepSize) {
        this.finiteDifferenceStepSize = finiteDifferenceStepSize;
    }

    /**
     * <p>
     * calculateDegreesOfFreedom.</p>
     *
     * @return a int.
     */
    public int calculateDegreesOfFreedom() {
        return openMMSystem.calculateDegreesOfFreedom();
    }

    /**
     * <p>
     * Getter for the field <code>numParticles</code>.</p>
     *
     * @return a int.
     */
    public int getNumParticles() {
        return numParticles;
    }

    /**
     * Compute the energy using the pure Java code path.
     *
     * @param x Input atomic coordinates
     * @return The energy (kcal/mol)
     */
    public double energyFFX(double[] x) {
        return energyFFX(x, false);
    }

    /**
     * Compute the energy using the pure Java code path.
     *
     * @param x       Input atomic coordinates
     * @param verbose Use verbose logging.
     * @return The energy (kcal/mol)
     */
    public double energyFFX(double[] x, boolean verbose) {
        return super.energy(x, verbose);
    }

    /**
     * Compute the energy and gradient using the pure Java code path.
     *
     * @param x Input atomic coordinates
     * @param g Storage for the gradient vector.
     * @return The energy (kcal/mol)
     */
    public double energyAndGradientFFX(double[] x, double[] g) {
        return energyAndGradientFFX(x, g, false);
    }

    /**
     * Compute the energy and gradient using the pure Java code path.
     *
     * @param x       Input atomic coordinates
     * @param g       Storage for the gradient vector.
     * @param verbose Use verbose logging.
     * @return The energy (kcal/mol)
     */
    public double energyAndGradientFFX(double[] x, double[] g, boolean verbose) {
        return super.energyAndGradient(x, g, verbose);
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

        // Make sure a context has been created.
        getContext();

        updateParameters(atoms);

        // Unscale the coordinates.
        unscaleCoordinates(x);

        setCoordinates(x);
        setOpenMMPositions(x, x.length);

        PointerByReference state = openMMContext.getState(OpenMM_State_Energy);

        double e = OpenMM_State_getPotentialEnergy(state) / OpenMM_KJPerKcal;
        if (!isFinite(e)) {
            String message = String.format(" Energy from OpenMM was a non-finite %8g", e);
            logger.warning(message);
            if (lambdaTerm) {
                openMMSystem.printLambdaValues();
            }
            throw new EnergyException(message);
        }
        OpenMM_State_destroy(state);

        if (verbose) {
            logger.log(Level.INFO, String.format(" OpenMM Energy: %14.10g", e));
        }

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
        if (lambdaBondedTerms) {
            return 0.0;
        }

        // Un-scale the coordinates.
        unscaleCoordinates(x);

        // Make sure a context has been created.
        getContext();

        setCoordinates(x);
        setOpenMMPositions(x, x.length);

        int infoMask = OpenMM_State_Energy + OpenMM_State_Forces;
        PointerByReference state = openMMContext.getState(infoMask);

        double e = OpenMM_State_getPotentialEnergy(state) / OpenMM_KJPerKcal;
        if (!isFinite(e)) {
            String message = String.format(" Energy from OpenMM was a non-finite %8g", e);
            logger.warning(message);
            if (lambdaTerm) {
                openMMSystem.printLambdaValues();
            }
            throw new EnergyException(message);
        }

        // if (vdwLambdaTerm) {
        //    PointerByReference parameterArray = OpenMM_State_getEnergyParameterDerivatives(state);
        //    int numDerives = OpenMM_ParameterArray_getSize(parameterArray);
        //    if (numDerives > 0) {
        //        double vdwdUdL = OpenMM_ParameterArray_get(parameterArray, pointerForString("vdw_lambda")) / OpenMM_KJPerKcal;
        //    }
        // }

        if (maxDebugGradient < Double.MAX_VALUE) {
            boolean extremeGrad = Arrays.stream(g).anyMatch((double gi)
                    -> (gi > maxDebugGradient || gi < -maxDebugGradient));
            if (extremeGrad) {
                File origFile = molecularAssembly.getFile();
                String timeString = LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyy_MM_dd-HH_mm_ss"));

                String filename = format("%s-LARGEGRAD-%s.pdb",
                        FilenameUtils.removeExtension(molecularAssembly.getFile().getName()), timeString);
                PotentialsFunctions ef = new PotentialsUtils();
                filename = ef.versionFile(filename);

                logger.warning(format(" Excessively large gradients detected; printing snapshot to file %s", filename));
                ef.saveAsPDB(molecularAssembly, new File(filename));
                molecularAssembly.setFile(origFile);
            }
        }

        if (verbose) {
            logger.log(Level.INFO, format(" OpenMM Energy: %14.10g", e));
        }

        PointerByReference forces = OpenMM_State_getForces(state);
        fillGradients(forces, g);

        // Scale the coordinates and gradients.
        scaleCoordinatesAndGradient(x, g);

        OpenMM_State_destroy(state);
        return e;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void setCrystal(Crystal crystal) {
        super.setCrystal(crystal);
        openMMContext.setPeriodicBoxVectors();
    }

    /**
     * {@inheritDoc}
     *
     * <p>
     * getGradient</p>
     */
    @Override
    public double[] getGradient(double[] g) {
        return openMMContext.getGradients(g);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Platform getPlatform() {
        return openMMContext.ffxPlatform;
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

        // Small optimization to only create the x array once.
        double[] x = new double[getNumberOfVariables()];
        getCoordinates(x);

        double currentLambda = lambda;
        double width = finiteDifferenceStepSize;
        double ePlus;
        double eMinus;

        if (twoSidedFiniteDifference) {
            if (currentLambda + finiteDifferenceStepSize > 1.0) {
                setLambda(currentLambda - finiteDifferenceStepSize);
                eMinus = energy(x);
                setLambda(currentLambda);
                ePlus = energy(x);
            } else if (currentLambda - finiteDifferenceStepSize < 0.0) {
                setLambda(currentLambda + finiteDifferenceStepSize);
                ePlus = energy(x);
                setLambda(currentLambda);
                eMinus = energy(x);
            } else {
                // Two sided finite difference estimate of dE/dL.
                setLambda(currentLambda + finiteDifferenceStepSize);
                ePlus = energy(x);
                setLambda(currentLambda - finiteDifferenceStepSize);
                eMinus = energy(x);
                width *= 2.0;
                setLambda(currentLambda);
            }
        } else {
            // One sided finite difference estimates of dE/dL
            if (currentLambda + finiteDifferenceStepSize > 1.0) {
                setLambda(currentLambda - finiteDifferenceStepSize);
                eMinus = energy(x);
                setLambda(currentLambda);
                ePlus = energy(x);
            } else {
                setLambda(currentLambda + finiteDifferenceStepSize);
                ePlus = energy(x);
                setLambda(currentLambda);
                eMinus = energy(x);
            }
        }

        // Compute the finite difference derivative.
        double dEdL = (ePlus - eMinus) / width;

        // logger.info(format(" getdEdL currentLambda: CL=%8.6f L=%8.6f dEdL=%12.6f", currentLambda, lambda, dEdL));
        return dEdL;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void getdEdXdL(double[] gradients) {
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
     * {@inheritDoc}
     */
    @Override
    public void setLambda(double lambda) {
        if (!lambdaTerm) {
            logger.fine(" Attempting to set a lambda value on a ForceFieldEnergyOpenMM with lambdaterm false.");
            return;
        }

        // Check for lambda outside the range [0 .. 1].
        if (lambda < 0.0 || lambda > 1.0) {
            String message = format(" Lambda value %8.3f is not in the range [0..1].", lambda);
            logger.warning(message);
            return;
        }

        super.setLambda(lambda);

        this.lambda = lambda;

        // Remove the beginning of the normal Lambda path.
        double mappedLambda = lambda;
        if (lambdaStart > 0) {
            double windowSize = 1.0 - lambdaStart;
            mappedLambda = lambdaStart + lambda * windowSize;
        }

        if (openMMSystem != null) {
            openMMSystem.setLambda(mappedLambda);
            if (atoms != null) {
                List<Atom> atomList = new ArrayList<>();
                for (Atom atom : atoms) {
                    if (atom.applyLambda()) {
                        atomList.add(atom);
                    }
                }
                // Update force field parameters based on defined lambda values.
                updateParameters(atomList.toArray(new Atom[0]));
            } else {
                updateParameters(null);
            }
        }
    }

    /**
     * Update parameters if the Use flags and/or Lambda value has changed.
     *
     * @param atoms Atoms in this list are considered.
     */
    public void updateParameters(Atom[] atoms) {
        openMMSystem.updateParameters(atoms);
    }

    public void setOpenMMPeriodicBoxVectors() {
        openMMContext.setPeriodicBoxVectors();
    }

    public void getPeriodicBoxVectors(PointerByReference state) {
        openMMSystem.getPeriodicBoxVectors(state);
    }

    /**
     * Creates and manage an OpenMM Context.
     */
    private class OpenMMContext {

        /**
         * OpenMM Integrator.
         */
        private OpenMMIntegrator openMMIntegrator;

        /**
         * OpenMM Context.
         */
        private PointerByReference context = null;
        /**
         * Requested Platform (i.e. Java or an OpenMM platform).
         */
        private final Platform ffxPlatform;
        /**
         * OpenMM Platform.
         */
        private PointerByReference platform = null;
        /**
         * Integrator string (default = VERLET).
         */
        private String integratorString = "VERLET";
        /**
         * Time step (default = 0.001 psec).
         */
        private double timeStep = 0.001;
        /**
         * Temperature (default = 298.15).
         */
        private double temperature = 298.15;
        /**
         * Constraint tolerance as a fraction of the constrained bond length.
         */
        private double constraintTolerance = ForceFieldEnergy.DEFAULT_CONSTRAINT_TOLERANCE;

        OpenMMContext(ForceField forceField, Platform requestedPlatform) {
            loadPlatform(requestedPlatform);
            ffxPlatform = requestedPlatform;
            openMMIntegrator = new OpenMMIntegrator(forceField, constraintTolerance);
        }

        public PointerByReference getContext() {
            if (context == null) {
                createContext(integratorString, timeStep, temperature);
            }
            return context;
        }

        public PointerByReference getIntegrator() {
            return openMMIntegrator.getIntegrator();
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
        double[] getOpenMMVelocities(PointerByReference velocities, int numberOfVariables, double[] v) {
            assert numberOfVariables == getNumberOfVariables();
            if (v == null || v.length < numberOfVariables) {
                v = new double[numberOfVariables];
            }
            int nAtoms = atoms.length;
            for (int i = 0; i < nAtoms; i++) {
                int offset = i * 3;
                OpenMM_Vec3 vel = OpenMM_Vec3Array_get(velocities, i);
                v[offset] = vel.x * OpenMM_AngstromsPerNm;
                v[offset + 1] = vel.y * OpenMM_AngstromsPerNm;
                v[offset + 2] = vel.z * OpenMM_AngstromsPerNm;
                Atom atom = atoms[i];
                double[] velocity = {v[offset], v[offset + 1], v[offset + 2]};
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
        double[] getOpenMMAccelerations(PointerByReference forces, int numberOfVariables, double[] mass, double[] a) {
            assert numberOfVariables == getNumberOfVariables();
            if (a == null || a.length < numberOfVariables) {
                a = new double[numberOfVariables];
            }
            int nAtoms = atoms.length;
            for (int i = 0; i < nAtoms; i++) {
                int offset = i * 3;
                OpenMM_Vec3 acc = OpenMM_Vec3Array_get(forces, i);
                a[offset] = (acc.x * OpenMM_AngstromsPerNm) / mass[i];
                a[offset + 1] = (acc.y * OpenMM_AngstromsPerNm) / mass[i + 1];
                a[offset + 2] = (acc.z * OpenMM_AngstromsPerNm) / mass[i + 2];
                Atom atom = atoms[i];
                double[] acceleration = {a[offset], a[offset + 1], a[offset + 2]};
                atom.setAcceleration(acceleration);
            }
            return a;
        }

        public double[] getGradients(double[] g) {
            PointerByReference state = OpenMM_Context_getState(context, OpenMM_State_Forces, enforcePBC);
            PointerByReference forces = OpenMM_State_getForces(state);
            g = fillGradients(forces, g);
            OpenMM_State_destroy(state);
            return g;
        }

        /**
         * createContext takes in a parameters to determine which integrator the
         * user requested during the start up of the simulation. A switch statement
         * is used with Strings as the variable to determine between Langevin,
         * Brownian, Custom, Compound and Verlet integrator
         *
         * @param integratorString a {@link java.lang.String} object.
         * @param timestep         Timestep in psec.
         * @param temperature      a double.
         */
        void createContext(String integratorString, double timestep, double temperature) {
            this.integratorString = integratorString;
            this.timeStep = timestep;
            this.temperature = temperature;

            if (context != null) {
                OpenMM_Context_destroy(context);
                context = null;
            }

            PointerByReference integrator = openMMIntegrator.createIntegrator(integratorString, this.timeStep, temperature);

            // Set lambda to 1.0 when creating a context to avoid OpenMM compiling out any terms.
            double currentLambda = lambda;
            if (lambdaTerm) {
                ForceFieldEnergyOpenMM.this.setLambda(1.0);
            }

            // Create a context.
            context = OpenMM_Context_create_2(openMMSystem.getOpenMMSystem(), integrator, platform);

            // Revert to the current lambda value.
            if (lambdaTerm) {
                ForceFieldEnergyOpenMM.this.setLambda(currentLambda);
            }

            // Set initial positions.
            double[] x = new double[numParticles * 3];
            int index = 0;
            for (int i = 0; i < numParticles; i++) {
                Atom atom = atoms[i];
                x[index] = atom.getX();
                x[index + 1] = atom.getY();
                x[index + 2] = atom.getZ();
                index += 3;
            }

            // Load current atomic positions.
            setOpenMMPositions(x, numParticles * 3);

            // Apply constraints starting from current atomic positions.
            OpenMM_Context_applyConstraints(context, constraintTolerance);

            // Get back constrained atomic coordinates for consistency.
            PointerByReference state = OpenMM_Context_getState(context, OpenMM_State_Positions, enforcePBC);
            PointerByReference positions = OpenMM_State_getPositions(state);
            getOpenMMPositions(positions, numParticles * 3, x);
            OpenMM_State_destroy(state);

            logger.info(format(" Context created (integrator=%s, time step=%6.2f fsec, temperature=%6.2f K).\n",
                    integratorString, this.timeStep * Constants.PSEC_TO_FSEC, temperature));
        }

        void destroyContext() {
            if (context != null) {
                OpenMM_Context_destroy(context);
            }
            context = null;
            openMMIntegrator.destroyIntegrator();
        }

        void reinitContext() {
            if (context != null) {
                int preserveState = 1;
                OpenMM_Context_reinitialize(context, preserveState);
            }
        }

        void setParameter(String name, double value) {
            if (context != null) {
                OpenMM_Context_setParameter(context, name, value);
            }
        }

        PointerByReference getState(int infoMask) {
            return OpenMM_Context_getState(context, infoMask, enforcePBC);
        }

        void setOpenMMCoordinates(double[] x, int numberOfVariables) {
            assert numberOfVariables == getNumberOfVariables();
            PointerByReference positions = OpenMM_Vec3Array_create(0);
            OpenMM_Vec3.ByValue pos = new OpenMM_Vec3.ByValue();
            for (int i = 0; i < numberOfVariables; i = i + 3) {
                pos.x = x[i] * OpenMM_NmPerAngstrom;
                pos.y = x[i + 1] * OpenMM_NmPerAngstrom;
                pos.z = x[i + 2] * OpenMM_NmPerAngstrom;
                OpenMM_Vec3Array_append(positions, pos);
            }
            OpenMM_Context_setPositions(context, positions);
            OpenMM_Vec3Array_destroy(positions);
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
        double[] getOpenMMPositions(PointerByReference positions, int numberOfVariables, double[] x) {
            assert numberOfVariables == getNumberOfVariables();
            if (x == null || x.length < numberOfVariables) {
                x = new double[numberOfVariables];
            }
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

        void setOpenMMVelocities(double[] v, int numberOfVariables) {
            assert numberOfVariables == getNumberOfVariables();
            PointerByReference velocities = OpenMM_Vec3Array_create(0);
            OpenMM_Vec3.ByValue vel = new OpenMM_Vec3.ByValue();
            for (int i = 0; i < numberOfVariables; i = i + 3) {
                vel.x = v[i] * OpenMM_NmPerAngstrom;
                vel.y = v[i + 1] * OpenMM_NmPerAngstrom;
                vel.z = v[i + 2] * OpenMM_NmPerAngstrom;
                OpenMM_Vec3Array_append(velocities, vel);
            }
            OpenMM_Context_setVelocities(context, velocities);
            OpenMM_Vec3Array_destroy(velocities);
        }

        private void setPeriodicBoxVectors() {
            Crystal crystal = getCrystal();
            if (!crystal.aperiodic()) {
                OpenMM_Vec3 a = new OpenMM_Vec3();
                OpenMM_Vec3 b = new OpenMM_Vec3();
                OpenMM_Vec3 c = new OpenMM_Vec3();
                double[][] Ai = crystal.Ai;
                a.x = Ai[0][0] * OpenMM_NmPerAngstrom;
                a.y = Ai[0][1] * OpenMM_NmPerAngstrom;
                a.z = Ai[0][2] * OpenMM_NmPerAngstrom;
                b.x = Ai[1][0] * OpenMM_NmPerAngstrom;
                b.y = Ai[1][1] * OpenMM_NmPerAngstrom;
                b.z = Ai[1][2] * OpenMM_NmPerAngstrom;
                c.x = Ai[2][0] * OpenMM_NmPerAngstrom;
                c.y = Ai[2][1] * OpenMM_NmPerAngstrom;
                c.z = Ai[2][2] * OpenMM_NmPerAngstrom;
                OpenMM_Context_setPeriodicBoxVectors(context, a, b, c);
            }
        }

        /**
         * Load an OpenMM Platform
         */
        private void loadPlatform(Platform requestedPlatform) {

            OpenMMUtils.init();

            // Print out the OpenMM Version.
            Pointer version = OpenMM_Platform_getOpenMMVersion();
            logger.log(Level.INFO, " Version: {0}", version.getString(0));

            // Print out the OpenMM lib directory.
            logger.log(Level.FINE, " Lib Directory:       {0}", OpenMMUtils.getLibDirectory());
            // Load platforms and print out their names.
            PointerByReference libs = OpenMM_Platform_loadPluginsFromDirectory(OpenMMUtils.getLibDirectory());
            int numLibs = OpenMM_StringArray_getSize(libs);
            logger.log(Level.FINE, " Number of libraries: {0}", numLibs);
            for (int i = 0; i < numLibs; i++) {
                String libString = stringFromArray(libs, i);
                logger.log(Level.FINE, "  Library: {0}", libString);
            }
            OpenMM_StringArray_destroy(libs);

            // Print out the OpenMM plugin directory.
            logger.log(Level.INFO, "\n Plugin Directory:  {0}", OpenMMUtils.getPluginDirectory());
            // Load plugins and print out their names.
            PointerByReference plugins = OpenMM_Platform_loadPluginsFromDirectory(OpenMMUtils.getPluginDirectory());
            int numPlugins = OpenMM_StringArray_getSize(plugins);
            logger.log(Level.INFO, " Number of Plugins: {0}", numPlugins);
            boolean cuda = false;
            for (int i = 0; i < numPlugins; i++) {
                String pluginString = stringFromArray(plugins, i);
                logger.log(Level.INFO, "  Plugin: {0}", pluginString);
                if (pluginString != null) {
                    pluginString = pluginString.toUpperCase();
                    boolean amoebaCudaAvailable = pluginString.contains("AMOEBACUDA");
                    if (amoebaCudaAvailable) {
                        cuda = true;
                    }
                }
            }
            OpenMM_StringArray_destroy(plugins);

            int numPlatforms = OpenMM_Platform_getNumPlatforms();
            logger.log(Level.INFO, " Number of Platforms: {0}", numPlatforms);

            if (requestedPlatform == Platform.OMM_CUDA && !cuda) {
                logger.severe(" The OMM_CUDA platform was requested, but is not available.");
            }

            // Extra logging to print out plugins that failed to load.
            if (logger.isLoggable(Level.FINE)) {
                PointerByReference pluginFailers = OpenMM_Platform_getPluginLoadFailures();
                int numFailures = OpenMM_StringArray_getSize(pluginFailers);
                for (int i = 0; i < numFailures; i++) {
                    String pluginString = stringFromArray(pluginFailers, i);
                    logger.log(Level.FINE, " Plugin load failure: {0}", pluginString);
                }
                OpenMM_StringArray_destroy(pluginFailers);
            }

            String defaultPrecision = "mixed";
            String precision = molecularAssembly.getForceField().getString("PRECISION", defaultPrecision).toLowerCase();
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
                int deviceID = molecularAssembly.getForceField().getInteger("CUDA_DEVICE", defaultDevice);
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

        @Override
        public String toString() {
            return String.format(" OpenMM context with integrator %s, timestep %9.3g fsec, temperature %9.3g K, constraintTolerance %9.3g", integratorString, timeStep, temperature, constraintTolerance);
        }
    }

    /**
     * Create and manage an OpenMM Integrator.
     */
    private class OpenMMIntegrator {

        /**
         * OpenMM Integrator.
         */
        private PointerByReference integrator = null;
        /**
         * Constraint tolerance as a fraction of the constrained bond length.
         */
        private double constraintTolerance;
        /**
         * Langevin friction coefficient.
         */
        private double frictionCoeff;
        /**
         * Andersen thermostat collision frequency.
         */
        private double collisionFreq;

        OpenMMIntegrator(ForceField forceField, double constraintTolerance) {
            this.constraintTolerance = constraintTolerance;

            frictionCoeff = forceField.getDouble("FRICTION_COEFF", 91.0);
            collisionFreq = forceField.getDouble("COLLISION_FREQ", 0.01);
        }

        public PointerByReference getIntegrator() {
            return integrator;
        }

        PointerByReference createIntegrator(String integratorString, double timeStep, double temperature) {
            double dt = timeStep;
            switch (integratorString) {
                case "LANGEVIN":
                    createLangevinIntegrator(temperature, frictionCoeff, dt);
                    break;
                case "RESPA":
                    // Read in the inner time step in psec.
                    int in = molecularAssembly.getProperties().getInt("respa-dt", 4);
                    if (in < 2) {
                        in = 2;
                    }
                    double inner = dt / in;
                    createRESPAIntegrator(inner, dt);
                    break;
            /*
            case "BROWNIAN":
                createBrownianIntegrator(temperature, frictionCoeff, dt);
                break;
            case "CUSTOM":
                createCustomIntegrator(dt);
                break;
            case "COMPOUND":
                createCompoundIntegrator();
                break;
             */
                case "VERLET":
                default:
                    createVerletIntegrator(dt);
            }

            return integrator;
        }

        /**
         * <p>
         * createRESPAIntegrator.</p>
         *
         * @param inner a double.
         * @param dt    a double.
         */
        private void createRESPAIntegrator(double inner, double dt) {
            createCustomIntegrator(dt);
            OpenMM_CustomIntegrator_addUpdateContextState(integrator);
            OpenMM_CustomIntegrator_setKineticEnergyExpression(integrator, "m*v*v/2");

            int n = (int) (round(dt / inner));
            StringBuilder e1 = new StringBuilder("v+0.5*(dt/" + n + ")*f0/m");
            StringBuilder e11 = new StringBuilder(n + "*(x-x1)/dt+" + e1);
            StringBuilder e2 = new StringBuilder("x+(dt/" + n + ")*v");

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
        }

        /**
         * Create a Langevin Integrator.
         *
         * @param temperature   Temperature (Kelvin).
         * @param frictionCoeff Frictional coefficient.
         * @param dt            Time step.
         */
        private void createLangevinIntegrator(double temperature, double frictionCoeff, double dt) {
            destroyIntegrator();

            integrator = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, dt);
            CompositeConfiguration properties = molecularAssembly.getProperties();
            if (properties.containsKey("randomseed")) {
                int randomSeed = properties.getInt("randomseed", 0);
                logger.info(String.format(" Setting random seed %d for Langevin dynamics", randomSeed));
                OpenMM_LangevinIntegrator_setRandomNumberSeed(integrator, randomSeed);
            }
            OpenMM_Integrator_setConstraintTolerance(integrator, constraintTolerance);
        }

        private void destroyIntegrator() {
            if (integrator != null) {
                OpenMM_Integrator_destroy(integrator);
                integrator = null;
            }
        }

        /**
         * Create a Verlet Integrator.
         *
         * @param dt Time step.
         */
        private void createVerletIntegrator(double dt) {
            destroyIntegrator();
            integrator = OpenMM_VerletIntegrator_create(dt);
            OpenMM_Integrator_setConstraintTolerance(integrator, constraintTolerance);
        }

        /**
         * Create a Custom Integrator.
         *
         * @param dt Time step.
         */
        private void createCustomIntegrator(double dt) {
            destroyIntegrator();
            integrator = OpenMM_CustomIntegrator_create(dt);
            OpenMM_Integrator_setConstraintTolerance(integrator, constraintTolerance);
        }

    }

    /**
     * Create and manage an OpenMM System.
     */
    private class OpenMMSystem {
        /**
         * The Force Field in use.
         */
        ForceField forceField;
        /**
         * Array of atoms in the sytem.
         */
        Atom[] atoms;
        /**
         * OpenMM System.
         */
        private PointerByReference system;
        /**
         * OpenMM thermostat. Currently an Andersen thermostat is supported.
         */
        private PointerByReference ommThermostat = null;
        /**
         * Barostat to be added if NPT (isothermal-isobaric) dynamics is requested.
         */
        private PointerByReference ommBarostat = null;
        /**
         * OpenMM center-of-mass motion remover.
         */
        private PointerByReference commRemover = null;
        /**
         * OpenMM AMOEBA Torsion Force.
         */
        private PointerByReference amoebaTorsionForce = null;
        /**
         * OpenMM Improper Torsion Force.
         */
        private PointerByReference amoebaImproperTorsionForce = null;
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
         * OpenMM Custom GB Force.
         */
        private PointerByReference customGBForce = null;
        /**
         * OpenMM Fixed Charge Non-Bonded Force.
         */
        private PointerByReference fixedChargeNonBondedForce = null;
        /**
         * Fixed charge softcore vdW force boolean.
         */
        private boolean softcoreCreated = false;
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
         * Lambda flag to indicate control of electrostatic scaling. If both elec
         * and vdW are being scaled, then vdW is scaled first, followed by elec.
         */
        private boolean elecLambdaTerm;
        /**
         * Lambda flag to indicate control of vdW scaling. If both elec and vdW are
         * being scaled, then vdW is scaled first, followed by elec.
         */
        private boolean vdwLambdaTerm;
        /**
         * Lambda flag to indicate control of torsional force constants (L=0
         * corresponds to torsions being off, and L=1 to torsions at full strength.
         */
        private boolean torsionLambdaTerm;
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
         * The lambda value that defines when the electrostatics will start to turn on for full path non-bonded term scaling.
         * <p>
         * A value of 0.6 works well for Chloride ion solvation, which is a difficult case due to
         * the ion having a formal negative charge and a large polarizability.
         */
        private double electrostaticStart = 0.6;
        /**
         * Electrostatics lambda is raised to this power.
         */
        private double electrostaticLambdaPower;
        /**
         * van der Waals softcore alpha.
         */
        private double vdWSoftcoreAlpha = 0.25;
        /**
         * van der Waals softcore beta.
         */
        private double vdwSoftcorePower = 3.0;
        /**
         * Torsional lambda power.
         */
        private double torsionalLambdaPower = 2.0;
        /**
         * When using MELD, our goal will be to scale down the potential by this factor.
         * A negative value indicates we're not using MELD.
         */
        private boolean useMeld;
        private static final double DEFAULT_MELD_SCALE_FACTOR = -1.0;
        private final double meldScaleFactor;

        OpenMMSystem(MolecularAssembly molecularAssembly) {
            // Create the OpenMM System
            system = OpenMM_System_create();
            logger.info(" System created.");

            forceField = molecularAssembly.getForceField();
            atoms = molecularAssembly.getAtomArray();

            // Load atoms.
            try {
                addAtoms();
            } catch (Exception e) {
                logger.severe(" Atom without mass encountered.");
            }

            boolean rigidHydrogen = forceField.getBoolean("RIGID_HYDROGEN", false);

            if (rigidHydrogen) {
                setUpHydrogenConstraints(openMMSystem.getOpenMMSystem());
            }

            boolean rigidBonds = forceField.getBoolean("RIGID_BONDS", false);
            if (rigidBonds) {
                setUpBondConstraints(openMMSystem.getOpenMMSystem());
            }

            boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
            if (rigidHydrogenAngles) {
                setUpHydrogenAngleConstraints(openMMSystem.getOpenMMSystem());
            }

            // Add Bond Force.
            if (rigidBonds) {
                logger.info(" Not creating AmoebaBondForce because bonds are constrained.");
            } else {
                addBondForce();
            }

            // Check for MELD use. If we're using MELD, set all lambda terms to true.
            meldScaleFactor = forceField.getDouble("MELD_SCALE_FACTOR", DEFAULT_MELD_SCALE_FACTOR);
            if (meldScaleFactor <= 1.0 && meldScaleFactor > 0.0) {
                useMeld = true;
                elecLambdaTerm = true;
                vdwLambdaTerm = true;
                torsionLambdaTerm = true;
            } else {
                useMeld = false;
                elecLambdaTerm = false;
                vdwLambdaTerm = false;
                torsionLambdaTerm = false;
            }

            // Read alchemical information -- this needs to be done before creating forces.
            elecLambdaTerm = forceField.getBoolean("ELEC_LAMBDATERM", elecLambdaTerm);
            vdwLambdaTerm = forceField.getBoolean("VDW_LAMBDATERM", vdwLambdaTerm);
            torsionLambdaTerm = forceField.getBoolean("TORSION_LAMBDATERM", torsionLambdaTerm);

            lambdaTerm = (elecLambdaTerm || vdwLambdaTerm || torsionLambdaTerm);

            VanDerWaals vdW = ForceFieldEnergyOpenMM.super.getVdwNode();
            if (vdW != null) {
                vdWSoftcoreAlpha = vdW.getAlpha();
                vdwSoftcorePower = (int) vdW.getBeta();
            }

            // Expand the path [lambda-start .. 1.0] to the interval [0.0 .. 1.0].
            lambdaStart = forceField.getDouble("LAMBDA_START", 0.0);
            if (lambdaStart > 1.0) {
                lambdaStart = 1.0;
            } else if (lambdaStart < 0.0) {
                lambdaStart = 0.0;
            }

            electrostaticStart = forceField.getDouble("ELEC_START", electrostaticStart);
            if (electrostaticStart > 1.0) {
                electrostaticStart = 1.0;
            } else if (electrostaticStart < 0.0) {
                electrostaticStart = 0.0;
            }
            electrostaticLambdaPower = forceField.getDouble("PERMANENT_LAMBDA_EXPONENT", 2.0);

            if (useMeld) {
                // lambda path starts at 0.0
                lambdaStart = 0.0;
                // electrostaticStart is ignored for MELD.
                electrostaticStart = 0.0;
                // electrostaticLambdaPower is ignored for MELD.
                electrostaticLambdaPower = 1.0;
                // vdW is linearly scaled for MELD.
                vdwSoftcorePower = 1;
                // No softcore offset for MELD.
                vdWSoftcoreAlpha = 0.0;
                // Torsions are linearly scaled for MELD.
                torsionalLambdaPower = 1.0;
                // Only need single-sided dU/dL
                twoSidedFiniteDifference = false;
            }

            // Add Angle Force.
            addAngleForce();
            addInPlaneAngleForce();

            // Add Stretch-Bend Force.
            addStretchBendForce();

            // Add Urey-Bradley Force.
            addUreyBradleyForce();

            // Out-of Plane Bend Force.
            addOutOfPlaneBendForce();

            // Add Torsion Force.
            addTorsionForce();

            // Add Improper Torsion Force.
            addImproperTorsionForce();

            // Add Pi-Torsion Force.
            addPiTorsionForce();

            // Add Torsion-Torsion Force.
            addTorsionTorsionForce();

            // Add coordinate restraints.
            addHarmonicRestraintForce();

            // Add bond restraints.
            addRestraintBonds();

            // Add stretch-torsion coupling terms.
            addStretchTorsionForce();

            // Add angle-torsion coupling terms.
            addAngleTorsionForce();

            if (vdW != null) {
                VanDerWaalsForm vdwForm = vdW.getVDWForm();
                if (vdwForm.vdwType == LENNARD_JONES) {
                    addFixedChargeNonBondedForce();
                } else {
                    // Add vdW Force.
                    addAmoebaVDWForce();

                    // Add Multipole Force.
                    addAmoebaMultipoleForce();
                }
            }

            if (lambdaTerm) {
                logger.info(format("\n Lambda path start:              %6.3f", lambdaStart));
                logger.info(format(" Lambda scales torsions:          %s", torsionLambdaTerm));
                if (torsionLambdaTerm) {
                    logger.info(format(" torsion lambda power:           %6.3f", torsionalLambdaPower));
                }
                logger.info(format(" Lambda scales vdW interactions:  %s", vdwLambdaTerm));
                if (vdwLambdaTerm) {
                    logger.info(format(" van Der Waals alpha:            %6.3f", vdWSoftcoreAlpha));
                    logger.info(format(" van Der Waals lambda power:     %6.3f", vdwSoftcorePower));
                }
                logger.info(format(" Lambda scales electrostatics:    %s", elecLambdaTerm));

                if (elecLambdaTerm) {
                    logger.info(format(" Electrostatics start:           %6.3f", electrostaticStart));
                    logger.info(format(" Electrostatics lambda power:    %6.3f", electrostaticLambdaPower));
                }
                logger.info(format(" Using Meld:                      %s", useMeld));
                if (useMeld) {
                    logger.info(format(" Meld scale factor:              %6.3f", meldScaleFactor));
                }
            }
        }

        public void printLambdaValues() {
            logger.info(format("\n Lambda Values\n Torsion: %6.3f vdW: %6.3f Elec: %6.3f ", lambdaTorsion, lambdaVDW, lambdaElec));
        }

        /**
         * {@inheritDoc}
         */
        public void setLambda(double lambda) {

            // Initially set all lambda values to 1.0.
            lambdaTorsion = 1.0;

            // Applied to softcore vdW forces.
            lambdaVDW = 1.0;

            // Applied to normal electrostatic parameters for alchemical atoms.
            lambdaElec = 1.0;

            if (torsionLambdaTerm) {
                // Multiply torsional potentials by L^2 (dU/dL = 0 at L=0).
                lambdaTorsion = pow(lambda, torsionalLambdaPower);
                if (useMeld) {
                    lambdaTorsion = meldScaleFactor + lambda * (1.0 - meldScaleFactor);
                }
            }

            if (elecLambdaTerm && vdwLambdaTerm) {
                // Lambda effects both vdW and electrostatics.
                if (lambda < electrostaticStart) {
                    // Begin turning vdW on with electrostatics off.
                    lambdaElec = 0.0;
                } else {
                    // Turn electrostatics on during the latter part of the path.
                    double elecWindow = 1.0 - electrostaticStart;
                    lambdaElec = (lambda - electrostaticStart) / elecWindow;
                    lambdaElec = pow(lambdaElec, electrostaticLambdaPower);
                }
                lambdaVDW = lambda;
                if (useMeld) {
                    lambdaElec = sqrt(meldScaleFactor + lambda * (1.0 - meldScaleFactor));
                    lambdaVDW = meldScaleFactor + lambda * (1.0 - meldScaleFactor);
                }
            } else if (vdwLambdaTerm) {
                // Lambda effects vdW, with electrostatics turned off.
                lambdaElec = 0.0;
                lambdaVDW = lambda;
                if (useMeld) {
                    lambdaVDW = meldScaleFactor + lambda * (1.0 - meldScaleFactor);
                }

            } else if (elecLambdaTerm) {
                // Lambda effects electrostatics, but not vdW.
                lambdaElec = lambda;
                if (useMeld) {
                    lambdaElec = sqrt(meldScaleFactor + lambda * (1.0 - meldScaleFactor));
                }
            }
        }

        public void destroy() {
            if (system != null) {
                OpenMM_System_destroy(system);
                system = null;
            }
        }

        PointerByReference getOpenMMSystem() {
            return system;
        }

        void addCOMMRemover(boolean addIfDuplicate) {
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
         * Add a Monte Carlo Barostat to the system.
         *
         * @param targetPressure The target pressure (in atm).
         * @param targetTemp     The target temperature.
         * @param frequency      The frequency to apply the barostat.
         */
        void addMonteCarloBarostat(double targetPressure, double targetTemp, int frequency) {
            if (ommBarostat == null) {
                double pressureInBar = targetPressure * Constants.ATM_TO_BAR;
                ommBarostat = OpenMM_MonteCarloBarostat_create(pressureInBar, targetTemp, frequency);

                CompositeConfiguration properties = molecularAssembly.getProperties();
                if (properties.containsKey("randomseed")) {
                    int randomSeed = properties.getInt("randomseed", 0);
                    logger.info(format(" Setting random seed %d for Monte Carlo Barostat", randomSeed));
                    OpenMM_MonteCarloBarostat_setRandomNumberSeed(ommBarostat, randomSeed);
                }
                OpenMM_System_addForce(system, ommBarostat);
                openMMContext.reinitContext();
                logger.info(format(" Added a Monte Carlo barostat with pressure %6.2f (atm), target temperature %6.2f K and MC move frequency %d.",
                        targetPressure, targetTemp, frequency));
            } else {
                logger.info(" Attempted to add a second barostat to an OpenMM force field!");
            }
        }

        /**
         * <p>
         * calculateDegreesOfFreedom.</p>
         *
         * @return a int.
         */
        int calculateDegreesOfFreedom() {
            int dof = numParticles * 3;
            dof = dof - OpenMM_System_getNumConstraints(system);
            if (commRemover != null) {
                dof -= 3;
            }
            return dof;
        }

        /**
         * Update parameters if the Use flags and/or Lambda value has changed.
         *
         * @param atoms Atoms in this list are considered.
         */
        void updateParameters(Atom[] atoms) {

            if (vdwLambdaTerm) {
                if (fixedChargeNonBondedForce != null) {
                    if (!softcoreCreated) {
                        addCustomNonbondedSoftcoreForce();
                        // Re-initialize the context.
                        openMMContext.reinitContext();
                        softcoreCreated = true;
                    }
                    openMMContext.setParameter("vdw_lambda", lambdaVDW);
                } else if (amoebaVDWForce != null) {
                    openMMContext.setParameter("AmoebaVdwLambda", lambdaVDW);
                    if (softcoreCreated) {
                        // Avoid any updateParametersInContext calls if vdwLambdaTerm is true, but not other alchemical terms.
                        if (!torsionLambdaTerm && !elecLambdaTerm) {
                            return;
                        }
                    } else {
                        softcoreCreated = true;
                    }
                }
            }

            if (torsionLambdaTerm) {
                if (amoebaTorsionForce != null) {
                    updateTorsionForce();
                }
                if (amoebaImproperTorsionForce != null) {
                    updateImproperTorsionForce();
                }
            }

            if (atoms == null || atoms.length == 0) {
                return;
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
         * Add an Andersen thermostat to the system.
         *
         * @param targetTemp    Target temperature in Kelvins.
         * @param collisionFreq Collision frequency in 1/psec
         */
        private void addAndersenThermostat(double targetTemp, double collisionFreq) {
            if (ommThermostat == null) {
                ommThermostat = OpenMM_AndersenThermostat_create(targetTemp, collisionFreq);
                OpenMM_System_addForce(system, ommThermostat);
                logger.info(format(" Added an Andersen thermostat at %6.2fK and collision frequency %6.2f.", targetTemp, collisionFreq));
            } else {
                OpenMM_AndersenThermostat_setDefaultTemperature(ommThermostat, targetTemp);
                OpenMM_AndersenThermostat_setDefaultCollisionFrequency(ommThermostat, collisionFreq);
                logger.info(format(" Updated existing Andersen thermostat at %6.2fK and collision frequency %6.2f.", targetTemp, collisionFreq));
            }
        }

        /**
         * Adds atoms from the molecular assembly to the OpenMM System and reports
         * to the user the number of particles added.
         */
        private void addAtoms() throws Exception {
            numParticles = 0;
            double totalMass = 0.0;
            for (Atom atom : atoms) {
                double mass = atom.getMass();
                totalMass += mass;
                if (mass < 0.0) {
                    throw new Exception(" Atom with mass less than 0.");
                }

                if (mass == 0.0) {
                    logger.info(format(" Atom %s has zero mass.", atom.toString()));
                }

                OpenMM_System_addParticle(system, mass);
                numParticles++;
            }
            logger.log(Level.INFO, format("  Atoms \t\t%6d\t%12.3f", atoms.length, totalMass));
        }

        private void updateAtomMass() {
            int index = 0;
            for (Atom atom : atoms) {
                double mass = 0.0;
                if (atom.isActive()) {
                    mass = atom.getMass();
                }
                OpenMM_System_setParticleMass(system, index++, mass);
            }
        }

        /**
         * Add a bond force to the OpenMM System.
         */
        private void addBondForce() {
            Bond[] bonds = getBonds();
            if (bonds == null || bonds.length < 1) {
                return;
            }
            PointerByReference amoebaBondForce = OpenMM_AmoebaBondForce_create();
            double kParameterConversion = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

            for (Bond bond : bonds) {
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

            int forceGroup = forceField.getInteger("BOND_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaBondForce, forceGroup);
            OpenMM_System_addForce(system, amoebaBondForce);
            logger.log(Level.INFO, format("  Bonds \t\t%6d\t\t%1d", bonds.length, forceGroup));
        }

        /**
         * Add an angle force to the OpenMM System.
         */
        private void addAngleForce() {
            Angle[] angles = getAngles();
            if (angles == null || angles.length < 1) {
                return;
            }
            int nAngles = angles.length;
            List<Angle> normalAngles = new ArrayList<>();
            // Sort all normal angles from in-plane angles
            for (int i = 0; i < nAngles; i++) {
                if (angles[i].getAngleMode() == AngleType.AngleMode.NORMAL) {
                    normalAngles.add(angles[i]);
                }
            }
            nAngles = normalAngles.size();
            if (nAngles < 1) {
                return;
            }

            boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);

            PointerByReference amoebaAngleForce = OpenMM_AmoebaAngleForce_create();
            for (Angle angle : normalAngles) {
                int i1 = angle.getAtom(0).getXyzIndex() - 1;
                int i2 = angle.getAtom(1).getXyzIndex() - 1;
                int i3 = angle.getAtom(2).getXyzIndex() - 1;
                int nh = angle.nh;
                if (isHydrogenAngle(angle) && rigidHydrogenAngles) {
                    logger.info("Not adding angle to AmoebaAngleForce because angle is constrained: " + angle);
                } else {
                    OpenMM_AmoebaAngleForce_addAngle(amoebaAngleForce, i1, i2, i3,
                            angle.angleType.angle[nh], OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant);
                }
            }

            if (angles[0].angleType.angleFunction == AngleFunction.SEXTIC) {
                OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleCubic(amoebaAngleForce, AngleType.cubic);
                OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleQuartic(amoebaAngleForce, AngleType.quartic);
                OpenMM_AmoebaAngleForce_setAmoebaGlobalAnglePentic(amoebaAngleForce, AngleType.quintic);
                OpenMM_AmoebaAngleForce_setAmoebaGlobalAngleSextic(amoebaAngleForce, AngleType.sextic);
            }


            int forceGroup = forceField.getInteger("ANGLE_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaAngleForce, forceGroup);
            OpenMM_System_addForce(system, amoebaAngleForce);
            logger.log(Level.INFO, format("  Angles \t\t%6d\t\t%1d", nAngles, forceGroup));
        }

        /**
         * Add an in-plane angle force to the OpenMM System.
         */
        private void addInPlaneAngleForce() {
            Angle[] angles = getAngles();
            if (angles == null || angles.length < 1) {
                return;
            }
            int nAngles = angles.length;
            List<Angle> inPlaneAngles = new ArrayList<>();
            //Sort all in-plane angles from normal angles
            for (int i = 0; i < nAngles; i++) {
                if (angles[i].getAngleMode() == AngleType.AngleMode.IN_PLANE) {
                    inPlaneAngles.add(angles[i]);
                }
            }
            nAngles = inPlaneAngles.size();
            if (nAngles < 1) {
                return;
            }

            PointerByReference amoebaInPlaneAngleForce = OpenMM_AmoebaInPlaneAngleForce_create();
            for (Angle angle : inPlaneAngles) {
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

            int forceGroup = forceField.getInteger("IN_PLANE_ANGLE_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaInPlaneAngleForce, forceGroup);
            OpenMM_System_addForce(system, amoebaInPlaneAngleForce);

            logger.log(Level.INFO, format("  In-plane Angles \t%6d\t\t%1d", nAngles, forceGroup));
        }

        /**
         * Add a Urey-Bradley force to the OpenMM System.
         */
        private void addUreyBradleyForce() {
            UreyBradley[] ureyBradleys = getUreyBradleys();
            if (ureyBradleys == null || ureyBradleys.length < 1) {
                return;
            }
            PointerByReference amoebaUreyBradleyForce = OpenMM_AmoebaBondForce_create();
            double kParameterConversion = UreyBradleyType.units * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

            for (UreyBradley ureyBradley : ureyBradleys) {
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

            int forceGroup = forceField.getInteger("UREY_BRADLEY_FORCE", 0);

            OpenMM_Force_setForceGroup(amoebaUreyBradleyForce, forceGroup);
            OpenMM_System_addForce(system, amoebaUreyBradleyForce);
            logger.log(Level.INFO, format("  Urey-Bradleys \t%6d\t\t%1d", ureyBradleys.length, forceGroup));
        }

        /**
         * Add an out-of-plane bend force to the OpenMM System.
         */
        private void addOutOfPlaneBendForce() {
            OutOfPlaneBend[] outOfPlaneBends = getOutOfPlaneBends();
            if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
                return;
            }
            PointerByReference amoebaOutOfPlaneBendForce = OpenMM_AmoebaOutOfPlaneBendForce_create();

            for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
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

            int forceGroup = forceField.getInteger("OUT_OF_PLANE_BEND_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaOutOfPlaneBendForce, forceGroup);
            OpenMM_System_addForce(system, amoebaOutOfPlaneBendForce);
            logger.log(Level.INFO, format("  Out-of-Plane Bends \t%6d\t\t%1d", outOfPlaneBends.length, forceGroup));
        }

        /**
         * Add a stretch-bend force to the OpenMM System.
         */
        private void addStretchBendForce() {
            StretchBend[] stretchBends = getStretchBends();
            if (stretchBends == null || stretchBends.length < 1) {
                return;
            }

            PointerByReference amoebaStretchBendForce = OpenMM_AmoebaStretchBendForce_create();
            for (StretchBend stretchBend : stretchBends) {
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

            int forceGroup = forceField.getInteger("STRETCH_BEND_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaStretchBendForce, forceGroup);
            OpenMM_System_addForce(system, amoebaStretchBendForce);
            logger.log(Level.INFO, format("  Stretch-Bends \t%6d\t\t%1d", stretchBends.length, forceGroup));
        }

        /**
         * Add a torsion force to the OpenMM System.
         */
        private void addTorsionForce() {
            Torsion[] torsions = getTorsions();
            if (torsions == null || torsions.length < 1) {
                return;
            }

            amoebaTorsionForce = OpenMM_PeriodicTorsionForce_create();
            for (Torsion torsion : torsions) {
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

            int fGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaTorsionForce, fGroup);
            OpenMM_System_addForce(system, amoebaTorsionForce);

            logger.log(Level.INFO, format("  Torsions \t\t%6d\t\t%1d", torsions.length, fGroup));
        }

        /**
         * Add an improper-torsion force to the OpenMM System.
         */
        private void addImproperTorsionForce() {
            ImproperTorsion[] improperTorsions = getImproperTorsions();
            if (improperTorsions == null || improperTorsions.length < 1) {
                return;
            }

            amoebaImproperTorsionForce = OpenMM_PeriodicTorsionForce_create();
            for (ImproperTorsion improperTorsion : improperTorsions) {
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

            int forceGroup = forceField.getInteger("IMPROPER_TORSION_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaImproperTorsionForce, forceGroup);
            OpenMM_System_addForce(system, amoebaImproperTorsionForce);

            logger.log(Level.INFO, format("  Improper Torsions \t%6d\t\t%1d", improperTorsions.length, forceGroup));
        }

        /**
         * Add a Pi-Torsion force to the OpenMM System.
         */
        private void addPiTorsionForce() {
            PiOrbitalTorsion[] piOrbitalTorsions = getPiOrbitalTorsions();
            if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
                return;
            }

            PointerByReference amoebaPiTorsionForce = OpenMM_AmoebaPiTorsionForce_create();
            double units = PiTorsionType.units;
            for (PiOrbitalTorsion piOrbitalTorsion : piOrbitalTorsions) {
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

            int forceGroup = forceField.getInteger("PI_ORBITAL_TORSION_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaPiTorsionForce, forceGroup);
            OpenMM_System_addForce(system, amoebaPiTorsionForce);
            logger.log(Level.INFO, format("  Pi-Orbital Torsions  \t%6d\t\t%1d", piOrbitalTorsions.length, forceGroup));
        }

        /**
         * Add a Torsion-Torsion force to the OpenMM System.
         */
        private void addTorsionTorsionForce() {
            TorsionTorsion[] torsionTorsions = getTorsionTorsions();
            if (torsionTorsions == null || torsionTorsions.length < 1) {
                return;
            }

            // Load the torsion-torsions.
            int nTypes = 0;
            LinkedHashMap<String, TorsionTorsionType> torTorTypes = new LinkedHashMap<>();
            PointerByReference amoebaTorsionTorsionForce = OpenMM_AmoebaTorsionTorsionForce_create();
            for (TorsionTorsion torsionTorsion : torsionTorsions) {
                int ia = torsionTorsion.getAtom(0).getXyzIndex() - 1;
                int ib = torsionTorsion.getAtom(1).getXyzIndex() - 1;
                int ic = torsionTorsion.getAtom(2).getXyzIndex() - 1;
                int id = torsionTorsion.getAtom(3).getXyzIndex() - 1;
                int ie = torsionTorsion.getAtom(4).getXyzIndex() - 1;

                TorsionTorsionType torsionTorsionType = torsionTorsion.torsionTorsionType;
                String key = torsionTorsionType.getKey();

                // Check if the TorTor parameters have already been added to the Hash.
                int gridIndex = 0;
                if (torTorTypes.containsKey(key)) {

                    // If the TorTor has been added, get its (ordered) index in the Hash.
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
                    // Add the new TorTor.
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

            // Load the Torsion-Torsion parameters.
            PointerByReference values = OpenMM_DoubleArray_create(6);
            int gridIndex = 0;
            for (String key : torTorTypes.keySet()) {
                TorsionTorsionType torTorType = torTorTypes.get(key);
                int nx = torTorType.nx;
                int ny = torTorType.ny;
                double[] tx = torTorType.tx;
                double[] ty = torTorType.ty;
                double[] f = torTorType.energy;
                double[] dx = torTorType.dx;
                double[] dy = torTorType.dy;
                double[] dxy = torTorType.dxy;

                // Create the 3D grid.
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
                    OpenMM_DoubleArray_set(values, addIndex, OpenMM_KJPerKcal * dxy[j]);
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

            int forceGroup = forceField.getInteger("TORSION_TORSION_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(amoebaTorsionTorsionForce, forceGroup);
            OpenMM_System_addForce(system, amoebaTorsionTorsionForce);
            logger.log(Level.INFO, format("  Torsion-Torsions  \t%6d\t\t%1d", torsionTorsions.length, forceGroup));
        }

        /**
         * Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic
         * acid model).
         */
        private void addStretchTorsionForce() {
            StretchTorsion[] stretchTorsions = getStretchTorsions();
            if (stretchTorsions == null || stretchTorsions.length < 1) {
                return;
            }

            PointerByReference stretchTorsionForce = OpenMM_CustomCompoundBondForce_create(4, StretchTorsion.stretchTorsionForm());
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

            for (StretchTorsion strTors : stretchTorsions) {
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

            int forceGroup = forceField.getInteger("STRETCH_TORSION_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(stretchTorsionForce, forceGroup);
            OpenMM_System_addForce(system, stretchTorsionForce);

            logger.log(Level.INFO, format("  Stretch-Torsions  \t%6d\t\t%1d", stretchTorsions.length, forceGroup));

        }

        /**
         * Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic
         * acid model).
         */
        private void addAngleTorsionForce() {
            AngleTorsion[] angleTorsions = getAngleTorsions();
            if (angleTorsions == null || angleTorsions.length < 1) {
                return;
            }

            PointerByReference angleTorsionForce = OpenMM_CustomCompoundBondForce_create(4, AngleTorsion.angleTorsionForm());
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
            for (AngleTorsion angleTorsion : angleTorsions) {
                double[] constants = angleTorsion.getConstants();
                PointerByReference atorsParams = OpenMM_DoubleArray_create(0);
                for (int m = 0; m < 2; m++) {
                    for (int n = 0; n < 3; n++) {
                        int index = (3 * m) + n;
                        double kmn = constants[index] * unitConv;
                        OpenMM_DoubleArray_append(atorsParams, kmn);
                    }
                }

                Atom[] atoms = angleTorsion.getAtomArray(true);

                // One thing that concerns me is whether it's correct to get angle[0] instead of angle[num hydrogens].
                // This is the way it is in FFX, but that may be a bug.

                OpenMM_DoubleArray_append(atorsParams, angleTorsion.angleType1.angle[0] * OpenMM_RadiansPerDegree);
                OpenMM_DoubleArray_append(atorsParams, angleTorsion.angleType2.angle[0] * OpenMM_RadiansPerDegree);

                PointerByReference atorsParticles = OpenMM_IntArray_create(0);
                for (int i = 0; i < 4; i++) {
                    OpenMM_IntArray_append(atorsParticles, atoms[i].getXyzIndex() - 1);
                }

                OpenMM_CustomCompoundBondForce_addBond(angleTorsionForce, atorsParticles, atorsParams);
                OpenMM_DoubleArray_destroy(atorsParams);
                OpenMM_IntArray_destroy(atorsParticles);
            }

            int forceGroup = forceField.getInteger("ANGLE_TORSION_FORCE_GROUP", 0);

            OpenMM_Force_setForceGroup(angleTorsionForce, forceGroup);
            OpenMM_System_addForce(system, angleTorsionForce);

            logger.log(Level.INFO, format("  Angle-Torsions  \t%6d\t\t%1d", angleTorsions.length, forceGroup));

        }

        /**
         * Uses arithmetic mean to define sigma and geometric mean for epsilon.
         */
        private void addFixedChargeNonBondedForce() {
            VanDerWaals vdW = getVdwNode();
            if (vdW == null) {
                return;
            }

        /*
          Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
          for epsilon is supported.
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

            // OpenMM vdW force requires a diameter (i.e. not radius).
            double radScale = 1.0;
            if (vdwForm.radiusSize == RADIUS) {
                radScale = 2.0;
            }

            // OpenMM vdw force requires atomic sigma values (i.e. not r-min).
            if (vdwForm.radiusType == R_MIN) {
                radScale /= 1.122462048309372981;
            }

            // Add particles.
            for (Atom atom : atoms) {
                VDWType vdwType = atom.getVDWType();
                double sigma = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
                double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
                double charge = 0.0;
                MultipoleType multipoleType = atom.getMultipoleType();
                if (multipoleType != null && atom.getElectrostatics()) {
                    charge = multipoleType.charge;
                }
                OpenMM_NonbondedForce_addParticle(fixedChargeNonBondedForce, charge, sigma, eps);
            }

            // Define 1-4 scale factors.
            double lj14Scale = vdwForm.getScale14();
            double coulomb14Scale = 1.0 / 1.2;
            ParticleMeshEwald pme = getPmeNode();
            if (pme != null) {
                coulomb14Scale = pme.getScale14();
            }
            Bond[] bonds = getBonds();
            if (bonds != null && bonds.length > 0) {
                PointerByReference bondArray = OpenMM_BondArray_create(0);
                for (Bond bond : bonds) {
                    int i1 = bond.getAtom(0).getXyzIndex() - 1;
                    int i2 = bond.getAtom(1).getXyzIndex() - 1;
                    OpenMM_BondArray_append(bondArray, i1, i2);
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

            Crystal crystal = getCrystal();
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

            int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);
            int pmeGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
            if (forceGroup != pmeGroup) {
                logger.severe(String.format(" ERROR: VDW-FORCE-GROUP is %d while PME-FORCE-GROUP is %d. "
                        + "This is invalid for fixed-charge force fields with combined nonbonded forces.", forceGroup, pmeGroup));
            }

            OpenMM_Force_setForceGroup(fixedChargeNonBondedForce, forceGroup);
            OpenMM_System_addForce(system, fixedChargeNonBondedForce);

            logger.log(Level.INFO, format("  Fixed charge non-bonded force \t%1d", forceGroup));

            GeneralizedKirkwood gk = getGK();
            if (gk != null) {
                addCustomGBForce();
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
        private void addCustomNonbondedSoftcoreForce() {
            VanDerWaals vdW = getVdwNode();
            if (vdW == null) {
                return;
            }

        /*
          Only 6-12 LJ with arithmetic mean to define sigma and geometric mean
          for epsilon is supported.
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

            // Sterics mixing rules.
            String stericsMixingRules = " epsilon = sqrt(epsilon1*epsilon2);";
            stericsMixingRules += " rmin = 0.5 * (sigma1 + sigma2) * 1.122462048309372981;";

            // Softcore Lennard-Jones, with a form equivalent to that used in FFX VanDerWaals class.
            String stericsEnergyExpression = "(vdw_lambda^beta)*epsilon*x*(x-2.0);";
            // Effective softcore distance for sterics.
            stericsEnergyExpression += " x = 1.0 / (alpha*(1.0-vdw_lambda)^2.0 + (r/rmin)^6.0);";
            // Define energy expression for sterics.
            String energyExpression = stericsEnergyExpression + stericsMixingRules;

            PointerByReference fixedChargeSoftcore = OpenMM_CustomNonbondedForce_create(energyExpression);

            // Get the Alpha and Beta constants from the VanDerWaals instance.
            double alpha = vdWSoftcoreAlpha;
            double beta = vdwSoftcorePower;

            OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "vdw_lambda", 1.0);
            OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "alpha", alpha);
            OpenMM_CustomNonbondedForce_addGlobalParameter(fixedChargeSoftcore, "beta", beta);
            OpenMM_CustomNonbondedForce_addPerParticleParameter(fixedChargeSoftcore, "sigma");
            OpenMM_CustomNonbondedForce_addPerParticleParameter(fixedChargeSoftcore, "epsilon");

            // Add particles.
            PointerByReference alchemicalGroup = OpenMM_IntSet_create();
            PointerByReference nonAlchemicalGroup = OpenMM_IntSet_create();
            DoubleByReference charge = new DoubleByReference();
            DoubleByReference sigma = new DoubleByReference();
            DoubleByReference eps = new DoubleByReference();

            int index = 0;
            for (Atom atom : atoms) {
                if (atom.applyLambda()) {
                    OpenMM_IntSet_insert(alchemicalGroup, index);
                } else {
                    OpenMM_IntSet_insert(nonAlchemicalGroup, index);
                }

                OpenMM_NonbondedForce_getParticleParameters(fixedChargeNonBondedForce, index, charge, sigma, eps);
                double sigmaValue = sigma.getValue();
                double epsValue = eps.getValue();

                // Handle cases where sigma is 0.0; for example Amber99 tyrosine hydrogen atoms.
                if (sigmaValue == 0.0) {
                    sigmaValue = 1.0;
                    epsValue = 0.0;
                }

                PointerByReference particleParameters = OpenMM_DoubleArray_create(0);
                OpenMM_DoubleArray_append(particleParameters, sigmaValue);
                OpenMM_DoubleArray_append(particleParameters, epsValue);
                OpenMM_CustomNonbondedForce_addParticle(fixedChargeSoftcore, particleParameters);
                OpenMM_DoubleArray_destroy(particleParameters);

                index++;
            }

            OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup, alchemicalGroup);
            OpenMM_CustomNonbondedForce_addInteractionGroup(fixedChargeSoftcore, alchemicalGroup, nonAlchemicalGroup);
            OpenMM_IntSet_destroy(alchemicalGroup);
            OpenMM_IntSet_destroy(nonAlchemicalGroup);

            Crystal crystal = getCrystal();
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
                    logger.info(format(" Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.", cut, off));
                    cut *= 0.99;
                }
            }

            OpenMM_CustomNonbondedForce_setCutoffDistance(fixedChargeSoftcore, OpenMM_NmPerAngstrom * off);
            OpenMM_CustomNonbondedForce_setUseSwitchingFunction(fixedChargeSoftcore, OpenMM_True);
            OpenMM_CustomNonbondedForce_setSwitchingDistance(fixedChargeSoftcore, OpenMM_NmPerAngstrom * cut);

            // Add energy parameter derivative
            // OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(fixedChargeSoftcore, "vdw_lambda");

            int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);

            OpenMM_Force_setForceGroup(fixedChargeSoftcore, forceGroup);
            OpenMM_System_addForce(system, fixedChargeSoftcore);

            // Alchemical with Alchemical could be either softcore or normal interactions (softcore here).
            PointerByReference alchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(stericsEnergyExpression);

            // Non-Alchemical with Alchemical is essentially always softcore.
            PointerByReference nonAlchemicalAlchemicalStericsForce = OpenMM_CustomBondForce_create(stericsEnergyExpression);

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
            int[][] torsionMask = vdW.getTorsionMask();

            for (int i = 0; i < range; i++) {
                OpenMM_NonbondedForce_getExceptionParameters(fixedChargeNonBondedForce, i, atomi, atomj, charge, sigma, eps);

                // Omit both Exclusions (1-2, 1-3) and Exceptions (scaled 1-4) from the CustomNonbondedForce.
                OpenMM_CustomNonbondedForce_addExclusion(fixedChargeSoftcore, atomi.getValue(), atomj.getValue());

                // Deal with scaled 1-4 torsions using the CustomBondForce
                int[] maskI = torsionMask[atomi.getValue()];
                int jID = atomj.getValue();
                boolean epsException = false;

                for (int mask : maskI) {
                    if (mask == jID) {
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

            OpenMM_Force_setForceGroup(alchemicalAlchemicalStericsForce, forceGroup);
            OpenMM_Force_setForceGroup(nonAlchemicalAlchemicalStericsForce, forceGroup);

            // OpenMM_CustomBondForce_addEnergyParameterDerivative(alchemicalAlchemicalStericsForce, "vdw_lambda");
            // OpenMM_CustomBondForce_addEnergyParameterDerivative(nonAlchemicalAlchemicalStericsForce, "vdw_lambda");

            OpenMM_System_addForce(system, alchemicalAlchemicalStericsForce);
            OpenMM_System_addForce(system, nonAlchemicalAlchemicalStericsForce);

            logger.log(Level.INFO, format("  Added fixed charge softcore force \t%d", forceGroup));
            logger.log(Level.INFO, format("   Alpha = %8.6f and beta = %8.6f", alpha, beta));
        }

        /**
         * Add a custom GB force to the OpenMM System.
         */
        private void addCustomGBForce() {
            GeneralizedKirkwood gk = getGK();
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

            // Particle pair term is the generalized Born cross term.
            OpenMM_CustomGBForce_addEnergyTerm(customGBForce,
                    "-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                            + "f=sqrt(r^2+B1*B2*exp(-r^2/(2.455*B1*B2)))",
                    OpenMM_CustomGBForce_ParticlePair);

            double[] baseRadii = gk.getBaseRadii();
            double[] overlapScale = gk.getOverlapScale();
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

            int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
            OpenMM_Force_setForceGroup(customGBForce, forceGroup);
            OpenMM_System_addForce(system, customGBForce);

            logger.log(Level.INFO, format("  Custom generalized Born force \t%d", forceGroup));
        }

        /**
         * Add an AMOEBA van der Waals force to the OpenMM System.
         */
        private void addAmoebaVDWForce() {
            VanDerWaals vdW = getVdwNode();
            if (vdW == null) {
                return;
            }

            amoebaVDWForce = OpenMM_AmoebaVdwForce_create();

            VanDerWaalsForm vdwForm = vdW.getVDWForm();
            NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
            Crystal crystal = getCrystal();

            double radScale = 1.0;
            if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
                radScale = 0.5;
            }

            int[] ired = vdW.getReductionIndex();
            int nAtoms = atoms.length;
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                VDWType vdwType = atom.getVDWType();
                int isAlchemical = atom.applyLambda() ? 1 : 0;
                OpenMM_AmoebaVdwForce_addParticle(amoebaVDWForce,
                        ired[i], OpenMM_NmPerAngstrom * vdwType.radius * radScale,
                        OpenMM_KJPerKcal * vdwType.wellDepth,
                        vdwType.reductionFactor, isAlchemical);
            }

            OpenMM_AmoebaVdwForce_setCutoffDistance(amoebaVDWForce, nonbondedCutoff.off * OpenMM_NmPerAngstrom);
            if (vdW.getDoLongRangeCorrection()) {
                OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_Boolean.OpenMM_True);
            } else {
                OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_Boolean.OpenMM_False);
            }

            if (crystal.aperiodic()) {
                OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
                        OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_NoCutoff);
            } else {
                OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
                        OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_CutoffPeriodic);
            }

            if (vdwLambdaTerm) {
                OpenMM_AmoebaVdwForce_setAlchemicalMethod(amoebaVDWForce,
                        AmoebaOpenMMLibrary.OpenMM_AmoebaVdwForce_AlchemicalMethod.OpenMM_AmoebaVdwForce_Decouple);
                OpenMM_AmoebaVdwForce_setSoftcoreAlpha(amoebaVDWForce, vdWSoftcoreAlpha);
                OpenMM_AmoebaVdwForce_setSoftcorePower(amoebaVDWForce, (int) vdwSoftcorePower);
            }

            int[][] bondMask = vdW.getBondMask();
            int[][] angleMask = vdW.getAngleMask();

            // Create exclusion lists.
            PointerByReference exclusions = OpenMM_IntArray_create(0);
            for (int i = 0; i < nAtoms; i++) {
                OpenMM_IntArray_append(exclusions, i);
                final int[] bondMaski = bondMask[i];
                for (int value : bondMaski) {
                    OpenMM_IntArray_append(exclusions, value);
                }
                final int[] angleMaski = angleMask[i];
                for (int value : angleMaski) {
                    OpenMM_IntArray_append(exclusions, value);
                }
                OpenMM_AmoebaVdwForce_setParticleExclusions(amoebaVDWForce, i, exclusions);
                OpenMM_IntArray_resize(exclusions, 0);
            }
            OpenMM_IntArray_destroy(exclusions);

            int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);
            OpenMM_Force_setForceGroup(amoebaVDWForce, forceGroup);
            OpenMM_System_addForce(system, amoebaVDWForce);

            logger.log(Level.INFO, format("  AMOEBA van der Waals force \t\t%d", forceGroup));
        }

        /**
         * Add an AMOEBA polarizable multipole force to the OpenMM System.
         */
        private void addAmoebaMultipoleForce() {
            ParticleMeshEwald pme = getPmeNode();
            if (pme == null) {
                return;
            }

            int[][] axisAtom = pme.getAxisAtoms();
            double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
            double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom
                    * OpenMM_NmPerAngstrom;
            double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

            amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();

            double polarScale = 1.0;
            ParticleMeshEwald.SCFAlgorithm scfAlgorithm = null;

            if (pme.getPolarizationType() != Polarization.MUTUAL) {
                OpenMM_AmoebaMultipoleForce_setPolarizationType(amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Direct);
                if (pme.getPolarizationType() == Polarization.NONE) {
                    polarScale = 0.0;
                }
            } else {
                String algorithm = forceField.getString("SCF_ALGORITHM", "CG");
                try {
                    algorithm = algorithm.replaceAll("-", "_").toUpperCase();
                    scfAlgorithm = ParticleMeshEwald.SCFAlgorithm.valueOf(algorithm);
                } catch (Exception e) {
                    scfAlgorithm = ParticleMeshEwald.SCFAlgorithm.CG;
                }

                switch (scfAlgorithm) {
                    case EPT:
                        /*
                         * Citation:
                         * Simmonett, A. C.;  Pickard, F. C. t.;  Shao, Y.;  Cheatham, T. E., 3rd; Brooks, B. R., Efficient treatment of induced dipoles. The Journal of chemical physics 2015, 143 (7), 074115-074115.
                         */
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

            int nAtoms = atoms.length;
            for (int i = 0; i < nAtoms; i++) {
                Atom atom = atoms[i];
                MultipoleType multipoleType = atom.getMultipoleType();
                PolarizeType polarType = atom.getPolarizeType();

                // Define the frame definition.
                int axisType;
                switch (multipoleType.frameDefinition) {
                    default:
                    case NONE:
                        axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
                        break;
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
                    case THREEFOLD:
                        axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
                        break;
                }

                double useFactor = 1.0;
                if (!atoms[i].getUse() || !atoms[i].getElectrostatics()) {
                    useFactor = 0.0;
                }

                double lambdaScale = lambdaElec; // Should be 1.0 at this point.
                if (!atom.applyLambda()) {
                    lambdaScale = 1.0;
                }

                useFactor *= lambdaScale;

                // Load local multipole coefficients.
                for (int j = 0; j < 3; j++) {
                    OpenMM_DoubleArray_set(dipoles, j,
                            multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
                }
                int l = 0;
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        OpenMM_DoubleArray_set(quadrupoles, l++,
                                multipoleType.quadrupole[j][k] * quadrupoleConversion * useFactor / 3.0);
                    }
                }

                int zaxis = 1;
                int xaxis = 1;
                int yaxis = 1;
                int[] refAtoms = axisAtom[i];
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

                double charge = multipoleType.charge * useFactor;

                // Add the multipole.
                OpenMM_AmoebaMultipoleForce_addMultipole(amoebaMultipoleForce,
                        charge, dipoles, quadrupoles,
                        axisType, zaxis, xaxis, yaxis,
                        polarType.thole,
                        polarType.pdamp * dampingFactorConversion,
                        polarType.polarizability * polarityConversion * polarScale);
            }
            OpenMM_DoubleArray_destroy(dipoles);
            OpenMM_DoubleArray_destroy(quadrupoles);

            Crystal crystal = getCrystal();
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

            int[][] ip11 = pme.getPolarization11();

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
                        if (!list12.contains(index) && !list13.contains(index)) {
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
                    if (!list12.contains(index) && !list13.contains(index) && !list14.contains(index)) {
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

                // AMOEBA does not scale between 1-2, 1-3, etc polarization groups.
            }

            OpenMM_IntArray_destroy(covalentMap);

            int forceGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
            OpenMM_Force_setForceGroup(amoebaMultipoleForce, forceGroup);
            OpenMM_System_addForce(system, amoebaMultipoleForce);

            logger.log(Level.INFO, format("  AMOEBA polarizable multipole force \t%d", forceGroup));

            GeneralizedKirkwood gk = getGK();
            if (gk != null) {
                addGKForce();
            }

            if (scfAlgorithm == ParticleMeshEwald.SCFAlgorithm.EPT) {
                logger.info("   Using extrapolated perturbation theory for polarization energy.");
            }
        }

        /**
         * Add a Generalized Kirkwood force to the OpenMM System.
         */
        private void addGKForce() {
            GeneralizedKirkwood gk = getGK();

            amoebaGeneralizedKirkwoodForce = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
            OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(amoebaGeneralizedKirkwoodForce, gk.getSolventPermittivity());
            OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(amoebaGeneralizedKirkwoodForce, 1.0);

            double[] overlapScale = gk.getOverlapScale();
            double[] baseRadii = gk.getBaseRadii();
            int nAtoms = atoms.length;
            for (int i = 0; i < nAtoms; i++) {
                MultipoleType multipoleType = atoms[i].getMultipoleType();
                OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle(amoebaGeneralizedKirkwoodForce,
                        multipoleType.charge, OpenMM_NmPerAngstrom * baseRadii[i], overlapScale[i]);
            }

            OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(amoebaGeneralizedKirkwoodForce,
                    gk.getProbeRadius() * OpenMM_NmPerAngstrom);

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

            int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
            OpenMM_Force_setForceGroup(amoebaGeneralizedKirkwoodForce, forceGroup);
            OpenMM_System_addForce(system, amoebaGeneralizedKirkwoodForce);

            logger.log(Level.INFO, format("  Generalized Kirkwood force \t\t%d", forceGroup));

            switch (nonpolar) {
                case CAV_DISP:
                case BORN_CAV_DISP:
                    addWCAForce();
                    break;
                case CAV:
                case HYDROPHOBIC_PMF:
                case BORN_SOLV:
                case NONE:
                default:
                    // WCA force is not being used.
            }
        }

        /**
         * Add a nonpolar Weeks-Chandler-Andersen dispersion force to the OpenMM System.
         */
        private void addWCAForce() {

            double epso = 0.1100;
            double epsh = 0.0135;
            double rmino = 1.7025;
            double rminh = 1.3275;
            double awater = 0.033428;
            double slevy = 1.0;
            double dispoff = 0.26;
            double shctd = 0.81;

            VanDerWaals vdW = getVdwNode();
            VanDerWaalsForm vdwForm = vdW.getVDWForm();
            double radScale = 1.0;
            if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
                radScale = 0.5;
            }

            amoebaWcaDispersionForce = OpenMM_AmoebaWcaDispersionForce_create();

            for (Atom atom : atoms) {
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

            int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);

            OpenMM_Force_setForceGroup(amoebaWcaDispersionForce, forceGroup);
            OpenMM_System_addForce(system, amoebaWcaDispersionForce);

            logger.log(Level.INFO, format("  WCA dispersion force \t\t\t%d", forceGroup));
        }

        /**
         * Adds harmonic restraints (CoordRestraint objects) to OpenMM as a custom
         * external force.
         * <p>
         * TODO: Make robust to flat-bottom restraints.
         */
        private void addHarmonicRestraintForce() {
            int nRestraints = getCoordRestraints().size();

            int forceGroup = forceField.getInteger("COORD_RESTRAINT_FORCE_GROUP", 0);

            for (CoordRestraint coordRestraint : getCoordRestraints()) {
                double forceConstant = coordRestraint.getForceConstant();
                forceConstant *= OpenMM_KJPerKcal;
                forceConstant *= (OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm);
                Atom[] restAtoms = coordRestraint.getAtoms();
                int nRestAts = coordRestraint.getNumAtoms();
                double[][] oCoords = coordRestraint.getOriginalCoordinates();
                for (int i = 0; i < nRestAts; i++) {
                    oCoords[i][0] *= OpenMM_NmPerAngstrom;
                    oCoords[i][1] *= OpenMM_NmPerAngstrom;
                    oCoords[i][2] *= OpenMM_NmPerAngstrom;
                }

                PointerByReference theRestraint = OpenMM_CustomExternalForce_create("k*periodicdistance(x,y,z,x0,y0,z0)^2");
                OpenMM_CustomExternalForce_addGlobalParameter(theRestraint, "k", forceConstant);
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

                OpenMM_Force_setForceGroup(theRestraint, forceGroup);
                OpenMM_System_addForce(system, theRestraint);
            }

            if (nRestraints > 0) {
                logger.log(Level.INFO, format("  Harmonic restraint force \t%6d\t%d", nRestraints, forceGroup));
            }
        }

        /**
         * Adds restraint bonds, if any.
         */
        private void addRestraintBonds() {
            List<RestraintBond> restraintBonds = getRestraintBonds();
            if (restraintBonds == null || restraintBonds.isEmpty()) {
                return;
            }

            int forceGroup = forceField.getInteger("BOND_RESTRAINT_FORCE_GROUP", 0);

            // OpenMM's HarmonicBondForce class uses k, not 1/2*k as does FFX.
            double kParameterConversion = BondType.units * 2.0 * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

            // Map from bond functional forms to the restraint-bonds using that functional form.
            Map<BondType.BondFunction, PointerByReference> restraintForces = new HashMap<>();

            for (RestraintBond restraintBond : getRestraintBonds()) {
                PointerByReference theForce;
                BondType bondType = restraintBond.bondType;
                BondType.BondFunction bondFunction = bondType.bondFunction;
                if (restraintForces.containsKey(bondFunction)) {
                    theForce = restraintForces.get(bondFunction);
                } else {
                    theForce = OpenMM_CustomBondForce_create(bondFunction.toMathematicalForm());
                    OpenMM_CustomBondForce_addPerBondParameter(theForce, "k");
                    OpenMM_CustomBondForce_addPerBondParameter(theForce, "r0");
                    if (bondFunction.hasFlatBottom()) {
                        OpenMM_CustomBondForce_addPerBondParameter(theForce, "fb");
                    }

                    // Wholly untested code.
                    switch (bondFunction) {
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

                    OpenMM_Force_setForceGroup(theForce, forceGroup);
                    OpenMM_System_addForce(system, theForce);
                    restraintForces.put(bondFunction, theForce);
                }

                double forceConst = bondType.forceConstant * kParameterConversion;
                double equilDist = bondType.distance * OpenMM_NmPerAngstrom;
                Atom[] ats = restraintBond.getAtomArray();
                int at1 = ats[0].getXyzIndex() - 1;
                int at2 = ats[1].getXyzIndex() - 1;

                PointerByReference bondParams = OpenMM_DoubleArray_create(0);
                OpenMM_DoubleArray_append(bondParams, forceConst);
                OpenMM_DoubleArray_append(bondParams, equilDist);
                if (bondFunction.hasFlatBottom()) {
                    OpenMM_DoubleArray_append(bondParams, bondType.flatBottomRadius * OpenMM_NmPerAngstrom);
                }
                OpenMM_CustomBondForce_addBond(theForce, at1, at2, bondParams);
                OpenMM_DoubleArray_destroy(bondParams);
            }

            logger.log(Level.INFO, format("  Restraint bonds force \t%6d\t%d", restraintBonds.size(), forceGroup));
        }

        /**
         * Updates the AMOEBA van der Waals force for changes in Use flags or Lambda.
         *
         * @param atoms Array of atoms to update.
         */
        private void updateAmoebaVDWForce(Atom[] atoms) {
            VanDerWaals vdW = getVdwNode();
            VanDerWaalsForm vdwForm = vdW.getVDWForm();
            double radScale = 1.0;
            if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
                radScale = 0.5;
            }

            int[] ired = vdW.getReductionIndex();
            for (Atom atom : atoms) {
                int index = atom.getXyzIndex() - 1;
                VDWType vdwType = atom.getVDWType();
                double useFactor = 1.0;
                if (!atom.getUse()) {
                    useFactor = 0.0;
                }

                int isAlchemical = atom.applyLambda() ? 1 : 0;
                double eps = OpenMM_KJPerKcal * vdwType.wellDepth * useFactor;
                OpenMM_AmoebaVdwForce_setParticleParameters(amoebaVDWForce, index, ired[index],
                        OpenMM_NmPerAngstrom * vdwType.radius * radScale, eps,
                        vdwType.reductionFactor, isAlchemical);
            }

            if (openMMContext.context != null) {
                OpenMM_AmoebaVdwForce_updateParametersInContext(amoebaVDWForce, openMMContext.context);
            }
        }

        /**
         * Updates the fixed-charge non-bonded force for changes in Use flags or Lambda.
         *
         * @param atoms Array of atoms to update.
         */
        private void updateFixedChargeNonBondedForce(Atom[] atoms) {
            VanDerWaals vdW = getVdwNode();
            // Only 6-12 LJ with arithmetic mean to define sigma and geometric mean for epsilon is supported.
            VanDerWaalsForm vdwForm = vdW.getVDWForm();
            if (vdwForm.vdwType != LENNARD_JONES
                    || vdwForm.radiusRule != ARITHMETIC
                    || vdwForm.epsilonRule != GEOMETRIC) {
                logger.log(Level.SEVERE, " Unsupported van der Waals functional form.");
                return;
            }

            // OpenMM vdW force requires a diameter (i.e. not radius).
            double radScale = 1.0;
            if (vdwForm.radiusSize == RADIUS) {
                radScale = 2.0;
            }

            // OpenMM vdw force requires atomic sigma values (i.e. not r-min).
            if (vdwForm.radiusType == R_MIN) {
                radScale /= 1.122462048309372981;
            }

            // Update parameters.
            for (Atom atom : atoms) {
                int index = atom.getXyzIndex() - 1;
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
                }

                if (!atom.getUse()) {
                    eps = 0.0;
                    charge = 0.0;
                }

                OpenMM_NonbondedForce_setParticleParameters(fixedChargeNonBondedForce, index, charge, sigma, eps);
            }

            // Update Exceptions.
            IntByReference particle1 = new IntByReference();
            IntByReference particle2 = new IntByReference();
            DoubleByReference chargeProd = new DoubleByReference();
            DoubleByReference sigma = new DoubleByReference();
            DoubleByReference eps = new DoubleByReference();

            int numExceptions = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);
            for (int i = 0; i < numExceptions; i++) {

                // Only update exceptions.
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

            /*
             Note that the minimum epsilon value cannot be zero, or OpenMM may
             report an error that the number of Exceptions has changed.
             */
                double minEpsilon = 1.0e-12;
                double lambdaValue = lambdaElec;

                if (lambdaValue < minEpsilon) {
                    lambdaValue = minEpsilon;
                }

                if (atom1.applyLambda()) {
                    qq *= lambdaValue;
                    if (vdwLambdaTerm) {
                        epsilon = minEpsilon;
                    }
                }
                if (atom2.applyLambda()) {
                    qq *= lambdaValue;
                    if (vdwLambdaTerm) {
                        epsilon = minEpsilon;
                    }
                }
                if (!atom1.getUse() || !atom2.getUse()) {
                    qq = minEpsilon;
                    epsilon = minEpsilon;
                }
                OpenMM_NonbondedForce_setExceptionParameters(fixedChargeNonBondedForce, i,
                        i1, i2, qq, sigma.getValue(), epsilon);
            }

            if (openMMContext.context != null) {
                OpenMM_NonbondedForce_updateParametersInContext(fixedChargeNonBondedForce, openMMContext.context);
            }
        }

        /**
         * Updates the custom GB force for changes in Use flags or Lambda.
         *
         * @param atoms Array of atoms to update.
         */
        private void updateCustomGBForce(Atom[] atoms) {
            GeneralizedKirkwood gk = getGK();
            double[] baseRadii = gk.getBaseRadii();
            double[] overlapScale = gk.getOverlapScale();
            boolean nea = gk.getNativeEnvironmentApproximation();

            double sTens = 0.0;
            if (gk.getNonPolarModel() == NonPolar.BORN_SOLV || gk.getNonPolarModel() == NonPolar.BORN_CAV_DISP) {
                sTens = gk.getSurfaceTension();
                sTens *= OpenMM_KJPerKcal;
                sTens *= 100.0; // 100 square Angstroms per square nanometer.
            }

            PointerByReference doubleArray = OpenMM_DoubleArray_create(0);
            for (Atom atom : atoms) {
                int index = atom.getXyzIndex() - 1;
                double chargeUseFactor = 1.0;
                if (!atom.getUse() || !atom.getElectrostatics()) {
                    chargeUseFactor = 0.0;
                }
                double lambdaScale = lambdaElec;
                if (!atom.applyLambda()) {
                    lambdaScale = 1.0;
                }

                chargeUseFactor *= lambdaScale;
                MultipoleType multipoleType = atom.getMultipoleType();
                double charge = multipoleType.charge * chargeUseFactor;
                double surfaceTension = sTens * chargeUseFactor;

                double overlapScaleUseFactor = nea ? 1.0 : chargeUseFactor;
                double oScale = overlapScale[index] * overlapScaleUseFactor;
                double baseRadius = baseRadii[index];

                OpenMM_DoubleArray_append(doubleArray, charge);
                OpenMM_DoubleArray_append(doubleArray, OpenMM_NmPerAngstrom * baseRadius);
                OpenMM_DoubleArray_append(doubleArray, oScale);
                OpenMM_DoubleArray_append(doubleArray, surfaceTension);
                OpenMM_CustomGBForce_setParticleParameters(customGBForce, index, doubleArray);

                // Reset the double array for the next atom.
                OpenMM_DoubleArray_resize(doubleArray, 0);
            }
            OpenMM_DoubleArray_destroy(doubleArray);

            if (openMMContext.context != null) {
                OpenMM_CustomGBForce_updateParametersInContext(customGBForce, openMMContext.context);
            }
        }

        /**
         * Updates the Amoeba electrostatic multipolar force for changes in Use flags or Lambda.
         *
         * @param atoms Array of atoms to update.
         */
        private void updateAmoebaMultipoleForce(Atom[] atoms) {
            ParticleMeshEwald pme = getPmeNode();
            int[][] axisAtom = pme.getAxisAtoms();
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

            for (Atom atom : atoms) {
                int index = atom.getXyzIndex() - 1;
                MultipoleType multipoleType = atom.getMultipoleType();
                PolarizeType polarType = atom.getPolarizeType();

                double useFactor = 1.0;
                if (!atom.getUse() || !atom.getElectrostatics()) {
                    useFactor = 0.0;
                }

                double lambdaScale = lambdaElec;
                if (!atom.applyLambda()) {
                    lambdaScale = 1.0;
                }
                useFactor *= lambdaScale;

                // Define the frame definition.
                int axisType;
                switch (multipoleType.frameDefinition) {
                    default:
                    case NONE:
                        axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
                        break;
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
                    case THREEFOLD:
                        axisType = OpenMM_AmoebaMultipoleForce_ThreeFold;
                        break;
                }

                // Load local multipole coefficients.
                for (int j = 0; j < 3; j++) {
                    OpenMM_DoubleArray_set(dipoles, j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
                }
                int l = 0;
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        OpenMM_DoubleArray_set(quadrupoles, l++, multipoleType.quadrupole[j][k]
                                * quadrupoleConversion / 3.0 * useFactor);
                    }
                }

                int zaxis = 1;
                int xaxis = 1;
                int yaxis = 1;
                int[] refAtoms = axisAtom[index];
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

                // Set the multipole parameters.
                OpenMM_AmoebaMultipoleForce_setMultipoleParameters(amoebaMultipoleForce, index,
                        multipoleType.charge * useFactor, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
                        polarType.thole, polarType.pdamp * dampingFactorConversion,
                        polarType.polarizability * polarityConversion * polarScale * useFactor);
            }

            OpenMM_DoubleArray_destroy(dipoles);
            OpenMM_DoubleArray_destroy(quadrupoles);

            if (openMMContext.context != null) {
                OpenMM_AmoebaMultipoleForce_updateParametersInContext(amoebaMultipoleForce, openMMContext.context);
            }
        }

        /**
         * Updates the AMOEBA Generalized Kirkwood force for changes in Use flags or Lambda.
         *
         * @param atoms Array of atoms to update.
         */
        private void updateAmoebaGeneralizedKirkwoodForce(Atom[] atoms) {
            GeneralizedKirkwood gk = getGK();
            double[] overlapScale = gk.getOverlapScale();
            double[] baseRadii = gk.getBaseRadii();
            boolean nea = gk.getNativeEnvironmentApproximation();

            for (Atom atom : atoms) {
                int index = atom.getXyzIndex() - 1;
                double chargeUseFactor = 1.0;
                if (!atom.getUse() || !atom.getElectrostatics()) {
                    chargeUseFactor = 0.0;
                }

                double lambdaScale = lambdaElec;
                if (!atom.applyLambda()) {
                    lambdaScale = 1.0;
                }

                chargeUseFactor *= lambdaScale;
                double overlapScaleUseFactor = nea ? 1.0 : chargeUseFactor;

                MultipoleType multipoleType = atom.getMultipoleType();
                OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters(amoebaGeneralizedKirkwoodForce, index,
                        multipoleType.charge * chargeUseFactor, OpenMM_NmPerAngstrom * baseRadii[index],
                        overlapScale[index] * overlapScaleUseFactor);
            }

//        OpenMM Bug: Surface Area is not Updated by "updateParametersInContext"
//
//        NonPolar nonpolar = gk.getNonPolarModel();
//        switch (nonpolar) {
//            case BORN_SOLV:
//            case BORN_CAV_DISP:
//            default:
//                // Configure a Born Radii based surface area term.
//                double surfaceTension = gk.getSurfaceTension() * OpenMM_KJPerKcal
//                        * OpenMM_AngstromsPerNm * OpenMM_AngstromsPerNm * lambdaElec;
//                OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(amoebaGeneralizedKirkwoodForce, -surfaceTension);
//                break;
//            case CAV:
//            case CAV_DISP:
//            case HYDROPHOBIC_PMF:
//            case NONE:
//                // This NonPolar model does not use a Born Radii based surface area term.
//                break;
//        }

            if (openMMContext.context != null) {
                OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(
                        amoebaGeneralizedKirkwoodForce, openMMContext.context);
            }

        }

        /**
         * Updates the WCA force for changes in Use flags or Lambda.
         *
         * @param atoms Array of atoms to update.
         */
        private void updateWCAForce(Atom[] atoms) {
            VanDerWaals vdW = getVdwNode();
            VanDerWaalsForm vdwForm = vdW.getVDWForm();
            double radScale = 1.0;
            if (vdwForm.radiusSize == VanDerWaalsForm.RADIUS_SIZE.DIAMETER) {
                radScale = 0.5;
            }

            for (Atom atom : atoms) {
                int index = atom.getXyzIndex() - 1;
                double useFactor = 1.0;
                if (!atom.getUse()) {
                    useFactor = 0.0;
                }

                // Scale all implicit solvent terms with the square of electrostatics lambda
                // (so dUdisp / dL is 0 at lambdaElec = 0).
                double lambdaScale = lambdaElec * lambdaElec;
                if (!atom.applyLambda()) {
                    lambdaScale = 1.0;
                }
                useFactor *= lambdaScale;

                VDWType vdwType = atom.getVDWType();
                double radius = vdwType.radius;
                double eps = vdwType.wellDepth;
                OpenMM_AmoebaWcaDispersionForce_setParticleParameters(amoebaWcaDispersionForce, index,
                        OpenMM_NmPerAngstrom * radius * radScale, OpenMM_KJPerKcal * eps * useFactor);
            }

            if (openMMContext.context != null) {
                OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(
                        amoebaWcaDispersionForce, openMMContext.context);
            }
        }

        /**
         * Updates the Torsion force for application of lambda scaling.
         */
        private void updateTorsionForce() {

            // Only update parameters if torsions are being scaled by lambda.
            if (!torsionLambdaTerm) {
                return;
            }

            // Check if this system has torsions.
            Torsion[] torsions = getTorsions();
            if (torsions == null || torsions.length < 1) {
                return;
            }

            int index = 0;
            for (Torsion torsion : torsions) {
                TorsionType torsionType = torsion.torsionType;
                int nTerms = torsionType.phase.length;
                if (!torsion.applyLambda()) {
                    index += nTerms;
                    continue;
                }

                int a1 = torsion.getAtom(0).getXyzIndex() - 1;
                int a2 = torsion.getAtom(1).getXyzIndex() - 1;
                int a3 = torsion.getAtom(2).getXyzIndex() - 1;
                int a4 = torsion.getAtom(3).getXyzIndex() - 1;

                for (int j = 0; j < nTerms; j++) {
                    double forceConstant = OpenMM_KJPerKcal * torsion.units * torsionType.amplitude[j] * lambdaTorsion;
                    OpenMM_PeriodicTorsionForce_setTorsionParameters(amoebaTorsionForce, index++, a1, a2, a3, a4,
                            j + 1, torsionType.phase[j] * OpenMM_RadiansPerDegree, forceConstant);
                }
            }

            if (openMMContext.context != null) {
                OpenMM_PeriodicTorsionForce_updateParametersInContext(amoebaTorsionForce, openMMContext.context);
            }
        }

        /**
         * Updates the Improper Torsion force for application of lambda scaling.
         */
        private void updateImproperTorsionForce() {
            // Only update parameters if torsions are being scaled by lambda.
            if (!torsionLambdaTerm) {
                return;
            }

            ImproperTorsion[] impropers = getImproperTorsions();
            if (impropers == null || impropers.length < 1) {
                return;
            }

            int nImpropers = impropers.length;
            for (int i = 0; i < nImpropers; i++) {
                ImproperTorsion improperTorsion = impropers[i];
                if (!improperTorsion.applyLambda()) {
                    continue;
                }

                int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
                int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
                int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
                int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;

                ImproperTorsionType improperTorsionType = improperTorsion.improperType;
                double forceConstant = OpenMM_KJPerKcal * improperTorsion.units
                        * improperTorsion.scaleFactor * improperTorsionType.k * lambdaTorsion;
                OpenMM_PeriodicTorsionForce_setTorsionParameters(amoebaImproperTorsionForce, i, a1, a2, a3, a4,
                        improperTorsionType.periodicity, improperTorsionType.phase * OpenMM_RadiansPerDegree, forceConstant);
            }

            if (openMMContext.context != null) {
                OpenMM_PeriodicTorsionForce_updateParametersInContext(amoebaImproperTorsionForce, openMMContext.context);
            }
        }

        /**
         * <p>
         * setUpHydrogenConstraints.</p>
         *
         * @param system a {@link com.sun.jna.ptr.PointerByReference} object.
         */
        private void setUpHydrogenConstraints(PointerByReference system) {
            Bond[] bonds = getBonds();
            if (bonds == null || bonds.length < 1) {
                return;
            }

            logger.info(" Setting up constraints for bonds to hydrogen atoms.");
            for (Bond bond : bonds) {
                Atom atom1 = bond.getAtom(0);
                Atom atom2 = bond.getAtom(1);
                if (atom1.isHydrogen() || atom2.isHydrogen()) {
                    BondType bondType = bond.bondType;
                    int iAtom1 = atom1.getXyzIndex() - 1;
                    int iAtom2 = atom2.getXyzIndex() - 1;
                    OpenMM_System_addConstraint(system, iAtom1, iAtom2, bondType.distance * OpenMM_NmPerAngstrom);
                }
            }
        }

        /**
         * <p>
         * setUpBondConstraints.</p>
         * Constrains all bonds in the system to a fixed length.
         *
         * @param system a {@link com.sun.jna.ptr.PointerByReference} object.
         */
        private void setUpBondConstraints(PointerByReference system) {
            Bond[] bonds = getBonds();
            if (bonds == null || bonds.length < 1) {
                return;
            }

            logger.info(" Setting up constraints for all bonds.");
            for (Bond bond : bonds) {
                Atom atom1 = bond.getAtom(0);
                Atom atom2 = bond.getAtom(1);
                int iAtom1 = atom1.getXyzIndex() - 1;
                int iAtom2 = atom2.getXyzIndex() - 1;
                OpenMM_System_addConstraint(system, iAtom1, iAtom2, bond.bondType.distance * OpenMM_NmPerAngstrom);
            }
        }

        /**
         * <p>setUpHydrogenAngleConstraints.</p>
         *
         * @param system a {@link com.sun.jna.ptr.PointerByReference} object.
         */
        private void setUpHydrogenAngleConstraints(PointerByReference system) {

            Angle[] angles = getAngles();
            if (angles == null || angles.length < 1) {
                return;
            }

            logger.info(" Setting up hydrogen angle constraints");

            for (Angle angle : angles) {
                if (isHydrogenAngle(angle)) {
                    Atom atom1 = angle.getAtom(0);
                    Atom atom3 = angle.getAtom(2);

                    // Calculate a "false bond" length between atoms 1 and 3 to constrain the angle using the law of cosines.
                    Bond bond1 = angle.getBond(0);
                    double distance1 = bond1.bondType.distance;

                    Bond bond2 = angle.getBond(1);
                    double distance2 = bond2.bondType.distance;

                    // Equilibrium angle value in degrees.
                    double angleVal = angle.angleType.angle[angle.nh];

                    // Law of cosines.
                    double falseBondLength = sqrt(distance1 * distance1 + distance2 * distance2
                            - 2.0 * distance1 * distance2 * cos(toRadians(angleVal)));

                    int iAtom1 = atom1.getXyzIndex() - 1;
                    int iAtom3 = atom3.getXyzIndex() - 1;
                    OpenMM_System_addConstraint(system, iAtom1, iAtom3, falseBondLength * OpenMM_NmPerAngstrom);
                }
            }
        }

        /**
         * Check to see if an angle is a hydrogen angle. This method only returns true for hydrogen angles that are less
         * than 160 degrees.
         *
         * @param angle Angle to check.
         * @return boolean indicating whether or not an angle is a hydrogen angle that is less than 160 degrees.
         */
        private boolean isHydrogenAngle(Angle angle) {
            if (angle.containsHydrogen()) {
                // Equilibrium angle value in degrees.
                double angleVal = angle.angleType.angle[angle.nh];
                // Make sure angle is less than 160 degrees because greater than 160 degrees will not be constrained
                // well using the law of cosines.
                if (angleVal < 160.0) {
                    Atom atom1 = angle.getAtom(0);
                    Atom atom2 = angle.getAtom(1);
                    Atom atom3 = angle.getAtom(2);
                    // Setting constraints only on angles where atom1 or atom3 is a hydrogen while atom2 is not a hydrogen.
                    return atom1.isHydrogen() && atom3.isHydrogen() && !atom2.isHydrogen();
                }
            }
            return false;
        }

        private void setDefaultPeriodicBoxVectors() {
            Crystal crystal = getCrystal();
            if (!crystal.aperiodic()) {
                OpenMM_Vec3 a = new OpenMM_Vec3();
                OpenMM_Vec3 b = new OpenMM_Vec3();
                OpenMM_Vec3 c = new OpenMM_Vec3();
                double[][] Ai = crystal.Ai;
                a.x = Ai[0][0] * OpenMM_NmPerAngstrom;
                a.y = Ai[0][1] * OpenMM_NmPerAngstrom;
                a.z = Ai[0][2] * OpenMM_NmPerAngstrom;
                b.x = Ai[1][0] * OpenMM_NmPerAngstrom;
                b.y = Ai[1][1] * OpenMM_NmPerAngstrom;
                b.z = Ai[1][2] * OpenMM_NmPerAngstrom;
                c.x = Ai[2][0] * OpenMM_NmPerAngstrom;
                c.y = Ai[2][1] * OpenMM_NmPerAngstrom;
                c.z = Ai[2][2] * OpenMM_NmPerAngstrom;
                OpenMM_System_setDefaultPeriodicBoxVectors(system, a, b, c);
            }
        }

        private void getPeriodicBoxVectors(PointerByReference state) {
            Crystal crystal = getCrystal();
            if (!crystal.aperiodic()) {
                OpenMM_Vec3 a = new OpenMM_Vec3();
                OpenMM_Vec3 b = new OpenMM_Vec3();
                OpenMM_Vec3 c = new OpenMM_Vec3();
                OpenMM_State_getPeriodicBoxVectors(state, a, b, c);
                double[][] latticeVectors = new double[3][3];
                latticeVectors[0][0] = a.x * OpenMM_AngstromsPerNm;
                latticeVectors[0][1] = a.y * OpenMM_AngstromsPerNm;
                latticeVectors[0][2] = a.z * OpenMM_AngstromsPerNm;
                latticeVectors[1][0] = b.x * OpenMM_AngstromsPerNm;
                latticeVectors[1][1] = b.y * OpenMM_AngstromsPerNm;
                latticeVectors[1][2] = b.z * OpenMM_AngstromsPerNm;
                latticeVectors[2][0] = c.x * OpenMM_AngstromsPerNm;
                latticeVectors[2][1] = c.y * OpenMM_AngstromsPerNm;
                latticeVectors[2][2] = c.z * OpenMM_AngstromsPerNm;
                crystal.setCellVectors(latticeVectors);
                ForceFieldEnergyOpenMM.super.setCrystal(crystal);
            }
        }
    }

    /**
     * Destroys pointer references to Context, Integrator and System to free up
     * memory.
     */
    private void freeOpenMM() {
        openMMContext.destroyContext();
        openMMSystem.destroy();
    }

    /**
     * Add an Andersen thermostat to the system.
     *
     * @param targetTemp    Target temperature in Kelvins.
     * @param collisionFreq Collision frequency in 1/psec
     */
    private void addAndersenThermostat(double targetTemp, double collisionFreq) {
        openMMSystem.addAndersenThermostat(targetTemp, collisionFreq);
    }

    /**
     * Private method for internal use, so we don't have subclasses calling
     * super.energy, and this class delegating to the subclass's getGradient
     * method.
     *
     * @param forces Reference to forces returned by OpenMM.
     * @param g      Gradient array to fill.
     * @return Gradient array.
     */
    private double[] fillGradients(PointerByReference forces, double[] g) {
        int n = getNumberOfVariables();
        if (g == null || g.length < n) {
            g = new double[n];
        }
        int index = 0;
        int nAtoms = atoms.length;
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i];
            if (a.isActive()) {
                OpenMM_Vec3 posInNm = OpenMM_Vec3Array_get(forces, i);

                // Convert OpenMM Forces in KJ/Nm into an FFX gradient in Kcal/A.
                double gx = -posInNm.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double gy = -posInNm.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                double gz = -posInNm.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
                if (isNaN(gx) || isInfinite(gx) || isNaN(gy) || isInfinite(gy) || isNaN(gz) || isInfinite(gz)) {
                    StringBuilder sb = new StringBuilder(
                            format(" The gradient of atom %s is (%8.3f,%8.3f,%8.3f).", a.toString(), gx, gy, gz));
                    double[] vals = new double[3];
                    a.getVelocity(vals);
                    sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getAcceleration(vals);
                    sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
                    a.getPreviousAcceleration(vals);
                    sb.append(format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
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

    private void setDefaultPeriodicBoxVectors() {
        openMMSystem.setDefaultPeriodicBoxVectors();
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

        logger.info(format(" Number of CUDA devices: %d.", nDevs));

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
     * Create a JNA Pointer to a String.
     *
     * @param string WARNING: assumes ascii-only string
     * @return pointer.
     */
    private static Pointer pointerForString(String string) {
        Pointer pointer = new Memory(string.length() + 1);
        pointer.setString(0, string);
        return pointer;
    }

    /**
     * Returns the platform array as a String
     *
     * @param stringArray The OpenMM String array.
     * @param i           The index of the String to return.
     * @return String The requested String.
     */
    private static String stringFromArray(PointerByReference stringArray, int i) {
        Pointer platformPtr = OpenMM_StringArray_get(stringArray, i);
        if (platformPtr == null) {
            return null;
        }
        return platformPtr.getString(0);
    }
}
