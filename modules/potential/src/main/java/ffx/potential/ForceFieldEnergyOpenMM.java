// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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

import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_destroy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_3D_DoubleArray_set;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_NonbondedMethod.OpenMM_AmoebaGKCavitationForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGKCavitationForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent12;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent13;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent14;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_Covalent15;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_CovalentType.OpenMM_AmoebaMultipoleForce_PolarizationCovalent11;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_Bisector;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_NoAxisType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ThreeFold;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZBisect;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZOnly;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_MultipoleAxisTypes.OpenMM_AmoebaMultipoleForce_ZThenX;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_NonbondedMethod.OpenMM_AmoebaMultipoleForce_PME;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Direct;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Extrapolated;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_PolarizationType.OpenMM_AmoebaMultipoleForce_Mutual;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_addMultipole;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setAEwald;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCovalentMap;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMultipoleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPmeGridDimensions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_setPolarizationType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaMultipoleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_AlchemicalMethod.OpenMM_AmoebaVdwForce_Decouple;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_CutoffPeriodic;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_NonbondedMethod.OpenMM_AmoebaVdwForce_NoCutoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticleType;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addParticle_1;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_addTypePair;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setAlchemicalMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleExclusions;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcoreAlpha;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setSoftcorePower;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_setUseDispersionCorrection;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaVdwForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_create;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setAwater;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setDispoff;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setEpsh;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setEpso;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setRminh;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setRmino;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setShctd;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_setSlevy;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AmoebaWcaDispersionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_AngstromsPerNm;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KJPerKcal;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_KcalPerKJ;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_RadiansPerDegree;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_setDefaultCollisionFrequency;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_AndersenThermostat_setDefaultTemperature;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_BondArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_False;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Boolean.OpenMM_True;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CMMotionRemover_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_applyConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_create_2;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_getState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_reinitialize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Context_setVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addAngle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_addPerAngleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_setAngleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomAngleForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_addPerBondParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomCompoundBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomExternalForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePair;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_ParticlePairNoExclusions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_ComputationType.OpenMM_CustomGBForce_SingleParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addComputedValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addEnergyTerm;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomGBForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addComputePerDof;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addConstrainVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addPerDofVariable;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_addUpdateContextState;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomIntegrator_setKineticEnergyExpression;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addExclusion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addGlobalParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addInteractionGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_addPerParticleParameter;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_setUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_resize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_DoubleArray_set;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setForceGroup;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Force_setName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_addBond;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_setBondParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_HarmonicBondForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_resize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntArray_set;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntSet_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntSet_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_IntSet_insert;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_setConstraintTolerance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Integrator_step;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinIntegrator_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LangevinIntegrator_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_LocalEnergyMinimizer_minimize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_MonteCarloBarostat_setRandomNumberSeed;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_createExceptionsFromBonds;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getNumExceptions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_getParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setCutoffDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setExceptionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setNonbondedMethod;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setPMEParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setParticleParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setSwitchingDistance;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setUseDispersionCorrection;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_setUseSwitchingFunction;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_addTorsion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_setTorsionParameters;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_PeriodicTorsionForce_updateParametersInContext;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getNumPlatforms;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getOpenMMVersion;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPlatformByName;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_getPluginLoadFailures;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_loadPluginsFromDirectory;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Platform_setPropertyDefaultValue;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Energy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Forces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Positions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_DataType.OpenMM_State_Velocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getForces;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getKineticEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPositions;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getPotentialEnergy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_State_getVelocities;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_get;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_StringArray_getSize;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addConstraint;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addForce;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_addParticle;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_getNumConstraints;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setDefaultPeriodicBoxVectors;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_System_setParticleMass;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_append;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_create;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_destroy;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_Vec3Array_get;
import static edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_VerletIntegrator_create;
import static ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar.GAUSS_DISP;
import static ffx.potential.nonbonded.VanDerWaalsForm.EPSILON_RULE.GEOMETRIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_RULE.ARITHMETIC;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_SIZE.RADIUS;
import static ffx.potential.nonbonded.VanDerWaalsForm.RADIUS_TYPE.R_MIN;
import static ffx.potential.nonbonded.VanDerWaalsForm.VDW_TYPE.LENNARD_JONES;
import static ffx.utilities.Constants.KCAL_TO_GRAM_ANG2_PER_PS2;
import static ffx.utilities.Constants.kB;
import static java.lang.Double.isFinite;
import static java.lang.Double.isInfinite;
import static java.lang.Double.isNaN;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.PI;
import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.pow;
import static org.apache.commons.math3.util.FastMath.round;
import static org.apache.commons.math3.util.FastMath.sqrt;
import static org.apache.commons.math3.util.FastMath.toRadians;

import com.sun.jna.Memory;
import com.sun.jna.Pointer;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import edu.rit.mp.CharacterBuf;
import edu.rit.pj.Comm;
import edu.uiowa.jopenmm.OpenMMLibrary;
import edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_CustomNonbondedForce_NonbondedMethod;
import edu.uiowa.jopenmm.OpenMMLibrary.OpenMM_NonbondedForce_NonbondedMethod;
import edu.uiowa.jopenmm.OpenMMUtils;
import edu.uiowa.jopenmm.OpenMM_Vec3;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.AngleTorsion;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.ImproperTorsion;
import ffx.potential.bonded.OutOfPlaneBend;
import ffx.potential.bonded.PiOrbitalTorsion;
import ffx.potential.bonded.RestraintBond;
import ffx.potential.bonded.RestraintTorsion;
import ffx.potential.bonded.StretchBend;
import ffx.potential.bonded.StretchTorsion;
import ffx.potential.bonded.Torsion;
import ffx.potential.bonded.TorsionTorsion;
import ffx.potential.bonded.UreyBradley;
import ffx.potential.extended.ExtendedSystem;
import ffx.potential.nonbonded.CoordRestraint;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar;
import ffx.potential.nonbonded.NonbondedCutoff;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ParticleMeshEwald.SCFAlgorithm;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.RestrainGroups;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.VanDerWaalsForm;
import ffx.potential.nonbonded.implicit.ChandlerCavitation;
import ffx.potential.nonbonded.implicit.DispersionRegion;
import ffx.potential.nonbonded.implicit.GaussVol;
import ffx.potential.parameters.AngleType;
import ffx.potential.parameters.AngleType.AngleFunction;
import ffx.potential.parameters.AngleType.AngleMode;
import ffx.potential.parameters.BondType;
import ffx.potential.parameters.BondType.BondFunction;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ImproperTorsionType;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.OutOfPlaneBendType;
import ffx.potential.parameters.PiOrbitalTorsionType;
import ffx.potential.parameters.PolarizeType;
import ffx.potential.parameters.TorsionTorsionType;
import ffx.potential.parameters.TorsionType;
import ffx.potential.parameters.UreyBradleyType;
import ffx.potential.parameters.VDWPairType;
import ffx.potential.parameters.VDWType;
import ffx.potential.utils.EnergyException;
import ffx.potential.utils.PotentialsFunctions;
import ffx.potential.utils.PotentialsUtils;
import ffx.utilities.Constants;
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
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * Compute the potential energy and derivatives using OpenMM.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
@SuppressWarnings("deprecation")
public class ForceFieldEnergyOpenMM extends ForceFieldEnergy {

  private static final Logger logger = Logger.getLogger(ForceFieldEnergyOpenMM.class.getName());
  /** Whether to enforce periodic boundary conditions when obtaining new States. */
  public final int enforcePBC;
  /** The atoms this ForceFieldEnergyOpenMM operates on. */
  private final Atom[] atoms;
  /** Number of particles. */
  private final int nAtoms;
  /** Lambda step size for finite difference dU/dL. */
  private final double finiteDifferenceStepSize;
  /** OpenMM Context. */
  private Context context;
  /** OpenMM System. */
  private System system;
  /**
   * Truncate the normal OpenMM Lambda Path from 0..1 to Lambda_Start..1. This is useful for
   * conformational optimization if full removal of vdW interactions is not desired (i.e. lambdaStart
   * = ~0.2).
   */
  private double lambdaStart = 0.0;
  /** Use two-sided finite difference dU/dL. */
  private boolean twoSidedFiniteDifference = true;

  private final boolean freeOpenMM;

  /**
   * ForceFieldEnergyOpenMM constructor; offloads heavy-duty computation to an OpenMM Platform while
   * keeping track of information locally.
   *
   * @param molecularAssembly Assembly to construct energy for.
   * @param requestedPlatform requested OpenMM platform to be used.
   * @param restraints Harmonic coordinate restraints.
   * @param nThreads Number of threads to use in the super class ForceFieldEnergy instance.
   */
  protected ForceFieldEnergyOpenMM(
      MolecularAssembly molecularAssembly,
      Platform requestedPlatform,
      List<CoordRestraint> restraints,
      int nThreads) {
    super(molecularAssembly, restraints, nThreads);

    Crystal crystal = getCrystal();
    int symOps = crystal.spaceGroup.getNumberOfSymOps();
    if (symOps > 1) {
      logger.info("");
      logger.severe(" OpenMM does not support symmetry operators.");
    }

    logger.info("\n Initializing OpenMM");

    ForceField forceField = molecularAssembly.getForceField();
    context = new Context(forceField, requestedPlatform);

    atoms = molecularAssembly.getAtomArray();
    nAtoms = atoms.length;
    system = new System(molecularAssembly);

    boolean aperiodic = super.getCrystal().aperiodic();
    boolean pbcEnforced = forceField.getBoolean("ENFORCE_PBC", !aperiodic);
    enforcePBC = pbcEnforced ? OpenMM_True : OpenMM_False;

    finiteDifferenceStepSize = forceField.getDouble("FD_DLAMBDA", 0.001);
    twoSidedFiniteDifference = forceField.getBoolean("FD_TWO_SIDED", twoSidedFiniteDifference);

    freeOpenMM = forceField.getBoolean("FREE_OPENMM", true);
  }

  /**
   * Gets the default co-processor device, ignoring any CUDA_DEVICE over-ride. This is either
   * determined by process rank and the availableDevices/CUDA_DEVICES property, or just 0 if neither
   * property is sets.
   *
   * @param props Properties in use.
   * @return Pre-override device index.
   */
  private static int getDefaultDevice(CompositeConfiguration props) {
    String availDeviceProp = props.getString("availableDevices", props.getString("CUDA_DEVICES"));
    if (availDeviceProp == null) {
      int nDevs = props.getInt("numCudaDevices", 1);
      availDeviceProp =
          IntStream.range(0, nDevs).mapToObj(Integer::toString).collect(Collectors.joining(" "));
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
          logger.severe(
              String.format(
                  " Failure at the allGather step for determining rank: %s\n%s",
                  ex, Utilities.stackTraceToString(ex)));
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
          logger.severe(
              String.format(
                  " Rank %d: Could not find any incoming host messages matching self %s!",
                  rank, host.trim()));
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
   * @param i The index of the String to return.
   * @return String The requested String.
   */
  private static String stringFromArray(PointerByReference stringArray, int i) {
    Pointer platformPtr = OpenMM_StringArray_get(stringArray, i);
    if (platformPtr == null) {
      return null;
    }
    return platformPtr.getString(0);
  }

  /**
   * Create an OpenMM Context.
   *
   * <p>Context.free() must be called to free OpenMM memory.
   *
   * @param integratorString Integrator to use.
   * @param timeStep Time step.
   * @param temperature Temperature (K).
   * @param forceCreation Force a new Context to be created, even if the existing one matches the
   *     request.
   */
  public void createContext(
      String integratorString, double timeStep, double temperature, boolean forceCreation) {
    context.create(integratorString, timeStep, temperature, forceCreation);
  }

  /**
   * Create an immutable OpenMM State.
   *
   * <p>State.free() must be called to free OpenMM memory.
   *
   * @param positions Retrieve positions.
   * @param energies Retrieve energies.
   * @param forces Retrieve forces.
   * @param velocities Retrieve velocities.
   * @return Returns the State.
   */
  public State createState(
      boolean positions, boolean energies, boolean forces, boolean velocities) {
    return new State(positions, energies, forces, velocities);
  }

  /** {@inheritDoc} */
  @Override
  public boolean destroy() {
    boolean ffxFFEDestroy = super.destroy();
    if (freeOpenMM) {
      free();
      logger.fine(" Destroyed the Context, Integrator, and OpenMMSystem.");
    }
    return ffxFFEDestroy;
  }

  /** {@inheritDoc} */
  @Override
  public double energy(double[] x) {
    return energy(x, false);
  }

  /** {@inheritDoc} */
  @Override
  public double energy(double[] x, boolean verbose) {

    if (lambdaBondedTerms) {
      return 0.0;
    }

    // Make sure a context has been created.
    context.getContextPointer();

    updateParameters(atoms);

    // Unscale the coordinates.
    unscaleCoordinates(x);

    setCoordinates(x);

    State state = new State(false, true, false, false);
    double e = state.potentialEnergy;
    state.free();

    if (!isFinite(e)) {
      String message = String.format(" Energy from OpenMM was a non-finite %8g", e);
      logger.warning(message);
      if (lambdaTerm) {
        system.printLambdaValues();
      }
      throw new EnergyException(message);
    }

    if (verbose) {
      logger.log(Level.INFO, String.format("\n OpenMM Energy: %14.10g", e));
    }

    // Rescale the coordinates.
    scaleCoordinates(x);

    return e;
  }

  /** {@inheritDoc} */
  @Override
  public double energyAndGradient(double[] x, double[] g) {
    return energyAndGradient(x, g, false);
  }

  /** {@inheritDoc} */
  @Override
  public double energyAndGradient(double[] x, double[] g, boolean verbose) {
    if (lambdaBondedTerms) {
      return 0.0;
    }

    // ZE BUG: updateParameters only gets called for energy(), not energyAndGradient().

    // Un-scale the coordinates.
    unscaleCoordinates(x);

    // Make sure a context has been created.
    context.getContextPointer();

    setCoordinates(x);

    State state = new State(false, true, true, false);
    double e = state.potentialEnergy;
    g = state.getGradient(g);
    state.free();

    if (!isFinite(e)) {
      String message = format(" Energy from OpenMM was a non-finite %8g", e);
      logger.warning(message);
      if (lambdaTerm) {
        system.printLambdaValues();
      }
      throw new EnergyException(message);
    }

    // if (vdwLambdaTerm) {
    //    PointerByReference parameterArray = OpenMM_State_getEnergyParameterDerivatives(state);
    //    int numDerives = OpenMM_ParameterArray_getSize(parameterArray);
    //    if (numDerives > 0) {
    //        double vdwdUdL = OpenMM_ParameterArray_get(parameterArray,
    // pointerForString("vdw_lambda")) / OpenMM_KJPerKcal;
    //    }
    // }

    if (maxDebugGradient < Double.MAX_VALUE) {
      boolean extremeGrad =
          Arrays.stream(g)
              .anyMatch((double gi) -> (gi > maxDebugGradient || gi < -maxDebugGradient));
      if (extremeGrad) {
        File origFile = molecularAssembly.getFile();
        String timeString =
            LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyy_MM_dd-HH_mm_ss"));

        String filename =
            format(
                "%s-LARGEGRAD-%s.pdb",
                FilenameUtils.removeExtension(molecularAssembly.getFile().getName()), timeString);
        PotentialsFunctions ef = new PotentialsUtils();
        filename = ef.versionFile(filename);

        logger.warning(
            format(
                " Excessively large gradients detected; printing snapshot to file %s", filename));
        ef.saveAsPDB(molecularAssembly, new File(filename));
        molecularAssembly.setFile(origFile);
      }
    }

    if (verbose) {
      logger.log(Level.INFO, format("\n OpenMM Energy: %14.10g", e));
    }

    // Scale the coordinates and gradients.
    scaleCoordinatesAndGradient(x, g);

    return e;
  }

  /**
   * Compute the energy and gradient using the pure Java code path.
   *
   * @param x Input atomic coordinates
   * @param g Storage for the gradient vector.
   * @return The energy (kcal/mol)
   */
  public double energyAndGradientFFX(double[] x, double[] g) {
    return super.energyAndGradient(x, g, false);
  }

  /**
   * Compute the energy and gradient using the pure Java code path.
   *
   * @param x Input atomic coordinates
   * @param g Storage for the gradient vector.
   * @param verbose Use verbose logging.
   * @return The energy (kcal/mol)
   */
  public double energyAndGradientFFX(double[] x, double[] g, boolean verbose) {
    return super.energyAndGradient(x, g, verbose);
  }

  /**
   * Compute the energy using the pure Java code path.
   *
   * @param x Atomic coordinates.
   * @return The energy (kcal/mol)
   */
  public double energyFFX(double[] x) {
    return super.energy(x, false);
  }

  /**
   * Compute the energy using the pure Java code path.
   *
   * @param x Input atomic coordinates
   * @param verbose Use verbose logging.
   * @return The energy (kcal/mol)
   */
  public double energyFFX(double[] x, boolean verbose) {
    return super.energy(x, verbose);
  }

  /** {@inheritDoc} */
  @Override
  public void finalize() throws Throwable {
    // Safer to leave super.finalize() in, even though right now that calls Object.finalize().
    logger.info(" ForceFieldEnergyOpenMM instance is being finalized.");
    super.finalize();
    if (destroyed) {
      logger.info(
          String.format(
              " Finalize called on a destroyed OpenMM ForceFieldEnergy %s", this.toString()));
    } else {
      destroy();
    }
  }

  /**
   * Returns the Context instance.
   *
   * @return context
   */
  public Context getContext() {
    return context;
  }

  /**
   * Re-compute the gradient using OpenMM and return it.
   *
   * @param g Gradient array.
   */
  @Override
  public double[] getGradient(double[] g) {
    State state = new State(false, false, true, false);
    g = state.getGradient(g);
    state.free();
    return g;
  }

  /** {@inheritDoc} */
  @Override
  public Platform getPlatform() {
    return context.platform;
  }

  /**
   * Get a reference to the System instance.
   *
   * @return Java wrapper to an OpenMM system.
   */
  public System getSystem() {
    return system;
  }

  /** {@inheritDoc} */
  @Override
  public double getd2EdL2() {
    return 0.0;
  }

  /** {@inheritDoc} */
  @Override
  public double getdEdL() {
    // No lambda dependence.
    if (!lambdaTerm) {
      return 0.0;
    }

    // Small optimization to only create the x array once.
    double[] x = new double[getNumberOfVariables()];
    getCoordinates(x);

    double currentLambda = getLambda();
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

    // logger.info(format(" getdEdL currentLambda: CL=%8.6f L=%8.6f dEdL=%12.6f", currentLambda,
    // lambda, dEdL));
    return dEdL;
  }

  /** {@inheritDoc} */
  @Override
  public void getdEdXdL(double[] gradients) {
    // Note for ForceFieldEnergyOpenMM this method is not implemented.
  }

  /** Update active atoms. */
  public void setActiveAtoms() {
    system.updateAtomMass();
    // Tests show reinitialization of the OpenMM Context is not necessary to pick up mass changes.
    // context.reinitContext();
  }

  /**
   * Set FFX and OpenMM coordinates for active atoms.
   *
   * @param x Atomic coordinates.
   */
  @Override
  public void setCoordinates(double[] x) {
    // Set both OpenMM and FFX coordinates to x.
    context.setOpenMMPositions(x);
  }

  /** {@inheritDoc} */
  @Override
  public void setCrystal(Crystal crystal) {
    super.setCrystal(crystal);
    context.setPeriodicBoxVectors();
  }

  /** {@inheritDoc} */
  @Override
  public void setLambda(double lambda) {

    if (!lambdaTerm) {
      logger.fine(" Attempting to set lambda for a ForceFieldEnergyOpenMM with lambdaterm false.");
      return;
    }

    // Check for lambda outside the range [0 .. 1].
    if (lambda < 0.0 || lambda > 1.0) {
      String message = format(" Lambda value %8.3f is not in the range [0..1].", lambda);
      logger.warning(message);
      return;
    }

    super.setLambda(lambda);

    // Remove the beginning of the normal Lambda path.
    double mappedLambda = lambda;
    if (lambdaStart > 0) {
      double windowSize = 1.0 - lambdaStart;
      mappedLambda = lambdaStart + lambda * windowSize;
    }

    if (system != null) {
      system.setLambda(mappedLambda);
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
    if (atoms == null) {
      atoms = this.atoms;
    }
    system.updateParameters(atoms);
  }

  /** Free OpenMM memory for the Context, Integrator and System. */
  private void free() {
    if (context != null) {
      context.free();
      context = null;
    }
    if (system != null) {
      system.free();
      system = null;
    }
  }

  /**
   * Creates and manage an OpenMM Context.
   *
   * <p>A Context stores the complete state of a simulation. More specifically, it includes: The
   * current time The position of each particle The velocity of each particle The values of
   * configurable parameters defined by Force objects in the System
   *
   * <p>You can retrieve a snapshot of the current state at any time by calling getState(). This
   * allows you to record the state of the simulation at various points, either for analysis or for
   * checkpointing. getState() can also be used to retrieve the current forces on each particle and
   * the current energy of the System.
   */
  public class Context {

    /** Requested Platform (i.e. Java or an OpenMM platform). */
    private final Platform platform;
    /** Instance of the OpenMM Integrator class. */
    private final Integrator integrator;
    /** Constraint tolerance as a fraction of the constrained bond length. */
    private final double constraintTolerance = ForceFieldEnergy.DEFAULT_CONSTRAINT_TOLERANCE;
    /** OpenMM Context pointer. */
    private PointerByReference contextPointer = null;
    /** Integrator string (default = VERLET). */
    private String integratorString = "VERLET";
    /** Time step (default = 0.001 psec). */
    private double timeStep = 0.001;
    /** OpenMM Platform pointer. */
    private PointerByReference platformPointer = null;
    /** Temperature (default = 298.15). */
    private double temperature = 298.15;

    /**
     * Create an OpenMM Context.
     *
     * @param forceField ForceField to used to create an integrator.
     * @param requestedPlatform Platform requested.
     */
    Context(ForceField forceField, Platform requestedPlatform) {
      loadPlatform(requestedPlatform);
      platform = requestedPlatform;
      integrator = new Integrator(forceField, constraintTolerance);
    }

    /**
     * Get a Pointer to the OpenMM Context. A Context is created if none has been instantiated yet.
     *
     * @return Context pointer.
     */
    public PointerByReference getContextPointer() {
      if (contextPointer == null) {
        create(integratorString, timeStep, temperature, true);
      }
      return contextPointer;
    }

    /**
     * Get a Pointer to the OpenMM Integrator.
     *
     * @return Integrator pointer.
     */
    public PointerByReference getIntegrator() {
      return integrator.getIntegratorPointer();
    }

    /**
     * Use the Context / Integrator combination to take the requested number of steps.
     *
     * @param numSteps Number of steps to take.
     */
    public void integrate(int numSteps) {
      OpenMM_Integrator_step(integrator.getIntegratorPointer(), numSteps);
    }

    /**
     * Use the Context to optimize the system to the requested tolerance.
     *
     * @param eps Convergence criteria (kcal/mole/A).
     * @param maxIterations Maximum number of iterations.
     */
    public void optimize(double eps, int maxIterations) {
      OpenMM_LocalEnergyMinimizer_minimize(
          contextPointer, eps / (OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ), maxIterations);
    }

    /**
     * The array x contains atomic coordinates only for active atoms.
     *
     * @param x Atomic coordinate array for only active atoms.
     */
    public void setOpenMMPositions(double[] x) {
      PointerByReference positions = OpenMM_Vec3Array_create(0);
      OpenMM_Vec3.ByValue coords = new OpenMM_Vec3.ByValue();
      double[] d = new double[3];
      int index = 0;
      for (Atom a : atoms) {
        if (a.isActive()) {
          a.moveTo(x[index++], x[index++], x[index++]);
          a.getXYZ(d);
          coords.x = d[0] * OpenMM_NmPerAngstrom;
          coords.y = d[1] * OpenMM_NmPerAngstrom;
          coords.z = d[2] * OpenMM_NmPerAngstrom;
          OpenMM_Vec3Array_append(positions, coords);
        } else {
          // OpenMM requires coordinates for even "inactive" atoms with mass of zero.
          coords.x = a.getX() * OpenMM_NmPerAngstrom;
          coords.y = a.getY() * OpenMM_NmPerAngstrom;
          coords.z = a.getZ() * OpenMM_NmPerAngstrom;
          OpenMM_Vec3Array_append(positions, coords);
        }
      }
      OpenMM_Context_setPositions(contextPointer, positions);
      if (freeOpenMM) {
        logger.finer(" Free OpenMM positions.");
        OpenMM_Vec3Array_destroy(positions);
        logger.finer(" Free OpenMM positions completed.");
      }
    }

    /**
     * The array v contains velocity values for active atomic coordinates.
     *
     * @param v Velocity array for active atoms.
     */
    public void setOpenMMVelocities(double[] v) {
      PointerByReference velocities = OpenMM_Vec3Array_create(0);
      OpenMM_Vec3.ByValue vel = new OpenMM_Vec3.ByValue();
      int index = 0;
      double[] velocity = new double[3];
      for (Atom a : atoms) {
        if (a.isActive()) {
          a.setVelocity(v[index++], v[index++], v[index++]);
          a.getVelocity(velocity);
          vel.x = velocity[0] * OpenMM_NmPerAngstrom;
          vel.y = velocity[1] * OpenMM_NmPerAngstrom;
          vel.z = velocity[2] * OpenMM_NmPerAngstrom;
          OpenMM_Vec3Array_append(velocities, vel);
        } else {
          // OpenMM requires velocities for even "inactive" atoms with mass of zero.
          a.setVelocity(0.0, 0.0, 0.0);
          vel.x = 0.0;
          vel.y = 0.0;
          vel.z = 0.0;
          OpenMM_Vec3Array_append(velocities, vel);
        }
      }
      OpenMM_Context_setVelocities(contextPointer, velocities);

      if (freeOpenMM) {
        logger.finer(" Free OpenMM velocities.");
        OpenMM_Vec3Array_destroy(velocities);
        logger.finer(" Free OpenMM velocities completed.");
      }
    }

    /** Set the periodic box vectors for a context based on the crystal instance. */
    public void setPeriodicBoxVectors() {
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
        OpenMM_Context_setPeriodicBoxVectors(contextPointer, a, b, c);
      }
    }

    /** {@inheritDoc} */
    @Override
    public String toString() {
      return String.format(
          " OpenMM context with integrator %s, timestep %9.3g fsec, temperature %9.3g K, constraintTolerance %9.3g",
          integratorString, timeStep, temperature, constraintTolerance);
    }

    /** Free OpenMM memory for the current Context and Integrator. */
    void free() {
      if (integrator != null) {
        integrator.free();
      }
      if (contextPointer != null && freeOpenMM) {
        logger.fine(" Free OpenMM Context.");
        OpenMM_Context_destroy(contextPointer);
        logger.fine(" Free OpenMM Context completed.");
        contextPointer = null;
      }
    }

    /**
     * Construct a new Context in which to run a simulation.
     *
     * @param integratorString Requested integrator.
     * @param timeStep Time step (psec).
     * @param temperature Temperature (K).
     * @param forceCreation Force creation of a new context, even if the current one matches.
     * @return Pointer to the created OpenMM context.
     */
    Context create(
        String integratorString, double timeStep, double temperature, boolean forceCreation) {
      // Check if the current context is consistent with the requested context.
      if (contextPointer != null && !forceCreation) {
        if (this.temperature == temperature
            && this.timeStep == timeStep
            && this.integratorString.equalsIgnoreCase(integratorString)) {
          // All requested features agree.
          return this;
        }
      }

      this.integratorString = integratorString;
      this.timeStep = timeStep;
      this.temperature = temperature;

      if (contextPointer != null) {
        logger.fine(" Free OpenMM Context.");
        OpenMM_Context_destroy(contextPointer);
        logger.fine(" Free OpenMM Context completed.");
        contextPointer = null;
      }

      logger.info("\n Creating OpenMM Context");

      PointerByReference integratorPointer =
          integrator.createIntegrator(integratorString, this.timeStep, temperature);

      // Set lambda to 1.0 when creating a context to avoid OpenMM compiling out any terms.
      double currentLambda = getLambda();

      if (lambdaTerm) {
        ForceFieldEnergyOpenMM.this.setLambda(1.0);
      }

      // Create a context.
      contextPointer = OpenMM_Context_create_2(system.getSystem(), integratorPointer,
          platformPointer);

      // Revert to the current lambda value.
      if (lambdaTerm) {
        ForceFieldEnergyOpenMM.this.setLambda(currentLambda);
      }

      // Get initial positions and velocities for active atoms.
      int nVar = ForceFieldEnergyOpenMM.super.getNumberOfVariables();
      double[] x = new double[nVar];
      double[] v = new double[nVar];
      double[] vel3 = new double[3];
      int index = 0;
      for (Atom a : atoms) {
        if (a.isActive()) {
          a.getVelocity(vel3);
          // X-axis
          x[index] = a.getX();
          v[index++] = vel3[0];
          // Y-axis
          x[index] = a.getY();
          v[index++] = vel3[1];
          // Z-axis
          x[index] = a.getZ();
          v[index++] = vel3[2];
        }
      }

      // Load the current periodic box vectors.
      setPeriodicBoxVectors();

      // Load current atomic positions.
      setOpenMMPositions(x);

      // Load current velocities.
      setOpenMMVelocities(v);

      // Apply constraints starting from current atomic positions.
      OpenMM_Context_applyConstraints(contextPointer, constraintTolerance);

      // Application of constraints can change coordinates and velocities.
      // Retrieve them for consistency.
      State state = new State(true, false, false, true);
      state.getPositions(x);
      state.getVelocities(v);
      state.free();

      return this;
    }

    /**
     * Reinitialize the context.
     *
     * <p>When a Context is created, it may cache information about the System being simulated and
     * the Force objects contained in it. This means that, if the System or Forces are then modified,
     * the Context might not see all of the changes. Call reinitialize() to force the Context to
     * rebuild its internal representation of the System and pick up any changes that have been
     * made.
     *
     * <p>This is an expensive operation, so you should try to avoid calling it too frequently.
     */
    void reinitContext() {
      if (contextPointer != null) {
        int preserveState = 1;
        OpenMM_Context_reinitialize(contextPointer, preserveState);
      }
    }

    /**
     * Set the value of an adjustable parameter defined by a Force object in the System.
     *
     * @param name the name of the parameter to set.
     * @param value the value of the parameter.
     */
    void setParameter(String name, double value) {
      if (contextPointer != null) {
        OpenMM_Context_setParameter(contextPointer, name, value);
      }
    }

    /** Load an OpenMM Platform */
    private void loadPlatform(Platform requestedPlatform) {

      OpenMMUtils.init();

      logger.log(Level.INFO, " Loaded from:\n {0}", OpenMMLibrary.JNA_NATIVE_LIB.toString());

      // Print out the OpenMM Version.
      Pointer version = OpenMM_Platform_getOpenMMVersion();
      logger.log(Level.INFO, " Version: {0}", version.getString(0));

      // Print out the OpenMM lib directory.
      logger.log(Level.FINE, " Lib Directory:       {0}", OpenMMUtils.getLibDirectory());
      // Load platforms and print out their names.
      PointerByReference libs =
          OpenMM_Platform_loadPluginsFromDirectory(OpenMMUtils.getLibDirectory());
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
      PointerByReference plugins =
          OpenMM_Platform_loadPluginsFromDirectory(OpenMMUtils.getPluginDirectory());
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
      String precision =
          molecularAssembly.getForceField().getString("PRECISION", defaultPrecision).toLowerCase();
      precision = precision.replace("-precision", "");
      switch (precision) {
        case "double":
        case "mixed":
        case "single":
          logger.info(String.format(" Precision level: %s", precision));
          break;
        default:
          logger.info(
              String.format(
                  " Could not interpret precision level %s, defaulting to %s",
                  precision, defaultPrecision));
          precision = defaultPrecision;
          break;
      }

      if (cuda && requestedPlatform != Platform.OMM_REF) {
        int defaultDevice = getDefaultDevice(molecularAssembly.getProperties());
        platformPointer = OpenMM_Platform_getPlatformByName("CUDA");
        int deviceID = molecularAssembly.getForceField().getInteger("CUDA_DEVICE", defaultDevice);
        String deviceIDString = Integer.toString(deviceID);

        OpenMM_Platform_setPropertyDefaultValue(
            platformPointer, pointerForString("CudaDeviceIndex"), pointerForString(deviceIDString));
        OpenMM_Platform_setPropertyDefaultValue(
            platformPointer, pointerForString("Precision"), pointerForString(precision));
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
        platformPointer = OpenMM_Platform_getPlatformByName("Reference");
        logger.info(" Platform: AMOEBA CPU Reference");
      }
    }
  }

  /**
   * Create and manage an OpenMM Integrator.
   *
   * <p>An Integrator defines a method for simulating a System by integrating the equations of
   * motion.
   *
   * <p>Each Integrator object is bound to a particular Context which it integrates. This connection
   * is specified by passing the Integrator as an argument to the constructor of the Context.
   */
  private class Integrator {

    /** Constraint tolerance as a fraction of the constrained bond length. */
    private final double constraintTolerance;
    /** Langevin friction coefficient. */
    private final double frictionCoeff;
    /** OpenMM Integrator pointer. */
    private PointerByReference integratorPointer = null;

    /**
     * Create an Integrator instance.
     *
     * @param forceField the ForceField instance containing integrator parameters.
     * @param constraintTolerance The integrator constraint tolerance.
     */
    Integrator(ForceField forceField, double constraintTolerance) {
      this.constraintTolerance = constraintTolerance;
      frictionCoeff = forceField.getDouble("FRICTION_COEFF", 91.0);
    }

    /**
     * Return a reference to the integrator.
     *
     * @return Integrator reference.
     */
    public PointerByReference getIntegratorPointer() {
      return integratorPointer;
    }

    /**
     * Create a integrator.
     *
     * @param integratorString Name of the integrator to use.
     * @param timeStep Time step (psec).
     * @param temperature Target temperature (kelvin).
     * @return Integrator reference.
     */
    PointerByReference createIntegrator(
        String integratorString, double timeStep, double temperature) {
      switch (integratorString) {
        case "LANGEVIN":
          createLangevinIntegrator(temperature, frictionCoeff, timeStep);
          break;
        case "RESPA":
          // Read in the inner time step in psec.
          int in = molecularAssembly.getProperties().getInt("respa-dt", 4);
          if (in < 2) {
            in = 2;
          }
          double inner = timeStep / in;
          createRESPAIntegrator(inner, timeStep);
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
          createVerletIntegrator(timeStep);
      }

      return integratorPointer;
    }

    /**
     * Create a RESPA integrator.
     *
     * @param inner Inner time step (psec).
     * @param dt Outer time step (psec).
     */
    private void createRESPAIntegrator(double inner, double dt) {
      createCustomIntegrator(dt);
      OpenMM_CustomIntegrator_addUpdateContextState(integratorPointer);
      OpenMM_CustomIntegrator_setKineticEnergyExpression(integratorPointer, "m*v*v/2");

      int n = (int) (round(dt / inner));
      StringBuilder e1 = new StringBuilder("v+0.5*(dt/" + n + ")*f0/m");
      StringBuilder e11 = new StringBuilder(n + "*(x-x1)/dt+" + e1);
      StringBuilder e2 = new StringBuilder("x+(dt/" + n + ")*v");

      OpenMM_CustomIntegrator_addPerDofVariable(integratorPointer, "x1", 0.0);
      OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", "v+0.5*dt*f1/m");
      for (int i = 0; i < n; i++) {
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", e1.toString());
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "x", e2.toString());
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "x1", "x");
        OpenMM_CustomIntegrator_addConstrainPositions(integratorPointer);
        OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", e11.toString());
        OpenMM_CustomIntegrator_addConstrainVelocities(integratorPointer);
      }
      OpenMM_CustomIntegrator_addComputePerDof(integratorPointer, "v", "v+0.5*dt*f1/m");
      OpenMM_CustomIntegrator_addConstrainVelocities(integratorPointer);
      logger.info("  Custom RESPA Integrator");
      logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
      logger.info(format("  Inner Time step:      %6.2f (fsec)", inner * 1000));
      logger.info(format("  Degrees of Freedom:   %6d", system.calculateDegreesOfFreedom()));
    }

    /**
     * Create a Langevin integrator.
     *
     * @param temperature Temperature (K).
     * @param frictionCoeff Frictional coefficient.
     * @param dt Time step (psec).
     */
    private void createLangevinIntegrator(double temperature, double frictionCoeff, double dt) {
      free();
      integratorPointer = OpenMM_LangevinIntegrator_create(temperature, frictionCoeff, dt);
      CompositeConfiguration properties = molecularAssembly.getProperties();
      if (properties.containsKey("integrator-seed")) {
        int randomSeed = properties.getInt("integrator-seed", 0);
        OpenMM_LangevinIntegrator_setRandomNumberSeed(integratorPointer, randomSeed);
      }

      OpenMM_Integrator_setConstraintTolerance(integratorPointer, constraintTolerance);
      logger.info("  Langevin Integrator");
      logger.info(format("  Target Temperature:   %6.2f (K)", temperature));
      logger.info(format("  Friction Coefficient: %6.2f", frictionCoeff));
      logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
      logger.info(format("  Degrees of Freedom:   %6d", system.calculateDegreesOfFreedom()));
    }

    /**
     * Create a Verlet integrator.
     *
     * @param dt Time step (psec).
     */
    private void createVerletIntegrator(double dt) {
      free();
      integratorPointer = OpenMM_VerletIntegrator_create(dt);
      OpenMM_Integrator_setConstraintTolerance(integratorPointer, constraintTolerance);
      logger.info("  Verlet Integrator");
      logger.info(format("  Time step:            %6.2f (fsec)", dt * 1000));
      logger.info(format("  Degrees of Freedom:   %6d", system.calculateDegreesOfFreedom()));
    }

    /**
     * Create a Custom integrator.
     *
     * @param dt Time step (psec).
     */
    private void createCustomIntegrator(double dt) {
      free();
      integratorPointer = OpenMM_CustomIntegrator_create(dt);
      OpenMM_Integrator_setConstraintTolerance(integratorPointer, constraintTolerance);
    }

    /** Destroy the integrator instance. */
    private void free() {
      if (integratorPointer != null && freeOpenMM) {
        logger.fine(" Free OpenMM Integrator.");
        OpenMM_Integrator_destroy(integratorPointer);
        logger.fine(" Free OpenMM Integrator completed.");
        integratorPointer = null;
      }
    }
  }

  /**
   * Create and manage an OpenMM System.
   *
   * <p>The definition of a System involves four elements:
   *
   * <p>The particles and constraints are defined directly by the System object, while forces are
   * defined by objects that extend the Force class. After creating a System, call addParticle() once
   * for each particle, addConstraint() for each constraint, and addForce() for each Force.
   *
   * <p>In addition, particles may be designated as "virtual sites". These are particles whose
   * positions are computed automatically based on the positions of other particles. To define a
   * virtual site, call setVirtualSite(), passing in a VirtualSite object that defines the rules for
   * computing its position.
   */
  public class System {

    private static final double DEFAULT_MELD_SCALE_FACTOR = -1.0;
    private final double meldScaleFactor;
    /** Andersen thermostat collision frequency. */
    private final double collisionFreq;
    /**
     * When using MELD, our goal will be to scale down the potential by this factor. A negative value
     * indicates we're not using MELD.
     */
    private final boolean useMeld;
    /** The Force Field in use. */
    ForceField forceField;
    /** Array of atoms in the sytem. */
    Atom[] atoms;
    /** OpenMM System. */
    private PointerByReference system;
    /** Barostat to be added if NPT (isothermal-isobaric) dynamics is requested. */
    private PointerByReference ommBarostat = null;
    /** OpenMM center-of-mass motion remover. */
    private PointerByReference commRemover = null;
    /**
     * This flag indicates bonded force constants and equilibria are updated (e.g. during ManyBody
     * titration).
     */
    private boolean updateBondedTerms = false;
    /**
     * If true, all torsions are treated as 6-fold, and all angles are treated as possibly changing
     * between normal and in-plane types.
     */
    private final boolean manyBodyTitration;
    /** OpenMM Custom Bond Force */
    private PointerByReference bondForce = null;
    /** OpenMM Custom Angle Force */
    private PointerByReference angleForce = null;
    /** OpenMM Custom Stretch-Bend Force */
    private PointerByReference stretchBendForce = null;
    /** OpenMM Custom In-Plane Angle Force */
    private PointerByReference inPlaneAngleForce = null;
    /** OpenMM Custom Urey-Bradley Force */
    private PointerByReference ureyBradleyForce = null;
    /** OpenMM Custom Out-of-Plane Bend Force */
    private PointerByReference outOfPlaneBendForce = null;
    /** OpenMM Custom Pi-Torsion Force */
    private PointerByReference piTorsionForce = null;
    /** OpenMM AMOEBA Torsion Force. */
    private PointerByReference torsionForce = null;
    private PointerByReference[] restraintTorsions = null;
    /** OpenMM Improper Torsion Force. */
    private PointerByReference improperTorsionForce = null;
    /** OpenMM AMOEBA van der Waals Force. */
    private PointerByReference amoebaVDWForce = null;
    /** OpenMM AMOEBA Multipole Force. */
    private PointerByReference amoebaMultipoleForce = null;
    /** OpenMM Generalized Kirkwood Force. */
    private PointerByReference amoebaGeneralizedKirkwoodForce = null;
    /** OpenMM AMOEBA WCA Dispersion Force. */
    private PointerByReference amoebaWcaDispersionForce = null;
    /** OpenMM AMOEBA WCA Cavitation Force. */
    private PointerByReference amoebaCavitationForce = null;
    /** OpenMM Custom GB Force. */
    private PointerByReference customGBForce = null;
    /** OpenMM Fixed Charge Non-Bonded Force. */
    private PointerByReference fixedChargeNonBondedForce = null;
    /** Fixed charge softcore vdW force boolean. */
    private boolean softcoreCreated = false;
    /** Boolean array, holds charge exclusion list. */
    private boolean[] chargeExclusion;
    /** Boolean array, holds van Der Waals exclusion list. */
    private boolean[] vdWExclusion;
    /** Double array, holds charge quantity value for exceptions. */
    private double[] exceptionChargeProd;
    /** Double array, holds epsilon quantity value for exceptions. */
    private double[] exceptionEps;
    /**
     * A map from vdW class values to OpenMM vdW types.
     */
    private Map<Integer, Integer> vdwClassToOpenMMType;
    /**
     * A class for a special vdW type that specifies zero energy (eps = 0.0; sigma = 1.0) for use
     * with the FFX "use" flag (e.g. use = false should give zero vdW energy for a many-body
     * expansion).
     */
    private int vdWClassForNoInteraction;
    /**
     * Lambda flag to indicate control of electrostatic scaling. If both elec and vdW are being
     * scaled, then vdW is scaled first, followed by elec.
     */
    private boolean elecLambdaTerm;
    /**
     * Lambda flag to indicate control of vdW scaling. If both elec and vdW are being scaled, then
     * vdW is scaled first, followed by elec.
     */
    private boolean vdwLambdaTerm;
    /**
     * Lambda flag to indicate control of torsional force constants (L=0 corresponds to torsions
     * being off, and L=1 to torsions at full strength.
     */
    private boolean torsionLambdaTerm;
    /** Value of the van der Waals lambda state variable. */
    private double lambdaVDW = 1.0;
    /** Value of the electrostatics lambda state variable. */
    private double lambdaElec = 1.0;
    /** Value of the electrostatics lambda state variable. */
    private double lambdaTorsion = 1.0;
    /**
     * The lambda value that defines when the electrostatics will start to turn on for full path
     * non-bonded term scaling.
     *
     * <p>A value of 0.6 works well for Chloride ion solvation, which is a difficult case due to the
     * ion having a formal negative charge and a large polarizability.
     */
    private double electrostaticStart = 0.6;
    /** Electrostatics lambda is raised to this power. */
    private double electrostaticLambdaPower;
    /** van der Waals softcore alpha. */
    private double vdWSoftcoreAlpha = 0.25;
    /** OpenMM thermostat. Currently an Andersen thermostat is supported. */
    private PointerByReference ommThermostat = null;
    /** van der Waals softcore beta. */
    private double vdwSoftcorePower = 3.0;
    /** Torsional lambda power. */
    private double torsionalLambdaPower = 2.0;

    /**
     * OpenMMSystem constructor.
     *
     * @param molecularAssembly MolecularAssembly used to construct the OpenMM system.
     */
    System(MolecularAssembly molecularAssembly) {
      // Create the OpenMM System
      system = OpenMM_System_create();
      logger.info("\n System created");

      forceField = molecularAssembly.getForceField();
      atoms = molecularAssembly.getAtomArray();

      // Load atoms.
      try {
        addAtoms();
      } catch (Exception e) {
        logger.severe(" Atom without mass encountered.");
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

      manyBodyTitration = forceField.getBoolean("MANYBODY_TITRATION", false);

      if (!forceField.getBoolean("LAMBDATERM", false)) {
        lambdaTerm = (elecLambdaTerm || vdwLambdaTerm || torsionLambdaTerm);
      } else {
        lambdaTerm = true;
      }

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

      electrostaticStart = forceField.getDouble("PERMANENT_LAMBDA_START", electrostaticStart);
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

      collisionFreq = forceField.getDouble("COLLISION_FREQ", 0.1);

      // Set up rigid constraints. These flags need to be set before bonds and angles are created
      // below.
      boolean rigidHydrogen = forceField.getBoolean("RIGID_HYDROGEN", false);
      boolean rigidBonds = forceField.getBoolean("RIGID_BONDS", false);
      boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
      if (rigidHydrogen) {
        addHydrogenConstraints();
      }
      if (rigidBonds) {
        addUpBondConstraints();
      }
      if (rigidHydrogenAngles) {
        setUpHydrogenAngleConstraints();
      }

      logger.info("\n Bonded Terms\n");

      // Add Bond Force.
      if (rigidBonds) {
        logger.info(" Not creating AmoebaBondForce because bonds are constrained.");
      } else {
        addBondForce();
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
      addRestraintBondForce();

      // Add Restrain Groups.
      addRestrainGroupsForce();

      // Add stretch-torsion coupling terms.
      addStretchTorsionForce();

      // Add angle-torsion coupling terms.
      addAngleTorsionForce();

      setDefaultPeriodicBoxVectors();

      addRestraintTorsions();

      if (vdW != null) {
        logger.info("\n Non-Bonded Terms\n");
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

    /**
     * Add an Andersen thermostat to the system.
     *
     * @param targetTemp Target temperature in Kelvins.
     */
    public void addAndersenThermostatForce(double targetTemp) {
      addAndersenThermostatForce(targetTemp, collisionFreq);
    }

    /**
     * Add an Andersen thermostat to the system.
     *
     * @param targetTemp Target temperature in Kelvins.
     * @param collisionFreq Collision frequency in 1/psec.
     */
    public void addAndersenThermostatForce(double targetTemp, double collisionFreq) {
      if (ommThermostat == null) {
        ommThermostat = OpenMM_AndersenThermostat_create(targetTemp, collisionFreq);
        OpenMM_System_addForce(system, ommThermostat);
        logger.info("\n Adding an Andersen thermostat");
        logger.info(format("  Target Temperature:   %6.2f (K)", targetTemp));
        logger.info(format("  Collision Frequency:  %6.2f (1/psec)", collisionFreq));
      } else {
        OpenMM_AndersenThermostat_setDefaultTemperature(ommThermostat, targetTemp);
        OpenMM_AndersenThermostat_setDefaultCollisionFrequency(ommThermostat, collisionFreq);
        logger.fine(" Updated the Andersen thermostat");
        logger.fine(format("  Target Temperature:   %6.2f (K)", targetTemp));
        logger.fine(format("  Collision Frequency:  %6.2f (1/psec)", collisionFreq));
      }
    }

    /** Adds a force that removes center-of-mass motion. */
    public void addCOMMRemoverForce() {
      int frequency = 100;
      if (commRemover == null) {
        commRemover = OpenMM_CMMotionRemover_create(frequency);
        OpenMM_System_addForce(system, commRemover);
        logger.info("\n Adding a center of mass motion remover");
        logger.info(format("  Frequency:            %6d", frequency));
      }
    }

    /**
     * Add a Monte Carlo Barostat to the system.
     *
     * @param targetPressure The target pressure (in atm).
     * @param targetTemp The target temperature.
     * @param frequency The frequency to apply the barostat.
     */
    public void addMonteCarloBarostatForce(
        double targetPressure, double targetTemp, int frequency) {
      if (ommBarostat == null) {
        double pressureInBar = targetPressure * Constants.ATM_TO_BAR;
        ommBarostat = OpenMM_MonteCarloBarostat_create(pressureInBar, targetTemp, frequency);
        CompositeConfiguration properties = molecularAssembly.getProperties();
        if (properties.containsKey("barostat-seed")) {
          int randomSeed = properties.getInt("barostat-seed", 0);
          logger.info(format(" Setting random seed %d for Monte Carlo Barostat", randomSeed));
          OpenMM_MonteCarloBarostat_setRandomNumberSeed(ommBarostat, randomSeed);
        }
        OpenMM_System_addForce(system, ommBarostat);
        logger.info("\n Adding a Monte Carlo barostat");
        logger.info(format("  Target Pressure:      %6.2f (atm)", targetPressure));
        logger.info(format("  Target Temperature:   %6.2f (K)", targetTemp));
        logger.info(format("  MC Move Frequency:    %6d", frequency));
      } else {
        logger.fine("\n Updating the Monte Carlo barostat");
        logger.fine(format("  Target Pressure:      %6.2f (atm)", targetPressure));
        logger.fine(format("  Target Temperature:   %6.2f (K)", targetTemp));
        logger.fine(format("  MC Move Frequency:    %6d", frequency));
      }
    }

    /**
     * Calculate the number of degrees of freedom.
     *
     * @return Number of degrees of freedom.
     */
    public int calculateDegreesOfFreedom() {
      // Begin from the 3 times the number of active atoms.
      int dof = getNumberOfVariables();
      // Remove OpenMM constraints.
      dof = dof - OpenMM_System_getNumConstraints(system);
      // Remove center of mass motion.
      if (commRemover != null) {
        dof -= 3;
      }
      return dof;
    }

    /** Destroy the system. */
    public void free() {
      if (system != null && freeOpenMM) {
        logger.fine(" Free OpenMM system.");
        OpenMM_System_destroy(system);
        logger.fine(" Free OpenMM system completed.");
        system = null;
      }
    }

    /** Print current lambda values. */
    public void printLambdaValues() {
      logger.info(
          format(
              "\n Lambda Values\n Torsion: %6.3f vdW: %6.3f Elec: %6.3f ",
              lambdaTorsion, lambdaVDW, lambdaElec));
    }

    /**
     * Set the overall lambda value for the system.
     *
     * @param lambda Current lambda value.
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

    public void setUpdateBondedTerms(boolean updateBondedTerms) {
      this.updateBondedTerms = updateBondedTerms;
    }

    /**
     * Return a reference to the System.
     *
     * @return System referenece.
     */
    PointerByReference getSystem() {
      return system;
    }

    /**
     * Set the default values of the vectors defining the axes of the periodic box (measured in nm).
     *
     * <p>Any newly created Context will have its box vectors set to these. They will affect any
     * Force added to the System that uses periodic boundary conditions.
     *
     * <p>Triclinic boxes are supported, but the vectors must satisfy certain requirements. In
     * particular, a must point in the x direction, b must point "mostly" in the y direction, and c
     * must point "mostly" in the z direction. See the documentation for details.
     */
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
            context.reinitContext();
            softcoreCreated = true;
          }
          context.setParameter("vdw_lambda", lambdaVDW);
        } else if (amoebaVDWForce != null) {
          context.setParameter("AmoebaVdwLambda", lambdaVDW);
          if (softcoreCreated) {
            // Avoid any updateParametersInContext calls if vdwLambdaTerm is true, but not other
            // alchemical terms.
            if (!torsionLambdaTerm && !elecLambdaTerm) {
              return;
            }
          } else {
            softcoreCreated = true;
          }
        }
      }

      // Note Stretch-Torsion and Angle-Torsion terms (for nucleic acids)
      // and Torsion-Torsion terms (for protein backbones) are not udpated yet.

      if (updateBondedTerms) {
        if (bondForce != null) {
          updateBondForce();
        }
        if (angleForce != null) {
          updateAngleForce();
        }
        if (stretchBendForce != null) {
          updateStretchBendForce();
        }
        if (inPlaneAngleForce != null) {
          updateInPlaneAngleForce();
        }
        if (ureyBradleyForce != null) {
          updateUreyBradleyForce();
        }
        if (outOfPlaneBendForce != null) {
          updateOutOfPlaneBendForce();
        }
        if (piTorsionForce != null) {
          updatePiTorsionForce();
        }
      }

      if (torsionLambdaTerm || updateBondedTerms) {
        if (torsionForce != null) {
          updateTorsionForce();
        }
        if (improperTorsionForce != null) {
          updateImproperTorsionForce();
        }
      }

      if (restraintTorsions != null && restraintTorsions.length > 0) {
        updateRestraintTorsions();
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
        updateGeneralizedKirkwoodForce(atoms);
      }

      // Update WCA Force.
      if (amoebaWcaDispersionForce != null) {
        updateWCAForce(atoms);
      }

      // Update WCA Force.
      if (amoebaCavitationForce != null) {
        updateCavitationForce(atoms);
      }
    }

    /**
     * Adds atoms from the molecular assembly to the OpenMM System and reports to the user the number
     * of particles added.
     */
    private void addAtoms() throws Exception {
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
      }
      logger.log(Level.INFO, format("  Atoms \t\t%6d", nAtoms));
      logger.log(Level.INFO, format("  Mass  \t\t%12.3f", totalMass));
    }

    /** This methods sets the mass of inactive atoms to zero. */
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

    /** Add a bond force to the OpenMM System. */
    private void addBondForce() {
      Bond[] bonds = getBonds();
      if (bonds == null || bonds.length < 1) {
        return;
      }

      String energy;
      if (bonds[0].bondType.bondFunction == BondFunction.QUARTIC) {
        energy = format("k*(d^2 + %.15g*d^3 + %.15g*d^4); d=r-r0",
            BondType.cubic / OpenMM_NmPerAngstrom,
            BondType.quartic / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
      } else {
        energy = "k*(d^2); d=r-r0";
      }
      bondForce = OpenMM_CustomBondForce_create(energy);
      OpenMM_CustomBondForce_addPerBondParameter(bondForce, "r0");
      OpenMM_CustomBondForce_addPerBondParameter(bondForce, "k");
      OpenMM_Force_setName(bondForce, "AmoebaBond");

      double kParameterConversion =
          OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getXyzIndex() - 1;
        int i2 = bond.getAtom(1).getXyzIndex() - 1;
        BondType bondType = bond.bondType;
        double r0 = bondType.distance * OpenMM_NmPerAngstrom;
        double k = kParameterConversion * bondType.forceConstant * BondType.units;
        OpenMM_DoubleArray_append(parameters, r0);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomBondForce_addBond(bondForce, i1, i2, parameters);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_DoubleArray_destroy(parameters);

      int forceGroup = forceField.getInteger("BOND_FORCE_GROUP", 0);
      OpenMM_Force_setForceGroup(bondForce, forceGroup);
      OpenMM_System_addForce(system, bondForce);
      logger.log(Level.INFO, format("  Bonds \t\t%6d\t\t%1d", bonds.length, forceGroup));
    }

    /** Update an existing bond force for the OpenMM System. */
    private void updateBondForce() {
      Bond[] bonds = getBonds();
      if (bonds == null || bonds.length < 1) {
        return;
      }

      double kParameterConversion =
          OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      int index = 0;
      for (Bond bond : bonds) {
        int i1 = bond.getAtom(0).getXyzIndex() - 1;
        int i2 = bond.getAtom(1).getXyzIndex() - 1;
        BondType bondType = bond.bondType;
        double r0 = bondType.distance * OpenMM_NmPerAngstrom;
        double k = kParameterConversion * bondType.forceConstant * BondType.units;
        OpenMM_DoubleArray_append(parameters, r0);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomBondForce_setBondParameters(bondForce, index++, i1, i2, parameters);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_DoubleArray_destroy(parameters);

      if (context.contextPointer != null) {
        OpenMM_CustomBondForce_updateParametersInContext(
            bondForce, context.contextPointer);
      }
    }

    /** Add an angle force to the OpenMM System. */
    private void addAngleForce() {
      Angle[] angles = getAngles();
      if (angles == null || angles.length < 1) {
        return;
      }
      boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
      String energy;
      if (angles[0].angleType.angleFunction == AngleFunction.SEXTIC) {
        energy = format(
            "k*(d^2 + %.15g*d^3 + %.15g*d^4 + %.15g*d^5 + %.15g*d^6); d=%.15g*theta-theta0",
            AngleType.cubic, AngleType.quartic, AngleType.quintic, AngleType.sextic, 180.0 / PI);
      } else {
        energy = format("k*(d^2); d=%.15g*theta-theta0", 180.0 / PI);
      }
      angleForce = OpenMM_CustomAngleForce_create(energy);
      OpenMM_CustomAngleForce_addPerAngleParameter(angleForce, "theta0");
      OpenMM_CustomAngleForce_addPerAngleParameter(angleForce, "k");
      OpenMM_Force_setName(angleForce, "Angle");

      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      int angleCount = 0;
      for (Angle angle : angles) {
        AngleMode angleMode = angle.angleType.angleMode;

        if (!manyBodyTitration && angleMode == AngleMode.IN_PLANE) {
          // Skip In-Plane angles unless this is ManyBody Titration.
          continue;
        } else if (isHydrogenAngle(angle) && rigidHydrogenAngles) {
          logger.log(Level.INFO, " Constrained angle %s was not added the AngleForce.", angle);
        } else {
          int i1 = angle.getAtom(0).getXyzIndex() - 1;
          int i2 = angle.getAtom(1).getXyzIndex() - 1;
          int i3 = angle.getAtom(2).getXyzIndex() - 1;
          double theta0 = angle.angleType.angle[angle.nh];
          double k = OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant;
          if (angleMode == AngleMode.IN_PLANE) {
            // This is a place-holder Angle, in case the In-Plane Angle is swtiched to a
            // Normal Angle during in the udpateAngleForce.
            k = 0.0;
          }
          OpenMM_DoubleArray_append(parameters, theta0);
          OpenMM_DoubleArray_append(parameters, k);
          OpenMM_CustomAngleForce_addAngle(angleForce, i1, i2, i3, parameters);
          angleCount++;
          OpenMM_DoubleArray_resize(parameters, 0);
        }
      }
      OpenMM_DoubleArray_destroy(parameters);

      int forceGroup = forceField.getInteger("ANGLE_FORCE_GROUP", 0);
      OpenMM_Force_setForceGroup(angleForce, forceGroup);
      OpenMM_System_addForce(system, angleForce);
      logger.log(Level.INFO, format("  Angles \t\t%6d\t\t%1d", angleCount, forceGroup));
    }

    /** Update the angle force. */
    private void updateAngleForce() {
      Angle[] angles = getAngles();
      if (angles == null || angles.length < 1) {
        return;
      }
      boolean rigidHydrogenAngles = forceField.getBoolean("RIGID_HYDROGEN_ANGLES", false);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      int index = 0;
      for (Angle angle : angles) {
        AngleMode angleMode = angle.angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleMode.IN_PLANE) {
          // Skip In-Plane angles unless this is ManyBody Titration.
          continue;
        }
        // Update angles that do not involve rigid hydrogen atoms.
        else if (!rigidHydrogenAngles || !isHydrogenAngle(angle)) {
          int i1 = angle.getAtom(0).getXyzIndex() - 1;
          int i2 = angle.getAtom(1).getXyzIndex() - 1;
          int i3 = angle.getAtom(2).getXyzIndex() - 1;
          double theta0 = angle.angleType.angle[angle.nh];
          double k = OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant;
          if (angleMode == AngleMode.IN_PLANE) {
            // Zero the force constant for In-Plane Angles.
            k = 0.0;
          }
          OpenMM_DoubleArray_append(parameters, theta0);
          OpenMM_DoubleArray_append(parameters, k);
          OpenMM_CustomAngleForce_setAngleParameters(angleForce, index++, i1, i2, i3, parameters);
          OpenMM_DoubleArray_resize(parameters, 0);
        }
      }
      OpenMM_DoubleArray_destroy(parameters);

      if (context.contextPointer != null) {
        OpenMM_CustomAngleForce_updateParametersInContext(
            angleForce, context.contextPointer);
      }
    }

    /** Add an in-plane angle force to the OpenMM System. */
    private void addInPlaneAngleForce() {
      Angle[] angles = getAngles();
      if (angles == null || angles.length < 1) {
        return;
      }

      String energy = format(
          "k*(d^2 + %.15g*d^3 + %.15g*d^4 + %.15g*d^5 + %.15g*d^6); d=theta-theta0; "
              + "theta = %.15g*pointangle(x1, y1, z1, projx, projy, projz, x3, y3, z3); "
              + "projx = x2-nx*dot; projy = y2-ny*dot; projz = z2-nz*dot; "
              + "dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3); "
              + "nx = px/norm; ny = py/norm; nz = pz/norm; "
              + "norm = sqrt(px*px + py*py + pz*pz); "
              + "px = (d1y*d2z-d1z*d2y); py = (d1z*d2x-d1x*d2z); pz = (d1x*d2y-d1y*d2x); "
              + "d1x = x1-x4; d1y = y1-y4; d1z = z1-z4; "
              + "d2x = x3-x4; d2y = y3-y4; d2z = z3-z4",
          AngleType.cubic, AngleType.quartic, AngleType.quintic, AngleType.sextic, 180.0 / PI);
      inPlaneAngleForce = OpenMM_CustomCompoundBondForce_create(4, energy);
      OpenMM_CustomCompoundBondForce_addPerBondParameter(inPlaneAngleForce, "theta0");
      OpenMM_CustomCompoundBondForce_addPerBondParameter(inPlaneAngleForce, "k");
      OpenMM_Force_setName(inPlaneAngleForce, "InPlaneAngle");

      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      for (Angle angle : angles) {
        AngleMode angleMode = angle.angleType.angleMode;

        if (!manyBodyTitration && angleMode == AngleMode.NORMAL) {
          // Skip Normal angles unless this is ManyBody Titration.
          continue;
        } else {
          double theta0 = angle.angleType.angle[angle.nh];
          double k = OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant;
          int i1 = angle.getAtom(0).getXyzIndex() - 1;
          int i2 = angle.getAtom(1).getXyzIndex() - 1;
          int i3 = angle.getAtom(2).getXyzIndex() - 1;
          int i4 = 0;
          if (angleMode == AngleMode.NORMAL) {
            // This is a place-holder Angle, in case the Normal Angle is switched to a
            // In-Plane Angle during in the udpateInPlaneAngleForce.
            k = 0.0;
            Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
            if (fourthAtom != null) {
              i4 = fourthAtom.getXyzIndex() - 1;
            } else {
              while (i1 == i4 || i2 == i4 || i3 == i4) {
                i4++;
              }
            }
          } else {
            i4 = angle.getAtom4().getXyzIndex() - 1;
          }
          OpenMM_IntArray_append(particles, i1);
          OpenMM_IntArray_append(particles, i2);
          OpenMM_IntArray_append(particles, i3);
          OpenMM_IntArray_append(particles, i4);
          OpenMM_DoubleArray_append(parameters, theta0);
          OpenMM_DoubleArray_append(parameters, k);
          OpenMM_CustomCompoundBondForce_addBond(inPlaneAngleForce, particles, parameters);
          OpenMM_IntArray_resize(particles, 0);
          OpenMM_DoubleArray_resize(parameters, 0);
        }
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);

      int forceGroup = forceField.getInteger("IN_PLANE_ANGLE_FORCE_GROUP", 0);
      OpenMM_Force_setForceGroup(inPlaneAngleForce, forceGroup);
      OpenMM_System_addForce(system, inPlaneAngleForce);
      logger.log(Level.INFO,
          format("  In-Plane Angles \t%6d\t\t%1d", angles.length, forceGroup));
    }

    /** Update the in-plane angle force. */
    private void updateInPlaneAngleForce() {
      Angle[] angles = getAngles();
      if (angles == null || angles.length < 1) {
        return;
      }
      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      int index = 0;
      for (Angle angle : angles) {
        AngleMode angleMode = angle.angleType.angleMode;
        if (!manyBodyTitration && angleMode == AngleMode.NORMAL) {
          // Skip Normal angles unless this is ManyBody Titration.
          continue;
        } else {
          double theta0 = angle.angleType.angle[angle.nh];
          double k = OpenMM_KJPerKcal * AngleType.units * angle.angleType.forceConstant;
          int i1 = angle.getAtom(0).getXyzIndex() - 1;
          int i2 = angle.getAtom(1).getXyzIndex() - 1;
          int i3 = angle.getAtom(2).getXyzIndex() - 1;
          // There is no 4th atom for normal angles, so set the index to first atom.
          int i4 = 0;
          if (angleMode == AngleMode.NORMAL) {
            // Zero the force constant for Normal Angles.
            k = 0.0;
            Atom fourthAtom = angle.getFourthAtomOfTrigonalCenter();
            if (fourthAtom != null) {
              i4 = fourthAtom.getXyzIndex() - 1;
            } else {
              while (i1 == i4 || i2 == i4 || i3 == i4) {
                i4++;
              }
            }
          } else {
            i4 = angle.getAtom4().getXyzIndex() - 1;
          }
          OpenMM_IntArray_append(particles, i1);
          OpenMM_IntArray_append(particles, i2);
          OpenMM_IntArray_append(particles, i3);
          OpenMM_IntArray_append(particles, i4);
          OpenMM_DoubleArray_append(parameters, theta0);
          OpenMM_DoubleArray_append(parameters, k);
          OpenMM_CustomCompoundBondForce_setBondParameters(inPlaneAngleForce, index++, particles,
              parameters);
          OpenMM_IntArray_resize(particles, 0);
          OpenMM_DoubleArray_resize(parameters, 0);
        }
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);

      if (context.contextPointer != null) {
        OpenMM_CustomCompoundBondForce_updateParametersInContext(
            inPlaneAngleForce, context.contextPointer);
      }
    }

    /** Add a Urey-Bradley force to the OpenMM System. */
    private void addUreyBradleyForce() {
      UreyBradley[] ureyBradleys = getUreyBradleys();
      if (ureyBradleys == null || ureyBradleys.length < 1) {
        return;
      }

      ureyBradleyForce = OpenMM_HarmonicBondForce_create();
      double kParameterConversion =
          UreyBradleyType.units * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

      for (UreyBradley ureyBradley : ureyBradleys) {
        int i1 = ureyBradley.getAtom(0).getXyzIndex() - 1;
        int i2 = ureyBradley.getAtom(2).getXyzIndex() - 1;
        UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
        double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
        // The implementation of UreyBradley in FFX & Tinker: k x^2
        // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
        double k = 2.0 * ureyBradleyType.forceConstant * kParameterConversion;
        OpenMM_HarmonicBondForce_addBond(ureyBradleyForce, i1, i2, length, k);
      }

      int forceGroup = forceField.getInteger("UREY_BRADLEY_FORCE", 0);
      OpenMM_Force_setForceGroup(ureyBradleyForce, forceGroup);
      OpenMM_System_addForce(system, ureyBradleyForce);
      logger.log(Level.INFO,
          format("  Urey-Bradleys \t%6d\t\t%1d", ureyBradleys.length, forceGroup));
    }

    /** Add a Urey-Bradley force to the OpenMM System. */
    private void updateUreyBradleyForce() {
      UreyBradley[] ureyBradleys = getUreyBradleys();
      if (ureyBradleys == null || ureyBradleys.length < 1) {
        return;
      }

      double kParameterConversion =
          UreyBradleyType.units * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

      int index = 0;
      for (UreyBradley ureyBradley : ureyBradleys) {
        int i1 = ureyBradley.getAtom(0).getXyzIndex() - 1;
        int i2 = ureyBradley.getAtom(2).getXyzIndex() - 1;
        UreyBradleyType ureyBradleyType = ureyBradley.ureyBradleyType;
        double length = ureyBradleyType.distance * OpenMM_NmPerAngstrom;
        // The implementation of UreyBradley in FFX & Tinker: k x^2
        // The implementation of Harmonic Bond Force in OpenMM:  k x^2 / 2
        double k = 2.0 * ureyBradleyType.forceConstant * kParameterConversion;
        OpenMM_HarmonicBondForce_setBondParameters(ureyBradleyForce, index++, i1, i2, length, k);
      }

      if (context.contextPointer != null) {
        OpenMM_HarmonicBondForce_updateParametersInContext(
            ureyBradleyForce, context.contextPointer);
      }
    }

    /** Add an out-of-plane bend force to the OpenMM System. */
    private void addOutOfPlaneBendForce() {
      OutOfPlaneBend[] outOfPlaneBends = getOutOfPlaneBends();
      if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
        return;
      }

      String energy = format(
          "k*(theta^2 + %.15g*theta^3 + %.15g*theta^4 + %.15g*theta^5 + %.15g*theta^6); "
              + "theta = %.15g*pointangle(x2, y2, z2, x4, y4, z4, projx, projy, projz); "
              + "projx = x2-nx*dot; projy = y2-ny*dot; projz = z2-nz*dot; "
              + "dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3); "
              + "nx = px/norm; ny = py/norm; nz = pz/norm; "
              + "norm = sqrt(px*px + py*py + pz*pz); "
              + "px = (d1y*d2z-d1z*d2y); py = (d1z*d2x-d1x*d2z); pz = (d1x*d2y-d1y*d2x); "
              + "d1x = x1-x4; d1y = y1-y4; d1z = z1-z4; "
              + "d2x = x3-x4; d2y = y3-y4; d2z = z3-z4",
          OutOfPlaneBendType.cubic, OutOfPlaneBendType.quartic, OutOfPlaneBendType.quintic,
          OutOfPlaneBendType.sextic, 180.0 / PI);
      outOfPlaneBendForce = OpenMM_CustomCompoundBondForce_create(4, energy);
      OpenMM_CustomCompoundBondForce_addPerBondParameter(outOfPlaneBendForce, "k");
      OpenMM_Force_setName(outOfPlaneBendForce, "OutOfPlaneBend");

      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
        OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
        int i1 = outOfPlaneBend.getAtom(0).getXyzIndex() - 1;
        int i2 = outOfPlaneBend.getAtom(1).getXyzIndex() - 1;
        int i3 = outOfPlaneBend.getAtom(2).getXyzIndex() - 1;
        int i4 = outOfPlaneBend.getAtom(3).getXyzIndex() - 1;
        double k = OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * OutOfPlaneBendType.units;
        OpenMM_IntArray_append(particles, i1);
        OpenMM_IntArray_append(particles, i2);
        OpenMM_IntArray_append(particles, i3);
        OpenMM_IntArray_append(particles, i4);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomCompoundBondForce_addBond(outOfPlaneBendForce, particles, parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);
      int forceGroup = forceField.getInteger("OUT_OF_PLANE_BEND_FORCE_GROUP", 0);
      OpenMM_Force_setForceGroup(outOfPlaneBendForce, forceGroup);
      OpenMM_System_addForce(system, outOfPlaneBendForce);
      logger.log(Level.INFO,
          format("  Out-of-Plane Bends \t%6d\t\t%1d", outOfPlaneBends.length, forceGroup));
    }

    /** Update the Out-of-Plane bend force. */
    private void updateOutOfPlaneBendForce() {
      OutOfPlaneBend[] outOfPlaneBends = getOutOfPlaneBends();
      if (outOfPlaneBends == null || outOfPlaneBends.length < 1) {
        return;
      }

      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      int index = 0;
      for (OutOfPlaneBend outOfPlaneBend : outOfPlaneBends) {
        OutOfPlaneBendType outOfPlaneBendType = outOfPlaneBend.outOfPlaneBendType;
        int i1 = outOfPlaneBend.getAtom(0).getXyzIndex() - 1;
        int i2 = outOfPlaneBend.getAtom(1).getXyzIndex() - 1;
        int i3 = outOfPlaneBend.getAtom(2).getXyzIndex() - 1;
        int i4 = outOfPlaneBend.getAtom(3).getXyzIndex() - 1;
        double k = OpenMM_KJPerKcal * outOfPlaneBendType.forceConstant * OutOfPlaneBendType.units;
        OpenMM_IntArray_append(particles, i1);
        OpenMM_IntArray_append(particles, i2);
        OpenMM_IntArray_append(particles, i3);
        OpenMM_IntArray_append(particles, i4);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomCompoundBondForce_setBondParameters(outOfPlaneBendForce, index++, particles,
            parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);

      if (context.contextPointer != null) {
        OpenMM_CustomCompoundBondForce_updateParametersInContext(
            outOfPlaneBendForce, context.contextPointer);
      }

    }

    /** Add a stretch-bend force to the OpenMM System. */
    private void addStretchBendForce() {
      StretchBend[] stretchBends = getStretchBends();
      if (stretchBends == null || stretchBends.length < 1) {
        return;
      }

      String energy = format(
          "(k1*(distance(p1,p2)-r12) + k2*(distance(p2,p3)-r23))*(%.15g*(angle(p1,p2,p3)-theta0))",
          180.0 / PI);
      stretchBendForce = OpenMM_CustomCompoundBondForce_create(3, energy);
      OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "r12");
      OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "r23");
      OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "theta0");
      OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "k1");
      OpenMM_CustomCompoundBondForce_addPerBondParameter(stretchBendForce, "k2");
      OpenMM_Force_setName(stretchBendForce, "AmoebaStretchBend");

      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      for (StretchBend stretchBend : stretchBends) {
        int i1 = stretchBend.getAtom(0).getXyzIndex() - 1;
        int i2 = stretchBend.getAtom(1).getXyzIndex() - 1;
        int i3 = stretchBend.getAtom(2).getXyzIndex() - 1;
        double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
        double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
        double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
        double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
        double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
        OpenMM_IntArray_append(particles, i1);
        OpenMM_IntArray_append(particles, i2);
        OpenMM_IntArray_append(particles, i3);
        OpenMM_DoubleArray_append(parameters, r12);
        OpenMM_DoubleArray_append(parameters, r23);
        OpenMM_DoubleArray_append(parameters, theta0);
        OpenMM_DoubleArray_append(parameters, k1);
        OpenMM_DoubleArray_append(parameters, k2);
        OpenMM_CustomCompoundBondForce_addBond(stretchBendForce, particles, parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);

      int forceGroup = forceField.getInteger("STRETCH_BEND_FORCE_GROUP", 0);
      OpenMM_Force_setForceGroup(stretchBendForce, forceGroup);
      OpenMM_System_addForce(system, stretchBendForce);
      logger.log(
          Level.INFO, format("  Stretch-Bends \t%6d\t\t%1d", stretchBends.length, forceGroup));
    }

    /** Update the Stretch-Bend force. */
    private void updateStretchBendForce() {
      StretchBend[] stretchBends = getStretchBends();
      if (stretchBends == null || stretchBends.length < 1) {
        return;
      }

      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      int index = 0;
      for (StretchBend stretchBend : stretchBends) {
        int i1 = stretchBend.getAtom(0).getXyzIndex() - 1;
        int i2 = stretchBend.getAtom(1).getXyzIndex() - 1;
        int i3 = stretchBend.getAtom(2).getXyzIndex() - 1;
        double r12 = stretchBend.bond0Eq * OpenMM_NmPerAngstrom;
        double r23 = stretchBend.bond1Eq * OpenMM_NmPerAngstrom;
        double theta0 = stretchBend.angleEq * OpenMM_RadiansPerDegree;
        double k1 = stretchBend.force0 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
        double k2 = stretchBend.force1 * OpenMM_KJPerKcal / OpenMM_NmPerAngstrom;
        OpenMM_IntArray_append(particles, i1);
        OpenMM_IntArray_append(particles, i2);
        OpenMM_IntArray_append(particles, i3);
        OpenMM_DoubleArray_append(parameters, r12);
        OpenMM_DoubleArray_append(parameters, r23);
        OpenMM_DoubleArray_append(parameters, theta0);
        OpenMM_DoubleArray_append(parameters, k1);
        OpenMM_DoubleArray_append(parameters, k2);
        OpenMM_CustomCompoundBondForce_setBondParameters(stretchBendForce, index++, particles,
            parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);

      if (context.contextPointer != null) {
        OpenMM_CustomCompoundBondForce_updateParametersInContext(stretchBendForce,
            context.contextPointer);
      }
    }

    /** Add a torsion force to the OpenMM System. */
    private void addTorsionForce() {
      Torsion[] torsions = getTorsions();
      if (torsions == null || torsions.length < 1) {
        return;
      }

      torsionForce = OpenMM_PeriodicTorsionForce_create();
      for (Torsion torsion : torsions) {
        int a1 = torsion.getAtom(0).getXyzIndex() - 1;
        int a2 = torsion.getAtom(1).getXyzIndex() - 1;
        int a3 = torsion.getAtom(2).getXyzIndex() - 1;
        int a4 = torsion.getAtom(3).getXyzIndex() - 1;
        TorsionType torsionType = torsion.torsionType;
        int nTerms = torsionType.phase.length;
        for (int j = 0; j < nTerms; j++) {
          OpenMM_PeriodicTorsionForce_addTorsion(torsionForce, a1, a2, a3, a4, j + 1,
              torsionType.phase[j] * OpenMM_RadiansPerDegree,
              OpenMM_KJPerKcal * torsion.units * torsionType.amplitude[j]);
        }
        // Enforce 6-fold torsions since TorsionType instances can have different lengths
        // when side-chain protonation changes.
        if (manyBodyTitration) {
          for (int j = nTerms; j < 6; j++) {
            OpenMM_PeriodicTorsionForce_addTorsion(torsionForce, a1, a2, a3, a4, j + 1,
                0.0, 0.0);
          }
        }
      }
      int fGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);

      OpenMM_Force_setForceGroup(torsionForce, fGroup);
      OpenMM_System_addForce(system, torsionForce);

      logger.log(Level.INFO, format("  Torsions \t\t%6d\t\t%1d", torsions.length, fGroup));
    }

    /** Update the Torsion force. */
    private void updateTorsionForce() {
      // Check if this system has torsions.
      Torsion[] torsions = getTorsions();
      if (torsions == null || torsions.length < 1) {
        return;
      }

      int index = 0;
      for (Torsion torsion : torsions) {
        TorsionType torsionType = torsion.torsionType;
        int nTerms = torsionType.phase.length;
        int a1 = torsion.getAtom(0).getXyzIndex() - 1;
        int a2 = torsion.getAtom(1).getXyzIndex() - 1;
        int a3 = torsion.getAtom(2).getXyzIndex() - 1;
        int a4 = torsion.getAtom(3).getXyzIndex() - 1;
        for (int j = 0; j < nTerms; j++) {
          double forceConstant =
              OpenMM_KJPerKcal * torsion.units * torsionType.amplitude[j] * lambdaTorsion;
          OpenMM_PeriodicTorsionForce_setTorsionParameters(torsionForce, index++, a1, a2, a3, a4,
              j + 1, torsionType.phase[j] * OpenMM_RadiansPerDegree, forceConstant);
        }
        // Enforce 6-fold torsions since TorsionType instances can have different lengths
        // when side-chain protonation changes.
        if (manyBodyTitration) {
          for (int j = nTerms; j < 6; j++) {
            OpenMM_PeriodicTorsionForce_setTorsionParameters(torsionForce, index++, a1, a2, a3, a4,
                j + 1,
                0.0, 0.0);
          }
        }
      }

      if (context.contextPointer != null) {
        OpenMM_PeriodicTorsionForce_updateParametersInContext(
            torsionForce, context.contextPointer);
      }
    }

    /** Add an improper-torsion force to the OpenMM System. */
    private void addImproperTorsionForce() {
      ImproperTorsion[] improperTorsions = getImproperTorsions();
      if (improperTorsions == null || improperTorsions.length < 1) {
        return;
      }

      improperTorsionForce = OpenMM_PeriodicTorsionForce_create();
      for (ImproperTorsion improperTorsion : improperTorsions) {
        int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
        int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
        int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
        int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
        ImproperTorsionType improperTorsionType = improperTorsion.improperType;
        OpenMM_PeriodicTorsionForce_addTorsion(
            improperTorsionForce,
            a1,
            a2,
            a3,
            a4,
            improperTorsionType.periodicity,
            improperTorsionType.phase * OpenMM_RadiansPerDegree,
            OpenMM_KJPerKcal
                * improperTorsion.units
                * improperTorsion.scaleFactor
                * improperTorsionType.k);
      }

      int forceGroup = forceField.getInteger("IMPROPER_TORSION_FORCE_GROUP", 0);

      OpenMM_Force_setForceGroup(improperTorsionForce, forceGroup);
      OpenMM_System_addForce(system, improperTorsionForce);

      logger.log(
          Level.INFO,
          format("  Improper Torsions \t%6d\t\t%1d", improperTorsions.length, forceGroup));
    }

    /** Update the Improper Torsion force. */
    private void updateImproperTorsionForce() {
      ImproperTorsion[] improperTorsions = getImproperTorsions();
      if (improperTorsions == null || improperTorsions.length < 1) {
        return;
      }

      int nImproperTorsions = improperTorsions.length;
      for (int i = 0; i < nImproperTorsions; i++) {
        ImproperTorsion improperTorsion = improperTorsions[i];
        int a1 = improperTorsion.getAtom(0).getXyzIndex() - 1;
        int a2 = improperTorsion.getAtom(1).getXyzIndex() - 1;
        int a3 = improperTorsion.getAtom(2).getXyzIndex() - 1;
        int a4 = improperTorsion.getAtom(3).getXyzIndex() - 1;
        ImproperTorsionType improperTorsionType = improperTorsion.improperType;
        double forceConstant = OpenMM_KJPerKcal * improperTorsion.units * improperTorsion.scaleFactor
            * improperTorsionType.k * lambdaTorsion;
        OpenMM_PeriodicTorsionForce_setTorsionParameters(improperTorsionForce, i, a1, a2, a3, a4,
            improperTorsionType.periodicity, improperTorsionType.phase * OpenMM_RadiansPerDegree,
            forceConstant);
      }

      if (context.contextPointer != null) {
        OpenMM_PeriodicTorsionForce_updateParametersInContext(
            improperTorsionForce, context.contextPointer);
      }
    }

    /** Add a Pi-Torsion force to the OpenMM System. */
    private void addPiTorsionForce() {
      PiOrbitalTorsion[] piOrbitalTorsions = getPiOrbitalTorsions();
      if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
        return;
      }

      String energy = "2*k*sin(phi)^2;"
          + "phi = pointdihedral(x3+c1x, y3+c1y, z3+c1z, x3, y3, z3, x4, y4, z4, x4+c2x, y4+c2y, z4+c2z); "
          + "c1x = (d14y*d24z-d14z*d24y); c1y = (d14z*d24x-d14x*d24z); c1z = (d14x*d24y-d14y*d24x); "
          + "c2x = (d53y*d63z-d53z*d63y); c2y = (d53z*d63x-d53x*d63z); c2z = (d53x*d63y-d53y*d63x); "
          + "d14x = x1-x4; d14y = y1-y4; d14z = z1-z4; "
          + "d24x = x2-x4; d24y = y2-y4; d24z = z2-z4; "
          + "d53x = x5-x3; d53y = y5-y3; d53z = z5-z3; "
          + "d63x = x6-x3; d63y = y6-y3; d63z = z6-z3";
      piTorsionForce = OpenMM_CustomCompoundBondForce_create(6, energy);
      OpenMM_CustomCompoundBondForce_addPerBondParameter(piTorsionForce, "k");
      OpenMM_Force_setName(piTorsionForce, "PiTorsion");

      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      for (PiOrbitalTorsion piOrbitalTorsion : piOrbitalTorsions) {
        int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
        int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
        int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
        int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
        int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
        int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
        PiOrbitalTorsionType type = piOrbitalTorsion.piOrbitalTorsionType;
        double k = OpenMM_KJPerKcal * type.forceConstant * PiOrbitalTorsionType.units;
        OpenMM_IntArray_append(particles, a1);
        OpenMM_IntArray_append(particles, a2);
        OpenMM_IntArray_append(particles, a3);
        OpenMM_IntArray_append(particles, a4);
        OpenMM_IntArray_append(particles, a5);
        OpenMM_IntArray_append(particles, a6);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomCompoundBondForce_addBond(piTorsionForce, particles, parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);

      int forceGroup = forceField.getInteger("PI_ORBITAL_TORSION_FORCE_GROUP", 0);
      OpenMM_Force_setForceGroup(piTorsionForce, forceGroup);
      OpenMM_System_addForce(system, piTorsionForce);
      logger.log(Level.INFO,
          format("  Pi-Orbital Torsions  \t%6d\t\t%1d", piOrbitalTorsions.length, forceGroup));
    }

    /** Update the Pi-Torsion force. */
    private void updatePiTorsionForce() {
      PiOrbitalTorsion[] piOrbitalTorsions = getPiOrbitalTorsions();
      if (piOrbitalTorsions == null || piOrbitalTorsions.length < 1) {
        return;
      }

      PointerByReference particles = OpenMM_IntArray_create(0);
      PointerByReference parameters = OpenMM_DoubleArray_create(0);
      int index = 0;
      for (PiOrbitalTorsion piOrbitalTorsion : piOrbitalTorsions) {
        int a1 = piOrbitalTorsion.getAtom(0).getXyzIndex() - 1;
        int a2 = piOrbitalTorsion.getAtom(1).getXyzIndex() - 1;
        int a3 = piOrbitalTorsion.getAtom(2).getXyzIndex() - 1;
        int a4 = piOrbitalTorsion.getAtom(3).getXyzIndex() - 1;
        int a5 = piOrbitalTorsion.getAtom(4).getXyzIndex() - 1;
        int a6 = piOrbitalTorsion.getAtom(5).getXyzIndex() - 1;
        PiOrbitalTorsionType type = piOrbitalTorsion.piOrbitalTorsionType;
        double k = OpenMM_KJPerKcal * type.forceConstant * PiOrbitalTorsionType.units;
        OpenMM_IntArray_append(particles, a1);
        OpenMM_IntArray_append(particles, a2);
        OpenMM_IntArray_append(particles, a3);
        OpenMM_IntArray_append(particles, a4);
        OpenMM_IntArray_append(particles, a5);
        OpenMM_IntArray_append(particles, a6);
        OpenMM_DoubleArray_append(parameters, k);
        OpenMM_CustomCompoundBondForce_setBondParameters(piTorsionForce, index++, particles,
            parameters);
        OpenMM_IntArray_resize(particles, 0);
        OpenMM_DoubleArray_resize(parameters, 0);
      }
      OpenMM_IntArray_destroy(particles);
      OpenMM_DoubleArray_destroy(parameters);

      if (context.contextPointer != null) {
        OpenMM_CustomCompoundBondForce_updateParametersInContext(
            piTorsionForce, context.contextPointer);
      }
    }

    /** Add a Torsion-Torsion force to the OpenMM System. */
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
        OpenMM_AmoebaTorsionTorsionForce_addTorsionTorsion(
            amoebaTorsionTorsionForce, ia, ib, ic, id, ie, iChiral, gridIndex);
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
        OpenMM_AmoebaTorsionTorsionForce_setTorsionTorsionGrid(
            amoebaTorsionTorsionForce, gridIndex++, grid3D);
        OpenMM_3D_DoubleArray_destroy(grid3D);
      }
      OpenMM_DoubleArray_destroy(values);

      int forceGroup = forceField.getInteger("TORSION_TORSION_FORCE_GROUP", 0);

      OpenMM_Force_setForceGroup(amoebaTorsionTorsionForce, forceGroup);
      OpenMM_System_addForce(system, amoebaTorsionTorsionForce);
      logger.log(
          Level.INFO,
          format("  Torsion-Torsions  \t%6d\t\t%1d", torsionTorsions.length, forceGroup));
    }

    /** Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic acid model). */
    private void addStretchTorsionForce() {
      StretchTorsion[] stretchTorsions = getStretchTorsions();
      if (stretchTorsions == null || stretchTorsions.length < 1) {
        return;
      }

      PointerByReference stretchTorsionForce =
          OpenMM_CustomCompoundBondForce_create(4, StretchTorsion.stretchTorsionForm());
      OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi1", 0);
      OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi2", Math.PI);
      OpenMM_CustomCompoundBondForce_addGlobalParameter(stretchTorsionForce, "phi3", 0);

      for (int m = 1; m < 4; m++) {
        for (int n = 1; n < 4; n++) {
          OpenMM_CustomCompoundBondForce_addPerBondParameter(
              stretchTorsionForce, String.format("k%d%d", m, n));
        }
      }

      for (int m = 1; m < 4; m++) {
        OpenMM_CustomCompoundBondForce_addPerBondParameter(
            stretchTorsionForce, String.format("b%d", m));
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

        OpenMM_CustomCompoundBondForce_addBond(
            stretchTorsionForce, strTorsParticles, strTorsParams);
        OpenMM_DoubleArray_destroy(strTorsParams);
        OpenMM_IntArray_destroy(strTorsParticles);
      }

      int forceGroup = forceField.getInteger("STRETCH_TORSION_FORCE_GROUP", 0);

      OpenMM_Force_setForceGroup(stretchTorsionForce, forceGroup);
      OpenMM_System_addForce(system, stretchTorsionForce);

      logger.log(
          Level.INFO,
          format("  Stretch-Torsions  \t%6d\t\t%1d", stretchTorsions.length, forceGroup));
    }

    /** Adds stretch-torsion couplings (as defined for the 2017 AMOEBA nucleic acid model). */
    private void addAngleTorsionForce() {
      AngleTorsion[] angleTorsions = getAngleTorsions();
      if (angleTorsions == null || angleTorsions.length < 1) {
        return;
      }

      PointerByReference angleTorsionForce =
          OpenMM_CustomCompoundBondForce_create(4, AngleTorsion.angleTorsionForm());
      OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi1", 0);
      OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi2", Math.PI);
      OpenMM_CustomCompoundBondForce_addGlobalParameter(angleTorsionForce, "phi3", 0);

      for (int m = 1; m < 3; m++) {
        for (int n = 1; n < 4; n++) {
          OpenMM_CustomCompoundBondForce_addPerBondParameter(
              angleTorsionForce, format("k%d%d", m, n));
        }
      }

      for (int m = 1; m < 3; m++) {
        OpenMM_CustomCompoundBondForce_addPerBondParameter(
            angleTorsionForce, format("a%d", m));
      }

      for (AngleTorsion angleTorsion : angleTorsions) {
        double[] constants = angleTorsion.getConstants();
        PointerByReference atorsParams = OpenMM_DoubleArray_create(0);
        for (int m = 0; m < 2; m++) {
          for (int n = 0; n < 3; n++) {
            int index = (3 * m) + n;
            double kmn = constants[index] * OpenMM_KJPerKcal;
            OpenMM_DoubleArray_append(atorsParams, kmn);
          }
        }

        Atom[] atoms = angleTorsion.getAtomArray(true);

        // One thing that concerns me is whether it's correct to get angle[0] instead of angle[num
        // hydrogens].
        // This is the way it is in FFX, but that may be a bug.

        OpenMM_DoubleArray_append(
            atorsParams, angleTorsion.angleType1.angle[0] * OpenMM_RadiansPerDegree);
        OpenMM_DoubleArray_append(
            atorsParams, angleTorsion.angleType2.angle[0] * OpenMM_RadiansPerDegree);

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

      logger.log(
          Level.INFO, format("  Angle-Torsions  \t%6d\t\t%1d", angleTorsions.length, forceGroup));
    }

    private void addRestraintTorsions() {
      if (rTors != null && rTors.length > 0) {
        int nRT = rTors.length;
        restraintTorsions = new PointerByReference[nRT];
        for (int i = 0; i < nRT; i++) {
          PointerByReference rtOMM = OpenMM_PeriodicTorsionForce_create();
          RestraintTorsion rt = rTors[i];
          int a1 = rt.getAtom(0).getXyzIndex() - 1;
          int a2 = rt.getAtom(1).getXyzIndex() - 1;
          int a3 = rt.getAtom(2).getXyzIndex() - 1;
          int a4 = rt.getAtom(3).getXyzIndex() - 1;
          int nTerms = rt.torsionType.terms;
          for (int j = 0; j < nTerms; j++) {
            OpenMM_PeriodicTorsionForce_addTorsion(
                rtOMM,
                a1,
                a2,
                a3,
                a4,
                j + 1,
                rt.torsionType.phase[j] * OpenMM_RadiansPerDegree,
                OpenMM_KJPerKcal * rt.units * rt.torsionType.amplitude[j]);
          }
          int fGroup = forceField.getInteger("TORSION_FORCE_GROUP", 0);

          OpenMM_Force_setForceGroup(rtOMM, fGroup);
          OpenMM_System_addForce(system, rtOMM);
          restraintTorsions[i] = rtOMM;
        }
        logger.info(format(" Added %d restraint torsions to OpenMM.", nRT));
      }
    }

    private void updateRestraintTorsions() {
      // update restraint torsions ONLY here.
      // Only update parameters if torsions are being scaled by lambda.

      // Check if this system has torsions.

      int nRT = restraintTorsions.length;
      for (int i = 0; i < nRT; i++) {
        RestraintTorsion rt = rTors[i];
        PointerByReference rtOMM = restraintTorsions[i];
        if (rt.applyLambda()) {
          int index = 0;
          TorsionType torsionType = rt.torsionType;
          int nTerms = torsionType.phase.length;
          int a1 = rt.getAtom(0).getXyzIndex() - 1;
          int a2 = rt.getAtom(1).getXyzIndex() - 1;
          int a3 = rt.getAtom(2).getXyzIndex() - 1;
          int a4 = rt.getAtom(3).getXyzIndex() - 1;
          for (int j = 0; j < nTerms; j++) {
            double forceConstant =
                OpenMM_KJPerKcal * rt.units * torsionType.amplitude[j] * rt.mapLambda(getLambda());
            OpenMM_PeriodicTorsionForce_setTorsionParameters(
                rtOMM,
                index++,
                a1,
                a2,
                a3,
                a4,
                j + 1,
                torsionType.phase[j] * OpenMM_RadiansPerDegree,
                forceConstant);
          }
        }

        if (context.contextPointer != null) {
          OpenMM_PeriodicTorsionForce_updateParametersInContext(
              rtOMM, context.contextPointer);
        }
      }
    }

    /** Uses arithmetic mean to define sigma and geometric mean for epsilon. */
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
        OpenMM_NonbondedForce_createExceptionsFromBonds(
            fixedChargeNonBondedForce, bondArray, coulomb14Scale, lj14Scale);
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
          OpenMM_NonbondedForce_getExceptionParameters(
              fixedChargeNonBondedForce, i, particle1, particle2, chargeProd, sigma, eps);
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
        OpenMM_NonbondedForce_setNonbondedMethod(
            fixedChargeNonBondedForce,
            OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_NoCutoff);
      } else {
        OpenMM_NonbondedForce_setNonbondedMethod(
            fixedChargeNonBondedForce,
            OpenMM_NonbondedForce_NonbondedMethod.OpenMM_NonbondedForce_PME);

        if (pme != null) {
          // Units of the Ewald coefficient are A^-1; Multiply by AngstromsPerNM to convert to
          // (Nm^-1).
          double aEwald = OpenMM_AngstromsPerNm * pme.getEwaldCoefficient();
          int nx = pme.getReciprocalSpace().getXDim();
          int ny = pme.getReciprocalSpace().getYDim();
          int nz = pme.getReciprocalSpace().getZDim();
          OpenMM_NonbondedForce_setPMEParameters(fixedChargeNonBondedForce, aEwald, nx, ny, nz);
        }

        NonbondedCutoff nonbondedCutoff = vdW.getNonbondedCutoff();
        double off = nonbondedCutoff.off;
        double cut = nonbondedCutoff.cut;
        OpenMM_NonbondedForce_setCutoffDistance(
            fixedChargeNonBondedForce, OpenMM_NmPerAngstrom * off);
        OpenMM_NonbondedForce_setUseSwitchingFunction(fixedChargeNonBondedForce, OpenMM_True);
        if (cut == off) {
          logger.warning(" OpenMM does not properly handle cutoffs where cut == off!");
          if (cut == Double.MAX_VALUE || cut == Double.POSITIVE_INFINITY) {
            logger.info(" Detected infinite or max-value cutoff; setting cut to 1E+40 for OpenMM.");
            cut = 1E40;
          } else {
            logger.info(
                String.format(
                    " Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.",
                    cut, off));
            cut *= 0.99;
          }
        }
        OpenMM_NonbondedForce_setSwitchingDistance(
            fixedChargeNonBondedForce, OpenMM_NmPerAngstrom * cut);
      }

      OpenMM_NonbondedForce_setUseDispersionCorrection(fixedChargeNonBondedForce, OpenMM_False);

      int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);
      int pmeGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
      if (forceGroup != pmeGroup) {
        logger.severe(
            String.format(
                " ERROR: VDW-FORCE-GROUP is %d while PME-FORCE-GROUP is %d. "
                    + "This is invalid for fixed-charge force fields with combined nonbonded forces.",
                forceGroup, pmeGroup));
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
     * Updates the fixed-charge non-bonded force for changes in Use flags or Lambda.
     *
     * @param atoms Array of atoms to update.
     */
    private void updateFixedChargeNonBondedForce(Atom[] atoms) {
      VanDerWaals vdW = getVdwNode();
      // Only 6-12 LJ with arithmetic mean to define sigma and geometric mean for epsilon is
      // supported.
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
          // If we're using vdwLambdaTerm, this atom's vdW interactions are handled by the Custom
          // Non-Bonded force.
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

        OpenMM_NonbondedForce_setParticleParameters(
            fixedChargeNonBondedForce, index, charge, sigma, eps);
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

        OpenMM_NonbondedForce_getExceptionParameters(
            fixedChargeNonBondedForce, i, particle1, particle2, chargeProd, sigma, eps);

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
        OpenMM_NonbondedForce_setExceptionParameters(
            fixedChargeNonBondedForce, i, i1, i2, qq, sigma.getValue(), epsilon);
      }

      if (context.contextPointer != null) {
        OpenMM_NonbondedForce_updateParametersInContext(
            fixedChargeNonBondedForce, context.contextPointer);
      }
    }

    /**
     * 1.) Handle interactions between non-alchemical atoms with our default OpenMM NonBondedForce.
     * Note that alchemical atoms must have eps=0 to turn them off in this force.
     *
     * <p>2.) Handle interactions between alchemical atoms and mixed non-alchemical <-> alchemical
     * interactions with an OpenMM CustomNonBondedForce.
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

        OpenMM_NonbondedForce_getParticleParameters(
            fixedChargeNonBondedForce, index, charge, sigma, eps);
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

      OpenMM_CustomNonbondedForce_addInteractionGroup(
          fixedChargeSoftcore, alchemicalGroup, alchemicalGroup);
      OpenMM_CustomNonbondedForce_addInteractionGroup(
          fixedChargeSoftcore, alchemicalGroup, nonAlchemicalGroup);
      OpenMM_IntSet_destroy(alchemicalGroup);
      OpenMM_IntSet_destroy(nonAlchemicalGroup);

      Crystal crystal = getCrystal();
      if (crystal.aperiodic()) {
        OpenMM_CustomNonbondedForce_setNonbondedMethod(
            fixedChargeSoftcore,
            OpenMM_CustomNonbondedForce_NonbondedMethod.OpenMM_CustomNonbondedForce_NoCutoff);
      } else {
        OpenMM_CustomNonbondedForce_setNonbondedMethod(
            fixedChargeSoftcore,
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
          logger.info(
              format(
                  " Detected cut %8.4g == off %8.4g; scaling cut to 0.99 of off for OpenMM.",
                  cut, off));
          cut *= 0.99;
        }
      }

      OpenMM_CustomNonbondedForce_setCutoffDistance(
          fixedChargeSoftcore, OpenMM_NmPerAngstrom * off);
      OpenMM_CustomNonbondedForce_setUseSwitchingFunction(fixedChargeSoftcore, OpenMM_True);
      OpenMM_CustomNonbondedForce_setSwitchingDistance(
          fixedChargeSoftcore, OpenMM_NmPerAngstrom * cut);

      // Add energy parameter derivative
      // OpenMM_CustomNonbondedForce_addEnergyParameterDerivative(fixedChargeSoftcore,
      // "vdw_lambda");

      int forceGroup = forceField.getInteger("VDW_FORCE_GROUP", 1);

      OpenMM_Force_setForceGroup(fixedChargeSoftcore, forceGroup);
      OpenMM_System_addForce(system, fixedChargeSoftcore);

      // Alchemical with Alchemical could be either softcore or normal interactions (softcore here).
      PointerByReference alchemicalAlchemicalStericsForce =
          OpenMM_CustomBondForce_create(stericsEnergyExpression);

      // Non-Alchemical with Alchemical is essentially always softcore.
      PointerByReference nonAlchemicalAlchemicalStericsForce =
          OpenMM_CustomBondForce_create(stericsEnergyExpression);

      // Currently both are treated the same (so we could condense the code below).
      OpenMM_CustomBondForce_addPerBondParameter(alchemicalAlchemicalStericsForce, "rmin");
      OpenMM_CustomBondForce_addPerBondParameter(alchemicalAlchemicalStericsForce, "epsilon");
      OpenMM_CustomBondForce_addGlobalParameter(
          alchemicalAlchemicalStericsForce, "vdw_lambda", 1.0);
      OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "alpha", alpha);
      OpenMM_CustomBondForce_addGlobalParameter(alchemicalAlchemicalStericsForce, "beta", beta);

      OpenMM_CustomBondForce_addPerBondParameter(nonAlchemicalAlchemicalStericsForce, "rmin");
      OpenMM_CustomBondForce_addPerBondParameter(nonAlchemicalAlchemicalStericsForce, "epsilon");
      OpenMM_CustomBondForce_addGlobalParameter(
          nonAlchemicalAlchemicalStericsForce, "vdw_lambda", 1.0);
      OpenMM_CustomBondForce_addGlobalParameter(
          nonAlchemicalAlchemicalStericsForce, "alpha", alpha);
      OpenMM_CustomBondForce_addGlobalParameter(nonAlchemicalAlchemicalStericsForce, "beta", beta);

      int range = OpenMM_NonbondedForce_getNumExceptions(fixedChargeNonBondedForce);

      IntByReference atomi = new IntByReference();
      IntByReference atomj = new IntByReference();
      int[][] torsionMask = vdW.getMask14();

      for (int i = 0; i < range; i++) {
        OpenMM_NonbondedForce_getExceptionParameters(
            fixedChargeNonBondedForce, i, atomi, atomj, charge, sigma, eps);

        // Omit both Exclusions (1-2, 1-3) and Exceptions (scaled 1-4) from the
        // CustomNonbondedForce.
        OpenMM_CustomNonbondedForce_addExclusion(
            fixedChargeSoftcore, atomi.getValue(), atomj.getValue());

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
          } else if ((atom1.applyLambda() && !atom2.applyLambda())
              || (!atom1.applyLambda() && atom2.applyLambda())) {
            oneAlchemical = true;
          }

          if (bothAlchemical) {
            PointerByReference bondParameters = OpenMM_DoubleArray_create(0);
            OpenMM_DoubleArray_append(bondParameters, sigma.getValue() * 1.122462048309372981);
            OpenMM_DoubleArray_append(bondParameters, eps.getValue());
            OpenMM_CustomBondForce_addBond(
                alchemicalAlchemicalStericsForce,
                atomi.getValue(),
                atomj.getValue(),
                bondParameters);
            OpenMM_DoubleArray_destroy(bondParameters);
          } else if (oneAlchemical) {
            PointerByReference bondParameters = OpenMM_DoubleArray_create(0);
            OpenMM_DoubleArray_append(bondParameters, sigma.getValue() * 1.122462048309372981);
            OpenMM_DoubleArray_append(bondParameters, eps.getValue());
            OpenMM_CustomBondForce_addBond(
                nonAlchemicalAlchemicalStericsForce,
                atomi.getValue(),
                atomj.getValue(),
                bondParameters);
            OpenMM_DoubleArray_destroy(bondParameters);
          }
        }
      }

      OpenMM_Force_setForceGroup(alchemicalAlchemicalStericsForce, forceGroup);
      OpenMM_Force_setForceGroup(nonAlchemicalAlchemicalStericsForce, forceGroup);

      // OpenMM_CustomBondForce_addEnergyParameterDerivative(alchemicalAlchemicalStericsForce,
      // "vdw_lambda");
      // OpenMM_CustomBondForce_addEnergyParameterDerivative(nonAlchemicalAlchemicalStericsForce,
      // "vdw_lambda");

      OpenMM_System_addForce(system, alchemicalAlchemicalStericsForce);
      OpenMM_System_addForce(system, nonAlchemicalAlchemicalStericsForce);

      logger.log(Level.INFO, format("  Added fixed charge softcore force \t%d", forceGroup));
      logger.log(Level.INFO, format("   Alpha = %8.6f and beta = %8.6f", alpha, beta));
    }

    /** Add a custom GB force to the OpenMM System. */
    private void addCustomGBForce() {
      GeneralizedKirkwood gk = getGK();
      if (gk == null) {
        return;
      }

      double sTens = 0.0;
      if (gk.getNonPolarModel() == NonPolar.BORN_SOLV
          || gk.getNonPolarModel() == NonPolar.BORN_CAV_DISP) {
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

      OpenMM_CustomGBForce_addGlobalParameter(
          customGBForce, "solventDielectric", gk.getSolventPermittivity());
      OpenMM_CustomGBForce_addGlobalParameter(customGBForce, "soluteDielectric", 1.0);
      OpenMM_CustomGBForce_addGlobalParameter(
          customGBForce,
          "dOffset",
          gk.getDielecOffset() * OpenMM_NmPerAngstrom); // Factor of 0.1 for Ang to nm.
      OpenMM_CustomGBForce_addGlobalParameter(
          customGBForce, "probeRadius", gk.getProbeRadius() * OpenMM_NmPerAngstrom);

      OpenMM_CustomGBForce_addComputedValue(
          customGBForce,
          "I",
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

      OpenMM_CustomGBForce_addComputedValue(
          customGBForce,
          "B",
          // "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
          // "psi=I*or; or=radius-0.009"
          "step(BB-radius)*BB + (1 - step(BB-radius))*radius;"
              + "BB = 1 / ( (3.0*III)^(1.0/3.0) );"
              + "III = step(II)*II + (1 - step(II))*1.0e-9/3.0;"
              + "II = maxI - I;"
              + "maxI = 1/(3.0*radius^3)",
          OpenMM_CustomGBForce_SingleParticle);

      OpenMM_CustomGBForce_addEnergyTerm(
          customGBForce,
          "surfaceTension*(radius+probeRadius+dOffset)^2*((radius+dOffset)/B)^6/6-0.5*138.935456*(1/soluteDielectric-1/solventDielectric)*q^2/B",
          OpenMM_CustomGBForce_SingleParticle);

      // Particle pair term is the generalized Born cross term.
      OpenMM_CustomGBForce_addEnergyTerm(
          customGBForce,
          "-138.935456*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
              + "f=sqrt(r^2+B1*B2*exp(-r^2/(2.455*B1*B2)))",
          OpenMM_CustomGBForce_ParticlePair);

      double[] baseRadii = gk.getBaseRadii();
      double[] overlapScale = gk.getOverlapScale();
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
      if (gk.getNonPolarModel() == NonPolar.BORN_SOLV
          || gk.getNonPolarModel() == NonPolar.BORN_CAV_DISP) {
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

      if (context.contextPointer != null) {
        OpenMM_CustomGBForce_updateParametersInContext(customGBForce, context.contextPointer);
      }
    }

    /** Add an AMOEBA van der Waals force to the OpenMM System. */
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

      /*
        The vdW class used to specify no vdW interactions for an atom will be Zero
        if all atom classes are greater than zero. Otherwise:
        vdWClassForNoInteraction = min(atomClass) - 1
       */
      vdWClassForNoInteraction = 0;
      // Add vdW parameters to the force and record their type.
      vdwClassToOpenMMType = new HashMap<>();

      Map<String, VDWType> vdwTypes = forceField.getVDWTypes();

      for (VDWType vdwType : vdwTypes.values()) {
        int atomClass = vdwType.atomClass;
        if (!vdwClassToOpenMMType.containsKey(atomClass)) {
          double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
          double rad = OpenMM_NmPerAngstrom * vdwType.radius * radScale;
          int type = OpenMM_AmoebaVdwForce_addParticleType(amoebaVDWForce, rad, eps);
          vdwClassToOpenMMType.put(atomClass, type);
          if (atomClass <= vdWClassForNoInteraction) {
            vdWClassForNoInteraction = atomClass - 1;
          }
        }
      }

      // Add a special vdW type for zero vdW energy and forces (e.g. to support the FFX "use" flag).
      int type = OpenMM_AmoebaVdwForce_addParticleType(amoebaVDWForce, OpenMM_NmPerAngstrom, 0.0);
      vdwClassToOpenMMType.put(vdWClassForNoInteraction, type);

      Map<String, VDWPairType> vdwPairTypeMap = forceField.getVDWPairTypes();
      for (VDWPairType vdwPairType : vdwPairTypeMap.values()) {
        int c1 = vdwPairType.atomClasses[0];
        int c2 = vdwPairType.atomClasses[1];
        int type1 = vdwClassToOpenMMType.get(c1);
        int type2 = vdwClassToOpenMMType.get(c2);
        double rMin = vdwPairType.radius * OpenMM_NmPerAngstrom;
        double eps = vdwPairType.wellDepth * OpenMM_KJPerKcal;
        OpenMM_AmoebaVdwForce_addTypePair(amoebaVDWForce, type1, type2, rMin, eps);
        OpenMM_AmoebaVdwForce_addTypePair(amoebaVDWForce, type2, type1, rMin, eps);
      }

      ExtendedSystem extendedSystem = vdW.getExtendedSystem();
      double[] vdwPrefactorAndDerivs = new double[3];

      int[] ired = vdW.getReductionIndex();
      for (int i = 0; i < nAtoms; i++) {
        Atom atom = atoms[i];
        VDWType vdwType = atom.getVDWType();
        int atomClass = vdwType.atomClass;
        type = vdwClassToOpenMMType.get(atomClass);
        int isAlchemical = atom.applyLambda() ? 1 : 0;
        double scaleFactor = 1.0;
        if (extendedSystem != null) {
          extendedSystem.getVdwPrefactor(i, vdwPrefactorAndDerivs);
          scaleFactor = vdwPrefactorAndDerivs[0];
        }

        OpenMM_AmoebaVdwForce_addParticle_1(amoebaVDWForce, ired[i], type, vdwType.reductionFactor,
            isAlchemical, scaleFactor);
      }

      double cutoff = nonbondedCutoff.off * OpenMM_NmPerAngstrom;
      OpenMM_AmoebaVdwForce_setCutoffDistance(amoebaVDWForce, cutoff);

      if (vdW.getDoLongRangeCorrection()) {
        OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_True);
      } else {
        OpenMM_AmoebaVdwForce_setUseDispersionCorrection(amoebaVDWForce, OpenMM_False);
      }

      if (crystal.aperiodic()) {
        OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce, OpenMM_AmoebaVdwForce_NoCutoff);
      } else {
        OpenMM_AmoebaVdwForce_setNonbondedMethod(amoebaVDWForce,
            OpenMM_AmoebaVdwForce_CutoffPeriodic);
      }

      if (vdwLambdaTerm) {
        OpenMM_AmoebaVdwForce_setAlchemicalMethod(amoebaVDWForce, OpenMM_AmoebaVdwForce_Decouple);
        OpenMM_AmoebaVdwForce_setSoftcoreAlpha(amoebaVDWForce, vdWSoftcoreAlpha);
        OpenMM_AmoebaVdwForce_setSoftcorePower(amoebaVDWForce, (int) vdwSoftcorePower);
      }

      int[][] bondMask = vdW.getMask12();
      int[][] angleMask = vdW.getMask13();

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

      ExtendedSystem extendedSystem = vdW.getExtendedSystem();
      double[] vdwPrefactorAndDerivs = new double[3];

      int[] ired = vdW.getReductionIndex();
      for (Atom atom : atoms) {
        int index = atom.getXyzIndex() - 1;
        VDWType vdwType = atom.getVDWType();

        // Get the OpenMM index for this vdW type.
        int type = vdwClassToOpenMMType.get(vdwType.atomClass);
        if (!atom.getUse()) {
          // Get the OpenMM index for a special vdW type that has no interactions.
          type = vdwClassToOpenMMType.get(vdWClassForNoInteraction);
        }
        int isAlchemical = atom.applyLambda() ? 1 : 0;
        double eps = OpenMM_KJPerKcal * vdwType.wellDepth;
        double rad = OpenMM_NmPerAngstrom * vdwType.radius * radScale;

        double scaleFactor = 1.0;
        if (extendedSystem != null) {
          extendedSystem.getVdwPrefactor(index, vdwPrefactorAndDerivs);
          scaleFactor = vdwPrefactorAndDerivs[0];
        }

        OpenMM_AmoebaVdwForce_setParticleParameters(
            amoebaVDWForce, index, ired[index],
            rad, eps, vdwType.reductionFactor, isAlchemical, type, scaleFactor);
      }

      if (context.contextPointer != null) {
        OpenMM_AmoebaVdwForce_updateParametersInContext(amoebaVDWForce, context.contextPointer);
      }
    }

    /** Add an AMOEBA polarizable multipole force to the OpenMM System. */
    private void addAmoebaMultipoleForce() {
      ParticleMeshEwald pme = getPmeNode();
      if (pme == null) {
        return;
      }

      int[][] axisAtom = pme.getAxisAtoms();
      double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
      double polarityConversion =
          OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
      double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

      amoebaMultipoleForce = OpenMM_AmoebaMultipoleForce_create();

      double polarScale = 1.0;
      SCFAlgorithm scfAlgorithm = null;

      if (pme.getPolarizationType() != ParticleMeshEwald.Polarization.MUTUAL) {
        OpenMM_AmoebaMultipoleForce_setPolarizationType(
            amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Direct);
        if (pme.getPolarizationType() == ParticleMeshEwald.Polarization.NONE) {
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
            OpenMM_AmoebaMultipoleForce_setPolarizationType(
                amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Extrapolated);
            PointerByReference exptCoefficients = OpenMM_DoubleArray_create(4);
            OpenMM_DoubleArray_set(exptCoefficients, 0, -0.154);
            OpenMM_DoubleArray_set(exptCoefficients, 1, 0.017);
            OpenMM_DoubleArray_set(exptCoefficients, 2, 0.657);
            OpenMM_DoubleArray_set(exptCoefficients, 3, 0.475);
            OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(
                amoebaMultipoleForce, exptCoefficients);
            OpenMM_DoubleArray_destroy(exptCoefficients);
            break;
          case CG:
          case SOR:
          default:
            OpenMM_AmoebaMultipoleForce_setPolarizationType(
                amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_Mutual);
            break;
        }
      }

      PointerByReference dipoles = OpenMM_DoubleArray_create(3);
      PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

      for (int i = 0; i < nAtoms; i++) {
        Atom atom = atoms[i];
        MultipoleType multipoleType = pme.getMultipoleType(i);
        PolarizeType polarType = pme.getPolarizeType(i);

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
          OpenMM_DoubleArray_set(
              dipoles, j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
        }
        int l = 0;
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            OpenMM_DoubleArray_set(
                quadrupoles,
                l++,
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
        OpenMM_AmoebaMultipoleForce_addMultipole(
            amoebaMultipoleForce,
            charge,
            dipoles,
            quadrupoles,
            axisType,
            zaxis,
            xaxis,
            yaxis,
            polarType.thole,
            polarType.pdamp * dampingFactorConversion,
            polarType.polarizability * polarityConversion * polarScale);
      }
      OpenMM_DoubleArray_destroy(dipoles);
      OpenMM_DoubleArray_destroy(quadrupoles);

      Crystal crystal = getCrystal();
      if (!crystal.aperiodic()) {
        OpenMM_AmoebaMultipoleForce_setNonbondedMethod(
            amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_PME);
        OpenMM_AmoebaMultipoleForce_setCutoffDistance(
            amoebaMultipoleForce, pme.getEwaldCutoff() * OpenMM_NmPerAngstrom);
        OpenMM_AmoebaMultipoleForce_setAEwald(
            amoebaMultipoleForce, pme.getEwaldCoefficient() / OpenMM_NmPerAngstrom);

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
        OpenMM_AmoebaMultipoleForce_setNonbondedMethod(
            amoebaMultipoleForce, OpenMM_AmoebaMultipoleForce_NoCutoff);
      }

      OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(amoebaMultipoleForce, 500);
      OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(
          amoebaMultipoleForce, pme.getPolarEps());

      int[][] ip11 = pme.getPolarization11();

      PointerByReference covalentMap = OpenMM_IntArray_create(0);
      for (int i = 0; i < nAtoms; i++) {
        Atom ai = atoms[i];

        // 1-2 Mask
        OpenMM_IntArray_resize(covalentMap, 0);
        for (Atom ak : ai.get12List()) {
          OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
        }
        OpenMM_AmoebaMultipoleForce_setCovalentMap(
            amoebaMultipoleForce, i, OpenMM_AmoebaMultipoleForce_Covalent12, covalentMap);

        // 1-3 Mask
        OpenMM_IntArray_resize(covalentMap, 0);
        for (Atom ak : ai.get13List()) {
          OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
        }
        OpenMM_AmoebaMultipoleForce_setCovalentMap(
            amoebaMultipoleForce, i, OpenMM_AmoebaMultipoleForce_Covalent13, covalentMap);

        // 1-4 Mask
        OpenMM_IntArray_resize(covalentMap, 0);
        for (Atom ak : ai.get14List()) {
          OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
        }
        OpenMM_AmoebaMultipoleForce_setCovalentMap(
            amoebaMultipoleForce, i, OpenMM_AmoebaMultipoleForce_Covalent14, covalentMap);

        // 1-5 Mask
        OpenMM_IntArray_resize(covalentMap, 0);
        for (Atom ak : ai.get15List()) {
          OpenMM_IntArray_append(covalentMap, ak.getIndex() - 1);
        }
        OpenMM_AmoebaMultipoleForce_setCovalentMap(
            amoebaMultipoleForce, i, OpenMM_AmoebaMultipoleForce_Covalent15, covalentMap);

        // 1-1 Polarization Groups.
        OpenMM_IntArray_resize(covalentMap, 0);
        for (int j = 0; j < ip11[i].length; j++) {
          OpenMM_IntArray_append(covalentMap, ip11[i][j]);
        }
        OpenMM_AmoebaMultipoleForce_setCovalentMap(
            amoebaMultipoleForce,
            i,
            OpenMM_AmoebaMultipoleForce_PolarizationCovalent11,
            covalentMap);

        // AMOEBA does not scale between 1-2, 1-3, etc polarization groups.
      }

      OpenMM_IntArray_destroy(covalentMap);

      int forceGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
      OpenMM_Force_setForceGroup(amoebaMultipoleForce, forceGroup);
      OpenMM_System_addForce(system, amoebaMultipoleForce);

      logger.log(Level.INFO, format("  AMOEBA polarizable multipole force \t%d", forceGroup));

      GeneralizedKirkwood gk = getGK();
      if (gk != null) {
        addGeneralizedKirkwoodForce();
      }

      if (scfAlgorithm == ParticleMeshEwald.SCFAlgorithm.EPT) {
        logger.info("   Using extrapolated perturbation theory for polarization energy.");
      }
    }

    /**
     * Updates the Amoeba electrostatic multipolar force for changes in Use flags or Lambda.
     *
     * @param atoms Array of atoms to update.
     */
    private void updateAmoebaMultipoleForce(Atom[] atoms) {
      ParticleMeshEwald pme = getPmeNode();
      double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
      double polarityConversion =
          OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
      double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

      double polarScale = 1.0;
      if (pme.getPolarizationType() == ParticleMeshEwald.Polarization.NONE) {
        polarScale = 0.0;
      }

      PointerByReference dipoles = OpenMM_DoubleArray_create(3);
      PointerByReference quadrupoles = OpenMM_DoubleArray_create(9);

      for (Atom atom : atoms) {
        int index = atom.getXyzIndex() - 1;
        MultipoleType multipoleType = pme.getMultipoleType(index);
        PolarizeType polarizeType = pme.getPolarizeType(index);
        int[] axisAtoms = atom.getAxisAtomIndices();

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
          OpenMM_DoubleArray_set(
              dipoles, j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
        }
        int l = 0;
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            OpenMM_DoubleArray_set(
                quadrupoles,
                l++,
                multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0 * useFactor);
          }
        }

        int zaxis = 1;
        int xaxis = 1;
        int yaxis = 1;

        if (axisAtoms != null) {
          zaxis = axisAtoms[0];
          if (axisAtoms.length > 1) {
            xaxis = axisAtoms[1];
            if (axisAtoms.length > 2) {
              yaxis = axisAtoms[2];
            }
          }
        } else {
          axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
        }

        // Set the multipole parameters.
        OpenMM_AmoebaMultipoleForce_setMultipoleParameters(
            amoebaMultipoleForce,
            index,
            multipoleType.charge * useFactor,
            dipoles,
            quadrupoles,
            axisType,
            zaxis,
            xaxis,
            yaxis,
            polarizeType.thole,
            polarizeType.pdamp * dampingFactorConversion,
            polarizeType.polarizability * polarityConversion * polarScale * useFactor);
      }

      OpenMM_DoubleArray_destroy(dipoles);
      OpenMM_DoubleArray_destroy(quadrupoles);

      if (context.contextPointer != null) {
        OpenMM_AmoebaMultipoleForce_updateParametersInContext(
            amoebaMultipoleForce, context.contextPointer);
      }
    }

    /** Add a Generalized Kirkwood force to the OpenMM System. */
    private void addGeneralizedKirkwoodForce() {
      GeneralizedKirkwood gk = getGK();

      amoebaGeneralizedKirkwoodForce = OpenMM_AmoebaGeneralizedKirkwoodForce_create();
      OpenMM_AmoebaGeneralizedKirkwoodForce_setSolventDielectric(
          amoebaGeneralizedKirkwoodForce, gk.getSolventPermittivity());
      OpenMM_AmoebaGeneralizedKirkwoodForce_setSoluteDielectric(
          amoebaGeneralizedKirkwoodForce, 1.0);
      OpenMM_AmoebaGeneralizedKirkwoodForce_setDielectricOffset(
          amoebaGeneralizedKirkwoodForce, gk.getDescreenOffset() * OpenMM_NmPerAngstrom);

      boolean usePerfectRadii = gk.getUsePerfectRadii();
      double perfectRadiiScale = 1.0;
      if (usePerfectRadii) {
        // No descreening when using perfect radii (OpenMM will just load the base radii).
        perfectRadiiScale = 0.0;
      }

      // Turn on tanh rescaling only when not using perfect radii.
      int tanhRescale = 0;
      if (gk.getTanhCorrection() && !usePerfectRadii) {
        tanhRescale = 1;
      }
      double[] betas = gk.getTanhBetas();
      OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhRescaling(amoebaGeneralizedKirkwoodForce,
          tanhRescale);
      OpenMM_AmoebaGeneralizedKirkwoodForce_setTanhParameters(amoebaGeneralizedKirkwoodForce,
          betas[0], betas[1], betas[2]);

      double[] baseRadius = gk.getBaseRadii();
      if (usePerfectRadii) {
        baseRadius = gk.getPerfectRadii();
      }

      double[] overlapScale = gk.getOverlapScale();
      double[] descreenRadius = gk.getDescreenRadii();
      double[] neckFactor = gk.getNeckScale();

      if (!usePerfectRadii && logger.isLoggable(Level.FINE)) {
        logger.fine("   GK Base Radii  Descreen Radius  Overlap Scale  Overlap");
      }

      for (int i = 0; i < nAtoms; i++) {
        MultipoleType multipoleType = atoms[i].getMultipoleType();
        double base = baseRadius[i] * OpenMM_NmPerAngstrom;
        double descreen = descreenRadius[i] * OpenMM_NmPerAngstrom * perfectRadiiScale;
        double overlap = overlapScale[i] * perfectRadiiScale;
        double neck = neckFactor[i] * perfectRadiiScale;
        OpenMM_AmoebaGeneralizedKirkwoodForce_addParticle_1(
            amoebaGeneralizedKirkwoodForce, multipoleType.charge, base, overlap, descreen, neck);

        if (!usePerfectRadii && logger.isLoggable(Level.FINE)) {
          logger.fine(format("   %s %8.6f %8.6f %5.3f",
              atoms[i].toString(), baseRadius[i], descreenRadius[i], overlapScale[i]));
        }
      }

      OpenMM_AmoebaGeneralizedKirkwoodForce_setProbeRadius(
          amoebaGeneralizedKirkwoodForce, gk.getProbeRadius() * OpenMM_NmPerAngstrom);

      NonPolar nonpolar = gk.getNonPolarModel();
      switch (nonpolar) {
        case BORN_SOLV:
        case BORN_CAV_DISP:
        default:
          // Configure a Born Radii based surface area term.
          double surfaceTension =
              gk.getSurfaceTension()
                  * OpenMM_KJPerKcal
                  * OpenMM_AngstromsPerNm
                  * OpenMM_AngstromsPerNm;
          OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(
              amoebaGeneralizedKirkwoodForce, OpenMM_True);
          OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(
              amoebaGeneralizedKirkwoodForce, -surfaceTension);
          break;
        case CAV:
        case CAV_DISP:
        case GAUSS_DISP:
        case SEV_DISP:
        case HYDROPHOBIC_PMF:
        case NONE:
          // This NonPolar model does not use a Born Radii based surface area term.
          OpenMM_AmoebaGeneralizedKirkwoodForce_setIncludeCavityTerm(
              amoebaGeneralizedKirkwoodForce, OpenMM_False);
          break;
      }

      int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
      OpenMM_Force_setForceGroup(amoebaGeneralizedKirkwoodForce, forceGroup);
      OpenMM_System_addForce(system, amoebaGeneralizedKirkwoodForce);

      logger.log(Level.INFO, format("  Generalized Kirkwood force \t\t%d", forceGroup));

      // Add dispersion
      switch (nonpolar) {
        case CAV_DISP:
        case GAUSS_DISP:
        case SEV_DISP:
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

      // Add cavitation
      if (nonpolar == GAUSS_DISP) {
        addCavitationForce();
      }
    }

    /**
     * Updates the AMOEBA Generalized Kirkwood force for changes in Use flags or Lambda.
     *
     * @param atoms Array of atoms to update.
     */
    private void updateGeneralizedKirkwoodForce(Atom[] atoms) {
      GeneralizedKirkwood gk = getGK();

      for (int i = 0; i < nAtoms; i++) {
        gk.udpateSoluteParameters(i);
      }

      boolean usePerfectRadii = gk.getUsePerfectRadii();
      double perfectRadiiScale = 1.0;
      if (usePerfectRadii) {
        // No descreening when using perfect radii (OpenMM will just load the base radii).
        perfectRadiiScale = 0.0;
      }

      double[] baseRadii = gk.getBaseRadii();
      if (usePerfectRadii) {
        baseRadii = gk.getPerfectRadii();
      }
      double[] overlapScale = gk.getOverlapScale();
      double[] descreenRadius = gk.getDescreenRadii();
      double[] neckFactors = gk.getNeckScale();

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

        double baseSize = baseRadii[index] * OpenMM_NmPerAngstrom;
        double descreenSize = descreenRadius[index] * OpenMM_NmPerAngstrom * perfectRadiiScale;

        chargeUseFactor *= lambdaScale;
        double overlapScaleUseFactor = nea ? 1.0 : chargeUseFactor;
        overlapScaleUseFactor = overlapScaleUseFactor * perfectRadiiScale;
        double overlap = overlapScale[index] * overlapScaleUseFactor;
        double neckFactor = neckFactors[index] * overlapScaleUseFactor;

        MultipoleType multipoleType = atom.getMultipoleType();
        OpenMM_AmoebaGeneralizedKirkwoodForce_setParticleParameters_1(
            amoebaGeneralizedKirkwoodForce,
            index, multipoleType.charge * chargeUseFactor, baseSize, overlap, descreenSize,
            neckFactor);
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
      //
      // OpenMM_AmoebaGeneralizedKirkwoodForce_setSurfaceAreaFactor(amoebaGeneralizedKirkwoodForce,
      // -surfaceTension);
      //                break;
      //            case CAV:
      //            case CAV_DISP:
      //            case HYDROPHOBIC_PMF:
      //            case NONE:
      //                // This NonPolar model does not use a Born Radii based surface area term.
      //                break;
      //        }

      if (context.contextPointer != null) {
        OpenMM_AmoebaGeneralizedKirkwoodForce_updateParametersInContext(
            amoebaGeneralizedKirkwoodForce, context.contextPointer);
      }
    }

    /** Add a nonpolar Weeks-Chandler-Andersen dispersion force to the OpenMM System. */
    private void addWCAForce() {

      GeneralizedKirkwood gk = getGK();
      DispersionRegion dispersionRegion = gk.getDispersionRegion();

      double epso = 0.1100;
      double epsh = 0.0135;
      double rmino = 1.7025;
      double rminh = 1.3275;
      double awater = 0.033428;
      double slevy = 1.0;
      double dispoff = dispersionRegion.getDispersionOffset();
      double shctd = dispersionRegion.getDispersionOverlapFactor();

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
        OpenMM_AmoebaWcaDispersionForce_addParticle(
            amoebaWcaDispersionForce,
            OpenMM_NmPerAngstrom * radius * radScale,
            OpenMM_KJPerKcal * eps);
      }

      OpenMM_AmoebaWcaDispersionForce_setEpso(amoebaWcaDispersionForce, epso * OpenMM_KJPerKcal);
      OpenMM_AmoebaWcaDispersionForce_setEpsh(amoebaWcaDispersionForce, epsh * OpenMM_KJPerKcal);
      OpenMM_AmoebaWcaDispersionForce_setRmino(
          amoebaWcaDispersionForce, rmino * OpenMM_NmPerAngstrom);
      OpenMM_AmoebaWcaDispersionForce_setRminh(
          amoebaWcaDispersionForce, rminh * OpenMM_NmPerAngstrom);
      OpenMM_AmoebaWcaDispersionForce_setDispoff(
          amoebaWcaDispersionForce, dispoff * OpenMM_NmPerAngstrom);
      OpenMM_AmoebaWcaDispersionForce_setAwater(
          amoebaWcaDispersionForce,
          awater / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom));
      OpenMM_AmoebaWcaDispersionForce_setSlevy(amoebaWcaDispersionForce, slevy);
      OpenMM_AmoebaWcaDispersionForce_setShctd(amoebaWcaDispersionForce, shctd);

      int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);

      OpenMM_Force_setForceGroup(amoebaWcaDispersionForce, forceGroup);
      OpenMM_System_addForce(system, amoebaWcaDispersionForce);

      logger.log(Level.INFO, format("  WCA dispersion force \t\t\t%d", forceGroup));
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
        OpenMM_AmoebaWcaDispersionForce_setParticleParameters(
            amoebaWcaDispersionForce,
            index,
            OpenMM_NmPerAngstrom * radius * radScale,
            OpenMM_KJPerKcal * eps * useFactor);
      }

      if (context.contextPointer != null) {
        OpenMM_AmoebaWcaDispersionForce_updateParametersInContext(
            amoebaWcaDispersionForce, context.contextPointer);
      }
    }

    /** Add a GaussVol cavitation force to the OpenMM System. */
    private void addCavitationForce() {

      GeneralizedKirkwood generalizedKirkwood = getGK();
      ChandlerCavitation chandlerCavitation = generalizedKirkwood.getChandlerCavitation();
      GaussVol gaussVol = chandlerCavitation.getGaussVol();
      if (gaussVol == null) {
        return;
      }

      amoebaCavitationForce = OpenMM_AmoebaGKCavitationForce_create();
      double surfaceTension =
          chandlerCavitation.getSurfaceTension()
              * OpenMM_KJPerKcal
              / OpenMM_NmPerAngstrom
              / OpenMM_NmPerAngstrom;
      double[] rad = gaussVol.getRadii();

      int index = 0;
      for (Atom atom : atoms) {
        int isHydrogen = OpenMM_False;
        double radius = rad[index++];
        if (atom.isHydrogen()) {
          isHydrogen = OpenMM_True;
          radius = 0.0;
        }
        OpenMM_AmoebaGKCavitationForce_addParticle(
            amoebaCavitationForce,
            radius * OpenMM_NmPerAngstrom,
            surfaceTension,
            isHydrogen);
      }

      OpenMM_AmoebaGKCavitationForce_setNonbondedMethod(
          amoebaCavitationForce, OpenMM_AmoebaGKCavitationForce_NoCutoff);

      int forceGroup = forceField.getInteger("GK_FORCE_GROUP", 1);
      OpenMM_Force_setForceGroup(amoebaCavitationForce, forceGroup);
      OpenMM_System_addForce(system, amoebaCavitationForce);

      logger.log(Level.INFO, format("  GaussVol cavitation force \t\t%d", forceGroup));
    }

    /**
     * Updates the Cavitation force for changes in Use flags or Lambda.
     *
     * @param atoms Array of atoms to update.
     */
    private void updateCavitationForce(Atom[] atoms) {
      GeneralizedKirkwood generalizedKirkwood = getGK();
      ChandlerCavitation chandlerCavitation = generalizedKirkwood.getChandlerCavitation();
      GaussVol gaussVol = chandlerCavitation.getGaussVol();
      if (gaussVol == null) {
        return;
      }

      double surfaceTension =
          chandlerCavitation.getSurfaceTension()
              * OpenMM_KJPerKcal
              / OpenMM_NmPerAngstrom
              / OpenMM_NmPerAngstrom;

      // Changing cavitation radii is not supported.
      // for (int i=0; i<nAtoms; i++) {
      //  gaussVol.updateAtom(i);
      // }
      double[] rad = gaussVol.getRadii();

      for (Atom atom : atoms) {
        int index = atom.getXyzIndex() - 1;
        double useFactor = 1.0;
        if (!atom.getUse()) {
          useFactor = 0.0;
        }
        // Scale all implicit solvent terms with the square of electrostatics lambda
        // (so dUcav / dL is 0 at lambdaElec = 0).
        double lambdaScale = lambdaElec * lambdaElec;
        if (!atom.applyLambda()) {
          lambdaScale = 1.0;
        }
        useFactor *= lambdaScale;

        double radius = rad[index];
        int isHydrogen = OpenMM_False;
        if (atom.isHydrogen()) {
          isHydrogen = OpenMM_True;
          radius = 0.0;
        }

        OpenMM_AmoebaGKCavitationForce_setParticleParameters(
            amoebaCavitationForce,
            index,
            radius * OpenMM_NmPerAngstrom,
            surfaceTension * useFactor,
            isHydrogen);
      }

      if (context.contextPointer != null) {
        OpenMM_AmoebaGKCavitationForce_updateParametersInContext(amoebaCavitationForce,
            context.contextPointer);
      }
    }

    /**
     * Adds harmonic restraints (CoordRestraint objects) to OpenMM as a custom external force.
     *
     * <p>TODO: Make robust to flat-bottom restraints.
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

        PointerByReference theRestraint = OpenMM_CustomExternalForce_create(
            "k*periodicdistance(x,y,z,x0,y0,z0)^2");
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
        logger.log(
            Level.INFO, format("  Harmonic restraint force \t%6d\t%d", nRestraints, forceGroup));
      }
    }

    /** Adds restraint bonds, if any. */
    private void addRestraintBondForce() {
      List<RestraintBond> restraintBonds = getRestraintBonds();
      if (restraintBonds == null || restraintBonds.isEmpty()) {
        return;
      }

      int forceGroup = forceField.getInteger("BOND_RESTRAINT_FORCE_GROUP", 0);

      // OpenMM's HarmonicBondForce class uses k, not 1/2*k as does FFX.
      double kParameterConversion =
          BondType.units * 2.0 * OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);

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
              OpenMM_CustomBondForce_addGlobalParameter(
                  theForce, "cubic", BondType.cubic / OpenMM_NmPerAngstrom);
              OpenMM_CustomBondForce_addGlobalParameter(
                  theForce,
                  "quartic",
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

      logger.log(
          Level.INFO,
          format("  Restraint bonds force \t%6d\t%d", restraintBonds.size(), forceGroup));
    }

    private void addRestrainGroupsForce() {
      RestrainGroups restrainGroups = getRestrainGroups();
      if (restrainGroups == null) {
        return;
      }

      // In the expression below, u and l are the upper and lower threshold
      PointerByReference force = OpenMM_CustomCentroidBondForce_create(2,
          "step(distance(g1,g2)-u)*k*(distance(g1,g2)-u)^2+step(l-distance(g1,g2))*k*(distance(g1,g2)-l)^2");
      OpenMM_CustomCentroidBondForce_addPerBondParameter(force, "k");
      OpenMM_CustomCentroidBondForce_addPerBondParameter(force, "l");
      OpenMM_CustomCentroidBondForce_addPerBondParameter(force, "u");

      // Create the Restrain Groups.
      int nGroups = restrainGroups.getNumberOfGroups();
      for (int j = 0; j < nGroups; j++) {
        PointerByReference group = OpenMM_IntArray_create(0);
        PointerByReference weight = OpenMM_DoubleArray_create(0);
        int[] groupMembers = restrainGroups.getGroupMembers(j);
        for (int i : groupMembers) {
          OpenMM_IntArray_append(group, i);
          OpenMM_DoubleArray_append(weight, atoms[i].getMass());
        }
        OpenMM_CustomCentroidBondForce_addGroup(force, group, weight);
        OpenMM_IntArray_destroy(group);
        OpenMM_DoubleArray_destroy(weight);
      }

      // Add the restraints between groups.
      double convert = OpenMM_KJPerKcal / (OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom);
      int nRestraints = restrainGroups.getNumberOfRestraints();
      int[] group1 = restrainGroups.getGroup1();
      int[] group2 = restrainGroups.getGroup2();
      double[] forceConstants = restrainGroups.getForceConstants();
      double[] smallerDistance = restrainGroups.getSmallerDistance();
      double[] largerDistance = restrainGroups.getLargerDistance();
      for (int i = 0; i < nRestraints; i++) {
        PointerByReference group = OpenMM_IntArray_create(0);
        OpenMM_IntArray_append(group, group1[i]);
        OpenMM_IntArray_append(group, group2[i]);
        PointerByReference params = OpenMM_DoubleArray_create(0);
        OpenMM_DoubleArray_append(params, forceConstants[i] * convert);
        OpenMM_DoubleArray_append(params, smallerDistance[i] * OpenMM_NmPerAngstrom);
        OpenMM_DoubleArray_append(params, largerDistance[i] * OpenMM_NmPerAngstrom);
        OpenMM_CustomCentroidBondForce_addBond(force, group, params);
        OpenMM_IntArray_destroy(group);
        OpenMM_DoubleArray_destroy(params);
      }

      // Add the constraint force.
      int forceGroup = forceField.getInteger("BOND_RESTRAINT_FORCE_GROUP", 0);
      OpenMM_Force_setForceGroup(force, forceGroup);

      if (getCrystal().aperiodic()) {
        OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(force, OpenMM_False);
      } else {
        OpenMM_CustomCentroidBondForce_setUsesPeriodicBoundaryConditions(force, OpenMM_True);
      }
      OpenMM_System_addForce(system, force);
      logger.log(Level.INFO, format("  Restrain Groups \t%6d\t\t%1d", nRestraints, forceGroup));
    }

    /** Add a constraint to every bond. */
    private void addUpBondConstraints() {
      Bond[] bonds = getBonds();
      if (bonds == null || bonds.length < 1) {
        return;
      }

      logger.info(" Adding constraints for all bonds.");
      for (Bond bond : bonds) {
        Atom atom1 = bond.getAtom(0);
        Atom atom2 = bond.getAtom(1);
        int iAtom1 = atom1.getXyzIndex() - 1;
        int iAtom2 = atom2.getXyzIndex() - 1;
        OpenMM_System_addConstraint(
            system, iAtom1, iAtom2, bond.bondType.distance * OpenMM_NmPerAngstrom);
      }
    }

    /** Add a constraint to every bond that includes a hydrogen atom. */
    private void addHydrogenConstraints() {
      Bond[] bonds = getBonds();
      if (bonds == null || bonds.length < 1) {
        return;
      }

      logger.info(" Adding constraints for hydrogen bonds.");
      for (Bond bond : bonds) {
        Atom atom1 = bond.getAtom(0);
        Atom atom2 = bond.getAtom(1);
        if (atom1.isHydrogen() || atom2.isHydrogen()) {
          BondType bondType = bond.bondType;
          int iAtom1 = atom1.getXyzIndex() - 1;
          int iAtom2 = atom2.getXyzIndex() - 1;
          OpenMM_System_addConstraint(
              system, iAtom1, iAtom2, bondType.distance * OpenMM_NmPerAngstrom);
        }
      }
    }

    /** Add a constraint to every angle that includes two hydrogen atoms. */
    private void setUpHydrogenAngleConstraints() {
      Angle[] angles = getAngles();
      if (angles == null || angles.length < 1) {
        return;
      }

      logger.info(" Adding hydrogen angle constraints.");

      for (Angle angle : angles) {
        if (isHydrogenAngle(angle)) {
          Atom atom1 = angle.getAtom(0);
          Atom atom3 = angle.getAtom(2);

          // Calculate a "false bond" length between atoms 1 and 3 to constrain the angle using the
          // law of cosines.
          Bond bond1 = angle.getBond(0);
          double distance1 = bond1.bondType.distance;

          Bond bond2 = angle.getBond(1);
          double distance2 = bond2.bondType.distance;

          // Equilibrium angle value in degrees.
          double angleVal = angle.angleType.angle[angle.nh];

          // Law of cosines.
          double falseBondLength =
              sqrt(
                  distance1 * distance1
                      + distance2 * distance2
                      - 2.0 * distance1 * distance2 * cos(toRadians(angleVal)));

          int iAtom1 = atom1.getXyzIndex() - 1;
          int iAtom3 = atom3.getXyzIndex() - 1;
          OpenMM_System_addConstraint(
              system, iAtom1, iAtom3, falseBondLength * OpenMM_NmPerAngstrom);
        }
      }
    }

    /**
     * Check to see if an angle is a hydrogen angle. This method only returns true for hydrogen
     * angles that are less than 160 degrees.
     *
     * @param angle Angle to check.
     * @return boolean indicating whether or not an angle is a hydrogen angle that is less than 160
     *     degrees.
     */
    private boolean isHydrogenAngle(Angle angle) {
      if (angle.containsHydrogen()) {
        // Equilibrium angle value in degrees.
        double angleVal = angle.angleType.angle[angle.nh];
        // Make sure angle is less than 160 degrees because greater than 160 degrees will not be
        // constrained
        // well using the law of cosines.
        if (angleVal < 160.0) {
          Atom atom1 = angle.getAtom(0);
          Atom atom2 = angle.getAtom(1);
          Atom atom3 = angle.getAtom(2);
          // Setting constraints only on angles where atom1 or atom3 is a hydrogen while atom2 is
          // not a hydrogen.
          return atom1.isHydrogen() && atom3.isHydrogen() && !atom2.isHydrogen();
        }
      }
      return false;
    }

  }

  /** Retrieve state information from an OpenMM Simulation. */
  public class State {

    /** Potential energy (kcal/mol). */
    public final double potentialEnergy;
    /** Kinetic energy (kcal/mol). */
    public final double kineticEnergy;
    /** Total energy (kcal/mol). */
    public final double totalEnergy;
    /** Temperature (K). */
    public final double temperature;
    /** Pointer to an OpenMM state. */
    private final PointerByReference state;
    /** If true, positions are available. */
    private final boolean positions;
    /** If true, energies are available. */
    private final boolean energies;
    /** If true, velocities are available. */
    private final boolean velocities;
    /** If true. forces are available. */
    private final boolean forces;

    /**
     * Construct an OpenMM State with the requested information.
     *
     * @param positions Retrieve positions.
     * @param energies Retrieve energies.
     * @param forces Retrieve forces.
     * @param velocities Retrieve velocities.
     */
    public State(boolean positions, boolean energies, boolean forces, boolean velocities) {
      this.positions = positions;
      this.energies = energies;
      this.forces = forces;
      this.velocities = velocities;

      // Construct the mask.
      int mask = 0;
      if (positions) {
        mask = OpenMM_State_Positions;
      }
      if (energies) {
        mask += OpenMM_State_Energy;
      }
      if (velocities) {
        mask += OpenMM_State_Velocities;
      }
      if (forces) {
        mask += OpenMM_State_Forces;
      }

      // Retrieve the state.
      state = OpenMM_Context_getState(context.getContextPointer(), mask, enforcePBC);

      if (energies) {
        potentialEnergy = OpenMM_State_getPotentialEnergy(state) * OpenMM_KcalPerKJ;
        kineticEnergy = OpenMM_State_getKineticEnergy(state) * OpenMM_KcalPerKJ;
        totalEnergy = potentialEnergy + kineticEnergy;
        double dof = system.calculateDegreesOfFreedom();
        temperature = 2.0 * kineticEnergy * KCAL_TO_GRAM_ANG2_PER_PS2 / (kB * dof);
      } else {
        potentialEnergy = 0.0;
        kineticEnergy = 0.0;
        totalEnergy = 0.0;
        temperature = 0.0;
      }
    }

    public void free() {
      if (state != null && freeOpenMM) {
        logger.fine(" Free OpenMM State.");
        OpenMM_State_destroy(state);
        logger.fine(" Free OpenMM State completed.");
      }
    }

    /**
     * The force array contains the OpenMM force information for all atoms. The returned array a
     * contains accelerations for active atoms only.
     *
     * @param a Acceleration components for only active atomic coordinates.
     * @return a The acceleration for each active atomic coordinate.
     */
    public double[] getAccelerations(double[] a) {
      if (!forces) {
        return a;
      }

      int n = getNumberOfVariables();
      if (a == null || a.length != n) {
        a = new double[n];
      }

      PointerByReference forcePointer = OpenMM_State_getForces(state);

      for (int i = 0, index = 0; i < nAtoms; i++) {
        Atom atom = atoms[i];
        if (atom.isActive()) {
          double mass = atom.getMass();
          OpenMM_Vec3 acc = OpenMM_Vec3Array_get(forcePointer, i);
          double xx = acc.x * OpenMM_AngstromsPerNm / mass;
          double yy = acc.y * OpenMM_AngstromsPerNm / mass;
          double zz = acc.z * OpenMM_AngstromsPerNm / mass;
          a[index++] = xx;
          a[index++] = yy;
          a[index++] = zz;
          atom.setAcceleration(xx, yy, zz);
        }
      }
      return a;
    }

    /**
     * The force array contains the OpenMM force information for all atoms. The returned array g
     * contains components for active atoms only.
     *
     * @param g Gradient components for only active atomic coordinates.
     * @return g The gradient includes only active atoms
     */
    public double[] getGradient(double[] g) {
      if (!forces) {
        return g;
      }

      int n = getNumberOfVariables();
      if (g == null || g.length != n) {
        g = new double[n];
      }

      PointerByReference forcePointer = OpenMM_State_getForces(state);
      for (int i = 0, index = 0; i < nAtoms; i++) {
        Atom atom = atoms[i];
        if (atom.isActive()) {
          OpenMM_Vec3 force = OpenMM_Vec3Array_get(forcePointer, i);
          double xx = -force.x * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
          double yy = -force.y * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
          double zz = -force.z * OpenMM_NmPerAngstrom * OpenMM_KcalPerKJ;
          if (isNaN(xx)
              || isInfinite(xx)
              || isNaN(yy)
              || isInfinite(yy)
              || isNaN(zz)
              || isInfinite(zz)) {
            StringBuilder sb =
                new StringBuilder(
                    format(
                        " The gradient of atom %s is (%8.3f,%8.3f,%8.3f).",
                        atom.toString(), xx, yy, zz));
            double[] vals = new double[3];
            atom.getVelocity(vals);
            sb.append(format("\n Velocities: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
            atom.getAcceleration(vals);
            sb.append(format("\n Accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
            atom.getPreviousAcceleration(vals);
            sb.append(
                format("\n Previous accelerations: %8.3g %8.3g %8.3g", vals[0], vals[1], vals[2]));
            throw new EnergyException(sb.toString());
          }
          g[index++] = xx;
          g[index++] = yy;
          g[index++] = zz;
          atom.setXYZGradient(xx, yy, zz);
        }
      }
      return g;
    }

    /**
     * Read the periodic lattice vectors from a state.
     *
     * <p>The crystal instance will be updated, and passed to the ForceFieldEnergy instance.
     */
    public void getPeriodicBoxVectors() {
      if (!positions) {
        return;
      }

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
        setCrystal(crystal);
      }
    }

    /**
     * The positions array contains the OpenMM atomic position information for all atoms. The
     * returned array x contains coordinates only for active atoms.
     *
     * @param x Atomic coordinates only for active atoms.
     * @return x The atomic coordinates for only active atoms.
     */
    public double[] getPositions(double[] x) {
      if (!positions) {
        return x;
      }

      int n = getNumberOfVariables();
      if (x == null || x.length != n) {
        x = new double[n];
      }

      PointerByReference positionsPointer = OpenMM_State_getPositions(state);
      for (int i = 0, index = 0; i < nAtoms; i++) {
        Atom atom = atoms[i];
        if (atom.isActive()) {
          OpenMM_Vec3 pos = OpenMM_Vec3Array_get(positionsPointer, i);
          double xx = pos.x * OpenMM_AngstromsPerNm;
          double yy = pos.y * OpenMM_AngstromsPerNm;
          double zz = pos.z * OpenMM_AngstromsPerNm;
          x[index++] = xx;
          x[index++] = yy;
          x[index++] = zz;
          atom.moveTo(xx, yy, zz);
        }
      }
      return x;
    }

    /**
     * The positions array contains the OpenMM atomic position information for all atoms. The
     * returned array x contains coordinates for active atoms only.
     *
     * @param v Velocity only for active atomic coordinates.
     * @return v The velocity for each active atomic coordinate.
     */
    public double[] getVelocities(double[] v) {
      if (!velocities) {
        return v;
      }

      int n = getNumberOfVariables();
      if (v == null || v.length != n) {
        v = new double[n];
      }

      PointerByReference velocitiesPointer = OpenMM_State_getVelocities(state);
      for (int i = 0, index = 0; i < nAtoms; i++) {
        Atom atom = atoms[i];
        if (atom.isActive()) {
          OpenMM_Vec3 vel = OpenMM_Vec3Array_get(velocitiesPointer, i);
          double xx = vel.x * OpenMM_AngstromsPerNm;
          double yy = vel.y * OpenMM_AngstromsPerNm;
          double zz = vel.z * OpenMM_AngstromsPerNm;
          v[index++] = xx;
          v[index++] = yy;
          v[index++] = zz;
          atom.setVelocity(xx, yy, zz);
        }
      }
      return v;
    }
  }
}
