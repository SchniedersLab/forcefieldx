// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2023.
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
package ffx.potential.openmm;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.pme.Polarization;
import ffx.potential.nonbonded.pme.SCFAlgorithm;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;

import java.util.logging.Level;
import java.util.logging.Logger;

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
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

public class AmoebaMultipoleForce extends OpenMMForce {

  private static final Logger logger = Logger.getLogger(AmoebaMultipoleForce.class.getName());

  public AmoebaMultipoleForce(OpenMMEnergy openMMEnergy) {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    if (pme == null) {
      return;
    }

    int[][] axisAtom = pme.getAxisAtoms();
    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    pointer = OpenMM_AmoebaMultipoleForce_create();

    double polarScale = 1.0;
    SCFAlgorithm scfAlgorithm = null;

    ForceField forceField = openMMEnergy.getMolecularAssembly().getForceField();

    if (pme.getPolarizationType() != Polarization.MUTUAL) {
      setPolarizationType(OpenMM_AmoebaMultipoleForce_Direct);
      if (pme.getPolarizationType() == Polarization.NONE) {
        polarScale = 0.0;
      }
    } else {
      String algorithm = forceField.getString("SCF_ALGORITHM", "CG");
      try {
        algorithm = algorithm.replaceAll("-", "_").toUpperCase();
        scfAlgorithm = SCFAlgorithm.valueOf(algorithm);
      } catch (Exception e) {
        scfAlgorithm = SCFAlgorithm.CG;
      }

      if (scfAlgorithm == SCFAlgorithm.EPT) {
        /*
         * Citation:
         * Simmonett, A. C.;  Pickard, F. C. t.;  Shao, Y.;  Cheatham, T. E., 3rd; Brooks, B. R.,
         * Efficient treatment of induced dipoles. The Journal of chemical physics 2015, 143 (7), 074115-074115.
         */
        setPolarizationType(OpenMM_AmoebaMultipoleForce_Extrapolated);
        OpenMMDoubleArray exptCoefficients = new OpenMMDoubleArray(4);
        exptCoefficients.set(0, -0.154);
        exptCoefficients.set(1, 0.017);
        exptCoefficients.set(2, 0.657);
        exptCoefficients.set(3, 0.475);
        setExtrapolationCoefficients(exptCoefficients);
        exptCoefficients.destroy();
      } else {
        setPolarizationType(OpenMM_AmoebaMultipoleForce_Mutual);
      }
    }

    OpenMMDoubleArray dipoles = new OpenMMDoubleArray(3);
    OpenMMDoubleArray quadrupoles = new OpenMMDoubleArray(9);

    OpenMMSystem system = openMMEnergy.getSystem();
    double lambdaElec = system.getLambdaElec();

    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    int nAtoms = atoms.length;
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = atoms[i];
      MultipoleType multipoleType = pme.getMultipoleType(i);
      PolarizeType polarType = pme.getPolarizeType(i);

      // Define the frame definition.
      int axisType = switch (multipoleType.frameDefinition) {
        case NONE -> OpenMM_AmoebaMultipoleForce_NoAxisType;
        case ZONLY -> OpenMM_AmoebaMultipoleForce_ZOnly;
        case ZTHENX -> OpenMM_AmoebaMultipoleForce_ZThenX;
        case BISECTOR -> OpenMM_AmoebaMultipoleForce_Bisector;
        case ZTHENBISECTOR -> OpenMM_AmoebaMultipoleForce_ZBisect;
        case THREEFOLD -> OpenMM_AmoebaMultipoleForce_ThreeFold;
      };

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
        dipoles.set(j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
      }
      int l = 0;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion * useFactor / 3.0);
        }
      }

      int zaxis = -1;
      int xaxis = -1;
      int yaxis = -1;
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
      addMultipole(charge, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis, polarType.thole,
          polarType.pdamp * dampingFactorConversion, polarType.polarizability * polarityConversion * polarScale);
    }
    dipoles.destroy();
    quadrupoles.destroy();

    Crystal crystal = openMMEnergy.getCrystal();
    if (!crystal.aperiodic()) {
      setNonbondedMethod(OpenMM_AmoebaMultipoleForce_PME);
      setCutoffDistance(pme.getEwaldCutoff() * OpenMM_NmPerAngstrom);
      setAEwald(pme.getEwaldCoefficient() / OpenMM_NmPerAngstrom);

      double ewaldTolerance = 1.0e-04;
      setEwaldErrorTolerance(ewaldTolerance);

      OpenMMIntArray gridDimensions = new OpenMMIntArray(3);
      ReciprocalSpace recip = pme.getReciprocalSpace();
      gridDimensions.set(0, recip.getXDim());
      gridDimensions.set(1, recip.getYDim());
      gridDimensions.set(2, recip.getZDim());
      setPmeGridDimensions(gridDimensions);
      gridDimensions.destroy();
    } else {
      setNonbondedMethod(OpenMM_AmoebaMultipoleForce_NoCutoff);
    }

    setMutualInducedMaxIterations(500);
    setMutualInducedTargetEpsilon(pme.getPolarEps());

    int[][] ip11 = pme.getPolarization11();

    OpenMMIntArray covalentMap = new OpenMMIntArray(0);
    for (int i = 0; i < nAtoms; i++) {
      Atom ai = atoms[i];

      // 1-2 Mask
      covalentMap.resize(0);
      for (Atom ak : ai.get12List()) {
        covalentMap.append(ak.getIndex() - 1);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent12, covalentMap);

      // 1-3 Mask
      covalentMap.resize(0);
      for (Atom ak : ai.get13List()) {
        covalentMap.append(ak.getIndex() - 1);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent13, covalentMap);

      // 1-4 Mask
      covalentMap.resize(0);
      for (Atom ak : ai.get14List()) {
        covalentMap.append(ak.getIndex() - 1);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent14, covalentMap);

      // 1-5 Mask
      covalentMap.resize(0);
      for (Atom ak : ai.get15List()) {
        covalentMap.append(ak.getIndex() - 1);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent15, covalentMap);

      // 1-1 Polarization Groups.
      covalentMap.resize(0);
      for (int j = 0; j < ip11[i].length; j++) {
        covalentMap.append(ip11[i][j]);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_PolarizationCovalent11, covalentMap);

      // AMOEBA does not scale between 1-2, 1-3, etc. polarization groups.
    }
    covalentMap.destroy();

    int forceGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
    setForceGroup(forceGroup);
    logger.log(Level.INFO, format("  AMOEBA polarizable multipole force \t%d", forceGroup));

    if (scfAlgorithm == SCFAlgorithm.EPT) {
      logger.info("   Using extrapolated perturbation theory for polarization energy.");
    }
  }

  /**
   * Convenience method to construct an AMOEBA Multipole Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the multipole information.
   * @return An AMOEBA Multipole Force, or null if there are no multipole interactions.
   */
  public static OpenMMForce constructForce(OpenMMEnergy openMMEnergy) {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    if (pme == null) {
      return null;
    }
    return new AmoebaMultipoleForce(openMMEnergy);
  }

  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    double polarScale = 1.0;
    if (pme.getPolarizationType() == Polarization.NONE) {
      polarScale = 0.0;
    }

    OpenMMDoubleArray dipoles = new OpenMMDoubleArray(3);
    OpenMMDoubleArray quadrupoles = new OpenMMDoubleArray(9);

    double lambdaElec = openMMEnergy.getSystem().getLambdaElec();

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
      int axisType = switch (multipoleType.frameDefinition) {
        case NONE -> OpenMM_AmoebaMultipoleForce_NoAxisType;
        case ZONLY -> OpenMM_AmoebaMultipoleForce_ZOnly;
        case ZTHENX -> OpenMM_AmoebaMultipoleForce_ZThenX;
        case BISECTOR -> OpenMM_AmoebaMultipoleForce_Bisector;
        case ZTHENBISECTOR -> OpenMM_AmoebaMultipoleForce_ZBisect;
        case THREEFOLD -> OpenMM_AmoebaMultipoleForce_ThreeFold;
      };

      // Load local multipole coefficients.
      for (int j = 0; j < 3; j++) {
        dipoles.set(j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * useFactor);
      }
      int l = 0;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0 * useFactor);
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
      setMultipoleParameters(index, multipoleType.charge * useFactor,
          dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
          polarizeType.thole, polarizeType.pdamp * dampingFactorConversion,
          polarizeType.polarizability * polarityConversion * polarScale * useFactor);
    }

    dipoles.destroy();
    quadrupoles.destroy();
    updateParametersInContext(openMMEnergy.getContext());
  }


  /**
   * Set the polarization method.
   *
   * @param method The polarization method.
   */
  public void setPolarizationType(int method) {
    OpenMM_AmoebaMultipoleForce_setPolarizationType(pointer, method);
  }

  /**
   * Set extrapolation coefficients.
   *
   * @param exptCoefficients The extrapolation coefficients.
   */
  public void setExtrapolationCoefficients(OpenMMDoubleArray exptCoefficients) {
    OpenMM_AmoebaMultipoleForce_setExtrapolationCoefficients(pointer, exptCoefficients.getPointer());
  }

  /**
   * Add a multipole.
   *
   * @param charge         The charge.
   * @param dipole         The dipole.
   * @param quadrupole     The quadrupole.
   * @param axisType       The axis type.
   * @param zaxis          The z-axis.
   * @param xaxis          The x-axis.
   * @param yaxis          The y-axis.
   * @param thole          The Thole parameter.
   * @param pdamp          The damping factor.
   * @param polarizability The polarizability.
   */
  public void addMultipole(double charge, OpenMMDoubleArray dipole, OpenMMDoubleArray quadrupole, int axisType,
                           int zaxis, int xaxis, int yaxis, double thole, double pdamp, double polarizability) {
    OpenMM_AmoebaMultipoleForce_addMultipole(pointer, charge, dipole.getPointer(), quadrupole.getPointer(),
        axisType, zaxis, xaxis, yaxis, thole, pdamp, polarizability);
  }

  /**
   * Set the multipole parameters.
   *
   * @param index          The atom index.
   * @param charge         The charge.
   * @param dipoles        The dipole.
   * @param quadrupoles    The quadrupole.
   * @param axisType       The axis type.
   * @param zaxis          The z-axis.
   * @param xaxis          The x-axis.
   * @param yaxis          The y-axis.
   * @param thole          The Thole parameter.
   * @param pdamp          The damping factor.
   * @param polarizability The polarizability.
   */
  public void setMultipoleParameters(int index, double charge, OpenMMDoubleArray dipoles, OpenMMDoubleArray quadrupoles,
                                     int axisType, int zaxis, int xaxis, int yaxis,
                                     double thole, double pdamp, double polarizability) {
    OpenMM_AmoebaMultipoleForce_setMultipoleParameters(pointer, index, charge,
        dipoles.getPointer(), quadrupoles.getPointer(), axisType, zaxis, xaxis, yaxis, thole, pdamp, polarizability);
  }

  /**
   * Set the nonbonded method for the multipole force.
   *
   * @param method The nonbonded method.
   */
  public void setNonbondedMethod(int method) {
    OpenMM_AmoebaMultipoleForce_setNonbondedMethod(pointer, method);
  }

  /**
   * Set the cutoff distance for the multipole force.
   *
   * @param cutoff The cutoff distance.
   */
  public void setCutoffDistance(double cutoff) {
    OpenMM_AmoebaMultipoleForce_setCutoffDistance(pointer, cutoff);
  }

  /**
   * Set the Ewald coefficient for the multipole force.
   *
   * @param aewald The Ewald coefficient.
   */
  public void setAEwald(double aewald) {
    OpenMM_AmoebaMultipoleForce_setAEwald(pointer, aewald);
  }

  /**
   * Set the Ewald error tolerance for the multipole force.
   *
   * @param ewaldTolerance The Ewald error tolerance.
   */
  public void setEwaldErrorTolerance(double ewaldTolerance) {
    OpenMM_AmoebaMultipoleForce_setEwaldErrorTolerance(pointer, ewaldTolerance);
  }

  /**
   * Set the PME grid dimensions for the multipole force.
   *
   * @param gridDimensions The PME grid dimensions.
   */
  public void setPmeGridDimensions(OpenMMIntArray gridDimensions) {
    OpenMM_AmoebaMultipoleForce_setPmeGridDimensions(pointer, gridDimensions.getPointer());
  }

  //    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(amoebaMultipoleForce, pme.getPolarEps());

  /**
   * Set the mutual induced target maximum number of iterations.
   *
   * @param iterations The mutual induced max iterations.
   */
  public void setMutualInducedMaxIterations(int iterations) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedMaxIterations(pointer, iterations);
  }

  /**
   * Set the mutual induced target epsilon.
   *
   * @param epsilon The mutual induced target epsilon.
   */
  public void setMutualInducedTargetEpsilon(double epsilon) {
    OpenMM_AmoebaMultipoleForce_setMutualInducedTargetEpsilon(pointer, epsilon);
  }

  /**
   * Set the covalent map.
   *
   * @param i            The atom index.
   * @param covalentType The covalent type.
   * @param covalentMap  The covalent map.
   */
  public void setCovalentMap(int i, int covalentType, OpenMMIntArray covalentMap) {
    OpenMM_AmoebaMultipoleForce_setCovalentMap(pointer, i, covalentType, covalentMap.getPointer());
  }

  /**
   * Update the parameters in the context.
   *
   * @param openMMContext The OpenMM context.
   */
  public void updateParametersInContext(OpenMMContext openMMContext) {
    if (openMMContext.hasContextPointer()) {
      OpenMM_AmoebaMultipoleForce_updateParametersInContext(pointer, openMMContext.getPointer());
    }
  }

}
