// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2026.
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
import ffx.openmm.DoubleArray;
import ffx.openmm.Force;
import ffx.openmm.IntArray;
import ffx.openmm.amoeba.MultipoleForce;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.ReciprocalSpace;
import ffx.potential.nonbonded.pme.AlchemicalParameters;
import ffx.potential.nonbonded.pme.Polarization;
import ffx.potential.nonbonded.pme.SCFAlgorithm;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.PolarizeType;

import java.util.HashSet;
import java.util.Set;
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
import static edu.uiowa.jopenmm.OpenMMAmoebaLibrary.OpenMM_NmPerAngstrom;
import static java.lang.String.format;
import static org.apache.commons.math3.util.FastMath.sqrt;

/**
 * AmoebaMultipoleForce.
 */
public class AmoebaMultipoleForce extends MultipoleForce {

  private static final Logger logger = Logger.getLogger(AmoebaMultipoleForce.class.getName());

  /**
   * Construct an AMOEBA Multipole Force.
   *
   * @param openMMEnergy The OpenMMEnergy instance that contains the multipole information.
   */
  public AmoebaMultipoleForce(OpenMMEnergy openMMEnergy) {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    if (pme == null) {
      destroy();
      return;
    }

    double doPolarization = configureForce(openMMEnergy);

    int[][] axisAtom = pme.getAxisAtoms();
    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    boolean lambdaTerm = pme.getLambdaTerm();
    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    double permLambda = alchemicalParameters.permLambda;
    double polarLambda = alchemicalParameters.polLambda;
    DoubleArray dipoles = new DoubleArray(3);
    DoubleArray quadrupoles = new DoubleArray(9);

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
      if (!atom.getUse() || !atom.getElectrostatics()) {
        useFactor = 0.0;
      }

      double permScale = useFactor;
      double polarScale = doPolarization * useFactor;
      if (lambdaTerm && atom.applyLambda()) {
        permScale *= permLambda;
        polarScale *= polarLambda;
      }

      // Load local multipole coefficients.
      for (int j = 0; j < 3; j++) {
        dipoles.set(j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * permScale);
      }
      int l = 0;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion * permScale / 3.0);
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

      double charge = multipoleType.charge * permScale;

      // Add the multipole.
      addMultipole(charge, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis, polarType.thole,
          polarType.pdamp * dampingFactorConversion, polarType.polarizability * polarityConversion * polarScale);
    }
    dipoles.destroy();
    quadrupoles.destroy();

    int[][] ip11 = pme.getPolarization11();
    IntArray covalentMap = new IntArray(0);
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
  }

  /**
   * Construct an AMOEBA Multipole Force.
   *
   * @param topology                 The topology index for the dual topology system.
   * @param openMMDualTopologyEnergy The OpenMM Dual-Topology Energy instance that contains the multipole parameters.
   */
  public AmoebaMultipoleForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    // Determine the other topology index.
    int otherTopology = 1 - topology;

    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ForceFieldEnergy otherForceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(otherTopology);

    ParticleMeshEwald pme = forceFieldEnergy.getPmeNode();
    ParticleMeshEwald otherPME = otherForceFieldEnergy.getPmeNode();
    if (pme == null || otherPME == null) {
      destroy();
      return;
    }

    double doPolarization = configureForce(forceFieldEnergy);

    // Dual-topology scale factor.
    // double scaleDT = Math.sqrt(openMMDualTopologyEnergy.getTopologyScale(topology));

    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    boolean lambdaTerm = pme.getLambdaTerm();
    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    double permLambda = alchemicalParameters.permLambda;
    double polarLambda = alchemicalParameters.polLambda;
    DoubleArray dipoles = new DoubleArray(3);
    DoubleArray quadrupoles = new DoubleArray(9);

    int nAtoms = openMMDualTopologyEnergy.getNumberOfAtoms();

    // Add a particle for each atom in the dual topology.
    for (int i = 0; i < nAtoms; i++) {
      Atom atom = openMMDualTopologyEnergy.getDualTopologyAtom(topology, i);
      if (atom.getTopologyIndex() == topology) {
        int index = atom.getArrayIndex();
        MultipoleType multipoleType = pme.getMultipoleType(index);
        PolarizeType polarType = pme.getPolarizeType(index);

        // Define the frame definition.
        int axisType = switch (multipoleType.frameDefinition) {
          case NONE -> OpenMM_AmoebaMultipoleForce_NoAxisType;
          case ZONLY -> OpenMM_AmoebaMultipoleForce_ZOnly;
          case ZTHENX -> OpenMM_AmoebaMultipoleForce_ZThenX;
          case BISECTOR -> OpenMM_AmoebaMultipoleForce_Bisector;
          case ZTHENBISECTOR -> OpenMM_AmoebaMultipoleForce_ZBisect;
          case THREEFOLD -> OpenMM_AmoebaMultipoleForce_ThreeFold;
        };

        // All multipoles are scaled by the dual topology scale factor.
        // double useFactor = scaleDT;
        double useFactor = 1.0;

        // Check if the atom is used and has electrostatics.
        if (!atom.getUse() || !atom.getElectrostatics()) {
          // This atom is not used or does not have electrostatics.
          useFactor = 0.0;
        }

        // Further scale the multipole coefficients for alchemical atoms.
        double permScale = useFactor;
        double polarScale = doPolarization * useFactor;
        if (lambdaTerm && atom.applyLambda()) {
          permScale *= permLambda;
          polarScale *= polarLambda;
        }

        // Load local multipole coefficients.
        double charge = multipoleType.charge * permScale;
        for (int j = 0; j < 3; j++) {
          dipoles.set(j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * permScale);
        }
        int l = 0;
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion * permScale / 3.0);
          }
        }

        int zaxis = -1;
        int xaxis = -1;
        int yaxis = -1;
        int[] refAtoms = atom.getAxisAtomIndices();
        if (refAtoms != null) {
          zaxis = refAtoms[0];
          zaxis = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, zaxis);
          if (refAtoms.length > 1) {
            xaxis = refAtoms[1];
            xaxis = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, xaxis);
            if (refAtoms.length > 2) {
              yaxis = refAtoms[2];
              yaxis = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, yaxis);
            }
          }
        } else {
          axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
        }

        // Add the multipole.
        addMultipole(charge, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
            polarType.thole, polarType.pdamp * dampingFactorConversion,
            polarType.polarizability * polarityConversion * polarScale);
      } else {
        // Add a fake multipole.
        // Define the frame definition.
        int axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
        // Load zero multipole coefficients.
        double charge = 0.0;
        for (int j = 0; j < 3; j++) {
          dipoles.set(j, 0.0);
        }
        int l = 0;
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            quadrupoles.set(l++, 0.0);
          }
        }
        // No frame defining atoms.
        int zaxis = -1;
        int xaxis = -1;
        int yaxis = -1;
        // Polarization parameters.
        double thole = 0.39;
        double pdamp = 1.0;
        // Add the fake site.
        addMultipole(charge, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
            thole, pdamp, 0.0);
      }
    }
    dipoles.destroy();
    quadrupoles.destroy();

    int[][] ip11 = pme.getPolarization11();
    int[][] ip11Other = otherPME.getPolarization11();

    IntArray covalentMap = new IntArray(0);
    // Use a set to ensure the same masking rules for both topologies.
    Set<Integer> covalentSet = new HashSet<>();

    for (int i = 0; i < nAtoms; i++) {
      Atom atom = openMMDualTopologyEnergy.getDualTopologyAtom(topology, i);
      Atom otherAtom = openMMDualTopologyEnergy.getDualTopologyAtom(otherTopology, i);
      int index = atom.getArrayIndex();
      int otherIndex = otherAtom.getArrayIndex();

      // 1-2 Mask
      covalentSet.clear();
      for (Atom ak : atom.get12List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      for (Atom ak : otherAtom.get12List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      covalentMap.resize(0);
      for (int ak : covalentSet) {
        covalentMap.append(ak);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent12, covalentMap);

      // 1-3 Mask
      covalentSet.clear();
      for (Atom ak : atom.get13List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      for (Atom ak : otherAtom.get13List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      covalentMap.resize(0);
      for (int ak : covalentSet) {
        covalentMap.append(ak);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent13, covalentMap);

      // 1-4 Mask
      covalentSet.clear();
      for (Atom ak : atom.get14List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      for (Atom ak : otherAtom.get14List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      covalentMap.resize(0);
      for (int ak : covalentSet) {
        covalentMap.append(ak);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent14, covalentMap);

      // 1-5 Mask
      covalentSet.clear();
      for (Atom ak : atom.get15List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      for (Atom ak : otherAtom.get15List()) {
        covalentSet.add(ak.getTopologyAtomIndex());
      }
      covalentMap.resize(0);
      for (int ak : covalentSet) {
        covalentMap.append(ak);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_Covalent15, covalentMap);

      // 1-1 Polarization Groups.
      covalentSet.clear();
      for (int k : ip11[index]) {
        int value = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, k);
        covalentSet.add(value);
      }
      for (int k : ip11Other[otherIndex]) {
        int value = openMMDualTopologyEnergy.mapToDualTopologyIndex(otherTopology, k);
        covalentSet.add(value);
      }
      covalentMap.resize(0);
      for (int k : covalentSet) {
        covalentMap.append(k);
      }
      setCovalentMap(i, OpenMM_AmoebaMultipoleForce_PolarizationCovalent11, covalentMap);
      // AMOEBA does not scale between 1-2, 1-3, etc. polarization groups.
    }
    covalentMap.destroy();
  }

  /**
   * Configure the AMOEBA Multipole Force based on the OpenMM Energy instance.
   *
   * @param forceFieldEnergy The ForceFieldEnergy instance that contains the multipole information.
   * @return The polarization factor for the force, which is 1.0 if polarization is enabled, or 0.0 if not.
   */
  private double configureForce(ForceFieldEnergy forceFieldEnergy) {
    ParticleMeshEwald pme = forceFieldEnergy.getPmeNode();
    ForceField forceField = forceFieldEnergy.getMolecularAssembly().getForceField();

    Polarization polarization = pme.getPolarizationType();
    double doPolarization = 1.0;
    SCFAlgorithm scfAlgorithm = null;
    if (polarization != Polarization.MUTUAL) {
      setPolarizationType(OpenMM_AmoebaMultipoleForce_Direct);
      if (pme.getPolarizationType() == Polarization.NONE) {
        doPolarization = 0.0;
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
        DoubleArray exptCoefficients = new DoubleArray(4);
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

    Crystal crystal = forceFieldEnergy.getCrystal();
    double cutoff = pme.getEwaldCutoff();
    double aewald = pme.getEwaldCoefficient();
    if (!crystal.aperiodic()) {
      setNonbondedMethod(OpenMM_AmoebaMultipoleForce_PME);
      setCutoffDistance(cutoff * OpenMM_NmPerAngstrom);
      double ewaldTolerance = 1.0e-04;
      setEwaldErrorTolerance(ewaldTolerance);
      ReciprocalSpace recip = pme.getReciprocalSpace();
      int nx = recip.getXDim();
      int ny = recip.getYDim();
      int nz = recip.getZDim();
      setPMEParameters(aewald / OpenMM_NmPerAngstrom, nx, ny, nz);
    } else {
      setNonbondedMethod(OpenMM_AmoebaMultipoleForce_NoCutoff);
    }

    setMutualInducedMaxIterations(500);
    double poleps = pme.getPolarEps();
    setMutualInducedTargetEpsilon(poleps);

    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    boolean lambdaTerm = pme.getLambdaTerm();

    int forceGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
    setForceGroup(forceGroup);
    if (logger.isLoggable(Level.INFO)) {
      StringBuilder sb = new StringBuilder();
      sb.append("\n  Electrostatics\n");
      sb.append(format("   Polarization:                       %8s\n", polarization.toString()));
      if (polarization == Polarization.MUTUAL) {
        sb.append(format("    SCF Convergence Criteria:         %8.3e\n", poleps));
      } else if (scfAlgorithm == SCFAlgorithm.EPT) {
        sb.append(format("    SCF Algorithm:                     %8s\n", scfAlgorithm));
      }
      if (aewald > 0.0) {
        sb.append("   Particle-mesh Ewald\n");
        sb.append(format("    Ewald Coefficient:                 %8.3f\n", aewald));
        sb.append(format("    Particle Cut-Off:                  %8.3f (A)", cutoff));
      } else if (cutoff != Double.POSITIVE_INFINITY) {
        sb.append(format("    Electrostatics Cut-Off:            %8.3f (A)", cutoff));
      } else {
        sb.append("    Electrostatics Cut-Off:                NONE");
      }
      logger.info(sb.toString());
      if (lambdaTerm) {
        sb = new StringBuilder("   Alchemical Parameters\n");
        double permLambdaStart = alchemicalParameters.permLambdaStart;
        double permLambdaEnd = alchemicalParameters.permLambdaEnd;
        sb.append(format("    Permanent Multipole Range:      %5.3f-%5.3f\n", permLambdaStart, permLambdaEnd));
        double permLambdaAlpha = alchemicalParameters.permLambdaAlpha;
        if (permLambdaAlpha != 0.0) {
          logger.severe(" Permanent multipole softcore not supported for OpenMM.");
        }
        double permLambdaExponent = alchemicalParameters.permLambdaExponent;
        sb.append(format("    Permanent Multipole Lambda Exponent:  %5.3f\n", permLambdaExponent));
        if (polarization != Polarization.NONE) {
          double polLambdaExponent = alchemicalParameters.polLambdaExponent;
          double polLambdaStart = alchemicalParameters.polLambdaStart;
          double polLambdaEnd = alchemicalParameters.polLambdaEnd;
          sb.append(format("    Polarization Lambda Exponent:         %5.3f\n", polLambdaExponent));
          sb.append(format("    Polarization Range:             %5.3f-%5.3f\n", polLambdaStart, polLambdaEnd));
          boolean doNoLigandCondensedSCF = alchemicalParameters.doNoLigandCondensedSCF;
          if (doNoLigandCondensedSCF) {
            logger.severe(" Condensed SCF without a ligand is not supported for OpenMM.");
          }
        }
        boolean doLigandGKElec = alchemicalParameters.doLigandGKElec;
        if (doLigandGKElec) {
          logger.severe(" Isolated ligand electrostatics are not supported for OpenMM.");
        }
        logger.info(sb.toString());
      }
      logger.log(Level.FINE, format("   Force group:\t\t%d\n", forceGroup));
    }
    return doPolarization;
  }

  /**
   * Convenience method to construct an AMOEBA Multipole Force.
   *
   * @param openMMEnergy The OpenMM Energy instance that contains the multipole information.
   * @return An AMOEBA Multipole Force, or null if there are no multipole interactions.
   */
  public static Force constructForce(OpenMMEnergy openMMEnergy) {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    if (pme == null) {
      return null;
    }
    return new AmoebaMultipoleForce(openMMEnergy);
  }

  /**
   * Convenience method to construct an AMOEBA Multipole Force.
   *
   * @param topology                 The topology index for the dual topology system.
   * @param openMMDualTopologyEnergy The OpenMM Dual-Topology Energy instance that contains the vdW information.
   * @return An AMOEBA Multipole Force, or null if there are no multipole interactions.
   */
  public static Force constructForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ParticleMeshEwald pme = forceFieldEnergy.getPmeNode();
    if (pme == null) {
      return null;
    }
    return new AmoebaMultipoleForce(topology, openMMDualTopologyEnergy);
  }

  /**
   * Update the force parameters for the AMOEBA Multipole Force.
   *
   * @param atoms        The array of Atoms for which the force parameters are to be updated.
   * @param openMMEnergy The OpenMMEnergy instance that contains the multipole information.
   */
  public void updateForce(Atom[] atoms, OpenMMEnergy openMMEnergy) {
    ParticleMeshEwald pme = openMMEnergy.getPmeNode();
    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    double doPolarization = 1.0;
    if (pme.getPolarizationType() == Polarization.NONE) {
      doPolarization = 0.0;
    }

    DoubleArray dipoles = new DoubleArray(3);
    DoubleArray quadrupoles = new DoubleArray(9);

    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    boolean lambdaTerm = pme.getLambdaTerm();
    if (lambdaTerm) {
      AlchemicalParameters.AlchemicalMode alchemicalMode = alchemicalParameters.mode;
      // Only scale mode is supported for OpenMM.
      if (alchemicalMode != AlchemicalParameters.AlchemicalMode.SCALE) {
        logger.severe(format(" Alchemical mode %s not supported for OpenMM.", alchemicalMode));
      }
      // Permanent multipole softcore is not supported for OpenMM.
      if (alchemicalParameters.permLambdaAlpha != 0.0) {
        logger.severe(" Permanent multipole softcore not supported for OpenMM.");
      }
      // Isolated ligand electrostatics are not supported for OpenMM.
      if (alchemicalParameters.doLigandGKElec || alchemicalParameters.doLigandVaporElec) {
        logger.severe(" Isolated ligand electrostatics are not supported for OpenMM.");
      }
      // Condensed SCF without a ligand is not supported for OpenMM.
      if (alchemicalParameters.doNoLigandCondensedSCF) {
        logger.severe(" Condensed SCF without a ligand is not supported for OpenMM.");
      }
    }

    double permLambda = alchemicalParameters.permLambda;
    double polarLambda = alchemicalParameters.polLambda;

    for (Atom atom : atoms) {
      int index = atom.getXyzIndex() - 1;
      MultipoleType multipoleType = pme.getMultipoleType(index);
      PolarizeType polarizeType = pme.getPolarizeType(index);
      int[] axisAtoms = atom.getAxisAtomIndices();

      double permScale = 1.0;
      double polarScale = doPolarization;

      if (!atom.getUse() || !atom.getElectrostatics()) {
        permScale = 0.0;
        polarScale = 0.0;
      }

      if (atom.applyLambda()) {
        permScale *= permLambda;
        polarScale *= polarLambda;
      }

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
        dipoles.set(j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * permScale);
      }
      int l = 0;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0 * permScale);
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
      setMultipoleParameters(index, multipoleType.charge * permScale,
          dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
          polarizeType.thole, polarizeType.pdamp * dampingFactorConversion,
          polarizeType.polarizability * polarityConversion * polarScale);
    }

    dipoles.destroy();
    quadrupoles.destroy();
    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update the force parameters for the AMOEBA Multipole Force in a dual topology system.
   *
   * @param atoms                    The array of Atoms for which the force parameters are to be updated.
   * @param topology                 The topology index for the dual topology system.
   * @param openMMDualTopologyEnergy The OpenMM Dual-Topology Energy instance that contains the multipole parameters.
   */
  public void updateForce(Atom[] atoms, int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ParticleMeshEwald pme = forceFieldEnergy.getPmeNode();
    double quadrupoleConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double polarityConversion = OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom * OpenMM_NmPerAngstrom;
    double dampingFactorConversion = sqrt(OpenMM_NmPerAngstrom);

    double doPolarization = 1.0;
    if (pme.getPolarizationType() == Polarization.NONE) {
      doPolarization = 0.0;
    }

    // Dual-topology scale factor.
    double scaleDT = Math.sqrt(openMMDualTopologyEnergy.getTopologyScale(topology));

    DoubleArray dipoles = new DoubleArray(3);
    DoubleArray quadrupoles = new DoubleArray(9);

    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    boolean lambdaTerm = pme.getLambdaTerm();
    if (lambdaTerm) {
      AlchemicalParameters.AlchemicalMode alchemicalMode = alchemicalParameters.mode;
      // Only scale mode is supported for OpenMM.
      if (alchemicalMode != AlchemicalParameters.AlchemicalMode.SCALE) {
        logger.severe(format(" Alchemical mode %s not supported for OpenMM.", alchemicalMode));
      }
      // Permanent multipole softcore is not supported for OpenMM.
      if (alchemicalParameters.permLambdaAlpha != 0.0) {
        logger.severe(" Permanent multipole softcore not supported for OpenMM.");
      }
      // Isolated ligand electrostatics are not supported for OpenMM.
      if (alchemicalParameters.doLigandGKElec || alchemicalParameters.doLigandVaporElec) {
        logger.severe(" Isolated ligand electrostatics are not supported for OpenMM.");
      }
      // Condensed SCF without a ligand is not supported for OpenMM.
      if (alchemicalParameters.doNoLigandCondensedSCF) {
        logger.severe(" Condensed SCF without a ligand is not supported for OpenMM.");
      }
    }

    double permLambda = alchemicalParameters.permLambda;
    double polarLambda = alchemicalParameters.polLambda;

    for (Atom atom : atoms) {
      if (atom.getTopologyIndex() != topology) {
        // Skip atoms that are not in this topology.
        continue;
      }

      int index = atom.getArrayIndex();
      MultipoleType multipoleType = pme.getMultipoleType(index);
      PolarizeType polarizeType = pme.getPolarizeType(index);
      int[] axisAtoms = atom.getAxisAtomIndices();

      double permScale = scaleDT;
      double polarScale = doPolarization;
      if (!atom.getUse() || !atom.getElectrostatics()) {
        permScale = 0.0;
        polarScale = 0.0;
      }

      if (atom.applyLambda()) {
        permScale *= permLambda;
        polarScale *= polarLambda;
      }

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
        dipoles.set(j, multipoleType.dipole[j] * OpenMM_NmPerAngstrom * permScale);
      }
      int l = 0;
      for (int j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
          quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0 * permScale);
        }
      }

      int zaxis = 1;
      int xaxis = 1;
      int yaxis = 1;

      if (axisAtoms != null) {
        zaxis = axisAtoms[0];
        zaxis = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, zaxis);
        if (axisAtoms.length > 1) {
          xaxis = axisAtoms[1];
          xaxis = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, xaxis);
          if (axisAtoms.length > 2) {
            yaxis = axisAtoms[2];
            yaxis = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, yaxis);
          }
        }
      } else {
        axisType = OpenMM_AmoebaMultipoleForce_NoAxisType;
      }

      // Set the multipole parameters.
      int dualTopologyIndex = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, index);
      setMultipoleParameters(dualTopologyIndex, multipoleType.charge * permScale,
          dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
          polarizeType.thole, polarizeType.pdamp * dampingFactorConversion,
          polarizeType.polarizability * polarityConversion * polarScale);
    }

    dipoles.destroy();
    quadrupoles.destroy();
    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
