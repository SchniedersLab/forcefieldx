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
package ffx.potential.ommffm;

import ffx.crystal.Crystal;
import ffx.openmm.ffm.DoubleArray;
import ffx.openmm.ffm.Force;
import ffx.openmm.ffm.IntArray;
import ffx.openmm.ffm.amoeba.MultipoleForce;
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

import static ffx.openmm.ffm.OpenMMUnits.NM_PER_ANGSTROM;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

/**
 * AMOEBA Multipole Force backed by FFM {@link MultipoleForce}.
 */
public class AmoebaMultipoleForce extends MultipoleForce {

  private static final Logger logger = Logger.getLogger(AmoebaMultipoleForce.class.getName());

  public static final int AXIS_TYPE_Z_THEN_X = 0;
  public static final int AXIS_TYPE_BISECTOR = 1;
  public static final int AXIS_TYPE_Z_BISECT = 2;
  public static final int AXIS_TYPE_THREE_FOLD = 3;
  public static final int AXIS_TYPE_Z_ONLY = 4;
  public static final int AXIS_TYPE_NO_AXIS_TYPE = 5;

  public static final int COVALENT_12 = 0;
  public static final int COVALENT_13 = 1;
  public static final int COVALENT_14 = 2;
  public static final int COVALENT_15 = 3;
  public static final int POLARIZATION_COVALENT_11 = 4;
  public static final int POLARIZATION_COVALENT_12 = 5;
  public static final int POLARIZATION_COVALENT_13 = 6;
  public static final int POLARIZATION_COVALENT_14 = 7;

  public static final int NONBONDED_METHOD_NO_CUTOFF = 0;
  public static final int NONBONDED_METHOD_PME = 1;

  public static final int POLARIZATION_TYPE_MUTUAL = 0;
  public static final int POLARIZATION_TYPE_DIRECT = 1;
  public static final int POLARIZATION_TYPE_EXTRAPOLATED = 2;

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
    double quadrupoleConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double polarityConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double dampingFactorConversion = sqrt(NM_PER_ANGSTROM);

    boolean lambdaTerm = pme.getLambdaTerm();
    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    double permLambda = alchemicalParameters.permLambda;
    double polarLambda = alchemicalParameters.polLambda;

    Atom[] atoms = openMMEnergy.getMolecularAssembly().getAtomArray();
    int nAtoms = atoms.length;

    try (DoubleArray dipoles = new DoubleArray(3);
         DoubleArray quadrupoles = new DoubleArray(9)) {
      for (int i = 0; i < nAtoms; i++) {
        Atom atom = atoms[i];
        MultipoleType multipoleType = pme.getMultipoleType(i);
        PolarizeType polarType = pme.getPolarizeType(i);

        // Define the frame definition.
        int axisType = switch (multipoleType.frameDefinition) {
          case NONE -> AXIS_TYPE_NO_AXIS_TYPE;
          case ZONLY -> AXIS_TYPE_Z_ONLY;
          case ZTHENX -> AXIS_TYPE_Z_THEN_X;
          case BISECTOR -> AXIS_TYPE_BISECTOR;
          case ZTHENBISECTOR -> AXIS_TYPE_Z_BISECT;
          case THREEFOLD -> AXIS_TYPE_THREE_FOLD;
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
          dipoles.set(j, multipoleType.dipole[j] * NM_PER_ANGSTROM * permScale);
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
          axisType = AXIS_TYPE_NO_AXIS_TYPE;
        }

        double charge = multipoleType.charge * permScale;

        // Add the multipole.
        addMultipole(charge, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis, polarType.thole,
            polarType.pdamp * dampingFactorConversion, polarType.polarizability * polarityConversion * polarScale);
      }
    }

    int[][] ip11 = pme.getPolarization11();
    try (IntArray covalentMap = new IntArray(0)) {
      for (int i = 0; i < nAtoms; i++) {
        Atom ai = atoms[i];

        // 1-2 Mask
        covalentMap.resize(0);
        for (Atom ak : ai.get12List()) {
          covalentMap.append(ak.getIndex() - 1);
        }
        setCovalentMap(i, COVALENT_12, covalentMap);

        // 1-3 Mask
        covalentMap.resize(0);
        for (Atom ak : ai.get13List()) {
          covalentMap.append(ak.getIndex() - 1);
        }
        setCovalentMap(i, COVALENT_13, covalentMap);

        // 1-4 Mask
        covalentMap.resize(0);
        for (Atom ak : ai.get14List()) {
          covalentMap.append(ak.getIndex() - 1);
        }
        setCovalentMap(i, COVALENT_14, covalentMap);

        // 1-5 Mask
        covalentMap.resize(0);
        for (Atom ak : ai.get15List()) {
          covalentMap.append(ak.getIndex() - 1);
        }
        setCovalentMap(i, COVALENT_15, covalentMap);

        // 1-1 Polarization Groups.
        covalentMap.resize(0);
        for (int j = 0; j < ip11[i].length; j++) {
          covalentMap.append(ip11[i][j]);
        }
        setCovalentMap(i, POLARIZATION_COVALENT_11, covalentMap);
      }
    }
  }

  /**
   * Construct a Dual Topology AMOEBA Multipole Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public AmoebaMultipoleForce(int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ParticleMeshEwald pme = forceFieldEnergy.getPmeNode();
    if (pme == null) {
      destroy();
      return;
    }

    double doPolarization = configureForce(forceFieldEnergy);
    boolean lambdaTerm = pme.getLambdaTerm();
    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    double permLambda = alchemicalParameters.permLambda;
    double polarLambda = alchemicalParameters.polLambda;

    int otherTopology = (topology == 0) ? 1 : 0;
    ForceFieldEnergy otherForceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(otherTopology);
    ParticleMeshEwald otherPME = otherForceFieldEnergy.getPmeNode();

    double quadrupoleConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double polarityConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double dampingFactorConversion = sqrt(NM_PER_ANGSTROM);

    try (DoubleArray dipoles = new DoubleArray(3);
         DoubleArray quadrupoles = new DoubleArray(9)) {
      int nAtoms = openMMDualTopologyEnergy.getNumberOfAtoms();
      for (int i = 0; i < nAtoms; i++) {
        Atom atom = openMMDualTopologyEnergy.getDualTopologyAtom(topology, i);
        int top = atom.getTopologyIndex();
        if (top == topology) {
          int index = atom.getArrayIndex();
          MultipoleType multipoleType = pme.getMultipoleType(index);
          PolarizeType polarType = pme.getPolarizeType(index);

          // Define the frame definition.
          int axisType = switch (multipoleType.frameDefinition) {
            case NONE -> AXIS_TYPE_NO_AXIS_TYPE;
            case ZONLY -> AXIS_TYPE_Z_ONLY;
            case ZTHENX -> AXIS_TYPE_Z_THEN_X;
            case BISECTOR -> AXIS_TYPE_BISECTOR;
            case ZTHENBISECTOR -> AXIS_TYPE_Z_BISECT;
            case THREEFOLD -> AXIS_TYPE_THREE_FOLD;
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
            dipoles.set(j, multipoleType.dipole[j] * NM_PER_ANGSTROM * permScale);
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
            axisType = AXIS_TYPE_NO_AXIS_TYPE;
          }

          double charge = multipoleType.charge * permScale;

          // Add the multipole.
          addMultipole(charge, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis, polarType.thole,
              polarType.pdamp * dampingFactorConversion, polarType.polarizability * polarityConversion * polarScale);
        } else {
          // Add a fake multipole for an atom not in this topology.
          int axisType = AXIS_TYPE_NO_AXIS_TYPE;
          double charge = 0.0;
          for (int j = 0; j < 3; j++) {
            dipoles.set(j, 0.0);
          }
          for (int j = 0; j < 9; j++) {
            quadrupoles.set(j, 0.0);
          }
          int zaxis = -1;
          int xaxis = -1;
          int yaxis = -1;
          double thole = 0.39;
          double pdamp = 1.0;
          double polarizability = 0.0;
          addMultipole(charge, dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis, thole,
              pdamp, polarizability);
        }
      }

      int[][] ip11 = pme.getPolarization11();
      int[][] ip11Other = (otherPME != null) ? otherPME.getPolarization11() : null;

      try (IntArray covalentMap = new IntArray(0)) {
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
          setCovalentMap(i, COVALENT_12, covalentMap);

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
          setCovalentMap(i, COVALENT_13, covalentMap);

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
          setCovalentMap(i, COVALENT_14, covalentMap);

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
          setCovalentMap(i, COVALENT_15, covalentMap);

          // 1-1 Polarization Groups.
          covalentSet.clear();
          if (ip11 != null && index < ip11.length && ip11[index] != null) {
            for (int k : ip11[index]) {
              int value = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, k);
              covalentSet.add(value);
            }
          }
          if (ip11Other != null && otherIndex < ip11Other.length && ip11Other[otherIndex] != null) {
            for (int k : ip11Other[otherIndex]) {
              int value = openMMDualTopologyEnergy.mapToDualTopologyIndex(otherTopology, k);
              covalentSet.add(value);
            }
          }
          covalentMap.resize(0);
          for (int k : covalentSet) {
            covalentMap.append(k);
          }
          setCovalentMap(i, POLARIZATION_COVALENT_11, covalentMap);
        }
      }
    }
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
    SCFAlgorithm scfAlgorithm;
    if (polarization != Polarization.MUTUAL) {
      setPolarizationType(POLARIZATION_TYPE_DIRECT);
      if (pme.getPolarizationType() == Polarization.NONE) {
        doPolarization = 0.0;
      }
    } else {
      String algorithm = forceField.getString("SCF_ALGORITHM", "CG");
      try {
        algorithm = algorithm.replace("-", "_").toUpperCase();
        scfAlgorithm = SCFAlgorithm.valueOf(algorithm);
      } catch (Exception e) {
        scfAlgorithm = SCFAlgorithm.CG;
      }
      if (scfAlgorithm == SCFAlgorithm.EPT) {
        setPolarizationType(POLARIZATION_TYPE_EXTRAPOLATED);
        try (DoubleArray exptCoefficients = new DoubleArray(4)) {
          exptCoefficients.set(0, -0.154);
          exptCoefficients.set(1, 0.017);
          exptCoefficients.set(2, 0.657);
          exptCoefficients.set(3, 0.475);
          setExtrapolationCoefficients(exptCoefficients);
        }
      } else {
        setPolarizationType(POLARIZATION_TYPE_MUTUAL);
      }
    }

    Crystal crystal = forceFieldEnergy.getCrystal();
    double cutoff = pme.getEwaldCutoff();
    double aewald = pme.getEwaldCoefficient();
    if (!crystal.aperiodic()) {
      setNonbondedMethod(NONBONDED_METHOD_PME);
      setCutoffDistance(cutoff * NM_PER_ANGSTROM);
      double ewaldTolerance = 1.0e-04;
      setEwaldErrorTolerance(ewaldTolerance);
      ReciprocalSpace recip = pme.getReciprocalSpace();
      int nx = recip.getXDim();
      int ny = recip.getYDim();
      int nz = recip.getZDim();
      setPMEParameters(aewald / NM_PER_ANGSTROM, nx, ny, nz);
    } else {
      setNonbondedMethod(NONBONDED_METHOD_NO_CUTOFF);
    }

    setMutualInducedMaxIterations(500);
    double poleps = pme.getPolarEps();
    setMutualInducedTargetEpsilon(poleps);

    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    boolean lambdaTerm = pme.getLambdaTerm();
    if (lambdaTerm) {
      AlchemicalParameters.AlchemicalMode alchemicalMode = alchemicalParameters.mode;
      if (alchemicalMode != AlchemicalParameters.AlchemicalMode.SCALE) {
        logger.severe(format(" Alchemical mode %s not supported for OpenMM.", alchemicalMode));
      }
      if (alchemicalParameters.permLambdaAlpha != 0.0) {
        logger.severe(" Permanent multipole softcore not supported for OpenMM.");
      }
      if (alchemicalParameters.doLigandGKElec || alchemicalParameters.doLigandVaporElec) {
        logger.severe(" Isolated ligand electrostatics are not supported for OpenMM.");
      }
      if (alchemicalParameters.doNoLigandCondensedSCF) {
        logger.severe(" Condensed SCF without a ligand is not supported for OpenMM.");
      }
    }

    int forceGroup = forceField.getInteger("PME_FORCE_GROUP", 1);
    setForceGroup(forceGroup);
    if (logger.isLoggable(Level.INFO)) {
      StringBuilder sb = new StringBuilder();
      sb.append(format("  Multipole Force \t\t\t%d\n", forceGroup));
      sb.append(format("   Polarization:                %10s\n", polarization));
      if (polarization == Polarization.MUTUAL) {
        sb.append(format("   Mutual Target Epsilon:       %10.2e", poleps));
      }
      logger.info(sb.toString());
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
   * Convenience method to construct a Dual Topology AMOEBA Multipole Force.
   *
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
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
    double quadrupoleConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double polarityConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double dampingFactorConversion = sqrt(NM_PER_ANGSTROM);

    double doPolarization = 1.0;
    if (pme.getPolarizationType() == Polarization.NONE) {
      doPolarization = 0.0;
    }

    AlchemicalParameters alchemicalParameters = pme.getAlchemicalParameters();
    boolean lambdaTerm = pme.getLambdaTerm();
    if (lambdaTerm) {
      AlchemicalParameters.AlchemicalMode alchemicalMode = alchemicalParameters.mode;
      if (alchemicalMode != AlchemicalParameters.AlchemicalMode.SCALE) {
        logger.severe(format(" Alchemical mode %s not supported for OpenMM.", alchemicalMode));
      }
      if (alchemicalParameters.permLambdaAlpha != 0.0) {
        logger.severe(" Permanent multipole softcore not supported for OpenMM.");
      }
      if (alchemicalParameters.doLigandGKElec || alchemicalParameters.doLigandVaporElec) {
        logger.severe(" Isolated ligand electrostatics are not supported for OpenMM.");
      }
      if (alchemicalParameters.doNoLigandCondensedSCF) {
        logger.severe(" Condensed SCF without a ligand is not supported for OpenMM.");
      }
    }

    double permLambda = alchemicalParameters.permLambda;
    double polarLambda = alchemicalParameters.polLambda;

    try (DoubleArray dipoles = new DoubleArray(3);
         DoubleArray quadrupoles = new DoubleArray(9)) {
      for (Atom atom : atoms) {
        int index = atom.getArrayIndex();
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
          case NONE -> AXIS_TYPE_NO_AXIS_TYPE;
          case ZONLY -> AXIS_TYPE_Z_ONLY;
          case ZTHENX -> AXIS_TYPE_Z_THEN_X;
          case BISECTOR -> AXIS_TYPE_BISECTOR;
          case ZTHENBISECTOR -> AXIS_TYPE_Z_BISECT;
          case THREEFOLD -> AXIS_TYPE_THREE_FOLD;
        };

        // Load local multipole coefficients.
        for (int j = 0; j < 3; j++) {
          dipoles.set(j, multipoleType.dipole[j] * NM_PER_ANGSTROM * permScale);
        }
        int l = 0;
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0 * permScale);
          }
        }

        int zaxis = -1;
        int xaxis = -1;
        int yaxis = -1;

        if (axisAtoms != null) {
          zaxis = axisAtoms[0];
          if (axisAtoms.length > 1) {
            xaxis = axisAtoms[1];
            if (axisAtoms.length > 2) {
              yaxis = axisAtoms[2];
            }
          }
        } else {
          axisType = AXIS_TYPE_NO_AXIS_TYPE;
        }

        // Set the multipole parameters.
        setMultipoleParameters(index, multipoleType.charge * permScale,
            dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
            polarizeType.thole, polarizeType.pdamp * dampingFactorConversion,
            polarizeType.polarizability * polarityConversion * polarScale);
      }
    }

    updateParametersInContext(openMMEnergy.getContext());
  }

  /**
   * Update existing AMOEBA multipole force for the Dual-Topology OpenMM System.
   *
   * @param atoms                    The array of Atoms for which the force parameters are to be updated.
   * @param topology                 The topology index for the OpenMM System.
   * @param openMMDualTopologyEnergy The OpenMMDualTopologyEnergy instance.
   */
  public void updateForce(Atom[] atoms, int topology, OpenMMDualTopologyEnergy openMMDualTopologyEnergy) {
    ForceFieldEnergy forceFieldEnergy = openMMDualTopologyEnergy.getForceFieldEnergy(topology);
    ParticleMeshEwald pme = forceFieldEnergy.getPmeNode();
    double quadrupoleConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double polarityConversion = NM_PER_ANGSTROM * NM_PER_ANGSTROM * NM_PER_ANGSTROM;
    double dampingFactorConversion = sqrt(NM_PER_ANGSTROM);

    double doPolarization = 1.0;
    if (pme.getPolarizationType() == Polarization.NONE) {
      doPolarization = 0.0;
    }

    // Dual-topology scale factor.
    double scaleDT = Math.sqrt(openMMDualTopologyEnergy.getTopologyScale(topology));

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

    try (DoubleArray dipoles = new DoubleArray(3);
         DoubleArray quadrupoles = new DoubleArray(9)) {
      for (Atom atom : atoms) {
        if (atom.getTopologyIndex() != topology) {
          // Skip atoms not in this topology.
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
          case NONE -> AXIS_TYPE_NO_AXIS_TYPE;
          case ZONLY -> AXIS_TYPE_Z_ONLY;
          case ZTHENX -> AXIS_TYPE_Z_THEN_X;
          case BISECTOR -> AXIS_TYPE_BISECTOR;
          case ZTHENBISECTOR -> AXIS_TYPE_Z_BISECT;
          case THREEFOLD -> AXIS_TYPE_THREE_FOLD;
        };

        // Load local multipole coefficients.
        for (int j = 0; j < 3; j++) {
          dipoles.set(j, multipoleType.dipole[j] * NM_PER_ANGSTROM * permScale);
        }
        int l = 0;
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 3; k++) {
            quadrupoles.set(l++, multipoleType.quadrupole[j][k] * quadrupoleConversion / 3.0 * permScale);
          }
        }

        int zaxis = -1;
        int xaxis = -1;
        int yaxis = -1;

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
          axisType = AXIS_TYPE_NO_AXIS_TYPE;
        }

        // Set the multipole parameters.
        int dualTopologyIndex = openMMDualTopologyEnergy.mapToDualTopologyIndex(topology, index);
        setMultipoleParameters(dualTopologyIndex, multipoleType.charge * permScale,
            dipoles, quadrupoles, axisType, zaxis, xaxis, yaxis,
            polarizeType.thole, polarizeType.pdamp * dampingFactorConversion,
            polarizeType.polarizability * polarityConversion * polarScale);
      }
    }

    updateParametersInContext(openMMDualTopologyEnergy.getContext());
  }

}
