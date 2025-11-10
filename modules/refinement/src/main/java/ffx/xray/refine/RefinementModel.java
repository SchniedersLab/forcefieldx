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
package ffx.xray.refine;

import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.MSNode;
import ffx.potential.bonded.Molecule;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import org.apache.commons.configuration2.CompositeConfiguration;

import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.numerics.math.ScalarMath.b2u;
import static java.lang.String.format;

/**
 * RefinementModel class.
 *
 * @author Timothy D. Fenn
 * @since 1.0
 */
public class RefinementModel {

  private static final Logger logger = Logger.getLogger(RefinementModel.class.getName());

  /**
   * An array of MolecularAssembly objects representing the different molecular assemblies
   * being processed.
   */
  private final MolecularAssembly[] molecularAssemblies;

  /**
   * The {@code RefinementMode} defines the approach used for refining model parameters.
   */
  private RefinementMode refinementMode;

  /**
   * Add a 6-parameter anisotropic b-factor to each heavy atom.
   */
  private final boolean addAnisou;

  /**
   * Refine b-factors by residue, where every atom in a residue has
   * the same b-factor.
   */
  private final boolean byResidue;

  /**
   * If byResidue is true, then nResiduePerBFactor groups bonded residues
   * together into residue groups.
   */
  private final int nResiduePerBFactor;

  /**
   * If true, hydrogen atom b-factors are taken from their heavy atom.
   * This also enforced using byResidue b-factor refinement (i.e., this flag is redundant).
   */
  private final boolean ridingHydrogen;

  /**
   * B-factor mass for extended Lagrangian.
   */
  private final double bMass;

  /**
   * A boolean flag that determines whether the molecular occupancies are subject
   * to refinement.
   */
  private final boolean refineMolOcc;

  /**
   * If true, H/D occupancy will each be set to 0.5.
   */
  private final boolean resetHDOccupancy;

  /**
   * Occupancy mass for extended Lagrangian.
   */
  private final double occMass;

  /**
   * All atoms that scatter are included in this array, including inactive atoms whose
   * atomic coordinates are frozen.
   */
  private final Atom[] scatteringAtoms;

  /**
   * This list has the same Atom instances as the refinedCoordinates in a defined order.
   */
  private final List<Atom> coordinateAtomList = new ArrayList<>();

  /**
   * All atomic coordinates that are being refined (inactive atoms are excluded).
   */
  private final Map<Atom, RefinedCoordinates> refinedCoordinates;

  /**
   * This list provides a defined order for the Atom instances in the refinedBFactors map.
   */
  private final List<Atom> bFactorAtomList = new ArrayList<>();

  /**
   * All active b-factors that are being refined.
   */
  private final Map<Atom, RefinedBFactor> refinedBFactors;

  /**
   * B-factors restraint list.
   */
  private final List<Atom[]> bFactorRestraints = new ArrayList<>();

  /**
   * This list has the same Atom instances as the refinedOccupancies in a defined order.
   */
  private final List<Atom> occupancyAtomList = new ArrayList<>();

  /**
   * All occupancies that are being refined.
   */
  private final Map<Atom, RefinedOccupancy> refinedOccupancies;

  /**
   * A list that holds all the refined parameters that are instances of {@code RefinedParameter},
   */
  private final List<RefinedParameter> allParametersList;

  /**
   * Each List should contain a Residue instance from each MolecularAssemblies (e.g., A, B and C).
   */
  private final List<List<Residue>> altResidues;

  /**
   * Each List should contain a Molecule instance from each MolecularAssemblies (e.g., A, B and C).
   */
  private final List<List<Molecule>> altMolecules;

  /**
   * Constructor for RefinementModel.
   *
   * @param molecularAssemblies an array of {@link ffx.potential.MolecularAssembly} objects.
   */
  public RefinementModel(MolecularAssembly[] molecularAssemblies) {
    this(RefinementMode.COORDINATES_AND_BFACTORS_AND_OCCUPANCIES, molecularAssemblies);
  }

  /**
   * Constructor for RefinementModel.
   *
   * @param refinementMode      The refinement mode in use.
   * @param molecularAssemblies an array of {@link ffx.potential.MolecularAssembly} objects.
   */
  public RefinementModel(RefinementMode refinementMode, MolecularAssembly[] molecularAssemblies) {

    this.refinementMode = refinementMode;
    this.molecularAssemblies = molecularAssemblies;
    MolecularAssembly rootAssembly = molecularAssemblies[0];

    // Load b-factor refinement properties.
    CompositeConfiguration properties = rootAssembly.getProperties();
    addAnisou = properties.getBoolean("add-anisou", false);
    byResidue = properties.getBoolean("residue-bfactor", false);
    nResiduePerBFactor = properties.getInt("n-residue-bfactor", 1);
    ridingHydrogen = properties.getBoolean("riding-hydrogen-bfactor", true);
    bMass = properties.getDouble("bfactor-mass", 5.0);
    occMass = properties.getDouble("occupancy-mass", 10.0);
    resetHDOccupancy = properties.getBoolean("reset-hd-occupancy", false);

    // Load occupancy refinement properties
    refineMolOcc = properties.getBoolean("refine-mol-occ", false);

    // Add anisotropic temperature factors
    if (addAnisou) {
      addAnisotropicBFactors();
    }

    // Regularize active atoms to be consistent with the b-factor refinement strategy.
    if (refinementMode.includesBFactors()) {
      regularizeActiveAtoms();
    }

    // Collect the atoms that scatter and atoms whose coordinates will be refined.
    List<Atom> scatteringList = new ArrayList<>();
    refinedCoordinates = createCoordinateModel(scatteringList);
    scatteringAtoms = scatteringList.toArray(new Atom[0]);

    // Collect atoms whose b-factors will be refined.
    refinedBFactors = createBFactorModel();

    // Collect atoms whose occupancy will be refined.
    altResidues = new ArrayList<>();
    altMolecules = new ArrayList<>();
    refinedOccupancies = createOccupancyModel();

    // Collect all refined parameters given the current RefinementMode
    allParametersList = new ArrayList<>();
    setRefinementMode(refinementMode);
  }

  /**
   * Sets the refinement mode and adjusts the refined parameters based on the provided mode.
   * Clears the current list of refined parameters and populates it with the relevant parameters
   * (coordinates, B-factors, or occupancies) if they are included in the specified refinement mode.
   * Finally, updates the parameter indices.
   *
   * @param mode the refinement mode to set, which determines the types of parameters
   *             to include (e.g., coordinates, B-factors, occupancies)
   */
  public void setRefinementMode(RefinementMode mode) {
    this.refinementMode = mode;
    allParametersList.clear();
    /*
     * The parameter maps do not maintain a defined order of the key-value pairs.
     * For this reason, the parameter lists are iterated over to create the overall list.
     */
    if (refinementMode.includesCoordinates()) {
      for (Atom atom : coordinateAtomList) {
        allParametersList.add(refinedCoordinates.get(atom));
      }
    }
    if (refinementMode.includesBFactors()) {
      for (Atom atom : bFactorAtomList) {
        allParametersList.add(refinedBFactors.get(atom));
      }
      collectBFactorRestraints();
    } else {
      bFactorRestraints.clear();
    }
    if (refinementMode.includesOccupancies()) {
      for (Atom atom : occupancyAtomList) {
        allParametersList.add(refinedOccupancies.get(atom));
      }
    }
    setParameterIndices();
  }

  /**
   * Getter for the field <code>totalAtomArray</code>.
   *
   * @return the totalAtomArray
   */
  public Atom[] getScatteringAtoms() {
    return scatteringAtoms;
  }

  /**
   * Getter for the field <code>activeAtomArray</code>.
   *
   * @return the activeAtomArray
   */
  public Atom[] getActiveAtoms() {
    return coordinateAtomList.toArray(new Atom[0]);
  }

  /**
   * Getter for the field <code>altMolecules</code>.
   *
   * @return the altMolecules
   */
  public List<List<Molecule>> getAltMolecules() {
    return altMolecules;
  }

  /**
   * Getter for the field <code>altResidues</code>.
   *
   * @return the altResidues
   */
  public List<List<Residue>> getAltResidues() {
    return altResidues;
  }

  /**
   * Retrieves the array of MolecularAssembly objects.
   *
   * @return An array of MolecularAssembly objects.
   */
  public MolecularAssembly[] getMolecularAssemblies() {
    return molecularAssemblies;
  }

  /**
   * List of Atom pairs that define B-factor restraints. There is a covalent bond between each pair.
   *
   * @return The bond list.
   */
  public List<Atom[]> getBFactorRestraints() {
    return bFactorRestraints;
  }

  /**
   * Adds the coordinate gradient from an alternate conformer to the overall coordinate gradient. Each
   * active atom from the alternate conformer stores its XYZ gradient and its index into the overall
   * refinement gradient array.
   *
   * @param assembly The index of the molecular assembly to retrieve active atoms from.
   * @param gradient The overall gradient array where the computed gradient data is aggregated.
   */
  public void addAssemblyGradient(int assembly, double[] gradient) {
    MolecularAssembly molecularAssembly = molecularAssemblies[assembly];
    Atom[] activeAtoms = molecularAssembly.getActiveAtomArray();
    double[] xyz = new double[3];
    for (Atom a : activeAtoms) {
      int index = a.getXrayCoordIndex() * 3;
      a.getXYZGradient(xyz);
      gradient[index] += xyz[0];
      gradient[index + 1] += xyz[1];
      gradient[index + 2] += xyz[2];
    }
  }

  /**
   * Provides a string representation of the refinement model, including details
   * about the number of atoms, active atoms, atoms in use, and the refinement variables.
   *
   * @return A string describing the refinement model, including the total
   * number of atoms, the number of atoms currently being used,
   * the number of active atoms, and the number of refinement variables
   * categorized by XYZ coordinates, B-Factors, and occupancies.
   */
  public String toString() {
    int nAtoms = scatteringAtoms.length;
    int nActive = coordinateAtomList.size();
    // Count the number of scatteringAtoms in use.
    int nUse = 0;
    for (Atom a : scatteringAtoms) {
      if (a.getUse()) {
        nUse++;
      }
    }
    int nXYZ = getNumCoordParameters();
    int nBFactors = getNumBFactorParameters();
    int nOccupancies = getNumOccupancyParameters();
    int n = nXYZ + nBFactors + nOccupancies;

    StringBuilder sb = new StringBuilder("\n Refinement Model\n");
    sb.append(format("  Number of atoms:        %d\n", nAtoms));
    sb.append(format("  Atoms being used:       %d\n", nUse));
    sb.append(format("  Number of active atoms: %d\n", nActive));
    sb.append(format("  Number of variables:    %d (nXYZ %d, nB %d, nOcc %d)\n",
        n, nXYZ, nBFactors, nOccupancies));

    return sb.toString();
  }

  /**
   * Retrieves the current refinement mode set in the object.
   *
   * @return the refinement mode of the object as a RefinementMode enum.
   */
  public RefinementMode getRefinementMode() {
    return refinementMode;
  }

  /**
   * Retrieves a list of refined parameters.
   *
   * @return a list containing the refined parameters.
   */
  public List<RefinedParameter> getRefinedParameters() {
    return allParametersList;
  }

  /**
   * Calculates the total number of parameters to be refined based on the specified refinement mode.
   * The parameters may include atomic coordinates, B-factors, and occupancies depending on the mode.
   *
   * @return The total number of parameters that need to be refined based on the specified mode.
   */
  public int getNumParameters() {
    return getNumCoordParameters()
        + getNumBFactorParameters()
        + getNumOccupancyParameters();
  }

  /**
   * Calculates the number of coordinate parameters to be refined based on the specified refinement mode.
   * If the mode includes coordinates, the number of parameters is determined by the size of the coordinate atom list.
   *
   * @return The number of coordinate parameters that need to be refined. Returns 0 if the mode does not include coordinates.
   */
  public int getNumCoordParameters() {
    // Coordinates
    if (refinementMode.includesCoordinates()) {
      return coordinateAtomList.size() * 3;
    }
    return 0;
  }

  /**
   * Calculates the total number of B-factor parameters based on the specified refinement mode.
   *
   * @return the total number of B-factor parameters, considering whether they are anisotropic or isotropic.
   */
  public int getNumBFactorParameters() {
    int num = 0;
    if (refinementMode.includesBFactors()) {
      for (RefinedBFactor bFactor : refinedBFactors.values()) {
        num += bFactor.getNumberOfParameters();
      }
    }
    return num;
  }

  /**
   * Calculates the number of occupancy parameters based on the provided refinement mode.
   *
   * @return the number of occupancy parameters if occupancies are included; otherwise, returns 0
   */
  public int getNumOccupancyParameters() {
    if (refinementMode.includesOccupancies()) {
      return occupancyAtomList.size();
    }
    return 0;
  }

  /**
   * Counts and returns the total number of ANISOU records present.
   *
   * @return the number of ANISOU records as an integer
   */
  public int getNumANISOU() {
    int numANISOU = 0;
    for (RefinedBFactor bFactor : refinedBFactors.values()) {
      if (bFactor.isAnisou()) {
        numANISOU++;
      }
    }
    return numANISOU;
  }

  /**
   * Get parameter values and store them into the provided array based on the refinement mode.
   * This method extracts and sets coordinates, B-factors, and occupancies for the atoms
   * depending on the specified refinement mode.
   *
   * @param x The array into which the parameter values are loaded. The array should be large
   *          enough to accommodate all required parameters based on the refinement mode.
   */
  public void getParameters(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.getParameters(x);
    }
  }

  /**
   * Set parameter values into the refinement model based on the provided refinement mode.
   * The method updates coordinates, B-factors, and occupancies for atoms, depending on what is
   * specified by the refinement mode.
   *
   * @param x An array of doubles containing the new parameter values. The values should be ordered as
   *          required by the refinement mode: first coordinates, then B-factors, and finally occupancies,
   *          if applicable.
   */
  public void setParameters(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.setParameters(x);
    }
  }

  /**
   * Get parameter velocities and store them into the provided array based on the refinement mode.
   * This method extracts and sets coordinates, B-factors, and occupancies for the atoms
   * depending on the specified refinement mode.
   *
   * @param x The array into which the parameter velocities are loaded. The array should be large
   *          enough to accommodate all required parameter velocities based on the refinement mode.
   */
  public void getVelocity(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.getVelocity(x);
    }
  }

  /**
   * Set parameter velocities into the refinement model based on the provided refinement mode.
   * The method updates coordinates, B-factors, and occupancies for atoms, depending on what is
   * specified by the refinement mode.
   *
   * @param x An array of doubles containing the new parameter velocities. The values should be ordered as
   *          required by the refinement mode: first coordinates, then B-factors, and finally occupancies,
   *          if applicable.
   */
  public void setVelocity(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.setVelocity(x);
    }
  }

  /**
   * Retrieves acceleration values for all refined parameters
   * in the list using the provided array.
   *
   * @param x an array of doubles to store acceleration values for each parameter
   */
  public void getAcceleration(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.getAcceleration(x);
    }
  }

  /**
   * Sets the acceleration values for all parameters in the parameters list.
   *
   * @param x an array of double values representing the acceleration to be set for each parameter
   */
  public void setAcceleration(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.setAcceleration(x);
    }
  }

  /**
   * Retrieves previous acceleration values for all refined parameters
   * in the list using the provided array.
   *
   * @param x an array of doubles to store previous acceleration values for each parameter
   */
  public void getPreviousAcceleration(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.getPreviousAcceleration(x);
    }
  }

  /**
   * Sets the previous acceleration values for all parameters in the parameters list.
   *
   * @param x an array of double values representing the previous acceleration to be set for each parameter
   */
  public void setPreviousAcceleration(double[] x) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.setPreviousAcceleration(x);
    }
  }


  /**
   * Populates the provided array with mass values retrieved from all refined parameters.
   *
   * @param mass an array where the mass values will be stored.
   */
  public void getMass(double[] mass) {
    for (RefinedParameter parameter : allParametersList) {
      if (parameter instanceof RefinedCoordinates) {
        parameter.getMass(mass, 5.0);
      } else if (parameter instanceof RefinedBFactor) {
        parameter.getMass(mass, bMass);
      } else if (parameter instanceof RefinedOccupancy) {
        parameter.getMass(mass, occMass);
      }
    }
  }

  /**
   * Loads the optimization scale factors into all refined parameters.
   *
   * @param optimizationScaling an array of scale factors to be applied to each refined parameter
   */
  public void loadOptimizationScaling(double[] optimizationScaling) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.setOptimizationScaling(optimizationScaling);
    }
  }

  /**
   * Zero out the gradient for all atoms being refined.
   */
  public void zeroGradient() {
    for (RefinedParameter parameter : allParametersList) {
      parameter.zeroGradient();
    }
  }

  /**
   * Loads the gradient values into the respective refined parameters based on the current refinement mode.
   *
   * @param gradient an array of double values representing the gradient to be loaded.
   */
  public void getGradient(double[] gradient) {
    for (RefinedParameter parameter : allParametersList) {
      parameter.getGradient(gradient);
    }
  }

  /**
   * Creates a coordinate refinement model by identifying scattering and active atoms
   * from the provided molecular assemblies and defining coordinate constraints for
   * specific atoms.
   *
   * @param scatteringList The list of atoms in the scattering model.
   * @return A map where the keys represent atoms with constrained coordinates, and the values
   * correspond to the indices of their paired atoms in the active atom array.
   */
  private Map<Atom, RefinedCoordinates> createCoordinateModel(List<Atom> scatteringList) {
    logger.fine("\n Creating Coordinate Refinement Model\n");

    MolecularAssembly rootAssembly = molecularAssemblies[0];
    // The keys are references to Atom instances.
    Map<Atom, RefinedCoordinates> coordinateMap = new IdentityHashMap<>();

    // Loop over all atoms in the root MolecularAssembly.
    Atom[] atomList = rootAssembly.getAtomArray();
    for (Atom a : atomList) {
      // All atoms from the root molecular assembly are added to the scattering list.
      scatteringList.add(a);
      if (a.isActive()) {
        // Active atoms can have their coordinates refined.
        a.setXrayCoordIndex(coordinateAtomList.size());
        coordinateAtomList.add(a);
        coordinateMap.put(a, new RefinedCoordinates(a));
        logger.fine(" Active: " + a);
      }
    }

    // Add scattering atoms from other topologies and create coordinate constraints for deuterium atoms.
    for (int i = 1; i < molecularAssemblies.length; i++) {
      MolecularAssembly molecularAssembly = molecularAssemblies[i];
      atomList = molecularAssembly.getAtomArray();
      for (Atom a : atomList) {
        Character altLoc = a.getAltLoc();
        Atom rootAtom = rootAssembly.findAtom(a, false);
        Atom deuteriumMatch = rootAssembly.findAtom(a, true);
        if (rootAtom != null && rootAtom.getAltLoc().equals(altLoc)) {
          // This atom is identical to an atom in Conformation A and does not scatter.
          if (rootAtom.isActive()) {
            // The coordinates are constrained to match those of the atom in the root topology.
            RefinedCoordinates refinedCoordinates = coordinateMap.get(rootAtom);
            refinedCoordinates.addConstrainedAtom(a);
            a.setXrayCoordIndex(rootAtom.getXrayCoordIndex());
          } else {
            // Ensure paired atoms active status concords.
            a.setActive(false);
          }
        } else if (deuteriumMatch != null) {
          // This is an H/D pair.
          scatteringList.add(a);
          if (deuteriumMatch.isActive()) {
            RefinedCoordinates refinedCoordinates = coordinateMap.get(deuteriumMatch);
            refinedCoordinates.addConstrainedAtomThatScatters(a);
            a.setXrayCoordIndex(deuteriumMatch.getXrayCoordIndex());
          } else {
            // Ensure paired atoms active status concords.
            a.setActive(false);
          }
        } else {
          // This atom is part of an alternate conformation not found in the root topology.
          scatteringList.add(a);
          if (a.isActive()) {
            a.setXrayCoordIndex(coordinateAtomList.size());
            coordinateAtomList.add(a);
            coordinateMap.put(a, new RefinedCoordinates(a));
          }
        }
      }
    }

    return coordinateMap;
  }

  /**
   * Creates a b-factor refinement model by identifying active atoms whose b-factors
   * should be refined based on the provided molecular assemblies and
   * defining b-factor constraints.
   *
   * @return A map where the keys represent atoms with constrained b-factors, and the values
   * correspond to the indices of their paired atoms in the b-factor array.
   */
  private Map<Atom, RefinedBFactor> createBFactorModel() {
    logger.fine("\n Creating B-Factor Refinement Model\n");
    MolecularAssembly rootAssembly = molecularAssemblies[0];
    Map<Atom, RefinedBFactor> bFactorMap = new IdentityHashMap<>();

    if (byResidue) {
      // Collect residue b-factors for the root topology.
      Polymer[] polymers = rootAssembly.getChains();
      for (Polymer polymer : polymers) {
        List<Residue> residues = polymer.getResidues();
        Atom heavyAtom = null;
        RefinedBFactor currentRefinedBFactor = null;
        for (int j = 0; j < residues.size(); j++) {
          Residue residue = residues.get(j);
          if (j % nResiduePerBFactor == 0 || heavyAtom == null) {
            heavyAtom = residue.getFirstActiveHeavyAtom();
            if (heavyAtom == null) {
              // This residue is inactive.
              continue;
            }
            bFactorAtomList.add(heavyAtom);
            currentRefinedBFactor = new RefinedBFactor(heavyAtom);
            bFactorMap.put(heavyAtom, currentRefinedBFactor);
          }
          // The rest of the atoms in this residue are constrained.
          for (Atom a : residue.getAtomList()) {
            if (a != heavyAtom) {
              currentRefinedBFactor.addConstrainedAtomThatScatters(a);
            }
          }
        }
      }
      List<MSNode> molecules = rootAssembly.getNodeList(true);
      for (MSNode m : molecules) {
        Atom heavyAtom = m.getFirstActiveHeavyAtom();
        if (heavyAtom == null) {
          // This molecule is inactive.
          continue;
        }
        bFactorAtomList.add(heavyAtom);
        RefinedBFactor currentRefinedBFactor = new RefinedBFactor(heavyAtom);
        bFactorMap.put(heavyAtom, currentRefinedBFactor);
        // The rest of the atoms in this residue are constrained.
        for (Atom a : m.getAtomList()) {
          if (a != heavyAtom) {
            currentRefinedBFactor.addConstrainedAtomThatScatters(a);
          }
        }
      }

      // Residue-based b-factors should generally not be used with alternative conformations.
      for (int i = 1; i < molecularAssemblies.length; i++) {
        MolecularAssembly molecularAssembly = molecularAssemblies[i];
        polymers = molecularAssembly.getChains();
        for (int j = 0; j < polymers.length; j++) {
          Polymer polymer = polymers[j];
          List<Residue> residues = polymer.getResidues();
          for (int k = 0; k < residues.size(); k++) {
            Residue residue = residues.get(k);
            Residue rootResidue = rootAssembly.getResidue(j, k);
            if (rootResidue == null) {
              logger.severe(format(" Residue %s not found in the root conformation.", residue));
              return null;
            }
            Atom rootAtom = rootResidue.getFirstActiveHeavyAtom();
            if (rootAtom == null) {
              // This residue is not active.
              continue;
            }
            RefinedBFactor refinedBFactor = bFactorMap.get(rootAtom);
            // The B-factors in this residue are constrained.
            for (Atom a : residue.getAtomList()) {
              Character altLoc = a.getAltLoc();
              if (!altLoc.equals(rootAtom.getAltLoc())) {
                refinedBFactor.addConstrainedAtomThatScatters(a);
              } else {
                refinedBFactor.addConstrainedAtom(a);
              }
            }
          }
        }

        molecules = molecularAssemblies[i].getNodeList(true);
        for (MSNode m : molecules) {
          Atom heavyAtom = m.getFirstActiveHeavyAtom();
          if (heavyAtom == null) {
            // This molecule is inactive.
            continue;
          }
          Character altLoc = heavyAtom.getAltLoc();
          if (!altLoc.equals(' ') && !altLoc.equals('A')) {
            bFactorAtomList.add(heavyAtom);
            RefinedBFactor refinedBFactor = new RefinedBFactor(heavyAtom);
            bFactorMap.put(heavyAtom, refinedBFactor);
            // The rest of the atoms in this molecule are constrained.
            for (Atom a : m.getAtomList()) {
              if (a != heavyAtom) {
                refinedBFactor.addConstrainedAtomThatScatters(a);
              }
            }
          }
        }
      }
    } else if (ridingHydrogen) {
      // Hydrogen will use the B-Factor of their heavy atom.
      // Add heavy atoms and deuterium get unique b-factors
      for (Atom atom : coordinateAtomList) {
        if (atom.isHydrogen() && !atom.isDeuterium()) {
          continue;
        }
        bFactorAtomList.add(atom);
        RefinedBFactor refinedBFactor = new RefinedBFactor(atom);
        bFactorMap.put(atom, refinedBFactor);
        // Non-scattering atoms constrained to this B-Factor
        RefinedCoordinates refinedCoords = refinedCoordinates.get(atom);
        for (Atom a : refinedCoords.constrainedAtoms) {
          refinedBFactor.addConstrainedAtom(a);
        }
        // Scattering atoms constrained to this B-factor
        for (Atom a : refinedCoords.constrainedAtomsThatScatter) {
          refinedBFactor.addConstrainedAtomThatScatters(a);
        }
      }
      // Add hydrogen constrained to their heavy atom.
      for (Atom atom : coordinateAtomList) {
        if (atom.isHydrogen() && !atom.isDeuterium()) {
          Atom heavy = atom.getBonds().getFirst().get1_2(atom);
          if (!heavy.isActive()) {
            // If the heavy atom is not active, then do not refine this hydrogen b-factor.
            continue;
          }
          if (bFactorMap.containsKey(heavy)) {
            RefinedBFactor refinedBFactor = bFactorMap.get(heavy);
            refinedBFactor.addConstrainedAtomThatScatters(atom);
            continue;
          }
          logger.info(" Could not locate a heavy atom B-factor for: " + atom);
        }
      }
    } else {
      // No special constraints.
      // The b-factors of all active atoms are refined.
      for (Atom atom : coordinateAtomList) {
        bFactorAtomList.add(atom);
        RefinedBFactor refinedBFactor = new RefinedBFactor(atom);
        bFactorMap.put(atom, refinedBFactor);
        // Non-scattering Atoms constrained to this BFactor
        RefinedCoordinates refinedCoords = refinedCoordinates.get(atom);
        for (Atom a : refinedCoords.constrainedAtoms) {
          refinedBFactor.addConstrainedAtom(a);
        }
        // Scattering Atoms constrained to this B-factor
        for (Atom a : refinedCoords.constrainedAtomsThatScatter) {
          refinedBFactor.addConstrainedAtomThatScatters(a);
        }
      }
    }

    return bFactorMap;
  }

  /**
   * Create a list of Bond restraints for B-factors being refined.
   */
  private void collectBFactorRestraints() {
    bFactorRestraints.clear();
    MolecularAssembly rootAssembly = molecularAssemblies[0];
    // Add each bond from the root assembly if at least one atom of the bond is active.
    List<Bond> rootBonds = rootAssembly.getBondList();
    for (Bond bond : rootBonds) {
      Atom a1 = bond.getAtom(0);
      Atom a2 = bond.getAtom(1);
      if (!a1.isActive() && !a2.isActive()) {
        continue;
      }
      bFactorRestraints.add(new Atom[]{a1, a2});
    }

    // Add each bond from alternate conformers is included if at least one atom is active and
    // at least one atom is from the alternate location.
    for (int i = 1; i < molecularAssemblies.length; i++) {
      MolecularAssembly molecularAssembly = molecularAssemblies[i];
      Character altLoc = molecularAssembly.getAlternateLocation();
      List<Bond> bonds = molecularAssembly.getBondList();
      for (Bond bond : bonds) {
        Atom a1 = bond.getAtom(0);
        Atom a2 = bond.getAtom(1);
        // One atom must be active.
        if (!a1.isActive() && !a2.isActive()) {
          continue;
        }
        // Both atoms are part of the alternate conformer.
        if (a1.getAltLoc().equals(altLoc) && a2.getAltLoc().equals(altLoc)) {
          bFactorRestraints.add(new Atom[]{a1, a2});
        } else if (a1.getAltLoc().equals(altLoc) && !a2.getAltLoc().equals(altLoc)) {
          // Atom 1 is part of the alternate conformer.
          a2 = rootAssembly.findAtom(a2);
          if (a2 != null) {
            bFactorRestraints.add(new Atom[]{a1, a2});
          }
        } else if (!a1.getAltLoc().equals(altLoc) && a2.getAltLoc().equals(altLoc)) {
          // Atom 2 is part of the alternate conformer.
          a1 = rootAssembly.findAtom(a1);
          if (a1 != null) {
            bFactorRestraints.add(new Atom[]{a1, a2});
          }
        }
      }
    }
  }

  /**
   * Creates an occupancy model by identifying and refining atoms within residues or molecules
   * that have alternate conformations or less-than-full occupancies. The method looks through
   * residues and molecules in molecular assemblies, evaluates their constituent atoms, and groups
   * those with alternate conformers into refined occupancy objects for further analysis or refinement.
   * <p>
   * Alternate residues and molecules are tracked separately in internal structures, and any atoms
   * that scatter due to constraints are also linked appropriately.
   *
   * @return A map of atoms to their corresponding refined occupancy objects. These objects
   * encapsulate information about atoms with alternate conformations or partial occupancies
   * and their constrained scattering atoms.
   */
  private Map<Atom, RefinedOccupancy> createOccupancyModel() {
    logger.fine("\n Creating Occupancy Refinement Model\n");
    Map<Atom, RefinedOccupancy> refinedOccupancies = new IdentityHashMap<>();

    boolean refineDeuterium = false;
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      if (molecularAssembly.hasDeuterium()) {
        refineDeuterium = true;
        break;
      }
    }

    if (refineDeuterium) {
      // Find polymer hydrogen / deuterium with occupancy less than one.
      MolecularAssembly rootAssembly = molecularAssemblies[0];
      Polymer[] polymers = rootAssembly.getChains();
      if (polymers != null) {
        for (Polymer polymer : polymers) {
          List<Residue> residues = polymer.getResidues();
          for (Residue residue : residues) {
            List<Atom> atoms = residue.getAtomList();
            for (Atom atom : atoms) {
              if (atom.isActive() && atom.isHydrogen()) {
                double occupancy = atom.getOccupancy();
                if (occupancy < 1.0) {
                  Atom heavy = atom.getBonds().getFirst().get1_2(atom);
                  if (refinedOccupancies.containsKey(heavy)) {
                    RefinedOccupancy refinedOccupancy = refinedOccupancies.get(heavy);
                    refinedOccupancy.addConstrainedAtomThatScatters(atom);
                  } else {
                    RefinedOccupancy refinedOccupancy = new RefinedOccupancy(atom);
                    refinedOccupancies.put(heavy, refinedOccupancy);
                    occupancyAtomList.add(heavy);
                  }
                }
              }
            }
          }
        }
      }
      // Find matching hydrogen / deuterium in the 2nd Assembly.
      if (molecularAssemblies.length > 1) {
        MolecularAssembly molecularAssembly = molecularAssemblies[1];
        for (Atom heavy : occupancyAtomList) {
          RefinedOccupancy refinedOccupancy = refinedOccupancies.get(heavy);
          // Collect the H/D atoms from conformation A.
          List<Atom> atoms = new ArrayList<>(refinedOccupancy.constrainedAtomsThatScatter);
          atoms.add(refinedOccupancy.atom);
          for (Atom atom : atoms) {
            // Find the matching atom in conformation B.
            Atom match = molecularAssembly.findAtom(atom, true);
            if (match != null) {
              double o1 = atom.getOccupancy();
              double o2 = match.getOccupancy();
              if (resetHDOccupancy) {
                logger.info(" Reset Occupancy for H/D Pair to 0.5/0.5:");
                atom.setOccupancy(0.5);
                match.setOccupancy(0.5);
                logger.info(format(" %s: %6.3f", atom, atom.getOccupancy()));
                logger.info(format(" %s: %6.3f", match, match.getOccupancy()));
              } else if (o1 + o2 != 1.0) {
                logger.info(" Occupancy Sum for H/D Pair is not 1.0:");
                logger.info(format(" %s: %6.3f", atom, atom.getOccupancy()));
                logger.info(format(" %s: %6.3f", match, match.getOccupancy()));
                double delta = (1.0 - o1 - o2) / 2.0;
                atom.setOccupancy(o1 + delta);
                match.setOccupancy(o2 + delta);
                logger.info(" Occupancy Sum for H/D Pair adjusted to 1.0:");
                logger.info(format(" %s: %6.3f", atom, atom.getOccupancy()));
                logger.info(format(" %s: %6.3f", match, match.getOccupancy()));
              }
              refinedOccupancy.addConstrainedAtomThatScattersComplement(match);
            }
          }
        }
      }
    } else {
      // Find Residues with alternate conformers.
      Polymer[] polymers = molecularAssemblies[0].getChains();
      if (polymers != null && polymers.length > 0) {
        for (int i = 0; i < polymers.length; i++) {
          List<Residue> residues = polymers[i].getResidues();
          for (int j = 0; j < residues.size(); j++) {
            List<Residue> list = getResidueConformers(i, j);
            if (list != null && !list.isEmpty()) {
              altResidues.add(list);
              for (Residue residue : list) {
                List<Atom> atomList = residue.getAtomList();
                RefinedOccupancy refinedOccupancy = null;
                Atom refinedAtom = null;
                for (Atom a : atomList) {
                  Character altLoc = a.getAltLoc();
                  double occupancy = a.getOccupancy();
                  if (!altLoc.equals(' ') || occupancy < 1.0) {
                    occupancyAtomList.add(a);
                    refinedOccupancy = new RefinedOccupancy(a);
                    refinedOccupancies.put(a, refinedOccupancy);
                    // logger.info(" Occupancy: " + a);
                    refinedAtom = a;
                    break;
                  }
                }
                if (refinedAtom != null) {
                  for (Atom a : atomList) {
                    if (a == refinedAtom) {
                      continue;
                    }
                    Character altLoc = a.getAltLoc();
                    double occupancy = a.getOccupancy();
                    if (!altLoc.equals(' ') || occupancy < 1.0) {
                      refinedOccupancy.addConstrainedAtomThatScatters(a);
                      // logger.info("  Constrained Occupancy: " + a);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Find Molecules with non-zero occupancies.
    if (refineMolOcc) {
      List<MSNode> molecules = molecularAssemblies[0].getMolecules();
      if (molecules != null && !molecules.isEmpty()) {
        for (int i = 0; i < molecules.size(); i++) {
          List<Molecule> list = getMoleculeConformers(i);
          if (list != null && !list.isEmpty()) {
            altMolecules.add(list);
            for (Molecule molecule : list) {
              List<Atom> atomList = molecule.getAtomList();
              RefinedOccupancy refinedOccupancy = null;
              Atom refinedAtom = null;
              for (Atom a : atomList) {
                Character altLoc = a.getAltLoc();
                double occupancy = a.getOccupancy();
                if (!altLoc.equals(' ') || occupancy < 1.0) {
                  occupancyAtomList.add(a);
                  refinedOccupancy = new RefinedOccupancy(a);
                  refinedOccupancies.put(a, refinedOccupancy);
                  logger.info(" Refined Occupancy: " + a);
                  refinedAtom = a;
                  break;
                }
              }
              if (refinedAtom != null) {
                for (Atom a : atomList) {
                  if (a == refinedAtom) {
                    continue;
                  }
                  Character altLoc = a.getAltLoc();
                  double occupancy = a.getOccupancy();
                  if (!altLoc.equals(' ') || occupancy < 1.0) {
                    refinedOccupancy.addConstrainedAtomThatScatters(a);
                    logger.info("  Constrained Occupancy: " + a);
                  }
                }
              }
            }
          }
        }
      }
    }

    return refinedOccupancies;
  }

  /**
   * Find residue-based alternate conformers.
   *
   * @param polymerID Polymer ID.
   * @param resID     The residue ID.
   * @return The constrained residues.
   */
  private List<Residue> getResidueConformers(int polymerID, int resID) {
    if (molecularAssemblies.length < 2) {
      return null;
    }

    double totalOccupancy = 0.0;
    List<Residue> residues = new ArrayList<>();
    // Check the residue from Conformer A.
    Residue residue = molecularAssemblies[0].getResidue(polymerID, resID);
    for (Atom a : residue.getAtomList()) {
      if (!a.getUse()) {
        continue;
      }
      Character altLoc = a.getAltLoc();
      double occupancy = a.getOccupancy();
      if (!altLoc.equals(' ') || occupancy < 1.0) {
        // Include this residue.
        logger.fine(format(" %s %c %5.3f", residue, altLoc, occupancy));
        totalOccupancy = occupancy;
        residues.add(residue);
        break;
      }
    }
    // No altLoc found for this residue.
    if (residues.isEmpty()) {
      return null;
    }

    // Find this residue in the other conformers.
    int numConformers = molecularAssemblies.length;
    for (int i = 1; i < numConformers; i++) {
      residue = molecularAssemblies[i].getResidue(polymerID, resID);
      for (Atom a : residue.getAtomList()) {
        if (!a.getUse()) {
          continue;
        }
        Character altLoc = a.getAltLoc();
        if (!altLoc.equals(' ') && !altLoc.equals('A')) {
          double occupancy = a.getOccupancy();
          totalOccupancy += occupancy;
          // Include this residue.
          residues.add(residue);
          logger.fine(format(" %s %c %5.3f", residue, altLoc, occupancy));
          break;
        }
      }
    }

    logger.fine("  Total occupancy: " + totalOccupancy);
    return residues;
  }

  /**
   * Find molecule-based alternate conformers.
   *
   * @param moleculeID The molecule ID.
   * @return The constrained molecules.
   */
  private List<Molecule> getMoleculeConformers(int moleculeID) {
    List<Molecule> molecules = new ArrayList<>();
    double totalOccupancy = 0.0;
    // Check the molecule from Conformer A.
    List<MSNode> molList = molecularAssemblies[0].getMolecules();
    if (molList != null && !molList.isEmpty()) {
      Molecule molecule = (Molecule) molList.get(moleculeID);
      for (Atom a : molecule.getAtomList()) {
        if (!a.getUse()) {
          continue;
        }
        Character altLoc = a.getAltLoc();
        double occupancy = a.getOccupancy();
        if (!altLoc.equals(' ') || occupancy < 1.0) {
          // Include this residue.
          totalOccupancy = occupancy;
          molecules.add(molecule);
          logger.fine(format(" %s %c %5.3f", molecule, altLoc, occupancy));
          break;
        }
      }
    }

    // No altLoc found for this residue.
    if (molecules.isEmpty()) {
      return null;
    }

    // Find this molecule in the other conformers.
    int numConformers = molecularAssemblies.length;
    for (int i = 1; i < numConformers; i++) {
      molList = molecularAssemblies[i].getMolecules();
      Molecule molecule = (Molecule) molList.get(moleculeID);
      for (Atom a : molecule.getAtomList()) {
        if (!a.getUse()) {
          continue;
        }
        Character altLoc = a.getAltLoc();
        if (!altLoc.equals(' ') && !altLoc.equals('A')) {
          // Include this residue.
          double occupancy = a.getOccupancy();
          totalOccupancy += occupancy;
          molecules.add(molecule);
          logger.fine(format(" %s %c %5.3f", molecule, altLoc, occupancy));
          break;
        }
      }
    }

    logger.fine("  Total occupancy: " + totalOccupancy);
    return molecules;
  }

  /**
   * Sets the parameter indices for the current refinement mode.
   * <p>
   * This method iterates over collections of refined parameters grouped by
   * coordinates, B-factors and occupancies to assign an incremental index to each parameter.
   * The indices are stored within the respective refined parameter objects.
   */
  private void setParameterIndices() {
    int index = 0;
    for (RefinedParameter parameter : allParametersList) {
      parameter.setIndex(index);
      index += parameter.getNumberOfParameters();
    }
  }

  /**
   * Adjusts the active status of atoms within each molecular assembly.
   * For riding hydrogen b-factor refinement, hydrogen atoms adopt the active status of
   * their associated heavy atom. For group-based b-factor refinement, if any heavy atom
   * in a residue or molecule is active, all its atoms are made active. Otherwise, all
   * atoms are inactive.
   */
  private void regularizeActiveAtoms() {
    for (MolecularAssembly molecularAssembly : molecularAssemblies) {
      // For riding hydrogen b-factor refinement, hydrogen atoms will adopt the active status of their
      // heavy atom.
      if (ridingHydrogen) {
        for (Atom atom : molecularAssembly.getAtomList()) {
          if (atom.isHydrogen()) {
            Atom other = atom.getBonds().getFirst().get1_2(atom);
            atom.setActive(other.isActive());
          }
        }
      }
      // For group-based b-factor refinement, if one heavy atom of a residue or molecule is active, then
      // all atoms will be active.
      if (byResidue) {
        List<MSNode> nodeList = molecularAssembly.getNodeList();
        for (MSNode node : nodeList) {
          Atom activeAtom = node.getFirstActiveHeavyAtom();
          if (activeAtom != null) {
            for (Atom atom : node.getAtomList()) {
              atom.setActive(true);
            }
          } else {
            // No active heavy atom -- set hydrogen to inactive.
            for (Atom atom : node.getAtomList()) {
              atom.setActive(false);
            }
          }
        }
      }
    }
  }

  /**
   * Adds anisotropic B-factors (ANISOU) to active heavy atoms in the provided molecular
   * assemblies that currently lack them.
   * Anisotropic B-factors for newly added atoms are initialized isotropically.
   */
  private void addAnisotropicBFactors() {
    if (addAnisou) {
      for (MolecularAssembly molecularAssembly : molecularAssemblies) {
        int count = 0;
        List<Atom> atomList = molecularAssembly.getAtomList();
        for (Atom a : atomList) {
          // Add an Anisou to each active heavy atom that lacks one.
          if (a.isHeavy() && a.isActive() && a.getAnisou(null) == null) {
            double[] anisou = new double[6];
            double u = b2u(a.getTempFactor());
            anisou[0] = u;
            anisou[1] = u;
            anisou[2] = u;
            anisou[3] = 0.0;
            anisou[4] = 0.0;
            anisou[5] = 0.0;
            a.setAnisou(anisou);
            count++;
          }
        }
        if (count > 0) {
          Character c = molecularAssembly.getAlternateLocation();
          logger.info(format(" %d anisotropic B-factors were added to conformer %c.", count, c));
        }
      }
    }
  }

}
