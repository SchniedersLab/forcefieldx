// ******************************************************************************
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
// ******************************************************************************
package ffx.algorithms.misc;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.GeneralizedKirkwood;
import ffx.potential.nonbonded.GeneralizedKirkwood.NonPolar;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.utils.PotentialsFunctions;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * Evaluate binding affinity from an MD trajectory using MM-GKSA/MM-GBSA/MM-PBSA approximation.
 * Binding affinity is mean value of (complex - protein - ligand), often weighted by component.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class MMgksa {
  private static final Logger logger = Logger.getLogger(MMgksa.class.getName());
  private final MolecularAssembly mola;
  private final PotentialsFunctions functs;
  private final SystemFilter filter;
  private final ForceFieldEnergy energy;
  private final Atom[] proteinAtoms;
  private final Atom[] ligandAtoms;
  private final boolean docking;
  private final EnergyTerm gksaTerm;
  private Atom[] ignoreAtoms;
  private double elecWt = 1.0;
  private double solvWt = 1.0;
  private double vdwWt = 1.0;
  private boolean decompose = false;
  /** Only used if decompose is true; if so, prints out cavitation and dispersion components. */
  private NonPolar nonPolar = NonPolar.NONE;

  private EnergyTerm[] terms;

  /**
   * Constructs an MMgksa instance on a single assembly (not for binding).
   *
   * @param mola a {@link ffx.potential.MolecularAssembly} object.
   * @param functs a {@link ffx.potential.utils.PotentialsFunctions} object.
   * @param filter a {@link ffx.potential.parsers.SystemFilter} object.
   */
  public MMgksa(MolecularAssembly mola, PotentialsFunctions functs, SystemFilter filter) {
    this.mola = mola;
    this.functs = functs;
    this.filter = filter;
    this.energy = mola.getPotentialEnergy();
    proteinAtoms = new Atom[0];
    ligandAtoms = new Atom[0];
    docking = false;
    gksaTerm =
        new EnergyTerm(
            "MM-GKSA",
            () -> {
              double e = (energy.getElectrostaticEnergy() * elecWt);
              e += (energy.getSolvationEnergy() * solvWt);
              e += (energy.getVanDerWaalsEnergy() * vdwWt);
              return e;
            });
  }

  /**
   * Constructs an MMgksa instance for MM-GKSA binding scoring.
   *
   * @param mola a {@link ffx.potential.MolecularAssembly} object.
   * @param functs a {@link ffx.potential.utils.PotentialsFunctions} object.
   * @param filter a {@link ffx.potential.parsers.SystemFilter} object.
   * @param proteinAtoms Atoms for first binding partner
   * @param ligandAtoms Atoms for second binding partner
   */
  public MMgksa(
      MolecularAssembly mola,
      PotentialsFunctions functs,
      SystemFilter filter,
      Atom[] proteinAtoms,
      Atom[] ligandAtoms) {
    this.mola = mola;
    this.functs = functs;
    this.filter = filter;
    this.energy = mola.getPotentialEnergy();
    this.proteinAtoms = new Atom[proteinAtoms.length];
    System.arraycopy(proteinAtoms, 0, this.proteinAtoms, 0, proteinAtoms.length);
    this.ligandAtoms = new Atom[ligandAtoms.length];
    System.arraycopy(ligandAtoms, 0, this.ligandAtoms, 0, ligandAtoms.length);
    docking = true;
    gksaTerm =
        new EnergyTerm(
            "MM-GKSA",
            () -> {
              double e = (energy.getElectrostaticEnergy() * elecWt);
              e += (energy.getSolvationEnergy() * solvWt);
              e += (energy.getVanDerWaalsEnergy() * vdwWt);
              return e;
            });
  }

  /**
   * Evaluates binding or average protein energy.
   *
   * @param frequency Evaluate energy every X frames
   */
  public void runMMgksa(int frequency) {
    runMMgksa(frequency, -1);
  }

  /**
   * runMMgksa.
   *
   * @param frequency a int.
   * @param maxEvals a int.
   */
  public void runMMgksa(int frequency, int maxEvals) {
    setEnergyTerms();
    int nTerms = terms.length;
    double[] meanE = new double[nTerms];
    double meanGKSA = 0.0;
    Arrays.fill(meanE, 0.0);
    int nEvals = 0;
    int counter = 0;
    do {
      if (counter++ % frequency != 0) {
        continue;
      }
      ++nEvals;
      double[] totE = currentEnergies();
      double totGKSA = gksaTerm.en();
      if (docking) {
        for (Atom atom : ligandAtoms) {
          atom.setUse(false);
        }
        double[] protE = currentEnergies();
        double protGKSA = gksaTerm.en();
        for (Atom atom : ligandAtoms) {
          atom.setUse(true);
        }

        for (Atom atom : proteinAtoms) {
          atom.setUse(false);
        }
        double[] ligE = currentEnergies();
        double ligGKSA = gksaTerm.en();
        for (Atom atom : proteinAtoms) {
          atom.setUse(true);
        }

        double[] snapE = new double[nTerms];
        for (int i = 0; i < nTerms; i++) {
          snapE[i] = totE[i] - (protE[i] + ligE[i]);
        }
        double snapGKSA = totGKSA - (protGKSA + ligGKSA);

        for (int i = 0; i < nTerms; i++) {
          meanE[i] += (snapE[i] - meanE[i]) / nEvals;
        }
        meanGKSA += (snapGKSA - meanGKSA) / nEvals;

        logger.info(
            String.format(" %10d frames read, %10d frames " + "evaluated", counter, nEvals));
        logger.info(formatHeader());
        logger.info(formatEnergy("Running Mean", meanGKSA, meanE));
        logger.info(formatEnergy("Snapshot", snapGKSA, snapE));
        logger.info(formatEnergy("Complex", totGKSA, totE));
        logger.info(formatEnergy("Protein", protGKSA, protE));
        logger.info(formatEnergy("Ligand", ligGKSA, ligE));
      }
    } while (filter.readNext() && (maxEvals < 0 || nEvals < maxEvals));

    logger.info(
        String.format(
            " MM-GKSA evaluation complete, %10d frames " + "read, %11.6f kcal/mol mean energy",
            nEvals, meanGKSA));
    filter.closeReader();
  }

  /**
   * Setter for the field <code>decompose</code>.
   *
   * @param decompose a boolean.
   */
  public void setDecompose(boolean decompose) {
    this.decompose = decompose;
    try {
      String cavModel = mola.getForceField().getString("CAVMODEL", "CAV_DISP").toUpperCase();
      nonPolar = GeneralizedKirkwood.getNonPolarModel(cavModel);
    } catch (Exception ex) {
      nonPolar = NonPolar.NONE;
    }
  }

  /**
   * Assigns weight to electrostatics.
   *
   * @param eWt a double.
   */
  public void setElectrostaticsWeight(double eWt) {
    this.elecWt = eWt;
  }

  /**
   * Sets atoms to be ignored entirely.
   *
   * @param ignored an array of {@link ffx.potential.bonded.Atom} objects.
   */
  public void setIgnoredAtoms(Atom[] ignored) {
    ignoreAtoms = new Atom[ignored.length];
    System.arraycopy(ignored, 0, ignoreAtoms, 0, ignored.length);
    for (Atom atom : ignoreAtoms) {
      atom.setUse(false);
    }
  }

  /**
   * Assigns weight to solvation.
   *
   * @param sWt a double.
   */
  public void setSolvationWeight(double sWt) {
    this.solvWt = sWt;
  }

  /**
   * Assigns weight to van der Waals energy.
   *
   * @param vWt a double.
   */
  public void setVdwWeight(double vWt) {
    this.vdwWt = vWt;
  }

  /**
   * Sets the appropriate solvation terms.
   *
   * @return
   */
  private List<EnergyTerm> getSolvationTerms() {
    List<EnergyTerm> termList = new ArrayList<>();
    EnergyTerm solvationTerm =
        new EnergyTerm("Solvation", energy::getSolvationEnergy); // May not be used.
    if (decompose) {
      EnergyTerm cavTerm = new EnergyTerm("Cavitation", energy::getCavitationEnergy);
      EnergyTerm dispTerm = new EnergyTerm("Dispersion", energy::getDispersionEnergy);
      switch (nonPolar) {
        case CAV_DISP:
          getE elecSolv =
              () -> {
                double e = energy.getSolvationEnergy();
                e -= energy.getCavitationEnergy();
                e -= energy.getDispersionEnergy();
                return e;
              };
          termList.add(new EnergyTerm("Elec. solv", elecSolv));
          termList.add(cavTerm);
          termList.add(dispTerm);
          break;
        case CAV:
          elecSolv = () -> (energy.getSolvationEnergy() - energy.getCavitationEnergy());
          termList.add(new EnergyTerm("Elec. solv", elecSolv));
          termList.add(cavTerm);
          break;
        case BORN_CAV_DISP:
          elecSolv = () -> (energy.getSolvationEnergy() - energy.getDispersionEnergy());
          termList.add(new EnergyTerm("Elec. solv", elecSolv));
          termList.add(dispTerm);
          break;
        default:
          termList.add(solvationTerm);
          break;
      }
    } else {
      termList.add(solvationTerm);
    }
    return termList;
  }

  /** Autodetects which energy terms are to be printed. */
  private void setEnergyTerms() {
    List<EnergyTerm> termSet = new ArrayList<>();
    // getE tot = this::weightEnergy;
    // termSet.add(new EnergyTerm("Total", tot));

    if (elecWt != 0.0) {
      if (decompose) {
        termSet.add(new EnergyTerm("Multipoles", energy::getPermanentMultipoleEnergy));
        termSet.add(new EnergyTerm("Polarization", energy::getPolarizationEnergy));
      } else {
        termSet.add(new EnergyTerm("Electrostatics", energy::getElectrostaticEnergy));
      }
    }

    if (solvWt != 0.0) {
      termSet.addAll(getSolvationTerms());
    }

    if (vdwWt != 0.0) {
      termSet.add(new EnergyTerm("van der Waals", energy::getVanDerWaalsEnergy));
    }
    terms = new EnergyTerm[termSet.size()];
    termSet.toArray(terms);
  }

  /**
   * Evaluates and decomposes energy of the system.
   *
   * @return MMGKSA energy components.
   */
  private double[] currentEnergies() {
    functs.energy(mola);
    int nTerms = terms.length;
    double[] energies = new double[nTerms];
    for (int i = 0; i < nTerms; i++) {
      energies[i] = terms[i].en();
    }
    return energies;
  }

  private String formatHeader() {
    StringBuilder sb = new StringBuilder(String.format("%15s%15s", " ", "GKSA"));
    for (EnergyTerm term : terms) {
      sb.append(String.format("%15s", term.toString()));
    }
    return sb.toString();
  }

  private String formatEnergy(String title, double gksa, double[] energies) {
    StringBuilder sb = new StringBuilder(String.format("%-15s", title));
    sb.append(String.format("  %13.5f", gksa));
    for (double e : energies) {
      sb.append(String.format("  %13.5f", e));
    }
    return sb.toString();
  }

  /** Functional interface for above. */
  private interface getE {
    double getEnergy();
  }

  /** Has a name, and some method to grab an energy. */
  private static class EnergyTerm {

    private final String name;
    private final getE method;

    public EnergyTerm(String name, getE method) {
      this.name = name;
      this.method = method;
    }

    public double en() {
      // return method.getEnergy(energy);
      return method.getEnergy();
    }

    public String getName() {
      return name;
    }

    @Override
    public String toString() {
      return getName();
    }
  }
}
