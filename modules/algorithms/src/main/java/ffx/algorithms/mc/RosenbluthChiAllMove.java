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
package ffx.algorithms.mc;

import static ffx.utilities.Constants.R;
import static java.lang.String.format;
import static java.lang.System.arraycopy;

import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueState;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.bonded.Torsion;
import ffx.potential.parameters.AtomType;
import ffx.potential.parsers.PDBFilter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Logger;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.util.FastMath;

/**
 * Represents a Boltzmann-drawn spin of all residue torsions for use with RosenbluthCBMC
 * (Configurational-Bias Monte Carlo).
 *
 * <p>Biases each torsion by drawing a test set from a Boltzmann distribution on torsion energy.
 * Selects from amongst each test set on Boltzmann weight of remaining energy.
 *
 * <p>Calculates Rosenbluth factors along the way; acceptance criterion = Wn/Wo.
 *
 * <p>Note: Contains much of the infrastructure necessary to generalize toward full
 * polymer-construction-type CBMC (i.e. changing angles and bonds as well). This implementation
 * leaves bonds/angles fixed and considers only torsion energy as dependent (Ubond in Frenkel/Smit
 * 13.2).
 *
 * @author Stephen D. LuCore
 */
public class RosenbluthChiAllMove implements MCMove {

  private static final Logger logger = Logger.getLogger(RosenbluthChiAllMove.class.getName());
  private static int numAccepted = 0;
  /** MolecularAssembly to operate on. */
  private final MolecularAssembly molecularAssembly;
  /** Target residue. */
  private final Residue target;
  /** Original residue state. */
  private final ResidueState origState;
  /** Force field energy to use. */
  private final ForceFieldEnergy forceFieldEnergy;

  private final double beta;
  private final int testSetSize;
  private final ThreadLocalRandom rand = ThreadLocalRandom.current();
  private final StringBuilder report = new StringBuilder();
  private final boolean[] doChi = new boolean[4];
  private final MODE mode;
  private final boolean torsionSampling;
  private final double CATASTROPHE_THRESHOLD = -10000;
  /** Convert nanoseconds to milliseconds. */
  private final double NS_TO_MS = 0.000001;
  /** Final energy. */
  double finalEnergy = 0.0;
  /** Proposed rotamer move. */
  private Rotamer proposedMove;
  /** PDBFilter to write out results. */
  private PDBFilter writer;
  /** Write out snapshots. */
  private SnapshotWriter snapshotWriter = null;

  private double Wn = 0.0;
  private double Wo = 0.0;
  private int moveNumber;
  private double[] proposedChis = new double[4];
  private boolean accepted = false;
  /** Original energy. */
  private double origEnergy;
  /** Start time. */
  private long startTime;
  /** End time. */
  private long endTime;

  private boolean verbose;
  private boolean randInts;
  private boolean noSnaps;
  private boolean printTestSets;
  private boolean logTimings;
  private boolean verboseEnergies;

  /**
   * Constructor for RosenbluthChiAllMove.
   *
   * @param molecularAssembly a {@link ffx.potential.MolecularAssembly} object.
   * @param target a {@link ffx.potential.bonded.Residue} object.
   * @param testSetSize a int.
   * @param forceFieldEnergy a {@link ffx.potential.ForceFieldEnergy} object.
   * @param temperature a double.
   * @param writeSnapshots a boolean.
   * @param moveNumber a int.
   * @param verbose a boolean.
   */
  RosenbluthChiAllMove(
      MolecularAssembly molecularAssembly,
      Residue target,
      int testSetSize,
      ForceFieldEnergy forceFieldEnergy,
      double temperature,
      boolean writeSnapshots,
      int moveNumber,
      boolean verbose) {

    CompositeConfiguration properties = molecularAssembly.getProperties();
    if (properties.containsKey("cbmc-type")) {
      mode = MODE.valueOf(properties.getString("cbmc-type"));
    } else {
      logger.severe("CBMC: must specify a bias type.");
      mode = MODE.CHEAP;
    }

    this.startTime = System.nanoTime();
    this.molecularAssembly = molecularAssembly;
    this.target = target;
    this.testSetSize = testSetSize;
    this.forceFieldEnergy = forceFieldEnergy;
    this.beta = 1 / (R * temperature);
    this.moveNumber = moveNumber;
    this.verbose = verbose;
    origState = target.storeState();

    randInts = properties.getBoolean("cbmc-randInts", false);
    noSnaps = properties.getBoolean("cbmc-noSnaps", false);
    printTestSets = properties.getBoolean("cbmc-printTestSets", false);
    logTimings = properties.getBoolean("cbmc-logTimings", false);
    verboseEnergies = properties.getBoolean("cbmc-verboseEnergies", false);

    if (writeSnapshots) {
      snapshotWriter = new SnapshotWriter(molecularAssembly, false);
    }

    doChi[0] = properties.getBoolean("cbmc-doChi0", false);
    doChi[1] = properties.getBoolean("cbmc-doChi1", false);
    doChi[2] = properties.getBoolean("cbmc-doChi2", false);
    doChi[3] = properties.getBoolean("cbmc-doChi3", false);
    if (!doChi[0] && !doChi[1] && !doChi[2] && !doChi[3]) {
      doChi[0] = true;
      doChi[1] = true;
      doChi[2] = true;
      doChi[3] = true;
    }
    updateAll();
    origEnergy = totalEnergy();

    torsionSampling = properties.getBoolean("cbmc-torsionSampler", false);
    if (torsionSampling) {
      logger.info(" Torsion Sampler engaged!");
      HashMap<Integer, BackBondedList> map =
          createBackBondedMap(AminoAcid3.valueOf(target.getName()));
      List<Torsion> allTors = new ArrayList<>();
      for (int i = 0; i < map.size(); i++) {
        Torsion tors = map.get(i).torsion;
        allTors.add(tors);
      }
      torsionSampler(allTors);
      System.exit(0);
    }
    try {
      switch (mode) {
        case EXPENSIVE:
          engageExpensive();
          break;
        case CHEAP:
          engageCheap();
          break;
        case CTRL_ALL:
          engageControlAll();
          break;
        default:
          logger.severe("CBMC: Unknown biasing type requested.");
          break;
      }
    } catch (ArithmeticException ex) {
      target.revertState(origState);
      accepted = false;
    }
  }

  /**
   * Getter for the field <code>mode</code>.
   *
   * @return a {@link ffx.algorithms.mc.RosenbluthChiAllMove.MODE} object.
   */
  public MODE getMode() {
    return mode;
  }

  /**
   * getWn.
   *
   * @return a double.
   */
  public double getWn() {
    return Wn;
  }

  /**
   * {@inheritDoc}
   *
   * <p>Performs the move associated with this MCMove. Also updates chi values in associated Torsion
   * objects.
   */
  @Override
  public void move() {
    RotamerLibrary.applyRotamer(target, proposedMove);
    updateAll();
  }

  /**
   * {@inheritDoc}
   *
   * <p>Reverts the last applied move() call.
   */
  @Override
  public void revertMove() {
    target.revertState(origState);
    updateAll();
  }

  /** {@inheritDoc} */
  @Override
  public String toString() {
    return format(
        "Rosenbluth Rotamer Move:\n   Res:   %s\n   Rota: %s",
        target.toString(), proposedMove.toString());
  }

  private void write() {
    if (noSnaps) {
      return;
    }
    if (writer == null) {
      writer = new PDBFilter(molecularAssembly.getFile(), molecularAssembly, null, null);
    }
    String filename = molecularAssembly.getFile().getAbsolutePath();
    if (!filename.contains("_mc")) {
      filename = FilenameUtils.removeExtension(filename) + "_mc.pdb";
    }
    File file = new File(filename);
    writer.writeFile(file, false);
  }

  /**
   * getWo.
   *
   * @return a double.
   */
  double getWo() {
    return Wo;
  }

  /**
   * This version foregoes doing a full energy eval (uExt) on each member of each chi test set.
   * Instead, each member of the test set is a full set of chi angles, the COMBINATION of which is
   * drawn from the Boltzmann.
   */
  private TrialSet cheapTorsionSet(List<Torsion> allTors, int setSize, String snapSuffix) {
    if (printTestSets) {
      report.append("    TrialSet_Cheap (uDep uExt)\n");
    }
    TrialSet trialSet = new TrialSet(setSize);
    double[] origChi = RotamerLibrary.measureRotamer(target, false);
    int i = 0;
    while (i < setSize) {
      double[] newChi = new double[origChi.length];
      arraycopy(origChi, 0, newChi, 0, origChi.length);
      for (int k = 0; k < origChi.length; k++) {
        if (doChi[k]) {
          if (randInts) {
            newChi[k] = rand.nextInt(360) - 180;
          } else {
            newChi[k] = rand.nextDouble(360.0) - 180;
          }
        }
      }
      Rotamer newState = createRotamer(target, newChi);
      RotamerLibrary.applyRotamer(target, newState);
      double uTors = 0;
      for (int j = 0; j < allTors.size(); j++) {
        if (doChi[j]) {
          uTors += allTors.get(j).energy(false);
          double offset = 0;
          try {
            offset = TORSION_OFFSET_AMPRO13.valueOf(target.getName() + j).offset;
          } catch (IllegalArgumentException ex) {
            logger.warning(ex.getMessage());
          }
          uTors += offset;
        }
      }
      double criterion = FastMath.exp(-beta * uTors);
      double rng = rand.nextDouble();
      if (rng < criterion) {
        trialSet.theta[i] = 0.0; // this cheap version does all thetas at once
        trialSet.rotamer[i] = newState;
        trialSet.uDep[i] = uTors;
        trialSet.uExt[i] = totalEnergy() - uTors;
        i++;
        writeSnapshot(snapSuffix, true);
        if (printTestSets) {
          if (i < 5 || i > setSize - 1) {
            report.append(
                format(
                    "       %3s %d:      %5.2f\t%5.2f\n",
                    snapSuffix, i, trialSet.uDep[i - 1], trialSet.uExt[i - 1]));
          } else if (i == 5) {
            report.append("       ...\n");
          }
        }
      }
    }
    target.revertState(origState);
    updateAll();
    return trialSet;
  }

  /**
   * Uses the accept-reject method (F/S Algorithm46) to draw new chi values for the given torsion.
   */
  private TrialSet expensiveTorsionSet(Torsion tors, int chiIndex, int setSize, String snapSuffix) {
    report.append(format("    TrialSet for Chi%d\t\t(Theta uDep uExt)\n", chiIndex));
    TrialSet trialSet = new TrialSet(setSize);
    double[] origChi = RotamerLibrary.measureRotamer(target, false);
    int i = 0;
    while (i < setSize) {
      double theta = rand.nextDouble(360.0) - 180;
      double[] newChi = new double[origChi.length];
      arraycopy(origChi, 0, newChi, 0, origChi.length);
      newChi[chiIndex] = theta;
      Rotamer newState = createRotamer(target, newChi);
      RotamerLibrary.applyRotamer(target, newState);
      double uTors = tors.energy(false);
      double offset = 0;
      try {
        offset = TORSION_OFFSET_AMPRO13.valueOf(target.getName() + chiIndex).offset;
      } catch (IllegalArgumentException ex) {
        logger.warning(ex.getMessage());
      }
      uTors += offset;
      double criterion = FastMath.exp(-beta * uTors);
      double rng = rand.nextDouble();
      report.append(format("    prop: %5.1f %.2g %.2g %.2g\n", theta, uTors, criterion, rng));
      if (rng < criterion) {
        report.append("    ^ Accepted!\n");
        writeSnapshot(snapSuffix, false);
        trialSet.theta[i] = theta;
        trialSet.rotamer[i] = newState;
        trialSet.uDep[i] = uTors;
        trialSet.uExt[i] = totalEnergy() - uTors; // Expensive!
        i++;
        writeSnapshot(snapSuffix, true);
        if (i < 4 || i > setSize - 1) {
          report.append(
              format(
                  "       %3s %d:      %5.2f\t%5.2f\t%5.2f\n",
                  snapSuffix, i, theta, trialSet.uDep[i - 1], trialSet.uExt[i - 1]));
        } else if (i == 4) {
          report.append("       ...\n");
        }
      }
    }
    target.revertState(origState);
    updateAll();
    return trialSet;
  }

  /**
   * So named, yet inadequate. Final stats on the ctrl vs bias algorithm for chis2,3 of LYS monomer
   * are: calc "4.23846e+07 / 22422" for the biased run calc "4.21623e+06 / 10000" for the control
   * run, where numerator = total timer; demon = accepted Possible accelerations are predicated on
   * the fact that the test above was performed using unbinned ie fully-continuous rotamers and
   * sampling was done with replacement. Rectifying either of these breaks balance, however, when
   * used with MD...
   */
  private void engageCheap() {
    report.append(format(" Rosenbluth CBMC Move: %4d  (%s)\n", moveNumber, target));

    AminoAcid3 name = AminoAcid3.valueOf(target.getName());
    double[] chi = RotamerLibrary.measureRotamer(target, false);
    HashMap<Integer, BackBondedList> map = createBackBondedMap(name);

    // For each chi, create a test set from Boltzmann distr on torsion energy (Ubond).
    // Select from among this set based on Boltzmann weight of REMAINING energy (Uext).
    // The Uext partition function of each test set (wn) becomes a factor of the overall Rosenbluth
    // (Wn).
    // ^ NOPE. Instead, Rosenbluth is going to be calculated just once for the whole combined set of
    // chis.
    List<Torsion> allTors = new ArrayList<>();
    for (int i = 0; i < chi.length; i++) {
      Torsion tors = map.get(i).torsion;
      allTors.add(tors);
    }
    TrialSet newTrialSet = cheapTorsionSet(allTors, testSetSize, "bkn");
    Wn = newTrialSet.sumExtBolt(); // yields uExt(1) + uExt(2) + ...
    if (Wn <= 0) {
      report.append("WARNING: Numerical instability in CMBC.");
      report.append("  Test set uExt values:  ");
      for (int i = 0; i < newTrialSet.uExt.length; i++) {
        report.append(format("%5.2f,  ", newTrialSet.uExt[i]));
      }
      report.append("  Discarding move.\n");
      target.revertState(origState);
      Wn = -1.0;
      Wo = 1000;
      logger.info(report.toString());
    }

    // Choose a proposal move from amongst this trial set (bn).
    double rng = rand.nextDouble(Wn);
    double running = 0.0;
    for (int j = 0; j < newTrialSet.uExt.length; j++) {
      double uExtBolt = FastMath.exp(-beta * newTrialSet.uExt[j]);
      running += uExtBolt;
      if (rng < running) {
        proposedMove = newTrialSet.rotamer[j];
        proposedChis = newTrialSet.theta;
        double prob = uExtBolt / Wn * 100;
        if (printTestSets) {
          report.append(
              format("       Chose %d   %7.4f\t  %4.1f%%\n", j + 1, newTrialSet.uExt[j], prob));
        }
        break;
      }
    }
    report.append(format("    Wn Total:  %g\n", Wn));

    // Reprise the above procedure for the old configuration.
    // The existing conformation forms first member of each test set (wo).
    double ouDep = 0.0;
    for (Torsion tors : allTors) {
      ouDep += tors.energy(false); // original-conf uDep
    }
    double ouExt = totalEnergy() - ouDep; // original-conf uExt
    double ouExtBolt = FastMath.exp(-beta * ouExt);
    if (printTestSets) {
      report.append(
          format("       %3s %d:  %9.5g  %9.5g  %9.5g\n", "bko", 0, ouDep, ouExt, ouExtBolt));
    }
    writeSnapshot("bko", true);
    TrialSet oldTrialSet = cheapTorsionSet(allTors, testSetSize - 1, "bko");
    Wo = ouExtBolt + oldTrialSet.sumExtBolt();

    report.append(format("    Wo Total:  %11.4g\n", Wo));
    report.append(format("    Wn/Wo:     %11.4g", Wn / Wo));

    RotamerLibrary.applyRotamer(target, proposedMove);
    proposedChis = RotamerLibrary.measureRotamer(target, false);
    finalEnergy = totalEnergy();
    if (finalEnergy < CATASTROPHE_THRESHOLD) {
      report.append("\nWARNING: Multipole catastrophe in CMBC.\n");
      report.append("  Discarding move.\n");
      target.revertState(origState);
      updateAll();
      Wn = -1.0;
      Wo = 1000;
      logger.info(report.toString());
    }

    double criterion = Math.min(1, Wn / Wo);
    rng = ThreadLocalRandom.current().nextDouble();
    report.append(format("    rng:    %5.2f\n", rng));
    if (rng < criterion) {
      accepted = true;
      numAccepted++;
      updateAll();
      write();
      report.append(format(" Accepted! %5d    NewEnergy: %.4f    Chi:", numAccepted, finalEnergy));
      for (int k = 0; k < proposedChis.length; k++) {
        report.append(format(" %7.2f", proposedChis[k]));
      }
      report.append(format("\n"));
    } else {
      accepted = false;
      target.revertState(origState);
      updateAll();
      report.append(format(" Denied.   %5d    OldEnergy: %.4f    Chi:", numAccepted, origEnergy));
      for (int k = 0; k < chi.length; k++) {
        report.append(format(" %7.2f", chi[k]));
      }
      report.append(format("\n"));
    }

    endTime = System.nanoTime();
    double took = (endTime - startTime) * NS_TO_MS;
    if (logTimings) {
      report.append(format("   Timing (ms): %.2f", took));
    }
    logger.info(report.toString());
  }

  /**
   * Follows Frenkel/Smit's derivation precisely, which for AMOEBA requires SCF calls in the inner
   * loops.
   */
  private void engageExpensive() {
    report.append(format(" Rosenbluth CBMC Move: %4d  %s\n", moveNumber, target));

    AminoAcid3 name = AminoAcid3.valueOf(target.getName());
    double chi[] = RotamerLibrary.measureRotamer(target, false);
    HashMap<Integer, BackBondedList> map = createBackBondedMap(name);

    // For each chi, create a test set from Boltzmann distr on torsion energy (Ubond).
    // Select from among this set based on Boltzmann weight of REMAINING energy (Uext).
    // The Uext partition function of each test set (wn) becomes a factor of the overall Rosenbluth
    // (Wn).
    double wn[] = new double[chi.length]; // factors of Wn
    double finalChi[] = new double[chi.length];
    for (int i = 0; i < chi.length; i++) {
      Torsion tors = map.get(i).torsion;
      //            TrialSet trialSet = boltzmannTorsionSet(tors, i, testSetSize, "bkn");
      TrialSet trialSet = expensiveTorsionSet(tors, i, testSetSize, "bkn");
      wn[i] = trialSet.sumExtBolt();
      if (i == 0) {
        Wn = wn[i];
      } else {
        Wn *= wn[i];
      }
      report.append(format(" wn,W running: %.2g %.2g", wn[i], Wn));
      logger.info(report.toString());

      // Choose a proposal move from amongst this trial set (bn).
      double rng = rand.nextDouble(wn[i]);
      double running = 0.0;
      for (int j = 0; j < trialSet.uExt.length; j++) {
        double uExtBolt = FastMath.exp(-beta * trialSet.uExt[j]);
        running += uExtBolt;
        if (rng < running) {
          finalChi[i] = trialSet.theta[j]; // Yes, I mean i then j.
          double prob = uExtBolt / wn[i] * 100;
          report.append(
              format(
                  "       Chose %d   %7.4f\t%7.4f\t  %4.1f%%\n",
                  j, trialSet.uExt[j], uExtBolt, prob));
          break;
        }
      }
    }
    report.append(format("    Wn Total:  %g\n", Wn));
    proposedMove = createRotamer(name, finalChi);

    // Reprise the above procedure for the old configuration.
    // The existing conformation forms first member of each test set (wo).
    // Overall Rosenbluth Wo is product of uExt partition functions.
    double[] wo = new double[chi.length]; // factors of Wo
    for (int i = 0; i < chi.length; i++) {
      Torsion tors = map.get(i).torsion;
      double ouDep = tors.energy(false); // original-conf uDep
      double ouExt = totalEnergy() - ouDep; // original-conf uExt
      double ouExtBolt = FastMath.exp(-beta * ouExt);
      TrialSet trialSet = expensiveTorsionSet(tors, i, testSetSize - 1, "bko");
      wo[i] = ouExtBolt + trialSet.sumExtBolt();
      if (i == 0) {
        Wo = wo[i];
      } else {
        Wo *= wo[i];
      }
    }
    report.append(format("    Wo Total:  %g\n", Wo));
    report.append(format("    Wn/Wo:     %g\n", Wn / Wo));

    target.revertState(origState);
    updateAll();
    if (verbose) {
      logger.info(report.toString());
    }
  }

  /**
   * For validation. Performs Monte Carlo chi moves WITHOUT biasing. Give ALL CHIs a random theta
   * simultaneously. Accept on the vanilla Metropolis criterion.
   */
  private void engageControlAll() {
    report.append(format(" Rosenbluth Control Move: %4d  %s\n", moveNumber, target));
    double origEnergy = totalEnergy();
    double[] origChi = RotamerLibrary.measureRotamer(target, false);
    double[] newChi = new double[origChi.length];
    arraycopy(origChi, 0, newChi, 0, origChi.length);
    for (int i = 0; i < origChi.length; i++) {
      if (doChi[i]) {
        double theta = rand.nextDouble(360.0) - 180;
        newChi[i] = theta;
      }
    }
    proposedChis = newChi;
    Rotamer newState = createRotamer(target, newChi);
    RotamerLibrary.applyRotamer(target, newState);

    finalEnergy = totalEnergy();
    if (this.finalEnergy < CATASTROPHE_THRESHOLD) {
      report.append("\nWARNING: Multipole catastrophe in CBMC.\n");
      report.append("  Discarding move.\n");
      target.revertState(origState);
      updateAll();
      Wn = -1.0;
      Wo = 1000;
      logger.info(report.toString());
    }

    double dU = finalEnergy - origEnergy;
    double criterion = FastMath.exp(-beta * dU);
    double rng = rand.nextDouble();
    report.append(format("    move (thetas):    "));
    for (int i = 0; i < newChi.length; i++) {
      report.append(format("%7.2f ", newChi[i]));
    }
    report.append(format("\n"));
    report.append(format("    orig, final, dU:  %.2g %.2g %.2g\n", origEnergy, finalEnergy, dU));
    report.append(format("    crit, rng:        %.2g %.2g\n", criterion, rng));
    if (rng < criterion) {
      accepted = true;
      numAccepted++;
      report.append(format(" Accepted! %5d    NewEnergy: %.4f    Chi:", numAccepted, finalEnergy));
      for (int k = 0; k < proposedChis.length; k++) {
        report.append(format(" %7.2f", proposedChis[k]));
      }
      report.append(format("\n"));
      updateAll();
      if (!noSnaps) {
        PDBFilter writer =
            new PDBFilter(molecularAssembly.getFile(), molecularAssembly, null, null);
        String filename = FilenameUtils.removeExtension(molecularAssembly.getFile().toString());
        filename = molecularAssembly.getFile().getAbsolutePath();
        if (!filename.contains("_mc")) {
          filename = FilenameUtils.removeExtension(filename) + "_mc.pdb";
        }
        File file = new File(filename);
        writer.writeFile(file, false);
      }
    } else {
      accepted = false;
      report.append(format(" Denied.   %5d    NewEnergy: %.4f    Chi:", numAccepted, origEnergy));
      for (int k = 0; k < origChi.length; k++) {
        report.append(format(" %7.2f", origChi[k]));
      }
      report.append(format("\n"));
      target.revertState(origState);
    }

    updateAll();
    endTime = System.nanoTime();
    double took = (endTime - startTime) * NS_TO_MS;
    if (logTimings) {
      report.append(format("   Timing (ms): %.2f", took));
    }
    logger.info(report.toString());
  }

  /**
   * wasAccepted.
   *
   * @return a boolean.
   */
  boolean wasAccepted() {
    return accepted;
  }

  private Rotamer createRotamer(AminoAcid3 name, double[] chi) {
    // Need to add sigma values to construct a new Rotamer with these chis.
    double[] values = new double[chi.length * 2];
    for (int k = 0; k < chi.length; k++) {
      int kk = 2 * k;
      values[kk] = chi[k];
      values[kk + 1] = 0.0;
    }
    return new Rotamer(name, values);
  }

  private Rotamer createRotamer(Residue res, double[] chi) {
    return createRotamer(AminoAcid3.valueOf(res.getName()), chi);
  }

  private void updateAll() {
    updateBonds();
    updateAngles();
    updateTorsions();
  }

  private void updateBonds() {
    for (Bond bond : target.getBondList()) {
      bond.update();
    }
  }

  private void updateAngles() {
    for (Angle angle : target.getAngleList()) {
      angle.update();
    }
  }

  private void updateTorsions() {
    for (Torsion torsion : target.getTorsionList()) {
      torsion.update();
    }
  }

  private double totalEnergy() {
    double x[] = new double[forceFieldEnergy.getNumberOfVariables() * 3];
    forceFieldEnergy.setEnergyTermState(Potential.STATE.BOTH);
    forceFieldEnergy.getCoordinates(x);
    return forceFieldEnergy.energy(false, verboseEnergies);
  }

  /**
   * Maps the back-bonded terms affected by key atoms in an amino acid. Here, 'key atom' refers to
   * each new rotamer-torsion-completing atom. e.g. VAL has 1 key atom (CG1), ARG has 4 key atoms
   * (CG,CD,NE,CZ). 'Back-bonded' means we only map terms that lead toward the backbone.
   */
  private HashMap<Integer, BackBondedList> createBackBondedMap(AminoAcid3 name) {
    HashMap<Integer, BackBondedList> map = new HashMap<>();
    List<Atom> chain = new ArrayList<>();
    Atom N = (Atom) target.getAtomNode("N");
    Atom CA = (Atom) target.getAtomNode("CA");
    Atom CB = (Atom) target.getAtomNode("CB");
    List<Atom> keyAtoms = new ArrayList<>();
    switch (name) {
      case VAL:
        {
          Atom CG1 = (Atom) target.getAtomNode("CG1");
          keyAtoms.add(CG1);
          keyAtoms.add(CB);
          break;
        }
      case LEU:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom CD1 = (Atom) target.getAtomNode("CD1");
          keyAtoms.add(CG);
          keyAtoms.add(CD1);
          break;
        }
      case ILE:
        {
          Atom CD1 = (Atom) target.getAtomNode("CD1");
          Atom CG1 = (Atom) target.getAtomNode("CG1");
          keyAtoms.add(CD1);
          keyAtoms.add(CG1);
          break;
        }
      case SER:
        {
          Atom OG = (Atom) target.getAtomNode("OG");
          Atom HG = (Atom) target.getAtomNode("HG");
          keyAtoms.add(OG);
          keyAtoms.add(HG);
          break;
        }
      case THR:
        {
          Atom OG1 = (Atom) target.getAtomNode("OG1");
          Atom HG1 = (Atom) target.getAtomNode("HG1");
          keyAtoms.add(OG1);
          keyAtoms.add(HG1);
          break;
        }
      case CYX:
      case CYD:
        {
          Atom SG = (Atom) target.getAtomNode("SG");
          keyAtoms.add(SG);
          break;
        }
      case PHE:
        {
          Atom CD1 = (Atom) target.getAtomNode("CD1");
          Atom CG = (Atom) target.getAtomNode("CG");
          keyAtoms.add(CG);
          break;
        }
      case PRO:
        {
          // Not allowed yet.
          Atom CD = (Atom) target.getAtomNode("CD");
          Atom CG = (Atom) target.getAtomNode("CG");
          keyAtoms.add(CG);
          keyAtoms.add(CD);
          break;
        }
      case TYR:
        {
          Atom CD1 = (Atom) target.getAtomNode("CD1");
          Atom CE2 = (Atom) target.getAtomNode("CE2");
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom CZ = (Atom) target.getAtomNode("CZ");
          Atom OH = (Atom) target.getAtomNode("OH");
          Atom HH = (Atom) target.getAtomNode("HH");
          // SPECIAL CASE: have to create map manualy.
          Bond b1 = CG.getBond(CB);
          Angle a1 = CG.getAngle(CB, CA);
          Torsion t1 = CG.getTorsion(CB, CA, N);
          Bond b2 = CD1.getBond(CG);
          Angle a2 = CD1.getAngle(CG, CB);
          Torsion t2 = CD1.getTorsion(CG, CB, CA);
          Bond b3 = HH.getBond(OH);
          Angle a3 = HH.getAngle(OH, CZ);
          Torsion t3 = HH.getTorsion(OH, CZ, CE2);
          BackBondedList bbl1 = new BackBondedList(b1, a1, t1);
          BackBondedList bbl2 = new BackBondedList(b2, a2, t2);
          BackBondedList bbl3 = new BackBondedList(b3, a3, t3);
          map.put(0, bbl1);
          map.put(1, bbl2);
          map.put(2, bbl3);
          return map; // Note the return here.
        }
      case TYD:
        {
          Atom CD1 = (Atom) target.getAtomNode("CD1");
          Atom CG = (Atom) target.getAtomNode("CG");
          keyAtoms.add(CG);
          keyAtoms.add(CD1);
          break;
        }
      case TRP:
        {
          Atom CD1 = (Atom) target.getAtomNode("CD1");
          Atom CG = (Atom) target.getAtomNode("CG");
          keyAtoms.add(CG);
          keyAtoms.add(CD1);
          break;
        }
      case HIS:
      case HID:
      case HIE:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom ND1 = (Atom) target.getAtomNode("ND1");
          keyAtoms.add(CG);
          keyAtoms.add(ND1);
          break;
        }
      case ASP:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          keyAtoms.add(CG);
          break;
        }
      case ASH:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom OD1 = (Atom) target.getAtomNode("OD1");
          keyAtoms.add(CG);
          keyAtoms.add(OD1);
          break;
        }
      case ASN:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom OD1 = (Atom) target.getAtomNode("OD1");
          keyAtoms.add(CG);
          keyAtoms.add(OD1);
          break;
        }
      case GLU:
      case GLH:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom CD = (Atom) target.getAtomNode("CD");
          Atom OE1 = (Atom) target.getAtomNode("OE1");
          keyAtoms.add(CG);
          keyAtoms.add(CD);
          keyAtoms.add(OE1);
          break;
        }
      case GLN:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom CD = (Atom) target.getAtomNode("CD");
          Atom OE1 = (Atom) target.getAtomNode("OE1");
          keyAtoms.add(CG);
          keyAtoms.add(CD);
          keyAtoms.add(OE1);
          break;
        }
      case MET:
        {
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom CE = (Atom) target.getAtomNode("CE");
          Atom SD = (Atom) target.getAtomNode("SD");
          keyAtoms.add(CG);
          keyAtoms.add(SD);
          keyAtoms.add(CE);
          break;
        }
      case LYS:
      case LYD:
        {
          Atom CD = (Atom) target.getAtomNode("CD");
          Atom CE = (Atom) target.getAtomNode("CE");
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom NZ = (Atom) target.getAtomNode("NZ");
          keyAtoms.add(CG);
          keyAtoms.add(CD);
          keyAtoms.add(CE);
          keyAtoms.add(NZ);
          break;
        }
      case ARG:
        {
          Atom CD = (Atom) target.getAtomNode("CD");
          Atom CG = (Atom) target.getAtomNode("CG");
          Atom CZ = (Atom) target.getAtomNode("CZ");
          Atom NE = (Atom) target.getAtomNode("NE");
          keyAtoms.add(CG);
          keyAtoms.add(CD);
          keyAtoms.add(NE);
          keyAtoms.add(CZ);
          break;
        }
      default:
        logger.severe(format("CBMC called on unsupported residue: %s", name.toString()));
    }
    // Build the chain and assign back-bonded terms.
    chain.add(N);
    chain.add(CA);
    chain.add(CB);
    chain.addAll(keyAtoms);
    for (int i = 3; i < chain.size(); i++) {
      Atom key = chain.get(i);
      Bond bond = key.getBond(chain.get(i - 1));
      Angle angle = key.getAngle(chain.get(i - 1), chain.get(i - 2));
      Torsion torsion = key.getTorsion(chain.get(i - 1), chain.get(i - 2), chain.get(i - 3));
      BackBondedList bbl = new BackBondedList(bond, angle, torsion);
      map.put(i - 3, bbl);
      //            report.append(String.format("    BackBondedMap: %s\t\t%s\n", key, torsion));
    }
    return map;
  }

  private void writeSnapshot(String suffix, boolean append) {
    if (snapshotWriter != null && !noSnaps) {
      snapshotWriter.write(suffix, append);
    }
  }

  /**
   * Collects data for plot of uTors vs. theta. This is for finding the appropriate offset to add to
   * each uTors such that the minimum is zero.
   *
   * @param allTors All torsions.
   */
  private void torsionSampler(List<Torsion> allTors) {
    updateAll();
    double[] origChi = measureLysine(target, false);
    StringBuilder sb = new StringBuilder();
    sb.append(format(" Torsion Sampling for %s\n", target.getName()));
    sb.append("   origChi:");
    for (double value : origChi) {
      sb.append(format(" %6.2f", value));
    }
    sb.append("\n");
    for (int k = 0; k < origChi.length; k++) {
      Torsion tors = allTors.get(k);
      AtomType type1 = tors.getAtomArray()[0].getAtomType();
      AtomType type2 = tors.getAtomArray()[1].getAtomType();
      AtomType type3 = tors.getAtomArray()[2].getAtomType();
      AtomType type4 = tors.getAtomArray()[3].getAtomType();
      sb.append(
          format(
              "   %d:    \"(%3d %3d %3s)  (%3d %3d %3s)  (%3d %3d %3s)  (%3d %3d %3s)\"\n",
              k,
              type1.type,
              type1.atomClass,
              type1.name,
              type2.type,
              type2.atomClass,
              type2.name,
              type3.type,
              type3.atomClass,
              type3.name,
              type4.type,
              type4.atomClass,
              type4.name));
    }
    logger.info(sb.toString());
    sb = new StringBuilder();
    sb.append(" Resi chi0 chi1 chi2 chi3 List<uTors> uTorsSum uTotal\n");
    logger.info(sb.toString());
    sb = new StringBuilder();
    boolean doChi0 = false, doChi1 = false;
    boolean doChi2 = true, doChi3 = true;
    double increment = 1.0;
    for (double chi0 = -180.0; chi0 < +180.0; chi0 += increment) {
      for (double chi1 = -180.0; chi1 <= +180.0; chi1 += increment) {
        for (double chi2 = -180.0; chi2 <= +180.0; chi2 += increment) {
          for (double chi3 = -180.0; chi3 <= +180.0; chi3 += increment) {
            sb.append(format(" %3s", target.getName()));
            double[] newChi = new double[4];
            if (doChi0) {
              newChi[0] = chi0;
            } else {
              newChi[0] = origChi[0];
            }
            if (doChi1) {
              newChi[1] = chi1;
            } else {
              newChi[1] = origChi[1];
            }
            if (doChi2) {
              newChi[2] = chi2;
            } else {
              newChi[2] = origChi[2];
            }
            if (doChi3) {
              newChi[3] = chi3;
            } else {
              newChi[3] = origChi[3];
            }
            Rotamer newState = createRotamer(target, newChi);
            RotamerLibrary.applyRotamer(target, newState);
            for (int wut = 0; wut < origChi.length; wut++) {
              sb.append(format(" %3.0f", newChi[wut]));
            }
            writeSnapshot(format("%.0f", chi1), false);
            double uTorsSum = 0.0;
            for (int k = 0; k < allTors.size(); k++) {
              double uTors = allTors.get(k).energy(false);
              sb.append(format(" %5.2f", uTors));
              uTorsSum += uTors;
            }
            sb.append(format(" %5.2f", uTorsSum));
            double totalE = 1000.0;
            try {
              totalE = totalEnergy();
            } catch (Exception ex) {
              //
            }
            if (totalE > 1000.0) {
              totalE = 1000.0;
            }
            sb.append(format(" %10.4f", totalE));
            sb.append(format("\n"));
            logger.info(sb.toString());
            sb = new StringBuilder();
            if (!doChi3) {
              break;
            }
          }
          if (!doChi2) {
            break;
          }
        }
        if (!doChi1) {
          break;
        }
      }
      if (!doChi0) {
        break;
      }
    }

    if (torsionSampling) {
      try {
        File output = new File("torsionSampler.log");
        BufferedWriter bw = new BufferedWriter(new FileWriter(output));
        bw.write(sb.toString());
        bw.close();
      } catch (IOException ex) {
        //
      }
    }
    System.exit(0);
  }

  /**
   * measureLysine.
   *
   * @param residue a {@link ffx.potential.bonded.Residue} object.
   * @param print a boolean.
   * @return an array of {@link double} objects.
   */
  private double[] measureLysine(Residue residue, boolean print) {
    if (!residue.getName().contains("LY")
        || (residue.getAminoAcid3() != AminoAcid3.LYS
            && residue.getAminoAcid3() != AminoAcid3.LYD)) {
      logger.severe("Yeah that ain't a lysine.");
    }
    double[] chi = new double[4];
    List<Torsion> torsions = residue.getTorsionList();
    Atom N = (Atom) residue.getAtomNode("N");
    Atom CA = (Atom) residue.getAtomNode("CA");
    Atom CB = (Atom) residue.getAtomNode("CB");
    Atom CD = (Atom) residue.getAtomNode("CD");
    Atom CE = (Atom) residue.getAtomNode("CE");
    Atom CG = (Atom) residue.getAtomNode("CG");
    Atom NZ = (Atom) residue.getAtomNode("NZ");
    logger.info(
        format(
            " Here's the atoms I found: \n%s\n%s\n%s\n%s\n%s\n%s\n%s", N, CA, CB, CD, CE, CG, NZ));
    logger.info(format(" Num torsions: %d", torsions.size()));
    int count = 0;
    for (Torsion torsion : torsions) {
      torsion.energy(false);
      logger.info(format(" Torsion numba %d: %s", count++, torsion));
      if (torsion.compare(N, CA, CB, CG)) {
        chi[0] = torsion.getValue();
        if (print) {
          logger.info(torsion.toString());
        }
      }
      if (torsion.compare(CA, CB, CG, CD)) {
        chi[1] = torsion.getValue();
        if (print) {
          logger.info(torsion.toString());
        }
      }
      if (torsion.compare(CB, CG, CD, CE)) {
        chi[2] = torsion.getValue();
        if (print) {
          logger.info(torsion.toString());
        }
      }
      if (torsion.compare(CG, CD, CE, NZ)) {
        chi[3] = torsion.getValue();
        if (print) {
          logger.info(torsion.toString());
        }
      }
    }
    return chi;
  }

  /** Mode of the RosenbluthChiAllMove instance. */
  public enum MODE {
    EXPENSIVE,
    CHEAP,
    CTRL_ALL;
  }

  /** Provides lookup values that make true the inequality: uTorsion + offset .ge. 0.0 */
  private enum TORSION_OFFSET_AMPRO13 {
    LYS0(1.610000),
    LYD0(1.610000),
    LYS1(0.939033),
    LYD1(0.939033),
    LYS2(1.000000),
    LYD2(1.000000),
    LYS3(0.800000),
    LYD3(0.800000);

    public final double offset;

    TORSION_OFFSET_AMPRO13(double offset) {
      this.offset = offset;
    }
  }

  private class TrialSet {

    public final Rotamer rotamer[];
    public final double uDep[];
    public final double uExt[];
    public final double theta[];

    public TrialSet(int setSize) {
      rotamer = new Rotamer[setSize];
      uDep = new double[setSize];
      uExt = new double[setSize];
      theta = new double[setSize];
    }

    public double prodExtBolt() {
      double prod = 0.0;
      for (int i = 0; i < uExt.length; i++) {
        prod *= FastMath.exp(-beta * uExt[i]);
      }
      return prod;
    }

    // We need to SUM over "members of the test set" i.e. ONLY do this for different TRIALS (b1 ...
    // bn)
    // AND THEN We want to MULTIPLY over "segments in the polymer chain" i.e. ONLY do prod over
    // different CHIS
    public double sumExtBolt() {
      double sum = 0.0;
      for (int i = 0; i < uExt.length; i++) {
        sum += FastMath.exp(-beta * uExt[i]);
      }
      return sum;
    }
  }

  private class BackBondedList {

    public final Bond bond;
    public final Angle angle;
    public final Torsion torsion;

    public BackBondedList(Bond bond, Angle angle, Torsion tors) {
      this.bond = bond;
      this.angle = angle;
      this.torsion = tors;
    }
  }

  private class SnapshotWriter {

    private final MolecularAssembly mola;
    private final PDBFilter filter;
    private final boolean interleaving;

    private SnapshotWriter(MolecularAssembly mola, boolean interleaving) {
      this.mola = mola;
      this.interleaving = interleaving;
      this.filter = new PDBFilter(mola.getFile(), mola, null, null);
      this.filter.setLogWrites(false);
    }

    private void write(String suffix, boolean append) {
      String filename =
          FilenameUtils.removeExtension(mola.getFile().toString())
              + "."
              + suffix
              + "-"
              + moveNumber;
      if (interleaving) {
        filename = mola.getFile().getAbsolutePath();
        if (!filename.contains("dyn")) {
          filename = FilenameUtils.removeExtension(filename) + "_dyn.pdb";
        }
      }
      File file = new File(filename);
      filter.setLogWrites(false);
      filter.writeFile(file, append);
    }
  }
}
