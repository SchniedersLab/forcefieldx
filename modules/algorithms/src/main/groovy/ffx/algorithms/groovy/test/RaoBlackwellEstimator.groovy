//******************************************************************************
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
//******************************************************************************
package ffx.algorithms.groovy.test

import edu.rit.mp.DoubleBuf
import edu.rit.pj.Comm
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.dynamics.MolecularDynamics
import ffx.algorithms.dynamics.MDEngine
import ffx.algorithms.dynamics.MolecularDynamicsOpenMM
import ffx.numerics.Potential
import ffx.potential.bonded.Residue
import ffx.potential.cli.WriteoutOptions
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import java.util.logging.Level
import java.util.logging.LogRecord

import static java.lang.String.format

/**
 * Use the BAR algorithm to estimate a free energy difference for a CpHMD system.
 * <br>
 * Usage: Umod calculation for model compounds
 * <br>
 * ffxc test.RaoBlackwellEstimator [options] &lt;filename&gt [file2...];
 */
@Command(description = " Use the Rao-Blackwell estimator to get a free energy difference for residues in a CpHMD system.", name = "test.RaoBlackwellEstimator")
class RaoBlackwellEstimator extends AlgorithmsScript {

  @Option(names = ['--aFi', '--arcFile'], paramLabel = "traj",
          description = 'A file containing the the PDB from which to build the ExtendedSystem. There is currently no default.')
  private String arcFileName = null

  @Option(names = ['--upState'], paramLabel = "1.0",
          description = 'State to perturb up to.')
  private double upState = 1.0

  @Option(names = ['--downState'], paramLabel = "0.0",
          description = 'State to perturb down to.')
  private double downState = 0.0


  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files",
          description = "PDB input file in the same directory as the ARC file.")
  private String filename

  private Potential forceFieldEnergy
  ArrayList<Double>[] oneZeroDeltaLists

  /**
   * Thermodynamics Constructor.
   */
  RaoBlackwellEstimator() {
    this(new Binding())
  }

  /**
   * Thermodynamics Constructor.
   * @param binding The Groovy Binding to use.
   */
  RaoBlackwellEstimator(Binding binding) {
    super(binding)
  }

  RaoBlackwellEstimator run() {

    if (!init()) {
      return this
    }

    if(arcFileName == null) {
      logger.severe("Archive file needs to be specified with --arcFile or --aFi.")
    }

    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }
    forceFieldEnergy = activeAssembly.getPotentialEnergy()

    // Set the filename.
    String filename = activeAssembly.getFile().getAbsolutePath()

    // Initialize and attach extended system first.
    ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, 7.0, null)
    int numESVs = esvSystem.extendedResidueList.size()

    // Set up the XPHFilter.
    File arcFile = new File(arcFileName)
    XPHFilter xphFilter = new XPHFilter(
            arcFile,
            activeAssembly,
            activeAssembly.getForceField(),
            activeAssembly.getProperties(),
            esvSystem)
    xphFilter.readFile()
    logger.info("Reading ESV lambdas from XPH file")

    esvSystem.setFixedTitrationState(true)
    esvSystem.setFixedTautomerState(true)
    forceFieldEnergy.attachExtendedSystem(esvSystem)
    logger.info(format(" Attached extended system with %d residues.", numESVs))

    double[] x = new double[forceFieldEnergy.getNumberOfVariables()]
    forceFieldEnergy.getCoordinates(x)
    forceFieldEnergy.energy(x, true)

    // Get pH from ARC file.
    double pH = 0.0
    // Read the first line of pHFind.
    try
    {
      BufferedReader br = new BufferedReader(new FileReader(arcFile))
      String line = br.readLine()
      String[] parts = line.split(" ")
      for(int i = 0; i < parts.length; i++) {
        if (parts[i].equals("pH")) {
          pH = Double.parseDouble(parts[i+1])
        }
      }
      br.close()
    } catch (IOException e) {
      e.printStackTrace()
    }
    logger.info("Setting ESV pH to " + pH)
    esvSystem.setConstantPh(pH)

    int index = 0
    while(xphFilter.readNext()) {
      logger.info("Reading frame " + index + 1)
      for (int i = 0; i < numESVs; i++) {
        Residue res = esvSystem.extendedResidueList.get(i)
        double titrationState = esvSystem.getTitrationLambda(res)
        esvSystem.setTitrationLambda(res, 0, false)
        forceFieldEnergy.getCoordinates(x)
        double zeroEnergy = forceFieldEnergy.energy(x, false)

        esvSystem.setTitrationLambda(res, 1, false)
        forceFieldEnergy.getCoordinates(x)
        double oneEnergy = forceFieldEnergy.energy(x, false)

        esvSystem.setTitrationLambda(res, titrationState, false)
        forceFieldEnergy.getCoordinates(x)
        double selfEnergy = forceFieldEnergy.energy(x, false)

        oneZeroDeltaLists[i].add(oneEnergy - zeroEnergy)
      }
      index++
    }

    // Calculate the Rao-Blackwell estimator for each residue.

    for(int i = 0; i < numESVs; i++) {
      // Calculate the free energy differences.
      ArrayList<Double> deltaU = selfLists[i]
      double temperature = 298.0
      double boltzmann = 0.001985875
      double beta = 1.0 / (temperature * boltzmann)

      ArrayList<Double> deltaExp = exp(mult(-beta,deltaU))
      ArrayList<Double> numerator = mult(beta,mult(deltaU,deltaExp)) / subtract(1.0,deltaExp)
      ArrayList<Double> denominator = mult(beta,deltaU) / subtract(1.0,deltaExp)
      double deltaG = -(1.0 / beta) * Math.log(average(numerator) / average(denominator))
    }
    return this
  }

  private static double average(ArrayList<Double> list) {
    double sum = 0.0
    for (Double d : list) {
      sum += d
    }
    return sum / list.size()
  }

    private static ArrayList<Double> mult(double a, ArrayList<Double> u) {
        ArrayList<Double> result = new ArrayList<Double>()
        for (Double d : u) {
            result.add(a * d)
        }
      return result
    }

  private static ArrayList<Double> mult(ArrayList<Double> v, ArrayList<Double> u) {
    if (v.size() != u.size()) {
      throw new IllegalArgumentException("Vector sizes must be equal.")
    }
    ArrayList<Double> result = new ArrayList<Double>()
    for (int i = 0; i < v.size(); i++) {
      result.add(v.get(i) * u.get(i))
    }
    return result
  }

  private static ArrayList<Double> subtract(double a, ArrayList<Double> u) {
    ArrayList<Double> result = new ArrayList<Double>()
    for (Double d : u) {
      result.add(a - d)
    }
    return result
  }

  private static ArrayList<Double> exp(ArrayList<Double> u) {
    ArrayList<Double> result = new ArrayList<Double>()
    for (Double d : u) {
      result.add(Math.exp(d))
    }
    return result
  }



}