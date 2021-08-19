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
package ffx.potential.groovy

import ffx.crystal.Crystal
import ffx.crystal.SpaceGroup
import ffx.crystal.SymOp
import ffx.numerics.Potential
import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static ffx.crystal.ReplicatesCrystal.replicatesCrystalFactory
import static ffx.crystal.SpaceGroupDefinitions.spaceGroupFactory
import static ffx.crystal.SpaceGroupInfo.getCCDCPercent
import static ffx.crystal.SpaceGroupInfo.getPDBRank
import static ffx.crystal.SymOp.randomSymOpFactory
import static java.lang.String.format

/**
 * The PrepareSpaceGroups creates sub-directories for the selected space groups.
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsXYZ [options] &lt;filename&gt;
 */
@Command(description = " Create sub-directories for selected space groups.", name = "ffxc PrepareSpaceGroups")
class PrepareSpaceGroups extends PotentialScript {

  /**
   * -r or --PDB only consider space groups ranked above the supplied cut-off in the Protein Databank.
   */
  @Option(names = ['-r', '--PDB'], paramLabel = '231',
      description = 'Only consider space groups populated in the PDB above the specified rank (231 excludes none).')
  int rank = 231

  /**
   * -p or --CSD only consider space groups populated above the specified percentage in the CSD.
   */
  @Option(names = ['-p', '--CSD'], paramLabel = '1.0',
      description = 'Only consider space groups populated above the specified percentage in the CSD.')
  double percent = 0.0

  /**
   * -c or --chiral to create directories only for chiral space groups.
   */
  @Option(names = ['-c', '--chiral'], paramLabel = 'false',
      description = 'Only consider chiral space groups.')
  boolean chiral = false

  /**
   *  -a or --achiral to create directories only for achiral space groups.
   */
  @Option(names = ['-a', '--achiral'], paramLabel = 'false',
      description = 'Only consider achiral space groups.')
  boolean achiral = false

  /**
   * -sym or --symOp random Cartesian symmetry operator will use the specified translation range -Arg .. Arg (default of 1.0 A).
   */
  @Option(names = ['--rsym', '--randomSymOp'], paramLabel = 'Arg',
      description = 'Random Cartesian sym op will choose a translation in the range -Arg .. Arg (default of 1.0 A).')
  double symScalar = 1.0

  /**
   * -d or --density random unit cell axes will be used achieve the specified density (g/cc).
   */
  @Option(names = ['-d', '--density'], paramLabel = '1.0',
      description = 'Random unit cell axes will be used to achieve the specified density (default of 1.0 g/cc).')
  double density = 1.0

  /**
   * -sg or --spacegroup prepare a directory for a single spacegroup (no default).
   */
  @Option(names = ['--sg', '--spacegroup'],
      description = 'Prepare a directory for a single spacegroup.')
  String sg

  /**
   * The final argument is a single filename in PDB or XYZ format.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in PDB or XYZ format.')
  String filename = null

  public int numberCreated = 0
  private ForceFieldEnergy energy

  /**
   * PrepareSpaceGroups Constructor.
   */
  PrepareSpaceGroups() {
    this(new Binding())
  }

  /**
   * PrepareSpaceGroups Constructor.
   * @param binding Groovy Binding to use.
   */
  PrepareSpaceGroups(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  PrepareSpaceGroups run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Turn off electrostatic interactions.
    System.setProperty("mpoleterm", "false")

    // Make sure an a-axis value is set.
    System.setProperty("a-axis", "10.0")

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    logger.info("\n Preparing space group directories for " + filename)

    energy = activeAssembly.getPotentialEnergy()
    CompositeConfiguration properties = activeAssembly.getProperties()

    File coordFile = activeAssembly.getFile()
    String coordName = FilenameUtils.getName(coordFile.getPath())

    File propertyFile = new File(properties.getString("propertyFile"))

    Atom[] atoms = activeAssembly.getAtomArray()
    double mass = activeAssembly.getMass()
    if (density < 0.0) {
      density = 1.0
    }

    for (int num = 1; num <= 230; num++) {
      SpaceGroup spacegroup = spaceGroupFactory(num)
      if (sg) {
        SpaceGroup spacegroup2 = spaceGroupFactory(sg)
        if (spacegroup2 == null) {
          logger.info(format("\n Space group %s was not recognized.\n", sg))
          return this
        }
        if (spacegroup2.number != spacegroup.number) {
          continue
        }
      }

      if (chiral && !spacegroup.respectsChirality()) {
        continue
      }

      if (achiral && spacegroup.respectsChirality()) {
        continue
      }

      if (getCCDCPercent(num) < percent) {
        continue
      }

      if (getPDBRank(spacegroup) > rank) {
        continue
      }

      logger.info(format("\n Preparing %s (CSD percent: %7.4f, PDB Rank: %d)",
          spacegroup.shortName, getCCDCPercent(num), getPDBRank(spacegroup)))

      // Create the directory.
      String sgDirName = spacegroup.shortName.replace('/', '_')

      File baseDir = baseDir
      File sgDir
      if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
        sgDir = new File(FilenameUtils.getFullPath(coordFile.getAbsolutePath()) + sgDirName)
      } else {
        sgDir = new File(baseDir.getAbsolutePath() + sgDirName)
      }

      if (!sgDir.exists()) {
        logger.info("\n Creating space group directory: " + sgDir.toString())
        sgDir.mkdir()
      }

      double[] abc = spacegroup.latticeSystem.resetUnitCellParams()
      Crystal crystal = new Crystal(abc[0], abc[1], abc[2], abc[3], abc[4], abc[5],
          spacegroup.shortName)
      crystal.setDensity(density, mass)
      double cutoff2 = energy.getCutoffPlusBuffer() * 2.0
      crystal = replicatesCrystalFactory(crystal, cutoff2)

      // Turn off special position checks.
      crystal.setSpecialPositionCutoff(0.0)
      crystal.getUnitCell().setSpecialPositionCutoff(0.0)
      energy.setCrystal(crystal)

      if (symScalar > 0.0) {
        SymOp symOp = randomSymOpFactory(symScalar)
        logger.info(format("\n Applying random Cartesian SymOp:\n %s", symOp.toString()))
        double[] xyz = new double[3]
        for (int i = 0; i < atoms.length; i++) {
          atoms[i].getXYZ(xyz)
          crystal.applyCartesianSymOp(xyz, xyz, symOp)
          atoms[i].setXYZ(xyz)
        }
      }

      // Save the coordinate file.
      File sgFile = new File(sgDir.getAbsolutePath() + File.separator + coordName)
      logger.info(" Saving " + sgDirName + " coordinates to: " + sgFile.toString())
      potentialFunctions.save(activeAssembly, sgFile)

      File keyFile = new File(
          FilenameUtils.removeExtension(sgFile.getAbsolutePath()) + ".properties")
      logger.info(" Saving " + sgDirName + " properties to: " + keyFile.toString())

      numberCreated++

      // Write the new properties file.
      new BufferedReader(new FileReader(propertyFile)).withCloseable {keyReader ->
        new PrintWriter(new BufferedWriter(new FileWriter(keyFile))).withCloseable {keyWriter ->
          try {
            while (keyReader.ready()) {
              String line = keyReader.readLine().trim()
              if (line != null && !line.equalsIgnoreCase("")) {
                String[] tokens = line.split(" +")
                if (tokens != null && tokens.length > 1) {
                  String first = tokens[0].toLowerCase()
                  if (first == "a-axis" || first == "b-axis" || first == "c-axis" ||
                      first == "alpha" || first == "beta" || first == "gamma" ||
                      first == "spacegroup") {
                    // This line is skipped, and updated below.
                    logger.fine(format(" Updating : %s", line))
                  } else if (first == "parameters") {
                    if (tokens.length > 1) {
                      keyWriter.println(format("parameters ../%s", tokens[1]))
                    } else {
                      logger.warning("Parameter file may not have been specified")
                    }
                  } else {
                    keyWriter.println(line)
                  }
                } else {
                  keyWriter.println(line)
                }
              }
            }
            // Update the space group and unit cell parameters.
            keyWriter.println(format("spacegroup %s", spacegroup.shortName))
            keyWriter.println(format("a-axis %12.8f", crystal.a))
            keyWriter.println(format("b-axis %12.8f", crystal.b))
            keyWriter.println(format("c-axis %12.8f", crystal.c))
            keyWriter.println(format("alpha  %12.8f", crystal.alpha))
            keyWriter.println(format("beta   %12.8f", crystal.beta))
            keyWriter.println(format("gamma  %12.8f", crystal.gamma))
          } catch (IOException ex) {
            logger.warning(" Exception writing keyfile." + ex.toString())
          }
        }
      }
    }

    System.clearProperty("mpoleterm")
    System.clearProperty("a-axis")

    return this
  }

  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (energy == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList((Potential) energy)
    }
    return potentials
  }

}
