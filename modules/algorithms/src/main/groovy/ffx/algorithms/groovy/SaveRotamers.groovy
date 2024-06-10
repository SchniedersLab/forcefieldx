//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2024.
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
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.optimize.TitrationManyBody
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.*
import ffx.potential.bonded.RotamerLibrary.NucleicSugarPucker
import ffx.potential.parameters.TitrationUtils
import ffx.potential.parsers.PDBFilter
import org.apache.commons.configuration2.CompositeConfiguration
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.getExtension
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * The SaveRotamers script saves out rotamers.
 * <br>
 * Usage:
 * <br>
 * ffxc SaveRotamers [options] &lt;filename&gt;
 */
@Command(description = " Save out rotamers.", name = "SaveRotamers")
class SaveRotamers extends AlgorithmsScript {

  /**
   * -c or --chain to specify chain
   */
  @Option(names = ['--chain', '-c'], paramLabel = "\" \"",
      description = 'Single character chain name.')
  char c = ' '

  /**
   * -l or --library to select rotamer library
   */
  @Option(names = ['--library', '-l'], paramLabel = "1", defaultValue = "1",
      description = 'Available rotamer libraries are (1) Ponder and Richards or (2) Richardson.')
  int library

  /**
   * -r or --resid to select residue number
   */
  @Option(names = ['--resid', '-r'], paramLabel = "1", defaultValue = "1",
      description = 'Residue number.')
  int resID

  /**
   * -i or --independent to draw nucleic acid rotamers independently of chain context
   */
  @Option(names = ['--independent', '-i'], paramLabel = 'false', defaultValue = 'false',
      description = 'Independent draws nucleic acid rotamers independently of chain context.')
  boolean independent

  /**
   * -s or --start to select first rotamer to draw. Indexed form rotamer 0
   */
  @Option(names = ['--start', '-s'], paramLabel = "0", defaultValue = "0",
      description = 'First rotamer to draw (indexed from rotamer 0).')
  int start

  /**
   * -f or --finish to select last rotamer to draw. Indexed from rotamer 0
   */
  @Option(names = ['--finish', '-f'], paramLabel = "-1", defaultValue = "-1",
      description = 'Last rotamer to draw (indexed from rotamer 0).')
  int finish

  /**
   * -x or --all to draw all rotamers beginning from the passed rotamer number (overrides other options). Indexed from rotamer 0.
   */
  @Option(names = ['--all', '-x'], paramLabel = "-1", defaultValue = "-1",
      description = 'Draw all rotamers beginning from the passed index (overrides other options).')
  int all

  /**
   * -u or --upstreamPucker to adjust the pucker of the 5\' residue to match the rotamer. Use flag to turn on
   */
  @Option(names = ['--upstreamPucker', '-u'], paramLabel = 'false', defaultValue = 'false',
      description = 'Adjusts the pucker of the 5\' residue to match the rotamer.')
  boolean upstreamPucker

  /**
   * --tR or --titrateResidue to allow for titration states to be saved
   */
  @Option(names = ['--tR', '--titrateResidue'], paramLabel = 'false', defaultValue = 'false',
      description = 'Titrate residues.')
  boolean titrateResidue

  /**
   * The final argument is an XYZ or PDB coordinate file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'The atomic coordinate file in XYZ or PDB format.')
  private String filename = null

  TitrationManyBody titrationManyBody
  TitrationUtils titrationUtils

  /**
   * SaveRotamers Constructor.
   */
  SaveRotamers() {
    this(new Binding())
  }

  /**
   * SaveRotamers Constructor.
   * @param binding Groovy Binding to use.
   */
  SaveRotamers(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  SaveRotamers run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }


    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()
    CompositeConfiguration properties = activeAssembly.getProperties()

    if (titrateResidue) {
      List<Residue> residues = activeAssembly.getResidueList()
      List<Integer> resNumberList = new ArrayList<>()
      for (Residue residue : residues) {
        resNumberList.add(residue.getResidueNumber())
      }
      titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
          resNumberList, 7.0)
      MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly()
      setActiveAssembly(protonatedAssembly)
    }
    RotamerLibrary rLib = new RotamerLibrary(
        RotamerLibrary.ProteinLibrary.intToProteinLibrary(library), true)
    String chain = Character.toString(c)

    boolean saveAllRotamers = false
    int allStart = 0
    // Start from to all
    if (all > -1) {
      allStart = all
      saveAllRotamers = true
    }

    logger.info("\n Saving rotamers for residue number " + resID + " of chain " + chain + ".")

    RotamerLibrary.initializeDefaultAtomicCoordinates(activeAssembly.getChains())
    Polymer polymer = activeAssembly.getChain(chain)
    if (polymer == null) {
      logger.info(" Polymer + " + chain + " does not exist.")
      return this
    }
    Residue residue = polymer.getResidue(resID)
    if (residue == null) {
      logger.info(" Residue + " + resID + " does not exist.")
      return this
    }

    residue.setRotamers(rLib)
    Rotamer[] rotamers = residue.getRotamers()
    if (rotamers == null) {
      logger.severe(" There are no rotamers for residue + " + residue.toString())
    }
    int nrotamers = rotamers.length

    boolean isDeoxy = false // Applies to prevResidue.
    Residue prevResidue = null
    // If upstreamPucker is specified true, ensure that it is valid to apply it
    // (residue is nucleic acid with an upstream partner), and set isDeoxy.
    // If it is invalid, set upstreamPucker false.
    if (upstreamPucker) {
      // Exception gets thrown if it's an amino acid, since "NA" is undefined.
      try {
        if (residue.getResidueType() == Residue.ResidueType.NA) {
          prevResidue = (Residue) residue.getPreviousResidue()
          // If no previous residue, set upstream pucker false.
          // The method used will ensure prevResidue is a nucleic acid.
          if (prevResidue == null) {
            upstreamPucker = false
          } else {
            Atom HOs = (Atom) prevResidue.getAtomNode("HO\'")
            if (HOs == null) {
              isDeoxy = true
            }
          }
        } else {
          upstreamPucker = false
        }
      } catch (Exception e) {
        upstreamPucker = false
      }
    }

    String ext = getExtension(filename)
    filename = removeExtension(filename)
    Set<Atom> excludeAtoms = new HashSet<>()
    List<Residue> removeAtomList = new ArrayList<>()
    if (saveAllRotamers) {
      if (allStart >= nrotamers) {
        logger.info(" Specified start range is outside of rotamer range. No action taken.")
      } else {
        for (int i = allStart; i < nrotamers; i++) {
          RotamerLibrary.applyRotamer(residue, rotamers[i], independent)
          logger.info(rotamers[i].getName())
          if (upstreamPucker) {
            double prevDelta = rotamers[i].chi1
            NucleicSugarPucker sugarPucker = NucleicSugarPucker.checkPucker(prevDelta, isDeoxy)
            RotamerLibrary.applySugarPucker(prevResidue, sugarPucker, isDeoxy, true)
          }
          if (ext.toUpperCase().contains("XYZ")) {
            logger.info(format(" Saving rotamer %d", i))
            potentialFunctions.saveAsXYZ(activeAssembly, new File(filename + ".xyz"))
          } else {
            if (titrateResidue) {
              excludeAtoms.clear()
              removeAtomList.clear()
              removeAtomList.add(residue)

              int[] currentRot = [i]
              boolean isTitrating = titrationManyBody.excludeExcessAtoms(excludeAtoms, currentRot,
                  removeAtomList)
              //properties.setProperty("standardizeAtomNames", "false")
              File modelFile = new File(filename + ".pdb")
              PDBFilter pdbFilter = new PDBFilter(modelFile, activeAssembly, activeAssembly.getForceField(),
                  properties)
              if (!pdbFilter.writeFile(modelFile, false, excludeAtoms, true, true)) {
                logger.info(format(" Save failed for %s", activeAssembly))
              }
            } else {
              logger.info(format(" Saving rotamer %d", i))
              potentialFunctions.saveAsPDB(activeAssembly, new File(filename + ".pdb"))
            }
          }
        }
      }
    } else {
      if (start >= nrotamers) {
        logger.info(" Specified start range is outside of rotamer range. No action taken.")
      } else {
        if (finish >= nrotamers) {
          finish = nrotamers - 1
        } else if (finish < start) {
          logger.info(
              " Specified finish point is before the start point drawing only rotamer " + start)
          finish = start
        }
        for (int i = start; i <= finish; i++) {
          RotamerLibrary.applyRotamer(residue, rotamers[i], independent)
          if (upstreamPucker) {
            double prevDelta = rotamers[i].chi1
            NucleicSugarPucker sugarPucker = NucleicSugarPucker.checkPucker(prevDelta, isDeoxy)
            RotamerLibrary.applySugarPucker(prevResidue, sugarPucker, isDeoxy, true)
          }
          if (ext.toUpperCase().contains("XYZ")) {
            logger.info(format(" Saving rotamer %d", i))
            potentialFunctions.saveAsXYZ(activeAssembly, new File(filename + ".xyz"))
          } else {
            logger.info(format(" Saving rotamer %d", i))
            potentialFunctions.saveAsPDB(activeAssembly, new File(filename + ".pdb"))
          }
        }
      }
    }

    return this
  }
}
