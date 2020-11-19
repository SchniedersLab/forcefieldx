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
package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.algorithms.optimize.Minimize
import ffx.algorithms.optimize.Minimize.MinimizationEngine
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.*
import ffx.potential.bonded.RotamerLibrary.ProteinLibrary
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

/**
 * The CreateRotamers script creates a set of conformation dependent rotamers.
 * <br>
 * Usage:
 * <br>
 * ffxc CreateRotamers [options] &lt;filename&gt;
 */
@Command(description = " Creates a set of conformation dependent rotamers.", name = "ffxc CreateRotamers")
class CreateRotamers extends AlgorithmsScript {

  @Mixin
  MinimizeOptions minimizeOptions

  // TODO: instead @Mixin a subset of current ManyBodyOptions.

  /**
   * -L or --library Choose either Ponder and Richards (1) or Richardson (2)
   * rotamer library.
   */
  @Option(names = ["-L", "--library"], paramLabel = "2",
      description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
  String library = 2

  /**
   * -d or --rmsd Set the RMSD cut off for rotamer exclusion
   * Any rotamer with an RMSD from any previous rotamer that is
   * less than or equal to the cut off will be thrown out.
   */
  @Option(names = ["-d", "--rmsd"], paramLabel = "0.1",
      description = "RMSD cut off for rotamer exclusion: only rotamers with an RMSD greater than this cut off value to all previously saved rotamers will be kept")
  double rmsdCutoff = 0.1

  /**
   * The final argument should be a filename.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = 'Atomic coordinate file in PDB or XYZ format.')
  List<String> filenames = null

  public ForceFieldEnergy forceFieldEnergy = null

  /**
   * CreateRotamers Constructor.
   */
  CreateRotamers() {
    this(new Binding())
  }

  /**
   * CreateRotamers Constructor.
   * @param binding The Groovy Binding to use.
   */
  CreateRotamers(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  CreateRotamers run() {

    if (!init()) {
      return this
    }

    if (filenames != null && filenames.size() > 0) {
      MolecularAssembly[] assemblies = [algorithmFunctions.open(filenames.get(0))]
      activeAssembly = assemblies[0]
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    String filename = activeAssembly.getFile().getAbsolutePath()
    logger.info(" Running CreateRotamers on " + filename)

    Atom[] atoms = activeAssembly.getAtomArray()
    int nAtoms = atoms.length

    // Set all atoms to be "inactive".
    for (int i = 0; i < nAtoms; i++) {
      atoms[i].setActive(false)
    }

    // For now, always use the original coordinates as a (fixed) rotamer.
    boolean useOriginalRotamers = true

    // AA Library
    RotamerLibrary rotamerLibrary = new RotamerLibrary(ProteinLibrary.getProteinLibrary(library),
        useOriginalRotamers)

    // Initialize Default NA Coordinates
    Polymer[] polymers = activeAssembly.getChains()
    RotamerLibrary.initializeDefaultAtomicCoordinates(polymers)

    // Get the residue list.
    List<Residue> residues = activeAssembly.getResidueList().stream().
        filter({
          Residue r ->
            Rotamer[] rots = r.getRotamers(rotamerLibrary)
            return rots != null && rots.length > 1
        }).collect(Collectors.toList())

    logger.info(String.format(" Number of residues: %d\n", residues.size()))

    // Loop over Residues and set side chain atoms to not be used.
    for (Residue residue : residues) {
      for (Atom atom : residue.getVariableAtoms()) {
        atom.setUse(false)
      }
    }

    // Create .rot file name: should match input file name and end in ".rot"
    String rotFileName = String.format("%s.rot", FilenameUtils.removeExtension(filename))

    BufferedWriter bw = null
    try {
      bw = new BufferedWriter(new FileWriter(new File(rotFileName)))

      // TODO: Make this ALGORITHM:[ALGORITHM]:[box/window number] instead of assuming global:1.
      bw.write("ALGORITHM:GLOBAL:1")
      bw.newLine()

      // Loop over Residues
      for (Residue residue : residues) {

        StringBuilder resLine = new StringBuilder(" RES:")
        resLine.append(residue.getChainID()).append(":")
        resLine.append(residue.getSegID()).append(":")
        resLine.append(residue.getName()).append(":")
        resLine.append(residue.getResidueNumber()).append("\n")
        bw.write(resLine.toString())

        // Get this residue's rotamers.
        Rotamer[] rotamers = residue.getRotamers(rotamerLibrary)

        assert rotamers != null && rotamers.length > 1

        // Configure "active" and "use" flags.
        // .getVariableAtoms returns backbone atoms for
        // nucleic acids, which is what we want
        List<Atom> sideChainAtoms = residue.getVariableAtoms()
        for (Atom atom : sideChainAtoms) {
          atom.setActive(true)
          atom.setUse(true)
        }

        // Define "all previously saved rotamers" arrayList to be used for RMSD comparison
        ArrayList<ResidueState> keptRotamers = new ArrayList<>()
        int keptRotamersCount = 0

        // Loop over rotamers for this Residue.
        for (int i = 0; i < rotamers.length; i++) {
          Rotamer rotamer = rotamers[i]

          // Apply the rotamer (i.e. amino acid side-chain or nucleic acid suite).
          RotamerLibrary.applyRotamer(residue, rotamer)

          if (i > 0 || !useOriginalRotamers) {
            // -Dplatform=omm
            MinimizationEngine engine =
                Minimize.defaultEngine(activeAssembly, activeAssembly.getPotentialEnergy())
            Minimize minimize = Minimize.minimizeFactory(activeAssembly,
                activeAssembly.getPotentialEnergy(), algorithmListener, engine)
            // Locally minimize.
            minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations())
          } else {
            logger.info(" Skipping minimization of original-coordinates rotamer.")
          }

          // Stores a copy of minimized residue coordinates
          // newResidueState will be added to keptRotamers if its RMSD
          // to all previously kept rotamers is above cut off (user defined or
          // default of 0.1 kcal/mol)
          ResidueState newResState = new ResidueState(residue)

          if (i == 0) {
            // Add 0th rotamer to the keptRotamers ArrayList
            keptRotamers.add(newResState)

            // Save out coordinates to a rotamer file (inputFileName.rot)
            bw.write(String.format("  ROT:%d\n", keptRotamersCount))
            for (Atom atom : sideChainAtoms) {
              double x = atom.getX()
              double y = atom.getY()
              double z = atom.getZ()
              logger.info(String.format(" %s %16.8f %16.8f %16.8f", atom.toString(), x, y, z))
              StringBuilder atomLine = new StringBuilder("   ATOM:")
              atomLine.append(atom.getName()).append(":")
              atomLine.append(x).append(":")
              atomLine.append(y).append(":")
              atomLine.append(z).append("\n")
              bw.write(atomLine.toString())
            }
            bw.write("  ENDROT\n")
            keptRotamersCount++
          } else {
            // For all but the 0th rotamer, do RMSD calculations to determine if the "new" rotamer
            // is within a cut off (default: 0.1 kcal/mol) of any other previously saved rotamer.
            logger.info("Number of rotamers kept for this residue: " + keptRotamers.size())

            // Define RMSD threshold value boolean
            boolean withinRange = false

            // Compare newResState (i.e.: residue with newly applied rotamer) to all
            // previously saved rotamers in keptRotamers
            for (int k = 0; k < keptRotamers.size(); k++) {
              double RMSD = newResState.compareTo(keptRotamers[k])
              logger.info("RMSD: " + RMSD + "\n")
              if (RMSD <= rmsdCutoff) {
                withinRange = true
              }
            }

            // Keep new rotamer if and only if it is different from all previously
            // kept rotamers by an RMSD of greater than the user-defined cut off
            // (or default cut off of 0.1 kcal/mol)
            if (withinRange) {
              // Rotamer not written because it's too energetically similar to another rotamer in the set
              // Keeping too many similar rotamers causes problems with Dead End Elimination
              logger.info("Rotamer not kept")
            } else {
              // Save out coordinates to a rotamer file.
              bw.write(String.format("  ROT:%d\n", keptRotamersCount))
              for (Atom atom : sideChainAtoms) {
                double x = atom.getX()
                double y = atom.getY()
                double z = atom.getZ()
                logger.info(String.format(" %s %16.8f %16.8f %16.8f", atom.toString(), x, y, z))
                StringBuilder atomLine = new StringBuilder("   ATOM:")
                atomLine.append(atom.getName()).append(":")
                atomLine.append(x).append(":")
                atomLine.append(y).append(":")
                atomLine.append(z).append("\n")
                bw.write(atomLine.toString())
              }
              bw.write("  ENDROT\n")

              // Add the new rotamer to keptRotamers list
              keptRotamers.add(newResState)
              keptRotamersCount++
            }

          }
        }

        // Set the Residue conformation back to rotamer 0.
        RotamerLibrary.applyRotamer(residue, rotamers[0])

        // Revert the active and use flags.
        for (Atom atom : sideChainAtoms) {
          atom.setActive(false)
          atom.setUse(false)
        }
      }

    } finally {
      bw?.flush()
      bw?.close()
    }

    return this
  }
}