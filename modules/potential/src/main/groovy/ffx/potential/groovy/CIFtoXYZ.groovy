//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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

import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.CIFFilter
import static java.lang.String.format;
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The CIFtoXYZ script converts a CIF file to an XYZ file including atom types.
 * TODO: Move CIF parsing into a parsers.CIFFilter class.
 * <br>
 * Usage:
 * <br>
 * ffxc CIFtoXYZ &lt;filename.cif&gt; &lt;filename.pdb&gt;
 */
@Command(description = " Convert a single molecule CIF file to XYZ format.", name = "ffxc CIFtoXYZ")
class CIFtoXYZ extends PotentialScript {

  /**
   * --sg or --spaceGroupNumber Override the CIF space group.
   */
  @Option(names = ['--sg', '--spaceGroupNumber'], paramLabel = "-1", defaultValue = "-1",
      description = 'Override the CIF space group.')
  private int sgNum

  /**
   * --name or --spaceGroupName Override the CIF space group.
   */
  @Option(names = ['--name', '--spaceGroupName'], paramLabel = "", defaultValue = "",
      description = 'Override the CIF space group.')
  private String sgName

  /**
   * --bt or --bondTolerance Tolerance added to covalent radius to bond atoms.
   */
  @Option(names = ['--bt', '--bondTolerance'], paramLabel = "0.2", defaultValue = "0.2",
      description = 'Tolerance added to covalent radius to determine if atoms should be bonded.')
  private double bondTolerance

  /**
   * --fl or --fixLattice Override CIF parameters to satisfy lattice conditions.
   */
  @Option(names = ['--fl', '--fixLattice'], paramLabel = "false", defaultValue = "false",
          description = 'Override CIF parameters to satisfy lattice conditions (Otherwise error).')
  private boolean fixLattice

  /**
   * The final argument(s) should be a CIF file and an XYZ file with atom types.
   */
  @Parameters(arity = "1..2", paramLabel = "files",
      description = "A CIF file and an XYZ file.")
  List<String> filenames = null

  /**
   * CIFtoXYZ Constructor.
   */
  CIFtoXYZ() {
    this(new Binding())
  }

  /**
   * CIFtoXYZ Constructor.
   * @param binding Groovy Binding to use.
   */
  CIFtoXYZ(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  CIFtoXYZ run() {

    // Turn off CDK logging.
    System.setProperty("cdk.logging.level", "fatal")

    // Turn off non-bonded interactions.
    System.setProperty("vdwterm", "false")

    if (!init()) {
      return this
    }

    if (filenames != null && filenames.size() == 2) {
      MolecularAssembly[] assemblies = potentialFunctions.openAll(filenames.get(1))
      System.clearProperty("mpoleterm")
      setActiveAssembly(assemblies[0])
      CIFFilter cifFilter = new CIFFilter(filenames.toArray() as String[], activeAssembly, sgNum, sgName, bondTolerance, fixLattice, baseDir);
      cifFilter.readFile()
    } else {
      logger.info(helpString())
      return this
    }
    return this
  }

  @Override
  List<Potential> getPotentials() {
    return new ArrayList<Potential>()
  }
}