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
package ffx.xray.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.xray.DiffractionData
import ffx.xray.cli.XrayOptions
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.removeExtension

/**
 * The X-ray ModelvsData script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.ModelvsData [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Compare the PDB model to the diffraction data.", name = "xray.ModelvsData")
class ModelvsData extends AlgorithmsScript {

  @Mixin
  XrayOptions xrayOptions

  /**
   * -m or --maps Output sigmaA weighted 2Fo-Fc and Fo-Fc electron density maps.
   */
  @Option(names = ['-p', '--maps'], paramLabel = 'false', description = 'Output sigmaA weighted 2Fo-Fc and Fo-Fc electron density maps.')
  boolean maps = false
  /**
   * -t or --timings Perform FFT timings.
   */
  @Option(names = ['-t', '--timings'], paramLabel = 'false', description = 'Perform FFT timings.')
  boolean timings = false
  /**
   * -w or --mtz Write out MTZ containing structure factor coefficients.
   */
  @Option(names = ['-w', '--mtz'], paramLabel = 'false',
      description = 'write out MTZ containing structure factor coefficients.')
  boolean mtz = false

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
  private List<String> filenames
  private DiffractionData diffractionData
  private MolecularAssembly[] molecularAssemblies

  /**
   * ModelvsData constructor.
   */
  ModelvsData() {
    this(new Binding())
  }

  /**
   * ModelvsData constructor.
   * @param binding The Groovy Binding to use.
   */
  ModelvsData(Binding binding) {
    super(binding)
  }

  @Override
  ModelvsData run() {

    if (!init()) {
      return this
    }

    xrayOptions.init()

    String filename
    if (filenames != null && filenames.size() > 0) {
      // Each alternate conformer is returned in a separate MolecularAssembly.
      molecularAssemblies = algorithmFunctions.openAll(filenames.get(0))
      activeAssembly = molecularAssemblies[0]
      filename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      molecularAssemblies = [activeAssembly]
      filename = activeAssembly.getFile().getAbsolutePath()
    }

    logger.info(format("\n Running xray.ModelvsData on %s", filename))

    // Combine script flags (in parseResult) with properties.
    CompositeConfiguration properties = activeAssembly.getProperties()
    xrayOptions.setProperties(parseResult, properties)

    // Set up diffraction data (can be multiple files)
    diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, parseResult)
    diffractionData.scaleBulkFit()
    diffractionData.printStats()
    algorithmFunctions.energy(molecularAssemblies)

    if (mtz) {
      diffractionData.writeData(removeExtension(filename) + "_ffx.mtz")
    }

    if (maps) {
      diffractionData.writeMaps(removeExtension(filename) + "_ffx")
    }

    if (timings) {
      diffractionData.timings()
    }

    return this
  }

  @Override
  List<Potential> getPotentials() {
    if (molecularAssemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(molecularAssemblies).filter {a -> a != null
      }.map {a -> a.getPotentialEnergy()
      }.filter {e -> e != null
      }.collect(Collectors.toList())
    }
  }

  @Override
  boolean destroyPotentials() {
    return diffractionData == null ? true : diffractionData.destroy()
  }
}
