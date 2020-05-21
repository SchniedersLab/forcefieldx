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
package ffx.xray.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.xray.DiffractionData
import ffx.xray.cli.XrayOptions
import ffx.xray.parsers.DiffractionFile
import ffx.xray.parsers.MTZWriter
import ffx.xray.parsers.MTZWriter.MTZType
import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import java.util.stream.Collectors

/**
 * The X-ray ComputeFc script.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.ComputeFc [options] &lt;filename [file2...]&gt;
 */
@Command(description = " Write out computed structure factors.", name = "ffxc ComputeFc")
class ComputeFc extends AlgorithmsScript {

  @Mixin
  XrayOptions xrayOptions

  /**
   * ComputeFc constructor.
   */
  ComputeFc() {
    this(new Binding())
  }

  /**
   * ComputeFc constructor.
   * @param binding The Groovy Binding to use.
   */
  ComputeFc(Binding binding) {
    super(binding)
  }

  /**
   * One or more filenames.
   */
  @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Diffraction input files.")
  private List<String> filenames

  private MolecularAssembly[] assemblies
  private DiffractionData diffractiondata

  @Override
  ComputeFc run() {

    if (!init()) {
      return this
    }

    xrayOptions.init()

    String modelfilename
    if (filenames != null && filenames.size() > 0) {
      assemblies = algorithmFunctions.openAll(filenames.get(0))
      activeAssembly = assemblies[0]
      modelfilename = filenames.get(0)
    } else if (activeAssembly == null) {
      logger.info(helpString())
      return this
    } else {
      modelfilename = activeAssembly.getFile().getAbsolutePath()
    }

    // Load parsed X-ray properties.
    CompositeConfiguration properties = activeAssembly.getProperties()
    xrayOptions.setProperties(parseResult, properties)

    // Set up diffraction data (can be multiple files)
    List<DiffractionData> diffractionfiles = xrayOptions.processData(filenames, assemblies)

    diffractiondata = new DiffractionData(assemblies, assemblies[0].getProperties(),
        xrayOptions.solventModel,
        diffractionfiles.toArray(new DiffractionFile[diffractionfiles.size()]))

    logger.info("\n Running xray.ComputeFc on " + modelfilename)

    // Compute structure factors
    diffractiondata.computeAtomicDensity()

    // Output Fcs
    MTZWriter mtzwriter = new MTZWriter(diffractiondata.reflectionlist[0],
        diffractiondata.refinementdata[0],
        FilenameUtils.getBaseName(modelfilename) + "_fc.mtz", MTZType.FCONLY)

    mtzwriter.write()

    return this
  }

  @Override
  List<Potential> getPotentials() {
    if (assemblies == null) {
      return new ArrayList<Potential>()
    } else {
      return Arrays.stream(assemblies).
          filter {a -> a != null
          }.
          map {a -> a.getPotentialEnergy()
          }.
          filter {e -> e != null
          }.
          collect(Collectors.toList())
    }
  }

  @Override
  boolean destroyPotentials() {
    return diffractiondata == null ? true : diffractiondata.destroy()
  }
}
