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

import ffx.crystal.Crystal
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.*
import static org.apache.commons.math3.util.FastMath.cos

/**
 * Automatic generation of Quantum Espresso input files based off of a XYZ file.
 *
 * @author Aaron J. Nessler
 * @author Michael J. Schnieders
 * <br>
 * Usage:
 * <br>
 * ffxc XYZtoQE [options] &lt;filename&gt;
 */
@Command(description = "Generate QE input from a XYZ file.", name = "ffxc XYZtoQE")
class XYZtoQE extends PotentialScript {

  /**
   * --ns or --numSteps Number of structural optimization steps performed in this run.
   */
  @Option(names = ['--ns', '--nstep'], paramLabel = "500", defaultValue = "500",
      description = 'Number of structural optimization steps performed in this run.')
  private int nStep

  /**
   * --ec or --etot_conv_thr Convergence threshold on total energy (a.u) for ionic minimization.
   */
  @Option(names = ['--ec', '--etot_conv_thr'], paramLabel = "1.0e-6", defaultValue = "1.0e-6",
      description = 'Convergence threshold on total energy (a.u) for ionic minimization.')
  private double etotConvThr

  /**
   * --ef or --forc_conv_thr Convergence threshold on forces (a.u) for ionic minimization.
   */
  @Option(names = ['--ef', '--forc_conv_thr'], paramLabel = "1.0e-4", defaultValue = "1.0e-4",
      description = 'Convergence threshold on forces (a.u) for ionic minimization.')
  private double forcConvThr

  /**
   * --ke or --ecutwfc Kinetic energy cutoff (Ry) for wavefunctions.
   */
  @Option(names = ['--ke', '--ecutwfc'], paramLabel = "50.0", defaultValue = "50.0",
      description = 'Kinetic energy cutoff (Ry) for wavefunctions.')
  private double ecutwfc

  /**
   * --rho or --ecutrho Kinetic energy cutoff (Ry) for charge density and potential.
   */
  @Option(names = ['--rho', '--ecutrho'], paramLabel = "500.0", defaultValue = "500.0",
      description = 'Kinetic energy cutoff (Ry) for charge density and potential.')
  private double ecutrho

  /**
   * --em or --electron_maxstep Maximum number of iterations in a scf step.
   */
  @Option(names = ['--em', '--electron_maxstep'], paramLabel = "1500", defaultValue = "1500",
      description = 'Maximum number of iterations in a scf step.')
  private int electronMaxstep

  /**
   * --ct or --convThr Convergence threshold for self consistency.
   */
  @Option(names = ['--ct', '--conv_thr'], paramLabel = "1.0e-8", defaultValue = "1.0e-8",
      description = 'Convergence threshold for self consistency.')
  private double convThr

  /**
   * --mb or --mixing_beta Mixing factor for self-consistency.
   */
  @Option(names = ['--mb', '--mixing_beta'], paramLabel = "0.5", defaultValue = "0.5",
      description = 'Mixing factor for self-consistency.')
  private double mixingBeta

  /**
   * --hx or --hexagonal Perform QE calculation on hexagonal rather than rhombohedral representation
   */
  @Option(names = ['--hx', '--hexagonal'], paramLabel = "true", defaultValue = "true",
      description = 'Perform QE on hexagonal system.')
  private boolean hexagonal

  /**
   * The final argument(s) should be one or more filenames.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = 'XYZ file to be converted.')
  String filename = null

  /**
   * XYZtoQE Constructor.
   */
  XYZtoQE() {
    this(new Binding())
  }

  /**
   * XYZtoQE Constructor.
   * @param binding Groovy Binding to use.
   */
  XYZtoQE(Binding binding) {
    super(binding)
  }

  /**
   * Execute the script.
   */
  @Override
  XYZtoQE run() {

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
    logger.info(format("\n Converting %s to QE format\n", filename))

    File saveDir = baseDir
    if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
      saveDir = new File(getFullPath(filename))
    }
    String dirName = saveDir.getAbsolutePath()
    String name = getName(filename)
    name = removeExtension(name)
    File modelFile = new File(dirName + File.separator + name + ".in")

    Crystal crystal = activeAssembly.getCrystal().getUnitCell()
    double xtalA = crystal.a
    double xtalB = crystal.b
    double xtalC = crystal.c

    HashMap<String, Double> atomTypes = new HashMap<>()
    String atomicPositions = ""

    Atom[] atoms = activeAssembly.getAtomArray()
    for (atom in atoms) {
      if (!atomTypes.containsKey(atom.name)) {
        atomTypes.put(atom.name, atom.getAtomType().atomicWeight)
      }
      double[] xyz = atom.getXYZ(null)
      crystal.toFractionalCoordinates(xyz, xyz)
      atomicPositions += format("%2s %16.12f %16.12f %16.12f\n", atom.name, xyz[0], xyz[1], xyz[2])
    }

    BufferedWriter bwQE = new BufferedWriter(new FileWriter(modelFile))
    // general variables controlling the run
    bwQE.write(format("&CONTROL\n" +
        "\tcalculation = 'vc-relax',\n" +
        "\trestart_mode = 'from_scratch',\n" +
        "\tprefix = '%s',\n" +
        "\tetot_conv_thr = %6.4E,\n" +
        "\tforc_conv_thr = %6.4E,\n" +
        "\tnstep = %d,\n" +
        "/\n", name, etotConvThr, forcConvThr, nStep))
    // structural information on the system under investigation
    bwQE.write(format("&SYSTEM\n" +
        "\tspace_group = " + crystal.spaceGroup.number + ",\n" +
        // Number of atoms
        "\tnat = " + activeAssembly.getAtomList().size() + ",\n" +
        // Number of different types of atoms (i.e. H, C, O = 3)
        "\tntyp = " + atomTypes.size() + ",\n" +
        "\ta = " + format("%16.12f", xtalA) + "\n" +
        "\tb = " + format("%16.12f", xtalB) + "\n" +
        "\tc = " + format("%16.12f", xtalC) + "\n" +
        "\tcosAB = " + format("%16.12f", cos(crystal.gamma)) + "\n" +
        "\tcosAC = " + format("%16.12f", cos(crystal.beta)) + "\n" +
        "\tcosBC = " + format("%16.12f", cos(crystal.alpha)) + "\n" +
        "\tecutwfc = %6.4f,\n" +
        "\tecutrho = %6.4f,\n" +
        "\tvdw_corr = 'XDM',\n" +
        // xdm_a1 and xdm_a2 do not need to be provided for B86bPBE
        // "\txdm_a1 = 0.6512,\n" +
        // "\txdm_a2 = 1.4633,\n" +
        "/\n", ecutwfc, ecutrho))
    if (hexagonal) {
      bwQE.write("\trhombohedral = .FALSE.,\n")
    }
    bwQE.write("/\n")
    // electronic variables: self-consistency, smearing
    bwQE.write(format("&ELECTRONS\n" +
        "\telectron_maxstep = %d,\n" +
        "\tconv_thr = %6.4E,\n" +
        "\tscf_must_converge = .TRUE.,\n" +
        "\tmixing_beta = %5.3f,\n" +
        "/\n", electronMaxstep, convThr, mixingBeta))
    // ionic variables: relaxation, dynamics
    bwQE.write("&IONS\n" +
        "\tion_dynamics = 'bfgs',\n" +
        "/\n")
    // variable-cell optimization or dynamics
    bwQE.write("&CELL\n" +
        "\tcell_dynamics = 'bfgs',\n" +
        "/\n")

    String line = ""
    Iterator it = atomTypes.entrySet().iterator()
    while (it.hasNext()) {
      Map.Entry pair = (Map.Entry) it.next()
      line += " " + pair.getKey() + " " + pair.getValue() + " " + pair.getKey() + ".b86bpbe.UPF\n"
    }
    bwQE.write("ATOMIC_SPECIES\n" + line + "\n")

    bwQE.write("ATOMIC_POSITIONS crystal_sg\n" +
        atomicPositions + "\n")

    // Set-Up K_Points Card
    // TODO: add explanations for k-point choices.
    int k1
    int k2
    int k3
    if (xtalA < 5) {
      k1 = 8
    } else if (xtalA <= 7) {
      k1 = 6
    } else if (xtalA <= 12) {
      k1 = 4
    } else {
      k1 = 2
    }

    if (xtalB < 5) {
      k2 = 8
    } else if (xtalB <= 7) {
      k2 = 6
    } else if (xtalB <= 12) {
      k2 = 4
    } else {
      k2 = 2
    }

    if (xtalC < 5) {
      k3 = 8
    } else if (xtalC <= 7) {
      k3 = 6
    } else if (xtalC <= 12) {
      k3 = 4
    } else {
      k3 = 2
    }

    bwQE.write("K_POINTS automatic\n" + k1 + " " + k2 + " " + k3 + " 1 1 1\n")
    bwQE.close()

    logger.info(format(" Saved QE file: %s", modelFile))

    return this
  }
}