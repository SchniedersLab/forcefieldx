//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2019.
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

package ffx.potential.groovy.test

import static java.lang.String.format

import org.apache.commons.io.FilenameUtils
import static org.apache.commons.math3.util.FastMath.cos

import ffx.crystal.Crystal
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

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
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'XYZ file to be converted.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    /**
     * Execute the script.
     */
    @Override
    XYZtoQE run() {

        if (!init()) {
            return this
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = [potentialFunctions.open(filenames.get(0))]
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()

        activeAssembly.computeFractionalCoordinates()

        Atom[] atoms = activeAssembly.getAtomArray()

        File saveDir = baseDir
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.getAbsolutePath()
        String fileName = FilenameUtils.getName(modelFilename)
        fileName = FilenameUtils.removeExtension(fileName) + ".in"
        File modelFile = new File(dirName + File.separator + fileName)

        Crystal crystal = activeAssembly.getCrystal().getUnitCell()
        double xtalA = crystal.a
        double xtalB = crystal.b
        double xtalC = crystal.c

        HashMap<String, Double> atomTypes = new HashMap<String, Double>()
        String atomicPositions = ""
        for (atom in atoms) {
            if (!atomTypes.containsKey(atom.name)) {
                atomTypes.put(atom.name, atom.getAtomType().atomicWeight)
            }
            double[] xyz = atom.getXYZ(null)
            crystal.toFractionalCoordinates(xyz, xyz)
            atomicPositions += format("%2s %16.12f %16.12f %16.12f\n", atom.name, xyz[0], xyz[1], xyz[2])
        }

        BufferedWriter bwQE = new BufferedWriter(new FileWriter(modelFile))
        bwQE.write("&CONTROL\n" +
                "\tcalculation = 'vc-relax',\n" +
                "\trestart_mode = 'from_scratch',\n" +
                "\tprefix = 'odh',\n" +
                "\tetot_conv_thr = 2.0D-6,\n" +
                "\tforc_conv_thr = 6.0D-4,\n" +
                "\tnstep = 500,\n" +
                "/\n")
        bwQE.write("&SYSTEM\n" +
                "\tspace_group = " + crystal.spaceGroup.number + ",\n" +
//                "\tibrav = 0,\n" +
                "\tnat = " + activeAssembly.getAtomList().size() + ",\n" + //Number of atoms
                "\tntyp = " + atomTypes.size() + ",\n" +  //Number of different types of atoms (i.e. H, C, O = 3)
                "\ta = " + format("%16.12f", xtalA) + "\n" +
                "\tb = " + format("%16.12f", xtalB) + "\n" +
                "\tc = " + format("%16.12f", xtalC) + "\n" +
                "\tcosAB = " + format("%16.12f", cos(crystal.gamma)) + "\n" +
                "\tcosAC = " + format("%16.12f", cos(crystal.beta)) + "\n" +
                "\tcosBC = " + format("%16.12f", cos(crystal.alpha)) + "\n" +
                "\tecutwfc = 50.0,\n" +
                "\tecutrho = 500.0,\n" +
                "\tvdw_corr = 'XDM',\n" +
                "\txdm_a1 = 0.6512,\n" +
                "\txdm_a2 = 1.4633,\n" +
                "/\n")
        bwQE.write("&ELECTRONS\n" +
                "\telectron_maxstep = 1500,\n" +
                "\tconv_thr = 1.D-8,\n" +
                "\tscf_must_converge = .TRUE.,\n" +
                "\tmixing_beta = 0.5D0,\n" +
                "/\n")
        bwQE.write("&IONS\n" +
                "\tion_dynamics = 'bfgs',\n" +
                "/\n")
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

//        bwQE.write("CELL_PARAMETERS angstrom\n")
//        bwQE.write(format("\t%16.12f %16.12f %16.12f\n", crystal.Ai[0][0], crystal.Ai[0][1], crystal.Ai[0][2]))
//        bwQE.write(format("\t%16.12f %16.12f %16.12f\n", crystal.Ai[1][0], crystal.Ai[1][1], crystal.Ai[1][2]))
//        bwQE.write(format("\t%16.12f %16.12f %16.12f\n", crystal.Ai[2][0], crystal.Ai[2][1], crystal.Ai[2][2]))

        bwQE.write("ATOMIC_POSITIONS crystal_sg\n" +
                atomicPositions + "\n")
        //Set-Up K_Points Card
        int k1
        int k2
        int k3
        if (xtalA < 5) {
            k1 = 8
        } else if (xtalA <= 8) {
            k1 = 6
        } else if (xtalA <= 12) {
            k1 = 3
        } else {
            k1 = 1
        }

        if (xtalB < 5) {
            k2 = 8
        } else if (xtalB <= 8) {
            k2 = 6
        } else if (xtalB <= 12) {
            k2 = 3
        } else {
            k2 = 1
        }

        if (xtalC < 5) {
            k3 = 8
        } else if (xtalC <= 8) {
            k3 = 6
        } else if (xtalC <= 12) {
            k3 = 3
        } else {
            k3 = 1
        }

        bwQE.write("K_POINTS automatic\n" + k1 + " " + k2 + " " + k3 + " 1 1 1\n")

        bwQE.close()
        return this
    }
}