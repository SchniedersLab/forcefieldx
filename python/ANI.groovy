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

import ffx.potential.bonded.Atom
import ffx.potential.cli.PotentialScript
import org.graalvm.polyglot.Context
import org.graalvm.polyglot.Value
import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import java.nio.file.Path
import java.nio.file.Paths

/**
 * The ANI script evaluates the ANI energy of a system.
 * This must be executed using GraalVM and a virtual environment with Torch installed.
 * <br>
 * Usage:
 * <br>
 * ffxc ANI &lt;filename&gt;
 */
@Command(description = " Compute the ANI-2x energy.", name = "ANI-2x")
class ANI extends PotentialScript {

    /**
     * The final argument is a PDB or XYZ coordinate file.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private String filename = null

    /**
     * Energy constructor.
     */
    ANI() {
        this(new Binding())
    }

    /**
     * Energy constructor.
     * @param binding The Groovy Binding to use.
     */
    ANI(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    ANI run() {
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

        logger.info("\n Running Energy on " + filename)

        String FFX_HOME = System.getProperty("basedir")
        Path graalpy = Paths.get(FFX_HOME, "ffx_venv", "bin", "graalpy")
        graalpyString = System.getProperty("graalpy", graalpy.toString())
        graalpy = Paths.get(graalpyString)
        logger.info(" graalpy (-Dgraalpy=path.to.graalpy):             " + graalpy)
        String torchScript = "ANI2x.pt"
        torchScript = System.getProperty("torchscript", torchScript)
        logger.info(" torchscript (-Dtorchscript=path.to.torchscript): " + torchScript)
        // Collect atomic number and coordinates for each atom.
        Atom[] atoms = activeAssembly.getAtomArray()
        int nAtoms = atoms.length
        int[] jspecies = new int[nAtoms]
        double[][] jcoords = new double[nAtoms][3]
        for (int i = 0; i < nAtoms; i++) {
            Atom a = atoms[i]
            jspecies[i] = a.getAtomicNumber()
            jcoords[i][0] = a.getX()
            jcoords[i][1] = a.getY()
            jcoords[i][2] = a.getZ()
        }
        double energy = 0.0
        double[] grad = new double[nAtoms * 3]
        // Construct a Polyglot Python environment.
        try (Context context = Context.newBuilder("python").allowAllAccess(true).
            option("python.Executable", graalpy.toString()).build()) {
            // Place the coords and species arrays into the context.
            Value polyglotBindings = context.getPolyglotBindings()
            polyglotBindings.putMember("jcoords", jcoords)
            polyglotBindings.putMember("jspecies", jspecies)

            // Construct the Python code to run ANI-2x using TorchScript.
            String torch = "import site\n"
            torch += "import polyglot\n"
            torch += "import torch\n"
            // Load the Java arrays into the Torch tensors.
            torch += "jspecies = polyglot.import_value('jspecies')\n"
            torch += "jcoords = polyglot.import_value('jcoords')\n"
            torch += "species = torch.tensor([jspecies], dtype=torch.int64)\n"
            torch += "coordinates = torch.tensor([jcoords], dtype=torch.double)\n"
            // Load and evaluate the ANI-2x forward method.
            torch += "ani = torch.jit.load('" + torchScript + "')\n"
            torch += "gradient = ani(species, coordinates)\n"

            // Evaluate ANI-2x and collect the energy and gradient.
            Value result = context.eval("python", torch)
            Value bindings = context.getBindings("python")
            Value ret = bindings.getMember("gradient")
            for (int i = 0; i < nAtoms * 3; i++) {
                grad[i] = ret.getArrayElement(i).asDouble()
            }
            energy = ret.getArrayElement(nAtoms * 3).asDouble()
            logger.info(" ANI-2x Energy (Hartree): " + energy)
            // Context close does not return.
            // context.close(true)
        } catch (Exception e) {
            logger.info(" Exception:\n" + e.toString())
        }

        return this
    }
}

