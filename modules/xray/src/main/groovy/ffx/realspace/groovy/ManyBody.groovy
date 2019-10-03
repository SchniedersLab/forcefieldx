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
package ffx.realspace.groovy

import org.apache.commons.configuration2.CompositeConfiguration
import org.apache.commons.io.FilenameUtils

import edu.rit.pj.Comm

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.numerics.Potential
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import ffx.realspace.RealSpaceData
import ffx.realspace.cli.RealSpaceOptions
import ffx.realspace.parsers.RealSpaceFile
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The ManyBody script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Discrete optimization using a many-body expansion and elimination expressions.", name = "ffxc realspace.ManyBody")
class ManyBody extends AlgorithmsScript {

    @Mixin
    RealSpaceOptions realSpace

    @Mixin
    ManyBodyOptions manyBody

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames
    private RefinementEnergy refinementEnergy;

    @Override
    ManyBody run() {

        if (!init()) {
            return this
        }

        String filename
        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
            filename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            filename = activeAssembly.getFile().getAbsolutePath()
        }

        CompositeConfiguration properties = activeAssembly.getProperties();
        if (!properties.containsKey("gk-suppressWarnings")) {
            properties.setProperty("gk-suppressWarnings", "true")
        }

        activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)

        List<RealSpaceFile> mapFiles = realSpace.processData(filenames, assemblies)
        RealSpaceData realSpaceData = new RealSpaceData(activeAssembly, activeAssembly.getProperties(),
                activeAssembly.getParallelTeam(), mapFiles.toArray(new RealSpaceFile[0]))

        refinementEnergy = new RefinementEnergy(realSpaceData, RefinementMinimize.RefinementMode.COORDINATES, null);

        RotamerOptimization rotamerOptimization = new RotamerOptimization(
                activeAssembly, refinementEnergy, algorithmListener)

        manyBody.initRotamerOptimization(rotamerOptimization, activeAssembly)

        double[] x = new double[refinementEnergy.getNumberOfVariables()]
        x = refinementEnergy.getCoordinates(x)
        refinementEnergy.energy(x, true)

        ArrayList<Residue> residueList = rotamerOptimization.getResidues()
        RotamerLibrary.measureRotamers(residueList, false);

        if (manyBody.algorithm == 1) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.INDEPENDENT);
        } else if (manyBody.algorithm == 2) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.ALL);
        } else if (manyBody.algorithm == 3) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BRUTE_FORCE);
        } else if (manyBody.algorithm == 4) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.WINDOW);
        } else if (manyBody.algorithm == 5) {
            rotamerOptimization.optimize(RotamerOptimization.Algorithm.BOX);
        }

        boolean master = true;
        if (Comm.world().size() > 1) {
            int rank = Comm.world().rank();
            if (rank != 0) {
                master = false;
            }
        }

        if (master) {
            logger.info(" Final Minimum Energy")

            algorithmFunctions.energy(activeAssembly)

            String ext = FilenameUtils.getExtension(filename);
            filename = FilenameUtils.removeExtension(filename);
            if (ext.toUpperCase().contains("XYZ")) {
                algorithmFunctions.saveAsXYZ(assemblies, new File(filename + ".xyz"))
            } else {
                algorithmFunctions.saveAsPDB(assemblies, new File(filename + ".pdb"))
            }
        }

        // manyBody.saveEliminatedRotamers()

        return this
    }

    @Override
    List<Potential> getPotentials() {
        return refinementEnergy == null ? Collections.emptyList() : Collections.singletonList(refinementEnergy);
    }
}
