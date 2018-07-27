package ffx.algorithms.groovy

import org.apache.commons.io.FilenameUtils

import edu.rit.pj.Comm

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.potential.ForceFieldEnergy
import ffx.algorithms.RotamerOptimization;
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary

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
@Command(description = " Run ManyBody algorithm on a system.", name = "ffxc ManyBody")
class ManyBody extends AlgorithmsScript {

    @Mixin
    ManyBodyOptions manyBody;

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = "PDB input file.")
    private List<String> filenames;

    private File baseDir = null;
    boolean testing = null;
    
    ForceFieldEnergy potentialEnergy;

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir;
    }
    
    @Override
    ManyBody run() {

        if (!init()) {
            return this;
        }

        String priorGKwarn = System.getProperty("gk-suppressWarnings");
        if (priorGKwarn == null || priorGKwarn.isEmpty()) {
            System.setProperty("gk-suppressWarnings", "true");
        }
        String modelFileName;
        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0));
            activeAssembly = assemblies[0];
            modelFileName = activeAssembly.getFile().getAbsolutePath();
        } else if (activeAssembly == null) {
            logger.info(helpString());
            return this;
        } else {
            logger.warning("Could not load file or active assembly.");
        }
        activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false);
        potentialEnergy = activeAssembly.getPotentialEnergy();

        RotamerOptimization rotamerOptimization = new RotamerOptimization(
            activeAssembly, activeAssembly.getPotentialEnergy(), algorithmListener);
        testing  = getTesting();

        if(testing == true){
            rotamerOptimization.turnRotamerSingleEliminationOff();
            rotamerOptimization.turnRotamerPairEliminationOff();
        }
        manyBody.initRotamerOptimization(rotamerOptimization, activeAssembly);

        ArrayList<Residue> residueList = rotamerOptimization.getResidues();

        boolean master = true;
        if (Comm.world().size() > 1) {
            int rank = Comm.world().rank();
            if (rank != 0) {
                master = false;
            }
        }

        algorithmFunctions.energy(activeAssembly)

        RotamerLibrary.measureRotamers(residueList, false);

        RotamerOptimization.Algorithm algo;
        switch (manyBody.getAlgorithmNumber()) {
            case 1:
                algo = RotamerOptimization.Algorithm.INDEPENDENT;
                break;
            case 2:
                algo = RotamerOptimization.Algorithm.ALL;
                break;
            case 3:
                algo = RotamerOptimization.Algorithm.BRUTE_FORCE;
                break;
            case 4:
                algo = RotamerOptimization.Algorithm.WINDOW;
                break;
            case 5:
                algo = RotamerOptimization.Algorithm.BOX;
                break;
            default:
                throw new IllegalArgumentException(String.format(" Algorithm choice was %d, not in range 1-5!", manyBody.getAlgorithmNumber()));
        }
        rotamerOptimization.optimize(algo);

        if (master) {
            logger.info(" Final Minimum Energy");

            algorithmFunctions.energy(activeAssembly);

            File saveDir = baseDir;
            if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
                saveDir = new File(FilenameUtils.getFullPath(modelFileName));
            }
            String dirName = saveDir.toString() + File.separator;
            String fileName = FilenameUtils.getName(modelFileName);
            fileName = FilenameUtils.removeExtension(fileName) + ".pdb";
            File modelFile = new File(dirName + fileName);

            algorithmFunctions.saveAsPDB(activeAssembly, modelFile);
        }

        manyBody.saveEliminatedRotamers();

        if (priorGKwarn == null) {
            System.clearProperty("gk-suppressWarnings");
        }

        return this;
    }

    /**
     * Returns the potential energy of the active assembly. Used during testing assertions.
     * @return potentialEnergy Potential energy of the active assembly.
     */
    ForceFieldEnergy getPotential(){
        return potentialEnergy;
    }


    /**
     * Set method for the testing boolean. When true, the testing boolean will shut off all elimination criteria forcing either a monte carlo or brute force search over all permutations.
     * @param testing A boolean flag that turns off elimination criteria for testing purposes.
     */
    void setTesting(boolean testing){
        this.testing = testing;
    }

    /**
     * Get method for the testing boolean. When true, the testing boolean will shut off all elimination criteria forcing either a monte carlo or brute force search over all permutations.
     * @return testing A boolean flag that turns off elimination criteria for testing purposes.
     */
    boolean getTesting(){
        return testing;
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
