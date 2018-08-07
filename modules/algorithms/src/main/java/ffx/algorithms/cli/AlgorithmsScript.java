/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms.cli;

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmUtils;
import ffx.algorithms.AlgorithmListener;
import ffx.utilities.BaseScript;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Base class for scripts in the Algorithms package, providing some key functions.
 *
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class AlgorithmsScript extends BaseScript {

    /**
     * An instance of AlgorithmFunctions passed into the current context.
     */
    public AlgorithmFunctions algorithmFunctions;

    /**
     * An active MolecularAssembly passed into the current context or loaded by
     * the Script from a file argument.
     */
    public MolecularAssembly activeAssembly;

    /**
     * An instance of the AlgorithmListener interface.
     */
    public AlgorithmListener algorithmListener;

    /**
     * The directory in which to place output files. Mostly for tests.
     */
    protected File saveDir;

    /**
     * Execute the BaseScript init method, then load algorithm functions.
     *
     * @return Returns true if the script should continue.
     */
    @Override
    public boolean init() {
        if (!super.init()) {
            return false;
        }

        if (context.hasVariable("functions")) {
            algorithmFunctions = (AlgorithmFunctions) context.getVariable("functions");
        } else {
            algorithmFunctions = new AlgorithmUtils();
        }

        activeAssembly = null;
        if (context.hasVariable("active")) {
            activeAssembly = (MolecularAssembly) context.getVariable("active");
        }

        algorithmListener = null;
        if (context.hasVariable("listener")) {
            algorithmListener = (AlgorithmListener) context.getVariable("listener");
        }

        return true;
    }

    /**
     * Returns a List of all Potential objects associated with this script.
     *
     * @return All Potentials. Sometimes empty, never null.
     */
    public List<Potential> getPotentials() {
        List<Potential> plist = new ArrayList<>();
        if (activeAssembly != null && activeAssembly.getPotentialEnergy() != null) {
            plist.add(activeAssembly.getPotentialEnergy());
        }
        return plist;
    }

    /**
     * Reclaims resources associated with all Potential objects associated with this script.
     *
     * @return If all Potentials had resources reclaimed.
     */
    public boolean destroyPotentials() {
        boolean allSucceeded = true;
        for (Potential potent : getPotentials()) {
            logger.fine(String.format(" Potential %s is being destroyed. ", potent));
            allSucceeded = allSucceeded && potent.destroy();
        }
        return allSucceeded;
    }

    /**
     * Sets the directory this script should save files to. Mostly used for tests.
     * @param saveDir Directory to save output to.
     */
    public void setSaveDir(File saveDir) {
        this.saveDir = saveDir;
    }

    /**
     * Gets a File in the save directory with the same name as the input file. Can just be the original
     * file if saveDir was never set, which is the case for production runs.'
     *
     * @param file File to find a save location for.
     * @return File to save to
     */
    protected File saveDirFile(File file) {
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            return file;
        } else {
            String baseName = file.getName();
            String newName = saveDir.getAbsolutePath() + File.separator + baseName;
            return new File(newName);
        }
    }
}
