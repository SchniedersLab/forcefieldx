/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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
package ffx.ui;

import java.io.File;

import ffx.algorithms.AlgorithmFunctions;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;

/**
 * The UIUtils class implements all Function interfaces, enabling lower modules
 * to be blind to the existence of User Interfaces but still update the GUI and
 * be in communication with User Interfaces.
 *
 * It presently gets away with simply implementing AlgorithmFunctions; if the
 * Automatic Parametrization method closures are similarly re-implemented, it
 * will probably be necessary to implement a new interface which extends both
 * AlgorithmFunctions and AutomaticParamFunctions interfaces.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class UIUtils implements AlgorithmFunctions {

    private final ModelingShell modelingShell;
    private final MainPanel mainPanel;

    public UIUtils(ModelingShell modelingShell, MainPanel mainPanel) {
        this.modelingShell = modelingShell;
        this.mainPanel = mainPanel;
    }

    @Override
    public void md(MolecularAssembly assembly, int nStep, double timeStep, double printInterval, double saveInterval, double temperature, boolean initVelocities, File dyn) {
        modelingShell.md(nStep, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
    }

    @Override
    public Potential minimize(MolecularAssembly assembly, double eps) {
        return modelingShell.minimize(eps);
    }

    @Override
    public boolean isLocal() {
        return false;
    }

    @Override
    public FFXSystem[] open(String file) {
        return mainPanel.openWait(file);
    }

    @Override
    public FFXSystem[] open(String[] files) {
        return mainPanel.openWait(files);
    }

    @Override
    public void close(MolecularAssembly assembly) {
        mainPanel.closeWait();
    }

    @Override
    public void closeAll(MolecularAssembly[] assemblies) {
        mainPanel.closeAll();
    }

    @Override
    public double time() {
        return modelingShell.time();
    }

    @Override
    public void save(MolecularAssembly assembly, File file) {
        mainPanel.saveAsXYZ(file);
    }

    @Override
    public void saveAsXYZ(MolecularAssembly assembly, File file) {
        mainPanel.saveAsXYZ(file);
    }

    @Override
    public void saveAsP1(MolecularAssembly assembly, File file) {
        mainPanel.saveAsP1(file);
    }

    @Override
    public void saveAsPDB(MolecularAssembly assembly, File file) {
        mainPanel.saveAsPDB(file);
    }

    @Override
    public void saveAsPDB(MolecularAssembly[] assemblies, File file) {
        mainPanel.saveAsPDB(assemblies, file);
    }
    
    @Override
    public void savePDBSymMates(MolecularAssembly assembly, File file) {
        mainPanel.savePDBSymMates(file, "_symMate");
    }
    
    @Override
    public void savePDBSymMates(MolecularAssembly assembly, File file, String suffix) {
        mainPanel.savePDBSymMates(file, suffix);
    }

    @Override
    public ForceFieldEnergy energy(MolecularAssembly assembly) {
        return modelingShell.energy();
    }

    @Override
    public double returnEnergy(MolecularAssembly assembly) {
        return modelingShell.returnEnergy();
    }

    @Override
    public MolecularAssembly[] convertDataStructure(Object data) {
        return mainPanel.convertWait(data, null);
    }

    @Override
    public MolecularAssembly[] convertDataStructure(Object data, File file) {
        return mainPanel.convertWait(data, file);
    }
    
    @Override
    public MolecularAssembly[] convertDataStructure(Object data, String filename) {
        File file = new File(filename);
        if (!file.exists() || file.isDirectory() || !file.canRead()) {
            throw new IllegalArgumentException(String.format("%s not a valid file name.", filename));
        }
        return mainPanel.convertWait(data, file);
    }
}
