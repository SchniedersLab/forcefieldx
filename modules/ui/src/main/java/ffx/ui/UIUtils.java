/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.ui;

import java.io.File;

import ffx.algorithms.AlgorithmFunctions;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.bonded.MolecularAssembly;

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
    public void close() {
        mainPanel.closeWait();
    }

    @Override
    public void closeAll() {
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
    public ForceFieldEnergy energy(MolecularAssembly assembly) {
        return modelingShell.energy();
    }

    @Override
    public double returnEnergy(MolecularAssembly assembly) {
        return modelingShell.returnEnergy();
    }
}
