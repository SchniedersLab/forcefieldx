// ******************************************************************************
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
// ******************************************************************************
package ffx.ui;

import static java.lang.String.format;

import ffx.algorithms.AlgorithmFunctions;
import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.AlgorithmUtils;
import ffx.numerics.Potential;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.PDBFilter.Mutation;
import ffx.potential.parsers.SystemFilter;
import ffx.utilities.Keyword;
import java.io.File;
import java.util.List;
import java.util.Optional;
import java.util.logging.Level;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * UIUtils implements core and extended functionality for many Force Field X algorithms and scripts,
 * such as opening and closing structure files, basic force field evaluations, molecular dynamics,
 * etc. This implementation performs additional tasks for the FFX graphical user interface, such as
 * updating the GUI and tree structure. This is widely used by our scripts. The AlgorithmUtils and
 * PotentialsUtils implementations of their respective interfaces lack the additional functionality
 * provided here.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class UIUtils extends AlgorithmUtils implements AlgorithmFunctions {

  private final ModelingShell modelingShell;
  private final MainPanel mainPanel;
  private SystemFilter lastFilter;

  public UIUtils(ModelingShell modelingShell, MainPanel mainPanel) {
    this.modelingShell = modelingShell;
    this.mainPanel = mainPanel;
  }

  @Override
  public void close(MolecularAssembly assembly) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    mainPanel.closeWait();
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public void closeAll(MolecularAssembly[] assemblies) {
    mainPanel.closeAll();
  }

  @Override
  public ForceFieldEnergy energy(MolecularAssembly assembly) {
    // TODO: Determine why this only runs energy on the last opened assembly, not the passed
    // assembly.
    Optional<FFXSystem> origSys = switchTo(assembly);
    ForceFieldEnergy ffe = modelingShell.energy();
    origSys.ifPresent(this::switchBack);
    return ffe;
  }

  @Override
  public MolecularAssembly getActiveAssembly() {
    return mainPanel.getHierarchy().getActive();
  }

  @Override
  public List<String> getArguments() {
    return modelingShell.getArgs();
  }

  @Override
  public AlgorithmListener getDefaultListener() {
    return modelingShell;
  }

  @Override
  public SystemFilter getFilter() {
    return lastFilter;
  }

  @Override
  public boolean isLocal() {
    return false;
  }

  @Override
  public void md(
      MolecularAssembly assembly,
      int nStep,
      double timeStep,
      double printInterval,
      double saveInterval,
      double temperature,
      boolean initVelocities,
      File dyn) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    modelingShell.md(
        nStep, timeStep, printInterval, saveInterval, temperature, initVelocities, dyn);
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public Potential minimize(MolecularAssembly assembly, double eps) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    Potential pot = modelingShell.minimize(eps);

    origSys.ifPresent(this::switchBack);
    return pot;
  }

  @Override
  public FFXSystem[] open(String[] files, int nThreads) {
    FFXSystem[] systems = mainPanel.openWait(files, nThreads);
    lastFilter = mainPanel.getFilter();
    return systems;
  }

  @Override
  public FFXSystem open(String file) {
    FFXSystem[] systems = mainPanel.openWait(file);
    lastFilter = mainPanel.getFilter();
    if (systems == null) {
      return null;
    }
    return systems[0];
  }

  @Override
  public FFXSystem[] openAll(String file) {
    FFXSystem[] systems = mainPanel.openWait(file);
    lastFilter = mainPanel.getFilter();
    return systems;
  }

  @Override
  public FFXSystem[] openAll(String[] files) {
    FFXSystem[] systems = mainPanel.openWait(files);
    lastFilter = mainPanel.getFilter();
    return systems;
  }

  @Override
  public FFXSystem[] openAll(String file, int nThreads) {
    FFXSystem[] systems = mainPanel.openWait(file, nThreads);
    lastFilter = mainPanel.getFilter();
    return systems;
  }

  /**
   * Mutates file on-the-fly as it is being opened. Used to open files for pHMD in fully-protonated
   * form.
   *
   * @param file a {@link java.io.File} object.
   * @param mutations a {@link java.util.List} object.
   * @return a {@link ffx.potential.MolecularAssembly} object.
   */
  @Override
  public MolecularAssembly openWithMutations(File file, List<Mutation> mutations) {
    file = new File(FilenameUtils.normalize(file.getAbsolutePath()));
    // Set the Current Working Directory based on this file.
    mainPanel.setCWD(file.getParentFile());

    // Create the CompositeConfiguration properties.
    CompositeConfiguration properties = Keyword.loadProperties(file);
    // Create an FFXSystem for this file.
    FFXSystem newSystem = new FFXSystem(file, "Open with mutations", properties);
    // Create a Force Field.
    ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
    ForceField forceField = forceFieldFilter.parse();
    String[] patches = properties.getStringArray("patch");
    for (String patch : patches) {
      logger.info(" Attempting to read force field patch from " + patch + ".");
      CompositeConfiguration patchConfiguration = new CompositeConfiguration();
      patchConfiguration.addProperty("parameters", patch);
      forceFieldFilter = new ForceFieldFilter(patchConfiguration);
      ForceField patchForceField = forceFieldFilter.parse();
      forceField.append(patchForceField);
      if (RotamerLibrary.addRotPatch(patch)) {
        logger.info(format(" Loaded rotamer definitions from patch %s.", patch));
      }
    }
    newSystem.setForceField(forceField);

    PDBFilter pdbFilter = new PDBFilter(file, newSystem, forceField, properties);
    pdbFilter.mutate(mutations);
    UIFileOpener opener = new UIFileOpener(pdbFilter, mainPanel);
    opener.run();
    lastFilter = pdbFilter;
    if (opener.getAllAssemblies().length > 1) {
      logger.log(
          Level.WARNING, "Found multiple assemblies in file {0}, opening first.", file.getName());
    }
    return opener.getAssembly();
  }

  @Override
  public double returnEnergy(MolecularAssembly assembly) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    double e = modelingShell.returnEnergy();
    origSys.ifPresent(this::switchBack);
    return e;
  }

  @Override
  public void save(MolecularAssembly assembly, File file) {
    saveAsXYZ(assembly, file);
  }

  @Override
  public void saveAsP1(MolecularAssembly assembly, File file) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    mainPanel.saveAsP1(file);
    lastFilter = mainPanel.getFilter();
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public void saveAsPDB(MolecularAssembly assembly, File file) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    mainPanel.saveAsPDB(file);
    lastFilter = mainPanel.getFilter();
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public void saveAsPDB(MolecularAssembly assembly, File file, boolean writeEnd, boolean append) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    mainPanel.saveAsPDB(file, false, append);
    lastFilter = mainPanel.getFilter();
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public void saveAsPDB(MolecularAssembly[] assemblies, File file) {
    mainPanel.saveAsPDB(assemblies, file);
    lastFilter = mainPanel.getFilter();
  }

  @Override
  public void saveAsXYZ(MolecularAssembly assembly, File file) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    mainPanel.saveAsXYZ(file);
    lastFilter = mainPanel.getFilter();
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public void savePDBSymMates(MolecularAssembly assembly, File file) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    mainPanel.savePDBSymMates(file, "_symMate");
    lastFilter = mainPanel.getFilter();
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public void savePDBSymMates(MolecularAssembly assembly, File file, String suffix) {
    Optional<FFXSystem> origSys = switchTo(assembly);
    mainPanel.savePDBSymMates(file, suffix);
    lastFilter = mainPanel.getFilter();
    origSys.ifPresent(this::switchBack);
  }

  @Override
  public double time() {
    return modelingShell.time();
  }

  /**
   * Switches the hierarchy's active system to assembly if assembly is present inside the hierarchy;
   * returns an Optional FFXSystem of the prior active system, or an empty Optional if no switch was
   * made
   *
   * @param assembly To switch to
   * @return Original system if switched, else empty Optional
   */
  private Optional<FFXSystem> switchTo(MolecularAssembly assembly) {
    Optional<FFXSystem> origSystem;
    if (!(assembly instanceof FFXSystem)) {
      origSystem = Optional.empty();
      return origSystem;
    }

    Hierarchy hierarchy = mainPanel.getHierarchy();
    FFXSystem activeSys = hierarchy.getActive();
    FFXSystem assemblySys = (FFXSystem) assembly;

    for (FFXSystem sys : hierarchy.getSystems()) {
      if (sys == assemblySys) {
        origSystem = Optional.of(activeSys);
        hierarchy.setActive(assemblySys);
        return origSystem;
      }
    }

    origSystem = Optional.empty();
    return origSystem;
  }

  /**
   * Switches the hierarchy's active system back to what it was.
   *
   * @param origSystem The active system.
   */
  private void switchBack(FFXSystem origSystem) {
    if (origSystem != null) {
      Hierarchy hierarchy = mainPanel.getHierarchy();
      hierarchy.setActive(origSystem);
    }
  }
}
