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
package ffx.ui;

import java.awt.Cursor;
import java.io.File;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.Utilities;
import ffx.potential.parsers.BioJavaFilter;
import ffx.potential.parsers.ConversionFilter;
import ffx.potential.parsers.FileOpener;

/**
 * The UIDataConverter class converts a data structure into a Force Field X data
 * structure using a filter from the ffx.potentials.parsers package. To avoid
 * freezing the FFX GUI, it implements the FileOpener interface, which extends Runnable.
 * <p>
 * Still need to finalize everything save the constructor.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class UIDataConverter implements FileOpener {

    private static final Logger logger = Logger.getLogger(UIDataConverter.class.getName());
    private static final long KB = 1024;
    private ConversionFilter conversionFilter;
    MainPanel mainPanel;
    private boolean timer = false;
    private boolean gc = false;
    private long occupiedMemory;
    private long time;
    private final Object dataStructure;

    /**
     * Constructs an object to convert a data structure to FFX MolecularAssembly.
     * <p>
     * If the data structure type is not recognized, an exception will be thrown when run.
     *
     * @param data             Data structure to load from.
     * @param file             Source file.
     * @param conversionFilter The conversionFilter.
     * @param mainPanel        The MainPanel instance.
     */
    public UIDataConverter(Object data, File file, ConversionFilter conversionFilter, MainPanel mainPanel) {
        if (conversionFilter instanceof BioJavaFilter) {
            this.dataStructure = data;
        } else {
            this.dataStructure = null;
        }
        this.conversionFilter = conversionFilter;
        this.mainPanel = mainPanel;
        if (System.getProperty("ffx.timer", "false").equalsIgnoreCase("true")) {
            timer = true;
            if (System.getProperty("ffx.timer.gc", "false").equalsIgnoreCase(
                    "true")) {
                gc = true;
            }
        }
    }

    /**
     * Converts the data structure to MolecularAssembly(s).
     */
    public void convert() {
        if (timer) {
            startTimer();
        }
        FFXSystem ffxSystem = null;
        // Continue if the file was read in successfully.
        if (conversionFilter.convert()) {
            ffxSystem = (FFXSystem) conversionFilter.getActiveMolecularSystem();
            if (!(conversionFilter instanceof BioJavaFilter)) {
                Utilities.biochemistry(ffxSystem, conversionFilter.getAtomList());
            }
            conversionFilter.applyAtomProperties();
            // Add the system to the multiscale hierarchy.
            mainPanel.getHierarchy().addSystemNode(ffxSystem);
            ForceFieldEnergy energy = ForceFieldEnergy.energyFactory(ffxSystem, conversionFilter.getCoordRestraints());
            ffxSystem.setPotential(energy);
            mainPanel.getHierarchy().setActive(ffxSystem);

            // Check if there are alternate conformers
            if (conversionFilter instanceof BioJavaFilter) {
                BioJavaFilter biojFilter = (BioJavaFilter) conversionFilter;
                List<Character> altLocs = biojFilter.getAltLocs();
                if (altLocs.size() > 1 || altLocs.get(0) != ' ') {
                    StringBuilder altLocString = new StringBuilder("\n Alternate locations found [ ");
                    for (Character c : altLocs) {
                        // Do not report the root conformer.
                        if (c == ' ') {
                            continue;
                        }
                        altLocString.append(String.format("(%s) ", c));
                    }
                    altLocString.append("]\n");
                    logger.info(altLocString.toString());
                }

                // Alternate conformers may have different chemistry, so they each need to be their own FFX system.
                for (Character c : altLocs) {
                    if (c.equals(' ') || c.equals('A')) {
                        continue;
                    }
                    FFXSystem newSystem = new FFXSystem(ffxSystem.getFile(),
                            "Alternate Location " + c, ffxSystem.getProperties());
                    newSystem.setForceField(ffxSystem.getForceField());
                    biojFilter.setAltID(newSystem, c);
                    biojFilter.clearSegIDs();
                    if (biojFilter.convert()) {
                        biojFilter.applyAtomProperties();
                        String fileName = ffxSystem.getFile().getAbsolutePath();
                        newSystem.setName(FilenameUtils.getBaseName(fileName) + " " + c);
                        mainPanel.getHierarchy().addSystemNode(newSystem);
                        energy = ForceFieldEnergy.energyFactory(newSystem, biojFilter.getCoordRestraints());
                        newSystem.setPotential(energy);
                    }
                }
            }
        } else {
            logger.warning(String.format(" Failed to convert structure %s", dataStructure.toString()));
        }
        mainPanel.setCursor(Cursor.getDefaultCursor());
        if (timer) {
            stopTimer(ffxSystem);
        }
    }

    /**
     * Returns the active MolecularAssembly from the user interface hierarchy.
     *
     * @return Active MolecularAssembly
     * @throws NullPointerException If no active MolecularAssembly
     */
    @Override
    public MolecularAssembly getAssembly() throws NullPointerException {
        MolecularAssembly assembly = mainPanel.getHierarchy().getActive();
        if (assembly == null) {
            throw new NullPointerException(" FFX hierarchy does not have an active assembly.");
        }
        return assembly;
    }

    /**
     * Returns all MolecularAssemblys in the user interface hierarchy.
     *
     * @return All MolecularAssembly objects stored by the hierarchy.
     * @throws NullPointerException If hierarchy has a null or empty list of
     *                              assemblies.
     */
    @Override
    public MolecularAssembly[] getAllAssemblies() throws NullPointerException {
        MolecularAssembly[] assemblies = mainPanel.getHierarchy().getSystems();
        if (assemblies == null) {
            throw new NullPointerException(" FFX hierarchy has a null list of assemblies.");
        } else if (assemblies.length == 0) {
            throw new NullPointerException(" FFX hierarchy has an empty list of assemblies.");
        } else {
            return assemblies;
        }
    }

    /**
     * Returns the properties of the hierarchy's active FFXSystem.
     *
     * @return Active properties
     */
    @Override
    public CompositeConfiguration getProperties() {
        return mainPanel.getHierarchy().getActive().getProperties();
    }

    /**
     * Returns the properties of all FFXSystems in the hierarchy.
     *
     * @return Properties for all systems.
     */
    @Override
    public CompositeConfiguration[] getAllProperties() {
        FFXSystem[] allSystems = mainPanel.getHierarchy().getSystems();
        int numSystems = allSystems.length;
        CompositeConfiguration[] allProperties = new CompositeConfiguration[numSystems];
        for (int i = 0; i < numSystems; i++) {
            allProperties[i] = allSystems[i].getProperties();
        }
        return allProperties;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        if (mainPanel != null && conversionFilter != null) {
            convert();
        }
    }

    /**
     * Rather verbose output for timed File Operations makes it easy to grep log
     * files for specific information.
     */
    private void startTimer() {
        Runtime runtime = Runtime.getRuntime();
        if (gc) {
            runtime.runFinalization();
            runtime.gc();
        }
        occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
        time -= System.nanoTime();
    }

    private void stopTimer(FFXSystem ffxSystem) {
        time += System.nanoTime();
        logger.log(Level.INFO, " Opened {0} with {1} atoms.\n File Op Time  (msec): {2}",
                new Object[]{ffxSystem.toString(), ffxSystem.getAtomList().size(), time * 1.0e-9});
        Runtime runtime = Runtime.getRuntime();
        if (gc) {
            runtime.runFinalization();
            runtime.gc();
            long moleculeMemory = (runtime.totalMemory() - runtime.freeMemory()) - occupiedMemory;
            logger.log(Level.INFO, " System Memory  (Kb): {0}", moleculeMemory / KB);
        }
        occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
        if (gc) {
            logger.log(Level.INFO, " Memory In Use  (Kb): {0}", occupiedMemory / KB);
        }
    }
}
