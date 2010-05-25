/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.ui;

import static java.lang.String.format;

import java.awt.Cursor;

import java.util.List;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;

import ffx.potential.bonded.Utilities;
import ffx.potential.PotentialEnergy;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.SystemFilter;

/**
 * The FileOpener class opens a file into Force Field X using a filter
 * from the ffe.parsers package. The OpenFile class implements the Runnable
 * interface so that opening a file does not freeze FFX.
 */
public class FileOpener
        implements Runnable {

    private static final Logger logger = Logger.getLogger(FileOpener.class.getName());
    private static final long KB = 1024;
    private static final long MB = KB * KB;
    SystemFilter systemFilter = null;
    MainPanel mainPanel = null;
    private boolean timer = false;
    private boolean gc = false;
    private long occupiedMemory;
    private long time;

    public FileOpener(SystemFilter systemFilter, MainPanel mainPanel) {
        this.systemFilter = systemFilter;
        this.mainPanel = mainPanel;
        if (System.getProperty("ffx.timer", "false").equalsIgnoreCase("true")) {
            timer = true;
            if (System.getProperty("ffx.timer.gc", "false").equalsIgnoreCase(
                    "true")) {
                gc = true;
            }
        }
    }

    private void open() {
        if (timer) {
            startTimer();
        }
        FFXSystem ffxSystem = null;
        // Continue if the file was read in successfully.
        if (systemFilter.readFile()) {
            ffxSystem = (FFXSystem) systemFilter.getActiveMolecularSystem();
            if (!(systemFilter instanceof PDBFilter)) {
                Utilities.biochemistry(ffxSystem, systemFilter.getAtomList());
            }
            // Add the system to the multiscale hierarchy.
            mainPanel.getHierarchy().addSystemNode(ffxSystem);
            PotentialEnergy energy = new PotentialEnergy(ffxSystem);
            ffxSystem.setPotential(energy);

            // Check if there are alternate conformers
            if (systemFilter instanceof PDBFilter) {
                PDBFilter pdbFilter = (PDBFilter) systemFilter;
                List<Character> altLocs = pdbFilter.getAltLocs();
                if (altLocs.size() > 1 || altLocs.get(0) != ' ') {
                    StringBuilder altLocString = new StringBuilder(" Alternate locations [ ");
                    for (Character c : altLocs) {
                        // Do not report the root conformer.
                        if (c == ' ') {
                            continue;
                        }
                        altLocString.append(format("(%s) ", c));
                    }
                    altLocString.append("]");
                    logger.info(altLocString.toString());
                } else {
                    logger.info(" No alternate conformers detected.");
                }

                /**
                 * Alternate conformers may have different chemistry, so they
                 * each need to be their own FFX system.
                 */
                for (Character c : altLocs) {
                    if (c.equals(' ') || c.equals('A')) {
                        continue;
                    }
                    FFXSystem newSystem = new FFXSystem(ffxSystem.getFile(),
                                                        "Alternate Location " + c, ffxSystem.getProperties());
                    newSystem.setForceField(ffxSystem.getForceField());
                    pdbFilter.setAltID(newSystem, c);
                    pdbFilter.clearSegIDs();
                    if (pdbFilter.readFile()) {
                        String fileName = ffxSystem.getFile().getAbsolutePath();
                        newSystem.setName(FilenameUtils.getBaseName(fileName) + " " + c);
                        mainPanel.getHierarchy().addSystemNode(newSystem);
                        energy = new PotentialEnergy(newSystem);
                        newSystem.setPotential(energy);
                    }
                }
            }
        }
        mainPanel.setCursor(Cursor.getDefaultCursor());
        if (timer) {
            stopTimer(ffxSystem);
        }
    }

    @Override
    public void run() {
        if (mainPanel != null && systemFilter != null) {
            open();
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
        logger.info("Opened " + ffxSystem.toString() + " with " + ffxSystem.getAtomList().size() + " atoms.\n" + "File Op Time  (msec): " + time * 1.0e-9);
        Runtime runtime = Runtime.getRuntime();
        if (gc) {
            runtime.runFinalization();
            runtime.gc();
            long moleculeMemory = (runtime.totalMemory() - runtime.freeMemory()) - occupiedMemory;
            logger.info("File Op Memory  (Kb): " + moleculeMemory / KB);
        }
        occupiedMemory = runtime.totalMemory() - runtime.freeMemory();
    }
}
