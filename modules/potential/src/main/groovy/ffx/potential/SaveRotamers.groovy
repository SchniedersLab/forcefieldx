package ffx.potentials

import org.apache.commons.io.FilenameUtils

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Polymer
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.cli.PotentialScript

import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * The SaveRotamers script saves a file as a PDB file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveRotamers [options] &lt;filename&gt;
 */
@Command(description = " Save out rotamers.", name = "ffxc SaveRotamers")
class SaveRotamers extends PotentialScript {

    /**
     * -c or --chain to specify chain
     */
    @Option(names = ['--chain', '-c'], paramLabel = "\" \"",
            description = 'Single character chain name.')
    char c = ' '

    /**
     * -l or --library to select rotamer library 
     */
    @Option(names = ['--library', '-l'], paramLabel = "1",
            description = 'Available rotamer libraries are (1) Ponder and Richards or (2) Richardson.')
    int library = 1

    /**
     * -r or --resid to select residue number 
     */
    @Option(names = ['--resid', '-r'], paramLabel = "1",
            description = 'Residue number.')
    int resID = 1

    /**
     * -i or --independent to draw nucleic acid rotamers independently of chain context
     */
    @Option(names = ['--independent', '-i'], paramLabel = 'false',
            description = 'Independent draws nucleic acid rotamers independently of chain context.')
    boolean independent = false

    /**
     * -s or --start to select first rotamer to draw. Indexed form rotamer 0
     */
    @Option(names = ['--start', '-s'], paramLabel = "0",
            description = 'First rotamer to draw (indexed from rotamer 0).')
    int start = 0

    /**
     * -f or --finish to select last rotamer to draw. Indexed from rotamer 0
     */
    @Option(names = ['--finish', '-f'], paramLabel = "-1",
            description = 'Last rotamer to draw (indexed from rotamer 0).')
    int finish = -1

    /**
     * -x or --all to draw all rotamers beginning from the passed rotamer number (overrides other options). Indexed from rotamer 0.
     */
    @Option(names = ['--all', '-x'], paramLabel = "-1",
            description = 'Draw all rotamers beginning from the passed index (overrides other options).')
    int all = -1

    /**
     * -u or --upstreamPucker to adjust the pucker of the 5\' residue to match the rotamer. Use flag to turn on
     */
    @Option(names = ['--upstreamPucker', '-u'], paramLabel = 'false',
            description = 'Adjusts the pucker of the 5\' residue to match the rotamer.')
    boolean upstreamPucker = false

    /**
     * One or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = "XYZ or PDB input file.")
    private List<String> filenames

    /**
     * Execute the script.
     */
    def run() {

        if (!init()) {
            return
        }

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            modelFilename = filenames.get(0)
            MolecularAssembly[] assemblies = potentialFunctions.open(modelFilename)
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        } else {
            modelFilename = activeAssembly.getFile().getAbsolutePath()
        }

        RotamerLibrary rLib = RotamerLibrary.getDefaultLibrary()
        String chain = Character.toString(c)

        boolean saveAllRotamers = false
        int allStart = 0
        // Start from to all
        if (all > -1) {
            allStart = all
            saveAllRotamers = true
        }

        logger.info("\n Saving rotamers for residue number " + resID + " of chain " + chain + ".")

        RotamerLibrary.initializeDefaultAtomicCoordinates(activeAssembly.getChains())
        Polymer polymer = activeAssembly.getChain(chain)
        if (polymer == null) {
            logger.info(" Polymer + " + chain + " does not exist.")
            return
        }
        Residue residue = polymer.getResidue(resID)
        if (residue == null) {
            logger.info(" Residue + " + resID + " does not exist.")
            return
        }

        if (library == 1) {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.PonderAndRichards)
        } else {
            rLib.setLibrary(RotamerLibrary.ProteinLibrary.Richardson)
        }

        Rotamer[] rotamers = residue.getRotamers(rLib)
        if (rotamers == null) {
            logger.severe(" There are no rotamers for residue + " + residue.toString())
        }
        int nrotamers = rotamers.length

        boolean isDeoxy = false // Applies to prevResidue.
        Residue prevResidue
        // If upstreamPucker is specified true, ensure that it is valid to apply it
        // (residue is nucleic acid with an upstream partner), and set isDeoxy.
        // If it is invalid, set upstreamPucker false.
        if (upstreamPucker) {
            // Exception gets thrown if it's an amino acid, since "NA" is undefined.
            try {
                if (residue.getResidueType() == NA) {
                    prevResidue = (Residue) residue.getPreviousResidue()
                    // If no previous residue, set upstream pucker false.
                    // The method used will ensure prevResidue is a nucleic acid.
                    if (prevResidue == null) {
                        upstreamPucker = false
                    } else {
                        Atom HOs = (Atom) prevResidue.getAtomNode("HO\'")
                        if (HOs == null) {
                            isDeoxy = true
                        }
                    }
                } else {
                    upstreamPucker = false
                }
            } catch (Exception e) {
                upstreamPucker = false
            }
        }

        String ext = FilenameUtils.getExtension(modelFilename)
        modelFilename = FilenameUtils.removeExtension(modelFilename)

        if (saveAllRotamers) {
            if (allStart >= nrotamers) {
                logger.info(" Specified start range is outside of rotamer range. No action taken.")
            } else {
                for (int i = allStart; i < nrotamers; i++) {
                    RotamerLibrary.applyRotamer(residue, rotamers[i], independent);
                    if (upstreamPucker) {
                        double prevDelta = rotamers[i].chi1;
                        if (RotamerLibrary.checkPucker(prevDelta) == 1) {
                            // North pucker
                            RotamerLibrary.applySugarPucker(prevResidue, 1, isDeoxy, true);
                        } else {
                            RotamerLibrary.applySugarPucker(prevResidue, 2, isDeoxy, true);
                        }
                    }
                    if (ext.toUpperCase().contains("XYZ")) {
                        logger.info(String.format("Saving rotamer %d", i))
                        potentialFunctions.saveAsXYZ(activeAssembly, new File(modelFilename + ".xyz"));
                    } else {
                        logger.info(String.format("Saving rotamer %d", i))
                        potentialFunctions.saveAsPDB(activeAssembly, new File(modelFilename + ".pdb"))
                    }
                }
            }
        } else {
            if (start >= nrotamers) {
                logger.info(" Specified start range is outside of rotamer range. No action taken.")
            } else {
                if (finish >= nrotamers) {
                    finish = nrotamers - 1
                } else if (finish < start) {
                    logger.info(" Specified finish point is before the start point drawing only rotamer " + start)
                    finish = start
                }
                for (int i = start; i <= finish; i++) {
                    RotamerLibrary.applyRotamer(residue, rotamers[i], independent)
                    if (upstreamPucker) {
                        double prevDelta = rotamers[i].chi1
                        if (RotamerLibrary.checkPucker(prevDelta) == 1) {
                            // North pucker
                            RotamerLibrary.applySugarPucker(prevResidue, 1, isDeoxy, true)
                        } else {
                            RotamerLibrary.applySugarPucker(prevResidue, 2, isDeoxy, true)
                        }
                    }
                    if (ext.toUpperCase().contains("XYZ")) {
                        logger.info(String.format("Saving rotamer %d", i))
                        potentialFunctions.saveAsXYZ(activeAssembly, new File(filename + ".xyz"))
                    } else {
                        logger.info(String.format("Saving rotamer %d", i))
                        potentialFunctions.saveAsPDB(activeAssembly, new File(modelFilename + ".pdb"))
                    }
                }
            }
        }
    }
}

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
