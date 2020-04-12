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
package ffx.potential.utils;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import static java.lang.String.format;

import org.apache.commons.io.FilenameUtils;

import ffx.crystal.Crystal;
import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.PDBFilter;
import ffx.potential.parsers.PDBFilter.Mutation;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.parsers.XYZFilter;

/**
 * <p>
 * PotentialsUtils implements core functionality for many Force Field X algorithms and
 * scripts, such as opening and closing structure files, basic force field evaluations,
 * etc. This implementation does not do anything on top of what is specified by the
 * interface, and is used primarily by tests. It is also potentially useful for third
 * parties who would like to use FFX without its graphical user interface.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 * @since 1.0
 */
public class PotentialsUtils implements PotentialsFunctions {

    private static final Logger logger = Logger.getLogger(PotentialsUtils.class.getName());
    private static final Logger potLog = Logger.getLogger("");
    private static Level levelBak = null;
    private final long initTime;
    private long interTime;
    private SystemFilter lastFilter;

    /**
     * <p>Constructor for PotentialsUtils.</p>
     */
    public PotentialsUtils() {
        initTime = System.nanoTime();
        interTime = initTime;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Shuts down parallel teams in the force field of the provided
     * MolecularAssembly. Kaminsky's ParallelTeamThreads' run() methods are
     * infinite loops, and because running threads are always GC roots, it is
     * necessary to send them a signal to shut down to enable garbage
     * collection.
     */
    @Override
    public void close(MolecularAssembly assembly) {
        assembly.destroy();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Shuts down parallel teams in the force fields of the provided
     * MolecularAssemblys.
     */
    @Override
    public void closeAll(MolecularAssembly[] assemblies) {
        for (MolecularAssembly assembly : assemblies) {
            assembly.destroy();
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Evaluates the energy of a MolecularAssembly and returns its
     * ForceFieldEnergy object.
     */
    @Override
    public ForceFieldEnergy energy(MolecularAssembly assembly) {
        if (assembly == null) {
            logger.info(" Molecular assembly was null - skipping energy");
            return null;
        } else {
            ForceFieldEnergy energy = assembly.getPotentialEnergy();
            if (energy == null) {
                energy = ForceFieldEnergy.energyFactory(assembly);
                assembly.setPotential(energy);
            }
            energy.energy(false, true);
            return energy;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public SystemFilter getFilter() {
        return lastFilter;
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns true (this is the local implementation).
     */
    @Override
    public boolean isLocal() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public MolecularAssembly open(String filename) {
        PotentialsFileOpener opener = new PotentialsFileOpener(filename);
        opener.run();
        lastFilter = opener.getFilter();
        if (opener.getAllAssemblies().length > 1) {
            logger.log(Level.WARNING, "Found multiple assemblies in file {0}, opening first.", filename);
        }
        return opener.getAssembly();
    }

    /**
     * One one file object.
     *
     * @param file a {@link java.io.File} object.
     * @return a {@link ffx.potential.MolecularAssembly} object.
     */
    public MolecularAssembly open(File file) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.run();
        lastFilter = opener.getFilter();
        if (opener.getAllAssemblies().length > 1) {
            logger.log(Level.WARNING, "Found multiple assemblies in file {0}, opening first.", file.getName());
        }
        return opener.getAssembly();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Opens an array of files and returns all created MolecularAssembly
     * objects, setting any underlying Potential to use a certain number of
     * threads.
     */
    @Override
    public MolecularAssembly[] open(String[] files, int nThreads) {
        PotentialsFileOpener opener = new PotentialsFileOpener(files);
        opener.setNThreads(nThreads);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Opens a file and returns all created MolecularAssembly objects.
     */
    @Override
    public MolecularAssembly[] openAll(String file) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Opens an array of files and returns the created MolecularAssembly
     * objects.
     */
    @Override
    public MolecularAssembly[] openAll(String[] files) {
        PotentialsFileOpener opener = new PotentialsFileOpener(files);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Opens a file and returns all created MolecularAssembly objects, setting
     * any underlying Potential to use a certain number of threads.
     */
    @Override
    public MolecularAssembly[] openAll(String file, int nThreads) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.setNThreads(nThreads);
        opener.run();
        lastFilter = opener.getFilter();
        return opener.getAllAssemblies();
    }

    /**
     * Open one filename string without printing all the header material.
     *
     * @param filename a {@link java.lang.String} object.
     * @return a {@link ffx.potential.MolecularAssembly} object.
     */
    public MolecularAssembly openQuietly(String filename) {
        setSilentPotential(true);
        MolecularAssembly mola = open(filename);
        setSilentPotential(false);
        return mola;
    }

    /**
     * Mutates file on-the-fly as it is being opened.
     * Used to open files for pHMD in fully-protonated form.
     *
     * @param file      a {@link java.io.File} object.
     * @param mutations a {@link java.util.List} object.
     * @return a {@link ffx.potential.MolecularAssembly} object.
     */
    public MolecularAssembly openWithMutations(File file, List<Mutation> mutations) {
        PotentialsFileOpener opener = new PotentialsFileOpener(file);
        opener.setMutations(mutations);
        opener.run();
        lastFilter = opener.getFilter();
        if (opener.getAllAssemblies().length > 1) {
            logger.log(Level.WARNING, "Found multiple assemblies in file {0}, opening first.", file.getName());
        }
        return opener.getAssembly();
    }

    /**
     * {@inheritDoc}
     * <p>
     * Returns the energy of a MolecularAssembly in kcal/mol (as a double) and
     * prints the energy evaluation
     */
    @Override
    public double returnEnergy(MolecularAssembly assembly) {
        if (assembly == null) {
            logger.info(" Molecular assembly was null - skipping energy");
            return 0.0;
        } else {
            ForceFieldEnergy energy = assembly.getPotentialEnergy();
            if (energy == null) {
                energy = ForceFieldEnergy.energyFactory(assembly);
                assembly.setPotential(energy);
            }
            return energy.energy(false, true);
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Saves the current state of a MolecularAssembly to an XYZ file.
     */
    @Override
    public void save(MolecularAssembly assembly, File file) {
        saveAsXYZ(assembly, file);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Saves the current state of a MolecularAssembly to an XYZ file as a P1
     * crystal.
     */
    @Override
    public void saveAsP1(MolecularAssembly assembly, File file) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            XYZFilter filter = new XYZFilter(file, assembly, null, null);
            ForceField forceField = assembly.getForceField();
            final double a = forceField.getDouble("A_AXIS", 10.0);
            final double b = forceField.getDouble("B_AXIS", a);
            final double c = forceField.getDouble("C_AXIS", a);
            final double alpha = forceField.getDouble("ALPHA", 90.0);
            final double beta = forceField.getDouble("BETA", 90.0);
            final double gamma = forceField.getDouble("GAMMA", 90.0);
            final String spacegroup = forceField.getString("SPACEGROUP", "P1");
            Crystal crystal = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
            if (!filter.writeFileAsP1(file, false, crystal)) {
                logger.info(format(" Save failed for %s", assembly.toString()));
            }
            lastFilter = filter;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Saves the current state of a MolecularAssembly to a PDB file.
     */
    @Override
    public void saveAsPDB(MolecularAssembly assembly, File file) {
        saveAsPDB(assembly, file, true, false);
    }

    /**
     * {@inheritDoc}
     * <p>
     * Saves the current state of a MolecularAssembly to a PDB file.
     */
    @Override
    public void saveAsPDB(MolecularAssembly assembly, File file, boolean writeEnd, boolean append) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file, assembly, null, null);
            if (!pdbFilter.writeFile(file, append, false, writeEnd)) {
                logger.info(format(" Save failed for %s", assembly.toString()));
            }
            lastFilter = pdbFilter;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Saves the current state of an array of MolecularAssemblys to a PDB file.
     */
    @Override
    public void saveAsPDB(MolecularAssembly[] assemblies, File file) {
        if (assemblies == null) {
            logger.info(" Null array of molecular assemblies to write.");
        } else if (assemblies.length == 0) {
            logger.info(" Zero-length array of molecular assemblies to write.");
        } else if (file == null) {
            logger.info(" No valid file to write to.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file,
                    Arrays.asList(assemblies), null, null);
            pdbFilter.writeFile(file, false);
            lastFilter = pdbFilter;
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Saves the current state of a MolecularAssembly to an XYZ file.
     */
    @Override
    public void saveAsXYZ(MolecularAssembly assembly, File file) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            XYZFilter xyzFilter = new XYZFilter(file, assembly, null, null);
            if (!xyzFilter.writeFile(file, false)) {
                logger.info(format(" Save failed for %s", assembly.toString()));
            }
            lastFilter = xyzFilter;
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void savePDBSymMates(MolecularAssembly assembly, File file) {
        savePDBSymMates(assembly, file, "_symMate");
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void savePDBSymMates(MolecularAssembly assembly, File file, String suffix) {
        if (assembly == null) {
            logger.info(" Assembly to save was null.");
        } else if (file == null) {
            logger.info(" No valid file provided to save assembly to.");
        } else {
            PDBFilter pdbFilter = new PDBFilter(file, assembly, null, null);
            lastFilter = pdbFilter;
            if (!pdbFilter.writeFile(file, false)) {
                logger.info(format(" Save failed for %s", assembly.toString()));
            } else {
                Crystal crystal = assembly.getCrystal();
                int nSymOps = crystal.spaceGroup.getNumberOfSymOps();
                String filename = FilenameUtils.removeExtension(file.getName());
                for (int i = 1; i < nSymOps; i++) {
                    pdbFilter.setSymOp(i);
                    String saveFileName = filename + suffix + "_" + i + ".pdb";
                    File saveFile = new File(saveFileName);
                    for (int j = 1; j < 1000; j++) {
                        if (!saveFile.exists()) {
                            break;
                        }
                        saveFile = new File(saveFileName + "_" + j);
                    }
                    StringBuilder symSb = new StringBuilder();
                    String[] symopLines = crystal.spaceGroup.getSymOp(i).toString().split("\\r?\\n");
                    int nLines = symopLines.length;
                    symSb.append("REMARK 350\nREMARK 350 SYMMETRY OPERATORS");
                    for (int j = 0; j < nLines; j++) {
                        symSb.append("\nREMARK 350 ").append(symopLines[j]);
                    }

                    symopLines = crystal.spaceGroup.getSymOp(i).toXYZString().split("\\r?\\n");
                    nLines = symopLines.length;
                    symSb.append("\nREMARK 350\nREMARK 350 SYMMETRY OPERATORS XYZ FORM");
                    for (int j = 0; j < nLines; j++) {
                        symSb.append("\nREMARK 350 ").append(symopLines[j]);
                    }

                    if (saveFile.exists()) {
                        logger.warning(format(" Could not successfully version file "
                                + "%s: appending to file %s", saveFileName, saveFile.getName()));
                        if (!pdbFilter.writeFileWithHeader(saveFile, symSb.toString(), true)) {
                            logger.info(format(" Save failed for %s", saveFile.getName()));
                        }
                    } else if (!pdbFilter.writeFileWithHeader(saveFile, symSb.toString(), false)) {
                        logger.info(format(" Save failed for %s", saveFile.getName()));
                    }
                }
            }
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Logs time since this interface was created and the last time this method
     * was called.
     */
    @Override
    public double time() {
        long currTime = System.nanoTime();
        logger.info(format(" Time since interface established: %f", (currTime - initTime) * 1.0E-9));
        double elapsed = (currTime - interTime) * 1.0E-9;
        interTime = currTime;
        logger.info(format(" Time since last timer call: %f", elapsed));
        return elapsed;
    }

    /**
     * <p>setSilentPotential.</p>
     *
     * @param silent a boolean.
     */
    void setSilentPotential(boolean silent) {
        if (silent && potLog.isLoggable(Level.INFO) && levelBak == null) {
            levelBak = potLog.getLevel();
            potLog.setLevel(Level.WARNING);
        }
        if (!silent && levelBak != null) {
            potLog.setLevel(levelBak);
            levelBak = null;
        }
    }

}
