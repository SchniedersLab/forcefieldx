package ffx.potential

// SAVE AS XYZ

// Apache Imports
import org.apache.commons.io.FilenameUtils
import org.apache.commons.configuration.CompositeConfiguration

// Groovy Imports
import groovy.cli.Option
import groovy.util.CliBuilder
import groovy.cli.Unparsed

// FFX Imports
import ffx.crystal.Crystal
import ffx.crystal.SpaceGroup
import ffx.crystal.SymOp
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.parameters.ForceField
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils

/**
 * The SaveAsXYZ script saves a file as an XYZ file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsXYZ [options] &lt;filename&gt;
 */
class PrepareSpaceGroups extends Script {
    /**
     * Options for the SaveAsXYZ script
     * <br>
     * Usage:
     * <br>
     * ffxc SaveAsXYZ [options] &lt;filename&gt;
     */
    public class Options{
        /**
         * -h or --help to print a help message
         */
        @Option(longName='help', shortName='h', defaultValue='false', description='Print this help message.') boolean help
        /**
         * -r or --PDB only consider space groups ranked above the supplied cut-off in the Protein Databank.
         */
        @Option(longName='PDB', shortName='r', defaultValue='231',
            description='Only consider space groups populated above the specified percentage in the Protein Databank.') int rank
        /**
         * -p or --CSD only consider space groups populated above the specified percentage in the CSD.
         */
        @Option(longName='CSD', shortName='p', defaultValue='0.0',
            description='Only consider space groups populated above the specified percentage in the CSD.') double percent
        /**
         * -c or --chiral to create directories only for chiral space groups.
         */
        @Option(longName='chiral', shortName='c', defaultValue='false',
            description='Only consider chiral space groups.') boolean chiral
        /**
         * -sym or --symOp random Cartesian symmetry operator will use the specified translation range -Arg .. Arg (default = 1.0 A).
         */
        @Option(shortName='rsym', longName='randomSymOp', defaultValue='1.0',
            description='Random Cartesian sym op will choose a translation in the range -Arg .. Arg (default = 1.0A).') double symScalar
        /**
         * -d or --density random unit cell axes will be used achieve the specified density (g/cc).
         */
        @Option(shortName='d', longName='density', defaultValue='1.0',
            description='Random unit cell axes will be used to achieve the specified density (default = 1.0 g/cc).') double density
        /**
         * -sg or --spacegroup prepare a directory for a single spacegroup (no default).
         */
        @Option(shortName='sg', longName='spacegroup', description='Prepare a directory for a single spacegroup.') String sg
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }


    /**
     * Execute the script.
     */
    def run() {
        int offset = 0;

        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc PrepareSpaceGroups [options] <filename>');
        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            return cli.usage()
        }

        List<String> arguments = options.filenames;

        String modelFilename = null
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            modelFilename = arguments.get(0)
            //open(modelFilename)
        } else if (active == null) {
            return cli.usage()
        } else {
            modelFilename = active.getFile()
        }

        logger.info("\n Preparing space group directories for " + modelFilename);

        System.setProperty("ewald-alpha", "0.0");

        PotentialsFunctions functions
        try {
            // Use a method closure to try to get an instance of UIUtils (the User Interfaces
            // implementation, which interfaces with the GUI, etc.).
            functions = getPotentialsUtils()
        } catch (MissingMethodException ex) {
            // If Groovy can't find the appropriate closure, catch the exception and build
            // an instance of the local implementation.
            functions = new PotentialsUtils()
        }

        MolecularAssembly[] assemblies = functions.open(modelFilename)
        MolecularAssembly activeAssembly = assemblies[0]
        ForceFieldEnergy energy = activeAssembly.getPotentialEnergy()
        CompositeConfiguration config = activeAssembly.getProperties()

        File coordFile = activeAssembly.getFile()
        String coordName = FilenameUtils.getName(coordFile.getPath())
        File propertyFile = new File(config.getProperty("propertyFile"))

        Atom[] atoms = activeAssembly.getAtomArray()
        double mass = activeAssembly.getMass()
        double density = options.density
        if (density < 0.0) {
            density = 1.0
        }

        for (int num = 1; num <=230; num++) {
            SpaceGroup spacegroup = SpaceGroup.spaceGroupFactory(num)
            if (options.sg) {
                SpaceGroup spacegroup2 = SpaceGroup.spaceGroupFactory(options.sg);
                if (spacegroup2 == null) {
                    logger.info(String.format("\n Space group %s was not recognized.\n", sg));
                    return;
                }
                if (spacegroup2.number != spacegroup.number) {
                    continue;
                }
            }

            if (options.chiral && !spacegroup.isChiral()) {
                continue
            }

            if (spacegroup.csdPercent[num - 1] < options.percent) {
                continue
            }

            if (SpaceGroup.getPDBRank(spacegroup) > options.rank) {
                continue
            }

            logger.info(String.format("\n Preparing %s (CSD percent: %7.4f, PDB Rank: %d)",
                    spacegroup.shortName, spacegroup.csdPercent[num - 1], SpaceGroup.getPDBRank(spacegroup)));

            // Create the directory.
            String sgDirName = spacegroup.shortName.replace('/', '_')
            File sgDir = new File(FilenameUtils.getFullPath(coordFile.getAbsolutePath()) + sgDirName)
            if (!sgDir.exists()) {
                logger.info("\n Creating space group directory: " + sgDir.toString())
                sgDir.mkdir()
            }

            double[] abc = spacegroup.randomUnitCellParams()
            Crystal crystal = new Crystal(abc[0], abc[1], abc[2], abc[3], abc[4], abc[5], spacegroup.shortName)
            crystal.setDensity(density, mass)
            energy.setCrystal(crystal);

            if (options.symScalar > 0.0) {
                SymOp symOp = SymOp.randomSymOpFactory(options.symScalar)
                logger.info(String.format("\n Applying random Cartesian SymOp:\n %s", symOp.toString()))
                double[] xyz = new double[3]
                for (int i=0; i<atoms.length; i++) {
                    atoms[i].getXYZ(xyz)
                    crystal.applyCartesianSymOp(xyz, xyz, symOp)
                    atoms[i].setXYZ(xyz)
                }
            }

            // Save the coordinate file.
            File sgFile = new File(sgDir.getAbsolutePath() + File.separator + coordName)
            logger.info(" Saving " + sgDirName + " coordinates to: " + sgFile.toString())
            functions.save(activeAssembly, sgFile)

            File keyFile = new File(FilenameUtils.removeExtension(sgFile.getAbsolutePath()) + ".properties")
            logger.info(" Saving " + sgDirName + " properties to: " + keyFile.toString())

            PrintWriter out = null;
            BufferedReader keyReader = null;
            try {
                out = new PrintWriter(new BufferedWriter(new FileWriter(keyFile)));
                keyReader = new BufferedReader(new FileReader(propertyFile));
                while (keyReader.ready()) {
                    String line = keyReader.readLine().trim();
                    if (line != null && !line.equalsIgnoreCase("")) {
                        String[] tokens = line.split(" +")
                        if (tokens != null && tokens.length > 1) {
                            if (tokens[0].equalsIgnoreCase("a-axis")) {
                                out.println(String.format("a-axis %12.8f", crystal.a));
                            } else if (tokens[0].equalsIgnoreCase("b-axis")) {
                                out.println(String.format("b-axis %12.8f", crystal.b));
                            } else if (tokens[0].equalsIgnoreCase("c-axis")) {
                                out.println(String.format("c-axis %12.8f", crystal.c));
                            } else if (tokens[0].equalsIgnoreCase("alpha")) {
                                out.println(String.format("alpha %12.8f", crystal.alpha));
                            } else if (tokens[0].equalsIgnoreCase("beta")) {
                                out.println(String.format("beta %12.8f", crystal.beta));
                            } else if (tokens[0].equalsIgnoreCase("gamma")) {
                                out.println(String.format("gamma %12.8f", crystal.gamma));
                            } else if (tokens[0].equalsIgnoreCase("spacegroup")) {
                                out.println(String.format("spacegroup %s", spacegroup.shortName));
                            } else {
                                out.println(line);
                            }
                        } else {
                            out.println(line);
                        }
                    }
                }
            } catch (IOException ex) {
                logger.warning(" Exception writing keyfile." + ex.toString());
            } finally {
                if (out != null) {
                    out.flush();
                    out.close();
                }
                if (keyReader != null) {
                    keyReader.close();
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
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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