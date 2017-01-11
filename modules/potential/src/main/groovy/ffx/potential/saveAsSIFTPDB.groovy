package ffx.potential


// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Groovy Imports
import groovy.cli.Option
import groovy.util.CliBuilder;
import groovy.cli.Unparsed

//FFX Imports
import ffx.potential.MolecularAssembly
import ffx.potential.utils.PotentialsFunctions
import ffx.potential.utils.PotentialsUtils;


/**
 * The SaveAsSIFTPDB script saves a file as a PDB file or as a SIFTPDB file
 * <br>
 * Usage:
 * <br>
 * ffxc SaveAsSIFTPDB [options] &lt;filename&gt;
 */
class SaveAsSIFTPDB extends Script {
    /**
     * The SaveAsSIFTPDB script saves a file as a PDB file or as a SIFTPDB file
     * <br>
     * Usage:
     * <br>
     * ffxc SaveAsSIFTPDB [options] &lt;filename&gt;
     */
    public class Options{
        /**
         * -h or --help to print a help message
         */
        @Option(longName='help', shortName='h', defaultValue='false', description='Print this help message.') boolean help
        /**
         * -f or --filename to input file of sift scores to enter
         */
        @Option(longName='fileName', shortName='f', description='File of sift scores to enter.') boolean filename
        /**
         * The final argument(s) should be one or more filenames.
         */
        @Unparsed List<String> filenames
    }
    
    /**
     * Execute the script.
     */
    def run() {
        String siftFileName = null;


        // Create the command line parser.
        def cli = new CliBuilder(usage:' ffxc SaveAsSIFTPDB [options] <filename>');
        def options = new Options()
        //cli.f(longOpt: 'fileName', args:1, argName:'GENENAME.out', 'File of sift scores to enter.');
        cli.parseFromInstance(options, args)
        
        if (options.help == true) {
            return cli.usage()
        }

        if (options.filename == true) {
            siftFileName = options.filename;
        }
        
        List<String> arguments = options.filenames;
        

        String[] data;

        File siftFile = new File(siftFileName);

        if (siftFile.exists()) {
            data = siftFile.text.split('\n')
        }

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

        logger.info("\n Writing out PDB for " + modelFilename);
        
        PotentialsFunctions functions
        try {
            // Use a method closure to try to get an instance of UIUtils (the User Interfaces
            // implementation, which interfaces with the GUI, etc.).
            functions = getPotentialsFunctions()
        } catch (MissingMethodException ex) {
            // If Groovy can't find the appropriate closure, catch the exception and build
            // an instance of the local implementation.
            functions = new PotentialsUtils()
        }

        MolecularAssembly[] assemblies = functions.open(modelFilename)
        MolecularAssembly activeAssembly = assemblies[0]
        modelFilename = FilenameUtils.removeExtension(modelFilename) + ".pdb";
        
        if (siftFileName == null) {
            functions.saveAsPDB(activeAssembly, new File(modelFilename));
        } else {
            functions.saveAsSIFTPDB(activeAssembly, new File(modelFilename), data)
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