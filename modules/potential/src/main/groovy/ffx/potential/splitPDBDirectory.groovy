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

package ffx.utilities

// Apache Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

// Force Field X Imports
import ffx.potential.parsers.CoordinateFileFilter;

// BioJava Imports

// PJ imports


File directory;
String newDirectoriesName = "ffx_splitDirectory_";
int numPerDirectory;
boolean parseDeep = true;
List<File> files;
boolean verbose = false;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc splitPDBDirectory [options] <directory>');
cli.h(longOpt:'help', 'Print this message.');
cli.n(longOpt:'numberPerDirectory', args:1, argName:'20', 'Split the directory into sub-directories of this many PDB files.');
cli.p(longOpt:'parseDeep', args:1, argName:'true', 'Split files with valid coordinate file format (true) or look only at extension (false).');
cli.d(longOpt:'directoryNames', args:1, argName:'ffx_splitDirectory_', 'New sub-directories begin with this prefix.');
cli.v(longOpt:'verbose', 'Flag to print progress messages.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

if (options.n) {
    numPerDirectory = Integer.parseInt(options.n);
}
if (options.p) {
    parseDeep = Boolean.parseBoolean(options.p);
}
if (options.d) {
    newDirectoriesName = options.d;
}
if (options.v) {
    verbose = true;
}

String directoryName = arguments.get(0);
//newDirectoriesName = directoryName.concat(FilenameUtils.EXTENSION_SEPARATOR_STR.concat(newDirectoriesName));
newDirectoriesName = FilenameUtils.separatorsToSystem(new StringBuilder(directoryName).append('/').append(newDirectoriesName).toString());
directory = new File(directoryName);
files = new ArrayList<>();
if (!directory.exists()) {
    return cli.usage();
}
if (directory.isFile()) {
    BufferedReader bw = new BufferedReader(new FileReader(directory));
    try {
        CoordinateFileFilter filter = new CoordinateFileFilter();
        String line = bw.readLine();
        while (line != null) {
            line = line.trim();
            File file = new File(line);
            if (file.exists() && file.isFile()) {
                if (parseDeep) {
                    if (filter.acceptDeep(file)) {
                        files.add(file);
                    }
                } else {
                    if (filter.accept(file)) {
                        files.add(file);
                    }
                }
            }
            line = bw.readLine();
        }
    } catch (IOException ex) {
        bw.close();
        logger.severe(String.format(" Error in reading file list: %s", ex.toString()));
    }
    bw.close();
} else {
    try {
        File[] fileArray = directory.listFiles();
        CoordinateFileFilter filter = new CoordinateFileFilter();
        if (parseDeep) {
            for (File file : fileArray) {
                if (filter.acceptDeep(file)) {
                    files.add(file);
                }
            }
        } else {
            for (File file : fileArray) {
                if (filter.accept(file)) {
                    files.add(file);
                }
            }
        }
    } catch (IOException ex) {
        logger.severe(String.format(" Error in reading files in directory: %s", ex.toString()));
    }
}

int numFiles = files.size();
if (verbose) {
    int numDirectories = ((numFiles % numPerDirectory) == 0 ? (numFiles / numPerDirectory) : ((numFiles / numPerDirectory) + 1));
    logger.info(String.format(" %d files found: %d directories of %d files will be made.",
            numFiles, numDirectories, numPerDirectory));
}

for (int i = 0; i < numFiles; i += numPerDirectory) {
    int endRange = i + numPerDirectory; // Exclusive when indexed for list, inclusive when indexed for user..
    endRange = (endRange < numFiles ? endRange : numFiles);
    int start = i + 1; // Indexed for user.
    if (verbose) {
        logger.info(String.format(" Moving files %d to %d", start, endRange));
    }
    try {
        File newDirectory = new File(newDirectoriesName + start + "_" + endRange);
        if (newDirectory.exists()) {
            for (int j = 2; j < 1000; j++) {
                newDirectory = new File(newDirectoriesName + "vers" + j + "_" + start + "_" + endRange);
                if (!newDirectory.exists()) {
                    break;
                }
            }
            if (newDirectory.exists()) {
                throw new IllegalArgumentException(String.format("Versioning failure: \n\
                    all directory names from %s_vers_2_%d_%d to %s_vers_999_%d_%d were taken.",
                        newDirectoriesName, start, end, newDirectoriesName, start, end));
            }
        }
        if (!newDirectory.mkdir()) { // Directory is made, and then false is returned if a failure occurred.
            throw new IOException("Error in creating directory.");
        }
        for (int j = i; j < endRange; j++) {
            File filej = files.get(j);
            FileUtils.moveFileToDirectory(filej, newDirectory, false);
        }
    } catch (IOException ex) {
        logger.severe(String.format(" Error in splitting files to directory: %s", ex.toString()));
    }
}
