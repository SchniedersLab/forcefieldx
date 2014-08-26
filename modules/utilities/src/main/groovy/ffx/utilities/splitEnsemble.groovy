/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2012.
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
 */

package ffx.utilities

// Apache Commons Imports
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.FileUtils;

// Groovy Imports
import groovy.util.CliBuilder;

int startModel = 1;
int endModel = -1;
boolean finished = false;
String suffix = "_split";
boolean addHeader = false;
boolean addListFile = true;

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc utilities.splitEnsemble [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
cli.x(longOpt:'all', args:1, argName:'1', 'Split out all PDB files starting from this model (over-rides other options.');
cli.s(longOpt:'start', args:1, argName:'-1', 'Starting model to split into separate PDB file.');
cli.f(longOpt:'finish', args:1, argName:'-1', 'Finishing model to split into separate PDB file.');
cli.X(longOpt:'suffix', args:1, argName: '_split', 'output suffix');
cli.d(longOpt:'header', args:1, argName: 'false', 'Whether to attach header lines to each PDB file.');
cli.dF(longOpt:'headerFileSource', args:1, argName: 'file', 'May use this file as source of header lines instead of ensemble file.');
cli.l(longOpt: 'listFiles', args:1, argName: 'true', 'Additionally output text file with the names of all produced PDB files as _suffix_modelList.txt');

def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

// Read in command line.
String filename = arguments.get(0);
String targetName = FilenameUtils.removeExtension(filename);

// Sets where we will start and stop
if (options.x) {
    startModel = Integer.parseInt(options.x);
    endModel = Integer.MAX_VALUE;
} else {
    if (options.s && options.f) {
        startModel = Integer.parseInt(options.s);
        endModel = Integer.parseInt(options.f);
        if (startModel > endModel) {
            logger.severe(" Start must be less than finish");
        }
    } else {
        logger.severe(" Either -x or both -s and -f must be specified.");
    }
}

// Sets the desired suffix for the output files.
if (options.X) {
    suffix = options.X;
}
targetName = targetName.concat(suffix);

if (options.d) {
    addHeader = Boolean.parseBoolean(options.d);
}


/*if (nextLine.startsWith("MODEL  ") || nextLine.startsWith("ATOM   ") || nextLine.startsWith("SSBOND  ") 
|| nextLine.startsWith("ANISOU  ") || nextLine.startsWith("CONECT ") || nextLine.startsWith("HET  ") 
|| nextLine.startsWith("HETATM  ") || nextLine.startsWith("FORMUL  ") || nextLine.startsWith("CRYST1")) {*/

try {
    File headerFile;
    BufferedReader fileReader = new BufferedReader(new FileReader(filename));
    int counter = 1;
    boolean reachedFirstModel = false;
    if (options.dF) {
        addHeader = true;
        headerFile = new File (targetName + "_head.tmp.txt");
        FileUtils.forceDeleteOnExit(headerFile);
        BufferedWriter headBw = new BufferedWriter(new FileWriter(headerFile, false));
        BufferedReader headReader = new BufferedReader(new FileReader(options.dF));
        boolean headerRead = false;
        while (!headerRead) {
            String nextLine = headReader.readLine();
            if (nextLine == null || nextLine.startsWith("CRYST1") || nextLine.startsWith("ATOM    ")
                || nextLine.contains("END OF HEADER") || nextLine.startsWith("MODEL  ")) {
                headerRead = true;
            } else {
                headBw.write(nextLine);
                headBw.newLine();
                headBw.flush();
            }
        }
        headBw.close();
        headReader.close();
        logger.info(" Header read successfully.");
    } else if (addHeader) {
        headerFile = new File(targetName + "_head.tmp.txt");
        FileUtils.forceDeleteOnExit(headerFile);
        BufferedWriter headBw = new BufferedWriter(new FileWriter(headerFile, false));
        while (!reachedFirstModel) {
            String nextLine = fileReader.readLine();
            if (nextLine.startsWith("MODEL   ")) {
                reachedFirstModel = true;
            } else {
                headBw.write(nextLine);
                headBw.newLine();
                headBw.flush();
            }
        }
        headBw.close();
        logger.info(" Header read successfully.");
    } else {
        while (!reachedFirstModel) {
            String nextLine = fileReader.readLine();
            reachedFirstModel = (nextLine.startsWith("MODEL  ") ? true : false);
        }
    }
    boolean startModelReached = (startModel < 2 ? true : false);
    while (!startModelReached) {
        String nextLine = fileReader.readLine();
        if (nextLine.startsWith("MODEL   ")) {
            ++counter;
            startModelReached = (counter >= startModel ? true : false);
        }
    }
    // get to appropriate model # (not coded)
    // While still in the range for writing sub-files:
    if (options.l) {
        addListFile = Boolean.parseBoolean(options.l);
    }
    File listOfFiles;
    BufferedWriter listWriter;
    if (addListFile) {
        listOfFiles = new File (targetName + "_modelList.txt");
        logger.info(String.format(" Creating file list (with the names of all created files): %s", targetName + "_modelList.txt"));
        listWriter = new BufferedWriter(new FileWriter(listOfFiles, false));
    }
    while (counter <= endModel && !finished) {
        File newModel = new File(targetName + "_" + counter + ".pdb");
        if (addListFile) {
            listWriter.write(targetName + "_" + counter + ".pdb");
            listWriter.newLine();
            listWriter.flush();
        }
        BufferedWriter bw;
        if (addHeader) {
            FileUtils.copyFile(headerFile, newModel);
            bw = new BufferedWriter(new FileWriter(newModel, true));
        } else {
            bw = new BufferedWriter(new FileWriter(newModel, false));
        }
        // Creates a file, creates a file writer which writes to the file, and creates a buffered writer which writes to the file writer.
        // I suspect the buffered writer helps in efficiency: FileWriter's methods can be bulky, so the buffered writer breaks it down in a buffer which is more efficient.
        boolean modelFinished = false;
        while (!modelFinished) {
            // While we have not finished this model's file, keep on reading in lines, skip MODEL lines, break ot if ENDMDL, otherwise write.
            String thisLine = fileReader.readLine();
            if (thisLine.startsWith("MODEL  ")) {
                // Skip. Should not be reached anyways, though.
            } else if (thisLine.contains("ENDMDL")) {
                modelFinished = true;
            } else {
                bw.write(thisLine);
                bw.newLine();
                bw.flush();
                // Write the line, write a new line character, and flush the stream. I'm not 100% sure on what flushing does.
            }
        }
        // Close the writers, increment the counter, check if there's more to the file.
        bw.close();
        logger.info(String.format(" Model %d written as file.", counter));
        ++counter;
        // Mark our place, check if there's anything more (it returns null if there's nothing), and reset.
        //fileReader.mark(5000);
        if (fileReader.readLine() == null) {
            finished = true;
        }
        //fileReader.reset();
    }
    listWriter.close();
    fileReader.close();
} catch (IOException ex) {
    logger.severe(String.format(" Exception writing to files: %s", ex.toString()));
}