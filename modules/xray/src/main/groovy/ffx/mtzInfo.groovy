// Groovy Imports
import groovy.util.CliBuilder;

import java.io.File;

import ffx.xray.MTZFilter;

// Things below this line normally do not need to be changed.
// ===============================================================================================

def today = new Date();
logger.info(" " + today);
logger.info(" command line variables:");
logger.info(" " + args + "\n");

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc mtzInfo [options] <mtzfilename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}

String mtzfile = arguments.get(0);

File file = new File(mtzfile);
if (!file.exists()){
  println("File: " + mtzfile + " not found!");
  return;
}

MTZFilter mtzfilter = new MTZFilter();
mtzfilter.getReflectionList(file);
mtzfilter.print_header();
