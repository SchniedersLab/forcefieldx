// Convert from PRM to Property

// Apache Imports
import org.apache.commons.configuration.CompositeConfiguration;

// Groovy Imports
import groovy.util.CliBuilder;

// FFX Imports
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.utilities.Keyword;

// Things below this line normally do not need to be changed.
// ===============================================================================================

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc prmToProperty <prm> [prm] ...');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() < 1) {
    return cli.usage();
}

// Read in the command line file.
String xyzname = arguments.get(0);
CompositeConfiguration properties = Keyword.loadProperties(null);
properties.setProperty("parameters", xyzname);
ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
ForceField forceField = forceFieldFilter.parse();

int prms = arguments.size();
for (int i=1; i<prms; i++) {
    xyzname = arguments.get(i);
    properties = Keyword.loadProperties(null);
    properties.setProperty("parameters", xyzname);
    forceFieldFilter = new ForceFieldFilter(properties);
    ForceField forceField2 = forceFieldFilter.parse();
    forceField.append(forceField2);
}

if (forceField != null) {
    forceField.print();
}
