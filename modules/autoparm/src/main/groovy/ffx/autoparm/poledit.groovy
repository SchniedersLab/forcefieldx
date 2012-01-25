import groovy.util.CliBuilder
//POLEDIT
//Give the gdmaoutfile and a peditinfile
//Poledit will create a .xyz file and .key file containing fixed multipole parameters and polarizabilities

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc poledit <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null) {
    return cli.usage();
}
String gdmaoutfname = arguments.get(0);
String peditinfname = arguments.get(1);
poledit(gdmaoutfname,peditinfname);



