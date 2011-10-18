import groovy.util.CliBuilder

//SUPERPOSE
//Give two xyzfiles
//Program will print out information about the distance between the two
def cli = new CliBuilder(usage:' ffxc poledit <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null || arguments.size() != 1) {
    return cli.usage();
}
String fname1 = arguments.get(0);
String fname2 = arguments.get(1);

superpose(fname1,fname2);
