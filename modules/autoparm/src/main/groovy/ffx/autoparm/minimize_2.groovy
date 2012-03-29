import groovy.util.CliBuilder

def cli = new CliBuilder(usage:' ffxc minimize_2 [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null) {
    return cli.usage();
}

String xyzfilename = arguments.get(0);
double eps = Double.parseDouble(arguments.get(1));
if(arguments.size() > 3 && arguments.get(2).equals("-k")){
    String keyfname = arguments.get(3);
    minimize_2(xyzfilename,eps,keyfname);
}
else{
    minimize_2(xyzfilename,eps,null);
}
