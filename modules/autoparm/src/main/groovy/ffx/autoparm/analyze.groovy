import groovy.util.CliBuilder

//ENERGY
//Gives system potential and prints out system multipoles
//filename = xyzfile
//use the -k flag to specify a keyfile
//use the -o flag to specify options
//options are: p -> prints out system multipole information; d-> prints out detailed information about the interactions (right now it only prints out information about tor-angles but that can be
//changed easily

// Create the command line parser.
def cli = new CliBuilder(usage:' ffxc analyze [options] <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);

List<String> arguments = options.arguments();
if (options.h || arguments == null) {
    return cli.usage();
}


String xyzfilename = arguments.get(0);
String keyfname = null;
String options1 = null;



if(arguments.size() == 5){
    if(arguments.get(1).equals("-k")){
        keyfname = arguments.get(2);
        options1 = arguments.get(4);
        analyze(xyzfilename,keyfname,options1)
    }
    else if(arguments.get(1).equals("-o")){
        keyfname = arguments.get(4);
        options1 = arguments.get(2);
        analyze(xyzfilename,keyfname,options1)
    }
}
else if(arguments.size() == 3){
    if(arguments.get(1).equals("-k")){
        keyfname = arguments.get(2);
        analyze(xyzfilename,keyfname,options1)
    }
    else if(arguments.get(1).equals("-o")){
        options1 = arguments.get(2);
        analyze(xyzfilename,keyfname,options1)
    }
}
else{
    analyze(xyzfilename,keyfname,options1);
}
