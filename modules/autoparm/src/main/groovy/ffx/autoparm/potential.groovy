import groovy.util.CliBuilder

//POTENTIAL
//Supply a choice (1 - 4)
//1. Get QM Potential from a Gaussian CUBE File
//2. Calculate the Model Potential for a System
//3. Compare a Model Potential to a Target Grid
//4. Fit Electrostatic Parameters to Target Grid
//if choice == 1, then supply the cube filename
//if choice == 2 || choice == 3, then supply the xyz filename
//if choice == 4, then supply the xyz filename and the eps

def cli = new CliBuilder(usage:' ffxc potential <filename>');
cli.h(longOpt:'help', 'Print this help message.');
def options = cli.parse(args);
List<String> arguments = options.arguments();
if (options.h || arguments == null) {
    return cli.usage();
}

//Get XYZ File
Double eps = null;
Integer choice = Integer.parseInt(arguments.get(0));
String filename = arguments.get(1);
int out_type = 5;
if(choice == 4){
	eps = Double.parseDouble(arguments.get(2));
}
potential(choice, filename, eps);
