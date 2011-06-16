//ENERGY
//Gives system potential and prints out system multipoles
//filename = xyzfile
//use the -k flag to specify a keyfile
//use the -o flag to specify options
//options are: p -> prints out system multipole information; d-> prints out detailed information about the interactions (right now it only prints out information about tor-angles but that can be
//changed easily
String xyzfilename = args[0];
String keyfname = null;
String options = null;
if(args.size() == 5){
    if(args[1].equals("-k")){
        keyfname = args[2];
        options = args[4];
        analyze(xyzfilename,keyfname,options)
    }
    else if(args[1].equals("-o")){
        keyfname = args[4];
        options = args[2];
        analyze(xyzfilename,keyfname,options)
    }
}
else if(args.size() == 3){
    if(args[1].equals("-k")){
        keyfname = args[2];
        analyze(xyzfilename,keyfname,options)
    }
    else if(args[1].equals("-o")){
        options = args[2];
        analyze(xyzfilename,keyfname,options)
    }
}
else{
    analyze(xyzfilename,keyfname,options);
}
