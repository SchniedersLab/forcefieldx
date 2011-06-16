String xyzfilename = args[0];
double eps = Double.parseDouble(args[1]);
if(args.size() > 3 && args[2].equals("-k")){
    String keyfname = args[3];
    minimize_2(xyzfilename,eps,keyfname);
}
else{
    minimize_2(xyzfilename,eps,null);
}
