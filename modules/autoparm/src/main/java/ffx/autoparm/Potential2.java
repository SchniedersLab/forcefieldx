package ffx.autoparm;
import static ffx.numerics.VectorMath.diff;
import static ffx.numerics.VectorMath.r;
//import static org.junit.Assert.*;

import java.io.*;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Random;
import java.util.logging.Logger;
import java.util.StringTokenizer;

//import org.junit.Test;
//import org.junit.runner.RunWith;
//import org.junit.runners.Parameterized;
//import org.junit.runners.Parameterized.Parameters;

import org.apache.commons.configuration.CompositeConfiguration;

import edu.rit.pj.ParallelTeam;
import ffx.crystal.Crystal;
import ffx.numerics.LBFGS;
import ffx.numerics.OptimizationListener;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Utilities;
import ffx.autoparm.PME_2;
import ffx.potential.nonbonded.ParticleMeshEwald;
import ffx.potential.nonbonded.VanDerWaals;
import ffx.potential.nonbonded.ParticleMeshEwald.Polarization;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.autoparm.ForceFieldFilter_2;
import ffx.potential.parsers.XYZFilter;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.ForceField.ForceFieldDouble;
import ffx.potential.parameters.ForceField.ForceFieldString;
import ffx.utilities.Keyword;
import ffx.numerics.OptimizationListener;
import ffx.numerics.LineSearch.LineSearchResult;

public class Potential2 implements OptimizationListener {
	private String info;
	private ArrayList<Integer> ipgrid = new ArrayList<Integer>();
	private double target_grid[][][];//nSymm, nAtoms, 4
	private double pot_grid[][][];//nSymm, nAtoms, 4
	private Atom atoms[];
	private int nAtoms;
	private File structure_key;
	private File structure_xyz;
	private File structure_cube;
	private File structure_prm;
	private MolecularAssembly molecularAssembly;
	private ArrayList<String> key = new ArrayList<String>();
	private ForceField forceField;
	private Crystal crystal;
	private int nSymm;
	private double x[];
	private double grad[];
	private double scaling[];
	private boolean done = false;
	private boolean terminate = false;
	private long time;
	private double grms;
	private int nSteps;
	private int nvars;
	private static final Logger logger = Logger.getLogger(Potential2.class.getName());
	public static final double BOHR = 0.52917720859;
	private PME_2 pme = null;

	InputStreamReader stdinput = new InputStreamReader(System.in);
	BufferedReader stdreader = new BufferedReader(stdinput);
	
	public Potential2() throws IOException{
		
		System.out.println("Would you like to:\n" +
				"1. Get QM Potential from a Gaussian CUBE File\n" +
				"2. Calculate the Model Potential for a System\n" +
				"3. Compare a Model Potential to a Target Grid\n" +
				"4. Fit Electrostatic Parameters to Target Grid");

		int choice = Integer.parseInt(stdreader.readLine());

		if(choice == 1){
			do_cube(null);
			System.exit(0);
		}

		boolean done = false;
		String xyzfname = null;
		System.out.println("Enter XYZ Filename: ");
		while(!done){
			xyzfname = stdreader.readLine();
			structure_xyz = new File(xyzfname);
			if(structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead()){
				done = true;
			}
			else{
				System.out.println("Couldn't find file. Please reenter filename: ");
			}
		}

		int n = 1;
		String oxyzfname = null;
		String old = xyzfname;
		while(structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead()){
			oxyzfname = xyzfname;
			n++;
			xyzfname = old;
			xyzfname = xyzfname+"_"+Integer.toString(n);
			structure_xyz = new File(xyzfname);
		}

		structure_xyz = new File(oxyzfname);
		int index = oxyzfname.lastIndexOf(".");
		String name = oxyzfname.substring(0, index);
		String keyfname = name+".key";
		structure_key = new File(keyfname);
		if(!(structure_key != null && structure_key.exists() && structure_key.canRead())){
			System.out.println("Enter Key Filename: ");
			done = false;
			while(!done){
				keyfname = stdreader.readLine();
				structure_key = new File(keyfname);
				if(structure_key != null && structure_key.exists() && structure_key.canRead()){
					done = true;
				}
				else{
					System.out.println("Couldn't find file. Please reenter filename: ");
				}
			}
			structure_key = new File(keyfname);
		}
		n = 1;
		String okeyfname = null;
		old = keyfname;
		while(structure_key != null && structure_key.exists() && structure_key.canRead()){
			okeyfname = keyfname;
			n++;
			keyfname = old;
			keyfname = keyfname+"_"+Integer.toString(n);
			structure_key = new File(keyfname);
		}

		structure_key = new File(okeyfname);

		molecularAssembly = new MolecularAssembly(name);
		molecularAssembly.setFile(structure_xyz);
		CompositeConfiguration properties = Keyword.loadProperties(structure_key);
		
		
		//ForceFieldFilter_2 forceFieldFilter = new ForceFieldFilter_2(properties, structure_key);
		//Decides difference between using prm file (below) or key file (above) for parameters
		ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
		
		
		forceField = forceFieldFilter.parse();
		//System.out.println(forceField.getMultipoleType("401 401 403"));
		molecularAssembly.setForceField(forceField);
		XYZFilter xyzFilter = new XYZFilter(structure_xyz,molecularAssembly,forceField,properties);
		boolean expectedReturn = true;
		boolean actualReturn = xyzFilter.readFile();
		//assertEquals(info, expectedReturn, actualReturn);
		Utilities.biochemistry(molecularAssembly, xyzFilter.getAtomList());
		//Had to comment out something in finalize method
		molecularAssembly.finalize(true);
		atoms = molecularAssembly.getAtomArray();
		nAtoms = atoms.length;
		ParallelTeam parallelTeam = new ParallelTeam();
		crystal = create_crystal(forceField, atoms);
		nSymm = crystal.spaceGroup.getNumberOfSymOps();
		VanDerWaals vanderWaals = new VanDerWaals(forceField, atoms, crystal, parallelTeam);
		store_key_file(structure_key);
		pme = new PME_2(forceField, atoms, crystal, parallelTeam, vanderWaals.getNeighborLists(), key);


		if(choice == 2){
			double pgrid[][][] = gen_pot_grid(null, atoms, 1);
			String potfilename = name+".pot";
			String ofname = potfilename;
			File outf = new File(potfilename);
			n = 1;
			while(outf != null && outf.exists() && outf.canRead()){
				potfilename = ofname;
				n++;
				potfilename = potfilename+"_"+Integer.toString(n);
				outf = new File(potfilename);
			}
			String gridfilename = name+".grid";
			String ofgname = gridfilename;
			File outfg = new File(gridfilename);
			n = 1;
			while(outfg != null && outfg.exists() && outfg.canRead()){
				gridfilename = ofgname;
				n++;
				gridfilename = gridfilename+"_"+Integer.toString(n);
				outfg = new File(gridfilename);
			}
			output_pgrid(outf, outfg, pgrid);
			System.exit(0);
		}
		
		String cubefname = name+".cube";
		if(!(structure_cube != null && structure_cube.exists() && structure_cube.canRead())){
			System.out.println("Enter Cube Filename: ");
			done = false;
			while(!done){
				cubefname = stdreader.readLine();
				structure_cube = new File(cubefname);
				if(structure_cube != null && structure_cube.exists() && structure_cube.canRead()){
					done = true;
				}
				else{
					System.out.println("Couldn't find file. Please reenter filename: ");
				}
			}
		}
		structure_cube = new File(cubefname);
		n = 1;
		String ocubefname = null;
		old = cubefname;
		while(structure_cube != null && structure_cube.exists() && structure_cube.canRead()){
			ocubefname = cubefname;
			n++;
			cubefname = old;
			cubefname = cubefname+"_"+Integer.toString(n);
			structure_cube = new File(cubefname);
		}

		structure_cube = new File(ocubefname);
		
		
		if(choice == 3){
			target_grid = gen_pot_grid(structure_cube, atoms, 0);
			pme.init_prms();
//			System.out.println("Output potential at every point? Y/N");
//			int ans = stdreader.read();
//			if(ans == 'Y'){
//				double pgrid[][][] = gen_pot_grid(structure_cube,atoms, 1);
//				String potfilename = name+".pot";
//				String ofname = potfilename;
//				File outf = new File(potfilename);
//				n = 1;
//				while(outf != null && outf.exists() && outf.canRead()){
//					potfilename = ofname;
//					n++;
//					potfilename = potfilename+"_"+Integer.toString(n);
//					outf = new File(potfilename);
//				}
//				output_pgrid(outf, pgrid);
//			}
			output_stats();
			System.exit(0);
		}
		else if(choice == 4){
			target_grid = gen_pot_grid(structure_cube, atoms, 0);
			pme.set_target_grid(target_grid);
			double a = avgrms();
			output_stats();
			int m = 7;
			double eps = .1;
			nvars = pme.getNumberOfVariables();
			x = new double[nvars];
			grad = new double[nvars];
			scaling = new double[nvars];
			pme.getCoordinates(x);
			double e = pme.energyAndGradient(x, grad);
			time = -System.nanoTime();
			int status = 0;
			long at = System.nanoTime();
			status = LBFGS.minimize(nvars, m, x, e, grad, eps, pme, this);
			long bt = System.nanoTime();
			System.out.println("TOTAL TIME "+(bt - at)*1e-9+"\n");
			switch (status) {
			case 0:
				logger.info(String.format("\n Optimization achieved convergence criteria: %8.5f\n", grms));
				break;
			case 1:
				logger.info(String.format("\n Optimization terminated at step %d.\n", nSteps));
				break;
			default:
				logger.warning("\n Optimization failed.\n");
			}
			pme.init_prms();
			output_stats();
			double b = avgrms();
			System.out.println("a = "+a+" b = "+b+" diff = "+(a - b));
		}
	}
	
	/**
	 * key_filename and cube_filename can be null if they are in the same directory as xyz_filename
	 * eps is only needed if choice == 4, otherwise leave as null
	 * cube_filename is needed if choice == 1
	 * @param choice
	 * @param xyz_filename
	 * @param key_filename
	 * @param cube_filename
	 * @param eps
	 * @throws IOException
	 */
	public Potential2(int choice, String xyz_filename, String cube_filename, Double eps) throws IOException{
		if(choice == 1){
			do_cube(cube_filename);
			System.exit(0);
		}
		structure_xyz = new File(xyz_filename);
		if(!(structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead())){
			System.out.println("Couldn't find xyz file");
			System.exit(1);
		}
		int n = 1;
		String oxyzfname = null;
		String old = xyz_filename;
		while(structure_xyz != null && structure_xyz.exists() && structure_xyz.canRead()){
			oxyzfname = xyz_filename;
			n++;
			xyz_filename = old;
			xyz_filename = xyz_filename+"_"+Integer.toString(n);
			structure_xyz = new File(xyz_filename);
		}
		structure_xyz = new File(oxyzfname);
		int index = oxyzfname.lastIndexOf(".");
		String name = oxyzfname.substring(0, index);
		String keyfname = name+".key";
		structure_key = new File(keyfname);
		while(!(structure_key != null && structure_key.exists() && structure_key.canRead())){
			System.out.println("Enter the Key Filename: ");
			keyfname = stdreader.readLine();
			structure_key = new File(keyfname);
			if(!(structure_key != null && structure_key.exists() && structure_key.canRead())){
				System.out.println("Couldn't find key file");
			}
		}

		n = 1;
		String okeyfname = null;
		old = keyfname;
		while(structure_key != null && structure_key.exists() && structure_key.canRead() && structure_key.length() != 0){
			okeyfname = keyfname;
			n++;
			keyfname = old;
			keyfname = keyfname+"_"+Integer.toString(n);
			structure_key = new File(keyfname);
		}

		structure_key = new File(okeyfname);
		
//		String prmfname = name+".prm";
//		structure_prm = new File(prmfname);
//		while(!(structure_prm != null && structure_prm.exists() && structure_prm.canRead())){
//			System.out.println("Enter the prm Filename: ");
//			prmfname = stdreader.readLine();
//			structure_prm = new File(prmfname);
//			if(!(structure_prm != null && structure_prm.exists() && structure_prm.canRead())){
//				System.out.println("Couldn't find prm file");
//			}
//		}
//
//		n = 1;
//		String oprmfname = null;
//		old = prmfname;
//		while(structure_prm != null && structure_prm.exists() && structure_prm.canRead() && structure_prm.length() != 0){
//			oprmfname = prmfname;
//			n++;
//			prmfname = old;
//			prmfname = prmfname+"_"+Integer.toString(n);
//			structure_prm = new File(prmfname);
//		}
//
//		structure_prm = new File(oprmfname);

		molecularAssembly = new MolecularAssembly(name);
		molecularAssembly.setFile(structure_xyz);
		CompositeConfiguration properties = Keyword.loadProperties(structure_key);
		ForceFieldFilter_2 forceFieldFilter = new ForceFieldFilter_2(properties, structure_key);
		//Decides difference between using prm file (below) or key file (above) for parameters
		//ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
		forceField = forceFieldFilter.parse();
		molecularAssembly.setForceField(forceField);
		XYZFilter xyzFilter = new XYZFilter(structure_xyz,molecularAssembly,forceField,properties);
		xyzFilter.readFile();
		Utilities.biochemistry(molecularAssembly, xyzFilter.getAtomList());
		molecularAssembly.finalize(true);
		atoms = molecularAssembly.getAtomArray();
		nAtoms = atoms.length;
		ParallelTeam parallelTeam = new ParallelTeam();
		crystal = create_crystal(forceField, atoms);
		nSymm = crystal.spaceGroup.getNumberOfSymOps();
		VanDerWaals vanderWaals = new VanDerWaals(forceField, atoms, crystal, parallelTeam);
		//RENAME
		//store_key_file(structure_prm);
		store_key_file(structure_key);
		pme = new PME_2(forceField, atoms, crystal, parallelTeam, vanderWaals.getNeighborLists(), key);

		if(choice == 2){
			double pgrid[][][] = gen_pot_grid(null, atoms, 1);
			String potfilename = name+".pot";
			String ofname = potfilename;
			File outf = new File(potfilename);
			n = 1;
			while(outf != null && outf.exists() && outf.canRead()){
				potfilename = ofname;
				n++;
				potfilename = potfilename+"_"+Integer.toString(n);
				outf = new File(potfilename);
			}
			String gridfilename = name+".grid";
			String ofgname = gridfilename;
			File outfg = new File(gridfilename);
			n = 1;
			while(outfg != null && outfg.exists() && outfg.canRead()){
				gridfilename = ofgname;
				n++;
				gridfilename = gridfilename+"_"+Integer.toString(n);
				outfg = new File(gridfilename);
			}
			output_pgrid(outf, outfg, pgrid);
			System.exit(0);
		}
		
		String cubefname = name+".cube";
		structure_cube = new File(cubefname);
		while(!(structure_cube != null && structure_cube.exists() && structure_cube.canRead())){
			System.out.println("Enter the Cube Filename: ");
			cubefname = stdreader.readLine();
			structure_cube = new File(cubefname);
			if(!(structure_cube != null && structure_cube.exists() && structure_cube.canRead())){
				System.out.println("Couldn't find cube file");
			}
		}
		n = 1;
		String ocubefname = null;
		old = cubefname;
		while(structure_cube != null && structure_cube.exists() && structure_cube.canRead()){
			ocubefname = cubefname;
			n++;
			cubefname = old;
			cubefname = cubefname+"_"+Integer.toString(n);
			structure_cube = new File(cubefname);
		}
		structure_cube = new File(ocubefname);
		
		if(choice == 3){
			target_grid = gen_pot_grid(structure_cube, atoms, 0);
			pme.init_prms();
			output_stats();
			System.exit(0);
		}
		else if(choice == 4){
			target_grid = gen_pot_grid(structure_cube, atoms, 0);
			pme.set_target_grid(target_grid);
			output_stats();
			int m = 7;
			nvars = pme.getNumberOfVariables();
			x = new double[nvars];
			grad = new double[nvars];
			scaling = new double[nvars];
			pme.getCoordinates(x);
			double e = pme.energyAndGradient(x, grad);
			time = -System.nanoTime();
			int status = 0;
			long at = System.nanoTime();
			status = LBFGS.minimize(nvars, m, x, e, grad, eps.doubleValue(), pme, this);
			long bt = System.nanoTime();
			System.out.println("TIME FOR LBFGS: "+(bt - at)*1e-9+" seconds\n");
			switch (status) {
			case 0:
				logger.info(String.format("\n Optimization achieved convergence criteria: %8.5f\n", grms));
				break;
			case 1:
				logger.info(String.format("\n Optimization terminated at step %d.\n", nSteps));
				break;
			default:
				logger.warning("\n Optimization failed.\n");
			}
			pme.init_prms();
			output_stats();
			//output_keyfile(pme.getallmpoles(x), prmfname);
			output_keyfile(pme.getallmpoles(x), keyfname);
		}
	}
	
	public void output_keyfile(double[] mpoles, String outfname) throws IOException{
		File outf = new File(outfname);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outf)));
		int r = 0;
		DecimalFormat myFormatter = new DecimalFormat(" ##########0.00000;-##########0.00000");
    	ArrayList<Integer> types = new ArrayList<Integer>();
    	int pos;
    	int polelen = 10;
    	//for traceless manipulations
    	int a = 4, b = 5, c = 6;
    	if(!pme.fitmpl){
    		polelen = polelen - 1;
    		a -= 1;
    		b -= 1;
    		c -= 1;
    	}
    	if(!pme.fitdpl){
    		polelen = polelen - 3;
    		a -= 3;
    		b -= 3;
    		c -= 3;
    	}
    	if(!pme.fitqdpl){
    		polelen = polelen - 5;
    	}

		//maintain traceless quadrupole at each multipole site
    	if(pme.fitqdpl){
    		double sum = 0;
    		double big = 0;
    		for(int i = 0; i < (mpoles.length)/polelen;i++){
    			sum = mpoles[(i*polelen) + a] + mpoles[(i*polelen) + b] + mpoles[(i*polelen) + c];
    			big = Math.max(Math.abs(mpoles[(i*polelen) + a]), Math.max(Math.abs(mpoles[(i*polelen) + b]), Math.abs(mpoles[(i*polelen) + c])));
    			int k = 0;
    			if(big == Math.abs(mpoles[(i*polelen) + a])) k = a;
    			else if(big == Math.abs(mpoles[(i*polelen) + b])) k = b;
    			else if(big == Math.abs(mpoles[(i*polelen) + c])) k = c;
    			if(k != 0){
    				mpoles[(i*polelen) + k] = mpoles[(i*polelen) + k] - sum;
    			}
    		}
    	}

    	for(int i = 0; i < nAtoms; i++){
    		Atom ai = atoms[i];
			if(!types.contains(ai.getType())){
				types.add(ai.getType());
			}
    	}
    	

		for(int i = 0; i < key.size(); i++){
			if(!key.get(i).toUpperCase().contains("MULTIPOLE")){
				bw.write(key.get(i)+"\n");
			}
			else if(key.get(i).toUpperCase().contains("MULTIPOLE")){
				pos = types.indexOf(Integer.parseInt(key.get(i).trim().split(" +")[1]));
				if(pme.fitmpl){
					bw.write(key.get(i).trim().split(" +")[0]+"   "+key.get(i).trim().split(" +")[1]+"   "+key.get(i).trim().split(" +")[2]+"   "+key.get(i).trim().split(" +")[3]+"\t"+myFormatter.format(mpoles[pos*polelen + r])+"\n");
					r = r + 1;
				}
				else{
					bw.write(key.get(i)+"\n");
				}
				if(pme.fitdpl){
					bw.write("                             \t\t\t\t"+myFormatter.format(mpoles[pos*polelen + r]/BOHR)+"   "+myFormatter.format(mpoles[pos*polelen + r+1]/BOHR)+"   "+myFormatter.format(mpoles[pos*polelen + r+2]/BOHR)+"\n");
					r = r + 3;
				}
				else{
					bw.write(key.get(i+1)+"\n");
				}
				if(pme.fitqdpl){
					bw.write("                             \t\t\t\t"+myFormatter.format(mpoles[pos*polelen + r]/(BOHR*BOHR))+"\n");
					bw.write("                             \t\t\t\t"+myFormatter.format(mpoles[pos*polelen + r+3]/(BOHR*BOHR))+"   "+myFormatter.format(mpoles[pos*polelen + r+1]/(BOHR*BOHR))+"\n");
					bw.write("                            \t\t\t\t"+myFormatter.format(mpoles[pos*polelen + r+4]/(BOHR*BOHR))+"   "+myFormatter.format(mpoles[pos*polelen + r+5]/(BOHR*BOHR))+"   "+myFormatter.format(mpoles[pos*polelen + r+2]/(BOHR*BOHR))+"\n");
					r = r + 6;
				}
				else{
					bw.write(key.get(i+2)+"\n");
					bw.write(key.get(i+3)+"\n");
					bw.write(key.get(i+4)+"\n");
				}
				i = i + 4;
				r = 0;
			}
		}
		bw.close();
		System.out.println("Keyfile written to "+outfname);
	}
	
	public void do_cube(String cfname) throws IOException{
		int n = 1;
		final double HARTREE = 627.5094688;
		String fname = cfname;
		if(fname == null){
			System.out.println("Enter the Gaussian CUBE Filename: ");
			fname = stdreader.readLine();
		}
		String or_fname = fname;
		File outf = new File(fname);

		while(outf != null && outf.exists() && outf.canRead()){
			fname = or_fname;
			n++;
			fname = fname+"_"+Integer.toString(n);
			outf = new File(fname);
		}
		File orf = new File(or_fname);
		String line;
		String ln[];
		double x, y, z, pot;
		int num = 0;
		int nAtoms;
		int nPoints;
		if (orf != null && orf.exists()
				&& orf.canRead()) {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(orf)));
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outf)));
			//read some initial stuff
			String title = br.readLine();
			br.readLine();
			nAtoms = Integer.parseInt(br.readLine().trim().split(" +")[0]);
			nPoints = Integer.parseInt(br.readLine().trim().split(" +")[0]);
			for(int i = 0; i < nAtoms + 2; i++){
				br.readLine();
			}
			double gridpoints[][] = new double[nPoints][4];
			while((line = br.readLine()) != null){
				ln = line.trim().split(" +");
				x = Double.parseDouble(ln[0]);
				gridpoints[num][0] = x;
				y = Double.parseDouble(ln[1]);
				gridpoints[num][1] = y;
				z = Double.parseDouble(ln[2]);
				gridpoints[num][2] = z;
				pot = Double.parseDouble(ln[3]) * HARTREE;
				gridpoints[num][3] = pot;
				num++;
			}
			br.close();
			bw.write("\t"+Integer.toString(nPoints)+" "+title+"\n");
			DecimalFormat myFormatter = new DecimalFormat("###########0.0000000");
			DecimalFormat myFormatter2 = new DecimalFormat("###########0.0000");
			String output;
			for(int i = 0; i < nPoints; i++){
				output = myFormatter.format(gridpoints[i][0])+"\t"+myFormatter.format(gridpoints[i][1])+"\t"+myFormatter.format(gridpoints[i][2])+"\t"+myFormatter2.format(gridpoints[i][3]);
				//bw.write("\t"+i+"\t"+gridpoints[i][0]%-12.6f+"\t"+gridpoints[i][1]%-12.6f+"\t"+gridpoints[i][2]%-12.6f+"\t"+gridpoints[i][3]%-12.6f+"\n");	
				bw.write("\t"+(i+1)+"\t"+output+"\n");
			}
			bw.close();
			System.out.println("Electrostatic Potential Written to File: " + fname);
		}	
	}

	public Crystal create_crystal(ForceField forceField, Atom atoms[]){
		// Determine the unit cell dimensions and Spacegroup
		String spacegroup;
		double a,b,c,alpha,beta,gamma;
		boolean aperiodic;
		int nAtoms = atoms.length;
		double ewaldOff = forceField.getDouble(ForceFieldDouble.EWALD_CUTOFF, 7.0);
		try {
			spacegroup = forceField.getString(ForceFieldString.SPACEGROUP);
			a = forceField.getDouble(ForceFieldDouble.A_AXIS);
			b = forceField.getDouble(ForceFieldDouble.B_AXIS, a);
			c = forceField.getDouble(ForceFieldDouble.C_AXIS, a);
			alpha = forceField.getDouble(ForceFieldDouble.ALPHA, 90.0);
			beta = forceField.getDouble(ForceFieldDouble.BETA, 90.0);
			gamma = forceField.getDouble(ForceFieldDouble.GAMMA, 90.0);
			aperiodic = false;
		} catch (Exception e) {
			//logger.info(" The system will be treated as aperiodic.");
			aperiodic = true;
			spacegroup = "P1";
			/**
			 * Search all atom pairs to find the largest pair-wise distance.
			 */
			double xr[] = new double[3];
			double maxr = 0.0;
			for (int i=0; i < nAtoms - 1; i++) {
				double[] xi =  atoms[i].getXYZ();
				for (int j=i+1; j < nAtoms; j++) {
					double[] xj = atoms[j].getXYZ();
					diff(xi, xj, xr);
					double r = r(xr);
					if (r > maxr) {
						maxr = r;
					}
				}
			}
			a = 2.0 * (maxr + ewaldOff);
			b = a;
			c = a;
			alpha = 90.0;
			beta = 90.0;
			gamma = 90.0;
		}
		Crystal unitCell = new Crystal(a, b, c, alpha, beta, gamma, spacegroup);
		unitCell.setAperiodic(aperiodic);
		return unitCell;
	}

	public double[][][] gen_pot_grid(File target_file, Atom atoms[], int type){
		double pot_grid[][][] = new double[nSymm][][];//nSymm, nPoints, 4
		ArrayList<Double[]> temp_grid = new ArrayList<Double[]>();
		if(target_file != null && target_file.exists() && target_file.canRead()){
			try {
				//won't work for nSymm stuff
				for(int i = 0; i < nSymm; i++){
					read_target_file(target_file, temp_grid, type);
					pot_grid[i] = new double[temp_grid.size()][4];
					for(int j = 0; j < temp_grid.size(); j++){
						Double xyzpot[] = temp_grid.get(j);
						pot_grid[i][j][0] = xyzpot[0];
						pot_grid[i][j][1] = xyzpot[1];
						pot_grid[i][j][2] = xyzpot[2];
						pot_grid[i][j][3] = xyzpot[3];
					}
					temp_grid.clear();
				}
			} catch (IOException e) {
				e.printStackTrace();
				System.out.println("Error opening or reading target/cube file");
			}
		}
		else{
			for(int i = 0; i < nSymm; i++){
				make_grid(temp_grid, atoms);
				pot_grid[i] = new double[temp_grid.size()][4];
				for(int j = 0; j < temp_grid.size(); j++){
					Double xyz[] = temp_grid.get(j);
					pot_grid[i][j][0] = xyz[0];
					pot_grid[i][j][1] = xyz[1];
					pot_grid[i][j][2] = xyz[2];
					pot_grid[i][j][3] = pme.potpoint(xyz);
				}
				temp_grid.clear();
			}
		}
		return pot_grid;
	}

	public void read_target_file(File target_file, ArrayList<Double[]> grid, int type) throws IOException{
		Double xyz[] = new Double[3];
		String line;
		double x, y, z, pot;
		int nPoints = 0;
		if (target_file != null && target_file.exists() && target_file.canRead()) {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(target_file)));
			line = br.readLine();
			StringTokenizer st1 = new StringTokenizer(line);
			nPoints = Integer.parseInt(st1.nextToken());
			while((line = br.readLine()) != null){
				Double xyzpot[] = new Double[4];
//				System.out.println(line.trim());
//				ln = line.split(" ");
				StringTokenizer st = new StringTokenizer(line);
				st.nextToken();
				x = Double.parseDouble(st.nextToken());
				y = Double.parseDouble(st.nextToken());
				z = Double.parseDouble(st.nextToken());
				xyzpot[0] = x;
				xyzpot[1] = y;
				xyzpot[2] = z;
				if(type == 0){
					pot = Double.parseDouble(st.nextToken());
					xyzpot[3] = pot;
				}
				else if(type == 1){
					xyz[0] = x;
					xyz[1] = y;
					xyz[2] = z;
					pot = pme.potpoint(xyz);
					xyzpot[3] = pot;
				}
				grid.add(xyzpot);
			}
		}
		int atmnum;
		double rad[] = new double[nAtoms];
		for(int i = 0; i < nAtoms; i++){
			rad[i] = 1.7;
			atmnum = atoms[i].getAtomicNumber();
			if (atmnum == 0)  rad[i] = 0.00;                                                      
			if (atmnum == 1)  rad[i] = 1.20;                                                      
			if (atmnum == 2)  rad[i] = 1.40;                                                      
			if (atmnum == 6)  rad[i] = 1.70;                                                      
			if (atmnum == 7)  rad[i] = 1.55;                                                      
			if (atmnum == 8)  rad[i] = 1.52;                                                      
			if (atmnum == 9)  rad[i] = 1.47;
			if (atmnum == 10)  rad[i] = 1.54;                                                     
			if (atmnum == 14)  rad[i] = 2.10;                                                     
			if (atmnum == 15)  rad[i] = 1.80;                                                     
			if (atmnum == 16)  rad[i] = 1.80;                                                     
			if (atmnum == 17)  rad[i] = 1.75;                                                     
			if (atmnum == 18)  rad[i] = 1.88;                                                     
			if (atmnum == 35)  rad[i] = 1.85;                                                     
			if (atmnum == 36)  rad[i] = 2.02;                                                     
			if (atmnum == 53)  rad[i] = 1.98;                                                     
			if (atmnum == 54)  rad[i] = 2.16;
		}
		double big = 1000;
		double small;
		for(int i = 0; i < nPoints; i++){
			small = big;
			Double xyzpot[] = grid.get(i);
			double xi = xyzpot[0];
			double yi = xyzpot[1];
			double zi = xyzpot[2];
			for(int k = 0; k < nAtoms; k++){
				double r2 = Math.pow((xi - atoms[k].getX()),2) + Math.pow((yi - atoms[k].getY()),2) + Math.pow((zi - atoms[k].getZ()),2);
				double dist = Math.sqrt(r2) - rad[k];
				if(dist < small){
					small = dist;
					ipgrid.add(k);
				}
			}
		}
	}

	public void make_grid(ArrayList<Double[]> grid, Atom atoms[]){
		
		int nAtoms = atoms.length;
		double rad[] = new double[nAtoms];
		double rad2[] = new double[nAtoms];
		double dot[][];
		boolean keep = true;
		int nshell = 4;
		double spacing = .35;
		double density = 4 * Math.PI / Math.pow(spacing, 2);
		double roffset = 1;

		int atmnum;
		for(int i = 0; i < nAtoms; i++){
			rad[i] = 1.7;
			atmnum = atoms[i].getAtomicNumber();
			if (atmnum == 0)  rad[i] = 0.00;                                                      
			if (atmnum == 1)  rad[i] = 1.20;                                                      
			if (atmnum == 2)  rad[i] = 1.40;                                                      
			if (atmnum == 6)  rad[i] = 1.70;                                                      
			if (atmnum == 7)  rad[i] = 1.55;                                                      
			if (atmnum == 8)  rad[i] = 1.52;                                                      
			if (atmnum == 9)  rad[i] = 1.47;
			if (atmnum == 10)  rad[i] = 1.54;                                                     
			if (atmnum == 14)  rad[i] = 2.10;                                                     
			if (atmnum == 15)  rad[i] = 1.80;                                                     
			if (atmnum == 16)  rad[i] = 1.80;                                                     
			if (atmnum == 17)  rad[i] = 1.75;                                                     
			if (atmnum == 18)  rad[i] = 1.88;                                                     
			if (atmnum == 35)  rad[i] = 1.85;                                                     
			if (atmnum == 36)  rad[i] = 2.02;                                                     
			if (atmnum == 53)  rad[i] = 1.98;                                                     
			if (atmnum == 54)  rad[i] = 2.16; 
			rad[i] += roffset;
			rad2[i] = Math.pow(rad[i], 2);
		}
		for(int m = 0; m < nshell; m++){
			if(m != 0){
				for(int i = 0; i < nAtoms; i++){
					rad[i] += spacing;
					rad2[i] = Math.pow(rad[i],2);
				}
			}
			for(int i = 0; i < nAtoms; i++){
				Atom ai = atoms[i];
				double xi = ai.getX();
				double yi = ai.getY();
				double zi = ai.getZ();
				int ndot = (int) (density * rad2[i]);
				dot = sphere(ndot);
				for(int j = 0; j < ndot; j++){
					Double xyzpot[] = new Double[4];
					double xj = xi + rad[i] * dot[j][0];
					double yj = yi + rad[i] * dot[j][1];
					double zj = zi + rad[i] * dot[j][2] ;
					for(int k = 0; k <= i-1; k++){
						double xr = xj - atoms[k].getX();
						double yr = yj - atoms[k].getY();
						double zr = zj - atoms[k].getZ();
						double r2 = xr*xr + yr*yr + zr*zr;
						if (r2 < rad2[k]){
							keep = false;
						}
					}
					if(keep){
						for(int k = i+1; k < nAtoms; k++){
							double xr = xj - atoms[k].getX();
							double yr = yj - atoms[k].getY();
							double zr = zj - atoms[k].getZ();
							double r2 = xr*xr + yr*yr + zr*zr;
							if (r2 < rad2[k]){
								keep = false;
							}
						}
						if(keep){
							xyzpot[0] = xj;
							xyzpot[1] = yj;
							xyzpot[2] = zj;
							xyzpot[3] = null;
							grid.add(xyzpot);
							ipgrid.add(i);
						}
					}
					keep = true;
				}
			}
		}
                System.out.println(grid.size());
	}

	public double[][] sphere(int ndot){
		double dot[][] = new double[ndot][3];
		double tot = (double) ndot;
		double tot1 = (double) (ndot - 1);
		double h; 
		double phi;
		double phiold = 0;
		double theta;

		for(int i = 0; i < ndot; i++){
			h = -1 + 2*i/tot1;
			theta = Math.acos(h);
			if(i == 0 || i == ndot - 1){
				phi = 0;
			}
			else{
				phi = (phiold+3.6/Math.sqrt(tot * (1 - h*h))) % (2*Math.PI);
			}
			dot[i][0] = Math.sin(theta) * Math.cos(phi);
			dot[i][1] = Math.sin(theta) * Math.sin(phi);
			dot[i][2] = Math.cos(theta);
			phiold = phi;
		}

		return dot;
	}

	public void output_pgrid(File outfpot, File outfgrid, double[][][] pgrid) throws IOException{
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfpot)));
		BufferedWriter bw2 = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outfgrid)));
		DecimalFormat myFormatter = new DecimalFormat("###########0.0000000");
		DecimalFormat myFormatter2 = new DecimalFormat("###########0.0000");
		String output;
		String output2;
		for(int i = 0; i < nSymm; i++){
			if(nSymm > 1){
				System.out.println(nSymm);
				bw.write(i+"\n");
			}
			bw.write(pgrid[i].length+"\n");
			for(int j = 0; j < pgrid[i].length; j++){
				output = myFormatter.format(pgrid[i][j][0])+"\t"+myFormatter.format(pgrid[i][j][1])+"\t"+myFormatter.format(pgrid[i][j][2])+"\t"+myFormatter2.format(pgrid[i][j][3]);
				output2 = myFormatter.format(pgrid[i][j][0])+"\t"+myFormatter.format(pgrid[i][j][1])+"\t"+myFormatter.format(pgrid[i][j][2]);
				bw.write("\t"+(j+1)+"\t"+output+"\n");
				bw2.write("\t"+output2+"\n");
			}
		}
		bw.close();
		bw2.close();
		System.out.println("Potential Grid Written to File: " + outfpot.getAbsolutePath());
		System.out.println("Gaussian CUBEGEN Input Written to File: " + outfgrid.getAbsolutePath());
	}

	public void output_stats(){
		ArrayList<Integer> natm = new ArrayList<Integer>();
		ArrayList<Double> rmsa = new ArrayList<Double>();
		ArrayList<Double> patm1 = new ArrayList<Double>();
		ArrayList<Double> patm2 = new ArrayList<Double>();
		Double xyz[] = new Double[3];
		double pot;
		double target;
		double a_pot = 0;
		double a_target = 0;
		double tave = 0;
		double uave = 0;
		double avgrms = 0;
		for(int i = 0; i < nAtoms; i++){
			natm.add(0);
			rmsa.add(0.0);
			patm1.add(0.0);
			patm2.add(0.0);
		}
		int k;
		int npoints = target_grid[0].length;
		for(int i = 0; i < nSymm; i++){
			for(int j = 0; j < npoints; j++){
    			xyz[0] = target_grid[i][j][0];
    			xyz[1] = target_grid[i][j][1];
    			xyz[2] = target_grid[i][j][2];
    			pot = pme.potpoint(xyz);
    			//System.out.println(pot);
    			target = target_grid[i][j][3];
    			a_pot += Math.abs(pot);
    			a_target += Math.abs(target);
    			tave += pot - target;
    			uave += Math.abs(pot - target);
				k = ipgrid.get(j);
				natm.set(k, natm.get(k) + 1);
				patm1.set(k, patm1.get(k) + pot);
				patm2.set(k, patm2.get(k) + target);
				rmsa.set(k, rmsa.get(k) + Math.pow(pot - target, 2));
			}
		}
		for(int i = 0; i < nAtoms; i++){
			patm1.set(i, patm1.get(i) / natm.get(i));
			patm2.set(i, patm2.get(i) / natm.get(i));
			rmsa.set(i, Math.sqrt(rmsa.get(i)/natm.get(i)));
		}
		a_pot = a_pot / npoints;
		a_target = a_target / npoints;
		tave = tave / npoints;
		uave = uave / npoints;
		avgrms = avgrms();
		System.out.println("\nAverage Electrostatic Potential over Atoms :\n (Kcal/mole per unit charge)\n");
		System.out.println("Atom \t Points \t Potential \t Target \t RMS Diff\n");
		for(int i = 0; i < nAtoms; i++){
			System.out.printf("%d \t %d \t %12.6g \t %12.6g \t %12.6g\n", i, natm.get(i), patm1.get(i), patm2.get(i), rmsa.get(i));
		}
		System.out.println("\nElectrostatic Potential over all Grid Points :\n");
		System.out.printf("Average Magnitude for Potential : \t\t %12.6g\n", a_pot);
		System.out.printf("Average Magnitude for Target : \t\t\t %12.6g\n", a_target);
		System.out.printf("Average Signed Potential Difference : \t\t %12.6g\n", tave);
		System.out.printf("Average Unsigned Potential Difference : \t %12.6g\n", uave);
		System.out.printf("Root Mean Square Potential Difference : \t %12.6g\n", avgrms);
	}
	
    public double avgrms(){
    	double er = 0;
    	pme.init_prms();
    	Double xyz[] = new Double[3];
    	for(int j = 0; j < target_grid[0].length; j++){
    		xyz[0] = target_grid[0][j][0];
    		xyz[1] = target_grid[0][j][1];
    		xyz[2] = target_grid[0][j][2];
    		er += Math.pow(pme.potpoint(xyz) - target_grid[0][j][3],2);
    	}
    	er = er/target_grid[0].length;
    	return Math.sqrt(er);
    }
	
	public void store_key_file(File keyfile) throws IOException{
		if (keyfile != null && keyfile.exists() && keyfile.canRead()) {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(keyfile)));
			String line;
			while((line = br.readLine()) != null){
				key.add(line);
			}
		}
	}
	

	@Override
	public boolean optimizationUpdate(int iter, int nfun, double grms,
			double xrms, double f, double df, double angle,
			LineSearchResult info) {
		long currentTime = System.nanoTime();
		Double seconds = (currentTime - time) * 1.0e-9;
		time = currentTime;
		this.grms = grms;
		this.nSteps = iter;

		if (iter == 0) {
			logger.info("\n Limited Memory BFGS Quasi-Newton Optimization: \n\n");
			logger.info(" Cycle       Energy      G RMS    Delta E   Delta X    Angle  Evals     Time      R  Rfree\n");
		}
		if (info == null) {
			logger.info(String.format("%6d %13.4g %11.4g\n",
					iter, f, grms));
		} else {
			if (info == LineSearchResult.Success) {
				                logger.info(String.format("%6d %13.4g %11.4g %11.4g %10.4g %9.2g %7d %8.3g %6.4f\n",
				                        iter, f, grms, df, xrms, angle, nfun, seconds,
				                        avgrms()));
			} else {
				logger.info(String.format("%6d %13.4g %11.4g %11.4g %10.4g %9.2g %7d %8s\n",
						iter, f, grms, df, xrms, angle, nfun, info.toString()));
			}
		}
		if (terminate) {
			logger.info(" The optimization recieved a termination request.");
			// Tell the L-BFGS optimizer to terminate.
			return false;
		}
		return true;
	}
	
	public static void main(String args[]) throws IOException{
		//Potential2 p1 = new Potential2(4, "/users/gchattree/Research/Compounds/test_compounds/phenobarbital-test/phenobarbital.xyz", null, .8);
		//Potential2 p2 = new Potential2(3, "/users/gchattree/Research/Compounds/test_compounds/phenobarbital-test/phenobarbital.xyz", null, null);
		Potential2 p1 = new Potential2(4, "/users/gchattree/Research/Compounds/test_compounds/12-ethanediol-poltypeffx/12-ethanediol.xyz", null, .1);
		//Potential2 p2 = new Potential2(2, "/users/gchattree/Research/Compounds/test_compounds/12-ethanediol-test/12-ethanediol.xyz", null, null);
	}
}


