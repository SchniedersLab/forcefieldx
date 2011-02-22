package ffx.autoparm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;

import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;

/**
 *  Poledit Provides Multipole Parameters from GDMA Output
 */

public class Poledit {
	
	private double localmultipole[][];
	private double globalmultipole[][];
	private ArrayList<Atom> atoms = new ArrayList<Atom>();
	private ArrayList<Double> polarity = new ArrayList<Double>();
	private ArrayList<Double> pdamp = new ArrayList<Double>();
	public static final double BOHR = 0.52917720859;

	
	/**
	 * This method takes data from GDMA and prints out multipole parameters
	 * 
	 * @param gdmaoutfname
	 * File location of multipole params output by GDMA
	 * 
	 * @param peditinfname
	 * File location of molecular polarization group information
	 */
	public Poledit(String gdmaoutfname, String peditinfname){
		readGDMA(gdmaoutfname);
		
		int index = gdmaoutfname.lastIndexOf(".");
		String name = gdmaoutfname.substring(0, index);
		
		setup_print_xyz(name);
		
		
	}
	
	public void readGDMA(String gdmaoutfname){
		ArrayList<double[]> sphmultipole = new ArrayList<double[]>();
		File gdmaoutf = new File(gdmaoutfname);
		try {
			if (gdmaoutf != null && gdmaoutf.exists() && gdmaoutf.canRead()) {
				BufferedReader br = new BufferedReader(new InputStreamReader(
						new FileInputStream(gdmaoutf)));
				String line;
				int i = 1;
				String element;
				while ((line = br.readLine()) != null) {
					if(line.contains("x =") && !line.contains("origin")){
						element = line.split(" +")[0];
						double xyz[] = new double[3];
						xyz[0] = Double.parseDouble(line.split(" +")[3]);
						xyz[1] = Double.parseDouble(line.split(" +")[6]);                       
						xyz[2] = Double.parseDouble(line.split(" +")[9]);
						line = br.readLine();
						//radius = Double.parseDouble(line.split(" +")[8]);
						line = br.readLine();
						double mp[] = new double[9];
						mp[0] = Double.parseDouble(line.split(" +")[3]); //Q00
						line = br.readLine();
						mp[1] = Double.parseDouble(line.split(" +")[5]); //Q10
						mp[2] = Double.parseDouble(line.split(" +")[8]); //Q11c
						mp[3] = Double.parseDouble(line.split(" +")[11]); //Q11s
						line = br.readLine();
						mp[4] = Double.parseDouble(line.split(" +")[5]); //Q20
						mp[5] = Double.parseDouble(line.split(" +")[8]); //Q21c
						mp[6] = Double.parseDouble(line.split(" +")[11]); //Q21s
						line = br.readLine();
						mp[7] = Double.parseDouble(line.split(" +")[3]); //Q22c
						mp[8] = Double.parseDouble(line.split(" +")[6]); //Q22s
						
						double radius[] = {0};
						double polart[] = {0};
						double pd[] = {0};
						AtomType at = get_atom_type(i, element, radius, polart, pd);

						Atom a = new Atom(Integer.toString(i));
						a.setXYZ(xyz);
						a.setAtomType(at);
						a.setBornRadius(radius[0]);
						polarity.add(polart[0]);
						pdamp.add(pd[0]);
						atoms.add(a);
						sphmultipole.add(mp);
						i++;
					}
				}
			}
		}
		catch(Exception e){
			e.printStackTrace();
			System.out.println("Error Reading File");
			System.exit(1);
		}
		
		//Convert multipoles from spherical harmonics to cartesian coordinates
		//should this be global or local?//
		globalmultipole = new double[atoms.size()][13];
		double term = Math.sqrt(.75);
		for(int i = 0; i < atoms.size(); i++)
		{
			globalmultipole[i][0] = sphmultipole.get(i)[0]; //Q00
			globalmultipole[i][1] = sphmultipole.get(i)[1]; //Q10
			globalmultipole[i][2] = sphmultipole.get(i)[2]; //Q11c
			globalmultipole[i][3] = sphmultipole.get(i)[3]; //Q11s
			globalmultipole[i][4] = -.5 * sphmultipole.get(i)[4] + term * sphmultipole.get(i)[7]; //-.5 * Q20 + term * Q22c
			globalmultipole[i][5] = term * sphmultipole.get(i)[8]; //term * Q22s
			globalmultipole[i][6] = term * sphmultipole.get(i)[5]; //term * Q21c
			globalmultipole[i][7] = term * sphmultipole.get(i)[8]; //term * Q22s (same as 5)
			globalmultipole[i][8] = -.5 * sphmultipole.get(i)[4] - term * sphmultipole.get(i)[7]; //-.5 * Q20 - term * Q22c
			globalmultipole[i][9] = term * sphmultipole.get(i)[6]; // term * Q21s
			globalmultipole[i][10] = term * sphmultipole.get(i)[5]; //term * Q21c (same as 6)
			globalmultipole[i][9] = term * sphmultipole.get(i)[6]; // term * Q21s (same as 9)
			globalmultipole[i][12] = sphmultipole.get(i)[4]; //Q20
		}
		
		//Maintain traceless values
		for(int i = 0; i < atoms.size(); i++){
			for(int j = 1; j < 4; j++){
				globalmultipole[i][j] = globalmultipole[i][j] * BOHR;
			}
			for(int j = 4; j < 13; j++){
				globalmultipole[i][j] = globalmultipole[i][j] * Math.pow(BOHR,2)/3;
			}
		}
	}
	
	/**
	 * Creates and returns an atomType object from an element name
	 * 
	 * To work on:
	 * What is 'thole' and is it necessary in this function.
	 * 
	 * @param type
	 * @param element
	 * @param radius
	 * @param polarity
	 * @param pdamp
	 * @return atomType
	 */
	
	public AtomType get_atom_type(int type, String element, double[] radius,
			double[] polarity, double[] pdamp) {
		int atmnum = 0;
		double atmmass;
		double rad, pol, pd;
		
		// Get atomic number from element name
		if (element.equalsIgnoreCase("SI")) {
			atmnum = 14;
		} else if (element.equalsIgnoreCase("CL")) {
			atmnum = 17;
		} else if (element.equalsIgnoreCase("BR")) {
			atmnum = 35;
		} else if (element.equalsIgnoreCase("H")) {
			atmnum = 1;
		} else if (element.equalsIgnoreCase("B")) {
			atmnum = 5;
		} else if (element.equalsIgnoreCase("C")) {
			atmnum = 6;
		} else if (element.equalsIgnoreCase("N")) {
			atmnum = 7;
		} else if (element.equalsIgnoreCase("O")) {
			atmnum = 8;
		} else if (element.equalsIgnoreCase("F")) {
			atmnum = 9;
		} else if (element.equalsIgnoreCase("P")) {
			atmnum = 15;
		} else if (element.equalsIgnoreCase("S")) {
			atmnum = 16;
		} else if (element.equalsIgnoreCase("I")) {
			atmnum = 53;
		}
		// Set atomic radius
		radius[0] = 0.77;
		if (atmnum == 0)
			radius[0] = 0.00;
		if (atmnum == 1)
			radius[0] = 0.37;
		if (atmnum == 2)
			radius[0] = 0.32;
		if (atmnum == 6)
			radius[0] = 0.77;
		if (atmnum == 7)
			radius[0] = 0.75;
		if (atmnum == 8)
			radius[0] = 0.73;
		if (atmnum == 9)
			radius[0] = 0.71;
		if (atmnum == 10)
			radius[0] = 0.69;
		if (atmnum == 14)
			radius[0] = 1.11;
		if (atmnum == 15)
			radius[0] = 1.06;
		if (atmnum == 16)
			radius[0] = 1.02;
		if (atmnum == 17)
			radius[0] = 0.99;
		if (atmnum == 18)
			radius[0] = 0.97;
		if (atmnum == 35)
			radius[0] = 1.14;
		if (atmnum == 36)
			radius[0] = 1.10;
		if (atmnum == 53)
			radius[0] = 1.33;
		if (atmnum == 54)
			radius[0] = 1.30;
		radius[0] = 1.1 * radius[0];
		
		// Set polarity
		atmmass = 1;
		polarity[0] = 0.0;
		if (atmnum == 1) {
			atmmass = 1.008;
			polarity[0] = 0.496;
		} else if (atmnum == 5) {
			atmmass = 10.810;
			polarity[0] = 1.600;
		} else if (atmnum == 6) {
			atmmass = 12.011;
			polarity[0] = 1.334;
		} else if (atmnum == 7) {
			atmmass = 14.007;
			polarity[0] = 1.073;
		} else if (atmnum == 8) {
			atmmass = 15.999;
			polarity[0] = 0.837;
		} else if (atmnum == 9) {
			atmmass = 18.998;
		} else if (atmnum == 14) {
			atmmass = 28.086;
		} else if (atmnum == 15) {
			atmmass = 30.974;
			polarity[0] = 1.828;
		} else if (atmnum == 16) {
			atmmass = 32.066;
			polarity[0] = 2.800;
		} else if (atmnum == 17) {
			atmmass = 35.453;
			polarity[0] = 4.000;
		} else if (atmnum == 35) {
			atmmass = 79.904;
			polarity[0] = 5.650;
		} else if (atmnum == 53) {
			atmmass = 126.904;
			polarity[0] = 7.250;
		}
        
		pdamp[0] = Math.pow(polarity[0], 1/6);
		
		AtomType at = new AtomType(type, type, element, null, atmnum, atmmass, 0);
		return at;
	}
	
	
	/**
	 * Sets up connectivities based on radii and then prints this information out to an xyz file
	 */
	public void setup_print_xyz(String name){
		//Set up connectivities based on radii
	    for(int i = 0; i < atoms.size() - 1; i++){
	    	Atom ai = atoms.get(i);
	    	for(int j = i+1; j < atoms.size(); j++){
		    	Atom aj = atoms.get(j);
	    		double xr = aj.getX() - ai.getX();
	    		double yr = aj.getY() - ai.getY(); 
	    		double zr = aj.getZ() - ai.getZ();
	    		double rij = ai.getBornRadius() + aj.getBornRadius();
	    		double dij = Math.sqrt(xr * xr + yr * yr + zr * zr);
	    		if(dij < rij){;
	    			Bond b = new Bond(ai,aj);
	    		}
	    	}
	    }
		
		File outf = new File(name + ".xyz");
		try {
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outf)));
			DecimalFormat myFormatter = new DecimalFormat(" ##########0.00000;-##########0.00000");
			
			bw.write(String.format("\t%5d\n",atoms.size()));
			for(Atom a : atoms){
				String output = String.format("%5d",a.getAtomType().type) + "  " + a.getAtomType().name + "\t" + myFormatter.format(a.getX()) + "   " + myFormatter.format(a.getY()) + "   " + myFormatter.format(a.getZ()) + "" + String.format("\t%-5d",a.getAtomType().atomClass);
				for(int i = 0; i < a.getBonds().size(); i++){
					output += "\t" + String.format("%-1d",a.getBonds().get(i).get1_2(a).getAtomType().type);
				}
				bw.write(output+"\n");
			}
			bw.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String args[]){
		Poledit p = new Poledit("/users/gchattree/Research/Compounds/test_compounds/12-ethanediol-test/12-ethanediol.gdmaout","/users/gchattree/Research/Compounds/test_compounds/12-ethanediol-test/12-ethanediol-peditin.txt");
		//Poledit p2 = new Poledit("/users/gchattree/Research/Compounds/test_compounds/phenobarbital-test/phenobarbital.gdmaout","/users/gchattree/Research/Compounds/test_compounds/phenobarbital-test/phenobarbital-peditin.txt");

	}

}
