// Force Field X Imports
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

// input PDB file
String pdbfilename = args[0];
if (pdbfilename == null){
  pdbfilename = "examples/1N7S.pdb";
}

// input data filename string (optional - if not given, data must be present as pdbfilename.[mtz/cif/ent/cns]
String xrayfilename = args[1];

String neutronfilename = args[2];

// output following scaling, sigmaA and likelihood calculation
// String outputfile = "examples/1N7S_ffx.mtz";
String outputfile = null;

// Things below this line normally do not need to be changed.
// ===============================================================================================

systems = open(pdbfilename);

DiffractionFile xrayfile = null;
if (xrayfilename != null) {
  xrayfile = new DiffractionFile(xrayfilename, 1.0, false);
} else {
  xrayfile = new DiffractionFile(systems, 1.0, false);
}

DiffractionFile neutronfile = null;
if (neutronfilename != null) {
  neutronfile = new DiffractionFile(neutronfilename, 1.0, true);
} else {
  neutronfile = new DiffractionFile(systems, 1.0, true);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, xrayfile, neutronfile);

diffractiondata.scalebulkfit();
diffractiondata.printstats();
energy();

if (outputfile != null){
  diffractiondata.writedata(outputfile);
}

// diffractiondata.timings();
