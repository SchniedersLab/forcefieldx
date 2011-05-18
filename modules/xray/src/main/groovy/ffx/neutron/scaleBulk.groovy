// Force Field X Imports
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

// input PDB file
String modelfilename = args[0];

// input data filename string (optional - if not given, data must be present as pdbfilename.[mtz/cif/ent/cns]
String xrayfilename = args[1];

String neutronfilename = args[2];

// Things below this line normally do not need to be changed.
// ===============================================================================================

if (modelfilename == null){
   modelfilename = "examples/1N7S.pdb";
}
systems = open(modelfilename);

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

diffractiondata.scaleBulkFit();
diffractiondata.printStats();
energy();

diffractiondata.writeData(FilenameUtils.removeExtension(modelfilename) + "_ffx.mtz");

// diffractiondata.timings();
