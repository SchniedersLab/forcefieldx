// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.xray.DiffractionData;
import ffx.xray.DiffractionFile;
import ffx.xray.CrystalReciprocalSpace.SolventModel;

// Name of the file (PDB or XYZ).
String modelfilename = args[0];

// input data filename string (optional - if not given, data must be present as modelfilename.[mtz/cif/ent/cns]
String datafilename = args[1];

// Things below this line normally do not need to be changed.
// ===============================================================================================

if (modelfilename == null){
  modelfilename = "examples/1N7S.pdb";
}
systems = open(modelfilename);

DiffractionFile diffractionfile = null;
if (datafilename != null) {
  diffractionfile = new DiffractionFile(datafilename, 1.0, false);
} else {
  diffractionfile = new DiffractionFile(systems, 1.0, false);
}

DiffractionData diffractiondata = new DiffractionData(systems, systems[0].getProperties(), SolventModel.POLYNOMIAL, diffractionfile);

diffractiondata.scalebulkfit();
diffractiondata.printstats();
energy();

diffractiondata.writedata(FilenameUtils.removeExtension(modelfilename) + "_ffx.mtz");

// diffractiondata.timings();
