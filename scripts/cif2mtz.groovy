// Apache Imports
import org.apache.commons.io.FilenameUtils;

// Force Field X Imports
import ffx.crystal.Crystal;
import ffx.crystal.ReflectionList;
import ffx.crystal.Resolution;
import ffx.xray.CIFFilter;
import ffx.xray.DiffractionRefinementData;
import ffx.xray.MTZWriter;

// Name of the PDB with crystal header information
String modelfilename = args[0];

// input CIF file
String datafilename = args[1];

// Things below this line normally do not need to be changed.
// ===============================================================================================

if (modelfilename == null
    || datafilename == null){
  println("\n Usage: ffxc cif2mtz PDBfilename CIFfilename");
  return;
}

systems = open(modelfilename);

CIFFilter ciffilter = new CIFFilter();
ReflectionList reflectionlist = ciffilter.getReflectionList(new File(datafilename), systems[0].getProperties());

if (reflectionlist == null) {
  println("Using crystal information from PDB to generate MTZ file");

  Crystal crystal = systems[0].getCrystal().getUnitCell();
  double res = ciffilter.getResolution(new File(datafilename), crystal);
  if (res < 0.0) {
    println("resolution could not be determined from PDB and CIF file");
    return;
  }

  Resolution resolution = new Resolution(res);
  reflectionlist = new ReflectionList(crystal, resolution, systems[0].getProperties());
}

DiffractionRefinementData refinementdata = new DiffractionRefinementData(systems[0].getProperties(), reflectionlist);
ciffilter.readFile(new File(datafilename), reflectionlist, refinementdata, systems[0].getProperties());

MTZWriter mtzwriter = new MTZWriter(reflectionlist, refinementdata, FilenameUtils.removeExtension(datafilename) + "_cif.mtz", true);
mtzwriter.write();
