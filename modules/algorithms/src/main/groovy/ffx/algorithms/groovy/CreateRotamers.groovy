package ffx.algorithms.groovy

import ffx.algorithms.Minimize
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.MinimizeOptions
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.bonded.RotamerLibrary.ProteinLibrary
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level
import java.util.stream.Collectors

/**
 * The CreateRotamers script creates a set of conformation dependent rotamers.
 * <br>
 * Usage:
 * <br>
 * ffxc CreateRotamers [options] &lt;filename&gt;
 */
@Command(description = " Creates a set of conformation dependent rotamers.", name = "ffxc CreateRotamers")
class CreateRotamers extends AlgorithmsScript {

    @Mixin
    MinimizeOptions minimizeOptions

    // TODO: instead @Mixin a subset of current ManyBodyOptions.

    /**
     * -L or --library Choose either Ponder and Richards (1) or Richardson (2)
     * rotamer library.
     */
    @Option(names = ["-L" , "--library"], paramLabel = "2",
            description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
    int library = 2

    /**
     * The final argument should be a filename.
     */
    @Parameters(arity = "1..*", paramLabel = "files", description = 'Atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    @Override
    CreateRotamers run() {

        if (!init()) {
            return this
        }

        if (filenames != null && filenames.size() > 0) {
            MolecularAssembly[] assemblies = algorithmFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return
        }

        String filename = activeAssembly.getFile().getAbsolutePath()
        logger.info(" Running CreateRotamers on " + filename)

        Atom[] atoms = activeAssembly.getAtomArray()
        int nAtoms = atoms.length

        // Set all atoms to be "inactive".
        for (int i = 0; i < nAtoms; i++) {
            atoms[i].setActive(false);
        }

        // For now, always use the original coordinates as a (fixed) rotamer.
        boolean useOriginalRotamers = true

        // TODO: handle getting AA vs. NA libraries.
        RotamerLibrary rotamerLibrary = new RotamerLibrary(ProteinLibrary.intToProteinLibrary(library), useOriginalRotamers)

        // Get the residue list.
        List<Residue> residues = activeAssembly.getResidueList();/*.stream().
                filter({ Residue r -> Rotamer[] rots = r.getRotamers(rotamerLibrary); return rots != null && rots.length > 1;}).
                collect(Collectors.toList());*/

        logger.info(String.format(" Number of residues: %d\n", residues.size()));

        // Loop over Residues and set sidechain atoms to not be used.
        for (Residue residue : residues) {
            for (Atom atom : residue.getVariableAtoms()) {
                atom.setUse(false);
            }
        }

        // Create .rot file name: should match input file name and end in ".rot"
        String rotFileName = String.format("%s.rot", FilenameUtils.removeExtension(filename));

        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(new File(rotFileName)));

            // Loop over Residues
            for (Residue residue : residues) {

                // Get this residue's rotamers.
                Rotamer[] rotamers = residue.getRotamers(rotamerLibrary);

                if (rotamers == null || rotamers.length < 2) {
                    continue;
                }

                // Configure "active" and "use" flags.
                List<Atom> sideChainAtoms = residue.getSideChainAtoms()
                for (Atom atom : sideChainAtoms) {
                    atom.setActive(true)
                    atom.setUse(true)
                }

                // Loop over rotamers for this Residue.
                for (int i=1; i< rotamers.length; i++) {
                    Rotamer rotamer = rotamers[i]

                    // Apply the rotamer (i.e. amino acid side-chain or nucleic acid suite).
                    RotamerLibrary.applyRotamer(residue, rotamer)

                    // Locally minimize.
                    Minimize minimize = new Minimize(activeAssembly, activeAssembly.getPotentialEnergy(), algorithmListener)
                    minimize.minimize(minimizeOptions.getEps(), minimizeOptions.getIterations())

                    // Save out coordinates to a rotamer file.
                    for (Atom atom : sideChainAtoms) {
                        double x = atom.getX()
                        double y = atom.getY()
                        double z = atom.getZ()
                        logger.info(String.format(" %s %16.8f %16.8f %16.8f", atom.toString(), x, y, z));

                        try{
                            /*chainID, chainName (segID),resName, resNum, rotamerNum, atomName, coordinates <x,y,z>*/
                            bw.write("[GLOBAL]:1:<residue>:"+residue.chainID+":"+residue.segID+":"+residue.name+":"+residue.getResidueNumber()+
                                    ":<rotamer>:"+rotamer.name+":"+i+":"+atom.name+
                                    ":<coordinates>:"+x+":"+y+":"+z+"\n")
                        } catch (Exception e){
                            String message = "Exception writing to file"
                            logger.log(Level.WARNING, message, e)
                            throw e
                        }

                    }
                }

                // Set the Residue conformation back to rotamer 0.
                RotamerLibrary.applyRotamer(residue, rotamers[0]);

                // Revert the active and use flags.
                for (Atom atom : sideChainAtoms) {
                    atom.setActive(false)
                    atom.setUse(false)
                }
            }

        } finally {
            bw?.flush();
            bw?.close();
        }

        return this;
    }
}