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

import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import java.util.logging.Level

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
            atoms[i].setActive(false)
            atoms[i].setUse(false)
        }

        // For now, always use the original coordinates as a (fixed) rotamer.
        boolean useOriginalRotamers = true

        // TODO: handle getting AA vs. NA libraries.
        RotamerLibrary rotamerLibrary = new RotamerLibrary(ProteinLibrary.Richardson, useOriginalRotamers)

        // Get the residue list.
        ArrayList<Residue> residues = activeAssembly.getResidueList()

        logger.info("Number of residues: "+residues.size().toString()+"\n")
        // Loop over Residues and set backbone atoms to being used.
        for (Residue residue : residues) {
            ArrayList<Atom> backBoneAtoms = residue.getBackboneAtoms()
            for (Atom atom : backBoneAtoms) {
                atom.setUse(true)
            }
        }

        // Define writers for .rot file
        FileWriter fw
        BufferedWriter bw

        // Create .rot file name: should match input file name and end in ".rot"
        int lenOfFile = filename.length()
        int stop = lenOfFile-4

        String rotFileName = filename.substring(0,stop)
        rotFileName = rotFileName.concat(".rot")

        try{
            fw = new FileWriter(rotFileName)
            bw = new BufferedWriter(fw)

        } catch (Exception e){
            String message = "Exception creating file writers"
            logger.log(Level.WARNING, message, e)
            throw e
        }

        // Loop over Residues
        for (Residue residue : residues) {

            // Get this residue's rotamers.
            Rotamer[] rotamers = residue.getRotamers(rotamerLibrary)

            // If there's not at least one rotamer beyond the original coordinates, continue.
            if (rotamers == null || rotamers.length < 2) {
                continue
            }

            // Configure "active" and "use" flags.
            ArrayList<Atom> sideChainAtoms = residue.getSideChainAtoms()
            for (Atom atom : sideChainAtoms) {
                atom.setActive(true)
                atom.setUse(true)
            }

            // Loop over rotamers for this Residue.
            for (int i=1; i<rotamers.length; i++) {
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
            RotamerLibrary.applyRotamer(residue, rotamers[0])

            // Revert the active and use flags.
            for (Atom atom : sideChainAtoms) {
                atom.setActive(false)
                atom.setUse(false)
            }
        }

        try {
            bw.flush()
            fw.flush()
            fw.close()
            //bw.flush()
            //bw.close()
        }catch (Exception e){
            String message = "Exception closing file writers"
            logger.log(Level.WARNING, message, e)
            throw e
        }

        return this
    }
}