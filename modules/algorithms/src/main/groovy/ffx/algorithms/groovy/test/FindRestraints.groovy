package ffx.algorithms.groovy.test

import ffx.algorithms.cli.AlgorithmsScript
import ffx.numerics.Potential
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Molecule
import picocli.CommandLine
import static java.lang.String.format

class FindRestraints extends AlgorithmsScript{
    /**
     * --hostName Molecule name of the host in the file.
     */
    @CommandLine.Option(names = ['--hostName'], paramLabel = 'None',
            description = 'Host molecule name in the file.')
    String hostName = "BCD"

    /**
     * --guestName Molecule name of the guest in the file.
     */
    @CommandLine.Option(names = ['--guestName'], paramLabel = 'LIG',
            description = 'Ligand molecule name in the file.')
    String guestName = "LIG"

    /**
     * --distanceCutoff Cutoff to use when selecting guest atoms near host COM.
     */
    @CommandLine.Option(names = ['--distanceCutoff'], paramLabel = '5',
            description = 'Cutoff to use when selecting guest atoms near host COM')
    double distanceCutoff = 5

    /**
     * One or more filenames.
     */
    @CommandLine.Parameters(arity = "1..*", paramLabel = "files",
            description = "XYZ or PDB input files.")
    private List<String> filenames

    /**
     * Creation of a public field to try and make the JUnit test work, original code does not declare this as a public field.
     * Originally it is declared in the run method
     */
    public Potential potential
    Potential getPotentialObject() {
        return potential
    }

    /**
     * Dynamics Constructor.
     */
    FindRestraints() {
        this(new Binding())
    }

    /**
     * Dynamics Constructor.
     * @param binding The Groovy Binding to use.
     */
    FindRestraints(Binding binding) {
        super(binding)
    }

    FindRestraints run() {
        if (!init()) {
            return this
        }

        activeAssembly = getActiveAssembly(filenames[0])
        if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        Molecule[] molArr = activeAssembly.getMoleculeArray()

        List<Atom> restrainList = new ArrayList<Atom>()
        double[] COM = new double[3]
        for (Molecule molecule: molArr) {
            logger.info(format("Molecule name: "+ molecule.getName()))
            if (molecule.getName().contains(hostName)) {
                Atom[] host_atoms = molecule.getAtomList()
                COM = getCOM(host_atoms)
                logger.info(format("Center of mass of host molecule: "+ COM))
            }
            else if (molecule.getName().contains(guestName)) {
                Atom[] guest_atoms = molecule.getAtomList()
                for (Atom atom : guest_atoms) {
                    double dist = Math.sqrt(Math.pow(atom.getXYZ().get()[0] - COM[0], 2) +
                            Math.pow(atom.getXYZ().get()[1] - COM[1], 2) +
                            Math.pow(atom.getXYZ().get()[2] - COM[2], 2))
                    logger.info(format("Atom: "+ atom))
                    logger.info(format("XYZ: "+ atom.getXYZ().get()))
                    logger.info(format("Distance from host COM: "+ dist))

                    if (dist < distanceCutoff && atom.isHeavy()) {
                        restrainList.add(atom)
                    }
                }
            }
        }
        logger.info(format("Number of atoms to restrain: "+ restrainList.size()))
        int[] restrainIndices = new int[restrainList.size()]
        restrainIndices = restrainList.collect { it.getIndex() }
        logger.info(format( "Restrain list indices: "+ restrainIndices))
    }

    /**
     * Gets the center of mass of a set of atoms
     * @param atoms
     * @return x,y,z coordinates of center of mass
     */
    private static double[] getCOM(Atom[] atoms){
        // Get center of mass of moleculeOneAtoms
        double[] COM = new double[3];
        double totalMass = 0.0;
        for(Atom s: atoms){
            double[] pos = s.getXYZ().get();
            COM[0] += pos[0] * s.getMass();
            COM[1] += pos[1] * s.getMass();
            COM[2] += pos[2] * s.getMass();
            totalMass += s.getMass();
        }
        totalMass = 1 / totalMass;
        COM[0] *= totalMass;
        COM[1] *= totalMass;
        COM[2] *= totalMass;

        return COM;
    }
}
