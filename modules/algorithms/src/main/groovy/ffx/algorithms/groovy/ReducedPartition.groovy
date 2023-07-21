package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.algorithms.optimize.TitrationManyBody
import ffx.numerics.math.DoubleMath
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.parsers.PDBFilter
import org.apache.commons.configuration2.CompositeConfiguration
import picocli.CommandLine
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard

/**
 * The ReductionPartition script performs a discrete optimization using a many-body expansion and elimination expressions.
 * <br>
 * Usage:
 * <br>
 * ffxc ManyBody [options] &lt;filename&gt;
 */
@Command(description = " Run ReducedPartition function for free energy change.", name = "ffxc ReducedPartition")
class ReducedPartition extends AlgorithmsScript {

    @Mixin
    ManyBodyOptions manyBodyOptions

    @Mixin
    AlchemicalOptions alchemicalOptions

    @CommandLine.Option(names = ["--mR", "--mutatingResidue"], paramLabel = "-1",
            description = "The residue that is mutating.")
    private int mutatingResidue = -1

    @CommandLine.Option(names = ["--resC", "--residueChain"], paramLabel = "A",
            description = "The chain that is mutating.")
    private String mutatingChain = 'A'

    @CommandLine.Option(names = ["--dC", "--distanceCutoff"], paramLabel = "-1",
            description = "Residues within the distance cutoff from the mutating residue will be optimized.")
    private double distanceCutoff = -1

    @CommandLine.Option(names = ["-n", "--residueName"], paramLabel = "ALA",
            description = "Mutant residue.")
    private String resName

    @CommandLine.Option(names = ["--rEE", "--ro-ensembleEnergy"], paramLabel = "0.0",
            description = "Keep permutations within ensemble Energy kcal/mol from the GMEC.")
    private String ensembleEnergy = "0.0"

    @CommandLine.Option(names = ["--un", "--unfolded"], paramLabel = "false",
            description = "Run the unfolded state tripeptide.")
    private boolean unfolded = false

    @CommandLine.Option(names = ["--pKa"], paramLabel = "false",
            description = "Calculating free energy change for pKa shift.")
    private boolean pKa = false

    @CommandLine.Option(names = ["--oT", "--onlyTitration"], paramLabel = "false",
            description = "Only include titratable residues in the residue selection.")
    private boolean onlyTitration = false

    @CommandLine.Option(names = ["--oP", "--onlyProtons"], paramLabel = "false",
            description = "Only allow proton movement of titratable reidues.")
    private boolean onlyProtons = false

    @CommandLine.Option(names = ["--pB", "--printBoltzmann"], paramLabel = "false",
            description = "Save the Boltzmann weights of protonated residue and total Boltzmann weights.")
    private boolean printBoltzmann = false

    @CommandLine.Option(names = ["--pF", "--printFiles"], paramLabel = "false",
            description = "Write to an energy restart file and ensemble file.")
    private boolean printFiles = false


    /**
     * An XYZ or PDB input file.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = "XYZ or PDB input file.")
    private List<String> filenames = null

    ForceFieldEnergy potentialEnergy
    TitrationManyBody titrationManyBody
    MolecularAssembly mutatedAssembly
    Binding mutatorBinding
    List<Residue> residues
    List<Residue> selectedResidues

    private String unfoldedFileName

    /**
     * ManyBody Constructor.
     */
    ReducedPartition() {
        this(new Binding())
    }

    /**
     * ManyBody Constructor.
     * @param binding The Groovy Binding to use.
     */
    ReducedPartition(Binding binding) {
        super(binding)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    ReducedPartition run() {
        if (!init()) {
            return this
        }

        double titrationPH = manyBodyOptions.getTitrationPH()
        if (manyBodyOptions.getTitration()) {
            System.setProperty("manybody-titration", "true")
        }
        System.setProperty("ro-ensembleEnergy", ensembleEnergy)
        activeAssembly = getActiveAssembly(filenames.get(0))

        //Make an unfolded state assembly when predicting folding free energy difference
        if (unfolded) {
            unfoldedFileName = "wt" + mutatingResidue + ".pdb"
            List<Atom> atoms = activeAssembly.getAtomList()
            Set<Atom> excludeAtoms = new HashSet<>()
            for (Atom atom : atoms) {
                if (atom.getResidueNumber() < mutatingResidue - 1 || atom.getResidueNumber() > mutatingResidue + 1) {
                    excludeAtoms.add(atom)
                } else if (atom.getResidueNumber() == mutatingResidue - 1 && atom.getName() == "H") {
                    excludeAtoms.add(atom)
                }
            }
            File file = new File(unfoldedFileName)
            PDBFilter pdbFilter = new PDBFilter(file, activeAssembly, activeAssembly.getForceField(),
                    activeAssembly.getProperties())
            pdbFilter.writeFile(file, false, excludeAtoms, true, true)
            setActiveAssembly(getActiveAssembly(unfoldedFileName))
        }

        String[] titratableResidues = ["HIS", "HIE", "HID", "GLU", "GLH", "ASP", "ASH", "LYS", "LYD"]
        List<String> titratableResiudesList = Arrays.asList(titratableResidues);
        double[] boltzmannWeights = new double[2]
        double[] offsets = new double[2]
        double[] titrateArray
        double[] titrateBoltzmann
        double totalBoltzmann = 0
        List<Residue> residueList = activeAssembly.getResidueList()

        List<Integer> residueNumber = new ArrayList<>()
        for (Residue residue : residueList) {
            residueNumber.add(residue.getResidueNumber())
        }


        String mutatedFileName = ""
        //Call the MutatePDB script and mutate the residue of interest
        if (filenames.size() == 1 && mutatingResidue != -1) {
            if (unfolded) {
                mutatorBinding = new Binding('-r', mutatingResidue.toString(), '-n', resName, unfoldedFileName)
            } else {
                mutatorBinding = new Binding('-r', mutatingResidue.toString(), '-n', resName, filenames.get(0), '--ch', mutatingChain)
            }

            MutatePDB mutatePDB = new MutatePDB(mutatorBinding)
            mutatePDB.run()
            mutatedFileName = mutatorBinding.getProperty('versionFileName')
        }

        String listResidues = ""
        //Select residues with alpha carbons within the distance cutoff
        if (mutatingResidue != -1 && distanceCutoff != -1) {
            double[] mutatingResCoor = new double[3]
            int index = residueNumber.indexOf(mutatingResidue)
            mutatingResCoor = residueList.get(index).getAtomByName("CA", true).getXYZ(mutatingResCoor)
            for (Residue residue: residueList) {
                double[] currentResCoor = new double[3]
                currentResCoor = residue.getAtomByName("CA", true).getXYZ(currentResCoor)
                double dist = DoubleMath.dist(mutatingResCoor, currentResCoor)
                if (dist < distanceCutoff) {
                    listResidues += "," + residue.getChainID() + residue.getResidueNumber()
                }
            }
            listResidues = listResidues.substring(1)
        }

        //Select only the titrating residues or the titrating residues and those within the distance cutoff
        if (onlyTitration || onlyProtons) {
            for (Residue residue : residueList) {
                if (titratableResiudesList.contains(residue.getName())) {
                    String titrateResNum = residue.getResidueNumber()
                    if(!listResidues.contains(titrateResNum)){
                       listResidues += "," + residue.getChainID() + titrateResNum
                    }

                    if (distanceCutoff != -1){
                        double[] titrationResCoor = new double[3]
                        //titrationResCoor = residue.getAtomByName("CA", true).getXYZ(titrationResCoor)
                        titrationResCoor = residue.getAtomList(true).last().getXYZ(titrationResCoor)
                        for (Residue residue2: residueList) {
                            double[] currentResCoor = new double[3]
                            currentResCoor = residue2.getAtomList(true).last().getXYZ(currentResCoor)
                            double dist = DoubleMath.dist(titrationResCoor, currentResCoor)
                            if (dist < distanceCutoff) {
                                String addedResidue = residue2.getChainID() + residue2.getResidueNumber()
                                if(!listResidues.contains(addedResidue)){
                                    listResidues += "," + residue2.getChainID() + addedResidue
                                }
                            }
                        }
                    }

                }

            }


            listResidues = listResidues.substring(1)
        }

        String filename = filenames.get(0)

        //Set the number of assemblies the partition function will be calculated for
        int numLoop = 1
        if (mutatingResidue != -1) {
            numLoop = 2
        }

        List<Residue> titrateResidues = new ArrayList<>()

        //Calculate all possible permutations for the number of assembles
        for (int j = 0; j < numLoop; j++) {

            // Load the MolecularAssembly second molecular assembly if applicable.
            if (j > 0) {
                if (filenames.size() == 1) {
                    mutatedAssembly = getActiveAssembly(mutatedFileName)
                    setActiveAssembly(mutatedAssembly)
                    logger.info(activeAssembly.getResidueList().toString())
                    activeAssembly.getPotentialEnergy().energy()
                    filename = mutatedFileName
                } else {
                    setActiveAssembly(getActiveAssembly(filenames.get(j)))
                    filename = filenames.get(j)
                }
            }

            if (activeAssembly == null) {
                logger.info(helpString())
                return this
            }

            CompositeConfiguration properties = activeAssembly.getProperties()

            // Application of rotamers uses side-chain atom naming from the PDB.
            if (properties.getBoolean("standardizeAtomNames", false)) {
                renameAtomsToPDBStandard(activeAssembly)
            }

            activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)
            potentialEnergy = activeAssembly.getPotentialEnergy()

            if (!pKa || onlyTitration || onlyProtons) {
                manyBodyOptions.setListResidues(listResidues)
            }


            // Collect residues to optimize.
            residues = manyBodyOptions.collectResidues(activeAssembly)
            if (residues == null || residues.isEmpty()) {
                logger.info(" There are no residues in the active system to optimize.")
                return this
            }

            // Handle rotamer optimization with titration.
            if (manyBodyOptions.getTitration()) {
                logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n")

                // Collect residue numbers.
                List<Integer> resNumberList = new ArrayList<>()
                for (Residue residue : residues) {
                    resNumberList.add(residue.getResidueNumber())
                    if (pKa) {
                        if (titratableResiudesList.contains(residue.getName())) {
                            titrateResidues.add(residue)
                        }
                    }
                }

                // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
                titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
                        resNumberList, titrationPH)
                MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly()
                setActiveAssembly(protonatedAssembly)
                potentialEnergy = protonatedAssembly.getPotentialEnergy()
            }

            //Run rotamer optimization with specified parameter
            RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly,
                    potentialEnergy, algorithmListener)
            rotamerOptimization.setPrintFiles(printFiles)
            rotamerOptimization.setWriteEnergyRestart(printFiles)
            rotamerOptimization.setOnlyProtons(onlyProtons)

            manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly)

            selectedResidues = rotamerOptimization.getResidues()

            logger.info("\n Initial Potential Energy:")
            potentialEnergy.energy(false, true)

            logger.info("\n Initial Rotamer Torsion Angles:")
            RotamerLibrary.measureRotamers(selectedResidues, false)

            // Run the optimization.
            rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(selectedResidues.size()))

            int[] currentRotamers = new int[selectedResidues.size()]

            //Keep track of the number of titrating residues
            if (pKa) {
                titrateArray = new double[selectedResidues.size()]
            }

            //Calculate possible permutations for assembly
            rotamerOptimization.checkPermutations(selectedResidues.toArray() as Residue[], 0, currentRotamers, titrateArray,
                    manyBodyOptions.getAlgorithm(selectedResidues.size()))

            //Collect the Bolztmann weights and calculated offset of each assembly
            boltzmannWeights[j] = rotamerOptimization.getTotalBoltzmann()
            offsets[j] = rotamerOptimization.getRefEnergy()

            //Calculate the fraction protonated for the titratable residues
            if (pKa) {
                titrateArray = rotamerOptimization.getFraction()
                if (printBoltzmann) {
                    titrateBoltzmann = rotamerOptimization.getTitrateBoltzmann()
                    totalBoltzmann = rotamerOptimization.getTotalBoltzmann()
                }

            }
        }

        //Print information from the fraction protonated calculations
        if (pKa) {
            int titrateCount = 0

            for (Residue residue : selectedResidues) {
                logger.info("Residue " + residue.getName() + residue.getResidueNumber() + " Fraction of Protonated: " +
                        titrateArray[titrateCount])
                if (printBoltzmann) {
                    logger.info("Residue " + residue.getName() + residue.getResidueNumber() + " Protonated Boltzmann: " +
                            titrateBoltzmann[titrateCount])
                    logger.info("Total Boltzmann: " + totalBoltzmann)
                }
                titrateCount += 1
            }
        } else {
            //Calculate Gibbs free energy change of mutating residues
            double gibbs = -(0.6) * (Math.log(boltzmannWeights[1] / boltzmannWeights[0]))
            logger.info("\n Gibbs Free Energy Change: " + gibbs)
        }


        return this
    }

    /**
     * Returns the potential energy of the active assembly. Used during testing assertions.
     * @return potentialEnergy Potential energy of the active assembly.
     */
    ForceFieldEnergy getPotential() {
        return potentialEnergy
    }


}
