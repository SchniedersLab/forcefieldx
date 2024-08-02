package ffx.realspace.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.algorithms.optimize.TitrationManyBody
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.LambdaInterface
import ffx.potential.bonded.Residue
import ffx.potential.bonded.Rotamer
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.parameters.ForceField
import ffx.potential.parsers.PDBFilter
import ffx.realspace.cli.RealSpaceOptions
import ffx.xray.RefinementEnergy
import org.apache.commons.configuration2.CompositeConfiguration
import picocli.CommandLine

import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard
import static java.lang.String.format

@CommandLine.Command(description = " Run GenZ function for free energy change.", name = "ffxc realspace.GenZ")
class GenZ extends AlgorithmsScript  {

    @CommandLine.Mixin
    ManyBodyOptions manyBodyOptions

    @CommandLine.Mixin
    AlchemicalOptions alchemicalOptions

    @CommandLine.Mixin
    private RealSpaceOptions realSpaceOptions

    @CommandLine.Option(names = ["--rEE", "--ro-ensembleEnergy"], paramLabel = "0.0",
            description = "Keep permutations within ensemble Energy kcal/mol from the GMEC.")
    private String ensembleEnergy = "0.0"

    @CommandLine.Option(names = ["--pF", "--printFiles"], paramLabel = "false",
            description = "Write to an energy restart file and ensemble file.")
    private boolean printFiles = false

    @CommandLine.Option(names = ["--pKa"], paramLabel = "false",
            description = "Calculating protonation populations for pKa shift.")
    private boolean pKa = false

    /**
     * One or more filenames.
     */
    @CommandLine.Parameters(arity = "1..*", paramLabel = "files", description = "PDB and Real Space input files.")
    private List<String> filenames
    private RefinementEnergy refinementEnergy

    ForceFieldEnergy potentialEnergy
    private MolecularAssembly[] molecularAssemblies
    private ForceField forceField
    TitrationManyBody titrationManyBody
    List<Residue> selectedResidues
    MolecularAssembly[] conformerAssemblies = new MolecularAssembly[3]

    /**
     * ManyBody Constructor.
     */
    GenZ() {
        this(new Binding())
    }

    /**
     * ManyBody Constructor.
     * @param binding The Groovy Binding to use.
     */
    GenZ(Binding binding) {
        super(binding)
    }

    /**
     * {@inheritDoc}
     */
    @Override
    GenZ run() {
        if (!init()) {
            return this
        }


        double titrationPH = manyBodyOptions.getTitrationPH()
        double inclusionCutoff = manyBodyOptions.getInclusionCutoff()
        int mutatingResidue = manyBodyOptions.getInterestedResidue()
        boolean onlyProtons = manyBodyOptions.getOnlyProtons()
        boolean onlyTitration = manyBodyOptions.getOnlyTitration()
        double pHRestraint = manyBodyOptions.getPHRestraint()
        if (manyBodyOptions.getTitration()) {
            System.setProperty("manybody-titration", "true")
        }


        // Many-Body expansion of the X-ray target converges much more quickly with the NEA.
        String nea = System.getProperty("native-environment-approximation", "true")
        System.setProperty("native-environment-approximation", nea)

        boolean lambdaTerm = alchemicalOptions.hasSoftcore()
        if (lambdaTerm) {
            // Turn on softcore van der Waals
            System.setProperty("lambdaterm", "true")
            // Turn of alchemical electrostatics
            System.setProperty("elec-lambdaterm", "false")
            // Turn on intra-molecular softcore
            System.setProperty("intramolecular-softcore", "true");
        }
        System.setProperty("ro-ensembleEnergy", ensembleEnergy)

        String modelFilename
        if (filenames != null && filenames.size() > 0) {
            molecularAssemblies = algorithmFunctions.openAll(filenames.get(0))
            activeAssembly = molecularAssemblies[0]
            modelFilename = filenames.get(0)
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        } else {
            molecularAssemblies = [activeAssembly]
            modelFilename = activeAssembly.getFile().getAbsolutePath()
        }

        CompositeConfiguration properties = activeAssembly.getProperties()
        activeAssembly.getPotentialEnergy().setPrintOnFailure(false, false)
        potentialEnergy = activeAssembly.getPotentialEnergy()
        forceField = activeAssembly.getForceField()
        if(forceField == null){
            logger.info("This force field is null")
        }

        String filename = filenames.get(0)

        List<Residue> titrateResidues = new ArrayList<>()

        //Prepare variables for saving out the highest population rotamers (optimal rotamers)
        Set<Atom> excludeAtoms = new HashSet<>()
        boolean isTitrating = false


        String[] titratableResidues = ["HIS", "HIE", "HID", "GLU", "GLH", "ASP", "ASH", "LYS", "LYD", "CYS", "CYD"]
        List<String> titratableResiudesList = Arrays.asList(titratableResidues);
        double[] boltzmannWeights = new double[1]
        double[] offsets = new double[1]
        double[][] populationArray = new double[1][55]
        double[][] titrateBoltzmann
        double[] protonationBoltzmannSums
        double totalBoltzmann = 0
        List<Residue> residueList = activeAssembly.getResidueList()

        // Collect residues to optimize.
        List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly)
        if (residues == null || residues.isEmpty()) {
            logger.info(" There are no residues in the active system to optimize.")
            return this
        }

        // Application of rotamers uses side-chain atom naming from the PDB.
        if (properties.getBoolean("standardizeAtomNames", false)) {
            renameAtomsToPDBStandard(activeAssembly)
        }

        // Handle rotamer optimization with titration.
        if (manyBodyOptions.getTitration()) {
            logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n")

            // Collect residue numbers.
            List<Integer> resNumberList = new ArrayList<>()
            for (Residue residue : residues) {
                resNumberList.add(residue.getResidueNumber())
            }

            // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
            titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
                    resNumberList, titrationPH, manyBodyOptions)
            MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly()
            setActiveAssembly(protonatedAssembly)
            potentialEnergy = protonatedAssembly.getPotentialEnergy()
        }
        molecularAssemblies = [activeAssembly]
        refinementEnergy = realSpaceOptions.toRealSpaceEnergy(filenames, molecularAssemblies)


        if (lambdaTerm) {
            alchemicalOptions.setFirstSystemAlchemistry(activeAssembly)
            LambdaInterface lambdaInterface = (LambdaInterface) potentialEnergy
            double lambda = alchemicalOptions.getInitialLambda()
            logger.info(format(" Setting ManyBody softcore lambda to: %5.3f", lambda))
            lambdaInterface.setLambda(lambda)
        }

        RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly, refinementEnergy, algorithmListener)
        rotamerOptimization.setPrintFiles(printFiles)
        rotamerOptimization.setWriteEnergyRestart(printFiles)
        rotamerOptimization.setPHRestraint(pHRestraint)
        rotamerOptimization.setOnlyProtons(onlyProtons)
        rotamerOptimization.setpH(titrationPH)


        manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly)

        double[] x = new double[refinementEnergy.getNumberOfVariables()]
        x = refinementEnergy.getCoordinates(x)
        double e = refinementEnergy.energy(x, true)
        logger.info(format("\n Initial target energy: %16.8f ", e))

        selectedResidues = rotamerOptimization.getResidues()
        rotamerOptimization.initFraction(selectedResidues)

        RotamerLibrary.measureRotamers(residueList, false)

        rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residueList.size()))



        int[] currentRotamers = new int[selectedResidues.size()]

        //Calculate possible permutations for assembly
        rotamerOptimization.getFractions(selectedResidues.toArray() as Residue[], 0, currentRotamers)
        if(pKa){
            rotamerOptimization.getProtonationPopulations(selectedResidues.toArray() as Residue[])
        }

        //Collect the Bolztmann weights and calculated offset of each assembly
        boltzmannWeights = rotamerOptimization.getTotalBoltzmann()
        offsets = rotamerOptimization.getRefEnergy()

        //Calculate the populations for the all residue rotamers
        populationArray = rotamerOptimization.getFraction()

        FileWriter fileWriter = new FileWriter("populations.txt")
        int residueIndex = 0
        for (Residue residue : selectedResidues) {
            fileWriter.write("\n")
            Rotamer[] rotamers = residue.getRotamers()
            for (int rotIndex=0; rotIndex < rotamers.length; rotIndex++) {
                String rotPop = format("%.6f", populationArray[residueIndex][rotIndex])
                fileWriter.write(residue.getName() + residue.getResidueNumber() + "\t" +
                        rotamers[rotIndex].toString() + "\t" + rotPop + "\n")

            }
            residueIndex += 1
        }
        fileWriter.close()
        System.out.println("\n Successfully wrote to the populations file.")

        int[][] conformers = rotamerOptimization.getConformers()
        char[] altLocs = new char[]{'C', 'B', 'A'}
        List<String> residueChainNum = new ArrayList<>()
        for(Residue residue: activeAssembly.getResidueList()){
            char chain = residue.getChainID()
            int resNum = residue.getResidueNumber()
            residueChainNum.add(String.valueOf(chain)+String.valueOf(resNum))
        }

        File structureFile = new File(filename)
        //List<String> rotNames = new ArrayList<>()
        int[] optimalRotamers = new int[selectedResidues.size()]
        int assemblyIndex = 0
        String[] rotNames = new String[selectedResidues.size()]
        for(int confIndex=2; confIndex > -1; confIndex--){
            List<Residue> conformerResidueList = new ArrayList<>()
            MolecularAssembly conformerAssembly = algorithmFunctions.open(filename)
            if(manyBodyOptions.getTitration()){
                logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n")

                // Collect residue numbers.
                List<Integer> resNumberList = new ArrayList<>()
                for (Residue residue : residues) {
                    resNumberList.add(residue.getResidueNumber())
                }

                // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
                titrationManyBody = new TitrationManyBody(filename, activeAssembly.getForceField(),
                        resNumberList, titrationPH, manyBodyOptions)
                conformerAssembly = titrationManyBody.getProtonatedAssembly()
                potentialEnergy = conformerAssembly.getPotentialEnergy()
            }
            for(int resIndex = 0; resIndex < selectedResidues.size(); resIndex++){
                Residue residueSelect = selectedResidues.get(resIndex)
                String resChainNum = String.valueOf(residueSelect.getChainID()) + String.valueOf(residueSelect.getResidueNumber())
                int index = residueChainNum.indexOf(resChainNum)
                Residue residue = conformerAssembly.getResidueList().get(index)
                conformerResidueList.add(residue)
                residue.setRotamers(manyBodyOptions.getRotamerLibrary(true))
                Rotamer[] rotamers = residue.getRotamers()
                int rotIndex = conformers[resIndex][confIndex]
                if(populationArray[resIndex][rotIndex]  != 0){
                    optimalRotamers[resIndex] = rotIndex
                    RotamerLibrary.applyRotamer(residue, rotamers[rotIndex])
                    double occupancy = populationArray[resIndex][rotIndex]
                    boolean diffStates = false
                    for (int i = 2; i > -1; i--) {
                        int rotamerInd = conformers[resIndex][i]
                        String rotName = rotamers[rotamerInd].getName()
                        double occupancyTest = populationArray[resIndex][rotamerInd]
                        logger.info("Residue index: " + resIndex)
                        if (i == 2 && rotNames[resIndex] != null) {
                            rotNames[resIndex] = rotName
                        } else if (i < 2 && occupancyTest != 0 && !rotNames[resIndex].contains(rotName)) {
                            diffStates = true
                            String newString = rotNames[resIndex] + rotName
                            logger.info(newString)
                            rotNames[resIndex] = newString
                        } else if (i == 2) {
                            rotNames[resIndex] = rotName
                        }
                    }
                    logger.info(rotNames.toString())
                    for(Atom atom: residue.getAtomList()){
                        if(!residue.getBackboneAtoms().contains(atom) || diffStates){
                            if(occupancy != 1){
                                atom.setAltLoc(altLocs[confIndex])
                                atom.setOccupancy(occupancy)
                            } else {
                                atom.setOccupancy(occupancy)
                                atom.setAltLoc(' ' as Character)
                            }
                        } else {
                            atom.setOccupancy(1.0)
                            atom.setAltLoc(' ' as Character)
                        }
                    }
                }

            }

            conformerAssemblies[assemblyIndex] = conformerAssembly
            titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, conformerAssemblies[assemblyIndex], conformerResidueList)
            excludeAtoms.addAll(titrationManyBody.getExcludeAtoms())
            assemblyIndex++
        }
        PDBFilter pdbFilter = new PDBFilter(structureFile, Arrays.asList(conformerAssemblies), forceField, properties)
        pdbFilter.writeFile(structureFile, false, excludeAtoms, true, true)
        System.setProperty("standardizeAtomNames", "false")

        return this
    }
}
