package ffx.xray.groovy

import edu.rit.pj.ParallelTeam
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.algorithms.optimize.TitrationManyBody
import ffx.numerics.Potential
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
import ffx.xray.DiffractionData
import ffx.xray.RefinementEnergy
import ffx.xray.RefinementMinimize
import ffx.xray.cli.XrayOptions
import org.apache.commons.configuration2.CompositeConfiguration
import picocli.CommandLine

import java.util.stream.Collectors

import static ffx.potential.bonded.NamingUtils.renameAtomsToPDBStandard
import static java.lang.String.format
import static java.lang.String.valueOf
import static org.apache.commons.io.FilenameUtils.removeExtension
import static org.apache.commons.io.FilenameUtils.removeExtension
import static org.apache.commons.io.FilenameUtils.removeExtension

@CommandLine.Command(description = " Run GenZ function for free energy change.", name = "ffxc GenZ")
class GenZ extends AlgorithmsScript {

    @CommandLine.Mixin
    ManyBodyOptions manyBodyOptions

    @CommandLine.Mixin
    AlchemicalOptions alchemicalOptions

    @CommandLine.Mixin
    XrayOptions xrayOptions

    @CommandLine.Option(names = ["--rEE", "--ro-ensembleEnergy"], paramLabel = "0.0",
            description = "Keep permutations within ensemble Energy kcal/mol from the GMEC.")
    private String ensembleEnergy = "0.0"

    @CommandLine.Option(names = ["--pF", "--printFiles"], paramLabel = "false",
            description = "Write to an energy restart file and ensemble file.")
    private boolean printFiles = false

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

        xrayOptions.init()
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
        int[] optimalRotamers
        Set<Atom> excludeAtoms = new HashSet<>()
        boolean isTitrating = false

        // The refinement mode must be coordinates.
        if (xrayOptions.refinementMode != RefinementMinimize.RefinementMode.COORDINATES) {
            logger.info(" Refinement mode set to COORDINATES.")
            xrayOptions.refinementMode = RefinementMinimize.RefinementMode.COORDINATES
        }

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
                    resNumberList, titrationPH)
            MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly()
            setActiveAssembly(protonatedAssembly)
            potentialEnergy = protonatedAssembly.getPotentialEnergy()
        }

        // Load parsed X-ray properties.
        xrayOptions.setProperties(parseResult, properties)

        // Set up the diffraction data, which could be multiple files.
        DiffractionData diffractionData = xrayOptions.getDiffractionData(filenames, molecularAssemblies, properties)
        refinementEnergy = xrayOptions.toXrayEnergy(diffractionData)
        refinementEnergy.setScaling(null)

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


        RotamerLibrary.measureRotamers(residueList, false)

        rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residueList.size()))

        selectedResidues = rotamerOptimization.getResidues()

        int[] currentRotamers = new int[selectedResidues.size()]

        //Calculate possible permutations for assembly
        rotamerOptimization.getPopulations(selectedResidues.toArray() as Residue[], 0, currentRotamers)

        //Collect the Bolztmann weights and calculated offset of each assembly
        boltzmannWeights = rotamerOptimization.getTotalBoltzmann()
        offsets = rotamerOptimization.getRefEnergy()

        //Calculate the populations for the all residue rotamers
        populationArray = rotamerOptimization.getFraction()
        logger.info("Population Array size: " + Arrays.toString(populationArray[0]))

        optimalRotamers = rotamerOptimization.getOptimumRotamers()
        if (manyBodyOptions.getTitration()) {
            isTitrating = titrationManyBody.excludeExcessAtoms(excludeAtoms, optimalRotamers, selectedResidues)
        }

        FileWriter fileWriter = new FileWriter("populations.txt")
        int residueIndex = 0
        for (Residue residue : selectedResidues) {
            fileWriter.write("\n")
            // Set sums for to protonated, deprotonated, and tautomer states of titratable residues
            double protSum = 0
            double deprotSum = 0
            double tautomerSum = 0
            Rotamer[] rotamers = residue.getRotamers()
            for (int rotIndex=0; rotIndex < rotamers.length; rotIndex++) {
                String rotPop = format("%.6f", populationArray[residueIndex][rotIndex])
                fileWriter.write(residue.getName() + residue.getResidueNumber() + "\t" +
                        rotamers[rotIndex].toString() + "\t" + rotPop + "\n")
                if (manyBodyOptions.titration) {
                    switch (rotamers[rotIndex].getName()) {
                        case "HIS":
                        case "LYS":
                        case "GLH":
                        case "ASH":
                        case "CYS":
                            protSum += populationArray[residueIndex][rotIndex]
                            break
                        case "HIE":
                        case "LYD":
                        case "GLU":
                        case "ASP":
                        case "CYD":
                            deprotSum += populationArray[residueIndex][rotIndex]
                            break
                        case "HID":
                            tautomerSum += populationArray[residueIndex][rotIndex]
                            break
                        default:
                            break
                    }
                }

            }
            if (manyBodyOptions.titration){
                String formatedProtSum = format("%.6f", protSum)
                String formatedDeprotSum = format("%.6f", deprotSum)
                String formatedTautomerSum = format("%.6f", tautomerSum)
                switch (residue.getName()) {
                    case "HIS":
                    case "HIE":
                    case "HID":
                        logger.info(residue.getResidueNumber() + "\tHIS" + "\t" + formatedProtSum + "\t" +
                                "HIE" + "\t" + formatedDeprotSum + "\t" +
                                "HID" + "\t" + formatedTautomerSum)
                        break
                    case "LYS":
                    case "LYD":
                        logger.info(residue.getResidueNumber() + "\tLYS" + "\t" + formatedProtSum + "\t" +
                                "LYD" + "\t" + formatedDeprotSum)
                        break
                    case "ASH":
                    case "ASP":
                        logger.info(residue.getResidueNumber() + "\tASP" + "\t" + formatedDeprotSum + "\t" +
                                "ASH" + "\t" + formatedProtSum)
                        break
                    case "GLH":
                    case "GLU":
                        logger.info(residue.getResidueNumber() + "\tGLU" + "\t" + formatedDeprotSum + "\t" +
                                "GLH" + "\t" + formatedProtSum)
                        break
                    case "CYS":
                    case "CYD":
                        logger.info(residue.getResidueNumber() + "\tCYS" + "\t" + formatedProtSum + "\t" +
                                "CYD" + "\t" + formatedDeprotSum)
                        break
                    default:
                        break
                }
            }
            residueIndex += 1
        }
        fileWriter.close()
        System.out.println("\n Successfully wrote to the populations file.")

        int[][] conformers = rotamerOptimization.getConformers()
        char[] altLocs = new char[]{'C', 'B', 'A'}
        List<Residue> residueAList = new ArrayList<>()
        List<double[]> atomListA = new ArrayList<>()
        List<String> residueChainNum = new ArrayList<>()
        for(Residue residue: activeAssembly.getResidueList()){
            char chain = residue.getChainID()
            int resNum = residue.getResidueNumber()
            residueChainNum.add(String.valueOf(chain)+String.valueOf(resNum))
        }

        File structureFile = new File(filename)
        int assemblyIndex = 0
        for(int confIndex=2; confIndex > -1; confIndex--){
            MolecularAssembly conformerAssembly = algorithmFunctions.open(filename)
            for(int resIndex = 0; resIndex < selectedResidues.size(); resIndex++){
                Residue residueSelect = selectedResidues.get(resIndex)
                String resChainNum = String.valueOf(residueSelect.getChainID()) + String.valueOf(residueSelect.getResidueNumber())
                int index = residueChainNum.indexOf(resChainNum)
                Residue residue = conformerAssembly.getResidueList().get(index)
                residue.setRotamers(manyBodyOptions.getRotamerLibrary(true))
                Rotamer[] rotamers = residue.getRotamers()
                if(conformers[resIndex][confIndex] != 0 || confIndex == 2){
                    int rotIndex = conformers[resIndex][confIndex]
                    RotamerLibrary.applyRotamer(residue, rotamers[rotIndex])
                        for(Atom atom: residue.getAtomList()){
                            String name = atom.getName()
                            if(!residue.getBackboneAtoms().contains(atom) && populationArray[resIndex][rotIndex] != 0 ||
                                    !residue.getBackboneAtoms().contains(atom) &&
                                    populationArray[resIndex][conformers[resIndex][1]] != 0 && confIndex == 2){
                                atom.setAltLoc(altLocs[confIndex])
                                double occupancy = populationArray[resIndex][rotIndex]
                                atom.setOccupancy(occupancy)
                                //double[] atomCoor = new double[]{atom.getResidueNumber(), atom.getX(), atom.getY(), atom.getZ()}
                                //atomListA.add(atomCoor)name != 'CA' && name != 'O' && name != 'C' && name != 'N'&& name != 'OXT'&& name != 'H'
                                //                                    && name != 'OT2' && name != 'H1' && name != 'H2'&& name != 'H3'&& name != 'HA'
                                //                                    && name != 'HA2' && name != 'HA3'
                            }
                        }
                    /*if(confIndex == 2){
                            residueAList.add(residue)
                    } else {
                        List<Atom> atomListB = residue.getAtomList()
                        double[][] atomACoor = new double[residue.getAtomList().size()][3]
                        int count = 0
                        for(int atomIndex = 0; atomIndex < atomListA.size(); atomIndex++){
                            if(atomListA.get(atomIndex)[0] == residue.getResidueNumber()){
                                for(int k = 1; k<4; k++){
                                    atomACoor[count][k-1] = atomListA.get(atomIndex)[k]
                                }
                                count++
                            }
                        }
                        for (int i = 0; i < residue.getAtomList().size(); i++) {
                            double coorAX = atomACoor[i][0]
                            double coorAY = atomACoor[i][1]
                            double coorAZ = atomACoor[i][2]
                            Atom atom = atomListB.get(i)
                            double coorBX = atom.getX()
                            double coorBY = atom.getY()
                            double coorBZ = atom.getZ()
                            if (coorAX == coorBX && coorAY == coorBY && coorAZ == coorBZ) {
                                atom.setAltLoc(' ' as Character)
                                int aResidueIndex = conformerAssemblies[0].getResidueList().indexOf(residue)
                                Residue aResidue = conformerAssemblies[0].getResidueList().get(aResidueIndex)
                                Atom atomA = aResidue.getAtomList()[i]
                                atomA.setAltLoc(' ' as Character)
                                atomA.setOccupancy(1.0)
                            } else {
                                atom.setAltLoc(altLocs[confIndex])
                                double occupancy = populationArray[resIndex][rotIndex]
                                atom.setOccupancy(occupancy)
                                int rotIndexA = conformers[resIndex][2]
                                int aResidueIndex = conformerAssemblies[0].getResidueList().indexOf(residue)
                                Residue aResidue = conformerAssemblies[0].getResidueList().get(aResidueIndex)
                                Atom atomA = aResidue.getAtomList()[i]
                                atomA.setAltLoc('A' as Character)
                                double occupancyA = populationArray[resIndex][rotIndexA]
                                atomA.setOccupancy(occupancyA)
                            }
                        }

                    }*/
                }

            }
            conformerAssemblies[assemblyIndex] = conformerAssembly
            assemblyIndex++
            logger.info("The assembly index " + assemblyIndex)
        }

        // Print the final energy of each conformer.
        /*algorithmFunctions.energy(conformerAssemblies)
        DiffractionData diffractionDataFinal = xrayOptions.getDiffractionData(filenames, conformerAssemblies, properties)
        refinementEnergy = xrayOptions.toXrayEnergy(diffractionDataFinal)
        refinementEnergy.setScaling(null)*/
        String remark = "String"
        PDBFilter pdbFilter = new PDBFilter(structureFile, Arrays.asList(conformerAssemblies), forceField, properties)
        pdbFilter.writeFile(structureFile, false, excludeAtoms, true, true)
        /*if (titrationPH > 0) {
            diffractionDataFinal.writeModel(removeExtension(filenames[0]) + ".pdb", excludeAtoms, titrationPH)
        } else {

        }*/
        System.setProperty("standardizeAtomNames", "false")


        return this
    }


    @Override
    List<Potential> getPotentials() {
        if (conformerAssemblies == null) {
            return new ArrayList<Potential>()
        } else {
            return Arrays.stream(conformerAssemblies).
                    filter {a -> a != null
                    }.map {a -> a.getPotentialEnergy()
            }.filter {e -> e != null
            }.collect(Collectors.toList())
        }
    }

}
