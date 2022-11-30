package ffx.algorithms.groovy

import edu.rit.pj.Comm
import ffx.algorithms.TitrationManyBody
import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.ManyBodyOptions
import ffx.algorithms.optimize.RotamerOptimization
import ffx.numerics.math.DoubleMath
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.bonded.RotamerLibrary
import ffx.potential.cli.AlchemicalOptions
import ffx.potential.parameters.ForceField
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
class ReducedPartition extends  AlgorithmsScript{

    @Mixin
    ManyBodyOptions manyBodyOptions

    @Mixin
    AlchemicalOptions alchemicalOptions

    @CommandLine.Option(names = ["--mR", "--mutatingResidue"], paramLabel = "1",
            description = "The residue that is mutating.")
    private int mutatingResidue = 1

    @CommandLine.Option(names = ["--dC", "--distanceCutoff"], paramLabel = "10.0",
            description = "Residues within the distance cutoff from the mutating residue will be optimized.")
    private double distanceCutoff = 10.0

    @CommandLine.Option(names = ["-n", "--residueName"], paramLabel = "ALA",
            description = "Mutant residue.")
    private String resName

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
        if (titrationPH > 0) {
            System.setProperty("manybody-titration", "true")
        }
        activeAssembly = getActiveAssembly(filenames.get(0))
        double[] boltzmannWeights = new double[2]
        double[] offsets = new double[2]
        List<Residue> residueList = activeAssembly.getResidueList()

        List<Integer> residueNumber = new ArrayList<>()
        for(Residue residue : residueList){
            residueNumber.add(residue.getResidueNumber())
        }


        String mutatedFileName = ""
        //Call the MutatePDB script
        if(filenames.size() == 1){
            mutatorBinding = new Binding('-r', mutatingResidue.toString(), '-n', resName, filenames.get(0))
            MutatePDB mutatePDB = new MutatePDB(mutatorBinding)
            mutatePDB.run()
            mutatedFileName = mutatorBinding.getProperty('versionFileName')
        }

        double[] mutatingResCoor = new double[3]
        int index = residueNumber.indexOf(mutatingResidue)
        mutatingResCoor = residueList.get(index).getAtomByName("CA", true).getXYZ(mutatingResCoor)
        String listResidues = ""
        int count = 0
        for(int k=0; k<residueList.size(); k++){
            double[] currentResCoor = new double[3]
            currentResCoor = residueList.get(k).getAtomByName("CA", true).getXYZ(currentResCoor)
            double dist = DoubleMath.dist(mutatingResCoor,currentResCoor)
            if(dist < distanceCutoff){
                if(count == 0){
                    listResidues += "A" + residueList.get(k).getResidueNumber()
                } else {
                    listResidues += ",A" + residueList.get(k).getResidueNumber()
                }
                count++
            }
        }

        for(int j=0; j<2; j++){
            // Load the MolecularAssembly.
            if(j>0){
                if(filenames.size() == 1){
                    mutatedAssembly = getActiveAssembly(mutatedFileName)
                    setActiveAssembly(mutatedAssembly)
                    logger.info(activeAssembly.getResidueList().toString())
                    activeAssembly.getPotentialEnergy().energy()
                } else{
                    setActiveAssembly(getActiveAssembly(filenames.get(j)))
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
            manyBodyOptions.setListResidues(listResidues)
            // Collect residues to optimize.
            List<Residue> residues = manyBodyOptions.collectResidues(activeAssembly)
            if (residues == null || residues.isEmpty()) {
                logger.info(" There are no residues in the active system to optimize.")
                return this
            }

            // Handle rotamer optimization with titration.
            if (titrationPH > 0) {
                logger.info("\n Adding titration hydrogen to : " + filenames.get(0) + "\n")

                // Collect residue numbers.
                List<Integer> resNumberList = new ArrayList<>()
                for (Residue residue : residues) {
                    resNumberList.add(residue.getResidueNumber())
                }

                // Create new MolecularAssembly with additional protons and update the ForceFieldEnergy
                titrationManyBody = new TitrationManyBody(filenames.get(0), activeAssembly.getForceField(),
                        resNumberList, titrationPH)
                MolecularAssembly protonatedAssembly = titrationManyBody.getProtonatedAssembly()
                setActiveAssembly(protonatedAssembly)
                potentialEnergy = protonatedAssembly.getPotentialEnergy()
            }



            RotamerOptimization rotamerOptimization = new RotamerOptimization(activeAssembly,
                    potentialEnergy, algorithmListener)

            manyBodyOptions.initRotamerOptimization(rotamerOptimization, activeAssembly)

            // rotamerOptimization.getResidues() returns a cached version of
            // manyBodyOptions.collectResidues(activeAssembly)
            List<Residue> residues1 = rotamerOptimization.getResidues()

            logger.info("\n Initial Potential Energy:")
            potentialEnergy.energy(false, true)

            logger.info("\n Initial Rotamer Torsion Angles:")
            RotamerLibrary.measureRotamers(residues1, false)

            // Run the optimization.
            rotamerOptimization.optimize(manyBodyOptions.getAlgorithm(residues1.size()))

            int[] currentRotamers = new int[residues1.size()]

            boltzmannWeights[j] = rotamerOptimization.partitionFunction(residues1.toArray() as Residue[], 0, currentRotamers)
            offsets[j] = rotamerOptimization.getRefEnergy()
        }

        double gibbs = -(0.6)*(Math.log(boltzmannWeights[1]/boltzmannWeights[0]))
        logger.info("\n Gibbs Free Energy Change: " + gibbs)


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
