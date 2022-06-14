//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.potential.groovy.test

import com.google.common.collect.MinMaxPriorityQueue
import ffx.crystal.Crystal
import ffx.numerics.Potential
import ffx.potential.AssemblyState
import ffx.potential.ForceFieldEnergy
import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Residue
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.PotentialScript
import ffx.potential.extended.ExtendedSystem
import ffx.potential.parsers.PDBFilter
import ffx.potential.parsers.SystemFilter
import ffx.potential.parsers.XPHFilter
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

import static ffx.potential.utils.StructureMetrics.momentsOfInertia
import static ffx.potential.utils.StructureMetrics.radiusOfGyration
import static java.lang.String.format
import static org.apache.commons.io.FilenameUtils.*

/**
 * The Energy script evaluates the energy of a system.
 * <br>
 * Usage:
 * <br>
 * ffxc Energy &lt;filename&gt;
 */
@Command(description = " Compute the force field potential energy.", name = "ffxc Energy")
class PhEnergy extends PotentialScript {

    @Mixin
    AtomSelectionOptions atomSelectionOptions

    /**
     * -m or --moments print out electrostatic moments.
     */
    @Option(names = ['-m', '--moments'], paramLabel = "false", defaultValue = "false",
            description = 'Print out electrostatic moments.')
    private boolean moments = false

    /**
     * --rg or --gyrate Print out the radius of gyration.
     */
    @Option(names = ['--rg', '--gyrate'], paramLabel = "false", defaultValue = "false",
            description = 'Print out the radius of gyration.')
    private boolean gyrate = false

    /**
     * --in or --inertia Print out the moments of inertia.
     */
    @Option(names = ['--in', '--inertia'], paramLabel = "false", defaultValue = "false",
            description = 'Print out the moments of inertia.')
    private boolean inertia = false

    /**
     * -g or --gradient to print out gradients.
     */
    @Option(names = ['-g', '--gradient'], paramLabel = "false", defaultValue = "false",
            description = 'Compute the atomic gradient as well as energy.')
    private boolean gradient = false

    /**
     * * --fl or --findLowest Return the n lowest energy structures from an ARC or PDB file.
     */
    @Option(names = ['--fl', '--findLowest'], paramLabel = "0", defaultValue = "0",
            description = 'Return the n lowest energies from an ARC/PDB file.')
    private int fl = 0

    /**
     * -v or --verbose enables printing out all energy components for multi-snapshot files (
     * the first snapshot is always printed verbosely).
     */
    @Option(names = ['-v', '--verbose'], paramLabel = "false", defaultValue = "false",
            description = "Print out all energy components for each snapshot.")
    private boolean verbose = false

    /**
     * --dc or --densityCutoff Collect structures above a specified density.
     */
    @Option(names = ['--dc', '--densityCutoff'], paramLabel = "0.0", defaultValue = "0.0",
            description = "Create ARC file of structures above a specified density.")
    private double dCutoff = 0.0

    /**
     * --ec or --energyCutOff Collect structures below a specified energy range from the minimum energy.
     */
    @Option(names = ['--ec', '--energyCutoff'], paramLabel = "0.0", defaultValue = "0.0",
            description = "Create ARC file of structures within a specified energy of the lowest energy structure.")
    private double eCutoff = 0.0

    /**
     * --pH or --constantPH Constant pH value for the test.
     */
    @Option(names = ['--pH', '--constantPH'], paramLabel = '7.4',
            description = 'pH value for the energy evaluation. (Only applies when esvTerm is true)')
    double pH = 7.4

    @Option(names = ['--aFi', '--arcFile'], paramLabel = "traj",
            description = 'A file containing the the PDB from which to build the ExtendedSystem. There is currently no default.')
    private String arcFileName = null

    /**
     * The final argument is a PDB or XPH coordinate file.
     */
    @Parameters(arity = "1", paramLabel = "file",
            description = 'The atomic coordinate file in PDB or XPH format.')
    private String filename = null

    public double energy = 0.0
    public ForceFieldEnergy forceFieldEnergy = null
    public MolecularAssembly mola = null

    /**
     * Energy constructor.
     */
    PhEnergy() {
        this(new Binding())
    }

    /**
     * Energy constructor.
     * @param binding The Groovy Binding to use.
     */
    PhEnergy(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    PhEnergy run() {
        // Init the context and bind variables.
        if (!init()) {
            return this
        }

        // Load the MolecularAssembly.
        activeAssembly = getActiveAssembly(filename)
        if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        // Set the filename.
        filename = activeAssembly.getFile().getAbsolutePath()

        logger.info("\n Running Energy on " + filename)
        forceFieldEnergy = activeAssembly.getPotentialEnergy()


        ExtendedSystem esvSystem = new ExtendedSystem(activeAssembly, null)
        esvSystem.setConstantPh(pH)
        int numESVs = esvSystem.extendedResidueList.size()
        forceFieldEnergy.attachExtendedSystem(esvSystem)
        logger.info(format(" Attached extended system with %d residues.", numESVs))

        // Apply atom selections
        atomSelectionOptions.setActiveAtoms(activeAssembly)

        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)

        if (gradient) {
            double[] g = new double[nVars]
            int nAts = (int) (nVars / 3)
            energy = forceFieldEnergy.energyAndGradient(x, g, true)
            logger.info(format("    Atom       X, Y and Z Gradient Components (kcal/mol/A)"))
            for (int i = 0; i < nAts; i++) {
                int i3 = 3 * i
                logger.info(format(" %7d %16.8f %16.8f %16.8f", i + 1, g[i3], g[i3 + 1], g[i3 + 2]))
            }
        } else {
            energy = forceFieldEnergy.energy(x, true)
        }

        if (moments) {
            forceFieldEnergy.getPmeNode().computeMoments(activeAssembly.getActiveAtomArray(), false)
        }

        if (gyrate) {
            double rg = radiusOfGyration(activeAssembly.getActiveAtomArray())
            logger.info(format(" Radius of gyration:           %10.5f A", rg))
        }

        if (inertia) {
            double[][] inertiaValue = momentsOfInertia(activeAssembly.getActiveAtomArray(), false, true, true)
        }

        SystemFilter systemFilter = null
        if(arcFileName != null){
            File arcFile = new File(arcFileName)
            systemFilter = new XPHFilter(arcFile, activeAssembly, activeAssembly.getForceField(), activeAssembly.getProperties(), esvSystem)
        } else{
            systemFilter = potentialFunctions.getFilter()
        }

        if (systemFilter instanceof XPHFilter || systemFilter instanceof PDBFilter) {
            int index = 1
            while (systemFilter.readNext()) {
                index++
                Crystal crystal = activeAssembly.getCrystal()
                forceFieldEnergy.setCrystal(crystal)
                forceFieldEnergy.getCoordinates(x)
                if (verbose) {
                    logger.info(format(" Snapshot %4d", index))
                    if (!crystal.aperiodic()) {
                        logger.info(format("\n Density:                                %6.3f (g/cc)",
                                crystal.getDensity(activeAssembly.getMass())))
                    }
                    energy = forceFieldEnergy.energy(x, true)
                } else {
                    energy = forceFieldEnergy.energy(x, false)
                    logger.info(format(" Snapshot %4d: %16.8f (kcal/mol)", index, energy))
                }
            }
        }

        return this
    }

    @Override
    List<Potential> getPotentials() {
        List<Potential> potentials
        if (forceFieldEnergy == null) {
            potentials = Collections.emptyList()
        } else {
            potentials = Collections.singletonList((Potential) forceFieldEnergy)
        }
        return potentials
    }

    private class StateContainer implements Comparable<StateContainer> {

        private final AssemblyState state
        private final double e

        StateContainer(AssemblyState state, double e) {
            this.state = state
            this.e = e
        }

        AssemblyState getState() {
            return state
        }

        double getEnergy() {
            return e
        }

        @Override
        int compareTo(StateContainer o) {
            return Double.compare(e, o.getEnergy())
        }
    }

}

