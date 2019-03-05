package ffx.potential.groovy

import ffx.potential.ForceFieldEnergy
import ffx.potential.bonded.Atom
import ffx.potential.parsers.XYZFilter
import org.apache.commons.io.FilenameUtils

import ffx.potential.MolecularAssembly
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.SystemFilter

import picocli.CommandLine.Command
import picocli.CommandLine.Parameters

import static java.lang.String.format

/**
 * The Superpose script superposes all molecules in an arc file to the first molecule in the arc file and reports the RMSD.
 * <br>
 * Usage:
 * <br>
 * ffxc Superpose [options] &lt;filename&gt;
 */
@Command(description = " Save the system as a PDB file.", name = "ffxc SaveAsPDB")
class Superpose extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    List<String> filenames = null

    private File baseDir = null

    void setBaseDir(File baseDir) {
        this.baseDir = baseDir
    }

    public ForceFieldEnergy forceFieldEnergy = null


    /**
     * Execute the script.
     */
    @Override
    Superpose run() {
        if (!init()) {
            return this
        }

        MolecularAssembly[] assemblies
        if (filenames != null && filenames.size() > 0) {
            assemblies = potentialFunctions.open(filenames.get(0))
            activeAssembly = assemblies[0]
        } else if (activeAssembly == null) {
            logger.info(helpString())
            return this
        }

        String modelFilename = activeAssembly.getFile().getAbsolutePath()

        forceFieldEnergy = activeAssembly.getPotentialEnergy()
        Atom[] atoms = activeAssembly.getAtomArray()

        int nVars = forceFieldEnergy.getNumberOfVariables()
        double[] x = new double[nVars]
        forceFieldEnergy.getCoordinates(x)

        SystemFilter systemFilter = potentialFunctions.getFilter()
        if (systemFilter instanceof XYZFilter) {
            XYZFilter xyzFilter = (XYZFilter) systemFilter

            double[] x2 = new double[nVars]
            double[] mass = new double[nVars / 3]

            int nAtoms = atoms.length;
            for (int i=0; i<nAtoms; i++) {
                mass[i] = atoms[i].getMass()
            }

            //Get heavy atom masses.
            int nHeavyVars = forceFieldEnergy.getNumberOfHeavyAtomVariables()
            double[] massHeavy = new double [nHeavyVars/3]
            for (int i=0; i<nHeavyVars/3; i++){
                if(!atoms[i].isHydrogen()) {
                    massHeavy[i] = atoms[i].getMass()
                }
            }

            //Array containing heavy atom indices.
            double[] heavyAtomPositions = new double[nHeavyVars/3];
            int j = 0;
            for(int i=0; i<nVars/3; i++){
                if(!atoms[i].isHydrogen()){
                    heavyAtomPositions[j]=i
                    j++
                }
            }

            while (xyzFilter.readNext()) {
                //Arrays for holding coordinates of heavy atoms after rotation and translation.
                double[] xHeavy = new double[nHeavyVars]
                double[] x2Heavy = new double[nHeavyVars]

                forceFieldEnergy.getCoordinates(x2)

                //Original RMSD.
                for(int i = 0; i < nHeavyVars/3; i++){
                    int positionOfHeavyAtom = heavyAtomPositions[i]
                    xHeavy[i*3]=x[positionOfHeavyAtom]
                    xHeavy[i*3+1]=x[positionOfHeavyAtom+1]
                    xHeavy[i*3+2]=x[positionOfHeavyAtom+2]
                    x2Heavy[i*3]=x2[positionOfHeavyAtom]
                    x2Heavy[i*3+1]=x2[positionOfHeavyAtom+1]
                    x2Heavy[i*3+2]=x2[positionOfHeavyAtom+2]
                }
                double origRMSDHeavy = ffx.potential.utils.Superpose.rmsd(xHeavy, x2Heavy, massHeavy);

                //Translated RMSD.
                ffx.potential.utils.Superpose.translate(x, mass, x2, mass)
                for(int i = 0; i < nHeavyVars/3; i++){
                    int positionOfHeavyAtom = heavyAtomPositions[i]
                    xHeavy[i*3]=x[positionOfHeavyAtom]
                    xHeavy[i*3+1]=x[positionOfHeavyAtom+1]
                    xHeavy[i*3+2]=x[positionOfHeavyAtom+2]
                    x2Heavy[i*3]=x2[positionOfHeavyAtom]
                    x2Heavy[i*3+1]=x2[positionOfHeavyAtom+1]
                    x2Heavy[i*3+2]=x2[positionOfHeavyAtom+2]
                }
                double transRMSDHeavy = ffx.potential.utils.Superpose.rmsd(xHeavy, x2Heavy, massHeavy)

                //Rotated RMSD.
                ffx.potential.utils.Superpose.rotate(x, x2, mass)
                for(int i = 0; i < nHeavyVars/3; i++){
                    int positionOfHeavyAtom = heavyAtomPositions[i]
                    xHeavy[i*3]=x[positionOfHeavyAtom]
                    xHeavy[i*3+1]=x[positionOfHeavyAtom+1]
                    xHeavy[i*3+2]=x[positionOfHeavyAtom+2]
                    x2Heavy[i*3]=x2[positionOfHeavyAtom]
                    x2Heavy[i*3+1]=x2[positionOfHeavyAtom+1]
                    x2Heavy[i*3+2]=x2[positionOfHeavyAtom+2]
                }
                double rotRMSDHeavy = ffx.potential.utils.Superpose.rmsd(xHeavy, x2Heavy, massHeavy)
                logger.info(format(
                        "\n Coordinate RMSD Based On Heavy Atoms (Angstroms)\n Original:\t\t%7.3f\n After Translation:\t%7.3f\n After Rotation:\t%7.3f\n",
                        origRMSDHeavy, transRMSDHeavy, rotRMSDHeavy))
            }
        }

        logger.info("\n Saving ARC for " + modelFilename)

        File saveDir = baseDir
        if (saveDir == null || !saveDir.exists() || !saveDir.isDirectory() || !saveDir.canWrite()) {
            saveDir = new File(FilenameUtils.getFullPath(modelFilename))
        }
        String dirName = saveDir.toString() + File.separator
        String fileName = FilenameUtils.getName(modelFilename)
        fileName = FilenameUtils.removeExtension(fileName) + ".pdb"
        File saveFile = potentialFunctions.versionFile(new File(dirName + fileName));
        potentialFunctions.saveAsPDB(activeAssembly, saveFile)

        return this
    }
}

/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2019.
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */