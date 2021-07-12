//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
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

package ffx.algorithms.groovy

import ffx.potential.MolecularAssembly
import ffx.potential.bonded.Atom
import ffx.potential.bonded.PolymerUtils
import ffx.potential.cli.PotentialScript
import ffx.potential.parsers.PDBFilter
import ffx.potential.utils.PotentialsUtils
import ffx.utilities.FFXScript
import ffx.potential.bonded.PolymerUtils
import ffx.potential.bonded.Residue
import picocli.CommandLine
import java.util.ArrayList

import static ffx.numerics.math.DoubleMath.dist
import static ffx.numerics.math.DoubleMath.log
import static java.lang.String.format
import org.apache.commons.io.FilenameUtils

@CommandLine.Command(description = " Fix chain breaks in a pdb file.", name = "ffxc ChainBreaks")
class ChainBreaks extends PotentialScript {

    /**
     * The final argument(s) should be one or more filenames.
     */
    @CommandLine.Parameters(arity = "1", paramLabel = "files",
            description = 'The atomic coordinate file in PDB or XYZ format.')
    private List<String> filenames = null

    MolecularAssembly molecularAssembly
    List<String> chainBreaks = new ArrayList<>()
    List<double[]> cCoor = new ArrayList<>()
    List<double[]> nCoor = new ArrayList<>()
    List<double[]> newCoordinates
    PotentialsUtils potentialsUtils = new PotentialsUtils()




    /**
     * ChainBreaks constructor.
     */
    ChainBreaks() {
        this(new Binding())
    }

    /**
     * ChainBreaks constructor.
     * @param binding The Groovy Binding to use.
     */
    ChainBreaks(Binding binding) {
        super(binding)
    }

    /**
     * Execute the script.
     */
    @Override
    ChainBreaks run() {

        if (!init()) {
            return this
        }

        molecularAssembly = potentialFunctions.open(filenames.get(0))
        List<Residue> residues = molecularAssembly.getResidueList()
        chainBreaks = findChainBreaks(residues, 3)
        logger.info(format(" Fixing Chain Breaks in %s", filenames.get(0)))
        String pdbName = FilenameUtils.getBaseName(filenames.get(0))
        String newPDBpath = FilenameUtils.getFullPath(filenames.get(0)).replace(filenames.get(0),"")  + pdbName + "_edited.pdb"
        File newPDBFile = new File (newPDBpath)
        logger.info(format(" Saving New Coordinates to %s", filenames.get(0) + "_edited.pdb"))
        potentialsUtils.saveAsPDB(molecularAssembly, newPDBFile, false, false)
        
        return this


    }

    private List<String> findChainBreaks(List<Residue> residues, double cutoff) {
        List<List<Residue>> subChains = new ArrayList<>()
        List<String> chainBreaks = new ArrayList<>()

        // Chain-start atom: N (amino)/O5* (nucleic)
        // Chain-end atom:   C (amino)/O3* (nucleic)
        Residue.ResidueType rType = residues.get(0).getResidueType();
        String startAtName;
        String endAtName;
        switch (rType) {
            case Residue.ResidueType.AA:
                startAtName = "N";
                endAtName = "C";
                break;
            case Residue.ResidueType.NA:
                boolean namedStar =
                        residues.stream()
                                .flatMap((Residue r) -> r.getAtomList().stream())
                                .anyMatch((Atom a) -> a.getName().equals("O5*"));
                if (namedStar) {
                    startAtName = "O5*";
                    endAtName = "O3*";
                } else {
                    startAtName = "O5\'";
                    endAtName = "O3\'";
                }
                break;
            case Residue.ResidueType.UNK:
            default:
                logger.fine(
                        " Not attempting to find chain breaks for chain with residue "
                                + residues.get(0).toString());
                List<List<Residue>> retList = new ArrayList<>();
                retList.add(residues);
                return retList;
        }

        List<Residue> subChain = null;
        Residue previousResidue = null;
        Atom priorEndAtom = null;

        for (Residue residue : residues) {
            List<Atom> resAtoms = residue.getAtomList();
            if (priorEndAtom == null) {
                // Initialization.
                subChain = new ArrayList<>();
                subChain.add(residue);
                subChains.add(subChain);
            } else {
                // Find the start atom of the current residue.
                Atom startAtom = null;
                for (Atom a : resAtoms) {
                    if (a.getName().equalsIgnoreCase(startAtName)) {
                        startAtom = a;
                        break;
                    }
                }
                if (startAtom == null) {
                    subChain.add(residue);
                    continue;
                }
                // Compute the distance between the previous carbonyl carbon and the current nitrogen.
                double r = dist(priorEndAtom.getXYZ(null), startAtom.getXYZ(null));

                if (r > cutoff) {
                    // Start a new chain.
                    subChain = new ArrayList<>()
                    subChain.add(residue)
                    subChains.add(subChain)
                    char ch1 = previousResidue.getChainID()
                    char ch2 = residue.getChainID();
                    chainBreaks.add("C " + previousResidue.toString() + " N " + residue.toString())

                    fixChainBreaks(priorEndAtom.getXYZ(null), startAtom.getXYZ(null))
                    priorEndAtom.setXYZ(newCoordinates.get(0))
                    startAtom.setXYZ(newCoordinates.get(1))

                } else {
                    // Continue the current chain.
                    subChain.add(residue)
                }
            }

            // Save the carbonyl carbon.
            for (Atom a : resAtoms) {
                if (a.getName().equalsIgnoreCase(endAtName)) {
                    priorEndAtom = a
                    break
                }
            }
            previousResidue = residue
        }

        return chainBreaks
    }

    private void fixChainBreaks(double[] cCoordinates, double[] nCoordinates){
        logger.info("Generating new coordinates")
        newCoordinates = new ArrayList<>()
        double distance = dist(cCoordinates,nCoordinates)
        while (distance > 3){
            double dx = Math.abs((cCoordinates[0] - nCoordinates[0]) / 4)
            double dy = Math.abs((cCoordinates[1] - nCoordinates[1]) / 4)
            double dz = Math.abs((cCoordinates[2] - nCoordinates[2]) / 4)

            if (cCoordinates[0] > nCoordinates[0]){
                cCoordinates[0] = cCoordinates[0] - dx
                nCoordinates[0] = nCoordinates[0] + dx
            } else if (cCoordinates[0] < nCoordinates[0]){
                cCoordinates[0] = cCoordinates[0] + dx
                nCoordinates[0] = nCoordinates[0] - dx
            }

            if (cCoordinates[1] > nCoordinates[1]){
                cCoordinates[1] = cCoordinates[1] - dy
                nCoordinates[1] = nCoordinates[1] + dy
            } else if (cCoordinates[1] < nCoordinates[1]){
                cCoordinates[1] = cCoordinates[1] + dy
                nCoordinates[1] = nCoordinates[1] - dy
            }

            if (cCoordinates[2] > nCoordinates[2]){
                cCoordinates[2] = cCoordinates[2] - dz
                nCoordinates[2] = nCoordinates[2] + dz
            } else if (cCoordinates[2] < nCoordinates[0]){
                cCoordinates[2] = cCoordinates[2] + dz
                nCoordinates[2] = nCoordinates[2] - dz
            }

            distance = dist(cCoordinates,nCoordinates)
        }
        newCoordinates.add(cCoordinates)
        newCoordinates.add(nCoordinates)

    }


}


