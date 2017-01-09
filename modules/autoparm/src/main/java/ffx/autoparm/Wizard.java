/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.autoparm;

// Java Imports
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * Convert SDF file data to SMILES strings, fragment
 *
 * @author Rae Ann Corrigan
 */
public class Wizard {

    protected Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
    protected String filename;
    private String fullSmiles;

    private final static Logger logger = Logger.getLogger(Wizard.class.getName());

    public Wizard(String sdffile) {
        this.filename = sdffile;
    }

    public String readSDF() throws FileNotFoundException, IOException {
        int recordsRead = 0;
        File file = new File(filename);
        BufferedReader bufferedReader = null;
        try {
            FileReader fileReader = new FileReader(file);
            bufferedReader = new BufferedReader(fileReader);
        } catch (IOException e) {
            logger.warning(String.format(" Could not open %s.\n %s\n", filename, e.toString()));
            return null;
        }

        IteratingSDFReader reader = null;
        try {
            reader = new IteratingSDFReader(bufferedReader, SilentChemObjectBuilder.getInstance());
            logger.info(String.format("\n Reading %s", file.getAbsoluteFile()));
            while (reader.hasNext()) {
                IAtomContainer molecule = reader.next();
                recordsRead++;
                try {
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                    CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(molecule);
                    aromaticity.apply(molecule);
                    /**
                     * Generate SMILES.
                     */
                    fullSmiles = assignSMILES(molecule);
                    /**
                     * Return the smiles string of the first molecule.
                     */
                    break;
                } catch (Exception x) {
                    logger.log(Level.SEVERE, String.format("[Record %d] Error %s\n", recordsRead, file.getAbsoluteFile()), x);
                }
            }
        } catch (Exception x) {
            logger.log(Level.SEVERE, String.format("[Record %d] Error %s\n", recordsRead, file.getAbsoluteFile()), x);
        } finally {
            try {
                reader.close();
            } catch (Exception x) {
                logger.warning(" Error closing file reader\n " + x.toString());
            }
        }
        return fullSmiles;
    }

    protected String assignSMILES(IAtomContainer molecule) throws Exception {
        //String absoluteSmiles = SmilesGenerator.absolute().create(molecule);
        //String uniqueSmiles = SmilesGenerator.unique().create(molecule);
        //String genericSmiles = SmilesGenerator.generic().create(molecule);
        String isomericSmiles = SmilesGenerator.isomeric().create(molecule);

        SmilesGenerator aromaticSmileGenerator = new SmilesGenerator().aromatic();
        String aromaticSmiles = aromaticSmileGenerator.create(molecule);

        //logger.info(String.format(" Absolute Smiles: %s", absoluteSmiles));
        //logger.info(String.format(" Unique Smiles: %s", uniqueSmiles));
        //logger.info(String.format(" Generic Smiles: %s", genericSmiles));
        logger.info(String.format(" Isomeric Smiles: %s", isomericSmiles));
        logger.info(String.format(" Aromatic Smiles: %s", aromaticSmiles));
        return aromaticSmiles;
    }

}
