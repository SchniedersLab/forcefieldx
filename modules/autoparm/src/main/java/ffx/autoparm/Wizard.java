/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

//Java Imports
//import java.io.BufferedReader;
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
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
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
    protected String sdffile;
    private final static Logger LOGGER = Logger.getLogger(Wizard.class.getName());

    public Wizard(String sdffile) {
        this.sdffile = sdffile;
        LOGGER.setLevel(Level.FINEST);
    }

    String fullSmi = new String();
    String fullSmi_aro = new String();
    
    public String readSDF() throws FileNotFoundException, IOException {
        int records_read = 0;
        int records_processed = 0;
        int records_error = 0;

        File file = new File(sdffile);

        BufferedReader read = null;

        try {
            FileReader fileReader = new FileReader(file);
            read = new BufferedReader(fileReader);

        } catch (IOException e) {
        }

        IteratingSDFReader reader = null;

        try {
            reader = new IteratingSDFReader(read, SilentChemObjectBuilder.getInstance());
            LOGGER.log(Level.INFO, String.format("\nReading %s", file.getAbsoluteFile()));
            System.out.println();
            while (reader.hasNext()) {

                IAtomContainer molecule = reader.next();

                records_read++;
                try {
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
                    CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(molecule);

                    aromaticity.apply(molecule);
                    /**
                     * Generate SMILES and assign as properties
                     */
                    fullSmi_aro = assignSMILES(molecule);
                    //System.out.println("\nFinished assignSMILES\n");

                    /**
                     * Descriptor calculation
                     */
                    //calculateLogP(molecule);
                    //System.out.println("Finished calculateLogP");
                    /**
                     * Write as plain comma delimited
                     */
                    /*for (FIELD field : FIELD.values()) {
                        Object value = molecule.getProperty(field.name());
                        writer.write(value == null ? "" : value.toString());
                        writer.write(",");
                    }
                    writer.write("\n");
                    //records_processed++;;
                    System.out.println("Finished writing");*/
                } catch (Exception x) {
                    System.err.println("*");
                    records_error++;
                    LOGGER.log(Level.SEVERE, String.format("[Record %d] Error %s\n", records_read, file.getAbsoluteFile()), x);
                }
            }
        } catch (Exception x) {
            x.printStackTrace();
            LOGGER.log(Level.SEVERE, String.format("[Record %d] Error %s\n", records_read, file.getAbsoluteFile()), x);
        } finally {
            try {
                reader.close();
            } catch (Exception x) {
            }
            /*try {

                writer.close();
                System.out.println("Closed writer\n");

            } catch (Exception x) {
                x.printStackTrace();
            }*/
        }
        return fullSmi_aro;
    }

    private XLogPDescriptor xlogp;
    private ALOGPDescriptor alogp;

    protected String assignSMILES(IAtomContainer molecule) throws Exception {
        //molecule.setProperty(FIELD.SMILES_Kekule.name(), SmilesGenerator.isomeric().create(molecule));
        //molecule.setProperty(FIELD.SMILES_Aromatic.name(), new SmilesGenerator().aromatic().create(molecule));

        SmilesGenerator sg1 = SmilesGenerator.isomeric();
        String smi1 = sg1.create(molecule);
        SmilesGenerator sg2 = new SmilesGenerator().aromatic();
        String smi2 = sg2.create(molecule);

        System.out.println("SMILES_Kekule,SMILES_Aromatic");
        System.out.print(smi1);
        System.out.print(",");
        System.out.print(smi2);
        System.out.println();
        
        return smi2;
    }

    protected void calculateLogP(IAtomContainer molecule) throws Exception {
        if (xlogp == null) {
            xlogp = new XLogPDescriptor();
        }
        if (alogp == null) {
            alogp = new ALOGPDescriptor();
        }

        DescriptorValue value = xlogp.calculate(molecule);
        String[] names = value.getNames();
        for (String name : names) {
            molecule.setProperty(name, value.getValue().toString());
        }

        value = alogp.calculate(molecule);
        names = value.getNames();
        if (value.getValue() instanceof DoubleArrayResult) {
            for (int i = 0; i < names.length; i++) {
                molecule.setProperty(names[i], ((DoubleArrayResult) value.getValue()).get(i));
            }
        }
    }

}
