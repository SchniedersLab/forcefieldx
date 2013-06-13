/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential;

import ffx.potential.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;
import ffx.potential.parameters.ForceField;
import ffx.potential.parsers.ForceFieldFilter;
import ffx.potential.parsers.PDBFilter;
import ffx.utilities.Keyword;
import java.io.File;
import java.util.ArrayList;
import org.apache.commons.configuration.CompositeConfiguration;
import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * @author William Tollefson
 */
public class RotamerLibraryTest {
    
    public RotamerLibraryTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of getRotamers method, of class RotamerLibrary.
     */
    @Test
    public void testGetRotamers() {
        ResidueEnumerations.AminoAcid3 name = AminoAcid3.ARG;
        int expResult = 14;
        Rotamer[] result = RotamerLibrary.getRotamers(name);
        assertEquals(expResult, result.length);
    }

    /**
     * Test of rotamerOptimization method, of class RotamerLibrary.
     */
    @Test
    public void testRotamerOptimization() {
        String filename = "ffx/potential/structures/peptide.pdb";
        ClassLoader cl = this.getClass().getClassLoader();
        File structure = new File(cl.getResource(filename).getPath());

        String name = structure.getName();
        int index = filename.lastIndexOf(".");
        name = filename.substring(0, index);
        MolecularAssembly molecularAssembly = new MolecularAssembly(name);
        molecularAssembly.setFile(structure);

        CompositeConfiguration properties = Keyword.loadProperties(structure);
        properties.setProperty("polarization", "direct");
        ForceFieldFilter forceFieldFilter = new ForceFieldFilter(properties);
        ForceField forceField = forceFieldFilter.parse();
        molecularAssembly.setForceField(forceField);
        PDBFilter pdbFilter = new PDBFilter(structure, molecularAssembly, forceField, properties);
        boolean read = pdbFilter.readFile();
        molecularAssembly.finalize(true);
        ForceFieldEnergy energy = new ForceFieldEnergy(molecularAssembly);
        molecularAssembly.setPotential(energy);
        assertEquals(read, true);
        ArrayList<Residue> residues = new ArrayList<Residue>();
        int permutations = 1;
        Polymer[] polymers = molecularAssembly.getChains();
        for (int i=14; i<=16; i++) {
            Residue residue = polymers[0].getResidue(i);
            residues.add(residue);
            ResidueEnumerations.AminoAcid3 aa = ResidueEnumerations.AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(aa);
            if (rotamers != null) {
                permutations *= rotamers.length;
            }
        }
        assertEquals(permutations, 108);
        ArrayList<Integer> optimum = new ArrayList<Integer>();
        double minEnergy = RotamerLibrary.rotamerOptimization(molecularAssembly, residues, Double.MAX_VALUE, optimum);
        for (int i=14; i<=16; i++) {
            Residue residue = polymers[0].getResidue(i);
            ResidueEnumerations.AminoAcid3 aa = ResidueEnumerations.AminoAcid3.valueOf(residue.getName());
            Rotamer[] rotamers = RotamerLibrary.getRotamers(aa);
            int j = optimum.remove(0);
            if (rotamers != null) {
                Rotamer rotamer = rotamers[j];
                RotamerLibrary.applyRotamer(aa, residue, rotamer);
            }
        }
        double finalEnergy = energy.energy(false, false);
        System.out.println(String.format(" Final energy %16.8f", finalEnergy));
        assertEquals(minEnergy, finalEnergy, 1.0e-8);
        assertEquals(finalEnergy, -586.514100, 1.0e-6);
    }
}