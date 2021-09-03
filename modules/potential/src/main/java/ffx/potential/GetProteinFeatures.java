package ffx.potential;

import ffx.potential.bonded.AminoAcidUtils;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.Torsion;
import org.openscience.cdk.AminoAcid;

import java.util.*;

public class GetProteinFeatures {
    private static final HashMap<AminoAcid3,String> polarityMap = new HashMap<>();
    private static final HashMap<AminoAcid3, String> acidityMap = new HashMap<>();
    private final String[] features = new String[5];
    private List<List<Torsion>> torsions;
    private int index;
    private static final NavigableMap<Double, String> phiToStructure = new TreeMap<>();
    private static final NavigableMap<Double, String> psiToStructure = new TreeMap<>();



    public GetProteinFeatures(List<List<Torsion>> torsions) {
        this.torsions = torsions;
    }



    static{
        //Map residue to its polarity
        polarityMap.put(AminoAcid3.ARG, "polar");
        polarityMap.put(AminoAcid3.ASN, "polar");
        polarityMap.put(AminoAcid3.ASP, "polar");
        polarityMap.put(AminoAcid3.ASH, "polar");
        polarityMap.put(AminoAcid3.CYS, "polar");
        polarityMap.put(AminoAcid3.GLN, "polar");
        polarityMap.put(AminoAcid3.GLU, "polar");
        polarityMap.put(AminoAcid3.GLH, "polar");
        polarityMap.put(AminoAcid3.HIS, "polar");
        polarityMap.put(AminoAcid3.HIE, "polar");
        polarityMap.put(AminoAcid3.HID, "polar");
        polarityMap.put(AminoAcid3.LYS, "polar");
        polarityMap.put(AminoAcid3.LYD, "polar");
        polarityMap.put(AminoAcid3.SER, "polar");
        polarityMap.put(AminoAcid3.THR, "polar");
        polarityMap.put(AminoAcid3.TYR, "polar");
        polarityMap.put(AminoAcid3.ALA, "nonpolar");
        polarityMap.put(AminoAcid3.GLY, "nonpolar");
        polarityMap.put(AminoAcid3.ILE, "nonpolar");
        polarityMap.put(AminoAcid3.LEU, "nonpolar");
        polarityMap.put(AminoAcid3.MET, "nonpolar");
        polarityMap.put(AminoAcid3.PHE, "nonpolar");
        polarityMap.put(AminoAcid3.PRO, "nonpolar");
        polarityMap.put(AminoAcid3.TRP, "nonpolar");
        polarityMap.put(AminoAcid3.VAL, "nonpolar");

        //Map residue to its acidity
        acidityMap.put(AminoAcid3.ASP, "acidic");
        acidityMap.put(AminoAcid3.ASN, "acidic");
        acidityMap.put(AminoAcid3.GLU, "acidic");
        acidityMap.put(AminoAcid3.GLN, "acidic");
        acidityMap.put(AminoAcid3.ARG, "basic");
        acidityMap.put(AminoAcid3.HIS, "basic");
        acidityMap.put(AminoAcid3.HIE, "basic");
        acidityMap.put(AminoAcid3.HID, "basic");
        acidityMap.put(AminoAcid3.LYS, "basic");
        acidityMap.put(AminoAcid3.LYD, "basic");
        acidityMap.put(AminoAcid3.GLY, "neutral");
        acidityMap.put(AminoAcid3.ALA, "neutral");
        acidityMap.put(AminoAcid3.VAL, "neutral");
        acidityMap.put(AminoAcid3.LEU, "neutral");
        acidityMap.put(AminoAcid3.ILE, "neutral");
        acidityMap.put(AminoAcid3.SER, "neutral");
        acidityMap.put(AminoAcid3.CYS, "neutral");
        acidityMap.put(AminoAcid3.THR, "neutral");
        acidityMap.put(AminoAcid3.MET, "neutral");
        acidityMap.put(AminoAcid3.PRO, "neutral");
        acidityMap.put(AminoAcid3.PHE, "neutral");
        acidityMap.put(AminoAcid3.TYR, "neutral");
        acidityMap.put(AminoAcid3.TRP, "neutral");

        //Map phi to the Alpha Helix (L) and (R) and Beta Sheet Range
        phiToStructure.put(-180.0, "Extended");
        phiToStructure.put(-150.0, "Structure");
        phiToStructure.put(-50.0, "Structure");
        phiToStructure.put(0.0, "Extended");
        phiToStructure.put(50.0, "Structure");
        phiToStructure.put(70.0, "Structure");
        phiToStructure.put(180.0, "Extended");

        //Map psi to the Alpha Helix (L) and (R) and Beta Sheet Range
        psiToStructure.put(-180.0, "Extended");
        psiToStructure.put(-70.0, "Alpha Helix (R)");
        psiToStructure.put(-50.0, "Alpha Helix (R)");
        psiToStructure.put(0.0, "Extended");
        psiToStructure.put(30.0, "Alpha Helix (L)");
        psiToStructure.put(70.0, "Alpha Helix (L)");
        psiToStructure.put(85.0, "Extended");
        psiToStructure.put(100.0, "Parallel Beta Sheet");
        psiToStructure.put(125.0, "Beta Sheet");
        psiToStructure.put(150.0, "Anti-Parallel Beta Sheet");
        psiToStructure.put(180.0, "Extended");
    }


    public void saveFeatures(Residue res){
        String name = res.getName();
        index = res.getResidueNumber() - 1;
        AminoAcid3 aa3 = res.getAminoAcid3();
        String acid = acidityMap.getOrDefault(aa3, "neutral");
        String structure = "";
        if (index == 0 || index == torsions.get(0).size() + 1){
            //Since phi, psi angles are determined between two residues, the first and last residue will not have values
            //and are default labeled as Extended Secondary Structure
            structure = "Extended";
        } else if (index > 0 && index <= torsions.get(0).size()) {
            structure = getSecondaryStructure();
        }
        features[0] = name;
        features[1] = String.valueOf(res.getResidueNumber());
        features[2] = polarityMap.get(aa3);
        features[3] = acid;
        features[4] = structure;
    }

    public String[] getFeatures(){
        return features;
    }

    //Data on phi,psi and secondary structure correlations is from Tamar Schlick's
    // Molecular Modeling and Simulation Textbook. Approx. page 108.
    //Alpha Helix: phi,psi of approximately -60,-50. Tighter helices are around -50,-25.
    //Looser helices are around -60,-70.

    //Collagen Helix: phi,psi of approximately -60,+125.

    //Parallel Beta Sheet: phi,psi of approximately -120,+115
    //Antiparallel Beta Sheet: phi,psi of approximately -140,+135
    public String getSecondaryStructure(){
        String secondaryStructure = "";
        String[] phiTorsion = torsions.get(0).get(index - 1).toString().split("\\s+");
        String[] psiTorsion = torsions.get(1).get(index - 1).toString().split("\\s+");
        Double phi = Double.parseDouble(phiTorsion[5].replace(",", ""));
        Double psi = Double.parseDouble(psiTorsion[5].replace(",", ""));
        //Use phi, psi ranges to determine which is the most likely secondary structure based on ramachandran values
        Double lowPhiKey = phiToStructure.floorKey(phi);
        String lowPhiStruct = phiToStructure.get(lowPhiKey);
        Double highPhiKey = phiToStructure.ceilingKey(phi);
        String highPhiStruct = phiToStructure.get(highPhiKey);
        if (lowPhiStruct.equals("Extended") || highPhiStruct.equals("Extended")){
            secondaryStructure = "Extended";
        } else {
            Double lowPsiKey = psiToStructure.floorKey(psi);
            String lowPsiStruct = psiToStructure.get(lowPsiKey);
            Double highPsiKey = psiToStructure.ceilingKey(psi);
            String highPsiStruct = psiToStructure.get(highPsiKey);
            if (lowPsiStruct.equals("Extended") || highPsiStruct.equals("Extended")){
                secondaryStructure = "Extended";
            } else if (lowPsiStruct.equals("Beta Sheet")){
                secondaryStructure = highPsiStruct;
            } else if (highPsiStruct.equals("Beta Sheet")){
                secondaryStructure = lowPsiStruct;
            } else {
                secondaryStructure = lowPsiStruct;
            }
        }
        return secondaryStructure;
    }


}
