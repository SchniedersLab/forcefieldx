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
    private final String[] features = new String[4];
    private List<List<Torsion>> torsions;
    private int index;


    public GetProteinFeatures(List<List<Torsion>> torsions) {
        this.torsions = torsions;
    }

    //Data on phi,psi and secondary structure correlations is from Tamar Schlick's
    // Molecular Modeling and Simulation Textbook. Approx. page 108.
    //Alpha Helix: phi,psi of approximately -60,-50. Tighter helices are around -50,-25.
    //Looser helices are around -60,-70.

    //Collagen Helix: phi,psi of approximately -60,+125.

    //Parallel Beta Sheet: phi,psi of approximately -120,+115
    //Antiparallel Beta Sheet: phi,psi of approximately -140,+135

    static{
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

    }


    public void saveFeatures(Residue res, int i){
        String name = res.getName();
        index = i;
        AminoAcid3 aa3 = res.getAminoAcid3();
        String acid = acidityMap.getOrDefault(aa3, "neutral");
        if (index < torsions.get(0).size() ) {
            getSecondaryStructure();
        }
        features[0] = name;
        features[1] = String.valueOf(res.getResidueNumber());
        features[2] = polarityMap.get(aa3);
        features[3] = acid;
    }

    public String[] getFeatures(){
        return features;
    }


    public String getSecondaryStructure(){
        String secondaryStructure = "";
        String[] phiTorsion = torsions.get(0).get(index).toString().split("\\s+");
        String[] psiTorsion = torsions.get(1).get(index).toString().split("\\s+");
        String phi = phiTorsion[5].replace(",", "");
        String psi = psiTorsion[5].replace(",", "");


        return secondaryStructure;
    }


}
