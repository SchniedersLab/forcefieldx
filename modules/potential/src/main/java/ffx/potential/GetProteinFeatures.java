package ffx.potential;

import ffx.potential.bonded.AminoAcidUtils;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Residue;

import java.util.HashMap;

public class GetProteinFeatures {
    private static final HashMap<AminoAcid3,String> polarityMap = new HashMap<>();
    private static final HashMap<AminoAcid3, String> acidityMap = new HashMap<>();
    private String[] features;


    static{
        // 1 is polar, 0 is nonpolar
        polarityMap.put(AminoAcid3.ARG, "polar");
        polarityMap.put(AminoAcid3.ASN, "polar");
        polarityMap.put(AminoAcid3.ASP, "polar");
        polarityMap.put(AminoAcid3.ASH, "polar");
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
        polarityMap.put(AminoAcid3.ALA, "nonpolar");
        polarityMap.put(AminoAcid3.CYS, "nonpolar");
        polarityMap.put(AminoAcid3.GLY, "nonpolar");
        polarityMap.put(AminoAcid3.ILE, "nonpolar");
        polarityMap.put(AminoAcid3.LEU, "nonpolar");
        polarityMap.put(AminoAcid3.MET, "nonpolar");
        polarityMap.put(AminoAcid3.PHE, "nonpolar");
        polarityMap.put(AminoAcid3.PRO, "nonpolar");
        polarityMap.put(AminoAcid3.TRP, "nonpolar");
        polarityMap.put(AminoAcid3.TYR, "nonpolar");
        polarityMap.put(AminoAcid3.VAL, "nonpolar");

        acidityMap.put(AminoAcid3.ASP, "+");
        acidityMap.put(AminoAcid3.ASH, "+");
        acidityMap.put(AminoAcid3.GLU, "+");
        acidityMap.put(AminoAcid3.GLH, "+");
        acidityMap.put(AminoAcid3.ARG, "-");
        acidityMap.put(AminoAcid3.HIS, "-");
        acidityMap.put(AminoAcid3.HIE, "-");
        acidityMap.put(AminoAcid3.HID, "-");
        acidityMap.put(AminoAcid3.LYS, "-");
        acidityMap.put(AminoAcid3.LYD, "-");

    }



    public void setFeatures(Residue res){
        String name = res.getName();
        AminoAcid3 aa3 = res.getAminoAcid3();
        String acid = acidityMap.getOrDefault(aa3, "null");
        features[0] = name;
        features[1] = String.valueOf(res.getResidueNumber());
        features[2] = polarityMap.get(aa3);
        features[3] = acid;
    }

    public String[] getFeatures(){
        System.out.println(features);
        return features;
    }
}
