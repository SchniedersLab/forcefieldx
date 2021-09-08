package ffx.potential;

import ffx.numerics.math.DoubleMath;
import ffx.potential.bonded.AminoAcidUtils.AminoAcid3;
import ffx.potential.bonded.Residue;

import java.util.*;

public class GetProteinFeatures {
    private static final HashMap<AminoAcid3,String> polarityMap = new HashMap<>();
    private static final HashMap<AminoAcid3, String> acidityMap = new HashMap<>();
    private final String[] features = new String[7];
    private static final NavigableMap<Double, String> phiToStructure = new TreeMap<>();
    private static final NavigableMap<Double, String> psiToStructure = new TreeMap<>();
    private double phi;
    private double psi;



    public GetProteinFeatures() {

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
        psiToStructure.put(-70.0, "Alpha Helix");
        psiToStructure.put(-50.0, "Alpha Helix");
        psiToStructure.put(0.0, "Extended");
        //psiToStructure.put(30.0, "Alpha Helix (L)");
        //psiToStructure.put(70.0, "Alpha Helix (L)");
        psiToStructure.put(85.0, "Extended");
        psiToStructure.put(100.0, "Beta Sheet");
        psiToStructure.put(150.0, "Beta Sheet");
        //psiToStructure.put(150.0, "Anti-Parallel Beta Sheet");
        psiToStructure.put(180.0, "Extended");
    }


    public void saveFeatures(Residue residue){
        String name = residue.getName();
        //index = res[0].getResidueNumber() - 1;
        AminoAcid3 aa3 = residue.getAminoAcid3();
        String acid = acidityMap.getOrDefault(aa3, "neutral");
        String structure = "";
        String phiString = "";
        String psiString = "";
        if (residue.getNextResidue() == null){
            //Since phi, psi angles are determined between two residues, the first and last residue will not have values
            //and are default labeled as Extended Secondary Structure
            structure = "Extended";
            getPhi(residue);
            phiString = String.valueOf(phi);
            psiString = null;
        } else if (residue.getPreviousResidue() == null) {
            structure = "Extended";
            psiString = String.valueOf(psi);
            phiString = null;
        } else {
            getPhi(residue);
            getPsi(residue);
            phiString = String.valueOf(phi);
            psiString = String.valueOf(psi);
            structure = getSecondaryStructure();
        }
        features[0] = name;
        features[1] = String.valueOf(residue.getResidueNumber());
        features[2] = polarityMap.get(aa3);
        features[3] = acid;
        features[4] = structure;
        features[5] = phiString;
        features[6] = psiString;
    }

    public String[] getFeatures(){
        return features;
    }


    public void getPhi(Residue currentRes){
        Residue previousRes = currentRes.getPreviousResidue();
        double[] cCoor = new double[3];
        double[] nCoor = new double[3];
        double[] caCoor = new double[3];
        double[] c2Coor = new double[3];
        phi = (DoubleMath.dihedralAngle(previousRes.getAtomByName("C", true).getXYZ(cCoor),
                    currentRes.getAtomByName("N", true).getXYZ(nCoor),
                    currentRes.getAtomByName("CA", true).getXYZ(caCoor),
                    currentRes.getAtomByName("C", true).getXYZ(c2Coor))) * 180/Math.PI;
    }

    public void getPsi(Residue currentRes){
        //res[0] is always current Res
        Residue nextRes = currentRes.getNextResidue();
        double[] nCoor = new double[3];
        double[] caCoor = new double[3];
        double[] cCoor = new double[3];
        double[] n2Coor = new double[3];
        psi = (DoubleMath.dihedralAngle(currentRes.getAtomByName("N", true).getXYZ(nCoor),
                currentRes.getAtomByName("CA", true).getXYZ(caCoor),
                currentRes.getAtomByName("C", true).getXYZ(cCoor),
                nextRes.getAtomByName("N", true).getXYZ(n2Coor))) * 180/Math.PI;
    }

    //Data on phi,psi and secondary structure correlations is from Tamar Schlick's
    // Molecular Modeling and Simulation Textbook. Approx. page 108.
    //Alpha Helix: phi,psi of approximately -60,-50. Tighter helices are around -50,-25.
    //Looser helices are around -60,-70.

    //Collagen Helix: phi,psi of approximately -60,+125.

    //Parallel Beta Sheet: phi,psi of approximately -120,+115
    //Antiparallel Beta Sheet: phi,psi of approximately -140,+135
    //Get surface area region to get surface area
    public String getSecondaryStructure(){
        String secondaryStructure = "";
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
        }        return secondaryStructure;
    }


}
