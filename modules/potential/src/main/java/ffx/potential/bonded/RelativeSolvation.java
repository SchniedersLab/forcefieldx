/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.bonded;

/**
 *
 * @author JacobLitman
 */
public class RelativeSolvation {
    
    private SolvationLibrary solvationLibrary;
    
    public RelativeSolvation() {
        this(SolvationLibrary.GK);
    }
    
    public RelativeSolvation(SolvationLibrary solvationLibrary) {
        this.solvationLibrary = solvationLibrary;
    }
    
    public void setSolvationLibrary(SolvationLibrary solvationLibrary) {
        this.solvationLibrary = solvationLibrary;
    }
    
    public SolvationLibrary getSolvationLibrary() {
        return solvationLibrary;
    }
    
    /**
     * Gets the solvation energy (desolvation penalty) for a given residue, allowing
     * for sequence optimization to include an estimate of energy relative to the 
     * unfolded state.
     * @param residue Residue to check
     * @param checkZeroes Throws an error if not in solvation energy library
     * @return Solvation energy
     */
    public double getSolvationEnergy(Residue residue, boolean checkZeroes) throws IllegalArgumentException {
        String resName = "";
        double energy = 0;
        switch (residue.getResidueType()) {
            case AA:
                if (residue instanceof MultiResidue) {
                    resName = ((MultiResidue) residue).getActive().getName();
                } else {
                    resName = residue.getName();
                }
                energy = getSolvationEnergy(ResidueEnumerations.AminoAcid3.valueOf(resName));
                break;
            case NA:
                if (residue instanceof MultiResidue) {
                    resName = ((MultiResidue) residue).getActive().getName();
                } else {
                    resName = residue.getName();
                }
                energy = getSolvationEnergy(ResidueEnumerations.NucleicAcid3.valueOf(resName));
                break;
            default:
                energy = 0;
        }
        if (checkZeroes && energy == 0) {
            throw new IllegalArgumentException(String.format(" Zero desolvation "
                    + "energy for residue %s: likely not in solvation library.", resName));
        }
        return energy;
    }
    
    /**
     * Returns amino acid solvation energy based on solvation library.
     * @param name Amino acid name
     * @return Solvation energy
     */
    public double getSolvationEnergy(ResidueEnumerations.AminoAcid3 name) {
        switch (solvationLibrary) {
            case WOLFENDEN:
                return getWolfendenSolvationEnergy(name);
            case CABANI:
                return getCabaniSolvationEnergy(name);
            case EXPLICIT:
                return getExplicitSolvationEnergy(name);
            case GK:
                return getGKSolvationEnergy(name);
            case MACCALLUM_SPC:
                return getMacCallumSPCSolvationEnergy(name);
            case MACCALLUM_TIP4P:
                return getMacCallumTIP4PSolvationEnergy(name);
            default:
                return 0;
        }
    }
    
    /**
     * Will return relative solvation energies for nucleic acids; currently returns
     * 0.
     * @param name Nucleic acid name
     * @return Relative solvation energy
     */
    public double getSolvationEnergy(ResidueEnumerations.NucleicAcid3 name) {
        return 0;
    }
    
    /**
     * Will return solvation energies relative to glycine for capped monomers
     * in GK solvent; currently wraps getExplicitSolvationEnergy.
     * @param name Amino acid name
     * @return Relative solvation energy
     */
    public double getGKSolvationEnergy(ResidueEnumerations.AminoAcid3 name) {
        return getExplicitSolvationEnergy(name);
    }
    
    /**
     * Will return solvation energies relative to glycine for capped monomers in
     * AMOEBA solvent; currently approximates this with charging energies in 
     * AMOEBA solvent.
     * @param name Amino acid name
     * @return Relative solvation energy
     */
    public double getExplicitSolvationEnergy(ResidueEnumerations.AminoAcid3 name) {
        switch (name) {
            case ALA:
                return 0.58;
            case CYS:
                return -0.85;
            case CYD:
                return -69.82;
            case ASP:
                return -69.45;
            case ASH:
                return -4.00;
            case GLU:
                return -71.40;
            case GLH:
                return -3.61;
            case PHE:
                return -1.59;
            case GLY:
                return 0.67;
            case HIS:
                return -45.36;
            case HID:
                return -8.42;
            case HIE:
                return -7.53;
            case ILE:
                return 0.14;
            case LYS:
                return -43.98;
            case LYD:
                return +0.35;
            case MET:
                return -3.48;
            case ASN:
                return -5.89;
            case PRO:
                return 7.82;
            case GLN:
                return -6.89;
            case ARG:
                return -42.57;
            case SER:
                return -2.14;
            case THR:
                return 0.58;
            case VAL:
                return 0.10;
            case TRP:
                return -4.64;
            case TYR:
                return 1.76;
            case TYD:
                return -41.71;
            default:
                return 0;
        }
    }
    
    /**
     * Returns absolute solvation energies for side chain analogs as calculated
     * by MacCallum for OPLS in TIP4P solvent.
     * 
     * @param name Amino acid name
     * @return Solvation energy
     */
    public double getMacCallumTIP4PSolvationEnergy(ResidueEnumerations.AminoAcid3 name) {
        switch (name) {
            case ALA:
                return 9.8;
            case CYS:
                return -0.5;
            case ASP:
                return -30.5;
            case GLU:
                return -19.0;
            case PHE:
                return -1.2;
            case HIS:
                return -28.0;
            case ILE:
                return 12.2;
            case LYS:
                return -13.6;
            case LEU:
                return 13.7;
            case MET:
                return -7.1;
            case ASN:
                return -34.5;
            case GLN:
                return -31.4;
            case ARG:
                return -43.9;
            case SER:
                return -20.2;
            case THR:
                return -20.3;
            case VAL:
                return 12.0;
            case TRP:
                return -16.2;
            case TYR:
                return -18.8;
            case GLY:
            case PRO:
            default:
                return 0;
        }
    }
    
    /**
     * Returns absolute solvation energies for side chain analogs as calculated
     * by MacCallum for OPLS in SPC solvent.
     * 
     * @param name Amino acid name
     * @return Solvation energy
     */
    public double getMacCallumSPCSolvationEnergy(ResidueEnumerations.AminoAcid3 name) {
        switch (name) {
            case ALA:
                return 9.3;
            case CYS:
                return -1.1;
            case ASP:
                return -30.1;
            case GLU:
                return -18.8;
            case PHE:
                return -1.4;
            case HIS:
                return -27.2;
            case ILE:
                return 11.9;
            case LYS:
                return -8.6;
            case LEU:
                return 12.6;
            case MET:
                return -5.1;
            case ASN:
                return -34.3;
            case GLN:
                return -30.8;
            case ARG:
                return -46.3;
            case SER:
                return -18.5;
            case THR:
                return -19.3;
            case VAL:
                return 11.3;
            case TRP:
                return -15.1;
            case TYR:
                return -18.2;
            case GLY:
            case PRO:
            default:
                return 0;
        }
    }
    
    /**
     * Returns absolute solvation energies for side chain analogs as experimentally
     * measured by Cabani et al.
     * 
     * @param name Amino acid name
     * @return Solvation energy
     */
    public double getCabaniSolvationEnergy(ResidueEnumerations.AminoAcid3 name) {
        switch (name) {
            case ALA:
                return 8.4;
            case CYS:
                return -5.2;
            case ASP:
                return -28.1;
            case GLU:
                return -27.1;
            case PHE:
                return -3.7;
            case HIS:
                return -27.4;
            case ILE:
                return 8.7;
            case LYS:
                return -15.5;
            case LEU:
                return 9.7;
            case MET:
                return 9.0;
            case ASN:
                return -40.6;
            case GLN:
                return -18.7;
            case ARG:
                return -30.1;
            case SER:
                return -21.4;
            case THR:
                return -21.0;
            case VAL:
                return 8.2;
            case TRP:
                return -12.3;
            case TYR:
                return -25.7;
            case GLY:
            case PRO:
            default:
                return 0;
        }
    }
    
    /**
     * Returns absolute solvation energies for side chain analogs as experimentally
     * measured by Wolfenden et al.
     * 
     * @param name Amino acid name
     * @return Solvation energy
     */
    public double getWolfendenSolvationEnergy(ResidueEnumerations.AminoAcid3 name) {
        switch (name) {
            case ALA:
                return 8.1;
            case CYS:
                return -5.1;
            case ASP:
                return -27.5;
            case GLU:
                return -26.6;
            case PHE:
                return -3.1;
            case HIS:
                return -42.1;
            case ILE:
                return 8.8;
            case LYS:
                return -18.0;
            case LEU:
                return 9.4;
            case MET:
                return -6.1;
            case ASN:
                return -39.9;
            case GLN:
                return -38.7;
            case ARG:
                return -44.8;
            case SER:
                return -20.8;
            case THR:
                return -20.1;
            case VAL:
                return 8.2;
            case TRP:
                return -24.3;
            case TYR:
                return -25.2;
            case GLY:
            case PRO:
            default:
                return 0; // Not listed.
        }
    }
    
    public String toString() {
        return "Relative solvation library: " + solvationLibrary.toString();
    }
    
    public enum SolvationLibrary {
        WOLFENDEN, CABANI, EXPLICIT, GK, MACCALLUM_SPC, MACCALLUM_TIP4P, OPLS_EXPLICIT,
        OPLS_GK
        /**
         * Citations:
         * Wolfenden et al:
         * Wolfenden, R., Andersson, L., Cullis, P. M. and Southgate, C. C. B. (1981) AFFINITIES OF AMINO-ACID SIDE-CHAINS FOR SOLVENT WATER. Biochemistry. 20, 849-855
         * 
         * Cabani et al:
         * Cabani, S., Gianni, P., Mollica, V. and Lepori, L. (1981) GROUP CONTRIBUTIONS TO THE THERMODYNAMIC PROPERTIES OF NON-IONIC ORGANIC SOLUTES IN DILUTE AQUEOUS-SOLUTION. Journal of Solution Chemistry. 10, 563-595
         * 
         * MacCallum OPLS libraries:
         * Maccallum, J. L. and Tieleman, D. P. (2003) Calculation of the water-cyclohexane transfer free energies of neutral amino acid side-chain analogs using the OPLS all-atom force field. Journal of Computational Chemistry. 24, 1930-1935
         */
    }
}
