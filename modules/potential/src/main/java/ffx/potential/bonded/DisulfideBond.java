/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.bonded;

import java.util.logging.Logger;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.SSBond;

/**
 *
 * @author jmlitman
 */
public class DisulfideBond implements SSBond {
    
    private static final Logger logger = Logger.getLogger(DisulfideBond.class.getName());
    private Group cys1;
    private Group cys2;

    int serNum; // serial number
    String chainID1;
    String chainID2;
    String resnum1;
    String resnum2;
    String insCode1;
    String insCode2;
    
    public DisulfideBond() {
        serNum = 0;
    }
    
    public DisulfideBond(Group cys1, Group cys2) {
        if (cys1.getPDBName().equalsIgnoreCase("CYS")) {
            this.cys1 = cys1;
            chainID1 = cys1.getChainId();
            ResidueNumber rn1 = cys1.getResidueNumber();
            resnum1 = rn1.getSeqNum().toString();
            insCode1 = "" + rn1.getInsCode();
        } else {
            logger.warning(String.format(" Attempt to construct disulfide bond with non-cysteine residue %s", cys1.toString()));
        }
        if (cys2.getPDBName().equalsIgnoreCase("CYS")) {
            this.cys2 = cys2;
            chainID2 = cys2.getChainId();
            ResidueNumber rn2 = cys2.getResidueNumber();
            resnum2 = rn2.getSeqNum().toString();
            insCode2 = "" + rn2.getInsCode();
        } else {
            logger.warning(String.format(" Attempt to construct disulfide bond with non-cysteine residue %s", cys2.toString()));
        }
    }
    
    @Override
    public String toPDB() {

        StringBuffer buf = new StringBuffer();
        toPDB(buf);
        return buf.toString();
    }

    /**
     * append the PDB representation of this SSBOND to the provided StringBuffer
     *
     * @param buf a StringBuffer to print the PDB representation to
     */
    @Override
    public void toPDB(StringBuffer buf) {

        /*12 - 14        LString(3)      "CYS"        Residue name.
         16             Character       chainID1     Chain identifier.
         18 - 21        Integer         seqNum1      Residue sequence number.
         22             AChar           icode1       Insertion code.
         26 - 28        LString(3)      "CYS"        Residue name.
         30             Character       chainID2     Chain identifier.
         32 - 35        Integer         seqNum2      Residue sequence number.
         36             AChar           icode2       Insertion code.
         60 - 65        SymOP           sym1         Symmetry oper for 1st resid
         67 - 72        SymOP           sym2         Symmetry oper for 2nd resid
         */
		//01234567890123456789012345678901234567890123456789012345678901234567890123456789
        //SSBOND   1 CYS      5    CYS     55                                     5PTI  67
        //SSBOND   2 CYS     14    CYS     38
        //SSBOND   3 CYS     30    CYS     51
        buf.append("SSBOND ");
        buf.append(String.format("%3d", serNum));
        buf.append(String.format(" CYS %s %4s%1s  ", chainID1, resnum1, insCode1));
        buf.append(String.format(" CYS %s %4s%1s  ", chainID2, resnum2, insCode2));
    }

    @Override
    public String getInsCode1() {
        return insCode1;
    }

    @Override
    public void setInsCode1(String insCode1) {
        this.insCode1 = insCode1;
    }

    @Override
    public String getInsCode2() {
        return insCode2;
    }

    @Override
    public void setInsCode2(String insCode2) {
        this.insCode2 = insCode2;
    }

    /**
     * set serial number of this SSBOND in PDB file
     *
     * @return the serial number
     */
    @Override
    public int getSerNum() {
        return serNum;
    }

    /**
     * get serial number of this SSBOND in PDB file
     *
     * @param serNum
     */
    @Override
    public void setSerNum(int serNum) {
        this.serNum = serNum;
    }

    @Override
    public DisulfideBond clone() {
        DisulfideBond nbond = new DisulfideBond();
        nbond.setChainID1(chainID1);
        nbond.setChainID2(chainID2);
        nbond.setResnum1(resnum1);
        nbond.setResnum2(resnum2);
        nbond.setCys1(cys1);
        nbond.setCys2(cys2);
        nbond.setInsCode1(insCode1);
        nbond.setInsCode2(insCode2);
        nbond.setSerNum(serNum);
        return nbond;
    }

    @Override
    public String getChainID1() {
        return chainID1;
    }

    @Override
    public void setChainID1(String chainID1) {
        this.chainID1 = chainID1;
    }

    @Override
    public String getChainID2() {
        return chainID2;
    }

    @Override
    public void setChainID2(String chainID2) {
        this.chainID2 = chainID2;
    }

    /**
     * get residue number for first CYS. number and insertion code are joint
     * together.
     *
     * @return the residue number of the first CYS.
     *
     */
    @Override
    public String getResnum1() {
        return resnum1;
    }

    @Override
    public void setResnum1(String resnum1) {
        this.resnum1 = resnum1;
    }

    /**
     * get residue number for second CYS. number and insertion code are joint
     * together.
     *
     * @return the residue number of the second CYS.
     *
     */
    @Override
    public String getResnum2() {
        return resnum2;
    }

    @Override
    public void setResnum2(String resnum2) {
        this.resnum2 = resnum2;
    }
    
    public void setCys1(Group cys1) {
        this.cys1 = cys1;
    }
    
    public void setCys2(Group cys2) {
        this.cys2 = cys2;
    }
    
    public Group getCys1() {
        return cys1;
    }
    
    public Group getCys2() {
        return cys2;
    }

    @Override
    public String toString() {
        String s = "[SSBOND:\n";

        s += "Atom 1:\n";
        s += "\tChain: " + chainID1 + "\n";
        s += "\tResidue #: " + resnum1 + "\n";
        s += "\tIns. Code: " + insCode1 + "\n";

        s += "Atom 2:\n";
        s += "\tChain: " + chainID2 + "\n";
        s += "\tResidue #: " + resnum2 + "\n";
        s += "\tIns. Code: " + insCode2 + "\n";

        s += "]";

        return s;
    }
    
}
