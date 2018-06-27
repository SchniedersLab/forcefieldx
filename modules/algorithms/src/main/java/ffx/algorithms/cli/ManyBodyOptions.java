/**
 * Title: Force Field X.
 * <p>
 * Description: Force Field X - Software for Molecular Biophysics.
 * <p>
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
 * <p>
 * This file is part of Force Field X.
 * <p>
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 * <p>
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * <p>
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 * <p>
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 * <p>
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
package ffx.algorithms.cli;

import picocli.CommandLine.Option;
import ffx.algorithms.RotamerOptimization;
import java.io.File;
import java.io.IOException;

/**
 * Represents command line options for scripts that use a many-body expansion for global optimization.
 *
 * @author Michael J. Schnieders
 * @author Mallory R. Tollefson
 * @since 1.0
 */
public class ManyBodyOptions {

    /**
     * PARAMETERS SHARED BY ALL ALGORITHMS. 
     */
    
    /**
     * -L or --library Choose either Ponder and Richards (1) or Richardson (2) rotamer library.
     */
    @Option(names = {"-L", "--library"}, paramLabel = "1",
            description = "Ponder and Richards (1) or Richardson (2) rotamer library.")
    int library = 1;
    
    /**
     * -a or --algorithm Choices are independent residues (1), all with rotamer elimination (2),
     * all brute force (3), sliding window (4), or box optimization (5).
     */
    @Option(names = {"-a", "--algorithm"}, paramLabel = "2",
            description = "Algorithm: independent residues (1), all with rotamer elimination (2), all brute force (3), sliding window (4), or box optimization (5)")
    int algorithm = 2;
    
    /**
     * --dee or --deadEnd Use dead-end elimination criteria instead of Goldstein criteria.
     */
    @Option(names = {"--dee", "--deadEnd"}, paramLabel = "false",
            description = "Use dead-end elimination criteria instead of Goldstein criteria.")
    boolean dee = false;
   
    /**
     * -c or --chain Single character chain ID of the residues to optimize.
     */
    @Option(names = {"-c", "--chain"}, paramLabel = "-1",
            description = "Single character chain ID of the residues to optimize.")
    String chain = "-1";
    
    /**
     * -s or --start Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.
     */
    @Option(names = {"-s", "--start"}, paramLabel = "-1",
            description = "Starting residue to perform the optimization on (-1 exits). For box optimization, first box to optimize.")
    int start = -1;
    
    /**
     * -f or --final Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.
     */
    @Option(names = {"-f", "--final"}, paramLabel = "-1",
            description = "Final residue to perform the optimization on (-1 exits). For box optimization, final box to optimize.")
    int finish = -1;
    
    /**
     * --tC or --twoBodyCutoff Cutoff distance for two-body interactions.
     */
    @Option(names = {"--tC", "--twoBodyCutoff"}, paramLabel = "-1.0",
            description = "Cutoff distance for two body interactions.")
    double twoBodyCutoff = -1.0;
    
    /**
     * -t or --threeBody Include 3-Body interactions in the elimination criteria.
     */
    @Option(names = {"-t", "--threeBody"}, 
            description = "Include 3-Body interactions in the elimination criteria.")
    boolean threeBody = false;
    
    /**
     * --thC or --threeBodyCutoff Cutoff distance for three-body interactions.
     */
    @Option(names = {"--thC", "--threeBodyCutoff"}, paramLabel = "9.0",
            description = "Cutoff distance for three-body interactions.")
    double threeBodyCutoff = 9.0;
         
    /**
     * -p or --prune Prune no clashes (0), only single clashes (1), or all clashes (2).
     */
    @Option(names = {"-p", "--prune"}, paramLabel = "2",
            description = "Prune no clashes (0), only single clashes (1), or all clashes (2)")
    int prune = 2;
       
    /**
     * -x or --all Optimize all residues beginning from the passed value (overrides other options);
     * for box optimization, optimizes all boxes beginning from the passed index.
     */
    @Option(names = {"-x", "--all"}, paramLabel = "-1",
            description = "Optimize all residues beginning from the passed value (overrides other options); for box optimization, optimizes all boxes beginning from the passed index.")
    int all = -1;
        
    /**
     * --eR or --energyRestart Load energy restart file from a previous run (requires that all parameters are the same).
     */
    @Option(names = {"--eR", "--energyRestart"}, paramLabel = "none",
            description = "Load energy restart file from a previous run (requires that all parameters are the same).")
    String energyRestart = "none";       
                
    /**
     * -m or --minimize Minimize the final structure to the given RMS gradient (Kcal/mole/A); the default is not
     * to minimize.
     */
    @Option(names = {"-m", "--minimize"}, paramLabel = "0.0",
            description = "Minimize the final structure to the given RMS gradient (Kcal/mole/A).")
    double minimize = 0.0;             

    /**
     * -v or --verbose Prints beginning and default-conformation energies.
     */
    @Option(names = {"-v", "--verbose"}, 
            description = "Prints beginning and default-conformation energies.")
    boolean verbose = false;     

    /**
     * -o or --original Do not include starting coordinates as their own rotamer.
     */
    @Option(names = {"-o", "--original"}, 
            description = "Do not include starting coordinates as their own rotamer.")
    boolean original = true;    

    /**
     * -d or --decompose Print energy decomposition for original-coordinates rotamers.
     */
    @Option(names = {"-d", "--decompose"}, 
            description = "Print energy decomposition for original-coordinates rotamers.")
    boolean decompose = false;    

    /**
     * --lR or --listResidues Choose a list of individual residues to optimize (eg. A11,A24,B40).
     */
    @Option(names = {"--lR", "--listResidues"}, paramLabel = "none", 
            description = "Choose a list of individual residues to optimize (eg. A11,A24,B40).")
    String listResidues = "none"; 
        
    /**
     * --sO or --sequence Choose a list of individual residues to sequence optimize (example: A2.A3.A5).
     */
    @Option(names = {"--sO", "--sequence"}, paramLabel = "none", 
            description = "Choose a list of individual residues to sequence optimize (example: A2.A3.A5)")
    String sequence = "none"; 

    /**
     * --tO or --titrationOptimization Optimize the titration states for a list of residues (example: H2.H3.H5).
     */
    @Option(names = {"--tO", "--titrationOptimization"}, paramLabel = "none", 
            description = "Optimize the titration states for a list of residues (example: H2.H3.H5).")
    String titrationOptimization = "none"; 
    
    /**
     * --mC or --monteCarlo Follow elimination criteria with 'n' Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.
     */
    @Option(names = {"--mC", "--monteCarlo"}, paramLabel = "-1", 
            description = "Follow elimination criteria with (n) Monte Carlo steps, or enumerate all remaining conformations, whichever is smaller.")
    int monteCarlo = -1; 
        
    /**
     * -z or --noRevert Do not revert an unfavorable change.
     */
    @Option(names = {"-z", "--revert"}, 
            description = "Undo an unfavorable change.")
    boolean noRevert = true;     
    
    /**
     * --fi or --file Choose the file type to write [PDB/XYZ].
     */
    @Option(names = {"--fi", "--file"}, paramLabel = "XYZ",
            description = "Choose file type to write [PDB/XYZ].")
    String fileType = "XYZ";

    /**
     * -out or --output Save eliminated singles and eliminated pairs to a text file (global and box optimization).
     */
    @Option(names = {"--out", "--output"}, paramLabel = "none",
            description = "Save eliminated singles and eliminated pairs to a text file.")
    boolean saveOutput = false;


    /**
     * PARAMETERS SPECIFIC TO SLIDING WINDOW.
     */
    
    /**
     * -w or --window Size of the sliding window with respect to adjacent residues (default = 7).
     */
    @Option(names = {"-w", "--window"}, paramLabel = "7",
            description = "Size of the sliding window with respect to adjacent residues.")
    int window = 7;
    
    /**
     * -i or --increment Sliding window increment (default = 3).
     */
    @Option(names = {"-I", "--INCREMENT"}, paramLabel = "3",
            description = "Sliding window increment.")
    int increment = 3;
    
    /**
     * -r or --cutoff The sliding window cutoff radius (Angstroms).
     */
    @Option(names = {"--r", "--cutoff"}, paramLabel = "2.0",
            description = "The sliding window cutoff radius (Angstroms).")
    double cutoff = 2.0;
    
    /**
     * -fR or --forceResidues Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.
     */
    @Option(names = {"--fR", "--forceResidues"}, paramLabel = "-1,-1",
            description = "Force residues in this range to be considered for sliding window radii, regardless of whether they lack rotamers.")
    String forceResidues = "-1,-1";
    
    
    
    /**
     * PARAMETERS SPECIFIC TO SLIDING BOX.
     */
    
    /**
     * -nB or --numBoxes Specify number of boxes along X, Y, and Z (default: '3,3,3').
     */
    @Option(names = {"--nB", "--numBoxes"}, paramLabel = "3,3,3",
            description = "Specify number of boxes along X, Y, and Z (default: 3,3,3)")
    String numBoxes = "3,3,3";
    
    /**
     * -bB or --boxBorderSize Extent of overlap between optimization boxes in Angstroms (default: 3.0).
     */
    @Option(names = {"--bB", "--boxBorderSize"}, paramLabel = "3.0",
            description = "Extent of overlap between optimization boxes in Angstroms.")
    double boxBorderSize = 3.0;
        
    /**
     * -bL or --approxBoxLength Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).
     * Box sizes are rounded up to make a whole number of boxes along each axis (default of 0 disables this function).
     */
    @Option(names = {"--bL", "--approxBoxLength"}, paramLabel = "0.0",
            description = "Approximate side lengths of boxes to be constructed (over-rides numXYZBoxes).")
    double approxBoxLength = 0.0;

    /**
     * -bC or --boxInclusionCriterion Criterion to use for adding residues to boxes.
     *      (1) uses C alpha only (N1/9 for nucleic acids)
     *      (2) uses any atom.
     *      (3) uses any rotamer
     */
    @Option(names = {"--bC", "--boxInclusionCriterion"}, paramLabel = "1",
            description = "Criterion to use for adding a residue to a box: (1) uses C alpha only (N1/9 for nucleic acids), (2) uses any atom, and (3) uses any rotamer")
    int boxInclusionCriterion = 1;    

        
 
    /**
     * PARAMETERS SPECIFIC TO NUCLEIC OPTIMIZATION.
     */
         
    /**
     * -nT or --nucleicCorrectionThreshold Nucleic acid rotamers adjusted by more than a threshold distance (A) are discarded (0 disables this function).
     */
    @Option(names = {"--nT", "--nucleicCorrectionThreshold"}, paramLabel = "0",
            description = "Nucleic acid rotamers adjusted by more than a threshold distance (A) are discarded (0 disables this function).")
    double nucleicCorrectionThreshold = 0;
    
    /**
     * -mN or --minimumAcceptedNARotamers Minimum number of NA rotamers to be accepted if a threshold distance is enabled.
     */
    @Option(names = {"--mN", "--minimumAcceptedNARotamers"}, paramLabel = "10",
            description = "Minimum number of NA rotamers to be accepted if a threshold distance is enabled.")
    int minimumAcceptedNARotamers = 10;
    
    RotamerOptimization rotamerOptimization;
    
    public void setRotamerOptimization(RotamerOptimization rotamerOptimization) {
        this.rotamerOptimization = rotamerOptimization;
    }
    
    public void setRotOptProperties(int [] numXYZBoxes, int forceResiduesStart, int forceResiduesEnd) {
        rotamerOptimization.setTwoBodyCutoff(twoBodyCutoff);
        rotamerOptimization.setThreeBodyCutoff(threeBodyCutoff);
        rotamerOptimization.setThreeBodyEnergy(threeBody);
        rotamerOptimization.setUseGoldstein(!dee);
        rotamerOptimization.setRevert(!noRevert);
        rotamerOptimization.setPruning(prune);
        rotamerOptimization.setWindowSize(window);
        rotamerOptimization.setIncrement(increment);
        rotamerOptimization.setDistanceCutoff(cutoff);
        rotamerOptimization.setNucleicCorrectionThreshold(nucleicCorrectionThreshold);
        rotamerOptimization.setMinimumNumberAcceptedNARotamers(minimumAcceptedNARotamers);
        rotamerOptimization.setVerboseEnergies(verbose);
        rotamerOptimization.setBoxBorderSize(boxBorderSize);
        rotamerOptimization.setApproxBoxLength(approxBoxLength);
        rotamerOptimization.setNumXYZBoxes(numXYZBoxes);
        rotamerOptimization.setBoxInclusionCriterion(boxInclusionCriterion);
        rotamerOptimization.setForcedResidues(forceResiduesStart, forceResiduesEnd);
        
        boolean monteCarloBool = false;
        if (monteCarlo > 1) {
            monteCarloBool = true;
        }
        rotamerOptimization.setMonteCarlo(monteCarloBool, monteCarlo);
        
        File energyRestartFile = null;
        if (!energyRestart.equalsIgnoreCase("none")) {
            energyRestartFile = new File(energyRestart);
            rotamerOptimization.setEnergyRestartFile(energyRestartFile);
        }
    }
    
    public void saveEliminatedRotamers () throws IOException {
        if (saveOutput){
            rotamerOptimization.outputEliminated();
        }
    }
}
