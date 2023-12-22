package ffx.algorithms.dynamics;

import edu.rit.mp.DoubleBuf;
import edu.rit.pj.Comm;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import ffx.algorithms.dynamics.integrators.IntegratorEnum;
import ffx.algorithms.dynamics.thermostats.ThermostatEnum;
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsUtils;
import ffx.potential.utils.Superpose;
import org.apache.commons.configuration2.CompositeConfiguration;

/**
 * This class implements the Weighted Ensemble algorithm. Simulations are ran in parallel with weights
 * assigned evenly to all simulations at the start. After a specified number of stochastic dynamics,
 * the parallel configurations are resampled (merged or split) and their weights are updated accordingly.
 *
 * Resampling is done by assigning each configuration a bin number. The bin number is determined by some
 * predetermined metric (e.g. RMSD, distance between residues, distance between complexes). The configurations
 * are then sorted by bin number and the configurations with the lowest weights are merged with other
 * configurations in the same bin. The configurations with the highest weights are split into two
 * configurations with half the weight of the parent.
 */
public class WeightedEnsemble {
    private static final Logger logger = Logger.getLogger(WeightedEnsemble.class.getName());


    private final int rank, worldSize;
    private final double[][] weightsBins, weightsBinsCopyRank, globalValues; // Local storage of all weights & bin numbers
    private final double[] myWeightBin, myGlobalValue;
    private final Comm world;
    private final DoubleBuf[] weightsBinsBuf, globalValuesBuf; // World storage of weights & bin numbers
    private final DoubleBuf myWeightBinBuf, myGlobalValueBuf;

    private boolean restart, staticBins;
    private boolean[] enteredNewBins;
    private int numBins, optNumPerBin;
    private long totalSteps, numStepsPerResample;
    private double weight, dt, temp, trajInterval;
    private double[] refCoords, binBounds, x;
    private File dynFile, dynTemp, trajFile, trajTemp, currentDir, refStructureFile;
    private MolecularDynamics molecularDynamics;
    private Potential potential;
    private OneDimMetric metric;


    public static enum OneDimMetric {RMSD, RESIDUE_DISTANCE, COM_DISTANCE, ATOM_DISTANCE, POTENTIAL_MIN, POTENTIAL_MAX};


    public WeightedEnsemble(OneDimMetric metric, int optNumPerBin,
                            MolecularDynamics md,
                            File refStructureFile) {
        // Init MPI stuff
        this.world = Comm.world();
        this.rank = world.rank();
        this.worldSize = world.size();
        /* Check if world size is valid
        if (worldSize < 2){
            logger.severe(" Weighted Ensemble requires at least 2 ranks.");
            System.exit(1);
        } else if (worldSize < 10){
            logger.warning(" Weighted Ensemble is not recommended for less than 10 parallel simulations.");
        }
         */

        // Init FFX stuff
        this.refStructureFile = refStructureFile;
        getRefCoords(refStructureFile, refCoords);
        this.dynFile = md.getDynFile();
        this.molecularDynamics = md;
        this.potential = null;
        if(molecularDynamics.molecularAssembly != null){
            this.potential = molecularDynamics.molecularAssembly[0].getPotentialEnergy();
        } else {
            logger.severe(" Molecular Assembly not set for Molecular Dynamics.");
            System.exit(1);
        }
        if (initFilesOrPrepRestart()){
            if (restart) {
                logger.info(" Restarting from previous Weighted Ensemble run.");
            } else {
                logger.info(" Initializing Weighted Ensemble run.");

            }
        } else {
            logger.severe(" Failed to initialize Weighted Ensemble run.");
            System.exit(1);
        }

        // Init Weighted Ensemble stuff
        CompositeConfiguration properties = molecularDynamics.molecularAssembly[0].getProperties();
        this.binBounds = getPropertyList(properties, "WE.BinBounds");
        this.numBins = binBounds.length + 1;
        this.staticBins = binBounds.length > 0;
        this.numBins = staticBins ? numBins : 500; // TODO: Implement dynamic binning & fix this hack for buffer size
        this.optNumPerBin = optNumPerBin;
        this.metric = metric;
        this.weight = 1.0 / worldSize;
        this.enteredNewBins = new boolean[worldSize];
        this.globalValues = new double[worldSize][1]; // [[globalValue (i.e. rmsd)] for each rank]
        this.weightsBins = new double[worldSize][2]; // [[weight, bin] for each rank]
        this.weightsBinsCopyRank = new double[worldSize][3]; // [[weight, bin, rankToCopyFrom] for each rank]
        this.weightsBinsBuf = new DoubleBuf[worldSize];
        this.globalValuesBuf = new DoubleBuf[worldSize];
        for (int i = 0; i < worldSize; i++){
            weightsBinsBuf[i] = DoubleBuf.buffer(weightsBins[i]);
            globalValuesBuf[i] = DoubleBuf.buffer(globalValues[i]);
        }
        this.myWeightBinBuf = weightsBinsBuf[rank];
        this.myGlobalValueBuf = globalValuesBuf[rank];
        myWeightBin = weightsBins[rank];
        myWeightBin[0] = weight;
        myGlobalValue = globalValues[rank];

        this.staticBins = false;
    }

    private double[] getPropertyList(CompositeConfiguration properties, String propertyName) {
        ArrayList<Double> list = new ArrayList<>();
        String[] split = properties.getString(propertyName, "").trim()
                .replace("[", "")
                .replace("]","")
                .replace(","," ")
                .split(" ");
        if (split[0].isEmpty()){
            return new double[0];
        }
        for (String s1 : split) {
            if (s1.isEmpty()) {
                continue;
            }
            list.add(Double.parseDouble(s1));
        }
        return list.stream().sorted().mapToDouble(Double::doubleValue).toArray();
    }

    private static void getRefCoords(File refStructureFile, double[] refCoords){
        PotentialsUtils utils = new PotentialsUtils();
        MolecularAssembly assembly = utils.open(refStructureFile);
        assembly.getPotentialEnergy().getCoordinates(refCoords);
    }

    private boolean initFilesOrPrepRestart() {
        // My directory info
        File parent = refStructureFile.getParentFile();
        currentDir = new File(parent + File.separator + rank);
        trajFile = new File(currentDir + File.separator + refStructureFile.getName() + ".arc");
        molecularDynamics.setArchiveFiles(new File[]{trajFile});
        dynFile = new File(currentDir + File.separator + refStructureFile.getName() + ".dyn");
        molecularDynamics.setFallbackDynFile(dynFile);
        restart = !currentDir.mkdir() && checkRestartFiles();
        molecularDynamics.setAutomaticWriteouts(false);
        return false;
    }

    private boolean checkRestartFiles(){
        // Check if restart files exist
        Path dynPath = dynFile.toPath();
        Path trajPath = trajFile.toPath();
        // TODO: Restart in place of where previous run left off instead of starting over
        // This could potentially be done by reading the number of traj frames
        return dynPath.toFile().exists() && trajPath.toFile().exists();
    }

    /**
     * Run the Weighted Ensemble algorithm.
     */
    public void run(long totalSteps, long numStepsPerResample, double temp, double dt) {
        this.totalSteps = totalSteps;
        this.numStepsPerResample = numStepsPerResample;
        this.temp = temp;
        this.dt = dt;
        if (!staticBins) {
            dynamics(numStepsPerResample * 3L);
        }

        int numCycles = (int) Math.ceil((double) totalSteps / numStepsPerResample);
        for (int i = 0; i < numCycles; i++) {
            dynamics(numStepsPerResample);
            molecularDynamics.writeRestart();
            calculateMyMetric();
            comms(); // Separate comms for dynamic binning
            binAssignment();
            comms();
            resample();
            commsSanityCheckAndFileMigration();
            logger.info(" Resampling cycle " + i + " complete.");
        }
    }

    private void dynamics(long numSteps){
        // Convert numStepsPerResample to picoseconds using dt (in femtoseconds)
        double dynamicTotalTime =  (long) (numSteps * dt / 1000.0);
        molecularDynamics.setCoordinates(x);
        potential.energy(x, false);
        molecularDynamics.dynamic(numSteps, dt, dynamicTotalTime/3.0,
                dynamicTotalTime/3.0, temp, true, dynFile);
        x = molecularDynamics.getCoordinates();
    }

    /**
     * Combine/merge configurations with the same bin number
     * <p>
     * I have to look at what everyone has to do --> merge(-1) / split(1) / nothing(0)
     * <p>
     * 'Ideal' particle weight is P/N where P is the total prob & N is number of particles in the bin
     * <p>
     * Merge: I have a low weight, my bin has > optNumPerBin
     *   | If I merge, I need to find something else to pick up
     *   | 'a' & 'b' are merged into c, c gets P(a) + P(b)
     *   | 'a' is chosen such that P(chose a) = P(a)/(P(a)+P(b))
     * <p>
     * Split: I have a high weight, my bin has < optNumPerBin
     *   | If I split, I need to find someone else to pick up my other half
     *   | If I have a weight > 2P/n, I should split such that each gets 1/m * P(parent) where P/n < m < 2P/n
     * <p>
     * Nothing: My bin has < optNumPerBin
     *   | If I do nothing, it's easy
     *   | Default to this if I can't find something to do or someone to split with
     *  <p>
     */
    private void resample(){
        // Initialize
        ArrayList<Integer>[] binRank = new ArrayList[numBins];
        PriorityQueue<Decision> merges = new PriorityQueue<>(); // Head is the lowest weight for priority queues
        PriorityQueue<Decision> splits = new PriorityQueue<>();
        for (int i = 0; i < numBins; i++){
            binRank[i] = new ArrayList<>();
        }

        // Sort ranks into bins
        for (int i = 0; i < worldSize; i++){
            int bin = (int) Math.round(weightsBins[i][1]);
            binRank[bin].add(i);
        }

        // Analyze each bin
        int m = 2;  // Number of particles to split into (default 2 -- must be > 1)
        for (int i = 0; i < numBins; i++){
            List<Integer> ranks = binRank[i]
                    .stream()
                    .sorted((a, b) -> Double.compare(weightsBins[a][0], weightsBins[b][0]))
                    .toList();
            if (ranks.isEmpty()){ continue; }
            double idealWeight = ranks.stream().mapToDouble(a -> weightsBins[a][0]).sum() / ranks.size();

            // Split Decisions
            for (int rank : ranks) {
                double weight = weightsBins[rank][0];
                if (weightsBins[rank][0] > 2.0 * idealWeight / m || enteredNewBins[rank]) {
                    splits.add(new Decision(new ArrayList<>(List.of(rank)), new ArrayList<>(List.of(weight))));
                }
            }

            // Merge Decisions
            Iterator<Integer> ranksIterator = ranks.stream().iterator();
            int rank = ranksIterator.next();
            double weight = weightsBins[rank][0]; // Smallest weight in bin
            double combinedWeight = 0.0;
            ArrayList<Decision> binMergeList = new ArrayList<>();
            binMergeList.add(new Decision(new ArrayList<>(), new ArrayList<>()));
            while(weight < idealWeight/2){
                if (!enteredNewBins[rank]) { // Don't merge if the rank just entered the bin?
                    if (combinedWeight + weight < idealWeight) { // Noraml case
                        binMergeList.getLast().ranks.add(rank);
                        binMergeList.getLast().weights.add(weight);
                        combinedWeight += weight;
                    } else if (combinedWeight + weight <= idealWeight * 1.5) { // Still allowed to merge with next rank
                        binMergeList.getLast().ranks.add(rank);
                        binMergeList.getLast().weights.add(weight);
                        combinedWeight = 0;
                        binMergeList.add(new Decision(new ArrayList<>(), new ArrayList<>()));
                    } else { // Can't merge with next rank, but still need to consider the current rank so skip the iteration
                        binMergeList.add(new Decision(new ArrayList<>(), new ArrayList<>()));
                        combinedWeight = 0;
                        continue;
                    }
                }

                if (ranksIterator.hasNext()){
                    rank = ranksIterator.next();
                    weight = weightsBins[rank][0];
                } else {
                    break;
                }
            }
            for (Decision decision : binMergeList){
                if (decision.ranks.size() > 1){
                    merges.add(decision);
                }
            }
        }

        // Apply decisions to ranks
        for (int i = 0; i < worldSize; i++){
            weightsBinsCopyRank[i][0] = weightsBins[i][0];
            weightsBinsCopyRank[i][1] = weightsBins[i][1];
            weightsBinsCopyRank[i][2] = -1;
        }
        // Merges
        ArrayList<Integer> freedRanks = new ArrayList<>();
        while (!merges.isEmpty()){
            Decision decision = merges.poll();
            // Sample a rank to merge into from all ranks in the decision to find target
            ArrayList<Double> weights = decision.weights;
            weights.replaceAll(a -> a / weights.stream().mapToDouble(Double::doubleValue).sum());
            double rand = Math.random();
            double cumulativeWeight = 0.0;
            int rankToMergeInto = -1;
            for (int i = 0; i < weights.size(); i++){
                cumulativeWeight += weights.get(i);
                if (rand <= cumulativeWeight){
                    rankToMergeInto = decision.ranks.get(i);
                    break;
                }
            }
            // Add non-target ranks to the freed ranks list & update the target rank
            for (int rank : decision.ranks){
                if (rank != rankToMergeInto){
                    freedRanks.add(rank);
                } else{
                    weightsBinsCopyRank[rank][0] = decision.getTotalWeight();
                }
            }
        }
        // Splits
        while(!splits.isEmpty() && !freedRanks.isEmpty()){
            if (freedRanks.size() < m-1){ // M will change from 2 eventually maybe
                m = freedRanks.size() + 1;
            }
            Decision decision = splits.poll();
            int parent = decision.ranks.getFirst();
            double weightToEach = decision.getTotalWeight() / m;
            for (int i = 0; i < m-1; i++){
                int childRank = freedRanks.removeFirst();
                weightsBinsCopyRank[parent][0] = weightToEach; // Update parent weight
                weightsBinsCopyRank[childRank][0] = weightToEach; // Update child weight
                weightsBinsCopyRank[childRank][1] = weightsBinsCopyRank[parent][1]; // Update child bin
                weightsBinsCopyRank[childRank][2] = parent; // Update where child needs to copy from
            }
        }

        // Apply my own decisions
        myWeightBin[0] = weightsBinsCopyRank[rank][0];
        myWeightBin[1] = weightsBinsCopyRank[rank][1];
    }

    private record Decision(ArrayList<Integer> ranks, ArrayList<Double> weights)
    implements Comparable<Decision> {
        //TODO: Add m as a parameter to separate decisions
        public double getTotalWeight(){
            return weights.stream().mapToDouble(Double::doubleValue).sum();
        }

        @Override
        public int compareTo(Decision decision) {
            return Double.compare(this.getTotalWeight(), decision.getTotalWeight());
        }

    }

    private void binAssignment() {
        double[] global = new double[worldSize];
        for (int i = 0; i < worldSize; i++){
            global[i] = globalValues[i][0];
        }
        double[] binBounds = getOneDimBinBounds(global);
        // TODO: Implement another communication to get global values to make bins dynamically
        // binBounds.length+1 is number of binds because we need to account for the two extra bins at the ends
        // 0 is the bin for < binBounds[0]
        // binBounds.length is the bin for >= binBounds[binBounds.length-1]
        for (int j = 0; j < worldSize; j++) {
            int oldBin = (int) Math.round(weightsBins[j][1]);
            for (int i = 0; i < numBins - 1; i++) {
                if (i == 0 && global[j] < binBounds[i]) {
                    weightsBins[j][1] = i;
                    break;
                }
                if (global[j] >= binBounds[i] && global[j] < binBounds[i + 1]) {
                    weightsBins[j][1] = i;
                    break;
                }
                if (i == numBins - 2 && global[j] >= binBounds[i + 1]) {
                    weightsBins[j][1] = i;
                    break;
                }
            }
            enteredNewBins[j] = oldBin != weightsBins[j][1];
        }
    }

    private void calculateMyMetric(){
        switch (metric){
            case RMSD:
                myGlobalValue[0] = Superpose.rmsd(refCoords, x,  potential.getMass());
                break;
            case RESIDUE_DISTANCE:
                break;
            case COM_DISTANCE:
                break;
            case ATOM_DISTANCE:
                break;
            case POTENTIAL_MIN:
                break;
            case POTENTIAL_MAX:
                break;
            default:
                break;
        }
    }

    private double[] getOneDimBinBounds(double[] globalValues){
        if (staticBins){
            return binBounds;
        }
        //TODO: Implement Voronoi binning for dynamic binning
        logger.severe(" Automatic binning not implemented yet.");
        switch (metric){
            case RMSD:
                logger.severe(" RMSD binning not implemented yet.");
                break;
            case RESIDUE_DISTANCE:
                logger.severe(" Residue distance binning not implemented yet.");
                break;
            case COM_DISTANCE:
                logger.severe(" COM distance binning not implemented yet.");
                break;
            case ATOM_DISTANCE:
                logger.severe(" Atom distance binning not implemented yet.");
                break;
            case POTENTIAL_MIN:
                logger.severe(" Potential min binning not implemented yet.");
                break;
            case POTENTIAL_MAX:
                logger.severe(" Potential max binning not implemented yet.");
                break;
            default:
                logger.severe(" Invalid metric for binning.");
                break;
        }
        return new double[0];
    }

    /**
     * Communicate weights and bin numbers between ranks
     */
    private void comms(){
        try{
            world.allGather(myWeightBinBuf, weightsBinsBuf);
            world.allGather(myGlobalValueBuf, globalValuesBuf);
        } catch (IOException e) {
            String message = " WeightedEnsemble allGather failed.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    private void commsSanityCheckAndFileMigration() {
        comms();
        // Check if the weightsBinsCopyRank & weightsBins are the same
        for (int i = 0; i < worldSize; i++) {
            if (weightsBins[i][0] != weightsBinsCopyRank[i][0] ||
                    weightsBins[i][1] != weightsBinsCopyRank[i][1]) {
                String message = " Rank " + i + " has mismatched weightsBins and weightsBinsCopyRank.";
                logger.severe(message);
                System.exit(1);
            }
        }
        // Get the files I need to copy into temp files then override my files
        if (weightsBinsCopyRank[rank][2] != -1 && weightsBinsCopyRank[rank][2] != rank){ // Second not strictly necessary
            copyOver((int) weightsBinsCopyRank[rank][2]);
            comms(); // Just so everyone is on the same page
            moveOnto();
        }
    }

    private void copyOver(int rank){
        // TODO: Figure out how to do this efficiently even as files get exceedingly large
        // Some graph algorithm w/ local storage of info in text file?
        File dyn = new File(currentDir.getParent() + File.separator + rank +
                File.separator + dynFile.getName());
        dynTemp = new File(currentDir.getParent() + File.separator + rank +
                File.separator + dynFile.getName() + ".temp");
        File traj = new File(currentDir.getParent() + File.separator + rank +
                File.separator + trajFile.getName());
        trajTemp = new File(currentDir.getParent() + File.separator + rank +
                File.separator + trajFile.getName() + ".temp");

        try {
            Files.copy(dyn.toPath(), dynTemp.toPath());
        } catch (IOException e) {
            String message = " Failed to copy dyn file from rank " + rank + " to rank " + this.rank;
            logger.log(Level.SEVERE, message, e);
        }
        try {
            Files.copy(traj.toPath(), trajTemp.toPath());
        } catch (IOException e) {
            String message = " Failed to copy traj file from rank " + rank + " to rank " + this.rank;
            logger.log(Level.SEVERE, message, e);
        }
    }

    private void moveOnto(){
        try {
            Files.move(dynTemp.toPath(), dynFile.toPath());
        } catch (IOException e) {
            String message = " Failed to move dyn file from rank " + rank + " to rank " + this.rank;
            logger.log(Level.SEVERE, message, e);
        }
        try {
            Files.move(trajTemp.toPath(), trajFile.toPath());
        } catch (IOException e) {
            String message = " Failed to move traj file from rank " + rank + " to rank " + this.rank;
            logger.log(Level.SEVERE, message, e);
        }
    }
}
