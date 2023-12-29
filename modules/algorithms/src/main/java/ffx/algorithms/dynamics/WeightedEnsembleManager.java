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
import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsUtils;
import ffx.potential.utils.ProgressiveAlignmentOfCrystals;
import ffx.potential.utils.StructureMetrics;
import ffx.potential.utils.Superpose;
import org.apache.commons.configuration2.CompositeConfiguration;
import org.apache.commons.io.FilenameUtils;

/**
 * This class implements the Weighted Ensemble algorithm. Simulations are ran in parallel with weights
 * assigned evenly to all simulations at the start. After a specified number of stochastic dynamics,
 * the parallel configurations are resampled (merged or split) and their weights are updated accordingly.
 * <p>
 * Resampling is done by assigning each configuration a bin number. The bin number is determined by some
 * predetermined metric (e.g. RMSD, distance between residues, distance between complexes). The configurations
 * are then sorted by bin number and the configurations with the lowest weights are merged with other
 * configurations in the same bin. The configurations with the highest weights are split into two
 * configurations with half the weight of the parent.
 * <p>
 * The algorithm is based on the Huber & Kim original paper and a following theory paper by Zuckerman:
 * <p>
 * Huber, G. A., & Kim, S. (1996). Weighted-ensemble Brownian dynamics simulations for protein association reactions.
 * Biophysical journal, 70, 97-110.
 * <p>
 * Zuckerman, D. M. (2010). The "Weighted Ensemble" path sampling method is statistically exact for a broad class of
 * stochastic processes and binning procedures. The Journal of Chemical Physics, 132.
 *
 */
public class WeightedEnsembleManager {
    private static final Logger logger = Logger.getLogger(WeightedEnsembleManager.class.getName());


    private final int rank, worldSize;
    private final double[][] weightsBins, weightsBinsCopyRank, globalValues; // Local storage of all weights & bin numbers
    //private final double[] myWeightBin, myGlobalValue;
    private final Comm world;
    private final DoubleBuf[] weightsBinsBuf, globalValuesBuf; // World storage of weights & bin numbers
    private final DoubleBuf myWeightBinBuf, myGlobalValueBuf;

    private boolean restart, staticBins;
    private final boolean[] enteredNewBins;
    private int numBins, optNumPerBin;
    private long totalSteps, numStepsPerResample, cycle;
    private double weight, dt, temp, trajInterval;
    private double[] refCoords, binBounds, x;
    private File dynFile, dynTemp, trajFile, trajTemp, currentDir, refStructureFile;
    private final Random random;
    private final MolecularDynamics molecularDynamics;
    private Potential potential;
    private final OneDimMetric metric;



    public enum OneDimMetric {RMSD, RESIDUE_DISTANCE, COM_DISTANCE, ATOM_DISTANCE, POTENTIAL, RADIUS_OF_GYRATION}

    public WeightedEnsembleManager(OneDimMetric metric, int optNumPerBin,
                                   MolecularDynamics md,
                                   File refStructureFile) {
        // Init MPI stuff
        this.world = Comm.world();
        this.rank = world.rank();
        this.worldSize = world.size();
        // Check if world size is valid
        if (worldSize < 2){
            logger.severe(" Weighted Ensemble requires at least 2 ranks.");
            System.exit(1);
        } else if (worldSize < 10){
            logger.warning(" Weighted Ensemble is not recommended for small scale parallel simulations.");
        }

        // Init FFX & file stuff
        this.refStructureFile = refStructureFile;
        refCoords = getRefCoords(refStructureFile);
        if(refCoords == null){
            logger.severe(" Failed to get reference coordinates.");
            System.exit(1);
        }
        this.dynFile = md.getDynFile();
        this.molecularDynamics = md;
        this.potential = null;
        if(molecularDynamics.molecularAssembly != null){
            this.potential = molecularDynamics.molecularAssembly[0].getPotentialEnergy();
            this.x = potential.getCoordinates(new double[potential.getNumberOfVariables()]);
        } else {
            logger.severe(" Molecular Assembly not set for Molecular Dynamics.");
            System.exit(1);
        }

        logger.info("\n\n ----------------------------- Initializing Weighted Ensemble Run -----------------------------");
        if (initFilesOrPrepRestart()){
            if (restart) {
                logger.info(" Restarting from previous Weighted Ensemble run.");
            }
        } else {
            logger.severe(" Failed to initialize Weighted Ensemble run.");
            System.exit(1);
        }

        // Init Weighted Ensemble stuff
        CompositeConfiguration properties = molecularDynamics.molecularAssembly[0].getProperties();
        this.random = new Random();
        random.setSeed(44); // TODO: Seed in a better way so that its the same for each rank but different each run
        this.binBounds = getPropertyList(properties, "WE.BinBounds");
        logger.info(" Bin bounds: " + Arrays.toString(binBounds));
        this.numBins = binBounds.length + 1;
        this.staticBins = binBounds.length > 0;
        if (staticBins){
            logger.info(" Using static binning with " + numBins + " bins.");
        } else {
            logger.info(" Using dynamic binning.");
        }
        this.numBins = staticBins ? numBins : worldSize / optNumPerBin;
        this.optNumPerBin = optNumPerBin;
        this.metric = metric;
        this.weight = 1.0 / worldSize;
        this.cycle = 0;
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
        weightsBins[rank][0] = weight;
        //myWeightBin = weightsBins[rank];
        //myWeightBin[0] = weight;
        //myGlobalValue = globalValues[rank];
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

    private static double[] getRefCoords(File refStructureFile){
        PotentialsUtils utils = new PotentialsUtils();
        MolecularAssembly assembly = utils.open(refStructureFile);
        return assembly.getPotentialEnergy().getCoordinates(new double[0]);
    }

    private boolean initFilesOrPrepRestart() {
        // My directory info
        logger.info("\n Rank " + rank + " structure is based on: " + refStructureFile.getAbsolutePath());
        File parent = refStructureFile.getParentFile();
        logger.info(" Rank " + rank + " is using parent directory: " + parent.getAbsolutePath());
        currentDir = new File(parent + File.separator + rank);
        logger.info(" Rank " + rank + " is using directory: " + currentDir.getAbsolutePath());
        trajFile = new File(currentDir + File.separator +
                FilenameUtils.getBaseName(refStructureFile.getName()) + ".arc");
        logger.info(" Rank " + rank + " is using trajectory file: " + trajFile.getAbsolutePath());
        molecularDynamics.setArchiveFiles(new File[]{trajFile});
        dynFile = new File(currentDir + File.separator +
                FilenameUtils.getBaseName(refStructureFile.getName()) + ".dyn");
        logger.info(" Rank " + rank + " is using dyn restart file: " + dynFile.getAbsolutePath());
        molecularDynamics.setFallbackDynFile(dynFile);

        restart = !currentDir.mkdir() && checkRestartFiles();
        logger.info(" \n");
        return currentDir.exists();
    }

    private boolean checkRestartFiles(){
        // Check if restart files exist
        Path dynPath = dynFile.toPath();
        Path trajPath = trajFile.toPath();
        //TODO: Read in ESV file counts
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
        if(restart) {
            long restartFrom = getRestartTime(dynFile);
            logger.info(" Restarting from " + restartFrom + " femtoseconds.");
            totalSteps -= restartFrom;
            logger.info(" Remaining cycles: " + (int) Math.ceil((double) totalSteps / numStepsPerResample));
        }

        // Initial information logging
        logger.info("\n ----------------------------- Start Weighted Ensemble Run -----------------------------");
        int numCycles = (int) Math.ceil((double) totalSteps / numStepsPerResample);

        for (int i = 0; i < numCycles; i++) {
            cycle = i;
            dynamics(numStepsPerResample);
            calculateMyMetric();
            comms(); // Prior comm call here for dynamic binning
            binAssignment();
            resample();
            comms();
            sanityCheckAndFileMigration(); // 1 comm call in here
            logger.info("\n ----------------------------- Resampling cycle #" + (i+1) +
                    " complete ----------------------------- ");
        }
        logger.info("\n\n\n ----------------------------- End Weighted Ensemble Run -----------------------------");
    }

    private long getRestartTime(File dynFile) {
        // TODO: Implement this (Time: in dyn file doesn't work)
        //done by reading in traj file (count num models)
        return 0;
    }

    private void dynamics(long numSteps){
        // Convert numStepsPerResample to picoseconds using dt (in femtoseconds)
        double dynamicTotalTime = (numSteps * dt / 1000.0);
        molecularDynamics.setCoordinates(x);
        potential.energy(x, false);
        molecularDynamics.dynamic(numSteps, dt, dynamicTotalTime/3.0,
                dynamicTotalTime/3.0, temp, true, dynFile);
        x = molecularDynamics.getCoordinates();
        molecularDynamics.writeRestart();
    }

    private void calculateMyMetric(){
        switch (metric){
            case RMSD:
                globalValues[rank][0] = Superpose.rmsd(refCoords, x,  potential.getMass());
                break;
            case RESIDUE_DISTANCE:
                break;
            case COM_DISTANCE:
                break;
            case ATOM_DISTANCE:
                break;
            case POTENTIAL:
                double refEnergy = potential.energy(refCoords, false);
                double myEnergy = potential.energy(x, false);
                globalValues[rank][0] = myEnergy - refEnergy;
                break;
            case RADIUS_OF_GYRATION:
                StructureMetrics.radiusOfGyration(x, potential.getMass());
                break;
            default:
                break;
        }
    }

    private void binAssignment() {
        double[] global = new double[worldSize];
        for (int i = 0; i < worldSize; i++){
            global[i] = globalValues[i][0]; // get global values from the comm() call
        }
        double[] binBounds = getOneDimBinBounds(global);
        numBins = binBounds.length + 1;
        // binBounds.length+1 is number of binds because we need to account for the two extra bins at the ends
        // 0 is the bin for < binBounds[0]
        // binBounds.length is the bin for >= binBounds[binBounds.length-1]
        for (int j = 0; j < worldSize; j++) {
            int oldBin = (int) Math.round(weightsBins[j][1]);
            for (int i = 0; i < numBins - 2; i++) { // -2 because there are 2 bins at the ends that get checked at same time
                if (i == 0 && global[j] < binBounds[i]) { // Lower than first bound
                    weightsBins[j][1] = i;
                    break;
                }
                if (global[j] >= binBounds[i] && global[j] < binBounds[i + 1]) { // Between bounds
                    weightsBins[j][1] = i + 1;
                    break;
                }
                if (i == numBins - 3 && global[j] >= binBounds[i + 1]) { // Higher than last bound
                    weightsBins[j][1] = i + 2;
                    break;
                }
            }
            enteredNewBins[j] = oldBin != weightsBins[j][1];
        }

        logger.info("\n Bin bounds:        " + Arrays.toString(binBounds));
        logger.info(" Rank global values with metric \"" + metricToString(metric) + "\": " + Arrays.toString(global));
        logger.info(" Entered new bins: " + Arrays.toString(enteredNewBins));
    }

    private String metricToString(OneDimMetric metric) {
        return switch (metric) {
            case RMSD -> "RMSD";
            case RESIDUE_DISTANCE -> "Residue Distance";
            case COM_DISTANCE -> "COM Distance";
            case ATOM_DISTANCE -> "Atom Distance";
            case POTENTIAL -> "Potential";
            default -> "Invalid Metric";
        };
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
        logger.info("\n\n ----------------------------- Resampling ----------------------------- ");
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
            double idealWeight = ranks.stream().mapToDouble(a -> weightsBins[a][0]).sum() / optNumPerBin;
            logger.info("\n Bin #" + i + " has " + ranks.size() + " ranks with ideal weight " + idealWeight + ".");

            // Split Options
            ArrayList<Integer> splitRanks = new ArrayList<>();
            for (int rank : ranks) {
                double weight = weightsBins[rank][0];
                if (weightsBins[rank][0] > 2.0 * idealWeight || enteredNewBins[rank]) {
                    splits.add(new Decision(new ArrayList<>(List.of(rank)), new ArrayList<>(List.of(weight))));
                    splitRanks.add(rank);
                }
            }

            // Merge Options
            Iterator<Integer> ranksIterator = ranks.stream().iterator();
            int rank = ranksIterator.next();
            double weight = weightsBins[rank][0]; // Smallest weight in bin
            double combinedWeight = 0.0;
            ArrayList<Decision> binMergeList = new ArrayList<>();
            binMergeList.add(new Decision(new ArrayList<>(), new ArrayList<>()));
            while(weight < idealWeight/2){ // Don't merge if the rank just entered the bin?
                boolean split = false;
                if (splitRanks.contains(rank) && ranks.size() > optNumPerBin){ // If the rank wants to split & merge, flip a coin
                    double choice = random.nextDouble();
                    if (choice < .5) {
                        int finalRank = rank;
                        splits.removeIf(decision -> decision.ranks.contains(finalRank));
                    } else {
                        split = true;
                    }
                }
                if (!split){
                    if (combinedWeight + weight < idealWeight) { // Normal case
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
        StringBuilder mergeString = new StringBuilder();
        for (int i = 0; i < worldSize; i++){
            weightsBinsCopyRank[i][0] = weightsBins[i][0];
            weightsBinsCopyRank[i][1] = weightsBins[i][1];
            weightsBinsCopyRank[i][2] = -1;
        }
        // Merges
        ArrayList<Integer> freedRanks = new ArrayList<>();
        int desiredFreeRanks = splits.size()*m - splits.size();
        while (!merges.isEmpty() && desiredFreeRanks > 0){
            Decision decision = merges.poll();
            if(decision.ranks.size()-1 > desiredFreeRanks){
                decision.ranks.subList(desiredFreeRanks+1, decision.ranks.size()).clear();
                decision.weights.subList(desiredFreeRanks+1, decision.weights.size()).clear();
            }
            desiredFreeRanks -= decision.ranks.size()-1;
            // Sample a rank to merge into from all ranks in the decision to find target
            ArrayList<Double> weights = (ArrayList<Double>) decision.weights.clone();
            double totalWeight = weights.stream().mapToDouble(Double::doubleValue).sum();
            weights.replaceAll(a -> a /totalWeight);
            double rand = random.nextDouble();
            double cumulativeWeight = 0.0;
            int rankToMergeInto = -1;
            for (int i = 0; i < weights.size(); i++){
                cumulativeWeight += weights.get(i);
                if (rand <= cumulativeWeight){
                    rankToMergeInto = decision.ranks.get(i);
                    break;
                }
            }
            mergeString.append("\t Ranks ").append(decision.ranks).append(" --> ").append(rankToMergeInto).append("\n");
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
        StringBuilder splitString = new StringBuilder();
        double[] global = new double[worldSize];
        for (int i = 0; i < worldSize; i++){
            global[i] = globalValues[i][0];
        }
        while(!splits.isEmpty() && !freedRanks.isEmpty()){
            if (freedRanks.size() < m-1){ // M will change from 2 eventually maybe
                m = freedRanks.size() + 1;
            }
            Decision decision = splits.poll();
            int parent = decision.ranks.getFirst();
            double weightToEach = decision.getTotalWeight() / m;
            ArrayList<Integer> ranksUsed = new ArrayList<>();
            ranksUsed.add(parent);
            weightsBinsCopyRank[parent][0] = weightToEach; // Update parent weight
            for (int i = 0; i < m-1; i++){
                int childRank = freedRanks.removeFirst();
                ranksUsed.add(childRank);
                weightsBinsCopyRank[childRank][0] = weightToEach; // Update child weight
                weightsBinsCopyRank[childRank][1] = weightsBinsCopyRank[parent][1]; // Update child bin
                weightsBinsCopyRank[childRank][2] = parent; // Update where child needs to copy from
                global[childRank] = globalValues[parent][0];
            }
            splitString.append("\t Rank ").append(parent).append(" --> ").append(ranksUsed).append("\n");
        }

        // Apply my own decisions
        weightsBins[rank][0] = weightsBinsCopyRank[rank][0];
        weightsBins[rank][1] = weightsBinsCopyRank[rank][1];
        globalValues[rank][0] = global[rank];

        // Log decisions
        logger.info("\n ----------------------------- Resampling Decisions ----------------------------- ");
        double[] weights = new double[worldSize];
        int[] bins = new int[worldSize];
        for (int i = 0; i < worldSize; i++) {
            weights[i] = weightsBinsCopyRank[i][0];
            bins[i] = (int) Math.round(weightsBinsCopyRank[i][1]);
        }
        logger.info("\n Rank bin numbers: " + Arrays.toString(bins));
        logger.info(" Rank weights:     " + Arrays.toString(weights));
        logger.info(" Weight sum:      " + Arrays.stream(weights).sum());
        logger.info(" Merges: \n" + mergeString + "\n");
        logger.info(" Splits: \n" + splitString + "\n");
        if ((Arrays.stream(weights).sum()-1) > 1e-6){
            logger.severe(" Weights do not sum to 1.0.");
        }
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

    private double[] getOneDimBinBounds(double[] globalValues){
        if (staticBins){
            return binBounds;
        }
        //TODO: Implement Voronoi binning for dynamic binning
        logger.severe(" Automatic binning not implemented yet.");
        return new double[0];
    }

    private void comms(){
        comms(false);
    } // This is for debugging

    /**
     * Communicate weights and bin numbers between ranks
     */
    private void comms(boolean log){
        if (log) {
            double[] weights = new double[worldSize];
            double[] global = new double[worldSize];
            for (int i = 0; i < worldSize; i++) {
                weights[i] = weightsBinsCopyRank[i][0];
                global[i] = globalValues[i][0];
            }
            logger.info(" Rank bin numbers Pre: " + Arrays.toString(weights));
            logger.info(" Rank global values Pre: " + Arrays.toString(global));
        }
        try{
            world.allGather(myWeightBinBuf, weightsBinsBuf);
            world.allGather(myGlobalValueBuf, globalValuesBuf);
        } catch (IOException e) {
            String message = " WeightedEnsemble allGather for weightsbins failed.";
            logger.severe(message);
        }
        if(log) {
            double[] weights = new double[worldSize];
            double[] global = new double[worldSize];
            for (int i = 0; i < worldSize; i++) {
                weights[i] = weightsBinsCopyRank[i][0];
                global[i] = globalValues[i][0];
            }
            logger.info(" Rank bin numbers Post: " + Arrays.toString(weights));
            logger.info(" Rank global values Post: " + Arrays.toString(global));
        }
    }

    private void sanityCheckAndFileMigration() {
        // Sanity check if the weightsBinsCopyRank & weightsBins are the same
        for (int i = 0; i < worldSize; i++) {
            if (weightsBins[i][0] != weightsBinsCopyRank[i][0] ||
                    weightsBins[i][1] != weightsBinsCopyRank[i][1]) {
                String message = " Rank " + i + " has mismatched weightsBins and weightsBinsCopyRank.";
                logger.info(" WeightsBins: " + Arrays.deepToString(weightsBins));
                logger.info(" WeightsBinsCopyRank: " + Arrays.deepToString(weightsBinsCopyRank));
                logger.severe(message);
                System.exit(1);
            }
        }

        // Get the files I need to copy into temp files then override my files
        if (weightsBinsCopyRank[rank][2] != -1 && weightsBinsCopyRank[rank][2] != rank) { // Second not strictly necessary
            copyOver((int) weightsBinsCopyRank[rank][2]);
        }
        comms(); // Force pause (I had this in the surrounding if statement earlier, and it's a horrible bug!)
        if (weightsBinsCopyRank[rank][2] != -1 && weightsBinsCopyRank[rank][2] != rank) {
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
            Files.move(dynTemp.toPath(), dynFile.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            String message = " Failed to move dyn file from rank " + rank + " to rank " + this.rank;
            logger.log(Level.SEVERE, message, e);
        }
        try {
            Files.move(trajTemp.toPath(), trajFile.toPath(), java.nio.file.StandardCopyOption.REPLACE_EXISTING);
        } catch (IOException e) {
            String message = " Failed to move traj file from rank " + rank + " to rank " + this.rank;
            logger.log(Level.SEVERE, message, e);
        }
    }
}
