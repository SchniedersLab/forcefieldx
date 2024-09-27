package ffx.algorithms.optimize;

import ffx.algorithms.optimize.manybody.ManyBodyCell;
import ffx.crystal.Crystal;
import ffx.crystal.SymOp;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueState;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.Double.parseDouble;
import static java.lang.String.format;

class BoxOptimization {

    private final RotamerOptimization rotamerOptimization;
    /**
     * Logger for this class.
     */
    private static final Logger logger = Logger.getLogger(RotamerOptimization.class.getName());
    /**
     * Number of boxes for box optimization in X, Y, Z.
     */

    public final int[] numXYZCells = {3, 3, 3};
    /**
     * Box border size.
     */
    public double cellBorderSize = 0;
    /**
     * Approximate box size.
     */
    public double approxBoxLength = 0;
    /**
     * Box optimization inclusion criteria.
     */
    public int boxInclusionCriterion = 1;
    /**
     * Index of the first box to optimize.
     */
    public int cellStart = 0;
    /**
     * Index of the last box to optimize.
     */
    public int cellEnd = -1;
    /**
     * Flag to indicate manual definition of a super box.
     */
    public boolean manualSuperbox = false;
    /**
     * Dimensions of the box.
     */
    public double[] boxDimensions;
    /**
     * Buffer size for the super box.
     */
    public double superboxBuffer = 8.0;
    /**
     * Box index loaded during a restart.
     */
    public int boxLoadIndex = -1;
    /**
     * Box indeces loaded during a restart.
     */
    public int[] boxLoadCellIndices;

    public BoxOptimization(RotamerOptimization rotamerOptimization) {
        this.rotamerOptimization = rotamerOptimization;
    }

    /**
     * Breaks down a structure into a number of overlapping boxes for optimization.
     *
     * @return Potential energy of final structure.
     */
    public double boxOptimization(List<Residue> residueList) throws Exception {
        rotamerOptimization.usingBoxOptimization = true;
        long beginTime = -System.nanoTime();
        Residue[] residues = residueList.toArray(new Residue[0]);
        /*
         * A new dummy Crystal will be constructed for an aperiodic system. The
         * purpose is to avoid using the overly large dummy Crystal used for
         * Ewald purposes. Atoms are not and should not be moved into the dummy
         * Crystal boundaries; to check if an Atom is inside a cell, use an
         * array of coordinates adjusted to be 0 < coordinate < 1.0.
         */
        Crystal crystal = generateSuperbox(residueList);

        // Cells indexed by x*(YZ divisions) + y*(Z divisions) + z.
        int totalCells = getTotalCellCount(crystal); // Also initializes cell count if using -bB
        if (cellStart > totalCells - 1) {
            rotamerOptimization.logIfRank0(format(" First cell out of range (%d) -- reset to first cell.", cellStart + 1));
            cellStart = 0;
        }
        if (cellEnd > totalCells - 1) {
            // Warn the user if the box end was explicitly set incorrectly.
            if (cellEnd != -1 && cellEnd != Integer.MAX_VALUE) {
                rotamerOptimization.logIfRank0(format(" Final cell out of range (%d) -- reset to last cell.", cellEnd + 1));
            }
            cellEnd = totalCells - 1;
        } else if (cellEnd < 0) {
            cellEnd = totalCells - 1;
        } else if(!crystal.aperiodic()){
            cellEnd = totalCells - 1;
        }
        ManyBodyCell[] cells = loadCells(crystal, residues);
        int numCells = cells.length;
        rotamerOptimization.logIfRank0(format(" Optimizing cells %d to %d", (cellStart + 1), (cellEnd + 1)));
        for (int i = 0; i < numCells; i++) {
            ManyBodyCell manyBodyCell = cells[i];
            List<Residue> residueSubsetList = manyBodyCell.getResiduesAsList();
            int[] cellIndices = manyBodyCell.getABCIndices();
            rotamerOptimization.logIfRank0(format("\n Iteration %d of cell based optimization.", (i + 1)));
            rotamerOptimization.logIfRank0(manyBodyCell.toString());
            int nResidueSubset = residueSubsetList.size();
            if (nResidueSubset > 0) {
                if (rotamerOptimization.rank0 && rotamerOptimization.writeEnergyRestart && rotamerOptimization.printFiles) {
                    String boxHeader = format(" Box %d: %d,%d,%d", i + 1, cellIndices[0], cellIndices[1], cellIndices[2]);
                    try {
                        rotamerOptimization.energyWriter.append(boxHeader);
                        rotamerOptimization.energyWriter.newLine();
                    } catch (IOException ex) {
                        logger.log(Level.SEVERE, " Exception writing box header to energy restart file.", ex);
                    }
                }
                if (rotamerOptimization.loadEnergyRestart) {
                    boxLoadIndex = i + 1;
                    boxLoadCellIndices = new int[3];
                    boxLoadCellIndices[0] = cellIndices[0];
                    boxLoadCellIndices[1] = cellIndices[1];
                    boxLoadCellIndices[2] = cellIndices[2];
                }

                long boxTime = -System.nanoTime();
                Residue firstResidue = residueSubsetList.getFirst();
                Residue lastResidue = residueSubsetList.get(nResidueSubset - 1);
                Residue[] residueSubsetArray = new Residue[residueSubsetList.size()];
                residueSubsetList.toArray(residueSubsetArray);
                if (rotamerOptimization.revert) {
                    ResidueState[] coordinates = ResidueState.storeAllCoordinates(residueSubsetList);
                    double startingEnergy = 0;
                    double finalEnergy = 0;
                    try {
                        startingEnergy = rotamerOptimization.currentEnergy(residueSubsetArray);
                    } catch (ArithmeticException ex) {
                        logger.severe(format(" Exception %s in calculating starting energy of a box; FFX shutting down", ex));
                    }
                    rotamerOptimization.globalOptimization(residueSubsetList);
                    try {
                        finalEnergy = rotamerOptimization.currentEnergy(residueSubsetArray);
                    } catch (ArithmeticException ex) {
                        logger.severe(format(" Exception %s in calculating starting energy of a box; FFX shutting down", ex));
                    }
                    if (startingEnergy <= finalEnergy) {
                        logger.info(
                                "Optimization did not yield a better energy. Reverting to original coordinates.");
                        ResidueState.revertAllCoordinates(residueSubsetList, coordinates);
                    } else {
                        // Copy sliding window optimal rotamers into the overall optimum array.
                        int r = 0;
                        for (Residue residue : residueSubsetList) {
                            int index = residueList.indexOf(residue);
                            rotamerOptimization.optimum[index] = rotamerOptimization.optimumSubset[r++];
                        }
                    }
                    long currentTime = System.nanoTime();
                    boxTime += currentTime;
                    if(rotamerOptimization.genZ){
                        int[] currentRotamers = new int[rotamerOptimization.optimumSubset.length];
                        rotamerOptimization.getFractions(residueSubsetArray,0,currentRotamers, true);
                        rotamerOptimization.getProtonationPopulations(residueSubsetArray);
                    }
                    rotamerOptimization.logIfRank0(format(" Time elapsed for this iteration: %11.3f sec", boxTime * 1.0E-9));
                    rotamerOptimization.logIfRank0(format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
                } else {
                    rotamerOptimization.globalOptimization(residueSubsetList);
                    // Copy sliding window optimal rotamers into the overall optimum array.
                    int r = 0;
                    for (Residue residue : residueSubsetList) {
                        int index = residueList.indexOf(residue);
                        rotamerOptimization.optimum[index] = rotamerOptimization.optimumSubset[r++];
                    }
                    long currentTime = System.nanoTime();
                    boxTime += currentTime;
                    if(rotamerOptimization.genZ){
                        int[] currentRotamers = new int[rotamerOptimization.optimumSubset.length];
                        rotamerOptimization.getFractions(residueSubsetArray,0,currentRotamers, true);
                        rotamerOptimization.getProtonationPopulations(residueSubsetArray);
                    }
                    rotamerOptimization.logIfRank0(format(" Time elapsed for this iteration: %11.3f sec", boxTime * 1.0E-9));
                    rotamerOptimization.logIfRank0(format(" Overall time elapsed: %11.3f sec", (currentTime + beginTime) * 1.0E-9));
                }
                if (rotamerOptimization.rank0 && rotamerOptimization.printFiles) {
                    // Don't write a file if it's the final iteration.
                    if (i == (numCells - 1)) {
                        continue;
                    }
                    try {
                        if (firstResidue != lastResidue) {
                            rotamerOptimization.logIfRank0(format(" File with residues %s ... %s in window written.", firstResidue, lastResidue));
                        } else {
                            rotamerOptimization.logIfRank0(format(" File with residue %s in window written.", firstResidue));
                        }
                    } catch (Exception e) {
                        logger.warning("Exception writing to file.");
                    }
                }
            } else {
                rotamerOptimization.logIfRank0(" Empty box: no residues found.");
            }
        }
        return 0.0;
    }

    public void update(String boxDim) {
        // String should be in format (buffer,xmin,xmax,ymin,ymax,zmin,zmax)
        try {
            String[] bdTokens = boxDim.split(",+");
            boxDimensions = new double[6];
            if (bdTokens.length != 7) {
                logger.warning(" Improper number of arguments to boxDimensions; default settings used.");
            } else {
                for (int i = 1; i < 7; i += 2) {
                    boxDimensions[i - 1] = parseDouble(bdTokens[i]);
                    boxDimensions[i] = parseDouble(bdTokens[i + 1]);
                    if (boxDimensions[i] < boxDimensions[i - 1]) {
                        logger.info(format(" Improper dimension min %8.5f > max %8.5f; max/min reversed.",
                                boxDimensions[i - 1], boxDimensions[i]));
                        double temp = boxDimensions[i];
                        boxDimensions[i] = boxDimensions[i - 1];
                        boxDimensions[i - 1] = temp;
                    }
                }
                superboxBuffer = parseDouble(bdTokens[0]);
                manualSuperbox = true;
            }
        } catch (Exception ex) {
            logger.warning(
                    format(" Error in parsing box dimensions: input discarded and defaults used: %s.", ex));
            manualSuperbox = false;
        }
    }

    private Crystal generateSuperbox(List<Residue> residueList) {
        double[] maxXYZ = new double[3];
        double[] minXYZ = new double[3];
        Crystal originalCrystal = rotamerOptimization.molecularAssembly.getCrystal();
        if (manualSuperbox) {
            for (int i = 0; i < maxXYZ.length; i++) {
                int ii = 2 * i;
                minXYZ[i] = boxDimensions[ii] - superboxBuffer;
                maxXYZ[i] = boxDimensions[ii + 1] + superboxBuffer;
            }
        } else if (originalCrystal.aperiodic()) {
            if (residueList == null || residueList.isEmpty()) {
                throw new IllegalArgumentException(
                        " Null or empty residue list when generating superbox.");
            }
            Atom initializerAtom = residueList.get(0).getReferenceAtom();
            initializerAtom.getXYZ(minXYZ);
            initializerAtom.getXYZ(maxXYZ);
            for (Residue residue : residueList) {
                Atom refAtom = residue.getReferenceAtom();
                double[] refAtomCoords = new double[3];
                refAtom.getXYZ(refAtomCoords);
                for (int i = 0; i < 3; i++) {
                    maxXYZ[i] = Math.max(refAtomCoords[i], maxXYZ[i]);
                    minXYZ[i] = Math.min(refAtomCoords[i], minXYZ[i]);
                }
            }
            for (int i = 0; i < 3; i++) {
                minXYZ[i] -= superboxBuffer;
                maxXYZ[i] += superboxBuffer;
            }
        } else {
            return originalCrystal.getUnitCell();
        }
        double newA = maxXYZ[0] - minXYZ[0];
        double newB = maxXYZ[1] - minXYZ[1];
        double newC = maxXYZ[2] - minXYZ[2];
        if (manualSuperbox) {
            logger.info(format(" Manual superbox set over (minX, maxX, minY, "
                            + "maxY, minZ, maxZ): %f, %f, %f, %f, %f, %f", minXYZ[0], maxXYZ[0], minXYZ[1],
                    maxXYZ[1], minXYZ[2], maxXYZ[2]));
        } else { // Crystal systems will have already returned.
            logger.info(
                    " System is aperiodic: protein box generated over these coordinates (minX, maxX, minY, maxY, minZ, maxZ):");
            String message = " Aperiodic box dimensions: ";
            for (int i = 0; i < minXYZ.length; i++) {
                message = message.concat(format("%f,%f,", minXYZ[i], maxXYZ[i]));
            }
            message = message.substring(0, message.length() - 1);
            logger.info(message);
        }
        logger.info(format(" Buffer size (included in dimensions): %f\n", superboxBuffer));
        return new Crystal(newA, newB, newC, 90.0, 90.0, 90.0, "P1");
    }

    /**
     * Returns the number of cells (boxes) for box optimization; if the -bB flag is set, sets the
     * final number of cells.
     *
     * @param crystal Crystal or dummy crystal being used to define boxes.
     * @return Total number of cells.
     */
    private int getTotalCellCount(Crystal crystal) {
        int numCells = 1;
        if (approxBoxLength > 0) {
            double[] boxes = new double[3];
            boxes[0] = crystal.a / approxBoxLength;
            boxes[1] = crystal.b / approxBoxLength;
            boxes[2] = crystal.c / approxBoxLength;
            for (int i = 0; i < boxes.length; i++) {
                if (boxes[i] < 1) {
                    numXYZCells[i] = 1;
                } else {
                    numXYZCells[i] = (int) boxes[i];
                }
            }
        }
        for (int numXYZBox : numXYZCells) {
            numCells *= numXYZBox;
        }
        return numCells;
    }

    /**
     * Creates and fills cells (boxes) for box optimization.
     *
     * @param crystal  Polymer crystal or dummy crystal
     * @param residues All residues to be optimized
     * @return Filled cells.
     */
    @SuppressWarnings("fallthrough")
    private ManyBodyCell[] loadCells(Crystal crystal, Residue[] residues) {
        double aCellBorderFracSize = (cellBorderSize / crystal.a);
        double bCellBorderFracSize = (cellBorderSize / crystal.b);
        double cCellBorderFracSize = (cellBorderSize / crystal.c);
        int numCells = cellEnd - cellStart + 1;
        rotamerOptimization.logIfRank0(format(" Number of fractional cells: %d = %d x %d x %d",
                numCells, numXYZCells[0], numXYZCells[1], numXYZCells[2]));

        ManyBodyCell[] cells = new ManyBodyCell[numCells];
        int currentIndex = 0;
        int filledCells = 0;
        int[] xyzIndices = new int[3];
        boolean doBreak = false; // Breaks the ijk loop if the last box passed.
        /*
         * Initializes coordinates for all the boxes, indexed linearly along z,
         * then y, then x (so the box with xyz indices 2,3,2 in a crystal with
         * 4, 5, and 3 boxes along xyz would be indexed 2*5*3 + 3*3 + 2 = 41).
         * The int[] indices stores separate x, y, and z indices.
         */
        for (int i = 0; i < numXYZCells[0]; i++) {
            if (doBreak) {
                break;
            }
            xyzIndices[0] = i;
            for (int j = 0; j < numXYZCells[1]; j++) {
                if (doBreak) {
                    break;
                }
                xyzIndices[1] = j;
                for (int k = 0; k < numXYZCells[2]; k++) {
                    if (currentIndex < cellStart) {
                        ++currentIndex;
                        continue;
                    } else if (currentIndex > cellEnd) {
                        doBreak = true;
                        break;
                    }
                    xyzIndices[2] = k;
                    double[] fracCoords = new double[6];
                    fracCoords[0] = (((1.0 * i) / numXYZCells[0]) - aCellBorderFracSize);
                    fracCoords[1] = (((1.0 * j) / numXYZCells[1]) - bCellBorderFracSize);
                    fracCoords[2] = (((1.0 * k) / numXYZCells[2]) - cCellBorderFracSize);
                    fracCoords[3] = (((1.0 + i) / numXYZCells[0]) + aCellBorderFracSize);
                    fracCoords[4] = (((1.0 + j) / numXYZCells[1]) + bCellBorderFracSize);
                    fracCoords[5] = (((1.0 + k) / numXYZCells[2]) + cCellBorderFracSize);
                    cells[filledCells++] = new ManyBodyCell(fracCoords, xyzIndices, currentIndex);
                    ++currentIndex;

                }
            }
        }
        assignResiduesToCells(crystal, residues, cells);
        for (ManyBodyCell cell : cells) {
            cell.sortCellResidues();
        }
        switch (rotamerOptimization.direction) {
            case BACKWARD:
                ManyBodyCell[] tempCells = new ManyBodyCell[numCells];
                for (int i = 0; i < numCells; i++) {
                    tempCells[i] = cells[numCells - (i + 1)];
                }
                cells = tempCells;
                // Fall through into forward case (for now).
            case FORWARD:
            default:
                return cells;
        }
    }

    /**
     * Constructs the cells for box optimization and assigns them residues, presently based on C
     * alpha fractional coordinates; by default, cells are sorted by global index. Presently,
     * specifying approxBoxLength over-rides numXYZBoxes, and always rounds the number of boxes down
     * (to ensure boxes are always at least the specified size).
     *
     * @param crystal  Crystal group.
     * @param residues List of residues to be optimized.
     * @param cells    BoxOptCell instance.
     */
    private void assignResiduesToCells(Crystal crystal, Residue[] residues, ManyBodyCell[] cells) {
        // Search through residues, add them to all boxes containing their
        // fractional coordinates.
        int nSymm = crystal.spaceGroup.getNumberOfSymOps();

        for (ManyBodyCell cell : cells) {
            Set<Residue> toAdd = new HashSet<>();
            for (int iSymm = 0; iSymm < nSymm; iSymm++) {
                SymOp symOp = crystal.spaceGroup.getSymOp(iSymm);
                for (Residue residue : residues) {
                    boolean contained;
                    switch (boxInclusionCriterion) {
                        default -> contained = cell.atomInsideCell(residue.getReferenceAtom(), crystal, symOp);
                        case 2 -> contained = cell.residueInsideCell(residue, crystal, symOp, true);
                        case 3 -> contained = cell.anyRotamerInsideCell(residue, crystal, symOp, true);
                    }
                    if (contained) {
                        toAdd.add(residue);
                    }
                }
                // If the identity symop produces nothing, skip checking other symops.
                if (toAdd.isEmpty()) {
                    break;
                }
            }
            toAdd.forEach(cell::addResidue);
        }
    }
}
