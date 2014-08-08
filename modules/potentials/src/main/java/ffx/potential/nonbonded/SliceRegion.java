/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.nonbonded;

import edu.rit.pj.IntegerForLoop;
import edu.rit.pj.IntegerSchedule;
import edu.rit.pj.ParallelRegion;
import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import static ffx.potential.nonbonded.SpatialDensityRegion.logger;
import java.nio.DoubleBuffer;
import java.util.logging.Level;

/**
 *
 * @author avdic
 */
public class SliceRegion extends ParallelRegion {

    // A list of atoms.
    Atom atoms[];
    int nAtoms;
    int gX, gY, gZ;
    int basisSize;
    int nSymm;
    int threadCount;
    double coordinates[][][];
    SliceLoop sliceLoop[];
    DoubleBuffer gridBuffer;
    GridInitLoop gridInitLoop[];
    Crystal crystal;
    double initValue = 0.0;
    int gridSize;
    double grid[];

    // Constructor
    public SliceRegion(int gX, int gY, int gZ, double grid[],
            int basisSize, int nSymm,
            int threadCount, Crystal crystal,
            Atom atoms[], double coordinates[][][]) {
        this.atoms = atoms;
        this.nAtoms = atoms.length;
        this.gX = gX;
        this.gY = gY;
        this.gZ = gZ;
        gridSize = gX * gY * gZ * 2;
        this.basisSize = basisSize;
        this.nSymm = nSymm;
        this.threadCount = threadCount;
        this.coordinates = coordinates;
        this.grid = grid;
        if (grid != null) {
            gridBuffer = DoubleBuffer.wrap(grid);
        }
        sliceLoop = new SliceLoop[threadCount];
        gridInitLoop = new GridInitLoop[threadCount];
        for (int i = 0; i < threadCount; i++) {
            gridInitLoop[i] = new GridInitLoop();
        }

    }

    public final void setCrystal(Crystal crystal, int gX, int gY, int gZ) {
        this.crystal = crystal.getUnitCell();
        //assert(this.crystal.spaceGroup.getNumberOfSymOps() == nSymm);
        this.gX = gX;
        this.gY = gY;
        this.gZ = gZ;
        gridSize = gX * gY * gZ * 2;
    }

    public void setGridBuffer(DoubleBuffer grid) {
        gridBuffer = grid;
    }

    public void setLoops(SliceLoop sliceLoops[]) {
        this.sliceLoop = sliceLoops;
    }

    @Override
    public void run() throws Exception {
        int threadIndex = getThreadIndex();
        try {
            execute(0, gridSize - 1, gridInitLoop[threadIndex]);
            execute(0, gZ - 1, sliceLoop[threadIndex]);
        } catch (Exception e) {
            String message = " Exception in SliceRegion.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * <p>
     * Setter for the field <code>initValue</code>.</p>
     *
     * @param initValue a double.
     */
    public void setInitValue(double initValue) {
        this.initValue = initValue;
    }

    private class GridInitLoop extends IntegerForLoop {

        private final IntegerSchedule schedule = IntegerSchedule.fixed();
        // Extra padding to avert cache interference.
        long pad0, pad1, pad2, pad3, pad4, pad5, pad6, pad7;
        long pad8, pad9, pada, padb, padc, padd, pade, padf;

        @Override
        public IntegerSchedule schedule() {
            return schedule;
        }

        @Override
        public void run(int lb, int ub) {
            if (gridBuffer != null) {
                //if (grid != null) {
                for (int i = lb; i <= ub; i++) {
                    //grid[i] = initValue;
                    gridBuffer.put(i, initValue);
                }
            }
        }
    }

}
