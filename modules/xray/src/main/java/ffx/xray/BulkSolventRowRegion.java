package ffx.xray;

import java.util.logging.Level;
import java.util.logging.Logger;

import static java.util.Arrays.fill;

import edu.rit.pj.ParallelTeam;

import ffx.crystal.Crystal;
import ffx.potential.bonded.Atom;
import ffx.potential.nonbonded.RowRegion;


public class BulkSolventRowRegion extends RowRegion {

    /**
     * Constant <code>logger</code>
     */
    protected static final Logger logger = Logger.getLogger(BulkSolventRowRegion.class.getName());

    private final BulkSolventList bulkSolventList;
    private final int gZ;
    private final int gY;

    /**
     * <p>
     * Constructor for BulkSolventDensityRegion.</p>
     *
     * @param gX a int.
     * @param gY a int.
     * @param gZ a int.
     * @param grid an array of double.
     * @param basisSize a int.
     * @param nSymm a int.
     * @param threadCount a int.
     * @param crystal a {@link ffx.crystal.Crystal} object.
     * @param atoms an array of {@link ffx.potential.bonded.Atom} objects.
     * @param coordinates an array of double.
     * @param cutoff a double.
     * @param parallelTeam a {@link edu.rit.pj.ParallelTeam} object.
     */
    public BulkSolventRowRegion(int gX, int gY, int gZ, double grid[],
            int basisSize, int nSymm, int threadCount, Crystal crystal,
            Atom atoms[], double coordinates[][][],
            double cutoff, ParallelTeam parallelTeam) {
        super(gX, gY, gZ, grid, basisSize, nSymm,
                threadCount, crystal, atoms, coordinates);

        this.gZ = gZ;
        this.gY = gY;
        // Asymmetric unit atoms never selected by this class.
        fill(select[0], false);
        bulkSolventList = new BulkSolventList(crystal, atoms, cutoff, parallelTeam);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void run() {
        try {
            execute(0, (gZ*gY) - 1, rowLoop[getThreadIndex()]);
        } catch (Exception e) {
            String message = " Exception in BulkSolventRowRegion.";
            logger.log(Level.SEVERE, message, e);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void selectAtoms() {
        bulkSolventList.buildList(coordinates, select, false);
    }
}
