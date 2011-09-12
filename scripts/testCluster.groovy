// TEST CLUSTER OPERATIONS

// Parallel Java Imports
import edu.rit.pj.Comm;
import edu.rit.mp.DoubleBuf;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   logger.info("\n Usage: ffxc -Dpj.nn=X testCluster filename");
   return;
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

logger.info("\n Running testCluster");
//systems = open(filename);

Comm comm = Comm.world();
int rank = comm.rank();
int size = comm.size();
logger.info(" My rank is " + rank + " out of " + size);

double[][] rootTempAndEnergy = null;
DoubleBuf[] rootBuf = null;

double[] myTempAndEnergy = new double[2];
DoubleBuf myBuf = DoubleBuf.buffer(myTempAndEnergy);

// Initialize the root copy of all temperatures
if (rank == 0) {
    rootTempAndEnergy = new double[size][2];
    rootBuf = new DoubleBuf[size];
    for (int i=0; i<size; i++) {
        rootBuf[i] = DoubleBuf.buffer(rootTempAndEnergy[i]);
        // Initialize temperature and energy
        rootTempAndEnergy[i][0] = 300.0 + i * 10.0;
        rootTempAndEnergy[i][1] = 0.0;
    }  
}

// Scatter the initial temperatures.
comm.scatter(0, rootBuf, myBuf);

// Load the temperture for this process.
double temperature = myTempAndEnergy[0];
double energy = 0.0;

logger.info(" Temperature for rank " + rank + " is " + temperature);

// Pick a random number for the energy of this process.
energy = Math.random() * 100.0;
myTempAndEnergy[1] = energy;

// Gather the energy from each node;
comm.gather(0, myBuf, rootBuf);

if (rank == 0) {
    // Could do RepEx here.
    for (int i=0; i<size; i++) {
        logger.info(" " + i + 
         " Temperature " + rootTempAndEnergy[i][0] +
         " Energy " + rootTempAndEnergy[i][1]);
    }
}






