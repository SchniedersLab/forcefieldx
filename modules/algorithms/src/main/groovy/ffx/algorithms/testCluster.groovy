/**
 * Title: Force Field X Description: Force Field X - Software for Molecular
 * Biophysics Copyright: Copyright (c) Michael J. Schnieders 2001-2012
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 3 as published by
 * the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */

// TEST CLUSTER OPERATIONS

// Parallel Java Imports
import edu.rit.pj.Comm;
import edu.rit.mp.DoubleBuf;

// Things below this line normally do not need to be changed.
// ===============================================================================================

logger.info("\n Running testCluster");

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

