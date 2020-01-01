//******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2020.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
//******************************************************************************
package ffx.utilities;

import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.ServerSocket;

/**
 * Port Utilities.
 *
 * @author Michael J. Schnieders
 */
public class PortUtils {

    /**
     * The maximum TCP Port (65535).
     */
    public static final int MAX_TCP_PORT = 65535;
    /**
     * The minimum TCP Port (0).
     */
    private static final int MIN_TCP_PORT = 0;


    /**
     * Check if a port is available.
     *
     * @param port The port id.
     * @return True if the port is available.
     */
    public static boolean isTcpPortAvailable(int port) {
        try (ServerSocket serverSocket = new ServerSocket()) {
            // setReuseAddress(false) is required only on OSX, otherwise the code will not work correctly on that platform
            serverSocket.setReuseAddress(false);

            // Try to bind the port.
            serverSocket.bind(new InetSocketAddress(InetAddress.getByName("localhost"), port), 1);
            return true;
        } catch (Exception ex) {
            return false;
        }
    }

    /**
     * Check if an int matches a valid TCP port (i.e. is a 16-bit unsigned integer).
     *
     * @param port A number to check.
     * @return If port &ge; 0 and &le; 65535
     */
    public static boolean tcpPortValid(int port) {
        return port >= MIN_TCP_PORT && port <= MAX_TCP_PORT;
    }

}
