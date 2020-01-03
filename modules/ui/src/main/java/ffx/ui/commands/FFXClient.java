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
package ffx.ui.commands;

import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.util.logging.Logger;

/**
 * The FFXClient class encapsulates a socket connection to an FFXServer
 * started by an executing FFX instance. FFXSystem and FFXUpdate objects
 * are sent by the FFXServer to the FFXClient on request.
 *
 * @author Michael J. Schnieders
 *
 */
public class FFXClient {

    private static final Logger logger = Logger.getLogger(FFXClient.class.getName());
    private Socket client; // Socket connection to the server
    private InetSocketAddress address; // Server address
    private InputStream in; // Input from the server
    private ObjectInputStream oin; // Input form the server
    private OutputStream out; // Output to the server
    private ObjectOutputStream oout; // Output to the server
    private SimulationDefinition system; // Tinker system definition
    private SimulationUpdate update; // Tinker update information
    private SimulationMessage message; // Message Passed Between Client & Server
    // Various connection status variables
    private int retryCount = 0; // Count of Client attempts to Connect
    private int retryLimit = 10000; // Maximum number of retries
    private boolean connectionMade = false; // True when a connection has been
    // made
    // closed is True if the server closes an open connection
    // or if the retryLimit is reached
    private boolean closed = false;

    /**
     * <p>
     * Constructor for FFXClient.</p>
     */
    public FFXClient() {
        address = new InetSocketAddress(2000);
    }

    /**
     * <p>
     * Constructor for FFXClient.</p>
     *
     * @param a a {@link java.net.InetSocketAddress} object.
     */
    public FFXClient(InetSocketAddress a) {
        address = a;
    }

    /**
     * <p>
     * Constructor for FFXClient.</p>
     *
     * @param port a int.
     */
    public FFXClient(int port) {
        address = new InetSocketAddress(port);
    }

    /**
     * Attempts to connect to a Tinker FServer. If this FClient is already
     * connected, the connection will be closed.
     */
    public void connect() {
        if (client != null && client.isConnected()) {
            release();
        }
        closed = false;
        client = new Socket();
        try {
            client.connect(address, 100);
            client.setTcpNoDelay(true);
            out = client.getOutputStream();
            oout = new ObjectOutputStream(out);
            in = client.getInputStream();
            oin = new ObjectInputStream(in);
            connectionMade = true;
            logger.info("Connected to FFX Server: " + client);
        } catch (Exception e) {
            client = null;
        } finally {
            if (client == null || !client.isConnected()) {
                release();
            }
        }
    }

    /**
     * <p>
     * Getter for the field <code>system</code>.</p>
     *
     * @return a {@link ffx.ui.commands.SimulationDefinition} object.
     */
    public SimulationDefinition getSystem() {
        readSocket();
        return system;
    }

    /**
     * <p>
     * Getter for the field <code>update</code>.</p>
     *
     * @return a {@link ffx.ui.commands.SimulationUpdate} object.
     */
    public SimulationUpdate getUpdate() {
        readSocket();
        return update;
    }

    /**
     * <p>
     * isClosed</p>
     *
     * @return a boolean.
     */
    public boolean isClosed() {
        return closed;
    }

    /**
     * <p>
     * isConnected</p>
     *
     * @return a boolean.
     */
    public boolean isConnected() {
        if (client != null && client.isConnected()) {
            return true;
        }
        return false;
    }

    /**
     * <p>
     * readSocket</p>
     */
    public void readSocket() {
        try {
            while (oin != null && in.available() > 0) {
                Object o = oin.readObject();
                if (o instanceof SimulationMessage) {
                    message = (SimulationMessage) o;
                    // logger.info(message.toString());
                    if (message.getMessage() == SimulationMessage.SYSTEM) {
                        system = (SimulationDefinition) oin.readObject();
                        system.read = false;
                    } else if (message.getMessage() == SimulationMessage.UPDATE) {
                        update = (SimulationUpdate) oin.readObject();
                        update.read = false;
                    } else if (message.getMessage() == SimulationMessage.CLOSING) {
                        closed = true;
                        release();
                    }
                }
            }
            if (system == null) {
                message = new SimulationMessage(SimulationMessage.SYSTEM);
                oout.reset();
                oout.writeObject(message);
                oout.flush();
            } else if (update == null || update.read) {
                message = new SimulationMessage(SimulationMessage.UPDATE);
                if (update != null) {
                    if (update.type == SimulationUpdate.SIMULATION) {
                        message.setTime(update.time);
                    } else {
                        message.setStep(update.step);
                    }
                }
                oout.reset();
                oout.writeObject(message);
                oout.flush();
            }
        } catch (Exception e) {
            logger.warning("Exception reading data from FFX\n"
                    + e.toString());
            release();
        }
    }

    /**
     * <p>
     * release</p>
     */
    public void release() {
        if (client == null) {
            return;
        }
        retryCount++;
        if (retryCount > retryLimit || connectionMade) {
            closed = true;
        }
        if (client != null && client.isConnected() && oout != null) {
            try {
                SimulationMessage close = new SimulationMessage(SimulationMessage.CLOSING);
                oout.reset();
                oout.writeObject(close);
                oout.flush();
            } catch (Exception e) {
                oout = null;
            }
        }
        try {
            if (oin != null) {
                oin.close();
            }
            if (in != null) {
                in.close();
            }
            if (oout != null) {
                oout.close();
            }
            if (out != null) {
                out.close();
            }
            if (client != null) {
                client.close();
            }
        } catch (Exception e) {
            client = null;
        } finally {
            in = null;
            oin = null;
            out = null;
            oout = null;
            client = null;
        }
    }
}
