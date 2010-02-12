/**
 * Title: Force Field X
 * Description: Force Field X - Software for Molecular Biophysics.
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2010
 *
 * This file is part of Force Field X.
 *
 * Force Field X is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as published
 * by the Free Software Foundation.
 *
 * Force Field X is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Force Field X; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 */
package ffx.ui.commands;

import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.net.InetSocketAddress;
import java.net.Socket;
import java.util.logging.Logger;

/**
 * The FFXClient class encapsulates a socket connection to an TinkerServer
 * started by an executing TINKER program. TinkerSystem and TinkerUpdate objects
 * are sent by the TinkerServer to the FFXClient on request.
 */
public class FFXClient {
	private static final Logger logger = Logger.getLogger(FFXClient.class.getName());
	private Socket client; // Socket connection to the server
	private InetSocketAddress address; // Server address
	private InputStream in; // Input from the server
	private ObjectInputStream oin; // Input form the server
	private OutputStream out; // Output to the server
	private ObjectOutputStream oout; // Output to the server
	private TinkerSystem system; // Tinker system definition
	private TinkerUpdate update; // Tinker update information
	private FFXMessage message; // Message Passed Between Client & Server
	// Various connection status variables
	private int retryCount = 0; // Count of Client attempts to Connect
	private int retryLimit = 10000; // Maximum number of retries
	private boolean connectionMade = false; // True when a connection has been
	// made
	// closed is True if the server closes an open connection
	// or if the retryLimit is reached
	private boolean closed = false;

	public FFXClient() {
		address = new InetSocketAddress(2000);
	}

	public FFXClient(InetSocketAddress a) {
		address = a;
	}

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
			logger.info("Connected to TINKER Server: " + client);
		} catch (Exception e) {
			client = null;
		} finally {
			if (client == null || !client.isConnected()) {
				release();
			}
		}
	}

	public TinkerSystem getSystem() {
		readSocket();
		return system;
	}

	public TinkerUpdate getUpdate() {
		readSocket();
		return update;
	}

	public boolean isClosed() {
		return closed;
	}

	public boolean isConnected() {
		if (client != null && client.isConnected()) {
			return true;
		}
		return false;
	}

	public void readSocket() {
		try {
			while (oin != null && in.available() > 0) {
				Object o = oin.readObject();
				if (o instanceof FFXMessage) {
					message = (FFXMessage) o;
					// logger.info(message.toString());
					if (message.getMessage() == FFXMessage.SYSTEM) {
						system = (TinkerSystem) oin.readObject();
						system.read = false;
					} else if (message.getMessage() == FFXMessage.UPDATE) {
						update = (TinkerUpdate) oin.readObject();
						update.read = false;
					} else if (message.getMessage() == FFXMessage.CLOSING) {
						closed = true;
						release();
					}
				}
			}
			if (system == null) {
				message = new FFXMessage(FFXMessage.SYSTEM);
				oout.reset();
				oout.writeObject(message);
				oout.flush();
			} else if (update == null || update.read) {
				message = new FFXMessage(FFXMessage.UPDATE);
				if (update != null) {
					if (update.type == TinkerUpdate.SIMULATION) {
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
			logger.warning("Exception reading data from TINKER\n"
					+ e.toString());
			release();
		}
	}

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
				FFXMessage close = new FFXMessage(FFXMessage.CLOSING);
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
