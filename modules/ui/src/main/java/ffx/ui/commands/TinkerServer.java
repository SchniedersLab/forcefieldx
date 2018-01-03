/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2018.
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
 *
 * Linking this library statically or dynamically with other modules is making a
 * combined work based on this library. Thus, the terms and conditions of the
 * GNU General Public License cover the whole combination.
 *
 * As a special exception, the copyright holders of this library give you
 * permission to link this library with independent modules to produce an
 * executable, regardless of the license terms of these independent modules, and
 * to copy and distribute the resulting executable under terms of your choice,
 * provided that you also meet, for each linked independent module, the terms
 * and conditions of the license of that module. An independent module is a
 * module which is not derived from or based on this library. If you modify this
 * library, you may extend this exception to your version of the library, but
 * you are not obligated to do so. If you do not wish to do so, delete this
 * exception statement from your version.
 */
package ffx.ui.commands;

import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.InetAddress;
import java.net.InetSocketAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.ListIterator;
import java.util.Vector;
import java.util.logging.Logger;

/**
 * The TinkerServer is launched by TINKER executables to allow Force Field X
 * Clients to connect.
 *
 * @author Michael J. Schnieders
 *
 */
@Deprecated
public class TinkerServer implements Runnable {

    private static final Logger logger = Logger.getLogger(TinkerServer.class.getName());
    private static FFXMessage closing = new FFXMessage(FFXMessage.CLOSING);
    private ServerSocket server;
    private int serverPort = 2000;
    private int serverTimeout = 100;
    private Thread thread;
    private int sleepTime = 1;
    private boolean init = true;
    private int cycle = 0;
    private boolean shutdown = false;
    private boolean request = false;
    private Vector<Socket> clients = new Vector<Socket>();
    private Vector<ObjectOutputStream> outputs = new Vector<ObjectOutputStream>();
    private Vector<ObjectInputStream> inputs = new Vector<ObjectInputStream>();
    private TinkerSystem system = null;
    private TinkerUpdate update = null;

    /**
     * <p>
     * Constructor for TinkerServer.</p>
     *
     * @param s a {@link ffx.ui.commands.TinkerSystem} object.
     */
    public TinkerServer(TinkerSystem s) {
        system = s;
    }

    private void accept() {
        Socket client = null;
        ObjectOutputStream oout = null;
        ObjectInputStream oin = null;
        if (server != null) {
            try {
                client = server.accept();
                if (client != null && client.isConnected()) {
                    client.setTcpNoDelay(true);
                    clients.add(client);
                    oout = new ObjectOutputStream(client.getOutputStream());
                    outputs.add(oout);
                    oin = new ObjectInputStream(client.getInputStream());
                    inputs.add(oin);
                    Logger.getLogger("ffx").info(
                            "Client connected\n" + client.toString());
                }
            } catch (Exception e) {
                if (client != null) {
                    clients.remove(client);
                }
                if (oout != null) {
                    outputs.remove(oout);
                }
                if (oin != null) {
                    inputs.remove(oin);
                }
            }
        }
    }

    private void closeClient(int index) {
        Socket client = clients.get(index);
        ObjectOutputStream oout = outputs.get(index);
        ObjectInputStream oin = inputs.get(index);
        try {
            oout.reset();
            oout.writeObject(closing);
            oout.flush();
        } catch (Exception e) {
            try {
                outputs.remove(index);
                inputs.remove(index);
                clients.remove(index);
                if (oout != null) {
                    oout.close();
                }
                if (oin != null) {
                    oin.close();
                }
                if (client != null) {
                    client.close();
                }
            } catch (Exception ex) {
                return;
            }
        }
    }

    private void closeServer() {
        for (int i = 0; i < clients.size(); i++) {
            lastUpdate(i);
        }
        while (!clients.isEmpty()) {
            closeClient(0);
            try {
                synchronized (this) {
                    wait(10);
                }
            } catch (Exception e) {
                Logger.getLogger("ffx").severe(e.toString());
            }
        }
        try {
            if (server != null) {
                server.close();
                server = null;
            }
        } catch (Exception e) {
            return;
        }
    }

    /**
     * <p>
     * isAlive</p>
     *
     * @return a boolean.
     */
    public boolean isAlive() {
        if (thread == null) {
            return false;
        } else if (thread.isAlive()) {
            return true;
        } else {
            return false;
        }
    }

    private void lastUpdate(int index) {
        try {
            ObjectOutputStream oout = outputs.get(index);
            Socket client = clients.get(index);
            if (client != null && client.isConnected() && oout != null) {
                FFXMessage last = new FFXMessage(FFXMessage.SYSTEM);
                if (system != null) {
                    oout.reset();
                    oout.writeObject(last);
                    oout.writeObject(system);
                    oout.flush();
                }
                if (update != null) {
                    oout.reset();
                    last.setMessage(FFXMessage.UPDATE);
                    oout.writeObject(last);
                    oout.writeObject(update);
                    oout.flush();
                }
                last.setMessage(FFXMessage.CLOSING);
                oout.reset();
                oout.writeObject(last);
                oout.flush();
            }
        } catch (Exception e) {
            logger.severe("" + e);
        }
    }

    /**
     * <p>
     * loadUpdate</p>
     *
     * @param u a {@link ffx.ui.commands.TinkerUpdate} object.
     */
    public void loadUpdate(TinkerUpdate u) {
        update = u;
    }

    /**
     * <p>
     * needUpdate</p>
     *
     * @return a boolean.
     */
    public boolean needUpdate() {
        if (clients.size() == 0) {
            sleepTime = 100;
            return false;
        }
        if (request != true) {
            return false;
        }
        sleepTime = 1;
        return true;
    }

    /**
     * <p>
     * run</p>
     */
    public void run() {
        startServer();
        while (!shutdown) {
            accept();
            send();
            try {
                Thread.sleep(sleepTime);
            } catch (Exception e) {
                // What to do here ??
                logger.severe("Server: thread sleep interrupted\n");
            }
            if (init) {
                cycle++;
                if (cycle >= 10) {
                    init = false;
                    serverTimeout = 1;
                }
            }
        }
        accept();
        send();
        closeServer();
    }

    private void send() {
        if (system == null) {
            return;
        }
        if (clients.size() == 0) {
            return;
        }
        ObjectOutputStream oout;
        ObjectInputStream oin;
        Socket client;
        ListIterator<ObjectOutputStream> lout;
        ListIterator<ObjectInputStream> lin;
        ListIterator<Socket> lclient;
        Vector<Socket> closed = new Vector<Socket>();
        for (lout = outputs.listIterator(), lin = inputs.listIterator(), lclient = clients
                .listIterator(); lout.hasNext();) {
            oout = lout.next();
            oin = lin.next();
            client = lclient.next();
            if (!client.isConnected() || client.isClosed()) {
                closed.add(client);
            } else {
                try {
                    FFXMessage message = null;
                    while (oin != null
                            && client.getInputStream().available() > 0) {
                        Object o = oin.readObject();
                        if (o instanceof FFXMessage) {
                            message = (FFXMessage) o;
                            if (message.getMessage() == FFXMessage.CLOSING) {
                                closed.add(client);
                                message = null;
                                break;
                            }
                        }
                    }
                    if (message != null) {
                        if (message.getMessage() == FFXMessage.SYSTEM) {
                            synchronized (system) {
                                oout.reset();
                                oout.writeObject(message);
                                oout.writeObject(system);
                                oout.flush();
                            }
                        } else if (message.getMessage() == FFXMessage.UPDATE) {
                            request = true;
                            if (update != null && update.isNewer(message)) {
                                synchronized (update) {
                                    // logger.info("Sending Update");
                                    oout.reset();
                                    oout.writeObject(message);
                                    oout.writeObject(update);
                                    oout.flush();
                                }
                            }
                        }
                    }
                } catch (Exception e) {
                    closed.add(client);
                }
            }
        }
        for (Socket s : closed) {
            int index = closed.indexOf(s);
            closeClient(index);
        }
    }

    /**
     * <p>
     * setUpdated</p>
     */
    public void setUpdated() {
        request = false;
    }

    /**
     * <p>
     * start</p>
     */
    public void start() {
        if (thread == null || !thread.isAlive()) {
            thread = new Thread(this);
            thread.setPriority(Thread.MAX_PRIORITY);
            thread.start();
        }
    }

    private void startServer() {
        try {
            server = new ServerSocket();
            server.setSoTimeout(serverTimeout);
            server.setReuseAddress(true);
            server.bind(new InetSocketAddress(InetAddress.getLocalHost(),
                    serverPort));
        } catch (Exception e) {
            try {
                server.bind(new InetSocketAddress(InetAddress.getByName(null),
                        serverPort));
            } catch (Exception ex) {
                Logger.getLogger("ffx").severe(
                        "SERVER -- Could not start\n" + e.toString());
                return;
            }
        }
        Logger.getLogger("ffx").info(
                "TINKER Server Address: " + server.getLocalSocketAddress());
        // JOptionPane.showMessageDialog(null, server.getLocalSocketAddress(),
        // "Server Address", JOptionPane.ERROR_MESSAGE);
    }

    /**
     * <p>
     * stop</p>
     */
    public void stop() {
        shutdown = true;
    }
}
