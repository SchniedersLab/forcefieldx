/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2015.
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
package ffx.ui;

import javax.swing.SwingUtilities;

/**
 * This is the 3rd version of SwingWorker (also known as SwingWorker 3), an
 * abstract class that you subclass to perform GUI-related work in a dedicated
 * thread. For instructions on using this class, see:
 * http://java.sun.com/docs/books/tutorial/uiswing/misc/threads.html Note that
 * the API changed slightly in the 3rd version: You must now invoke start() on
 * the SwingWorker after creating it.
 *
 * @author Michael J. Schnieders
 *
 */
public abstract class SwingWorker {

    /**
     * Class to maintain reference to current worker thread under separate
     * synchronization control.
     */
    private static class ThreadVar {

        private Thread thread;

        ThreadVar(Thread t) {
            thread = t;
        }

        synchronized void clear() {
            thread = null;
        }

        synchronized Thread get() {
            return thread;
        }
    }
    // private Thread thread;
    private Object value; // see getValue(), setValue()
    private ThreadVar threadVar;

    /**
     * Start a thread that will call the <code>construct</code> method and then
     * exit.
     */
    public SwingWorker() {
        final Runnable doFinished = new Runnable() {
            public void run() {
                finished();
            }
        };
        Runnable doConstruct = new Runnable() {
            public void run() {
                try {
                    setValue(construct());
                } finally {
                    threadVar.clear();
                }
                SwingUtilities.invokeLater(doFinished);
            }
        };
        Thread t = new Thread(doConstruct);
        threadVar = new ThreadVar(t);
    }

    /**
     * Compute the value to be returned by the <code>get</code> method.
     *
     * @return a {@link java.lang.Object} object.
     */
    public abstract Object construct();

    /**
     * Called on the event dispatching thread (not on the worker thread) after
     * the <code>construct</code> method has returned.
     */
    public void finished() {
    }

    /**
     * Return the value created by the <code>construct</code> method. Returns
     * null if either the constructing thread or the current thread was
     * interrupted before a value was produced.
     *
     * @return the value created by the <code>construct</code> method
     */
    public Object get() {
        while (true) {
            Thread t = threadVar.get();
            if (t == null) {
                return getValue();
            }
            try {
                t.join();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt(); // propagate
                return null;
            }
        }
    }

    /**
     * Get the value produced by the worker thread, or null if it hasn't been
     * constructed yet.
     *
     * @return a {@link java.lang.Object} object.
     */
    protected synchronized Object getValue() {
        return value;
    }

    /**
     * A new method that interrupts the worker thread. Call this method to force
     * the worker to stop what it's doing.
     */
    public void interrupt() {
        Thread t = threadVar.get();
        if (t != null) {
            t.interrupt();
        }
        threadVar.clear();
    }

    /**
     * Set the value produced by worker thread
     */
    private synchronized void setValue(Object x) {
        value = x;
    }

    /**
     * Start the worker thread.
     */
    public void start() {
        Thread t = threadVar.get();
        if (t != null) {
            t.start();
        }
    }
}
