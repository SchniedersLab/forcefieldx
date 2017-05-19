/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2017.
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
package ffx.potential.extended;

import java.util.IllegalFormatException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static ffx.potential.extended.ExtUtils.arrayToStrings;
import static ffx.potential.extended.ExtUtils.prop;

/**
 * The SB utility:
 *  1) relieves ".append(.format(" fatigue,
 *  2) handles newlines succinctly,
 *  3) cuts memory allocation,
 *  4) and prints through the logger of the calling class.

 * Usage:
 *      import static SB.SB;    // That's all the setup you need!
 *      SB.logfn("format", foo, ...);
 *      SB.print();
 * 
 * Of course, you can create your own instances if you need to inter-weave
 * creation of distinct log messages (or need a threaded context).
 * Stop wasting your life away with Logger and switch to SB today!
 */
public final class SBLogger {
    /**
     * A singleton-like shared instantiation.
     */
    public static final SBLogger SB = DefaultInstance.SB;
    private static class DefaultInstance {
        private DefaultInstance() {}
        private static final SBLogger SB = new SBLogger();
    }
	public SBLogger() {}
	private static final Logger fallback = Logger.getLogger(SBLogger.class.getName());
    /**
     * Allows avoidance of resizing.
     */
    private static final int initialCapacity = prop("sys.sbCapacity", 500000);
    /**
     * Identifies the class of method calls so that original loggers can
     * be used whenever possible.
     */
    private static final CallerID cid = new CallerID();
    private StringBuffer sb = new StringBuffer(initialCapacity);
    /**
     * (log) with (f)ormat
     */
    public synchronized void logf(String msg, Object... args) {
        sb.append(format(msg, args));
    }
    /**
     * (log) with (f)ormat, (n)ewline
     */
    public synchronized void logfn(String msg, Object... args) {
        sb.append(format(msg, args)).append("\n");
    }
    /**
     * (n)ewline, (log) with (f)ormat
     */
    public synchronized void nlogf(String msg, Object... args) {
        sb.append("\n").append(format(msg, args));
    }
    /**
     * (n)ewline, (log) with (f)ormat, (n)ewline
     */
    public synchronized void nlogfn(String msg, Object... args) {
        sb.append("\n").append(format(msg, args)).append("\n");
    }
    /**
     * (log) with (f)ormat, (p)rint
     */
    public synchronized void logfp(String msg, Object... args) {
		sb.append(format(msg, args));
		write(Level.INFO);
    }
    /**
     * (n)ew(l)ine
     */
    public synchronized void nl() {
        sb.append("\n");
    }
    /**
     * Send to calling logger.
     */
    public synchronized void print() {
		write(Level.INFO);
    }
    /**
     * Kick it old school to thwart Levels, Filters, Handlers, LogManagers,
     * ResourceBundles, and other such travesties of the modern age.
     */
    public synchronized void force() {
		Logger.getAnonymousLogger().log(Level.INFO, sb.toString());
    }
    /**
     * @see SBLogger::forceToConsole
     */
    public synchronized void force(String msg, Object... args) {
		sb.append(format(msg, args));
		force();
    }
    /**
     * Send to calling logger as warning.
     */
    public synchronized void warning() {
		write(Level.WARNING);
    }
    /**
     * @see SBLogger::warning
     */
    public synchronized void warning(String msg, Object... args) {
		headern(msg, args);
		write(Level.WARNING);
	}
    /**
     * Send to calling logger as severe.
	 * Ensures termination in the event of caught exceptions via System call.
     */
    public synchronized void crash() {
		try		{ write(Level.SEVERE);	}
		finally { System.exit(1);		}
    }
	private final boolean synchronous		 = true;
	private final boolean suppressTrailingNL = true;
	private final boolean useCallerID		 = true;
	private synchronized void write(Level level) {
		String msg = sb.toString();
		if (!msg.isEmpty()) {
			if (suppressTrailingNL && msg.endsWith("%n"))
				msg = msg.substring(0, msg.length()-1);
			if (useCallerID)
				 cid.getCallingLogger().log(level, msg);
			else fallback.log(level, msg);
		}
		clear();
	}
    /**
     * @see SBLogger::crash
     */
    public synchronized void crash(String msg, Object... args) {
		headern(msg, args);
		crash();
    }
	/**
	 * Insert a formatted string preceding current message.
	 */
	public synchronized void header(String msg, Object... args) {
		sb.insert(0, format(msg, args));
	}
    /**
     * Insert a formatted string preceding current message; newline included.
     */
    public synchronized void headern(String msg, Object... args) {
        sb.insert(0, format(msg, args) + format("\n"));
    }
    /**
     * Send to calling logger or discard message.
     */
    public synchronized void printIf(boolean print) {
        if (print) {
            print();
        } else {
            clear();
        }
    }
	public synchronized boolean isEmpty() {
		return sb.toString().isEmpty();
	}
    /**
     * Discard message.
     */
    public synchronized void clear() {
		sb = new StringBuffer();
    }
    /**
     * Attempt to automatically coerce toString() requirements.
     */
    private  synchronized static String format(String msg, Object... args) {
        try {
            return String.format(msg, args);
        } catch (IllegalFormatException ex) {
            // NR: search for indexes of %s rather than replacing all
            return String.format(msg, arrayToStrings(args));
        }
    }
    
    private static class CallerID extends SecurityManager {
        public synchronized Class<?> getCallingClass() {
            Class<?>[] callStack = getClassContext();
            for (Class<?> caller : callStack) {
                if (caller != this.getClass()) {
                    return caller;
                }
            }
            return getClass();      // fallback
        }
        public synchronized Logger getCallingLogger() {
            return Logger.getLogger(getCallingClass().getName());
        }
    }
}
