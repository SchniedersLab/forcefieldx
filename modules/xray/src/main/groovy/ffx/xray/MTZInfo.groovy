package ffx.xray

import groovy.cli.Option
import groovy.cli.Unparsed
import groovy.cli.picocli.CliBuilder

import ffx.xray.parsers.MTZFilter

/**
 * The MTZInfo script prints out MTZ reflection file header info.
 * <br>
 * Usage:
 * <br>
 * ffxc xray.MTZInfo &lt;filename&gt;
 */
class MTZInfo extends Script {

    /**
     * Options for the MTZInfo Script.
     * <br>
     * Usage:
     * <br>
     * ffxc xray.MTZInfo &lt;filename&gt;
     */
    class Options {
        /**
         * -h or --help to print a help message
         */
        @Option(shortName = 'h', defaultValue = 'false', description = 'Print this help message.')
        boolean help
        /**
         * The final argument should be an MTZ filename.
         */
        @Unparsed(description = 'The MTZ file.')
        List<String> filename
    }

    /**
     * Execute the script.
     */
    def run() {

        def cli = new CliBuilder()
        cli.name = "ffxc xray.MTZInfo"

        def options = new Options()
        cli.parseFromInstance(options, args)

        if (options.help == true) {
            return cli.usage()
        }

        List<String> arguments = options.filename
        String mtzfile
        if (arguments != null && arguments.size() > 0) {
            // Read in command line.
            mtzfile = arguments.get(0)
        } else {
            return cli.usage()
        }

        logger.info("\n Running MTZInfo on " + mtzfile)

        File file = new File(mtzfile)
        if (!file.exists()) {
            println(" File " + mtzfile + " was not found.")
            return
        }

        MTZFilter mtzfilter = new MTZFilter()
        mtzfilter.getReflectionList(file)
        mtzfilter.printHeader()
    }
}

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