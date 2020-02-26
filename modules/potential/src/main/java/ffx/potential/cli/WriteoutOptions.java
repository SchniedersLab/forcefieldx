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
package ffx.potential.cli;

import ffx.potential.MolecularAssembly;
import ffx.potential.utils.PotentialsFunctions;
import picocli.CommandLine;

import java.io.File;

/**
 * Represents command line options for scripts that periodically write out structures.
 *
 * @author Michael J. Schnieders
 * @author Soham Ali
 * @since 1.0
 */
public class WriteoutOptions {

    /**
     * -F or --fileFormat Choose the file type to write [PDB/XYZ].
     */
    @CommandLine.Option(names = {"-F", "--fileFormat"}, paramLabel = "XYZ",
            description = "Choose file type to write [PDB/XYZ].")
    private String fileType = "XYZ";

    /**
     * <p>Getter for the field <code>fileType</code>.</p>
     *
     * @return a {@link java.lang.String} object.
     */
    public String getFileType() {
        return fileType;
    }

    public static String toArchiveExtension(String fileType) {
        return Extensions.nameToExt(fileType).archive;
    }

    /**
     * Saves a single-snapshot file to either .xyz or .pdb, depending on the value of fileType.
     * @param baseFileName Basic file name without extension.
     * @param pFuncts      A PotentialFunctions object.
     * @param ma           MolecularAssembly to save.
     * @return             File written to.
     */
    public File saveFile(String baseFileName, PotentialsFunctions pFuncts, MolecularAssembly ma) {
        String outFileName = baseFileName;
        File outFile;
        if (fileType.equalsIgnoreCase("XYZ")) {
            outFileName = outFileName + ".xyz";
            outFile = pFuncts.versionFile(new File(outFileName));

            pFuncts.saveAsXYZ(ma, outFile);
        } else {
            outFileName = outFileName + ".pdb";
            outFile = pFuncts.versionFile(new File(outFileName));

            pFuncts.saveAsPDB(ma, outFile);
        }
        return outFile;
    }

    private enum Extensions {
        XYZ("xyz", "arc"), PDB("pdb", "pdb");

        private final String single;
        private final String archive;

        Extensions(String single, String archive) {
            this.single = single;
            this.archive = archive;
        }

        static Extensions nameToExt(String name) {
            switch(name.toLowerCase()) {
                case "pdb":
                    return PDB;
                case "xyz":
                case "arc":
                default:
                    return XYZ;
            }
        }
    }
}
