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
package ffx.algorithms;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Residue;
import ffx.potential.bonded.ResidueEnumerations.AminoAcid3;
import ffx.potential.bonded.ResidueEnumerations.NucleicAcid3;
import ffx.potential.bonded.Rotamer;
import ffx.potential.bonded.RotamerLibrary;
import ffx.potential.parsers.PDBFilter;

/**
 * The GenerateRotamers class helps generate a rotamer library (particularly for 
 * nonstandard amino acids) for a Residue.
 *
 * @author Jacob M. Litman
 * @author Michael J. Schnieders
 */
public class GenerateRotamers {
    private static final Logger logger = Logger.getLogger(GenerateRotamers.class.getName());
    
    private final Residue residue;
    private final int nChi;
    private final double[] currentChi;
    
    private final File outFile;
    private final MolecularAssembly mola;
    private final AlgorithmListener listener;
    private final RotamerLibrary library;
    
    private final Potential potential;
    private double[] x;
    
    private AminoAcid3 baselineAAres;
    private Rotamer[] baselineRotamers;
    
    private double incr = 10.0; // Degrees to rotate each torsion by.
    private boolean aroundLibrary = false;
    private double width = 180.0; // Left 180 + right 180 = 360 degree coverage.
    
    private int startDepth = 0;
    private int endDepth;
    
    private boolean print = false; // If true, log energy at each torsion.
    private int nEval = 0;
    
    private boolean writeVideo = false;
    private File videoFile;
    private PDBFilter videoFilter;
    
    private static final Pattern atRangePatt = Pattern.compile("(\\d+)-(\\d+)");
    
    /**
     * Intended to create rotamer sets for nonstandard amino acids.
     * @param mola
     * @param potential
     * @param residue
     * @param file Output file
     * @param nChi Number of rotameric torsions in the residue
     * @param listener
     */
    public GenerateRotamers(MolecularAssembly mola, Potential potential, 
            Residue residue, File file, int nChi, AlgorithmListener listener) {
        this(mola, potential, residue, file, nChi, listener, RotamerLibrary.getDefaultLibrary());
    }
    
    /**
     * Intended to create rotamer sets for nonstandard amino acids.
     * @param mola
     * @param potential
     * @param residue
     * @param file Output file
     * @param nChi Number of rotameric torsions in the residue
     * @param listener
     * @param library Rotamer library to use
     */
    public GenerateRotamers(MolecularAssembly mola, Potential potential, 
            Residue residue, File file, int nChi, AlgorithmListener listener, 
            RotamerLibrary library) {
        this.residue = residue;
        this.nChi = nChi;
        endDepth = nChi - 1;
        this.potential = potential;
        this.mola = mola;
        this.listener = listener;
        this.currentChi = new double[nChi];
        this.library = library;
        Arrays.fill(currentChi, 0.0);
        File outputFile = file;
        if (outputFile.exists()) {
            String outName = outputFile.getName();
            for (int i = 1; i < 1000; i++) {
                outputFile = new File(String.format("%s_%d", outName, i));
                if (!outputFile.exists()) {
                    break;
                }
            }
            if (outputFile.exists()) {
                logger.severe(String.format(" Could not version file %s", outName));
            }
        }
        outFile = outputFile;
        this.baselineAAres = residue.getAminoAcid3();
        baselineRotamers = library.getRotamers(baselineAAres);
    }
    
    /**
     * Sets a standard amino acid to be the baseline for rotamer generation. For
     * example, use TYR as a baseline for phosphotyrosine.
     * @param aa3 
     */
    public void setBaselineAARes(AminoAcid3 aa3) {
        this.baselineAAres = aa3;
        if (aa3 == aa3.UNK) {
            boolean orig = library.getUsingOrigCoordsRotamer();
            library.setUseOrigCoordsRotamer(false);
            baselineRotamers = residue.getRotamers(library);
            library.setUseOrigCoordsRotamer(orig);
        } else {
            baselineRotamers = library.getRotamers(aa3);
        }
        aroundLibrary = true;
    }
    
    /**
     * Sets algorithm to log all torsions/energies (not just to file).
     * @param print 
     */
    public void setPrint(boolean print) {
        this.print = print;
    }
    
    /**
     * Set which torsions to work on. Negative end values set the final depth
     * to be the total number of torsions.
     * @param start
     * @param end
     */
    public void setDepth(int start, int end) {
        startDepth = start;
        if (end < 1 || end > nChi) {
            endDepth = nChi - 1;
        } else {
            endDepth = end;
        }
    }
    
    /**
     * Sets the width around each torsion to search (+/-, so 10 degree width will
     * search a 20 degree arc).
     * @param width 
     */
    public void setSearchWidth(double width) {
        this.width = width;
    }
    
    /**
     * Sets the angle to change torsions by.
     * @param incr 
     */
    public void setIncrement(double incr) {
        this.incr = incr;
    }
    
    /**
     * Null file indicates to not write a video.
     * @param videoFile Filename for video or null
     */
    public void setVideo(String videoFile) {
        if (videoFile != null) {
            File vidFile = new File(videoFile);
            if (vidFile.exists()) {
                for (int i = 0; i < 1000; i++) {
                    vidFile = new File(String.format("%s_%d", videoFile, i));
                    if (!vidFile.exists()) {
                        this.videoFile = vidFile;
                        writeVideo = true;
                        videoFilter = new PDBFilter(this.videoFile, mola, mola.getForceField(), null);
                        videoFilter.setLogWrites(false);
                        break;
                    }
                }
                if (vidFile.exists()) {
                    logger.warning(String.format(" Could not version video file %s", videoFile));
                }
            } else {
                this.videoFile = vidFile;
                writeVideo = true;
                videoFilter = new PDBFilter(this.videoFile, mola, mola.getForceField(), null);
                videoFilter.setLogWrites(false);
            }
        } else {
            writeVideo = false;
        }
    }
    
    /**
     * Inactivates electrostatics for atom sets defined by 'start-end,start-end,...'.
     * @param electrostatics Input string
     */
    public void setElectrostatics(String electrostatics) {
        if (electrostatics != null) {
            String[] toks = electrostatics.split(",");
            Atom[] atoms = mola.getAtomArray();
            for (String tok : toks) {
                Matcher m = atRangePatt.matcher(tok);
                if (m.matches()) {
                    int begin = Integer.parseInt(m.group(1));
                    int end = Integer.parseInt(m.group(2));
                    logger.info(String.format(" Inactivating electrostatics for atoms %d-%d", begin, end));
                    for (int i = begin; i <= end; i++) {
                        Atom ai = atoms[i - 1];
                        ai.setElectrostatics(false);
                        ai.print();
                    }
                } else {
                    logger.info(String.format(" Discarding electrostatics input %s", tok));
                }
            }
        }
    }
    
    /**
     * Inactivates atom sets defined by 'start-end,start-end,...'.
     *
     * @param iatoms Input string
     */
    public void setInactiveAtoms(String iatoms) {
        if (iatoms != null) {
            String[] toks = iatoms.split(",");
            Atom[] atoms = mola.getAtomArray();
            for (String tok : toks) {
                Matcher m = atRangePatt.matcher(tok);
                if (m.matches()) {
                    int begin = Integer.parseInt(m.group(1));
                    int end = Integer.parseInt(m.group(2));
                    logger.info(String.format(" Inactivating atoms %d-%d", begin, end));
                    for (int i = begin; i <= end; i++) {
                        Atom ai = atoms[i - 1];
                        ai.setUse(false);
                        ai.print();
                    }
                } else {
                    logger.info(String.format(" Discarding inactive atoms input %s", tok));
                }
            }
        }
    }
    
    /**
     * Main driver method; spins torsions, evaluates energy, and prints to file.
     */
    public void tryRotamers() {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outFile))) {
            if (aroundLibrary) {
                for (Rotamer rotamer : baselineRotamers) {
                    applyLibRotamer(rotamer);
                    turnChi(startDepth, bw);
                }
            } else {
                RotamerLibrary.measureRotamer(residue, currentChi, false);
                turnChi(startDepth, bw);
            }
        } catch (IOException ex) {
            logger.warning(String.format(" IO exception in rotamer generation: %s", ex.toString()));
        }
    }
    
    /**
     * Applies a library rotamer, filling in the end with 0s as necessary. For 
     * example, PTR requires more torsions than TYR, so one cannot simply apply
     * a TYR rotamer to a PTR residue.
     * @param rotamer 
     */
    private void applyLibRotamer(Rotamer rotamer) {
        double[] rotValues = new double[nChi*2];
        Arrays.fill(rotValues, 0.0);
        Arrays.fill(currentChi, 0.0);
        int fillTo = Math.min(startDepth, (rotamer.length - 1));
        for (int i = 0; i < fillTo; i++) {
            int ii = 2*i;
            rotValues[ii] = rotamer.angles[i];
            currentChi[i] = rotamer.angles[i];
            rotValues[ii+1] = rotamer.sigmas[i];
        }
        Rotamer newRot = generateRotamer(rotValues);
        RotamerLibrary.applyRotamer(residue, newRot);
    }
    
    /**
     * Generates an aa/na/unk Rotamer for the selected residue given values.
     * @param values
     * @return 
     */
    private Rotamer generateRotamer(double[] values) {
        switch (residue.getResidueType()) {
            case AA:
                AminoAcid3 aa3 = residue.getAminoAcid3();
                return new Rotamer(aa3, values);
            case NA:
                NucleicAcid3 na3 = residue.getNucleicAcid3();
                return new Rotamer(na3, values);
            case UNK:
            default:
                return new Rotamer(values);
        }
    }
    
    /**
     * Recursive method for turning the torsions.
     * @param depth Current depth of the recursion
     * @param bw
     * @throws IOException 
     */
    private void turnChi(int depth, BufferedWriter bw) throws IOException {
        double chi = currentChi[depth];
        for (double i = chi - width; i <= chi + width; i+= incr) {
            currentChi[depth] = i;
            if (depth == endDepth) {
                evaluateChi(bw);
                if (listener != null) {
                    listener.algorithmUpdate(mola);
                }
                if (writeVideo) {
                    writeSnapshot();
                }
            } else {
                turnChi(depth+1, bw);
            }
        }
        currentChi[depth] = chi; // Reset the chi value to where it was.
    }
    
    private void writeSnapshot() {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(videoFile, true))) {
            bw.write(String.format("MODEL %d", nEval));
            bw.newLine();
            bw.write(String.format("REMARK 301 TORSIONS %s", formatChi()));
            bw.newLine();
            bw.flush();
            videoFilter.writeFile(videoFile, true);
            bw.write("ENDMDL");
            bw.newLine();
            bw.flush();
        } catch (IOException ex) {
            logger.warning(String.format(" Exception writing to video file %s: "
                    + "%s", videoFile.getName(), ex.toString()));
        }
    }
    
    /**
     * Called at a leaf of the recursion.
     * @param bw
     * @throws IOException 
     */
    private void evaluateChi(BufferedWriter bw) throws IOException {
        try {
            applyChi();
            double e = currentEnergy();
            String result = String.format("%s,%12f", formatChi(), e);
            ++nEval;
            if (print) {
                logger.info(String.format(" Evaluation %10d %s, energy %10.5f kcal/mol", nEval, formatChi(), e));
            }
            bw.write(result);
            bw.newLine();
            if (nEval % 1000 == 0) {
                logger.info(String.format(" %12.7e states evaluated", (double) nEval));
            }
        } catch (ArithmeticException ex) {
            logger.info(String.format(" Force field exception at chi %s", formatChi()));
        }
    }
    
    /**
     * Converts the chi array to a Rotamer and applies it.
     */
    private void applyChi() {
        double[] rotValues = new double[nChi*2];
        for (int i = 0; i < nChi; i++) {
            int ii = 2*i;
            rotValues[ii] = currentChi[i];
            rotValues[ii+1] = 0.0;
        }
        RotamerLibrary.applyRotamer(residue, generateRotamer(rotValues));
    }
    
    /**
     * Returns a formatted String with chi values.
     * @return 
     */
    private String formatChi() {
        StringBuilder sb = new StringBuilder(String.format("%8f", currentChi[0]));
        for (int i = 1; i < nChi; i++) {
            //sb.append(",").append(currentChi[i]);
            sb.append(String.format(",%8f", currentChi[i]));
        }
        return sb.toString();
    }
    
    /**
     * Accessory method for more simplistic saving of specific torsion states.
     * @param torSets 
     */
    public void applyAndSaveTorsions(String[] torSets) {
        for (String torSet : torSets) {
            String[] torsions = torSet.split(",");
            double[] values = new double[nChi*2];
            Arrays.fill(values, 0.0);
            for (int i = 0; i < (Math.min(torsions.length, nChi)); i++) {
                double chival = Double.parseDouble(torsions[i]);
                currentChi[i] = chival;
                values[2*i] = chival;
            }
            Rotamer newRot = generateRotamer(values);
            RotamerLibrary.applyRotamer(residue, newRot);
            writeSnapshot();
        }
    }

    /**
     * Calculates the energy at the current state.
     *
     * @return Energy of the current state.
     */
    private double currentEnergy() throws ArithmeticException {
        if (x == null) {
            int nVar = potential.getNumberOfVariables();
            x = new double[nVar * 3];
            
        }
        potential.getCoordinates(x);
        return potential.energy(x);
    }
}
