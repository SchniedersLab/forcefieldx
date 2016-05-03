/**
 * Title: Force Field X.
 *
 * Description: Force Field X - Software for Molecular Biophysics.
 *
 * Copyright: Copyright (c) Michael J. Schnieders 2001-2016.
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

import ffx.potential.ForceFieldEnergy;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Atom;
import ffx.potential.parsers.SystemFilter;
import ffx.potential.utils.PotentialsFunctions;
import java.util.logging.Logger;

/**
 * Evaluate binding affinity from an MD trajectory using MM-GKSA/MM-GBSA/MM-PBSA
 * approximation. Binding affinity is mean value of (complex - protein - ligand),
 * often weighted by component.
 * 
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class MMgksa {
    private static final Logger logger = Logger.getLogger(MMgksa.class.getName());
    private final MolecularAssembly mola;
    private final PotentialsFunctions functs;
    private final SystemFilter filter;
    private final ForceFieldEnergy energy;
    private final Atom[] proteinAtoms;
    private final Atom[] ligandAtoms;
    private final boolean docking;
    private Atom[] ignoreAtoms;
    private double elecWt = 1.0;
    private double solvWt = 1.0;
    private double vdwWt = 1.0;
    
    public MMgksa(MolecularAssembly mola, PotentialsFunctions functs, SystemFilter filter) {
        this.mola = mola;
        this.functs = functs;
        this.filter = filter;
        this.energy = mola.getPotentialEnergy();
        proteinAtoms = new Atom[0];
        ligandAtoms = new Atom[0];
        docking = false;
    }
    
    public MMgksa(MolecularAssembly mola, PotentialsFunctions functs, SystemFilter filter, 
            Atom[] proteinAtoms, Atom[] ligandAtoms) {
        this.mola = mola;
        this.functs = functs;
        this.filter = filter;
        this.energy = mola.getPotentialEnergy();
        this.proteinAtoms = new Atom[proteinAtoms.length];
        System.arraycopy(proteinAtoms, 0, this.proteinAtoms, 0, proteinAtoms.length);
        this.ligandAtoms = new Atom[ligandAtoms.length];
        System.arraycopy(ligandAtoms, 0, this.ligandAtoms, 0, ligandAtoms.length);
        docking = true;
    }
    
    public void setIgnoredAtoms(Atom[] ignored) {
        ignoreAtoms = new Atom[ignored.length];
        System.arraycopy(ignored, 0, ignoreAtoms, 0, ignored.length);
        for (Atom atom : ignoreAtoms) {
            atom.setUse(false);
        }
    }
    
    public void setElectrostaticsWeight(double eWt) {
        this.elecWt = eWt;
    }
    
    public void setSolvationWeight(double sWt) {
        this.solvWt = sWt;
    }
    
    public void setVdwWeight(double vWt) {
        this.vdwWt = vWt;
    }
    
    /**
     * Evaluates binding or average protein energy.
     * @param frequency Evaluate energy every X frames
     */
    public void runMMgksa(int frequency) {
        runMMgksa(frequency, -1);
    }
    
    /**
     * Evaluates binding or average protein energy.
     * @param frequency Evaluate energy every X frames
     * @param maxEvals Negative value runs to end of trajectory.
     */
    public void runMMgksa(int frequency, int maxEvals) {
        double meanElecE = 0.0;
        double meanSolvE = 0.0;
        double meanVdwE = 0.0;
        double meanE = 0.0;
        int nEvals = 0;
        // Total frames read.
        int counter = 0;
        
        do {
            if (counter++ % frequency != 0) {
                continue;
            }
            ++nEvals;
            functs.energy(mola);
            double totElecE = energy.getElectrostaticEnergy();
            double totSolvE = energy.getSolvationEnergy();
            double totVdwE = energy.getVanDerWaalsEnergy();
            double totE = weightEnergy(totElecE, totSolvE, totVdwE);
            
            if (docking) {
                for (Atom atom : proteinAtoms) {
                    atom.setUse(false);
                }
                functs.energy(mola);
                double ligElecE = energy.getElectrostaticEnergy();
                double ligSolvE = energy.getSolvationEnergy();
                double ligVdwE = energy.getVanDerWaalsEnergy();
                double ligE = weightEnergy(ligElecE, ligSolvE, ligVdwE);
                for (Atom atom : proteinAtoms) {
                    atom.setUse(true);
                }
                
                for (Atom atom : ligandAtoms) {
                    atom.setUse(false);
                }
                functs.energy(mola);
                double protElecE = energy.getElectrostaticEnergy();
                double protSolvE = energy.getSolvationEnergy();
                double protVdwE = energy.getVanDerWaalsEnergy();
                double protE = weightEnergy(protElecE, protSolvE, protVdwE);
                for (Atom atom : ligandAtoms) {
                    atom.setUse(true);
                }
                
                double elec = totElecE - ligElecE - protElecE;
                double solv = totSolvE - ligSolvE - protSolvE;
                double vdw = totVdwE - ligVdwE - protVdwE;
                double e = totE - ligE - protE;
                
                meanElecE += (elec - meanElecE) / nEvals;
                meanSolvE += (solv - meanSolvE) / nEvals;
                meanVdwE += (vdw - meanVdwE) / nEvals;
                meanE += (e - meanE) / nEvals;
                
                logger.info(String.format(" %-10d frames read, %-10d frames "
                        + "evaluated", counter, nEvals));
                logger.info(formatHeader());
                logger.info(formatEnergy("Running Mean", meanE, meanElecE, meanSolvE, meanVdwE));
                logger.info(formatEnergy("Snapshot", e, elec, solv, vdw));
                logger.info(formatEnergy("Complex", totE, totElecE, totSolvE, totVdwE));
                logger.info(formatEnergy("Protein", protE, protElecE, protSolvE, protVdwE));
                logger.info(formatEnergy("Ligand", ligE, ligElecE, ligSolvE, ligVdwE));
                /*logger.info(String.format("%15s%15s%15s%15s%15s", " ", "Energy", 
                        "Electrostatic", "Solvation", "van der Waals"));
                logger.info(String.format("%-15s  %13.5f  %13.5f  %13.5f  %13.5f", 
                        "Running Mean", meanE, meanElecE, meanSolvE, meanVdwE));
                logger.info(String.format("%-15s  %13.5f  %13.5f  %13.5f  %13.5f", 
                        "Snapshot", e, elec, solv, vdw));
                logger.info(String.format("%-15s  %13.5f  %13.5f  %13.5f  %13.5f", 
                        "Complex", totE, totElecE, totSolvE, totVdwE));
                logger.info(String.format("%-15s  %13.5f  %13.5f  %13.5f  %13.5f", 
                        "Protein", protE, protElecE, protSolvE, protVdwE));
                logger.info(String.format("%-15s  %13.5f  %13.5f  %13.5f  %13.5f", 
                        "Ligand", ligE, ligElecE, ligSolvE, ligVdwE));*/
            } else {
                meanElecE += (totElecE - meanElecE) / nEvals;
                meanSolvE += (totSolvE - meanSolvE) / nEvals;
                meanVdwE += (totVdwE - meanVdwE) / nEvals;
                meanE += (totE - meanE) / nEvals;
                
                logger.info(String.format(" %10d frames read, %10d frames "
                        + "evaluated", counter, nEvals));
                logger.info(String.format("%15c%15s%15s%15s%15s", ' ', "Energy",
                        "Electrostatic", "Solvation", "van der Waals"));
                logger.info(String.format("%-15s  %13.5f  %13.5f  %13.5f  %13.5f", 
                        "Running Mean", meanE, meanElecE, meanSolvE, meanVdwE));
                logger.info(String.format("%-15s  %13.5f  %13.5f  %13.5f  %13.5f", 
                        "Snapshot", totE, totElecE, totSolvE, totVdwE));
            }
        } while (filter.readNext() && (maxEvals < 0 || nEvals < maxEvals));
        
        logger.info(String.format(" MM-GKSA evaluation complete, %10d frames "
                + "read, %11.6f kcal/mol mean energy", nEvals, meanE));
    }
    
    private String formatHeader() {
        StringBuilder sb = new StringBuilder(String.format("%15s%15s", " ", "Energy"));
        if (elecWt != 0.0) {
            sb.append(String.format("%15s", "Electrostatic"));
        }
        if (solvWt != 0.0) {
            sb.append(String.format("%15s", "Solvation"));
        }
        if (vdwWt != 0.0) {
            sb.append(String.format("%15s", "van der Waals"));
        }
        return sb.toString();
    }
    
    private String formatEnergy(String title, double total, double electro, double solv, double vdW) {
        StringBuilder sb = new StringBuilder(String.format("%-15s", title));
        sb.append(String.format("  %13.5f", total));
        if (elecWt != 0.0) {
            sb.append(String.format("  %13.5f", electro));
        }
        if (solvWt != 0.0) {
            sb.append(String.format("  %13.5f", solv));
        }
        if (vdwWt != 0.0) {
            sb.append(String.format("  %13.5f", vdW));
        }
        return sb.toString();
    }
    
    public double weightEnergy(double elec, double solv, double vdW) {
        return (elecWt * elec) + (solvWt * solv) + (vdwWt * vdW);
    }
}
