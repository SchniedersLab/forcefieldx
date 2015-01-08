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
package ffx.algorithms;

import java.util.ArrayList;
import java.util.Random;
import java.util.logging.Logger;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.random;

import ffx.numerics.Potential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;

/**
 *
 * @author Jordan
 */
public class Protonate {

    private static final Logger logger = Logger.getLogger(Protonate.class.getName());

    private int stepCount; //the current step
    static int howOften; //how often an MC switch should be attempted

    static int numTitratable;  //the number of titratable residues

    private int numChains; //the number of chains
    private int numResidue; //the number of residue in the current chain

    private final Polymer chains[]; //an array of the chains of the species
    private ArrayList<Residue> residues; //an arraylist of the residues of the chain currently being checked for titratable residues
    private ArrayList<Residue> titratable; //an arraylist of titratable residues

    private Random randomGenerator = new Random(); //a random number generator

    private int numAccepted; //number of accepted MC moves

    private final double boltzmann = 0.0019872041; //kcal/mol/kelvin
    private final double temperature = 298.15;
    private final double kT = boltzmann * temperature;

    public enum titratableList {

        CYS, TYR, HIS, ASP, GLU, LYS, ARG;

    }; //list of residues that have titratable R groups

    //initializes variables and searches species for titratable resiudes
    Protonate(MolecularAssembly molecularassembly, int howoften) {

        //initialize stepcount and the number of accepted moves
        stepCount = 0;
        numAccepted = 0;

        //set how often a MC switch should be attempted
        howOften = howoften;

        //retrieve chains of species from molecular assembly and count the number of chains
        chains = molecularassembly.getChains();
        numChains = chains.length;

        //retrieve residues from a given chain
        for (int i = 0; i != numChains; i++) {
            residues = chains[i].getResidues();
            numResidue = residues.size();

            //check titratability of retrieved residues
            for (int j = 0; j != numResidue; j++) {

                if (CheckTitratability(residues.get(j).getName())) {
                    titratable.add(residues.get(j));
                    numTitratable++;
                }
            }
        }
    }

    //checks if a given residue is contained in the list of titratable residues
    private boolean CheckTitratability(String name) {
        for (titratableList names : titratableList.values()) {
            if (names.name().equals(name)) {
                return true;
            }
        }

        return false;
    }

    public void TryMCSwitch(MolecularAssembly molecularassembly, Potential potential, double InitialEnergy, double pH, double elecref) {
        if (CheckStep()) {
            return;
        }
        //determine which titratable residue should have a protonation state switch attempt
        int Random = randomGenerator.nextInt(numTitratable);
        //DOTHEMIKEYSWITCH(titratable.get(Random));

        //calculate the free energy to be used as the MC acceptance criteria
        double elec = potential.getTotalEnergy() - InitialEnergy;

        double pKaref;

        String name = (titratable.get(Random)).getName();
        switch (name) {
            case "ARG":
                pKaref = 12.48;
                break;
            case "ASP":
                pKaref = 4.0;
                break;
            case "CYS":
                pKaref = 8.18;
                break;
            case "GLU":
                pKaref = 4.25;
                break;
            case "HIS":
                pKaref = 6.00;
                break;
            case "LYS":
                pKaref = 10.53;
                break;
            case "TYR":
                pKaref = 10.07;
                break;
            default:
                pKaref = 7.0;

        }

        double deltaG = kT * (pH - pKaref) * Math.log(10) + elec - elecref;

        //accept the change if the energy change is favorable
        if (deltaG < 0) {
            numAccepted++;
            return;
        }

        //determine if an unenergetically favorable switch should be made
        double boltzmann = exp(-deltaG / kT);
        double metropolis = random();

        if (metropolis < boltzmann) {
            numAccepted++;
            return;
        }

        //DOTHEMIKEYSWITCH(titratable.get(Random));
    }

    //increments stepCount and checks if a monte carlo switch should be attempted
    private boolean CheckStep() {
        stepCount++;
        return (stepCount % howOften == 0);
    }

    //finds the
    public double getAcceptanceRate() {
        double numTries = stepCount / howOften;
        return numAccepted / numTries;
    }
}
