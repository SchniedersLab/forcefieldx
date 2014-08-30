/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.algorithms;

import ffx.numerics.Potential;
import ffx.potential.bonded.MolecularAssembly;
import ffx.potential.bonded.Polymer;
import ffx.potential.bonded.Residue;

import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.random;

import java.util.ArrayList;
import java.util.Random;
import java.util.logging.Logger;

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
