// TIMER
import ffx.potential.ForceFieldEnergy;

// Name of the file (PDB or XYZ).
String filename = args[0];
if (filename == null) {
   println("\n Usage: ffxc timer filename nEvals");
   return;
}

nEvals = 5;
if (args.size() > 1) {
    nEvals = Integer.parseInt(args[1]);
}

// Things below this line normally do not need to be changed.
// ===============================================================================================

println("\n Timing energy and gradient for " + filename);

open(filename);

ForceFieldEnergy energy = active.getPotentialEnergy();
for (int i=0; i<nEvals; i++) {
    energy.energy(true, true);
}