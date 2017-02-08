package ffx.potential.nonbonded;

import java.util.List;
import java.util.logging.Logger;

import ffx.potential.bonded.Angle;
import ffx.potential.bonded.Atom;
import ffx.potential.bonded.Bond;
import ffx.potential.parameters.AtomType;
import ffx.potential.parameters.ForceField;
import ffx.potential.parameters.MultipoleType;
import ffx.potential.parameters.MultipoleType.MultipoleFrameDefinition;
import ffx.potential.parameters.PolarizeType;

/**
 * Bundles an Atom's MultipoleType together with the indices of its axis atoms
 * (ie: those neighbors which feel its torque).
 */
public class Multipole {
    
    private static final Logger logger = Logger.getLogger(Multipole.class.getName());
    
    public final MultipoleType multipoleType;
    public final int[] axisAtom;
    public final MultipoleFrameDefinition frameDefinition;
    
    private Multipole(MultipoleType multipoleType, int[] axisAtom,
            MultipoleFrameDefinition frameDefinition) {
        this.multipoleType = multipoleType;
        this.axisAtom = axisAtom;
        this.frameDefinition = frameDefinition;
    }

    public static Multipole buildMultipole(Atom atom, ForceField forceField) {
        double[] multipole = new double[10];
        int[] axisAtom;
        MultipoleFrameDefinition frame;
        AtomType atomType = atom.getAtomType();
        if (atomType == null) {
            String message = " Multipoles can only be assigned to atoms that have been typed.";
            logger.severe(message);
            return null;
        }

        PolarizeType polarizeType = forceField.getPolarizeType(atomType.getKey());
        if (polarizeType != null) {
            atom.setPolarizeType(polarizeType);
        } else {
            String message = " No polarization type was found for " + atom.toString();
            logger.fine(message);
            double polarizability = 0.0;
            double thole = 0.0;
            int polarizationGroup[] = null;
            polarizeType = new PolarizeType(atomType.type,
                    polarizability, thole, polarizationGroup);
            forceField.addForceFieldType(polarizeType);
            atom.setPolarizeType(polarizeType);
        }

        String key;
        // No reference atoms.
        key = atomType.getKey() + " 0 0";
        MultipoleType multipoleType = forceField.getMultipoleType(key);
        if (multipoleType != null) {
            atom.setMultipoleType(multipoleType);
            return new Multipole(multipoleType, null, multipoleType.frameDefinition);
        }

        // No bonds.
        List<Bond> bonds = atom.getBonds();
        if (bonds == null || bonds.size() < 1) {
            String message = "Multipoles can only be assigned after bonded relationships are defined.\n";
            logger.severe(message);
        }

        // 1 reference atom.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            key = atomType.getKey() + " " + atom2.getAtomType().getKey() + " 0";
            multipoleType = multipoleType = forceField.getMultipoleType(key);
            if (multipoleType != null) {
                int multipoleReferenceAtoms[] = new int[1];
                multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                atom.setMultipoleType(multipoleType);
                axisAtom = multipoleReferenceAtoms;
                frame = multipoleType.frameDefinition;
                return new Multipole(multipoleType, axisAtom, frame);
            }
        }

        // 2 reference atoms.
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                key = atomType.getKey() + " " + key2 + " " + key3;
                multipoleType = forceField.getMultipoleType(key);
                if (multipoleType != null) {
                    int multipoleReferenceAtoms[] = new int[2];
                    multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                    multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                    atom.setMultipoleType(multipoleType);
                    axisAtom = multipoleReferenceAtoms;
                    frame = multipoleType.frameDefinition;
                    return new Multipole(multipoleType, axisAtom, frame);
                }
            }
        }

        /**
         * 3 reference atoms.
         */
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            for (Bond b2 : bonds) {
                if (b == b2) {
                    continue;
                }
                Atom atom3 = b2.get1_2(atom);
                String key3 = atom3.getAtomType().getKey();
                for (Bond b3 : bonds) {
                    if (b == b3 || b2 == b3) {
                        continue;
                    }
                    Atom atom4 = b3.get1_2(atom);
                    String key4 = atom4.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        int multipoleReferenceAtoms[] = new int[3];
                        multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                        multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                        multipoleReferenceAtoms[2] = atom4.getIndex() - 1;
                        atom.setMultipoleType(multipoleType);
                        axisAtom = multipoleReferenceAtoms;
                        frame = multipoleType.frameDefinition;
                        return new Multipole(multipoleType, axisAtom, frame);
                    }
                }
                List<Angle> angles = atom.getAngles();
                for (Angle angle : angles) {
                    Atom atom4 = angle.get1_3(atom);
                    if (atom4 != null) {
                        String key4 = atom4.getAtomType().getKey();
                        key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                        multipoleType = forceField.getMultipoleType(key);
                        if (multipoleType != null) {
                            int multipoleReferenceAtoms[] = new int[3];
                            multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                            multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                            multipoleReferenceAtoms[2] = atom4.getIndex() - 1;
                            atom.setMultipoleType(multipoleType);
                            axisAtom = multipoleReferenceAtoms;
                            frame = multipoleType.frameDefinition;
                            return new Multipole(multipoleType, axisAtom, frame);
                        }
                    }
                }
            }
        }

        /**
         * Revert to a 2 reference atom definition that may include a 1-3 site.
         * For example a hydrogen on water.
         */
        for (Bond b : bonds) {
            Atom atom2 = b.get1_2(atom);
            String key2 = atom2.getAtomType().getKey();
            List<Angle> angles = atom.getAngles();
            for (Angle angle : angles) {
                Atom atom3 = angle.get1_3(atom);
                if (atom3 != null) {
                    String key3 = atom3.getAtomType().getKey();
                    key = atomType.getKey() + " " + key2 + " " + key3;
                    multipoleType = forceField.getMultipoleType(key);
                    if (multipoleType != null) {
                        int multipoleReferenceAtoms[] = new int[2];
                        multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                        multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                        atom.setMultipoleType(multipoleType);
                        axisAtom = multipoleReferenceAtoms;
                        frame = multipoleType.frameDefinition;
                        return new Multipole(multipoleType, axisAtom, frame);
                    }
                    for (Angle angle2 : angles) {
                        Atom atom4 = angle2.get1_3(atom);
                        if (atom4 != null && atom4 != atom3) {
                            String key4 = atom4.getAtomType().getKey();
                            key = atomType.getKey() + " " + key2 + " " + key3 + " " + key4;
                            multipoleType = forceField.getMultipoleType(key);
                            if (multipoleType != null) {
                                int multipoleReferenceAtoms[] = new int[3];
                                multipoleReferenceAtoms[0] = atom2.getIndex() - 1;
                                multipoleReferenceAtoms[1] = atom3.getIndex() - 1;
                                multipoleReferenceAtoms[2] = atom4.getIndex() - 1;
                                atom.setMultipoleType(multipoleType);
                                axisAtom = multipoleReferenceAtoms;
                                frame = multipoleType.frameDefinition;
                                return new Multipole(multipoleType, axisAtom, frame);
                            }
                        }
                    }
                }
            }
        }
        return null;
    }
    
}
