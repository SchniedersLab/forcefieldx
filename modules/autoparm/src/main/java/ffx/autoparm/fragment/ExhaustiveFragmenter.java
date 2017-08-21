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
package ffx.autoparm.fragment;

/* Copyright (C) 2010  Rajarshi Guha <rajarshi.guha@gmail.com>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.SpanningTree;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * Generate fragments exhaustively.
 * <p/>
 * This fragmentation scheme simply breaks single non-ring bonds. By default
 * fragments smaller than 6 atoms in size are not considered, but this can be
 * changed by the user. Side chains are retained.
 *
 * @author Rajarshi Guha
 * @cdk.module fragment
 * @cdk.githash
 * @cdk.keyword fragment
 */
/**
 * Potential edits made to fit FFX specifications
 *
 * @author rcorrigan
 */
public class ExhaustiveFragmenter implements IFragmenter {

    private static final int DEFAULT_MIN_FRAG_SIZE = 6;

    Map<String, IAtomContainer> fragMap;
    SmilesGenerator smilesGenerator;
    String[] fragments = null;
    int minFragSize = 6;
    private static ILoggingTool logger = LoggingToolFactory
            .createLoggingTool(ExhaustiveFragmenter.class);

    /**
     * Instantiate fragmenter with default minimum fragment size.
     */
    public ExhaustiveFragmenter() {
        this(DEFAULT_MIN_FRAG_SIZE);
    }

    /**
     * Instantiate fragmenter with user specified minimum fragment size.
     *
     * @param minFragSize the minimum fragment size desired
     */
    public ExhaustiveFragmenter(int minFragSize) {
        this.minFragSize = minFragSize;
        fragMap = new HashMap<>();
        smilesGenerator = SmilesGenerator.unique().aromatic();
    }

    /**
     * Set the minimum fragment size.
     *
     * @param minFragSize the smallest size fragment that will be returned
     */
    public void setMinimumFragmentSize(int minFragSize) {
        this.minFragSize = minFragSize;
    }

    /**
     * Generate fragments for the input molecule.
     *
     * @param atomContainer The input molecule.
     */
    @Override
    public void generateFragments(IAtomContainer atomContainer) throws CDKException {
        fragMap.clear();
        run(atomContainer);
    }

    private List<IAtomContainer> run(IAtomContainer atomContainer) throws CDKException {

        ArrayList<IAtomContainer> fragments = new ArrayList<>();

        if (atomContainer.getBondCount() < 3) {
            return fragments;
        }
        List<IBond> splitableBonds = getSplitableBonds(atomContainer);
        if (splitableBonds.isEmpty()) {
            return fragments;
        }
        logger.debug("Got " + splitableBonds.size() + " splittable bonds");

        String tmpSmiles;
        String[] uniqueAtomNames;
        String[] uniqueAtomNames2;
        boolean flag = true;
        for (IBond bond : splitableBonds) {
            List<IAtomContainer> parts = FragmentUtils.splitMolecule(atomContainer, bond);
            // make sure we don't add the same fragment twice
            //try using regular for loop to ensure atoms/fragments are in correct order
            //for (IAtomContainer partContainer : parts) {
            for (int i = 0; i < parts.size(); i++) {
                IAtomContainer partContainer = parts.get(i);
                uniqueAtomNames = new String[partContainer.getAtomCount()];
                //go thru partContainer, get all atoms, record and store atom IDs for later reassignment
                // HashMap<IAtom, String> hashMap = new HashMap<IAtom, String>();
                for (int j = 0; j < partContainer.getAtomCount(); j++) {
                    //atoms names recorded
                    IAtom atom = partContainer.getAtom(j);
                    String name = atom.getID();
                    uniqueAtomNames[j] = atom.getID();
                }
                AtomContainerManipulator.clearAtomConfigurations(partContainer);
                for (IAtom atom : partContainer.atoms()) {
                    atom.setImplicitHydrogenCount(null);
                }
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(partContainer);
                CDKHydrogenAdder.getInstance(partContainer.getBuilder()).addImplicitHydrogens(partContainer);
                Aromaticity.cdkLegacy().apply(partContainer);
                tmpSmiles = smilesGenerator.create(partContainer);
                if (partContainer.getAtomCount() >= minFragSize && !fragMap.containsKey(tmpSmiles)) {
                    fragments.add(partContainer);
                    fragMap.put(tmpSmiles, partContainer);
                }
            }
            flag = false;
        }

        // try and partition the fragments
        List<IAtomContainer> tmp = new ArrayList<>(fragments);
        for (IAtomContainer fragment : fragments) {
            if (fragment.getBondCount() < 3 || fragment.getAtomCount() < minFragSize) {
                continue;
            }
            if (getSplitableBonds(fragment).isEmpty()) {
                continue;
            }

            List<IAtomContainer> frags = run(fragment);
            if (frags.isEmpty()) {
                continue;
            }

            //for (IAtomContainer frag : frags) {
            for (int x = 0; x < frags.size(); x++) {
                IAtomContainer frag = frags.get(x);
                //sets up uniqueAtomsNames2 array and fills it with unique atom ID's from fragment atoms
                uniqueAtomNames2 = new String[frag.getAtomCount()];
                for (int y = 0; y < frag.getAtomCount(); y++) {
                    uniqueAtomNames2[y] = frag.getAtom(y).getID();
                    //System.out.println("UniqueAtomNames2: " + uniqueAtomNames2[y]);
                }
                if (frag.getBondCount() < 3) {
                    continue;
                }
                AtomContainerManipulator.clearAtomConfigurations(frag);
                for (IAtom atom : frag.atoms()) {
                    atom.setImplicitHydrogenCount(null);
                }
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(frag);
                CDKHydrogenAdder.getInstance(frag.getBuilder()).addImplicitHydrogens(frag);
                Aromaticity.cdkLegacy().apply(frag);
                tmpSmiles = smilesGenerator.create(frag);
                for (int z = 0; z < frag.getAtomCount(); z++) {
                    IAtom test = frag.getAtom(z);
                }
                //***Atom Types and unique names/ID's match here!***

                //System.out.println("\nFinished with fragment x = "+x);
                if (frag.getAtomCount() >= minFragSize && !fragMap.containsKey(tmpSmiles)) {
                    tmp.add(frag);
                    fragMap.put(tmpSmiles, frag);
                }
            }
        }
        fragments = new ArrayList<>(tmp);
        return fragments;
    }

    private List<IBond> getSplitableBonds(IAtomContainer atomContainer) throws CDKException {
        // do ring detection
        SpanningTree spanningTree = new SpanningTree(atomContainer);
        IRingSet allRings = spanningTree.getAllRings();

        // find the splitable bonds
        ArrayList<IBond> splitableBonds = new ArrayList<>();

        for (IBond bond : atomContainer.bonds()) {
            boolean isInRing = false;
            boolean isTerminal = false;

            // lets see if it's in a ring
            IRingSet rings = allRings.getRings(bond);
            if (rings.getAtomContainerCount() != 0) {
                isInRing = true;
            }

            // lets see if it is a terminal bond
            for (IAtom atom : bond.atoms()) {
                if (atomContainer.getConnectedAtomsCount(atom) == 1) {
                    isTerminal = true;
                    break;
                }
            }

            if (!(isInRing || isTerminal)) {
                splitableBonds.add(bond);
            }
        }
        return splitableBonds;
    }

    /**
     * Get the fragments generated as SMILES strings.
     *
     * @return a String[] of the fragments.
     */
    @Override
    public String[] getFragments() {
        return (new ArrayList<>(fragMap.keySet())).toArray(new String[0]);
    }

    /**
     * Get the fragments generated as {@link IAtomContainer} objects..
     *
     * @return a IAtomContainer[] of the fragments.
     */
    @Override
    public IAtomContainer[] getFragmentsAsContainers() {
        return (new ArrayList<>(fragMap.values())).toArray(new IAtomContainer[0]);
    }

}
