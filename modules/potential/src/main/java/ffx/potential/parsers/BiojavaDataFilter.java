/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.parsers;

import org.biojava.bio.structure.Structure;

/**
 *
 * @author JacobLitman
 */
public class BiojavaDataFilter implements DataFilter {

    @Override
    public boolean accept(Object ob) {
        return ob instanceof Structure;
    }

    @Override
    public String getDescription() {
        return "Biojava 3 Structure";
    }
    
}
