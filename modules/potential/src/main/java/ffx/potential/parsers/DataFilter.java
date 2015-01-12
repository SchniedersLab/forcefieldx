/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ffx.potential.parsers;

/**
 *
 * @author JacobLitman
 */
public interface DataFilter {
    public boolean accept(Object ob);
    public String getDescription();
}
