// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2021.
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
// ******************************************************************************
package ffx.utilities;

import org.apache.commons.configuration2.CompositeConfiguration;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * This represents the context of an FFX MolecularAssembly.
 * <p>
 * CompositeConfiguration is extended by adding the API used by the Groovy Binding class.
 * <p>
 * The Binding API represents the variable bindings of a script which can be altered from outside the
 * script object or created outside of a script and passed into it.
 * <p>
 * The Binding API is not supposed to be used in a multithreaded context.
 */
public class FFXContext extends CompositeConfiguration {

  private Map<String, Object> variables;

  public FFXContext() {
    super();
  }

  public FFXContext(Map<String, Object> variables) {
    this();
    this.variables = variables;
  }

  /**
   * A helper constructor used in main(String[]) method calls
   *
   * @param args are the command line arguments from a main()
   */
  public FFXContext(String[] args) {
    this();
    setVariable("args", args);
  }

  /**
   * @param name the name of the variable to lookup
   * @return the variable value
   */
  public Object getVariable(String name) throws Exception {
    if (variables == null) {
      throw new Exception(name);
    }

    Object result = variables.get(name);

    if (result == null && !variables.containsKey(name)) {
      throw new Exception(name);
    }

    return result;
  }

  /**
   * Sets the value of the given variable
   *
   * @param name the name of the variable to set
   * @param value the new value for the given variable
   */
  public void setVariable(String name, Object value) {
    if (variables == null) {
      variables = new LinkedHashMap<>();
    }
    variables.put(name, value);
  }

  /**
   * Remove the variable with the specified name.
   *
   * @param name the name of the variable to remove
   */
  public void removeVariable(String name) {
    if (null == variables) {
      return;
    }

    variables.remove(name);
  }

  /**
   * Simple check for whether the context contains a particular variable or not.
   *
   * @param name the name of the variable to check for
   */
  public boolean hasVariable(String name) {
    return variables != null && variables.containsKey(name);
  }

  /**
   * Get the Map of all variables. This returns a reference and not a copy.
   *
   * @return a Map of the variables.
   */
  public Map<String, Object> getVariables() {
    if (variables == null) {
      variables = new LinkedHashMap<>();
    }
    return variables;
  }

  /**
   * Set the Map of all variables. The reference is kept without making a copy.
   *
   * @param variables a Map of the variables.
   */
  public void setVariables(Map<String, Object> variables) {
    this.variables = variables;
  }

}

