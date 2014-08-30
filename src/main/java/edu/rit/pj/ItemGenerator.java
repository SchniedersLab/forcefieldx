//******************************************************************************
//
// File:    ItemGenerator.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.ItemGenerator
//
// This Java source file is copyright (C) 2007 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj;

/**
 * Class ItemGenerator is the abstract base class for an object that generates
 * items for a {@linkplain ParallelIteration}.
 *
 * @param <T> Data type of the items iterated over.
 *
 * @author Alan Kaminsky
 * @version 04-Jun-2007
 */
abstract class ItemGenerator<T> {

// Hidden data members.
    // Current item sequence number for filling in the item holder.
    int myCurrentSequenceNumber;

    // Current iteration index for use by the ordered() construct.
    int myOrderedIndex;

    // True to break out of the parallel iteration.
    boolean myBreak;

// Exported constructors.
    /**
     * Construct a new item generator.
     */
    public ItemGenerator() {
    }

// Exported operations.
    /**
     * Return an item holder containing the next item to be processed plus
     * associated information.
     *
     * @return Item holder, or null if no more items.
     */
    public abstract ItemHolder<T> nextItem();

}
