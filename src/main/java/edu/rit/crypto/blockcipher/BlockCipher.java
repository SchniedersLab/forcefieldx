//******************************************************************************
//
// File:    BlockCipher.java
// Package: edu.rit.crypto.blockcipher
// Unit:    Class edu.rit.crypto.blockcipher.BlockCipher
//
// This Java source file is copyright (C) 2005 by Alan Kaminsky. All rights
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

package edu.rit.crypto.blockcipher;

/**
 * Class BlockCipher is the abstract base class for a block cipher. A block
 * cipher object can encrypt a block of plaintext, yielding a block of
 * ciphertext. (Decryption is defined in another class.) The encryption key is
 * specified when the block cipher object is constructed and may also be
 * specified by calling the <TT>setKey()</TT> method. A block cipher object can
 * report the number of bytes in the block (the block length) and the number of
 * bytes in the key (the key length).
 * <P>
 * The actual encryption algorithm is implemented in a subclass.
 *
 * @author  Alan Kaminsky
 * @version 06-Dec-2005
 */
public abstract class BlockCipher
	{

// Hidden data members.

	/**
	 * Block length (bytes).
	 */
	private int myBlockLength;

	/**
	 * Key length (bytes).
	 */
	private int myKeyLength;

// Exported constructors.

	/**
	 * Construct a new block cipher object.
	 *
	 * @param  theBlockLength
	 *     Number of bytes in the plaintext or ciphertext block.
	 * @param  theKeyLength
	 *     Number of bytes in the encryption key.
	 *
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theBlockLength</TT> &lt;= 0 or
	 *     <TT>theKeyLength</TT> &lt;= 0.
	 */
	public BlockCipher
		(int theBlockLength,
		 int theKeyLength)
		{
		if (theBlockLength <= 0 || theKeyLength <= 0)
			{
			throw new IllegalArgumentException();
			}
		myBlockLength = theBlockLength;
		myKeyLength = theKeyLength;
		}

// Exported operations.

	/**
	 * Determine the block length of this block cipher.
	 *
	 * @return  Number of bytes in the plaintext or ciphertext block.
	 */
	public int getBlockLength()
		{
		return myBlockLength;
		}

	/**
	 * Determine the key length of this block cipher.
	 *
	 * @return  Number of bytes in the key.
	 */
	public int getKeyLength()
		{
		return myKeyLength;
		}

	/**
	 * Set the key to be used for all subsequent encryptions and decryptions.
	 * After this method returns, <TT>theKey</TT> is no longer needed and may be
	 * erased (set to all 0s). The length of <TT>theKey</TT> must be at least
	 * <TT>getKeyLength()</TT>. Only the first <TT>getKeyLength()</TT> bytes are
	 * used.
	 * <P>
	 * A subclass must override this method to implement the key setting
	 * algorithm.
	 *
	 * @param  theKey  Key (byte array).
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theKey</TT> is null. 
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theKey.length</TT> &lt;
	 *     <TT>getKeyLength()</TT>.
	 */
	public abstract void setKey
		(byte[] theKey);

	/**
	 * Encrypt the given plaintext block in place. On input, <TT>theBlock</TT>
	 * contains the plaintext. On output, the contents of <TT>theBlock</TT> have
	 * been replaced by the ciphertext. The length of <TT>theBlock</TT> must be
	 * at least <TT>getBlockLength()</TT>. Only the first
	 * <TT>getBlockLength()</TT> bytes are encrypted.
	 * <P>
	 * This is a convenience method that simply calls <TT>encrypt (theBlock,
	 * theBlock)</TT>.
	 *
	 * @param  theBlock  Block to be encrypted.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>theBlock</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>theBlock.length</TT> &lt;
	 *     <TT>getBlockLength()</TT>.
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if the key has been erased and not
	 *     re-set.
	 */
	public void encrypt
		(byte[] theBlock)
		{
		encrypt (theBlock, theBlock);
		}

	/**
	 * Encrypt the given plaintext block. On input, <TT>thePlaintext</TT>
	 * contains the plaintext. On output, the contents of <TT>theCiphertext</TT>
	 * have been replaced by the ciphertext. <TT>thePlaintext</TT> and
	 * <TT>theCiphertext</TT> may be the same block. The length of
	 * <TT>thePlaintext</TT> must be at least <TT>getBlockLength()</TT>. Only
	 * the first <TT>getBlockLength()</TT> bytes are read. The length of
	 * <TT>theCiphertext</TT> must be at least <TT>getBlockLength()</TT>. Only
	 * the first <TT>getBlockLength()</TT> bytes are written.
	 * <P>
	 * A subclass must override this method to implement the encryption
	 * algorithm.
	 *
	 * @param  thePlaintext   Input plaintext block to be encrypted.
	 * @param  theCiphertext  Output ciphertext block.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if <TT>thePlaintext</TT> is null or
	 *     <TT>theCiphertext</TT> is null.
	 * @exception  IllegalArgumentException
	 *     (unchecked exception) Thrown if <TT>thePlaintext.length</TT> &lt;
	 *     <TT>getBlockLength()</TT> or <TT>theCiphertext.length</TT> &lt;
	 *     <TT>getBlockLength()</TT>.
	 * @exception  IllegalStateException
	 *     (unchecked exception) Thrown if the key has been erased and not
	 *     re-set.
	 */
	public abstract void encrypt
		(byte[] thePlaintext,
		 byte[] theCiphertext);

	/**
	 * Erase this block cipher object's key material.
	 * <P>
	 * A subclass must override this method to set all of this block cipher
	 * object's key material to innocuous values (like 0). It is not acceptable
	 * merely to null out references to the key material, as that leaves the key
	 * material still stored in memory.
	 */
	public abstract void erase();

	}
