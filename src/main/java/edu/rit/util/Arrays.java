//******************************************************************************
//
// File:    Arrays.java
// Package: edu.rit.util
// Unit:    Class edu.rit.util.Arrays
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

package edu.rit.util;

import java.lang.reflect.Array;

/**
 * Class Arrays provides static methods for various operations on arrays and
 * matrices of primitive types and object types.
 * <P>
 * <I>Note:</I> The operations in class Arrays are not multiple thread safe.
 *
 * @author  Alan Kaminsky
 * @version 19-Nov-2007
 */
public class Arrays
	{

// Prevent construction.

	private Arrays()
		{
		}

// Exported array allocation operations.

	/**
	 * Allocate all elements in the given object array. Each array element
	 * (object) is created using the given class's no-argument constructor.
	 *
	 * @param  <T>    Array element data type.
	 * @param  <ST>   Data type of the objects to be created.
	 * @param  array  Array.
	 * @param  type   Class for the objects to be created.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  InstantiationException
	 *     Thrown if an instance of the given class cannot be instantiated.
	 * @exception  IllegalAccessException
	 *     Thrown if the given class or its no-argument constructor is not
	 *     accessible.
	 */
	public static <T, ST extends T> void allocate
		(T[] array,
		 Class<ST> type)
		throws InstantiationException, IllegalAccessException
		{
		allocate (array, new Range (0, array.length-1), type);
		}

	/**
	 * Allocate the elements within the given index range in the given object
	 * array. Each array element (object) is created using the given class's
	 * no-argument constructor.
	 *
	 * @param  <T>    Array element data type.
	 * @param  <ST>   Data type of the objects to be created.
	 * @param  array  Array.
	 * @param  range  Range of indexes to allocate.
	 * @param  type   Class for the objects to be created.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>range</TT> is
	 *     outside the bounds of the <TT>array</TT>.
	 * @exception  InstantiationException
	 *     Thrown if an instance of the given class cannot be instantiated.
	 * @exception  IllegalAccessException
	 *     Thrown if the given class or its no-argument constructor is not
	 *     accessible.
	 */
	public static <T, ST extends T> void allocate
		(T[] array,
		 Range range,
		 Class<ST> type)
		throws InstantiationException, IllegalAccessException
		{
		int lb = range.lb();
		int ub = range.ub();
		int stride = range.stride();
		if (0 > lb || ub >= array.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): array indexes = 0.." + (array.length-1) +
				 ", range = " + range);
			}
		for (int i = lb; i <= ub; i += stride)
			{
			array[i] = type.newInstance();
			}
		}

// Exported matrix allocation operations.

	/**
	 * Allocate the elements in the given Boolean matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is allocated with the
	 * given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(boolean[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given
	 * Boolean matrix. Each matrix row in the given row index range is allocated
	 * with the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(boolean[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new boolean [ncols];
			}
		}

	/**
	 * Allocate the elements in the given byte matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is allocated with the
	 * given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(byte[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given byte
	 * matrix. Each matrix row in the given row index range is allocated with
	 * the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(byte[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new byte [ncols];
			}
		}

	/**
	 * Allocate the elements in the given character matrix. Each matrix row in
	 * the full row index range (0..<TT>matrix.length</TT>-1) is allocated with
	 * the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(char[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given
	 * character matrix. Each matrix row in the given row index range is
	 * allocated with the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(char[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new char [ncols];
			}
		}

	/**
	 * Allocate the elements in the given double matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is allocated with the
	 * given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(double[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given
	 * double matrix. Each matrix row in the given row index range is allocated
	 * with the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(double[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new double [ncols];
			}
		}

	/**
	 * Allocate the elements in the given float matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is allocated with the
	 * given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(float[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given
	 * float matrix. Each matrix row in the given row index range is allocated
	 * with the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(float[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new float [ncols];
			}
		}

	/**
	 * Allocate the elements in the given integer matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is allocated with the
	 * given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(int[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given
	 * integer matrix. Each matrix row in the given row index range is allocated
	 * with the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(int[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new int [ncols];
			}
		}

	/**
	 * Allocate the elements in the given long matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is allocated with the
	 * given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(long[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given long
	 * matrix. Each matrix row in the given row index range is allocated with
	 * the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(long[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new long [ncols];
			}
		}

	/**
	 * Allocate the elements in the given short matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is allocated with the
	 * given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(short[][] matrix,
		 int ncols)
		{
		allocate (matrix, new Range (0, matrix.length-1), ncols);
		}

	/**
	 * Allocate the elements within the given row index range in the given short
	 * matrix. Each matrix row in the given row index range is allocated with
	 * the given number of columns.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 */
	public static void allocate
		(short[][] matrix,
		 Range rowRange,
		 int ncols)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = new short [ncols];
			}
		}

	/**
	 * Allocate the elements in the given object matrix. Each matrix row in the
	 * full row range (0..<TT>matrix.length</TT>-1) is allocated with the given
	 * number of columns. Within each row, the full column range
	 * (0..<TT>ncols</TT>-1) of elements is allocated. Each matrix element
	 * (object) is created using the given class's no-argument constructor.
	 *
	 * @param  <T>       Matrix element data type.
	 * @param  <ST>      Data type of the objects to be created.
	 * @param  matrix    Matrix.
	 * @param  ncols     Number of columns in each row.
	 * @param  type      Class for the objects to be created.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 * @exception  InstantiationException
	 *     Thrown if an instance of the given class cannot be instantiated.
	 * @exception  IllegalAccessException
	 *     Thrown if the given class or its no-argument constructor is not
	 *     accessible.
	 */
	public static <T, ST extends T> void allocate
		(T[][] matrix,
		 int ncols,
		 Class<ST> type)
		throws InstantiationException, IllegalAccessException
		{
		allocate
			(matrix,
			 new Range (0, matrix.length-1),
			 ncols,
			 new Range (0, ncols-1),
			 type);
		}

	/**
	 * Allocate the elements within the given row index range in the given
	 * object matrix. Each matrix row in the given row index range is allocated
	 * with the given number of columns. Within each row, the full column range
	 * (0..<TT>ncols</TT>-1) of elements is allocated. Each matrix element
	 * (object) is created using the given class's no-argument constructor.
	 *
	 * @param  <T>       Matrix element data type.
	 * @param  <ST>      Data type of the objects to be created.
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 * @param  type      Class for the objects to be created.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 * @exception  InstantiationException
	 *     Thrown if an instance of the given class cannot be instantiated.
	 * @exception  IllegalAccessException
	 *     Thrown if the given class or its no-argument constructor is not
	 *     accessible.
	 */
	public static <T, ST extends T> void allocate
		(T[][] matrix,
		 Range rowRange,
		 int ncols,
		 Class<ST> type)
		throws InstantiationException, IllegalAccessException
		{
		allocate
			(matrix,
			 rowRange,
			 ncols,
			 new Range (0, ncols-1),
			 type);
		}

	/**
	 * Allocate the elements within the given row and column index ranges in the
	 * given object matrix. Each matrix row is allocated with the given number
	 * of columns. Each matrix element (object) is created using the given
	 * class's no-argument constructor.
	 *
	 * @param  <T>       Matrix element data type.
	 * @param  <ST>      Data type of the objects to be created.
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to allocate.
	 * @param  ncols     Number of columns in each row.
	 * @param  colRange  Range of column indexes to allocate.
	 * @param  type      Class for the objects to be created.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>. Thrown if any index in
	 *     the <TT>colRange</TT> is outside the bounds 0 .. <TT>ncols</TT>-1.
	 * @exception  NegativeArraySizeException
	 *     (unchecked exception) Thrown if <TT>ncols</TT> &lt; 0.
	 * @exception  InstantiationException
	 *     Thrown if an instance of the given class cannot be instantiated.
	 * @exception  IllegalAccessException
	 *     Thrown if the given class or its no-argument constructor is not
	 *     accessible.
	 */
	public static <T, ST extends T> void allocate
		(T[][] matrix,
		 Range rowRange,
		 int ncols,
		 Range colRange,
		 Class<ST> type)
		throws InstantiationException, IllegalAccessException
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		int colLb = colRange.lb();
		int colUb = colRange.ub();
		int colStride = colRange.stride();
		if (0 > colLb || colUb >= ncols)
			{
			throw new IndexOutOfBoundsException
				("Arrays.allocate(): matrix column indexes = 0.." + (ncols-1) +
				 ", column range = " + colRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			T[] matrix_i = (T[]) Array.newInstance (type, ncols);
			matrix[i] = matrix_i;
			for (int j = colLb; j <= colUb; j += colStride)
				{
				matrix_i[j] = type.newInstance();
				}
			}
		}

// Exported array deallocation operations.

	/**
	 * Deallocate all elements in the given object array. Each array element
	 * is set to null.
	 *
	 * @param  <T>    Array element data type.
	 * @param  array  Array.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static <T> void deallocate
		(T[] array)
		{
		deallocate (array, new Range (0, array.length-1));
		}

	/**
	 * Deallocate the elements within the given index range in the given object
	 * array. Each array element is set to null.
	 *
	 * @param  <T>    Array element data type.
	 * @param  array  Array.
	 * @param  range  Range of indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>range</TT> is
	 *     outside the bounds of the <TT>array</TT>.
	 */
	public static <T> void deallocate
		(T[] array,
		 Range range)
		{
		int lb = range.lb();
		int ub = range.ub();
		int stride = range.stride();
		if (0 > lb || ub >= array.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): array indexes = 0.." + (array.length-1) +
				 ", range = " + range);
			}
		for (int i = lb; i <= ub; i += stride)
			{
			array[i] = null;
			}
		}

// Exported matrix deallocation operations.

	/**
	 * Deallocate the elements in the given Boolean matrix. Each matrix row in
	 * the full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(boolean[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * Boolean matrix. Each matrix row in the given row index range is set to
	 * null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(boolean[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given byte matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(byte[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * byte matrix. Each matrix row in the given row index range is set to null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(byte[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given character matrix. Each matrix row in
	 * the full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(char[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * character matrix. Each matrix row in the given row index range is set to
	 * null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(char[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given double matrix. Each matrix row in
	 * the full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(double[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * double matrix. Each matrix row in the given row index range is set to
	 * null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(double[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given float matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(float[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * float matrix. Each matrix row in the given row index range is set to
	 * null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(float[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given integer matrix. Each matrix row in
	 * the full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(int[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * integer matrix. Each matrix row in the given row index range is set to
	 * null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(int[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given long matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(long[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * long matrix. Each matrix row in the given row index range is set to null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(long[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given short matrix. Each matrix row in the
	 * full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static void deallocate
		(short[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * short matrix. Each matrix row in the given row index range is set to
	 * null.
	 *
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static void deallocate
		(short[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements in the given object matrix. Each matrix row in
	 * the full row index range (0..<TT>matrix.length</TT>-1) is set to null.
	 *
	 * @param  <T>       Matrix element data type.
	 * @param  matrix    Matrix.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 */
	public static <T> void deallocate
		(T[][] matrix)
		{
		deallocate (matrix, new Range (0, matrix.length-1));
		}

	/**
	 * Deallocate the elements within the given row index range in the given
	 * object matrix. Each matrix row in the given row index range is set to
	 * null.
	 *
	 * @param  <T>       Matrix element data type.
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>.
	 */
	public static <T> void deallocate
		(T[][] matrix,
		 Range rowRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			matrix[i] = null;
			}
		}

	/**
	 * Deallocate the elements within the given row and column index ranges in
	 * the given object matrix. Each matrix element within the given row and
	 * column index ranges is set to null.
	 *
	 * @param  <T>       Matrix element data type.
	 * @param  matrix    Matrix.
	 * @param  rowRange  Range of row indexes to deallocate.
	 * @param  colRange  Range of column indexes to deallocate.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in the <TT>rowRange</TT> is
	 *     outside the row bounds of the <TT>matrix</TT>. Thrown if any index in
	 *     the <TT>colRange</TT> is outside the bounds 0 .. <TT>ncols</TT>-1.
	 */
	public static <T> void deallocate
		(T[][] matrix,
		 Range rowRange,
		 Range colRange)
		{
		int rowLb = rowRange.lb();
		int rowUb = rowRange.ub();
		int rowStride = rowRange.stride();
		if (0 > rowLb || rowUb >= matrix.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.deallocate(): matrix row indexes = 0.." +
				 (matrix.length-1) + ", row range = " + rowRange);
			}
		int colLb = colRange.lb();
		int colUb = colRange.ub();
		int colStride = colRange.stride();
		for (int i = rowLb; i <= rowUb; i += rowStride)
			{
			T[] matrix_i = matrix[i];
			int ncols = matrix_i.length;
			if (0 > colLb || colUb >= ncols)
				{
				throw new IndexOutOfBoundsException
					("Arrays.allocate(): matrix column indexes = 0.." +
					 (ncols-1) + ", column range = " + colRange);
				}
			for (int j = colLb; j <= colUb; j += colStride)
				{
				matrix_i[j] = null;
				}
			}
		}

// Exported array copying operations.

	/**
	 * Copy a range of elements from one Boolean array to a range of elements in
	 * another Boolean array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(boolean[] src,
		 Range srcRange,
		 boolean[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one byte array to a range of elements in
	 * another byte array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(byte[] src,
		 Range srcRange,
		 byte[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one character array to a range of elements
	 * in another character array. The number of elements copied is the smaller
	 * of <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or
	 * both of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be
	 * greater than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(char[] src,
		 Range srcRange,
		 char[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one double array to a range of elements in
	 * another double array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(double[] src,
		 Range srcRange,
		 double[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one float array to a range of elements in
	 * another float array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(float[] src,
		 Range srcRange,
		 float[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one integer array to a range of elements in
	 * another integer array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(int[] src,
		 Range srcRange,
		 int[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one long array to a range of elements in
	 * another long array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(long[] src,
		 Range srcRange,
		 long[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one short array to a range of elements in
	 * another short array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 *
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static void copy
		(short[] src,
		 Range srcRange,
		 short[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

	/**
	 * Copy a range of elements from one object array to a range of elements in
	 * another object array. The number of elements copied is the smaller of
	 * <TT>srcRange</TT>'s length and <TT>dstRange</TT>'s length. Either or both
	 * of <TT>srcRange</TT>'s and <TT>dstRange</TT>'s strides may be greater
	 * than 1.
	 * <P>
	 * <I>Note:</I> This method does a <I>shallow copy.</I> Only the references
	 * to the array elements are copied from the source array to the destination
	 * array.
	 *
	 * @param  <ST>      Source array element data type.
	 * @param  <DT>      Destination array element data type.
	 * @param  src       Source array.
	 * @param  srcRange  Range of source elements.
	 * @param  dst       Destination array.
	 * @param  dstRange  Range of destination elements.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRange</TT> is
	 *     outside the bounds of the source array. Thrown if any index in
	 *     <TT>dstRange</TT> is outside the bounds of the destination array.
	 */
	public static <DT, ST extends DT> void copy
		(ST[] src,
		 Range srcRange,
		 DT[] dst,
		 Range dstRange)
		{
		int len = Math.min (srcRange.length(), dstRange.length());
		if (len == 0) return;
		int srcLower = srcRange.lb();
		int dstLower = dstRange.lb();
		if (src == dst && srcLower == dstLower) return;
		int srcStride = srcRange.stride();
		int dstStride = dstRange.stride();
		if (srcStride == 1 && dstStride == 1)
			{
			System.arraycopy (src, srcLower, dst, dstLower, len);
			}
		else
			{
			int srcUpper = srcRange.ub();
			int dstUpper = dstRange.ub();
			if (0 > srcLower || srcUpper >= src.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): src indexes = 0.." + (src.length-1) +
					 ", srcRange = " + srcRange);
				}
			if (0 > dstLower || dstUpper >= dst.length)
				{
				throw new IndexOutOfBoundsException
					("Arrays.copy(): dst indexes = 0.." + (dst.length-1) +
					 ", dstRange = " + dstRange);
				}
			if (src != dst || srcLower > dstLower)
				{
				for
					(int i = srcLower, j = dstLower;
					 i <= srcUpper;
					 i += srcStride, j += dstStride)
					{
					dst[j] = src[i];
					}
				}
			else if (srcLower < dstLower)
				{
				for
					(int i = srcUpper, j = dstUpper;
					 i >= srcLower;
					 i -= srcStride, j -= dstStride)
					{
					dst[j] = src[i];
					}
				}
			}
		}

// Exported matrix copying operations.

	/**
	 * Copy a range of elements from one Boolean matrix to a range of elements
	 * in another Boolean matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(boolean[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 boolean[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one byte matrix to a range of elements
	 * in another byte matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(byte[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 byte[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one character matrix to a range of elements
	 * in another character matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(char[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 char[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one double matrix to a range of elements
	 * in another double matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(double[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 double[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one float matrix to a range of elements
	 * in another float matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(float[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 float[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one integer matrix to a range of elements
	 * in another integer matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(int[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 int[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one long matrix to a range of elements
	 * in another long matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(long[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 long[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one short matrix to a range of elements
	 * in another short matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 *
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static void copy
		(short[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 short[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

	/**
	 * Copy a range of elements from one object matrix to a range of elements
	 * in another object matrix. The number of rows copied is the smaller of
	 * <TT>srcRowRange</TT>'s length and <TT>dstRowRange</TT>'s length. Within
	 * each row, the number of columns copied is the smaller of
	 * <TT>srcColRange</TT>'s length and <TT>dstColRange</TT>'s length. Any of
	 * <TT>srcRowRange</TT>'s, <TT>srcColRange</TT>'s, <TT>dstRowRange</TT>'s,
	 * and <TT>dstColRange</TT>'s strides may be greater than 1. It is assumed
	 * that the source matrix is fully allocated; each row in the source matrix
	 * is the same length; the destination matrix is fully allocated; and each
	 * row in the destination matrix is the same length.
	 * <P>
	 * <I>Note:</I> This method does a <I>shallow copy.</I> Only the references
	 * to the matrix elements are copied from the source matrix to the
	 * destination matrix.
	 *
	 * @param  <ST>         Source matrix element data type.
	 * @param  <DT>         Destination matrix element data type.
	 * @param  src          Source matrix.
	 * @param  srcRowRange  Range of source rows.
	 * @param  srcColRange  Range of source columns.
	 * @param  dst          Destination matrix.
	 * @param  dstRowRange  Range of destination rows.
	 * @param  dstColRange  Range of destination columns.
	 *
	 * @exception  NullPointerException
	 *     (unchecked exception) Thrown if any argument is null.
	 * @exception  IndexOutOfBoundsException
	 *     (unchecked exception) Thrown if any index in <TT>srcRowRange</TT> is
	 *     outside the row bounds of the source matrix. Thrown if any index in
	 *     <TT>srcColRange</TT> is outside the column bounds of the source
	 *     matrix. Thrown if any index in <TT>dstRowRange</TT> is outside the
	 *     row bounds of the destination matrix. Thrown if any index in
	 *     <TT>dstColRange</TT> is outside the column bounds of the destination
	 *     matrix.
	 */
	public static <DT, ST extends DT> void copy
		(ST[][] src,
		 Range srcRowRange,
		 Range srcColRange,
		 DT[][] dst,
		 Range dstRowRange,
		 Range dstColRange)
		{
		int len = Math.min (srcRowRange.length(), dstRowRange.length());
		if (len == 0) return;
		int srcRowLower = srcRowRange.lb();
		int dstRowLower = dstRowRange.lb();
		int srcRowStride = srcRowRange.stride();
		int dstRowStride = dstRowRange.stride();
		int srcRowUpper = srcRowRange.ub();
		int dstRowUpper = dstRowRange.ub();
		if (0 > srcRowLower || srcRowUpper >= src.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): src row indexes = 0.." + (src.length-1) +
				 ", srcRowRange = " + srcRowRange);
			}
		if (0 > dstRowLower || dstRowUpper >= dst.length)
			{
			throw new IndexOutOfBoundsException
				("Arrays.copy(): dst row indexes = 0.." + (dst.length-1) +
				 ", dstRowRange = " + dstRowRange);
			}
		if (src != dst || srcRowLower > dstRowLower)
			{
			for
				(int i = srcRowLower, j = dstRowLower;
				 i <= srcRowUpper;
				 i += srcRowStride, j += dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		else if (srcRowLower < dstRowLower)
			{
			for
				(int i = srcRowUpper, j = dstRowUpper;
				 i >= srcRowLower;
				 i -= srcRowStride, j -= dstRowStride)
				{
				copy (src[i], srcColRange, dst[j], dstColRange);
				}
			}
		}

// Exported array length operations.

	/**
	 * Determine the number of elements in the given Boolean array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(boolean[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given byte array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(byte[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given character array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(char[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given double array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(double[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given float array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(float[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given integer array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(int[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given long array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(long[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given short array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static int length
		(short[] array)
		{
		return array == null ? 0 : array.length;
		}

	/**
	 * Determine the number of elements in the given object array. If
	 * <TT>array</TT> is null, 0 is returned.
	 *
	 * @param  <T>    Array element data type.
	 * @param  array  Array.
	 *
	 * @return  Number of elements in <TT>array</TT>.
	 */
	public static <T> int length
		(T[] array)
		{
		return array == null ? 0 : array.length;
		}

// Exported matrix length operations.

	/**
	 * Determine the number of rows in the given Boolean matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(boolean[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given byte matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(byte[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given character matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(char[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given double matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(double[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given float matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(float[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given integer matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(int[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given long matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(long[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given short matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static int rowLength
		(short[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of rows in the given object matrix. If
	 * <TT>matrix</TT> is null, 0 is returned.
	 *
	 * @param  <T>     Matrix element data type.
	 * @param  matrix  Matrix.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 */
	public static <T> int rowLength
		(T[][] matrix)
		{
		return matrix == null ? 0 : matrix.length;
		}

	/**
	 * Determine the number of columns in the given row of the given Boolean
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(boolean[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given byte
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(byte[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given character
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(char[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given double
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(double[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given float
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(float[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given integer
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(int[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given long
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(long[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given short
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static int colLength
		(short[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	/**
	 * Determine the number of columns in the given row of the given object
	 * matrix. If <TT>matrix</TT> is null or the given row is not allocated, 0
	 * is returned.
	 *
	 * @param  <T>     Matrix element data type.
	 * @param  matrix  Matrix.
	 * @param  i       Row index.
	 *
	 * @return  Number of rows in <TT>matrix</TT>.
	 *
	 * @exception  ArrayIndexOutOfBoundsException
	 *     (unchecked exception) Thrown if <TT>i</TT> is out of bounds.
	 */
	public static <T> int colLength
		(T[][] matrix,
		 int i)
		{
		return matrix == null ? 0 : length (matrix[i]);
		}

	}
