import edu.gwu.lintool.ComplexNumber;
import edu.gwu.lintool.LinResult;
import edu.gwu.lintool.LinToolEmpty;

/**
 * Rohan Sunder
 * CSCI 6342 - Linear Algebra: A Computational Approach
 * @author rsunder
 *
 */
public class RohanA3_LinTool extends LinToolEmpty {
    
	//Vector Operations
	// Assignment 1 Methods
	///////////////////////////////////////////////
	
	public double[] add (double[] u, double[] v) {
		
		//Null check
		if (u == null || v == null) {
			System.out.println("Vectors are null. Aborting.");
			return null;
		}
    	
		//Size Check
		if (u.length != v.length) {
			System.out.println("Vector parameters have unequal sizes. Not adding together.");
			return null;
		}
		
		double[] sumVector = new double[u.length];
		
		for (int i = 0; i < u.length; i++) {
			sumVector[i] = u[i] + v[i];
		}
		
		return sumVector;
    }
	
    public double norm (double[] v) {
    	
    	//Null Check
    	if (v == null) {
    		System.out.println("Vector is null. Aborting");
    		return -1;
    	}
    	
    	double val = 0.0;
    	
    	for (int i = 0; i < v.length; i++) {
    		val += v[i] * v[i];
    	}
    	
    	 return Math.sqrt(val);
    }
    
    public double dotProduct (double[] u, double[] v) {
    	
		// Null check
    	if (u == null || v == null) {
			System.out.println("Vectors are null. Aborting.");
			return -1;
		}
		
    	//Size check
		if (u.length != v.length) {
			System.out.println("Vector parameters have unequal sizes. Not adding together.");
			return -1;
		}
    	
    	double val = 0.0;
    	
    	for (int i = 0; i < u.length; i++) {
    		val += u[i] * v[i];
    	}
    	
    	return val;
    }
    
    public double[] scalarProduct (double alpha, double[] v) {
    	
    	//Null check
    	if (v == null) {
    		System.out.println("Vector is null. Aborting");
    		return null;
    	}
    	
    	double ret[] = new double[v.length];
    	
    	for (int i = 0; i < v.length; i++) {
    		ret[i] = v[i] * alpha;
    	}
    	
    	return ret;
    	
    }
    public boolean approxEquals (double[] u, double[] v, double errorTolerance) {
    	
    	if (u == null || v == null) {
    		System.out.println("Vector(s) are null. Aborting.");
    		return false;
    	}
    	
    	//Length check
    	if (u.length != v.length) {
    		return false;
    	}
    	
    	for (int i = 0; i < u.length; i++) {
    		if (Math.abs(u[i] - v[i]) > errorTolerance) {
    			return false;
    		}
    	}
    	
    	return true;
    }
    
    public double cosine (double[] u, double[] v) {
    	
    	if (u == null || v == null) {
    		System.out.println("Vector(s) are null. Aborting.");
    		return -1;
    	}
    	
    	RohanA3_LinTool tool = new RohanA3_LinTool();
    	
    	return (tool.dotProduct(u, v) / (tool.norm(v)) * (tool.norm(u)));
    	
    }
    
    // Matrix operations:
    
    public double[][] add (double[][] A, double[][] B) {
    	
    	if (A == null || B == null) {
    		System.out.println("Matrices are null. Aborting");
    		return null;
    	}
    	
    	if (A.length != B.length || A[0].length != B[0].length) {
    		System.out.println("Matrices are of differing dimensions and cannot be added.");
    		return null;
    	}
    	
    	double retMatrix[][] = new double[A.length][A[0].length];
    	
    	for (int i = 0; i < A.length; i++) {
    		for (int j = 0; j < A[i].length; j++) {
    			retMatrix[i][j] = A[i][j] + B[i][j];
    		}
    	}
    	
    	return retMatrix;
    }
    
    public double[][] scalarProduct (double alpha, double[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null.");
    		return null;
    	}
    	
    	double retMatrix[][] = new double[A.length][A[0].length];
    
    	for (int i = 0; i < A.length; i++) {
    		for (int j = 0; j < A[i].length; j++) {
    			retMatrix[i][j] = A[i][j] * alpha;
    		}
    	}
    	
    	return retMatrix;
    	    	
    }
    
    public double[][] mult (double[][] A, double[][] B) {
    	
		if (A == null || B == null) {
			System.out.println("Matrices are empty");
			return null;
		}
		
		int Arows = A.length;
		int Acolumns = A[0].length; 
		
		int Brows = B.length;
		int Bcolumns = B[0].length;
		
		if (Acolumns != Brows) {
			System.out.println("Invalid matrix dimensions. Aborting program.");
			return null;
		}
		
		double[][] ret = new double[Arows][Bcolumns];
		double tempValue = 0.0;
		
		for (int i = 0; i < Arows; i++) {
			for (int j = 0; j < Bcolumns; j++) {
				for (int k = 0; k < Brows; k++) {
					tempValue += A[i][k] * B[k][j];
				}
				ret[i][j] = tempValue;
				tempValue = 0.0;
			}
		}

		return ret;
    	
    }
    public double[] matrixVectorMult (double[][] A, double[] v) {
    	
    	if (A == null || v == null) {
    		System.out.println("Matrix or vector is empty. Aborting.");
    		return null;
    	}
    	
    	if (v.length != A[0].length) {
    		return null;
    	}
    	
		double[] ret = new double[A.length];
		double tempValue = 0.0;
		
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[i].length; j++) {
				tempValue += A[i][j] * v[j];
			}
			ret[i] = tempValue;
			tempValue = 0.0;
		}
		
		return ret;
    	
    }
    public double[] vectorMatrixMult (double[] v, double[][] A) {
    	
    	if (A == null || v == null) {
    		System.out.println("Matrix or vector is empty. Aborting.");
    		return null;
    	}
    	
    	int n = A[0].length;
    	int m = A.length;
    	
    	if (v.length != m) {
    		
    		return null;
    	}
    	
		double[] ret = new double[n];
		double tempValue = 0.0;
		
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				tempValue += A[j][i] * v[j];
			}
			ret[i] = tempValue;
			tempValue = 0.0;
		}
		
		return ret;
    }
    
    public double[][] vectorLeftMult (double[] u, double[] v) {
    
    	if (u == null || v == null) {
    		System.out.println("Vectors are null.");
    		return null;
    	}
    	
    	if (u.length != v.length) {
    		System.out.println("Vectors are of different sizes.");
    		return null;
    	}
    	
    	double ret[][] = new double[v.length][v.length];
    	for (int i = 0; i < ret.length; i++) {
    		for (int j = 0; j < ret[i].length; j++) {
    			ret[i][j] = u[i] * v[j];
    		}
    	}
    	
    	return ret;
    }
    
    public double[][] transpose (double[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null. Aborting.");
    		return null;
    	}
    	
    	double[][] B = new double[A[0].length][A.length];
    	
    	for (int i = 0; i < B.length; i++) {
    		for (int j = 0; j < B[i].length; j++) {
    			B[i][j] = A[j][i];
    		}
    	}
    	
    	return B;
    }
    
    public double frobeniusNorm (double[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null. Aborting.");
    		return -1;
    	}
    	
    	double sum = 0.0;
    	
    	for (int i = 0; i < A.length; i++) {
    		for (int j = 0; j < A[0].length; j++) {
    			sum = sum + (A[i][j] * A[i][j]);
    		}
    	}
    	
    	return Math.sqrt(sum);
    }
    
    public boolean approxEquals (double[][] A, double[][] B, double errorTolerance) {
    	
    	if (A.length != B.length || A[0].length != B[0].length) {
    		System.out.println("Matrices are of differing dimensions and cannot be equal.");
    		return false;
    	}
    	
    	for (int i = 0; i < A.length; i++) {
    		for (int j = 0; j < A[i].length; j++) {
    			if (Math.abs(A[i][j] - B[i][j]) > errorTolerance) {
        			return false;
        		}
    		}
    	}
    	
    	return true;
    }
    
    public double[] getColumnAsVector (double[][] A, int col) {
    	
    	double[] vectorColumn = new double[A[0].length];
    	
    	if (col >= A.length || col < 0) {
    		System.out.println("Col number is out of the matrix");
    		return null;
    	}
    	
    	for (int i = 0; i < A[col].length; i++) {
    		vectorColumn[i] = A[col][i];
    	}
    	
    	return vectorColumn;
    }
    
    public double[] getRowAsVector (double[][] A, int row) {
    	
    	if (row >= A.length || row < 0) {
    		System.out.println("Row number is out of the matrix");
    		return null;
    	}
    	
    	return A[row];
    }
    
    // Complex vector operations:
    public ComplexNumber[] add (ComplexNumber[] u, ComplexNumber[] v) {
    	
    	if (u == null || v == null) {
    		System.out.println("Complex vectors are null.");
    		return null;
    	}
    	
    	if (u.length != v.length) {
    		System.out.println("Complex vectors are not of the same size.");
    		return null;
    	}
    	
    	ComplexNumber[] ret = new ComplexNumber[u.length];
    	
    	for (int i = 0; i < u.length; i++) {
    		ret[i] = u[i].add(v[i]);
    	}
    	
    	return ret;
    }
    
    public double norm (ComplexNumber[] v) {
    	
    	if (v == null) {
    		System.out.println("Complex vector is null. Aborting.");
    		return -1;
    	}
    	
    	double norm = 0.0;
    	
    	for (int i = 0; i < v.length; i++) {
    		norm += v[i].pow(2).magnitude();
    	}
    	
    	return Math.sqrt(norm);
    }
    
    public ComplexNumber[] scalarProduct (ComplexNumber alpha, ComplexNumber[] v) {
    	
    	if (alpha == null || v == null) {
    		System.out.println("Alpha or the vector is null. Aborting.");
    		return null;
    	}
    	
    	ComplexNumber[] ret = new ComplexNumber[v.length];
    	
    	for (int i = 0; i < v.length; i++) {
    		ret[i] = v[i].mult(alpha);
    	}
    	
    	return ret;
    }
    
    public ComplexNumber dotProduct (ComplexNumber[] u, ComplexNumber[] v) {
    	
    	if (u == null || v == null) {
    		System.out.println("Complex vectors are null.");
    		return null;
    	}
    	
    	if (u.length != v.length) {
    		System.out.println("Complex vectors are not of the same size.");
    		return null;
    	}
    	
    	double re = 0.0;
    	double im = 0.0;
    	
    	for (int i = 0; i < u.length; i++) {
    		re += (u[i].re * v[i].re) + (u[i].im * v[i].im);
    		im += (v[i].re * u[i].im) - (u[i].re * v[i].im);
    	}
    	
    	return ComplexNumber.makeComplexNumber(re, im);
    	
    }
    // Note: dot product is the Hermitian dot product:
    
    
    // Assignment 2 Methods
    /////////////////////////////////////////////////////////////////////////////////
    // Solving equations:
    public LinResult computeREF (double[][] A, double[] b) {
    	
    	if (A == null || b == null) {
    		System.out.println("Matrix or vector is empty. Aborting.");
    		return null;
    	}
    	
    	if (b.length != A.length) {
    		System.out.println("Matrix rows do no match vector rows.");
    		return null;
    	}
    	
    	System.out.println("Beginning matrix:");
    	this.printMatrix(A);
    	
    	double[][] Aplus = this.computeAugmentedMatrix(A, b);
    	
    	LinResult lin = this.privateComputeRREF(Aplus);

    	return lin;
    }
    
    public LinResult computeRREF (double[][] A, double[] b) {
    	
    	if (A == null || b == null) {
    		System.out.println("Matrix or vector is empty. Aborting.");
    		return null;
    	}
    	
    	if (b.length != A.length) {
    		System.out.println("Matrix rows do no match vector rows.");
    		return null;
    	}
    	
    	double[][] Aplus = this.computeAugmentedMatrix(A, b);
    	
    	LinResult lin = this.privateComputeRREF(Aplus);
    	
    	return lin;
    	
    }
    
    public LinResult solveFromREF (double[][] A, double[] b) {
    	
    	if (A == null || b == null) {
    		System.out.println("Matrix or vector is empty. Aborting.");
    		return null;
    	}
    	
    	if (b.length != A.length) {
    		System.out.println("Matrix rows do no match vector rows.");
    		return null;
    	}
    	
    	double[][] Aplus = this.computeAugmentedMatrix(A, b);
    	
    	LinResult lin = this.privateComputeRREF(Aplus);
    	
    	LinResult ret = lin;
//    	double[][] refM = ret.ref;
//    	
//    	int end = refM[0].length;
    	
//    	double[] x = new double[ret.rank];
//    	for (int i = b.length-1; i >= 0; i--) {
//    		double sum = 0.0;
//    		for (int j = i+1; j < b.length; j++) {
//    			sum += refM[i][j] * x[j];
//    		}
//    		x[i] = (refM[i][end-1] - sum) / refM[i][i];
//    	}
    	
    	//ret.x = x;
    	
    	return ret;
    	
    }
    
    public LinResult solveFromRREF (double[][] A, double[] b) {
    	
    	if (A == null || b == null) {
    		System.out.println("Matrix or vector is empty. Aborting.");
    		return null;
    	}
    	
    	if (b.length != A.length) {
    		System.out.println("Matrix rows do no match vector rows.");
    		return null;
    	}
    	
    	double[][] Aplus = this.computeAugmentedMatrix(A, b);
    	
    	LinResult lin = this.privateComputeRREF(Aplus);
    	
    	return lin;
    }
    
    public LinResult inverse (double[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is empty. Aborting.");
    		return null;
    	}
    	
    	double[] b = new double[A.length];
    	
    	for (int i = 0; i < b.length; i++) {
    		b[i] = 0.0;
    	}
    	
    	double[][] Aplus = this.computeAugmentedMatrix(A, b);
    	
    	LinResult lin = this.privateComputeRREF(Aplus);
    	
    	return lin;
    }
    
    /**
     * Given the original matrix A and the result vector b, produce an augmented matrix that
     * combines the vector with 
     * @param A
     * @param b
     * @return
     */
    private double[][] computeAugmentedMatrix(double[][] A, double[] b) {
    	
    	double[][] ret = new double[A.length][A[0].length + 1];
    	
    	for (int i = 0; i < A.length; i++) {
    		
    		for (int j = 0; j < A[i].length; j++) {
    			ret[i][j] = A[i][j];
    		}
    		ret[i][A[i].length] = b[i];
    	}
    	
    	return ret;
    }
    
    /**
     * Finds the next pivot given the currentRow, currentColumn and the augmentedMatrix.
     * @param row
     * @param col
     * @param Aplus
     * @return
     */
    private int[] findNextPivot (int row, int col, double[][] Aplus) {
    	
    	int retIndexes[] = new int[2];
    	
    	if (row >= Aplus.length || col >= Aplus[0].length-1) {
    		return null;
    	}
    	
    	if (Aplus[row][col] != 0) {
    		
    		retIndexes[0] = row;
    		retIndexes[1] = col;
    		return retIndexes;
    		
    	} else {
    		
    		for (int j = col; j < Aplus[row].length-1; j++) {
        		
    			for (int i = row; i < Aplus.length; i++) {
        			
        			if (Aplus[i][j] != 0) {
        				
        				retIndexes[0] = i;
        				retIndexes[1] = j;
        				
        				System.out.println("retIndexes row and column: " + retIndexes[0] + " " + retIndexes[1]);
        				return retIndexes;
        				
        			}
        		}
    		}
    	}
    	
    	return null;
    }
    
    /**
     * Swaps the row that has the desired pivot with the original row in the augmented matrix.
     * @param swapRow
     * @param originalRow
     * @param aMatrix
     */
    private void swap (int swapRow, int originalRow, double[][] aMatrix) {
    	
    	double[] swappedRow = aMatrix[swapRow];
    	double[] origRow = aMatrix[originalRow];
    	
    	aMatrix[originalRow] = swappedRow;
    	aMatrix[swapRow] = origRow;
    	
    }
    
    /**
     * Private version of computing the RREF. Will be used in pretty much all of the main methods.
     * @param A
     * @param b
     * @return
     */
    private LinResult privateComputeRREF(double[][] Aplus) {
    	
    	LinResult lin = new LinResult();
    	int currentRow = 0;
    	int currentCol = 0;
    	
    	boolean inverseAvailable = false;

    	int numOfRows = Aplus.length;
    	int numOfCols = Aplus[0].length;
    	
    	if (numOfRows == numOfCols-1) {
    		inverseAvailable = true;
    	}
    	
    	//Identity matrix for inverse.
    	double[][] I = this.generateIdentityMatrix(numOfCols-1);
    	
    	lin.rank = 0;
    	lin.isPivotColumn = new boolean[numOfCols-1];
    	lin.pivotRow = new int[numOfRows];
    	
    	// Index 0 = Row. Index 0 = Column.
    	int pivotIndexes[] = new int[2];
    	
    	//Invokes finding the next Pivot.
    	pivotIndexes = this.findNextPivot(currentRow, currentCol, Aplus);
    	
    	System.out.println("Original Matrix:");
    	this.printMatrix(Aplus);
    	
    	// The Magic
    	while ((pivotIndexes != null) && (pivotIndexes[0] >= 0 && pivotIndexes[1] >= 0)) {
    		
    		//Swap into currentRow
    		if (pivotIndexes[0] > currentRow) {
    			this.swap (pivotIndexes[0], currentRow, Aplus);
    			System.out.println("Swapped Matrix:");
    			this.printMatrix(Aplus);
    		}
    		
    		currentCol = pivotIndexes[1];
    		
    		lin.isPivotColumn[currentCol] = true;
    		lin.rank = lin.rank + 1;
    		lin.pivotRow[lin.rank-1] = currentRow;

    		
    		//Make coefficient of the pivot equal 1.
    		double alpha = Aplus[currentRow][currentCol];
    		for (int j = 0; j < numOfCols; j++) {
    			Aplus[currentRow][j] = Aplus[currentRow][j] / alpha;
    			
    			if (inverseAvailable && j < numOfCols-1) {
    				I[currentRow][j] = I[currentRow][j] / alpha;
    			}
    		}
    		
    		// Zero out the coefficients below the pivot.
    		for (int r = currentRow+1; r < numOfRows; r++) {
    			
    			double beta = Aplus[r][currentCol] / Aplus[currentRow][currentCol];
    			
    			for (int j = 0; j < numOfCols; j++) {
    				Aplus[r][j] = Aplus[r][j] - (Aplus[currentRow][j] * beta);
    				
    				if (inverseAvailable && j < numOfCols-1) {
        				I[r][j] = I[r][j] - (I[currentRow][j] * beta);
    				}

    			}
    		}
    		
    		// Increment counters
    		currentRow = currentRow + 1;
    		currentCol = currentCol + 1;
    		pivotIndexes = this.findNextPivot(currentRow, currentCol, Aplus);
    	}
    	
    	System.out.println("REF Matrix:");
    	this.printMatrix(Aplus);
    	lin.ref = this.makeDeepCopy(Aplus);
    	lin.rref = null;
    	
    	int pivotIndex[];
    	
    	for (int k = 1; k <= lin.rank; k++) {
    		
    		pivotIndex = this.findPivot(lin, k, Aplus);
    		
    		int pRow = pivotIndex[0];
    		int pCol = pivotIndex[1];
    		
    		// Zero out the coefficients above the pivot.
    		for (int r = pRow-1; r >= 0; r--) {
    			
    			double beta = Aplus[r][pCol] / Aplus[pRow][pCol];
    			
    			for (int j = 0; j < numOfCols; j++) {
    				Aplus[r][j] = Aplus[r][j] - (Aplus[pRow][j] * beta);
    				
    				//Add row operations to inverse.
    				if (inverseAvailable && j < numOfCols-1) {
    					I[r][j] = I[r][j] - (I[pRow][j] * beta);	
    				}
    				
    			}
    		}
    	}
    	
    	System.out.println("RREF Matrix:");
    	lin.rref = Aplus;
    	this.printMatrix(Aplus);
    	
    	if (this.existsContradictionRow(Aplus, lin)) {
    		
    		lin.solutionExists = false;
    		lin.isUniqueSolution = false;
    		lin.x = null;
    		lin.Ainv = null;
    		
    	} else if (this.existsAllZeroRows(Aplus, lin)) {
    		
    		lin.solutionExists = true;
    		lin.isUniqueSolution = false;
    		lin.x = null;
    		lin.Ainv = null;
    		
    	} else {
    		lin.isUniqueSolution = true;
    		lin.solutionExists = true;
    		
    		lin.x = new double[lin.rank];
    		
    		for (int i = 0; i < lin.x.length; i++) {
    			lin.x[i] = Aplus[i][numOfCols-1];
    		}
    		
    		//Include the inverse, given a unique solution.
    		if (lin.rank != Aplus.length) {
    			lin.Ainv = null;
    		} else if (!inverseAvailable) {
    			lin.Ainv = null;
    		} else {
    			lin.Ainv = I;
    		}
    	}
    	
    	return lin;
    }
    
    /**
     * Finds the pivot after computing the REF.
     * @param l
     * @param k
     * @return
     */
    private int[] findPivot(LinResult l, int k, double[][] Aplus) {
    	
    	int[] pivotIndex = new int[2];
    	pivotIndex[0] = l.pivotRow[k-1];
    	
    	for (int i = pivotIndex[0]; i < Aplus[pivotIndex[0]].length-1; i++) {
    		
    		if (l.isPivotColumn[i]) {
    			if (Aplus[pivotIndex[0]][i] == 1) {
    				pivotIndex[1] = i;
    				return pivotIndex;
    			}
    		}
    	}
    	
    	return pivotIndex;
    }
    
    /**
     * Returns a boolean indicating there are any contradictions in the RREF.
     * @param Aplus
     * @return
     */
    private boolean existsContradictionRow(double[][] Aplus, LinResult l) {
    	
    	double temp = 0.0;
    	boolean pivotCheck = false;
    	
    	for (int i = 0; i < Aplus.length; i++) {
    		
    		if (Aplus[i][Aplus[i].length-1] != 0) {
        		
    			for (int j = 0; j < Aplus[i].length-1; j++) {
        			temp += Aplus[i][j];
        			if (temp >= 1) {
        				pivotCheck = true;
        			}
        		}
    			if (temp == 0.0 && !pivotCheck) {
    				return true;
    			} else {
    				temp = 0.0;
    				pivotCheck = false;
    			}
    		}
    	}
    	
    	return false;
    }
    
    /**
     * Returns a boolean indicating whether there are multiple solutions in the RREF.
     * @param Aplus
     * @return
     */
    private boolean existsAllZeroRows(double[][] Aplus, LinResult l) {
    	
    	double temp = 0.0;
    	boolean pivotCheck = false;
    	
    	for (int i = 0; i < Aplus.length; i++) {
    		
    		for (int j = 0; j < Aplus[i].length-1; j++) {
    			
    			temp += Aplus[i][j];
    			if (temp >= 1) {
    				pivotCheck = true;
    			}
    		}
    		
    		if ((temp > 0 && ( temp != 1)) || temp < 0) {
    			return true;
    		} else if ((temp == 0.0) && (Aplus.length <= 2 || pivotCheck)) {
    			return true;
    		} else {
    			temp = 0.0;
    			pivotCheck = false;
    		}
    	}
    	
    	return false;
    }
    
    /**
     * Generates identity matrix of size n x n.
     * @param n
     * @return
     */
    private double[][] generateIdentityMatrix(int n) {
    	
    	double[][] I = new double[n][n];
    	
    	int row = 0;
    	int col = 0;
    	
    	for (int i = 0; i < n; i++) {
    		for (int j = 0; j < n; j++) {
    			if (i == row && j == col) {
    				I[i][j] = 1.0;
    				row = row + 1;
    				col = col + 1;
    			} else {
    				I[i][j] = 0.0;
    			}
    		}
    	}
    	
    	return I;
    }
    
    /**
     * Utility method used for debugging.
     * @param A
     */
    private void printMatrix(double[][] A) {
    	
    	for (int i = 0; i < A.length; i++) {
    		System.out.print("[");
    		for (int j = 0; j < A[i].length; j++) {
    			System.out.print(A[i][j] + ", ");
    		}
    		System.out.print("]\n");
    	}
    }
    
    /**
     * Utility method for copying.
     * @param A
     * @return
     */
    private double[][] makeDeepCopy(double[][] A) {
    	
    	int rows = A.length;
    	int cols = A[0].length;
    	
    	double[][] copy = new double[rows][cols];
    	
    	for (int i = 0; i < rows; i++) {
    		for (int j = 0; j < cols; j++) {
    			copy[i][j] = A[i][j];
    		}
    	}
    	
    	return copy;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Assignment 3 Code for LinTool
    
    // Orthogonality:
    public boolean areColumnsLI (double[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null. Aborting.");
    		return false;
    	}
    	
//    	if (Arows != Acols) {
//    		System.out.println("Matrix is not n x n dimensions.");
//    		return false;
//    	}
    	
    	double[] b = new double[A.length];
    	
    	for (int i = 0; i < b.length; i++) {
    		b[i] = 0.0;
    	}
    	
    	double[][] Aplus = this.computeAugmentedMatrix(A, b);
    	
    	LinResult lin = this.privateComputeRREF(Aplus);
    	
    	if (lin.isUniqueSolution) {
    		return true;
    	}
    	
    	return false;
    }
    
    
    public LinResult gramSchmidt (double[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null. Aborting.");
    		return null;
    	}
    	
    	//Transpose of A so that each row can be easily treated like a vector.
    	double[][] aT = this.transpose(A);
    	LinResult ret = new LinResult();
    	
    	// Transposed matrix that contains the orthogonal vectors from G-S. 
    	// This will be re-transposed to produce the proper return matrix.
    	double[][] bT = new double[aT.length][aT[0].length];
    	
    	// Transposed matrix that contains the normalized orthogonal vectors from G-S.
    	// This will be re-transposed to produce the proper return matrix.
    	double[][] cT = new double[aT.length][aT[0].length];
    	
    	// Temporary vectors to be used in the algorithm.
    	double[] s = new double[aT[0].length];
    	double[] sJ = null;
    	
    	double projVal = 0.0;
    	
    	// v_1 = w_1
    	bT[0] = aT[0];
    	cT[0] = this.scalarProduct((1 / this.norm(bT[0])), bT[0]);
    	
    	for (int i = 1; i < aT.length; i++) {
    		
        	for (int a = 0; a < s.length; a++) {
        		s[a] = 0.0;
        	}
    		
    		for (int j = 0; j < i; j++) {
    			projVal = this.proj(aT[i], bT[j]);
    			sJ = this.scalarProduct(projVal, bT[j]);
    			s = this.add(s, sJ);
    		}
    		
    		bT[i] = this.add(aT[i], this.makeVectorNegative(s));
    		cT[i] = this.scalarProduct((1 / this.norm(bT[i])), bT[i]);
    	}
    	
    	double[][] b = this.transpose(bT);
    	ret.V = b;
    	
    	double[][] c = this.transpose(cT);
    	ret.Q = c;
    	
    	return ret;
    }
    
    public LinResult computeQR (double[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null. Aborting.");
    		return null;
    	}
    	
    	//this.printMatrix(A);
    	
    	//Transpose of A so that each row can be easily treated like a vector.
    	double[][] aT = this.transpose(A);
    	LinResult ret = new LinResult();
    	
    	// Transposed matrix that contains the orthogonal vectors from G-S. 
    	// This will be re-transposed to produce the proper return matrix.
    	double[][] bT = new double[aT.length][aT[0].length];
    	
    	// Transposed matrix that contains the normalized orthogonal vectors from G-S.
    	// This will be re-transposed to produce the proper return matrix.
    	double[][] cT = new double[aT.length][aT[0].length];
    	
    	// Transposed matrix that contains the r values from QR decomposition algorithm.
    	// This will be re-transposed to produce the proper return matrix.
    	double[][] dT = new double[aT.length][aT[0].length];
    	
    	// Temporary vectors to be used in the algorithm.
    	double[] s = new double[aT[0].length];
    	double[] sJ = null;
    	
    	double projVal = 0.0;
    	
    	// v_1 = w_1
    	bT[0] = aT[0];
    	cT[0] = this.scalarProduct((1 / this.norm(bT[0])), bT[0]);
    	dT[0][0] = this.norm(aT[0]);
    	
    	for (int i = 1; i < aT.length; i++) {
    		
        	for (int a = 0; a < s.length; a++) {
        		s[a] = 0.0;
        	}
    		
    		for (int j = 0; j < i; j++) {
    			projVal = this.proj(aT[i], bT[j]);
    			
    			System.out.println(projVal + "");
    			sJ = this.scalarProduct(projVal, bT[j]);
    			s = this.add(s, sJ);
    			dT[i][j] = projVal * this.norm(bT[j]); 
    		}
    		
    		bT[i] = this.add(aT[i], this.makeVectorNegative(s));
    		cT[i] = this.scalarProduct((1 / this.norm(bT[i])), bT[i]);
    		dT[i][i] = this.norm(bT[i]);
    	}
    	
    	// Debugging Code:
//    	System.out.println("END");
//    	this.printMatrix(bT);
//    	System.out.println("END");
//    	this.printMatrix(cT);
//    	System.out.println("END");
//    	this.printMatrix(dT);
    	
    	double[][] b = this.transpose(bT);
    	ret.V = b;
    	
    	double[][] c = this.transpose(cT);
    	ret.Q = c;
    	
    	double[][] d = this.transpose(dT);
    	ret.R = d;
    	
    	return ret;
    }

    // Complex matrix operations:
    public ComplexNumber[][] add (ComplexNumber[][] A, ComplexNumber[][] B) {
    	
    	if (A == null || B == null) {
    		System.out.println("Matrices are null. Aborting");
    		return null;
    	}
    	
    	if (A.length != B.length || A[0].length != B[0].length) {
    		System.out.println("Matrices are of differing dimensions and cannot be added.");
    		return null;
    	}
    	
    	ComplexNumber retMatrix[][] = new ComplexNumber[A.length][A[0].length];
    	
    	for (int i = 0; i < A.length; i++) {
    		for (int j = 0; j < A[i].length; j++) {
    			retMatrix[i][j] = A[i][j].add(B[i][j]);
    		}
    	}
    	
    	return retMatrix;
    }
    public ComplexNumber[][] scalarProduct (ComplexNumber alpha, ComplexNumber[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null.");
    		return null;
    	}
    	
    	ComplexNumber retMatrix[][] = new ComplexNumber[A.length][A[0].length];
        
    	for (int i = 0; i < A.length; i++) {
    		for (int j = 0; j < A[i].length; j++) {
    			retMatrix[i][j] = A[i][j].mult(alpha);
    		}
    	}
    	
    	return retMatrix;
    }
    
    public ComplexNumber[][] mult (ComplexNumber[][] A, ComplexNumber[][] B) {
    	
		if (A == null || B == null) {
			System.out.println("Matrices are empty");
			return null;
		}
		
		int Arows = A.length;
		int Acolumns = A[0].length; 
		
		int Brows = B.length;
		int Bcolumns = B[0].length;
		
		if (Acolumns != Brows) {
			System.out.println("Invalid matrix dimensions. Aborting program.");
			return null;
		}
		
		this.printMatrix(A);
		this.printMatrix(B);
		ComplexNumber[][] ret = new ComplexNumber[Arows][Bcolumns];
		ComplexNumber tempValue = ComplexNumber.makeComplexNumber(0, 0);
		
		for (int i = 0; i < Arows; i++) {
			for (int j = 0; j < Bcolumns; j++) {
				for (int k = 0; k < Brows; k++) {
					tempValue = tempValue.add(A[i][k].mult(B[k][j]));
				}
				ret[i][j] = tempValue;
				tempValue = ComplexNumber.makeComplexNumber(0, 0);
			}
		}
		this.printMatrix(ret);
		return ret;
    }
    
    public ComplexNumber[] matrixVectorMult (ComplexNumber[][] A, ComplexNumber[] v) {
    	
    	if (A == null || v == null) {
    		System.out.println("Matrix or vector is empty. Aborting.");
    		return null;
    	}
    	
		ComplexNumber[] ret = new ComplexNumber[A.length];
		ComplexNumber tempValue = ComplexNumber.makeComplexNumber(0, 0);
		
		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < A[i].length; j++) {
				tempValue = tempValue.add(A[i][j].mult(v[j]));
			}
			ret[i] = tempValue;
			tempValue = ComplexNumber.makeComplexNumber(0, 0);
		}
		return ret;
    }
    
    public ComplexNumber[][] hermitianTranspose (ComplexNumber[][] A) {
    	
    	if (A == null) {
    		System.out.println("Matrix is null. Aborting.");
    		return null;
    	}
    	
    	ComplexNumber[][] B = new ComplexNumber[A[0].length][A.length];
    	
    	for (int i = 0; i < B.length; i++) {
    		for (int j = 0; j < B[i].length; j++) {
    			B[i][j] = A[j][i];
    			B[i][j].im = B[i][j].im * -1;
    		}
    	}
    	
    	return B;
    }

    /**
     * Utility method used for debugging.
     * @param A
     */
    private void printMatrix(ComplexNumber[][] A) {
    	
    	for (int i = 0; i < A.length; i++) {
    		System.out.print("[");
    		for (int j = 0; j < A[i].length; j++) {
    			System.out.print(A[i][j].re + "+" + A[i][j].im + "i, ");
    		}
    		System.out.print("]\n");
    	}
    	System.out.println("END");
    }
    
    /**
     * Helper method to calculate projected value for G-S algorithm.
     * @param w
     * @param v
     * @return
     */
    private double proj (double[] w, double[] v) {
		
    	double dpNum = this.dotProduct(w, v);
    	double dpDen = this.dotProduct(v, v);
	
		double scalarValue = dpNum / dpDen;

		return scalarValue; // Temporarily
    }
    
    private double[] makeVectorNegative(double[] v) {
    	
    	for (int i = 0; i < v.length; i++) {
    		v[i] = v[i] * -1;
    	}
    	
    	return v;
    }
    
}
