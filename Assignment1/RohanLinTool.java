import edu.gwu.lintool.LinToolEmpty;

/**
 * Rohan Sunder
 * CSCI 6342 - Linear Algebra: A Computational Approach
 * Assignment 1
 * @author rsunder
 *
 */
public class RohanLinTool extends LinToolEmpty {

    //Vector Operations
	
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
    	
		//
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
    	
    	RohanLinTool tool = new RohanLinTool();
    	
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
    	
		// Add the square of every value in the matrix. Then take the square root.
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
}
