import edu.gwu.lintool.*;

import java.text.*;
import java.util.*;
import java.io.*;


/**
 * Rohan Sunder
 * CSCI 6342 - Linear Algebra: A Computational Approach
 * @author rsunder
 *
 */
public class TestA2 {

    public static void main (String[] argv)
    {
    	
    RohanA2_LinTool lin = new RohanA2_LinTool();
	// Test complex vector operations.
	LinTest.testComplexVectors (lin);

	// REPLACE the line below with your tool.
	//LinTool lin = new LinToolImpl ();

	//LinTest.testREF (lin);
	//LinTest.testRREF (lin);
	//LinTest.testSolveFromREF (lin);
	//LinTest.testSolveFromRREF (lin);
	LinTest.testInverse (lin);
    }

}