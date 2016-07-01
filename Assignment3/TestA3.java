import edu.gwu.lintool.*;

import java.text.*;
import java.util.*;
import java.io.*;

public class TestA3 {

    public static void main (String[] argv)
    {
    	
    	// REPLACE the line below with your tool.
    	LinTool lin = new RohanA3_LinTool ();
    	// Test complex matrix operations.
    	LinTest.testComplexMatrices (lin);

    	LinTest.testAreColumnsLI (lin);
    	LinTest.testGramSchmidt (lin);
    	LinTest.testQR (lin);
    }

}