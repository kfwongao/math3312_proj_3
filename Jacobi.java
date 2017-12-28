/**
 * MATH3312 Project 3
 * name: WONG, Kin Fat Huanqiang
 * this code has some modifications from the
 * source code found in GitHub which
 * uploaded by Rafael Rendón Pablo <smart.rendon@gmail.com>
 */

/* * This class provides a simple implementation of the Jacobi method for solving
 * systems of linear equations. */

/*
  How to use:
  The program reads an augmented matrix from standard input,
  for example:
 
   3
   5 -2  3 -1
  -3  9  1  2
   2 -1 -7  3
  The number in the first line is the number of equations
  and number of variables. You can put this values in a file
  and then execute the program as follows:
  $ java Jacobi < equations.txt
  If the matrix isn't diagonally dominant the program tries
  to convert it(if possible) by rearranging the rows.
*/


import java.io.*;
import java.util.Arrays;
import java.util.StringTokenizer;

public class Jacobi {

  public static int maxNumOfIteration = 30;
  public static double tol = 0.001;
  private double[][] M;
  
  public Jacobi(double [][] matrix) { M = matrix; }

  public void print()
  {
	System.out.println("This is what you have inputed");
	System.out.println("maxNumOfIteration = " + maxNumOfIteration + "\ntol = " + tol);
	System.out.println("Ax=b in augement form: ");
    int n = M.length;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n + 1; j++)
        System.out.print(M[i][j] + " ");
      System.out.println();
    }
  }

  public boolean transformToDominant(int r, boolean[] V, int[] R)
  {
    int n = M.length;
    if (r == M.length) {
      double[][] T = new double[n][n+1];
      for (int i = 0; i < R.length; i++) {
        for (int j = 0; j < n + 1; j++)
          T[i][j] = M[R[i]][j];
      }

      M = T;
      
      return true;
    }

    for (int i = 0; i < n; i++) {
      if (V[i]) continue;

      double sum = 0;
      
      for (int j = 0; j < n; j++)
        sum += Math.abs(M[i][j]);

      if (2 * Math.abs(M[i][r]) > sum) { // diagonally dominant?
        V[i] = true;
        R[r] = i;

        if (transformToDominant(r + 1, V, R))
          return true;

        V[i] = false;
      }
    }

    return false;
  }
  
  
  /**
   * Returns true if is possible to transform M(data member) to a diagonally
   * dominant matrix, false otherwise.
  */
  public boolean makeDominant()
  {
    boolean[] visited = new boolean[M.length];
    int[] rows = new int[M.length];

    Arrays.fill(visited, false);

    return transformToDominant(0, visited, rows);
  }


  /**
   * Applies Jacobi method to find the solution of the system
   * of linear equations represented in matrix M.
   * M is a matrix with the following form:
   * a_11 * x_1 + a_12 * x_2 + ... + a_1n * x_n = b_1
   * a_21 * x_1 + a_22 * x_2 + ... + a_2n * x_n = b_2
   * .                 .                  .        .
   * .                 .                  .        .
   * .                 .                  .        .
   * a_n1 * x_n + a_n2 * x_2 + ... + a_nn * x_n = b_n
  */
  public void solve()
  {
    int iterations = 0;
    int n = M.length;
    double epsilon = tol;
    double[] X = new double[n]; // Approximations
    double[] P = new double[n]; // Prev
    Arrays.fill(X, 0);
    Arrays.fill(P, 0);

    System.out.println("The solution in each iteration is:\n");
    System.out.print("The starting x = (");
    for (int i = 0; i < n; i++)
      System.out.print("0" + " ");
    System.out.println(")");
    
    while (true) {
      for (int i = 0; i < n; i++) {
        double sum = M[i][n]; // b_n

        for (int j = 0; j < n; j++)
          if (j != i)
            sum -= M[i][j] * P[j];

        // Update x_i but it's no used in the next row calculation
        // but up to de next iteration of the method
        X[i] = 1/M[i][i] * sum;
      }
      
      System.out.print("The " + (iterations+1) + "th iteration x = (");
      for (int i = 0; i < n; i++)
        System.out.print(X[i] + " ");
      System.out.println(")");

      iterations++;
//      if (iterations == 1) continue;

      boolean stop = true;
      for (int i = 0; i < n && stop; i++)
        if (Math.abs(X[i] - P[i]) > epsilon)
          stop = false;

      if (stop || iterations == maxNumOfIteration) break;
      P = (double[])X.clone();
    }
  }

  public static void main(String[] args) throws IOException
  {
    int n;
    double[][] M;

    BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
    PrintWriter writer = new PrintWriter(System.out, true);
    System.out.println("Please input the maximum no. of iteration you need to solve the system Ax=b");
    maxNumOfIteration = Integer.parseInt(reader.readLine());
    
    System.out.println("Please input the tolerance");
    tol = Double.parseDouble(reader.readLine());
    
    System.out.println("Please input n, the size of the nxn matrix A");
    n = Integer.parseInt(reader.readLine());
    
    M = new double[n][n+1];
    System.out.println("Please input the i th linear equation in Ax=b "
			+ "\nYou only need to input the coefficients "
			+ "\nExample: −x1 + 11x2 − x3 + 3x4 = 25, you only need to input -1 11 -1 3 25\n");
    
    for (int i = 0; i < n; i++) {
    	System.out.println("Please input the " + (i+1) + "th linear equation in Ax=b ");
    		StringTokenizer strtk = new StringTokenizer(reader.readLine());

    		while (strtk.hasMoreTokens())
    	  for (int j = 0; j < n + 1 && strtk.hasMoreTokens(); j++)
    		  M[i][j] = Integer.parseInt(strtk.nextToken());
    }


    
    
    
    
    Jacobi jacobi = new Jacobi(M);

    if (!jacobi.makeDominant()) {
      writer.println("The system isn't diagonally dominant: " + 
                     "The method cannot guarantee convergence.");
    }

    writer.println();

    jacobi.print();

    jacobi.solve();

  }
}

