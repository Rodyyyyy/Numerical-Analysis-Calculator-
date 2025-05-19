package com.example;

public class System_Linear_Non {
    public void jacobi(int n, double[][] A, double[] B, double[] x, int iterations) {
        double[] xNew = new double[n];

        for (int itr = 1; itr <= iterations; itr++) {
            for (int i = 0; i < n; i++) {
                double sum = B[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= A[i][j] * x[j];
                    }
                }
                xNew[i] = sum / A[i][i];
            }

            System.out.print("Iteration " + itr + ": ");
            for (int i = 0; i < n; i++) {
                System.out.printf("x[%d] = %.6f ", i, xNew[i]);
                x[i] = xNew[i];
            }
            System.out.println();
        }
    }

    public void gauss(int n, double[][] A, double[] B, double[] x, int iterations) {
        for (int itr = 1; itr <= iterations; itr++) {
            for (int i = 0; i < n; i++) {
                double sum = B[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= A[i][j] * x[j];
                    }
                }
                x[i] = sum / A[i][i];
            }

            System.out.print("Iteration " + itr + ": ");
            for (int i = 0; i < n; i++) {
                System.out.printf("x[%d] = %.6f ", i, x[i]);
            }
            System.out.println();
        }
    }

    public void newton() {
        System.out.println("Non-linear system not supported in GUI yet.");
    }
}